//! Minimal internal radix-2 FFT (power-of-two sizes, 3-D separable).
//!
//! The survey's first choice was crates/fs-fft (`FftNd`), but at build
//! time its dependency closure pulled fs-exec -> asupersync, an EXTERNAL
//! sibling repository that is not present in this tree — so fs-fft was
//! not usable from a standalone out-of-workspace crate, and the build
//! spec's fallback applied: a minimal internal radix-2 FFT.
//!
//! PHASE F3 UPDATE: fs-fft IS now importable from here — depend on it as
//! `fs-fft = { path = "../../../../crates/fs-fft", default-features =
//! false }` (its default-ON `exec` feature carries the fs-exec/asupersync
//! coupling; disabling it keeps the whole serial `Fft`/`FftNd`/`RealFft`
//! surface). The swap is deliberately NOT performed in F3: this internal
//! FFT is golden-gated, and swapping kernels risks output changes that
//! would force a golden re-freeze. Keep the structured `NotPow2`
//! rejection when swapping — see gum-core/PORTABILITY.md.
//!
//! Deterministic by
//! construction: precomputed twiddles, fixed butterfly and axis order.
//! Forward is unnormalized (X_m = sum_j x_j e^{-2 pi i m j / n});
//! `inverse` scales by 1/n so inverse(forward(x)) = x.

use core::f64::consts::PI;

/// A complex number (f64 re/im).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl C64 {
    /// Construct.
    #[inline]
    #[must_use]
    pub const fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }

    #[inline]
    fn mul(self, o: C64) -> C64 {
        C64::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }

    #[inline]
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }

    #[inline]
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}

fn is_pow2(x: usize) -> bool {
    x >= 1 && (x & (x - 1)) == 0
}

/// In-place iterative radix-2 Cooley–Tukey with bit-reversal permutation.
fn fft1d(data: &mut [C64], twiddles: &[C64]) {
    let n = data.len();
    if n <= 1 {
        return;
    }
    // Bit-reversal permutation (fixed order).
    let mut j = 0usize;
    for i in 0..n - 1 {
        if i < j {
            data.swap(i, j);
        }
        let mut m = n >> 1;
        while m >= 1 && (j & m) != 0 {
            j ^= m;
            m >>= 1;
        }
        j |= m;
    }
    // Butterflies, fixed stage/group order. twiddles[k] = e^{-2 pi i k / n},
    // k < n/2; stage with half-size `hs` uses stride n / (2 hs).
    let mut hs = 1usize;
    while hs < n {
        let stride = n / (2 * hs);
        let mut base = 0usize;
        while base < n {
            for k in 0..hs {
                let w = twiddles[k * stride];
                let a = data[base + k];
                let b = data[base + k + hs].mul(w);
                data[base + k] = a.add(b);
                data[base + k + hs] = a.sub(b);
            }
            base += 2 * hs;
        }
        hs *= 2;
    }
}

/// A planned 3-D complex FFT over a row-major buffer (last axis
/// contiguous). Every axis length must be a power of two.
pub struct Fft3 {
    dims: [usize; 3],
    /// Forward twiddles per axis: e^{-2 pi i k / n}, k < n/2.
    tw_f: [Vec<C64>; 3],
    /// Inverse twiddles per axis (conjugates).
    tw_i: [Vec<C64>; 3],
}

impl Fft3 {
    /// Plan a transform.
    ///
    /// # Panics
    /// If any axis is not a power of two.
    #[must_use]
    pub fn new(dims: [usize; 3]) -> Fft3 {
        assert!(dims.iter().all(|&d| is_pow2(d)), "axes must be powers of two: {dims:?}");
        let mk = |n: usize, sign: f64| -> Vec<C64> {
            (0..n / 2)
                .map(|k| {
                    let ang = sign * 2.0 * PI * k as f64 / n as f64;
                    C64::new(ang.cos(), ang.sin())
                })
                .collect()
        };
        Fft3 {
            dims,
            tw_f: [mk(dims[0], -1.0), mk(dims[1], -1.0), mk(dims[2], -1.0)],
            tw_i: [mk(dims[0], 1.0), mk(dims[1], 1.0), mk(dims[2], 1.0)],
        }
    }

    fn run(&self, data: &mut [C64], inverse: bool) {
        let [n0, n1, n2] = self.dims;
        assert_eq!(data.len(), n0 * n1 * n2, "buffer length mismatch");
        let tw = if inverse { &self.tw_i } else { &self.tw_f };
        // Axis 2 (contiguous pencils).
        for row in data.chunks_exact_mut(n2) {
            fft1d(row, &tw[2]);
        }
        // Axis 1: gather/scatter pencils, fixed order.
        let mut pencil = vec![C64::default(); n1.max(n0)];
        for i in 0..n0 {
            for k in 0..n2 {
                for j in 0..n1 {
                    pencil[j] = data[(i * n1 + j) * n2 + k];
                }
                fft1d(&mut pencil[..n1], &tw[1]);
                for j in 0..n1 {
                    data[(i * n1 + j) * n2 + k] = pencil[j];
                }
            }
        }
        // Axis 0.
        for j in 0..n1 {
            for k in 0..n2 {
                for i in 0..n0 {
                    pencil[i] = data[(i * n1 + j) * n2 + k];
                }
                fft1d(&mut pencil[..n0], &tw[0]);
                for i in 0..n0 {
                    data[(i * n1 + j) * n2 + k] = pencil[i];
                }
            }
        }
        if inverse {
            let s = 1.0 / (n0 * n1 * n2) as f64;
            for v in data.iter_mut() {
                v.re *= s;
                v.im *= s;
            }
        }
    }

    /// Forward 3-D DFT, unnormalized.
    pub fn forward(&self, data: &mut [C64]) {
        self.run(data, false);
    }

    /// Inverse 3-D DFT with 1/total normalization.
    pub fn inverse(&self, data: &mut [C64]) {
        self.run(data, true);
    }
}
