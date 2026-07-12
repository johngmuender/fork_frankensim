//! Diagnostics kit and the deterministic perturbation generator.
//!
//! Auxiliary cell fields (velocities, gradients, perturbations, tangent
//! directions) share one flat layout with `engine::Grad`:
//! `idx(a, i, j, k) = ((a n + i) n + j) n + k` over the REAL cells.

use fs_gum_field::Field3;
use fs_math::det;

/// Flat index into an auxiliary cell field (Grad layout).
#[inline]
#[must_use]
pub fn aidx(n: usize, a: usize, i: usize, j: usize, k: usize) -> usize {
    ((a * n + i) * n + j) * n + k
}

// ---------------------------------------------------------------------------
// deterministic RNG (hand-rolled LCG; documented constants)
// ---------------------------------------------------------------------------

/// 64-bit linear congruential generator, Knuth MMIX constants
/// a = 6364136223846793005, c = 1442695040888963407 (mod 2^64).  Uniform
/// doubles take the top 53 bits.  Four warm-up steps decorrelate small
/// seeds.  Chosen per the survey's dictum: determinism matters more than
/// generator pedigree — the stream is part of the frozen protocol, not a
/// statistical claim.
pub struct Lcg(u64);

impl Lcg {
    #[must_use]
    pub fn new(seed: u64) -> Self {
        let mut s = Lcg(seed);
        for _ in 0..4 {
            s.next_u64();
        }
        s
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    /// Uniform in [0, 1) with 53-bit resolution.
    pub fn u01(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Uniform in [lo, hi).
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.u01()
    }
}

// ---------------------------------------------------------------------------
// perturbation fields (field3d perturbation() / halo_seed(), LCG-driven)
// ---------------------------------------------------------------------------

/// Fixed-seed smooth non-axisymmetric bump field (field3d
/// `Engine.perturbation` semantics): `kbumps` Gaussian bumps, centres
/// uniform in [-1.25 R*, 1.25 R*]^3, widths in [0.5, 0.9), 4-component
/// amplitudes in [-amp, amp).  Draw order per bump (frozen): cx, cy, cz,
/// w, A0, A1, A2, A3.  Returns an auxiliary cell field (Grad layout).
#[must_use]
pub fn bump_field(n: usize, half: f64, seed: u64, kbumps: usize, amp: f64) -> Vec<f64> {
    let rstar = fs_gum_field::rstar();
    let h = 2.0 * half / (n as f64);
    let x = |i: usize| -half + (i as f64 + 0.5) * h;
    let mut rng = Lcg::new(seed);
    let mut dq = vec![0.0_f64; 4 * n * n * n];
    for _ in 0..kbumps {
        let cx = rng.uniform(-1.25 * rstar, 1.25 * rstar);
        let cy = rng.uniform(-1.25 * rstar, 1.25 * rstar);
        let cz = rng.uniform(-1.25 * rstar, 1.25 * rstar);
        let w = rng.uniform(0.5, 0.9);
        let mut a = [0.0_f64; 4];
        for av in a.iter_mut() {
            *av = rng.uniform(-amp, amp);
        }
        let i2w2 = 1.0 / (2.0 * w * w);
        for i in 0..n {
            let dx = x(i) - cx;
            for j in 0..n {
                let dy = x(j) - cy;
                for k in 0..n {
                    let dz = x(k) - cz;
                    let gau = det::exp(-(dx * dx + dy * dy + dz * dz) * i2w2);
                    for (c, &av) in a.iter().enumerate() {
                        dq[aidx(n, c, i, j, k)] += av * gau;
                    }
                }
            }
        }
    }
    dq
}

/// Explicit small tilt-halo shell in q1 (field3d `Engine.halo_seed`):
/// dq1 = eta exp(-(r - r0)^2 / (2 w^2)).
#[must_use]
pub fn halo_seed_field(n: usize, half: f64, eta: f64, r0: f64, w: f64) -> Vec<f64> {
    let h = 2.0 * half / (n as f64);
    let x = |i: usize| -half + (i as f64 + 0.5) * h;
    let i2w2 = 1.0 / (2.0 * w * w);
    let mut dq = vec![0.0_f64; 4 * n * n * n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let r = (x(i) * x(i) + x(j) * x(j) + x(k) * x(k)).sqrt();
                dq[aidx(n, 1, i, j, k)] = eta * det::exp(-(r - r0) * (r - r0) * i2w2);
            }
        }
    }
    dq
}

// ---------------------------------------------------------------------------
// field <-> auxiliary-buffer plumbing
// ---------------------------------------------------------------------------

/// Deep copy of a field (real cells; ghosts are vacuum in both).
#[must_use]
pub fn clone_field(f: &Field3) -> Field3 {
    let n = f.n();
    let mut c = Field3::new(n, f.half());
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                c.set(i, j, k, f.get(i, j, k));
            }
        }
    }
    c
}

/// Snapshot the real cells into an auxiliary buffer (Grad layout).
pub fn save_cells(f: &Field3, buf: &mut [f64]) {
    let n = f.n();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let q = f.get(i, j, k);
                for (a, &v) in q.iter().enumerate() {
                    buf[aidx(n, a, i, j, k)] = v;
                }
            }
        }
    }
}

/// Restore the real cells from a snapshot.
pub fn restore_cells(f: &mut Field3, buf: &[f64]) {
    let n = f.n();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let mut q = [0.0_f64; 4];
                for (a, qv) in q.iter_mut().enumerate() {
                    *qv = buf[aidx(n, a, i, j, k)];
                }
                f.set(i, j, k, q);
            }
        }
    }
}

/// q += scale * dq on the real cells (no renormalisation; call
/// `f.renormalize()` afterwards, matching Python `normalize(q + dq)`).
pub fn add_scaled(f: &mut Field3, dq: &[f64], scale: f64) {
    let n = f.n();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let mut q = f.get(i, j, k);
                for (a, qv) in q.iter_mut().enumerate() {
                    *qv += scale * dq[aidx(n, a, i, j, k)];
                }
                f.set(i, j, k, q);
            }
        }
    }
}

/// Pointwise tangent projection of an auxiliary field: u -= (u.q) q.
pub fn tangent_project(f: &Field3, u: &mut [f64]) {
    let n = f.n();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let q = f.get(i, j, k);
                let mut dot = 0.0_f64;
                for (a, &qv) in q.iter().enumerate() {
                    dot += u[aidx(n, a, i, j, k)] * qv;
                }
                for (a, &qv) in q.iter().enumerate() {
                    u[aidx(n, a, i, j, k)] -= dot * qv;
                }
            }
        }
    }
}

/// Normalise a flat buffer to unit Euclidean norm; returns the pre-fix norm.
pub fn normalize_l2(u: &mut [f64]) -> f64 {
    let mut s = 0.0_f64;
    for &v in u.iter() {
        s += v * v;
    }
    let nrm = s.sqrt();
    for v in u.iter_mut() {
        *v /= nrm;
    }
    nrm
}

/// Fixed-order dot product of two flat buffers.
#[must_use]
pub fn dot_flat(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0_f64;
    for (&x, &y) in a.iter().zip(b.iter()) {
        s += x * y;
    }
    s
}

// ---------------------------------------------------------------------------
// halo fraction + clock (field3d halo_fraction / Section K)
// ---------------------------------------------------------------------------

/// I-fraction outside r > 1.5 R* around the (1 - q0) centroid
/// (field3d `Engine.halo_fraction`).  Returns (fraction, centroid).
#[must_use]
pub fn halo_fraction(f: &Field3) -> (f64, [f64; 3]) {
    let n = f.n();
    let (mut wsum, mut cx, mut cy, mut cz) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let q = f.get(i, j, k);
                let w = 1.0 - q[0];
                wsum += w;
                cx += w * f.x(i);
                cy += w * f.x(j);
                cz += w * f.x(k);
            }
        }
    }
    cx /= wsum;
    cy /= wsum;
    cz /= wsum;
    let hr = crate::halo_radius();
    let hr2 = hr * hr;
    let (mut tot, mut halo) = (0.0_f64, 0.0_f64);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let q = f.get(i, j, k);
                let di = q[1] * q[1] + q[2] * q[2];
                tot += di;
                let (dx, dy, dz) = (f.x(i) - cx, f.x(j) - cy, f.x(k) - cz);
                if dx * dx + dy * dy + dz * dz > hr2 {
                    halo += di;
                }
            }
        }
    }
    (halo / tot, [cx, cy, cz])
}

/// Clock bisection at frozen fields: the L solving L^2 = (2/3) I E_static
/// (200 bisection steps on [0, 20], the field3d Section-K procedure; the
/// closed form sqrt((2/3) I E) is the cross-check, not the primary).
#[must_use]
pub fn clock_bisect(i_val: f64, estat: f64) -> f64 {
    let (mut lo, mut hi) = (0.0_f64, 20.0_f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid * mid - (2.0 / 3.0) * i_val * estat > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Rotational energy fraction E_rot / E_tot = (L^2/2I) / (E_static + L^2/2I).
#[must_use]
pub fn erot_frac(l: f64, i_val: f64, estat: f64) -> f64 {
    let erot = l * l / (2.0 * i_val);
    erot / (estat + erot)
}
