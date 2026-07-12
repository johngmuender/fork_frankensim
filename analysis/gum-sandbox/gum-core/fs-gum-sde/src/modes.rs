//! Mode-sum wavefunction machinery — port of born.py's `ModeSystem`.
//!
//! 2-D infinite square well [0, π]², ħ = m = 1:
//!
//! ```text
//! φ_mn(x, y) = (2/π) sin(m x) sin(n y),   E_mn = (m² + n²)/2
//! ψ(x, y, t) = M^{-1/2} Σ_{(m,n)∈S} e^{i θ_mn} φ_mn e^{-i E_mn t}
//! ```
//!
//! Exact evolution (no PDE solve): only the complex coefficients rotate.
//! ψ and ∇ψ are evaluated through the separable bilinear form
//! ψ = sxᵀ C sy with per-axis sine/cosine tables, exactly as born.py's
//! `ModeSystem.velocity` (cost O(kmax) per point, O(kmax²) per coefficient
//! refresh).
//!
//! The fixed random phases of the Tier-3 pilot (numpy `default_rng(42)`,
//! `uniform(0, 2π, 16)`) are embedded as [`BORN_THETA_16`] so the Rust ψ is
//! the SAME wavefunction born.py evolved — the constants are shortest
//! round-trip decimal reprs, which parse to bit-identical f64s.

use fs_math::det;

/// Box side L = π.
pub const BOX_L: f64 = std::f64::consts::PI;

/// born.py's 16 mode phases (numpy seed 42), index `(m-1)*4 + (n-1)`
/// matching `meshgrid(ks, ks, indexing="ij").ravel()` order. Verified
/// against `np.random.default_rng(42).uniform(0, 2π, 16)` bit patterns.
pub const BORN_THETA_16: [f64; 16] = [
    4.862909272689599,
    2.757554564287996,
    5.3947298351621535,
    4.381692553882582,
    0.5917337285168199,
    6.13001602516006,
    4.782381792256834,
    4.938987693414485,
    0.8049616944763924,
    2.829858307545725,
    2.3297926977893746,
    5.823036161141988,
    4.0455238622962515,
    5.169563679814652,
    2.7860535790666945,
    1.42778299794038,
];

/// Density floor used by born.py when dividing by |ψ|² (guards nodal lines).
pub const DENS_FLOOR: f64 = 1e-30;

/// Superposition of the first kmax × kmax box eigenmodes, equal weights,
/// fixed phases.
#[derive(Debug, Clone)]
pub struct ModeSystem {
    /// Modes per axis (M = kmax²).
    pub kmax: usize,
    /// Energies E[(m-1)*kmax + (n-1)] = (m² + n²)/2.
    pub e: Vec<f64>,
    /// Phases θ[(m-1)*kmax + (n-1)].
    pub theta: Vec<f64>,
}

/// Time-dependent complex coefficient matrix C[m-1][n-1] =
/// exp(i(θ − E t))/√M, row-major, split into real/imaginary parts.
/// The (2/π) eigenmode norm is NOT included (it cancels in ∇ψ/ψ);
/// density evaluations multiply it back in explicitly.
#[derive(Debug, Clone)]
pub struct ModeCoeffs {
    /// Modes per axis.
    pub kmax: usize,
    /// Re C, row-major (m-1)*kmax + (n-1).
    pub cr: Vec<f64>,
    /// Im C, same layout.
    pub ci: Vec<f64>,
}

/// The complex log-derivative W = ∇ψ/ψ at a point, plus the (un-normalized,
/// floored) density |ψ|². In ħ = m = 1 units the Nelson velocities are
/// v = Im W (current) and u = Re W (osmotic) — see [`crate::nelson`].
#[derive(Debug, Clone, Copy)]
pub struct LogDeriv {
    /// Re(∂ₓψ/ψ).
    pub wx_re: f64,
    /// Im(∂ₓψ/ψ).
    pub wx_im: f64,
    /// Re(∂ᵧψ/ψ).
    pub wy_re: f64,
    /// Im(∂ᵧψ/ψ).
    pub wy_im: f64,
    /// max(|ψ|², 1e-30) WITHOUT the (2/π)² norm (born.py convention).
    pub dens: f64,
}

/// Hard cap on kmax so per-point tables live on the stack.
pub const KMAX_MAX: usize = 8;

impl ModeSystem {
    /// New system from per-mode phases (`theta.len() == kmax²`).
    #[must_use]
    pub fn new(kmax: usize, theta: &[f64]) -> Self {
        assert!(kmax >= 1 && kmax <= KMAX_MAX, "kmax out of range");
        assert_eq!(theta.len(), kmax * kmax, "need kmax² phases");
        let mut e = Vec::with_capacity(kmax * kmax);
        for m in 1..=kmax {
            for n in 1..=kmax {
                e.push(0.5 * ((m * m + n * n) as f64));
            }
        }
        ModeSystem {
            kmax,
            e,
            theta: theta.to_vec(),
        }
    }

    /// The Tier-3 M = 16 system (kmax = 4, born.py phases, seed 42).
    #[must_use]
    pub fn born_m16() -> Self {
        Self::new(4, &BORN_THETA_16)
    }

    /// The single stationary ground mode φ₁₁ (M = 1; the phase is a global
    /// factor and drops out of every observable).
    #[must_use]
    pub fn phi11() -> Self {
        Self::new(1, &[0.0])
    }

    /// Number of modes M = kmax².
    #[must_use]
    pub fn m_modes(&self) -> usize {
        self.kmax * self.kmax
    }

    /// Coefficient matrix at time t (pure function of t — safe to cache and
    /// share across trajectories; see `nelson::CoeffCache`).
    #[must_use]
    pub fn coeffs(&self, t: f64) -> ModeCoeffs {
        let m = self.m_modes();
        let inv_sqrt_m = 1.0 / det::sqrt(m as f64);
        let mut cr = Vec::with_capacity(m);
        let mut ci = Vec::with_capacity(m);
        for i in 0..m {
            let phase = self.theta[i] - self.e[i] * t;
            cr.push(det::cos(phase) * inv_sqrt_m);
            ci.push(det::sin(phase) * inv_sqrt_m);
        }
        ModeCoeffs {
            kmax: self.kmax,
            cr,
            ci,
        }
    }

    /// W = ∇ψ/ψ and |ψ|² at one point via the separable bilinear form
    /// (born.py `velocity`, with the real part — the osmotic direction —
    /// returned as well; born.py only formed the imaginary part).
    #[must_use]
    pub fn log_derivative(&self, c: &ModeCoeffs, x: f64, y: f64) -> LogDeriv {
        debug_assert_eq!(c.kmax, self.kmax);
        let kmax = self.kmax;
        let mut sx = [0.0f64; KMAX_MAX];
        let mut cxm = [0.0f64; KMAX_MAX];
        let mut sy = [0.0f64; KMAX_MAX];
        let mut cym = [0.0f64; KMAX_MAX];
        for k in 0..kmax {
            let a = (k + 1) as f64;
            sx[k] = det::sin(a * x);
            cxm[k] = a * det::cos(a * x); // d/dx of sin(a x)
            sy[k] = det::sin(a * y);
            cym[k] = a * det::cos(a * y);
        }
        let (mut pr, mut pi) = (0.0f64, 0.0f64);
        let (mut gxr, mut gxi) = (0.0f64, 0.0f64);
        let (mut gyr, mut gyi) = (0.0f64, 0.0f64);
        for j in 0..kmax {
            // a[j] = Σ_i sx[i] C[i][j] (ψ and ∂ᵧψ share it);
            // bx[j] = Σ_i cxm[i] C[i][j] (∂ₓψ).
            let (mut arj, mut aij, mut bxr, mut bxi) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            for i in 0..kmax {
                let re = c.cr[i * kmax + j];
                let im = c.ci[i * kmax + j];
                arj += sx[i] * re;
                aij += sx[i] * im;
                bxr += cxm[i] * re;
                bxi += cxm[i] * im;
            }
            pr += arj * sy[j];
            pi += aij * sy[j];
            gxr += bxr * sy[j];
            gxi += bxi * sy[j];
            gyr += arj * cym[j];
            gyi += aij * cym[j];
        }
        let dens = (pr * pr + pi * pi).max(DENS_FLOOR);
        LogDeriv {
            // W = ∇ψ · conj(ψ) / |ψ|².
            wx_re: (gxr * pr + gxi * pi) / dens,
            wx_im: (gxi * pr - gxr * pi) / dens, // born.py's vx
            wy_re: (gyr * pr + gyi * pi) / dens,
            wy_im: (gyi * pr - gyr * pi) / dens, // born.py's vy
            dens,
        }
    }

    /// |ψ(t)|² at the centers of an n × n fine grid, row-major
    /// `out[i*n + j] = |ψ(x_i, y_j)|²`, x_i = (i + ½)π/n — includes the
    /// (2/π) eigenmode norm (born.py `psi2_grid`).
    #[must_use]
    pub fn psi2_grid(&self, t: f64, n: usize) -> Vec<f64> {
        let c = self.coeffs(t);
        let kmax = self.kmax;
        let norm = 2.0 / std::f64::consts::PI;
        // sxg[i*kmax + k] = sin((k+1) x_i).
        let mut sxg = vec![0.0f64; n * kmax];
        for i in 0..n {
            let x = (i as f64 + 0.5) * BOX_L / n as f64;
            for k in 0..kmax {
                sxg[i * kmax + k] = det::sin((k + 1) as f64 * x);
            }
        }
        // tmp[i*kmax + l] = Σ_k sxg[i][k] C[k][l] * norm.
        let mut tr = vec![0.0f64; n * kmax];
        let mut ti = vec![0.0f64; n * kmax];
        for i in 0..n {
            for l in 0..kmax {
                let (mut ar, mut ai) = (0.0f64, 0.0f64);
                for k in 0..kmax {
                    let s = sxg[i * kmax + k];
                    ar += s * c.cr[k * kmax + l];
                    ai += s * c.ci[k * kmax + l];
                }
                tr[i * kmax + l] = ar * norm;
                ti[i * kmax + l] = ai * norm;
            }
        }
        let mut out = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let (mut pr, mut pi) = (0.0f64, 0.0f64);
                for l in 0..kmax {
                    let s = sxg[j * kmax + l];
                    pr += tr[i * kmax + l] * s;
                    pi += ti[i * kmax + l] * s;
                }
                out[i * n + j] = pr * pr + pi * pi;
            }
        }
        out
    }

    /// |ψ(t)|² at one point, WITH the (2/π) norm (born.py `psi2_points`;
    /// rejection-sampling kernel).
    #[must_use]
    pub fn psi2_point(&self, c: &ModeCoeffs, x: f64, y: f64) -> f64 {
        let kmax = self.kmax;
        let norm = 2.0 / std::f64::consts::PI;
        let mut sx = [0.0f64; KMAX_MAX];
        let mut sy = [0.0f64; KMAX_MAX];
        for k in 0..kmax {
            let a = (k + 1) as f64;
            sx[k] = det::sin(a * x);
            sy[k] = det::sin(a * y);
        }
        let (mut pr, mut pi) = (0.0f64, 0.0f64);
        for i in 0..kmax {
            for j in 0..kmax {
                let w = sx[i] * sy[j];
                pr += w * c.cr[i * kmax + j];
                pi += w * c.ci[i * kmax + j];
            }
        }
        let (pr, pi) = (pr * norm, pi * norm);
        pr * pr + pi * pi
    }

    /// Rejection-sampling envelope for |ψ(t)|²: 1.2 × the max over a
    /// 512 × 512 grid (born.py `sample_psi2` bound).
    #[must_use]
    pub fn psi2_max_bound(&self, t: f64) -> f64 {
        let g = self.psi2_grid(t, 512);
        let mut mx = 0.0f64;
        for &v in &g {
            if v > mx {
                mx = v;
            }
        }
        1.2 * mx
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn psi2_grid_is_normalized() {
        // ∫|ψ|² over the box = 1 (equal weights, orthonormal modes);
        // midpoint rule on 256² is exact-ish for this trig polynomial.
        let sys = ModeSystem::born_m16();
        for &t in &[0.0, 1.7] {
            let n = 256;
            let g = sys.psi2_grid(t, n);
            let cell = (BOX_L / n as f64) * (BOX_L / n as f64);
            let total: f64 = g.iter().sum::<f64>() * cell;
            assert!((total - 1.0).abs() < 1e-9, "norm {total} at t={t}");
        }
    }

    #[test]
    fn phi11_log_derivative_is_real_cotangent() {
        // ψ = φ₁₁: W = (cot x, cot y), purely real (v = 0, u = cot).
        let sys = ModeSystem::phi11();
        let c = sys.coeffs(3.21);
        let (x, y) = (0.7, 2.1);
        let ld = sys.log_derivative(&c, x, y);
        assert!((ld.wx_re - det::cos(x) / det::sin(x)).abs() < 1e-12);
        assert!((ld.wy_re - det::cos(y) / det::sin(y)).abs() < 1e-12);
        assert!(ld.wx_im.abs() < 1e-12);
        assert!(ld.wy_im.abs() < 1e-12);
    }

    #[test]
    fn psi2_point_matches_grid() {
        let sys = ModeSystem::born_m16();
        let t = 0.9;
        let c = sys.coeffs(t);
        let n = 64;
        let g = sys.psi2_grid(t, n);
        let (i, j) = (17, 41);
        let x = (i as f64 + 0.5) * BOX_L / n as f64;
        let y = (j as f64 + 0.5) * BOX_L / n as f64;
        let p = sys.psi2_point(&c, x, y);
        assert!((p - g[i * n + j]).abs() < 1e-12 * p.max(1.0));
    }
}
