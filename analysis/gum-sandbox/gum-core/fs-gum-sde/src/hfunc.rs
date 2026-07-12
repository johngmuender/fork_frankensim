//! Coarse-grained H-function and decay fit — port of born.py's
//! `h_function` and `fit_tau`.
//!
//! H̄ = Σ_cells P̄ ln(P̄/Q̄) with P̄ from particle counts and Q̄ from
//! box-averaging the fine |ψ|² grid; empty cells (P̄ = 0) are SKIPPED and
//! Q̄ is floored at 1e-300 (born.py conventions, verbatim). The finite-N
//! sampling noise floor is `cg²/(2N)` (occupied cells / 2N).
//!
//! All reductions are straight sequential sums in fixed (row-major) order
//! and all logs go through fs-math `det` — bit-identical replay.

use fs_math::det;

use crate::modes::BOX_L;

/// Normalized coarse-grained particle histogram P̄ on a cg × cg grid over
/// [0, π]², row-major `p[i*cg + j]` (i indexes x).
#[must_use]
pub fn histogram_p(xs: &[f64], ys: &[f64], cg: usize) -> Vec<f64> {
    assert_eq!(xs.len(), ys.len());
    let mut counts = vec![0u64; cg * cg];
    let scale = cg as f64 / BOX_L;
    for (&x, &y) in xs.iter().zip(ys) {
        let i = ((x * scale) as usize).min(cg - 1);
        let j = ((y * scale) as usize).min(cg - 1);
        counts[i * cg + j] += 1;
    }
    let n = xs.len() as f64;
    counts.iter().map(|&c| c as f64 / n).collect()
}

/// Q̄: block-mean of the fine |ψ|² grid (nfine × nfine, row-major) onto
/// cg × cg, normalized to Σ Q̄ = 1. `nfine` must be a multiple of `cg`.
#[must_use]
pub fn coarse_q(q_fine: &[f64], nfine: usize, cg: usize) -> Vec<f64> {
    assert_eq!(q_fine.len(), nfine * nfine);
    assert_eq!(nfine % cg, 0, "nfine must be a multiple of cg");
    let b = nfine / cg;
    let mut q = vec![0.0f64; cg * cg];
    for bi in 0..cg {
        for bj in 0..cg {
            let mut sum = 0.0;
            for i in 0..b {
                let row = (bi * b + i) * nfine + bj * b;
                for j in 0..b {
                    sum += q_fine[row + j];
                }
            }
            q[bi * cg + bj] = sum; // block mean = sum/b²; normalization absorbs b²
        }
    }
    let tot: f64 = q.iter().sum();
    for v in &mut q {
        *v /= tot;
    }
    q
}

/// H̄ = Σ_{P̄>0} P̄ ln(P̄ / max(Q̄, 1e-300)).
#[must_use]
pub fn h_bar(p: &[f64], q: &[f64]) -> f64 {
    assert_eq!(p.len(), q.len());
    let mut h = 0.0;
    for (&pv, &qv) in p.iter().zip(q) {
        if pv > 0.0 {
            h += pv * det::ln(pv / qv.max(1e-300));
        }
    }
    h
}

/// Finite-N noise floor ≈ occupied cells / (2N) = cg²/(2N).
#[must_use]
pub fn noise_floor(cg: usize, n_part: usize) -> f64 {
    (cg * cg) as f64 / (2.0 * n_part as f64)
}

/// Exponential-decay fit result (born.py `fit_tau`).
#[derive(Debug, Clone, Copy)]
pub struct TauFit {
    /// Decay time from the log-space linear fit H ≈ H0_fit e^{−t/τ}.
    pub tau: f64,
    /// exp(intercept).
    pub h0_fit: f64,
    /// Fit-window end actually used.
    pub t_end: f64,
    /// R² of the linear fit in log space.
    pub r2: f64,
}

/// Fit H ~ H0 exp(−t/τ) over the initial decay: up to where H first drops
/// to 10% of H(0), or t = 2π, whichever comes first (born.py `fit_tau`,
/// verbatim port including the ≥ 4-point requirement).
#[must_use]
pub fn fit_tau(ts: &[f64], hs: &[f64]) -> TauFit {
    assert_eq!(ts.len(), hs.len());
    let h0 = hs[0];
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut t_end = two_pi;
    for (i, &h) in hs.iter().enumerate() {
        if h <= 0.1 * h0 {
            t_end = ts[i].min(two_pi);
            break;
        }
    }
    let sel: Vec<(f64, f64)> = ts
        .iter()
        .zip(hs)
        .filter(|&(&t, &h)| t <= t_end && h > 0.0)
        .map(|(&t, &h)| (t, det::ln(h)))
        .collect();
    if sel.len() < 4 {
        return TauFit {
            tau: f64::NAN,
            h0_fit: f64::NAN,
            t_end,
            r2: f64::NAN,
        };
    }
    let n = sel.len() as f64;
    let (mut sx, mut sy, mut sxx, mut sxy) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for &(t, ly) in &sel {
        sx += t;
        sy += ly;
        sxx += t * t;
        sxy += t * ly;
    }
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let intercept = (sy - slope * sx) / n;
    let mean_y = sy / n;
    let (mut ss_res, mut ss_tot) = (0.0f64, 0.0f64);
    for &(t, ly) in &sel {
        let r = ly - (slope * t + intercept);
        ss_res += r * r;
        ss_tot += (ly - mean_y) * (ly - mean_y);
    }
    let r2 = if ss_tot > 0.0 {
        1.0 - ss_res / ss_tot
    } else {
        f64::NAN
    };
    TauFit {
        tau: -1.0 / slope,
        h0_fit: det::exp(intercept),
        t_end,
        r2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn h_bar_zero_when_p_matches_q() {
        let p = vec![0.25; 4];
        let q = vec![0.25; 4];
        assert_eq!(h_bar(&p, &q), 0.0);
    }

    #[test]
    fn h_bar_positive_and_empty_cells_skipped() {
        let p = vec![0.5, 0.5, 0.0, 0.0];
        let q = vec![0.25; 4];
        let h = h_bar(&p, &q);
        assert!((h - det::ln(2.0)).abs() < 1e-14);
    }

    #[test]
    fn fit_tau_recovers_synthetic_decay() {
        let tau_true = 3.0;
        let ts: Vec<f64> = (0..80).map(|i| 0.1 * i as f64).collect();
        let hs: Vec<f64> = ts.iter().map(|&t| 0.5 * det::exp(-t / tau_true)).collect();
        let f = fit_tau(&ts, &hs);
        assert!((f.tau - tau_true).abs() < 1e-6, "tau {}", f.tau);
        assert!((f.h0_fit - 0.5).abs() < 1e-9);
        assert!(f.r2 > 0.999999);
    }

    #[test]
    fn noise_floor_formula() {
        assert!((noise_floor(32, 20000) - 0.0256).abs() < 1e-15);
        assert!((noise_floor(16, 20000) - 0.0064).abs() < 1e-15);
    }
}
