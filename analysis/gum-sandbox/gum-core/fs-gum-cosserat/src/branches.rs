//! Branch extraction, classification, co-motion fidelity, and small-k fits
//! (spectrum.py `branches` / `classify_full` / `fit_w2`, ported 1:1 on the
//! k ∥ z gauge where the exact decoupling u_z | φ_z | transverse holds; the
//! full 6×6 solve is cross-validated against the sub-block union at every k).

use crate::eigh::solve;
use crate::moduli::Moduli;
use crate::symbol::{mass_matrix, stiffness};
use fs_math::c64::C64;

/// Exactly-decoupled sub-blocks for k ∥ z (spectrum.py order).
pub const IDX_UZ: [usize; 1] = [2];
/// Longitudinal φ_z block.
pub const IDX_PZ: [usize; 1] = [5];
/// Transverse block order: u_x, u_y, φ_x, φ_y.
pub const IDX_T: [usize; 4] = [0, 1, 3, 4];

/// Branch arrays over a k grid plus diagnostics (spectrum.py `branches`).
pub struct Branches {
    /// Longitudinal acoustic ω²(k) (u_z block).
    pub b1: Vec<f64>,
    /// Light transverse doublet ω² (both members, ascending).
    pub b2: Vec<[f64; 2]>,
    /// Heavy transverse doublet ω² (both members, ascending).
    pub b3: Vec<[f64; 2]>,
    /// Longitudinal twist ω²(k) (φ_z block).
    pub b4: Vec<f64>,
    /// Co-motion fidelity r = |φ_⊥| / ((k/2)|u_⊥|) on the lowest transverse mode.
    pub r_light: Vec<f64>,
    /// Max |Δω²| between the full 6×6 spectrum and the sub-block union.
    pub xval: f64,
    /// Most negative ω² encountered.
    pub min_w2: f64,
    /// Worst embedding pair residual encountered.
    pub pair_resid: f64,
}

/// Evaluate the four branches over the k grid, k ∥ z.
#[must_use]
pub fn branches(m: &Moduli, kk: &[f64]) -> Branches {
    let md = mass_matrix(m);
    let n = kk.len();
    let mut out = Branches {
        b1: vec![0.0; n],
        b2: vec![[0.0; 2]; n],
        b3: vec![[0.0; 2]; n],
        b4: vec![0.0; n],
        r_light: vec![0.0; n],
        xval: 0.0,
        min_w2: f64::INFINITY,
        pair_resid: 0.0,
    };
    for (i, &k) in kk.iter().enumerate() {
        let kmat = stiffness([0.0, 0.0, k], m);
        let full = solve(&kmat, &md, None);
        let w1 = solve(&kmat, &md, Some(&IDX_UZ));
        let w4 = solve(&kmat, &md, Some(&IDX_PZ));
        let wt = solve(&kmat, &md, Some(&IDX_T));
        // Union of sub-block eigenvalues, sorted, vs the full spectrum.
        let mut union: Vec<f64> = Vec::with_capacity(6);
        union.extend_from_slice(&w1.w2);
        union.extend_from_slice(&w4.w2);
        union.extend_from_slice(&wt.w2);
        union.sort_by(f64::total_cmp);
        for (a, b) in union.iter().zip(full.w2.iter()) {
            let d = (a - b).abs();
            if d > out.xval {
                out.xval = d;
            }
        }
        for &w in full.w2.iter().chain(wt.w2.iter()) {
            if w < out.min_w2 {
                out.min_w2 = w;
            }
        }
        for pr in [full.pair_resid, w1.pair_resid, w4.pair_resid, wt.pair_resid] {
            if pr > out.pair_resid {
                out.pair_resid = pr;
            }
        }
        out.b1[i] = w1.w2[0];
        out.b4[i] = w4.w2[0];
        out.b2[i] = [wt.w2[0], wt.w2[1]];
        out.b3[i] = [wt.w2[2], wt.w2[3]];
        // Lowest transverse mode: r = |φ_⊥| / ((k/2)|u_⊥|).
        let v = &wt.vecs[0]; // order u_x, u_y, φ_x, φ_y
        let u_amp = (v[0].norm_sq() + v[1].norm_sq()).sqrt();
        let p_amp = (v[2].norm_sq() + v[3].norm_sq()).sqrt();
        out.r_light[i] = if u_amp > 0.0 { p_amp / (0.5 * k * u_amp) } else { f64::NAN };
    }
    out
}

/// Least-squares fit ω² = gap + slope·k² over k ≤ kmax (centered normal
/// equations — numerically equivalent to spectrum.py's polyfit here).
#[must_use]
pub fn fit_w2(kk: &[f64], w2: &[f64], kmax: f64) -> (f64, f64) {
    let mut n = 0.0f64;
    let mut sx = 0.0f64;
    let mut sy = 0.0f64;
    for (i, &k) in kk.iter().enumerate() {
        if k <= kmax {
            n += 1.0;
            sx += k * k;
            sy += w2[i];
        }
    }
    let (mx, my) = (sx / n, sy / n);
    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    for (i, &k) in kk.iter().enumerate() {
        if k <= kmax {
            let dx = k * k - mx;
            sxx += dx * dx;
            sxy += dx * (w2[i] - my);
        }
    }
    let slope = sxy / sxx;
    let gap = my - slope * mx;
    (gap, slope)
}

/// Least-squares fit ln y = amp + p·ln x over the masked points
/// (for the χ₃ splitting power law). Returns (exponent p, amplitude).
#[must_use]
pub fn fit_loglog(x: &[f64], y: &[f64], mask: impl Fn(usize) -> bool) -> (f64, f64) {
    let mut n = 0.0f64;
    let mut sx = 0.0f64;
    let mut sy = 0.0f64;
    for i in 0..x.len() {
        if mask(i) {
            n += 1.0;
            sx += x[i].ln();
            sy += y[i].ln();
        }
    }
    let (mx, my) = (sx / n, sy / n);
    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    for i in 0..x.len() {
        if mask(i) {
            let dx = x[i].ln() - mx;
            sxx += dx * dx;
            sxy += dx * (y[i].ln() - my);
        }
    }
    let p = sxy / sxx;
    (p, my - p * mx)
}

/// Classify the 6 eigenmodes of the full 6×6 problem at wavevector k·k̂:
/// longitudinal-u / longitudinal-φ / transverse dominant-subspace weights
/// (reduces to spectrum.py `classify_full` for k̂ = e_z).
#[must_use]
pub fn classify_full(k: f64, khat: [f64; 3], m: &Moduli) -> Vec<(f64, &'static str, f64)> {
    let kmat = stiffness([k * khat[0], k * khat[1], k * khat[2]], m);
    let md = mass_matrix(m);
    let modes = solve(&kmat, &md, None);
    let mut rows = Vec::with_capacity(6);
    for n in 0..6 {
        let v = &modes.vecs[n];
        // Complex longitudinal projections u·k̂, φ·k̂.
        let mut ul = C64::ZERO;
        let mut pl = C64::ZERO;
        let mut tot = 0.0f64;
        for a in 0..3 {
            ul = ul + v[a].scale(khat[a]);
            pl = pl + v[3 + a].scale(khat[a]);
            tot += v[a].norm_sq() + v[3 + a].norm_sq();
        }
        let wl_u = ul.norm_sq();
        let wl_p = pl.norm_sq();
        let wt = tot - wl_u - wl_p;
        let fr = [wl_u / tot, wl_p / tot, wt / tot];
        let (mut best, mut lab) = (fr[0], "long-u");
        if fr[1] > best {
            best = fr[1];
            lab = "long-phi";
        }
        if fr[2] > best {
            best = fr[2];
            lab = "transverse";
        }
        rows.push((modes.w2[n], lab, best));
    }
    rows
}
