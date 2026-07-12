//! 1-D radial hedgehog solver at fixed t — the Rust port of
//! `tier2-closure/radial_solve.py` (graded grid, midpoint assembly,
//! diagonally preconditioned energy-monotone flow, damped (Levenberg)
//! Newton with the exact symmetric-tridiagonal Hessian and Thomas solves,
//! eps-dial outer iteration on t).
//!
//! Sector energies (radial integrals), total E = t (E2 + E4) + E6 + E0
//! (a6 = a0 = m = 1):
//! ```text
//! E2 = INT (f'^2 r^2 + 2 sin^2 f) dr
//! E4 = INT sin^2 f (2 f'^2 + sin^2 f / r^2) dr
//! E6 = INT sin^4 f f'^2 / r^2 dr
//! E0 = INT (1 - cos f) r^2 dr
//! ```
//! Deterministic: no RNG, fixed sweep orders, elementary functions from
//! `fs_math::det`.

use fs_math::det;

/// Solved radial profile and its sector energies.
pub struct RadialProfile {
    /// Grid nodes r_0 = 0 .. r_N = rmax (graded).
    pub r: Vec<f64>,
    /// Profile values f(r_i), f(0) = pi, f(rmax) = 0 (Dirichlet).
    pub f: Vec<f64>,
    /// (E2, E4, E6, E0) at the solution.
    pub sectors: [f64; 4],
    /// Coupling t = a2 = a4 the profile was solved at.
    pub t: f64,
    /// Achieved max |grad| (Newton floor).
    pub gmax: f64,
}

/// np.interp semantics: piecewise-linear, clamped to the end values outside
/// the (monotonically increasing) node range.
#[must_use]
pub fn interp(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    let n = xp.len();
    if x <= xp[0] {
        return fp[0];
    }
    if x >= xp[n - 1] {
        return fp[n - 1];
    }
    // last index i with xp[i] <= x (binary search on the sorted nodes)
    let i = xp.partition_point(|&v| v <= x) - 1;
    fp[i] + (fp[i + 1] - fp[i]) * (x - xp[i]) / (xp[i + 1] - xp[i])
}

/// Graded grid on [0, rmax]: N cells, N+1 nodes.  Node density boosted near
/// r = 0 (core), near the compacton edge r ~ R*, and through the Yukawa-tail
/// window; deterministic inverse-CDF construction (200001-sample trapezoid,
/// exactly the Python weight function).
#[must_use]
pub fn make_grid(n: usize, rmax: f64) -> Vec<f64> {
    const M: usize = 200_001;
    let rstar = crate::rstar();
    let mut rr = vec![0.0_f64; M];
    let mut w = vec![0.0_f64; M];
    for (i, (rv, wv)) in rr.iter_mut().zip(w.iter_mut()).enumerate() {
        let s = i as f64 / (M - 1) as f64;
        let r = s * rmax;
        *rv = r;
        let a = r / 0.5;
        let b = (r - rstar) / 0.45;
        let c = (r - (rstar + 1.2)) / 0.9;
        *wv = 0.55 + 2.0 * det::exp(-(a * a)) + 3.0 * det::exp(-(b * b)) + 0.8 * det::exp(-(c * c));
    }
    let mut cw = vec![0.0_f64; M];
    for i in 1..M {
        cw[i] = cw[i - 1] + 0.5 * (w[i] + w[i - 1]) * (rr[i] - rr[i - 1]);
    }
    let total = cw[M - 1];
    for v in cw.iter_mut() {
        *v /= total;
    }
    let mut r = vec![0.0_f64; n + 1];
    for (j, rv) in r.iter_mut().enumerate() {
        let u = j as f64 / n as f64;
        *rv = interp(u, &cw, &rr);
    }
    r[0] = 0.0;
    r[n] = rmax;
    r
}

/// Smoothed-compacton initial profile: f = 2 arccos(u / (1 + u^8)^(1/8)),
/// u = r/R* ((1+u^8)^(1/8) via three IEEE square roots).
#[must_use]
pub fn init_profile(r: &[f64]) -> Vec<f64> {
    let rstar = crate::rstar();
    let n = r.len() - 1;
    let mut f = vec![0.0_f64; n + 1];
    for (fv, &rv) in f.iter_mut().zip(r.iter()) {
        let u = rv / rstar;
        let u2 = u * u;
        let u8 = u2 * u2 * u2 * u2;
        let x = u / (1.0 + u8).sqrt().sqrt().sqrt();
        *fv = 2.0 * det::acos(x.clamp(0.0, 1.0));
    }
    f[0] = core::f64::consts::PI;
    f[n] = 0.0;
    f
}

/// Output of one midpoint assembly.
pub struct Assembly {
    /// Total energy t (E2 + E4) + E6 + E0.
    pub e: f64,
    /// (E2, E4, E6, E0).
    pub sectors: [f64; 4],
    /// Gradient at the interior nodes 1..N-1 (length N-1).
    pub g: Vec<f64>,
    /// Hessian diagonal at the interior nodes (length N-1; empty if
    /// `need_hess` was false).
    pub diag: Vec<f64>,
    /// Hessian off-diagonal couplings (length N-2; empty if `need_hess`
    /// was false).
    pub off: Vec<f64>,
}

/// Midpoint (cell-centred) assembly of the discrete energy, its exact
/// gradient, and (optionally) the exact symmetric-tridiagonal Hessian.
/// Per cell: e_c = h_c L(r_m, f_m, d), f_m = (f_i + f_{i+1})/2,
/// d = (f_{i+1} - f_i)/h_c (radial_solve.py `assemble`, same expressions).
#[must_use]
pub fn assemble(r: &[f64], f: &[f64], t: f64, need_hess: bool) -> Assembly {
    let nc = r.len() - 1; // cells
    let mut e2 = 0.0_f64;
    let mut e4 = 0.0_f64;
    let mut e6 = 0.0_f64;
    let mut e0 = 0.0_f64;
    let mut gl = vec![0.0_f64; nc];
    let mut gr = vec![0.0_f64; nc];
    let (mut hll, mut hrr, mut hlr) = if need_hess {
        (vec![0.0_f64; nc], vec![0.0_f64; nc], vec![0.0_f64; nc])
    } else {
        (Vec::new(), Vec::new(), Vec::new())
    };
    for c in 0..nc {
        let h = r[c + 1] - r[c];
        let rm = 0.5 * (r[c + 1] + r[c]);
        let rm2 = rm * rm;
        let fm = 0.5 * (f[c + 1] + f[c]);
        let d = (f[c + 1] - f[c]) / h;
        let d2 = d * d;
        let s = det::sin(fm);
        let cc = det::cos(fm);
        let s2 = s * s;
        let s3 = s2 * s;
        let s4 = s2 * s2;
        let sin2f = 2.0 * s * cc;
        let cos2f = cc * cc - s2;
        e2 += h * (d2 * rm2 + 2.0 * s2);
        e4 += h * (s2 * (2.0 * d2 + s2 / rm2));
        e6 += h * (s4 * d2 / rm2);
        e0 += h * ((1.0 - cc) * rm2);
        let l_d = 2.0 * t * d * rm2 + 4.0 * t * s2 * d + 2.0 * s4 * d / rm2;
        let l_f = 2.0 * t * sin2f * (1.0 + d2) + 4.0 * (t + d2) * s3 * cc / rm2 + s * rm2;
        gl[c] = 0.5 * h * l_f - l_d;
        gr[c] = 0.5 * h * l_f + l_d;
        if need_hess {
            let l_dd = 2.0 * t * rm2 + 4.0 * t * s2 + 2.0 * s4 / rm2;
            let l_fd = 4.0 * t * d * sin2f + 8.0 * s3 * cc * d / rm2;
            let l_ff = 4.0 * t * cos2f * (1.0 + d2)
                + 4.0 * (t + d2) * (3.0 * s2 * cc * cc - s4) / rm2
                + cc * rm2;
            hll[c] = 0.25 * h * l_ff - l_fd + l_dd / h;
            hrr[c] = 0.25 * h * l_ff + l_fd + l_dd / h;
            hlr[c] = 0.25 * h * l_ff - l_dd / h;
        }
    }
    let mut g = vec![0.0_f64; nc - 1];
    for i in 0..nc - 1 {
        g[i] = gr[i] + gl[i + 1];
    }
    let (mut diag, mut off) = (Vec::new(), Vec::new());
    if need_hess {
        diag = vec![0.0_f64; nc - 1];
        for i in 0..nc - 1 {
            diag[i] = hrr[i] + hll[i + 1];
        }
        off = vec![0.0_f64; nc.saturating_sub(2)];
        for i in 0..nc.saturating_sub(2) {
            off[i] = hlr[i + 1];
        }
    }
    let e = t * (e2 + e4) + e6 + e0;
    Assembly { e, sectors: [e2, e4, e6, e0], g, diag, off }
}

/// Solve the symmetric tridiagonal system (diag; off) x = b (Thomas
/// algorithm, sequential, deterministic).
#[must_use]
pub fn thomas(diag: &[f64], off: &[f64], b: &[f64]) -> Vec<f64> {
    let n = diag.len();
    let mut dd = diag.to_vec();
    let mut bb = b.to_vec();
    for i in 1..n {
        let m = off[i - 1] / dd[i - 1];
        dd[i] -= m * off[i - 1];
        bb[i] -= m * bb[i - 1];
    }
    let mut x = vec![0.0_f64; n];
    x[n - 1] = bb[n - 1] / dd[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (bb[i] - off[i] * x[i + 1]) / dd[i];
    }
    x
}

fn gmax_of(g: &[f64]) -> f64 {
    let mut m = 0.0_f64;
    for &v in g {
        let a = v.abs();
        if a > m {
            m = a;
        }
    }
    m
}

/// Diagonally preconditioned, energy-monotone gradient flow to the Newton
/// basin (|g|_max < `gmax_target`).
pub fn flow(r: &[f64], f: &mut Vec<f64>, t: f64, gmax_target: f64, maxsteps: usize) {
    let n = f.len() - 1;
    let mut asm = assemble(r, f, t, true);
    let mut dt = 0.5;
    for _ in 0..maxsteps {
        if gmax_of(&asm.g) < gmax_target {
            break;
        }
        let mut dmax = f64::NEG_INFINITY;
        for &v in &asm.diag {
            if v > dmax {
                dmax = v;
            }
        }
        let mut fnew = f.clone();
        for i in 1..n {
            let prec = asm.diag[i - 1].max(1.0e-3 * dmax);
            fnew[i] = (f[i] - dt * asm.g[i - 1] / prec).clamp(0.0, core::f64::consts::PI);
        }
        let asm_new = assemble(r, &fnew, t, true);
        if asm_new.e < asm.e {
            *f = fnew;
            asm = asm_new;
            dt = (dt * 1.1).min(5.0);
        } else {
            dt *= 0.5;
            if dt < 1.0e-9 {
                break;
            }
        }
    }
}

/// Damped (Levenberg) Newton on grad(E) = 0 with the exact tridiagonal
/// Hessian; a step is accepted only if max |grad| decreases.  Returns
/// (sectors, gmax, |dE| of the last accepted step, iterations).
pub fn newton(
    r: &[f64],
    f: &mut Vec<f64>,
    t: f64,
    gtol: f64,
    maxit: usize,
) -> ([f64; 4], f64, f64, usize) {
    let n = f.len() - 1;
    let mut asm = assemble(r, f, t, true);
    let mut lam = 1.0e-3;
    let mut de_last = f64::INFINITY;
    let mut it = 0;
    while it < maxit {
        let gmax = gmax_of(&asm.g);
        if gmax < gtol {
            break;
        }
        let mut accepted = false;
        while lam < 1.0e14 {
            let mut damped = asm.diag.clone();
            for (dv, &d0) in damped.iter_mut().zip(asm.diag.iter()) {
                *dv = d0 + lam * (d0.abs() + 1.0e-30);
            }
            let rhs: Vec<f64> = asm.g.iter().map(|&v| -v).collect();
            let step = thomas(&damped, &asm.off, &rhs);
            if step.iter().all(|v| v.is_finite()) {
                let mut fnew = f.clone();
                for i in 1..n {
                    fnew[i] += step[i - 1];
                }
                let asm_new = assemble(r, &fnew, t, true);
                let gmax_new = gmax_of(&asm_new.g);
                if asm_new.e.is_finite() && gmax_new < gmax {
                    de_last = (asm.e - asm_new.e).abs();
                    *f = fnew;
                    asm = asm_new;
                    lam = (lam * 0.25).max(1.0e-14);
                    accepted = true;
                    break;
                }
            }
            lam *= 10.0;
        }
        if !accepted {
            break; // stalled at the machine floor
        }
        it += 1;
    }
    (asm.sectors, gmax_of(&asm.g), de_last, it)
}

/// Full profile solve at fixed t: graded grid, smoothed-compacton init,
/// flow to the basin, Newton to max |grad| < 1e-12.
#[must_use]
pub fn radial_solve(t: f64, n: usize, rmax: f64) -> RadialProfile {
    let r = make_grid(n, rmax);
    let mut f = init_profile(&r);
    flow(&r, &mut f, t, 1.0e-2, 6000);
    let (sectors, gmax, _de, _it) = newton(&r, &mut f, t, 1.0e-12, 300);
    RadialProfile { r, f, sectors, t, gmax }
}

/// Outer eps-dial iteration on t: re-solve and update
/// t <- t * target / ratio (ratio = t (E2+E4)/(E6+E0) is ~linear in t, so
/// this fixed-point secant converges in a few steps; radial_solve.py outer
/// loop, verbatim).  Returns the dialled t, the final profile, and the
/// (t, ratio) history.
#[must_use]
pub fn eps_dial(
    t0: f64,
    target: f64,
    tol: f64,
    n: usize,
    rmax: f64,
    max_outer: usize,
) -> (f64, RadialProfile, Vec<(f64, f64)>) {
    let r = make_grid(n, rmax);
    let mut f = init_profile(&r);
    let mut t = t0;
    let mut history = Vec::new();
    let mut sectors = [0.0_f64; 4];
    let mut gmax = 0.0_f64;
    for _ in 0..max_outer {
        flow(&r, &mut f, t, 1.0e-2, 6000);
        let (sect, gm, _de, _it) = newton(&r, &mut f, t, 1.0e-12, 300);
        sectors = sect;
        gmax = gm;
        let ratio = t * (sect[0] + sect[1]) / (sect[2] + sect[3]);
        history.push((t, ratio));
        if (ratio - target).abs() < tol {
            break;
        }
        t *= target / ratio;
    }
    (t, RadialProfile { r, f, sectors, t, gmax }, history)
}

/// Signed 1-D degree K = -(2/pi) INT sin^2 f f' dr (midpoint quadrature).
#[must_use]
pub fn degree(r: &[f64], f: &[f64]) -> f64 {
    let nc = r.len() - 1;
    let mut s = 0.0_f64;
    for c in 0..nc {
        let fm = 0.5 * (f[c + 1] + f[c]);
        let sn = det::sin(fm);
        s += sn * sn * (f[c + 1] - f[c]);
    }
    -(2.0 / core::f64::consts::PI) * s
}

/// Yukawa tail fit result.
pub struct TailFit {
    /// Fitted decay constant (compare 1/sqrt(2t)).
    pub mu: f64,
    /// Fitted amplitude A.
    pub a: f64,
    /// RMS of the log-space residual over the fit window.
    pub rms: f64,
    /// Number of nodes in the window.
    pub npts: usize,
    /// (r_first, r_last) of the window.
    pub window: (f64, f64),
}

/// Fit f ~ A (1/(mu r) + 1/(mu r)^2) e^{-mu r} on the window
/// f in (1e-6, 1e-2), r > R* (log-space least squares; 120-step golden
/// section over mu, analytic solve for log A) — radial_solve.py `tail_fit`.
#[must_use]
pub fn tail_fit(r: &[f64], f: &[f64], mu_lo: f64, mu_hi: f64) -> TailFit {
    let rstar = crate::rstar();
    let mut rw = Vec::new();
    let mut y = Vec::new();
    for (&rv, &fv) in r.iter().zip(f.iter()) {
        if fv > 1.0e-6 && fv < 1.0e-2 && rv > rstar {
            rw.push(rv);
            y.push(det::ln(fv));
        }
    }
    let misfit = |mu: f64| -> (f64, f64) {
        let m = rw.len() as f64;
        let mut sum = 0.0_f64;
        let mut model = vec![0.0_f64; rw.len()];
        for (i, &rv) in rw.iter().enumerate() {
            let x = mu * rv;
            model[i] = -x + det::ln(1.0 / x + 1.0 / (x * x));
            sum += y[i] - model[i];
        }
        let lg_a = sum / m;
        let mut ss = 0.0_f64;
        for (i, &mv) in model.iter().enumerate() {
            let res = y[i] - mv - lg_a;
            ss += res * res;
        }
        ((ss / m).sqrt(), lg_a)
    };
    let invphi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let (mut a, mut b) = (mu_lo, mu_hi);
    let mut c = b - invphi * (b - a);
    let mut d = a + invphi * (b - a);
    let (mut fc, _) = misfit(c);
    let (mut fd, _) = misfit(d);
    for _ in 0..120 {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - invphi * (b - a);
            fc = misfit(c).0;
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + invphi * (b - a);
            fd = misfit(d).0;
        }
    }
    let mu = 0.5 * (a + b);
    let (rms, lg_a) = misfit(mu);
    TailFit {
        mu,
        a: det::exp(lg_a),
        rms,
        npts: rw.len(),
        window: (rw[0], rw[rw.len() - 1]),
    }
}

/// True if f is monotone non-increasing (tolerance 1e-12, Python check).
#[must_use]
pub fn monotone_decreasing(f: &[f64]) -> bool {
    f.windows(2).all(|w| w[1] - w[0] <= 1.0e-12)
}
