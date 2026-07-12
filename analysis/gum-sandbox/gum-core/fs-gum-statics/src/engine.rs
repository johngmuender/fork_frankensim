//! Forward objective + analytic variational gradient of the campaign
//! functional on the stored field, for both field3d discretisations
//! (`Scheme::Central4` measurement engine, `Scheme::Corner` descent engine).
//!
//! Port of `field3d_solve.py Engine.energy` (lines 410-641): each vectorized
//! numpy line became a loop nest over evaluation points; the layout-heavy
//! buffer choreography of the Python version collapses into per-point local
//! arrays plus a deterministic scatter through the scheme adjoint.
//! The gradient is the EXACT derivative of the discrete objective (checked
//! against central finite differences at 1e-5 by gate G-A):
//!
//!  * pass 1 (forward): sector sums, penalties, objective;
//!  * pass 2 (gradient): per evaluation point, the (2+4) flux
//!    `P[ax][a] = (2/4pi)[(c2 + c4(nn_o1+nn_o2)) D_ax q_a
//!                        - c4 (D_ax.D_o1) D_o1 q_a - c4 (D_ax.D_o2) D_o2 q_a]`,
//!    the sextic cofactor terms from the 6-term pair expansion with
//!    det-weight `w6 = (c6 + MU_F devf) det/(2pi) + cst sgn6/(2pi^2)`
//!    (`cst = MU dev - MU_F devf FLOOR/deg_ref`, the one-sided guards), and
//!    the density derivative w.r.t. the evaluation-point value `q_E`; then
//!    the scheme adjoint scatters flux and density terms back to the REAL
//!    cells (ghost flux = zero: contributions to/from ghost cells are
//!    dropped, exactly the transpose of reading fixed vacuum ghosts);
//!  * cell terms: E0 density `-(c0 + MU_F devf)/(4pi)` on q0 and the
//!    Routhian term `-(L^2/2I^2) 4 q_{1,2}`;
//!  * finally `g *= h^3` and the pointwise tangent projection
//!    `g -= (g.q) q`.
//!
//! Sector weights `c = [c2, c4, c6, c0]` generalise the frozen objective
//! (`c = [t, t, 1, 1]` = E_static) so that gate G-A can check every sector's
//! gradient separately with the SAME code path.

use fs_gum_field::{bps_floor, stencil, Field3, Scheme, T_FROZEN};

use core::f64::consts::PI;

use crate::diag::aidx;
use crate::{DEG_BAND, MU_DEG, MU_FLOOR};

const FOURPI: f64 = 4.0 * PI;

/// Determinant sign convention (field3d `sgn6`): measured on the hedgehog
/// the raw det[q, Dx q, Dy q, Dz q] integrates to a negative number, so
/// sgn6 = -1 fixes degree = +1 for the f: pi -> 0 hedgehog.  Identical to
/// fs-gum-field's b = -det/(2 pi^2).
pub const SGN6: f64 = -1.0;

/// 6-term pair expansion of det[q, Dx q, Dy q, Dz q] (field3d PAIRS):
/// det = sum_p sg_p (qE_a Dx_b - qE_b Dx_a)(Dy_c Dz_d - Dy_d Dz_c).
const PAIRS: [(usize, usize, usize, usize, f64); 6] = [
    (0, 1, 2, 3, 1.0),
    (0, 2, 1, 3, -1.0),
    (0, 3, 1, 2, 1.0),
    (1, 2, 0, 3, 1.0),
    (1, 3, 0, 2, -1.0),
    (2, 3, 0, 1, 1.0),
];

/// One-sided degree anchor E_pen = (mu/2) min(0, deg - (deg_ref - band))^2.
#[derive(Clone, Copy, Debug)]
pub struct Anchor {
    pub mu: f64,
    pub band: f64,
}

/// One-sided Bogomolny-floor wall E_fpen = (mu/2) min(0, fgap - fgap_ref)^2,
/// fgap = E6 + E0 - (32 sqrt2/15) deg / deg_ref.
#[derive(Clone, Copy, Debug)]
pub struct Wall {
    pub mu: f64,
    pub fgap_ref: f64,
}

/// Objective specification.  The objective differentiated by [`eval`] is
/// `obj = c[0] E2 + c[1] E4 + c[2] E6 + c[3] E0 + L^2/(2I) + E_pen + E_fpen`
/// (rot/penalty terms only when configured).  `t` is carried separately so
/// the reported `estat` is always the bare frozen-dial static energy.
#[derive(Clone, Copy, Debug)]
pub struct Opts {
    pub scheme: Scheme,
    pub t: f64,
    /// Objective weights on (E2, E4, E6, E0); `[t, t, 1, 1]` = E_static.
    pub c: [f64; 4],
    /// Fixed isospin L of the Routhian term L^2/(2I); `None` = statics.
    pub l: Option<f64>,
    /// Degree-anchor / floor-gap reference degree (engine's own hedgehog
    /// degree in the campaign protocol; 1.0 before calibration).
    pub deg_ref: f64,
    pub anchor: Option<Anchor>,
    pub wall: Option<Wall>,
}

impl Opts {
    /// Bare E_static objective at frozen t (no guards, no rotation).
    #[must_use]
    pub fn estatic(scheme: Scheme) -> Self {
        let t = T_FROZEN;
        Opts { scheme, t, c: [t, t, 1.0, 1.0], l: None, deg_ref: 1.0, anchor: None, wall: None }
    }

    /// Bare fixed-L Routhian objective at frozen t.
    #[must_use]
    pub fn routhian(scheme: Scheme, l: f64) -> Self {
        Opts { l: Some(l), ..Opts::estatic(scheme) }
    }

    /// Add the campaign guards with their frozen strengths (MU = 5000,
    /// band = 0.005, MU_F = 400) referenced to the engine's own hedgehog
    /// (`deg_ref` = hedgehog degree, `fgap_ref` = hedgehog gap - 0.01).
    #[must_use]
    pub fn with_guards(mut self, deg_ref: f64, fgap_ref: f64) -> Self {
        self.deg_ref = deg_ref;
        self.anchor = Some(Anchor { mu: MU_DEG, band: DEG_BAND });
        self.wall = Some(Wall { mu: MU_FLOOR, fgap_ref });
        self
    }
}

/// One evaluation of the objective: bare sectors + Routhian + guard
/// penalties.  Reported energies are ALWAYS the bare functionals; the
/// penalties appear only in `epen`/`efpen`/`obj`.
#[derive(Clone, Copy, Debug)]
pub struct Out {
    pub e2: f64,
    pub e4: f64,
    pub e6: f64,
    pub e0: f64,
    pub i: f64,
    pub deg: f64,
    /// t (E2 + E4) + E6 + E0 at `Opts::t` (bare).
    pub estat: f64,
    /// estat + L^2/(2I) (bare Routhian; = estat when L is None).
    pub r: f64,
    /// E6 + E0 - (32 sqrt2/15) deg / deg_ref (drift diagnostic).
    pub floor_gap: f64,
    pub epen: f64,
    pub efpen: f64,
    /// The descent objective: weighted sectors + rotation + penalties.
    pub obj: f64,
}

/// Gradient on the REAL cells, component-major unpadded layout
/// `idx(a, i, j, k) = ((a n + i) n + j) n + k`, h^3-weighted and
/// tangent-projected (the exact convention of the Python `g` buffer).
pub struct Grad {
    pub n: usize,
    pub data: Vec<f64>,
}

impl Grad {
    fn zeros(n: usize) -> Self {
        Grad { n, data: vec![0.0; 4 * n * n * n] }
    }

    /// Flat index of component `a` at real cell (i, j, k).
    #[inline]
    #[must_use]
    pub fn idx(&self, a: usize, i: usize, j: usize, k: usize) -> usize {
        ((a * self.n + i) * self.n + j) * self.n + k
    }

    /// Euclidean norm of the flat gradient vector (fixed sequential sum).
    #[must_use]
    pub fn norm(&self) -> f64 {
        let mut s = 0.0_f64;
        for &v in &self.data {
            s += v * v;
        }
        s.sqrt()
    }
}

#[inline]
fn nn_dd(d: &[[f64; 4]; 3]) -> ([f64; 3], [f64; 3]) {
    let mut nn = [0.0_f64; 3];
    for ax in 0..3 {
        let v = &d[ax];
        nn[ax] = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
    }
    let dot = |i: usize, j: usize| -> f64 {
        d[i][0] * d[j][0] + d[i][1] * d[j][1] + d[i][2] * d[j][2] + d[i][3] * d[j][3]
    };
    (nn, [dot(0, 1), dot(0, 2), dot(1, 2)])
}

/// det[qE, Dx, Dy, Dz] via the 6-term pair expansion, fixed term order
/// (field3d PAIRS accumulation).
#[inline]
fn det_pairs(q: &[f64; 4], d: &[[f64; 4]; 3]) -> f64 {
    let mut det = 0.0_f64;
    for &(a, b, c, e, sg) in &PAIRS {
        let m = q[a] * d[0][b] - q[b] * d[0][a];
        let pc = d[1][c] * d[2][e] - d[1][e] * d[2][c];
        det += sg * (m * pc);
    }
    det
}

/// Visit every derivative-evaluation point of the scheme in fixed ascending
/// order, passing (q at the point, D_i q).  Central4: the N^3 cells.
/// Corner: the (N+1)^3 corners with the 8-cell average.
fn for_each_point<F: FnMut(usize, usize, usize, &[f64; 4], &[[f64; 4]; 3])>(
    f: &Field3,
    scheme: Scheme,
    mut cb: F,
) {
    let n = f.n();
    match scheme {
        Scheme::Central4 => {
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let q = f.get(i, j, k);
                        let d = stencil::derivs4(f, i, j, k);
                        cb(i, j, k, &q, &d);
                    }
                }
            }
        }
        Scheme::Corner => {
            for ci in 0..=n {
                for cj in 0..=n {
                    for ck in 0..=n {
                        let (qe, d) = stencil::corner(f, ci, cj, ck);
                        cb(ci, cj, ck, &qe, &d);
                    }
                }
            }
        }
        Scheme::Central2 => {
            panic!("fs-gum-statics: gradients are defined for Central4/Corner only")
        }
    }
}

/// Per-point flux and density derivative.  Returns
/// (`p[ax][a]` = d(dens)/d(D_ax q_a), `gq[a]` = d(dens)/d(qE_a) sextic part).
#[inline]
fn point_flux(
    q: &[f64; 4],
    d: &[[f64; 4]; 3],
    c2: f64,
    c4: f64,
    w6s: f64,
    cstw: f64,
) -> ([[f64; 4]; 3], [f64; 4]) {
    let (nn, dd) = nn_dd(d);
    // (o1, o2, dd index of (ax,o1), dd index of (ax,o2)) per axis
    const OTHERS: [(usize, usize, usize, usize); 3] = [(1, 2, 0, 1), (0, 2, 0, 2), (0, 1, 1, 2)];
    let two4pi = 2.0 / FOURPI;
    let mut p = [[0.0_f64; 4]; 3];
    for ax in 0..3 {
        let (o1, o2, k1, k2) = OTHERS[ax];
        let t1 = c2 + c4 * (nn[o1] + nn[o2]);
        for a in 0..4 {
            p[ax][a] =
                two4pi * (t1 * d[ax][a] - c4 * (dd[k1] * d[o1][a] + dd[k2] * d[o2][a]));
        }
    }
    // sextic cofactor terms through the det weight
    let det = det_pairs(q, d);
    let w6 = w6s * det + cstw;
    let mut gq = [0.0_f64; 4];
    for &(a, b, c, e, sg) in &PAIRS {
        let m = q[a] * d[0][b] - q[b] * d[0][a];
        let pc = d[1][c] * d[2][e] - d[1][e] * d[2][c];
        let t2 = sg * w6 * pc;
        gq[a] += t2 * d[0][b];
        gq[b] -= t2 * d[0][a];
        p[0][b] += t2 * q[a];
        p[0][a] -= t2 * q[b];
        let t3 = sg * w6 * m;
        p[1][c] += t3 * d[2][e];
        p[1][e] -= t3 * d[2][c];
        p[2][e] += t3 * d[1][c];
        p[2][c] -= t3 * d[1][e];
    }
    (p, gq)
}

/// Evaluate the objective (and, if `need_grad`, its exact tangent-projected
/// analytic gradient) of the stored field under `o`.  See the module docs
/// for the algorithm; every constant and guard form is the frozen field3d
/// protocol.
#[must_use]
pub fn eval(f: &Field3, o: &Opts, need_grad: bool) -> (Out, Option<Grad>) {
    let n = f.n();
    let h = f.h();
    let h3 = h * h * h;

    // ---- pass 1: forward sums -------------------------------------------
    let (mut s_e2, mut s_e4, mut s_det, mut s_det2) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for_each_point(f, o.scheme, |_i, _j, _k, q, d| {
        let (nn, dd) = nn_dd(d);
        s_e2 += nn[0] + nn[1] + nn[2];
        s_e4 += (nn[0] * nn[1] - dd[0] * dd[0])
            + (nn[0] * nn[2] - dd[1] * dd[1])
            + (nn[1] * nn[2] - dd[2] * dd[2]);
        let det = det_pairs(q, d);
        s_det += det;
        s_det2 += det * det;
    });
    let (mut s_e0, mut s_i) = (0.0_f64, 0.0_f64);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let q = f.get(i, j, k);
                s_e0 += 1.0 - q[0];
                s_i += q[1] * q[1] + q[2] * q[2];
            }
        }
    }
    let e2 = s_e2 * h3 / FOURPI;
    let e4 = s_e4 * h3 / FOURPI;
    let e6 = s_det2 * h3 / FOURPI; // = pi^3 INT b^2, b = -det/(2 pi^2)
    let e0 = s_e0 * h3 / FOURPI;
    let i_val = 2.0 * s_i * h3;
    let deg = SGN6 * s_det * h3 / (2.0 * PI * PI);
    let estat = o.t * (e2 + e4) + e6 + e0;
    let rot = match o.l {
        Some(l) => l * l / (2.0 * i_val),
        None => 0.0,
    };
    let r = estat + rot;
    let floor_gap = e6 + e0 - bps_floor() * deg / o.deg_ref;
    let (dev, epen) = match o.anchor {
        Some(a) => {
            let dev = (deg - (o.deg_ref - a.band)).min(0.0);
            (dev, 0.5 * a.mu * dev * dev)
        }
        None => (0.0, 0.0),
    };
    let (devf, efpen, wall_mu) = match o.wall {
        Some(w) => {
            let devf = (floor_gap - w.fgap_ref).min(0.0);
            (devf, 0.5 * w.mu * devf * devf, w.mu)
        }
        None => (0.0, 0.0, 0.0),
    };
    let obj = o.c[0] * e2 + o.c[1] * e4 + o.c[2] * e6 + o.c[3] * e0 + rot + epen + efpen;
    let out =
        Out { e2, e4, e6, e0, i: i_val, deg, estat, r, floor_gap, epen, efpen, obj };
    if !need_grad {
        return (out, None);
    }

    // ---- pass 2: gradient ------------------------------------------------
    // guard-dressed weights (field3d energy() gradient path, verbatim):
    //   d obj / d E6  = c6 + MU_F devf          -> det weight w6s * det
    //   d obj / d E0  = c0 + MU_F devf          -> e0c on q0
    //   d obj / d deg = MU dev - MU_F devf FLOOR/deg_ref  -> cstw on det
    let w6s = (o.c[2] + wall_mu * devf) / (2.0 * PI);
    let mut cst = 0.0_f64;
    if dev < 0.0 {
        // anchor is Some here by construction (dev != 0 only when active)
        cst += o.anchor.map_or(0.0, |a| a.mu) * dev;
    }
    if devf < 0.0 {
        cst -= wall_mu * devf * bps_floor() / o.deg_ref;
    }
    let cstw = cst * SGN6 / (2.0 * PI * PI);
    let e0c = -(o.c[3] + wall_mu * devf) / FOURPI;
    let rotfac = match o.l {
        Some(l) => -(l * l) / (2.0 * i_val * i_val) * 4.0,
        None => 0.0,
    };
    let mut g = Grad::zeros(n);
    match o.scheme {
        Scheme::Central4 => {
            // eval points = cells; adjoint of the 4th-order stencil with
            // zero ghost flux (contributions to/from ghost cells dropped).
            let c1 = 8.0 / (12.0 * h);
            let c2s = 1.0 / (12.0 * h);
            for_each_point(f, o.scheme, |i, j, k, q, d| {
                let (p, gq) = point_flux(q, d, o.c[0], o.c[1], w6s, cstw);
                for a in 0..4 {
                    g.data[aidx(n, a, i, j, k)] += gq[a];
                    // axis x
                    let pv = p[0][a];
                    if i + 1 < n {
                        g.data[aidx(n, a, i + 1, j, k)] += c1 * pv;
                    }
                    if i >= 1 {
                        g.data[aidx(n, a, i - 1, j, k)] -= c1 * pv;
                    }
                    if i + 2 < n {
                        g.data[aidx(n, a, i + 2, j, k)] -= c2s * pv;
                    }
                    if i >= 2 {
                        g.data[aidx(n, a, i - 2, j, k)] += c2s * pv;
                    }
                    // axis y
                    let pv = p[1][a];
                    if j + 1 < n {
                        g.data[aidx(n, a, i, j + 1, k)] += c1 * pv;
                    }
                    if j >= 1 {
                        g.data[aidx(n, a, i, j - 1, k)] -= c1 * pv;
                    }
                    if j + 2 < n {
                        g.data[aidx(n, a, i, j + 2, k)] -= c2s * pv;
                    }
                    if j >= 2 {
                        g.data[aidx(n, a, i, j - 2, k)] += c2s * pv;
                    }
                    // axis z
                    let pv = p[2][a];
                    if k + 1 < n {
                        g.data[aidx(n, a, i, j, k + 1)] += c1 * pv;
                    }
                    if k >= 1 {
                        g.data[aidx(n, a, i, j, k - 1)] -= c1 * pv;
                    }
                    if k + 2 < n {
                        g.data[aidx(n, a, i, j, k + 2)] -= c2s * pv;
                    }
                    if k >= 2 {
                        g.data[aidx(n, a, i, j, k - 2)] += c2s * pv;
                    }
                }
            });
        }
        Scheme::Corner => {
            // adjoints of the corner operators: corner (ci,cj,ck) touches
            // the 8 cells (ci-1+di, cj-1+dj, ck-1+dk); ghosts receive
            // nothing.  dD_x/dcell = (1/h) sx/4, dqe/dcell = 1/8 with
            // sx = +1 on the di = 1 face, -1 on di = 0 (and cyclically).
            let ih = 1.0 / h;
            for_each_point(f, o.scheme, |ci, cj, ck, qe, d| {
                let (p, gq) = point_flux(qe, d, o.c[0], o.c[1], w6s, cstw);
                for di in 0..2usize {
                    if ci + di == 0 || ci + di > n {
                        continue;
                    }
                    let i = ci + di - 1;
                    let sx = if di == 1 { 1.0 } else { -1.0 };
                    for dj in 0..2usize {
                        if cj + dj == 0 || cj + dj > n {
                            continue;
                        }
                        let j = cj + dj - 1;
                        let sy = if dj == 1 { 1.0 } else { -1.0 };
                        for dk in 0..2usize {
                            if ck + dk == 0 || ck + dk > n {
                                continue;
                            }
                            let k = ck + dk - 1;
                            let sz = if dk == 1 { 1.0 } else { -1.0 };
                            for a in 0..4 {
                                g.data[aidx(n, a, i, j, k)] += 0.25
                                    * ih
                                    * (sx * p[0][a] + sy * p[1][a] + sz * p[2][a])
                                    + 0.125 * gq[a];
                            }
                        }
                    }
                }
            });
        }
        Scheme::Central2 => unreachable!(),
    }
    // cell-density terms (E0 with the wall's dE0 part; Routhian dI part)
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                g.data[aidx(n, 0, i, j, k)] += e0c;
                if o.l.is_some() {
                    let q = f.get(i, j, k);
                    g.data[aidx(n, 1, i, j, k)] += rotfac * q[1];
                    g.data[aidx(n, 2, i, j, k)] += rotfac * q[2];
                }
            }
        }
    }
    // h^3 weight, then pointwise tangent projection g -= (g.q) q
    for v in g.data.iter_mut() {
        *v *= h3;
    }
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let q = f.get(i, j, k);
                let mut dot = 0.0_f64;
                for a in 0..4 {
                    dot += g.data[aidx(n, a, i, j, k)] * q[a];
                }
                for a in 0..4 {
                    let o = g.idx(a, i, j, k);
                    g.data[o] -= dot * q[a];
                }
            }
        }
    }
    (out, Some(g))
}
