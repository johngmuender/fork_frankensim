//! Objective + analytic gradient on the anisotropic field — the exact
//! generalization of `fs_gum_statics::engine::eval` to (nx, ny, nz), with
//! one deliberate extension for the two-charge sector: the Bogomolny-floor
//! normalization degree `deg_unit` is carried SEPARATELY from the anchor
//! reference `deg_ref`.
//!
//! In the single-knot protocol fgap = E6 + E0 - FLOOR * deg / deg_ref uses
//! deg_ref (the engine hedgehog degree, ~1) as "one unit of degree".  For a
//! degree-2 configuration the same formula with deg_ref ~ 2 would halve the
//! per-degree floor; here `fgap = E6 + E0 - FLOOR * deg / deg_unit` with
//! `deg_unit = deg_ref / B` restores the intended per-unit-degree wall
//! (B = 2 for the pair, B = 1 reproduces fs-gum-statics verbatim — asserted
//! bitwise by the twin gate).  Everything else (PAIRS expansion, guard
//! forms, adjoints, projection, accumulation order) is the frozen field3d
//! protocol, transcribed with anisotropic bounds.

use fs_gum_field::{bps_floor, Scheme};

use core::f64::consts::PI;

use crate::field_a::{aidx_a, corner_a, derivs4_a, FieldA};

const FOURPI: f64 = 4.0 * PI;

/// Determinant sign convention (fs-gum-statics SGN6).
pub const SGN6: f64 = -1.0;

/// 6-term pair expansion of det[q, Dx q, Dy q, Dz q] (field3d PAIRS).
const PAIRS: [(usize, usize, usize, usize, f64); 6] = [
    (0, 1, 2, 3, 1.0),
    (0, 2, 1, 3, -1.0),
    (0, 3, 1, 2, 1.0),
    (1, 2, 0, 3, 1.0),
    (1, 3, 0, 2, -1.0),
    (2, 3, 0, 1, 1.0),
];

/// One-sided degree anchor (fs-gum-statics `Anchor`).
#[derive(Clone, Copy, Debug)]
pub struct AnchorA {
    pub mu: f64,
    pub band: f64,
}

/// One-sided Bogomolny-floor wall (fs-gum-statics `Wall`).
#[derive(Clone, Copy, Debug)]
pub struct WallA {
    pub mu: f64,
    pub fgap_ref: f64,
}

/// Objective specification (fs-gum-statics `Opts` + `deg_unit`).
#[derive(Clone, Copy, Debug)]
pub struct OptsA {
    pub scheme: Scheme,
    pub t: f64,
    pub c: [f64; 4],
    /// Degree-anchor reference (the seed's own degree in this protocol).
    pub deg_ref: f64,
    /// Floor-normalization degree ("one unit"); = deg_ref / B.
    pub deg_unit: f64,
    pub anchor: Option<AnchorA>,
    pub wall: Option<WallA>,
}

impl OptsA {
    /// Bare E_static objective at frozen t (no guards).
    #[must_use]
    pub fn estatic(scheme: Scheme) -> Self {
        let t = fs_gum_field::T_FROZEN;
        OptsA {
            scheme,
            t,
            c: [t, t, 1.0, 1.0],
            deg_ref: 1.0,
            deg_unit: 1.0,
            anchor: None,
            wall: None,
        }
    }

    /// Campaign guards at frozen strengths, referenced to the seed's own
    /// degree and floor gap (`deg_unit = deg_ref / b_charge`).
    #[must_use]
    pub fn with_guards(mut self, deg_ref: f64, deg_unit: f64, fgap_ref: f64) -> Self {
        self.deg_ref = deg_ref;
        self.deg_unit = deg_unit;
        self.anchor = Some(AnchorA { mu: fs_gum_statics::MU_DEG, band: fs_gum_statics::DEG_BAND });
        self.wall = Some(WallA { mu: fs_gum_statics::MU_FLOOR, fgap_ref });
        self
    }
}

/// One evaluation (fs-gum-statics `Out`, minus the Routhian which the
/// statics protocol here never uses).
#[derive(Clone, Copy, Debug)]
pub struct OutA {
    pub e2: f64,
    pub e4: f64,
    pub e6: f64,
    pub e0: f64,
    pub i: f64,
    pub deg: f64,
    pub estat: f64,
    pub floor_gap: f64,
    pub epen: f64,
    pub efpen: f64,
    pub obj: f64,
}

/// Gradient on the REAL cells, anisotropic Grad layout, h^3-weighted and
/// tangent-projected.
pub struct GradA {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub data: Vec<f64>,
}

impl GradA {
    fn zeros(nx: usize, ny: usize, nz: usize) -> Self {
        GradA { nx, ny, nz, data: vec![0.0; 4 * nx * ny * nz] }
    }

    /// Euclidean norm (fixed sequential sum).
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

fn for_each_point<F: FnMut(usize, usize, usize, &[f64; 4], &[[f64; 4]; 3])>(
    f: &FieldA,
    scheme: Scheme,
    mut cb: F,
) {
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
    match scheme {
        Scheme::Central4 => {
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let q = f.get(i, j, k);
                        let d = derivs4_a(f, i, j, k);
                        cb(i, j, k, &q, &d);
                    }
                }
            }
        }
        Scheme::Corner => {
            for ci in 0..=nx {
                for cj in 0..=ny {
                    for ck in 0..=nz {
                        let (qe, d) = corner_a(f, ci, cj, ck);
                        cb(ci, cj, ck, &qe, &d);
                    }
                }
            }
        }
        Scheme::Central2 => {
            panic!("fs-gum-twoknot: gradients are defined for Central4/Corner only")
        }
    }
}

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
    const OTHERS: [(usize, usize, usize, usize); 3] = [(1, 2, 0, 1), (0, 2, 0, 2), (0, 1, 1, 2)];
    let two4pi = 2.0 / FOURPI;
    let mut p = [[0.0_f64; 4]; 3];
    for ax in 0..3 {
        let (o1, o2, k1, k2) = OTHERS[ax];
        let t1 = c2 + c4 * (nn[o1] + nn[o2]);
        for a in 0..4 {
            p[ax][a] = two4pi * (t1 * d[ax][a] - c4 * (dd[k1] * d[o1][a] + dd[k2] * d[o2][a]));
        }
    }
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

/// Evaluate objective (and optionally the exact tangent-projected analytic
/// gradient) — fs-gum-statics `eval`, anisotropic bounds, `deg_unit` floor.
#[must_use]
#[allow(clippy::too_many_lines)]
pub fn eval_a(f: &FieldA, o: &OptsA, need_grad: bool) -> (OutA, Option<GradA>) {
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
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
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let q = f.get(i, j, k);
                s_e0 += 1.0 - q[0];
                s_i += q[1] * q[1] + q[2] * q[2];
            }
        }
    }
    let e2 = s_e2 * h3 / FOURPI;
    let e4 = s_e4 * h3 / FOURPI;
    let e6 = s_det2 * h3 / FOURPI;
    let e0 = s_e0 * h3 / FOURPI;
    let i_val = 2.0 * s_i * h3;
    let deg = SGN6 * s_det * h3 / (2.0 * PI * PI);
    let estat = o.t * (e2 + e4) + e6 + e0;
    let floor_gap = e6 + e0 - bps_floor() * deg / o.deg_unit;
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
    let obj = o.c[0] * e2 + o.c[1] * e4 + o.c[2] * e6 + o.c[3] * e0 + epen + efpen;
    let out = OutA { e2, e4, e6, e0, i: i_val, deg, estat, floor_gap, epen, efpen, obj };
    if !need_grad {
        return (out, None);
    }

    // ---- pass 2: gradient --------------------------------------------------
    let w6s = (o.c[2] + wall_mu * devf) / (2.0 * PI);
    let mut cst = 0.0_f64;
    if dev < 0.0 {
        cst += o.anchor.map_or(0.0, |a| a.mu) * dev;
    }
    if devf < 0.0 {
        cst -= wall_mu * devf * bps_floor() / o.deg_unit;
    }
    let cstw = cst * SGN6 / (2.0 * PI * PI);
    let e0c = -(o.c[3] + wall_mu * devf) / FOURPI;
    let mut g = GradA::zeros(nx, ny, nz);
    match o.scheme {
        Scheme::Central4 => {
            let c1 = 8.0 / (12.0 * h);
            let c2s = 1.0 / (12.0 * h);
            for_each_point(f, o.scheme, |i, j, k, q, d| {
                let (p, gq) = point_flux(q, d, o.c[0], o.c[1], w6s, cstw);
                for a in 0..4 {
                    g.data[aidx_a(nx, ny, nz, a, i, j, k)] += gq[a];
                    // axis x
                    let pv = p[0][a];
                    if i + 1 < nx {
                        g.data[aidx_a(nx, ny, nz, a, i + 1, j, k)] += c1 * pv;
                    }
                    if i >= 1 {
                        g.data[aidx_a(nx, ny, nz, a, i - 1, j, k)] -= c1 * pv;
                    }
                    if i + 2 < nx {
                        g.data[aidx_a(nx, ny, nz, a, i + 2, j, k)] -= c2s * pv;
                    }
                    if i >= 2 {
                        g.data[aidx_a(nx, ny, nz, a, i - 2, j, k)] += c2s * pv;
                    }
                    // axis y
                    let pv = p[1][a];
                    if j + 1 < ny {
                        g.data[aidx_a(nx, ny, nz, a, i, j + 1, k)] += c1 * pv;
                    }
                    if j >= 1 {
                        g.data[aidx_a(nx, ny, nz, a, i, j - 1, k)] -= c1 * pv;
                    }
                    if j + 2 < ny {
                        g.data[aidx_a(nx, ny, nz, a, i, j + 2, k)] -= c2s * pv;
                    }
                    if j >= 2 {
                        g.data[aidx_a(nx, ny, nz, a, i, j - 2, k)] += c2s * pv;
                    }
                    // axis z
                    let pv = p[2][a];
                    if k + 1 < nz {
                        g.data[aidx_a(nx, ny, nz, a, i, j, k + 1)] += c1 * pv;
                    }
                    if k >= 1 {
                        g.data[aidx_a(nx, ny, nz, a, i, j, k - 1)] -= c1 * pv;
                    }
                    if k + 2 < nz {
                        g.data[aidx_a(nx, ny, nz, a, i, j, k + 2)] -= c2s * pv;
                    }
                    if k >= 2 {
                        g.data[aidx_a(nx, ny, nz, a, i, j, k - 2)] += c2s * pv;
                    }
                }
            });
        }
        Scheme::Corner => {
            let ih = 1.0 / h;
            for_each_point(f, o.scheme, |ci, cj, ck, qe, d| {
                let (p, gq) = point_flux(qe, d, o.c[0], o.c[1], w6s, cstw);
                for di in 0..2usize {
                    if ci + di == 0 || ci + di > nx {
                        continue;
                    }
                    let i = ci + di - 1;
                    let sx = if di == 1 { 1.0 } else { -1.0 };
                    for dj in 0..2usize {
                        if cj + dj == 0 || cj + dj > ny {
                            continue;
                        }
                        let j = cj + dj - 1;
                        let sy = if dj == 1 { 1.0 } else { -1.0 };
                        for dk in 0..2usize {
                            if ck + dk == 0 || ck + dk > nz {
                                continue;
                            }
                            let k = ck + dk - 1;
                            let sz = if dk == 1 { 1.0 } else { -1.0 };
                            for a in 0..4 {
                                g.data[aidx_a(nx, ny, nz, a, i, j, k)] += 0.25
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
    // cell-density terms (E0 with the wall's dE0 part)
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                g.data[aidx_a(nx, ny, nz, 0, i, j, k)] += e0c;
            }
        }
    }
    // h^3 weight, then pointwise tangent projection g -= (g.q) q
    for v in g.data.iter_mut() {
        *v *= h3;
    }
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let q = f.get(i, j, k);
                let mut dot = 0.0_f64;
                for a in 0..4 {
                    dot += g.data[aidx_a(nx, ny, nz, a, i, j, k)] * q[a];
                }
                for a in 0..4 {
                    let o2 = aidx_a(nx, ny, nz, a, i, j, k);
                    g.data[o2] -= dot * q[a];
                }
            }
        }
    }
    (out, Some(g))
}
