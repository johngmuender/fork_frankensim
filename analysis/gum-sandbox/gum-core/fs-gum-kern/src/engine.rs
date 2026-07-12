//! The gum sweep kernels re-hosted on the tiled deterministic layer:
//! forward sector sums as tiled reductions, the gradient pass with the
//! adjoint INVERTED from scatter to gather (two parallel-map passes), and
//! renormalize / vector ops / norms as tiled sweeps.
//!
//! Drop-in contract: [`eval`] takes and returns the SAME types as
//! `fs_gum_statics::engine::eval` (`Opts`, `Out`, `Grad`) plus a thread
//! count, and computes the SAME discrete objective and its exact analytic
//! gradient.  The per-point physics (`nn_dd`, `det_pairs`, `point_flux`,
//! the PAIRS table, guard dressing) is copied VERBATIM from
//! `fs-gum-statics/src/engine.rs` so per-point values are bit-identical;
//! what changes — deliberately, once — is the ACCUMULATION ORDER:
//!
//!  * reductions run per-tile (ascending (i,j,k) within the i-slab) and
//!    combine up the fixed pairwise tree — vs the old global flat sum;
//!  * the adjoint runs as a GATHER: each real cell collects its 13
//!    (central4: 4 flux terms per axis + its own gq) or 8 (corner)
//!    contributions in the pinned order stated in
//!    [`crate::GUM_KERN_BIT_SEMANTICS`] — vs the old data-race-prone
//!    scatter whose per-cell arrival order was the global loop order.
//!
//! Both changes are the documented Phase-G1 golden bump (see RESULTS.md);
//! the gates binary measures the old-vs-new deltas once at machine-eps
//! class and freezes NEW bit goldens for the tiled order.
//!
//! GATHER STENCIL DERIVATION (verified against the scatter by gate K-A
//! and the in-crate test): the central4 scatter at eval point i sends
//! `+c1 p` to i+1, `-c1 p` to i-1, `-c2 p` to i+2, `+c2 p` to i-2
//! (ghost flux dropped), so cell i GATHERS `+c1 P[i-1] - c1 P[i+1]
//! - c2 P[i-2] + c2 P[i+2]` per axis, out-of-range eval points dropped,
//! plus its own gq.  The corner scatter sends
//! `0.25/h (sx p0 + sy p1 + sz p2) + 0.125 gq` from corner c to the 8
//! cells `c - 1 + d` (d in {0,1}^3, sx = +1 on the d_x = 1 face), so cell
//! i gathers from its 8 surrounding corners `i + e` (e in {0,1}^3) with
//! `sx = +1` for `e_x = 0` (the cell sits on the upper face of that
//! corner) and `-1` for `e_x = 1`; no boundary dropping arises because
//! all 8 corners of a real cell exist.

use fs_gum_field::{bps_floor, stencil, Field3, Scheme, GHOST};
use fs_gum_statics::engine::{Grad, Opts, Out, SGN6};

use core::f64::consts::PI;
use std::ops::Range;

use crate::tile::{split_tiles, sweep_reduce, tile_ranges, tiled_run, MaxP, SumP};
use crate::reduce::Reduce;

const FOURPI: f64 = 4.0 * PI;

/// 6-term pair expansion table (fs-gum-statics engine.rs, verbatim).
const PAIRS: [(usize, usize, usize, usize, f64); 6] = [
    (0, 1, 2, 3, 1.0),
    (0, 2, 1, 3, -1.0),
    (0, 3, 1, 2, 1.0),
    (1, 2, 0, 3, 1.0),
    (1, 3, 0, 2, -1.0),
    (2, 3, 0, 1, 1.0),
];

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
/// (fs-gum-statics engine.rs, verbatim).
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

/// Per-point flux and density derivative (fs-gum-statics engine.rs,
/// verbatim): (`p[ax][a]` = d(dens)/d(D_ax q_a), `gq[a]` = sextic
/// d(dens)/d(qE_a)).
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

// ---------------------------------------------------------------------------
// tiled forward reductions
// ---------------------------------------------------------------------------

/// Forward stencil-sum partial (s_e2, s_e4, s_det, s_det2).
#[derive(Clone, Copy)]
struct FwdP {
    e2: f64,
    e4: f64,
    det: f64,
    det2: f64,
}

impl Reduce for FwdP {
    fn identity() -> Self {
        FwdP { e2: 0.0, e4: 0.0, det: 0.0, det2: 0.0 }
    }
    fn merge(self, o: Self) -> Self {
        FwdP {
            e2: self.e2 + o.e2,
            e4: self.e4 + o.e4,
            det: self.det + o.det,
            det2: self.det2 + o.det2,
        }
    }
}

/// Cell-density partial (s_e0, s_i).
#[derive(Clone, Copy)]
struct CellP {
    e0: f64,
    i: f64,
}

impl Reduce for CellP {
    fn identity() -> Self {
        CellP { e0: 0.0, i: 0.0 }
    }
    fn merge(self, o: Self) -> Self {
        CellP { e0: self.e0 + o.e0, i: self.i + o.i }
    }
}

/// One forward-accumulation step — the exact statements of the old
/// engine's pass-1 closure.
#[inline]
fn fwd_point(p: &mut FwdP, q: &[f64; 4], d: &[[f64; 4]; 3]) {
    let (nn, dd) = nn_dd(d);
    p.e2 += nn[0] + nn[1] + nn[2];
    p.e4 += (nn[0] * nn[1] - dd[0] * dd[0])
        + (nn[0] * nn[2] - dd[1] * dd[1])
        + (nn[1] * nn[2] - dd[2] * dd[2]);
    let det = det_pairs(q, d);
    p.det += det;
    p.det2 += det * det;
}

fn forward_sums(f: &Field3, scheme: Scheme, threads: usize) -> FwdP {
    let n = f.n();
    match scheme {
        Scheme::Central4 => sweep_reduce(n, threads, |r: Range<usize>| {
            let mut p = FwdP::identity();
            for i in r {
                for j in 0..n {
                    for k in 0..n {
                        let q = f.get(i, j, k);
                        let d = stencil::derivs4(f, i, j, k);
                        fwd_point(&mut p, &q, &d);
                    }
                }
            }
            p
        }),
        Scheme::Corner => sweep_reduce(n + 1, threads, |r: Range<usize>| {
            let mut p = FwdP::identity();
            for ci in r {
                for cj in 0..=n {
                    for ck in 0..=n {
                        let (qe, d) = stencil::corner(f, ci, cj, ck);
                        fwd_point(&mut p, &qe, &d);
                    }
                }
            }
            p
        }),
        Scheme::Central2 => {
            panic!("fs-gum-kern: gradients are defined for Central4/Corner only")
        }
    }
}

fn cell_sums(f: &Field3, threads: usize) -> CellP {
    let n = f.n();
    sweep_reduce(n, threads, |r: Range<usize>| {
        let mut p = CellP::identity();
        for i in r {
            for j in 0..n {
                for k in 0..n {
                    let q = f.get(i, j, k);
                    p.e0 += 1.0 - q[0];
                    p.i += q[1] * q[1] + q[2] * q[2];
                }
            }
        }
        p
    })
}

// ---------------------------------------------------------------------------
// eval: forward + gather-form gradient
// ---------------------------------------------------------------------------

/// Reusable scratch for [`eval_ws`] — the pass-A flux record buffer
/// (16 f64 per eval point: P[3][4] then gq[4]; ~117 MB at N = 96 corner).
/// Reusing it across ANF iterations avoids a large per-eval allocation.
#[derive(Default)]
pub struct Ws {
    flux: Vec<f64>,
}

impl Ws {
    #[must_use]
    pub fn new() -> Self {
        Ws::default()
    }
}

/// Flux-record layout: 16 doubles per eval point, point-major
/// (k fastest); slots 0..11 = P[ax][a] (ax-major), 12..15 = gq[a].
#[inline]
fn flux_base(i: usize, j: usize, k: usize, mj: usize, mk: usize) -> usize {
    ((i * mj + j) * mk + k) * 16
}

/// Evaluate the objective (and, if `need_grad`, its exact tangent-
/// projected analytic gradient) on `threads` threads through the tiled
/// deterministic layer.  Bit-identical for every `threads >= 1` by
/// construction; allocation-free on the flux path when `ws` is reused.
#[must_use]
pub fn eval_ws(
    f: &Field3,
    o: &Opts,
    need_grad: bool,
    threads: usize,
    ws: &mut Ws,
) -> (Out, Option<Grad>) {
    let n = f.n();
    let h = f.h();
    let h3 = h * h * h;

    // ---- pass 1: forward sums (tiled reductions) --------------------------
    let fwd = forward_sums(f, o.scheme, threads);
    let cells = cell_sums(f, threads);

    // ---- Out assembly (fs-gum-statics engine.rs, verbatim) ----------------
    let e2 = fwd.e2 * h3 / FOURPI;
    let e4 = fwd.e4 * h3 / FOURPI;
    let e6 = fwd.det2 * h3 / FOURPI;
    let e0 = cells.e0 * h3 / FOURPI;
    let i_val = 2.0 * cells.i * h3;
    let deg = SGN6 * fwd.det * h3 / (2.0 * PI * PI);
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
    let out = Out { e2, e4, e6, e0, i: i_val, deg, estat, r, floor_gap, epen, efpen, obj };
    if !need_grad {
        return (out, None);
    }

    // ---- guard-dressed weights (verbatim) ----------------------------------
    let w6s = (o.c[2] + wall_mu * devf) / (2.0 * PI);
    let mut cst = 0.0_f64;
    if dev < 0.0 {
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

    // ---- pass A: per-point flux records (parallel map, disjoint chunks) ---
    let mi = match o.scheme {
        Scheme::Central4 => n,
        Scheme::Corner => n + 1,
        Scheme::Central2 => unreachable!(),
    };
    let (mj, mk) = (mi, mi);
    let len = 16 * mi * mj * mk;
    if ws.flux.len() < len {
        ws.flux.resize(len, 0.0);
    }
    let flux_all = &mut ws.flux[..len];
    {
        let bounds = tile_ranges(mi);
        let jobs: Vec<Vec<&mut [f64]>> = split_tiles(flux_all, 1, mj * mk * 16, &bounds);
        let scheme = o.scheme;
        let (c2, c4) = (o.c[0], o.c[1]);
        tiled_run(jobs, threads, |t, mut job| {
            let chunk = job.pop().expect("one component");
            let r = bounds[t].clone();
            let mut off = 0usize;
            for i in r {
                for j in 0..mj {
                    for k in 0..mk {
                        let (q, d) = match scheme {
                            Scheme::Central4 => (f.get(i, j, k), stencil::derivs4(f, i, j, k)),
                            Scheme::Corner => stencil::corner(f, i, j, k),
                            Scheme::Central2 => unreachable!(),
                        };
                        let (p, gq) = point_flux(&q, &d, c2, c4, w6s, cstw);
                        for (ax, pax) in p.iter().enumerate() {
                            chunk[off + 4 * ax..off + 4 * ax + 4].copy_from_slice(pax);
                        }
                        chunk[off + 12..off + 16].copy_from_slice(&gq);
                        off += 16;
                    }
                }
            }
        });
    }
    let flux: &[f64] = flux_all;

    // ---- pass B: gather adjoint + cell terms + h^3 + projection -----------
    // (parallel map over output-cell tiles, disjoint writes; per-cell term
    // order is pinned — see GUM_KERN_BIT_SEMANTICS.)
    let mut g = Grad { n, data: vec![0.0_f64; 4 * n * n * n] };
    let bounds = tile_ranges(n);
    let jobs = split_tiles(&mut g.data, 4, n * n, &bounds);
    let rotate = o.l.is_some();
    let finish = |q: &[f64; 4], ga: &mut [f64; 4]| {
        ga[0] += e0c;
        if rotate {
            ga[1] += rotfac * q[1];
            ga[2] += rotfac * q[2];
        }
        for v in ga.iter_mut() {
            *v *= h3;
        }
        let mut dot = 0.0_f64;
        for a in 0..4 {
            dot += ga[a] * q[a];
        }
        for a in 0..4 {
            ga[a] -= dot * q[a];
        }
    };
    match o.scheme {
        Scheme::Central4 => {
            let c1 = 8.0 / (12.0 * h);
            let c2s = 1.0 / (12.0 * h);
            tiled_run(jobs, threads, |t, mut job| {
                let r = bounds[t].clone();
                let i0 = r.start;
                for i in r {
                    for j in 0..n {
                        for k in 0..n {
                            let q = f.get(i, j, k);
                            let mut ga = [0.0_f64; 4];
                            for (a, gv) in ga.iter_mut().enumerate() {
                                let mut acc = 0.0_f64;
                                // axis x (slot 0): eval points i-1, i+1, i-2, i+2
                                if i >= 1 {
                                    acc += c1 * flux[flux_base(i - 1, j, k, mj, mk) + a];
                                }
                                if i + 1 < n {
                                    acc -= c1 * flux[flux_base(i + 1, j, k, mj, mk) + a];
                                }
                                if i >= 2 {
                                    acc -= c2s * flux[flux_base(i - 2, j, k, mj, mk) + a];
                                }
                                if i + 2 < n {
                                    acc += c2s * flux[flux_base(i + 2, j, k, mj, mk) + a];
                                }
                                // axis y (slot 1)
                                if j >= 1 {
                                    acc += c1 * flux[flux_base(i, j - 1, k, mj, mk) + 4 + a];
                                }
                                if j + 1 < n {
                                    acc -= c1 * flux[flux_base(i, j + 1, k, mj, mk) + 4 + a];
                                }
                                if j >= 2 {
                                    acc -= c2s * flux[flux_base(i, j - 2, k, mj, mk) + 4 + a];
                                }
                                if j + 2 < n {
                                    acc += c2s * flux[flux_base(i, j + 2, k, mj, mk) + 4 + a];
                                }
                                // axis z (slot 2)
                                if k >= 1 {
                                    acc += c1 * flux[flux_base(i, j, k - 1, mj, mk) + 8 + a];
                                }
                                if k + 1 < n {
                                    acc -= c1 * flux[flux_base(i, j, k + 1, mj, mk) + 8 + a];
                                }
                                if k >= 2 {
                                    acc -= c2s * flux[flux_base(i, j, k - 2, mj, mk) + 8 + a];
                                }
                                if k + 2 < n {
                                    acc += c2s * flux[flux_base(i, j, k + 2, mj, mk) + 8 + a];
                                }
                                // own gq, last
                                acc += flux[flux_base(i, j, k, mj, mk) + 12 + a];
                                *gv = acc;
                            }
                            finish(&q, &mut ga);
                            let local = ((i - i0) * n + j) * n + k;
                            for (a, chunk) in job.iter_mut().enumerate() {
                                chunk[local] = ga[a];
                            }
                        }
                    }
                }
            });
        }
        Scheme::Corner => {
            let ih = 1.0 / h;
            tiled_run(jobs, threads, |t, mut job| {
                let r = bounds[t].clone();
                let i0 = r.start;
                for i in r {
                    for j in 0..n {
                        for k in 0..n {
                            let q = f.get(i, j, k);
                            let mut ga = [0.0_f64; 4];
                            for (a, gv) in ga.iter_mut().enumerate() {
                                let mut acc = 0.0_f64;
                                // the 8 surrounding corners, ascending
                                // (ei, ej, ek); sx = +1 on the ei = 0
                                // corner (the cell is its upper-x face)
                                for ei in 0..2usize {
                                    let sx = if ei == 0 { 1.0 } else { -1.0 };
                                    for ej in 0..2usize {
                                        let sy = if ej == 0 { 1.0 } else { -1.0 };
                                        for ek in 0..2usize {
                                            let sz = if ek == 0 { 1.0 } else { -1.0 };
                                            let ob =
                                                flux_base(i + ei, j + ej, k + ek, mj, mk);
                                            acc += 0.25
                                                * ih
                                                * (sx * flux[ob + a]
                                                    + sy * flux[ob + 4 + a]
                                                    + sz * flux[ob + 8 + a])
                                                + 0.125 * flux[ob + 12 + a];
                                        }
                                    }
                                }
                                *gv = acc;
                            }
                            finish(&q, &mut ga);
                            let local = ((i - i0) * n + j) * n + k;
                            for (a, chunk) in job.iter_mut().enumerate() {
                                chunk[local] = ga[a];
                            }
                        }
                    }
                }
            });
        }
        Scheme::Central2 => unreachable!(),
    }
    (out, Some(g))
}

/// [`eval_ws`] with a throwaway workspace.
#[must_use]
pub fn eval(f: &Field3, o: &Opts, need_grad: bool, threads: usize) -> (Out, Option<Grad>) {
    eval_ws(f, o, need_grad, threads, &mut Ws::new())
}

// ---------------------------------------------------------------------------
// parallel-map sweeps: renormalize, vector ops, norms
// ---------------------------------------------------------------------------

/// Pointwise unit-norm renormalisation of all real cells as a tiled
/// parallel map (disjoint writes; per-cell arithmetic identical to
/// `Field3::renormalize`, so the updated FIELD is bit-identical to the
/// serial original at every thread count).  Returns max |nrm - 1| before
/// the fix via a tiled max-reduction (order-insensitive over a fixed
/// value set, so this too matches the old serial value bitwise).
pub fn renormalize(f: &mut Field3, threads: usize) -> f64 {
    let n = f.n();
    let p = f.padded();
    let p2 = p * p;
    // padded plane ranges of the real i-slabs
    let bounds: Vec<Range<usize>> =
        tile_ranges(n).into_iter().map(|r| r.start + GHOST..r.end + GHOST).collect();
    let jobs = split_tiles(f.data_mut(), 4, p2, &bounds);
    let drifts = tiled_run(jobs, threads, |t, mut job| {
        let r = bounds[t].clone();
        let ip0 = r.start;
        let mut drift = 0.0_f64;
        for ip in r {
            for j in 0..n {
                for k in 0..n {
                    let local = ((ip - ip0) * p + (j + GHOST)) * p + (k + GHOST);
                    let q = [job[0][local], job[1][local], job[2][local], job[3][local]];
                    let nrm =
                        (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
                    let d = (nrm - 1.0).abs();
                    if d > drift {
                        drift = d;
                    }
                    for (a, chunk) in job.iter_mut().enumerate() {
                        chunk[local] = q[a] / nrm;
                    }
                }
            }
        }
        MaxP(drift)
    });
    crate::reduce::pairwise_fold(drifts).0
}

/// Euclidean norm of a flat auxiliary field (Grad layout, `4 n^3`) as a
/// tiled reduction: per-tile plain sums in pinned (a, i, j, k) ascending
/// order, combined up the fixed tree, sqrt at the end.  This is the NEW
/// canonical g.norm order (the old one was a flat a-major running sum —
/// part of the documented golden bump).
#[must_use]
pub fn norm_flat(data: &[f64], n: usize, threads: usize) -> f64 {
    let n2 = n * n;
    let n3 = n2 * n;
    let s = sweep_reduce(n, threads, |r: Range<usize>| {
        let mut s = 0.0_f64;
        for a in 0..4 {
            for i in r.clone() {
                for &v in &data[a * n3 + i * n2..a * n3 + (i + 1) * n2] {
                    s += v * v;
                }
            }
        }
        SumP(s)
    });
    s.0.sqrt()
}

/// `y[i] -= a * x[i]` over equal-length flat buffers as a parallel map
/// (elementwise, so bit-identical at any thread count and any tiling; the
/// expression matches the ANF velocity update `*vm -= a * gm`).
pub fn axpy_sub(y: &mut [f64], a: f64, x: &[f64], threads: usize) {
    let unit = 4096usize; // bit-neutral chunk grain (pure elementwise map)
    let m = y.len().div_ceil(unit);
    let bounds: Vec<Range<usize>> = tile_ranges(m)
        .into_iter()
        .map(|r| r.start * unit..(r.end * unit).min(y.len()))
        .collect();
    let mut jobs: Vec<(usize, &mut [f64])> = Vec::with_capacity(bounds.len());
    let mut rest = y;
    let mut pos = 0usize;
    for b in &bounds {
        let (chunk, tail) = rest.split_at_mut(b.end - b.start);
        jobs.push((pos, chunk));
        pos = b.end;
        rest = tail;
    }
    tiled_run(jobs, threads, |_t, (start, chunk)| {
        for (o, ym) in chunk.iter_mut().enumerate() {
            *ym -= a * x[start + o];
        }
    });
}

/// `q += dt * v` on all real cells as a parallel map (the ANF position
/// step; per-cell arithmetic identical to the old get/set loop).
pub fn step_q(f: &mut Field3, v: &[f64], dt: f64, threads: usize) {
    let n = f.n();
    let p = f.padded();
    let p2 = p * p;
    let n2 = n * n;
    let n3 = n2 * n;
    let bounds: Vec<Range<usize>> =
        tile_ranges(n).into_iter().map(|r| r.start + GHOST..r.end + GHOST).collect();
    let jobs = split_tiles(f.data_mut(), 4, p2, &bounds);
    tiled_run(jobs, threads, |t, mut job| {
        let r = bounds[t].clone();
        let ip0 = r.start;
        for ip in r {
            let i = ip - GHOST;
            for j in 0..n {
                for k in 0..n {
                    let local = ((ip - ip0) * p + (j + GHOST)) * p + (k + GHOST);
                    for (a, chunk) in job.iter_mut().enumerate() {
                        chunk[local] += dt * v[a * n3 + i * n2 + j * n + k];
                    }
                }
            }
        }
    });
}

/// Pointwise tangent re-projection of a flat auxiliary field:
/// `u -= (u . q) q` per cell, as a parallel map (the ANF velocity
/// re-projection; per-cell arithmetic identical to the old loop).
pub fn tangent_project(f: &Field3, u: &mut [f64], threads: usize) {
    let n = f.n();
    let n2 = n * n;
    let bounds = tile_ranges(n);
    let jobs = split_tiles(u, 4, n2, &bounds);
    tiled_run(jobs, threads, |t, mut job| {
        let r = bounds[t].clone();
        let i0 = r.start;
        for i in r {
            for j in 0..n {
                for k in 0..n {
                    let q = f.get(i, j, k);
                    let local = ((i - i0) * n + j) * n + k;
                    let mut dot = 0.0_f64;
                    for (a, &qv) in q.iter().enumerate() {
                        dot += job[a][local] * qv;
                    }
                    for (a, &qv) in q.iter().enumerate() {
                        job[a][local] -= dot * qv;
                    }
                }
            }
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_gum_statics::diag::{add_scaled, bump_field};
    use fs_gum_statics::Opts as SOpts;

    const N: usize = 12;
    const LBOX: f64 = 4.5;

    fn test_field() -> Field3 {
        let mut f = Field3::new(N, LBOX);
        f.sample_compacton();
        let bump = bump_field(N, LBOX, 77, 4, 0.05);
        add_scaled(&mut f, &bump, 1.0);
        let _ = f.renormalize();
        f
    }

    fn test_opts(scheme: Scheme) -> Opts {
        // rotor + both guards forced active: exercises every gradient term
        SOpts::routhian(scheme, 3.0).with_guards(1.0, 10.0)
    }

    #[test]
    fn matches_the_old_engine_at_machine_eps_class() {
        let f = test_field();
        for scheme in [Scheme::Central4, Scheme::Corner] {
            let o = test_opts(scheme);
            let (old, gold) = fs_gum_statics::eval(&f, &o, true);
            let (new, gnew) = eval(&f, &o, true, 2);
            for (a, b) in [
                (old.e2, new.e2),
                (old.e4, new.e4),
                (old.e6, new.e6),
                (old.e0, new.e0),
                (old.i, new.i),
                (old.deg, new.deg),
                (old.obj, new.obj),
            ] {
                assert!((a / b - 1.0).abs() < 1.0e-11, "{scheme:?}: {a} vs {b}");
            }
            let (gold, gnew) = (gold.unwrap(), gnew.unwrap());
            let gmax = gold.data.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
            let dmax = gold
                .data
                .iter()
                .zip(gnew.data.iter())
                .fold(0.0_f64, |m, (&a, &b)| m.max((a - b).abs()));
            assert!(dmax / gmax < 1.0e-10, "{scheme:?}: gather-vs-scatter {dmax} / {gmax}");
        }
    }

    #[test]
    fn eval_is_bit_identical_across_thread_counts() {
        let f = test_field();
        for scheme in [Scheme::Central4, Scheme::Corner] {
            let o = test_opts(scheme);
            let (out1, g1) = eval(&f, &o, true, 1);
            let g1 = g1.unwrap();
            for threads in 2..=4 {
                let (outt, gt) = eval(&f, &o, true, threads);
                let gt = gt.unwrap();
                for (a, b) in [
                    (out1.obj, outt.obj),
                    (out1.deg, outt.deg),
                    (out1.estat, outt.estat),
                    (out1.i, outt.i),
                ] {
                    assert_eq!(a.to_bits(), b.to_bits(), "{scheme:?} t={threads}");
                }
                assert!(
                    g1.data.iter().zip(gt.data.iter()).all(|(a, b)| a.to_bits() == b.to_bits()),
                    "{scheme:?} grad t={threads}"
                );
            }
        }
    }

    #[test]
    fn map_sweeps_match_their_serial_references_bitwise() {
        let f0 = test_field();
        // renormalize: new tiled vs Field3::renormalize
        let mut fa = fs_gum_statics::diag::clone_field(&f0);
        let mut fb = fs_gum_statics::diag::clone_field(&f0);
        let da = fa.renormalize();
        let db = renormalize(&mut fb, 3);
        assert_eq!(da.to_bits(), db.to_bits());
        assert!(fa.data().iter().zip(fb.data().iter()).all(|(a, b)| a.to_bits() == b.to_bits()));
        // axpy / step_q / tangent_project vs the old ANF loops
        let n = N;
        let n3 = n * n * n;
        let g = bump_field(n, LBOX, 5, 3, 1.0);
        let mut v1: Vec<f64> = bump_field(n, LBOX, 6, 3, 0.5);
        let mut v2 = v1.clone();
        for (vm, &gm) in v1.iter_mut().zip(g.iter()) {
            *vm -= 0.37 * gm;
        }
        axpy_sub(&mut v2, 0.37, &g, 3);
        assert!(v1.iter().zip(v2.iter()).all(|(a, b)| a.to_bits() == b.to_bits()));
        let mut fc = fs_gum_statics::diag::clone_field(&f0);
        let mut fd = fs_gum_statics::diag::clone_field(&f0);
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let mut q = fc.get(i, j, k);
                    for (c, qv) in q.iter_mut().enumerate() {
                        *qv += 0.01 * v1[fs_gum_statics::diag::aidx(n, c, i, j, k)];
                    }
                    fc.set(i, j, k, q);
                }
            }
        }
        step_q(&mut fd, &v1, 0.01, 3);
        assert!(fc.data().iter().zip(fd.data().iter()).all(|(a, b)| a.to_bits() == b.to_bits()));
        let mut u1 = v1.clone();
        let mut u2 = v1.clone();
        fs_gum_statics::diag::tangent_project(&f0, &mut u1);
        tangent_project(&f0, &mut u2, 3);
        assert!(u1.iter().zip(u2.iter()).all(|(a, b)| a.to_bits() == b.to_bits()));
        // norm: same value class as the old flat sum (order differs)
        let old = {
            let mut s = 0.0_f64;
            for &x in &g[..4 * n3] {
                s += x * x;
            }
            s.sqrt()
        };
        let newn = norm_flat(&g, n, 3);
        assert!((old / newn - 1.0).abs() < 1.0e-13);
        for threads in 2..=4 {
            assert_eq!(norm_flat(&g, n, 1).to_bits(), norm_flat(&g, n, threads).to_bits());
        }
    }
}
