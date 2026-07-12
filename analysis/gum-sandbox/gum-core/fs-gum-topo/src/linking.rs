//! Hopf invariant, Method B: preimage linking (the integer referee).
//!
//! Pipeline (survey build_spec, Algorithm 2B):
//!   1. Pick two regular values y1 (a z-hat tilted by fixed
//!      irrational-looking angles — dodges lattice-symmetric
//!      non-transversality) and y2 = -y1.
//!   2. Extract n^{-1}(y_a) by marching tetrahedra: each grid cube is
//!      split into the 6 Freudenthal/Kuhn tets (face-diagonal choices
//!      consistent across neighbouring cubes); in each tet the two
//!      frame coordinates u = n.e1, v = n.e2 about y_a are linearly
//!      interpolated, so their common zero is a unique straight segment
//!      (the reason for tets over cubes); the hemisphere check
//!      w = n.y_a > 0 discards the antipodal sheet. Segments are
//!      oriented by grad u x grad v (the pullback of the target area
//!      form, e1 x e2 = y_a) and chained into closed polylines by their
//!      EXACT shared-face keys (sorted global node triples), so chaining
//!      needs no floating-point matching.
//!   3. Lk(L1, L2) by the exact polyline Gauss-linking sum: per segment
//!      pair, the signed spherical area of the Gauss-map quadrilateral,
//!      evaluated as two van Oosterom-Strackee triangle solid angles —
//!      no quadrature error, only the O(h^2) geometric error in the
//!      curve positions, and the total is an exact integer up to
//!      floating point. H = SIGN_LINK * sum, sign fixed once against
//!      the analytic Hopf-1 referee.
//!
//! Pitfalls handled (survey list): non-regular values are detected as
//! interpolation degeneracies / edge-grazing crossings / odd crossing
//! counts and answered by a deterministic re-tilt-and-retry schedule;
//! curves that fail to close (torus wrap or boundary exit) return
//! [`TopoError::OpenCurve`] — the same obstruction Method A rejects as
//! net flux.

use crate::field::{DirectorField, Field3};
use crate::{solid_angle_origin, v3, TopoError, SIGN_LINK};
use core::f64::consts::PI;
use std::collections::HashMap;

/// Result of the preimage-linking route.
#[derive(Debug, Clone)]
pub struct LinkingOut {
    /// The Hopf invariant (sign-fixed Gauss-linking sum; an integer up
    /// to floating point).
    pub hopf: f64,
    /// The raw linking sum before the sign fix.
    pub raw: f64,
    /// Closed preimage loops found for y1 and y2.
    pub loops: [usize; 2],
    /// Total polyline segments for y1 and y2.
    pub segments: [usize; 2],
    /// Which retry attempt succeeded (0 = first tilt).
    pub attempt: usize,
    /// The regular value y1 actually used.
    pub y1: [f64; 3],
}

/// The 6 Freudenthal/Kuhn tetrahedra of the unit cube: paths from
/// (0,0,0) to (1,1,1) adding one axis step at a time, one per
/// permutation. Face diagonals agree across neighbouring cubes.
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// One extracted preimage segment, endpoints tagged by exact face keys.
struct Seg {
    p: [f64; 3],
    q: [f64; 3],
    f_start: [usize; 3],
    f_end: [usize; 3],
}

enum FaceRes {
    Miss,
    Hit([f64; 3]),
    Bad,
}

/// Extract the closed preimage polylines n^{-1}(y) of a director field
/// by marching tetrahedra over the non-wrapping cubes of the grid.
///
/// # Errors
/// [`TopoError::Degenerate`] on non-transversal geometry (caller
/// re-tilts and retries) and [`TopoError::OpenCurve`] when a curve
/// fails to close inside the marching region.
pub fn extract_preimage_curves(
    f: &DirectorField,
    y: [f64; 3],
) -> Result<Vec<Vec<[f64; 3]>>, TopoError> {
    let [n0, n1, n2] = f.n;
    let e3 = v3::normalize(y);
    let a = if e3[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] };
    let e1 = v3::normalize(v3::cross(a, e3));
    let e2 = v3::cross(e3, e1); // e1 x e2 = e3.

    // Frame coordinates at every node (normalized defensively).
    let total = n0 * n1 * n2;
    let mut u = vec![0.0_f64; total];
    let mut v = vec![0.0_f64; total];
    let mut w = vec![0.0_f64; total];
    for (idx, nv) in f.data.iter().enumerate() {
        let m = v3::normalize(*nv);
        u[idx] = v3::dot(m, e1);
        v[idx] = v3::dot(m, e2);
        w[idx] = v3::dot(m, e3);
    }
    let nid = |i: usize, j: usize, k: usize| (i * n1 + j) * n2 + k;
    let decode = |id: usize| -> (usize, usize, usize) {
        (id / (n1 * n2), (id / n2) % n1, id % n2)
    };
    let node_pos = |id: usize| -> [f64; 3] {
        let (i, j, k) = decode(id);
        f.pos(i, j, k)
    };

    // Zero-crossing of the two linear frame coordinates on a triangular
    // face, keyed by the sorted global node triple.
    let mut cache: HashMap<[usize; 3], FaceRes> = HashMap::new();
    let mut degenerate: Option<String> = None;
    let mut face_hit = |key: [usize; 3],
                        u: &[f64],
                        v: &[f64],
                        w: &[f64],
                        degenerate: &mut Option<String>|
     -> Option<[f64; 3]> {
        let res = cache.entry(key).or_insert_with(|| {
            let (u0, u1, u2) = (u[key[0]], u[key[1]], u[key[2]]);
            let (v0, v1, v2) = (v[key[0]], v[key[1]], v[key[2]]);
            let (du1, du2) = (u1 - u0, u2 - u0);
            let (dv1, dv2) = (v1 - v0, v2 - v0);
            let det = du1 * dv2 - du2 * dv1;
            let scale = du1.abs().max(du2.abs()).max(dv1.abs()).max(dv2.abs());
            if det.abs() < 1.0e-14 * scale * scale {
                // Coordinates nearly collinear on the face: only a
                // problem when zero is actually reachable there.
                let umin = u0.min(u1).min(u2);
                let umax = u0.max(u1).max(u2);
                let vmin = v0.min(v1).min(v2);
                let vmax = v0.max(v1).max(v2);
                if umin <= 0.0 && umax >= 0.0 && vmin <= 0.0 && vmax >= 0.0 {
                    return FaceRes::Bad;
                }
                return FaceRes::Miss;
            }
            let l1 = (-u0 * dv2 + du2 * v0) / det;
            let l2 = (-du1 * v0 + u0 * dv1) / det;
            let l0 = 1.0 - l1 - l2;
            let eps = 1.0e-9;
            if l1.abs() < eps || l2.abs() < eps || l0.abs() < eps {
                return FaceRes::Bad; // Edge/vertex grazing: retry tilt.
            }
            if l1 < 0.0 || l2 < 0.0 || l0 < 0.0 {
                return FaceRes::Miss;
            }
            let wv = w[key[0]] + l1 * (w[key[1]] - w[key[0]]) + l2 * (w[key[2]] - w[key[0]]);
            if wv <= 0.0 {
                return FaceRes::Miss; // Antipodal sheet (n = -y).
            }
            if wv < 0.5 {
                return FaceRes::Bad; // Should be ~ +1 for unit data.
            }
            let (p0, p1, p2) = (node_pos(key[0]), node_pos(key[1]), node_pos(key[2]));
            FaceRes::Hit([
                p0[0] + l1 * (p1[0] - p0[0]) + l2 * (p2[0] - p0[0]),
                p0[1] + l1 * (p1[1] - p0[1]) + l2 * (p2[1] - p0[1]),
                p0[2] + l1 * (p1[2] - p0[2]) + l2 * (p2[2] - p0[2]),
            ])
        });
        match res {
            FaceRes::Miss => None,
            FaceRes::Hit(p) => Some(*p),
            FaceRes::Bad => {
                *degenerate = Some("edge-grazing or non-transversal face crossing".into());
                None
            }
        }
    };

    let mut segs: Vec<Seg> = Vec::new();
    for ci in 0..n0.saturating_sub(1) {
        for cj in 0..n1.saturating_sub(1) {
            for ck in 0..n2.saturating_sub(1) {
                // Cube prefilter: both frame coordinates must change
                // sign somewhere in the cube, on the w > 0 sheet.
                let mut umin = f64::INFINITY;
                let mut umax = f64::NEG_INFINITY;
                let mut vmin = f64::INFINITY;
                let mut vmax = f64::NEG_INFINITY;
                let mut wmax = f64::NEG_INFINITY;
                for di in 0..2 {
                    for dj in 0..2 {
                        for dk in 0..2 {
                            let id = nid(ci + di, cj + dj, ck + dk);
                            umin = umin.min(u[id]);
                            umax = umax.max(u[id]);
                            vmin = vmin.min(v[id]);
                            vmax = vmax.max(v[id]);
                            wmax = wmax.max(w[id]);
                        }
                    }
                }
                if !(umin <= 0.0 && umax >= 0.0 && vmin <= 0.0 && vmax >= 0.0 && wmax > 0.0) {
                    continue;
                }
                for perm in &PERMS {
                    let mut offs = [[0usize; 3]; 4];
                    let mut cur = [0usize; 3];
                    for (s, &ax) in perm.iter().enumerate() {
                        cur[ax] += 1;
                        offs[s + 1] = cur;
                    }
                    let ids: [usize; 4] = core::array::from_fn(|m| {
                        nid(ci + offs[m][0], cj + offs[m][1], ck + offs[m][2])
                    });
                    // Face m omits vertex m; sort for the exact key.
                    let mut hits: Vec<([usize; 3], [f64; 3])> = Vec::new();
                    for m in 0..4 {
                        let mut key = [0usize; 3];
                        let mut t = 0;
                        for (s, &id) in ids.iter().enumerate() {
                            if s != m {
                                key[t] = id;
                                t += 1;
                            }
                        }
                        key.sort_unstable();
                        if let Some(p) = face_hit(key, &u, &v, &w, &mut degenerate) {
                            hits.push((key, p));
                        }
                        if let Some(msg) = degenerate.take() {
                            return Err(TopoError::Degenerate(msg));
                        }
                    }
                    match hits.len() {
                        0 => {}
                        2 => {
                            // Orient by grad u x grad v inside the tet.
                            let p0 = node_pos(ids[0]);
                            let rows: [[f64; 3]; 3] = core::array::from_fn(|r| {
                                v3::sub(node_pos(ids[r + 1]), p0)
                            });
                            let du: [f64; 3] =
                                core::array::from_fn(|r| u[ids[r + 1]] - u[ids[0]]);
                            let dv: [f64; 3] =
                                core::array::from_fn(|r| v[ids[r + 1]] - v[ids[0]]);
                            let detm = crate::det3(rows[0], rows[1], rows[2]);
                            if detm.abs() < 1.0e-300 {
                                return Err(TopoError::Degenerate(
                                    "degenerate tetrahedron".into(),
                                ));
                            }
                            let solve = |d: [f64; 3]| -> [f64; 3] {
                                // Cramer on M g = d with M rows = tet edges.
                                let gx = crate::det3(
                                    [d[0], rows[0][1], rows[0][2]],
                                    [d[1], rows[1][1], rows[1][2]],
                                    [d[2], rows[2][1], rows[2][2]],
                                ) / detm;
                                let gy = crate::det3(
                                    [rows[0][0], d[0], rows[0][2]],
                                    [rows[1][0], d[1], rows[1][2]],
                                    [rows[2][0], d[2], rows[2][2]],
                                ) / detm;
                                let gz = crate::det3(
                                    [rows[0][0], rows[0][1], d[0]],
                                    [rows[1][0], rows[1][1], d[1]],
                                    [rows[2][0], rows[2][1], d[2]],
                                ) / detm;
                                [gx, gy, gz]
                            };
                            let gu = solve(du);
                            let gv = solve(dv);
                            let tdir = v3::cross(gu, gv);
                            let d = v3::dot(v3::sub(hits[1].1, hits[0].1), tdir);
                            if d == 0.0 {
                                return Err(TopoError::Degenerate(
                                    "zero-length or tangent-degenerate segment".into(),
                                ));
                            }
                            let (s0, s1) = if d > 0.0 { (0, 1) } else { (1, 0) };
                            segs.push(Seg {
                                p: hits[s0].1,
                                q: hits[s1].1,
                                f_start: hits[s0].0,
                                f_end: hits[s1].0,
                            });
                        }
                        _ => {
                            return Err(TopoError::Degenerate(format!(
                                "{} face crossings in one tetrahedron",
                                hits.len()
                            )));
                        }
                    }
                }
            }
        }
    }

    // Chain by exact face keys: each face is the end of exactly one
    // segment and the start of exactly one.
    let mut start_map: HashMap<[usize; 3], usize> = HashMap::new();
    for (idx, s) in segs.iter().enumerate() {
        if start_map.insert(s.f_start, idx).is_some() {
            return Err(TopoError::Degenerate("face starts two segments".into()));
        }
    }
    let mut visited = vec![false; segs.len()];
    let mut loops: Vec<Vec<[f64; 3]>> = Vec::new();
    for s0 in 0..segs.len() {
        if visited[s0] {
            continue;
        }
        let mut poly: Vec<[f64; 3]> = Vec::new();
        let mut cur = s0;
        loop {
            visited[cur] = true;
            poly.push(segs[cur].p);
            let Some(&nxt) = start_map.get(&segs[cur].f_end) else {
                return Err(TopoError::OpenCurve);
            };
            if nxt == s0 {
                break;
            }
            if visited[nxt] {
                return Err(TopoError::Degenerate("chain re-entered a segment".into()));
            }
            cur = nxt;
        }
        loops.push(poly);
    }
    Ok(loops)
}

/// Signed spherical area of the Gauss-map quadrilateral of one segment
/// pair (a1 -> a2, b1 -> b2): corners (b1-a1, b1-a2, b2-a2, b2-a1) in
/// (s, t) orientation, split into two VOS triangles.
#[inline]
fn quad_omega(a1: [f64; 3], a2: [f64; 3], b1: [f64; 3], b2: [f64; 3]) -> f64 {
    let r1 = v3::sub(b1, a1);
    let r2 = v3::sub(b1, a2);
    let r3 = v3::sub(b2, a2);
    let r4 = v3::sub(b2, a1);
    solid_angle_origin(r1, r2, r3) + solid_angle_origin(r1, r3, r4)
}

/// Exact Gauss-linking number of two closed polylines (an integer up to
/// floating point for disjoint curves).
#[must_use]
pub fn gauss_linking(la: &[[f64; 3]], lb: &[[f64; 3]]) -> f64 {
    let ma = la.len();
    let mb = lb.len();
    let mut s = 0.0;
    for i in 0..ma {
        let a1 = la[i];
        let a2 = la[(i + 1) % ma];
        for j in 0..mb {
            let b1 = lb[j];
            let b2 = lb[(j + 1) % mb];
            s += quad_omega(a1, a2, b1, b2);
        }
    }
    s / (4.0 * PI)
}

/// Deterministic regular-value schedule: a z-hat tilted by fixed
/// irrational-looking angles, re-tilted on retry `t`.
fn regular_value(t: usize) -> [f64; 3] {
    let alpha = 1.9106332362490188 + 0.0533 * t as f64;
    let beta = 0.9086036519234123 + 0.0771 * t as f64;
    [alpha.sin() * beta.cos(), alpha.sin() * beta.sin(), alpha.cos()]
}

/// The preimage-linking Hopf invariant of a director field (Method B,
/// the integer referee): H = SIGN_LINK * sum of Gauss-linking numbers
/// between the n^{-1}(y1) and n^{-1}(y2) loop families, y2 = -y1, with
/// a deterministic re-tilt-and-retry schedule (up to 5 attempts) for
/// non-regular values.
///
/// # Errors
/// Propagates [`TopoError::OpenCurve`] (torus wrap / boundary exit) and
/// returns the last [`TopoError::Degenerate`] if every retry fails.
pub fn hopf_preimage_linking(f: &Field3<3>) -> Result<LinkingOut, TopoError> {
    let mut last: Option<TopoError> = None;
    for t in 0..5 {
        let y1 = regular_value(t);
        let y2 = [-y1[0], -y1[1], -y1[2]];
        let c1 = match extract_preimage_curves(f, y1) {
            Ok(c) => c,
            Err(e @ TopoError::Degenerate(_)) => {
                last = Some(e);
                continue;
            }
            Err(e) => return Err(e),
        };
        let c2 = match extract_preimage_curves(f, y2) {
            Ok(c) => c,
            Err(e @ TopoError::Degenerate(_)) => {
                last = Some(e);
                continue;
            }
            Err(e) => return Err(e),
        };
        let mut raw = 0.0;
        for la in &c1 {
            for lb in &c2 {
                raw += gauss_linking(la, lb);
            }
        }
        return Ok(LinkingOut {
            hopf: SIGN_LINK * raw,
            raw,
            loops: [c1.len(), c2.len()],
            segments: [
                c1.iter().map(Vec::len).sum(),
                c2.iter().map(Vec::len).sum(),
            ],
            attempt: t,
            y1,
        });
    }
    Err(last.unwrap_or_else(|| TopoError::Degenerate("no regular value found".into())))
}
