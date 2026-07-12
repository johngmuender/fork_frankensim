//! S3 degree (b_P) on STORED quaternion fields, plus the geometric
//! solid-angle hedgehog charge.
//!
//! Three density schemes:
//!   - `Central2`: the fs-cosserat-pilot b-density, ported operation-for-
//!     operation so the stored-field run reproduces the pilot's certified
//!     digits (B(48) = 0.9830385118051, B(96) = 0.9956367241585) on the
//!     compacton referee;
//!   - `Central4`: field3d_solve.py's `_derivs_central4` (4th-order
//!     central differences, densities at cells);
//!   - `Corner`: field3d_solve.py's `_derivs_corner` midpoint scheme —
//!     D_i q = edge difference of face averages at the (N+1)^3 cell
//!     corners, q at the corner = 8-cell average; width-1 stencils leave
//!     no sub-grid null modes (robust, not the most accurate quadrature —
//!     see field3d_RESULTS.md section G).
//!
//! Convention (stated, identical to the pilot): b = -(1/2 pi^2)
//! det(q, dq/dx, dq/dy, dq/dz), sign fixed so the f: pi -> 0 hedgehog
//! has degree +1. `degree_b_density` also returns the b field and
//! E6 = pi^3 INT b^2 (free consistency check against the pilot's D3).

use crate::field::{Bc, Field3, QField};
use crate::{det4, solid_angle_origin, v3};
use core::f64::consts::PI;

/// Finite-difference scheme for the b-density.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scheme {
    /// O(h^2) central differences (the pilot's scheme).
    Central2,
    /// O(h^4) central differences (field3d_solve.py `central4`).
    Central4,
    /// Corner/midpoint scheme at (N+1)^3 corners (field3d_solve.py
    /// `corner`).
    Corner,
}

/// Output of [`degree_b_density`].
pub struct DegreeOut {
    /// deg = INT b d^3x.
    pub deg: f64,
    /// E6 = pi^3 INT b^2 d^3x (the sextic-sector consistency number).
    pub e6: f64,
    /// The b-density field (at cells for the central schemes, at the
    /// (N+1)^3 corners for `Corner`).
    pub b: Field3<1>,
}

/// Topological degree of a stored quaternion field via the b-density.
#[must_use]
pub fn degree_b_density(q: &QField, scheme: Scheme) -> DegreeOut {
    match scheme {
        Scheme::Central2 => central(q, false),
        Scheme::Central4 => central(q, true),
        Scheme::Corner => corner(q),
    }
}

/// Central-difference schemes, densities at the stored nodes.
/// `fourth = false` reproduces the pilot's arithmetic exactly:
/// dq = (q(+1) - q(-1)) * inv2h, b = -det4/(2 pi^2), accumulated in
/// fixed (i, j, k) order and scaled by h^3 once at the end.
fn central(q: &QField, fourth: bool) -> DegreeOut {
    let [n0, n1, n2] = q.n;
    let h = q.h;
    let inv2h = 1.0 / (2.0 * h);
    let c1 = 8.0 / (12.0 * h);
    let c2 = 1.0 / (12.0 * h);
    let two_pi2 = 2.0 * PI * PI;
    let mut bdata = Vec::with_capacity(n0 * n1 * n2);
    let (mut s_deg, mut s_b2) = (0.0_f64, 0.0_f64);
    for i in 0..n0 {
        let ii = i as isize;
        for j in 0..n1 {
            let jj = j as isize;
            for k in 0..n2 {
                let kk = k as isize;
                let qc = q.at(i, j, k);
                let mut dx = [0.0; 4];
                let mut dy = [0.0; 4];
                let mut dz = [0.0; 4];
                if fourth {
                    let (xp1, xm1) = (q.get(ii + 1, jj, kk), q.get(ii - 1, jj, kk));
                    let (xp2, xm2) = (q.get(ii + 2, jj, kk), q.get(ii - 2, jj, kk));
                    let (yp1, ym1) = (q.get(ii, jj + 1, kk), q.get(ii, jj - 1, kk));
                    let (yp2, ym2) = (q.get(ii, jj + 2, kk), q.get(ii, jj - 2, kk));
                    let (zp1, zm1) = (q.get(ii, jj, kk + 1), q.get(ii, jj, kk - 1));
                    let (zp2, zm2) = (q.get(ii, jj, kk + 2), q.get(ii, jj, kk - 2));
                    for c in 0..4 {
                        dx[c] = (xp1[c] - xm1[c]) * c1 - (xp2[c] - xm2[c]) * c2;
                        dy[c] = (yp1[c] - ym1[c]) * c1 - (yp2[c] - ym2[c]) * c2;
                        dz[c] = (zp1[c] - zm1[c]) * c1 - (zp2[c] - zm2[c]) * c2;
                    }
                } else {
                    let (qxp, qxm) = (q.get(ii + 1, jj, kk), q.get(ii - 1, jj, kk));
                    let (qyp, qym) = (q.get(ii, jj + 1, kk), q.get(ii, jj - 1, kk));
                    let (qzp, qzm) = (q.get(ii, jj, kk + 1), q.get(ii, jj, kk - 1));
                    for c in 0..4 {
                        dx[c] = (qxp[c] - qxm[c]) * inv2h;
                        dy[c] = (qyp[c] - qym[c]) * inv2h;
                        dz[c] = (qzp[c] - qzm[c]) * inv2h;
                    }
                }
                let b = -det4(&qc, &dx, &dy, &dz) / two_pi2;
                s_deg += b;
                s_b2 += b * b;
                bdata.push([b]);
            }
        }
    }
    let vol = h * h * h;
    DegreeOut {
        deg: s_deg * vol,
        e6: PI.powi(3) * s_b2 * vol,
        b: Field3::from_parts(q.n, h, q.origin, Bc::Vacuum([0.0]), bdata),
    }
}

/// Corner (midpoint) scheme: densities at the (N+1)^3 cell corners.
/// For corner (ic, jc, kc) between cell layers ic-1 and ic (etc.):
///   F_x(i) = quarter-sum of q over the yz-face cells,
///   D_x q  = (F_x(ic) - F_x(ic-1)) / h,   qe = (F_x(ic) + F_x(ic-1)) / 2,
/// and cyclically — the exact port of `_derivs_corner`.
fn corner(q: &QField) -> DegreeOut {
    let [n0, n1, n2] = q.n;
    let h = q.h;
    let two_pi2 = 2.0 * PI * PI;
    let m = [n0 + 1, n1 + 1, n2 + 1];
    let mut bdata = Vec::with_capacity(m[0] * m[1] * m[2]);
    let (mut s_deg, mut s_b2) = (0.0_f64, 0.0_f64);
    // Face average of q over the 4 cells at fixed layer `a` of axis
    // `ax`, spanned by (b0-1..b0, c0-1..c0) on the other two axes.
    let face = |ax: usize, a: isize, b0: isize, c0: isize| -> [f64; 4] {
        let mut s = [0.0; 4];
        for db in -1..=0 {
            for dc in -1..=0 {
                let (i, j, k) = match ax {
                    0 => (a, b0 + db, c0 + dc),
                    1 => (b0 + db, a, c0 + dc),
                    _ => (b0 + db, c0 + dc, a),
                };
                let v = q.get(i, j, k);
                for c in 0..4 {
                    s[c] += v[c];
                }
            }
        }
        for c in &mut s {
            *c *= 0.25;
        }
        s
    };
    for ic in 0..m[0] {
        let i = ic as isize;
        for jc in 0..m[1] {
            let j = jc as isize;
            for kc in 0..m[2] {
                let k = kc as isize;
                let fx1 = face(0, i, j, k);
                let fx0 = face(0, i - 1, j, k);
                let fy1 = face(1, j, i, k);
                let fy0 = face(1, j - 1, i, k);
                let fz1 = face(2, k, i, j);
                let fz0 = face(2, k - 1, i, j);
                let mut qe = [0.0; 4];
                let mut dx = [0.0; 4];
                let mut dy = [0.0; 4];
                let mut dz = [0.0; 4];
                for c in 0..4 {
                    qe[c] = 0.5 * (fx1[c] + fx0[c]);
                    dx[c] = (fx1[c] - fx0[c]) / h;
                    dy[c] = (fy1[c] - fy0[c]) / h;
                    dz[c] = (fz1[c] - fz0[c]) / h;
                }
                let b = -det4(&qe, &dx, &dy, &dz) / two_pi2;
                s_deg += b;
                s_b2 += b * b;
                bdata.push([b]);
            }
        }
    }
    let vol = h * h * h;
    let origin = [
        q.origin[0] - 0.5 * h,
        q.origin[1] - 0.5 * h,
        q.origin[2] - 0.5 * h,
    ];
    DegreeOut {
        deg: s_deg * vol,
        e6: PI.powi(3) * s_b2 * vol,
        b: Field3::from_parts(m, h, origin, Bc::Vacuum([0.0]), bdata),
    }
}

/// Geometric (VOS-kernel, integer-snapping) hedgehog charge of a
/// director field through a closed oriented triangulated surface:
/// sample n at the triangle corners (trilinear + renormalize), sum the
/// signed spherical-triangle solid angles, Q = sum / 4 pi. Lattice-
/// geometric and integer-snapping — the complement of the O(h^2)
/// density route.
#[must_use]
pub fn hedgehog_charge_geometric(n: &Field3<3>, tris: &[[[f64; 3]; 3]]) -> f64 {
    let mut omega = 0.0;
    for t in tris {
        let a = v3::normalize(n.sample_trilinear(t[0]));
        let b = v3::normalize(n.sample_trilinear(t[1]));
        let c = v3::normalize(n.sample_trilinear(t[2]));
        omega += solid_angle_origin(a, b, c);
    }
    omega / (4.0 * PI)
}

/// Outward-oriented triangulation of the axis-aligned cube shell of
/// half-width `s` around `center`, `m x m` quads per face (12 m^2
/// triangles total) — the standard closed probe surface around a
/// suspected point defect.
#[must_use]
pub fn cube_shell(center: [f64; 3], s: f64, m: usize) -> Vec<[[f64; 3]; 3]> {
    let mut tris = Vec::with_capacity(12 * m * m);
    for ax in 0..3usize {
        let bx = (ax + 1) % 3;
        let cx = (ax + 2) % 3;
        for side in [-1.0_f64, 1.0] {
            let pt = |u: usize, w: usize| -> [f64; 3] {
                let mut p = center;
                p[ax] += side * s;
                p[bx] += -s + 2.0 * s * u as f64 / m as f64;
                p[cx] += -s + 2.0 * s * w as f64 / m as f64;
                p
            };
            for u in 0..m {
                for w in 0..m {
                    let (p00, p10, p11, p01) = (pt(u, w), pt(u + 1, w), pt(u + 1, w + 1), pt(u, w + 1));
                    // e_b x e_c = +e_a: (p00,p10,p11) is outward on the
                    // +side face; swap winding on the -side face.
                    if side > 0.0 {
                        tris.push([p00, p10, p11]);
                        tris.push([p00, p11, p01]);
                    } else {
                        tris.push([p00, p11, p10]);
                        tris.push([p00, p01, p11]);
                    }
                }
            }
        }
    }
    tris
}
