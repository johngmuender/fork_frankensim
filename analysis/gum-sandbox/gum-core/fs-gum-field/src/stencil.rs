//! Finite-difference current stencils on STORED data.
//!
//! Three engines, all reading through the 2-layer vacuum ghost rind:
//!  * `derivs2` — 2nd-order central differences at cell centres
//!    (fs-cosserat-pilot's stencil, now on stored data);
//!  * `derivs4` — 4th-order central differences at cell centres, the
//!    field3d `central4` measurement scheme, coefficients c1 = 8/(12h),
//!    c2 = 1/(12h);
//!  * `corner` — the compact corner (midpoint) scheme at the (N+1)^3 cell
//!    corners: D_i q = edge difference of the perpendicular-face 4-cell
//!    average, q at the corner = 8-cell average (field3d `_derivs_corner`,
//!    with the same nested-average evaluation order).

use crate::field::{Field3, GHOST};

/// 2nd-order central derivatives D_i q at real cell (i, j, k):
/// D_x q = (q[i+1] - q[i-1]) / (2h), reading ghosts at the rim.
/// Returns `d[axis][component]`.
#[must_use]
pub fn derivs2(f: &Field3, i: usize, j: usize, k: usize) -> [[f64; 4]; 3] {
    let inv2h = 1.0 / (2.0 * f.h());
    let (ip, jp, kp) = (i + GHOST, j + GHOST, k + GHOST);
    let mut d = [[0.0; 4]; 3];
    for c in 0..4 {
        d[0][c] = (f.get_pad(c, ip + 1, jp, kp) - f.get_pad(c, ip - 1, jp, kp)) * inv2h;
        d[1][c] = (f.get_pad(c, ip, jp + 1, kp) - f.get_pad(c, ip, jp - 1, kp)) * inv2h;
        d[2][c] = (f.get_pad(c, ip, jp, kp + 1) - f.get_pad(c, ip, jp, kp - 1)) * inv2h;
    }
    d
}

/// 4th-order central derivatives D_i q at real cell (i, j, k):
/// D = c1 (q[+1] - q[-1]) - c2 (q[+2] - q[-2]), c1 = 8/(12h), c2 = 1/(12h)
/// (field3d `_derivs_central4`, same expression order).
#[must_use]
pub fn derivs4(f: &Field3, i: usize, j: usize, k: usize) -> [[f64; 4]; 3] {
    let c1 = 8.0 / (12.0 * f.h());
    let c2 = 1.0 / (12.0 * f.h());
    let (ip, jp, kp) = (i + GHOST, j + GHOST, k + GHOST);
    let mut d = [[0.0; 4]; 3];
    for c in 0..4 {
        d[0][c] = (f.get_pad(c, ip + 1, jp, kp) - f.get_pad(c, ip - 1, jp, kp)) * c1
            - (f.get_pad(c, ip + 2, jp, kp) - f.get_pad(c, ip - 2, jp, kp)) * c2;
        d[1][c] = (f.get_pad(c, ip, jp + 1, kp) - f.get_pad(c, ip, jp - 1, kp)) * c1
            - (f.get_pad(c, ip, jp + 2, kp) - f.get_pad(c, ip, jp - 2, kp)) * c2;
        d[2][c] = (f.get_pad(c, ip, jp, kp + 1) - f.get_pad(c, ip, jp, kp - 1)) * c1
            - (f.get_pad(c, ip, jp, kp + 2) - f.get_pad(c, ip, jp, kp - 2)) * c2;
    }
    d
}

/// Compact corner (midpoint) scheme at corner (ci, cj, ck), indices
/// 0..=n per axis.  Corner (0,0,0) touches one layer of ghost cells (the
/// scheme needs a 1-layer rind; our 2-layer rind covers it).
///
/// Returns `(qe, d)` where `qe` is the 8-cell average of q at the corner
/// and `d[axis][component]` the edge-difference derivatives.  The nested
/// 0.5*(a+b) averaging order matches field3d `_derivs_corner`.
#[must_use]
pub fn corner(f: &Field3, ci: usize, cj: usize, ck: usize) -> ([f64; 4], [[f64; 4]; 3]) {
    let ih = 1.0 / f.h();
    // padded cell indices of the 8 cells around this corner: the corner
    // sits between real cells (c-1) and (c) per axis, i.e. padded
    // (c + GHOST - 1) and (c + GHOST).
    let (i0, j0, k0) = (ci + GHOST - 1, cj + GHOST - 1, ck + GHOST - 1);
    let mut qe = [0.0; 4];
    let mut d = [[0.0; 4]; 3];
    for a in 0..4 {
        let v = |di: usize, dj: usize, dk: usize| f.get_pad(a, i0 + di, j0 + dj, k0 + dk);
        // branch ay -> ayz -> {D_x, qe}
        let ay = |di: usize, dk: usize| 0.5 * (v(di, 1, dk) + v(di, 0, dk));
        let ayz = |di: usize| 0.5 * (ay(di, 1) + ay(di, 0));
        d[0][a] = (ayz(1) - ayz(0)) * ih;
        qe[a] = 0.5 * (ayz(1) + ayz(0));
        // branch ax -> {axz -> D_y, axy -> D_z}
        let ax = |dj: usize, dk: usize| 0.5 * (v(1, dj, dk) + v(0, dj, dk));
        let axz = |dj: usize| 0.5 * (ax(dj, 1) + ax(dj, 0));
        d[1][a] = (axz(1) - axz(0)) * ih;
        let axy = |dk: usize| 0.5 * (ax(1, dk) + ax(0, dk));
        d[2][a] = (axy(1) - axy(0)) * ih;
    }
    (qe, d)
}
