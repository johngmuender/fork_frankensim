//! Disclination diagnostics for RP2 (line-field) directors:
//! pi1(RP2) = Z/2 loop holonomy by sign-gauge transport, half-integer
//! planar winding for framed loops, and dual-plaquette defect-line
//! tracing with the lines-cannot-end structural check.
//!
//! The stored data are ordinary 3-component vectors; the diagnostics
//! never assume a global vector lift — sign-gauge transport IS the RP2
//! semantics (survey build_spec, Algorithm 3).

use crate::field::Field3;
use crate::{v3, TopoError};
use core::f64::consts::PI;
use std::collections::HashMap;

/// Default core tolerance: transport aborts when |n_a . n_b| drops
/// below this (the loop grazes a defect core / the grid is too coarse).
pub const CORE_TOL: f64 = 0.1;

/// The two elements of pi1(RP2) = Z/2.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Z2 {
    /// Contractible holonomy class.
    Trivial,
    /// The nontrivial class (the loop encircles a disclination).
    NonTrivial,
}

/// A closed lattice loop given by its grid-node vertices in order
/// (the closing edge back to `verts[0]` is implicit).
#[derive(Clone, Debug)]
pub struct GridLoop {
    /// Node indices (i, j, k) in traversal order.
    pub verts: Vec<[usize; 3]>,
}

impl GridLoop {
    /// Axis-aligned rectangle loop at fixed k, counter-clockwise seen
    /// from +z: (i0..i1, j0) -> (i1, j0..j1) -> (i1..i0, j1) ->
    /// (i0, j1..j0). Requires i1 > i0 and j1 > j0.
    #[must_use]
    pub fn rectangle_xy(i0: usize, i1: usize, j0: usize, j1: usize, k: usize) -> GridLoop {
        let mut verts = Vec::new();
        for i in i0..i1 {
            verts.push([i, j0, k]);
        }
        for j in j0..j1 {
            verts.push([i1, j, k]);
        }
        for i in (i0 + 1..=i1).rev() {
            verts.push([i, j1, k]);
        }
        for j in (j0 + 1..=j1).rev() {
            verts.push([i0, j, k]);
        }
        GridLoop { verts }
    }
}

/// Z/2 loop holonomy by sign-gauge transport: lift the first director,
/// then n_{s+1} := sign(n_{s+1} . n_s) n_{s+1}; the class is nontrivial
/// iff the lift closes to -n_0. The only frame-independent invariant
/// for a general 3-D loop.
///
/// # Errors
/// [`TopoError::CoreHit`] when |n_{s+1} . n_s| < `tol` at some step.
pub fn z2_holonomy(n: &Field3<3>, lp: &GridLoop, tol: f64) -> Result<Z2, TopoError> {
    let read = |v: [usize; 3]| v3::normalize(n.at(v[0], v[1], v[2]));
    let mut cur = read(lp.verts[0]);
    let mut sign = 1.0_f64;
    for s in 1..=lp.verts.len() {
        let mut nx = read(lp.verts[s % lp.verts.len()]);
        let d = v3::dot(cur, nx);
        if d.abs() < tol {
            return Err(TopoError::CoreHit { step: s, dot: d.abs() });
        }
        if d < 0.0 {
            sign = -sign;
            nx = v3::scale(nx, -1.0);
        }
        cur = nx;
    }
    Ok(if sign < 0.0 { Z2::NonTrivial } else { Z2::Trivial })
}

/// Half-integer winding of a FRAMED PLANAR loop with normal `nu`:
/// project the LIFTED director onto the plane, accumulate the atan2
/// angle increments, k = sum / 2 pi. Half-integer iff the Z/2 holonomy
/// is nontrivial.
///
/// # Errors
/// [`TopoError::Degenerate`] when the in-plane component drops below
/// `tol` (director parallel to `nu` — the projection pitfall), and
/// [`TopoError::CoreHit`] as in [`z2_holonomy`].
pub fn planar_winding(
    n: &Field3<3>,
    lp: &GridLoop,
    nu: [f64; 3],
    tol: f64,
) -> Result<f64, TopoError> {
    let nuh = v3::normalize(nu);
    let a = if nuh[2].abs() < 0.9 { [0.0, 0.0, 1.0] } else { [1.0, 0.0, 0.0] };
    let u1 = v3::normalize(v3::cross(a, nuh));
    let u2 = v3::cross(nuh, u1);
    let read = |v: [usize; 3]| v3::normalize(n.at(v[0], v[1], v[2]));
    let angle_of = |m: [f64; 3]| -> Result<f64, TopoError> {
        let (x, y) = (v3::dot(m, u1), v3::dot(m, u2));
        if (x * x + y * y).sqrt() < tol {
            return Err(TopoError::Degenerate(
                "director parallel to the loop normal".into(),
            ));
        }
        Ok(y.atan2(x))
    };
    let mut cur = read(lp.verts[0]);
    let mut theta_prev = angle_of(cur)?;
    let mut total = 0.0_f64;
    for s in 1..=lp.verts.len() {
        let mut nx = read(lp.verts[s % lp.verts.len()]);
        let d = v3::dot(cur, nx);
        if d.abs() < CORE_TOL {
            return Err(TopoError::CoreHit { step: s, dot: d.abs() });
        }
        if d < 0.0 {
            nx = v3::scale(nx, -1.0);
        }
        let theta = angle_of(nx)?;
        let mut dth = theta - theta_prev;
        while dth > PI {
            dth -= 2.0 * PI;
        }
        while dth < -PI {
            dth += 2.0 * PI;
        }
        total += dth;
        theta_prev = theta;
        cur = nx;
    }
    Ok(total / (2.0 * PI))
}

/// A traced disclination line on the dual lattice: the ordered centers
/// of the pierced plaquettes it threads. `closed` distinguishes loops
/// from boundary-terminated lines (Z/2 lines cannot end in the bulk).
#[derive(Clone, Debug)]
pub struct DefectLine {
    /// Pierced-plaquette centers along the line.
    pub points: Vec<[f64; 3]>,
    /// Closed loop (true) or boundary-terminated (false).
    pub closed: bool,
}

/// Sign-gauge holonomy around one plaquette given its four node values.
fn plaquette_sign(vals: [[f64; 3]; 4], tol: f64) -> Result<f64, TopoError> {
    let mut s = 1.0_f64;
    for e in 0..4 {
        let d = v3::dot(
            v3::normalize(vals[e]),
            v3::normalize(vals[(e + 1) % 4]),
        );
        if d.abs() < tol {
            return Err(TopoError::CoreHit { step: e, dot: d.abs() });
        }
        if d < 0.0 {
            s = -s;
        }
    }
    Ok(s)
}

/// Trace all disclination lines of an RP2 director field: evaluate the
/// Z/2 holonomy of every elementary plaquette; pierced plaquettes'
/// dual edges are chained through the cubes into lines. Structural
/// self-check (dd = 0 analogue): every cube must contain an EVEN number
/// of pierced faces — Z/2 lines cannot end in a cell — and every chain
/// must close or terminate on the domain boundary.
///
/// # Errors
/// [`TopoError::CoreHit`] if any plaquette transport degenerates,
/// [`TopoError::Degenerate`] on an odd pierced-face count (structural
/// violation) or a line junction (>= 4 pierced faces in one cube,
/// unresolvable without finer data).
pub fn disclination_lines(n: &Field3<3>, tol: f64) -> Result<Vec<DefectLine>, TopoError> {
    let [n0, n1, n2] = n.n;
    // A pierced face: the two cubes it separates (None = outside), and
    // its center.
    struct Face {
        below: Option<[usize; 3]>,
        above: Option<[usize; 3]>,
        center: [f64; 3],
    }
    let mut faces: Vec<Face> = Vec::new();
    // Plaquettes with normal along `ax`, spanned by the other two axes.
    for ax in 0..3usize {
        let bx = (ax + 1) % 3;
        let cx = (ax + 2) % 3;
        let mut dims = [0usize; 3];
        dims[ax] = n.n[ax];
        dims[bx] = n.n[bx] - 1;
        dims[cx] = n.n[cx] - 1;
        for a in 0..dims[ax] {
            for b in 0..dims[bx] {
                for c in 0..dims[cx] {
                    let mk = |db: usize, dc: usize| -> [usize; 3] {
                        let mut v = [0usize; 3];
                        v[ax] = a;
                        v[bx] = b + db;
                        v[cx] = c + dc;
                        v
                    };
                    let corners = [mk(0, 0), mk(1, 0), mk(1, 1), mk(0, 1)];
                    let vals: [[f64; 3]; 4] =
                        core::array::from_fn(|e| n.at(corners[e][0], corners[e][1], corners[e][2]));
                    if plaquette_sign(vals, tol)? > 0.0 {
                        continue;
                    }
                    // Cube (I, J, K) spans nodes I..I+1 etc.; the
                    // plaquette at layer `a` of axis `ax` separates the
                    // cubes with ax-index a-1 and a.
                    let cube = |aa: isize| -> Option<[usize; 3]> {
                        if aa < 0 || aa as usize >= n.n[ax] - 1 {
                            return None;
                        }
                        let mut v = [0usize; 3];
                        v[ax] = aa as usize;
                        v[bx] = b;
                        v[cx] = c;
                        Some(v)
                    };
                    let p0 = n.pos(corners[0][0], corners[0][1], corners[0][2]);
                    let mut center = p0;
                    center[bx] += 0.5 * n.h;
                    center[cx] += 0.5 * n.h;
                    faces.push(Face {
                        below: cube(a as isize - 1),
                        above: cube(a as isize),
                        center,
                    });
                }
            }
        }
    }
    // Cube -> pierced-face incidence; even-count structural check.
    let mut incidence: HashMap<[usize; 3], Vec<usize>> = HashMap::new();
    for (fi, f) in faces.iter().enumerate() {
        for cube in [f.below, f.above].into_iter().flatten() {
            incidence.entry(cube).or_default().push(fi);
        }
    }
    for (cube, fl) in &incidence {
        if fl.len() % 2 != 0 {
            return Err(TopoError::Degenerate(format!(
                "odd pierced-face count {} in cube {cube:?} (a Z/2 line would end there)",
                fl.len()
            )));
        }
        if fl.len() > 2 {
            return Err(TopoError::Degenerate(format!(
                "line junction: {} pierced faces in cube {cube:?}",
                fl.len()
            )));
        }
    }
    // Walk chains: boundary-started lines first, then interior cycles.
    let mut used = vec![false; faces.len()];
    let mut lines: Vec<DefectLine> = Vec::new();
    let walk = |start: usize, mut cube: Option<[usize; 3]>, used: &mut Vec<bool>| -> DefectLine {
        let mut points = vec![faces[start].center];
        used[start] = true;
        let mut cur = start;
        loop {
            let Some(cb) = cube else {
                return DefectLine { points, closed: false };
            };
            let fl = &incidence[&cb];
            let nxt = fl[usize::from(fl[0] == cur)];
            if nxt == start {
                return DefectLine { points, closed: true };
            }
            used[nxt] = true;
            points.push(faces[nxt].center);
            let f = &faces[nxt];
            cube = if f.below == Some(cb) { f.above } else { f.below };
            cur = nxt;
        }
    };
    for fi in 0..faces.len() {
        if used[fi] || (faces[fi].below.is_some() && faces[fi].above.is_some()) {
            continue;
        }
        // Start from the boundary side and walk inward.
        let into = if faces[fi].below.is_none() { faces[fi].above } else { faces[fi].below };
        lines.push(walk(fi, into, &mut used));
    }
    for fi in 0..faces.len() {
        if used[fi] {
            continue;
        }
        lines.push(walk(fi, faces[fi].above, &mut used));
    }
    Ok(lines)
}
