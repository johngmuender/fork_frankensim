//! Analytic referee fixtures (survey validation_referees): the exact
//! eps = 0 BPS compacton hedgehog and its Hopf projection, the +-1/2
//! disclination line fields, the controls, and the linked-circle
//! calibration curves.

use crate::field::{hopf_project_one, Bc, DirectorField, Field3, QField};
use crate::v3;
use core::f64::consts::PI;

/// R*^2 = 2^(5/3) (R* = 2^(5/6)), evaluated once, deterministically —
/// the pilot's convention.
#[must_use]
pub fn rstar_sq() -> f64 {
    32.0_f64.cbrt()
}

/// The exact eps = 0 compacton hedgehog (fs-cosserat-pilot closed form):
/// f0 = 2 arccos(r/R*) gives q0 = 2 r^2/R*^2 - 1 and
/// q_i = (2/R*^2) sqrt(R*^2 - r^2) x_i inside R*; vacuum (1,0,0,0)
/// outside. |q| = 1 identically; degree +1.
#[must_use]
pub fn q_compacton(x: f64, y: f64, z: f64, rs2: f64) -> [f64; 4] {
    let r2 = x * x + y * y + z * z;
    if r2 >= rs2 {
        return [1.0, 0.0, 0.0, 0.0];
    }
    let k = 2.0 * (rs2 - r2).sqrt() / rs2;
    [2.0 * r2 / rs2 - 1.0, k * x, k * y, k * z]
}

/// The compacton stored on the campaign's cell-centred nn^3 grid over
/// [-half, half]^3 with the vacuum ghost (1,0,0,0).
#[must_use]
pub fn compacton_qfield(nn: usize, half: f64) -> QField {
    let rs2 = rstar_sq();
    Field3::cell_centered(nn, half, Bc::Vacuum([1.0, 0.0, 0.0, 0.0]), |x, y, z| {
        q_compacton(x, y, z, rs2)
    })
}

/// The analytic Hopf-1 texture: n = hopf_project(q_compacton), stored
/// directly (H = deg q = 1 exactly; n = z-hat outside R*, so the
/// periodic net-flux obstruction vanishes).
#[must_use]
pub fn hopfion_director(nn: usize, half: f64, bc: Bc<3>) -> DirectorField {
    let rs2 = rstar_sq();
    Field3::cell_centered(nn, half, bc, |x, y, z| {
        hopf_project_one(q_compacton(x, y, z, rs2))
    })
}

/// Normalized vector part of a quaternion field as a director (the S2
/// hedgehog map v-hat = (q1, q2, q3)/|v|; `fallback` where |v| ~ 0,
/// i.e. at the vacuum and at the compacton center — callers must only
/// sample where the map is defined).
#[must_use]
pub fn vector_director(q: &QField, fallback: [f64; 3]) -> DirectorField {
    let conv = |v: [f64; 4]| -> [f64; 3] {
        let w = [v[1], v[2], v[3]];
        let n = v3::norm(w);
        if n < 1.0e-14 {
            fallback
        } else {
            v3::scale(w, 1.0 / n)
        }
    };
    let bc = match q.bc {
        Bc::Periodic => Bc::Periodic,
        Bc::Vacuum(g) => Bc::Vacuum(conv(g)),
    };
    let data = q.data.iter().map(|&v| conv(v)).collect();
    Field3::from_parts(q.n, q.h, q.origin, bc, data)
}

/// Globally target-rotate a director field: n'(x) = R n(x), R the
/// Rodrigues rotation about `axis` by `angle`. H is invariant.
#[must_use]
pub fn rotate_target(f: &DirectorField, axis: [f64; 3], angle: f64) -> DirectorField {
    let k = v3::normalize(axis);
    let (s, c) = angle.sin_cos();
    let rot = |m: [f64; 3]| -> [f64; 3] {
        v3::add(
            v3::add(v3::scale(m, c), v3::scale(v3::cross(k, m), s)),
            v3::scale(k, v3::dot(k, m) * (1.0 - c)),
        )
    };
    let bc = match f.bc {
        Bc::Periodic => Bc::Periodic,
        Bc::Vacuum(g) => Bc::Vacuum(rot(g)),
    };
    let data = f.data.iter().map(|&m| rot(m)).collect();
    Field3::from_parts(f.n, f.h, f.origin, bc, data)
}

/// Domain parity mirror z -> -z (exact on the cell-centred grid):
/// n'(i, j, k) = n(i, j, N-1-k). Orientation-reversing, so H -> -H.
#[must_use]
pub fn mirror_z(f: &DirectorField) -> DirectorField {
    let [n0, n1, n2] = f.n;
    let mut data = Vec::with_capacity(f.data.len());
    for i in 0..n0 {
        for j in 0..n1 {
            for k in 0..n2 {
                data.push(f.at(i, j, n2 - 1 - k));
            }
        }
    }
    Field3::from_parts(f.n, f.h, f.origin, f.bc, data)
}

/// Uniform director (all z-hat): H = 0, all loops trivial.
#[must_use]
pub fn uniform_director(nn: usize, half: f64, bc: Bc<3>) -> DirectorField {
    Field3::cell_centered(nn, half, bc, |_, _, _| [0.0, 0.0, 1.0])
}

/// The +-1/2 disclination line field n = (cos(phi/2), s sin(phi/2), 0),
/// phi = atan2(y, x), s = +-1 — an RP2 line field (discontinuous as a
/// vector across phi = pi, continuous as a line field), with a straight
/// defect line along the z-axis.
#[must_use]
pub fn disclination_field(nn: usize, half: f64, s: f64) -> DirectorField {
    Field3::cell_centered(nn, half, Bc::Vacuum([1.0, 0.0, 0.0]), |x, y, _| {
        let phi = y.atan2(x);
        [(0.5 * phi).cos(), s * (0.5 * phi).sin(), 0.0]
    })
}

/// Periodic slab of 2-D baby skyrmions (a degree-1 (x, y)-texture at
/// every z): carries net B_z flux through the torus 2-cycle, so the
/// Whitehead route MUST reject it with the net-flux obstruction.
#[must_use]
pub fn baby_skyrmion_slab(nn: usize, half: f64) -> DirectorField {
    let rho0 = 0.6 * half;
    Field3::cell_centered(nn, half, Bc::Periodic, |x, y, _| {
        let rho = (x * x + y * y).sqrt();
        if rho >= rho0 {
            return [0.0, 0.0, 1.0];
        }
        let g = PI * (1.0 - rho / rho0);
        let (sg, cg) = g.sin_cos();
        let phi = y.atan2(x);
        [sg * phi.cos(), sg * phi.sin(), cg]
    })
}

/// The standard Hopf-link calibration pair: a unit circle in the
/// xy-plane about the origin and a unit circle in the xz-plane about
/// (1, 0, 0) — |Lk| = 1 exactly. `m` segments per circle.
#[must_use]
pub fn linked_circles(m: usize) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let c1 = (0..m)
        .map(|s| {
            let t = 2.0 * PI * s as f64 / m as f64;
            [t.cos(), t.sin(), 0.0]
        })
        .collect();
    let c2 = (0..m)
        .map(|s| {
            let t = 2.0 * PI * s as f64 / m as f64;
            [1.0 + t.cos(), 0.0, t.sin()]
        })
        .collect();
    (c1, c2)
}

/// An unlinked control pair: the same circles pulled far apart (Lk = 0).
#[must_use]
pub fn unlinked_circles(m: usize) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let (c1, mut c2) = linked_circles(m);
    for p in &mut c2 {
        p[0] += 4.0;
    }
    (c1, c2)
}
