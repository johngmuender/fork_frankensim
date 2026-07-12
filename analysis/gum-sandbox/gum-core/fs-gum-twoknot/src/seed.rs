//! Two-knot product-ansatz seeding and the App I.2 relative orientations.
//!
//! Product ansatz (paper App C.3-C.4 "product-ansatz bookkeeping"; the
//! coordinator-fixed form): q(x) = q1(x - d/2 X) * q2'(x + d/2 X)
//! (quaternion Hamilton product; degree adds -> B = 2, verified by the
//! degree diagnostic).  Knot 2 is isorotated by the channel-theorem
//! orientation R: q2' = a (x) q2 (x) conj(a) with a the unit quaternion of
//! the isorotation (q0 invariant, vector part rotated) — the standard
//! U -> A U A^dagger Skyrme isorotation.
//!
//! Channels (Thm VII.2 channel theorem, Sec VII.D):
//!  * ATTRACTIVE:  theta = pi, axis n perpendicular to the separation X
//!    (the pi/perp rule) -> a = (0, 0, 0, 1) for X = x-hat, n = z-hat;
//!  * REPULSIVE (strong): theta = pi, n parallel X -> a = (0, 1, 0, 0);
//!  * ALIGNED (weak repulsion, the intermediate channel): theta = 0 ->
//!    a = identity.
//!
//! Separation axis = grid x (the long axis).  Knot centres +- d/2 with
//! d/2 = m h an integer multiple of h, so both centres sit on cell corners
//! for every separation (and for the single-knot reference at the origin):
//! one sub-grid alignment class for the whole campaign, cancelling the
//! per-knot discretisation offset in every energy difference.

use fs_gum_field::radial::{interp, RadialProfile};
use fs_gum_field::{qconj, qmul};
use fs_math::det;

use crate::field_a::FieldA;

/// Relative orientation channel (App I.2's three orientations).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Channel {
    /// theta = pi about z-hat (perpendicular to separation): attractive.
    Attract,
    /// theta = pi about x-hat (parallel to separation): strong repulsion.
    Repulse,
    /// Identity (aligned): weak repulsion.
    Align,
}

impl Channel {
    #[must_use]
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "attract" => Some(Channel::Attract),
            "repulse" => Some(Channel::Repulse),
            "align" => Some(Channel::Align),
            _ => None,
        }
    }

    #[must_use]
    pub fn as_str(&self) -> &'static str {
        match self {
            Channel::Attract => "attract",
            Channel::Repulse => "repulse",
            Channel::Align => "align",
        }
    }

    /// Isorotation unit quaternion a (q2' = a q2 conj(a)).
    #[must_use]
    pub fn quat(&self) -> [f64; 4] {
        match self {
            Channel::Attract => [0.0, 0.0, 0.0, 1.0],
            Channel::Repulse => [0.0, 1.0, 0.0, 0.0],
            Channel::Align => [1.0, 0.0, 0.0, 0.0],
        }
    }
}

/// Hedgehog quaternion at position (x, y, z) relative to the knot centre,
/// from the stored radial profile (Field3::sample_hedgehog pointwise form:
/// same interp clamp, same det trig, same r guard).
#[inline]
#[must_use]
pub fn hedgehog_q(prof: &RadialProfile, x: f64, y: f64, z: f64) -> [f64; 4] {
    let r = (x * x + y * y + z * z).sqrt();
    let f = interp(r, &prof.r, &prof.f);
    let q0 = det::cos(f);
    let sfr = if r > 1.0e-300 { det::sin(f) / r } else { 0.0 };
    [q0, sfr * x, sfr * y, sfr * z]
}

/// Seed the two-knot product ansatz: centres at (-s, 0, 0) and (+s, 0, 0),
/// knot 2 (at +s) isorotated by the channel quaternion.  The product of two
/// unit quaternions is unit; `renormalize()` afterwards is hygiene only.
pub fn seed_two_knot(f: &mut FieldA, prof: &RadialProfile, s: f64, ch: Channel) {
    let a = ch.quat();
    let ac = qconj(a);
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
    for i in 0..nx {
        let x = f.x(i);
        for j in 0..ny {
            let y = f.y(j);
            for k in 0..nz {
                let z = f.z(k);
                // q1 centred at -s (evaluated at x - (-s)), q2 at +s
                let q1 = hedgehog_q(prof, x + s, y, z);
                let q2 = hedgehog_q(prof, x - s, y, z);
                let q2r = qmul(qmul(a, q2), ac);
                f.set(i, j, k, qmul(q1, q2r));
            }
        }
    }
}

/// Seed a single centred knot (the reference run, same grid class).
pub fn seed_single_knot(f: &mut FieldA, prof: &RadialProfile) {
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
    for i in 0..nx {
        let x = f.x(i);
        for j in 0..ny {
            let y = f.y(j);
            for k in 0..nz {
                let z = f.z(k);
                f.set(i, j, k, hedgehog_q(prof, x, y, z));
            }
        }
    }
}

/// Effective separation diagnostic: (1 - q0)-weighted x-centroid of the
/// x > 0 half minus that of the x < 0 half.  Tracks whether the relaxed
/// pair drifted from the seeded separation (merging / sliding caveat).
#[must_use]
pub fn effective_separation(f: &FieldA) -> f64 {
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
    let (mut wp, mut xp, mut wm, mut xm) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    for i in 0..nx {
        let x = f.x(i);
        for j in 0..ny {
            for k in 0..nz {
                let q = f.get(i, j, k);
                let w = 1.0 - q[0];
                if x > 0.0 {
                    wp += w;
                    xp += w * x;
                } else {
                    wm += w;
                    xm += w * x;
                }
            }
        }
    }
    xp / wp - xm / wm
}
