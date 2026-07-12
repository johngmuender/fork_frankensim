//! GUM core Phase B1 — `fs-gum-field`: the stored SU(2)/quaternion field type
//! and its measurement kernels.
//!
//! This crate is the Rust port of the *forward* (measurement) path of the
//! Tier-4A Python engine `analysis/gum-sandbox/tier4-field/field3d_solve.py`
//! plus the Tier-2b 1-D radial hedgehog solver
//! `analysis/gum-sandbox/tier2-closure/radial_solve.py`, grown from the
//! certified Tier-4B skeleton `fs-cosserat-pilot`.  Descent machinery
//! (gradients, arrested Newton flow, penalties) is deliberately OUT of scope
//! here — the data layout keeps the padded-ghost SoA structure of the Python
//! engine so the adjoint stencils port 1:1 in a later phase.
//!
//! Conventions (Tier-2 unit map, `axi_RESULTS.md`, a6 = a0 = m = 1):
//!
//! ```text
//! field      q = (q0,q1,q2,q3), |q| = 1, vacuum q = (1,0,0,0)
//! currents   a_i = conj(q) * D_i q   (Maurer-Cartan; |a_i|^2 = |D_i q|^2
//!            for unit q; the scalar part of a_i is a discretisation
//!            diagnostic, identically 0 in the continuum)
//! E2 = (1/4pi) INT sum_i |D_i q|^2
//! E4 = (1/4pi) INT sum_{i<j} (|D_i q|^2 |D_j q|^2 - (D_i q . D_j q)^2)
//! E6 = pi^3 INT b^2,   b = -(1/2 pi^2) det[q, D_x q, D_y q, D_z q]
//!      (first-row-expansion det4, the fs-cosserat-pilot form; analytically
//!      identical to field3d's 6-term PAIRS expansion; sign fixed so the
//!      f: pi -> 0 hedgehog has degree +1)
//! E0 = (1/4pi) INT (1 - q0)
//! I  = INT 2 (q1^2 + q2^2)
//! degree = INT b
//! E_static = t (E2 + E4) + E6 + E0,   t = 0.008276434949296802 (eps = 0.05)
//! weighted Derrick virial: t E2 - t E4 - 3 E6 + 3 E0 = 0
//! ```
//!
//! Determinism contract: all sweeps and reductions are plain sequential f64
//! loops in fixed (component-outer where stated, i, j, k ascending) order —
//! no tree reductions, no parallelism, no platform libm in kernels
//! (elementary functions come from `fs_math::det`; IEEE `sqrt` and the
//! `cbrt`-for-R* precedent of fs-cosserat-pilot are the only intrinsics).
//! Two runs of any pipeline built on this crate must agree bit-for-bit.
//!
//! Epistemic notice (binding): everything here is within-model numerical
//! engineering on a speculative theory's functional.  A passing gate
//! certifies the discretisation and the port — never anything about nature.

pub mod field;
pub mod radial;
pub mod sector;
pub mod stencil;

pub use field::{Field3, GHOST};
pub use sector::{measure, virial, Scheme, Sectors};

/// eps = 0.05 dial value of t = a2 = a4 (Tier-2 step 1, frozen).
pub const T_FROZEN: f64 = 0.008276434949296802;

/// Compacton radius R* = 2^(5/6) for a6 = a0 = m = 1, evaluated as
/// sqrt(cbrt(32)) (R*^2 = 2^(5/3) = cbrt(32); fs-cosserat-pilot precedent).
#[must_use]
pub fn rstar() -> f64 {
    rstar_sq().sqrt()
}

/// R*^2 = 2^(5/3), evaluated once, deterministically (pilot convention).
#[must_use]
pub fn rstar_sq() -> f64 {
    32.0_f64.cbrt()
}

/// Continuum Bogomolny floor of E6 + E0 per unit degree: 32 sqrt(2) / 15.
#[must_use]
pub fn bps_floor() -> f64 {
    32.0 * core::f64::consts::SQRT_2 / 15.0
}

/// Hamilton product of two quaternions stored as [w, x, y, z].
#[must_use]
pub fn qmul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [
        a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
        a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
        a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1],
        a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0],
    ]
}

/// Quaternion conjugate (= inverse for unit quaternions).
#[must_use]
pub fn qconj(q: [f64; 4]) -> [f64; 4] {
    [q[0], -q[1], -q[2], -q[3]]
}

/// Maurer-Cartan currents a_i = conj(q) * D_i q for the three spatial
/// directions.  Full 4-component quaternions are returned: entries 1..4 are
/// the pure-vector su(2) parts, entry 0 (the scalar part) is a
/// discretisation diagnostic (identically 0 in the continuum for unit q).
/// The identity |a_i|^2 = |D_i q|^2 (unit q) is what lets E2 be summed from
/// raw derivatives.
#[must_use]
pub fn mc_currents(q: [f64; 4], d: &[[f64; 4]; 3]) -> [[f64; 4]; 3] {
    let qc = qconj(q);
    [qmul(qc, d[0]), qmul(qc, d[1]), qmul(qc, d[2])]
}
