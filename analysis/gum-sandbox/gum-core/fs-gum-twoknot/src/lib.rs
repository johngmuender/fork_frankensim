//! GUM core Phase G4 — `fs-gum-twoknot`: the App I.2 two-knot
//! bond-equation flagship (the corpus protocol never executed anywhere).
//!
//! Built ON the finished crates (`fs-gum-field` measurement path,
//! `fs-gum-statics` guard/ANF protocol), neither of which is modified.
//! `Field3` is cubic-only, so this crate carries a minimal anisotropic
//! (nx, ny, nz) twin of the field + engine + ANF, gated BITWISE against
//! fs-gum-statics on a cubic grid by `twoknot_gates` before any campaign
//! run is trusted.
//!
//! Protocol (paper App I.2 reconstructed; grid per scale_survey_v4
//! DIMENSION 4, long axis extended for the far-point anchor):
//!
//! ```text
//! grid        160 x 96 x 96 at the frozen h = 0.09375 (box 15 x 9 x 9;
//!             transverse box = the certified LBOX = 4.5 referee geometry);
//!             deviation from the survey's 144x96x96: +16 long-axis cells so
//!             the x ~ 4 asymptote anchor fits with the full ln(1e6)/mu =
//!             1.78 Yukawa tail margin — one uniform grid class for ALL runs
//! seeding     product ansatz q = q1(x - d/2) q2'(x + d/2), knot 2
//!             isorotated per channel (attract: pi about z-hat perp X;
//!             repulse: pi about x-hat par X; align: identity); centres at
//!             +-(m h) on cell corners for every d — one sub-grid alignment
//!             class, so per-knot discretisation offsets cancel in every
//!             energy difference
//! relax       guarded ANF, corner objective, eps = 0.05
//!             (t = 0.0082764349...), deg_ref = the seed's own degree (~2),
//!             deg_unit = deg_ref / 2 (per-unit-degree Bogomolny wall),
//!             fgap_ref = seed floor gap - 0.01, frozen MU = 5000 / 400,
//!             fixed iteration cap (uniform across runs; tail flatness
//!             recorded)
//! observable  E_int(d) = E(two-knot, d) - 2 E(single knot, same grid),
//!             cross-anchored on E(far point x ~ 4); fitted downstream to
//!             the pair-law screened-dipole form (7.5)
//! ```
//!
//! Determinism: each run is one OS process, plain sequential loops, no RNG
//! anywhere in the seeding or descent (the product ansatz is analytic) —
//! per-run bit-identity is inherited from the fs-gum-statics contract.  The
//! 4-process task farm changes scheduling only, never any number.
//!
//! Epistemic notice (binding): everything here is within-model numerical
//! engineering on a speculative theory's functional.  Measured facts only;
//! adjudication is the coordinator's.

pub mod anf_a;
pub mod engine_a;
pub mod field_a;
pub mod seed;

pub use anf_a::{anf_a, AnfParamsA, AnfResultA, RecordA, StatusA};
pub use engine_a::{eval_a, AnchorA, GradA, OptsA, OutA, WallA, SGN6};
pub use field_a::{aidx_a, corner_a, derivs4_a, FieldA, GHOST};
pub use seed::{effective_separation, hedgehog_q, seed_single_knot, seed_two_knot, Channel};

/// Frozen campaign grid spacing (the certified single-knot referee h:
/// 2 * 4.5 / 96).
pub const H_FROZEN: f64 = 0.09375;

/// Campaign grid: long axis (separation axis).
pub const NX: usize = 160;
/// Campaign grid: transverse axes (the certified LBOX = 4.5 geometry).
pub const NT: usize = 96;

/// Radial-profile protocol (the fs-gum-statics gate values).
pub const RADIAL_N: usize = 4000;
pub const RADIAL_RMAX: f64 = 6.0;

/// Yukawa screening mass mu = 1/sqrt(2t) at the frozen dial (radial referee
/// R4: mu_fit 7.7726 vs this closed form 7.7725).
#[must_use]
pub fn mu_frozen() -> f64 {
    1.0 / (2.0 * fs_gum_field::T_FROZEN).sqrt()
}
