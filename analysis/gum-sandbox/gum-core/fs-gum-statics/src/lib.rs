//! GUM core Phase E1 — `fs-gum-statics`: the SOLVER half of the
//! Skyrme/chiral statics dimension (gap_survey_v3.json), built on the
//! finished measurement crate `fs-gum-field`.
//!
//! This crate is the Rust port of the gradient / adjoint / arrested-Newton-
//! flow machinery of the Tier-4A Python engine
//! `analysis/gum-sandbox/tier4-field/field3d_solve.py` (`Engine.energy`
//! gradient path, `anf`, the near-BPS lattice guards, `halo_fraction`, the
//! clock bisection).  The forward measurement path (stored field, stencils,
//! sector energies, 1-D radial solver) lives in `fs-gum-field` and is reused
//! unchanged; the padded-ghost SoA layout was kept identical to the Python
//! engine there precisely so the adjoint stencils here port 1:1.
//!
//! Frozen algorithms and parameters (survey "Skyrme/chiral" entry):
//!
//! ```text
//! (2+4) flux    P[i,a] = (2/4pi) [ (c2 + c4 (nn_o1 + nn_o2)) D_i q_a
//!                                  - c4 (D_i.D_o1) D_o1 q_a
//!                                  - c4 (D_i.D_o2) D_o2 q_a ]
//!               (c2 = c4 = t reproduces field3d's w24 = 2t/4pi form)
//! sextic        6-term pair expansion of det[q, Dx q, Dy q, Dz q]
//!               (PAIRS table), det-weight w6 = det/(2 pi); the one-sided
//!               guards enter ONLY through this weight and the E0 density
//! E0 density    -1/(4 pi) on q0
//! Routhian      -(L^2 / 2 I^2) dI/dq = -(L^2 / 2 I^2) 4 q_{1,2}
//! adjoints      4th-order stencil transpose with zero ghost flux;
//!               corner average/difference operator transposes
//! projection    g -= (g . q) q   (pointwise tangent projection)
//! guards        E_pen  = (MU/2)   min(0, deg  - (deg_ref - 0.005))^2, MU  = 5000
//!               E_fpen = (MU_F/2) min(0, fgap - fgap_ref)^2,          MU_F = 400
//!               fgap = E6 + E0 - (32 sqrt2 / 15) deg / deg_ref
//!               penalties enter ONLY the descent objective; all reported
//!               energies are the bare functionals with the penalties logged
//! ANF           velocity field, acceleration = -grad, revert + zero-velocity
//!               arrest on objective increase, adaptive dt (grow 1.01x to
//!               dt_max = 0.05, shrink 0.6x on arrest, stall below 1e-7),
//!               per-step renormalisation + velocity tangent re-projection
//! halo          fraction of I outside 1.5 R* about the (1 - q0) centroid;
//!               threshold kappa_c = 1/sqrt(8 pi) = 0.19947 (F-R5)
//! clock         L^2 = (2/3) I E_static  ->  E_rot/E_tot = 1/4
//! ```
//!
//! Determinism contract (inherited from fs-gum-field): all sweeps and
//! reductions are plain sequential f64 loops in fixed ascending
//! (point, component) order — no tree reductions, no parallelism, no
//! platform libm in kernels.  The perturbation generator is a hand-rolled
//! 64-bit LCG with documented constants (`diag::Lcg`); determinism matters
//! more than generator pedigree.  Two runs of any pipeline built on this
//! crate agree bit-for-bit (asserted by the gate binary, gate G-E).
//!
//! Epistemic notice (binding): everything here is within-model numerical
//! engineering on a speculative theory's functional.  The halo-descent gate
//! replicates the campaign's F-R5 mechanism as an ENGINE CAPABILITY — a
//! passing gate certifies the discretisation, the gradients, and the port,
//! never anything about nature.

pub mod anf;
pub mod diag;
pub mod engine;

pub use anf::{anf, AnfParams, AnfResult, Record, Status};
pub use engine::{eval, Anchor, Grad, Opts, Out, Wall, SGN6};

/// One-sided degree-anchor strength MU (field3d `MU_DEG`, frozen).
pub const MU_DEG: f64 = 5000.0;

/// Degree-anchor dead band (field3d `DEG_BAND`, frozen).
pub const DEG_BAND: f64 = 0.005;

/// One-sided Bogomolny-floor wall strength MU_F (field3d `MU_FLOOR`, frozen).
pub const MU_FLOOR: f64 = 400.0;

/// Floor-wall dead band below the engine hedgehog's own floor gap
/// (field3d `FLOOR_BAND`, frozen): fgap_ref = hedgehog gap - FLOOR_BAND.
pub const FLOOR_BAND: f64 = 0.01;

/// Halo threshold kappa_c = 1/sqrt(8 pi) = 0.19947... (Step 3, F-R5).
#[must_use]
pub fn kappa_threshold() -> f64 {
    1.0 / (8.0 * core::f64::consts::PI).sqrt()
}

/// Halo diagnostic radius 1.5 R* (field3d `HALO_RADIUS`).
#[must_use]
pub fn halo_radius() -> f64 {
    1.5 * fs_gum_field::rstar()
}
