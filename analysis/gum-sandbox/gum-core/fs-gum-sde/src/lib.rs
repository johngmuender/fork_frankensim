//! fs-gum-sde — Phase B4 of the GUM physics core: the stochastic
//! (Nelson / stochastic-QM) sector.
//!
//! Spec: `analysis/gum-sandbox/gum-core/gap_survey_v3.json`, dimension
//! "Nelson/stochastic-QM". Port source for the wavefunction machinery:
//! `analysis/gum-sandbox/tier3-born/born.py` (Tier-3 referee:
//! `tier3-born/REPORT.md`).
//!
//! Modules:
//! - [`wiener`]: `WienerStream` — the (trajectory, step-slot, component) →
//!   draw-index contract on fs-rand's counter-based Philox streams. Fixed
//!   consumption per slot makes ensembles order-independent and every single
//!   trajectory replayable in isolation.
//! - [`sde`]: the [`sde::Sde`] problem trait and the [`sde::Stepper`]
//!   integrator trait with Euler–Maruyama and the additive-noise SRA1-class
//!   order-1.5 stochastic Runge–Kutta (with the correlated ΔZ integral).
//! - [`modes`]: born.py's `ModeSystem` — exact mode-sum ψ(t) on the 2-D box
//!   [0,π]² (ħ = m = 1), its log-derivative ∇ψ/ψ, |ψ|² grids and point
//!   evaluations, and the embedded born.py phase constants (numpy seed 42).
//! - [`nelson`]: the Nelson SDE dx = (v + u) dt + √(ħ/m) dW with the current
//!   velocity v = Im(∇ψ/ψ) AND the osmotic velocity u = Re(∇ψ/ψ),
//!   reflecting box walls, and the substepped macro-step with drift
//!   displacement cap (born.py's guard structure adapted to noise).
//! - [`hfunc`]: the coarse-grained H̄ = Σ P̄ ln(P̄/Q̄) functional
//!   (empty-cell skip), the noise-floor formula cg²/(2N), and born.py's
//!   log-space exponential decay fit.
//!
//! Everything is built on fs-math `det` functions and fs-rand Philox
//! streams: bit-identical replay across runs and ISAs by construction.
//! No external crates; std only (campaign closure doctrine).
//!
//! Epistemic scope: this crate validates the ENGINE (SDE integrators +
//! relaxation diagnostics) within the model. The H-theorem import is the
//! corpus's [IM] ingredient — replicating it here says nothing about the
//! substrate ontology (see Tier-3 REPORT.md's epistemic notice).

pub mod hfunc;
pub mod modes;
pub mod nelson;
pub mod sde;
pub mod wiener;
