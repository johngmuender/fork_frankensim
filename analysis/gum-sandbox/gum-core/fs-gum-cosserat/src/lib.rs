//! GUM physics core, Phase B3: the linearized micropolar/Cosserat continuum.
//!
//! Line-by-line port of the validated Tier-1 pilot
//! `analysis/gum-sandbox/tier1-spectrum/spectrum.py`, generalized from
//! k ∥ z to arbitrary k ∈ R³, with two genuinely new pieces the survey
//! (gap_survey_v3.json, dimension "Micropolar/Cosserat continuum") called
//! for:
//!
//!   1. a complex-Hermitian eigensolve WITH eigenvectors, built as the
//!      real-symmetric embedding  H = A + iB  →  S = [[A, −B], [B, A]]
//!      on fs-la's deterministic `jacobi_eigh` (eigenvalues come in exact
//!      duplicate pairs; any real eigenvector (x; y) of S recovers the
//!      complex eigenvector x + iy of H), cross-checked against fs-la's
//!      complex QR `eig` as a values-only oracle — both kernels consumed
//!      directly from fs-la (`default-features = false` disables its
//!      fs-exec/asupersync `exec` coupling; Phase F3). Phases B3–E2 used
//!      verbatim vendored copies (src/fsla_vendored.rs, now deleted); the
//!      swap back was gate-verified bit-identical (unchanged Merkle root);
//!   2. a mass-matrix-aware symplectic Verlet (kick–drift–kick with
//!      q̇ = M⁻¹p) evolving plane-wave amplitudes of each branch in the
//!      TIME domain — the repo's first time-domain validation of the
//!      coupled second-order system  ρ₀ü = −K_uu u − K_uφ φ,
//!      J φ̈ = −K_φu u − K_φφ φ.
//!
//! Physics (spectrum.py conventions, exactly):
//!   fields u (displacement) and φ (micro-rotation); plane-wave amplitude
//!   q̂ = (u_x, u_y, u_z, φ_x, φ_y, φ_z) ∈ C⁶;
//!   W₂ = (λ/2)|e_kk|² + μ e_(ij)e_(ij) + μ_c e_[ij]e_[ij]
//!      + (α/2)|Γ_kk|² + (β/2)Γ_(ij)² + (γ/2)Γ_[ij]²
//!      + mass term (Case A: (m_V²/2)|ψ|², ψ = φ − ½ curl u;
//!                   Case B: (m_V²/2)|φ|²)
//!      + χ₃ Re(e_[ij]* Γ_[ij]);
//!   e_ij = ∂_i u_j − ε_ijl φ_l,  Γ_ij = ∂_i φ_j  (∂_i → i k_i);
//!   T = (ρ₀/2)|u̇|² + (J/2)|φ̇|²,  M = diag(ρ₀×3, J×3);
//!   each energy term c·Σ|T_ij|² contributes K += 2c R†R so that
//!   W = ½ q†Kq, and the spectrum solves ω² M v = K v via
//!   M^(−1/2) K M^(−1/2).
//!
//! Phase E2 extension (dimension "Skyrme/chiral", the W_χ part): the full
//! chiral-transduction term of corpus II.C eq (2.2),
//! W_χ = χ₁e_kkΓ_ll + χ₂e_(ij)Γ_(ij) + χ₃e_[ij]Γ_[ij], as
//!
//!   1. symbol assembly for all three couplings (this file's K(k) path,
//!      previously χ₃-only), with the text-unpinned trace-vs-deviatoric
//!      χ₂ convention as a documented switch;
//!   2. a REAL-SPACE energy/gradient module [`wchi`] for stored (u, φ)
//!      fields on a periodic grid, cross-validated against the symbol
//!      path on plane waves (exact discrete identity, gate G22);
//!   3. the Dzyaloshinskii soft-sector condensate functional
//!      [`softsector`] (corpus II.H, the E-H1 audit): margin
//!      𝔪 = χ²/(γΔ²), threshold at 𝔪 = 1, pitch q* = χ/γ, and the
//!      Thm H-3 tilt sinθ_c = √(1 − 1/𝔪).
//!
//! All χ couplings are FREE runtime parameters defaulting to 0 — the
//! χ = 0 default reproduces the Phase B3 battery bit-for-bit and leaves
//! the campaign's χ-free knot-statics functional unchanged. See the
//! module docs of [`wchi`] and [`softsector`] for the precise
//! pinned-vs-free convention ledger.
//!
//! Epistemic notice (binding): every PASS in the gates binary certifies a
//! WITHIN-MODEL property of a speculative theory's linearized dispersion.
//! It validates the port, the eigensolve path, and the integrator — never
//! the physics.

pub mod branches;
pub mod eigh;
pub mod moduli;
pub mod softsector;
pub mod symbol;
pub mod verlet;
pub mod wchi;

pub use branches::{branches, classify_full, fit_loglog, fit_w2, Branches};
pub use eigh::{eigh_hermitian, oracle_check, solve, Modes};
pub use moduli::{MassCase, Moduli};
pub use softsector::SoftSector;
pub use symbol::{build_ops, deviatoric, mass_matrix, stiffness, tensor_parts, Ops};
pub use verlet::{evolve_mode, verlet_mass_step, ModeRun};
pub use wchi::{energy_displaced, grad_dot, wchi_energy, wchi_grad, MicroField};
