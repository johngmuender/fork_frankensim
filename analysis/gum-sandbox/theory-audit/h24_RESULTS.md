# H2.4 — A2b-1 helical counterexample (the null-transport "exactness", broken constructively)

**Campaign:** replication Phase H2, follow-on №2 of the theory audit (`T1_NR_gates.md` §1.2/§6).
**Date:** 2026-07-16. **Code:** `h24_helical.py` (python3 + sympy/numpy/scipy, deterministic, runtime < 1 s).
**Raw output:** `h24_results.json`. **Units:** c = ω = 1, so ƛ = c/ω = 1.
**Status: measured/derived facts below; §6 marked as analysis; promotion is the coordinator's call.**

## 0. Verdict in one line

**The counterexample exists, is fully explicit, and SURVIVES every textually-plausible
premise-completion we could construct — zero net axial momentum, zero net axial current,
compact support, closed streamlines, exact z-parity of the current field, every carrier
orbiting at the entrained ω, and Route B's current counting I = qω/2π — while giving
⟨r_⊥²⟩_j = 0.33 ƛ² (annulus), 0.57 ƛ² (knotted torus), or a² for any a ≤ ƛ (thin shell).
The only completions that restore ⟨r_⊥²⟩_j = ƛ² are restatements of the contested clause
itself (pointwise j_z ≡ 0 / v = the pattern's rigid-rotation field). Theorem A2b-1's "exact"
is therefore underdetermined by its printed premises: the audit's GAP is upgraded to a
demonstrated underdetermination, with the needed repair pinned.**

## 1. The premises, read precisely from NR-A2b §1 and used as stated

- **P1 (null transport):** "exact masslessness … leaves no sub-luminal branch for a stationary
  transported structure built of channel content; **|v| = c on every occupied streamline**."
- **P2 (single-frequency stationarity):** "Stationarity at one frequency (the entrained clock,
  ω = E/ħ exact, IV.J) means the transported pattern co-rotates" — which the instrument
  continues "**: v = ωr_⊥ azimuthally**". The premise is the co-rotation of the *pattern*;
  the continuation "v = ωr_⊥ azimuthally" is the step under audit (T1 §1.2), and we do NOT
  assume it.
- **Implicit premises honored:** charge conservation (continuity), smooth single-valued
  velocity field on the support, nonnegative density.
- **Route B (§1):** "the winding-counter/flux argument … I = qω/2π; Φ_eff = full co-rotation
  flux" — the I-counting is honored by the counterexamples (measured below); the Φ_eff = πƛ²B
  step presumes the cylinder radius (NR-A1 §3.3: "the routes share the cylinder lemma").

## 2. CE-1 — the balanced helical annulus (measured facts)

Two coaxial annular families, ρ = f(r)(1 + cos(φ−ωt)), **v = ωr φ̂ ± √(c²−ω²r²) ẑ**:
family A r ∈ [0.55, 0.75] drifting +z; family B r ∈ [0.35, 0.55] drifting −z, weight 1.2242
chosen to balance axial momentum.

| check | result |
|---|---|
| \|v\| = c pointwise | symbolic 0; numeric ≤ 1.1×10⁻¹⁶ |
| continuity ∂_tρ + ∇·(ρv) = 0 | symbolic residual 0 (any radial profile f, either drift sign) |
| pattern co-rotation | j(r,φ,z,t) = J(r, φ−ωt, z) rigidly; v itself rotation-equivariant |
| lab-frame spectrum at a probe point | power at {0, ω} only; all other bins < 10⁻²⁴ (relative) |
| streamlines | helices: r_⊥ drift 1.8×10⁻¹¹ over 3 turns (RK45); fitted v_z exact to 10⁻⁹; **dφ/dt = ω on every streamline** |
| net axial momentum / current | 0.0 (balanced; quadrature exact) |
| Route B current: I/(qω/2π) | 1.000000000000 (time-averaged) |
| **⟨r_⊥²⟩_ρ = ⟨r_⊥²⟩_{\|j\|}** | **0.331580 ƛ² — a 67% deficit against the theorem's "= ƛ² exactly"** |

(ρ- and |j|-weighting coincide because |j| = ρc pointwise under P1.)

## 3. CE-2 / CE-3 — the compact, closed-streamline, parity-symmetric version (measured facts)

Support = torus r_⊥ = R₀ + b cosΘ, z = b sinΘ (R₀ = 0.7); **v = ωr_⊥ φ̂ + √(c²−ω²r_⊥²) Θ̂**;
σ ∝ (1+cos(φ−ωt))/(r_⊥ v_p). With b tuned so ωT_poloidal = 2π/3 (b* = 0.2121037079, brentq to
10⁻¹⁴), **every streamline is a closed (1,3) torus knot traversed in exactly one clock period
2π/ω** (closure |X(2π/ω) − X₀| = 7.6×10⁻¹²; poloidal turns = 3.000000000; speed = 1 to 10⁻¹²
along the orbit; on-surface deviation 5×10⁻¹³).

- Surface continuity: symbolic residual 0. Spectrum {0, ω} (same rigid co-rotation).
- Net axial momentum: 0 exactly (integrand ∝ cosΘ). I = qω/2π exactly (both integrals carry
  the same ∫dΘ/v_p).
- **⟨r_⊥²⟩_j = 0.567690 ƛ² — a 43% deficit.**
- **CE-3 (parity):** CE-2 alone is not z-mirror symmetric (mirror reverses the poloidal
  circulation — parity *does* constrain a single drifting tube). Two CE-2 tori at z = ±h with
  opposite poloidal circulation restore exact parity: v(Mx) = Mv(x) verified to 2.1×10⁻¹⁶ at
  91 sample points; also invariant under full inversion; reversed-torus closure 1.2×10⁻¹¹;
  ⟨r_⊥²⟩ unchanged. So the parity completion does not rescue the theorem either.
- Tunability: a single thin shell at radius a gives ⟨r_⊥²⟩ = a² = (c²−v_z²)/ω² exactly for any
  a ∈ (0, ƛ] — the audit's √(c²−v_z²)/ω-class; a = ƛ recovers the corpus's null ring as the
  boundary case of the family.

## 4. The premise-completion scoreboard (the fairest version of the audit)

| completion (unstated, textually plausible) | kills the counterexample? |
|---|---|
| zero net axial momentum | NO (CE-1 balanced; CE-2/3 exact) |
| zero net axial current | NO (same) |
| compact support | NO (CE-2/3) |
| closed streamlines / periodic orbits | NO (CE-2/3: (1,3) knots, period = one clock period) |
| every carrier orbits at the entrained ω | NO (dφ/dt = ω on every streamline, all CEs) |
| z-parity of the current field | NO (CE-3 mirror pair; it does kill a *single* drifting tube) |
| Route B's I = qω/2π (winding counter) | NO (satisfied exactly, all CEs) |
| pointwise j_z ≡ 0 ("current support is planar") | YES — but this *is* the contested clause |
| v = pattern's rigid-rotation field ωẑ×x | YES — equivalent to the above |

Structural reason the gap is real (derived, not assumed): co-rotating stationarity constrains
the velocity field only up to flows **tangent to the pattern's level sets** — axial drift on a
z-uniform pattern, poloidal recirculation on a poloidally-uniform one are pattern-invisible.
Continuity then fixes the density weighting (σ ∝ 1/(r_⊥v_p)) but never the drift itself.

## 5. What breaks downstream (measured consequence)

With ⟨r_⊥²⟩_j free in (0, ƛ²], NR-A2b's w1 "exact" fails as printed, and everything calibrated
on the null ring inherits the freedom: the moment reading (𝔐 ∝ I·π⟨r_⊥²⟩ shrinks by the same
factor — 0.33× to 0.57× in our examples), the w3 protection quartet's "all current support at
r_⊥ = ƛ identically" condition (T1 §1.2 already graded w3 conditional on A2b-1), and the ¾
calibration's "projection extent = ƛ" cross-reference. Route B does not independently rescue
w1: its current counting is satisfied by every counterexample, so the ƛ² conclusion rests
entirely on the flux step that assumes the cylinder.

## 6. Status and the repair (analysis)

**Counterexample status: SURVIVES all plausible premise-completions.** No completion we could
construct from the corpus's own text kills it short of asserting the azimuthal-only clause
itself. The corpus's repair options, pinned: (i) a new lemma that the channel rule/entrainment
identity forbids drift along pattern-symmetry directions on occupied streamlines (the audit's
"what would close it" — nothing of the sort is printed); or (ii) a variational selection
argument that the stationary solution's current support is planar — which would need an energy
functional for occupied streamlines that NR-A2b does not define. Until one lands, A2b-1 is
correctly graded a demonstrated underdetermination, not a theorem, and NR-A1 §2.4's own honest
earlier grade ("the exactness of (iii) beyond the ideal limit is not yet a theorem") was the
right one. Caveat, printed: "textually plausible" is bounded by our reading of NR-A2b §1 and
NR-A1's cited clauses; a completion living elsewhere in the corpus (e.g. the unreproduced VII.F
operator construction) is not excluded — that is exactly the corpus's self-declared replication
gate, and the gate is the right one.

## 7. Deliverables

- `h24_helical.py` — constructions + all verifications (symbolic and quadrature; < 1 s).
- `h24_results.json` — premise checks, counterexample data, scoreboard, tunability scan.
- This memo.

— end of memo —
