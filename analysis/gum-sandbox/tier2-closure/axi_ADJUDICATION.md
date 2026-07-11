# Tier 2b Step 2 — Spheroidal Isorotating Closure: Coordinator's Adjudication

**The headline: FINDING F-R4 — the first substantive within-model stress on
the corpus's centerpiece.** Everything below is grounded in `axi_results.json`
(validation gates, two deformation families, ε→0 endpoints, ε-scan) computed
with **zero fitted constants** from real field-level quadratures of the
Maurer–Cartan currents.

## Validation gates (all load-bearing ones pass)
- Inertia 2-D vs 1-D radial formula: agreement to **1.6×10⁻⁷** (V2).
- Field-level g(0.42) vs the paper's spheroid formula: **1.285547004 vs
  1.285547005 — 9 digits** (V4). The unit map validated independently: the
  analytic scaling-rung endpoint reproduces `2√(ê₀𝔦₀) = 2.0532877` to 7
  digits (endA_analytic).
- SDiff-flatness of E₆+E₀ under pullback: exact to ~10⁻¹⁵ per sector at the
  analytic level (V3_cross final: E₀ drift 1.6×10⁻¹⁵, I drift 1.2×10⁻¹⁵).
- V1 (2-D↔1-D sector integrals): worst sector −0.62% vs the 0.5% target —
  marginal, discretization-level, accepted with note.

## FINDING F-R4 (substantive): the SDiff inertia-enhancement mechanism
fails in its literal reading — measured at machine precision

The paper's shape-exact endpoint (`g = 3/2`, `𝔠₀ = 128√42/105π = 2.515`)
rests on App. G.5's argument: the isorotation inertia, taken "over volume-V
SDiff images," can approach `sup g = 3/2` in the extreme-oblate limit while
the (6+0) energy stays flat (B.3). The direct numerical fact from this
solve: **under volume-preserving (SDiff) pullback, the isorotation inertia
is invariant identically** — measured drift 1.2×10⁻¹⁵ at λ = 0.6 — for the
same reason ∫sin²f is invariant: the inertia density `2(q₁²+q₂²)` is a
function of the field *value*, and pullbacks preserve value distributions.
**Both the BPS energy and the inertia are frozen on the SDiff orbit; the
g-dial cannot be turned there.** The deformation class that *does* turn it
(contour deformation with the hedgehog direction held — family B, which
reproduces the paper's own g(λ) formula to 9 digits) is **not** SDiff and
is **not** energy-flat.

Consequences, measured:
- **Family A (honest SDiff):** g ≡ 1 exactly; the ε→0 closure lands on the
  paper's *scaling rung* — `𝔠 = 2.0533`, `κ = 0.9354`, `V = √2` exact — not
  the claimed shape-exact endpoint (2.515, 0.764).
- **Family B (g-enhancing, non-SDiff):** ε→0 closes at `𝔠 = 2.169`,
  `g = 1.078`, `λ = 0.817`; at ε = 0.05: `𝔠 = 2.279`, `κ = 0.907`,
  `λ* = 0.818`, `g* = 1.078`, `V = 1.461`.
- **Against the ⟨r1⟩ benchmark (2.37 ± 0.09, 0.802 ± 0.018, λ* = 0.42,
  g* = 1.31, V = 1.409):** 𝔠 within 1.0σ (family B), but **κ +5.8σ,
  V +5σ, λ*/g* far off** — the benchmark tuple is NOT reproduced in either
  family with zero fitted constants.
- Invariants that DO hold everywhere: `E_rot/E = ¼` to 10⁻⁹ (automatic, as
  the paper claims); `V(ε→0) = √2` exact in family A; clock residuals ≤ 10⁻⁸.

## What F-R4 does and does not establish (the honest boundary)
- It **refutes the stated mechanism** (G.5's SDiff supremum) in its literal
  compositional reading, at machine precision, by the corpus's own
  deformation vocabulary.
- It does **not** refute the benchmark *number*: my families are restricted;
  the corpus's ⟨r1⟩ came from an unrestricted 2-D field solve, which may
  find inertia-gaining modes cheaper than family B's (they just cannot be
  SDiff modes). Under the corpus's own rules this is a **named doubt with a
  discharge path**, not a kill: the unrestricted solve must now show which
  non-SDiff modes supply g* = 1.31 at ε = 0.05 and g → 3/2 as ε → 0, and at
  what energy cost — or the deep-BPS endpoint 𝔠₀ = 2.515 reverts toward the
  scaling rung and the ħ-closure's physical constants shift by ~20%.
- Alternative reading kept open: G.5's "SDiff images" may be intended
  non-compositionally (the corpus's text is one paragraph); but its own
  proof sketch ("∫sin²f is SDiff-invariant and sin²θ ≤ 1") uses exactly the
  invariance that freezes the inertia in the compositional reading.

## Status
Tier 2b Step 2 **complete** as a restricted-family closure: unit map
validated to 7–9 digits, ¼ and √2 anchors exact, the scaling rung
reproduced exactly — and the shape-exact endpoint and the ⟨r1⟩ fine
structure now carry **F-R4**, the campaign's first finding that bears on
the centerpiece rather than on print-level arithmetic. The unrestricted
2-D solve (Tier 2b Step 3) is now not just the benchmark's replication but
**F-R4's adjudicator**.
