# Tier 2b Step 1 — Radial Field Solve: Coordinator's Adjudication

The pilot (`radial_RESULTS.md`, `radial_solve.py`) executed the campaign's
**first field-level PDE solve**: the 1-D hedgehog profile at ε = 0.05
(gradient flow + damped Newton to max-gradient ~4×10⁻¹³, graded grid,
N vs 2N refinement deltas as error bars).

## Gate verdicts (against the corpus's own App. I.1 acceptance criteria)

| Gate | Corpus criterion | Measured | Verdict |
|------|-----------------|----------|---------|
| Degree | ∫b = 1.000 | K = 1.0000001 | **PASS** |
| Derrick virial | residual ≤ 3×10⁻⁴ | **5.96×10⁻⁷** | **PASS** (500× inside the gate) |
| ε dial | sector ratio = 0.05 | 0.049985 | **PASS** |
| Tail form + mass | Yukawa dipole, μ² = m̃²/(2a_ψ) | fitted μ = 7.772564 vs derived 7.772547 (2×10⁻⁶) | **PASS** |
| Tail amplitude 𝔭 | 0.84 ± 0.03 | nearest candidate A/2π = 0.805 (~1.2σ low) | **CONVENTION-LIMITED** |
| Near-BPS structure | small positive deficit | +0.48% above the Bogomolny bound | **PASS** |

## Adjudication rulings

1. **The paper's tail-mass formula is independently confirmed — blind.**
   My task sheet mis-stated the linearization as μ = m√(a₀/a₂), dropping the
   ½ from 𝒱 = m²(1−cos f) ≈ (m²/2)f². The pilot caught the error and
   re-derived μ = m/√(2a₂) from scratch — which is **exactly the corpus's
   own App. B.5 statement μ² = m̃²/(2a_ψ)**. Since the pilot did not have
   the paper's formula, this is a blind external re-derivation landing on
   the corpus's value, then verified numerically to 2×10⁻⁶. The "√2 defect"
   the pilot attributed to conventions was *my spec's* defect; **the corpus
   is vindicated on this point.** Filed per the print-your-own-defects rule.
2. **The virial identity holds in coefficient-absorbed convention.** The
   paper's App. C.1 "E₂ − E₄ − 3E₆ + 3E₀ = 0" reads with sector energies
   *including* their moduli (the paper's normal usage); the pilot's weighted
   re-derivation (tE₂ − tE₄ − 3E₆ + 3E₀ = 0 with bare integrals) is the
   same identity, satisfied at 6×10⁻⁷. Not a corpus defect — a notation
   ambiguity in my spec, now pinned.
3. **The 𝔭 = 0.84 ± 0.03 amplitude remains unreproducible in normalization**
   — five natural candidates computed (0.365–5.06), the nearest at 0.805
   (~1.2σ low). This is the expected consequence of the frozen ⟨r1⟩
   pipeline's absence from the archive: the tail *form* and *mass* replicate
   exactly; the *amplitude convention* cannot be recovered from the text.
   Recommendation for the corpus: pin 𝔭's definition (this is precisely the
   class of spec-underdetermination the v2 assessment predicted).
4. **The BPS scaffolding behaves as claimed at the field level**: the ε=0.05
   solution sits +0.48% above the (re-derived, normalization-matched)
   Bogomolny bound, with the compacton R* recovered in-convention and the
   edge smoothed into the Yukawa tail — the boundary-layer hazard resolved
   by the graded grid exactly as anticipated.

## What this milestone establishes
The 1-D rung of the ⟨r1⟩ pipeline is now externally replicated at the field
level: profile, topology, virial, near-BPS deficit, and tail all check
against the corpus's own gates. **What remains of Tier 2b** is the
axisymmetric isorotating 2-D solve — the (κ, V, c_g, δ, over-spin onset)
five-target test on which the ħ-derivation's benchmark rests.
