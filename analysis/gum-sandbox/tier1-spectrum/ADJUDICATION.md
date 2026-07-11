# Tier 1 — Coordinator's Adjudication

The pilot (`REPORT.md`) reported two literal FAILs against the test spec's
machine-precision expectations. Both are adjudicated here against what the
corpus's theorems **actually claim**. Ruling: **Tier 1 replicates every
claimed result; the two literal failures were over-strict tests in my spec,
not corpus defects.** Filed per the program's print-your-own-defects rule.

## Ruling 1 — D6 / Theorem II.2 (exact masslessness): REPLICATED
Theorem II.2 states the potential sector "contributes **zero mass** to B2 to
all orders" — a statement about the **gap** ω(k→0), not the full dispersion.
Measured: light-branch gap invariant under m_V ∈ {0,1,5,20} to **8×10⁻¹⁵**;
IR speed invariant to 2.7×10⁻⁹. The finite-k deviation (5.3×10⁻³ at k=3,
falling as O(k⁴/A_c); closed form derived in REPORT.md) is the heavy-field
elimination renormalizing the effective stiffness — precisely the paper's own
Theorem V.C.2 behavior ("cone-preserving pieces renormalize the **value** of
c without dispersion" at leading order). My spec demanded pointwise ω(k)
invariance to 1e-12; the theorem never claimed that. **PASS as claimed.**

## Ruling 2 — D7 / Theorem II.3 (tree achirality): REPLICATED (IR)
Measured circular splitting of the light doublet: Δω² ∝ k^4.994 (analytic:
χ₃k⁵/A_c — the χ₃ coupling vanishes identically on the ψ=0 locked wave),
i.e. **exact achirality as k→0** with steep polynomial suppression; the
ψ-carrying heavy doublet splits at O(χ₃k) as it should. The paper's specific
"(ka)²-suppressed regeneration" involves the lattice/Floquet scale absent
from this continuum quadratic model, so the *power* is not comparable — the
*claim* (no tree-level splitting; suppressed regeneration only) is
reproduced. **PASS as claimed, with the scale caveat recorded.**

## Confirmed without qualification
- **The factor 4:** ω₀² = (m_V² + 4μ_c)/J with **d(ω₀²)/dμ_c = 4.000000000**
  measured — the paper's m̃² = m_V² + 4μ̃_c convention confirmed numerically.
- B1 speed = √((λ+2μ)/ρ₀) to 1.3×10⁻¹⁶; mode count 3+3 exact; B4 gap and
  slope as constructed; no negative ω² anywhere.

## The M-1 demonstration — sharper than the pedagogy
Case B (mass on absolute φ, the M-1 defect) produced the defect's true
linearized signature, which **refines** the corpus's teaching-case telling:
- **No gap opens on the light branch** (translation invariance protects a
  gapless transverse wave regardless of the defect — the naive "photon mass"
  reading is subtler than the Course implies).
- Instead the light wave's **rotational content dies**: co-motion fidelity
  r = 4μ_c/(4μ_c + m_V²) — measured 0.9524 / 0.4444 / 0.0476 at m_V = 1/5/20,
  matching the analytic form to 9 digits — and the branch speed shifts
  (c² = μ + μ_c·m_V²/(4μ_c + m_V²), matched to 1e-9).
- The **rotationally-dominant branch is gapped** — that is where the "photon
  mass" lives: the wave that carries E ∝ φ̇ and B ∝ ∇×φ acquires the gap.
Either reading is phenomenologically fatal (light stripped of its EM content,
or EM waves gapped) — so the M-1 defect kills the model as the corpus says,
by a mechanism one level more precise than the corpus states it.

## Net Tier-1 verdict
Four branches with the claimed structure; exact masslessness (gap) at machine
precision; the (m_V² + 4μ_c) convention confirmed; tree achirality exact in
the IR; the M-1 defect exhibited with its mechanism derived in closed form.
**No corpus defect found. Two spec-level over-strictness findings (mine),
adjudicated and filed.**
