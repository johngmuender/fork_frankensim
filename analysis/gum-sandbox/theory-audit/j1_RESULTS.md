# J1 — the (B.4) 𝔭-grid, executed (App I.1's never-run protocol step)

**ANALYSIS LAYER.** Every claim below is WITHIN-MODEL only: does the corpus's
printed number solve the corpus's printed equations under a stated convention.
Nothing here validates physics/nature.

Reproduce: `python3 j1_pgrid.py` (deterministic, no RNG; writes
`j1_results.json` stage-by-stage and prints the gate table).
Companion files: `j1_pgrid.py`, `j1_results.json`.

STATUS: RUNNING — sections below fill in as stages land.

## What was run

App I.1 (corpus2/01-GUM-Omega-Paper-v3.0-ext.md, App. I): "Profile
shooting/relaxation: virial ≤ 3×10⁻⁴; ∫b = 1.000; **𝔭-grid over the (B.4)
window**." App B.4: 𝒱 = m̃²(1−σ) + c₂(1−σ)², positivity + saturability window
**c₂ ∈ [−m̃²/2, m̃²/6)**. The Tier-2b pilot (tier2-closure/radial_solve.py,
all gates PASS) solved c₂ = 0 only (i3_RESULTS.md caveat 4). J1 extends the
certified machinery with the c₂ term and scans the window.

* Solver: identical midpoint/graded-grid energy minimization (flow + exact
  tridiagonal Newton), ε-dial t(E₂+E₄)/(E₆+E₀) = 0.05 re-tuned per c₂.
* E₀ generalized: ∫ [m²(1−cos f) + c₂(1−cos f)²] r² dr (m = a₆ = a₀ = 1).
* R\*(c₂): ε = 0 compacton radius from the Bogomolny saturation quadrature
  R\*³ = 3∫₀^π sin²f /√𝒱 df (reproduces 2^{5/6} at c₂ = 0).
* Conventions extracted at each grid point (μ_fit, A_fit from the tail fit):
  raw/naive **A e^{−μR\*}/2π** (pilot's nearest, 0.805 at c₂ = 0) and I3's
  edge-referenced **A e^{−μR\*}/6** (0.843 at c₂ = 0); also /(15π/8), the
  literal A/μ, and f_fit(R\*) for the record.

## (a) Grid table

(pending)

## (b) Where 0.84 ± 0.03 lands

(pending)

## (c) Joint 𝔟_eff check (G5b relaxed cap-1200 C_d)

(pending)

## (d) Grid convergence

(pending)

## Gate table

(pending)

## Verdict (analysis layer)

(pending)

## Honest caveats

(pending)
