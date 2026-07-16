# J1 — the (B.4) 𝔭-grid, executed (App I.1's never-run protocol step)

**ANALYSIS LAYER.** Every claim below is WITHIN-MODEL only: does the corpus's
printed number solve the corpus's printed equations under a stated convention.
Nothing here validates physics/nature. Measured facts only; adjudication is
the coordinator's.

Reproduce: `python3 j1_pgrid.py` (deterministic, no RNG, 53 s; writes
`j1_results.json` stage-flushed, prints the gate table — `j1_run.log`).

## What was run

App I.1 (corpus2/01-GUM-Omega-Paper-v3.0-ext.md, App. I): "Profile
shooting/relaxation: virial ≤ 3×10⁻⁴; ∫b = 1.000; **𝔭-grid over the (B.4)
window**." App B.4: 𝒱 = m̃²(1−σ) + c₂(1−σ)², window **c₂ ∈ [−m̃²/2, m̃²/6)**.
The certified Tier-2b pilot (tier2-closure/radial_solve.py, all gates PASS)
solved c₂ = 0 only (i3_RESULTS.md caveat 4). J1 adds the c₂ term to the SAME
machinery (midpoint energy, graded grid, flow + exact tridiagonal Newton,
ε-dial t(E₂+E₄)/(E₆+E₀) = 0.05 re-tuned per c₂; m = a₆ = a₀ = 1) and scans
9 c₂ points across the window at N = 4000 and 2N = 8000 cells (4N = 16000 at
sentinels −0.5, 0, +0.16). R\*(c₂) is the ε = 0 compacton radius from the
Bogomolny saturation quadrature R\*³ = 3∫₀^π sin²f/√𝒱 df, verified against
closed forms at both checkable points (2^{5/6} at c₂ = 0; (6√2)^{1/3} at
c₂ = −1/2; both to 1×10⁻¹³).

## (a) Grid table (2N values; all solver gates PASS at every point/resolution)

| c₂ | t = a₂ = a₄ | μ_fit | R\*(c₂) | A_fit | **A e^{−μR\*}/2π** (naive) | **A e^{−μR\*}/6** (edge/I3) | /(15π/8) |
|---|---|---|---|---|---|---|---|
| −0.500 | 0.0036324 | 11.7324 | 2.039649 | 1.556×10¹¹ | 1.0024 | 1.0497 | 1.0692 |
| −0.400 | 0.0049331 | 10.0677 | 1.929433 | 1.804×10⁹ | 1.0518 | 1.1015 | 1.1220 |
| −0.300 | 0.0059168 | 9.1927 | 1.876490 | 1.889×10⁸ | 0.9692 | 1.0149 | 1.0338 |
| −0.200 | 0.0067762 | 8.5900 | 1.838111 | 4.085×10⁷ | 0.9031 | 0.9457 | 0.9633 |
| −0.100 | 0.0075634 | 8.1307 | 1.807479 | 1.287×10⁷ | 0.8494 | 0.8895 | 0.9061 |
| **0.000** | 0.0082764 | 7.7726 | 1.781797 | 5.233×10⁶ | **0.8054** | **0.8434** | 0.8591 |
| +0.060 | 0.0086868 | 7.5867 | 1.768126 | 3.290×10⁶ | 0.7821 | 0.8190 | 0.8342 |
| +0.120 | 0.0090844 | 7.4189 | 1.755495 | 2.166×10⁶ | 0.7607 | 0.7966 | 0.8114 |
| +0.160 | 0.0093434 | 7.3153 | 1.747574 | 1.674×10⁶ | 0.7474 | 0.7827 | 0.7972 |

Measured facts about the grid itself:

* **𝔭 is strongly c₂-dependent**: 34.8% spread across the window under every
  amplitude convention; d𝔭/dc₂ ≈ −0.44 per unit c₂ near the benchmark
  (≈ −1.5σ_corpus per 0.1 of c₂). The printed ±0.03 is therefore NOT the
  grid spread over the full window (that would be ±0.16); ±0.03 corresponds
  to a c₂ neighborhood of half-width ≈ 0.07 around a central point.
* **B.5's μ² = m̃²/(2a_ψ) re-confirms across the whole window**, not just at
  c₂ = 0: |μ_fit/μ_derived − 1| ≤ 1×10⁻⁴ gate PASS at all 9 points (μ runs
  7.315 → 11.732); c₂ drops out of the linearization ((1−σ)² = O(f⁴)) and
  enters μ only through the ε-dial's t. Degree ∫b = 1.000 and
  virial ≤ 3×10⁻⁴ (measured ≤ 10⁻⁵ class) PASS everywhere.
* The curve is non-monotone at the low edge (𝔭 peaks at c₂ = −0.4, dips at
  −0.5). This is resolution-stable at N/2N/4N (deltas ≤ 7×10⁻⁵), i.e. a real
  feature of the model at the massless-core edge 𝒱(π) → 0, not numerics.

## (b) Where 𝔭 = 0.84 ± 0.03 lands — the discrimination

In-band c₂ intervals (𝔭 ∈ [0.81, 0.87], linear interpolation on the 2N grid)
and the exact 0.84 crossing:

| convention | 𝔭(c₂=0) | 0.84 crossed at | in-band c₂ interval | contains benchmark c₂ = 0? | μ at crossing vs frozen 7.77255 |
|---|---|---|---|---|---|
| naive A e^{−μR\*}/2π | 0.8054 (−1.15σ) | c₂ = −0.0786 | [−0.138, −0.010] | **NO** | 8.054 (+3.6%) |
| **edge-ref. A e^{−μR\*}/6** | **0.8434 (+0.11σ)** | **c₂ = +0.0083** | **[−0.058, +0.084]** | **YES** | 7.747 (−0.33%) |
| edge-ref. A e^{−μR\*}/(15π/8) | 0.8591 (+0.64σ) | c₂ = +0.0460 | [−0.023, +0.124] | YES | 7.630 (−1.8%) |

Measured facts:

1. **The edge-referenced /6 convention puts the corpus's 0.84 at
   c₂ = +0.008 — the benchmark potential c₂ = 0 to 0.4% of 𝔭** (0.28σ of
   the corpus band), and its in-band interval brackets c₂ = 0 roughly
   symmetrically. Under this convention the corpus's printed 𝔭 is the grid
   value AT the corpus's own calibration point.
2. **The naive /2π convention reaches 0.84 only off-benchmark**
   (c₂ = −0.079); its entire in-band interval excludes c₂ = 0. To hold
   /2π AND 𝔭 = 0.84 one must also move μ to 8.05 — 3.6% above the frozen
   campaign referee μ = 7.77255 = m/√(2a₂), a mismatch ~10³× the μ-fit
   precision (2×10⁻⁶). At the benchmark, /2π stays 1.15σ low, exactly the
   pilot's 0.805.
3. The grid supplies the discrimination I3 asked for (caveat 4): the /6-vs-/2π
   gap (0.038 = 1.26σ) maps to a c₂ offset of 0.086, well resolved by the
   scan. **CONFIRM for the edge-referencing element** in the specific sense:
   the 0.84 is recovered in-window, and it is recovered AT the benchmark only
   under the edge-referenced convention class.
4. Residual degeneracy, printed honestly: /(15π/8) = 0.859 at c₂ = 0 is also
   in-band (+0.64σ) with crossing at c₂ = +0.046. The B.4 grid does NOT
   separate /6 from /(15π/8) (0.7σ apart at benchmark); it separates the
   edge-referenced class from /2π-at-benchmark, i.e. it breaks exactly the
   1.26σ ambiguity I3 flagged, no more.

## (c) Joint 𝔟_eff check (G5b relaxed cap-1200 C_d = 4.261×10⁸ ± 19%)

𝔟 = C_d μ e^{−μR\*}/(F_b 𝔭²) at the c₂ = 0 benchmark (where the two-knot runs
were done), vs 42 ± 6 with the C_d fit error folded in:

| 𝔭 convention | F_b = 2π | F_b = 4π² | F_b = 8π² |
|---|---|---|---|
| naive /2π (0.8054) | 785.7 (+4.9σ) | 125.1 (+3.3σ) | 62.5 (+1.52σ) |
| edge /6 (0.8434) | 716.5 (+4.9σ) | 114.0 (+3.1σ) | **57.0 (+1.19σ)** |
| edge /(15π/8) (0.8591) | 690.6 (+4.8σ) | 109.9 (+3.1σ) | **55.0 (+1.06σ)** |

* **Yes — the same edge-referenced convention keeps 𝔟_eff = 42 ± 6
  recoverable at ≤ 1.2σ with the best-converged C_d**, with one named
  factor F_b = 8π² (matching i3's per-estimator best for G5b, 54.95/+1.06
  reproduced here to 0.02%). The naive-𝔭 row is marginally outside 1.5σ.
* **Convergence caveat, printed honestly (binding):** C_d[G5b cap-1200] is
  NOT a converged measurement — every two-knot run ended at iter_cap, the
  relaxed amplitude moved ×54 from cap-400 to cap-1200 and was still rising,
  ν pinned at fit bounds, dof = 1, ±19% fit error (i3 caveat 2). The
  superseded cap-400 C_d closes instead at F_b = 2 (i3's 41.55/−0.02); the
  factor identity (2 vs 4π² vs 8π²) rides entirely on the unconverged C_d
  systematic and is NOT pinned by J1. What J1 adds is only: the c₂-grid does
  not disturb the joint edge-referenced solution (the 𝔭 entering 𝔟 is the
  benchmark value, which the grid confirms as the corpus's own).

## (d) Grid convergence (decisive numbers, two–three resolutions)

* 𝔭 conventions, max |N − 2N| over the whole grid: ≤ 7.3×10⁻⁵ absolute
  (0.009%, i.e. 0.0024σ of the corpus band). Sentinels at 4N: 𝔭/6 at
  c₂ = 0 is 0.843408 / 0.843379 / 0.843372 (N/2N/4N) — shifts ~10× smaller
  again, O(h²) behavior as certified for the pilot.
* R\*(c₂) quadrature: 2000 vs 4000 Gauss–Legendre nodes agree to ≤ 1.4×10⁻¹³.
* The 0.84-crossing location moves < 10⁻³ in c₂ between N and 2N.
* Tail-fit rms ≤ 1.1×10⁻⁵ (log-space) at N, halving at 2N, at every point;
  fit windows hold 663–3356 points.

## Gate table

All 21 solver runs (9 points × {N, 2N} + 3 × 4N) PASS all seven gates:
virial ≤ 3×10⁻⁴ (corpus I.1 gate), degree K = 1 ± 10⁻³ (∫b = 1.000),
ε-dial 0.05 ± 10⁻³, μ vs B.5 formula ≤ 10⁻⁴, tail rms ≤ 10⁻⁴, containment
f(> 0.97 R_max) ≤ 10⁻⁸, monotonicity. Full table in `j1_run.log` and
per-point booleans in `j1_results.json`.

## Verdict (analysis layer — measured facts; adjudication is the coordinator's)

**The B.4 𝔭-grid CONFIRMS the I3 edge-referencing element.** Executed for the
first time, the printed protocol step yields a strongly c₂-dependent 𝔭
(35% across the window), so the corpus's single 𝔭 = 0.84 ± 0.03 pins a
convention-plus-c₂ pair. Under I3's edge-referenced convention
A e^{−μR\*}/6, the corpus value sits at c₂ = +0.008 ≈ the benchmark
potential (0.28 corpus-σ from the c₂ = 0 value 0.8434); under the
raw/naive /2π normalization it is reachable only at c₂ = −0.079, where the
frozen referee μ is violated by 3.6%, and at the benchmark it remains the
pilot's 1.15σ-low 0.805. The joint 𝔟_eff = 42 ± 6 stays recovered (+1.06 to
+1.19σ) under the same edge convention with G5b's C_d and one named factor
8π² — subject, unchanged, to the C_d non-convergence caveat. Within the
edge-referenced class the /6 vs /(15π/8) residual (0.7σ) remains open.

## Honest caveats

1. **Benchmark identification.** "0.84 is recovered at the corpus's
   calibration point" assumes the ⟨r1⟩ benchmark potential is c₂ = 0. This is
   supported within-model (the certified pilot at c₂ = 0 passes every printed
   I.1 gate and reproduces the corpus-corrected B.5 μ to 2×10⁻⁶) but the
   corpus never prints its benchmark c₂ explicitly.
2. **The frozen μ = 7.77255 used to score the /2π crossing is the campaign's
   own referee value** (radial_results.json), not a corpus-printed number;
   the corpus prints the formula (B.5), not a μ value. The 3.6% mismatch
   statement is a within-campaign consistency fact.
3. **𝔟's factor identity is untouched**: F_b = 8π² closes G5b, F_b = 2 closes
   the superseded cap-400 C_d; the ×54 unconverged C_d systematic spans the
   named-factor family (i3 caveats 1–2 inherited verbatim).
4. **Low-edge behavior**: at c₂ = −0.5 the core vacuum is massless
   (𝒱(π) = 0) and A reaches 1.6×10¹¹ (μR\* = 23.9); gates still pass and the
   non-monotone dip is resolution-stable, but the k₁-normalized A is an
   extreme extrapolation there — the edge-referenced products remain O(1),
   which is itself evidence for the edge-referenced bookkeeping.
5. The 9-point grid is a coarsening of a continuum scan (runtime trivially
   allowed more; 9 points already resolve the band structure to Δc₂ ≈ 0.01
   by interpolation, deltas checked at 3 resolutions).
6. Everything here is WITHIN-MODEL: whether printed numbers solve printed
   equations under stated conventions. No physics/nature validation.

## Files

```
theory-audit/
  j1_pgrid.py       the scan (deterministic, 53 s; PASS/FAIL gate table)
  j1_results.json   per-point records (N/2N/4N), gates, band analysis, joint b_eff
  j1_run.log        console transcript incl. full gate table
  j1_RESULTS.md     this file
```
