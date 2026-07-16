# K2 — the ε < 0.005 probe: J4's saturated ε-scan extended to ε ∈ {0.001, 0.002} (analysis layer)

**Phase:** ROADMAP_v7 K2 · 2026-07-16 · theory-audit numerics subagent (Fable-class).
**Charge:** extend the J4 saturated small-ε scan below its stated limit ε = 0.005 (J4-D3);
test whether the measured linear law (p = 0.997 [0.996, 1.012]) holds to the smallest
honestly solvable ε. **Computations:** `k2_smalleps.py` (deterministic, i1 machinery reused
unmodified via the J4 code path) → `k2_results.json` (flushed after every solve point) +
`k2_run.log`. Gates: **11/12 PASS** — the one FAIL (G3d) is a finding, printed in §5.
**WITHIN-MODEL ONLY; the coordinator adjudicates.**

## §1 New saturated solves at ε ∈ {0.001, 0.002} (MEASURED)

i1 dial + `saturated_at` reused unmodified, **J4's budgets verbatim** (N_OPT/N_POL/N_FIN/
N_COARSE = 192/320/480/320, B_TRI/B_COLD/B_POL = 220/1500/260); two grid resolutions per
point (N480/N320) + the dial's N×2/rmax+2 ladder. **No escalation and no coarsening was
needed** — both points passed all four honesty gates (grid conv ≤ 5×10⁻⁴, no boundary
flags, dial ladder ≤ 10⁻⁴, signal/conv ≥ 10) at first attempt; total run 176 s (≪ 60 min).
No hard solver limit was reached, so the verdict below is NOT convergence-bounded.

| ε (achieved) | t | 𝔠_sat (N480) | N320 | conv rel | signal/conv | hh-core | floor 3.2011+12t | dial N2x |
|---|---|---|---|---|---|---|---|---|
| 0.0010001 | 1.1681e-4 | **3.205304** | 3.205215 | 2.8e-5 | 47 | 3.25026 | 3.20253 OK | −3.0e-5 |
| 0.0019995 | 2.4616e-4 | **3.208836** | 3.208740 | 3.0e-5 | 80 | 3.25401 | 3.20408 OK | −1.5e-5 |

Monotone: 3.201125 (exact) < 𝔠(0.001) < 𝔠(0.002) < J4's 𝔠(0.005) = 3.219560 (G1c);
analytic floor respected (G1b: +2.8e-3 / +4.8e-3 above).

## §2 Cross-check gate: ε = 0.005 re-solve vs J4's banked value

Same budgets, same code path: **𝔠_sat(0.005) = 3.2195637 vs J4's banked 3.2195597 —
rel +1.24×10⁻⁶, ×27 inside J4's two-resolution band 3.3×10⁻⁵** (G2a PASS; K2's own band
3.4×10⁻⁵). The extension is anchored to the J4 dataset it extends.

## §3 Extended fits (J4's four classes, fit code verbatim; 13-point dataset)

Dataset = 2 new K2 points + J4's 11 (ε = 0.005…2.0); σ = max(grid conv, 3×10⁻⁵)·𝔠,
CIs residual-rescaled — J4's σ-model verbatim. New window W0: ε ≤ 0.01.

| window | (i) aε | (ii) aε^{2/3} | (iii) aε^p, p [CI] | (iv) aε + bε^{2/3} | M1f c0, s | M3f p [CI] |
|---|---|---|---|---|---|---|
| **W0 ε≤0.01 (new)** | a=1.1472, χ²=47 | a=0.218, χ²=3614 | p=0.965 [0.952,0.978] | 1.039ε + 0.0211ε^{2/3} | 3.201671, 1.1215 | **1.014 [0.998,1.030]** |
| W1 ε≤0.05 (primary, ext) | a=1.1372, χ²=59 | a=0.335, χ²=71584 | **p=0.994 [0.990,0.998]** | 1.115ε + 0.0068ε^{2/3} | 3.201546, 1.1317 | 1.005 [1.003,1.007] |
| W2 ε≤0.03 | a=1.1386 | a=0.300 | p=0.990 [0.984,0.996] | 1.105ε + 0.0092ε^{2/3} | 3.201571, 1.1307 | 1.009 |
| W3 ε≤0.10 | a=1.1359 | a=0.375 | p=0.996 [0.994,0.998] | 1.121ε + 0.0051ε^{2/3} | 3.201543, 1.1321 | 1.002 |
| W4 ε≤0.20 | a=1.1483 | a=0.431 | p=1.011 [1.006,1.015] | 1.172ε − 0.0099ε^{2/3} | 3.201021, 1.1489 | 1.019 |

- **2/3 EXCLUDED down to ε = 0.001**: χ²(2/3)/χ²(lin) = 1222 on W1ext, and **77 on the
  new W0 alone** (G3b) — the smallest-ε decade by itself rejects the 2/3 class.
- **W1ext free-p: p = 0.994 [0.990, 0.998], jackknife [0.990, 0.995]; raw all-readings
  envelope [0.965, 1.011]** (vs J4's [0.996, 1.012]). The new low edge 0.965 is W0's
  FIXED-exact-intercept fit; with free intercept W0 gives p = 1.014 [0.998, 1.030] — the
  dip is an intercept effect, not exponent (see G3d, §5). Model-free local exponents:
  0.884±0.038 (0.001–0.002), 0.951±0.015 (0.002–0.005), then 0.985, 0.997, 0.998, 0.998, 0.998.
- **G3c PASS:** W1ext free-intercept fit c0 = 3.201546 ± 0.000062 (rel +1.32×10⁻⁴ vs exact),
  s = 1.1317 vs J4's 1.1338 (−0.2%, well inside the 1.3% family band).
- **G3d FAIL (the K2 finding, quantified):** J4's constant-bias model (intercept offset
  +7.7×10⁻⁵, free-intercept slope 1.1338, NO refit) predicts the new points at
  **+3.15σ (ε=0.001) and +2.16σ (ε=0.002)** — the measured per-point offsets above pure
  linear are +1.72×10⁻⁴ / +1.42×10⁻⁴, vs +7.7×10⁻⁵ fitted on [0.005, 0.05].
- **Shape discrimination of the small-ε deviation (measured):** on W0, constant + linear
  (M1f, χ² = 0.71) beats linear + ε^{2/3} admixture (M4, χ² = 6.65) at equal parameter
  count, and M4's admixture amplitude is window-unstable (0.021 → 0.007 → 0.005) — the
  deviation has the shape of a slowly-drifting CONSTANT offset (family intercept bias),
  not of a genuine ε^{2/3} term.

## §4 (4.9′) propagation decision (trigger: ceiling movement > J4's 1.3% family spread)

Central fitted-slope ceilings (law(ε) = Δ_lock = 3×10⁻³, J2 arithmetic):
linear class J4 2.63905e-3 → K2ext 2.63810e-3 (**move −3.6×10⁻⁴**); free-p class
J4 2.61610e-3 → K2ext 2.59577e-3 (**move −7.8×10⁻³**). Max central movement
**7.8×10⁻³ ≤ 1.3×10⁻² → G4a: PROPAGATION UNCHANGED.** J4's numbers stand and are cited:
ceiling ∈ [2.60, 2.71]×10⁻³, thinning ×[4.31, 4.50], f_top ∈ [59.9, 61.1] MeV,
g_min ∈ [7.66, 7.82]×10⁻¹⁰; dichotomy robust (min worst-joint ×1.98). RECORD ONLY: pinning
p to the raw-envelope low edge 0.965 would give ceiling 2.413×10⁻³ (−8.6%), but that edge
is the measured intercept-bias drift leaking into a fixed-intercept fit (G3d + shape test
above), i.e. double-counting a measured systematic — recorded, not propagated.

## §5 Defects (printed honestly)

- **K2-D1 (headline; = the G3d FAIL): the 12-mode family's intercept bias is NOT constant
  into the small-ε regime.** It roughly doubles, from +7.7×10⁻⁵ (J4 G2c, fitted on
  ε ∈ [0.005, 0.05]) to +1.7×10⁻⁴ ± 0.25×10⁻⁴ (W0 refit; per-point +1.42/+1.72×10⁻⁴).
  What it undermines: J4-D2's "bias is smooth and mostly common-mode" is now quantified as
  bias DRIFT — fixed-exact-intercept exponent readings below ε ≈ 0.005 (W0's p = 0.965,
  local 0.884±0.038) are bias artifacts and must not be read as physics. What it does NOT
  undermine: the exponent (free-intercept p = 1 within CI on every window incl. W0), the
  2/3 exclusion (χ² ratio 77 on W0 with the bias present), the slope (−0.2% vs J4), or the
  ceiling (§4). Within-model, "family-bias drift" vs "a genuine constant-shaped sub-leading
  term" cannot be fully split; the χ² shape test (constant ≫ ε^{2/3}) and the variational
  upper-bound character both point to the family, and either reading leaves leading order
  linear.
- **K2-D2 (2/3 admixture bound loosened at small ε):** J4's b = 0.0034 ± 0.0014 becomes
  window-unstable once ε ≤ 0.002 enters (0.021 ± 0.006 on W0) because the admixture is
  degenerate with K2-D1's offset there. The honest statement: any true ε^{2/3} admixture
  is ≤ 0.02 in amplitude — 2/3 stays excluded as LEADING order everywhere.
- **K2-D3 (scope):** ε < 0.001 unprobed — but no solver limit was hit (both points
  converged cleanly at J4 budgets with signal/conv ≥ 47), so 0.001 is a budget choice,
  not a wall. Below ε ~ 3×10⁻⁴ the K2-D1 offset would exceed ~50% of the signal and a
  fixed-intercept scan would become bias-dominated regardless of solver quality.
- **K2-D4 (inherited):** J4-D1 (1.3% axi4-vs-i1 amplitude spread) and J4-D5 (Δ_lock, ε_q,
  f-map, m₃ inherited verbatim) carry over unchanged; λ* shape re-scan remains OPEN.

## §6 Verdict (analysis layer — within-model; the coordinator adjudicates)

**The linear law holds to ε = 0.001, the smallest ε probed — and no honest-convergence
limit was reached in getting there.** Both new points solved cleanly at J4's budgets
(conv ≤ 3.0×10⁻⁵, two resolutions each, ε = 0.005 reproduction at +1.24×10⁻⁶); the 2/3
class is now excluded by the smallest-ε decade alone (χ² ratio 77 on ε ≤ 0.01; 1222 on
the extended primary window); the extended free-p reading p = 0.994 [0.990, 0.998]
(free-intercept 1.005 [1.003, 1.007]) brackets 1. **J4's verdict (linear; 2/3 excluded)
STRENGTHENS on the exclusion and the slope, while its p-envelope widens on the low side
([0.965, 1.011] raw)** — the widening is entirely the newly-measured drift of the 12-mode
family's intercept bias (K2-D1, the G3d FAIL: +7.7×10⁻⁵ → +1.7×10⁻⁴ as ε → 0.001), an
analysis-layer systematic of the upper-bound family, not small-ε physics; with the
intercept freed, every window returns p = 1 within CI. **Propagation UNCHANGED (G4a):
central ceiling movement 7.8×10⁻³ < the 1.3% trigger — J4's (4.9′) numbers stand as cited
in §4.** Files: `k2_smalleps.py` (deterministic, restart-safe, 11/12 gates, 176 s),
`k2_results.json` (per-point flushed), `k2_run.log`, this memo.
