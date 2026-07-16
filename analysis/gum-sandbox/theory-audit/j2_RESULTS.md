# J2 RESULTS — repaired (4.8′)/(4.9′) ε-scan arithmetic (analysis layer)

**Status: COMPLETE — 44/44 gates PASS** (`j2_epsscan.py` → `j2_results.json`, stage-flushed).
All claims WITHIN-MODEL only; nothing here bears on nature. Workstream J2 (first
arithmetic obligation of the I2 memo). Sources: i2_RESULTS.md/i2_propagate.py;
T3_repaired_closure.md/t3_closures.py; tier2-closure/axi4_locked_results.json
(the only saturated-branch ε data); the two corpus editions (pins below).

## 0. Provenance pins (11 verbatim text pins, all PASS)
- **(4.8) v2.0.1:** 𝔠(ε) = 𝔠₀[1 − c_g ε^{2/3}(1+O(ε^{1/3}lnε))], 𝔠₀ = 128√42/(105π) = 2.5147,
  c_g = 0.42 ± 0.04; scan domain ε ∈ [10⁻³, 0.1]; anchor 2.515·[1−0.42·0.05^{2/3}] = 2.371.
- **(4.9) v2.0.1:** ceiling ε_e ≤ (Δ_lock/c_g)^{3/2} ≈ 6×10⁻⁴; floor 1.5×10⁻⁶; window
  [1×10⁻⁵, 3×10⁻³]₆₈; Δ_lock = λ_field/ω ≈ 3×10⁻³ (IV.H.3; λ_field = (3.1±0.7)×10⁻³ω₀).
- **(4.8′) v3.0-ext:** "3.2011 → 3.3855 at ε = 0.05 (+5.76%; ≈ +1.15ε class; the rung's law
  is −1.03ε^{1.0})"; **(4.9′):** "ceiling ≈ 2.6×10⁻³", window "[1.5×10⁻⁶, ~3×10⁻³]".
- **Q-6′ v3.0-ext:** ε_dress = (f_q/𝔪_Sk)² ≈ 10⁻², "twenty-fold above the entrainment
  ceiling" — retained verbatim (defect D1); **V15.7:** ε_q = (0.14–0.20/1.7)² = 0.7–1.4×10⁻²;
  **Q-5:** f_q ≈ 0.14–0.20 GeV; **F-A15-3:** old operative range ε_e ∈ [1.0×10⁻⁵, 6×10⁻⁴]₆₈.
- **B.6 record:** "saturated: +1.15ε class"; **Repaired-Closure v1.0:** "margin thins ×4.3",
  obligation "(4.8′)/(4.9′)". Exact constants: 𝔠₀ = 2.5147536, 𝔠_sat = 64√2/(9π) = 3.2011247.

## 1. (4.8′) — ε-scan under the saturated closure (MEASURED FACTS)
The saturated branch (fixed oblate compacton + pinned marginal halo) has exactly
**two** measured ε points (axi4 ring-saturated): 𝔠(t→0) = 3.203459 (engine; exact
𝔠_sat = 3.2011247, bias +7.3×10⁻⁴ = D4) and 𝔠(0.05) = 3.385470.

- **The law:** 𝔠_sat(ε) = 3.201125·(1 + s·ε) leading order, **s = +1.152** (exact-base;
  engine-base +1.136, spread 1.4% = D4). Rise at ε = 0.05: **+5.76%** — the printed (4.8′)
  values reproduced exactly. **Sign flips** (falling → rising) and the **exponent moves
  2/3 → 1**: the old fractional power was the λ*²∝ε^{2/3} pancake-geometry narrative; the
  saturated minimiser's oblateness is fixed O(1), so the leading dependence is analytic.
  Linear class is measured at only two points — exponent unmeasured (D3).
- **Which constants run:** only 𝔠 (and g_tot: engine +2.84% at ε = 0.05, ≈ +0.57ε; hence
  **κ²g_tot(0.05) ≈ 1.5008 vs the S4′ stake 35/24 = 1.45833** — benches at ε ≈ 0.05 should
  expect ≈1.50; the discriminating stake value is the ε→0 limit).
- **Which constants are ε-FROZEN exactly** (saturation/T3.1's KKT condition is ε-free;
  engine bit-identical at both points): κ = 1/√2, binding depth 29.29%, ω_th = √2ω₀,
  E_rot/E = ¼. Under (4.8) these ran with ε (⟨r1⟩'s 0.802 / 19.8% at ε = 0.05).
  The ε = 0.05 cap-ladder (𝔠 3.05→3.29 rising, κ 0.752→0.722 falling, no interior optimum)
  confirms the ring channel still binds at ε > 0.
- **Repaired one-sided bands (F-A15-1 redo, side flips):** drift strictly positive ⇒
  𝔠_phys = 3.2011⁺⁰·⁰⁰⁹⁶₋₀ (i.e. [3.2011, 3.2107] at Δ_lock); κ/depth/ω_th: **zero ε-band**
  (exact). Old bands were negative-side (2.515₋₀.₀₀₈⁺⁰).
- **Scan table over the printed domain [10⁻³, 0.1]** (old falling vs new rising):
  ε = 10⁻³: 2.5042 | 3.2048; ε = 10⁻²: 2.4657 | 3.2380; ε = 0.05: 2.3714 | **3.3855 (MEASURED)**;
  ε = 0.1: 2.2876 | 3.5698. All non-0.05 saturated rows are interp/extrap grade (D5) —
  the multi-ε engine re-scan itself remains OPEN; J2 discharges the arithmetic only.

## 2. The forfeited ε-anchor: survive / shift / undefined census of (4.8)'s numbers
| printed item | verdict |
|---|---|
| 𝔠₀ = 2.5147 (intercept) | STRUCK → 𝔠_sat = 3.201125 [T3.1] |
| c_g = 0.42 ± 0.04 (amplitude) | UNDEFINED on repaired branch; → s ≈ +1.15 (two-point, no band yet) |
| exponent 2/3 + λ*² ∝ ε^{2/3} narrative | UNDEFINED (fixed-oblateness minimiser ⇒ linear) |
| sign (falling) | FLIPS (rising) |
| anchor 2.371 vs ⟨r1⟩ 2.37±0.09 (0.016σ) | arithmetic SURVIVES as restricted-branch code-validation record; **FORFEITED as corroboration** — repaired 𝔠(0.05) = 3.38547 is **11.283σ** from ⟨r1⟩, correctly (F-R5/F-R15); repaired branch has NO measured scan |
| 3-way loop δ𝔠/𝔠₀ = 0.057 ⇒ g* = 1.329 vs 1.31±0.04 | SURVIVES as restricted-branch record (B.6); UNDEFINED on saturated branch (no ε-running λ*) |
| κ(ε), depth(ε), ω_th(ε) running + one-sided bands | SHIFT to exact ε-frozen constants; bands collapse to zero (only 𝔠 keeps a +side band) |
| λ*(ε) = 0.42ε^{1/3}; λ*(ε_q) = 0.246 | SUPERSEDED (V15.7 row); shape re-scan obligation stays OPEN |
| ⟨r1⟩ "benchmark-at-ε = 0.05" | SURVIVES with restricted-branch/transit label |
| domain ε ∈ [10⁻³, 0.1] | SURVIVES as domain; saturated branch measured at 2 points only (D5) |
| λ_e ∼ ε_e^{1/3} ∼ 0.01–0.1 (IV.J consequences) | UNDEFINED/superseded — **still printed in v3.0-ext IV.J′** (D2) |
| ε-suppression leg of point-likeness (≤6×10⁻⁴) | SHIFTS to ≤2.6×10⁻³ (×4.3 weaker, still small) |

## 3. (4.9′) — censorship ceiling, Q-6′ done properly
(pending)

## 4. Confinement dichotomy verdict per reading
(pending)

## 5. Cross-check vs i2_propagate.py
(pending)

## 6. Defects
(pending)

## 7. Verdict (analysis layer)
(pending)
