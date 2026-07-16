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

## 3. (4.9′) — the censorship ceiling, Q-6′ done properly
- **Old ceiling:** (Δ_lock/c_g)^{3/2} = **6.037×10⁻⁴** (printed "≈6×10⁻⁴" ✓); c_g ± 0.04
  band [5.27, 7.01]×10⁻⁴.
- **Repaired ceiling:** |δ𝔠/𝔠_sat| = s·ε ≤ Δ_lock ⇒ ε ≤ Δ_lock/s = **2.605×10⁻³**
  (printed "≈2.6×10⁻³" ✓; engine-base 2.640×10⁻³; Δ_lock band [2.4, 3.8]×10⁻³ ⇒
  ceiling [2.08, 3.30]×10⁻³).
- **The ×4.3 with provenance:** ceil_new/ceil_old = (Δ/s)/((Δ/c_g)^{3/2}) = **×4.315** —
  the exact number behind the memo-chain's "×4.3".
- **Exponent sensitivity (D3):** reading the same two saturated points at exponent p,
  ceiling(p) is monotone increasing on p ∈ [2/3, 1]: p = 2/3 → 5.95×10⁻⁴ (amplitude
  0.4243 — numerologically ≈ old c_g 0.42±0.04; flagged, not charged), p = 1 → 2.60×10⁻³.
  **The linear reading is the loosest ceiling**; every alternative reading tightens it and
  *widens* the censorship margin.
- **ε_q per reading:** (f_q/𝔪_Sk)², f_q ∈ [0.14, 0.20] GeV, 𝔪_Sk = 1.7 GeV (𝔠-blind) ⇒
  band **[6.78×10⁻³, 1.384×10⁻²]**, central 1.000×10⁻² ("≈10⁻²" ✓).
- **Margins:** old [11.23, 22.93], central 16.57 — "twenty-fold" was fair for v2.0.1.
  Repaired: **[2.60, 5.31], central 3.84**. The v3.0-ext Q-6′ sentence retains
  "twenty-fold" verbatim against its own (4.9′) — **stale (D1)**; it should read
  "≈×4 above (×2.6–5.3)".

## 4. Confinement dichotomy — verdict at every reading
| reading | margin ε_q/ceiling |
|---|---|
| R1 central ε_q / central ceiling | ×3.84 |
| R2 low-edge ε_q (6.78×10⁻³) / central ceiling | ×2.60 |
| R3 printed ε_dress ≈ 10⁻² / central ceiling | ×3.84 |
| R4 low ε_q / Δ_lock-high ceiling (3.8×10⁻³) | ×2.06 |
| R5 low ε_q / engine-base-slope ceiling | ×2.57 |
| R6 low ε_q / exponent-2/3 ceiling | ×11.41 |
| R7 worst joint (Δ hi + engine slope) | **×2.03** |

**Verdict: the dichotomy SURVIVES AT EVERY READING** — all six quarks censored
(ε_q > ceiling, min margin ×2.03), all leptons free (floor 1.5×10⁻⁶ < operative bottom
1.0×10⁻⁵ < ceiling 2.6×10⁻³; termination threshold 8.4×10⁻⁶ untouched; window
[1.5×10⁻⁶, ~3×10⁻³] survives numerically, top now ceiling-set at 86.8% of 3×10⁻³).
**Breach conditions:** a further ×2.60 ceiling rise (central; ×2.03 at the worst joint
reading) — equivalently Δ_lock ≥ 7.8×10⁻³ (λ_field **+6.7σ** above its measured
(3.1±0.7)×10⁻³ω₀), or slope s ≤ 0.442 (×2.6 below +1.15), or σ-inversion low edge
f_q ≤ 0.087 GeV (vs 0.14). No reading is within a factor 2 of breach.

## 5. Cross-check vs i2_propagate.py (task c)
All 11 overlaps AGREE (rel ≤ 10⁻⁴, most ≤ 10⁻⁶): ceiling old/new, censorship margins
old/new (both edges), anchor pulls 0.0156σ / 11.283σ, operative-range top, f-window top
59.87 MeV, slope "+1.15". **No discrepancy found.** One print-vs-computed note (D6):
computed f-window bottom 3.71 MeV vs F-A15-3's printed "[4, 29] MeV" — corpus rounding,
substance identical (pre-existing, shared with i2).

## 6. Defects (printed honestly; includes our own)
- **D1 (corpus v3.0-ext, Q-6′):** stale "twenty-fold above the entrainment ceiling"
  contradicts its own (4.9′); consistent figure ×2.6–5.3 (central ≈×3.8).
- **D2 (corpus v3.0-ext, IV.J′):** retains λ_e ∼ ε_e^{1/3} ∼ 0.01–0.1 — the superseded
  restricted-branch shape law inside the repaired section.
- **D3 (the repaired law / this memo):** "+1.15ε" linear class rests on TWO ε points;
  exponent unmeasured. Dichotomy is exponent-robust (§3), but the ×4.3 loosening — and
  the downstream 60 MeV f-top and 7.8×10⁻¹⁰ Majoron edge — are linear-reading-specific;
  a 2/3 reading reverts the ceiling to ≈6×10⁻⁴ (amplitude 0.424 ≈ old c_g: numerology flag).
- **D4 (engine):** axi4 ring_t0 base 3.203459 vs exact 3.201125 (+7.3×10⁻⁴, ×5 its own
  N-ladder estimate 1.4×10⁻⁴) ⇒ 1.4% slope-reading spread (1.152 vs 1.136), ceiling
  2.60 vs 2.64×10⁻³. Not material to any verdict.
- **D5 (scope):** J2 discharges the (4.8′)/(4.9′) *arithmetic*; the multi-ε saturated
  re-scan and the λ* shape re-scan remain OPEN obligations.
- **D6:** f-window bottom 3.71 vs printed "4" MeV (rounding).

## 7. Verdict (analysis layer)
**(4.8′) re-derived:** on the saturated branch only 𝔠 (and g_tot) run with ε —
𝔠_sat(ε) = 3.201125·(1 + 1.152ε) leading order (sign flipped, exponent 2/3 → 1, both
traceable to the fixed-oblateness minimiser); κ = 1/√2, depth 29.29%, ω_th = √2ω₀,
E_rot/E = ¼ are ε-frozen exactly, and the one-sided physical-constant bands flip to the
+𝔠 side with zero band on κ/depth/ω_th. The ⟨r1⟩ ε-anchor is confirmed FORFEITED
(repaired prediction 11.283σ from ⟨r1⟩, correctly — restricted-branch code validation);
the repaired branch has no measured scan (2 points; the re-scan proper stays open).
**(4.9′) re-derived:** ceiling Δ_lock/s = 2.605×10⁻³ (×4.315 above the old 6.037×10⁻⁴),
Q-6′ margins thin to [2.60, 5.31] (central 3.84), and the confinement dichotomy —
all six quarks censored, all leptons free — **survives at every reading examined**
(worst joint ×2.03); breach requires a further ×2.6 ceiling rise ≈ a +6.7σ λ_field.
The v3.0-ext Q-6′ "twenty-fold" sentence and the IV.J′ ε^{1/3} shape sentence are stale
(D1/D2) and should be amended in the next revision pass. All i2 overlaps agree; 44/44
gates PASS. Within-model only; the coordinator adjudicates.

*Files: j2_epsscan.py (44/44), j2_results.json (stage-flushed), this memo.*
