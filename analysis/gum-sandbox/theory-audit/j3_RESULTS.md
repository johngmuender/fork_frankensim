# J3 — MAJORON BATTERY RE-RUN AT THE PROPAGATED LOW EDGE g ≈ 7.8×10⁻¹⁰

**Phase:** ROADMAP_v6 J3 · 2026-07-16 · theory-audit numerics subagent (Fable-class).
**Charge:** the second arithmetic obligation filed by I2 (corpus2/10 §Obligations: "the
Majoron battery re-run at g ≈ 7.8×10⁻¹⁰"): re-run the corpus's OWN battery arithmetic on
the repaired phason window f ∈ [3.7, ≈60] MeV, old-window control first.
**Computations:** `j3_majoron.py` (deterministic, standalone, stdlib-only; **17/17 gates
PASS**) → `j3_results.json` (stage-flushed).
**Standing notice — ANALYSIS LAYER, WITHIN-MODEL ONLY.** This re-executes the corpus's
printed arithmetic at shifted inputs; it validates nothing about nature.

## §0 The battery as the corpus actually prints it (source census — measured fact)

There is **no row-by-row numerical battery table anywhere in the archive.** App K.3 says
"Majoron dictionary g = m_ν/f with the SN/BBN/streaming battery table," but the table's
content is only ever printed as four stamped rows (Ω VII.I, both v2.0.1 and v3.0-ext,
verbatim identical; Substrate Course 19.3 "imports that literature's constraint battery
wholesale and passes it"):
**R1** free-streaming: f ≥ B-ν1 = (T·M_Pl·m_ν²)^{1/4} = 1.57 → **1.4 MeV** "with thermal
factors" (the only row with a printed formula+number); inverted: floor B-ν1′ ε_e ≥ 1.5×10⁻⁶.
**R2** BBN ✓ — quantitative anchor ΔN_eff = (4/7)(10.75/106.75)^{4/3} = 0.0268 [CJ-thermal];
no BBN bound value printed. **R3** "SN-cooling band evaded from below ✓" — **band edges
never printed**. **R4** "late recoupling harmless ✓" — **no formula or number printed**.
Perimeter pricing (WS-K2 b/c): static PV/spin forces "at g² ≲ 10⁻¹⁶ — negligible";
Cs-PNC folding distortion "(λZ/a₀)² = 10⁻⁴–3×10⁻³ per unit admixture x" over f ∈ [4, 29].
Window versions on record: **body dictionary** (VII.I/II.H): f ∈ [4, 60] MeV, g ∈ [8×10⁻¹⁰,
1.3×10⁻⁸]; **F-A15-3 operative** (window ∩ old ceiling 6×10⁻⁴): f ∈ [4, 29] MeV, g ∈
[1.7×10⁻⁹, 1.3×10⁻⁸] "(battery re-run: passes throughout)"; worksheets (WS-K row 9, WS-T
row 6, Course 19.3) carry [4, 29]; WS-nu-P4 r1 re-certified the V15.5 arithmetic only.

## §1 Control: the old window (f ∈ [3.7, 29] MeV) — corpus prints reproduced (A1–A8)

| item | recomputed | corpus print | gate |
|---|---|---|---|
| f(ε=1×10⁻⁵) = 0.69√ε·1.7 GeV | 3.7094 MeV | "4" (generous round; i2 "3.7") | A1 PASS |
| f(ε=6×10⁻⁴) | 28.73 MeV | "29" | A2 PASS |
| g low edge, m₃=0.0468/28.7 | 1.63×10⁻⁹ | i2 "1.6×10⁻⁹" | A3 PASS |
| g endpoints, m=0.0503 (osc. floor)/[29, 4] | 1.73×10⁻⁹, 1.26×10⁻⁸ | F-A15-3 "[1.7×10⁻⁹, 1.3×10⁻⁸]" | A3 PASS |
| ΔN_eff | 0.026772 | 0.0268 | A4 PASS |
| ε-floor = (1.4/(0.69·1700))² | 1.424×10⁻⁶ | "1.4 → 1.5×10⁻⁶" | A5 PASS |
| B-ν1 = 1.57 MeV needs T = 0.227 eV | recombination-class | convention-limited (tier5 C5c) | A6 PASS |
| body low edge 0.047 eV/60 MeV | 7.83×10⁻¹⁰ | "8×10⁻¹⁰" | A7 PASS |

Old-window battery verdict reproduced (A8): R1 margin f_lo/B-ν1 = 3.71/1.4 = **2.65** at
central m₃ (worst m-band corner m=0.115: B-ν1∝√m → 2.19 MeV, margin **1.69**; raw-1.57
scaling 1.51); ε-floor slack (max 3.5×10⁻⁶ < 1×10⁻⁵); R2 ✓; R3/R4 carry the corpus's
stamp with bounds unprinted. **Control = corpus's "battery passed": reproduced.**

## §2 Re-run: repaired window (f ∈ [3.7, ≈60] MeV, g → 7.8×10⁻¹⁰), m_ν swept (B1–B9)

Repaired ceiling 3×10⁻³/1.152 = 2.604×10⁻³ → f_top = **59.86 MeV** (coincident with the
II.H window top "60"); new low edge g = 0.0468 eV/59.86 MeV = **7.82×10⁻¹⁰** — the I2
obligation number reproduces exactly (B1, B2). Full sweep: 121-point f-grid × 25-point
m-grid over the band [0.019, 0.115] eV (plus live window [0.050, 0.057] in the JSON);
all extrema verified at corners, no interior reversal (B8). g spans [3.17×10⁻¹⁰ (m=0.019,
f=60), 3.10×10⁻⁸ (m=0.115, f=3.71)].

| row | old window | repaired window | verdict |
|---|---|---|---|
| R1 free-streaming f ≥ B-ν1(m) | margin 2.65 (m₃) / 1.69 (m=0.115) | **bit-identical** — binding edge f_lo unchanged | **PASS, margin unchanged** (B3) |
| R1b ε-floor B-ν1′ | 1.4×10⁻⁶ (m₃), 3.5×10⁻⁶ (m=0.115) | bit-identical, < operative bottom 1×10⁻⁵ | **PASS** (B4) |
| R2 BBN / ΔN_eff | 0.026772 | bit-identical (f,g-independent) | **PASS** (B5) |
| R3 SN-cooling evaded from below | binding edge g_max = m/f_lo: 1.26×10⁻⁸ (m₃), 3.10×10⁻⁸ (m=0.115) | **g_max unchanged for every m**; g_min drops 1.6×10⁻⁹ → 7.8×10⁻¹⁰ (band corner 6.6×10⁻¹⁰ → 3.2×10⁻¹⁰) — evasion deepens | **INHERITED-PASS** + spec-recovery (B6) |
| R4 late recoupling harmless | binds at g_max | g_max unchanged; new strip monotonically safer | **INHERITED-PASS** + spec-recovery (B7) |
| lab perimeter (WS-K2 b/c) | g² ≤ 1.6×10⁻¹⁶ (dictionary endpoints); PNC per-x [5×10⁻⁵, 3×10⁻³], λ 6.9–53 fm | g² top unchanged; PNC per-x extends DOWN to 1.2×10⁻⁵ (λ 3.3–53 fm) — compliance improves | **PASS, margin-shift favorable** (B9) |

**The one physics direction that could have bitten — checked:** a row binding from above
on f (or below on g) would newly engage at the 60-MeV edge. The printed battery contains
none: every coupling-driven row binds at g_max = m_ν/f_lo, and f_lo (= operative ε bottom
1×10⁻⁵, repair-invariant per F-A15-3/i2) does not move. The newly opened strip is exactly
g ∈ [7.8×10⁻¹⁰, 1.6×10⁻⁹) — monotonically deeper evasion on every printed constraint.
Signature rows (not gates): P-ν4 "Majoron-mode 0νββ at g ∼ 10⁻⁹" — the new strip sits
at/below 10⁻⁹, so the far-horizon signature retreats further at the new edge; ΔN_eff
[CJ-thermal] conditionality: smaller g makes primordial thermalization less likely (grade
unchanged; a null does not kill).

## §3 Spec-recovery limits and defects (printed, per discipline)

1. **R3/R4 pass by inheritance, not recomputation.** The SN-cooling band edges and the
   recoupling bound are nowhere printed in the archive; "PASS" for those rows means: the
   corpus stamped ✓ on a domain whose binding edge (g_max) is bit-identical before and
   after the repair, and the repair moves the free edge away. If the corpus's imported
   "wholesale literature battery" contains any row binding at high f or low g, it was
   never printed and is outside this re-run's reach.
2. **B-ν1 is convention-limited** (tier5 C5c standing): 1.57 MeV recovers with T = 0.227 eV
   (full M_Pl); the → 1.4 MeV "thermal factors" (×0.89) are unprinted. J3 uses the printed
   1.4 (and 1.57 as cross-check), scaled by the formula's own ∝ √m_ν.
3. **m-convention mixing in the printed endpoints** (tier5 C5d standing): low edges use
   m₃ ≈ 0.047, F-A15-3's 1.7×10⁻⁹ needs m = 0.0503; the band-top corner (m = 0.115,
   f = 3.71) gives g = 3.1×10⁻⁸, g² = 9.6×10⁻¹⁶ — ×6 above WS-K2(b)'s "g² ≲ 10⁻¹⁶" class.
   Not repair-driven (that corner exists identically in the old window): the corpus prices
   its perimeter at dictionary-central m, not at the m-band corner. Filed as an
   observation on the corpus's print, not a J3 failure.
4. The corpus's own body dictionary (VII.I both versions) already stamped all four ✓ on
   f ∈ [4, 60] — i.e. on the II.H admission window BEFORE intersecting the old ceiling.
   The repair makes the operative window (ceiling 2.6×10⁻³) coincide with that domain, so
   "battery passed on [4, 60]" was already a corpus claim; J3's contribution is the control
   + edge-to-edge arithmetic showing nothing in the printed battery distinguishes the two.
5. f(1×10⁻⁵) = 3.71 MeV printed as "4" is a generous round (tier5 C7b standing); i2's
   "3.7" is the honest value; J3 uses 3.7094.

## §4 Verdict (analysis layer, within-model)

**I2's "SHIFTS-HARMLESSLY" is CONFIRMED at the battery level: no row newly binds or fails
anywhere in the repaired window f ∈ [3.7, 59.9] MeV, across the full m_ν band [0.019,
0.115] eV and the live window [0.050, 0.057] eV.** Old-window control reproduces every
printed number and the "passed" verdict (A1–A8); the repaired re-run leaves R1/R1b/R2
bit-identical (their binding edges are repair-invariant), deepens R3/R4's evasion
monotonically (inherited-PASS; bounds unprinted — spec-recovery limits 1), and improves
the lab-perimeter margins. Free-streaming margin at the worst swept corner: **1.69**
(thermal convention) / 1.51 (raw), identical old vs new — the only margin under 2×
anywhere in the battery, and it is m-driven, not repair-driven. The corpus's "battery
passed" claim survives the T3.1 propagation over the whole repaired window — with the
honest caveat that two of its four rows were only ever stamps, and stamps can only be
inherited, not re-run.

*Files: `j3_majoron.py` (17/17 gates), `j3_results.json`, this memo. Cross-links:
i2_RESULTS.md §1.2/V15.5 (the propagated window); tier5-family C5c/C5d/C7b (the
convention limits, reused); the (4.8′)/(4.9′) ε-scan re-derivation remains the other
open I2 obligation — if it moves the ceiling again, g_min moves with it but the binding
edges of this battery still do not.*
