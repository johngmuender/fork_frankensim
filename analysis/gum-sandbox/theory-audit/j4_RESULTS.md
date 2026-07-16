# J4 — multi-ε saturated re-scan: the (4.8′) exponent measured (analysis layer)

**Phase:** ROADMAP_v6 J4 · 2026-07-16 · theory-audit numerics subagent (Fable-class).
**Charge:** J2 defect D3 — the repaired law 𝔠_sat(ε) = 3.2011·(1 + 1.152ε) is a LINEAR
class fitted on TWO measured saturated points. Densify the small-ε saturated branch,
measure the exponent, and re-propagate (4.9′) per surviving reading.
**Computations:** `j4_epslaw.py` (deterministic, i1 machinery reused unmodified) →
`j4_results.json` (stage-flushed). **WITHIN-MODEL ONLY; the coordinator adjudicates.**

## §0 Provenance: which two points J2 used, and how i1's curve relates

**J2's linear fit used exactly two points, both from `tier2-closure/axi4_locked_results.json`
(ring-saturated engine):** `saturated.ring_t0` 𝔠 = 3.203459 (t → 0; exact 𝔠_sat = 3.2011247,
engine bias +7.3×10⁻⁴ = J2's D4) and `saturated.ring_eps` 𝔠 = 3.385470 (ε = 0.05) →
s_exact = 1.1518, s_engine = 1.1359 (G0b/G0c reproduce both).
**i1's curve is a different family, same objective:** the 12-mode axi3-basis saturated
closure G_t on the ε-dial hedgehog base (i1_r5test.py `saturated_at`), measured at
ε ∈ {0.02, …, 2.0}. At ε = 0.05: i1 gives 3.382843 — 7.8×10⁻⁴ rel BELOW axi4's ring value,
i.e. the tighter of the two upper bounds on the true G_t*; its exact-base slope reading at
0.05 is 1.135 vs axi4's 1.152 (≈1.5% family-amplitude spread, same class as J2's D4).
Consequence for J4: **within-family consistency (i1 machinery reused unmodified at every
new ε) decides the EXPONENT; the axi4-vs-i1 family spread is an AMPLITUDE band** carried
through §3 as separate readings. Both families are variational upper bounds; the exact
ε → 0 anchor 64√2/9π is the only point known exactly.

## §1 The densified small-ε dataset (MEASURED)

New solves: i1 dial + `saturated_at` reused unmodified, full i1 budgets (N_OPT/N_POL/
N_FIN/N_COARSE = 192/320/480/320, B_TRI/B_COLD/B_POL = 220/1500/260 — **no coarsening
needed**; 249 s total, well inside budget; no non-convergent point, no boundary flag).
Two grid resolutions per point (N=480 / N=320) + the dial's own N×2 / rmax+2 ladder.

| ε (achieved) | t | 𝔠_sat (N480) | N320 | conv rel | hh-core | floor 3.2011+12t | dial N2x |
|---|---|---|---|---|---|---|---|
| 0.0050005 | 6.6331e-4 | **3.219560** | 3.219453 | 3.3e-5 | 3.26538 | 3.20908 OK | −5.7e-6 |
| 0.0099997 | 1.4097e-3 | **3.237616** | 3.237497 | 3.6e-5 | 3.28447 | 3.21804 OK | −2.8e-6 |
| 0.0199991 | 3.0101e-3 | **3.273947** | 3.273806 | 4.3e-5 | 3.32299 | 3.23725 OK | −1.4e-6 |
| 0.0300056 | 4.7046e-3 | **3.310313** | 3.310152 | 4.9e-5 | 3.36187 | 3.25758 OK | −9.0e-7 |

Full fit dataset = these four + i1's rows at ε = 0.05…2.0 (achieved-ε from i1's dial,
conv flags respected; all `floor_ok`). **Reproduction gate G1a: 𝔠_sat(0.02) = 3.2739466
vs i1's 3.2739555 — rel −2.7×10⁻⁶, ×16 inside the 4.3×10⁻⁵ convergence band** (G1b: dial
ε agrees to 1×10⁻⁴). Monotone rising (G1d); analytic floor respected everywhere (G1c).

## §2 Exponent fits (window, classes, exclusions)

y = 𝔠_sat(ε)/𝔠₀ − 1, 𝔠₀ = 64√2/9π exact; per-point σ = max(grid conv, 3×10⁻⁵)·𝔠;
CIs residual-rescaled (χ²/dof → 1). Primary window **W1: ε ≤ 0.05 (n = 5)** — leading
order defensible there because the measured quadratic correction is ≲ 2% of the linear
term at 0.05 (and W4 shows curvature onset by ε ~ 0.2). Fit table (per window):

| window | (i) aε | (ii) aε^{2/3} | (iii) aε^p, p [CI] | (iv) aε + bε^{2/3} |
|---|---|---|---|---|
| W1 ε≤0.05 | a=1.1368, χ²=6.8 | a=0.341, χ²=57186 | **p=0.997 [0.996, 0.998]**, a=1.123 | 1.126ε + 0.0034ε^{2/3}, χ²=1.2 |
| W2 ε≤0.03 | a=1.1379, χ²=5.0 | a=0.306, χ²=23124 | p=0.996 [0.993, 0.998] | 1.123ε + 0.0041ε^{2/3} |
| W3 ε≤0.10 | a=1.1356, χ²=10.8 | a=0.381, χ²=152170 | p=0.998 [0.997, 0.998] | 1.127ε + 0.0031ε^{2/3} |
| W4 ε≤0.20 | a=1.1482, χ²=793 | a=0.436, χ²=426476 | p=1.012 [1.007, 1.016] | 1.177ε − 0.0123ε^{2/3} |

- **Best-fit exponent p = 0.997, honest envelope [0.996, 1.012]** (profile CI ⊕ jackknife
  [0.996, 0.998] ⊕ free-intercept fits M3f ⊕ all windows; W4's 1.012 is the curvature of
  the full law leaking in, not small-ε physics). Model-free local exponents:
  0.985±0.010 (0.005–0.01), 0.997±0.005, 0.998±0.006, 0.998±0.004, 0.998±0.002 (0.05–0.1).
- **(ii) pure 2/3 is EXCLUDED as leading order**: χ²(2/3)/χ²(lin) = 8445 on W1; its
  amplitude is not even stable (0.31→0.44 across windows). **(iv) measures the 2/3
  admixture directly: b = 0.0034 ± 0.0014 — consistent with zero at ~2.4σ(fit), 0.3% of
  the linear term at ε = 0.05.** (i) linear survives everywhere; (iii) ⇒ p = 1 within CI.
- **Intercept gate G2c:** free-intercept linear fit gives 𝔠₀_fit = 3.201371 ± 0.000097 vs
  exact 3.201125 (rel +7.7×10⁻⁵, +side) — consistent with the exact anchor given the
  12-mode family's upper-bound bias; ×9 smaller than axi4's t→0 engine bias (+7.3×10⁻⁴).
- Amplitude (i1 family, W1): s = 1.137 (fixed intercept) / 1.134 (free intercept) vs
  J2's axi4 reading 1.152 — the ~1.3% two-family spread of §0, carried below as readings.

## §3 (4.9′) propagation per reading (J2/J3 formulas verbatim)

Ceiling solves law(ε) = Δ_lock = 3×10⁻³ (band [2.4, 3.8]×10⁻³); old ceiling
(3×10⁻³/0.42)^{3/2} = 6.037×10⁻⁴; margins ε_q/ceiling, ε_q = (f_q/1.7)², f_q ∈ [0.14,
0.20] GeV; f_top = 0.69·√ceiling·1700 MeV (f_lo = 3.709 MeV fixed); g_min = 0.0468 eV/f_top.

| reading | ceiling | ×old | margins [lo, c, hi] | worst-joint | f_top MeV | g_min |
|---|---|---|---|---|---|---|
| RL J2 axi4 s=1.152 (control) | 2.605e-3 | 4.31 | [2.60, 3.84, 5.31] | ×2.06 | 59.87 | 7.82e-10 |
| RL fit, i1-family s=1.137 | 2.639e-3 | 4.37 | [2.57, 3.79, 5.24] | ×2.03 | 60.26 | 7.77e-10 |
| RP free-p best (1.123ε^0.997) | 2.616e-3 | 4.33 | [2.59, 3.82, 5.29] | ×2.04 | 60.00 | 7.80e-10 |
| RP envelope p=0.996 | 2.610e-3 | 4.32 | [2.60, 3.83, 5.30] | ×2.05 | 59.92 | 7.81e-10 |
| RP envelope p=1.012 | 2.715e-3 | 4.50 | [2.50, 3.68, 5.10] | ×1.98 | 61.12 | 7.66e-10 |
| RC combo lin+2/3 | 2.608e-3 | 4.32 | [2.60, 3.84, 5.31] | ×2.05 | 59.90 | 7.81e-10 |
| ~~R23 pure 2/3~~ (EXCLUDED, ref) | 8.26e-4 | 1.37 | [8.21, 12.10, 16.75] | ×5.76 | 33.72 | 1.39e-9 |

- **G3a control PASS:** the linear reading reproduces J2/J3 exactly (2.6047×10⁻³, ×4.315,
  [2.60, 5.31] c 3.84, 59.87 MeV, 7.82×10⁻¹⁰) — formula-reuse verified.
- **ROBUST across surviving readings:** ceiling ∈ [2.60, 2.71]×10⁻³ (window top now
  87–90% of 3×10⁻³); thinning ×[4.31, 4.50]; f_top ∈ [59.9, 61.1] MeV; g_min ∈ [7.66,
  7.82]×10⁻¹⁰; margins central 3.7–3.8; breach needs a further ×2.50–2.60 (central).
  **The ×4.3 loosening, the ≈60 MeV f-top and the ≈7.8×10⁻¹⁰ Majoron edge are CONFIRMED
  and promoted from linear-reading-specific (J2 D3) to measured-exponent-backed.**
- **READING-SPECIFIC (now excluded):** the 2/3 ceiling ≈ 6–8×10⁻⁴ class and its ×11-wide
  Q-6′ margin (J2's R6, i2's "2/3 reading widens to ×11") die with the exponent. NOTE:
  R6's ×11.41 used amplitude 0.4243 pinned at the axi4 0.05 point; the honest W1 2/3 fit
  gives 0.341 → ceiling 8.3×10⁻⁴, margin ×8.2 — same verdict, different arithmetic.
- **Dichotomy (G3b): survives at EVERY reading including the excluded reference** — all
  six quarks censored, min worst-joint margin ×1.98 (RP hi-envelope + Δ_lock-high);
  leptons free (floor 1.5×10⁻⁶ < 1×10⁻⁵ < every ceiling). J2's "no reading within a
  factor 2 of breach" thins to ×1.98 at the new worst joint — grazing, not crossed.

## §4 Cross-checks and gates (14/14 PASS, one invocation, restart-fresh)

G0a exact anchor; G0b/G0c J2's two axi4 points + slopes 1.152/1.136 reproduced;
G1a ε=0.02 reproduction vs i1 (−2.7×10⁻⁶ rel); G1b dial agreement; G1c analytic floors;
G1d monotonicity; G1e in-band convergence ≤ 4.9×10⁻⁵; G2a exponent envelope excludes 2/3;
G2b χ² ratio 8445; G2c intercept vs 64√2/9π (+7.7×10⁻⁵, within fit σ + upper-bound bias);
G3a J2/J3 control reproduction; G3b dichotomy at every reading; G3c ×4.3 confirmed.
No disagreement found ⇒ no cross-check defect filed.

## §5 Defects (printed honestly)

- **J4-D1 (two-family amplitude spread, inherits J2-D4):** axi4 ring vs i1 12-mode give
  s = 1.152 vs 1.137 (~1.3%); axi4's t→0 base carries +7.3×10⁻⁴ bias vs the i1 family's
  +7.7×10⁻⁵. The ceiling band [2.60, 2.64]×10⁻³ across the two linear readings is real
  method spread, not statistics. Not material to any verdict.
- **J4-D2 (upper-bound family):** every 𝔠_sat(ε>0) is a 12-mode variational upper bound
  (downward-open ≤ 2.5% per i1's mode ladder). The exponent is trusted because the bias
  is smooth and mostly common-mode (free-intercept fits move p by < 0.01); a family
  enlargement could still shift the amplitude a few ×10⁻³ downward → ceiling slightly up.
- **J4-D3 (smallest-ε point):** the local exponent 0.985±0.010 on [0.005, 0.01] sits 1.5σ
  below 1 — consistent with the intercept-bias leakage of D2, but ε < 0.005 was not
  probed; the claim "p = 1" is measured on [0.005, 0.2], asymptotics below that inferred.
- **J4-D4 (W4 drift):** p(W4 ε≤0.2) = 1.012 — quadratic curvature of the full law, why
  the leading-order window must stop at ε ≈ 0.05–0.1. Envelope includes it anyway.
- **J4-D5 (scope):** Δ_lock = 3×10⁻³, ε_q, f-map 0.69√ε·1.7 GeV, m₃ = 0.0468 eV inherited
  verbatim from J2/J3 — J4 re-measures the ε-law only. λ* shape re-scan remains OPEN.

## §6 Verdict (analysis layer — within-model; the coordinator adjudicates)

**The saturated branch's small-ε law is measured LINEAR: p = 0.997, profile CI [0.996,
0.998], all-readings envelope [0.996, 1.012] — the 2/3 class is excluded (χ² ratio 8445;
direct 2/3-admixture 0.0034 ± 0.0014, ≤0.3% of linear at ε = 0.05), and the ε→0 intercept
reproduces 64√2/9π to +7.7×10⁻⁵.** J2's defect D3 is DISCHARGED: the "+1.15ε" class now
rests on a 7-point measured scan (ε = 0.005–0.2, two resolutions each), not two points;
amplitude s = 1.137 (i1 family) vs 1.152 (axi4) is a 1.3% method spread. Downstream, all
surviving readings give ceiling 2.60–2.71×10⁻³: **the ×4.3 loosening (now ×4.31–4.50),
the ≈60 MeV f-window top (59.9–61.1) and the 7.8×10⁻¹⁰ Majoron low edge (7.66–7.82) are
confirmed as measured, no longer linear-assumption-specific; the 2/3 reversion scenario
(ceiling ≈ 6×10⁻⁴, ×11 margin) is excluded.** The confinement dichotomy survives at every
reading (min worst-joint ×1.98 — a hair inside J2's "factor 2 of breach" phrase, printed
honestly). Files: `j4_epslaw.py` (14/14 gates, deterministic, 249 s), `j4_results.json`
(stage-flushed per solve point), this memo.
