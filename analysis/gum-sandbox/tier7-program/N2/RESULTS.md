# N2 — F-R9 route (a): the honest second-moment treatment (F-T7-N2)

**Goal.** F-R9's route (a): instead of rescuing the corpus's ⟨A⟩²-based mean-rate claim
(NR-K3a Step 2), compute what an HONEST second-moment treatment gives at the observable
level — spatial self-averaging of the sample-total rate, the rate's Lorenz/burstiness
structure — and adjudicate the consequence for the corpus's 0νββ-insert (K3a)
phenomenology (P-ν4, the g-window, the J3 battery re-run) under a pre-registered
one-of-three decision frame, no tuning.

**Campaign:** Tier 7 Phase N (ROADMAP_v9_PROGRAM.md addendum, N2).
**Date:** 2026-07-18. **Code:** `n2_secondmoment.py` (seed 20260718, runtime 9 s),
`n2_fig.py`. **Raw output:** `n2_results.json`. **Figure:** `n2_fig.png`.
**Prior art (mandatory reads, honored):** theory-audit/h21_RESULTS.md + h21_mc.py
(parametrization reused verbatim), M3/RESULTS.md (route (b): NOT-VIABLE),
theory-audit/j3_RESULTS.md (Majoron battery re-run at g ≈ 7.8×10⁻¹⁰).

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**The honest sample-total rate is R × the corpus's claimed rate (R = 1.1×10¹⁴–8.9×10¹⁷
at the four anchors), self-averaging to relative width ≤ 5×10⁻⁴ for any macroscopic
sample, carried almost entirely by the ~10⁻¹⁸–10⁻¹³ hot-site fraction (on-tube nuclei
alone: 94.3–99.6%); propagated at rate ∝ g², the printed g-window maps to
g_equiv = √R·g = 8.5×10⁻³–12.3, i.e. 6.9–10.1 orders above the corpus's own P-ν4
far-horizon line — but the archive prints NO experimental Majoron-mode rate bound to
compare against, so the pre-registered verdict is (iii) INDETERMINATE-FROM-ARCHIVE, with
the blocking unprinted number named (g_lim, the Majoron-mode 0νββ detectability/limit
line) and the conditional kill filed: any archive-external g_lim < 8×10⁻³ excludes the
honest rate at every anchor and every g in the printed window. Route (i)
(survive-by-repricing) is closed by arithmetic: the required retreat (f_req = √R·f =
4.9×10²–4.4×10⁴ TeV, ε_req = 1.8×10¹¹–1.4×10¹⁵) overshoots the corpus's own printed
ε-ceiling 2.604×10⁻³ by 14–18 orders. The J3 battery itself is untouched (its rows bind
on the actual g). S1 is untouched either way.**

## 1. Gates

| id | spec | measured | verdict |
|---|---|---|---|
| G1 | h21 anchor reproduction: R within ×3 at ≥ 2 anchors; mean 1.000-class | 4/4 anchors at ratio-to-memo 0.984–1.014 (×3 met with ~200× margin); ⟨A⟩/m̄ = 0.99988–1.000059; direct importance-sampled MC re-anchors: R = (1.1170 ± 0.0026)×10¹⁴ vs h21 leg-3b (1.115 ± 0.003)×10¹⁴ (audit ref) and R = (8.879 ± 0.019)×10¹⁷ vs Campbell 8.90×10¹⁷ (thin extreme) | PASS |
| G2 | burstiness measured: Lorenz fractions (top 10⁻ᵏ of sites, k = 2..12) with MC errors at ≥ 2 anchors | measured at all 4 anchors, quadrature + 2×10⁶-environment importance-sampled MC with 20-slice jackknife: share(top 10⁻ᵏ) = 1.000000 for every k = 2..12 at every anchor (the whole rate sits deeper: kernel-reach site fraction 1.9×10⁻¹⁵–2.0×10⁻¹³); knee resolved at k = 14–16 (audit ref: 0.999987 at k=14, 0.1186 at k=15); on-tube share 0.9426–0.9960 (MC agrees to ≤ 6×10⁻⁴, SE ≤ 2×10⁻⁴); q₅₀ = 2.96×10⁻¹⁹–4.21×10⁻¹⁵ | PASS |
| G3 | propagation table (printed rate → honest rate → battery consequence) + ONE pre-registered verdict with arithmetic | table in §4 (per anchor: R → √R → g_equiv at the window edges, the J3 edge, the P-ν4 line → consequence); verdict: **(iii) INDETERMINATE-FROM-ARCHIVE** — blocking unprinted number: g_lim (equivalently the mode's rate normalization T₁/₂(g)); route (i) closed by the ε-ceiling overshoot 6.8×10¹³–5.4×10¹⁷; conditional kill filed at F-T7-N2 (fires for any external g_lim < 8×10⁻³) | PASS (verdict delivered) |
| G4 | honest-limits paragraph: upstream K3a physics; S1 untouched; what remains corpus-side | printed in §6 | PASS |

## 2. Key numbers

**G1 — anchor reproduction (h21 parametrization verbatim):**

| anchor (a, r0, L) | R this run (Campbell) | R h21 memo | ratio | ⟨A⟩/m̄ |
|---|---|---|---|---|
| 52.6 fm, 2 fm, 1.0 μm (audit ref) | 1.116×10¹⁴ | 1.1×10¹⁴ | 1.014 | 1.000038 |
| 52.6 fm, 2 fm, 4.2 μm (corpus central L) | 1.968×10¹⁵ | 2.0×10¹⁵ | 0.984 | 1.000038 |
| 6.8 fm, 2 fm, 4.2 μm (corpus central) | 9.375×10¹⁶ | 9.4×10¹⁶ | 0.997 | 1.000059 |
| 2.0 fm, 1 fm, 4.2 μm (thin extreme) | 8.90×10¹⁷ | 8.9×10¹⁷ | 1.000 | 0.999880 |

**Spatial self-averaging (physics point 1).** The web is static on experiment
timescales — h21 §6 priced motional narrowing at web transport 1.5×10⁸–1.2×10⁹ c
(campaign ledger prints 5×10⁷–1.2×10⁹ c); each decay event samples a frozen snapshot,
so the ensemble average over events IS the spatial average over nucleus sites. For a
sample of N = 10²⁷ nuclei the total rate T = Σᵢ A(xᵢ)² has
RSD²(T) = (K₄−1)/N + (9/8)/(ρ_L π R_b²) (point-sampling + web-realization terms,
Campbell-quadrature): K₄ = ⟨A⁴⟩/⟨A²⟩² = 1.15×10¹⁴–1.36×10¹⁸, and the worst total RSD
across sample diameters 1 cm–1 km is **5.0×10⁻⁴** (dominated by the tube-length
fluctuation term at 1 cm; ≤ 3.7×10⁻⁵ at ≥ 1 m). **The honest sample-total rate is
N·⟨A²⟩ = R × the corpus's claimed N·m̄² rate, to sub-10⁻³ relative width — R is not a
tail curiosity; it IS the observable.**

**G2 — the hot-site picture (burstiness):**

| anchor | φ (on-tube fraction) | kernel-reach fraction (A > 0) | on-tube rate share (quad / MC) | q₅₀ | q₉₉ | per-hot-site rate ÷ claimed mean rate | hot nuclei in N = 10²⁷ |
|---|---|---|---|---|---|---|---|
| 52.6 fm, 1.0 μm | 8.7×10⁻¹⁵ | 2.0×10⁻¹³ | 0.9960 / 0.9962 ± 6×10⁻⁵ | 4.2×10⁻¹⁵ | 8.6×10⁻¹⁵ | 1.3×10²⁸ | 8.7×10¹² |
| 52.6 fm, 4.2 μm | 4.9×10⁻¹⁶ | 1.1×10⁻¹⁴ | 0.9960 / 0.9962 ± 8×10⁻⁵ | 2.4×10⁻¹⁶ | 4.9×10⁻¹⁶ | 4.1×10³⁰ | 4.9×10¹¹ |
| 6.8 fm, 4.2 μm | 8.2×10⁻¹⁸ | 7.6×10⁻¹⁵ | 0.9674 / 0.9675 ± 2×10⁻⁴ | 3.5×10⁻¹⁸ | 9.8×10⁻¹⁸ | 1.5×10³⁴ | 8.2×10⁹ |
| 2.0 fm, 4.2 μm | 7.1×10⁻¹⁹ | 1.9×10⁻¹⁵ | 0.9426 / 0.9432 ± 2×10⁻⁴ | 3.0×10⁻¹⁹ | 1.1×10⁻¹⁸ | 2.0×10³⁶ | 7.1×10⁸ |

Share(top 10⁻ᵏ), k = 2..12: **1.000000 at every anchor** (in-model exact for
q > kernel-reach fraction; MC confirms, jackknife SE ≤ 10⁻¹⁵ where nontrivial). The
knee sits at k ≈ 14–16 (audit ref: 0.999987 at k = 14, 0.1186 at k = 15, 0.0119 at
k = 16 — inside the tube the rate is near-uniform, so share ≈ q/φ there). The honest
rate is carried by near-tube shells: ~10⁹–10¹³ nuclei per 10²⁷ decay at 10²⁸–10³⁶ ×
the corpus's claimed per-nucleus rate; the other (1 − 10⁻¹³)-fraction is quiescent.

## 3. Method

- **Parametrization reused verbatim from `theory-audit/h21_mc.py`** (gate condition):
  isotropic Poisson line process, ρ_L = 1/L², φ = πa²/L², m̂ = 1 on tubes, m̄ = φ;
  kernel w(r) = exp(−r/r0)/(4π r0 r²); identical quadrature code for the one-tube
  response g(d); h21's leg-3b importance-sampled physical-scale MC reused for the
  direct re-anchors and extended (same sampling law) to the Lorenz estimator.
- **Lorenz structure:** sparse limit (encounter probability λ ~ πρ_L d_max² ~ 10⁻¹³ —
  h21's own regime): a site's amplitude is the nearest-tube response g(d), d with
  measure 2πρ_L d dd; g verified monotone nonincreasing on the table (worst relative
  increase ≤ 1.4×10⁻⁵), so the top-q site set is a distance ball: share(q) =
  ∫₀^{d_q} 2πt g² dt / J₂, q = πρ_L d_q². Quadrature on the g-table grid + MC with
  20-slice jackknife (2×10⁶ environments/anchor).
- **Self-averaging:** Campbell moments J₂ = ∫2πd g² dd, J₄ = ∫2πd g⁴ dd →
  K₄ = (ρ_L J₄ + 3(ρ_L J₂)²)/(ρ_L J₂)²; web-realization term from the total tube
  length in a ball (Poisson line process: Var/mean² = (9/8)/(ρ_L πR_b²)).
- **Propagation:** Majoron-mode rate ∝ g² (corpus dictionary g = m_ν/f); honest rate
  = R × printed rate ⇔ g_equiv = √R·g. Comparators, all corpus-printed: window
  g ∈ [8×10⁻¹⁰, 1.3×10⁻⁸] (Ω VII.I), P-ν4 line g ~ 10⁻⁹ (far-horizon), J3 repaired
  low edge 7.8×10⁻¹⁰, f-dictionary f = 1.173√ε GeV (J3 A1), ε-bracket
  [1.5×10⁻⁶, 2.604×10⁻³] (B-ν1′ floor; repaired ceiling).

## 4. G3 — the propagation table and the pre-registered adjudication

R → √R → g_equiv = √R·g at the four corpus-printed g values:

| anchor | R | √R | g_equiv @ 8×10⁻¹⁰ (window low) | @ 7.8×10⁻¹⁰ (J3 edge) | @ 10⁻⁹ (P-ν4) | @ 1.3×10⁻⁸ (window top) | retreat to keep P-ν4 printed: f_req (ε_req) | ε-ceiling overshoot |
|---|---|---|---|---|---|---|---|---|
| 52.6 fm, 1.0 μm | 1.12×10¹⁴ | 1.06×10⁷ | 8.5×10⁻³ | 8.2×10⁻³ | 1.1×10⁻² | 0.137 | 494 TeV (1.8×10¹¹) | 6.8×10¹³ |
| 52.6 fm, 4.2 μm | 1.97×10¹⁵ | 4.44×10⁷ | 3.5×10⁻² | 3.5×10⁻² | 4.4×10⁻² | 0.577 | 2.1×10³ TeV (3.1×10¹²) | 1.2×10¹⁵ |
| 6.8 fm, 4.2 μm | 9.37×10¹⁶ | 3.06×10⁸ | 0.245 | 0.239 | 0.306 | 3.98 | 1.4×10⁴ TeV (1.5×10¹⁴) | 5.7×10¹⁶ |
| 2.0 fm, 4.2 μm | 8.90×10¹⁷ | 9.43×10⁸ | 0.755 | 0.736 | 0.943 | 12.3 | 4.4×10⁴ TeV (1.4×10¹⁵) | 5.4×10¹⁷ |

**Where the printed window lands.** In g-equivalent units the whole printed window
[8×10⁻¹⁰, 1.3×10⁻⁸] maps to [8.5×10⁻³, 12.3] — **6.9–10.1 orders above the corpus's
own P-ν4 far-horizon line (g ~ 10⁻⁹) and 5.8–9.0 orders above the corpus's own window
top**; at 2 of 4 anchors g_equiv reaches or passes O(1) (the honest observable mimics a
nonperturbative coupling). Against the J3 battery: rows R1/R1b/R2/R3/R4 and the lab
perimeter all bind on the ACTUAL coupling g through astro/lab emission; R multiplies
only the 0νββ insertion's spatial second moment, so **"SHIFTS-HARMLESSLY at
7.8×10⁻¹⁰" survives untouched as a battery statement** — what inverts is J3's
signature-row remark: the new strip's far-horizon "retreat" becomes a 14–18-order
rate ADVANCE.

**The pre-registered one-of-three (arithmetic, no tuning):**

- **(i) Survives with re-priced parameters — CLOSED.** Keeping the printed
  far-horizon observable requires g′ = g/√R, i.e. f_req = √R·f = 4.9×10²–4.4×10⁴ TeV
  against the printed f-window [4, 60] MeV (7–9 orders outside), i.e.
  ε_req = (f_req/1.173 GeV)² = 1.8×10¹¹–1.4×10¹⁵ against the printed two-sided
  ε-bracket [1.5×10⁻⁶, 2.604×10⁻³] — ceiling overshoot 6.8×10¹³–5.4×10¹⁷
  (≈ 0.6R, as it must be). No admissible retreat exists inside the corpus's own
  printed parameter space; g = m_ν/f is pinned two-sided by the corpus's own
  cosmological floor (B-ν1′) and ħ-entrainment ceiling.
- **(ii) Overshoots a corpus-printed bound — NOT AVAILABLE AS PRINTED.** J3 §0's
  source census stands: the archive prints NO experimental Majoron-mode 0νββ bound,
  no T₁/₂, no SN band edges, no detectability line — the battery's only numeric rows
  bind on actual g and are untouched. There is no printed number for the honest rate
  to overshoot.
- **(iii) INDETERMINATE-FROM-ARCHIVE — the verdict.** The blocking unprinted number
  (the L3-style isolated residue): **g_lim, the experimental Majoron-mode 0νββ
  limit/detectability line in g units (equivalently the mode's rate normalization
  T₁/₂(g))**. The conditional statement is filed at F-T7-N2 as a kill sharpening:
  **for ANY archive-external g_lim < 8×10⁻³, the honest second moment converts P-ν4
  from far-horizon signature to already-excluded prediction at EVERY anchor and EVERY
  g in the printed window** (min over the table: g_equiv = 8.2×10⁻³ at the audit-ref
  anchor and the J3 edge). The corpus can only keep the far-horizon stamp by printing
  a detectability line above 8×10⁻³ — five-plus orders above its own entire coupling
  window.

## 5. What route (a) establishes (the lemma-level statement)

For a static line-supported web, the sample-total 0νββ-insert rate of a macroscopic
sample is a SPATIALLY SELF-AVERAGED quantity: T = N⟨A²⟩ = N·R·m̄² to relative width
≤ 5×10⁻⁴ (N = 10²⁷, any cm–km sample), with the rate carried by the ~φ–10⁻¹³ hot-site
fraction at per-site enhancements 1/φ² = 10²⁸–10³⁶. There is no "rare-fluctuation"
escape: R is not the tail of a distribution an experiment might miss — it is the mean
of the true quadratic observable, and burstiness only concentrates WHERE the decays
happen, not how many. The corpus's choice is binary: adopt the R-enhanced rate (and
face §4), or retreat to route (c) — a smooth bulk m-component — at the cost h21
already priced: Step 1's own line-support theorem (Theorem K-2).

## 6. Honest limits (G4)

1. **Upstream K3a physics stays corpus-side.** The kernel is unprinted (h21 caveat 3;
   r0 = 1, 2 fm variants change nothing in the thick regime); the insertion's
   linearity in m(y) is the corpus's own model (route (b)'s saturation escape is
   M3-closed: NOT-VIABLE); straight static tubes and uniform m̂ per h21 (geometry
   class shifts R < 3%; non-uniform |m| only raises the second moment).
2. **The g² scaling and the g_equiv unit** are the workstream's pre-registered frame
   (rate ∝ g² for the Majoron mode, per the corpus's own dictionary); R itself is
   g-independent — expressing R as √R in g units changes presentation, not content.
3. **S1 (Σm_ν) is untouched either way** — F-R9 wounds the K3a 0νββ MODE-rate claim
   (P-ν4's grading and the V-K1 four-cell phenomenology), not the mass calibration
   (h21's first-moment verification 1.000 ± 0.001 stands). S6's occurrence stake,
   hash, funnels and kill conditions are verbatim-unchanged (Watch memo).
4. **Route (a) is a recomputation, not a repair.** It does not discharge F-R9; it
   prices what honesty costs. The corpus-side residue: print g_lim (or the mode's
   T₁/₂(g) normalization) and adjudicate §4's conditional kill, or adopt route (c)
   at the line-support cost. The battery (J3) needs no re-run: its rows do not see R.
5. **Sparse-limit Lorenz estimator** carries O(λ) ~ 10⁻¹³ multi-tube corrections
   (h21-validated regime); the k = 2..12 shares equal to 1.000000 are in-model exact
   statements (zero amplitude beyond kernel reach; physical exponential tails are
   e⁻¹⁰⁰-class there).

## 7. Deliverables

- `n2_secondmoment.py` — computation (deterministic, seed 20260718; 9 s).
- `n2_results.json` — G1 anchors + direct-MC re-anchors, self-averaging table,
  full Lorenz tables (quadrature + MC with jackknife errors, k = 2..16), propagation
  table, verdict block.
- `n2_fig.py` / `n2_fig.png` — panel A: rate Lorenz curves at the 4 anchors with the
  on-tube fraction marked; panel B: the printed g-window mapped to g-equivalent units
  per anchor vs the corpus's printed window, P-ν4 line, and g_equiv = 1.
- This memo. **F-R9 route (a) status: executed — verdict (iii)
  INDETERMINATE-FROM-ARCHIVE with the conditional kill sharpening filed as F-T7-N2
  (blocking number: g_lim); route (i) repricing closed by arithmetic; battery
  untouched; S1 untouched.**

— end of memo —
