# M3 — F-R9 route (b) priced: on-tube saturation nonlinearity (F-T7-M3)

**Goal.** Quantitatively price the one undischarged rescue route for F-R9 (the 0νββ-insert
second-moment defect, h21/T-H2): does a saturation nonlinearity of the on-tube insertion
amplitude exist that (i) tames R = ⟨A²⟩/⟨A⟩² from 10¹⁴–10¹⁸ to O(1), while (ii) preserving
the mean-rate phenomenology the corpus claims (h21's verified ⟨A⟩/m̄ = 1.000 ± 0.001,
rate ∝ m̄²), and (iii) staying inside the corpus's own scale window?

**Campaign:** Tier 7 program extensions, workstream M3 (ROADMAP_v9_PROGRAM.md).
**Date:** 2026-07-18. **Code:** `m3_saturation.py` (seed 20260718, runtime 38 s),
`m3_fig.py`. **Raw output:** `m3_results.json`. **Figure:** `m3_fig.png`.

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**Route (b) is NOT-VIABLE.** On h21's line-supported web, saturation A → A_sat·tanh(A/A_sat)
never brings R below 2.2×10¹³ anywhere in a 20-decade scan of A_sat/⟨A⟩ (target: R ≤ 10);
R falls only logarithmically in A_sat (measured d ln R/d ln A_sat = 0.04–0.08 in the clip
regime) while the claimed mean rate collapses quadratically (d ln D/d ln A_sat = 1.92–1.96).
Reaching R = 10 would require A_sat ≈ 10^(−4×10⁷) – 10^(−3.4×10⁸) of the on-tube amplitude
— tens to hundreds of millions of orders of magnitude below the corpus's own on-tube
amplitude class — and the mean rate there would be suppressed by 10^(8×10⁷)-class factors.
The 0νββ insert's repair menu gets a closed door, not an entry.

## 1. Gates

| id | spec | measured | verdict |
|---|---|---|---|
| G1 | unsaturated control reproduces h21: ≥ 1 anchor R within ×3; mean 1.000 ± 0.001-class | 4 anchors reproduced at ratio-to-memo 0.984–1.014 (×3 satisfied with margin ~200×); direct-MC re-anchor R = (1.117 ± 0.026)×10¹⁴ vs h21's (1.115 ± 0.003)×10¹⁴ same cell; mean check ⟨A⟩/m̄ = 0.99988–1.00006 on all anchors | PASS |
| G2 | pricing curves R(A_sat), D(A_sat) measured over ≥ 4 decades; report A_sat* (R ≤ 10) and D there | measured over 20 decades (A_sat/⟨A⟩ ∈ [10⁻², 10¹⁸], 81 points, 4 anchors spanning the h21 window R_unsat = 1.1×10¹⁴–8.9×10¹⁷); **A_sat\* does not exist in the scan** — min R = 2.2×10¹³ (audit-ref anchor, scan edge); analytic extension: R = 10 requires log₁₀(A_sat/m̂) = −4.1×10⁷ to −3.4×10⁸, where log₁₀ D_claim = −8.2×10⁷ to −6.9×10⁸ | PASS (curves measured; the physics answer is "no crossing") |
| G3 | pre-registered rule: VIABLE iff \|D−1\| ≤ 0.5 at A_sat* AND A_sat* inside h21's on-tube amplitude range; CONDITIONAL if range not extractable | A_sat* does not exist for any A_sat > 0 reachable in (or analytically beyond) the scan; on-tube amplitude range IS extractable from h21's parametrization (A_on = g(0)·m̂ = 0.91–1.00 m̂, i.e. A_on/⟨A⟩ = 1.2×10¹⁴–1.3×10¹⁸), so no CONDITIONAL grade is needed; required A_sat is ~10⁷·⁶-orders below that range | **NOT-VIABLE** |
| G4 | lemma-level account + the trade quantified | printed in §5: ⟨A_s²⟩ ≤ A_sat⟨A_s⟩ ⇒ R ≤ A_sat/⟨A_s⟩ (taming only at scale A_sat), but support bound R ≥ 1/P(A > 0-class): saturation changes values, not support; measured trade: 1.92–1.96 decades of mean rate lost per 0.04–0.08 decades of R gained (trade ratio 25–46) | PASS |

## 2. Key numbers

**G1 — h21 reproduction (unsaturated control; parametrization reused verbatim):**

| anchor (a, r0, L) | R this run | R h21 memo | ratio | ⟨A⟩/m̄ |
|---|---|---|---|---|
| 52.6 fm, 2 fm, 1.0 μm (audit ref) | 1.116×10¹⁴ | 1.1×10¹⁴ | 1.014 | 1.000038 |
| 52.6 fm, 2 fm, 4.2 μm (corpus central L) | 1.968×10¹⁵ | 2.0×10¹⁵ | 0.984 | 1.000038 |
| 6.8 fm, 2 fm, 4.2 μm (corpus central) | 9.375×10¹⁶ | 9.4×10¹⁶ | 0.997 | 1.000059 |
| 2.0 fm, 1 fm, 4.2 μm (thin extreme) | 8.90×10¹⁷ | 8.9×10¹⁷ | 1.000 | 0.999880 |

Direct importance-sampled MC re-anchor (a = 52.6 fm, r0 = 2 fm, L = 1 μm):
R = (1.117 ± 0.026)×10¹⁴ vs h21 leg-3b (1.115 ± 0.003)×10¹⁴ — 0.2% agreement.

**G2 — the pricing scan** (A_sat/⟨A⟩ ∈ [10⁻², 10¹⁸], per anchor):

| anchor | R_unsat | min R in scan (at x = 10⁻²) | x where D_claim = 1 | R there | ⟨A_s⟩/m̄ there | analytic log₁₀(A_sat/m̂) for R = 10 | log₁₀ D_claim there |
|---|---|---|---|---|---|---|---|
| a=52.6, L=1.0 | 1.12×10¹⁴ | 2.21×10¹³ | 7.1×10⁶ | 4.8×10¹³ | 1.4×10⁻⁷ | −4.1×10⁷ | −8.2×10⁷ |
| a=52.6, L=4.2 | 1.97×10¹⁵ | 3.56×10¹⁴ | 2.9×10⁷ | 8.0×10¹⁴ | 3.6×10⁻⁸ | −1.7×10⁸ | −3.4×10⁸ |
| a=6.8, L=4.2 | 9.37×10¹⁶ | 7.6×10¹⁴ | 5.8×10⁷ | 3.0×10¹⁵ | 1.8×10⁻⁸ | −1.7×10⁸ | −3.4×10⁸ |
| a=2.0, L=4.2 | 8.90×10¹⁷ | 2.96×10¹⁵ | 1.2×10⁸ | 1.2×10¹⁶ | 9.1×10⁻⁹ | −3.4×10⁸ | −6.9×10⁸ |

- Clip-regime log-slopes (x ∈ [10², 10¹⁰]): d ln D_claim/d ln A_sat = 1.92–1.96;
  d ln R/d ln A_sat = 0.042–0.076. **Trade ratio: 25–46 decades of mean rate per decade of R.**
- The other direction of the vise: at the D_true = 0.5 boundary (halving the *unsaturated*
  true rate), A_sat sits at 0.85–0.98 of the on-tube amplitude and R has fallen only 2–9%
  from its unsaturated value.
- Mean-amplitude check under saturation: at the only point where the claimed rate is
  preserved (D_claim = 1), ⟨A_s⟩/m̄ ≈ 10⁻⁸–10⁻⁷ — h21's verified 1.000 ± 0.001 first-moment
  identity is destroyed by 7–8 orders. The corpus cannot keep both its Step-2 mean statement
  and a saturation-preserved rate simultaneously at any A_sat.

**D definitions (stated precisely, per roadmap G2 and h21's mean claim):**
- D_true(A_sat) = ⟨A_s²⟩/⟨A²⟩ — roadmap's ⟨rate⟩_sat/⟨rate⟩, denominator = the unsaturated
  model's true mean rate.
- D_claim(A_sat) = ⟨A_s²⟩/m̄² — vs the corpus's *claimed* mean rate: NR-K3a asserts the
  amplitude is governed by the volume mean m̄ (h21 verified ⟨A⟩/m̄ = 1.000 ± 0.001), so the
  claimed rate is m̄². The G3 verdict is identical under either definition (both D's are
  ≪ 0.5 at any A_sat approaching the R = 10 requirement; no A_sat* exists regardless).

## 3. Method

- **Parametrization: reused verbatim from `theory-audit/h21_mc.py`** (gate condition):
  isotropic Poisson line process, ρ_L = 1/L_web² (φ = πa²/L²), tube radius a, uniform
  amplitude m̂ = 1 on tubes, m̄ = φ; kernel w(r) = exp(−r/r0)/(4π r0 r²); one-tube response
  g(d) by the identical quadrature code; the corpus scale window L ∈ {1…10} μm,
  a ∈ {2…200} fm. Four anchors span h21's R = 1.1×10¹⁴ → 8.9×10¹⁷.
- **Saturation:** A → A_s = A_sat·tanh(A/A_sat) applied to the point amplitude (roadmap
  spec). Saturated moments by generalized Campbell in the sparse limit:
  ⟨S(A)⟩ = ρ_L∫2πd·S(g(d))dd, ⟨S(A)²⟩ = ρ_L∫2πd·S(g(d))²dd + ⟨S(A)⟩² (exact for linear S;
  multi-tube corrections O(λ) relative, λ ~ πρ_L d_max² ~ 10⁻¹³ at physical scales).
- **Validation of the nonlinear formula** where multi-tube effects are *largest*
  (φ = 1.5×10⁻³–3×10⁻³, explicit Poisson and jittered-lattice networks, saturation applied
  to the full multi-tube sum): MC/prediction ratios 0.955–0.983 across all A_sat, identical
  to the offset h21's own *linear* MC shows at those φ (0.967/0.977) — the
  saturation-specific deviation is ≤ 1.5%, and it scales down with φ. At physical scale
  (importance-sampled direct MC, 2×10⁶ environments): agreement 0.03–0.3% at four A_sat
  values spanning 12 decades.
- **Scan:** 81-point log grid, A_sat/⟨A⟩ ∈ [10⁻², 10¹⁸] (20 decades; spec required ≥ 4);
  crossings by log-log interpolation; R = 10 requirement beyond the scan by
  exponential-tail extrapolation of g(d) (decay scale r0, d_req = L/√(9π)).

## 4. Why the door closes (G4 lemma, printed)

**Lemma (what saturation does to a quadratic Campbell functional).** S(A) ≤ A_sat pointwise
gives ⟨A_s²⟩ ≤ A_sat·⟨A_s⟩, hence R = ⟨A_s²⟩/⟨A_s⟩² ≤ A_sat/⟨A_s⟩: the second moment is
bounded and R is tamed *at the scale A_sat* — this is the intuition that made route (b)
worth pricing. **But** a support bound runs the other way: S(0) = 0 and S is monotone, so
saturation changes the *values* of A, never its *support*. By Cauchy–Schwarz,
R ≥ 1/P(A ≳ A_sat-class), and on a sparse line-supported web that probability is
πρ_L·d_eff(A_sat)² with d_eff = a + r0·ln(m̂g(0)/A_sat) — it grows only *logarithmically*
as A_sat shrinks, because off-tube amplitude dies exponentially on the kernel scale r0 ≪ L.
So R(A_sat) ≈ [πρ_L(a + r0 ln(m̂/A_sat))²]⁻¹: 14 decades of clipping buy a factor ~5–300 in
R. Meanwhile the mean and the mean rate are carried by the same rare on-tube configurations
(the φ-fraction with A ≈ m̂ = m̄/φ), so once A_sat < m̂ both collapse ∝ A_sat and A_sat²
respectively. **The competition, measured:** second moment dominated by rare close-tube
configurations (h21 §5), mean not — so the scan shows d ln D/d ln A_sat ≈ 1.94 against
d ln R/d ln A_sat ≈ 0.05: every decade of R relief costs 25–46 decades of the claimed rate.
The h21 motional-narrowing exclusion (web speeds 10⁸–10⁹ c) is untouched by any of this.

## 5. Caveats

1. **The tanh form is one saturation shape.** The closure argument does not depend on it:
   the support bound (§4) holds for *any* monotone S with S(0) = 0 and any ceiling —
   R ≥ [πρ_L d_eff²]⁻¹ with d_eff growing at most logarithmically for any exponentially
   decaying kernel. A hard clip or power-law compression moves the measured curves by O(1).
2. **Sparse-limit formula for the saturated moments** carries O(λ) ~ 10⁻¹³ multi-tube
   corrections at physical scales; validated at moderate φ where the corrections are
   ~10¹¹ times larger and found ≤ 1.5% beyond the linear case's own finite-box offset.
3. **Kernel choice inherited from h21** (corpus prints no 0νββ kernel): r0 = 1, 2 fm
   variants; the closure worsens for shorter-range kernels (d_eff grows slower) and could
   only be undone by a kernel with range ≳ L_web ~ μm, which h21 already priced as
   contradicting the corpus's own fm-scale virtuality.
4. **Saturation applied to the point amplitude A** per the roadmap spec. Applying it
   instead to the tube density m(y) before convolution is the same operation up to the
   kernel smoothing (on-tube values m̂ ↦ A_sat-class) and cannot evade the support bound
   for the same reason.
5. What this run does **not** adjudicate: other rescue routes not named in F-R9 (e.g. a
   smooth bulk component of m — already priced by h21 as costing the corpus its Step-1
   line-support premise), and the promotion question for the four-cell table, which stays
   with the coordinator.

## 6. Deliverables

- `m3_saturation.py` — simulation (deterministic, seed 20260718; 38 s single-threaded).
- `m3_results.json` — G1 anchors, both validation legs, full 81-point scans for 4 anchors,
  crossings, analytic requirements, verdict block.
- `m3_fig.py` / `m3_fig.png` — panel A: R(A_sat/⟨A⟩) for the 4 anchors with the R ≤ 10
  target band (never reached); panel B: D_claim and D_true for the corpus-central anchor
  with the |D−1| ≤ 0.5 band and the on-tube amplitude marked.
- This memo. **F-R9 route (b) status: priced and closed — NOT-VIABLE.**

— end of memo —
