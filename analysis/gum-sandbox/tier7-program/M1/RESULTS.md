# M1 — WS4-M4c toy: w(z) from the Γ(H) relaxation profile [F-T7-M1]

**Phase:** Tier 7 program extensions (ROADMAP_v9_PROGRAM.md item M1;
Sec. XI open #4; F-R12 anchor).  **Date:** 2026-07-18.
**Goal.** First independent instantiation of VI.D's relaxation dynamics:
the slow mode δq displaced from its minimum by Hubble drag,
ρ_Λ = ½χ⁻¹δq², w(z) from the dynamics — testing (i) the exact tracker
identity w(z) = −1 + Ω_m(z), (ii) the printed memory-term 𝔴 pairs
(β = −0.3 → 1.28→1.03; β = −1.0 → 1.72→1.05), and (iii) the CPL
projection w_a = +3𝔴Ω_m0(1−Ω_m0) with its always-positive sign — the
F-R12 strengthening verified on dynamics rather than algebra.
**Code:** `m1_wz.py` (algebraic + ODE legs, all gates; ~11 s,
deterministic), `m1_fig.py`.  **Outputs:** `m1_results.json`,
`m1_wz_fig.png`.  **Units:** H₀ = 1; test-field background
E²(x) = Ω_m0 x³ + Ω_Λ0, x = 1+z, Ω_m0 = 0.3.

**Epistemic frame (binding).** Within-model; nothing here bears on
nature.  The whole sector additionally inherits the P1 conditionality
(grand-potential sourcing premise NOT-AUDITABLE-FROM-TEXT; T4 §G4) —
this workstream tests the printed *shape family*, not its premise.

---

## Derivations (printed per gate spec)

**G2 — the exact tracker identity (two lines).**
1. Flat matter+Λ: H²(x) = H₀²[Ω_m0 x³ + (1−Ω_m0)] ⇒
   d ln H²/d ln x = 3Ω_m0 x³ H₀²/H² = **3 Ω_m(z)**.
2. ρ_Λ ∝ H² and continuity (w = −1 + ⅓ d ln ρ/d ln x) ⇒
   **w(z) = −1 + Ω_m(z), exactly, at every z.** ∎

**Memory form, closed.**  Ḣ = −(3/2)Ω_m H² in the same background, so
ρ = αH²M\*² + βḢM\*² = H²M\*²(α + bΩ_m) with **b ≡ −(3/2)β**; with
dΩ_m/d ln x = 3Ω_m(1−Ω_m):
**𝔴(z) = (w+1)/Ω_m = 1 + b(1−Ω_m)/(α + bΩ_m)** — for β < 0: 𝔴 ≥ 1
everywhere, monotone in Ω_m, → 1 in the matter era, maximal at z = 0.

**ODE leg.**  δq̇ = −Γδq + S·H (Γ const, the fast-relaxation/tracking
regime the task fixes).  Linear ⇒ any two solutions differ by Ce^{−Γt}
(the attractor approach rate is Γ *exactly* — measured below).
Gradient expansion of the particular solution:
δq = (S/Γ)[H − Ḣ/Γ + O(Γ⁻²)] ⇒ ρ = (S²/2χΓ²)[H² − 2HḢ/Γ + …]:
the leading term **is** the corpus's αH² with α = S²/2χΓ² > 0, and the
first correction is a Ḣ term with coefficient −2α(H\*/Γ) < 0 — **the
corpus's β < 0 is forced by the dynamics**, not assumed.  To first
order in H/Γ: **𝔴_dyn(z) = 1 + 3(H/Γ)(1 − Ω_m/2)** (verified below).

**Γ(H)-profile family (the workstream title).**  Γ = γ₀H₀(H/H₀)ⁿ ⇒
quasi-static δq ∝ H^{1−n} ⇒ ρ ∝ H^{2(1−n)} ⇒ **𝔴 = 1 − n**: the
relaxation profile sets 𝔴 directly.  n = 0 → pure tracker (𝔴 = 1);
n = 1 (the ROADMAP's literal "Γ = γH") → δq = S/γ = const → w ≡ −1
(Λ, 𝔴 = 0): the tracker class requires Γ ≁ H, and the corpus's αH²
form corresponds to n = 0.  Verified dynamically at n = ½ (below).

**CPL tangent.**  Ω_m(a) = Ω_m0/(Ω_m0+(1−Ω_m0)a³) ⇒
−dΩ_m/da|₁ = 3Ω_m0(1−Ω_m0) ⇒ w_a|_{a=1 tangent} = 3𝔴Ω_m0(1−Ω_m0):
the corpus's projection is exactly the derivative at a = 1. ∎

---

## Gates

| id | spec | measured | verdict |
|---|---|---|---|
| G1 | numerics to 10⁻⁸; ρ_Λ > 0 throughout (q-variable sign theorem's content) | ODE cross-integrator (DOP853 1e−12 vs Radau 1e−10) rel. diff **4.9×10⁻¹¹**; FD-vs-analytic w **1.8×10⁻¹²**; min ρ over all dynamics runs 7.8×10⁻⁷ > 0; memory shape factor min(1+bΩ_m) = 1.045 > 0; ρ = ½δq² ≥ 0 by construction | **PASS** |
| G2 | pure tracker: w(z) = −1+Ω_m(z) to ≤ 10⁻⁶, z ∈ [0,10]; derivation printed | max\|w−(−1+Ω_m)\| = **1.8×10⁻¹²** (8001-pt grid, 4th-order FD); two-line derivation above | **PASS** |
| G3 | CPL fit around a = 0.7: w_a = +3𝔴Ω_m0(1−Ω_m0) within 5% at 𝔴 = 1; sign STRICTLY positive across the family | tangent at a = 1 (the formula's own definition): **w_a = 0.630000** (dev 6×10⁻¹⁰ %) ✓; literal LS windows: [0.6,0.8] → **+1.051** (+67%), [0.5,1] → +0.969, [0.7,1] → +0.853; density-level fits +0.83…+0.96 — no finite window centered below a = 1 meets 5%. Sign: min w_a over 15-member family (𝔴 = 0.25…1.5; β = −0.1…−2; n = ±½; ODE Γ = 50/200/800) × two conventions = **+0.1575 > 0**, strictly | **CONDITIONAL** (5% clause holds only under the a = 1-tangent convention — the corpus's own algebraic definition; sign clause PASS unconditionally) |
| G4 | memory 𝔴 vs printed pairs (1.28→1.03; 1.72→1.05) within 15%; else FAIL + convention diagnosis | task-spec reading (z = 0.5→0): (1.145, 1.278) and (1.325, 1.724) — devs 10.5/**24**% and 23/**64**% → fails.  **Recovered convention (α = 1, β = printed, endpoints z = 0 → z = 2):** (1.2775, 1.0253) and (1.7241, 1.0501) — all four printed values are the 2-dp roundings, devs **0.19/0.46% and 0.24/0.01%**; scan confirms z = 2.00 is where W hits 1.05 for β = −1 (1.82 for the 1.03 target) | **PASS-MATCH** (under recovered z-window; spec-assumed window FAILs — diagnosis below) |
| G5 | Branch-A operational content: (w+1) ∝ Ω_m(z) exactly in tracker limit; memory deviation shape reported | max\|(w+1)/Ω_m − 1\| = **1.8×10⁻¹²** (tracker); memory: (w+1)−Ω_m = bΩ_m(1−Ω_m)/(1+bΩ_m) ≥ 0, 𝔴(z) ≥ 1 and monotone-decreasing in z (measured true for both β); deviation peaks at z = 0.245 (0.093, β = −0.3) / z = 0.14 (0.225, β = −1) and dies as Ω_m → 1 | **PASS** |

## Key numbers

- **Tracker identity:** max deviation 1.8×10⁻¹² over z ∈ [0,10] — exact
  to machine-FD floor; the F-R12 strengthening ("𝔴 = 1 EXACTLY", not
  "O(1)") is verified on independent code.
- **Attractor (the dynamics underwriting the algebra):** approach rate
  measured/Γ − 1 = **−5.4×10⁻⁸** (Γ = 50, two-solution log-slope fit,
  window Γt ∈ [1,15]); quasi-static δq ∝ H confirmed with first-order
  lag 𝔴_dyn(0) − 1 = 2.55/Γ verified: measured/predicted =
  1.067 / 1.016 / **1.004** at Γ/H₀ = 50 / 200 / 800; lag → 0 with
  log-slope −1.02 in Γ.  Tracking entry is late (Γ = H at z = 19/50/128
  for the three Γ): at Γ = 50 the z ≳ 3 range still shows the physical
  entry transient (w spikes to +0.97 at z ≈ 5; δq/δq_qs = 3.5 at z = 10),
  fully relaxed by z ≲ 0.5.
- **Sign structure from dynamics (F-R12's "always-positive w_a" on
  dynamical footing):** every Γ = const run gives 𝔴_dyn(z) ≥ 1 across
  z ∈ [0,10] (min 1.0032 at Γ = 800) — the relaxation ODE *cannot*
  produce 𝔴 < 1, hence w+1 ≥ Ω_m > 0 and w_a > 0, with β_eff < 0 forced
  (derivation above).  Γ(H)-profile: 𝔴 = 1−n; measured 𝔴 = 0.5029 vs
  0.5 at n = ½ (γ₀ = 400).
- **Memory pairs:** 𝔴(0) = 1.2775 (β = −0.3), 1.7241 (β = −1.0);
  𝔴(2) = 1.0253, 1.0501 — the corpus's four printed values to 2 dp.
- **CPL:** tangent w_a = 0.6300 (𝔴 = 1); memory tangents +0.709
  (β = −0.3), +0.749 (β = −1.0); broad-window fits up to +1.05; **all
  positive, every member, every convention.**
- **Self-consistency (one iteration, ρ_Λ(z) fed back into H):** pure
  tracker w(0): −0.700 → −0.490 (Δ = +0.21); β = −0.3: Δw(0) = +0.27;
  β = −1.0: Δw(0) = +0.31.  The pure αH² form has **no self-consistent
  accelerating fixed point** (iteration flows toward EdS) — the corpus's
  own printed pathology ("a pure αH² term tracks and cannot
  accelerate"), now measured; the exact identity is a *test-field*
  statement on the ΛCDM background.

## Method

Test-field: background fixed to ΛCDM (Ω_m0 = 0.3); ρ_Λ read off each
form; w from the continuity equation, w = −1 + ⅓ d ln ρ/d ln(1+z),
via 4th-order central differences on an 8001-point ln(1+z) grid
(cross-checked against closed forms).  ODE leg: dδq/dN = S − (Γ/E)δq
integrated z = 3000 → 0 from δq = 0 (DOP853, rtol 1e−12; Radau
cross-check), w from the exact RHS (no FD noise); rate from the decay
of the difference of two solutions.  CPL: 5-point one-sided derivative
at a = 1 plus LS fits over four windows (plus density-level fits).
Self-consistency: one Friedmann iteration on the grid, ρ_de
renormalized to Ω_Λ0 at z = 0, w re-derived.

## Caveats

1. **G4's convention is our inference.**  The corpus prints
   "𝔴 = 1.28→1.03 / 1.72→1.05" with no z-window; the recovery (all four
   values = 2-dp roundings of 𝔴(z = 0) and 𝔴(z = 2) at α = 1,
   β = printed) is strong (two independent β's, four numbers, ≤ 0.5%)
   but the corpus's β-convention remains not fully specified — the
   roadmap's own honesty flag stands.  Under the task-sheet's assumed
   window (z = 0.5 → 0) the values disagree by up to 64%.
2. **The corpus's "+0.36 at β = −0.3" CPL print is NOT recovered** under
   any convention tried (tangent +0.709; windows +0.78…+0.92): sign
   agrees, magnitude does not.  Flagged as an unresolved corpus-internal
   arithmetic inconsistency with its own 𝔴 > 1 memory values
   (report-only; not gated).
3. **G3's 5% clause is convention-bound:** the printed formula is exactly
   the a = 1 tangent; any finite fit window reaching a ≈ 0.7 measures
   the curvature of Ω_m(a) and overshoots by 35–67%.  The *sign* claim —
   the physics content of F-R12 — is convention-robust.
4. **Test-field is not self-consistent** (Δw(0) = +0.21…+0.31 after one
   feedback iteration; pure-tracker fixed point is EdS).  The corpus
   prints this pathology openly and hangs acceleration on the memory
   term; within this toy the memory term does not cure it (its ρ is
   still ∝ H² up to the Ω_m-modulation).  All gate statements are
   background-ΛCDM statements.
5. **Γ-profile scope:** the ROADMAP's literal "Γ = γH" gives δq = const
   (w ≡ −1), not a tracker; the task-sheet's Γ = const is what
   reproduces the corpus's αH² residual (𝔴 = 1−n, n = 0).  Reported,
   not hidden.
6. Sector-level: P1 (grand-potential sourcing) conditionality inherited
   from T4; this workstream neither tests nor discharges it.

**Promotion note for the coordinator.**  WS4-M4c now has a first
independent toy instantiation: the exact tracker identity holds at
machine precision on independent dynamics; the attractor, the sign
β < 0, and 𝔴 ≥ 1 (hence w_a > 0, freezing class) all *emerge* from the
relaxation ODE rather than being imposed — the F-R12 strengthening is
underwritten at toy grade.  The printed memory pairs are recovered
exactly under a diagnosed z-window convention (z = 0 → 2).  Two honest
debits: the CPL "+0.36" print does not reproduce, and the test-field
reading is load-bearing.  File as **F-T7-M1**.

## Files

- `m1_wz.py` — all measurements (reruns in ~11 s, deterministic).
- `m1_results.json` — every measured number (gates G1–G5, dynamics
  block, self-consistency block).
- `m1_fig.py` / `m1_wz_fig.png` — w(z) and 𝔴(z) with the corpus prints
  marked at the recovered (z = 0, z = 2) endpoints.
