# T3 — THEORY AUDIT: THE REPAIRED CLOSURE (the constructive companion)

**Auditor:** theory-audit subagent (Fable-class), 2026-07-16.
**Charge:** F-R4 proved the SDiff inertia dial frozen (`tier2-closure/axi_FR4_NOTE.md`);
F-R5 proved the corpus's stated variational problem selects halo saturation and that
every printed tuple is a halo-unstable saddle (`axi3_ADJUDICATION.md`, tier4); G3
closed the saturated objective in exact form: G\* = 16√2/9, attained by the oblate
compacton, saturated closure 𝔠 = 64√2/(9π) at κ = 1/√2 exactly
(`gstar_RESULTS.md`). This memo is the **constructive** piece: formalize the
candidate WELL-POSED closures the corpus could adopt, compute what each selects,
propagate every downstream constant, and steelman any reading under which the
printed 𝔠₀ = 128√42/(105π) = 2.514754 is recoverable. Argument-level, within-model
only; nothing here bears on nature.

**Sources:** paper Secs. IV.A–J, IX (S4′), Apps. A.4, B, G, I.4;
`tier2-closure/{axi,axi3,axi4_locked,gstar}_RESULTS.md` + adjudications;
`SESSION_HANDOFF.md` §§2–3 (unit maps, exact thresholds);
`theory-audit/T2_closure_pillars.md` (the surviving kinematic pillars);
`tier5-family/family_ADJUDICATION.md` (P-O1 blindness);
`DISCHARGE_PACKAGE/README.md` §§7–8.
**Computations:** `theory-audit/t3_closures.py` (deterministic, numpy/scipy quad;
all numbers below reproduced by it; every closed form checked to ≥ 9 digits, the
three new exact identities to 12–14 digits).

**Verdict vocabulary:** SOUND / GAP / DEFECT-CANDIDATE / NOT-AUDITABLE-FROM-TEXT,
as in T1/T2.

---

## 0. What a repaired closure must specify, and the one identity that organizes all of them

T2 established that the closure's *kinematic skeleton* is sound and
configuration-space-free: E_rot/E_tot = j/2 (uses only L = j𝔠 and the clock),
w·j(1−j) = 1 selection, the M.0 loop κ = √2ê/𝔠, the ±½ error-signal exponents.
What F-R4/F-R5 broke is the *specification of the configuration space*: the paper's
(V, g)-family presupposes an energy-flat inertia dial that its own SDiff invariance
forbids (F-R4), and its enlarged space contains a tilt/ring halo channel that makes
every printed tuple a saddle (F-R5). A **well-posed repaired closure** therefore
consists of the unchanged two conditions

- (C-spin) L = 𝔠/2, (C-clock) E_tot = 𝔠ω — equivalently the engine fixed point
  L² = (2/3)·𝕀·E_static —

**plus an explicit, declared configuration class 𝒦 (and/or stability constraint)**
on which the fixed-L Routhian minimization is bounded and attained. The candidates
audited here are the four honest choices of 𝒦.

**The organizing identity (new, this audit; verified to 10⁻⁹).** From the clock
fixed point alone, κ_ours² = (2/3)E_static/𝕀 = E_tot/(2𝕀); pushing through the
validated unit map (ê = √2E_engine/π, 𝔦 = I_engine/(2√2π²), both checked against
ê₀ = 64/15π and 𝔦₀ = 256/105π):

**κ_paper² = ê_tot / (2 𝔦_tot)** — for ANY closure solution, any 𝒦. (T3.0)

In-family this reproduces κ² = 7/8g exactly (ê_tot = √2ê₀ g-free, 𝔦 = 𝔦₀gV,
V = √2: computed 0.875/0.583 at g = 1, 3/2 ✓). It is the honest generalization of
the S4′ stake, and it says at a glance what each candidate does to κ.

---

## 1. C1 — the hedgehog/rigid-direction-restricted closure (the scaling rung)

### 1.1 Formalization

𝒦₁ = {hedgehog ansatz: q = (cos f(r), sin f(r)·x̂), isotropic dilation}. The
closure on 𝒦₁ is Step 1/2's radial problem — well-posed (coercive in the dilation
variable, 1-D profile ODE), solution reproduced by the campaign to 7 digits.

### 1.2 Selected constants (all exact)

**𝔠 = 2√(ê₀𝔦₀) = 256/(15√7 π) = 2.053288; κ = √(7/8) = 0.935414; g = 1;
V = √2; E_rot/E = ¼; S4′ κ²g = 7/8 holds; mass law E_tot = √2ê₀Λm̃ unchanged.**

### 1.3 The strongest in-family extension — the corpus's own G.5 spheroidal ansatz, honestly costed (new computation)

The natural objection to C1 is "the corpus never claimed the hedgehog; it claimed
the spheroidal g(λ) family." So this audit built that family *with its true energy
cost*: F(x) = f₀(R\*s), s the ellipsoid radius (a = R\*λ^{−1/3}, c = R\*λ^{2/3},
volume-preserving), direction re-locked Θ = θ (the only reading under which g ≠ 1
is reachable at all — F-R4). For this locked transport family the sectors reduce in
closed form (derived in `t3_closures.py` header; b = sin²F·F_r/(2π²r²) because
Θ = θ makes dΘ∧dθ ≡ 0):

- **E₀(λ) = E₀(1) exactly** (value function + volume preservation);
- **E₆(λ) = E₆(1)·P(λ)**, P(λ) = ½∫sinθ·m̂³dθ ≥ 1 (Jensen; sphere strict min),
  m̂² = γsin²θ + cos²θ/γ², γ = λ^{2/3} — the honest cost the paper's
  "energy-flat at ε = 0" premise sets to zero;
- **𝕀(λ) = 𝕀(1)·Q(λ)**, Q = ¾∫sin³θ·m̂⁻³dθ — and Q(λ) **reproduces the corpus's
  corrected g(λ) (AUD-15 F-A15-7 oblate branch) to 6+ digits at every tested λ**.
  This is an independent argument-level validation that the corpus's g(λ) is
  correct *geometry* attached to a false *energetics*.

Closure on the family (fixed-L min over (λ, d), then clock) collapses to:
u\* = √(2A/B), L² = 2AC, and λ-stationarity ⟺ **(E₆/𝕀)′(λ) = 0** — the closure
selects the minimum of the cost-to-inertia ratio, not the inertia supremum.
Computed optimum:

| λ\* | g\* = Q | E₆ cost P | 𝔠\* | κ\* | V\* | E_rot/E | κ²g |
|---|---|---|---|---|---|---|---|
| **0.817257** | **1.078096** | 1.034796 | **2.168732** | **0.900897** | 1.438608 | 0.250000 | **0.875000000000** |

Three facts of record: (i) the honest G.5-family endpoint is **2.17, not 2.515** —
13.8 % short, −2.2σ even against the ε = 0.05 benchmark; the g-dial buys 5.6 % of
𝔠, not 22.5 %. (ii) **κ²g = 7/8 is an exact identity on the entire locked
transport family** (proof: κ_paper²·g = 4πE₀(1)/𝕀(1) = 7/8 — E₀ invariance makes
g cancel; verified to 12 digits) — so S4′ is *family-robust* but thereby also
*value-free within C1*: it cannot distinguish g\* = 1.08 from g = 3/2. (iii) the
per-λ fixed-point curve 𝔠(λ) does pass through 2.515 (near λ ≈ 0.59), but that
member is not λ-stationary — 𝔠₀ is not selected by anything.

### 1.4 Verdict

**SOUND as a variational problem** (bounded, attained, replicated), with **two
GAPs the corpus would have to own**: (a) no printed dynamical principle justifies
excluding the tilt/F-halo channels — the restriction is a solver artifact promoted
to physics; (b) every C1 solution is halo-unstable in the enlarged space —
measured: the rung is above even the *uniform locked* threshold (axi4 §4,
x = 1.080), and the honest spheroidal optimum κ\* = 0.9009 is above the ring
threshold 1/√2. C1 is internally consistent and externally indefensible.

---

## 2. C2 — the stability-constrained closure: which channel binds, and what it selects

### 2.1 Formalization

Minimize the fixed-L Routhian over the **unrestricted** degree-1 class subject to
the stability side condition κ_paper ≤ κ_crit(channel), then close with the clock.
The channel menu is the campaign's exact family (axi4 §3, engine-verified to
≤ 3×10⁻⁵): **κ_crit(w) = 1/√(2⟨sin²θ⟩_w)** — polar s²=|cos|: 1 exactly; uniform:
√3/2; s = sin²: √(7/12) exactly; equatorial-ring limit: 1/√2; unlocked tilt: 1/√2.

### 2.2 Which channel binds

The binding constraint is the **infimum of the menu, κ_crit = 1/√2**, reached by
two independent channels (unlocked far-field tilt, F-R5; and the locked equatorial
ring limit, 4C — the lock does not remove it). Imposing only a *weaker* channel's
constraint does not produce a well-posed problem: the functional stays unbounded
below through the stronger channel — measured fact, axi4 §5 ("the G_L minimisation
over the full locked family is not well-posed — finite-amplitude equatorial rings
drive it to −inf"). So the only self-consistent single constraint is κ ≤ 1/√2.

### 2.3 Theorem T3.1: the constrained optimum IS the saturation point, exactly

At fixed L, a marginal far halo in channel w adds statics dE₀ = dI/(16π⟨sin²θ⟩_w)
and rotation d(L²/2I) — so with the binding channel (⟨s²⟩ = 1) the halo ledger is
*linear*: E_static = G + I/(16π) with G ≡ E_static − I/(16π) halo-blind. The
fixed-L objective R(I) = G + I/(16π) + L²/(2I) is strictly convex in the halo-fed
inertia I with interior stationarity at

  1/(16π) = L²/(2I²) ⟺ κ_ours = 1/√(8π) ⟺ **κ_paper = 1/√2 exactly.**

The *unconstrained* fixed-L minimum therefore already sits exactly ON the
constraint boundary: the KKT multiplier is zero, the constraint is marginally
active, and adding it changes nothing. (With a strict inequality κ < 1/√2 the
infimum is not attained; its closure is the same point.) Substituting I\* and the
clock L² = (2/3)IE_static gives, for any core, **L = √(8π)·G and
𝔠_paper = (4/π)·G**; minimizing over cores gives G → G\* = 16√2/9 (G3's sharp
bound + attainer). The general-channel version 𝔠 = (4/π)√⟨s²⟩·G_w,
κ = 1/√(2⟨s²⟩) — derived here — is validated against the engine: it reproduces
axi4's printed uniform-channel pair (G_L\* = 2.423818 → 𝔠 = 2.519791 vs printed
2.51979, rel. 3×10⁻⁷).

**Verdict: SOUND, and C2 ≡ C3.** The stability-constrained closure selects
*exactly* the saturation point — not approximately. Within this model there is no
intermediate "stabilized" tuple between the rung and saturation for the corpus to
land on.

---

## 3. C3 — the saturated closure (the tuple the corpus's own equations select)

### 3.1 The complete closed-form solution

Configuration: the **oblate BPS compacton** g(cos(F/2)sinθ) = ρ³/(6√2)
(direction-locked, reducing to the campaign's spherical compacton on the axis,
G3-verified at 10⁻¹¹–10⁻¹⁴) **plus a marginal equatorial far halo** whose
amplitude is pinned jointly by saturation (κ_ours = 1/√(8π)) and the clock
(L² = (2/3)IE_static) — verified: both identities hold on the G3 tuple to 10⁻¹⁵.

| quantity | closed form | value | corpus (4.10) |
|---|---|---|---|
| 𝔠_phys | 64√2/(9π) | **3.201125** | 2.514754 |
| κ_phys | 1/√2 | **0.707107** | 0.763763 |
| binding depth 1−κ | 1−1/√2 | **29.29 %** | 23.61 % |
| ω_th/ω₀ = 2κ | √2 | **1.414214** | √(7/3) = 1.527525 |
| E_rot/E | 1/4 | exact | exact (survives) |
| g_core = 𝕀_core/𝕀₀ | **35/24** | 1.458333 | (g\* = 3/2 claimed) |
| g_tot = 𝕀_tot/𝕀₀ | **35/12** | 2.916667 | — |
| ê_core, 𝔦_core | 40/(9π), 32/(9π) | ratio **5/4** | ê₀/𝔦₀ = 7/4 |
| ê_tot = 𝔦_tot | 64/(9π) | 2.263537 | — |
| mass law E_tot/(Λm̃) | √2ê₀ × **5/(3√2)** = 64/(9π) | +17.85 % | √2ê₀ = 1.920675 |
| support volume ratio | **3/2 exactly** | 1.500000 | V = √2 = 1.414 |
| S4′ replacement | κ²g_tot = **35/24**; κ²g_core = 35/48 | 1.4583 / 0.7292 | κ²g = 7/8 |

New exact identities established by this audit (each verified to 12–14 digits in
`t3_closures.py`):

1. **Saturation ⟺ ê_tot = 𝔦_tot.** By (T3.0), κ² = ê_tot/2𝔦_tot; κ = 1/√2 is
   precisely the point where the total (paper-unit) energy and inertia numbers
   coincide, both = 64/(9π). This is the cleanest possible restatement of the
   repaired closure for the corpus's algebra chapter.
2. **g_core = 35/24, g_tot = 35/12.** The true minimiser's core realizes 97.2 % of
   the G.5 supremum 3/2 — but with ê_core/𝔦_core = 5/4, not 7/4: the "7" that
   generates √(7/8g), √(7/12) and √(7/3) does not survive off the spherical
   backbone (see §4(v)).
3. **The support volume of the oblate compacton is exactly (3/2)× the spherical
   compacton's** (equivalently ∫₀^π g(sinθ)/sin²θ dθ = 2; quad-verified to
   1.3×10⁻¹⁵). The in-family invariant V = √2 does **not** survive C3 (the
   solution is d-stationary at d = 1); the natural volume diagnostic is 3/2. A
   bench measuring "V" would be testing a different number — 1.500 vs 1.414 —
   and the ⟨r1⟩ measured 1.409 ± 0.010 sits 9σ from it.

### 3.2 Honest caveats (the price of C3)

- **Marginal stability, not strict.** The halo direction is exactly flat at the
  solution: the minimum is degenerate along a manifold (core unique; halo pinned
  only in its ∫η²-moments, its shape/radius free). The corpus would have to sell
  this zero mode as a modulus; radiative "silence" at threshold is marginal, and
  the consistency ring's third link (T2 §d.3) is repaired only in the weak sense
  "no strict descent," not "gapped."
- **The ε-law flips sign.** Measured: saturated 𝔠 rises with ε
  (3.201125 → 3.38547 at ε = 0.05, +5.76 %), against (4.8)'s falling
  𝔠₀[1 − 0.42ε^{2/3}]. Sec. IV.J's anchor arithmetic dies with the old law
  (see §5).
- **The convention question travels with it.** κ = 1/√2 is stated in the closure
  algebra's own normalization (κ_paper = ω/(√2μ), pinned by the campaign's unit
  map); the corpus's onset-quote convention must be reconciled per axi3 §3 —
  under either resolution C3 is the consistent tuple, but the App. I.4 control-run
  sentence must be re-worded to name its channel (s² = |cos|, whose exact onset is
  1.00000).

**Verdict: SOUND and well-posed — the unique reading under which the corpus's
stated variational problem has its stated infimum attained, now fully closed-form.
This is the tuple the corpus's own equations select; adopting it is the minimal
mathematically honest repair.**

---

## 4. C4 — the steelman: is 2.515 recoverable under ANY reading?

Every route this audit could construct or find in the record, each tested:

**(i) 𝔠₀ = G\* (the value coincidence).** 2.5147536 vs 2.5141574: **disproven
exactly** (G3): algebraic vs transcendental, gap 5.96×10⁻⁴ pinned; the near-miss
is 24√21 ≈ 35π (2.37×10⁻⁴). Also a category error: G\* is an objective value; the
closure maps it to 𝔠 by (4/π).

**(ii) The uniform-channel saturated reference.** (4/π)√(2/3)·G_L: at ε = 0.05 it
lands 2.51979 (+0.20 % of 𝔠₀) — and at t = 0 it lands 2.37096 (the benchmark's
2.37 ± 0.09; DISCHARGE §7 identity 3). Tempting, but: the uniform channel is not
binding (§2.2 — the functional is unbounded below through rings, so this
"closure" is a d-relaxed fixed point, not a minimum), and the 𝔠₀-proximity is
ε-specific. Coincidence-class, twice over.

**(iii) The zero-cost spheroidal reading (the corpus's actual derivation).**
Taking the G.5 family with E(λ) ≡ E(sphere) and g → sup = 3/2 gives
𝔠 = 2.053288·√(3/2) = **2.514754 = 𝔠₀ exactly** (verified). This is the honest
formalization of what Sec. IV.C–E does. It fails as mathematics on two
independent counts: (a) the supremum is an **open endpoint** — "approached, not
attained" by the corpus's own G.5 parenthesis — so 𝔠₀ is the sup of a family of
non-solutions, not the optimum of any functional; (b) the zero-cost premise is
**refuted by the F-R4 lemma at argument level**: 𝕀 = ∫2(q₁²+q₂²) is a value
function, so *no measure-preserving pullback of any kind* can change it — reaching
g > 1 requires changing field values, which costs E₆ (this audit's P(λ) ≥ 1 with
equality only at the sphere is the quantitative form). The honest cost turns
2.515 into **2.169** (§1.3). There is no limit, scaling, or convention in which
the cost vanishes while the gain survives: both live in the same ε-independent
sectors (E₆ vs 𝕀).

**(iv) The κ-side.** κ₀ = √(7/12) IS recoverable exactly — as the locked halo
onset of the s = sin² channel (⟨sin⁶⟩/⟨sin⁴⟩ = 6/7; axi4 exact identity). But
that is a *threshold*, not a closure solution, and its coincidence with κ(g = 3/2)
has a visible common arithmetic root (both are seventh-order sine-moment ratios);
nothing propagates from it to 𝔠₀.

**(v) The g = 3/2 shape-exact algebra composed with the true (oblate-compacton)
geometry — the commissioned test.** It does not compose: plugging the true
minimiser's g_core = 35/24 into the family algebra gives
𝔠 = 2√(g ê₀𝔦₀) = 2.479581 and κ = √(7/8g) = √(3/5) = 0.774597 — neither 𝔠₀ nor
the C3 values. The family algebra presupposes ê/𝔦 = 7/4 (the spherical Haar
moments); the true minimiser has ê_core/𝔦_core = **5/4** and ê_tot/𝔦_tot = **1**.
The √42 = √(2·3·7) in 𝔠₀ is inherited from the 7 in 105 = 3·5·7 of the spherical
𝔦₀; the saturated geometry replaces that 7 by 1 (whence the pure √2 of G\*). **No
closed-form identity connects the two constants** (G3 §4.3, confirmed here at
argument level): the g = 3/2 algebra and the oblate geometry are two different
theories of the same functional, and the functional agrees with the second.

**(vi) Numerology registry (flagged, not charged):** R_eq = (3√2π)^{1/3} =
2.370984 vs the benchmark 𝔠 = 2.37(9) — dimensional mismatch; the per-λ curve of
§1.3 passes through 2.515 at a non-stationary λ ≈ 0.59.

**Verdict C4: DEFECT-CANDIDATE for the claim that any well-posed within-model
functional selects 𝔠₀ = 2.5147** — every candidate reading is an idealization
whose load-bearing premise is refuted (iii), an exact-form disproof (i), a
non-binding channel accident (ii), or a threshold masquerading as a solution (iv).
One reading remains NOT-AUDITABLE-FROM-TEXT: that the frozen ⟨r1⟩ code contains an
explicit constraint term realizing (iii)'s idealization as an honest penalty whose
ε → 0 limit lands 𝔠₀; the code is absent from the archive (DISCHARGE §8 route (a)
stands open). Absent that, 2.515 is not recoverable.

---

## 5. Downstream propagation (everything below 𝔠₀), per candidate

### 5.1 The constants table

| quantity | corpus (4.10) | C1 (rung) | C1′ (honest G.5 family) | C2 = C3 (saturated) |
|---|---|---|---|---|
| 𝔠_phys | 2.514754 | 2.053288 (−18.4 %) | 2.168732 (−13.8 %) | 3.201125 (+27.3 %) |
| Λ√J = ħ/𝔠 (recalibration ×) | 1 | ×√(3/2) = 1.224745 | ×1.159548 | ×6√21/35 = 0.785584 |
| κ_phys | 0.763763 | 0.935414 | 0.900897 | 0.707107 (−7.42 %) |
| binding depth 1−κ | 23.62 % | 6.46 % | 9.91 % | 29.29 % |
| ω_th = 2κ·ω₀ | √(7/3) = 1.5275 | √(7/2) = 1.8708 | 1.8018 | **√2 = 1.4142** |
| interaction range (×ƛ_C) | 0.7638 | 0.9354 | 0.9009 | 0.7071 |
| S4′ κ²g | 7/8 | 7/8 ✓ | 7/8 ✓ (exact identity, §1.3) | **35/24** (g_tot) / 35/48 (g_core) |
| V | √2 | √2 ✓ | 1.4386 | d\* = 1; support ratio **3/2** |
| E_rot/E | ¼ | ¼ | ¼ | ¼ (survives everywhere) |
| mass-law factor E_tot/(Λm̃) | √2ê₀ | √2ê₀ (unchanged — E_tot is g-free, Cor. IV.3′) | √2ê₀·√P·(…) ≈ +1.7 % | ×5/(3√2) = +17.85 %, = 64/(9π) |
| ε-law (4.8) | 𝔠₀[1 − 0.42ε^{2/3}] | sign FLIPS: −1.03ε^{1.0} deficit (Step 2) | flips (same family class) | flips: +1.15ε class (3.2011 → 3.3855 at ε = 0.05) |
| B-U1′ ceiling on ε_e | 6×10⁻⁴ | ≈ 2.9×10⁻³ | ≈ 2.9×10⁻³ | ≈ 2.6×10⁻³ |

Notes: ħ = 𝔠Λ√J calibrates only the *product* Λ√J (the Λ-vs-J split is fixed
elsewhere); the entrainment ceiling uses |δ𝔠/𝔠| ≤ Δ_lock = 3×10⁻³ with the
repaired (linear, sign-flipped) laws — the admission window survives numerically
([1.5×10⁻⁶, ~3×10⁻³]) but its exponent-3/2 derivation and the "𝔠 falls toward
𝔠₀" narrative do not. ω_th = √2ω₀ under C3 also displaces every ω_th-inheriting
guidance number (e.g. WS-A8's line-spectra envelope).

### 5.2 The family sector is 𝔠₀-blind (P-O1's argument, formalized and verified)

Claim to verify: no candidate shifts the tier-5-replicated family/quark/neutrino
arithmetic. Formalized: the corpus's mass law is M_kc² = 𝒢·ê_k·Λm̃_k with 𝒢 a
**species-universal geometric factor** (√2 in-family; 5·√2/3·(ê₀-normalized) —
i.e. E_tot = 2G\*[k] — under C3): universality holds because the dimensionless
closure problem is m̃-free (Theorems IV.1/IV.3, whose m̃-cancellation the campaign
verified and F-R4/F-R5 do not touch). Hence **M_k/M_j = (ê_k m̃_k)/(ê_j m̃_j)** —
ratios of potential-sector Haar numbers and m̃-generator outputs, with 𝒢 and Λ
cancelling. Every tier-5 quantity (⟨r7⟩/⟨r8⟩ Σ(p) structure, ⟨r11⟩ belts, the
0.05σ/0.43σ pulls, the S1 pipeline) is built from such ratios: **blind to all of
C1/C2/C3/C4** — confirming `family_ADJUDICATION.md`'s statement at argument level.
The *non-blind* set, for the record: everything in the table above, plus any
absolute Λ-unit quantity and the κ-dependent edge-stiffening exponent's numerical
anchor. The α-drift/μ-drift *correlation structure* of P-O1 survives (it needs
only "mass ratios 𝔥-blind", which is exactly what was just verified); its
magnitude ~10⁻² (binding fraction) would rescale with the new binding depths.

### 5.3 What Sec. IV / App. G would need to change

**Under C1** (not recommended — unstable and unprincipled, §1.4): IV.C–E rewritten
to define the closure on the declared hedgehog(+dilation) class; Theorem IV.4 and
the bracket (4.6) demoted to "supremum of an inertia diagnostic, closure-inert";
G.5 struck (F-R4); rows 5–6 of the IV.D table struck or re-derived; (4.8)–(4.10)
replaced by rung values with a sign-flipped ε-law; ⟨r1⟩ re-explained as a saddle
artifact; App. I.4 must state the restriction and re-name its onset channel.

**Under C2 = C3** (the minimal honest repair): IV.C–E replaced by: the halo
threshold theorem κ_crit(w) = 1/√(2⟨sin²θ⟩_w) (new App. G content — it absorbs
the corpus's own onset measurement as the s²=|cos| channel, exactly); Theorem
T3.1 (constrained optimum = saturation, KKT-marginal); G3's oblate compacton as
the exact minimiser with the tuple of §3.1; (4.10) becomes 𝔠 = 64√2/(9π),
κ = 1/√2, depth 29.3 %, ω_th = √2ω₀, S4′ → κ²g_tot = 35/24 (or, cleanest,
"ê_tot = 𝔦_tot at saturation"); the IV.D table gains a final row "saturated
(well-posed limit)"; IV.F's cloak/ring paragraph re-argued at marginal (not
gapped) stability; the (4.8) ε-law re-measured with the sign flip and (4.9)
re-derived (ceiling ≈ 2.6×10⁻³); App. I.4's protocol amended to declare the halo
channel and the benchmark re-labelled "restricted-branch code-validation point."
Sec. VII.B's mass normalization picks up 5/(3√2) (absorbed by recalibrating Λm̃;
no family observable moves, §5.2).

---

## 6. Cross-cutting summary

| candidate | selects (𝔠, κ) | well-posed? | verdict |
|---|---|---|---|
| C1 hedgehog-restricted | (2.053288, 0.935414) | yes (attained) | SOUND as math; GAP — restriction undeclared/unjustified; solution halo-unstable in every enlarged space |
| C1′ honest G.5 spheroidal family (new) | (2.168732, 0.900897), λ\* = 0.8173, g\* = 1.0781 | yes | SOUND as math; same GAPs; proves the honest g-dial buys 5.6 %, not 22.5 % |
| C2 stability-constrained | = C3 exactly (Theorem T3.1; binding channel ring/tilt at 1/√2; weaker channels ⟹ ill-posed) | yes | SOUND; collapses to C3 |
| C3 saturated | (64√2/(9π), 1/√2) = (3.201125, 0.707107); oblate compacton + pinned marginal halo; g_tot = 35/12; ê_tot = 𝔦_tot = 64/(9π); support ratio 3/2 | yes (attained; degenerate flat halo manifold) | **SOUND — the unique well-posed reading of the corpus's own problem; the repair to adopt** |
| C4 any 𝔠₀-recovering reading | 2.514754 | no route found | DEFECT-CANDIDATE (routes (i)–(vi) each disproven/accidental); NOT-AUDITABLE only via the absent frozen ⟨r1⟩ code (discharge route (a)) |
| P-O1 family blindness | — | — | SOUND (mass ratios are (ê·m̃)-ratios; 𝒢, Λ universal ⟹ cancel) |

**Uncertainty, printed honestly.** (i) All C2/C3 statements inherit the campaign's
halo algebra and compacton unit map — validated at 7–9 digits against three
independent corpus quantities and engine-confirmed (G3 ladder budget 1.5×10⁻⁶),
but a convention error there would move §§2–3 wholesale; the places to attack are
listed in axi3 §5(b). (ii) The C1′ computation is this audit's own construction
(closed-form reduction + scipy quadrature, self-checked: volume constraint = 1 to
10⁻¹⁵, Q ≡ corpus g(λ), κ²g = 7/8 to 10⁻¹², rung recovered to 9 digits at λ = 1);
it has not been run on the campaign's frozen engine — a cheap follow-on. (iii) The
support-volume identity ∫g(sinθ)/sin²θ dθ = 2 is verified to 1.3×10⁻¹⁵ but not
symbolically proven here. (iv) The C4 verdict is about the *auditable record*; the
corpus's frozen ⟨r1⟩ code could still discharge via route (a), which would convert
C4's DEFECT-CANDIDATE into a spec-omission GAP without changing which closure the
*stated* equations select. (v) Everything is within-model; no statement here bears
on nature, and "repair" means only: make the printed mathematics select its
printed numbers.
