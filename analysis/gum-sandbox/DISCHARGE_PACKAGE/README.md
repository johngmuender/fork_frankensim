# To the authors of the GUM corpus — Replication Campaign Findings Package

**Contents:** this cover document; `REPRODUCE.md` (exact reproduction
instructions for every load-bearing number); `CLAIMS.json` (the findings
ledger, machine-readable). All artifacts are cited by repository path and
stay where the campaign produced them — nothing is copied, everything is
checkable in place.

**Purpose:** to give you everything needed to either **discharge** or
**concede** findings **F-R4** and **F-R5** under your own rules — your
ledger classes, your named-doubt discipline, your print-your-own-defects
norm. The same rules bind this package: every claim below is
**within-model** (it concerns whether your printed numbers solve your
printed equations); **nothing here bears on nature**, on GUM's physical
viability, or on any empirical stake — S1's neutrino window is untouched
(its pipeline replicated cleanly).

---

## 1. What the campaign was

A replication-first campaign: your computational claims re-derived
independently, on an independent codebase (frankensim's certified interval
arithmetic and evidence machinery), from your printed equations and
conventions, with **zero fitted constants** wherever you claim zero
(one calibrated constant where you calibrate one). Tiers: the constants
gauntlet (Tier 0), the linear spectrum (Tier 1), the reduced and
field-level ħ-closure (Tier 2a/2b Steps 1–3), Born-rule relaxation
(Tier 3), the unrestricted 3-D field and certified diagnostics (Tier 4),
and the family/neutrino/color arithmetic perimeter (Tier 5a). Index and
per-tier adjudications: `analysis/gum-sandbox/REPLICATION_CAMPAIGN_STATUS.md`.

Everything deterministic, most of it certified: 52/52 gauntlet checks as
fail-closed `Certified<f64>` claims under a content-addressed BLAKE3
Merkle root with bit-identical replay; 13/13 certified 3-D diagnostics
(Tier 4B); 91 machine-checked gates across the four ported physics crates
plus 18 solver gates including the F-R5 halo referee. Golden roots and
commands in `REPRODUCE.md`.

## 2. The findings ledger (F-R1 – F-R8)

| ID | Severity | One-line statement | Discharge path |
|----|----------|--------------------|----------------|
| F-R1 | minor (print) | NR-D2 §1's printed KE coefficient is ×25 too large at its stated v = 10⁻³c (2×10⁻⁴f printed vs 8×10⁻⁶f computed) — conservative direction; the murk kill fires ~25× harder | Correct the printed coefficient; the kill is unaffected |
| F-R2 | minor (print) | AUD-15 V15.2's quoted g(0.37–0.47) span [1.28, 1.32] is inconsistent with its own corrected formula (true: [1.257, 1.315]) | Reprint the span; conclusion unaffected |
| F-R3 | minor (print) | ⟨r10⟩'s ±0.90 band corresponds to ρ_eff = −0.36, not the printed ρ = +0.45 (which gives ±1.344) | Pin the covariance sign convention; S1 window unaffected (physical floors) |
| **F-R4** | **substantive** | App. G.5's SDiff inertia-enhancement mechanism fails in its literal reading: the isorotation inertia is invariant on the SDiff orbit (two-line proof + 10⁻¹⁵ measurement); the free g-dial to 3/2 does not exist there | §5 below |
| **F-R5** | **substantive** | The closure's stated variational problem does not select 𝔠₀ = 2.5147 or the ⟨r1⟩ tuple: the fixed-L Routhian is unstable to a far-field halo whenever κ_paper > 1/√2; every quoted solution sits above the threshold; the saturated infimum obeys the rigorous bound 𝔠_paper ≥ 2√2 | §6–7 below |
| F-R6 | moderate (provenance) | The ⟨r7⟩ → E-F1 → ⟨r8⟩ narrative chain is arithmetically false as written (2.9 → ×2 → +14% = 6.61 ≠ 5.6; printed errors shrink where the chain forces growth); the coherent reading is that ⟨r8⟩ is an independent re-evaluation | State ⟨r8⟩'s actual pipeline and correct §VII.C/K.1's provenance sentence. Direction-neutral: the honest chain would *worsen* data agreement |
| F-R7 | minor (print) | App. K.2's "qR\* ∼ 10⁻¹⁷" is six orders from its own inversion (ln(1/qR\*) = 24.36 ⇒ qR\* = 2.6×10⁻¹¹); the inversion side is coherent | Fix the stray prose magnitude; no downstream number moves |
| F-R8 | minor (print) | The locked-bond comparison's printed errors give a 1.54–1.59σ pull where the body claims 1.3σ (and F-A15-5 says 1.4σ) | Print the unprinted denominator uncertainty or amend the grade |

Artifact paths and key numbers for each: `CLAIMS.json`.

## 3. What replicated cleanly — the credit side, stated first

The ledger above is short because most of what we tested **replicated**,
often to more digits than you quote. In fairness, and because your own
discipline demands it, the credit side in full:

- **Every arithmetic gate in the swept archive: 52/52 PASS** (Tier 0,
  `analysis/gum-sandbox/tier0-gauntlet/RESULTS.md`). Every exactness claim
  is exact as claimed (7/4, ¼, √2, 7/3, 8/35, the 𝔠₀ square identity);
  "two computations one number" confirmed (Haar ∩ compacton = 64/15π,
  width 6×10⁻⁵); the ⟨r10⟩ bridge lands m₃ = 0.0468 eV; the g_s = 2 → μ_B
  identity and the Unruh recovery check; the WS-A4 demand tables reproduce
  from one δ = 10⁻¹² window including the "thin corner ∼40" at 38.6;
  NR-D1b's W-eternal floor replicates with your G→8 rounding
  conservative-direction. **No arithmetic defect was found in the corpus
  at Tier 0.**
- **The four-branch linear spectrum, exactly** (Tier 1): B3's factor 4
  measured at 4.000000000; Theorem II.2's exact masslessness invariant
  under m_V → 8×10⁻¹⁵; and **M-1 exhibited with a sharper mechanism than
  your own telling** — your fatal call on the defect stands, by a more
  precise route.
- **The ħ-closure's kinematics, everywhere, at every tier**: E_rot/E = ¼
  to 10⁻⁹ (automatic, as Cor. IV.2 claims), V(ε→0) = √2 exact, the spin
  selection family (j = 1 pole, j = 3/2 impossibility), the scaling rung
  𝔠 = 2√(ê₀𝔦₀) = 2.0533 reproduced to 7 digits, your spheroid formula
  g(0.42) reproduced to **9 digits**, clock residuals ≤ 4×10⁻⁹. These
  survive every finding below.
- **The radial rung** (Tier 2b Step 1): all App. I.1 gates pass (Derrick
  virial 6×10⁻⁷, 500× inside your gate); and your tail-mass formula
  μ² = m̃²/(2a_ψ) was **blind-re-derived and confirmed to 2×10⁻⁶ against
  our own spec error** — the corpus vindicated against its replicator.
- **Born-rule relaxation** (Tier 3): near-exponential H̄ decay for
  M ≳ 4 (τ = 83/10.2/4.9 for M = 4/9/16), no revival of the particle
  distribution, the equivariance control flat at the noise floor; the one
  tension (pocketed M = 4) resolved **in your favor** by the generic-set
  rerun (τ = 4.48, r² = 0.987).
- **⟨r10⟩'s central value** and the ⟨r11⟩ **quark belts vindicated by full
  reconstruction** (Tier 5a): R = 1 + B_geo/A − B_tube/A reproduces your
  printed ratios to the error digits, the down-type sign flip is forced,
  PDG lands exactly with direct m_t = 172.5 GeV (MS-bar measurably fails —
  your scheme choice pinned), four pulls at 0.17–0.32σ as claimed. This is
  *stronger* than your own summary sentence.
- **The antifragility record.** Your fail-closed instincts worked against
  us too: of the defects found en route, several were **ours** and are
  printed (a quote-precision test bug; a provenance-charset schema
  violation caught by fs-checker; two over-strict Tier-1 tests; the
  Step-1 spec error above). Your "external hands must close" Watch-Mode
  gates (WS-A4, NR-D1b, NR-D2 arithmetic) are now all externally closed —
  in your direction.

**Certification totals:** 52/52 (gauntlet, golden root `e19d0bcb…`) +
13/13 (Tier-4B certified 3-D diagnostics, golden root `f88731af…`) + 91
machine-checked crate gates (fs-gum-field 25, fs-gum-topo 25,
fs-gum-cosserat 20, fs-gum-sde 21) + 18 solver gates with the halo
referee (fs-gum-statics) — all bit-replay verified.

## 4. Where the two substantive findings live

Both concern one object: **App. G.5 / App. I.4's variational closure**
— the mechanism and the minimization that produce 𝔠₀ = 128√42/105π =
2.5147, κ₀ = √(7/12), and the ⟨r1⟩ benchmark tuple (𝔠 = 2.37 ± 0.09,
κ = 0.802 ± 0.018, λ\* = 0.42, g\* = 1.31, V = 1.409). Nothing else in
the corpus is implicated: the family sector is 𝔠₀-blind (your own P-O1
argument), and the selection theorems replicate.

## 5. F-R4 — the SDiff invariance lemma

Your G.5 argument: the isorotation inertia, taken over volume-preserving
(SDiff) images of the static solution, approaches sup g = 3/2 in the
extreme-oblate limit while the (6+0) energy stays flat — the basis of the
deep-BPS endpoint 𝔠₀ = 2.515. The campaign's finding, first measured at
machine precision (drift 1.2×10⁻¹⁵, `tier2-closure/axi_results.json`,
V3_cross) and then proved, verbatim from
`analysis/gum-sandbox/tier2-closure/axi_FR4_NOTE.md`:

> **The invariance lemma (two lines).**
> Let T be a volume-preserving diffeomorphism (det DT ≡ 1) and Q′ = Q ∘ T
> the pulled-back texture. For **any** density that is a function of the
> field *value* alone — the potential density `1 − σ_P(Q)`, and the
> **isorotation inertia density** `2(q₁² + q₂²)` — the change of variables
> y = T(x) gives
>
> ∫ Φ(Q′(x)) d³x = ∫ Φ(Q(T(x))) d³x = ∫ Φ(Q(y)) d³y.
>
> The sextic sector obeys the same because the topological density
> transforms as `b_{Q∘T} = (b_Q ∘ T)·det DT = b_Q ∘ T`, so `∫ b² d³x` is
> likewise invariant. **Hence E₀, E₆, and 𝕀 are all constant on the SDiff
> orbit of any base configuration.** Only the ε-suppressed kinetic sectors
> (E₂, E₄) vary.

Note that G.5's own proof sketch — "∫sin²f is SDiff-invariant and
sin²θ ≤ 1, so sup g = 3/2" — uses precisely the invariance that freezes
𝕀: the weight sin²f·sin²θ *is* q₁² + q₂², a value function. The supremum
is real as an inequality over all configurations; the claim that SDiff
images *approach* it is what fails. The deformation class that does turn
the g-dial (contour deformation at held hedgehog direction — which
reproduces your own g(λ) formula to 9 digits) is not SDiff and not
energy-flat: transverse to the orbit the static Hessian is positive (your
own Derrick stability), so inertia gain carries an ε-independent
quadratic cost, and the compacton's full 12-direction marginal spectrum
was measured (every dE/dp = 0, every d²E/dp² > 0, every dI/dp ≠ 0 —
inertia is never free and never refused; `tier2-closure/axi3_RESULTS.md`
§6).

**Hardening (Step 3).** The natural repair — "G.5 is a wording defect;
non-SDiff modes supply the enhancement at finite cost and the numbers
survive" — was tested and fails: adding the direction-tilt modes (the
only modes that can turn the g-dial, by the lemma) drives the ε→0 closure
from the scaling rung (2.0533) **past** your endpoint 2.5147 without
stopping (§6). No wording repair of G.5 survives Step 3. In 3-D the
supremum reappears in its lawful role: the below-threshold control's
inertia self-limits at exactly the (3/2)·I ceiling — a bound, not a dial.

**What would discharge F-R4** (from `axi_FR4_NOTE.md`): (a) exhibit
non-SDiff modes whose cost/gain balance lands g\* = 1.31 at ε = 0.05 and
drives 𝔠(ε→0) → 2.515 — then G.5 has a wording defect and the numbers
survive (Step 3 measured the balance; it does not land there — see F-R5);
or (b) produce a non-compositional reading of "SDiff images" under which
G.5's derivation is sound. Absent both, the deep-BPS endpoint reverts
toward the scaling rung and every constant downstream of 𝔠₀ shifts at the
~20% level within the model.

## 6. F-R5 — the halo instability and the saturated closure

**The mechanism (derived, then verified, then confirmed in 3-D).** For
any configuration, add a low-amplitude far-field halo of amplitude η in
the tilt channel. Exact small-η energetics in the closure's own
conventions (E₀ = (1/4π)∫(1−cosF), I = ∫2sin²F sin²Θ; gradient costs
vanish in the soft-halo limit):

```
dR = dE_static − (L²/2I²)·dI = [∫η² dV] · ( 1/(8π) − κ_ours² )
```

The fixed-L Routhian **decreases** under halo growth whenever
κ_ours > 1/√(8π), i.e. **κ_paper > 1/√2 = 0.70711**. This is the
standard over-spin criterion ω > μ(tilt) of your own imported literature
([8]/[24]) — nothing exotic. Verified by engine-level shell probes:
κ_crit measured 0.2163–0.2169 (finite width) vs ideal 1/√(8π) = 0.19947,
dI quadratic (ratio 4.000), position-independent to 5×10⁻⁵.

**Placement.** Every solution the corpus quotes sits above the
threshold: the scaling rung (κ = 0.935), the ⟨r1⟩ benchmark (0.802), the
deep-BPS endpoint (√(7/12) = 0.764). The ⟨r1⟩ benchmark is a
**halo-unstable saddle of its own functional** — a stationary point with
a measured, strictly downhill fixed-L direction.

**Saturation, not runaway.** Growing the halo raises I and lowers
κ = L/I; the minimizing sequence self-limits at **κ_paper = 1/√2
exactly**, where the closure factorizes: 𝔠_paper = (4/π)·G\*,
G\* = min(E_static − I/16π), with the Bogomolny-type bound G ≥ π/√2
(quadrature-verified to 2×10⁻¹⁶) giving the **rigorous**
𝔠_paper ≥ 2√2 = 2.82843 — which **excludes 𝔠₀ = 2.5147 from below** in
the enlarged space. Family values: 𝔠(t→0) ∈ [2.828, 3.201]. The G3
precision push subsequently closed this in exact form
(`tier2-closure/gstar_RESULTS.md`): **G\* = 16√2/9 exactly** (a sharp
Bogomolny bound attained by a closed-form oblate compacton), so the
saturated closure value is **𝔠_paper = 64√2/(9π) = 3.2011247 exactly**
at κ_paper = 1/√2.

**3-D, assumption-free** (Tier 4A, `tier4-field/field3d_RESULTS.md`): on
a free 96³ Cartesian grid with no symmetry, no ansatz, the fixed-L
descent fell monotonically 0.560 below the axisymmetric stationary value;
κ passed through the threshold; the halo was box-limited (a 1.3× box
descends deeper — the runaway is physical). The mechanism fingerprint:
**88% of the inertia gain is priced in E₀ at exactly the
threshold-defining rate 1/(16π)**. Two-sided: the below-threshold control
formed no condensate and self-limited at the (3/2)·I tilt ceiling. And
**the clock condition forces the over-spun regime**: imposing your clock
on a B = 1 solution lands κ at 1.245–1.474× threshold in every reading
measured.

**Direction-locking does not rescue it** (Tier 4C,
`tier2-closure/axi4_locked_RESULTS.md`): hedgehog-locking the direction
field removes the tilt channel but leaves a one-parameter family of
profile-halo onsets κ_crit_paper = 1/√(2⟨sin²θ⟩_w) sliding from √3/2
(angle-uniform) down to the unrestricted 1/√2 (equatorial ring); the
converged locked closure runs away accordingly (𝔠 = 3.26 at ε = 0.05,
+9.9σ; the locked ε-scan also has the opposite sign to your
0.42·ε^(2/3) deficit law).

## 7. The exact identities — where your numbers appear to come from

The 4C channel analysis produced three identities, stated here as
measured facts (the reading after them is a flagged hypothesis):

1. **Your measured over-spin onset is a halo threshold, exactly.** The
   polar-weighted locked channel (s² ∝ |cosθ|) has ⟨sin²θ⟩_w = 1/2
   exactly, hence κ_crit_paper = **1.00000 exactly** — against your
   App. I.4 control-run quote "emission onset at κ = 1.000 ± 0.004"
   (smoothed-probe verification 0.98810–0.99089).
2. **Your deep-BPS κ₀ is a halo threshold, exactly.** The s = sin²
   channel has ⟨sin⁶θ⟩/⟨sin⁴θ⟩ = 6/7 exactly, hence κ_crit_paper =
   **√(7/12) exactly** — the same closed form as your κ₀ = √(7/12) =
   0.7638 (which you derive from κ²g = 7/8 at g = 3/2).
3. **Your benchmark 𝔠 sits at a saturation value.** The locked
   uniform-channel saturated closure on the hedgehog core lands at
   𝔠_paper = **2.37096** — against your ⟨r1⟩ benchmark 2.37 ± 0.09 (a
   one-part-in-300 match of a solve nobody tuned).

A fourth proximity was recorded here without a claim — G\* converging
to within ~10⁻⁴ of 𝔠₀ — and has since been **resolved as an accident**
(`tier2-closure/gstar_RESULTS.md`): G\* = 16√2/9 = 2.5141574 exactly
(algebraic), 𝔠₀ = 128√42/(105π) = 2.5147536 (transcendental); the gap
|G\* − 𝔠₀| = 5.96×10⁻⁴ is pinned in closed form, and the 0.024%
near-miss is the coincidence 24√21 ≈ 35π. We print this correction
against our own earlier speculation per the shared discipline.

Separately, a convention tension that is a finding regardless of the
hypothesis: your closure-algebra κ's (√(7/8), 0.802, √(7/12)) and your
onset quote (κ = 1.000) **cannot both hold in one convention** for the
tilt channel. Under the compacton-exact unit map (validated by
reproducing your scaling-rung tuple to 7 digits and g(0.42) to 9),
κ_paper = ω/(√2μ), so the physical onset ω = μ sits at κ = 1/√2 — and if
instead your κ = ω/μ, translating the closure-algebra values consistently
still places rung/benchmark/endpoint at ω/μ = 1.08–1.32 > 1. Either way
at least one printed claim fails within-model.

**The coherent reading (hypothesis, flagged as such):** the ⟨r1⟩
pipeline's reported numbers are thresholds and saturation values of the
halo family its own solver could not represent — the solver converged to
(or drifted along) the edge of an instability it was blind to, its onset
diagnostic measured the one channel it could see (polar), and its
benchmark landed at the uniform-channel saturation point. Identities 1–2
are rational-moment facts of the sphere and could in principle be
coincidences of sin²-moment algebra; identity 3 is numerical. Only your
frozen code can decide.

## 8. Precise discharge conditions for F-R5

Any one of the following discharges the finding:

**(a) Produce the frozen ⟨r1⟩ code.** Show the term, constraint, or
restriction in the *actually solved* functional that gaps the halo
channel above the clock frequency (candidates the campaign could not
distinguish: a direction lock **plus** a radial-only or
equator-suppressed far field; a compact computational domain acting as a
box constraint; a stability side condition; or Newton convergence to the
saddle). Then F-R5 becomes a **spec-omission finding**: the benchmark
survives as the solution of the *constrained* problem, App. I.4 is
amended to state the constraint, and G.5's ε→0 endpoint still falls
separately under F-R4.

**(b) Identify the stabilizing constraint independently** — i.e., state
and justify, from the corpus's own text or an amendment, the physical
side condition under which the printed tuple is the true minimum, and
show it excludes *both* the tilt halo and the equatorial profile halo
(direction-locking alone was tested and does not suffice — Tier 4C).

**(c) Exhibit an error in the halo algebra or the unit map.** The three
named attack points, in decreasing order of leverage:
  1. **The E₀ coefficient**: 1/(8π) potential cost per unit ∫η² (from
     1 − cosη ≈ η²/2 with the 1/4π sector constant);
  2. **The inertia normalization**: KE = ∫(q̇₁² + q̇₂²), fixed by the
     same convention as I itself (dI = 2∫η² in the tilt channel);
  3. **The compacton calibration** (the unit map 𝔠_paper = 𝔠_ours/√(2π³),
     κ_paper = 2√π·κ_ours) — validated at 7–9 digits against three
     independent paper quantities, but it carries the √2 that the
     convention question in §7 turns on.

Every number entering these three is reproducible from the commands in
`REPRODUCE.md`; the small-η algebra is six lines and printed in full in
`tier2-closure/axi3_RESULTS.md` §3 and (locked channels)
`tier2-closure/axi4_locked_RESULTS.md` §3.

**If none of (a)–(c) is available:** within the model, the ħ-derivation's
value 𝔠₀ has no derivation from the printed equations; the honest answers
are the scaling rung with a stability constraint (𝔠 = 2.053, −18%) or the
saturated closure (≥ 2.828, +12% to +27%), and every constant downstream
of 𝔠₀ (the ħ calibration 𝔠Λ√J, κ_phys = √(7/12), the S4′ ε-run anchor,
binding depths, ω_th) inherits that shift. The selection theorems, the ¼
and √2 invariants, the clock algebra, the spectrum, the radial rung, the
Born-rule results, the family sector, and the S1 stake are untouched.

## 9. Epistemic notice (binding on this package)

Every claim above is within-model, made under your own rules: named
doubts with pre-stated discharge paths, defects printed on both sides
(including ours), severity graded, nothing adjudicated that your frozen
code can adjudicate better. Reproducible arithmetic is not evidence about
nature — in either direction. The campaign found your archive's
arithmetic essentially flawless (52/52), your instrument design
antifragile, and two substantive within-model defects in one appendix's
variational argument. We would consider the package a success under
exactly one condition: that you run it.
