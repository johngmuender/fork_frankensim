# H3.3 — PINNING THE UNDERDETERMINED: Lemma II.1 convexity (T-H8) and the E.6 response exponents (T-H10)

**Campaign:** theory Phase H3, workstream H3.3 (`ROADMAP_v5_THEORY.md`).
**Date:** 2026-07-16. **Code:** `h33_convexity.py` (14/14 internal checks PASS,
8 s) and `h33_pq.py` (21/21 PASS, 1.4 s); python3 + sympy/numpy, deterministic.
**Machine-readable results:** `h33_summary.json`.
**Sources:** `substrate-suite/01-GUM-Omega-Paper-v2.0.1.md` (Lemma II.1;
App. F.1; Thm VIII′.1; Secs. II.A–C, VI.A–C; App. E.6);
`theory-audit/T4_gravity_oneworld.md` §§G2, W7; `theory-audit/h25_RESULTS.md`.
**Epistemic rules:** within-model only; every algebraic claim below is either
machine-checked (check ID cited) or explicitly grade-flagged. Verdict
vocabulary as in T4. **Status: RESULTS-layer analysis; promotion is the
coordinator's call.**

## 0. Verdicts in two lines

- **T-H8:** App. F.1's convexity claim is **FALSE** — strict machine-verified
  counterexample (det F_t < 0 on an open interval, not merely a boundary
  touch); the lemma as printed is **ill-posed as to its admissible class**;
  the conclusion is true for the global-diffeo class (by a deep import the
  corpus did not cite for the job) and, as actually *used* by Theorem VIII′.1,
  true under an elementary **history-path repair** stated in §4. Pairs-only
  survives; the printed proof does not.
- **T-H10:** the corpus's own action **pins the sign after all — against
  E.6**. The displacement kinetic term forces **p = q − 1** for the photon
  doublet B2; the h2.5 positivity wedge then requires q ∈ (2/3, 5/6), and
  every textually consistent reading puts q on the half-integer lattice —
  so **𝔞₁^(B2) < 0 unconditionally** (−1/3, −2/15, or −5/33 by reading;
  −1/12 cone-only). E.6's premise w_s > 0 ∀s is false on the corpus's own
  action under *every* permitted reading: T-H10's conditional defect becomes
  **unconditional**. The rescue exists but is a new quantization-measure
  postulate (window of width 1/18 in density weight, containing no natural
  convention the corpus uses), priced in §8.

---

# PART A — T-H8: Lemma II.1 / App. F.1

## 1. The claim as printed

Main text (Sec. II.B):

> **Lemma II.1 (charge transfer and protection)** [DF]. For all admissible u,
> deg R̃[u] = 0 (Diff_c(ℝ³) contractible; degree homotopy-invariant), hence
> deg P̃ = deg Q̃ ≡ K: topological charge is measured in the relative texture
> and can change only where det F → 0 — *particle identity is guarded by the
> integrity of matter*, and (Sec. VIII′) creation is thereby localized to the
> coincidence stratum. ∎

App. F.1 (the printed proof):

> **F.1 Lemma II.1.** Admissible deformations (det F > 0, identity at ∞) form
> a convex, hence contractible, set; R̃[u] is a continuous map into
> Maps(ℝ³, SU(2)) from a contractible domain ⟹ deg R̃ ≡ deg R̃[0] = 0;
> deg P̃ = deg Q̃. Contrapositive: ΔΣK ≠ 0 ⟹ det F → 0 — the pairs-only
> theorem's engine.

## 2. Its role in Theorem VIII′.1

Sec. VIII′.B rests the nucleation sector's architecture on exactly the
contrapositive:

> **Theorem VIII′.1 (creation is pairwise)** [DF]. ΔΣK ≠ 0 requires det F → 0
> somewhere (Lemma II.1 contrapositive): **single-knot creation is forbidden
> in tear-free matter** … creation lives on the **coincidence/small-pair
> boundary stratum** of configuration space. ∎

Downstream of it: the IBC formalism's boundary stratum, UV-finite creation
rates, the selection rules, and (via the threshold 8′.2) the spectator-pair
phenomenology. The logical chain is: admissible set contractible ⟹ deg R̃[u]
constant = 0 ⟹ deg P̃ = deg Q̃ (since P̃ = R̃†Q̃ and degree is additive under
the product) ⟹ total charge can change only through inadmissibility
(det F → 0), which is the tear/coincidence stratum.

## 3. Formalization, readings, and the counterexample

**Formalized claim (as printed).** Let
𝒜 = {u ∈ C¹(ℝ³,ℝ³) : det(𝟙+∇u)(x) > 0 ∀x, u → 0 at ∞}. F.1 asserts:
u₀, u₁ ∈ 𝒜 ⟹ (1−t)u₀ + tu₁ ∈ 𝒜 for all t ∈ [0,1].

**Counterexample 1 (boundary touch — T4's W7, re-verified).** The twist map
φ(x) = R_z(θ(r))x with θ ≡ π for r ≤ 1, smoothstepped to 0 on [1,2]: each
sphere rotates rigidly, and det ∇φ = 1 identically for a *generic* profile
θ(r) — proved symbolically (check C1; the determinant of the rank-one update
R_z(θ)(𝟙 + θ′ J_z x r̂ᵀ-term) collapses because r̂ᵀJ_z x = 0). φ is a global
diffeomorphism, identity at ∞ — admissible on any reading. The segment to the
identity has, in the core, det((1−t)𝟙 + tR_z(π)) = (1−2t)², which **vanishes
at t = ½** (checks C2, C2b; T4's numeric row +1, ¼, 0, ¼, +1 reproduced,
C2c). The set is not convex — and since u₀ = 0 is one endpoint, it is **not
even star-shaped about the reference configuration** (the weakest reading of
"convex" that could serve F.1).

**Counterexample 2 (strict — new in this run).** The touch at det = 0 could
tempt a "measure-zero, perturb it away" defense. It cannot be perturbed away:
compose the twist with the stretch ψ(x) = x + e^{−r²/9}x₁ê₁. Then
(i) ψ is admissible: det ∇ψ = 1 + η + η′x₁²/r symbolically (C3a), with global
minimum 1 − 2e^{−3/2} = 0.5537 > 0 (C3b/C3c, grid + dense scan + analytic);
(ii) φ∘ψ is admissible by the chain rule (C3g), with
F(0) = R_z(π)·diag(2,1,1) = diag(−2,−1,1) (verified by finite differences,
C3d); (iii) the segment to the identity has det F_t(0) = (1−3t)(1−2t)
(C3e), **strictly negative on the open interval t ∈ (1/3, 1/2)** (C3f;
det = −0.04 at t = 0.4). The convex combination of two admissible
deformations is *orientation-reversing* on an open set — an open-condition
failure, robust to any perturbation of the endpoints. (Pointwise content:
GL⁺(3) is not convex; F.1 imported a convexity that fails already at a
single material point.)

**Reading enumeration (the audit pattern).**

| # | reading of "admissible … convex, hence contractible" | verdict |
|---|---|---|
| L1 | 𝒜 (pointwise det F > 0, identity at ∞) is convex — F.1 as printed | **FALSE** — strict counterexample above; not star-shaped either |
| L1′ | same 𝒜 is contractible (drop "convex", keep the conclusion) | **NOT-PROVEN, plausibly FALSE** — for the local-diffeo class (interpenetration allowed) the h-principle heuristic gives π₀ ≅ π₃(SO(3)) = ℤ (T4 W7 repair-2, sketch grade): a disconnected set on which deg R̃ plausibly labels components |
| L2 | admissible = *global* orientation-preserving diffeos decaying to identity (non-interpenetration, the physical reading); main text's "(Diff_c(ℝ³) contractible)" | **conclusion TRUE, printed proof still false** — the twist map *is* a global diffeo and the segment still exits; what carries the conclusion is π₀ Diff_c(ℝ³) = 0 (Cerf) / full contractibility (Hatcher 1983, the Smale conjecture) — a deep theorem, not one-line convexity. Note the corpus's two texts disagree: the main text cites the deep fact, F.1 "proves" it by false convexity |
| L4 | history form: along a continuous physical history u(t) with det F(t,·) > 0, deg R̃[u(t)] is constant | **TRUE (elementary)** — and it is all VIII′.1 uses |

The lemma is therefore **ill-posed as printed**: the text never fixes whether
"admissible" means local (det F > 0 only) or global (embedding) deformations;
the printed proof is false on both; the conclusion is true on the global
reading (deep import) and false-in-danger on the local one.

**Degree machinery sanity (C4).** The numeric degree integrator validates on
a hedgehog control (|deg| = 0.989 ≈ 1 at grid resolution, C4a), and along the
*admissible* path s ↦ R_z(sθ(r)) the SU(2)-lift degree integrand is
identically zero (C4b): the history-path argument operating exactly where the
convex path fails.

## 4. The repair

**Lemma II.1′ (corrected).** *Let t ↦ u(t) be a continuous history with
u(t) → 0 at ∞ (uniformly) and det F(t,x) ≥ δ > 0 on compacts. Then
deg R̃[u(t)] is constant in t; if u(0) = 0 (or u(0) is connected to 0 within
the admissible class), deg R̃[u(t)] ≡ 0 and deg P̃(t) = deg Q̃(t) throughout.
Contrapositive: ΔΣK ≠ 0 along a history ⟹ inf_x det F → 0 along it.*

*Proof.* Polar decomposition F ↦ R is continuous on GL⁺(3), uniformly on
det F ≥ δ; the SU(2) lift of a map ℝ³ → SO(3) trivial at ∞ exists and is
unique up to global sign (ℝ³ simply connected); Brouwer degree of the
compactified map S³ → S³ is locally constant under uniform convergence. No
convexity, contractibility, or even connectedness of the admissible set is
used. ∎

This is strictly weaker than what F.1 tried to prove and strictly sufficient
for Theorem VIII′.1, which quantifies over histories (nucleation events), not
over abstract pairs of configurations. **Pairs-only, the coincidence stratum,
the IBC construction, and the threshold phenomenology all survive under
II.1′.** Recommended erratum: replace F.1's convexity sentence with the path
argument (two lines); if the corpus wants the set-level statement too, it
must (a) define admissible = global embeddings and (b) cite Cerf/Hatcher at
import grade [IM], not [DF].

**T-H8 verdict: DEFECT-CANDIDATE confirmed and sharpened (printed proof
false, strict counterexample; statement ill-posed between two classes);
conclusion SOUND under the stated repair.** Concurs with T4-W7 and upgrades
its boundary-touch counterexample to an open-set orientation reversal.

---

# PART B — T-H10: deriving (p, q) from the corpus's own action

## 5. The question, and every constraint the text actually places

h2.5 established: every summed sector of (6.1) reduces to
D = −w_t(x)∂_τ² − ∂_i w_x(x)∂_i + M(x) with flat functional measure, and the
R-coefficient is the moduli-independent master form (independently
re-derived here from the Gilkey/Vassilevich formulas, checks M1a–M1e,
including the sign convention anchor M1b and the mass-response separation
M1d):

**𝔞₁(p, q) = −(5p + q) / (12(2p + q))**, w_x ∝ 1 + p h, w_t ∝ 1 + q h,

positive iff (5p+q)(2p+q) < 0 — an 8.5%-measure wedge of opposite-sign
response rays containing the unimodular ray q = −3p (anchors M2a–M2d: on
that ray E = O(h²), i.e. the flat measure *is* the covariant measure, and
𝔞₁ = +1/6). h2.5's verdict was conditional: "the (p,q) response is genuinely
underdetermined by the corpus." This run does the derivation the corpus
should have done. What the text fixes:

1. **II.C:** the moduli λ, μ, μ_c, α, β, γ, and the densities ρ₀, J are
   *constants* — numbers of the medium, campaign-validated at flat
   background; no defect-response function for any of them is printed
   anywhere in the paper or archive.
2. **VI.A [IM: KBKK]:** "disclination density **is** the curvature of an
   effective Riemann–Cartan connection." A weak static disclination
   background therefore enters as a non-flat **material metric** g₃ (and
   frame) — not as a new field the moduli could depend on. Probe class used
   here: conformal g₃_ij = (1+h)δ_ij (sufficient: a universal positivity
   claim that fails on conformal backgrounds fails).
3. **VI.B:** the background is "slowly varying defect-geometric"; the
   *displayed* sector perturbation is g_s = ḡ + 2δ_s c² n⊗n — a **cone-only**
   (time-time) modulation, i.e. (p,q) = (0,1).
4. **II.B (F10′, the M-1 lesson):** stored energy may depend on rotations
   only through the micro-frame's rotation *relative to the material frame* —
   objectivity ties the micro-rotation's algebra index to the material
   **frame** (soldered), not to coordinates.
5. **Flat functional measure** on the medium's fields u, φ (h2.5, forced by
   the Sec. III construction: the Fisher/Madelung chain is written in flat
   field coordinates).

What the text does **not** fix: the index variance with which W₂'s
contractions covariantize (u and the soldering map as covector / frame /
vector objects). That is the entire residual freedom — and it turns out not
to matter for the verdict.

**A preliminary pincer (reading "strain-only").** If one refuses
covariantization altogether — the defect enters only as a background strain
inside flat-space W₂ — then, W₂ being exactly quadratic with constant
coefficients, the fluctuation operator is *unmodulated*: (p,q) = (0,0), no
induced R-term exists, and (6.1) is empty — Sakharov gravity, Theorem VI.1,
and Sec. VI die structurally. The corpus therefore *needs* the covariant
reading for its gravity sector to exist at all. The covariant reading is what
we now evaluate.

## 6. The derivation

Covariantize each sector's quadratic action on g₃ = (1+h)δ with constant
moduli and measure √g₃ d³x (this much is forced by items 1–2; time is
untouched — the disclination background is static spatial geometry in the
medium rest frame). Every coefficient of the reduced branch operator is then
a product of exactly three kinds of factors (sympy series extraction, M3a):

| factor | power of (1+h) |
|---|---|
| √g₃ (measure) | +3/2 |
| each derivative-index contraction g₃^{ij} | −1 |
| each field-index pair: covector g₃^{ij} / frame δ_ab / vector g₃_{ij} | −1 / 0 / +1 |
| group-trace kinetic Tr[(Q̃⁻¹∂_tQ̃)†(Q̃⁻¹∂_tQ̃)] (intrinsic Killing form, **no spatial index**) | 0 |

Sector structures (from II.C/II.E, cf. the h2.5 sector table): B2's inertia
is ρ₀u̇u̇ (one field pair, no derivative pair) and its stiffness μe_(ij)e_(ij)
(one derivative pair + one field pair); B3/B4 and the bosonic knot band-edge
carry the Q̃ kinetic term (group trace — **metric-blind under every reading**;
the Killing form is basis-intrinsic) and wryness stiffness (one derivative
pair + one soldered/algebra pair). Results (all machine-extracted, M3a/M3b):

| reading of the spatial indices | B2 (p, q) | 𝔞₁^(B2) | B3/B4/knot_b (p, q) | 𝔞₁^(Q-sectors) |
|---|---|---|---|---|
| covector (unsoldered) | (−1/2, 1/2) | **−1/3** | (−1/2, 3/2) | **+1/6** |
| **frame / soldered (favored by F10′)** | (1/2, 3/2) | **−2/15** | (1/2, 3/2) | **−2/15** |
| vector | (3/2, 5/2) | **−5/33** | (3/2, 3/2) | **−1/6** |
| cone-only (VI.B's displayed g_s) | (0, 1) | **−1/12** | (0, 1) | **−1/12** |

**The two structural invariants (the new result):**

- **B2: p = q − 1, always.** The stiffness differs from the inertia by
  exactly one derivative-index contraction (−1), for *every* uniform variance
  choice and — since a field density-weight redefinition (det g₃)^w shifts p
  and q together by 3w — for **every density weight as well** (M4). The
  positivity wedge with p = q − 1 is q ∈ (2/3, 5/6) (M5a), an interval of
  width 1/6 that contains **no half-integer**; the reachable q are
  {1/2, 3/2, 5/2} (M5b). Hence **𝔞₁^(B2) < 0 in every textually consistent
  reading** — the photon doublet, the one sector whose operator the campaign
  validated at Tier 1, carries a negative weight no matter how the residual
  ambiguity is resolved.
- **Q̃-sectors: q = 3/2 is group-protected** (the kinetic trace has no
  spatial index), and only p moves with the reading: {−1/2, 1/2, 3/2}. The
  wedge at q = 3/2 is p ∈ (−3/4, −3/10) (M5c), which contains exactly one
  reachable value: p = −1/2, the **unsoldered covector** reading (M5d),
  landing precisely on the unimodular ray with 𝔞₁ = +1/6. But that reading
  identifies the algebra index with a *coordinate* covector rather than the
  material frame — violating F10′'s frame-relative discipline (the corpus's
  own cured defect M-1), abandoning the KBKK soldering that VI.A itself
  imports, and (declared limit, §9) generating soldering-gradient connection
  terms the two-exponent master form does not capture. And it cannot rescue
  B2 in any case.

## 7. Verdict: the action pins the signs — against E.6

Assembled sign map (with the h2.5 knot-statistics readings carried along):

| reading | B2 | B3/B4/knots_b | Σ_s w_s | consequence for Thm VI.1 / E.6 |
|---|---|---|---|---|
| frame/soldered (textually favored: F10′ + KBKK) | −2/15 | −2/15 | < 0 | hull formally survives (common sign divides out of (6.2)) but **1/16πG = Σw_s < 0: inverted Newton constant**, contradicting VI.B's G ∼ +c³/(ħN_effΛ_UV²) |
| vector | −5/33 | −1/6 | < 0 | same: G_ind < 0 |
| covector (violates F10′) | −1/3 | +1/6 | indefinite | **mixed signs: c_GW² can exit the hull; the slaving bound (6.2)/(6.3) fails as an inequality** |
| cone-only (VI.B displayed) | −1/12 | −1/12 | < 0 | G_ind < 0 (h2.5's row, confirmed) |
| any of the above + statistics-signed knots (h2.5 K3) | − | knots + | indefinite | mixed signs: hull exit |
| strain-only (no covariantization) | — | — | — | (6.1) empty: no induced gravity at all |

**T-H10 verdict: the underdetermination is NOT genuinely in the action — the
action decides, and it decides against App. E.6.** Every reading the text
permits yields at least one negative weight (B2 always; all sectors under
the favored frame reading), so the premise "w_s > 0 for every summed sector"
of Theorem VI.1 is **false on the corpus's own action**: the h2.5 conditional
DEFECT-CANDIDATE is hereby made **unconditional** (promotable on the F-R
standard: quantitative, machine-verified, no reading escape — the 8.5% wedge
h2.5 left open is closed by the half-integer lattice of the action's
reachable rays). Which way it fails is reading-dependent but the disjunction
is exhaustive: *either* all weights are negative (favored reading) and the
induced Newton constant inverts — the slaving bound survives arithmetically
but gravity comes out repulsive, contradicting VI.B — *or* the signs are
mixed and c_GW² can exit the convex hull, voiding the GW170817 "automatic
pass" (6.3) as an inequality (T4's order-of-magnitude fallback is then all
that remains). P-M3's "definite sign" prediction loses its premise in either
branch.

## 8. The minimal constitutive postulate that would settle it the other way, and its price

**P-acoustic (the unique rescue family, quantified).** Declare that on a
defected background the one-loop functional measure is the
**acoustic-covariant** one — equivalently: the flat-measure field is the
(det g₃)^w-weighted field with **w ∈ (−5/18, −2/9)** (M6a; the *entire*
rescue family — any repair whatsoever must place the effective weight in
this width-1/18 window), of which the distinguished point is **w = −1/4**,
the exact-acoustic/unimodular measure, landing every bosonic sector on
q = −3p with 𝔞₁ = +1/6 (M6b). In moduli language the same operator is
reached by dressing the inertial densities ρ₀, J ∝ (det g₃)^{−1} (M6d).
Costs, priced:

1. **The Sec. III premise.** The corpus's quantum chain (Fisher/Madelung,
   the Wallstrom discharge, the Born rates of III.D) is built in flat field
   coordinates — the flat measure is not a convenience but a premise.
   P-acoustic makes the measure defect-dependent, re-opening those [DF]
   derivations on curved backgrounds. No campaign *numeric* breaks (every
   tier is flat-background, and all readings coincide at h = 0 — Tier-1
   dispersion is blind to (p,q)); the cost is derivational, at the theory's
   foundations rather than its arithmetic.
2. **The medium ontology.** ρ₀, J ∝ 1/det g₃ says defect-inserted material
   carries no extra inertia — in tension with VI.A's own KBKK picture (a
   disclination *inserts a material wedge*) and the F5/F8 homogeneous-lattice
   bookkeeping behind the constant moduli.
3. **Unforcedness.** Nothing printed selects the window except the desired
   E.6 conclusion; the natural conventions — w = 0 (the corpus's actual
   choice) and w = −1/2 (half-density) — are both outside it (M6c). By the
   corpus's own VI.D standard ("without fine-tuning"), a width-1/18 window
   hit only by construction is tuning.
4. **What it does buy, honestly:** P-acoustic is Unruh-*aligned* — the
   exact-acoustic ray is precisely the Unruh case (h2.5), so the corpus
   could present it as "the medium's excitations gravitate acoustically."
   The honest trade is E.6 saved at the price of the Sec. III measure
   premise. And it is **not sufficient alone**: the fermionic knot sector
   still needs the separate supertrace erratum (h2.5 K3) — (p,q) cannot fix
   statistics bookkeeping. Two independent repairs remain the minimum bill
   for one appendix line.

## 9. Honest limits of this run

- Linear order in h, conformal probe class: sufficient to refute a universal
  positivity premise; not a computation of 𝔞₁ on general anisotropic
  disclination backgrounds (there the response is tensorial and h2.5's
  non-universality pathology — position-dependent ray — applies on top).
- B3's composite cone (the (β+γ)/2J + A_c/2ρ₀ mix, h2.5) was exponent-checked
  to track the pure-wryness scaling under the frame and covector readings
  (both components respond with the same exponent); the vector-reading mix
  was not fully traced — immaterial, since the vector reading is negative
  through its clean components already.
- Unsoldered readings generate additional connection terms (gradients of the
  algebra-to-space identification) beyond the two-exponent master form;
  these could move the covector reading's +1/6 by O(1) — flagged, and it
  only *strengthens* the conclusion (the sole positive cell is also the
  least stable one).
- Torsion corrections (the EC sector's own claim) enter E by contact terms
  and cannot restore a positive R-coefficient at linear order; not computed.
- The h-principle remark under reading L1′ remains at sketch grade (as in
  T4); it is not load-bearing — L1's falsity is machine-verified and L4's
  truth is elementary.

## 10. Findings routed to the campaign ledger (proposed; adjudicator's call)

- **F-H33-1 (moderate, proof-integrity; sharpens F-T4-1):** App. F.1's
  convexity is false with a *strict* counterexample (open-interval
  orientation reversal, det F_t(0) = (1−3t)(1−2t)); the lemma is ill-posed
  between the local- and global-diffeo classes; main text and F.1 cite
  different (one deep-true, one false) justifications. Erratum: Lemma II.1′
  (§4). Pairs-only theorem survives.
- **F-H33-2 (serious, argument-integrity; promotes T-H10 and completes
  F-T4-4/h2.5):** the E.6 positivity premise is false on the corpus's own
  action under every textually permitted reading — 𝔞₁^(B2) < 0
  unconditionally (p = q − 1 invariant; wedge q ∈ (2/3,5/6) unreachable on
  the action's half-integer ray lattice). Theorem VI.1 as printed fails
  either to inverted G or to hull exit; (6.3) demotes to an
  order-of-magnitude estimate. Discharge requires P-acoustic (§8, priced)
  plus the knot supertrace erratum — both new physics postulates, not
  errata.
- **F-H33-3 (minor, structural):** the strain-only reading self-destructs
  (empty Sakharov term) — worth printing, since it forces the covariant
  reading on which F-H33-2 then bites; the corpus cannot retreat from
  covariantization without losing Sec. VI entirely.

## Deliverables

- `h33_convexity.py` — Part A: symbolic twist-map determinant, segment
  determinants, strict counterexample realizability (grid + analytic
  minimum), degree-integrator control + admissible-path demonstration;
  14/14 checks.
- `h33_pq.py` — Part B: independent Gilkey/Vassilevich re-derivation of the
  master form with anchors, sector (p,q) extraction per reading, the
  p = q − 1 and q = 3/2 invariants, wedge-closure lemmas, rescue-window
  quantification; 21/21 checks.
- `h33_summary.json` — machine-readable: both verdicts, counterexample data,
  reading tables, sign map, minimal postulate with costs, full check lists.
- `h33_RESULTS.md` — this memo.
