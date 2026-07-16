# T2 — THEORY AUDIT: THE SURVIVING PILLARS OF THE ħ-CLOSURE

**Auditor:** theory-audit subagent (Fable-class), 2026-07-16.
**Charge:** the campaign refuted the closure's *value* (F-R4: frozen SDiff dial;
F-R5: halo-unstable saddle / saturated infimum ≥ 2√2), but its
selection/kinematics theorems replicated numerically (tier0 gate A8, ⟨r1⟩
invariants ¼ and √2, clock residual ≤ 4×10⁻⁹ even at halo saturation). This memo
audits those theorems' PROOFS AS WRITTEN — argument-level, within-model only.
Nothing here asserts anything about nature.

**Sources:** `substrate-suite/01-GUM-Omega-Paper-v2.0.1.md` (Secs. III, IV,
VII.F, VIII.D; Apps. A, E, F, G; AUD-15 §V15.1), `tier2-closure/axi3_ADJUDICATION.md`
§3 (the √2 convention tension), `tier2-closure/axi4_locked_RESULTS.md` §§3–4
(the locked threshold family), `SESSION_HANDOFF.md` §§2–3 (unit maps, exact
thresholds). Computations in
`scratchpad/t2_checks.py`, `scratchpad/t2_convention.py` (numbers reproduced
inline below; all re-runnable with numpy/scipy).

**Verdict vocabulary:** SOUND (valid given premises) / GAP (missing step, stated
precisely) / DEFECT-CANDIDATE (step fails; counterexample or derivation given) /
NOT-AUDITABLE-FROM-TEXT (premises unrecoverable; what is missing stated).

---

## 0. The closure system, formalized once

Units Λ√J (action) and Λm̃ (energy); ω in units ω₀ = m̃/√J; 𝔠 ≡ ħ/(Λ√J);
ê₀ = 64/15π = 1.35812; 𝔦₀ = 256/105π = 0.77607; ê₀/𝔦₀ = 7/4 exactly (verified).
Configuration family: BPS compacton backbone × dilation V × SDiff shape entering
only via g ≡ 𝕀/(𝕀₀V) ∈ [1, 3/2] (App. G.5). The two closure conditions:

- **(C-spin)** L ≡ 𝕀ω = j𝔠 with j the rotor number (j = ½ claimed forced);
- **(C-clock)** E_tot = 𝔠ω (de Broglie, promoted to identity by IV.H.1).

Energy functional on the family (the object every pillar below lives on):

  E(V; L) = (ê/2)(V⁻¹ + V) + L²/(2𝔦₀gV).

Envelope over V gives the exact one-parameter family
**E(L) = ê·√(1 + x)**, x ≡ L²/(ê𝔦₀g), with ω_mech = dE/dL = L/(𝔦₀gV),
V² = 1 + x. The clock condition is equivalent to E = 2L·dE/dL ⟺ **x = 1** at
j = ½. (Derived in §e below; verified numerically to 4×10⁻⁹ relative.)

The campaign's standing context, which several verdicts must carry: **F-R5
established that this family is not exhaustive of the corpus's own stated
variational problem** — the honest minimizer runs to a far-field tilt halo,
saturating at κ_paper = 1/√2. Which pillar conclusions survive outside the
family, and which are family-bound, is a recurring axis of this audit.

---

## (a) The spin-selection theorem T-B1 (Sec. IV.I)

### a.1 The argument, formalized

Premises:
- **P1** (rotor kinematics): the state is a stationary isorotating texture with
  angular momentum **L = j𝔠** (the *semiclassical* identification: L linear
  in j, not L² = j(j+1)𝔠²), E_rot = L²/(2𝕀), 𝕀 = 𝔦₀gV.
- **P2** (field energy): E_field(V) = (ê/2)(V⁻¹+V) on the BPS locus
  (E₆ ∝ V⁻¹, E₀ ∝ V — dimensional scaling under volume dilation, sound).
- **P3** (clock): E_tot = 𝔠ω.
- **P4** (quantization): j ∈ ½ℤ, j > 0 admissible ("half-integrality,
  Theorem c″").
- **P5** (positivity): w ≡ 𝔠²/(ê𝔦₀g) > 0.

Derivation chain (each step re-derived by this audit, matching AUD-15 V15.1):
stationarity ∂_V E = 0 ⟹ **V² = 1 + j²w**; the clock plus L = j𝔠 gives
ω = j𝔠/𝕀, E_tot = j𝔠²/𝕀, hence **E_rot/E_tot = j/2** (this step is purely
kinematic — it uses only P1 and P3, no potential structure) and
E_field = E_rot(2−j)/j, i.e. **1 + V² = j(2−j)w**. Subtracting:
2 = w·2j(1−j), i.e.

  **w·j(1−j) = 1, V² = 1/(1−j).**

Solvable with w > 0 iff **0 < j < 1**; the j = 1 pole (V² → ∞) and the j = 3/2
sign obstruction (w = −4/3 < 0) fall out exactly as printed. Anchors at j = ½:
w = 4 ⟹ 𝔠² = 4ê𝔦₀g (the admission normalization), V = √2, E_rot/E = ¼ — all
three recomputed ✓ (𝔠(g=1) = 2.053288, 𝔠(g=3/2) = 2.514754 = 128√42/105π,
κ = √(7/8g): 0.935414 / 0.763763). **The algebra is SOUND.** This confirms and
extends the campaign's tier0 A8 gate.

### a.2 What w is, and why it is positive (the commissioned premise question)

**w is not an integer and nothing requires it to be one.** w ≡ 𝔠²/(ê𝔦₀g) is a
*positive real* — the squared spin-to-energy-scale ratio of the branch,
dimensionless because every factor is a pure number in Λ-units. Its positivity
is inherited from: 𝔠² > 0 (ħ, Λ, J real nonzero — Λ, J are stiffness/inertia
constitutive constants, positive by the medium's stability, Sec. II); ê ≥
ê₀ > 0 (Bogomolny bound, App. B); 𝔦₀ > 0 (manifestly, Eq. 4.2 is an integral of
a square); g ∈ [1, 3/2] (Theorem IV.4, two-sided). So P5 is SOUND and cheap.
The *discreteness* input is entirely in j (P4), not in w: the theorem's
structure is "continuous solvability window (0,1) ∩ discrete spectrum {½, 1,
3/2, …} = {½}". If the task-sheet phrase "w must be a positive integer" reflects
a reading of the corpus, that reading is wrong but harmlessly so — the corpus
never claims it; integer-w plays no role anywhere in Sec. IV.

### a.3 The load-bearing unstated premise: L = jħ vs L² = j(j+1)ħ²

P1 identifies the isorotating solution's classical angular momentum with jħ.
Standard rigid-rotor quantization in the corpus's own imported literature
(Skyrme practice, refs [6–8]/[24]: Adkins–Nappi–Witten-class) uses
E_rot = j(j+1)ħ²/2𝕀. The theorem is **silently convention-dependent in its
numbers but — this audit computed it — robust in its selection**. Setting
a ≡ j (semiclassical) or a ≡ √(j(j+1)) (quantum rotor), the identical algebra
gives w·a(1−a) = 1, solvable iff 0 < a < 1, V² = 1/(1−a), E_rot/E = a/2:

| convention | j = ½ | j = 1 | j = 3/2 | anchors at j = ½ |
|---|---|---|---|---|
| L = jħ | a = 0.5, solves | a = 1, pole | a = 1.5, w < 0 | V = 1.4142, ¼, 𝔠(3/2) = 2.5148 |
| L² = j(j+1)ħ² | a = 0.8660, solves | a = 1.4142, none | a = 1.9365, none | V = 2.7321, 0.4330, 𝔠(3/2) = 3.6914 |

(`t2_convention.py`.) So: **j = ½ is the unique half-integer solution under
both conventions** — the selection conclusion is convention-robust. But the
*anchors* (V = √2, E_rot/E = ¼) and the closure *value* 𝔠₀ = 2.5147 hold only
under the semiclassical identification. The ⟨r1⟩ code's measured invariants
(E_rot/E = 0.2500 ± 0.0002, V = 1.409 ± 0.010) confirm the code implements
L = jħ — internally consistent — but the paper nowhere argues the choice, and
the corpus's own reference tradition uses the other one. **Verdict: GAP** (one
missing sentence of justification for P1; consequences confined to the value
side, which F-R4/F-R5 already own, not to the selection side).

### a.4 Scope: which exclusions survive outside the (V, g) family

The exclusions have two different logical grades, which the printed theorem
does not distinguish:

- **j ≥ 2 (semiclassical; a ≥ 2 generally): family-free.** E_rot/E = j/2 uses
  only P1 + P3; j ≥ 2 ⟹ E_field ≤ 0, impossible for any configuration with
  positive field energy. This kills those rungs in ANY enlarged space,
  including the halo-saturated branch.
- **1 ≤ j < 2 (in particular j = 1 and j = 3/2): family-bound.** Kinematically
  E_field = E(2−j)/2 > 0 is fine (at j = 3/2, E_field = E/4). The
  contradiction (2 = 0 at j = 1; w < 0 at j = 3/2) is produced by P2's
  specific dilation-virial structure V² = 1+j²w vs 1+V² = j(2−j)w. F-R5
  established that the corpus's own minimizer leaves this family (tilt halo,
  g_tot ≈ 2.9 halo-carried, saturated closure of a different functional form
  𝔠 = (4/π)G*). **The printed proof therefore does not establish "no j = 1 or
  j = 3/2 solution" for the enlarged configuration space the corpus's own
  dynamics explores.** Whether the exclusion survives there is an open
  computation (checkable: redo the two-condition system on the saturated
  branch, where E_rot/E = j/2 still holds kinematically but the E(V) profile
  is halo-modified).

**Verdict: the theorem as scoped (restricted family) is SOUND; the universally
quantified conclusion "every elementary massive knot species is spin-½" and the
armed kill clause "elementary spin-3/2 falsifies the closure" carry a GAP** —
they quantify over states the proof's family does not contain, and the corpus
itself (via the campaign's F-R5, within-model) is on record that the family is
not exhaustive. Honest grading: j = ½ uniqueness is proven-in-family and
kinematically favored globally (j ≥ 2 dead everywhere; j = 0 has no clock);
j = 1, 3/2 exclusion is in-family only.

### a.5 Circularity audit of the closure system

The commissioned question: is the closure circular anywhere? Map of the
dependency graph as printed:

1. (C-spin) L = ħ/2 cites Theorem c″ (IV.B). But c″ is a *topology* result
   (π₁ = ℤ₂, exchange = rotation sign): it makes half-odd-integer j
   **consistent**, it does not produce L = ħ/2. The j = ½ *selection* is
   T-B1 — which presupposes the closure equations. As printed this is a
   citation circle: IV.B ← c″ for what only IV.I proves, and IV.I ← IV.B's
   system.
2. IV.H.1(iii) states "Theorem c″ forces j = ½". **Misattribution — c″ forces
   nothing of the kind** (it constrains j's parity class at most; see §b.3 on
   whether it even does that). The half of IV.H.1 that matters ((i) phase ≡
   mechanical angle; (ii) Nelson-stationary ω_phase = E/𝔥_stat) does not use
   (iii).
3. **The circle is repairable, and the repair is available from the corpus's
   own parts** — this audit's proposed non-circular ordering:
   (P-i) quantization on the rotor's configuration space forces j ∈ ½ℤ
   (either parity sector); (P-ii) IV.H.1(i)+(ii) give the clock condition
   without any j input; (P-iii) T-B1's window (0,1) then selects j = ½
   *regardless of sector* (the integer candidates 0 and 1 die by no-clock and
   the pole); (P-iv) j = ½ requires the double-valued (fermionic) sector,
   which c″ certifies as consistent and whose exchange sign c″ ties to the
   rotation sign. Under this ordering nothing is circular, and c″'s actual
   content (consistency + spin-statistics correlation) is used exactly where
   it is valid.
4. Residual unrepaired premise: IV.H.1(i) — "the quantum phase is a real
   texture angle, and for the isorotating knot that angle IS the internal
   rotation" — is an ontological identification asserted from III.C, not
   derived from the field equations (no computation is offered that the wake
   phase gradient equals the isorotation angle's gradient). It is the sole
   support for C-clock being "an identity, not an imposed equation". **GAP**
   (a definitional/postulate step presented with theorem typography).

**Verdicts (a):** algebra SOUND; w-premise SOUND; L = jħ convention GAP
(selection robust, values not); universal quantification GAP (family-bound
exclusions for 1 ≤ j < 2); as-printed c″ citation structure DEFECT-CANDIDATE
(misattribution in IV.H.1(iii)/IV.B) with a valid repair exhibited; C-clock's
"identity" status GAP.

---

## (b) Theorem c″ (App. F.5 / VII.F): spin-statistics, both routes

### b.1 What is actually established vs asserted, Route 1 (fibration/obstruction)

Claimed chain: stratified one-knot space 𝒞_strat (degree-1 textures + one
disclination line ℓ) fibers over defect moduli with fiber F_ℓ = {degree-1 maps
ℝ³∖ℓ → S³ with meridian winding}; restriction to Maps(T²) is a fibration; the
π₁-contributing group sits "in the top cell": π₄(S³) = ℤ₂; meridian data add
"only abelian winding A"; hence π₁(𝒞_strat) = ℤ₂ × A with the 2π-rotation loop
generating the ℤ₂ at sign −1 ("Williams-class computation").

The *unstratified* skeleton is classical and correct (Finkelstein–Rubinstein
1968; Williams 1970: π₁(Maps₁(S³,S³)) ≅ π₄(S³) = ℤ₂ and the 2π-rotation loop is
the generator) — the corpus's "undressed corollary" is sound as imported
mathematics [IM]. The *stratified extension* — the theorem's actual novelty and
the thing the closure needs (the physical knot is dressed) — rests on four
steps none of which is more than named in the text:

- **(G-b1)** local triviality of the map 𝒞_strat → (defect moduli): asserted by
  the word "fibering"; no slice/neighborhood construction. Without it the
  homotopy exact sequence used implicitly does not exist.
- **(G-b2)** the base's own topology: the space of defect configurations
  (lines in ℝ³) has nontrivial homotopy (e.g. unoriented-line direction
  spaces carry ℤ₂'s of their own). The text folds everything non-top-cell
  into "abelian winding A" without computing the base's π₁ or the sequence
  π₁(F) → π₁(𝒞) → π₁(B); the direct-product structure ℤ₂ × A (rather than a
  possibly non-split extension) is asserted. A non-split extension would not
  change the *existence* of a ℤ₂ character but could change whether
  exchange and rotation loops are forced to the same character value.
- **(G-b3)** "the constraint lives on lower skeleta" — the obstruction-theory
  step that the meridian-winding condition cannot move the π₄(S³)
  contribution. Plausible (the meridian condition is a condition on a
  2-skeleton restriction), but the cell decomposition and the obstruction
  cocycle comparison are not given.
- **(G-b4)** "rotation-covariant meridian" — the statement that the 2π-rotation
  loop still hits the ℤ₂ generator after dressing. This is the theorem's
  crux and is one parenthesis.

App. F opens by saying appendices give "the checkable spine and the manifest
pointer" where the full computation lives in the archived TN series. **The TN
series is not in the corpus archive** (verified: `corpus/` is the WS/MIP/NR
worksheet set; no TN files). So the spine is all there is to audit.
**Verdict Route 1: GAP** — a plausible sketch whose four load-bearing steps
(G-b1…G-b4) are named, not proven, and whose full version is cited to a
document outside the auditable record. (Not DEFECT: no step is false on its
face; the unstratified core is correct mathematics.)

### b.2 Route 2 (Mayer–Vietoris)

As printed: "decompose 𝒞 = N(defect stratum) ∪ (smooth stratum); the connecting
homomorphism lands the rotation class in the same ℤ₂, sign −1."

**Mayer–Vietoris computes homology; π₁ is not a homology functor.** The
π₁-analogue is Seifert–van Kampen, which has no connecting homomorphism. There
is a legitimate reading: a sign character χ: π₁ → {±1} factors through
H₁ = π₁^ab, so an H₁ computation *can* in principle locate a ℤ₂ and evaluate
the rotation loop's class in it — but then (i) the MV decomposition of the
(infinite-dimensional, stratified) configuration space into open sets with
computable H₁(A), H₁(B), H₁(A∩B) must be exhibited, (ii) the rotation loop's
image under the maps of the sequence must be traced, and (iii) "sign −1" must
mean "the loop's H₁ class is the nonzero element", with the representation
value still a separate matter (§b.3). None of (i)–(iii) appears; the route is
one sentence. **Verdict Route 2: NOT-AUDITABLE-FROM-TEXT** — the decomposition,
the groups, and the connecting map are unrecoverable; as literally worded
("Mayer–Vietoris … connecting homomorphism … π₁-class") the sentence is not
even type-correct, so no counterexample can be constructed against it either.
What is missing, exactly: the two strata as open sets, their H₁'s, the MV
segment, and the rotation class's image.

### b.3 The forced-vs-consistent defect (the corpus's own F12 flag, applied to itself)

Suppose both routes succeed and π₁(𝒞_strat) = ℤ₂ × A with the 2π-rotation loop
generating the ℤ₂. **That still does not make the knots fermions.** Quantization
on a multiply-connected configuration space is classified by characters (more
generally unitary irreps) of π₁: the ℤ₂ factor admits *two* consistent sectors,
χ_rot = +1 (bosonic) and χ_rot = −1 (fermionic). Topology forces the
*correlation* χ_exch = χ_rot (given the exchange–rotation homotopy, F.3), never
the *value*. This is not an external standard being imposed: **the corpus
itself prints exactly this distinction one paragraph later** — App. F.6 / flag
F12, for the Hopf sector: "π₁(config) = ℤ₂ (Krusch–Speight): fermionic
quantization *consistent, not forced*", and VII.I re-flags it ("one priced
discrete vacuum choice, contrast Theorem c″"). No relevant disanalogy between
the two cases is stated or apparent: both are ℤ₂ fundamental groups generated
(in the relevant part) by the rotation loop. The phrase "maps to the generator
with sign −1" conflates two different statements — "the loop is homotopically
nontrivial" (a theorem, Williams-class) and "the loop is represented by −1"
(a sector choice).

**Verdict: DEFECT-CANDIDATE for the conclusion as worded** ("identical dressed
knots ARE fermions, spin-statistics satisfied", grade [DF]): the printed
premises support "fermionic quantization is consistent, and in the fermionic
sector spin-statistics is automatic" — precisely the F12 grade the corpus
assigns the Hopf sector, no stronger. The asymmetry between c″'s [DF-forced]
and F12's [consistent-not-forced] is an internal inconsistency of grading.

**The repair exists and is the same one as §a.5:** the *dynamics* selects the
sector — T-B1 has no solution at integer j, and a j = ½ rotor state exists only
in the χ_rot = −1 sector; therefore a substrate that hosts stationary knot
particles at all hosts them as fermions. This makes spin-statistics a
*closure* corollary (dynamical selection + topological correlation), not a
topology theorem — a demotion in provenance but a valid argument, and it is
visibly the argument the numbers verified. Note the repair inherits (a)'s
family-scope GAP through T-B1, and c″'s remaining role (the correlation
χ_exch = χ_rot) still needs F.3's lift.

### b.4 F.3's exchange–rotation lift (used by both routes)

"The exchange–rotation homotopy lifts to the stratified space (defect
separations bounded below along the interpolation)" — the parenthesis is the
proof. For the standard FR homotopy the interpolating configurations are
explicit and one can plausibly track the dressing lines, but "bounded below" is
asserted for an unspecified interpolation of line defects that must move with
their knots and could in principle be forced through each other or to infinity.
**GAP** (one construction missing: exhibit the interpolation and the bound).

---

## (c) IV.G.2 core self-averaging and IV.G.3 carrier-locked exclusion

### c.1 Theorem IV.G.1 (context, audited because G.2/G.3 stand on it)

The Bogomolny bound E₆₊₀ ≥ C₆|K| uses Λ, 𝒱_tot, deg only — sound, and
band-blind *as a bound*. The stated conclusion "ê = ê₀ in any band environment,
up to O(ε)" additionally needs *attainment* (the compacton saturates using
potential-sector structure only; kinetic modulation enters at O(ε·t)) — the
"up to O(ε)" qualifier covers this. **SOUND as scoped**, with the scope being
constitutive flag F5′ (modulation confined to the kinetic sector; Λ, m_V
uniform). F5′ is a *choice* the corpus itself flags as such and elevates by
noting F5″ would refracture ħ; the theorem is conditional on it and honestly
so.

### c.2 Theorem IV.G.2 (self-averaging)

Claim: 𝕀 = 2∫J(x)sin²f sin²θ with J(x) = J̄[1 + δ_J m(x/a)], m zero-mean
isotropic ⟹ 𝔦 = 𝔦₀[1 + O((a/R\*)²)], corrections ≤ 10⁻²⁶.

Audit: the correction term is δ_J∫m(x/a)F(x)d³x with F = 2J̄sin²f sin²θ smooth
on scale R\*. For zero-mean m this is the classic rapidly-oscillating-integral
estimate; O((a/R\*)²) is the standard two-scale homogenization order
(and for smooth periodic m it is even better). Unstated regularity premises
(m of bounded variation / finite correlation length; f smooth to the edge —
note the compacton edge is only C¹-class, which could in principle degrade the
order at the boundary layer, a measure-O(a) region contributing O(a/R\*)·(edge
weight) — at these magnitudes irrelevant). Magnitude check (this audit):
(a/R\*)² with the corpus's own scales — a from the Bragg audit's ħΩ = πħc/a ≳
10¹⁹ eV gives a ≲ 6.2×10⁻²⁶ m, and taking R\* at the electron scale,
(0.511 MeV/10¹⁹ eV)² = **2.6×10⁻²⁷ ≤ 10⁻²⁶** ✓ — the printed ceiling is
consistent with its own inputs. The corollary ("neither factor renormalizes"
⟹ one ħ across the smooth-core sector) follows from G.1 + G.2 + Theorem
IV.3's determinacy. **Verdict: SOUND as scoped** (standard-technique step with
standard unstated regularity premises; arithmetic verified). The 10¹⁹ eV input
itself is imported from elsewhere in the corpus (UV/Lorentz sector) and is not
re-derived here — flagged, not charged.

### c.3 Proposition IV.G.3 (carrier-locked exclusion)

Formalized: (i) any texture not carrier-locked self-averages (G.2) and closes
at the universal 𝔠; (ii) a carrier-locked texture has an O(1)-different
effective inertia 𝔧_n, so its closure returns a different pure number
("𝔠√𝔧_n ≠ 𝔠" — telegraphic; read: 𝔠_n = 2√(ê_n g_n 𝔧_n) ≠ 𝔠); (iii) one
vacuum has one ħ (IV.H.1), so a species with 𝔠_n ≠ 𝔠 has no stationary state
(entrainment bandwidth, IV.J); (iv) therefore no elementary species is
carrier-locked.

Two audit points:

- **The exclusion is generic, not universal.** Step (ii) needs
  𝔠_n ≠ 𝔠 *strictly*. Nothing forbids a coincidence 4ê_n g_n 𝔧_n = 𝔠²
  (indeed IV.G.4 defines admission by exactly this equation — a
  carrier-locked branch that *happened* to satisfy admission would be a
  species by the corpus's own criterion). The proposition is thus a
  genericity statement; the printed text half-acknowledges this ("GUM
  therefore **asserts** (and the published C2 null independently
  corroborates)"). **GAP**, precisely: no lemma excludes the measure-zero
  admission coincidence for locked branches, and no bound on |𝔠_n/𝔠 − 1| for
  locked states is derived anywhere in the auditable record.
- **The exclusion's teeth come from entrainment**, i.e. from IV.H.3/IV.J's
  bandwidth Δ_lock ≈ 3×10⁻³ — so (c) inherits (e)'s verdicts: the *error
  signal* is sound, the *no-stationary-state-outside-bandwidth* dynamical
  claim is the GAP-graded half (see §e).

The Bragg arithmetic in the same paragraph (lattice gaps ħΩ ∼ πħc/a ≳ 10¹⁹ eV,
thirteen orders above lepton masses) is dimensional-consistent (checked:
a ≈ 6×10⁻²⁶ m ⟹ 0.511 MeV is 13.3 orders below 10¹⁹ eV ✓) and supports the
retirement of "Floquet-gap family" language (erratum E-Ω-2) — that part is
SOUND given the imported a.

**Verdict (c): G.1, G.2 SOUND as scoped (F5′ conditional, flagged); G.3 GAP
(generic-not-universal + dependence on the entrainment half of (e)); severity
low — the corpus's own label ("asserts … corroborates") is nearly the right
grade already, and the C2 null is real corroboration within-model.**

---

## (d) IV.D radiative self-consistency: does the superradiant rejection survive the √2 tension?

### d.1 The rigid-rung algebra (verified)

Rigid closure: E(κ) = ê₀ + ½𝔦₀κ², spin 𝔠 = 2𝔦₀κ, clock ê₀ + ½𝔦₀κ² = 𝔠κ ⟹
κ² = 2ê₀/3𝔦₀ = (2/3)(7/4) = **7/6 exactly**; κ = 1.08012, 𝔠 = 2𝔦₀κ = 1.67650 —
matches the printed row (1.676, 1.080) and the "8% above the band edge" gloss
(8.012%). The other rungs recompute as printed (𝔠 = ê₀ = 1.358; 𝔠 = 2𝔦₀ =
1.552 at κ := 1; 2.0533/0.9354 at g = 1; 2.5148/0.7638 at g = 3/2). SOUND.

### d.2 The rejection inequality under every convention on the table

The campaign established (axi3 §3; axi4_locked §§3–4; exact closed forms in
SESSION_HANDOFF §2) that the corpus's κ-conventions are in √2 tension: the
closure algebra's κ is exactly ω/(√2μ) under the validated unit map (μ = tilt
channel = Yukawa-tail gap), while the corpus's over-spin control run quotes
onset at κ = 1.000 ± 0.004; and that under a direction lock the onset is not
one number but a channel family κ_crit = 1/√(2⟨sin²θ⟩_w). The full threshold
menu, against the rigid rung:

| threshold (convention/channel) | value | rigid rung 1.0801 |
|---|---|---|
| unlocked tilt/ring channel (closure-algebra convention) | 1/√2 = 0.7071 | ABOVE (+53%) |
| locked, s = sin² channel | √(7/12) = 0.7638 | ABOVE |
| locked, angle-uniform | √3/2 = 0.8660 | ABOVE |
| locked, s² = \|cos\| — **the corpus's own measured onset, 1.000 ± 0.004** | 1.0000 | ABOVE (+8%, > 3σ of the control run) |
| locked, s = \|cos\| (polar extreme) | √(5/4) = 1.1180 | below |
| translated to ω/μ (if the onset convention is taken instead) | rigid ⟹ ω/μ = √2·√(7/6) = √(7/3) = 1.5275 | ABOVE any unit threshold |

(`t2_checks.py`; threshold values are the campaign's closed forms, engine-
verified to ≤ 3×10⁻⁵ in axi4_locked.) **The rejection of the rigid rung
survives both conventions and every measured channel except one**, and that one
(the polar s = \|cos\| extreme at 1.118) is excluded as the operative onset by
the corpus's *own* control measurement (1.000 ± 0.004 sits 30σ below it).
Amusing exact identity found en route: under the ω/μ translation the rigid rung
lands at √(7/3) — numerically identical to the corpus's pair-emission threshold
ω_th = √(7/3)ω₀ (√2·√(7/6) = √(7/3) exactly). Coincidence of the algebra, not
used by either side.

**Verdict: the superradiant rejection of the rigid rung is SOUND and
convention-robust.** Also premise-checked: the rejection leans on Thm
VIII.D.2/E.3 (above-threshold √-onset radiation — standard gapped-channel
kinematics, sound as imported) and VIII.D.1/E.2 (below-cut even-analyticity ⟹
no runaways — sound as a stated property of the gapped kernel; the
even-analyticity itself is asserted at App-E spine level, minor GAP not
affecting the rejection, which needs only "above threshold radiates").

### d.3 What does NOT survive: the argument's discriminating power

The same criterion applied uniformly must judge all rungs. Under the closure
algebra's own convention (threshold 1/√2 for the tilt channel, which the lock
does not remove — axi4_locked: the ring limit survives locking, it only
requires equator-concentrated halos):

- scaling rung κ = 0.9354: ABOVE 0.7071 — over-spun;
- ⟨r1⟩ benchmark κ = 0.802: ABOVE — over-spun (measured directly by the
  campaign: halo-unstable, F-R5);
- deep-BPS endpoint κ = 0.7638: ABOVE the ring limit; sits *exactly at* the
  s = sin² locked threshold (√(7/12) — the same closed form, an exact
  marginality the corpus never mentions).

So IV.D's narrative — "self-consistency forces centrifugal deformation **until
κ < 1**: the in-gap character of matter is a derived necessity" — is
**DEFECT-CANDIDATE as a self-consistency selection principle**: the criterion
that rejects the rigid rung, applied in the convention its own algebra fixes,
also rejects every "consistent" rung in the table (this is precisely F-R5
restated at argument level; axi3 §3 proved the two corpus statements cannot
cohere in one convention, so at least one fails within-model). The table's
"consistent" column is coherent only under an unstated direction-lock plus an
unstated channel restriction (onset read in the s² = \|cos\| channel where the
control run measured 1.000). The derivation chain "spin-½ ⟹ oblateness ⟹
in-gap clock ⟹ radiative silence" (the consistency ring, IV.F) therefore has a
broken third link *as written*; it is repairable only by importing the lock as
an explicit constraint, which App. I.4's protocol does not state (F-R5's
discharge path (a)).

**Split verdict (d): rejection of the rigid rung SOUND (robust under both κ
conventions — the commissioned question, answered affirmatively); the
surrounding self-consistency/selection argument DEFECT-CANDIDATE (inherits
F-R5; asymmetric application of the radiation criterion; exact marginality of
the deep-BPS endpoint against the sin²-channel threshold left unremarked).**

---

## (e) Theorem IV.H.3: the error signal and entrainment

### e.1 The half-power law (verified exactly, with a stronger general form)

Claim: displacing 𝔥 = 𝔥\*(1+δ) along the constrained family, with the envelope
theorem dE/dL = ω: d ln ω_mech/d ln 𝔥 = +½ and d ln ω_phase/d ln 𝔥 = −½.

Audit derivation (family of §0): L = 𝔥/2 maintained under displacement;
E(L) = ê√(1+x), x = L²/(ê𝔦₀g);

- ln ω_mech = ln L − ½ln(1+x) + const ⟹ d ln ω_mech/d ln 𝔥 = 1 − x/(1+x) =
  **1/(1+x)**;
- ω_phase = E/𝔥 ⟹ d ln ω_phase/d ln 𝔥 = x/(1+x) − 1 = **−1/(1+x)**.

The closure fixed point is exactly x = 1 (E = 2L·dE/dL ⟺ 1+x = 2x; this
recovers u = g, V = √2, E = √2ê — Corollary IV.3′'s loop, all reproduced), so
the exponents are **±½ exactly at the fixed point** — the corpus's claim is
correct, and this audit adds the off-fixed-point form ±1/(1+x), which shows the
error signal is monotone and sign-faithful along the *whole* family (the two
derivatives never change sign for x > 0), a slightly stronger statement than
printed. Numerical confirmation by central differences on the exact family:
+0.500000/−0.500000 at δ = 10⁻³ (`t2_checks.py`; fixed-point identity
E = 𝔥ω verified to 4×10⁻¹⁰ relative). Premise check: the envelope theorem
needs only stationarity in the shape variables (V, g), not minimality, so this
half is untouched by F-R5's saddle finding. **Verdict: SOUND** (and
independently re-derived).

### e.2 The entrainment half

The remaining chain: beat at Δω radiates into the gapless B2/B4 channels with
an "Adler locking torque ∝ −sin Δφ" ⟹ each knot is a phase-locked
self-oscillator ⟹ the vacuum's 𝔥 is an attractor with rate λ.

- The *existence of a coupling* to gapless channels is fine (B2/B4 are gapless
  per Sec. II.E — mode count audited by the campaign's Tier 1).
- The *Adler form* and, critically, the *sign* (restoring, not anti-restoring)
  are not derived in the paper: the sign is established numerically in a toy
  class (⟨r2⟩: λ = (0.021 ± 0.004)g_χ²ω\* — "restoring in the
  thermalizing-bath class, **marginal in the integrable limit** — the claim's
  honest boundary") and at field level with one coupling input (⟨r3⟩:
  λ_field = (3.1 ± 0.7)×10⁻³ω₀). Numerics of that kind are outside this
  memo's argument-level scope; as an *argument* the restoring sign is
  premise-by-computation, and the corpus prints its own boundary honestly.
  **GAP** (precisely: no analytic sign theorem for the locking torque; the
  integrable-limit marginality means the conclusion is class-conditional).
- **F-R5 coupling (new, this audit):** "attractor" requires the displaced
  states to relax *toward the fixed point* rather than into any lower-lying
  configuration. The campaign proved the fixed-point states are halo-unstable
  saddles of the stated functional (κ_paper = 0.802 > 1/√2). Within-model,
  a self-oscillator perched on a saddle has a competing relaxation channel
  (halo shedding at rate set by the unstable mode) absent from IV.H.3's
  analysis; the entrainment conclusion is therefore conditional on the same
  unstated lock as (d). The half-power *error signal* survives regardless
  (it needs stationarity only, and the campaign measured the clock identity
  holding even at halo saturation).
- Bookkeeping: the printed τ_𝔥 ∼ 10⁻¹² s / "40 orders" in IV.H.3 was already
  corrected by the corpus's own AUD-15 (F-A15-3: 10⁻¹⁷–10⁻¹⁶ s, ≈34 orders,
  conclusion unchanged) — noted, no new charge.
- The edge-stiffening exponent λ ∝ (1−κ)^{−0.5±0.15} is [DW]-graded
  measured-scaling; nothing to audit at argument level.

**Split verdict (e): error-signal half SOUND (verified, generalized);
entrainment/attractor half GAP (sign class-conditional by the corpus's own
print; attractor premise undermined within-model by F-R5's saddle unless the
lock is imported explicitly).**

---

## Cross-cutting summary

| Pillar | Component | Verdict |
|---|---|---|
| (a) T-B1 | w·j(1−j) = 1 algebra, pole, j = 3/2 sign, anchors | SOUND (verified) |
| (a) | w a positive real (not integer); positivity | SOUND |
| (a) | L = jħ (vs j(j+1)) identification | GAP — selection robust either way (computed), anchors/value convention-pinned, choice unargued |
| (a) | "every species is spin-½" / spin-3/2 kill clause | GAP — 1 ≤ j < 2 exclusions are family-bound; F-R5 shows the family is not exhaustive; j ≥ 2 exclusion is family-free |
| (a) | closure circularity | DEFECT-CANDIDATE as printed (IV.H.1(iii)/IV.B misattribute j = ½ to c″); repairable non-circular ordering exhibited; C-clock "identity" premise remains a GAP |
| (b) c″ | unstratified FR core (π₄(S³) = ℤ₂; rotation loop nontrivial) | SOUND [IM] |
| (b) | Route 1 stratified extension | GAP (G-b1–G-b4; TN proofs absent from archive) |
| (b) | Route 2 (Mayer–Vietoris) | NOT-AUDITABLE-FROM-TEXT (MV is homology; no decomposition given) |
| (b) | "knots ARE fermions" [DF-forced] | DEFECT-CANDIDATE — topology gives consistent-not-forced (corpus's own F12 grade for the isomorphic Hopf case); repair via T-B1 dynamical sector selection exhibited, inheriting (a)'s scope GAP |
| (b) | F.3 exchange–rotation lift | GAP (bound asserted, construction missing) |
| (c) | IV.G.1 band-blind bound | SOUND as scoped (conditional on F5′, flagged) |
| (c) | IV.G.2 self-averaging ≤ 10⁻²⁶ | SOUND as scoped (arithmetic verified: 2.6×10⁻²⁷) |
| (c) | IV.G.3 carrier-locked exclusion | GAP (generic, not universal; admission coincidence unexcluded; teeth borrowed from (e)'s GAP half) |
| (d) | rigid rung κ² = 7/6, superradiant rejection | SOUND — **robust under both κ conventions** (1.0801 above every operative threshold; corpus's own onset 1.000 ± 0.004 suffices) |
| (d) | "deform until κ < 1 ⟹ in-gap derived" selection narrative | DEFECT-CANDIDATE (criterion applied asymmetrically; under the closure algebra's own convention all accepted rungs are over-spun — F-R5 at argument level; deep-BPS endpoint exactly marginal against √(7/12) channel) |
| (e) | ±½ error-signal exponents | SOUND (re-derived exactly; general form ±1/(1+x); numeric ±0.500000) |
| (e) | Adler entrainment / attractor | GAP (no analytic sign theorem; class-conditional by corpus's own boundary print; attractor premise vs F-R5 saddle) |

**The campaign's framing survives this audit with sharpened edges:** what
"replicated" is, in every case, the *kinematic/algebraic* skeleton (a-algebra,
e-error-signal, d-rigid-rejection, c-self-averaging) — and each skeleton is
genuinely sound, several verified here beyond the corpus's own statements. What
does *not* upgrade to theorem grade is exactly the set of steps that quantify
over configuration space (a-scope, b-stratified, d-selection, e-attractor),
and all four fail or gap for the same root cause the campaign's F-R5 isolated
numerically: **the closure's proofs live on a restricted family/lock that its
prose does not declare.** One pillar (b) has an additional, independent defect
class — sector-forcing — that the corpus's own F12 flag shows it knows how to
grade correctly, and doesn't, plus a repair this memo makes explicit.

**Uncertainty, printed honestly:** (i) the c″ verdicts are limited by the
absence of the TN-series proofs — Route 1's gaps may all be discharged in
documents we cannot see; the verdict is about the auditable record, not about
the (unseen) full proofs. (ii) The (a)/(d) family-scope charges assume the
campaign's F-R5 mechanics (halo algebra + unit map); those were verified at
7–9 digits and 3-D-confirmed, but if the corpus's discharge path (a) of
axi3 §5 succeeds (a lock term exhibited in the frozen ⟨r1⟩ code), the
DEFECT-CANDIDATE in (d) downgrades to a spec-omission GAP and the attractor
caveat in (e) closes; the misattribution in (a.5), the sector-forcing issue in
(b.3), and the convention-dependence in (a.3) are unaffected by any such
discharge. (iii) The quantum-rotor sensitivity analysis (a.3) uses the standard
rigid-body quantization map; if the corpus intends a genuinely classical
spinning texture with L literally = ħ/2 (a reading its own IV.H.1 supports),
the convention GAP shrinks to a definitional flag. No physical-truth claim is
made anywhere in this memo.
