# I4 — TOPOLOGY COMPLETION PACKAGE (exchange lift / curved moduli / even w)

**Auditor:** theory-audit subagent (Fable-class), 2026-07-16. Phase I4 of
ROADMAP_v6_SYNTHESIS: the three scoped opens of h4_completion.md §8 worked to
ground. **Sources:** h4_completion.md (∂ = ±2[rot], corrected group
π₁(𝒞_strat) = ℤ₂ × A/⟨2a_rot⟩, offset homomorphism ρ₊); h31_structure.md Part
B (MH, N1/N2, completions I/II/III); corpus App F (F.3 exchange lift, F.4
no-braiding, F.5 Theorem c″ + worked B = 2). **Computations:** i4_compute.py
(reuses h4_compute.py machinery). Verdicts: i4_summary.json.

---

## VERDICTS (headline)

> **A — the exchange lift (N2/F.3): FAILS as stated — obstruction exhibited.**
> In the stratified two-defect space with F.3's own proviso enforced
> (defect separations bounded below), the ordered-pair double cover exists,
> its swap monodromy σ: π₁(𝒞₂,strat) → ℤ₂ is a homomorphism, and
> σ(exchange) = −1 ≠ +1 = σ(rotation). No homotopy between them exists in
> that space; the proviso the corpus added to keep the interpolation inside
> the stratum is precisely what obstructs the lift. This is the geon
> phenomenon [IM: Friedman–Sorkin]: hard (non-mergeable) defect cores break
> the forced spin–statistics correlation. Repaired statement: χ_exch =
> χ(σ)·χ_rot with χ(σ) = ±1 a NEW free discrete choice. "ARE fermions" is
> not earned; the minimal hypothesis gains a clause (χ(σ) = +1, priced), or
> the space must admit line-merging strata — outside F.3's own proviso.
> F.4 (no braiding) SURVIVES: E² = r₁r₂ machine-verified factor by factor;
> χ(E²) = 1 for every character of the separated model.
>
> **B — curved-line moduli: ∂ = ±2[rot] SURVIVES — theorem-grade.**
> Naturality under the straight-line inclusion forces 2[rot] ∈ im ∂ on any
> moduli containing straight lines, and a Stiefel–Whitney parity theorem
> (w₂(N) = w₁(T)², which pulls back to 0 on every S²) forces
> ρ₊(∂q) ∈ 2wℤ for EVERY sphere family q of embedded lines — curved,
> dodging, or knotted-component alike. Hence [rot] (ρ₊ = −w) is never in
> im ∂: completion III is structurally impossible on any embedded-line
> moduli. The naive deformation retract to straight lines is REFUTED
> (modulo one Hatcher-input caveat, §B.3): π₂ of the asymptotic-ray data is
> larger than π₂(ℝP²); but nothing load-bearing needed the retract.
>
> **C — even w ≥ 2: χ(rot) = −1 EXISTS for all w ≠ 0, splitness-free —**
> H4 §5.5's even-w caveat is DISCHARGED and sharpened. i₊[rot] has order
> exactly 2 in π₁(𝒞_strat) for every w ≠ 0 (only ρ₊ and exactness needed),
> and an order-2 element of any abelian group admits a U(1)-character
> taking −1 (divisibility of U(1)). The dichotomy H4 saw is real but is
> about ℤ₂-VALUED characters: for even w every χ(rot) = −1 sector is
> theta-twisted (the offset loop carries a nontrivial 2w-th root of unity)
> unless the unproven FR splitting is granted. The B = 2 worked example
> survives as a sector-consistency check, NOT as a derivation: "J = 0
> forbidden" holds exactly in the (χ(σ) = +1, χ_rot = −1) sectors — 2 of 8
> characters at w = 1 — and is contingent on Part A's new free sign.

---

## PART A — THE EXCHANGE LIFT (N2 / F.3)

### A.1 The two-defect model M₂, declared

Extend H4's model M to two defects. Two parallel straight disclination
lines ℓ₁, ℓ₂ along ẑ through a₁ = (+d/2, 0), a₂ = (−d/2, 0), each with a
closed tube N_i of radius r₀ ≪ d; domain V₂ = ℝ³∖(int N₁ ∪ int N₂)
(compactified as in H4). Fiber F₂ = {U: V₂ → S³, U|∂N_i valued in U(1)_k
with meridian winding w (longitude 0) on each tube, relative degree 2}.
Base ℳ₂ = unordered pairs of disjoint straight unoriented lines. Total
space 𝒞₂,strat = the corpus's dressed two-knot sector (F.3's home).

**The based configuration** (machine check A1): U₀ = m·H₁·H₂ with
m(x) = exp(k·w·[θ₁(x) + θ₂(x)]), θ_i = arg((x−a_i)_⊥), and H₂ := H₁∘R_z(−π)
(compacton hedgehogs, disjoint supports clearing both tubes, deg 1 each).
Then U₀ is exactly R_z(π)-invariant: θ_i(R_z(−π)x) = θ_{3−i}(x) − π, so m
picks up e^{−2πkw} = 1; the hedgehogs swap and commute (disjoint supports).

### A.2 The exchange loop E and the rotation loop R₁

**E (exchange):** E(t) = R_z(πt)·U₀, t ∈ [0,1] — the rigid half-turn about
the midpoint axis. Closed (A.1's invariance). Line separation is constant
= d throughout (rigid motion): **the corpus's proviso "defect separations
bounded below along the interpolation" is realized exactly.** Meridian
winding w is preserved at every t (covariance, H4 C3; re-checked as A2).

**R₁ (rotation of one dressed knot):** the 2π rotation of knot 1 about its
own line axis, dressing included — the two-defect version of H4's [rot]:
R₁ = [2π rotation of H₁ about ℓ₁] · [offset loop β₁^{−w} on tube 1]
(the rigid one-line decomposition of H4 §5.6/C6: rigid own-axis rotation =
hedgehog frame rotation × boundary-offset drift −w). A pure fiber loop:
lines fixed, invariants (FR, off₁, off₂) = (1, −w, 0).

### A.3 The obstruction: swap monodromy

In the separations-bounded-below space, every configuration carries exactly
two disjoint lines, so the **ordered double cover** 𝒞̃₂ = {(config,
ordering of its two lines)} is an honest 2-sheeted covering (the ℤ₂ swap
action is free: the lines are distinct). Covering-space theory gives a
homomorphism σ: π₁(𝒞₂,strat, U₀) → ℤ₂ (which sheet the lift ends on),
constant on homotopy classes.

- σ(E) = −1: along E the lines are continuously labeled (positions
  e^{iπt}a₁, e^{iπt}a₂); at t = 1 line 1 occupies ℓ₂ and vice versa — the
  lift ends on the swapped sheet. (Machine check A3: continuous label
  transport, endpoint permutation = transposition.)
- σ(R₁) = +1: the lines never move.

**Therefore [E] ≠ [R₁] in π₁(𝒞₂,strat) — no exchange–rotation homotopy
exists in the space F.3's proviso defines.** The argument is basepoint- and
representative-independent: ANY exchange loop has σ = −1 by definition of
"exchange"; any rotation loop has σ = +1. The only escape is to leave the
two-disjoint-line stratum (let separations → 0, lines merge/reconnect
through a winding-2w stratum), where the double cover degenerates and σ
dies — exactly the move the corpus's parenthesis forbids. This is the known
mechanism by which spin–statistics fails for gravitational geons and for
hard-core defects [IM: Friedman–Sorkin; Sorkin's kink/particle dichotomy]:
the FR null-homotopy in the maps space necessarily passes through
configurations where the two lumps lose their separate identities (§A.6
verifies this on the undressed model), and defect lines cannot do that
while staying separated.

### A.4 Meridian/winding data tracked through the exchange

Boundary data on the moving tube 1 at loop time t, in the fixed-frame
identification (translate the tube back; meridian angle θ from x̂):
φ_t(θ) = w(θ − πt) + w·arg((a₁ − a₂)e^{iπt} + r₀e^{i(θ−πt)}) + const.
Machine checks (A2, A4, A5):
- meridian winding = w at every t (constraint never violated);
- at t = 1 the arriving data coincides with tube 2's original data exactly
  (analytically: the two arg-terms differ by π and the −wπ meridian shift
  supplies the other half of −2πw; quadrature confirms closure);
- **offset drift of E² = −w per tube** (quadrature; the meridian
  re-angling term −2πw is the whole effect — in the pulled-back frame the
  spectator puncture sits at its original position, so there is no
  relative-direction pickup). This is EXACTLY the per-tube drift of
  r₁·r₂ (own-axis rotations, H4 C6): **E² = r₁r₂ with no offset
  correction** — the undressed framed-pair relation holds factor by
  factor in the dressed model. (An earlier hand-analysis expected a
  cancellation to 0; the machine corrected it, and the corrected value
  makes the relation tighter, not looser.)
- The pure-translation braid loop B₁₂ (carry line 1 once around line 2
  without rotating) has offset drift +w per tube: [B₁₂] = β₁^w β₂^w. The
  dressing holonomy is topologically an INTEGER offset winding (loop
  closes; class possibly nonzero), machine check A6.

### A.5 F.4's claim (U(1) dressing holonomy dynamical-only; double
exchanges contractible) — machine-checked topological part

Every computable invariant of E²·(r₁r₂)⁻¹ vanishes: σ = 0; FR part = 0
(rigid 2π rotation of a degree-2 configuration; Williams parity
B mod 2 = 0 matches FR(r₁) + FR(r₂), certified via the SU(2)-lift endpoint
machinery, A7); per-tube offsets of E² equal those of r₁r₂ exactly (A5);
base loop contractible (the relative direction sweeps a great circle,
π₁(S²) = 0). In the
separated-model group of §A.6 the relation e² = r₁r₂ makes χ(E²) =
χ(r₁)χ(r₂) = (χ_F ζ^{−w})² = ζ^{−2w} = 1 for EVERY character: **no
braiding phase in any abelian sector — F.4's no-braiding content is
VINDICATED in-model.** Caveat, priced: the braid class [B₁₂] = β₁^wβ₂^w is
offset-nontrivial unless killed by the two-line connecting map im ∂₂ or by
A-blind characters — the same "A-blind" clause h31's MH already carries
(clause 3). F.4's "dynamical phase only" is exactly true for the
continuous part (integer winding ⇒ the U(1) endpoint holonomy is trivial;
what survives is at most a discrete character phase, and χ(B₁₂)² =
χ(β₁β₂)^{2w} = 1 in the corrected torsion group: at most a sign).

### A.6 The class of E in the corrected group

Separated-regime model 𝒮 (maps into 𝒞₂,strat; relations proven there hold
in π₁(𝒞₂,strat); distinctions certified by σ hold there too): fiber
invariant group ℤ₂^{FR} × ℤ_{β₁} × ℤ_{β₂}, quotiented by each line's own
direction-sweep boundary im ∂ ⊇ ⟨β₁^{2w}⟩, ⟨β₂^{2w}⟩ (H4's ∂ = ±2[rot]
applied per line; the sweep dodges the spectator line through the
contractible translation directions), extended by the exchange:

> G₂ = [ ℤ₂^{FR} × (ℤ/2w)_{β₁} × (ℤ/2w)_{β₂} ] ⋊ ⟨e⟩,
> e β₁ e⁻¹ = β₂, e r₁ e⁻¹ = r₂, e² = r₁r₂ = (FR-even, −w, −w),
> r_i = (1, −w δ_{i1}, −w δ_{i2}), σ(e) = −1.

Characters (abelianization forces χ(β₁) = χ(β₂) =: ζ, ζ^{2w} = 1):
χ(r₁) = χ(r₂) = ρ := χ_F·ζ^{−w} (χ_F = ±1), and χ(e)² = ρ² ⇒
**χ(e) = s·ρ with s = ±1 free** — the swap sign σ detects, characters
realize. Machine check A8 enumerates the w = 1 character table (8 sectors)
and verifies: χ_exch = χ_rot holds in exactly half the sectors; nothing in
the topology ties s.

**Undressed sanity check (A9):** in the undressed maps space π₁(Maps_B) =
π₄(S³) = ℤ₂ for all B [IM], both e and r₁ map to the generator — FR7
holds — while in the separated framed-pair model e·r₁⁻¹ is the nontrivial
swap class. So the classical FR null-homotopy PROVABLY exits the
separated regime; the stratified space with separations bounded below
cannot follow it. F.3's proviso is not a technicality — it is the
obstruction.

### A.7 Verdict A

**FAILS (obstruction exhibited).** F.3 as printed — "the exchange–rotation
homotopy lifts to the stratified space (defect separations bounded below
along the interpolation)" — is false in model M₂: the parenthetical
hypothesis implies the swap homomorphism σ exists, and σ separates the two
loops. χ_exch = χ_rot is NOT earned; the corrected statement is χ_exch =
χ(σ)·χ_rot, χ(σ) a new free ℤ₂ choice (geon loophole). Scoped opens: (i)
if the corpus's 𝒞_strat includes line-merging/reconnection strata
(winding-2w collision lines), σ dies there and the lift may be restorable
— that is a DIFFERENT claim than F.3's and needs a new interpolation
through the deeper stratum plus a redo of H4's ∂ on that base; (ii) the
two-line im ∂₂ beyond the per-line sweeps was not exhausted (affects only
how much of the offset factor survives, not σ). Confidence: HIGH on the
obstruction (covering-space argument + machine-checked model), MEDIUM-HIGH
that M₂ is the corpus's dressed two-knot space (same scope caveats as H4
§8, plus the two-line base).

---

## PART B — CURVED-LINE MODULI: DOES ∂ = ±2[rot] SURVIVE?

### B.1 What survives automatically: naturality

Let ℰ be any moduli of properly-embedded unknotted lines containing the
straight lines, ι: ℳ_str ≃ ℝP² ↪ ℰ, and let 𝒞^c → ℰ be the stratified
space over it (same fiber F_ℓ at the straight basepoint). The straight-line
bundle is the restriction, so the LES is natural under ι and
∂_c(ι₊ gen) = ∂_str(gen) = ±2[rot] (H4). Hence **im ∂_c ⊇ ⟨2[rot]⟩ always**
— the belt-trick mechanism does not need to be re-derived on the bigger
base. Free corollary (machine check B1): since ρ₊(2[rot]) = −2w ≠ 0,
ι₊(gen) ≠ 0 in π₂(ℰ) — the direction-sweep sphere stays essential in every
enlargement; the retraction question cannot erase it.

The only risk to H4's conclusions is therefore the CONVERSE: a new
generator q ∈ π₂(ℰ) with ∂_c(q) hitting an ODD multiple of [rot] (or any
class of odd ρ₊/w) would kill the rotation class (completion-III
resurrection). H4 §8 conjectured a "structural parity protection". §B.2
proves it.

### B.2 The parity protection theorem (the load-bearing result)

**Setup.** As in H4 §1.4, restriction-to-boundary is a fibration
ρ: 𝒞^c → 𝔅^c over ℰ, where the fiber of 𝔅^c → ℰ over ℓ is the winding-w
component of Maps(∂N(ℓ), U(1)_k) ≃ U(1) (the offset circle; canonical π₁
orientation from the fixed k-direction, so the fiberwise π₁ = ℤ local
system is trivial over each w-component). Naturality of the two LES gives
ρ₊ ∘ ∂_c = ∂_𝔅: π₂(ℰ) → π₁(U(1)) = ℤ, and for a U(1)-fibration ∂_𝔅 is the
Euler pairing: ρ₊(∂_c q) = ⟨e(L), q⟩, L = the offset circle bundle over ℰ.

**Weight-w reduction.** The offset coordinate of the boundary datum
e^{k·w·θ} transforms under a re-framing of the meridian origin θ ↦ θ + α
by offset ↦ offset + wα: the structure circle acts with weight w. Hence
L ≅ K^{⊗w}, where K is the meridian-framing circle bundle = the circle
bundle of the normal 2-plane bundle N of the line at a marked point, and
⟨e(L), q⟩ = w·⟨e(K), q⟩. (Machine check B2: the drift computations of H4
C5/C6 are exactly linear in w — re-verified at w = 1, 2, 3; and on the
straight generator ⟨e(K), gen⟩ = ±2 is H4's C2 clutching number.)

**Parity.** At the marked point, N ⊕ T = ℝ³ (trivial), T = the tautological
tangent line bundle over ℰ. Whitney: (1 + w₁N + w₂N)(1 + w₁T) = 1, so
w₁N = w₁T and w₂N = w₁N·w₁T = w₁T². For ANY sphere family q: S² → ℰ,
q*w₁T ∈ H¹(S²; ℤ₂) = 0, hence q*w₂N = (q*w₁T)² = 0, i.e. ⟨e(K), q⟩ is
EVEN. (Machine check B3: the Whitney algebra verified symbolically mod 2.)

> **Theorem (parity protection).** On any moduli ℰ of embedded lines —
> curved, translated, or a knotted component alike — ρ₊(∂_c q) ∈ 2wℤ for
> every q ∈ π₂(ℰ). Since ρ₊[rot] = −w, the equation [rot] = ∂_c(q) forces
> −w ∈ 2wℤ: impossible for w ≠ 0. **[rot] ∉ im ∂_c; the rotation class
> survives with order exactly 2 (2[rot] ∈ im ∂_c by B.1) on every
> embedded-line moduli, for every w ≠ 0.**

This upgrades H4 §8's mitigation ("any completion-III mechanism must
produce an odd multiple of [rot] ... structural, not accidental") from an
observed pattern to a proof: the obstruction is w₂ of the line's normal
bundle, which no S²-family of embedded lines can carry.

### B.3 The retraction question, answered honestly

Does ℰ deformation-retract to straight lines? **Load-bearing answer: it no
longer matters** (B.2 covers all of π₂(ℰ)). Direct answer, reconstructed:
model ℰ as smooth proper lines straight outside a compact set. Sending the
compact part to its asymptotic data fibers ℰ over the ray-pair space 𝒜
with fiber the space of long unknotted arcs with fixed ends, which is
contractible [IM: Hatcher's theorem on the unknot component of long-knot
spaces — the one external input, flagged]. So ℰ ≃ 𝒜 = {unordered pairs of
disjoint rays}/(end swap). The two end-directions can be swept
independently (dodging the codim-1 intersecting-pair wall through the
contractible translation directions), giving π₂(𝒜) ⊇ ℤ² against
π₂(ℝP²) = ℤ: **the naive retract to straight lines FAILS** — e.g. the
one-end sweep (bend an elbow, sweep one arm over S²) is a genuinely new
sphere family. MEDIUM confidence on the ℤ² lower bound (the dodge's
closure over S² was constructed only up to the standard obstruction
vanishing in contractible fibers, not machine-checked); HIGH that no
retract exists is not needed anywhere downstream. ∂_c on the new
generator: not computed individually — B.2 makes every value even in the
offset and the FR factor cannot be hit either (same argument as H4 §4.3(ii)
applies to any family of ambient bendings: their boundary loops are
isotopy-induced; recorded as a remark, not a proof — the theorem in B.2
does not need it).

**Knotted components:** knotted proper lines form different components of
the embedding space; the corpus's defect sector is the unknotted one. If
one insists: π₂ of knotted long-knot components is generically nontrivial
[IM: Budney–Hatcher-type results], so the ∂-domain grows again, and the
required computation is: evaluate ⟨e(K), q⟩ on Gramain-type generators —
but B.2's parity theorem already applies verbatim (it never used
unknottedness), so the rotation class survives there too. Stated, not
executed, per scope.

### B.4 Verdict B

**SURVIVES — theorem-grade.** ∂ = ±2[rot] on the straight generator
persists by naturality; no π₂ generator of any embedded-line moduli can
reach [rot] (parity protection, w₂-argument); order exactly 2 and
χ(rot) = −1 availability are moduli-independent. The one scoped caveat:
the U(1)-tube constraint model (H4 §8's declared choice) is still the
model in which all of this lives.

---

## PART C — EVEN w ≥ 2: DOES χ(rot) = −1 EXIST?

### C.1 The order-2 theorem for ALL w ≠ 0 (splitness-free)

H4 §8 left even w scoped out: "χ(rot) = −1 needs the FR factor split off".
Re-examining the algebra shows the caveat was about ℤ₂-VALUED characters
only. Three steps, each consuming nothing beyond H4's computed data:

1. **[rot] has infinite order in π₁(F_ℓ)^ab** for every w ≠ 0: ρ₊[rot] =
   −w (H4 C6), and ρ₊ factors through the abelianization (ℤ target).
2. **i₊[rot] ≠ 0 in π₁(𝒞_strat), and even in its abelianization:**
   ker i₊ = im ∂ = ⟨2[rot]⟩ (straight moduli; ⊆ even-ρ₊ classes on any
   moduli by B.2). [rot] = 2k[rot] would give (2k−1)[rot] = 0, contradicting
   step 1 — for EVERY w, even or odd. In G^ab the same ρ₊ argument runs
   unchanged. (H4's §5.5 parity argument needed odd w only to exclude
   [rot] ∈ 2·π₁(F)^ab + im ∂, which is the ℤ₂-character criterion.)
3. **Order exactly 2, hence χ(rot) = −1 exists:** 2[rot] ∈ im ∂ kills the
   square; an element of order exactly 2 in any abelian group admits a
   U(1)-character taking −1 on it (define χ = −1 on ⟨g⟩ ≅ ℤ₂ and extend:
   U(1) is divisible, hence injective as a ℤ-module).

> **Theorem.** For every w ≠ 0 (even included) and every embedded-line
> moduli, i₊[rot] has order exactly 2 in π₁(𝒞_strat) and U(1)-characters
> with χ(rot) = −1 exist, with no assumption on the structure of π₁(F_ℓ).

Machine check C1: in the corrected group model G = (ℤ₂ × ℤ)/⟨(0, 2w)⟩ =
ℤ₂ × ℤ/2w, exhaustive at w = 1..6: i[rot] = (1, −w) has order 2 and a
character with χ(rot) = −1 exists at every w; and in the WORST-case
non-split model (drop the ℤ₂ factor entirely: π₁(F) = ℤ_β, [rot] = −w·β...
realized as [rot] = (w/2)·(2β) only if w even — model: G = ℤ/2w with
rot = w) the order-2 statement and the −1 character still hold (C2).

### C.2 The precise structural parity protection, and its even-w price

The sharp dichotomy (machine check C2 enumerates it):

- **odd w:** χ(rot) = −1 is achievable by a ℤ₂-VALUED character with the
  offset factor blind (H4 §5.5's result); no theta-structure needed.
- **even w:** every character with χ(rot) = −1 satisfies χ_F·ζ^{−w} = −1
  with ζ = χ(β), ζ^{2w} = 1. If the FR factor is not split (χ_F absent /
  not independently addressable), −1 must come from ζ^w = −1, i.e. ζ is a
  primitive 2w-th-root-type: **the sector is theta-twisted — the offset/
  winding loop β carries a nontrivial discrete phase.** Physically: at
  even w the spinor sector exists but is necessarily dressed with a
  topological theta-angle on the U(1) winding loop; at odd w it exists
  untwisted. This is the precise content behind H4's "structural parity
  protection": the protection ([rot] never dies) is w-independent (C.1 +
  B.2); only the CHARACTER TYPE realizing −1 depends on w parity, and the
  corpus's unproven splitting is needed only to have an untwisted even-w
  fermionic sector.

### C.3 The corpus's worked B = 2 under the corrected group

Corpus F.5: "Worked B = 2: bond symmetries; J = 0 forbidden; J = 1
internal-antisymmetric ground pattern — the deuteron analog validating the
pipeline." Its topological inputs, audited against Parts A + C:

- constituents at w = ±1 (odd): χ_rot = −1 sectors exist untwisted — H4 +
  C.1 underwrite this (fine);
- 2π rotation of the B = 2 composite: FR-even (Williams parity), so
  integer J for the composite — consistent, machine check A7;
- **"J = 0 forbidden" consumes χ_exch = −1**, i.e. exchange antisymmetry
  of the two identical dressed constituents. Under Part A, χ_exch =
  χ(σ)·χ_rot = s·ρ with s free: the selection rule holds iff (s, ρ) =
  (+1, −1). Machine check C3 (character table of G₂, w = 1, 8 sectors):
  J = 0 is forbidden in exactly 2 of 8 sectors — precisely those with
  s = +1, ρ = −1; in the s = −1, ρ = −1 sectors the constituents are
  spinors whose exchange is EVEN (geon-like: spin-statistics violated) and
  J = 0 is allowed.

**Verdict on B = 2:** the worked example is internally consistent and
VALIDATES the pipeline *within the FR sector*, but under the corrected
group it is a sector-conditional statement, not a derivation: the corpus's
"J = 0 forbidden" silently consumes N2/F.3, which Part A refutes as
stated. The deuteron analog survives as: "in the (χ(σ) = +1, χ_rot = −1)
quantization sector, J = 0 is forbidden and J = 1 is the ground pattern."

### C.4 Verdict C

**EXISTS — for all w ≠ 0, splitness-free (theta-twisted at even w unless
the FR splitting is granted).** H4's even-w scope row can be closed with
the theorem of C.1; the corpus's top-cell splitting claim is now needed
only for the untwisted-sector refinement, not for QUANT. The B = 2 example
survives conditionally (C.3), inheriting Part A's new discrete choice.

---

## CONSEQUENCES, SCOPE, CONFIDENCE

**Theorem c″ / (7.6).** The three-part regrade of h31/h4 is now fully
adjudicated: (i) rotation ℤ₂ alive with χ(rot) = −1 available — COMPUTED,
and now moduli-independent and w-independent (B, C); (ii) χ_exch = χ_rot —
**REFUTED as a forced correlation in the proviso-enforced space** (A): it
is a choice, χ(σ) = +1, of the same epistemic type as the F12 vacuum
choice; (iii) the fermionic value — dynamical selection via T-B1,
unchanged. "Identical dressed knots ARE fermions" therefore needs THREE
inputs: QUANT (now theorem-grade in-model), T-B1's sector selection, and
the new discrete choice χ(σ) = +1 (or a line-merging extension of
𝒞_strat that F.3's own proviso excludes). Honest grade for c″: closure
corollary with one additional priced discrete choice — a strict demotion
of F.3, a strict promotion of F.5's N1 content.

**The ħ-closure / QUANT / T-B1 / T-H3.** Unaffected and strengthened:
QUANT consumed only N1, which Parts B + C make moduli- and w-robust. The
solvability verdict of H4 stands on strictly larger ground. Part A's
refutation costs the closure nothing — it costs the STATISTICS claim.

**F.4.** Vindicated at the topological level in-model (E² = r₁r₂ verified
factor by factor; χ(E²) = 1 in every abelian sector), with the
A-blindness caveat on single-braid offset classes matching MH clause 3.

**B = 2 / F.5 worked example.** Sector-conditional (C.3): valid pipeline
check inside (χ(σ), χ_rot) = (+1, −1); not a derivation of "J = 0
forbidden".

| item | status |
|---|---|
| N2/F.3 exchange lift, separations bounded below | **REFUTED** (swap monodromy σ; A.3) |
| N2 with line-merging strata admitted | OPEN — different claim; needs interpolation through the 2w-collision stratum + ∂ on that base |
| ∂ = ±2[rot] on curved/unknotted embedded-line moduli | **PROVEN to persist** (naturality + w₂ parity theorem; B.2) |
| deformation retract of curved moduli to straight lines | REFUTED (π₂ grows), MEDIUM confidence, not load-bearing; Hatcher input [IM] flagged |
| knotted-line components | parity theorem applies verbatim; individual ∂ values stated-as-computable, not computed |
| χ(rot) = −1 at even w | **EXISTS unconditionally** (order-2 + U(1) divisibility; C.1); untwisted version still needs the FR splitting |
| two-line im ∂₂ beyond per-line sweeps | not exhausted (affects offset factor of G₂ only, not σ, not the verdicts) |

**Confidence.** A: HIGH on the obstruction (covering-space argument is
elementary and machine-modeled; the geon precedent is standard [IM]),
MEDIUM-HIGH on model fidelity to the corpus's two-knot space. B: HIGH
(the parity theorem's only inputs are Whitney duality and H4's computed
drifts; the retract refutation MEDIUM but non-load-bearing). C: HIGH
(finite algebra, machine-enumerated; the U(1)-divisibility step is
classical). All three parts are argument-level, within-model audits, as
throughout the H/I series — nothing here is a claim about nature.

**Files:** i4_topology.md (this memo), i4_compute.py (machine checks,
reusing h4_compute.py machinery), i4_summary.json (structured verdicts).
Reproduce: `python3 i4_compute.py` (deterministic; numpy + sympy).
