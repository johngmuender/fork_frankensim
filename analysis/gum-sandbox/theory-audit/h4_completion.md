# H4 — T-H4 STRATIFIED COMPLETION: THE DISCHARGE COMPUTATION

**Auditor:** theory-audit subagent (Fable-class), 2026-07-16. Phase H4 of
ROADMAP_v5_THEORY — the last open item of the theory-audit program.
**Charge:** h31_structure.md B.3.4 identified the determinate open problem on
which Theorem c″'s stratified extension (and, through QUANT, the solvability
of the entire ħ-closure) hangs: *evaluate the connecting homomorphism
∂: π₂(line-moduli base) → π₁(F_ℓ) on the generator of π₂(ℝP²) = ℤ, and
determine which of the three exact completions (I split / II ℤ₄ / III fatal)
is realized.* This memo performs that computation.
**Sources:** h31_structure.md (Part B, esp. B.3.2–B.3.4); h31_fr.py (certified
FR machinery, reused: SU(2)-lift certificates B3, degree engine B1/B2, finite
exactness models B5); 01-GUM-Omega-Paper-v2.0.1.md App F.5, VII.F, V.D.
**Computations:** h4_compute.py — 28/28 checks pass (numpy quadrature/phase
integration; sympy exact; finite exactness models). Verdicts: h4_summary.json.

---

## VERDICT (headline)

> **∂(generator) = ±2·[rot_ℓ]** — the connecting homomorphism sends the
> direction-sweep generator of π₂(ℝP²) to (plus or minus, an orientation
> convention) **twice the 2π-rotation loop about the line's own axis**, i.e.
> to the 4π-rotation loop. Consequently, in the declared model:
>
> - the rotation class **survives** in π₁(𝒞_strat) with order **exactly 2**;
> - **completion III (fatal: ∂ kills the rotation class) is REFUTED** — by
>   two independent routes;
> - **completion II (ℤ₄) is REFUTED** for odd meridian winding w (the
>   physical case: the elementary charge is w = ±1);
> - **completion I — the corpus's claim — holds in corrected form**: the
>   corpus's implicit ∂ = 0 is *false* (∂ is injective), but im ∂ lies in the
>   winding factor, so MH's load-bearing clause (the FR ℤ₂ injects; χ(rot) =
>   −1 characters exist) is **DISCHARGED in-model**. The printed
>   "π₁(𝒞_strat) = ℤ₂ × A" needs two corrections: A → A/⟨2a_rot⟩ (the winding
>   factor acquires torsion) and the base π₁(ℝP²) = ℤ₂ does **not** survive
>   into π₁(𝒞_strat) for w ≠ 0 (p₊ = 0: a charged disclination cannot be
>   transported around the orientation-reversing loop and return to its own
>   component).
>
> **The ħ-closure's solvability is safe**: QUANT (j ∈ ½ℤ admissible) is
> underwritten, T-B1's premise P4 stands, and T-H3's benign fixed-point
> verdict no longer consumes an undischarged topological hypothesis — *within
> the declared model and its declared scope* (straight-line moduli; U(1)-tube
> meridian constraint; N2/F.3 exchange-lift untouched). Confidence: HIGH on
> the model computation (two independent routes agree, every discrete
> invariant machine-verified); MEDIUM-HIGH that the model is the corpus's
> 𝒞_strat (scope table in §8).

An unexpected corollary with physical content (§5.4): **in the stratified
space the 4π rotation of the dressed knot is null-homotopic *through line
sweeps*** — the defect line is a belt that performs the belt trick. The ℤ₂ is
preserved not because ∂ vanishes (it does not) but because ∂'s image is
*exactly* the even part of the rotation cyclic group.

---

## 1. The model, declared (task step 1)

Everything below is proven for the following declared model **M**; §8 states
honestly what separates M from the corpus's 𝒞_strat.

**1.1 Domain and fiber.** Compactify: ℝ³ ∪ {∞} = S³; a straight disclination
line ℓ through the origin closes to a great circle ℓ̂ ∋ ∞. Its closed tubular
neighborhood N(ℓ̂) is a solid torus; the domain is the complementary solid
torus V = S³ ∖ int N(ℓ̂), ∂V = T². The fiber over ℓ is the
meridian-constrained mapping space (App F.5's object, made precise):

> F_ℓ = { U: V → S³ continuous : U|_∂V takes values in the fixed circle
> U(1)_k = {cos θ + k sin θ} ⊂ S³ with meridian winding w ≠ 0 (longitudinal
> winding 0), and relative degree 1 }.

Relative degree is well-defined because the boundary values lie in a fixed
1-complex U(1)_k ⊂ S³ (H₃(V, ∂V) = ℤ evaluated against the pulled-back
volume class); it is the corpus's "degree-1 texture". Nonemptiness at degree
1 for every w: U₀ = m_w · H with m_w(x) = exp(k·w·φ_c(x)) (φ_c = meridian
angle about the z-axis) and H a compacton hedgehog centered off-axis clearing
the tube — on ∂V, U₀ = m_w is U(1)_k-valued with winding w; rel-degree =
deg H = 1 (h31 B1's exact radial formula). Machine check **C3**.

*Why U(1)-valuedness and not exact Dirichlet data:* the constraint must be
rotation-covariant for the moduli fibration to exist (F.5's own phrase:
"rotation-covariant meridian"). Exact Dirichlet data is not
rotation-invariant (a rotation shifts the boundary phase); the U(1)-valued
winding-w condition is, and it is the standard homotopy-theoretic model of a
defect core with fixed Frank class. Covariance machine-checked (**C3**):
R_\*m_w restricted to a meridian of the rotated tube is U(1)_k-valued with
winding w for generic R ∈ SO(3), all w tested.

**1.2 Base.** ℳ = affine unoriented straight lines in ℝ³ = {(±n, b): n ∈ S²,
b ⊥ n} — the tautological quotient of the rank-2 bundle γ^⊥ over ℝP².
**Deformation retract** (task step 1's explicit argument): H_t(±n, b) =
(±n, (1−t)b) is fiberwise-linear scaling; it preserves b ⊥ n (sympy exact,
**C11**), is polynomial in t (continuous), is ℤ₂-equivariant under
(n, b) ↦ (−n, b) hence descends to the unoriented quotient, and H₁ is the
zero section. So ℳ ≃ ℝP² and π₂(ℳ) = π₂(ℝP²) ≅ π₂(S²) = ℤ (covering-space
theorem: p: S² → ℝP² induces an isomorphism on π₂ [IM]). Restricting to
lines meeting a neighborhood of the knot restricts b to an open disk
sub-bundle: same retract, same homotopy type.

**1.3 The fibration and its lifting property.** Total space
𝒞_strat = {(ℓ, U): ℓ ∈ ℳ, U ∈ F_ℓ}, projection π(ℓ, U) = ℓ. After the
retract we may take ℳ = ℝP² (lines through the origin). Two structures:

*(a) Bundle structure.* SO(3) acts by pushforward R·U = U ∘ R⁻¹ carrying
F_ℓ homeomorphically onto F_{Rℓ} (constraint covariance, C3). The map
SO(3) × F_{ℓ₀} → 𝒞_strat, (R, U) ↦ (Rℓ₀, R·U) exhibits

> 𝒞_strat ≅ SO(3) ×_{O(2)_ℓ} F_{ℓ₀}  over ℝP² = SO(3)/O(2)_ℓ,

(O(2)_ℓ = stabilizer of the unoriented z-axis line), and its pullback along
p: S² → ℝP² is the associated bundle of the principal SO(2)-bundle
SO(3) → S² (SO(2) = stabilizer of the *oriented* ẑ):

> 𝒞̃ := p\*𝒞_strat ≅ SO(3) ×_{SO(2)} F_{ℓ₀}  over S².

Well-definedness: [Rg, g⁻¹·U] ↦ (Rgẑ-line, R·U) is independent of g in the
stabilizer since gẑ = ẑ. Local triviality: local sections of SO(3) → S²
(e.g. n ↦ rotation about ẑ×n by the polar angle, on the complement of the
south pole) trivialize the associated bundle in the standard way.

*(b) Homotopy lifting.* A fiber bundle over a paracompact base is a Hurewicz
fibration (Huebsch–Hurewicz; a CW base and Serre fibration would suffice for
everything used here). This is the HLP invoked below — **not assumed** as the
corpus's G-b1 was, but *constructed* for the straight-line moduli by (a).
This is one genuine upgrade over h31's starting position: for straight lines,
Route 1's fibration hypothesis is a theorem, not a grant.

**1.4 The restriction fibration and the offset invariant.** Restriction
ρ: F_ℓ → 𝔅, U ↦ U|_∂V, where 𝔅 = {T² → U(1)_k, winding (w, 0)} — a Hurewicz
fibration since ∂V ↪ V is a cofibration (standard). Since U(1)_k = K(ℤ,1),
each component of Maps(T², U(1)) is homotopy-equivalent to U(1) (global
phase offset) × contractible, so

> π₁(𝔅, m_w) = ℤ, generated by the **offset loop** β: the constant-in-space
> phase rotation m_w ↦ e^{kθ}m_w, θ: 0 → 2π.

Hence a homomorphism ρ₊: π₁(F_ℓ) → ℤ ("offset winding"), the machine-readable
invariant used throughout; it factors through abelianization automatically
(ℤ target). This is the "evaluation/restriction fibration" of the task's
step 1, and it is also the corpus's own architecture ("restriction fibration
to Maps(T²)", F.5 Route 1).

---

## 2. The generator, pinned (task step 2)

The **tautological direction family** is F: S² → ℳ, n ↦ span(n) — sweep the
line's direction over the sphere. As a map to ℝP² this *is* the covering
projection p: S² → ℝP². Identification of its class: p₊: π₂(S²) → π₂(ℝP²) is
an isomorphism (covering-space theorem), and π₂(S²) = ℤ is generated by
[id_{S²}] (Hopf degree); F = p ∘ id, so

> [F] = p₊[id_{S²}] = the generator of π₂(ℝP², ℓ₀) = ℤ.

(Bookkeeping **C11**.) The based representative used in Route B: in polar
coordinates (ρ, φ) on D², q(ρ, φ) = line with direction
n(πρ, φ) = (sin πρ cos φ, sin πρ sin φ, cos πρ); q maps ∂D² (ρ = 1, n = −ẑ)
to the basepoint line ℓ₀ (the z-axis), and q = p ∘ (standard degree-1
collapse D²/∂D² → S²).

Because π₂(S²) → π₂(ℝP²) is an isomorphism and fibers of 𝒞̃ = p*𝒞_strat are
canonically identified with fibers of 𝒞_strat, naturality of the LES under
the pullback square gives ∂_{ℝP²}(gen) = ∂̃_{S²}(gen) **on the nose** — the
computation may be done over the oriented base S², where the bundle is the
associated bundle of SO(3) → S² (§1.3a). All π₂ of the base is exhausted by
this single generator: nothing else to evaluate ∂ on (straight-line scope).

---

## 3. Route A — the clutching computation: ∂(gen) = ±2[rot]

**3.1 Naturality reduction.** The bundle map κ: SO(3) → 𝒞̃ = SO(3) ×_{SO(2)}
F_{ℓ₀}, R ↦ [R, U₀], covers id_{S²} and restricts on fibers to the **orbit
map** ω: SO(2) → F_{ℓ₀}, R_z(θ) ↦ R_z(θ)·U₀. Naturality of the LES under
(κ, id) gives the commuting square

    π₂(S²) --∂_P--> π₁(SO(2))
       ‖                | ω₊
    π₂(S²) --∂̃--->  π₁(F_ℓ₀)

so ∂̃ = ω₊ ∘ ∂_P: the whole problem reduces to the principal bundle
SO(2) → SO(3) → S², a classical object.

**3.2 ∂_P = ±2, twice over.** *(i) Exact-sequence argument:* the segment
π₂(SO(3)) → π₂(S²) --∂_P--> π₁(SO(2)) --i₊--> π₁(SO(3)) → π₁(S²) reads
0 → ℤ → ℤ → ℤ₂ → 0. The premise that i₊ is onto with kernel 2ℤ is
machine-checked (**C1**): the SU(2) lift of the 2π fiber loop ends at −1
(nontrivial in π₁(SO(3))), the 4π lift closes (h31 B3's certificates,
reproduced). Exactness forces im ∂_P = ker i₊ = 2ℤ, and π₂(SO(3)) = 0 [IM:
π₂ of any Lie group vanishes] forces ∂_P injective; hence **∂_P(1) = ±2**.
*(ii) Clutching argument, fully numeric (C2):* with the explicit north
section s_N(n) = rotation ẑ→n about ẑ×n and south section s_S(n) =
(rotation −ẑ→n about −(ẑ×n)) ∘ R_π(ŷ), the transition τ(φ) = s_S⁻¹s_N on the
equator is verified to land in the fiber SO(2) (fixes ẑ to 10⁻¹²) and its
SO(3)-winding integrates to **−2.000000** (4001-point phase quadrature).
This is the Euler-number-2 of the unit tangent bundle of S², recomputed from
scratch rather than cited.

**3.3 Push into the fiber.** ω₊ sends the generator of π₁(SO(2)) (the 2π
z-rotation of SO(2)) to [θ ↦ R_z(θ)·U₀, θ: 0→2π] =: **[rot]** — the
2π-rotation loop of the dressed configuration about its own line axis, which
is a loop in F_{ℓ₀} because R_z preserves ℓ₀ and the fiber constraint
(covariance, C3). Therefore

> **∂(gen) = ω₊(±2) = ±2[rot]** — the 4π-rotation loop about the line axis.

Sign = orientation convention for the π₂ generator; the invariant statement
is: **im ∂ = ⟨2[rot]⟩ = the even rotation classes.** No step in 3.1–3.3 is
"it can be shown": 3.1 is LES naturality applied to an explicitly constructed
bundle map, 3.2 is machine-verified at both ends, 3.3 is the definition of
the orbit map.

---

## 4. Route B — the explicit lift (task step 3), independent of Route A

**4.1 The quaternion-field family over S² minus a disk.** Over the based
representative q(ρ, φ) of §2, define the lift

> Λ(ρ, φ) = 𝔮(ρ, φ)·U₀, 𝔮(ρ, φ) = (cos(πρ/2), sin(πρ/2)·û(φ)) ∈ SU(2),
> û(φ) = (−sin φ, cos φ, 0),

i.e. the configuration U₀ pushed forward by the rotation with quaternion
𝔮(ρ, φ) — an explicit continuous two-parameter quaternion family (it is a
hemisphere of S³). Machine checks (**C4**): continuity at both delicate ends
(ρ = 0: 𝔮 ≡ 1, φ-independent; ρ = 1: 𝔮 = (0, û(φ)), the flip circle);
projection: 𝔮(ρ, φ) rotates ẑ to n(πρ, φ) exactly — so π ∘ Λ = q, a genuine
lift of the tautological family; the line-defect boundary data stays
admissible throughout (covariance, C3). This is the "quaternion-field family
over S² minus a disk, continuous, with line-defect boundary data" the task
specifies: the disk removed is an arbitrarily small cap at n = −ẑ, where the
family closes up onto the boundary loop.

**4.2 The boundary loop.** At ρ = 1 the direction is −ẑ — the *same
unoriented line* ℓ₀ — and 𝔮(1, φ) covers R_π(û(φ)), the π-rotation about the
horizontal axis û(φ), which maps the z-axis to itself (ẑ ↦ −ẑ, checked to
10⁻¹², C4). So the boundary of the lifted disk is the loop **in the fiber**

> γ(φ) = R_π(û(φ))·U₀, φ ∈ [0, 2π] (closed: û(2π) = û(0)),

and by the definition of the connecting homomorphism, [γ] represents ∂(gen)
up to basepoint-translation along the flip path (the translation is the
standard π₁-conjugation ambiguity and does not affect any invariant used
below; Route A's based computation is the clean version, and §4.4 confirms
they agree).

**4.3 The two factors of ∂, measured on γ.**

*(i) The winding factor (explicit integration, C5).* Restricted to the tube
boundary, γ(φ)|_∂V (x) = m_w(R_π(û(φ))⁻¹x). Numerical phase quadrature
(4001-point, unwrapped, w = 1, 2, 3):
- **meridian winding of each γ(φ) = −w**: the flip reverses the meridian
  orientation in fixed coordinates — γ lives in the winding-(−w) component
  (this is the π₀ monodromy fact used in §5.3);
- **offset winding of the loop = +2w**: ρ₊-invariant of the boundary loop is
  2w times the offset generator β. Analytically the datum is
  exp(k·w·(2φ + π − φ_c)): offset 2wφ + const, two full turns per circuit —
  the machine confirms the exact analytic form.

Route A cross-check: ρ₊(∂gen) = ρ₊(±2[rot]) = ∓2w, since **ρ₊[rot] = −w**
(the 2π rotation about the line shifts every meridian phase by one full turn
times w; explicit quadrature, **C6**). Magnitudes agree exactly; the sign is
the same π₂-orientation convention as in §3.3.

*(ii) The FR ℤ₂ factor (certified SU(2) lift machinery, C8).* The boundary
loop is rotation-induced: γ = ω ∘ g with g(φ) = R_π(û(φ)) =
R_z(φ)R_π(ŷ)R_z(−φ) a loop in SO(3). By the Williams/FR machinery (h31 FR6,
certified in h31 B3), a rotation-induced loop of a degree-1 configuration
carries the FR ℤ₂ class iff the rotation loop is noncontractible in SO(3),
iff its continuous SU(2) lift is open. Machine check: the explicit lift
g̃(φ) = q_z(φ) · j · q_z(φ)⁻¹ = (−sin φ)i + (cos φ)j is continuous and
**CLOSED** (g̃(2π) = g̃(0) = j, to 10⁻¹²), and covers g (verified pointwise).
So the flip loop is null-homotopic in SO(3): **the FR ℤ₂ factor of ∂(gen) is
trivial (0)**. Consistency with Route A: ∂gen = ±2[rot] and [rot] has FR
part 1 (the 2π lift is OPEN — endpoint −1, h31 B3(ii) reproduced in C8/C1),
so 2[rot] has FR part 0; and the 4π lift closes with h31 B3(iii)'s explicit
certified null-homotopy (max deviation from S³: 3×10⁻¹⁶) — the two routes
agree factor by factor.

**4.4 Agreement.** Route A: ∂gen = ±2[rot], with (FR, offset) = (0, ∓2w).
Route B: (FR, offset) of the boundary loop = (0, +2w). Identical up to the
single orientation sign. Two genuinely different derivations — one through
the classical Euler-class-2 clutching of SO(3) → S², one by direct
construction and quadrature on an explicit quaternion family — give the same
answer. This is the "second independent route" discipline the task demands,
applied to the *positive* result (it was commissioned for a FATAL verdict;
a benign verdict earned by computation deserves the same standard).

---

## 5. The discharge: what ∂ = 2[rot] does to the exact sequence

Throughout, w ≠ 0 (the corpus's dressed knot is charged; charge IS the
disclination, V.D — the elementary case is w = ±1). Write the LES of
𝒞_strat → ℳ ≃ ℝP² at the relevant segment:

    π₂(𝒞) → π₂(ℝP²) --∂--> π₁(F_ℓ) --i₊--> π₁(𝒞_strat) --p₊--> π₁(ℝP²) --∂₀--> π₀(F_ℓ)
                ℤ                                                  ℤ₂

**5.1 [rot] ∉ im ∂: the rotation class survives (completion III refuted).**
im ∂ = ℤ·(2[rot]). If [rot] = 2n[rot] for some n, apply ρ₊:
−w = 2n(−w) ⟹ 2n = 1 in ℤ — impossible (**C7**, all w ≠ 0 tested; no
splitness of π₁(F) assumed, no knowledge of π₁(F) beyond the homomorphism
ρ₊ needed). By exactness ker i₊ = im ∂, so **i₊[rot] ≠ 0**.

**5.2 Order exactly 2.** 2[rot] ∈ im ∂ = ker i₊, so i₊[rot]² = 0: the
rotation class has order exactly 2 in π₁(𝒞_strat) — precisely what a spinor
sector needs, and now with im ∂ known rather than assumed.

**5.3 p₊ = 0 for w ≠ 0 (the base π₁ does not survive).** Transporting U₀
around the orientation-reversing generator of π₁(ℝP²) ends at a flipped
configuration, whose meridian winding in fixed coordinates is −w (C9/C5):
for w ≠ 0 this is a *different path component* of F_ℓ (meridian winding
w.r.t. fixed orientation is locally constant). So the boundary map
∂₀: π₁(ℝP²) → π₀(F_ℓ) is injective, and exactness gives im p₊ = ker ∂₀ = 0.
Hence i₊ is onto and

> **π₁(𝒞_strat) ≅ π₁(F_ℓ) / ⟨2[rot]⟩.**

(Also: 𝒞_strat stays connected — the monodromy identifies the ±w components
of the fiber at π₀ level. And exactness upstream gives im(π₂(𝒞_strat) →
π₂(ℝP²)) = ker ∂ = 0, a free corollary: the direction-sweep sphere never
lifts closed.)

**5.4 The belt-trick reading.** i₊(2[rot]) = 0 says: *the 4π rotation of the
dressed knot about its line axis is null-homotopic in the stratified space,
and the null-homotopy consists of sweeping the disclination direction over
the full sphere of lines* — the defect line performs the belt trick. The ℤ₂
survives not because the line-moduli topology decouples (∂ = 0, the corpus's
implicit reading) but because ∂'s image is exactly the even rotation
subgroup, which is what a ℤ₂ quotient needs. This is the standard escape
pattern for order parameters with defects (cf. the rigid-rotor SU(2) double
cover itself: same 2, same mechanism), and it is pleasing that the stratified
space rederives it through the Euler number of SO(3) → S².

**5.5 χ(rot) = −1 characters exist — splitness-free for odd w (C7).** A
U(1)-character of G = π₁(𝒞_strat) with χ(i₊[rot]) = −1 exists iff [rot] has
nonzero image in G^ab/2G^ab, i.e. iff [rot] ∉ 2π₁(F)^ab + ⟨2[rot]⟩. Suppose
[rot] = 2x + 2n[rot] in π₁(F)^ab; apply ρ₊: −w = 2ρ₊(x) − 2nw, i.e.
x-solvability forces w even (sympy: x = w(2n−1)/2, never an integer for odd
w). **For odd w — the elementary charge w = ±1 — the fermionic sector exists
with no assumption on the structure of π₁(F_ℓ) whatsoever**: no ℤ₂×A
splitting, no top-cell/π₄ input, only the offset homomorphism. The
disclination dressing itself protects the spinor sector: [rot] carries odd
meridian offset and can never become even. (For even w the FR ℤ₂ factor must
be split off to conclude — scope note in §8; finite character model at
w = 1: G ⊇ ℤ₂ × ℤ₂ with i₊[rot] = (1,1), exactly 2 of 4 characters give −1.)

**5.6 The corrected structure.** Granting additionally the corpus's own
structure claim π₁(F_ℓ) = ℤ₂^{FR} × A with A ∋ a_rot := ρ-projection of
[rot] (A ⊇ ℤ·β-offset; a_rot = −w in the offset coordinate):

> π₁(𝒞_strat) = (ℤ₂ × A)/⟨(0, 2a_rot)⟩ = **ℤ₂ × (A/2a_rotℤ)** —
> e.g. A = ℤ: π₁ = ℤ₂ × ℤ/2w; at w = 1: ℤ₂ × ℤ₂, i₊[rot] = (1, 1).

Corrections to F.5's printed sentence, itemized: (a) ∂ ≠ 0 — it is
*injective*; "the constraint lives on lower skeleta, adding only abelian
winding factors" is right about **im ∂ ⊆ (FR-trivial part)** and wrong as an
implicit π₂-decoupling; (b) the winding factor is torsioned by the quotient
(A → A/2a_rot: the offset winding is only defined mod 2w in the total
space); (c) π₁(ℝP²) = ℤ₂ contributes nothing for w ≠ 0 (p₊ = 0) — the
corpus's A, if it was meant to contain the base classes, loses them; (d) the
2π-rotation loop does not generate a *direct factor* on the nose (it is the
diagonal-ish (1, ā_rot)); it generates a ℤ₂ and admits χ(rot) = −1
characters, which is everything N1/QUANT needs.

---

## 6. Verdict among the three completions (task step 4)

| completion (h31 B.3.4) | status after computation | ground |
|---|---|---|
| **I. ∂ = 0, split: π₁ = ℤ₂ × A** | **holds in corrected form I′**: ∂ = injective with im ∂ = ⟨2[rot]⟩ ⊆ FR-trivial part; π₁(𝒞_strat) = ℤ₂ × A/⟨2a_rot⟩; rotation class order exactly 2; χ(rot) = −1 exists | §§3–5; C1–C10 |
| **II. ∂ = 0, non-split (ℤ₄; every ℤ₂-character gives χ(rot) = +1)** | **REFUTED for odd w** (elementary charge): χ(rot) = −1 exhibited splitness-free; also the extension question dissolves — p₊ = 0 for all w ≠ 0, so there is no base-extension to be non-split | §5.5, §5.3; C7, C9 |
| **III. ∂ onto the ℤ₂ (FATAL: [rot] dies; integer j only; closure unsolvable)** | **REFUTED — twice**: Route A (∂gen = 2[rot], and [rot] ∉ ⟨2[rot]⟩ by ρ₊-parity) and Route B (the boundary flip loop's SU(2) lift closes: FR part of ∂gen is 0, so ∂ never reaches the odd rotation class) | §§3, 4, 5.1; C2, C5–C8 |

**MH (h31 B.3.2), clause by clause.** Clause 1 (im ∂ ⊆ A / FR ℤ₂ injects):
**DISCHARGED in-model** — proven, not assumed, with the correction that
∂ ≠ 0. Clause 2 (a χ(rot) = −1 character on the resulting group):
**DISCHARGED in-model**, splitness-free for odd w. Clause 3 (N2: the dressed
exchange loop / F.3's interpolation lift): **NOT TOUCHED** by this
computation — it remains the GAP h31/T2 recorded (it lives in the two-knot
space, a different base). The tidy sufficient version of MH ("dressing is
π₁-inert") is *false as stated* (the dressing does act: it torsions A and
kills the base π₁) but its load-bearing consequence (the rotation ℤ₂
survives with fermionic characters available) is true and now computed.

**Honesty item — what the corpus wrote vs what is true.** F.5 Route 1's
*conclusion* (rotation ℤ₂ alive at sign −1, winding factors abelian) is
**vindicated at every point the closure consumes**; its *proof sketch* was
not only unproven (h31's finding) but its implicit mechanism (∂ = 0,
π₂-decoupling, "top cell exactly as unstratified") is **false in the
model** — the true mechanism is ∂ = 2[rot] with even image. The corpus was
right for a wrong-in-detail reason. Route 2 (Mayer–Vietoris) remains as
h31 B.4 graded it: type-incorrect, cannot be repaired into independence;
"verified by two independent routes" stays unsupported *as printed* — though
this memo now supplies, ironically, the two honest routes.

---

## 7. Consequences (task step 5) — measured

**T-B1 / spin selection.** T-B1's premise P4 (j ∈ ½ℤ admissible, half-odd
sector available) consumed QUANT, and QUANT consumed MH (h31's escalation:
"the stakes are not statistics but solvability"). MH clauses 1–2 are now
discharged in-model: **the half-odd-integer sector exists for the dressed,
charged (odd-w) knot; the closure system w·j(1−j) = 1 retains its admissible
spectrum; the ħ-closure fixed point (j, w, V) = (½, 4, √2) stands on
topology that has now been computed rather than hypothesized** — within the
declared model. The fatal scenario (empty admissible spectrum, no stationary
knot at any j) is off the table for straight-line moduli.

**T-H3 (closure cycle).** Part A's benign fixed-point verdict consumed QUANT
as an external premise (h31 cross-cutting table: "T-H4 is load-bearing for
the closure's solvability"). That load is now carried: the benign verdict's
last consumed topological hypothesis is discharged in-model. The T-H3 × T-H4
coupling row of h31's table can be closed accordingly.

**Theorem c″ / (7.6).** The regrade path is unchanged and now better
grounded: c″'s honest content = (i) rotation ℤ₂ alive in 𝒞_strat with
χ(rot) = −1 available [now COMPUTED in-model, was GAP/MH]; (ii) χ_exch =
χ_rot correlation [needs N2/F.3's lift — still GAP]; (iii) the fermionic
*value* [dynamical selection via T-B1 — h31's repair, unchanged; "ARE
fermions" remains overreach independent of this memo]. The printed
π₁(𝒞_strat) = ℤ₂ × A should be corrected to ℤ₂ × A/⟨2a_rot⟩ with the p₊ = 0
remark; the "two independent routes" sentence should be replaced by the
actual two routes (§§3–4) or deleted.

**What this does NOT do.** It does not make the dressed knot a fermion (value
still dynamically selected); it does not close N2; it does not extend to
curved/moving line moduli (below); it says nothing about nature — this is
argument-level, within-model auditing, as throughout the H-series.

---

## 8. Scope, residual opens, and confidence

| item | status |
|---|---|
| straight-line moduli (≃ ℝP², after retract C11) | computation COMPLETE: ∂ evaluated on the full π₂ (one generator) |
| curved/shaped line moduli (Emb(ℝ, ℝ³)-type; the TN-series' claimed home) | OPEN: π₂ is larger; additional generators' ∂ not computed here. Mitigation: every extra generator maps into the fiber through loops of ambient isotopies; a repeat of §3's naturality on the larger base is well-posed, finite work |
| U(1)-tube meridian constraint as the model of F.5's "meridian-constrained mapping space" | DECLARED model choice; rotation-covariant (C3), matches F.5's own words; other homotopy-equivalent core models would not change π₁-level conclusions |
| even w ≥ 2 | χ(rot) = −1 needs the FR factor split off (corpus's top-cell claim, unproven); odd w (elementary charge) needs nothing |
| w = 0 (uncharged dressed knot) | outside this memo's w ≠ 0 hypotheses (π₀ argument and ρ₊-parity both die); but w = 0 is the undressed sector, already SOUND classically (h31 B.1–B.2) |
| N2 (exchange–rotation lift in the dressed two-knot space, F.3) | UNTOUCHED — remains the GAP of record |
| π₁(F_ℓ) full structure (the ℤ₂ × A claim itself) | not needed for §§5.1–5.5 (that is the point of the ρ₊-parity argument); consumed only in §5.6's corrected-form statement |

**Confidence: HIGH** that ∂(gen) = ±2[rot] in model M (two independent
derivations, one classical-topology, one constructive-numerical; every
discrete invariant machine-verified; 28/28). **HIGH** on the downstream
algebra (§5: exactness bookkeeping, finite models, sympy parity).
**MEDIUM-HIGH** that model M is the corpus's 𝒞_strat: the two declared gaps
are line-shape moduli and the constraint-model choice; neither has an
identified mechanism to resurrect completion III (any completion-III
mechanism must produce an *odd* multiple of [rot] in im ∂, which the
ρ₊-parity obstruction blocks for every rotation-induced sphere family — a
structural, not accidental, protection).

**Files:** `h4_completion.md` (this memo), `h4_compute.py` (28/28 checks),
`h4_summary.json` (structured verdicts). Reproduce: `python3 h4_compute.py`
(deterministic; numpy + sympy; ~5 s).
