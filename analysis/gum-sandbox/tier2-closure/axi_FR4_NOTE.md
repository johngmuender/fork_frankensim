# F-R4 — The Analytical Note (making the numerical finding self-contained)

## The invariance lemma (two lines)
Let T be a volume-preserving diffeomorphism (det DT ≡ 1) and Q′ = Q ∘ T the
pulled-back texture. For **any** density that is a function of the field
*value* alone — the potential density `1 − σ_P(Q)`, and the **isorotation
inertia density** `2(q₁² + q₂²)` — the change of variables y = T(x) gives

∫ Φ(Q′(x)) d³x = ∫ Φ(Q(T(x))) d³x = ∫ Φ(Q(y)) d³y.

The sextic sector obeys the same because the topological density transforms
as `b_{Q∘T} = (b_Q ∘ T)·det DT = b_Q ∘ T`, so `∫ b² d³x` is likewise
invariant. **Hence E₀, E₆, and 𝕀 are all constant on the SDiff orbit of any
base configuration.** Only the ε-suppressed kinetic sectors (E₂, E₄) vary.
(Numerically confirmed at 10⁻¹⁵ in `axi_results.json` V3_cross.)

## The corollary that stresses the centerpiece
At ε = 0, over the paper's own moduli — the SDiff orbit plus the (non-SDiff)
dilation V — the closure has **no shape dependence at all**: statics and
inertia are both frozen on the orbit, and dilation alone gives exactly the
**scaling closure**, `𝔠 = 2√(ê₀𝔦₀) = 2.0533`, `κ = √(7/8)`, `V = √2`
(reproduced to 7 digits by the Step-2 endpoint). For the paper's
shape-exact endpoint (`g = 3/2`, `𝔠₀ = 2.515`) to exist, the closure must
recruit configurations **off** the SDiff orbit; but transverse to the orbit
the static Hessian is positive (the paper's own Derrick stability), so
inertia gain carries an ε-**independent** quadratic energy cost. The
optimum therefore sits at a finite cost/gain balance — **not** at the free
geometric supremum — and the true ε → 0 endpoint is whatever that balance
yields (Step 3's run R1 measures it), not 𝔠₀ = 2.515 by construction.

Note also that G.5's own proof sketch — "∫sin²f is SDiff-invariant and
sin²θ ≤ 1, so sup g = 3/2, approached in the extreme-oblate limit" — uses
precisely the invariance that freezes 𝕀: the weight `sin²f·sin²θ` *is*
`q₁² + q₂²`, a value function, so no pullback can redistribute it. The
supremum is real as an inequality over *all* configurations; the claim that
SDiff images *approach* it is what fails.

## Relation to the standard literature
This is consistent with standard practice in the BPS-Skyrme literature the
paper imports ([8]/[24]): spinning solutions there are obtained by solving
the isorotation-deformed field equations, with the moment of inertia
computed *on the solution* — not by maximizing 𝕀 over a flat moduli space.
The paper's innovation (the free inertia dial g ∈ [1, 3/2]) is exactly the
step the lemma forbids.

## What would discharge F-R4
Either (a) Step 3's enlarged-basis solve (or the corpus's own unrestricted
⟨r1⟩ code, if produced) exhibits non-SDiff modes whose cost/gain balance
lands g* = 1.31 at ε = 0.05 and drives 𝔠(ε→0) → 2.515 — in which case G.5
has a wording defect and the numbers survive; or (b) a non-compositional
reading of "SDiff images" is produced under which G.5's derivation is
sound. Absent both, the deep-BPS endpoint reverts toward the scaling rung
and every constant downstream of 𝔠₀ (κ_phys, binding depth, ω_th, the
ε-scan anchor) shifts at the ~20% level **within the model**.
