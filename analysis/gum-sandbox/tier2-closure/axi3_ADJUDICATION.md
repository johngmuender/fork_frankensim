# Tier 2b Step 3 — Enlarged-Basis Closure: Coordinator's Adjudication

**Headline: F-R4 is HARDENED, and the run surfaces a second, sharper finding
— F-R5: the corpus's closure problem, minimized honestly over its own
enlarged configuration space, does not select 𝔠₀ = 2.515 or the ⟨r1⟩
benchmark at all. The fixed-L Routhian is variationally unstable to a
far-field direction-tilt halo whenever κ_paper > 1/√2; every solution the
corpus quotes sits above that threshold; the true infimum is
halo-saturated at κ_paper = 1/√2 exactly, with the rigorous bound
𝔠_paper ≥ 2√2 = 2.828 — which excludes the deep-BPS centerpiece
𝔠₀ = 2.5147 from below.** All numbers from `axi3_results.json` /
`axi3_run.log` (gates G0–G5 pass; the solver reproduces Step 2's solA
closure to 8×10⁻⁶ before any basis enlargement).

## 1. What was asked and what came back

Step 3 was commissioned as F-R4's adjudicator: enlarge the basis beyond
Step 2's two restricted families with direction-tilt Θ modes (the only
modes that can turn the g-dial, per the F-R4 lemma) and see whether the
ε→0 closure moves from the scaling rung (2.0533) toward the paper's
shape-exact endpoint (2.5147). The answer is neither:

| R1 (t→0) closure | 𝔠_paper | κ_paper | g* | verdict |
|---|---|---|---|---|
| hedgehog (scaling rung) | 2.0533 | 0.9354 | 1.000 | the g=1 anchor |
| b-only (tilt modes) | 2.1901 | 0.8946 | 1.093 | tilt alone: modest |
| a-only (contour) | 2.5906 | 0.8022 | 0.536 | boundary-limited (amax cap) |
| joint / extended (18 modes) | 2.775 / 2.767 | 0.772 / 0.773 | 0.62 / 0.65 | boundary-limited, still climbing |
| **saturated (well-posed limit)** | **[2.828, 3.201]** | **1/√2 exact** | g_tot ≈ 2.9 (halo-carried) | the true infimum bracket |

The direct minimizations hit their parameter caps while still descending
(`!! BOUNDARY` flags) — the family-too-small signal — and the saturated
closure identifies where they are going. At ε = 0.05 (R2) the same
happens: direct joint 𝔠 = 3.187 (+9.1σ from the benchmark 2.37±0.09),
saturated 3.381; **the ⟨r1⟩ tuple is not the minimum of its own
variational problem in the enlarged space.**

## 2. FINDING F-R5 (substantive): the tilt-halo instability and the
saturated closure

**The mechanism (derived, then verified).** A far-field tilt halo of
amplitude η changes the Routhian by dR = [∫η² dV]·(1/(8π) − κ_ours²):
the E₀ potential cost per unit ∫η² is 1/(8π) (from 1−cos η ≈ η²/2 with
the 1/4π sector constant), the rotational gain is κ_ours², and the
gradient cost per unit ∫η² can be made arbitrarily small by widening the
halo. Explicit shell probes on the full engine confirm: κ_crit measured
0.2163–0.2169 (finite-width) vs ideal 1/√(8π) = 0.19947, quadratic dI
scaling 3.999 ≈ 4. Every Step-2 solution (κ_ours 0.254–0.266) is above
threshold — as is the ⟨r1⟩ benchmark itself (κ_paper = 0.802 ↦
κ_ours = 0.226) under the compacton-exact unit map.

**The threshold is the standard over-spin criterion.** In our units the
tilt channel's temporal gap is μ_ours = √((1/8π)/1) = 1/√(8π) — identical
to κ_crit. So dR < 0 is exactly ω > μ: the classic isospinning-soliton
statement that above the meson mass the rotating solution cannot remain
localized (the corpus's own imported literature, [8]/[24], works below
this threshold or solves the deformed field equations with the shed
radiation accounted). Nothing exotic is happening; the finding is that
**the corpus's closure conditions land its solutions above the
threshold.**

**Saturation, not runaway.** Growing the halo raises I, which at fixed L
lowers κ = L/I; the minimizing sequence self-limits at κ_ours = 1/√(8π),
i.e. **κ_paper = 1/√2 exactly**, where halo addition is marginal (a flat
direction — halo_D ≈ 16–21 in the solved configurations). There the
closure reduces to 𝔠_paper = (4/π)·G*, G* = min(E_static − I/16π), with
the Bogomolny-type bound G ≥ π/√2 (verified by quadrature to 2×10⁻¹⁶)
giving the **rigorous** 𝔠_paper ≥ 2√2 = 2.82843; the family value is
G* = 2.5144 → 𝔠 = 3.2015 (R1) and 3.3812 (R2). The ¼ invariant and the
clock residual (≤ 4×10⁻⁹) hold throughout — those anchors are kinematic,
as the corpus says, and survive even at saturation.

## 3. The convention question — handled before anyone else asks it

The corpus reports (§IV.D, App. I.4) a "deliberate over-spin control run:
emission onset at κ = 1.000 ± 0.004." Our threshold sits at κ_paper =
1/√2 = 0.707. These are reconcilable only by a √2 of normalization, so I
pin it: the compacton-exact unit map (validated by reproducing the
paper's own scaling-rung tuple 2√(ê₀𝔦₀) = 2.0533, κ = √(7/8) to 7
digits, and g(0.42) to 9 digits) makes our κ_paper **exactly ω/(√2 μ)**,
where μ is the tilt-channel (= Yukawa tail) mass μ² = m̃²/(2a_ψ) that
Step 1 blind-confirmed. That is: the paper's closure-algebra κ normalizes
the clock by √2·μ. Under that convention the physical over-spin
threshold ω = μ sits at κ = 1/√2, **not** at κ = 1. The corpus's two
statements — "the closure algebra gives κ = √(7/8), 0.802, √(7/12)" and
"onset is at κ = 1" — can therefore not both hold in one convention for
the tilt channel:

- If its κ = ω/(√2μ) (the closure algebra's convention, which our exact
  map fixes), the onset quote κ = 1.000 is a √2 too high for the tilt
  channel, and every quoted solution (0.764–0.935) is over-spun — F-R5
  as stated.
- If its κ = ω/μ (onset at 1 correct), then its closure-algebra values
  are in a different convention than its onset value, and translating
  them consistently still places rung/benchmark/endpoint at
  ω/μ = 1.08–1.32 > 1 — over-spun again.

**Either way, at least one of the corpus's claims fails within-model.**
The natural resolution — and my working hypothesis — is that the ⟨r1⟩
code restricts the direction field (hedgehog-locked n̂), which (a)
removes the tilt-halo channel entirely, so its Newton–Krylov converges
to what the enlarged space reveals as a **saddle**, and (b) makes its
over-spin control run measure a different (direction-locked) channel
whose onset genuinely is κ = 1. That reading makes the corpus's numbers
internally reproducible — Step 2 reproduced g(λ) to 9 digits inside
exactly such a restricted family — but it means the ⟨r1⟩ benchmark and
the deep-BPS endpoint are artifacts of an unstated restriction, not
solutions of the stated variational problem.

## 4. What this does to F-R4 and to the corpus's centerpiece

- **F-R4 is hardened, in an unexpected direction.** The question posed
  ("do non-SDiff tilt modes drive 𝔠(ε→0) → 2.515?") is answered: tilt
  modes drive 𝔠 **past** 2.515 without stopping — through it, to a
  saturated bracket [2.828, 3.201]. The G.5 mechanism (free inertia dial
  to g = 3/2 at flat energy) fails exactly as the lemma said; what
  exists instead is an unbounded-below-until-saturation dial that no
  wording repair of G.5 can rescue. 𝔠₀ = 128√42/105π = 2.5147 is now
  **excluded from below by a rigorous bound in the enlarged space** and
  unreachable from above in every restricted family (2.05–2.28).
- **The ⟨r1⟩ benchmark tuple is a halo-unstable saddle** of its own
  functional (direct measurement: κ above threshold; shell probes give
  descent directions). Its numbers can only be recovered by adding a
  constraint (direction-locking, a stability side condition, or an
  explicit halo-exclusion) that appears nowhere in App. I.4's protocol.
- **The ε-scan sign flip (Step 2 addendum) is explained**: the corpus's
  +0.42·ε^(2/3) law is a property of the restricted (saddle) branch; on
  the honest minimum the ε-sectors are a positive penalty (exponent
  ~1.0, negative sign in both restricted families).
- **Downstream, within-model**: every constant the corpus derives
  through 𝔠₀ (ħ = 𝔠Λ√J calibration, κ_phys = √(7/12), the S4′ relation
  κ²g = 7/8's ε-run anchor, binding depths, ω_th) inherits the
  ambiguity: scaling rung −18%, saturated closure +12% to +27%, versus
  the printed 2.5147. The ħ-closure's *selection* theorems (spin-½ only,
  E_rot/E = ¼, V = √2 at the rung) are untouched — they replicated to
  10⁻⁹ throughout; it is the *value* 𝔠₀ that is now unsupported by the
  stated variational problem.

## 5. What would discharge F-R5

(a) The corpus produces the frozen ⟨r1⟩ code showing a term or constraint
in the actually-solved functional that gaps the tilt-halo channel above
the clock frequency (and amends App. I.4 to state it) — then F-R5
becomes a spec-omission finding and the benchmark survives as the
solution of the *constrained* problem, with G.5's ε→0 endpoint still
dead per F-R4; or (b) an error is exhibited in the halo cost/gain
algebra or the unit map — the places to attack are the E₀ coefficient
(1/8π per ∫η²), the inertia normalization (KE = ∫(q̇₁²+q̇₂²), fixed by
the same convention as I itself), and the compacton calibration
(validated at 7–9 digits against three independent paper quantities).
Absent both: within the model, the ħ-derivation's numerical centerpiece
rests on a restricted stationary point, its stated deep-BPS value is
excluded by its own enlarged dynamics, and the honest answers are
"scaling rung with a stability constraint" (2.053) or "saturated
closure" (≥ 2.828) — the printed 2.515 is neither.

## 6. Status and epistemics

Tier 2b Step 3 **complete**; the campaign's field-level program
(Steps 1–3) is complete. Every claim above is **within-model**: it
concerns whether the corpus's printed numbers solve the corpus's printed
equations, and says nothing about nature. The boundary-limited direct
solves are labeled as such (caps hit; saturated closure supplies the
well-posed limit); the single-solution shell probes, the N = 336/480/720
convergence ladders, and the G0–G5 gate table are the error budget. Per
the print-your-own-defects rule: my Step-3 task sheet did not anticipate
the halo channel — the agent found it, derived the threshold, and built
the saturated formulation; the adjudication frame (convention pinning,
saddle hypothesis, F-R4/F-R5 split) is mine.
