# Tier 4 — Coordinator's Adjudication (4A: 3-D unrestricted / 4B: certified
diagnostics / 4C: direction-locked closure)

**Headline: F-R5 is CONFIRMED in unrestricted 3-D with the mechanism
quantitatively fingerprinted, and the direction-locking working hypothesis
of the Step-3 adjudication is REFUTED as a rescue — but 4C's channel
analysis produced two exact identities that, for the first time, explain
*where the corpus's numbers come from*: its measured over-spin onset
(1.000 ± 0.004) equals the polar-halo channel threshold (exactly 1), and
its deep-BPS κ₀ = √(7/12) equals the sin²-weighted halo threshold
(exactly √(7/12), a rational-moment identity 6/7). The corpus's numbers
are real numbers — they are thresholds and saturation values of the halo
family, not minima of the stated closure.**

Sources: `field3d_RESULTS.md`/`field3d_results.json` (4A),
`fs-cosserat-pilot/RESULTS.md` (4B, golden root `f88731af…`),
`../tier2-closure/axi4_locked_RESULTS.md`/`axi4_locked_results.json` (4C).

## 1. 4A — the assumption-free 3-D verdict on F-R5

The fixed-L Routhian descent on a full 3-D Cartesian grid, no symmetry
imposed, from the perturbed static solution at the Step-2 spin:

- R fell **monotonically 0.560 below** the axisymmetric stationary value
  (same-grid comparison), crossing it at iteration ≈110.
- κ fell 0.374 → through solA's 0.2656 → **through the threshold
  1/√(8π)** → 0.182 at the iteration cap, still falling; inertia rose
  through the soft-halo saturation √(8π)·L; halo fraction rising and
  **box-limited** (the 1.3× box descends deeper at matched iterations —
  the runaway is physical, not discretization).
- **The mechanism fingerprint**: 88% of the inertia gain is priced in E₀
  at exactly the threshold-defining rate 1/(16π); the surplus (including
  the sub-threshold κ overshoot) is the bounded tilt dial of Step 3
  riding on top. This is the halo channel, not lattice mush.
- **Two-sided**: the L = 4.0 control (below threshold) formed no
  far-field condensate; its inertia self-limited at the tilt ceiling
  (3/2)·I_hedgehog — the F-R4 supremum appearing in its lawful role as a
  *bound*, not a free dial.
- **The clock forces the over-spun regime**: L_clock (bisection,
  E_rot/E = 0.250000 exactly) gives κ = 1.245× the halo threshold on the
  3-D solution (1.331× on solA, 1.474× on the clean axisymmetric
  control). The instability is not an edge case of the closure — it is
  where the clock condition lives.

Standing caveats (all documented in the report): near-BPS lattice
pathologies (charge unwinding; quadrature descent through the Bogomolny
floor) required a corner-scheme objective, a degree anchor and a one-sided
BPS wall (measured footprint ≤ 2×10⁻³ at every endpoint); descents ended
at iteration caps — the verdicts are monotone trends with the crossing
points resolved, not converged endpoints.

## 2. 4C — the working hypothesis dies, and something better replaces it

Step 3's adjudication hypothesized that the corpus's frozen ⟨r1⟩ code
direction-locks the hedgehog, which would remove the tilt-halo channel
and make the benchmark a genuine constrained minimum. Measured outcome:

- **The lock does not close the halo family.** A direction-locked
  profile halo F = η·s(θ) has threshold κ_crit_paper = 1/√(2⟨sin²θ⟩_w)
  (w = s²): √3/2 = 0.866 for angle-uniform halos (my prediction —
  confirmed), but **sliding down to the unrestricted 1/√2 as the halo
  concentrates at the equator**. The benchmark (κ = 0.802) is stable
  against uniform halos and **unstable against ring halos**; the fully
  converged locked closure runs away accordingly: at ε = 0.05 the
  interior-converged tuple is 𝔠 = 3.26 (+9.9σ), κ = 0.725, g* = 0.40
  (−23σ), V = 7.1 (+570σ). **Direction-locking is not the corpus's
  missing constraint.** Whatever stabilized ⟨r1⟩ must also suppress
  equatorial profile halos — a radial-only far field, a compact
  computational domain, or plain Newton-to-a-saddle convergence.
- The locked ε-scan still contradicts the corpus's deficit law (fit
  −0.447·ε^0.418: wrong sign again; cap-limited rows flagged).

**The two exact identities (the run's permanent contribution):**

1. **The corpus's over-spin onset is a halo threshold.** For a
  polar-weighted channel (s² ∝ |cosθ|), ⟨sin²θ⟩_w = 1/2 *exactly*, so
  κ_crit_paper = **1.00000 exactly** — against the corpus's measured
  "emission onset at κ = 1.000 ± 0.004." Verified numerically at
  0.98810–0.99089 for the smoothed variant (finite-width correction).
2. **The deep-BPS κ₀ is a halo threshold.** For s = sin²,
  ⟨sin⁶⟩/⟨sin⁴⟩ = 6/7 *exactly*, so κ_crit_paper = **√(7/12) exactly** —
  the corpus's κ₀ = √(7/12) = 0.7638, which it derives from closure
  algebra (κ²g = 7/8 at g = 3/2).
3. (Looser, numerical:) the locked **uniform-channel saturated closure at
  the hedgehog core** gives 𝔠_paper = **2.37096** at κ = √3/2 — against
  the benchmark 𝔠 = 2.37 ± 0.09.

Identities 1–2 are rational-moment facts of the sphere and could in
principle be coincidences of sin²-moment algebra (the paper's κ₀ also
comes from sin²-moments); identity 3 is a one-part-in-300 numerical match
of a solve nobody tuned. A fourth proximity, recorded without a claim:
the saturated objective G* = min(E_static − I/16π) converges downward
with basis size to 2.51442–2.51446 (4C's independent engine cross-checks
axi3's value at 6×10⁻⁵) — within ~1.2×10⁻⁴ of 𝔠₀ = 2.514735 itself,
though measurably *below* it and still decreasing; whether the halo-
saturated G* and the paper's 𝔠₀ coincide exactly in some limit is open. Taken together they support a coherent reading:
**the ⟨r1⟩ pipeline's reported numbers are threshold/saturation values of
the halo family its own solver was blind to** — the solver sat at (or
drifted along) the edge of an instability it could not represent, its
"onset" diagnostic measured the one channel it could represent (polar),
and its benchmark landed at the uniform-channel saturation point. This
reading is a hypothesis, flagged as such; discharging it requires the
frozen ⟨r1⟩ code.

## 3. 4B — the certified rung (context)

13/13 certified claims PASS (golden root `f88731af…`, bit-identical
replay): degree with measured O(h²) convergence, sector energies checked
two ways against fs-ivl enclosures, the F-R4 SDiff-invariance and the ¼
invariant as certified interval claims. The diagnostics that every solver
above relies on are now in frankensim's evidence discipline, ready for
the full fs-cosserat port when the solver spec freezes.

## 4. Consolidated within-model status of the ħ-closure chain after Tier 4

- **Kinematics (replicates, certified):** spin selection, E_rot/E = ¼,
  V = √2 at the rung, the scaling rung 𝔠 = 2√(ê₀𝔦₀) = 2.0533, the clock
  algebra, the four-branch spectrum. Untouched by every finding.
- **The stated variational closure (fails, three independent ways):**
  F-R4 (SDiff dial frozen — proof + machine precision), F-R5 (halo
  instability above κ_paper = 1/√2 — 2-D enlarged basis, now 3-D
  assumption-free, two-sided, clock-forced), 4C (direction-locking does
  not rescue it).
- **The corpus's printed numbers (𝔠₀ = 2.515, κ₀ = 0.764, benchmark
  2.37/0.802/1.31/0.42, onset 1.000):** not minima of the stated problem;
  consistent with thresholds and saturation points of the halo family
  (two exact identities + one 0.3% match), i.e. artifacts of a solver
  restriction nowhere stated in App. I.4.
- **What the corpus must produce to discharge F-R4/F-R5:** the frozen
  ⟨r1⟩ code, with the stabilizing restriction or term identified — and an
  amended App. G.5/I.4 stating it. Until then, within the model, the
  ħ-derivation's value 𝔠₀ has no derivation from the printed equations.

## 5. Epistemics

Everything above is within-model. No statement here bears on nature, on
GUM's physical viability, or on any empirical stake (S1's neutrino window
is untouched — its pipeline replicated cleanly in Tier 0). The campaign's
findings are about the relationship between the corpus's printed
equations and its printed numbers, established with the corpus's own
tools: derived thresholds, converged solves, certified arithmetic, and
pre-stated discharge paths.
