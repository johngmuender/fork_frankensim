# M4 — Nelson non-equilibrium relaxation: τ(ν) (F-T7-M4)

**Goal.** First measurement of whether GUM's *actual* kinematics — Nelson
diffusion at ν = ħ/2m — relaxes to Born equilibrium faster than the
deterministic de Broglie–Bohm (dBB) guidance used in every rates anchor
(⟨V-SIM-4⟩, tier3-born). Every III.D relaxation-rate statement in the corpus
was calibrated on dBB; at equilibrium the statistics are ν-blind (tier6 L2),
but out of equilibrium the osmotic term u = ν∇ln|ψ|² drives ρ toward |ψ|²
directly. Pre-registered direction: τ(ν) strictly decreasing.

**Campaign:** Tier 7 program extensions, workstream M4 (ROADMAP_v9_PROGRAM.md).
**Date:** 2026-07-18. **Code:** `m4_nelson_relax.py` (fully seeded; resumable
per-run cache in `cache/`; total compute ≈ 56 min across staged invocations,
log `m4_run.log`). **Raw output:** `m4_results.json`, `m4_hdata.csv`.
**Figure:** `m4_relax.png`. **Units:** ħ = m = 1 (so GUM's ν = ħ/2m = 0.5).

**Within-model; nothing here bears on nature.**

## 0. Verdict in one line

Nelson noise **accelerates** Born relaxation at every ν tested, monotonically
in ν for both mode sets: at GUM's own ν = 0.5 the relaxation time is
**17× (M = 9) / 41× (M = 16) shorter than dBB** — so the corpus's dBB-calibrated
III.D rates are direction-signed conservative bounds — but the small-ν approach
to the dBB limit is *logarithmically slow* (τ ≈ a + b·ln(1/ν), dBB rate
recovered only near ν ~ 10⁻⁸–10⁻⁹), so the pre-registered continuity threshold
τ(0.1)/τ(0) ≥ 0.7 **fails honestly** (measured 0.13 / 0.07).

## 1. Setup (tier3-born conventions, reused)

Everything tier3: 2-D box [0, π]², ψ(t) = M^(−1/2) Σ e^(iθ_mn) φ_mn e^(−iE_mn t)
(equal weights, phase seed 42, exact mode-sum evolution), M = 9 and 16
(first k×k modes; M = 4 skipped per spec); ρ₀ = |φ₁₁|² non-equilibrium start
(rejection-sampled, seed 12345 + M); coarse-grained
H̄(t) = Σ P̄ ln(P̄/Q̄) at 32×32 (headline) and 16×16, Q̄ from box-averaging |ψ|²
on a 256² fine grid; exponential fit window up to H̄ ≤ 0.1·H̄(0) or t = 2π.
tier3's exact ρ₀ and phase conventions were fully recoverable from
`born.py`/REPORT.md — no substitutions needed.

New here: **dX = (v + u)dt + √(2ν)dW**, v = Im(∇ψ/ψ), u = ν∇ln|ψ|²;
ν ∈ {0 (dBB control), 0.1, 0.5, 1.0} + small-ν diagnostics {0.02, 0.05};
N = 10,000 (roadmap spec; tier3 used 20,000 — noise floor 32²: 0.0512 vs
0.0256). ν = 0 uses tier3's RK4 + substep integrator verbatim (dt = 2×10⁻³);
ν > 0 uses Euler–Maruyama (dt = 10⁻³) with the same per-particle drift-substep
guard and displacement cap, noise √(2νh) per substep. **Walls:** reflecting
(positions folded back into the box each substep); the osmotic term
~ 2ν·cot(x) already repels where |ψ|² → 0. No-flux verified numerically:
all particles in-box at all times; reflection rate ≤ 1.26 events/particle
over t = 4π at ν = 1; equilibrium-control wall-band occupancy (band π/32)
matches the Born value within Poisson error (7 of 8 configs; worst
0.0065 vs 0.0046, ~2.8σ single-snapshot; the flat-H̄ control excludes any
coarse-grained wall flux). Two deviations from tier3, both stated: (i) the
fit window additionally ends where H̄ ≤ 2× the finite-N floor — inactive for
every ν = 0 run, so the G1 comparison is convention-identical to tier3;
(ii) output cadence 0.05/0.02 (instead of 0.1) for ν > 0, resolution only.

**Error model (stated honestly):** 200-resample particle bootstrap; the
bootstrap τ distribution's *location* is biased upward (resampling doubles
the multinomial noise, inflating the effective floor — worst for the fastest
decays), so quoted errors are the bootstrap *relative* spread applied to the
point estimate; cross-checked against 3-seed replicate scatter where
available (bootstrap is the conservative i.e. larger of the two).

## 2. Gates

| id | spec | measured | verdict |
|---|---|---|---|
| G1 | dBB control (ν = 0) within 25% of tier3's τ (10.2 / 4.9 at M = 9/16) | τ = 10.66 (dev 4.5%), 5.33 (dev 9.9%); fit r² 0.90/0.98; comparison caveat: N = 10k vs tier3's 20k (floor 2×), same integrator/binning/window | PASS |
| G2 | equivariance control at every ν: equilibrium-born ensemble stays at floor | 8/8 controls (both M × 4 ν) flat over full t = 4π: max H̄/floor = 1.12–1.18 (tier3's own control peaked at 1.18); ν-blindness of equilibrium confirmed at every ν | PASS |
| G3 | τ(ν) strictly decreasing in ν for both M (pre-registered) | M9: 10.66 → 1.391 → 0.635 → 0.412; M16: 5.33 → 0.395 → 0.130 → 0.085; strict also with ν = 0.02, 0.05 inserted; every step decisive except the last (0.5 → 1.0): paired-bootstrap P(τ(1.0) < τ(0.5)) = 0.74 (M9) / 0.87 (M16), direction consistent at ~1σ | PASS |
| G4 | τ(0.5)/τ(dBB) both M with bootstrap errors; continuity τ(0.1)/τ(0) ≥ 0.7 | ratios measured: **0.060 ± 0.005 (M9), 0.024 ± 0.004 (M16)**; continuity clause: τ(0.1)/τ(0) = **0.131 ± 0.011 / 0.074 ± 0.005 — both ≪ 0.7 → clause FAILS**; diagnosis in §4: singular (logarithmic) advective-mixing limit, not a numerical defect (dt-converged, controls flat, seed-replicated) | PARTIAL (quantification PASS; pre-registered continuity criterion FAIL, diagnosed) |
| G5 | consequence paragraph: which corpus claims move which way | printed in §5 (report-only) | PASS |
| dt | halving changes τ by < 5% | worst single-seed 4.6%; three configs exceeded 5% single-seed (8.3–11.1%) and resolved to 1.5–3.1% under 3-seed averaging at both dt — scatter, not bias | PASS |

## 3. Key numbers (32² grain; errors as per §1 error model)

| ν | τ (M = 9) | fit r² | τ (M = 16) | fit r² |
|---|---|---|---|---|
| 0 (dBB) | 10.66 ± 0.26 | 0.90 | 5.334 ± 0.113 | 0.98 |
| 0.02 † | 2.662 ± 0.042 | 0.96 | 1.122 ± 0.102 | 0.96 |
| 0.05 † | 1.825 ± 0.101 | 0.95 | 0.786 ± 0.073 | 0.91 |
| 0.1 | 1.391 ± 0.116 | 0.91 | 0.395 ± 0.028 | 0.94 |
| 0.5 (= ħ/2m) | 0.635 ± 0.048 | 0.79 | 0.130 ± 0.022 | 0.93 |
| 1.0 | 0.412 ± 0.017 | 0.75 | 0.085 ± 0.015 | 0.91 |

† diagnostic points (added to map the small-ν approach after the continuity
clause failed; not part of the pre-registered grid).

16² grain, same ordering (grain-robustness): M9 τ = 7.94 / 2.21 / 1.49 /
1.03 / 0.46 / 0.28; M16 τ = 4.15 / 1.06 / 0.68 / 0.385 / 0.119 / 0.074.

Headline ratios: **τ(0.5)/τ(dBB) = 0.060 ± 0.005 (M9), 0.024 ± 0.004 (M16)**
— a 17×/41× acceleration at GUM's own diffusion constant.

## 4. The failed continuity clause, diagnosed

τ(0.1)/τ(0) = 0.13/0.07 ≪ 0.7. The small-ν diagnostic points show the
mechanism: over ν ∈ [0.02, 1.0] the measured law is close to
**τ ≈ 0.28 + 0.56·ln(1/ν)** (M9, r² = 0.97) and **τ ≈ 0.26·ln(1/ν)**
(M16, r² = 0.92) — the classic high-Péclet advective-mixing form. The dBB
flow stirs f = ρ/|ψ|² into exponentially finer filaments; any finite ν cuts
the cascade at the Batchelor scale, so the coarse-grained relaxation rate
depends on ν only logarithmically. Extrapolating the log law to τ = τ_dBB
puts the crossover near **ν ~ 10⁻⁸ (M9) / 10⁻⁹ (M16)**: the ν → 0 limit is
continuous in principle but singular in practice — no ν on any laboratory-like
scale behaves like dBB out of equilibrium. This is physics of the model, not
numerics: dt-halving converged (< 5% seed-averaged), the equivariance
controls are flat at every ν, and the effect is 20× the quoted errors.
Pre-registration is pre-registration: the clause is graded FAIL as written.

## 5. Consequences for the corpus (G5 — the honest double edge)

**Strengthened — laboratory/equilibrium quality (III.D).** Every corpus
statement of the form "residual non-equilibrium relaxes on system dynamical
times" was calibrated with dBB mixing (⟨V-SIM-4⟩, tier3: τ ≈ 10.2/4.9).
Under GUM's actual Nelson kinematics at ν = ħ/2m the same observable relaxes
17–41× faster in this box. The III.D [CAL] rates therefore become
direction-signed *conservative bounds*: whatever equilibrium quality the
corpus claimed with dBB rates, the model's own kinematics does strictly
better at every ν tested, and equivariance (the padlock premise) stays exact
at every ν (G2). This is the first ν-dependence measurement of the Born
H-theorem rates under GUM's stated kinematics.

**Weakened — relic non-equilibrium survival windows.** The same factor cuts
the other way: any corpus arithmetic in which primordial or relic quantum
non-equilibrium survives long enough to be observable inherited dBB
survival times. Under Nelson kinematics the survival window shrinks by the
same ~17–41× (box-toy scale), and §4's log law says even ν orders of
magnitude below ħ/2m still relaxes several-fold faster than dBB. Relic
non-equilibrium signatures are correspondingly harder to keep alive in this
model class. Both edges must be carried together: the corpus may not quote
the acceleration for equilibrium quality and the dBB rates for relic windows.

## 6. Method notes

Exact ψ(t) by mode sum (no PDE error); drift evaluated analytically from the
separable bilinear form (tier3's O(N·k) contraction, extended with
u = 2ν·Re(∇ψ/ψ)); per-particle substepping keyed to |v + u| with tier3's
displacement cap (guards node and wall spikes); reflection applied after
every substep with out-of-box events counted; H̄ and fits exactly tier3's
estimator code paths; bootstrap resamples whole trajectories (multinomial
over particles) and refits within the same window convention; dt-halving
run for all 8 pre-registered configs (escalated to 3 seeds per dt where the
single-seed comparison exceeded 5%). All seeds fixed and printed in
`m4_results.json`; every run cached (`cache/`) so the pipeline is resumable
and every number reproducible bit-for-bit.

## 7. Caveats (standing)

1. N = 10,000 (spec) → 32² floor 0.0512: the fast-config fits reach only
   ~2× floor; residual floor contamination biases fast τ *upward* — i.e.
   the quoted accelerations are, if anything, understated (conservative for
   the headline direction).
2. Bootstrap location bias documented in §1; errors are the relative spread
   (conservative vs 3-seed scatter). The M16 ν ∈ {0.5, 1.0} τ's carry
   17–54% relative bootstrap spreads; their difference is only ~1σ (G3
   caveat).
3. G1 compares N = 10k to tier3's N = 20k (floors differ 2×); deviation
   4.5–9.9% is comfortably inside the 25% gate but is not a same-N rerun.
4. Single phase seed (42) and single ρ₀ convention (|φ₁₁|²), inherited from
   tier3; τ(ν) values are for this box/mode family — the *monotonicity* and
   the log-law shape are the transportable content, not the absolute τ's.
5. Euler–Maruyama is weak-order-1; the dt gate bounds the residual bias at
   < 5% on τ. Wall reflection + drift cap is a standard but not
   boundary-exact SDE scheme; the wall-band check bounds the residual at
   the Poisson level of one config (~2.8σ single snapshot).
6. The ν ~ 10⁻⁸–10⁻⁹ crossover in §4 is a log-law extrapolation two-plus
   decades beyond the measured ν range — shape statement, not a measured
   point.

## 8. Files

`m4_nelson_relax.py` (script, self-contained), `m4_results.json` (every
measured number incl. gates), `m4_hdata.csv` (all H̄(t) curves, long format),
`m4_relax.png` (H̄ decay both M + τ(ν)), `m4_run.log` (run log),
`cache/` (per-run resumable cache). Filed as **F-T7-M4**.
