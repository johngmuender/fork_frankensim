# L2 — Momentum time-of-flight: dBB exact + Nelson leg (F-T6-2b / F-T6-6-EXEC)

**Goal.** Execute Bricmont's Appendix-1 momentum time-of-flight computation
exactly for the 1D Gaussian ground-state packet under free evolution in
deterministic de Broglie-Bohm dynamics, then repeat the identical ToF
statistic under Nelson stochastic dynamics (GUM Sec. III's actual kinematics)
at finite diffusion, including a diffusion-constant dial. Shows that the
"measured momentum" delivered by the ToF protocol is an apparatus-active
statistic — not the instantaneous particle momentum, which is exactly zero at
t = 0 — and that at quantum equilibrium the statistic is diffusion-blind.

Within-model; nothing here bears on nature.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | dBB trajectory law: RMS of X_num(T) − X(0)√(1+T²) over 1000 trajectories ≤ 1e-6·√(1+T²) = 1.0000e-4 | RMS = 2.64e-12 (max abs 1.18e-11) | PASS |
| G2 | dBB ToF: KS of p = X(T)/T vs N(0, 1/√2), p-value > 0.1, N = 20000 | D = 0.00467, p = 0.773 | PASS |
| G3 | Nelson equivariance at t ∈ {1, 10, 100}: sup-norm of (binned − analytic ρ) over bins with ρ > 1e-3·max, within the 4/√(N_per_bin) floor | max \|diff\|/floor = 0.587 (t=1), 0.631 (t=10), 0.464 (t=100); sup-norms 1.50e-2, 2.29e-3, 2.93e-4 | PASS |
| G4 | Nelson ToF (ν = 1/2): KS of X(T)/T vs N(0, 1/√2), p > 0.1 | D = 0.00464, p = 0.781 | PASS |
| G5 | ν-dial: pairwise two-sample KS across ν ∈ {0.1, 0.5, 1.0}, all p > 0.05 | p = 0.894 (0.1 vs 0.5), 0.338 (0.1 vs 1.0), 0.384 (0.5 vs 1.0) | PASS |

## Key numbers

- T = 100, N = 20000 particles, seed 20260717; target σ_p = 1/√2 ≈ 0.70711.
- dBB ToF sample std: 0.70475; Nelson (ν = 1/2) ToF sample std: 0.71161.
  (Exact equilibrium prediction at finite T: σ = √((1+T²)/2)/T = 0.70714.)
- Instantaneous dBB velocity at t = 0: exactly 0 for every particle (Ψ real),
  measured max |v(X₀, 0)| = 0.0 — yet the ToF statistic reproduces the full
  |FT Ψ|² = π^(−1/2) e^(−p²) distribution.
- Nelson one-sample KS vs target: p = 0.988 (ν = 0.1), 0.781 (ν = 0.5),
  0.233 (ν = 1.0).
- G3 statistical floor used: absolute per-bin floor 4·√(ρ/(N·h)) =
  4ρ/√(N_per_bin) with N_per_bin = N·ρ·h (61 bins over ±4.5σ(t); 51 bins
  survive the ρ > 1e-3·max cut; smallest expected count in a kept bin ≈ 1.3).

## Method

Closed forms only (no PDE solve): Ψ(x,t) = (1+it)^(−1/2) π^(−1/4)
exp(−x²/(2(1+it))), ρ(x,t) = π^(−1/2)(1+t²)^(−1/2) exp(−x²/(1+t²)),
current velocity v = tx/(1+t²), osmotic velocity u = −2νx/(1+t²). Units
ħ = m = 1.

- **dBB leg:** dX/dt = v integrated by RK4, dt = 1e-3 for t < 10 then
  dt = 1e-2 to T = 100; compared against the analytic solution
  X(t) = X(0)√(1+t²) (G1). ToF estimator p = X(T)/T tested against
  Gaussian(0, 1/√2) by one-sample KS (G2).
- **Nelson leg:** dX = (v+u)dt + √(2ν)dW by Euler-Maruyama; dt = 1e-3 for
  t < 10, then geometric growth (×1.02 per step) capped at dt = 0.05, landing
  exactly on t = 100; X(0) ~ ρ(·,0) = N(0, 1/2); N = 20000, fixed seeds.
  Equivariance checked by binned density vs analytic ρ at t = 1, 10, 100
  (G3); ToF KS at ν = 1/2 (G4); ν-dial with independent ensembles at
  ν ∈ {0.1, 0.5, 1.0} and pairwise two-sample KS (G5).

## Connection to the GUM corpus

GUM Sec. III is Nelson kinematics with ħ = 2mν; the corpus asserts but never
exhibits that its stochastic kinematics reproduces dBB measurement
phenomenology at equilibrium. This is that exhibit (F-T6-2b / F-T6-6-EXEC):
the momentum "measurement" is apparatus-active in exactly Bell's sense — the
ground-state particle has zero instantaneous (dBB) momentum, and the Nelson
particle has zero current velocity at t = 0, yet the free-flight protocol
manufactures the Gaussian(0, 1/√2) momentum distribution |FT Ψ|² in both
dynamics. Because equivariance holds for every ν (the osmotic term rescales
with ν in lockstep with the noise), the ToF statistic is diffusion-blind at
equilibrium: ν ∈ {0.1, 0.5, 1.0} yield statistically indistinguishable ToF
distributions (all pairwise KS p > 0.33). The diffusion constant is invisible
to this class of measurement.

## Caveats

- KS p-values are seed-dependent single draws (uniform under H0); the fixed
  seed 20260717 was chosen once and not shopped.
- G3's 4/√(N_per_bin) floor treats bin counts as Gaussian; the outermost kept
  bins have expected counts ~1.3 where the Poisson tail is skewed. All
  measured ratios sit well below 1 (≤ 0.63), so the conclusion does not lean
  on the approximation.
- Euler-Maruyama is weak order 1; late-time dt = 0.05 is small relative to
  the local dynamical timescale (~t), and residual discretization bias is
  bounded above by the KS distances actually observed (D < 0.01 ≈ 1/√N
  noise floor).
- Diffusion-blindness is shown at quantum equilibrium only; out-of-equilibrium
  initial conditions would expose ν (relaxation-rate dependence) and are out
  of scope here.
- Statement about "instantaneous momentum = 0" is the dBB/Nelson current
  velocity of the real initial wavefunction — a within-model bookkeeping
  fact, not an operational claim.

## Files

- `l2_tof.py` — simulation + gates + figure.
- `L2_results.json` — all measured numbers.
- `L2_tof.png` — 4-panel figure: (a) dBB ToF vs |FT Ψ|², (b) Nelson ToF,
  (c) equivariance at t = 1/10/100 in scaled units, (d) ν-dial CDF overlay.
- `l2_samples.npz` — ToF samples and snapshots.
