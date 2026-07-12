# fs-gum-sde — Phase B4 RESULTS (Nelson / stochastic-QM sector)

**Status: COMPLETE — 21/21 gates pass.** Spec: `gap_survey_v3.json` dimension
"Nelson/stochastic-QM"; port source `tier3-born/born.py`; Tier-3 referee
`tier3-born/REPORT.md`. Standalone crate (empty `[workspace]`, path deps
`fs-rand`, `fs-math` only), `cargo build --release` warning-free, 13/13 unit
tests. Gate/demo binary: `cargo run --release` (wall time ≈ 309 s,
single-threaded — vs born.py's 829 s for the deterministic pilot).
Outputs: `hdata.csv`, `hdata_fine.csv`, `results.json` (all deterministic —
no wall-clock content).

## Gate table

| id | gate | value | target | pass |
|----|------|-------|--------|------|
| G1.a | EM/SRA1 share ΔW bits at every (traj, slot, comp) | true | true | PASS |
| G1.b | single-trajectory isolation replay (bitwise) | true | true | PASS |
| G1.c | ensemble order-independence (reverse sweep, bitwise) | true | true | PASS |
| G1.d | OU ensemble in-process replay (bitwise) | true | true | PASS |
| G2.a | EM OU \|mean err\| at T=1 (exact 0.367879, N=10⁵) | 1.88e-3 | < 5e-3 | PASS |
| G2.b | EM OU \|var rel err\| (exact 0.108083) | 2.33e-3 | < 2e-2 | PASS |
| G2.c | SRA1 OU \|mean err\| | 2.11e-4 | < 5e-3 | PASS |
| G2.d | SRA1 OU \|var rel err\| | 1.03e-2 | < 2e-2 | PASS |
| G2.e | strong order, EM on OU (additive ⇒ 1.0) | 1.000 | [0.85, 1.15] | PASS |
| G2.f | strong order, SRA1 on OU (linear ⇒ superconverges to 2.0) | 2.000 | [1.85, 2.15] | PASS |
| G2.g | strong order, EM on cubic well (additive ⇒ 1.0) | 1.032 | [0.85, 1.15] | PASS |
| G2.h | strong order, SRA1 on cubic well (additive ⇒ 1.5) | 1.609 | [1.35, 1.65] | PASS |
| G2.i | strong order, EM on GBM (multiplicative ⇒ 0.5) | 0.491 | [0.35, 0.65] | PASS |
| G2.j | weak order, EM on OU (⇒ 1) | 1.119 | [0.8, 1.3] | PASS |
| G2.k | weak order, SRA1 on OU (⇒ 2) | 2.182 | [1.6, 2.5] | PASS |
| G3.a | stationary φ₁₁: max H̄₃₂ over t ∈ [0, 2π] | 0.0275 | < 2× floor (0.0512) | PASS |
| G3.b | stationary φ₁₁: mean H̄₃₂ | 0.0245 | < 1.5× floor (0.0384) | PASS |
| G4.a | M=16 cg32 floor-subtracted exponential-fit r² | 0.9778 | ≥ 0.95 | PASS |
| G4.b | M=16 H̄₃₂(4π) | 0.0261 | < 2× floor (0.0512) | PASS |
| G4.c | control max H̄₃₂ (equilibrium stays at floor) | 0.0295 | < 2× floor (0.0512) | PASS |
| G4.d | Nelson τ (cg32, born.py protocol) | 0.132 | (0, 6.07] | PASS |

Convergence-order measurement protocol: dt-halving ladders on FIXED Wiener
paths — the finest-level (ΔW, ΔZ) pairs are generated once through the
WienerStream addressing and coarser levels aggregate them exactly (ΔW by
summation, ΔZ by the composition rule ΔZ_coarse = Σⱼ[ΔZⱼ + h_f·(W(tⱼ)−W(t₀))]).
Strong error = E|X_k(T) − X_{k+1}(T)| (2000 coupled paths, dt = 2⁻⁴…2⁻¹⁰,
reference ladder to 2⁻¹²); weak error = |E[X_k(T) − X_{k+1}(T)]| with common
random numbers (20 000 paths). OU: θ=1, σ=0.5; cubic well: dX = (X−X³)dt +
0.5 dW; GBM: μ=σ=1 against the exact pathwise solution.

**SRA1 order fine print** (measured, documented in `src/sde.rs`): the generic
additive-noise strong order 1.5 shows up on the NONLINEAR cubic well (1.609
measured, the h^1.5-limiting residuals carry f''); on LINEAR drift (OU) those
residuals vanish identically and SRA1 superconverges to strong order 2.0
(2.000 measured). Milstein was NOT implemented: for additive noise its
correction term ½σ(∂σ/∂x)((ΔW)²−h) vanishes identically and the scheme
degenerates to Euler–Maruyama exactly — SRA1 is the genuinely higher-order
scheme for this sector (survey dimension 6, missing-piece 6).

## The Nelson H-theorem demo (G4)

Setup: N = 20 000 trajectories, M = 16 modes (kmax = 4, born.py's exact
phases — numpy seed 42 embedded as round-trip decimal constants), ρ₀ =
|φ₁₁|² by per-trajectory rejection sampling (streams seeded 12345),
Nelson SDE dX = (v + u)dt + √(ħ/m) dW with v = Im ∇ψ/ψ, u = Re ∇ψ/ψ
(ħ = m = 1, ν = ħ/2m = ½), EM per substep with born.py's guard structure
(|b|h ≤ 0.05, ≤ 32 substeps, drift-displacement cap), reflecting walls,
dt = 2×10⁻³, T = 4π (the exact ψ-revival period). H̄ = Σ P̄ ln(P̄/Q̄) at
32×32 and 16×16 with empty-cell skip; noise floors cg²/2N = 0.0256 / 0.0064.

### H̄(t) — Nelson M=16 and equilibrium-born control

| t | H̄₃₂ Nelson | H̄₁₆ Nelson | H̄₃₂ control | H̄₁₆ control |
|------|--------|--------|--------|--------|
| 0.0000 | 0.5645 | 0.4866 | 0.0238 | 0.0062 |
| 0.1000 | 0.1552 | 0.1328 | 0.0261 | 0.0067 |
| 0.2000 | 0.0778 | 0.0575 | 0.0257 | 0.0076 |
| 0.3000 | 0.0494 | 0.0297 | 0.0261 | 0.0061 |
| 0.4000 | 0.0371 | 0.0172 | 0.0268 | 0.0064 |
| 0.5000 | 0.0328 | 0.0148 | 0.0249 | 0.0061 |
| 1.0000 | 0.0298 | 0.0100 | 0.0282 | 0.0067 |
| 2.0000 | 0.0269 | 0.0076 | 0.0274 | 0.0070 |
| 3.1000 | 0.0272 | 0.0067 | 0.0266 | 0.0074 |
| 6.3000 | 0.0262 | 0.0067 | 0.0269 | 0.0071 |
| 9.4000 | 0.0257 | 0.0070 | 0.0275 | 0.0070 |
| 12.5000 | 0.0248 | 0.0060 | 0.0271 | 0.0064 |
| 12.5664 (=4π) | 0.0261 | 0.0054 | 0.0270 | 0.0062 |

(Full curves: `hdata.csv` every 0.1; `hdata_fine.csv` every 0.01 for
t ≤ 1.5.) Both Tier-3 qualitative behaviors reproduce: H̄ decays
near-exponentially to the finite-N noise floor and STAYS there through the
ψ revival at t = 4π (no revival of the particle distribution: H̄₃₂(4π) =
0.0261 ≈ floor 0.0256, vs Bohm's 0.0457); the equilibrium-born control is
pinned at the floor for the whole run (mean 0.0262, max 0.0295 — matching
Tier-3's control mean 0.0261, max 0.0303). The stationary-state gate (G3,
ψ = φ₁₁, where v = 0 and u = (cot x, cot y)) holds the floor for 2π
(mean 0.0245) — the sharpest test of the σ² = ħ/m vs u = (ħ/2m)∇ln ρ
factor bookkeeping: a factor-2 error in either fails it immediately.

### Nelson vs de Broglie–Bohm relaxation — the new measurement

| quantity | de Broglie–Bohm (Tier-3, deterministic) | Nelson (this crate) |
|----------|------------------------|---------------------|
| τ (cg32, born.py fit protocol) | 4.853 (r² 0.987) | **0.132** (r² 0.962) |
| τ (cg32, floor-subtracted fit) | — | 0.0958 (r² 0.978) |
| τ (cg16, born.py protocol) | 3.996 | 0.107 (r² 0.983) |
| H̄₃₂ at t = 4π | 0.0457 | 0.0261 (= floor) |
| ratio τ_Nelson/τ_Bohm (cg32) | — | **0.0272 (≈ 37× faster)** |

No campaign artifact had measured this pair on the same system before.
Nelson dynamics relaxes ρ → |ψ|² roughly **37× faster** than deterministic
de Broglie–Bohm guidance on the identical wavefunction: the osmotic drift +
diffusion (ν = ½) erases coarse-grained structure at rate ~ν k² for mode
number k, so the |φ₁₁|²-vs-|ψ|² disequilibrium (structure up to k ≈ 4) dies
on t ~ 0.1, while Bohm relaxation must wait for chaotic advection to fold
trajectories (t ~ 5). This satisfies the survey's expectation τ_Nelson ≤
τ_Bohm = 4.85 with two orders of magnitude to spare.

Protocol caveats for the τ comparison: (1) the Nelson decay is far too fast
for born.py's 0.1 diagnostic cadence, so the fits use a fine early-time
series (every 0.01, t ≤ 1.5; `hdata_fine.csv`) — the coarse cadence would
leave ~4 points in the fit window; (2) the headline 0.132 uses born.py's
raw-log fit protocol for parity with the Tier-3 τ = 4.853, but the Nelson
curve REACHES the noise floor inside its fit window (Bohm's never gets
close within t ≤ 2π), which bends the raw log fit — the floor-subtracted
fit (τ = 0.0958, r² = 0.978, the G4.a gate) is the cleaner exponentiality
test; (3) samplers differ (fs-rand Philox vs numpy PCG), so H̄(0) = 0.5645
vs Tier-3's 0.5572 — resampling-level difference, consistent with the ~5%
fit tolerance Tier-3 quotes.

## Determinism statement

Every draw is a pure function of (seed, kernel, trajectory, slot, component)
through fs-rand Philox4x32-10 counter streams; all transcendentals go
through fs-math `det` (cross-ISA bitwise contracts); all reductions are
fixed-order sequential sums; no threads, no platform libm in any kernel.
Verified in-run: single-trajectory isolation replay and reverse-order
ensemble sweep are bitwise identical (G1.b/c), an ensemble re-run
in-process is bitwise identical (G1.d). Verified across runs: two
consecutive full executions of the gate binary produce byte-identical
`hdata.csv`, `hdata_fine.csv`, and `results.json`; golden fingerprint
(FNV-1a 64 over all H̄ curves): **cca041e8558ce556**. The WienerStream
draw-index contract (4 reserved draws per (slot, component): Z1 for ΔW,
Z2 for the SRA1 ΔZ partner, reserved whether consumed or not) is documented
in `src/wiener.rs` and makes ensembles order-independent and every
trajectory replayable in isolation by construction.

## Deviations from spec (with reasons)

1. **SRA1 strong-order gate split in two.** Measured SRA1 on OU is 2.0, not
   1.5 — linear-drift superconvergence (the h^1.5-limiting Itô–Taylor
   residuals carry f'' and vanish for OU). The generic 1.5 is gated on a
   cubic-well SDE (measured 1.609); the OU superconvergence is gated at 2.0
   and documented. The spec's "SRA1 ≈ 1.5" is confirmed where it
   mathematically applies.
2. **Displacement cap applies to the drift contribution only** (born.py
   capped the full RK4 displacement): capping the Wiener increment would
   bias the diffusion coefficient itself; walls handle noise overshoot by
   specular reflection (which born.py's clip cannot, once noise can cross a
   wall in one substep).
3. **EM (not RK4) per substep for the Nelson demo**: born.py's RK4 is a
   deterministic-ODE integrator inapplicable to the SDE; EM is strong order
   1.0 for this additive-noise system, and the H̄ diagnostics are
   statistical. SRA1 is available behind the same trait but buys nothing
   near nodal lines where the drift is only piecewise smooth.
4. **Fine early-time τ sampling + floor-subtracted fit** (see protocol
   caveats above): the born.py cadence and raw-log fit under-resolve a decay
   37× faster than the one they were designed for. Both protocols are
   reported; the headline comparison uses born.py's.
5. **N kept at 20 000** (no reduction needed): demo wall time ≈ 309 s total,
   well under the 15-minute budget.
6. **Phases embedded, sampler not emulated**: born.py's 16 phase values
   (numpy `default_rng(42)`) are embedded as exact round-trip decimal
   constants so ψ is the SAME wavefunction; the ensemble samplers use
   fs-rand streams (per-trajectory rejection with deterministic
   consumption), so ensembles differ at resampling level from numpy's.
7. **Evidence-stack certification not wired** (spec marks it optional): the
   determinism story is carried by the G1 gates + the cross-run golden
   fingerprint above; the fs-evidence/fs-package/fs-checker wrapper can be
   added later without touching the physics.

## Epistemic notice

Within-model validation only. This crate validates the ENGINE — the SDE
integrators against analytic referees, and the H-theorem machinery against
Tier-3's qualitative gate. The dynamical-relaxation H-theorem is the
corpus's [IM] imported ingredient (Valentini-class); reproducing it under
Nelson dynamics, and measuring the Nelson-vs-Bohm relaxation-rate ratio,
says nothing about GUM's substrate ontology (see Tier-3 REPORT.md). The
τ_Nelson ≪ τ_Bohm result is a property of the two borrowed dynamics on this
system, offered as engine-level information for the corpus comparison.
