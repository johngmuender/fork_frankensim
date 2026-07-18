# O1 — the flow-map first-crossing proof: backward-reachability criterion on the kill bins (F-T7-O1)

**Goal.** Execute the terminal proof shape identified by F-T7-N3 (which proved
the line criterion cannot work and that the L7b cutoff is a transport
phenomenon): elevate the L7b kill-bin zeros (t in [6.4, 6.8] and [6.8, 7.2],
|n_y| >= 0.75, 0/2000 arrivals) to a FIELD/FLOW statement by backward-
integrating the guidance ODE from EVERY kill-bin detector point and showing
each one already crossed the detector line earlier — so none can carry a
first crossing (arrival) in the kill bins.

**Campaign:** Tier 7 Phase O (ROADMAP_v9_PROGRAM.md addendum, O1).
**Date:** 2026-07-18. **Code:** `o1_run.py` (runtime 585 s), `o1_fig.py`.
**Raw output:** `o1_results.json`, `o1_raw.npz`. **Figure:** `o1_fig.png`.
**Prior art (mandatory reads, honored):** tier7-program/N3/RESULTS.md +
n3_field.py (the refuted line criterion; engine-reuse pattern),
tier6-foundations/L7b/RESULTS.md + l7b_run.py, L5prime/l5prime_run.py —
engine reused verbatim (exact-in-z factorized evolution, C2 smoothstep
taper, k0 = 2, waveguide z in [0, 1], box x in [-15, 45], Nx = 2048
reference, dt = 2.5e-4 <= 5e-4, clamp ±500, density floor 1e-14 relative),
including the exact velocity decomposition
v(x, z, t; n_y) = a(x, z, t) + n_y b(x, z, t).

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**ALL FIVE GATES PASS — the kill-bin zeros are elevated to a FIELD/FLOW
statement at proof-grade-modulo-discretization:** of the 64,000 kill-bin
detector points (200 z x 80 t x 4 kill n_y), 63,984 back-trace across
x = d at a certified earlier epoch (class (a), min margin 0.621 =
22,068 x the p95 backward-path discretization error), the remaining 16
are extreme wall cells whose total crossing flux is <= 2.8e-7 of the
|phi_0|^2 measure (3.5 orders below the 1e-3 gate), and there are ZERO
violations — no backward path reaches t = 0 inside supp(phi_0) without
pre-crossing — identically at reference, doubled Nx, and halved dt
(classification agreement 80,000/80,000 in both refinements). The
first-crossing (arrival) measure of both kill bins is therefore
<= 2.8e-7 per kill member GIVEN the computed field — a ~5,000-fold,
field-level sharpening of L7b's sampled 95% CL bound 1.5e-3. The
discrimination control confirms the criterion sees physics, not an
artifact: at n_y = +0.25 the same backward map classifies 74% of the
same bins as legitimately fresh, and the fresh cells' flux (0.0187)
reproduces L7b's measured bin probability (0.0235) to 20%. The fully
analytic proof (T4-W5) remains open; its terminal shape — this
backward-reachability argument — is now demonstrated end to end.

## The theorem (stated and proved)

**Backward-reachability first-crossing criterion.** Let v be the guidance
velocity of the bench member n_y,

    v_x = Im(f'/f)(x, t) − n_y π cot(π z),   v_z = n_y Re(f'/f)(x, t),

on U = {ρ > 0}, ρ = 2|f|² sin²(π z). For every t > 0, f(·, t) is the free
evolution of a compactly supported L² ∩ L¹ datum, hence (Paley–Wiener,
applied to the Fresnel representation f(x,t) = (2πit)^{-1/2} e^{ix²/2t}
∫ e^{-ixy/t} [e^{iy²/2t} f₀(y)] dy) real-analytic in x; v is therefore C^∞,
in particular locally Lipschitz, on U.

**(1) Two-sided uniqueness.** By Picard–Lindelöf, through every point
(p, t) ∈ U there is exactly one maximal integral curve of dX/ds = v(X, s) —
unique FORWARD and BACKWARD in time. Let X(s) = Φ_{s,t}(p) denote it.

**(2) Pre-crossing kills the arrival.** Fix the detector line x = d and a
detector point p = (d, z) at time t. Suppose the backward segment satisfies
max_{s<t} X_x(s) ≥ d. If X_x(0) ≥ d the unique trajectory through (p, t)
starts at or beyond the line and has no first crossing at t at all.
Otherwise X_x(0) < d ≤ X_x(s₀) for some s₀ < t, and by continuity and the
intermediate value theorem the trajectory reaches x = d at some
s* ≤ s₀ < t: its FIRST crossing time is ≤ s* < t. Either way, by (1) this
is the ONLY trajectory through (p, t), so the point (z, t) carries no
first-crossing (arrival) density. ∎

**(3) Flux domination of the exceptional set.** The first-crossing density
is dominated by the upward crossing flux: a first crossing at (z, t) is in
particular an upward crossing of the line (v_x ≥ 0 there — N3's crossing
lemma), and the equivariant transported ensemble density is ρ, so

    dP_first/(dz dt) ≤ ρ(d, z, t) · max(v_x(d, z, t), 0).

**(4) Conclusion.** If every detector point of B = (0,1) × (kill bins) is
of type (2) except a subset N, then the first-crossing measure of B is
≤ ∫_N ρ (v_x)₊ dz dt. If the bound is zero (negligible), the arrival
density in the kill bins vanishes (to that bound) at continuum/field
level. ∎

This is exactly the flow-map/transport argument N3 showed is required: it
certifies where the not-yet-crossed set IS (nowhere near the kill bins),
not the sign of v_x on the line (which N3 proved has the wrong sign).

## Execution

From every detector point (d = 1, z_i, t_j) — z on the 200-point cell-centre
grid of (0, 1), t on 40 cell centres per kill bin (spacing 0.01), per
n_y ∈ {+1, −1, +0.75, −0.75} and the control +0.25 (80,000 backward
trajectories total per resolution) — the SAME field's guidance ODE is
integrated backward (RK4 on fields at s, s − dt/2, s − dt from the exact
propagator, engine clamps and z-clip identical to l7b_run.py) from s = t
down to s = 0, and the start point is classified:

- **(a) PRE-CROSSED** — max_{s<t} X_x ≥ d + 10⁻³ before any ρ-floor event
  (margin and certified crossing epoch s_pc recorded);
- **(b) RHO-FLOOR** — the path enters ρ < 10⁻¹⁰ × ρ_peak(s) (near-node /
  near-wall zones where the discrete field is ill-conditioned); its
  possible first-crossing measure is bounded by its detector-cell flux
  ρ (v_x)₊ Δz Δt;
- **(c) VIOLATION** — the path reaches s = 0 inside supp φ₀ with
  max X_x < d + 10⁻⁹ throughout — a genuine fresh kill-bin arrival, which
  would REFUTE the L7b zeros;
- plus the bookkeeping classes GRAZE-UNCERTAIN (10⁻⁹ < margin < 10⁻³),
  FRESH-OUT (reaches s = 0 outside supp φ₀, where ρ₀ = 0), BOX-EXIT
  (reaches the box rim) — all three added to the (b)-style flux bound.

Refinements over the FULL grid: Nx = 4096 (doubled) and dt = 1.25e-4
(halved backward step).

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | engine bars, 1e-6 class or better | 2-D wall diagnostic: norm dev 2.2e-16, wall max 3.4e-16, energy drift 1.2e-16 (rel), z-mode leakage 0.0; backward field evolutions (all three resolutions): norm dev <= 2.0e-12, energy drift <= 3.2e-13, and backward REVERSIBILITY — F accumulated from s = 7.195 back to 0 vs F0 — <= 6.1e-12 relative | **PASS** |
| G2 | zero VIOLATIONS at reference AND under Nx-doubling + dt-halving; classification stable | 0 violations for every kill n_y at ref, Nx = 4096, and dt = 1.25e-4 (12 zero cells); classification agreement 80,000/80,000 points in BOTH refinements (not a single label moved, control included); fresh count at +0.25 identical (11,833) at all three resolutions | **PASS** |
| G3 | min class-(a) margin (max X_x − d) >> backward-path discretization error, measured via refinement differences of the paths | min certified margin 0.6207 (at \|n_y\| = 0.75; the \|n_y\| = 1 families are entirely at the cap 2.0). Per-path margin differences over the 9,996 uncapped class-(a) paths: Nx-doubling median 8.1e-6 / p95 2.8e-5 / max 6.0e-3; dt-halving median 2.6e-9 / p95 4.0e-8 / max 4.8e-5. Separation: 0.6207 / 2.8e-5 = **22,068** (vs required 10) | **PASS** |
| G4 | measure reachable through class-(b) (+ graze / out-supp / box-exit) <= 1e-3 of \|phi_0\|^2 | ZERO rho-floor events before classification (class (b) empty at 200-z resolution: z_min = 0.0025 keeps start cells 5 orders above the 1e-10 floor); the entire bounded budget is 16 box-exit paths (extreme wall cells z = 0.0025/0.9975, backward exit via the LEFT rim x = −14.5): flux 2.79e-7 (\|n_y\| = 1) / 1.53e-7 (0.75) per family — 3.5 orders below the gate | **PASS** |
| G5 | discrimination: n_y = +0.25 backward runs from its POPULATED bins show a substantial legitimately-fresh fraction | 11,833 / 16,000 grid points (74.0%) legitimately fresh (reach t = 0 inside supp phi_0 with max X_x < d + 1e-9); fresh-cell flux per bin 0.0104 / 0.0083 vs L7b measured Pi 0.0120 / 0.0115 — total ratio 0.795 (>= 0.3 gate) — the backward map quantitatively reproduces the mechanism that populates those bins for the no-cutoff member | **PASS** |

## Key numbers

- Classification, reference run (per kill n_y; 16,000 detector points each;
  200 z x 40 t per bin, bins [6.4, 6.8] and [6.8, 7.2], d = 1):

  | n_y | (a) PRE-CROSSED | (b) RHO-FLOOR | (c) VIOLATION | graze | box-exit | bounded flux |
  |---|---|---|---|---|---|---|
  | +1.00 | 15,995 | 0 | **0** | 0 | 5 | 2.79e-7 |
  | −1.00 | 15,995 | 0 | **0** | 0 | 5 | 2.79e-7 |
  | +0.75 | 15,997 | 0 | **0** | 0 | 3 | 1.53e-7 |
  | −0.75 | 15,997 | 0 | **0** | 0 | 3 | 1.53e-7 |
  | +0.25 (control) | 4,165 | 0 | 11,833 fresh (legit) | 2 | 0 | 3.68e-7 |

  The +/− mirror pairs are numerically identical (z → 1 − z symmetry of the
  guidance law with the z-symmetric preparation) — an engine self-check, as
  in N3's covering.
- Certified pre-crossing epochs s_pc (the time the backward path is first
  a full margin 1e-3 beyond the line): |n_y| = 1: min 3.599 / median 5.700 /
  max 7.195; |n_y| = 0.75: min 2.383 / median 4.876 / max 7.195 — all deep
  in the trusted mid-evolution regime; the late-epoch points are those with
  v_x(d, z, t) < 0 (N3's z < z*(t) branch), whose backward paths sit above
  the line immediately.
- Margins: |n_y| = 1: every class-(a) path reaches the cap 2.0 (min 2.0000);
  |n_y| = 0.75: 72.9% capped, uncapped minimum 0.6207. Margin histogram and
  the discretization comparison in the figure (panel C).
- Field-level bound: first-crossing measure in the two kill bins
  <= 2.8e-7 per kill member (the bounded-class flux), vs L7b's ensemble
  95% CL of 1.5e-3 from 0/2000 — a ~5,400x sharpening, and a statement
  about the CONTINUUM flow, not a sample.
- The zeros are not for lack of line flux: total upward flux rho (v_x)+
  through the kill bins is 0.0302 (|n_y| = 1) / 0.0257 (0.75) per family —
  ALL of it re-crossing flux of trajectories that first crossed earlier
  (backward-verified), none of it first-crossing flux.
- Control quantitative consistency: fresh region is the band z >= 0.2425
  (mixed near the far wall), pre-crossed band z <= 0.24 — the slow lane
  that n_y = +0.25's weak z-advection has not yet emptied; fresh flux
  0.010401 + 0.008285 = 0.018685 vs L7b Pi = 0.0120 + 0.0115 = 0.0235
  (ratio 0.795; the deficit is the finite-cell discretization of the
  flux integral and L7b's N = 2000 sampling). The 2 graze-uncertain
  points (margins 2.8e-4, 9.9e-4) sit exactly on the fresh/pre-crossed
  boundary z ≈ 0.26, as they must.
- Runtime: ref 145 s + Nx-doubled 148 s + dt-halved 291 s = 585 s total;
  80,000 backward trajectories per resolution, all five n_y families in
  one field evolution each.

## Method

Engine reused verbatim from L7b (`l7b_run.py`) as in N3: factorized
exact-in-z state, f propagated exactly in 1-D Fourier space on [-15, 45];
a = Im(f'/f), b = Re(f'/f) with the engine's density floor (1e-14
relative) and clamp (±500); v_x re-clamped after the spin term; z clipped
to [1e-9, 1 − 1e-9]. Backward fields by half-step accumulation with the
conjugate phase factor from the exact F(s_start); backward reversibility
(F returned to F0 at s = 0) is a gate bar. Backward RK4 mirrors the
forward integrator stage-for-stage (fields at s, s − dt/2, s − dt). All
five n_y families advance in one field evolution per resolution.
Activation times are exact multiples of dt (cell centres at spacing 0.01;
s_start = 7.195). Trajectories freeze at: a ρ-floor event, the box rim, or
once the pre-crossing margin exceeds the cap 2.0 (classification settled;
recorded margin then understates the true max — conservative for G3).
Classification precedence: a certified pre-crossing only counts if NO
ρ-floor event precedes it on the backward path, so the certified segment
[s_pc, t] lies entirely in the trusted region. Detector-cell flux weights
ρ(d, z, t) max(v_x, 0) Δz Δt use the same interpolated, clamped fields the
integrator sees. 2-D FFT wall diagnostic as in L7b. hbar = m = 1.

## Honest scope

This is the pre-registered verdict target achieved at its stated grade and
no more: the kill-bin zeros are now a FIELD/FLOW statement — every kill-bin
detector point of the continuum flow, up to a quantified 2.8e-7 flux
residue, is backward-reachable across the detector line and therefore
carries no first crossing — GIVEN the computed field and RK4 flow map at
the printed resolutions, with stability certified by a full-grid
Nx-doubling and dt-halving in which not one of 80,000 classifications
moved. It is proof-grade MODULO DISCRETIZATION: the field is the exact
spectral propagator (error at the 1e-12 reversibility bar), but the
backward trajectories are RK4 on interpolated fields, not certified
(interval) integration, and the pre-crossing margins are compared against
a refinement-difference error ESTIMATE (x22,068 separation), not a
rigorous enclosure. The fully analytic T4-W5 proof remains open — what
this workstream settles is its SHAPE: N3 proved no line-local certificate
exists; O1 demonstrates that the backward-reachability/transport argument
closes the gap, and that nothing about the flow resists it (no rho-floor
quarantine was even needed). Within-model; nothing here bears on nature.

## Caveats

- Not certified integration: RK4 + linear-in-x interpolation on the
  discrete field; the guarantee is refinement stability (identical
  classification; per-path margin differences median 8.1e-6, p95 2.8e-5),
  not interval arithmetic. The max margin difference 6.0e-3 (one
  quasi-nodal passage) is still 100x below the min margin.
- The margin cap 2.0 freezes certified class-(a) paths early; recorded
  margins understate the true backward max X_x (conservative for G3), and
  frozen paths are not integrated to s = 0 (irrelevant to the theorem:
  the certified segment [s_pc, t] is complete and floor-free).
- Class (b) is EMPTY at this detector grid because the 200-point z-grid's
  wall cells (z = 0.0025) sit ~5 orders above the 1e-10 rho floor and no
  interior backward path entered a floored zone before classification; at
  a much finer z-grid the extreme wall cells would populate class (b)
  instead of box-exit — same accounting, same negligible flux.
- The 16 box-exit paths leave via the LEFT rim (x = −14.5, behind the
  source) under the near-wall clamped velocity; their detector cells carry
  sin^2(pi z) ~ 6e-5 weight and are bounded by flux, not classified.
  Box-exit via the right rim is impossible without prior certified
  pre-crossing (any path at x > d + 1e-3 is already class (a)).
- The engine's velocity regularization (clamp +-500, density floor 1e-14
  relative) is inherited verbatim from L5prime/L7b; near the walls the
  clamped field is not the exact j/rho, which is why wall-adjacent
  bookkeeping is flux-bounded rather than trusted. The flux weights
  themselves use the same clamped field the integrator sees.
- The within-step ordering (pre-cross checked before floor at the same
  step) is a one-step ambiguity; with zero floor events it never fired.
- G5's flux-vs-Pi agreement (0.795) is a consistency check between a
  cell-discretized flux integral and an N = 2000 binomial estimate
  (Pi 95% intervals ~ +-40%); it is corroboration, not a calibration.
- Configuration-specific: d = 1 ahead of the source, k0 = 2, literal
  amplitude, box [-15, 45]; the kill-bin window t <= 7.2 is far inside
  the wrap-free regime established by L5prime's box control.
