# N3 — L7 continuum proof attempt of the cutoff zeros: the field-level crossing criterion (F-T7-N3)

**Goal.** Elevate L7b's sampled hard-cutoff zeros (tau_max = 5.130950 at
n_y = +1; two exactly-empty kill bins [6.4, 6.8], [6.8, 7.2] at |n_y| >= 0.75)
from ensemble statements to a FIELD statement via the pre-registered crossing
criterion: measure v_x(d, z, t) on the detector line d = d_near = 1 over the
full (z, t) grid, find tau* = sup{t : max_z v_x >= 0}, and verify the margin
sup_z v_x < 0 on (tau*, T] at proof-grade-modulo-discretization.

**Campaign:** Tier 7 Phase N (ROADMAP_v9_PROGRAM.md addendum, N3).
**Date:** 2026-07-18. **Code:** `n3_field.py` (runtime 427 s), `n3_fig.py`.
**Raw output:** `n3_results.json`, `n3_raw.npz`. **Figure:** `n3_fig.png`.
**Prior art (mandatory reads, honored):** tier6-foundations/L5prime/RESULTS.md +
l5prime_run.py, tier6-foundations/L7b/RESULTS.md + l7b_run.py — engine reused
verbatim (exact-in-z factorized evolution, C2 smoothstep taper, k0 = 2,
waveguide z in [0, 1], box x in [-15, 45], Nx = 2048 reference, dt = 2.5e-4
<= 5e-4, both raw archives loaded for the consistency checks), including the
exact L7b velocity decomposition v(x, z, t; n_y) = a(x, z, t) + n_y b(x, z, t).

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**NEGATIVE RESULT — theorem intact, hypothesis refuted:** the crossing
criterion is true and two lines long, but its hypothesis
sup_z v_x(d, z, t) < 0 is FALSE at every t in (0, T] — on the full line and
on every rho-floored region M(t) = max_z v_x >= 0 at ALL 1601 output times
(tau* = T = 16, no negative window exists), and even the single most
favorable point z = 1/2, where the spin term vanishes and
v_x(d, 1/2, t) = A(t) is the SAME for every n_y, stays positive through both
kill bins (A in [0.857, 0.964]) and recurs >= 0 until t = 15.87 — forced
independently by L5prime's own 219 sampled axial arrivals up to 15.787, each
of which requires A >= 0 at the line. The wrong sign exceeds the measured
v_x discretization bound by ~3.5e4. The L7b cutoff is therefore a TRANSPORT
phenomenon (the not-yet-crossed set empties by tau_max), not a
local-velocity-sign phenomenon; no line-local certificate exists in this
geometry. The kill-bin zeros themselves are STRENGTHENED: a 60,000-trajectory
deterministic covering gives 0 kill-bin crossings, with grid support edges
5.142 (|n_y| = 1) and 6.017 (|n_y| = 0.75), far below the first kill bin at
6.4. The fully analytic proof (T4-W5 residual) remains open and is now known
to require a flow-map argument.

## The theorem (stated and proved — it is two lines)

**Crossing criterion.** Let (X(t), Z(t)) solve dX/dt = v_x(X, Z, t),
dZ/dt = v_z(X, Z, t). If X first crosses x = d at t0 (X(t0) = d, X(t) < d
for t < t0), then v_x(d, Z(t0), t0) >= 0: for every h > 0,
(X(t0) − X(t0 − h))/h >= 0, and the h → 0+ limit is dX/dt(t0) =
v_x(d, Z(t0), t0). ∎

**Contrapositive.** If sup_z v_x(d, z, t) < 0 for all t in (tau*, T], no
trajectory first-crosses x = d in (tau*, T]; the continuum first-crossing
density vanishes there and its support within [0, T] lies in [0, tau*]. ∎

Density weighting, handled as pre-registered: v_x is defined wherever
rho > 0; rho(d, z, t) = 2|f(d, t)|^2 sin^2(pi z) vanishes only at the walls
(plus possible isolated zeros of f), so the sup is evaluated on the full
open line AND on the rho-floored regions {sin^2(pi z) >= eps},
eps = 1e-1, 1e-2, 1e-3 (wall-complement |phi|^2-measure ∝ (4pi^2/3)z_fl^3,
negligible). The subtlety turns out immaterial: the hypothesis fails ON the
floored regions and at mid-channel, not in the low-rho zones.

## The structural identity that decides the workstream

For the factorized state phi = f(x, t) sin(pi z) e^{−i pi^2 t/2}, the L7b
velocity decomposition evaluated ON the line x = d separates exactly:

    v_x(d, z, t; n_y) = A(t) + n_y S(z),
    A(t) = Im(f'/f)(d, t)  [z-independent],  S(z) = −pi cot(pi z)  [t-independent].

Two consequences, both fatal to the criterion in this configuration:

1. **Mid-channel:** S(1/2) = 0, so v_x(d, 1/2, t) = A(t) for EVERY n_y.
   A(t) is the axial (n_y = 0) line velocity, and L5prime's G3 measured 219
   axial first crossings of this very line at times up to 15.787, each
   forcing A >= 0 there. The present field measurement confirms it directly:
   A(t) > 0 on ALL of (0, 14.0] — including the entire kill-bin band
   (min 0.857) — nonnegative on 96.9% of [0, 16], last nonnegative time
   t = 15.87.
2. **Wall branch:** S(z) → +∞ as z → 1 (n_y > 0) resp. z → 0 (n_y < 0), so
   on any rho-floored region max_z v_x = A(t) + |n_y| pi cot(pi z_floor),
   with additive constants 9.42|n_y| at eps = 1e-1, 31.1|n_y| at 1e-2,
   95.2|n_y| at 1e-3 (measured kill-bin minima of M at n_y = +1: 10.27 /
   32.00 / 96.06, i.e. A_min + the constant) — positive at every t.
   tau* = T identically on every admissible region and for every n_y in
   {+1, +0.75, −0.75, −1, +0.25}.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | field quality at the L5prime bars (norm 1e-6; walls < 1e-10; energy 1e-6) | line evolution (exact propagator): norm dev 6.7e-16, energy drift 3.9e-16 (rel); 2-D wall diagnostic: norm dev 2.2e-16, wall max 3.4e-16, energy drift 1.2e-16, z-mode leakage 0.0; E_2D = 7.18520 = E_x (2.25040) + pi^2/2 | **PASS** |
| G2 | tau* exists for n_y = +1 and is stable <= 1% under Nx-doubling AND dt-halving | tau* does NOT exist below T: M(t) >= 0 at every one of the 1601 output times on the full line (M pinned at the +500 clamp by the z → 1 wall branch) and on every rho-floored region (M >= A + 10.3 at the mildest floor); tau* = T = 16.0, identically at Nx x 2 and dt/2 (rel 0.0 — stable but vacuous). Most favorable measure-zero restriction (z = 1/2 alone): last-nonnegative time 15.87, identical at both refinements — 3.1x the sampled cutoff, not an estimate of it | **FAIL** |
| G3 | consistency: tau* >= 5.130950 with the gap reported and explained (sampled max <= continuum edge); tau* < sampled tau_max = CONTRADICTION-FAIL | No contradiction in any variant: tau* = 16.0 (full/floored) resp. 15.87 (mid-only) >= 5.130950; gap 10.869 resp. 10.739. But the gap is NOT the anticipated sampled-max-below-continuum-edge sliver: tau* is pinned near T by the n_y-independent mid-channel velocity and the positive wall branch, so it does not estimate the support edge — the criterion is vacuous above 5.131. The real sampled-vs-continuum gap appears in the covering diagnostic instead: grid support edge 5.1422 vs sampled 5.13095 (+0.011, +0.22%), the expected direction, and still 1.26 below the first kill bin | **PARTIAL** |
| G4 | margin: delta_min on (tau* + 0.1, T] >= 10x the v_x discretization bound | Window empty (tau* = T). Where the gate wants delta = −M > 0, the field gives M >= A >= 0.857 in the kill bins — the WRONG SIGN by ~3.5e4 x the measured v_x line bound (Nx-doubling + dt-halving, pointwise sum: kill-bin median 9.2e-6, p95 2.2e-5, max 2.5e-5; over all t > 5.131: median 1.7e-5, max 3.4e-3 at quasi-nodal passages). The hypothesis is refuted, not unverified, at proof-grade-modulo-discretization | **FAIL** |
| G5 | discrimination: n_y = +0.25 shows no tau* below T | Literal clause holds: at +0.25, M(t) >= 0 recurring to T (tau* = T full-line; mid-restriction 15.87). But n_y = +1 behaves IDENTICALLY under the criterion (same tau* = T, same mid value 15.87, since v_x(d, 1/2, t) is n_y-independent): the criterion does not separate the cutoff member from the no-cutoff member, which was this gate's purpose | **PARTIAL** |

## Key numbers

- Line velocity A(t) = v_x(d = 1, z = 1/2, t) (reference Nx = 2048,
  dt_out = 0.01, 1601 samples): A > 0 on all of (0, 14.0]; nonnegative on
  96.9% of [0, 16]; kill-bin range [0.857, 0.964]; last nonnegative time
  15.87 (identical at Nx = 4096 and dt = 1.25e-4); A(16) = −1.90;
  rho(d, t = 16) = 9.6e-4.
- tau* (n_y = +1): full line 16.0 = T; floors 1e-1/1e-2/1e-3: 16.0/16.0/16.0;
  mid-channel-only 15.87. Same values for n_y = +0.75, −0.75, −1, +0.25.
  Kill-bin minima of M(t), n_y = +1: full 500 (clamp); floors 10.27 / 32.00 /
  96.06; mid 0.857.
- v_x discretization bound on the line (pointwise |A_2048 − A_4096| +
  |A_dt − A_dt/2|): kill bins median 9.2e-6 / p95 2.2e-5 / max 2.5e-5
  (Nx term dominant; dt term < 8e-12 — the propagator is exact, dt enters
  only through float accumulation); t > 5.131 overall: median 1.7e-5,
  max 3.4e-3 (isolated quasi-nodal spikes). Wrong-sign ratio in the kill
  bins: 0.857 / 2.5e-5 = 3.5e4.
- Consistency with sampled statistics: sampled tau_max(+1) = 5.130950
  (L5prime) / 5.130946 (L7b) <= tau* in every variant. A interpolated at
  the 219 sampled late axial arrival times: 98.6% >= 0 (99.85% >= −0.05
  over all 1976 axial arrivals; min −0.29 at a quasi-nodal oscillation,
  within the RK4-crossing vs output-grid interpolation tolerance) — field
  and trajectories agree. Sampled axial last arrival 15.787 <= 15.87 =
  tau*_mid: for n_y = 0 the criterion IS valid (v_z = 0 freezes z, so
  sup_z v_x = A) and delivers a true field statement — axial continuum
  first-crossing support within [0, 16] ⊂ [0, 15.87], margin
  delta_min = 1.87 on (15.97, 16] vs bound < 3.4e-3 — sharp to 0.5% of the
  sampled endpoint, exactly where the spin term is absent.
- Covering diagnostic (deterministic 150 x 100 cell-centre grid on the
  compact support, 15,000 trajectories per n_y, dt = 2.5e-4, 424 s):

  | n_y | crossed | grid tau_max | beyond sampled tau_max | kill-bin crossings | non-crossers (weight) |
  |---|---|---|---|---|---|
  | +1 | 14,857 | 5.1422 | 1 pt (w = 2.0e-4): (x0, z0) = (−5.21, 0.215), tau = 5.1422 | **0** | 143 (2.2e-9) |
  | −1 | 14,857 | 5.1422 | 1 pt (mirror z0 = 0.785) | **0** | 143 (2.2e-9) |
  | +0.75 | 14,699 | 6.0174 | 0 | **0** | 301 (4.8e-7) |
  | −0.75 | 14,699 | 6.0174 | 0 | **0** | 301 (4.8e-7) |

  Weighted 99.9th-percentile crossing times 5.105 (|n_y| = 1) / 5.959
  (|n_y| = 0.75). The +/− mirror pairs are numerically identical
  (z → 1 − z symmetry) — an engine self-check. Non-crossers are thin
  near-wall / support-edge strips of negligible |phi0|^2 weight that drift
  backward and never cross (they produce no late arrivals and cannot touch
  the kill-bin statement).
- Engine: total runtime 427 s; criterion z-grid 1999 points; all six n_y
  evaluated from one A(t) per resolution (the decomposition is exact).

## Method

Engine reused verbatim from L7b (`l7b_run.py`): factorized exact-in-z
state, f propagated exactly in 1-D Fourier space on [-15, 45], Nx = 2048
reference; a = Im(f'/f), b = Re(f'/f) with the engine's density floor
(1e-14 relative) and clamp (±500); the line field A(t) is a(x, t) linearly
interpolated to x = d exactly as the RK4 integrator sees it. Three error
probes for v_x on the line: exact spectral evaluation at each of the 1601
output times (the propagator is exact — no time-stepping error exists for
the field), the engine's cumulative half-step accumulation at dt = 2.5e-4
and dt = 1.25e-4 (float-accumulation bound, i.e. the dt-halving probe),
and Nx = 4096 (spatial + interpolation bound, the Nx-doubling probe).
v_x(d, z, t; n_y) assembled from the exact decomposition A(t) + n_y S(z)
on a 1999-point open z-grid, clamped as in the engine; M(t) maximized over
the full line, the rho-floored regions {sin^2(pi z) >= eps}, and z = 1/2
alone; tau* = last output time with M >= 0; margins delta = −M where a
negative window exists. Covering diagnostic: the L7b trajectory integrator
(RK4 on stored fields at t, t + dt/2, t + dt; crossing interpolated within
the step; escapees frozen at the box rim as in L7b) run from a
deterministic 150 x 100 cell-centre grid instead of the |phi0|^2 sample,
n_y in {−1, −0.75, +0.75, +1} advanced in one field evolution. 2-D FFT
wall diagnostic as in L7b (Nx = 2048, Nz = 128, 17 checkpoints).
Consistency data: L5prime_raw.npz and l7b_raw.npz loaded unmodified.
hbar = m = 1.

## Diagnosis: why the criterion cannot work here — and what would

1. **The detector line is ahead of the source.** For free spreading from
   compact support around x0 = −5, the late-time line velocity approaches
   the ballistic field A(t) ≈ (d − x0)/t = 6/t > 0; measured A stays
   O(0.1–1) and strictly positive to t = 14, with taper-lobe interference
   producing brief negative excursions only beyond ~14. At ANY line ahead
   of the source, sup_z v_x < 0 at late times is impossible in this
   geometry. A hard first-crossing cutoff at a forward detector is
   therefore necessarily a statement about where the not-yet-crossed
   trajectories ARE (transport), never about the sign of v_x on the line.
2. **The spin term cannot help; it hurts.** On the line it is the
   t-independent profile n_y S(z), odd about mid-channel with a
   positive-divergent branch at one wall for every n_y ≠ 0: it can only
   raise sup_z v_x, and it vanishes at z = 1/2 where the shared axial
   velocity already recurs positive. The criterion is structurally blind
   to what actually produces the cutoff: the z-advection
   v_z = n_y Re(f'/f), which lifts laggards into the fast wall lane and
   empties the not-yet-crossed set by tau_max (the Das–Dürr mechanism as
   realized in L5prime/L7b).
3. **Where the criterion IS valid it is sharp** (n_y = 0: z frozen,
   sup_z v_x = A, proved support ⊂ [0, 15.87] vs sampled 15.787, margin
   1.87 vs bound 3.4e-3) — the field computation is sound; the failure at
   |n_y| > 0 is physics, not numerics.
4. **What a continuum proof now requires:** a flow-map / characteristics
   argument — certified integration of the not-yet-crossed set's boundary,
   or a Lipschitz covering of the initial support. The 60,000-trajectory
   deterministic covering here is the uncertified numerical version and
   confirms the zeros (0 kill-bin crossings; support edges 5.142 / 6.017
   vs first kill bin at 6.4). The T4-W5 residual is thereby sharpened from
   "prove the zeros" to "prove them by transport — the line criterion
   provably cannot".

## Honest scope

This workstream proves the pre-registered continuum support statement only
for the n_y = 0 family member (GIVEN the computed field; bounds printed).
For the kill set |n_y| >= 0.75 it establishes the opposite of the
pre-registered hope, at the same evidentiary grade: the criterion's
hypothesis is false at every t in (0, T] (wrong sign by ~3.5e4 x the
discretization bound), so no line-local certificate of the cutoff exists in
this configuration. The kill-bin zeros stand as ensemble statements now
corroborated by a deterministic covering — a field-like but uncertified
upgrade, not the theorem the workstream sought. The fully analytic proof
remains the named open item (T4-W5), with its required shape identified.
Within-model; nothing here bears on nature.

## Caveats

- All statements are about the computed field of the validated engine at
  the printed resolutions; the discretization bound is an Nx-doubling +
  dt-halving difference estimate, not a rigorous interval enclosure.
  Isolated quasi-nodal spikes of A(t) inflate the bound's max (3.4e-3);
  medians/p95 are quoted alongside, and no conclusion is within 3 orders
  of the bound.
- tau* values are read on the dt_out = 0.01 output grid (±0.01
  granularity); the quoted refinement stability (identical values) is at
  that granularity.
- The 96.9%/15.87 recurrence structure of A(t) and the covering's grid
  support edges are configuration facts of this bench (k0 = 2, literal
  amplitude, box [-15, 45] periodic); L5prime's box-control showed the
  wrap does not move gate-relevant near-field quantities, and the kill-bin
  window t <= 7.2 is far inside the wrap-free regime.
- The covering is reference-resolution RK4, not certified integration; its
  non-crossers (143–301 per n_y, |phi0|^2 weight <= 4.8e-7) are near-wall /
  support-edge strips outside the sampled ensemble's range — they produce
  no crossings and cannot affect the kill-bin statement, but they mark the
  covering's resolution limit at the walls.
- The single grid point beyond the sampled tau_max (tau = 5.1422 at
  |n_y| = 1, weight 2.0e-4) means the continuum support edge exceeds the
  sampled max by at least 0.011 — consistent with L7b's caveat that
  tau_max is an order statistic, and irrelevant to the kill bins (edge
  still 1.26 below them).
- The criterion's failure is specific to this geometry (detector ahead of
  a forward-boosted source in a straight guide); a detector line behind
  the source is not covered by this diagnosis and remains a conceivable
  home for a line-local certificate.
