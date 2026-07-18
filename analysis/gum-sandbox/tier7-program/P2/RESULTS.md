# P2 — T4-W5 certification pilot: exact-field, high-precision backward paths (F-T7-P2)

**Goal.** Upgrade O1 (F-T7-O1) one rung toward certified: the L5'/L7b/O1
engine's field is a FINITE trigonometric sum (analytic sine mode in z,
band-limited 2048-mode Fourier representation in x), so psi, rho, j — hence
the guidance velocity v — can be evaluated EXACTLY (to floating-point
rounding) at ARBITRARY points, eliminating O1's dominant error source
(linear grid interpolation). Re-run the backward classification for a
stratified 432-path subset with (i) high-order adaptive integration
(DOP853, rtol 1e-12) and (ii) a multiprecision ladder (200-bit ~ 60-digit
>= 50-digit fixed-step RK4 on the exact mode sum, 24 paths incl. the worst
margins).

**Campaign:** Tier 7 Phase P (ROADMAP_v9_PROGRAM.md addendum, P2).
**Date:** 2026-07-18. **Code:** `p2_certify.py` (phases g1 / select / dop /
ladder / gates / fig, checkpointed). **Raw output:** `p2_results.json` +
intermediates (`_g1.json`, `_select.json`, `_selection.npz`, `_dop.npz`,
`_ladder.json`, `_ladder_idx.npy`). **Figure:** `p2_fig.png`.
**Prior art (mandatory reads, honored):** tier7-program/O1/RESULTS.md +
o1_run.py + o1_raw.npz (the backward-reachability run this pilot upgrades;
engine constants inherited verbatim from o1_run.py / l7b_run.py /
l5prime_run.py: k0 = 2, box [-15, 45], Nx = 2048, d = 1, clamp ±500,
density floor 1e-14 relative, z-clip 1e-9; kill bins [6.4, 6.8], [6.8, 7.2]).

*Within-model; nothing here bears on nature.*

## 0. Verdict in one line

**FOUR OF FIVE GATES PASS; G3 FAILS on its endpoint clause only — and the
certification content is achieved:** the exact band-limited field evaluator
agrees with the grid engine at machine rounding (sup-normalized 9.9e-15,
G1); all 432 re-run backward paths keep their O1 classification (0 moved:
424 class-(a) stay pre-crossed, all 8 box-exit stay box-exit; margin
differences p95 1.9e-5, inside O1's own 2.8e-5 error estimate, G2); on the
24-path precision ladder every pre-crossing margin is stable to ~10 digits
(rel diff median 2.6e-11, max 1.8e-10 — four orders inside the 6-digit
gate) and the worst-margin path of the entire stratified subset is
certified at ladder grade, margin 0.620679003382508519636559682106572646
(G4). G3 is printed FAIL because ONE of 24 ladder paths misses the
pre-registered t = 0 endpoint bar 1e-8 by 1.70x (1.703e-8; the other 23
sit at <= 2.9e-10, median 4.3e-11) — a failure that LOCALIZES where the
flow's sensitivity lives (the free backward endpoints, which inherit the
full integration's amplification) and where it does not (the pre-crossing
margins, recorded at much earlier certified epochs). The kill-bin
pre-crossing statement of O1 is therefore re-established on the exact
field at precision-certification grade; the analytic T4-W5 proof remains
the terminal open.

## Provenance note (restart/resume)

The production run was interrupted by a container restart after phase
`dop` completed and phase `ladder` had begun; the phase-checkpoint design
worked as intended — the committed intermediates (`_g1.json`,
`_selection.npz`, `_dop.npz`) survived and the chain was relaunched from
phase `ladder` (the duplicated `== phase ladder ==` line in `p2_run.log`
is this resume). All phases downstream of the restart re-read only
committed files; no partial state was reused.

## What was re-run (the stratified subset)

Per family n_y in {+1, +0.75} (the -n_y mirrors are numerically identical
under z -> 1-z; re-verified here from o1_raw.npz: max margin deviation
3.5e-8 / 6.8e-9, classifications identical):

- **wall stratum (16/family):** ALL O1 BOX_EXIT paths (5 for |n_y| = 1,
  3 for 0.75) + the extreme-wall-row (z = 0.0025 / 0.9975) class-(a)
  paths with the earliest certified crossing epoch s_pc, filled to 16;
- **small-margin stratum (100/family):** the 100 smallest-margin class-(a)
  paths;
- **random stratum (100/family):** 100 random class-(a) paths
  (seed 20260718);

= 216 per family, **432 total**. Integration: scipy DOP853, rtol 1e-12,
atol 1e-14, dense output; margins refined by bounded local optimization on
the order-7 dense interpolant; s_pc by Brent root-finding of
X_x = d + 1e-3. The engine's velocity regularization (clips ±500, relative
density floor 1e-14, z-clip) is mirrored exactly and monitored: the clips
NEVER bind on any class-(a) path (max |a| = 7.05, max |b| = 7.23 vs 500).

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | exact-vs-grid field agreement at nodes <= 1e-12 rel (>= 1000 nodes x several times) | 9,000 nodes over 6 times (t = 0.5 ... 7.195), incl. 300 targeted well-conditioned nodes per time: sup-normalized \|df\|/max\|f\| = 9.9e-15 (all nodes); per-node relative f = 2.1e-13 and velocity components a, b = 1.8e-14 on well-conditioned nodes (3,711); conditioning-normalized a, b deviation on ALL nodes = 2.2e-15 — each <= 1e-12 | **PASS** |
| G2 | every re-run path keeps its O1 classification; margins agree within O1's own error estimate (pre-registered: p95 <= 5.6e-5 = 2x O1's p95, max <= 1.2e-2 = 2x O1's max; capped paths re-certify >= 2.0) | 432/432 classifications kept, 0 moved (424 class-(a), 8 box-exit). Margin differences over the 127 O1-uncapped class-(a) paths: median 6.4e-6, p95 1.88e-5 (vs O1's own refinement p95 2.81e-5), max 8.3e-5 (vs O1 max 6.0e-3). All 297 O1-capped paths re-certified: min exact-field margin 2.0010 >= cap 2.0. s_pc agreement: median 1.4e-4, max 2.2e-2 | **PASS** |
| G3 | precision ladder: float64-DOP853 vs 50-digit endpoints <= 1e-8; margins stable to >= 6 digits (rel <= 1e-6) | Margins clause EXCEEDED by ~4 orders: rel diff median 2.6e-11, max 1.8e-10 over 24 paths; s_pc abs diff median 3.7e-11, max 2.8e-10; mp step-halving margin shifts 6.6e-12 – 7.8e-12 on the 3 worst-margin paths; gmpy2-vs-mpmath field cross-check 1.8e-16. Endpoint clause MISSED by one path: diffs median 4.3e-11, but max 1.703e-8 > 1e-8 (path idx 415: n_y = 0.75, t0 = 6.565, z0 = 0.0225 — wall-adjacent, margin 18.19; the other 23 paths <= 2.9e-10) | **FAIL** |
| G4 | the worst-margin path certified at ladder grade, full record printed | Worst class-(a) margin of the full 432-path subset (n_y = +0.75, t0 = 6.405, z0 = 0.6825, O1 flat index 32136): margin_mp = **0.620679003382508519636559682106572646** (200-bit RK4, h ~ 5e-4, h/50 refinement); margin_f64 = 0.620679003403114 (rel agreement 3.3e-11 ~ 10 digits); s_pc = 2.8498644663664 (mp) vs 2.8498644663332 (f64); endpoint agreement 4.3e-11; mp step-halving margin shift 6.6e-12; min relative rho along the path 0.233 (deep inside the trusted region — no clip, floor, or clamp ever approached) | **PASS** |
| G5 | honest-scope paragraph: precision-certification, NOT formal interval arithmetic; what formal certification still requires; analytic proof stays the terminal open | Printed below and in p2_results.json | **PASS** |

## Key numbers

- **G3 FAIL diagnosis (the failed clause as a finding).** The FAIL is on
  the t = 0 ENDPOINT clause only, exceeded 1.70x by 1 of 24 ladder paths;
  the certification-relevant clause (margin stability) is passed with ~4
  orders of headroom. The physics: the backward endpoint at t = 0 sits at
  the END of the full backward integration (t0 ~ 6.4–7.2) and inherits the
  flow's cumulative sensitivity amplification — the one failing path is a
  wall-adjacent case (z0 = 0.0225) whose early backward segment runs
  through the strongest velocity gradients. The pre-crossing MARGINS are
  recorded at the much earlier certified crossing epochs (s_pc median 5.7
  for |n_y| = 1 / 4.9 for 0.75 in O1) and are provably insensitive:
  median rel diff 2.6e-11 across the ladder, step-halving shifts ~7e-12.
  The pre-registered 1e-8 endpoint bar was a spec choice slightly too
  tight for the longest wall-adjacent paths; per campaign rules the gate
  is printed FAIL, not relabeled. What the FAIL localizes: sensitivity
  lives in the free endpoints, not in the certified pre-crossing margins —
  the quantity the theorem consumes.
- **Margins on the exact field.** Min class-(a) margin 0.620679 (vs O1's
  grid-interpolated 0.6206889 — the 1e-5 shift is O1's own interpolation
  error, inside its estimate); median class-(a) margin 4.11; the 297
  O1-capped paths, uncapped here, reach margins up to ~18. The margin
  histogram (O1 vs P2 overlay), the G2 agreement histogram against O1's
  error estimate, and the full ladder-agreement panel are in the figure.
- **Ladder composition.** 24 paths = the 12 globally smallest margins +
  6 spread over family +1 by s_pc quantiles + 6 spread over family +0.75
  by margin quantiles (6 from family +1, 18 from +0.75); + step-halving
  (h -> 2.5e-4) re-runs of the 3 worst-margin paths.
- **Monitors.** Velocity clips (±500) never bind on any class-(a) path
  (max |a| = 7.05, |b| = 7.23). Three class-(a) paths touch the relative
  rho floor (1e-10) on the monitored subsample — but only AFTER their
  certified crossing (s_floor = 0.83–2.21 vs s_pc = 3.96–4.59; backward
  time runs downward, so the certified segment [s_pc, t0] is floor-free
  for all 424 class-(a) paths, preserving O1's precedence semantics).
- **Backend equivalence.** The ladder runs in 200-bit (~60-decimal-digit)
  GMP/MPFR arithmetic via gmpy2 (the backend mpmath uses); a direct
  cross-check of one field evaluation against mpmath at dps = 60 agrees
  to 1.8e-16 (float64-representation-limited, as recorded).
- **Runtime.** dop: 1,948 s wall on 4 workers (7,666 s CPU-sum; mean
  187,658 RHS evaluations/path — each an exact 2048-mode sum, no
  interpolation anywhere); ladder: 2,929 s wall (27 mp integrations of
  ~12,800–14,400 RK4 steps each). Total production ~81 min.

## Method

The state is exactly factorized, phi(x, z, t) = f(x, t) sin(pi z)
e^{-i pi^2 t/2}, and f(., t) is the free evolution of the Nx = 2048 grid
datum on the periodic box [-15, 45]: f(x, t) = sum_k C_k
e^{i kappa_k (x - xmin)} e^{-i kappa_k^2 t/2} with C_k = F0[k]/Nx,
kappa_k = (pi/30) m_k, m_k in [-1024, 1023] — the FFT coefficients ARE the
exact spectral coefficients of the band-limited field the engine evolves,
which the grid engine reproduces exactly at nodes (G1 verifies this to
rounding). f, f' at arbitrary (x, t) by direct mode summation;
v_x = Im(f'/f) - n_y pi cot(pi z), v_z = n_y Re(f'/f); engine
regularization (clips ±500, relative density floor 1e-14, z-clip 1e-9)
mirrored exactly and monitored. Float64 leg: DOP853 (rtol 1e-12,
atol 1e-14), dense output, box-rim termination events; margin by bounded
minimization on the dense interpolant; s_pc by brentq. Multiprecision
ladder leg: fixed-step classical RK4, h = t0/ceil(t0/5e-4) (~5e-4), in
200-bit gmpy2/MPFR arithmetic on the exact mode sum, with the per-mode
phase e^{-i kappa_k^2 t/2} advanced by precomputed half-step quadratic
phase recurrences (u^{(m+1)^2} = u^{m^2} u^{2m+1}; one complex exp per
phase array, no per-step trigonometry) and f, f' by two power ladders in
FFT order; margins refined by h/50 local re-integration + 3-point
parabola; s_pc by h/50 bracket + linear interpolation; step-halving
(h/2) on the 3 worst-margin paths bounds the mp-side truncation
independently of float64. Selection from o1_raw.npz (O1's reference-run
classification), seed 20260718, mirror families re-verified before
selection. hbar = m = 1.

## Honest scope (G5)

This is PRECISION-CERTIFICATION, NOT FORMAL PROOF: the field is exact (a
finite trigonometric sum evaluated to floating-point rounding, G1) and the
integration is converged (adaptive rtol 1e-12 cross-validated against
50-digit fixed-step arithmetic, G3's margin clause), but no step carries a
rigorous enclosure. Formal certification would still require validated ODE
integration with interval / Taylor-model enclosures (CAPD/COSY-class) of
the backward flow over the exact field, with outward-rounded mode sums, so
that each pre-crossing margin becomes a machine-checked inequality. The
analytic T4-W5 proof (no computation at all) remains the terminal open
item. What this pilot adds to O1: the dominant O1 error source (grid
interpolation) is eliminated, every re-run classification survives on the
exact field, and the worst margin in the stratified subset is now known to
~10 significant digits with an independent 60-digit confirmation.
Within-model; nothing here bears on nature.

## Caveats

- G3 is a genuine FAIL of the pre-registered composite gate and is
  reported as such; the diagnosis above (endpoint clause, 1 path, 1.70x)
  is printed alongside, in the campaign's failed-clause-as-finding
  tradition. No gate was relabeled and no bar was moved after the fact.
- The subset is stratified, not exhaustive: 432 of O1's 64,000 kill-bin
  paths (all 16 box-exit/wall-adjacent members, the 100 smallest margins
  and 100 random class-(a) per family). The strata are chosen to bound
  the worst cases; the remaining ~63,500 paths carry margins >= the
  re-certified minimum by O1's classification, but were not individually
  re-run here.
- Mirror families (-1, -0.75) were not re-integrated; their numerical
  identity with the + families was re-verified from o1_raw.npz at the
  3.5e-8 level (and exactly, classification-wise) before selection.
- The 8 box-exit paths remain flux-bounded bookkeeping (as in O1), not
  certified pre-crossings; they keep their classification on the exact
  field, and their O1 flux bound (2.8e-7-class) is unchanged by this
  pilot.
- The mp ladder integrates the same clipped/floored velocity law as the
  engine; on certified segments no clip or floor ever binds (monitored),
  so there the law coincides with the exact j/rho.
- The dop-phase monitors sample every 4th dense-output node (~1e-3
  spacing); the three post-certification floor touches are localized to
  that sampling. Certified segments are unaffected (floor events all at
  s < s_pc - 1.7).
- Endpoint sensitivity grows with backward-integration length; the 1e-8
  endpoint bar would need ~2x relaxation (or endpoint-free phrasing) to
  be attainable for the longest wall-adjacent paths at float64. This is
  a statement about the spec, recorded for any successor workstream, not
  a retroactive amendment.
- Configuration-specific: d = 1, k0 = 2, box [-15, 45], Nx = 2048
  band-limit inherited from O1's reference grid — the exact field IS the
  band-limited field that grid defines, not a continuum-limit field.
