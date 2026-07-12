# Tier-4A — symmetry-free 3-D Cartesian fixed-L Routhian descent at eps = 0.05

Reproduce with `python3 field3d_solve.py` (deterministic — the only RNG uses
fixed seeds 424242/11/12; 4461 s wall; writes `field3d_results.json`,
`field3d_decay.png`, `field3d_run.log`; `--smoke` runs a small-N pipeline
test with separate outputs).  This run adjudicates the Step-3 fixed-L
claims in fully 3-D Cartesian space with NO symmetry assumption: no
axisymmetry, no mode basis, no profile ansatz — the field is a free
4-component unit quaternion on an N = 96 cell-centred grid.  Construction
and measured numbers only; adjudication is the coordinator's layer.

## 1. Construction

- Grid: cell-centred 96^3, box half-width 4.5 = 2.53 R* (R* = 2^(5/6));
  h = 0.09375; fixed-vacuum ghost cells.  Box-sanity rerun at 1.3x box
  (half-width 5.85, same N, h = 0.12188).
- Field q = (q0,q1,q2,q3), |q| = 1 enforced by renormalization after every
  update.  Sector conventions exactly as Tier-2 (`axi_RESULTS.md`):
  E2, E4 with 1/4pi, E6 = pi^3 INT b^2 via the 4x4 determinant density
  det[q, Dq] (deg = INT b, sign fixed on the hedgehog), E0 = (1/4pi)
  INT (1-q0), I = INT 2(q1^2+q2^2); E_static = t(E2+E4)+E6+E0 at the
  frozen Step-1 dial t = 0.0082764349; R(q; L) = E_static + L^2/(2I).
- 1-D hedgehog profile re-solved in-script with the Step-1 method (graded
  midpoint grid, flow + damped Newton, t frozen): E2 = 12.117459,
  E4 = 6.191126, E6 = 1.523936, E0 = 1.507588, eps-ratio 0.049985,
  gmax = 4e-13 — the Step-1 numbers to 7 s.f.  Boundary tail of the 3-D
  hedgehog: |sin f|_max/pi = 5.5e-11 on the box faces (gate 1e-6, PASS).
- Two discretizations, both gated: `central4` (4th-order central
  differences, the MEASUREMENT engine) and `corner` (compact midpoint
  scheme — derivative sectors at the (N+1)^3 cell corners, D_i = edge
  difference of face averages, the 3-D analogue of Step-1's midpoint rule
  — the DESCENT objective).  Analytic gradients; FD-vs-analytic gradient
  gate over both schemes, E_static and Routhian objectives: worst rel
  3.2e-6 (tol 1e-5, PASS).
- Minimizer: arrested Newton flow (velocity on q, acceleration = minus the
  tangent-projected functional gradient; step reverted and velocity zeroed
  on any objective increase; adaptive dt), objective monotone by
  construction.

### Regularization protocol (numerical obstacles, measured before being walled)

At eps = 0.05 the functional is 95% BPS and the BPS sector (E6+E0) does
not control either the degree or the shape on a lattice.  Two failure
channels were measured in unregularized probe runs and then excluded:

1. **Charge leak.** The unconstrained flow unwinds the soliton: degree
   0.99991 -> 0.9994 -> 0.95 -> 0.32 with exponentially growing rate; the
   continuum anti-unwinding barrier of the eps-suppressed (2+4) sectors is
   O(t) and sits below the discretization error.  (Step 1 pins f(0) = pi
   as a Dirichlet condition; the Step-2/3 ansaetze build it in; a free 3-D
   lattice has no such pin.)  Excluded by a ONE-SIDED degree anchor
   E_pen = (MU/2) min(0, deg-(deg_ref-0.005))^2, MU = 5000, deg_ref = the
   engine's own hedgehog degree.  Its gradient acts only through the
   topological density; the far-field tilt channel (q1,q2 at large r,
   det ~ 0) is untouched.
2. **Quadrature-cheat drift.** With the degree held, the flow still lowers
   the DISCRETE (6+0) sectors below the continuum Bogomolny floor
   E6+E0 >= (32 sqrt2/15) B = 3.016989 B by roughening the core (measured:
   central4-objective static descent reached E_static = 3.0094 <
   3.01699 deg — impossible in continuum; the corner scheme drifts too,
   more slowly).  Excluded by a ONE-SIDED Bogomolny-floor wall
   E_fpen = (MU_F/2) min(0, fgap - fgap_ref)^2, MU_F = 400,
   fgap = E6+E0-3.016989 deg/deg_ref, fgap_ref = the engine's own
   hedgehog fgap minus 0.01.  Every continuum configuration satisfies the
   floor, so the wall only excludes lattice-artifact territory.
3. **Halo seed.** The descent perturbation adds, to the fixed-seed random
   bumps (amplitude 0.02), an explicit small tilt-halo shell (eta = 0.05
   in q1 at r ~ 3.2) — an accelerator applied IDENTICALLY to the over-spun
   and control runs, so stability/instability of the seeded direction is
   tested two-sidedly.

All reported energies and Routhians are the BARE functionals; the anchors'
values are logged at every record (final |E_pen|, |E_fpen| <= 1.1e-3 in
every run).  The verdict-carrying instruments — I, kappa = L/I, halo
fraction, degree, E0 — contain no derivatives, are quadrature-exact on the
cells, and are identical in both schemes.

## 2. Section G — gates (hedgehog vs Step-1 references)

| engine | N | E2 rel | E4 rel | E6 rel | E0 rel | I rel | degree | virial rel |
|---|---|---|---|---|---|---|---|---|
| central4 | 48 | -2.2e-2 | -2.2e-3 | +3.9e-3 | +2.4e-6 | -7.7e-6 | 0.998575 | -6.3e-3 |
| central4 | 64 | -9.7e-3 | -9.3e-4 | +1.7e-3 | -2.8e-6 | -8.4e-6 | 0.999541 | -2.7e-3 |
| **central4** | **96** | **-2.7e-3** | **-2.7e-4** | **+4.1e-4** | **+1.0e-8** | **-5.1e-8** | **0.999911** | **-6.6e-4** |
| corner | 64 | -4.1e-2 | -3.2e-2 | -4.9e-2 | -2.8e-6 | -8.4e-6 | 0.968958 | +7.2e-2 |
| corner | 96 | -1.9e-2 | -1.4e-2 | -2.2e-2 | +1.0e-8 | -5.1e-8 | 0.986012 | +3.2e-2 |

- central4 at N = 96: worst sector 2.7e-3 (gate 1e-2 **PASS**), |deg-1| =
  8.9e-5 (gate 5e-3 **PASS**); N64/N96 error ratios 3.4-4.1 (better than
  the O(h^2) value 2.25).
- corner at N = 96: worst sector 2.2e-2 (gate 1e-2 **FAIL** — reported per
  protocol and proceeded with the caveat), |deg-1| = 1.4e-2 (**FAIL**);
  N64/N96 ratios 2.11/2.18/2.22 = clean O(h^2).  The corner engine is the
  descent objective because it is robust (compact stencils, no sub-grid
  null modes), not because it is the more accurate quadrature; every
  descent conclusion below is a SAME-FUNCTIONAL comparison (descent
  endpoint vs axisymmetric reference evaluated in the same corner
  functional), which cancels the scheme offset.
- E0 and I are quadrature-exact (1e-8 level) in both schemes at both N.

## 3. Section A — on-grid axisymmetric family-A references (min over dilation d of R(hedgehog(x/d); L))

| box, L | d* | R_A (corner) | R_A (central4) | I | kappa |
|---|---|---|---|---|---|
| main, L = 8.4979 | 1.12569 | 4.492979 | 4.514038 | 31.764 | 0.26753 |
| big, L = 8.4979 | 1.12379 | 4.478743 | 4.514401 | 31.604 | 0.26888 |
| main, L = 4.0 | 1.03089 | 3.492102 | 3.523918 | 24.397 | 0.16396 |

The central4 on-grid reference reproduces the Step-2 spectral solA
(R = 4.513718, d = 1.12844, kappa = 0.265578) to 7e-5 relative in R — the
two tiers' engines agree on the axisymmetric stationary point.  (The task
sheet quoted "4.615-ish" for this value; the recorded `axi_results.json`
solA Routhian is 4.513718, and that is the number used and logged here.)

## 4. Section S — static minimization (hedgehog + fixed-seed perturbation)

Start: corner-hedgehog E_static = 3.146624 (deg 0.986012); perturbed
+4.6e-4.  ANF, 1200 iterations (cap): E_static 3.14709 -> 3.10316
(-1.4%), gnorm/gnorm0 = 0.54 — the 1e-6 gradient-reduction target was
NOT reached; the run ends in the iteration cap with the degree pinned at
the anchor band edge (0.98043) and the floor gap pinned at the wall
(-0.0314).  Measured behaviour: the flow does NOT return to the exact
hedgehog and does not converge to an isolated discrete minimum; it
wanders down a quasi-flat internal rearrangement valley of the 95%-BPS
functional (E0 1.5077 -> 1.4032 and E6 1.4906 -> 1.5653 with E6+E0 moving
only -0.014, I 22.264 -> 22.034, centroid drift < 1.2e-3), staying in the
hedgehog's neighbourhood in every quadrature-exact metric while the
corner objective drifts by the wall-bounded -1.4%.  Virial residual
t E2 - t E4 - 3E6 + 3E0 = -0.455 (rel -0.147) at the endpoint, against
-2.1e-3 (rel -6.6e-4) for the clean hedgehog on the central4 engine and
6e-7 in Step 1 — the endpoint is not a continuum-quality stationary
point.  central4 re-measurement of the endpoint gives E_static = 6.83
(2.2x the corner value): the configuration carries grid-scale structure.
This is the measured 3-D lattice landscape of the near-BPS functional at
this resolution: below the ~1.5% level the static problem is
quadrature-degenerate, and "the static minimum" is defined only up to
that mush.  (This is the documented caveat inherited by everything
below; the descent conclusions are therefore stated only through the
quadrature-exact instruments and same-functional differentials, which
are 10-40x larger than the mush.)

## 5. Sections M, C, B — the fixed-L Routhian descents

All three start from the static endpoint (big box: its own perturbed
hedgehog) + the same fixed-seed perturbation + the same halo seed.
Corner objective; every quantity below is the bare functional.

### L = 8.4979 (Step-2 solA charge; kappa_start = 0.374, 1.9x threshold)

| it | R | E_static | E0 | E6 | I | kappa | halo frac | deg |
|---|---|---|---|---|---|---|---|---|
| 0 | 4.7062 | 3.1171 | 1.4169 | 1.5652 | 22.722 | 0.3740 | 0.028 | 0.9803 |
| 400 | 4.2582 | 3.1293 | 1.5638 | 1.4047 | 31.984 | 0.2657 | 0.046 | 0.9804 |
| 800 | 4.1235 | 3.1423 | 1.6508 | 1.3175 | 36.797 | 0.2309 | 0.060 | 0.9804 |
| 1600 | 4.0047 | 3.1541 | 1.7559 | 1.2126 | 42.449 | 0.2002 | 0.081 | 0.9804 |
| 2400 (cap) | 3.9330 | 3.1591 | 1.8346 | 1.1338 | 46.655 | 0.1821 | 0.100 | 0.9804 |

Measured facts:

- **R descended monotonically to 0.560 BELOW the axisymmetric family-A
  stationary value measured in the same functional on the same grid**
  (3.9330 vs 4.4930), crossing it between iterations 100 and 125, and
  was still descending at the iteration cap (status `iter_cap`, slope
  ~ -8e-5/it).  For scale: the wall-bounded static mush is -0.044 and the
  control's excess descent is -0.145; the main channel is 4-13x those.
- **kappa fell from 0.374 through the axisymmetric stationary value
  0.26558 (at it ~ 400) and through 1/sqrt(8 pi) = 0.19947 (at it ~ 1600)
  to 0.1821 at the cap, still falling.**  I rose from 22.72 through the
  soft-halo saturation value sqrt(8 pi) L = 42.60 (at it ~ 1620) to 46.65.
- Inertia migrated outward monotonically: halo fraction (I outside
  1.5 R*) 0.028 -> 0.100 and rising.  Degree pinned at the anchor
  (0.9804) throughout; |E_pen|, |E_fpen| <= 1.1e-3.
- Sector accounting at the endpoint: E6 fell and E0 rose to almost
  exactly the solA (dilated-core) proportions (corner 1.134/1.835 vs
  continuum solA 1.060/2.166 with the -2% corner offset), i.e. the core
  dilated as in family A; but I = 46.65 = 1.46x solA's 31.998 with LESS
  E0 than solA — the surplus inertia is carried by internal-direction
  tilt (sin^2 Theta -> 1), the Step-3 tilt channel.  Of the total
  dI = +23.9, 16 pi dE0 = +21.0 (88%) is priced at exactly the
  soft-condensate rate 1/(16 pi) that defines the threshold; the
  remainder is E0-free tilt.  The tilt surplus is also why kappa passes
  below the pure-halo saturation value: the halo reservoir sets
  I -> sqrt(8 pi) L while the bounded tilt dial adds on top at fixed L.

### L = 4.0 control (kappa_start = 0.176, 0.88x threshold)

| it | R | E_static | E0 | E6 | I | kappa | halo frac |
|---|---|---|---|---|---|---|---|
| 0 | 3.4692 | 3.1171 | 1.4169 | 1.5652 | 22.722 | 0.1760 | 0.028 |
| 800 | 3.3713 | 3.1092 | 1.5336 | 1.4347 | 30.520 | 0.1311 | 0.041 |
| 1400 (cap) | 3.3473 | 3.1094 | 1.5836 | 1.3847 | 33.622 | 0.1190 | 0.050 |

Measured facts: no far-field condensate forms.  E_static is FLAT
(-0.008 net) — the control never buys inertia at the E0-paying condensate
rate; its I growth saturates at the BOUNDED internal tilt ceiling
(I -> 33.62 ~ (3/2) I_hedgehog = 33.40, the maximal g-dial value), with
halo fraction <= 0.050 (vs 0.100/0.163 in the over-spun runs) — the
r > 1.5 R* fraction it does gain tracks the d* = 1.03 dilation and skirt,
not a growing far field.  R descends 0.145 below its axisymmetric
reference — the bounded tilt gain plus the static-mush drift, an order
of magnitude short of the over-spun runs' channel and self-limiting in I.
(In the continuum the tilt directions carry small positive stiffness —
axi3 Section 6 measured d2E > 0 for every tilt mode — so the full extent
of the control's ride up the tilt dial is lattice-mush-assisted; its
BOUNDEDNESS, and the absence of the E0-paying far-field channel below
threshold, are the control facts this run establishes.)

### Box-size sanity (L = 8.4979, box 1.3x, same N, h = 0.12188)

| it | R | I | kappa | halo frac |
|---|---|---|---|---|
| 0 | 4.7098 | 22.950 | 0.3703 | 0.028 |
| 800 | 4.0158 | 42.297 | 0.2009 | 0.130 |
| 1600 (cap) | 3.9209 | 48.218 | 0.1762 | 0.163 |

At MATCHED iteration count (it = 1600) the 1.3x box sits 0.070 DEEPER
relative to its own same-grid axisymmetric reference (dR = -0.558 vs
-0.488 in the main box), with more inertia (48.2 vs 42.4), lower kappa
(0.176 vs 0.200) and a substantially larger halo fraction (0.163 vs
0.081); both runs are still descending at their caps.  The descent
deepens and the inertia migrates further outward when the box grows —
the far-field growth is room-limited, not scheme-limited.

## 6. Section K — clock invariant on the control solution

On the L = 4.0 endpoint (E_static = 3.10941, I = 33.6221):

- E_rot/E_tot at L = 4.0: **0.071083** — NOT 1/4: the invariant does not
  hold at arbitrary spin.
- Clock L^2 = (2/3) I E_static solved by bisection at frozen fields
  (identical to the closed form): **L_clock = 8.34845**; at L_clock,
  E_rot/E_tot = **0.250000** — the 1/4 invariant holds exactly and only
  when the clock is imposed (algebraic identity, verified numerically).
- **kappa_ours(L_clock) = 0.248302 = 1.245x the halo threshold
  1/sqrt(8 pi)** — the clock charge places the solution IN the over-spun
  regime.  Cross-check on the clean on-grid axisymmetric control
  reference (d* = 1.0309, E_static = 3.16419, I = 24.3966, i.e. reading
  the control's continuum-expected endpoint instead of the
  tilt-saturated lattice one): L_clock = 7.1738, kappa = 0.29405 =
  1.474x threshold.  Step 2's solA itself sits at 1.331x.  Under every
  reading measured here, imposing the corpus's clock on a B = 1 solution
  of this functional lands the charge above the tilt-halo threshold.
  (The central4 re-measure of the rough control endpoint gives
  E_static = 22.9 and is not a usable clock input; it is listed in the
  JSON as the roughness diagnostic it is.)

## 7. Convergence, caveats, determinism

- Statuses: every ANF run ended at its iteration cap (static 1200, main
  2400, control 1400, big box 1600), monotone in the objective
  throughout; gnorm/gnorm0 at the caps: 0.54 / 0.11 / 0.05 / 0.29.  None
  is a converged stationary point; the over-spun descents are minimizing
  SEQUENCES still in motion — every verdict-relevant statement above is a
  trend or a crossed reference value, not an endpoint claim.
- The descent endpoints carry grid-scale structure: central4
  re-measurement of the corner endpoints gives 6.8 (static), 23.2
  (control), 77.0 / 74.7 (over-spun) — recorded in the JSON as
  `final_central4`.  Conclusions rest exclusively on (a) quadrature-exact
  instruments (I, kappa, E0, halo fraction, degree) and (b)
  same-functional corner-vs-corner differentials.
- The corner engine misses the 1% sector gate at N = 96 (2.2% worst,
  O(h^2) trend); the central4 measurement engine passes all gates.  The
  main descent's crossing margins (-0.560 in R, -0.084 in kappa below
  solA at threshold-crossing time) are 10-25x the corner scheme's
  absolute offsets on smooth configurations.
- Regularization footprint at every reported endpoint: |E_pen| +
  |E_fpen| <= 2.1e-3, degree within 0.006 of the engine's hedgehog
  value, floor gap within 0.012 of the hedgehog's own.  Runs are
  deterministic (fixed seeds; reruns reproduce the log except timings).
- Runtime 4461 s at N = 96 on 4 cores.

## 8. Files

`field3d_solve.py` (engine + all sections, self-contained),
`field3d_results.json` (gates, references, full instrumented time series
of all four flows, clock block), `field3d_decay.png` (R, kappa, halo
fraction vs iteration for both L values with the big-box overlay and the
solA / threshold reference lines), `field3d_run.log` (printed tables).
Smoke-test artifacts (`*_smoke.*`) are pipeline checks at N = 32, not
physics.

## 9. Epistemic notice

Everything above is a within-model computation on the corpus's defined
functional, sector conventions and closure definitions (Tier-2 unit map,
zero calibrated constants), executed on a finite lattice with the two
documented one-sided regularizations and an explicit seeded perturbation.
The near-BPS lattice pathologies (charge leak, quadrature mush,
grid-scale roughening of minimizing sequences) are properties of the
discretized problem and are reported as such; nothing here is a statement
about nature, about the corpus's own numerical methods, or about any
reading of its text beyond the compositional definitions in force since
Step 1.
