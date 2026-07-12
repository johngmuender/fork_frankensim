# Tier-4C — direction-locked (hedgehog-locked) isorotating closure

Reproduce with `python3 axi4_locked_solve.py` (deterministic, no RNG; 2527.3 s
wall at full resolution; writes `axi4_locked_results.json`,
`axi4_locked_landscape.png`, `axi4_locked_epsscan.png`; every printed table in
`axi4_locked_run.log`). Construction and measured numbers only; adjudication
is the coordinator's layer.

## 1. The locked problem, basis, quadrature, cap protocol

Hypothesis under test (axi3_ADJUDICATION.md Section 3): the corpus's frozen
`<r1>` benchmark code plausibly direction-locks the field — hedgehog
direction fixed, `Theta(r,theta) == theta` exactly, only the profile
F(r,theta) free — removing Step 3's tilt-halo channel, so that the benchmark
could be a genuine constrained minimum. This step builds that locked problem
and measures it.

```
q = ( cos F,  sin F sin theta,  0,  sin F cos theta )          (the lock)
F(r,theta) = f_base( r * exp( SUM_kl a_kl B_k(r) C_l(mu) ) )
```

- `B_k`: K = 4 (default) / K = 5 (basis ladder) normalised cubic-B-spline
  bumps spanning the base support; `C_l`: even Legendre {1, P2, P4} —
  12 / 15 parameters, richer than axi3's contour group (K = 3, {1, P2}).
- Overall dilation d carried analytically (`E2 ~ d, E4 ~ 1/d, E6 ~ 1/d^3,
  E0 ~ d^3, I ~ d^3`; golden in ln d inside every objective evaluation).
- Quadrature engine, sector constants, unit map (compacton-anchored,
  `c_paper = c_ours/sqrt(2 pi^3)`, `kappa_paper = 2 sqrt(pi) kappa_ours`,
  zero calibrated constants) and clock: REUSED unchanged from
  `axi_solve.py` / `axi3_solve.py`.
- Closure: minimise `R(a; L) = E_static + L^2/(2I)` at fixed L
  (deterministic Nelder-Mead, warm-started outers), update L by the
  closed-form clock fixed point `L^2 = (2/3) I E_static`, iterate; after
  convergence/termination the (d, L) fixed point is re-solved exactly at
  frozen a*, so the reported CLOCK residuals are `|2L^2/I - R|/R <= 7.9e-9`
  at every reported solution (most <= 6e-9; the 1e-8 residual target is met
  in that sense at every row). The outer |dL/L| record is in Section 8 —
  for cap-limited runs it stalls at 1e-4..6e-4, the boundary-creep signal.
- Caps are soft quartic walls; ACCEPTANCE is convergence, not caps: a level
  is accepted only if no optimum sits within 80% of any cap, otherwise every
  cap is widened and the solve repeats (ladder):

| cap level | domain | shrink cap (max A) | per-coeff box | stretch cap | N_opt |
|---|---|---|---|---|---|
| L0 | 1.5 r_sup | 0.45 | 0.80 | 0.3755 | 200 |
| L1 | 2.0 r_sup | 0.80 | 1.40 | 0.6631 | 240 |
| L2 | 2.7 r_sup | 1.40 | 2.40 | 0.9633 | 288 |

Polish grid N = 400; final evaluation and gates N = 720; coarse rerun of
every final answer N = 360. Seeds: hedgehog, spheroidal-shell projections
onto {1, P2, P4} at lam = 0.7 / 0.45, outer-bump skirt seed; scan and ladder
runs warm-started.

## 2. Validation gates (all PASS, measured at N = 720)

| gate | what it checks | measured |
|---|---|---|
| G1 | p = 0 sectors vs step-1 radial integrals | E2 -3.4e-4, E4 -1.7e-4, E6 -1.1e-4, E0 +6.3e-6; degree 0.99991 |
| G2 | p = 0 inertia vs radial formula | +2.0e-7 |
| G3 | compacton p = 0 vs exact | E6 +1.9e-5, E0 +8.9e-7, I -2.4e-7; degree 1.00003 |
| G4L | locked deformed config (K = 4, P4 active) vs independent spherical-coordinate quadrature | E0 +7.8e-7, I +2.9e-7; degree 1.00004 |
| G5 | p = 0 closure at eps = 0.05 vs step-2 solA | c_paper 2.15824 vs 2.15825 (rel -7.5e-6) |

## 3. Task A — the locked halo threshold (derived independently, then measured)

Independent derivation (this step's algebra). For the locked field, add a
far halo `F = eta chi(r) s(theta)` disjoint from the core. Exactly at
O(eta^2), in these conventions (`E0 = (1/4pi) INT (1-cos F) dV`,
`I = INT 2 sin^2 F sin^2 theta dV`; the t-weighted gradient costs per unit
halo norm scale as 1/r_c^2 and vanish in the far-field limit; dE6 = O(eta^4)):

```
dE0 = (1/8pi) INT eta^2 chi^2 s^2 dV
dI  = 2 <sin^2 theta>_w INT eta^2 chi^2 s^2 dV        (w = s^2 weight)
dR  = [INT eta^2 chi^2 s^2 dV] ( 1/(8pi) - kappa_ours^2 <sin^2 theta>_w )

=>  kappa_crit_ours(s)  = 1/sqrt(8 pi <sin^2 theta>_w)
    kappa_crit_paper(s) = 1/sqrt(2 <sin^2 theta>_w)
```

For the angle-uniform halo s = 1 (`<sin^2 theta> = 2/3`) this reproduces the
task sheet's expected locked threshold, CONFIRMED:
`kappa_crit_ours = sqrt(3/(16 pi)) = 0.24430`, `kappa_crit_paper = sqrt(3)/2
= 0.86603`.

**Flagged discrepancy vs the single-threshold expectation.** The lock
freezes Theta but not the angular shape of F itself, so the locked problem
retains a one-parameter FAMILY of halo onsets labelled by the halo's angular
profile s. At t = 0 the angular concentration is exactly free at O(eta^2)
(E0 and I are pointwise densities — no gradient enters), and at t > 0 the
angular-gradient cost per unit norm falls off as 1/r_c^2. The infimum over s
is the equatorial-ring limit `s^2 -> delta(theta - pi/2)`, which reproduces
the Step-3 tilt threshold `1/sqrt(8 pi)` (paper `1/sqrt2`) exactly: the lock
does not remove that channel, it only requires equator-concentrated far
halos to reach it. Three members of the family are closed-form identities:

- s = 1: `kappa_paper = sqrt(3)/2 = 0.86603` (exact algebra);
- s = sin^2(theta): `<sin^2 theta>_w = 6/7` exactly, so `kappa_paper =
  sqrt(7/12) = 0.76376` EXACTLY — the same closed form as the corpus's
  deep-BPS kappa_0;
- s^2 = |cos(theta)|: `<sin^2 theta>_w = 1/2` exactly, so `kappa_paper =
  1.00000` EXACTLY.

Engine-level probes (explicit shell at r_c = 2.1 r_sup, width 0.45 r_sup,
eta = 0.02 / 0.04, ALL sectors measured, N = 480):

| angular profile | <sin^2 th>_w | derived paper | measured t=0 | measured eps=0.05 | dI ratio eta->2eta |
|---|---|---|---|---|---|
| s = |cos| (polar) | 2/5 | 1.11803 | 1.11805 | 1.12132 | 4.000 |
| s^2 ~ |cos| (smoothed) | 0.512129 | 0.98809 | 0.98810 | 0.99089 | 3.999 |
| s^2 = |cos| (exact limit) | 1/2 | **1.00000** | — | — | — |
| s = 1 (uniform) | 2/3 | **0.86603** | 0.86606 | 0.86844 | 3.999 |
| s = sin | 4/5 | 0.79057 | 0.79059 | 0.79280 | 3.999 |
| s = sin^2 | 6/7 | **0.76376 = sqrt(7/12)** | 0.76378 | 0.76595 | 3.999 |
| s = sin^4 | 10/11 | 0.74162 | 0.74164 | 0.74381 | 3.999 |
| ring limit (s^2 -> delta) | 1 | 0.70711 | — | — | — |

(The smoothed polar row uses `s = (cos^2 + 0.01)^{1/4}`; its own
`<sin^2 th>_w = 0.512129` is what its derived value reflects. The tilt-halo
channel that the lock EXCLUDES sits at 0.70711 ideal / 0.76659 as measured
by axi3's partial-tilt shell.)

- Quadratic scaling: dI(2 eta)/dI(eta) = 3.999-4.000 on every row, both
  profiles — the channel is quadratic as derived.
- Radius scan (position independence): at t = 0, moving the shell from
  r_c = 1.5 r_sup (w = 0.22) to 3.2 r_sup changes the measured threshold by
  < 5e-5 (0.24431 at both radii; sin^2 row 0.21546 vs ideal 0.21545) —
  position-independent, as the pointwise-density derivation requires. At
  eps = 0.05 the measured values sit 0.3-1.1% above ideal (0.24700 at
  1.5 r_sup, 0.24495 at 3.2 r_sup vs 0.24430) — the finite-t, finite-r_c
  gradient cost, decaying with radius as derived.
- Edge shells (mid-field overlap probes, Task D): a shell straddling the
  profile (r_c = R* for the compacton; r_c = 0.40 r_sup for the eps profile,
  where f ~ 0.1) has dI LINEAR in eta — dI(2 eta)/dI(eta) = 2.059 (t = 0) /
  2.077 (eps) with dE quadratic (ratios 3.955 / 3.992) — i.e. it is a
  first-order profile-relaxation direction already inside the minimised
  family, not a distinct onset channel; no additional threshold exists at
  the compacton edge.

## 4. Task A — placement of the corpus tuples (measured, not assumed)

`kappa_ours = kappa_paper / (2 sqrt(pi))`; x columns are ratios to the
uniform threshold 0.24430 and the ring limit 0.19947; `req <s^2>` is the
minimal halo concentration `<sin^2 theta>_w = 1/(8 pi kappa_ours^2)` that
opens a descent direction against the tuple.

| tuple | kappa_paper | kappa_ours | x uniform | x ring | vs uniform halo | vs ring halo | req <s^2> |
|---|---|---|---|---|---|---|---|
| scaling rung | 0.9354 | 0.26387 | 1.080 | 1.323 | UNSTABLE | unstable | 0.5714 (< 2/3) |
| `<r1>` benchmark | 0.8020 | 0.22624 | 0.926 | 1.134 | stable | UNSTABLE | 0.7774 |
| deep-BPS endpoint | 0.7638 | 0.21546 | 0.882 | 1.080 | stable | UNSTABLE | 0.85706 = 6/7 to 1e-4 |

Measured placement: (i) the scaling rung remains halo-unstable even under
the lock (1.080x the uniform threshold — even angle-uniform halos open);
(ii) the benchmark and the deep-BPS endpoint sit BELOW the uniform locked
threshold and ABOVE the ring limit — the task-sheet placements (0.2262
below, 0.2639 above, 0.2154 below vs 0.24430) are confirmed as stated, with
the addition that descent against the two "below" tuples reopens for halos
concentrated beyond `<sin^2 theta>_w` = 0.7774 / 0.857 respectively;
(iii) the deep-BPS kappa equals the s = sin^2 member of the locked onset
family exactly (both are sqrt(7/12); its required concentration is that
profile's 6/7).

## 5. Task B — locked closure at eps = 0.05 (t = 0.0082764349): cap ladder

Direct fixed-L + clock closure over the locked family. NO cap level yields
an interior optimum at eps = 0.05 — every level terminates on a wall with c
still rising and kappa still falling toward the ring values; the ladder (and
the K = 5 basis enlargement) is the measured no-interior-optimum record:

| config | c_paper | kappa_paper | g* | lam_eff | V* | V_E0 | walls at >= 80% | interior |
|---|---|---|---|---|---|---|---|---|
| K4/L0 | 3.05361 | 0.75248 | 0.6160 | 0.5579 | 4.129 | 1.643 | amax 100%, stretch 95% | NO |
| K4/L1 | 3.23535 | 0.72827 | 0.4183 | 0.3791 | 6.657 | 1.685 | stretch 97% | NO |
| K4/L2 (accepted widest) | **3.26271** | **0.72483** | 0.3968 | 0.3405 | 7.112 | 1.691 | amax 94% | NO |
| K5/L2 (15 params) | **3.28533** | **0.72230** | 0.3113 | 0.3074 | 9.160 | 1.697 | amax 100% | NO |
| ring-saturated locked closure (well-posed limit) | 3.38547 | 0.70711 exact | g_core 0.991 | — | — | — | interior (34%) | — |

Every row: E_rot/E = 0.250000, degree B = 0.9990-0.9998, clock residual
<= 7.9e-9. The achieved sector ratio at the cap-limited eps = 0.05 solutions
is eps_at_solution = 0.0364 (the skirt dilutes the (2+4) share; same
phenomenon as axi3's R2). The ring-saturated G-minimisation over the same
locked basis is interior and its internal identity `E_core + E0_halo =
(3/2) G*` holds to 1.1e-16.

Sigma pulls vs the `<r1>` benchmark (c 2.37 +/- 0.09, kappa 0.802 +/- 0.018,
g* 1.31 +/- 0.04, lam 0.42 +/- 0.05, V 1.409 +/- 0.010, E_rot/E 0.2500):

| quantity | K4 (widest caps) | pull | K5 | pull |
|---|---|---|---|---|
| c_paper | 3.2627 | **+9.9 sigma** | 3.2853 | +10.2 sigma |
| kappa_paper | 0.7248 | **-4.3 sigma** | 0.7223 | -4.4 sigma |
| g* | 0.397 | -22.8 sigma | 0.311 | -25.0 sigma |
| lam_eff | 0.341 | -1.6 sigma | 0.307 | -2.3 sigma |
| V* | 7.112 | +570 sigma | 9.160 | +775 sigma |
| E_rot/E | 0.250000 | exact | 0.250000 | exact |

(lam_eff = sqrt(2<z^2>/<rho^2>) with (1-q0) weight, = 1 for a sphere; it is
NOT the paper's ansatz lambda — reported for shape reference only, and its
pull should be read with that caveat. g* = I(a*)/I(0) at d = 1;
V* = (d*/d0)^3 against the static hedgehog optimum d0.)

Stability placement of the locked solutions against the Section-3 channels
(the task's saddle-signature measurement):

- K4: kappa_ours = 0.20447 — BELOW the uniform-halo threshold 0.24430,
  ABOVE the ring/tilt threshold 0.19947 (1.025x); required destabilising
  concentration <sin^2 theta>_w = 0.9517.
- K5: kappa_ours = 0.20376 — same placement (1.021x ring); required
  concentration 0.9584.

Both solutions sit BETWEEN the two thresholds: stable against every
angular profile the K = 4/5 basis can build (the family's skirt has
exhausted its reachable concentration), unstable in principle to halos with
`<sin^2 theta>_w > 0.95` that the basis cannot represent. The cap-ladder
trend (kappa 0.752 -> 0.728 -> 0.725 -> 0.722, c 3.05 -> 3.29, monotone
toward 0.70711 / 3.385) is the family-size dial of that residual channel.

Saturated / reference values at eps = 0.05: ring-saturated locked
c_paper = 3.38547 at kappa_paper = 1/sqrt2 (G* = 2.658942);
uniform-channel hedgehog-core saturated reference c_paper = 2.51979 at
kappa_paper = sqrt(3)/2 (G_L* = 2.423818, d-relaxed hedgehog core; the G_L
minimisation over the full locked family is not well-posed — finite-
amplitude equatorial rings drive it to -inf — so only this dilation-relaxed
fixed point is quoted for the uniform channel).

## 6. Task C — locked eps-scan (t re-dialled per point) and t = 0 endpoint

All points at the final caps (L2), K = 4, warm-started, same discipline.
Achieved ratios 0.02001, 0.03501, 0.04998, 0.07497, 0.09991.

| eps | c_paper | kappa_paper | g* | lam_eff | V* | max cap | interior |
|---|---|---|---|---|---|---|---|
| 0.02001 | 3.17242 | 0.72353 | 0.2168 | 0.3452 | 12.883 | 82% | NO |
| 0.03501 | 3.21459 | 0.72435 | 0.3000 | 0.3368 | 9.348 | 96% | NO |
| 0.04998 | 3.26271 | 0.72483 | 0.3968 | 0.3405 | 7.112 | 94% | NO |
| 0.07497 | 3.34415 | 0.72592 | 0.5275 | 0.3379 | 5.406 | 100% | NO |
| 0.09991 | 3.42565 | 0.72727 | 0.6469 | 0.3378 | 4.454 | 100% | NO |
| **t = 0 endpoint** | **2.90721** | **0.75034** | 0.2053 | 0.5336 | 12.160 | **60%** | **YES (interior)** |

- The t = 0 locked endpoint is the ONE interior optimum of the whole run:
  max cap usage 60%, outer |dL/L| down to 6.5e-6 (polish 7.2e-6), K = 5
  rerun shifts c by 7.8e-3 relative (2.90721 -> 2.92996), coarse-grid shift
  3.0e-4. Its kappa_ours = 0.21167 sits between the ring (1.061x) and
  uniform (0.866x) thresholds, above the s = sin^4 member (0.20921);
  required destabilising concentration 0.888 — beyond the family's reach,
  consistent with the interior status.
- Every eps > 0 scan row is CAP-LIMITED (interior = NO). The deficit fit
  below therefore compares boundary-pinned scan values against an interior
  endpoint and is a family-/cap-conditioned measurement, not a clean law.

Deficit fit against the locked t = 0 endpoint:

```
1 - c/c0 = -0.4469 eps^0.4181        (locked, K = 4, caps L2)
```

Sign: NEGATIVE deficit — c INCREASES with eps, as in both Step-2 restricted
families (which gave -1.03 eps^1.00). The corpus law is c(eps) =
c0 [1 - 0.42 eps^(2/3)]: positive deficit, c decreasing. The locked
converged solve does NOT restore the corpus's sign; the fitted exponent
0.418 (cap-conditioned, see caveat above) matches neither the corpus's 2/3
nor Step 2's ~1.0.

## 7. Task D — locked over-spin onsets in paper units

All measured onsets of the locked problem, next to the corpus's control-run
quote (App. I.4: "emission onset at kappa = 1.000 +/- 0.004") and the
uniform locked value 0.866:

| channel | derived | measured t=0 | measured eps=0.05 |
|---|---|---|---|
| F-halo, s = |cos| (polar) | 1.11803 | 1.11805 | 1.12132 |
| **F-halo, s^2 = |cos| (ideal limit)** | **1.00000 EXACT** | — | — |
| F-halo, s^2 ~ |cos| (smoothed probe) | 0.98809 | 0.98810 | 0.99089 |
| **corpus onset quote** | **1.000 +/- 0.004** | | |
| **F-halo, s = 1 (uniform)** | **0.86603 = sqrt(3)/2** | 0.86606 | 0.86844 |
| F-halo, s = sin | 0.79057 | 0.79059 | 0.79280 |
| F-halo, s = sin^2 | 0.76376 = sqrt(7/12) EXACT | 0.76378 | 0.76595 |
| F-halo, s = sin^4 | 0.74162 | 0.74164 | 0.74381 |
| F-halo, ring limit | 0.70711 | — | — |
| tilt halo (EXCLUDED by the lock; axi3 values for reference) | 0.70711 | 0.76659 | 0.76882 |
| compacton-edge breathing/overlap shell | no threshold (dI linear: first-order relaxation) | ratio 2.059 | ratio 2.077 |

Measured/derived facts stated without interpretation:

- The locked onset family contains a member whose ideal threshold is
  kappa_paper = 1.00000 EXACTLY (`s^2 = |cos theta|`, `<sin^2 theta>_w = 1/2`
  exactly); the corpus quotes its onset at 1.000 +/- 0.004.
- The angle-uniform locked onset is sqrt(3)/2 = 0.86603 — the task sheet's
  0.866.
- The s = sin^2 locked onset is sqrt(7/12) = 0.76376 exactly (closed-form
  identity, <sin^2 th>_w = 6/7) — the same closed form as the corpus's
  deep-BPS kappa_0 = sqrt(7/12).
- The benchmark kappa = 0.802 +/- 0.018 lies 0.6 sigma from the s = sin
  locked onset 0.79057.
- The uniform-channel saturated reference on the hedgehog core at t = 0
  lands at c_paper = 2.37096, kappa_paper = 0.86603 (Section 5 definition);
  the `<r1>` benchmark tuple is c = 2.37 +/- 0.09, kappa = 0.802 +/- 0.018.
  At eps = 0.05 the same reference gives c_paper = 2.51979.
- There is no locked channel at the compacton edge: the mid-field shells
  are first-order relaxation directions (dI linear), not onsets.

## 8. Convergence, boundary and determinism evidence

- Clock residual `|2L^2/I - R|/R` at every reported solution: <= 7.9e-9
  (K4 eps=0.05: 3.6e-11; t=0: 6.3e-9; scan rows <= 6.2e-9) — the 1e-8
  residual target, met at the exactly re-solved (d, L) fixed point.
- Outer |dL/L| (p-consistency): t = 0 interior optimum 6.5e-6 (outer) /
  7.2e-6 (polish); cap-limited runs stall at 1.4e-4..1.1e-3 while c is
  still rising — the same still-climbing-at-the-wall signature axi3
  recorded for its a-mode runs (this IS the boundary-limited result, not a
  solver failure; the interior t = 0 run shows the iteration converges when
  an interior optimum exists).
- Grid ladder (N = 720 vs 360 re-evaluation of final answers): |dc|/c =
  2.1e-3 (K4 eps=0.05), 2.5e-3 (K5), 3.0e-4 (t=0 K4), 3.5e-4 (t=0 K5),
  1.8e-3..2.4e-3 (scan rows), 6.7e-4 / 1.4e-4 (ring-saturated eps / t0).
- Basis ladder (12 -> 15 params at fixed caps): c shifts by 6.9e-3 (eps =
  0.05) and 7.8e-3 (t = 0) relative; kappa by 3.5e-3 / 5.2e-3; the shape
  diagnostics g*, V* move by 20-30% (the skirt, not the clock outputs,
  absorbs the extra freedom) — closure outputs (c, kappa) are basis-stable
  at the 1e-2 level, shape diagnostics are not.
- Ring-saturated G* over the locked basis at t = 0: 2.515991, vs axi3's
  a-only (6-param, tighter caps) 2.51583 — consistent at 6e-5 relative;
  saturated internal identity residual <= 2.2e-16.
- Halo probes: dI quadratic-scaling ratio 3.999-4.000 on all 12 far-shell
  rows; t = 0 position scan flat to < 5e-5.
- Degree B = 1 measured at every accepted optimum (worst 0.99905).
- Deterministic throughout (fixed seeds, no RNG); nfev per closure 7.7k
  (full) / 5.8k (warm scan rows).

## 9. Figures

- `axi4_locked_landscape.png` — (a) F contours of the locked eps = 0.05
  optimum vs the hedgehog (d = 1 frame); (b) the contour-warp field A(r,mu)
  (equatorial stretch / polar shrink skirt); (c) the locked onset family:
  derived curve `1/sqrt(2 <sin^2 th>_w)` with all engine probe points
  (t = 0 and eps = 0.05) and the corpus lines 1.000 / 0.9354 / 0.866 /
  0.802 / sqrt(7/12) / 1/sqrt2; (d, e) c and kappa across the cap/basis
  ladder against the benchmark bands and thresholds; (f) mode amplitudes.
- `axi4_locked_epsscan.png` — (a) c(eps) vs the locked t = 0 endpoint, the
  corpus law c0[1 - 0.42 eps^(2/3)], and the benchmark point; (b) log-log
  deficit with the fitted exponent 0.418 and the corpus 0.42 eps^(2/3)
  slope; (c) kappa(eps) vs the two locked thresholds; (d) g*(eps).

## 10. Epistemic notice

Everything above is a *within-model* computation: the defined closure
(fixed-L Routhian minimisation + clock), the frozen step-1/2 conventions and
compacton unit map (zero calibrated constants), one 12/15-parameter
direction-locked deformation family with an explicit cap ladder, and exact
small-amplitude algebra for the halo channels verified by engine probes.
The cap-limited eps > 0 numbers are family-artifacts by construction (that
is what the interior = NO column records); the t = 0 endpoint is the one
interior optimum; the onset identities (1.00000, sqrt(3)/2, sqrt(7/12)) are
exact algebra of the defined functional. Nothing here is a statement about
nature, about the corpus's numerical methods, or about readings of the
corpus's text beyond the compositional definitions in force since Step 1;
which reading of the `<r1>` code the measurements support is left to the
adjudication layer.
