# Tier-2b STEP 3 — enlarged-basis (contour + tilt) isorotating closure

Reproduce with `python3 axi3_solve.py` (deterministic, no RNG; 3240.9 s wall
at full resolution; writes `axi3_results.json`, `axi3_modes.png`;
`axi3_run.log` holds every printed table). Companion files: `axi3_solve.py`,
`axi3_results.json`, `axi3_modes.png`, `axi3_run.log`. Construction and
measured numbers only; adjudication is the coordinator's layer.

## 1. Ansatz, mode basis, quadrature

General axisymmetric unit-winding field, on the phi = 0 half-plane in
spherical coordinates (r, theta), mu = cos(theta):

```
q = ( cos F,  sin F sin Theta,  0,  sin F cos Theta )
F(r,theta)     = f_base( r * exp( SUM_kl a_kl B_k(r) C_l(mu) ) )
Theta(r,theta) = theta + SUM_kl b_kl B_k(r) S_l(theta)
```

- `B_k`: K = 3 normalised cubic-B-spline bumps spanning the base support
  (centres `r_sup (k+1/2)/3`, width `r_sup/3`, peak 1).
- `C_l`: even Legendre polynomials {1, P2} (default), {1, P2, P4} (extended).
- `S_l`: odd-in-mu tilt functions {sin t cos t, sin t cos^3 t}
  (extended adds sin t cos^5 t), all vanishing at theta = 0, pi/2, pi.
- Default basis: **6 contour (a) + 6 tilt (b) modes**; extended: 9 + 9.
- Axis and equator boundary conditions (`Theta(r,0)=0, Theta(r,pi)=pi,
  Theta(r,pi/2)=pi/2`) and the hedgehog-class z -> -z symmetry are exact by
  construction; every deformation is continuous from the hedgehog, so degree
  1 is expected and is *measured* at every accepted optimum (worst reported
  |B-1| = 1.5e-4, table entries below).
- Overall dilation `d` is carried analytically (`E2 ~ d, E4 ~ 1/d,
  E6 ~ 1/d^3, E0 ~ d^3, I ~ d^3`) and minimised in closed form (golden in
  ln d) inside every objective evaluation.
- Quadrature: the Step-2 engine `axi_solve.sector_pieces` REUSED unchanged
  (cell-centred quarter-plane grid, ghost ring with exact reflection
  parities, analytic L_phi, `dV = 2 pi rho drho dz` doubled for z < 0). The
  builder uses signed `sin theta = rho/r` and `theta = arctan2(rho, z)`,
  which extends it to the ghost rings with the exact parities.
- Grids: optimisation N = 240, polish N = 480, final evaluation and gates
  N = 720, coarse rerun of every final answer N = 360 (Section 8). Fixed
  quadrature domain `1.5 * r_sup` per run so the objective is a smooth
  deterministic function of the parameters.
- Parameter caps (soft quartic walls, proximity reported per optimum):
  |a_kl| <= 0.60, |b_kl| <= 1.50, shrink cap max_A <= 0.45 on
  `A = ln(r_eff/r)`, stretch cap `-min_A <= ln(1.5) - 0.03 = 0.3755`
  (support must stay inside the domain).
- Optimiser: hand-rolled deterministic Nelder-Mead (triage over fixed seeds:
  hedgehog, spheroidal-shell P2 projections at lam = 0.7/0.45, tilt seeds,
  Step-2-family-B-like combos; then a long run + a 0.25-step restart).
  Total objective evaluations, R1 direct: none 0 / a-only 11 782 / b-only
  6 242 / joint 16 169 / extended 14 668.

Closure (direct, task-literal): from `E_tot = E_static + L^2/(2I)` and the
clock `c kappa = E_tot` (c = 2L, kappa = L/I), the solution point obeys
`L^2 = (2/3) I E_static`. This is **not** substituted into R; the iteration
is: minimise `R(p; L) = E_static(p) + L^2/(2 I(p))` over p at fixed L
(Nelder-Mead), then update L by the closed-form clock fixed point, repeat to
|dL/L| < 1e-7; after convergence the (d, L) fixed point is re-solved exactly
at frozen p*, so reported clock residuals are <= 6.7e-9 (most exactly 0.0).
`E_rot/E = 0.250000` at every reported solution (automatic identity,
verified).

## 2. Validation gates (all PASS, measured at N = 720)

| gate | what it checks | measured |
|---|---|---|
| G0 | clock closed form `L = sqrt((2/3) I E_stat)` against all four Step-2 solutions | max rel 4.4e-9 |
| G0b | Bogomolny normalisation: quadrature vs `32 sqrt2/15` (known) and vs `pi/sqrt2` (Section 4) | 1.1e-16, 2.2e-16 |
| G1 | p = 0 sectors vs step-1 radial integrals (eps = 0.05 profile) | E2 −3.4e-4, E4 −1.7e-4, E6 −1.1e-4, E0 +6.3e-6; degree 0.99991 |
| G2 | p = 0 inertia vs radial formula `(16 pi/3) INT r^2 sin^2 f dr` | +2.0e-7 |
| G3 | compacton p = 0 vs exact (`E6 = E0 = 16 sqrt2/15`, I exact) | E6 +1.9e-5, E0 +8.9e-7, I −2.4e-7; degree 1.00003 |
| G4 | deformed test config (a, b nonzero): cylindrical engine vs independent spherical-coordinate quadrature of E0 and I | E0 +9.7e-7, I +2.3e-7; degree 1.00002 |
| G5 | p = 0 closure at eps = 0.05 reproduces Step-2 solA | c_paper 2.15824 vs 2.15825 (rel −7.6e-6) |

Unit map to (e0, i0) paper units: identical to Step 2 (compacton-anchored,
zero calibrated constants); the p = 0, t = 0 closure reproduces the analytic
scaling rung `c = 2.0533, kappa = 0.9354, V = sqrt2` to 4 decimal places
(R1-none row below).

## 3. The tilted-halo channel (measured property of this energy functional)

For any configuration, add a low-amplitude tail ("halo") `F = eta` far from
the core with the internal direction tilted equatorward (`Theta ~ pi/2`,
allowed off-axis). Exact small-eta energetics in these conventions
(`E0 = (1/4pi) INT (1-cos F) dV`, `I = INT 2 sin^2 F sin^2 Theta dV`;
dE2, dE4 vanish in the soft-halo limit, dE6 is O(eta^6)):

```
dE_static ~ (1/8pi) INT eta^2 dV      dI ~ 2 INT eta^2 dV
dR = dE_static - (L^2/2I^2) dI = [INT eta^2 dV] ( 1/(8pi) - kappa_ours^2 )
```

so the fixed-L Routhian decreases under halo growth whenever
`kappa_ours > 1/sqrt(8 pi) = 0.19947`, i.e. `kappa_paper > 1/sqrt2 =
0.70711` (the isorotation channel's mass threshold in paper units, m = 1).

**Engine-level probe** (explicit shell `chi(r)` at `r_c = 2.1 r_sup`, width
`0.45 r_sup`, partial parity-exact tilt `Theta = theta + sin t cos t chi`;
all sectors measured, N = 480):

| profile | kappa_crit measured (eta = 0.02 / 0.04) | ideal (full tilt) | dI ratio eta -> 2 eta | compared solution | above threshold |
|---|---|---|---|---|---|
| t = 0 compacton | 0.21625 / 0.21628 | 0.19947 | 4.000 (quadratic) | endA kappa_ours = 0.26388 | yes |
| t = 0.0082764 (eps = 0.05) | 0.21688 / 0.21690 | 0.19947 | 4.000 (quadratic) | solA kappa_ours = 0.26558 | yes |

Measured placement of all Step-2 closure solutions relative to the
threshold: solA 1.33x, solB 1.28x, endA 1.32x, endB 1.27x critical. The
paper-unit values of the corpus tuples (deep-BPS kappa = sqrt(7/12) =
0.7638, `<r1>` kappa = 0.802) also lie above 1/sqrt2 = 0.70711.

## 4. Saturated closure: formulation, what is proven vs computed

At `kappa_ours > 1/sqrt(8pi)` the halo grows until `L^2/I^2 = 1/(8 pi)`
exactly; the halo contribution to `G := E_static - I/(16 pi)` is zero at
threshold, so the closure factorises (derivation in the script header,
algebra exact):

```
kappa_ours* = 1/sqrt(8pi)   (kappa_paper* = 1/sqrt2, EXACT)
core minimises G;  L = sqrt(8pi) G*,  E_stat_tot = (3/2) G*,
I_tot = 8pi G*,  R = 2 G*,  E_rot/E = 1/4,  c_paper = (4/pi) G*.
```

Status of each statement:
- **Proven (exact algebra, given the closure definitions):** the threshold
  formula for dR; the saturation values `kappa_paper = 1/sqrt2` and
  `c_paper = (4/pi) G*`; halo G-neutrality. Verified numerically: the
  internal identity `E_stat_core + E0_halo = (3/2) G*` holds to 2.2e-16 at
  every reported saturated solution; the halo channel coefficients are
  verified by the Section-3 probe; the saturation is approached only in the
  soft/delocalised halo limit (at t > 0 the E2 cost of a finite halo is
  positive and -> 0 with its wavelength), so the saturated tuple is the
  *limit* of the direct closure, not a normalisable interior optimum.
- **Proven (Bogomolny, normalisation verified in G0b):** pointwise
  `sin^2 Theta <= 1` gives `G >= E6 + (1/8pi) INT (1-cos F)^2 dV`, and the
  sextic+potential bound with `U_eff = (1-cos f)^2/2` gives
  `G >= 2 INT_0^pi sin^2 f sqrt(U_eff) df = pi/sqrt2` (closed form), hence
  `c_paper >= 2 sqrt2 = 2.82843` for every t >= 0.
- **Computed (family-limited upper value):** `G*` minimised over the
  enlarged basis (well-posed: interior in all but near-neutral halo-like
  directions), Section 6.

## 5. R1 (t -> 0, compacton base): the decisive run

### 5a. Direct closure (task-literal), 2x2 mode-group ablation, N = 720

| config | c_paper | kappa_paper | g* | lam_eff | V | E_rot/E | degree B | boundary walls hit |
|---|---|---|---|---|---|---|---|---|
| none (hedgehog + dilation) | 2.0533 | 0.9354 | 1.0000 | 1.0000 | 1.4142 | 0.250000 | 1.00003 | — |
| a-only (contour) | 2.5906 | 0.8022 | 0.5363 | 0.7407 | 3.880 | 0.250000 | 1.00002 | amax (100% of cap) |
| b-only (tilt) | 2.1901 | 0.8946 | 1.0933 | 1.0000 | 1.4426 | 0.250000 | 1.00003 | — |
| joint | 2.7752 | 0.7718 | 0.6229 | 1.0168 | 3.719 | 0.250000 | 0.99999 | amax (100% of cap) |
| extended (18 modes) | 2.7673 | 0.7727 | 0.6501 | 1.1083 | 3.549 | 0.250000 | 0.99998 | amax (100% of cap) |

(g* = I_1(p*)/I_1(0) at d = 1; V = (d*/d0)^3; V_E0 = E0-ratio definition
also recorded in the JSON: 1.530/1.443/1.577/1.574 for a/b/joint/ext.)

Fixed-L 2x2 Routhian table at L* = 10.9271 (joint solution's spin):
`R(none) = 5.0726`, `R(a-only) = 4.7820` (dR = −0.2906), `R(b-only) =
4.9591` (dR = −0.1135), `R(joint) = 4.7582` (dR = −0.3144). Tilt modes
lower R on their own and beyond what contour modes alone reach; the joint
gain exceeds either group.

**Boundary behaviour (measured):** every configuration containing a-modes
terminates ON the shrink-cap wall (amax at 100% of its cap), with the outer
L-iteration still creeping upward at the 12-iteration cap (last |dL/L| ~
1e-4..8e-5); the b-only run converges to an interior optimum (|dL/L|
2.2e-8). Family-size scan with both A-walls dialled together, joint config:

| cap on max A | stretch cap | c_paper | kappa_paper | V | g* | wall hit |
|---|---|---|---|---|---|---|
| 0.30 | 0.2324 | 2.6050 | 0.7996 | 2.684 | 0.782 | amax 100% |
| 0.45 | 0.3755 | 2.7752 | 0.7718 | 3.719 | 0.623 | amax 100% |
| 0.65 | 0.5296 | 2.9405 | 0.7454 | 5.894 | 0.431 | amax 100% |

c rises and kappa falls monotonically with every enlargement, with V
growing without sign of settling — the direct closure has no interior
optimum in any of these families; the trend runs toward the Section-4
saturated values.

### 5b. Saturated closure (G-minimisation over the same basis), N = 720

| config | G* | c_paper = (4/pi) G* | g_tot | g_core (at scale) | V_E0 (incl. halo) | halo Delta | degree B |
|---|---|---|---|---|---|---|---|
| hedgehog core | 2.54985 | 3.2466 | 2.958 | 1.183 | 1.690 | 19.23 | 1.00003 |
| a-only core | 2.51583 | 3.2033 | 2.919 | 1.428 | 1.668 | 16.15 | 1.00002 |
| b-only core | 2.53767 | 3.2311 | 2.944 | 1.255 | 1.682 | 18.29 | 1.00003 |
| joint core | 2.51445 | 3.2015 | 2.917 | 1.447 | 1.667 | 15.92 | 1.00002 |
| extended core | 2.51442 | 3.2015 | 2.917 | 1.449 | 1.667 | 15.90 | 1.00003 |

kappa_paper = 1/sqrt2 = 0.70711 exact at saturation; rigorous lower bound
G >= pi/sqrt2 = 2.22144. Measured t -> 0 bracket:

```
c_paper(t -> 0)  in  [ 2.82843 (rigorous bound) , 3.20146 (family value) ]
```

The extended basis moves G* by 2.7e-5 (c by 3.5e-5) relative to the default
basis — converged within this family; the joint/extended core carries a mild
amax-wall contact (halo-like skirt directions are G-flat to G-slightly
positive, so the wall costs G nothing measurable: a-only, interior, differs
by 3.8e-4 in G). For reference: the paper's shape-exact endpoint values are
(c, kappa) = (2.5147, 0.7638) and the scaling rung is (2.0533, 0.9354);
both c-values lie below the measured bracket's lower edge, and neither
kappa equals 1/sqrt2.

## 6. R3: marginal-mode spectrum (d2 E_static/dp^2, dI/dp), N = 480

At the **compacton** (p = 0, d = 1; E_static = 3.017046, I = 21.6643):

| mode | dE/dp | d2E/dp2 | dI/dp | dI/d2E |
|---|---|---|---|---|
| a[0,0] | +0.0007 | 1.693 | +2.172 | +1.28 |
| a[0,2] | +0.0001 | 0.338 | −0.434 | −1.28 |
| a[1,0] | +0.0049 | 11.571 | −3.326 | −0.29 |
| a[1,2] | +0.0003 | 2.308 | +0.664 | +0.29 |
| a[2,0] | +0.0036 | 25.276 | −62.911 | −2.49 |
| a[2,2] | +0.0002 | 5.047 | +12.610 | +2.50 |
| b[0,1] | +0.0001 | 0.286 | +0.655 | +2.29 |
| b[0,2] | +0.0001 | 0.105 | +0.281 | +2.69 |
| b[1,1] | +0.0004 | 1.192 | +4.581 | +3.85 |
| b[1,2] | +0.0002 | 0.435 | +1.964 | +4.51 |
| b[2,1] | +0.0003 | 1.061 | +6.437 | +6.07 |
| b[2,2] | +0.0002 | 0.387 | +2.759 | +7.12 |

Measured facts: (i) every dE/dp is zero to the finite-difference floor
(<= 5e-3 of the mode scale) — the compacton is E-stationary in all 12
directions; (ii) every mode carries dI/dp != 0, i.e. inertia is bought at
*quadratic* energy cost, never free, and never refused; (iii) the tilt (b)
modes have the highest inertia-per-stiffness ratios (dI/d2E up to 7.1, and
up to 29.5 at the R1 direct optimum — outer-region tilt is the cheapest
inertia in the family, the finite-basis shadow of the Section-3 halo
channel); (iv) the largest raw dI belongs to a[2,0]/a[2,2] (edge contour
modes — skirt formation). Full table at the R1 direct-joint optimum
(d* = 1.5493) is in `axi3_run.log` and `axi3_results.json`.

## 7. R2 (eps = 0.05, t frozen from step 1)

| solution | c_paper | kappa_paper | g* | V | V_E0 | boundary |
|---|---|---|---|---|---|---|
| direct joint (boundary-limited) | 3.1871 | 0.7325 | 0.658 | 4.147 | 1.676 | stretch 95%, amax 100%, bbox |
| saturated (unrestricted limit) | 3.3812 | 0.70711 (exact) | g_tot 2.997 (core 1.561) | V_core 2.091 | 1.713 | amax, bbox (core) |
| saturated, hedgehog core (reference) | 3.4402 | 0.70711 | 3.050 | — | 1.735 | — |

Against the `<r1>` benchmark (c = 2.37 ± 0.09, kappa = 0.802 ± 0.018,
g* = 1.31 ± 0.04, V = 1.409 ± 0.010):

| quantity | direct | pull | saturated | pull |
|---|---|---|---|---|
| c | 3.1871 | +9.1 sigma | 3.3812 | +11.2 sigma |
| kappa | 0.7325 | −3.9 sigma | 0.70711 | −5.3 sigma |
| g* | 0.658 (d=1 def.) | — | 2.997 (total) | +42 sigma |
| V | 4.147 / V_E0 1.676 | — | V_E0 1.713 | +30 sigma (V_E0 def.) |

Measured facts: the direct R2 closure hits three walls at once and its
outer iteration is still climbing at the cap (final |dL/L| = 2.7e-4,
`outer_converged = false` recorded); the achieved sector ratio at its
boundary-pinned solution is eps = 0.0319 (the halo-like skirt dilutes the
(2+4) share). The `<r1>` tuple sits at kappa_ours = 0.802/KMAP = 0.226 >
1/sqrt(8pi): by the Section-3 measured threshold it is a stationary point
with a strictly downhill fixed-L direction in this family (a saddle of the
Routhian), not its minimum — that placement is the measured fact; what it
means for the corpus is left to adjudication. The rigorous floor
`c_paper >= 2.828` applies at eps = 0.05 as well (the bound drops
t E2 + t E4 >= 0).

## 8. Convergence and determinism evidence

- Final answers re-evaluated on the coarse grid (N = 360 vs 720):
  |delta c|/c = 4.1e-5 (R1 direct joint), 4.8e-5 (R1 ext), 1.0e-6 (R1
  saturated joint, G to 2.5e-6), 4.4e-4 (R2 direct), 2.4e-4 (R2 saturated).
- Gate quantities at N = 720 sit at the 1e-4..1e-7 level against the exact
  and step-1 references (Section 2 table).
- Saturated-closure internal identity `E_core + E0_halo = (3/2) G*`:
  residual 2.2e-16 (R1), 1.1e-16 (R2).
- Halo probe quadratic scaling: dI(2 eta)/dI(eta) = 4.000 at both t values.
- Everything is deterministic (fixed seeds, no RNG); rerunning reproduces
  the log byte-for-byte apart from timings.

## 9. Figure

`axi3_modes.png`: (a) F-contours of the R1 saturated-core optimum vs the
compacton (oblate-waisted core); (b) the internal-direction tilt field
Delta Theta (equatorward everywhere off-axis, peak ~0.39 rad); (c) optimal
mode amplitudes R1 vs R2 (R2 tilt coefficients b[2,1], b[2,2] at the 1.5
cap); (d) the compacton marginal spectrum (Section 6) — inertia gain vs
quadratic cost per mode; (e, f) the family-size scan of the direct closure
with the saturated asymptotes `(4/pi) G*_fam` and `1/sqrt2` and the corpus
reference lines.

## 10. Epistemic notice

Everything above is a *within-model* computation: the defined closure
(fixed-L Routhian minimisation + clock), the defined sector conventions
(a6 = a0 = m = 1, step-1/step-2 unit map, zero calibrated constants), and a
specific 12/18-parameter deformation family plus one exact bound over the
unrestricted field space. The boundary-limited direct numbers are
family-artifacts by construction (that is what the walls-hit columns
record); the saturated numbers are limits of the defined closure, exact in
kappa and bracketed in c. Nothing here is a statement about nature, about
the paper's numerical methods, or about readings of the paper's text other
than the compositional definitions used since Step 1.
