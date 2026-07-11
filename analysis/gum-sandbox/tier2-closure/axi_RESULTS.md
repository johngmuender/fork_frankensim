# Tier-2b STEP 2 — spheroidal-ansatz isorotating closure at eps = 0.05 (field level, zero calibrated constants)

Reproduce with `python3 axi_solve.py` (deterministic, no RNG; ~170 s; writes
`axi_results.json`, `axi_landscape.png`, `axi_epsscan.png`; printed table in
`axi_run.log`). Grid: quarter-plane cell-centred 512x512 (= half-plane
512x1024; z-parity `P(rho,-z) = sigma1 P(rho,z) sigma1` derived and used),
L_rho/L_z by central differences of P, L_phi analytic from the winding;
coarse rerun at 336x336 for convergence. Adjudication of what these results
mean for the corpus is filed separately (`axi_ADJUDICATION.md`, F-R4); this
report records the construction and the measured numbers.

## Constants (derived, then verified — nothing fitted)

Evaluating the specified densities on the hedgehog and matching step 1's
radial-integral convention (E_rad = E_3D/4pi) fixes every sector constant
analytically **before** any closure:

| sector | density (orthonormal frame D_i) | constant |
|---|---|---|
| E2 | `-(1/2) Tr(D_i D_i) = sum |a_i|^2` | 1/4pi |
| E4 | `-(1/8) Tr([D_i,D_j]^2) = |a_i x a_j|^2` | 1/4pi |
| E6 | `b^2`, `b = -(1/2pi^2) a_rho.(a_z x a_phi/rho)` | pi^3 |
| E0 | `1 - q0` | 1/4pi |
| I  | `2 (q1^2 + q2^2)` | 1 (full 3-D integral, task convention) |

Gates V1/V2 confirm all five at the quadrature-error level (<= 3.0e-4).

**Unit map to the paper's (e0, i0) convention** — fixed with zero calibration
by the exact eps = 0 compacton, which both conventions possess in closed form:
`u_E = e0/(2 E6_compacton) = sqrt2/pi`, `u_I = i0/I_compacton = 2^(-3/2)/pi^2`,
hence `c_paper = c_ours/sqrt(2 pi^3)`, `kappa_paper = 2 sqrt(pi) kappa_ours`.
Identity check: the family-A t = 0 closure must land on the reduced model's
g = 1 endpoint `c = 2 sqrt(e0 i0) = 2.0533`, `kappa = 0.9354`, `V = sqrt2` —
it does (below, to 4-5 digits). Note the paper's advertised eps->0 endpoint
(2.5147, 0.7638) is its g = 3/2 pancake endpoint, which presumes the inertia
enhancement; see FINDING.

## Deformation families

- **A — the task's literal family (SDiff pullback + dilation):**
  `q_(lam,d)(x) = q_base(x/(ed), y/(ed), z/(pd))`, `e = lam^(-1/3)`,
  `p = lam^(2/3)`. Pullback to base coordinates gives exact closed-form
  (lam, d) dependence of all sectors from one base evaluation (verified
  against an independent, incommensurate lab-grid run in V3).
- **B — the paper's spheroidal-shell ansatz (NOT an SDiff image):**
  `q = (cos f(m), sin f(m) xhat_lab)`, `m^2 = (rho/e)^2 + (z/p)^2`, then
  `x -> x/d`. This is the minimal field realisation of the paper's
  `I = 2J INT sin^2 f sin^2 theta_lab`; it reproduces `g_oblate(lam)` exactly.

**FINDING (measured; exact by change of variables):** any functional with no
spatial gradients — E0 *and* the inertia density `2(q1^2+q2^2)` — is invariant
under every volume-preserving diffeomorphism. In family A, `g(lam) == 1`
identically (measured 1 - g < 2e-15 at lam = 0.6 and 0.42 on independent
grids). The paper's g(lam) enhancement exists only in family B, whose E6 is
then no longer flat (table below). At t = 0 the family-A lambda direction is
exactly degenerate — no pancake runaway, no gradient in lambda at all
(the reported lam* there is an arbitrary point of the flat valley).

## Validation gates (N = 512)

| gate | quantity | measured | tolerance | verdict |
|---|---|---|---|---|
| V1 | E2 vs 12.117459 | -3.00e-4 | 5e-3 | **PASS** |
| V1 | E4 vs 6.191126 | -1.54e-4 | 5e-3 | **PASS** |
| V1 | E6 vs 1.523936 | -9.89e-5 | 5e-3 | **PASS** |
| V1 | E0 vs 1.507588 | +5.55e-6 | 5e-3 | **PASS** |
| V1 | degree B | +0.999924 | — | PASS (sign fixes b-orientation) |
| V2 | I_2D vs (16pi/3) INT r^2 sin^2 f dr = 22.26826 | +2.0e-7 | 5e-3 | **PASS** |
| V3 | (E6+E0) shift, lam = 0.6, independent iso-grid | -5.1e-5 | 2e-3 | **PASS** (SDiff-flat BPS sectors) |
| V3 | lab vs closed-form A: E2, E4, E6, E0, I | -2.7e-4, -1.4e-4, -1.0e-4, 2e-15, 1e-15 | — | cross-check PASS |
| V4 | I(0.42)/I(1), family A | 1.000000 | expected 1.2855 | **= 1 exactly** (FINDING; 28.6% below the paper curve) |
| V4 | I(0.42)/I(1), family B | 1.285547 | g_oblate = 1.285547 | matches to 3e-10 (paper's curve is exactly the shell-ansatz sin^2-average) |

(The V3/V4 family-A lab runs use square cells incommensurate with the
pullback of the base grid; a deformation-commensurate grid would reproduce
the closed forms to machine precision by construction and test nothing.)

## Closure at eps = 0.05 (t = 0.0082764349 frozen) and t = 0 endpoints

`E_static = t E2 + t E4 + E6 + E0`; `R = E_static + L^2/(2I)`; minimise over
(lam, d); clock `2L^2/I* = R*` by bisection (clock residual <= 1e-8).

| quantity | A, t=0 | **A, eps=0.05** | B, t=0 | **B, eps=0.05** | target |
|---|---|---|---|---|---|
| L | 8.0847 | 8.4979 | 8.5392 | 8.9728 | — |
| c_ours = 2L | 16.1694 | 16.9958 | 17.0785 | 17.9456 | — |
| c_paper | 2.0533 | 2.1583 | 2.1688 | **2.2789** | 2.37 +/- 0.09 |
| kappa_ours = L/I* | 0.26388 | 0.26558 | 0.25414 | 0.25583 | — |
| kappa_paper | 0.9354 | 0.9414 | 0.9009 | **0.9069** | 0.802 +/- 0.018 |
| lambda* | (flat) | 1.0000 | 0.8173 | **0.8180** | 0.42 +/- 0.05 |
| g* = I(lam*,d*)/I(1,d*) | 1 (exact) | 1.0000 | 1.0781 | **1.0777** | 1.31 +/- 0.04 |
| V* = (d*/d0)^3 | 1.41421 | 1.4370 | 1.4386 | **1.4615** | 1.409 +/- 0.010 |
| E_rot/E_tot | 0.250000 | 0.250000 | 0.250000 | **0.250000** | 0.2500 |
| eps at solution | 0 | 0.04913 | 0 | 0.04920 | (0.05 at sphere, d=1) |

Family-A t = 0 closed form (`s = L^2/2I = E6`, `V = sqrt2`,
`c_paper = 2 sqrt(e0 i0)`): matched to 6e-5 relative — internal minimiser +
quadrature check. Compacton-engine accuracy: E6 +1.7e-5, E0 +8.6e-7, I +1.1e-7
vs exact values.

## Why lambda* stops at 0.82: the field-level cost/gain (family B, eps = 0.05, d = 1)

| lam | E6(lam)/E6(1) | I(lam)/I(1) [= g_oblate] | E0(lam)/E0(1) |
|---|---|---|---|
| 1.000 | 1.0000 | 1.0000 | 1 (exact at all lam) |
| 0.825 | 1.0315 | 1.0746 | 1 |
| 0.650 | 1.1744 | 1.1590 | 1 |
| 0.425 | 1.8716 | 1.2826 | 1 |
| 0.250 | 4.4263 | 1.3869 | 1 |
| 0.150 | 11.516 | 1.4443 | 1 |

The sextic cost is O(1) and steep (at lam = 0.42 it is DE6 = +1.33 in our
units against a rotational gain of ~0.25 at the solved spin), while the
inertia gain saturates at 3/2. The balance lands at lam* = 0.818,
g* = 1.078 — far short of the paper's (0.42, 1.31), and lam*(eps) is nearly
eps-independent (0.8177 -> 0.8183 over the whole scan): the shape optimum is
set by the BPS-sector curvature, not by the eps-suppressed sectors.

## eps-scan (t re-dialled per point; achieved ratios 0.02001, 0.03501, 0.04998, 0.07497, 0.09991)

| eps | c_paper A | c_paper B | lam*_B | g*_B | V*_A | V*_B | kappa_paper B |
|---|---|---|---|---|---|---|---|
| 0.020 | 2.0952 | 2.2126 | 0.8177 | 1.0779 | 1.4233 | 1.4477 | 0.9032 |
| 0.035 | 2.1267 | 2.2457 | 0.8179 | 1.0778 | 1.4302 | 1.4546 | 0.9050 |
| 0.050 | 2.1583 | 2.2789 | 0.8180 | 1.0777 | 1.4370 | 1.4615 | 0.9069 |
| 0.075 | 2.2108 | 2.3342 | 0.8182 | 1.0777 | 1.4484 | 1.4730 | 0.9101 |
| 0.100 | 2.2634 | 2.3896 | 0.8183 | 1.0776 | 1.4599 | 1.4846 | 0.9134 |

**Deficit exponent fits** (vs the same-family t = 0 endpoint):

- family A: `1 - c/c0 = -1.0305 eps^1.0026`
- family B: `1 - c/c0 = -1.0316 eps^1.0052`

Both are **exponent ~ 1.00 with NEGATIVE deficit** (c grows with eps): the
eps-sectors are a positive penalty that raises R* and therefore the clock
charge, while the shape optimum barely moves. The corpus law
(+0.42 eps^(2/3), c *decreasing*) requires the shape optimum to migrate with
eps along a statically flat direction (as in the reduced model, where
g(lam*(eps)) falls from 3/2); at field level that flat direction with
g-variation does not exist in either family, so both the exponent (1.00 vs
2/3) and the sign disagree. Tier-2a's reduced model gave 0.75 with the
correct sign — the field level does not reproduce it.

## Grid convergence (N = 336 vs 512, eps = 0.05)

All closure outputs move by <= 8.1e-5 relative (c: 7.2e-5 A / 8.1e-5 B;
lam*: 1.6e-5 / 5.1e-5; g*: 0 / 1.8e-5; V*: 4.6e-6 / 3.8e-6; R*: 7.5e-5);
gate quantities improve as O(h^2) (V1 max err 7.0e-4 -> 3.0e-4; V3
-1.2e-4 -> -5.1e-5; g_B(0.42) unchanged to 5 dp). Quadrature error is far
below every stated tolerance and far below all physics discrepancies.

## Figures

- `axi_landscape.png` — Routhian over (lam, d) at the solved L for both
  families; family A's valley is flat in lambda with the minimum at
  lam = 1.000 (marginally selected by the eps-sectors), family B's minimum
  sits at (0.818, 1.135); static d0 marked.
- `axi_epsscan.png` — c(eps) in paper units vs both endpoints and the
  benchmark; log-log deficit with fitted slopes (A and B nearly coincide);
  lam*(eps); g*(eps).

## Verdicts vs the <r1> benchmark tuple (family B = the paper's own ansatz)

| target | family A (literal SDiff family) | family B (shell ansatz) | verdict |
|---|---|---|---|
| lambda* = 0.42 +/- 0.05 | 1.000 (sphere) | 0.818 (+8.0 sigma) | **FAIL both** |
| g* = 1.31 +/- 0.04 | 1.000 (exactly; theorem) | 1.078 (-5.8 sigma) | **FAIL both** |
| V* = 1.409 +/- 0.010 | 1.437 (+2.8 sigma, above sqrt2) | 1.462 (+5.2 sigma) | **FAIL both** (benchmark sits *below* sqrt2; field level lands *above*) |
| E_rot/E = 0.2500 | 0.250000 | 0.250000 | **PASS** (automatic spin+clock kinematics, as claimed) |
| c = 2.37 +/- 0.09 (paper units) | 2.158 (-2.4 sigma) | 2.279 (-1.0 sigma) | A FAIL / **B marginal PASS** |
| kappa = 0.802 +/- 0.018 (paper units) | 0.941 (+7.7 sigma) | 0.907 (+5.8 sigma) | **FAIL both** (benchmark kappa is below even the g = 1 eps->0 endpoint 0.9354) |
| deficit exponent 2/3-class | 1.003, wrong sign | 1.005, wrong sign | **FAIL both** |

c and kappa carry the compacton unit map (u_E, u_I above) — the map itself is
exact, but any unstated convention offset in the corpus would shift both;
lambda*, g*, V*, E_rot/E and the exponent are convention-free as reported.

## Epistemic notice

Within-model replication of a speculative corpus's internal consistency;
nothing here bears on nature. Interpretation of the failures (which reading
of the corpus's G.5 mechanism they refute) is in `axi_ADJUDICATION.md`.
