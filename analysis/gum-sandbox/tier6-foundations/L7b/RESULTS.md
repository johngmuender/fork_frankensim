# L7b — spin-family affinity test: POVM exclusion on the A3 observable class (F-T6-L7b)

**Goal:** Instantiate the missing T4-W5 exclusion step of the corpus's
Sec. VIII.F gate logic, at testbed grade, on the A3-relevant observable:
test whether the model's Bohmian first-crossing statistics, taken over a
FAMILY of spin preparations chi(n-hat) (experimentally: a spin rotator
before the waveguide), are realizable by any single fixed POVM. Field
configuration and engine are inherited unchanged from the validated
L5prime bench (k0 = 2, C2 smoothstep taper, hard-wall waveguide,
d_near = 1, N = 2000, seed 20260717, T = 16).

Within-model; nothing here bears on nature.

## The operational schema

A fixed physical detector induces, via the Naimark construction, ONE
preparation-independent POVM {E(dt)} on the system Hilbert space. So over
a FAMILY of preparations psi, POVM-realizable arrival statistics must be
a quadratic (sesquilinear) functional Pi_psi(B) = <psi|E(B)|psi>. A
single preparation excludes nothing (any one distribution is trivially a
POVM); all content lives in the family. A measured violation over the
family IS the exclusion step the corpus's Sec. VIII.F gate logic is
missing (T4-W5) — conditional, exactly as VIII.F requires, on these
statistics being what the experiment observes.

## Key reduction (why the test is one-dimensional)

With no magnetic field, psi = phi(x,z,t) chi(n-hat) with phi
spin-independent: ONE field evolution total. The Pauli guidance current is
j = Im(phi* grad phi) + (n_y/2)(-d rho/dz, +d rho/dx) — the spin enters
ONLY through the out-of-plane Bloch component n_y, so the velocity field
is affine at the field level, v(x,z,t; n_y) = a(x,z,t) + n_y b(x,z,t)
(the trajectory map and its first-crossing statistics need not be, and
are not). A POVM gives per bin Pi_B(n-hat) = A + B . n-hat (affine on the
Bloch sphere, from the 2x2 effective operator <phi|E_B|phi>); pure states
with the same n_y but different (n_x, n_z) have IDENTICAL dynamics,
forcing B_x = B_z = 0; hence Pi_B must be AFFINE IN n_y:
Pi_B(n_y) = A + B n_y. Any measured nonlinearity in n_y, or any
interval-vanishing-with-positivity-elsewhere pattern, excludes every POVM
for this family.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | norm 1e-6; wall amplitude < 1e-10; energy drift 1e-6; n_y = +1 and n_y = 0 reproduce L5prime transverse/axial (tau_max within 1e-3 rel; KS p > 0.5) | ensemble bit-identical to L5prime (max abs diff 0.0); 1-D engine: norm dev 9.0e-13, energy drift 2.9e-13 (rel); 2-D diagnostic: norm dev 2.2e-16, wall max 3.4e-16, energy drift 1.2e-16, z-mode leakage 0.0; tau_max(n_y=+1) = 5.1309462 vs L5prime 5.1309505 → rel 8.3e-7; KS vs L5prime: transverse stat 0.001, p = 1.000; axial stat 0.0035, p = 1.000 | **PASS** |
| G2 | affine fit Pi_B(n_y) = A + B n_y per bin (40 bins on [0,16] + non-arrival bin; >= 20 expected counts in populated mid-range); pre-registered exclusion: >= 3 independent bins with max residual > 5 sigma_boot (2000 resamples) | **8 bins > 5 sigma** (bins 12–19, t in [4.8, 8.0]); max violation **40.2 sigma** at t in [5.2, 5.6]; mid-range bins 3–12 have >= 20 counts at every n_y; quadratic diagnostic: C-term significant (|C|/sigma_C > 5) in 9 bins, max 12.1 sigma, C < 0 throughout the late-time violation band | **PASS** |
| G3 | >= 1 time bin with ZERO arrivals (0/2000, Pi < 1.9e-3 at 95% CL) for every n_y in an interval while Pi > 1e-2 measured outside it; zeros verified at one refined resolution | **2 kill bins**: t in [6.4, 6.8] and [6.8, 7.2] have 0/2000 arrivals for EVERY n_y in {-1, -0.75, +0.75, +1} (each Pi < 1.50e-3 at 95% CL) while Pi = 0.0120 resp. 0.0115 at n_y = +0.25 (counts 22/15/11/24/17 resp. 18/9/14/23/16 across the five inner n_y); zeros confirmed at double Nx (4096) AND on the extended box [-15, 105] (wrap-around control) | **PASS** |
| G4 | tau_max(n_y) varies by > 20% across the family and is monotone in \|n_y\|, or reported otherwise with diagnosis | tau_max: 5.131 (\|n_y\| = 1) → 6.01 (0.75) → 7.94 (0.5) → 15.69 (0.25) → 15.84 (0): variation **209%**, strictly monotone decreasing in \|n_y\| (both branches separately and averaged); empty-gap check: gap >= 2 IQR holds for all \|n_y\| >= 0.5 (gap 8.02–10.87 vs 2 IQR 2.35–3.07); for \|n_y\| <= 0.25 the support runs to ~15.8 — no cutoff inside the window (n_y = 0: 20 non-arrivals), as in L5prime's axial case | **PASS** |
| G5 | (report-only) effect sizes | max linearity violation 40.2 sigma (bin [5.2, 5.6]); mean arrival time 2.938 (n_y = -1) → 3.546 (n_y = +0.25), a 20.7% swing, even-dominated (maximum near small \|n_y\|, minima at \|n_y\| = 1); max quadratic-term significance 12.1 sigma | REPORT-ONLY |

## The affine-vanishing lemma (G3 kill)

An affine function Pi_B(n_y) = A + B n_y that vanishes at >= 2 distinct
n_y values (a fortiori on an interval of n_y values) is identically zero
— contradicting the measured positivity at interior n_y. Quantitatively:
any affine Pi_B with Pi_B(-1) and Pi_B(+1) both < 1.50e-3 (the 95% CL
bound from 0/2000) forces Pi_B(0) = (Pi_B(-1) + Pi_B(+1))/2 < 1.50e-3
and by convexity Pi_B < 1.5e-3 on all of [-1, 1] — a factor 8 below the
measured Pi = 0.0120 at n_y = +0.25. No single POVM reproduces the
family.

## Key numbers

- tau_max(n_y), reference run (Nx = 2048, dt = 2.5e-4, box [-15, 45]):

  | n_y | -1 | -0.75 | -0.5 | -0.25 | 0 | +0.25 | +0.5 | +0.75 | +1 |
  |---|---|---|---|---|---|---|---|---|---|
  | tau_max | 5.1289 | 6.0263 | 7.9753 | 15.758 | 15.836 | 15.627 | 7.9134 | 5.9993 | 5.1309 |
  | mean tau | 2.9376 | 3.1019 | 3.2927 | 3.4859 | 3.4759 | 3.5461 | 3.3105 | 3.0755 | 2.9589 |
  | non-arrivals | 0 | 0 | 0 | 0 | 20 | 0 | 0 | 0 | 0 |

- Affine-fit violation profile (max over the 9 n_y of |residual|/sigma_boot),
  bins over 5 sigma: bin 12 [4.8, 5.2] 6.2; bin 13 [5.2, 5.6] **40.2**;
  bin 14 [5.6, 6.0] 24.4; bin 15 [6.0, 6.4] 20.5; bin 16 [6.4, 6.8] 10.0;
  bin 17 [6.8, 7.2] 9.6; bin 18 [7.2, 7.6] 7.2; bin 19 [7.6, 8.0] 6.6.
  Next largest: 4.1 (non-arrival bin), 4.0 (bins [1.6, 2.0], [2.8, 3.2]).
  Full 41-bin profile in l7b_results.json.
- Statistics are even in n_y to sampling accuracy (mirror symmetry
  z -> 1 - z, n_y -> -n_y of guidance law + z-symmetric preparation):
  +/- pair differences are O(0.01-0.03) in mean tau, consistent with the
  finite-sample z-mirror; consequently the best affine fit has B ~ 0 and
  the violation is carried by the even part — the quadratic C-term is
  negative and 5-12 sigma significant across the whole violation band
  t in [4.8, 8.4] (support shrinks with n_y^2), 9 bins total > 5 sigma.
- Kill-bin counts (bin [6.4, 6.8]), n_y = -1 ... +1:
  0, 0, 22, 15, 11, 24, 17, 0, 0 (of 2000 each). 95% CL upper bound for
  0/2000: p < 1.497e-3 (< 1.9e-3 spec). Confirmed 0 at Nx = 4096 and on
  the extended box.
- Zero-of-2000 stability: tau_max per n_y changes by < 3e-6 (|n_y| >= 0.5)
  under Nx doubling; at \|n_y\| <= 0.25 the endpoint (a slow-tail order
  statistic near the window edge) moves by up to 0.19 between the
  reference and the extended box — it does not enter any gate margin.
- Engine: one field evolution per resolution drives all 9 spin ensembles
  (18,000 trajectories); reference 147 s, refine (Nx x 2) 157 s, extended
  box 159 s. Frozen box escapees (trajectories reaching within 0.5 of a
  box edge, frozen there): 110-155 per n_y, all long after crossing d_near.
- G1 cross-checks: tau_max(n_y = +1) = 5.1309462 = L5prime's
  pre-registered-box control value to all printed digits; L5prime
  extended-box reference 5.1309505 (rel diff 8.3e-7); KS(our n_y = 0 vs
  L5prime axial): stat 0.0035, p ~ 1.0 (1980 vs 1976 window-16 arrivals —
  the 4-trajectory difference is the documented wrap effect).

## Method

Factorized exact-in-z engine reused verbatim from L5prime:
phi = f(x, t) sin(pi z) e^{-i pi^2 t / 2}, f propagated exactly in 1-D
Fourier space (Nx = 2048 on [-15, 45], dt = 2.5e-4, T = 16). Initial
state: Gaussian exp(-(x-x0)^2/(2 sigma_x^2)) x S2-smoothstep taper
(W = 1 for |x-x0| <= 2.5 sigma, 0 at 3.5 sigma) x e^{i k0 x}, sigma_x = 1,
x0 = -5, k0 = 2, ground mode sin(pi z). N = 2000 positions sampled from
|phi(0)|^2 (analytic inverse-CDF, seed 20260717) — verified bit-identical
to L5prime's ensemble — and reused for every n_y. Guidance:
v_x = Im(f'/f) - n_y pi cot(pi z), v_z = n_y Re(f'/f) (the exact
factorized form of j/rho with the Pauli spin term; n_y = +1 and n_y = 0
reproduce L5prime's transverse and axial branches bit-for-bit). RK4 on
stored fields at t, t + dt/2, t + dt; arrival = first crossing of x = 1;
n_y in {-1, -0.75, ..., +1}, all 9 ensembles advanced inside one field
evolution. Statistics: 40 uniform bins on [0, 16] + non-arrival bin;
per-bin LSQ fits (affine and quadratic in n_y) over the 9 family points;
per-(bin, n_y) bootstrap sigma from 2000 resamples with shared resample
indices (paired ensembles), sigma floored at the 1-count level 3.5e-4 to
keep zero-count bins finite. Refinements: double Nx (4096, same box) and
extended box [-15, 105] at identical dx (wrap-around control). 2-D FFT
wall diagnostic (Nx = 2048, Nz = 128, 17 spectral checkpoints) for the
wall-amplitude gate.

## Realization

The family's first-crossing statistics are maximally far from the
POVM-required affine form: they are (to sampling accuracy) EVEN in n_y —
the linear coefficient a POVM would supply is consistent with zero, and
the entire spin dependence sits in components an affine function cannot
have. The violation is not a tail effect of marginal bins: eight
consecutive pre-registered bins exceed 5 sigma, peaking at 40 sigma where
the |n_y| = 1 hard cutoff (tau_max = 5.13) has emptied the bin that
n_y ~ 0 still fills at the percent level. And the exclusion does not
lean on error bars at all in the G3 form: two bins are exactly empty
(0/2000) for the four outermost n_y while measurably populated (>1%)
inside — an affine function vanishing at two distinct n_y values is
identically zero. The tau_max(n_y) profile (209% variation, strictly
monotone in |n_y|) shows the mechanism: the n_y-weighted spin current
advances the slow half of the guide, and the support endpoint contracts
continuously — but non-affinely — as |n_y| grows.

**Connection.** This is the missing T4-W5 step instantiated on the
corpus's own A3 observable class (Sec. VIII.F / IX.A): the spin-covariant
cutoff statistic, taken across the spin-preparation family
(experimentally: a spin rotator before the waveguide), is not realizable
by ANY fixed POVM — not merely absent from proposals. Per the operational
schema above the exclusion is conditional, exactly as VIII.F requires, on
the model's statistics being what the experiment observes. Filed as
F-T6-L7b.

## Caveats

- "Within-model and configuration-specific: this shows the model's
  first-crossing statistics over this preparation family are not
  realizable by any single POVM; it assumes the idealized non-perturbing
  screen (Prop. VIII.2's premise). A physical absorbing detector
  back-reacts and is itself POVM-describable — that dichotomy is what the
  A3 experiment adjudicates."
- Testbed grade: N = 2000 sampled trajectories per n_y, 9 family points;
  Pi-hat and tau_max are ensemble estimates, not continuum densities.
  The G3 zeros are exact at this N and at both refinements, but "zero
  arrivals" means Pi < 1.5e-3 at 95% CL, not Pi = 0.
- The reference box [-15, 45] is periodic (FFT) and cannot contain the
  fast tail to T = 16; the extended-box control shows the wrap touches
  only the late-tail order statistics at |n_y| <= 0.25 (tau_max shifts
  up to 0.19; 20 vs 24 non-arrivals at n_y = 0) and no gate quantity:
  kill-bin zeros, the >5-sigma set, and the G4 ordering are identical on
  [-15, 105].
- The bootstrap sigma is per-(bin, n_y) with shared resample indices;
  residual correlations induced by the 2-parameter fit across the 9
  points are ignored (they shrink, not inflate, the quoted
  significances' denominator only through the fit's 7 dof — with
  violations at 40 sigma and exact zeros, no gate is near this margin).
- The +/- n_y evenness is a property of this z-symmetric preparation in
  a symmetric waveguide, not of the exclusion logic; a z-asymmetric
  preparation would add an odd component without restoring affinity.
- Bins outside the populated mid-range (early bins < t = 1.2 and
  t > 13) carry small counts; they contribute no >5-sigma bins and the
  pre-registered criterion was evaluated on all 41 bins regardless.
