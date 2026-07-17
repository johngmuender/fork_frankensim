# F-T6-L7a — POVM-exclusion testbed: quadraticity/polarization test of Bohmian first-arrival statistics (1-D spatial family, exact)

**Goal.** Test whether the model's exact Bohmian first-crossing statistics over a
two-component preparation family are realizable as a quadratic (sesquilinear)
functional of the state — i.e., by any single preparation-independent POVM. A
measured violation is the exclusion step the corpus's Sec. VIII.F gate logic is
missing (T4-W5), at testbed grade, conditional (exactly as VIII.F requires) on
these statistics being what the experiment observes.

Within-model; nothing here bears on nature.

## The operational schema

A fixed physical detector induces, via the Naimark construction, ONE
preparation-independent POVM {E(dt)} on the system Hilbert space. So over a
FAMILY of preparations psi, POVM-realizable arrival statistics must be a
quadratic (sesquilinear) functional Pi_psi(B) = <psi|E(B)|psi>. A single
preparation excludes nothing (any one distribution is trivially a POVM); all
content lives in the family.

## Gates

| id | spec | measured | verdict |
|----|------|----------|---------|
| G1 | component norms conserved <= 1e-10; at a backflow-heavy family point, 4000-trajectory RK4 fan first-crossing ECDF within 2x DKW-95 band of the running-max CDF; no trajectory crossing | norm dev 4.4e-16; star prep (theta,phi) = (0.229pi, 0.0625pi), backflow 2.387e-4; sup-norm ECDF vs running-max CDF = 1.250e-4 (= 1/2n stratification limit) vs band 4.295e-2; ordering inversions: 0 (min gap 5.6e-4 > 0). Note: the backflow-heavy point exists only on the refined grid — the registered 7x8 grid contains none | PASS |
| G2 | comparator floor: same POVM form fitted to the manifestly quadratic N_psi(t_i); floor <= 1e-6 | floor = 2.22e-15 (registered 7x8 grid), 2.89e-15 (refined 25x32 grid) | PASS |
| G3 | Bohmian per-bin residual > 1e3 x floor AND > 1e-4 absolute in >= 1 backflow bin | registered 7x8 grid: **no exclusion** (zero backflow on grid, residual 3.9e-16 = floor). After sanctioned resolution refinement of the preparation grid (13x16: residual 4.2e-5, below threshold; 25x32: residual **1.843e-4** in backflow bins 15-17, t in [3.75, 4.5]; threshold max(1e3 x 2.89e-15, 1e-4) = 1e-4) | PASS (after grid refinement; registered sampling alone excludes nothing) |
| G4 | no-backflow control (k1 = k2 = 1.8, centers -8/-11): N monotone for all family members (one separation adjustment allowed); Bohmian residual <= 10x floor | centers -8/-11: NOT monotone (grid backflow 4.27e-4, residual 3.71e-4 — itself a clean backflow-borne violation). **Adjusted once: -11 -> -9.5** (separation 1.5): all 56 family members monotone (runmax - N = 0.0 in FP); residual 3.72e-16 <= 10x floor = 2.22e-14 | PASS (after printed adjustment; continuum caveat below) |

## Key numbers

- Overlap <psi_1|psi_2> = exp(-sigma^2(k2-k1)^2/4) e^{i(k2-k1)x0} =
  -0.68699 + 0.12162i, |.| = 0.6977 — **not** exponentially negligible.
  Numeric matches analytic to 3.6e-15; constant in t to 3.6e-15
  (unitarity check). Family normalizations n^2 in [0.313, 1.687].
- Component norm conservation (quadrature, composite 12-pt Gauss-Legendre,
  panels of 0.5 over |x| <= 160): max |∫|psi_j|^2 - 1| = 4.4e-16 over all t
  (spec 1e-12/1e-10). Right-half integral verified against
  1 - ∫_{-inf}^0 via the same closure (identity to machine precision).
- Registered 7x8 (theta,phi) grid, 56 preparations: max backflow
  max_s[max_{s'<=s}N(s') - N(s)] = 3.9e-31 (i.e. zero); Bohmian per-bin
  residual 3.92e-16 — quadratic at floor. **No exclusion at the registered
  sampling.**
- Continuum diagnosis (2x2 detector-current matrix J(t), j(0,t) = c†J(t)c):
  lambda_min(J(t)) = -4.92e-4 at t = 3.82 < 0 — negative-current states
  exist; fine 121x128 scan: backflow pocket max 2.83e-4 at
  (theta, phi) = (0.772, 0.245), 221/15488 grid points > 1e-5. The 7x8 grid
  misses the pocket entirely -> G3's registered failure is a resolution
  failure, refined per campaign rules (2 retries: 13x16, then 25x32).
- Refined 25x32 grid (800 preparations, still 4 fit params/bin): max
  backflow 2.387e-4; per-bin residual profile peaks at 1.843e-4 (bin 16,
  t in [4.0, 4.25]); residual > 1e-4 in bins 15-17; backflow-contaminated
  bins (running max pinned above N at a bin edge for some prep):
  10-18, 23-48. Max residual in **zero-backflow bins: 5.3e-16** (= floor).
- G1 fan (RK4, dt = 0.001, 4000 stratified-quantile initial points from
  |psi_0|^2): sup|ECDF - runmax CDF| = 1.250e-4 = 1/(2n); the fan resolves
  the backflow signal (sup|ECDF - N| = 1.54e-4; CDF - N peaks at 2.39e-4):
  trajectories follow the running max, not N. Zero ordering violations in
  12000 steps. A second fan at the control family's backflow point
  (sep 3.0, dt = 0.001): sup-norm 1.249e-4, 0 inversions — running-max CDF
  validated where backflow is genuinely present.
- Continuity cross-check dN/dt vs closed-form j(0,t): max dev 4.7e-6
  (finite-difference truncation limited); min control current -2.7e-3.
- G4 control: sep 3.0 grid backflow 4.27e-4, residual 3.71e-4 confined to
  its backflow bins {17-19, 37-39}, 2.7e-16 elsewhere; separation scan
  (continuum, fine grid): max backflow 3.1e-4 / 3.3e-4 / 3.7e-4 / 4.8e-4 /
  6.4e-4 at sep 1.0/1.5/2.0/3.0/4.0 — **no separation in [1,4] gives a
  continuum-backflow-free family**; after the sanctioned single adjustment
  to sep 1.5 the 56 registered family members are exactly monotone and the
  residual sits at 3.72e-16.
- Backflow-borne correlation (both families): residual at floor (<= 5.3e-16)
  in every bin without backflow contamination; violation (up to 1.8e-4 /
  3.7e-4) exclusively in backflow bins.

## Method

- Components: psi_j = free evolution of (pi sigma^2)^{-1/4}
  exp(-(x-x0)^2/(2 sigma^2)) e^{i k_j x}, sigma = 1, x0 = -8, k1 = 1.2,
  k2 = 2.4; closed form psi_j(x,t) = (pi sigma^2)^{-1/4} sqrt(sigma^2/w)
  exp(-(x - x0 - k_j t)^2/(2w) + i k_j x - i k_j^2 t/2), w = sigma^2 + it
  (derived by Fourier transform + exact Gaussian integration; verified by
  the norm/overlap unitarity checks above). Detector d = 0, horizon T = 12,
  dt = 0.005.
- Family psi_thetaphi = cos(theta) psi_1 + sin(theta) e^{i phi} psi_2,
  **explicitly normalized** because <psi_1|psi_2> = 0.698 != 0. For
  unnormalized psi the POVM prediction <psi|E|psi>/<psi|psi> is a ratio of
  quadratics; normalizing makes the test exact: data(theta,phi) must equal
  [a cos^2 theta + b sin^2 theta + sin 2theta (c cos phi + d sin phi)] /
  n^2(theta,phi), which is linear in (a,b,c,d) once the basis functions are
  divided by n^2 — this is the least-squares design used (4 params/bin;
  the 2x2 restriction of any Hermitian E(B) to span{psi_1, psi_2} has
  exactly these 4 real parameters, so the fit is fully general).
- N_psi(s) = ∫_0^inf |psi_s|^2 dx assembled exactly from the three pair
  integrals ∫_0^inf psi_i^* psi_j dx (composite Gauss-Legendre, 12-pt
  panels of width 0.5 to x = 160), cross-checked against the closed-form
  erfc expression for the diagonal terms and against 1 - ∫_{-inf}^0.
- Exact first-arrival CDF (1-D no-crossing / running-max theorem): Bohmian
  trajectories solve dx/dt = Im(psi'/psi), a first-order ODE with a single-
  valued velocity field, so in 1-D trajectories cannot cross and positions
  at equal times preserve the initial ordering. Hence the set of initial
  conditions whose trajectory has ever crossed d by time t is a right
  interval in the ordering, and its mass — the first-crossing CDF — equals
  the running maximum of the presence probability:
  CDF_psi(t) = max_{s<=t} N_psi(s). (Ordering argument as in Leavens's
  1-D arrival-time analyses; Daumer-Dürr-Goldstein-Zanghì, J. Stat. Phys.
  88, 967 (1997); Vona, Hinrichs, Dürr, PRL 111, 220404 (2013). For
  current-positive dynamics this reduces to the quantum-flux/Kijowski
  arrival distribution — the content of G4.) Independently verified here by
  the explicit RK4 trajectory fan (G1), including at a genuinely
  backflow-carrying preparation.
- Bins: 48 uniform CDF-increment bins on [0, 12] plus the non-arrival mass
  at T (bin 48). Per-bin max-abs residual of the least-squares POVM-form
  fit; comparator floor from the same fit applied to the manifestly
  quadratic N_psi(t_i) (G2).
- Backflow-bin flag: a bin's increment is contaminated iff the running max
  is pinned above N at either bin edge for some preparation (an interior
  dip that recovers before the edge leaves the increment quadratic).
- Diagnostics: detector-current matrix J_jl(t) = [psi_j^* psi_l' -
  psi_j'^* psi_l]/(2i) at x = 0; lambda_min(J(t)) >= 0 for all t would
  imply zero backflow for the entire continuum family. Fine 121x128
  (theta,phi) scans map the backflow pocket. Physics of the pocket: the
  components' local phase-gradient mismatch at the detector is
  dk_loc(t) = (k2-k1) sigma^4/(sigma^4+t^2) (0.099 at t = 3.33, 0.046 at
  t = 5 — an order below the mean local wavenumber ~1.6), because
  co-located packets' local momenta converge to (x-x0)/t; backflow survives
  only in a narrow (theta,phi) pocket via the envelope-gradient
  interference term. The displaced-center control family keeps a
  center-offset-driven mismatch, which is why it carries robust backflow
  near destructive phi ~ pi at every separation scanned.

## Connection (T4-W5)

This is the polarization-identity leg of the T4-W5 exclusion: Bohmian
first-arrival statistics over spatial superposition families fail the
quadratic-form (single-POVM) constraint exactly where backflow lives — the
per-bin residual rises from the 1e-15 comparator floor to 1.8e-4 precisely
in the bins where the running-max CDF detaches from the presence
probability, and returns to the floor in every backflow-free bin of both
families. Kijowski-class and quantum-flux proposals are quadratic in the
state by construction and therefore cannot reproduce these statistics; the
measured violation is the operational exclusion step Sec. VIII.F's gate
logic is missing, at testbed grade, conditional on the model's
first-crossing statistics being what the experiment observes. File as
F-T6-L7a.

## Caveats

- "Within-model and configuration-specific: this shows the model's
  first-crossing statistics over this preparation family are not realizable
  by any single POVM; it assumes the idealized non-perturbing screen
  (Prop. VIII.2's premise). A physical absorbing detector back-reacts and
  is itself POVM-describable — that dichotomy is what the A3 experiment
  adjudicates."
- Pre-registration deviation (G3): the registered 7x8 (theta,phi) grid
  contains no backflow-carrying preparation — at that sampling the
  statistics are quadratic to 3.9e-16 and exclude nothing. The violation
  was recovered only after refining the preparation grid (13x16, then
  25x32) of the *same* registered family; physics (k's, sigma, x0, d, T,
  bins) untouched. The 13x16 refinement shows a violation (4.2e-5 > 1e3 x
  floor) that misses the 1e-4 absolute threshold; the exclusion magnitude
  is grid-dependent because the backflow pocket is narrow (~1.4% of the
  (theta,phi) torus above 1e-5).
- G4's design premise ("same momentum, different centers => no backflow")
  is not exact: the continuum family carries a backflow pocket
  (3.1e-4-6.4e-4) at every separation in [1, 4]. Monotonicity after the
  sanctioned adjustment (sep 1.5) holds for the 56 registered family
  members (exactly, in floating point) — which is what the fit consumes —
  but not for the continuum family; symmetrically, the sep-3.0 control run
  itself furnishes a second, clean backflow-borne violation. The
  backflow-borne conclusion rests on the bin-level on/off correlation in
  both families, which is exact (floor vs 1e-4-scale).
- Exclusion magnitude (1.8e-4) is bounded by the model's backflow
  (2.4e-4 on the 25x32 grid; 2.8e-4 continuum max): a real experiment
  would need arrival-probability resolution at the few-1e-4 level per
  0.25-wide time bin for this configuration.
- The fan's DKW comparison uses stratified quantile sampling, whose
  deviation (1/2n = 1.25e-4) is far inside the iid DKW band (4.3e-2); the
  band is therefore a loose bound, and the real agreement is at the
  stratification limit.
- The G2 floor (2.2e-15) is a same-quadrature comparator: the fitted data
  and the design share the pair-integral quadrature, so quadrature error
  cancels in the floor; absolute quadrature accuracy is separately
  bounded by the unitarity checks (4.4e-16) and the continuity check
  (4.7e-6, finite-difference limited).

## Files

- `l7a_povm_exclusion.py` — closed-form model, registered pipeline (stage 1)
- `l7a_diagnostics.py` — current-matrix lambda_min scans, separation scan,
  control fan, continuity check
- `l7a_refine.py` — pocket map, grid refinements, final gate logic
- `l7a_figure.py` — figure + merged `L7a_results.json`
- `L7a_results.json` — every measured number
- `L7a_povm_exclusion.png` — 4-panel summary figure
