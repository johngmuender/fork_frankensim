export const meta = {
  name: 'tier6-L7-povm-exclusion',
  description: 'Execute L7a/L7b POVM-exclusion testbeds for the T4-W5 gap (quadraticity + spin-affinity tests)',
  phases: [
    { title: 'Execute', detail: 'Two independent testbeds: 1-D quadraticity (exact) and spin-family affinity (L5-prime engine)' },
  ],
}

const RESULT_SCHEMA = {
  type: 'object',
  required: ['workstream', 'gates', 'key_numbers', 'files', 'summary'],
  properties: {
    workstream: { type: 'string' },
    gates: {
      type: 'array',
      items: {
        type: 'object',
        required: ['id', 'verdict', 'measured'],
        properties: {
          id: { type: 'string' },
          verdict: { type: 'string', enum: ['PASS', 'FAIL', 'PARTIAL'] },
          measured: { type: 'string' },
        },
      },
    },
    key_numbers: { type: 'object' },
    files: { type: 'array', items: { type: 'string' } },
    caveats: { type: 'array', items: { type: 'string' } },
    summary: { type: 'string' },
  },
}

const COMMON = `
You are executing one workstream of the GUM replication campaign's Tier 6
(foundations), Phase L7 — the T4-W5 POVM-exclusion testbed. Campaign
rules, binding:
- Work entirely inside your assigned directory (create it).
- python3 with numpy/scipy/matplotlib is installed.
- Write clean self-contained code; run it; iterate until every gate has an
  honest verdict. Do NOT tune physics to force a PASS; if a gate fails for
  resolution reasons refine (~2 retries); if for physics reasons, report
  FAIL with diagnosis. Honesty over success.
- Deliverables: script(s), <name>_results.json with every measured number,
  one PNG figure, RESULTS.md in campaign style (goal header; the line
  "Within-model; nothing here bears on nature."; gates table
  id | spec | measured | verdict; key numbers; method; caveats).
- Units hbar = m = 1.
- Your final structured output must faithfully mirror RESULTS.md.

THE OPERATIONAL SCHEMA (print it in RESULTS.md): a fixed physical
detector induces, via the Naimark construction, ONE
preparation-independent POVM {E(dt)} on the system Hilbert space. So over
a FAMILY of preparations psi, POVM-realizable arrival statistics must be
a quadratic (sesquilinear) functional Pi_psi(B) = <psi|E(B)|psi>. A
single preparation excludes nothing (any one distribution is trivially a
POVM); all content lives in the family. Your job: test whether the
model's Bohmian first-crossing statistics over a preparation family are
quadratic-form-realizable — a measured violation IS the exclusion step
the corpus's Sec. VIII.F gate logic is missing (T4-W5), at testbed grade,
conditional (exactly as VIII.F requires) on these statistics being what
the experiment observes. Scope caveat to print verbatim: "Within-model
and configuration-specific: this shows the model's first-crossing
statistics over this preparation family are not realizable by any single
POVM; it assumes the idealized non-perturbing screen (Prop. VIII.2's
premise). A physical absorbing detector back-reacts and is itself
POVM-describable — that dichotomy is what the A3 experiment adjudicates."
`

phase('Execute')

const [l7a, l7b] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM L7a — quadraticity/polarization test (1-D spatial family, exact).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L7a/

Model (all closed-form; no PDE solve):
- Components: psi_j(x,t) = free evolution of Gaussian envelope
  exp(-(x-x0)^2/(2 sigma^2)) * exp(i k_j x), sigma = 1, x0 = -8, with
  k1 = 1.2, k2 = 2.4. Use the analytic free-Gaussian propagation formula
  for each component (a boosted Gaussian stays a closed-form Gaussian:
  center x0 + k t, spreading (1 + i t/sigma^2)-type factor — derive
  carefully and verify norm conservation to 1e-12 by quadrature).
- Family: psi_thetaphi = cos(theta) psi_1 + sin(theta) e^{i phi} psi_2,
  automatically normalized iff <psi_1|psi_2> = 0 — it is NOT exactly zero
  (overlap exp-suppressed in (k2-k1)); compute the exact overlap, print
  it, and NORMALIZE each family member explicitly (the POVM form
  <psi|E|psi>/<psi|psi> for unnormalized psi is a ratio of quadratics —
  handle by normalizing states; note this in RESULTS.md).
- Detector at d = 0; horizon T = 12; time grid dt <= 0.005.
- N_psi(s) = integral_0^inf |psi_s|^2 dx by high-order quadrature on a
  fine grid (verify against 1 - integral_{-inf}^0).
- EXACT Bohmian first-arrival CDF (1-D no-crossing theorem):
  CDF_psi(t) = max_{s<=t} N_psi(s). State + cite the ordering argument in
  RESULTS.md (trajectories cannot cross in 1-D; the ever-crossed set is a
  right-interval whose mass is the running max of the presence
  probability N).
- Grid: theta in {0, pi/12, ..., pi/2} (7 values), phi in {0, pi/4, ...,
  7pi/4} (8 values) — 56 preparations (well-determined: 4 fit params/bin).
- Bins: ~48 uniform time bins over [0, T] in CDF increments (arrival
  probability per bin), plus the non-arrival mass at T.
- Per-bin least-squares fit of the POVM form
  F(theta,phi) = a cos^2(theta) + b sin^2(theta)
               + sin(2 theta) (c cos(phi) + d sin(phi));
  residual = max |data - fit| per bin, and global max.

Gates:
- G1 (numerics + running-max control): component norms conserved 1e-10;
  at one backflow-heavy family point (find one where N(s) is
  non-monotone; report its backflow magnitude max_s [max_{s'<=s} N(s') -
  N(s)]), integrate an explicit trajectory fan (RK4, >= 4000 initial
  points sampled from |psi_0|^2) and compare its empirical first-crossing
  CDF to the running-max CDF: sup-norm <= 2 x Dvoretzky-Kiefer-Wolfowitz
  95% band. Also verify no trajectory crossing (ordering preserved) in
  the fan.
- G2 (comparator floor): fit the SAME POVM form to the manifestly
  quadratic functional N_psi(t_i) (presence probability at bin edges,
  computed from normalized states): report max residual = the harness
  floor. Spec: floor <= 1e-6.
- G3 (exclusion, pre-registered): Bohmian per-bin residual > 1e3 x the G2
  floor AND > 1e-4 absolute, in at least one bin where backflow is
  present; report the full residual-vs-bin profile and which bins carry
  backflow.
- G4 (no-backflow control): sub-family with k1 = k2 = 1.8 and envelope
  centers x0 = -8, -11 (positions differing, same momentum — verify
  numerically that N(s) is monotone for all family members i.e. zero
  backflow; if not, adjust separation once and print the adjustment):
  Bohmian residual at the G2 floor (<= 10x floor). This shows the G3
  violation is backflow-borne, matching the flux = first-crossing theorem
  for current-positive dynamics.

RESULTS.md connection paragraph: this is the polarization-identity leg of
the T4-W5 exclusion — Bohmian arrival statistics over spatial
superposition families fail the quadratic-form constraint exactly where
backflow lives; Kijowski-class and quantum-flux proposals are quadratic
by construction and cannot reproduce them. File as F-T6-L7a.`, { label: 'L7a:quadraticity', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM L7b — spin-family affinity test (the A3-relevant leg).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L7b/

CONTEXT — read these first:
- /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L5prime/RESULTS.md
- /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L5prime/l5prime_run.py
  (REUSE this validated engine — exact-in-z factorized evolution, C2
  smoothstep taper, k0 = 2, hard-wall waveguide z in [0,1], ground sine
  mode, box x in [-15,45], reference Nx = 2048 acceptable here, dt <=
  5e-4 reference; N = 2000 trajectories, seed 20260717, d_near = 1,
  T = 16.)

Physics: spin state chi(n-hat) on the Bloch sphere, NO magnetic field, so
psi = phi(x,z,t) chi with phi spin-independent: ONE field evolution
total. The Pauli guidance current is j = Im(phi* grad phi) +
(n_y/2)(-d(rho)/dz, +d(rho)/dx) — the spin term enters ONLY through n_y
(out-of-plane Bloch component). So the velocity field is
v(x,z,t; n_y) = a(x,z,t) + n_y b(x,z,t): precompute/store the two vector
fields once, then run the trajectory ensemble per n_y cheaply.
KEY REDUCTION (print in RESULTS.md): a POVM gives per bin
Pi_B(n-hat) = A + B . n-hat (affine on the Bloch sphere, from the 2x2
effective operator <phi|E_B|phi>); pure states with the same n_y but
different (n_x, n_z) have IDENTICAL dynamics, forcing B_x = B_z = 0;
hence Pi_B must be AFFINE IN n_y: Pi_B(n_y) = A + B n_y. Any measured
nonlinearity in n_y, or any interval-vanishing-with-positivity-elsewhere
pattern, excludes every POVM for this family.

Runs: n_y in {-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1} (9 values),
identical initial ensemble (N = 2000, seed 20260717) for every n_y.
Record first-crossing times at d_near = 1 (and non-arrivals by T = 16).

Gates:
- G1 (engine): norm 1e-6; wall amplitude < 1e-10; energy drift 1e-6
  (expect ~1e-16 from the L5prime engine); cross-check: the n_y = +1 and
  n_y = 0 runs must reproduce L5prime's transverse/axial results
  (tau_max within 1e-3 relative; KS p > 0.5 between arrival samples).
- G2 (linearity test): bin arrivals (~40 bins over [0, 16] chosen so
  expected counts >= 20 in the populated mid-range, plus the non-arrival
  bin); per bin, least-squares fit Pi_B(n_y) = A + B n_y over the 9
  points; residuals compared to per-bin bootstrap sigma (2000 resamples).
  Exclusion criterion (pre-registered): >= 3 independent bins with max
  residual > 5 sigma_boot. Report the residual/sigma profile for all
  bins. ALSO run the same fit with a quadratic term A + B n_y + C n_y^2
  and report the C significances — diagnostic of the violation shape
  (not gated).
- G3 (affine-vanishing kill): measure tau_max(n_y) for each n_y (last
  arrival + the empty-gap check as in L5prime). Exhibit >= 1 time bin
  B with ZERO arrivals (0/2000, i.e. Pi_B < 1.9e-3 at 95% CL) for every
  n_y in an interval at one end (e.g. |n_y| >= n*) while Pi_B > 1e-2
  measured for some n_y outside it. Print the lemma: an affine function
  of n_y that vanishes on an interval of n_y values is identically zero,
  contradicting the measured positivity — no single POVM reproduces the
  family. Verify the zero-arrival claim for the critical bins at one
  refined resolution (double Nx or halve dt).
- G4 (structure, report + gate on continuity): tau_max(n_y) profile —
  spec: varies by > 20% across the family and is monotone in |n_y| or
  reported otherwise with diagnosis.
- G5 (report-only): effect sizes — max linearity violation in sigma
  units; the n_y-dependence of mean arrival time.

RESULTS.md connection paragraph: this is the missing T4-W5 step
instantiated on the corpus's own A3 observable class (Sec. VIII.F/IX.A):
the spin-covariant cutoff statistic, taken across the spin-preparation
family (experimentally: a spin rotator before the waveguide), is not
realizable by ANY fixed POVM — not merely absent from proposals. File as
F-T6-L7b.`, { label: 'L7b:spin-affinity', phase: 'Execute', schema: RESULT_SCHEMA }),
])

return { l7a, l7b }