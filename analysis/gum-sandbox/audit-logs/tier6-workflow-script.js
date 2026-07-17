export const meta = {
  name: 'tier6-foundations-execute',
  description: 'Execute Tier-6 foundations workstreams L1/L2/L3/L5 (contextuality, ToF, leakage lemma, arrival times) in parallel',
  phases: [
    { title: 'Execute', detail: 'Four independent numerical workstreams, each with pre-registered gates' },
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
(foundations). Campaign style rules, binding:
- Work entirely inside your assigned directory (create it).
- python3 is available with numpy/scipy/matplotlib freshly installed.
- Write clean, self-contained code; run it; iterate until every gate has an
  honest verdict. Do NOT tune physics to force a PASS: if a gate fails for
  resolution reasons, refine resolution (up to ~2 retries); if it fails for
  physics reasons, report FAIL with a diagnosis. Honesty over success.
- Deliverables: the script(s), a machine-readable <name>_results.json with
  every measured number, one PNG figure, and RESULTS.md in campaign style:
  a header stating the goal, an epistemic notice line "Within-model;
  nothing here bears on nature.", a gates table (id | spec | measured |
  verdict), key numbers, method summary, and caveats.
- Units hbar = m = 1 unless stated.
- Your final structured output must faithfully mirror what RESULTS.md says.
`

phase('Execute')

const [l1, l2, l3, l5] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM L1 — Contextuality exhibit: the flipped Stern-Gerlach.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L1/

Goal: first in-repo execution of the de Broglie-Bohm guidance-sector
demonstration (Bell 1982 / Bricmont 2019) that "measurements don't
measure": same initial wave function, same initial particle position,
reversed Stern-Gerlach gradient => identical spatial trajectory, flipped
outcome LABEL.

Model (1-D, two-component spinor):
- Initial spinor Psi(z,0) = phi(z) * (|up> + |down>)/sqrt(2), phi = Gaussian
  centered 0, sigma = 1, on a periodic grid z in [-40, 40], N_grid >= 4096.
- Impulsive SG coupling: multiply components by exp(+i k z) (up) and
  exp(-i k z) (down) with kick k = 2.5 (this is the impulsive limit of
  H_int = -lambda z sigma_z). "Reversed gradient" = k -> -k.
- Free split-step FFT evolution afterward (H = p^2/2), dt <= 0.005, evolve
  until the two packets are separated by >= 8 sigma.
- Guidance: dZ/dt = j_total/rho_total with j_total = sum over components of
  Im(psi_a* dz psi_a), rho_total = sum |psi_a|^2. Integrate trajectories
  (RK4, linear interpolation of the velocity field between grid points and
  time steps). Ensemble: N = 4000 initial positions Z(0) sampled from
  |phi|^2 (fixed RNG seed).

Gates:
- G1 (nodal-line theorem): for the symmetric initial state, count median
  (z=0) crossings across all trajectories and all time steps. Spec: 0.
- G2 (outcome flip): rerun with k -> -k, identical Psi(0) and identical
  Z(0) ensemble. Spec: for EVERY trajectory the spatial exit side (sign of
  Z(T)) is unchanged, while the spin label attached to that side (which
  component dominates locally at Z(T): up-packet vs down-packet) is flipped
  for every trajectory. Report max |Z_flip(T) - Z_orig(T)| too (should be
  at machine/integration precision, since the total density and current are
  invariant under the flip - verify and state this).
- G3 (Born check): upper-exit fraction = 0.5 within 3 sigma binomial for
  N = 4000.
- G4 (determinism control): halve dt and re-integrate 100 trajectories;
  spec max position deviation <= 1e-6 * packet separation.

RESULTS.md must include one paragraph connecting the result to GUM:
this is the Gleason-Jauch/non-contextuality violation Bell says any
pilot-wave theory must exhibit; GUM Sec. III inherits it verbatim, and it
is the conceptual basis of the corpus's claim that its S3/A3 arrival-time
gate statistics are apparatus-active (Sec. VIII.F). File as finding
exhibit F-T6-2a-EXEC.`, { label: 'L1:flipped-SG', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM L2 — Momentum time-of-flight: dBB exact + Nelson leg.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L2/

Goal: execute Bricmont's Appendix-1 momentum ToF computation exactly
(deterministic de Broglie-Bohm), then the GUM-specific extension: the same
ToF statistic under NELSON stochastic dynamics (GUM Sec. III's actual
kinematics) at finite diffusion. Shows "measured momentum" is not
instantaneous momentum (which is zero), and that the result is
diffusion-blind.

Setup: Psi(x,0) = pi^(-1/4) exp(-x^2/2); free evolution has the closed
form Psi(x,t) = (1+it)^(-1/2) pi^(-1/4) exp(-x^2/(2(1+it))), giving
rho(x,t) = pi^(-1/2) (1+t^2)^(-1/2) exp(-x^2/(1+t^2)),
current velocity v(x,t) = t x/(1+t^2),
osmotic velocity u(x,t) = -2 nu x/(1+t^2) (u = nu * d/dx ln rho).
Use these closed forms; no PDE solve needed.
- dBB leg: dX/dt = v; analytic solution X(t) = X(0) sqrt(1+t^2). Integrate
  numerically (RK4, dt adaptive or <= 1e-3 early) to T = 100 and compare.
  ToF estimator p = X(T)/T. Target distribution |FT Psi|^2 =
  pi^(-1/2) exp(-p^2), i.e. Gaussian(0, sigma^2 = 1/2).
- Nelson leg: dX = (v + u) dt + sqrt(2 nu) dW, Euler-Maruyama, dt <= 1e-3
  for t < 10 then may grow geometrically; nu = 1/2 baseline; N = 20000
  particles, X(0) ~ rho(.,0), fixed seed.

Gates:
- G1: dBB trajectory law - RMS of X_num(T) - X(0) sqrt(1+T^2) over 1000
  trajectories <= 1e-6 * sqrt(1+T^2).
- G2: dBB ToF - KS test of p = X(T)/T against Gaussian(0, 1/sqrt(2)):
  p-value > 0.1 at N = 20000.
- G3: Nelson equivariance - at t in {1, 10, 100}, binned ensemble density
  vs analytic rho(x,t): report sup-norm of (binned - analytic) over bins
  with rho > 1e-3 max, spec within the ~4/sqrt(N_per_bin) statistical
  floor (state the floor you compute).
- G4: Nelson ToF - KS of X(T)/T vs Gaussian(0, 1/sqrt(2)), p > 0.1.
- G5: nu-dial - repeat Nelson ToF for nu in {0.1, 0.5, 1.0}; pairwise
  two-sample KS p > 0.05 (ToF statistic is diffusion-blind).

RESULTS.md connection paragraph: GUM Sec. III is Nelson kinematics with
hbar = 2 m nu; the corpus asserts but never exhibits that its stochastic
kinematics reproduces dBB measurement phenomenology at equilibrium. This
is that exhibit (F-T6-2b/F-T6-6-EXEC): apparatus-active momentum
"measurement" (Bell's morals), diffusion-blind at equilibrium.`, { label: 'L2:momentum-ToF', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM L3 — Lock-2 conservation lemma: the leakage-exponent ladder.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L3/

Goal: toy-grade discharge of the GUM corpus's named debt T4-W3: the
paper's Lock-2 leakage exponent (c/c_L)^5 for the stiff elliptic
constraint sector at finite c_L presupposes vanishing monopole and dipole
channels of the constraint source; the corpus never proves this. Show by
exact computation that source conservation supplies exactly that: the
radiated-power exponent ladder in 1/c_L is 1 (non-conserved monopole
source), 3 (dipole, momentum-violating), 5 (quadrupole, conserved-source
class).

Method (exact retarded solution of the 3-D scalar wave equation
(1/c_L^2) d2p/dt2 - Laplacian p = s, no PDE grid):
- Represent sources as small arrays of point sources with harmonic
  amplitude cos(omega t), array half-size a, and evaluate the EXACT
  retarded solution p(x,t) = (1/4pi) sum_i q_i cos(omega(t - |x - x_i|/c_L))
  / |x - x_i| at far-field radius R (R >= 50 * c_L,max/omega... choose R
  >> wavelength for the largest c_L used; state your choice).
- Arrays: (a) monopole: single source q = 1 at origin (total "charge"
  oscillates: non-conserved); (b) dipole: q = +1 at +a z, q = -1 at -a z
  (sum q = 0, dipole moment oscillates); (c) linear quadrupole: q = +1 at
  +2a z... use the clean pattern (+1 at +a, -2 at 0, +1 at -a) (sum q = 0
  and sum q z = 0).
- Compute time-averaged radiated power P = integral over the sphere of
  <dp/dt * (-dp/dr)> * c_L-appropriate flux ... measure instead the
  simplest correct proxy: far-field time-averaged intensity
  I = <(dp/dt)^2>/c_L integrated over the sphere at radius R (use
  Gauss-Legendre quadrature in cos(theta), symmetry in phi). Verify
  numerically that P is R-independent (two radii, <=1% relative
  difference) - this is your flux-consistency gate.
- Sweep c_L over at least 1.5 decades (e.g. c_L in logspace so that
  k a = omega a / c_L <= 0.05 at the SMALLEST c_L; pick omega = 1, a small
  accordingly). Fit slope of ln P vs ln(1/c_L).

Gates:
- G0: flux consistency - P at two radii agrees <= 1% for every case.
- G1: monopole slope = 1 +/- 0.1.
- G2: dipole slope = 3 +/- 0.15.
- G3: quadrupole slope = 5 +/- 0.25.
- G4: lemma printed - RESULTS.md states and proof-sketches the lemma: if
  the constraint source s is built from a locally conserved density
  (s = d(rho)/dt with total rho and total momentum conserved), the l = 0
  and l = 1 radiation channels cancel identically (show via the retarded
  multipole expansion: far-field amplitude at order l scales as
  (omega a/c_L)^l, and conservation kills l = 0, 1), leaving quadrupole
  leakage P proportional to (c/c_L)^5.

RESULTS.md connection paragraph: GUM VIII.D Lock 2 prints leakage
(c/c_L)^5 <= 1e-20; the T4 audit flagged the exponent as presupposing the
vanishing of monopole ((c/c_L)^1) and dipole ((c/c_L)^3) channels, with
honest fallback (c/c_L)^3. This workstream isolates the missing step as a
conservation lemma and discharges it at toy grade: IF the constraint
source is built from conserved knot mass/momentum densities (the residual
named assumption, now isolated), the printed exponent 5 is restored.
File as F-T6-5-EXEC. Also print the fallback arithmetic: with c_L/c >=
1e4 (the corpus's Bell-baseline floor), exponents 1/3/5 give leakage
1e-4 / 1e-12 / 1e-20.`, { label: 'L3:leakage-ladder', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM L5 — A3 arrival-time phenomenon on an independent engine.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier6-foundations/L5/

Goal: first in-repo instantiation of the GUM Layer-1 core phenomenon
(Das-Durr class, Sci. Rep. 9, 2242 (2019)): spin-dependent Bohmian
arrival times in a waveguide. Transverse (out-of-plane) spin => the
Pauli-current spin term modifies trajectories and the near-field
arrival-time density acquires a HARD SUPPORT CUTOFF tau_max; in-plane
spin => pure convective (flux) trajectories with no such cutoff; far
field => the two distributions converge.

Model (2-D, coordinates (x, z); spin state constant since there is no
magnetic field, so psi = phi(x,z,t) chi with chi fixed):
- Waveguide: hard walls at z = 0 and z = L (L = 1): expand in
  sin(n pi z / L) modes (sine basis, e.g. via DST or odd-extension FFT);
  initial z-profile = ground mode sin(pi z/L) (compact support
  automatically). x: TRUNCATED Gaussian, sigma_x = 1, truncated at
  +/- 3 sigma (set to zero outside, renormalize; compact support),
  centered x0 = -5, boost exp(i k0 x) with k0 = 4. Grid: x in [-15, 45],
  Nx >= 1024; Nz >= 64. Split-step: kinetic via FFT in x and sine
  transform in z (or full 2-D FFT on the odd extension in z - your
  choice, verify walls). dt <= 2e-3. Evolve to T ~ 8 (packet crosses the
  near detector).
- Pauli guidance current: j = Im(phi* grad phi) + j_spin, where for spin
  chi = eigenstate of sigma_y (out-of-plane, "transverse"):
  M_y = +rho/2 and j_spin = (-(1/2) d(rho)/dz, +(1/2) d(rho)/dx)
  [the 2-D restriction of (1/2) curl(psi-dagger sigma psi)]; for in-plane
  spin (sigma_x or sigma_z eigenstate, "axial"): j_spin contributes
  nothing in-plane, so j = Im(phi* grad phi). Velocity v = j/rho.
- Trajectories: N = 2000 initial points sampled from |phi(0)|^2 (fixed
  seed), RK4 with bilinear space interpolation and stored velocity fields
  (or on-the-fly). Arrival time = first crossing of the detector line
  x = d. Near field d_near = 1 (i.e. ~1 sigma_x past the leading edge
  region); far field d_far = 25. Record non-arrivals within T honestly.
- Run both spin cases with IDENTICAL initial ensemble.

Gates:
- G1: evolution quality - norm conserved to 1e-6; hard-wall condition
  |phi| at walls < 1e-10; energy drift <= 1e-6.
- G2: transverse near-field hard cutoff - the arrival-time density at
  d_near has empirical support ending at a measurable tau_max: ZERO
  arrivals in (tau_max, T] with T - tau_max >= 2 * (interquartile range
  of arrivals), AND tau_max stable to <= 2% under one refinement
  (double Nx or halve dt).
- G3: axial contrast - the axial-spin arrival density at d_near has
  arrivals beyond the transverse tau_max: fraction > 5% (or, if the
  contrast realizes differently, e.g. the cutoff appears in a different
  form: report the measured contrast honestly and grade PARTIAL with
  diagnosis).
- G4: far-field null - at d_far, two-sample KS between transverse and
  axial arrival-time samples: p > 0.05 (distributions converge).
- G5 (report-only, ungated): ratio of mean arrival times
  transverse/axial at d_near with bootstrap error.

Scope caveat to print verbatim in RESULTS.md: "This reproduces the
phenomenon; it neither proves nor tests POVM-exclusion (the corpus's
T4-W5 theorem gap) - it builds the testbed any such argument must
clear." Connection paragraph: this is the S3/A3 observable of GUM
Sec. IX.A (Layer 1), whose gate logic (VIII.F) the T4 audit repaired;
the hard cutoff is the non-POVM candidate statistic. File as
F-T6-L5-EXEC.`, { label: 'L5:arrival-times', phase: 'Execute', schema: RESULT_SCHEMA }),
])

return { l1, l2, l3, l5 }