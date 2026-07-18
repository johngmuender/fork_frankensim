export const meta = {
  name: 'tier7-program-execute',
  description: 'Execute Tier-7 program workstreams M1-M4 (w(z) tracker dynamics, servo coefficient, F-R9 saturation pricing, Nelson relaxation)',
  phases: [
    { title: 'Execute', detail: 'Four independent workstreams from ROADMAP_v9_PROGRAM with pre-registered gates' },
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
          verdict: { type: 'string', enum: ['PASS', 'FAIL', 'PARTIAL', 'PASS-MATCH', 'PASS-MEASURED', 'CONDITIONAL'] },
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
You are executing one workstream of the GUM replication campaign's
Tier 7 (program extensions). Campaign rules, binding:
- Work entirely inside your assigned directory (create it).
- python3 with numpy/scipy/matplotlib is installed.
- The full workstream spec INCLUDING GATES is in
  /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/ROADMAP_v9_PROGRAM.md
  — read it first; your prompt below restates the essentials but the
  roadmap text is authoritative on gate wording.
- Write clean self-contained code; run it; iterate until every gate has
  an honest verdict. Do NOT tune physics to force a PASS; resolution
  refinements allowed (~2 retries); physics failures reported as FAIL
  with diagnosis. Honesty over success.
- Deliverables: script(s), <name>_results.json with every measured
  number, one PNG figure, RESULTS.md in campaign style (goal header;
  "Within-model; nothing here bears on nature."; gates table
  id | spec | measured | verdict; key numbers; method; caveats).
- Units hbar = m = 1 unless stated.
- Your final structured output must faithfully mirror RESULTS.md.
`

phase('Execute')

const [m1, m2, m3, m4] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM M1 — WS4-M4c toy: w(z) from the Gamma(H) relaxation profile.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M1/

Context: GUM VI.D's dark-energy mechanism — a slow mode q displaced from
its minimum by Hubble drag, rho_Lambda = (1/2) chi^{-1} dq^2, with the
corpus's printed claims: exact tracker identity w(z) = -1 + Omega_m(z)
(pure tracker, W = 1); memory forms beta = -0.3 -> W: 1.28->1.03 and
beta = -1.0 -> W: 1.72->1.05 (W defined by w(z) = -1 + W*Omega_m(z),
evaluated over z = 0.5 -> 0); CPL projection w_a = +3 W Omega_m0 (1 -
Omega_m0) (~ +0.63 at W = 1, Omega_m0 = 0.3, always positive).

Model (derive before coding; print derivations in RESULTS.md):
- Flat FRW background: for gates use the test-field approximation
  (background = LCDM, Omega_m0 = 0.3, Omega_L0 = 0.7) AND check
  self-consistency by one iteration (feed rho_Lambda(z) back into H);
  report the difference.
- Pure tracker: the corpus's residual is rho_Lambda,eff = alpha H^2
  M*^2 — implement rho(z) proportional to H^2(z) and compute
  w(z) = -1 + (1/3) dln rho/dln(1+z). Verify ANALYTICALLY +
  numerically the exact identity: in a flat matter+Lambda universe,
  dln H^2/dln(1+z) = 3 Omega_m(z), hence w(z) = -1 + Omega_m(z)
  EXACTLY (show the two-line derivation).
- Dynamics leg: also implement the relaxation ODE
  d(dq)/dt = -Gamma dq + S H with Gamma = const in the fast-relaxation
  (tracking) regime, verify the quasi-static attractor gives
  dq proportional to H hence rho proportional to H^2 (the ODE
  underwrites the algebraic form; measure the attractor approach rate).
- Memory forms: use the corpus's own printed form
  rho_Lambda,eff = alpha H^2 M*^2 + beta Hdot M*^2 directly. Test the
  reading beta = printed value with alpha = 1 (i.e. beta/alpha = -0.3,
  -1.0): compute W(z) = (w(z)+1)/Omega_m(z) at z = 0.5 and z = 0 and
  compare to the printed pairs (1.28->1.03; 1.72->1.05). If that
  reading fails, try alternative natural normalizations and report
  which (if any) matches — this is the G4 convention-recovery task; an
  honest no-match is a legitimate outcome with diagnosis.

Gates (roadmap authoritative): G1 numerics + rho_Lambda positivity
(the q-variable sign theorem's content); G2 exact tracker identity to
1e-6 over z in [0,10] + printed derivation; G3 CPL fit around a = 0.7:
w_a = +3 W Omega_m0(1-Omega_m0) within 5% at W = 1 and w_a sign
STRICTLY positive across the whole family tested; G4 memory-value
recovery vs printed pairs within 15% (else FAIL with convention
diagnosis); G5 Branch-A operational content: (w+1) proportional to
Omega_m(z) exactly in the tracker limit; deviation shape with memory
reported.
RESULTS.md connection: first independent instantiation of WS4-M4c
(Sec. XI open #4) at toy grade; underwrites the F-R12 strengthening
(exact identity + always-positive w_a) on dynamics rather than
algebra. File as F-T7-M1.`, { label: 'M1:wz-tracker', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM M2 — the bootstrap servo coefficient (r2 toy law).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M2/

Context: GUM IV.H.3: each knot is a phase-locked self-oscillator of the
medium; displacing the vacuum order parameter hbar = h*(1+delta) gives
clock mismatch with EXACT exponents d ln omega_mech/d ln h = +1/2,
d ln omega_phase/d ln h = -1/2 (h27-confirmed), locking torque
proportional to -sin(Dphi), and the corpus's UNREPLICATED toy law:
relaxation rate lambda = (0.021 +/- 0.004) g_chi^2 omega* in the
thermalizing-bath class; marginal (no relaxation) in the integrable
limit (h27: gamma = 0 exactly).

Model: build the minimal faithful toy and DERIVE your reduced equations
in RESULTS.md before presenting numbers. Recommended architecture:
- Pair (delta, Dphi): Dphi_dot = omega* [(1+delta)^{1/2} -
  (1+delta)^{-1/2}] (exact exponents built in; ~ omega* delta small-
  delta), delta_dot = -g_chi^2 A sin(Dphi) + bath coupling. The
  deterministic pair is conservative-oscillatory (marginal); the BATH
  supplies the irreversibility that converts the beat into secular
  relaxation (the corpus's 'thermalizing-bath class').
- Bath: N >= 100 harmonic modes bilinearly coupled (to Dphi's conjugate
  or to delta — your choice, state it), Ohmic spectral density
  J(omega) = eta omega exp(-omega/omega_c), omega_c = 5 omega*, small
  temperature. Average over >= 20 realizations.
- Measure: envelope decay rate lambda of delta(t) from delta0 = 0.05.

Gates (roadmap authoritative): G1 integrator quality + bath statistics
verification; G2 lambda proportional to g_chi^2 across >= 1 decade
(log-slope 2.0 +/- 0.1); G3 lambda/omega* invariant under omega*
rescaling x4 (+/- 10%); G4 coefficient c = lambda/(g_chi^2 omega*)
with error vs 0.021 +/- 0.004 — within 2 sigma_combined = PASS-MATCH;
else PASS-MEASURED with the honest caveat (r2's bath spec unreleased;
coefficient bath-dependent) and c reported for >= 2 bath spectral
shapes (Ohmic + super-Ohmic J proportional to omega^3) to bound the
convention spread; G5 integrable-limit control: bath off => marginal
(fitted |lambda| consistent with 0 at the fit-noise floor) — h27's
gamma = 0 class boundary reproduced.
RESULTS.md connection: the last unreplicated [CAL] anchor of the
quantum-adjacent core (IV.H.3 r2); h27 covered exponents/torque/
boundary — this covers the LAW. File as F-T7-M2.`, { label: 'M2:servo-coefficient', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM M3 — F-R9 route (b) priced: on-tube saturation nonlinearity.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M3/

MANDATORY CONTEXT (gate condition, not suggestion): read
/home/user/fork_frankensim/analysis/gum-sandbox/theory-audit/h21_RESULTS.md
and /home/user/fork_frankensim/analysis/gum-sandbox/theory-audit/h21_mc.py
first, and REUSE h21's parametrization of the line-supported web
(Campbell/shot-noise statistics of the amplitude A at a point, tube
density/geometry parameters, the scale window that produced
R = <A^2>/<A>^2 = 1e14-1e18, and the mean-rate check 1.000 +/- 0.001).
Do not invent a new parametrization.

Context: F-R9 (substantive): NR-K3a's 0nubb-insert rate claim drops the
second moment; rate is quadratic in A; the line-supported web gives
R = <A^2>/<A>^2 = 1e14-1e18 across the corpus's own scale window; the
motional-narrowing rescue is causally excluded (5e7-1.2e9 c). Discharge
route (b), never priced: an on-tube saturation nonlinearity. Your job:
price it.

Method:
- Reproduce h21's unsaturated R at >= 1 anchor point of its scale
  window (agreement within x3) and the mean check (1.000-class).
- Saturation model: A -> A_sat tanh(A/A_sat). Scan A_sat/<A> over >= 4
  decades (log grid). For each: R(A_sat) = <A_s^2>/<A_s>^2 and the
  mean-rate distortion D(A_sat) = <A_s^2>/<A^2-based claimed mean rate>
  — use h21's mean-claim definition and state it precisely. The point:
  the rescue must tame the fluctuation excess WITHOUT destroying the
  mean-rate phenomenology the corpus claims.
- Determine A_sat* where R <= 10; report D(A_sat*).
- Verdict per the pre-registered decision rule (roadmap G3): VIABLE iff
  |D-1| <= 0.5 at A_sat* AND A_sat* is inside h21's on-tube amplitude
  range; if the range is not extractable from h21's parametrization,
  grade CONDITIONAL and print the required range as the isolated
  assumption (the L3-style residue).
- G4: print the lemma-level account: saturation bounds <A_s^2> <=
  A_sat^2 so R is tamed at scale A_sat; the competition that makes
  route (b) hard: on a sparse line-supported web the SECOND moment is
  dominated by rare close-tube configurations while the mean is not —
  quantify the trade in your scan.

Gates: G1 h21 reproduction (anchor R within x3; mean 1.000-class);
G2 pricing curves R(A_sat), D(A_sat) measured over >= 4 decades;
G3 verdict per decision rule (VIABLE / NOT-VIABLE / CONDITIONAL with
the isolated assumption printed); G4 lemma + trade quantified.
RESULTS.md connection: F-R9's route (b) moves from 'nowhere in the
text' to priced — whichever way it lands, the 0nubb insert's repair
menu gets one concrete entry or one closed door. File as F-T7-M3.`, { label: 'M3:saturation-pricing', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM M4 — Nelson non-equilibrium relaxation: tau(nu).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M4/

MANDATORY CONTEXT (gate condition): read
/home/user/fork_frankensim/analysis/gum-sandbox/tier3-born/REPORT.md
first and REUSE its conventions: the 2-D box setup, mode sets, the
non-equilibrium initial density rho_0 choice, the coarse-grained
H-function estimator H(t) = integral rho_bar ln(rho_bar/|psi|^2)
(binning included), and its tau values (M = 9: tau = 10.2; M = 16:
tau = 4.9) for the dBB control comparison. If tier3's exact phase or
rho_0 conventions are not fully recoverable from the REPORT, use fixed-
seed equivalents and state the substitution.

Context: GUM Sec. III's actual kinematics is Nelson diffusion at
nu = hbar/2m, but every rates anchor (V-SIM-4, tier3) used
deterministic dBB guidance. At equilibrium the statistics are nu-blind
(proved in tier6 L2). Out of equilibrium they cannot be: the osmotic
term drives rho toward |psi|^2 directly. Pre-registered direction:
tau(nu) strictly decreasing.

Model: 2-D box, psi = equal-weight superposition of the first M modes
(M in {9, 16}; skip M = 4), exact psi(t) by mode sum (hard-wall box:
sine modes). Particles: N = 10,000, X(0) from tier3's rho_0 convention.
Dynamics: dX = (v + u)dt + sqrt(2 nu)dW with v = grad S (from the mode
sum), u = nu grad ln|psi|^2; nu in {0 (dBB control), 0.1, 0.5, 1.0};
reflecting boundary for the SDE at the walls (state your treatment; at
the walls |psi|^2 -> 0 so the osmotic term repels — verify no-flux
numerically). dt small enough that halving changes tau by < 5%.
Measure H(t) with tier3's binning; fit tau in the near-exponential
window.

Gates (roadmap authoritative): G1 dBB control within 25% of tier3's
tau (state comparison caveats honestly); G2 equivariance control at
every nu (equilibrium-born ensemble stays at floor); G3 monotonicity:
tau(nu) strictly decreasing in nu for both M (pre-registered; if
non-monotone, FAIL + diagnosis); G4 quantification:
tau(0.5)/tau(dBB) both M with bootstrap errors; continuity check
tau(0.1)/tau(0) >= 0.7; G5 consequence paragraph (report-only): which
corpus claims move which way (laboratory equilibrium quality:
strengthened; relic non-equilibrium survival window: weakened — print
the honest double edge).
RESULTS.md connection: first nu-dependence measurement of the Born
H-theorem rates under GUM's actual kinematics; the corpus's III.D
[CAL] rates become direction-signed conservative statements. File as
F-T7-M4.`, { label: 'M4:nelson-relaxation', phase: 'Execute', schema: RESULT_SCHEMA }),
])

return { m1, m2, m3, m4 }