export const meta = {
  name: 'tier8-S-execute',
  description: 'Execute ROADMAP v10 workstreams S1-S4 (alternatives + extensions)',
  phases: [
    { title: 'Execute', detail: 'four parallel gated workstreams' },
  ],
}

const ROOT = '/home/user/fork_frankensim/analysis/gum-sandbox'
const T8 = `${ROOT}/tier8-assessment`

const DISC = `EXECUTION DISCIPLINE (binding, from ROADMAP_v10_ALTERNATIVES.md — read it first at ${T8}/ROADMAP_v10_ALTERNATIVES.md, your workstream section is the spec; gates are PRE-REGISTERED and may not be moved):
- Within-model; nothing bears on nature. Honest FAILs are findings. Grade language never rises.
- Work in your workstream directory (create it). Deliver: the script(s), machine-readable results JSON, a figure (PNG via matplotlib, Agg backend) where meaningful, and RESULTS.md in the campaign house style (header with date/phase/code/outputs, epistemic frame, numbered sections, a gates table with requirement/measured/verdict, key-numbers table, honest caveats section).
- ANY long-running computation (> ~10 min) MUST be launched with nohup and MUST checkpoint per-phase/per-item to JSON/npz files so a finisher agent can resume after a container restart without recomputation (the P2 pattern: script takes a phase/item CLI arg, reads committed intermediates). Never leave a long run un-nohup'd.
- If the pre-registered spec proves defective or too expensive, do NOT silently change it: implement the closest feasible variant, and PRINT the deviation in RESULTS.md as a coordinator-owned spec amendment (this is the campaign's standing rule; deviations are findings).
- Budget your wall-clock: target the whole workstream <= ~2.5 hours. Benchmark unit cost FIRST, then scale the run to fit, printing any scale-down as a spec amendment.
- python3 with numpy/scipy/matplotlib/sympy/mpmath (incl. mpmath.iv) and gmpy2 are installed.
- Your final structured output must faithfully report gate verdicts as graded — a FAIL is reported as FAIL.`

const SCHEMA = {
  type: 'object',
  required: ['workstream', 'gates', 'key_numbers', 'files', 'honest_notes', 'status'],
  properties: {
    workstream: { type: 'string' },
    gates: { type: 'array', items: { type: 'object', required: ['gate', 'requirement', 'measured', 'verdict'], properties: { gate: { type: 'string' }, requirement: { type: 'string' }, measured: { type: 'string' }, verdict: { type: 'string', enum: ['PASS', 'FAIL', 'PARTIAL'] } } } },
    key_numbers: { type: 'object' },
    files: { type: 'array', items: { type: 'string' } },
    honest_notes: { type: 'string' },
    status: { type: 'string', enum: ['COMPLETE', 'COMPLETE-WITH-AMENDMENTS', 'INCOMPLETE'] },
    finding_summary: { type: 'string' },
  },
}

const S1 = `${DISC}

WORKSTREAM S1 — P-acoustic adoption pilot: pricing cost (1) of the F-R14 sound cell. Directory: ${T8}/S1/.

CONTEXT (read these): ${ROOT}/tier7-program/O3/RESULTS.md (the decision matrix; the sound cell = P-acoustic density weight w = -1/4 + supertrace; SS4 lists the three printed costs — you are pricing cost (1): "reopens Sec. III's flat-measure [DF] derivations (Fisher/Madelung, Wallstrom, Born rates) on curved backgrounds"). ${ROOT}/corpus3/01-GUM-Omega-Paper-v4.3-ext.md Sec. III (the flat-measure Fisher/Madelung route to the quantum potential — find the actual printed derivation and its conventions; grep for Fisher, Madelung, quantum potential). ${ROOT}/tier7-program/O3/o3_matrix.py (conventions).

TASK: symbolic, sympy, nothing fitted.
(a) Reproduce the flat-measure identity exactly: variation of the Fisher information functional I_F[rho] = integral |grad rho|^2 / rho yields the Bohm quantum potential Q = -(hbar^2/2m) lap(sqrt rho)/sqrt rho with the standard coefficient linkage (state the campaign-units convention you adopt and cite where the paper prints it).
(b) Repeat on a curved 3-metric g_3 = (1+h) delta with quantization-measure weight rho_0, J proportional to (det g_3)^w, general symbolic w, to LINEAR order in h. Treat two cases separately: spatially constant h, and h(x) with nonzero gradient. Compute the corrected Euler-Lagrange/quantum-potential expression and isolate the h-dependent correction terms (which break the flat [DF] chain).
(c) Evaluate the correction at w = 0, w = -1/4 (P-acoustic), w = -1/2; determine whether ANY w cancels the corrections (solve symbolically); print the correction's structure and coefficients at w = -1/4.

GATES: S1-G1 flat identity exact; S1-G2 general-w linear-order expressions printed symbolically; S1-G3 the cost priced (cancellation w found exactly, or correction structure+coefficients at w=-1/4 printed, with honest statement of the effect on the [DF] label); S1-G4 decision-support-only framing, no adoption recommended, F-R14 regrade language unchanged. Deliver s1_price.py (deterministic, self-checking), s1_results.json, RESULTS.md (figure optional for symbolic work).`

const S2 = `${DISC}

WORKSTREAM S2 — Chirality-interpolation family: is the saturated closure constant distinctive? Directory: ${T8}/S2/.

CONTEXT (read these): ${ROOT}/tier2-closure/gstar_RESULTS.md IN FULL (the t->0 saturated objective G = integral [pi^3 b^2 + W(q)] dV, W = (1/8pi)[(1-cos psi)^2 + sin^2 psi cos^2 Theta]; sharp bound G* = 16 sqrt2 / 9 via pointwise AM-GM + degree pullback; closed-form direction-locked minimiser via g(x) = arcsin x - x sqrt(1-x^2); sector decomposition E6/E0/I via pushforward dV = dOmega/(2 sqrt(pi W)); identities E6 = E0 - I/16pi at d=1 and E_stat_tot = (3/2) G*). ${ROOT}/tier2-closure/gstar_solve.py (REUSE its conventions and quadrature machinery wherever possible — parametrization-reuse mandate; nothing refitted).

TASK: embed W in W_lambda = (1/8pi)[(1-cos psi)^2 + lambda sin^2 psi cos^2 Theta], lambda in [0,1] (lambda=1 GUM chiral, lambda=0 achiral).
(a) Sharp bound B(lambda) = (1/sqrt2) integral integral sin^2 psi sin Theta sqrt[(1-cos psi)^2 + lambda sin^2 psi cos^2 Theta] dpsi dTheta by high-order quadrature (target 1e-12); verify closed forms at both endpoints (B(1) = 16 sqrt2/9, B(0) = pi/sqrt2 — derive B(0) yourself from the integral to confirm).
(b) Attainment at >= 9 interior lambda values: the same per-theta BPS reduction (pi^3 b^2 = W_lambda pointwise, direction-locked Theta = theta) — determine whether the closed-form construction generalizes (the per-theta ODE integrates for general lambda?); measure the BPS saturation residual; if attainment fails off the endpoints, report that honestly (B(lambda) then only a bound).
(c) Sector decomposition on the lambda-minimiser via the pushforward: E6(lambda), E0(lambda), I(lambda); test E6 = E0 - I/16pi and E_stat = (3/2) G* at every sampled lambda, residuals printed; verdict per identity: GENERIC (holds for all lambda) / DISTINCTIVE (lambda=1 only) / OTHER.
(d) Anchor-hit scan: does c(lambda) = (4/pi) G*(lambda) cross any printed corpus constant (3.201124679..., 2.514753626*(4/pi)=3.201883... careful, list the printed constants you scan against from gstar_RESULTS.md) at lambda != 1? Report crossings with lambda values.

GATES: S2-G1 endpoints <= 1e-10; S2-G2 attainment residual <= 1e-8 at every sampled lambda OR honest attainment-fails statement; S2-G3 identity persistence table with GENERIC/DISTINCTIVE/OTHER verdicts; S2-G4 honest bottom line (distinctiveness is within-model structure, confers no physical support). Deliver s2_chiral.py, s2_results.json, s2_fig.png (G*(lambda), sector curves, identity residuals), RESULTS.md.`

const S3 = `${DISC}

WORKSTREAM S3 — O-1 pilot: validated forward flow-map covering (Q1 technology extension). Directory: ${T8}/S3/.

CONTEXT (read these): ${ROOT}/tier7-program/Q2/ANALYTIC_RECON.md SS5 (obligation O-1 — the dominant obligation: early-time flow control on [0,T0] x supp rho_0; "a validated-enclosure covering of the flow map... certifying where the not-yet-crossed set IS at T0"). ${ROOT}/tier7-program/Q1/RESULTS.md and ${ROOT}/tier7-program/Q1/q1_enclose.py IN FULL (the mpmath.iv Lohner QR engine, Picard boxes, certified band-limited field — REUSE this engine; adapt from single backward paths to forward box integration). ${ROOT}/tier6-foundations/ L5/L5prime RESULTS for the physical setup (arrival-time testbed, detector at x=d, tau_max = 5.130950 at n_y=1).

TASK (pilot; NOT an O-1 discharge — the deliverable is the feasibility datum):
(a) Engine validation: forward-integrate the time-reverse of one Q1-certified backward path; the endpoint enclosure must CONTAIN the Q1 anchor point (containment, not distance).
(b) Covering: n_y = 1; initial boxes at t0 = 1.0 covering the pre-registered subdomain (the Q2 T_ball band mapped to the testbed's (x, z) coordinates — read Q2/Q1 to fix the mapping honestly; if the roadmap's (w, z) rectangle w in [-1.2,-0.6], z in [0.05,0.95] maps awkwardly, print the coordinate choice as a spec clarification): M = 24 boxes (6 x 4), integrate each forward to T0 = 5.14. Per-box nohup + JSON checkpoint; benchmark ONE box first and scale (fewer boxes = printed amendment) to fit the wall-clock budget.
(c) Classification per box: (i) CROSSED (enclosure wholly past x = d before T0), (ii) LOCALIZED (wholly x < d at T0, enclosure printed), (iii) BLOWUP (width > 1 before T0) -> subdivide once 2x2, re-run children, then stop.
(d) Feasibility datum: fraction (i)/(ii)/FAIL, median validated-step count and wall-cost per box, width-growth factors, and the honest extrapolated cost of a full O-1 covering.

GATES: S3-G1 containment PASS; S3-G2 covering executed with zero unaccounted subdomain; S3-G3 feasibility datum printed incl. full-O-1 extrapolation; S3-G4 scope honesty (pilot != discharge; kill-bin theorem remains open; register unchanged). Deliver s3_cover.py, per-box checkpoints, s3_results.json, s3_fig.png (subdomain map colored by classification), RESULTS.md.`

const S4 = `${DISC}

WORKSTREAM S4 — Population-scale certification screener (P2 extension). Directory: ${T8}/S4/.

CONTEXT (read these): ${ROOT}/tier7-program/P2/RESULTS.md IN FULL and ${ROOT}/tier7-program/P2/p2_certify.py (the 432-configuration exact-field certification: parametrization grid, margin definition, the mpmath precision ladder, the G3 endpoint FAIL carried honestly). Your population MUST be drawn from the SAME parametrization, refined (parametrization-reuse mandate).

TASK:
(a) Screener validation: a float64 re-implementation of the P2 margin pipeline; on 24 randomly chosen (seeded, printed seed) P2 configs, reproduce the P2 margins to <= 1e-6 relative (if the P2 pipeline's own precision floor makes 1e-6 unattainable in float64, measure and print the actual agreement floor and its cause as a finding, gate graded against the honest achievable bar with the deviation printed as a spec amendment).
(b) Population screen: refine the P2 grid toward N = 10,000 configs; benchmark per-config float cost FIRST and scale N to fit the budget (N >= 2,000 minimum; any scale-down from 10,000 printed as amendment). nohup + chunked checkpoints (resume-safe). Deliver the margin distribution, the minimum, and the location of the worst configs in parameter space.
(c) Worst-K certification: K = 24 lowest-margin configs re-run through the P2 high-precision ladder (200-bit mpmath, nohup + per-config checkpoints); verdict per config STABLE/UNSTABLE with certified digits; any UNSTABLE is a FINDING (report it, do not soften).
(d) Honest statement: screened coverage != formally-verified kernel; that open remains open.

GATES: S4-G1 screener agreement at the honest bar; S4-G2 screen complete at final N with distribution + minimum; S4-G3 worst-K certified with per-config verdicts; S4-G4 honest scope statement. Deliver s4_screen.py, s4_results.json, s4_fig.png (margin distribution + worst-config map), RESULTS.md.`

phase('Execute')
const results = await parallel([
  () => agent(S1, { label: 'S1:p-acoustic-cost', phase: 'Execute', schema: SCHEMA }),
  () => agent(S2, { label: 'S2:chirality-family', phase: 'Execute', schema: SCHEMA }),
  () => agent(S3, { label: 'S3:o1-pilot-covering', phase: 'Execute', schema: SCHEMA }),
  () => agent(S4, { label: 'S4:population-screen', phase: 'Execute', schema: SCHEMA }),
])

return { S1: results[0], S2: results[1], S3: results[2], S4: results[3] }