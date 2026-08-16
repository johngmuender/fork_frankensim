export const meta = {
  name: 'tier10-U-execute',
  description: 'Execute ROADMAP v12 junction workstreams U1-U3',
  phases: [
    { title: 'Execute', detail: 'three parallel gated workstreams' },
  ],
}

const ROOT = '/home/user/fork_frankensim/analysis/gum-sandbox'
const T10 = `${ROOT}/tier10-junction`

const DISC = `EXECUTION DISCIPLINE (binding). REQUIRED FIRST READS, in order: (1) ${T10}/ROADMAP_v12_JUNCTION.md — your workstream section IS the frozen spec; gates are PRE-REGISTERED and may not move; (2) ${T10}/JUNCTION_CONTEXT_ANALYSIS.md; (3) ${T10}/JUNCTION_PAPER_DIGEST.md. Also: ${T10}/u_context_agents.json (raw context outputs with file:line evidence — mine it rather than re-searching).
- Within-model; nothing bears on nature. STAR/lattice/Regge values are [IM] anchors. Honest FAILs are findings. Grade language never rises. All corpus-side changes are OFFERS.
- SEAL q-THETA: any hadron-mass-confronting number reportable ONLY at V.F grade ('dimensionally secure / structurally plausible / quantitatively unclaimed'). NO FITTING to STAR observables (1.84, 0.64, 1.04) anywhere — they may appear only as [IM] anchors in confrontation tables, never as targets. The junction interpretation of STAR data is itself contested (CGC-saturation rival, neutron-skin, strangeness confounds) — carry rivals at equal weight.
- FORBIDDEN SENTENCES: 'GUM predicted the junction result', 'STAR confirms/refutes GUM', any grade rise, any stake issuance. Numerology hazard (dimensionless c-frak = 2.37 vs GeV digit coincidences) pre-empted in every deliverable.
- Work in your workstream directory (create it). Deliver: script(s) where computational, machine-readable results JSON, figure where meaningful (matplotlib Agg), RESULTS.md in campaign house style (header phase/date/code/outputs/sources; epistemic frame; numbered sections; gates table requirement/measured/verdict; key-numbers table; honest caveats).
- python3 has numpy/scipy/matplotlib/sympy/mpmath/gmpy2. Long runs (>10 min): nohup + checkpoints (none expected).
- Spec defects: closest feasible variant + PRINT the coordinator-owned amendment. Structured output reports gates as graded — FAIL as FAIL.`

const SCHEMA = {
  type: 'object',
  required: ['workstream', 'gates', 'key_numbers', 'files', 'honest_notes', 'status', 'finding_summary'],
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

const U1 = `${DISC}

WORKSTREAM U1 — Junction-antijunction (J-Jbar) spectrum from printed parameters. Directory: ${T10}/U1/.
Your spec: ROADMAP_v12_JUNCTION.md section U1 (governing). Key elements:
- REUSE ${ROOT}/tier9-glueball/T1/t1_spectrum.py (read it fully): the variable-mass radial Schrodinger solver, frozen cutoff f* = 1.9724, Zhu-Kroemer/BDD orderings. Extend to the dumbbell: separation coordinate L between two junctions joined by three parallel tubes; leading energy 3*sigma*L (+ 2*m_J); configuration inertia mu(L) from tube mass 3*sigma*L + 2*m_J (imitate T1's mu(rho) = 2*pi*sigma*rho pattern); transverse-phonon zero-point content: three tubes x 2 polarizations with the T1 census conventions and the frozen cutoff as junction-end smoothing (transfer choices = printed amendments). Y-law end geometry note: at leading order the three tubes are parallel (junction-antijunction axis); state this geometric idealization as an amendment.
- VALIDATION FIRST (U1-G1): (a) two-body linear-potential spectrum vs exact Airy zeros to <= 1e-6 (set the reduced-mass/tension parameters to a textbook case); (b) reproduce a T1 limit (e.g. your solver with T1's loop potential must reproduce T1's 0++ epsilon = 3.583 at sigma = 0.19 to <= 1e-3 rel).
- Spectrum (U1-G2): m_J/sqrt(sigma) in {0, 0.1355 [IM 2+1D lattice anchor - dimensional-transfer caveat printed], 0.39 [= f_q central 0.17/0.436]}; lowest 3+ states each; scheme envelope as uncertainty; zero fitted parameters.
- Confrontation (U1-G3, V.F): m(J-Jbar ground)/sqrt(sigma) vs Tier-9 anchors (X 5.41-5.45; lattice 0-+ 5.5-5.9; 0++ 3.5) and vs T1's loop states (IP 0++ 3.583, 0-+ 8.275; NG 4.80): does the junction sector interleave/duplicate/extend the loop sector? Tetraquark-kill note (n=3 only). Either answer is the finding.
Deliver u1_dumbbell.py, u1_results.json, u1_fig.png, RESULTS.md. File as F-T10-U1.`

const U2 = `${DISC}

WORKSTREAM U2 — Two-carrier transport sign test. Directory: ${T10}/U2/.
Your spec: ROADMAP_v12_JUNCTION.md section U2 (governing). Key elements:
- The toy (pre-registered by the roadmap, sign/monotonicity ONLY): channel V (valence-led) moves B by transporting three knots of mass scale M_knot across a rapidity interval; channel J (junction-led) migrates the junction+tube structure (scales: sigma = 0.19 GeV^2, m_J from U1's parameter set {0, 0.1355, 0.39}*sqrt(sigma)) and pays 2x pair-minting cost at the printed ~100-MeV-class snap scale (Primer ex. 15.1 chain: ~1e5 N x 1e-15 m — cite it). Q rides only on knots (V.D winding; K-1). Build the simplest cost/suppression model that lets you compute the SIGN of the B-vs-Q transport asymmetry (e.g. Boltzmann/exponential suppression in transported inertia across a rapidity gap — state your suppression ansatz as the [CJ-new] construction it is, print it, and show the sign is ansatz-robust for any monotone-decreasing cost function).
- Scan the FULL printed windows: M_knot in [m_Sk = 1.7 GeV, m-tilde_t-class 172.5 GeV]; f_q in [0.14, 0.20]; m_J over the U1 set; snap scale 100-140 MeV. NO point tuning; NO use of 1.84/0.64/1.04 as targets (mechanical self-scan required in G4).
- Validation limits (U2-G1): minting disabled -> junction channel closed -> asymmetry exactly 0 (valence-only); tube+junction+minting cost -> 0 -> maximal junction dominance. Both exact in the code.
- Output (U2-G2/G3): the asymmetry sign over the whole window (robustness map figure); the reconciliation statement as computed: is K-1 (knot-carried B) compatible with junction-led transport — i.e. does B move without any original knot moving, with net core count conserved? State the carrier/enforcer/transporter trichotomy verdict. Magnitudes UNCLAIMABLE (the corpus prints no transport law — this is a [CJ-new] structural-genericity test, not a prediction).
Deliver u2_transport.py, u2_results.json, u2_fig.png, RESULTS.md. File as F-T10-U2.`

const U3 = `${DISC}

WORKSTREAM U3 — B-dictionary audit + scorecard + priced alternative + annotations. Directory: ${T10}/U3/.
Your spec: ROADMAP_v12_JUNCTION.md section U3 (governing). Key elements:
- (a) Dictionary memo (U3-G1): Theorem K-1's full chain from the archaeology (census crisis C-K1 -> residue lock rho = w_em mod 1 -> K-1/K-2/K-3; K-frak = 3B + L, F-K1-2), every claim file:line-cited (mine u_context_agents.json junctionArch section, then VERIFY the load-bearing lines yourself in corpus/WS-K-*.md); the worksheet-stratum caveat explicit (zero operative-edition hits for 'baryon number' — the unaffixed delta-pack stubs = FOLD DEBT, file it); the F-K0-1 junction-order seam (paper-conditional vs Primer-asserted n=3) carried with all edition sites.
- (b) Formal conservation audit (U3-G2): for each printed operation — (i) taping-rule termination, (ii) web reconnection (net-conserving, census column), (iii) snap-minting (pair creation) — plus (iv) the junction-migration scenario (U2's mode): does B = (1/3)*Sigma(signed residues) change? Work each case explicitly from the printed rules; any failure is a finding.
- (c) Scorecard (U3-G3): rows = STAR's three measurements + the compact-B imaging result (arXiv:2603.03730) + the junction-exotics class; columns = corpus-printed content (file:line) / U1-U2 computed content (POLL for ${T10}/U1/RESULTS.md and ${T10}/U2/RESULTS.md with a watcher loop, up to ~40 min; on timeout fill PENDING + print amendment) / verdict in {SUPPORTED-STRUCTURALLY, V.F-CLASS, SILENT-NEEDS-NEW-WORK, TENSION} / MANDATORY RIVALS column (CGC saturation PRC 111 024912; neutron skin; strangeness AMPT) so no row overstates the junction interpretation. Verdict language within the Tier-9 stakes constraints (no stake, no clock, no grade motion — reuse tier9-glueball/t_context_agents.json stakes section).
- (d) Priced alternative (in G3's RESULTS section): 'K-1-as-printed only' (carrier answered, transporter silent — the STAR-shaped question unanswerable in-model) vs 'K-1 + junction-led transport completion' (costs: [CJ-new] transport layer + junction energy + unprinted dynamics; buys: the trichotomy dissolving the apparent tension). Symmetric, neither recommended.
- (e) Annotations (U3-G4): T10 offer drafts modeled on corpus3/01's T6 boxes: (i) VII.J carrier/enforcer/transporter note; (ii) K-stub affixation reminder (fold debt); (iii) F-K0-1 seam cross-reference; (iv) junction-network state class -> third entry of Tier-9's two-sector dictionary (T4/t4_annotations.md - read it); (v) optional WS1-H3-discharge suggestion citing Komargodski-Zhong PRD 110 056018 + arXiv:2508.00608 lattice junction mass as the adjudicator class. All strike-able, cited, zero grade motion; mechanical forbidden-sentence scan.
Deliver RESULTS.md, u3_scorecard.json, u3_annotations.md. File as F-T10-U3.`

phase('Execute')
const results = await parallel([
  () => agent(U1, { label: 'U1:jjbar-spectrum', phase: 'Execute', schema: SCHEMA }),
  () => agent(U2, { label: 'U2:transport-sign', phase: 'Execute', schema: SCHEMA }),
  () => agent(U3, { label: 'U3:dictionary-scorecard', phase: 'Execute', schema: SCHEMA }),
])

return { U1: results[0], U2: results[1], U3: results[2] }