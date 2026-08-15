export const meta = {
  name: 'tier9-T-execute',
  description: 'Execute ROADMAP v11 glueball workstreams T1-T4',
  phases: [
    { title: 'Execute', detail: 'four parallel gated workstreams' },
  ],
}

const ROOT = '/home/user/fork_frankensim/analysis/gum-sandbox'
const T9 = `${ROOT}/tier9-glueball`

const DISC = `EXECUTION DISCIPLINE (binding). REQUIRED FIRST READS, in order: (1) ${T9}/ROADMAP_v11_GLUEBALL.md — your workstream section IS the frozen spec; gates are PRE-REGISTERED and may not move; (2) ${T9}/GLUEBALL_CONTEXT_ANALYSIS.md; (3) ${T9}/GLUEBALL_PAPER_DIGEST.md. Also available: ${T9}/t_context_agents.json (raw context-agent outputs with file:line evidence — mine it rather than re-searching).
- Within-model; nothing bears on nature. BESIII/lattice values are [IM] imported anchors. Honest FAILs are findings. Grade language never rises. All corpus-side changes are OFFERS.
- SEAL q-THETA discipline: any within-model mass number confronting a hadron mass is reportable ONLY at V.F grade ('dimensionally secure / structurally plausible / quantitatively unclaimed'). The pipeline is frozen by the roadmap — no tuning toward 2359/2376 MeV. The numerology hazard (dimensionless c-frak = 2.37 vs m_X = 2.37 GeV digit coincidence) must be pre-empted explicitly in your RESULTS.md and never cited as structure.
- FORBIDDEN SENTENCES everywhere: 'GUM predicted X(2370)', 'X(2370) confirms/refutes GUM', any grade rise, any stake issuance.
- Work in your workstream directory (create it). Deliver: script(s) where computational, machine-readable results JSON, figure where meaningful (matplotlib Agg), and RESULTS.md in campaign house style (header: phase/date/code/outputs/sources; epistemic frame; numbered sections; gates table requirement/measured/verdict; key-numbers table; honest caveats).
- Long runs (>10 min): nohup + per-item JSON checkpoints (none expected here). python3 has numpy/scipy/matplotlib/sympy/mpmath/gmpy2.
- Spec defects: implement closest feasible variant + PRINT the deviation as a coordinator-owned amendment.
- Structured output must report gate verdicts as graded — FAIL reported as FAIL.`

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

const T1 = `${DISC}

WORKSTREAM T1 — Closed-tube spectrum from printed corpus parameters. Directory: ${T9}/T1/.
Your spec: ROADMAP_v11_GLUEBALL.md section T1. Key elements (the roadmap text governs):
- VALIDATION-FIRST: implement the Isgur-Paton adiabatic loop quantization (Isgur & Paton, PRD 31, 2910 (1985): transverse phonons of frequency m/rho on a circular flux loop of radius rho, each phonon carrying angular momentum +-m about the loop axis; adiabatic effective potential V_eff(rho) = 2 pi rho sigma + phonon zero-point/occupation energies with the IP prescription incl. their short-distance cutoff/fudge treatment; radial Schroedinger equation in rho). It MUST reproduce IP's published lowest 0++ ~ 1.52 GeV at their sigma to <= 3% BEFORE any GUM number enters (T1-G1). Document the IP parameter conventions you adopt with citations (their sigma ~ 0.18 GeV^2 via b = 0.18 GeV^2? — pin exactly what reproduces 1.52 and print the provenance; if the exact IP prescription cannot be pinned from accessible literature, implement the closest published variant and print the deviation as a coordinator-owned amendment with your achieved validation number).
- SECOND ROUTE: exact free Nambu-Goto closed-string levels; for glueball reading use the standard closed-string glueball treatments (zero winding circular/folded closed strings; E^2 = 8 pi sigma [(N_L+N_R)/2 - (D-2)/24] with N_L = N_R level matching at zero momentum, D = 4; also the Sonnenschein-Weissman closed-slope alpha'/2 Regge cross-check). Machine-precision identity checks on the formula.
- Evaluate BOTH routes at the corpus tension sigma = 0.19 GeV^2 [IM]; lowest 3+ states per route, J^PC as the models assign (IP phonon content m=2 doublet -> 0++/2++ etc.; carry the 0- worldsheet-axion caveat from the context analysis); uncertainty bands from model systematics ONLY (IP cutoff freedom; NG short-string validity: expansion parameter 1/(sigma l^2) printed per state), zero fitted parameters.
- Confrontation table (T1-G3, V.F-graded): m/sqrt(sigma) per state vs [IM] anchors: X(2370) m/sqrt(sigma) = 5.41-5.45; lattice 0-+/0++ = 1.50; lattice m(0++)/sqrt(sigma) ~ 3.5, m(0-+)/sqrt(sigma) ~ 5.5-5.9. Answer the pre-registered question honestly: does the corpus-parametrized closed tube put its lightest pseudoscalar in the X(2370)/lattice-0-+ class or not? EITHER answer is the finding.
Deliver t1_spectrum.py (deterministic, self-checking), t1_results.json, t1_fig.png (spectrum vs anchors), RESULTS.md. File as F-T9-T1.`

const T2 = `${DISC}

WORKSTREAM T2 — Worldsheet-mode census on the GUM tube (the chirality hook). Directory: ${T9}/T2/.
Your spec: ROADMAP_v11_GLUEBALL.md section T2. Key elements:
- Background: REUSE the campaign's straight-static-tube machinery — read ${ROOT}/theory-audit/h21_RESULTS.md (and h21 script if present), ${ROOT}/tier7-program/M3/RESULTS.md, ${ROOT}/tier7-program/N2/RESULTS.md for the tube parametrization (mandate: cite what you reuse). Also read the printed substrate energy structure: corpus3/01-GUM-Omega-Paper-v4.3-ext.md Sec. II (micropolar + chiral chi_1/chi_2/chi_3 classes — grep 'chi_1\\|chi_2\\|chi_3\\|wryness\\|chiral gate'), Sec. VII.J, and the h25 branch-spectrum table (${ROOT}/theory-audit/h25_RESULTS.md).
- Mode enumeration (T2-G1): table of candidate worldsheet degrees of freedom on a straight T2 mismatch tube: 2 transverse translation Goldstones; locking-stratum orientation modes on the core (Q-2/Q-5 orientation field, stiffness f_q); phason (T3) admixtures via printed cross-stratum couplings; anything else the printed field content supports. Each row: printed evidence (file:line) or [CJ-new] flag.
- Symbolic reduction (T2-G2, sympy or careful hand-derivation shown step-by-step): quadratic worldsheet action for each candidate mode from the printed energy terms. VALIDATION: the 2 transverse modes MUST come out massless (Goldstone check) — if they do not, your reduction is wrong. For internal modes: mass scale symbolic (M_gap, f_q — no invented numbers; M_gap is symbol-only per the context analysis).
- Parity classification (T2-G3): for each mode, worldsheet parity; the question: does ANY printed chi-class chiral coupling induce a parity-odd (pseudoscalar) worldsheet mode on the tube — GUM's candidate for the worldsheet axion that 4D closed-string glueball spectra require? Every unprinted assumption flagged [CJ-new] inline. Verdict: plain-bosonic-string vs string-with-internal-modes (+ pseudoscalar or not). If NO pseudoscalar emerges from printed structure, THAT is the finding (GUM tube inherits the known 3+1D IP/NG J^PC problems).
Deliver t2_census.py (symbolic derivations executable/self-checking where feasible), t2_results.json, RESULTS.md (figure optional). File as F-T9-T2.`

const T3 = `${DISC}

WORKSTREAM T3 — Seven-criteria scorecard + the closure obligation. Directory: ${T9}/T3/.
Your spec: ROADMAP_v11_GLUEBALL.md section T3. Key elements:
- The scorecard (T3-G1): seven rows = BESIII criteria (mass; spin-parity 0-+; production rate; eta_c decay-pattern similarity; flavor-singlet; narrow partial widths; gamma-omega/gamma-phi suppression). Columns: (a) what the corpus-as-printed says (file:line from t_context_agents.json + your own verification reads), (b) what T1/T2 computed (READ ${T9}/T1/RESULTS.md and ${T9}/T2/RESULTS.md — they run in parallel with you, so POLL for them: check every few minutes with a watcher loop until both exist or ~40 min elapse; if either is missing after that, fill those cells 'PENDING-T1/T2' and print the amendment), (c) verdict per row from the four pre-registered categories: SUPPORTED-STRUCTURALLY / V.F-CLASS (number, unclaimed) / SILENT-NEEDS-NEW-WORK / TENSION.
- The closure obligation (T3-G2): formalize as a named open problem in F-Q2-prime style. From print: Q-6-prime language lock 'only asymptotic states owe closure' (corpus3/01 line ~615); the closure machinery Sec. IV (spin-clock, w*j(1-j)=1, for knots); 'line-neutral composites close at the emergent level'. State the requirement; enumerate the three branches (emergent-level closure mechanism for a knotless state / principled excusal / cannot-close => contradiction with the WS-D relic census that presupposes the class); adjudicate what print answers (expected NOTHING — NEW-OPEN verdict); draft the proposed obligation text in offer language.
- Relic-census consistency note (T3-G3): WS-D D-0 'stratum composites (glueball-analogs) -> decay to hadrons/phasons [check]' vs the [IM] measured Gamma_X = 170 MeV hadronic width: qualitative census-grade consistency note, no overclaim either direction.
- T3-G4: stakes-machinery constraints verbatim (read the stakes agent's output in t_context_agents.json): no stake, no board clock, forbidden sentences absent, V.F on all numbers.
Deliver RESULTS.md (+ scorecard also as t3_scorecard.json). File as F-T9-T3.`

const T4 = `${DISC}

WORKSTREAM T4 — The two-sector dictionary (alternative-component analysis). Directory: ${T9}/T4/.
Your spec: ROADMAP_v11_GLUEBALL.md section T4. Key elements:
- Dictionary table (T4-G1): state classes (baryons; K-Kbar mesons; glueball-analogs; hybrids-analogs if any; the deuteron-analog B=2) vs owning sector (Q-5 emergent Skyrme orientation field vs T2 closed-string/tube sector) vs printed evidence (file:line) vs literature precedent (Skyrme model has no glueballs — cite standard references; flux-tube/Skyrme complementarity in real QCD).
- Double-counting audit: states both sectors could claim (scalar glueball vs Skyrmion breathing/monopole modes; 0-+ vs any pseudoscalar mode of the orientation field), with the standard literature resolutions.
- Priced alternative (T4-G2), symmetric honesty: 'Q-5 as printed (Skyrme-only)' — cost: the glueball-dominant X(2370)-class state is HOMELESS (the corpus's own relic census presupposes the class the identification cannot house); vs 'Q-5 + closed-string completion' — cost: new structure (loop quantization, closure obligation, worldsheet content per T2), itemized. Neither recommended — priced.
- Proposed T9 annotation set (T4-G3), charter-compliant offer language (strike-able, cites finding + artifact path, grades never rise): (i) the VII.J two-sector note; (ii) the closure obligation (READ ${T9}/T3/RESULTS.md — poll with a watcher for up to ~40 min; if missing, draft from the roadmap's T3 spec and print the amendment); (iii) the relic-census cross-reference (WS-D D-0 <-> VII.J); (iv) OPTIONALLY a T-N5-adjacent monitoring-row suggestion (offer only). Model the annotation format on existing T6/T7 boxes in corpus3/01 (grep 'Proposed annotation').
Deliver RESULTS.md (+ t4_annotations.md with the draft annotation texts). File as F-T9-T4.`

phase('Execute')
const results = await parallel([
  () => agent(T1, { label: 'T1:closed-tube-spectrum', phase: 'Execute', schema: SCHEMA }),
  () => agent(T2, { label: 'T2:worldsheet-census', phase: 'Execute', schema: SCHEMA }),
  () => agent(T3, { label: 'T3:scorecard-closure', phase: 'Execute', schema: SCHEMA }),
  () => agent(T4, { label: 'T4:two-sector-dictionary', phase: 'Execute', schema: SCHEMA }),
])

return { T1: results[0], T2: results[1], T3: results[2], T4: results[3] }