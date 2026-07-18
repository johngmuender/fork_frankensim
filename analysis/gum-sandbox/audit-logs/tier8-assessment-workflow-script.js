export const meta = {
  name: 'tier8-S-assessment',
  description: 'Program-wide GUM assessment: sector surveys, pros/cons briefs, adversarial audit',
  phases: [
    { title: 'Survey', detail: 'six parallel sector surveyors' },
    { title: 'Briefs', detail: 'pros-advocate and cons-critic' },
    { title: 'Audit', detail: 'adversarial evidence audit of both briefs' },
  ],
}

const ROOT = '/home/user/fork_frankensim/analysis/gum-sandbox'

const COMMON = `You are a survey agent in a multi-agent assessment of the GUM (Geometrische Umdeutung Mechanik) speculative-physics replication campaign in ${ROOT}. GUM is a chiral micropolar (Cosserat) substrate theory; the campaign has run Tiers 0-7 plus Phases N/O/P/Q/R of numerical replication, adversarial testing, and document folds. Register discipline: everything is "within-model; nothing bears on nature" — the campaign tests internal consistency and replication of the corpus's claims, not physical truth. Ledger classes: [DF]=derived-fixed, [DW]=derived-weak, [CAL]=calibrated, [IM]=imported, [CJ]=conjecture.

Your job: read the files listed below (use Read/Grep; be selective — read RESULTS.md / adjudication / analysis files fully, skim scripts) and return a STRUCTURED assessment of your sector. For every claim you make, cite a specific file path (and section/line if possible) as evidence. Distinguish: (a) genuine strengths (things the program does well — methodology, honesty, reproducibility, derivational economy), (b) genuine weaknesses (circularity, calibration disguised as derivation, unfalsifiable elements, gaps, failed gates, scope limits), (c) open questions. Be adversarially honest in BOTH directions: do not inflate strengths, do not manufacture weaknesses. FAILed gates that were honestly reported are a methodological strength AND possibly a content weakness — record both sides where true.`

const SCHEMA = {
  type: 'object',
  required: ['sector', 'strengths', 'weaknesses', 'opens'],
  properties: {
    sector: { type: 'string' },
    strengths: { type: 'array', items: { type: 'object', required: ['claim', 'evidence'], properties: { claim: { type: 'string' }, evidence: { type: 'string' }, weight: { type: 'string', enum: ['major', 'minor'] } } } },
    weaknesses: { type: 'array', items: { type: 'object', required: ['claim', 'evidence'], properties: { claim: { type: 'string' }, evidence: { type: 'string' }, weight: { type: 'string', enum: ['major', 'minor'] } } } },
    opens: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
}

const SECTORS = [
  { key: 'theory-core', files: `${ROOT}/corpus3/01-GUM-Omega-Paper-v4.3-ext.md (read fully - ~970 lines), ${ROOT}/corpus3/REVISION_CHARTER_v2.md, ${ROOT}/theory-audit/ (all files)`, focus: 'The theory itself as presented in the v4.3 paper: derivational structure, hbar-closure, invariants, ledger-class honesty, what is derived vs calibrated vs imported, internal consistency, the T6/T7 fold blocks.' },
  { key: 'early-tiers', files: `${ROOT}/tier0-gauntlet/, ${ROOT}/tier1-spectrum/, ${ROOT}/tier2-closure/, ${ROOT}/tier3-born/ (RESULTS.md, ADJUDICATION files, ANALYSIS files in each)`, focus: 'Tiers 0-3: the gauntlet, spectrum replication, closure saturation, Born-rule emergence. What replicated, what failed, quality of gates.' },
  { key: 'field-family', files: `${ROOT}/tier4-field/, ${ROOT}/tier5-family/ (RESULTS.md, adjudication/analysis files)`, focus: 'Tiers 4-5: field-theoretic tests and family/parameter structure. Includes T4-W3, T4-W5 lineage. What replicated, what failed.' },
  { key: 'foundations', files: `${ROOT}/tier6-foundations/ (FOUNDATIONS_ANALYSIS.md, TIER6_ADJUDICATION.md, ROADMAP_v8_FOUNDATIONS.md, L1/L2/L3/L5/L5prime/L7a/L7b RESULTS.md)`, focus: 'Tier 6: GUM vs quantum-foundations literature (Bell, dBB, Nelson, GRWm). The discriminator ledger, POVM-exclusion testbed, arrival-time work. Strength of the foundations positioning.' },
  { key: 'program-phases', files: `${ROOT}/tier7-program/ (PROGRAM_ANALYSIS_T7.md, TIER7_ADJUDICATION.md, ROADMAP_v9_PROGRAM.md, M1-M4/N2/N3/O1-O3/P1-P2/Q1-Q2 RESULTS.md)`, focus: 'Tier 7 + Phases N-Q: cosmology sector, servo, F-R9 closure, relaxation, second-moment rate honesty, backward-reachability proof, g_lim archaeology, F-R14 matrix, exact-field certification, validated enclosures, Lemma A. The certification chain quality.' },
  { key: 'campaign-meta', files: `${ROOT}/REPLICATION_CAMPAIGN_STATUS.md, ${ROOT}/SESSION_HANDOFF.md, ${ROOT}/DISCHARGE_PACKAGE/ , ${ROOT}/corpus3/FOLD_RECORD.md, ${ROOT}/corpus3/09-External-Replication-Record-v2.0.md, ${ROOT}/audit-logs/MANIFEST_tier6.json (skim)`, focus: 'The campaign as a scientific process: pre-registration discipline, adversarial honesty, reproducibility infrastructure, audit trail, fold protocol, restoration docs. Meta-level strengths/weaknesses of the replication methodology itself.' },
]

phase('Survey')
const surveys = await parallel(SECTORS.map(s => () =>
  agent(`${COMMON}\n\nYOUR SECTOR: ${s.key}\nFILES: ${s.files}\nFOCUS: ${s.focus}\n\nReturn 5-12 strengths and 5-12 weaknesses (each with file evidence and weight major/minor), plus open questions.`,
    { label: `survey:${s.key}`, phase: 'Survey', schema: SCHEMA })
))

const okSurveys = surveys.filter(Boolean)
log(`surveys complete: ${okSurveys.length}/6 sectors`)

const surveyDigest = JSON.stringify(okSurveys)

phase('Briefs')
const BRIEF_SCHEMA = {
  type: 'object',
  required: ['thesis', 'points'],
  properties: {
    thesis: { type: 'string' },
    points: { type: 'array', items: { type: 'object', required: ['title', 'argument', 'evidence'], properties: { title: { type: 'string' }, argument: { type: 'string' }, evidence: { type: 'string' }, weight: { type: 'string', enum: ['major', 'minor'] } } } },
  },
}
const [prosBrief, consBrief] = await parallel([
  () => agent(`You are the PROS ADVOCATE in an assessment of the GUM program at ${ROOT}. Six sector surveys (JSON below) give strengths/weaknesses with file evidence. Build the STRONGEST honest case FOR the program: 8-14 points, each grounded in survey evidence (verify key claims yourself with Read/Grep against the cited files where feasible). Cover both content (what the theory achieves within-model) and process (methodology). Do not overstate: a point that inflates evidence will be struck by an adversarial auditor. Register: within-model; nothing bears on nature.\n\nSURVEYS:\n${surveyDigest}`, { label: 'brief:pros', phase: 'Briefs', schema: BRIEF_SCHEMA }),
  () => agent(`You are the CONS CRITIC in an assessment of the GUM program at ${ROOT}. Six sector surveys (JSON below) give strengths/weaknesses with file evidence. Build the STRONGEST honest case AGAINST the program: 8-14 points, each grounded in survey evidence (verify key claims yourself with Read/Grep against the cited files where feasible). Cover content (circularity, calibration, unfalsifiability, failed gates, scope) and process (any methodological weaknesses). Do not manufacture: a point that misreads evidence will be struck by an adversarial auditor. Register: within-model; nothing bears on nature.\n\nSURVEYS:\n${surveyDigest}`, { label: 'brief:cons', phase: 'Briefs', schema: BRIEF_SCHEMA }),
])

phase('Audit')
const AUDIT_SCHEMA = {
  type: 'object',
  required: ['verdicts'],
  properties: {
    verdicts: { type: 'array', items: { type: 'object', required: ['point_title', 'verdict', 'reason'], properties: { point_title: { type: 'string' }, verdict: { type: 'string', enum: ['SUPPORTED', 'OVERSTATED', 'UNSUPPORTED'] }, reason: { type: 'string' }, correction: { type: 'string' } } } },
  },
}
const audits = await parallel([
  () => agent(`You are an adversarial AUDITOR. Below is a PROS brief about the GUM program at ${ROOT}. For EACH point, check its evidence against the actual files (Read/Grep). Verdict SUPPORTED only if the file evidence genuinely backs the claim at the stated strength; OVERSTATED if the kernel is real but inflated (give correction); UNSUPPORTED if the evidence does not back it. Default skeptical.\n\nBRIEF:\n${JSON.stringify(prosBrief)}`, { label: 'audit:pros', phase: 'Audit', schema: AUDIT_SCHEMA }),
  () => agent(`You are an adversarial AUDITOR. Below is a CONS brief about the GUM program at ${ROOT}. For EACH point, check its evidence against the actual files (Read/Grep). Verdict SUPPORTED only if the file evidence genuinely backs the criticism at the stated strength; OVERSTATED if the kernel is real but inflated (give correction); UNSUPPORTED if the evidence does not back it (e.g. the criticism ignores an honest disclosure or a ledger-class label that already concedes the point). Default skeptical.\n\nBRIEF:\n${JSON.stringify(consBrief)}`, { label: 'audit:cons', phase: 'Audit', schema: AUDIT_SCHEMA }),
])

return { surveys: okSurveys, prosBrief, consBrief, prosAudit: audits[0], consAudit: audits[1] }