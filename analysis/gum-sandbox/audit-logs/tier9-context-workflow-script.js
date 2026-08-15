export const meta = {
  name: 'tier9-T-context',
  description: 'Glueball context: corpus color-sector archaeology + field literature sweep',
  phases: [
    { title: 'Context', detail: 'three corpus readers + two literature agents + stakes reader' },
  ],
}

const ROOT = '/home/user/fork_frankensim/analysis/gum-sandbox'

const COMMON = `You are a context agent in Tier 9 (Phase T) of the GUM replication campaign at ${ROOT}. The campaign is analyzing the GUM corpus (a speculative chiral-micropolar-substrate theory; within-model register: nothing bears on nature) in the context of BESIII's 2026 identification of X(2370) as dominantly the lightest 0-+ glueball (arXiv:2607.20366). REQUIRED FIRST READ: ${ROOT}/tier9-glueball/GLUEBALL_PAPER_DIGEST.md (the coordinator's extreme-detail digest of that paper, including the corpus's color-sector hooks in Sec. 6). Cite every corpus claim with file path + line/section. Distinguish [DF]/[DW]/[CAL]/[IM]/[CJ] ledger classes faithfully. Honest gaps are findings: if the corpus is silent on something, say SILENT, do not infer.`

const CORPUS_SCHEMA = {
  type: 'object',
  required: ['agent', 'findings', 'key_numbers', 'silences', 'notes'],
  properties: {
    agent: { type: 'string' },
    findings: { type: 'array', items: { type: 'object', required: ['topic', 'content', 'evidence', 'ledger_class'], properties: { topic: { type: 'string' }, content: { type: 'string' }, evidence: { type: 'string' }, ledger_class: { type: 'string' } } } },
    key_numbers: { type: 'object' },
    silences: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
}

const LIT_SCHEMA = {
  type: 'object',
  required: ['agent', 'entries', 'unverified', 'notes'],
  properties: {
    agent: { type: 'string' },
    entries: { type: 'array', items: { type: 'object', required: ['topic', 'content', 'source'], properties: { topic: { type: 'string' }, content: { type: 'string' }, source: { type: 'string' }, confidence: { type: 'string', enum: ['VERIFIED', 'SNIPPET-ONLY'] } } } },
    unverified: { type: 'array', items: { type: 'string' } },
    notes: { type: 'string' },
  },
}

phase('Context')
const results = await parallel([
  () => agent(`${COMMON}

ROLE: COLOR-SECTOR ARCHAEOLOGIST (the O2-archaeology pattern: exhaustive, classified, file:line catalog).
Sweep: ${ROOT}/corpus3/01-GUM-Omega-Paper-v4.3-ext.md Sec. VII.J + everywhere color/confinement appears (grep for: color, confin, tube, disclination, junction, quark, Skyrme, mismatch, winding, stratum, phason, censorship, hadron); ${ROOT}/corpus2/ (same greps — the dated record may print more detail than v4.3 carries); ${ROOT}/theory-audit/ (esp. any files auditing VII.J, h25 knot rows, sector content); ${ROOT}/tier5-family/ (the r11 belt/tube reconstruction).
DELIVER: (1) the complete catalog of printed color-sector claims with ledger classes (Q-1, Q-2, Q-3, Q-5, Q-6-prime, sigma = 0.19 GeV^2, f_q, M_gap value if printed, epsilon_dress, f-top ~60 MeV, phason-string window, Y-law junction, WS1-H3 tetraquark kill, F-Q1, F-Q2-prime, quark-spacing rule); (2) THE CENTRAL QUESTION: does the corpus anywhere print, imply, accommodate, or forbid a KNOT-FREE bound state of the pure tube/locking-stratum sector (a closed disclination loop, closed tube, glueball analog, pure-glue state, gluon analog)? Search hard: closed loop, loop, torus, ring, breather, vacuum-web excitation, stratum excitation, gap mode. Report every near-miss with evidence, and SILENT verdicts where silent. (3) What the corpus says the GLUON maps to (the mediator of the gapped stratum — is there an explicit gluon dictionary entry?). (4) What 'dual-constrained by sigma and the hadronic spectrum' (F-Q2-prime) concretely references.`,
    { label: 'ctx:color-archaeology', phase: 'Context', schema: CORPUS_SCHEMA }),

  () => agent(`${COMMON}

ROLE: MASS-ANCHOR READER.
Read: ${ROOT}/corpus3/01-GUM-Omega-Paper-v4.3-ext.md Secs. IV (hbar closure, Lambda calibration), V (if scales appear), VII.A-VII.E (band-edge mass theorem, knot profile, virial, pair law), VII.J; grep for: M_gap, m_Sk (or frak-m_Sk), Lambda, band edge, mass scale, GeV, calibration.
DELIVER: (1) the complete chain of how ABSOLUTE mass/energy scales are set in the model (what is calibrated [CAL], what imported [IM], what derived [DF/DW]); (2) every printed number with GeV units relevant to the strong sector (sigma = 0.19 GeV^2 - where does it come from, is it imported from lattice/Regge phenomenology?; f_q = 0.14-0.20 GeV; m_Sk value; M_gap; f-top ~ 60 MeV; the phason window); (3) whether the model has any printed statement about the HADRONIC spectrum it must reproduce (the F-Q2-prime dual constraint) - which hadron masses are anchors; (4) what a 2.37 GeV state would be in the model's natural units (ratios to sqrt(sigma) = 436 MeV, to m_Sk, to f_q) - pure arithmetic from printed numbers, no new physics.`,
    { label: 'ctx:mass-anchors', phase: 'Context', schema: CORPUS_SCHEMA }),

  () => agent(`${COMMON}

ROLE: QUANTUM-NUMBERS / SELECTION-RULES READER.
Read: ${ROOT}/corpus3/01-GUM-Omega-Paper-v4.3-ext.md VII.E-VII.G (knot matter, quantization, Theorem c-doubleprime, Finkelstein-Rubinstein), VII.J, IV.I (spin-selection), plus ${ROOT}/theory-audit/i4_topology.md and h4_completion.md (the exchange/rotation monodromy machinery).
DELIVER: (1) how composite states get quantum numbers in the model (J^PC of line-neutral composites; is there a printed parity/C-parity operation for tube states?); (2) the flavor structure: what distinguishes K-Kbar mesons from other composites - do knots carry flavor labels, and is there a printed flavor symmetry (SU(3)-flavor analog) of the emergent Skyrme identification?; (3) is there any printed G-parity or generalized-G-parity analog, or any selection rule forbidding specific composite decays?; (4) for a hypothetical KNOT-FREE closed tube state: from the printed quantization machinery (FR constraints, chirality of the substrate, the locking-stratum orientation field), what J^PC would the model's own rules assign to the lightest such state - answer ONLY from printed machinery, flag every step that requires new assumptions as [CJ-new]; (5) does the substrate's CHIRALITY print any consequence for parity assignments of emergent states (the chiral condensate, double-twist network handedness)?`,
    { label: 'ctx:selection-rules', phase: 'Context', schema: CORPUS_SCHEMA }),

  () => agent(`You are a literature agent in Tier 9 of the GUM campaign. FIRST READ ${ROOT}/tier9-glueball/GLUEBALL_PAPER_DIGEST.md Sec. 7 (your verification targets). Your tools: WebSearch (works). WebFetch is BLOCKED for arxiv.org (EGRESS_BLOCKED) and possibly other domains - try non-arxiv domains (journals.aps.org, phys.org, inspirehep.net, en.wikipedia.org) but fall back to WebSearch snippets; label each entry VERIFIED (you fetched/read a full page confirming it) or SNIPPET-ONLY (search-result snippet). Any number you cannot pin gets listed in 'unverified'. Cite source URLs.

ROLE: GLUEBALL THEORY LITERATURE.
Targets: (1) lattice QCD glueball spectrum - lightest 0++ mass (~1.5-1.8 GeV), 2++ (~2.4), 0-+ (~2.3-2.6 GeV): get the standard values from Morningstar-Peardon 1999, Chen et al 2006, Gui et al 2019 (quenched 0-+ ~2.395 GeV?), any unquenched/full-QCD values, and the RATIO 0-+/0++; (2) the Isgur-Paton flux-tube model of glueballs as closed flux tubes (Phys Rev D 31, 2910 (1985)): its predicted spectrum, what J^PC the lowest closed-tube states have, the phonon quantization on the loop, and any modern closed-flux-tube glueball spectra (e.g. Johnson-Teper, torelon studies, effective-string corrections); the closed-string mass formula m^2 or m vs sqrt(sigma) relations; (3) glueball-qqbar and glueball-ccbar mixing (Zhang et al PLB 827 136960; Chen-Gui-Li-Sun 2026 mixing angle 2-5 degrees); (4) the sqrt-OZI rule (Robson 1977); (5) scalar glueball candidates f0(1710)/f0(1500) status and the 2025-2026 'lattice evidence that scalar glueballs are small' (arXiv:2508.21821) if findable; (6) string tension sigma ~ 0.18-0.19 GeV^2 standard value provenance (Regge slope / lattice).`,
    { label: 'lit:glueball-theory', phase: 'Context', schema: LIT_SCHEMA }),

  () => agent(`You are a literature agent in Tier 9 of the GUM campaign. FIRST READ ${ROOT}/tier9-glueball/GLUEBALL_PAPER_DIGEST.md Secs. 3 and 5. Tools: WebSearch (works); WebFetch blocked for arxiv.org, try other domains, label VERIFIED vs SNIPPET-ONLY; unverifiable numbers to 'unverified'; cite URLs.

ROLE: X(2370) EXPERIMENT + ALTERNATIVES LITERATURE.
Targets: (1) the X(2370) measurement chain: 2011 discovery (PRL 106 072002), 2020 confirmation (EPJC 80 746), 2024 spin-parity (PRL 132 181901), the 2026 combined mass/width 2359 MeV & 170 MeV source (arXiv:2605.26495 - K_S K_S pi0 / pi0 pi0 eta observation paper), and the August 2026 ICHEP announcement (phys.org, IHEP press release); (2) alternative interpretations in the literature and their current status: eta-eta-prime radial excitation (Yu-Sun-Liu-Zhao PRD 83 114021), the 15-200 MeV K*Kbar width prediction (Wang-Luo-Sun-Liu PRD 96 034013), baryonium, hybrid, 4-quark (Su-Chen PRD 106 014023, Dong et al EPJC 80 749, Wang PRD 112 034010); (3) any published skepticism or caveats about the glueball identification (2026 reactions, Science/AAAS coverage, community commentary - e.g. mixing-fraction unknowns, need for partial-wave analysis); (4) the eta(1405)/eta(1475) situation (one-state vs two-state, K*(892)K as major mode) as the contrast case; (5) prior glueball-candidate history for X(1835), f0(1500), f0(1710) - why earlier candidates stayed unconfirmed.`,
    { label: 'lit:x2370-alternatives', phase: 'Context', schema: LIT_SCHEMA }),

  () => agent(`${COMMON}

ROLE: STAKES-MACHINERY READER.
Read: ${ROOT}/corpus3/01-GUM-Omega-Paper-v4.3-ext.md Sec. IX (experimental program + stakes table v2), Sec. XII (How GUM Dies - the kill registry), Sec. I.D (ledgered results), plus grep for 'banked consistency', 'stake', 'S1' through 'S7', 'wager'.
DELIVER: (1) the exact machinery: what qualifies as a STAKE vs a BANKED CONSISTENCY vs a KILL CONDITION in the corpus's own language - quote the definitions; (2) the current stakes table content (what S1-S7 are); (3) how the corpus handles NEW external experimental results arriving after the edition (is there a printed protocol? the SWP-1 Standing Watch - what does it say?); (4) precedents: how were previous external-data confrontations graded (the family sector's PDG pulls, the Y-law 'banked consistency') - what grade language was used; (5) your assessment (flagged as coordinator-input, not corpus text): under the corpus's own rules, what is the CORRECT grade-language for a hypothetical Tier-9 finding that the corpus's closed-tube sector is/is-not consistent with the X(2370) glueball evidence - which categories are available without violating 'grades never rise by replication'?`,
    { label: 'ctx:stakes-machinery', phase: 'Context', schema: CORPUS_SCHEMA }),
])

return { colorArch: results[0], massAnchors: results[1], selectionRules: results[2], litTheory: results[3], litX2370: results[4], stakes: results[5] }