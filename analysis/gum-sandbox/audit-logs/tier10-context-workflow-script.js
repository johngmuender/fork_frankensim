export const meta = {
  name: 'tier10-U-context',
  description: 'Junction context: corpus archaeology + tube dynamics + field literature',
  phases: [
    { title: 'Context', detail: 'two corpus readers + two literature agents' },
  ],
}

const ROOT = '/home/user/fork_frankensim/analysis/gum-sandbox'

const COMMON = `You are a context agent in Tier 10 (Phase U) of the GUM replication campaign at ${ROOT}. The campaign is analyzing the GUM corpus (speculative chiral-micropolar-substrate theory; within-model register: nothing bears on nature) in the context of STAR's baryon-junction evidence (arXiv:2408.15441: baryon number traced by the Y-shaped gluonic junction, not valence quarks). REQUIRED FIRST READS: (1) ${ROOT}/tier10-junction/JUNCTION_PAPER_DIGEST.md (the coordinator's full digest incl. GUM hooks in Sec. 5); (2) ${ROOT}/tier9-glueball/GLUEBALL_CONTEXT_ANALYSIS.md (the Tier-9 baseline: color-sector archaeology, stakes-machinery constraints — REUSE its stakes conclusions, do not re-derive; seal q-theta binds all mass numbers to V.F grade). Cite every corpus claim with file path + line/section. Ledger classes [DF]/[DW]/[CAL]/[IM]/[CJ] faithful. Corpus silence = SILENT verdict, never inferred content.`

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

ROLE: JUNCTION ARCHAEOLOGIST (O2-archaeology pattern: exhaustive, classified, file:line catalog).
Sweep corpus/, corpus2/, corpus3/, theory-audit/, substrate-suite/ for: junction, Y-law, Y-junction, n-junction, WS1-H3, baryon, baryon number, reconnection, disclination line/network, mismatch tube meeting points, tetraquark, junction order.
DELIVER:
(1) THE CENTRAL QUESTION — the corpus-internal counterpart of the STAR title question: in GUM, do KNOTS or JUNCTIONS carry baryon number? Catalog everything printed about baryon number as a conserved quantity: is there a dictionary entry (B = ? in substrate variables)? What does the WS-K charge census say exactly (web-line charges, pi_1-classes, junction constraints, reconnection net-conserving)? Is baryon-number conservation derived, asserted, or SILENT? If B is only implicit, state precisely what IS printed and what reading it supports.
(2) Everything printed about the junction itself: junction order n (WS1-H3 — the Tier-9 sweep found the computation is NOT in the archive; confirm and catalog every reference to it), the Y-law banked consistency (exact language, all sites), junction energy/geometry (any printed junction-energy term?), the tetraquark kill's exact wording.
(3) Junction-antijunction / knot-free junction states: does the corpus anywhere print, imply, or forbid a J-Jbar composite (three tubes between two junctions, no knots) or larger junction networks (the STAR paper's 'baryonium glueballs, gluonic graphene, buckyballs' analogs)? Same ACCOMMODATED/FORBIDDEN/CONSTRUCTED/SILENT verdict discipline as Tier-9's closed-loop question. Note how Course 20.1's 'exactly' inventory sentence bears on this class.
(4) Reconnection: every printed statement about tube/line reconnection (WS-K 'reconnection events (net-conserving)' — net-conserving of WHAT exactly?); jets/snapping (Primer's cord-snapping = quark-pair minting) — the dynamics closest to junction transport.
(5) Anything printed about baryon asymmetry/baryogenesis in the corpus (the relic census? cosmology sector?) — the STAR paper frames B-conservation since baryogenesis; does GUM print any B-genesis story? Likely SILENT — verify.`,
    { label: 'ctx:junction-archaeology', phase: 'Context', schema: CORPUS_SCHEMA }),

  () => agent(`${COMMON}

ROLE: TUBE/JUNCTION DYNAMICS + TRANSPORT READER.
Read: corpus3/01-GUM-Omega-Paper-v4.3-ext.md VII.J + Course Ch. 20 + Primer Ch. 15 (the color sector's dynamics content); the campaign tube machinery (${ROOT}/theory-audit/h21_RESULTS.md, ${ROOT}/tier7-program/M3/RESULTS.md, ${ROOT}/tier7-program/N2/RESULTS.md); tier9-glueball/T1/RESULTS.md and T2/RESULTS.md (the validated closed-tube spectrum machinery and worldsheet census — Tier-10 workstreams will extend these); the mass-anchor inventory in tier9-glueball/t_context_agents.json (mass-anchor reader section).
DELIVER:
(1) What the corpus prints about DYNAMICS of tubes/junctions: propagation, mobility, tension vs mass scales, anything bearing on 'junctions are light/slow-momentum, knots are heavy' — the GUM analog of the STAR transport asymmetry. Expected mostly SILENT; catalog near-misses (e.g. string-breaking/minting, censorship dressing floor, the knot mass scales m-tilde/m_Sk vs tube scale sqrt(sigma)/f_q).
(2) The printed mass-scale asymmetry: knot masses (m-tilde-heavy species; m_Sk = 1.7 GeV floor) vs tube/junction scales (sqrt(sigma) = 436 MeV, f_q = 0.14-0.20 GeV) — assemble the complete printed-number basis for a two-carrier mobility toy (Tier-10 U2 will pre-register it; your job is the INPUT inventory with citations, not the model).
(3) What T1's validated machinery provides for a junction-network extension (U1): the IP/NG implementations, the frozen f* = 1.9724 cutoff, the Y-law geometry needed for a J-Jbar dumbbell (three tubes + 2 junctions) — what additional structure is needed and what the corpus prints for it (junction energy: printed or [CJ-new]?).
(4) Any printed statement on EIC-adjacent observables (u-channel, backward production, baryon distribution in nuclei) — expected SILENT, verify.`,
    { label: 'ctx:dynamics-transport', phase: 'Context', schema: CORPUS_SCHEMA }),

  () => agent(`You are a literature agent in Tier 10 of the GUM campaign. FIRST READ ${ROOT}/tier10-junction/JUNCTION_PAPER_DIGEST.md Sec. 6 (your verification targets). Tools: WebSearch (works); WebFetch is BLOCKED for arxiv.org — try non-arxiv domains (journals.aps.org, inspirehep.net, phys.org, link.springer.com), label each entry VERIFIED (full page read) or SNIPPET-ONLY; unpinnable numbers into 'unverified'; cite URLs.

ROLE: JUNCTION THEORY LITERATURE.
Targets: (1) the founding papers: Artru NPB 85, 442 (1975) (string model with baryons: topology); Rossi-Veneziano NPB 123, 507 (1977) and the 2016 update JHEP 06, 041 (string-junction picture of multiquark states); (2) Kharzeev PLB 378, 238 (1996) 'Can gluons trace baryon number?' — pin the Regge intercepts behind 0.42 < alpha_B < 1 (J+J vs J+Pomeron) and the exp(-alpha_B Delta-y) transport law; (3) Csorgo-Gyulassy-Kharzeev J. Phys. G 30, L17 (2004) 'Buckyballs and gluon junction networks on the femtometer scale' — WHAT states do they construct (junction-network buckyballs/graphene), what masses/scales do they print (this is the [IM] anchor set for GUM's U1 junction-network spectrum); (4) lattice Y-law: Takahashi et al. PRL 86, 18 (2001), Suganuma et al. AIP Conf. Proc. 756 (2005), Bissey et al. PRD 76, 114512 (2007) — Y vs Delta ansatz status, junction-energy terms if quantified; (5) PYTHIA junction machinery: Christiansen-Skands JHEP 08, 003 (2015) (color reconnection beyond leading colour, junction formation), Lonnblad-Shah EPJC 83, 1105 (2023) (baryon correlations); (6) any theoretical estimates of junction-antijunction bound states (baryonium glueball masses) beyond CGK 2004.`,
    { label: 'lit:junction-theory', phase: 'Context', schema: LIT_SCHEMA }),

  () => agent(`You are a literature agent in Tier 10 of the GUM campaign. FIRST READ ${ROOT}/tier10-junction/JUNCTION_PAPER_DIGEST.md Secs. 2-4. Tools: WebSearch (works); WebFetch blocked for arxiv.org, try other domains; VERIFIED vs SNIPPET-ONLY discipline; 'unverified' list; URLs.

ROLE: EXPERIMENT + ALTERNATIVES + RECEPTION LITERATURE.
Targets: (1) the STAR measurement chain: arXiv:2408.15441 publication status (journal? it was submitted 2024 — where did it land? Science? PRX?), the companion methodology paper Lewis et al. EPJC 84, 590 (2024) (arXiv:2205.05685), the isobar dataset provenance (PRC 105, 014901 CME search), prior baryon-stopping measurements (NA49 PRL 82, 2471; BRAHMS PLB 677, 267; STAR PRC 79, 034909); (2) RECEPTION AND CRITICISM 2024-2026: any published skepticism of the junction interpretation — valence-quark-side defenses, alternative explanations of the isobar B/Delta-Q (neutron-skin-only readings? hadronic rescattering?), responses to the gamma+Au slope; any theory papers post-2408.15441 testing or challenging it (e.g. Kovchegov, Vitev, Lv-et-al works on baryon stopping at small-x; 'baryon stopping as a probe' literature 2024-2026); (3) the EIC u-channel program (Cebra et al. PRC 106, 015204 (2022), arXiv:2204.07915) and the ePIC/EIC baryon-junction observables plans; (4) LHC-side evidence: LHCb/ALICE baryon-to-charge or baryon-stopping analyses bearing on junctions (e.g. ALICE antimatter/baryon transport at mid-rapidity at 13 TeV); (5) status of junction exotics searches (baryonium tetraquark candidates via Rossi-Veneziano 2016; any experimental limits).`,
    { label: 'lit:experiment-reception', phase: 'Context', schema: LIT_SCHEMA }),
])

return { junctionArch: results[0], dynamics: results[1], litTheory: results[2], litExp: results[3] }