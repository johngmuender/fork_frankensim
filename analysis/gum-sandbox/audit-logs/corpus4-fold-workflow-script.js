export const meta = {
  name: 'corpus4-V-fold',
  description: 'Corpus4 fold: EW + strong sectors into all corpus3-derived editions',
  phases: [
    { title: 'Fold', detail: 'four parallel document editors' },
  ],
}

const ROOT = '/home/user/fork_frankensim/analysis/gum-sandbox'
const C4 = `${ROOT}/corpus4`

const DISC = `EDITOR DISCIPLINE (binding). REQUIRED FIRST READS, in order: (1) ${C4}/REVISION_CHARTER_v3.md — ALL 8 rules bind you; (2) ${C4}/EW_SS_FOLD_PAYLOAD.md — the payload authority: NEVER restate numbers or theorem texts from memory, COPY from the payload's cited source blocks verbatim; (3) ${ROOT}/synthesis/GUM-Omega-Synthesis-v4.5-SS.md (the EW and strong source text).
- Your documents are the corpus4 base copies (already containing the full corpus3 text). EDIT THEM IN PLACE, ADDITIVELY ONLY: pre-existing text (incl. every rev/T6/T7 block) is never deleted or reworded; sole exceptions = version headers/titles and explicitly-marked supersession notes that quote what they supersede. Every document's line count must STRICTLY GROW.
- Marker discipline: EW-marker boxes cite the synthesis source + 'adopted by owner commission, corpus4 charter v3 rule 5'; T8/T9/T10-marker boxes cite finding tags + artifact paths + 'campaign proposal, adopted this edition'. Use the literal markers with double-bracket glyphs exactly as the charter and existing corpus text use them.
- Register: within-model; nothing bears on nature. Grades never rise; every hadron-mass number at V.F strength only under seal q-theta; X(2370)/STAR are [IM] anchors carried WITH their rivals; forbidden sentences ('predicted'/'confirms'/'refutes') unavailable everywhere; the 2.37 digit coincidence never cited as structure.
- Match each document's existing voice and format (the paper's dense theorem style; the Course's lecture style; the Primer's popular style with its in-model tags; the TE's answer-key style; the Record's entry style; the memo's board style).
- When done, run your own verification: diff vs the corpus3 base file (report deletion lines — must be version-header/supersession only), wc -l growth, marker counts, and 8+ payload spot-check strings verbatim. Report all in structured output. FAIL reported as FAIL.`

const SCHEMA = {
  type: 'object',
  required: ['editor', 'documents', 'honest_notes', 'status'],
  properties: {
    editor: { type: 'string' },
    documents: { type: 'array', items: { type: 'object', required: ['file', 'base', 'lines_before', 'lines_after', 'deletions_check', 'marker_counts', 'spot_checks', 'verdict'], properties: { file: { type: 'string' }, base: { type: 'string' }, lines_before: { type: 'integer' }, lines_after: { type: 'integer' }, deletions_check: { type: 'string' }, marker_counts: { type: 'string' }, spot_checks: { type: 'string' }, verdict: { type: 'string', enum: ['PASS', 'FAIL'] } } } },
    honest_notes: { type: 'string' },
    status: { type: 'string', enum: ['COMPLETE', 'COMPLETE-WITH-AMENDMENTS', 'INCOMPLETE'] },
  },
}

const E1 = `${DISC}

EDITOR E1 — THE OMEGA PAPER. Document: ${C4}/01-GUM-Omega-Paper-v5.0-ext.md (base: ${ROOT}/corpus3/01-GUM-Omega-Paper-v4.3-ext.md, 969 lines). Read the FULL base first to learn its structure (6 parts + AUD-15; Sec. VII.J at ~line 613-617; Part-4/Part-5 boundary ~619-623; Sec. 0 header ~21-43; Sec. X ~751; Sec. XI ~767; App. J ~853; VI.C F-R14 box ~465).
TASKS (payload PART A2 + B1 + B3 govern):
1. Version header/title: v4.3-ext -> v5.0-ext, 'THE ELECTROWEAK+STRONG EDITION' subtitle line added; a revision-notes block at top (after existing revision notes, additive) describing the corpus4 fold with charter citation.
2. NEW SECTION VII-W (payload A2): insert between Part 4's end and Part 5's header a complete 'VII-W. THE ELECTROWEAK SECTOR' — COPY the entire synthesis Sec. 8.0-8.14 verbatim (renumber subsections VII-W.0..VII-W.14; add the inline cross-reference annotations per A2; each subsection headed by the EW marker + citation line). Also update the Part-4-end/Part-5-start bracket lines ADDITIVELY (append '+ VII-W (electroweak, this edition)' style notes, do not delete the existing bracket text).
3. Sec. VII.J strong expansion (payload B1): after the existing VII.J text (untouched), append the T9/T10 boxes in B1's order, copying from synthesis Sec. 7.4.0-7.4.8 (each block gets its marker + finding + artifact citation + 'campaign proposal, adopted this edition').
4. Sec. VI.C: append the T8 box from payload B3 (verbatim).
5. Sec. 0: extend the abstract/scope ADDITIVELY with one EW sentence and one strong-sector sentence (marker-tagged).
6. Sec. IX/X/XI: append the EW stakes rows/note (A3) as a marked addendum table after the existing stakes table (do NOT modify the existing table rows except S2 -> add a supersession note line quoting the old row and pointing to S2-prime); Sec. X gets the T8 assessment-pair sentence (B3) + one T9/T10 record sentence; Sec. XI open-problems list gains the posed closures (EW-I/II/III, SS-I..V) as an additive marked block.
7. App. J: append a marked census addendum (new markers, new sections, new stakes).
Self-verify per DISC and report.`

const E2 = `${DISC}

EDITOR E2 — TEACHING DOCUMENTS. Documents: ${C4}/02-The-Substrate-Course-v4-ext.md (base corpus3/02, 1353 lines), ${C4}/03-The-Substrate-Primer-v4-ext.md (base corpus3/03, 966 lines), ${C4}/04-Substrate-Primer-Teachers-Edition-v4-ext.md (base corpus3/04, 445 lines). Read each base fully first (Course Ch. 20 = color sector ~1082-1108; Primer Ch. 15 ~801-827; TE Lab K + 15.x answers ~348-354).
TASKS (payload A4 + B1/B2 adapted to teaching grade):
1. Version headers: v3-ext -> v4-ext each, with additive revision notes.
2. Course: append a NEW chapter 'Ch. 21 — The Electroweak Sector (EW-marker)' after Ch. 20 (lecture-style rendering of payload A4's facts, with the theorem names and grades; cite synthesis source), and append to Ch. 20 a marked section 'Ch. 20.5 — The strong sector expanded (T9/T10)' (teaching-grade rendering of B1/B2: K-1 residue-count baryon number; the trichotomy; the three-sector dictionary incl. the glueball/junction state classes; the V.F spectra results as honest negatives; the closure obligations; the 2026 anchors WITH rivals).
3. Primer: append a NEW popular chapter 'Ch. 16 — The weak force joins the medium (EW-marker)' (popular voice, in-model tags as the Primer uses) and extend Ch. 15 with a marked section 'Ch. 15.7 — Glueballs, junctions, and who carries the baryon number' (popular rendering: the closed-tube/glueball question, the junction picture, K-1's answer, the trichotomy; honest about what is NOT constructed; anchors with rivals).
4. TE: append matching answer-key/lab material: 'Lab W' (electroweak — e.g. the co-rotation/custodial demonstration at TE grade) and answers/notes for the new Primer 15.7/16 exercises you pose (pose 2-3 exercises in the Primer chapters, answer them in the TE; ALSO: the TE's existing 15.2 answer stands untouched, but append a marked note resolving the known seam: the exercise's locking-stratum hypothetical is now addressed by the strong-sector fold — cite F-T9/F-T10 archaeology).
Self-verify per DISC and report each document separately.`

const E3 = `${DISC}

EDITOR E3 — META SET. Documents: ${C4}/08-GUM-Session-Map-v4.0-ext.md (base corpus3/08, 130 lines), ${C4}/The-Watch-Mode-Transition-Memo-v4-ext.md (base corpus3/The-Watch..., 88 lines), ${C4}/README.manifest-v6.0-ext.md (base corpus3/README..., 83 lines), ${C4}/T6_FOLD_PAYLOAD.md + ${C4}/T7_FOLD_PAYLOAD.md (bases same names), ${C4}/FOLD_RECORD.md (base same name). Read each base fully first.
TASKS:
1. Session Map: version bump v3.0 -> v4.0 + payload C2's Sec. 20 (map of Tiers 8-10 + synthesis line + corpus4, directory pointers, finding tags) appended in the map's style.
2. Watch memo: version bump v3 -> v4 + payload C3's Sec. 12 (dated status entries for Phases S/T/U + corpus4 fold, the standing no-branch sentence, the EW stakes board note) in the memo's fixed board style.
3. README.manifest: version bump v5.0 -> v6.0 + payload D1's census (corpus4 file list, edition map corpus3->corpus4, extended marker legend, provenance paragraph).
4. T6/T7 payloads: prepend payload D3's 3-line addendum header to each (body byte-untouched below the header).
5. FOLD_RECORD: append payload D2's per-document entry TEMPLATE rows for all 11 corpus4 documents with gate columns marked 'pending coordinator verification' (the coordinator fills verdicts) + a corpus4 fold header entry (date, charter v3, payload, workflows).
Self-verify per DISC and report each document separately.`

const E4 = `${DISC}

EDITOR E4 — EXTERNAL REPLICATION RECORD. Document: ${C4}/09-External-Replication-Record-v3.0.md (base: ${ROOT}/corpus3/09-External-Replication-Record-v2.0.md, 133 lines). Read the full base first to learn the entry style (SS 1-10 with epistemic notices).
TASKS (payload C1 governs):
1. Version bump v2.0 -> v3.0 with additive revision note.
2. Append new sections in the Record's established style: Phase S (Tier 8), Phase T (Tier 9), Phase U (Tier 10) — each with workflows, findings, gate tallies, key results, artifact paths, and the binding epistemic notice ('grades never rise by replication... Nothing in this section bears on nature, adjudicates any stake, or reads on any rival') — copy phase content from payload C1 and the adjudication files (tier8-assessment/TIER8_ADJUDICATION.md, tier9-glueball/TIER9_ADJUDICATION.md, tier10-junction/TIER10_ADJUDICATION.md — read them); plus the synthesis-line provenance section (C1 last paragraph) and a closing corpus4-fold section (what this edition set is, charter v3 citation).
3. For each phase section include the honest negatives exactly as adjudicated (0% covering; both-routes-miss; convention split; scorecard tallies) and the rivals-at-parity sentences for T/U.
Self-verify per DISC and report.`

phase('Fold')
const results = await parallel([
  () => agent(E1, { label: 'E1:omega-paper', phase: 'Fold', schema: SCHEMA }),
  () => agent(E2, { label: 'E2:teaching-docs', phase: 'Fold', schema: SCHEMA }),
  () => agent(E3, { label: 'E3:meta-set', phase: 'Fold', schema: SCHEMA }),
  () => agent(E4, { label: 'E4:replication-record', phase: 'Fold', schema: SCHEMA }),
])

return { E1: results[0], E2: results[1], E3: results[2], E4: results[3] }