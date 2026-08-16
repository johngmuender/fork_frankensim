# CORPUS3 FOLD RECORD — coordinator verification and disposition

Run 2026-07-17, workflow wf_2c33ba17-f9b (four parallel document
editors) + two coordinator fixes. Within-model; nothing bears on nature.

## Verification (coordinator re-run, final state)

Additive-only check — `diff corpus2/<orig> corpus3/<edition>`:

| Edition | Deleted lines | All deletions are version headers? | Added lines | ⟦T6⟧ markers |
|---|---|---|---|---|
| 01 Omega paper v4.0-ext | 6 | ✓ (main title + 5 PART-N running heads) | 38 | 30 |
| 02 Course v3-ext | 2 | ✓ (title + subtitle) | 24 | 15 |
| 03 Primer v3-ext | 2 | ✓ | 14 | 8 |
| 04 Teachers' Ed. v3-ext | 3 | ✓ (title, subtitle, colophon line) | 19 | 6 |
| 08 Session Map v3.0-ext | 1 | ✓ | 24 | 12 |
| Watch-Mode memo v3-ext | 1 | ✓ | 10 | 5 |
| 09 Replication Record v2.0 | 1 | ✓ | 24 | 13 |
| README.manifest v5.0-ext | new file | n/a | 68 | 18 |

All pre-existing corpus2 markers (⟦rev⟧, ⟦H4 resolved⟧, ⟦I4 reported⟧,
⟦J…⟧) verbatim (agent diffs; count deltas explained by mentions inside
new ⟦T6⟧ text only). Every folded number sourced from
T6_FOLD_PAYLOAD.md; L7 carries [DW, testbed-grade,
configuration-specific] + the verbatim scope caveat at every citation;
L5′ G4 reported FAIL-AT-d=25/null-approached everywhere; all five
annotations labeled proposals.

## Coordinator fixes (post-agent, pre-commit)

1. **Gate-tally harmonization**: the paper's Sec. 0 tally ("25 gates /
   23 PASS", from the workflow instruction) and Sec. X tally ("16
   original gates…", from the payload) were reconcilable but not
   face-identical — both replaced by the unambiguous per-workstream
   form (L1 4/4, L2 5/5, L3 5/5, L5 3+2 PARTIAL→converted, L7a 4/4,
   L7b 5/5; sole gated non-PASS = L5′ G4 FAIL-AT-d=25). Payload P7
   harmonized identically. [The instruction-vs-payload mismatch was the
   coordinator's — filed in the own-defects tradition.]
2. **Running heads**: the five PART-N-OF-6 heads bumped to v4.0-ext
   (Tier-6 fold) under the version-header exception (paper agent's
   caveat #1 accepted).

## Disposition

corpus3 is the replicators' proposed Tier-6 fold edition set, offered
to the corpus's authors per REVISION_CHARTER_v2; adoption is theirs;
corpus2 stands untouched as the dated record. Entry point:
README.manifest-v5.0-ext.md.

---

## ⟦EW⟧⟦T8⟧⟦T9⟧⟦T10⟧ CORPUS4 FOLD — header entry (extended additively; the corpus3 record above stands verbatim)

**The fold.** Run 2026-08-16, **commissioned by the corpus's owner**
(session goal, 2026-08-16). Governing charter:
corpus4/REVISION_CHARTER_v3.md (binding rules 1–8, inheriting the
corpus3 charter v2). Payload of record: corpus4/EW_SS_FOLD_PAYLOAD.md
(Parts A–D; charter rule 7 — editors copy from the payload's cited
source blocks, never restate from memory). Sources folded: the
synthesis line — synthesis/GUM-Omega-Synthesis-v4.4-EW.md §8, carried
verbatim through synthesis/GUM-Omega-Synthesis-v4.5-SS.md (the EW and
strong source text; ⟦EW⟧ content adopted by owner commission, corpus4
charter v3 rule 5) — plus the Tier-8/9/10 campaign phases (findings
F-T8-S1…S4, F-T9-T1…T4, F-T10-U1…U3; adjudications in
tier8-assessment/, tier9-glueball/, tier10-junction/; each ⟦T8⟧/⟦T9⟧/
⟦T10⟧ box a campaign proposal, adopted this edition). Executed as the
corpus4 fold workflow (parallel document editors over the eleven
corpus4 documents; workflow IDs, editor assignments, and audit-archive
paths recorded by the coordinator at verification). Within-model;
nothing bears on nature.

## Per-document verification templates (charter v3 rule 8, gates (a)–(e); verdicts filled by the coordinator at verification)

Gate key: (a) diff vs corpus3 base shows deletions ONLY in version
headers/supersession notes · (b) line count strictly grows · (c) marker
counts reported · (d) payload spot-checks verbatim · (e) register scan
(no forbidden sentences; V.F language on masses).

| corpus4 document | corpus3 base | (a) additive-only diff | (b) line growth | (c) marker counts | (d) payload spot-checks | (e) register scan |
|---|---|---|---|---|---|---|
| 01-GUM-Omega-Paper-v5.0-ext.md | 01-GUM-Omega-Paper-v4.3-ext.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| 02-The-Substrate-Course-v4-ext.md | 02-The-Substrate-Course-v3-ext.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| 03-The-Substrate-Primer-v4-ext.md | 03-The-Substrate-Primer-v3-ext.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| 04-Substrate-Primer-Teachers-Edition-v4-ext.md | 04-Substrate-Primer-Teachers-Edition-v3-ext.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| 08-GUM-Session-Map-v4.0-ext.md | 08-GUM-Session-Map-v3.0-ext.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| 09-External-Replication-Record-v3.0.md | 09-External-Replication-Record-v2.0.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| The-Watch-Mode-Transition-Memo-v4-ext.md | The-Watch-Mode-Transition-Memo-v3-ext.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| README.manifest-v6.0-ext.md | README.manifest-v5.0-ext.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| T6_FOLD_PAYLOAD.md (3-line ⟦T10-era note⟧ header; body byte-untouched) | T6_FOLD_PAYLOAD.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| T7_FOLD_PAYLOAD.md (3-line ⟦T10-era note⟧ header; body byte-untouched) | T7_FOLD_PAYLOAD.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |
| FOLD_RECORD.md (this document, extended additively) | FOLD_RECORD.md | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) | PASS (coordinator-verified 2026-08-16) |

**Disposition (corpus4).** corpus4 is the owner-commissioned EW+SS fold
edition set: the corpus3 files stand untouched as the dated record;
gate verdicts (a)–(e) above are the coordinator's to fill; entry point:
README.manifest-v6.0-ext.md. Within-model; nothing bears on nature.

**Coordinator verification disposition (2026-08-16).**  Independent
re-verification of all 11 documents (diff vs corpus3 base; growth;
deletion-class audit): every document strictly grew (total 4,678 →
5,300 lines, +622); deletions are version-header/title class only
(6/2/2/3/1/1/3/1/0/0/0 per document, zero suspicious); T6/T7 payload
bodies byte-identical below their 3-line headers; corpus4 census = 13
files (11 editions + charter v3 + EW_SS payload).  Editor gate reports
(a)–(e) accepted as filed: E1 paper 98/98 verbatim-copy audit, 14 spot
strings; E2 teaching docs 20/20 + 9+ + 13 spot strings with
negative-control discipline, chapter-collision amendments printed
(Ch. 20-W / Ch. 15½ / §15.7-kept, charter-rule-3-forced); E3 meta set
including the workflow-ID and phase-date disclosures; E4 record 19/19
spot strings with rivals-at-parity language verified.  Fold workflow:
wf_1cadd65b-be9 (4 editors).  The corpus4 edition set is complete,
additive, and internally consistent; corpus3 stands untouched as the
dated record.  Within-model; nothing here bears on nature.
