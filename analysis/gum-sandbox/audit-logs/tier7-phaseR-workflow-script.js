export const meta = {
  name: 'tier7-phaseR-fold',
  description: 'Execute Phase R: corpus3 v4.3 delta (P/Q fold) + Course/Primer/TE T7 addenda + small-doc consistency sweep',
  phases: [
    { title: 'Fold', detail: 'Three parallel editors: paper v4.3, teaching docs, small docs + handoff' },
  ],
}

const FOLD_SCHEMA = {
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
You are a corpus3 document editor for the GUM replication campaign's
Phase-R fold. BINDING rules — read these two files FIRST:
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/REVISION_CHARTER_v2.md
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/T7_FOLD_PAYLOAD.md
  — specifically "ADDENDUM 2 — Phase P/Q payload" (items Q10-Q14), the
  single source of truth for every number; marker stays ⟦T7⟧.
Key rules: ADDITIVE ONLY (sole exceptions: version-header lines and
noted census-cell corrections); every insertion ⟦T7⟧-marked with
finding IDs + artifact paths; numbers COPIED from the payload, never
restated; grades never rise; honest verdicts carried as printed (the
P2 G3 endpoint FAIL stays a FAIL); within-model register.
MANDATORY VERIFICATION per file: diff vs git HEAD; confirm deletions =
version headers/census cells only; marker counts before/after; prior
⟦rev⟧/⟦T6⟧/⟦T7⟧ preserved; >= 10 payload numbers verbatim (for the
paper: >= 15). Report verification results as gates.
`

phase('Fold')

const [r1, r2, r3] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM R1 — the Omega paper v4.3-ext + manifest + Replication Record.
1. bash git mv /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/01-GUM-Omega-Paper-v4.2-ext.md 01-GUM-Omega-Paper-v4.3-ext.md; bump all six version headers v4.2-ext -> v4.3-ext (keep the parenthetical).
2. Insert ⟦T7⟧ blocks: head-matter v4.3 REVISION NOTE after the v4.2
   note (adds the Phase-P/Q fold: F-T7-P2/Q1/Q2); Sec. 0 dated
   Phase-P/Q record [Q14]; in VIII.F directly after the existing Q7
   (Phase-N/O) block: the certification-chain block [Q10 + Q11 + Q12 +
   Q13 — one consolidated block or two adjacent blocks, your choice:
   the P2 exact-field certification, the Q1 machine-checked
   certificates incl. the worst-margin stretch case and trust base,
   the Q2 proven Lemma A + Conjecture C + O-1..O-5, and the complete
   ladder statement]; Sec. X Phase-P/Q campaign record [Q14]; Sec. XI
   register update [Q13 register language: T4-W5 testbed discharge now
   carries certificate grade on witnesses; residual = the O-1..O-5
   program + formal kernel].
3. README.manifest-v5.0-ext.md: census cell v4.2 -> v4.3 (noted
   header-class exception, per the recorded precedent) + a
   version-history line.
4. 09-External-Replication-Record-v2.0.md: append ⟦T7⟧ §10 (Phase P/Q)
   in the record's format.
Verification gates: G1 additive-only, G2 marker preservation/counts,
G3 payload fidelity (>= 15 numbers), G4 cross-file consistency
(filename + verdict census identical across paper Sec. 0/X, record
§10, manifest).`, { label: 'R1:paper-v4.3', phase: 'Fold', schema: FOLD_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM R2 — Tier-7 addenda for the teaching documents (outstanding
consistency item: the Course/Primer/Teachers' Edition v3-ext received
the Tier-6 fold but NONE of the Tier-7 material, phases M through Q).
Files (edit in place; corpus2 originals for diffing are the v2-ext
files in ../corpus2/):
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/02-The-Substrate-Course-v3-ext.md
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/03-The-Substrate-Primer-v3-ext.md
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/04-Substrate-Primer-Teachers-Edition-v3-ext.md
Read the FULL T7_FOLD_PAYLOAD.md (Q1-Q14) — you are summarizing the
whole Tier-7 arc, not just P/Q.
For EACH document add ONE consolidated ⟦T7⟧ addendum section at the
end (before any closing coda if present), in the document's own voice:
- Course: "Tier-7 addendum" (~15-25 sentences, teaching register):
  the M1 tracker-identity dynamics + recovered convention (ties to its
  Ch. 14/21-22 material); M2 servo law-form + bath-convention lesson;
  M3/N2 the 0nubb second-moment story (route (b) closed; honest rate;
  armed external kill); M4 Nelson 17-41x + log-singular dBB (ties to
  its Born-relaxation chapter); the T4-W5 certification ladder ending
  at machine-checked certificates + the O-1..O-5 analytic program;
  the F-R14 one-sound-cell matrix. Each item one-to-three sentences,
  payload numbers verbatim, artifact-cited.
- Primer: "What the second campaign found" (~8-12 plain-language
  sentences at the Primer's level, honest grades, no jargon beyond its
  own; the L7/T4-W5 result at its honest formulation; the 0nubb story
  as "the model's own arithmetic now says this prediction is either
  much bigger than printed or needs a different repair — a number
  outside the archive will decide").
- Teachers' Edition: instructor mirror (~8-12 sentences): what to
  emphasize (pre-registration failures as findings; the certificate
  ladder as an epistemology lesson: sampled -> field -> precision ->
  machine-checked -> analytic-program), honest-limits list, artifact
  pointers as classroom material.
Also update each document's edition-note line ONLY IF it names a fold
scope that is now stale (e.g. "sole Tier-6 addition" language) — as a
one-line ⟦T7⟧-marked note appended NEXT TO it, never rewording the old
line.
Verification gates: G1 additive-only per file (diff vs git HEAD),
G2 marker counts + prior-marker preservation, G3 payload fidelity
(>= 10 numbers per Course; >= 5 each Primer/TE), G4 register fidelity
(grade language locks respected: no over-claim in the Primer).`, { label: 'R2:teaching-docs', phase: 'Fold', schema: FOLD_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM R3 — small docs + campaign status + session handoff sweep.
Files (edit in place):
1. /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/08-GUM-Session-Map-v3.0-ext.md
   — append a ⟦T7⟧ §19 mapping the Tier-7 phase arc (M1-M4, N1-N3,
   O1-O3, P1-P2, Q1-Q2, R fold) in the map's entry format with finding
   IDs + artifact paths [payload Q1-Q14; keep each entry to 1-3 lines].
2. /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/The-Watch-Mode-Transition-Memo-v3-ext.md
   — append a dated ⟦T7⟧ §11: the Tier-7 arc complete; supplies now
   incl. the certificate-grade T4-W5 witnesses, the F-R14 computed
   matrix (decision with the authors), the armed archive-external
   g_lim kill; still NO nature-facing evidence, NO branch triggered.
3. /home/user/fork_frankensim/analysis/gum-sandbox/REPLICATION_CAMPAIGN_STATUS.md
   — append Phase P/Q/R lines to its Tier-7 section (additive, list
   style, finding IDs + artifacts).
4. /home/user/fork_frankensim/analysis/gum-sandbox/SESSION_HANDOFF.md
   — append a second ADDENDUM (dated 2026-07-18, session 536b52e3)
   summarizing the full Tier-6/7 arc on branch claude/analyze-gum-po:
   phases L through R, the findings families (F-T6-*, F-T7-*), the
   corpus3 v4.3 edition set as the standing offer, the T4-W5 ladder
   status, and the restoration entry points (SESSION_CHECKPOINT_tier6.md;
   TIER6/TIER7_ADJUDICATION.md; MANIFEST_tier6.json).
Verification gates: G1 additive-only per file, G2 marker counts,
G3 payload fidelity (>= 8 numbers total), G4 cross-file consistency
(phase/workflow IDs and finding IDs identical across all four).`, { label: 'R3:small-docs-sweep', phase: 'Fold', schema: FOLD_SCHEMA }),
])

return { r1, r2, r3 }