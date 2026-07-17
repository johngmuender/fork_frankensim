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
