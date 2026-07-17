export const meta = {
  name: 'corpus3-t6-fold',
  description: 'Produce corpus3 editions: Omega paper v4.0-ext + Course/Primer/Teachers + small docs, folding the Tier-6 findings',
  phases: [
    { title: 'Fold', detail: 'Four parallel document editors applying T6 payload with additive-only marker discipline' },
  ],
}

const FOLD_SCHEMA = {
  type: 'object',
  required: ['document', 'sites', 'verification', 'summary'],
  properties: {
    document: { type: 'string' },
    sites: { type: 'array', items: { type: 'object', required: ['location', 'payload_items'], properties: { location: { type: 'string' }, payload_items: { type: 'string' }, note: { type: 'string' } } } },
    verification: { type: 'object', required: ['additive_only', 'marker_count', 'details'], properties: { additive_only: { type: 'string' }, marker_count: { type: 'number' }, details: { type: 'string' } } },
    files: { type: 'array', items: { type: 'string' } },
    caveats: { type: 'array', items: { type: 'string' } },
    summary: { type: 'string' },
  },
}

const COMMON = `
You are a corpus3 document editor for the GUM replication campaign's
Tier-6 fold. BINDING rules — read these two files FIRST and follow them
exactly:
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/REVISION_CHARTER_v2.md
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/T6_FOLD_PAYLOAD.md
Key rules restated: (1) ADDITIVE ONLY — never delete or reword existing
text (sole exceptions: the version header/title lines of your document,
which you update to the new edition name with a one-line edition note;
and explicitly-marked supersession notes that QUOTE what they
supersede). Pre-existing corpus2 markers (like the rev marker in double
white brackets) stay verbatim. (2) Every insertion is a block or inline
note beginning with the marker T6 in white double brackets: ⟦T6⟧ — and
cites its [F-T6-…] finding and tier6-foundations artifact path.
(3) COPY numbers from the payload file; never restate from memory.
(4) Grades never rise; proposals are labeled proposals; L7 carries its
testbed-grade/configuration-specific labels and scope caveat; L5-prime
G4 is always FAIL-AT-d25/null-approached. (5) The register: within-model;
nothing bears on nature.
VERIFICATION (mandatory, report results): after editing, run a
diff-based check against the corpus2 original file to confirm (a) every
changed hunk is an insertion of a ⟦T6⟧-marked block, a version-header
change, or a quoted supersession note — no other deletions; (b) count of
⟦T6⟧ markers; report both. Use bash diff/grep.
Your final structured output lists each insertion site (location =
section/anchor; payload_items = which P-items from the payload it
carries) and the verification results.
`

phase('Fold')

const [paper, course, primer, small] = await parallel([
  () => agent(`${COMMON}
DOCUMENT: the Omega paper v4.0-ext.
File (already copied from corpus2's v3.0-ext — edit in place):
/home/user/fork_frankensim/analysis/gum-sandbox/corpus3/01-GUM-Omega-Paper-v4.0-ext.md

This is the flagship. The file is large (~900 lines, very long lines);
Read it in windows and Grep for anchors. Make these insertions (each a
compact ⟦T6⟧ block in the corpus's own prose style — dense, ledgered,
artifact-cited; match the existing rev-block voice):

1. HEADER (line 1 region): retitle v3.0-ext → v4.0-ext, add "THE
   TIER-6 FOUNDATIONS FOLD" note + a REVISION NOTE (v4.0-ext) block
   directly after the existing v3.0-ext revision note (do not touch that
   note): states this edition = corpus3 fold of the Tier-6 foundations
   campaign (F-T6-1…6, L1/L2/L3/L5/L5-prime/L7a/L7b), binding charter =
   corpus3/REVISION_CHARTER_v2.md, marker = ⟦T6⟧, S-30 compliance
   (corpus2 files untouched; this is a NEW document), adoption is the
   authors'.
2. Sec. 0 (Public Record): dated ⟦T6⟧ paragraph: second external
   campaign phase (2026-07-17, Tier 6): foundations analysis vs the
   Bell/Bricmont/Goldstein/GRWm canon + 25 pre-registered gates across
   L1–L7 (23 PASS incl. both L5 PARTIAL conversions; one honest
   FAIL-AT-d25 with measured approach; three coordinator-owned spec
   defects printed). [P7]
3. Sec. I.D results table region: a ⟦T6⟧ note (do not restructure the
   table): rows 1–2's quantum-sector chain now carries executed
   canonical exhibits (L1 contextuality, L2 ToF incl. Nelson leg);
   row 19 Lock-2 flag upgraded to conservation-lemma-conditional
   (proposal #4); row 20's gate quantifier flag now carries the L7
   testbed discharge (proposal #5). [P2, P3, P4, P6]
4. Sec. III.A/III.C/III.D area: ⟦T6⟧ annotation: the two canonical
   "measurements don't measure" exhibits executed in-repo (L1 flipped
   SG: side preserved, label flipped, rho/j invariant; L2 ToF: exact
   Bricmont App.-1 + Nelson diffusion-blindness, first verification of
   Sec. III's actual stochastic kinematics against dBB phenomenology at
   equilibrium); F-T6-1 credit note (von Neumann-immunity nontrivial;
   T-B1 selection). [P1 credit, P2, P3]
5. Sec. VI.E: ⟦T6⟧ proposal #3 box: DGNSZ 2013 anchor + the two deltas
   (matter-sourced F; scheduled detection program). [P1/F-T6-4]
6. Sec. VIII.A: ⟦T6⟧ proposal #1 box: the architectural credit vs GRWm
   (medium-first PO hypersurface-independent trivially; the GRWm
   problem never arises) + the F-T6-3 discriminator ledger pointer +
   the NEW mutual kill (A3-positive + confirmed collapse contradict
   jointly). [P1/F-T6-3]
7. Sec. VIII.D Lock 2: ⟦T6⟧ proposal #4 box: the L3 conservation lemma
   with the measured exponent ladder (1.000000/2.999844/4.999869), the
   isolated source-form assumption, the proposed regrade
   [DF core + conservation-lemma-conditional]. Also annotate App. E.5
   (same lemma, one-line cross-pointer). [P4]
8. Sec. VIII.F: ⟦T6⟧ proposal #5 box: the L7 testbed discharge — the
   operational schema (fixed detector = one Naimark POVM), L7b's
   40.2-sigma non-affinity + affine-vanishing kill, L7a's
   backflow-borne quadraticity violation (1.84e-4 vs 2.9e-15 floor),
   status [DW, testbed-grade, configuration-specific], residual
   continuum-proof obligation, verbatim scope caveat. Also the
   collapse-class sharpening (A3-positive kills the (A4)/(A5) collapse
   family, not just orthodox proposals). [P6, P1/F-T6-3]
9. Sec. IX.A + IX.B: ⟦T6⟧ notes: A3 phenomenon instantiated in-repo
   (L5/L5-prime: tau_max exists, spin-dependent, tau_max(n_y) monotone
   209% contraction); the (k0, d_far) JOINT bench-design constraint;
   IX.B gate text now carries the L7 testbed citation with its grade
   labels. [P5, P6]
10. Sec. X (audit trail): ⟦T6⟧ dated paragraph: the Tier-6 campaign
    record incl. the three coordinator-owned spec defects (own-defects
    discipline maintained). [P7]
11. Sec. XI (open problems): ⟦T6⟧ register update: T4-W3
    discharged-at-toy-grade (residual: source-form verification);
    T4-W5 testbed-discharged [DW, configuration-specific] (new named
    open: continuum proof of the cutoff zeros + generality); the
    (k0, d_far) design note filed. [P4, P6, P5]
12. Sec. XII (How GUM dies): one ⟦T6⟧ line adding the mutual-kill
    cross-index (A3-positive + confirmed objective collapse cannot
    both stand). [P1/F-T6-3]
13. References section: ⟦T6⟧ reference-group extension per P7 (extend,
    do not duplicate — check which of Bell 1982/Das-Durr 2019 etc. are
    already present; GRWm 2014, DGNSZ 2013, Goldstein "Bell on Bohm",
    Bricmont 2019 are the likely additions) + the Tier-6 artifact-set
    pointer alongside the existing frankensim artifact-set entry.
Then run the mandatory diff verification vs
corpus2/01-GUM-Omega-Paper-v3.0-ext.md.`, { label: 'fold:paper-v4', phase: 'Fold', schema: FOLD_SCHEMA }),

  () => agent(`${COMMON}
DOCUMENT: The Substrate Course v3-ext.
File (already copied from corpus2's v2-ext — edit in place):
/home/user/fork_frankensim/analysis/gum-sandbox/corpus3/02-The-Substrate-Course-v3-ext.md
Corpus2 original for diffing:
/home/user/fork_frankensim/analysis/gum-sandbox/corpus2/02-The-Substrate-Course-v2-ext.md

The Course is the corpus's long-form teaching document (~1300 lines).
First skim its structure (grep for headings/section markers) to locate:
the quantum-sector/Nelson/Born chapters, any measurement or
Stern-Gerlach discussion, the one-world/foliation/nonlocality chapters,
the arrival-time/experimental-program material, and any Lock-2/elliptic
silence discussion. Then:
1. Update the title/header to v3-ext with a one-line edition note +
   a compact ⟦T6⟧ REVISION NOTE block near the top (charter pointer,
   marker discipline, adoption-is-theirs).
2. Insert ⟦T6⟧ teaching boxes at the matched sections — the Course is
   pedagogical, so each box should be a short "now executed:" exhibit
   note (2-6 sentences) pointing students at the runnable artifact:
   L1 flipped SG at the measurement/contextuality material [P2]; L2 ToF
   + Nelson diffusion-blindness at the Nelson/Born material [P3]; the
   GRWm comparison + mutual kill at the one-world/nonlocality material
   [P1/F-T6-3, F-T6-4]; L5/L5-prime + the (k0, d_far) design note at
   the arrival-time/experimental material [P5]; the L7 testbed at the
   gate/POVM material with grade labels + scope caveat [P6]; the L3
   lemma at any Lock-2/elliptic-silence material [P4]. If a natural
   anchor for an item does not exist in the Course, add at most ONE
   consolidated ⟦T6⟧ appendix-style section at the end carrying the
   remaining items, titled as a Tier-6 addendum.
3. Mandatory diff verification vs the corpus2 original.`, { label: 'fold:course-v3', phase: 'Fold', schema: FOLD_SCHEMA }),

  () => agent(`${COMMON}
DOCUMENTS (two): The Substrate Primer v3-ext AND the Teachers' Edition v3-ext.
Files (already copied from corpus2 — edit in place):
/home/user/fork_frankensim/analysis/gum-sandbox/corpus3/03-The-Substrate-Primer-v3-ext.md
/home/user/fork_frankensim/analysis/gum-sandbox/corpus3/04-Substrate-Primer-Teachers-Edition-v3-ext.md
Corpus2 originals for diffing:
/home/user/fork_frankensim/analysis/gum-sandbox/corpus2/03-The-Substrate-Primer-v2-ext.md
/home/user/fork_frankensim/analysis/gum-sandbox/corpus2/04-Substrate-Primer-Teachers-Edition-v2-ext.md

The Primer is the accessible exposition; the Teachers' Edition is its
instructor companion. Skim structure first (grep headings). For EACH:
1. Update title/header to v3-ext + compact ⟦T6⟧ REVISION NOTE.
2. Primer: insert brief reader-level ⟦T6⟧ notes (1-3 sentences each,
   plain language, no jargon beyond the Primer's own) at the matched
   spots: the measurement/observer discussion (L1: "the same particle,
   the same starting point, the flipped magnet — the answer label
   flips" [P2]); the momentum/uncertainty discussion (L2 [P3]); the
   nonlocality/other-theories discussion (GRWm comparison + the
   discriminating experiments [P1/F-T6-3]); the arrival-time experiment
   discussion (L5/L7 status at its honest grade [P5, P6]). Respect the
   Primer's grade-language locks: never over-claim; the L7 result is
   "in the model's own simulations, no fixed detector-response rule of
   the standard kind can reproduce these spin-dependent arrival
   patterns" — with the caveat that a real detector's back-action is
   exactly what the experiment tests.
3. Teachers' Edition: mirror the Primer insertions with instructor
   notes: what to emphasize, the honest-limits list (testbed-grade,
   configuration-specific, FAIL-AT-d25 far-field lesson), and pointers
   to the runnable artifacts as classroom material.
4. Mandatory diff verification vs both corpus2 originals.`, { label: 'fold:primer-teachers-v3', phase: 'Fold', schema: FOLD_SCHEMA }),

  () => agent(`${COMMON}
DOCUMENTS (four small ones):
1. /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/08-GUM-Session-Map-v3.0-ext.md
   (copied from corpus2 v2.0-ext; original at corpus2/08-GUM-Session-Map-v2.0-ext.md)
2. /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/The-Watch-Mode-Transition-Memo-v3-ext.md
   (copied from corpus2 v2-ext; original at corpus2/The-Watch-Mode-Transition-Memo-v2-ext.md)
3. /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/09-External-Replication-Record-v2.0.md
   (copied from corpus2 v1.0; original at corpus2/09-External-Replication-Record-v1.0.md)
4. NEW FILE to write from scratch:
   /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/README.manifest-v5.0-ext.md
   — based on reading corpus2/README.manifest-v4.0-ext.md (do NOT edit
   that file).

These are short documents — read each fully first.
1. Session Map v3.0-ext: header update + ⟦T6⟧ entries extending the map
   with the Tier-6 session (the L1-L7 workstreams, the two workflows,
   L5-prime, their artifact paths, findings F-T6-1…6 + F-T6-L5 +
   F-T6-L7a/b). Match the map's existing entry format.
2. Watch-Mode memo v3-ext: header update + a dated ⟦T6⟧ status
   paragraph: second external campaign phase complete before any
   adjudicator reported; supplies (i) the canonical quantum-sector
   exhibits, (ii) the Lock-2 lemma, (iii) the A3 testbed + POVM
   exclusion at testbed grade — still NO nature-facing evidence,
   triggers NO branch; gate-register consequences: the Layer-2 funding
   gate's POVM-exclusion requirement now has a constructive candidate
   [DW, testbed-grade, configuration-specific] with the continuum proof
   as the named residual; the (k0, d_far) joint design constraint filed
   for the Layer-1 envelope. [P5, P6, P7]
3. External Replication Record v2.0: header update + a ⟦T6⟧ Tier-6
   section appended in the record's own format summarizing the campaign
   phase: gate counts, headline numbers (from the payload only), the
   three coordinator-owned spec defects, artifact + audit-log paths.
4. README.manifest v5.0-ext (new): mirror the v4.0-ext manifest's
   structure; list every corpus3 file with its edition lineage
   (corpus2 source -> corpus3 edition), the charter + payload files,
   the carried-forward-by-reference rule for unlisted corpus2 files
   (S-30), and a short Tier-6 fold summary with finding IDs and
   artifact pointers. Begin it with the standard epistemic notice.
5. Mandatory diff verification for the three edited files vs their
   corpus2 originals (the manifest is new — no diff needed, but state
   that).`, { label: 'fold:small-docs', phase: 'Fold', schema: FOLD_SCHEMA }),
])

return { paper, course, primer, small }