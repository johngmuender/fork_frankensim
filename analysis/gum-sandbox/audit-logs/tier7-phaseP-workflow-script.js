export const meta = {
  name: 'tier7-phaseP-execute',
  description: 'Execute Phase P: corpus3 v4.2 delta (N/O fold) + T4-W5 certification pilot',
  phases: [
    { title: 'Execute', detail: 'Two independent workstreams: document fold and exact-field precision certification' },
  ],
}

const RESULT_SCHEMA = {
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
          verdict: { type: 'string', enum: ['PASS', 'FAIL', 'PARTIAL', 'CONDITIONAL'] },
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
You are executing one workstream of the GUM replication campaign's
Tier 7 Phase P. Campaign rules, binding:
- python3 with numpy/scipy/matplotlib is installed (mpmath: pip install
  if needed).
- The pre-registered spec is the Phase-P addendum of
  /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/ROADMAP_v9_PROGRAM.md
  — read it first; it is authoritative on gate wording.
- Honesty over success; resolution refinements allowed (~2); failures
  reported with diagnosis.
- Deliverables per spec; RESULTS.md (or verification record) in campaign
  style with "Within-model; nothing here bears on nature.", gates table,
  key numbers, method, caveats. Units hbar = m = 1.
- Your final structured output must faithfully mirror your deliverable.
`

phase('Execute')

const [p1, p2] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM P1 — corpus3 v4.2 delta: the N/O fold (document editor).
Read FIRST, binding:
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/REVISION_CHARTER_v2.md
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/T7_FOLD_PAYLOAD.md
  — specifically the "ADDENDUM — Phase N/O payload" (items Q6-Q9), the
  single source of truth for every number; the marker stays ⟦T7⟧.
Also skim tier7-program/N1's prior fold pattern by inspecting the
current paper (grep for the v4.1 REVISION NOTE).

Tasks:
1. bash git mv /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/01-GUM-Omega-Paper-v4.1-ext.md 01-GUM-Omega-Paper-v4.2-ext.md; bump ALL six version-header lines v4.1-ext -> v4.2-ext (keep the "(Tier-6/7 fold)" parenthetical; version-header exception; everything else additive-only).
2. Insert compact ⟦T7⟧ blocks (corpus prose style, Q-payload-sourced,
   artifact-cited):
   - Head matter: a ⟦T7⟧ REVISION NOTE (v4.2-ext) directly after the
     v4.1 note: adds the Phase-N/O fold (F-T7-N2/N3/O1/O2/O3); charter
     + payload-addendum pointers; adoption is the authors'.
   - Sec. 0: dated Phase-N/O record [Q9].
   - Sec. VII.I or IX.E area (find the Majoron dictionary / P-nu4
     text): the Q6 block — honest second-moment rate, ABSENT-CERTIFIED
     archaeology, the armed conditional kill, the F-R9 terminal menu,
     S1/S6 untouched.
   - Sec. VIII.F area (after the existing T6 proposal #5 box): the Q7
     block — kill-bin zeros at field/flow grade (O1 numbers), the N3
     line-criterion refutation, T4-W5 residual narrowed to the
     certified/analytic form.
   - Sec. VI.C dispute box or App. E.6 area (find the F-R14 text): the
     Q8 block — the computed decision matrix, ONE sound cell
     (P-acoustic + supertrace), totals, tuning fraction 1/4, minimum
     bill computed; decision support only.
   - Sec. X: Phase-N/O campaign record [Q9].
   - Sec. XI: register update — T4-W5 residual narrowed (line route
     closed by proof; field/flow grade attained; certified/analytic
     form the terminal open); F-R9 terminal menu + the armed external
     kill; F-R14 carries the matrix pointer.
3. Update README.manifest-v5.0-ext.md (census cell v4.1 -> v4.2 as the
   noted header-class exception + a version-history line) and
   09-External-Replication-Record-v2.0.md (append ⟦T7⟧ §9 for phases
   N/O in the record's format).
4. VERIFICATION (mandatory): diff the renamed paper vs git HEAD's
   v4.1-ext content (git show HEAD:analysis/gum-sandbox/corpus3/01-GUM-Omega-Paper-v4.1-ext.md);
   confirm deletions = version headers only; count ⟦T7⟧ markers before/
   after; confirm ⟦rev⟧/⟦T6⟧ and the eight prior ⟦T7⟧ blocks preserved
   verbatim; spot-check >= 15 Q6-Q9 payload numbers verbatim.
Gates: G1 additive-only, G2 marker preservation + counts, G3 payload
fidelity, G4 cross-file consistency. File as the v4.2 fold record.`, { label: 'P1:v4.2-fold', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM P2 — T4-W5 certification pilot (exact-field, high-precision backward paths).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/P2/ (create).
MANDATORY context: read tier7-program/O1/RESULTS.md + o1_run.py (the
backward-reachability run you are upgrading; its classification data is
in O1/o1_raw.npz — inspect its arrays), and the engine in
tier6-foundations/L5prime/l5prime_run.py.

KEY IDEA: the engine's field is a FINITE trigonometric sum — analytic
sine-mode factors in z (the initial state is the pure ground mode, so
z-structure is exactly known) and a band-limited Fourier representation
in x (the FFT grid coefficients ARE the exact spectral coefficients of
the band-limited field). Therefore psi (hence rho, j, v) can be
evaluated EXACTLY (to floating-point rounding) at ARBITRARY (x, z, t)
by direct mode summation: x-part = sum over the Nx Fourier modes with
their exact e^{i k x} factors and exact free-evolution phases; z-part
analytic. This removes O1's dominant error source (grid interpolation)
entirely. NOTE the L5prime engine is exact-in-z factorized — recover
its representation faithfully; verify your direct evaluator against
the grid engine AT GRID NODES (should agree to ~1e-13 rel).

Runs (n_y in {+1, +0.75} — mirrors are numerically identical per O1):
- Stratified path subset per the roadmap: ALL 16 box-exit/wall-adjacent
  paths + the 100 smallest-margin class-(a) paths + 100 random
  class-(a) paths per n_y (>= 432 total; read O1's raw npz to select).
- Integrate backward with scipy DOP853 (rtol 1e-12, atol 1e-14) on the
  exact-field velocity (v = a + n_y b computed from direct mode sums).
- Precision ladder: for >= 20 paths (the worst margins + a spread),
  re-integrate with mpmath at 50 digits (fixed-step high-order RK or
  mpmath odefun; state the method) and compare endpoints/margins.

Gates (roadmap authoritative):
- G1: exact-vs-grid field agreement at nodes <= 1e-12 rel (sampled over
  >= 1000 random nodes x several times).
- G2: every re-run path keeps its O1 classification; margins agree with
  O1 within O1's own error estimate (p95 2.8e-5 class).
- G3: precision ladder — float64-DOP853 vs 50-digit endpoints <= 1e-8;
  margins stable to >= 6 digits on the ladder subset.
- G4: the worst-margin path certified at ladder grade (print its full
  record: start point, margin, crossing epoch, ladder agreement).
- G5: honest scope printed — precision-certification, NOT formal
  interval arithmetic; name what formal certification still requires
  (validated ODE enclosures over the exact field) and that the analytic
  proof remains the terminal open.
Deliverables: P2/p2_certify.py, p2_results.json, one figure (margin
distribution + ladder agreement), RESULTS.md. File as F-T7-P2.`, { label: 'P2:certification-pilot', phase: 'Execute', schema: RESULT_SCHEMA }),
])

return { p1, p2 }