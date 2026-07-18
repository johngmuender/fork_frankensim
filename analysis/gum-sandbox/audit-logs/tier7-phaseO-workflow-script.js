export const meta = {
  name: 'tier7-phaseO-execute',
  description: 'Execute Phase O: flow-map first-crossing proof, g_lim archaeology, F-R14 rescue decision matrix',
  phases: [
    { title: 'Execute', detail: 'Three independent workstreams closing the named opens' },
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
          verdict: { type: 'string', enum: ['PASS', 'FAIL', 'PARTIAL', 'CONDITIONAL', 'VERDICT'] },
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
Tier 7 Phase O. Campaign rules, binding:
- python3 with numpy/scipy/matplotlib is installed.
- The pre-registered spec is the Phase-O addendum of
  /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/ROADMAP_v9_PROGRAM.md
  — read it first; it is authoritative on gate wording.
- Honesty over success: no tuning; resolution refinements allowed (~2);
  physics/verdict failures reported honestly with diagnosis.
- Deliverables: script(s) where computational, <name>_results.json, one
  PNG figure (computational workstreams), RESULTS.md (or the named memo)
  in campaign style with the "Within-model; nothing here bears on
  nature." notice, gates table, key numbers, method, caveats.
- Units hbar = m = 1 unless stated.
- Your final structured output must faithfully mirror your deliverable.
`

phase('Execute')

const [o1, o2, o3] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM O1 — the flow-map first-crossing proof (T4-W5 terminal shape).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/O1/ (create).
MANDATORY context: read tier7-program/N3/RESULTS.md + n3_field.py (the
refuted line criterion + the engine reuse pattern), tier6-foundations/
L7b/RESULTS.md + l7b_run.py, and L5prime/l5prime_run.py. REUSE the
validated engine (exact-in-z factorized evolution, C2 smoothstep,
k0 = 2, waveguide z in [0,1], box x in [-15,45], Nx = 2048 reference,
dt <= 5e-4) and the exact velocity decomposition
v(x,z,t; n_y) = a(x,z,t) + n_y b(x,z,t).

Theorem to state + prove (short): guidance trajectories are unique
forward AND backward (the velocity field is smooth wherever rho > 0);
therefore if the backward trajectory from a detector point (x = d, z,
t) — integrated with the SAME field — satisfies max_{s < t} X_x(s) >= d,
the unique trajectory through that point already crossed the line
earlier, so it contributes NO first crossing at time t. If EVERY
kill-bin detector point (up to a quantified rho-negligible set) is of
this type, the first-crossing (arrival) density in the kill bins is
zero at field level.

Execute for n_y in {+1, -1, +0.75, -0.75}, d = 1, kill bins
t in [6.4, 6.8] and [6.8, 7.2]:
- Backward-integrate (RK4, the stored/reconstructed velocity fields;
  integrate from t down to 0) from a grid >= 200 z-points x 40 t-points
  per bin per n_y.
- Classify each start point: (a) PRE-CROSSED: the backward path has
  max X_x >= d + margin at some s < t (record the margin and the
  crossing epoch); (b) RHO-FLOOR: the backward path enters regions with
  rho < rho_floor (choose rho_floor ~ 1e-10 of peak; these are
  near-node/near-wall zones where the field is ill-conditioned) —
  count them and bound the |phi_0|^2-measure they can carry;
  (c) VIOLATION: the backward path reaches t = 0 inside supp(phi_0)
  with max X_x < d throughout — a genuine fresh kill-bin arrival, which
  would REFUTE the zeros. Also record paths that leave the box.
- Refinement: repeat classification at doubled Nx and halved backward
  dt for at least the boundary-adjacent subset (or all, if cheap);
  classifications must be stable.
- Discrimination control: n_y = +0.25, its POPULATED bins near t ~ 6.5:
  backward runs must show a substantial legitimately-fresh fraction
  (paths reaching supp(phi_0) without pre-crossing) — the mechanism
  that populates those bins for the no-cutoff member.

Gates (roadmap authoritative): G1 engine bars (1e-6 class or better);
G2 zero VIOLATIONS at reference + stable under refinement;
G3 pre-crossing margin: min over class-(a) of (max X_x - d) quantified
and >> the backward-trajectory discretization error (measure via
refinement differences of the paths themselves); G4 rho-floor
accounting <= 1e-3 of |phi_0|^2-measure (else PARTIAL with the bound);
G5 discrimination at n_y = +0.25. Verdict target: kill-bin zeros
elevated to FIELD/FLOW statement, proof-grade modulo discretization;
print the honest scope (analytic proof still open; its demonstrated
shape is this backward-reachability argument). File as F-T7-O1.`, { label: 'O1:flow-map-proof', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM O2 — g_lim archaeology (the F-T7-N2 blocking number).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/O2/ (create).
Context: read tier7-program/N2/RESULTS.md first (the conditional kill:
any Majoron-mode 0nubb rate-level bound g_lim < 8e-3 would exclude the
honest second-moment rate at every anchor and every printed g; the
archive appeared to print NO such bound — your job is to certify or
refute that at archaeology grade).

Sweep the ENTIRE archive exhaustively:
- Directories: /home/user/fork_frankensim/analysis/gum-sandbox/corpus/
  (the 104-file archive), corpus2/, corpus3/, theory-audit/,
  DISCHARGE_PACKAGE/, tier0-gauntlet ... tier7-program, gum-core/,
  and the top-level *.md files of analysis/gum-sandbox/.
- Use >= 2 independent pattern families (grep -ri):
  family A (mode/observable): '0nubb|0\\u03bd\\u03b2\\u03b2|neutrinoless|double.beta|T_?1/2|half.life|Majoron.mode|P-nu4|P-\\u03bd4'
  family B (bounds/limits): 'g_?lim|g <|g<|g \\u2264|limit on g|bound on g|KamLAND|EXO|GERDA|LEGEND|nEXO|NEMO|CUORE|detectab'
  plus a numeric family C: '1e-9|10\\^-9|10\\u207b\\u2079' near Majoron
  contexts. Iterate patterns as needed; log every pattern actually run.
- For EVERY hit: record file path, line, verbatim quote, and classify:
  (1) Majoron-mode RATE-level bound (g_lim or T_1/2-equivalent or a
  quantitative detectability line) — the blocking number; (2) a
  non-rate g-constraint (SN cooling, free-streaming, BBN, recoupling —
  constrains g itself, not the mode-rate observable); (3) stake/kill
  language about 0nubb OCCURRENCE (S6 — occurrence is not a rate
  number); (4) other/incidental.
- Decision rule (pre-registered): if class-(1) found -> execute the
  F-T7-N2 conditional-kill arithmetic against it (does the honest rate
  R x printed exceed it? at which anchors?) and report which way it
  fires; if only class-(2)/(3) -> catalog them + state precisely why
  they cannot adjudicate the honest-rate question; if nothing -> certify
  ABSENT with the full search protocol printed (reproducible).

Gates: G1 sweep coverage — file counts per directory + the full pattern
log printed; G2 the catalog — every hit quoted + classified (expect the
corpus's own battery items: B-nu1 free-streaming f >= 1.4 MeV, SN
band, BBN, the g-window [8e-10, 1.3e-8], P-nu4 'far-horizon' language,
S6/LEGEND-1000/nEXO stake rows — classify each); G3 the verdict per
decision rule. Deliverables: O2/GLIM_ARCHAEOLOGY.md + o2_results.json
(the machine-readable catalog) + the grep driver script o2_sweep.sh or
.py. No figure needed. File as F-T7-O2.`, { label: 'O2:glim-archaeology', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM O3 — F-R14 rescue decision matrix (decision support).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/O3/ (create).
MANDATORY context: read theory-audit/h33_RESULTS.md (Part B and the
rescue-pricing section 8; also h33_pq.py if present) and
theory-audit/h25_RESULTS.md (the K3 supertrace row), plus the F-R14
dispute box and App. E.6 rewrite in
corpus3/01-GUM-Omega-Paper-v4.1-ext.md (grep 'F-R14' and 'E.6').

Goal: turn the corpus's one unrepaired structural defect into a fully
computed decision matrix — every textually available reading priced,
nothing recommended (the constitutive choice is the authors').

Implement:
- The master formula a1(p,q) = -(5p+q)/(12(2p+q)) (h33's derived form)
  where (p,q) encode the sector's index-variance/density-weight
  structure; the B2 photon doublet has the structural invariant
  p = q - 1; Qtilde-sectors have q = 3/2 group-protected; the knot
  band-edge (Dirac-type) sector contributes a1 = -1/3 per field BEFORE
  the loop sign, with the supertrace-signed count per h25 K3 (show the
  row both ways: naive count vs supertrace-signed).
- Readings to evaluate (rows of the matrix): frame/soldered (favored),
  covector, vector, cone-only — extract each reading's (p,q) per
  sector from h33's anchor values and VERIFY: your computed a1^(B2)
  must reproduce h33's printed anchors exactly: -1/3 (covector),
  -2/15 (frame/soldered), -5/33 (vector), -1/12 (cone-only).
- The P-acoustic family: effective density weight w in (-5/18, -2/9)
  (width 1/18), distinguished point w = -1/4 on the q = -3p locus with
  a1 = +1/6 — reproduce that anchor too; scan w across the window
  (>= 50 points): per-sector a1(w), weights w_s = N_s Lambda_s^2 a1,
  sign of Sigma w_s (1/16 pi G), convexity/hull status for c_GW^2
  (with a representative sector content: B2 doublet + B3 + B4 +
  knot band-edges; use the corpus's N_s multiplicities where printed,
  else parametrize and say so).
- Matrix columns per reading/window-point: G-sign, hull status
  (convex combination? can c_GW^2 exit?), (6.3) margin status
  (the '4.0 orders boundary-exact' arithmetic under that reading),
  and the priced costs (for P-acoustic: the three printed costs +
  the tuning fraction = window width / distance to nearest natural
  convention w = 0 or -1/2).

Gates: G1 reproduce h33's five anchors EXACTLY (the four readings +
the P-acoustic +1/6); G2 the decision matrix fully computed (every
cell a number or a computed status, none asserted); G3 P-acoustic
pricing: tuning fraction quantified; supertrace erratum shown both
ways on the knot row (does the G-sign conclusion flip with it in any
reading? report); G4 the honest bottom line printed: decision support
only; adoption is the authors'; no reading recommended; the F-R14
regrade language unchanged. Deliverables: O3/o3_matrix.py,
o3_results.json, one figure (weights vs w across the P-acoustic
window), RESULTS.md with the matrix as its centerpiece table. File as
F-T7-O3.`, { label: 'O3:FR14-matrix', phase: 'Execute', schema: RESULT_SCHEMA }),
])

return { o1, o2, o3 }