export const meta = {
  name: 'tier7-phaseN-execute',
  description: 'Execute Phase N: corpus3 v4.1 fold (F-T7), F-R9 route (a) honest second moment, L7 continuum proof',
  phases: [
    { title: 'Execute', detail: 'Three independent workstreams: document fold, second-moment propagation, field-level crossing criterion' },
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
Tier 7 Phase N. Campaign rules, binding:
- python3 with numpy/scipy/matplotlib is installed.
- The pre-registered spec is the Phase-N addendum of
  /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/ROADMAP_v9_PROGRAM.md
  — read it first; it is authoritative on gate wording.
- Honesty over success: no tuning to pass; resolution refinements
  allowed (~2); physics failures = FAIL with diagnosis.
- Numerical workstreams deliver: script(s), <name>_results.json, one
  PNG figure, RESULTS.md in campaign style (goal header; "Within-model;
  nothing here bears on nature."; gates table; key numbers; method;
  caveats). Units hbar = m = 1 unless stated.
- Your final structured output must faithfully mirror your RESULTS.md
  (or, for the document workstream, your verification record).
`

phase('Execute')

const [n1, n2, n3] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM N1 — corpus3 v4.1 delta: the F-T7 fold (document editor).
Read FIRST, binding:
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/REVISION_CHARTER_v2.md
- /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/T7_FOLD_PAYLOAD.md (single source of truth for every number; marker for this pass is ⟦T7⟧, distinct from ⟦rev⟧ and ⟦T6⟧ which stay verbatim)

Tasks:
1. git mv /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/01-GUM-Omega-Paper-v4.0-ext.md to 01-GUM-Omega-Paper-v4.1-ext.md (use bash git mv). Bump ALL version-header lines (main title + the five PART-N running heads) v4.0-ext -> v4.1-ext (version-header exception; everything else additive-only).
2. Insert compact ⟦T7⟧ blocks (corpus prose style, payload-sourced, artifact-cited):
   - Head matter: a ⟦T7⟧ REVISION NOTE (v4.1-ext) directly after the ⟦T6⟧ revision note: this edition adds the Tier-7 program-extension fold (F-T7-M1..M4 + F-T7-1); charter + payload pointers; adoption is the authors'.
   - Sec. 0: dated third-phase paragraph [payload Q5].
   - Sec. III.D area (after the existing ⟦rev⟧/⟦T6⟧ annotations): M4 — Nelson relaxation 17-41x faster, dBB log-singular, rates = direction-signed conservative bounds, double edge [Q4].
   - Sec. IV.H.3 area: M2 — servo law-form replicated; c = eta/2 bath functional; corpus 0.021 <-> eta = 0.042; constant-c g^2 law pins the r2 bath to the Ohmic class [Q2].
   - Sec. VI.D area: M1 — tracker identity dynamical (1.8e-12); Gamma ~ H^n => W = 1-n; memory-pair convention RECOVERED (W(z=0)->W(z=2)); EdS pathology measured; PLUS the F-T7-1 minor flag on the "(+0.36 at beta=-0.3)" print [Q1].
   - Sec. X: T7 campaign record incl. the two coordinator-owned spec defects [Q5].
   - Sec. XI: register update — WS4-M4c (open #4) partially discharged at toy grade (Gamma(H)->W law + convention recovered; residual: field-level Gamma and the P1 premise); F-R9 repair menu narrowed by M3 (route (b) closed; routes (a)/(c) remain) [Q1, Q3].
3. Update /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/README.manifest-v5.0-ext.md: version-history line for v4.1 + the paper's new filename in the census (additive edits + the filename cell correction, which you may treat as a version-header-class exception, noting it).
4. Append a ⟦T7⟧ section to /home/user/fork_frankensim/analysis/gum-sandbox/corpus3/09-External-Replication-Record-v2.0.md (a §8 Tier-7 phase record in the record's own format; header untouched).
5. Append a Tier-6/7 section to /home/user/fork_frankensim/analysis/gum-sandbox/REPLICATION_CAMPAIGN_STATUS.md (additive; matches that file's list style; summarizes tiers 6-7 with finding IDs and artifact paths; do not modify existing lines).
6. VERIFICATION (mandatory, report): diff the renamed paper against git HEAD's v4.0-ext content (git show HEAD:analysis/gum-sandbox/corpus3/01-GUM-Omega-Paper-v4.0-ext.md > /tmp/v40.md then diff) — confirm every deletion is a version-header line and every insertion is ⟦T7⟧-marked; count ⟦T7⟧ markers in every touched file; confirm all ⟦rev⟧/⟦T6⟧ markers preserved; spot-check >= 10 payload numbers appear verbatim.
No PNG required for this workstream; your structured output's gates array should encode the verification checks as gates (G1 additive-only, G2 marker counts, G3 payload fidelity, G4 cross-file consistency).`, { label: 'N1:v4.1-fold', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM N2 — F-R9 route (a): the honest second-moment treatment.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/N2/ (create).
MANDATORY: read theory-audit/h21_RESULTS.md + h21_mc.py (reuse the
parametrization verbatim) AND tier7-program/M3/RESULTS.md (the route-b
closure — your sibling result) AND theory-audit's J3 references in
REPLICATION_CAMPAIGN_STATUS.md (the Majoron battery re-run context;
also read theory-audit/j3_RESULTS.md if present).

Goal: F-R9's route (a) — instead of rescuing the corpus's <A>^2-based
mean-rate claim, compute what an HONEST second-moment treatment gives
at the observable level, and adjudicate the consequence for the
corpus's 0nubb-insert (K3a) phenomenology.

Physics to establish and compute:
1. Spatial self-averaging: for a macroscopic sample (N_nuclei ~ 1e27
   over cm^3-km^3 scales vs the web's um-scale correlation), the
   sample-total rate = N * <A^2> (the spatial mean of the true local
   quadratic rate). The web is static on experiment timescales (cite
   h21's motional-narrowing pricing: narrowing would need 5e7-1.2e9 c).
   Hence the honest total rate is R x the corpus's claimed rate with
   R = <A^2>/<A>^2 = 1e14-1e18 across the corpus scale window
   (recompute R at the 4 h21/M3 anchors — G1 reproduction).
2. Burstiness/heterogeneity: measure the rate Lorenz structure — the
   fraction of total rate carried by the top 10^-k fraction of sites
   (k = 2..12) at >= 2 anchors: importance-sampled MC on h21's
   parametrization. Report the effective 'hot-site' picture (rate
   concentrated on near-tube shells).
3. Propagation to the corpus's phenomenology: the Majoron-mode 0nubb
   rate scales as g^2 (state the corpus's P-nu4 g ~ 1e-9 target and the
   g-window [8e-10, 1.3e-8]); honest rate = R x printed rate means the
   effective g_equiv = sqrt(R) x g = 1e7-1e9 x g. Build the propagation
   table: per anchor -> R -> g_equiv shift -> where the printed window
   lands relative to (a) the corpus's own far-horizon P-nu4
   detectability line and (b) the J3 battery's 'SHIFTS-HARMLESSLY'
   verdict at g = 7.8e-10. Pre-registered decision frame (no tuning):
   conclude ONE of: (i) the battery/insert survives with re-priced
   parameters (state the required retreat), (ii) the honest rate
   overshoots a corpus-printed bound => a NEW candidate defect/kill
   sharpening filed at F-T7-N2, or (iii) indeterminate-from-archive
   (name exactly which unprinted number blocks the verdict — the
   L3-style isolated residue). IMPORTANT honesty note: the 0nubb
   MODE rate claim (K3a insert) is what F-R9 wounds; S1 (Sigma m_nu)
   is untouched either way — print this.
4. Also state what route (a) does NOT settle: it is an honest
   recomputation, not a repair — the corpus must either adopt the
   R-enhanced rate (with its phenomenological consequences) or retreat
   to route (c) at the line-support cost.

Gates: G1 h21 anchor reproduction (R within x3 at >= 2 anchors; mean
1.000-class); G2 burstiness measured (Lorenz fractions with MC errors);
G3 the propagation table + ONE pre-registered verdict with arithmetic;
G4 honest-limits paragraph (upstream K3a physics; S1 untouched; what
remains corpus-side). File as F-T7-N2.`, { label: 'N2:second-moment', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM N3 — L7 continuum proof of the cutoff zeros (T4-W5 residual).
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/N3/ (create).
MANDATORY: read tier6-foundations/L5prime/RESULTS.md + l5prime_run.py
and tier6-foundations/L7b/RESULTS.md + l7b_run.py — REUSE the validated
engine (exact-in-z factorized evolution, C2 smoothstep, k0 = 2,
waveguide z in [0,1], box x in [-15,45], Nx = 2048 reference,
dt <= 5e-4) and the exact L7b velocity decomposition
v(x,z,t; n_y) = a(x,z,t) + n_y b(x,z,t).

Goal: elevate the sampled hard-cutoff zeros (tau_max = 5.130950 at
n_y = +1; two empty kill bins) from ensemble statements to a FIELD
statement via the crossing criterion.

Theorem (state + prove in RESULTS.md — it is two lines): if a
trajectory crosses the line x = d at time t moving rightward (a
first crossing from x < d requires X(t) = d with dX/dt >= 0 there),
then v_x(d, Z(t), t) >= 0. Contrapositive: if
sup_z v_x(d, z, t) < 0 for ALL t in (tau*, T], then NO trajectory
crosses x = d in (tau*, T] — the continuum first-crossing density
vanishes there; support subset [0, tau*]. (Also note the density-
weighting subtlety honestly: v_x is defined wherever rho > 0; handle
rho -> 0 regions by reporting the criterion on the region
rho >= rho_floor AND showing the complement carries negligible
not-yet-crossed measure — or better, show v_x < 0 holds on the whole
line including low-rho zones within the field's resolution.)

Compute (n_y = +1, d = d_near = 1, T = 16):
- v_x(d, z, t) on the full (z, t) grid (fine t sampling, dt_out <=
  0.01); M(t) = max_z v_x(d, z, t); tau* = sup{t : M(t) >= 0}.
- Margin delta(t) = -M(t) on (tau*, T]; delta_min; compare to a
  discretization error bound for v_x on the line (estimate via
  Nx-doubling and dt-halving differences of v_x itself).
- Consistency with the sampled tau_max = 5.130950 (and the kill bins
  [6.4,6.8], [6.8,7.2] at |n_y| >= 0.75: also compute tau*(n_y) for
  n_y in {0.75, -0.75, -1} to cover the kill-bin set).
- Discrimination control at n_y = +0.25: M(t) must recur >= 0 to late
  times (no tau* below T).

Gates: G1 field quality (norm/walls/energy at the L5prime bars);
G2 tau* exists for n_y = +1 and is stable <= 1% under Nx-doubling AND
dt-halving; G3 consistency: tau* >= 5.130950 with the gap reported and
explained (sampled max <= continuum edge; if tau* < sampled tau_max
that is a CONTRADICTION - FAIL with diagnosis); G4 margin:
delta_min on (tau* + 0.1, T] exceeds the v_x discretization bound by
>= 10x (quantify both) — the theorem's hypothesis verified at
proof-grade-modulo-discretization, upgrading the kill-bin zeros to a
field statement for every |n_y| >= 0.75 checked; G5 discrimination:
n_y = +0.25 shows no tau* below T. Honest scope paragraph: proves the
continuum support statement GIVEN the computed field (bounds printed);
the fully analytic proof remains the named open. File as F-T7-N3.`, { label: 'N3:continuum-proof', phase: 'Execute', schema: RESULT_SCHEMA }),
])

return { n1, n2, n3 }