export const meta = {
  name: 'tier7-phaseQ-execute',
  description: 'Execute Phase Q: validated interval enclosures for pre-crossing certificates + analytic reconnaissance of the terminal T4-W5 proof',
  phases: [
    { title: 'Execute', detail: 'Two independent workstreams: formal-rung interval certification and the analytic memo' },
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
Tier 7 Phase Q. Campaign rules, binding:
- python3 with numpy/scipy/matplotlib/mpmath (incl. mpmath.iv) is
  installed and verified working.
- The pre-registered spec is the Phase-Q addendum of
  /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/ROADMAP_v9_PROGRAM.md
  — read it first; authoritative on gate wording.
- Honesty over success; no tuning; failures reported with diagnosis;
  long runs MUST be launched with nohup in the background so they
  survive turn boundaries and container restarts, with per-phase or
  per-path checkpoint files so work is resumable.
- Deliverables: script(s)/memo per spec, <name>_results.json where
  computational, RESULTS.md (or the named memo) in campaign style with
  "Within-model; nothing here bears on nature.", gates table, key
  numbers, method, caveats. Units hbar = m = 1.
- Your final structured output must faithfully mirror your deliverable.
`

phase('Execute')

const [q1, q2] = await parallel([
  () => agent(`${COMMON}
WORKSTREAM Q1 — validated interval enclosures for the pre-crossing certificates.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/Q1/ (create).
MANDATORY context: read tier7-program/P2/RESULTS.md + p2_certify.py
(the exact band-limited field representation and the path data; the
selection/margins live in P2/_selection.npz, _dop.npz, _ladder.json,
p2_results.json — inspect them; the worst-margin path record is in
p2_results.json G4: selection_index 232, n_y = 0.75, start t = 6.405,
z = 0.6825, margin ~0.6207, crossing epoch s_pc ~ 2.85). Also read
tier7-program/O1/RESULTS.md for the theorem being certified.

CRITICAL SCOPE INSIGHT (from the roadmap): the certificate needs ONLY
the backward segment from the detector point (d=1, z, t) until X_x
exceeds d + delta at some earlier s — for most paths that happens
within ~0.5-1.5 time units backward (s_pc median 5.7/4.9 vs t in
[6.4,7.2]); WAIT — note s_pc for the worst path is 2.85, i.e. ~3.5
time units backward. Certify the segment [s_stop, t] where s_stop is
where the enclosure's X_x lower bound first exceeds d + delta — you do
NOT need to reach s_pc if the tube crosses d earlier... careful: the
backward path X_x DIPS below d first (it is in the not-yet-... no: it
was AT d at time t and pre-crossed earlier, meaning backward-in-time
the path goes left then comes back above d at s_pc. You must integrate
backward until the tube's X_x lower bound > d. Check the float path
first (from P2 data or by quick float integration) to know each
path's actual excursion structure and expected certification depth;
choose candidate paths with SHORT certification depth for the >= 5
extra paths (many class-(a) paths cross d + margin quickly), and
attempt the worst-margin path with its ~3.5-unit depth as the
stretch case (report honestly if the tube wraps before certifying —
that is the G4 failure-accounting gate, not a scandal).

Implementation:
- Interval field: the exact representation from p2_certify.py (finite
  mode sum; analytic z-factor; free-evolution phases). Reimplement its
  evaluation in mpmath.iv: interval x, z, t -> interval psi (real +
  imag), rho, j, v = a + n_y b. Spectral coefficients: float64 values
  declared exact-model data (rationals). Precision: iv.mp.prec ~ 80-120
  bits; keep the mode count manageable — if the full Nx = 2048 mode sum
  is too slow in interval arithmetic, you may certify on a truncated
  spectral model ONLY IF you add a rigorous truncation error bound
  (sum of |dropped coefficients| bounds the field error; propagate it
  as an interval inflation) — state the choice.
- Validated integrator: stepwise a priori box (Picard: B = X + [0,h] *
  F(B_trial), verify containment F(B) * h inside), then a first-order
  (or midpoint) interval step with rigorous remainder via an interval
  bound on the derivative over B (finite differences are NOT rigorous —
  bound dF via interval evaluation of F over the box and a Lipschitz
  estimate from interval arithmetic on the analytic derivative if
  implementable, else use the crude but rigorous |F(B)| * h + width
  inflation with the second-order term bounded by interval evaluation
  of F at box corners... choose a defensible construction and PRINT it;
  the construction must be genuinely rigorous given the trust
  statement, not heuristic).
- Stop when lower(X_x) > d + delta (delta = 1e-3 target); record s_stop,
  the certified lower bound, and max interval width reached.

Gates (roadmap authoritative): G1 containment checks (interval field
encloses >= 1e4 random float evaluations, 0 violations; outward
rounding sanity on >= 100 cases); G2 the worst-margin path certified
(or honestly reported as wrap-failed under G4 accounting — but then
certify instead the deepest path that succeeds and report the achieved
depth); G3 >= 5 additional paths certified across both kill bins and
both |n_y| families incl. >= 1 small-margin path; G4 failure
accounting (any wraps/blow-ups with widths); G5 trust + scope
statement (mpmath.iv correct rounding assumed; floats-as-exact-model;
CAPD-grade = the remaining formal step; analytic proof = terminal).
Deliverables: Q1/q1_enclose.py, q1_results.json, one figure (tube
widths + certified bounds), RESULTS.md. File as F-T7-Q1.`, { label: 'Q1:interval-enclosures', phase: 'Execute', schema: RESULT_SCHEMA }),

  () => agent(`${COMMON}
WORKSTREAM Q2 — analytic reconnaissance of the terminal T4-W5 proof.
Directory: /home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/Q2/ (create).
MANDATORY context: read tier7-program/N3/RESULTS.md (the refuted line
criterion + the transport mechanism diagnosis), O1/RESULTS.md (the
backward-reachability theorem + measured crossing epochs), and the
field construction in tier6-foundations/L5prime/l5prime_run.py (initial
data: C2-smoothstep-truncated Gaussian x-profile, boost k0 = 2, ground
sine mode in z; free evolution; hard walls z in [0,1]).

Deliverable: Q2/ANALYTIC_RECON.md + supporting q2_check.py — a scoped
analytic memo, campaign register, with:

(i) A PROVEN LEMMA: the late-time ballistic asymptotic for the
band-limited field class. For free 1-D evolution of a band-limited
packet (finite Fourier sum, compact spectral support [k_min, k_max]
with k > 0 after the boost — check the actual spectrum of the L5prime
initial data first; the C2 truncation makes the spectral support
technically all of k-space but with rapidly decaying tails — handle
honestly: either (a) prove the lemma for a strictly band-limited
proxy and bound the proxy-vs-actual difference numerically, or
(b) prove it with explicit tail-dependent error terms). Target
statement class: for x/t inside the group-velocity range,
v_x(x, t) = x/t + E with |E| <= C(x,t) explicit and -> 0; derive C
exactly via the Fresnel/stationary-phase decomposition of the finite
mode sum (this is exact computable analysis, not asymptotics folklore
— show the derivation). Since the 2-D problem factorizes (z-factor =
ground mode alone for the AXIAL part; the transverse spin term adds
the known b-term), state clearly which part of v the lemma controls.
Verify the lemma numerically at >= 100 (x, t) points against the
engine's field (tier7-program/P2's exact evaluator pattern).

(ii) THE CONJECTURE, precisely: the transport/no-slow-lane mechanism —
for |n_y| >= n*, every not-yet-crossed trajectory is z-advected into
the fast wall lane and crosses by a finite T_max(n_y); write the
needed structure as explicit proof obligations (e.g., a monotone
functional along trajectories; an invariant-region argument for the
(z, x - d) system; what the b-term's sign structure gives near the
walls: b_z ~ +n_y d(rho)/dx / (2 rho)-class — derive the actual
near-wall asymptotics of the spin drift from the ground-mode z-factor
and state the candidate Lyapunov functional).

(iii) WHAT THE LEMMA ALREADY YIELDS: an explicit rigorous statement of
the form "no first crossings at the detector for t > T_ball under
conditions C" — derive T_ball for this configuration's numbers and
state honestly whether the kill bins [6.4, 7.2] fall inside or outside
its reach (likely outside — the lemma controls the far tail, not the
kill window; say so plainly).

Gates: G1 lemma proven + numerically verified (>= 100 points, error
within the derived bound, ZERO violations); G2 conjecture + proof
obligations printed precisely (a competent analyst could pick them up);
G3 the honest delta stated (what the lemma covers vs what the kill-bin
theorem needs). FAIL on any part is permitted and must be diagnosed.
File as F-T7-Q2.`, { label: 'Q2:analytic-recon', phase: 'Execute', schema: RESULT_SCHEMA }),
])

return { q1, q2 }