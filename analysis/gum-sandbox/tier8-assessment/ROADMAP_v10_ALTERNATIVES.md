# ROADMAP v10 — Phase S: Program Assessment + Alternatives & Extensions

**Date:** 2026-07-18.  **Phase:** Tier 8 (assessment + alternatives).
**Branch:** `claude/analyze-gum-po`.  **Coordinator:** session 536b52e3
(Claude Fable, remote).  **Status at pre-registration:** Phases L–R closed
(see `tier7-program/TIER7_ADJUDICATION.md`, `corpus3/` v4.3 edition set);
remaining opens per the Phase-R closeout: (i) formally-verified kernel +
population-scale certification [engineering], (ii) analytic theorem via
O-1…O-5 [mathematics], (iii) F-R14 adoption + archive-external g_lim
[corpus-side].

**Epistemic frame (binding, unchanged).** Within-model; nothing here bears
on nature.  All workstreams test internal consistency, distinctiveness, or
robustness of the corpus's own structures.  Ledger classes [DF]/[DW]/[CAL]/
[IM]/[CJ] as in the charter.  Grade language may never rise as a result of
Phase-S work; honest FAILs are findings.

---

## 0. Phase-S deliverable structure

1. **S-Assess** (workflow `wf_f8099700-bc8`): six parallel sector surveys
   (theory-core, early-tiers, field-family, foundations, program-phases,
   campaign-meta) → pros-advocate + cons-critic briefs → independent
   adversarial evidence audits of both briefs.  Coordinator synthesis into
   two separate documents:
   - `PROS_ASSESSMENT.md` — the strongest audited case FOR the program.
   - `CONS_ASSESSMENT.md` — the strongest audited case AGAINST the program.
   Only SUPPORTED points (or OVERSTATED points at their audited-corrected
   strength) may appear; every point carries file evidence.
2. **S1–S4** (below): four pre-registered executable workstreams, each an
   extension or an *alternative to a GUM component*, each with gates fixed
   before execution.
3. Adjudication (`TIER8_ADJUDICATION.md`), checkpoint addendum, audit-log
   archives, full commit/push (credit-resilience protocol).

---

## S1 — P-acoustic adoption pilot: pricing cost (1) of the F-R14 sound cell

**Component being alternatived:** the flat quantization measure underlying
Sec. III's [DF] derivations (Fisher/Madelung → quantum potential; the Born
machinery), replaced by the P-acoustic density weight w = −1/4
(ρ₀, J ∝ (det g₃)^w-family) — the *only computationally sound cell* of the
O3 decision matrix (`tier7-program/O3/RESULTS.md` §3: P-acoustic +
supertrace, tuning fraction 1/4).

**Question.** O3 carried three printed costs of adoption; cost (1) —
"reopens Sec. III's flat-measure [DF] derivations on curved backgrounds" —
is qualitative in the corpus.  Price it *quantitatively*: does the
Fisher/Madelung derivation of the quantum-potential form survive the
measure change at linear order in the metric perturbation h, and if not,
what is the correction term and its coefficient?

**Method.** Symbolic (sympy), no fitting.  (a) Reproduce the flat-measure
Fisher → quantum-potential identity exactly (the corpus's Sec. III route:
variation of the Fisher information functional I_F[ρ] yields the Bohm
quantum potential with coefficient ħ²/8m, in campaign units).  (b) Repeat
on g₃ = (1+h)δ with measure weight (det g₃)^w, general w, to linear order
in h (constant and gradient h separately).  (c) Evaluate at w = 0 (flat
convention), w = −1/4 (P-acoustic), w = −1/2 (the other natural
convention); identify which terms break covariance of the [DF] chain and
whether any w cancels them.

**Gates (pre-registered).**
- **S1-G1** Flat identity reproduced symbolically, coefficient exact.
- **S1-G2** Curved-measure variation computed at linear order for general
  w; results as explicit symbolic expressions (no truncation beyond O(h²)).
- **S1-G3** The cost is *priced*: either the correction vanishes for some
  printed w (state which, exactly), or the correction term's structure and
  coefficient at w = −1/4 is printed — with the honest statement of what
  this does to the [DF] label on curved backgrounds.
- **S1-G4** Decision-support-only framing: no adoption recommended; F-R14
  regrade language unchanged.

## S2 — Chirality-interpolation family: is the saturated closure constant distinctive?

**Component being alternatived:** the chiral (micropolar, direction-locked)
structure of the substrate energy at closure.  Per
`tier2-closure/gstar_RESULTS.md`, the t → 0 saturated objective is the
BPS-type functional G = ∫[π³b² + W(q)]dV with
W = (1/8π)[(1−cos ψ)² + sin²ψ cos²Θ], whose sharp bound 16√2/9 is attained
by the oblate compacton; the paper's constant is c = (4/π)G* = 64√2/(9π).
Dropping the chiral cos²Θ term weakens the bound to π/√2 (axi3).

**Question.** Embed the corpus's W in the one-parameter family
W_λ = (1/8π)[(1−cos ψ)² + λ sin²ψ cos²Θ], λ ∈ [0, 1] (λ = 1 = GUM;
λ = 0 = achiral alternative).  Is λ = 1 *distinguished* — i.e. do the
closure identities the corpus leans on (E6 = E0 − I/16π at d = 1;
E_stat_tot = (3/2)G*; the compacton support law) persist across λ (generic
→ not distinctive) or hold only at λ = 1 (distinctive)?  What is G*(λ) and
does any *other* λ reproduce a corpus anchor?

**Method.** Numerical, reusing the frozen gstar parametrization (mandate:
same conventions as `gstar_solve.py`; nothing refitted).  (a) The sharp
lower bound B(λ) = (1/√2)∫∫ sin²ψ sinΘ √[(1−cosψ)² + λ sin²ψ cos²Θ] dψ dΘ
by high-order quadrature (target 1e-12); closed-form checks at λ = 0, 1.
(b) Attainment at sampled λ via the same per-θ BPS reduction (π³b² = W_λ
pointwise, direction-locked ansatz); verify saturation residual.
(c) Sector decomposition E6(λ), E0(λ), I(λ) on the λ-minimiser via the
pushforward quadrature; test the two identities and the support law at
each λ.  (d) Scan for accidental anchor hits: does c(λ) = (4/π)G*(λ) cross
any printed corpus constant (3.2011…, 𝔠₀ = 2.5147…·(4/π), π, etc.) at
λ ≠ 1?

**Gates (pre-registered).**
- **S2-G1** Endpoints exact: B(1) = 16√2/9 to ≤ 1e-10; B(0) = π/√2 to
  ≤ 1e-10 (both closed forms known independently).
- **S2-G2** Attainment: BPS saturation residual ≤ 1e-8 at every sampled λ
  (≥ 9 interior points), else the honest statement that attainment fails
  off λ ∈ {0, 1} and B(λ) is only a bound there.
- **S2-G3** Identity persistence table: E6 = E0 − I/16π and
  E_stat = (3/2)G* tested at every sampled λ with residuals printed;
  verdict per identity: GENERIC (holds ∀λ) vs DISTINCTIVE (λ = 1 only)
  vs OTHER (structure found).
- **S2-G4** Honest bottom line: what distinctiveness does / does not mean
  within-model; no claim that distinctiveness confers physical support.

## S3 — O-1 pilot: validated forward flow-map covering (Q1 technology extension)

**Component being extended:** obligation O-1 of Conjecture C
(`tier7-program/Q2/ANALYTIC_RECON.md` §5) — *the dominant obligation*:
early-time flow control on [0, T0] × supp ρ₀, for which "the only
identified technology is the Q1 validated-enclosure covering."  Q1
certified 7 single backward paths; O-1 needs a *covering* by validated
boxes, forward in time.

**Question (pilot scope, honest).** Can the Q1 Lohner/interval engine be
extended from single paths to a box covering at all, at what cost, and
what fraction of a pre-registered subdomain can be certified
crossed-or-localized by T0?  A pilot ≠ an O-1 discharge; the pilot's
deliverable is the *feasibility datum* (cost per box, wrapping-effect
growth, fraction certified) that decides whether full O-1 is tractable.

**Method.** Reuse `tier7-program/Q1/q1_enclose.py` (mpmath.iv Lohner QR,
Picard boxes, the certified band-limited field) with forward integration.
Subdomain (pre-registered): n_y = 1; the lag-band × mid-channel rectangle
w ∈ [−1.2, −0.6], z ∈ [0.05, 0.95] at t₀ = 1.0 (the Q2 T_ball band),
mapped to (x, z) boxes; M = 24 initial boxes (6 × 4); integrate to
T0 = 5.14 (the measured T_max(1)).  Per-box nohup + JSON checkpoint
(restart-resilient, per the P2 protocol).  A box certifies if its
enclosure either (i) wholly crosses x = d before T0, or (ii) remains
wholly localized (x < d) with printed enclosure at T0; a box FAILS if the
enclosure blows up (width > 1) before T0 — subdivide once (2 × 2) on
failure, then stop.

**Gates (pre-registered).**
- **S3-G1** Engine validation: forward-integrate the reverse of one Q1
  certified backward path; endpoint enclosure must contain the Q1 anchor
  (containment check, not distance).
- **S3-G2** Covering executed: all 24 boxes (plus ≤ 1 round of
  subdivisions) run to T0 or terminal blow-up, all checkpointed, zero
  uncertified gaps in the subdomain accounting.
- **S3-G3** Feasibility datum printed: fraction certified (i)/(ii)/FAIL,
  median validated-step cost per box, width-growth factor distribution —
  and the extrapolated cost of a full O-1 covering, honest.
- **S3-G4** Scope honesty: pilot explicitly NOT an O-1 discharge; the
  kill-bin theorem remains open; register unchanged.

## S4 — Population-scale certification screener (P2 extension)

**Component being extended:** the P2 exact-field certification (432/432
configurations, `tier7-program/P2/RESULTS.md`) toward the engineering open
"population-scale certification."

**Question.** Does the P2 stability margin survive at population scale?
Strategy: cheap float64 pre-screen of a large configuration population
drawn from the P2 parametrization grid refined ~23× (target N = 10,000),
rank by margin, then re-certify the worst decile boundary with the P2
high-precision ladder (K = 24 configs at 200-bit).

**Gates (pre-registered).**
- **S4-G1** Screener validation: float64 margins reproduce the P2 pipeline
  on 24 randomly chosen P2 configs to ≤ 1e-6 relative.
- **S4-G2** Population screen complete: N ≥ 10,000 configs, margin
  distribution + minimum located; nohup + chunk checkpoints.
- **S4-G3** Worst-K certification: the K = 24 lowest-margin configs
  re-run at 200-bit; verdict per config (STABLE/UNSTABLE) with certified
  digits; any UNSTABLE is a finding, not a failure of the workstream.
- **S4-G4** Honest statement: this is *screened* population coverage, not
  a formally-verified kernel; the formal-kernel open remains open.

---

## Execution discipline (all workstreams)

- Long runs: `nohup` + per-phase/per-item CLI checkpoints (P2 protocol);
  background watcher loops for monitoring; container-restart recovery from
  committed intermediates.
- Every workstream: `RESULTS.md` + machine-readable JSON + figure where
  meaningful; gates adjudicated in `TIER8_ADJUDICATION.md`.
- Spec defects discovered mid-run are coordinator-owned and printed.
- Commit/push after each major stage (credit-resilience).
