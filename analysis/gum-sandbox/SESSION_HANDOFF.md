# SESSION HANDOFF — GUM Replication Campaign + Physics-Core Program

**Audience: a future Claude (Fable-class) session resuming this work.**
Checkpoint timestamp: 2026-07-12, branch `claude/analyze-gum-ivm-viz` on
`johngmuender/fork_frankensim`. Everything below is committed; two agents
were still executing at checkpoint (§6 — resume there first).

## 1. What this is

Two consecutive programs on one branch, both complete except §6:

1. **The replication campaign** (Tiers 0–5a): adversarial, within-model
   replication of the GUM speculative-physics corpus's tagged numerical
   claims on frankensim's certified numerics. Index:
   `REPLICATION_CAMPAIGN_STATUS.md`. Product: findings ledger F-R1–F-R8.
2. **The physics-core engineering program** (ROADMAP_v3 executed in
   full; ROADMAP_v4_SCALE executing): the fs-gum-* crates under
   `gum-core/`, closing the five gaps of GAP_ANALYSIS_v3 and now the
   scale (CPU→GPU) frontier.

**The binding epistemic rule of everything here**: every claim is
WITHIN-MODEL (does the corpus's printed number solve its printed
equations?); nothing asserts anything about nature; every PASS/FAIL is
gated against frozen referees; defects are printed on both sides
(including our own — see the G* correction in
`DISCHARGE_PACKAGE/README.md` §7 and `tier4-field/TIER4_ADJUDICATION.md`).

## 2. The findings ledger (the campaign's product)

| ID | Severity | One line | Where |
|---|---|---|---|
| F-R1 | conservative slip | NR-D2 KE coefficient ×25 too large at stated v; kill fires harder | tier0-gauntlet/RESULTS.md |
| F-R2 | harmless | AUD-15 g-span inconsistent with its own formula | tier2-closure (Tier 2a) |
| F-R3 | band arithmetic | ⟨r10⟩ ±0.90 band implies ρ = −0.36 not +0.45 | tier0-gauntlet/RESULTS.md |
| F-R4 | **substantive** | G.5's SDiff inertia dial is frozen — proof + machine precision | axi_FR4_NOTE.md, axi_ADJUDICATION.md |
| F-R5 | **substantive** | The ħ-closure's variational problem selects halo saturation, not 𝔠₀ = 2.5147; the ⟨r1⟩ tuple is a halo-unstable saddle; confirmed unrestricted-3-D; clock-forced | axi3_*, tier4-field/*, DISCHARGE_PACKAGE/ |
| F-R6 | moderate | ⟨r7⟩→E-F1→⟨r8⟩ provenance chain arithmetically false as written | tier5-family/ |
| F-R7 | minor | K.2 "qR* ~ 1e-17" six orders from its own inversion | tier5-family/ |
| F-R8 | minor | locked-bond pull 1.55σ printed as 1.3σ | tier5-family/ |

**Exact results discovered along the way** (all in closed form):
- Locked halo thresholds κ_crit_paper = 1/√(2⟨sin²θ⟩_w): polar channel
  = 1 exactly (matches corpus onset 1.000±0.004); s=sin² channel
  = √(7/12) exactly (= corpus κ₀); ring limit = 1/√2.
- **G\* = 16√2/9 exactly** (G3, `tier2-closure/gstar_RESULTS.md`):
  the halo-saturated objective, attained by the newly-found **oblate
  compacton** g(cos(F/2)sinθ) = ρ³/(6√2); saturated closure
  𝔠_paper = 64√2/(9π) = 3.2011247 at κ = 1/√2 exactly. The 𝔠₀
  proximity is the accident 24√21 ≈ 35π (0.024%).
- NB corrected digit transposition: 𝔠₀ = 128√42/(105π) = **2.514754**
  (some earlier docs wrote 2.514735).

**What replicated cleanly (credit side, always state it)**: 52/52
certified archive arithmetic (root e19d0bcb…), four-branch spectrum +
M-1, ¼ and √2 invariants, spin selection, radial rung (tail mass
blind-confirmed), Born relaxation + equivariance + generic-M≳4, ⟨r10⟩
central value, ⟨r11⟩ quark belts fully reconstructed (0.17–0.32σ),
family-sector pulls (0.05σ, 0.43σ), soft-sector tilt.

## 3. Key tacit knowledge (unit maps & conventions — memorize before touching closure numbers)

- Sector conventions (axi_RESULTS.md): E2 = (1/4π)∫Σ|D_i q|²,
  E4 = (1/4π)∫Σ_{i<j}(|D_i|²|D_j|²−(D_i·D_j)²), E6 = π³∫b²,
  E0 = (1/4π)∫(1−q0), I = ∫2(q1²+q2²), b = −(1/2π²)det(q,∂x,∂y,∂z q).
  E_static = t(E2+E4)+E6+E0; ε=0.05 ⇔ t = 0.0082764349. Compacton
  R* = 2^(5/6); BPS floor E6+E0 ≥ 32√2/15 per unit degree.
- Unit map (compacton-exact, validated 7–9 digits):
  𝔠_paper = 2L/√(2π³); κ_paper = 2√π·κ_ours, κ_ours = L/I.
  **κ_paper is exactly ω/(√2μ)** where μ = tilt-channel mass = 1/√(8π)
  in engine units — hence the halo threshold at κ_paper = 1/√2.
- Clock: L² = (2/3)·I·E_static ⇒ E_rot/E = ¼ identically.
  κ(L_clock)/threshold ≈ 1.25 (over-spun) on every B=1 solution.
- Halo instability: dR = [∫η²dV](1/(8π) − κ_ours²) for tilt halos;
  direction-locked profile halos get the ⟨sin²θ⟩_w generalization.
- Benchmark targets (⟨r1⟩): 𝔠 = 2.37±0.09, κ = 0.802±0.018,
  g* = 1.31±0.04, λ* = 0.42±0.05, V = 1.409±0.010, onset 1.000±0.004.

## 4. The engineering program state (gum-core/)

All crates standalone (empty `[workspace]`, path deps), all gates
bit-replay verified. ROADMAP_v3: **executed in full** (Phases A, B1–B4,
E1–E2, F1–F3). ROADMAP_v4_SCALE: G-A, G2, G3 done; G1, G4 in flight
(§6); G5 queued.

| Crate | Gates | Golden root / note |
|---|---|---|
| fs-gum-field | 25/25 | 6c1e7850… — Field3, stencils, sectors, radial solver |
| fs-gum-topo | 25/25 | 9cfbe9a6… — degree, Hopf ×2 methods, disclinations |
| fs-gum-cosserat | 26/26 | 9ddc0fe4… — dispersion, Verlet, W_χ, Dzyaloshinskii |
| fs-gum-sde | 21/21 | Nelson ≈37× faster than Bohm (τ 0.132 vs 4.853) |
| fs-gum-statics | 18/18 | analytic gradients (sextic 1e-8), guarded ANF, halo referee |
| fs-gum-e2e | 22/22 | d30de830… — the composed pipeline |
| fs-gum-gpu | all | WGSL f64 on llvmpipe; bands 100× margins; K4b identical descent decisions |
| fs-cosserat-pilot (4B) | 13/13 | f88731af… — analytic-config diagnostics |
| tier0-gauntlet | 52/52 | e19d0bcb… |

Environment facts a future session needs: 4-core Cascade-Lake-class
Xeon, AVX-512 present but `-C target-cpu=native` is a MEASURED 10–25%
regression on the real kernels; multiprocess scaling 3.96×; llvmpipe
software Vulkan at /usr/share/vulkan/icd.d/lvp_icd.json (wgpu/naga
30.0.0 resolve through the proxy); scipy 1.17.1 + numba installed;
**asupersync sibling repo ABSENT** → fs-la/fs-fft need
`default-features = false` from out-of-tree (fixed in F3; see
gum-core/PORTABILITY.md); the fork's workspace root has never resolved
in this checkout.

## 5. Working practices (keep these)

- **Two-layer authorship**: agents write RESULTS.md (measured facts,
  no adjudication language); the coordinator writes *_ADJUDICATION.md
  (findings, severities, discharge paths). Never blur them.
- **Print-your-own-defects**: task-sheet errors, harness bugs, and
  wrong hypotheses get documented in the committed record (examples:
  the Tier-1 over-strict specs, the μ spec error, the direction-locking
  hypothesis refuted by 4C, the G*≈𝔠₀ speculation corrected by G3).
- Gates before physics; referees before solvers; one deliberate golden
  bump only with a documented reason (fs-gum-kern's tiled order).
- Commit style: descriptive multi-paragraph messages; trailer
  `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>` + session
  link; push to `claude/analyze-gum-ivm-viz` only; no PRs unasked.
- Background agents' completion waiters are FLAKY: if a run's outputs
  exist but the agent went quiet, send it a message (it resumes from
  transcript). Snapshot WIP commits keep the stop hook satisfied.
- Determinism doctrine: bit-identical goldens within a backend;
  tolerance-band/metric goldens across backends and at ANF endpoints
  (arrest cascades amplify ulps).

## 6. IN-FLIGHT AT CHECKPOINT — RESOLVED before session end

**Both in-flight builds completed and are committed.** SUBSEQUENT
SESSIONS ALSO COMPLETED: ROADMAP_v4 G5a/G5b/G5c (N=192 ladder with the
crossing OBSERVED at it 2070-2080; two-knot well 5.0σ, x₀ = 1.95 ±
0.07 vs corpus 1.92 ± 0.08); NATIVE_REMEASURE (host-specific
regression resolved); Phase H theory audit (48 arguments; T-ledger
T-H1..T-H11; the repaired-closure Theorem T3.1) and the FULL H2
numeric program (F-R9 promoted — NR-K3a second moment, 13-18 orders;
F-R4 floor confirmed on the collision problem; E.6 positivity
T-H10; D1b honest floor T-H11; A2b counterexamples survive; the
anchor discriminator and torque-sign gates close IN THE CORPUS'S
FAVOR). Ledger now F-R1-F-R9 + T-H1-T-H11. FINAL ADDITIONS (2026-07-16 close): H4 executed — the fatal branch
REFUTED (∂ = ±2·[rot], two routes); T-H4 discharged with corrections
(F-R16); ledger FINAL at F-R1–F-R16, T-ledger EMPTY. And corpus2/ —
the complete externally-revised edition: Omega paper v3.0-ext (67/67
map entries), Course/Primer/Teachers-Ed v2-ext (35/35 + 14), manifest
v4.0-ext (105-file census: 3 revised/67 delta/35 carried), Session Map
v2.0-ext, Watch-Mode v2, 09-External-Replication-Record,
10-The-Repaired-Closure, 70 S-30 delta memos (98 annotations), all
H4-pending boxes patched and verified. The program is COMPLETE; any
new session starts from corpus2/README.manifest-v4.0-ext.md and the
corpus's response. G1: 13/13 gates,
serial↔4-thread bit-identical, 3.61× speedup. G4: full protocol run
(23 relaxations, 3 orientations), bond loop IN BAND — closed-loop
x₀ = 1.88/1.97 vs predicted 1.90 ± 0.05; direct 2.0 ± 0.15 vs corpus
1.92 ± 0.08; 𝔟_eff absolute normalization not recoverable (documented
spec-underdetermination, 𝔭-class). See ROADMAP_v4_SCALE execution
record and fs-gum-twoknot/twoknot_RESULTS.md.

**UPDATE (2026-07-16): G5 also complete** — G5a (N=192 F-R5 ladder:
state-space path grid-robust, clock ladder resolved to a state
function with floor ≈1.245–1.2525, fixed-h box axis a NULL correcting
4A's plank; fs-gum-kern/n192_RESULTS.md) and G5b (two-knot refinement:
well 5.0σ, x₀ = 1.95 ± 0.07 vs corpus 1.92 ± 0.08 at 0.3σ;
fs-gum-twoknot/twoknot_REFINE.md). DISCHARGE_PACKAGE carries the G5
refinements + corrections. ROADMAP v4 is EXECUTED in full; the next
session starts from ROADMAP_v4_SCALE's "Next (unscheduled)" list
(the ~2500-cap N=192 crossing run is the cheapest open item; then
explicit SIMD, real-GPU execution, blue-fog, the corpus's response).
Ops note: container restarts kill background farms/monitors — bank
anchors first, use stage-flushed JSONs, and expect to wake agents
manually (this bit us three times today). The subsections below are retained as the
historical checkpoint state.

### (historical) in-flight state at the checkpoint commit

**G1 `fs-gum-kern`** (deterministic parallel kernel layer): all source
modules committed (tile.rs TILE_I=4 canonical order, reduce.rs pairwise
tree, engine.rs gather-form sweeps, anf.rs, gates binary); the
freeze/gates run was executing at checkpoint; RESULTS.md is a
PLACEHOLDER. To finish: run
`cd gum-core/fs-gum-kern && cargo run --release --bin gum_kern_gates`,
verify (a) old-engine equivalence at machine-eps class (the documented
golden bump ~12th digit), (b) serial↔threaded bit-identity at 1/2/3/4
threads, (c) speedup ≥3× at N=96 corner-gradient, (d) the N=48 descent
referee; then write RESULTS.md (gate table, bit-contract, golden-bump
numbers, speedup table), adjudicate, commit. Spec: the G1 agent prompt
(recoverable from ROADMAP_v4_SCALE Phase G1 + scale_survey_v4.json
dims 1–2).

**G4 `fs-gum-twoknot`** (the I.2 bond-equation flagship — the corpus
protocol never executed anywhere): crate committed (field_a.rs
anisotropic wrapper, seed.rs product ansatz, engine_a/anf_a, runner +
gates binaries, farm.sh, analyze_twoknot.py). At checkpoint: single-knot
reference + attractive-channel separations m14/16/18 DONE (runs/*.json),
m20/23/26/29 EXECUTING as a 4-process farm (~10¹³ FLOP each). To
finish: wait/rerun `./farm.sh` for missing separations (binary:
`./target/release/twoknot_run attract <m>` where m = 10×x roughly —
check farm.sh), then `python3 analyze_twoknot.py` → E_int(d) table,
screened-dipole fit (μ = 7.7725 engine units), bond-equation closure →
**x₀ vs predicted 1.90 ± 0.05 (𝔟 = 42) and corpus-measured
1.92 ± 0.08**; write twoknot_RESULTS.md (+ orientations 2–3 as budget
allows), adjudicate (this is a NEW corpus replication — F-R-ledger
grade if it fails, credit-side if it lands), commit.

**G5 (queued)**: N=192 F-R5 halo endpoint on fs-gum-kern's threaded
sweeps (8× the campaign's cells) — hardens DISCHARGE_PACKAGE.

## 7. Remaining roadmap beyond G5 (unscheduled)

Blue-fog condensate/disclination boxes + tower-depth at scale (the
fs-gum-gpu backend is the on-ramp; needs a real GPU); gravitation
sector; nucleation; the corpus's response to `DISCHARGE_PACKAGE/`
(watch for the frozen ⟨r1⟩ code — the discharge routes are in its
README §8).

## 8. File map (fast orientation)

- `REPLICATION_CAMPAIGN_STATUS.md` — campaign index (Tiers 0–5a).
- `GUM_SIMULATION_ASSESSMENT_v2_FABLE.md` — the original tiered plan.
- `substrate-suite/01-GUM-Omega-Paper-v2.0.1.md` — the corpus's flagship
  paper; `corpus/` — the 105-file archive (no ⟨r1⟩–⟨r11⟩ code in it).
- `tier2-closure/` — the closure chain: radial → axi (Step 2, F-R4) →
  axi3 (Step 3, F-R5) → axi4_locked (4C) → gstar (G3 exact).
- `tier4-field/` — 3-D pilots + TIER4_ADJUDICATION.md (F-R5 in 3-D).
- `tier5-family/` — family/quark sector audit (F-R6/7/8).
- `DISCHARGE_PACKAGE/` — the corpus-facing deliverable (verified
  reproduction commands in REPRODUCE.md).
- `gum-core/` — ROADMAP_v3/v4 + GAP_ANALYSIS_v3/v4 + surveys (JSON) +
  all fs-gum-* crates + PORTABILITY.md.

Every RESULTS.md is self-describing; every gates binary re-verifies its
crate in minutes. Trust the committed roots; recompute before extending.

## Addendum (2026-07-16, end of session-X continuation): Phases I and J closed

Phase I executed and folded (DISCHARGE_PACKAGE v2.2 with F-R16/F-R17 +
phase_I block; corpus2 carries ⟦I4 reported⟧ boxes at every residual-opens
site). Phase J (the filed obligations) executed and folded: J1 recovered
the 𝔭/𝔟_eff convention (edge-referenced class; OPEN → RECOVERED-AS-CLASS
in corpus2), J2 executed the repaired (4.8′)/(4.9′) arithmetic (ceiling
2.605×10⁻³; ε-anchor forfeit 11.283σ; dichotomy survives; two stale
sentences in OUR corpus2 revision corrected in print), J3 confirmed the
Majoron battery SHIFTS-HARMLESSLY (17/17). Ledger remains F-R1–F-R17;
T-ledger EMPTY. The former named residual — the multi-ε saturated re-scan — was then
executed as J4: exponent measured LINEAR (p = 0.997 [0.996, 1.012]; 2/3
excluded at χ² ratio 8445; ×4.3/60 MeV/7.8×10⁻¹⁰ all MEASURED; dichotomy
reading-robust). NO named opens remain; only the unscheduled Beyond items. Artifacts:
theory-audit/j{1,2,3}_* + ROADMAP_v6_SYNTHESIS.md Phase-J verdicts;
audit-logs/ archive current through the Phase-J sweep (35 mapped agents).

## Addendum 2 (2026-07-16, credits-exhaustion checkpoint #2): Phase K state — RESTORATION GUIDE

You are (probably) a fresh Claude Fable session resuming this campaign
after credit exhaustion. Read this section, then REPLICATION_CAMPAIGN_
STATUS.md, then theory-audit/ROADMAP_v7_HARDENING.md. Branch:
claude/analyze-gum-ivm-viz (develop and push ONLY here; no PRs).

### Where Phase K stood at checkpoint time
- **K2 ✅ DONE, coordinator-verified, committed** (5c9f47a + re-verify):
  linear ε-law holds to ε = 0.001; 2/3 excluded by the new window alone
  (χ² 77); J4 STRENGTHENED; propagation UNCHANGED. Headline defect K2-D1:
  the 12-mode family's intercept bias DRIFTS (~+1.7e-4 at ε = 0.001,
  +3.15σ vs the constant model) — intercept-shaped, not exponent physics.
  k2_smalleps.py resumes from banked points in ~2 s (verified).
- **K1 ▶ RUNNING at checkpoint**: the /6 vs /(15π/8) discriminator.
  Skeleton k1_RESULTS.md (5 sections PENDING); k1_discriminator.py +
  k1_results.json (stage-flushed) committed as WIP; its run log was at
  "STAGE B5 — C6 and mu (amplitude-free anchors, re-verified)". If the
  agent died: read k1_results.json's banked stages, re-run
  k1_discriminator.py (deterministic), and finish the memo from its
  output. Remember: a verified NULL ("indistinguishable within the
  archive") is a legitimate terminal verdict — do not force a winner.
- **K4 ▶ RUNNING at checkpoint**: the REPRODUCE freshness gauntlet.
  DISCHARGE_PACKAGE/FRESHNESS.md had 24 rows PASS / 25 PENDING; queue
  farm (13 k4run.sh jobs; queues A cheap/B mid/C mid/D monsters) with
  per-command logs + timings.csv in the session scratchpad k4/ dir
  (NOT in the repo — logs are gitignored and the scratchpad dies with
  the container). One finding so far: K36 CMD-DEFECT (a REPRODUCE
  command as printed fails; exit 2 at 0.5 s in timings.csv) — the
  corrected K36b (fs-gum-kern n192_fr5 crossing, cap 2500) was running
  (the multi-hour long pole). If the farm died with the container:
  FRESHNESS.md's PENDING rows are exactly the unfinished work-list;
  re-run rows by hand from REPRODUCE.md (all paths relative to repo
  root), update rows, and write the findings section + freshness
  verdict. Rust rows: baseline flags, NOT -C target-cpu=native.
- **K3 NOT STARTED** (by design, after K1/K2/K4): WGSL f64
  sector-measure + E_static gradient vs fs-gum-statics CPU goldens on
  llvmpipe (tolerance-band, NOT bit-identical across backends).

### Adjudication debts a future session owes
1. K1 verdict → fold into ROADMAP_v7 (verdict section), campaign
   status, and IF decisive: corpus2's J1 box (upgrade RECOVERED-AS-CLASS
   → a named convention) + DISCHARGE_PACKAGE CLAIMS.json.
2. K4 FRESHNESS.md finalization → commit; CMD-DEFECT findings also
   belong in DISCHARGE_PACKAGE README (a defects-of-the-package note).
3. K2 fold: ROADMAP_v7 verdict section + a one-line amendment where
   J4 is cited (K2-D1 nuance), CLAIMS.json phase_K block when the
   phase closes, REPRODUCE rows for k1/k2 scripts.
4. Audit sweep after each fold (tar/xz commands in §5 of this file;
   MANIFEST agent map — add the three K agents, phase "Phase K").
5. Task tracker: #53 (K1) and #55 (K4) in_progress at checkpoint.

### Working knowledge worth keeping (beyond §5's playbook)
- Re-verify discipline unchanged: coordinator re-runs every gate script
  before final commit; physics fields must be identical (timing fields
  exempt). COMMIT THE AGENT'S FINAL JSON BEFORE re-running scripts that
  rewrite their results file (the J4 mid-rewrite lesson, f35c1ae).
- Agents park with nohup children alive; the harness notification fires
  anyway. Check pgrep, set a `while kill -0 PID` background watcher,
  and SendMessage the agent when its process exits (J4/K2/K4 pattern).
- K4's re-runs legitimately rewrite committed result JSONs (timing
  fields) — verify timing-only before shrugging; h24/h25/gstar/tier1
  REPORT and j1_results.json were confirmed such cases.
- Scripts with per-point flushing resume cheaply; prefer wiring resume
  into new long scripts (k2_smalleps.py is the template).
