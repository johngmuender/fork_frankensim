# SESSION CHECKPOINT — Tier 6 Foundations (session 536b52e3)

**Audience: a future Claude (Fable-class) session restoring this line of
work.** Written 2026-07-17 on branch `claude/analyze-gum-po`
(`johngmuender/fork_frankensim`). Everything referenced is committed and
pushed. Read this file, then `FOUNDATIONS_ANALYSIS.md`,
`TIER6_ADJUDICATION.md`, and the top-level `../SESSION_HANDOFF.md`
(the ivm-viz line's handoff) — that set restores the full working state.

## 0. Branch topology (read first — two lines, deliberately kept in sync)

- `claude/analyze-gum-ivm-viz`: the original campaign line (Tiers 0–5a,
  theory-audit Phases H/I/J/K, gum-core crates, DISCHARGE_PACKAGE,
  audit-logs for session d0b81976). Phase K is COMPLETE there
  (K3 folded, K28 corrected 5/5 bit-identical, gauntlet 46/46 zero
  skips); only the unscheduled "Beyond" pair remains.
- `claude/analyze-gum-po` (THIS branch): forked from ivm-viz, adds
  **Tier 6 (foundations)** — this directory — and has ivm-viz's Phase-K
  completion **merged in additively** (merge commit `0fcd532`; no
  collisions; both lines' artifacts coexist).
- A legacy branch `claude/analyze-gum-po-5x5cne` exists from session
  bootstrap; it is NOT the working branch. Develop and push on
  `claude/analyze-gum-po` unless the user redirects.

Commit spine of this session: `e6caebf` (fork point) → `97c9afa`
(analysis + roadmap) → `71bbd1c` (L1/L2/L3/L5 execution + adjudication)
→ `c167971` (L1 trajectory dump) → `0fcd532` (Phase-K merge) →
checkpoint commit (this file + audit-log archives).

## 1. What this session did (chronological)

1. **Read the corpus in extreme detail**: corpus2/01-GUM-Omega-Paper
   v3.0-ext (all 6 parts + embedded AUD-15), REPLICATION_CAMPAIGN_STATUS
   (F-R1–F-R17, Phases J/K), REVISION_CHARTER, corpus2 layout.
2. **Read four foundations texts in extreme detail** (user-uploaded PDFs,
   ephemeral — see §3 for restoration-grade summaries):
   Bell 1982 "On the Impossible Pilot Wave"; Goldstein "Bell on Bohm";
   Bricmont 2019 Hvar lectures; Bedingham–Dürr–Ghirardi–Goldstein–
   Tumulka–Zanghì GRWm (arXiv:1111.1425).
3. **Wrote `FOUNDATIONS_ANALYSIS.md`** — GUM vs the canon; findings
   F-T6-1…6 (see §4).
4. **Wrote `ROADMAP_v8_FOUNDATIONS.md`** — workstreams L1/L2/L3/L5 with
   pre-registered gates (L4/L6 analysis-only, folded into the memo).
5. **Executed the roadmap** as a 4-agent parallel workflow
   (`wf_5852a821-c22`, archived in `../audit-logs/`): all four agents
   wrote + ran independent Python (numpy/scipy/matplotlib, installed
   fresh via pip — NOT preinstalled in this container image).
6. **Adjudicated** (`TIER6_ADJUDICATION.md`): 14/16 gates PASS, 2 honest
   PARTIALs assigned to roadmap-spec defects (coordinator-owned).
7. **Merged** ivm-viz Phase-K completion additively; **checkpointed**
   (this file + `../audit-logs/MANIFEST_tier6.json` + xz archives).

## 2. Tier-6 results at restoration grade

| WS | Verdict | Key numbers to trust |
|---|---|---|
| L1 flipped SG | 4/4 PASS | 4000/4000 side-preserved + label-flipped; ρ,j bitwise invariant under k→−k (max dev 0.0); 0 median crossings; Born 0.92σ; dt-halving 5.9e-9 |
| L2 ToF | 5/5 PASS | dBB law RMS 2.6e-12 vs X(0)√(1+t²); ToF KS p=0.77 vs |Ψ̂|² with v(t=0)≡0; Nelson equivariance ≤0.63× floor; Nelson ToF p=0.78; ν-blind (pairwise KS ≥0.34 over ν∈{0.1,0.5,1}) |
| L3 leakage ladder | 5/5 PASS | slopes 1.000000 / 2.999844 / 4.999869 (ln P vs ln 1/c_L); flux R-indep 2.5e-4; lemma: conservation kills ℓ=0,1 ⇒ (c/c_L)⁵ restored; residual assumption isolated (source = conserved knot densities) |
| L5 arrival times | CONFIRMED (G2/G3 PARTIAL) | transverse hard cutoff τ_max = 2.28503 (stable 7e-7 dt-halving, 1.6e-4 Nx-doubling); axial support to 3.2245 (41% past cutoff, zero transverse events there); far-field KS p=0.31; engines cross-validate at KS p=1.0, norms/walls/energy ~1e-16 |

L5's PARTIAL causes (both spec artifacts, mine): (i) one trajectory
(0.05%) born on the discontinuous ±3σ truncation edge is resolution-
non-convergent — spec should have used a C¹-smoothed truncation;
(ii) G3's 5% tail-fraction bar was set blind against the k₀=4 boost
(measured 2.8%; direction + separation confirmed).

**Named follow-ups (the natural next roadmap items):**
- L5′: rerun with C¹-smoothed truncation and slower packet / farther
  near-detector → convert both PARTIALs.
- The four proposed corpus annotations (TIER6_ADJUDICATION §"Proposed"):
  VIII.A architectural credit; VIII.F A3-vs-collapse-class sharpening +
  mutual kill; VI.E/VIII.B DGNSZ-2013 anchor; VIII.D Lock-2 lemma
  adoption. These are OFFERS per the corpus2 charter — adoption is the
  authors'; a corpus2 fold (v3.1-ext delta or tier6 pointer notes) is a
  candidate next session task.
- L6 full transcription (PLC selection principle applied to VI.D's P1
  premise) — deferred, needs the induced-gravity sector.
- POVM-exclusion theorem for the L5 cutoff statistic (T4-W5 gap) — the
  hard open; L5 is its testbed.

## 3. Foundations texts — restoration summaries (PDFs were ephemeral uploads)

**Bell 1982, "On the Impossible Pilot Wave" (Found. Phys. 12, 989–999).**
Minimal dBB model: Ψ(a,x,t), measurement Hamiltonian gO(ħ/i)∂/∂x
(infinite-mass ⇒ it is the whole H), packets Φₙ(x−gOₙt) separate;
add X(t), dX/dt = j/ρ; probability enters once (initial conditions).
Autopsies: von Neumann's additivity postulate absurd for individual
results (±½ ≠ (±½±½)/√2); Gleason–Jauch/Kochen–Specker (anticipated by
Specker 1960) hides non-contextuality — dBB violates it because
different completions are different Hamiltonians ("denial of Bohr's
insight"); Jost's identical-particles proof dissolves since Ψ carries
the symmetry. Morals: test on simple models; only position observations;
ban "measurement" (quantum logic = artifact of the word); only
local-causality impossibility results survive, and dBB is what started
that line.

**Goldstein, "Bell on Bohm" (memorial chapter).** One Bell, not two:
1966 autopsy → noticed dBB's explicit nonlocal mechanism → 1964 theorem
("even smarter than Bohm won't remove nonlocality"; publication order
inverted by delay). Nonlocality's origin = ψ on configuration space.
Wigner/Hawking misread documented. The Big Question (Bell's last word):
which of dBB/GRW can be made Lorentz invariant. Progress ledger: Tumulka
flash-GRW 2006; Bedingham et al. 2014 (GRWm, below); DGNSZ 2013
covariant ψ→foliation map for Bohmian mechanics.

**Bricmont 2019 Hvar slides (84 pp).** State = (X, Ψ); Schrödinger
(never collapses) + guidance dXₖ/dt = ∇ₖS = Im(Ψ*·∇ₖΨ)/Ψ*·Ψ;
equivariance from continuity ⇒ quantum equilibrium ρ = |Ψ|². Flipped
Stern–Gerlach: same Ψ, same X(0), reversed gradient ⇒ same spatial
motion, flipped outcome — "measurements don't measure." Momentum ToF
(App. 1): Ψ real Gaussian ⇒ S=0, all particles at rest; free evolution
gives X(t) = X(0)√(1+t²), p = lim X/t = X(0) distributed as
|Ψ̂(p)|² = π^{-1/2}e^{-p²}. Effective collapse: empty branch never
re-overlaps for N~10²³ (no environment needed); decoherence alone puts
the observer back (no fact distinguishes branches without X).

**Bedingham et al. 2014, GRWm (arXiv:1111.1425).** PO problem for
collapse theories: flashes vs matter density m(x,t) = Σᵢmᵢ∫δ³(x−qᵢ)|ψ|².
Eq. (1) not Lorentz invariant. Their law: m(x) = ⟨ψ_PLC(x)|M(x)|ψ_PLC(x)⟩
— evaluate on the past light cone of x (Eqs. 5–7; Dirac 4-vector form
Eqs. 8–9). Assumptions (A1) Tomonaga–Schwinger ψ_Σ + conditional law;
(A2) lightlike limits; (A3) no dependence on spacelike external fields;
(A4) macroscopic superpositions collapse at Born weights; (A5)
hypersurface-consistent outcomes. Results: m(x) hypersurface-independent
(required of matter); no-micro-signalling (10); pointer/ψ outcome
agreement (11); empirical adequacy (12) ⇒ nonlocal by Bell. Examples:
detection updates m at distance only after light delay (∫m d³x not
conserved); EPR: naive Eq.-(1) reading holds on post-measurement Σ₁ but
fails on wing-threading Σ₂/Σ₃. FLC variant rejected only by
no-micro-signalling + counter-intuitive retro-dependence. Applied to
Tumulka's process, empirically equivalent to flash PO. Answers Maudlin:
he assumed Eq. (1) in a frame, never a law of form (5).

## 4. Tier-6 findings register (final state)

- F-T6-1 credit: GUM's von Neumann-immunity nontrivial (T-B1 selection).
- F-T6-2a/b: canonical exhibits missing → **discharged by L1/L2**.
- F-T6-3: GUM↔GRWm discriminator ledger; all lines already armed in
  Sec. XII; NEW mutual kill: A3-positive + confirmed collapse
  contradict jointly.
- F-T6-4: GUM = matter-sourced member of covariant-foliation family;
  DGNSZ 2013 anchor missing at VI.E/VIII.B (documentation-grade).
- F-T6-5: T4-W3 lemma → **discharged at toy grade by L3**; residual
  assumption isolated.
- F-T6-6: Nelson-reproduces-dBB never exhibited → **discharged by L2**.
- F-T6-L5: A3 phenomenon confirmed on independent engine; τ_max = 2.285
  (this geometry); two spec artifacts, follow-ups named.

## 5. Tacit knowledge for restoration

- **Epistemic register is binding**: every claim within-model; "nothing
  bears on nature" notice on every RESULTS.md; credit side always
  stated; defects printed on both sides including our own (L5 spec
  artifacts are coordinator-owned, in the campaign's two-harness-defects
  tradition). Grade-language locks (RT-8) apply if touching corpus text.
- **Environment**: remote container; python numpy/scipy/matplotlib must
  be pip-installed fresh each session; poppler-utils (pdftotext) needed
  apt-get (mirror 404s → run `apt-get update` first). No `gh` CLI — use
  GitHub MCP tools. Scratchpad is ephemeral.
- **Corpus geography**: corpus2/ = replicators' proposed v3.0-ext
  edition (charter rules binding, esp. S-30 annotate-don't-rewrite and
  rule 4 grades-never-rise-by-replication); theory-audit/ = h2x/h3x/T3/
  T4/H4/I4/J/K artifacts; DISCHARGE_PACKAGE/ = certified claims +
  REPRODUCE + FRESHNESS; tier0–tier5 = campaign tiers; tier6-foundations
  = this session; audit-logs/ = session archives (MANIFEST.json for
  d0b81976, MANIFEST_tier6.json for this one).
- **The paper's current center of mass** (if resuming corpus work):
  saturated closure 𝔠 = 64√2/(9π) = 3.2011 at κ = 1/√2 [T3.1];
  restricted-branch record (4.5); ¼ and √2 invariants exact; F-R14
  gravitation dispute OPEN (P-acoustic rescue priced); F-R17 exchange
  lift FAILS (χ_exch = χ(σ)·χ_rot, new free ℤ₂); S1 stake live
  (Σm_ν ∈ [0.058, 0.11] eV, NO).
- **User constraints**: Claude credits may be limited — checkpoint
  early/often, document aggressively, push everything to GitHub
  (`git push -u origin claude/analyze-gum-po`); commits carry the
  Co-Authored-By + Claude-Session trailers.

## 6. Log archives (this session)

In `../audit-logs/`:
- `tier6-coordinator-session.jsonl.xz` — main session transcript
  (snapshotted during checkpointing; final commit exchanges absent by
  construction — git history carries the residue).
- `tier6-workflow-wf_5852a821-c22.tar.xz` — workflow journal + all four
  agent transcripts + meta (L1/L2/L3/L5).
- `tier6-workflow-script.js` — the exact orchestration script (meta,
  schema, four agent prompts) as persisted by the runtime.
- `MANIFEST_tier6.json` — index of the above with agent IDs/verdicts.

To restore deepest context: read this file + FOUNDATIONS_ANALYSIS.md +
TIER6_ADJUDICATION.md + the four L*/RESULTS.md; decompress the archives
only if the precise agent-level reasoning is needed.

---

## ADDENDUM (same session, post-checkpoint): L5′ EXECUTED

The first named follow-up of §2 is done (user-directed). `L5prime/`
(F-T6-L5-EXEC-2): pre-registered fixes (C² smoothstep truncation;
k₀ = 2), full gate battery, no exclusions. **Both L5 PARTIALs
converted**: G2 PASS (τ_max = 5.130950, stability 8.3e-8/5.6e-10,
all 2000 trajectories converge — edge-outlier mode eliminated;
energy grid-converged); G3 PASS (axial-beyond fraction 11.08% > 5%).
**G4 = FAIL-AT-d=25, honest**: coverage-balanced KS p = 9.9e-8 at the
spec'd d = 25, with a measured monotone KS(d) decay (0.126 → 0.057,
p → 3.7e-3 at d = 35) — the far-field null is approached, not
violated; at slow boost d = 25 is not asymptotic. New bench-design
record: (k₀, d_far) is a joint Layer-1 design constraint. Full
adjudication: TIER6_ADJUDICATION.md addendum. Updated open items:
POVM-exclusion theorem (T4-W5) stays the hard open; optional KS(d)
extension to the measured asymptote; corpus annotations 1–4 still
offered. Session archives updated: `../audit-logs/`
tier6-agent-L5prime.jsonl.xz + refreshed coordinator archive;
MANIFEST_tier6.json amended.

---

## ADDENDUM (same session): L7 POVM-EXCLUSION TESTBED EXECUTED (T4-W5)

The hard open of §2 now has its testbed-grade discharge (user-directed;
pre-registered at commit e21a5a7; workflow wf_1fc03d62-b33). Schema: a
fixed detector = one Naimark POVM ⇒ over a preparation family the
statistics must be quadratic in ψ; for the spin family, affine in n_y.
- **L7b (headline)**: on the L5′ engine, arrival statistics across
  n_y ∈ {−1…+1} (9 pts × 2000 identical-ensemble trajectories, one
  field solve, v = a + n_y·b exact) violate affinity at up to **40.2σ**
  (8 bins > 5σ; statistics even in n_y); **affine-vanishing kill**: two
  bins empty (Π < 1.5e-3, 95% CL) for all |n_y| ≥ 0.75 but Π = 0.012
  at n_y = +0.25 (affine bound ~1.9e-3 — contradiction ×7), zeros
  confirmed at double Nx + extended box; τ_max(n_y) 15.84 → 5.13
  monotone. **No single POVM reproduces the family.** F-T6-L7b.
- **L7a**: exact running-max CDF over spatial superpositions violates
  the quadratic form by 1.84e-4 (6×10¹⁰ × the 2.9e-15 comparator
  floor), exactly and only in backflow bins; backflow-free bins at
  floor. Deviations (honest, printed): registered 7×8 grid missed the
  backflow pocket (coordinator spec defect; refined 25×32); G4 control
  adjusted once, continuum has pockets at every separation. F-T6-L7a.
- Adjudication + proposed annotation #5 (VIII.F/IX.B may cite the
  testbed at [DW, testbed-grade, configuration-specific]):
  TIER6_ADJUDICATION.md L7 addendum. Residual to theorem grade:
  continuum proof of the cutoff zeros; generality beyond this
  configuration.

---

## ADDENDUM (same session): CORPUS3 — THE TIER-6 FOLD EDITIONS

User-directed. `../corpus3/` now contains the replicators' proposed
Tier-6 fold: **Omega paper v4.0-ext** (15 ⟦T6⟧ insertions: revision
note; Sec. 0 record; I.D note; Sec. III exhibits; proposals #1–#5 at
VIII.A/VI.E/VIII.D/VIII.F; IX.A/IX.B; X record; XI register; XII
mutual-kill; App. E.5; References), **Course v3-ext** (9 teaching
boxes), **Primer + Teachers' Edition v3-ext**, **Session Map
v3.0-ext** (§18), **Watch-Mode memo v3-ext** (§10 dated status),
**Replication Record v2.0** (§7), **README.manifest-v5.0-ext** (new),
under **REVISION_CHARTER_v2** (additive-only, ⟦T6⟧ marker, payload
authority = T6_FOLD_PAYLOAD.md, grades-never-rise, offers-not-
adoptions). Verification: every deletion across all seven editions is
a version-header line (FOLD_RECORD.md table); all numbers
payload-sourced. Two coordinator fixes post-agents: gate-tally
harmonization to per-workstream form (instruction-vs-payload mismatch
was mine); running heads bumped. Workflow wf_2c33ba17-f9b (4 agents,
399,795 tokens), archived in ../audit-logs/. corpus2 untouched (S-30).

---

## ADDENDUM (same session): TIER 7 — PROGRAM EXTENSIONS EXECUTED

User-directed. `../tier7-program/`: PROGRAM_ANALYSIS_T7.md (sector
scorecard; ranked open register; critical path), ROADMAP_v9_PROGRAM.md
(M1–M4, pre-registered), executed via workflow wf_95da92d5-375,
adjudicated in TIER7_ADJUDICATION.md. Headlines:
- **M1** (WS4-M4c toy): tracker identity w = −1 + Ω_m dynamical
  (1.8e-12); Γ ∝ Hⁿ ⟹ W = 1−n; the corpus's memory pairs recovered as
  W(z=0)→W(z=2) roundings; NEW minor flag **F-T7-1**: the "(+0.36 at
  β=−0.3)" CPL print unrecoverable + inconsistent with W ≥ 1; EdS
  fixed-point pathology measured (Δw(0) +0.21..+0.31).
- **M2** (servo law): form replicated (g² slope 1.9932, ω*-linear,
  integrable marginal); coefficient c = η/2 = pure bath convention;
  corpus's 0.021 ⇔ η = 0.042; super-Ohmic destroys the g² law ⇒ the
  printed law pins the r2 bath to the Ohmic class.
- **M3** (F-R9 route b): **CLOSED, NOT-VIABLE** — saturation can't
  change sparse support; R log-slow vs mean quadratic collapse (25–46
  decades traded per decade); repair menu narrows to routes (a)/(c).
- **M4** (Nelson relaxation): **17×/41× faster** at ν = ħ/2m than the
  dBB-calibrated tier3 rates; τ ~ ln(1/ν) — dBB a log-singular limit;
  III.D rates = direction-signed conservative bounds (double edge).
- Two coordinator spec defects printed (M1 Γ-profile; M4 continuity
  clause — the wrong pre-registration that produced the discovery).
Natural next candidates: corpus3 delta folding F-T7-* (a v4.1 pass);
F-R9 route (a) second-moment treatment; L7 continuum proof.

---

## ADDENDUM (same session): PHASE N EXECUTED (v4.1 fold; F-R9 route (a); L7 continuum attempt)

User-directed; workflow wf_9558d09c-007; adjudicated in
TIER7_ADJUDICATION.md Phase-N addendum.
- **N1**: corpus3 paper now `01-GUM-Omega-Paper-v4.1-ext.md` (8 ⟦T7⟧
  blocks; additive-only verified; 26/26 payload numbers verbatim);
  manifest/record/campaign-status updated additively.
- **N2 (F-T7-N2)**: honest second-moment rate = R × claimed (spatial
  self-averaging, RSD ≤ 5e-4); extreme hot-site Lorenz (100% of rate
  on ≤ 1e-12 site fraction); g_equiv = √R·g = 8.5e-3–12.3, 6.9–10.1
  orders above the corpus's own P-ν4 line; repricing retreat CLOSED by
  arithmetic (f_req 490–44,000 TeV vs [4,60] MeV). Verdict
  INDETERMINATE-FROM-ARCHIVE: the archive prints no Majoron-mode 0νββ
  bound; **conditional kill filed: any external g_lim < 8e-3 excludes
  the honest rate everywhere**. F-R9 terminal menu: adopt-R (face
  g_lim) or route (c) at line-support cost. S1/S6 untouched.
- **N3 (F-T7-N3)**: the line-local crossing criterion REFUTED in this
  geometry (v_x = A(t) + n_y S(z); S(1/2) = 0 and positive wall
  branch ⇒ no sign certificate ahead of a forward-boosted source);
  T4-W5 residual sharpened: the cutoff is a TRANSPORT phenomenon —
  prove via flow-map. Banked: valid n_y = 0 field statement (axial
  support ⊂ [0,15.87] vs sampled 15.787); 60,000-trajectory covering
  with 0 kill-bin crossings (support edges 5.142/6.017 vs kill bin
  6.4). Third productive coordinator spec artifact of the tier.
Next natural opens: the flow-map proof (T4-W5 terminal shape); g_lim
archaeology (corpus-side); F-R14 (authors' choice); corpus3 remains
the standing offer.

---

## ADDENDUM (same session): PHASE O EXECUTED (flow-map proof; g_lim archaeology; F-R14 matrix)

Workflow wf_c362f0c0-baa; adjudicated in TIER7_ADJUDICATION.md Phase-O
addendum.
- **O1 (F-T7-O1)**: backward-reachability proof executed — 64,000
  kill-bin backward trajectories, zero violations at three resolutions,
  min margin 0.621 (22,068× p95 error), residue flux ≤ 2.8e-7;
  kill-bin zeros now FIELD/FLOW grade (~5,400× sharper than L7b);
  T4-W5's remaining open is only the certified/analytic form.
- **O2 (F-T7-O2)**: g_lim ABSENT-CERTIFIED (658 files, 14 patterns,
  370 hits catalogued; nearest miss triply disqualified); the F-T7-N2
  conditional kill premise certified, armed archive-externally at
  g_lim < 8e-3.
- **O3 (F-T7-O3)**: F-R14 matrix computed exactly — ONE sound cell of
  ten (P-acoustic + supertrace); supertrace alone trades inverted-G
  for hull-exit; tuning fraction 1/4; h33's minimum bill computed.
  Decision remains the authors'.
Remaining opens after Phase O: certified/analytic T4-W5 form;
archive-external g_lim; the authors' F-R14 adoption; corpus3 v4.1
standing offer. A corpus3 v4.2 delta folding F-T7-N*/O* is the natural
next fold if directed.

---

## ADDENDUM (same session): PHASE P EXECUTED (v4.2 fold + T4-W5 certification pilot)

Workflow wf_d3b1ed78-3d0 + finisher agent; adjudicated in
TIER7_ADJUDICATION.md Phase-P addendum.
- **P1**: corpus3 paper now `01-GUM-Omega-Paper-v4.2-ext.md` — the
  Phase-N/O fold (7 new ⟦T7⟧ blocks; additive-only verified 39/39).
  The standing offer to the authors = the v4.2 edition set.
- **P2 (F-T7-P2)**: exact-field certification — grid interpolation
  eliminated (9.9e-15 field agreement); 432/432 classifications kept;
  24-path 200-bit ladder: margins stable ~10 digits; worst margin
  certified to 36 digits (0.620679003382508519636…). G3 FAIL honest:
  endpoint clause only (1 path, 1.70×), margins 4 orders inside —
  sensitivity localized to free endpoints. Survived one process death
  + one container restart via committed phase intermediates.
- **T4-W5 terminal ladder now**: validated ODE enclosures (formal) →
  analytic proof (terminal). All other named opens are corpus-side
  (F-R14 adoption; archive-external g_lim).

---

## ADDENDUM (same session): PHASE Q EXECUTED (formal rung + analytic program)

Workflow wf_74d72548-c07; adjudicated in TIER7_ADJUDICATION.md Phase-Q
addendum.
- **Q1 (F-T7-Q1)**: FORMAL RUNG ATTAINED — 7/7 kill-bin pre-crossing
  certificates machine-checked by validated interval enclosures
  (Lohner QR + Picard boxes over the exact 2048-mode field with
  outward rounding); worst-margin path certified over its full
  3.555-unit excursion (X_x ≥ 1.0010188761, width 7.8e-4, 233k
  validated steps, zero retries); 16,000/16,000 containment; engine
  regularization provably never binds. Trust base one rung below
  CAPD-grade (stated).
- **Q2 (F-T7-Q2)**: Lemma A proven exactly (Fresnel identity + moment
  bounds, explicit constants; 0/610 verification violations);
  Conjecture C with proof obligations O-1…O-5 each carrying a
  sufficient condition; honest lemma-vs-kill-window delta printed.
- **T4-W5 ladder now**: L7b zeros → O1 field/flow → P2 exact-field
  precision → Q1 machine-checked witnesses + Q2 analytic program.
  Remaining: formally-verified kernel + population certification
  (engineering); the analytic theorem via O-1…O-5 (mathematics).

---

## ADDENDUM (same session): PHASE R — THE v4.3 EDITION SET (session fold complete)

Workflow wf_8444d72f-484, 12/12 verification gates PASS.
- Paper → `01-GUM-Omega-Paper-v4.3-ext.md` (P2/Q1/Q2 + the complete
  T4-W5 ladder folded at VIII.F; Sec. 0/X/XI updated).
- Course/Primer/TE received their consolidated Tier-7 addenda (the
  outstanding consistency item — closed; zero deletions).
- Session Map §19, Watch memo §11, campaign status P/Q/R entries,
  SESSION_HANDOFF ADDENDUM 2 (restoration entry points) — all additive,
  85 payload strings verbatim, workflow IDs consistent.
**corpus3 v4.3 = the complete, internally consistent standing offer.**
Session state: all measurement phases (L, L5′, L7, M, N, O, P, Q) and
all folds (T6, v4.1, v4.2, v4.3) executed and adjudicated. Remaining
opens are exactly: the formally-verified kernel + population-scale
certification (engineering), the analytic theorem via O-1…O-5
(mathematics), F-R14 adoption and archive-external g_lim (both
corpus-side).

---

## ADDENDUM (same session): PHASE S — TIER 8 ASSESSMENT + ALTERNATIVES

New directory: `tier8-assessment/`. Two workflows + one nohup chain.

**S-Assess** (wf_f8099700-bc8): the audited assessment pair —
`PROS_ASSESSMENT.md` (14 points: kinematic skeleton genuinely derived
and replicated; closed-form repairability; blind re-derivation
convergence; honest-failure carriage; steelmanned adversarial process)
and `CONS_ASSESSMENT.md` (14 points: centerpiece failed at value level,
survives as replicators' repair; calibration-shielded ħ; gravity/
statistics/0νββ broken or conditional as printed; family flagship
unauditable at core; single self-adjudicating actor system; terminally
undischargeable findings). Audit symmetric: 11 SUPPORTED + 3
OVERSTATED-corrected + 0 UNSUPPORTED on EACH side. Raw briefs/verdicts:
`s_assess_briefs_audits.json`. The pair is the program's standing
summary; neither document may be quoted without the other.

**ROADMAP v10 executed** (wf_871291ab-d94; 16/16 gates):
- **F-T8-S1** (`S1/`): F-R14 sound-cell cost (1) priced — no measure
  weight both rescues F-R14 and preserves Sec. III's [DF] chain
  (cancelling w = 1/2 outside window); w = −1/4 correction
  (9ħ²/16m)∇h·∇ρ/ρ (scheme A). Third quantified IOU of adoption.
- **F-T8-S2** (`S2/`): closed-form G*(λ) for the chirality family
  (midpoint (4+π)/3); attainment generic; corpus identities are λ = 1
  faces of generic BPS identities (I_λ = λI + 2(1−λ)J); chirality
  selects the VALUE, not the STRUCTURE; c(0.81778) = π numerological.
- **F-T8-S3** (`S3/`): forward Lohner covering pilot — 0% certified,
  wrapping-bound (not flow-bound); O-1 needs Taylor-model technology;
  costs priced both ends (1.2e13–9.4e16 boxes this engine vs ~3.5e6
  wrapping-free).
- **F-T8-S4** (`S4/`): N = 10,000 population screen (hybrid RK4/DOP853
  by measured necessity; attempt-1 FAIL preserved); min margin 0.6216 =
  P2 reference + 0.14%; worst-24 all STABLE at 200-bit. Engineering
  open narrows to the formal kernel exactly.

Adjudication: `tier8-assessment/TIER8_ADJUDICATION.md` (incl. the
reflexivity note: CONS C12–C13 apply to Phase S itself; only outside
re-derivation discharges). Archives: MANIFEST workflow10/11 +
tier8-*.tar.xz + both workflow scripts. Campaign status has the ⟦T8⟧
Phase-S entry.

**Restoration note for future sessions:** the Phase-S findings are NOT
yet folded into corpus3 (no v4.4). A future fold phase would draw its
payload from the four RESULTS.md files + TIER8_ADJUDICATION.md; the
assessment pair is campaign-side documentation and does not fold into
corpus editions. Tacit knowledge unchanged (pip3 install numpy scipy
matplotlib sympy mpmath gmpy2 in fresh containers; git from repo root;
nohup + per-item checkpoints for long runs; workflow subagents not
resumable — use fresh finisher agents with full context).

---

## ADDENDUM (session continuation, 2026-08-15): PHASE T — TIER 9, THE GLUEBALL CONFRONTATION

New directory: `tier9-glueball/`. Trigger: user-supplied BESIII Letter
(arXiv:2607.20366, X(2370) = lightest 0⁻⁺ glueball-dominant; read in
full, digest at `GLUEBALL_PAPER_DIGEST.md`).

**Context** (wf_14a465fb-5ae): `GLUEBALL_CONTEXT_ANALYSIS.md` +
`t_context_agents.json` (raw). Core facts for restoration: the corpus's
color sector is VII.J (Q-1 confinement, Q-2 trichotomy, Q-3 dichotomy,
Q-5 emergent-Skyrme with σ = 0.19 GeV² [IM] → f_q 0.14–0.20 GeV
inversion, Q-6′ censorship); glueball-analogs printed exactly once (WS-D
relic census, v1 record); tube = plain-string question was OPEN; closure
obligation (knot-free asymptotic states owe closure) unadjudicated;
seal q-θ forces V.F grading on all hadron-mass numbers.

**ROADMAP v11 executed** (wf_f1635ff3-d44, 16/16 gates, findings
F-T9-T1…T4): T1 both axion-less closed-string routes MISS the
X/lattice-0⁻⁺ window (IP 8.27 high / NG 4.80 low, ratio-failing both);
T2 printed tube = plain bosonic string, no pseudoscalar mode; W_χ =
unique parity-odd class → "tube-core axion" [CJ-new] offer; T3 scorecard
1/5/1, F-Q8 closure obligation NEW-OPEN; T4 two-sector dictionary +
four ⟦T9⟧ annotations drafted (`T4/t4_annotations.md`).

**Restoration note:** Tier-9 findings are NOT folded into corpus3 (no
v4.4). A future fold would draw from the four RESULTS.md + the
annotation drafts. The single highest-leverage next computation: the
tube-core-axion existence/gap derivation (W_χ-sourced parity-odd core
mode) — T1+T2 jointly show the corpus's ability to house an X-class 0⁻⁺
hinges on it. Also open: F-Q8 discharge; the h25 missing mediator
branch; kappa_q valuation. Tacit knowledge: arxiv.org egress BLOCKED
(WebFetch and curl; WebSearch works; user-uploaded PDFs readable);
Phase-S restoration notes still apply.

---

## ADDENDUM (session continuation, 2026-08-16): PHASE U — TIER 10, THE BARYON-JUNCTION CONFRONTATION

New directory: `tier10-junction/`. Trigger: user-supplied STAR paper
(arXiv:2408.15441 = Science 10.1126/science.ads5962, 2026-08-13;
baryon number traced by the gluonic Y-junction, not valence quarks;
digest at `JUNCTION_PAPER_DIGEST.md`).

**Core restoration facts:** the corpus prints B = (1/3)·(signed count
of fractional-residue knot cores) — Theorem K-1, WS-K worksheet stratum
ONLY (operative editions B-silent = fold debt; delta-pack stubs
K-Δ-1..9 unaffixed). Census object 𝕂 = 3B + L. Junction = enforcer via
taping rule; transport SILENT. WS1-H3 (junction order n = 3) absent
from archive — ansatz, [CAL]-at-best, F-K0-1 Primer-vs-paper seam.
Findings F-T10-U1..U3 (12/12 gates): U1 J–J̄ dumbbell spectrum
(convention split ring-analog 2.6–2.9 vs stretch 5.5–5.8 in m/√σ —
priced; interleaving robust; T1 machinery extended, Airy-validated);
U2 junction-led transport strictly cheaper across the whole printed
window, exact ledger proves K-1-compatible (trichotomy: knot-carried,
taping-enforced, junction-transported); U3 ΔB = 0 exact in all four
printed operations, scorecard 4 SILENT/1 V.F/0 TENSION, five ⟦T10⟧
offer annotations (`U3/u3_annotations.md`).

**Restoration note:** Tiers 9–10 findings NOT folded into corpus3 (no
v4.4). A future fold draws from tier9 T1–T4 + tier10 U1–U3 RESULTS +
both annotation files. Highest-leverage next computations: (1) the
tube-core-axion existence/gap derivation (Tier-9); (2) the U1
inertia-convention resolution (derive the dumbbell's true kinetic term
from the substrate action); (3) WS1-H3 discharge via junction-EFT +
lattice junction-mass technology; (4) the K-stub affixation fold.
New [IM] anchors this phase: M_J/√σ = 0.1355(36) (2+1D); STAR
1.84/0.64/1.04 set (contested — CGC-saturation rival at parity).
Tacit knowledge: ALL scholarly domains egress-blocked this session
(WebSearch only; user-uploaded PDFs readable); prior notes apply.

---

## ADDENDUM (session continuation, 2026-08-16): PHASE V — THE CORPUS4 FOLD

`corpus4/` = complete new edition set of all corpus3 documents with EW
(synthesis v4.4-EW/v4.5-SS §8, owner-commissioned adoption) + strong
sector (Tiers 8–10) folded in. Governing docs:
`corpus4/REVISION_CHARTER_v3.md` (8 rules) +
`corpus4/EW_SS_FOLD_PAYLOAD.md` (Parts A–D, payload authority).
Editions: paper v5.0-ext (1,275 lines), Course v4-ext, Primer v4-ext,
TE v4-ext, Session Map v4.0-ext, Record v3.0, manifest v6.0-ext, Watch
memo v4-ext, T6/T7 payloads (3-line headers, bodies byte-identical),
FOLD_RECORD (verdicts filled). Fold workflow wf_1cadd65b-be9; 11/11
coordinator-verified PASS; +622 lines total; deletions header-class
only. Notable editor amendments (all printed in edition notes):
chapter-collision resolutions (Course Ch. 20-W not Ch. 21; Primer
Ch. 15½ not Ch. 16; §15.7 number kept); paper PART-header version
bumps; E1 self-caught SS-Closure-V verbatim defect repaired pre-gates.

**Restoration note:** the session also survived a container restart at
Phase-V start (fresh clone checked out the ORIGINAL branch
claude/analyze-gum-po-5x5cne at a Phase-K commit — fix: git fetch +
checkout claude/analyze-gum-po; reinstall python deps). The synthesis
line lives in `synthesis/` (v4.4-EW preserved + v4.5-SS, published as
artifact 62369470-e9fa-449e-b340-7f35f3fc16e2). Remaining opens
unchanged in kind: SS-Closure-I…V, EW-Closure-I…III (now printed in the
paper's Sec. XI), the Tier-8 engineering/math opens, F-R14 corpus-side.
