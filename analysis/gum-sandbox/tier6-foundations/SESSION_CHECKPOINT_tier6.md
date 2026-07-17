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
