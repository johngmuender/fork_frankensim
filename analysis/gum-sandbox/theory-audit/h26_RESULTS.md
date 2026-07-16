# H2.6 — the j(j+1)-convention NON-RIGID closure: the T2 §a.3 GAP, measured

**Phase:** H2.6 (ROADMAP_v5_THEORY item 6; numeric follow-on of
`T2_closure_pillars.md` §a.3).  **Date:** 2026-07-16.
**Engines:** `fs-gum-field` (measurement path, radial solver),
`fs-gum-statics` (guarded ANF protocol), `fs-gum-kern` (deterministic tiled
threaded layer, bit-identical at any thread count) — none modified.
**Solver:** `h26_solve/` (this directory), bins `h26_convention`
(descent tier + H2.7a displacements) and `h26_dilation` (exact-family tier).
**Compute:** 14.1 min (descent tier, 4 threads) + 16 s (exact tier).
**Outputs:** `h26_results.json` (descent tier + rigid level),
`h26_dilation.json` (exact-family tier), `h27a_field_exponents.json`
(the H2.7(a) block, consumed by `h27_entrain.py`).
**Determinism:** frozen LCG perturbation stream; fs-gum-kern bit contract;
in-process bitwise replay of the final corpus-anchor relaxation asserted
PASS; same binary + args ⇒ bit-identical JSON.

**Epistemic frame (binding).** Within-model numerical engineering on a
speculative theory's functional.  Everything below is a measured fact about
the campaign's certified functional and the corpus's argument structure —
never about nature; corpus-level adjudication is the coordinator's.

---

## 1. The question and the modified clock condition (derived once)

T2 §a.3 (GAP): the corpus anchors spin semiclassically, **L = j𝔠**, where
its own imported rotor literature uses **L² = j(j+1)𝔠²**; the j = ½
*selection* is convention-robust, but the closure *tuple* is not — and the
paper never argues the choice.  Writing a = j (semiclassical) or
a = √(j(j+1)) (quantum rotor), the T-B1 family relations w·a(1−a) = 1,
V² = 1/(1−a), E_rot/E = a/2 give at j = ½ the rigid-level predictions
(regenerated here, the `t2_convention.py` lines, all reproduced):

| convention | a | w | V | E_rot/E | 𝔠(g=1) | 𝔠(g=3/2) |
|---|---|---|---|---|---|---|
| L = j𝔠 | 0.5 | 4 | 1.414214 | 0.250000 | 2.053288 | 2.514754 |
| L² = j(j+1)𝔠² | 0.8660254 | 8.618802 | 2.732051 | 0.433013 | 3.014000 | 3.691381 |

The solver-form clock condition under anchor a, derived from the same two
closure premises only (L = a𝔠, E_tot = 𝔠ω, ω = L/𝕀):
E_tot = L²/(a𝕀) ⟹ E_static = L²(2−a)/(2a𝕀) ⟹

**L² = [2a/(2−a)] · I · E_static,  E_rot/E_tot = a/2 identically at any root.**

a = ½ recovers the campaign's frozen Section-K form L² = (2/3)·I·E_static;
a = √3/2 gives the coefficient 2√3/(4−√3) = **1.527415881581** (amusing
near-coincidence, checked and NOT exact: √(7/3) = 1.527525, rel 7.2×10⁻⁵
away — unrelated to the ω_th identity of T2 §d.2).  Because a/2 is
kinematically forced *at* a root, the field-level discriminating content is
(i) existence/character of the root, (ii) the rest of the tuple (κ against
the halo and radiation thresholds, the paper-unit 𝔠-analogue, the inertia
dilation V), (iii) whether the solved E_rot/E stay separated at the ⟨r1⟩
precision (0.2500 ± 0.0002).

**Protocol (both tiers share the anchor state):** N = 48, LBOX = 4.5,
ε = 0.05 (t = 0.0082764349), corner scheme, frozen guards (MU = 5000 at the
engine hedgehog degree, MU_F = 400 at hedgehog gap − 0.01); static stage =
hedgehog + bump(424242, K = 6, amp 0.02), 900 guarded ANF iterations —
identical to the certified gates G-B stage.  Unit maps (SESSION_HANDOFF §3):
𝔠_paper = (L/a)/√(2π³); κ_ours = L/I, κ_paper = 2√π·κ_ours (halo threshold
κ_paper = 1/√2; corpus measured radiation onset κ_paper = 1.000 ± 0.004).

## 2. Tier 2 — the exact dilation-family solve (budget-independent)

Under x → sx the sector functionals scale exactly (E2 → sE2, E4 → E4/s,
E6 → E6/s³, E0 → s³E0, I → s³I; V ≡ s³ is App G.5's dilation dial with the
ε-suppressed E2/E4 corrections the pure corpus family drops), so anchoring
the family at the relaxed static solution (E2 = 8.183790, E4 = 5.814824,
E6 = 1.620982, E0 = 1.249693, I = 22.60711, E_stat = 2.9865333) makes the
constrained stationary family **closed form**: envelope s\*(L) and the
Section-K bisection both run to machine precision (envelope-theorem check
|dE/dL / ω − 1| ≤ 1.8×10⁻¹⁰).  Static optimum s₀ = 1.043310 (V₀ = 1.135637)
is the V-reference.  Solved tuples:

| anchor | L\* | **E_rot/E** | V/V₀ (rigid pred.) | κ_ours (/thr) | **κ_paper** | **𝔠_paper** |
|---|---|---|---|---|---|---|
| corpus L = j𝔠 | 8.789304 | **0.250000000000** | 1.432772 (√2 = 1.414214, +1.3%) | 0.238942 (1.198×) | 0.847028 | 2.232260 |
| alt L² = j(j+1)𝔠² | 22.936161 | **0.433012701892** | 2.859953 (2.732051, +4.7%) | 0.312376 (1.566×) | **1.107343** | 3.363181 |

The V/V₀ values land on the rigid-level convention predictions with the
expected ε-class corrections (+1.3% / +4.7%) — the T2 §a.3 sensitivity
table is reproduced at field-anchored level.

## 3. Tier 3 — the descent-relaxed solve (the campaign's 4A machinery)

At each candidate L the field relaxes under the guarded Routhian (cold
restart from the common static state, 350 iterations per bisection
evaluation, 700 at reporting points; scanned bracket + 16-step bisection +
3-step fixed-point polish):

| anchor | L\* | E_rot/E (resid) | I\*/I_static | halo | κ_ours (/thr) | κ_paper | 𝔠_paper |
|---|---|---|---|---|---|---|---|
| corpus | 10.285957 | 0.249450 (−2.9×10⁻³) | 2.305 | 0.0026 | 0.197394 (0.990×) | **0.699745** | 2.612371 |
| alternative | 17.220874 | 0.430771 (−9.2×10⁻³) | 2.777 | 0.0037 | 0.274269 (1.375×) | 0.972257 | 2.525136 |

**Measured structural fact (the F-R5 mechanism transplanted to the clock):
the descent tier has no budget-independent root.**  At the clock condition
the state is over-spun (κ above the halo threshold at every rigid-level
anchor), i.e. exactly the F-R5 halo-unstable regime: I(L) keeps growing
along the slow halo/dilation flow, so the self-consistent root drifts
upward with iteration budget (corpus anchor: L\* ≈ 6.9 at 60-iteration
evaluations, 10.29 at the 350/700 protocol; adjacent-evaluation F-noise
gives a ±2%-class root band, visible in the logged bisection tail).  The
drift has a *direction*: the corpus-anchor endpoint sits at
κ_ours/threshold = 0.990, κ_paper = 0.6997 ≈ **1/√2** — the descent walks
the clock-charged state toward the campaign's halo-saturation locus
(κ_paper = 1/√2 exactly, F-R5/G\*), while E_rot/E stays pinned at a/2
throughout.  The alternative anchor, needing 2.3× the corpus L², lags the
same flow (κ_paper 0.972, still shedding at budget).  This is the
within-model reason the tuple side of the closure is owned by F-R4/F-R5
regardless of convention — and it is measured here, not asserted.

## 4. The convention discriminator — verdict: SEPARATED

- Exact-family tier: E_rot/E = 0.250000000000 vs 0.433012701892 —
  separation 0.183013 = **915× the ⟨r1⟩ quoted σ** (0.2500 ± 0.0002).
- Descent tier (all budgets): 0.249450 vs 0.430771 — separation 0.181321 =
  **907σ**.  The a/2 identity holds at every budget; the two anchors NEVER
  converge — the halo flow moves κ, V, 𝔠, but cannot move E_rot/E off a/2.

**The corpus's ⟨r1⟩ measurement E_rot/E = 0.2500 ± 0.0002 therefore SELECTS
the semiclassical anchor L = jħ at ~900σ** — a convention-selection
*observable*, closing T2's §a.3 gap in the corpus's favor at the level of
internal consistency: the code and the theorem use the same convention, and
the measurement discriminates it.  (What this does NOT repair: the paper
still never *argues* the choice — the gap remains a one-sentence
documentation gap, no longer a live ambiguity.)

Two supplementary discriminators, same direction:

1. **Radiative admissibility.**  The alternative-anchor closure state sits
   at κ_paper = 1.107 (exact tier) — above the corpus's own measured
   radiation onset 1.000 ± 0.004 (27σ) and above every locked-channel
   threshold except the polar extreme (1.118), which the corpus's control
   run itself excludes (T2 §d.2).  Under the j(j+1) anchor the closure
   would place every elementary fermion in the superradiant regime: the
   corpus's radiative-self-consistency sector is only coherent under its
   semiclassical anchor.  (The corpus anchor's 0.847 is sub-onset, in-gap
   by the corpus's convention; its own tension with the ring limit 1/√2 is
   the known F-R5 charge, not new.)
2. **V-anchor.**  Exact-tier V/V₀ = 1.433 vs 2.860: the ⟨r1⟩ measured
   V = 1.409 ± 0.010 likewise selects the semiclassical anchor
   (the alternative is ~145σ away).

## 5. Verdict table

| item | result |
|---|---|
| modified clock condition under a = √(j(j+1)) | **L² = 1.527416·I·E_static** (derived §1, implemented, solved) |
| rigid-level convention table (T2 §a.3) | regenerated exactly (§1) |
| non-rigid solve, both anchors, exact constrained family | roots exist; tuples §2; E_rot/E = a/2 to machine precision |
| non-rigid solve, both anchors, guarded descent | roots exist per budget but drift (F-R5 halo flow, measured §3); E_rot/E pinned at a/2 |
| discriminator | **SEPARATED — 907–915× the ⟨r1⟩ σ; the ¼ measurement selects L = jħ** |
| alternative-anchor radiative status | κ_paper = 1.107 > onset 1.000 ± 0.004: superradiant under the corpus's own criterion |
| determinism | bitwise replay PASS; kern bit contract; frozen seeds |

**Promotion note for the coordinator:** T2 §a.3's GAP verdict
("selection robust, anchors convention-pinned, choice unargued") is
numerically confirmed and *sharpened*: the anchor choice is now an
internally-measured selection (¼ and V both discriminate at huge
significance, and the alternative anchor is radiatively inadmissible
within-model).  The residual charge is documentation-grade.

## Files

- `h26_solve/Cargo.toml`, `h26_solve/src/bin/h26_convention.rs` — descent
  tier + H2.7(a) displacement runs (usage: `h26_convention <outdir>
  [maxit_l] [maxit_final] [bisect_steps] [threads]`, run at 350 700 16 4).
- `h26_solve/src/bin/h26_dilation.rs` — exact dilation-family tier
  (usage: `h26_dilation <outdir> [threads]`).
- `h26_results.json` — rigid level, descent tier (full evaluation log),
  discriminator block, H2.7(a) descent-tier exponents, replay flag.
- `h26_dilation.json` — exact-family anchors + H2.7(a) family exponents.
- `h27a_field_exponents.json` — the H2.7(a) block for `h27_entrain.py`.
