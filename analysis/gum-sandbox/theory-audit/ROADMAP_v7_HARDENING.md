# GUM Program — Roadmap v7: Hardening & Frontier (post-Phase-J)

Successor to ROADMAP_v6 (executed: Phases I and J; ledger F-R1–F-R17;
T-ledger empty; no named opens). What remains is HARDENING: the three
micro-opens the J memos printed honestly, plus the strongest guarantee
an external replication package can carry — a fresh end-to-end re-run
of every command it promises.

## Phase K — four workstreams ▶ EXECUTING

**K1 — the /6 vs /(15π/8) discriminator (J1's remaining split).** The
edge-referenced amplitude convention is recovered AS A CLASS; the two
in-class candidates differ by 1.87% in 𝔭 (0.7σ against ±0.03 — archive
error bars alone cannot split them). Find and execute a discriminating
computation: an independent archive quantity sensitive to the amplitude
convention at better than ~2% (candidates: the I.2/G4 two-knot tail
amplitude, the B.4 band edges themselves, the I.1 gate constants). A
verified null — no such quantity exists in the archive — is a
legitimate terminal state and closes the item as "indistinguishable
within the archive"; print it as such. → k1_* artifacts.

**K2 — the ε < 0.005 probe (J4's stated limit).** Extend the J4 scan
to ε ∈ {0.001, 0.002} (two resolutions each, honest convergence gates;
if the solver cannot converge at a point, print the limit rather than
forcing it). Does the linear law hold to the smallest honestly
solvable ε? → k2_* artifacts extending j4_results.json's dataset.

**K3 — GPU pipeline extension (the "Beyond" on-ramp made concrete).**
Port the sector-measure + E_static gradient pipeline to WGSL f64 and
verify tolerance-band against the fs-gum-statics CPU goldens on
llvmpipe — the full readiness demonstration for real-GPU physics
(scheduled after K1/K2/K4 land; the heavy item).

**K4 — the REPRODUCE gauntlet (package freshness).** Execute every
command in DISCHARGE_PACKAGE/REPRODUCE.md fresh on this host, top to
bottom: record PASS/FAIL, runtime, and any drift from the printed
expectations (goldens must reproduce bit-identically within-backend;
tolerance-band items within band). Deliverable:
DISCHARGE_PACKAGE/FRESHNESS.md — the dated re-run record an external
replicator sees first. Any failure is a finding, not an embarrassment:
print it. → k4 artifacts.

## Phase K verdicts (K1/K2 adjudicated; K4 running at fold time)
- **K1: VERIFIED NULL — closed as "indistinguishable within the
  archive"** (17/17 gates; coordinator re-verified; ~2 s determinism).
  The archive prints exactly two amplitude-bearing calibrands; every
  channel falls short of 2σ: the 𝔭 band 0.52σ (NOTE: the 0.7σ quoted
  in J1/ROADMAP_v6/this file's K1 brief OVERSTATED — measured 0.52σ,
  the null is stronger than advertised); 𝔟_eff ≤ 0.37σ and doubly
  voided structurally (F_b relabel absorbs the switch to 0.06σ; the
  C_d ×54 systematic); the bond loop 0.075σ; the B.4 window edges
  re-derive EXACTLY but are 𝒱-functionals with zero amplitude
  sensitivity (the edge signs newly closed: −m̃²/2 anti-vacuum
  positivity, +m̃²/6 the a₄ sign change); C₆, μ, both x₀ measurements,
  the I.1 gate constants all amplitude-free; no field snapshot exists
  in the frozen two-knot records for a 3-D re-fit. One conditional
  pointer printed, NOT claimed: "0.84" rounds from /6's 0.8434 and not
  /(15π/8)'s 0.8591 ("0.86") — a 6.6σ-equivalent split IF three
  unprinted assumptions hold. Defect acted on: the on-disk
  j1_results.json partial stage-flush (regenerated to the full 9
  points in this fold). k1_RESULTS.md.
- **K2: the linear law holds to ε = 0.001 — J4 STRENGTHENED** (11/12;
  coordinator re-verified; both points at unchanged budgets, conv
  3×10⁻⁵; ε = 0.005 reproduction +1.24×10⁻⁶; 2/3 excluded by the new
  window alone, χ² 77; free-intercept p = 1.005 [1.003, 1.007]).
  Headline defect K2-D1 (the one FAIL, printed): the 12-mode family's
  intercept bias DRIFTS (~+1.7×10⁻⁴ at ε = 0.001, +3.15σ vs constant)
  — intercept-shaped by the χ² shape test, so small-ε fixed-intercept
  exponents are bias artifacts; undermines neither exponent, slope,
  nor ceiling. Propagation UNCHANGED (movement 7.8×10⁻³ ≤ the 1.3%
  trigger); J4's (4.9′) envelopes stand. k2_RESULTS.md.
- **K4: THE PACKAGE REPRODUCES — zero physics drift.** All 46
  REPRODUCE rows executed fresh on this host (sole skip: K28's ~24-min
  corrected form, documented); every physics comparison bit-identical —
  7 Merkle roots, 3 fingerprints, every results-JSON physics field,
  every PNG, and the FULL 2500-iteration N=192 crossing series on the
  4-thread solver (window [2070, 2080] identical). All ten findings
  process-level: F1/F5/F7 two printed-command defects (one with a
  clobber hazard; corrected forms verified or documented), F6/F8/F10
  runtime-column drift with F8 ROOT-CAUSED — h23 unpinned hung 5.6 h in
  sympy's Wang/Hensel factorization; PYTHONHASHSEED=0 reproduces the
  primary physics-identically in 184.5 s; the pin-fix and a re-printed
  runtime column are written into FRESHNESS.md — F2/F3 stale test-count
  expectations, F9 the shipped field3d log was itself truncated (fresh
  run completes it). DISCHARGE_PACKAGE/FRESHNESS.md is the dated
  record. Phase K complete; K3 (GPU pipeline extension) is the one
  remaining scheduled workstream.

## Beyond (unscheduled)
The corpus's response to corpus2/DISCHARGE_PACKAGE; real-GPU hardware.
