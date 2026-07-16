# GUM Program — Roadmap v6: Synthesis (post-corpus2)

Successor to ROADMAP_v5 (executed: ledger F-R1–F-R16, T-ledger empty,
corpus2 built). The remaining frontier is SYNTHESIS: does the repaired
model cohere — and does the corpus's own frozen archive already
corroborate the repair?

## Phase I — four workstreams ▶ EXECUTING

**I1 — the ⟨r5⟩ corroboration test (the headline question).** The
corpus's frozen quartic-regime closure ⟨r5⟩ reports 𝔠_q = 3.1 ± 0.2 —
0.5σ from the repaired saturated value 64√2/9π = 3.2011. The halo
threshold κ_ours = 1/√(8π) is t-INDEPENDENT (F-R5's derivation), so the
saturated closure exists at every ε. Compute: (a) the full 𝔠_sat(ε)
curve from ε = 0.05 through the quartic regime (ε ≈ 0.5–2) on the
axi/gstar machinery; (b) the direct (restricted) quartic closure for
comparison with the corpus's method-class; (c) the verdict: does
⟨r5⟩'s 3.1 ± 0.2 match OUR quartic-regime value — i.e., was the
corpus's own ⟨r5⟩ measuring the saturated branch all along? Either
outcome is major: corroboration (the archive contains the repair's
fingerprint) or a pinned quartic-regime discrepancy.

**I2 — repaired-closure phenomenology propagation.** With T3.1 adopted
(as corpus2 does), propagate downstream systematically: the ħ = 𝔠Λ√J
calibration (Λ√J refit; what is OBSERVABLE within-model vs absorbed),
κ_phys = 1/√2 consequences (range 29.3% shorter; the P8-7 bench
number), ω_th = √2ω₀, the S4′ retarget, binding depth, the AUD-15
V15.3–V15.9 phenomenological confrontations re-run under the repaired
constants (which survive, which shift, which improve), P-O1 blindness
verified sector by sector. Deliverable: the propagation memo corpus2's
§IX needs.

**I3 — the joint convention solve (𝔭 × 𝔟_eff).** Two spec-
underdeterminations resisted single-constraint recovery: 𝔭 = 0.84 ±
0.03 (tail amplitude; nearest candidate 0.805) and 𝔟_eff = 42 ± 6
(bond constant; ~7 orders unrecovered). Enumerate the normalization-
convention space programmatically (amplitude definitions × density
weights × 2π/4π factors × per-degree vs total × the App-A candidates)
and solve JOINTLY: a single convention satisfying both plus the
already-recovered quantities (μ, R*, the unit map) would reconstruct
the frozen pipeline's conventions — discharging two open items at
once. A null (no joint solution) pins the incompatibility.

**I4 — the topology completion package.** The three scoped opens from
H4/H3.1 worked to ground: (a) the N2/F.3 exchange lift attempted
constructively on H4's explicit model (the two-line moduli exchange
loop; success ⟹ "ARE fermions" earned at last, failure ⟹ the minimal
hypothesis pinned sharper); (b) curved-line moduli (does ∂ = ±2[rot]
survive the retraction's relaxation?); (c) even w ≥ 2 (the composite
sectors; the structural parity protection tested).

## Phase I: ✅ EXECUTED — the synthesis verdicts
- **I1**: the ⟨r5⟩ test is SADDLE-EXCLUDING and repair-consistent —
  under every ε-stable reading the archive's 3.1 ± 0.2 selects the
  saturated branch (~7σ separation; anchor 64√2/9π at 0.51σ) and
  cannot come from the corpus's printed chain. Two new rigorous
  floors (𝔠_sat ≥ 64√2/9π at every ε; Faddeev-type +12t). Filed as
  archive-side evidence for the repair (normalization caveat printed).
- **I2**: propagation memo — every data confrontation SURVIVES (S1
  verified untouched; bridge ≤0.12σ); the bench sharpens (S4′ now
  discriminating at 13.4σ; the four-rung over-spin ladder all exact);
  one honest cost (the forfeited ε-anchor) + Q-6′ margin ×4.3 thinner.
- **I3**: joint convention solve — ONE edge-referencing element
  resolves both 𝔭 and the entire "7 orders" of 𝔟_eff simultaneously
  (B.4/I.2 consistent); O(1) factor degenerate; the discriminating
  follow-on named (the never-run B.4-window 𝔭-grid).
- **I4 → F-R17**: the exchange lift FAILS with the obstruction
  exhibited (geon phenomenon; F.3's own proviso constructs it; F.4
  vindicated; χ(σ) the new priced ℤ₂); curved/knotted-line moduli
  verdict THEOREM-GRADE (the fatal completion structurally dead
  everywhere); χ(rot) = −1 extends to ALL w ≠ 0 splitness-free
  (even-w fermionic sectors theta-twisted); B=2 sector-conditional.

## Phase J — the filed obligations ✅ EXECUTED (post-fold)

Phase I's fold is in (DISCHARGE_PACKAGE v2.2; corpus2 ⟦I4 reported⟧
boxes). What I1–I4 left open is exactly three named computations, now
scheduled as parallel agents:

**J1 — the B.4-window 𝔭-grid + converged C_d.** I3's discriminator:
App I.1's never-run "𝔭-grid over the (B.4) window", extracted under
both the raw (0.805) and edge-referenced conventions; joint 𝔟_eff
check on the best converged C_d (G5b cap-1200), convergence caveats
printed. → j1_pgrid.py / j1_results.json / j1_RESULTS.md.

**J2 — the repaired (4.8′)/(4.9′) ε-scan arithmetic.** I2's first
obligation: the ε-scan law under the saturated closure (fixed oblate
compacton + halo; B.6 sign flip; the forfeited ε-anchor stated), and
the Q-6′ censorship ceiling re-derivation with the ×4.3 thinned margin
and the dichotomy verdict. → j2_epsscan.py / j2_results.json /
j2_RESULTS.md.

**J3 — the Majoron battery at g ≈ 7.8×10⁻¹⁰.** I2's second
obligation: the corpus's SN/BBN/streaming battery table reconstructed
(old-window control must reproduce its printed PASS), then re-run over
the repaired window f ∈ [3.7, ≈60] MeV — is "SHIFTS-HARMLESSLY"
confirmed, or does a row newly bind? → j3_majoron.py /
j3_results.json / j3_RESULTS.md.

## Phase J: ✅ EXECUTED — the obligation verdicts
- **J1**: the B.4-window 𝔭-grid (first execution; 21 runs, 7/7 gates,
  N/2N/4N ≤ 7×10⁻⁵) **CONFIRMS the I3 edge-referencing element** —
  𝔭 is 35% c₂-dependent across the window; edge-referenced /6 puts the
  corpus's 0.84 at the benchmark potential (0.28σ); the naive /2π only
  reaches it off-benchmark with the referee μ violated 3.6%. 𝔭 upgrades
  OPEN → RECOVERED-AS-CLASS; 𝔟_eff stays in band (+1.06σ) jointly
  (C_d convergence caveat printed). j1_RESULTS.md.
- **J2**: (4.8′) re-derived — 𝔠_sat(ε) = 3.2011·(1 + 1.152ε) leading
  order; κ, depth, ω_th, ¼ all ε-frozen exactly; the ⟨r1⟩ ε-anchor
  forfeit CONFIRMED at 11.283σ (correctly, as restricted-branch code
  validation). (4.9′): ceiling 2.605×10⁻³ — the exact provenance of
  "×4.3"; Q-6′ margins [2.60, 5.31]; the confinement dichotomy SURVIVES
  every reading (worst ×2.03; exponent-robust). All 11 i2 overlaps agree
  ≤10⁻⁴. Two defects of OUR corpus2 revision found and corrected in
  print (the stale "twenty-fold"; the superseded ε^{1/3} shape sentence);
  the multi-ε saturated re-scan named as the residual open. j2_RESULTS.md.
- **J3**: the Majoron battery re-run — I2's "SHIFTS-HARMLESSLY"
  **CONFIRMED** (17/17; old-window control reproduces every corpus
  print; the binding edges are repair-invariant; the moving edge is
  monotonically safer; SN/recoupling rows pass by inheritance — their
  bounds are never printed: spec-recovery limits). The corpus's own
  VII.I dictionary had already stamped ✓ on f ∈ [4, 60]; the re-run
  makes that arithmetic. j3_RESULTS.md.
- **J4 (the Phase-J residual, executed post-fold)**: the multi-ε
  saturated re-scan — dense small-ε solves on the i1 machinery (ε =
  0.005/0.01/0.02/0.03 × two resolutions; ε = 0.02 reproduces i1 to
  −2.7×10⁻⁶; intercept hits 64√2/9π at +7.7×10⁻⁵) measure the leading
  exponent **LINEAR: p = 0.997 [0.996, 1.012]** — the 2/3 class
  EXCLUDED (χ² ratio 8445; admixture 0.0034 ± 0.0014). J2's D3
  discharged in the strong direction: the ×4.3 loosening, 60 MeV f-top
  and 7.8×10⁻¹⁰ Majoron edge are MEASURED (envelopes ×4.31–4.50 /
  59.9–61.1 / 7.66–7.82×10⁻¹⁰); the 2/3 reversion dies with the
  exponent; the dichotomy survives every reading incl. the excluded
  reference (min worst-joint ×1.98 — grazing J2's "factor 2", printed).
  Defects: 1.3% axi4-vs-i1 amplitude spread; upper-bound family bias;
  ε < 0.005 unprobed; W4 curvature drift. 14/14 gates. j4_RESULTS.md.

## Beyond (unscheduled)
Real-GPU physics; the corpus's response to corpus2/DISCHARGE_PACKAGE.
(The former Phase-J residual — the multi-ε saturated re-scan — was
executed as J4: exponent linear, 2/3 excluded, ×4.3 measured.)
