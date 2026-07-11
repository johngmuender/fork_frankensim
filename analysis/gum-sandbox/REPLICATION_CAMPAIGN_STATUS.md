# GUM Replication Campaign — Execution Status

Executing the tiered roadmap of `GUM_SIMULATION_ASSESSMENT_v2_FABLE.md`.
Method: replication-first — independently re-derive the corpus's own
computational claims on an independent codebase (frankensim's certified
numerics), under its own epistemic rules. Every result is within-model;
nothing here bears on nature.

## Tier 0 — Constants Gauntlet ✅ COMPLETE
**`tier0-gauntlet/`** (Rust, frankensim crates; std-only closure).
- **37/37 checks PASS** with fs-ivl certified interval enclosures: closure
  algebra (7/4 exact; 𝔠₀ two-routes + exact square identity; superradiant
  rung certified >1; spin-selection family incl. j=1 pole, j=3/2
  impossibility), quadratures ("two computations one number": Haar ∩
  compacton = 64/15π; hedgehog degree K=1; 8/35), the ⟨r10⟩ bridge
  (m₃ = 0.0468 eV), AUD-15 V15.3–V15.9 phenomenological arithmetic, the
  g_s = 2 → μ_B identity, and the Unruh recovery.
- **Certification loop closed:** 37 `Certified<f64>` → EvidencePackage with
  content-addressed BLAKE3 certificates → fs-checker three ways (deny-all
  refuses unauthenticated Verified claims; recompute-and-compare capability
  passes; tampered root fails). **Bit-identical replay**; golden Merkle root
  `e7a8e897c55984d1a11eaea70ffd7af02a4d528a9d45150d53b3602c21c9df91`.
- **Kill-content: no arithmetic defect found in the corpus.** Two defects
  found were the harness's own (quote-precision test bug; provenance-charset
  schema violation caught by fs-checker) — filed in `RESULTS.md`.

## Tier 1 — Linear Spectrum ✅ COMPLETE
**`tier1-spectrum/`** (Python pilot + coordinator adjudication).
- Four-branch structure replicated exactly: B1 acoustic (c_L to 1e-16),
  B2 gapless light doublet, B3 Klein–Gordon with **ω₀² = (m_V² + 4μ_c)/J and
  the factor 4 measured at 4.000000000**, B4 twist.
- **Theorem II.2 (exact masslessness) replicated:** B2 gap invariant under
  m_V to 8×10⁻¹⁵. Tree achirality exact in the IR (splitting ∝ χ₃k⁵).
- **The M-1 defect exhibited with a sharper mechanism than the corpus's own
  telling:** no gap opens on the gapless branch (translation Goldstone);
  instead the wave's rotational (EM) content dies — co-motion fidelity
  r = 4μ_c/(4μ_c + m_V²), measured = analytic to 9 digits — and the
  rotationally-dominant branch is gapped. Either reading is fatal, as the
  corpus says, by a more precise route (`ADJUDICATION.md`).
- Two literal spec FAILs adjudicated as over-strict tests (mine), with the
  deviations derived in closed form. **No corpus defect found.**

## Tier 3 — Born-Rule Relaxation ✅ COMPLETE
**`tier3-born/`** (Python pilot + adjudication). 2-D box, exact ψ evolution,
N = 20,000 de Broglie trajectories to t = 4π.
- **Relaxation to Born equilibrium replicated**: τ = 83 / 10.2 / 4.9 for
  M = 4 / 9 / 16; M = 16 decays near-exponentially (r² = 0.99) to the
  finite-N noise floor within one revival period; the particle distribution
  does NOT revive when ψ does.
- **Equivariance control verified**: an equilibrium-born ensemble stays flat
  at the noise floor for the whole run — the H-theorem fixed point.
- **Mild tension adjudicated**: M = 4 is pocketed/non-exponential (r² 0.43) —
  attributable to the first-2×2 mode set's degenerate, commensurate spectrum
  (the corpus's own "near-integrable small-M" caveat case); generic-set M = 4
  replication filed as follow-up (`tier3-born/REPORT.md`).

## Corpus survey (supporting) ✅ COMPLETE
The full 104-file archive (`corpus/`) contains **no ⟨r1⟩–⟨r11⟩ pipeline
code**; ⟨r10⟩ is the only release whose arithmetic appears end-to-end
(WS-nu-P2 — replicated in Tier 0/D5). The archive's own Watch-Mode memo
designates the open NR-* gates (NR-A2b ¾-Casimir, NR-D1b core-energy floor,
NR-D2 murk kill, NR-K2b/K3a) as what "external hands must close" —
candidate Tier-2+ targets alongside the ⟨r1⟩ axisymmetric closure.

## Tier 2a — Reduced Closure Loop ✅ COMPLETE (2b remains)
**`tier2-closure/`** (Python + adjudication). One constant calibrated from
one datum (λ* = 0.42); everything else predicted.
- **Coarse structure replicates**: ε→0 endpoint to 6 decimals; ¼ invariant
  automatic (Cor. IV.2 confirmed); g* −0.6σ and 𝔠 +0.9σ inside the ⟨r1⟩
  bands; the ε^{2/3}-class scan law from pure spheroid geometry.
- **Fine structure does not**: κ +2.4σ, V +4.2σ (model pushes V above √2;
  benchmark sits below), c_g ~2× low — a quantified statement that the ⟨r1⟩
  tuple's fine structure requires the field-level (V-coupled, binding)
  ε-sector. Tier 2b's falsification targets sharpened to five numbers:
  κ = 0.802, V = 1.409, c_g = 0.42, δ = 0.29, onset κ = 1.000.
- **FINDING F-R2 (minor)**: AUD-15 V15.2's quoted g(0.37–0.47) span
  [1.28, 1.32] is inconsistent with its own corrected formula
  (true: [1.257, 1.315]); conclusion unaffected.

## Tier 2b Step 1 — Radial Field Solve ✅ COMPLETE
**`tier2-closure/radial_*`** — the campaign's first field-level PDE solve:
the 1-D hedgehog at ε = 0.05 against the corpus's own App. I.1 gates.
- **All gates PASS**: degree 1.0000001; Derrick virial 6×10⁻⁷ (500× inside
  the 3×10⁻⁴ gate); ε dial 0.049985; near-BPS deficit +0.48%; the tail
  form + mass replicate with the paper's B.5 formula μ² = m̃²/(2a_ψ)
  **blind-re-derived and confirmed to 2×10⁻⁶** (against my own spec error —
  the corpus vindicated).
- **𝔭 = 0.84 amplitude: CONVENTION-LIMITED** — nearest natural candidate
  0.805 (~1.2σ low); the frozen pipeline's normalization is not recoverable
  from the text (the predicted spec-underdetermination, now concrete).

## Next (not yet executed)
- **Tier 2b Step 2 (the kill-shot proper):** the axisymmetric isorotating
  2-D field solve
  with the real ε-sector (log-divergent compacton-edge boundary layer
  resolved), targeting κ = 0.802 ± 0.018, V = 1.409 ± 0.010, c_g = 0.42 ±
  0.04, δ = 0.29 ± 0.05, onset κ = 1.000 ± 0.004.
- ~~WS-A4 demand tables~~ ✅ closed (gauntlet Group G: all five rows from
  one δ = 1.0×10⁻¹² window; the "thin corner ∼40" confirmed at 38.6).
- ~~NR-D1b core-energy floor~~ ✅ closed (gauntlet Group H: W-eternal
  replicated; the corpus's G→8 rounding is conservative-direction).
- Remaining: generic-set M = 4 Born rerun (minor).

**Gauntlet now at 50/50** — golden Merkle root
`5640e4fe1e2e52e112af0fac1e7420a15894c384744210be492bbdd50d926b1e`,
bit-deterministic. Every arithmetic-recomputable gate the archive's
Watch-Mode memo left open is now externally recomputed; what remains open
is exactly the *argument-level* physics (NR-A2b's ¾ operator-law selection,
NR-D2's SDiff-reachability, NR-D1b's c_h/r_min inputs) and the field-level
Tier 2b solve.
