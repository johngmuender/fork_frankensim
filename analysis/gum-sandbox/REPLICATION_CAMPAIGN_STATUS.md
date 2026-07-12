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

## Tier 2b Step 2 — Spheroidal Isorotating Closure ✅ COMPLETE → FINDING F-R4
**`tier2-closure/axi_*`** — zero-fit field-level closure over two deformation
families, unit map validated to 7–9 digits (g(0.42) matches the paper's
formula to 9 digits; the scaling rung 2√(ê₀𝔦₀) reproduced to 7).
- **F-R4 (substantive):** under honest SDiff pullback the isorotation
  inertia is invariant to 10⁻¹⁵ — the paper's G.5 mechanism (SDiff shape
  modes supplying g → 3/2 at zero energy cost, the basis of the deep-BPS
  endpoint 𝔠₀ = 2.515) **fails in its literal reading**. The g-enhancing
  family is non-SDiff and not energy-flat; restricted-family closures land
  at 𝔠 = 2.05–2.28, κ = 0.90–0.94 — the ⟨r1⟩ tuple NOT reproduced (κ +5.8σ)
  though 𝔠 is within ~1σ. ¼ and V(ε→0) = √2 hold exactly.
- Not a kill (families are restricted; the corpus's unrestricted solve may
  find cheaper non-SDiff modes) — but the deep-BPS endpoint now carries a
  named, machine-precision-grounded doubt whose adjudicator is the
  unrestricted 2-D solve.

## Tier 2b Step 3 — Enlarged-Basis Closure ✅ COMPLETE → F-R4 HARDENED + FINDING F-R5
**`tier2-closure/axi3_*`** — the F-R4 adjudicator: direction-tilt Θ modes
added to the axisymmetric basis (the only modes that can turn the g-dial,
per the F-R4 lemma), gates G0–G5 all PASS, Step-2's closure reproduced to
8×10⁻⁶ before enlargement.
- **The answer to F-R4's question is "neither"**: tilt modes do not carry
  the ε→0 closure from the scaling rung (2.0533) to the paper's endpoint
  (2.5147) — they drive it **past** 2.5147 without stopping. Direct solves
  hit their caps still descending (𝔠 = 2.78 and climbing); the well-posed
  saturated limit gives **𝔠_paper ∈ [2.828 (rigorous bound), 3.201]**.
- **FINDING F-R5 (substantive):** the fixed-L Routhian is variationally
  unstable to a far-field tilt halo whenever κ_paper > 1/√2 — the standard
  over-spin criterion ω > μ(tilt), derived (dR = [∫η²](1/8π − κ_ours²)) and
  verified by shell probes (κ_crit 0.216 measured vs 0.199 ideal, dI
  quadratic). **Every solution the corpus quotes (rung 0.935, benchmark
  0.802, deep-BPS 0.764) sits above the threshold**; the true infimum is
  halo-saturated at κ_paper = 1/√2 exactly. The ⟨r1⟩ benchmark is a
  **halo-unstable saddle** of its own functional; 𝔠₀ = 2.5147 is excluded
  from below by a rigorous bound in the enlarged space. The corpus's
  "onset at κ = 1.000" claim is normalized inconsistently with its
  closure-algebra κ's (a pinned √2; both readings fail within-model) —
  working hypothesis: the frozen ⟨r1⟩ code direction-locks the hedgehog,
  an unstated restriction. Full analysis: `axi3_ADJUDICATION.md`;
  measured facts: `axi3_RESULTS.md`.
- Untouched throughout: E_rot/E = ¼ (10⁻⁹), clock residuals ≤ 4×10⁻⁹, the
  spin-selection and rung algebra — the *selection* theorems survive; the
  *value* 𝔠₀ does not, within the stated variational problem.

## Findings ledger (the campaign's product)
- **F-R1** (Tier 0/F): NR-D2's printed KE coefficient ×25 too large at its
  stated v — conservative direction; the murk kill fires harder.
- **F-R2** (Tier 2a): AUD-15 V15.2's g-span [1.28,1.32] inconsistent with
  its own corrected formula (true [1.257,1.315]) — harmless.
- **F-R3** (Tier 0/I): ⟨r10⟩'s ±0.90 band implies ρ_eff = −0.36, not the
  printed +0.45 — band-arithmetic slip; the S1 stake unaffected.
- **F-R4** (Step 2, hardened by Step 3): G.5's SDiff inertia-enhancement
  mechanism fails at machine precision — inertia is invariant on the SDiff
  orbit (proof: `axi_FR4_NOTE.md`); no wording repair survives Step 3.
- **F-R5** (Step 3): the closure's stated variational problem does not
  select 𝔠₀ = 2.5147 or the ⟨r1⟩ tuple — tilt-halo instability above
  κ_paper = 1/√2, saturated infimum ≥ 2√2. Discharge path: the corpus
  produces the frozen ⟨r1⟩ code exhibiting the stabilizing
  constraint/term absent from App. I.4's printed protocol.

Within-model consequence if F-R4/F-R5 stand: every constant downstream of
𝔠₀ (the ħ calibration 𝔠Λ√J, κ_phys = √(7/12), the S4′ ε-run anchor,
binding depths) shifts by −18% (rung, with a stability constraint) to
+12…+27% (saturated) against the printed values. The replicated successes
stand alongside: the four-branch spectrum, M-1, the ¼ and √2 invariants,
spin selection, the radial rung (tail mass blind-confirmed), Born-rule
relaxation incl. equivariance and the generic-M ≳ 4 onset, ⟨r10⟩'s
central value, and all archive arithmetic (52/52).

## Closed side-gates
- ~~WS-A4 demand tables~~ ✅ (gauntlet Group G: all five rows from one
  δ = 1.0×10⁻¹² window; the "thin corner ∼40" confirmed at 38.6).
- ~~NR-D1b core-energy floor~~ ✅ (gauntlet Group H: W-eternal replicated;
  the corpus's G→8 rounding is conservative-direction).
- ~~Generic-set M = 4 Born rerun~~ ✅ (τ = 4.48, r² = 0.987; the corpus's
  "M ≳ 4" holds for generic sets — Tier-3 tension resolved in the
  corpus's favor).

**Gauntlet at 52/52** — golden Merkle root
`e19d0bcb3e97e430d9cca8d0823705f7b6b0807bf3601d5f40923c5f010a85f3`,
bit-deterministic. Every arithmetic-recomputable gate the archive's
Watch-Mode memo left open is externally recomputed; what remains open is
the *argument-level* physics (NR-A2b's ¾ operator-law selection, NR-D2's
SDiff-reachability, NR-D1b's c_h/r_min inputs) and, now, the corpus's
response to F-R4/F-R5.

## Campaign status: COMPLETE (Tiers 0–3, Steps 1–3, archive gates)
Possible extensions (not scheduled): a fully unrestricted 2-D
Newton–Krylov solve reproducing the corpus's saddle under explicit
direction-locking (to test the F-R5 working hypothesis); the fs-cosserat
certified-Rust port of the closure pipeline from the feasibility plan.
