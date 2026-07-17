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
- **F-R5** (Step 3; confirmed in 3-D by Tier 4A, sharpened by 4C): the
  closure's stated variational problem does not select 𝔠₀ = 2.5147 or
  the ⟨r1⟩ tuple — halo instability above κ_paper = 1/√2 (tilt channel;
  locked profile channels between √3/2 and 1/√2), saturated infimum
  ≥ 2√2, clock-forced. Direction-locking does not rescue it. Discharge
  path: the corpus produces the frozen ⟨r1⟩ code exhibiting the
  stabilizing constraint/term absent from App. I.4's printed protocol.

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

## Tier 4 — 3-D Texture Field ✅ COMPLETE (three workstreams)
**`tier4-field/` + `tier2-closure/axi4_locked_*`** — adjudication:
`tier4-field/TIER4_ADJUDICATION.md`.
- **4A (unrestricted 3-D pilot)**: F-R5's mechanism **confirmed in full
  3-D with no symmetry assumption** — fixed-L descent fell 0.56 below the
  axisymmetric stationary value, κ through the threshold, halo box-limited;
  88% of the inertia gain priced in E₀ at exactly the threshold rate
  1/(16π); the below-threshold control self-limited at the (3/2)·I tilt
  ceiling (two-sided verification). **The clock condition forces the
  over-spun regime** (κ(L_clock) = 1.25–1.47× threshold in every reading).
- **4B (fs-cosserat-pilot)**: 13/13 certified 3-D diagnostics, golden root
  `f88731af…`, bit-identical replay — degree O(h²), sector energies vs
  fs-ivl enclosures two ways, F-R4's SDiff invariance and the ¼ invariant
  as certified claims.
- **4C (direction-locked closure)**: the Step-3 working hypothesis
  (direction-locking rescues the benchmark) **refuted** — locked profile
  halos have channel-dependent thresholds κ_crit = 1/√(2⟨sin²θ⟩_w)
  sliding from √3/2 (uniform) to the unrestricted 1/√2 (equatorial ring);
  the converged locked closure runs away (𝔠 = 3.26, V = 7.1 at ε = 0.05).
  **Two exact identities found**: the corpus's measured over-spin onset
  1.000 ± 0.004 = the polar-channel threshold (exactly 1, ⟨sin²θ⟩ = 1/2);
  its deep-BPS κ₀ = √(7/12) = the sin²-channel threshold (exactly 6/7
  moment identity); plus the locked uniform-saturated closure landing at
  𝔠_paper = 2.37096 vs the benchmark 2.37 ± 0.09. Reading (hypothesis,
  flagged): **the ⟨r1⟩ numbers are thresholds/saturation values of the
  halo family the corpus's solver was blind to.**

## Tier 5a — Family/Neutrino/Color Arithmetic ✅ COMPLETE
**`tier5-family/`** — the arithmetic perimeter of the "no fitted
parameters" sector (the Σ(p) integrals themselves are not recoverable
from the text; everything downstream is audited).
- **Replicates cleanly**: lepton data logs and both spacing pulls (0.05σ,
  0.43σ); shape ratios; ΔN_eff = 0.026772; termination arithmetic; the
  soft-sector window (the "Cabibbo" symbol collision is the corpus's own
  declared and sealed non-issue).
- **⟨r11⟩ quark belts VINDICATED by full reconstruction**: the rule
  R = 1 + B_geo/A − B_tube/A reproduces the printed ratios to the error
  digits, the down-type sign flip is forced, PDG data lands exactly with
  direct m_t = 172.5 GeV (MS-bar fails — scheme pinned), four pulls at
  0.17–0.32σ as claimed.
- **FINDING F-R6 (moderate)**: the ⟨r7⟩→E-F1→⟨r8⟩ narrative chain is
  arithmetically false as written (2.9 ×2 +14% = 6.61 ≠ 5.6; printed
  errors shrink where the chain forces growth). Coherent reading: ⟨r8⟩
  is an independent re-evaluation and the provenance sentence is wrong.
  Direction-neutral (the honest chain would worsen data agreement).
- **F-R7 (minor)**: K.2's "qR* ∼ 10⁻¹⁷" is six orders from its own
  inversion (2.6×10⁻¹¹); inversion side coherent; S1 unaffected.
- **F-R8 (minor)**: locked-bond pull prints 1.3σ where its own errors
  give 1.54–1.59σ.

## Campaign status: COMPLETE (Tiers 0–5a, Steps 1–3, archive gates)
- **F-R9 (Phase H2.1; substantive for the 0νββ insert)**: promoted from
  theory-audit T-H2 — NR-K3a's "fluctuations negligible" clause fails
  by 13–18 orders (rate is quadratic in the amplitude; line-supported
  web gives R = ⟨A²⟩/⟨A⟩² = 10¹⁴–10¹⁸ across the corpus's own scale
  window, MC + Campbell-exact; the MEAN claim simultaneously verified
  at 1.000 ± 0.001 — the defect is purely the dropped second moment).
  Motional-narrowing rescue priced at 5×10⁷–1.2×10⁹ c — causally
  excluded in-model; both V-K1 fork branches unsupported. Discharge:
  a second-moment treatment, an on-tube saturation nonlinearity
  (nowhere in the text), or a smooth bulk m-component (which costs
  Step 1's line-support theorem). S1 (Σm_ν) untouched; the K3a 0νββ
  rate phenomenology is what's affected. theory-audit/h21_RESULTS.md.

- **F-R10 (moderate; Phase H2.4 + T-audit)**: NR-A2b as printed — the
  ¾-Casimir normalization mixing and the null-transport exactness's
  explicit counterexample classes. theory-audit/h24_RESULTS.md.
- **F-R11 (minor; Phase H2.3)**: the D1b floor overstated 8% (honest
  relaxed floor 7.7024); Frank-convention resolution printed with our
  own mis-audit correction. theory-audit/h23_RESULTS.md.
- **F-R12 (moderate; Phase H T4)**: VI.D — the exact tracker identity
  w(z) = −1 + Ω_m(z) contradicts the printed direction-check sentence.
  theory-audit/T4_gravity_oneworld.md.

- **F-R13 (moderate; Phase H3.3-A)**: Lemma II.1/App F.1's convexity
  proof is FALSE by strict counterexample (det F_t < 0 on an open
  interval; the admissible set is not even star-shaped) and the
  admissible class is ill-posed; the CONCLUSION survives under the
  supplied repair Lemma II.1′ (degree locally constant along
  det F ≥ δ histories) — pairs-only and the IBC construction stand.
  theory-audit/h33_RESULTS.md.
- **F-R14 (substantive; Phase H3.3-B, promoted from T-H10)**: the
  corpus's OWN action pins the E.6 response exponents and closes the
  positivity escape wedge — the photon doublet has p = q − 1 under
  every index-variance reading and density weight, forcing
  𝔞₁^(B2) < 0 in every textually permitted reading; Theorem VI.1's
  convexity core fails UNCONDITIONALLY as printed (favored reading →
  inverted Newton constant vs VI.B; covector reading → mixed signs,
  c_GW² can exit the hull). Minimal rescue quantified ("P-acoustic"
  measure, weight window 1/18 wide) with its priced costs against the
  Sec III quantum chain and F5/KBKK. theory-audit/h33_RESULTS.md.

- **F-R15 (moderate, F-R5-corroborating; Phase H3.2)**: IV.D's
  selection narrative fails as a dynamical claim — the honest
  quasi-static flow from the rigid rung crosses κ = 1 at 4% of the
  flow with no distinguished feature, crosses the benchmark 0.802
  without pause at cycle 130, and converges to the saturation-class
  locus (attractors at/below 1/√2) exactly as Theorem T3.1 predicts;
  even the tightest compact-support rescue converges 5.5σ below the
  benchmark. NEW measured fact: the (κ, 𝔠) flow passes within the
  ⟨r1⟩ tuple's own error bars at cycles 130–170 — the benchmark reads
  as a budget-stamped snapshot of the transit. The superradiant
  REJECTION itself survives untouched. theory-audit/h32_RESULTS.md.

- **F-R16 (minor; Phase H4)**: the printed stratified statistics group
  "ℤ₂ × A" is false as written — the connecting homomorphism is
  nonzero (∂ = ±2·[rot], two independent routes, 28/28 machine
  checks), so π₁(𝒞_strat) = ℤ₂ × A/⟨2a_rot⟩ (A torsioned; base π₁
  absent; charged lines cannot traverse the orientation-reversing
  loop). Every consequence the ħ-closure consumes SURVIVES: the fatal
  branch is refuted, the rotation class has order exactly 2,
  χ(rot) = −1 exists splitness-free at w = ±1 — T-H4 is DISCHARGED
  largely in the corpus's favor; T-H3's benign verdict is now
  unconditional. Residual GAP of record: the N2/F.3 exchange lift
  ("ARE fermions" value-forcing). theory-audit/h4_completion.md.

- **F-R17 (moderate; Phase I4)**: App F.3's exchange-lift claim is
  refuted as stated — in the two-defect space with F.3's OWN proviso
  (separations bounded below), the ordered-line double cover's swap
  monodromy gives σ(exchange) = −1 ≠ +1 = σ(rotation): the proviso
  constructs the obstruction (the geon phenomenon; the classical FR
  null-homotopy provably exits the separated regime). Repair:
  χ_exch = χ(σ)·χ_rot with χ(σ) a new free ℤ₂ the corpus must price;
  "ARE fermions" demoted further; the B=2 worked example is
  sector-conditional (J=0 forbidden in exactly 2 of 8 sectors). F.4
  (no braiding) VINDICATED factor by factor; QUANT unaffected and
  STRENGTHENED (I4-B's parity theorem kills the fatal completion on
  ANY embedded-line moduli; I4-C extends χ(rot) = −1 to all w ≠ 0
  splitness-free). theory-audit/i4_topology.md.

Findings ledger: F-R1–F-R17 (four substantive: F-R4/F-R5 on the ħ-closure, F-R9 on the 0νββ insert, F-R14 on the gravitation sector;
one moderate: F-R6 provenance-chain; five print-level/minor).

Phase J (the filed obligations, 2026-07-16, post-Phase-I fold — no new
F-R numbers; one recovery, two confirmations, two self-corrections):
J1 ran App I.1's never-run B.4-window 𝔭-grid (first execution): the I3
edge-referencing element CONFIRMED — 𝔭 = 0.84 recovered at the benchmark
to 0.28σ under /6; the 𝔭/𝔟_eff convention item upgrades OPEN →
RECOVERED-AS-CLASS (j1_RESULTS.md). J2 executed the repaired
(4.8′)/(4.9′) arithmetic: 𝔠_sat(ε) = 3.2011(1 + 1.152ε), the ε-anchor
forfeit confirmed at 11.283σ, the Q-6′ ceiling 2.605×10⁻³ (the "×4.3"
provenance), dichotomy survives every reading — and found two stale
sentences in OUR corpus2 revision (the Q-6′ "twenty-fold"; the ε^{1/3}
shape law), both corrected in print per the shared discipline
(j2_RESULTS.md). J3 re-ran the Majoron battery at g ≈ 7.8×10⁻¹⁰:
"SHIFTS-HARMLESSLY" confirmed 17/17 with the old-window control exact;
SN/recoupling bounds never printed — spec-recovery limits stated
(j3_RESULTS.md). The residual open was then executed as J4 (2026-07-16): dense small-ε
solves measure the saturated exponent LINEAR — p = 0.997 [0.996, 1.012],
2/3 EXCLUDED (χ² ratio 8445), intercept at 64√2/9π to +7.7×10⁻⁵; the
×4.3 ceiling loosening, 60 MeV f-top and 7.8×10⁻¹⁰ Majoron edge are
MEASURED rather than linear-assumption-specific; dichotomy survives
every reading (j4_RESULTS.md). No named opens remain in Phase J.

Phase K (hardening, 2026-07-16, ROADMAP_v7 — K1/K2 adjudicated, K4's
freshness gauntlet running at fold time): K1 worked J1's remaining
/6-vs-/(15π/8) split to a VERIFIED NULL — every amplitude-bearing
archive channel < 2σ (𝔭 band 0.52σ, correcting the 0.7σ previously
quoted; 𝔟_eff ≤ 0.37σ structurally voided; bond loop 0.075σ; B.4 edges
exactly re-derived yet amplitude-insensitive) — the item closes as
"indistinguishable within the archive," with a conditional rounding
pointer printed-not-claimed (k1_RESULTS.md). K2 extended the saturated
ε-scan to ε = 0.001: the linear law HOLDS (2/3 excluded by the new
window alone, χ² 77), J4 strengthened, propagation unchanged; headline
defect K2-D1 — the 12-mode family's intercept bias drifts at small ε
(bias artifact, not exponent physics) (k2_RESULTS.md). K4's REPRODUCE
gauntlet COMPLETED: all 46 rows fresh, ZERO PHYSICS DRIFT — every root,
fingerprint, JSON physics field, PNG, and the full N=192 crossing
series bit-identical; ten process-level findings (two printed-command
defects; runtime-column drift root-caused to PYTHONHASHSEED with the
pin-fix written in; stale expectations; the shipped field3d log itself
truncated). The package REPRODUCES: DISCHARGE_PACKAGE/FRESHNESS.md is
the dated record. Phase K closed; K3 (GPU extension) remains scheduled. Out of
envelope (documented, not scheduled): blue-fog/disclination boxes and
nucleation (GPU-scale/theorem-level, per the v2 assessment's Tier-5
flag), the Σ(p) integrals from geometry, and the full fs-cosserat solver
port — the natural next engineering step once the corpus responds to
F-R4/F-R5 (the 4B pilot proved the packaging).
