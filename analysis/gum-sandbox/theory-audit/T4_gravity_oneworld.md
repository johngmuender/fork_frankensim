# T4 — THEORY AUDIT: THE GRAVITATION AND ONE-WORLD SECTORS

**Auditor:** theory-audit subagent (Fable-class), 2026-07-16.
**Charge:** first argument-level audit of Omega-paper Sec. VI (Einstein–Cartan
by theorem; Theorem VI.1 / App. E.6; the UHECR pincer; Gibbs–Duhem vacuum
energy, the sign theorem, and w(z)) and Secs. VIII/VIII′ (the triple padlock:
Lock 1 Valentini, Lock 2 elliptic no-pole, Theorem VIII.D.1, the T-KILL gate,
QFT-1/2, and the pairs-only theorem's Lemma II.1). These sectors were never
touched by the numerical campaign (Tiers 0–5a); nothing here has a frozen
referee — this memo is pure argument-and-arithmetic.

**Epistemic rules:** WITHIN-MODEL ONLY. Every verdict is about whether the
corpus's printed argument is valid given its printed premises. Nothing below
asserts anything about nature. Where a small computation settles a step, it
was run (`scratchpad/t4_checks.py`, python3/numpy; all numbers reproduced
inline).

**Sources:** `substrate-suite/01-GUM-Omega-Paper-v2.0.1.md` Secs. II.F–H, III.D,
VI, VIII, VIII′, IX, Apps. A–F, K, AUD-15; corpus files
`WS-B-Q1-Theorem-B0-Trilemma-v0_1.md`, `WS-B-Q2-Reading-Circularity-Audit-v0_1.md`,
`CA-1-Conformal-Gravity-Comparison-Audit-v0_1.md`; `SESSION_HANDOFF.md` §§1–3
(F-R1–F-R8 context, esp. F-R5).

**Verdict vocabulary:** SOUND (valid given premises) / GAP (missing step,
stated precisely) / DEFECT-CANDIDATE (the step fails; counterexample or
derivation given) / NOT-AUDITABLE-FROM-TEXT (premises unrecoverable; what is
missing stated).

---

## 0. Summary table

| # | Argument | Verdict | One line |
|---|---|---|---|
| G1 | VI.A Einstein–Cartan by theorem | GAP | kinematic dictionary is theorem-grade import; the *dynamics* step (EC field equations) is delegated to a sketch that carries G2's gap |
| G2 | Thm VI.1 slaved cone + App. E.6 | SOUND given premises / **GAP in E.6** | convexity algebra is right; the premise 𝔞₁^(s) > 0 is asserted, not derived, and is false in known field classes |
| G3 | The UHECR pincer (6.3) | GAP (arithmetic ✓) | every number checks; the max-over-*all*-sectors step is uncovered by the species-specific bounds |
| G4 | VI.D Gibbs–Duhem cancellation | NOT-AUDITABLE-FROM-TEXT (core premise) / SOUND (identity + B.2 down-payment) | the identity ω = −P is exact; *why the gravitating source is ω, not ε* is imported with no in-model derivation |
| G5 | VI.D sign theorem (α > 0) | **DEFECT-CANDIDATE as printed** (repairable) | the printed expansion drops the linear term −n₀δμ; conclusion recoverable via the Volovik q-variable expansion |
| G6 | VI.D residual shape (6.4) | SOUND given premises (verified) | w = −1 + 𝔴Ω_m(z) follows *exactly* (𝔴 = 1) for the pure tracker — stronger than the paper prints |
| G7 | VI.D "DESI-era direction check" | **DEFECT-CANDIDATE** | the family's own CPL projection gives w_a > 0 (+0.63𝔴); the cited preference is w_a < 0 — opposite sign, and the cited (w₀, w_a) pair is growth-*anti*-correlated at z ≳ 0.4, the family's own kill direction |
| W1 | VIII.A conserved exponential | NOT-AUDITABLE-FROM-TEXT (as theorem) | components are real cited imports; "the exponential is conserved" is never formalized across ontology classes |
| W2 | Lock 1 (Valentini premises as used) | SOUND given imports | premise chain complete and honestly grade-limited [IM+CAL] at the equilibrium hinge |
| W3 | Lock 2 (−1/k² no-pole; (c/c_L)⁵) | SOUND core / GAP-minor | no-pole argument exact at c_L = ∞; leakage exponent 5 presupposes vanishing monopole+dipole, unproven in-text; arithmetic ✓ 10⁻²⁰ |
| W4 | Thm VIII.D.1 (no runaways) | SOUND (scoped) | reality-below-gap ⟹ even-analyticity ⟹ no ω³ term is a correct chain; valid only for all-gapped channels, which the paper respects (E.4 keeps Larmor) |
| W5 | Thm VIII.1 + the S3/A3 gate | SOUND (VIII.1) / **GAP (gate forward implication)** | the commutation theorem is trivially valid; "S3 observed ⟹ non-POVM" needs an all-POVM exclusion, but only an enumeration over *proposals* is offered |
| W6 | QFT-1 / QFT-2 (nucleation covariance) | GAP (two debts, both printed) | proof skeleton valid; conditional on all-orders spacelike commutativity (WS-QFT-2, open) and the POVM-record scoping (F-A15-9) |
| W7 | VIII′.1 pairs-only via Lemma II.1 / App. F.1 | **DEFECT-CANDIDATE (printed proof)** / SOUND (conclusion, repaired) | the convexity claim is false — explicit counterexample, verified numerically; the conclusion survives by a history-path argument the corpus did not make |
| W8 | Threshold (8′.2) ω_th = √(7/3)ω₀ | SOUND arithmetic; F-R5-contingent value | 2√(7/12) = √(7/3) = 1.527525 ✓; under F-R5's saturated closure the same formula gives √2 ω₀ |

Credit side, stated per campaign rules: every piece of pure arithmetic in these
two sectors that this audit could check, checked (12/12 — §5); three of the
four load-bearing *identities* (Gibbs–Duhem, B.2 dust, VIII.1 commutation,
VI.1 convexity-given-weights) are exactly right; and the corpus itself printed
the two largest one-world gaps (WS-QFT-2, the equilibrium hinge) before this
audit existed.

---

## PART I — GRAVITATION (Sec. VI)

### G1. "Einstein–Cartan by theorem" (VI.A)

**Skeleton.**
- P1 [IM]: KBKK dictionary — in a defected continuum, dislocation density =
  torsion, disclination density = curvature of an effective Riemann–Cartan
  connection. (True at theorem grade in the continuum-defect literature;
  legitimate import.)
- P2: GUM's coarse-grained defect distributions therefore define (g, T^λ_μν);
  knots carry Cosserat spin density conjugate to antisymmetric stress.
- C: "The stiff-medium field equations are **Einstein–Cartan**: torsion
  algebraic, GR in torsion-free vacuum."

**Analysis.** P1+P2 deliver *kinematics*: a defected micropolar medium has a
natural RC-geometric description, and its spin density is the natural torsion
source. They do not deliver *dynamics*: which action governs (g, T) is a
separate question, and "Einstein–Cartan" is a specific answer (Einstein–Hilbert
in the RC connection, no independent torsion kinetic term). The paper's route
to the dynamics is VI.B's Sakharov induced action (6.1) — a one-loop sketch
whose R-coefficient positivity is exactly the E.6 gap (G2). Torsion
non-propagation ("algebraic") holds *iff* the induced action contains no
torsion kinetic term; nothing in the text shows the substrate loop expansion
suppresses ∂T terms — plausible (they are higher-dimension), not proven.

**Verdict: GAP.** Missing step, stated precisely: a derivation that the
substrate's induced action for the defect geometry is EC-form — specifically
(i) sign/positivity of the induced Λ_s²𝔞₁ R term (delegated to E.6, itself a
gap), (ii) absence/suppression of independent torsion kinetics. The audited
compliance riders are arithmetic-consistent: spin–torsion corrections at
∼10⁵⁴ kg/m³ vs nuclear 4×10¹⁷ kg/m³ = **36.4 orders** (paper "∼36" ✓; the
10⁵⁴ figure itself matches the standard Hehl-class EC critical density for
fermions — order-of-magnitude consistent as an import).

### G2. Theorem VI.1 (slaved graviton cone) and App. E.6

**Skeleton.**
- P1: induced TT graviton kinetic operator = sector sum with weights
  w_s = N_s Λ_s² 𝔞₁^(s) (from expanding Σ_s √(−g_s)R[g_s], g_s = ḡ + 2δ_s c²n⊗n).
- P2 (App. E.6, in full): "w_s > 0 ⟹ the TT kinetic operator is a convex
  combination." Justification offered: "heat-kernel positivity."
- C1: c_GW² = Σw_s c_s²/Σw_s ∈ [min_s c_s², max_s c_s²].
- C2: |Δc_GW/c| ≤ max_s |δ_s(E ∼ Λ_UV)| ≡ δ_UV.

**Analysis of the algebra (SOUND given P1, P2).** Given positive weights, the
weighted mean lies in the hull; with c_s² = c²(1 + 2δ_s), Δc_GW/c ≈ Σwδ/Σw
and |Δc_GW/c| ≤ max|δ_s|. Measuring against the photon cone (B2 is itself one
of the summed sectors) costs at most a factor 2 in the bound (spread vs
deviation) — harmless at the quoted order. The step is valid.

**Analysis of P2 (GAP — the load-bearing line).** "Heat-kernel positivity" is
a true statement about 𝔞₀ (the heat kernel trace and its leading Weyl
coefficient are positive). It is **not** a theorem about 𝔞₁. The Seeley–DeWitt
𝔞₁ integrand is tr(R/6 − E), whose sign is field-content-dependent:

- non-minimally coupled scalar: coefficient (1/6 − ξ) — **negative for ξ > 1/6**;
- Dirac-type operator (E = R/4·𝟙₄): tr gives 4(1/6 − 1/4) = **−1/3** before
  the fermion-loop sign, i.e. sign depends on the statistics bookkeeping;
- the sign-indefiniteness of the Sakharov-induced Newton constant is a
  standard, documented feature of induced gravity (Visser's modern-perspective
  review states it as the central caveat).

Within GUM this matters twice: (i) the sector sum includes "knot band-edges",
and knots are *fermionic by the corpus's own Theorem c″* — the one sector class
where the naive bosonic 𝔞₁ > 0 intuition demonstrably cannot be assumed;
(ii) acoustic metrics for micropolar branches generically arise with
non-minimal couplings (disformal/ξ-type terms), and no computation of the
effective ξ_s is anywhere in the text or archive. If any w_s < 0, the mean is
no longer convex, c_GW² can exit the hull, and the slaving bound fails as an
inequality (it survives only as an order-of-magnitude estimate when no large
weight cancellation occurs).

**Verdict: Theorem VI.1 SOUND given premises; App. E.6 GAP** — the premise
𝔞₁^(s) > 0 for every summed sector is asserted at one line under a
positivity label that does not cover it; a per-sector 𝔞₁ computation (a
finite, doable calculation — see §5 checkable items) is the missing step.
This gap is *upstream of the sector's headline theorem* and also feeds G1.

### G3. The UHECR pincer (6.3)

**Skeleton.** (1) Loop weights w_s ∝ Λ_s² peak at the cutoff, where the C1
flow makes δ largest. (2) ħc/a ≳ 2×10¹⁹ eV for a ≲ 10⁻²⁶ m ⟹ UHECRs probe δ
essentially at the cutoff. (3) Observed UHECR propagation (no vacuum
Cherenkov, no photon decay) pins δ_UV ≲ 10⁻¹⁹. (4) Hence |Δc_GW/c| ≲ 10⁻¹⁹,
beating GW170817's 10⁻¹⁵ with ≥ 4 orders of margin.

**Arithmetic (all checked, python3):**
- ħc/a at a = 10⁻²⁶ m = **1.97×10¹⁹ eV** — the "≳ 2×10¹⁹ eV" line is exact. ✓
- Threshold-scaling magnitudes for the bounds behind step (3):
  proton vacuum-Cherenkov class δ ≲ (m_p/E)² = 3.5×10⁻²² at E = 5×10¹⁹ eV
  (2.2×10⁻²³ at 2×10²⁰); photon-decay class δ ≲ (2m_e/E_γ)² = 1.0×10⁻¹⁶ at
  E_γ = 10¹⁴ eV. The quoted δ_UV ≲ 10⁻¹⁹ sits *inside the span* of the
  species-pair bounds — defensible as a conservative summary for
  cutoff-adjacent species pairs, and the paper's "species-pair-dependent and
  partially one-sided" rider is honest. ✓
- Margin: log₁₀(10⁻¹⁵/10⁻¹⁹) = **4.0 exactly** — "≥ 4 orders" is
  boundary-exact against the round 10⁻¹⁵; against the tight side of the actual
  GW170817 band (−3×10⁻¹⁵, +7×10⁻¹⁶) it is 3.85 orders. Cosmetic. ✓

**The logical gap.** Step (3) bounds the cone splittings of the sectors UHECR
observables couple to (proton-constituent strata, photon). Step (4) needs
**max over all sectors s in the sum** — including B3, B4, and the knot
band-edge sectors, none of which is directly probed by UHECR kinematics. The
bridge would be C1-flow universality (all splittings run together toward the
common cone), but the C1 flow as printed (II.G) gives a *rate*, not a UV
anchor: two sectors can share the flow yet have different δ_UV. The pincer as
stated silently substitutes "the δ's UHECRs see" for "max_s δ_s". Also, the
partial one-sidedness the paper flags matters more than it lets on: a
one-sided bound on δ for one species pair does not bound |δ| in the convex
mean, and P-M3's "definite sign" prediction leans on the very loop-weight
signs G2 leaves unproven.

**Verdict: GAP** (coverage step from probed species to all summed sectors),
with all printed arithmetic verified correct.

### G4. Gibbs–Duhem vacuum-energy cancellation (VI.D)

**Skeleton.**
- P1 [IM, Volovik-class]: for a self-sustained medium, the vacuum source of
  the induced gravitational equations is the grand-potential density
  ω_vac = ε_vac − Σμ_i n_i, not ε_vac.
- P2 (thermodynamic identity): Gibbs–Duhem ⟹ ω = −P_ext.
- P3: self-sustained ⟹ equilibrium at P_ext = 0.
- C: ω_vac = 0 without fine-tuning; only departures gravitate.
- Down-payment (App. B.2): on the Bogomolny locus P = ½Λ²b² − 𝒱 = 0 pointwise.

**Analysis.** P2 is an exact identity (ε + P = Ts + μn at T = 0 gives
ω ≡ ε − μn = −P); P3 is the definition of self-sustainment. Given P1 the
conclusion follows — the argument is *valid*. The entire weight is on P1, and
P1 has **no in-model derivation anywhere in the paper or archive**: nothing in
GUM's induced-action construction (6.1) is shown to produce ω rather than ε as
the source. In Volovik's own q-theory this step is a specific dynamical
argument about a conserved vacuum variable entering the action; the paper
imports the conclusion, labels it [IM] honestly, but then builds a "sign
theorem" and a w(z) family on it as if the identification were secured.

The down-payment is genuinely sound: at Bogomolny saturation Λb = √(2𝒱), so
½Λ²b² − 𝒱 = 0 identically — the compensation pattern does operate natively in
the (6+0) sector. (Checked symbolically; trivial.) Note its limitation: it
shows *P = 0 for the knot sector on-shell*, which is the EOS statement (dust),
not the ω-sourcing statement P1.

Also checked, the "nightmare" sizing: ħc/a⁴ at the section's own a ≲ 10⁻²⁶ m
gives **3.2×10⁷⁸ J/m³**, not the printed "≳ 10⁹⁰ J/m³" — the 10⁹⁰ figure
corresponds to a = 1.33×10⁻²⁹ m, which appears nowhere. Rhetoric-only slip
(the catastrophe is 88 orders instead of 99; the argument is unaffected), but
it is a checkable number printed wrong from its own stated input.
**Finding F-T4-a (minor, arithmetic).**

**Verdict: core premise P1 NOT-AUDITABLE-FROM-TEXT** (what is missing: the
derivation that GUM's induced-gravity source is the grand potential — the
in-model analog of Volovik's q-theory argument); **the identity chain P2–P3–C
SOUND given P1; App. B.2 down-payment SOUND.**

### G5. The sign theorem (α > 0)

**As printed:** "Expanding the grand potential about the self-sustained point,
ω(μ) ≈ ½χ⁻¹(μ−μ₀)² with χ > 0 by substrate stability; Hubble dilution drives
a lag δμ ∝ H/Γ; hence ρ_Λ = ½χ⁻¹δμ² > 0."

**The defect.** ∂ω/∂μ = −n. An expansion of ω in μ about μ₀ has leading term
**−n(μ₀)δμ**, which is *linear and of indefinite sign* unless n(μ₀) = 0. But
n(μ₀) ≠ 0 is the very heart of the mechanism (the cancellation is between
large ε and large μ₀n₀ — if the vacuum's conserved density vanished there
would be nothing to cancel). As printed, the expansion is therefore wrong: the
quadratic form is claimed in the one variable whose linear response does not
vanish, and the sign conclusion does not follow. **DEFECT-CANDIDATE, by
derivation.**

**The repair (available, not printed).** Volovik's own expansion is in the
conserved vacuum variable q at fixed μ = μ₀: ρ_Λ(q) = ε(q) − μ₀q satisfies
dρ_Λ/dq|_{q₀} = ε′(q₀) − μ₀ = 0 *by the equilibrium condition itself*, and
d²ρ_Λ/dq² = ε″(q₀) = χ⁻¹ > 0 by stability, so a lag δq gives
ρ_Λ ≈ ½χ⁻¹δq² > 0. The theorem's conclusion is recoverable by transcribing
μ → q; the printed proof is not. (Whether the Hubble lag lives in q with
δq ∝ H/Γ is then a further [DW] modeling step, as the paper elsewhere admits
via the WS4-M4c debt.)

**Verdict: DEFECT-CANDIDATE as printed; conclusion recoverable under the
q-variable repair, at which point the grade is GAP** (the lag mechanism and
the χ, Γ ↦ M\*² dimensional bridge remain asserted). Magnitude check run:
with M\*² = c⁴/8πG and H₀ = 70 km/s/Mpc, matching ρ_Λ = 6×10⁻¹⁰ J/m³ needs
α = **0.73** — the O(1) claim is arithmetically fair. ✓

### G6. The residual shape (6.4): w(z) ≈ −1 + 𝔴·Ω_m(z)

**Checked by direct computation** (ρ_DE = M\*²(αH² + βḢ), separately
conserved, ΛCDM-like background, w = −1 − (1/3)dlnρ/dlna):

- Pure tracker (β = 0): the computation gives **w(z) = −1 + Ω_m(z) exactly**
  (numerical max deviation 1.4×10⁻⁴ = grid error; analytically,
  dlnH²/dlna = −3Ω_m on a matter+Λ background). So (6.4) holds with
  **𝔴 = 1 exactly** for the pure αH² term — a *stronger and cleaner* result
  than the paper's "𝔴 = O(1)" with its value left as a debt. The paper
  underclaims its own shape here.
- With the memory term: β = −0.3 gives 𝔴(z) = 1.28 → 1.03 over z = 0 → 2;
  β = −1.0 gives 1.72 → 1.05: 𝔴 = O(1), mildly running, w > −1 throughout,
  approaching −1 as Ω_m → 0. The printed shape and O(1) claim are **SOUND
  given the premises** (the ρ_DE form and separate conservation).
- The tracking-pathology honesty ("a pure αH² term cannot accelerate") is
  correct in its intended self-consistent sense (ρ_DE ∝ ρ_tot ⟹ constant
  fraction, no domination) — properly printed as a pathology, not hidden.

**Verdict: SOUND given premises** (β < 0 itself remains underived — the
paper's own WS4-M4c debt — so the *full* (6.4) package is GAP at exactly the
point the paper already flags).

### G7. The "direction check" sentence — DEFECT-CANDIDATE

The paper: "the DESI-era preference (w₀ > −1, w_a < 0) sits on exactly this
side of Λ."

**Computation.** The family (6.4) has w increasing with redshift (w → −1 from
above as Ω_m → 0, i.e. toward the future). Its own CPL projection at a = 1:
w_a ≡ −dw/da|₁ = −𝔴·dΩ_m/da|₁ = **+3𝔴Ω_m0(1−Ω_m0) = +0.63𝔴 > 0**
(numerically +0.36 at β = −0.3 with the memory term; +0.63 for the pure
tracker). The DESI-era CPL preference cited is w_a < 0 — **opposite sign**.
Worse for the sentence: the cited (w₀ ≈ −0.7, w_a ≈ −1) pair has
w(z) crossing below −1 at z ≈ 0.4 and (w+1) *anti*-correlated with Ω_m(z) at
higher z — which is precisely the family's own printed "kills outright"
direction. The only part of the DESI preference the family matches is
w₀ > −1 (the low-z side of Λ).

**Verdict: DEFECT-CANDIDATE** on the direction-check sentence as written: the
family does not sit on the side of the *cited pair*; it sits on the side of
w₀ alone. Honest uncertainty printed: CPL is an extrapolating parametrization,
and low-z data dominate the w₀ pull, so a defender can argue the physical
content of the DESI preference is compatible with thawing-from-above families;
but then the sentence should cite w₀ only. As printed, half the cited evidence
points at the family's own kill clause, unflagged. (Within-model: this also
tightens SWP Branch A — the pre-bound escape hatch requires "the
growth-correlated sign of (6.4)", which the current CPL-preferred sign is
*not*; the corpus's own integrity clause is working against it here, which the
corpus should be told.)

---

## PART II — ONE WORLD (Secs. VIII, VIII′)

### W1. VIII.A — the conserved exponential

The three cost results (Spreeuw's mode count; Deutsch–Hayden/Bédard descriptor
cost; configuration-space dimension) are genuine, correctly attributed
imports. The connective claim — "the exponential is conserved; only its
address changes" — is a meta-theorem over *all* realist readings and is never
formalized: no common cost measure is defined under which the three exponents
are the same conserved quantity, and no exhaustiveness argument is given in
the paper itself (the archive's WS-B Q1 trilemma supplies an exhaustiveness
*map* at [DF-structural] grade, but over ontology classes, not over cost
measures). The paper's *operational* claim — GUM pays via Norsen towers, depth
3–5 sufficing at 10⁻³ for area-law entanglement — is a ⟨r1⟩ [CAL] numeric
(campaign credit side; not re-run here).

**Verdict: NOT-AUDITABLE-FROM-TEXT as a theorem** (missing: a defined cost
functional and a proof it is invariant across the three addresses); the
individual imports and the tower numerics stand at their printed grades.

### W2. Lock 1 — Valentini's theorem, premises as used

**Skeleton.** Valentini's signal-locality theorem: marginals at B independent
of spacelike operations at A **iff** ρ = |ψ|². Premises as used by the paper:
(i) leaf-wise equivariance (from HBD dynamics on the supplied foliation —
standard import, valid); (ii) equilibrium ρ = |ψ|², carried by III.D =
imported H-theorem [IM: Valentini] + simulation-anchored rates [CAL];
(iii) operations at A are local (Hamiltonian/CP on A's factor).

**Analysis.** The theorem itself is standard and its leaf-wise transcription
is legitimate given (i). The load-bearing premise is (ii), and the paper is
honest about its grade: the kinematic H-theorem is imported, the rates are
[CAL], and "equilibrium is coarse-grained — always violated below scale a" is
printed. The archive's own circularity audit (WS-B Q2, S8) stamps this chain
"DERIVED-AT-EQUILIBRIUM" with the hinge inherited from S7's [IM+CAL] — this
audit concurs with that stamping. The remaining unclosed premise, also flagged
in-corpus (WS-B Q2, D1): *why the post-bounce initial condition is
near-equilibrium* — without it, Lock 1 is a fixed-point stability statement,
not an unconditional lock. The self-application flourish ("agents are
equilibrium knots") adds nothing load-bearing.

**Verdict: SOUND given imports**, with the equilibrium hinge correctly
grade-limited by the corpus itself; residual GAP = the initial-condition
premise (already on the corpus's own debt list).

### W3. Lock 2 — the −1/k² no-pole argument and the (c/c_L)⁵ leakage

**Skeleton.** (i) At c_L = ∞ the constraint sector is elliptic:
∇²p = −s, propagator −1/k², no ω-dependence, hence no pole surface
ω = ω(k): no radiation for *any* source motion — "a solver, not a channel."
(ii) At finite c_L (Bell-baseline floor c_L/c ≳ 10⁴ [IM: Salart-class]),
leakage ≤ (c/c_L)⁵ ≤ 10⁻²⁰ (App. E.5: "quadrupole leakage").

**Analysis.** (i) is exactly right as mathematics: an elliptic constraint has
no propagating on-shell surface, its Green function does no far-field work in
the wave-zone sense. The energy face of the corollary follows. (ii) arithmetic:
(10⁻⁴)⁵ = **10⁻²⁰ exactly** ✓. The exponent 5 is the *quadrupole* power ratio
(P_quad ∝ ⃛Q²/c⁵ channel-by-channel, so the c_L-channel emission is (c/c_L)⁵
of the reference). This presupposes the monopole and dipole channels vanish:
monopole radiation would leak at (c/c_L)¹ ∼ 10⁻⁴ and dipole at (c/c_L)³ ∼
10⁻¹², both fatally larger than the printed bound. Vanishing requires the
pressure-channel source to carry conserved total strength (monopole) and
conserved first moment (dipole) — physically plausible if s is built from
conserved knot mass/momentum densities, but **nowhere exhibited**; App. E.5
just says "quadrupole leakage" as if the multipole order were established.

**Verdict: SOUND core (the c_L = ∞ statement and the corollary); GAP-minor**,
stated precisely: a proof that the constraint source's monopole and dipole
moments are non-radiating (conservation laws of the B1 channel), which fixes
the leading leakage exponent at 5. Absent it, the honest bound from the text's
own premises is (c/c_L)³ ≤ 10⁻¹² (momentum conservation being the easier
half), still small but eight orders off the printed number.

### W4. Theorem VIII.D.1 — no runaways by even-analyticity

**Skeleton.** The self-force response kernel χ(ω) of a *gapped* channel:
(i) below the cut (|ω| < ω₀) there is no absorption, so χ is real there;
(ii) reality of the time-domain kernel gives χ(−ω) = χ*(ω); (i)+(ii) ⟹ χ even
and real-analytic below the cut ⟹ Taylor series in ω² only ⟹ the derivative
expansion of the self-force contains only even time-derivatives (conservative
mass-renormalization-class terms), and **the ω³ Abraham–Lorentz term cannot
appear** ⟹ no radiation-reaction runaways/preacceleration in that channel.

**Analysis.** The chain is correct and each link is standard: below-threshold
non-absorption ⟹ Im χ = 0 there; reality + Im = 0 ⟹ evenness; gap ⟹
analyticity at ω = 0. Two scope conditions, both of which the paper respects:
(a) it holds per-channel and only for channels with a gap — the gapless B2
(EM) channel keeps its Larmor term, which App. E.4 states explicitly, so the
theorem does not overclaim against its own electrodynamics; (b) "no runaways"
is here the statement that no odd low-order term generates the classic
third-derivative instability — the full nonperturbative statement (the
resummed kernel defines stable dynamics) is stronger, but within the derivative
expansion the printed claim is exactly what evenness buys. Also consistent
in-context: the same evenness is what rejected the superradiant closure rung
(IV.D), i.e. the corpus uses the theorem where it bites its own numbers —
a good sign for within-model coherence.

**Verdict: SOUND (scoped as printed).**

### W5. Theorem VIII.1 and the S3/A3 gate

**Theorem VIII.1 (POVM order-independence).** If Bob's detector is a POVM
E_B(dt) on H_B, then for spacelike Alice operations, local CP maps on distinct
tensor factors commute, so Tr[ρ(E_A^± ⊗ E_B(dt))] is order-independent.
**SOUND** — this is elementary and exact ((E_A ⊗ id) and (id ⊗ E_B) commute
as superoperators). The paper's use of it is also honest: it *proves the
hiding* (why the foliation is invisible to POVM statistics) rather than
asserting it.

**Prop. VIII.2** (first-crossing statistics are not POVMs on H_B; 9.6%
near-field order-dependence) — the structural point is right (a functional of
the conditional wave function depends on Alice's F-time and is not of the form
Tr[ρ_B E_B]); the 9.6% is ⟨r1⟩ [CAL], not re-run here (campaign's tier-3/T-KILL
replication touched the 2×10⁻⁵ ABR side).

**The gate's forward implication — GAP.** The gate is stated as a
biconditional: "*S3 observed* ⟺ the physical detector implements a non-POVM
arrival observable." The reverse direction (¬S3 + ABR ⟹ T-KILL) is a clean
modus tollens on GUM's own prediction — fine. The forward direction requires:
**no POVM on H_B can produce the observed statistic**. What the text offers
(IX.A) is "a hard support cutoff τ_max — absent from every POVM *proposal*."
Enumeration over proposals is not exclusion over the POVM class: POVMs as a
class can realize essentially arbitrary outcome distributions, including
spin-dependent compactly-supported ones (nothing in positivity or normalization
forbids hard support; e.g. a time-gated detection model is a POVM with a
cutoff). The known analytic results (Kijowski-class arrival POVMs are
support-unbounded) cover natural families, not all of them. So "S3 observed"
would establish "orthodox proposals fail," not "Theorem VIII.1's premise fails
in nature" — and Layer-2's license (the whole point of the gate) needs the
latter. Missing step, stated precisely: a theorem that the *specific
spin-covariant cutoff statistic* (jointly with the far-field null and the
X-Vb1 consistency value) is unrealizable by any POVM on H_B compatible with
the experiment's preparation class.

**Verdict: VIII.1 SOUND; Prop. VIII.2 structure SOUND / numeric [CAL];
the gate's forward implication GAP** (quantifier slip from "every proposal" to
"every POVM").

### W6. QFT-1 / QFT-2 — nucleation covariance

**Skeleton (QFT-1).** For admissible foliations F, F′: (i) leaf-wise
equivariance of the jump process (minimal-jump-rate construction — a standard
Bell-type-QFT result, legitimate import) makes record statistics functionals
of Ψ alone; (ii) Ψ evolves by multi-time Tomonaga–Schwinger dynamics whose
integrability = spacelike commutativity of H_int densities — **verified at
effective order only**, all-orders open (WS-QFT-2, printed); (iii) records are
configuration facts, POVM-representable on a completion leaf (scoping
tightened by the corpus's own F-A15-9). Conclusion: P_F[O] = P_{F′}[O].

**Analysis.** The proof skeleton is the correct one (it mirrors the
hypersurface-Bohm–Dirac / Bell-type-QFT equivariance results in the cited
literature), and the paper's honesty block is genuinely load-bearing: it
states at equal size that the jump *micro-history* is foliation-relative, and
it scopes records to POVM-representables — exactly the two places where a
stronger claim would be false (Prop. VIII.2's non-POVM functionals are
order-dependent, so unscoped "records" would contradict the corpus's own
proposition; AUD-15 caught this pre-release, F-A15-9). What remains is exactly
what the corpus prints as debt: (a) all-orders spacelike commutativity with
IBC boundary terms (WS-QFT-2) — without it QFT-1 is proven only at effective
order; (b) the admissibility class of foliations is defined by fiat
("spacelike, asymptotically the substrate rest slicing").

**QFT-2 (triple padlock)** = QFT-1 + Lock 1 + Lock 2 conjoined; it inherits
their statuses (W2's equilibrium hinge, W3's multipole gap, W6a's all-orders
debt). No new step to audit.

**Verdict: GAP** — the two conditions are real and both are printed by the
corpus itself; given them, the argument is valid. (Within-model, this is the
best-behaved theorem in the audited set: every soft spot is self-flagged.)

### W7. VIII′.1 pairs-only theorem and Lemma II.1 (App. F.1) — the printed proof fails

**The printed proof (App. F.1):** "Admissible deformations (det F > 0,
identity at ∞) form a **convex, hence contractible, set**; R̃[u] is continuous
into Maps(ℝ³, SU(2)) from a contractible domain ⟹ deg R̃ ≡ deg R̃[0] = 0;
deg P̃ = deg Q̃. Contrapositive: ΔΣK ≠ 0 ⟹ det F → 0" — the pairs-only
theorem's engine.

**DEFECT-CANDIDATE — the convexity claim is false.** Counterexample
(constructed, then verified numerically): let φ be a twist map — rotation by
angle θ(r) about the x₃-axis with θ = π inside a ball, θ → 0 outside — a
smooth deformation, identity at ∞, with det F = 1 > 0 everywhere (each sphere
rotates rigidly). φ is admissible; so is the identity. Their convex
combination at t = ½ has, in the core region,
F_½ = ½R_z(π) + ½I = diag(0, 0, 1):

  det((1−t)I + tR_z(π)) at t = 0, ¼, ½, ¾, 1: **+1.0000, +0.2500, 0.0000,
  +0.2500, +1.0000** (numpy, exact zero at t = ½).

The straight-line segment between two admissible deformations exits the
admissible set: the set is **not convex**, and "convex, hence contractible"
collapses. (This is just the classical fact that GL⁺(3) is not convex,
imported pointwise.)

**Is the conclusion still true? Yes — by an argument the corpus did not make.**
Two repairs, in increasing strength:

1. **History-path repair (sufficient for the theorem, elementary).** The
   pairs-only theorem quantifies over *physical histories*: the medium's
   configuration evolves continuously from its initial state, so u(t) is a
   continuous path with det F(t, x) > 0 along it. Degree is a locally constant
   functional along any continuous path *within* the admissible set — no
   convexity, contractibility, or even connectedness of the whole set is
   needed. The contrapositive (ΔΣK ≠ 0 along a history ⟹ admissibility fails
   somewhere ⟹ det F → 0) follows immediately. **The theorem's conclusion is
   SOUND under this repair**, and the repair is strictly weaker than what the
   corpus tried to prove.
2. **Set-level repair (what F.1 actually claimed).** If "admissible" means
   global orientation-preserving embeddings decaying to identity
   (non-interpenetration of matter — the physically motivated class), the
   space *is* contractible, but by a deep theorem (the Smale-conjecture
   corollary: Diff_c(ℝ³) is contractible, Hatcher 1983), not by one-line
   convexity. If "admissible" means only pointwise det F > 0 (local
   diffeomorphisms, interpenetration allowed), the set-level claim is in real
   danger: by Gromov's h-principle for open Diff-invariant relations the
   solution space is weakly equivalent to Maps_c(ℝ³, GL⁺(3)) ≃ based
   Maps(S³, SO(3)), whose components are π₃(SO(3)) = ℤ — a *disconnected*
   admissible set, on which deg R̃ plausibly labels components and the lemma's
   set-level statement would be false. The corpus's text does not specify
   which class it means.

**Verdict: printed proof DEFECT-CANDIDATE (counterexample above); theorem
VIII′.1's conclusion SOUND after the history-path repair.** Downstream users
(the IBC/coincidence-stratum construction, "particles are created in pairs")
survive, because they only ever use the history form. Recommended erratum: replace
F.1's convexity sentence with the path argument (two lines, and it is the
argument the physics actually licenses).

### W8. Threshold (8′.2) — arithmetic and an F-R5 cross-link

- 2κ_phys = 2√(7/12) = √(7/3) = **1.527525** ✓ (printed 1.5275(15); AUD-15's
  one-sided band correction applies).
- "born already 23.6% bound": 1 − √(7/12) = **0.2362** ✓.
- **Cross-link (new, within-model):** κ_phys = √(7/12) is the closure value
  the campaign's F-R5 showed is *not* selected by the corpus's own variational
  problem (the unrestricted minimizer runs to halo saturation at
  κ_paper = 1/√2). If the saturated closure is used instead, the same formula
  gives ω_th = 2·(1/√2)ω₀ = **√2 ω₀ = 1.4142**, binding 29.3%. The threshold
  is therefore an F-R5-contingent number: the *argument* (bound-pair mass sets
  the spectator-assisted threshold; Proposition Ω-1 forbids vacuum pair
  creation; the √-onset from Thm VIII.D.2) is unaffected, but the numeric
  value inherits the campaign's substantive finding. The nucleation sector
  nowhere flags this dependence.

**Verdict: arithmetic SOUND; value carries an unflagged F-R5 dependency**
(inherits DISCHARGE_PACKAGE routing, not a new defect class).

---

## 5. Checkable-arithmetic ledger (all run, `scratchpad/t4_checks.py`)

| Item | Computed | Paper | Status |
|---|---|---|---|
| (c/c_L)⁵ at c_L/c = 10⁴ | 1.000×10⁻²⁰ | ≤ 10⁻²⁰ | ✓ |
| ħc/a at a = 10⁻²⁶ m | 1.973×10¹⁹ eV | ≳ 2×10¹⁹ eV | ✓ |
| GW margin log₁₀(10⁻¹⁵/10⁻¹⁹) | 4.0 (3.85 vs the +7×10⁻¹⁶ side) | ≥ 4 orders | ✓ boundary-exact |
| UHECR species bounds (m/E)² | p: 2.2×10⁻²³–3.5×10⁻²²; γ: 1.0×10⁻¹⁶ | δ_UV ≲ 10⁻¹⁹ | ✓ in-span (coverage gap separate) |
| spin–torsion 10⁵⁴ vs 4×10¹⁷ kg/m³ | 36.4 orders | ∼36 | ✓ |
| zero-point ħc/a⁴ at a = 10⁻²⁶ m | 3.2×10⁷⁸ J/m³ (10⁹⁰ ⟹ a = 1.33×10⁻²⁹) | ≳ 10⁹⁰ J/m³ | ✗ **F-T4-a** (rhetoric-only) |
| α for ρ_Λ = 6×10⁻¹⁰ J/m³, M\*² = c⁴/8πG | 0.73 | O(1) | ✓ |
| pure-tracker w(z) | −1 + Ω_m(z) exactly (dev 1.4×10⁻⁴ = grid) | ≈ −1 + 𝔴Ω_m, 𝔴 = O(1) | ✓ (stronger: 𝔴 = 1) |
| 𝔴 with β = −0.3 / −1.0 | 1.28→1.03 / 1.72→1.05 over z = 0→2 | O(1) | ✓ |
| family's CPL w_a | **+0.63𝔴** (+0.36 at β = −0.3) | cited DESI w_a < 0 as favorable | ✗ **G7 defect** |
| det((1−t)I + tR_z(π)) | 0.0000 at t = ½ | F.1 "convex" | ✗ **W7 defect** |
| 2√(7/12) = √(7/3); 1 − κ | 1.527525; 0.2362 | 1.5275(15); 23.6% | ✓ |
| sidereal τ\* = vd/c² (30 m, 370 km/s) | 0.1235 ns | 0.123 ns | ✓ (bonus) |
| ΔN_eff (4/7)(10.75/106.75)^{4/3} | 0.0268 | 0.0268 | ✓ (bonus) |

## 6. Uncertainty statement

- G2/G1: the 𝔞₁-sign gap is a *gap*, not a demonstrated failure — a per-sector
  computation could come out all-positive; this audit asserts only that the
  positivity is currently unproven and that the offered justification
  ("heat-kernel positivity") does not cover it. Confidence the criticism is
  correct as stated: high (standard Seeley–DeWitt facts).
- G5: the defect is in the printed derivation; confidence the *intended*
  Volovik argument goes through under the q-transcription: high, but that
  argument's own premise (P1, G4) remains an import.
- G7: confidence in the sign computation: high (analytic + numeric agree).
  Whether it *matters* depends on how much weight "direction check (no
  detection claimed)" was meant to carry; the SWP Branch A interaction makes
  it more than cosmetic.
- W5: the gap claim rests on "POVMs can have hard-support distributions,"
  which is elementary; the possibility remains that the *full* S3 statistic
  (cutoff + covariance structure jointly) is provably non-POVM — that proof
  would discharge the gap and should be requested from the corpus.
- W7: the counterexample is machine-verified; the h-principle remark in
  repair 2 is offered at sketch grade (flagged as such) and is not load-bearing
  for the verdict — the history-path repair is elementary and suffices.
- Prop. VIII.2's 9.6%, the 8×10⁻⁴ D3 RMS, and tower-depth numerics are ⟨r1⟩
  [CAL] values outside this memo's charge; they are neither confirmed nor
  challenged here.

## 7. Findings routed to the campaign ledger (proposed, adjudicator's call)

- **F-T4-1 (moderate, proof-integrity):** App. F.1's convexity proof of
  Lemma II.1 is false (verified counterexample); conclusion recoverable by the
  history-path argument. Erratum-class fix; pairs-only theorem survives.
- **F-T4-2 (moderate, argument-integrity):** VI.D sign theorem's printed
  expansion drops the −n₀δμ linear term; valid only after transcription to the
  q-variable. Erratum-class; conclusion survives under the repair.
- **F-T4-3 (moderate, evidence-integrity):** the (6.4) "DESI direction check"
  cites (w₀ > −1, w_a < 0) as favorable while the family's own CPL projection
  gives w_a > 0; the cited pair is growth-anti-correlated at z ≳ 0.4 — the
  family's kill direction. Interacts with SWP Branch A's pre-bound clause.
- **F-T4-4 (minor, gap-registration):** App. E.6's 𝔞₁^(s) > 0 is asserted
  under a positivity label that does not cover it; per-sector 𝔞₁ computation
  should join the open register beside R-M3b (it is "arithmetic" in exactly
  the corpus's own sense).
- **F-T4-5 (minor, arithmetic):** VI.D's "≳ 10⁹⁰ J/m³" does not follow from
  a ≲ 10⁻²⁶ m (gives 3.2×10⁷⁸); rhetoric-only.
- **F-T4-6 (minor, scope):** Lock 2's (c/c_L)⁵ presupposes non-radiating
  monopole/dipole of the constraint source (unproven in-text; conservative
  fallback exponent 3 ⟹ 10⁻¹²).
- **F-T4-7 (minor, quantifier):** the A3/S3 gate's forward implication needs
  all-POVM exclusion; only proposal-enumeration is in the text.
- **F-T4-8 (cross-link):** ω_th = √(7/3)ω₀ is F-R5-contingent (saturated
  closure ⟹ √2 ω₀); nucleation sector should carry the flag.

*Two-layer authorship note: this memo is RESULTS-layer argument analysis by a
subagent; severities above are proposals for the coordinator's adjudication
pass, not adjudications.*
