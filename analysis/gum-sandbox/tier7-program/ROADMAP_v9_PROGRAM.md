# ROADMAP v9 — PROGRAM EXTENSIONS (Tier 7)

Successor to ROADMAP_v8_FOUNDATIONS. Source analysis:
`PROGRAM_ANALYSIS_T7.md` (ranked open register §2). Method unchanged:
within-model, independent code, pre-registered gates, RESULTS.md per
workstream, honesty over success. Units ħ = m = 1 unless stated.

## M1 — WS4-M4c toy: 𝔴 from the Γ(H) relaxation profile [open #4; F-R12 anchor]
**Goal.** First independent instantiation of VI.D's relaxation
dynamics: the q-variable displaced from its minimum by Hubble drag,
ρ_Λ,eff = ½χ⁻¹δq², w(z) from the dynamics — testing the exact tracker
identity, the memory-term 𝔴 values, and the CPL projection the corpus
prints.
**Model.** Flat FRW, Ω_m0 = 0.3, matter + the relaxation field. Slow
mode: δq̇ = −Γ(H)δq + S·H (drive ∝ expansion rate; S sets the
amplitude, cancels in w), Γ = γH (tracker class). Effective DE:
ρ_Λ(z) = ½χ⁻¹δq²; w(z) from the continuity equation
w = −1 + (1/3)dln ρ_Λ/dln(1+z). Memory variant: Γ(H) with lag β (the
corpus's β = −0.3, −1.0 forms — implement as the printed linear-response
correction to the pure tracker and report the mapping used).
**Gates.** G1: numerics — background solved to 10⁻⁸; ρ_Λ positive
throughout (the q-variable sign theorem's content). G2 (exact tracker
identity): in the pure-tracker limit, measured w(z) = −1 + Ω_m(z) to
≤ 10⁻⁶ across z ∈ [0, 10] — the F-R12 strengthening verified on
independent dynamics, incl. the analytic derivation printed in
RESULTS.md. G3 (CPL projection): fitted (w₀, w_a) around a = 0.7 gives
w_a = +3𝔴Ω_m0(1−Ω_m0) within 5% for the pure tracker (𝔴 = 1 ⇒
w_a ≈ +0.63); sign STRICTLY positive across every profile in the
family. G4 (memory values): report measured 𝔴(z=0.5→0) for the two
corpus β-forms against the printed pairs (1.28→1.03; 1.72→1.05);
agreement within 15% = PASS, else FAIL with the mapping ambiguity
diagnosed (the corpus's β-convention is not fully specified — say so
honestly). G5 (Branch-A operational check): confirm the pre-bound
clause's operational content on this dynamics — (w+1) ∝ Ω_m(z) exactly
in the tracker limit; report the deviation shape when memory is on.
**Deliverables.** `M1/m1_wz.py`, `M1/m1_results.json`, `M1/RESULTS.md`,
w(z) figure.

## M2 — The bootstrap servo coefficient [⟨r2⟩; last quantum-core [CAL] anchor]
**Goal.** Independently measure the toy-law λ = (0.021 ± 0.004) g_χ²ω\*
(IV.H.3 ⟨r2⟩) that h27 did NOT cover (h27 verified exponents/torque/
class boundary only).
**Model.** A phase-locked self-oscillator: knot clock φ with
ω_mech(𝔥) and phase-lock torque −λ₀ sin Δφ against the medium's
ω_phase(𝔥); 𝔥 displaced by δ; coupling g_χ to a thermalizing bath of
N ≥ 100 modes (Langevin, Ohmic-class, temperature small); measure the
relaxation rate λ of δ(t) → 0 (fit exponential envelope).
Use the h27-confirmed structure: d ln ω_mech/d ln 𝔥 = +½,
d ln ω_phase/d ln 𝔥 = −½ (build these in exactly).
**Gates.** G1: integrator quality (energy drift of the deterministic
limit ≤ 10⁻⁶; bath statistics verified Gaussian). G2 (scaling):
λ ∝ g_χ² measured across ≥ 1 decade of g_χ (log-slope 2.0 ± 0.1).
G3 (linearity in ω\*): λ/ω\* invariant under ω\* rescaling ×4 (± 10%).
G4 (coefficient): measured c ≡ λ/(g_χ²ω\*) reported with error;
compared to 0.021 ± 0.004. Agreement within 2σ_combined = PASS-MATCH;
disagreement = reported honestly as PASS-MEASURED with the caveat that
the corpus's bath spec (r2) is unreleased and the coefficient is
bath-spectrum-dependent — in that case report c across ≥ 2 bath
spectral shapes to bound the convention spread. G5 (class boundary
control): integrable limit (bath off) ⇒ marginal (no relaxation) —
reproducing the corpus's honest boundary and h27's γ = 0 result.
**Deliverables.** `M2/m2_servo.py`, `M2/m2_results.json`,
`M2/RESULTS.md`, λ-vs-g² figure.

## M3 — F-R9 route (b) priced: on-tube saturation nonlinearity [0νββ insert]
**Goal.** Quantitatively price the one undischarged rescue route for
F-R9: does an on-tube amplitude saturation exist that (i) tames
R = ⟨A²⟩/⟨A⟩² from 10¹⁴–10¹⁸ to O(1) while (ii) preserving the MEAN
rate (verified at 1.000 ± 0.001) and (iii) staying inside the corpus's
own scale window?
**Model.** Reproduce h21's Campbell-process setup in reduced form:
line-supported web, rate ∝ A², A = superposition of tube contributions
(Campbell/shot-noise statistics; match h21's parametrization — read
theory-audit/h21_RESULTS.md + h21_mc.py first and REUSE its
parametrization; do not invent a new one). Add saturation
A → A_sat·tanh(A/A_sat); scan A_sat/⟨A⟩ over ≥ 4 decades; compute
R(A_sat) and the mean-rate distortion D(A_sat) = ⟨rate⟩_sat/⟨rate⟩.
**Gates.** G1: unsaturated control reproduces h21's R within its own
window (10¹⁴–10¹⁸ across the corpus scale window — reproduce at least
one anchor point to ×3) and mean 1.000 ± 0.001-class. G2 (the pricing
curve): R(A_sat) and D(A_sat) measured; report A_sat\* where R ≤ 10
and the mean-distortion there. G3 (the verdict, pre-registered
decision rule): route (b) is VIABLE iff at A_sat\* the mean distortion
|D − 1| ≤ 0.5 AND A_sat\* lies within the corpus's printed on-tube
amplitude range (extract the range from h21's parametrization; if the
range is not extractable, grade CONDITIONAL and print the required
range as the isolated assumption — the honest analog of L3's
source-form residue). G4: the lemma-level statement printed: what
saturation does to a quadratic Campbell functional (⟨A²⟩ → bounded;
R → O(A_sat²/⟨A⟩²) class), with the h21 motional-narrowing exclusion
left untouched.
**Deliverables.** `M3/m3_saturation.py`, `M3/m3_results.json`,
`M3/RESULTS.md`, R(A_sat) figure.

## M4 — Nelson non-equilibrium relaxation: τ(ν) [III.D rates; unlock windows]
**Goal.** First test of whether GUM's actual (Nelson) kinematics
relaxes to Born equilibrium faster than the dBB dynamics used in every
rates anchor — a definite-direction conservatism statement for the
corpus's equilibrium claims and relic windows.
**Model.** tier3-born's setup, reduced for cost: 2-D box, superposition
of M = 9 and M = 16 modes (the clean near-exponential cases; skip
M = 4), exact ψ(t) by mode sum; N = 10,000 particles started from a
non-equilibrium ρ₀ (tier3's choice: |ψ|² of a SUBSET of modes — read
tier3-born/REPORT.md and reuse its ρ₀ convention and H̄ estimator:
coarse-grained H̄(t) = ∫ρ̄ ln(ρ̄/|ψ|²)). Dynamics: dX = (v + u)dt +
√(2ν)dW at ν ∈ {0 (dBB control), 0.1, 0.5, 1.0}; measure τ(ν, M) from
the exponential fit window.
**Gates.** G1: dBB control consistent with tier3's τ (10.2/4.9 for
M = 9/16) within 25% (different N and binning — state the comparison
honestly). G2 (equivariance control at every ν): equilibrium-born
ensemble stays at the floor. G3 (monotonicity): τ(ν) strictly
decreasing in ν for both M — the pre-registered direction; if
non-monotone, FAIL with diagnosis. G4 (quantification): report
τ(ν=0.5)/τ(dBB) for both M with bootstrap errors; and the small-ν
behavior (does ν → 0 recover dBB continuously — τ(0.1)/τ(0) ≥ 0.7?).
G5 (consequence paragraph, report-only): the corpus's III.D rates and
relic-window arithmetic are conservative if τ_Nelson < τ_dBB — state
which corpus claims move in which direction (equilibrium quality:
strengthened; relic non-equilibrium survival: weakened — the honest
double edge).
**Deliverables.** `M4/m4_nelson_relax.py`, `M4/m4_results.json`,
`M4/RESULTS.md`, H̄(t) figure.

## Execution
M1–M4 are independent → one parallel workflow fan-out; coordinator
adjudication in `TIER7_ADJUDICATION.md`; artifacts to branch
`claude/analyze-gum-po`. M3 and M4 must READ the prior artifacts they
extend (h21, tier3-born) before writing code — parametrization reuse
is a gate condition, not a suggestion.
