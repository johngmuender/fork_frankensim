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

---

## ADDENDUM — Phase N (pre-registered): v4.1 fold, F-R9 route (a), L7 continuum proof

### N1 — corpus3 v4.1 delta (the F-T7 fold)
Editor task under REVISION_CHARTER_v2 (marker ⟦T7⟧; payload authority
corpus3/T7_FOLD_PAYLOAD.md): bump the paper to v4.1-ext (git mv + all
version headers), insert ⟦T7⟧ blocks at Sec. 0 (third-phase record),
III.D (M4), IV.H.3 (M2), VI.D (M1 + F-T7-1 flag), Sec. X (T7 record),
Sec. XI (register update: WS4-M4c toy-discharged; F-R9 menu narrowed);
update README.manifest v5 version history + Replication Record (§8) +
append a Tier-6/7 section to ../REPLICATION_CAMPAIGN_STATUS.md
(additive). Verification: additive-only diff (headers exempt), payload-
sourced numbers, grades never rise.

### N2 — F-R9 route (a): the honest second-moment treatment
On h21's parametrization (mandatory reuse): compute the observable-level
consequence of including the second moment. Key physics to establish:
(i) for a macroscopic sample the total rate self-averages SPATIALLY to
N·⟨A²⟩ ⇒ enhancement factor R = 1e14–1e18 over the corpus's ⟨A⟩²-based
claim (staticity of the web on experiment timescales follows from h21's
own motional-narrowing pricing — cite it); (ii) rate-distribution
statistics across nuclei (burstiness: fraction of rate carried by the
top 1e-x of sites); (iii) propagate through the corpus's Majoron/P-ν4
phenomenology: with rate ∝ g²-class scaling, the honest rate = R × the
printed rate ⇒ effective sensitivity shift √R = 1e7–1e9 in g-equivalent
units; check against the corpus's own printed g-window [8e-10, 1.3e-8]
and the J3 battery re-run ("SHIFTS-HARMLESSLY at 7.8e-10") — does the
honest second moment break the battery, force a parameter retreat, or
strengthen a kill? Pre-registered decision frame: report which of the
three, with the arithmetic; no tuning. Gates: G1 h21 anchor
reproduction; G2 the enhancement + burstiness measured; G3 the
propagation table (printed rate → honest rate → battery consequence)
with an explicit verdict; G4 honest-limits paragraph (what route (a)
does NOT settle: the K3a insert's upstream physics; S1 untouched).

### N3 — L7 continuum proof of the cutoff zeros (T4-W5 residual)
Elevate L7b's sampled zeros toward theorem grade via a FIELD-LEVEL
sufficient condition: Theorem (crossing criterion): if
sup_z v_x(d, z, t) < 0 for all t ∈ (τ*, T], then no trajectory crosses
x = d in (τ*, T] (a crossing at time t requires v_x ≥ 0 at the crossing
point), hence the first-crossing density vanishes there — continuum
support ⊂ [0, τ*]. Execute on the L5′/L7b engine (n_y = +1): compute
v_x(d_near, z, t) on the full (z, t) grid; measure τ* = sup{t :
max_z v_x(d, z, t) ≥ 0}; margin δ(t) = −max_z v_x on (τ*, T].
Gates: G1 field quality (machine precision as before); G2 τ* exists,
stable ≤ 1% under Nx-doubling + dt-halving; G3 consistency:
τ* ≥ τ_max(sampled) = 5.130950 and the gap τ* − τ_max reported +
explained (the sampled max is ≤ the continuum edge); G4 margin: δ(t) ≥
δ_min > 0 on (τ* + ε, T] with δ_min ≫ discretization error bound
(quantify both) — the hypothesis of the theorem verified at
proof-grade-modulo-discretization; G5 discrimination control: same
computation at n_y = +0.25 must show max_z v_x ≥ 0 recurring out to
late times (no τ* below T) — the criterion distinguishes the cutoff
from the no-cutoff family member. Honest scope: this proves the
continuum support statement GIVEN the computed field (discretization
bounds printed); it upgrades the L7b kill-bin zeros from ensemble
statement to field statement; the fully analytic proof remains open.

---

## ADDENDUM — Phase O (pre-registered): flow-map proof, g_lim archaeology, F-R14 decision matrix

### O1 — the flow-map first-crossing proof (T4-W5 terminal shape, per F-T7-N3)
N3 proved the line criterion cannot work; the transport argument can.
**Theorem (backward-reachability criterion; state + prove):** trajectories
of the guidance ODE are unique both ways; if the backward trajectory
from every kill-bin detector point (d, z, t), z ∈ (0,1), t ∈ the kill
bins, satisfies max_{s<t} X_x(s) ≥ d (it already crossed the line
earlier), then NO first crossing occurs in the kill bins — the arrival
density there is exactly zero at continuum (field) level.
Execute on the validated L5′/L7b engine, n_y ∈ {+1, −1, +0.75, −0.75}:
backward-integrate from a fine grid over the kill bins (≥ 200 z × 40 t
per bin per n_y); classify each point: (a) PRE-CROSSED (max_{s<t} X_x ≥
d + margin), (b) ρ-floor (backward path enters ρ < ρ_floor near-node/
wall zones — report the ρ-measure it carries), (c) VIOLATION (backward
path reaches t = 0 inside supp ρ₀ without pre-crossing = a genuine
fresh arrival — would refute the zeros). Gates: G1 engine bars;
G2 zero VIOLATIONS at reference resolution AND under Nx-doubling +
dt-halving of the backward integration; G3 pre-crossing margin
quantified (min over class-(a) of max X_x − d) ≫ trajectory
discretization error (quantified by refinement differences);
G4 ρ-floor accounting: total |φ₀|²-measure reachable through class-(b)
≤ 10⁻³ (else PARTIAL with the honest bound); G5 discrimination:
n_y = +0.25 backward runs from its POPULATED late bins must show a
substantial fresh-arrival (class-c-like, legitimately arriving)
fraction — the criterion distinguishes. Verdict target: the kill-bin
zeros elevated to FIELD/FLOW statement (proof-grade modulo
discretization); the fully analytic proof remains open but its shape
is now demonstrated. File as F-T7-O1.

### O2 — g_lim archaeology (the N2 blocking number)
Exhaustive archive sweep (corpus/ all 104 files, corpus2/, corpus3/,
theory-audit/, DISCHARGE_PACKAGE/, tier*/) for ANY printed: Majoron-mode
0νββ bound, T_1/2 limit, experimental rate number, P-ν4 detectability
line or rate normalization, LEGEND/nEXO sensitivity figures, or any
g-constraint of any channel (catalog all, with verbatim quotes + file
paths). Decision rule (pre-registered): if a Majoron-mode rate-level
bound (g_lim or T_1/2-equivalent) is found → execute the F-T7-N2
conditional-kill arithmetic against it and report which way it fires;
if only non-rate g-constraints exist → catalog them and state precisely
why they do not adjudicate the honest-rate question (the observable
distinction); if nothing → certify ABSENT with the search protocol
printed (patterns, file counts) so the certification is reproducible.
Gates: G1 sweep coverage (file count + pattern list printed; ≥ 2
independent pattern families); G2 the catalog (every hit quoted);
G3 the verdict per decision rule. Deliverable: memo
tier7-program/O2/GLIM_ARCHAEOLOGY.md (+ results JSON). File as F-T7-O2.

### O3 — F-R14 rescue decision matrix (decision support for the authors)
Implement h33's master formula 𝔞₁(p,q) = −(5p+q)/(12(2p+q)) and the
corpus sector table (B2 doublet p = q − 1 structural; Q̃-sectors
q = 3/2 group-protected; knot band-edges with the supertrace-signed
count per h25 K3). Compute, for each textually available reading
(frame/soldered, covector, vector, cone-only) AND across the
P-acoustic window w ∈ (−5/18, −2/9) (incl. the distinguished
w = −1/4): per-sector weights w_s, the sign of Σw_s (1/16πG), whether
c_GW² is a convex combination (hull check with the printed sector
cones), and the (6.3) margin status. Gates: G1 reproduce h33's five
printed anchors exactly (−1/3, −2/15, −5/33, −1/12; +1/6 at w = −1/4
on q = −3p); G2 the decision matrix (readings × outcomes) with every
cell computed, not asserted; G3 the P-acoustic pricing quantified:
window width vs distance to natural conventions (w = 0, −1/2) as a
tuning fraction; the supertrace erratum's effect on the knot-sector
row shown both ways; G4 the honest bottom line: this is DECISION
SUPPORT — the constitutive choice is the authors' (charter language);
no reading is recommended, each is priced. File as F-T7-O3.

---

## ADDENDUM — Phase P (pre-registered): v4.2 delta + T4-W5 certification pilot

### P1 — corpus3 v4.2 delta (the N/O fold)
Editor task under REVISION_CHARTER_v2 (marker ⟦T7⟧; payload authority =
T7_FOLD_PAYLOAD.md ADDENDUM Q6–Q9): git mv the paper to
01-GUM-Omega-Paper-v4.2-ext.md; bump all six version headers (v4.1 →
v4.2, parenthetical unchanged); add a ⟦T7⟧ v4.2 REVISION NOTE after the
v4.1 note; insert ⟦T7⟧ blocks: Sec. 0 (Q9 phases record), the VII.I/
IX.E Majoron area (Q6 — honest rate + certified-absent + the armed
conditional kill), VIII.F/IX.B T4-W5 box area (Q7 — field/flow-grade
zeros; residual narrowed to certified/analytic form), VI.C dispute-box/
App. E.6 area (Q8 — the decision matrix, one sound cell, priced),
Sec. X (Q9), Sec. XI (register: T4-W5 residual narrowed; F-R9 terminal
menu; F-R14 now carries the computed matrix pointer). Update
README.manifest v5 (census cell + version history) and Replication
Record (§9). Verification identical to N1's (additive-only diff vs git
HEAD v4.1; marker counts; ≥ 15 payload numbers verbatim; ⟦rev⟧/⟦T6⟧/
prior-⟦T7⟧ preserved).

### P2 — T4-W5 certification pilot (exact-field, high-precision backward paths)
Upgrade O1 one rung toward certified: the engine's field is a FINITE
trigonometric sum (analytic sine-mode factors in z; band-limited FFT
representation in x), so ψ, ρ, j — hence v — can be evaluated EXACTLY
(to fp rounding) at ARBITRARY points, eliminating O1's dominant error
source (grid interpolation). Implement exact band-limited evaluation
(direct mode sums at query points; verify against the grid engine at
grid nodes to ~1e-13) and re-run the backward classification for a
STRATIFIED subset: all 16 box-exit/wall-adjacent paths + the 100
smallest-margin class-(a) paths + 100 random class-(a) paths per
|n_y| ∈ {1, 0.75} (≥ 432 paths total), with (i) high-order integration
(DOP853 rtol 1e-12) and (ii) a precision ladder (float64 vs mpmath
50-digit on ≥ 20 paths incl. the worst margins). Gates: G1 exact-vs-
grid field agreement at nodes ≤ 1e-12 rel; G2 every re-run path keeps
its O1 classification (pre-crossed stays pre-crossed; margins agree
with O1 within O1's own error estimate); G3 precision ladder: float64
vs 50-digit path endpoints agree ≤ 1e-8 (integration-truncation
bounded), margins stable to ≥ 6 digits on the ladder subset;
G4 the worst-margin path's pre-crossing certified at ladder grade with
the margin printed; G5 honest scope: this is precision-certification
(exact field + converged integration), NOT formal interval arithmetic;
state what a formal certification would still require (validated ODE
enclosures) and that the analytic proof remains the terminal open.
File as F-T7-P2.
