# S1 — P-acoustic adoption pilot: pricing cost (1) of the F-R14 sound cell [F-T8-S1]

**Phase:** Tier 8 (Phase S), ROADMAP_v10_ALTERNATIVES.md workstream S1.
**Date:** 2026-07-18.  **Code:** `s1_price.py` (pure sympy, deterministic,
nothing fitted; 14/14 internal checks PASS, ~11 s).
**Outputs:** `s1_results.json`, `s1_fig.png`.
**Sources:** `tier7-program/O3/RESULTS.md` §4 (the three printed costs of the
P-acoustic w = −1/4 cell; this workstream prices cost (1)); `o3_matrix.py`
(measure bookkeeping: det g₃ = (1+h)³, weight (det g₃)^w = (1+h)^{3w});
`corpus3/01-GUM-Omega-Paper-v4.3-ext.md` Sec. III.B–III.C (the flat-measure
Fisher → quantum-potential derivation, eqs. (3.2)–(3.4), Theorem III.1).

**Epistemic frame (binding).  Within-model; nothing here bears on nature.**
This workstream is **DECISION SUPPORT ONLY**: it turns O3's qualitative cost
(1) — "reopens Sec. III's flat-measure [DF] derivations (Fisher/Madelung,
Wallstrom, Born rates) on curved backgrounds" — into an exact symbolic price.
**No adoption is recommended; the F-R14 regrade language is unchanged.**

---

## 1. What is computed

The corpus's Theorem III.1 chain, in campaign units as printed in Sec. III.B:

- E_Q[ρ] = (ħ²/8m) I_F[ρ], I_F = ∫ |∇ρ|²/ρ d³x — eq. (3.2) (Fisher
  stiffness **ħ²/8m**);
- δE_Q/δρ = (ħ²/8m)[|∇ρ|²/ρ² − 2∇²ρ/ρ] — eq. (3.3);
- = −(ħ²/2m) ∇²√ρ/√ρ = U_Q, Bohm's quantum potential — eq. (3.4),
  Theorem III.1 [DF].

is (a) reproduced exactly, then (b) repeated on the curved 3-metric
g₃ = (1+h)δ_ij with quantization-measure weight (det g₃)^w = (1+h)^{3w}
(the o3_matrix.py bookkeeping), general symbolic w, to linear order in h
(h → εh(x); |∇h|² terms are O(ε²) and are dropped by the series expansion —
no truncation below O(h²)).  The **correction** is defined as

  C = (curved chemical potential) − U_cov,  U_cov = −(ħ²/2m) Δ_g√f/√f,

Δ_g the Laplace–Beltrami operator of g₃ — i.e. the failure of the printed
identity to close in its natural covariant form.

**Formalization amendment (coordinator-owned, printed per the standing
rule).** The corpus prints *no* curved Fisher functional; "repeat the
variation with measure weight w" therefore requires a formalization choice
for where (det g₃)^w enters.  Rather than silently picking one, **both**
natural readings are computed and gated together (the spec's single
"curved-measure variation" is delivered as a two-scheme robustness pair;
their conclusions agree):

- **Scheme A (covariant kinetic term, weighted quantization measure):**
  E_A[ρ] = (ħ²/8m) ∫ g₃^{ij} ∂_iρ ∂_jρ/ρ (det g₃)^w d³x, chemical potential
  conjugate to ρ under dμ_w = (det g₃)^w d³x (so w = 1/2 is the Riemannian
  measure).
- **Scheme B (the printed flat functional verbatim on the weighted
  density):** ρ = (det g₃)^w σ with σ scalar; E_B[σ] = (ħ²/8m) ∫ |∇ρ|²/ρ d³x
  exactly as printed in (3.2) (flat contraction, coordinate measure);
  chemical potential conjugate to σ; target U_cov in σ.

## 2. Gate S1-G1 — the flat identity, exact

Both printed steps verified as exact symbolic zeros (3-D, general ρ):
δE_Q/δρ equals the intermediate form (3.3) with coefficient ħ²/8m, and (3.3)
equals −(ħ²/2m)∇²√ρ/√ρ identically.  Theorem III.1's flat [DF] status is
reproduced, coefficient linkage ħ²/8m ↔ ħ²/2m exact.

## 3. Gate S1-G2 — the curved variation at general w (the centerpiece)

**Scheme A** (exact symbolic result, then linearized):

  C_A = (ħ²/8m) **(3 − 6w)** ∇h·∇ρ/ρ + O(h²).

- **Constant h:** C_A ≡ 0 *exactly, to all orders, for every w* — the
  measure factor cancels between the functional and the conjugate pairing;
  the surviving (1+h)^{−1} is precisely the covariant g^{ij} contraction.
  The chain survives constant conformal rescaling untouched.
- **Gradient h:** a single new structure, ∇h·∇ρ/ρ — a genuine h-gradient
  force term absent from the flat chain.

**Scheme B** (linearized; core ≡ |∇σ|²/σ² − 2∇²σ/σ, the flat U_Q form):

  C_B = (ħ²/8m)[ **(3w+1)** h·core + **(1−6w)** ∇h·∇σ/σ **− 6w** ∇²h ] + O(h²).

- **Constant h:** only the multiplicative piece (3w+1)h × (flat U_Q)
  survives — a w-dependent rescaling of the stiffness ħ²/8m, absorbable into
  ħ²/m only if h is position-independent; zero only at w = −1/3.
- **Gradient h:** two independent new structures (∇h·∇σ/σ and ∇²h) with
  incompatible zeros.

## 4. Gate S1-G3 — the price (cost (1) quantified)

**Cancellation, solved exactly.**

| scheme | cancelling w | status |
|---|---|---|
| A | **w = 1/2, unique** (the Riemannian measure; correction vanishes to *all orders* there) | **outside** the F-R14 rescue window (−5/18, −2/9); not among printed conventions {0, −1/4, −1/2} |
| B | **none** — the three structures' zeros are w = −1/3 (h·core), 1/6 (∇h·∇σ/σ), 0 (∇²h): simultaneous solve empty | no w cancels; in particular each printed w leaves ≥ 2 structures alive |

**No weight both rescues F-R14 and preserves the Fisher → Bohm identity:**
the only cancelling weight in either scheme, w = 1/2, lies outside the O3
rescue window — the geometry of the two requirements is disjoint.

**Correction at the printed values** (coefficients × ħ²/8m):

| w | scheme A: ∇h·∇ρ/ρ | scheme B: (h·core, ∇h·∇σ/σ, ∇²h) |
|---|---|---|
| 0 (flat convention) | 3 | (1, 1, 0) |
| **−1/4 (P-acoustic)** | **9/2** | **(1/4, 5/2, 3/2)** |
| −1/2 (other convention) | 6 | (−1/2, 4, 3) |

**The P-acoustic price, explicitly.**  At w = −1/4 the curved chain fails to
close with correction

  **C_A(w=−1/4) = (9ħ²/16m) ∇h·∇ρ/ρ**  (scheme A),
  C_B(w=−1/4) = (ħ²/8m)[¼ h·core + (5/2) ∇h·∇σ/σ + (3/2) ∇²h]  (scheme B),

i.e. a new osmotic-type force ∝ ∇h with O(1) coefficient (9/2 × the printed
stiffness in scheme A), plus, in scheme B, a σ-independent potential term
(3/2)(ħ²/8m)∇²h sourced purely by the metric.

**Effect on the [DF] label (honest statement).**  Theorem III.1 remains
[DF] **on flat backgrounds** — nothing here touches the h = 0 chain (and in
scheme A even spatially constant h is exactly harmless).  On curved
backgrounds with the P-acoustic weight, the Fisher → Bohm identity **does
not close as printed**: the h-gradient corrections above are not derived
anywhere in the corpus, so the curved extension of Sec. III.B–III.C
(Madelung closure inherits U_Q; Wallstrom and Born-rate machinery sit
downstream of the closure) would require a **new derivation** and is not
[DF] as it stands.  O3's cost (1) is thereby confirmed and priced: adopting
w = −1/4 buys the sound F-R14 cell at the cost of an open curved-background
re-derivation obligation with the explicit correction terms above.  (This
prices the reopening; it does not assert the re-derivation is impossible.)

## 5. Gates

| gate | requirement | measured | verdict |
|---|---|---|---|
| S1-G1 | flat identity reproduced symbolically, coefficient exact | (3.2)→(3.3)→(3.4) both steps exact symbolic zeros; stiffness ħ²/8m, Bohm ħ²/2m | **PASS** |
| S1-G2 | curved variation at linear order, general w, explicit symbolic expressions | C_A = (ħ²/8m)(3−6w)∇h·∇ρ/ρ; C_B = (ħ²/8m)[(3w+1)h·core + (1−6w)∇h·∇σ/σ − 6w∇²h]; constant-h and gradient-h cases separated; O(ε²) truncation only; numeric ε²-scaling check ratio 4.000 | **PASS** |
| S1-G3 | cost priced: cancelling w exact, or structure+coefficients at w=−1/4 + honest [DF] statement | both delivered: unique scheme-A cancellation w = 1/2 (exact, all orders), outside rescue window; no scheme-B cancellation; w=−1/4 coefficients 9/2 and (1/4, 5/2, 3/2); [DF] statement §4 | **PASS** |
| S1-G4 | decision-support-only; no adoption recommended; F-R14 regrade language unchanged | frame printed in header and §6; no recommendation anywhere; regrade language untouched | **PASS** |

## 6. The honest bottom line (Gate S1-G4)

**DECISION SUPPORT ONLY.  No adoption is recommended; the constitutive
choice is the authors' (charter language).**  The F-R14 regrade language is
unchanged: Theorem VI.1 remains OPEN-DEFECT as printed; (6.3) remains
demoted to order-of-magnitude; P-M3 remains premise-suspended.  What this
workstream adds is the exact price of O3's cost (1): the P-acoustic weight
w = −1/4 — the only computationally sound F-R14 cell — makes the Sec. III
Fisher/Madelung [DF] chain fail to close on curved backgrounds at linear
order in h, with correction (9ħ²/16m)∇h·∇ρ/ρ (scheme A) resp. three-term
structure (1/4, 5/2, 3/2)×(ħ²/8m) (scheme B), and **no** measure weight
cancels the corrections inside (or anywhere near) the rescue window.  The
two repairs h33 called "the minimum bill" therefore carry a third, now
quantified, IOU: a curved-background re-derivation of Sec. III.

## 7. Key numbers

| quantity | value |
|---|---|
| flat identity | exact ((3.2)→(3.3)→(3.4); stiffness ħ²/8m) |
| scheme A correction (general w) | (ħ²/8m)(3−6w) ∇h·∇ρ/ρ |
| scheme B correction (general w) | (ħ²/8m)[(3w+1)h·core + (1−6w)∇h·∇σ/σ − 6w∇²h] |
| cancelling w, scheme A | 1/2 exactly (unique; Riemannian measure; all orders) |
| cancelling w, scheme B | none (zeros −1/3, 1/6, 0 mutually incompatible) |
| rescue window vs cancellation | w = 1/2 ∉ (−5/18, −2/9): disjoint |
| correction at w = −1/4 | scheme A: **9ħ²/16m** on ∇h·∇ρ/ρ; scheme B: (1/4, 5/2, 3/2) × ħ²/8m |
| corrections at w = 0 / −1/2 (scheme A) | 3 / 6 × ħ²/8m |
| constant-h case | scheme A: zero exactly, all w, all orders; scheme B: (3w+1)h rescaling |
| checks | 14/14 PASS, ~11 s |

## 8. Caveats (honest limits)

1. **The formalization choice is ours, printed as an amendment.**  The corpus
   prints no curved Fisher functional; schemes A and B bracket the natural
   readings, and their shared conclusion (no cancelling w in the printed set
   or the rescue window; explicit O(1) corrections at w = −1/4) is the robust
   content.  A third reading engineered ad hoc to cancel is not excluded by
   this computation — but it would itself be a new constitutive postulate,
   i.e. it would add to, not reduce, the F-R14 bill.
2. **Linear order in h, conformal metric class** g₃ = (1+h)δ — inherited
   deliberately from h33/O3 (same probe class).  Anisotropic perturbations
   and O(h²) terms are not computed; they can only add structures, not
   restore the closed identity.
3. **Scope:** cost (1) names Fisher/Madelung, Wallstrom, and Born rates; this
   workstream prices the Fisher/Madelung link (the head of the chain, whose
   U_Q feeds Theorem III.2's closure).  Wallstrom (phase single-valuedness)
   and the Born-rate machinery sit downstream and are not separately
   re-derived here; the priced failure at the head already reopens them, as
   O3's cost language states.
4. Time-independent, single-particle setting, static h — matching Sec. III's
   own derivation scope.
5. This memo prices; it does not adjudicate.  Adoption is the authors' call.

## Deliverables

- `s1_price.py` — the computation (14/14 checks, deterministic, self-checking).
- `s1_results.json` — machine-readable: conventions, both schemes' symbolic
  corrections, cancellation solves, per-w tables, gates, check list.
- `s1_fig.png` — correction coefficients vs w, rescue window and printed
  conventions marked (scheme A zero at w = 1/2 visibly outside the window).
- `RESULTS.md` — this memo.  Filed as **F-T8-S1**.
