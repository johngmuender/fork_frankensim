# T1 — Theory audit of the open NR gates
## NR-A2b (¾-Casimir) · NR-D2 (collision kill) · NR-D1b (core-energy floor) · NR-K2b/K3a (two-scale pair)

**Auditor:** theory-audit subagent, replication campaign Phase H (T1).
**Date:** 2026-07-16. **Sources:** `analysis/gum-sandbox/corpus/NR-A2b-…`, `NR-D2-…`,
`NR-D1b-…`, `NR-K2b-K3a-…` (all v0_1), plus upstream corpus instruments
(WS-A1/NR-A1, WS-A4, WS-A v0.2, WS-D3, WS-K delta pack) and the campaign record
(tier0-gauntlet/RESULTS.md Groups F/G/H = F-R1; tier2-closure/axi_FR4_NOTE.md = F-R4).

**Epistemic frame (binding).** Everything below is WITHIN-MODEL: I audit whether the
corpus's printed arguments are valid given their own printed premises, never whether
the physics is true of nature. Classification scheme: **SOUND** (valid given premises) /
**GAP** (a missing step, stated precisely) / **DEFECT-CANDIDATE** (the step fails; the
counterexample or derivation is shown) / **NOT-AUDITABLE-FROM-TEXT** (premises or
definitions unrecoverable; what is missing is stated). Where a small computation
settles a step I ran it (python3; numbers printed inline). These four instruments are
the archive's own list of "what external hands must close"; each one *self-declares*
its replication-gate inputs, and part of this audit is checking whether the
self-declared gate list is complete. Verdict in one line: **it is not complete in two
places** (NR-A2b's ¾ calibration; NR-K3a's fluctuation step), **and it is
direction-safe in the one place that matters most** (NR-D2's kill).

---

## 1. NR-A2b — the transport-profile theorem (the ¾-Casimir)

### 1.1 Argument skeleton

- **A2b-1 (w1, null transport):**
  P1: charge transport lives in the strictly gapless channel; exact masslessness
  leaves "no sub-luminal branch for a stationary transported structure" ⟹ |v| = c on
  every occupied streamline.
  P2: stationarity at the single entrained frequency ⟹ the pattern co-rotates:
  **v = ωr_⊥ azimuthally**.
  C: |v| = c ∧ v = ωr_⊥ ⟹ r_⊥ = c/ω = ƛ on every streamline ⟹ ⟨r_⊥²⟩_j = ƛ² exactly.
  Route B is cited as independent (NR-A1's flux counting, I = qω/2π, Φ_eff = πBƛ²).
- **A2b-2 (w2, the Casimir–projection split):**
  P1: at point level the displacement about the guidance center is operator-valued;
  Wigner–Eckart leaves only Ŝ to build a vector from.
  P2 (projection face): 𝔐 = μ_B(L+2S)/ħ; spin eigenvalues ±μ_B; "the stretched-axis
  reading, geometric radius ƛ"; "in natural units the projection face carries
  (2s)² = 1 at s = ½".
  P3 (Casimir face): ⟨δr̂²⟩ is a rotational invariant; the only invariant is the
  Casimir; "calibrating the operator loop on the stretched axis (its projection
  extent = ƛ, from (a)/§1)" gives **⟨r_⊥²⟩_ρ/ƛ² = s(s+1)|_{s=½} = ¾ exactly**.
  Caveat F-A2b-2′ (printed by the corpus): at s = ½ the candidate laws s(s+1) and
  (s+1)/2 coincide; the operator law is NOT selected; harmless because Thm 12.4
  forbids elementary s ≥ 1.
- **w3 (seal opening):** the classical shell law δ𝒥 = −c_sh·ε^{2/3}, c_sh = O(1),
  read as a moment correction, violates WS-A4's demand table (|c| ≤ 1.4×10⁻¹⁰ at
  p = 2/3) by ~10 orders; rescue: the moment is assembled from exact structures
  (n, ħ/L = 2, c, ƛ) so the ε^{2/3} content "lands in the Casimir face … and cancels
  identically in the projection face"; surviving corrections ≲ 3×10⁻²¹.
- **w4:** re-reading of the frozen classical profile class as the operator's
  silhouette; annotation only, no new mathematics.

### 1.2 Verdicts

**A2b-1 — GAP.** The step from P1∧P2 to r_⊥ = ƛ *pointwise* needs an additional
premise nowhere stated or proved: **the occupied streamlines have no axial (or
radial-circulation) velocity component.** A stationary single-frequency pattern can
carry helical transport: v = ωr_⊥ φ̂ + v_z ẑ with the pattern stationary in the
co-rotating frame; then |v| = c gives ω²r_⊥² + v_z² = c², i.e. **r_⊥ ≤ ƛ with equality
only for purely azimuthal flow**, and ⟨r_⊥²⟩_j < ƛ² whenever any streamline is
helical. The corpus's own earlier version of this argument (NR-A1 §2.4) was graded
honestly — "the exactness of (iii) beyond the ideal limit is not yet a theorem",
δ𝒥 sealed — and NR-A2b's upgrade to "exact" is carried entirely by the unproven
azimuthal-only clause. Note also that "Route B (independent)" is overstated by the
corpus's own record: NR-A1 §3.3 states in print that "the routes share the cylinder
lemma," so the two routes are not independent on exactly the load-bearing point.
*What would close it:* a lemma that stationarity + the channel rule forbid axial
components on occupied streamlines (e.g., that any v_z ≠ 0 breaks single-valuedness
or the entrainment identity), or an explicit variational demonstration that the
stationary solution's current support is planar.

**A2b-2(b), the ¾ — DEFECT-CANDIDATE (calibration inconsistency).** Formalize the
"operator loop": the displacement is a vector operator built from Ŝ, so δr̂ = c·Ŝ for
some scalar c (this is exactly what Wigner–Eckart licenses). Then for **any** c,

  ⟨δr̂²⟩ / (max ⟨δr̂_z⟩)² = s(s+1)/s² = **3** at s = ½   (computed; = 3.000),

not ¾. The two determinate readings of the document's calibration both fail to give ¾:
(i) calibrate the stretched projection to the §1 ring **radius** ƛ (the document's own
cross-reference: "its projection extent = ƛ, from (a)/§1", where §1 puts every
streamline at r_⊥ = ƛ): then ⟨δr̂²⟩ = 3ƛ²; (ii) calibrate ⟨δr̂²⟩ itself to the ring:
then the ratio is 1, and there is no ¾ either. The printed ¾ is obtained only by
comparing **two differently normalized operators across the two faces**: the
projection face is computed with 2Ŝ/ħ ("(2s)² = 1" — the g = 2-dressed winding
counter), while the Casimir face is computed with Ŝ/ħ (s(s+1) = ¾). Equivalently, ¾
requires reading "projection extent = ƛ" as *peak-to-peak* (eigenvalues ±ƛ/2), which
contradicts §1's radius-ƛ null ring (peak-to-peak 2ƛ). No single linear-in-Ŝ
displacement operator consistent with §1's geometry yields ¾ƛ².
  *Uncertainty, printed:* the underlying operator construction (VII.F FR-quantized
moduli) is not reproduced in this instrument; if the corpus elsewhere defines the
loop's calibration so that the stretched projection is ƛ/2 *and* reconciles that with
§1's r_⊥ = ƛ, the verdict demotes to GAP (unstated, load-bearing convention). As the
instrument is written, the algebra does not close. Context that raises the bar: the ¾
was **inverted from hydrogen data first** (WS-A v0.2 §2 row: "¾, inverted
[⬛-inversion; kill attached to NR-A2b]"), so the derivation is a postdiction whose
one free convention lands exactly on the known target — the calibration step is
therefore precisely where a within-model auditor must demand a theorem, and it isn't
one. This is the audit's sharpest new finding: **the corpus's self-declared gate for
NR-A2b ("the operator construction from the FR quantization proper") is the right
gate, but the corpus grades the ¾ as "exact, as algebra" when the algebra as printed
gives 3, not ¾.**

**A2b-2′ (the degeneracy caveat) — SOUND.** s(s+1) = (s+1)/2 ⟺ s = ½ (checked:
0.750 = 0.750 at ½; 2 vs 1 at s = 1; 3.75 vs 1.25 at 3/2). Given Thm 12.4 (spin
selection — replicated by the campaign on the credit side), the law is only ever
evaluated at s = ½ and the ambiguity is indeed empty *within the model*. The caveat's
logic is valid; note only that it concerns *which* invariant law, not the calibration
defect above, which it does not cure.

**w3 — GAP (conditional argument, inherits A2b-1).** The "protection quartet" has
real logical force *conditional on A2b-1 being exact*: if all current support sits at
r_⊥ = ƛ identically (zero radial spread of j), then a shell-thickness correction
cannot move the current-weighted flux loop, and the ε^{2/3} content can only appear
in the charge-sampling (ρ) distribution. But (a) the "cancels identically in the
projection face" clause is asserted, not derived — no computation decomposes the
sealed δ𝒥_shell law into faces; and (b) the condition is A2b-1's exact null ring,
which is itself a GAP (above). The demand-table arithmetic checks: WS-A4's p = 2/3 row
reads |c| ≤ 1.4×10⁻¹⁰ (verified in the file), so c_sh = O(1) is a ~10-order violation
as stated, and the campaign's tier-0 Group G certified the WS-A4 tables. The imported
bounds 3×10⁻²¹ and 10⁻²⁷ are not derivable from this text (inputs from the lattice/
cone and transduction instruments) — arithmetic-plausible, not audited here.

**w4 — SOUND-as-annotation.** No new mathematical content to audit; it re-labels the
classical exclusions, which stand independently.

---

## 2. NR-D2 — the collision problem (the fired kill)

### 2.1 Argument skeleton

P1 (classicality): λ_dB ∼ 10² fm ≪ p = 4 μm ⟹ classical encounter at tether scale.
P2 (criterion): if V_int > KE = ½m_𝔪v² ≈ 2×10⁻⁴f (v = 10⁻³c, m = 16f), the objects
bounce at tether contact; σ ≈ π(2p)² geometric.
P3 (naive wall): director mismatch cost K/p² × p³ = f²p per crossing, K = f².
P4 (screening): the mismatch "**is** SDiff-reachable (a smooth interpolating
reorientation lives on the soft manifold)" ⟹ cost is ε-lifted: V_int ≈ ε·f²p =
(7.6×10⁷ε)·f.
P5 (band, imported): ε ∈ [1.0×10⁻⁵, 6×10⁻⁴], frozen by the frustration sector.
C1 (F-D2b-1): transparency needs ε ≤ 3×10⁻¹²; shortfall ×3×10⁶ to ×10¹⁰·⁵ ⟹ hard
μm sphere; σ/m ≈ 5×10¹⁸ cm²/g vs Bullet ≲ 1 ⟹ KILL, over by 18+ orders.
Escape hunt (α)–(ε) closes attractive overlap, velocity structure, unlifted soft
directions, tether shedding, and the LC anchor. Channel 2 (web pinning) is filed as
over-determination.

### 2.2 Arithmetic spine, recomputed (independent of tier-0, agreeing with it)

| Quantity | Corpus | Recomputed | Note |
|---|---|---|---|
| f·p | 7.6×10⁷ | 7.602×10⁷ | PASS |
| KE(v=10⁻³c, m=16f) | 2×10⁻⁴ f | **8.00×10⁻⁶ f** | **F-R1 confirmed: ×25; printed value = v = 5×10⁻³c** |
| ε bound (printed KE) | 3×10⁻¹² | 2.63×10⁻¹² | PASS as rounding |
| ε bound (true KE) | — | 1.05×10⁻¹³ | kill fires ~25× harder |
| shortfall, favorable corner | ×3×10⁶ | 3.8×10⁶ (printed KE); 9.5×10⁷ (true) | PASS |
| shortfall, far corner (+Bullet ×10²) | ×10¹⁰·⁵ | 2.3×10¹⁰ = 10¹⁰·⁴ | PASS |
| σ/m | ≈5×10¹⁸ cm²/g | 5.0×10¹⁸ at m = 60f; 1.9×10¹⁹ at m = 16f; 4.7×10¹⁸ at 64f (tier-0's mid-band) | PASS in order; corpus figure corresponds to a mid-band mass, not the m = 16f floor used for KE — a benign mixed-normalization inside one section |
| λ_dB | "∼10² fm" | **3.29×10³ fm** (ħ/mv, m = 16f, v = 10⁻³c); 3.3×10² fm at v = 10⁻²c | **minor slip, new**: the printed number matches the Bullet velocity, not the stated v; anti-conservative in direction (flatters classicality) but harmless in magnitude — λ_dB/p = 8×10⁻⁷ ≪ 1 regardless |

### 2.3 Verdicts

**The hard-wall criterion + kill arithmetic (F-D2b-1) — SOUND given premises.**
Given P1 (true with six orders to spare even after correcting the λ_dB print), P4-or-
its-negation (see below), and P5 (an import this audit takes as frozen), the
comparison V_int vs KE and the geometric σ/m are valid and the margins are as printed
or larger. The known F-R1 defect (×25 in KE) is **conservative**: every reading
tightens the kill. The one-line merger closure (attractive ⟹ same geometric σ) is
valid as an inequality argument.

**P4, the SDiff-reachability/ε-lift step — GAP, and this is where F-R4 bears.**
The reachability claim is a single parenthesis; no construction of the interpolating
reorientation on the soft manifold is given, and the campaign's tier-0 addendum
already noted the two physics gates (reachability; hard-wall applicability) are
argument-level. What F-R4 adds, precisely, in both directions:

1. **Against the mechanism as stated (methodological precedent + a lemma-shaped
   obstruction).** F-R4 (axi_FR4_NOTE.md) proved, for the flagship's functional, that
   every value-function integral (potential, sextic via b∘T·det DT = b∘T, inertia) is
   **frozen on the compositional SDiff orbit**: pullback Q∘T rearranges field values
   but cannot change their distribution. The G.5 defect it exposed was exactly a
   claimed "SDiff images approach X" where no pullback can redistribute a value
   function. NR-D2's P4 is the same move-class: a *reconciliation* of two mismatched
   double-twist textures must **change director values** in the overlap region, and a
   changed value distribution is not reachable by composition with any
   volume-preserving diffeomorphism. So if "the soft manifold" means the SDiff orbit
   in F-R4's (compositional) sense, the reconciled configuration is strictly
   off-orbit and the ε-lifted cost is a modeling estimate, not a consequence of
   reachability. (Caveat, printed: F-R4's lemma is proven for the flagship sector
   conventions; NR-D2's medium is the Frank-elastic chiral sector with its own
   frustration-manifold ε-lifting — the transfer is structural, not literal. The
   corpus never defines the murk-sector soft manifold precisely enough to decide;
   that definition is itself missing text.)
2. **For the kill (direction-safety).** The gap is harmless to the verdict, twice
   over. If reachability fails, the cost reverts toward the unscreened f²p — the kill
   fires ~1/ε ≈ 10³·²–10⁵ harder. If it holds, V_int ≈ ε f²p ≥ 7.6×10² f versus
   KE = 8×10⁻⁶ f — still eight orders short of transparency at the corpus's most
   favorable corner. And if the murk-sector soft manifold *is* SDiff-like, F-R4's
   frozen-value lemma **rigorously closes escape (γ)** — an exactly-flat direction
   cannot perform a value-changing reorientation, so no zero-cost reconciliation
   channel can exist — which is stronger than NR-D2's own phenomenological closure
   ("a zero-cost channel would have falsified VII.I"). Net: **F-R4 makes the kill
   more robust while making the corpus's stated screening *mechanism* less
   theorem-like.** Every disposition of the gap lands on KILL.

**Channel 2 (the web) — NOT-AUDITABLE-FROM-TEXT.** The per-crossing reconnection
cost "∼ε·f-class" is not derived anywhere in the instrument or its cited upstream;
the comparison to KE is arithmetic-consistent (ε·f ∈ [10⁻⁵, 6×10⁻⁴]f vs
KE = 8×10⁻⁶f, so "comparable to pinning" holds and strengthens post-F-R1), but the
coefficient's provenance is unrecoverable. The corpus itself files this channel as
moot-but-printed; the classification matches that self-assessment.

**The escape hunt (α, β, δ) — SOUND given premises** (each is a one-step inequality
or a citation of NR-D1's Hessian/size-minimum); **(ε) NOT-AUDITABLE-FROM-TEXT** (an
[IM] literature anchor, consistent in direction, not checkable here).

---

## 3. NR-D1b — the core-energy floor

### 3.1 Argument skeleton

E_core = E₁ + E₂ + E₃ with signs tracked:
(i) E₁ = c_h·K·r_geo, c_h ∈ [8π/3, 8π] [IM], r_geo ≥ r_min ≈ 1/f (the flagged
assumption) ⟹ E₁/f ≥ 8.4;
(ii) E₂ (soft-rerouted far field) ∼ c·ε·(f·p) ∈ O(10²)–O(10⁴): "large, positive,
wildly uncertain";
(iii) binding/overlap corrections "O(1)-fractional of E₁ … cannot flip the sign of
the total".
F-D1b-1: "Every contribution to E_core is positive; the bare near-core term alone
gives 𝒢_c ≥ 8" ⟹ clear of eternal threshold 5.6 by ×1.5; S = π·64 ≈ 200 vs 92 ⟹
W-eternal. F-D1b-2: a one-sided bound answers the three-fate question.

### 3.2 Arithmetic (verified; agrees with tier-0 Group H)

8π/3 = 8.378, 8π = 25.13 (corpus [8.4, 25] PASS); 8.378/5.6 = 1.496 ("×1.5" PASS);
π·8² = 201.1 vs 92, ratio 2.19 ("more than ×2" PASS); exponent margin 109 → "e^{100+}"
PASS (tier-0 adds that the G→8 rounding is conservative-direction). E₂ band:
ε·f·p = 7.6×10² (ε = 10⁻⁵) to 4.6×10⁴ (ε = 6×10⁻⁴), matching "O(10²)–O(10⁴)".

### 3.3 Verdicts

**F-D1b-1 — GAP, in three named pieces, one of them NOT self-declared.**
- (a) *c_h,min = 8π/3* — imported [IM]; the corpus flags it as replication-gate
  input 1. Note a checkable subtlety for the follow-on: in one-constant elasticity the
  radial and hyperbolic hedgehogs are **degenerate** (n_hyp = M·r̂ with M ∈ O(3), and
  |∇(Mr̂)|² = |∇r̂|², giving c_h = 8π for both under E = c_h K r); the [8π/3, 8π] band
  therefore encodes elastic anisotropy and/or interior escape structure that the text
  does not specify. The *bound* uses the band's minimum, which is the safe end — but
  the band itself is unrecoverable from the text.
- (b) *r_min ≈ 1/f* — the corpus's own flagged single load-bearing clause
  (gate input 2), with the correct fail-soft note (a sub-core absorption mechanism
  reopens W-open, not W-decayed). The flag is honest; the logic of the fail-soft
  clause is valid.
- (c) **The binding-correction clause — the unflagged gap.** §1(iii) concedes the
  binding/overlap terms may be negative ("cannot flip the sign of the total" — a
  sign statement only), yet F-D1b-1's first clause asserts "every contribution …
  positive," and the floor is taken as E₁ alone. With clearance ×1.496, the eternal
  verdict survives only if the net negative correction is **< 33% of E₁**
  (1 − 5.6/8.378 = 0.332, computed) — and "O(1)-fractional" is exactly the regime in
  which 33% is not a theorem. The Hessian-minimum citation bounds nothing
  quantitatively (a minimum can sit at any depth). The rescue inside the registered
  ε-band is E₂ ≥ 7.6×10²·c, which dwarfs the threshold — but the corpus itself says
  the full 𝒢_c spans down to O(10) "across the ε-tail's uncertainty," so at the
  tail's bottom the ×1.5 clearance is the whole margin and the binding fraction is
  live. *Precise missing step:* a quantitative bound |E₃| ≤ η·E₁ with η < 1/3, or a
  floor on E₂ valid across the whole tail.

**F-D1b-2 — SOUND conditional on F-D1b-1.** The logic "a one-sided bound collapses a
three-window question when every point of the uncertain range is on one side" is
valid; it inherits (a)–(c) above and nothing else. Consequence (5)'s arithmetic
(m ≳ 2E_core ⟹ 𝒢_m ≳ 16; σ/m eased ×2–8) is consistent with the floor.

---

## 4. NR-K2b / NR-K3a — the two-scale pair

### 4.1 NR-K2b (factorization) — skeleton and verdict

P1 (imported, WS-K2): the weak vertex couples through V2, a single-integer-winding
line reconnection — a line event, stratum-local. C: A(Q) = [line vertex](Q;
stratum-scale only) × Z_shadow^{1/2}(Q-independent) — so ν-DIS sees a point and the
μm texture hides in the walled normalization.

**Verdict: GAP.** The factorization is an ansatz with a mechanism-shaped
justification, not a derivation: nothing in the text shows the in/out texture overlap
Z is Q-independent (the outgoing lepton's μm texture is boosted by Q-dependent
kinematics; why the overlap of differently-boosted textures is constant in Q is
exactly the content a factorization theorem would have to prove — compare the labor a
QCD factorization proof does). The WS-A precedent is an analogy, and Route B's
counter-channel checks are census assertions (⬛-anchored, not derivable here). The
corpus itself names "the factorization theorem in the line-dynamics formalism" as the
replication window — i.e., the corpus's self-declared gate is exactly this gap, and
the self-declaration is accurate and complete for K2b.

### 4.2 NR-K3a (zero-momentum theorem) — skeleton and verdict

Step 1: the Majorana insertion is line-supported (m(z) on the web); the physical
propagating m_ν reads the same insertion at k ≈ 0 with O(1) engagement (the measured
mass calibrates m̄). Step 2: in the second-order amplitude the insertion transfers
k = q₁ − q₂ to the background; external kinematics dominate at k = 0;
**m̃(k) = m̄·δ³(k) + fluctuations at k ∼ 1/L**; hence the amplitude is governed by the
volume-averaged m̄ — "line-support and bulk-support are amplitude-identical in the
mean"; fluctuations are "(fm/μm)-class, negligible"; fork V-K1 → branch (a).

**Step 1 — SOUND given premises.** For a delocalized propagating state over a
statistically homogeneous background, the first-order mass eigenvalue is the spatial
average of the mixing density — standard perturbation theory; the "one insertion, two
momenta" identification is definitionally coherent within the framework.

**Step 2 — DEFECT-CANDIDATE, two-part, with the derivation:**

*(i) The printed spectral formula is false for line-supported m.* A density supported
on sparse tubes (radius a, spacing L, a ≪ L) does not have power only "at k ∼ 1/L":
its power spectrum extends flat-ish out to k ∼ 1/a. Sparse support is broadband by
construction — the δ³(k) + (k∼1/L)-fluctuations form printed in §2 describes a
*smooth* density with long-wavelength ripples, i.e., it assumes away exactly the
line-support the theorem is about.

*(ii) The observable is quadratic in the amplitude, and the mean-amplitude step drops
the dominant term.* The event amplitude for a nucleus at x₀ is
A(x₀) = ∫d³y m(y)·K(x₀ − y) with K the double-propagator kernel of range
r ∼ 1/q̄ ∼ 2 fm (q̄ ∼ 100 MeV). The theorem computes ⟨A⟩ ∝ m̄ (true). But the decay
rate of a sample is Σᵢ|A(xᵢ)|² ∝ ⟨|A|²⟩ over the nucleus–web geometry, and for
kernel range ≪ tube spacing this is variance-dominated:
  ⟨|A|²⟩/|⟨A⟩|² ≈ 1/φ, φ = πa²/L² (tube volume fraction) when a ≫ r;
  ≈ L²/(π²r²) when a ≪ r.
Computed at the corpus's own scales (L = μm; a = 1/f = 52.6 fm; r = 2 fm):
φ = 8.7×10⁻¹⁵, **enhancement ∼10¹⁴**; thin-tube regime gives ∼2.5×10¹⁶. In plain
terms: the φ-fraction of nuclei sitting on tubes each carry amplitude m̂ = m̄/φ, so
the sample rate is φ·m̂² = m̄²/φ, not m̄². The fluctuation term the theorem calls
"(fm/μm)-class, negligible" is, on this estimate, **larger than the retained term by
fourteen-plus orders.** Notably this fails in the *opposite direction* from the
audit's branch-(b) (collapse ×10⁻¹⁸): the second moment says *enhancement*, i.e.,
**both** the fork's original branches may be wrong, and the four-cell table's collapse
to the top row (§3) is unsupported either way.

*Uncertainty, printed at full size.* This is an estimate-backed challenge, not a
proven defect: it assumes (1) a static (or slowly varying) web on the decay's
virtuality timescale, (2) the insertion linear in the local m(y) (the corpus's own
model), (3) kernel range set by the standard 0νββ virtuality, (4) uniform m̂ on the
support. A corpus rescue would need one of: a demonstrated smooth (bulk) component of
m(y) carrying most of m̄; motional narrowing (web dynamics fast enough to average m
within a single event's coherence time); or saturation/nonlinearity of the insertion
on tubes. None of these is in the text; the text's own caveat ("fluctuation
corrections estimated, not derived") marks the joint but grades it negligible without
the second-moment check. **The corpus's self-declared gate ("the zero-momentum
theorem with the fluctuation bound") names the right theorem but the instrument's
verdict (branch (a), adverse cells unreached) already leans on the bound it does not
have.** Downstream: §3's four-cell collapse and the basis-review stand-down are
conditional on Step 2 and should be graded conditional, not landed.

---

## 5. Cross-cutting findings

1. **The corpus's self-gating is mostly accurate, with two omissions.** Of the eight
   load-bearing steps audited, six are exactly the steps the instruments themselves
   flag as replication-gate inputs (A2b's operator construction; D2's reachability
   and hard-wall clauses; D1b's c_h and r_min; K2b's factorization; K3a's fluctuation
   bound — flagged but pre-graded "negligible"). The two defects the gates do NOT
   name: **A2b-2's calibration inconsistency** (the ¾ as printed is not the algebra
   of any single operator consistent with §1) and **K3a's second-moment omission**
   (flagged as a caveat but the verdict is landed as if bounded). D1b additionally
   has one unflagged gap (the binding fraction vs the ×1.5 clearance).
2. **Direction analysis.** NR-D2 is the only gate whose conclusion is robust to every
   disposition of its gaps (both P4 dispositions tighten or preserve the kill; F-R1
   tightens ×25; the λ_dB slip is 6 orders inside its margin). NR-D1b's conclusion
   survives its flagged gaps by design (fail-soft) but is exposed to the unflagged
   binding clause only at the ε-tail's bottom. NR-A2b's and NR-K3a's headline
   conclusions rest directly on their defect-candidate steps.
3. **F-R4's bearing on NR-D2, stated once, both directions:** it undermines the
   *screening mechanism's* theorem-status (reconciliation is value-changing, hence
   off the compositional SDiff orbit — the same defect-class F-R4 proved against
   G.5), while *strengthening the kill* (no exactly-flat channel can reorient values,
   so the ε-lifted scale is a floor, closing escape (γ) more rigorously than the
   corpus's own argument). Transfer caveat: F-R4's lemma is proven in the flagship's
   sector conventions; the murk sector's soft manifold is never defined precisely
   enough in the text to make the transfer literal.
4. **Two new minor numeric findings** (F-R-ledger candidates, both sub-severity of
   F-R1): NR-D2's λ_dB print ("∼10² fm" vs computed 3.3×10³ fm at the stated v — the
   printed value matches v = 10⁻²c) and NR-D2's σ/m mass normalization (5×10¹⁸
   corresponds to m ≈ 60f, not the m = 16f used in the same section's KE; benign,
   order unchanged).

## 6. What a numeric follow-on could check

1. **A2b ¾-calibration exhaust (cheap, symbolic/numeric):** enumerate all calibration
   conventions δr̂ = cŜ with c fixed by each of {max projection = ƛ, ƛ/2;
   ⟨δr̂²⟩ = ƛ²; moment face via (L+2S)}, and tabulate ⟨r²⟩_ρ/ƛ². Confirms/refutes
   that no single convention consistent with §1 yields ¾. (Partially done here:
   3, 1, ¾-only-via-mixed-normalization.)
2. **A2b-1 helical counterexample (cheap):** construct an explicit stationary
   single-frequency null vector field with v_z ≠ 0 on a compact support satisfying
   the instrument's stated premises (i)–(ii), exhibiting ⟨r_⊥²⟩_j < ƛ²; this makes
   the GAP a demonstrated underdetermination.
3. **NR-D2 reconciliation-field solve (moderate):** on fs-gum-field/statics, prepare
   two overlapping director textures with mismatched double-twist axes and relax
   under (a) full functional, (b) SDiff-projected moves only; measure the actual
   interpolation cost vs ε·f²p. Decides P4's GAP within the engine's sector
   conventions, and tests the F-R4-transfer claim (value-distribution frozen ⟹
   cost floor) directly.
4. **NR-D1b c_h (cheap):** numerically integrate the one-constant Frank energy for
   radial vs hyperbolic hedgehogs (confirm degeneracy at 8π under E = c_h K r), then
   with two/three elastic constants scan the anisotropy needed to reach 8π/3 —
   recovers or refutes the [IM] band's lower end.
5. **NR-D1b binding fraction (moderate):** on the engine's NR-D1-class configuration,
   measure E₃/E₁ at pitch separation; the eternal verdict needs < 1/3 at the ε-tail
   bottom.
6. **NR-K3a second moment (cheap, decisive for the estimate):** Monte-Carlo a static
   line network (spacing L, tube radius a), convolve with a 2-fm double-propagator
   kernel, and compute ⟨|A|²⟩/|⟨A⟩|² across (a, L, r); then repeat with a
   time-averaged (moving-web) kernel to price the motional-narrowing rescue. This
   converts the DEFECT-CANDIDATE into either a confirmed defect or a discharged one.
7. **λ_dB / σ/m prints (trivial):** add the two NR-D2 minor slips to the tier-0
   gauntlet as exhibits (both directions priced).

— end of memo —
