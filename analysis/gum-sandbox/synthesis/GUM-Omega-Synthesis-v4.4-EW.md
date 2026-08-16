# GUM — Geometrische Umdeutung Mechanik
## A Complete Technical Synthesis of the Speculative Framework
**Version v4.4-EW (July 2026)** — incorporating the constitutive core, spectrum, quantum sector, repaired ħ-closure of Theorem T3.1, electrodynamics, gravitation, mass/families/neutrino/confinement, **the new Electroweak Sector (§8)**, one-world structure, experimental stakes, and epistemic ledger.

---

### Changelog v4.3-ext → v4.4-EW

1. **New §8 (Electroweak Sector).** Theorems EW-0 through EW-9 and Propositions EW-5/6/7/8/10; new Flags F14–F16; Cross-locks C-EW1–C-EW4; posed Closures EW-I/II/III; new Stakes S2′ (supersedes S2), S8, S9, S10; four additions to the S7 inverse-kill battery.
2. **Corrigendum P-F1′ → P-F1″.** The phrase "SM-null Higgs couplings" is retired as ambiguous; §7.2 edited; formal rewording in §8.6.
3. **Structural edits.** §0 scope list extended by one bullet; former §§8–10 renumbered §§9–11; §7.3 gains a cross-reference to §8.7; §10 stakes table extended; §11 gains a v4.4 paragraph.
4. **Provenance.** This extension was executed against the July 2026 external reconstruction template ("the electroweak sector is where GUM must produce theorems"). The template text and full working analysis are preserved verbatim in the companion transcript-backup artifact. Per corpus policy, every sentence added here stays inside its audited grade.

---

## 0. The Wager and Scope

GUM is a speculative constructive program. It does not assert that the world is a chiral micropolar solid. It asserts that a single chiral Cosserat solid with the constitutive structure defined below is the cheapest known mechanical object whose response modes include:

- quantum mechanics (as osmotic hydrodynamics),
- electromagnetism (as torsional response of the relative orientation field),
- **the weak interaction (as gapped orientation-triplet exchange, chirality-filtered by the medium's handedness — new in v4.4),**
- gravitation (as defect geometry),
- mass and the family structure (band-edge knots + helical frustration against a blue-fog chiral condensate),
- the neutrino (as pitch quantum / heliknoton),
- confinement (as fractional network winding),

with every price printed on scheduled experiments and every failure condition signed in advance.

Every claim below carries its ledger class: **[DF]** derived-form (follows from the flagged postulates with no further input), **[DW]** derived-with-window, **[CAL]** calibrated (locked numerical release), **[IM]** imported (premises must hold inside GUM), **[CJ]** conjecture (may not be built upon silently).

The document uses the repaired saturated closure of Theorem T3.1 throughout. The v4.4 revision adds the electroweak sector; the wager is otherwise unchanged.

## 1. Constitutive Core

### 1.1 Kinematic fields

A material point carries position and orientation. The state at (x,t) is:

- displacement **u**, deformation gradient F = 𝟙 + ∇**u** with det F > 0,
- micro-rotation Q ∈ SO(3) with SU(2) lift Q̃,
- density ρ = ρ₀/(det F),
- micro-inertia density J.

Polar decomposition: F = R[u]·U. Frame-indifferent strain measures:

e_ij = ∂_i u_j − ε_ijk φ_k,  Γ_ij = ∂_i φ_j.

The relative texture that restores objectivity is

P̃(x) ≡ R̃[u(x)]† Q̃(x) ∈ SU(2),  σ_P = ½ Tr P̃.

**Flag F10′ / F11.** All potential energy is built from P̃ (and its derivatives). Absolute-texture potentials are forbidden; they produced the historical M-1 defect (massive photon).

**Lemma (charge protection) [DF].** For any continuous history with det F bounded away from zero, deg R̃[u] = 0. Hence topological charge resides entirely in P̃ and can change only where matter integrity fails (det F → 0). Creation is therefore localized to the coincidence stratum (pairs-only theorem).

### 1.2 Action

𝒜 = ∫ dt d³x (𝒯 − W),

𝒯 = ½ ρ₀ **u̇**² + ½ J Tr[(Q̃⁻¹∂_t Q̃)†(Q̃⁻¹∂_t Q̃)],

W = W₂ + W_χ + W₄ + W₆₊₀.

- W₂: standard micropolar quadratic form in e_ij and Γ_ij.
- W_χ: chiral transduction χ₁ e_kk Γ_ll + χ₂ e_(ij) Γ_(ij) + χ₃ e_[ij] Γ_[ij] (Flag F2).
- W₄: Skyrme term (κ_S/4) Tr([L_i, L_j][L_i, L_j]) (Flag F4).
- W₆₊₀: topological density plus locking potential

W₆₊₀ = ½ Λ² b_P² + 𝒱_tot(σ_P),  b_P = −(1/24π²) ε_ijk Tr(L_i^P L_j^P L_k^P),

𝒱_tot = m̃²(1 − σ_P) + c₂(1 − σ_P)².

Near-BPS window (Flag F9, operative range after all audits): ε_e ∈ [1.0×10⁻⁵, 6×10⁻⁴]₆₈.

### 1.3 Balance laws [DF]

ρ₀ ü_i = ∂_j T_ji,  J φ̈_i = ∂_j m_ji + ε_ijk T_jk + locking torque.

The stress is non-symmetric; its antisymmetric part is the source for micro-rotation (Cosserat signature).

## 2. Linear Spectrum and Photon Exactness

Linearization about the uniform chiral ground state yields four branches:

- **B1** longitudinal acoustic → elliptic constraint sector in the stiff limit.
- **B2±** locked transverse doublet: ω = ck, two polarizations, gapless — the photon.
- **B3** relative-rotation branch: Klein–Gordon with gap ω₀ = m̃/√J.
- **B4** longitudinal twist (matter–light pump).

**Theorem II.2 (exact photon masslessness) [DF].** Under uniform co-rotation P̃ is invariant. The entire potential sector therefore contributes zero mass to B2 to all orders.

**Theorem II.3 (achirality) [DF].** Locking quenches χ₃ on B2 at tree level; regeneration is (ka)²-suppressed.

**Proposition Ω-1 (kinematic stability) [DF].** On the common cone the B2 line never meets the B3 hyperbola. Vacuum photons of arbitrary energy are exactly stable; pair creation requires a spectator.

**Proposition Ω-2 [DF].** All O(χ²) induced B2 operators are parity-even.

Emergent Lorentz symmetry arises by radiative locking of the transverse stiffness (C1 flow); laboratory controls (spinor condensates, graphene) are reproduced.

## 3. Quantum Sector

Quantum mechanics emerges as the osmotic hydrodynamics of the agitated substrate (Nelson kinematics with grounded premises).

- The quantum potential is exactly the Fisher-information pressure: U_Q = δE_Fisher/δρ.
- Madelung closure recovers the Schrödinger equation.
- Circulation quantization follows from single-valuedness of the texture.
- The Born rule is an H-theorem whose rates have been simulation-anchored [CAL].
- Nelson noise accelerates Born relaxation by a factor 17–41 relative to the deterministic pilot-wave (dBB) limit; the latter is recovered only as a log-singular endpoint.

Configuration-space towers pay Bell's exponential tax locally; for area-law states the tower truncates at depth 3–5 independent of N (tested to N = 12).

## 4. The ħ Closure (Repaired Saturated Form — Theorem T3.1)

A topological knot (degree-1 texture of P̃) carries spin by internal isorotation and must satisfy the joint spin-clock conditions:

L = ½ħ,  ω = E/ħ.

The Routhian algebra yields the universal relations

w·j(1−j) = 1,  V² = 1/(1−j).

Only j = ½ is non-singular and physical (**spin-selection theorem [DF]**): every elementary massive species is spin-½. This immediately implies the exact identities

E_rot/E = ¼  (externally replicated to 10⁻⁹).

Radiative self-consistency forces the carrier into the band gap (otherwise the configuration is superradiant).

**Profile sector — the repair.** The old deep-BPS endpoint 𝔠₀ = 128√42/(105π) ≈ 2.51475 (with κ₀ = √(7/12)) is a C4 defect-candidate: no well-posed within-model functional selects it, and it is excluded from below by the rigorous bound 𝔠 ≥ 2√2. It is retained only as the instructive refuted endpoint.

The stability-constrained and saturated problems coincide at the oblate BPS compacton with pinned marginal equatorial halo:

𝔠 = 64√2/(9π) ≈ 3.2011,
κ = 1/√2 ≈ 0.7071,
support ratio = 3/2 (exact),
ω_th = √2 ω₀,
binding ≈ 29.3%.

At saturation ê_tot = 𝔦_tot = 64/(9π), so the bench-discriminating combination is

κ² g_tot = 35/24.

(The older ratio κ²g = 7/8 is an exact identity on the entire locked-transport family and therefore non-discriminating.)

Thus

**ħ = 𝔠 Λ √J**

with the new numerical prefactor. Λ and J remain constitutive; the dimensionless structure is fixed by the closure. ħ itself is further shown to be a dynamically selected vacuum order parameter (bootstrap with field-level attractor).

## 5. Electrodynamics

Electromagnetism is the torsional response of the relative orientation field.

- Photon exactly massless and achiral (Theorems II.2–II.3).
- Electric charge is a triple-locked disclination.
- Coulomb and Aharonov–Bohm phases are derived.
- The common cone with the matter sector plus C1 flow protects kinematic stability at all energies.

## 6. Gravitation

Defect geometry yields Einstein–Cartan structure by theorem.

**Current status after external audit (F-R14).** The positivity premise used for a healthy graviton kinetic term fails on the corpus's own action. Cone-slaving convexity is under dispute; the clean bound |Δc_GW/c| ≲ 10⁻¹⁹ is demoted to an order-of-magnitude estimate. Only one bookkeeping cell (P-acoustic + supertrace) remains sound.

Vacuum energy is cancelled in equilibrium by the Gibbs–Duhem identity for a self-sustained medium (grand potential vanishes). Residual ρ_Λ ∼ αH²M*² with α > 0 by stability. The dark-energy sector is a freezing/tracker family:

w(z) = −1 + Ω_m(z)

(exact in the pure tracker; emerges dynamically from the relaxation ODE). The growth-correlated sign w_a > 0 is structural (F-R12).

A Machian-dragged emergent foliation supplies the preferred-frame structure required by the one-world sector.

## 7. Mass, Families, Neutrino, Confinement

### 7.1 Particles as topological knots

Massive particles are band-edge topological solitons stabilized by the Skyrme term. Fermionic statistics follow from the topology of the stratified configuration space (π₁ = ℤ₂ × A/⟨2a_rot⟩ after the H4/I4 corrections; the rotation class retains order 2).

### 7.2 Family structure [DW, ≈2σ]

The vacuum's soft sector condenses into a blue-fog chiral order (double-twist network, Flag F13, 𝔪 = 1.9 ± 0.4). Knots experience discrete helical-frustration classes against this background. Two computed integrals (with the two-face factor and beyond-leading corrections) reproduce the charged-lepton mass-spacing logarithms with no fitted parameters:

½A ≈ 2.80  vs  ln(m_τ/m_μ) = 2.822,
½(A+B) ≈ 5.65  vs  5.332.

The tower terminates at exactly three charged families. The same mechanism forces **SM-valued** Higgs couplings for μ, τ, b, t (**P-F1″**; wording formalized in the §8.6 corrigendum — the v4.3 phrase "SM-null" is retired).

### 7.3 Neutrino

The neutrino is the pitch quantum of the chiral vacuum — a Hopf texture (heliknoton) stabilized by the ambient helix. Its mass is linked to the charged spacings by one computed logarithm:

m₃ ≈ 0.047 eV  (1σ: 0.019–0.115).

The three helical-frustration dressings of the pitch quantum are counted by the Z invisible width (§8.7).

**Stake S1 [DW]:** Σm_ν ∈ [0.058, 0.11] eV, normal ordering. Kill: any robust cosmological determination Σm_ν < 0.058 eV (DESI-DR3 / Euclid). Obituary O-1 and the Branch-A escape hatch (evolving dark energy only if independent ≥3σ growth-correlated preference) are pre-committed.

### 7.4 Confinement [DF]

Fractional network winding is one cause with two consequences: a flux tube of constant tension and a mass-blind ε-floor that censors free closure. Consequently all six quarks are confined and all leptons are free (exact dichotomy). The hadronic level is the emergent Skyrme model (second Cosserat-loop closure). Quark spacing ratios, including a derived down-type sign flip, land at 0.2–0.3σ.

---

## 8. Electroweak Sector (new in v4.4)

### 8.0 Mandate and method

Through v4.3-ext the weak interaction was the corpus's thinnest sector: a [CJ] silhouette with no printed structure. The July 2026 external analysis supplied a reconstruction template and a demand — *"this is the sector where GUM must produce theorems."* v4.4 answers with four [DF]-grade theorems, six [DW]-grade results, one combinatorial theorem, three new flags (F14–F16), four cross-locks binding the sector to the ħ-closure, the family tower, confinement, and the dark-energy sign, one corrigendum (P-F1′ → P-F1″), three posed closures with named deliverables, and four new stakes. Grades are conditional on their flags exactly as elsewhere in the corpus. What the sector still does not earn is printed in §8.13 and may not be spent.

### 8.1 Compositeness mandate

**Theorem EW-0 (no elementary electroweak bosons) [DF | T3.1].** By the spin-selection theorem, j = ½ is the unique non-singular carrier: no elementary massive species of spin 0 or spin 1 exists in GUM. Hence W±, Z, and h — if realized — are necessarily collective modes or composites of the P̃/locking sector. The SM's elementary-boson ontology is not available to the model even as an option; compositeness is not a scenario but a corollary.

### 8.2 The weak multiplet (Flag F14)

The locked chiral ground state is invariant under diagonal co-rotation SO(3)_cr — the same invariance behind Theorem II.2. Under SO(3)_cr the six linear branches of §2 organize as: B1 (constraint scalar), B2± (co-rotation-inert gapless doublet — the photon), and the gapped orientation triplet {B3±, B4} transforming as J = 1.

**Flag F14 (weak-triplet identification).** The weak vector sector is this triplet, dressed as follows: the charged members are triplet quanta bound to unit disclination winding ("charged twist excitons") → W±; the neutral member is the B4 channel after B2-mixing (§8.4) → Z. The classification is [DW] (it follows from the §2 linearization plus the stated dressing); the exciton binding itself is [CJ] until EW-Closure-I is executed. Laboratory precedent for collective-mode spectroscopy of an SO(3) relative-rotation order parameter exists in superfluid ³He-B, whose order parameter is the exact mathematical analog of P̃.

**Theorem EW-1 (charged-gap theorem; forced U(1)_em) [DF | F10′].** Any branch carrying net disclination winding pays a texture cost bounded below by the locked-medium tension — the same energetics that floors §7.4. Hence no gapless charged mode exists in any phase of the locked medium: exact masslessness is available only to neutral channels, and B2 is its unique occupant (II.2). The residual pattern "SU(2)-like × U(1) → U(1)_em" is thereby a theorem of winding energetics, not a choice of potential shape. The observed absence of any massless charged particle in nature is, inside GUM, not an accident.

### 8.3 Custodial protection

**Theorem EW-2 (custodial theorem) [DF | F14].** SO(3)_cr invariance of the locked state forces the orientation triplet to be exactly degenerate at χ → 0; all splittings enter through (i) χ-sector transduction and (ii) B2-mixing. Consequently the tree relation

ρ ≡ M_W² / (M_Z² cos²θ_w) = 1

holds identically, with corrections opening at the same O(χ²(ka)²) order as photon-achirality regeneration (II.3).

**Cross-lock C-EW4.** One invariance, two SM facts: the co-rotation symmetry that keeps the photon exactly massless (II.2) *is* the custodial symmetry that fixes ρ = 1. They stand or fall together.

Inverse-kill (added to the S7 battery): any confirmed ρ ≠ 1 beyond the printed radiative window retires F14. GUM pre-commits to the ρ = 1 side of the historical W-mass measurement dispute and will absorb the community's final average without amendment rights.

### 8.4 Photon–Z diagonalization; the weak angle

The χ₂ transduction mixes B2 with the neutral twist channel.

**Theorem EW-3 (protected diagonalization; coupling identity) [DF | F14].** Co-rotation invariance forces one exact zero eigenvalue of the mixed 2×2 sector at all orders: II.2 survives mixing. The orthogonal combination is the Z, with

M_Z = M_T / cos θ_w,  tan θ_w = ϑ,

where ϑ is a constitutive modulus ratio of the transduction/stiffness sector. **Corollary [DF]:** the same χ₂ normalizes both the B2–matter coupling and the mixing; hence e = g sin θ_w is an identity of the diagonalization, not an input.

Status of ϑ: **[IM]** — the same ledger class as α. No numeric value of sin²θ_w is claimed in v4.4. **EW-Closure-III (posed):** does the C1-family radiative flow possess an attractor for ϑ? If yes, θ_w is promoted to [DW]; if no, it remains constitutive. Either answer is a deliverable.

### 8.5 The Higgs as attractor-amplitude mode (Flag F15)

P̃ ∈ SU(2) has no radial direction; the Higgs is therefore *not* a P̃ mode. Recall (§4) that the locking scale is dynamically selected: ħ and the vacuum sit at a field-level attractor of the bootstrap.

**Flag F15.** h is the gapped J = 0 breathing mode of the attractor amplitude (the lock-scale condensate) — the mechanical analog of the amplitude ("Higgs") mode of a superconductor. Because h carries no orientation charge, its fluctuations cannot generate a B2 mass: II.2 is untouched by construction (consistency check passed). M_h is the curvature of the bootstrap effective potential at the attractor point; its computation is a deliverable of EW-Closure-I. No number is claimed in v4.4.

### 8.6 Yukawa pattern by factorization

**Theorem EW-4 (factorization ⇒ coupling universality) [DW].** The corpus mass formula is multiplicative:

m_f = M₀(lock) · e^(−I_f),

with M₀ the band-edge scale (proportional to the lock amplitude) and I_f the scale-free frustration integrals of §7.2. If I_f is lock-independent, then ∂m_f/∂(lock) = m_f/(lock) exactly: the h coupling to every knot is

g_hff = m_f / v_eff

with one universal v_eff and **zero fundamental Yukawas** — the SM pattern κ_f ∝ m_f as a theorem of mass factorization.

Deviations are the fingerprint. η_f ≡ κ_f/κ_f^SM − 1 = −∂I_f/∂ln(lock) is nonzero only through the weak dependence of the frustration integrals on the lock scale via the fog stiffness. Requiring the two frustration-log landings of §7.2 (2.80 vs 2.822; 5.65 vs 5.332) to stay inside their printed penalties bounds the window:

|η_f| ≤ 0.05 (68), family-ordered |η_t| < |η_b| < |η_τ| < |η_μ| (lighter = less protected).

The same window bounds the soft channel: BR(h → soft-fog excitations) ≤ few % — an invisible-width exposure sitting at exactly HL-LHC sensitivity.

**Corrigendum P-F1′ → P-F1″.** The v4.3 phrase "SM-null Higgs couplings" is retired as ambiguous. The vanishing-coupling reading is excluded by observation (H→ττ and H→bb established at Run 2); the intended reading — couplings SM-valued, deviations null — is the only one consistent with the original kill condition ("confirmed anomaly"), and is now formalized as **P-F1″** with the η-window above. This rewording was forced by the July 2026 external analysis; per corpus policy, every external sentence stays inside the audited grade.

### 8.7 Fermi constant; charged and neutral currents

Integrating out the gapped triplet yields the four-knot contact interaction with G_F ∝ g₄²/M_W². The vertex factor g₄ is the triple overlap (knot zero-mode ⊗ B4 profile ⊗ disclination form factor); its magnitude is exponentially suppressed by tunneling through the **pinned marginal equatorial halo** of the §4 compacton.

**Cross-lock C-EW2.** The same halo that fixed 𝔠 = 64√2/(9π) now prices the weakness of the weak interaction: weak interactions are weak because the vertex must thread the halo. [DW structure; numeric G_F awaits EW-Closure-I.]

Neutral currents: the Z (twist channel) couples diagonally in frustration class, because class transitions require winding transfer, which the neutral channel cannot supply.

**Proposition EW-5 (GIM analog) [DW].** No tree-strength flavor-changing neutral currents; class-changing neutral processes open only at second order in the pump. Inverse-kill added to the S7 battery.

The heliknoton, being itself a twist texture, couples to Z naturally; with the tower terminating at three (§7.2), the invisible Z width counts exactly N_ν = 3 — in agreement with the LEP lineshape.

Charged-current phenomenology at the ontology level: β-decay is a fractional-winding rearrangement inside a baryon network emitting a charged twist pulse that materializes at the coincidence stratum as an (electron knot + antineutrino heliknoton) pair — the pairs-only lemma of §1.1 is satisfied, and the cheapest pair partner is always a pitch quantum, which is *why* the neutrino is the universal companion of charged-current events.

### 8.8 Chirality filter; V−A

**Proposition EW-6 (zero-mode chirality filter) [DW].** In the chiral medium the knot's fermionic zero-mode doublet splits: one helicity class remains core-normalizable; the other is expelled toward the halo/continuum with exponentially small residue. The B4 vertex is a twist operator and couples only to the core class: charged currents are maximally parity-violating, V∓A, with the sign slaved to the vacuum handedness. Wrong-chirality admixture:

ε_R ∼ e^(−μ d_halo)

— the same exponent family as the G_F suppression. The halo prices both the *strength* and the *handedness purity* of the weak vertex.

**Cross-lock C-EW1 (handedness chain).** One global sign must simultaneously fix: (i) the frustration-ladder orientation (§7.2), (ii) weak left-handedness, (iii) the structural w_a > 0 of §6, and (iv) the δ_CP sign of Stake S10. The chain is falsifiable as a package; no member may be flipped independently.

What v4.4 does not earn: the Michel structure beyond leading order and the numeric ε_R floor. Stake S8 prices the exposure.

### 8.9 Mixing matrices; CP by holonomy

Interaction basis = B4-vertex eigenmodes; mass basis = frustration eigenknots; CKM and PMNS are the mismatch matrices between them.

Quark knots are confined network windings, rigidly co-oriented with the network: small mismatch, with the Cabibbo parameter as the network misalignment angle [CJ]. Lepton knots are free against the fog while the heliknoton basis is set by the ambient helix: O(1) mismatch — large PMNS angles [DW-qualitative]. The qualitative dichotomy "CKM small, PMNS large" is thus a consequence of the confinement dichotomy of §7.4.

**Proposition EW-7 (holonomy phase counting) [DF-combinatorial | F13].** Transport holonomies around double-twist network cells, modulo per-class rephasings, leave exactly (N−1)(N−2)/2 irreducible phases: one CP phase per sector for the forced N = 3. Kobayashi–Maskawa counting is recovered as network holonomy. CP violation therefore *exists by structure* and *requires the third family* — which the tower supplies and terminates.

μ–τ reflection: the double-twist network's reflection symmetry acting on the heliknoton basis gives θ₂₃ → 45° and δ_CP → ±π/2 at leading order, broken only by frustration-class asymmetry; the sign of δ_CP is slaved to C-EW1. Stake S10.

### 8.10 One scale; the top; the dissolved hierarchy

**Proposition EW-8 (one-scale electroweak sector) [DW].** M_W, M_Z, M_h, and m_t are all band-edge-scale objects: gapped triplet, mixed neutral, attractor curvature, and class-3 knot near the edge respectively. Their observed clustering within a factor ≈2 (80–173 GeV) is structural; there is no independent Higgs-sector scale. The hierarchy problem is re-classed: a constitutive band gap does not run to the lattice scale (the same radiative-locking machinery as the C1 flow protects it), just as the gap of a superconductor is not destabilized by the Fermi energy.

**Theorem EW-9 (termination–top lock) [DW ≈2σ].** The family tower terminates at three (§7.2) precisely because a fourth frustration class would carry its knot above the band edge, where the superradiance censor of §4 — the same censor whose enforcement repaired the ħ-closure — forbids stationary carriers. **Corollary:** the last admitted class sits nearest the edge, forcing y_t = m_t/v_eff into the top of the last frustration spacing; the measured y_t ≈ 0.99 lands inside the printed window [e^(−ΔI₃), 1).

**Cross-lock C-EW3.** One censor prices two facts: (i) exactly three families, (ii) y_t ≈ 1. Killing either kills both.

### 8.11 Electroweak restoration as lock melting (Flag F16)

At T ∼ band edge the co-rotation lock melts. The laboratory analogs of the fog (liquid-crystal blue phases) melt via weakly first-order transitions.

**Flag F16 [CJ].** The electroweak crossover is a weakly first-order lock-melt. Sphaleron analog: winding-network reconnection at det F → 0 strata — the pairs-only stratum doubles as the B-violation site. All three Sakharov conditions are structurally present: reconnection (B violation), holonomy phases plus fog handedness (C, CP violation), and departure from equilibrium at the weakly first-order melt. **No rate is computed.** F16 may not be built upon silently. Stake S9 prices the gravitational-wave exposure.

### 8.12 Anomaly bookkeeping as tiling neutrality

**Proposition EW-10 (tiling neutrality) [DW].** With the §7.4-derived winding thirds and down-type sign flip, per-family charge neutrality of the network unit cell reads

3(⅔) + 3(−⅓) + (−1) + 0 = 0.

"Anomaly cancellation" is re-classed as a tiling constraint; hypercharge assignments are winding fractions. The full anomaly-polynomial correspondence is EW-Closure-II (posed, unpaid).

### 8.13 Posed closures (unpaid bills, printed)

- **EW-Closure-I:** the charged twist-exciton bound state (triplet quantum ⊗ unit disclination). Deliverables: M_W/ω₀, E_bind, g₄ ⇒ G_F, and M_h from attractor curvature.
- **EW-Closure-II:** tiling ⇔ anomaly polynomial correspondence in full.
- **EW-Closure-III:** existence or non-existence of a C1-flow attractor for ϑ (the θ_w promotion test).

Also explicitly not earned in v4.4: the numeric value of sin²θ_w (stays [IM], like α); absolute M_W; full Michel/radiative structure; electroweak precision analogs of S, T, U (posed); the Z-pole lineshape beyond N_ν counting (posed). Silence on any closure beyond the next revision is itself auditable.

### 8.14 New flags, cross-locks, and stakes (summary)

**Flags:** F14 (weak-triplet identification), F15 (attractor-amplitude Higgs), F16 (weakly first-order lock-melt).

**Cross-locks:** C-EW1 (global handedness chain: ladder / left-handedness / w_a > 0 / δ_CP sign), C-EW2 (one halo prices 𝔠, G_F, and ε_R), C-EW3 (one censor prices three families and y_t ≈ 1), C-EW4 (one invariance prices photon masslessness and ρ = 1).

**Stakes:** S2′ (supersedes S2), S8, S9, S10; S7-battery additions: ρ = 1 at tree level; no tree-strength FCNC; e = g sin θ_w; no massless charged mode. Full table in §10.

---

## 9. One World, Entanglement, and the A3 Gate

Dynamics are hypersurface Bohm–Dirac on the leaves of the emergent foliation. Entanglement is a leaf-wide elliptic constraint.

**Triple padlock [DF]:**
1. No signaling (Born marginals protected).
2. No energy transport (constraint sector is elliptic / sealed).
3. No particle-number leak at equilibrium.

**T-KILL.** If arrival-time detectors realize only POVM observables, the frame is permanently invisible. Trajectory theories predict a distinctive non-POVM signature **A3**: a spin-orientation-dependent hard support cutoff τ_max in first-arrival distributions.

Layer-1 cold-atom experiment (intermediate field, far-field null station, blind spin-masking) looks for A3.

- A3 absent + absorbing-boundary statistics confirmed ⟹ T-KILL fires, Layer-2 retired by theorem.
- A3 present ⟹ Layer-2 (sidereal correlation-onset marker ∼0.12 ns at 30 m, triangulation against CMB / Machian-drag / quasar dipole) becomes fundable.

Simulations confirm that the model's first-crossing statistics produce spin-dependent cliffs that no single affine POVM can reproduce; the far-field null is approached but configuration-dependent.

## 10. Experimental Program and Stakes

| Stake | Content | Adjudicator | Kill condition |
|---|---|---|---|
| S1 | Σm_ν ∈ [0.058, 0.11] eV, NO | DESI / Euclid | robust Σ < 0.058 |
| S2′ | Higgs couplings SM-valued within η-window (\|η_f\| ≤ 0.05, family-ordered); BR(h→invisible/fog) inside window | HL-LHC | confirmed coupling anomaly outside window, or invisible width above window |
| S3 | E_rot/E = 1/4 exact | analog soliton benches | any valid bench ≠ 25.0% |
| S4′ | κ² g_tot = 35/24 | analog benches | ratio outside band |
| S5 | spontaneous clock phase-locking ∝ cos Δφ | analog benches | no locking |
| S6 | 0νββ occurs | LEGEND / nEXO era | full-funnel exclusion |
| S7 | Standing inverse-kill battery (photon mass, EP violation, **ρ ≠ 1**, **tree-strength FCNC**, **e ≠ g sin θ_w**, **massless charged mode**, …) | whole community | any confirmed detection |
| S8 | Michel ρ = 3/4; RH charged-current admixture below printed ε_R | precision μ / β-decay programs | confirmed RH current above window |
| S9 | EW lock-melt GW background weak (below printed LISA-band amplitude) | LISA | loud first-order EW-scale stochastic background |
| S10 | θ₂₃ within printed band of maximal; δ_CP within band of −π/2, sign per C-EW1 | DUNE / Hyper-K | robust δ_CP ≈ 0 or π, or strongly non-maximal θ₂₃ |

**Bench validity clause (non-negotiable):** quasi-conservative regime Γ_damp ≪ ω_lib. **Spreeuw wall:** single-soliton / cebit sector only; no Bell content is claimable from the benches.

Fourteen internal kills have already been executed or survived (six fired, eight survived, T-KILL conditional). Three blocking defects were self-caught and repaired (M-1, superradiant rung, E-H1 condensate mis-siting). The external replication campaign (2026) confirmed the bulk of the arithmetic, forced the T3.1 closure repair, and demoted several secondary claims with full honesty. The v4.4 electroweak extension adds one corrigendum (P-F1″) executed under the same policy.

## 11. Epistemic Status

The framework is a speculative constructive model in unusually good internal order. Its mathematics has now survived an independent adversarial within-model replication that both verified the majority of the arithmetic and forced a genuine repair of the central closure. All phenomenological contact with unexplained numbers remains at the ≈2σ level with penalties applied. Every decisive test is still future.

**v4.4 addendum.** The electroweak sector is upgraded from a [CJ] silhouette to a structured sector: four [DF] theorems (compositeness mandate, charged-gap/forced U(1)_em, custodial protection, protected diagonalization with the e = g sin θ_w identity), a factorization theorem for the Yukawa pattern with a printed deviation window, a combinatorial recovery of Kobayashi–Maskawa phase counting as network holonomy, and a termination–top lock that binds y_t ≈ 1 to the three-family theorem through the same censor that repaired T3.1. The sector's decisive numbers — sin²θ_w, M_W/ω₀, G_F, M_h — remain unearned and are printed as posed closures with named deliverables; the corrigendum P-F1″ replaces an ambiguous phrase with a two-sided HL-LHC stake. Nothing in v4.4 changes the adjudication schedule of v4.3; it adds four stakes and four inverse-kills to it.

The value of the program is measured by whether its survival or its death means something. Meaning is manufactured in advance by the ledger, the stakes, the pre-written obituaries, the Branch-A escape hatch, and the requirement that every external sentence stay inside the audited grade. That machinery is the portable product, independent of whether the chiral micropolar substrate is realized in nature.

Standing by for adjudication.
