# Digest: "Tracking the baryon number with nuclear collisions"
## STAR Collaboration, arXiv:2408.15441 [nucl-ex], 27 Aug 2024 — read in extreme detail
### (Short title: "Do quarks or gluons carry baryon number?")

**Phase:** Tier 10 (Phase U).  **Date:** 2026-08-16.  **Source:**
user-uploaded PDF (24 pp: main text pp. 1–9, acknowledgments p. 10,
Materials & Methods pp. 11–19, figures + references pp. 20–24); also
inspirehep.net/literature/2822399.  **Register note (binding for all
Tier-10 use):** [IM]-class imported experimental anchor.  Within-model;
nothing in the campaign's use of it bears on nature.

---

## 1. The question and the two scenarios

Baryon number B is conserved since baryogenesis.  Two candidate carriers:
1. **Valence quarks** (conventional): each carries B = 1/3, Q ≠ 0.
2. **The baryon junction** (Artru 1975 [ref 3]; Rossi–Veneziano 1977
   [4], 2016 update [10]; Kharzeev 1996 "Can gluons trace baryon
   number?" [5]): a **non-perturbative Y-shaped topology of neutral
   gluons** carrying B = 1 and Q = 0; the three ends connect to valence
   quarks with few exceptions.  Studied in lattice QCD [6–8:
   Suganuma/Takahashi three-quark potential; Bissey gluon-flux
   distribution].  Junction-without-quarks exceptions named: baryonium
   tetraquarks [10], **exotic gluonic states of junction and
   anti-junction such as baryonium glueballs, gluonic graphene and
   buckyballs** [9: Csörgő–Gyulassy–Kharzeev, J. Phys. G 30, L17
   (2004)].  None discovered.  Neither carrier scenario experimentally
   verified before this work.

**The discriminating physics:** valence quarks carry a large fraction of
the colliding baryon's momentum → hard to transport from beam rapidity
to mid-rapidity; junctions are made of **low-momentum gluons → easy to
transport**.  Junction picture ⇒ enhanced baryon transport ("baryon
stopping") relative to charge transport.

## 2. The three independent STAR measurements

**(a) Isobar collisions (Ru+Ru vs Zr+Zr, √s_NN = 200 GeV, ~2×10⁹ events
each, 2018 dataset).**  Observable: net-baryon ⟨B⟩ over the isobar
net-charge difference ΔQ (Eq. 3, via double ratios R2_π, R2_K, R2_p —
precision ~10⁻³ since detector conditions cancel).  Valence-quark
expectation: ⟨B⟩/ΔQ ≈ A/ΔZ (=96/4); junction expectation: > A/ΔZ (the
junction carries B but no Q).  **Result: ⟨B⟩/ΔQ × ΔZ/A = 1.84 ± 0.02
(stat) ± 0.09 (syst) ± 0.16 (feed-down) in 0–10% central** — far above
the valence expectation 1; decreases monotonically toward peripheral
(partly a neutron-skin effect: Zr's larger neutron skin, captured in
trend by TRENTO).  Model failures: UrQMD (no junctions) 0.5–0.7; AMPT
even smaller; HERWIG 7.2 / PYTHIA 8.3 default (p+p) 0.5–0.6.  **PYTHIA
8.3 with Color Reconnection Mode 2 — which dynamically forms junctions —
rises to 0.99 ± 0.03.**  Neutron yields estimated (not measured) via
d̄/d ratios + hyperon feed-down (Eqs. 5–12, THERMUS-fitted; ~30% of
neutron yield is feed-down); ⟨B⟩/ΔQ approximates the true ratio to ≤1%.

**(b) Photonuclear γ+Au (from Au+Au at √s_NN = 54.4 GeV; ⟨E_γ⟩ ≈ 0.8
GeV, W_γN ~ 9 GeV; ~2M events selected by 1n/Xn ZDC asymmetry + BBC/VPD
cuts).**  The photon carries no baryon number → clean baryon-transport
probe.  Net-proton rapidity density fit with f(y) ∝ exp(−α_B Δy),
Δy = Y_beam − y (Eq. 4).  **Result: α_B = 1.04 ± 0.22**, consistent with
the Regge-theory junction prediction **0.42 < α_B < 1** (J+ℙ = 0.42
lower limit, J+J = 1 upper limit) [5].  Antiproton slope 0.02 ± 0.05
(flat — pure pair production, as expected).  PYTHIA 8.3 γ*+p (both
default and CR tunes) **overpredicts the slope at ~2.6** — valence
transport too steep.

**(c) Au+Au exponential scaling (√s_NN = 7.7–200 GeV, published
net-proton mid-rapidity yields).**  Fit vs Δy across beam energies:
**α_B = 0.64 ± 0.05, independent of centrality** — inside the junction
Regge window; consistent with γ+Au within 1.7σ.  UrQMD gives the
*opposite* centrality trend (α_B decreasing peripheral→central by ~0.1,
the valence multiple-scattering expectation [35: Itakura–Kovchegov–
McLerran–Teaney]).  Au+Au data consistent with PYTHIA 8.3 CR, not with
default PYTHIA/HERWIG.

## 3. Conclusion as printed

"These findings, corroborated by previous measurements in Au+Au
collisions, **disfavor the valence quark picture**… for [alternative
theories] to be viable, they should simultaneously account for all the
phenomena, an achievement that currently appears to be possible with the
**baryon junction picture**."  Future: EIC u-channel backward production
[36, 37] to probe the junction distribution in nuclei directly.

## 4. Honest caveats the paper carries

- Neutrons are estimated, not measured (feed-down uncertainty is the
  dominant systematic, quoted separately: ±0.16).
- The peripheral isobar trend is contaminated by the neutron-skin
  geometry effect (TRENTO reproduces the *trend* without any junction).
- PYTHIA 8.3 CR is only "a partial imitation of the baryon junction
  mechanism"; both PYTHIA tunes overpredict the γ+Au slope — "possibly
  due to a lack of certain physics ingredients".
- α_B could depend on collision kinematics [3, 21]; current data
  insufficient to resolve.
- γ+Au absolute yields carry an unaccounted selection-efficiency
  normalization (slopes unaffected).

## 5. GUM hooks (campaign-side; why this digest exists)

The corpus's color sector (VII.J) is *architecturally* a junction
theory:
- **Q-1** [DF-structural]: finite-energy states are line-neutral
  composites — "K–K̄ mesons; **n-junction baryons**; junction order n
  owned by WS1-H3, with the tetraquark kill armed."  The baryon IS a
  Y-junction of ℤ₃ mismatch tubes (disclination lines) terminated by
  knots.
- **The Y-law junction geometry is a banked structural consistency**
  ("genuine central junction, favored by lattice QCD at long distance")
  — the same lattice three-quark-potential literature this STAR paper
  cites as refs [6–8].
- Knots carry the **electric charge** (winding) and the momentum/mass
  (m̃-heavy); the junction+tubes are light network topology — the GUM
  cartoon of "junction = low-momentum gluonic structure, quarks = hard
  valence carriers" is structurally the STAR cartoon.
- **The corpus-internal counterpart of the paper's title question is
  sharp and (probably) unadjudicated: in GUM, do knots or junctions
  carry baryon number?**  WS-K's charge census enumerates "web-line
  charges — π₁-classes…; junction constraints; reconnection events
  (net-conserving)" without, to current knowledge, printing a baryon-
  number dictionary entry (B = junction count?).  Tier-10 archaeology
  question #1.
- The junction-without-quarks exotics named by the paper (junction–
  anti-junction "baryonium glueballs", gluonic graphene, buckyballs
  [9]) are, in GUM language, **knot-free tube-network composites with
  junctions** — a state class adjacent to Tier-9's closed loop, equally
  NOT-CONSTRUCTED, and computable in the same validated closed-tube
  machinery (T1) extended by the printed Y-law junction: the J–J̄
  three-tube dumbbell is the minimal case.
- [IM] anchors for Tier-10: ⟨B⟩/ΔQ × ΔZ/A = 1.84 ± 0.02 ± 0.09 ± 0.16
  (0–10%); α_B(γ+Au) = 1.04 ± 0.22; α_B(Au+Au) = 0.64 ± 0.05,
  centrality-independent; Regge junction window 0.42 < α_B < 1;
  PYTHIA-CR 0.99 ± 0.03.

## 6. Verification targets for the Tier-10 literature sweep

Artru NPB 85, 442 (1975); Rossi–Veneziano NPB 123, 507 (1977) + JHEP 06,
041 (2016); Kharzeev PLB 378, 238 (1996) (the α_B intercepts);
Csörgő–Gyulassy–Kharzeev J. Phys. G 30, L17 (2004) (buckyball/graphene
state masses if printed); Takahashi et al. PRL 86, 18 (2001) + Suganuma
(Y-law lattice); Bissey PRD 76, 114512 (2007); Lewis et al. EPJC 84, 590
(2024) (the companion junction-search methodology paper); PYTHIA CR Mode
2 junction formation (Christiansen–Skands JHEP 08, 003 (2015));
post-2024 reactions/criticism of the junction interpretation; EIC
u-channel program (Cebra et al. PRC 106, 015204 (2022)).
