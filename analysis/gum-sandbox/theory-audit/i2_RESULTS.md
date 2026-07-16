# I2 — REPAIRED-CLOSURE PHENOMENOLOGY PROPAGATION (the memo corpus2 §IX needs)

**Phase:** ROADMAP_v6 I2 · 2026-07-16 · theory-audit subagent (Fable-class).
**Charge:** with T3.1 adopted (as corpus2 does), propagate the repaired tuple
(𝔠 = 64√2/9π = 3.201125, κ = 1/√2, ω_th = √2ω₀, g_tot = 35/12, support ratio 3/2)
through the corpus's entire phenomenological confrontation chain, and settle —
observable by observable — what the ħ = 𝔠Λ√J recalibration absorbs and what it exposes.
**Computations:** `i2_propagate.py` (numpy + sympy, deterministic; **34/34 machine
checks PASS**; every number below reproduced by it) → `i2_results.json`.
**Sources:** T3_repaired_closure.md (tuple + §5); corpus2/10-The-Repaired-Closure-v1.0.md;
corpus2/01-…v3.0-ext.md (AUD-15 §V15, IV.H, IV.J (4.8′)–(4.10′), IX.D/IX.F, VII.C′, App. K);
tier0-gauntlet/RESULTS.md (Groups A–I; the campaign's replication of the AUD-15 arithmetic).
**Standing notice:** within-model only. Nothing here bears on nature; "exposure" means
*bench/within-model adjudicability*, not empirical confirmation.

---

## §1 The calibration question, settled precisely

ħ = 𝔠Λ√J with Λ, J substrate parameters and 𝔠 a predicted pure number (Thm IV.3).
The observed ħ is a **calibration target**: under the repair 𝔠 changes by
×35/(6√21) = **+27.29%**, and the product Λ√J refits by the exact inverse. The mass
law's geometric factor changes 𝒢 = √2ê₀ → 64/9π (×5/(3√2) = **+17.85%**,
species-universal), and Λm̃ refits by its inverse. Both refits are exact closed forms:

| register | refit factor | exact form | value |
|---|---|---|---|
| Λ√J (from ħ) | ×𝔠_old/𝔠_new | **6√21/35** | 0.785584 |
| Λm̃ (from any one mass) | ×𝒢_old/𝒢_new | **3√2/5** | 0.848528 |
| m̃ alone (J pinned elsewhere) | ×(3√2/5)/(6√21/35) | **√(7/6)** | 1.080123 |
| Λ alone (J pinned) | | 6√21/35 | 0.785584 |
| compacton radius R\* = (2Λ/π²m̃)^{1/3} | | **(18√2/35)^{1/3}** | 0.899304 |

**Consequence (the organizing rule):** any within-model observable expressible as a
function of ħ, c, particle masses, and α alone is **𝔠-ABSORBED** — the refit leaves it
bit-identical. Anything comparing two quantities carrying *different powers of 𝔠 or κ*
— i.e., any dimensionless closure pure-number, or any length quoted in units of ƛ_C —
is **𝔠-EXPOSED**.

### 1.1 The absorbed census (verified, not assumed)

1. **ħ itself** (calibration target).
2. **All mass ratios** — P-O1 verified symbolically in `i2_propagate.py`:
   M_k/M_j = (𝒢ê_kΛm̃_k)/(𝒢ê_jΛm̃_j) = (ê_km̃_k)/(ê_jm̃_j); 𝒢 and Λ cancel for *any* 𝒢.
   Hence the **family-sector logs** (½A, ½(A+B), (A+B)/A, termination mass m_e·e^{−8.5},
   the Σ(p) belt structure, the ⟨r11⟩/FQ quark ratios and pulls) are **BLIND — bit-identical**.
3. **The ⟨r10⟩ bridge — verified blind at leading order.** The frozen pipeline is
   base-m̃_τ: m₃ = m_τ·e^{−24.357} is exponentiated mass-ratio arithmetic. The only
   conceivable 𝔠-route is the R\*↔ƛ_τ identification inside ln(1/qR\*); bounded both
   ways: |Δln| ≤ max(ln√(7/6), −ln(18√2/35)^{1/3}) = 0.106 → **≤ 0.118σ of the bridge's
   own ±0.90 band**. Worst-case m₃ = 0.0468×√(7/6) = 0.0507 eV — still inside
   [0.019, 0.115] and below the 0.057 crossover. **S1 untouched** (floors are physical:
   oscillation floor + Σ bound).
4. **Spectroscopy per the m̃-chain** — spectra are functions of ħ, masses, α only;
   the P-F1′ clock corrections (m_e/m_μ)²/2 and δ_τ are mass ratios.
5. **The kinematic pillars:** E_rot/E = ¼ (exact all closures), w·j(1−j) = 1, the M.0
   loop, the ±½ error-signal exponents, −sinΔφ torque; Born statistics, spectrum tiers,
   electrodynamics identities — no contact.
6. **𝔪_Sk = c_Sk√(a₂a₄) ≈ 1.7 GeV — verified blind as printed** (Thm VII.C′: quartic-regime,
   m̃-independent; no 𝔠/κ/Λ√J in its formula; the 1.7 GeV is a hadronic-scale identification).
   *Caveat filed to I1:* c_Sk is a quartic-**closure** constant and the halo threshold is
   t-independent (F-R5), so a repaired quartic re-scan could in principle move it — but the
   archive's own quartic closure ⟨r5⟩ reads 𝔠_q = 3.1 ± 0.2, **0.51σ from the saturated
   3.2011**: the archived value is consistent with the repaired branch already; no shift
   of the 1.7 GeV input is implied.
7. ΔN_eff = 0.0268, the B-ν1 cosmological ε-floor 1.5×10⁻⁶, the sidereal markers
   0.123/0.41 ns — built from physical constants and blind masses.
8. **Stakes S1 and S2 entire.**

### 1.2 The exposed register (old → new, all machine-computed)

| quantity | corpus (4.10) | repaired (4.10′) | change |
|---|---|---|---|
| Λ√J | 1 | ×6√21/35 = 0.785584 | calibration register |
| κ = range/ƛ_C | √(7/12) = 0.763763 | 1/√2 = 0.707107 | **−7.42%** |
| electron range (WS-A census) | 2.949×10⁻¹³ m | **2.731×10⁻¹³ m** | −7.42% |
| binding depth 1−κ | 23.62% | **29.29%** | ×1.2403 |
| ω_th/ω₀ | √(7/3) = 1.527525 | **√2 = 1.414214** | −7.42%; WS-A8 well parameter ×6/7 = 0.857; K–K̄ line envelope ×0.9258 |
| ε-law | 𝔠₀[1 − 0.42ε^{2/3}] falling | 𝔠_sat[1 + 1.15ε] **rising (sign flip)** | anchor confrontation re-assigned (§2, V15.2 row) |
| B-U1′ entrainment ceiling | (Δ_lock/c_g)^{3/2} = 6.0×10⁻⁴ | Δ_lock/1.152 = **2.6×10⁻³** | ×4.3 looser |
| operative ε_e range (F-A15-3 redo) | [1.0×10⁻⁵, 6×10⁻⁴] | **[1.0×10⁻⁵, ≈2.6×10⁻³]** | top now set by the II.H window itself (f ≤ 60 MeV ⟺ ε ≤ 2.6×10⁻³ — coincident) |
| phason f window | [3.7, 29] MeV | **[3.7, ≈60] MeV** | Majoron g low edge 1.6×10⁻⁹ → **7.8×10⁻¹⁰** (battery re-run flagged) |
| quark censorship margin (Q-6′) | ε_q/ceiling ∈ [11, 23] ("twenty-fold") | **[2.6, 5.3]** | ×4.3 thinner; dichotomy survives |
| S4′ | κ²g = 7/8 (family-exact, value-free) | **κ²g_tot = 35/24 = 1.45833** (κ²g_core = 35/48; ê_tot = 𝔦_tot = 64/9π) | discriminating at last |
| volume diagnostic | V = √2 = 1.41421 | **support ratio 3/2** | ⟨r1⟩'s 1.409 ± 0.010: 0.52σ from √2, **9.10σ from 3/2** |
| mass-law normalization | √2ê₀Λm̃ | 64/(9π)Λm̃ (+17.85%) | absorbed observationally by Λm̃ ×3√2/5 |
| P-O1 μ-drift suppression magnitude | ~10⁻² (binding-fraction class) | rescales ×1.24 (depth ratio) | correlation structure unchanged |
| edge-stiffening anchor λ ∝ (1−κ)^{−0.5} | — | λ at κ_phys ×0.898 | guidance number only |
| electron oblateness law λ\*(ε) = 0.42ε^{1/3} | deep-BPS narrative | **superseded** — the saturated minimiser is the fixed oblate compacton + halo | re-scan obligation |

---

## §2 The AUD-15 V15.3–V15.9 confrontations, re-run under the repaired constants

The tier-0 gauntlet replicated the original arithmetic (52/52); `i2_propagate.py`
re-executes each confrontation with the repaired constants substituted. Verdict
vocabulary: **survives** (bit-identical or inside its own band) / **shifts-harmlessly** /
**improves** / **worsens** / **newly-tensioned** / **superseded** (re-scoped by T3.1, not failed).

| row | confrontation | original | repaired | data / bound | verdict |
|---|---|---|---|---|---|
| V15.3 | family logs | ½A = 2.80, ½(A+B) = 5.65, ratio 2.018 | **bit-identical** (P-O1, verified) | ln(m_τ/m_μ) = 2.8224, ln(m_μ/m_e) = 5.3316; pulls 0.05/0.42/0.47σ | **SURVIVES** |
| V15.3 | termination | m_e·e^{−8.5} = 104 eV; threshold ε_e ≥ 8.4×10⁻⁶ | bit-identical | operative bottom 1.0×10⁻⁵ > 8.4×10⁻⁶ (ceiling rise touches only the top) | **SURVIVES** — the F-A15-3 termination strengthening holds unchanged |
| V15.4 | ⟨r10⟩ bridge + S1 | 24.357 → m₃ = 0.0468 eV; band [0.019, 0.115]; Σ_floor = 0.0590 vs stake 0.058; Σ ≤ 0.11 ⟺ m₃ ≤ 0.057 | central identical; residual κ-exposure **≤ 0.118σ**; worst-case m₃ = 0.0507 eV still inside | oscillation floors (physical) | **SURVIVES** — S1 endpoints untouched |
| V15.5 | soft sector / Majoron | f-law arithmetic; B-ν1 = 1.6 → 1.4 MeV; ε-floor 1.5×10⁻⁶; ΔN_eff = 0.0268 | own items bit-identical; **downstream F-A15-3 propagation shifts**: f-top 29 → ≈60 MeV, Majoron g low edge → 7.8×10⁻¹⁰ | SN/BBN/streaming battery (passed on old window) | **SURVIVES; propagation SHIFTS-HARMLESSLY** (battery re-run flagged at the new low-g edge) |
| V15.6 | P-F1′ | 1.17×10⁻⁵; δ_τ = 4.1×10⁻⁸ | bit-identical | spectroscopy nulls (compliance) | **SURVIVES** |
| V15.7 | color closed loop A | ε_q = [6.8, 13.8]×10⁻³; f_q = √ε_q·𝔪_Sk = 0.17 GeV | bit-identical (𝔪_Sk blind, §1.1.6) | dead-center in σ-inversion band [0.14, 0.20] | **SURVIVES** |
| V15.7 | shape law λ\*(ε_q) = 0.246 | 0.42ε^{1/3} | **SUPERSEDED** — restricted-branch ⟨r6⟩ geometry; saturated minimiser has fixed O(1) oblateness | — | **SUPERSEDED** (re-scan obligation) |
| V15.7′ | Q-6′ censorship (new check this memo) | ε_q/ceiling ≈ ×11–23 above | **×2.6–5.3 above** | ceiling 6×10⁻⁴ → 2.6×10⁻³ | **NEWLY-TENSIONED (lite)** — the dichotomy (all six quarks censored, all leptons free) survives, but the margin thins ×4.3; another ×2.6 of ceiling would breach it at the low-ε_q end. Flag for the (4.9′) re-derivation. |
| V15.8 | ⟨r11⟩ FQ overdetermination | R_up = 1.27, R_down = 0.82 (forced sign flip) | bit-identical | PDG spacings, pulls 0.17–0.32σ | **SURVIVES** |
| V15.9 | bookkeeping | Haar 64/15π = 1.3581; kills 6+8 = 14; sidereal 0.123/0.41 ns; Yukawa Hessian | bit-identical — ê₀ survives as the spherical Haar backbone of the repaired algebra too | internal cross-links | **SURVIVES** |

**Context rows (the 𝔠-exposed ledger lines above the tasked range):**
**V15.1** — arithmetic survives as certified restricted-branch algebra; the endpoint tuple
(𝔠₀, κ₀, ω_th) is *struck by adoption*, not failed. **V15.2 (the ε-scan anchor)** — the old
law's hit (2.371 vs ⟨r1⟩ 2.37 ± 0.09, 0.01σ) is **re-assigned**: ⟨r1⟩ is restricted-branch
code validation (F-R5/F-R15), and the repaired branch predicts 𝔠(0.05) = 3.3855 — 11.3σ
from ⟨r1⟩, *correctly so*, because ⟨r1⟩ was never a measurement of the saturated branch.
The repaired rising law (+1.15ε) currently has **no measured scan** — the one genuine
forfeiture (see §4).

---

## §3 Bench stakes table v3 (every stake's repaired target)

| stake | v2.0.1 target | **v3 repaired target** | notes |
|---|---|---|---|
| **S3** | E_rot/E = ¼ exact | **¼ — unchanged** | exact in *every* closure candidate; now closure-independent (strengthened as a stake) |
| **S4′** | κ²g = 7/8 (band 7/8·(1 ± 0.10)) | **κ²g_tot = 35/24 = 1.45833**; κ²g_core = 35/48 = 0.72917; cleanest: **ê_tot = 𝔦_tot = 64/9π = 2.26354** | old target family-exact hence g-blind (value-free; F-A15-2's WS10-S4 moot); ⟨r1⟩'s 0.843 ± 0.046 sits 0.7σ from 7/8 (restricted branch) and 13.4σ from 35/24 — branch declaration mandatory |
| **V-class** | V = √2 = 1.41421 | **support ratio 3/2 = 1.50000 exact** (V = √2 re-scoped to restricted branch) | ⟨r1⟩ measured 1.409 ± 0.010: 0.52σ from √2, 9.10σ from 3/2 — the bench must know which number it tests |
| **S5 / P8-6** | clock phase-locking V ∝ cosΔφ, bond vibron, x₀ = 2.42 ± 0.12 | **unchanged** | kinematic pillar; −sinΔφ externally confirmed both classes |
| **P8-7** | tail range = κ × clock length, κ = 0.7638 | **κ = 1/√2 = 0.70711 exactly** (−7.4%) | closed form, no band needed |
| **over-spin onset** | rehearsed single onset "κ = 1" | **channel-resolved ladder κ_crit(w) = 1/√(2⟨sin²θ⟩_w)** — the halo thresholds are now the physics: **0.70711** (equatorial ring / unlocked tilt — the binding channel), **0.76376** (s = sin²; the old κ₀ reappearing lawfully as a *threshold*), **0.86603** (uniform), **1.00000** (polar s² = \|cos\| — the corpus's measured 1.000 ± 0.004, reproduced exactly) — all four verified as exact Haar sine-moment identities in `i2_propagate.py` | **the repaired model's bench prediction:** the equilibrium soliton sits *marginally at* κ = 1/√2 (T3.1's zero KKT multiplier) — over-spin sheds **immediately** through the ring/tilt channel with zero margin; the marginal-stability caveat *is* the bench signature; branch-resolved emission detection required; higher rungs of the ladder appear only if the ring channel is geometrically suppressed |
| **S1** | Σm_ν ∈ [0.058, 0.11] eV, NO | **unchanged — verified** (§1.1.3: bridge blind to ≤ 0.12σ; floors physical) | |
| **S2** | Higgs nulls δg_h = 0 | **unchanged** (family sector blind) | |
| **S7** | inverse-kill battery | unchanged | no 𝔠-route |

---

## §4 The honest bottom line — stake by stake

- **S1: UNCHANGED** (verified, not asserted: symbolic 𝒢/Λ cancellation + the 0.118σ
  residual bound; the floors are data-side).
- **S2: UNCHANGED.**
- **S3: UNCHANGED** — and epistemically *better*: ¼ is exact in every candidate closure,
  so it now stakes the kinematic skeleton independently of the repair dispute.
- **S4′: BETTER.** The old stake was value-free (κ²g = 7/8 is an exact identity on the
  whole locked-transport family — it could not discriminate anything). 35/24 is a real,
  falsifiable, closed-form number. Price: the bench must declare its branch (restricted-ε
  platforms will and should see ≈7/8·[1 − c₄ε^{2/3}]).
- **V-diagnostic: BETTER (more discriminating), with branch-ID risk.** 3/2 is 9σ from
  the archived restricted-branch measurement — by design, since it diagnoses the *other*
  branch; a bench that confuses branches manufactures a false kill.
- **S5/P8-6: UNCHANGED.**
- **P8-7: SHIFTED, BETTER** — target moves −7.4% to an exact closed form (1/√2).
- **Over-spin: BETTER and riskier** — one rehearsed number becomes a four-rung exact
  ladder, the corpus's own measured onset 1.000 ± 0.004 is reproduced exactly (polar
  channel), and the marginal saturation point makes "immediate shedding at κ = 1/√2"
  a sharp new prediction.
- **The one forfeiture (WORSE):** the ε-scan anchor. v2.0.1 owned a 0.01σ hit
  (2.371 vs 2.37 ± 0.09); under the repair that hit is re-scoped to code validation of a
  saddle, and the replacement law (𝔠 rising, +1.15ε; ceiling 2.6×10⁻³) is **unmeasured** —
  a standing obligation (mitigated only by F-R15's transit reading).
- **The one thinned margin (NEWLY-TENSIONED, lite):** Q-6′ quark censorship — the
  ε_q-to-ceiling margin drops from ~×17 ("twenty-fold") to **×2.6–5.3**. The exact
  confinement dichotomy survives, but the (4.9′) ceiling re-derivation should print this.

**Net, within-model:** the repaired closure leaves the corpus's empirical exposure
**BETTER on the bench program** (three stakes become exact and discriminating, one gains
a richer signature), **UNCHANGED across the entire data-confronting sector** (every
V15.3–V15.9 number bit-identical or ≤ 0.12σ; S1/S2 verified untouched), and **WORSE at
exactly one point** — the forfeited ⟨r1⟩ ε-anchor corroboration — plus one thinned
censorship margin and the standing marginal-stability caveat. That is the trade §IX
should print: the repair costs the corpus its only closure-side *data hit* and pays it
back in adjudicability.

---

*Files: `i2_propagate.py` (34/34 checks), `i2_results.json`, this memo.
Cross-links: I1 (⟨r5⟩ 𝔠_q = 3.1 ± 0.2 vs saturated 3.2011, 0.51σ — the 𝔪_Sk caveat and
the headline corroboration question are the same question); the (4.8′)/(4.9′) re-scan
and the Majoron-battery re-run at g ≈ 7.8×10⁻¹⁰ are the two new arithmetic obligations
this propagation surfaces.*
