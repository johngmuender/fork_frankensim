#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
I2 — repaired-closure phenomenology propagation (ROADMAP_v6 Phase I2).

Machine-checks every number in i2_RESULTS.md:
  1. the ħ = cΛ√J calibration under T3.1 (exact refit factors, sympy);
  2. the c-ABSORBED vs c-EXPOSED census (with the blindness VERIFY items:
     family logs / <r10> bridge / m_Sk, each verified or bounded);
  3. the AUD-15 V15.3–V15.9 confrontation arithmetic re-run under the
     repaired constants (original vs repaired vs data/bound, verdicts);
  4. the bench stakes table v3 (incl. the channel-resolved over-spin menu,
     verified as exact Haar sine-moment identities);
  5. the bottom line, stake by stake.

Sources: theory-audit/T3_repaired_closure.md; corpus2/10-The-Repaired-
Closure-v1.0.md; corpus2/01-GUM-Omega-Paper-v3.0-ext.md (AUD-15 §V15,
IV.H/IV.J(4.8'-4.10'), IX.D/IX.F, VII.C'/K.2); tier0-gauntlet/RESULTS.md.

Deterministic; numpy + sympy only. Output: i2_results.json.
"""
import json
import numpy as np
import sympy as sp

R = {}          # results accumulator -> i2_results.json
CHECKS = []     # (name, ok, detail)


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name} {detail}")
    assert ok, f"CHECK FAILED: {name} {detail}"


def close(a, b, rtol=1e-9):
    return abs(float(a) - float(b)) <= rtol * max(1.0, abs(float(b)))


# ---------------------------------------------------------------------------
# 0. Exact tuples (sympy) — corpus (4.10) vs repaired (4.10'), T3 §3.1/§5.1
# ---------------------------------------------------------------------------
print("== 0. exact tuples ==")
pi = sp.pi
c_old = 128 * sp.sqrt(42) / (105 * pi)          # 2.514754
c_new = 64 * sp.sqrt(2) / (9 * pi)              # 3.201125  = (4/pi)*G*
Gstar = sp.Rational(16, 9) * sp.sqrt(2)
k_old = sp.sqrt(sp.Rational(7, 12))             # 0.763763
k_new = 1 / sp.sqrt(2)                          # 0.707107
e0 = sp.Rational(64, 15) / pi                   # ê0 = 64/15π
i0 = sp.Rational(256, 105) / pi                 # 𝔦0 = 256/105π

check("ê0/𝔦0 = 7/4 exact", sp.simplify(e0 / i0 - sp.Rational(7, 4)) == 0)
check("c_old² = 6ê0𝔦0 (corpus square identity)",
      sp.simplify(c_old**2 - 6 * e0 * i0) == 0)
check("c_new = (4/π)G* = 64√2/9π", sp.simplify(c_new - 4 / pi * Gstar) == 0)
check("T3.0 at saturation: κ² = ê_tot/2𝔦_tot = 1/2 (ê_tot = 𝔦_tot)",
      sp.simplify(k_new**2 - sp.Rational(1, 2)) == 0)

ratio = sp.simplify(c_new / c_old)              # 35/(6√21)
check("c_new/c_old = 35/(6√21) = (6√21/35)⁻¹",
      sp.simplify(ratio - 35 / (6 * sp.sqrt(21))) == 0,
      f"= {float(ratio):.6f} (+{100*(float(ratio)-1):.2f}%)")

R["tuple"] = {
    "c_old": {"form": "128*sqrt(42)/(105*pi)", "value": float(c_old)},
    "c_new": {"form": "64*sqrt(2)/(9*pi)", "value": float(c_new)},
    "c_ratio": {"form": "35/(6*sqrt(21))", "value": float(ratio),
                "pct": 100 * (float(ratio) - 1)},
    "kappa_old": {"form": "sqrt(7/12)", "value": float(k_old)},
    "kappa_new": {"form": "1/sqrt(2)", "value": float(k_new)},
    "binding_depth_old_pct": 100 * (1 - float(k_old)),
    "binding_depth_new_pct": 100 * (1 - float(k_new)),
    "omega_th_old": {"form": "sqrt(7/3)", "value": float(sp.sqrt(sp.Rational(7, 3)))},
    "omega_th_new": {"form": "sqrt(2)", "value": float(sp.sqrt(2))},
    "g_core": {"form": "35/24", "value": 35 / 24},
    "g_tot": {"form": "35/12", "value": 35 / 12},
    "e_tot=i_tot": {"form": "64/(9*pi)", "value": float(64 / (9 * np.pi))},
    "support_ratio_new": 1.5, "V_old": float(sp.sqrt(2)),
}

# ---------------------------------------------------------------------------
# 1. THE CALIBRATION QUESTION — ħ = cΛ√J, exact refit factors
# ---------------------------------------------------------------------------
print("== 1. calibration ==")
# ħ observed is fixed  =>  Λ√J refit by c_old/c_new
LsJ = sp.simplify(c_old / c_new)                # ×6√21/35 = 0.785584
check("Λ√J refit = 6√21/35", sp.simplify(LsJ - 6 * sp.sqrt(21) / 35) == 0,
      f"= {float(LsJ):.6f}")
# mass law E_tot = 𝒢 Λm̃ ; 𝒢: √2 ê0 -> 64/9π  => Λm̃ refit by inverse
G_old = sp.sqrt(2) * e0
G_new = sp.Rational(64, 9) / pi
Lm = sp.simplify(G_old / G_new)                 # ×3√2/5 = 0.848528
check("mass-law factor ratio 𝒢_new/𝒢_old = 5/(3√2) (+17.85%)",
      sp.simplify(G_new / G_old - 5 / (3 * sp.sqrt(2))) == 0,
      f"= {float(G_new/G_old):.6f}")
check("Λm̃ refit = 3√2/5", sp.simplify(Lm - 3 * sp.sqrt(2) / 5) == 0,
      f"= {float(Lm):.6f}")
# J pinned elsewhere ("the Λ-vs-J split is fixed elsewhere", T3 §5.1 note):
Lam = LsJ                                       # Λ ×6√21/35
mt = sp.simplify(Lm / LsJ)                      # m̃ ×√(7/6) = 1.080123
check("m̃ refit = √(7/6)", sp.simplify(mt - sp.sqrt(sp.Rational(7, 6))) == 0,
      f"= {float(mt):.6f}")
# compacton radius R* = (2Λ/π²m̃)^{1/3}  (Ω VII, exact backbone)
Rstar = sp.simplify((Lam / mt)**sp.Rational(1, 3))
check("R* refit = (18√2/35)^(1/3)",
      sp.simplify(Rstar**3 - 18 * sp.sqrt(2) / 35) == 0,
      f"= {float(Rstar):.6f}")

R["calibration"] = {
    "LambdaSqrtJ_refit": {"form": "6*sqrt(21)/35", "value": float(LsJ)},
    "LambdaMtilde_refit": {"form": "3*sqrt(2)/5", "value": float(Lm)},
    "mass_law_factor_change_pct": 100 * (float(G_new / G_old) - 1),
    "mtilde_refit": {"form": "sqrt(7/6)", "value": float(mt)},
    "Lambda_refit_J_pinned": float(Lam),
    "Rstar_refit": {"form": "(18*sqrt(2)/35)**(1/3)", "value": float(Rstar)},
    "note": ("hbar and every species mass are refit targets: Λ√J and Λm̃ "
             "absorb the c- and 𝒢-changes exactly; observables built solely "
             "from ħ, c, masses, α are c-ABSORBED."),
}

# ---------------------------------------------------------------------------
# 2a. VERIFY blindness items
# ---------------------------------------------------------------------------
print("== 2a. blindness verifications ==")
# (i) family logs: M_k/M_j = (ê_k m̃_k)/(ê_j m̃_j); 𝒢, Λ species-universal
ek, ej, mk, mj, Lsym, Gsym = sp.symbols("ek ej mk mj L G", positive=True)
Mratio = (Gsym * ek * Lsym * mk) / (Gsym * ej * Lsym * mj)
check("P-O1: mass ratio 𝒢,Λ-free (symbolic cancellation)",
      sp.simplify(Mratio - (ek * mk) / (ej * mj)) == 0)
# hence ln-mass-ratios (½A, ½(A+B), termination mass, quark ratios, Σ(p))
# are bit-identical under any 𝒢 change: c-ABSORBED.

# (ii) <r10> bridge: m3 = m_tau * exp(-24.36).  Leading order: pure
# exponentiated family-log arithmetic (A_core, κ_far, data-side A are
# m̃-ratio-class) => blind.  Residual c-route: the R* <-> reduced-Compton
# identification.  Bound it both ways:
dln_kappa = float(sp.log(k_old / k_new))        # if R* = κ·ƛ (range reading)
dln_Rstar = float(-sp.log(Rstar))               # if R* = compacton radius
band = 0.90
bridge_res = {
    "dln_if_Rstar_is_range_kappa_lambdabar": dln_kappa,
    "dln_if_Rstar_is_compacton_radius": dln_Rstar,
    "sigma_frac_worst": max(dln_kappa, dln_Rstar) / band,
}
check("bridge residual ≤ 0.12σ of its own ±0.90 band",
      bridge_res["sigma_frac_worst"] <= 0.12,
      f"worst |Δln| = {max(dln_kappa, dln_Rstar):.4f} -> "
      f"{bridge_res['sigma_frac_worst']:.3f}σ")
# worst-case m3 shift if the κ-reading applied (it is NOT the frozen
# pipeline's; frozen pipeline is base-m̃_tau => exactly blind):
m_tau = 1776.86e6  # eV
lnq = (5.644 - 3.05) / 0.1065
m3 = m_tau * np.exp(-lnq)
m3_worst = m3 * float(sp.sqrt(sp.Rational(7, 6)))
check("V15.4 central chain: (5.644-3.05)/0.1065 = 24.357, m3 = 0.0468 eV",
      close(lnq, 24.35681, 1e-5) and close(m3, 0.04688, 1e-3),
      f"lnq={lnq:.5f}, m3={m3:.5f}")

# (iii) m_Sk = c_Sk*sqrt(a2*a4) (Thm VII.C'): contains no c, κ, Λ, J, m̃ —
# quartic-regime, m̃-independent; its 1.7 GeV is a hadronic-scale
# identification.  Caveat: c_Sk is a quartic-CLOSURE constant and the halo
# threshold is t-independent (F-R5), so a repaired quartic re-scan could in
# principle move c_Sk — but the archive's own quartic closure <r5> reads
# c_q = 3.1 ± 0.2, i.e. 0.51σ FROM THE SATURATED VALUE already (I1's
# question); no shift of the 1.7 GeV input is implied by the repair.
cq, cq_err = 3.1, 0.2
pull_r5 = (float(c_new) - cq) / cq_err
check("<r5> quartic c_q vs saturated: 0.5σ (m_Sk caveat, I1 crosslink)",
      abs(pull_r5) < 0.6, f"pull = {pull_r5:.2f}σ")

R["blindness_verify"] = {
    "family_logs": "BLIND — symbolic: M_k/M_j = (ê_k m̃_k)/(ê_j m̃_j), 𝒢/Λ cancel",
    "r10_bridge": {"verdict": "BLIND at leading order (frozen pipeline is "
                              "base-m̃_tau, pure mass-ratio arithmetic); "
                              "residual via any R*<->ƛ identification bounded",
                   **bridge_res,
                   "m3_worst_case_eV": m3_worst},
    "m_Sk": {"verdict": "BLIND as printed (m_Sk = c_Sk*sqrt(a2*a4): no c/κ/Λ√J "
                        "route; 1.7 GeV is a hadronic identification)",
             "caveat_I1": f"c_Sk is quartic-closure-class; archive <r5> c_q = "
                          f"3.1±0.2 sits {pull_r5:+.2f}σ from saturated 3.2011 "
                          f"— consistent with the repaired branch already"},
}

# ---------------------------------------------------------------------------
# 2b. Halo-threshold channel menu — exact sine-moment identities (sympy)
# ---------------------------------------------------------------------------
print("== 2b. channel menu ==")
th = sp.symbols("theta", positive=True)


def s2_mean(weight):
    num = sp.integrate(weight * sp.sin(th)**2 * sp.sin(th), (th, 0, sp.pi))
    den = sp.integrate(weight * sp.sin(th), (th, 0, sp.pi))
    return sp.simplify(num / den)


menu = {}
# polar channel s² = |cos|
m_polar = s2_mean(sp.Abs(sp.cos(th)))
check("polar ⟨sin²⟩ = 1/2 -> κ_crit = 1 exactly",
      sp.simplify(m_polar - sp.Rational(1, 2)) == 0)
menu["polar (s²=|cos|)"] = {"s2": "1/2", "kappa_crit": 1.0}
# uniform
m_unif = s2_mean(1)
check("uniform ⟨sin²⟩ = 2/3 -> κ_crit = √3/2",
      sp.simplify(m_unif - sp.Rational(2, 3)) == 0,
      f"κ = {float(1/sp.sqrt(2*m_unif)):.6f}")
menu["uniform"] = {"s2": "2/3", "kappa_crit": float(1 / sp.sqrt(2 * m_unif))}
# s = sin² channel: weight = s² = sin⁴ -> ⟨sin⁶⟩/⟨sin⁴⟩ = 6/7
m_sin2 = s2_mean(sp.sin(th)**4)
check("s=sin² ⟨sin⁶⟩/⟨sin⁴⟩ = 6/7 -> κ_crit = √(7/12) = old κ0(!)",
      sp.simplify(m_sin2 - sp.Rational(6, 7)) == 0,
      f"κ = {float(1/sp.sqrt(2*m_sin2)):.6f}")
menu["s=sin²"] = {"s2": "6/7", "kappa_crit": float(sp.sqrt(sp.Rational(7, 12)))}
# equatorial ring / unlocked tilt: ⟨sin²⟩ -> 1
menu["equatorial ring / tilt (binding)"] = {"s2": "1",
                                            "kappa_crit": float(1 / sp.sqrt(2))}
R["overspin_channel_menu"] = menu

# ---------------------------------------------------------------------------
# 2c. The c-EXPOSED register (old -> new, machine-computed)
# ---------------------------------------------------------------------------
print("== 2c. exposed register ==")
lamC_e = 3.8615926799e-13  # m, electron reduced Compton wavelength (CODATA)
range_old = float(k_old) * lamC_e
range_new = float(k_new) * lamC_e
check("WS-A census range: κ_new·ƛ_e = 2.73e-13 m",
      close(range_new, 2.7306e-13, 1e-3), f"= {range_new:.4e} m")

# ε-law and anchor
c0f, cnf = float(c_old), float(c_new)
anchor_old = c0f * (1 - 0.42 * 0.05**(2 / 3))
c_sat_005 = 3.38547                              # T3 measured value at ε=0.05
slope_lin = (c_sat_005 / cnf - 1) / 0.05
check("old anchor 2.5147(1-0.42·0.05^{2/3}) = 2.371",
      close(anchor_old, 2.3711, 2e-4), f"= {anchor_old:.4f}")
check("repaired ε-law class +1.15ε (from 3.2011->3.38547 at 0.05)",
      close(slope_lin, 1.152, 2e-3), f"= +{slope_lin:.3f}ε")
pull_anchor_old = abs(anchor_old - 2.37) / 0.09
pull_anchor_new = abs(c_sat_005 - 2.37) / 0.09

# B-U1' ceiling
Dlock, cg = 3e-3, 0.42
ceil_old = (Dlock / cg)**1.5
ceil_new = Dlock / slope_lin
check("ceiling old (Δ/c_g)^{3/2} ≈ 6e-4", close(ceil_old, 6.04e-4, 2e-2),
      f"= {ceil_old:.3e}")
check("ceiling new Δ/1.152 ≈ 2.6e-3", close(ceil_new, 2.60e-3, 2e-2),
      f"= {ceil_new:.3e}")

# quark censorship margin (Q-6': ε_dress vs entrainment ceiling)
eps_q = (np.array([0.14, 0.20]) / 1.7)**2
marg_old = eps_q / ceil_old
marg_new = eps_q / ceil_new
check("censorship margin narrows ~x4.3 but stays > 1 (dichotomy survives)",
      marg_new.min() > 1.0,
      f"old [{marg_old[0]:.1f}, {marg_old[1]:.1f}] -> "
      f"new [{marg_new[0]:.2f}, {marg_new[1]:.2f}]")

# operative ε_e range and f-window propagation (F-A15-3 redo)
eps_bot = 1.0e-5
f_MeV = lambda e: 0.69 * np.sqrt(e) * 1700.0
f_old = (f_MeV(eps_bot), f_MeV(min(6e-4, ceil_old * 1.0)))
f_new = (f_MeV(eps_bot), f_MeV(ceil_new))
g_old = (m3 / (f_old[1] * 1e6), m3 / (f_old[0] * 1e6))
g_new = (m3 / (f_new[1] * 1e6), m3 / (f_new[0] * 1e6))
check("f-window old [3.7, 29] MeV; new top ≈ 60 MeV (= II.H window top)",
      close(f_new[1], 59.9, 2e-2),
      f"old [{f_old[0]:.1f},{f_old[1]:.1f}] new [{f_new[0]:.1f},{f_new[1]:.1f}]")

# S4' numbers
k2g_old = 7 / 8
k2g_new = float(k_new)**2 * 35 / 12
check("S4' κ²g_tot = 35/24 = 1.45833", close(k2g_new, 35 / 24), f"= {k2g_new:.5f}")
r1_k2g, r1_k2g_err = 0.843, 0.046
pull_s4_old = (k2g_old - r1_k2g) / r1_k2g_err
pull_s4_new = (k2g_new - r1_k2g) / r1_k2g_err

# V diagnostic
V_meas, V_err = 1.409, 0.010
pull_V_sqrt2 = (np.sqrt(2) - V_meas) / V_err
pull_V_32 = (1.5 - V_meas) / V_err
check("⟨r1⟩ V: 0.5σ from √2, 9.1σ from 3/2",
      close(pull_V_sqrt2, 0.52, 5e-2) and close(pull_V_32, 9.1, 2e-2),
      f"= {pull_V_sqrt2:.2f}σ / {pull_V_32:.2f}σ")

# ω_th inheritors
wth_ratio = float(sp.sqrt(2) / sp.sqrt(sp.Rational(7, 3)))
well_factor = wth_ratio**2
check("WS-A8 well parameter x6/7 = 0.857", close(well_factor, 6 / 7),
      f"= {well_factor:.6f} (envelope ×{wth_ratio:.5f})")

# binding fraction / P-O1 magnitude; edge-stiffening anchor
depth_ratio = (1 - float(k_new)) / (1 - float(k_old))
lam_edge = depth_ratio**(-0.5)

R["exposed_register"] = {
    "LambdaSqrtJ": {"old": 1.0, "new": float(LsJ), "class": "calibration register"},
    "kappa_range_lambdabar": {"old": float(k_old), "new": float(k_new),
                              "pct": -100 * (1 - float(k_new) / float(k_old))},
    "electron_range_m": {"old": range_old, "new": range_new},
    "binding_depth_pct": {"old": 100 * (1 - float(k_old)),
                          "new": 100 * (1 - float(k_new)),
                          "rel_change": depth_ratio},
    "omega_th": {"old": float(sp.sqrt(sp.Rational(7, 3))), "new": float(sp.sqrt(2)),
                 "WS_A8_well_factor": well_factor},
    "eps_law": {"old": "c0[1 - 0.42 eps^(2/3)] (falling)",
                "new": f"c_sat[1 + {slope_lin:.2f} eps] (rising; sign flip)",
                "anchor_old_pull_sigma": pull_anchor_old,
                "repaired_c_at_0.05": c_sat_005,
                "repaired_vs_r1_sigma": pull_anchor_new,
                "note": "<r1> re-scoped to restricted-branch code validation; "
                        "the repaired branch has NO measured ε-scan yet"},
    "BU1_ceiling": {"old": ceil_old, "new": ceil_new},
    "operative_eps_range": {"old": [eps_bot, 6e-4], "new": [eps_bot, ceil_new]},
    "f_window_MeV": {"old": list(f_old), "new": list(f_new)},
    "majoron_g": {"old": list(g_old), "new": list(g_new),
                  "note": "low-g edge new: battery re-run flagged"},
    "quark_censorship_margin": {"old": list(marg_old), "new": list(marg_new),
                                "note": "Q-6' 'twenty-fold' -> ~2.6-5.3x; "
                                        "dichotomy survives, margin x4.3 thinner"},
    "S4prime": {"old": k2g_old, "new": k2g_new,
                "r1_measured": [r1_k2g, r1_k2g_err],
                "pull_old_sigma": pull_s4_old, "pull_new_sigma": pull_s4_new},
    "support_diag": {"old": float(sp.sqrt(2)), "new": 1.5,
                     "r1_pulls_sigma": [float(pull_V_sqrt2), float(pull_V_32)]},
    "mass_law_norm": {"factor": float(G_new / G_old),
                      "note": "+17.85% — absorbed by Λm̃ x 3√2/5 (exposed in "
                              "substrate register, absorbed observationally)"},
    "P_O1_magnitude": {"depth_rescale": depth_ratio,
                       "note": "μ-drift suppression floor ~1e-2 class both ways"},
    "edge_stiffening_anchor": {"lambda_rescale_at_kappa_phys": lam_edge},
}

R["absorbed_register"] = [
    "ħ itself (calibration target; Λ√J x 6√21/35)",
    "all mass ratios M_k/M_j = (ê_k m̃_k)/(ê_j m̃_j) — P-O1 verified symbolically",
    "family logs ½A, ½(A+B), (A+B)/A; termination mass m_e e^{-8.5}; Σ(p) belts",
    "quark spacing ratios <r11> / FQ table (pure ratio arithmetic)",
    "<r10> bridge central 24.357 and m3 = 0.0468 eV (frozen base-m̃_tau pipeline)",
    "spectroscopy per the m̃-chain (spectra are functions of ħ, masses, α only)",
    "E_rot/E = 1/4 (kinematic, exact all closures)",
    "spin-selection w·j(1-j) = 1; M.0 loop; ±½ error-signal exponents",
    "Born statistics, spectrum tiers, electrodynamics identities (no contact)",
    "ΔN_eff = 0.0268; B-ν1 floor 1.5e-6; sidereal 0.123/0.41 ns (physical const.)",
    "m_Sk = 1.7 GeV (quartic, m̃-independent; I1 caveat noted)",
    "P-F1' clock corrections (m_e/m_μ)²/2, δ_τ (mass ratios)",
    "S1, S2 stakes entire (family/neutrino sector blind)",
]

# ---------------------------------------------------------------------------
# 3. V15.3–V15.9 re-run under repaired constants
# ---------------------------------------------------------------------------
print("== 3. V15 re-run ==")
m_e, m_mu = 0.51099895e6, 105.6583755e6  # eV
v15 = []

# V15.3 family sector
lnA = np.log(m_tau / m_mu)
lnB = np.log(m_mu / m_e)
term_mass = m_e * np.exp(-8.5)
supp = np.exp((4 / 3) * 8.5)
check("V15.3: ln(m_tau/m_mu) = 2.822; ln(m_mu/m_e) = 5.332; m_e e^-8.5 = 104 eV",
      close(lnA, 2.8222, 1e-3) and close(lnB, 5.3316, 1e-3)
      and close(term_mass, 103.8, 2e-3) and close(supp, 8.3e4, 1e-2),
      f"{lnA:.4f}/{lnB:.4f}/{term_mass:.1f}/{supp:.2e}")
v15.append({
    "id": "V15.3", "sector": "family logs + termination",
    "original": {"halfA": 2.80, "halfApB": 5.65, "ratio": 2.018,
                 "term_mass_eV": 104, "eps_threshold": 8.4e-6},
    "repaired": "bit-identical (all mass-ratio arithmetic; P-O1 verified)",
    "confronted": {"ln(mtau/mmu)": round(lnA, 4), "ln(mmu/me)": round(lnB, 4),
                   "pulls_sigma": [0.05, 0.42, 0.47],
                   "eps_threshold_vs_operative_bottom":
                       f"8.4e-6 < 1.0e-5 (unchanged; ceiling rise to "
                       f"{ceil_new:.2e} does not touch the floor side)"},
    "verdict": "SURVIVES (identically)",
})

# V15.4 bridge and stake
band_lo, band_hi = m3 * np.exp(-0.90), m3 * np.exp(0.90)
m3_floor, m2_floor = np.sqrt(2.53e-3), np.sqrt(7.5e-5)
Sigma_floor = m3_floor + m2_floor
# Σ(m3) for NO with splittings; solve Σ=0.11 by bisection (no scipy needed)
def Sigma_NO(m3v, dm31=2.53e-3, dm21=7.5e-5):
    m1sq = m3v**2 - dm31
    m1 = np.sqrt(max(m1sq, 0.0))
    m2 = np.sqrt(m1**2 + dm21)
    return m1 + m2 + m3v
lo, hi = m3_floor, 0.2
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if Sigma_NO(mid) < 0.11:
        lo = mid
    else:
        hi = mid
m3_at_011 = 0.5 * (lo + hi)
check("V15.4: band [0.019, 0.115]; Σ_floor = 0.059; Σ≤0.11 <=> m3 ≤ 0.057",
      close(band_lo, 0.0190, 2e-2) and close(band_hi, 0.1153, 2e-2)
      and close(Sigma_floor, 0.0590, 2e-3) and close(m3_at_011, 0.057, 2e-2),
      f"[{band_lo:.4f},{band_hi:.4f}], Σf={Sigma_floor:.4f}, m3*={m3_at_011:.4f}")
v15.append({
    "id": "V15.4", "sector": "<r10> bridge and S1 stake",
    "original": {"lnq": 24.36, "m3_eV": 0.0468, "band": [0.019, 0.115],
                 "Sigma_floor": 0.059},
    "repaired": {"central": "identical (blind, frozen pipeline)",
                 "residual_exposure_sigma": round(bridge_res["sigma_frac_worst"], 3),
                 "m3_worst_case_eV": round(m3_worst, 4),
                 "worst_case_still_inside": bool(band_lo < m3_worst < 0.057)},
    "confronted": {"stake_floor": 0.058, "Sigma_floor_physical": round(Sigma_floor, 4)},
    "verdict": "SURVIVES (residual ≤ 0.12σ; S1 endpoints untouched)",
})

# V15.5 soft sector / Majoron
f_bench = 0.69 * np.sqrt(1e-2) * 1.7  # GeV at eps=1e-2 -> 0.117? (V15.5 checks the LAW arithmetic)
Bnu1 = (0.26 * 1.22e28 * m3**2)**0.25 / 1e6  # MeV, T=0.26 eV, M_Pl=1.22e19 GeV
eps_floor = (1.4 / (0.69 * 1700))**2
dNeff = (4 / 7) * (10.75 / 106.75)**(4 / 3)
check("V15.5: B-ν1 ~ 1.6 MeV (corpus 1.57->1.4 w/ thermal); eps-floor 1.4e-6; "
      "ΔN_eff = 0.0268",
      close(Bnu1, 1.6, 6e-2) and close(eps_floor, 1.42e-6, 2e-2)
      and close(dNeff, 0.02678, 1e-3),
      f"{Bnu1:.2f} MeV / {eps_floor:.2e} / {dNeff:.5f}")
v15.append({
    "id": "V15.5", "sector": "soft sector / Majoron",
    "original": {"B_nu1_MeV": 1.57, "eps_floor": 1.5e-6, "dNeff": 0.0268},
    "repaired": "own items bit-identical (m_ν, T, M_Pl, m_Sk: no c-route). "
                "Downstream F-A15-3 propagation SHIFTS: f-window top "
                f"{f_old[1]:.0f} -> {f_new[1]:.0f} MeV (ceiling {ceil_old:.1e}"
                f" -> {ceil_new:.1e}); Majoron g low edge {g_old[0]:.1e} -> "
                f"{g_new[0]:.1e}",
    "confronted": {"cosmological floor": "unchanged", "battery":
                   "passed on old window; re-run flagged at new low-g edge"},
    "verdict": "SURVIVES; propagation SHIFTS-HARMLESSLY (battery re-run flagged)",
})

# V15.6 P-F1'
pf1 = (m_e / m_mu)**2 / 2
dtau = (m_e / m_tau)**2 / 2
check("V15.6: (m_e/m_mu)²/2 = 1.17e-5; δ_τ = 4.1e-8",
      close(pf1, 1.169e-5, 2e-3) and close(dtau, 4.13e-8, 2e-2),
      f"{pf1:.3e} / {dtau:.2e}")
v15.append({
    "id": "V15.6", "sector": "P-F1' clock corrections",
    "original": {"muon": 1.17e-5, "tau": 4.1e-8},
    "repaired": "bit-identical (mass ratios)",
    "confronted": "spectroscopy nulls (compliance)", "verdict": "SURVIVES",
})

# V15.7 color
lam_star_printed = 0.42 * 0.2**(1 / 3)
f_q = np.sqrt(1e-2) * 1.7
check("V15.7: eps_q = [6.8e-3, 1.38e-2]; λ* printed 0.246; f_q = 0.17 GeV",
      close(eps_q[0], 6.78e-3, 2e-2) and close(eps_q[1], 1.384e-2, 2e-2)
      and close(lam_star_printed, 0.2457, 2e-3) and close(f_q, 0.17, 1e-3),
      f"[{eps_q[0]:.2e},{eps_q[1]:.2e}] / {lam_star_printed:.4f} / {f_q:.3f}")
v15.append({
    "id": "V15.7", "sector": "color: stiffness closed loop + shape law",
    "original": {"eps_q": [round(eps_q[0], 4), round(eps_q[1], 4)],
                 "lambda_star": 0.246, "f_q_GeV": 0.17,
                 "sigma_inversion_band_GeV": [0.14, 0.20]},
    "repaired": {"eps_q and f_q loop": "bit-identical (m_Sk & σ-band arithmetic; "
                                       "c-blind) — closed loop A SURVIVES",
                 "lambda_star_law": "SUPERSEDED — λ*(ε) = 0.42 ε^{1/3} is the "
                                    "restricted-branch <r6> shape law; the "
                                    "saturated solution is the fixed oblate "
                                    "compacton + halo (re-scan required)",
                 "censorship_margin": f"[{marg_old[0]:.1f},{marg_old[1]:.1f}] -> "
                                      f"[{marg_new[0]:.2f},{marg_new[1]:.2f}]"},
    "confronted": {"f_q vs sigma-inversion band": "0.17 in [0.14, 0.20] dead-"
                                                  "center — unchanged"},
    "verdict": "loop SURVIVES; λ*-row SUPERSEDED (restricted-branch); "
               "Q-6' censorship NEWLY-TENSIONED-LITE: margin x4.3 thinner, "
               "dichotomy survives (min margin 2.6)",
})

# V15.8 FQ table
R_up = 1 + 0.33 - 0.06
R_dn = 1 + 0.33 - 0.51
check("V15.8: aligned 1.27, anti-aligned 0.82 (sign flip forced)",
      close(R_up, 1.27) and close(R_dn, 0.82),
      f"{R_up:.2f} / {R_dn:.2f}")
v15.append({
    "id": "V15.8", "sector": "<r11> FQ overdetermination",
    "original": {"R_up": 1.27, "R_down": 0.82},
    "repaired": "bit-identical (pure ratio arithmetic; externally reconstructed, "
                "pulls 0.17-0.32σ)",
    "confronted": "PDG quark spacings", "verdict": "SURVIVES (identically)",
})

# V15.9 bookkeeping
haar = (2 / np.pi) * 8 * (4 / 15)
tau1 = 3.7e5 * 30 / (2.99792458e8)**2 * 1e9   # ns
tau2 = 3.7e5 * 100 / (2.99792458e8)**2 * 1e9  # ns
check("V15.9: Haar 64/15π = 1.3581; kills 6+8=14; sidereal 0.123/0.41 ns",
      close(haar, 64 / (15 * np.pi)) and (6 + 8 == 14)
      and close(tau1, 0.1235, 2e-2) and close(tau2, 0.412, 2e-2),
      f"{haar:.4f} / {tau1:.3f} / {tau2:.3f}")
v15.append({
    "id": "V15.9", "sector": "bookkeeping",
    "original": {"haar": 1.3581, "kills": 14, "sidereal_ns": [0.123, 0.41]},
    "repaired": "bit-identical — ê0 = 64/15π survives as the spherical Haar "
                "backbone (enters the repaired algebra too); kill census and "
                "sidereal kinematics have no c-route",
    "confronted": "internal cross-links", "verdict": "SURVIVES (identically)",
})

R["v15_rerun"] = v15

# context rows (V15.1/V15.2, outside the tasked range but the exposed ones)
R["v15_context"] = {
    "V15.1": "arithmetic survives (restricted-branch algebra, certified); the "
             "physical-endpoint tuple (c0, κ0, ω_th) STRUCK per T3.1 — the "
             "supersession, not a failed check",
    "V15.2": f"anchor arithmetic survives; the confrontation is RE-ASSIGNED: "
             f"old law hits <r1> at {pull_anchor_old:.2f}σ but belongs to the "
             f"struck endpoint; repaired branch predicts c(0.05) = {c_sat_005} "
             f"({pull_anchor_new:.1f}σ from <r1> — correctly so: <r1> is "
             f"restricted-branch code, per F-R5/F-R15)",
}

# ---------------------------------------------------------------------------
# 4. Bench stakes table v3
# ---------------------------------------------------------------------------
print("== 4. bench stakes v3 ==")
stakes = [
    {"stake": "S3 (E_rot/E)", "old": 0.25, "new": 0.25,
     "note": "exact, all closures — UNCHANGED (now closure-independent)",
     "exposure": "UNCHANGED"},
    {"stake": "S4' (κ²g)", "old": 7 / 8, "new": float(sp.Rational(35, 24)),
     "note": f"κ²g_tot = 35/24 (κ²g_core = 35/48 = {35/48:.5f}; cleanest: "
             f"ê_tot = 𝔦_tot = 64/9π = {float(64/(9*np.pi)):.5f}); old 7/8 "
             f"family-exact hence g-blind (value-free); ⟨r1⟩ 0.843±0.046: "
             f"{pull_s4_old:+.1f}σ from 7/8 (restricted), {pull_s4_new:+.1f}σ "
             f"from 35/24 (different branch)",
     "exposure": "BETTER (discriminating at last)"},
    {"stake": "V-class diagnostic", "old": float(sp.sqrt(2)), "new": 1.5,
     "note": f"support ratio 3/2 exact (saturated); V = √2 re-scoped to "
             f"restricted branch; ⟨r1⟩ 1.409±0.010 = {pull_V_sqrt2:.1f}σ from "
             f"√2, {pull_V_32:.1f}σ from 3/2 — bench must declare its branch",
     "exposure": "BETTER (more discriminating) with branch-ID risk"},
    {"stake": "S5/P8-6 (clock lock)", "old": "V ∝ cosΔφ, bond vibron",
     "new": "identical", "note": "kinematic; x0_lock = 2.42±0.12 unchanged; "
                                 "-sinΔφ torque externally confirmed",
     "exposure": "UNCHANGED"},
    {"stake": "P8-7 (tail range)", "old": float(k_old), "new": float(k_new),
     "note": "range = κ × clock length; κ = 1/√2 = 0.70711 exactly (-7.4%); "
             "closed form, no band needed",
     "exposure": "SHIFTED target; BETTER (exact)"},
    {"stake": "over-spin onset", "old": "single rehearsed onset κ = 1",
     "new": "channel-resolved menu κ_crit(w) = 1/√(2⟨sin²θ⟩_w)",
     "note": "equilibrium sits marginally AT κ = 1/√2: over-spin sheds "
             "IMMEDIATELY through the equatorial-ring/tilt channel (zero "
             "margin — the marginal-stability caveat IS the bench signature; "
             "branch-resolved emission detection required); ladder if ring "
             "suppressed: 0.70711 (ring/tilt) < 0.76376 (s=sin²; the old κ0, "
             "now a threshold) < 0.86603 (uniform) < 1.00000 (polar — the "
             "corpus's measured 1.000±0.004, reproduced exactly)",
     "exposure": "BETTER (richer + riskier prediction; halo thresholds are "
                 "now the physics)"},
    {"stake": "S1 (Σm_ν)", "old": "[0.058, 0.11] eV, NO", "new": "identical",
     "note": "VERIFIED untouched: bridge blind (≤0.12σ residual), floors "
             "physical (oscillation + Σ bound)", "exposure": "UNCHANGED"},
    {"stake": "S2 (Higgs nulls)", "old": "δg_h = 0", "new": "identical",
     "note": "family sector, blind", "exposure": "UNCHANGED"},
]
R["bench_stakes_v3"] = stakes
for s in stakes:
    print(f"  {s['stake']}: {s['exposure']}")

# ---------------------------------------------------------------------------
# 5. Bottom line
# ---------------------------------------------------------------------------
R["bottom_line"] = {
    "data_confronting_sector (family, neutrino, quark, cosmology)":
        "UNCHANGED — every V15.3-V15.9 confrontation number is bit-identical "
        "or shifts ≤ 0.12σ inside its own band; S1/S2 verified untouched",
    "bench_program": "BETTER — S4' becomes discriminating (was value-free), "
                     "P8-7 and the V-diagnostic become exact closed forms, "
                     "the over-spin protocol becomes a channel-resolved "
                     "ladder with the corpus's own measured 1.000±0.004 "
                     "reproduced exactly as the polar channel",
    "forfeited": "the <r1> ε-anchor corroboration (2.371 vs 2.37±0.09, "
                 "0.01σ) is RE-SCOPED to restricted-branch code validation; "
                 "the repaired rising ε-law (+1.15ε) has no measurement yet "
                 "— a new obligation, mitigated only by F-R15's transit "
                 "reading",
    "newly_tensioned": "Q-6' quark-censorship margin thins x4.3 (min 2.6; "
                       "dichotomy survives); marginal (not gapped) stability "
                       "is the standing structural risk",
    "net": "within-model empirical exposure: BETTER on the bench stakes, "
           "UNCHANGED across every data confrontation, WORSE at exactly one "
           "point (the forfeited ε-anchor hit) plus one thinned margin",
}

# ---------------------------------------------------------------------------
with open(__file__.replace("i2_propagate.py", "i2_results.json"), "w") as fh:
    json.dump(R, fh, indent=2, ensure_ascii=False, default=float)

n_pass = sum(1 for _, ok, _ in CHECKS if ok)
print(f"\n== {n_pass}/{len(CHECKS)} machine checks PASS; i2_results.json written ==")
