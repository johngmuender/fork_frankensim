#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
J2 — the repaired (4.8')/(4.9') eps-scan arithmetic (first arithmetic
obligation filed by the I2 propagation memo).

(a) (4.8'): the eps-dependence of the closure constants on the saturated
    branch (fixed oblate compacton + pinned marginal halo), across the
    corpus's printed eps range [1e-3, 0.1]; survive/shift/undefined census
    of every (4.8)-derived printed number (the forfeited eps-anchor).
(b) (4.9'): the censorship ceiling under repaired constants; the Q-6'
    margin per reading; the x4.3 thinning with provenance; the confinement
    dichotomy verdict at every reading + breach conditions.
(c) cross-check vs i2_propagate.py / i2_results.json where they overlap.

Deterministic, standalone (numpy + sympy + json + stdlib). Stage-flushed
output: j2_results.json. PASS/FAIL gate table printed at the end.
WITHIN-MODEL ONLY: nothing here bears on nature.
"""
import json
import os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
SBX = os.path.dirname(HERE)
OUT = os.path.join(HERE, "j2_results.json")

R = {"_notice": "WITHIN-MODEL ONLY; analysis layer; J2 obligation of I2"}
GATES = []


def gate(name, ok, detail=""):
    GATES.append({"gate": name, "pass": bool(ok), "detail": str(detail)})
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")


def flush():
    with open(OUT, "w") as fh:
        json.dump(R, fh, indent=1, ensure_ascii=False, default=float)


def close(a, b, rtol=1e-9):
    return abs(float(a) - float(b)) <= rtol * max(1.0, abs(float(b)))


def pin(path, substrings, name):
    """Text pin: every substring must occur verbatim in the file."""
    with open(path, encoding="utf-8") as fh:
        txt = fh.read()
    missing = [s for s in substrings if s not in txt]
    gate(f"PIN {name}", not missing,
         ("all substrings found" if not missing else f"MISSING: {missing}"))
    return not missing


# ===========================================================================
# 0. Exact constants + printed-text pins (provenance layer)
# ===========================================================================
print("== 0. exact constants and printed-text pins ==")
pi = sp.pi
c_old = 128 * sp.sqrt(42) / (105 * pi)      # corpus (4.10) endpoint, struck
c_sat = 64 * sp.sqrt(2) / (9 * pi)          # T3.1 saturated, (4.10')
Gstar = sp.Rational(16, 9) * sp.sqrt(2)
k_sat = 1 / sp.sqrt(2)
c_oldf, c_satf = float(c_old), float(c_sat)

gate("c_old = 128*sqrt(42)/(105*pi) = 2.514754 (printed '2.5147')",
     close(c_oldf, 2.5147536, 1e-6), f"= {c_oldf:.7f}")
gate("c_sat = (4/pi)*G* = 64*sqrt(2)/(9*pi) = 3.201125 (printed '3.2011')",
     sp.simplify(c_sat - 4 / pi * Gstar) == 0 and close(c_satf, 3.2011253, 1e-6),
     f"= {c_satf:.7f}")

P = os.path.join
v201 = P(SBX, "substrate-suite", "01-GUM-Omega-Paper-v2.0.1.md")
v30 = P(SBX, "corpus2", "01-GUM-Omega-Paper-v3.0-ext.md")
rep = P(SBX, "corpus2", "10-The-Repaired-Closure-v1.0.md")

pin(v201, ["𝔠(ε) = 𝔠₀[1 − c_g ε^{2/3}(1 + O(ε^{1/3}ln ε))]",
           "𝔠₀ = 128√42/(105π) = 2.5147, c_g = 0.42 ± 0.04",
           "ε ∈ [10⁻³, 0.1]",
           "2.515·[1 − 0.42·0.05^{2/3}] = 2.371"], "(4.8) as printed, v2.0.1")
pin(v201, ["ε_e ≤ (Δ_lock/c_g)^{3/2} ≈ 6×10⁻⁴",
           "ε_e ≥ 1.5×10⁻⁶", "admission window [1×10⁻⁵, 3×10⁻³]₆₈",
           "Δ_lock = λ_field/ω ≈ 3×10⁻³"], "(4.9)+Δ_lock as printed, v2.0.1")
pin(v30, ["3.2011 → 3.3855 at ε = 0.05 (+5.76%; ≈ +1.15ε class",
          "the rung's law is −1.03ε^{1.0}"], "(4.8') as printed, v3.0-ext")
pin(v30, ["the repaired linear, sign-flipped law gives ceiling ≈ 2.6×10⁻³",
          "the admission window survives numerically: [1.5×10⁻⁶, ~3×10⁻³]"],
    "(4.9') as printed, v3.0-ext")
pin(v30, ["ε_dress = (f_q/𝔪_Sk)² ≈ 10⁻²",
          "twenty-fold above the entrainment ceiling"],
    "Q-6' as printed, v3.0-ext (stale 'twenty-fold' retained — defect D1)")
pin(v30, ["λ_e ∼ ε_e^{1/3} ∼ 0.01–0.1"],
    "IV.J' consequences retains superseded shape law (defect D2)")
pin(v30, ["ε_q = (0.14–0.20/1.7)² = 0.7–1.4×10⁻²"], "V15.7 eps_q as printed")
pin(v30, ["f_q ≈ 0.14–0.20 GeV"], "Q-5 sigma-inversion band as printed")
pin(v30, ["ε_e ∈ [1.0×10⁻⁵, 6×10⁻⁴]₆₈"], "F-A15-3 old operative range")
pin(v30, ["saturated: +1.15ε class"], "App B.6 record, v3.0-ext")
pin(rep, ["censorship margin thins ×4.3", "the repaired (4.8′)/(4.9′)"],
    "obligation sentence, Repaired-Closure v1.0")

R["s0_constants"] = {"c_old": c_oldf, "c_sat": c_satf,
                     "kappa_sat": float(k_sat), "Gstar": float(Gstar)}
flush()

# ===========================================================================
# A. (4.8') — the eps-scan under the saturated closure
# ===========================================================================
print("== A. (4.8') saturated-branch eps arithmetic ==")
with open(P(SBX, "tier2-closure", "axi4_locked_results.json")) as fh:
    axi4 = json.load(fh)
ring0 = axi4["saturated"]["ring_t0"]      # eps -> 0 (t = 0)
ring5 = axi4["saturated"]["ring_eps"]     # eps = 0.05 (t = 0.0082764349)

# A.1 the two measured saturated points (the ONLY saturated-branch data)
c0_eng, c5_eng = ring0["c_paper"], ring5["c_paper"]
gate("A.1 saturated points: c(t0) = 3.203459, c(0.05) = 3.385470 (axi4)",
     close(c0_eng, 3.2034586, 1e-6) and close(c5_eng, 3.3854702, 1e-6),
     f"{c0_eng:.6f} / {c5_eng:.6f}")
bias0 = c0_eng / c_satf - 1
gate("A.2 engine t0 base vs exact c_sat: bias +7.3e-4 (defect D4 if > N-ladder)",
     abs(bias0) < 1e-3, f"bias = {bias0:+.2e} (N-ladder est. 1.4e-4 -> excess x5)")

# A.3 kappa / depth / omega_th / E_rot are eps-FROZEN on the saturated branch
k0, k5 = ring0["kappa_paper"], ring5["kappa_paper"]
gate("A.3 kappa_paper(t0) == kappa_paper(0.05) == 1/sqrt(2) (bit-identical)",
     k0 == k5 == float(k_sat) and ring0["erot_frac"] == ring5["erot_frac"] == 0.25,
     f"kappa = {k0!r} both; E_rot/E = 1/4 both. Saturation (T3.1 KKT) is "
     "eps-free analytically; axi4's eps=0.05 cap-ladder (c 3.05->3.29 rising, "
     "kappa 0.752->0.722 falling, NO interior optimum) shows the ring channel "
     "still binds at eps > 0")

# A.4 the slope: two readings (exact base vs engine base), nominal eps = 0.05
eps5 = 0.05           # nominal; locked-scan achieved 0.0499847 (0.03% effect)
rise_exact = c5_eng / c_satf - 1
rise_engine = c5_eng / c0_eng - 1
s_exact = rise_exact / eps5
s_engine = rise_engine / eps5
gate("A.4 rise at eps=0.05: +5.76% (printed '+5.76%'); slope +1.15 class",
     close(100 * rise_exact, 5.759, 2e-3) and close(s_exact, 1.152, 2e-3),
     f"rise = +{100*rise_exact:.3f}%, s_exact = {s_exact:.4f}, "
     f"s_engine = {s_engine:.4f} (1.4% spread = D4)")

# A.5 sign flip vs the printed falling law + the old anchor arithmetic
anchor_old = c_oldf * (1 - 0.42 * 0.05 ** (2 / 3))
r1_c, r1_e = 2.37, 0.09
pull_old = abs(anchor_old - r1_c) / r1_e
pull_new = abs(c5_eng - r1_c) / r1_e
gate("A.5 old anchor 2.5147*(1-0.42*0.05^{2/3}) = 2.3714, 0.016 sigma from <r1>",
     close(anchor_old, 2.371404, 1e-5) and close(pull_old, 0.0156, 2e-2),
     f"anchor = {anchor_old:.6f}, pull = {pull_old:.4f} sigma")
gate("A.6 forfeited anchor: repaired c(0.05) = 3.38547 is 11.28 sigma from <r1>",
     close(pull_new, 11.283, 1e-3),
     f"pull = {pull_new:.3f} sigma — correctly so (<r1> = restricted-branch "
     "code validation, F-R5/F-R15); the repaired branch has NO measured scan")

# A.7 the scan table across the printed range [1e-3, 0.1]
def c_old_law(e):
    return c_oldf * (1 - 0.42 * e ** (2 / 3))

def c_new_law(e):
    return c_satf * (1 + s_exact * e)

scan_rows = []
for e, grade in [(1e-3, "linear-interp (unmeasured)"),
                 (1e-2, "linear-interp (unmeasured)"),
                 (0.05, "MEASURED (axi4 ring-saturated)"),
                 (0.1, "linear-EXTRAP (unmeasured)")]:
    scan_rows.append({"eps": e, "c_old_falling": round(c_old_law(e), 5),
                      "c_sat_rising": round(c_new_law(e), 5),
                      "kappa_sat": float(k_sat), "grade": grade})
gate("A.7 scan-table endpoints: c_old(0.1) = 2.288 falling; c_sat(0.1) ~ 3.570",
     close(c_old_law(0.1), 2.2876, 1e-3) and close(c_new_law(0.1), 3.5698, 1e-3),
     "only eps={t->0, 0.05} are measured on the saturated branch — "
     "the multi-eps engine re-scan itself remains OPEN (D5)")

# A.8 g_tot runs with eps (the one other running constant); S4' drift
gt0, gt5 = ring0["g_total"], ring5["g_total"]
gate("A.8 g_tot: engine t0 = 2.91873 vs exact 35/12 = 2.91667 (+0.07%)",
     close(gt0, 35 / 12, 8e-4), f"g_tot(t0) = {gt0:.6f}")
s_gtot = (gt5 / gt0 - 1) / eps5
k2g_05 = 0.5 * gt5
gate("A.9 S4' eps-drift: kappa^2*g_tot(0.05) = 1.5008 vs stake 35/24 = 1.45833",
     close(k2g_05, 1.50076, 1e-3),
     f"g_tot(0.05) = {gt5:.5f} (+{100*(gt5/gt0-1):.2f}%; ~+{s_gtot:.2f}*eps); "
     "the discriminating stake is the eps->0 value — benches at eps~0.05 "
     "should expect ~1.50")

# A.10 repaired one-sided physical-constant bands (F-A15-1 redo, sign flipped)
D_lock = 3e-3
c_band_hi = c_satf * (1 + D_lock)
gate("A.10 (4.10') bands: drift strictly POSITIVE -> c_phys = 3.2011 +0.0096/-0;"
     " kappa/depth/omega_th EXACT (zero eps-band)",
     close(c_band_hi - c_satf, 0.009603, 1e-3),
     f"c_phys in [{c_satf:.4f}, {c_band_hi:.4f}]; old F-A15-1 bands were "
     "negative-side (2.515 -0.008/+0) — side flips with the law's sign")

R["sA_48prime"] = {
    "saturated_points": {"t0": {"c": c0_eng, "kappa": k0, "g_tot": gt0},
                         "eps0.05": {"c": c5_eng, "kappa": k5, "g_tot": gt5},
                         "engine_bias_t0": bias0,
                         "grid_ladder_rel": [1.4e-4, 6.7e-4]},
    "law": {"form": "c_sat*(1 + s*eps), leading order (sign FLIPPED vs (4.8))",
            "s_exact_base": s_exact, "s_engine_base": s_engine,
            "rise_at_0.05_pct": 100 * rise_exact,
            "exponent_status": "linear class ASSUMED (analytic leading order: "
                               "the saturated minimiser's oblateness is fixed "
                               "O(1), so the 2/3 fractional-power geometry "
                               "dies); measured at TWO eps points only (D3)"},
    "frozen_constants": {"kappa": float(k_sat), "binding_depth_pct": 29.289,
                         "omega_th_over_omega0": float(sp.sqrt(2)),
                         "erot_over_e": 0.25,
                         "note": "eps-frozen exactly on the saturated branch "
                                 "(saturation condition is eps-free); under "
                                 "(4.8) these ran with eps (<r1>: 0.802, 19.8%)"},
    "running_constants": {"c": f"+{s_exact:.3f}*eps",
                          "g_tot": f"+{s_gtot:.2f}*eps (kappa^2*g_tot(0.05) = "
                                   f"{k2g_05:.4f} vs stake 1.45833)"},
    "scan_table": scan_rows,
    "anchor": {"old_2.371_pull_sigma": pull_old,
               "repaired_c005_pull_sigma": pull_new,
               "verdict": "FORFEITED as physical-branch corroboration; "
                          "survives as restricted-branch code-validation record"},
    "bands_4.10prime": {"c_phys": [c_satf, c_band_hi],
                        "kappa_phys": [float(k_sat), float(k_sat)],
                        "note": "one-sided band flips to the + side of c; "
                                "kappa/depth/omega_th have NO eps-band"},
    "survive_shift_undefined": [
        {"item": "c0 = 128sqrt42/105pi = 2.5147 (intercept)",
         "verdict": "STRUCK/superseded -> c_sat = 3.201125 [T3.1]"},
        {"item": "c_g = 0.42 +/- 0.04 (amplitude)",
         "verdict": "UNDEFINED on the repaired branch (dies with the law); "
                    "replaced by s ~ +1.15 (two-point, no band yet)"},
        {"item": "exponent 2/3 + lambda*^2 ~ eps^{2/3} geometric narrative",
         "verdict": "UNDEFINED (fixed-oblateness minimiser -> analytic/linear)"},
        {"item": "sign (falling)", "verdict": "FLIPS (rising)"},
        {"item": "anchor 2.371 vs <r1> 2.37+/-0.09 (0.016 sigma)",
         "verdict": "arithmetic SURVIVES as restricted-branch record; "
                    "FORFEITED as corroboration (repaired: 11.28 sigma, correctly)"},
        {"item": "three-way loop dc/c0 = 0.057 -> g* = 1.329 vs 1.31+/-0.04",
         "verdict": "SURVIVES as restricted-branch record (B.6); UNDEFINED on "
                    "the saturated branch (no eps-running lambda*)"},
        {"item": "kappa(eps), depth(eps), omega_th(eps) running + one-sided bands",
         "verdict": "SHIFT to exact eps-frozen constants (1/sqrt2, 29.29%, "
                    "sqrt2*omega0); bands collapse to zero; only c keeps a "
                    "(+)side band"},
        {"item": "lambda*(eps) = 0.42 eps^{1/3}; lambda*(eps_q) = 0.246",
         "verdict": "SUPERSEDED (V15.7 row, i2); shape re-scan obligation OPEN"},
        {"item": "<r1> 'benchmark-at-eps=0.05' designation",
         "verdict": "SURVIVES with restricted-branch/transit label"},
        {"item": "eps range [1e-3, 0.1] as scan domain",
         "verdict": "SURVIVES as domain; saturated branch measured at 2 points "
                    "only — full re-scan OPEN (D5)"},
        {"item": "electron oblateness sentence lambda_e ~ eps^{1/3} ~ 0.01-0.1",
         "verdict": "UNDEFINED/superseded — still printed in v3.0-ext IV.J' (D2)"},
        {"item": "eps-suppression third leg of point-likeness (<= 6e-4)",
         "verdict": "SHIFTS to <= 2.6e-3 (x4.3 weaker, still small)"},
    ],
}
flush()

# ===========================================================================
# B. (4.9') — the censorship ceiling and the Q-6' margins, done properly
# ===========================================================================
print("== B. (4.9') ceiling + Q-6' censorship ==")
c_g, dc_g = 0.42, 0.04
ceil_old = (D_lock / c_g) ** 1.5
ceil_old_band = [(D_lock / (c_g + dc_g)) ** 1.5, (D_lock / (c_g - dc_g)) ** 1.5]
gate("B.1 old ceiling (3e-3/0.42)^{3/2} = 6.04e-4 (printed '~6e-4')",
     close(ceil_old, 6.0368e-4, 1e-3),
     f"= {ceil_old:.4e}; c_g+/-0.04 band [{ceil_old_band[0]:.2e}, "
     f"{ceil_old_band[1]:.2e}]")

ceil_new = D_lock / s_exact
ceil_new_eng = D_lock / s_engine
lam_f, dlam_f = 3.1e-3, 0.7e-3           # lambda_field = (3.1+/-0.7)e-3 omega0
ceil_new_Dband = [2.4e-3 / s_exact, 3.8e-3 / s_exact]
gate("B.2 repaired ceiling Delta_lock/s = 2.60e-3 (printed '~2.6e-3')",
     close(ceil_new, 2.6047e-3, 1e-3),
     f"= {ceil_new:.4e} (engine-base {ceil_new_eng:.4e}); Delta_lock band "
     f"[2.4,3.8]e-3 -> ceiling [{ceil_new_Dband[0]:.2e}, {ceil_new_Dband[1]:.2e}]")

thin = ceil_new / ceil_old
gate("B.3 thinning provenance: ceil_new/ceil_old = x4.31 (the printed 'x4.3')",
     close(thin, 4.315, 2e-3), f"= x{thin:.3f} = (D/s)/((D/c_g)^{{3/2}})")

# B.4 exponent sensitivity: same two saturated points read at exponent p
def ceiling_p(p, rise=rise_exact):
    amp = rise / eps5 ** p
    return (D_lock / amp) ** (1 / p), amp

ceil_23, amp_23 = ceiling_p(2 / 3)
ceil_1, amp_1 = ceiling_p(1.0)
mono = all(ceiling_p(p)[0] <= ceiling_p(q)[0] + 1e-15
           for p, q in zip(np.linspace(2 / 3, 1, 8)[:-1],
                           np.linspace(2 / 3, 1, 8)[1:]))
gate("B.4 exponent sensitivity: ceiling(p) monotone up on [2/3,1]; "
     "p=2/3 -> 5.9e-4, p=1 -> 2.60e-3",
     mono and close(ceil_23, 5.94e-4, 2e-2) and close(ceil_1, ceil_new, 1e-9),
     f"p=2/3: amp = {amp_23:.4f} (numerology flag: ~ old c_g 0.42+/-0.04), "
     f"ceiling = {ceil_23:.2e}; the linear reading is the LOOSEST ceiling (D3)")

# B.5 eps_q per reading
mSk = 1.7
fq_lo, fq_c, fq_hi = 0.14, 0.17, 0.20
eq_lo, eq_c, eq_hi = (fq_lo / mSk) ** 2, (fq_c / mSk) ** 2, (fq_hi / mSk) ** 2
gate("B.5 eps_q = (0.14-0.20/1.7)^2 = [6.78e-3, 1.384e-2], central 1.00e-2 "
     "(printed '0.7-1.4e-2', '~1e-2')",
     close(eq_lo, 6.782e-3, 1e-3) and close(eq_hi, 1.3841e-2, 1e-3)
     and close(eq_c, 1.0e-2, 1e-6),
     f"[{eq_lo:.4e}, {eq_hi:.4e}], central {eq_c:.4e}")

marg_old = [eq_lo / ceil_old, eq_c / ceil_old, eq_hi / ceil_old]
marg_new = [eq_lo / ceil_new, eq_c / ceil_new, eq_hi / ceil_new]
gate("B.6 OLD margins [11.2, 22.9], central 16.6 — 'twenty-fold' was fair for "
     "v2.0.1; STALE in v3.0-ext Q-6' (D1)",
     close(marg_old[0], 11.234, 1e-3) and close(marg_old[2], 22.927, 1e-3)
     and close(marg_old[1], 16.565, 1e-3),
     f"old = [{marg_old[0]:.2f}, {marg_old[1]:.2f}, {marg_old[2]:.2f}]")
gate("B.7 REPAIRED margins [2.60, 5.31], central 3.84 — Q-6' should read "
     "'~x4 above' (x2.6-5.3), not 'twenty-fold'",
     close(marg_new[0], 2.6037, 1e-3) and close(marg_new[2], 5.3137, 1e-3)
     and close(marg_new[1], 3.8392, 1e-3),
     f"new = [{marg_new[0]:.2f}, {marg_new[1]:.2f}, {marg_new[2]:.2f}]")

# B.8 the dichotomy at every reading (quarks censored <=> eps_q > ceiling;
#     leptons free <=> admission window below ceiling non-empty)
readings = [
    ("R1 central: eps_q central / ceiling central", eq_c / ceil_new),
    ("R2 band-edge low eps_q / ceiling central", eq_lo / ceil_new),
    ("R3 printed 'eps_dress ~ 1e-2' / ceiling central", 1e-2 / ceil_new),
    ("R4 low eps_q / Delta_lock-high ceiling (3.8e-3)", eq_lo / (3.8e-3 / s_exact)),
    ("R5 low eps_q / engine-base-slope ceiling", eq_lo / ceil_new_eng),
    ("R6 low eps_q / exponent-2/3 ceiling (tightest)", eq_lo / ceil_23),
    ("R7 WORST JOINT: low eps_q / (Delta hi & engine slope)",
     eq_lo / (3.8e-3 / s_engine)),
]
min_marg = min(m for _, m in readings)
for nm, m in readings:
    print(f"      {nm}: margin x{m:.2f}")
gate("B.8 dichotomy verdict: eps_q > ceiling at EVERY reading (all six quarks "
     "censored); lepton window [1e-5, ceiling] non-empty (all leptons free)",
     min_marg > 1.0 and 1.5e-6 < 1.0e-5 < ceil_new,
     f"min margin x{min_marg:.2f} (worst joint reading); floor 1.5e-6 < "
     f"operative bottom 1e-5 < ceiling {ceil_new:.2e}; termination threshold "
     "8.4e-6 < 1e-5 unchanged — SURVIVES AT EVERY READING")

# B.9 breach conditions (what further ceiling shift kills the dichotomy)
breach_factor_central = eq_lo / ceil_new
breach_D = s_exact * eq_lo                    # Delta_lock needed to breach
breach_D_sigma = (breach_D - lam_f) / dlam_f  # in lambda_field sigmas
breach_s = D_lock / eq_lo                     # slope needed to breach
breach_fq = np.sqrt(ceil_new) * mSk           # f_q low edge needed to breach
gate("B.9 breach conditions: ceiling must rise a FURTHER x2.60 (central; "
     "x2.03 at the worst joint reading)",
     close(breach_factor_central, 2.604, 1e-3) and close(min_marg, 2.028, 2e-3),
     f"equivalently Delta_lock >= {breach_D:.2e} (lambda_field "
     f"{breach_D_sigma:+.1f} sigma above its measured 3.1+/-0.7e-3), or slope "
     f"s <= {breach_s:.3f} (vs +1.15), or sigma-inversion low edge f_q <= "
     f"{breach_fq:.4f} GeV (vs 0.14)")

# B.10 window survival
gate("B.10 admission window: ceiling 2.60e-3 vs II.H window top 3e-3 — "
     "window survives numerically, top now ceiling-set (86.8% of 3e-3)",
     0.8 < ceil_new / 3e-3 < 1.0, f"ratio = {ceil_new/3e-3:.3f}")

R["sB_49prime"] = {
    "Delta_lock": {"value": D_lock, "provenance":
                   "IV.H.3 lambda_field/omega ~ 3e-3 (lambda_field = "
                   "(3.1+/-0.7)e-3 omega0)"},
    "ceiling_old": {"form": "(Delta/c_g)^{3/2}", "value": ceil_old,
                    "c_g_band": ceil_old_band},
    "ceiling_new": {"form": "Delta/s (linear sign-flipped law)",
                    "value": ceil_new, "engine_base": ceil_new_eng,
                    "Delta_band": ceil_new_Dband,
                    "exponent_sensitivity": {"p=2/3": ceil_23, "p=1": ceil_1,
                                             "amp_at_2/3": amp_23}},
    "thinning_factor": thin,
    "eps_q": {"band": [eq_lo, eq_hi], "central": eq_c,
              "provenance": "(f_q/m_Sk)^2, f_q in [0.14, 0.20] GeV "
                            "(sigma-inversion, Q-5), m_Sk = 1.7 GeV (blind)"},
    "margins": {"old": marg_old, "new": marg_new,
                "readings": [{"reading": nm, "margin": m} for nm, m in readings],
                "min_margin": min_marg},
    "dichotomy": {"verdict": "SURVIVES AT EVERY READING (min x2.03)",
                  "quarks": "all six censored (eps_q > ceiling everywhere)",
                  "leptons": "all free (window [1e-5, 2.6e-3] non-empty; "
                             "floor and termination threshold untouched)"},
    "breach": {"further_ceiling_factor_central": breach_factor_central,
               "further_ceiling_factor_worst": min_marg,
               "Delta_lock_needed": breach_D,
               "lambda_field_sigma": breach_D_sigma,
               "slope_needed": breach_s, "fq_low_edge_needed_GeV": breach_fq},
    "defect_D1": "v3.0-ext Q-6' retains 'twenty-fold above the entrainment "
                 "ceiling' verbatim while its own (4.9') prints 2.6e-3 — "
                 "internally consistent figure is x2.6-5.3 (central ~x3.8)",
}
flush()

# ===========================================================================
# C. Cross-check vs i2_propagate.py / i2_results.json
# ===========================================================================
print("== C. cross-check vs i2 ==")
with open(P(HERE, "i2_results.json")) as fh:
    i2 = json.load(fh)
er = i2["exposed_register"]
pairs = [
    ("ceiling_old", er["BU1_ceiling"]["old"], ceil_old, 1e-9),
    ("ceiling_new", er["BU1_ceiling"]["new"], ceil_new, 1e-4),
    ("margin_new_lo", er["quark_censorship_margin"]["new"][0], marg_new[0], 1e-4),
    ("margin_new_hi", er["quark_censorship_margin"]["new"][1], marg_new[2], 1e-4),
    ("margin_old_lo", er["quark_censorship_margin"]["old"][0], marg_old[0], 1e-9),
    ("margin_old_hi", er["quark_censorship_margin"]["old"][1], marg_old[2], 1e-9),
    ("anchor_old_pull", er["eps_law"]["anchor_old_pull_sigma"], pull_old, 1e-6),
    ("repaired_pull_r1", er["eps_law"]["repaired_vs_r1_sigma"], pull_new, 1e-4),
    ("operative_top_new", er["operative_eps_range"]["new"][1], ceil_new, 1e-4),
    ("f_window_top_MeV", er["f_window_MeV"]["new"][1],
     0.69 * np.sqrt(ceil_new) * 1700.0, 1e-4),
]
disc = []
for nm, iv, jv, tol in pairs:
    ok = close(iv, jv, tol)
    if not ok:
        disc.append({"item": nm, "i2": iv, "j2": jv})
    gate(f"C i2-overlap {nm}", ok, f"i2 = {iv:.6e}, j2 = {jv:.6e}")
slope_i2 = float(er["eps_law"]["new"].split("+")[1].strip().split(" ")[0])
gate("C i2 slope string '+1.15 eps' vs j2 s_exact",
     abs(slope_i2 - round(s_exact, 2)) < 5e-3,
     f"i2 prints +{slope_i2}, j2 = {s_exact:.4f} (exact-base reading)")
note_f = ("i2/j2 f-window bottom 3.71 MeV vs corpus F-A15-3's printed "
          "'f in [4, 29] MeV' — corpus rounds 3.7->4; substance identical (D6)")
print("   note:", note_f)
R["sC_crosscheck"] = {"pairs": [{"item": nm, "i2": iv, "j2": jv,
                                 "agree": close(iv, jv, tol)}
                                for nm, iv, jv, tol in pairs],
                      "discrepancies": disc, "note_f_window": note_f}
flush()

# ===========================================================================
# D. Defect register (printed honestly) + gate table
# ===========================================================================
R["sD_defects"] = [
    {"id": "D1", "against": "corpus v3.0-ext Q-6'",
     "text": "retains 'twenty-fold above the entrainment ceiling' while its "
             "own (4.9') prints ceiling 2.6e-3; consistent figure x2.6-5.3 "
             "(central ~x3.8). Stale sentence — the (4.9') re-derivation "
             "obligation includes this fix."},
    {"id": "D2", "against": "corpus v3.0-ext IV.J' consequences",
     "text": "retains 'lambda_e ~ eps_e^{1/3} ~ 0.01-0.1' — the superseded "
             "restricted-branch shape law inside the repaired section."},
    {"id": "D3", "against": "the repaired law itself (and this memo)",
     "text": "the '+1.15 eps' LINEAR class rests on TWO saturated points "
             "(t->0, eps=0.05); the exponent is unmeasured. Within-model "
             "justification: fixed-oblateness minimiser -> analytic leading "
             "order. Sensitivity: ceiling(p) rises monotonically to the "
             "linear reading's 2.60e-3; a 2/3 reading gives 5.9e-4 with "
             "amplitude 0.424 ~ old c_g (numerology flag, not charged). The "
             "dichotomy is exponent-robust; the x4.3 loosening (and the 60 "
             "MeV f-top / 7.8e-10 Majoron edge downstream) is linear-specific."},
    {"id": "D4", "against": "axi4 engine base",
     "text": "ring_t0 c_paper = 3.203459 vs exact 3.201125 (+7.3e-4), x5 its "
             "own N-ladder estimate 1.4e-4; propagates a 1.4% slope-reading "
             "spread (1.152 exact-base vs 1.136 engine-base) and 2.60 vs "
             "2.64e-3 in the ceiling. Not material to any verdict."},
    {"id": "D5", "against": "scope of this memo",
     "text": "J2 discharges the (4.8')/(4.9') ARITHMETIC; the multi-eps "
             "engine re-scan of the saturated branch (and the lambda* shape "
             "re-scan) remain OPEN — rows at eps != 0.05 in the scan table "
             "are interpolation/extrapolation grade."},
    {"id": "D6", "against": "i2/j2 vs corpus print",
     "text": "f-window bottom: computed 3.71 MeV vs F-A15-3's printed "
             "'[4, 29] MeV' (corpus rounding; substance identical)."},
]

n_pass = sum(1 for g in GATES if g["pass"])
R["gate_table"] = GATES
R["summary"] = {
    "gates": f"{n_pass}/{len(GATES)} PASS",
    "48prime": f"c_sat(eps) = 3.201125*(1 + {s_exact:.3f}*eps) leading order "
               "(sign flipped, exponent 2/3 -> 1); kappa/depth/omega_th/Erot "
               "eps-FROZEN exactly; g_tot runs +0.57*eps; anchor forfeited "
               "(11.28 sigma from <r1>, correctly)",
    "49prime": f"ceiling = Delta_lock/s = {ceil_new:.3e} (x{thin:.2f} looser); "
               f"Q-6' margins [2.60, 5.31] central 3.84; dichotomy SURVIVES "
               f"at every reading (min x{min_marg:.2f}); breach needs a "
               f"further x2.60 (central) ceiling rise",
    "crosscheck": "all i2 overlaps agree" if not disc else f"DISCREPANCIES: {disc}",
}
flush()

print(f"\n== J2 gates: {n_pass}/{len(GATES)} PASS; j2_results.json written ==")
if n_pass != len(GATES):
    raise SystemExit(1)
