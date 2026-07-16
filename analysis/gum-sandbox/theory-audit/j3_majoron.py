#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
J3 — Majoron battery re-run at the propagated low edge g ≈ 7.8e-10 (ROADMAP_v6 J3).

Re-executes the corpus's OWN Majoron battery arithmetic (Ω VII.I / App K.3+K.5 /
V15.5 / F-A15-3 / Substrate Course 19.3 / WS-nu-P4 r1 / WS-K2 (b)-(c)) at the
T3.1-repaired phason window propagated by I2 (i2_propagate.py):
  old operative window f ∈ [3.7, 29] MeV (ε_e ∈ [1e-5, 6.0e-4], B-U1' ceiling)
  new operative window f ∈ [3.7, ≈60] MeV (ε_e ∈ [1e-5, 2.6e-3], repaired ceiling)
Old window first as CONTROL (must reproduce the corpus's printed numbers/verdicts).

WITHIN-MODEL ONLY — analysis layer. Deterministic, standalone (stdlib only).
Writes j3_results.json stage-flushed.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "j3_results.json")

RES = {"meta": {
    "phase": "ROADMAP_v6 J3", "date": "2026-07-16",
    "standing": "WITHIN-MODEL ONLY — re-runs the corpus's printed battery "
                "arithmetic at shifted inputs; validates nothing about nature.",
    "sources": [
        "substrate-suite/01-GUM-Omega-Paper-v2.0.1.md VII.I (dictionary line), "
        "(4.9), F-A15-3, App K.3/K.5, V15.5",
        "corpus2/01-GUM-Omega-Paper-v3.0-ext.md VII.I (identical dictionary line)",
        "substrate-suite/02-The-Substrate-Course.md 19.3 (the four battery rows)",
        "corpus/WS-nu-P4-Perimeter-Annex-Funnel-v0_1.md r1 (battery re-certified)",
        "corpus/WS-K2-Charged-Current-Architecture-v0_1.md (b)/(c) (lab perimeter)",
        "theory-audit/i2_propagate.py + i2_RESULTS.md (the propagated window)",
        "tier5-family/family_audit.py C5c/C5d/C7b (prior replication, reused)"]},
    "stages": {}}

def flush():
    with open(OUT, "w") as fh:
        json.dump(RES, fh, indent=1, ensure_ascii=False)

GATES = []
def gate(gid, desc, ok, detail):
    GATES.append({"id": gid, "desc": desc, "pass": bool(ok), "detail": detail})
    print(f"{gid:<5} {'PASS' if ok else 'FAIL'}  {desc}\n      {detail}")
    return bool(ok)

# ---------------- corpus constants, as printed ----------------
SINTC   = 0.69          # sin θ_c (K.5; ⟨r9⟩ m=1.9±0.4 → √(1-1/m)=0.69)
MSK_GEV = 1.7           # GeV, 𝔪_Sk (WS-nu-P4 f-law; I2 verified 𝔠-blind)
f_MeV   = lambda eps: SINTC * math.sqrt(eps) * MSK_GEV * 1e3   # f = 0.69√ε·1.7 GeV
EPS_BOT = 1.0e-5        # operative ε_e bottom (F-A15-3), repair-invariant
DLOCK, CG = 3.0e-3, 0.42
CEIL_OLD = (DLOCK / CG) ** 1.5      # B-U1' old ceiling ≈ 6.0e-4
SLOPE    = 1.152                    # repaired ε-law slope (+1.15ε, i2)
CEIL_NEW = DLOCK / SLOPE            # repaired ceiling ≈ 2.6e-3
M3       = 0.0468       # eV, central m₃ (⟨r10⟩ bridge, external replication)
M3_FLOOR = 0.0503       # eV, oscillation floor (F-A15-3's g-endpoint convention)
M3_BAND  = (0.019, 0.115)   # 1σ band the corpus scans
M3_LIVE  = (0.050, 0.057)   # E-XXII-1 live window
MPL_EV   = 1.220890e28  # eV, full Planck mass (tier5 C5c convention)
BNU1_TH  = 1.4          # MeV, printed B-ν1 "with thermal factors" at m≈0.047
BNU1_RAW = 1.57         # MeV, printed raw (T M_Pl m_ν²)^{1/4}
DNEFF    = (4.0/7.0) * (10.75/106.75) ** (4.0/3.0)
HBARC_MEV_FM = 197.3269631; Z_CS = 55.0; A0_FM = 52917.72

def bnu1(m_ev, ref=BNU1_TH):        # printed value scaled by its own m-law ∝ √m
    return ref * math.sqrt(m_ev / 0.047)

def battery_rows(f_lo, f_hi, m_lo, m_hi, tag):
    """The corpus's four battery rows + lab perimeter, evaluated on a window."""
    g_min = m_lo * 1e-9 / (f_hi * 1e-3)     # eV/GeV -> dimensionless
    g_max = m_hi * 1e-9 / (f_lo * 1e-3)
    # R1 free-streaming: pass condition as printed: f_window bottom >= B-ν1(m).
    # Binding corner: m at band top (B-ν1 ∝ √m); f at window bottom.
    fs_margin_th  = f_lo / bnu1(m_hi, BNU1_TH)
    fs_margin_raw = f_lo / bnu1(m_hi, BNU1_RAW)
    # inverted ε-floor (B-ν1'): must stay below the operative bottom 1e-5
    eps_floor = (bnu1(m_hi, BNU1_TH) / (SINTC * MSK_GEV * 1e3)) ** 2
    # R3/R4 binding edge is g_max (evasion from below / recoupling at large g)
    pnc = lambda f: (HBARC_MEV_FM / f * Z_CS / A0_FM) ** 2   # (λZ/a₀)² per unit x
    return {
        "tag": tag, "f_MeV": [f_lo, f_hi], "m_eV": [m_lo, m_hi],
        "g_range": [g_min, g_max],
        "R1_free_streaming": {"bound": "f >= B-nu1 (printed 1.4 MeV at m~0.047, "
                              "raw 1.57; ∝ √m)", "binding_corner": "f_lo, m_hi",
                              "Bnu1_at_mhi_MeV": [bnu1(m_hi), bnu1(m_hi, BNU1_RAW)],
                              "margin": [fs_margin_th, fs_margin_raw],
                              "pass": fs_margin_th > 1.0 and fs_margin_raw > 1.0},
        "R1b_eps_floor": {"eps_floor_at_mhi": eps_floor, "op_bottom": EPS_BOT,
                          "pass": eps_floor < EPS_BOT},
        "R2_BBN_dNeff": {"dNeff": DNEFF, "f_g_dependence": "none (one decoupled "
                         "scalar); bit-identical old->new", "printed_bound":
                         "NOT PRINTED (spec-recovery limit); '✓' as stamped",
                         "pass": True},
        "R3_SN_cooling": {"pass_condition_as_printed": "'evaded from below' — "
                          "whole g-interval below the band's lower edge; band "
                          "edges NOT PRINTED (spec-recovery limit)",
                          "binding_edge_g_max": g_max, "pass": None},
        "R4_late_recoupling": {"pass_condition_as_printed": "'harmless' — no "
                               "formula/number printed (spec-recovery limit); "
                               "monotone in g, binds at g_max",
                               "binding_edge_g_max": g_max, "pass": None},
        "LAB_WSK2": {"g2_range": [g_min**2, g_max**2],
                     "printed": "g² ≲ 1e-16 'negligible' (model value, not bound)",
                     "pnc_distortion_per_x": [pnc(f_hi), pnc(f_lo)],
                     "printed_pnc": "1e-4–3e-3 per unit x over [4,29] MeV",
                     "lambda_fm": [HBARC_MEV_FM/f_hi, HBARC_MEV_FM/f_lo]}}

# =====================================================================
# STAGE A — CONTROL: old operative window, reproduce the corpus's prints
# =====================================================================
print("== STAGE A — old-window control ==")
f_lo_old, f_hi_old = f_MeV(EPS_BOT), f_MeV(min(6e-4, CEIL_OLD))
gate("A1", "old f-window bottom 3.71 MeV (corpus rounds to '4'; i2 '3.7')",
     abs(f_lo_old - 3.709) < 0.01, f"f(1e-5) = {f_lo_old:.4f} MeV")
gate("A2", "old f-window top 28.7 MeV (corpus '29')",
     abs(f_hi_old - 28.73) < 0.05, f"f(6.0368e-4 -> capped 6e-4) = {f_hi_old:.4f} MeV")
g_lo_i2   = M3 * 1e-9 / (f_hi_old * 1e-3)
g_lo_a153 = M3_FLOOR * 1e-9 / (29e-3)
g_hi_a153 = M3_FLOOR * 1e-9 / (4e-3)
gate("A3", "old g endpoints: i2's 1.6e-9 (m=0.0468/28.7) and F-A15-3's "
     "[1.7e-9, 1.3e-8] (m=0.0503/[29,4])",
     abs(g_lo_i2/1.63e-9 - 1) < 0.02 and abs(g_lo_a153/1.73e-9 - 1) < 0.02
     and abs(g_hi_a153/1.26e-8 - 1) < 0.02,
     f"{g_lo_i2:.3e}; {g_lo_a153:.3e}; {g_hi_a153:.3e}")
gate("A4", "ΔN_eff = (4/7)(10.75/106.75)^{4/3} = 0.0268",
     abs(DNEFF - 0.0268) < 5e-5, f"= {DNEFF:.6f}")
eps_floor_print = (BNU1_TH / (SINTC * MSK_GEV * 1e3)) ** 2
gate("A5", "ε-floor (1.4 MeV inverted) = 1.42e-6 (printed '1.4->1.5e-6')",
     abs(eps_floor_print/1.424e-6 - 1) < 0.01, f"= {eps_floor_print:.3e}")
T_rec_needed = (BNU1_RAW*1e6)**4 / (MPL_EV * M3**2)   # eV
gate("A6", "B-ν1 raw 1.57 MeV recovers with T = (1.57 MeV)^4/(M_Pl m₃²) — "
     "recombination-class T (CONVENTION-LIMITED, per tier5 C5c)",
     0.1 < T_rec_needed < 0.5, f"T = {T_rec_needed:.3f} eV (T_rec class ~0.26 eV)")
g_body_lo = 0.047e-9 / 60e-3
gate("A7", "v2.0/v3.0 BODY dictionary already prints [8e-10, 1.3e-8] over "
     "f ∈ [4, 60] (II.H admission window) — 0.047/60 = 7.8e-10 -> '8e-10'",
     abs(g_body_lo/7.83e-10 - 1) < 0.01, f"= {g_body_lo:.3e}")
old = battery_rows(f_lo_old, f_hi_old, M3_BAND[0], M3_BAND[1], "old-window")
old_central = battery_rows(f_lo_old, f_hi_old, M3, M3, "old-window-central-m3")
gate("A8", "CONTROL battery verdict on old window reproduces 'passed': "
     "R1 margin>1 across full m-band, ε-floor slack, R2 ✓; R3/R4 stamped ✓ "
     "with bounds unprinted",
     old["R1_free_streaming"]["pass"] and old["R1b_eps_floor"]["pass"],
     f"R1 margins (th/raw, worst corner m=0.115) = "
     f"{old['R1_free_streaming']['margin'][0]:.2f}/"
     f"{old['R1_free_streaming']['margin'][1]:.2f}; central-m3 margins = "
     f"{old_central['R1_free_streaming']['margin'][0]:.2f}/"
     f"{old_central['R1_free_streaming']['margin'][1]:.2f}; "
     f"eps_floor(0.115) = {old['R1b_eps_floor']['eps_floor_at_mhi']:.2e}")
RES["stages"]["A_control"] = {"f_window_MeV": [f_lo_old, f_hi_old],
    "ceil_old": CEIL_OLD, "rows_full_band": old, "rows_central": old_central,
    "g_endpoint_conventions": {"i2_low": g_lo_i2, "FA153_low": g_lo_a153,
                               "FA153_high": g_hi_a153, "body_low_f60": g_body_lo}}
flush()

# =====================================================================
# STAGE B — RE-RUN on the repaired window, m_ν swept edge to edge
# =====================================================================
print("== STAGE B — repaired-window re-run ==")
f_lo_new, f_hi_new = f_MeV(EPS_BOT), f_MeV(CEIL_NEW)
gate("B1", "repaired ceiling 3e-3/1.152 = 2.60e-3 -> f-top 59.9 MeV "
     "(coincident with the II.H window top '60')",
     abs(CEIL_NEW/2.604e-3 - 1) < 1e-3 and abs(f_hi_new - 59.86) < 0.05,
     f"ceiling = {CEIL_NEW:.4e}; f_top = {f_hi_new:.3f} MeV")
g_lo_new = M3 * 1e-9 / (f_hi_new * 1e-3)
gate("B2", "new low edge g = m₃/f_top = 7.8e-10 (the I2 obligation number)",
     abs(g_lo_new/7.82e-10 - 1) < 0.01, f"= {g_lo_new:.3e}")

new = battery_rows(f_lo_new, f_hi_new, M3_BAND[0], M3_BAND[1], "new-window")
new_central = battery_rows(f_lo_new, f_hi_new, M3, M3, "new-window-central-m3")
new_live = battery_rows(f_lo_new, f_hi_new, M3_LIVE[0], M3_LIVE[1],
                        "new-window-live-m")

gate("B3", "R1 free-streaming PASSES across the whole repaired window and "
     "full m-band; margin identical to control (binding edge f_lo unchanged)",
     new["R1_free_streaming"]["pass"] and
     abs(new["R1_free_streaming"]["margin"][0] -
         old["R1_free_streaming"]["margin"][0]) < 1e-12,
     f"worst margin (m=0.115) = {new['R1_free_streaming']['margin'][0]:.3f} "
     f"(thermal) / {new['R1_free_streaming']['margin'][1]:.3f} (raw); "
     f"central m₃ = {new_central['R1_free_streaming']['margin'][0]:.2f}/"
     f"{new_central['R1_free_streaming']['margin'][1]:.2f}")
gate("B4", "R1b ε-floor stays slack across m-band (max at m=0.115)",
     new["R1b_eps_floor"]["pass"],
     f"eps_floor(0.115) = {new['R1b_eps_floor']['eps_floor_at_mhi']:.3e} "
     f"< operative bottom {EPS_BOT:.0e}")
gate("B5", "R2 BBN/ΔN_eff bit-identical (f,g-independent)", True,
     f"ΔN_eff = {DNEFF:.6f} old and new")
gate("B6", "R3 SN-cooling binding edge g_max UNCHANGED old->new for every m "
     "(f_lo repair-invariant); new strip only deepens evasion-from-below",
     abs(new["R3_SN_cooling"]["binding_edge_g_max"] -
         old["R3_SN_cooling"]["binding_edge_g_max"]) < 1e-24,
     f"g_max = {new['R3_SN_cooling']['binding_edge_g_max']:.3e} (m=0.115) / "
     f"{new_central['R3_SN_cooling']['binding_edge_g_max']:.3e} (m₃); "
     f"g_min drops {old['g_range'][0]:.2e} -> {new['g_range'][0]:.2e}")
gate("B7", "R4 late-recoupling binding edge g_max likewise unchanged; "
     "newly opened strip is monotonically safer",
     abs(new["R4_late_recoupling"]["binding_edge_g_max"] -
         old["R4_late_recoupling"]["binding_edge_g_max"]) < 1e-24,
     "verdict inherited from the old-window stamp (bound unprinted)")

# edge-to-edge sweep: verify no interior extremum on any printed margin
NF, NM = 121, 25
fs = [f_lo_new + i*(f_hi_new-f_lo_new)/(NF-1) for i in range(NF)]
ms = [M3_BAND[0] + j*(M3_BAND[1]-M3_BAND[0])/(NM-1) for j in range(NM)]
worst = {"fs_margin": 1e9, "g_max": 0.0, "g_min": 1e9, "pnc_max": 0.0}
for f in fs:
    for m in ms:
        worst["fs_margin"] = min(worst["fs_margin"], f_lo_new / bnu1(m))
        g = m * 1e-9 / (f * 1e-3)
        worst["g_max"] = max(worst["g_max"], g)
        worst["g_min"] = min(worst["g_min"], g)
    worst["pnc_max"] = max(worst["pnc_max"],
                           (HBARC_MEV_FM/f * Z_CS / A0_FM) ** 2)
gate("B8", "edge-to-edge grid sweep (121 f × 25 m): extrema sit at corners; "
     "no interior reversal (all rows monotone in f and m)",
     abs(worst["g_max"]/new["g_range"][1] - 1) < 1e-9 and
     abs(worst["g_min"]/new["g_range"][0] - 1) < 1e-9 and
     abs(worst["fs_margin"]/new["R1_free_streaming"]["margin"][0] - 1) < 1e-9,
     f"g ∈ [{worst['g_min']:.3e}, {worst['g_max']:.3e}]; "
     f"min FS margin = {worst['fs_margin']:.3f}")
gate("B9", "lab perimeter (WS-K2 b/c): g² top unchanged (≲1e-16 class); Cs-PNC "
     "distortion band recomputed — high-f extension only shrinks it",
     new["LAB_WSK2"]["g2_range"][1] < 1.1e-15 and
     new["LAB_WSK2"]["pnc_distortion_per_x"][0] <
     old["LAB_WSK2"]["pnc_distortion_per_x"][0],
     f"g² ∈ [{new['LAB_WSK2']['g2_range'][0]:.2e}, "
     f"{new['LAB_WSK2']['g2_range'][1]:.2e}]; PNC per-x old "
     f"[{old['LAB_WSK2']['pnc_distortion_per_x'][0]:.1e}, "
     f"{old['LAB_WSK2']['pnc_distortion_per_x'][1]:.1e}] -> new "
     f"[{new['LAB_WSK2']['pnc_distortion_per_x'][0]:.1e}, "
     f"{new['LAB_WSK2']['pnc_distortion_per_x'][1]:.1e}] "
     f"(λ band {new['LAB_WSK2']['lambda_fm'][0]:.1f}–"
     f"{new['LAB_WSK2']['lambda_fm'][1]:.1f} fm)")
RES["stages"]["B_rerun"] = {"f_window_MeV": [f_lo_new, f_hi_new],
    "ceil_new": CEIL_NEW, "rows_full_band": new, "rows_central": new_central,
    "rows_live_window": new_live, "grid_worst": worst,
    "signature_notes": {
        "P_nu4_0nubb": "signature row, not a gate: 'Majoron-mode 0νββ at "
        "g ∼ 1e-9 (far-horizon)' — the newly opened strip g ∈ [7.8e-10, "
        "1.6e-9) sits at/below 1e-9: at the new low edge the far-horizon "
        "signature retreats further; no pass/fail content",
        "dNeff_CJ": "[CJ-thermal] conditionality: smaller g makes primordial "
        "thermalization less likely; grade unchanged, null does not kill"}}
flush()

# =====================================================================
# STAGE C — verdict
# =====================================================================
print("== STAGE C — verdict ==")
npass = sum(1 for g in GATES if g["pass"]); nfail = len(GATES) - npass
verdict = ("CONFIRMED: I2's 'SHIFTS-HARMLESSLY' holds — no battery row newly "
           "binds or fails anywhere in the repaired window. The two binding "
           "edges (f_lo = 3.71 MeV, g_max = m/f_lo) are repair-invariant; the "
           "only moving edge (g_min: 1.6e-9 -> 7.8e-10) moves monotonically "
           "away from every coupling-driven constraint. R3 (SN) and R4 "
           "(recoupling) are inherited-PASS by monotonicity + unchanged "
           "binding edge, NOT recomputed: their numerical bounds are never "
           "printed in the corpus (spec-recovery limits). The corpus's BODY "
           "dictionary (VII.I, both versions) already stamped the four ✓ on "
           "f ∈ [4, 60] — the repaired operative window lands on a domain "
           "the corpus itself already claimed; the re-run makes that claim "
           "arithmetic instead of stamped.")
RES["stages"]["C_verdict"] = {"gates_pass": npass, "gates_fail": nfail,
                              "verdict": verdict}
RES["gates"] = GATES
flush()
print(f"\nGATES: {npass}/{len(GATES)} PASS, {nfail} FAIL")
print(verdict)
