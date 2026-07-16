#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K2 — the eps < 0.005 probe (ROADMAP_v7 Phase K; extends J4 = J2 defect D3).

QUESTION. J4 measured the saturated branch's small-eps law as LINEAR
(p = 0.997, profile CI [0.996, 0.998], envelope [0.996, 1.012]) on
eps in [0.005, 0.2] and printed defect J4-D3: eps < 0.005 was not probed.
Does the linear law hold at eps in {0.001, 0.002} — the smallest honestly
solvable eps — or does the exponent drift?

METHOD (J4 code path REUSED VERBATIM where marked):
  (a) new saturated solves at eps in {0.001, 0.002}: i1's dial +
      saturated_at unmodified, i1 full budgets (N_OPT/N_POL/N_FIN/N_COARSE
      = 192/320/480/320, B_TRI/B_COLD/B_POL = 220/1500/260), two grid
      resolutions per point (N_FIN=480 / N_COARSE=320) + the dial's own
      N2x/rmax+ ladder.  If a point fails the honesty gates (grid conv
      <= 5e-4, no boundary flags, dial ladder <= 1e-4, signal/conv >= 10)
      it is re-solved ONCE with escalated budgets (printed); if it still
      fails it is declared UNCONVERGED and excluded from fits — a bounded
      verdict is acceptable, forcing is not.
  (b) cross-check gate: re-solve eps = 0.005 (same budgets) and require
      agreement with J4's banked 3.2195597 within J4's two-resolution band.
  (c) J4's four fit classes on the EXTENDED dataset (K2 points + J4's
      s2_fits dataset), windows extended down (new W0: eps <= 0.01);
      free-intercept variants carry the family-bias question — J4's G2c
      measured the 12-mode family intercept bias +7.7e-5 rel, which is
      ~6.7% of the SIGNAL at eps = 0.001, so fixed-exact-intercept
      exponents are expected to dip; a constant-bias model is TESTED.
  (d) propagation decision: re-fit the slope; propagate through (4.9')
      only if the ceiling moves > J4's family spread (1.3%); otherwise
      state propagation unchanged and cite J4's numbers.

Deterministic (no RNG; fixed NM seeds inherited from i1/axi3).
k2_results.json is flushed after EVERY solve point; on restart, banked
points are reloaded and not re-solved (same code path => same numbers).
WITHIN-MODEL ONLY: nothing here bears on nature.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "k2_results.json")

import i1_r5test as i1          # noqa: E402  (i1 machinery REUSED unmodified)
ax = i1.ax
a3 = i1.a3

SQ2 = math.sqrt(2.0)
PI = math.pi
C00 = 64.0 * SQ2 / (9.0 * PI)            # exact saturated anchor 3.2011247
D_LOCK = 3.0e-3                          # Delta_lock (IV.H.3), j2 verbatim
SIG_FLOOR = 3.0e-5                       # per-point rel sigma floor (j4 verbatim)
EPS_NEW = [0.001, 0.002]                 # the K2 extension
EPS_GATE = 0.005                         # J4 cross-check re-solve
FAM_SPREAD = 0.013                       # J4's 1.3% family spread = propagation trigger
CONV_GATE = 5.0e-4                       # honest grid-convergence requirement (j4 G1e)
LADDER_GATE = 1.0e-4                     # dial-ladder honesty band
SNR_GATE = 10.0                          # signal (c-C00) must be >= 10x grid conv

BUDGET_KEYS = ("N_OPT", "N_POL", "N_FIN", "N_COARSE", "B_TRI", "B_COLD", "B_POL")
BUDGET_DEFAULT = {k: getattr(i1, k) for k in BUDGET_KEYS}
BUDGET_ESC = dict(N_OPT=256, N_POL=416, N_FIN=640, N_COARSE=480,
                  B_TRI=300, B_COLD=2000, B_POL=340)

GATES = []
T0 = time.time()

# ---- restart-safe result store: reload banked points if present ----------
R = {"_notice": "WITHIN-MODEL ONLY; analysis layer; K2 = ROADMAP_v7 "
                "eps<0.005 probe (extends J4)"}
BANK = {}
if os.path.exists(OUT):
    try:
        with open(OUT) as fh:
            old = json.load(fh)
        for p in old.get("s1_new_points", []) + \
                ([old["s2_gate_point"]] if old.get("s2_gate_point") else []):
            BANK[round(p["eps_nominal"], 6)] = p
        if BANK:
            print("RESTART: reloaded %d banked solve point(s) from "
                  "k2_results.json: %s" % (len(BANK), sorted(BANK)))
    except Exception as e:                                   # noqa: BLE001
        print("RESTART: k2_results.json unreadable (%s) — starting fresh" % e)


def gate(name, ok, detail=""):
    GATES.append({"gate": name, "pass": bool(ok), "detail": str(detail)})
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")


def flush():
    R["gate_table"] = GATES
    R["runtime_s"] = time.time() - T0
    with open(OUT, "w") as fh:
        json.dump(R, fh, indent=1, ensure_ascii=False, default=float)


def close(a, b, rtol=1e-9):
    return abs(float(a) - float(b)) <= rtol * max(1.0, abs(float(b)))


# ===========================================================================
# 0. provenance: J4's dataset + its banked eps = 0.005 point
# ===========================================================================
print("== 0. provenance: exact anchor + J4 dataset ==")
gate("G0a exact anchor c_sat(0) = 64*sqrt(2)/(9*pi) = 3.2011247",
     close(C00, 3.2011247, 1e-7), f"= {C00:.7f}")

with open(os.path.join(HERE, "j4_results.json")) as fh:
    J4 = json.load(fh)
j4_pts = {round(p["eps_nominal"], 6): p for p in J4["s1_new_points"]}
j4_005 = j4_pts[0.005]
j4_data = J4["s2_fits"]["dataset"]        # 11 rows, eps 0.005 .. 2.0
j4_W1 = J4["s2_fits"]["windows"]["W1 eps<=0.05 (primary)"]
b_j4 = j4_W1["M1f_free_intercept"]["c0"] / C00 - 1.0     # measured family bias
s_j4f = j4_W1["M1f_free_intercept"]["s_over_c0"]         # free-intercept slope
gate("G0b J4 dataset loaded: banked c_sat(0.005) = 3.2195597, 11-point "
     "fit dataset, measured family intercept bias +7.7e-5",
     close(j4_005["c_sat"], 3.2195597, 1e-6) and len(j4_data) == 11
     and close(b_j4, 7.71e-5, 2e-2),
     f"c(0.005) = {j4_005['c_sat']:.7f}; bias = {b_j4:+.2e}; "
     f"free-intercept slope {s_j4f:.4f}")
gate("G0c i1 budgets are J4's budgets (192/320/480/320, 220/1500/260)",
     tuple(BUDGET_DEFAULT[k] for k in BUDGET_KEYS)
     == (192, 320, 480, 320, 220, 1500, 260),
     str(BUDGET_DEFAULT))
R["s0_provenance"] = {"exact_anchor": C00, "j4_c_005": j4_005["c_sat"],
                      "j4_conv_005": j4_005["conv_c_rel"],
                      "j4_family_bias_rel": b_j4,
                      "j4_free_intercept_slope": s_j4f,
                      "budgets_default": BUDGET_DEFAULT,
                      "budgets_escalated_if_needed": BUDGET_ESC}
flush()


# ===========================================================================
# 1. new saturated solves (J4 code path verbatim, wrapped for escalation)
# ===========================================================================
def solve_point(eps, t_guess, budgets=None, tag=""):
    """One saturated solve point — the j4_epslaw.py stage-1 body verbatim,
    with optional budget override (printed) and honesty verdicts attached."""
    tA = time.time()
    if budgets:
        for k, v in budgets.items():
            setattr(i1, k, v)
        print("   BUDGET CHANGE (printed per K2 charge): %s" % budgets)
    try:
        sol = i1.dial(eps, t_guess, N=4000, rmax=i1.rmax_for(eps))
        lad = i1.dial_ladder(eps, sol)
        R["s1_dial_in_flight"] = {"eps": eps, "t": sol["t"],
                                  "ratio": sol["ratio"],
                                  "ladder": {k: v["d_ratio_rel"]
                                             for k, v in lad.items()}}
        flush()                       # bank the dial before the ~50 s solve
        r, f = sol["r"], sol["f"]
        fspl = ax.CubicSpline1D(r, f)
        idx = np.nonzero(f > 2e-5)[0]
        dom_tail = float(min(r[idx[-1]] * 1.05, r[-1]))
        r_sup_b = float(r[np.nonzero(f > 0.02)[0][-1]])
        dom_sat = float(min(max(1.5 * r_sup_b, dom_tail), r[-1]))
        rep, rep0 = i1.saturated_at(fspl, sol["t"], r_sup_b, dom_sat,
                                    "eps=%.4f%s" % (eps, tag))
    finally:
        if budgets:
            for k, v in BUDGET_DEFAULT.items():
                setattr(i1, k, v)
    floor_c = C00 + 12.0 * sol["t"]
    sig_over_conv = ((rep["c_paper"] - C00)
                     / max(rep["conv_c_rel"] * rep["c_paper"], 1e-300))
    lad_max = max(abs(v["d_ratio_rel"]) for v in lad.values())
    reasons = []
    if rep["conv_c_rel"] > CONV_GATE:
        reasons.append("grid conv %.1e > %.0e" % (rep["conv_c_rel"], CONV_GATE))
    if any(rep["flags"].values()):
        reasons.append("boundary flags %s"
                       % [k for k, v in rep["flags"].items() if v])
    if lad_max > LADDER_GATE:
        reasons.append("dial ladder %.1e > %.0e" % (lad_max, LADDER_GATE))
    if sig_over_conv < SNR_GATE:
        reasons.append("signal/conv %.1f < %.0f" % (sig_over_conv, SNR_GATE))
    pt = dict(eps_nominal=eps, eps_achieved=sol["ratio"], t=sol["t"],
              c_sat=rep["c_paper"], c_sat_coarse=rep["coarse"]["c_paper"],
              conv_c_rel=rep["conv_c_rel"],
              c_sat_hedgehog_core=rep0["c_paper"],
              d=rep["d"], degree=rep["degree"],
              flags={k: bool(v) for k, v in rep["flags"].items()},
              floor_analytic=floor_c,
              floor_ok=bool(rep["c_paper"] >= floor_c - 5e-3),
              dial_ladder={k: {"d_ratio_rel": v["d_ratio_rel"]}
                           for k, v in lad.items()},
              signal_over_conv=sig_over_conv,
              budgets=(budgets or BUDGET_DEFAULT),
              escalated=bool(budgets),
              converged=(not reasons), unconverged_reasons=reasons,
              nfev=rep["nfev"], secs=time.time() - tA)
    print("   eps=%.4f%s: t=%.6g  c_sat=%.6f (coarse %.6f, conv %.1e)  "
          "floor %.5f %s  sig/conv %.0f  ladder %.1e  flags %s  (%.0f s)%s"
          % (eps, tag, sol["t"], pt["c_sat"], pt["c_sat_coarse"],
             pt["conv_c_rel"], floor_c, "OK" if pt["floor_ok"] else "VIOL",
             sig_over_conv, lad_max,
             [k for k, v in pt["flags"].items() if v], pt["secs"],
             "" if pt["converged"] else "  ** UNCONVERGED: %s" % reasons))
    return pt


print("\n== 1. new saturated solves at eps in {0.001, 0.002} + gate re-solve "
      "at 0.005 (i1 dial + saturated_at, J4 budgets; two grid resolutions "
      "N_FIN=480 / N_COARSE=320 per point) ==")
new_pts = []
R["s1_new_points"] = new_pts
R["s2_gate_point"] = None
plan = EPS_NEW + [EPS_GATE]
t_guess = ax.TFROZEN * plan[0] / 0.05
for j, eps in enumerate(plan):
    key = round(eps, 6)
    if key in BANK:
        pt = BANK[key]
        print("   eps=%.4f: BANKED (restart) c_sat=%.6f conv=%.1e %s"
              % (eps, pt["c_sat"], pt["conv_c_rel"],
                 "converged" if pt["converged"] else "UNCONVERGED"))
    else:
        pt = solve_point(eps, t_guess)
        if not pt["converged"]:
            print("   eps=%.4f FAILED honesty gates with J4 budgets — "
                  "escalating ONCE (printed) and re-solving" % eps)
            pt2 = solve_point(eps, t_guess, budgets=BUDGET_ESC, tag="/esc")
            if pt2["converged"] or pt2["conv_c_rel"] < pt["conv_c_rel"]:
                pt = pt2
        if not pt["converged"]:
            print("   eps=%.4f remains UNCONVERGED after escalation — "
                  "declared a LIMIT, excluded from fits (bounded verdict)"
                  % eps)
    if eps == EPS_GATE:
        R["s2_gate_point"] = pt
    else:
        if pt not in new_pts:
            new_pts.append(pt)
    flush()
    if j + 1 < len(plan):
        t_guess = pt["t"] * plan[j + 1] / eps

conv_new = [p for p in new_pts if p["converged"]]
gate("G1a honest convergence at the new points (grid conv <= 5e-4, no "
     "flags, ladder <= 1e-4, signal/conv >= 10)",
     len(conv_new) == len(new_pts),
     "; ".join("eps=%g: conv %.1e, sig/conv %.0f, %s"
               % (p["eps_nominal"], p["conv_c_rel"], p["signal_over_conv"],
                  "OK" if p["converged"] else "UNCONV " +
                  str(p["unconverged_reasons"])) for p in new_pts))
gate("G1b analytic floor c >= 3.2011 + 12 t at every new point",
     all(p["floor_ok"] for p in new_pts),
     "; ".join("eps=%g: c-floor=%+.1e"
               % (p["eps_nominal"], p["c_sat"] - p["floor_analytic"])
               for p in new_pts))
chain = [p["c_sat"] for p in sorted(new_pts, key=lambda q: q["eps_nominal"])]
gate("G1c monotone rising: C00 < c(0.001) < c(0.002) < J4's c(0.005)",
     all(a < b for a, b in zip([C00] + chain, chain + [j4_005["c_sat"]])),
     " < ".join("%.6f" % v for v in [C00] + chain + [j4_005["c_sat"]]))
flush()

# ===========================================================================
# 2. cross-check gate: eps = 0.005 re-solve vs J4's banked value
# ===========================================================================
print("\n== 2. cross-check gate: eps = 0.005 re-solve vs J4 ==")
g5 = R["s2_gate_point"]
rel5 = g5["c_sat"] / j4_005["c_sat"] - 1.0
band5 = j4_005["conv_c_rel"]                       # J4's two-resolution band
gate("G2a eps=0.005 re-solve agrees with J4's banked 3.2195597 within "
     "J4's two-resolution band",
     abs(rel5) <= band5,
     f"k2 = {g5['c_sat']:.7f}, j4 = {j4_005['c_sat']:.7f}, rel = "
     f"{rel5:+.2e} (J4 band {band5:.1e}; k2 band {g5['conv_c_rel']:.1e})")
R["s2_crosscheck"] = {"c_k2": g5["c_sat"], "c_j4": j4_005["c_sat"],
                      "rel_diff": rel5, "band_j4": band5,
                      "band_k2": g5["conv_c_rel"]}
flush()

# ===========================================================================
# 3. J4's four fit classes on the EXTENDED dataset (fit code j4 verbatim)
# ===========================================================================
print("\n== 3. exponent fits on the EXTENDED dataset ==")
data = [(p["eps_achieved"], p["c_sat"], max(p["conv_c_rel"], SIG_FLOOR))
        for p in conv_new]
for row in j4_data:
    data.append((row["eps"], row["c_sat"], max(row["sigma_rel"], SIG_FLOOR)))
data.sort()
E = np.array([d[0] for d in data])
C = np.array([d[1] for d in data])
SC = np.array([d[2] for d in data]) * C
Y = C / C00 - 1.0
SY = SC / C00
print("   dataset: %d points, eps in [%.4g, %.4g] (%d new K2 + %d J4)"
      % (len(E), E[0], E[-1], len(conv_new), len(j4_data)))

WINDOWS = {"W0 eps<=0.01 (K2 small-eps)": E <= 0.0105,
           "W1 eps<=0.05 (J4 primary, extended)": E <= 0.055,
           "W2 eps<=0.03": E <= 0.035,
           "W3 eps<=0.10": E <= 0.105,
           "W4 eps<=0.20": E <= 0.205}


# ---- fit machinery: j4_epslaw.py VERBATIM -------------------------------
def wlsq_amp(x, y, s):
    """weighted LSQ amplitude for y = a*x (fixed zero intercept)."""
    w = 1.0 / s ** 2
    a = float(np.sum(w * x * y) / np.sum(w * x * x))
    chi2 = float(np.sum(w * (y - a * x) ** 2))
    return a, chi2


def wlsq_lin(X, y, s):
    """weighted LSQ for y = X @ beta; returns beta, chi2, cov."""
    w = 1.0 / s
    A = X * w[:, None]
    b = y * w
    beta, *_ = np.linalg.lstsq(A, b, rcond=None)
    chi2 = float(np.sum((b - A @ beta) ** 2))
    cov = np.linalg.inv(A.T @ A)
    return beta, chi2, cov


def fit_free_p(e, y, s, free_intercept=False, pgrid=None):
    """profile chi2 over p for y = a*e^p (+ const); Delta-chi2=1 CI after
    rescaling to chi2/dof = 1 (honest residual-based errors)."""
    if pgrid is None:
        pgrid = np.arange(0.40, 1.6001, 0.0005)
    best, prof = None, []
    for p in pgrid:
        if free_intercept:
            X = np.column_stack([np.ones_like(e), e ** p])
            beta, chi2, _ = wlsq_lin(X, y, s)
            a = beta[1]
        else:
            a, chi2 = wlsq_amp(e ** p, y, s)
        prof.append(chi2)
        if best is None or chi2 < best[1]:
            best = (float(p), chi2, float(a))
    prof = np.array(prof)
    dof = len(e) - (3 if free_intercept else 2)
    scale = max(best[1] / max(dof, 1), 1.0)      # rescale if chi2/dof > 1
    lo = hi = best[0]
    for p, c2 in zip(pgrid, prof):
        if c2 <= best[1] + scale:
            lo = min(lo, float(p))
            hi = max(hi, float(p))
    return dict(p=best[0], p_lo=lo, p_hi=hi, amp=best[2], chi2=best[1],
                dof=dof, chi2_scale=scale)
# ---- end verbatim block --------------------------------------------------


fits = {}
for wname, mask in WINDOWS.items():
    e, y, s = E[mask], Y[mask], SY[mask]
    n = int(mask.sum())
    a1, chi1 = wlsq_amp(e, y, s)
    a23, chi23 = wlsq_amp(e ** (2.0 / 3.0), y, s)
    fp = fit_free_p(e, y, s)
    Xc = np.column_stack([e, e ** (2.0 / 3.0)])
    bc, chic, covc = wlsq_lin(Xc, y, s)
    ci, si = C[mask], SC[mask]
    Xl = np.column_stack([np.ones_like(e), e])
    bl, chil, covl = wlsq_lin(Xl, ci, si)
    fpi = fit_free_p(e, ci, si, free_intercept=True)
    fits[wname] = dict(
        n=n, eps_range=[float(e[0]), float(e[-1])],
        M1_linear={"amp": a1, "chi2": chi1, "chi2_dof": chi1 / max(n - 1, 1)},
        M2_twothirds={"amp": a23, "chi2": chi23,
                      "chi2_dof": chi23 / max(n - 1, 1)},
        M3_free_p=fp,
        M4_combo={"amp_lin": float(bc[0]), "amp_23": float(bc[1]),
                  "chi2": chic,
                  "amp_23_sigma": float(math.sqrt(covc[1, 1])
                                        * math.sqrt(max(chic / max(n-2, 1), 1.0)))},
        M1f_free_intercept={"c0": float(bl[0]), "slope_abs": float(bl[1]),
                            "s_over_c0": float(bl[1] / bl[0]), "chi2": chil,
                            "c0_sigma": float(math.sqrt(covl[0, 0])
                                              * math.sqrt(max(chil / max(n-2, 1), 1.0)))},
        M3f_free_p_free_intercept=fpi)
    print(f"   {wname} (n={n}): M1 a={a1:.4f} chi2={chi1:.1f} | "
          f"M2(2/3) a={a23:.4f} chi2={chi23:.1f} | "
          f"M3 p={fp['p']:.4f} [{fp['p_lo']:.3f},{fp['p_hi']:.3f}] "
          f"a={fp['amp']:.4f} chi2={fp['chi2']:.1f} | "
          f"M4 lin={bc[0]:.3f} 2/3={bc[1]:+.4f} chi2={chic:.1f} | "
          f"M1f c0={bl[0]:.5f} s={bl[1]/bl[0]:.4f} | "
          f"M3f p={fpi['p']:.4f} [{fpi['p_lo']:.3f},{fpi['p_hi']:.3f}]")

# local (pairwise) exponents — model-free (j4 pattern, extended down)
loc = []
for i in range(len(E) - 1):
    if E[i + 1] <= 0.105:
        pl = math.log(Y[i + 1] / Y[i]) / math.log(E[i + 1] / E[i])
        spl = math.sqrt((SY[i] / Y[i]) ** 2 + (SY[i + 1] / Y[i + 1]) ** 2) \
            / abs(math.log(E[i + 1] / E[i]))
        loc.append({"pair": [float(E[i]), float(E[i + 1])], "p_local": pl,
                    "sigma": spl})
        print(f"   local exponent [{E[i]:.4f},{E[i+1]:.4f}]: "
              f"p = {pl:.4f} +/- {spl:.4f}")

# jackknife on the extended primary window (W1) free-p fit
mW1 = WINDOWS["W1 eps<=0.05 (J4 primary, extended)"]
eW, yW, sW = E[mW1], Y[mW1], SY[mW1]
jk = []
for i in range(len(eW)):
    sel = np.arange(len(eW)) != i
    jk.append(fit_free_p(eW[sel], yW[sel], sW[sel])["p"])
jk_spread = (min(jk), max(jk))

fw = fits["W1 eps<=0.05 (J4 primary, extended)"]
f0 = fits["W0 eps<=0.01 (K2 small-eps)"]
p_best = fw["M3_free_p"]["p"]
p_lo, p_hi = fw["M3_free_p"]["p_lo"], fw["M3_free_p"]["p_hi"]
p_all = [p_best, p_lo, p_hi, fw["M3f_free_p_free_intercept"]["p"],
         fw["M3f_free_p_free_intercept"]["p_lo"],
         fw["M3f_free_p_free_intercept"]["p_hi"], *jk_spread]
for wname in fits:
    p_all += [fits[wname]["M3_free_p"]["p"]]
p_env = (min(p_all), max(p_all))

# ---- the family-bias (constant-offset) model, TESTED not assumed --------
# J4's G2c measured the 12-mode family intercept bias b = +7.7e-5 rel and
# the free-intercept slope s = 1.1338 on eps in [0.005, 0.05].  If the
# small-eps dip of fixed-intercept exponents is that bias and nothing else,
# then y_meas(eps) = s*eps + b must hold at 0.001/0.002 with NO refit.
bias_rows = []
for p in conv_new:
    y_m = p["c_sat"] / C00 - 1.0
    y_pred = s_j4f * p["eps_achieved"] + b_j4
    sig = max(p["conv_c_rel"], SIG_FLOOR) * p["c_sat"] / C00
    bias_rows.append({"eps": p["eps_achieved"], "y_meas": y_m,
                      "y_pred_j4_bias_model": y_pred,
                      "resid_sigma": (y_m - y_pred) / sig})
    print("   bias-model test eps=%.4g: y_meas=%.6e vs y_pred=%.6e "
          "(J4 free-intercept law, NO refit) -> resid %+.2f sigma"
          % (p["eps_achieved"], y_m, y_pred, bias_rows[-1]["resid_sigma"]))

gate("G3a extended free-p fit (W1ext): linear class holds; envelope "
     "EXCLUDES 2/3",
     p_env[0] > 2.0 / 3.0 + 0.1,
     f"p = {p_best:.3f}, profile CI [{p_lo:.3f}, {p_hi:.3f}], jackknife "
     f"[{jk_spread[0]:.3f}, {jk_spread[1]:.3f}], all-readings envelope "
     f"[{p_env[0]:.3f}, {p_env[1]:.3f}] vs 2/3 = 0.667")
chi_ratio_W1 = fw["M2_twothirds"]["chi2"] / fw["M1_linear"]["chi2"]
chi_ratio_W0 = f0["M2_twothirds"]["chi2"] / f0["M1_linear"]["chi2"]
gate("G3b pure-2/3 vs pure-linear chi2 ratio, W1ext AND the new W0",
     chi_ratio_W1 > 25.0 and chi_ratio_W0 > 25.0,
     f"chi2(2/3)/chi2(lin): W1ext = {chi_ratio_W1:.0f}, W0 = "
     f"{chi_ratio_W0:.0f} — 2/3 EXCLUDED down to eps = "
     f"{E[0]:.4g}")
c0_f = fw["M1f_free_intercept"]["c0"]
c0_s = fw["M1f_free_intercept"]["c0_sigma"]
gate("G3c free-intercept fit on W1ext consistent with J4's (c0, s) and "
     "the exact anchor",
     abs(c0_f / C00 - 1.0) < 1.5e-3
     and abs(fw["M1f_free_intercept"]["s_over_c0"] / s_j4f - 1.0) < FAM_SPREAD,
     f"c0 = {c0_f:.6f} (rel {c0_f/C00-1:+.2e} vs exact, sigma {c0_s:.1e}); "
     f"s = {fw['M1f_free_intercept']['s_over_c0']:.4f} vs J4 {s_j4f:.4f}")
gate("G3d constant-bias model (J4's measured +7.7e-5, NO refit) matches "
     "the new points within 3 sigma",
     all(abs(r["resid_sigma"]) <= 3.0 for r in bias_rows),
     "; ".join("eps=%.4g: %+.2f sigma" % (r["eps"], r["resid_sigma"])
               for r in bias_rows))

R["s3_fits"] = {"dataset": [{"eps": float(a), "c_sat": float(b),
                             "sigma_rel": float(c / b)}
                            for a, b, c in zip(E, C, SC)],
                "windows": fits, "local_exponents": loc,
                "jackknife_p_W1ext": jk, "p_envelope_raw": list(p_env),
                "bias_model_test": bias_rows,
                "sigma_model": f"per-point sigma = max(grid conv, {SIG_FLOOR}) "
                               "* c; residual-rescaled CIs (j4 verbatim)"}
flush()

# ===========================================================================
# 4. propagation decision (trigger: ceiling moves > J4's 1.3% family spread)
# ===========================================================================
print("\n== 4. (4.9') propagation decision ==")
a_j4 = j4_W1["M1_linear"]["amp"]
fp_j4 = j4_W1["M3_free_p"]
a_ext = fw["M1_linear"]["amp"]
fp_ext = fw["M3_free_p"]
ceil_lin_j4 = D_LOCK / a_j4
ceil_lin_ext = D_LOCK / a_ext
ceil_fp_j4 = (D_LOCK / fp_j4["amp"]) ** (1.0 / fp_j4["p"])
ceil_fp_ext = (D_LOCK / fp_ext["amp"]) ** (1.0 / fp_ext["p"])
mv_lin = ceil_lin_ext / ceil_lin_j4 - 1.0
mv_fp = ceil_fp_ext / ceil_fp_j4 - 1.0
mv_max = max(abs(mv_lin), abs(mv_fp))
# raw-envelope-low-edge ceiling, FOR THE RECORD (see G3d: that edge is the
# measured family intercept bias leaking into fixed-intercept small-eps
# fits, not small-eps physics)
amp_lo, _ = wlsq_amp(eW ** p_env[0], yW, sW)
ceil_env_lo = (D_LOCK / amp_lo) ** (1.0 / p_env[0])
print("   ceiling(linear):  J4 %.6e -> K2ext %.6e  (move %+.2e)"
      % (ceil_lin_j4, ceil_lin_ext, mv_lin))
print("   ceiling(free-p):  J4 %.6e -> K2ext %.6e  (move %+.2e)"
      % (ceil_fp_j4, ceil_fp_ext, mv_fp))
print("   [record only] raw-envelope low edge p=%.3f would give ceiling "
      "%.4e (%+.1f%% vs J4 linear) — an artifact of the measured intercept "
      "bias in fixed-intercept small-eps windows (G3d), not propagated"
      % (p_env[0], ceil_env_lo, 100 * (ceil_env_lo / ceil_lin_j4 - 1.0)))
propagate_needed = mv_max > FAM_SPREAD
gate("G4a propagation trigger: central fitted-slope ceiling movement vs "
     "J4's 1.3% family spread",
     True,
     f"max central movement {mv_max:.2e} "
     f"{'>' if propagate_needed else '<='} {FAM_SPREAD} -> "
     f"{'PROPAGATE' if propagate_needed else 'PROPAGATION UNCHANGED'}")
R["s4_propagation"] = {
    "ceiling_linear_j4": ceil_lin_j4, "ceiling_linear_ext": ceil_lin_ext,
    "ceiling_freep_j4": ceil_fp_j4, "ceiling_freep_ext": ceil_fp_ext,
    "movement_linear": mv_lin, "movement_freep": mv_fp,
    "movement_max": mv_max, "trigger": FAM_SPREAD,
    "propagate_needed": bool(propagate_needed),
    "raw_envelope_low_edge_ceiling_record_only": ceil_env_lo,
    "j4_surviving_ranges_cited": J4["s3_propagation"]["surviving_ranges"]}
if not propagate_needed:
    sr = J4["s3_propagation"]["surviving_ranges"]
    print("   PROPAGATION UNCHANGED — citing J4: ceiling [%.4e, %.4e], "
          "thinning x[%.2f, %.2f], f_top [%.1f, %.1f] MeV, g_min "
          "[%.2e, %.2e]"
          % (sr["ceiling"][0], sr["ceiling"][1], sr["thinning"][0],
             sr["thinning"][1], sr["f_top_MeV"][0], sr["f_top_MeV"][1],
             sr["g_low_edge"][0], sr["g_low_edge"][1]))
flush()

# ===========================================================================
# 5. summary
# ===========================================================================
n_pass = sum(1 for g in GATES if g["pass"])
verdict_scope = ("linear holds to eps = %.4g" % min(
    (p["eps_achieved"] for p in conv_new), default=EPS_GATE))
if len(conv_new) < len(new_pts):
    bad = [p for p in new_pts if not p["converged"]]
    verdict_scope += "; UNCONVERGED at " + ", ".join(
        "eps=%g (%s)" % (p["eps_nominal"], "; ".join(p["unconverged_reasons"]))
        for p in bad)
R["summary"] = {
    "gates": f"{n_pass}/{len(GATES)} PASS",
    "exponent_extended": f"p(W1ext) = {p_best:.3f} (CI [{p_lo:.3f}, "
                         f"{p_hi:.3f}]; raw envelope [{p_env[0]:.3f}, "
                         f"{p_env[1]:.3f}])",
    "verdict_scope": verdict_scope,
    "j4_verdict_status": ("STRENGTHENS (2/3 excluded to smaller eps; "
                          "bias-model gate G3d PASS)" if
                          (n_pass == len(GATES)) else "SEE GATE TABLE"),
    "propagation": ("UNCHANGED — J4's numbers cited"
                    if not propagate_needed else "RE-PROPAGATED (moved)"),
}
flush()
print(f"\n== K2 gates: {n_pass}/{len(GATES)} PASS; k2_results.json written "
      f"({time.time()-T0:.0f} s) ==")
if n_pass != len(GATES):
    raise SystemExit(1)
