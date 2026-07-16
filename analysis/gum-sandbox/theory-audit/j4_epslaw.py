#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
J4 — multi-eps saturated re-scan (the one residual open item of Phase J:
J2 defect D3).

QUESTION. J2's repaired (4.8') law c_sat(eps) = 3.2011*(1 + 1.152*eps) is a
LINEAR leading-order class fitted on exactly TWO measured saturated points
(axi4 ring_t0 and ring_eps at eps = 0.05). Does a dense small-eps scan of
the saturated branch support the linear exponent, a 2/3 exponent (the shape
of the corpus's superseded (4.8) law), or something else — and what does
each reading do to the (4.9') ceiling and its downstream edges?

METHOD.
  (a) assemble the measured c_sat(eps) curve from i1_results.json (12-mode
      saturated closure, i1_r5test machinery) and DENSIFY at small eps:
      new solves at eps in {0.005, 0.01, 0.03} + eps = 0.02 as a
      reproduction gate, using i1's dial + saturated_at REUSED UNMODIFIED
      (same budgets: N_OPT/N_POL/N_FIN/N_COARSE = 192/320/480/320,
      B_TRI/B_COLD/B_POL = 220/1500/260; two grid resolutions per point:
      N_FIN = 480 and N_COARSE = 320, plus the dial's own N2x/rmax+ ladder).
  (b) fit (c_sat(eps)/c_sat(0) - 1) against: (i) a*eps, (ii) a*eps^{2/3},
      (iii) a*eps^p free p, (iv) a*eps + b*eps^{2/3}; fixed exact intercept
      c_sat(0) = 64 sqrt2/(9 pi) AND free intercept; window sensitivity.
  (c) propagate each reading through J2's (4.9') arithmetic — formulas
      reused verbatim from j2_epsscan.py: ceiling, Q-6' margins, f-window
      top (0.69*sqrt(ceiling)*1700 MeV), Majoron g low edge (m3/f_top,
      m3 = 0.0468 eV, j3 convention).
  (d) cross-checks: eps = 0.02 vs i1 (must reproduce); eps -> 0 intercept
      vs 64 sqrt2/(9 pi); linear-reading control vs J2/J3 printed numbers.

Deterministic (no RNG; fixed NM seeds inherited from i1/axi3). Stage-flushed
j4_results.json. WITHIN-MODEL ONLY: nothing here bears on nature.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "j4_results.json")

import i1_r5test as i1          # noqa: E402  (i1 machinery REUSED unmodified)
ax = i1.ax
a3 = i1.a3

SQ2 = math.sqrt(2.0)
PI = math.pi
C00 = 64.0 * SQ2 / (9.0 * PI)           # exact saturated anchor 3.2011247
D_LOCK = 3.0e-3                          # Delta_lock (IV.H.3), j2 verbatim
D_BAND = (2.4e-3, 3.8e-3)                # lambda_field (3.1+/-0.7)e-3 omega0
C_G, DC_G = 0.42, 0.04                   # old (4.8) amplitude, j2 verbatim
M_SK = 1.7                               # GeV, c-blind
FQ_LO, FQ_C, FQ_HI = 0.14, 0.17, 0.20    # GeV, sigma-inversion band
M3_EV = 0.0468                           # eV, j3's m3 (dictionary central)
EPS_NEW = [0.005, 0.01, 0.02, 0.03]      # 0.02 = reproduction point
SIG_FLOOR = 3.0e-5                       # per-point rel sigma floor (~grid ladder)

R = {"_notice": "WITHIN-MODEL ONLY; analysis layer; J4 = J2 defect D3 "
                "(multi-eps saturated re-scan)"}
GATES = []
T0 = time.time()


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
# 0. provenance: J2's two points, i1's curve, exact anchor
# ===========================================================================
print("== 0. provenance: J2's two axi4 points vs i1's 12-mode curve ==")
gate("G0a exact anchor c_sat(0) = 64*sqrt(2)/(9*pi) = 3.2011247",
     close(C00, 3.2011247, 1e-7), f"= {C00:.7f}")

with open(os.path.join(i1.T2DIR, "axi4_locked_results.json")) as fh:
    axi4 = json.load(fh)
c0_a4 = axi4["saturated"]["ring_t0"]["c_paper"]
c5_a4 = axi4["saturated"]["ring_eps"]["c_paper"]
gate("G0b J2's two points are axi4 ring_t0/ring_eps = 3.203459 / 3.385470",
     close(c0_a4, 3.2034586, 1e-6) and close(c5_a4, 3.3854702, 1e-6),
     f"{c0_a4:.6f} / {c5_a4:.6f}")
s_j2_exact = (c5_a4 / C00 - 1.0) / 0.05
s_j2_engine = (c5_a4 / c0_a4 - 1.0) / 0.05
gate("G0c J2 slope reproduction: s_exact = 1.152 (engine-base 1.136)",
     close(s_j2_exact, 1.1518, 1e-3) and close(s_j2_engine, 1.1359, 1e-3),
     f"s_exact = {s_j2_exact:.4f}, s_engine = {s_j2_engine:.4f}")

with open(os.path.join(HERE, "i1_results.json")) as fh:
    I1 = json.load(fh)
i1_curve = I1["curve"]
i1_dial = {d["eps"]: d for d in I1["dial"]}
c5_i1 = next(c for c in i1_curve if c["eps"] == 0.05)["c_sat"]
s_i1_05 = (c5_i1 / C00 - 1.0) / i1_dial[0.05]["ratio"]
method_note = (
    "METHOD-CLASS NOTE (measured): J2's linear fit used the two axi4 "
    "ring-saturated engine points (t->0 and eps=0.05). i1's curve is the "
    "12-mode axi3-basis saturated closure on the eps-dial hedgehog base — "
    "same objective G_t, different family. At eps=0.05: axi4 3.385470 vs "
    f"i1 {c5_i1:.6f} (i1 sits {(c5_i1/c5_a4-1):+.2e} rel, the TIGHTER upper "
    "bound); exact-base slope readings 1.152 (axi4) vs "
    f"{s_i1_05:.4f} (i1). Both are family upper bounds on the true G_t*; "
    "within-family consistency is what decides the EXPONENT, family spread "
    "decides the AMPLITUDE band.")
print("   " + method_note)
R["s0_provenance"] = {
    "exact_anchor": C00,
    "j2_points": {"axi4_ring_t0": c0_a4, "axi4_ring_eps0.05": c5_a4,
                  "s_exact_base": s_j2_exact, "s_engine_base": s_j2_engine},
    "i1_at_0.05": {"c_sat": c5_i1, "s_exact_base": s_i1_05},
    "method_class_note": method_note,
    "i1_curve_conv_flags": [{"eps": c["eps"], "conv": c["c_sat_conv"],
                             "floor_ok": c["floor_ok"]} for c in i1_curve]}
flush()

# ===========================================================================
# 1. densify: new saturated-branch solves at small eps (i1 machinery REUSED)
# ===========================================================================
print("\n== 1. new saturated solves (i1 dial + saturated_at, full budgets; "
      "two grid resolutions per point: N_FIN=480 / N_COARSE=320) ==")
print(f"   budgets: N_OPT={i1.N_OPT} N_POL={i1.N_POL} N_FIN={i1.N_FIN} "
      f"N_COARSE={i1.N_COARSE} B_TRI={i1.B_TRI} B_COLD={i1.B_COLD} "
      f"B_POL={i1.B_POL} (no coarsening needed — printed for the record)")
new_pts = []
R["s1_new_points"] = new_pts
t_guess = ax.TFROZEN * EPS_NEW[0] / 0.05
for eps in EPS_NEW:
    tA = time.time()
    sol = i1.dial(eps, t_guess, N=4000, rmax=i1.rmax_for(eps))
    t_guess = sol["t"] * 2.0          # next point is ~2x this eps
    lad = i1.dial_ladder(eps, sol)
    r, f = sol["r"], sol["f"]
    fspl = ax.CubicSpline1D(r, f)
    idx = np.nonzero(f > 2e-5)[0]
    dom_tail = float(min(r[idx[-1]] * 1.05, r[-1]))
    r_sup_b = float(r[np.nonzero(f > 0.02)[0][-1]])
    dom_sat = float(min(max(1.5 * r_sup_b, dom_tail), r[-1]))
    rep, rep0 = i1.saturated_at(fspl, sol["t"], r_sup_b, dom_sat,
                                "eps=%.3f" % eps)
    floor_c = C00 + 12.0 * sol["t"]
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
              nfev=rep["nfev"], secs=time.time() - tA)
    new_pts.append(pt)
    print("   eps=%.4f: t=%.6g  c_sat=%.6f (coarse %.6f, conv %.1e)  "
          "floor %.5f %s  flags %s  (%.0f s)"
          % (eps, sol["t"], pt["c_sat"], pt["c_sat_coarse"],
             pt["conv_c_rel"], floor_c, "OK" if pt["floor_ok"] else "VIOL",
             [k for k, v in pt["flags"].items() if v], pt["secs"]))
    flush()

# gates on the new points
p02 = next(p for p in new_pts if p["eps_nominal"] == 0.02)
c02_i1 = next(c for c in i1_curve if c["eps"] == 0.02)["c_sat"]
rep_rel = p02["c_sat"] / c02_i1 - 1.0
conv02 = max(p02["conv_c_rel"],
             next(c for c in i1_curve if c["eps"] == 0.02)["c_sat_conv"])
gate("G1a REPRODUCTION: c_sat(0.02) vs i1's 3.2739555 within stated "
     "convergence", abs(rep_rel) <= 3.0 * max(conv02, 1e-5),
     f"j4 = {p02['c_sat']:.7f}, i1 = {c02_i1:.7f}, rel = {rep_rel:+.2e} "
     f"(conv band {conv02:.1e})")
gate("G1b dial(0.02) achieved eps vs i1",
     close(p02["eps_achieved"], i1_dial[0.02]["ratio"], 1e-4),
     f"{p02['eps_achieved']:.7f} vs {i1_dial[0.02]['ratio']:.7f}")
gate("G1c all new points respect the analytic floor c >= 3.2011 + 12 t",
     all(p["floor_ok"] for p in new_pts),
     "; ".join(f"eps={p['eps_nominal']}: c-floor={p['c_sat']-p['floor_analytic']:+.1e}"
               for p in new_pts))
gate("G1d monotone rising in eps (new + i1 small-eps rows)",
     all(a["c_sat"] < b["c_sat"] for a, b in zip(new_pts, new_pts[1:]))
     and new_pts[-1]["c_sat"] < c5_i1, "strictly increasing")
gate("G1e in-band grid convergence at every new point (<= 5e-4 rel)",
     all(p["conv_c_rel"] <= 5e-4 for p in new_pts),
     "; ".join(f"{p['eps_nominal']}: {p['conv_c_rel']:.1e}" for p in new_pts))
flush()

# ===========================================================================
# 2. the fit: exponent classes over the small-eps window
# ===========================================================================
print("\n== 2. exponent fits on (c_sat/c_sat(0) - 1) ==")
# dataset: new points (achieved eps) + i1 rows (achieved eps from its dial);
# drop the duplicate 0.02 (keep the j4 re-solve; identical to i1 anyway)
data = [(p["eps_achieved"], p["c_sat"], max(p["conv_c_rel"], SIG_FLOOR))
        for p in new_pts]
for c in i1_curve:
    if c["eps"] >= 0.05:
        data.append((i1_dial[c["eps"]]["ratio"], c["c_sat"],
                     max(c["c_sat_conv"], SIG_FLOOR)))
data.sort()
E = np.array([d[0] for d in data])
C = np.array([d[1] for d in data])
SC = np.array([d[2] for d in data]) * C          # abs sigma on c
Y = C / C00 - 1.0
SY = SC / C00

WINDOWS = {"W1 eps<=0.05 (primary)": E <= 0.055,
           "W2 eps<=0.03": E <= 0.035,
           "W3 eps<=0.10": E <= 0.105,
           "W4 eps<=0.20": E <= 0.205}


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


fits = {}
for wname, mask in WINDOWS.items():
    e, y, s = E[mask], Y[mask], SY[mask]
    n = int(mask.sum())
    a1, chi1 = wlsq_amp(e, y, s)
    a23, chi23 = wlsq_amp(e ** (2.0 / 3.0), y, s)
    fp = fit_free_p(e, y, s)
    Xc = np.column_stack([e, e ** (2.0 / 3.0)])
    bc, chic, covc = wlsq_lin(Xc, y, s)
    # free-intercept variants (family-bias-tolerant): fit c = A + B*e^p
    ci, si = C[mask], SC[mask]
    Xl = np.column_stack([np.ones_like(e), e])
    bl, chil, covl = wlsq_lin(Xl, ci, si)
    fpi = fit_free_p(e, ci, si, free_intercept=True)
    rms = lambda c2, m: math.sqrt(c2 / m)        # chi-units RMS
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

# local (pairwise) exponents — model-free
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

# jackknife on the primary window free-p fit
m1 = WINDOWS["W1 eps<=0.05 (primary)"]
eW, yW, sW = E[m1], Y[m1], SY[m1]
jk = []
for i in range(len(eW)):
    sel = np.arange(len(eW)) != i
    jk.append(fit_free_p(eW[sel], yW[sel], sW[sel])["p"])
jk_spread = (min(jk), max(jk))

fw = fits["W1 eps<=0.05 (primary)"]
p_best = fw["M3_free_p"]["p"]
p_lo, p_hi = fw["M3_free_p"]["p_lo"], fw["M3_free_p"]["p_hi"]
# widen the honest uncertainty by window + jackknife + free-intercept spread
p_all = [p_best, p_lo, p_hi, fw["M3f_free_p_free_intercept"]["p"],
         fw["M3f_free_p_free_intercept"]["p_lo"],
         fw["M3f_free_p_free_intercept"]["p_hi"], *jk_spread]
for wname in fits:
    p_all += [fits[wname]["M3_free_p"]["p"]]
p_env = (min(p_all), max(p_all))

gate("G2a exponent measured: free-p fit on W1, CI + jackknife + window "
     "envelope EXCLUDES 2/3",
     p_env[0] > 2.0 / 3.0 + 0.1,
     f"p = {p_best:.3f}, profile CI [{p_lo:.3f}, {p_hi:.3f}], jackknife "
     f"[{jk_spread[0]:.3f}, {jk_spread[1]:.3f}], all-readings envelope "
     f"[{p_env[0]:.3f}, {p_env[1]:.3f}] vs 2/3 = 0.667")
chi_ratio = fw["M2_twothirds"]["chi2"] / fw["M1_linear"]["chi2"]
gate("G2b pure-2/3 class vs pure-linear class on W1 (chi2 ratio)",
     chi_ratio > 25.0,
     f"chi2(2/3)/chi2(lin) = {chi_ratio:.0f} — 2/3 EXCLUDED as leading order")
c0_fit = fw["M1f_free_intercept"]["c0"]
c0_sig = fw["M1f_free_intercept"]["c0_sigma"]
gate("G2c intercept: free-intercept linear fit eps->0 vs exact 64sqrt2/9pi",
     abs(c0_fit / C00 - 1.0) < 1.5e-3,
     f"c0_fit = {c0_fit:.6f} vs {C00:.6f} (rel {c0_fit/C00-1:+.2e}, "
     f"fit sigma {c0_sig:.1e}); 12-mode family bias at t->0 "
     f"{'+' if c0_fit>C00 else '-'}side, consistent with upper-bound family")

R["s2_fits"] = {"dataset": [{"eps": float(a), "c_sat": float(b),
                             "sigma_rel": float(c / b)}
                            for a, b, c in zip(E, C, SC)],
                "windows": fits, "local_exponents": loc,
                "jackknife_p_W1": jk, "p_envelope": list(p_env),
                "sigma_model": f"per-point sigma = max(grid conv, {SIG_FLOOR}) "
                               "* c; residual-rescaled CIs (chi2/dof -> 1)"}
flush()

# ===========================================================================
# 3. (4.9') propagation per reading — j2/j3 formulas VERBATIM
# ===========================================================================
print("\n== 3. (4.9') propagation per reading (j2_epsscan/j3 formulas) ==")
ceil_old = (D_LOCK / C_G) ** 1.5
eq_lo, eq_c, eq_hi = (FQ_LO / M_SK) ** 2, (FQ_C / M_SK) ** 2, (FQ_HI / M_SK) ** 2
f_lo_MeV = 0.69 * math.sqrt(1.0e-5) * 1700.0


def bisect_root(fun, lo, hi, it=200):
    flo = fun(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if fun(mid) * flo <= 0.0:
            hi = mid
        else:
            lo = mid
            flo = fun(lo)
    return 0.5 * (lo + hi)


def propagate(law, lo_hint=1e-7, hi_hint=0.5):
    """ceiling = eps solving law(eps) = D_lock; then j2/j3 downstream."""
    ceil = bisect_root(lambda e: law(e) - D_LOCK, lo_hint, hi_hint)
    ceil_b = [bisect_root(lambda e: law(e) - D, lo_hint, hi_hint)
              for D in D_BAND]
    marg = [eq_lo / ceil, eq_c / ceil, eq_hi / ceil]
    f_top = 0.69 * math.sqrt(ceil) * 1700.0            # MeV (j2/j3 verbatim)
    g_min = M3_EV / (f_top * 1.0e6)                    # m3/f_top (j3)
    return dict(
        ceiling=ceil, ceiling_Dband=ceil_b, thinning_vs_old=ceil / ceil_old,
        margins_lo_c_hi=marg,
        worst_joint_margin=eq_lo / ceil_b[1],
        breach_further_factor_central=eq_lo / ceil,
        f_window_MeV=[f_lo_MeV, f_top], g_low_edge=g_min,
        window_top_frac_of_3em3=ceil / 3.0e-3)


aW = fw["M1_linear"]["amp"]                 # i1-family linear amp, W1
a23W = fw["M2_twothirds"]["amp"]
fpW = fw["M3_free_p"]
combo = fw["M4_combo"]
readings = {
    "RL_J2_axi4_linear_s1.152 (J2 control)":
        (lambda e: s_j2_exact * e, "SURVIVES (linear class confirmed)"),
    "RL_fit_i1family_linear":
        (lambda e: aW * e, "SURVIVES (measured amplitude, i1 family)"),
    "RP_free_p_best":
        (lambda e: fpW["amp"] * e ** fpW["p"], "SURVIVES (best fit)"),
    "R23_pure_twothirds (EXCLUDED class, reference only)":
        (lambda e: a23W * e ** (2.0 / 3.0), "EXCLUDED by G2a/G2b"),
    "RC_combo_lin_plus_23":
        (lambda e: combo["amp_lin"] * e + combo["amp_23"] * e ** (2.0 / 3.0),
         "combination fit (2/3 admixture measured small)"),
}
# envelope-edge free-p readings (amp refit at pinned p)
amp_lo, _ = wlsq_amp(eW ** p_env[0], yW, sW)
amp_hi, _ = wlsq_amp(eW ** p_env[1], yW, sW)
readings["RP_free_p_envelope_low_p=%.3f" % p_env[0]] = (
    lambda e: amp_lo * e ** p_env[0], "SURVIVES (envelope edge)")
readings["RP_free_p_envelope_high_p=%.3f" % p_env[1]] = (
    lambda e: amp_hi * e ** p_env[1], "SURVIVES (envelope edge)")

prop = {}
for nm, (law, status) in readings.items():
    pr = propagate(law)
    pr["status"] = status
    prop[nm] = pr
    print(f"   {nm}:")
    print(f"      ceiling {pr['ceiling']:.4e} (x{pr['thinning_vs_old']:.2f} "
          f"vs old 6.04e-4); margins [{pr['margins_lo_c_hi'][0]:.2f}, "
          f"{pr['margins_lo_c_hi'][2]:.2f}] c {pr['margins_lo_c_hi'][1]:.2f}; "
          f"worst-joint x{pr['worst_joint_margin']:.2f}; f_top "
          f"{pr['f_window_MeV'][1]:.2f} MeV; g_min {pr['g_low_edge']:.2e}")

ctrl = prop["RL_J2_axi4_linear_s1.152 (J2 control)"]
gate("G3a linear-reading control reproduces J2/J3: ceiling 2.6047e-3, "
     "thinning x4.315, margins [2.60, 5.31] c 3.84, f_top 59.87 MeV, "
     "g_min 7.8e-10",
     close(ctrl["ceiling"], 2.6047e-3, 1e-3)
     and close(ctrl["thinning_vs_old"], 4.315, 2e-3)
     and close(ctrl["margins_lo_c_hi"][0], 2.6037, 1e-3)
     and close(ctrl["margins_lo_c_hi"][2], 5.3137, 1e-3)
     and close(ctrl["f_window_MeV"][1], 59.87, 1e-3)
     and close(ctrl["g_low_edge"], 7.82e-10, 1e-2),
     f"ceiling {ctrl['ceiling']:.4e}, x{ctrl['thinning_vs_old']:.3f}, "
     f"f_top {ctrl['f_window_MeV'][1]:.2f}, g_min {ctrl['g_low_edge']:.2e}")

surv = {k: v for k, v in prop.items() if "EXCLUDED" not in k}
ceil_rng = (min(v["ceiling"] for v in surv.values()),
            max(v["ceiling"] for v in surv.values()))
ftop_rng = (min(v["f_window_MeV"][1] for v in surv.values()),
            max(v["f_window_MeV"][1] for v in surv.values()))
gmin_rng = (min(v["g_low_edge"] for v in surv.values()),
            max(v["g_low_edge"] for v in surv.values()))
thin_rng = (min(v["thinning_vs_old"] for v in surv.values()),
            max(v["thinning_vs_old"] for v in surv.values()))
min_marg_all = min(v["worst_joint_margin"] for v in prop.values())
gate("G3b dichotomy robust: eps_q > ceiling at EVERY reading incl. the "
     "excluded 2/3 reference and every worst-joint",
     min_marg_all > 1.0 and all(1.0e-5 < v["ceiling"] for v in prop.values()),
     f"min worst-joint margin x{min_marg_all:.2f}; lepton window "
     "[1e-5, ceiling] non-empty at every reading")
gate("G3c x4.3 loosening CONFIRMED across surviving readings",
     thin_rng[0] > 3.5 and thin_rng[1] < 5.5,
     f"thinning range x[{thin_rng[0]:.2f}, {thin_rng[1]:.2f}] "
     f"(J2's x4.315 inside)")

R["s3_propagation"] = {
    "formulas": "j2_epsscan.py verbatim: ceiling solves law(eps)=Delta_lock"
                " (=3e-3, band [2.4,3.8]e-3); old ceiling (3e-3/0.42)^1.5;"
                " margins eps_q/ceiling with eps_q=(f_q/1.7)^2,"
                " f_q in [0.14,0.20]; f_top=0.69*sqrt(ceiling)*1700 MeV;"
                " g_min=0.0468 eV/f_top (j3)",
    "ceiling_old": ceil_old, "readings": prop,
    "surviving_ranges": {"ceiling": ceil_rng, "thinning": thin_rng,
                         "f_top_MeV": ftop_rng, "g_low_edge": gmin_rng},
    "robust_vs_reading_specific": {
        "robust": ["confinement dichotomy (all quarks censored, leptons "
                   "free) — holds at every reading incl. excluded 2/3",
                   "lepton floor 1.5e-6, operative bottom 1e-5, f_lo 3.71 "
                   "MeV, all Majoron binding edges (g_max) — ceiling-blind"],
        "now_measured_not_reading_specific": [
            f"x4.3 loosening: x[{thin_rng[0]:.2f}, {thin_rng[1]:.2f}]",
            f"f_top ~60 MeV: [{ftop_rng[0]:.1f}, {ftop_rng[1]:.1f}] MeV",
            f"g_min ~7.8e-10: [{gmin_rng[0]:.2e}, {gmin_rng[1]:.2e}]"],
        "excluded": ["2/3 reading (ceiling 6e-4 class, x11 margin) — "
                     "measured exponent excludes it"]}}
flush()

# ===========================================================================
# 4. summary
# ===========================================================================
n_pass = sum(1 for g in GATES if g["pass"])
R["summary"] = {
    "gates": f"{n_pass}/{len(GATES)} PASS",
    "exponent": f"p = {p_best:.3f} (profile CI [{p_lo:.3f}, {p_hi:.3f}]; "
                f"envelope [{p_env[0]:.3f}, {p_env[1]:.3f}]) — LINEAR class "
                "supported; 2/3 EXCLUDED",
    "law": f"c_sat(eps) = 3.201125*(1 + s*eps), s(i1-family W1) = {aW:.3f} "
           f"(axi4 reading 1.152; family amplitude spread ~"
           f"{abs(s_j2_exact-aW)/s_j2_exact*100:.1f}%)",
    "ceiling_range_surviving": list(ceil_rng),
    "verdict_one_line": "linear exponent measured and confirmed; 2/3 "
                        "excluded; x4.3 loosening + 60 MeV f-top + 7.8e-10 "
                        "Majoron edge promoted from linear-reading-specific "
                        "to measured (with small amplitude-band spread); "
                        "dichotomy robust at every reading",
}
flush()
print(f"\n== J4 gates: {n_pass}/{len(GATES)} PASS; j4_results.json written "
      f"({time.time()-T0:.0f} s) ==")
if n_pass != len(GATES):
    raise SystemExit(1)
