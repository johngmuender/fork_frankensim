#!/usr/bin/env python3
"""Phase G5b refine analysis: cap-1600 attractive-channel reruns vs the G4
cap-400 campaign.  Deterministic (no RNG).  Reads runs/ (G4, cap 400) and
runs1600/ (G5b, cap 1600); writes twoknot_refine_results.json and
twoknot_refine_fig.png.  Does NOT touch twoknot_results.json.

What it computes (roadmap G5b step 2-3):

  * per-run tail extrapolation of E_static from the instrumented series:
    Aitken delta^2 on three consecutive block means of the tail (the
    exponential-approach estimator), with a deterministic error bar =
    max(|E_inf - E_last|, block scatter).  Falls back to (E_last, tail
    variation) when the block means are not a contraction (oscillatory
    arrest-dominated tails).
  * anchored E_int at cap 1600 (far anchor m = 38, same cap) and its
    extrapolated twin, each with propagated error bars.
  * THE CAP SYSTEMATIC: E_int(cap 1600) - E_int(cap 400) per separation.
  * a refreshed noise model for the 1600 runs and the well-depth
    significance at x = 1.894.
  * direct parabolic x0 with band; pair-law + core fit on the 1600 points;
    b = 42 running-estimator crossing; all vs corpus targets
    (predicted 1.90 +- 0.05, measured 1.92 +- 0.08).
"""

import glob
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS400 = os.path.join(HERE, "runs")
RUNS1600 = os.path.join(HERE, "runs1600")

RSTAR = 2.0 ** (5.0 / 6.0)
T_FROZEN = 0.008276434949296802
MU = 1.0 / math.sqrt(2.0 * T_FROZEN)
P_DIP = 0.84
P_DIP_ERR = 0.03
B_CORPUS = 42.0
X0_PRED = (1.90, 0.05)
X0_MEAS = (1.92, 0.08)
SCHEME = "corner"  # the protocol instrument (G4 caveat 4)
FAR_M = 38
WELL_M = 18
MS_SCAN = (16, 18, 20, 23)  # the well window


def h_par(d):
    d = np.asarray(d, dtype=float)
    x = MU * d
    return (2.0 + 2.0 * x + x * x) * np.exp(-x) / d**3


def load_runs(dirname):
    runs = {}
    for fn in sorted(glob.glob(os.path.join(dirname, "run_*.json"))):
        with open(fn) as f:
            r = json.load(f)
        runs[(r["protocol"]["channel"], r["protocol"]["m"])] = r
    return runs


def tail_extrapolate(run, frac=0.5):
    """Exponential-approach (Aitken delta^2) extrapolation of E_static from
    the instrumented series.

    The tail (last `frac` of the recorded series) is split into 3 equal
    consecutive blocks with means m1, m2, m3.  For a clean exponential
    approach the block means are a geometric contraction; then
    E_inf = m3 - (m3 - m2)^2 / ((m3 - m2) - (m2 - m1))  [Aitken].
    Validity requires 0 < r < 1 with r = (m3-m2)/(m2-m1).
    Error bar: |E_inf - E_last| + within-block scatter of the last block.
    Fallback (oscillatory / non-contracting tail): E_inf = E_last,
    error = peak-to-peak of the last quarter of the series.
    Returns dict(e_inf, sigma, method, r, e_last, resid_slope)."""
    ser = run["series"]
    es = np.array([e["estat"] for e in ser])
    its = np.array([e["it"] for e in ser])
    e_last = float(es[-1])
    n = len(es)
    i0 = max(0, int(round(n * (1.0 - frac))))
    tail = es[i0:]
    nt = len(tail)
    # residual slope of the tail (per iteration), always reported
    it_t = its[i0:].astype(float)
    slope = float(np.polyfit(it_t, tail, 1)[0]) if nt >= 3 else 0.0
    k = nt // 3
    out = {"e_last": e_last, "resid_slope_per_iter": slope}
    if k >= 2:
        m1 = float(np.mean(tail[nt - 3 * k : nt - 2 * k]))
        m2 = float(np.mean(tail[nt - 2 * k : nt - k]))
        m3 = float(np.mean(tail[nt - k :]))
        d1, d2 = m2 - m1, m3 - m2
        if d1 != 0.0:
            r = d2 / d1
            if 0.0 < r < 0.9:
                e_inf = m3 - d2 * d2 / (d2 - d1) if d2 != d1 else m3
                scat = float(np.std(tail[nt - k :]))
                out.update(
                    e_inf=float(e_inf),
                    sigma=float(abs(e_inf - e_last) + scat),
                    method="aitken",
                    ratio=float(r),
                    blocks=[m1, m2, m3],
                )
                return out
            out["ratio"] = float(r)
    # fallback: oscillatory tail — no net exponential approach resolvable
    q = es[max(0, n - 1 - max(1, n // 4)) :]
    out.update(
        e_inf=e_last,
        sigma=float(np.max(q) - np.min(q)),
        method="last_value_p2p",
    )
    return out


def diff_series(run, far):
    """Differential E_static series E_m(it) - E_far(it) on the common
    instrumentation grid.  The absolute descent (the shared cap systematic)
    cancels iteration-by-iteration; what remains is the physical
    interaction plus arrest-desynchronization noise."""
    fes = {e["it"]: e["estat"] for e in far["series"]}
    its, dv = [], []
    for e in run["series"]:
        if e["it"] in fes:
            its.append(e["it"])
            dv.append(e["estat"] - fes[e["it"]])
    return np.array(its, dtype=float), np.array(dv, dtype=float)


def tail_extrapolate_diff(its, dv, frac=0.5):
    """Cap-extrapolated E_int from the differential series.

    Aitken delta^2 on three consecutive block MEDIANS of the tail (medians,
    not means: arrest-desync spikes — measured at cap 400, e.g. m=23's
    endpoint -2.25e-4 vs a tail plateau of -6e-6 — are outliers, not
    signal).  Validity: contraction 0 < r < 0.9.  Error bar =
    |E_inf - med3| + half the peak-to-peak of the last block.
    Fallback (no contraction): median of the last quarter, error = half
    peak-to-peak of the last quarter + |linear trend over the quarter|."""
    n = len(dv)
    i0 = max(0, int(round(n * (1.0 - frac))))
    tail = dv[i0:]
    it_t = its[i0:]
    nt = len(tail)
    slope = float(np.polyfit(it_t, tail, 1)[0]) if nt >= 3 else 0.0
    out = {"resid_slope_per_iter": slope, "endpoint": float(dv[-1])}
    k = nt // 3
    if k >= 4:
        m1 = float(np.median(tail[nt - 3 * k : nt - 2 * k]))
        m2 = float(np.median(tail[nt - 2 * k : nt - k]))
        m3 = float(np.median(tail[nt - k :]))
        d1, d2 = m2 - m1, m3 - m2
        p2p3 = float(np.max(tail[nt - k :]) - np.min(tail[nt - k :]))
        if d1 != 0.0:
            r = d2 / d1
            out["ratio"] = float(r)
            if 0.0 < r < 0.9:
                e_inf = m3 - d2 * d2 / (d2 - d1)
                out.update(
                    e_inf=float(e_inf),
                    sigma=float(abs(e_inf - m3) + 0.5 * p2p3),
                    method="aitken_diff",
                    blocks=[m1, m2, m3],
                )
                return out
    q0 = max(0, n - max(2, n // 4))
    quarter = dv[q0:]
    med = float(np.median(quarter))
    p2p = float(np.max(quarter) - np.min(quarter))
    trend = abs(float(np.polyfit(its[q0:], quarter, 1)[0])) * (its[-1] - its[q0])
    out.update(
        e_inf=med,
        sigma=float(0.5 * p2p + trend),
        method="median_diff_p2p",
    )
    return out


def wls_two(dd, ee, ss, nu):
    b1 = -h_par(dd)
    b2 = np.exp(-nu * dd)
    s1 = np.max(np.abs(b1))
    s2 = np.max(np.abs(b2))
    A = np.column_stack([b1 / s1, b2 / s2])
    w = 1.0 / ss
    Aw = A * w[:, None]
    ATA = Aw.T @ Aw
    try:
        coef_s = np.linalg.solve(ATA, Aw.T @ (ee * w))
        cov_s = np.linalg.inv(ATA)
    except np.linalg.LinAlgError:
        return None
    scale = np.array([1.0 / s1, 1.0 / s2])
    coef = coef_s * scale
    if coef[0] <= 0.0 or coef[1] <= 0.0:
        return None
    cov = cov_s * np.outer(scale, scale)
    chi2 = float(np.sum(((ee - A @ coef_s) * w) ** 2))
    return float(coef[0]), float(coef[1]), cov, chi2


def fit_model(dd, ee, ss):
    def chi_at(nu):
        r = wls_two(dd, ee, ss, nu)
        return (math.inf if r is None else r[3]), r

    grid = np.linspace(1.2 * MU, 2.6 * MU, 561)
    chis = [chi_at(nu)[0] for nu in grid]
    i = int(np.argmin(chis))
    if not math.isfinite(chis[i]):
        return None
    a, b = grid[max(0, i - 1)], grid[min(len(grid) - 1, i + 1)]
    g = (math.sqrt(5.0) - 1.0) / 2.0
    c, d = b - g * (b - a), a + g * (b - a)
    for _ in range(120):
        if chi_at(c)[0] < chi_at(d)[0]:
            b = d
        else:
            a = c
        c, d = b - g * (b - a), a + g * (b - a)
    nu = 0.5 * (a + b)
    chi2, r = chi_at(nu)
    if r is None:
        best = None
        for nug in grid:
            c2, rg = chi_at(nug)
            if rg is not None and (best is None or c2 < best[0]):
                best = (c2, nug, rg)
        if best is None:
            return None
        chi2, nu, r = best
    C, B, cov, _ = r
    return {"C_d": C, "B_core": B, "nu": nu, "cov_CB": cov.tolist(),
            "chi2": chi2, "dof": max(1, len(dd) - 3)}


def model_v(d, C, B, nu):
    return -C * h_par(d) + B * np.exp(-nu * np.asarray(d, dtype=float))


def golden_min(fun, lo, hi, iters=240):
    g = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    for _ in range(iters):
        if fun(c) < fun(d):
            b = d
        else:
            a = c
        c, d = b - g * (b - a), a + g * (b - a)
    m = 0.5 * (a + b)
    if m - lo < 1e-5 or hi - m < 1e-5:
        return None
    return m


def parabola_vertex(x3, y3):
    (x0, x1, x2), (y0, y1, y2) = x3, y3
    denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
    b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom
    if a <= 0:
        return None
    return -b / (2.0 * a)


def parabola_min_with_band(xs, ys, ss):
    i = int(np.argmin(ys))
    if i == 0 or i == len(xs) - 1:
        return None, None, None
    x3 = (xs[i - 1], xs[i], xs[i + 1])
    y3 = np.array([ys[i - 1], ys[i], ys[i + 1]])
    s3 = np.array([ss[i - 1], ss[i], ss[i + 1]])
    v0 = parabola_vertex(x3, y3)
    if v0 is None:
        return None, None, None
    vs = []
    for sgn in range(8):
        pert = np.array([(1 if sgn >> k & 1 else -1) for k in range(3)]) * s3
        v = parabola_vertex(x3, y3 + pert)
        if v is not None and x3[0] <= v <= x3[2]:
            vs.append(v)
    if vs:
        return v0, min(vs), max(vs)
    return v0, None, None


def beff_running(eint, d):
    return -np.asarray(eint) * MU / (2.0 * math.pi * P_DIP**2 * h_par(d))


def crossing_x(xs, vals, target):
    for i in range(len(xs) - 1):
        a, b = vals[i] - target, vals[i + 1] - target
        if a == 0.0:
            return float(xs[i])
        if a * b < 0.0:
            return float(xs[i] + (xs[i + 1] - xs[i]) * (-a) / (b - a))
    return None


def main():
    r400 = load_runs(RUNS400)
    r1600 = load_runs(RUNS1600)
    far16 = r1600[("attract", FAR_M)]
    far4 = r400[("attract", FAR_M)]
    single16 = r1600[("single", 0)]

    results = {
        "constants": {
            "rstar": RSTAR, "mu": MU, "t_frozen": T_FROZEN, "p_dipole": P_DIP,
            "b_corpus": B_CORPUS,
            "x0_predicted_corpus": list(X0_PRED),
            "x0_measured_corpus": list(X0_MEAS),
            "cap_g4": 400, "cap_g5b": 1600,
        },
        "per_run": [], "tails": {}, "cap_systematic": {},
    }

    # ---- tail extrapolation for every 1600 run (incl. anchor + single) ----
    tails = {}
    for key, run in sorted(r1600.items()):
        t = tail_extrapolate(run)
        tails[key] = t
        results["tails"]["%s_m%02d" % key] = t

    # ---- per-run table: E_int anchored, at both caps + extrapolated -------
    xs, dd = [], []
    e400, e1600, einf, sinf = [], [], [], []
    for m in MS_SCAN:
        run16 = r1600[("attract", m)]
        run4 = r400[("attract", m)]
        d = run16["protocol"]["d"]
        x = run16["protocol"]["x_sep"]
        ei400 = run4["final"][SCHEME]["estat"] - far4["final"][SCHEME]["estat"]
        ei1600 = run16["final"][SCHEME]["estat"] - far16["final"][SCHEME]["estat"]
        t = tails[("attract", m)]
        dits, dv = diff_series(run16, far16)
        td = tail_extrapolate_diff(dits, dv)
        ei_inf = td["e_inf"]
        s_inf = td["sigma"]
        results["tails"]["diff_attract_m%02d" % m] = td
        row = {
            "m": m, "d": d, "x": x,
            "E_1600": run16["final"][SCHEME]["estat"],
            "Eint_400": ei400, "Eint_1600": ei1600,
            "Eint_extrap": ei_inf, "sigma_extrap": s_inf,
            "cap_shift_1600_minus_400": ei1600 - ei400,
            "Eint_seed": (run16["seed"][SCHEME]["estat"]
                          - 2.0 * single16["seed"][SCHEME]["estat"]),
            "deg_final": run16["final"][SCHEME]["deg"],
            "d_eff_final": run16["final"]["d_eff"],
            "arrests": run16["anf"]["arrests"],
            "status": run16["anf"]["status"],
            "seconds": run16["anf"]["seconds"],
            "tail_method": td["method"],
            "resid_slope_per_iter_abs": t["resid_slope_per_iter"],
            "resid_slope_per_iter_diff": td["resid_slope_per_iter"],
        }
        xs.append(x)
        dd.append(d)
        e400.append(ei400)
        e1600.append(ei1600)
        einf.append(ei_inf)
        sinf.append(s_inf)
        results["per_run"].append(row)
    xs = np.array(xs)
    dd = np.array(dd)
    e400 = np.array(e400)
    e1600 = np.array(e1600)
    einf = np.array(einf)
    sinf = np.array(sinf)

    # anchor + single cross-check at 1600
    results["anchor_crosscheck_1600"] = {
        "E_far_minus_2E_single_final": (
            far16["final"][SCHEME]["estat"]
            - 2.0 * single16["final"][SCHEME]["estat"]),
        "value_at_400": (far4["final"][SCHEME]["estat"]
                         - 2.0 * r400[("single", 0)]["final"][SCHEME]["estat"]),
    }

    # ---- the cap systematic ------------------------------------------------
    shift = e1600 - e400
    results["cap_systematic"] = {
        "per_point": {int(m): float(s) for m, s in zip(MS_SCAN, shift)},
        "rms": float(np.sqrt(np.mean(shift**2))),
        "max_abs": float(np.max(np.abs(shift))),
        "note": "Eint(cap1600) - Eint(cap400), far-anchored, corner scheme; "
                "this IS the measured cap systematic of the G4 numbers.",
    }

    # ---- refreshed noise model at cap 1600 ---------------------------------
    # Only the attractive channel exists at 1600, so the G4 cross-channel
    # scatter is not re-measurable.  Three independent handles:
    #   (a) the m=23 point: the seed estimator bounds the physical
    #       interaction there at |E| <= 5e-6, so the measured anchored
    #       value at m=23 is (physical + noise) ~ noise;
    #   (b) the per-run tail-extrapolation error bars (differential,
    #       propagated through the anchor);
    #   (c) the G4->G5b cap shift at m=20,23 (where physics is <= 5e-6):
    #       the shift there is almost pure differential cap noise moving
    #       between the two caps, an upper bound on the 1600 noise.
    m23 = list(MS_SCAN).index(23)
    m20 = list(MS_SCAN).index(20)
    seed23 = abs(results["per_run"][m23]["Eint_seed"])
    noise_a = abs(float(e1600[m23]))
    noise_b = float(np.median(sinf))
    noise_c = float(max(abs(shift[m20]), abs(shift[m23])))
    sigma_ref = max(noise_a, min(noise_b, noise_c))
    d23 = float(dd[m23])
    sig_dc = sigma_ref * np.minimum(1.0, np.exp(-MU * (dd - d23)))
    ss = np.maximum(sig_dc, 0.10 * np.abs(e1600))
    results["noise_model_1600"] = {
        "handle_a_abs_Eint_m23": noise_a,
        "handle_a_seed_bound_m23": seed23,
        "handle_b_median_sigma_extrap": noise_b,
        "handle_c_cap_shift_quiet_zone": noise_c,
        "sigma_ref_adopted": float(sigma_ref),
        "form": "sigma_dc(d) = sigma_ref * min(1, exp(-mu (d - d_23))), "
                "sigma_i = max(sigma_dc, 0.10 |E_int|)",
        "per_point_sigma": {int(m): float(s) for m, s in zip(MS_SCAN, ss)},
    }

    # ---- the well ----------------------------------------------------------
    iw = list(MS_SCAN).index(WELL_M)
    well = {
        "x": float(xs[iw]),
        "Eint_400": float(e400[iw]),
        "Eint_1600": float(e1600[iw]),
        "Eint_extrap": float(einf[iw]),
        "sigma_1600": float(ss[iw]),
        "sigma_extrap": float(sinf[iw]),
        "depth_over_sigma_1600": float(-e1600[iw] / ss[iw]),
        "depth_over_sigma_extrap": float(-einf[iw] / sinf[iw]) if sinf[iw] > 0 else None,
        "g4_significance": 2.1,
    }
    results["well"] = well

    # ---- direct x0 with band (three estimators) ----------------------------
    bond = {}
    x0_16, lo_16, hi_16 = parabola_min_with_band(xs, e1600, ss)
    x0_inf, lo_inf, hi_inf = parabola_min_with_band(xs, einf, np.maximum(sinf, ss))
    bond["x0_direct_1600"] = x0_16
    bond["x0_direct_1600_band"] = [lo_16, hi_16]
    bond["x0_direct_extrap"] = x0_inf
    bond["x0_direct_extrap_band"] = [lo_inf, hi_inf]

    # pair-law + core fit on the four 1600 points (dof = 1)
    fit = fit_model(dd, e1600, ss)
    if fit is not None:
        fit["b_eff"] = fit["C_d"] * MU / (2.0 * math.pi * P_DIP**2)
        dmin = golden_min(lambda d_, f=fit: float(
            model_v(d_, f["C_d"], f["B_core"], f["nu"])), 2.2, 7.4)
        fit["x0_fit"] = dmin / RSTAR if dmin else None
    results["fit_1600"] = fit if fit is not None else {"failed": True}
    bond["x0_fit_1600"] = None if fit is None else fit.get("x0_fit")

    # b = 42 running-estimator crossing (reading B of G4)
    beff16 = beff_running(e1600, dd)
    beffinf = beff_running(einf, dd)
    bond["x0_b42_running_1600"] = crossing_x(xs, beff16, B_CORPUS)
    bond["x0_b42_running_extrap"] = crossing_x(xs, beffinf, B_CORPUS)
    bond["beff_running_1600"] = beff16.tolist()
    bond["beff_running_extrap"] = beffinf.tolist()

    # quoted x0: centre of the estimator ensemble, band = envelope of the
    # direct-parabola sigma bands and the estimator spread
    ests = [v for v in (x0_16, x0_inf, bond["x0_fit_1600"]) if v is not None]
    if ests:
        x0_q = float(np.mean(ests))
        lo_all = [v for v in (lo_16, lo_inf, *ests) if v is not None]
        hi_all = [v for v in (hi_16, hi_inf, *ests) if v is not None]
        bond["x0_quoted"] = x0_q
        bond["x0_quoted_band"] = [float(min(lo_all)), float(max(hi_all))]
    results["bond"] = bond

    with open(os.path.join(HERE, "twoknot_refine_results.json"), "w") as f:
        json.dump(results, f, indent=1, sort_keys=True, default=float)

    # ---- figure -------------------------------------------------------------
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
        ax = axes[0]
        ax.errorbar(xs, e400, fmt="s", ms=4, color="#999999", label="cap 400 (G4)")
        ax.errorbar(xs, e1600, yerr=ss, fmt="o", ms=4, color="#1c6fb8",
                    label="cap 1600")
        ax.errorbar(xs + 0.012, einf, yerr=sinf, fmt="^", ms=4, color="#e08214",
                    label="cap-extrapolated")
        if fit is not None:
            xg = np.linspace(1.6, 2.5, 300)
            ax.plot(xg, model_v(xg * RSTAR, fit["C_d"], fit["B_core"], fit["nu"]),
                    "-", lw=1, color="#e08214", alpha=0.7, label="pair-law+core fit (1600)")
        ax.axhline(0, color="k", lw=0.6)
        ax.axvspan(X0_MEAS[0] - X0_MEAS[1], X0_MEAS[0] + X0_MEAS[1], alpha=0.12,
                   color="purple", label="corpus 1.92 +- 0.08")
        ax.axvline(X0_PRED[0], color="purple", ls=":", lw=1)
        if bond["x0_direct_1600"]:
            ax.axvline(bond["x0_direct_1600"], color="#1c6fb8", ls="--", lw=0.8)
        ax.set_xlabel("x = d / R*")
        ax.set_ylabel("E_int (corner, far-anchored)")
        ax.set_yscale("symlog", linthresh=1e-5)
        ax.set_title("G5b: attractive channel, cap 400 vs 1600 vs extrapolated")
        ax.legend(fontsize=7)

        ax = axes[1]
        for m in MS_SCAN + (FAR_M,):
            run = r1600[("attract", m)]
            ser = run["series"]
            it = [e["it"] for e in ser]
            es = [e["estat"] for e in ser]
            ax.plot(it, np.array(es) - es[-1], lw=0.9,
                    label=f"attract m={m} (x={run['protocol']['x_sep']:.2f})")
        ax.set_yscale("symlog", linthresh=1e-6)
        ax.axvline(400, color="k", ls=":", lw=0.8)
        ax.text(400, ax.get_ylim()[1], " G4 cap", fontsize=7, va="top")
        ax.set_xlabel("ANF iteration")
        ax.set_ylabel("E_static(it) - E_static(1600)")
        ax.set_title("Descent tails (cap 1600)")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "twoknot_refine_fig.png"), dpi=140)
    except Exception as ex:
        print("figure generation skipped:", ex)

    # ---- console summary ----------------------------------------------------
    print("== G5b refine: cap 400 -> 1600, attractive channel ==")
    print(f"anchor cross-check (E_far - 2 E_single): "
          f"{results['anchor_crosscheck_1600']['E_far_minus_2E_single_final']:+.6e} at 1600 "
          f"vs {results['anchor_crosscheck_1600']['value_at_400']:+.6e} at 400")
    print(f"{'m':>3s} {'x':>7s} {'Eint400':>12s} {'Eint1600':>12s} {'shift':>10s} "
          f"{'Eint_extrap':>12s} {'sig_ex':>9s} {'sig1600':>9s} {'tail':>14s}")
    for i, m in enumerate(MS_SCAN):
        row = results["per_run"][i]
        print(f"{m:3d} {xs[i]:7.4f} {e400[i]:12.4e} {e1600[i]:12.4e} "
              f"{shift[i]:+10.2e} {einf[i]:12.4e} {sinf[i]:9.2e} {ss[i]:9.2e} "
              f"{row['tail_method']:>14s}")
    cs = results["cap_systematic"]
    print(f"cap systematic: rms {cs['rms']:.2e}, max |shift| {cs['max_abs']:.2e}")
    nm = results["noise_model_1600"]
    print(f"noise handles: a(m23)={nm['handle_a_abs_Eint_m23']:.2e} "
          f"b(extrap)={nm['handle_b_median_sigma_extrap']:.2e} "
          f"c(shift quiet)={nm['handle_c_cap_shift_quiet_zone']:.2e} "
          f"-> sigma_ref={nm['sigma_ref_adopted']:.2e}")
    w = results["well"]
    print(f"WELL x={w['x']:.4f}: Eint 400={w['Eint_400']:.3e} "
          f"1600={w['Eint_1600']:.3e} extrap={w['Eint_extrap']:.3e}")
    print(f"     depth/sigma: 1600 {w['depth_over_sigma_1600']:.2f} "
          f"extrap {w['depth_over_sigma_extrap']} (G4: 2.1)")
    b = results["bond"]
    print(f"x0 direct 1600 = {b['x0_direct_1600']} band {b['x0_direct_1600_band']}")
    print(f"x0 direct extrap = {b['x0_direct_extrap']} band {b['x0_direct_extrap_band']}")
    print(f"x0 fit 1600 = {b['x0_fit_1600']}")
    print(f"x0 b42 running: 1600 {b['x0_b42_running_1600']} "
          f"extrap {b['x0_b42_running_extrap']}")
    if "x0_quoted" in b:
        print(f"x0 QUOTED = {b['x0_quoted']:.3f} band {b['x0_quoted_band']} "
              f"vs corpus {X0_MEAS[0]} +- {X0_MEAS[1]} (pred {X0_PRED[0]} +- {X0_PRED[1]})")
    if fit is not None:
        print(f"fit 1600: C_d={fit['C_d']:.3e} nu={fit['nu']:.2f} "
              f"(nu/mu={fit['nu']/MU:.2f}) chi2/dof={fit['chi2']:.2f}/{fit['dof']} "
              f"b_eff={fit['b_eff']:.3e}")


if __name__ == "__main__":
    main()
