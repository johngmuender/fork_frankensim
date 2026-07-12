#!/usr/bin/env python3
"""Phase G4 analysis: E_int(d), pair-law fit, bond-equation closure.

Deterministic (no RNG; fixed grids and golden-section refinements).
Reads runs/*.json written by twoknot_run; writes twoknot_results.json and
twoknot_fig.png.

Corpus forms used (verbatim from 01-GUM-Omega-Paper-v2.0.1.md):

  (7.5)  V_lin = C_d Q_ij H_ij(X), H = dd(e^{-mu r}/r):
             h_par  = (2 + 2x + x^2) e^{-x} / r^3     [x = mu r]
             h_perp = -(1 + x) e^{-x} / r^3
  Channel theorem (VII.D): V = C_d [cos th TrH + (1 - cos th) n^T H n];
         attractive channel (th = pi, n perp X): V = -C_d h_par;
         strong repulsion (th = pi, n par X); weak repulsion aligned.
  App A units lock: b := B / (2 pi p^2), [B] = E L^2, [C_d] = E L^3,
         p = 0.84 +- 0.03  =>  B_eff = C_d mu is the unique E L^2
         combination; b_eff = C_d mu / (2 pi p^2).
  App H: pair-law domain x >= 1.5, MANDATORY calibration window [1.5, 3]
         (x read in the campaign's bond coordinate x = d / R*).
  VII.D/App I.2: "bond against derived core repulsion"; closed loop: the
         bond equation at b = 42 returns x0 = 1.90 +- 0.05 vs measured
         x0 = 1.92 +- 0.08.

Fit model (attractive channel, calibration window):

    E_int(d) = -C_d h_par(mu d) + B_core exp(-nu d)

with mu FROZEN at 7.7725 (radial referee) and the core exponent nu fitted
(the corpus's "derived core repulsion" closed form is not printed in the
text; a Yukawa-class exponential is assumed and nu is reported).  (C_d,
B_core) enter linearly -> WLS at each nu on a deterministic grid + golden
refinement.

Two E_int estimators, both reported:
  * relaxed  (protocol): E_int(d) = E_relaxed(d) - E_relaxed(d_far),
    far anchor x ~ 4; cross-check E(d_far) - 2 E_single printed as the
    cap-convergence systematic;
  * seed (product ansatz, iteration-free, deterministic): E_seed(d) -
    2 E_seed(single) — the clean tail instrument (no cap noise), carrying
    the product-ansatz bias instead.

Bond closures reported:
  * x0_direct: parabolic minimum of the relaxed anchored E_int samples;
  * x0_fit: argmin of the fitted model;
  * x0(b=42), reading 1 (amplitude substitution): argmin of
    -C_d(42) h_par + fitted core, C_d(42) = 2 pi p^2 * 42 / mu;
  * x0(b=42), reading 2 (running estimator): the x where
    b_eff(x) = -E_int(x) mu / (2 pi p^2 h_par(mu d)) crosses +42.
"""

import glob
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs")

RSTAR = 2.0 ** (5.0 / 6.0)
T_FROZEN = 0.008276434949296802
MU = 1.0 / math.sqrt(2.0 * T_FROZEN)
P_DIP = 0.84
P_DIP_ERR = 0.03
B_CORPUS = 42.0
X0_PRED = (1.90, 0.05)
X0_MEAS = (1.92, 0.08)
WINDOW = (1.5, 3.06)  # mandatory [1.5, 3] + the m=29 point at x = 3.0517
SCHEMES = ("corner", "central4")
FAR_M = 38


def h_par(d):
    d = np.asarray(d, dtype=float)
    x = MU * d
    return (2.0 + 2.0 * x + x * x) * np.exp(-x) / d**3


def load_runs():
    runs = {}
    for fn in sorted(glob.glob(os.path.join(RUNS, "run_*.json"))):
        with open(fn) as f:
            r = json.load(f)
        runs[(r["protocol"]["channel"], r["protocol"]["m"])] = r
    return runs


def conv_sigma(run):
    """Cap-tail systematic proxy: objective drop over the last quarter of
    the recorded descent + last-sample estat drift."""
    ser = run["series"]
    n = len(ser)
    i0 = max(0, n - 1 - max(1, n // 4))
    drop = abs(ser[i0]["obj"] - ser[-1]["obj"])
    est = abs(ser[-2]["estat"] - ser[-1]["estat"]) if n >= 2 else 0.0
    return max(drop, est)


def wls_two(dd, ee, ss, nu):
    """WLS of E = -C h_par + B exp(-nu d) with C > 0, B > 0 enforced
    (attraction amplitude and core repulsion are positive by construction;
    an unconstrained solution violating the signs is rejected as inf).
    Returns (C, B, cov, chi2) or None."""
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
    """Deterministic nu-grid + golden refinement of the 3-parameter model.
    nu is restricted to [1.2 mu, 2.6 mu]: the core must decay strictly
    faster than the pair-law tail, and nu -> mu makes the two basis
    functions collinear (degenerate fit, observed and rejected)."""

    def chi_at(nu):
        r = wls_two(dd, ee, ss, nu)
        return (math.inf if r is None else r[3]), r

    grid = np.linspace(1.2 * MU, 2.6 * MU, 561)
    chis = [chi_at(nu)[0] for nu in grid]
    i = int(np.argmin(chis))
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
        # fall back to the best valid grid point
        best = None
        for nug in grid:
            c2, rg = chi_at(nug)
            if rg is not None and (best is None or c2 < best[0]):
                best = (c2, nug, rg)
        if best is None:
            return None
        chi2, nu, r = best
    C, B, cov, _ = r
    return {"C_d": C, "B_core": B, "nu": nu, "cov_CB": cov.tolist(), "chi2": chi2,
            "dof": max(1, len(dd) - 3)}


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
    """Vertex of the parabola through the 3 points bracketing the sample
    minimum, with a deterministic 1-sigma band from per-point sigmas
    (vertex recomputed under +-sigma shifts of each point, extremes taken;
    2^3 x sign patterns = 8 deterministic evaluations)."""
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
    return -eint * MU / (2.0 * math.pi * P_DIP**2 * h_par(d))


def crossing_x(xs, vals, target):
    """First x where vals crosses target from below-x side (linear interp
    between adjacent samples)."""
    for i in range(len(xs) - 1):
        a, b = vals[i] - target, vals[i + 1] - target
        if a == 0.0:
            return float(xs[i])
        if a * b < 0.0:
            return float(xs[i] + (xs[i + 1] - xs[i]) * (-a) / (b - a))
    return None


def main():
    runs = load_runs()
    single = runs[("single", 0)]
    channels = {ch: sorted(m for (c, m) in runs if c == ch)
                for ch in ("attract", "align", "repulse")}
    far = runs.get(("attract", FAR_M))

    results = {
        "constants": {
            "rstar": RSTAR, "mu": MU, "t_frozen": T_FROZEN, "p_dipole": P_DIP,
            "b_corpus": B_CORPUS, "window_x": [1.5, 3.0],
            "x0_predicted_corpus": list(X0_PRED),
            "x0_measured_corpus": list(X0_MEAS),
        },
        "single_knot": {"sigma_conv": conv_sigma(single)},
        "per_run": [], "fits": {}, "bond": {}, "anchor_crosscheck": {},
    }
    for scheme in SCHEMES:
        results["single_knot"][scheme] = {
            "estat_final": single["final"][scheme]["estat"],
            "estat_seed": single["seed"][scheme]["estat"],
        }
    results["single_knot"]["deg_final"] = single["final"]["corner"]["deg"]
    results["single_knot"]["anf"] = {k: single["anf"][k]
                                     for k in ("status", "iters", "arrests", "seconds")}

    sig_far = conv_sigma(far) if far else 0.0
    table = {}
    for ch, ms in channels.items():
        for m in ms:
            r = runs[(ch, m)]
            row = {
                "channel": ch, "m": m,
                "d": r["protocol"]["d"], "x": r["protocol"]["x_sep"],
                "deg_seed": r["seed"]["corner"]["deg"],
                "deg_final": r["final"]["corner"]["deg"],
                "d_eff_seed": r["seed"]["d_eff"],
                "d_eff_final": r["final"]["d_eff"],
                "arrests": r["anf"]["arrests"], "status": r["anf"]["status"],
                "seconds": r["anf"]["seconds"],
                "sigma_conv": conv_sigma(r),
                "epen_final": r["final"]["corner"]["epen"],
                "efpen_final": r["final"]["corner"]["efpen"],
            }
            for scheme in SCHEMES:
                e1s = single["seed"][scheme]["estat"]
                e1f = single["final"][scheme]["estat"]
                ef = r["final"][scheme]["estat"]
                row[f"E_{scheme}"] = ef
                row[f"Eint2s_{scheme}"] = ef - 2.0 * e1f
                row[f"Eint_seed_{scheme}"] = r["seed"][scheme]["estat"] - 2.0 * e1s
                if far is not None:
                    row[f"Eint_{scheme}"] = ef - far["final"][scheme]["estat"]
            row["sigma_anch"] = math.hypot(row["sigma_conv"], sig_far)
            table[(ch, m)] = row
            results["per_run"].append(row)

    if far is not None:
        for scheme in SCHEMES:
            results["anchor_crosscheck"][scheme] = (
                far["final"][scheme]["estat"]
                - 2.0 * single["final"][scheme]["estat"]
            )

    # ---- channel-theorem tail ratios (seed estimator, corner) --------------
    # Theorem VII.2 large-separation predictions relative to the attractive
    # channel: align/(-attract) = x^2/(2+2x+x^2) [V = C_d TrH],
    # repulse/(-attract) = (2+x)^2/(2+2x+x^2) [V = C_d(-TrH + 2 h_par-part)],
    # x = mu d.
    ratios = []
    for m in (23, 26, 29):
        row_a = table[("attract", m)]
        xmu = MU * row_a["d"]
        pred_align = xmu * xmu / (2.0 + 2.0 * xmu + xmu * xmu)
        pred_rep = (2.0 + xmu) ** 2 / (2.0 + 2.0 * xmu + xmu * xmu)
        ea = row_a["Eint_seed_corner"]
        ratios.append({
            "m": m, "x": row_a["x"],
            "align_over_minus_attract": table[("align", m)]["Eint_seed_corner"] / (-ea),
            "align_predicted": pred_align,
            "repulse_over_minus_attract": table[("repulse", m)]["Eint_seed_corner"] / (-ea),
            "repulse_predicted": pred_rep,
        })
    results["channel_theorem_tail_ratios"] = ratios

    # ---- differential-noise calibration (measured, not assumed) ------------
    # In the mid zone the three channels' far-anchored E_int scatter by
    # ~2.5e-4 with random signs at m = 20..23 while the seed estimator says
    # the physical interaction there is <= 5e-6: that scatter IS the
    # differential cap-convergence noise.  It decays with separation like
    # the interaction itself (all-channel magnitudes track e^{-mu d}), so
    # the noise model is  sigma_dc(d) = RMS_23 * exp(-mu (d - d_23))
    # floored below by 10% relative.
    m_ref = 23
    d_ref = 2.0 * m_ref * 0.09375
    rms_ref = {}
    for scheme in SCHEMES:
        vals = [table[(ch, m_ref)][f"Eint_{scheme}"]
                for ch in ("attract", "align", "repulse")]
        rms_ref[scheme] = float(np.sqrt(np.mean(np.square(vals))))
    results["noise_model"] = {
        "rms_at_m23": rms_ref,
        "form": "sigma_dc(d) = rms_23 * min(1, exp(-mu (d - d_23))), "
                "sigma_i = max(sigma_dc, 0.10 |E_int|)",
    }

    for scheme in SCHEMES:
        att_ms = [m for m in channels["attract"] if m != FAR_M]
        att = [table[("attract", m)] for m in att_ms]
        xx = np.array([r["x"] for r in att])
        dd = np.array([r["d"] for r in att])
        ee_rel = np.array([r[f"Eint_{scheme}"] for r in att])
        sig_dc = rms_ref[scheme] * np.minimum(1.0, np.exp(-MU * (dd - d_ref)))
        ss_rel = np.maximum(sig_dc, 0.10 * np.abs(ee_rel))
        ee_seed = np.array([r[f"Eint_seed_{scheme}"] for r in att])
        ss_seed = np.array([max(0.01 * abs(v), 1e-11) for v in ee_seed])

        inw = (xx >= WINDOW[0] - 1e-9) & (xx <= WINDOW[1])
        fits = {}
        for name, ee, ss in (("relaxed_anchored", ee_rel, ss_rel),
                             ("seed_ansatz", ee_seed, ss_seed)):
            fit = fit_model(dd[inw], ee[inw], ss[inw])
            if fit is None:
                fits[name] = {"failed": True}
                continue
            C = fit["C_d"]
            fit["b_eff"] = C * MU / (2.0 * math.pi * P_DIP**2)
            cerr = math.sqrt(max(fit["cov_CB"][0][0], 0.0))
            fit["C_d_err"] = cerr
            fit["b_eff_err"] = (abs(fit["b_eff"]) * math.hypot(
                cerr / abs(C) if C else math.inf, 2.0 * P_DIP_ERR / P_DIP))
            dmin = golden_min(lambda d_, f=fit: float(
                model_v(d_, f["C_d"], f["B_core"], f["nu"])), 2.2, 7.4)
            fit["x0_fit"] = dmin / RSTAR if dmin else None
            fits[name] = fit

        # direct minima with deterministic bands
        x0_rel, x0_rel_lo, x0_rel_hi = parabola_min_with_band(xx, ee_rel, ss_rel)
        x0_seed, _, _ = parabola_min_with_band(xx, ee_seed, ss_seed)

        # closed loop at b = 42
        cd42 = 2.0 * math.pi * P_DIP**2 * B_CORPUS / MU
        bond42 = {}
        for name, fit in fits.items():
            if fit.get("failed"):
                bond42[name] = None
                continue
            dmin = golden_min(lambda d_, f=fit: float(
                model_v(d_, cd42, f["B_core"], f["nu"])), 1.8, 12.0)
            bond42[name] = dmin / RSTAR if dmin else None
        # reading 2: running-estimator crossing of +42 (relaxed + seed)
        beff_rel = beff_running(ee_rel, dd)
        beff_seed = beff_running(ee_seed, dd)
        x42_rel = crossing_x(xx, beff_rel, B_CORPUS)
        x42_seed = crossing_x(xx, beff_seed, B_CORPUS)

        # free-mu tail check on the seed estimator (pure-attraction points)
        tail = ee_seed < 0
        mu_eff = None
        if int(np.sum(tail)) >= 3:
            dt_, et_ = dd[tail], ee_seed[tail]

            def chi_mu(muv):
                xv = muv * dt_
                b1 = -(2.0 + 2.0 * xv + xv * xv) * np.exp(-xv) / dt_**3
                a = float(np.sum(b1 * et_) / np.sum(b1 * b1))
                return float(np.sum((et_ - a * b1) ** 2 / et_**2)), a

            gridm = np.linspace(4.0, 12.0, 641)
            cm = [chi_mu(m)[0] for m in gridm]
            i = int(np.argmin(cm))
            a, b = gridm[max(0, i - 1)], gridm[min(len(gridm) - 1, i + 1)]
            g = (math.sqrt(5.0) - 1.0) / 2.0
            c, d = b - g * (b - a), a + g * (b - a)
            for _ in range(120):
                if chi_mu(c)[0] < chi_mu(d)[0]:
                    b = d
                else:
                    a = c
                c, d = b - g * (b - a), a + g * (b - a)
            mu_eff = 0.5 * (a + b)

        results["fits"][scheme] = fits
        results["fits"][scheme]["mu_eff_seed_tail"] = mu_eff
        results["bond"][scheme] = {
            "x0_direct_relaxed": x0_rel,
            "x0_direct_relaxed_band": [x0_rel_lo, x0_rel_hi],
            "x0_direct_seed": x0_seed,
            "x0_fit_relaxed": fits["relaxed_anchored"].get("x0_fit"),
            "x0_fit_seed": fits["seed_ansatz"].get("x0_fit"),
            "C_d_at_b42": cd42,
            "x0_b42_amplitude_relaxed": bond42["relaxed_anchored"],
            "x0_b42_amplitude_seed": bond42["seed_ansatz"],
            "x0_b42_running_relaxed": x42_rel,
            "x0_b42_running_seed": x42_seed,
            "beff_running_relaxed": beff_rel.tolist(),
            "beff_running_seed": beff_seed.tolist(),
        }

    with open(os.path.join(HERE, "twoknot_results.json"), "w") as f:
        json.dump(results, f, indent=1, sort_keys=True, default=float)

    # ---- figure --------------------------------------------------------------
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
        colors = {"attract": "#1c6fb8", "align": "#5f9e6e", "repulse": "#c23b22"}
        ax = axes[0]
        for ch, ms in channels.items():
            rows = [table[(ch, m)] for m in ms if m != FAR_M or ch != "attract"]
            xs = [r["x"] for r in rows]
            es = [r["Eint_corner"] for r in rows]
            ax.plot(xs, es, "o-", ms=4, lw=1, color=colors[ch], label=ch)
        ax.axhline(0, color="k", lw=0.6)
        ax.axvspan(1.5, 3.0, alpha=0.08, color="green", label="calibration window")
        ax.set_xlabel("x = d / R*")
        ax.set_ylabel("E_int (corner, far-anchored)")
        ax.set_yscale("symlog", linthresh=1e-5)
        ax.legend(fontsize=8)
        ax.set_title("Two-knot interaction energy by channel")

        ax = axes[1]
        att_ms = [m for m in channels["attract"] if m != FAR_M]
        att = [table[("attract", m)] for m in att_ms]
        xs = np.array([r["x"] for r in att])
        er = np.array([r["Eint_corner"] for r in att])
        ds = np.array([r["d"] for r in att])
        sr = np.maximum(
            rms_ref["corner"] * np.minimum(1.0, np.exp(-MU * (ds - d_ref))),
            0.10 * np.abs(er),
        )
        es = np.array([r["Eint_seed_corner"] for r in att])
        ax.errorbar(xs, er, yerr=sr, fmt="o", ms=4, color="#1c6fb8",
                    label="relaxed, far-anchored")
        ax.plot(xs, es, "s", ms=4, color="#7b3294", label="seed product ansatz")
        fit = results["fits"]["corner"]["relaxed_anchored"]
        xg = np.linspace(1.55, 3.3, 400)
        if not fit.get("failed"):
            ax.plot(xg, model_v(xg * RSTAR, fit["C_d"], fit["B_core"], fit["nu"]),
                    "-", color="#e08214", lw=1.2, label="pair-law (7.5) + core fit")
        b = results["bond"]["corner"]
        if b["x0_direct_relaxed"]:
            ax.axvline(b["x0_direct_relaxed"], color="#1c6fb8", ls="--", lw=0.8)
        ax.axvspan(X0_MEAS[0] - X0_MEAS[1], X0_MEAS[0] + X0_MEAS[1], alpha=0.12,
                   color="purple", label="corpus 1.92 +- 0.08")
        ax.axvline(X0_PRED[0], color="purple", ls=":", lw=1, label="corpus bond eq. 1.90")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlim(1.55, 3.3)
        span = max(1e-4, 3 * abs(min(er.min(), es.min())))
        ax.set_ylim(-span, span * 4)
        ax.set_xlabel("x = d / R*")
        ax.set_ylabel("E_int (corner)")
        ax.set_title("Attractive channel: the bond region")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "twoknot_fig.png"), dpi=140)
    except Exception as ex:
        print("figure generation skipped:", ex)

    # ---- console summary ------------------------------------------------------
    print(f"single: E_corner {single['final']['corner']['estat']:.9f} "
          f"(seed {single['seed']['corner']['estat']:.9f}) "
          f"sigma_conv {results['single_knot']['sigma_conv']:.2e}")
    for scheme in SCHEMES:
        if results["anchor_crosscheck"]:
            print(f"[{scheme}] anchor cross-check E(x~4) - 2 E_single = "
                  f"{results['anchor_crosscheck'][scheme]:+.6e} (cap systematic)")
    hdr = (f"{'ch':8s} {'m':>3s} {'x':>7s} {'Eint_anch':>12s} {'Eint_seed':>12s} "
           f"{'sig':>9s} {'deg_f':>8s} {'deff_f':>7s} {'arr':>4s} {'sec':>5s}")
    print(hdr)
    for row in results["per_run"]:
        print(f"{row['channel']:8s} {row['m']:3d} {row['x']:7.4f} "
              f"{row.get('Eint_corner', float('nan')):12.4e} "
              f"{row['Eint_seed_corner']:12.4e} {row['sigma_anch']:9.2e} "
              f"{row['deg_final']:8.5f} {row['d_eff_final']:7.4f} "
              f"{row['arrests']:4d} {row['seconds']:5.0f}")
    for scheme in SCHEMES:
        for name in ("relaxed_anchored", "seed_ansatz"):
            f_ = results["fits"][scheme][name]
            if f_.get("failed"):
                print(f"[{scheme}/{name}] fit FAILED (no valid nu/positivity)")
                continue
            print(f"[{scheme}/{name}] C_d={f_['C_d']:.4e}+-{f_['C_d_err']:.1e} "
                  f"B_core={f_['B_core']:.4e} nu={f_['nu']:.4f} "
                  f"(nu/mu={f_['nu'] / MU:.3f}) chi2/dof={f_['chi2']:.2f}/{f_['dof']} "
                  f"b_eff={f_['b_eff']:.4g} x0_fit={f_['x0_fit']}")
        print(f"[{scheme}] mu_eff(seed tail) = "
              f"{results['fits'][scheme]['mu_eff_seed_tail']} (frozen {MU:.4f})")
        b_ = results["bond"][scheme]
        print(f"[{scheme}] x0_direct = {b_['x0_direct_relaxed']} "
              f"band {b_['x0_direct_relaxed_band']} | x0_seed {b_['x0_direct_seed']} | "
              f"x0(b=42): ampl {b_['x0_b42_amplitude_relaxed']}/{b_['x0_b42_amplitude_seed']}, "
              f"running {b_['x0_b42_running_relaxed']}/{b_['x0_b42_running_seed']}")


if __name__ == "__main__":
    main()
