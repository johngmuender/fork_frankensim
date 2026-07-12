#!/usr/bin/env python3
"""Phase G4 analysis: E_int(d), pair-law fit, bond-equation closure.

Deterministic (no RNG; fixed algorithms). Reads runs/*.json written by
twoknot_run, writes twoknot_results.json + figures.

Corpus forms used (stated verbatim from 01-GUM-Omega-Paper-v2.0.1.md):

  (7.5)  V_lin = C_d Q_ij H_ij(X),  H = dd(e^{-mu r}/r):
         h_par  = (2 + 2x + x^2) e^{-x} / r^3          [x = mu r]
         h_perp = -(1 + x) e^{-x} / r^3
  Channel theorem: V = C_d [cos(th) TrH + (1 - cos(th)) n^T H n];
         attractive channel (th = pi, n perp X): V = -C_d h_par.
  App A: b := B / (2 pi p^2), [B] = E L^2, [C_d] = E L^3, p = 0.84 +- 0.03
         => B_eff = C_d mu  (the unique EL^2 combination of C_d and mu).
  App H: pair-law domain x >= 1.5, MANDATORY calibration window [1.5, 3]
         (read here, per the campaign's operating convention, in the bond
         coordinate x = d / R*); alpha-Josephson class e^{-2 mu r} taken as
         the core-repulsion form (App C's "derived core repulsion"
         expression is not printed in the corpus text).
  VII.D / App I.2: bond against derived core repulsion; closed loop: the
         bond equation at b = 42 returns x0 = 1.90 +- 0.05 vs measured
         x0 = 1.92 +- 0.08.

Fit model on the calibration window (attractive channel):

  E_int(d) = -C_d * h_par(mu d) + B_rep * exp(-2 mu d)

linear in (C_d, B_rep) -> weighted least squares; x0 = argmin of the model;
closed loop: replace C_d by C_d(b=42) = 2 pi p^2 b / mu, keep the fitted
B_rep, re-solve the bond equation dE/dd = 0.
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
MU = 1.0 / math.sqrt(2.0 * T_FROZEN)  # 7.77254686... (radial referee R4)
P_DIP = 0.84  # p = 0.84 +- 0.03 (Sec VII.C [CAL])
P_DIP_ERR = 0.03
B_CORPUS = 42.0
X0_PRED = (1.90, 0.05)   # corpus bond equation at b = 42
X0_MEAS = (1.92, 0.08)   # corpus measured
WINDOW = (1.5, 3.0)      # mandatory calibration window (App H)
SCHEMES = ("corner", "central4")


def h_par(d):
    x = MU * d
    return (2.0 + 2.0 * x + x * x) * np.exp(-x) / d**3


def rep(d):
    return np.exp(-2.0 * MU * d)


def load_runs():
    runs = {}
    for fn in sorted(glob.glob(os.path.join(RUNS, "run_*.json"))):
        with open(fn) as f:
            r = json.load(f)
        key = (r["protocol"]["channel"], r["protocol"]["m"])
        runs[key] = r
    return runs


def conv_sigma(run, scheme):
    """Tail-flatness proxy: objective drop over the last quarter of the
    descent (upper-bound-flavoured estimate of the remaining cap error),
    plus the drift between the last two recorded estat samples."""
    ser = run["series"]
    n = len(ser)
    i0 = max(0, n - 1 - max(1, n // 4))
    drop = abs(ser[i0]["obj"] - ser[-1]["obj"])
    est_drift = abs(ser[-2]["estat"] - ser[-1]["estat"]) if n >= 2 else 0.0
    return max(drop, est_drift)


def e_final(run, scheme):
    return run["final"][scheme]["estat"]


def weighted_lsq(dd, ee, ss):
    """WLS for E = -C_d h_par(d) + B_rep rep(d); returns coef, cov (2x2)."""
    # scaled basis for conditioning
    b1 = -h_par(dd)
    b2 = rep(dd)
    s1, s2 = np.max(np.abs(b1)), np.max(np.abs(b2))
    A = np.column_stack([b1 / s1, b2 / s2])
    w = 1.0 / ss
    Aw = A * w[:, None]
    yw = ee * w
    ATA = Aw.T @ Aw
    ATy = Aw.T @ yw
    coef_s = np.linalg.solve(ATA, ATy)
    cov_s = np.linalg.inv(ATA)
    scale = np.array([1.0 / s1, 1.0 / s2])
    coef = coef_s * scale
    cov = cov_s * np.outer(scale, scale)
    resid = ee - (A @ coef_s)
    chi2 = float(np.sum((resid * w) ** 2))
    dof = max(1, len(dd) - 2)
    return coef, cov, chi2, dof


def model_e(d, cd, br):
    return -cd * h_par(np.asarray(d, dtype=float)) + br * rep(np.asarray(d, dtype=float))


def argmin_model(cd, br, lo=2.0, hi=7.5):
    """Deterministic golden-section minimization of the model on [lo, hi]."""
    g = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c = b - g * (b - a)
    dd = a + g * (b - a)
    for _ in range(200):
        if model_e(c, cd, br) < model_e(dd, cd, br):
            b = dd
        else:
            a = c
        c = b - g * (b - a)
        dd = a + g * (b - a)
    dm = 0.5 * (a + b)
    # interior check
    if dm - lo < 1e-6 or hi - dm < 1e-6:
        return None
    return dm


def x0_band(coef, cov, lo=2.0, hi=7.5, npts=720):
    """x0 range over the deterministic 1-sigma ellipse of (C_d, B_rep)."""
    vals, vecs = np.linalg.eigh(cov)
    vals = np.clip(vals, 0.0, None)
    x0s = []
    for i in range(npts):
        th = 2.0 * math.pi * i / npts
        dp = vecs @ (np.sqrt(vals) * np.array([math.cos(th), math.sin(th)]))
        cd, br = coef + dp
        if cd <= 0 or br <= 0:
            continue
        dm = argmin_model(cd, br, lo, hi)
        if dm is not None:
            x0s.append(dm / RSTAR)
    if not x0s:
        return None, None
    return float(min(x0s)), float(max(x0s))


def fit_free_mu(dd, ee, ss):
    """One-parameter scan over mu' with WLS amplitude of -h_par(mu' d):
    checks the measured decay constant against the frozen mu."""

    def chi2_of(muv):
        x = muv * dd
        b1 = -(2.0 + 2.0 * x + x * x) * np.exp(-x) / dd**3
        w = 1.0 / ss
        num = np.sum(b1 * ee * w * w)
        den = np.sum(b1 * b1 * w * w)
        a = num / den
        r = (ee - a * b1) * w
        return float(np.sum(r * r)), float(a)

    grid = np.linspace(2.0, 14.0, 481)
    chis = [chi2_of(m)[0] for m in grid]
    i = int(np.argmin(chis))
    lo, hi = grid[max(0, i - 1)], grid[min(len(grid) - 1, i + 1)]
    g = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    for _ in range(120):
        if chi2_of(c)[0] < chi2_of(d)[0]:
            b = d
        else:
            a = c
        c, d = b - g * (b - a), a + g * (b - a)
    mubest = 0.5 * (a + b)
    chi, amp = chi2_of(mubest)
    return mubest, amp, chi


def parabola_min(xs, ys):
    """Quadratic through the 3 points bracketing the minimum sample."""
    i = int(np.argmin(ys))
    if i == 0 or i == len(xs) - 1:
        return None
    x0, x1, x2 = xs[i - 1], xs[i], xs[i + 1]
    y0, y1, y2 = ys[i - 1], ys[i], ys[i + 1]
    denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
    b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / denom
    if a <= 0:
        return None
    return -b / (2.0 * a)


def main():
    runs = load_runs()
    single = runs[("single", 0)]
    far_key = ("attract", 38)

    channels = {}
    for ch in ("attract", "align", "repulse"):
        ms = sorted(m for (c, m) in runs if c == ch)
        channels[ch] = ms

    results = {
        "constants": {
            "rstar": RSTAR,
            "mu": MU,
            "t_frozen": T_FROZEN,
            "p_dipole": P_DIP,
            "b_corpus": B_CORPUS,
            "window_x": list(WINDOW),
            "x0_predicted_corpus": list(X0_PRED),
            "x0_measured_corpus": list(X0_MEAS),
        },
        "single_knot": {},
        "per_run": [],
        "fits": {},
        "bond": {},
    }

    for scheme in SCHEMES:
        results["single_knot"][scheme] = {
            "estat_final": e_final(single, scheme),
            "estat_seed": single["seed"][scheme]["estat"],
        }
    results["single_knot"]["deg_final"] = single["final"]["corner"]["deg"]
    results["single_knot"]["sigma_conv"] = conv_sigma(single, "corner")
    results["single_knot"]["anf"] = {
        k: single["anf"][k] for k in ("status", "iters", "arrests", "seconds")
    }

    table = {}
    for ch, ms in channels.items():
        for m in ms:
            r = runs[(ch, m)]
            d = r["protocol"]["d"]
            x = r["protocol"]["x_sep"]
            row = {
                "channel": ch,
                "m": m,
                "d": d,
                "x": x,
                "deg_seed": r["seed"]["corner"]["deg"],
                "deg_final": r["final"]["corner"]["deg"],
                "d_eff_seed": r["seed"]["d_eff"],
                "d_eff_final": r["final"]["d_eff"],
                "arrests": r["anf"]["arrests"],
                "status": r["anf"]["status"],
                "seconds": r["anf"]["seconds"],
                "sigma_conv": conv_sigma(r, "corner"),
                "epen_final": r["final"]["corner"]["epen"],
                "efpen_final": r["final"]["corner"]["efpen"],
            }
            for scheme in SCHEMES:
                e1 = e_final(single, scheme)
                e2k = e_final(r, scheme)
                row[f"E_{scheme}"] = e2k
                row[f"Eint_{scheme}"] = e2k - 2.0 * e1
                row[f"Eint_seed_{scheme}"] = (
                    r["seed"][scheme]["estat"] - 2.0 * single["seed"][scheme]["estat"]
                )
            table[(ch, m)] = row
            results["per_run"].append(row)

    # anchor cross-check: E(far attract) - 2 E_single
    if far_key in runs:
        results["anchor_crosscheck"] = {
            scheme: table[far_key][f"Eint_{scheme}"] for scheme in SCHEMES
        }

    # ---- fits + bond equation (attractive channel) -------------------------
    for scheme in SCHEMES:
        att = [table[("attract", m)] for m in channels["attract"]]
        dd_all = np.array([r["d"] for r in att])
        xx_all = dd_all / RSTAR
        ee_all = np.array([r[f"Eint_{scheme}"] for r in att])
        sig_single = results["single_knot"]["sigma_conv"]
        ss_all = np.array(
            [
                max(
                    math.hypot(r["sigma_conv"], 2.0 * sig_single),
                    0.02 * abs(r[f"Eint_{scheme}"]),
                    1e-12,
                )
                for r in att
            ]
        )
        inw = (xx_all >= WINDOW[0] - 1e-9) & (xx_all <= WINDOW[1] + 0.06)
        dd, ee, ss = dd_all[inw], ee_all[inw], ss_all[inw]

        coef, cov, chi2, dof = weighted_lsq(dd, ee, ss)
        cd_fit, br_fit = float(coef[0]), float(coef[1])
        cd_err = float(math.sqrt(max(cov[0, 0], 0.0)))
        br_err = float(math.sqrt(max(cov[1, 1], 0.0)))

        beff = cd_fit * MU  # B_eff = C_d mu (App A dimensional lock)
        b_small = beff / (2.0 * math.pi * P_DIP**2)
        b_small_err = b_small * math.sqrt(
            (cd_err / abs(cd_fit)) ** 2 + (2.0 * P_DIP_ERR / P_DIP) ** 2
        ) if cd_fit != 0 else float("nan")

        d0 = argmin_model(cd_fit, br_fit)
        x0_fit = d0 / RSTAR if d0 else None
        x0_lo, x0_hi = x0_band(coef, cov) if d0 else (None, None)

        # direct minimum from measured points (all attract points incl. x<1.5)
        x0_direct = parabola_min(xx_all, ee_all)

        # closed loop at b = 42
        cd_42 = 2.0 * math.pi * P_DIP**2 * B_CORPUS / MU
        d0_42 = argmin_model(cd_42, br_fit)
        x0_42 = d0_42 / RSTAR if d0_42 else None

        mu_free, amp_free, chi_free = fit_free_mu(dd, ee, ss)

        results["fits"][scheme] = {
            "window_points_x": [float(v) for v in xx_all[inw]],
            "C_d": cd_fit,
            "C_d_err": cd_err,
            "B_rep": br_fit,
            "B_rep_err": br_err,
            "chi2": chi2,
            "dof": dof,
            "B_eff=C_d*mu": beff,
            "b_eff=B_eff/(2 pi p^2)": b_small,
            "b_eff_err": b_small_err,
            "mu_free_fit": mu_free,
            "mu_free_amp": amp_free,
            "mu_free_chi2": chi_free,
        }
        results["bond"][scheme] = {
            "x0_fit_min": x0_fit,
            "x0_fit_band": [x0_lo, x0_hi],
            "x0_direct_parabola": x0_direct,
            "x0_closed_loop_b42": x0_42,
            "C_d_at_b42": cd_42,
        }

    with open(os.path.join(HERE, "twoknot_results.json"), "w") as f:
        json.dump(results, f, indent=1, sort_keys=True)

    # ---- figures ------------------------------------------------------------
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
        colors = {"attract": "#1c6fb8", "align": "#8a8a8a", "repulse": "#c23b22"}
        ax = axes[0]
        for ch, ms in channels.items():
            xs = [table[(ch, m)]["x"] for m in ms]
            es = [table[(ch, m)]["Eint_corner"] for m in ms]
            ax.plot(xs, es, "o-", ms=4, lw=1, color=colors[ch], label=ch)
        ax.axhline(0, color="k", lw=0.6)
        ax.axvspan(*WINDOW, alpha=0.08, color="green", label="calibration window")
        ax.set_xlabel("x = d / R*")
        ax.set_ylabel("E_int (corner)")
        ax.set_yscale("symlog", linthresh=1e-6)
        ax.legend(fontsize=8)
        ax.set_title("Two-knot interaction energy by channel")

        ax = axes[1]
        att = [table[("attract", m)] for m in channels["attract"]]
        xs = np.array([r["x"] for r in att])
        es = np.array([r["Eint_corner"] for r in att])
        ax.plot(xs, es, "o", color="#1c6fb8", label="measured (attract)")
        fit = results["fits"]["corner"]
        xg = np.linspace(1.3, 4.1, 400)
        ax.plot(
            xg,
            model_e(xg * RSTAR, fit["C_d"], fit["B_rep"]),
            "-",
            color="#e08214",
            lw=1.2,
            label="pair-law fit (7.5) + e^{-2mu d}",
        )
        b = results["bond"]["corner"]
        if b["x0_fit_min"]:
            ax.axvline(b["x0_fit_min"], color="#e08214", ls="--", lw=0.8)
        ax.axvspan(X0_MEAS[0] - X0_MEAS[1], X0_MEAS[0] + X0_MEAS[1], alpha=0.12,
                   color="purple", label="corpus x0 = 1.92 +- 0.08")
        ax.axvline(X0_PRED[0], color="purple", ls=":", lw=1,
                   label="corpus bond eq. 1.90")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel("x = d / R*")
        ax.set_ylabel("E_int (corner)")
        ax.set_title("Attractive channel: bond")
        ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, "twoknot_fig.png"), dpi=140)
    except Exception as ex:  # figures are convenience, not protocol
        print("figure generation skipped:", ex)

    # ---- console summary -----------------------------------------------------
    print(f"single knot: E_corner = {e_final(single, 'corner'):.9f}  "
          f"E_central4 = {e_final(single, 'central4'):.9f}  "
          f"sigma_conv = {results['single_knot']['sigma_conv']:.2e}")
    if "anchor_crosscheck" in results:
        print("anchor cross-check E(x~4) - 2 E_single:",
              {k: f"{v:.3e}" for k, v in results["anchor_crosscheck"].items()})
    print(f"{'ch':8s} {'m':>3s} {'x':>7s} {'Eint_corner':>13s} {'Eint_c4':>13s} "
          f"{'deg_f':>8s} {'deff_f':>7s} {'sig':>9s} {'arr':>4s}")
    for row in results["per_run"]:
        print(f"{row['channel']:8s} {row['m']:3d} {row['x']:7.4f} "
              f"{row['Eint_corner']:13.4e} {row['Eint_central4']:13.4e} "
              f"{row['deg_final']:8.5f} {row['d_eff_final']:7.4f} "
              f"{row['sigma_conv']:9.2e} {row['arrests']:4d}")
    for scheme in SCHEMES:
        f_ = results["fits"][scheme]
        b_ = results["bond"][scheme]
        print(f"[{scheme}] C_d = {f_['C_d']:.4e} +- {f_['C_d_err']:.1e}, "
              f"B_rep = {f_['B_rep']:.4e} +- {f_['B_rep_err']:.1e}, "
              f"chi2/dof = {f_['chi2']:.2f}/{f_['dof']}")
        print(f"[{scheme}] b_eff = {f_['b_eff=B_eff/(2 pi p^2)']:.4g} "
              f"(corpus 42 +- 6); mu_free = {f_['mu_free_fit']:.4f} "
              f"(frozen {MU:.4f})")
        print(f"[{scheme}] x0_fit = {b_['x0_fit_min']}, band = {b_['x0_fit_band']}, "
              f"x0_direct = {b_['x0_direct_parabola']}, "
              f"x0(b=42) = {b_['x0_closed_loop_b42']}")


if __name__ == "__main__":
    main()
