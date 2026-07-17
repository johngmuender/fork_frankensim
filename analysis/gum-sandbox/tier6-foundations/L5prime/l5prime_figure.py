#!/usr/bin/env python3
"""L5' summary figure (single PNG, light surface, validated palette).

Panels:
  A  smoothed initial profile vs L5's hard truncation (fix 1, log scale)
  B  near-field arrival histograms, tau_max, axial-beyond region (G2/G3)
  C  same data log-y: hard transverse edge vs smooth axial tail
  D  far-field arrival CDFs at d = 25, T = 40 retry (G4)
  E  KS distance vs detector position (G4 diagnosis: null approached)
  F  sample-maximum trajectory vs resolution ladder (L5 failure mode, gone)
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# validated palette (dataviz reference instance, light mode)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"
C_T = "#2a78d6"   # transverse (slot 1, blue)
C_A = "#008300"   # axial      (slot 2, green)

res = json.load(open(os.path.join(OUTDIR, "L5prime_results.json")))
retry = json.load(open(os.path.join(OUTDIR, "L5prime_g4retry.json")))
raw = np.load(os.path.join(OUTDIR, "L5prime_raw.npz"))
rawR = np.load(os.path.join(OUTDIR, "L5prime_g4retry_raw.npz"))

T_NEAR = res["params"]["T_near_window"]
tau_max = res["G2"]["tau_max"]


def fin(a):
    return a[np.isfinite(a)]


tT = fin(raw["tau_T_near"])
tA = fin(raw["tau_A_near"])
tT16 = tT[tT <= T_NEAR]
tA16 = tA[tA <= T_NEAR]
fT = fin(rawR["tau_T_far"])       # retry: T = 40, balanced coverage
fA = fin(rawR["tau_A_far"])

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8.5,
    "text.color": INK, "axes.edgecolor": BASE, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})

fig, axs = plt.subplots(2, 3, figsize=(13.2, 7.4))
fig.subplots_adjust(hspace=0.44, wspace=0.30, top=0.87, bottom=0.08,
                    left=0.055, right=0.985)

# ---- A: initial profile, smoothed vs hard ----
ax = axs[0, 0]
x = np.linspace(-9.2, -0.8, 2001)
d2 = (x + 5.0) ** 2
hard = np.exp(-d2 / 2.0)
hard[np.abs(x + 5.0) > 3.0] = 0.0
u = np.abs(x + 5.0)
W = np.ones_like(u)
W[u >= 3.5] = 0.0
mseg = (u > 2.5) & (u < 3.5)
tt = (3.5 - u[mseg]) / 1.0
W[mseg] = tt ** 3 * (10 + tt * (-15 + 6 * tt))
smooth = np.exp(-d2 / 2.0) * W
ax.semilogy(x, hard, color=MUTED, lw=1.4, ls="--")
ax.semilogy(x, smooth, color=INK, lw=1.8)
ax.annotate("L5: hard cut at 3σ\n(discontinuous)", xy=(-7.9, 2e-2),
            color=MUTED, fontsize=8, ha="center")
ax.annotate("L5′: S2 taper,\nW = 1 to 2.5σ, 0 at 3.5σ\n(C² edge)",
            xy=(-2.6, 2e-4), color=INK, fontsize=8, ha="center")
ax.set_ylim(1e-8, 3)
ax.set_xlabel("x")
ax.set_ylabel(r"$|f_0(x)|$ (log)")
ax.set_title("A  fix 1: smoothed compact truncation", loc="left",
             fontsize=9.5, color=INK)

# ---- B: near-field histograms ----
ax = axs[0, 1]
bins = np.linspace(0.0, T_NEAR, 129)
ax.hist(tA16, bins=bins, color=C_A, alpha=0.55, label="axial spin",
        edgecolor="none")
ax.hist(tT16, bins=bins, color=C_T, alpha=0.55, label="transverse spin",
        edgecolor="none")
ax.axvline(tau_max, color=INK, lw=1.1, ls="--")
ax.axvspan(tau_max, T_NEAR, color=C_A, alpha=0.08)
ax.annotate(r"$\tau_{max}$ = %.4f" % tau_max, xy=(tau_max, 0.96),
            xycoords=("data", "axes fraction"), xytext=(5, 0),
            textcoords="offset points", color=INK, fontsize=8.5)
ax.annotate("axial beyond cutoff:\n%.1f%% (gate: > 5%%)"
            % (100 * res["G3"]["frac_axial_beyond_tau_max"]),
            xy=(0.97, 0.58), xycoords="axes fraction", ha="right",
            color=C_A, fontsize=8)
ax.set_xlim(0, T_NEAR)
ax.set_xlabel("arrival time at $x = d_{near} = 1$")
ax.set_ylabel("count / bin")
ax.set_title("B  near-field arrivals, N = 2000 each (G2, G3)", loc="left",
             fontsize=9.5, color=INK)
ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.92))

# ---- C: log-y tail ----
ax = axs[0, 2]
ax.hist(tA16, bins=bins, color=C_A, alpha=0.55, edgecolor="none")
ax.hist(tT16, bins=bins, color=C_T, alpha=0.55, edgecolor="none")
ax.set_yscale("log")
ax.axvline(tau_max, color=INK, lw=1.1, ls="--")
ax.set_xlim(0, T_NEAR)
ax.set_ylim(0.7, None)
ax.set_xlabel("arrival time at $x = d_{near} = 1$")
ax.set_ylabel("count / bin (log)")
ax.set_title("C  hard transverse edge vs smooth axial tail", loc="left",
             fontsize=9.5, color=INK)
ax.annotate("zero transverse arrivals in\n(%.3f, 40] on retry horizon"
            % tau_max, xy=(0.97, 0.82), xycoords="axes fraction", ha="right",
            color=C_T, fontsize=8)

# ---- D: far-field CDFs (retry, balanced coverage) ----
ax = axs[1, 0]
for data, c, lab in ((fT, C_T, "transverse"), (fA, C_A, "axial")):
    s = np.sort(data)
    ax.step(s, np.arange(1, s.size + 1) / s.size, where="post", color=c,
            lw=1.8, label="%s (n = %d)" % (lab, s.size))
ax.set_xlabel("arrival time at $x = d_{far} = 25$  (T = 40 retry)")
ax.set_ylabel("empirical CDF")
ax.set_title("D  far field at d = 25: G4 FAIL-AT-d25", loc="left",
             fontsize=9.5, color=INK)
ax.annotate("KS stat %.4f, p = %.1e\ncoverage balanced:\n0 / 88 non-arrivals"
            % (retry["gate"]["ks_stat"], retry["gate"]["ks_p"]),
            xy=(0.97, 0.06), xycoords="axes fraction", ha="right",
            color=INK2, fontsize=8)
ax.legend(frameon=False, loc="upper left")

# ---- E: KS distance vs detector position ----
ax = axs[1, 1]
tab = retry["ks_vs_d_T40"]
items = sorted(tab.values(), key=lambda v: v["d"])
ds = [v["d"] for v in items]
ss = [v["ks_stat"] for v in items]
nT0, nA0 = items[-3]["n_T"], items[-3]["n_A"]
crit = 1.358 * np.sqrt((nT0 + nA0) / (nT0 * nA0))
ax.plot(ds, ss, "-", color=C_T, lw=1.8, zorder=3)
ax.plot(ds, ss, "o", color=C_T, ms=7, zorder=4, markeredgecolor=SURFACE,
        markeredgewidth=1.5)
ax.axhline(crit, color=MUTED, lw=1.0, ls="--")
ax.annotate("p = 0.05 threshold (≈ %.3f)" % crit, xy=(2, crit),
            xytext=(0, -11), textcoords="offset points", color=MUTED,
            fontsize=8)
for v in (items[0], items[4], items[-1]):
    ax.annotate(("p = %.0e" % v["ks_p"]) if v["ks_p"] < 1e-3
                else ("p = %.4f" % v["ks_p"]),
                xy=(v["d"], v["ks_stat"]), xytext=(5, 5),
                textcoords="offset points", color=INK2, fontsize=7.5)
ax.set_ylim(0, 0.148)
ax.set_xlabel("detector position d")
ax.set_ylabel("two-sample KS distance")
ax.set_title("E  KS(d): far-field null approached, not violated",
             loc="left", fontsize=9.5, color=INK)

# ---- F: sample-max trajectory convergence ladder ----
ax = axs[1, 2]
mt = res["G2"]["max_trajectory"]
order = ["Nx2048_dt2e-3", "Nx4096_dt1e-3", "Nx4096_dt5e-4",
         "Nx4096_dt2.5e-4", "Nx8192_dt2.5e-4", "Nx4096_dt1.25e-4"]
labels = ["2048\n2e-3", "4096\n1e-3", "4096\n5e-4", "4096\n2.5e-4\n(ref)",
          "8192\n2.5e-4", "4096\n1.25e-4"]
vals = [mt["tau_by_run"].get(k) for k in order]
xs = np.arange(len(order))
ax.plot(xs, vals, "-", color=C_T, lw=1.8, zorder=3)
ax.plot(xs, vals, "o", color=C_T, ms=7, zorder=4,
        markeredgecolor=SURFACE, markeredgewidth=1.5)
vmid = float(np.mean(vals))
ax.set_ylim(vmid - 4e-6, vmid + 4e-6)
ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=7)
ax.set_xlabel("resolution ladder (Nx / dt)")
ax.set_ylabel("latest-trajectory arrival time")
ax.set_title("F  sample-max trajectory: converged (G2)", loc="left",
             fontsize=9.5, color=INK)
ax.annotate("idx %d: rel spread %.1e\nL5's edge outlier (idx 1249,\n"
            "1.67 ... 7.34 / never): gone" % (mt["index"], mt["spread_rel"]),
            xy=(0.04, 0.30), xycoords="axes fraction", va="top",
            color=INK2, fontsize=8)

fig.suptitle("L5$'$ — Pauli-guidance arrival times with smoothed truncation "
             "and $k_0$ = 2  (F-T6-L5-EXEC-2)", fontsize=12, color=INK,
             x=0.055, ha="left")
fig.text(0.055, 0.905, "Within-model; nothing here bears on nature.",
         fontsize=9, color=INK2)

out = os.path.join(OUTDIR, "L5prime_figure.png")
fig.savefig(out, dpi=150)
print("wrote", out)
