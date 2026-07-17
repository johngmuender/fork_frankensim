#!/usr/bin/env python3
"""L5' summary figure (single PNG, light surface).

Panels:
  A  near-field arrival histograms (transverse vs axial), tau_max line,
     axial-beyond-cutoff region shaded.
  B  same data, log-y: the hard transverse edge vs the smooth axial tail.
  C  far-field arrival CDFs (transverse vs axial) + KS annotation.
  D  sample-maximum trajectory arrival vs resolution ladder (the L5 failure
     mode, shown converged in L5').
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
raw = np.load(os.path.join(OUTDIR, "L5prime_raw.npz"))

T_NEAR = res["params"]["T_near_window"]
tau_max = res["G2"]["tau_max"]

tT = raw["tau_T_near"]; tT = tT[np.isfinite(tT)]
tA = raw["tau_A_near"]; tA = tA[np.isfinite(tA)]
fT = raw["tau_T_far"]; fT = fT[np.isfinite(fT)]
fA = raw["tau_A_far"]; fA = fA[np.isfinite(fA)]
tT16 = tT[tT <= T_NEAR]
tA16 = tA[tA <= T_NEAR]

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9,
    "text.color": INK, "axes.edgecolor": BASE, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
})

fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.2))
fig.subplots_adjust(hspace=0.42, wspace=0.28, top=0.90, bottom=0.08,
                    left=0.075, right=0.975)

bins = np.linspace(0.0, T_NEAR, 129)

# ---- A: near-field histograms ----
ax = axs[0, 0]
ax.hist(tA16, bins=bins, color=C_A, alpha=0.55, label="axial spin",
        edgecolor="none")
ax.hist(tT16, bins=bins, color=C_T, alpha=0.55, label="transverse spin",
        edgecolor="none")
ax.axvline(tau_max, color=INK, lw=1.2, ls="--")
ax.axvspan(tau_max, T_NEAR, color=C_A, alpha=0.08)
ax.annotate(r"$\tau_{max}$ = %.3f" % tau_max, xy=(tau_max, 0.95),
            xycoords=("data", "axes fraction"), xytext=(6, 0),
            textcoords="offset points", color=INK, fontsize=9)
ax.annotate("axial arrivals beyond cutoff: %.1f%%"
            % (100 * res["G3"]["frac_axial_beyond_tau_max"]),
            xy=(0.98, 0.72), xycoords="axes fraction", ha="right",
            color=INK2, fontsize=8.5)
ax.set_xlim(0, T_NEAR)
ax.set_xlabel("arrival time at $x = d_{near} = 1$")
ax.set_ylabel("count / bin")
ax.set_title("A  near-field arrival times (N = 2000 each)", loc="left",
             fontsize=10, color=INK)
ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.68))

# ---- B: log-y tail ----
ax = axs[0, 1]
ax.hist(tA16, bins=bins, color=C_A, alpha=0.55, label="axial",
        edgecolor="none")
ax.hist(tT16, bins=bins, color=C_T, alpha=0.55, label="transverse",
        edgecolor="none")
ax.set_yscale("log")
ax.axvline(tau_max, color=INK, lw=1.2, ls="--")
ax.set_xlim(0, T_NEAR)
ax.set_ylim(0.7, None)
ax.set_xlabel("arrival time at $x = d_{near} = 1$")
ax.set_ylabel("count / bin (log)")
ax.set_title("B  hard transverse edge vs smooth axial tail (log scale)",
             loc="left", fontsize=10, color=INK)
ax.annotate("zero transverse arrivals\nin (%.3f, %g]" % (tau_max, T_NEAR),
            xy=(0.97, 0.80), xycoords="axes fraction", ha="right",
            color=C_T, fontsize=8.5)

# ---- C: far-field CDFs ----
ax = axs[1, 0]
for data, c, lab in ((fT, C_T, "transverse"), (fA, C_A, "axial")):
    s = np.sort(data)
    ax.step(s, np.arange(1, s.size + 1) / s.size, where="post", color=c,
            lw=2.0, label="%s (n = %d)" % (lab, s.size))
ax.set_xlabel("arrival time at $x = d_{far} = 25$")
ax.set_ylabel("empirical CDF")
ax.set_title("C  far-field null: distributions converge", loc="left",
             fontsize=10, color=INK)
ax.annotate("KS stat %.4f,  p = %.3f" % (res["G4"]["ks_stat"],
                                         res["G4"]["ks_p"]),
            xy=(0.97, 0.10), xycoords="axes fraction", ha="right",
            color=INK2, fontsize=9)
ax.legend(frameon=False, loc="upper left")

# ---- D: sample-max trajectory convergence ladder ----
ax = axs[1, 1]
mt = res["G2"]["max_trajectory"]
order = ["Nx2048_dt2e-3", "Nx4096_dt1e-3", "Nx4096_dt5e-4",
         "Nx4096_dt2.5e-4", "Nx8192_dt2.5e-4", "Nx4096_dt1.25e-4"]
labels = ["2048\n2e-3", "4096\n1e-3", "4096\n5e-4", "4096\n2.5e-4\n(ref)",
          "8192\n2.5e-4", "4096\n1.25e-4"]
vals = [mt["tau_by_run"].get(k) for k in order]
xs = np.arange(len(order))
ax.plot(xs, vals, "-", color=C_T, lw=2.0, zorder=3)
ax.plot(xs, vals, "o", color=C_T, ms=8, zorder=4,
        markeredgecolor=SURFACE, markeredgewidth=2)
if mt.get("tau_2d") is not None:
    ax.plot([0], [mt["tau_2d"]], "s", color=INK2, ms=8, zorder=4,
            markeredgecolor=SURFACE, markeredgewidth=2)
    ax.annotate("2-D engine", xy=(0, mt["tau_2d"]), xytext=(8, -12),
                textcoords="offset points", color=INK2, fontsize=8.5)
vmid = np.mean([v for v in vals if v is not None])
ax.set_ylim(vmid - 0.05, vmid + 0.05)
ax.set_xticks(xs)
ax.set_xticklabels(labels, fontsize=7.5)
ax.set_xlabel("resolution (Nx / dt)")
ax.set_ylabel(r"arrival time of latest trajectory")
ax.set_title("D  sample-maximum trajectory: converged across ladder",
             loc="left", fontsize=10, color=INK)
ax.annotate("idx %d, spread %.1e (rel)\nL5's edge outlier: non-convergent\n"
            "(1.67 ... 7.34 / never) — gone" % (mt["index"], mt["spread_rel"]),
            xy=(0.03, 0.96), xycoords="axes fraction", va="top",
            color=INK2, fontsize=8.5)

fig.suptitle("L5$'$ — Pauli-guidance arrival times, smoothed truncation, "
             "$k_0$ = 2  (F-T6-L5-EXEC-2)", fontsize=12, color=INK, x=0.075,
             ha="left")
fig.text(0.075, 0.925, "Within-model; nothing here bears on nature.",
         fontsize=9, color=INK2)

out = os.path.join(OUTDIR, "L5prime_figure.png")
fig.savefig(out, dpi=160)
print("wrote", out)
