#!/usr/bin/env python3
"""Figure for workstream L5: spin-dependent Bohmian arrival times."""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = os.path.dirname(os.path.abspath(__file__))

# palette (dataviz reference instance, light mode)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"
C_T = "#2a78d6"   # categorical slot 1 (blue)  -> transverse spin
C_A = "#008300"   # categorical slot 2 (green) -> axial spin

with open(os.path.join(OUTDIR, "L5_results.json")) as fh:
    R = json.load(fh)
raw = np.load(os.path.join(OUTDIR, "L5_refine_raw.npz"))   # reference resolution
tr = np.load(os.path.join(OUTDIR, "L5_traj.npz"))          # 2-D engine traces

tT = raw["tau_T_near"]; tT = tT[np.isfinite(tT)]
tA = raw["tau_A_near"]; tA = tA[np.isfinite(tA)]
fT = raw["tau_T_far"]; fT = fT[np.isfinite(fT)]
fA = raw["tau_A_far"]; fA = fA[np.isfinite(fA)]
tau_max = R["G2"]["tau_max"]
axial_max = R["G3"]["axial_last_arrival"]
T_FINAL = R["params"]["T"]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "text.color": INK,
    "axes.edgecolor": BASE,
    "axes.labelcolor": INK2,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.7,
    "axes.axisbelow": True,
})

fig = plt.figure(figsize=(10.4, 7.6), facecolor=SURFACE)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05],
                      hspace=0.42, wspace=0.26,
                      left=0.075, right=0.97, top=0.87, bottom=0.08)
axN = fig.add_subplot(gs[0, 0], facecolor=SURFACE)
axF = fig.add_subplot(gs[0, 1], facecolor=SURFACE)
axX = fig.add_subplot(gs[1, :], facecolor=SURFACE)
for ax in (axN, axF, axX):
    ax.spines[["top", "right"]].set_visible(False)

# ---- (a) near-field arrival-time densities ----
bins = np.linspace(0.0, T_FINAL, 121)
axN.hist(tA, bins=bins, density=True, histtype="stepfilled",
         color=C_A, alpha=0.25, lw=0)
axN.hist(tA, bins=bins, density=True, histtype="step", color=C_A, lw=1.6)
axN.hist(tT, bins=bins, density=True, histtype="stepfilled",
         color=C_T, alpha=0.25, lw=0)
axN.hist(tT, bins=bins, density=True, histtype="step", color=C_T, lw=1.6)
axN.axvline(tau_max, color=INK, lw=1.1, ls=(0, (4, 3)))
axN.set_xlim(0, 6.0)
axN.set_yscale("log")
axN.set_ylim(3e-3, 6)
axN.annotate(r"transverse cutoff $\tau_{\max}$ = %.3f" % tau_max,
             (tau_max, 4.5), xytext=(6, 0),
             textcoords="offset points", ha="left", va="top",
             color=INK, fontsize=9)
axN.annotate("zero transverse arrivals\nbeyond the cutoff;\n"
             "axial tail extends to %.3f" % axial_max,
             (3.55, 0.15), ha="left", color=INK2, fontsize=8.5)
axN.set_xlabel(r"arrival time $\tau$ at $x = d_{\rm near} = 1$")
axN.set_ylabel("density")
axN.set_title("(a)  Near field: hard support cutoff", loc="left",
              color=INK, fontsize=10.5)

# ---- (b) far-field arrival-time densities ----
lo = min(fT.min(), fA.min()) - 0.1
binsF = np.linspace(lo, T_FINAL, 61)
axF.hist(fA, bins=binsF, density=True, histtype="stepfilled",
         color=C_A, alpha=0.25, lw=0)
axF.hist(fA, bins=binsF, density=True, histtype="step", color=C_A, lw=1.6)
axF.hist(fT, bins=binsF, density=True, histtype="stepfilled",
         color=C_T, alpha=0.25, lw=0)
axF.hist(fT, bins=binsF, density=True, histtype="step", color=C_T, lw=1.6)
axF.set_xlabel(r"arrival time $\tau$ at $x = d_{\rm far} = 25$")
axF.set_ylabel("density")
axF.set_title("(b)  Far field: distributions converge", loc="left",
              color=INK, fontsize=10.5)
axF.annotate("KS $p$ = %.3f\n($n_T$=%d, $n_A$=%d)"
             % (R["G4"]["ks_p"], R["G4"]["n_far_T"], R["G4"]["n_far_A"]),
             (0.03, 0.95), xycoords="axes fraction", va="top",
             color=INK2, fontsize=9)

# ---- (c) trajectories x(t) ----
ts = tr["ts"]; PT = tr["T"]; PA = tr["A"]   # (nt, np, 2)
for k in range(PT.shape[1]):
    axX.plot(ts, PA[:, k, 0], color=C_A, lw=0.55, alpha=0.45, zorder=2)
for k in range(PT.shape[1]):
    axX.plot(ts, PT[:, k, 0], color=C_T, lw=0.55, alpha=0.45, zorder=3)
axX.axhline(1.0, color=INK, lw=1.0, ls=(0, (4, 3)))
axX.axhline(25.0, color=INK, lw=1.0, ls=(0, (2, 2)))
axX.axvline(tau_max, color=INK, lw=0.9, ls=(0, (4, 3)), alpha=0.6)
axX.annotate(r"$d_{\rm near}=1$", (0.06, 1.0), xytext=(0, 5),
             textcoords="offset points", color=INK, fontsize=9)
axX.annotate(r"$d_{\rm far}=25$", (0.06, 25.0), xytext=(0, 5),
             textcoords="offset points", color=INK, fontsize=9)
axX.annotate(r"$\tau_{\max}$", (tau_max, axX.get_ylim()[0]),
             xytext=(4, 10), textcoords="offset points",
             color=INK2, fontsize=9)
axX.set_xlabel(r"$t$")
axX.set_ylabel(r"$x(t)$")
axX.set_xlim(0, T_FINAL)
axX.set_title("(c)  Bohmian trajectories (80 of 2000 shown, identical initial points)",
              loc="left", color=INK, fontsize=10.5)

# legend (fixed hue order: transverse slot 1, axial slot 2)
handles = [plt.Line2D([], [], color=C_T, lw=2.4,
                      label="transverse spin ($\\sigma_y$): Pauli spin current on"),
           plt.Line2D([], [], color=C_A, lw=2.4,
                      label="axial spin ($\\sigma_x/\\sigma_z$): convective flow only")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.075, 0.955),
           frameon=False, ncol=2, fontsize=9.5, labelcolor=INK2)
fig.suptitle("L5 - Spin-dependent Bohmian arrival times in a hard-wall waveguide "
             "(Das-Dürr class)", x=0.075, y=0.985, ha="left",
             color=INK, fontsize=12, fontweight="bold")

fig.savefig(os.path.join(OUTDIR, "L5_arrival_times.png"), dpi=170,
            facecolor=SURFACE)
print("wrote L5_arrival_times.png")
