#!/usr/bin/env python3
"""L7b figure: spin-family affinity test (4 panels)."""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

OUTDIR = os.path.dirname(os.path.abspath(__file__))
res = json.load(open(os.path.join(OUTDIR, "l7b_results.json")))
raw = np.load(os.path.join(OUTDIR, "l7b_raw.npz"))

NY = raw["ny_values"]
tau = raw["tau_ref"]
Phat = raw["Phat"]
sigma_eff = raw["sigma_eff"]
coef1 = raw["affine_coef"]
ratio = raw["ratio"]
edges = raw["bin_edges"]
T = 16.0
N = 2000

# diverging n_y colormap: blue (-1) -> neutral gray (0) -> orange (+1)
div = LinearSegmentedColormap.from_list(
    "nydiv", ["#2a78d6", "#9a9a97", "#eb6834"])
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3e0"
SURF = "#fcfcfb"

plt.rcParams.update({
    "font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "figure.facecolor": SURF, "axes.facecolor": SURF,
})

fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6))
fig.subplots_adjust(hspace=0.34, wspace=0.26, left=0.07, right=0.97,
                    top=0.90, bottom=0.08)
(axA, axB), (axC, axD) = axes
for ax in axes.flat:
    ax.grid(True, color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

# ---------------- Panel A: arrival-time CDFs by n_y --------------------
for i, ny in enumerate(NY):
    t = tau[i]
    fin = np.sort(t[np.isfinite(t) & (t <= T)])
    y = np.arange(1, fin.size + 1) / N
    c = div((ny + 1) / 2)
    axA.step(np.concatenate([[0], fin, [T]]),
             np.concatenate([[0], y, [y[-1]]]),
             where="post", color=c, lw=1.6)
axA.set_title("A  First-crossing CDF at $x=1$, one curve per $n_y$")
axA.set_xlabel("arrival time $t$")
axA.set_ylabel("cumulative fraction of 2000")
axA.set_xlim(0, T)
axA.set_ylim(0, 1.02)
sm = plt.cm.ScalarMappable(cmap=div, norm=plt.Normalize(-1, 1))
cb = fig.colorbar(sm, ax=axA, pad=0.01, aspect=28)
cb.set_label("$n_y$", rotation=0, labelpad=8)
cb.outline.set_visible(False)
axA.annotate("$\\pm n_y$ pairs coincide:\nstatistics are even in $n_y$",
             xy=(9.5, 0.86), fontsize=8.5, color=INK2)

# ---------------- Panel B: Pi_B(n_y) for selected bins + affine fits ----
# early bin (shape contrast), max-violation bin, kill bin, non-arrival bin
sel_bins = [7, int(res["G2"]["max_violation_bin"]),
            int(res["G3"]["kill_bins"][0]["bin"]), 40]
cat = ["#2a78d6", "#008300", "#eda100", "#eb6834", "#4a3aa7"]
nyf = np.linspace(-1, 1, 100)
for k, b in enumerate(sel_bins):
    c = cat[k % len(cat)]
    lab = ("non-arrival by $T=16$" if b == 40 else
           "$t \\in [%.1f, %.1f]$" % (edges[b], edges[b + 1]))
    axB.errorbar(NY, Phat[:, b], yerr=sigma_eff[:, b], fmt="o", ms=4.5,
                 color=c, lw=1.2, capsize=2, zorder=3, label=lab)
    axB.plot(nyf, coef1[0, b] + coef1[1, b] * nyf, "--", color=c, lw=1.4,
             alpha=0.75, zorder=2)
axB.legend(loc="center right", fontsize=8.5, frameon=False,
           bbox_to_anchor=(1.0, 0.62))
axB.set_title("B  $\\Pi_B(n_y)$ vs the best affine fit $A + B\\,n_y$ (dashed)")
axB.set_xlabel("$n_y$")
axB.set_ylabel("bin probability $\\Pi_B$")
axB.set_xticks(NY)

# ---------------- Panel C: residual/sigma profile ----------------------
mx = ratio.max(axis=0)
centers = np.concatenate([0.5 * (edges[:-1] + edges[1:]), [16.6]])
width = 0.34
over = mx > 5
axC.bar(centers[~over], mx[~over], width=width, color="#2a78d6")
axC.bar(centers[over], mx[over], width=width, color="#eb6834")
axC.axhline(5, color=INK2, lw=1.1, ls=":")
axC.annotate("5$\\sigma$ pre-registered criterion", xy=(0.2, 5.4),
             fontsize=8.5, color=INK2)
axC.annotate("N/A", xy=(16.6, -0.02), xycoords=("data", "axes fraction"),
             ha="center", fontsize=7.5, color=INK2)
axC.set_title("C  Affine-fit violation per bin: "
              "$\\max_{n_y} |r|/\\sigma_{\\rm boot}$ "
              "(%d bins $>5\\sigma$)" % int(res["G2"]["n_bins_over_5sigma"]))
axC.set_xlabel("bin time (last bar = non-arrival)")
axC.set_ylabel("violation ($\\sigma$ units)")
axC.set_yscale("log")
axC.set_ylim(0.05, max(mx) * 2)

# ---------------- Panel D: tau_max(n_y) + kill bin ---------------------
prof = res["G4"]["profile"]
nys = [p["ny"] for p in prof]
tmx = [p["tau_max"] for p in prof]
if res["G3"]["kill_bins"]:
    kb = res["G3"]["kill_bins"][0]
    t0k, t1k = kb["t_range"]
    axD.axhspan(t0k, t1k, color="#eda100", alpha=0.28, zorder=1)
    axD.annotate("kill bin $[%.1f, %.1f]$ (shaded):\n0/2000 arrivals for "
                 "every $|n_y| \\geq %.2f$,\nyet $\\Pi_B = %.3f$ at "
                 "$n_y = %+g$"
                 % (t0k, t1k, kb["nstar"], kb["max_inner_p"],
                    kb["max_inner_ny"]),
                 xy=(-1.0, 10.2), fontsize=8.5, color=INK)
axD.plot(nys, tmx, "-", color="#2a78d6", lw=1.8, zorder=2)
axD.plot(nys, tmx, "o", color="#2a78d6", ms=5.5, zorder=3,
         markeredgecolor=SURF, markeredgewidth=1.2)
axD.set_title("D  Support endpoint $\\tau_{\\max}(n_y)$ "
              "(last arrival of 2000 by $T=16$)")
axD.set_xlabel("$n_y$")
axD.set_ylabel("$\\tau_{\\max}$")
axD.set_xticks(NY)

fig.suptitle("L7b — spin-family affinity test: any single POVM requires "
             "$\\Pi_B(n_y) = A + B\\,n_y$; measured: even, non-affine, "
             "interval-vanishing", fontsize=11.5, fontweight="bold", y=0.97)
fig.text(0.07, 0.925,
         "waveguide $z\\in[0,1]$, $k_0=2$, $N=2000$ per $n_y$, identical "
         "ensemble, seed 20260717; detector $x=1$; within-model",
         fontsize=8.5, color=INK2)
fig.savefig(os.path.join(OUTDIR, "l7b_figure.png"), dpi=160)
print("saved l7b_figure.png")
