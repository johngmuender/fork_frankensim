#!/usr/bin/env python3
"""N2 figure: panel A — rate Lorenz structure (cumulative share of the sample rate
vs top site fraction, 4 h21 anchors); panel B — the honest-rate propagation in
g-equivalent units vs the corpus's printed window and P-nu4 line.
Palette: dataviz reference categorical slots 1-4, validated (validate_palette.js,
light mode: PASS; contrast WARN discharged by direct labels)."""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from n2_secondmoment import (ANCHORS, FM_PER_UM, make_h_table, build_g_table,
                             lorenz_quadrature, G_WINDOW, G_J3_EDGE, G_PNU4)

COLORS = ["#2a78d6", "#008300", "#e87ba4", "#eda100"]
INK = "#0b0b0b"
INK2 = "#52514e"

res = json.load(open("/home/user/fork_frankensim/analysis/gum-sandbox/"
                     "tier7-program/N2/n2_results.json"))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.6, 5.0), dpi=160)
fig.patch.set_facecolor("#fcfcfb")

# ---------------- panel A: Lorenz curves -------------------------------------------
h_tabs = {r0: make_h_table(r0) for r0 in (1.0, 2.0)}
labels_short = ["a=52.6 fm, L=1.0 um", "a=52.6 fm, L=4.2 um",
                "a=6.8 fm, L=4.2 um", "a=2.0 fm, L=4.2 um"]
for i, (a, r0, L_um, _, label) in enumerate(ANCHORS):
    sg, hg = h_tabs[r0]
    dg, gg, dmax = build_g_table(a, r0, sg, hg, extent_r0=100.0)
    L = L_um * FM_PER_UM
    rho_L = 1.0 / L ** 2
    phi = np.pi * a ** 2 / L ** 2
    C, q, J2, _ = lorenz_quadrature(dg, gg, rho_L, a)
    sel = q > 1e-22
    axA.plot(q[sel], C[sel], color=COLORS[i], lw=2, label=labels_short[i])
    axA.plot([phi], [np.interp(a, dg, C)], marker="o", ms=6, color=COLORS[i],
             mec="#fcfcfb", mew=1.2, zorder=5)
axA.set_xscale("log")
axA.set_xlim(1e-20, 1e-11)
axA.set_ylim(0, 1.02)
axA.set_xlabel("top site fraction q (nuclei ranked by rate)", color=INK)
axA.set_ylabel("fraction of total sample rate carried", color=INK)
axA.set_title("A — rate Lorenz structure: the hot-site picture\n"
              "(dots: on-tube fraction $\\phi$; spec window k=2..12 all at 1.000)",
              fontsize=10, color=INK, loc="left")
axA.legend(frameon=False, fontsize=8, loc="upper left")
axA.grid(True, which="major", color="#e8e7e3", lw=0.7)
axA.set_axisbelow(True)
axA.annotate("all of the rate sits on ~$10^{-18}$–$10^{-13}$ of sites;\n"
             "on-tube nuclei alone carry 94–99.6%",
             xy=(3e-15, 0.5), fontsize=8, color=INK2)

# ---------------- panel B: propagation to g-equivalent -----------------------------
y = np.arange(len(ANCHORS))[::-1]
axB.axvspan(G_WINDOW[0], G_WINDOW[1], color="#e8e7e3", zorder=0)
axB.text(1.6e-8, 3.62, "printed g-window\n[8e-10, 1.3e-8]",
         ha="left", fontsize=8, color=INK2)
axB.axvline(G_PNU4, color=INK2, lw=1.2, ls="--")
axB.text(1.3e-9, 1.38, "P-$\\nu$4 line\n$g\\sim10^{-9}$\n(far-horizon)",
         ha="left", fontsize=8, color=INK2)
axB.axvline(1.0, color=INK2, lw=1.2, ls=":")
axB.text(1.45, 3.72, "$g_{equiv}=1$", fontsize=8, color=INK2)
for i, entry in enumerate(res["G3_propagation"]):
    lo = entry["g_equiv"]["window_low_8e-10"]["g_equiv"]
    hi = entry["g_equiv"]["window_top_1.3e-8"]["g_equiv"]
    j3 = entry["g_equiv"]["J3_edge_7.8e-10"]["g_equiv"]
    axB.plot([lo, hi], [y[i], y[i]], color=COLORS[i], lw=6,
             solid_capstyle="round", zorder=3)
    axB.plot([j3], [y[i]], marker="D", ms=6, color=COLORS[i], mec="#fcfcfb",
             mew=1.2, zorder=5)
    axB.annotate(f"$\\sqrt{{R}}$ = {entry['sqrtR']:.2g}", (hi * 2.2, y[i] - 0.06),
                 fontsize=8, color=INK, va="center")
    # connector from printed window to honest band
    axB.annotate("", xy=(lo, y[i] + 0.18), xytext=(G_WINDOW[0], y[i] + 0.18),
                 arrowprops=dict(arrowstyle="->", color=COLORS[i], lw=0.9, alpha=0.6))
axB.set_yticks(y)
axB.set_yticklabels(labels_short, fontsize=8, color=INK)
axB.set_xscale("log")
axB.set_xlim(2e-10, 3e3)
axB.set_ylim(-0.5, 4.1)
axB.set_xlabel("Majoron-mode coupling, g-equivalent units "
               "($g_{equiv}=\\sqrt{R}\\,g$; rate $\\propto g^2$)", color=INK)
axB.set_title("B — honest rate in g-equivalent units: the printed window maps\n"
              "6.9–10.1 orders above the corpus's own P-$\\nu$4 far-horizon line",
              fontsize=10, color=INK, loc="left")
axB.grid(True, which="major", axis="x", color="#e8e7e3", lw=0.7)
axB.set_axisbelow(True)

for ax in (axA, axB):
    ax.set_facecolor("#fcfcfb")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c9c8c2")
    ax.tick_params(colors=INK2, labelsize=8)

fig.suptitle("N2 — F-R9 route (a): the honest second-moment treatment "
             "(diamonds: J3 edge g = 7.8e-10)", fontsize=11, color=INK, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = ("/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/N2/n2_fig.png")
fig.savefig(out, facecolor=fig.get_facecolor())
print("wrote", out)
