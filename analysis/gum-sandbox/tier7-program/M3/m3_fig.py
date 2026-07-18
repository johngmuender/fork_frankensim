#!/usr/bin/env python3
"""M3 figure: the pricing of F-R9 route (b). Reads m3_results.json."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M3/"
r = json.load(open(BASE + "m3_results.json"))
scans = r["G2_scans"]

# fixed categorical order (Okabe-Ito, CVD-safe)
COLORS = ["#0072B2", "#E69F00", "#009E73", "#CC79A7"]
INK = "#1a1a1a"
MUTED = "#666666"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
    "axes.grid": True, "grid.color": "#d9d9d9", "grid.linewidth": 0.6,
    "grid.alpha": 0.6, "axes.axisbelow": True})

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.4, 8.2), sharex=True,
                               gridspec_kw={"hspace": 0.12})

# ---- panel A: R(A_sat) ------------------------------------------------------
for s, c in zip(scans, COLORS):
    x = np.array(s["x_grid_A_sat_over_mean"])
    R = np.array(s["R"])
    lab = (f"a={s['a_fm']:g} fm, r0={s['r0_fm']:g} fm, L={s['L_um']:g} um")
    ax1.plot(x, R, color=c, lw=2, label=lab)
ax1.axhspan(1, 10, color="#bbbbbb", alpha=0.35, lw=0)
ax1.axhline(10, color=MUTED, lw=1, ls="--")
ax1.text(2e-2, 14, "rescue target  R <= 10  (never reached)", color=MUTED,
         fontsize=9, va="bottom")
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_ylim(1, 1e19)
ax1.set_ylabel(r"R = $\langle A_s^2\rangle/\langle A_s\rangle^2$")
ax1.set_title("F-R9 route (b): on-tube saturation, priced "
              "(h21 web, corpus scale window)", fontsize=11, color=INK)
ax1.legend(loc="lower right", bbox_to_anchor=(1.0, 0.10), fontsize=8.5,
           frameon=False)
# direct label of the flat plateaus
ax1.text(3e13, 4e18, "unsaturated R (h21 window)", fontsize=8.5,
         color=MUTED, ha="center")

# ---- panel B: mean-rate distortion, corpus-central anchor -------------------
s = scans[2]  # a=6.8, L=4.2 um (corpus central)
c = COLORS[2]
x = np.array(s["x_grid_A_sat_over_mean"])
Dc = np.array(s["D_claim"])
Dt = np.array(s["D_true"])
ax2.plot(x, Dc, color=c, lw=2,
         label=r"$D_{\rm claim}=\langle A_s^2\rangle/\bar m^2$  (vs corpus's claimed rate)")
ax2.plot(x, Dt, color=c, lw=2, ls="--",
         label=r"$D_{\rm true}=\langle A_s^2\rangle/\langle A^2\rangle$  (vs unsaturated rate)")
ax2.axhspan(0.5, 1.5, color="#bbbbbb", alpha=0.35, lw=0)
ax2.text(2e-2, 2.2, "|D-1| <= 0.5 band", color=MUTED, fontsize=9)
xc = s["x_at_D_claim_1"]
ax2.axvline(xc, color=MUTED, lw=1, ls=":")
ax2.text(xc * 2.2, 1e-22, f"$D_{{\\rm claim}}=1$ here,\nbut R = "
         f"{s['R_at_D_claim_1']:.1e}\nand $\\langle A_s\\rangle/\\bar m\\approx 2\\times10^{{-8}}$",
         fontsize=8.5, color=MUTED)
xon = s["g0"] / s["phi"]
ax2.axvline(xon, color=MUTED, lw=1, ls=":")
ax2.text(xon * 2.2, 1e-10, "on-tube\namplitude", fontsize=8.5, color=MUTED)
ax2.set_yscale("log")
ax2.set_xscale("log")
ax2.set_ylim(1e-34, 1e19)
ax2.set_xlim(1e-2, 1e18)
ax2.set_xlabel(r"$A_{\rm sat}/\langle A\rangle$")
ax2.set_ylabel("mean-rate distortion D")
ax2.legend(loc="upper left", fontsize=8.5, frameon=False,
           title="corpus-central anchor (a=6.8 fm, L=4.2 um)", title_fontsize=8.5,
           alignment="left")

for ax in (ax1, ax2):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

fig.savefig(BASE + "m3_fig.png", dpi=160, bbox_inches="tight",
            facecolor="white")
print("wrote m3_fig.png")
