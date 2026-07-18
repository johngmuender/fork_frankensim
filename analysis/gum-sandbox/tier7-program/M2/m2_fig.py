#!/usr/bin/env python3
"""M2 figure: lambda-vs-g scaling law + envelope decay + coefficient spread
+ bath FDT check.  Reads m2_results.json, writes m2_fig.png."""

import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, GREEN, MAGENTA = "#2a78d6", "#008300", "#e87ba4"  # validated categorical
INK, INK2, GRID = "#333333", "#666666", "#e3e2de"

with open("m2_results.json") as fh:
    R = json.load(fh)

plt.rcParams.update({
    "font.size": 9.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "text.color": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
    "axes.axisbelow": True, "figure.facecolor": "white"})

fig, axes = plt.subplots(2, 2, figsize=(10.2, 8.0))
(axA, axB), (axC, axD) = axes

# ---- (a) scaling law: lambda vs g, log-log --------------------------------
gs = np.array(R["G2"]["gs"])
lam = np.array(R["G2"]["lam"])
sig = np.array(R["G2"]["sig"])
axA.errorbar(gs, lam, yerr=sig, fmt="o", ms=6, color=BLUE, mec="white",
             mew=0.8, capsize=2.5, lw=1.4, zorder=4,
             label="Ohmic bath, $\\eta=1$ (measured)")
lin = np.array([p["lam_lin_exact"] for p in R["scan_ohmic"]])
axA.plot(gs, lin, "-", color=BLUE, lw=1.2, alpha=0.55, zorder=3,
         label="exact linear response")
gg = np.geomspace(gs[0], gs[-1], 50)
cref = lam[3] / gs[3] ** 2
axA.plot(gg, cref * gg ** 2, "--", color=INK2, lw=1.4, zorder=2,
         label="slope-2 guide $\\propto g_\\chi^2$")
sop = R["superohmic"]["points"]
gso = np.array([p["g"] for p in sop])
lso = np.array([p["lam"] for p in sop])
axA.errorbar(gso, lso, yerr=[p["sig"] for p in sop], fmt="D", ms=6,
             color=MAGENTA, mec="#a04a6e", mew=0.8, capsize=2.5, zorder=4)
axA.plot(gso, [p["lam_lin_exact"] for p in sop], "-", color=MAGENTA, lw=1.2,
         alpha=0.6, zorder=3)
axA.annotate(f"super-Ohmic $J\\propto\\omega^3$, $\\eta_3={R['superohmic']['eta3']:.0f}$\n"
             f"slope $\\approx$ {R['superohmic']['slope']:.2f} — $g^2$ law fails",
             xy=(gso[1], lso[1]), xytext=(0.35, 0.06), textcoords="axes fraction",
             fontsize=8.5, color="#a04a6e",
             arrowprops=dict(arrowstyle="-", color="#a04a6e", lw=0.8))
axA.set_xscale("log"); axA.set_yscale("log")
axA.set_xlabel("$g_\\chi$"); axA.set_ylabel("$\\lambda/\\omega_*$")
axA.set_title(f"(a) G2: $\\lambda\\propto g_\\chi^2$ over one decade — "
              f"slope {R['G2']['slope']:.3f} $\\pm$ {R['G2']['sig_slope']:.3f}",
              fontsize=10)
axA.legend(frameon=False, fontsize=8.5, loc="upper left")

# ---- (b) envelope decay + integrable control ------------------------------
ex = R["example_envelope"]
t = np.array(ex["t"]); q = np.array(ex["Qbar"])
E0 = q[0]
axB.semilogy(t * ex["lam"], q / E0, color=BLUE, lw=1.6,
             label=f"bath on: $\\bar Q(t)$, $g_\\chi={ex['g']:.3f}$")
axB.semilogy(t * ex["lam"], np.exp(-2 * ex["lam"] * t), "--", color=INK2,
             lw=1.3, label="$e^{-2\\lambda t}$ fit")
t5 = np.array(R["G5"]["env_t"]); e5 = np.array(R["G5"]["env_E"])
axB.semilogy(t5 * ex["lam"], e5 / e5[0], color=GREEN, lw=1.6,
             label="bath off (integrable): marginal")
axB.axhline(R["meta"]["T_over_wstar"] / E0, color=INK2, lw=0.9, ls=":",
            label="thermal floor $T$")
axB.set_xlim(0, 2.8)
axB.set_ylim(bottom=0.35 * R["meta"]["T_over_wstar"] / E0)
axB.set_xlabel("$\\lambda t$")
axB.set_ylabel("energy envelope $\\bar Q(t)/\\bar Q(0)$")
axB.set_title("(b) relaxation of $\\delta_0=0.05$: thermalizing vs integrable",
              fontsize=10)
axB.legend(frameon=False, fontsize=8.5, loc="lower left")

# ---- (c) coefficient c = lambda/(g^2 w*) ----------------------------------
ci = np.array(R["G4_ohmic_eta1"]["c_per_g"])
axC.plot(gs, ci, "o-", ms=6, color=BLUE, mec="white", mew=0.8, lw=1.2,
         label="Ohmic $\\eta=1$: $c=\\lambda/(g_\\chi^2\\omega_*)$")
axC.axhline(0.5, color=BLUE, lw=0.9, ls="--", alpha=0.6)
axC.text(gs[0], 0.52, "$\\eta/2$", color=BLUE, fontsize=8.5)
el = R["eta_linearity"]
axC.plot([0.2237], [el["c"]], "s", ms=7, color=GREEN, mec="white", mew=0.8,
         label="Ohmic $\\eta=0.1$")
axC.axhline(0.05, color=GREEN, lw=0.9, ls="--", alpha=0.6)
cso = np.array([p["c"] for p in sop])
axC.plot(gso, cso, "D-", ms=6, color=MAGENTA, mec="#a04a6e", mew=0.8, lw=1.2,
         label=f"super-Ohmic $\\eta_3={R['superohmic']['eta3']:.0f}$ ($g$-dependent)")
axC.axhspan(0.021 - 0.004, 0.021 + 0.004, color="#c9c8c2", alpha=0.55, zorder=1)
axC.text(0.44, 0.0235, "corpus $\\langle r2\\rangle$:\n$0.021\\pm0.004$",
         fontsize=8.5, color=INK)
axC.set_yscale("log"); axC.set_xscale("log")
axC.set_xlabel("$g_\\chi$")
axC.set_ylabel("$c=\\lambda/(g_\\chi^2\\omega_*)$")
axC.set_title("(c) G4: coefficient is a bath-convention number ($c=\\eta/2$)",
              fontsize=10)
axC.legend(frameon=False, fontsize=8, loc="center left")

# ---- (d) bath statistics: FDT ---------------------------------------------
bs = R["G1_bath_stats"]
tg = np.array(bs["tgrid"])
axD.plot(tg, np.array(bs["C_meas"]) * 1e8, "o", ms=5, color=BLUE, mec="white",
         mew=0.7, label="$\\langle\\xi(t)\\xi(0)\\rangle$ sampled, $M=5000$")
axD.plot(tg, np.array(bs["Tgam_cont"]) * 1e8, "-", color=INK2, lw=1.5,
         label="classical FDT: $T\\gamma(t)$ (analytic)")
g0 = bs["gauss"]["t0"]
axD.text(0.35, 0.55,
         f"Gaussianity of $\\xi$:\nskew $= {g0['skew']:+.3f}\\ (\\pm{g0['se_skew']:.3f})$"
         f"\nex. kurt $= {g0['kurtosis_excess']:+.3f}\\ (\\pm{g0['se_kurt']:.3f})$"
         f"\nKS $p = {g0['ks_p']:.2f}$",
         transform=axD.transAxes, fontsize=8.5,
         bbox=dict(fc="white", ec=GRID))
axD.set_xlabel("$t\\,\\omega_*$")
axD.set_ylabel("$\\langle\\xi(t)\\xi(0)\\rangle \\times 10^{8}$")
axD.set_title("(d) G1: bath noise is Gaussian and obeys the FDT", fontsize=10)
axD.legend(frameon=False, fontsize=8.5)

fig.suptitle("M2 — the bootstrap servo coefficient (IV.H.3 $\\langle r2\\rangle$ toy law): "
             "$\\lambda = c\\,g_\\chi^2\\omega_*$ in the thermalizing-bath class",
             fontsize=11.5, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig("m2_fig.png", dpi=160)
print("wrote m2_fig.png")
