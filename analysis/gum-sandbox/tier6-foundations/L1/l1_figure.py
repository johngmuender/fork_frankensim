#!/usr/bin/env python3
"""L1 figure: Bohmian trajectory fan for the flipped Stern-Gerlach exhibit."""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- palette (validated reference instance, light mode) ---------------------
SURFACE = "#fcfcfb"
PAGE = "#f9f9f7"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
BLUE = "#2a78d6"   # slot 1: upper-exit trajectories
GREEN = "#008300"  # slot 2: lower-exit trajectories

d = np.load("L1_traj.npz")
t, Z = d["traj_t"], d["traj_Z"]          # (n_snap,), (n_snap, N_traj)
Zf = d["Zf_o"]
z = d["z"]
rho_up, rho_dn = d["rho_up_o"], d["rho_dn_o"]

rng = np.random.default_rng(7)
sel = rng.choice(Z.shape[1], size=160, replace=False)

fig = plt.figure(figsize=(9.6, 5.4), dpi=200)
fig.patch.set_facecolor(PAGE)
gs = fig.add_gridspec(1, 2, width_ratios=[3.1, 1.0], wspace=0.06,
                      left=0.075, right=0.965, top=0.845, bottom=0.115)
ax = fig.add_subplot(gs[0, 0])
axd = fig.add_subplot(gs[0, 1], sharey=ax)

for a in (ax, axd):
    a.set_facecolor(SURFACE)
    for s in a.spines.values():
        s.set_color(BASELINE)
        s.set_linewidth(0.8)
    a.tick_params(colors=MUTED, labelsize=8.5, length=3)
    for lbl in a.get_xticklabels() + a.get_yticklabels():
        lbl.set_color(INK2)

# ---- trajectory fan ---------------------------------------------------------
ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
ax.axhline(0.0, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
for i in sel:
    c = BLUE if Zf[i] > 0 else GREEN
    ax.plot(t, Z[:, i], color=c, linewidth=0.55, alpha=0.55, zorder=3)
ax.set_xlim(0.0, t[-1])
ax.set_ylim(-9.5, 9.5)
ax.set_xlabel("t", color=INK2, fontsize=10)
ax.set_ylabel("z", color=INK2, fontsize=10)
ax.text(0.03, 0.036, "median z = 0 : never crossed (0 crossings in 4000 trajectories)",
        transform=ax.transAxes, fontsize=8, color=MUTED)

# direct labels on the two branches (identity not by color alone)
box = dict(boxstyle="round,pad=0.45", facecolor=SURFACE, edgecolor=GRID, alpha=0.92)
ax.text(0.86, 6.9, "upper branch\nk = +2.5 labels it “up”\nk = −2.5 labels it “down”",
        fontsize=8.5, color=BLUE, ha="left", va="center", linespacing=1.35,
        bbox=box, zorder=6)
ax.text(0.86, -6.9, "lower branch\nk = +2.5 labels it “down”\nk = −2.5 labels it “up”",
        fontsize=8.5, color=GREEN, ha="left", va="center", linespacing=1.35,
        bbox=box, zorder=6)

# ---- final-density marginal -------------------------------------------------
axd.grid(True, color=GRID, linewidth=0.6, zorder=0)
axd.axhline(0.0, color=MUTED, linewidth=0.9, linestyle=(0, (4, 3)), zorder=2)
axd.fill_betweenx(z, 0, rho_up, color=BLUE, alpha=0.28, linewidth=0, zorder=3)
axd.plot(rho_up, z, color=BLUE, linewidth=1.6, zorder=4)
axd.fill_betweenx(z, 0, rho_dn, color=GREEN, alpha=0.28, linewidth=0, zorder=3)
axd.plot(rho_dn, z, color=GREEN, linewidth=1.6, zorder=4)
axd.set_xlim(0, 1.12 * max(rho_up.max(), rho_dn.max()))
axd.set_xlabel(r"$|\psi_a(z,T)|^2$", color=INK2, fontsize=10)
axd.tick_params(labelleft=False)
axd.text(0.96, 0.86, "up-kicked\npacket", fontsize=8, color=BLUE,
         ha="right", va="center", linespacing=1.3, transform=axd.transAxes)
axd.text(0.96, 0.14, "down-kicked\npacket", fontsize=8, color=GREEN,
         ha="right", va="center", linespacing=1.3, transform=axd.transAxes)

# ---- titles & legend --------------------------------------------------------
fig.text(0.075, 0.955, "Measurements don’t measure: the flipped Stern–Gerlach",
         fontsize=13, color=INK, fontweight="bold")
fig.text(0.075, 0.905,
         "de Broglie–Bohm guidance, 160 of 4000 trajectories shown. Reversing the SG gradient "
         "(k → −k) leaves every trajectory\nbitwise identical (max |ΔZ(T)| = 0) and flips every "
         "spin label — the outcome names the apparatus, not a preexisting spin.",
         fontsize=8.7, color=INK2, va="top", linespacing=1.4)

from matplotlib.lines import Line2D
leg = ax.legend(
    handles=[Line2D([], [], color=BLUE, lw=1.8, label="exit z > 0 (1971/4000)"),
             Line2D([], [], color=GREEN, lw=1.8, label="exit z < 0 (2029/4000)")],
    loc="upper left", fontsize=8, frameon=True, borderpad=0.7, handlelength=1.6)
leg.get_frame().set_facecolor(SURFACE)
leg.get_frame().set_edgecolor(GRID)
for txt in leg.get_texts():
    txt.set_color(INK2)

fig.savefig("L1_flipped_sg.png")
print("wrote L1_flipped_sg.png")
