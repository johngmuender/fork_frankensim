#!/usr/bin/env python3
"""M1 figure: w(z) and W(z)=(w+1)/Om(z) for the tracker family.
Palette: dataviz reference categorical slots 1-4 (validated, light mode);
sub-3:1 slots carry direct labels (relief rule)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from m1_wz import Om, E2, W_memory, run_dq, w_dyn_of, Om_a

SURF, INK, MUT, GRID = "#fcfcfb", "#0b0b0b", "#898781", "#e1e0d9"
C1, C2, C3, C4 = "#2a78d6", "#008300", "#e87ba4", "#eda100"

z = np.linspace(0.0, 3.0, 601)
x = 1.0 + z
a = 1.0 / x
om = Om(x)

w_tr = -1.0 + om
w_m03 = -1.0 + W_memory(x, -0.3) * om
w_m10 = -1.0 + W_memory(x, -1.0) * om
sol = run_dq(200.0)
w_dy = w_dyn_of(sol, 200.0)(a)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.4), dpi=150)
fig.patch.set_facecolor(SURF)
for ax in (ax1, ax2):
    ax.set_facecolor(SURF)
    ax.grid(color=GRID, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=MUT, labelsize=9)

# ---- panel 1: w(z)
ax1.axhline(-1.0, color=MUT, lw=1.0, ls=":")
ax1.plot(z, w_tr, color=C1, lw=2)
ax1.plot(z, w_m03, color=C2, lw=2)
ax1.plot(z, w_m10, color=C3, lw=2)
ax1.plot(z, w_dy, color=C4, lw=2, ls="--")
ax1.text(2.06, -0.115, "pure tracker\nw = −1 + Ωm(z)  [exact]",
         color=C1, fontsize=9, ha="left", va="top")
ax1.text(0.78, -0.175, "β = −0.3", color=C2, fontsize=9)
ax1.text(0.13, -0.30, "β = −1.0", color="#b13a68", fontsize=9)
ax1.text(1.25, -0.44, "ODE, Γ = 200 H₀", color="#a87200",
         fontsize=9, rotation=-8)
ax1.text(2.55, -0.955, "w = −1", color=MUT, fontsize=8.5)
ax1.set_xlabel("z", color=INK, fontsize=10)
ax1.set_ylabel("w(z)", color=INK, fontsize=10)
ax1.set_title("Equation of state: tracker identity, memory forms, dynamics",
              color=INK, fontsize=10.5, loc="left")

# ---- panel 2: W(z) with corpus prints
W03 = W_memory(x, -0.3)
W10 = W_memory(x, -1.0)
W_dyz = (w_dy + 1.0) / om
ax2.axhline(1.0, color=C1, lw=2)
ax2.plot(z, W03, color=C2, lw=2)
ax2.plot(z, W10, color=C3, lw=2)
ax2.plot(z, W_dyz, color=C4, lw=2, ls="--")
pr = [(0.0, 1.28, C2), (2.0, 1.03, C2), (0.0, 1.72, C3), (2.0, 1.05, C3)]
for zz, Wv, cc in pr:
    ax2.scatter([zz], [Wv], s=52, facecolors=SURF, edgecolors=INK,
                zorder=5, linewidths=1.3)
ax2.text(0.09, 1.72, "corpus print 1.72", color=INK, fontsize=8.5, va="center")
ax2.text(0.09, 1.285, "corpus print 1.28", color=INK, fontsize=8.5, va="center")
ax2.text(1.55, 1.115, "1.03 / 1.05 at z = 2", color=INK, fontsize=8.5)
ax2.text(2.35, 0.945, "W = 1 (tracker)", color=C1, fontsize=9)
ax2.text(1.02, 1.21, "β = −0.3", color=C2, fontsize=9)
ax2.text(0.62, 1.50, "β = −1.0", color="#b13a68", fontsize=9)
ax2.text(0.25, 1.035, "ODE Γ = 200 H₀", color="#a87200", fontsize=8.5)
ax2.set_xlabel("z", color=INK, fontsize=10)
ax2.set_ylabel("W(z) = (w+1)/Ωm", color=INK, fontsize=10)
ax2.set_title("W(z): printed pairs land at z = 0 and z = 2 "
              "(recovered convention)", color=INK, fontsize=10.5, loc="left")
ax2.set_ylim(0.90, 1.82)

fig.suptitle("M1 — w(z) from the Γ(H) relaxation profile "
             "(test-field, Ωm0 = 0.3; within-model)",
             color=INK, fontsize=11, x=0.008, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig("/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M1/"
            "m1_wz_fig.png", facecolor=SURF)
print("figure written")
