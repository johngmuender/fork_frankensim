#!/usr/bin/env python3
"""O1 figure: (A) the pre-crossing epoch s_pc(z, t) of the backward paths at
n_y = +1 — every kill-bin detector point traces back across x = d well before
its start time; (B) sample backward trajectories X(s) showing the pre-crossing
above x = d and the return to the initial support; (C) the class-(a)
pre-crossing margin distribution vs the per-path refinement error — the G3
separation; (D) the discrimination control n_y = +0.25: the same backward map
classifies a large fresh region whose detector-cell flux reproduces L7b's
measured bin probabilities."""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

BASE = "/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/O1"
INK = "#0b0b0b"
INK2 = "#52514e"
SURF = "#fcfcfb"
BLUES = ["#b7cff0", "#7fa8e0", "#2a78d6", "#12448c", "#0a2a5c"]
ACC = "#e87ba4"
GRN = "#3e8f6b"

raw = np.load(f"{BASE}/o1_raw.npz")
res = json.load(open(f"{BASE}/o1_results.json"))

ny_flat = raw["ny_flat"]
t_flat = raw["t_flat"]
z_flat = raw["z_flat"]
zc = raw["zc"]
tc = raw["tc"]
cls = raw["cls_ref"]
margin = raw["margin_ref"]
s_pc = raw["s_pc_ref"]
NZ = zc.size
NT = tc.size
KILL = [(6.4, 6.8), (6.8, 7.2)]
D = 1.0

def fam(ny):
    return np.isclose(ny_flat, ny)

fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6), dpi=160)
fig.patch.set_facecolor(SURF)
(axA, axB), (axC, axD) = axes

# ---------------- A: pre-crossing epoch map, n_y = +1 ------------------------
m1 = fam(1.0)
spc_map = s_pc[m1].reshape(NT, NZ)
cls_map = cls[m1].reshape(NT, NZ)
im = axA.pcolormesh(tc, zc, spc_map.T, cmap="viridis", shading="auto",
                    rasterized=True)
bad = cls_map.T != 0
if bad.any():
    jj, ii = np.where(bad)
    axA.plot(tc[ii], zc[jj], "x", color=ACC, ms=4)
axA.axvline(6.8, color="white", lw=0.8, ls=":")
ep = res["G3"]["crossing_epochs"]["+1.00"]
axA.set_xlabel("detector time t", color=INK)
axA.set_ylabel("z", color=INK)
axA.set_title("A — certified pre-crossing epoch $s_{pc}(z,t)$ at $n_y{=}+1$: every\n"
              "kill-bin detector point back-traces across $x{=}d$ at $s_{pc} < t$"
              " (median %.2f)" % ep["median"],
              fontsize=10, color=INK, loc="left")
cb = fig.colorbar(im, ax=axA, pad=0.01, aspect=28)
cb.set_label("$s_{pc}$", color=INK2, fontsize=8)
cb.ax.tick_params(labelsize=7, colors=INK2)

# ---------------- B: sample backward trajectories, n_y = +1 ------------------
hist_idx = raw["hist_idx"]
histS = raw["histS"]
histX = raw["histX"]
sel = np.where(np.isclose(ny_flat[hist_idx], 1.0))[0]
axB.axhline(D, color=INK, lw=1.2, ls="--")
axB.text(6.9, 1.12, "$x = d$", fontsize=8, color=INK, ha="right")
axB.axhspan(-8.5, -1.5, color="#e8e7e3", zorder=0)
axB.text(0.15, -5.3, "supp $\\phi_0$", fontsize=8, color=INK2)
n_shown = 0
for k, j in enumerate(sel):
    x = histX[:, j]
    ok = np.isfinite(x)
    xs, ss = x[ok], histS[ok]
    # trim the frozen tail (history keeps the last value after freezing)
    chg = np.where(np.abs(np.diff(xs)) > 0)[0]
    if chg.size:
        last = chg[-1] + 2
        xs, ss = xs[:last], ss[:last]
    axB.plot(ss, xs, lw=1.0, color=BLUES[k % len(BLUES)], alpha=0.85)
    n_shown += 1
axB.set_xlim(0, 7.3)
axB.set_ylim(-9.5, 3.6)
axB.set_xlabel("backward time s", color=INK)
axB.set_ylabel("$X_x(s)$", color=INK)
axB.set_title("B — backward paths from kill-bin points ($n_y{=}+1$, %d shown):\n"
              "each exceeds $x=d$ earlier (pre-crossing); frozen at margin cap 2"
              % n_shown,
              fontsize=10, color=INK, loc="left")

# ---------------- C: margin distribution vs refinement error -----------------
kill_sel = np.isin(np.round(ny_flat * 100).astype(int),
                   [100, -100, 75, -75]) & (cls == 0)
mar = margin[kill_sel]
bins = np.logspace(-4, 0.45, 60)
axC.hist(np.clip(mar, 1e-4, None), bins=bins, color=BLUES[2],
         edgecolor="none")
axC.set_xscale("log")
axC.set_yscale("log")
g3 = res["G3"]
axC.axvline(g3["min_margin_class_a_kill"], color=INK, lw=1.2, ls="--")
axC.text(g3["min_margin_class_a_kill"] * 1.15, 2e3,
         "min margin\n%.3f" % g3["min_margin_class_a_kill"],
         fontsize=8, color=INK)
axC.axvline(g3["margin_err_p95"], color=ACC, lw=1.2, ls="-.")
axC.text(g3["margin_err_p95"] * 1.15, 3e2,
         "p95 refinement\npath error %.1e" % g3["margin_err_p95"],
         fontsize=8, color="#a63d68")
axC.text(2e-3, 30.0, "separation $\\times$%.0f" %
         g3["ratio_min_margin_over_p95err"], fontsize=9, color=INK2)
axC.set_ylim(bottom=10.0)
axC.set_xlabel("class-(a) pre-crossing margin  max$_s X_x - d$  (capped at 2)",
               color=INK)
axC.set_ylabel("detector points", color=INK)
axC.set_title("C — G3: the smallest certified margin sits far above the\n"
              "backward-trajectory discretization error",
              fontsize=10, color=INK, loc="left")

# ---------------- D: control n_y = +0.25 classification map ------------------
mc = fam(0.25)
cmap_cls = ListedColormap([BLUES[1], "#d9b13b", GRN, "#999999", "#777777",
                           ACC])
cm = cls[mc].reshape(NT, NZ)
imD = axD.pcolormesh(tc, zc, cm.T, cmap=cmap_cls, vmin=-0.5, vmax=5.5,
                     shading="auto", rasterized=True)
axD.axvline(6.8, color="white", lw=0.8, ls=":")
g5 = res["G5"]
axD.set_xlabel("detector time t", color=INK)
axD.set_ylabel("z", color=INK)
axD.set_title("D — discrimination ($n_y{=}+0.25$, populated bins): fresh region\n"
              "(green, %d pts) carries flux %.4f vs L7b measured $\\Pi$ = %.4f"
              % (g5["n_fresh"], g5["flux_fresh_total"], g5["l7b_pi_total"]),
              fontsize=10, color=INK, loc="left")
axD.text(6.42, 0.04, "blue = pre-crossed", fontsize=8, color="white")
axD.text(6.42, 0.90, "green = legitimately fresh (first crossings live here)",
         fontsize=8, color="white")

for ax in (axA, axB, axC, axD):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c9c8c2")
    ax.tick_params(colors=INK2, labelsize=8)
    if ax in (axB, axC):
        ax.grid(True, which="major", color="#e8e7e3", lw=0.7)
        ax.set_axisbelow(True)

fig.suptitle("O1 — backward-reachability proof of the kill-bin zeros: every "
             "detector point pre-crossed, zero violations at three resolutions "
             "(F-T7-O1)", fontsize=11.5, color=INK, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.965))
out = f"{BASE}/o1_fig.png"
fig.savefig(out, facecolor=fig.get_facecolor())
print("wrote", out)
