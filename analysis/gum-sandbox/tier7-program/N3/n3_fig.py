#!/usr/bin/env python3
"""N3 figure: (A) v_x(d, z, t) at n_y = +1 — the field the criterion must make
negative (diverging map, neutral at 0); (B) the mid-channel line velocity
A(t) = v_x(d, 1/2, t) (n_y-independent) with the sampled axial arrivals it
must contradict; (C) M(t) = sup of v_x over the region ladder — never
negative below t ~ 15.9; (D) the v_x discretization bounds vs the signal.
Diverging + sequential color per dataviz color-formula; single-series panels
carry no legend box where the title names the series."""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

BASE = "/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/N3"
INK = "#0b0b0b"
INK2 = "#52514e"
SURF = "#fcfcfb"
BLUES = ["#b7cff0", "#7fa8e0", "#2a78d6", "#12448c", "#0a2a5c"]  # sequential
ACC = "#e87ba4"

raw = np.load(f"{BASE}/n3_raw.npz")
res = json.load(open(f"{BASE}/n3_results.json"))
tg = raw["tgrid"]
A = raw["A_ref"]
zg = raw["zg_line"]
heat = raw["heat_p1"]           # (nt, nz) v_x at n_y = +1, clamped
bNx = raw["b_Nx"]
bdt = raw["b_dt"]
tauax = raw["tau_axial_sampled"]

TAUS = 5.130950
KILL = [(6.4, 6.8), (6.8, 7.2)]
TSTARMID = res["G2"]["tau_star_mid"]

fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.6), dpi=160)
fig.patch.set_facecolor(SURF)
(axA, axB), (axC, axD) = axes

# ---------------- A: heatmap of v_x(d,z,t), n_y = +1 -------------------------
Vc = np.clip(heat, -8, 8)
norm = TwoSlopeNorm(vmin=-8, vcenter=0.0, vmax=8)
im = axA.pcolormesh(tg, zg, Vc.T, norm=norm, cmap="RdBu_r", shading="auto",
                    rasterized=True)
axA.axvline(TAUS, color=INK, lw=1.2, ls="--")
axA.text(TAUS + 0.15, 0.06, "sampled $\\tau_{max}$ = 5.131", fontsize=8,
         color=INK, rotation=90, va="bottom")
for b0, b1 in KILL:
    axA.axvline(b0, color=INK2, lw=0.8, ls=":")
axA.axvline(KILL[1][1], color=INK2, lw=0.8, ls=":")
axA.text(6.42, 0.03, "kill bins", fontsize=8, color="white")
axA.annotate("$v_x>0$ wall branch: $-\\pi\\cot(\\pi z)\\to+\\infty$, $z\\to1$",
             xy=(8.2, 0.94), fontsize=8, color="white")
axA.annotate("mid-channel: $v_x = A(t) > 0$ through the kill bins",
             xy=(5.8, 0.47), fontsize=8, color="#7a1f3d")
axA.set_xlabel("t", color=INK)
axA.set_ylabel("z", color=INK)
axA.set_title("A — $v_x(d{=}1,z,t)$ at $n_y{=}+1$: the criterion needs an all-blue\n"
              "column after $\\tau^*$; red persists at every t (clip $\\pm$8)",
              fontsize=10, color=INK, loc="left")
cb = fig.colorbar(im, ax=axA, pad=0.01, aspect=28)
cb.set_label("$v_x$", color=INK2, fontsize=8)
cb.ax.tick_params(labelsize=7, colors=INK2)

# ---------------- B: A(t) with sampled axial arrivals ------------------------
axB.axhline(0.0, color=INK2, lw=1.0)
for b0, b1 in KILL:
    axB.axvspan(b0, b1, color="#e8e7e3", zorder=0)
axB.plot(tg, np.clip(A, -12, 90), color=BLUES[2], lw=1.6)
late = tauax[tauax > TAUS]
axB.plot(late, np.full(late.size, -10.2), "|", ms=9, color=ACC, alpha=0.8)
axB.text(4.6, -3.9, "sampled axial arrivals ($n_y{=}0$) beyond $\\tau_{max}$\n"
         "(219 of them, to t = 15.79; each requires $A\\geq0$)",
         fontsize=8, color="#a63d68")
axB.axvline(TAUS, color=INK, lw=1.2, ls="--")
axB.axvline(TSTARMID, color=INK, lw=1.2, ls="-.")
axB.text(TSTARMID - 0.45, 2.6, "last $A\\geq0$:\nt = %.2f" % TSTARMID,
         fontsize=8, color=INK, ha="right")
axB.text(6.45, 12, "$A$ = 0.86–0.96 > 0\nin the kill bins", fontsize=8,
         color=INK2)
axB.set_xlim(0, 16)
axB.set_ylim(-12, 90)
axB.set_yscale("symlog", linthresh=2.0)
axB.set_xlabel("t", color=INK)
axB.set_ylabel("$A(t) = v_x(d, z{=}1/2, t)$   [$n_y$-independent]", color=INK)
axB.set_title("B — the mid-channel line velocity every family member shares:\n"
              "positive through the kill bins, recurring $\\geq$0 until 15.87",
              fontsize=10, color=INK, loc="left")

# ---------------- C: M(t) over the region ladder, n_y = +1 -------------------
S = -np.pi / np.tan(np.pi * zg)
sin2 = np.sin(np.pi * zg) ** 2
ladder = [("full line (engine clamp)", np.ones(zg.size, bool), BLUES[0]),
          ("$\\rho$-floor $10^{-3}$", sin2 >= 1e-3, BLUES[1]),
          ("$\\rho$-floor $10^{-2}$", sin2 >= 1e-2, BLUES[2]),
          ("$\\rho$-floor $10^{-1}$", sin2 >= 1e-1, BLUES[3]),
          ("mid-channel only $z{=}1/2$", None, BLUES[4])]
for lab, mask, c in ladder:
    if mask is None:
        M = A.copy()
    else:
        M = np.clip(A[:, None] + 1.0 * S[None, :], -500, 500)[:, mask].max(axis=1)
    axC.plot(tg, M, color=c, lw=1.6, label=lab)
axC.axhline(0.0, color=INK2, lw=1.0)
axC.set_yscale("symlog", linthresh=1.0)
axC.set_xlim(0, 16)
axC.axvline(TAUS, color=INK, lw=1.2, ls="--")
for b0, b1 in KILL:
    axC.axvspan(b0, b1, color="#e8e7e3", zorder=0)
axC.text(0.3, -0.55, "criterion region: $M(t)<0$ — never reached below t = 15.87",
         fontsize=8, color=INK2)
axC.set_xlabel("t", color=INK)
axC.set_ylabel("$M(t) = \\sup_z v_x(d,z,t)$ on region  ($n_y{=}+1$)", color=INK)
axC.set_title("C — the criterion never bites: M(t) $\\geq$ 0 at every t on every\n"
              "$\\rho$-floored region; $\\tau^* = T$ (mid-only: $\\tau^*$ = 15.87)",
              fontsize=10, color=INK, loc="left")
axC.legend(frameon=False, fontsize=8, loc="upper right")

# ---------------- D: discretization bounds vs signal -------------------------
axD.semilogy(tg, np.abs(A) + 1e-18, color="#c9c8c2", lw=1.2,
             label="$|A(t)|$ (signal scale)")
axD.semilogy(tg, bNx + 1e-18, color=BLUES[2], lw=1.2,
             label="$|A_{Nx=2048}-A_{Nx=4096}|$")
axD.semilogy(tg, bdt + 1e-18, color=ACC, lw=1.2,
             label="$|A_{dt}-A_{dt/2}|$ (engine accumulation)")
kb = res["discretization_bounds"]["total_pointwise"]["kill_bins"]
axD.text(0.3, 3e-16, "kill-bin bound: median %.1e, p95 %.1e\n"
         "vs required sign flip of $A\\sim0.9$: deficit $\\sim10^{4}\\times$"
         % (kb["median"], kb["p95"]), fontsize=8, color=INK2)
axD.set_xlim(0, 16)
axD.set_ylim(1e-17, 3e2)
axD.set_xlabel("t", color=INK)
axD.set_ylabel("magnitude", color=INK)
axD.set_title("D — the failure is not discretization: the bound on $v_x$ at the\n"
              "line sits 4–13 orders below the wrong-signed field",
              fontsize=10, color=INK, loc="left")
axD.legend(frameon=False, fontsize=8, loc="upper right")

for ax in (axA, axB, axC, axD):
    ax.set_facecolor(SURF)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c9c8c2")
    ax.tick_params(colors=INK2, labelsize=8)
    if ax is not axA:
        ax.grid(True, which="major", color="#e8e7e3", lw=0.7)
        ax.set_axisbelow(True)

fig.suptitle("N3 — field-level crossing criterion at d = 1: hypothesis refuted, "
             "the L7b cutoff is a transport phenomenon (F-T7-N3)",
             fontsize=11.5, color=INK, y=0.995)
fig.tight_layout(rect=(0, 0, 1, 0.965))
out = f"{BASE}/n3_fig.png"
fig.savefig(out, facecolor=fig.get_facecolor())
print("wrote", out)
