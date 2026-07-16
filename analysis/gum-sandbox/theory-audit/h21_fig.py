#!/usr/bin/env python3
"""Phase H2.1 figure: (A) R = <A^2>/<A>^2 vs tube radius a across web spacings,
Campbell-exact curves + direct-MC points; (B) motional-narrowing rescue price R(delta/L).
Reads h21_results.json (produced by h21_mc.py)."""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/user/fork_frankensim/analysis/gum-sandbox/theory-audit"
res = json.load(open(f"{BASE}/h21_results.json"))

# palette (validated reference palette, light mode)
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
C = {1.0: "#2a78d6", 4.2: "#008300", 10.0: "#e87ba4", 1.7: "#eda100"}  # per L_um

grid = res["physical_grid"]
narrow = res["motional_narrowing"]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), facecolor=SURFACE)
for ax in axes:
    ax.set_facecolor(SURFACE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, which="major", color="#e8e7e3", lw=0.7)
    ax.set_axisbelow(True)

# ---- panel A: R vs a, r0 = 2 fm ------------------------------------------------
ax = axes[0]
r0_show = 2.0
for L_um in (1.0, 4.2, 10.0):
    rows = sorted([r for r in grid if r["r0_fm"] == r0_show and r["L_um"] == L_um],
                  key=lambda r: r["a_fm"])
    a = np.array([r["a_fm"] for r in rows])
    R = np.array([r["R_campbell"] for r in rows])
    Raud = np.array([r["R_audit_thick"] for r in rows])
    ax.plot(a, R, "-", color=C[L_um], lw=2, zorder=3, label=f"L = {L_um:g} μm")
    ax.plot(a, Raud, "--", color=C[L_um], lw=1.1, alpha=0.55, zorder=2)
# direct-MC points (any r0/L; plot only r0=2 rows on this panel)
for r in grid:
    if "mc_physical" in r and r["r0_fm"] == r0_show:
        mc = r["mc_physical"]
        ax.errorbar(r["a_fm"], mc["R"], yerr=mc["R_se"], fmt="o", ms=7,
                    mfc=C[r["L_um"]], mec=SURFACE, mew=1.5, ecolor=INK2,
                    capsize=3, zorder=4)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("tube radius a  [fm]", color=INK, fontsize=10)
ax.set_ylabel(r"R = $\langle A^2\rangle\,/\,\langle A\rangle^2$", color=INK, fontsize=10)
ax.set_title("A.  Static web: second moment vs mean  (kernel r0 = 2 fm)\n"
             "solid = Campbell-exact (Poisson web) · dots = direct MC · dashed = audit 1/φ",
             loc="left", fontsize=9.5, color=INK)
ax.axhline(1.0, color=INK2, lw=0.8)
ax.annotate("corpus Step-2 claim: R ≈ 1", xy=(2.2, 1.0), xytext=(0, 5),
            textcoords="offset points", fontsize=8.5, color=INK2)
leg = ax.legend(loc="upper right", frameon=False, fontsize=9, labelcolor=INK)

# ---- panel B: motional narrowing ------------------------------------------------
ax = axes[1]
for entry in narrow:
    key = (entry["a_fm"], entry["r0_fm"], entry["L_um"])
    if key not in [(52.6, 2.0, 4.2), (52.6, 2.0, 1.0)]:
        continue
    x = np.array([c["delta_over_L"] for c in entry["curve"]])
    R = np.array([c["R"] for c in entry["curve"]])
    col = C[entry["L_um"]]
    ls = "-" if entry["L_um"] == 4.2 else (0, (5, 3))
    ax.plot(x, R, ls=ls, color=col, lw=2, label=f"L = {entry['L_um']:g} μm")
ax.legend(loc="upper right", frameon=False, fontsize=9, labelcolor=INK,
          title="curves coincide in δ/L", title_fontsize=8, alignment="right")
ax.axhline(2.0, color=INK2, lw=1.0, ls=":")
ax.annotate("R = 2 (rescue)", xy=(2e-8, 2.0), xytext=(0, 5),
            textcoords="offset points", fontsize=8.5, color=INK2)
e = [n for n in narrow if (n["a_fm"], n["r0_fm"], n["L_um"]) == (52.6, 2.0, 4.2)][0]
d2 = e["delta_over_L_at_R2"]
if d2:
    ax.axvline(d2, color=INK2, lw=0.8, ls=":")
    ax.annotate(f"δ ≈ {d2:.2f} L needed inside one\nevent coherence time τ ≈ r0/c\n→ web speed ≈ "
                f"{e['required_speed_over_c_at_R2']:.1e} c",
                xy=(d2, 30), xytext=(-8, 0), textcoords="offset points",
                ha="right", fontsize=8.5, color=INK)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("web displacement per coherence time  δ / L", color=INK, fontsize=10)
ax.set_ylabel("R (time-averaged kernel)", color=INK, fontsize=10)
ax.set_title("B.  Motional-narrowing rescue price  (a = 52.6 fm, r0 = 2 fm)",
             loc="left", fontsize=9.5, color=INK)

fig.suptitle("H2.1 — NR-K3a second moment: line-supported 0νββ insertion, "
             "R = ⟨A²⟩/⟨A⟩² (T-H2)", fontsize=11, color=INK, x=0.02, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig.savefig(f"{BASE}/h21_fig.png", dpi=160, facecolor=SURFACE)
print("wrote h21_fig.png")
