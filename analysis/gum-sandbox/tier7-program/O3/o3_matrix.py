#!/usr/bin/env python3
"""
O3 -- F-R14 rescue decision matrix [F-T7-O3].

Turns the corpus's one unrepaired structural defect (F-R14: the E.6
positivity premise of Theorem VI.1, FALSE as printed on the corpus's own
action) into a fully computed decision matrix: every textually available
reading priced, nothing recommended.  DECISION SUPPORT ONLY -- the
constitutive choice is the authors' (charter language).

Sources (authoritative inputs, not re-derived here except where checked):
  * theory-audit/h33_RESULTS.md Part B + section 8 (rescue pricing);
    theory-audit/h33_pq.py (master form, factor dictionary, anchors M1-M6).
  * theory-audit/h25_RESULTS.md section 5 (knot rows K2/K3: Dirac-type
    a1 = 4(1/6 - 1/4) = -1/3 per field BEFORE the loop sign; K3
    statistics-signed loop = +2 minimal-scalar equivalents per field),
    section 3 (sector table: N_s = 2, 2, 1 for B2+-, B3, B4).
  * corpus3/01-GUM-Omega-Paper-v4.1-ext.md: F-R14 dispute box (VI.C),
    App. E.6 rewrite, (6.3) margin arithmetic ("the margin is 4.0,
    boundary-exact, IF the bound held").

Implements:
  master form   a1(p,q) = -(5p+q)/(12(2p+q))          [per scalar mode]
  factor dict   sqrt(g3) -> +3/2; derivative pair -> -1;
                field pair: covector -1 / frame 0 / vector +1;
                group-trace kinetic -> 0 (Killing form, no spatial index)
  B2 invariant  p = q - 1 (structural, all variances + density weights)
  Q~ sectors    q = 3/2 group-protected
  knot sector   Dirac-type: -1/3 per field naive (K2) vs +1/3 per field
                supertrace-signed (K3); shown BOTH ways in every row
  P-acoustic    density weight w shifts (p,q) by (3w,3w) off the frame
                baseline; window w in (-5/18, -2/9); distinguished
                w = -1/4 on q = -3p with a1 = +1/6

Sector content for the sums (Sigma w_s = 1/16piG), representative:
  B2+- doublet N=2, B3 doublet N=2, B4 N=1 (printed via the corpus's own
  operator, h25 sec.3 table); knot band-edges: N_knot Dirac species --
  NOT printed anywhere in the corpus -> parametrized (default 3), with
  all thresholds reported as functions of N_knot.  Lambda_s: no
  per-sector cutoffs are printed ("loop weights peak at the cutoff") ->
  equal Lambda_s = Lambda assumed, stated; weights quoted in units of
  Lambda^2.

Exact rational arithmetic (fractions.Fraction) for the matrix; float
scan only for the w-window figure.  Deterministic, no tuning.
Outputs: o3_results.json, o3_fig.png.
"""

import json
import math
import pathlib
from fractions import Fraction as Fr

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent

checks = []


def check(name, passed, detail=""):
    checks.append({"name": name, "pass": bool(passed), "detail": str(detail)})
    print(f"[{'PASS' if passed else 'FAIL'}] {name}  {detail}")


def a1(p, q):
    """h33 master form, per scalar mode. None on the pole ray q = -2p."""
    p, q = Fr(p), Fr(q)
    den = 12 * (2 * p + q)
    if den == 0:
        return None
    return Fr(-(5 * p + q), 1) / den


# ----------------------------------------------------------------------
# 1. Sector (p,q) per reading -- h33's factor dictionary, re-implemented
#    independently in exact arithmetic (exponents are additive).
# ----------------------------------------------------------------------
SQRTG = Fr(3, 2)          # sqrt(g3) measure
DPAIR = Fr(-1)            # each derivative-index contraction g3^{ij}
FPAIR = {"covector": Fr(-1), "frame": Fr(0), "vector": Fr(1)}
# group-trace kinetic (Q~^-1 dt Q~): intrinsic Killing form -> exponent 0


def sector_pq(reading, w=Fr(0)):
    """(p,q) for B2 and for the Q~-sectors (B3/B4) under a reading,
    with an optional density weight w (shifts both by 3w; det g3=(1+h)^3,
    weight enters squared in the quadratic form -> per h33 M4/M6 the
    effective shift on (p,q) is +3w on the frame baseline)."""
    if reading == "cone-only":
        base = {"B2": (Fr(0), Fr(1)), "Q": (Fr(0), Fr(1))}
    else:
        f = FPAIR[reading]
        base = {
            "B2": (SQRTG + DPAIR + f, SQRTG + f),   # (p, q)
            "Q":  (SQRTG + DPAIR + f, SQRTG),       # q = 3/2 group-protected
        }
    return {k: (pv + 3 * w, qv + 3 * w) for k, (pv, qv) in base.items()}


# ----------------------------------------------------------------------
# 2. GATE G1 -- the five h33 anchors, exactly.
# ----------------------------------------------------------------------
anchor_expect = {"covector": Fr(-1, 3), "frame": Fr(-2, 15),
                 "vector": Fr(-5, 33), "cone-only": Fr(-1, 12)}
anchors = {}
for rd, expect in anchor_expect.items():
    pq = sector_pq(rd)["B2"]
    val = a1(*pq)
    anchors[rd] = {"pq_B2": [str(pq[0]), str(pq[1])], "a1_B2": str(val),
                   "expected": str(expect), "match": val == expect}
    check(f"G1 anchor a1_B2({rd}) = {expect}", val == expect,
          f"(p,q)=({pq[0]},{pq[1]}) -> {val}")

# P-acoustic distinguished point: frame baseline + w = -1/4 -> q = -3p, +1/6
wstar = Fr(-1, 4)
pq_star = sector_pq("frame", wstar)["B2"]
a1_star = a1(*pq_star)
on_ray = (pq_star[1] + 3 * pq_star[0] == 0)
anchors["P-acoustic w=-1/4"] = {
    "pq": [str(pq_star[0]), str(pq_star[1])], "on_q_eq_-3p": on_ray,
    "a1": str(a1_star), "expected": "1/6", "match": a1_star == Fr(1, 6)}
check("G1 anchor P-acoustic w=-1/4: on q=-3p, a1 = +1/6",
      on_ray and a1_star == Fr(1, 6),
      f"(p,q)=({pq_star[0]},{pq_star[1]}) -> {a1_star}")

G1_pass = all(v["match"] for v in anchors.values())
check("G1 ALL FIVE ANCHORS EXACT", G1_pass)

# structural invariants (re-verified in exact arithmetic)
inv_ok = all(sector_pq(rd, wv)["B2"][0] == sector_pq(rd, wv)["B2"][1] - 1
             for rd in ("covector", "frame", "vector")
             for wv in (Fr(0), Fr(-1, 4), Fr(1, 7), Fr(-1, 2)))
check("B2 invariant p = q - 1 (all variances, sample density weights)",
      inv_ok)
check("Q~ invariant q = 3/2 group-protected (variance readings, w=0)",
      all(sector_pq(rd)["Q"][1] == Fr(3, 2)
          for rd in ("covector", "frame", "vector")))

# ----------------------------------------------------------------------
# 3. Knot band-edge row (Dirac-type), both bookkeepings (h25 sec. 5).
#    K2 naive:      a1_eff = 4*(1/6 - 1/4) = -1/3 per field (no loop sign)
#    K3 supertrace: Gamma_F = -Tr ln -> +2 minimal-scalar equivalents
#                   = +2*(1/6) = +1/3 per field
# ----------------------------------------------------------------------
A1_KNOT_NAIVE = 4 * (Fr(1, 6) - Fr(1, 4))
A1_KNOT_SUPER = 2 * Fr(1, 6)
check("knot Dirac-type a1 = -1/3 per field before loop sign (h25 K2)",
      A1_KNOT_NAIVE == Fr(-1, 3), str(A1_KNOT_NAIVE))
check("knot supertrace-signed = +2 minimal-scalar equiv = +1/3 (h25 K3)",
      A1_KNOT_SUPER == Fr(1, 3), str(A1_KNOT_SUPER))

# ----------------------------------------------------------------------
# 4. The decision matrix.
# ----------------------------------------------------------------------
N_B2, N_B3, N_B4 = 2, 2, 1          # printed via the corpus operator (h25)
N_KNOT_DEFAULT = 3                   # NOT printed -> parametrized
MARGIN_63 = math.log10(1e-15) - math.log10(1e-19)   # 4.0, boundary-exact

COSTS = {
    "frame": ("favored: consistent with F10' objectivity (M-1 lesson) + "
              "VI.A's KBKK soldering; zero new postulates"),
    "covector": ("violates F10' frame-relative discipline (the corpus's own "
                 "cured defect M-1); abandons the KBKK soldering VI.A "
                 "imports; unmodeled soldering-gradient connection terms "
                 "(h33 sec.9) make its +1/6 the least stable cell"),
    "vector": ("no printed text selects the vector variance; same F10' "
               "tension as covector without the +1/6 payoff"),
    "cone-only": ("reads VI.B's displayed g_s literally (time-time only); "
                  "ignores that spatial covariantization is forced by VI.A "
                  "(the strain-only retreat self-destructs, F-H33-3)"),
    "P-acoustic w=-1/4": ("NEW quantization-measure postulate ([CJ]-class, "
                          "not an erratum), three printed costs: (1) reopens "
                          "Sec. III's flat-measure [DF] derivations on "
                          "curved backgrounds; (2) rho0, J ~ 1/det g3 "
                          "against KBKK wedge-insertion + F5/F8 bookkeeping; "
                          "(3) unforced width-1/18 window excluding both "
                          "natural conventions (tuning by VI.D's own "
                          "standard) -- tuning fraction quantified below"),
}


def row(reading, w=Fr(0), n_knot=N_KNOT_DEFAULT):
    """One fully computed matrix row (both knot bookkeepings)."""
    pq = sector_pq("frame" if reading.startswith("P-acoustic") else reading,
                   w)
    a_B2, a_Q = a1(*pq["B2"]), a1(*pq["Q"])
    sig_bos = N_B2 * a_B2 + (N_B3 + N_B4) * a_Q      # units Lambda^2
    out = {"reading": reading,
           "pq_B2": [str(pq["B2"][0]), str(pq["B2"][1])],
           "pq_Q": [str(pq["Q"][0]), str(pq["Q"][1])],
           "a1_B2": str(a_B2), "a1_B3B4": str(a_Q),
           "sum_bosonic_over_Lambda2": str(sig_bos),
           "N_s": {"B2": N_B2, "B3": N_B3, "B4": N_B4,
                   "knot": f"{n_knot} (parametrized, not printed)"},
           "cost": COSTS[reading]}
    for tag, a_k in (("naive_K2", A1_KNOT_NAIVE),
                     ("supertrace_K3", A1_KNOT_SUPER)):
        total = sig_bos + n_knot * a_k
        signs = {s for s in
                 (int(np.sign(a_B2)), int(np.sign(a_Q)),
                  int(np.sign(a_k)))}
        if signs == {1}:
            hull = ("CONVEX: all w_s > 0; c_GW^2 cannot exit the sector-"
                    "cone hull; (6.2) holds as an inequality")
            m63 = (f"(6.3) restored as inequality; margin = "
                   f"{MARGIN_63:.1f} orders, boundary-exact")
        elif signs == {-1}:
            hull = ("FORMALLY CONVEX (all w_s < 0: common sign divides out "
                    "of the ratio (6.2)) -- but Sigma w_s < 0")
            m63 = (f"(6.2) ratio arithmetic retained (margin {MARGIN_63:.1f} "
                   "orders IF the bound held) but 1/16piG < 0 voids the "
                   "sector upstream (contradicts VI.B)")
        else:
            hull = ("MIXED SIGNS: weights w_s/Sigma w_s not in (0,1); "
                    "c_GW^2 CAN exit the hull; (6.2) fails as an inequality")
            m63 = (f"(6.3) demoted to order-of-magnitude estimate; the "
                   f"verified arithmetic (margin {MARGIN_63:.1f} orders, "
                   "boundary-exact) holds only IF the bound held")
        # supertrace G-flip threshold: total > 0 iff n*|a_k| > |sig_bos|
        thr = None
        if a_k > 0 and sig_bos < 0:
            thr_frac = -sig_bos / a_k          # N_knot must EXCEED this
            thr = {"N_knot_exceeds": str(thr_frac),
                   "min_integer_N_knot": int(thr_frac) + 1,
                   "note": "equal-Lambda assumption"}
        out[tag] = {
            "a1_knot_per_field": str(a_k),
            "w_knot_over_Lambda2": str(n_knot * a_k),
            "sum_ws_over_Lambda2": str(total),
            "G_sign": "+" if total > 0 else ("-" if total < 0 else "0"),
            "G_status": ("G_ind > 0" if total > 0 else
                         "1/16piG = Sigma w_s < 0: INVERTED Newton constant"
                         if total < 0 else "degenerate Sigma w_s = 0"),
            "hull_status": hull, "margin_6p3_status": m63,
        }
        if thr:
            out[tag]["G_flip_threshold"] = thr
    return out


readings = ["frame", "covector", "vector", "cone-only", "P-acoustic w=-1/4"]
matrix = [row(rd, w=(wstar if rd.startswith("P-acoustic") else Fr(0)))
          for rd in readings]

# G2: every cell computed
n_cells = sum(len(r["naive_K2"]) + len(r["supertrace_K3"]) + 6
              for r in matrix)
check("G2 decision matrix fully computed",
      len(matrix) == 5 and all("naive_K2" in r and "supertrace_K3" in r
                               for r in matrix),
      f"5 rows x 2 bookkeepings, ~{n_cells} computed cells")

# G3a: supertrace both ways -- does the G-sign conclusion flip?
flips = {}
for r in matrix:
    gn = r["naive_K2"]["G_sign"]
    gs = r["supertrace_K3"]["G_sign"]
    flips[r["reading"]] = {
        "naive": gn, "supertrace": gs,
        "flips_at_default_N": gn != gs,
        "threshold": r["supertrace_K3"].get("G_flip_threshold")}
check("G3a supertrace-flip analysis computed for every reading",
      len(flips) == 5, str({k: v["flips_at_default_N"]
                            for k, v in flips.items()}))

# ----------------------------------------------------------------------
# 5. P-acoustic pricing: window, tuning fraction, w-scan.
# ----------------------------------------------------------------------
W_LO, W_HI = Fr(-5, 18), Fr(-2, 9)
WIDTH = W_HI - W_LO
check("window width = 1/18", WIDTH == Fr(1, 18), str(WIDTH))

# window-edge verification from the wedge condition on the frame baseline
p_of_w = lambda w: Fr(1, 2) + 3 * w
q_of_w = lambda w: Fr(3, 2) + 3 * w
check("edge w=-2/9 is the 5p+q=0 zero; w=-5/18 the 2p+q=0 pole",
      5 * p_of_w(W_HI) + q_of_w(W_HI) == 0
      and 2 * p_of_w(W_LO) + q_of_w(W_LO) == 0)

NATURAL = {"w=0 (corpus flat measure)": Fr(0),
           "w=-1/2 (half-density)": Fr(-1, 2)}


def set_dist(x):
    if W_LO < x < W_HI:
        return Fr(0)
    return min(abs(x - W_LO), abs(x - W_HI))


dists_set = {k: set_dist(v) for k, v in NATURAL.items()}
center = (W_LO + W_HI) / 2
dists_ctr = {k: abs(v - center) for k, v in NATURAL.items()}
d_set = min(dists_set.values())
d_ctr = min(dists_ctr.values())
tuning_set = WIDTH / d_set          # window width / set-distance
tuning_ctr = WIDTH / d_ctr          # window width / center-distance
check("G3b tuning fraction (set-distance) = (1/18)/(2/9) = 1/4",
      d_set == Fr(2, 9) and tuning_set == Fr(1, 4),
      f"d(w=0)={dists_set['w=0 (corpus flat measure)']}, "
      f"d(w=-1/2)={dists_set['w=-1/2 (half-density)']} (equal)")
check("G3b tuning fraction (center-distance variant) = (1/18)/(1/4) = 2/9",
      center == Fr(-1, 4) and tuning_ctr == Fr(2, 9))

# exact crossover of the NAIVE-knot total inside the window:
# 5*a1(w) = N/3  ->  w_cross = -(20 + 5N/ ... solve exactly with sympy-free
# algebra: a1(w) = -(4+18w)/(30+108w); 5*a1 = N/3
# -> -15(4+18w) = N(30+108w) -> w = -(60+30N)/(270+108N)
N = N_KNOT_DEFAULT
w_cross = Fr(-(60 + 30 * N), (270 + 108 * N))
a1_at_cross = a1(p_of_w(w_cross), q_of_w(w_cross))
check("naive-knot G>0 subwindow edge exact (N_knot=3): w = -25/99",
      w_cross == Fr(-25, 99) and a1_at_cross == Fr(N, 15),
      f"w_cross={w_cross}, a1={a1_at_cross}")
sub_lo, sub_hi = W_LO, w_cross
sub_width = sub_hi - sub_lo
check("naive G>0 subwindow width = 5/198 (of the 1/18 window)",
      sub_width == Fr(5, 198),
      f"fraction of window = {sub_width / WIDTH} = "
      f"{float(sub_width / WIDTH):.4f}")

# ---- float scan across the window (>= 50 interior points) ----
NPTS = 61
lo, hi = float(W_LO), float(W_HI)
ws = lo + (hi - lo) * (np.arange(1, NPTS + 1) / (NPTS + 1))   # open interval
a1w = np.array([-(4 + 18 * w) / (30 + 108 * w) for w in ws])
w_B2 = N_B2 * a1w
w_B3 = N_B3 * a1w
w_B4 = N_B4 * a1w
sig_bos_w = w_B2 + w_B3 + w_B4                       # = 5 a1(w)
knot_naive = N * float(A1_KNOT_NAIVE)
knot_super = N * float(A1_KNOT_SUPER)
tot_naive = sig_bos_w + knot_naive
tot_super = sig_bos_w + knot_super
hull_naive = np.where(a1w > 0, "MIXED (bos+, knot-): can exit",
                      "degenerate")                  # a1w>0 throughout window
hull_super = np.where(a1w > 0, "CONVEX (all +)", "degenerate")
check("scan: a1(w) > 0 at every interior point (window = positivity wedge)",
      bool((a1w > 0).all()),
      f"min a1 = {a1w.min():.4f}, max a1 = {a1w.max():.2f} "
      "(diverges toward the pole edge w=-5/18)")
check("scan: supertrace total > 0 and convex across the ENTIRE window",
      bool((tot_super > 0).all()))
frac_naive_pos = float((tot_naive > 0).mean())
check("scan: naive-knot total positive only on the pole-side subwindow",
      abs(frac_naive_pos - float(sub_width / WIDTH)) < 0.02,
      f"scan fraction {frac_naive_pos:.3f} vs exact {float(sub_width/WIDTH):.3f}"
      " -- and there the hull is still MIXED (knots negative)")

# G3 aggregate
G3_pass = (checks[-1]["pass"] and tuning_set == Fr(1, 4)
           and len(flips) == 5)

# ----------------------------------------------------------------------
# 6. Figure -- weights vs w across the P-acoustic window.
# ----------------------------------------------------------------------
# Okabe-Ito (colorblind-safe) fixed assignment
C_BOS, C_NAIVE, C_SUPER, C_KNOT = "#0072B2", "#D55E00", "#009E73", "#888888"

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.2, 8.4))

# panel A: a1(w) on an extended range, window + conventions in context
wext = np.linspace(-0.55, 0.05, 1200)
a1ext = np.where(np.abs(30 + 108 * wext) > 1e-9,
                 -(4 + 18 * wext) / (30 + 108 * wext), np.nan)
a1ext[np.abs(30 + 108 * wext) < 0.15] = np.nan       # mask pole for plotting
ax1.axvspan(lo, hi, color="#E69F00", alpha=0.25, lw=0,
            label="P-acoustic window (width 1/18)")
ax1.plot(wext, a1ext, color=C_BOS, lw=2,
         label=r"$\mathfrak{a}_1(w)$, frame baseline $+3w$")
ax1.axhline(0, color="#999999", lw=0.8)
ax1.axvline(0, color="#444444", lw=1.2, ls=":",)
ax1.axvline(-0.5, color="#444444", lw=1.2, ls=":")
ax1.axvline(float(W_LO), color="#C05050", lw=1.0, ls="--")
ax1.plot([-0.25], [1 / 6], marker="o", ms=8, color=C_SUPER, zorder=5)
ax1.annotate(r"$w=-\frac{1}{4}$: $\mathfrak{a}_1=+\frac{1}{6}$ (unimodular ray)",
             (-0.25, 1 / 6), xytext=(-0.21, 0.42), fontsize=9,
             arrowprops=dict(arrowstyle="-", lw=0.8, color="#555555"))
ax1.plot([0.0], [-2 / 15], marker="s", ms=7, color="#444444", zorder=5)
ax1.annotate(r"$w=0$ (corpus): $-\frac{2}{15}$", (0.0, -2 / 15),
             xytext=(-0.115, -0.42), fontsize=9,
             arrowprops=dict(arrowstyle="-", lw=0.8, color="#555555"))
ax1.plot([-0.5], [-5 / 24], marker="s", ms=7, color="#444444", zorder=5)
ax1.annotate(r"$w=-\frac{1}{2}$ (half-density): $-\frac{5}{24}$", (-0.5, -5 / 24),
             xytext=(-0.545, -0.46), fontsize=9,
             arrowprops=dict(arrowstyle="-", lw=0.8, color="#555555"))
ax1.annotate("pole edge\n$2p+q=0$\n($w=-5/18$)", (float(W_LO), 0.75),
             fontsize=8, ha="right", color="#C05050",
             xytext=(float(W_LO) - 0.012, 0.55))
ax1.set_ylim(-0.62, 1.05)
ax1.set_xlim(-0.55, 0.05)
ax1.set_xlabel("density weight $w$")
ax1.set_ylabel(r"$\mathfrak{a}_1$ per bosonic mode")
ax1.set_title("A. The rescue window in context: both natural conventions "
              "lie outside", fontsize=10, loc="left")
ax1.legend(loc="upper right", fontsize=8, frameon=False)

# panel B: sector weights and totals inside the window
ax2.axvspan(lo, hi, color="#E69F00", alpha=0.15, lw=0)
ax2.plot(ws, sig_bos_w, color=C_BOS, lw=2,
         label=r"$\Sigma_{\rm bos}\, w_s/\Lambda^2 = 5\,\mathfrak{a}_1(w)$"
               "  (B2$\\pm$+B3+B4)")
ax2.axhline(knot_naive, color=C_KNOT, lw=2, ls="--",
            label=rf"knots naive K2: ${N}\times(-1/3) = -1$")
ax2.axhline(knot_super, color=C_KNOT, lw=2, ls="-.",
            label=rf"knots supertrace K3: $+{N}\times(1/3) = +1$")
ax2.plot(ws, tot_naive, color=C_NAIVE, lw=2.4,
         label=r"$\Sigma w_s$ naive (G-sign: $-$ then $+$; hull MIXED)")
ax2.plot(ws, tot_super, color=C_SUPER, lw=2.4,
         label=r"$\Sigma w_s$ supertrace (G $>0$; hull CONVEX)")
ax2.axhline(0, color="#999999", lw=0.8)
ax2.axvline(-0.25, color=C_SUPER, lw=1.0, ls=":")
ax2.axvline(float(w_cross), color=C_NAIVE, lw=1.0, ls=":")
ax2.annotate(r"$w=-\frac{25}{99}$: naive $\Sigma w_s$ crosses 0",
             (float(w_cross), 0.05), xytext=(float(w_cross) + 0.002, 1.9),
             fontsize=8, color=C_NAIVE)
ax2.annotate(r"$w=-\frac{1}{4}$", (-0.25, 0), xytext=(-0.2495, -1.6),
             fontsize=8, color=C_SUPER)
ax2.set_xlim(lo - 0.002, hi + 0.002)
ax2.set_ylim(-2.1, 6.5)
ax2.set_xlabel("density weight $w$  (window $-5/18 < w < -2/9$, open)")
ax2.set_ylabel(r"weights $/\ \Lambda^2$  ($N_{\rm knot}=3$, equal $\Lambda_s$)")
ax2.set_title("B. Weights across the window: even P-acoustic needs the "
              "supertrace erratum (two repairs minimum)",
              fontsize=10, loc="left")
ax2.legend(loc="upper left", fontsize=8, frameon=False)

fig.suptitle("O3 / F-T7-O3 -- F-R14 rescue decision support: the P-acoustic "
             "family, priced (within-model)", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.97))
figpath = HERE / "o3_fig.png"
fig.savefig(figpath, dpi=160)
print(f"wrote {figpath}")

# ----------------------------------------------------------------------
# 7. Gates + JSON.
# ----------------------------------------------------------------------
G2_pass = checks[[c["name"] for c in checks].index(
    "G2 decision matrix fully computed")]["pass"]
G4_bottom_line = (
    "DECISION SUPPORT ONLY. The constitutive choice is the authors' "
    "(charter language). No reading is recommended here; each is priced. "
    "The F-R14 regrade language is UNCHANGED by this workstream: Theorem "
    "VI.1 remains OPEN-DEFECT (unconditional as printed); (6.3) remains "
    "demoted to an order-of-magnitude estimate; P-M3 remains "
    "premise-suspended pending the E.6 discharge; adoption of P-acoustic "
    "(+ the knot supertrace erratum -- two independent repairs is the "
    "minimum bill) would be a new [CJ]-class constitutive postulate, "
    "not an erratum.")
print("\n" + G4_bottom_line)

n_pass = sum(c["pass"] for c in checks)
print(f"\ninternal checks: {n_pass}/{len(checks)} PASS")

gates = {
    "G1_anchors_exact": {"pass": G1_pass,
                         "measured": {k: v["a1" if "a1" in v else "a1_B2"]
                                      for k, v in anchors.items()}},
    "G2_matrix_fully_computed": {"pass": bool(G2_pass),
                                 "rows": len(matrix),
                                 "bookkeepings_per_row": 2},
    "G3_pricing_and_supertrace": {
        "pass": bool(G3_pass),
        "tuning_fraction_set_distance": str(tuning_set),
        "tuning_fraction_center_distance": str(tuning_ctr),
        "supertrace_flips": flips},
    "G4_honest_bottom_line": {"pass": True, "text": G4_bottom_line},
}

results = {
    "workstream": "O3 -- F-R14 rescue decision matrix [F-T7-O3]",
    "date": "2026-07-18",
    "epistemic_frame": "Within-model; nothing here bears on nature. "
                       "Decision support only; the constitutive choice is "
                       "the authors'.",
    "master_form": "a1(p,q) = -(5p+q)/(12(2p+q)) per scalar mode "
                   "(h33 M1, independently re-implemented in exact "
                   "rational arithmetic here)",
    "anchors": anchors,
    "structural_invariants": {
        "B2": "p = q - 1 (all variances + density weights; re-verified)",
        "Q_sectors": "q = 3/2 group-protected (re-verified)",
        "knot_Dirac": "a1 = 4(1/6 - 1/4) = -1/3 per field before loop "
                      "sign; supertrace-signed = +2*(1/6) = +1/3"},
    "sector_content": {
        "N_s_printed_via_h25_operator_table": {"B2": 2, "B3": 2, "B4": 1},
        "N_knot": f"{N_KNOT_DEFAULT} Dirac species -- NOT printed in the "
                  "corpus; parametrized (thresholds reported in N_knot)",
        "Lambda_s": "no per-sector cutoffs printed; equal Lambda assumed; "
                    "all weights in units Lambda^2"},
    "decision_matrix": matrix,
    "p_acoustic": {
        "window": ["-5/18", "-2/9"], "width": "1/18",
        "edges": {"w=-2/9": "5p+q = 0 (a1 -> 0)",
                  "w=-5/18": "2p+q = 0 POLE (a1 diverges; h25's 'ansatz "
                             "fails structurally' ray is the window edge)"},
        "distinguished_point": {"w": "-1/4", "ray": "q = -3p (unimodular/"
                                "exact-acoustic, Unruh-aligned)",
                                "a1": "1/6"},
        "natural_conventions": {
            k: {"w": str(v), "set_distance_to_window": str(dists_set[k]),
                "center_distance": str(dists_ctr[k])}
            for k, v in NATURAL.items()},
        "tuning_fraction": {
            "definition_primary": "window width / distance from the window "
                                  "(as a set) to the nearest natural "
                                  "convention",
            "value_primary": str(tuning_set),
            "value_center_variant": str(tuning_ctr),
            "note": "both natural conventions are EQUIDISTANT (2/9 set / "
                    "1/4 center) -- the window sits exactly midway"},
        "three_printed_costs": COSTS["P-acoustic w=-1/4"],
        "scan": {
            "n_points": NPTS, "w": ws.tolist(),
            "a1_per_bosonic_mode": a1w.tolist(),
            "sum_bosonic_over_Lambda2": sig_bos_w.tolist(),
            "knot_naive_over_Lambda2": knot_naive,
            "knot_supertrace_over_Lambda2": knot_super,
            "total_naive": tot_naive.tolist(),
            "total_supertrace": tot_super.tolist(),
            "hull_naive": "MIXED at every scanned point (bosonic +, knots "
                          "-): c_GW^2 can exit",
            "hull_supertrace": "CONVEX at every scanned point",
            "naive_G_positive_subwindow": {
                "exact": ["-5/18", str(w_cross)],
                "width": str(sub_width),
                "fraction_of_window": float(sub_width / WIDTH)}},
    },
    "margin_6p3": {"orders": MARGIN_63,
                   "arithmetic": "log10(1e-15) - log10(1e-19) = 4.0, "
                                 "boundary-exact (corpus language retained: "
                                 "'IF the bound held')"},
    "gates": gates,
    "checks": checks, "n_pass": n_pass, "n_total": len(checks),
}
(HERE / "o3_results.json").write_text(json.dumps(results, indent=1))
print(f"wrote {HERE / 'o3_results.json'}")
