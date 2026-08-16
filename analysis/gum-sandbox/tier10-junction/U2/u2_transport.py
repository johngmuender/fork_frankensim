#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
U2 — Two-carrier transport sign test (Tier 10, Phase U; F-T10-U2).

Pre-registered by ROADMAP_v12_JUNCTION.md section U2 (frozen spec).
Within-model; nothing here bears on nature.  Sign/monotonicity ONLY —
magnitudes are UNCLAIMABLE (the corpus prints no transport law; context
analysis section 2: transporter SILENT).  The dynamics below are a
[CJ-new] within-model construction on printed scale inputs; every
unprinted step is flagged inline.

THE TOY (as pre-registered):
  Channel V (valence-led): transports B by moving the three original
    knots, mass scale M_knot each, across a rapidity interval.
    Transported inertia  W_V = 3 * M_knot.
  Channel J (junction-led): transports B by migrating the junction+tube
    structure (scales: sigma = 0.19 GeV^2 [IM]; junction mass m_J from
    U1's pre-registered set {0, 0.1355, 0.39}*sqrt(sigma)) and paying
    2x pair-minting cost at destination, at the printed ~100-MeV-class
    snap scale (Primer ex. 15.1 chain: ~1e5 N x 1e-15 m ~ 100 MeV;
    corpus3/03:807,811,827).
    Transported inertia  W_J = m_J + 3*sigma*ell_tube + 2*E_snap.
  Q rides ONLY on knots (V.D winding; K-1 v1-census dictionary).
  Junction carries no charge of any kind (WS-K census rows).

[CJ-new] CONSTRUCTION CHOICES (printed here, none from the corpus):
  (c1) Suppression ansatz: the probability of transporting inertia W
       across a rapidity gap dy is P(W) = g(W), with g ANY normalized
       monotone-decreasing cost function (g(0)=1).  Three families are
       run to demonstrate ansatz-robustness of the sign:
         exponential  g(W) = exp(-W*dy/Lambda)   (Boltzmann-like)
         power-law    g(W) = (1 + W/Lambda)^(-p)
         Gaussian     g(W) = exp(-(W*dy/Lambda)^2)
       Lambda is an arbitrary [CJ-new] scale: NO value is claimable, so
       it is scanned only to show the SIGN never depends on it.
  (c2) Migrating tube-structure inertia: the junction drags three tube
       stubs of length ell_tube; corpus prints no size, so BOTH printed
       hadronic-length readings are scanned:
         compact reading      ell_tube = 1/sqrt(sigma)
         conservative reading ell_tube = 1/f_q   (longest printed
                              hadronic length; makes W_J LARGEST =
                              hardest case for the junction channel)
  (c3) Junction order n = 3 tubes (Primer-asserted 03:811;
       paper-CONDITIONAL on WS1-H3, 01:615/779; execution record absent
       — the F-K0-1 seam is carried, not resolved).
  (c4) Channel weights: equal prior weight, no tuning.  Per unit of
       baryon transported, channel V co-transports one unit of its knot
       charge; channel J transports zero NET charge (minting is local
       and pairwise net-zero at destination; C-blind ensemble).

WHAT IS COMPUTED:
  G1  Validation limits, EXACT in code:
      (a) minting disabled -> junction channel closed -> B-vs-Q
          asymmetry == 0.0 exactly (valence-only);
      (b) tube+junction+minting cost -> 0 -> W_J = 0 -> P_J = g(0) = 1.0
          exactly (its maximum), and in the sharp-suppression limit the
          junction share -> 1.0 exactly (maximal junction dominance).
  G2  The SIGN of the B-vs-Q transport asymmetry and of junction
      dominance over the FULL printed windows:
        M_knot in [1.7, 172.5] GeV;  f_q in [0.14, 0.20] GeV;
        m_J/sqrt(sigma) in {0, 0.1355, 0.39};  E_snap in [0.100, 0.140]
      with the robustness map (margin W_V/W_J) and the analytic
      flip-boundary location.
  G3  The reconciliation ledger: exact rational-arithmetic bookkeeping
      of the junction-led event showing B moves with ZERO original
      knots moving, net signed core count conserved, junction itself
      never carrying B or Q -> the carrier/enforcer/transporter
      trichotomy verdict, as computed.
  G4  Mechanical self-scan: the three STAR anchor values may not appear
      as targets anywhere in the deliverables (regex scan; in
      RESULTS.md they are permitted only on [IM]-anchor-labelled lines).

NO FITTING anywhere.  No STAR observable enters any formula.
"""

import json
import math
import os
import re
import sys
from fractions import Fraction
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# Printed inputs ([IM] / printed-window status as inventoried by the
# context analysis and u_context_agents.json; no other numbers enter).
# ----------------------------------------------------------------------
SIGMA = 0.19                       # GeV^2, corpus tension [IM import]
SQRT_SIGMA = math.sqrt(SIGMA)      # 0.43589 GeV
M_KNOT_WINDOW = (1.7, 172.5)       # GeV: m_Sk floor .. m-tilde_t-class
F_Q_WINDOW = (0.14, 0.20)          # GeV: printed f_q band
MJ_OVER_SQRT_SIGMA = (0.0, 0.1355, 0.39)   # U1 pre-registered set
E_SNAP_WINDOW = (0.100, 0.140)     # GeV: printed 100-MeV-class snap
                                   # scale (Primer ex. 15.1 chain
                                   # ~1e5 N x 1e-15 m ~ 100 MeV; pion
                                   # 140 MeV comparison printed there)
N_TUBE = 3                         # junction order (c3 caveat above)

# Grid resolution (scan, not tuning: results reported over the WHOLE grid)
N_MKNOT, N_FQ, N_SNAP = 241, 61, 9


# ----------------------------------------------------------------------
# Channel inertias (the toy's only "physics"; [CJ-new] construction)
# ----------------------------------------------------------------------
def w_valence(m_knot):
    """Channel V transported inertia: three knots move."""
    return 3.0 * m_knot


def w_junction(m_j, f_q, e_snap, reading):
    """Channel J transported inertia: junction + three tube stubs +
    2x pair-minting cost at destination."""
    if reading == "compact":
        ell = 1.0 / SQRT_SIGMA
    elif reading == "conservative":
        ell = 1.0 / f_q
    else:
        raise ValueError(reading)
    return m_j + N_TUBE * SIGMA * ell + 2.0 * e_snap


# ----------------------------------------------------------------------
# Suppression ansatz families ([CJ-new] c1) — all normalized g(0)=1,
# all strictly monotone decreasing on W >= 0.
# ----------------------------------------------------------------------
def g_exp(w, dy, lam):
    return math.exp(-w * dy / lam)


def g_pow(w, dy, lam, p=2.0):
    # dy folded into lam for the power family (dimensionless form)
    return (1.0 + w * dy / lam) ** (-p)


def g_gauss(w, dy, lam):
    return math.exp(-((w * dy / lam) ** 2))


G_FAMILIES = {"exponential": g_exp, "power_law": g_pow, "gaussian": g_gauss}
# Arbitrary [CJ-new] scale choices, scanned ONLY to show sign-invariance:
LAMBDA_SCAN = (0.5, 1.0, 5.0, 50.0)   # GeV
DY_SCAN = (0.5, 1.0, 3.0)             # rapidity gaps (dimensionless)


def channel_transport(m_knot, m_j, f_q, e_snap, reading, gfun, dy, lam,
                      mint_enabled=True, zero_junction_cost=False):
    """Return (B_transported, Q_transported, P_V, P_J) in toy units.

    Channel V transports 1 unit of B and 1 unit of knot charge per
    success; channel J transports 1 unit of B and 0 net charge, and is
    CLOSED if minting is disabled (a snapped tube must terminate on a
    fresh fractional core: taping rule — no minting, no junction-led
    delivery)."""
    p_v = gfun(w_valence(m_knot), dy, lam)
    if not mint_enabled:
        p_j = 0.0                       # exact closure of the channel
    else:
        w_j = 0.0 if zero_junction_cost else w_junction(m_j, f_q, e_snap,
                                                        reading)
        p_j = gfun(w_j, dy, lam)
    b_t = p_v + p_j
    q_t = p_v
    return b_t, q_t, p_v, p_j


# ----------------------------------------------------------------------
# G1 — validation limits, exact
# ----------------------------------------------------------------------
def validate_limits():
    res = {"checks": [], "all_pass": True}

    def record(name, ok, detail):
        res["checks"].append({"check": name, "pass": bool(ok),
                              "detail": detail})
        if not ok:
            res["all_pass"] = False

    corner = dict(m_knot=M_KNOT_WINDOW[0],
                  m_j=MJ_OVER_SQRT_SIGMA[2] * SQRT_SIGMA,
                  f_q=F_Q_WINDOW[0], e_snap=E_SNAP_WINDOW[1],
                  reading="conservative")

    # (a) minting disabled -> asymmetry EXACTLY zero, all families/scales
    exact_a = True
    for gname, gfun in G_FAMILIES.items():
        for dy in DY_SCAN:
            for lam in LAMBDA_SCAN:
                b, q, _, pj = channel_transport(
                    gfun=gfun, dy=dy, lam=lam, mint_enabled=False, **corner)
                ratio_ok = (q == 0.0) or (b / q == 1.0)  # q may underflow
                if not (pj == 0.0 and (b - q) == 0.0 and ratio_ok):
                    exact_a = False
    record("G1a_minting_disabled_asymmetry_exactly_zero", exact_a,
           "P_J == 0.0, B - Q == 0.0, B/Q == 1.0 exact float equality "
           "for all 3 g-families x %d dy x %d Lambda"
           % (len(DY_SCAN), len(LAMBDA_SCAN)))

    # (b) junction-channel cost -> 0: P_J = g(0) = 1.0 exactly (the
    # maximum of every normalized monotone-decreasing g), hence the
    # asymmetry attains its supremum over W_J at that point ...
    exact_b1 = True
    wj_probe = np.linspace(0.0, 10.0, 101)
    for gname, gfun in G_FAMILIES.items():
        for dy in DY_SCAN:
            for lam in LAMBDA_SCAN:
                b, q, pv, pj = channel_transport(
                    gfun=gfun, dy=dy, lam=lam, zero_junction_cost=True,
                    **corner)
                if pj != 1.0:
                    exact_b1 = False
                # supremum check: g(0) >= g(w) for every probed w
                if any(gfun(w, dy, lam) > pj for w in wj_probe):
                    exact_b1 = False
    record("G1b_zero_cost_PJ_is_exact_maximum", exact_b1,
           "P_J == 1.0 exactly at W_J = 0 and g(0) >= g(w) on a 101-pt "
           "probe, all families/scales")

    # ... and in the sharp-suppression limit the junction share -> 1.0
    # EXACTLY (g(W_V) underflows to 0.0 while g(0) = 1.0): maximal
    # junction dominance.
    b, q, pv, pj = channel_transport(
        gfun=g_exp, dy=1.0, lam=1e-3, zero_junction_cost=True, **corner)
    share = pj / (pv + pj)
    exact_b2 = (share == 1.0 and pv == 0.0 and pj == 1.0)
    record("G1b_sharp_limit_junction_share_exactly_one", exact_b2,
           "exp family, Lambda = 1e-3 GeV: g(W_V) == 0.0 (underflow), "
           "P_J == 1.0, share == 1.0 exact")
    return res


# ----------------------------------------------------------------------
# G2 — sign scan over the full printed windows
# ----------------------------------------------------------------------
def scan_windows():
    m_knot = np.geomspace(*M_KNOT_WINDOW, N_MKNOT)
    f_q = np.linspace(*F_Q_WINDOW, N_FQ)
    e_snap = np.linspace(*E_SNAP_WINDOW, N_SNAP)
    m_j = np.array(MJ_OVER_SQRT_SIGMA) * SQRT_SIGMA

    out = {"grid": {"n_m_knot": N_MKNOT, "n_f_q": N_FQ,
                    "n_e_snap": N_SNAP, "n_m_j": len(m_j),
                    "m_knot_GeV": [float(m_knot[0]), float(m_knot[-1])],
                    "f_q_GeV": list(F_Q_WINDOW),
                    "e_snap_GeV": list(E_SNAP_WINDOW),
                    "m_j_GeV": [float(x) for x in m_j]},
           "readings": {}}

    # margin(W_V / W_J) over the full 4D grid, per tube-length reading
    maps = {}
    for reading in ("compact", "conservative"):
        MK, FQ, MJ, ES = np.meshgrid(m_knot, f_q, m_j, e_snap,
                                     indexing="ij")
        if reading == "compact":
            ell = 1.0 / SQRT_SIGMA
            WJ = MJ + N_TUBE * SIGMA * ell + 2.0 * ES
        else:
            WJ = MJ + N_TUBE * SIGMA / FQ + 2.0 * ES
        WV = 3.0 * MK
        margin = WV / WJ
        i_min = np.unravel_index(np.argmin(margin), margin.shape)
        i_max = np.unravel_index(np.argmax(margin), margin.shape)
        frac_dom = float(np.mean(margin > 1.0))
        out["readings"][reading] = {
            "margin_min": float(margin[i_min]),
            "margin_min_at": {"m_knot_GeV": float(MK[i_min]),
                              "f_q_GeV": float(FQ[i_min]),
                              "m_j_GeV": float(MJ[i_min]),
                              "e_snap_GeV": float(ES[i_min])},
            "margin_max": float(margin[i_max]),
            "W_J_range_GeV": [float(WJ.min()), float(WJ.max())],
            "W_V_range_GeV": [float(WV.min()), float(WV.max())],
            "fraction_of_grid_with_junction_dominance": frac_dom,
            "sign_flip_anywhere_in_window": bool(frac_dom < 1.0),
        }
        # worst-case 2D map over (M_knot, f_q): min margin over m_J, E_snap
        maps[reading] = {"m_knot": m_knot, "f_q": f_q,
                         "margin_worst": margin.min(axis=(2, 3))}

    # analytic flip boundary for the tube-length ansatz (c2 fragility):
    # sign of (W_V - W_J) flips when ell >= (3 M_knot - m_J - 2 E_snap)/(3 sigma)
    mk0 = M_KNOT_WINDOW[0]
    mj_max = MJ_OVER_SQRT_SIGMA[2] * SQRT_SIGMA
    es_max = E_SNAP_WINDOW[1]
    ell_flip = (3.0 * mk0 - mj_max - 2.0 * es_max) / (3.0 * SIGMA)  # GeV^-1
    GEVM1_TO_FM = 0.1973269804
    out["flip_boundary"] = {
        "statement": "at the hardest corner (M_knot = m_Sk floor, m_J max, "
                     "E_snap max) the junction channel stays cheaper until "
                     "the migrating tube-stub length reaches ell_flip",
        "ell_flip_GeV^-1": ell_flip,
        "ell_flip_fm": ell_flip * GEVM1_TO_FM,
        "ell_used_conservative_GeV^-1": 1.0 / F_Q_WINDOW[0],
        "ell_used_conservative_fm": GEVM1_TO_FM / F_Q_WINDOW[0],
        "ell_used_compact_GeV^-1": 1.0 / SQRT_SIGMA,
        "ell_used_compact_fm": GEVM1_TO_FM / SQRT_SIGMA,
        "headroom_factor_conservative": ell_flip * F_Q_WINDOW[0],
        "headroom_factor_compact": ell_flip * SQRT_SIGMA,
    }

    # ansatz-robustness: for EVERY monotone-decreasing g the dominance
    # sign equals sign(W_V - W_J); verified mechanically on a subsampled
    # grid for the three families x Lambda x dy, plus the strict
    # B - Q > 0 sign wherever minting is on.
    sub_mk = m_knot[::24]
    sub_fq = f_q[::12]
    sub_es = e_snap[::4]
    robust = True
    strict_sign = True
    n_pts = 0
    for reading in ("compact", "conservative"):
        for mkv, fqv, mjv, esv in product(sub_mk, sub_fq, m_j, sub_es):
            wv, wj = w_valence(mkv), w_junction(mjv, fqv, esv, reading)
            for gname, gfun in G_FAMILIES.items():
                for dy in DY_SCAN:
                    for lam in LAMBDA_SCAN:
                        pv, pj = gfun(wv, dy, lam), gfun(wj, dy, lam)
                        n_pts += 1
                        # dominance sign == cost sign (monotonicity);
                        # underflow-safe comparison
                        if wj < wv and pj < pv:
                            robust = False
                        if wj > wv and pj > pv:
                            robust = False
                        # B - Q = P_J > 0 whenever the channel is open
                        # (up to float underflow, reported as such)
                        if pj < 0.0:
                            strict_sign = False
    out["ansatz_robustness"] = {
        "n_points_checked": n_pts,
        "dominance_sign_equals_cost_sign_everywhere": bool(robust),
        "B_minus_Q_nonnegative_everywhere": bool(strict_sign),
        "analytic_statement": "for ANY monotone-decreasing g: "
            "W_J < W_V  <=>  g(W_J) >= g(W_V)  (junction dominance), and "
            "B - Q = P_J = g(W_J) > 0 whenever minting is enabled — the "
            "B >= Q sign needs NO cost model at all, only channel-J "
            "openness plus Q-blind minting; the cost model enters only "
            "the dominance statement.",
    }
    return out, maps


# ----------------------------------------------------------------------
# G3 — reconciliation ledger (exact rational arithmetic)
# ----------------------------------------------------------------------
def reconciliation_ledger():
    """Junction-led event, bookkept exactly.

    Initial: baryon = 3 original knots (signed core +1 each) + junction
    (constraint point, no charge, no core) at beam rapidity.
    Event: junction + three tube stubs migrate to destination; each of
    the three tubes snaps AT the destination, minting one knot-antiknot
    pair per tube (locally, pairwise).  The junction constraint (taping
    rule: thirds sum to integers; fractional line ends only on
    fractional core) forces the three destination tube ends to
    terminate on the three minted KNOTS -> signed assembly B = +1 at
    destination.  The minted antiknots stay in the destination
    neighbourhood (they cap the snapped remainders) -> net minted
    charge in the destination region is 0 for EVERY flavor assignment.
    Original knots: never move.
    """
    third = Fraction(1, 3)
    u_charge, d_charge = Fraction(2, 3), Fraction(-1, 3)
    flavors = {"u": u_charge, "d": d_charge}

    results = {"assignments": [], "all_exact": True}
    for combo in product(flavors, repeat=3):
        minted_knot_charges = [flavors[f] for f in combo]
        # per snapped tube: knot(q) + antiknot(-q), both at destination
        minted_pairs_net_q = sum(q + (-q) for q in minted_knot_charges)
        # signed core count: original +3 (beam), minted +3 -3 (dest)
        cores_before = 3
        cores_after = 3 + 3 - 3
        b_before = Fraction(cores_before, 3)
        b_after = Fraction(cores_after, 3)
        b_destination = 3 * third          # three minted knots on junction
        original_knots_moved = 0
        ok = (minted_pairs_net_q == 0 and b_before == b_after == 1
              and b_destination == 1 and original_knots_moved == 0)
        if not ok:
            results["all_exact"] = False
        results["assignments"].append({
            "minted_flavors": "".join(combo),
            "destination_baryon_charge": str(sum(minted_knot_charges)),
            "net_minted_charge_in_destination_region": str(
                minted_pairs_net_q),
            "B_total_before": str(b_before), "B_total_after": str(b_after),
            "B_at_destination": str(b_destination),
            "original_knots_moved": original_knots_moved,
            "exact": ok})

    results["trichotomy_verdict"] = {
        "carrier": "KNOTS — K-1 intact at every step: B is at all times "
                   "(1/3) x signed core count; the junction never carries "
                   "a unit of B (or any charge) in the ledger",
        "enforcer": "the web/taping rule + junction constraint — it is "
                    "exactly this constraint that FORCES the signed "
                    "assembly B = +1 at destination (minted thirds must "
                    "terminate on the junction and sum to an integer)",
        "transporter": "the junction+tube migration mode — B moves across "
                       "the rapidity interval with ZERO original knots "
                       "moving and net signed core count exactly "
                       "conserved (pairs net-zero)",
        "compatibility": "COMPATIBLE, as computed: knot-carried B "
                         "(K-1) and junction-led transport are not in "
                         "tension — carrier, enforcer and transporter are "
                         "three different structural roles, and the "
                         "junction-led mode realizes the third without "
                         "touching the first",
        "status": "[CJ-new] within-model construction: the corpus prints "
                  "the legality of every ingredient (migration-free "
                  "statics, snap+minting existence, reconnection "
                  "net-conservation, taping rule) but no dynamics; this "
                  "ledger tests structural consistency, not a prediction",
    }
    return results


# ----------------------------------------------------------------------
# G4 — mechanical self-scan for the STAR anchor values used as targets
# ----------------------------------------------------------------------
ANCHOR_PATTERNS = [r"(?<![\d.])1\.84(?!\d)", r"(?<![\d.])0\.64(?!\d)",
                   r"(?<![\d.])1\.04(?!\d)"]


def scan_deliverables():
    files = ["u2_transport.py", "u2_results.json", "RESULTS.md"]
    report = {"patterns": ANCHOR_PATTERNS, "files": {}, "pass": True}
    for fn in files:
        path = os.path.join(HERE, fn)
        entry = {"exists": os.path.exists(path), "hits": []}
        if entry["exists"]:
            with open(path, "r", encoding="utf-8") as fh:
                lines = fh.readlines()
            for i, line in enumerate(lines, 1):
                for pat in ANCHOR_PATTERNS:
                    if re.search(pat, line):
                        anchor_ctx = ("[IM]" in line or "anchor" in
                                      line.lower())
                        ok = (fn == "RESULTS.md") and anchor_ctx
                        # redact the anchor digits in the stored report so
                        # u2_results.json itself never contains the values
                        # (keeps the scan idempotent on re-run)
                        red = line.strip()[:160]
                        for p2 in ANCHOR_PATTERNS:
                            red = re.sub(p2, "[ANCHOR-REDACTED]", red)
                        entry["hits"].append(
                            {"line": i, "pattern": pat,
                             "anchor_context": anchor_ctx,
                             "legal": ok, "text": red})
                        if not ok:
                            report["pass"] = False
        report["files"][fn] = entry
    # code and JSON must be completely clean; RESULTS.md hits legal only
    # in [IM]-anchor context
    return report


# ----------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------
def make_figure(maps, scan):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, LogNorm

    surface, ink, ink2 = "#fcfcfb", "#0b0b0b", "#52514e"
    blue, orange = "#2a78d6", "#eb6834"
    seq = LinearSegmentedColormap.from_list(
        "blues1", ["#dbe9fa", "#7fb0e6", "#2a78d6", "#123c73"])

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), dpi=160)
    fig.patch.set_facecolor(surface)

    # Panel A: worst-case margin map, conservative reading
    ax = axes[0]
    ax.set_facecolor(surface)
    mm = maps["conservative"]
    X, Y = np.meshgrid(mm["m_knot"], mm["f_q"], indexing="ij")
    Z = mm["margin_worst"]
    pc = ax.pcolormesh(X, Y, Z, cmap=seq, norm=LogNorm(vmin=1.0,
                       vmax=Z.max()), shading="auto", rasterized=True)
    ax.set_xscale("log")
    cb = fig.colorbar(pc, ax=ax, pad=0.02)
    cb.set_label("margin  $W_V/W_J$  (worst case over $m_J$, $E_{snap}$)",
                 color=ink2, fontsize=9)
    cb.ax.tick_params(colors=ink2, labelsize=8)
    i0 = np.unravel_index(np.argmin(Z), Z.shape)
    ax.plot(X[i0], Y[i0], marker="o", ms=8, mfc="none", mec=ink, mew=1.4)
    ax.annotate("min margin %.3f\n(sign unchanged)" % Z[i0],
                xy=(X[i0], Y[i0]), xytext=(3.2, 0.147), color=ink,
                fontsize=9, arrowprops=dict(arrowstyle="-", color=ink2,
                                            lw=0.8))
    ax.set_xlabel(r"$M_{knot}$  [GeV]  (printed window)", color=ink)
    ax.set_ylabel(r"$f_q$  [GeV]", color=ink)
    ax.set_title("A — junction-dominance margin, conservative reading\n"
                 r"($\ell_{tube}=1/f_q$; flip contour $W_V=W_J$ lies "
                 "outside the window)", fontsize=10, color=ink, loc="left")

    # Panel B: margin vs M_knot at the hardest light-sector corner
    ax = axes[1]
    ax.set_facecolor(surface)
    mk = maps["conservative"]["m_knot"]
    mj_max = MJ_OVER_SQRT_SIGMA[2] * SQRT_SIGMA
    es_max = E_SNAP_WINDOW[1]
    for reading, color, lab in (
            ("conservative", blue,
             r"conservative  $\ell=1/f_q$ ($f_q=0.14$)"),
            ("compact", orange, r"compact  $\ell=1/\sqrt{\sigma}$")):
        fqv = F_Q_WINDOW[0]
        wj = w_junction(mj_max, fqv, es_max, reading)
        ax.plot(mk, 3.0 * mk / wj, color=color, lw=2.0)
    ax.axhline(1.0, color=ink2, lw=1.0, ls="--")
    ax.text(2.0, 1.08, "sign flip boundary  $W_V=W_J$", color=ink2,
            fontsize=8.5)
    ax.text(22.0, 5.6, r"conservative  $\ell=1/f_q$", color=blue,
            fontsize=9)
    ax.text(2.4, 14.0, r"compact  $\ell=1/\sqrt{\sigma}$", color=orange,
            fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$M_{knot}$  [GeV]", color=ink)
    ax.set_ylabel(r"margin  $W_V/W_J$", color=ink)
    ax.set_title("B — margin at hardest corner ($m_J$, $E_{snap}$ max)\n"
                 "junction channel cheaper over the entire printed window",
                 fontsize=10, color=ink, loc="left")

    for ax in axes:
        for s in ax.spines.values():
            s.set_color(ink2)
            s.set_linewidth(0.6)
        ax.tick_params(colors=ink2, labelsize=8)
        ax.grid(True, which="major", color=ink2, alpha=0.15, lw=0.5)

    fig.suptitle("U2 — two-carrier transport sign test: sign/monotonicity "
                 "only, magnitudes unclaimable  [within-model; CJ-new toy]",
                 fontsize=10.5, color=ink, y=1.02)
    fig.tight_layout()
    out = os.path.join(HERE, "u2_fig.png")
    fig.savefig(out, bbox_inches="tight", facecolor=surface)
    return out


# ----------------------------------------------------------------------
def main():
    only_scan = "--scan-only" in sys.argv
    json_path = os.path.join(HERE, "u2_results.json")

    if only_scan and os.path.exists(json_path):
        with open(json_path, "r", encoding="utf-8") as fh:
            results = json.load(fh)
        results["G4_anchor_scan"] = scan_deliverables()
        with open(json_path, "w", encoding="utf-8") as fh:
            json.dump(results, fh, indent=1)
        print("scan-only: G4 scan refreshed; pass =",
              results["G4_anchor_scan"]["pass"])
        return

    results = {
        "workstream": "U2 — two-carrier transport sign test (F-T10-U2)",
        "date": "2026-08-16",
        "register": "within-model; nothing bears on nature; sign/"
                    "monotonicity only; magnitudes UNCLAIMABLE (corpus "
                    "prints no transport law); [CJ-new] construction on "
                    "printed scale inputs; no fitting; STAR values are "
                    "[IM] anchors only and appear nowhere in this file",
        "printed_inputs": {
            "sigma_GeV2": SIGMA, "sqrt_sigma_GeV": SQRT_SIGMA,
            "M_knot_window_GeV": list(M_KNOT_WINDOW),
            "f_q_window_GeV": list(F_Q_WINDOW),
            "m_J_over_sqrt_sigma_set": list(MJ_OVER_SQRT_SIGMA),
            "E_snap_window_GeV": list(E_SNAP_WINDOW),
            "snap_scale_provenance": "Primer ex. 15.1 chain: ~1e5 N x "
                                     "1e-15 m ~ 100 MeV-class "
                                     "(corpus3/03:807,811,827)",
            "junction_order_n": N_TUBE,
            "junction_order_caveat": "Primer-asserted (03:811), "
                                     "paper-conditional on WS1-H3 "
                                     "(01:615,779); F-K0-1 seam carried",
        },
        "construction_CJ_new": {
            "W_V": "3 * M_knot (three knots move; Q rides only on knots "
                   "by V.D winding / K-1)",
            "W_J": "m_J + 3*sigma*ell_tube + 2*E_snap (junction + tube "
                   "stubs migrate; 2x pair-minting paid at destination)",
            "ell_tube_readings": {"compact": "1/sqrt(sigma)",
                                  "conservative": "1/f_q (hardest case)"},
            "suppression_ansatz": "P(W) = g(W), g normalized monotone "
                                  "decreasing; three families run "
                                  "(exponential, power-law, Gaussian); "
                                  "the SIGN is g-independent",
        },
    }

    results["G1_validation_limits"] = validate_limits()
    scan_out, maps = scan_windows()
    results["G2_window_scan"] = scan_out
    results["G3_reconciliation_ledger"] = reconciliation_ledger()

    fig_path = make_figure(maps, None)
    results["figure"] = os.path.basename(fig_path)

    # headline
    cons = scan_out["readings"]["conservative"]
    comp = scan_out["readings"]["compact"]
    results["headline"] = {
        "B_vs_Q_asymmetry_sign": "+1 (B transport >= Q transport) at "
            "every grid point with minting enabled, for every cost "
            "family — the sign requires only channel-J openness plus "
            "Q-blind minting",
        "junction_dominance": "junction channel strictly cheaper (W_J < "
            "W_V) at 100% of the 4D grid in BOTH tube-length readings",
        "min_margin_conservative": cons["margin_min"],
        "min_margin_compact": comp["margin_min"],
        "max_margin": max(cons["margin_max"], comp["margin_max"]),
        "fragility_note": "the conservative-reading margin at the "
            "hardest corner is thin; see flip_boundary — an ell_tube "
            "only ~x%.2f above 1/f_q would flip the corner sign"
            % scan_out["flip_boundary"]["headroom_factor_conservative"],
    }

    results["G4_anchor_scan"] = scan_deliverables()

    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=1)

    print("G1 all_pass:", results["G1_validation_limits"]["all_pass"])
    print("G2 conservative: min margin %.4f at" % cons["margin_min"],
          cons["margin_min_at"], "| dominance fraction",
          cons["fraction_of_grid_with_junction_dominance"])
    print("G2 compact:      min margin %.4f" % comp["margin_min"],
          "| dominance fraction",
          comp["fraction_of_grid_with_junction_dominance"])
    print("G2 robustness:", scan_out["ansatz_robustness"][
        "dominance_sign_equals_cost_sign_everywhere"],
        "on", scan_out["ansatz_robustness"]["n_points_checked"], "points")
    print("flip boundary:", {k: round(v, 4) for k, v in
          scan_out["flip_boundary"].items() if isinstance(v, float)})
    print("G3 ledger all_exact:",
          results["G3_reconciliation_ledger"]["all_exact"])
    print("G4 scan pass:", results["G4_anchor_scan"]["pass"])
    print("figure:", fig_path)
    print("json:", json_path)


if __name__ == "__main__":
    main()
