#!/usr/bin/env python3
"""
Tier-7 N2 — F-R9 route (a): the honest second-moment treatment (0nubb insert).

Question. h21 (theory-audit, F-R9) established that the corpus's NR-K3a Step-2 claim
(0nubb rate governed by the volume mean m_bar; "fluctuations negligible") drops the
dominant term: the observable is quadratic in the amplitude and the line-supported web
gives R = <A^2>/<A>^2 = 1e14-1e18 across the corpus's own scale window (mean claim
simultaneously verified at 1.000 +/- 0.001). M3 closed rescue route (b) (on-tube
saturation: NOT-VIABLE). This workstream executes route (a): not a rescue — the honest
recomputation of the observable, and its propagation into the corpus's Majoron/P-nu4
phenomenology.

Physics computed here:
 1. SPATIAL SELF-AVERAGING. For a macroscopic sample of N nuclei the total rate is
    T = sum_i A(x_i)^2 -> N <A^2> with computed relative fluctuation: the web is static
    on experiment timescales (h21's motional-narrowing pricing: narrowing would need web
    transport at 1.5e8-1.2e9 c — causally excluded; campaign ledger prints 5e7-1.2e9 c),
    so each nucleus samples a frozen snapshot and the ensemble average is a SPATIAL
    average. RSD^2(T) = (K4 - 1)/N + Var(Lambda)/<Lambda>^2, with K4 = <A^4>/<A^2>^2
    (point-sampling term, Campbell: <A^4> ~= rho_L*J4 in the sparse limit) and the
    web-realization term from the total tube length Lambda in the sample volume
    (Poisson line process in a ball of radius R_b: Var/mean^2 = (9/8)/(rho_L pi R_b^2)).
 2. BURSTINESS / LORENZ STRUCTURE. Fraction of the total rate carried by the top 10^-k
    fraction of sites (k = 2..16), by quadrature AND importance-sampled MC with
    jackknife errors; plus the site-fraction quantiles q50/q90/q99 and the on-tube
    (d <= a) share — the 'hot-site' picture.
 3. PROPAGATION. Majoron-mode 0nubb rate ~ g^2 (corpus P-nu4: g ~ 1e-9, far-horizon;
    printed window g in [8e-10, 1.3e-8]); honest rate = R x printed rate <=>
    g_equiv = sqrt(R) x g. Table per anchor; retreat arithmetic against the corpus's
    own printed epsilon-bracket; pre-registered one-of-three verdict (no tuning).

Parametrization: REUSED from theory-audit/h21_mc.py verbatim (gate condition).
  * Static isotropic Poisson line process, tube radius a, uniform pitch amplitude
    m_hat = 1 inside, 0 outside; rho_L = 1/L_web^2 (phi = pi a^2/L^2); m_bar = phi.
  * Kernel w(r) = exp(-r/r0)/(4 pi r0 r^2), r0 in {1, 2} fm (h21's acknowledged choice).
  * Anchors (as M3): (a, r0, L) = (52.6, 2, 1.0), (52.6, 2, 4.2), (6.8, 2, 4.2),
    (2.0, 1, 4.2) [fm, fm, um] — spanning h21's R = 1.1e14 .. 8.9e17.

Deterministic seed. Output: n2_results.json.
"""

import json
import time
import numpy as np

RNG_SEED = 20260718
FM_PER_UM = 1.0e9

T_START = time.time()

# ============================================================================
# Section 1 — REUSED h21 machinery (copied verbatim from theory-audit/h21_mc.py)
# ============================================================================

def make_h_table(r0, s_min_frac=1e-9, s_max_over_r0=120.0, n_s=1400, n_t=3000):
    """Line-integrated kernel  h(s) = int_-inf^inf w(sqrt(s^2+z^2)) dz."""
    s = np.geomspace(s_min_frac * r0, s_max_over_r0 * r0, n_s)
    t = np.linspace(0.0, 32.0, n_t)
    ch = np.cosh(t)
    x = s[:, None] / r0
    expo = -x * ch[None, :]
    integ = np.where(expo > -745.0, np.exp(np.maximum(expo, -745.0)), 0.0) / ch[None, :]
    I = np.trapezoid(integ, t, axis=1)
    h = I / (2.0 * np.pi * r0 * s)
    return s, h


def h_interp(s, s_grid, h_grid):
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    inside = (s > s_grid[0]) & (s < s_grid[-1])
    out[inside] = np.exp(np.interp(np.log(s[inside]), np.log(s_grid),
                                   np.log(np.maximum(h_grid, 1e-300))))
    out[s <= s_grid[0]] = h_grid[0]
    return out


def g_of_d(d, a, r0, s_grid, h_grid, n_nodes=900):
    """g(d) = int over disc(center distance d, radius a) of h(s) d^2s."""
    s_cut = s_grid[-1]
    s_lo = max(d - a, 0.0)
    s_hi = min(d + a, s_cut)
    if s_lo >= s_hi:
        return 0.0
    nodes = [np.geomspace(max(s_lo, 1e-10 * r0), s_hi, n_nodes)]
    if s_lo == 0.0:
        nodes.append(np.geomspace(1e-10 * r0, min(1e-2 * r0, s_hi), 120))
    kink = abs(d - a)
    if s_lo < kink < s_hi:
        nodes.append(np.linspace(max(s_lo, kink - 6 * r0), min(s_hi, kink + 6 * r0), 300))
    nodes.append(np.linspace(max(s_lo, 1e-10 * r0), s_hi, 400))
    s = np.unique(np.concatenate(nodes))
    s = s[(s > 0) & (s <= s_hi)]
    h = h_interp(s, s_grid, h_grid)
    if d < 1e-12:
        psi = np.where(s <= a, 2 * np.pi, 0.0)
    else:
        psi = np.zeros_like(s)
        full = s <= (a - d)
        psi[full] = 2 * np.pi
        part = (~full) & (s > abs(d - a)) & (s < d + a)
        arg = (s[part] ** 2 + d ** 2 - a ** 2) / (2 * s[part] * d)
        psi[part] = 2 * np.arccos(np.clip(arg, -1.0, 1.0))
    return float(np.trapezoid(h * s * psi, s))


def build_g_table(a, r0, s_grid, h_grid, extent_r0=100.0, n_extra=500):
    d_max = a + extent_r0 * r0
    grid = [np.array([0.0]),
            np.geomspace(1e-3 * min(a, r0), d_max, 500),
            np.linspace(max(0.0, a - extent_r0 * r0), d_max, n_extra),
            np.linspace(max(0.0, a - 8 * r0), min(d_max, a + 8 * r0), 300)]
    d = np.unique(np.clip(np.concatenate(grid), 0.0, d_max))
    g = np.array([g_of_d(x, a, r0, s_grid, h_grid) for x in d])
    return d, g, d_max


def campbell_moments(d, g, rho_L):
    mean = rho_L * np.trapezoid(2 * np.pi * d * g, d)
    var = rho_L * np.trapezoid(2 * np.pi * d * g ** 2, d)
    return mean, var


def physical_mc_R(a, r0, L_web, d_grid, g_grid, d_max, rng, n=2_000_000):
    """h21 leg 3b verbatim: importance-sampled direct MC of R at physical scales."""
    rho_L = 1.0 / L_web ** 2
    d_lo = 1e-3 * min(a, r0)
    lnZ = np.log(d_max / d_lo)
    d = d_lo * np.exp(rng.random(n) * lnZ)          # q(d) = 1/(d lnZ)
    gg = np.interp(d, d_grid, g_grid)
    w1 = 2 * np.pi * d * gg * (d * lnZ)             # integrand / proposal
    w2 = 2 * np.pi * d * gg ** 2 * (d * lnZ)
    I1 = w1.mean()
    I2 = w2.mean()
    mean = rho_L * I1
    var = rho_L * I2
    R = 1.0 + var / mean ** 2
    k = 20
    slices = np.array_split(np.arange(n), k)
    Rj = np.empty(k)
    S1, S2 = w1.sum(), w2.sum()
    for i, s in enumerate(slices):
        n_i = n - s.size
        I1_i = (S1 - w1[s].sum()) / n_i
        I2_i = (S2 - w2[s].sum()) / n_i
        Rj[i] = 1.0 + (rho_L * I2_i) / (rho_L * I1_i) ** 2
    R_se = np.sqrt((k - 1) / k * np.sum((Rj - Rj.mean()) ** 2))
    return {"R": float(R), "R_se": float(R_se), "mean": float(mean),
            "mean_over_phi": float(mean / (rho_L * np.pi * a ** 2))}


# ============================================================================
# Section 2 — NEW: burstiness (Lorenz structure) by quadrature and MC
# ============================================================================

def lorenz_quadrature(d_grid, g_grid, rho_L, a):
    """Cumulative rate share C(d) = int_0^d 2 pi t g^2 dt / J2 on the g-table grid,
    site fraction q(d) = pi rho_L d^2 (sparse limit: nearest-tube distance measure).
    Requires g monotone nonincreasing (checked; report worst violation)."""
    integ = 2 * np.pi * d_grid * g_grid ** 2
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1])
                                           * np.diff(d_grid))])
    J2 = cum[-1]
    C = cum / J2
    q = np.pi * rho_L * d_grid ** 2
    # monotonicity of g (allow numerical wiggle)
    dg = np.diff(g_grid)
    worst = float(dg.max() / max(g_grid[0], 1e-300))
    return C, q, J2, worst


def share_at_site_fraction(C, q, q_target):
    """Rate share carried by the top q_target fraction of sites (in-model exact 1.0
    when q_target exceeds the kernel-reach fraction q[-1])."""
    if q_target >= q[-1]:
        return 1.0
    return float(np.interp(q_target, q, C))


def site_fraction_at_share(C, q, share_target):
    """Site fraction q needed to carry share_target of the rate."""
    return float(np.interp(share_target, C, q))


def lorenz_mc(a, r0, L_web, d_grid, g_grid, d_max, rng, q_targets, a_share=True,
              n=2_000_000, n_jk=20):
    """Importance-sampled MC (h21 leg-3b sampling) of the Lorenz shares with jackknife
    errors. Sites at nearest-tube distance d with measure 2 pi rho_L d dd; rate weight
    g(d)^2. Share(top q) = sum(w2[d <= d_q]) / sum(w2), d_q = sqrt(q/(pi rho_L))."""
    rho_L = 1.0 / L_web ** 2
    d_lo = 1e-3 * min(a, r0)
    lnZ = np.log(d_max / d_lo)
    d = d_lo * np.exp(rng.random(n) * lnZ)
    gg = np.interp(d, d_grid, g_grid)
    w2 = 2 * np.pi * d * gg ** 2 * (d * lnZ)
    thresholds = {}
    for q in q_targets:
        d_q = np.sqrt(q / (np.pi * rho_L))
        thresholds[q] = d_q
    if a_share:
        thresholds["on_tube"] = a
    slices = np.array_split(np.arange(n), n_jk)
    out = {}
    S_tot = w2.sum()
    for key, d_q in thresholds.items():
        if d_q >= d_max:
            out[key] = {"share": 1.0, "se": 0.0, "exact_in_model": True,
                        "d_q_fm": float(d_q)}
            continue
        mask = d <= d_q
        S_in = w2[mask].sum()
        share = S_in / S_tot
        jk = np.empty(n_jk)
        for i, s in enumerate(slices):
            m_s = mask[s]
            jk[i] = (S_in - w2[s][m_s].sum()) / (S_tot - w2[s].sum())
        se = np.sqrt((n_jk - 1) / n_jk * np.sum((jk - jk.mean()) ** 2))
        out[key] = {"share": float(share), "se": float(se),
                    "exact_in_model": False, "d_q_fm": float(d_q)}
    # MC quantiles q50/q90/q99 with jackknife
    order = np.argsort(d)
    d_sorted = d[order]
    w_sorted = w2[order]
    cw = np.cumsum(w_sorted)
    quantiles = {}
    for lvl in (0.5, 0.9, 0.99):
        i_c = int(np.searchsorted(cw, lvl * cw[-1]))
        d_lvl = d_sorted[min(i_c, n - 1)]
        # jackknife on d_lvl -> q_lvl
        jk = np.empty(n_jk)
        slice_of = np.empty(n, dtype=int)
        for i, s in enumerate(slices):
            slice_of[s] = i
        so_sorted = slice_of[order]
        for i in range(n_jk):
            keep = so_sorted != i
            cwi = np.cumsum(w_sorted[keep])
            ii = int(np.searchsorted(cwi, lvl * cwi[-1]))
            jk[i] = d_sorted[keep][min(ii, len(cwi) - 1)]
        d_se = np.sqrt((n_jk - 1) / n_jk * np.sum((jk - jk.mean()) ** 2))
        q_lvl = np.pi * rho_L * d_lvl ** 2
        q_se = 2 * np.pi * rho_L * d_lvl * d_se
        quantiles[f"q{int(lvl*100)}"] = {"d_fm": float(d_lvl),
                                         "site_fraction": float(q_lvl),
                                         "site_fraction_se": float(q_se)}
    return out, quantiles


# ============================================================================
# Section 3 — self-averaging (Campbell-quadrature; the web is static, h21 sec 6)
# ============================================================================

def self_averaging(d_grid, g_grid, rho_L, N_nuclei, D_samples_m, L_web_fm):
    """RSD of the sample-total rate T = sum_i A_i^2 over N nuclei in a ball of
    diameter D: RSD^2 = (K4 - 1)/N  +  (9/8)/(rho_L pi R_b^2).
    K4 = <A^4>/<A^2>^2 with <A^4> ~= rho_L J4 (sparse limit; subleading + 3(rho_L J2)^2
    included), <A^2> ~= rho_L J2 (variance-dominated).  Lengths in fm."""
    J2 = np.trapezoid(2 * np.pi * d_grid * g_grid ** 2, d_grid)
    J4 = np.trapezoid(2 * np.pi * d_grid * g_grid ** 4, d_grid)
    m2 = rho_L * J2                       # <A^2> (variance part dominates mean^2)
    m4 = rho_L * J4 + 3.0 * (rho_L * J2) ** 2
    K4 = m4 / m2 ** 2
    rows = []
    for D_m in D_samples_m:
        R_b_fm = (D_m / 2.0) * 1.0e15
        web_term = (9.0 / 8.0) / (rho_L * np.pi * R_b_fm ** 2)
        point_term = (K4 - 1.0) / N_nuclei
        rows.append({"D_sample_m": D_m,
                     "web_realization_RSD": float(np.sqrt(web_term)),
                     "point_sampling_RSD": float(np.sqrt(point_term)),
                     "total_RSD": float(np.sqrt(web_term + point_term))})
    return {"K4": float(K4), "J2": float(J2), "J4": float(J4),
            "N_nuclei": N_nuclei, "samples": rows}


# ============================================================================
# main
# ============================================================================

ANCHORS = [
    # (a [fm], r0 [fm], L [um], h21 memo R (Campbell), label)
    (52.6, 2.0, 1.0, 1.1e14, "a=52.6, L=1.0um (audit ref)"),
    (52.6, 2.0, 4.2, 2.0e15, "a=52.6, L=4.2um (corpus central L)"),
    (6.8, 2.0, 4.2, 9.4e16, "a=6.8,  L=4.2um (corpus central)"),
    (2.0, 1.0, 4.2, 8.9e17, "a=2.0,  L=4.2um (thin extreme)"),
]

# Corpus-printed comparators (sources cited in RESULTS.md):
G_WINDOW = (8.0e-10, 1.3e-8)      # Omega VII.I body dictionary
G_J3_EDGE = 7.8e-10               # J3 repaired low edge (I2 obligation number)
G_PNU4 = 1.0e-9                   # P-nu4: "Majoron-mode 0nubb at g ~ 1e-9 (far-horizon)"
M3_EV = 0.0468                    # m3 [eV] (J3's central m3)
F_DICT_GEV_PER_SQRTEPS = 1.173    # f = 0.69*sqrt(eps)*1.7 GeV = 1.173*sqrt(eps) GeV (J3 A1)
EPS_BRACKET = (1.5e-6, 2.604e-3)  # B-nu1' floor; repaired ceiling (J2/J3)
N_NUCLEI = 1.0e27                 # spec: macroscopic sample class (~ tonne-scale)


def main():
    rng = np.random.default_rng(RNG_SEED)
    results = {
        "seed": RNG_SEED,
        "parametrization": "REUSED from theory-audit/h21_mc.py: Poisson line process "
                           "rho_L=1/L^2 (phi=pi a^2/L^2), m_hat=1 on tubes, m_bar=phi; "
                           "kernel w(r)=exp(-r/r0)/(4 pi r0 r^2)",
        "staticity": "web static on experiment timescales: h21 sec 6 motional-narrowing "
                     "pricing — narrowing needs web transport 1.5e8-1.2e9 c "
                     "(ledger prints 5e7-1.2e9 c); slow web => frozen snapshots, "
                     "spatial second moment is the observable",
        "comparators": {"g_window": G_WINDOW, "g_J3_edge": G_J3_EDGE,
                        "g_Pnu4": G_PNU4, "m3_eV": M3_EV,
                        "eps_bracket": EPS_BRACKET,
                        "f_dictionary": "f = 1.173*sqrt(eps) GeV"},
    }

    print("== kernel tables ==", flush=True)
    h_tabs = {r0: make_h_table(r0) for r0 in (1.0, 2.0)}
    g_tabs = {}
    for (a, r0, L_um, _, _) in ANCHORS:
        if (a, r0) not in g_tabs:
            sg, hg = h_tabs[r0]
            g_tabs[(a, r0)] = build_g_table(a, r0, sg, hg, extent_r0=100.0)

    # ---------------- G1: h21 anchor reproduction -------------------------------------
    print("== G1: h21 anchor reproduction ==", flush=True)
    g1_rows = []
    for (a, r0, L_um, R_memo, label) in ANCHORS:
        dg, gg, dmax = g_tabs[(a, r0)]
        L = L_um * FM_PER_UM
        rho_L = 1.0 / L ** 2
        phi = np.pi * a ** 2 / L ** 2
        mean_c, var_c = campbell_moments(dg, gg, rho_L)
        R_c = 1.0 + var_c / mean_c ** 2
        row = {"a_fm": a, "r0_fm": r0, "L_um": L_um, "phi": float(phi),
               "R_campbell": float(R_c), "R_h21_memo": R_memo,
               "ratio_to_memo": float(R_c / R_memo),
               "mean_check": float(mean_c / phi)}
        g1_rows.append(row)
        print(f"  {label}: R={R_c:.4g} (memo {R_memo:.3g}, ratio {R_c/R_memo:.3f}); "
              f"<A>/m_bar={mean_c/phi:.6f}", flush=True)
    # direct importance-sampled MC re-anchors at 2 cells (h21 leg 3b)
    mc_anchors = {}
    for (a, r0, L_um) in [(52.6, 2.0, 1.0), (2.0, 1.0, 4.2)]:
        dg, gg, dmax = g_tabs[(a, r0)]
        mc = physical_mc_R(a, r0, L_um * FM_PER_UM, dg, gg, dmax, rng)
        mc_anchors[f"a{a}_r{r0}_L{L_um}"] = mc
        print(f"  direct-MC a={a} r0={r0} L={L_um}um: R={mc['R']:.4g}+-{mc['R_se']:.2g}; "
              f"<A>/m_bar={mc['mean_over_phi']:.4f}", flush=True)
    results["G1"] = {"anchors": g1_rows, "direct_mc": mc_anchors}

    # ---------------- self-averaging (physics point 1) --------------------------------
    print("== self-averaging of the sample total ==", flush=True)
    sa_rows = []
    for (a, r0, L_um, _, label) in ANCHORS:
        dg, gg, dmax = g_tabs[(a, r0)]
        L = L_um * FM_PER_UM
        rho_L = 1.0 / L ** 2
        sa = self_averaging(dg, gg, rho_L, N_NUCLEI, [0.01, 1.0, 1000.0], L)
        sa["anchor"] = label
        sa["a_fm"], sa["r0_fm"], sa["L_um"] = a, r0, L_um
        sa_rows.append(sa)
        worst = max(r["total_RSD"] for r in sa["samples"])
        print(f"  {label}: K4={sa['K4']:.3g}; worst total RSD (D=1cm..1km, N=1e27) = "
              f"{worst:.3g}", flush=True)
    results["self_averaging"] = sa_rows

    # ---------------- G2: burstiness / Lorenz -----------------------------------------
    print("== G2: Lorenz structure (quadrature + MC with jackknife) ==", flush=True)
    k_list = list(range(2, 17))
    lorenz_rows = []
    for (a, r0, L_um, _, label) in ANCHORS:
        dg, gg, dmax = g_tabs[(a, r0)]
        L = L_um * FM_PER_UM
        rho_L = 1.0 / L ** 2
        phi = np.pi * a ** 2 / L ** 2
        C, q, J2, worst_mono = lorenz_quadrature(dg, gg, rho_L, a)
        quad = {}
        for k in k_list:
            quad[k] = share_at_site_fraction(C, q, 10.0 ** (-k))
        quad_ontube = float(np.interp(a, dg, C))
        q_reach = float(q[-1])
        quantq = {f"q{int(l*100)}": site_fraction_at_share(C, q, l)
                  for l in (0.5, 0.9, 0.99)}
        mc_shares, mc_quant = lorenz_mc(a, r0, L, dg, gg, dmax, rng,
                                        [10.0 ** (-k) for k in k_list])
        row = {"a_fm": a, "r0_fm": r0, "L_um": L_um, "label": label,
               "phi": float(phi), "kernel_reach_site_fraction": q_reach,
               "g_monotone_worst_relative_increase": worst_mono,
               "quad_share_top_10^-k": {str(k): quad[k] for k in k_list},
               "quad_on_tube_share": quad_ontube,
               "quad_site_fraction_quantiles": quantq,
               "mc_share_top_10^-k": {str(k): mc_shares[10.0 ** (-k)]
                                      for k in k_list},
               "mc_on_tube_share": mc_shares["on_tube"],
               "mc_site_fraction_quantiles": mc_quant,
               "per_hot_site_rate_over_claimed": float(1.0 / phi ** 2),
               "hot_nuclei_in_N": float(phi * N_NUCLEI)}
        lorenz_rows.append(row)
        print(f"  {label}: phi={phi:.3g}; reach fraction={q_reach:.3g}; "
              f"on-tube share quad={quad_ontube:.4f} "
              f"mc={mc_shares['on_tube']['share']:.4f}+-{mc_shares['on_tube']['se']:.1g}; "
              f"q50={quantq['q50']:.3g} (mc {mc_quant['q50']['site_fraction']:.3g}"
              f"+-{mc_quant['q50']['site_fraction_se']:.1g})", flush=True)
        print(f"    top-1e-14 share: quad={quad[14]:.4f} "
              f"mc={mc_shares[1e-14]['share']:.4f}+-{mc_shares[1e-14]['se']:.1g}; "
              f"top-1e-12 share = {quad[12]:.6f}", flush=True)
    results["G2_lorenz"] = lorenz_rows

    # ---------------- G3: propagation table + verdict ---------------------------------
    print("== G3: propagation table ==", flush=True)
    prop_rows = []
    for row1, (a, r0, L_um, _, label) in zip(g1_rows, ANCHORS):
        R = row1["R_campbell"]
        sqrtR = float(np.sqrt(R))
        entry = {"anchor": label, "a_fm": a, "r0_fm": r0, "L_um": L_um,
                 "R": R, "sqrtR": sqrtR, "g_equiv": {}}
        for name, g in [("window_low_8e-10", G_WINDOW[0]),
                        ("J3_edge_7.8e-10", G_J3_EDGE),
                        ("Pnu4_1e-9", G_PNU4),
                        ("window_top_1.3e-8", G_WINDOW[1])]:
            ge = sqrtR * g
            entry["g_equiv"][name] = {
                "g": g, "g_equiv": float(ge),
                "ratio_to_Pnu4_line": float(ge / G_PNU4),
                "ratio_to_window_top": float(ge / G_WINDOW[1])}
        # retreat arithmetic: hold the observable at the printed P-nu4 class
        # => g_req = g/sqrt(R), f_req = sqrt(R) * f, eps_req = (f_req/1.173 GeV)^2
        f_pnu4_MeV = M3_EV / G_PNU4 * 1e-6      # m3/g in eV -> MeV
        f_req_MeV = sqrtR * f_pnu4_MeV
        eps_req = (f_req_MeV * 1e-3 / F_DICT_GEV_PER_SQRTEPS) ** 2
        entry["retreat"] = {
            "g_req_for_printed_observable": float(G_PNU4 / sqrtR),
            "f_req_MeV": float(f_req_MeV),
            "f_req_TeV": float(f_req_MeV * 1e-6),
            "eps_req": float(eps_req),
            "eps_ceiling_printed": EPS_BRACKET[1],
            "eps_overshoot_factor": float(eps_req / EPS_BRACKET[1]),
            "f_window_printed_MeV": [4.0, 60.0]}
        prop_rows.append(entry)
        print(f"  {label}: R={R:.3g}, sqrtR={sqrtR:.3g}; "
              f"g_equiv(J3 7.8e-10)={sqrtR*G_J3_EDGE:.3g}; "
              f"g_equiv(top 1.3e-8)={sqrtR*G_WINDOW[1]:.3g}; "
              f"retreat: f_req={f_req_MeV*1e-6:.3g} TeV, eps_req={eps_req:.3g} "
              f"(ceiling overshoot {eps_req/EPS_BRACKET[1]:.2g})", flush=True)
    results["G3_propagation"] = prop_rows

    g_equiv_min = min(e["g_equiv"]["J3_edge_7.8e-10"]["g_equiv"] for e in prop_rows)
    g_equiv_max = max(e["g_equiv"]["window_top_1.3e-8"]["g_equiv"] for e in prop_rows)
    results["G3_verdict"] = {
        "frame": "(i) survives-with-repricing / (ii) overshoots-a-printed-bound "
                 "/ (iii) indeterminate-from-archive — pre-registered, no tuning",
        "battery_rows": "R1-R4 + perimeter constrain the ACTUAL coupling g via "
                        "astro/lab emission; R multiplies only the 0nubb insertion's "
                        "second moment => the J3 battery verdict itself is untouched",
        "route_i_closed": "holding the printed far-horizon observable requires "
                          "f_req = sqrt(R) x f = 4.9e2-4.4e4 TeV, i.e. eps_req "
                          "1.8e11-1.4e15 — overshoots the corpus's printed epsilon "
                          "ceiling 2.604e-3 by 6.8e13-5.5e17 and leaves the printed "
                          "f-window [4,60] MeV by 7-9 orders: no admissible retreat "
                          "inside the corpus's own printed parameter space",
        "g_equiv_range_over_table": [float(g_equiv_min), float(g_equiv_max)],
        "verdict": "(iii) INDETERMINATE-FROM-ARCHIVE",
        "blocking_unprinted_number": "g_lim — the experimental Majoron-mode 0nubb "
                                     "limit / detectability line in g units "
                                     "(equivalently the mode's rate normalization "
                                     "T_1/2(g)); nowhere printed (J3 sec 0 census)",
        "conditional_kill": "for ANY archive-external g_lim < 8.2e-3, the honest rate "
                            "is excluded at EVERY anchor and EVERY g in the printed "
                            "window [8e-10, 1.3e-8] — P-nu4 inverts from far-horizon "
                            "to already-excluded; filed as F-T7-N2",
        "honesty_note": "F-R9 wounds the K3a 0nubb MODE-rate claim only; S1 (Sigma "
                        "m_nu) is untouched either way; S6's occurrence stake, hash, "
                        "funnels, kill conditions verbatim-unchanged (Watch memo)"}

    results["runtime_s"] = time.time() - T_START
    out = ("/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/N2/"
           "n2_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {out} ({results['runtime_s']:.0f}s)")


if __name__ == "__main__":
    main()
