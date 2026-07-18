#!/usr/bin/env python3
"""
Tier-7 M3 — F-R9 route (b) priced: on-tube saturation nonlinearity (0nubb insert).

Question. h21 (theory-audit) established that with line-supported mixing density the
0nubb rate observable (quadratic in the amplitude A) carries a second-moment excess
R = <A^2>/<A>^2 = 1e14-1e18 across the corpus's own scale window, while the corpus's
FIRST-moment claim <A>/m_bar = 1.000 +/- 0.001 is verified. The motional-narrowing
rescue is causally excluded (needs web speeds 1e8-1e9 c). The one undischarged rescue
(h21 caveat #2, F-R9 route (b)): a saturation nonlinearity of the insertion on tubes.
This script prices it.

Saturation model (per ROADMAP_v9 M3): A -> A_s = A_sat * tanh(A / A_sat), applied to
the point amplitude. Scan A_sat/<A> over >= 4 decades (we do 20); measure
    R(A_sat)       = <A_s^2>/<A_s>^2                       (fluctuation excess)
    D_true(A_sat)  = <A_s^2>/<A^2>                         (roadmap G2: <rate>_sat/<rate>)
    D_claim(A_sat) = <A_s^2>/m_bar^2                       (vs the corpus's CLAIMED mean
                                                            rate; h21's mean claim is
                                                            <A> = m_bar, rate ~ m_bar^2)
    M1(A_sat)      = <A_s>/m_bar                           (h21's verified 1.000-check)

Parametrization: REUSED from theory-audit/h21_mc.py verbatim (gate condition).
  * Static isotropic Poisson line process, tube radius a, uniform pitch amplitude
    m_hat = 1 inside, 0 outside; rho_L = 1/L_web^2 (phi = pi a^2/L^2); m_bar = phi.
  * Kernel w(r) = exp(-r/r0)/(4 pi r0 r^2), r0 in {1, 2} fm (h21's acknowledged choice).
  * Scale window: L in {1, 1.7, 4.2, 10} um; a in {2, 6.8, 20, 52.6, 200} fm.
  * One-tube response g(d) tables: identical code (make_h_table / g_of_d / build_g_table).

Saturated moments. At h21's physical scales the encounter probability is
lam ~ pi rho_L d_max^2 ~ 1e-13, so A at a point is (up to O(lam) relative corrections)
the single nearest-tube contribution m_hat*g(d), d distributed with intensity
2 pi rho_L d dd. Generalized (nonlinear) Campbell, sparse limit:
    <S(A)>   = rho_L * int 2 pi d S(g(d)) dd
    <S(A)^2> = rho_L * int 2 pi d S(g(d))^2 dd  +  <S(A)>^2
Exact for linear S (Campbell's theorem, h21 leg 2); for tanh the multi-tube error is
O(lam) relative — VALIDATED here against explicit-network MC at moderate separation
(phi ~ 1e-3..1e-2, where the multi-tube error is at its LARGEST) and against
importance-sampled direct MC at physical scale.

Deterministic seed. Output: m3_results.json, m3_fig.png.
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


def sample_poisson_lines(rho_L, R_b, rng):
    n = rng.poisson(rho_L * np.pi * R_b ** 2)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    ref = np.where(np.abs(u[:, [2]]) < 0.9, [[0, 0, 1.0]], [[1.0, 0, 0]])
    e1 = np.cross(u, ref)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(u, e1)
    rad = R_b * np.sqrt(rng.random(n))
    th = 2 * np.pi * rng.random(n)
    p = rad[:, None] * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
    return p, u


def sample_lattice_lines(rho_L, R_b, rng, jitter_frac=0.25):
    L_lat = np.sqrt(3.0 / rho_L)
    m = int(np.ceil(2 * R_b / L_lat)) + 2
    base = (np.arange(m) - m / 2) * L_lat
    ps, us = [], []
    for axis in range(3):
        off = rng.uniform(0, L_lat, 2)
        t1, t2 = np.meshgrid(base + off[0], base + off[1], indexing="ij")
        t1 = t1.ravel() + rng.uniform(-jitter_frac, jitter_frac, t1.size) * L_lat
        t2 = t2.ravel() + rng.uniform(-jitter_frac, jitter_frac, t2.size) * L_lat
        keep = (t1 ** 2 + t2 ** 2) <= (R_b + L_lat) ** 2
        t1, t2 = t1[keep], t2[keep]
        p = np.zeros((t1.size, 3))
        u = np.zeros((t1.size, 3))
        u[:, axis] = 1.0
        tr = [i for i in range(3) if i != axis]
        p[:, tr[0]] = t1
        p[:, tr[1]] = t2
        ps.append(p)
        us.append(u)
    return np.vstack(ps), np.vstack(us)


# ============================================================================
# Section 2 — NEW: saturation
# ============================================================================

def S_of(x, A_sat):
    """Saturation nonlinearity A -> A_sat tanh(A/A_sat); identity for A_sat = inf."""
    x = np.asarray(x, dtype=float)
    if np.isinf(A_sat):
        return x.copy()
    return A_sat * np.tanh(x / A_sat)


def sat_moments(d, g, rho_L, A_sat, m_hat=1.0):
    """Sparse-limit generalized Campbell for the saturated amplitude.
    Returns (<A_s>, <A_s^2>, R) with <A_s^2> = kappa2 + kappa1^2."""
    Sg = S_of(m_hat * g, A_sat)
    k1 = rho_L * np.trapezoid(2 * np.pi * d * Sg, d)
    k2 = rho_L * np.trapezoid(2 * np.pi * d * Sg ** 2, d)
    m1 = k1
    m2 = k2 + k1 ** 2
    return m1, m2, m2 / m1 ** 2


def network_mc_sat(kind, a, r0, phi, d_grid, g_grid, d_max, A_sat_list, rng,
                   n_blocks=4, n_per_block=500_000, chunk=100_000):
    """Explicit-network MC (h21 leg 3a, extended): moments of S(A) for each A_sat,
    saturation applied to the FULL multi-tube sum A — this is the ground truth the
    sparse-limit formula is checked against."""
    rho_L = phi / (np.pi * a ** 2)
    L_eff = 1.0 / np.sqrt(rho_L)
    R_b = 5.0 * L_eff + 2.0 * d_max
    R_in = R_b - d_max
    nA = len(A_sat_list)
    blocks = np.zeros((n_blocks, nA, 2))          # m1, m2 per block per A_sat
    for b in range(n_blocks):
        if kind == "poisson":
            p, u = sample_poisson_lines(rho_L, R_b, rng)
        else:
            p, u = sample_lattice_lines(rho_L, R_b, rng)
        s1 = np.zeros(nA)
        s2 = np.zeros(nA)
        done = 0
        while done < n_per_block:
            n = min(chunk, n_per_block - done)
            x = rng.normal(size=(n, 3))
            x /= np.linalg.norm(x, axis=1)[:, None]
            x *= R_in * rng.random(n)[:, None] ** (1.0 / 3.0)
            A = np.zeros(n)
            for i in range(p.shape[0]):
                dvec = x - p[i]
                proj = dvec @ u[i]
                per2 = np.einsum("ij,ij->i", dvec, dvec) - proj ** 2
                near = per2 < d_max ** 2
                if near.any():
                    A[near] += np.interp(np.sqrt(per2[near]), d_grid, g_grid)
            for k, A_sat in enumerate(A_sat_list):
                SA = S_of(A, A_sat)
                s1[k] += SA.sum()
                s2[k] += (SA ** 2).sum()
            done += n
        blocks[b, :, 0] = s1 / n_per_block
        blocks[b, :, 1] = s2 / n_per_block
    out = []
    for k, A_sat in enumerate(A_sat_list):
        m1 = blocks[:, k, 0]
        m2 = blocks[:, k, 1]
        Rb = m2 / m1 ** 2
        out.append({"A_sat": float(A_sat), "m1": float(m1.mean()),
                    "m1_se": float(m1.std(ddof=1) / np.sqrt(n_blocks)),
                    "m2": float(m2.mean()),
                    "m2_se": float(m2.std(ddof=1) / np.sqrt(n_blocks)),
                    "R": float(Rb.mean()),
                    "R_se": float(Rb.std(ddof=1) / np.sqrt(n_blocks))})
    return out


def physical_mc_sat(a, r0, L_web, d_grid, g_grid, d_max, A_sat, rng, n=2_000_000):
    """h21's importance-sampled physical-scale MC (leg 3b), with saturation applied
    to the single-tube amplitude (sparse limit; lam ~ 1e-13 corrections)."""
    rho_L = 1.0 / L_web ** 2
    d_lo = 1e-3 * min(a, r0)
    lnZ = np.log(d_max / d_lo)
    d = d_lo * np.exp(rng.random(n) * lnZ)
    gg = np.interp(d, d_grid, g_grid)
    Sg = S_of(gg, A_sat)
    w1 = 2 * np.pi * d * Sg * (d * lnZ)
    w2 = 2 * np.pi * d * Sg ** 2 * (d * lnZ)
    I1 = w1.mean()
    I2 = w2.mean()
    m1 = rho_L * I1
    m2 = rho_L * I2 + m1 ** 2
    R = m2 / m1 ** 2
    # jackknife over 20 slices
    k = 20
    slices = np.array_split(np.arange(n), k)
    Rj = np.empty(k)
    S1, S2 = w1.sum(), w2.sum()
    for i, s in enumerate(slices):
        n_i = n - s.size
        I1_i = (S1 - w1[s].sum()) / n_i
        I2_i = (S2 - w2[s].sum()) / n_i
        m1_i = rho_L * I1_i
        Rj[i] = (rho_L * I2_i + m1_i ** 2) / m1_i ** 2
    R_se = np.sqrt((k - 1) / k * np.sum((Rj - Rj.mean()) ** 2))
    return {"R": float(R), "R_se": float(R_se), "m1": float(m1), "m2": float(m2)}


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


def main():
    rng = np.random.default_rng(RNG_SEED)
    results = {
        "seed": RNG_SEED,
        "parametrization": "REUSED from theory-audit/h21_mc.py: Poisson line process "
                           "rho_L=1/L^2 (phi=pi a^2/L^2), m_hat=1 on tubes, m_bar=phi; "
                           "kernel w(r)=exp(-r/r0)/(4 pi r0 r^2); scale window "
                           "L in {1,1.7,4.2,10} um, a in {2,6.8,20,52.6,200} fm",
        "saturation_model": "A -> A_sat * tanh(A/A_sat) applied to the point amplitude",
        "D_definitions": {
            "D_true": "<A_s^2>/<A^2>  (roadmap G2: <rate>_sat/<rate>, rate ~ A^2, "
                      "denominator = the unsaturated model's true mean rate)",
            "D_claim": "<A_s^2>/m_bar^2  (vs the corpus's CLAIMED mean rate: NR-K3a says "
                       "the amplitude is governed by the volume mean m_bar, verified by "
                       "h21 as <A>/m_bar = 1.000+-0.001; the claimed rate is m_bar^2)",
            "M1": "<A_s>/m_bar  (h21's first-moment 1.000-check under saturation)"},
    }

    print("== kernel tables ==", flush=True)
    h_tabs = {r0: make_h_table(r0) for r0 in (1.0, 2.0)}
    g_tabs = {}
    for (a, r0, L_um, _, _) in ANCHORS:
        if (a, r0) not in g_tabs:
            sg, hg = h_tabs[r0]
            g_tabs[(a, r0)] = build_g_table(a, r0, sg, hg, extent_r0=100.0)

    # ---------------- G1: unsaturated reproduction of h21 -----------------------------
    print("== G1: h21 reproduction (unsaturated control) ==", flush=True)
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
               "mean_check": float(mean_c / phi),
               "g0": float(gg[0]),
               "A_on_tube_over_mean": float(gg[0] / phi)}
        g1_rows.append(row)
        print(f"  {label}: R={R_c:.4g} (memo {R_memo:.3g}, ratio {R_c/R_memo:.3f}); "
              f"<A>/m_bar={mean_c/phi:.6f}; g(0)={gg[0]:.4f}", flush=True)
    # physical-scale direct MC re-anchor (h21 leg 3b) at the audit reference cell
    a, r0, L_um = 52.6, 2.0, 1.0
    dg, gg, dmax = g_tabs[(a, r0)]
    mc = physical_mc_sat(a, r0, L_um * FM_PER_UM, dg, gg, dmax, np.inf, rng)
    print(f"  direct-MC re-anchor (a=52.6, r0=2, L=1um, unsat): "
          f"R={mc['R']:.4g}+-{mc['R_se']:.2g} (h21 leg-3b same cell: "
          f"1.115e14+-0.003e14)", flush=True)
    results["G1"] = {"anchors": g1_rows, "direct_mc_unsat_a52.6_L1": mc,
                     "h21_direct_mc_value_same_cell": 1.115e14}

    # ---------------- validation: sparse-limit sat formula vs network MC --------------
    print("== validation: saturated sparse-limit formula vs explicit-network MC ==",
          flush=True)
    val = []
    for (a_v, phi_v, kinds) in [(0.5, 1.5e-3, ("poisson",)),
                                (2.0, 3.0e-3, ("poisson", "lattice"))]:
        r0_v = 1.0
        sg, hg = h_tabs[r0_v]
        dgv, ggv, dmaxv = build_g_table(a_v, r0_v, sg, hg, extent_r0=60.0)
        rho_Lv = phi_v / (np.pi * a_v ** 2)
        g0 = ggv[0]
        A_sat_list = [0.05 * g0, 0.2 * g0, 1.0 * g0, np.inf]
        preds = []
        for A_sat in A_sat_list:
            m1, m2, R = sat_moments(dgv, ggv, rho_Lv, A_sat)
            preds.append({"A_sat": float(A_sat), "m1": m1, "m2": m2, "R": R})
        cell = {"a": a_v, "r0": r0_v, "phi": phi_v, "g0": float(g0),
                "sparse_limit_pred": preds, "mc": {}}
        for kind in kinds:
            t0 = time.time()
            mcrows = network_mc_sat(kind, a_v, r0_v, phi_v, dgv, ggv, dmaxv,
                                    A_sat_list, rng)
            cell["mc"][kind] = mcrows
            for pr, mr in zip(preds, mcrows):
                print(f"  a={a_v} phi={phi_v} {kind} A_sat={mr['A_sat']:.3g}: "
                      f"R_mc={mr['R']:.1f}+-{mr['R_se']:.1f} vs pred {pr['R']:.1f} "
                      f"(ratio {mr['R']/pr['R']:.3f})", flush=True)
            print(f"    [{time.time()-t0:.0f}s]", flush=True)
        val.append(cell)
    results["validation_network_mc"] = val

    # physical-scale saturated MC vs quadrature (audit ref cell)
    print("== validation: saturated physical-scale MC vs quadrature ==", flush=True)
    a, r0, L_um = 52.6, 2.0, 1.0
    dg, gg, dmax = g_tabs[(a, r0)]
    L = L_um * FM_PER_UM
    rho_L = 1.0 / L ** 2
    phi = np.pi * a ** 2 / L ** 2
    val_phys = []
    for x in (1e3, 1e7, 1e11, 1e15):
        A_sat = x * phi
        m1q, m2q, Rq = sat_moments(dg, gg, rho_L, A_sat)
        mc = physical_mc_sat(a, r0, L, dg, gg, dmax, A_sat, rng)
        val_phys.append({"A_sat_over_mean": x, "R_quad": Rq, "R_mc": mc["R"],
                         "R_mc_se": mc["R_se"], "ratio": mc["R"] / Rq})
        print(f"  x=A_sat/<A>={x:.0e}: R_quad={Rq:.4g} R_mc={mc['R']:.4g}"
              f"+-{mc['R_se']:.2g} (ratio {mc['R']/Rq:.4f})", flush=True)
    results["validation_physical_mc"] = val_phys

    # ---------------- G2: the pricing scan --------------------------------------------
    print("== G2: pricing scan R(A_sat), D(A_sat) ==", flush=True)
    x_grid = np.geomspace(1e-2, 1e18, 81)     # A_sat/<A>, 20 decades
    scans = []
    for (a, r0, L_um, R_memo, label) in ANCHORS:
        dg, gg, dmax = g_tabs[(a, r0)]
        L = L_um * FM_PER_UM
        rho_L = 1.0 / L ** 2
        phi = np.pi * a ** 2 / L ** 2
        mean_u, var_u = campbell_moments(dg, gg, rho_L)
        m2_u = var_u + mean_u ** 2            # unsaturated <A^2>
        R_u = m2_u / mean_u ** 2
        R_arr, Dt_arr, Dc_arr, M1_arr = [], [], [], []
        for x in x_grid:
            A_sat = x * phi                    # A_sat in units of m_hat; <A> = m_bar = phi
            m1, m2, R = sat_moments(dg, gg, rho_L, A_sat)
            R_arr.append(R)
            Dt_arr.append(m2 / m2_u)
            Dc_arr.append(m2 / phi ** 2)
            M1_arr.append(m1 / phi)
        R_arr = np.array(R_arr)
        Dt_arr = np.array(Dt_arr)
        Dc_arr = np.array(Dc_arr)
        M1_arr = np.array(M1_arr)

        def crossing_down(y, target):
            """largest x where y crosses down through target (y decreasing with x?
            here R and D increase with x; find x where y first reaches target from
            below as x decreases — i.e. crossing on the decreasing-x direction)."""
            for i in range(len(x_grid) - 1, 0, -1):
                y1, y0 = y[i], y[i - 1]
                if (y0 - target) * (y1 - target) <= 0 and y0 != y1:
                    t = (np.log(target) - np.log(y0)) / (np.log(y1) - np.log(y0))
                    return float(np.exp(np.log(x_grid[i - 1]) +
                                        t * (np.log(x_grid[i]) - np.log(x_grid[i - 1]))))
            return None

        # A_sat* where R <= 10 (scan)
        idx_ok = np.where(R_arr <= 10.0)[0]
        A_sat_star = float(x_grid[idx_ok[-1]]) if idx_ok.size else None
        # D_claim = 1 crossing and R there
        x_Dc1 = crossing_down(Dc_arr, 1.0)
        R_at_Dc1 = (float(np.exp(np.interp(np.log(x_Dc1), np.log(x_grid),
                                           np.log(R_arr)))) if x_Dc1 else None)
        # D_true = 0.5 crossing (the |D-1|<=0.5 boundary, roadmap definition) and R there
        x_Dt05 = crossing_down(Dt_arr, 0.5)
        R_at_Dt05 = (float(np.exp(np.interp(np.log(x_Dt05), np.log(x_grid),
                                            np.log(R_arr)))) if x_Dt05 else None)
        # trade slopes in the deep-clip regime (x in [1e2, 1e10])
        sel = (x_grid >= 1e2) & (x_grid <= 1e10)
        lx = np.log(x_grid[sel])
        slope_D = float(np.polyfit(lx, np.log(Dc_arr[sel]), 1)[0])
        slope_R = float(np.polyfit(lx, np.log(R_arr[sel]), 1)[0])

        # analytic requirement for R = 10 (exponential-tail extrapolation of g)
        d_req = 1.0 / np.sqrt((10.0 - 1.0) * np.pi * rho_L)
        i_tail = np.where(gg > 0)[0][-1]
        log10_A_req_mhat = (np.log10(gg[i_tail])
                            - (d_req - dg[i_tail]) / (r0 * np.log(10.0)))
        log10_A_req_over_mean = log10_A_req_mhat - np.log10(phi)
        # D_claim at that point: <A_s^2> ~ rho_L pi d_req^2 A_req^2
        log10_D_req = (np.log10(np.pi * rho_L * d_req ** 2)
                       + 2 * log10_A_req_mhat - 2 * np.log10(phi))

        scan = {"a_fm": a, "r0_fm": r0, "L_um": L_um, "phi": float(phi),
                "label": label,
                "R_unsat": float(R_u), "g0": float(gg[0]),
                "x_grid_A_sat_over_mean": [float(v) for v in x_grid],
                "R": [float(v) for v in R_arr],
                "D_true": [float(v) for v in Dt_arr],
                "D_claim": [float(v) for v in Dc_arr],
                "M1": [float(v) for v in M1_arr],
                "A_sat_star_R_le_10": A_sat_star,
                "min_R_in_scan": float(R_arr.min()),
                "min_R_at_x": float(x_grid[int(np.argmin(R_arr))]),
                "x_at_D_claim_1": x_Dc1, "R_at_D_claim_1": R_at_Dc1,
                "x_at_D_true_05": x_Dt05, "R_at_D_true_05": R_at_Dt05,
                "dlnD_claim_dlnAsat_clip": slope_D, "dlnR_dlnAsat_clip": slope_R,
                "analytic_R10_requirement": {
                    "d_req_fm": float(d_req),
                    "log10_A_sat_req_over_mhat": float(log10_A_req_mhat),
                    "log10_A_sat_req_over_mean": float(log10_A_req_over_mean),
                    "log10_D_claim_at_req": float(log10_D_req)}}
        scans.append(scan)
        print(f"  {label}: R_unsat={R_u:.3g}; min R over 20-decade scan = "
              f"{R_arr.min():.3g} (at x={x_grid[np.argmin(R_arr)]:.1g}); "
              f"A_sat*(R<=10): {A_sat_star}; "
              f"D_claim=1 at x={x_Dc1:.3g} where R={R_at_Dc1:.3g}; "
              f"slopes clip regime: dlnD/dlnAsat={slope_D:.3f}, "
              f"dlnR/dlnAsat={slope_R:.4f}", flush=True)
        print(f"    analytic R=10 requirement: d_req={d_req:.3g} fm, "
              f"log10(A_sat/m_hat)={log10_A_req_mhat:.3g}, "
              f"log10 D_claim there = {log10_D_req:.3g}", flush=True)
    results["G2_scans"] = scans

    # ---------------- G3: verdict ------------------------------------------------------
    any_star = any(s["A_sat_star_R_le_10"] is not None for s in scans)
    all_D_fail = True  # evaluated only if a star existed
    verdict = ("NOT-VIABLE" if not any_star else "see per-anchor")
    results["G3"] = {
        "decision_rule": "VIABLE iff |D-1| <= 0.5 at A_sat* (R<=10) AND A_sat* inside "
                         "h21's on-tube amplitude range; CONDITIONAL if range not "
                         "extractable; else NOT-VIABLE",
        "A_sat_star_exists_anywhere_in_scan": bool(any_star),
        "on_tube_amplitude_range_extractable": True,
        "on_tube_amplitude_over_mean_range": [
            min(r["A_on_tube_over_mean"] for r in g1_rows),
            max(r["A_on_tube_over_mean"] for r in g1_rows)],
        "verdict": verdict}
    print(f"== G3 verdict: {verdict} ==", flush=True)

    results["runtime_s"] = time.time() - T_START
    out = ("/home/user/fork_frankensim/analysis/gum-sandbox/tier7-program/M3/"
           "m3_results.json")
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {out} ({results['runtime_s']:.0f}s)")


if __name__ == "__main__":
    main()
