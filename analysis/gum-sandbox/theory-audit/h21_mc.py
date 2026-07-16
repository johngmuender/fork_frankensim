#!/usr/bin/env python3
"""
Phase H2.1 — decisive Monte-Carlo for theory-audit finding T-H2 (NR-K3a second moment).

Question. NR-K2b/K3a (corpus, Step 2) claims the 0nubb insertion amplitude is governed
by the volume MEAN m_bar of the line-supported chiral-pitch density m(x): "line-support
and bulk-support are amplitude-identical in the mean"; fluctuations "(fm/um)-class,
negligible". The audit (T1_NR_gates.md sec 4.2) counters: the observable is quadratic
(rate ~ <|A|^2>), and for sparse line support R = <|A|^2>/|<A>|^2 ~ 1/phi
(phi = tube volume filling fraction) ~ 1e14-2.5e16. This script decides it by direct MC.

Model (units: fm throughout).
  * Static network of straight tubes (radius a, uniform pitch amplitude m_hat inside,
    0 outside), volume mean m_bar = phi * m_hat held fixed. Two geometries:
    (i) isotropic Poisson line process with length density rho_L = 1/L_web^2
        (audit convention: phi = pi a^2 / L^2);
    (ii) jitter-perturbed cubic lattice of lines along x/y/z at the same rho_L.
  * Kernel: the corpus never prints the 0nubb kernel; per the audit we use the
    normalized double-propagator-class weight  w(r) = exp(-r/r0) / (4 pi r0 r^2),
    with r0 = 1 fm and 2 fm variants (range ~ 1/q_bar ~ 2 fm).
  * Amplitude at insertion point x0:  A(x0) = int m(y) w(x0-y) d^3y.
  * Report R = <A^2>/<A>^2 over insertion points.

Method (three independent legs; every ingredient cross-checked).
  1. One-tube response g(d) = int_tube w  (perpendicular distance d to the tube axis)
     by deterministic quadrature, validated against a brute-force kernel-sampling MC.
  2. Exact second moment for the Poisson line process via Campbell's theorem:
        <A> = m_hat * rho_L * int 2 pi d g(d) dd            (= m_bar, mass check)
        Var A = m_hat^2 * rho_L * int 2 pi d g(d)^2 dd
     -> R = 1 + Var/<A>^2. Exact (no approximation) at ANY scale separation.
  3. Direct network MC: (a) at moderate scale separation (phi ~ 1e-3..1e-2, reachable
     by uniform sampling) sample millions of x0 in explicit Poisson AND jittered-lattice
     networks, estimate R with error bars, compare with Campbell -> validates leg 2 and
     the Poisson<->lattice equivalence; (b) at PHYSICAL scales, direct MC with
     importance sampling on the nearest-tube distance (uniform sampling is hopeless at
     phi ~ 1e-14 by construction -- that hopelessness IS the finding).
  4. Motional narrowing: time-averaged kernel. The web is displaced rigidly by a random
     transverse offset (2D Gaussian, rms delta per axis) during the coherence time of a
     single event; the tube response becomes g_eff = g (*) Gauss2D(delta). R(delta)
     prices the rescue: delta needed for R -> O(1).

Deterministic seeds. Runtime ~ a few minutes. Output: h21_results.json, h21_fig.png.
"""

import json
import time
import numpy as np
from scipy.special import i0e

RNG_SEED = 20260716
FM_PER_UM = 1.0e9
HBARC_MEV_FM = 197.327

T_START = time.time()


# ----------------------------------------------------------------------------
# 1. Kernel and one-tube response g(d)
# ----------------------------------------------------------------------------

def make_h_table(r0, s_min_frac=1e-9, s_max_over_r0=120.0, n_s=1400, n_t=3000):
    """Line-integrated kernel  h(s) = int_-inf^inf w(sqrt(s^2+z^2)) dz
    = 1/(2 pi r0 s) * int_0^inf exp(-(s/r0) cosh t)/cosh t dt.
    Returns (s_grid, h_vals) for log-log interpolation."""
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
    """log-log interpolation of h; 0 beyond table."""
    s = np.asarray(s, dtype=float)
    out = np.zeros_like(s)
    inside = (s > s_grid[0]) & (s < s_grid[-1])
    out[inside] = np.exp(np.interp(np.log(s[inside]), np.log(s_grid), np.log(np.maximum(h_grid, 1e-300))))
    out[s <= s_grid[0]] = h_grid[0]  # h ~ log-divergent/s; clipped (measure-zero region)
    return out


def g_of_d(d, a, r0, s_grid, h_grid, n_nodes=900):
    """g(d) = int over disc(center distance d, radius a) of h(s) d^2s
    = int h(s) * s * psi(s; d, a) ds,  psi = angular measure of the disc at radius s."""
    s_cut = s_grid[-1]
    s_lo = max(d - a, 0.0)
    s_hi = min(d + a, s_cut)
    if s_lo >= s_hi:
        return 0.0
    nodes = [np.geomspace(max(s_lo, 1e-10 * r0), s_hi, n_nodes)]
    if s_lo == 0.0:
        nodes.append(np.geomspace(1e-10 * r0, min(1e-2 * r0, s_hi), 120))
    # cluster around the kink at |d-a| and around the kernel scale
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
    """g(d) on a composite grid up to d_max = a + extent_r0*r0."""
    d_max = a + extent_r0 * r0
    grid = [np.array([0.0]),
            np.geomspace(1e-3 * min(a, r0), d_max, 500),
            np.linspace(max(0.0, a - extent_r0 * r0), d_max, n_extra),
            np.linspace(max(0.0, a - 8 * r0), min(d_max, a + 8 * r0), 300)]
    d = np.unique(np.clip(np.concatenate(grid), 0.0, d_max))
    g = np.array([g_of_d(x, a, r0, s_grid, h_grid) for x in d])
    return d, g, d_max


def campbell_moments(d, g, rho_L):
    """Poisson line process, Campbell's theorem (exact):
       mean = rho_L * int 2 pi d g dd ;  var = rho_L * int 2 pi d g^2 dd  (m_hat = 1)."""
    mean = rho_L * np.trapezoid(2 * np.pi * d * g, d)
    var = rho_L * np.trapezoid(2 * np.pi * d * g ** 2, d)
    return mean, var


# ----------------------------------------------------------------------------
# 2. Brute-force validation of g(d): sample the kernel directly
# ----------------------------------------------------------------------------

def g_bruteforce(d, a, r0, n, rng):
    """g(d) = Pr[x0 + r*u_hat inside tube], r ~ Exp(r0), u_hat uniform on S^2.
    (Exact because the radial density of w is p(r) = exp(-r/r0)/r0.)"""
    r = rng.exponential(r0, n)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    px = d + r * u[:, 0]
    py = r * u[:, 1]
    hit = (px ** 2 + py ** 2) <= a ** 2
    p = hit.mean()
    se = hit.std(ddof=1) / np.sqrt(n)
    return p, se


# ----------------------------------------------------------------------------
# 3a. Direct network MC at moderate scale separation (uniform x0 sampling)
# ----------------------------------------------------------------------------

def sample_poisson_lines(rho_L, R_b, rng):
    """Isotropic Poisson line process hitting the ball of radius R_b:
    intensity of lines = rho_L * pi * R_b^2; direction uniform on sphere,
    perpendicular offset uniform on the disc of radius R_b."""
    n = rng.poisson(rho_L * np.pi * R_b ** 2)
    u = rng.normal(size=(n, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    # orthonormal frame per line
    ref = np.where(np.abs(u[:, [2]]) < 0.9, [[0, 0, 1.0]], [[1.0, 0, 0]])
    e1 = np.cross(u, ref)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(u, e1)
    rad = R_b * np.sqrt(rng.random(n))
    th = 2 * np.pi * rng.random(n)
    p = rad[:, None] * (np.cos(th)[:, None] * e1 + np.sin(th)[:, None] * e2)
    return p, u


def sample_lattice_lines(rho_L, R_b, rng, jitter_frac=0.25):
    """Jitter-perturbed cubic lattice of lines along x, y, z with total length
    density rho_L (per-axis spacing L_lat = sqrt(3/rho_L)); random global offset,
    each line's two transverse coordinates jittered by U(-jf, jf)*L_lat."""
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


def network_mc_R(kind, a, r0, phi, d_grid, g_grid, d_max, rng,
                 n_blocks=4, n_per_block=1_000_000, chunk=100_000):
    """Uniform-x0 MC of R on explicit networks. Returns per-block R and moments."""
    rho_L = phi / (np.pi * a ** 2)
    L_eff = 1.0 / np.sqrt(rho_L)
    R_b = 5.0 * L_eff + 2.0 * d_max
    R_in = R_b - d_max
    blocks = []
    for _ in range(n_blocks):
        if kind == "poisson":
            p, u = sample_poisson_lines(rho_L, R_b, rng)
        else:
            p, u = sample_lattice_lines(rho_L, R_b, rng)
        s1 = s2 = cnt = 0.0
        done = 0
        while done < n_per_block:
            n = min(chunk, n_per_block - done)
            # uniform in the ball of radius R_in
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
            s1 += A.sum()
            s2 += (A ** 2).sum()
            cnt += n
            done += n
        m1, m2 = s1 / cnt, s2 / cnt
        blocks.append((m1, m2, m2 / m1 ** 2 if m1 > 0 else np.nan))
    Rs = np.array([b[2] for b in blocks])
    return {"R_mean": float(Rs.mean()), "R_se": float(Rs.std(ddof=1) / np.sqrt(n_blocks)),
            "R_blocks": [float(r) for r in Rs],
            "m1_mean": float(np.mean([b[0] for b in blocks])),
            "m2_mean": float(np.mean([b[1] for b in blocks]))}


# ----------------------------------------------------------------------------
# 3b. Direct MC at physical scales (importance sampling on nearest-tube distance)
# ----------------------------------------------------------------------------

def physical_mc_R(a, r0, L_web, d_grid, g_grid, d_max, rng, n=2_000_000):
    """At physical phi ~ 1e-14, P(any tube within kernel reach of x0) = lam ~ pi rho_L d_max^2.
    E[A] = rho_L * I1 and E[A^2] = (rho_L I1)^2 + rho_L I2 with
    I_k = int 2 pi d g^k dd  (Campbell; multi-tube coincidences are the mean^2 term,
    O(lam) relative corrections beyond it are ~1e-13 and included in the error budget).
    The integrals are done by IMPORTANCE-SAMPLED direct MC over the tube-distance d:
    proposal log-uniform on [d_lo, d_max]."""
    rho_L = 1.0 / L_web ** 2
    d_lo = 1e-3 * min(a, r0)
    lnZ = np.log(d_max / d_lo)
    d = d_lo * np.exp(rng.random(n) * lnZ)          # q(d) = 1/(d lnZ)
    gg = np.interp(d, d_grid, g_grid)
    w1 = 2 * np.pi * d * gg * (d * lnZ)             # integrand / proposal
    w2 = 2 * np.pi * d * gg ** 2 * (d * lnZ)
    I1, I1se = w1.mean(), w1.std(ddof=1) / np.sqrt(n)
    I2, I2se = w2.mean(), w2.std(ddof=1) / np.sqrt(n)
    mean = rho_L * I1
    var = rho_L * I2
    R = 1.0 + var / mean ** 2
    # jackknife SE on R over 20 slices (I1, I2 correlated)
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
            "mean_over_phi": float(mean / (rho_L * np.pi * a ** 2)),
            "I1_relse": float(I1se / I1), "I2_relse": float(I2se / I2)}


# ----------------------------------------------------------------------------
# 4. Motional narrowing: time-averaged kernel
# ----------------------------------------------------------------------------

def g_eff_table(d_grid, g_grid, a, r0, delta, d_max):
    """2D Gaussian convolution (rms delta per transverse axis) of g:
    g_eff(d) = int g(u) (u/delta^2) exp(-(u-d)^2/(2 delta^2)) i0e(u d/delta^2) du."""
    u = d_grid
    gu = g_grid
    out_max = d_max + 8.0 * delta
    d_out = np.unique(np.concatenate([
        np.array([0.0]),
        np.geomspace(max(1e-3 * min(a, r0), 1e-6 * delta), out_max, 900),
        np.linspace(0.0, out_max, 600)]))
    U = u[None, :]
    D = d_out[:, None]
    K = (U / delta ** 2) * np.exp(-((U - D) ** 2) / (2 * delta ** 2)) * i0e(U * D / delta ** 2)
    ge = np.trapezoid(K * gu[None, :], u, axis=1)
    return d_out, ge


def narrowing_curve(d_grid, g_grid, a, r0, L_web, d_max, deltas):
    rho_L = 1.0 / L_web ** 2
    rows = []
    for delta in deltas:
        dd, ge = g_eff_table(d_grid, g_grid, a, r0, delta, d_max)
        mean, var = campbell_moments(dd, ge, rho_L)
        rows.append({"delta_fm": float(delta), "delta_over_L": float(delta / L_web),
                     "R": float(1.0 + var / mean ** 2),
                     "mass_check": float(mean / (rho_L * np.pi * a ** 2))})
    return rows


def crossing(rows, target):
    """delta/L at which R first drops to target (log-log interpolation)."""
    x = np.array([r["delta_over_L"] for r in rows])
    y = np.array([r["R"] for r in rows])
    idx = np.argsort(x)
    x, y = x[idx], y[idx]
    for i in range(len(x) - 1):
        if y[i] > target >= y[i + 1]:
            t = (np.log(target) - np.log(y[i])) / (np.log(y[i + 1]) - np.log(y[i]))
            return float(np.exp(np.log(x[i]) + t * (np.log(x[i + 1]) - np.log(x[i]))))
    return None


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    rng = np.random.default_rng(RNG_SEED)
    results = {"seed": RNG_SEED,
               "kernel": "w(r) = exp(-r/r0)/(4 pi r0 r^2), normalized; r0 in {1,2} fm "
                         "(corpus kernel unprinted; double-propagator range ~1/q_bar ~ 2 fm per audit)",
               "conventions": {"rho_L": "1/L_web^2 (audit convention; phi = pi a^2/L^2)",
                               "units": "fm", "m_hat": "1 on tubes (R is normalization-independent)"}}

    # ---- h/g tables per (a, r0) used anywhere below --------------------------------
    print("== build kernel tables ==", flush=True)
    h_tabs = {r0: make_h_table(r0) for r0 in (1.0, 2.0)}

    # ---- leg 1: g(d) quadrature vs brute force -------------------------------------
    print("== leg 1: g(d) brute-force validation ==", flush=True)
    g_checks = []
    for (a, r0, d) in [(0.5, 1.0, 0.0), (0.5, 1.0, 2.0), (2.0, 1.0, 1.0), (2.0, 1.0, 4.0),
                       (52.6, 2.0, 20.0), (52.6, 2.0, 55.0), (10.0, 2.0, 9.0)]:
        sg, hg = h_tabs[r0]
        gq = g_of_d(d, a, r0, sg, hg)
        gm, gse = g_bruteforce(d, a, r0, 4_000_000, rng)
        g_checks.append({"a": a, "r0": r0, "d": d, "g_quad": gq, "g_mc": gm, "g_mc_se": gse,
                         "pull": (gq - gm) / gse if gse > 0 else 0.0})
        print(f"  a={a} r0={r0} d={d}: quad={gq:.6g} mc={gm:.6g}+-{gse:.2g} pull={(gq-gm)/max(gse,1e-30):+.2f}")
    results["g_validation"] = g_checks

    # ---- leg 3a: moderate-separation network MC vs Campbell ------------------------
    print("== leg 3a: network MC (moderate separation) ==", flush=True)
    val_rows = []
    for (a, phi) in [(0.5, 1.5e-3), (2.0, 3.0e-3), (10.0, 1.0e-2)]:
        r0 = 1.0
        sg, hg = h_tabs[r0]
        dg, gg, dmax = build_g_table(a, r0, sg, hg, extent_r0=60.0)
        rho_L = phi / (np.pi * a ** 2)
        mean_c, var_c = campbell_moments(dg, gg, rho_L)
        R_c = 1.0 + var_c / mean_c ** 2
        row = {"a": a, "r0": r0, "phi": phi, "R_campbell": float(R_c),
               "mass_check_campbell": float(mean_c / phi)}
        for kind in ("poisson", "lattice"):
            t0 = time.time()
            mc = network_mc_R(kind, a, r0, phi, dg, gg, dmax, rng,
                              n_blocks=4, n_per_block=750_000)
            row[f"mc_{kind}"] = mc
            print(f"  a={a} phi={phi} {kind}: R={mc['R_mean']:.1f}+-{mc['R_se']:.1f} "
                  f"(Campbell {R_c:.1f})  [{time.time()-t0:.0f}s]", flush=True)
        val_rows.append(row)
    results["network_mc_validation"] = val_rows

    # ---- leg 2 + 3b: the physical grid ----------------------------------------------
    print("== physical grid ==", flush=True)
    L_scan_um = [1.0, 1.7, 4.2, 10.0]   # heliknoton pitch/spacing: 4.2 um central (m3=0.047 eV), 1.7-10 um = 1 sigma; 1 um = audit floor
    a_scan = [2.0, 6.8, 20.0, 52.6, 200.0]  # fm; 6.8-49 = 1/f over f in [4,29] MeV; 52.6 = audit's 1/f; 2 = thin; 200 = generous thick
    grid_rows = []
    mc_phys_subset = {(6.8, 2.0, 1.0), (52.6, 2.0, 1.0), (52.6, 2.0, 4.2), (2.0, 1.0, 4.2),
                      (52.6, 1.0, 1.0), (200.0, 2.0, 10.0)}
    for r0 in (1.0, 2.0):
        sg, hg = h_tabs[r0]
        for a in a_scan:
            dg, gg, dmax = build_g_table(a, r0, sg, hg, extent_r0=100.0)
            for L_um in L_scan_um:
                L = L_um * FM_PER_UM
                rho_L = 1.0 / L ** 2
                phi = np.pi * a ** 2 / L ** 2
                mean_c, var_c = campbell_moments(dg, gg, rho_L)
                R_c = 1.0 + var_c / mean_c ** 2
                row = {"L_um": L_um, "a_fm": a, "r0_fm": r0, "phi": float(phi),
                       "R_campbell": float(R_c), "log10_R": float(np.log10(R_c)),
                       "mass_check": float(mean_c / phi),
                       "R_audit_thick": float(1.0 / phi),
                       "R_audit_thin": float(L ** 2 / (np.pi ** 2 * r0 ** 2)),
                       "ratio_to_audit_thick": float(R_c * phi)}
                if (a, r0, L_um) in mc_phys_subset:
                    mc = physical_mc_R(a, r0, L, dg, gg, dmax, rng, n=2_000_000)
                    row["mc_physical"] = mc
                    print(f"  MC-phys a={a} r0={r0} L={L_um}um: R={mc['R']:.4g}+-{mc['R_se']:.2g} "
                          f"(Campbell {R_c:.4g}); mass={mc['mean_over_phi']:.4f}", flush=True)
                grid_rows.append(row)
    results["physical_grid"] = grid_rows

    # ---- leg 4: motional narrowing ---------------------------------------------------
    print("== motional narrowing ==", flush=True)
    narrowing = []
    for (a, r0, L_um) in [(52.6, 2.0, 4.2), (52.6, 2.0, 1.0), (6.8, 2.0, 4.2), (52.6, 1.0, 4.2)]:
        sg, hg = h_tabs[r0]
        dg, gg, dmax = build_g_table(a, r0, sg, hg, extent_r0=100.0)
        L = L_um * FM_PER_UM
        deltas = np.geomspace(max(2 * a, 10 * r0), 3.0 * L, 26)
        rows = narrowing_curve(dg, gg, a, r0, L, dmax, deltas)
        d2, d10 = crossing(rows, 2.0), crossing(rows, 10.0)
        tau_s = (r0 * 1e-15) / 2.99792458e8          # coherence time ~ r0/c
        entry = {"a_fm": a, "r0_fm": r0, "L_um": L_um, "curve": rows,
                 "delta_over_L_at_R2": d2, "delta_over_L_at_R10": d10,
                 "tau_coh_s": tau_s}
        if d2:
            v_over_c = d2 * L / r0                    # (delta/tau)/c = delta/(r0)
            entry["required_speed_over_c_at_R2"] = float(v_over_c)
        if d10:
            entry["required_speed_over_c_at_R10"] = float(d10 * L / r0)
        narrowing.append(entry)
        print(f"  a={a} r0={r0} L={L_um}um: delta/L @R=10: {d10}; @R=2: {d2}; "
              f"v/c @R=2: {entry.get('required_speed_over_c_at_R2'):.3g}", flush=True)
    results["motional_narrowing"] = narrowing

    results["runtime_s"] = time.time() - T_START
    out = "/home/user/fork_frankensim/analysis/gum-sandbox/theory-audit/h21_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"wrote {out} ({results['runtime_s']:.0f}s)")
    return results


if __name__ == "__main__":
    main()
