#!/usr/bin/env python3
"""
Tier-7 M4: Nelson non-equilibrium relaxation tau(nu) — F-T7-M4.

Extends tier3-born (analysis/gum-sandbox/tier3-born/born.py, REPORT.md) from
deterministic de Broglie-Bohm guidance to Nelson stochastic mechanics:

    dX = (v + u) dt + sqrt(2 nu) dW,
    v = Im(grad psi / psi)   (dBB current velocity),
    u = nu * grad ln|psi|^2  (osmotic velocity),

with nu in {0 (dBB control), 0.1, 0.5, 1.0}; hbar = m = 1 (GUM Sec. III has
nu = hbar/2m = 1/2 in these units; we scan through it).

Conventions REUSED from tier3 verbatim:
  * 2-D box [0, pi]^2, psi = M^{-1/2} sum e^{i theta_mn} phi_mn e^{-i E_mn t},
    equal weights, phases seeded 42; M = 9 (kmax 3) and M = 16 (kmax 4).
  * rho_0 = |phi_11|^2 (non-equilibrium start), rejection-sampled with
    seed 12345 + M (tier3's SEED_SAMPLE + m convention).
  * Coarse-grained H(t) = sum Pbar ln(Pbar/Qbar) at 32x32 (headline) and
    16x16, Qbar from box-averaging |psi|^2 on a 256^2 fine grid.
  * Exponential fit window: up to where H first drops to 10% of H(0),
    or t = 2*pi, whichever is first (tier3 fit_tau, unchanged).
  * dBB control integrator: tier3's RK4 with per-particle substepping and
    displacement cap, dt = 2e-3.

New for nu > 0:
  * Euler-Maruyama with the same per-particle substep guard on the DRIFT
    displacement (|v+u| h <= 0.05 per substep, capped), noise sqrt(2 nu h)
    per substep, dt = 1e-3 (halving check at 5e-4).
  * Reflecting walls: positions folded back into the box after each substep
    (z -> fold(z mod 2L)); the osmotic term ~ 2 nu cot(x) already repels
    from the walls where |psi|^2 -> 0. No-flux is verified numerically:
    (a) reflection-event rate logged, (b) all particles remain in the box,
    (c) the equilibrium control's wall-band occupancy matches the Born value,
    (d) the equilibrium control H stays at the finite-N floor.

N = 10,000 particles (roadmap spec; tier3 used 20,000 — noise floors double).
Bootstrap (200 resamples of particle trajectories) gives tau errors.

Outputs (next to this script): m4_results.json, m4_hdata.csv, m4_relax.png
"""

import json
import os
import sys as _sys
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Parameters (tier3 values kept where the convention is reused)
# ----------------------------------------------------------------------------
OUTDIR = os.path.dirname(os.path.abspath(__file__))

L = np.pi                     # box side
DT_DBB = 2e-3                 # tier3 RK4 step (nu = 0 control)
DT_SDE = 1e-3                 # Euler-Maruyama base step (nu > 0)
T_FINAL = 4.0 * np.pi         # tier3 horizon (= exact psi revival period)
N_PART = int(os.environ.get("M4_N", 10000))   # roadmap spec (tier3: 20,000)
DISP_MAX = 0.05               # max |drift| * h per substep (tier3 value)
MAX_SUB = 32                  # cap on per-particle substeps (tier3 value)
EPS_WALL = 1e-4               # dBB clip guard (tier3); SDE uses 1e-12 nudge
EPS_REFL = 1e-12              # post-reflection nudge off the exact wall
NFINE = 256                   # fine grid for |psi|^2 cell integrals (tier3)
CGS = (32, 16)                # coarse-graining levels (tier3)
SEED_PHASES = 42              # tier3 mode phases
SEED_SAMPLE = 12345           # tier3 ensemble sampling base
SEED_NOISE = 20260718         # Wiener noise base (new; deterministic per run)
NU_LIST = (0.0, 0.1, 0.5, 1.0)
NU_DIAG = (0.02, 0.05)        # small-nu diagnostic runs (continuity map)
KMAX_LIST = (3, 4)            # M = 9, 16 (roadmap: skip M = 4)
N_BOOT = 200                  # bootstrap resamples for tau errors
TWO_PI = 2.0 * np.pi

SMOKE = "--smoke" in _sys.argv  # short timing run


def dt_out_for(nu):
    """Diagnostic cadence: tier3's 0.1 for the dBB control; finer for noisy
    runs whose tau is O(0.1-1) (resolution refinement, not a physics change)."""
    if nu == 0.0:
        return 0.1
    return 0.05 if nu < 0.5 else 0.02


# ----------------------------------------------------------------------------
# Mode system (tier3 ModeSystem + osmotic drift)
# ----------------------------------------------------------------------------
class ModeSystem:
    """Superposition of the first kmax x kmax box eigenmodes, equal weights,
    fixed random phases (seed 42). Identical to tier3; drift() added."""

    def __init__(self, kmax):
        self.kmax = kmax
        ks = np.arange(1, kmax + 1)
        mm, nn = np.meshgrid(ks, ks, indexing="ij")
        self.mm = mm.ravel().astype(float)
        self.nn = nn.ravel().astype(float)
        self.M = kmax * kmax
        self.E = 0.5 * (self.mm ** 2 + self.nn ** 2)
        rng = np.random.default_rng(SEED_PHASES)
        self.theta = rng.uniform(0.0, 2.0 * np.pi, self.M)
        self.ks = ks.astype(float)

    def coeff_matrix(self, t):
        c = np.exp(1j * (self.theta - self.E * t)) / np.sqrt(self.M)
        return c.reshape(self.kmax, self.kmax)

    def _fields(self, x, y, t):
        """psi and grad psi contractions at particle positions (norm
        prefactor omitted: cancels in v and u)."""
        C = self.coeff_matrix(t)
        Cr = np.ascontiguousarray(C.real)
        Ci = np.ascontiguousarray(C.imag)
        xk = x[:, None] * self.ks
        yk = y[:, None] * self.ks
        sx = np.sin(xk)
        cxm = np.cos(xk) * self.ks
        sy = np.sin(yk)
        cym = np.cos(yk) * self.ks
        ar, ai = sx @ Cr, sx @ Ci
        pr = np.einsum("ij,ij->i", ar, sy)
        pi_ = np.einsum("ij,ij->i", ai, sy)
        gxr = np.einsum("ij,ij->i", cxm @ Cr, sy)
        gxi = np.einsum("ij,ij->i", cxm @ Ci, sy)
        gyr = np.einsum("ij,ij->i", ar, cym)
        gyi = np.einsum("ij,ij->i", ai, cym)
        dens = np.maximum(pr * pr + pi_ * pi_, 1e-30)
        return pr, pi_, gxr, gxi, gyr, gyi, dens

    def velocity(self, x, y, t):
        """dBB velocity v = Im(grad psi / psi) (tier3)."""
        pr, pi_, gxr, gxi, gyr, gyi, dens = self._fields(x, y, t)
        vx = (gxi * pr - gxr * pi_) / dens
        vy = (gyi * pr - gyr * pi_) / dens
        return vx, vy

    def drift(self, x, y, t, nu):
        """Nelson forward drift b = v + u,  u = nu grad ln|psi|^2
        = 2 nu Re(grad psi / psi)."""
        pr, pi_, gxr, gxi, gyr, gyi, dens = self._fields(x, y, t)
        vx = (gxi * pr - gxr * pi_) / dens
        vy = (gyi * pr - gyr * pi_) / dens
        ux = 2.0 * nu * (gxr * pr + gxi * pi_) / dens
        uy = 2.0 * nu * (gyr * pr + gyi * pi_) / dens
        return vx + ux, vy + uy

    def psi2_grid(self, t, n=NFINE):
        xg = (np.arange(n) + 0.5) * L / n
        sxg = np.sin(np.outer(xg, self.ks))
        C = self.coeff_matrix(t) * (2.0 / np.pi)
        psi = sxg @ C @ sxg.T
        return psi.real ** 2 + psi.imag ** 2

    def psi2_points(self, x, y, t):
        C = self.coeff_matrix(t) * (2.0 / np.pi)
        sx = np.sin(x[:, None] * self.ks)
        sy = np.sin(y[:, None] * self.ks)
        pr = np.einsum("ij,ij->i", sx @ np.ascontiguousarray(C.real), sy)
        pi_ = np.einsum("ij,ij->i", sx @ np.ascontiguousarray(C.imag), sy)
        return pr ** 2 + pi_ ** 2


# ----------------------------------------------------------------------------
# Sampling (tier3, verbatim)
# ----------------------------------------------------------------------------
def sample_phi11(n, rng):
    """rho_0 = |phi_11|^2 (tier3 non-equilibrium initial ensemble)."""
    xs, ys, got = [], [], 0
    while got < n:
        prop = rng.uniform(0.0, L, size=(2 * n, 2))
        u = rng.uniform(0.0, 1.0, 2 * n)
        acc = u < np.sin(prop[:, 0]) ** 2 * np.sin(prop[:, 1]) ** 2
        xs.append(prop[acc, 0])
        ys.append(prop[acc, 1])
        got += int(acc.sum())
    x = np.concatenate(xs)[:n]
    y = np.concatenate(ys)[:n]
    return np.clip(x, EPS_WALL, L - EPS_WALL), np.clip(y, EPS_WALL, L - EPS_WALL)


def sample_psi2(sysm, n, rng):
    """Equilibrium control ensemble from |psi(0)|^2 (tier3)."""
    bound = 1.2 * sysm.psi2_grid(0.0, n=512).max()
    xs, ys, got = [], [], 0
    while got < n:
        prop = rng.uniform(0.0, L, size=(2 * n, 2))
        u = rng.uniform(0.0, bound, 2 * n)
        acc = u < sysm.psi2_points(prop[:, 0], prop[:, 1], 0.0)
        xs.append(prop[acc, 0])
        ys.append(prop[acc, 1])
        got += int(acc.sum())
    x = np.concatenate(xs)[:n]
    y = np.concatenate(ys)[:n]
    return np.clip(x, EPS_WALL, L - EPS_WALL), np.clip(y, EPS_WALL, L - EPS_WALL)


# ----------------------------------------------------------------------------
# Integrators
# ----------------------------------------------------------------------------
def _rk4_substep(sysm, x, y, t, h, k1=None):
    """tier3 RK4 substep with displacement cap and wall clip (nu = 0)."""
    k1x, k1y = sysm.velocity(x, y, t) if k1 is None else k1
    k2x, k2y = sysm.velocity(x + 0.5 * h * k1x, y + 0.5 * h * k1y, t + 0.5 * h)
    k3x, k3y = sysm.velocity(x + 0.5 * h * k2x, y + 0.5 * h * k2y, t + 0.5 * h)
    k4x, k4y = sysm.velocity(x + h * k3x, y + h * k3y, t + h)
    dx = h * (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0
    dy = h * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0
    disp = np.hypot(dx, dy)
    big = disp > DISP_MAX
    if big.any():
        f = DISP_MAX / disp[big]
        dx[big] *= f
        dy[big] *= f
    x = np.clip(x + dx, EPS_WALL, L - EPS_WALL)
    y = np.clip(y + dy, EPS_WALL, L - EPS_WALL)
    return x, y


def macro_step_dbb(sysm, x, y, t, dt):
    """tier3 macro step (RK4 + adaptive per-particle substepping)."""
    vx, vy = sysm.velocity(x, y, t)
    speed = np.hypot(vx, vy)
    nsub = np.clip(np.ceil(speed * dt / DISP_MAX).astype(int), 1, MAX_SUB)
    uniq = np.unique(nsub)
    if uniq.size == 1 and uniq[0] == 1:
        return _rk4_substep(sysm, x, y, t, dt, k1=(vx, vy))
    for n in uniq:
        idx = np.where(nsub == n)[0]
        xs, ys = x[idx], y[idx]
        h = dt / n
        k1 = (vx[idx], vy[idx])
        for j in range(n):
            xs, ys = _rk4_substep(sysm, xs, ys, t + j * h, h, k1=k1)
            k1 = None
        x[idx], y[idx] = xs, ys
    return x, y


def _reflect(z, stats):
    """Reflecting boundary: fold back into [0, L]; count wall events."""
    out = (z < 0.0) | (z > L)
    n_out = int(out.sum())
    if n_out:
        stats["reflections"] += n_out
        z = np.mod(z, 2.0 * L)
        z = np.where(z > L, 2.0 * L - z, z)
    return np.clip(z, EPS_REFL, L - EPS_REFL)


def macro_step_sde(sysm, x, y, t, dt, nu, rng, stats):
    """One Euler-Maruyama macro step with per-particle drift substepping.

    Substep count from the initial drift speed (|b| h <= DISP_MAX), drift
    displacement hard-capped per substep (node/wall spikes), Gaussian noise
    sqrt(2 nu h) per substep, reflection after every substep."""
    bx, by = sysm.drift(x, y, t, nu)
    speed = np.hypot(bx, by)
    nsub = np.clip(np.ceil(speed * dt / DISP_MAX).astype(int), 1, MAX_SUB)
    for n in np.unique(nsub):
        idx = np.where(nsub == n)[0]
        xs, ys = x[idx], y[idx]
        h = dt / n
        sig = np.sqrt(2.0 * nu * h)
        cbx, cby = bx[idx], by[idx]
        for j in range(n):
            if j > 0:
                cbx, cby = sysm.drift(xs, ys, t + j * h, nu)
            dx = h * cbx
            dy = h * cby
            disp = np.hypot(dx, dy)
            big = disp > DISP_MAX
            if big.any():
                f = DISP_MAX / disp[big]
                dx[big] *= f
                dy[big] *= f
            xs = _reflect(xs + dx + sig * rng.standard_normal(xs.size), stats)
            ys = _reflect(ys + dy + sig * rng.standard_normal(ys.size), stats)
        x[idx], y[idx] = xs, ys
    return x, y


# ----------------------------------------------------------------------------
# Coarse-grained H-function (tier3, verbatim) + bootstrap helper
# ----------------------------------------------------------------------------
def coarse_q(q_fine, ncell):
    b = NFINE // ncell
    q = q_fine.reshape(ncell, b, ncell, b).mean(axis=(1, 3))
    return q / q.sum()


def h_function(x, y, q_fine, ncell):
    counts, _, _ = np.histogram2d(x, y, bins=ncell, range=[[0.0, L], [0.0, L]])
    p = counts / counts.sum()
    q = coarse_q(q_fine, ncell)
    mask = p > 0
    qm = np.maximum(q[mask], 1e-300)
    return float(np.sum(p[mask] * np.log(p[mask] / qm)))


def h_from_counts(counts, qbar):
    p = counts / counts.sum()
    mask = p > 0
    qm = np.maximum(qbar[mask], 1e-300)
    return float(np.sum(p[mask] * np.log(p[mask] / qm)))


def fit_tau(ts, hs, floor=0.0):
    """tier3 fit: window up to H <= 0.1 H(0) or t = 2 pi, whichever first.

    One guard added for the fast noisy runs: the window also ends where H
    first drops below 2x the finite-N noise floor (floor = ncells/2N), so
    the exponential is never fitted through the floor. For the nu = 0 dBB
    runs this guard is inactive (H stays above 2x floor for t <= 2 pi), so
    the tier3 convention is exactly preserved where G1 compares to tier3."""
    h0 = hs[0]
    thresh = max(0.1 * h0, 2.0 * floor)
    below = np.where(hs <= thresh)[0]
    t_end = ts[below[0]] if below.size else TWO_PI
    t_end = min(t_end, TWO_PI)
    sel = (ts <= t_end) & (hs > 0)
    if sel.sum() < 4:
        return np.nan, np.nan, float(t_end), np.nan
    slope, intercept = np.polyfit(ts[sel], np.log(hs[sel]), 1)
    tau = -1.0 / slope
    resid = np.log(hs[sel]) - (slope * ts[sel] + intercept)
    ss_tot = np.sum((np.log(hs[sel]) - np.log(hs[sel]).mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan
    return float(tau), float(np.exp(intercept)), float(t_end), float(r2)


# ----------------------------------------------------------------------------
# Run driver
# ----------------------------------------------------------------------------
def run_ensemble(sysm, x, y, nu, dt, label, noise_seed,
                 store_positions=False, t_final=T_FINAL, dt_out=None):
    """Evolve to t_final, sample H every dt_out (tier3 cadence for nu = 0;
    last sample exactly at t_final). Positions stored only for t <= 2 pi
    (the bootstrap window). Returns ts, {cg: H}, positions, stats."""
    t0 = time.time()
    if dt_out is None:
        dt_out = dt_out_for(nu)
    rng = np.random.default_rng(noise_seed)
    stats = {"reflections": 0}
    n_full = int(np.floor(t_final / dt_out))
    steps_per_out = int(round(dt_out / dt))
    ts, h_by_cg = [], {cg: [] for cg in CGS}
    positions = []

    def record(t):
        qf = sysm.psi2_grid(t)
        ts.append(t)
        for cg in CGS:
            h_by_cg[cg].append(h_function(x, y, qf, cg))
        if store_positions and t <= TWO_PI + 1e-9:
            positions.append((x.copy(), y.copy()))

    t = 0.0
    record(t)
    for _k in range(1, n_full + 1):
        for _ in range(steps_per_out):
            if nu == 0.0:
                x, y = macro_step_dbb(sysm, x, y, t, dt)
            else:
                x, y = macro_step_sde(sysm, x, y, t, dt, nu, rng, stats)
            t += dt
        record(t)
    rem = t_final - t
    if rem > 1e-12:
        n_extra = max(1, int(np.ceil(rem / dt)))
        h = rem / n_extra
        for _ in range(n_extra):
            if nu == 0.0:
                x, y = macro_step_dbb(sysm, x, y, t, h)
            else:
                x, y = macro_step_sde(sysm, x, y, t, h, nu, rng, stats)
            t += h
        t = t_final
        record(t)
    stats["in_box"] = bool(np.all((x >= 0) & (x <= L) & (y >= 0) & (y <= L)))
    stats["runtime_s"] = time.time() - t0
    stats["final_xy"] = (x, y)
    print(f"[{label}] {stats['runtime_s']:.1f} s  H32(final) = "
          f"{h_by_cg[32][-1]:.4f}  reflections = {stats['reflections']}",
          flush=True)
    return (np.array(ts), {cg: np.array(v) for cg, v in h_by_cg.items()},
            positions, stats)


def bootstrap_taus(positions, ts_w, qbars_w, n_boot, seed, floor):
    """Bootstrap tau (32^2 grain): resample particle trajectories with
    replacement, recompute H(t) over the stored window (t <= 2 pi), refit."""
    rng = np.random.default_rng(seed)
    n = positions[0][0].size
    edges = np.linspace(0.0, L, 33)
    taus = np.empty(n_boot)
    for b in range(n_boot):
        ii = rng.integers(0, n, n)
        hs = np.empty(len(ts_w))
        for k, (px, py) in enumerate(positions):
            counts, _, _ = np.histogram2d(px[ii], py[ii], bins=[edges, edges])
            hs[k] = h_from_counts(counts, qbars_w[k])
        taus[b], _, _, _ = fit_tau(np.asarray(ts_w), hs, floor=floor)
    return taus


def wall_band_check(sysm, x, y, t, band=L / 32):
    """No-flux check: fraction of particles within `band` of any wall vs the
    Born (|psi(t)|^2) value computed on the fine grid."""
    d = np.minimum(np.minimum(x, L - x), np.minimum(y, L - y))
    frac = float((d < band).mean())
    q = sysm.psi2_grid(t)
    q = q / q.sum()
    n = q.shape[0]
    xg = (np.arange(n) + 0.5) * L / n
    dg = np.minimum(xg, L - xg)
    in_band = (dg[:, None] < band) | (dg[None, :] < band)
    born = float(q[in_band].sum())
    return frac, born


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
FLOORS = {cg: cg * cg / (2.0 * N_PART) for cg in CGS}
CACHE_DIR = os.path.join(OUTDIR, "cache")


def _cache_path(key):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, key + ".json")


def _cache_load(key):
    p = _cache_path(key)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return None


def _cache_save(key, obj):
    with open(_cache_path(key), "w") as f:
        json.dump(obj, f)


def nonequilibrium_run(sysm, M, nu, dt, nseed, t_final, n_boot,
                       dt_out=None, do_boot=True):
    """One non-equilibrium run from tier3's rho_0 = |phi_11|^2 (seeded
    12345 + M): evolve, floor-guarded tier3 fit at both grains, bootstrap.
    Cached to disk (seeded, so a cache hit is bit-identical work)."""
    key = (f"run_N{N_PART}_M{M}_nu{nu:g}_dt{dt:g}_s{nseed}_T{t_final:g}"
           f"_b{n_boot if do_boot else 0}")
    hit = _cache_load(key)
    if hit is not None:
        print(f"[cache hit] {key}", flush=True)
        return (hit["rr"], np.array(hit["ts"]),
                {32: np.array(hit["h32"]), 16: np.array(hit["h16"])},
                np.array(hit["taus_b"]) if hit["taus_b"] is not None
                else None)
    x0, y0 = sample_phi11(N_PART, np.random.default_rng(SEED_SAMPLE + M))
    ts, hcg, pos, stats = run_ensemble(
        sysm, x0, y0, nu, dt, f"M={M} nu={nu} dt={dt:g}", nseed,
        store_positions=do_boot, t_final=t_final, dt_out=dt_out)
    rr = {"nu": nu, "dt": dt, "dt_out": dt_out or dt_out_for(nu),
          "reflections_per_particle": stats["reflections"] / N_PART,
          "all_in_box": stats["in_box"], "runtime_s": stats["runtime_s"]}
    for cg in CGS:
        tau, h0f, t_end, r2 = fit_tau(ts, hcg[cg], floor=FLOORS[cg])
        rr[f"cg{cg}"] = {"H_initial": float(hcg[cg][0]),
                         "H_final": float(hcg[cg][-1]),
                         "tau": tau, "H0_fit": h0f,
                         "fit_window_t_end": t_end, "fit_r2": r2}
    taus_b = None
    if do_boot:
        nw = len(pos)
        ts_w = ts[:nw]
        qb_w = [coarse_q(sysm.psi2_grid(t), 32) for t in ts_w]
        taus_b = bootstrap_taus(pos, ts_w, qb_w, n_boot,
                                seed=nseed + 5000, floor=FLOORS[32])
        rr["tau32_boot_mean"] = float(np.nanmean(taus_b))
        rr["tau32_boot_std"] = float(np.nanstd(taus_b))
    _cache_save(key, {"rr": rr, "ts": list(map(float, ts)),
                      "h32": list(map(float, hcg[32])),
                      "h16": list(map(float, hcg[16])),
                      "taus_b": (list(map(float, taus_b))
                                 if taus_b is not None else None)})
    return rr, ts, hcg, taus_b


def main():
    t_start = time.time()
    t_final = 2.0 if SMOKE else T_FINAL
    n_boot = 10 if SMOKE else N_BOOT
    results = {"N": N_PART, "dt_dbb": DT_DBB, "dt_sde": DT_SDE,
               "t_final": float(t_final), "seed_phases": SEED_PHASES,
               "seed_sample": SEED_SAMPLE, "seed_noise": SEED_NOISE,
               "nu_list": list(NU_LIST), "nu_diag": list(NU_DIAG),
               "n_boot": n_boot,
               "noise_floor": {f"cg{cg}": FLOORS[cg] for cg in CGS},
               "tier3_reference": {"M9_tau32": 10.201107001665816,
                                   "M16_tau32": 4.852765879922653},
               "runs": {}, "controls": {}, "dt_halving": {}}

    h_store = {}      # (M, nu) -> (ts, H32, H16); ("ctrl", M, nu) -> (ts, H32)
    boot_store = {}   # (M, nu) -> bootstrap tau array

    for kmax in KMAX_LIST:
        sysm = ModeSystem(kmax)
        M = sysm.M
        all_nus = list(NU_LIST) + list(NU_DIAG)
        for inu, nu in enumerate(all_nus):
            is_diag = nu in NU_DIAG
            dt = DT_DBB if nu == 0.0 else DT_SDE
            nseed = SEED_NOISE + 100 * M + inu
            rr, ts, hcg, taus_b = nonequilibrium_run(
                sysm, M, nu, dt, nseed, t_final, n_boot)
            rr["diagnostic"] = is_diag
            boot_store[(M, nu)] = taus_b
            if taus_b is not None:
                # Error model: the bootstrap tau distribution's LOCATION is
                # biased upward (resampling with replacement doubles the
                # multinomial noise, inflating the effective floor — worst
                # for fast decays), so we quote the RELATIVE spread applied
                # to the point estimate; cross-checked against independent
                # seed-replicate scatter (conservative there).
                rel = float(np.nanstd(taus_b) / np.nanmean(taus_b))
                rr["tau32_boot_rel_spread"] = rel
                rr["tau32_err"] = rr["cg32"]["tau"] * rel
            results["runs"][f"M={M},nu={nu}"] = rr
            h_store[(M, nu)] = (ts, hcg[32], hcg[16])

            if is_diag:
                continue   # diagnostics: no halving check, no control

            # --- dt-halving convergence check ---
            rrh, _, _, _ = nonequilibrium_run(
                sysm, M, nu, dt / 2.0, nseed + 9000, t_final, n_boot,
                dt_out=dt_out_for(nu), do_boot=False)
            tau_full, tau_h = rr["cg32"]["tau"], rrh["cg32"]["tau"]
            hv = {"tau_dt": tau_full, "tau_dt_half": tau_h,
                  "rel_change": abs(tau_h - tau_full) / tau_full}
            if hv["rel_change"] > 0.05 and not SMOKE:
                # stochastic scatter vs discretization bias: average taus
                # over 3 independent noise seeds at each dt level
                t_a, t_b = [tau_full], [tau_h]
                for k in range(2):
                    ra, _, _, _ = nonequilibrium_run(
                        sysm, M, nu, dt, nseed + 11000 + k, t_final,
                        n_boot, do_boot=False)
                    rb, _, _, _ = nonequilibrium_run(
                        sysm, M, nu, dt / 2.0, nseed + 13000 + k, t_final,
                        n_boot, dt_out=dt_out_for(nu), do_boot=False)
                    t_a.append(ra["cg32"]["tau"])
                    t_b.append(rb["cg32"]["tau"])
                hv["seed_avg_tau_dt"] = float(np.mean(t_a))
                hv["seed_avg_tau_dt_half"] = float(np.mean(t_b))
                hv["seed_taus_dt"] = t_a
                hv["seed_taus_dt_half"] = t_b
                hv["rel_change_seed_avg"] = float(
                    abs(np.mean(t_b) - np.mean(t_a)) / np.mean(t_a))
            results["dt_halving"][f"M={M},nu={nu}"] = hv

            # --- equivariance control at every nu, both M (cached) ---
            ckey = f"ctrl_N{N_PART}_M{M}_nu{nu:g}_s{nseed + 7000}_T{t_final:g}"
            chit = _cache_load(ckey)
            if chit is not None:
                print(f"[cache hit] {ckey}", flush=True)
                results["controls"][f"M={M},nu={nu}"] = chit["cc"]
                h_store[("ctrl", M, nu)] = (np.array(chit["ts"]),
                                            np.array(chit["h32"]))
            else:
                rngc = np.random.default_rng(SEED_SAMPLE + 999)  # tier3 seed
                xc, yc = sample_psi2(sysm, N_PART, rngc)
                tsc, hcgc, _, statsc = run_ensemble(
                    sysm, xc, yc, nu, dt, f"ctrl M={M} nu={nu}",
                    nseed + 7000, store_positions=False, t_final=t_final)
                xcf, ycf = statsc["final_xy"]
                frac, born = wall_band_check(sysm, xcf, ycf, t_final)
                cc = {"nu": nu,
                      "wall_band_frac": frac, "wall_band_born": born,
                      "all_in_box": statsc["in_box"],
                      "reflections_per_particle":
                          statsc["reflections"] / N_PART}
                for cg in CGS:
                    cc[f"cg{cg}"] = {"H_initial": float(hcgc[cg][0]),
                                     "H_final": float(hcgc[cg][-1]),
                                     "H_max": float(hcgc[cg].max()),
                                     "H_mean": float(hcgc[cg].mean())}
                _cache_save(ckey, {"cc": cc, "ts": list(map(float, tsc)),
                                   "h32": list(map(float, hcgc[32]))})
                results["controls"][f"M={M},nu={nu}"] = cc
                h_store[("ctrl", M, nu)] = (tsc, hcgc[32])

    # ------------------------------------------------------------------
    # Gate arithmetic
    # ------------------------------------------------------------------
    g = {}
    # G1: dBB control vs tier3
    for M, ref in ((9, results["tier3_reference"]["M9_tau32"]),
                   (16, results["tier3_reference"]["M16_tau32"])):
        tau0 = results["runs"][f"M={M},nu=0.0"]["cg32"]["tau"]
        g[f"G1_M{M}"] = {"tau_dbb": tau0, "tier3_tau": ref,
                         "rel_dev": abs(tau0 - ref) / ref,
                         "pass": bool(abs(tau0 - ref) / ref <= 0.25)}
    # G2: controls at floor (criterion: H_max < 2x tier3-style floor,
    # floor = ncells/2N; tier3's control peaked at 1.18x floor)
    floor32 = results["noise_floor"]["cg32"]
    g2_pass = True
    g2_rows = {}
    for key, cc in results["controls"].items():
        ratio = cc["cg32"]["H_max"] / floor32
        ok = ratio < 2.0
        g2_rows[key] = {"H_max": cc["cg32"]["H_max"],
                        "Hmax_over_floor": ratio, "pass": bool(ok)}
        g2_pass = g2_pass and ok
    g["G2"] = {"floor32": floor32, "rows": g2_rows, "pass": bool(g2_pass)}
    # G3: strict monotonicity of tau(nu) for both M (pre-registered grid);
    # supplementary: with the small-nu diagnostic points inserted
    g3 = {}
    for kmax in KMAX_LIST:
        M = kmax * kmax
        taus = [results["runs"][f"M={M},nu={nu}"]["cg32"]["tau"]
                for nu in NU_LIST]
        dec = all(taus[i + 1] < taus[i] for i in range(len(taus) - 1))
        nus_full = sorted(list(NU_LIST) + list(NU_DIAG))
        taus_full = [results["runs"][f"M={M},nu={nu}"]["cg32"]["tau"]
                     for nu in nus_full]
        dec_full = all(taus_full[i + 1] < taus_full[i]
                       for i in range(len(taus_full) - 1))
        g3[f"M{M}"] = {"nus": list(NU_LIST), "taus": taus,
                       "strictly_decreasing": bool(dec),
                       "nus_with_diag": nus_full,
                       "taus_with_diag": taus_full,
                       "strictly_decreasing_with_diag": bool(dec_full)}
    g["G3"] = {**g3, "pass": bool(all(v["strictly_decreasing"]
                                      for v in g3.values()))}
    # G4: tau(0.5)/tau(0) with bootstrap errors; continuity tau(0.1)/tau(0)
    # >= 0.7 (pre-registered); small-nu diagnostic ratios for the approach
    g4 = {}

    def boot_ratio_err(M, nu_a, nu_b, ratio):
        """Delta-method ratio error from bootstrap RELATIVE spreads
        (independent runs; see error-model note above)."""
        ta, tb = boot_store[(M, nu_a)], boot_store[(M, nu_b)]
        ra = np.nanstd(ta) / np.nanmean(ta)
        rb = np.nanstd(tb) / np.nanmean(tb)
        return float(ratio * np.sqrt(ra * ra + rb * rb))

    for kmax in KMAX_LIST:
        M = kmax * kmax
        tau = {nu: results["runs"][f"M={M},nu={nu}"]["cg32"]["tau"]
               for nu in list(NU_LIST) + list(NU_DIAG)}
        r05 = tau[0.5] / tau[0.0]
        r01 = tau[0.1] / tau[0.0]
        g4[f"M{M}"] = {
            "ratio_05_over_0": r05,
            "ratio_05_err": boot_ratio_err(M, 0.5, 0.0, r05),
            "ratio_01_over_0": r01,
            "ratio_01_err": boot_ratio_err(M, 0.1, 0.0, r01),
            "diag_ratios": {str(nu): tau[nu] / tau[0.0] for nu in NU_DIAG},
            "continuity_pass": bool(r01 >= 0.7)}
    g["G4"] = {**g4, "pass": bool(all(v["continuity_pass"]
                                      for v in g4.values()))}
    # dt-halving summary (seed-averaged value supersedes single-seed where
    # the escalation ran)
    eff = {k: v.get("rel_change_seed_avg", v["rel_change"])
           for k, v in results["dt_halving"].items()}
    worst_key = max(eff, key=eff.get)
    g["dt_convergence"] = {"worst_config": worst_key,
                           "worst_rel_change": eff[worst_key],
                           "per_config": eff,
                           "pass": bool(eff[worst_key] < 0.05)}
    results["gates"] = g

    # ------------------------------------------------------------------
    # m4_hdata.csv (long format: cadence differs across runs)
    # ------------------------------------------------------------------
    with open(os.path.join(OUTDIR, "m4_hdata.csv"), "w") as f:
        f.write("kind,M,nu,t,H32,H16\n")
        for key, val in h_store.items():
            if key[0] == "ctrl":
                _, M, nu = key
                ts_k, h32 = val
                for t, a in zip(ts_k, h32):
                    f.write(f"ctrl,{M},{nu},{t:.6g},{a:.8g},\n")
            else:
                M, nu = key
                ts_k, h32, h16 = val
                for t, a, b in zip(ts_k, h32, h16):
                    f.write(f"noneq,{M},{nu},{t:.6g},{a:.8g},{b:.8g}\n")

    # ------------------------------------------------------------------
    # Figure: H(t) per M (ordinal blue ramp over nu) + tau(nu) panel
    # ------------------------------------------------------------------
    nu_col = {0.0: "#86b6ef", 0.1: "#3987e5", 0.5: "#1c5cab", 1.0: "#0d366b"}
    m_col = {9: "#eda100", 16: "#e34948"}      # tier3 M colors
    gray = "#898781"
    ink, muted = "#0b0b0b", "#52514e"

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
    for ax, kmax in zip(axes[:2], KMAX_LIST):
        M = kmax * kmax
        for nu in NU_LIST:
            lab = ("$\\nu = 0$ (dBB)" if nu == 0.0 else f"$\\nu$ = {nu}")
            ts_k, h32, _h16 = h_store[(M, nu)]
            ax.plot(ts_k, np.maximum(h32, 1e-6),
                    color=nu_col[nu], lw=1.8, label=lab)
            rr = results["runs"][f"M={M},nu={nu}"]["cg32"]
            tt = np.linspace(0, rr["fit_window_t_end"], 50)
            ax.plot(tt, rr["H0_fit"] * np.exp(-tt / rr["tau"]),
                    color=nu_col[nu], lw=1.0, ls=":", alpha=0.9)
        ts_c, h32_c = h_store[("ctrl", M, 1.0)]
        ax.plot(ts_c, np.maximum(h32_c, 1e-6),
                color=gray, lw=1.4, ls="--",
                label="equilib. control ($\\nu$ = 1)")
        ax.axhline(floor32, color=muted, lw=1.0, ls="-.",
                   label=f"noise floor $\\approx$ {floor32:.3g}")
        ax.set_yscale("log")
        ax.set_xlabel("t", color=ink)
        ax.set_title(f"M = {M}", fontsize=11, color=ink)
        ax.grid(True, which="major", color="#DDDDDD", lw=0.6)
        ax.tick_params(colors=muted)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel(r"$\bar{H}(t)=\sum \bar{P}\,\ln(\bar{P}/\bar{Q})$"
                       "  (32$\\times$32)", color=ink)
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")

    ax = axes[2]
    nus_full = sorted(list(NU_LIST) + list(NU_DIAG))
    for kmax in KMAX_LIST:
        M = kmax * kmax
        taus = [results["runs"][f"M={M},nu={nu}"]["cg32"]["tau"]
                for nu in nus_full]
        errs = [results["runs"][f"M={M},nu={nu}"]["tau32_err"]
                for nu in nus_full]
        ax.errorbar(nus_full, taus, yerr=errs, color=m_col[M], lw=1.8,
                    marker="o", ms=5, capsize=3, label=f"M = {M}")
        ax.annotate(f"M = {M}", xy=(nus_full[2], taus[2]),
                    xytext=(6, 6), textcoords="offset points",
                    fontsize=9, color=muted)
    ax.axhline(results["tier3_reference"]["M9_tau32"], color=m_col[9],
               lw=0.9, ls="--", alpha=0.6)
    ax.axhline(results["tier3_reference"]["M16_tau32"], color=m_col[16],
               lw=0.9, ls="--", alpha=0.6)
    ax.text(0.62, results["tier3_reference"]["M9_tau32"] * 1.04,
            "tier3 dBB (N = 20k)", fontsize=8, color=muted)
    ax.set_yscale("log")
    ax.set_xlabel(r"Nelson diffusion $\nu$", color=ink)
    ax.set_ylabel(r"relaxation time $\tau$ (32$^2$ grain)", color=ink)
    ax.set_title(r"$\tau(\nu)$: pre-registered direction = decreasing",
                 fontsize=11, color=ink)
    ax.grid(True, which="major", color="#DDDDDD", lw=0.6)
    ax.tick_params(colors=muted)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    fig.suptitle("Nelson non-equilibrium relaxation in the tier3 2-D box: "
                 f"H-decay vs diffusion $\\nu$ (N = {N_PART:,}; "
                 "dotted: exponential fits; dX = (v+u)dt + $\\sqrt{2\\nu}$dW)",
                 fontsize=12, color=ink)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(OUTDIR, "m4_relax.png"), dpi=160)
    plt.close(fig)

    results["runtime_s"] = time.time() - t_start
    with open(os.path.join(OUTDIR, "m4_results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results["gates"], indent=2))
    print(f"TOTAL RUNTIME: {results['runtime_s']:.1f} s", flush=True)


if __name__ == "__main__":
    main()
