#!/usr/bin/env python3
"""
Tier-3 pilot: Born-rule relaxation (Valentini-Westman class) for the GUM
replication program.

2-D infinite square well [0, pi]^2, hbar = m = 1.
  phi_mn(x, y) = (2/pi) sin(m x) sin(n y),  E_mn = (m^2 + n^2)/2
  psi(x, y, t) = (1/sqrt(M)) sum_{(m,n) in S} exp(i theta_mn) phi_mn exp(-i E_mn t)
with fixed random phases theta_mn (numpy seed 42).

de Broglie guidance v = Im(grad psi / psi), grad psi analytic (exact time
evolution, no PDE solve). Initial ensemble rho_0 = |phi_11|^2 (non-equilibrium).
RK4 with per-particle substep guard (|v| dt <= 0.05 per substage) and a final
displacement cap; positions clipped strictly inside the box.

Diagnostics: coarse-grained H(t) = sum P ln(P/Q) at 32x32 and 16x16,
exponential fit of the initial decay, M=16 snapshots, and an equivariance
control sampled from |psi(0)|^2.

Fully seeded / reproducible. Outputs (written next to this script):
  hdata.csv, hbar_decay.png, snapshots_M16.png, results.json
"""

import json
import os
import time

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
OUTDIR = os.path.dirname(os.path.abspath(__file__))

L = np.pi                     # box side
DT = 2e-3                     # base integrator step
DT_OUT = 0.1                  # diagnostic sampling interval
T_FINAL = 4.0 * np.pi         # = exact revival period of psi (all E multiples of 1/2)
N_PART = 20000                # ensemble size
DISP_MAX = 0.05               # max allowed |v|*dt per substage
MAX_SUB = 32                  # cap on per-particle substeps
EPS_WALL = 1e-4               # keep particles strictly inside (eps, pi - eps)
NFINE = 256                   # fine grid for |psi|^2 cell integrals
CGS = (32, 16)                # coarse-graining levels
SEED_PHASES = 42              # fixed mode phases
SEED_SAMPLE = 12345           # ensemble sampling
KMAX_LIST = (2, 3, 4)         # -> M = 4, 9, 16

# ----------------------------------------------------------------------------
# Mode system
# ----------------------------------------------------------------------------
class ModeSystem:
    """Superposition of the first kmax x kmax box eigenmodes, equal weights,
    fixed random phases (seed 42)."""

    def __init__(self, kmax):
        self.kmax = kmax
        ks = np.arange(1, kmax + 1)
        mm, nn = np.meshgrid(ks, ks, indexing="ij")
        self.mm = mm.ravel().astype(float)
        self.nn = nn.ravel().astype(float)
        self.im = mm.ravel() - 1          # index into per-k trig tables
        self.jn = nn.ravel() - 1
        self.M = kmax * kmax
        self.E = 0.5 * (self.mm ** 2 + self.nn ** 2)
        rng = np.random.default_rng(SEED_PHASES)
        self.theta = rng.uniform(0.0, 2.0 * np.pi, self.M)
        self.ks = ks.astype(float)

    def coeffs(self, t):
        """Time-dependent complex coefficients (norm prefactor 2/pi omitted:
        it cancels in v = Im(grad psi / psi))."""
        return np.exp(1j * (self.theta - self.E * t)) / np.sqrt(self.M)

    def coeff_matrix(self, t):
        """Coefficients arranged as a (kmax, kmax) matrix C[m-1, n-1]."""
        return self.coeffs(t).reshape(self.kmax, self.kmax)

    def velocity(self, x, y, t):
        """de Broglie velocity v = Im(grad psi / psi) at particle positions.

        Uses the separable bilinear form psi_i = sx_i^T C sy_i with sx, sy the
        per-axis sine tables, so cost is O(N * kmax) rather than O(N * M).
        """
        C = self.coeff_matrix(t)
        Cr = np.ascontiguousarray(C.real)
        Ci = np.ascontiguousarray(C.imag)
        xk = x[:, None] * self.ks
        yk = y[:, None] * self.ks
        sx = np.sin(xk)                        # (N, kmax)
        cxm = np.cos(xk) * self.ks             # d/dx factors: m cos(m x)
        sy = np.sin(yk)
        cym = np.cos(yk) * self.ks
        ar, ai = sx @ Cr, sx @ Ci              # (N, kmax); reused for psi & d/dy
        pr = np.einsum("ij,ij->i", ar, sy)
        pi_ = np.einsum("ij,ij->i", ai, sy)
        gxr = np.einsum("ij,ij->i", cxm @ Cr, sy)
        gxi = np.einsum("ij,ij->i", cxm @ Ci, sy)
        gyr = np.einsum("ij,ij->i", ar, cym)
        gyi = np.einsum("ij,ij->i", ai, cym)
        dens = np.maximum(pr * pr + pi_ * pi_, 1e-30)
        vx = (gxi * pr - gxr * pi_) / dens     # Im(dpsi_x conj(psi)) / |psi|^2
        vy = (gyi * pr - gyr * pi_) / dens
        return vx, vy

    def psi2_grid(self, t, n=NFINE):
        """|psi(t)|^2 at the centers of an n x n fine grid (normalized density)."""
        xg = (np.arange(n) + 0.5) * L / n
        sxg = np.sin(np.outer(xg, self.ks))    # (n, kmax)
        C = self.coeff_matrix(t) * (2.0 / np.pi)
        psi = sxg @ C @ sxg.T                  # psi[i, j] at (x_i, y_j)
        return psi.real ** 2 + psi.imag ** 2

    def psi2_points(self, x, y, t):
        """|psi(t)|^2 at arbitrary points (for rejection sampling)."""
        C = self.coeff_matrix(t) * (2.0 / np.pi)
        sx = np.sin(x[:, None] * self.ks)
        sy = np.sin(y[:, None] * self.ks)
        pr = np.einsum("ij,ij->i", sx @ np.ascontiguousarray(C.real), sy)
        pi_ = np.einsum("ij,ij->i", sx @ np.ascontiguousarray(C.imag), sy)
        return pr ** 2 + pi_ ** 2


# ----------------------------------------------------------------------------
# Sampling (rejection, seeded)
# ----------------------------------------------------------------------------
def sample_phi11(n, rng):
    """rho_0 = |phi_11|^2 = (2/pi)^2 sin^2 x sin^2 y via rejection sampling."""
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


def sample_psi2(sys, n, rng):
    """Equilibrium control: sample from |psi(0)|^2 via rejection sampling."""
    bound = 1.2 * sys.psi2_grid(0.0, n=512).max()
    xs, ys, got = [], [], 0
    while got < n:
        prop = rng.uniform(0.0, L, size=(2 * n, 2))
        u = rng.uniform(0.0, bound, 2 * n)
        acc = u < sys.psi2_points(prop[:, 0], prop[:, 1], 0.0)
        xs.append(prop[acc, 0])
        ys.append(prop[acc, 1])
        got += int(acc.sum())
    x = np.concatenate(xs)[:n]
    y = np.concatenate(ys)[:n]
    return np.clip(x, EPS_WALL, L - EPS_WALL), np.clip(y, EPS_WALL, L - EPS_WALL)


# ----------------------------------------------------------------------------
# Integrator: RK4 + adaptive per-particle substepping + displacement cap
# ----------------------------------------------------------------------------
def _rk4_substep(sys, x, y, t, h, k1=None):
    k1x, k1y = sys.velocity(x, y, t) if k1 is None else k1
    k2x, k2y = sys.velocity(x + 0.5 * h * k1x, y + 0.5 * h * k1y, t + 0.5 * h)
    k3x, k3y = sys.velocity(x + 0.5 * h * k2x, y + 0.5 * h * k2y, t + 0.5 * h)
    k4x, k4y = sys.velocity(x + h * k3x, y + h * k3y, t + h)
    dx = h * (k1x + 2.0 * k2x + 2.0 * k3x + k4x) / 6.0
    dy = h * (k1y + 2.0 * k2y + 2.0 * k3y + k4y) / 6.0
    # hard guard: cap displacement magnitude (catches mid-substage node spikes)
    disp = np.hypot(dx, dy)
    big = disp > DISP_MAX
    if big.any():
        f = DISP_MAX / disp[big]
        dx[big] *= f
        dy[big] *= f
    x = np.clip(x + dx, EPS_WALL, L - EPS_WALL)
    y = np.clip(y + dy, EPS_WALL, L - EPS_WALL)
    return x, y


def macro_step(sys, x, y, t, dt):
    """One base step dt; particles whose initial |v| dt > DISP_MAX get their
    step subdivided (grouped by required substep count to stay vectorized)."""
    vx, vy = sys.velocity(x, y, t)
    speed = np.hypot(vx, vy)
    nsub = np.clip(np.ceil(speed * dt / DISP_MAX).astype(int), 1, MAX_SUB)
    uniq = np.unique(nsub)
    if uniq.size == 1 and uniq[0] == 1:
        # common case: no particle needs subdividing; reuse guard eval as k1
        return _rk4_substep(sys, x, y, t, dt, k1=(vx, vy))
    for n in uniq:
        idx = np.where(nsub == n)[0]
        xs, ys = x[idx], y[idx]
        h = dt / n
        k1 = (vx[idx], vy[idx])
        for j in range(n):
            xs, ys = _rk4_substep(sys, xs, ys, t + j * h, h, k1=k1)
            k1 = None
        x[idx], y[idx] = xs, ys
    return x, y


# ----------------------------------------------------------------------------
# Coarse-grained H-function
# ----------------------------------------------------------------------------
def h_function(x, y, q_fine, ncell):
    """H = sum_cells Pbar ln(Pbar/Qbar); Pbar from particle counts, Qbar from
    box-averaging |psi|^2 on the fine grid. Empty cells are skipped."""
    counts, _, _ = np.histogram2d(x, y, bins=ncell, range=[[0.0, L], [0.0, L]])
    p = counts / counts.sum()
    b = NFINE // ncell
    q = q_fine.reshape(ncell, b, ncell, b).mean(axis=(1, 3))
    q = q / q.sum()
    mask = p > 0
    qm = np.maximum(q[mask], 1e-300)
    return float(np.sum(p[mask] * np.log(p[mask] / qm)))


# ----------------------------------------------------------------------------
# Run driver
# ----------------------------------------------------------------------------
def run_ensemble(sys, x, y, label, keep_snapshots=False):
    """Evolve ensemble to T_FINAL, sampling H every DT_OUT (last sample lands
    exactly on t = 4*pi, the psi revival time)."""
    t0 = time.time()
    n_full = int(np.floor(T_FINAL / DT_OUT))          # full 0.1 output blocks
    steps_per_out = int(round(DT_OUT / DT))           # 50
    ts, h_by_cg = [], {cg: [] for cg in CGS}
    snaps = {}

    def record(t):
        qf = sys.psi2_grid(t)
        ts.append(t)
        for cg in CGS:
            h_by_cg[cg].append(h_function(x, y, qf, cg))

    if keep_snapshots:
        snaps[0.0] = (x.copy(), y.copy())
    t = 0.0
    record(t)
    for k in range(1, n_full + 1):
        for _ in range(steps_per_out):
            x, y = macro_step(sys, x, y, t, DT)
            t += DT
        record(t)
    # partial final block up to exactly t = 4*pi
    rem = T_FINAL - t
    if rem > 1e-12:
        n_extra = max(1, int(np.ceil(rem / DT)))
        h = rem / n_extra
        for _ in range(n_extra):
            x, y = macro_step(sys, x, y, t, h)
            t += h
        t = T_FINAL
        record(t)
    if keep_snapshots:
        snaps[T_FINAL] = (x.copy(), y.copy())
    print(f"[{label}] done in {time.time() - t0:.1f} s "
          f"(final H32 = {h_by_cg[32][-1]:.4f})", flush=True)
    out = {cg: np.array(v) for cg, v in h_by_cg.items()}
    return np.array(ts), out, snaps


def fit_tau(ts, hs):
    """Fit H ~ H0 exp(-t/tau) over the initial decay: up to where H first
    drops to 10% of H(0), or t = 2*pi, whichever comes first."""
    h0 = hs[0]
    below = np.where(hs <= 0.1 * h0)[0]
    t_end = ts[below[0]] if below.size else 2.0 * np.pi
    t_end = min(t_end, 2.0 * np.pi)
    sel = (ts <= t_end) & (hs > 0)
    if sel.sum() < 4:
        return np.nan, np.nan, t_end, np.nan
    slope, intercept = np.polyfit(ts[sel], np.log(hs[sel]), 1)
    tau = -1.0 / slope
    # goodness: R^2 of the linear fit in log space
    resid = np.log(hs[sel]) - (slope * ts[sel] + intercept)
    ss_tot = np.sum((np.log(hs[sel]) - np.log(hs[sel]).mean()) ** 2)
    r2 = 1.0 - np.sum(resid ** 2) / ss_tot if ss_tot > 0 else np.nan
    return tau, float(np.exp(intercept)), float(t_end), float(r2)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    t_start = time.time()
    results = {"N": N_PART, "dt": DT, "t_final": float(T_FINAL),
               "seed_phases": SEED_PHASES, "seed_sample": SEED_SAMPLE,
               "runs": {}}

    all_ts = None
    h_curves = {}          # (label, cg) -> array
    snaps_m16 = None
    systems = {}

    # --- non-equilibrium runs: M = 4, 9, 16 ---
    for kmax in KMAX_LIST:
        sys = ModeSystem(kmax)
        systems[kmax] = sys
        m = sys.M
        rng = np.random.default_rng(SEED_SAMPLE + m)
        x, y = sample_phi11(N_PART, rng)
        ts, hcg, snaps = run_ensemble(sys, x, y, f"M={m}",
                                      keep_snapshots=(m == 16))
        all_ts = ts
        if m == 16:
            snaps_m16 = snaps
        rr = {}
        for cg in CGS:
            h_curves[(f"M{m}", cg)] = hcg[cg]
            tau, h0fit, t_end, r2 = fit_tau(ts, hcg[cg])
            rr[f"cg{cg}"] = {"H_initial": float(hcg[cg][0]),
                             "H_final": float(hcg[cg][-1]),
                             "tau": float(tau), "H0_fit": h0fit,
                             "fit_window_t_end": t_end, "fit_r2": r2}
        results["runs"][f"M={m}"] = rr

    # --- equivariance control: ensemble drawn from |psi(0)|^2, M = 16 ---
    sys16 = systems[4]
    rng = np.random.default_rng(SEED_SAMPLE + 999)
    xc, yc = sample_psi2(sys16, N_PART, rng)
    ts, hcg, _ = run_ensemble(sys16, xc, yc, "control(M=16,equilibrium)")
    rr = {}
    for cg in CGS:
        h_curves[("ctrl", cg)] = hcg[cg]
        rr[f"cg{cg}"] = {"H_initial": float(hcg[cg][0]),
                         "H_final": float(hcg[cg][-1]),
                         "H_max": float(hcg[cg].max()),
                         "H_mean": float(hcg[cg].mean())}
    results["runs"]["control"] = rr

    # noise floors ~ n_occupied_cells / (2N)
    results["noise_floor"] = {f"cg{cg}": cg * cg / (2.0 * N_PART) for cg in CGS}

    # ------------------------------------------------------------------
    # hdata.csv
    # ------------------------------------------------------------------
    cols = ["t"]
    data = [all_ts]
    for m in (4, 9, 16):
        for cg in CGS:
            cols.append(f"H{cg}_M{m}")
            data.append(h_curves[(f"M{m}", cg)])
    for cg in CGS:
        cols.append(f"H{cg}_ctrl")
        data.append(h_curves[("ctrl", cg)])
    arr = np.column_stack(data)
    np.savetxt(os.path.join(OUTDIR, "hdata.csv"), arr,
               delimiter=",", header=",".join(cols), comments="",
               fmt="%.8g")

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    # categorical colors (fixed order, CVD-safe), control in neutral gray
    # validated categorical trio (CVD-safe, fixed order); control is a neutral
    # reference line distinguished by dash style, not a categorical slot
    col = {4: "#2a78d6", 9: "#eda100", 16: "#e34948"}
    gray = "#898781"
    ink, muted = "#0b0b0b", "#52514e"

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    for ax, cg in zip(axes, CGS):
        floor = results["noise_floor"][f"cg{cg}"]
        for m in (4, 9, 16):
            ax.plot(all_ts, h_curves[(f"M{m}", cg)], color=col[m],
                    lw=1.8, label=f"M = {m}")
            tau = results["runs"][f"M={m}"][f"cg{cg}"]["tau"]
            h0f = results["runs"][f"M={m}"][f"cg{cg}"]["H0_fit"]
            tend = results["runs"][f"M={m}"][f"cg{cg}"]["fit_window_t_end"]
            tt = np.linspace(0, tend, 50)
            ax.plot(tt, h0f * np.exp(-tt / tau), color=col[m], lw=1.0,
                    ls=":", alpha=0.9)
        ax.plot(all_ts, np.maximum(h_curves[("ctrl", cg)], 1e-6), color=gray,
                lw=1.4, ls="--", label="control ($|\\psi(0)|^2$, M = 16)")
        ax.axhline(floor, color=muted, lw=1.0, ls="-.",
                   label=f"noise floor $\\approx$ {floor:.3g}")
        ax.axvline(T_FINAL, color=muted, lw=1.0, ls=":")
        ax.text(T_FINAL, ax.get_ylim()[0], " ", fontsize=8)
        ax.set_yscale("log")
        ax.set_xlabel("t", color=ink)
        ax.set_title(f"{cg}$\\times${cg} coarse-graining", fontsize=11,
                     color=ink)
        ax.grid(True, which="major", color="#DDDDDD", lw=0.6)
        ax.tick_params(colors=muted)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel(r"$\bar{H}(t)=\sum \bar{P}\,\ln(\bar{P}/\bar{Q})$",
                       color=ink)
    # revival marker annotation on left panel
    for ax in axes:
        ymin, ymax = ax.get_ylim()
        ax.annotate("$\\psi$ revival\n$t = 4\\pi$", xy=(T_FINAL, ymax),
                    xytext=(T_FINAL - 0.15, ymax * 0.55), ha="right",
                    fontsize=8, color=muted)
    axes[0].legend(frameon=False, fontsize=9, loc="lower left")
    fig.suptitle("Coarse-grained H-function decay, de Broglie ensemble in a 2-D box "
                 f"(N = {N_PART:,}; dotted: exponential fits)",
                 fontsize=12, color=ink)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(OUTDIR, "hbar_decay.png"), dpi=160)
    plt.close(fig)

    # snapshots for M = 16
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 9.2))
    times = [0.0, T_FINAL]
    labels = ["t = 0", "t = 4$\\pi$ ($\\psi$ revival)"]
    for j, (tt, lab) in enumerate(zip(times, labels)):
        xs, ys = snaps_m16[tt]
        ax = axes[0, j]
        sub = slice(0, 8000)  # subsample for legibility
        ax.scatter(xs[sub], ys[sub], s=1.2, c="#2a78d6", alpha=0.35,
                   linewidths=0)
        ax.set_title(f"particles, {lab}", fontsize=10, color=ink)
        q = systems[4].psi2_grid(tt)
        ax2 = axes[1, j]
        ax2.imshow(q.T, origin="lower", extent=(0, L, 0, L), cmap="viridis",
                   aspect="equal")
        ax2.set_title(f"$|\\psi|^2$, {lab}", fontsize=10, color=ink)
        for a in (ax, ax2):
            a.set_xlim(0, L)
            a.set_ylim(0, L)
            a.set_aspect("equal")
            a.set_xlabel("x", color=muted)
            a.set_ylabel("y", color=muted)
            a.tick_params(colors=muted)
    fig.suptitle("M = 16: particle ensemble vs quantum density "
                 "(ensemble starts as $|\\phi_{11}|^2$, relaxes toward $|\\psi|^2$)",
                 fontsize=11, color=ink)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(OUTDIR, "snapshots_M16.png"), dpi=160)
    plt.close(fig)

    results["runtime_s"] = time.time() - t_start
    with open(os.path.join(OUTDIR, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(json.dumps(results, indent=2))
    print(f"TOTAL RUNTIME: {results['runtime_s']:.1f} s", flush=True)


if __name__ == "__main__":
    main()
