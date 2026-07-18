#!/usr/bin/env python3
"""
WORKSTREAM O1 -- the flow-map first-crossing proof (T4-W5 terminal shape,
per F-T7-N3).  File id: F-T7-O1.

Theorem (backward-reachability first-crossing criterion; the proof N3
showed is required).  Let v(x, z, t) be the guidance velocity field of the
L5prime/L7b bench member n_y,
    v_x = Im(f'/f)(x,t) - n_y pi cot(pi z),   v_z = n_y Re(f'/f)(x,t),
defined on U = {rho > 0}.  For t > 0 the free evolution of the compactly
supported initial state is real-analytic in x (Paley-Wiener: f(.,t) is the
Fourier transform of a compactly-supported-in-x, L1 datum times a Gaussian
phase), so v is C^infty (locally Lipschitz) on U; by Picard-Lindelof the
integral curves of dX/ds = v(X, s) through any (p, t) in U are UNIQUE both
FORWARD and BACKWARD.  Let X(s) = Phi_{s,t}(p) be the backward segment,
s in [0, t], from a detector point p = (d, z).
  (i) If max_{s<t} X_x(s) >= d, then by continuity and the intermediate
      value theorem the unique trajectory through p already reached the
      line x = d at some s0 < t; its first crossing time is <= s0 < t, so
      it contributes NO first-crossing (arrival) density at time t.
  (ii) The first-crossing density at (z, t) is dominated by the upward
      crossing flux: dP_first/(dz dt) <= rho(d,z,t) max(v_x(d,z,t), 0)
      (a first crossing at (z,t) is in particular an upward crossing, and
      the equivariant transported density is rho).
  (iii) Hence if EVERY detector point of a set B in (line x kill bins) is
      of type (i) except a subset N, the first-crossing measure of B is
      <= integral over N of rho (v_x)_+ dz dt.  If that bound is 0 (or
      negligible), the arrival density in B vanishes at field level.  QED
(The forward statement used in (i): a trajectory with X_x(s0) >= d and
X_x(0) < d has a first crossing <= s0.  If X_x(0) >= d the trajectory
STARTS at/beyond the line and again contributes no first-crossing at t.)

This workstream EXECUTES the criterion on the validated engine: from every
kill-bin detector point (d = 1, z, t), z on a 200-point grid, t on 40
points per bin per n_y in {+1, -1, +0.75, -0.75} (+0.25 discrimination
control), integrate the SAME field's guidance ODE BACKWARD (RK4, fields
from the exact propagator) from t down to 0 and classify:
  (a) PRE-CROSSED : max_{s<t} X_x >= d + 1e-3 before any rho-floor event
      (crossing epoch + margin recorded);
  (b) RHO-FLOOR   : the backward path enters rho < 1e-10 x rho_peak(s)
      (near-node / near-wall zones; the field is ill-conditioned there) --
      counted, and the first-crossing measure they can carry bounded by
      the detector-cell flux rho (v_x)_+ dz dt;
  (c) VIOLATION   : the path reaches t = 0 inside supp(phi_0) with
      max X_x < d + 1e-9 throughout -- a genuine fresh kill-bin arrival,
      which would REFUTE the L7b zeros;
plus GRAZE-UNCERTAIN (1e-9 < margin < 1e-3: not certified either way;
flux added to the (b) bound), FRESH-OUT (reaches t = 0 outside
supp(phi_0): rho_0 = 0 there, zero measure; flux bounded like (b)) and
BOX-EXIT (leaves [-15,45] rim; flux bounded like (b)).

Engine inherited verbatim from L7b/N3: factorized exact-in-z state
phi = f(x,t) sin(pi z) e^{-i pi^2 t/2}, f propagated exactly in 1-D
Fourier space, k0 = 2, C2 smoothstep taper, literal amplitude, box
x in [-15, 45], Nx = 2048 reference, dt = 2.5e-4 <= 5e-4, clamp +-500,
density floor 1e-14 relative, z clipped to [1e-9, 1-1e-9]; velocity
decomposition v = a + n_y b exactly as in l7b_run.py.  Backward fields by
the same half-step accumulation with the conjugate phase (dt-halving
probes the accumulation).  Refinements: Nx = 4096 and dt = 1.25e-4, both
over the FULL grid (not just boundary-adjacent).  hbar = m = 1.

Gates (pre-registered, ROADMAP_v9 Phase O / O1):
  G1 engine bars 1e-6 class or better (2-D wall diagnostic + 1-D norm/
     energy + backward reversibility F(0) == F0).
  G2 zero VIOLATIONS at reference AND at both refinements; classification
     stable.
  G3 min over certified class-(a) of (max X_x - d) >> per-path refinement
     difference (ratio >= 10 vs p95).
  G4 flux bound over {RHO-FLOOR, GRAZE, FRESH-OUT, BOX-EXIT} <= 1e-3 of
     the |phi_0|^2 measure per kill n_y (else PARTIAL with the bound).
  G5 discrimination: n_y = +0.25 backward runs from the same (populated)
     bins show a substantial legitimately-fresh fraction: >= 50 fresh
     points and fresh-cell flux >= 0.3 x the L7b measured bin probability
     (report the full comparison).
"""
import json
import os
import time

import numpy as np
from scipy import fft as sfft

# ----------------------------------------------------------------------
# Parameters (inherited verbatim from L7b / N3)
# ----------------------------------------------------------------------
L = 1.0
SIGMA_X = 1.0
X0 = -5.0
K0 = 2.0
W_IN, W_OUT = 2.5, 3.5
D_NEAR = 1.0
VCLAMP = 500.0
RHO_EPS_REL = 1e-14
XMIN, XMAX = -15.0, 45.0
NX_REF = 2048
DT_REF = 2.5e-4

KILL_BINS = [(6.4, 6.8), (6.8, 7.2)]
NY_KILL = [1.0, -1.0, 0.75, -0.75]
NY_CTRL = 0.25
NY_SET = NY_KILL + [NY_CTRL]

NZ_DET = 200                     # z-points on the line (cell centres)
NT_PER_BIN = 40                  # t-points per bin (cell centres)

RHO_FLOOR_REL = 1e-10            # class-(b) floor: rho < 1e-10 x peak(s)
MARGIN_EPS = 1e-9                # "touched the line" threshold
MARGIN_CLASS = 1e-3              # certified pre-crossing threshold
MARGIN_CAP = 2.0                 # freeze once margin >= cap (class (a) settled)

ZCLIP = 1e-9                     # engine z-clip

OUTDIR = os.path.dirname(os.path.abspath(__file__))
T6 = os.path.join(os.path.dirname(os.path.dirname(OUTDIR)), "tier6-foundations")
L7B_JSON = os.path.join(T6, "L7b", "l7b_results.json")

CLS_NAME = {0: "A_PRECROSSED", 1: "B_RHOFLOOR", 2: "C_FRESH_IN_SUPP",
            3: "FRESH_OUT_SUPP", 4: "BOX_EXIT", 5: "GRAZE_UNCERTAIN"}
BOUNDED_CLS = (1, 3, 4, 5)       # classes whose flux enters the G4 bound


# ----------------------------------------------------------------------
# Initial profile + engine pieces (verbatim from l7b_run.py)
# ----------------------------------------------------------------------
def smoothstep(t):
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))


def taper(x):
    u = np.abs(np.asarray(x, float) - X0) / SIGMA_X
    W = np.ones_like(u)
    W[u >= W_OUT] = 0.0
    m = (u > W_IN) & (u < W_OUT)
    W[m] = smoothstep((W_OUT - u[m]) / (W_OUT - W_IN))
    return W


def amp(x):
    d2 = (np.asarray(x, float) - X0) ** 2 / SIGMA_X ** 2
    return np.exp(-d2 / 2.0) * taper(x)


def build_f(Nx, xmin, xmax):
    dx = (xmax - xmin) / Nx
    xg = xmin + dx * np.arange(Nx)
    f = amp(xg).astype(complex)
    f *= np.exp(1j * K0 * xg)
    f /= np.sqrt(np.sum(np.abs(f) ** 2) * dx)
    kx = 2.0 * np.pi * sfft.fftfreq(Nx, dx)
    return xg, dx, f, kx


def lin_interp(V, px, xg0, dx, Nx):
    g = np.clip((px - xg0) / dx, 0.0, Nx - 1.000001)
    i = g.astype(np.int64)
    fr = g - i
    return V[i] * (1 - fr) + V[i + 1] * fr


# ----------------------------------------------------------------------
# 2-D wall-amplitude diagnostic (verbatim from l7b_run.py / n3_field.py)
# ----------------------------------------------------------------------
def diag_2d(Nx=2048, Nz=128, xmin=XMIN, xmax=XMAX, t_final=16.0, n_checks=17):
    dx = (xmax - xmin) / Nx
    xg = xmin + dx * np.arange(Nx)
    dz = L / Nz
    zg_ext = dz * np.arange(2 * Nz)
    fx = amp(xg).astype(complex) * np.exp(1j * K0 * xg)
    gz = np.sin(np.pi * zg_ext / L)
    phi0 = fx[:, None] * gz[None, :]
    dA = dx * dz
    phi0 /= np.sqrt(np.sum(np.abs(phi0[:, :Nz]) ** 2) * dA)
    kx = 2.0 * np.pi * sfft.fftfreq(Nx, dx)
    kz = 2.0 * np.pi * sfft.fftfreq(2 * Nz, dz)
    K2 = kx[:, None] ** 2 + kz[None, :] ** 2
    F0 = sfft.fft2(phi0)
    norm0 = float(np.sum(np.abs(phi0[:, :Nz]) ** 2) * dA)
    Pspec = np.sum(np.abs(F0) ** 2)
    E0 = float(np.sum(0.5 * K2 * np.abs(F0) ** 2) / Pspec)
    rowpow = np.sum(np.abs(F0) ** 2, axis=0)
    mode1 = float(rowpow[1] + rowpow[2 * Nz - 1])
    leak = float((rowpow.sum() - mode1) / rowpow.sum())
    norms, walls, energies = [], [], []
    for tchk in np.linspace(0.0, t_final, n_checks):
        Ft = F0 * np.exp(-0.5j * K2 * tchk)
        phit = sfft.ifft2(Ft)
        norms.append(float(np.sum(np.abs(phit[:, :Nz]) ** 2) * dA))
        walls.append(float(max(np.abs(phit[:, 0]).max(),
                               np.abs(phit[:, Nz]).max())))
        energies.append(float(np.sum(0.5 * K2 * np.abs(Ft) ** 2)
                              / np.sum(np.abs(Ft) ** 2)))
    return dict(norm_deviation_max=float(np.max(np.abs(np.array(norms)
                                                       - norm0))),
                wall_amplitude_max=float(np.max(walls)),
                energy_drift_rel=float(np.max(np.abs(np.array(energies)
                                                     - E0)) / abs(E0)),
                energy_mean=E0, z_mode_leakage=leak)


# ----------------------------------------------------------------------
# Detector grid (shared by all runs)
# ----------------------------------------------------------------------
def detector_grid():
    zc = (np.arange(NZ_DET) + 0.5) / NZ_DET
    tcs = []
    for (b0, b1) in KILL_BINS:
        w = (b1 - b0) / NT_PER_BIN
        tcs.append(b0 + (np.arange(NT_PER_BIN) + 0.5) * w)
    tc = np.concatenate(tcs)                       # 80 start times
    NT_TOT = tc.size
    # layout: index = (i_ny * NT_TOT + j_t) * NZ_DET + i_z
    ny_flat = np.repeat(np.array(NY_SET), NT_TOT * NZ_DET)
    t_flat = np.tile(np.repeat(tc, NZ_DET), len(NY_SET))
    z_flat = np.tile(zc, NT_TOT * len(NY_SET))
    return zc, tc, ny_flat, t_flat, z_flat


ZC, TC, NY_FLAT, T_FLAT, Z_FLAT = detector_grid()
NTOT = NY_FLAT.size
S_START = float(TC.max())                          # 7.195
DZ_CELL = 1.0 / NZ_DET
DT_CELL = (KILL_BINS[0][1] - KILL_BINS[0][0]) / NT_PER_BIN


# ----------------------------------------------------------------------
# Backward run: activate at start times, RK4 backward to s = 0
# ----------------------------------------------------------------------
def backward_run(Nx, dt, label, hist_idx=None, hist_every=25):
    xg, dx, f0, kx = build_f(Nx, XMIN, XMAX)
    F0 = sfft.fft(f0)
    nsteps = int(round(S_START / dt))
    assert abs(nsteps * dt - S_START) < 1e-12
    act_step = np.round((S_START - T_FLAT) / dt).astype(np.int64)
    assert np.max(np.abs((S_START - T_FLAT) / dt - act_step)) < 1e-9
    Ehb = np.exp(+0.5j * kx ** 2 * (dt / 2.0))     # backward half step
    F = F0 * np.exp(-0.5j * kx ** 2 * S_START)     # exact field at S_START
    norm0 = float(np.sum(np.abs(f0) ** 2) * dx)
    p0 = np.abs(F0) ** 2
    E0 = float(np.sum(0.5 * kx ** 2 * p0) / np.sum(p0))

    def fields(Fc):
        f = sfft.ifft(Fc)
        fp = sfft.ifft(1j * kx * Fc)
        rho = np.abs(f) ** 2
        den = np.maximum(rho, RHO_EPS_REL * rho.max())
        w = np.conj(f) * fp
        a = np.clip(np.imag(w) / den, -VCLAMP, VCLAMP)
        b = np.clip(np.real(w) / den, -VCLAMP, VCLAMP)
        return a, b, rho, rho.max()

    def vel(ab, px, pz, ny):
        a, b = ab[0], ab[1]
        va = lin_interp(a, px, xg[0], dx, Nx)
        vb = lin_interp(b, px, xg[0], dx, Nx)
        vx = np.clip(va - ny * np.pi / np.tan(np.pi * pz / L),
                     -VCLAMP, VCLAMP)
        return vx, ny * vb

    X = np.full(NTOT, np.nan)
    Z = np.full(NTOT, np.nan)
    started = np.zeros(NTOT, bool)
    frozen = np.zeros(NTOT, bool)
    maxX = np.full(NTOT, -np.inf)
    pc_cert = np.zeros(NTOT, bool)                 # margin >= MARGIN_CLASS
    s_pc = np.full(NTOT, np.nan)                   # certified crossing epoch
    floored = np.zeros(NTOT, bool)
    s_floor = np.full(NTOT, np.nan)
    boxexit = np.zeros(NTOT, bool)
    vx_det = np.full(NTOT, np.nan)
    rho_det = np.full(NTOT, np.nan)

    if hist_idx is not None:
        n_h = nsteps // hist_every + 2
        histX = np.full((n_h, hist_idx.size), np.nan)
        histS = np.full(n_h, np.nan)
        jh = 0

    check_every = max(1, nsteps // 16)
    norm_dev = 0.0
    e_drift = 0.0
    ab = fields(F)
    t_start_wall = time.time()
    for n in range(nsteps):
        s = S_START - n * dt
        if n % check_every == 0:
            f = sfft.ifft(F)
            norm_dev = max(norm_dev, abs(float(np.sum(np.abs(f) ** 2) * dx)
                                         - norm0))
            p = np.abs(F) ** 2
            e_drift = max(e_drift, abs(float(np.sum(0.5 * kx ** 2 * p)
                                             / np.sum(p)) - E0) / abs(E0))
        # activate trajectories starting at this clock time
        act = act_step == n
        if act.any():
            started[act] = True
            X[act] = D_NEAR
            Z[act] = Z_FLAT[act]
            vx0, _ = vel(ab, X[act], Z[act], NY_FLAT[act])
            vx_det[act] = vx0
            rho_det[act] = (lin_interp(ab[2], X[act], xg[0], dx, Nx)
                            * 2.0 * np.sin(np.pi * Z_FLAT[act]) ** 2)
        if hist_idx is not None and n % hist_every == 0:
            histX[jh] = X[hist_idx]
            histS[jh] = s
            jh += 1
        # backward fields at s - dt/2 and s - dt
        F = F * Ehb
        ab_h = fields(F)
        F = F * Ehb
        ab_n = fields(F)
        idx = np.where(started & ~frozen)[0]
        if idx.size:
            px, pz = X[idx], Z[idx]
            nyv = NY_FLAT[idx]
            k1x, k1z = vel(ab, px, pz, nyv)
            k2x, k2z = vel(ab_h, px - 0.5 * dt * k1x,
                           np.clip(pz - 0.5 * dt * k1z, ZCLIP, L - ZCLIP), nyv)
            k3x, k3z = vel(ab_h, px - 0.5 * dt * k2x,
                           np.clip(pz - 0.5 * dt * k2z, ZCLIP, L - ZCLIP), nyv)
            k4x, k4z = vel(ab_n, px - dt * k3x,
                           np.clip(pz - dt * k3z, ZCLIP, L - ZCLIP), nyv)
            nx_ = px - dt / 6.0 * (k1x + 2 * k2x + 2 * k3x + k4x)
            nz_ = np.clip(pz - dt / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z),
                          ZCLIP, L - ZCLIP)
            X[idx] = nx_
            Z[idx] = nz_
            maxX[idx] = np.maximum(maxX[idx], nx_)
            # certified pre-crossing (checked BEFORE the floor: within-step
            # ties resolve to (a); one-step ambiguity is refinement-probed)
            newpc = (~pc_cert[idx]) & (~floored[idx]) \
                & (maxX[idx] >= D_NEAR + MARGIN_CLASS)
            if newpc.any():
                pc_cert[idx[newpc]] = True
                s_pc[idx[newpc]] = s - dt
            # rho-floor at the new position/time
            rho_rel = (lin_interp(ab_n[2], nx_, xg[0], dx, Nx)
                       * np.sin(np.pi * nz_) ** 2) / ab_n[3]
            flr = rho_rel < RHO_FLOOR_REL
            if flr.any():
                newf = idx[flr]
                floored[newf] = True
                s_floor[newf] = s - dt
                frozen[newf] = True
            # freeze once the (a) classification is settled with margin
            capd = maxX[idx] - D_NEAR >= MARGIN_CAP
            if capd.any():
                frozen[idx[capd]] = True
            # box rim (engine freeze rule)
            esc = (nx_ < XMIN + 0.5) | (nx_ > XMAX - 0.5)
            if esc.any():
                ei = idx[esc]
                boxexit[ei] = True
                frozen[ei] = True
        ab = ab_n
        if n % 4000 == 0:
            print("  [%s] step %d/%d s=%.3f active=%d elapsed=%.0fs"
                  % (label, n, nsteps, s, int((started & ~frozen).sum()),
                     time.time() - t_start_wall), flush=True)
    if hist_idx is not None:
        histX[jh] = X[hist_idx]
        histS[jh] = 0.0
        jh += 1
    # backward reversibility: F should have returned to F0 exactly
    rev = float(np.max(np.abs(F - F0)) / np.max(np.abs(F0)))

    # classification
    margin = maxX - D_NEAR
    reached0 = started & ~floored & ~boxexit & ~pc_cert & ~frozen
    graze = (margin > MARGIN_EPS) & (margin < MARGIN_CLASS)
    insupp = np.abs(X - X0) <= W_OUT * SIGMA_X
    cls = np.full(NTOT, -1, np.int64)
    cls[pc_cert] = 0
    cls[~pc_cert & floored] = 1
    cls[~pc_cert & ~floored & boxexit] = 4
    cls[reached0 & graze] = 5
    cls[reached0 & ~graze & insupp] = 2
    cls[reached0 & ~graze & ~insupp] = 3
    assert not (cls == -1).any()

    flux = rho_det * np.maximum(vx_det, 0.0) * DZ_CELL * DT_CELL
    out = dict(label=label, Nx=Nx, dt=dt, nsteps=nsteps,
               elapsed=time.time() - t_start_wall,
               engine=dict(norm_dev=norm_dev, energy_drift_rel=e_drift,
                           reversibility_rel=rev, E0=E0))
    if hist_idx is not None:
        out["hist"] = (histS[:jh], histX[:jh])
    return out, dict(cls=cls, margin=margin, s_pc=s_pc, s_floor=s_floor,
                     Xend=X.copy(), Zend=Z.copy(), vx_det=vx_det,
                     rho_det=rho_det, flux=flux)


# ----------------------------------------------------------------------
def per_ny_tables(rec):
    """Per-n_y, per-bin classification counts and flux accounting."""
    cls, flux = rec["cls"], rec["flux"]
    tables = {}
    for ny in NY_SET:
        sel_ny = NY_FLAT == ny
        entry = {}
        for (b0, b1) in KILL_BINS:
            m = sel_ny & (T_FLAT > b0) & (T_FLAT < b1)
            e = {}
            for c, name in CLS_NAME.items():
                e[name] = int(np.sum(cls[m] == c))
            e["flux_total_pos"] = float(flux[m].sum())
            e["flux_bounded"] = float(flux[m & np.isin(cls, BOUNDED_CLS)].sum())
            e["flux_violation"] = float(flux[m & (cls == 2)].sum())
            e["flux_fresh"] = e["flux_violation"]
            entry["bin_%g_%g" % (b0, b1)] = e
        m = sel_ny
        entry["all_bins"] = {CLS_NAME[c]: int(np.sum(cls[m] == c))
                             for c in CLS_NAME}
        entry["all_bins"]["flux_bounded"] = float(
            flux[m & np.isin(cls, BOUNDED_CLS)].sum())
        entry["all_bins"]["flux_violation"] = float(
            flux[m & (cls == 2)].sum())
        entry["all_bins"]["n_violation"] = int(np.sum((cls == 2) & m)) \
            if ny != NY_CTRL else 0
        entry["all_bins"]["n_fresh"] = int(np.sum((cls == 2) & m))
        # class-(a) margin statistics
        ma = rec["margin"][m & (cls == 0)]
        if ma.size:
            entry["class_a_margin"] = dict(
                n=int(ma.size), min=float(ma.min()),
                median=float(np.median(ma)), max=float(ma.max()),
                frac_capped=float(np.mean(ma >= MARGIN_CAP)))
            spc = rec["s_pc"][m & (cls == 0)]
            entry["class_a_crossing_epoch"] = dict(
                min=float(np.nanmin(spc)), median=float(np.nanmedian(spc)),
                max=float(np.nanmax(spc)))
        tables["%+.2f" % ny] = entry
    return tables


def main():
    t0_wall = time.time()
    print("O1 backward-reachability run: NTOT=%d start points, S_START=%.3f"
          % (NTOT, S_START), flush=True)

    # history subset for the figure: 10 z x 2 t per n_y
    NT_TOT = TC.size
    hist_idx = []
    for i_ny in range(len(NY_SET)):
        for j_t in (20, 60):
            for i_z in range(10, NZ_DET, 20):
                hist_idx.append((i_ny * NT_TOT + j_t) * NZ_DET + i_z)
    hist_idx = np.array(sorted(set(hist_idx)))

    print("2-D wall diagnostic ...", flush=True)
    d2 = diag_2d()
    print("  ", json.dumps(d2), flush=True)

    runs = {}
    recs = {}
    for label, Nx, dt in [("ref", NX_REF, DT_REF),
                          ("nx2", 2 * NX_REF, DT_REF),
                          ("dt2", NX_REF, DT_REF / 2.0)]:
        print("run %s: Nx=%d dt=%g ..." % (label, Nx, dt), flush=True)
        runs[label], recs[label] = backward_run(
            Nx, dt, label, hist_idx=hist_idx if label == "ref" else None)
        r = runs[label]
        print("  elapsed %.1f s  norm_dev %.2e  E_drift %.2e  rev %.2e"
              % (r["elapsed"], r["engine"]["norm_dev"],
                 r["engine"]["energy_drift_rel"],
                 r["engine"]["reversibility_rel"]), flush=True)
        cls = recs[label]["cls"]
        for c in CLS_NAME:
            print("    %-16s %d" % (CLS_NAME[c], int(np.sum(cls == c))),
                  flush=True)

    ref = recs["ref"]
    tab_ref = per_ny_tables(ref)
    tab_nx2 = per_ny_tables(recs["nx2"])
    tab_dt2 = per_ny_tables(recs["dt2"])

    # ------------------------------------------------------------------
    # classification stability + per-path refinement differences
    # ------------------------------------------------------------------
    stab = {}
    for lab, rec in (("nx2", recs["nx2"]), ("dt2", recs["dt2"])):
        agree = ref["cls"] == rec["cls"]
        moved = np.where(~agree)[0]
        pairs = {}
        for i in moved:
            key = "%s->%s" % (CLS_NAME[ref["cls"][i]], CLS_NAME[rec["cls"][i]])
            pairs[key] = pairs.get(key, 0) + 1
        # margin differences over stable class-(a), uncapped in both
        both_a = (ref["cls"] == 0) & (rec["cls"] == 0)
        uncap = both_a & (ref["margin"] < MARGIN_CAP) \
            & (rec["margin"] < MARGIN_CAP)
        dmar = np.abs(ref["margin"][uncap] - rec["margin"][uncap])
        stab[lab] = dict(
            n_agree=int(agree.sum()), n_moved=int(moved.size),
            agree_frac=float(agree.mean()), moved_pairs=pairs,
            moved_flux_ref=float(ref["flux"][moved].sum()),
            margin_diff_uncapped=dict(
                n=int(dmar.size), median=float(np.median(dmar)),
                p95=float(np.percentile(dmar, 95)), max=float(dmar.max()))
            if dmar.size else None)

    # ------------------------------------------------------------------
    # G1: engine bars
    # ------------------------------------------------------------------
    g1 = dict(diag_2d=d2,
              runs={lab: runs[lab]["engine"] for lab in runs})
    g1_ok = (d2["norm_deviation_max"] < 1e-6
             and d2["wall_amplitude_max"] < 1e-10
             and d2["energy_drift_rel"] <= 1e-6
             and all(runs[lab]["engine"]["norm_dev"] < 1e-6
                     and runs[lab]["engine"]["energy_drift_rel"] <= 1e-6
                     and runs[lab]["engine"]["reversibility_rel"] < 1e-9
                     for lab in runs))
    g1["verdict"] = "PASS" if g1_ok else "FAIL"

    # ------------------------------------------------------------------
    # G2: zero violations, all kill n_y, all three resolutions
    # ------------------------------------------------------------------
    nv = {}
    for lab, tab in (("ref", tab_ref), ("nx2", tab_nx2), ("dt2", tab_dt2)):
        nv[lab] = {("%+.2f" % ny): tab["%+.2f" % ny]["all_bins"]["n_fresh"]
                   for ny in NY_KILL}
    zero_all = all(v == 0 for lab in nv for v in nv[lab].values())
    g2 = dict(n_violations=nv, zero_at_all_resolutions=zero_all,
              classification_stability=stab)
    g2["verdict"] = "PASS" if zero_all else "FAIL"

    # ------------------------------------------------------------------
    # G3: certified margin vs per-path discretization error
    # ------------------------------------------------------------------
    kill_sel = np.isin(NY_FLAT, NY_KILL)
    a_sel = kill_sel & (ref["cls"] == 0)
    min_margin = float(ref["margin"][a_sel].min())
    err95 = max(stab["nx2"]["margin_diff_uncapped"]["p95"]
                if stab["nx2"]["margin_diff_uncapped"] else 0.0,
                stab["dt2"]["margin_diff_uncapped"]["p95"]
                if stab["dt2"]["margin_diff_uncapped"] else 0.0)
    errmed = max(stab["nx2"]["margin_diff_uncapped"]["median"]
                 if stab["nx2"]["margin_diff_uncapped"] else 0.0,
                 stab["dt2"]["margin_diff_uncapped"]["median"]
                 if stab["dt2"]["margin_diff_uncapped"] else 0.0)
    ratio = min_margin / err95 if err95 > 0 else np.inf
    g3 = dict(min_margin_class_a_kill=min_margin,
              margin_err_p95=err95, margin_err_median=errmed,
              ratio_min_margin_over_p95err=float(ratio),
              per_ny_margin={("%+.2f" % ny):
                             tab_ref["%+.2f" % ny].get("class_a_margin")
                             for ny in NY_KILL},
              crossing_epochs={("%+.2f" % ny):
                               tab_ref["%+.2f" % ny].get(
                                   "class_a_crossing_epoch")
                               for ny in NY_KILL})
    g3["verdict"] = "PASS" if ratio >= 10.0 else "FAIL"

    # ------------------------------------------------------------------
    # G4: bounded-class flux accounting
    # ------------------------------------------------------------------
    fb = {("%+.2f" % ny): tab_ref["%+.2f" % ny]["all_bins"]["flux_bounded"]
          for ny in NY_KILL}
    fbmax = max(fb.values())
    g4 = dict(flux_bounded_per_ny=fb, flux_bounded_max=fbmax,
              spec="total |phi0|^2-measure reachable through class-(b) "
                   "(+ graze/out-supp/box-exit) <= 1e-3 per kill n_y",
              per_bin={("%+.2f" % ny): {
                  k: dict(flux_bounded=tab_ref["%+.2f" % ny][k]["flux_bounded"],
                          flux_total_pos=tab_ref["%+.2f" % ny][k][
                              "flux_total_pos"])
                  for k in tab_ref["%+.2f" % ny] if k.startswith("bin_")}
                  for ny in NY_KILL})
    g4["verdict"] = "PASS" if fbmax <= 1e-3 else "PARTIAL"

    # ------------------------------------------------------------------
    # G5: discrimination control at n_y = +0.25
    # ------------------------------------------------------------------
    l7b = json.load(open(L7B_JSON))
    c025 = l7b["G2"]["counts"]["+0.25"]
    pi_l7b = {"bin_6.4_6.8": c025[16] / 2000.0, "bin_6.8_7.2": c025[17] / 2000.0}
    tctrl = tab_ref["%+.2f" % NY_CTRL]
    n_fresh = tctrl["all_bins"]["n_fresh"]
    fresh_frac = n_fresh / float(NT_TOT * NZ_DET)
    flux_fresh = {k: tctrl[k]["flux_fresh"]
                  for k in tctrl if k.startswith("bin_")}
    flux_fresh_tot = float(sum(flux_fresh.values()))
    pi_tot = float(sum(pi_l7b.values()))
    g5 = dict(n_fresh=n_fresh, fresh_frac_of_grid=fresh_frac,
              flux_fresh_per_bin=flux_fresh,
              l7b_measured_pi=pi_l7b,
              flux_fresh_total=flux_fresh_tot,
              l7b_pi_total=pi_tot,
              flux_over_pi=float(flux_fresh_tot / pi_tot),
              n_fresh_nx2=tab_nx2["%+.2f" % NY_CTRL]["all_bins"]["n_fresh"],
              n_fresh_dt2=tab_dt2["%+.2f" % NY_CTRL]["all_bins"]["n_fresh"],
              classification=tctrl["all_bins"])
    g5["verdict"] = "PASS" if (n_fresh >= 50
                               and flux_fresh_tot >= 0.3 * pi_tot) else "FAIL"

    # ------------------------------------------------------------------
    results = dict(
        workstream="O1", file_id="F-T7-O1",
        params=dict(L=L, sigma_x=SIGMA_X, x0=X0, k0=K0,
                    taper="S2 smoothstep, W=1 for |x-x0|<=2.5 sigma, 0 at 3.5",
                    amplitude="literal exp(-(x-x0)^2/(2 sigma^2)) * W",
                    box=[XMIN, XMAX], Nx_ref=NX_REF, dt_ref=DT_REF,
                    d=D_NEAR, kill_bins=KILL_BINS,
                    ny_kill=NY_KILL, ny_control=NY_CTRL,
                    nz_det=NZ_DET, nt_per_bin=NT_PER_BIN,
                    rho_floor_rel=RHO_FLOOR_REL, margin_eps=MARGIN_EPS,
                    margin_class=MARGIN_CLASS, margin_cap=MARGIN_CAP,
                    vclamp=VCLAMP, rho_eps_rel=RHO_EPS_REL,
                    s_start=S_START, dz_cell=DZ_CELL, dt_cell=DT_CELL),
        classification_ref=tab_ref,
        classification_nx2=tab_nx2,
        classification_dt2=tab_dt2,
        stability=stab,
        G1=g1, G2=g2, G3=g3, G4=g4, G5=g5,
        run_elapsed={lab: runs[lab]["elapsed"] for lab in runs},
        elapsed_total=time.time() - t0_wall,
    )
    with open(os.path.join(OUTDIR, "o1_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    histS, histX = runs["ref"]["hist"]
    np.savez_compressed(
        os.path.join(OUTDIR, "o1_raw.npz"),
        ny_flat=NY_FLAT, t_flat=T_FLAT, z_flat=Z_FLAT, zc=ZC, tc=TC,
        cls_ref=ref["cls"], cls_nx2=recs["nx2"]["cls"],
        cls_dt2=recs["dt2"]["cls"],
        margin_ref=ref["margin"], margin_nx2=recs["nx2"]["margin"],
        margin_dt2=recs["dt2"]["margin"],
        s_pc_ref=ref["s_pc"], s_floor_ref=ref["s_floor"],
        Xend_ref=ref["Xend"], Zend_ref=ref["Zend"],
        vx_det_ref=ref["vx_det"], rho_det_ref=ref["rho_det"],
        flux_ref=ref["flux"],
        hist_idx=hist_idx, histS=histS, histX=histX)
    print("verdicts:", {g: results[g]["verdict"]
                        for g in ("G1", "G2", "G3", "G4", "G5")})
    print("total elapsed %.1f s" % results["elapsed_total"])


if __name__ == "__main__":
    main()
