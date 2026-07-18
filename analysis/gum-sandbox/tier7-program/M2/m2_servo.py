#!/usr/bin/env python3
"""
M2 — the bootstrap servo coefficient (IV.H.3 <r2> toy law).

Minimal faithful toy of GUM IV.H.3's phase-locked knot self-oscillator:

  Reduced pair (Dphi, delta), canonical under
      H_S = omega* U(delta) + g^2 omega* (1 - cos Dphi),
      U'(delta) = f(delta) = (1+delta)^{1/2} - (1+delta)^{-1/2}
  so that
      Dphi_dot = omega* f(delta)          (exact +-1/2 clock exponents)
      delta_dot = -g^2 omega* sin Dphi + F_bath        (Adler torque, -sin)

  Bath: N harmonic modes (x_k, p_k) bilinearly coupled to the ANGLE Dphi
  (standard Caldeira-Leggett with counterterm), so the bath force appears
  in the delta_dot equation.  Couplings c_k ~ g_chi: both the coherent
  lock torque (g^2) and the incoherent friction kernel (c_k^2 ~ g^2) are
  second order in the knot-medium coupling.
      J(omega) = (pi/2) sum_k (c_k^2/omega_k) delta(omega - omega_k)
      Ohmic:      J = eta  g^2 omega            exp(-omega/omega_c)
      super-Ohmic:J = eta3 g^2 omega^3/omega_c^2 exp(-omega/omega_c)
  omega_c = 5 omega* (spec-fixed).  Classical thermal initial conditions
  at small T; the "noise" is bath initial data (closed Hamiltonian
  system -> total-H conservation is an integrator gate).

  Weak-coupling linear response (derivation in RESULTS.md):
      lambda = omega* J(Omega_lib) / (2 Omega_lib),  Omega_lib = g omega*
      Ohmic  -> lambda = (eta/2) g^2 omega* exp(-g omega*/omega_c)
      cubic  -> lambda = (eta3/50) g^4 omega* exp(-g omega*/omega_c)
  so c = lambda/(g^2 omega*) = eta/2 for Ohmic (pure bath convention),
  and the g^2 law itself FAILS for the super-Ohmic shape (lambda ~ g^4).

Units hbar = m = 1.  Deterministic given the fixed seeds.
Writes m2_results.json; figure produced by m2_fig() (m2_fig.png).
"""

import json
import time

import numpy as np

# ----------------------------------------------------------------------
# model pieces
# ----------------------------------------------------------------------

DELTA0 = 0.05          # spec-fixed initial displacement of the order parameter
WMIN_ABS = 1.0e-3      # absolute low end of bath grid (kept absolute on purpose:
                       # the omega* rescaling gate then probes a real asymmetry)


def f_beat(delta):
    """Exact clock mismatch: (1+d)^{1/2} - (1+d)^{-1/2}  (the +-1/2 exponents)."""
    s = np.sqrt(1.0 + delta)
    return s - 1.0 / s


def U_pot(delta):
    """Potential with U' = f_beat, U(0) = 0;  U ~ delta^2/2 for small delta."""
    s = np.sqrt(1.0 + delta)
    return (2.0 / 3.0) * (s ** 3 - 1.0) - 2.0 * (s - 1.0)


def J_spec(w, shape, eta, g, wc):
    if shape == "ohmic":
        return eta * g * g * w * np.exp(-w / wc)
    if shape == "cubic":
        return eta * g * g * (w ** 3 / wc ** 2) * np.exp(-w / wc)
    raise ValueError(shape)


def lam_pred(g, wstar, eta, shape):
    """Weak-coupling linear-response prediction lambda = w* J(Om)/(2 Om)."""
    om = g * wstar
    wc = 5.0 * wstar
    if shape == "ohmic":
        return 0.5 * eta * g * g * wstar * np.exp(-om / wc)
    if shape == "cubic":
        return 0.5 * eta * g * g * wstar * (om / wc) ** 2 * np.exp(-om / wc)
    raise ValueError(shape)


def lam_lin_exact(g, wstar, eta, shape, tol=1e-12):
    """Exact linear-response decay rate for the CONTINUUM bath: complex root of
        D(s) = s^2 + Omega^2 + wstar * s^2 * K(s) = 0,
        K(s) = (2/pi) int_0^inf J(w) / (w (s^2 + w^2)) dw
    near s = -lam + i Omega_R.  Captures cutoff suppression AND the
    counterterm-included frequency (Lamb) renormalization, i.e. everything the
    linearized model predicts beyond the weak-coupling formula."""
    from scipy.integrate import quad
    om = g * wstar
    wc = 5.0 * wstar

    def K(s):
        def integ(u, part):
            val = J_spec(u, shape, eta, g, wc) / (u * (s * s + u * u))
            return val.real if part == 0 else val.imag
        re = quad(integ, 0, 60 * wc, args=(0,), limit=400,
                  points=[abs(s.imag), wc])[0]
        im = quad(integ, 0, 60 * wc, args=(1,), limit=400,
                  points=[abs(s.imag), wc])[0]
        return (2.0 / np.pi) * (re + 1j * im)

    def D(s):
        return s * s + om * om + wstar * s * s * K(s)

    s = -lam_pred(g, wstar, eta, shape) + 1j * om
    for _ in range(60):
        h = 1e-7 * max(abs(s), 1.0)
        dD = (D(s + h) - D(s - h)) / (2 * h)
        step = D(s) / dD
        s = s - step
        if abs(step) < tol * max(abs(s), 1.0):
            break
    return float(-s.real), float(abs(s.imag))


def build_bath(R, N, wmin, wmax, shape, eta, g, wc, rng):
    """Log grid, one mode per cell, frequency jittered per realization.
    c_k^2 = (2/pi) omega_k J(omega_k) Delta_k  (J-locked couplings)."""
    edges = np.geomspace(wmin, wmax, N + 1)
    lo, hi = edges[:-1], edges[1:]
    w = lo + rng.uniform(size=(R, N)) * (hi - lo)
    c2 = (2.0 / np.pi) * w * J_spec(w, shape, eta, g, wc) * (hi - lo)
    return w, np.sqrt(c2)


# ----------------------------------------------------------------------
# integrators
# ----------------------------------------------------------------------

def run_config(g, wstar, eta, shape, dt, t_end, R, N, T, seed,
               bath_on=True, nrec_target=1500):
    """Velocity-Verlet (kick-drift-kick) on the full closed system.
    Returns records of (t, Dphi, delta, E_S) for all R realizations plus
    total-H drift diagnostics."""
    rng = np.random.default_rng(seed)
    wc = 5.0 * wstar
    g2ws = g * g * wstar

    if bath_on:
        w, c = build_bath(R, N, WMIN_ABS, 25.0 * wstar, shape, eta, g, wc, rng)
        w2 = w * w
        S2 = (c * c / w2).sum(axis=1)                       # counterterm curvature
        x = rng.standard_normal((R, N)) * (np.sqrt(T) / w)  # thermal ICs about phi=0
        p = rng.standard_normal((R, N)) * np.sqrt(T)
    phi = np.zeros(R)
    dlt = np.full(R, DELTA0)

    nsteps = int(round(t_end / dt))
    stride = max(1, nsteps // nrec_target)
    nrec = nsteps // stride
    t_rec = np.empty(nrec)
    phi_r = np.empty((nrec, R))
    dlt_r = np.empty((nrec, R))
    es_r = np.empty((nrec, R))
    h_r = np.empty((nrec, R))

    def E_sys(phi, dlt):
        return wstar * U_pot(dlt) + g2ws * (1.0 - np.cos(phi))

    def H_tot(phi, dlt, x, p):
        h = E_sys(phi, dlt)
        h = h + 0.5 * (p * p + w2 * x * x).sum(axis=1)
        h = h - phi * (c * x).sum(axis=1) + 0.5 * S2 * phi * phi
        return h

    # cached forces
    if bath_on:
        Fphi = -g2ws * np.sin(phi) + (c * x).sum(axis=1) - phi * S2
        Fx = -w2 * x + c * phi[:, None]
    else:
        Fphi = -g2ws * np.sin(phi)

    irec = 0
    for step in range(nsteps):
        # half kick
        dlt += 0.5 * dt * Fphi
        if bath_on:
            p += 0.5 * dt * Fx
        # drift
        phi += dt * wstar * f_beat(dlt)
        if bath_on:
            x += dt * p
        # new forces + half kick
        if bath_on:
            cx = (c * x).sum(axis=1)
            Fphi = -g2ws * np.sin(phi) + cx - phi * S2
            Fx = -w2 * x + c * phi[:, None]
            p += 0.5 * dt * Fx
        else:
            Fphi = -g2ws * np.sin(phi)
        dlt += 0.5 * dt * Fphi

        if (step + 1) % stride == 0 and irec < nrec:
            t_rec[irec] = (step + 1) * dt
            phi_r[irec] = phi
            dlt_r[irec] = dlt
            es_r[irec] = E_sys(phi, dlt)
            h_r[irec] = H_tot(phi, dlt, x, p) if bath_on else es_r[irec]
            irec += 1

    h0 = h_r[0]
    drift_max = float(np.abs((h_r - h0) / h0).max())
    # secular trend of H (linear fit, mean over realizations)
    hm = h_r.mean(axis=1)
    sec = np.polyfit(t_rec, (hm - hm[0]) / abs(hm[0]), 1)[0] * (t_rec[-1] - t_rec[0])
    return dict(t=t_rec, phi=phi_r, dlt=dlt_r, es=es_r,
                H_drift_max=drift_max, H_drift_secular=float(sec))


def run_pair_rk4(g, wstar, dt, t_end, nrec_target=2000):
    """Deterministic (bath-off) pair, classic RK4 — G1 energy gate + G5."""
    g2ws = g * g * wstar

    def rhs(y):
        phi, dlt = y
        return np.array([wstar * f_beat(dlt), -g2ws * np.sin(phi)])

    nsteps = int(round(t_end / dt))
    stride = max(1, nsteps // nrec_target)
    y = np.array([0.0, DELTA0])
    nrec = nsteps // stride
    t_rec = np.empty(nrec)
    phi_r = np.empty(nrec)
    dlt_r = np.empty(nrec)
    irec = 0
    for step in range(nsteps):
        k1 = rhs(y)
        k2 = rhs(y + 0.5 * dt * k1)
        k3 = rhs(y + 0.5 * dt * k2)
        k4 = rhs(y + dt * k3)
        y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        if (step + 1) % stride == 0 and irec < nrec:
            t_rec[irec] = (step + 1) * dt
            phi_r[irec], dlt_r[irec] = y
            irec += 1
    es = wstar * U_pot(dlt_r) + g2ws * (1.0 - np.cos(phi_r))
    e0 = wstar * U_pot(DELTA0)
    drift_max = float(np.abs(es / e0 - 1.0).max())
    sec = np.polyfit(t_rec, es / e0 - 1.0, 1)[0] * (t_rec[-1] - t_rec[0])
    return dict(t=t_rec, phi=phi_r[:, None], dlt=dlt_r[:, None], es=es[:, None],
                E_drift_max=drift_max, E_drift_secular=float(sec))


# ----------------------------------------------------------------------
# envelope fit
# ----------------------------------------------------------------------

def _fit_once(t, Q, floor, lam, wlo, whi):
    mask = (t >= wlo / lam) & (t <= whi / lam) & (Q - floor > 0)
    if mask.sum() < 8:
        return None
    y = np.log(Q[mask] - floor)
    s, b = np.polyfit(t[mask], y, 1)
    return -0.5 * s


def fit_lambda(rec, wstar, g, lam_guess, floor, wlo=0.2, whi=2.2, nblocks=6):
    """Envelope rate of delta(t): fit ln(Qbar - floor) with the ripple-reduced
    envelope Q = E_S + lambda * Dphi * delta (exact e^{-2 lambda t} envelope of
    the damped linear pair).  Iterated window; block errors over realizations;
    window-variation systematic."""
    t = rec["t"]

    def qbar(lam, idx=None):
        phi = rec["phi"] if idx is None else rec["phi"][:, idx]
        dlt = rec["dlt"] if idx is None else rec["dlt"][:, idx]
        es = rec["es"] if idx is None else rec["es"][:, idx]
        return (es + lam * phi * dlt).mean(axis=1)

    lam_win = lam_guess          # window scale (always the last POSITIVE rate)
    lam_c = None                 # measured rate (may be ~0 or negative:
    for it in range(4):          #  marginal/integrable case — report as is)
        lam_q = lam_win if lam_c is None else lam_c   # ripple-correction rate
        out = _fit_once(t, qbar(lam_q), floor, lam_win, wlo, whi)
        if out is None:
            break
        lam_c = out
        if out > 0.1 * lam_guess:
            lam_win = out
    if lam_c is None:
        lam_c = float("nan")

    # block statistical error
    R = rec["phi"].shape[1]
    lam_b = []
    if R >= nblocks:
        for j in range(nblocks):
            idx = np.arange(R)[j::nblocks]
            lb = _fit_once(t, qbar(lam_c, idx), floor, lam_win, wlo, whi)
            if lb is not None:
                lam_b.append(lb)
    sig_stat = float(np.std(lam_b) / np.sqrt(len(lam_b))) if len(lam_b) > 2 else 0.0

    # window systematic
    alts = []
    for fac in (0.7, 1.4):
        la = _fit_once(t, qbar(lam_c), floor, lam_win, wlo * fac, whi * fac)
        if la is not None:
            alts.append(la)
    sig_sys = 0.5 * (max(alts + [lam_c]) - min(alts + [lam_c])) if alts else 0.0
    return dict(lam=float(lam_c), sig_stat=sig_stat, sig_sys=float(sig_sys),
                sig=float(np.hypot(sig_stat, sig_sys)))


# ----------------------------------------------------------------------
# bath statistics / FDT verification (gate G1b)
# ----------------------------------------------------------------------

def bath_stats(seed=901, g=0.2, wstar=1.0, eta=1.0, T=1e-7, N=2000, M=5000):
    """Sample the free-bath force xi(t) = sum c_k [x_k(0) cos w_k t
    + p_k(0) sin(w_k t)/w_k] over M thermal draws; verify Gaussianity and the
    classical FDT  <xi(t) xi(0)> = T * gamma(t),
    gamma(t) = (2/pi) int J(w)/w cos(wt) dw = (2 eta g^2/pi) wc/(1+wc^2 t^2)."""
    from scipy import stats as sst
    rng = np.random.default_rng(seed)
    wc = 5.0 * wstar
    # fixed grid (no jitter) so the discrete kernel is exact and comparable
    edges = np.geomspace(WMIN_ABS, 25.0 * wstar, N + 1)
    w = np.sqrt(edges[:-1] * edges[1:])
    c = np.sqrt((2.0 / np.pi) * w * J_spec(w, "ohmic", eta, g, wc)
                * (edges[1:] - edges[:-1]))
    x0 = rng.standard_normal((M, N)) * (np.sqrt(T) / w)
    p0 = rng.standard_normal((M, N)) * np.sqrt(T)

    tgrid = np.linspace(0.0, 4.0 / wstar, 41)
    xi = np.empty((len(tgrid), M))
    for i, tt in enumerate(tgrid):
        xi[i] = x0 @ (c * np.cos(w * tt)) + p0 @ (c * np.sin(w * tt) / w)
    xi0 = xi[0]
    C_meas = (xi * xi0).mean(axis=1)
    gam_cont = (2.0 * eta * g * g / np.pi) * wc / (1.0 + (wc * tgrid) ** 2)
    gam_disc = np.array([((c * c / (w * w)) * np.cos(w * tt)).sum() for tt in tgrid])

    # Gaussianity of xi at t=0 and one interior time
    tests = {}
    for lab, s in (("t0", xi0), ("t_mid", xi[len(tgrid) // 2])):
        z = (s - s.mean()) / s.std()
        tests[lab] = dict(
            skew=float(sst.skew(s)), kurtosis_excess=float(sst.kurtosis(s)),
            se_skew=float(np.sqrt(6.0 / M)), se_kurt=float(np.sqrt(24.0 / M)),
            ks_p=float(sst.kstest(z, "norm").pvalue))
    # FDT match, normalized by peak
    dev_cont = float(np.abs(C_meas - T * gam_cont).max() / (T * gam_cont[0]))
    dev_disc = float(np.abs(C_meas - T * gam_disc).max() / (T * gam_disc[0]))
    samp_err = float(np.sqrt(2.0 / M))  # expected relative sampling error at peak
    return dict(tgrid=tgrid.tolist(), C_meas=C_meas.tolist(),
                Tgam_cont=(T * gam_cont).tolist(), Tgam_disc=(T * gam_disc).tolist(),
                gauss=tests, fdt_maxdev_over_peak_cont=dev_cont,
                fdt_maxdev_over_peak_disc=dev_disc,
                expected_sampling_err=samp_err,
                gamma0_disc_over_cont=float(gam_disc[0] / gam_cont[0]), M=M, N=N)


# ----------------------------------------------------------------------
# campaign driver
# ----------------------------------------------------------------------

def main(quick=False):
    t0 = time.time()
    res = {"meta": dict(
        delta0=DELTA0, N_modes=2000, R_realizations=24, T_over_wstar=1e-7,
        omega_c_over_wstar=5.0, wmin_abs=WMIN_ABS, wmax_over_wstar=25.0,
        dt_main=0.02, seed_base=1000,
        architecture="canonical pair (Dphi, delta); CL bath on Dphi (counterterm); c_k ~ g_chi")}
    N, R, T = 2000, 24, 1e-7
    if quick:
        N, R = 400, 8

    def one(g, wstar, eta, shape, dt, seed, span=2.6, Tloc=None,
            Nloc=None, Rloc=None):
        lp = lam_pred(g, wstar, eta, shape)
        lam_lin, om_r = lam_lin_exact(g, wstar, eta, shape)
        Tl = (T * wstar) if Tloc is None else Tloc
        rec = run_config(g, wstar, eta, shape, dt, span / lam_lin,
                         Rloc or R, Nloc or N, Tl, seed)
        fit = fit_lambda(rec, wstar, g, lam_lin, floor=Tl)
        fit.update(lam_pred=float(lp), lam_lin_exact=float(lam_lin),
                   Omega_R=float(om_r),
                   H_drift_max=rec["H_drift_max"],
                   H_drift_secular=rec["H_drift_secular"], g=g, wstar=wstar,
                   eta=eta, shape=shape, dt=dt, t_end=float(span / lam_lin))
        return rec, fit

    # ---- A. Ohmic g-scan (G2, G4) --------------------------------------
    print("== A: Ohmic scan eta=1 ==", flush=True)
    gs = np.geomspace(0.06, 0.6, 8)
    scan = []
    rec_example = None
    for i, g in enumerate(gs):
        rec, fit = one(float(g), 1.0, 1.0, "ohmic", 0.02, 1000 + i)
        scan.append(fit)
        if i == 4:                  # keep one mid-scan envelope for the figure
            rec_example = dict(t=rec["t"].tolist(),
                               Qbar=(rec["es"] + fit["lam"] * rec["phi"] * rec["dlt"])
                               .mean(axis=1).tolist(),
                               dlt0=rec["dlt"][:, 0].tolist(), g=float(g),
                               lam=fit["lam"])
        print(f"  g={g:.4f}  lam={fit['lam']:.6e} (lin {fit['lam_lin_exact']:.3e},"
              f" wk {fit['lam_pred']:.3e}) +- {fit['sig']:.1e} "
              f" Hdrift={fit['H_drift_max']:.2e}", flush=True)
    res["scan_ohmic"] = scan
    res["example_envelope"] = rec_example

    # G2: log-log slope
    lg = np.log(gs)
    ll = np.log([s["lam"] for s in scan])
    sig_ll = np.array([max(s["sig"], 1e-12) / s["lam"] for s in scan])
    wts = 1.0 / sig_ll ** 2
    Sw, Sx, Sy = wts.sum(), (wts * lg).sum(), (wts * ll).sum()
    Sxx, Sxy = (wts * lg * lg).sum(), (wts * lg * ll).sum()
    den = Sw * Sxx - Sx * Sx
    slope = (Sw * Sxy - Sx * Sy) / den
    sig_slope = np.sqrt(Sw / den)
    # cutoff-corrected slope: remove the known exp(-g/5) of the spec-fixed cutoff
    ll_c = ll + gs / 5.0
    slope_c = ((Sw * (wts * lg * ll_c).sum() - Sx * (wts * ll_c).sum()) / den)
    res["G2"] = dict(slope=float(slope), sig_slope=float(sig_slope),
                     slope_cutoff_corrected=float(slope_c),
                     predicted_raw_slope="2 - <g> w*/w_c (cutoff drift)",
                     gs=gs.tolist(), lam=[s["lam"] for s in scan],
                     sig=[s["sig"] for s in scan])
    print(f"G2 slope = {slope:.4f} +- {sig_slope:.4f} (cutoff-corr {slope_c:.4f})")

    # G4: coefficient c = lam/(g^2 w*)
    c_i = np.array([s["lam"] / (s["g"] ** 2) for s in scan])
    c0_i = c_i * np.exp(gs / 5.0)     # cutoff-corrected -> should be eta/2
    cbar = float(c_i.mean())
    sig_c = float(np.hypot(c_i.std(ddof=1) / np.sqrt(len(c_i)),
                           0.5 * (c_i.max() - c_i.min()) / np.sqrt(3)))
    res["G4_ohmic_eta1"] = dict(c_per_g=c_i.tolist(), c_cutoffcorr_per_g=c0_i.tolist(),
                                c_mean=cbar, c_sig=sig_c,
                                c_cutoffcorr_mean=float(c0_i.mean()),
                                eta_over_2=0.5)
    print(f"G4 c = {cbar:.4f} +- {sig_c:.4f}  (cutoff-corr {c0_i.mean():.4f}; eta/2 = 0.5)")

    # ---- B. omega* rescaling x4 (G3) + dt cross-check ------------------
    print("== B: omega* invariance ==", flush=True)
    g_inv = float(gs[4])            # shares g with the scan -> dt cross-check
    _, fit_w1 = one(g_inv, 1.0, 1.0, "ohmic", 0.005, 2001)
    _, fit_w4 = one(g_inv, 4.0, 1.0, "ohmic", 0.005, 2002)
    r_inv = (fit_w4["lam"] / 4.0) / fit_w1["lam"]
    sig_r = r_inv * np.hypot(fit_w1["sig"] / fit_w1["lam"], fit_w4["sig"] / fit_w4["lam"])
    res["G3"] = dict(g=g_inv, lam_w1=fit_w1["lam"], sig_w1=fit_w1["sig"],
                     lam_w4=fit_w4["lam"], sig_w4=fit_w4["sig"],
                     ratio_lam_over_wstar=float(r_inv), sig_ratio=float(sig_r))
    res["dt_check"] = dict(g=g_inv, lam_dt02=scan[4]["lam"], sig_dt02=scan[4]["sig"],
                           lam_dt005=fit_w1["lam"], sig_dt005=fit_w1["sig"],
                           rel_diff=float(fit_w1["lam"] / scan[4]["lam"] - 1.0))
    print(f"G3 ratio (lam/w*)_4 / (lam/w*)_1 = {r_inv:.4f} +- {sig_r:.4f}; "
          f"dt check rel diff {res['dt_check']['rel_diff']:+.3f}")

    # ---- C. eta linearity (coefficient is pure bath convention) --------
    # smaller eta -> narrower linewidth -> denser bath needed (modes within
    # the width lam*rho(Omega) must stay >= a few)
    print("== C: eta=0.1 ==", flush=True)
    _, fit_eta = one(0.2, 1.0, 0.1, "ohmic", 0.02, 3001,
                     Nloc=(N if quick else 8000), Rloc=(R if quick else 12))
    res["eta_linearity"] = dict(eta=0.1, lam=fit_eta["lam"], sig=fit_eta["sig"],
                                lam_lin=fit_eta["lam_lin_exact"],
                                c=fit_eta["lam"] / 0.04,
                                c_over_eta=fit_eta["lam"] / 0.04 / 0.1)
    print(f"  c(eta=0.1) = {fit_eta['lam']/0.04:.4f}  (eta/2 = 0.05; "
          f"lin {fit_eta['lam_lin_exact']/0.04:.4f})")

    # ---- D. super-Ohmic shape (G4 spread) ------------------------------
    print("== D: super-Ohmic ==", flush=True)
    ETA3 = 6.0
    so = []
    for i, g in enumerate([0.35, 0.475, 0.65]):
        _, fit_so = one(g, 1.0, ETA3, "cubic", 0.02, 4001 + i,
                        Nloc=(N if quick else 8000), Rloc=(R if quick else 12))
        fit_so["c"] = fit_so["lam"] / (g * g)
        so.append(fit_so)
        print(f"  g={g}  lam={fit_so['lam']:.4e} (lin {fit_so['lam_lin_exact']:.3e},"
              f" wk {fit_so['lam_pred']:.3e})  c={fit_so['c']:.4f}", flush=True)
    lgso = np.log([s["g"] for s in so])
    llso = np.log([s["lam"] for s in so])
    slope_so = float(np.polyfit(lgso, llso, 1)[0])
    slope_so_lin = float(np.polyfit(lgso, np.log([s["lam_lin_exact"] for s in so]), 1)[0])
    res["superohmic"] = dict(eta3=ETA3, points=so, slope=slope_so,
                             slope_lin_exact=slope_so_lin,
                             note="lambda ~ g^4-class: the g^2 law fails off the Ohmic class")
    print(f"  super-Ohmic slope = {slope_so:.3f} (exact-linear {slope_so_lin:.3f})")

    # ---- E/F. integrable limit (G5) + deterministic energy gate (G1a) --
    print("== E: integrable control ==", flush=True)
    det_rk4 = run_pair_rk4(0.2, 1.0, 0.01, 2000.0)
    lp02 = lam_pred(0.2, 1.0, 1.0, "ohmic")
    fit_int = fit_lambda(det_rk4, 1.0, 0.2, lp02, floor=0.0)
    # whole-run envelope rate as well
    fit_int_long = _fit_once(det_rk4["t"], det_rk4["es"][:, 0], 0.0,
                             2.2 / det_rk4["t"][-1], 0.2, 2.2)
    rec_v = run_config(0.2, 1.0, 1.0, "ohmic", 0.02, 2000.0, 1, 0, 0.0, 5001,
                       bath_on=False)
    fit_int_v = fit_lambda(rec_v, 1.0, 0.2, lp02, floor=0.0)
    res["G5"] = dict(
        lam_rk4_windowfit=fit_int["lam"], sig_rk4=fit_int["sig"],
        lam_rk4_wholerun=float(fit_int_long) if fit_int_long else None,
        lam_verlet_windowfit=fit_int_v["lam"], sig_verlet=fit_int_v["sig"],
        lam_bath_same_g=lp02,
        ratio_floor=abs(fit_int["lam"]) / lp02,
        env_t=rec_v["t"][::4].tolist(), env_E=rec_v["es"][::4, 0].tolist())
    res["G1_energy"] = dict(
        rk4_drift_max=det_rk4["E_drift_max"], rk4_drift_secular=det_rk4["E_drift_secular"],
        verlet_drift_max=rec_v["H_drift_max"], verlet_drift_secular=rec_v["H_drift_secular"],
        bath_run_H_drift_max=max(s["H_drift_max"] for s in scan),
        bath_run_H_drift_secular=max(abs(s["H_drift_secular"]) for s in scan))
    print(f"G5 |lam_integrable| = {abs(fit_int['lam']):.2e} (RK4), "
          f"{abs(fit_int_v['lam']):.2e} (Verlet) vs bath {lp02:.2e}")
    print(f"G1a RK4 drift max {det_rk4['E_drift_max']:.2e} "
          f"secular {det_rk4['E_drift_secular']:.2e}")

    # ---- G. bath statistics (G1b) --------------------------------------
    print("== G: bath statistics ==", flush=True)
    res["G1_bath_stats"] = bath_stats()
    bs = res["G1_bath_stats"]
    print(f"  FDT maxdev/peak = {bs['fdt_maxdev_over_peak_disc']:.3f} "
          f"(sampling {bs['expected_sampling_err']:.3f}); "
          f"skew {bs['gauss']['t0']['skew']:+.3f} kurt "
          f"{bs['gauss']['t0']['kurtosis_excess']:+.3f} KS p {bs['gauss']['t0']['ks_p']:.2f}")

    # ---- G4 verdict arithmetic -----------------------------------------
    sig_comb = float(np.hypot(sig_c, 0.004))
    z = abs(cbar - 0.021) / sig_comb
    res["G4_verdict"] = dict(c_corpus=0.021, sig_corpus=0.004,
                             c_measured=cbar, sig_measured=sig_c,
                             sigma_combined=sig_comb, z=z,
                             eta_required_for_match=2 * 0.021,
                             verdict="PASS-MATCH" if z <= 2 else "PASS-MEASURED")
    print(f"G4: z = {z:.1f} sigma_combined -> {res['G4_verdict']['verdict']}")

    res["meta"]["wall_seconds"] = time.time() - t0
    with open("m2_results.json", "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"wall {time.time()-t0:.1f}s -> m2_results.json")
    return res


if __name__ == "__main__":
    import sys
    main(quick="--quick" in sys.argv)
