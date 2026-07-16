#!/usr/bin/env python3
"""Phase H2.7(b) — the <r2>-class Adler locking-torque sign in the
thermalizing-bath class vs the integrable limit (ROADMAP_v5 item 7),
plus assembly of h27_results.json (merging the H2.7(a) field-level
+-1/2 exponents produced by h26_solve/h26_convention).

THE TOY (minimal, per task: "a 10-line ODE with the same symbol-matrix
coefficients — simpler is better"), built from the corpus's own reduced
objects:

  * The knot is the closure family's OWN nonlinear oscillator in
    action-angle form: E(J) = e*sqrt(1 + J^2/(e*i)) — exactly the exact
    one-parameter family of Omega-paper Sec. IV / T2 §0 with L -> J,
    e = e0 = 64/15pi, i = i0 = 256/105pi (shape factor g = 1; the choice
    only scales J*, omega*).  omega(J) = dE/dJ; closure fixed point
    x = J^2/(e*i) = 1: J* = sqrt(e*i), omega* = sqrt(e/2i) = sqrt(7/8)
    (e0/i0 = 7/4 exact).  The action-dependent frequency IS the
    nonlinearity that makes phase-locking a pendulum problem.
  * The vacuum channels are TWO explicit modes with fs-gum-cosserat's
    VALIDATED benchmark eigenfrequencies (RESULTS.md, gate G19 table,
    moduli rho0 = J = lambda = mu = 1, mu_c = 5, alpha = beta = gamma_m
    = 0.5, m_V = 1): B2 light transverse at k = 1, omega_2 =
    0.947870595456; B4 longitudinal twist at k = 1, omega_4 =
    4.690415759823 (gap^2 = 21 + (alpha+beta)k^2).  These are the
    gapless/twist channels IV.H.3 names.
  * The vacuum consensus is an external reference oscillation
    B cos(omega_ref t) at omega_ref = omega* (the medium's hbar sets the
    clock; the knot is displaced off it via J(0)).
  * Coupling (bilinear, strength g_chi): H_c = -g_chi cos(theta) *
    (y2 + y4 + B cos(omega_ref t)).

  thetadot = omega(J)
  Jdot     = -g_chi sin(theta) (y2 + y4 + B cos(omega_ref t))
  y_k''    = -omega_k^2 y_k - gamma y_k' + g_chi cos(theta)

  gamma > 0  = the thermalizing-bath variant (the two channels shed the
               beat irreversibly — "add weak damping");
  gamma = 0  = the integrable limit (fully conservative flow).

MEASUREMENTS
  M0  RK4 dissipation control (g = 0, gamma = 0): the artificial damping
      floor of the integrator, so "marginal" is separated from numerics.
  M1  Open-loop locking torque: clamp theta = omega_ref t + dphi, average
      <Jdot>(dphi) over many periods, fit T0 + Ts sin + Tc cos.  The
      Adler form and its SIGN (corpus: torque ∝ -sin dphi).
  M2  Closed-loop: full flow from J(0) = J*(1 + dJ); measure the
      locking rate lambda as the exponential decay rate of the beat
      u(t) = omega(J) - omega_ref (windowed-RMS envelope fit), for a
      gamma scan including gamma = 0.  Power-law fit lambda ∝ gamma^p
      pins the class boundary.
  M3  g_chi scan at fixed gamma: lambda ∝ g_chi^2 (the <r2> toy law
      lambda = 0.021 g_chi^2 omega* class).

Deterministic: no RNG anywhere; fixed dt, fixed windows; same
interpreter => identical JSON to the printed digits.

Epistemic notice (binding): a within-model toy of a speculative theory's
locking argument.  Nothing here says anything about nature.
"""

import json
import math
import os

E0 = 64.0 / (15.0 * math.pi)   # e0hat
I0 = 256.0 / (105.0 * math.pi)  # i0hat
JSTAR = math.sqrt(E0 * I0)
WSTAR = math.sqrt(E0 / (2.0 * I0))          # = sqrt(7/8), fixed point x = 1
WPRIME = WSTAR / (2.0 * JSTAR)               # d omega/dJ at x = 1 (dlnw/dlnJ = 1/2)
G_CHI = 0.05                                 # default coupling
B_REF = 1.0
C_B2 = 1.0                                   # B2 speed, mu/rho0 = 1 (validated G05)


def olib(g):
    """Small-libration frequency of the Adler pendulum at coupling g."""
    return math.sqrt(WPRIME * g * B_REF / 2.0)


# The B2 channel is GAPLESS (validated: gap <= 1e-12, G04; speed c^2 =
# mu/rho0, G05), so the beat radiates on-branch at the clock's libration
# sideband: omega_2 = c_B2 * k_2 with k_2 = (omega* + Olib)/c_B2.  The B4
# twist channel is taken at the validated k = 1 benchmark eigenfrequency
# (gap^2 = 21 + (alpha+beta)k^2, G08/G19) — far off-resonant, reactive.
W2_BASE = WSTAR + olib(G_CHI)
W4 = 4.690415759823                          # fs-gum-cosserat B4 @ k=1 (validated)
WREF = WSTAR
DT = 0.025

HERE = os.path.dirname(os.path.abspath(__file__))


def omega_of(j):
    x = j * j / (E0 * I0)
    return j / (I0 * math.sqrt(1.0 + x))


def rhs(t, s, g, gam, w2):
    th, j, y2, p2, y4, p4 = s
    c = math.cos(th)
    drive = B_REF * math.cos(WREF * t)
    return (
        omega_of(j),
        -g * math.sin(th) * (y2 + y4 + drive),
        p2,
        -w2 * w2 * y2 - gam * p2 + g * c,
        p4,
        -W4 * W4 * y4 - gam * p4 + g * c,
    )


def rk4_run(s0, g, gam, tmax, sample_every, clamp=False, w2=W2_BASE):
    """Fixed-step RK4.  clamp=True freezes (theta = WREF t + dphi0, J)
    — the open-loop torque protocol (dphi0 is encoded in s0[0])."""
    s = list(s0)
    dphi0 = s0[0]
    t = 0.0
    n = int(round(tmax / DT))
    stride = max(1, int(round(sample_every / DT)))
    out_t, out_j, out_jdot, out_th = [], [], [], []
    for it in range(n + 1):
        if clamp:
            s[0] = WREF * t + dphi0
        if it % stride == 0:
            k = rhs(t, s, g, gam, w2)
            out_t.append(t)
            out_j.append(s[1])
            out_jdot.append(k[1])
            out_th.append(s[0])
        k1 = rhs(t, s, g, gam, w2)
        s1 = [s[i] + 0.5 * DT * k1[i] for i in range(6)]
        if clamp:
            s1[0] = WREF * (t + 0.5 * DT) + dphi0
            s1[1] = s[1]
        k2 = rhs(t + 0.5 * DT, s1, g, gam, w2)
        s2 = [s[i] + 0.5 * DT * k2[i] for i in range(6)]
        if clamp:
            s2[0] = WREF * (t + 0.5 * DT) + dphi0
            s2[1] = s[1]
        k3 = rhs(t + 0.5 * DT, s2, g, gam, w2)
        s3 = [s[i] + DT * k3[i] for i in range(6)]
        if clamp:
            s3[0] = WREF * (t + DT) + dphi0
            s3[1] = s[1]
        k4 = rhs(t + DT, s3, g, gam, w2)
        for i in range(6):
            s[i] += (DT / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        if clamp:
            s[1] = s0[1]
        t += DT
    return out_t, out_j, out_jdot, out_th


def lstsq_3(ts, ys, w):
    """Least squares of ys on [1, sin(w x), cos(w x)] with x = ts."""
    a00 = len(ts)
    s1 = sum(math.sin(w * x) for x in ts)
    c1 = sum(math.cos(w * x) for x in ts)
    ss = sum(math.sin(w * x) ** 2 for x in ts)
    cc = sum(math.cos(w * x) ** 2 for x in ts)
    sc = sum(math.sin(w * x) * math.cos(w * x) for x in ts)
    b0 = sum(ys)
    b1 = sum(y * math.sin(w * x) for x, y in zip(ts, ys))
    b2 = sum(y * math.cos(w * x) for x, y in zip(ts, ys))
    # solve 3x3 (Cramer)
    m = [[a00, s1, c1], [s1, ss, sc], [c1, sc, cc]]
    b = [b0, b1, b2]

    def det3(mm):
        return (mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
                - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
                + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0]))

    d = det3(m)
    sol = []
    for col in range(3):
        mm = [row[:] for row in m]
        for r in range(3):
            mm[r][col] = b[r]
        sol.append(det3(mm) / d)
    return sol  # [T0, Ts, Tc]


def linfit(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    return slope, my - slope * mx


def m0_dissipation_control():
    """g = 0, gamma = 0, both bath modes rung with unit amplitude:
    measured RK4 energy decay rate (the numerics floor for 'marginal')."""
    s = [0.0, JSTAR, 1.0, 0.0, 1.0, 0.0]
    t = 0.0
    n = int(round(20000.0 / DT))
    e_start = None
    for it in range(n + 1):
        if it == 0:
            e_start = 0.5 * (s[3] ** 2 + W2_BASE ** 2 * s[2] ** 2 + s[5] ** 2 + W4 ** 2 * s[4] ** 2)
        k1 = rhs(t, s, 0.0, 0.0, W2_BASE)
        s1 = [s[i] + 0.5 * DT * k1[i] for i in range(6)]
        k2 = rhs(t + 0.5 * DT, s1, 0.0, 0.0, W2_BASE)
        s2 = [s[i] + 0.5 * DT * k2[i] for i in range(6)]
        k3 = rhs(t + 0.5 * DT, s2, 0.0, 0.0, W2_BASE)
        s3 = [s[i] + DT * k3[i] for i in range(6)]
        k4 = rhs(t + DT, s3, 0.0, 0.0, W2_BASE)
        for i in range(6):
            s[i] += (DT / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
        t += DT
    e_end = 0.5 * (s[3] ** 2 + W2_BASE ** 2 * s[2] ** 2 + s[5] ** 2 + W4 ** 2 * s[4] ** 2)
    gam_num = -0.5 * math.log(e_end / e_start) / 20000.0  # amplitude rate
    return {"tmax": 20000.0, "dt": DT, "E_start": e_start, "E_end": e_end,
            "gamma_numerical_amplitude_rate": gam_num}


def m1_torque_curve(gam, g, w2=W2_BASE):
    """Open-loop: clamp theta = WREF t + dphi, measure <Jdot> over the
    last 2/3 of the window, 16 phases, fit T0 + Ts sin + Tc cos."""
    tmax, warm = 3000.0, 1000.0
    phases, torques = [], []
    for idx in range(16):
        dphi = 2.0 * math.pi * idx / 16.0
        s0 = (dphi, JSTAR, 0.0, 0.0, 0.0, 0.0)
        ts, _js, jdots, _th = rk4_run(s0, g, gam, tmax, 0.5, clamp=True, w2=w2)
        acc = [jd for t, jd in zip(ts, jdots) if t >= warm]
        phases.append(dphi)
        torques.append(sum(acc) / len(acc))
    t0, tsin, tcos = lstsq_3(phases, torques, 1.0)
    resid = max(abs(y - (t0 + tsin * math.sin(x) + tcos * math.cos(x)))
                for x, y in zip(phases, torques))
    return {"gamma": gam, "g_chi": g, "phases": phases, "mean_torque": torques,
            "fit_T0": t0, "fit_Tsin": tsin, "fit_Tcos": tcos,
            "fit_max_resid": resid, "adler_prediction_Tsin": -g * B_REF / 2.0}


def m2_locking_run(gam, g, dj=0.02, tmax=60000.0, w2=W2_BASE):
    """Closed loop from J(0) = J*(1+dj); lambda from the windowed-RMS
    envelope of u(t) = omega(J) - omega_ref."""
    s0 = (0.0, JSTAR * (1.0 + dj), 0.0, 0.0, 0.0, 0.0)
    ts, js, _jd, ths = rk4_run(s0, g, gam, tmax, 0.5, w2=w2)
    us = [omega_of(j) - WREF for j in js]

    def windows(wlen):
        per = int(round(wlen / 0.5))
        wt, wr = [], []
        for k in range(int(tmax / wlen)):
            seg = us[k * per:(k + 1) * per]
            if not seg:
                break
            wt.append((k + 0.5) * wlen)
            wr.append(math.sqrt(sum(v * v for v in seg) / len(seg)))
        return wt, wr

    def fit(wlen):
        """-slope of ln(RMS envelope), on the decaying section if the
        envelope reaches its floor, else on all windows (marginal/slow)."""
        wt, wr = windows(wlen)
        ntail = max(2, len(wr) // 10)
        floor = sum(wr[-ntail:]) / ntail
        xs, ys, crossed = [], [], False
        for x, r in zip(wt[1:], wr[1:]):
            if r > 2.0 * floor and not crossed:
                xs.append(x)
                ys.append(math.log(r))
            elif xs:
                crossed = True
        if len(xs) < 5:
            if crossed and wlen > 100.0:
                return fit(wlen / 4.0)  # fast lock: refine windows
            xs = list(wt[1:])
            ys = [math.log(r) for r in wr[1:]]
        if len(xs) < 5:
            return None, wt, wr
        slope, _ = linfit(xs, ys)
        return -slope, wt, wr

    lam, wt, wr = fit(400.0)  # base window ~6.8 libration periods at g=0.05
    floor = sum(wr[-max(2, len(wr) // 10):]) / max(2, len(wr) // 10)
    env_ratio = wr[-1] / wr[0]
    # model-free effective rate over the whole run (robust when the decay
    # is non-exponential, e.g. the overdamped-absorber regime)
    lam_env = -math.log(env_ratio) / (wt[-1] - wt[0])
    # locked phase (circular mean of theta - WREF t over the last 10%)
    ntail = max(1, len(ts) // 10)
    cs = sum(math.cos(th - WREF * t) for th, t in zip(ths[-ntail:], ts[-ntail:]))
    sn = sum(math.sin(th - WREF * t) for th, t in zip(ths[-ntail:], ts[-ntail:]))
    dphi_lock = math.atan2(sn, cs)
    return {"gamma": gam, "g_chi": g, "dJ0": dj, "tmax": tmax, "omega_2": w2,
            "lambda_lock": lam, "lambda_envelope": lam_env, "rms_floor": floor,
            "envelope_last_over_first": env_ratio,
            "dphi_locked_tail": dphi_lock,
            "rms_windows_t": wt[:: max(1, len(wt) // 40)],
            "rms_windows": wr[:: max(1, len(wt) // 40)]}


def drag_analytic(gam, g, w2=W2_BASE):
    """Time-averaged bath drag on the clamped carrier (steady state):
    T0 = -sum_k g^2 gamma omega_ref / (2[(omega_k^2-omega_ref^2)^2
         + gamma^2 omega_ref^2])  (linear-response, exact for the toy)."""
    tot = 0.0
    for wk in (w2, W4):
        d = wk * wk - WREF * WREF
        tot += g * g * gam * WREF / (2.0 * (d * d + gam * gam * WREF * WREF))
    return -tot


def declock_boundaries(g):
    """The toy's locking window: locked iff |T0(gamma)| < gB/2.  The
    near-resonant B2 drag is non-monotone in gamma (rises ~ gamma, falls
    ~ 1/gamma once the linewidth exceeds the detuning), so there can be a
    declocked interval [g_lo, g_hi]; bisect both crossings."""
    target = g * B_REF / 2.0
    # peak of |T0| at gamma* ~ |D2|/omega_ref
    d2 = abs(W2_BASE * W2_BASE - WREF * WREF) / WREF
    if abs(drag_analytic(d2, g)) <= target:
        return None  # never declocks
    lo, hi = 1.0e-6, d2
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if abs(drag_analytic(mid, g)) > target:
            hi = mid
        else:
            lo = mid
    g_lo = math.sqrt(lo * hi)
    lo, hi = d2, 10.0
    for _ in range(80):
        mid = math.sqrt(lo * hi)
        if abs(drag_analytic(mid, g)) > target:
            lo = mid
        else:
            hi = mid
    g_hi = math.sqrt(lo * hi)
    return (g_lo, g_hi)


def main():
    out = {"phase": "H2.7 — entrainment +-1/2 (field level, part a) + Adler locking-torque sign / class boundary (toy, part b)",
           "date": "2026-07-16",
           "model": {
               "knot": "E(J) = e0*sqrt(1 + J^2/(e0*i0)), action-angle; fixed point x = 1",
               "e0": E0, "i0": I0, "J_star": JSTAR, "omega_star": WSTAR,
               "omega_star_closed_form": "sqrt(7/8) (e0/i0 = 7/4 exact)",
               "bath_modes": {"B2_sideband": W2_BASE,
                              "B2_choice": "gapless branch omega = c k (validated G04/G05) at the upper libration sideband k_2 = (omega* + Olib(g))/c, Olib = sqrt(omega'(J*) g B/2)",
                              "Olib_default": olib(G_CHI),
                              "B4_k1": W4,
                              "source": "fs-gum-cosserat RESULTS.md G05/G08/G19 (validated benchmark moduli)"},
               "coupling": "H_c = -g_chi cos(theta) (y2 + y4 + B cos(omega* t)), B = 1",
               "dt": DT, "integrator": "fixed-step RK4, no RNG"}}

    print("H2.7(b) Adler toy: J*=%.6f omega*=%.6f (sqrt(7/8)=%.6f), B2=%.6f B4=%.6f"
          % (JSTAR, WSTAR, math.sqrt(7.0 / 8.0), W2_BASE, W4))

    # M0
    m0 = m0_dissipation_control()
    print("M0 integrator floor: gamma_num = %.3e (amplitude rate, gamma=0,g=0)"
          % m0["gamma_numerical_amplitude_rate"])
    out["M0_dissipation_control"] = m0

    # M1 torque curves
    g = G_CHI
    m1 = []
    for gam in [0.0, 1.0e-2, 3.0e-2]:
        r = m1_torque_curve(gam, g)
        r["drag_analytic_T0"] = drag_analytic(gam, g)
        m1.append(r)
        print("M1 gamma=%g: torque fit T0=%+.3e Ts=%+.6e Tc=%+.3e (Adler pred Ts=%+.6e, analytic T0=%+.3e, max resid %.1e)"
              % (gam, r["fit_T0"], r["fit_Tsin"], r["fit_Tcos"],
                 r["adler_prediction_Tsin"], r["drag_analytic_T0"],
                 r["fit_max_resid"]))
    out["M1_openloop_torque"] = m1
    bounds = declock_boundaries(g)
    out["M1_locking_window"] = {
        "criterion": "|T0(gamma)| < g B/2 (drag vs maximal Adler torque)",
        "declocked_interval": list(bounds) if bounds else None,
        "drag_curve_note": "B2 sits at the libration sideband (omega_2 - omega* = %.4f); drag is small and the window criterion may never fire" % (W2_BASE - WREF)}
    if bounds:
        print("M1 locking window: declocked (drag > gB/2) for gamma in [%.4e, %.4e]" % bounds)

    # M2 gamma scan
    gammas = [0.0, 3.0e-5, 1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3, 1.0e-2, 3.0e-2, 1.0e-1]
    m2 = []
    for gam in gammas:
        r = m2_locking_run(gam, g)
        m2.append(r)
        lam = r["lambda_lock"]
        print("M2 gamma=%-8g lambda=%s lambda_env=%+.3e  env(last/first)=%.4f  dphi_lock=%+.3f"
              % (gam, ("%.3e" % lam) if lam is not None else "none",
                 r["lambda_envelope"],
                 r["envelope_last_over_first"], r["dphi_locked_tail"]))
    out["M2_gamma_scan"] = m2
    # power law on the small-gamma asymptotic decades (below the declocking
    # window), where lambda was measured and positive
    xs = [math.log(r["gamma"]) for r in m2
          if 0 < r["gamma"] <= 3.0e-3 and r["lambda_lock"] and r["lambda_lock"] > 0]
    ys = [math.log(r["lambda_lock"]) for r in m2
          if 0 < r["gamma"] <= 3.0e-3 and r["lambda_lock"] and r["lambda_lock"] > 0]
    if len(xs) >= 3:
        p, c = linfit(xs, ys)
        out["M2_powerlaw"] = {"p": p, "prefactor": math.exp(c),
                              "gammas_used": [math.exp(x) for x in xs],
                              "scope": "small-gamma asymptote, gamma <= 3e-3"}
        print("M2 power law (gamma <= 3e-3): lambda ~ %.3e * gamma^%.3f" % (math.exp(c), p))
    # integrable-limit marginality bound
    r0 = m2[0]
    out["M2_integrable_limit"] = {
        "lambda_measured": r0["lambda_lock"],
        "envelope_last_over_first": r0["envelope_last_over_first"],
        "bound": "envelope drift over the full run vs the M0 numerics floor"}

    # M3 g scan at fixed gamma, with the B2 mode retuned to each coupling's
    # own libration sideband (physical: the gapless branch has a mode at
    # every omega; the relevant one tracks the beat)
    m3 = []
    for gg in [0.05, 0.1]:
        r = m2_locking_run(1.0e-3, gg, w2=WSTAR + olib(gg))
        m3.append(r)
        lam = r["lambda_lock"]
        print("M3 g_chi=%g gamma=1e-3: lambda=%s lambda/(g^2 omega*)=%s"
              % (gg, ("%.3e" % lam) if lam else "none",
                 ("%.4f" % (lam / (gg * gg * WSTAR))) if lam else "none"))
    out["M3_g_scan"] = [{**{k: r[k] for k in ("gamma", "g_chi", "lambda_lock")},
                         "lambda_over_g2_omega_star":
                             (r["lambda_lock"] / (r["g_chi"] ** 2 * WSTAR)
                              if r["lambda_lock"] else None)} for r in m3]
    out["r2_reference"] = {"lambda_toy": "0.021 +- 0.004 g_chi^2 omega* (thermalizing class)",
                           "integrable": "marginal (corpus print, IV.H.3)"}

    # merge part (a) from the Rust field-level run
    pa = os.path.join(HERE, "h27a_field_exponents.json")
    if os.path.exists(pa):
        with open(pa) as fh:
            out["parta_field_exponents"] = json.load(fh)
    else:
        out["parta_field_exponents"] = "h27a_field_exponents.json not present (run h26_convention first)"

    with open(os.path.join(HERE, "h27_results.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote h27_results.json")


if __name__ == "__main__":
    main()
