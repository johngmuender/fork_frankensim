#!/usr/bin/env python3
"""
WORKSTREAM Q2 -- analytic reconnaissance of the terminal T4-W5 proof
(F-T7-Q2): numerical companion to ANALYTIC_RECON.md.

THE PROVEN LEMMA (derivation in ANALYTIC_RECON.md; this script computes its
constants and verifies it):

For the continuum free evolution of the compactly supported L5' datum
    f0(y) = A0(y - x0) e^{i k0 y},  A0(u) = e^{-u^2/2} W(|u|)  (C2, even,
    supp A0 = [-3.5, 3.5]),  x0 = -5, k0 = 2,
the Fresnel kernel factorizes EXACTLY (no asymptotics) after recentering
at the source:  with  xi = x - x0,  vtil = xi/t,  w = vtil - k0,

    f(x,t)   = (2 pi i t)^{-1/2} e^{i(k0 x0 + xi^2/2t)}  Phi_t(w),
    Phi_t(w) = INT A0(u) e^{-i w u} e^{i u^2/2t} du,

    f'/f = i vtil + Phi_t'(w) / (t Phi_t(w)),
    =>  a := Im(f'/f) = (x - x0)/t + Im[Phi_t'/Phi_t](w)/t,
        b := Re(f'/f) =              Re[Phi_t'/Phi_t](w)/t.

Since |e^{i u^2/2t} - 1| <= u^2/2t on supp A0:
    |Phi_t(w)  - Ahat0(w) | <= m2/(2t),    m2 = INT u^2  A0 du,
    |Phi_t'(w) - Ahat0'(w)| <= m3/(2t),    m3 = INT |u|^3 A0 du,
with Ahat0 (real, even) the Fourier transform of A0 and Ahat0' its
derivative (real, odd).  Whenever  Ahat0(w) > m2/(2t)  (the admissibility
condition), writing beta(w) = Ahat0'(w)/Ahat0(w):

    |Phi_t'/Phi_t - beta(w)| <= Dbeta(w,t)
        := [ m3/(2t) + |beta(w)| m2/(2t) ] / [ Ahat0(w) - m2/(2t) ],

hence the three rigorous bounds (all constants computable):

    (B1)  | a(x,t) - (x - x0)/t |  <=  Dbeta(w,t)/t   =: C_lem(x,t),
    (B2)  | b(x,t) - beta(w)/t  |  <=  C_lem(x,t),
    (B3)  | a(x,t) - x/t |        <=  (|x0| + Dbeta)/t   (spec's x/t form).

C_lem = O(1/t^2) at fixed w.  The lemma controls the AXIAL part a and the
spin-drift generator b (v_z = n_y b entirely); the transverse spin term
n_y S(z) = -n_y pi cot(pi z) in v_x is known exactly and is NOT part of
the error term.

WHAT THIS SCRIPT DOES
  1. Computes m0, m2, m3, Ahat0, Ahat0' by piecewise Gauss-Legendre
     quadrature (A0 analytic on each piece), with node-doubling
     convergence checks and an mpmath (dps 30) cross-check of Phi_t.
  2. Verifies the decomposition identity against direct Fresnel
     quadrature of f, f' at 8 points (pure identity check, 1e-12 class).
  3. Verifies (B1)-(B3) at >= 100 (x,t) points BOTH for the continuum
     field (quadrature) AND for the ENGINE's exact band-limited field
     (P2's evaluator pattern: 2048-mode direct sum of the Nx = 2048
     FFT spectrum on [-15, 45]); ZERO violations required.  The
     engine-vs-continuum delta (band-limitation + periodization layer)
     is measured and reported at every point.
  4. Kill-window detector rows: rigorous continuum bounds for a(d, t),
     b(d, t) on t in [6.4, 7.2] at d = 1, compared with N3's measured
     line field A(t).
  5. Late-time layer diagnostic: a_engine vs a_continuum at (1, 12),
     (1, 14), (1, 15.9) -- where the periodized field departs from the
     continuum (the wrap), quantifying why the certified object is the
     discretized field.
  6. T_ball numbers for the memo's part (iii): the rigorous
     no-slow-lane time T_ball(n_y; band, T0) = T0 exp(Dz/(|n_y|
     beta_eff)) with beta_eff = min over the lag band of
     (beta - Dbeta)|_{t=T0}, scanned over T0; plus the (flagged)
     heuristic T_max consistency estimate against the measured tau_max.

hbar = m = 1.  Output: q2_results.json.
"""
import json
import os
import time

import numpy as np
from scipy import fft as sfft

# ----------------------------------------------------------------------
# Parameters (inherited verbatim from L5prime / L7b / O1 / P2)
# ----------------------------------------------------------------------
L = 1.0
SIGMA_X = 1.0
X0 = -5.0
K0 = 2.0
W_IN, W_OUT = 2.5, 3.5
D_NEAR = 1.0
XMIN, XMAX = -15.0, 45.0
NX = 2048
RHO_EPS_REL = 1e-14
KILL_BINS = [(6.4, 6.8), (6.8, 7.2)]

OUTDIR = os.path.dirname(os.path.abspath(__file__))
N3DIR = os.path.join(os.path.dirname(OUTDIR), "N3")


# ----------------------------------------------------------------------
# Envelope A0 and piecewise Gauss-Legendre machinery
# ----------------------------------------------------------------------
def smoothstep(t):
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))


def A0(u):
    u = np.abs(np.asarray(u, float))
    g = np.exp(-0.5 * u * u)
    W = np.ones_like(u)
    W[u >= W_OUT] = 0.0
    m = (u > W_IN) & (u < W_OUT)
    W[m] = smoothstep((W_OUT - u[m]) / (W_OUT - W_IN))
    return g * W


PIECES = [(-3.5, -2.5), (-2.5, 0.0), (0.0, 2.5), (2.5, 3.5)]


def gl_nodes(n_per):
    x, w = np.polynomial.legendre.leggauss(n_per)
    xs, ws = [], []
    for a, b in PIECES:
        xs.append(0.5 * (b - a) * x + 0.5 * (a + b))
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(xs), np.concatenate(ws)


NGL = 200
UN1, WN1 = gl_nodes(NGL)
A01 = A0(UN1)
UN2, WN2 = gl_nodes(2 * NGL)
A02 = A0(UN2)


def moment(p, lvl=1):
    u, w, a0 = (UN1, WN1, A01) if lvl == 1 else (UN2, WN2, A02)
    return float(np.sum(w * np.abs(u) ** p * a0))


def ahat0(wv, lvl=1):
    """Fourier transform of A0 (complex return; imag ~ 0 checked)."""
    u, w, a0 = (UN1, WN1, A01) if lvl == 1 else (UN2, WN2, A02)
    wv = np.atleast_1d(np.asarray(wv, float))
    ph = np.exp(-1j * wv[:, None] * u[None, :])
    return (ph * (w * a0)[None, :]).sum(axis=1)


def ahat0p(wv, lvl=1):
    u, w, a0 = (UN1, WN1, A01) if lvl == 1 else (UN2, WN2, A02)
    wv = np.atleast_1d(np.asarray(wv, float))
    ph = np.exp(-1j * wv[:, None] * u[None, :])
    return (ph * (w * a0 * (-1j * u))[None, :]).sum(axis=1)


def phi_t(wv, t, lvl=1):
    """Phi_t(w), Phi_t'(w) by piecewise GL quadrature."""
    u, w, a0 = (UN1, WN1, A01) if lvl == 1 else (UN2, WN2, A02)
    wv = np.atleast_1d(np.asarray(wv, float))
    ph = np.exp(1j * (-wv[:, None] * u[None, :]
                      + (u * u)[None, :] / (2.0 * t)))
    P = (ph * (w * a0)[None, :]).sum(axis=1)
    Pp = (ph * (w * a0 * (-1j * u))[None, :]).sum(axis=1)
    return P, Pp


# ----------------------------------------------------------------------
# Engine field: exact band-limited mode sum (P2 evaluator pattern)
# ----------------------------------------------------------------------
def build_spectrum():
    dx = (XMAX - XMIN) / NX
    xg = XMIN + dx * np.arange(NX)
    f0 = (A0(xg - X0)).astype(complex) * np.exp(1j * K0 * xg)
    f0 /= np.sqrt(np.sum(np.abs(f0) ** 2) * dx)
    kx = 2.0 * np.pi * sfft.fftfreq(NX, dx)
    F0 = sfft.fft(f0)
    return xg, dx, kx, F0


XG, DX, KX, F0 = build_spectrum()
CFT = F0 / NX
IKC = 1j * KX * CFT


def engine_ab(x, t):
    """Engine a, b at arbitrary (x, t): exact 2048-mode direct sum.
    No clamps applied (well-conditioned points; rho_rel monitored)."""
    ph = np.exp(1j * (KX * (x - XMIN) - 0.5 * KX ** 2 * t))
    f = (CFT * ph).sum()
    fp = (IKC * ph).sum()
    r = f.real ** 2 + f.imag ** 2
    q = np.conj(f) * fp
    return q.imag / r, q.real / r, r


def grid_rho_max(t):
    Ft = F0 * np.exp(-0.5j * KX ** 2 * t)
    return float((np.abs(sfft.ifft(Ft)) ** 2).max())


# ----------------------------------------------------------------------
# Lemma constants + bound evaluation at a point
# ----------------------------------------------------------------------
M0 = moment(0)
M2 = moment(2)
M3 = moment(3)
MOM_CONV = dict(m0=abs(M0 - moment(0, 2)), m2=abs(M2 - moment(2, 2)),
                m3=abs(M3 - moment(3, 2)))


def lemma_at(x, t):
    """Returns dict of lemma quantities at (x, t); admissibility flag."""
    xi = x - X0
    vt = xi / t
    w = vt - K0
    Ah = float(ahat0(w).real[0])
    Ahp = float(ahat0p(w).real[0])
    adm = Ah > M2 / (2.0 * t)
    out = dict(x=float(x), t=float(t), w=float(w), vtil=float(vt),
               Ahat=Ah, Ahatp=Ahp, admissible=bool(adm))
    if not adm:
        return out
    beta = Ahp / Ah
    Dbeta = (M3 / (2 * t) + abs(beta) * M2 / (2 * t)) / (Ah - M2 / (2 * t))
    out.update(beta=float(beta), Dbeta=float(Dbeta),
               C_lem=float(Dbeta / t))
    return out


def continuum_ab(x, t, lvl=1):
    xi = x - X0
    w = xi / t - K0
    P, Pp = phi_t(w, t, lvl)
    ratio = (Pp / P)[0]
    return xi / t + ratio.imag / t, ratio.real / t


# ----------------------------------------------------------------------
# Phase 1: constants, Ahat0 window, zeros
# ----------------------------------------------------------------------
def phase_constants():
    # imag-part sanity of Ahat0 (must vanish: A0 real even)
    wscan = np.linspace(-8.0, 8.0, 8001)
    Av = ahat0(wscan)
    imag_max = float(np.abs(Av.imag).max())
    Ar = Av.real
    # first zeros of Ahat0 on w >= 0 (even function)
    pos = wscan >= 0
    wp, Ap = wscan[pos], Ar[pos]
    sgn = np.sign(Ap)
    flips = np.where(sgn[1:] * sgn[:-1] < 0)[0]
    zeros = [float(0.5 * (wp[i] + wp[i + 1])) for i in flips[:3]]
    # quadrature convergence of Ahat0 on the verification band
    band = np.linspace(-1.8, 1.6, 69)
    d_ah = float(np.abs(ahat0(band) - ahat0(band, 2)).max())
    d_ahp = float(np.abs(ahat0p(band) - ahat0p(band, 2)).max())
    Ab = ahat0(band).real
    return dict(
        m0=M0, m2=M2, m3=M3, moment_convergence=MOM_CONV,
        ahat0_imag_max=imag_max,
        ahat0_first_zeros_pos_w=zeros,
        ahat0_min_on_band=dict(band=[-1.8, 1.6], value=float(Ab.min()),
                               at_w=float(band[int(np.argmin(Ab))])),
        ahat0_quadrature_conv=dict(ahat=d_ah, ahatp=d_ahp),
        gaussian_reference=dict(
            note="pure (untapered) Gaussian would give Ahat0 = "
                 "sqrt(2 pi) e^{-w^2/2}, m2/m0 = 1; taper corrections "
                 "are the differences below",
            m0_gauss=float(np.sqrt(2 * np.pi)),
            ahat_at_0_over_gauss=float(ahat0(0.0).real[0]
                                       / np.sqrt(2 * np.pi))))


# ----------------------------------------------------------------------
# Phase 2: identity check (decomposition vs direct Fresnel quadrature)
# ----------------------------------------------------------------------
def fresnel_direct(x, t):
    xi = x - X0
    ker = np.exp(1j * ((xi - UN2) ** 2 / (2.0 * t) + K0 * (X0 + UN2)))
    pref = 1.0 / np.sqrt(2.0 * np.pi * 1j * t)
    f = pref * np.sum(WN2 * A02 * ker)
    fp = pref * np.sum(WN2 * A02 * ker * 1j * (xi - UN2) / t)
    q = np.conj(f) * fp
    r = f.real ** 2 + f.imag ** 2
    return q.imag / r, q.real / r


def phase_identity():
    pts = [(1.0, 6.8), (1.0, 6.4), (0.0, 5.0), (5.0, 5.0), (9.0, 7.0),
           (15.0, 8.0), (-2.0, 6.0), (20.0, 9.5)]
    rows, worst = [], 0.0
    for x, t in pts:
        ad, bd = fresnel_direct(x, t)
        ac, bc = continuum_ab(x, t, lvl=2)
        d = max(abs(ad - ac), abs(bd - bc))
        worst = max(worst, d)
        rows.append(dict(x=x, t=t, a_direct=float(ad), a_decomp=float(ac),
                         b_direct=float(bd), b_decomp=float(bc),
                         max_diff=float(d)))
    return dict(points=rows, max_abs_diff=float(worst))


def phase_mpmath():
    import mpmath as mp
    mp.mp.dps = 30

    def a0_mp(u):
        uu = abs(u)
        if uu >= W_OUT:
            return mp.mpf(0)
        g = mp.e ** (-uu * uu / 2)
        if uu > W_IN:
            s = (W_OUT - uu) / (W_OUT - W_IN)
            g *= s ** 3 * (10 - 15 * s + 6 * s * s)
        return g

    rows, worst = [], 0.0
    for w, t in [(-1.1176470588, 6.8), (-0.5, 5.0), (0.0, 7.0),
                 (1.2, 9.0)]:
        f = mp.quad(lambda u: a0_mp(u)
                    * mp.e ** (mp.mpc(0, 1) * (-w * u + u * u / (2 * t))),
                    [-3.5, -2.5, 0, 2.5, 3.5])
        P, _ = phi_t(w, t, lvl=2)
        d = abs(complex(f) - complex(P[0]))
        worst = max(worst, d)
        rows.append(dict(w=w, t=t, gl=str(complex(P[0])),
                         mpmath=str(complex(f))[:48], abs_diff=float(d)))
    return dict(points=rows, max_abs_diff=float(worst))


# ----------------------------------------------------------------------
# Phase 3: >= 100-point verification (continuum AND engine field)
# ----------------------------------------------------------------------
def phase_verify():
    tgrid = [4.8, 5.4, 6.0, 6.4, 6.8, 7.2, 8.0, 9.0, 10.0]
    wgrid = np.linspace(-1.7, 1.5, 13)
    pts = [(X0 + (K0 + w) * t, t) for t in tgrid for w in wgrid]
    pts += [(D_NEAR, t) for t in (6.4, 6.6, 6.8, 7.0, 7.2)]

    rho_max = {t: grid_rho_max(t) for t in sorted(set(t for _, t in pts))}

    rows = []
    nviol = dict(a_cont=0, b_cont=0, a_eng=0, b_eng=0, xt_eng=0)
    hr_a_eng, hr_b_eng, hr_a_c, hr_b_c = [], [], [], []
    d_a, d_b, rho_rels = [], [], []
    for x, t in pts:
        lem = lemma_at(x, t)
        assert lem["admissible"], "inadmissible verification point (%f,%f)" \
            % (x, t)
        C = lem["C_lem"]
        vt, beta = lem["vtil"], lem["beta"]
        ac, bc = continuum_ab(x, t)
        aN, bN, r = engine_ab(x, t)
        rho_rel = r / rho_max[t]
        e_ac = abs(ac - vt)
        e_bc = abs(bc - beta / t)
        e_aN = abs(aN - vt)
        e_bN = abs(bN - beta / t)
        e_xt = abs(aN - x / t)
        Cxt = (abs(X0) + lem["Dbeta"]) / t
        v = dict(a_cont=e_ac > C, b_cont=e_bc > C,
                 a_eng=e_aN > C, b_eng=e_bN > C, xt_eng=e_xt > Cxt)
        for k in nviol:
            nviol[k] += int(v[k])
        hr_a_c.append(C - e_ac)
        hr_b_c.append(C - e_bc)
        hr_a_eng.append(C - e_aN)
        hr_b_eng.append(C - e_bN)
        d_a.append(abs(aN - ac))
        d_b.append(abs(bN - bc))
        rho_rels.append(rho_rel)
        rows.append(dict(x=round(float(x), 6), t=float(t),
                         w=round(lem["w"], 6), vtil=round(vt, 6),
                         C_lem=float(C), beta=float(beta),
                         a_cont=float(ac), b_cont=float(bc),
                         a_eng=float(aN), b_eng=float(bN),
                         err_a_eng=float(e_aN), err_b_eng=float(e_bN),
                         delta_engine_a=float(abs(aN - ac)),
                         rho_rel=float(rho_rel),
                         violation=bool(any(v.values()))))
    def stats(a):
        a = np.array(a)
        return dict(min=float(a.min()), median=float(np.median(a)),
                    max=float(a.max()))
    return dict(
        n_points=len(pts),
        n_detector_points=5,
        t_range=[min(tgrid), max(tgrid)],
        w_range=[float(wgrid[0]), float(wgrid[-1])],
        violations=nviol,
        zero_violations=bool(sum(nviol.values()) == 0),
        headroom=dict(a_cont=stats(hr_a_c), b_cont=stats(hr_b_c),
                      a_eng=stats(hr_a_eng), b_eng=stats(hr_b_eng)),
        engine_vs_continuum_delta=dict(a=stats(d_a), b=stats(d_b)),
        rho_rel=stats(rho_rels),
        rows=rows)


# ----------------------------------------------------------------------
# Phase 4: kill-window detector bounds + N3 comparison
# ----------------------------------------------------------------------
def phase_killwindow():
    try:
        raw = np.load(os.path.join(N3DIR, "n3_raw.npz"))
        tg, Aref = raw["tgrid"], raw["A_ref"]
        have_n3 = True
    except Exception:
        have_n3 = False
    rows = []
    for t in (6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2):
        lem = lemma_at(D_NEAR, t)
        C = lem["C_lem"]
        lo_a, hi_a = lem["vtil"] - C, lem["vtil"] + C
        lo_b, hi_b = lem["beta"] / t - C, lem["beta"] / t + C
        row = dict(t=t, w=round(lem["w"], 6),
                   a_bounds=[float(lo_a), float(hi_a)],
                   b_bounds=[float(lo_b), float(hi_b)],
                   a_positive_certified=bool(lo_a > 0),
                   b_positive_certified=bool(lo_b > 0))
        if have_n3:
            row["A_n3_measured"] = float(np.interp(t, tg, Aref))
            row["A_n3_inside_bounds"] = bool(lo_a - 2e-3
                                             <= row["A_n3_measured"]
                                             <= hi_a + 2e-3)
        rows.append(row)
    return dict(rows=rows,
                all_a_positive=bool(all(r["a_positive_certified"]
                                        for r in rows)),
                all_b_positive=bool(all(r["b_positive_certified"]
                                       for r in rows)))


# ----------------------------------------------------------------------
# Phase 5: late-time periodization layer diagnostic
# ----------------------------------------------------------------------
def phase_late_layer():
    rows = []
    for t in (8.0, 10.0, 12.0, 14.0, 15.9):
        ac, bc = continuum_ab(D_NEAR, t, lvl=2)
        aN, bN, r = engine_ab(D_NEAR, t)
        lem = lemma_at(D_NEAR, t)
        rows.append(dict(t=t, a_cont=float(ac), a_eng=float(aN),
                         delta=float(abs(aN - ac)),
                         a_lower_bound_cont=float(lem["vtil"]
                                                  - lem["C_lem"]),
                         C_lem=float(lem["C_lem"])))
    return dict(rows=rows)


# ----------------------------------------------------------------------
# Phase 6: T_ball + heuristic T_max consistency
# ----------------------------------------------------------------------
def phase_tball():
    wband = np.linspace(-1.2, -0.6, 61)
    Ah = ahat0(wband).real
    Ahp = ahat0p(wband).real
    beta = Ahp / Ah
    zband = [0.05, 0.95]
    Dz = zband[1] - zband[0]
    best = None
    scan = []
    for T0 in np.arange(3.0, 12.01, 0.5):
        m2t, m3t = M2 / (2 * T0), M3 / (2 * T0)
        ok = Ah > m2t
        if not ok.all():
            continue
        Dbeta = (m3t + np.abs(beta) * m2t) / (Ah - m2t)
        beff = float((beta - Dbeta).min())
        if beff <= 0:
            continue
        row = dict(T0=float(T0), beta_eff=beff,
                   T_ball_ny1=float(T0 * np.exp(Dz / (1.0 * beff))),
                   T_ball_ny075=float(T0 * np.exp(Dz / (0.75 * beff))))
        scan.append(row)
        if best is None or row["T_ball_ny1"] < best["T_ball_ny1"]:
            best = row
    # heuristic (NOT rigorous): raw beta, early T0, no Dbeta
    bbar = float(np.interp(-1.1, wband, beta))
    heur = {ny: dict(T0=2.0, beta_bar=bbar,
                     T_max_est=float(2.0 * np.exp(0.9 / (ny * bbar))))
            for ny in (1.0, 0.75)}
    return dict(
        w_band=[float(wband[0]), float(wband[-1])],
        z_band=zband,
        beta_on_band=dict(min=float(beta.min()), max=float(beta.max())),
        scan=scan, best=best,
        kill_bins=KILL_BINS,
        kill_bins_inside_reach=bool(best is not None
                                    and best["T_ball_ny1"] <= 7.2),
        heuristic_consistency=dict(
            note="NON-RIGOROUS ballpark: T_max ~ T0 exp(Dz/(|n_y| "
                 "beta_bar)) with raw beta at w = -1.1, T0 = 2; "
                 "compare measured tau_max 5.1309 (|n_y|=1), "
                 "6.017 (0.75); flagged heuristic, not evidence",
            estimates={"%.2f" % ny: heur[ny] for ny in heur},
            measured=dict(ny1=5.130950, ny075=6.0174)))


# ----------------------------------------------------------------------
def main():
    t0 = time.time()
    print("== constants ==", flush=True)
    constants = phase_constants()
    print(json.dumps({k: constants[k] for k in
                      ("m0", "m2", "m3", "ahat0_first_zeros_pos_w")}),
          flush=True)
    print("== identity check ==", flush=True)
    ident = phase_identity()
    print("  max |decomp - direct| = %.3e" % ident["max_abs_diff"],
          flush=True)
    print("== mpmath cross-check ==", flush=True)
    mpc = phase_mpmath()
    print("  max |GL - mpmath| = %.3e" % mpc["max_abs_diff"], flush=True)
    print("== verification (>=100 points) ==", flush=True)
    ver = phase_verify()
    print("  n=%d  violations=%s  min headroom (engine a) = %.4f"
          % (ver["n_points"], ver["violations"],
             ver["headroom"]["a_eng"]["min"]), flush=True)
    print("  engine-vs-continuum delta a: median %.2e max %.2e"
          % (ver["engine_vs_continuum_delta"]["a"]["median"],
             ver["engine_vs_continuum_delta"]["a"]["max"]), flush=True)
    print("== kill window ==", flush=True)
    kw = phase_killwindow()
    print("  all a>0 certified: %s   all b>0 certified: %s"
          % (kw["all_a_positive"], kw["all_b_positive"]), flush=True)
    print("== late layer ==", flush=True)
    late = phase_late_layer()
    for r in late["rows"]:
        print("  t=%5.1f  a_cont=%+.4f  a_eng=%+.4f  delta=%.2e"
              % (r["t"], r["a_cont"], r["a_eng"], r["delta"]), flush=True)
    print("== T_ball ==", flush=True)
    tb = phase_tball()
    print("  best: %s" % json.dumps(tb["best"]), flush=True)
    print("  kill bins inside reach: %s" % tb["kill_bins_inside_reach"],
          flush=True)

    results = dict(
        workstream="Q2", file_id="F-T7-Q2",
        params=dict(L=L, sigma_x=SIGMA_X, x0=X0, k0=K0,
                    taper=[W_IN, W_OUT], box=[XMIN, XMAX], Nx=NX,
                    d=D_NEAR, kill_bins=KILL_BINS, ngl_per_piece=NGL),
        lemma_constants=constants,
        identity_check=ident,
        mpmath_crosscheck=mpc,
        verification=ver,
        kill_window=kw,
        late_layer=late,
        t_ball=tb,
        gates=dict(
            G1=dict(
                spec="lemma proven + verified at >= 100 points, error "
                     "within derived bound, ZERO violations",
                n_points=ver["n_points"],
                zero_violations=ver["zero_violations"],
                identity_max_diff=ident["max_abs_diff"],
                verdict="PASS" if (ver["zero_violations"]
                                   and ver["n_points"] >= 100
                                   and ident["max_abs_diff"] < 1e-10)
                        else "FAIL"),
            G2=dict(spec="conjecture + proof obligations printed "
                         "precisely", where="ANALYTIC_RECON.md section 4",
                    verdict="PASS (memo)"),
            G3=dict(spec="honest delta stated",
                    T_ball_best=tb["best"],
                    kill_bins_inside_reach=tb["kill_bins_inside_reach"],
                    where="ANALYTIC_RECON.md section 5",
                    verdict="PASS (memo)")),
        elapsed_s=float(time.time() - t0))
    with open(os.path.join(OUTDIR, "q2_results.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    print("total elapsed %.1f s -> q2_results.json" % (time.time() - t0),
          flush=True)


if __name__ == "__main__":
    main()
