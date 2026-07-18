#!/usr/bin/env python3
"""
WORKSTREAM Q1 -- validated interval enclosures for the pre-crossing
certificates (F-T7-Q1).  Upgrades P2 (F-T7-P2) one rung: the pre-crossing
margins of selected kill-bin backward paths become MACHINE-CHECKED
INEQUALITIES via rigorous interval enclosures of the backward flow over the
exact band-limited field.

THE CERTIFIED OBJECT (identical to P2/O1): the finite trigonometric sum
    f(x, t)   = sum_k C_k e^{i kappa_k (x - xmin)} e^{-i kappa_k^2 t / 2},
    C_k = FFT(f0)/Nx (float64 values declared EXACT model data),
    kappa_k = (pi/30) m_k, m_k in [-1024, 1023],
with the analytic z-factor sin(pi z) e^{-i pi^2 t/2}; guidance velocity
    v_x = Im(f'/f) - n_y pi cot(pi z),  v_z = n_y Re(f'/f).
A path's pre-crossing certificate: starting at the detector point
(x, z) = (d, z0) at s = t0 and integrating BACKWARD, the enclosure tube's
X_x LOWER bound exceeds d + delta (delta = 1e-3) at some s_stop < t0.
By O1's theorem (backward-reachability first-crossing criterion) that kills
the arrival at (z0, t0).  The engine's regularization (clips +-500,
relative rho floor, z-clip) is NOT mirrored; instead the enclosure PROVES
it never binds on the tube (rho, |a|, |b|, z monitored with rigorous
bounds), so the enclosed flow is that of BOTH the ideal and the engine law.

RIGOROUS-CONSTRUCTION STATEMENT (printed; the trust base is in G5):
(1) Scalar real intervals: IEEE-754 float64 pairs [lo, hi]; every +,-,*,/
    is the machine op (correctly rounded to nearest per IEEE 754) followed
    by one math.nextafter step outward on each endpoint -- a rigorous
    outward rounding since the true result lies strictly within the
    1-ulp neighbours of a correctly-rounded result.
(2) Complex quantities: midpoint-radius ("disc") intervals (center complex
    float64, radius float64 upper bound).  Disc products/quotients use the
    Gargantini-Henrici formulas; center magnitudes are bounded above by
    sqrt(re^2+im^2)(1+2^-40) (sqrt correctly rounded per IEEE 754, so the
    bound is rigorous AND tight -- looseness would compound through the
    phase-power chains); all radius arithmetic is padded outward by
    relative 2^-46 (= 128 ulp, >= 20x the worst-case accumulated IEEE-754
    rounding of the few-op radius formulas) plus 5e-324 absolute.  Disc
    arithmetic is rotation-tight (no box wrapping), which is what makes
    the mode-phase products below width-faithful.
(3) Mode phases: e^{i theta_k}, theta_k = m_k b1 - m_k^2 b2, are built by
    binary exponentiation (11 + 21 doubling levels) of TWO base rotations
    W = e^{i b1}, U = e^{-i b2}; the only transcendental evaluations are
    cos/sin of b1, b2 midpoints and the z-factor cot/csc, all via
    mpmath.iv (prec 120, outward rounding assumed correct), converted
    outward to float64.  Disc powers reproduce the true phase-interval
    widths |m| w(b1) + m^2 w(b2) up to the declared pads.  The per-mode
    sums f^(n) = sum_k C_k (i kappa_k)^n e^{i theta_k} are evaluated
    vectorized in numpy float64 (IEEE 754 ops), radii carried as arrays
    with the same outward pads; the coefficient discs C_k (i kappa_k)^n
    are precomputed once in mpmath.iv prec 120.
(4) Validated integrator (backward step s0 -> s1 = s0 - h, exact identity)
        X(s1) = X0 - h v(X0, s0) + int_{s1}^{s0} (u - s1) D(X(u), u) du,
    D = d/du[v] = dt_v + (v . grad) v  (the exact total derivative), so
        X(s1)  in  X0 - h v(X0, s0) + (h^2/2) hull(D(B, [s1, s0]))
    for ANY a-priori box B containing the step's solution segment.  B is
    obtained by a verified Picard step:  Y0 + [0,h].(-v(B, I_s)) subset B
    (checked in interval arithmetic; inflated and re-verified on failure).
    D needs only f..f''' :  dt_vx = Re(G3)/2, dt_vz = -(n_y/2) Im(G3),
    G3 = f'''/f - (f'/f)(f''/f); grad v from G2 = f''/f - (f'/f)^2 and
    the analytic csc^2 term.  The X0-dependence is propagated in Lohner
    QR mean-value form: X0 - c in A q  (A float 2x2, q interval box),
        X(s1) in [c - h v(c,s0) + R2] + (I - h J(B, I_s)) A q,
    J the interval Jacobian (componentwise mean-value theorem, valid since
    Y0 subset B), re-factored each step through the orthogonal factor of
    the midpoint matrix (rigorous interval inverse via 2x2 adjugate).
    The step-enclosure box is intersected with the direct evaluation.
(5) Stop rule: certify when lower(X_x) > d + delta, delta = 1e-3.  Abort
    honestly (G4) if the tube width exceeds WIDTH_CAP, if Picard fails
    after 6 inflations, or if the f-disc reaches 0.

GATES (pre-registered, ROADMAP_v9 Phase Q / Q1): G1 containment (>= 1e4
random exact-point evaluations enclosed, 0 violations; outward-rounding
sanity on >= 100 cases); G2 the P2 worst-margin path certified (widths
<= 1e-3 at the certificate point) or honestly wrap-failed under G4;
G3 >= 5 additional paths across both kill bins and both |n_y| families
incl. >= 1 small-margin path; G4 failure accounting; G5 trust statement.

Usage: q1_enclose.py [test|g1|certify|gates|fig|all]   (checkpointed)
hbar = m = 1.  File as F-T7-Q1.
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy import fft as sfft

# ----------------------------------------------------------------------
# Engine constants (verbatim O1/P2)
# ----------------------------------------------------------------------
L = 1.0
SIGMA_X = 1.0
X0C = -5.0
K0 = 2.0
W_IN, W_OUT = 2.5, 3.5
D_NEAR = 1.0
XMIN, XMAX = -15.0, 45.0
NX = 2048
VCLAMP = 500.0
RHO_EPS_REL = 1e-14
ZCLIP = 1e-9
DELTA = 1e-3                       # certificate clearance (= P2 s_pc class)
WIDTH_CAP = 2e-2                   # honest wrap-abort threshold
OUTDIR = os.path.dirname(os.path.abspath(__file__))
P2DIR = os.path.join(os.path.dirname(OUTDIR), "P2")

# path roster: selection_index in P2/_selection.npz -> step size (exact
# binary floats so s_k = t0 - k h is a single rounding); chosen from the
# float recon amplification integrals (see RESULTS.md method)
PATHS = [
    dict(idx=134, h=2.0 ** -14, role="extra: family +1.00, bin [6.4,6.8]"),
    dict(idx=194, h=2.0 ** -14, role="extra: family +1.00, bin [6.8,7.2]"),
    dict(idx=372, h=2.0 ** -14, role="extra: family +0.75, bin [6.4,6.8]"),
    dict(idx=375, h=2.0 ** -14, role="extra: family +0.75, bin [6.8,7.2]"),
    dict(idx=212, h=2.0 ** -19, role="extra: wall-adjacent immediate-exit"),
    dict(idx=395, h=2.0 ** -16, role="extra: SMALL-MARGIN (1.305), depth 2.7"),
    dict(idx=232, h=2.0 ** -16, role="G2 stretch: worst margin 0.6207, depth 3.6"),
]

# ----------------------------------------------------------------------
# Spectrum (exact model data)
# ----------------------------------------------------------------------
def build_spectrum():
    def smoothstep(t):
        return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))
    dx = (XMAX - XMIN) / NX
    xg = XMIN + dx * np.arange(NX)
    u = np.abs(xg - X0C) / SIGMA_X
    W = np.ones_like(u)
    W[u >= W_OUT] = 0.0
    m = (u > W_IN) & (u < W_OUT)
    W[m] = smoothstep((W_OUT - u[m]) / (W_OUT - W_IN))
    f0 = np.exp(-(xg - X0C) ** 2 / (2 * SIGMA_X ** 2)) * W \
        * np.exp(1j * K0 * xg)
    f0 = f0.astype(complex)
    f0 /= np.sqrt(np.sum(np.abs(f0) ** 2) * dx)
    F0 = sfft.fft(f0)
    return F0 / NX


CFT = build_spectrum()                       # exact model data (float64)
MVALS = np.concatenate([np.arange(0, NX // 2), np.arange(-NX // 2, 0)])
M2VALS = (MVALS.astype(np.int64)) ** 2       # <= 1048576, exact in float64

# ----------------------------------------------------------------------
# (1) scalar real interval arithmetic  (tuples (lo, hi), outward 1-ulp)
# ----------------------------------------------------------------------
INF = math.inf


def _lo(a):
    return math.nextafter(a, -INF)


def _hi(a):
    return math.nextafter(a, INF)


def ipt(x):
    return (x, x)


def iadd(a, b):
    return (_lo(a[0] + b[0]), _hi(a[1] + b[1]))


def isub(a, b):
    return (_lo(a[0] - b[1]), _hi(a[1] - b[0]))


def imul(a, b):
    p = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (_lo(min(p)), _hi(max(p)))


def idiv(a, b):
    if b[0] <= 0.0 <= b[1]:
        raise ZeroDivisionError("interval division by zero-containing interval")
    p = (a[0] / b[0], a[0] / b[1], a[1] / b[0], a[1] / b[1])
    return (_lo(min(p)), _hi(max(p)))


def iscale(a, c):
    p = (a[0] * c, a[1] * c)
    return (_lo(min(p)), _hi(max(p)))


def ineg(a):
    return (-a[1], -a[0])


def ihull(a, b):
    return (min(a[0], b[0]), max(a[1], b[1]))


def iwid(a):
    return _hi(a[1] - a[0])


def imid(a):
    return 0.5 * (a[0] + a[1])


def imag_up(a):
    return max(abs(a[0]), abs(a[1]))


def isubset(a, b):
    """a subset of b."""
    return b[0] <= a[0] and a[1] <= b[1]


def i_isect(a, b):
    lo, hi = max(a[0], b[0]), min(a[1], b[1])
    if lo > hi:
        raise RuntimeError("empty intersection (enclosure logic error)")
    return (lo, hi)


# ----------------------------------------------------------------------
# (2) scalar complex disc arithmetic  ((center complex, radius >= 0))
# ----------------------------------------------------------------------
RPAD = 2.0 ** -46                  # blanket relative pad (= 128 ulp >= worst
TINY = 5e-324                      # accumulated rounding of the formulas)


def _au(z):
    """rigorous TIGHT upper bound of |z|: sqrt is IEEE-754 correctly
    rounded, so sqrt(fl(re^2+im^2)) <= |z| (1+3u); pad by 2^-40 >> 3u.
    Tightness matters: a sqrt(2)-loose bound would compound per disc
    multiplication through the 21-level phase-power chains."""
    return math.sqrt(z.real * z.real + z.imag * z.imag) \
        * (1.0 + 2.0 ** -40) + TINY


def _al(z):
    """rigorous lower bound of |z| via correctly-rounded sqrt."""
    s = math.sqrt(z.real * z.real + z.imag * z.imag)
    return max(0.0, s * (1.0 - 2.0 ** -40))


def dmul(a, b):
    c = a[0] * b[0]
    aa, ab = _au(a[0]), _au(b[0])
    r = (aa * b[1] + ab * a[1] + a[1] * b[1] + aa * ab * RPAD) \
        * (1.0 + RPAD) + TINY
    return (c, r)


def ddiv(a, b):
    bl = _al(b[0])
    if bl - b[1] <= 0.0:
        raise ZeroDivisionError("disc division: denominator contains 0")
    c = a[0] / b[0]
    r = ((_au(a[0]) * b[1] + _au(b[0]) * a[1]) / (bl * (bl - b[1]))
         + _au(c) * RPAD) * (1.0 + RPAD) + TINY
    return (c, r)


def dsub(a, b):
    return (a[0] - b[0], (a[1] + b[1] + _au(a[0] - b[0]) * RPAD)
            * (1.0 + RPAD) + TINY)


def d_re(a):
    return (_lo(a[0].real - a[1]), _hi(a[0].real + a[1]))


def d_im(a):
    return (_lo(a[0].imag - a[1]), _hi(a[0].imag + a[1]))


# ----------------------------------------------------------------------
# transcendental scalar enclosures via mpmath.iv (prec 120)
# ----------------------------------------------------------------------
from mpmath import iv as _iv                                   # noqa: E402
_iv.prec = 120


def _iv2i(x):
    """mpmath.iv interval -> outward float64 interval."""
    return (_lo(_lo(float(x.a))), _hi(_hi(float(x.b))))


PI_I = _iv2i(_iv.pi)
KB_I = _iv2i(_iv.pi / 30)                       # kappa spacing pi/30
B2C_I = _iv2i((_iv.pi / 30) ** 2 / 2)           # (pi/30)^2/2
PISQ_I = _iv2i(_iv.pi ** 2)


def iv_sincos(phi_i):
    """rigorous (sin, cos) enclosures of a float interval argument."""
    a = _iv.mpf([phi_i[0], phi_i[1]])
    return _iv2i(_iv.sin(a)), _iv2i(_iv.cos(a))


def expi_disc(phi_i):
    """disc enclosure of e^{i phi}, phi in the float interval phi_i."""
    m = imid(phi_i)
    sm, cm = iv_sincos((m, m))
    hw = _hi(max(m - phi_i[0], phi_i[1] - m))    # arc half-width (chord<=arc)
    c = complex(imid(cm), imid(sm))
    r = (hw + 0.5 * iwid(cm) + 0.5 * iwid(sm) + _au(c) * RPAD) \
        * (1.0 + RPAD) + TINY
    return (c, r)


# ----------------------------------------------------------------------
# (3) vectorized disc mode sums
# ----------------------------------------------------------------------
NBIT_M, NBIT_M2 = 11, 21
ABSM = np.abs(MVALS)
MNEG = MVALS < 0
MBITS = [(ABSM >> j) & 1 == 1 for j in range(NBIT_M)]
M2BITS = [(M2VALS >> j) & 1 == 1 for j in range(NBIT_M2)]
RPV = 2.0 ** -46


def _build_coeff_discs(nmax=3):
    """A_n = C_k (i kappa_k)^n as disc arrays, via mpmath.iv prec 120."""
    AC, AR, AABS = [], [], []
    kap_lo = np.empty(NX)
    kap_hi = np.empty(NX)
    for n in range(nmax + 1):
        cc = np.empty(NX, complex)
        rr = np.empty(NX)
        for k in range(NX):
            m = int(MVALS[k])
            kn = (_iv.pi / 30 * m) ** n
            # (i kappa)^n = i^n kappa^n
            i_n = 1j ** (n % 4)
            klo, khi = _iv2i(kn)
            kmid = 0.5 * (klo + khi)
            khw = max(kmid - klo, khi - kmid)
            base = complex(CFT[k]) * i_n          # exact float64 * unit i^n
            cc[k] = base * kmid
            rr[k] = (abs(base) * khw + _au(cc[k]) * RPAD) * (1 + RPAD) + TINY
            if n == 1:
                kap_lo[k], kap_hi[k] = klo, khi
        AC.append(cc)
        AR.append(rr)
        AABS.append(np.sqrt(cc.real ** 2 + cc.imag ** 2) * (1 + 2.0 ** -40)
                    + rr)
    return AC, AR, AABS


_T0 = time.time()
AC, AR, AABS = _build_coeff_discs()
SUMABS = [float(np.sum(a)) * (1 + RPV) for a in AABS]   # sum_k |A_n_k| upper


def _vau(c):
    """vector rigorous upper |c| via correctly-rounded sqrt."""
    return np.sqrt(c.real ** 2 + c.imag ** 2) * (1.0 + 2.0 ** -40)


def vec_pow(base, bits, n):
    """disc array base^e for exponent bit arrays; base scalar disc."""
    ec = np.ones(n, complex)
    er = np.zeros(n)
    b = base
    for j, mask in enumerate(bits):
        if mask.any():
            aa = _vau(ec[mask])
            ab = _au(b[0])
            er[mask] = (aa * b[1] + ab * er[mask] + er[mask] * b[1]
                        + aa * ab * RPV) * (1.0 + RPV) + TINY
            ec[mask] = ec[mask] * b[0]
        if j + 1 < len(bits):
            b = dmul(b, b)
    return ec, er


def fields_disc(x_i, t_i, nmax):
    """disc enclosures of f^(0..nmax) at interval (x, t): full 2048 modes."""
    b1 = imul(KB_I, isub(x_i, ipt(XMIN)))
    b2 = imul(B2C_I, t_i)
    W = expi_disc(b1)
    U = expi_disc(ineg(b2))
    ec, er = vec_pow(W, MBITS, NX)
    ec = np.where(MNEG, np.conj(ec), ec)
    uc, ur = vec_pow(U, M2BITS, NX)
    aa, ab = _vau(ec), _vau(uc)
    er2 = (aa * ur + ab * er + er * ur + aa * ab * RPV) * (1 + RPV) + TINY
    ec2 = ec * uc
    ae = _vau(ec2) + er2
    out = []
    for nn in range(nmax + 1):
        tc = AC[nn] * ec2
        tr = AABS[nn] * er2 + AR[nn] * ae + AABS[nn] * _vau(ec2) * RPV
        c = complex(np.sum(tc))
        r = (float(np.sum(tr))
             + float(np.sum(np.abs(tc.real) + np.abs(tc.imag))) * 64
             * 2.0 ** -53) * (1.0 + RPV) + TINY
        out.append((c, r))
    return out


# ----------------------------------------------------------------------
# velocity / Jacobian / D assembly (scalar intervals + discs)
# ----------------------------------------------------------------------
def vel_pack(x_i, z_i, t_i, ny, nmax):
    """returns dict with v (and J, D if nmax >= 3) as scalar intervals,
    plus monitors.  nmax = 1 (velocity only) or 3 (full step pack)."""
    F = fields_disc(x_i, t_i, nmax)
    f, f1 = F[0], F[1]
    if _al(f[0]) - f[1] <= 0.0:
        raise ZeroDivisionError("f-disc contains 0")
    g = ddiv(f1, f)
    a_i, b_i = d_im(g), d_re(g)
    if not (0.0 < z_i[0] and z_i[1] < 1.0):
        raise RuntimeError("z interval left (0,1)")
    pz = imul(PI_I, z_i)
    s_i, c_i = iv_sincos(pz)
    if s_i[0] <= 0.0:
        raise RuntimeError("sin(pi z) lower bound <= 0")
    cot = idiv(c_i, s_i)
    vx = isub(a_i, iscale(imul(PI_I, cot), ny))
    vz = iscale(b_i, ny)
    rho_x_lo = max(0.0, _al(f[0]) - f[1]) ** 2       # |f|^2 lower bound
    out = dict(vx=vx, vz=vz, a=a_i, b=b_i,
               rho_x_lo=rho_x_lo, g_up=_au(g[0]) + g[1])
    if nmax >= 3:
        q2 = ddiv(F[2], f)
        G2 = dsub(q2, dmul(g, g))
        G3 = dsub(ddiv(F[3], f), dmul(g, q2))
        Jxx = d_im(G2)
        Jzx = iscale(d_re(G2), ny)
        csc2 = idiv(PISQ_I, imul(s_i, s_i))
        Jxz = iscale(csc2, ny)
        dtvx = iscale(d_re(G3), 0.5)
        dtvz = iscale(d_im(G3), -0.5 * ny)
        Dx = iadd(dtvx, iadd(imul(vx, Jxx), imul(vz, Jxz)))
        Dz = iadd(dtvz, imul(vx, Jzx))
        out.update(Jxx=Jxx, Jxz=Jxz, Jzx=Jzx, Dx=Dx, Dz=Dz)
    return out


# ----------------------------------------------------------------------
# (4) Lohner QR mean-value step
# ----------------------------------------------------------------------
def mat_ivmul_float(S, A):
    """interval 2x2 S times float 2x2 A -> interval 2x2."""
    return [[iadd(iscale(S[i][0], A[0][j]), iscale(S[i][1], A[1][j]))
             for j in range(2)] for i in range(2)]


def mat_ivmul_iv(P, Q):
    return [[iadd(imul(P[i][0], Q[0][j]), imul(P[i][1], Q[1][j]))
             for j in range(2)] for i in range(2)]


def matvec_iv(P, v):
    return [iadd(imul(P[i][0], v[0]), imul(P[i][1], v[1])) for i in range(2)]


def matvec_float_iv(A, v):
    return [iadd(iscale(v[0], A[i][0]), iscale(v[1], A[i][1]))
            for i in range(2)]


def inv2_iv(A):
    """rigorous interval inverse of a float 2x2 (adjugate / det)."""
    det = isub(imul(ipt(A[0][0]), ipt(A[1][1])),
               imul(ipt(A[0][1]), ipt(A[1][0])))
    return [[idiv(ipt(A[1][1]), det), idiv(ipt(-A[0][1]), det)],
            [idiv(ipt(-A[1][0]), det), idiv(ipt(A[0][0]), det)]]


def qr2(M):
    """orthogonal factor (float) of a 2x2 with column pivoting."""
    n0 = math.hypot(M[0][0], M[1][0])
    n1 = math.hypot(M[0][1], M[1][1])
    j0, j1 = (0, 1) if n0 >= n1 else (1, 0)
    a, b = M[0][j0], M[1][j0]
    n = math.hypot(a, b)
    if n < 1e-300:
        return [[1.0, 0.0], [0.0, 1.0]]
    q1 = (a / n, b / n)
    q2 = (-q1[1], q1[0])
    return [[q1[0], q2[0]], [q1[1], q2[1]]]


def certify_path(cfg, sel, quiet=False, max_hours=6.0):
    """run one path's backward enclosure until certified / failed.
    Checkpoint every CK_EVERY steps to _q1_ck_<idx>.json."""
    idx = cfg["idx"]
    h = cfg["h"]
    ny = float(sel["ny"][idx])
    t0 = float(sel["t0"][idx])
    z0 = float(sel["z0"][idx])
    ckfile = os.path.join(OUTDIR, "_q1_ck_%d.json" % idx)
    CK_EVERY = 20000
    REC_EVERY = max(1, int(round(2.0 ** -8 / h)))    # ~256 samples per unit s

    # state
    c = [D_NEAR, z0]
    A = [[1.0, 0.0], [0.0, 1.0]]
    q = [ipt(0.0), ipt(0.0)]
    box = [ipt(D_NEAR), ipt(z0)]
    k = 0
    rec = []
    mon = dict(min_rho_x_lo=INF, max_g_up=0.0, z_lo=z0, z_hi=z0,
               max_wid=0.0, picard_retries=0, picard_max_B_wid=0.0)
    wall0 = time.time()
    if os.path.exists(ckfile):
        ck = json.load(open(ckfile))
        k = ck["k"]
        c = [float.fromhex(v) for v in ck["c"]]
        A = [[float.fromhex(v) for v in row] for row in ck["A"]]
        q = [tuple(float.fromhex(v) for v in pair) for pair in ck["q"]]
        box = [tuple(float.fromhex(v) for v in pair) for pair in ck["box"]]
        rec = ck["rec"]
        mon = ck["mon"]
        if not quiet:
            print("  [%d] resumed at step %d" % (idx, k), flush=True)

    status, s_stop, x_lo_stop, wid_stop = None, None, None, None
    while True:
        # times: s0 = t0 - k h, s1 = s0 - h (k h exact: h binary, k < 2^24)
        s0 = t0 - k * h
        s1 = t0 - (k + 1) * h
        s0_i = (_lo(s0), _hi(s0))
        s1_i = (_lo(s1), _hi(s1))
        is_i = ihull(s1_i, s0_i)
        if s1 <= 0.05:
            status = "FAIL_REACHED_T0"
            break
        if (time.time() - wall0) > max_hours * 3600:
            status = "FAIL_TIMEOUT"
            break

        # thin velocity at the center point
        pc = vel_pack(ipt(c[0]), ipt(c[1]), s0_i, ny, 1)

        # a-priori box via verified Picard, then fat pack over (B, I_s)
        vg = [pc["vx"], pc["vz"]]
        ok = False
        infl = 1.0
        for attempt in range(6):
            B = []
            for d_ in range(2):
                mv = ihull(ipt(0.0), iscale(vg[d_], -h))       # [0,h](-v)
                w_ = iwid(box[d_]) + iwid(mv)
                Bd = iadd(box[d_], mv)
                pad = infl * (0.55 * w_ + 4.0 * h * 1e-3 + 1e-14)
                B.append((_lo(Bd[0] - pad), _hi(Bd[1] + pad)))
            try:
                pf = vel_pack(B[0], B[1], is_i, ny, 3)
            except (ZeroDivisionError, RuntimeError) as e:
                status = "FAIL_" + str(e)[:40].replace(" ", "_")
                break
            V = [pf["vx"], pf["vz"]]
            cont = all(isubset(iadd(box[d_],
                                    ihull(ipt(0.0), iscale(V[d_], -h))),
                               B[d_]) for d_ in range(2))
            if cont:
                vg = V
                ok = True
                break
            vg = V
            infl *= 2.5
            mon["picard_retries"] += 1
        if status is not None:
            break
        if not ok:
            status = "FAIL_PICARD"
            break
        mon["picard_max_B_wid"] = max(mon["picard_max_B_wid"],
                                      iwid(B[0]), iwid(B[1]))

        # step:  u_new = c - h v(c,s0) + (h^2/2) D(B, I_s)
        h22 = 0.5 * h * h
        u_new = [
            iadd(isub(ipt(c[0]), iscale(pc["vx"], h)),
                 iscale(pf["Dx"], h22)),
            iadd(isub(ipt(c[1]), iscale(pc["vz"], h)),
                 iscale(pf["Dz"], h22)),
        ]
        # S = I - h J(B, I_s)
        S = [[isub(ipt(1.0), iscale(pf["Jxx"], h)),
              ineg(iscale(pf["Jxz"], h))],
             [ineg(iscale(pf["Jzx"], h)), ipt(1.0)]]
        M = mat_ivmul_float(S, A)
        Aq = matvec_float_iv(A, q)
        direct = [iadd(u_new[d_], matvec_iv(S, Aq)[d_]) for d_ in range(2)]
        # QR refactor
        Mm = [[imid(M[i][j]) for j in range(2)] for i in range(2)]
        A2 = qr2(Mm)
        Ainv = inv2_iv(A2)
        c2 = [imid(u_new[0]), imid(u_new[1])]
        e = [isub(u_new[d_], ipt(c2[d_])) for d_ in range(2)]
        BQ = mat_ivmul_iv(Ainv, M)
        q2v = [iadd(matvec_iv(BQ, q)[d_], matvec_iv(Ainv, e)[d_])
               for d_ in range(2)]
        lohner = [iadd(ipt(c2[d_]), matvec_float_iv(A2, q2v)[d_])
                  for d_ in range(2)]
        box2 = [i_isect(lohner[d_], direct[d_]) for d_ in range(2)]

        c, A, q, box = c2, A2, q2v, box2
        k += 1

        # monitors
        w0, w1 = iwid(box[0]), iwid(box[1])
        mon["max_wid"] = max(mon["max_wid"], w0, w1)
        mon["min_rho_x_lo"] = min(mon["min_rho_x_lo"], pf["rho_x_lo"])
        mon["max_g_up"] = max(mon["max_g_up"], pf["g_up"])
        mon["z_lo"] = min(mon["z_lo"], box[1][0])
        mon["z_hi"] = max(mon["z_hi"], box[1][1])
        if k % REC_EVERY == 0 or box[0][0] > D_NEAR + DELTA:
            rec.append([s1, box[0][0], box[0][1], w0, w1])

        # certificate check
        if box[0][0] > D_NEAR + DELTA:
            status = "CERTIFIED"
            s_stop, x_lo_stop = s1, box[0][0]
            wid_stop = max(w0, w1)
            break
        if max(w0, w1) > WIDTH_CAP:
            status = "FAIL_WIDTH_CAP"
            s_stop, x_lo_stop, wid_stop = s1, box[0][0], max(w0, w1)
            break

        if k % CK_EVERY == 0:
            _ck = dict(k=k, c=[v.hex() for v in c],
                       A=[[v.hex() for v in row] for row in A],
                       q=[[v[0].hex(), v[1].hex()] for v in q],
                       box=[[v[0].hex(), v[1].hex()] for v in box],
                       rec=rec, mon=mon)
            json.dump(_ck, open(ckfile + ".tmp", "w"))
            os.replace(ckfile + ".tmp", ckfile)
            if not quiet:
                print("  [%d] step %d  s=%.4f  x=[%.6f,%.6f] w=%.2e  %.0fs"
                      % (idx, k, s1, box[0][0], box[0][1],
                         max(w0, w1), time.time() - wall0), flush=True)

    out = dict(
        idx=idx, ny=ny, t0=t0, z0=z0, h=h, role=cfg["role"],
        status=status, n_steps=k,
        s_stop=s_stop, x_lo_stop=x_lo_stop,
        certified=bool(status == "CERTIFIED"),
        delta_achieved=(None if x_lo_stop is None
                        else x_lo_stop - D_NEAR),
        wid_stop=wid_stop, depth=(None if s_stop is None else t0 - s_stop),
        final_box=[[box[0][0], box[0][1]], [box[1][0], box[1][1]]],
        mon=mon, rec=rec, elapsed=time.time() - wall0)
    json.dump(out, open(os.path.join(OUTDIR, "_q1_path_%d.json" % idx),
                        "w"))
    if os.path.exists(ckfile):
        os.remove(ckfile)
    if not quiet:
        print("  [%d] %s  steps=%d  s_stop=%s  x_lo=%s  wid=%s  (%.0f s)"
              % (idx, status, k, s_stop, x_lo_stop, wid_stop,
                 time.time() - wall0), flush=True)
    return out


# ----------------------------------------------------------------------
# gmpy2 high-precision reference point evaluator (for G1 checks only)
# ----------------------------------------------------------------------
def ref_eval(x, t, nmax=1, bits=200):
    import gmpy2
    from gmpy2 import mpfr, mpc
    gmpy2.get_context().precision = bits
    pi = gmpy2.const_pi()
    base = pi / 30
    th0 = base * (mpfr(x) - mpfr(XMIN))
    th1 = base * base * mpfr(t) / 2
    out = [mpc(0) for _ in range(nmax + 1)]
    for kk in range(NX):
        m = int(MVALS[kk])
        ph = gmpy2.exp(mpc(mpfr(0), th0 * m - th1 * (m * m)))
        cterm = mpc(mpfr(CFT[kk].real), mpfr(CFT[kk].imag)) * ph
        ik = mpc(mpfr(0), base * m)
        f = cterm
        for n in range(nmax + 1):
            out[n] += f
            f = f * ik
    return [complex(o) for o in out]


def in_disc(val, disc, slack=1.0):
    return abs(val - disc[0]) <= disc[1] * slack


# ----------------------------------------------------------------------
# phase G1: containment + outward-rounding sanity
# ----------------------------------------------------------------------
def phase_g1():
    rng = np.random.default_rng(20260718)
    n_pts = 10000
    xs = rng.uniform(-12.0, 25.0, n_pts)
    ts = rng.uniform(0.05, 7.2, n_pts)
    zs = rng.uniform(0.01, 0.99, n_pts)
    nys = rng.choice([0.75, 1.0], n_pts)
    n_checks = viol = 0
    n_vel = 2000
    wmax_f = 0.0
    t_start = time.time()
    for i in range(n_pts):
        nmax = 1 if i < n_vel else 0
        F = fields_disc(ipt(xs[i]), ipt(ts[i]), nmax)
        ref = ref_eval(xs[i], ts[i], nmax)
        for n in range(nmax + 1):
            n_checks += 1
            if not in_disc(ref[n], F[n]):
                viol += 1
        wmax_f = max(wmax_f, F[0][1])
        if i < n_vel:
            # velocity containment: a, b intervals must contain ref a, b
            pk = vel_pack(ipt(xs[i]), ipt(zs[i]), ipt(ts[i]), nys[i], 1)
            gr = ref[1] / ref[0]
            n_checks += 2
            if not (pk["a"][0] <= gr.imag <= pk["a"][1]):
                viol += 1
            if not (pk["b"][0] <= gr.real <= pk["b"][1]):
                viol += 1
        if (i + 1) % 2000 == 0:
            print("  g1 containment %d/%d  (viol %d)  %.0f s"
                  % (i + 1, n_pts, viol, time.time() - t_start), flush=True)

    # outward-rounding sanity, >= 100 cases
    sane = 0
    n_sane = 0
    fails = []
    # (a) 60 fat-box cases: box enclosure must contain interior point evals
    for j in range(60):
        xc = rng.uniform(-8, 20)
        tc = rng.uniform(0.2, 7.0)
        dx_, dt_ = 10.0 ** rng.uniform(-6, -3), 10.0 ** rng.uniform(-7, -4)
        Fb = fields_disc((xc - dx_, xc + dx_), (tc - dt_, tc + dt_), 0)[0]
        okc = True
        for _ in range(2):
            xi = rng.uniform(xc - dx_, xc + dx_)
            ti = rng.uniform(tc - dt_, tc + dt_)
            Fp = fields_disc(ipt(xi), ipt(ti), 0)[0]
            # interior thin disc must sit inside fat disc
            if abs(Fp[0] - Fb[0]) + Fp[1] > Fb[1]:
                okc = False
        n_sane += 1
        sane += okc
        if not okc:
            fails.append(("fatbox", j))
    # (b) 40 width-positivity + ref containment + width scale checks
    for j in range(40):
        xi = rng.uniform(-8, 20)
        ti = rng.uniform(0.2, 7.0)
        Fp = fields_disc(ipt(xi), ipt(ti), 0)[0]
        ref = ref_eval(xi, ti, 0)[0]
        okc = (Fp[1] > 0.0) and (Fp[1] < 1e-9) and in_disc(ref, Fp)
        n_sane += 1
        sane += okc
        if not okc:
            fails.append(("thin", j))
    # (c) 20 mpmath.iv cross-checks on the scalar interval layer
    for j in range(20):
        a = rng.uniform(-3, 3)
        b = rng.uniform(0.1, 2)
        r1 = imul(ipt(a), ipt(b))
        ex = _iv.mpf(a) * _iv.mpf(b)
        okc = (r1[0] <= float(ex.a) and float(ex.b) <= r1[1])
        s_i, c_i = iv_sincos((a, a + 1e-9))
        okc = okc and (s_i[0] <= math.sin(a + 5e-10) <= s_i[1])
        n_sane += 1
        sane += okc
        if not okc:
            fails.append(("scalar", j))

    ok = (viol == 0) and (sane == n_sane)
    out = dict(n_points=n_pts, n_containment_checks=n_checks,
               violations=viol, n_velocity_pts=n_vel,
               max_thin_f_width=wmax_f,
               sanity_cases=n_sane, sanity_passed=sane,
               sanity_fails=fails,
               reference=("gmpy2/MPFR 200-bit full 2048-mode point sums "
                          "(error <= 1e-55, far below enclosure widths)"),
               verdict="PASS" if ok else "FAIL")
    json.dump(out, open(os.path.join(OUTDIR, "_q1_g1.json"), "w"), indent=2)
    print("G1:", out["verdict"],
          json.dumps({k: v for k, v in out.items() if k != "sanity_fails"}),
          flush=True)
    return out


# ----------------------------------------------------------------------
# phase TEST: quick validation before production
# ----------------------------------------------------------------------
def phase_test():
    print("== construction ==")
    print(__doc__[__doc__.find("RIGOROUS-CONSTRUCTION"):
                  __doc__.find("GATES")])
    sel = np.load(os.path.join(P2DIR, "_selection.npz"))
    print("coeff discs built in %.1f s" % (time.time() - _T0))

    # (i) containment of a few reference evals
    for (x, t) in [(1.0, 6.405), (0.5, 3.3), (-4.0, 1.0), (12.0, 7.0)]:
        F = fields_disc(ipt(x), ipt(t), 3)
        ref = ref_eval(x, t, 3)
        stat = [(abs(ref[n] - F[n][0]) <= F[n][1]) for n in range(4)]
        print("  x=%5.1f t=%.3f  contain f..f''' %s  widths %s"
              % (x, t, stat, ["%.1e" % F[n][1] for n in range(4)]))

    # (ii) FD validation of J and D formulas (float, sign check only)
    def vflt(x, z, t, ny):
        ph = np.exp(1j * ((np.pi / 30) * MVALS * (x - XMIN)
                          - 0.5 * ((np.pi / 30) * MVALS) ** 2 * t))
        f = np.sum(CFT * ph)
        f1 = np.sum(CFT * 1j * (np.pi / 30) * MVALS * ph)
        g = f1 / f
        return np.array([g.imag - ny * np.pi / np.tan(np.pi * z),
                         ny * g.real])
    x, z, t, ny = 0.83, 0.61, 5.2, 0.75
    pk = vel_pack(ipt(x), ipt(z), ipt(t), ny, 3)
    eps = 1e-6
    Jfd = np.zeros((2, 2))
    Jfd[:, 0] = (vflt(x + eps, z, t, ny) - vflt(x - eps, z, t, ny)) / (2 * eps)
    Jfd[:, 1] = (vflt(x, z + eps, t, ny) - vflt(x, z - eps, t, ny)) / (2 * eps)
    dtv = (vflt(x, z, t + eps, ny) - vflt(x, z, t - eps, ny)) / (2 * eps)
    v0 = vflt(x, z, t, ny)
    Dfd = dtv + Jfd @ v0
    print("  J analytic vs FD: Jxx %.6f/%.6f Jxz %.6f/%.6f Jzx %.6f/%.6f"
          % (imid(pk["Jxx"]), Jfd[0, 0], imid(pk["Jxz"]), Jfd[0, 1],
             imid(pk["Jzx"]), Jfd[1, 0]))
    print("  D analytic vs FD: Dx %.6f/%.6f  Dz %.6f/%.6f"
          % (imid(pk["Dx"]), Dfd[0], imid(pk["Dz"]), Dfd[1]))
    assert abs(imid(pk["Dx"]) - Dfd[0]) < 1e-4 * max(1, abs(Dfd[0]))
    assert abs(imid(pk["Dz"]) - Dfd[1]) < 1e-4 * max(1, abs(Dfd[1]))
    print("  J, D formula FD-validation PASS")

    # (iii) 400-step trial on path 134, tube must contain gmpy2-free float
    cfg = dict(idx=134, h=2.0 ** -14, role="test")
    t_w = time.time()
    idx = cfg["idx"]
    ny = float(sel["ny"][idx])
    t0 = float(sel["t0"][idx])
    z0 = float(sel["z0"][idx])
    # temporary shallow run: hack via reduced WIDTH_CAP? just run 400 steps
    # by monkey-limit: use certify_path but stop early via max_hours tiny?
    # simpler: inline mini-run
    c = [D_NEAR, z0]
    A = [[1.0, 0.0], [0.0, 1.0]]
    q = [ipt(0.0), ipt(0.0)]
    box = [ipt(D_NEAR), ipt(z0)]
    h = cfg["h"]
    from scipy.integrate import solve_ivp

    def rhsf(s, y):
        return vflt(y[0], y[1], s, ny)
    ref = solve_ivp(rhsf, (t0, t0 - 400 * h), [D_NEAR, z0], method="DOP853",
                    rtol=1e-12, atol=1e-14, dense_output=True)
    ok_all = True
    for k in range(400):
        s0 = t0 - k * h
        s1 = s0 - h
        s0_i, s1_i = (_lo(s0), _hi(s0)), (_lo(s1), _hi(s1))
        is_i = ihull(s1_i, s0_i)
        pc = vel_pack(ipt(c[0]), ipt(c[1]), s0_i, ny, 1)
        vg = [pc["vx"], pc["vz"]]
        for attempt in range(6):
            B = []
            for d_ in range(2):
                mv = ihull(ipt(0.0), iscale(vg[d_], -h))
                Bd = iadd(box[d_], mv)
                pad = (2.5 ** attempt) * (0.55 * (iwid(box[d_]) + iwid(mv))
                                          + 4.0 * h * 1e-3 + 1e-14)
                B.append((_lo(Bd[0] - pad), _hi(Bd[1] + pad)))
            pf = vel_pack(B[0], B[1], is_i, ny, 3)
            V = [pf["vx"], pf["vz"]]
            if all(isubset(iadd(box[d_], ihull(ipt(0.0),
                                               iscale(V[d_], -h))), B[d_])
                   for d_ in range(2)):
                break
            vg = V
        h22 = 0.5 * h * h
        u_new = [iadd(isub(ipt(c[0]), iscale(pc["vx"], h)),
                      iscale(pf["Dx"], h22)),
                 iadd(isub(ipt(c[1]), iscale(pc["vz"], h)),
                      iscale(pf["Dz"], h22))]
        S = [[isub(ipt(1.0), iscale(pf["Jxx"], h)),
              ineg(iscale(pf["Jxz"], h))],
             [ineg(iscale(pf["Jzx"], h)), ipt(1.0)]]
        M = mat_ivmul_float(S, A)
        Aq = matvec_float_iv(A, q)
        direct = [iadd(u_new[d_], matvec_iv(S, Aq)[d_]) for d_ in range(2)]
        Mm = [[imid(M[i][j]) for j in range(2)] for i in range(2)]
        A2 = qr2(Mm)
        Ainv = inv2_iv(A2)
        c2 = [imid(u_new[0]), imid(u_new[1])]
        e = [isub(u_new[d_], ipt(c2[d_])) for d_ in range(2)]
        BQ = mat_ivmul_iv(Ainv, M)
        q = [iadd(matvec_iv(BQ, q)[d_], matvec_iv(Ainv, e)[d_])
             for d_ in range(2)]
        lohner = [iadd(ipt(c2[d_]), matvec_float_iv(A2, q)[d_])
                  for d_ in range(2)]
        box = [i_isect(lohner[d_], direct[d_]) for d_ in range(2)]
        c, A = c2, A2
        if k % 100 == 99:
            yr = ref.sol(s1)
            inb = (box[0][0] <= yr[0] <= box[0][1]
                   and box[1][0] <= yr[1] <= box[1][1])
            ok_all = ok_all and inb
            print("  step %d: s=%.4f  x=[%.8f, %.8f] (w %.2e)  ref-in %s"
                  % (k + 1, s1, box[0][0], box[0][1], iwid(box[0]), inb))
    dt = (time.time() - t_w) / 400
    print("  trial: per-step %.1f ms -> est deep path (233k steps) %.0f min"
          % (dt * 1e3, dt * 233000 / 60))
    print("  float-reference containment over trial:", ok_all)
    assert ok_all
    print("TEST PASS")


# ----------------------------------------------------------------------
# phase CERTIFY (parallel across paths)
# ----------------------------------------------------------------------
def _worker(cfg):
    sel = np.load(os.path.join(P2DIR, "_selection.npz"))
    try:
        return certify_path(cfg, sel)
    except Exception as exc:                       # noqa: BLE001
        import traceback
        traceback.print_exc()
        return dict(idx=cfg["idx"], status="FAIL_EXC_" + str(exc)[:60],
                    certified=False)


def phase_certify():
    from multiprocessing import Pool
    todo = [cfg for cfg in PATHS
            if not os.path.exists(os.path.join(
                OUTDIR, "_q1_path_%d.json" % cfg["idx"]))]
    done = [cfg["idx"] for cfg in PATHS if cfg not in todo]
    if done:
        print("already done:", done, flush=True)
    if not todo:
        return
    # deep paths first so they overlap with the shallow ones
    todo.sort(key=lambda cf: -0.0 if cf["idx"] in (232, 395) else 1.0)
    with Pool(min(4, len(todo))) as pool:
        for r in pool.imap_unordered(_worker, todo):
            print("done: path %d -> %s" % (r["idx"], r["status"]),
                  flush=True)


# ----------------------------------------------------------------------
# phase GATES
# ----------------------------------------------------------------------
def phase_gates():
    g1 = json.load(open(os.path.join(OUTDIR, "_q1_g1.json")))
    sel = np.load(os.path.join(P2DIR, "_selection.npz"))
    p2 = json.load(open(os.path.join(P2DIR, "p2_results.json")))
    runs = {}
    for cfg in PATHS:
        fn = os.path.join(OUTDIR, "_q1_path_%d.json" % cfg["idx"])
        if os.path.exists(fn):
            runs[cfg["idx"]] = json.load(open(fn))

    worst = runs.get(232)
    g2_ok = (worst is not None and worst["certified"]
             and worst["wid_stop"] <= 1e-3)
    g2 = dict(path=worst, spec=("worst-margin path (P2 G4: idx 232, "
                                "margin 0.6207, s_pc 2.8499) certified "
                                "with widths <= 1e-3 at certificate"),
              p2_margin_mp=p2["G4"]["margin_mp_str"],
              verdict="PASS" if g2_ok else "FAIL")

    extras = [runs[c["idx"]] for c in PATHS if c["idx"] != 232
              and c["idx"] in runs]
    cert = [r for r in extras if r["certified"]]
    fams = set(round(abs(r["ny"]), 2) for r in cert)
    bins = set("bin1" if r["t0"] < 6.8 else "bin2" for r in cert)
    small = [r for r in cert if r["idx"] == 395]
    g3_ok = (len(cert) >= 5 and fams == {0.75, 1.0}
             and bins == {"bin1", "bin2"} and len(small) >= 1)
    g3 = dict(n_extra_certified=len(cert),
              families=sorted(fams), bins=sorted(bins),
              small_margin_certified=bool(small),
              detail=[dict(idx=r["idx"], ny=r["ny"], t0=r["t0"],
                           z0=r["z0"], role=r["role"], h=r["h"],
                           s_stop=r["s_stop"], depth=r["depth"],
                           x_lo_stop=r["x_lo_stop"],
                           delta_achieved=r["delta_achieved"],
                           wid_stop=r["wid_stop"],
                           max_wid=r["mon"]["max_wid"],
                           n_steps=r["n_steps"]) for r in extras],
              verdict="PASS" if g3_ok else
              ("PARTIAL" if len(cert) >= 3 else "FAIL"))

    failures = [dict(idx=r["idx"], status=r["status"],
                     s_reached=(r["rec"][-1][0] if r.get("rec") else None),
                     max_wid=r.get("mon", {}).get("max_wid"))
                for r in runs.values() if not r["certified"]]
    g4 = dict(n_failures=len(failures), failures=failures,
              note=("failure accounting: every wrap/blow-up reported with "
                    "widths; absence of failures also recorded here"),
              verdict="PASS")

    # engine-law coincidence: on every certified tube the engine's
    # regularization provably never binds, so the enclosed flow is that of
    # BOTH the ideal velocity field and the engine's regularized law.
    # rmax(t) <= (sum_k |C_k|)^2 (triangle inequality, outward-rounded).
    s0u = SUMABS[0]
    floor_ub = _hi(RHO_EPS_REL * s0u * s0u)
    coin = dict(
        rho_floor_upper_bound=floor_ub,
        checks=[dict(idx=r["idx"],
                     min_rho_x_lo=r["mon"]["min_rho_x_lo"],
                     rho_headroom=r["mon"]["min_rho_x_lo"] / floor_ub,
                     max_g_up=r["mon"]["max_g_up"],
                     clamp_headroom=VCLAMP / max(r["mon"]["max_g_up"],
                                                 1e-30),
                     z_range=[r["mon"]["z_lo"], r["mon"]["z_hi"]],
                     z_clip_ok=bool(r["mon"]["z_lo"] > ZCLIP
                                    and r["mon"]["z_hi"] < 1 - ZCLIP))
                for r in runs.values()],
    )
    coin["all_ok"] = all(c["rho_headroom"] > 1e3 and c["clamp_headroom"] > 1
                         and c["z_clip_ok"] for c in coin["checks"])

    g5 = dict(
        trust=("TRUST BASE: (i) IEEE-754 float64 semantics of numpy/CPython "
               "arithmetic (+,-,*,/ and sqrt correctly rounded to nearest); "
               "(ii) math.nextafter correctness; (iii) mpmath.iv outward "
               "rounding correctness for pi, sin, cos at prec 120 (used "
               "only for base rotations and the z-factor); (iv) gmpy2/MPFR "
               "correct rounding (used ONLY in G1 reference checks, not in "
               "any certificate); (v) the float64 spectral coefficients "
               "C_k and path data (t0, z0, n_y, d) are EXACT model data."),
        scope=("Within-model: the certified object is the Nx=2048 "
               "band-limited discretized field O1/P2 computed, not the "
               "continuum field; the enclosures make the pre-crossing "
               "margins of the certified paths machine-checked "
               "inequalities GIVEN the trust base. CAPD/COSY-grade work "
               "would replace (i)-(iii) by a formally verified interval "
               "kernel and extend certification to all 64,000 kill-bin "
               "paths; the analytic T4-W5 proof remains the terminal "
               "open item."),
        verdict="PASS")

    results = dict(
        workstream="Q1", file_id="F-T7-Q1",
        params=dict(Nx=NX, box=[XMIN, XMAX], d=D_NEAR, delta=DELTA,
                    width_cap=WIDTH_CAP, rpad=RPAD,
                    paths=[dict(idx=c["idx"], h=c["h"], role=c["role"])
                           for c in PATHS]),
        construction=__doc__[__doc__.find("RIGOROUS-CONSTRUCTION"):
                             __doc__.find("GATES")].strip(),
        engine_coincidence=coin,
        G1=g1, G2=g2, G3=g3, G4=g4, G5=g5,
        runs={str(k): {kk: vv for kk, vv in v.items() if kk != "rec"}
              for k, v in runs.items()},
    )
    json.dump(results, open(os.path.join(OUTDIR, "q1_results.json"), "w"),
              indent=1)
    print("verdicts:", {g: results[g]["verdict"]
                        for g in ("G1", "G2", "G3", "G4", "G5")}, flush=True)
    return results


# ----------------------------------------------------------------------
# phase FIG
# ----------------------------------------------------------------------
def phase_fig():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    INK, INK2, SURF = "#0b0b0b", "#52514e", "#fcfcfb"
    BLU, BLUD, ACC, GRN = "#2a78d6", "#12448c", "#e87ba4", "#3e8f6b"

    runs = {}
    for cfg in PATHS:
        fn = os.path.join(OUTDIR, "_q1_path_%d.json" % cfg["idx"])
        if os.path.exists(fn):
            runs[cfg["idx"]] = json.load(open(fn))
    res = json.load(open(os.path.join(OUTDIR, "q1_results.json")))

    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.4), dpi=160)
    fig.patch.set_facecolor(SURF)
    axA, axB, axC = axes

    # A: X_x tube for the two deep paths
    for idx, col, lab in ((232, BLUD, "worst margin (idx 232)"),
                          (395, ACC, "small margin (idx 395)")):
        if idx not in runs:
            continue
        r = np.array(runs[idx]["rec"])
        axA.fill_between(r[:, 0], r[:, 1], r[:, 2], color=col, alpha=0.35,
                         lw=0)
        axA.plot(r[:, 0], r[:, 1], color=col, lw=1.0, label=lab)
    axA.axhline(D_NEAR + DELTA, color=GRN, lw=1.2, ls="--")
    axA.text(3.1, D_NEAR + DELTA + 0.05, "d + $\\delta$", color=GRN,
             fontsize=8)
    axA.set_xlabel("backward epoch s")
    axA.set_ylabel("X$_x$ enclosure")
    axA.invert_xaxis()
    axA.set_title("A  validated X$_x$ tubes (deep paths)\n"
                  "backward from the detector point", fontsize=10,
                  loc="left", color=INK)
    axA.legend(fontsize=8, frameon=False)

    # B: tube widths vs s
    for cfg in PATHS:
        idx = cfg["idx"]
        if idx not in runs or not runs[idx].get("rec"):
            continue
        r = np.array(runs[idx]["rec"])
        axB.semilogy(r[:, 0], np.maximum(r[:, 3], 1e-18), lw=1.1,
                     label="%d (h=2$^{%d}$)" % (idx,
                                                round(math.log2(cfg["h"]))))
    axB.axhline(1e-3, color=INK2, lw=1.0, ls="--")
    axB.text(4.5, 1.4e-3, "G2 width gate 1e-3", fontsize=7, color=INK2)
    axB.set_xlabel("backward epoch s")
    axB.set_ylabel("X$_x$ tube width")
    axB.invert_xaxis()
    axB.set_title("B  enclosure widths along the backward flow",
                  fontsize=10, loc="left", color=INK)
    axB.legend(fontsize=6.5, frameon=False, ncol=2)

    # C: certified clearances
    labs, vals, wids, cols = [], [], [], []
    for cfg in PATHS:
        idx = cfg["idx"]
        if idx not in runs:
            continue
        r = runs[idx]
        labs.append(str(idx))
        vals.append(r["delta_achieved"] if r["certified"] else 0.0)
        wids.append(r["wid_stop"] if r["wid_stop"] else np.nan)
        cols.append(BLU if r["certified"] else ACC)
    xpos = np.arange(len(labs))
    axC.bar(xpos, vals, color=cols, width=0.62)
    for xp, w in zip(xpos, wids):
        axC.text(xp, max(vals) * 0.03, "w=%.0e" % w, rotation=90,
                 fontsize=6.5, ha="center", color=INK)
    axC.axhline(DELTA, color=GRN, lw=1.2, ls="--")
    axC.set_xticks(xpos, labs)
    axC.set_xlabel("path (selection index)")
    axC.set_ylabel("certified X$_x$ lower bound $-$ d at s$_{stop}$")
    axC.set_title("C  certified pre-crossing clearances\n"
                  "(rigorous lower bounds; w = tube width at certificate)",
                  fontsize=10, loc="left", color=INK)

    for ax in axes:
        ax.set_facecolor(SURF)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
    fig.suptitle("Q1 — validated interval enclosures for the pre-crossing "
                 "certificates (F-T7-Q1)", fontsize=11, color=INK, x=0.01,
                 ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(OUTDIR, "q1_fig.png"), facecolor=SURF,
                bbox_inches="tight")
    print("figure written", flush=True)


# ----------------------------------------------------------------------
def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if phase == "test":
        phase_test()
        return
    if phase in ("g1", "all") and not os.path.exists(
            os.path.join(OUTDIR, "_q1_g1.json")):
        print("== phase g1 ==", flush=True)
        phase_g1()
    if phase in ("certify", "all"):
        print("== phase certify ==", flush=True)
        phase_certify()
    if phase in ("gates", "all"):
        print("== phase gates ==", flush=True)
        phase_gates()
    if phase in ("fig", "all"):
        print("== phase fig ==", flush=True)
        phase_fig()
    print("total elapsed %.1f s" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
