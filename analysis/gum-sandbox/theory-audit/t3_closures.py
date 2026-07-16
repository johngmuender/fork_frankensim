# T3 repaired-closure computations (within-model; campaign engine conventions)
# Engine sector constants (axi_solve conventions, t->0):
#   E6 = pi^3 INT b^2, E0 = (1/4pi) INT (1-q0), I = INT 2(q1^2+q2^2)
# Unit map (compacton-exact): c_paper = 2L/sqrt(2 pi^3), kappa_paper = 2 sqrt(pi) kappa_ours,
#   kappa_ours = L/I.  Paper-unit maps derived here: e_hat = sqrt(2) E_engine / pi,
#   i_frak = I_engine / (2 sqrt(2) pi^2).
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, brentq

pi = np.pi; s2 = np.sqrt(2.0)

print("="*78)
print("SECTION 0: exact constants and closed-form cross-checks")
print("="*78)
e0h = 64/(15*pi)          # e-hat_0
i0h = 256/(105*pi)        # i-frak_0
c0  = 128*np.sqrt(42)/(105*pi)   # corpus deep-BPS constant
c1  = 2*np.sqrt(e0h*i0h)         # scaling rung
Gs  = 16*s2/9                    # G* (G3, exact)
c3  = 64*s2/(9*pi)               # saturated closure (4/pi)G*
print(f"e0_hat={e0h:.9f} i0_hat={i0h:.9f} ratio={e0h/i0h:.9f} (7/4={7/4})")
print(f"c0   = 128 sqrt42/105pi     = {c0:.9f}")
print(f"c1   = 2 sqrt(e0 i0) (rung) = {c1:.9f}   c0/c1 = {c0/c1:.9f}  sqrt(3/2)={np.sqrt(1.5):.9f}")
print(f"G*   = 16 sqrt2/9           = {Gs:.9f}")
print(f"c3   = 64 sqrt2/9pi         = {c3:.9f}   c3/c0 = {c3/c0:.9f}  35/(6 sqrt21)={35/(6*np.sqrt(21)):.9f}")
print(f"accident: 24 sqrt21 = {24*np.sqrt(21):.6f} vs 35 pi = {35*pi:.6f}  rel={24*np.sqrt(21)/(35*pi)-1:.3e}")

# G3 tuple in engine units
R3   = 2**(5/6)
E6c  = 8*s2/9; E0c = 4*s2/3; Icore = 64*s2*pi/9
halo = 32*s2*pi/9; E0h = 4*s2/9
Itot = 8*pi*Gs; Estat = 24*s2/9; L3 = np.sqrt(8*pi)*Gs
I0eng = 128*pi*(2**2.5)/105
print(f"\nengine: Icore={Icore:.6f} halo_D={halo:.6f} Itot={Itot:.6f} (=8 pi G* {8*pi*Gs:.6f})")
print(f"clock check L^2=(2/3)I E_stat: {L3**2:.6f} vs {(2/3)*Itot*Estat:.6f}")
print(f"kappa_ours=L/I={L3/Itot:.9f} vs 1/sqrt(8pi)={1/np.sqrt(8*pi):.9f}")
# paper-unit maps
def ehat(E): return s2*E/pi
def ifrak(I): return I/(2*s2*pi**2)
print(f"map check: ehat(32sqrt2/15)={ehat(32*s2/15):.9f} vs e0h={e0h:.9f}")
print(f"map check: ifrak(I0eng)={ifrak(I0eng):.9f} vs i0h={i0h:.9f}")
print(f"\noblate core: ehat={ehat(E6c+E0c):.9f} (40/9pi={40/(9*pi):.9f})  "
      f"ifrak={ifrak(Icore):.9f} (32/9pi={32/(9*pi):.9f})  ratio={(E6c+E0c)/ (Icore/(2*pi**2)) /1:.6f}")
print(f"oblate core ehat/ifrak = {ehat(E6c+E0c)/ifrak(Icore):.9f}  (5/4? {5/4})")
g_core = Icore/I0eng; g_tot = Itot/I0eng
print(f"g_core = Icore/I0 = {g_core:.9f}  (35/24={35/24:.9f})")
print(f"g_tot  = Itot /I0 = {g_tot:.9f}  (35/12={35/12:.9f})")
# universal clock-virial identity kappa_paper^2 = ehat_tot/(2 ifrak_tot)
Etot3 = 2*Gs
et3 = ehat(Etot3); it3 = ifrak(Itot)
print(f"ehat_tot={et3:.9f} ifrak_tot={it3:.9f}  (64/9pi={64/(9*pi):.9f})  kappa^2={et3/(2*it3):.9f}")
etR = s2*e0h  # family E_tot = sqrt2 e0 (g-free); ifrak_fam = i0*g*V
for g in (1.0, 1.5):
    itR = i0h*g*s2
    print(f"family g={g}: kappa^2 = ehat_tot/(2 ifrak) = {etR/(2*itR):.9f}  (7/8g = {7/(8*g):.9f})")

print()
print("="*78)
print("SECTION 1: general-channel saturated closure  c=(4/pi) sqrt(<s2>) G_w ; k=1/sqrt(2<s2>)")
print("="*78)
# derivation:  E_static = G_w + alpha I, alpha = 1/(16 pi <sin^2>_w); fixed-L min: I*=L/sqrt(2 alpha)
# clock L^2=(2/3) I E_static  =>  L = G_w/sqrt(2 alpha);  c_paper = 2L/sqrt(2 pi^3) = (4/pi) sqrt(<s2>) G_w
for name, s2w in [("ring/tilt (binding)",1.0),("s=sin^2",6/7),("uniform",2/3),("s^2=|cos|",0.5),("polar s=|cos|",2/5)]:
    print(f"  {name:22s} <s2>={s2w:.6f}  kappa={1/np.sqrt(2*s2w):.6f}  c=(4/pi)sqrt(<s2>)*G_w")
# engine validation vs axi4 printed pair (uniform channel, eps=0.05): G_L*=2.423818 -> c=2.51979
c_check = (4/pi)*np.sqrt(2/3)*2.423818
print(f"  axi4 validation (uniform, eps=.05): formula c={c_check:.6f} vs printed 2.51979  rel={c_check/2.51979-1:.2e}")
print(f"  implied uniform-channel G_w at t=0 from printed c=2.37096: G_w={2.37096*pi/(4*np.sqrt(2/3)):.6f}")
print(f"  ring channel: c=(4/pi)*G* = {(4/pi)*Gs:.9f}")

print()
print("="*78)
print("SECTION 2 (C4 test): the corpus's OWN G.5 spheroidal locked family, honestly costed, t->0")
print("="*78)
# family: F(x)=f0(R* s), s^2 = rho^2/a^2 + z^2/c^2, a=R*/sqrt(gam), c=R* gam, lam=c/a=gam^(3/2)
# direction locked Theta=theta => b = sin^2F F_r/(2 pi^2 r^2)  (only F_r enters)
# closed-form reduction (derived in memo):
#   E6(lam) = E6(1) * P(lam),  P = (1/2) INT sin th m^3 dth,  m^2 = gam sin^2 + cos^2/gam^2
#   E0(lam) = E0(1)            (exact: value-function + volume-preserving)
#   I(lam)  = I(1) * Q(lam),   Q = (3/4) INT sin^3 th m^-3 dth
E6_1 = 16*s2/15; E0_1 = 16*s2/15; I_1 = I0eng
def mhat2(th, gam): return gam*np.sin(th)**2 + np.cos(th)**2/gam**2
def P(lam):
    gam = lam**(2/3)
    return 0.5*quad(lambda th: np.sin(th)*mhat2(th,gam)**1.5, 0, pi, limit=200)[0]
def Q(lam):
    gam = lam**(2/3)
    return 0.75*quad(lambda th: np.sin(th)**3*mhat2(th,gam)**-1.5, 0, pi, limit=200)[0]
def volchk(lam):
    gam = lam**(2/3)
    return 0.5*quad(lambda th: np.sin(th)*mhat2(th,gam)**-1.5, 0, pi, limit=200)[0]
def g_paper(lam):  # corrected F-A15-7 oblate branch / prolate arctan branch
    if lam < 1:
        b = np.sqrt(1-lam**2); return 1.5*(1/b**2 - (1-b**2)*np.arctanh(b)/b**3)
    if lam == 1: return 1.0
    a = np.sqrt(lam**2-1); return 1.5*((1+a**2)*np.arctan(a)/a**3 - 1/a**2)
print("lam    volchk    Q(lam)    g_paper   P(lam)")
for lam in (0.1,0.3,0.5,0.7,0.9,1.0,1.5):
    print(f"{lam:4.2f}  {volchk(lam):.6f}  {Q(lam):.6f}  {g_paper(lam):.6f}  {P(lam):.6f}")
# closure on the family: per-lam fixed point L^2=2AC, u=sqrt(2A/B); selection: min over lam of A/C ~ P/Q
h = lambda lam: P(lam)/Q(lam)
res = minimize_scalar(h, bounds=(1e-3, 1.0), method="bounded",
                      options={"xatol":1e-10})
lam_s = res.x
A = E6_1*P(lam_s); B = E0_1; C = I_1*Q(lam_s)
c_sph = 2*np.sqrt(A*C/pi**3)
kap_sph = 2*np.sqrt(pi*B/C)
u = np.sqrt(2*A/B)
print(f"\nHONEST spheroidal-locked closure optimum:")
print(f"  lam* = {lam_s:.6f}   g* = Q = {Q(lam_s):.6f}   E6 cost factor P = {P(lam_s):.6f}")
print(f"  c_paper* = {c_sph:.6f}   kappa_paper* = {kap_sph:.6f}   V* = u = {u:.6f}")
print(f"  vs c0 = {c0:.6f}  (short by {(1-c_sph/c0)*100:.2f}%)   vs rung {c1:.6f}  vs benchmark 2.37+-0.09")
print(f"  sanity lam=1: c = {2*np.sqrt(E6_1*I_1/pi**3):.9f} (rung {c1:.9f}), kappa = {2*np.sqrt(pi*E0_1/I_1):.9f} (sqrt(7/8)={np.sqrt(7/8):.9f})")
# zero-cost reading (corpus's G.5 idealization): c(lam) = c1 sqrt(Q); sup at lam->0 = c1*sqrt(3/2) = c0
print(f"  zero-cost reading: sup c = c1*sqrt(3/2) = {c1*np.sqrt(1.5):.9f}  (= c0 exactly: {c0:.9f}) — open endpoint, not attained")
# where the honest family crosses kappa thresholds
print(f"  kappa*(honest) vs ring 1/sqrt2={1/s2:.4f}: over-spun? {kap_sph>1/s2}")
# pull of honest spheroidal tuple vs r1 benchmark
print(f"  pulls vs <r1>: c {(c_sph-2.37)/0.09:+.1f} sigma, kappa {(kap_sph-0.802)/0.018:+.1f} sigma")

print()
print("="*78)
print("SECTION 3: C2 = C3 theorem (fixed-L ring-channel algebra)")
print("="*78)
# R(I) = G + I/(16 pi) + L^2/(2I): convex in I, minimum at I* = L sqrt(8 pi) i.e. kappa_ours = 1/sqrt(8pi)
# => unconstrained fixed-L min sits EXACTLY at the constraint boundary kappa<=kappa_crit(ring):
Ltest = L3
Igrid = np.linspace(0.2,4,200)*Ltest*np.sqrt(8*pi)
Rg = Gs + Igrid/(16*pi) + Ltest**2/(2*Igrid)
imin = np.argmin(Rg)
print(f"argmin_I R: kappa_ours = {Ltest/Igrid[imin]:.6f} vs 1/sqrt(8pi) = {1/np.sqrt(8*pi):.6f}")
print(f"=> KKT: constraint kappa<=1/sqrt2 marginally active, multiplier 0; C2 selects the C3 tuple exactly.")

print()
print("="*78)
print("SECTION 4: downstream table (within-model shifts)")
print("="*78)
rows = [("c_phys", c0, c1, c3),
        ("Lambda*sqrtJ factor (hbar/c: x c0/c)", 1.0, c0/c1, c0/c3),
        ("kappa_phys", np.sqrt(7/12), np.sqrt(7/8), 1/s2),
        ("binding depth 1-kappa", 1-np.sqrt(7/12), 1-np.sqrt(7/8), 1-1/s2),
        ("omega_th/omega0 = 2 kappa", np.sqrt(7/3), np.sqrt(7/2), s2),
        ("S4' kappa^2 g (core)", 7/8, 7/8, 0.5*35/24),
        ("S4' kappa^2 g (tot)",  7/8, 7/8, 0.5*35/12),
        ("mass-law factor E_tot/(Lam mt)", s2*e0h, s2*e0h, ehat(2*Gs)),
        ]
print(f"{'quantity':38s} {'corpus':>10s} {'C1 rung':>10s} {'C2=C3 sat':>10s}")
for nm, a, b, c in rows:
    print(f"{nm:38s} {a:10.6f} {b:10.6f} {c:10.6f}")
print(f"mass-law shift under C3: {ehat(2*Gs)/(s2*e0h):.9f}  (5/(3 sqrt2) = {5/(3*s2):.9f})")
print(f"Mc^2 under C3 = (64/9pi) Lam m~: {64/(9*pi):.6f} = ehat_tot")
# eps-law / entrainment ceiling
dc_C3 = 3.38547/3.201125 - 1   # axi4 ring-saturated eps=.05 vs exact t->0
print(f"\neps-law: corpus  c(eps)=c0[1-0.42 eps^(2/3)] (c falls). C3 measured: dc/c(+.05)={dc_C3:+.4f} (c RISES)")
print(f"   linear coeff ~ {dc_C3/0.05:.3f}; C1-family (Step2): deficit -1.03 eps^1.00 (c rises)")
for nm,k in (("C1",1.03),("C3",dc_C3/0.05)):
    print(f"   B-U1' repaired ceiling ({nm}, |dc/c|<=3e-3, exp 1): eps_e <= {3e-3/k:.2e}  (corpus: 6e-4 w/ exp 3/2)")

print()
print("="*78)
print("SECTION 5: C3 solution geometry (V-invariant fate; support volume)")
print("="*78)
gfun = lambda x: np.arcsin(x) - x*np.sqrt(1-x**2)
integ = quad(lambda th: gfun(np.sin(th))/np.sin(th)**2, 1e-9, pi-1e-9, limit=400)[0]
Vol_obl = (2*pi/3)*6*s2*integ
Vol_sph = (4*pi/3)*2**2.5
print(f"oblate compacton support volume = {Vol_obl:.6f}; sphere = {Vol_sph:.6f}; ratio = {Vol_obl/Vol_sph:.6f}")
print(f"equatorial radius (3 sqrt2 pi)^(1/3) = {(3*s2*pi)**(1/3):.6f}; axis radius 2^(5/6) = {2**(5/6):.6f}")
print(f"axis ratio lam_eff = {2**(5/6)/(3*s2*pi)**(1/3):.6f}")
print(f"(V=sqrt2 in-family identity does NOT carry over: C3 is d-stationary at d=1)")

print()
print("="*78)
print("SECTION 6: C4 coincidence registry (each tested)")
print("="*78)
print(f"(i)  G* vs c0: {Gs:.7f} vs {c0:.7f}  rel={(c0-Gs)/c0:.3e}  -> disproven exact (G3)")
c_u05 = 2.51979
print(f"(ii) uniform-channel saturated ref eps=.05: {c_u05} vs c0 {c0:.6f} rel={(c_u05-c0)/c0:+.3e} (eps-specific, channel ill-posed)")
print(f"(iii) zero-cost spheroidal sup = c0 EXACT but open-endpoint + zero-cost premise refuted (F-R4 lemma)")
print(f"(iv) kappa side: s=sin^2 locked threshold = sqrt(7/12) = corpus kappa0 EXACTLY (threshold, not solution)")
print(f"(v)  shape-exact algebra on oblate geometry: c(g_core=35/24)=2 sqrt(g e0 i0)={2*np.sqrt((35/24)*e0h*i0h):.6f} != c0, != c3")
print(f"     kappa(g=35/24)=sqrt(7/(8g))={np.sqrt(7/(8*35/24)):.6f} = sqrt(3/5) != 1/sqrt2; e/i ratio 5/4 breaks the '7'")
print(f"(vi) R_eq = {(3*s2*pi)**(1/3):.5f} vs benchmark c=2.37(9): dimensional mismatch, numerology only")
