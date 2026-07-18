# S2 — Chirality-interpolation family: is the saturated closure constant distinctive? [F-T8-S2]

**Phase:** Tier 8 (Phase S), ROADMAP_v10_ALTERNATIVES.md workstream S2.
**Date:** 2026-07-18.  **Code:** `s2_chiral.py` (deterministic, no RNG;
reuses the frozen `gstar_solve.py` inversion machinery `gfun`/`u_of_y`
unmodified — parametrization-reuse mandate, nothing refitted; ~250 s wall).
**Outputs:** `s2_results.json`, `s2_fig.png`, `s2_checkpoint.json`
(per-phase checkpoints), `s2_run.log`.
**Sources:** `tier2-closure/gstar_RESULTS.md` (the t → 0 saturated objective
G = ∫[π³b² + W]dV, sharp bound 16√2/9, the closed-form direction-locked
minimiser via g(x) = arcsin x − x√(1−x²), sector pushforward
dV = dΩ/(2√(πW))); `tier2-closure/axi3_solve.py` (the halo-saturation
fixed-point algebra: κ* = 1/√(8π), L = √(8π)G*, I_tot = 8πG*,
E_stat_tot = (3/2)G*).

**Epistemic frame (binding).  Within-model; nothing here bears on nature.**
This workstream alternatives ONE component — the chiral (micropolar,
direction-locked) cos²Θ term of the closure potential — inside a
one-parameter family, and asks which of the corpus's closure identities are
distinctive to the chiral point λ = 1 and which are generic family
structure.  Distinctiveness, where found, is a statement about the
internal architecture of the model family only; it confers no physical
support on the corpus (§7).

---

## 1. Question and answer

**Question.** Embed the corpus's closure potential in
W_λ = (1/8π)[(1−cos ψ)² + λ sin²ψ cos²Θ], λ ∈ [0,1] (λ = 1 GUM chiral,
λ = 0 achiral).  (a) What is the sharp bound B(λ)?  (b) Is it attained off
the endpoints?  (c) Do the closure identities E6 = E0 − I/16π and
E_stat = (3/2)G* persist across λ?  (d) Does c(λ) = (4/π)G*(λ) hit any
printed corpus constant at λ ≠ 1?

**Answer.**
1. **B(λ) is attained at every λ** — the corpus's per-θ BPS reduction
   integrates in closed form for the whole family (§2.2), so
   G*(λ) = B(λ) exactly; the BPS saturation residual is ≤ 1.0e-11 at all
   9 interior sample points (gate S2-G2 PASS).
2. **The whole curve is closed form** (found and verified here):

   ```
   G*(λ) = (√2/3) [ √λ(3−2λ)/(1−λ) + (3−4λ)·arccos(√λ)/(1−λ)^{3/2} ]
   G*(0) = π/√2,   G*(1/2) = (4+π)/3,   G*(1) = 16√2/9.
   ```

   The midpoint is π-rational: at λ = 1/2 the entire sector tuple is exact
   (E0 = 2, I = 32π/3, J = 16π − 8π²/3).
3. **The literal identities are DISTINCTIVE to λ = 1; their structure is
   GENERIC.**  E6 = E0 − I/16π fails off λ = 1 with residual
   −(1−λ)(2J−I)/16π (verified to ≤ 2.8e-16 against quadrature), but it is
   the λ = 1 face of the generic d-stationarity identity
   E6 = E0 − I_λ/16π with the interpolated inertia I_λ = λI + 2(1−λ)J.
   Likewise E_stat = (3/2)G* holds at every λ under λ-consistent
   bookkeeping (halo fixed-point algebra is λ-free) and fails under literal
   chiral-I bookkeeping.  κ_paper = 1/√2 is λ-independent.  Verdict table
   in §4.
4. **One genuine interior anchor hit:** c(λ*) = π at
   λ* = 0.817779781174653.  No other printed corpus constant is crossed at
   interior λ in its own dimension class (§5).

## 2. Derivation (frozen conventions; nothing refitted)

### 2.1 Sector embedding of the family

Using (1−c)² = 2(1−c) − s² and s²cos²Θ = s² − s²sin²Θ,

```
W_λ = (1/4π)(1−q0) − i_λ(q)/(16π),
i_λ = 2λ(q1²+q2²) + 2(1−λ)|q⃗|²,     I_λ := ∫i_λ dV = λI + 2(1−λ)J,
J := ∫ sin²ψ dV,
```

so the family interpolates the **inertia sector** from the corpus's chiral
transverse I (λ = 1) to an isotropic 2∫|q⃗|²dV (λ = 0), with
G_λ = E6 + E0 − I_λ/16π = ∫[π³b² + W_λ]dV.  This decomposition is unique
given the E0 sector and is the coherent reading of the pre-registered
family (spec-interpretation note, §7).

### 2.2 The per-θ BPS ODE integrates in closed form for every λ

Locked ansatz Θ = θ.  With half-angle c = cos(F/2) and
**k(θ) := √(1 − λcos²θ)**, the identity
(1−cosF)² + λsin²F cos²θ = 4sin²(F/2)(1 − k²c²) turns the BPS condition
π³b² = W_λ into 4c² c_r = (r²/√2)√(1−k²c²), which integrates with the
campaign's own g(x) = arcsin x − x√(1−x²):

```
g( k·cos(F/2) ) = k³ r³ / (6√2),      support edge  g(k) = k³R(θ)³/(6√2).
```

At λ = 1 (k = sinθ) this is EXACTLY the corpus law
g(cos(F/2)sinθ) = ρ³/(6√2); at λ = 0 (k ≡ 1) it is a spherical profile of
radius (3√2π)^{1/3}.  Saturation reduces to the algebraic identity
√[(1−cosF)² + λsin²Fcos²θ] = 2sin(F/2)√(1−x²) at x = k cos(F/2) — exact,
b single-signed, degree 1 (deg = 1.00000000 measured at every λ).  Hence
**attainment is generic**: G*(λ) = B(λ) for all λ ∈ [0,1].
Geometry: the equatorial radius R(π/2) = (3√2π)^{1/3} = 2.370984 is
**λ-independent**; the axis radius runs from 2.370984 (λ = 0, sphere) to
2^{5/6} = 1.781797 (λ = 1, the campaign compacton) — the spherical
compacton axis profile F = 2arccos(r/2^{5/6}) appears ONLY at λ = 1.

### 2.3 Closed-form bound curve and sectors

Inner c-integral (sympy, exact):
P(m) = √(1−m)(4m²−4m+3)/(48m²) + (2m−1)asin(√m)/(16m^{5/2}); outer
u-integral by parts (m = 1−λu²; the log terms cancel exactly) gives the
B(λ) closed form of §1.  Sectors by the same route, with a = √(1−λ) and
M_k := ∫₀¹ v^{2k}/√(1−a²v²)dv (note M₁ = g(a)/(2a³) — the campaign g
again):

```
E6 = B/2  (BPS half-split, generic),      E0 = 2√2 (M0 − M1),
J  = (16√2π/3)(M0 + M1 − 2M2),            I_λ = 16π(E0 − B/2),
I  = (I_λ − 2(1−λ)J)/λ   (λ > 0;  I(0) = 8√2π²/3).
Special points:  λ=0:  E0 = π/√2,  I = 8√2π²/3,  J = 2√2π²
                 λ=1/2: B = (4+π)/3,  E0 = 2,  I = 32π/3,  J = 16π−8π²/3
                 λ=1:  E0 = 4√2/3,  I = 64√2π/9,  J = 224√2π/45.
```

All closed forms verified against independent tanh-sinh quadrature (dps 30)
to ≤ 9e-29 at every sampled λ; the λ = 1/2 specials to ≤ 7e-40 at dps 40.

### 2.4 Saturated fixed point across the family

The axi3 halo algebra (clock L² = (2/3)I E_stat; halo growth threshold) is
λ-independent: the marginal halo is small-amplitude equatorial (Θ = π/2),
where i_λ = 2η² for EVERY λ, so κ* = 1/√(8π) (κ_paper = 1/√2) exactly,
and with the λ-consistent inertia the whole tuple follows:
E_stat_tot = (3/2)G*(λ), I_λ,tot = 8πG*(λ), L = √(8π)G*(λ),
c(λ) = (4/π)G*(λ).

## 3. Numerical verification — pre-registered gates

| gate | requirement | measured | verdict |
|---|---|---|---|
| S2-G1 | B(1) = 16√2/9 and B(0) = π/√2 to ≤ 1e-10 | 2.7e-16 and 7.3e-17 (dps 40); float 2-D/1-D routes ≤ 5.3e-15 | **PASS** |
| S2-G2 | BPS saturation residual ≤ 1e-8 at every sampled λ (≥ 9 interior) | worst interior max-pointwise residual 1.0e-11 (9 interior λ = 0.1…0.9 + endpoints); attainment holds everywhere | **PASS** |
| S2-G3 | identity persistence table with GENERIC/DISTINCTIVE/OTHER verdicts | table in §4, residuals printed per λ | **PASS** |
| S2-G4 | honest bottom line; no physical-support claim | §7 | **PASS** |

Internal cross-checks: three independent B routes (2-D target GL, 1-D
reduced GL, mpmath tanh-sinh) agree to ≤ 5.3e-15; closed forms vs
quadrature ≤ 9e-29; field-level spectral G vs B ≤ 4.5e-16; degree
= 1.00000000 at all λ; support-law construction residual ≤ 1.5e-11
(3.6e-9 at λ = 0 only, an arcsin-conditioning artifact at the x → 1 edge,
noted in §7).

## 4. Identity persistence table (S2-G3)

Residuals from mpmath sector quadrature (dps 30), all λ sampled at
0, 0.1, …, 1.0.  id1 = E6 − (E0 − I/16π) (literal corpus form, chiral I);
pred = −(1−λ)(2J−I)/16π (derived); id2lit = E_stat − (3/2)G* under literal
chiral-I halo bookkeeping; the generalized forms use I_λ = λI + 2(1−λ)J.

| λ | E6 | E0 | I | J | id1 literal | id1 gen | id2 lit | id2 gen |
|---|---|---|---|---|---|---|---|---|
| 0.0 | 1.110720735 | 2.221441469 | 37.220609 | 27.915457 | −0.370240245 | 0.0e+00 | +0.370240 | 0.0e+00 |
| 0.1 | 1.128391403 | 2.151986220 | 36.056875 | 26.581000 | −0.306266077 | −2.2e-16 | +0.306266 | +4.4e-16 |
| 0.2 | 1.144918168 | 2.103486286 | 35.244251 | 25.708774 | −0.257406023 | +2.2e-16 | +0.257406 | 0.0e+00 |
| 0.3 | 1.160638273 | 2.063915477 | 34.581236 | 25.020924 | −0.215305378 | −2.2e-16 | +0.215305 | 0.0e+00 |
| 0.4 | 1.175718654 | 2.029961350 | 34.012329 | 24.444991 | −0.177588913 | +2.2e-16 | +0.177589 | +4.4e-16 |
| 0.5 | 1.190265442 | 2.000000000 | 33.510322 | 23.946537 | −0.143067891 | 0.0e+00 | +0.143068 | −4.4e-16 |
| 0.6 | 1.204353507 | 1.973073446 | 33.059163 | 23.505726 | −0.111028790 | 0.0e+00 | +0.111029 | 0.0e+00 |
| 0.7 | 1.218039003 | 1.948555357 | 32.648358 | 23.109843 | −0.080997901 | +2.2e-16 | +0.080998 | −4.4e-16 |
| 0.8 | 1.231365829 | 1.926007892 | 32.270572 | 22.750152 | −0.052639433 | 0.0e+00 | +0.052639 | 0.0e+00 |
| 0.9 | 1.244369350 | 1.905110082 | 31.920426 | 22.420342 | −0.025704038 | 0.0e+00 | +0.025704 | 0.0e+00 |
| 1.0 | 1.257078722 | 1.885618083 | 31.593834 | 22.115684 | −2.2e-16 | −2.2e-16 | 0.0e+00 | 0.0e+00 |

The literal residual matches its derived prediction −(1−λ)(2J−I)/16π to
≤ 2.8e-16 at every λ (the identity's failure off λ = 1 is itself exactly
understood).

**Verdicts:**

| identity / structure | verdict | comment |
|---|---|---|
| E6 = E0 − I/16π (literal, chiral I) | **DISTINCTIVE** (λ = 1 only) | residual −(1−λ)(2J−I)/16π, up to −0.37 at λ = 0 |
| E6 = E0 − I_λ/16π (d-stationarity) | **GENERIC** | ≤ 4.4e-16 at all λ; the literal identity is its λ = 1 face |
| E_stat = (3/2)G* (literal chiral-I bookkeeping) | **DISTINCTIVE** (λ = 1 only) | residual +(1−λ)(2J−I)/16π |
| E_stat = (3/2)G* (λ-consistent bookkeeping) | **GENERIC** | halo fixed-point algebra is λ-free; ≤ 4.4e-16 |
| support law, literal (axis compacton 2^{5/6}; g(sinθ) = ρ³/6√2) | **DISTINCTIVE** (λ = 1 only) | R_axis(λ) runs 2.371 → 1.782; sinθ = k only at λ = 1 |
| support law, generalized g(k cos(F/2)) = k³r³/6√2 | **GENERIC** | construction residual ≤ 1.5e-11 (λ > 0) |
| equatorial radius R_eq = (3√2π)^{1/3} | **GENERIC** | λ-independent exactly |
| κ_paper = 1/√2 | **GENERIC** | halo threshold sits at Θ = π/2 where i_λ = 2η² for every λ |

Net structural reading (OTHER-type finding): **chirality enters the closure
in exactly one place — which inertia sector is discounted** (transverse I
vs isotropic 2J) — and through the VALUE of G*(λ); every identity the
corpus leans on is a generic consequence of the BPS + halo fixed-point
architecture, not of chirality itself.

## 5. Anchor-hit scan (S2 part d)

Ranges: c(λ) = (4/π)G*(λ) ∈ [2√2, 64√2/(9π)] = [2.828427, 3.201125];
G*(λ) ∈ [π/√2, 16√2/9] = [2.221441, 2.514157].  Constants scanned (every
printed constant of gstar_RESULTS.md in or near range, exact values):
c_paper = 64√2/(9π) = 3.201124680; (4/π)𝔠₀ = 3.201883…; π; the axi3 c
bracket prints 3.20146 and 2.82843; κ_paper = 1/√2; G* = 16√2/9;
𝔠₀ = 2.514753626; axi3 family values 2.5144238 / 2.5144510; axi4 locked
2.515991; axi3 bound π/√2; R_eq = (3√2π)^{1/3}; R_axis = 2^{5/6};
E_stat_tot(1), E6(1), E0(1).

| curve | anchor | λ | classification |
|---|---|---|---|
| c(λ) | c_paper = 64√2/(9π) | 1.0 | endpoint (by construction) |
| c(λ) | **π** | **0.817779781174653** | **INTERIOR HIT** |
| c(λ) | 2.82843 (axi3 bracket lo print) | 6.1e-6 | endpoint: rounded print of 2√2 = c(0) exactly |
| c(λ) | (4/π)𝔠₀ = 3.201883, bracket hi 3.20146 | — | no crossing (both above c(1) = 3.201125) |
| G*(λ) | 16√2/9 | 1.0 | endpoint |
| G*(λ) | π/√2 | 0.0 | endpoint (achiral bound, by construction) |
| G*(λ) | 𝔠₀ = 2.514754, axi3 2.5144238/2.5144510, axi4 2.515991 | — | no crossing (all ABOVE max G*(1) = 2.514157 — consistent with the bound) |
| G*(λ) | R_eq = (3√2π)^{1/3} = 2.370984 | 0.466815480353222 | numerological only: LENGTH-class constant vs energy-class curve (dimension mismatch, flagged) |

The single genuine same-class interior hit is **c(λ*) = π at
λ* = 0.8178** — i.e. a hypothetical 82%-chiral variant of the model would
print c = π exactly.  λ* solves a transcendental equation
((4/π)B(λ*) = π); it shows no simple closed form and matches no printed
corpus constant.  This is numerology *within* the family, reported as the
scan demanded; it carries no significance beyond illustrating how easily
one-parameter families manufacture constant hits — a useful calibration
for the corpus's own coincidence claims (cf. the 𝔠₀ ≈ G* accident, 2.4e-4,
already disproven in gstar).

## 6. Key numbers

| quantity | closed form | value |
|---|---|---|
| G*(1) (GUM) | 16√2/9 | 2.514157444218836 |
| G*(0) (achiral) | π/√2 | 2.221441469079183 |
| G*(1/2) | (4+π)/3 | 2.380530884529931 |
| G*(λ) curve | (√2/3)[√λ(3−2λ)/(1−λ) + (3−4λ)arccos√λ/(1−λ)^{3/2}] | verified ≤ 9e-29 vs quadrature |
| c(λ) range | [2√2, 64√2/(9π)] | [2.828427, 3.201125] |
| λ* : c(λ*) = π | transcendental root | 0.817779781174653 |
| E0(1/2), I(1/2), J(1/2) | 2, 32π/3, 16π − 8π²/3 | 2.000000, 33.510322, 23.946537 |
| J(1) (new closed form) | 224√2π/45 | 22.115684 |
| R_eq (all λ) | (3√2π)^{1/3} | 2.370984 |
| R_axis(0) → R_axis(1) | (3√2π)^{1/3} → 2^{5/6} | 2.370984 → 1.781797 |
| worst interior BPS residual | — | 1.0e-11 |
| worst literal-identity residual | −(1−0)(2J−I)/16π at λ=0 | −0.370240245 |
| κ_paper | 1/√2 (all λ) | 0.707106781 |

## 7. Honest bottom line and caveats (S2-G4)

1. **What distinctiveness means here — and what it does not.**  The
   corpus's printed closure identities in their literal form (chiral-I
   discount; the 2^{5/6} axis compacton; the g(sinθ) support law) hold
   ONLY at λ = 1: within the family, the chiral point is where the
   discounted inertia sector coincides with the transverse inertia of the
   direction-locked target.  But every one of these identities generalizes
   verbatim across the family (I → I_λ, sinθ → k), and the fixed-point
   tuple (κ = 1/√2, E_stat = (3/2)G*, I_tot = 8πG*, c = (4/π)G*) is pure
   family structure.  So the corpus's closure architecture is NOT special
   to chirality; only the constant's VALUE (16√2/9 vs π/√2 vs (4+π)/3) and
   the oblateness of the minimiser change with λ.  **This is a
   within-model statement about a family of functionals.  It confers no
   physical support on the corpus: nothing here tests nature, and a
   distinctive-or-generic verdict cannot upgrade any ledger grade.**
2. **Identity 2 required an embedding choice.**  E_stat = (3/2)G* at λ ≠ 1
   is only defined once one says which inertia enters the rotational
   closure; the λ-consistent choice (I_λ, forced by W_λ itself, §2.1) makes
   it GENERIC, the literal chiral-I choice makes it DISTINCTIVE with
   exactly the id1 residual.  Both bookkeepings are computed and printed;
   this dual reporting is the printed resolution of a genuine spec
   underdetermination (coordinator-owned interpretation, flagged here).
3. **The λ = 0 support-law residual is 3.6e-9** (vs ≤ 1.5e-11 elsewhere):
   an arcsin-conditioning artifact of the residual EVALUATION at the
   x → 1 support edge (g′ → ∞), not a defect of the minimiser — the BPS
   residual at λ = 0 is 4.5e-12 and G − B = 0 to 8e-15.
4. **The anchor scan is deliberately numerological** and is reported as
   such; the one genuine hit (c = π at λ* ≈ 0.8178) is an illustration of
   family-tuning capacity, not a finding about nature.  The
   dimension-mismatched G* = R_eq coincidence at λ ≈ 0.4668 is flagged and
   should never be quoted without that flag.
5. **Bonus finding beyond the pre-registered scope** (printed, not
   gate-bearing): the closed forms of §2.3 — G*(λ), the π-rational
   midpoint (4+π)/3, and J(1) = 224√2π/45 — extend the gstar closed-form
   tuple to the whole family.  They were derived after (and verified
   against) the pre-registered quadrature program; no gate was moved.
6.  All computations sit in the campaign's frozen conventions
   (a6 = a0 = m = 1, sector constants of axi_solve); the λ-family is an
   artificial construction of this workstream, not a corpus object.  The
   corpus text is untouched by any statement here.
