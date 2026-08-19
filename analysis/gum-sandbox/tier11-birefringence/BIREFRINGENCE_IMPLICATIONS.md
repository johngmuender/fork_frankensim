# Digest + Implications Map: "Cosmic birefringence from a joint analysis of ACT and Planck"
## Johannes R. Eskilt (Oslo), arXiv:2608.06480v1 [astro-ph.CO], 6 Aug 2026 — read in full (8 pp)

**Phase:** Tier 11 reconnaissance (no /goal; implications map only — no
workstreams executed, no adjudication issued).  **Source:** user-uploaded
PDF (arXiv egress-blocked).  **Register (binding):** the measurement is a
nature-facing experimental result; every GUM-side statement below is
within-model bookkeeping about what the corpus's own machinery says, and
issues no verdict.

---

## 1. What the paper measures

**Claim:** joint analysis of ACT DR6 and Planck PR4 (NPIPE) polarization
gives an isotropic cosmic-birefringence angle
**β = 0.277° ± 0.057°, excluding β = 0 at 4.8σ** — the strongest
birefringence significance to date.

**Method.**  A non-zero β mixes CMB E/B modes, producing
C_ℓ^{EB,o} = tan(4β)/2·(C_ℓ^{EE,o} − C_ℓ^{BB,o}) + C_ℓ^{EB,CMB}/cos(4β)
(intrinsic CMB EB set to zero).  The killer systematic is the
miscalibration angle α (β → α + β, degenerate per instrument).  The
degeneracy is broken **independently on the two legs**:
- *Planck:* the Minami–Komatsu foreground method — Galactic emission is
  rotated by α only (it is local), the CMB by α + β; a filamentary dust
  EB model (C_ℓ^{EB,dust} = A_ℓ C_ℓ^{EE,dust} sin 4ψ_ℓ, ψ_ℓ from 353 GHz
  TB/TE; three A_ℓ amplitudes, flat positive priors; LFI synchrotron EB
  assumed zero).
- *ACT:* optics-model instrumental priors on α_i (σ ≈ 0.09–0.11° per
  band, 90% correlated within a polarization array; the Galactic plane is
  masked and only ℓ > 600 is probed).
Jointly sampled: 17 miscalibration angles + 3 dust-EB parameters + β,
over 347 EB spectra (25 ACT×ACT, 140 ACT×Planck, 182 Planck×Planck);
HFI pre-flight α priors σ_αi = 0.31° (ρ = 0.31); PA4 f150 discarded
(failed null tests).

**Key numbers and internal consistency.**
| subset | β | note |
|---|---|---|
| joint | **0.277° ± 0.057° (4.8σ)** | PTE 9.7% |
| ACT×ACT | 0.207° ± 0.073° | PTE 0.1% (poor fit ℓ ≈ 1800–1900, ACT-known); matches Diego-Palazuelos–Komatsu 0.215° ± 0.074° |
| ACT×Planck | 0.178° ± 0.078° (2.3σ) | PTE 58.6% |
| Planck×Planck | 0.387° ± 0.094° | PTE 89.1%; 1.5σ above ACT leg |
| foreground-robustness run (30% Galactic cut, drop 30/44/353 GHz, ℓ > 600) | **0.236° ± 0.067° (3.5σ)** | "cannot be explained by foreground EB correlations" |
| no Planck α priors | 0.280° ± 0.059° | foreground still constrains the Planck leg |
| no α priors at all | 0.442° ± 0.098° | prior-sensitivity exposure |

**The author's own honesty:** "unresolved systematics in the data must be
understood before we can draw strong cosmological conclusions"; β = 0
would require the dust EB model AND the instrumental priors to fail
*simultaneously*; and — crucially — **the ACT team itself did not
interpret its own ⟨ψ⟩ = 0.20° ± 0.08° as birefringence**, reporting it
as an unexplained polarization-angle discrepancy of PA5 f090/f150.
Ten-year context (Fig. 4): Minami–Komatsu 2020 (2.4σ) → Diego-Palazuelos
2022 (0.30° ± 0.11°) → Eskilt 2022 PR4 (0.33° ± 0.10°) → Eskilt–Komatsu
2022 Planck+WMAP (0.342°, 3.6σ) → ACT DR6 2025 → this work (4.8σ).
The signal is **frequency-independent** (shown in earlier work) — the
signature of an axion-like Chern–Simons coupling, and the discriminant
that excludes Faraday rotation (∝ ν⁻²).

## 2. Standard-physics implications (field-level)

1. **If cosmological, ΛCDM is parity-violated.**  The natural carrier is
   an ultralight pseudoscalar φ coupled via (g/4)φFF̃: β = (g/2)Δφ over
   the line of sight; β ≈ 0.28° ⇒ g·Δφ ≈ 9.7 mrad.  Dark-energy-like
   (m ≲ H₀), ALP dark matter (with washout constraints), and
   early-dark-energy realizations are all live; tomography (the
   reionization-bump EB at low ℓ; LiteBIRD-class) distinguishes them.
2. **The measurement is systematics-limited, not statistics-limited.**
   Both degeneracy-breaking legs rest on assumptions (dust filament
   model; pre-flight optics priors), and the two legs sit in mild 1.5σ
   tension.  Independent confirmation (Simons Observatory with its
   calibration hardware, SPT, LiteBIRD) is the actual adjudicator.
3. **If instead systematics:** the "signal" would be the first case of
   two independent experiments' calibration chains conspiring at the
   0.2–0.4° level — itself consequential for all future CMB polarization
   science (r searches inherit the same angle calibration problem).

## 3. GUM implications (within-model bookkeeping; no verdict issued)

**3.1 This touches the armed board — a first for this arc.**  The
corpus's Watch-Mode board carries **optical-activity null rows (T-N2/
T-N3, "breathing")**, and the paper registers a kill from the photon
achirality theorems: **"any tree-level vacuum optical activity"**
(Prop. Ω-2 context: induced photon operators parity-even through O(χ²),
first helicity-odd effect at O(χ³)(ka); Theorem II.3: χ₃ quenched at
tree level).  Unlike X(2370) and the STAR junction result — which
touched no stake and no row — a confirmed cosmic-birefringence detection
is an adjudicating release *touching a board row class*, which under the
Watch constitution engages the **30-day law** (a filed response
instrument, corpus-side).  At minimum, this paper warrants a dated
watch-memo response line rather than a quarterly-memo mention.

**3.2 The kill-wording question (the sharp seam).**  Two readings of
"tree-level vacuum optical activity":
- *Narrow (mechanism-specific):* the kill targets **intrinsic
  photon-sector chirality** — rotation from the substrate's own χ₃-class
  couplings, present in laboratory vacuum, frequency-flat per unit
  length.  GUM's suppression theorems predict this is O(χ³)(ka) —
  utterly negligible — so a lab null keeps the row breathing, and cosmic
  birefringence sourced by an **evolving pseudoscalar background** is a
  different mechanism that does not fire the kill as written.
- *Broad (phenomenon-level):* any confirmed rotation of linear
  polarization across cosmological vacuum is "vacuum optical activity";
  a robust β ≠ 0 fires the kill.
The corpus does not adjudicate between readings (the kill wording was
written against the lab/tree-level context of II.3/Ω-2).  **This is an
F-K0-1-class wording seam, now load-bearing**: it should be resolved
*before* the community resolves β, not after — otherwise the registry's
"failure conditions signed in advance" discipline is compromised
exactly where it may first be needed.

**3.3 The paradox of the chiral substrate.**  GUM is built on a chiral
medium, yet its printed photon-sector results are all *suppression*
theorems (II.3 quenching; Ω-2 parity-even; the knot-matter achirality
gate χ* = 0).  If nature shows cosmological parity violation in photon
propagation, the model's engineered achirality flips from a success
(matching lab nulls) to a liability: **GUM would need a mechanism for
β ≠ 0 that its own theorems don't obviously permit — or it must predict
β = 0 and stake it.**  Conversely, GUM owns exactly the ingredient the
ALP interpretation needs: a cosmological pseudoscalar sector (the T3
phasons, f ∈ [4, 60] MeV; the blue-fog chiral condensate with a global
handedness).  What is UNPRINTED: any photon–phason/fog Chern–Simons
coupling, the phason mass/evolution on cosmological timescales, and any
birefringence budget.  Whether the printed structure can produce
β ~ 0.3° — or is bounded far below it — is a well-posed, currently
unanswered within-model computation.

**3.4 The χ-sector triangle (extends Tier-9's C-SS1).**  W_χ is the
corpus's unique parity-odd energy class, and it now carries **three**
opposite-pulling duties:
1. *Suppressed* on the photon (II.3/Ω-2 — lab nulls, and the kill);
2. *Expressed* on the tube core (the Tier-9 tube-core axion, needed if
   the strong sector is to house an X-class 0⁻⁺);
3. *Bounded or expressed cosmologically* (this paper: either GUM
   predicts β compatible with 0.277° ± 0.057° via a fog/phason
   background, or it predicts β = 0 and inherits the detection as a
   potential kill under the broad reading).
One sector, three sign-and-magnitude obligations — a consistency
triangle that is falsifiable as a package and would be the natural
sixth member of the handedness chain C-EW1 (whose members already
include weak left-handedness, w_a > 0, and the δ_CP sign).  Note the
suggestive alignment the corpus may NOT claim credit for: C-EW1 already
commits GUM to one global vacuum handedness; a measured, signed,
sky-uniform β is precisely the kind of observable such a chain would
feed — *if* a coupling exists and *if* the sign comes out right, neither
of which is printed.

**3.5 Robustness-bar bookkeeping (why nothing fires today).**  The
corpus's SWP robustness definition is "≥3σ in ≥2 independent pipelines,
prior-stable."  Current status: two independent datasets, yes; joint
4.8σ and foreground-robust 3.5σ, yes; **prior-stable, no** — β moves
0.236° → 0.442° under prior removal, the two legs sit in 1.5σ tension
with different degeneracy-breaking assumptions, the ACT instrument team
itself declines the birefringence interpretation of its own angle, and
the author prints "unresolved systematics… before we can draw strong
cosmological conclusions."  Under the corpus's own bar this is a
**watch-row event, not an adjudication**: the row moves from "breathing
null" to "contested detection pending systematics," no kill fires, no
grade moves.

**3.6 What a commissioned Phase W would execute** (posed, not run):
- **W1 — the kill-wording adjudication memo**: resolve narrow-vs-broad
  from printed evidence, and (offer-class) draft the sharpened kill
  wording with the SWP bar attached, before the community resolves β.
- **W2 — the GUM birefringence budget**: compute the intrinsic
  photon-sector β from the printed O(χ³)(ka) suppression (expected:
  absurdly below observable — GUM's intrinsic prediction is β ≈ 0), and
  formalize the fog/phason-background channel's status (expected:
  UNPRINTED coupling → [CJ-new]; the W_χ triangle as the constraint).
- **W3 — the sign-chain extension**: if any within-model β mechanism
  exists, its sign is slaved to the C-EW1 global handedness — an
  additional package-falsifiability member (offer-class).
- **W4 — watch-memo response instrument**: the 30-day-law dated response
  line for the T-N3-class row (corpus-side, offer).

## 4. Bottom line

For the field: the strongest-yet (4.8σ joint, 3.5σ foreground-robust)
hint that the Universe rotates CMB polarization — either new
parity-violating physics in the dark sector or a two-experiment
calibration conspiracy; systematics, not statistics, now decide, and
the author says so plainly.  For GUM: the first external release of
this arc that touches its armed board; it exposes a kill-wording seam
that should be resolved in advance, sets up a three-duty consistency
triangle on the corpus's unique parity-odd sector, and forces a choice
the corpus has so far never had to make — predict β = 0 and stake it,
or exhibit the unprinted coupling that could make its chiral vacuum
show up in the sky.
