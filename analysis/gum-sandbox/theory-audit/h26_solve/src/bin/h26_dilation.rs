//! Phase H2.6 tier-2 — the ENGINE-ANCHORED EXACT DILATION FAMILY solve of
//! both spin-anchor closures, plus the H2.7(a) exponents on the same family.
//!
//! WHY THIS TIER EXISTS: the descent-relaxed solve (h26_convention) measures
//! the guarded ANF endpoint at fixed iteration budget, and at the clock root
//! the state is over-spun (kappa above the halo threshold), i.e. F-R5
//! halo-unstable: the descent never becomes stationary and the root drifts
//! with budget.  The corpus's IV.H.3 claims (and T2's +-1/2 proof) live on
//! the CONSTRAINED STATIONARY family — envelope over shape at fixed L —
//! which a descent engine cannot sit on (saddle).  This tier realizes that
//! family EXACTLY on the engine's own relaxed solution:
//!
//! Under the spatial dilation x -> s x of any stored configuration the
//! sector functionals scale EXACTLY (2 derivatives: s^{-2} * volume s^3;
//! etc.):
//!
//!     E2 -> s E2h,  E4 -> E4h/s,  E6 -> E6h/s^3,  E0 -> s^3 E0h,
//!     I  -> s^3 Ih          (V := s^3 is the corpus's dilation dial:
//!                            E6 ~ V^{-1}, E0 ~ V, I ~ V — App G.5's family,
//!                            with the eps-suppressed E2/E4 corrections the
//!                            pure corpus family drops)
//!
//! so, anchoring (E2h, E4h, E6h, E0h, Ih) at the campaign's N = 48 guarded
//! static solution (eps = 0.05 protocol, seed 424242 — identical to
//! h26_convention's static stage), the constrained family
//!
//!     E(s; L) = t (s E2h + E4h/s) + E6h/s^3 + s^3 E0h + L^2/(2 Ih s^3)
//!
//! is CLOSED FORM.  Envelope: s*(L) = argmin_s E(s; L) (1-D bisection on
//! dE/ds, machine precision).  Clock under anchor a (derivation in
//! h26_convention.rs): L^2 = [2a/(2-a)] I(s*) E_static(s*), bisected to
//! machine precision.  Everything is budget-independent and deterministic.
//!
//! Outputs: the solved tuple per anchor (E_rot/E, V = (s*/s0)^3, kappa vs
//! thresholds, c_paper), and the H2.7(a) exponents d ln omega_mech / d ln h,
//! d ln omega_phase / d ln h at delta = 1e-3 (central differences on the
//! exact family; hbar = 2L).  T2 predicts +1/2 / -1/2 exactly on the PURE
//! family at x = 1; the measured deviations here are the eps-corrections
//! (E2/E4) plus the engine anchoring — the field-anchored empirical values.
//!
//! Usage: h26_dilation <outdir> [threads]     (~2 min: radial + static ANF)
//!
//! Epistemic notice (binding): within-model numerical engineering on a
//! speculative theory's functional; nothing about nature.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{add_scaled, bump_field};
use fs_gum_statics::{kappa_threshold, AnfParams, Opts, FLOOR_BAND};

use std::time::Instant;

const LBOX: f64 = 4.5;
const N: usize = 48;
const MAXIT_STATIC: usize = 900;
const SEED_PERT: u64 = 424242;
const PI: f64 = std::f64::consts::PI;

/// The exact dilation family anchored at (e2h, e4h, e6h, e0h, ih).
#[derive(Clone, Copy)]
struct Fam {
    t: f64,
    e2h: f64,
    e4h: f64,
    e6h: f64,
    e0h: f64,
    ih: f64,
}

impl Fam {
    fn estat(&self, s: f64) -> f64 {
        self.t * (s * self.e2h + self.e4h / s) + self.e6h / (s * s * s)
            + s * s * s * self.e0h
    }
    fn i(&self, s: f64) -> f64 {
        s * s * s * self.ih
    }
    fn e(&self, s: f64, l: f64) -> f64 {
        self.estat(s) + l * l / (2.0 * self.i(s))
    }
    fn deds(&self, s: f64, l: f64) -> f64 {
        let s2 = s * s;
        let s4 = s2 * s2;
        self.t * (self.e2h - self.e4h / s2) - 3.0 * self.e6h / s4
            + 3.0 * s2 * self.e0h
            - 3.0 * l * l / (2.0 * self.ih * s4)
    }
    /// Envelope minimizer s*(L): bisection on dE/ds (dE/ds is increasing
    /// through the minimum; bracket [0.3, 6] covers V in [0.027, 216]).
    fn sstar(&self, l: f64) -> f64 {
        let (mut lo, mut hi) = (0.3_f64, 6.0_f64);
        assert!(self.deds(lo, l) < 0.0 && self.deds(hi, l) > 0.0, "envelope bracket");
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if self.deds(mid, l) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }
    /// Clock mismatch F(L) = L^2 - c_a I(s*) E_static(s*) on the envelope.
    fn clock_f(&self, l: f64, c_a: f64) -> f64 {
        let s = self.sstar(l);
        l * l - c_a * self.i(s) * self.estat(s)
    }
    /// Bisect the clock root (Section-K bisection on the exact family).
    fn clock_root(&self, c_a: f64) -> f64 {
        let (mut lo, mut hi) = (1.0e-6_f64, 40.0_f64);
        assert!(self.clock_f(lo, c_a) < 0.0 && self.clock_f(hi, c_a) > 0.0, "clock bracket");
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if self.clock_f(mid, c_a) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

fn je(v: f64) -> String {
    format!("{v:.17e}")
}

fn main() {
    let t0 = Instant::now();
    let args: Vec<String> = std::env::args().collect();
    let outdir = args.get(1).cloned().unwrap_or_else(|| ".".to_string());
    let threads: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(4);

    println!("h26_dilation: N={N} LBOX={LBOX} t={T_FROZEN:.15} threads={threads}");

    // ---- the SAME static stage as h26_convention ---------------------------
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    let mut f = Field3::new(N, LBOX);
    f.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let (out_h, _) = fs_gum_kern::eval(&f, &Opts::estatic(Scheme::Corner), false, threads);
    let deg_ref = out_h.deg;
    let fgap_ref = out_h.e6 + out_h.e0 - fs_gum_field::bps_floor() - FLOOR_BAND;
    let bump = bump_field(N, LBOX, SEED_PERT, 6, 0.02);
    add_scaled(&mut f, &bump, 1.0);
    let _ = f.renormalize();
    let opts_static = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
    let res_s = fs_gum_kern::anf(&mut f, &opts_static, &AnfParams::new(MAXIT_STATIC), "static", threads);
    let o = res_s.out;
    println!(
        "static anchor: E2={:.6} E4={:.6} E6={:.6} E0={:.6} I={:.5} Estat={:.7} ({:.1}s)",
        o.e2, o.e4, o.e6, o.e0, o.i, o.estat, t0.elapsed().as_secs_f64()
    );

    let fam = Fam { t: T_FROZEN, e2h: o.e2, e4h: o.e4, e6h: o.e6, e0h: o.e0, ih: o.i };

    // the family's own static optimum (L = 0) — the V reference
    let s0 = fam.sstar(0.0);
    println!("family static optimum: s0={s0:.9} (V0={:.9}), Estat(s0)={:.7}", s0 * s0 * s0, fam.estat(s0));

    let a_semi = 0.5_f64;
    let a_qr = 0.75_f64.sqrt();
    let c_of = |a: f64| 2.0 * a / (2.0 - a);
    let sqrt2pi3 = (2.0 * PI * PI * PI).sqrt();
    let kthr = kappa_threshold();

    let mut anchors_json = Vec::new();
    let mut l_corpus = 0.0_f64;
    for (name, a) in [("corpus_L_eq_j", a_semi), ("alt_L2_eq_jjp1", a_qr)] {
        let c_a = c_of(a);
        let l = fam.clock_root(c_a);
        let s = fam.sstar(l);
        let (iv, es) = (fam.i(s), fam.estat(s));
        let erot = l * l / (2.0 * iv);
        let etot = es + erot;
        let erot_frac = erot / etot;
        let v_rel = (s / s0).powi(3);
        let kappa = l / iv;
        let kappa_paper = 2.0 * PI.sqrt() * kappa;
        let c_paper = (l / a) / sqrt2pi3;
        // envelope-theorem check: dE/dL along the envelope vs omega = L/I
        let dl = 1.0e-6 * l;
        let dedl = (fam.e(fam.sstar(l + dl), l + dl) - fam.e(fam.sstar(l - dl), l - dl)) / (2.0 * dl);
        let env_check = (dedl / (l / iv) - 1.0).abs();
        println!(
            "[{name}] a={a:.10}: L*={l:.9} s*={s:.9} V/V0={v_rel:.9} Erot/E={erot_frac:.12} (a/2={:.12})",
            a / 2.0
        );
        println!(
            "         kappa_ours={kappa:.9} = {:.6}x halo thr; kappa_paper={kappa_paper:.9}; c_paper={c_paper:.9}; envelope check |dE/dL / omega - 1|={env_check:.2e}",
            kappa / kthr
        );
        if name == "corpus_L_eq_j" {
            l_corpus = l;
        }
        anchors_json.push(format!(
            "\"{name}\": {{\"a\": {}, \"c_a\": {}, \"L\": {}, \"s_star\": {}, \"V_over_V0\": {}, \"I\": {}, \"Estat\": {}, \"Etot\": {}, \"Erot_over_Etot\": {}, \"kappa_ours\": {}, \"kappa_over_halo_threshold\": {}, \"kappa_paper\": {}, \"c_paper\": {}, \"envelope_theorem_dev\": {}}}",
            je(a), je(c_a), je(l), je(s), je(v_rel), je(iv), je(es), je(etot), je(erot_frac),
            je(kappa), je(kappa / kthr), je(kappa_paper), je(c_paper), je(env_check)
        ));
    }

    // ---- H2.7(a) exponents on the exact family -----------------------------
    println!("\n[H2.7a/dilation] exponents at the corpus fixed point L*={l_corpus:.9} (hbar = 2L)");
    let mut exps_json = Vec::new();
    for delta in [1.0e-3_f64, 1.0e-5] {
        let lp = l_corpus * (1.0 + delta);
        let lm = l_corpus * (1.0 - delta);
        let (sp, sm) = (fam.sstar(lp), fam.sstar(lm));
        let dlnh = ((1.0 + delta) / (1.0 - delta)).ln();
        let ex_mech = ((lp / fam.i(sp)) / (lm / fam.i(sm))).ln() / dlnh;
        let ex_phase = ((fam.e(sp, lp) / (2.0 * lp)) / (fam.e(sm, lm) / (2.0 * lm))).ln() / dlnh;
        let ex_i = (fam.i(sp) / fam.i(sm)).ln() / dlnh;
        println!(
            "  delta={delta:.0e}: dln omega_mech/dln hbar = {ex_mech:+.9}   dln omega_phase/dln hbar = {ex_phase:+.9}   (dln I = {ex_i:+.9})"
        );
        exps_json.push(format!(
            "{{\"delta\": {}, \"dln_omega_mech_dln_h\": {}, \"dln_omega_phase_dln_h\": {}, \"dln_I_dln_h\": {}}}",
            je(delta), je(ex_mech), je(ex_phase), je(ex_i)
        ));
    }

    let body = format!(
        "{{\n\"phase\": \"H2.6 tier-2 — exact dilation-family closure (engine-anchored), + H2.7a exponents on the family\",\n\"date\": \"2026-07-16\",\n\"anchor_state\": {{\"protocol\": \"N=48 LBOX=4.5 eps=0.05 guarded static ANF, seed 424242, maxit 900 (identical to h26_convention static stage)\", \"E2\": {}, \"E4\": {}, \"E6\": {}, \"E0\": {}, \"I\": {}, \"Estat\": {}, \"deg\": {}}},\n\"family\": \"E(s;L) = t(s E2h + E4h/s) + E6h/s^3 + s^3 E0h + L^2/(2 Ih s^3); V = s^3 (exact sector scaling laws)\",\n\"s0_static_optimum\": {},\n\"anchors\": {{{}, {}}},\n\"h27a_exponents_dilation\": [{}],\n\"expected\": {{\"mech\": 0.5, \"phase\": -0.5, \"note\": \"T2 (e.1) exact on the pure V-family at x=1; deviations here are the eps = 0.05 E2/E4 corrections, measured\"}},\n\"runtime_seconds\": {}\n}}\n",
        je(o.e2), je(o.e4), je(o.e6), je(o.e0), je(o.i), je(o.estat), je(o.deg),
        je(s0), anchors_json[0], anchors_json[1], exps_json.join(", "),
        je(t0.elapsed().as_secs_f64())
    );
    let path = format!("{outdir}/h26_dilation.json");
    std::fs::write(&path, &body).expect("write json");
    println!("\nwrote {path}");
    println!("total {:.1}s", t0.elapsed().as_secs_f64());
}
