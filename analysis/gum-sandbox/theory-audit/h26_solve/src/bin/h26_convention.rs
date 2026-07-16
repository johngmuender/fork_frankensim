//! Phase H2.6 — the j(j+1)-convention NON-RIGID closure (ROADMAP_v5 item 6),
//! plus the Phase H2.7(a) field-level entrainment-exponent measurement.
//!
//! BACKGROUND (T2_closure_pillars.md §a.3, the GAP under audit): the corpus
//! anchors the spin condition semiclassically, L = j*c (c = the closure
//! number in engine action units), and never argues the choice against the
//! quantum-rotor alternative L^2 = j(j+1) c^2 used by its own imported
//! Skyrme literature.  T2 proved at the RIGID level that the j = 1/2
//! SELECTION is convention-robust but the closure TUPLE is not:
//! writing a = j (semiclassical) or a = sqrt(j(j+1)) (quantum rotor), the
//! identical T-B1 algebra gives w*a(1-a) = 1, V^2 = 1/(1-a), E_rot/E = a/2,
//! so at j = 1/2: a = 1/2 -> (V = 1.4142, E_rot/E = 1/4) vs
//! a = sqrt(3)/2 -> (V = 2.7321, E_rot/E = 0.4330).
//!
//! THE MODIFIED CLOCK CONDITION, derived once (this is the whole algebra):
//! spin anchor L = a*c and clock E_tot = c*omega with omega = L/I give
//! E_tot = (L/a)(L/I) = L^2/(a I); with E_tot = E_static + L^2/(2I),
//!
//!     E_static = L^2 (2 - a) / (2 a I)
//!     =>  L^2 = [2a/(2-a)] * I * E_static           (clock condition C(a))
//!     =>  E_rot/E_tot = a/2  IDENTICALLY at any root of C(a).
//!
//! a = 1/2 recovers the campaign's frozen Section-K form
//! L^2 = (2/3) I E_static (E_rot/E = 1/4); a = sqrt(3)/2 gives
//! L^2 = [2*sqrt(3)/(4-sqrt(3))] I E_static = 1.5274601... I E_static
//! (E_rot/E = sqrt(3)/4 = 0.4330127).  Because a/2 is kinematically forced
//! AT an exact root, the field-level discriminator is (i) whether the
//! non-rigid solve under each anchor HAS a root inside the guarded family
//! at all, (ii) the rest of the solved tuple (kappa vs the halo/onset
//! thresholds, the paper-unit c-analogue, the inertia dilation), and
//! (iii) whether the two anchors' solved E_rot/E stay separated at the
//! corpus's own measurement precision (<r1>: 0.2500 +- 0.0002).
//!
//! NON-RIGID means: at every candidate L the FIELD relaxes under the
//! guarded Routhian objective (E_static + L^2/2I + frozen guards) before
//! (I, E_static) are read — the field3d/fs-gum-statics 4A machinery, not
//! the rigid V-dilation family.  Every relaxation cold-starts from the SAME
//! static solution, so (I(L), E_static(L)) is a pure deterministic function
//! of L and the bisection is well-posed.
//!
//! PROTOCOL (the campaign's frozen N = 48 gate protocol):
//! LBOX = 4.5, eps = 0.05 (t = T_FROZEN), corner scheme, guards
//! MU = 5000/band 0.005 anchored at the engine hedgehog degree, wall
//! MU_F = 400 at hedgehog gap - 0.01; static stage = hedgehog +
//! bump(424242, K = 6, amp = 0.02), 900 guarded ANF iterations (gates G-B);
//! clock stage = cold restart from the static solution, MAXIT_L guarded
//! Routhian iterations per bisection evaluation, MAXIT_FINAL for the
//! reported endpoints; bisection = the Section-K bisection transplanted to
//! the self-consistent target (fixed step count on a scanned bracket),
//! then a short fixed-point polish L <- sqrt(c_a I(L) E(L)) at
//! MAXIT_FINAL to shrink the finite-relaxation residual (recorded).
//! Sweeps run on the fs-gum-kern deterministic tiled layer (bit-identical
//! at any thread count by construction; a replay of the final corpus-anchor
//! relaxation is asserted bitwise in-process).
//!
//! H2.7(a) (IV.H.3 error-signal exponents, field level): displace the
//! effective hbar-anchor by delta — realized exactly as the task states,
//! clock target L -> L*(1+delta) at fixed shape protocol, hbar = L/j = 2L —
//! re-relax, and measure d ln omega_mech / d ln hbar (omega_mech = L/I on
//! the solved state) and d ln omega_phase / d ln hbar (omega_phase =
//! E_tot/hbar, the IV.H.1(ii) Nelson-stationary phase rate) by central
//! log-differences at delta = 1e-3 and 3e-3.  T2 proved +-1/2 exactly on
//! the reduced family E(L) = e*sqrt(1+x) at the fixed point x = 1; this is
//! the field-level empirical check.
//!
//! Unit maps (SESSION_HANDOFF §3, compacton-exact, validated 7-9 digits):
//! c_paper = (L/a)/sqrt(2 pi^3); kappa_ours = L/I; kappa_paper =
//! 2 sqrt(pi) kappa_ours (halo threshold kappa_paper = 1/sqrt(2), i.e.
//! kappa_ours = 1/sqrt(8 pi); corpus measured radiation onset
//! kappa_paper = 1.000 +- 0.004).
//!
//! Usage: h26_convention <outdir> [maxit_l] [maxit_final] [bisect_steps] [threads]
//! (defaults 350 700 16 4).  Same binary + args => bit-identical JSON.
//!
//! Epistemic notice (binding): within-model numerical engineering on a
//! speculative theory's functional.  Measured facts about the campaign's
//! own certified functional; nothing here says anything about nature.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{add_scaled, bump_field, halo_fraction, restore_cells, save_cells};
use fs_gum_statics::{kappa_threshold, AnfParams, Opts, FLOOR_BAND};

use std::fmt::Write as _;
use std::time::Instant;

const LBOX: f64 = 4.5;
const N: usize = 48;
const MAXIT_STATIC: usize = 900;
const SEED_PERT: u64 = 424242;
const PI: f64 = std::f64::consts::PI;

/// One relaxed-at-fixed-L measurement (endpoint of the guarded Routhian ANF).
#[derive(Clone, Copy, Debug)]
struct Solved {
    l: f64,
    i: f64,
    estat: f64,
    etot: f64,
    erot_frac: f64,
    kappa: f64,
    halo: f64,
    deg: f64,
    floor_gap: f64,
    epen: f64,
    efpen: f64,
    gn_ratio: f64,
    arrests: usize,
    iters: usize,
}

fn je(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.17e}")
    } else {
        format!("\"{v}\"")
    }
}

fn solved_json(s: &Solved) -> String {
    format!(
        "{{\"L\": {}, \"I\": {}, \"Estat\": {}, \"Etot\": {}, \"Erot_over_Etot\": {}, \"kappa_ours\": {}, \"halo\": {}, \"deg\": {}, \"floor_gap\": {}, \"epen\": {}, \"efpen\": {}, \"gn_ratio\": {}, \"arrests\": {}, \"iters\": {}}}",
        je(s.l), je(s.i), je(s.estat), je(s.etot), je(s.erot_frac), je(s.kappa), je(s.halo),
        je(s.deg), je(s.floor_gap), je(s.epen), je(s.efpen), je(s.gn_ratio), s.arrests, s.iters
    )
}

struct Ctx {
    q_static: Vec<f64>,
    deg_ref: f64,
    fgap_ref: f64,
    threads: usize,
}

impl Ctx {
    /// Cold-restart from the static solution, relax under the guarded
    /// Routhian at fixed L for `maxit` iterations, measure the endpoint.
    fn relax_at(&self, f: &mut Field3, l: f64, maxit: usize, label: &str) -> Solved {
        self.relax_from(f, &self.q_static, l, maxit, label)
    }

    /// Same, but warm-started from an arbitrary saved state (used for the
    /// H2.7(a) +-delta displacements about the solved fixed point: both
    /// signs inherit the SAME residual unconvergence, which cancels to
    /// first order in the central log-difference).
    fn relax_from(&self, f: &mut Field3, start: &[f64], l: f64, maxit: usize, label: &str) -> Solved {
        restore_cells(f, start);
        let opts = Opts::routhian(Scheme::Corner, l).with_guards(self.deg_ref, self.fgap_ref);
        let prm = AnfParams::new(maxit);
        let res = fs_gum_kern::anf(f, &opts, &prm, label, self.threads);
        let out = res.out;
        let (halo, _c) = halo_fraction(f);
        let erot = l * l / (2.0 * out.i);
        Solved {
            l,
            i: out.i,
            estat: out.estat,
            etot: out.estat + erot,
            erot_frac: erot / (out.estat + erot),
            kappa: l / out.i,
            halo,
            deg: out.deg,
            floor_gap: out.floor_gap,
            epen: out.epen,
            efpen: out.efpen,
            gn_ratio: res.gn / res.gn0,
            arrests: res.arrests,
            iters: res.iters,
        }
    }
}

/// Clock mismatch under anchor a at the relaxed state:
/// F(L) = L^2 - c_a * I(L) * E_static(L)  (root <=> clock condition C(a)).
fn clock_f(s: &Solved, c_a: f64) -> f64 {
    s.l * s.l - c_a * s.i * s.estat
}

fn main() {
    let t0 = Instant::now();
    let args: Vec<String> = std::env::args().collect();
    let outdir = args.get(1).cloned().unwrap_or_else(|| ".".to_string());
    let maxit_l: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(350);
    let maxit_final: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(700);
    let bisect_steps: usize = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(16);
    let threads: usize = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(4);

    println!(
        "h26_convention: N={N} LBOX={LBOX} t={T_FROZEN:.15} maxit_l={maxit_l} maxit_final={maxit_final} bisect={bisect_steps} threads={threads}"
    );

    // ---- anchors ----------------------------------------------------------
    // a = j (semiclassical, corpus) and a = sqrt(j(j+1)) (quantum rotor), j = 1/2.
    let a_semi = 0.5_f64;
    let a_qr = 0.75_f64.sqrt(); // sqrt(3)/2
    let c_of = |a: f64| 2.0 * a / (2.0 - a);
    // rigid-level regenerated lines (t2_convention.py class, T2 §a.3):
    let e0hat = 64.0 / (15.0 * PI);
    let i0hat = 256.0 / (105.0 * PI);
    let rigid = |a: f64| -> (f64, f64, f64, f64, f64) {
        let w = 1.0 / (a * (1.0 - a));
        let v = (1.0 / (1.0 - a)).sqrt();
        let c_g1 = (w * e0hat * i0hat).sqrt();
        let c_g32 = (w * e0hat * i0hat * 1.5).sqrt();
        (w, v, a / 2.0, c_g1, c_g32)
    };
    let (w_s, v_s, er_s, c1_s, c32_s) = rigid(a_semi);
    let (w_q, v_q, er_q, c1_q, c32_q) = rigid(a_qr);
    println!("rigid-level (regenerated t2_convention lines):");
    println!(
        "  L=j c      : w={w_s:.6} V={v_s:.6} Erot/E={er_s:.6} c(g=1)={c1_s:.6} c(g=3/2)={c32_s:.6}"
    );
    println!(
        "  L^2=j(j+1) : w={w_q:.6} V={v_q:.6} Erot/E={er_q:.6} c(g=1)={c1_q:.6} c(g=3/2)={c32_q:.6}"
    );
    println!(
        "  clock coefficients c_a = 2a/(2-a): corpus {:.15} (=2/3), alternative {:.15}",
        c_of(a_semi),
        c_of(a_qr)
    );

    // ---- radial profile + hedgehog references (gates G-B protocol) -------
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    let mut f = Field3::new(N, LBOX);
    f.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let (out_h, _) = fs_gum_kern::eval(&f, &Opts::estatic(Scheme::Corner), false, threads);
    let deg_ref = out_h.deg;
    let fgap_h = out_h.e6 + out_h.e0 - fs_gum_field::bps_floor();
    let fgap_ref = fgap_h - FLOOR_BAND;
    println!(
        "hedgehog (corner N={N}): Estat={:.7} deg={:.6} I={:.5} -> deg_ref={deg_ref:.6} fgap_ref={fgap_ref:+.6}",
        out_h.estat, out_h.deg, out_h.i
    );

    // ---- static stage ------------------------------------------------------
    let bump = bump_field(N, LBOX, SEED_PERT, 6, 0.02);
    add_scaled(&mut f, &bump, 1.0);
    let _ = f.renormalize();
    let opts_static = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
    let prm_s = AnfParams::new(MAXIT_STATIC);
    let res_s = fs_gum_kern::anf(&mut f, &opts_static, &prm_s, "static", threads);
    let mut q_static = vec![0.0_f64; 4 * N * N * N];
    save_cells(&f, &mut q_static);
    let (i_s, e_s) = (res_s.out.i, res_s.out.estat);
    println!(
        "static: status={} iters={} arrests={} Estat={:.7} I={:.5} deg={:.6} ({:.1}s)",
        res_s.status.as_str(),
        res_s.iters,
        res_s.arrests,
        e_s,
        i_s,
        res_s.out.deg,
        t0.elapsed().as_secs_f64()
    );

    let ctx = Ctx { q_static, deg_ref, fgap_ref, threads };
    let mut log: Vec<(String, Solved, f64)> = Vec::new(); // (tag, state, F)

    // ---- per-anchor non-rigid clock solve ---------------------------------
    let sqrt2pi3 = (2.0 * PI * PI * PI).sqrt();
    let kthr = kappa_threshold();
    let mut anchor_json = Vec::new();
    let mut roots = Vec::new();
    let mut q_root: Option<Vec<f64>> = None;
    for (name, a) in [("corpus_L_eq_j", a_semi), ("alt_L2_eq_jjp1", a_qr)] {
        let ta = Instant::now();
        let c_a = c_of(a);
        let l0 = (c_a * i_s * e_s).sqrt();
        println!("\n[{name}] a={a:.12} c_a={c_a:.12}  rigid-from-static L0={l0:.6}");
        let mut eval_at = |l: f64, maxit: usize, tag: &str| -> (Solved, f64) {
            let s = ctx.relax_at(&mut f, l, maxit, tag);
            let fv = clock_f(&s, c_a);
            println!(
                "  {tag}: L={:.6} -> I={:.5} Estat={:.6} F={:+.5} kappa={:.5} halo={:.4} deg={:.5} gn/gn0={:.2e} arr={}",
                s.l, s.i, s.estat, fv, s.kappa, s.halo, s.deg, s.gn_ratio, s.arrests
            );
            log.push((format!("{name}/{tag}"), s, fv));
            (s, fv)
        };
        // scan bracket around the rigid-from-static estimate
        let (mut lo, mut hi) = (0.75 * l0, 1.25 * l0);
        let (_slo0, mut flo) = eval_at(lo, maxit_l, "scan_lo");
        let (_smid, _fmid) = eval_at(l0, maxit_l, "scan_mid");
        let (_shi0, mut fhi) = eval_at(hi, maxit_l, "scan_hi");
        let mut widen = 0;
        while flo > 0.0 && widen < 4 {
            lo *= 0.8;
            let (_s, fv) = eval_at(lo, maxit_l, "widen_lo");
            flo = fv;
            widen += 1;
        }
        while fhi < 0.0 && widen < 8 {
            hi *= 1.25;
            let (_s, fv) = eval_at(hi, maxit_l, "widen_hi");
            fhi = fv;
            widen += 1;
        }
        assert!(flo < 0.0 && fhi > 0.0, "clock bracket failed for {name}: F({lo})={flo}, F({hi})={fhi}");
        // the campaign's standard bisection, transplanted to the
        // self-consistent target (fixed step count, monotone bracket)
        for step in 0..bisect_steps {
            let mid = 0.5 * (lo + hi);
            let (_s, fv) = eval_at(mid, maxit_l, &format!("bisect_{step:02}"));
            if fv > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let l_bisect = 0.5 * (lo + hi);
        // fixed-point polish at the reporting budget
        let mut l_cur = l_bisect;
        let mut fin = None;
        for p in 0..3 {
            let (s, _fv) = eval_at(l_cur, maxit_final, &format!("polish_{p}"));
            l_cur = (c_a * s.i * s.estat).sqrt();
            fin = Some(s);
        }
        let (s_fin, f_fin) = eval_at(l_cur, maxit_final, "final");
        let _ = fin; // last polish state superseded by the final evaluation
        let resid_rel = f_fin / (s_fin.l * s_fin.l);
        let c_paper = (s_fin.l / a) / sqrt2pi3;
        let kappa_paper = 2.0 * PI.sqrt() * s_fin.kappa;
        println!(
            "[{name}] SOLVED: L*={:.6} Erot/E={:.9} (target a/2={:.9}, resid_rel={:+.2e})",
            s_fin.l,
            s_fin.erot_frac,
            a / 2.0,
            resid_rel
        );
        println!(
            "         kappa_ours={:.6} = {:.4}x halo threshold; kappa_paper={:.6} (onset 1.000+-0.004, ring 1/sqrt2=0.7071)",
            s_fin.kappa,
            s_fin.kappa / kthr,
            kappa_paper
        );
        println!(
            "         c_paper=(L/a)/sqrt(2 pi^3)={:.6}; I*/I_static={:.5}; halo={:.4}  ({:.1}s)",
            c_paper,
            s_fin.i / i_s,
            s_fin.halo,
            ta.elapsed().as_secs_f64()
        );
        anchor_json.push(format!(
            "\"{name}\": {{\"a\": {}, \"c_a\": {}, \"L_bisect\": {}, \"solved\": {}, \"clock_resid_rel\": {}, \"erot_target_a_over_2\": {}, \"kappa_over_halo_threshold\": {}, \"kappa_paper\": {}, \"c_paper\": {}, \"I_over_I_static\": {}}}",
            je(a), je(c_a), je(l_bisect), solved_json(&s_fin), je(resid_rel), je(a / 2.0),
            je(s_fin.kappa / kthr), je(kappa_paper), je(c_paper), je(s_fin.i / i_s)
        ));
        // keep the endpoint field of the corpus anchor for the H2.7(a) warm starts
        if name == "corpus_L_eq_j" {
            // f still holds the "final" relaxation endpoint
            let mut buf = vec![0.0_f64; 4 * N * N * N];
            save_cells(&f, &mut buf);
            q_root = Some(buf);
        }
        roots.push((name, a, s_fin));
    }
    let q_root = q_root.expect("corpus root state saved");

    // ---- H2.7(a): +-1/2 exponents at the corpus fixed point ----------------
    // Warm-started from the solved fixed-point state (common residual
    // unconvergence cancels to first order in the central difference);
    // one cold-start pair kept as a protocol cross-check, plus a
    // relaxation-budget scan so the exponent's convergence trend is itself
    // a measured fact.
    let (_, _, s_star) = roots[0];
    let l_star = s_star.l;
    println!("\n[H2.7a] entrainment exponents at the corpus fixed point L*={l_star:.6} (hbar = 2L)");
    let mut h27 = Vec::new();
    let mut measure_pair = |start: &[f64], delta: f64, maxit: usize, tag: &str| -> String {
        let sp = ctx.relax_from(&mut f, start, l_star * (1.0 + delta), maxit, "delta+");
        let sm = ctx.relax_from(&mut f, start, l_star * (1.0 - delta), maxit, "delta-");
        let dlnh = ((1.0 + delta) / (1.0 - delta)).ln();
        let ex_mech = ((sp.l / sp.i) / (sm.l / sm.i)).ln() / dlnh;
        let ex_phase = ((sp.etot / (2.0 * sp.l)) / (sm.etot / (2.0 * sm.l))).ln() / dlnh;
        let ex_i = (sp.i / sm.i).ln() / dlnh;
        let ex_etot = (sp.etot / sm.etot).ln() / dlnh;
        println!(
            "  {tag} delta={delta:.0e} maxit={maxit}: dln omega_mech/dln hbar = {ex_mech:+.6}   dln omega_phase/dln hbar = {ex_phase:+.6}   (dln I = {ex_i:+.6}, dln Etot = {ex_etot:+.6})"
        );
        format!(
            "{{\"variant\": \"{tag}\", \"delta\": {}, \"maxit\": {maxit}, \"plus\": {}, \"minus\": {}, \"dln_omega_mech_dln_h\": {}, \"dln_omega_phase_dln_h\": {}, \"dln_I_dln_h\": {}, \"dln_Etot_dln_h\": {}}}",
            je(delta), solved_json(&sp), solved_json(&sm), je(ex_mech), je(ex_phase), je(ex_i), je(ex_etot)
        )
    };
    // primary: warm-started, both deltas, reporting budget
    for delta in [1.0e-3_f64, 3.0e-3] {
        h27.push(measure_pair(&q_root, delta, maxit_final, "warm"));
    }
    // budget scan (warm), delta = 3e-3
    for m in [maxit_l, 2 * maxit_final, 4 * maxit_final] {
        h27.push(measure_pair(&q_root, 3.0e-3, m, "warm_budget"));
    }
    // protocol cross-check: cold-started from the static state
    h27.push(measure_pair(&ctx.q_static.clone(), 3.0e-3, maxit_final, "cold"));

    // ---- determinism replay -------------------------------------------------
    let s_replay = ctx.relax_at(&mut f, l_star, maxit_final, "replay");
    let replay_ok = s_replay.i.to_bits() == s_star.i.to_bits()
        && s_replay.estat.to_bits() == s_star.estat.to_bits()
        && s_replay.deg.to_bits() == s_star.deg.to_bits();
    println!(
        "\nreplay of the corpus final relaxation: bitwise (I, Estat, deg) identical = {replay_ok}"
    );
    assert!(replay_ok, "determinism replay failed");

    // ---- discriminator verdict ---------------------------------------------
    let er_c = roots[0].2.erot_frac;
    let er_a = roots[1].2.erot_frac;
    let sep_sigma = (er_a - er_c).abs() / 2.0e-4; // <r1> quoted sigma on Erot/E
    println!(
        "\nDISCRIMINATOR: corpus-anchor solved Erot/E = {er_c:.6}; alternative-anchor solved Erot/E = {er_a:.6}"
    );
    println!(
        "  separation = {:.4} = {sep_sigma:.0} x the <r1> quoted sigma (0.0002); <r1> measured 0.2500 +- 0.0002",
        er_a - er_c
    );

    // ---- JSON ---------------------------------------------------------------
    let mut evals_json = String::new();
    for (i, (tag, s, fv)) in log.iter().enumerate() {
        let _ = write!(
            evals_json,
            "{}\n  {{\"tag\": \"{tag}\", \"F\": {}, \"state\": {}}}",
            if i == 0 { "" } else { "," },
            je(*fv),
            solved_json(s)
        );
    }
    let mut exps_json = String::new();
    for (i, r) in h27.iter().enumerate() {
        let _ = write!(exps_json, "{}\n  {r}", if i == 0 { "" } else { "," });
    }
    let body = format!(
        "{{\n\"phase\": \"H2.6 j(j+1)-convention non-rigid closure (+ H2.7a field exponents)\",\n\"date\": \"2026-07-16\",\n\"protocol\": {{\"N\": {N}, \"LBOX\": {LBOX}, \"t_frozen\": {}, \"eps\": 0.05, \"scheme\": \"corner\", \"guards\": {{\"deg_ref\": {}, \"fgap_ref\": {}, \"mu_deg\": 5000.0, \"mu_floor\": 400.0}}, \"static\": {{\"seed\": {SEED_PERT}, \"kbumps\": 6, \"amp\": 0.02, \"maxit\": {MAXIT_STATIC}, \"I\": {}, \"Estat\": {}}}, \"maxit_l\": {maxit_l}, \"maxit_final\": {maxit_final}, \"bisect_steps\": {bisect_steps}, \"threads\": {threads}, \"kernel\": \"fs-gum-kern tiled deterministic layer (bit-identical at any thread count)\", \"clock_family\": \"L^2 = [2a/(2-a)] I(L) E_static(L), field relaxed at fixed L from the common static state\"}},\n\"rigid_level_regenerated\": {{\"e0hat\": {}, \"i0hat\": {}, \"semiclassical\": {{\"a\": {}, \"w\": {}, \"V\": {}, \"Erot_over_E\": {}, \"c_g1\": {}, \"c_g32\": {}}}, \"quantum_rotor\": {{\"a\": {}, \"w\": {}, \"V\": {}, \"Erot_over_E\": {}, \"c_g1\": {}, \"c_g32\": {}}}}},\n\"anchors\": {{{}, {}}},\n\"discriminator\": {{\"erot_corpus\": {}, \"erot_alternative\": {}, \"separation\": {}, \"r1_quoted\": [0.2500, 0.0002], \"separation_in_r1_sigma\": {}}},\n\"h27a_exponents\": [{}\n],\n\"unit_maps\": {{\"c_paper\": \"(L/a)/sqrt(2 pi^3)\", \"kappa_paper\": \"2 sqrt(pi) L/I\", \"halo_threshold_ours\": {}, \"halo_threshold_paper\": 0.7071067811865475}},\n\"determinism\": {{\"replay_bitwise_identical\": {replay_ok}}},\n\"evaluations\": [{}\n],\n\"runtime_seconds\": {}\n}}\n",
        je(T_FROZEN), je(deg_ref), je(fgap_ref), je(i_s), je(e_s),
        je(e0hat), je(i0hat),
        je(a_semi), je(w_s), je(v_s), je(er_s), je(c1_s), je(c32_s),
        je(a_qr), je(w_q), je(v_q), je(er_q), je(c1_q), je(c32_q),
        anchor_json[0], anchor_json[1],
        je(er_c), je(er_a), je(er_a - er_c), je(sep_sigma),
        exps_json, je(kappa_threshold()), evals_json,
        je(t0.elapsed().as_secs_f64())
    );
    let path = format!("{outdir}/h26_results.json");
    std::fs::write(&path, &body).expect("write json");
    println!("\nwrote {path}");
    // the h27(a) block again under the h27 prefix for the H2.7 assembler
    let mut h27_body = String::from("{\n\"phase\": \"H2.7a field-level entrainment exponents (produced by h26_convention)\",\n\"L_star\": ");
    let _ = write!(
        h27_body,
        "{},\n\"fixed_point\": {},\n\"exponents\": [{}\n],\n\"expected\": {{\"mech\": 0.5, \"phase\": -0.5, \"proof_scope\": \"T2 (e.1): exact on the reduced family at x = 1; this is the field-level check\"}}\n}}\n",
        je(l_star),
        solved_json(&s_star),
        exps_json
    );
    let path27 = format!("{outdir}/h27a_field_exponents.json");
    std::fs::write(&path27, &h27_body).expect("write json");
    println!("wrote {path27}");
    println!("total {:.1}s", t0.elapsed().as_secs_f64());
}
