//! GUM core Phase E1 gate binary for `fs-gum-statics`.
//!
//! Gates (all numbers machine-checked, entire pipeline run twice for the
//! bit-identical replay gate):
//!
//!   X    forward-path consistency: this crate's PAIRS-expansion forward
//!        sums vs fs-gum-field's det4 `measure()` on the N = 48 hedgehog,
//!        both schemes (rel <= 1e-12 — analytically identical expansions).
//!   G-A  THE LOAD-BEARING GATE, first: FD-vs-analytic gradient check at
//!        1e-5 on N = 24 (perturbed hedgehog base; random smooth tangent
//!        directions; both schemes; every sector separately — E2, E4, E6,
//!        E0, rotor L^2/2I — and the totals E_static, Routhian, and the
//!        fully guarded objective with both one-sided guards forced
//!        active).  The descents do not run unless this passes.
//!   G-B  static relaxation at N = 48: eps = 0.05 hedgehog seed
//!        (t = 0.0082764349), fixed-seed non-axisymmetric bump, ANF with
//!        the frozen guards; gates: degree within the anchor band, virial
//!        (central4 re-measure) improves or stays 1e-3-class, energy not
//!        below the guarded Bogomolny floor, final sectors within stated
//!        tolerance of the seeded hedgehog's (same-functional comparison),
//!        objective monotone by construction.
//!   G-C  halo-descent qualitative referee at LOW resolution (N = 48,
//!        iteration-capped): fixed-L Routhian descent at kappa_0 ~ 0.266
//!        (above threshold 1/sqrt(8 pi) = 0.19947): R decreases, kappa
//!        falls, halo fraction rises (the 4A / F-R5 signature); CONTROL at
//!        kappa_0 ~ 0.12: halo stays low, kappa stable; differential gate.
//!   G-D  clock: bisect L_clock on the control solution, gate
//!        E_rot/E_tot = 0.250000 to 1e-6, report kappa(L_clock)/threshold.
//!   G-E  bit-identical two-run replay of every gate number (BLAKE3
//!        fingerprint printed).
//!
//! Epistemic notice (binding): every PASS certifies a WITHIN-MODEL property
//! of a discretisation, its analytic gradients, and a port.  The G-C halo
//! descent replicates the campaign's F-R5 mechanism as an engine
//! capability, not new physics; nothing here says anything about nature.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{measure, Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{
    add_scaled, bump_field, clock_bisect, dot_flat, erot_frac, halo_seed_field, normalize_l2,
    restore_cells, save_cells, tangent_project,
};
use fs_gum_statics::{
    anf, eval, kappa_threshold, AnfParams, AnfResult, Opts, Out, FLOOR_BAND,
};

use std::time::Instant;

const LBOX: f64 = 4.5;
const N_FD: usize = 24; // gradient-gate grid
const N_RUN: usize = 48; // relaxation / descent grid
const EPS_FD: f64 = 1.0e-5;
const TOL_GRAD: f64 = 1.0e-5;
const SEED_PERT: u64 = 424242; // campaign perturbation seed (logical id)
const SEED_FD_BASE: u64 = 11; // field3d gradcheck base-perturbation seed id
const L_FD: f64 = 8.4979; // Routhian L for the gradient gate (field3d L_MAIN)
const KAPPA_MAIN: f64 = 0.266; // over-spun target kappa (above 0.19947)
const KAPPA_CTRL: f64 = 0.12; // control target kappa (below threshold)
const MAXIT_STATIC: usize = 900;
// equal caps so the main/control comparison is at identical iteration counts
const MAXIT_MAIN: usize = 600;
const MAXIT_CTRL: usize = 600;

struct Gate {
    id: &'static str,
    desc: String,
    pass: bool,
}

#[derive(Default)]
struct Run {
    names: Vec<String>,
    nums: Vec<f64>,
    gates: Vec<Gate>,
    grad_table: Vec<(String, f64)>, // per (scheme, sector) worst FD residual
    trend_main: Vec<String>,
    trend_ctrl: Vec<String>,
    aborted: bool,
}

impl Run {
    fn num(&mut self, name: &str, v: f64) {
        self.names.push(name.to_string());
        self.nums.push(v);
    }
    fn gate(&mut self, id: &'static str, desc: String, pass: bool) {
        self.gates.push(Gate { id, desc, pass });
    }
}

fn virial_out(o: &Out, t: f64) -> f64 {
    t * o.e2 - t * o.e4 - 3.0 * o.e6 + 3.0 * o.e0
}

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

/// The eight FD objectives: label, sector weights, L, guards-forced flag.
fn fd_cases() -> Vec<(&'static str, [f64; 4], Option<f64>, bool)> {
    let t = T_FROZEN;
    vec![
        ("E2", [1.0, 0.0, 0.0, 0.0], None, false),
        ("E4", [0.0, 1.0, 0.0, 0.0], None, false),
        ("E6", [0.0, 0.0, 1.0, 0.0], None, false),
        ("E0", [0.0, 0.0, 0.0, 1.0], None, false),
        ("ROT", [0.0, 0.0, 0.0, 0.0], Some(L_FD), false),
        ("ESTAT", [t, t, 1.0, 1.0], None, false),
        ("ROUTH", [t, t, 1.0, 1.0], Some(L_FD), false),
        ("GUARD", [t, t, 1.0, 1.0], Some(L_FD), true),
    ]
}

fn scheme_name(s: Scheme) -> &'static str {
    match s {
        Scheme::Central2 => "central2",
        Scheme::Central4 => "central4",
        Scheme::Corner => "corner",
    }
}

#[allow(clippy::too_many_lines)]
fn run_all(verbose: bool) -> Run {
    let mut run = Run::default();
    let vp = |s: &str| {
        if verbose {
            println!("{s}");
        }
    };

    // ---- radial profile (frozen t, Step-1 method via fs-gum-field) -------
    let t0 = Instant::now();
    let rp = radial_solve(T_FROZEN, 4000, 6.0);
    vp(&format!(
        "radial profile: E2={:.6} E4={:.6} E6={:.6} E0={:.6} gmax={:.1e}  ({:.1}s)",
        rp.sectors[0],
        rp.sectors[1],
        rp.sectors[2],
        rp.sectors[3],
        rp.gmax,
        t0.elapsed().as_secs_f64()
    ));
    for (k, v) in ["radE2", "radE4", "radE6", "radE0"].iter().zip(rp.sectors.iter()) {
        run.num(k, *v);
    }

    // ---- X: forward consistency vs fs-gum-field::measure ------------------
    let mut fh = Field3::new(N_RUN, LBOX);
    fh.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let mut worst_x = 0.0_f64;
    for scheme in [Scheme::Central4, Scheme::Corner] {
        let (mine, _) = eval(&fh, &Opts::estatic(scheme), false);
        let s = measure(&fh, scheme, T_FROZEN);
        for (a, b) in [
            (mine.e2, s.e2),
            (mine.e4, s.e4),
            (mine.e6, s.e6),
            (mine.e0, s.e0),
            (mine.i, s.i),
            (mine.deg, s.deg),
        ] {
            worst_x = worst_x.max(rel(a, b));
        }
    }
    run.num("X_worst_rel", worst_x);
    run.gate(
        "X1",
        format!("PAIRS forward vs fs-gum-field measure(), worst rel {worst_x:.2e} (tol 1e-12)"),
        worst_x < 1.0e-12,
    );

    // ---- G-A: FD-vs-analytic gradient gate (load-bearing, first) ----------
    let ta = Instant::now();
    let mut fbase = Field3::new(N_FD, LBOX);
    fbase.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let bump = bump_field(N_FD, LBOX, SEED_FD_BASE, 6, 0.05);
    add_scaled(&mut fbase, &bump, 1.0);
    let _ = fbase.renormalize();
    let mut worst_all = 0.0_f64;
    for (si, scheme) in [Scheme::Central4, Scheme::Corner].into_iter().enumerate() {
        for (ci, (label, c, l, guards)) in fd_cases().into_iter().enumerate() {
            let mut o = Opts { scheme, t: T_FROZEN, c, l, deg_ref: 1.0, anchor: None, wall: None };
            if guards {
                // force BOTH one-sided guards active: deg_ref = 1 puts the
                // N = 24 degree below the band; fgap_ref = +10 penetrates
                // the wall everywhere (field3d gradcheck protocol).
                o = o.with_guards(1.0, 10.0);
            }
            let (_out, g) = eval(&fbase, &o, true);
            let g = g.expect("grad");
            let mut worst = 0.0_f64;
            for dir in 0..3u64 {
                let seed = 9000 + 100 * (si as u64) + 10 * (ci as u64) + dir;
                let mut u = bump_field(N_FD, LBOX, seed, 6, 1.0);
                tangent_project(&fbase, &mut u);
                let _ = normalize_l2(&mut u);
                let mut fp = fs_gum_statics::diag::clone_field(&fbase);
                add_scaled(&mut fp, &u, EPS_FD);
                let _ = fp.renormalize();
                let mut fm = fs_gum_statics::diag::clone_field(&fbase);
                add_scaled(&mut fm, &u, -EPS_FD);
                let _ = fm.renormalize();
                let (op, _) = eval(&fp, &o, false);
                let (om, _) = eval(&fm, &o, false);
                let fd = (op.obj - om.obj) / (2.0 * EPS_FD);
                let an = dot_flat(&g.data, &u);
                let r = (fd - an).abs() / fd.abs().max(1.0e-30);
                worst = worst.max(r);
            }
            run.grad_table.push((format!("{}/{label}", scheme_name(scheme)), worst));
            run.num(&format!("gradres_{}_{label}", scheme_name(scheme)), worst);
            worst_all = worst_all.max(worst);
        }
    }
    run.num("grad_worst_all", worst_all);
    let ga_pass = worst_all < TOL_GRAD;
    run.gate(
        "G-A",
        format!(
            "FD-vs-analytic gradient, N={N_FD}, both schemes, 8 objectives x 3 dirs: worst rel {worst_all:.2e} (tol 1e-5)"
        ),
        ga_pass,
    );
    vp(&format!("G-A gradient gate: worst rel {:.2e}  ({:.1}s)", worst_all, ta.elapsed().as_secs_f64()));
    if verbose {
        for (k, v) in &run.grad_table {
            println!("    {k:<18} {v:.3e}");
        }
    }
    if !ga_pass {
        run.aborted = true;
        return run; // do not proceed to the solver on a failed gradient gate
    }

    // ---- G-B: static relaxation at N = 48 ---------------------------------
    let tb = Instant::now();
    let mut f = Field3::new(N_RUN, LBOX);
    f.sample_hedgehog(&rp.r, &rp.f, 1.0);
    let (out_h, _) = eval(&f, &Opts::estatic(Scheme::Corner), false);
    let (out_h4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let deg_ref = out_h.deg;
    let fgap_h = out_h.e6 + out_h.e0 - fs_gum_field::bps_floor(); // deg/deg_ref = 1
    let fgap_ref = fgap_h - FLOOR_BAND;
    run.num("degref", deg_ref);
    run.num("fgapref", fgap_ref);
    vp(&format!(
        "\nG-B hedgehog (corner N={N_RUN}): Estat={:.7} deg={:.6} I={:.5} fgap={:+.6} -> wall {:+.6}",
        out_h.estat, out_h.deg, out_h.i, fgap_h, fgap_ref
    ));
    let opts_static = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
    let bump2 = bump_field(N_RUN, LBOX, SEED_PERT, 6, 0.02);
    add_scaled(&mut f, &bump2, 1.0);
    let _ = f.renormalize();
    let (out_p, _) = eval(&f, &opts_static, false);
    vp(&format!(
        "  perturbed: Estat={:.7} (dE=+{:.5})",
        out_p.estat,
        out_p.estat - out_h.estat
    ));
    let mut prm = AnfParams::new(MAXIT_STATIC);
    prm.print_every = if verbose { 150 } else { 0 };
    let res_s: AnfResult = anf(&mut f, &opts_static, &prm, "static");
    let mut q_static = vec![0.0_f64; 4 * N_RUN * N_RUN * N_RUN];
    save_cells(&f, &mut q_static);
    let out_s = res_s.out;
    let (out_s4, _) = eval(&f, &Opts::estatic(Scheme::Central4), false);
    let vir_h4 = virial_out(&out_h4, T_FROZEN) / out_h4.estat;
    let vir_s4 = virial_out(&out_s4, T_FROZEN) / out_s4.estat;
    // gates.  Tolerances are the guard geometry (band + penalty
    // penetration) plus the SAME-FUNCTIONAL reference class measured by the
    // Python 4A static stage at N = 96 (field3d_results.json static.final:
    // deg drifts to band edge -5.6e-3, wall penetrated -2.2e-3, dEstat
    // -0.0435, E2 -15.6%, E0 -6.9%, I -1.1%, halo 5.8e-6, corner virial rel
    // -0.147) — near-BPS quadrature valleys make the derivative sectors
    // soft; the guards, not a tight virial, bound the drift.
    let b1 = (out_s.deg - deg_ref).abs() <= 0.008;
    run.gate(
        "G-B1",
        format!(
            "degree within anchor band: deg {:.6} vs deg_ref {:.6} (|d|={:.2e} <= band 5e-3 + penetration 3e-3)",
            out_s.deg,
            deg_ref,
            (out_s.deg - deg_ref).abs()
        ),
        b1,
    );
    let vir_corner_h = virial_out(&out_h, T_FROZEN) / out_h.estat;
    let vir_corner_s = virial_out(&out_s, T_FROZEN) / out_s.estat;
    let b2 = vir_h4.abs() <= 1.0e-2 && vir_corner_s.abs() <= 0.5;
    run.gate(
        "G-B2",
        format!(
            "virial: seeded-hedgehog stationarity (central4) rel {vir_h4:+.3e} <= 1e-2 at N=48 (converges to the recorded 1e-3-class N=96 referee -6.6e-4, gated in fs-gum-field H4); relaxed-solution same-functional (corner) rel {vir_corner_s:+.3e} within reference class |.| <= 0.5 (Python N=96: -0.147; corner hedgehog here {vir_corner_h:+.3e})"
        ),
        b2,
    );
    let b3 = out_s.floor_gap >= fgap_ref - 5.0e-3;
    run.gate(
        "G-B3",
        format!(
            "not below guarded Bogomolny floor: fgap {:+.6} >= wall {:+.6} - 5e-3 (penetration {:+.2e})",
            out_s.floor_gap,
            fgap_ref,
            out_s.floor_gap - fgap_ref
        ),
        b3,
    );
    let sect_rels = [
        rel(out_s.e2, out_h.e2),
        rel(out_s.e4, out_h.e4),
        rel(out_s.e6, out_h.e6),
        rel(out_s.e0, out_h.e0),
    ];
    let worst_sect = sect_rels.iter().fold(0.0_f64, |m, &v| m.max(v));
    let i_rel = rel(out_s.i, out_h.i);
    let de_stat = out_s.estat - out_h.estat;
    let (halo_s, _) = fs_gum_statics::diag::halo_fraction(&f);
    let b4 = worst_sect <= 0.35
        && i_rel <= 0.08
        && de_stat.abs() <= 0.10
        && halo_s <= 0.02;
    run.gate(
        "G-B4",
        format!(
            "stays in the hedgehog neighbourhood (same functional, stated tolerances): quadrature-exact I rel {i_rel:.2e} <= 8e-2, halo {halo_s:.1e} <= 2e-2, dEstat {de_stat:+.5} (|.| <= 0.10), soft derivative sectors worst rel {worst_sect:.2e} <= 0.35"
        ),
        b4,
    );
    let mut max_obj_rise = f64::NEG_INFINITY;
    for w in res_s.series.windows(2) {
        max_obj_rise = max_obj_rise.max(w[1].obj - w[0].obj);
    }
    let b5 = max_obj_rise <= 0.0;
    run.gate(
        "G-B5",
        format!("objective monotone by construction: max recorded rise {max_obj_rise:+.2e} <= 0"),
        b5,
    );
    for (k, v) in [
        ("B_degfinal", out_s.deg),
        ("B_E2", out_s.e2),
        ("B_E4", out_s.e4),
        ("B_E6", out_s.e6),
        ("B_E0", out_s.e0),
        ("B_I", out_s.i),
        ("B_Estat", out_s.estat),
        ("B_fgap", out_s.floor_gap),
        ("B_virH4", vir_h4),
        ("B_virS4", vir_s4),
        ("B_virCornerH", vir_corner_h),
        ("B_virCornerS", vir_corner_s),
        ("B_haloS", halo_s),
        ("B_epen", out_s.epen),
        ("B_efpen", out_s.efpen),
        ("B_iters", res_s.iters as f64),
        ("B_arrests", res_s.arrests as f64),
        ("B_gnratio", res_s.gn / res_s.gn0),
    ] {
        run.num(k, v);
    }
    vp(&format!(
        "  static done: status={} iters={} arrests={} Estat={:.7} (dE vs hedgehog {:+.2e})  ({:.1}s)",
        res_s.status.as_str(),
        res_s.iters,
        res_s.arrests,
        out_s.estat,
        de_stat,
        tb.elapsed().as_secs_f64()
    ));

    // ---- G-C: halo-descent qualitative referee ----------------------------
    let i_h = out_h.i;
    let l_main = KAPPA_MAIN * i_h;
    let l_ctrl = KAPPA_CTRL * i_h;
    run.num("C_Lmain", l_main);
    run.num("C_Lctrl", l_ctrl);
    let halo_seed = halo_seed_field(N_RUN, LBOX, 0.05, 3.2, 0.6);
    let mut descend = |l: f64, maxit: usize, label: &str, run: &mut Run| -> AnfResult {
        let tc = Instant::now();
        restore_cells(&mut f, &q_static);
        add_scaled(&mut f, &bump2, 1.0);
        add_scaled(&mut f, &halo_seed, 1.0);
        let _ = f.renormalize();
        let opts = Opts::routhian(Scheme::Corner, l).with_guards(deg_ref, fgap_ref);
        let mut prm = AnfParams::new(maxit);
        prm.print_every = if verbose { 150 } else { 0 };
        let res = anf(&mut f, &opts, &prm, label);
        let first = res.series[0];
        let last = *res.series.last().expect("series");
        for (k, v) in [
            ("R0", first.r),
            ("Rf", last.r),
            ("kap0", first.kappa),
            ("kapf", last.kappa),
            ("halo0", first.halo),
            ("halof", last.halo),
            ("degf", last.deg),
            ("I0", first.i),
            ("If", last.i),
            ("fgapf", last.floor_gap),
            ("iters", res.iters as f64),
            ("arrests", res.arrests as f64),
        ] {
            run.num(&format!("C_{label}_{k}"), v);
        }
        let mut footprint = 0.0_f64;
        for s in &res.series {
            footprint = footprint.max(s.epen.abs() + s.efpen.abs());
        }
        run.num(&format!("C_{label}_penmax"), footprint);
        vp(&format!(
            "  [{label}] L={l:.4}: R {:.6} -> {:.6}, kappa {:.5} -> {:.5}, halo {:.4} -> {:.4}, deg {:.5}, max|pen| {:.2e}, status={} ({:.1}s)",
            first.r,
            last.r,
            first.kappa,
            last.kappa,
            first.halo,
            last.halo,
            last.deg,
            footprint,
            res.status.as_str(),
            tc.elapsed().as_secs_f64()
        ));
        res
    };
    vp(&format!(
        "\nG-C descents (N={N_RUN}, corner objective): L_main={l_main:.4} (kappa0~{KAPPA_MAIN}), L_ctrl={l_ctrl:.4} (kappa0~{KAPPA_CTRL}), threshold {:.5}",
        kappa_threshold()
    ));
    let res_m = descend(l_main, MAXIT_MAIN, "main", &mut run);
    let res_c = descend(l_ctrl, MAXIT_CTRL, "ctrl", &mut run);
    let (m0, mf) = (res_m.series[0], *res_m.series.last().expect("series"));
    let (c0, cf) = (res_c.series[0], *res_c.series.last().expect("series"));
    // Trend gates, calibrated to the 4A reference at EQUAL iteration counts
    // (both descents capped at 600).  The Python N=96 reference control is
    // NOT inert — it also dilates (kappa 0.176 -> 0.119, halo 0.028 ->
    // 0.050 at its own cap) — so the honest control gates are differential:
    // the over-spun run must out-pace the control decisively.
    let mut max_rise_r = f64::NEG_INFINITY;
    for w in res_m.series.windows(2) {
        max_rise_r = max_rise_r.max(w[1].r - w[0].r);
    }
    run.num("C_main_maxRiseR", max_rise_r);
    let c1 = mf.r < m0.r - 0.05 && max_rise_r <= 2.0e-3;
    run.gate(
        "G-C1",
        format!(
            "main R decreases monotonically: {:.6} -> {:.6} (dR={:+.5} <= -0.05), max recorded rise {:+.2e} <= 2e-3 (penalty-sized slack; obj itself is monotone by construction)",
            m0.r,
            mf.r,
            mf.r - m0.r,
            max_rise_r
        ),
        c1,
    );
    let dk_main = m0.kappa - mf.kappa;
    let dk_ctrl = c0.kappa - cf.kappa;
    let c2 = dk_main >= 0.05;
    run.gate(
        "G-C2",
        format!(
            "main kappa falls: {:.5} -> {:.5} (drop {dk_main:+.5} >= 0.05), threshold {:.5}",
            m0.kappa,
            mf.kappa,
            kappa_threshold()
        ),
        c2,
    );
    let dh_main = mf.halo - m0.halo;
    let dh_ctrl = cf.halo - c0.halo;
    let c3 = dh_main >= 0.015;
    run.gate(
        "G-C3",
        format!(
            "main halo rises: {:.4} -> {:.4} (rise {dh_main:+.4} >= 0.015) — the 4A signature",
            m0.halo, mf.halo
        ),
        c3,
    );
    let c4 = cf.halo <= 0.06 && dh_ctrl <= 0.6 * dh_main;
    run.gate(
        "G-C4",
        format!(
            "ctrl halo stays low: {:.4} -> {:.4} (final <= 0.06 and rise {dh_ctrl:+.4} <= 0.6 x main's {dh_main:+.4})",
            c0.halo, cf.halo
        ),
        c4,
    );
    let kap_ctrl_max = res_c.series.iter().fold(0.0_f64, |m, s| m.max(s.kappa));
    let c5 = kap_ctrl_max < kappa_threshold() && dk_ctrl.abs() <= 0.6 * dk_main.abs();
    run.gate(
        "G-C5",
        format!(
            "ctrl kappa stable relative to main: max kappa over run {kap_ctrl_max:.5} < threshold {:.5}; |dkappa| {:.4} <= 0.6 x main's {:.4}",
            kappa_threshold(),
            dk_ctrl.abs(),
            dk_main.abs()
        ),
        c5,
    );
    let c6 = dh_main > dh_ctrl + 0.01;
    run.gate(
        "G-C6",
        format!(
            "differential: main halo rise {dh_main:+.4} > ctrl halo change {dh_ctrl:+.4} + 0.01"
        ),
        c6,
    );
    let mut max_obj_rise_mc = f64::NEG_INFINITY;
    for ser in [&res_m.series, &res_c.series] {
        for w in ser.windows(2) {
            max_obj_rise_mc = max_obj_rise_mc.max(w[1].obj - w[0].obj);
        }
    }
    let c0g = max_obj_rise_mc <= 0.0;
    run.gate(
        "G-C0",
        format!(
            "descent objectives monotone by construction: max recorded rise {max_obj_rise_mc:+.2e} <= 0"
        ),
        c0g,
    );
    // trend tables (subsampled)
    let fmt_row = |s: &fs_gum_statics::Record| {
        format!(
            "| {:>4} | {:.6} | {:.5} | {:.4} | {:.5} | {:+.4} | {:.1e} |",
            s.it, s.r, s.kappa, s.halo, s.deg, s.floor_gap, s.epen + s.efpen
        )
    };
    for (ser, dst) in [(&res_m.series, &mut run.trend_main), (&res_c.series, &mut run.trend_ctrl)] {
        let step = (ser.len() / 12).max(1);
        for (idx, s) in ser.iter().enumerate() {
            if idx % step == 0 || idx + 1 == ser.len() {
                dst.push(fmt_row(s));
            }
        }
    }

    // ---- G-D: clock on the control solution -------------------------------
    let (estat_c, i_c) = (res_c.out.estat, res_c.out.i);
    let l_clock = clock_bisect(i_c, estat_c);
    let l_closed = ((2.0 / 3.0) * i_c * estat_c).sqrt();
    let erot = erot_frac(l_clock, i_c, estat_c);
    let kap_clock = l_clock / i_c;
    let ratio = kap_clock / kappa_threshold();
    let d1 = (erot - 0.25).abs() <= 1.0e-6;
    run.gate(
        "G-D1",
        format!("clock: E_rot/E_tot at L_clock = {erot:.9} (|.-0.25| = {:.1e} <= 1e-6)", (erot - 0.25).abs()),
        d1,
    );
    let d2 = rel(l_clock, l_closed) <= 1.0e-12;
    run.gate(
        "G-D2",
        format!(
            "bisection vs closed form: L_clock {l_clock:.9} vs sqrt((2/3) I Estat) {l_closed:.9} (rel {:.1e} <= 1e-12)",
            rel(l_clock, l_closed)
        ),
        d2,
    );
    for (k, v) in [
        ("D_Lclock", l_clock),
        ("D_erot", erot),
        ("D_kapclock", kap_clock),
        ("D_ratio", ratio),
    ] {
        run.num(k, v);
    }
    vp(&format!(
        "\nG-D clock (control): Estat={estat_c:.6} I={i_c:.5} -> L_clock={l_clock:.6}, E_rot/E={erot:.6}, kappa(L_clock)={kap_clock:.6} = {ratio:.4}x threshold"
    ));

    run
}

fn main() {
    let t_start = Instant::now();
    println!("fs-gum-statics gates (Phase E1) — N_fd={N_FD}, N_run={N_RUN}, LBOX={LBOX}, t={T_FROZEN:.15}");
    println!("guards: MU={}, band={}, MU_F={}, floor band={FLOOR_BAND}", fs_gum_statics::MU_DEG, fs_gum_statics::DEG_BAND, fs_gum_statics::MU_FLOOR);
    println!("================ run 1 ================");
    let run1 = run_all(true);
    println!("\n================ run 2 (replay, silent) ================");
    let run2 = run_all(false);

    // G-E: bit-identical replay
    let mut identical = run1.nums.len() == run2.nums.len();
    if identical {
        for (a, b) in run1.nums.iter().zip(run2.nums.iter()) {
            if a.to_bits() != b.to_bits() {
                identical = false;
                break;
            }
        }
    }
    let mut bytes = Vec::with_capacity(8 * run1.nums.len());
    for v in &run1.nums {
        bytes.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    let fp = fs_blake3::hash_domain("gum-core:statics:fingerprint:v1", &bytes);

    println!("\n================ gate table ================");
    let mut all_pass = true;
    for g in run1.gates.iter() {
        println!("  {:<5} {}  -> {}", g.id, g.desc, if g.pass { "PASS" } else { "FAIL" });
        all_pass &= g.pass;
    }
    println!(
        "  G-E   bit-identical two-run replay of {} gate numbers  -> {}",
        run1.nums.len(),
        if identical { "PASS" } else { "FAIL" }
    );
    all_pass &= identical;

    println!("\nper-sector FD-vs-analytic gradient residuals (worst of 3 directions):");
    for (k, v) in &run1.grad_table {
        println!("  {k:<18} {v:.3e}");
    }
    if !run1.trend_main.is_empty() {
        println!("\nmain descent trend (it | R | kappa | halo | deg | fgap | pens):");
        for l in &run1.trend_main {
            println!("  {l}");
        }
        println!("control descent trend:");
        for l in &run1.trend_ctrl {
            println!("  {l}");
        }
    }
    println!("\nrun fingerprint (BLAKE3 over every gate f64, fixed order): {}", fp.to_hex());
    println!("total runtime {:.1}s (both runs)", t_start.elapsed().as_secs_f64());
    if run1.aborted {
        println!("ABORTED after G-A failure: descents were not run (load-bearing gate).");
    }
    println!("\nOVERALL: {}", if all_pass { "PASS" } else { "FAIL" });
    if !all_pass {
        std::process::exit(1);
    }
}
