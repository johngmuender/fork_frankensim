//! GUM core Phase G5a — the F-R5 halo-descent referee at N = 192 on
//! fs-gum-kern's threaded sweeps (ROADMAP_v4_SCALE Phase G5).
//!
//! This binary re-runs the frozen fs-gum-statics G-C protocol (Phase E1,
//! itself the Rust port of the tier4-field 4A protocol) at larger grids
//! through `fs_gum_kern::{eval, anf}` — the campaign's F-R5 referee at up
//! to 8x the 4A cells and 64x the E1 gate-suite cells.  Everything is the
//! FROZEN protocol: hedgehog seed from the Step-1 radial solve at
//! t = T_FROZEN, LCG bump seed 424242 (K = 6, amp 0.02), tilt-halo shell
//! (eta = 0.05, r0 = 3.2, w = 0.6), guards MU = 5000 / band 0.005 /
//! MU_F = 400 / floor band 0.01 referenced to the engine's own hedgehog,
//! L_main = 0.266 I_hedgehog (over-spun), L_ctrl = 0.12 I_hedgehog
//! (sub-threshold), corner objective, arrested Newton flow.  Only the
//! grid (N, LBOX) and the iteration caps vary — both are command-line
//! arguments so the runs and their budgets are explicit in the launch
//! commands recorded in n192_RESULTS.md.
//!
//! Subcommands:
//!
//!   budget <N> <LBOX> <iters> [threads]
//!       Timing probe: seeded main-descent configuration, `iters` ANF
//!       iterations, prints s/eval — run BEFORE committing to caps.
//!
//!   proto <out.json> <N> <LBOX> <cap_static> <cap_main> <cap_ctrl> [threads]
//!       The full protocol: hedgehog refs -> R1 static relaxation ->
//!       R2 over-spun main descent -> R3 sub-threshold control
//!       (cap_ctrl = 0 skips R3) -> R4 clock bisection on the control
//!       endpoint.  Full instrumented series + measured-fact gates to
//!       `out.json`.
//!
//!   bitcheck <N> <LBOX> <iters> <threads_a> <threads_b>
//!       Thread bit-identity spot check at scale: the K-D construction
//!       (BLAKE3 over the full instrumented series, run metadata and the
//!       final field bytes of a short seeded main descent) at two thread
//!       counts; prints both hashes and PASS/FAIL on equality.
//!
//! Epistemic notice (binding, inherited): everything here is a
//! within-model computation on a speculative theory's functional.  The
//! descents replicate and refine the campaign's F-R5 mechanism as grid
//! measurements; PASS lines are measured facts against the frozen E1
//! thresholds, never adjudication and never statements about nature.

use fs_gum_field::radial::{radial_solve, RadialProfile};
use fs_gum_field::{bps_floor, Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{
    add_scaled, bump_field, clock_bisect, erot_frac, halo_seed_field, restore_cells, save_cells,
};
use fs_gum_statics::{kappa_threshold, AnfParams, AnfResult, Opts, Out, Record, FLOOR_BAND};

use std::fmt::Write as _;
use std::time::Instant;

const SEED_PERT: u64 = 424242; // the campaign perturbation seed
const BUMP_K: usize = 6;
const BUMP_AMP: f64 = 0.02;
const HALO_ETA: f64 = 0.05;
const HALO_R0: f64 = 3.2;
const HALO_W: f64 = 0.6;
const KAPPA_MAIN: f64 = 0.266; // over-spun target kappa_0 (E1 G-C)
const KAPPA_CTRL: f64 = 0.12; // control target kappa_0 (E1 G-C)
const RADIAL_N: usize = 4000;
const RADIAL_RMAX: f64 = 6.0;
const DEFAULT_THREADS: usize = 4;

fn hedgehog(rp: &RadialProfile, n: usize, lbox: f64) -> Field3 {
    let mut f = Field3::new(n, lbox);
    f.sample_hedgehog(&rp.r, &rp.f, 1.0);
    f
}

/// The seeded over-spun start configuration used by `budget`/`bitcheck`:
/// hedgehog + bump + halo shell (no static stage — the K-D construction).
fn seeded_main(rp: &RadialProfile, n: usize, lbox: f64) -> (Field3, Opts) {
    let mut f = hedgehog(rp, n, lbox);
    let (out_h, _) = fs_gum_kern::eval(&f, &Opts::estatic(Scheme::Corner), false, DEFAULT_THREADS);
    let deg_ref = out_h.deg;
    let fgap_ref = out_h.e6 + out_h.e0 - bps_floor() - FLOOR_BAND;
    let bump = bump_field(n, lbox, SEED_PERT, BUMP_K, BUMP_AMP);
    let halo = halo_seed_field(n, lbox, HALO_ETA, HALO_R0, HALO_W);
    add_scaled(&mut f, &bump, 1.0);
    add_scaled(&mut f, &halo, 1.0);
    let _ = f.renormalize();
    let o = Opts::routhian(Scheme::Corner, KAPPA_MAIN * out_h.i).with_guards(deg_ref, fgap_ref);
    (f, o)
}

// ---------------------------------------------------------------------------
// JSON (hand-rolled: fixed keys, f64s at full round-trip precision)
// ---------------------------------------------------------------------------

fn jf(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.17e}")
    } else {
        "null".to_string()
    }
}

fn json_record(r: &Record) -> String {
    format!(
        "{{\"it\":{},\"r\":{},\"obj\":{},\"epen\":{},\"efpen\":{},\"floor_gap\":{},\"estat\":{},\"e2\":{},\"e4\":{},\"e6\":{},\"e0\":{},\"i\":{},\"kappa\":{},\"deg\":{},\"halo\":{},\"gnorm\":{},\"dt\":{},\"arrests\":{}}}",
        r.it,
        jf(r.r),
        jf(r.obj),
        jf(r.epen),
        jf(r.efpen),
        jf(r.floor_gap),
        jf(r.estat),
        jf(r.e2),
        jf(r.e4),
        jf(r.e6),
        jf(r.e0),
        jf(r.i),
        jf(r.kappa),
        jf(r.deg),
        jf(r.halo),
        jf(r.gnorm),
        jf(r.dt),
        r.arrests
    )
}

fn json_run(res: &AnfResult, l: Option<f64>, wall: f64) -> String {
    let mut s = String::new();
    let _ = write!(
        s,
        "{{\"l\":{},\"status\":\"{}\",\"iters\":{},\"arrests\":{},\"gn0\":{},\"gn\":{},\"wall_s\":{},\"series\":[",
        l.map_or("null".to_string(), jf),
        res.status.as_str(),
        res.iters,
        res.arrests,
        jf(res.gn0),
        jf(res.gn),
        jf(wall)
    );
    for (i, r) in res.series.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&json_record(r));
    }
    s.push_str("]}");
    s
}

// ---------------------------------------------------------------------------
// measured-fact gates (the frozen E1 thresholds, restated verbatim)
// ---------------------------------------------------------------------------

struct Gate {
    id: &'static str,
    desc: String,
    pass: bool,
}

fn crossing(series: &[Record], thr: f64) -> Option<(usize, usize)> {
    for w in series.windows(2) {
        if w[0].kappa >= thr && w[1].kappa < thr {
            return Some((w[0].it, w[1].it));
        }
    }
    None
}

fn max_rise(series: &[Record], get: impl Fn(&Record) -> f64) -> f64 {
    let mut m = f64::NEG_INFINITY;
    for w in series.windows(2) {
        m = m.max(get(&w[1]) - get(&w[0]));
    }
    m
}

// ---------------------------------------------------------------------------
// proto
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_lines)]
fn proto(
    out_path: &str,
    n: usize,
    lbox: f64,
    cap_static: usize,
    cap_main: usize,
    cap_ctrl: usize,
    threads: usize,
) {
    let t_all = Instant::now();
    let thr = kappa_threshold();
    println!(
        "G5a proto: N={n} LBOX={lbox} h={:.9} caps static/main/ctrl {cap_static}/{cap_main}/{cap_ctrl} threads={threads}",
        2.0 * lbox / (n as f64)
    );
    println!("bit contract: {}", fs_gum_kern::GUM_KERN_BIT_SEMANTICS);
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);

    // ---- hedgehog references (the engine's own referee) --------------------
    let t0 = Instant::now();
    let fh = hedgehog(&rp, n, lbox);
    let (out_h, _) = fs_gum_kern::eval(&fh, &Opts::estatic(Scheme::Corner), false, threads);
    let deg_ref = out_h.deg;
    let fgap_h = out_h.e6 + out_h.e0 - bps_floor(); // deg/deg_ref = 1 on the reference itself
    let fgap_ref = fgap_h - FLOOR_BAND;
    let (halo_h, _) = fs_gum_statics::diag::halo_fraction(&fh);
    println!(
        "hedgehog (corner): Estat={:.7} deg={:.6} I={:.5} fgap={:+.6} -> wall {:+.6} halo={:.2e}  ({:.1}s)",
        out_h.estat,
        out_h.deg,
        out_h.i,
        fgap_h,
        fgap_ref,
        halo_h,
        t0.elapsed().as_secs_f64()
    );

    let bump = bump_field(n, lbox, SEED_PERT, BUMP_K, BUMP_AMP);
    let halo_seed = halo_seed_field(n, lbox, HALO_ETA, HALO_R0, HALO_W);
    let mut gates: Vec<Gate> = Vec::new();

    // ---- R1: static relaxation ---------------------------------------------
    let opts_static = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
    let mut f = fs_gum_statics::diag::clone_field(&fh);
    add_scaled(&mut f, &bump, 1.0);
    let _ = f.renormalize();
    let mut prm = AnfParams::new(cap_static);
    prm.print_every = 100;
    let ts = Instant::now();
    let res_s = fs_gum_kern::anf(&mut f, &opts_static, &prm, "static", threads);
    let wall_s = ts.elapsed().as_secs_f64();
    let out_s = res_s.out;
    let (halo_s, _) = fs_gum_statics::diag::halo_fraction(&f);
    let mut q_static = vec![0.0_f64; 4 * n * n * n];
    save_cells(&f, &mut q_static);

    let rel = |a: f64, b: f64| (a / b - 1.0).abs();
    gates.push(Gate {
        id: "R1-deg",
        desc: format!(
            "static degree within anchor band: deg {:.6} vs deg_ref {:.6} (|d|={:.2e} <= 0.008)",
            out_s.deg,
            deg_ref,
            (out_s.deg - deg_ref).abs()
        ),
        pass: (out_s.deg - deg_ref).abs() <= 0.008,
    });
    gates.push(Gate {
        id: "R1-floor",
        desc: format!(
            "floor held: fgap {:+.6} >= wall {:+.6} - 5e-3 (penetration {:+.2e})",
            out_s.floor_gap,
            fgap_ref,
            out_s.floor_gap - fgap_ref
        ),
        pass: out_s.floor_gap >= fgap_ref - 5.0e-3,
    });
    let sect_worst = [
        rel(out_s.e2, out_h.e2),
        rel(out_s.e4, out_h.e4),
        rel(out_s.e6, out_h.e6),
        rel(out_s.e0, out_h.e0),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max);
    let i_rel = rel(out_s.i, out_h.i);
    let de_stat = out_s.estat - out_h.estat;
    gates.push(Gate {
        id: "R1-hood",
        desc: format!(
            "sector neighbourhood vs seeded hedgehog: I rel {i_rel:.2e} <= 8e-2, halo {halo_s:.1e} <= 2e-2, dEstat {de_stat:+.5} (|.| <= 0.10), soft sectors worst rel {sect_worst:.2e} <= 0.35"
        ),
        pass: sect_worst <= 0.35 && i_rel <= 0.08 && de_stat.abs() <= 0.10 && halo_s <= 0.02,
    });
    gates.push(Gate {
        id: "R1-mono",
        desc: format!(
            "static objective monotone: max recorded rise {:+.2e} <= 0",
            max_rise(&res_s.series, |r| r.obj)
        ),
        pass: max_rise(&res_s.series, |r| r.obj) <= 0.0,
    });
    println!(
        "R1 static done: status={} iters={} arrests={} Estat={:.7} (dE {de_stat:+.4}) deg={:.6} ({:.1}s)",
        res_s.status.as_str(),
        res_s.iters,
        res_s.arrests,
        out_s.estat,
        out_s.deg,
        wall_s
    );

    // ---- R2: over-spun main descent ----------------------------------------
    let l_main = KAPPA_MAIN * out_h.i;
    restore_cells(&mut f, &q_static);
    add_scaled(&mut f, &bump, 1.0);
    add_scaled(&mut f, &halo_seed, 1.0);
    let _ = f.renormalize();
    let opts_main = Opts::routhian(Scheme::Corner, l_main).with_guards(deg_ref, fgap_ref);
    let mut prm = AnfParams::new(cap_main);
    prm.print_every = 100;
    let tm = Instant::now();
    let res_m = fs_gum_kern::anf(&mut f, &opts_main, &prm, "main", threads);
    let wall_m = tm.elapsed().as_secs_f64();
    let (m0, mf) = (res_m.series[0], *res_m.series.last().expect("series"));
    let cross_m = crossing(&res_m.series, thr);
    let rise_r = max_rise(&res_m.series, |r| r.r);
    gates.push(Gate {
        id: "R2-desc",
        desc: format!(
            "main R decreases monotonically: {:.6} -> {:.6} (dR {:+.5} <= -0.05), max recorded rise {rise_r:+.2e} <= 2e-3",
            m0.r,
            mf.r,
            mf.r - m0.r
        ),
        pass: mf.r < m0.r - 0.05 && rise_r <= 2.0e-3,
    });
    gates.push(Gate {
        id: "R2-kappa",
        desc: format!(
            "main kappa falls through threshold {thr:.5}: {:.5} -> {:.5} (drop {:+.5} >= 0.05); crossing {}",
            m0.kappa,
            mf.kappa,
            m0.kappa - mf.kappa,
            cross_m.map_or("NONE".to_string(), |(a, b)| format!("between it {a} and {b}")),
        ),
        pass: m0.kappa - mf.kappa >= 0.05 && cross_m.is_some(),
    });
    gates.push(Gate {
        id: "R2-halo",
        desc: format!(
            "main halo rises (the 4A signature): {:.4} -> {:.4} (rise {:+.4} >= 0.015)",
            m0.halo,
            mf.halo,
            mf.halo - m0.halo
        ),
        pass: mf.halo - m0.halo >= 0.015,
    });
    gates.push(Gate {
        id: "R2-mono",
        desc: format!(
            "main objective monotone: max recorded rise {:+.2e} <= 0",
            max_rise(&res_m.series, |r| r.obj)
        ),
        pass: max_rise(&res_m.series, |r| r.obj) <= 0.0,
    });
    println!(
        "R2 main done: L={l_main:.4} kappa {:.5} -> {:.5} (threshold crossing {}), halo {:.4} -> {:.4}, R {:.6} -> {:.6} ({:.1}s)",
        m0.kappa,
        mf.kappa,
        cross_m.map_or("NONE".to_string(), |(a, b)| format!("it {a}..{b}")),
        m0.halo,
        mf.halo,
        m0.r,
        mf.r,
        wall_m
    );

    // ---- R3: sub-threshold control + R4 clock ------------------------------
    let mut ctrl_json = "null".to_string();
    let mut clock_json = "null".to_string();
    if cap_ctrl > 0 {
        let l_ctrl = KAPPA_CTRL * out_h.i;
        restore_cells(&mut f, &q_static);
        add_scaled(&mut f, &bump, 1.0);
        add_scaled(&mut f, &halo_seed, 1.0);
        let _ = f.renormalize();
        let opts_ctrl = Opts::routhian(Scheme::Corner, l_ctrl).with_guards(deg_ref, fgap_ref);
        let mut prm = AnfParams::new(cap_ctrl);
        prm.print_every = 100;
        let tc = Instant::now();
        let res_c = fs_gum_kern::anf(&mut f, &opts_ctrl, &prm, "ctrl", threads);
        let wall_c = tc.elapsed().as_secs_f64();
        let (c0, cf) = (res_c.series[0], *res_c.series.last().expect("series"));
        let kap_ctrl_max = res_c.series.iter().fold(0.0_f64, |m, s| m.max(s.kappa));
        let dh_main = mf.halo - m0.halo;
        let dh_ctrl = cf.halo - c0.halo;
        let dk_main = m0.kappa - mf.kappa;
        let dk_ctrl = c0.kappa - cf.kappa;
        gates.push(Gate {
            id: "R3-halo",
            desc: format!(
                "ctrl halo stays low: {:.4} -> {:.4} (final <= 0.06 and rise {dh_ctrl:+.4} <= 0.6 x main's {dh_main:+.4})",
                c0.halo, cf.halo
            ),
            pass: cf.halo <= 0.06 && dh_ctrl <= 0.6 * dh_main,
        });
        gates.push(Gate {
            id: "R3-kappa",
            desc: format!(
                "ctrl self-limiting: max kappa {kap_ctrl_max:.5} < threshold {thr:.5}; |dkappa| {:.4} <= 0.6 x main's {:.4}",
                dk_ctrl.abs(),
                dk_main.abs()
            ),
            pass: kap_ctrl_max < thr && dk_ctrl.abs() <= 0.6 * dk_main.abs(),
        });
        // tilt ceiling: I_final vs (3/2) I_hedgehog (reported, F-R4 bound)
        let tilt_ceiling = 1.5 * out_h.i;
        println!(
            "R3 ctrl done: L={l_ctrl:.4} kappa {:.5} -> {:.5}, halo {:.4} -> {:.4}, I {:.4} -> {:.4} (tilt ceiling 1.5 I_h = {tilt_ceiling:.4}) ({wall_c:.1}s)",
            c0.kappa, cf.kappa, c0.halo, cf.halo, c0.i, cf.i
        );
        // R4: the clock on the control endpoint (the E1 G-D / 4A Section-K referee)
        let (estat_c, i_c) = (res_c.out.estat, res_c.out.i);
        let l_clock = clock_bisect(i_c, estat_c);
        let l_closed = ((2.0 / 3.0) * i_c * estat_c).sqrt();
        let erot = erot_frac(l_clock, i_c, estat_c);
        let kap_clock = l_clock / i_c;
        let ratio = kap_clock / thr;
        gates.push(Gate {
            id: "R4-clock",
            desc: format!(
                "clock: E_rot/E_tot at L_clock = {erot:.9} (|.-0.25| = {:.1e} <= 1e-6); bisection vs closed form rel {:.1e} <= 1e-12",
                (erot - 0.25).abs(),
                (l_clock / l_closed - 1.0).abs()
            ),
            pass: (erot - 0.25).abs() <= 1.0e-6 && (l_clock / l_closed - 1.0).abs() <= 1.0e-12,
        });
        println!(
            "R4 clock (ctrl endpoint): Estat={estat_c:.6} I={i_c:.5} -> L_clock={l_clock:.6}, E_rot/E={erot:.9}, kappa(L_clock)={kap_clock:.6} = {ratio:.4}x threshold"
        );
        ctrl_json = json_run(&res_c, Some(l_ctrl), wall_c);
        clock_json = format!(
            "{{\"estat\":{},\"i\":{},\"l_clock\":{},\"erot\":{},\"kappa_clock\":{},\"ratio_to_threshold\":{}}}",
            jf(estat_c),
            jf(i_c),
            jf(l_clock),
            jf(erot),
            jf(kap_clock),
            jf(ratio)
        );
    }

    // ---- gate table + JSON --------------------------------------------------
    println!("\n---- measured-fact gate table (frozen E1 thresholds) ----");
    let mut all_pass = true;
    for g in &gates {
        println!("  {:<9} {}  -> {}", g.id, g.desc, if g.pass { "PASS" } else { "FAIL" });
        all_pass &= g.pass;
    }

    let mut j = String::new();
    let _ = write!(
        j,
        "{{\n\"phase\":\"G5a\",\"bit_contract\":\"{}\",\n\"n\":{n},\"lbox\":{},\"h\":{},\"threads\":{threads},\n\"protocol\":{{\"seed\":{SEED_PERT},\"bump_k\":{BUMP_K},\"bump_amp\":{BUMP_AMP},\"halo_seed\":[{HALO_ETA},{HALO_R0},{HALO_W}],\"kappa_main\":{KAPPA_MAIN},\"kappa_ctrl\":{KAPPA_CTRL},\"t_frozen\":{},\"caps\":[{cap_static},{cap_main},{cap_ctrl}]}},\n\"kappa_threshold\":{},\n\"hedgehog\":{{\"estat\":{},\"e2\":{},\"e4\":{},\"e6\":{},\"e0\":{},\"i\":{},\"deg\":{},\"fgap\":{},\"deg_ref\":{},\"fgap_ref\":{},\"halo\":{}}},\n",
        fs_gum_kern::GUM_KERN_BIT_SEMANTICS,
        jf(lbox),
        jf(2.0 * lbox / (n as f64)),
        jf(T_FROZEN),
        jf(thr),
        jf(out_h.estat),
        jf(out_h.e2),
        jf(out_h.e4),
        jf(out_h.e6),
        jf(out_h.e0),
        jf(out_h.i),
        jf(out_h.deg),
        jf(fgap_h),
        jf(deg_ref),
        jf(fgap_ref),
        jf(halo_h)
    );
    let _ = write!(
        j,
        "\"static\":{},\n\"main\":{},\n\"main_kappa_crossing\":{},\n\"ctrl\":{ctrl_json},\n\"clock\":{clock_json},\n\"gates\":[",
        json_run(&res_s, None, wall_s),
        json_run(&res_m, Some(l_main), wall_m),
        cross_m.map_or("null".to_string(), |(a, b)| format!("[{a},{b}]")),
    );
    for (i, g) in gates.iter().enumerate() {
        if i > 0 {
            j.push(',');
        }
        let _ = write!(
            j,
            "{{\"id\":\"{}\",\"desc\":\"{}\",\"pass\":{}}}",
            g.id,
            g.desc.replace('"', "'"),
            g.pass
        );
    }
    let _ = write!(
        j,
        "],\n\"all_pass\":{all_pass},\n\"total_wall_s\":{}\n}}\n",
        jf(t_all.elapsed().as_secs_f64())
    );
    std::fs::write(out_path, j).expect("write json");
    println!(
        "\nwrote {out_path}; total wall {:.1}s; OVERALL: {}",
        t_all.elapsed().as_secs_f64(),
        if all_pass { "PASS" } else { "FAIL" }
    );
    if !all_pass {
        std::process::exit(1);
    }
}

// ---------------------------------------------------------------------------
// budget + bitcheck
// ---------------------------------------------------------------------------

fn budget(n: usize, lbox: f64, iters: usize, threads: usize) {
    println!("G5a budget probe: N={n} LBOX={lbox} iters={iters} threads={threads}");
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);
    let t0 = Instant::now();
    let (mut f, o) = seeded_main(&rp, n, lbox);
    println!("  setup (hedgehog + refs + seeds) {:.1}s", t0.elapsed().as_secs_f64());
    let prm = AnfParams::new(iters);
    let t1 = Instant::now();
    let res = fs_gum_kern::anf(&mut f, &o, &prm, "budget", threads);
    let wall = t1.elapsed().as_secs_f64();
    // anf performs iters+1 gradient evals (initial + one per iteration)
    let per = wall / (iters as f64 + 1.0);
    println!(
        "  {iters} ANF iterations in {wall:.2}s -> {per:.3} s/eval (iters+1 evals incl. maps/norms; arrests {})",
        res.arrests
    );
    println!(
        "  projections: 100 it = {:.1} min; 600 it = {:.1} min; 900 it = {:.1} min",
        per * 101.0 / 60.0,
        per * 601.0 / 60.0,
        per * 901.0 / 60.0
    );
}

fn hash_descent(res: &AnfResult, f: &Field3) -> String {
    let mut h = fs_blake3::Blake3::new();
    h.update(fs_gum_kern::GUM_KERN_BIT_SEMANTICS.as_bytes());
    h.update(b"gum-kern:g5a-bitcheck:v1");
    for s in &res.series {
        for v in [
            s.r, s.obj, s.epen, s.efpen, s.floor_gap, s.estat, s.e2, s.e4, s.e6, s.e0, s.i,
            s.kappa, s.deg, s.halo, s.gnorm, s.dt,
        ] {
            h.update(&v.to_bits().to_le_bytes());
        }
    }
    for v in [res.gn0, res.gn] {
        h.update(&v.to_bits().to_le_bytes());
    }
    h.update(&(res.iters as u64).to_le_bytes());
    h.update(&(res.arrests as u64).to_le_bytes());
    for v in f.data() {
        h.update(&v.to_bits().to_le_bytes());
    }
    h.finalize().to_hex()
}

fn bitcheck(n: usize, lbox: f64, iters: usize, ta: usize, tb: usize) {
    println!("G5a thread bit-identity spot check: N={n} LBOX={lbox} iters={iters} threads {ta} vs {tb}");
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);
    let mut hashes = Vec::new();
    for threads in [ta, tb] {
        let (mut f, o) = seeded_main(&rp, n, lbox);
        let prm = AnfParams::new(iters);
        let t0 = Instant::now();
        let res = fs_gum_kern::anf(&mut f, &o, &prm, "bitcheck", threads);
        let hh = hash_descent(&res, &f);
        println!(
            "  threads={threads}: {hh}  ({:.1}s, {} arrests)",
            t0.elapsed().as_secs_f64(),
            res.arrests
        );
        hashes.push(hh);
    }
    let pass = hashes[0] == hashes[1];
    println!(
        "bit-identity ({iters}-iter seeded main descent at N={n}, series + endpoint field bytes): {}",
        if pass { "PASS" } else { "FAIL" }
    );
    if !pass {
        std::process::exit(1);
    }
}

// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let usage = "usage: n192_fr5 budget <N> <LBOX> <iters> [threads]\n       n192_fr5 proto <out.json> <N> <LBOX> <cap_static> <cap_main> <cap_ctrl> [threads]\n       n192_fr5 bitcheck <N> <LBOX> <iters> <threads_a> <threads_b>";
    let cmd = args.get(1).map(String::as_str).unwrap_or("");
    let p = |i: usize| -> usize { args[i].parse().expect("integer arg") };
    let pf = |i: usize| -> f64 { args[i].parse().expect("float arg") };
    match cmd {
        "budget" if args.len() >= 5 => {
            let threads = if args.len() > 5 { p(5) } else { DEFAULT_THREADS };
            budget(p(2), pf(3), p(4), threads);
        }
        "proto" if args.len() >= 8 => {
            let threads = if args.len() > 8 { p(8) } else { DEFAULT_THREADS };
            proto(&args[2], p(3), pf(4), p(5), p(6), p(7), threads);
        }
        "bitcheck" if args.len() >= 7 => {
            bitcheck(p(2), pf(3), p(4), p(5), p(6));
        }
        _ => {
            eprintln!("{usage}");
            std::process::exit(2);
        }
    }
}
