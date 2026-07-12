//! GUM core Phase G1 gate binary for `fs-gum-kern`.
//!
//! Gates (all numbers machine-checked):
//!
//!   K-A  ONE-TIME CROSS-CHECK (the golden-bump audit): the new tiled
//!        serial path vs the OLD fs-gum-statics engine on identical
//!        inputs (N = 48 hedgehog + frozen LCG bump, both schemes, bare
//!        and fully-guarded rotor objectives) at machine-eps-class
//!        tolerance.  The measured deltas ARE the deliberate golden bump
//!        (tiled/gather order vs old flat/scatter order) and are printed
//!        in full for RESULTS.md — the fs-la KC-contract precedent.
//!   K-B  FD-vs-analytic gradient gate re-run THROUGH the new layer
//!        (N = 24, both schemes, 8 objectives x 3 tangent directions,
//!        tol 1e-5) — the gather adjoint is exact, not just close to the
//!        old scatter.
//!   K-C  NEW bit goldens, frozen for the tiled order: BLAKE3 over the
//!        new-serial eval outputs (Out + gradient bits + norms), domain-
//!        tagged by GUM_KERN_BIT_SEMANTICS.
//!   K-D  THE HEADLINE GATE: serial vs threaded at 1/2/3/4 threads —
//!        BLAKE3 hash of all outputs (eval fwd+grad, renormalize field
//!        bytes + drift, norms, and an 80-iteration arrested-Newton
//!        descent: full instrumented series + final field bytes)
//!        bit-identical across thread counts, BY CONSTRUCTION.
//!   K-E  the fs-gum-statics G-C main halo-descent protocol at N = 48
//!        run through the new layer (static relax 900 + over-spun main
//!        600 iters): same qualitative referee (monotone R, kappa
//!        falling through threshold, halo growth), plus the endpoint
//!        compared against the OLD engine within stated tolerance bands
//!        (it need not bit-match: the arrest cascade amplifies the
//!        ~1e-12-class order change into diverging trajectories —
//!        documented, the sanctioned cross-backend golden policy).
//!   K-F  measured performance: ns/point for forward and gradient
//!        sweeps at N = 48/96/128, serial-old vs serial-new vs
//!        threaded-4; the speedup gate is >= 3.0x threaded-4 over
//!        serial-new at N = 96 on the corner gradient path (the ANF hot
//!        path).
//!   K-G  bit-identical two-run replay of the K-A..K-D sections (BLAKE3
//!        fingerprint printed).
//!
//! Epistemic notice (binding): every PASS certifies a WITHIN-MODEL
//! property of a kernel layer — determinism by construction, equivalence
//! to a frozen reference, measured speedup.  Nothing here says anything
//! about nature.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{Field3, Scheme, T_FROZEN};
use fs_gum_statics::diag::{
    add_scaled, bump_field, clone_field, dot_flat, halo_seed_field, normalize_l2, restore_cells,
    save_cells, tangent_project,
};
use fs_gum_statics::{AnfParams, AnfResult, Opts, Out, Record, FLOOR_BAND};

use std::hint::black_box;
use std::time::Instant;

const LBOX: f64 = 4.5;
const N_FD: usize = 24; // gradient-gate grid
const N_HASH: usize = 32; // thread-identity mini-descent grid
const N_RUN: usize = 48; // G-C protocol grid
const EPS_FD: f64 = 1.0e-5;
const TOL_GRAD: f64 = 1.0e-5;
const SEED_PERT: u64 = 424242;
const SEED_FD_BASE: u64 = 11;
const L_FD: f64 = 8.4979;
const KAPPA_MAIN: f64 = 0.266;
const MAXIT_STATIC: usize = 900;
const MAXIT_MAIN: usize = 600;
const MAXIT_HASH: usize = 80;
const THREADS: usize = 4;

/// K-A tolerance: machine-eps-class agreement of the tiled/gather order
/// with the old flat/scatter order (the measured deltas are the golden
/// bump; anything past this bound would be a real bug, not order noise).
const TOL_XCHECK: f64 = 5.0e-11;

/// K-C frozen goldens for the NEW canonical (tiled + gather) order,
/// domain-tagged by GUM_KERN_BIT_SEMANTICS.  Deliberately re-frozen in
/// this phase — the documented golden bump (see RESULTS.md).
const GOLDEN_EVAL: &str = "5a4cc174d094b0e878ee1f46c50f764202bc5cdc6827130a8d1c1a190d07fefa";

/// K-E endpoint tolerance bands (new-layer vs old-engine G-C main
/// endpoint after 900 + 600 iterations).  These are PHYSICS-SIZED bands
/// for two legitimately different fixed-order backends whose arrest
/// cascades diverge from a ~1e-12 seed difference; measured deltas are
/// printed alongside.  R/estat: absolute; kappa/halo/deg: absolute.
const BAND_R: f64 = 0.02;
const BAND_KAPPA: f64 = 0.02;
const BAND_HALO: f64 = 0.02;
const BAND_DEG: f64 = 0.005;

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

fn rel(a: f64, b: f64) -> f64 {
    (a / b - 1.0).abs()
}

fn out_fields(o: &Out) -> [f64; 12] {
    [o.e2, o.e4, o.e6, o.e0, o.i, o.deg, o.estat, o.r, o.floor_gap, o.epen, o.efpen, o.obj]
}

fn hash_stream(domain: &str, parts: &[&[f64]], extra: &[u64]) -> fs_blake3::ContentHash {
    let mut h = fs_blake3::Blake3::new();
    h.update(fs_gum_kern::GUM_KERN_BIT_SEMANTICS.as_bytes());
    h.update(domain.as_bytes());
    for p in parts {
        for v in *p {
            h.update(&v.to_bits().to_le_bytes());
        }
    }
    for v in extra {
        h.update(&v.to_le_bytes());
    }
    h.finalize()
}

fn hedgehog(rp: &fs_gum_field::radial::RadialProfile, n: usize) -> Field3 {
    let mut f = Field3::new(n, LBOX);
    f.sample_hedgehog(&rp.r, &rp.f, 1.0);
    f
}

/// The eight FD objectives of the statics G-A gate (copied verbatim).
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

/// Serialize an ANF result for hashing: every instrumented record plus
/// the final Out; iteration/arrest counts ride as u64s.
fn anf_stream(res: &AnfResult) -> (Vec<f64>, Vec<u64>) {
    let mut fs = Vec::new();
    for s in &res.series {
        let rec: &Record = s;
        fs.extend_from_slice(&[
            rec.r, rec.obj, rec.epen, rec.efpen, rec.floor_gap, rec.estat, rec.e2, rec.e4,
            rec.e6, rec.e0, rec.i, rec.kappa, rec.deg, rec.halo, rec.gnorm, rec.dt,
        ]);
    }
    fs.extend_from_slice(&out_fields(&res.out));
    fs.push(res.gn0);
    fs.push(res.gn);
    let us = vec![res.iters as u64, res.arrests as u64];
    (fs, us)
}

// ---------------------------------------------------------------------------
// K-A .. K-D (the replayed core)
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_lines)]
fn core(rp: &fs_gum_field::radial::RadialProfile, verbose: bool) -> Run {
    let mut run = Run::default();
    let vp = |s: &str| {
        if verbose {
            println!("{s}");
        }
    };

    // shared canonical input: N = 48 hedgehog + frozen bump, renormalised
    let mut f48 = hedgehog(rp, N_RUN);
    let bump = bump_field(N_RUN, LBOX, SEED_PERT, 6, 0.02);
    add_scaled(&mut f48, &bump, 1.0);
    let _ = f48.renormalize();

    // ---- K-A: one-time cross-check vs the OLD engine ----------------------
    let ta = Instant::now();
    let mut worst_fwd = 0.0_f64;
    let mut worst_grad = 0.0_f64;
    let mut worst_norm = 0.0_f64;
    for scheme in [Scheme::Central4, Scheme::Corner] {
        for (case, o) in [
            ("estatic", Opts::estatic(scheme)),
            ("rotor+guards", Opts::routhian(scheme, L_FD).with_guards(1.0, 10.0)),
        ] {
            let (oo, og) = fs_gum_statics::eval(&f48, &o, true);
            let (no, ng) = fs_gum_kern::eval(&f48, &o, true, 1);
            let mut wf = 0.0_f64;
            for (a, b) in out_fields(&oo).iter().zip(out_fields(&no).iter()) {
                if *b != 0.0 || *a != 0.0 {
                    wf = wf.max(rel(*a, *b));
                }
            }
            let (og, ng) = (og.expect("grad"), ng.expect("grad"));
            let gmax = og.data.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
            let dmax = og
                .data
                .iter()
                .zip(ng.data.iter())
                .fold(0.0_f64, |m, (&a, &b)| m.max((a - b).abs()));
            let wg = dmax / gmax;
            let nrel = rel(og.norm(), fs_gum_kern::norm_flat(&ng.data, N_RUN, 1));
            run.num(&format!("KA_fwd_{}_{case}", scheme_name(scheme)), wf);
            run.num(&format!("KA_grad_{}_{case}", scheme_name(scheme)), wg);
            worst_fwd = worst_fwd.max(wf);
            worst_grad = worst_grad.max(wg);
            worst_norm = worst_norm.max(nrel);
            if case == "estatic" {
                vp(&format!(
                    "  K-A {}: obj old {:.17e} | new {:.17e} (rel {:.2e}); grad max|d|/max|g| {:.2e}; gnorm rel {:.2e}",
                    scheme_name(scheme),
                    oo.obj,
                    no.obj,
                    rel(oo.obj, no.obj),
                    wg,
                    nrel
                ));
            }
        }
    }
    run.num("KA_worst_fwd", worst_fwd);
    run.num("KA_worst_grad", worst_grad);
    run.num("KA_worst_gnorm", worst_norm);
    run.gate(
        "K-A1",
        format!(
            "old-engine cross-check, forward sums (tiled vs flat order): worst rel {worst_fwd:.2e} (tol {TOL_XCHECK:.0e}) — the measured golden bump"
        ),
        worst_fwd < TOL_XCHECK,
    );
    run.gate(
        "K-A2",
        format!(
            "old-engine cross-check, gradient (gather vs scatter order): worst max|d|/max|g| {worst_grad:.2e}, gnorm rel {worst_norm:.2e} (tol {TOL_XCHECK:.0e})"
        ),
        worst_grad < TOL_XCHECK && worst_norm < TOL_XCHECK,
    );
    vp(&format!("K-A cross-check done ({:.1}s)", ta.elapsed().as_secs_f64()));

    // ---- K-B: FD-vs-analytic gradient through the new layer ----------------
    let tb = Instant::now();
    let mut fbase = hedgehog(rp, N_FD);
    let bump_fd = bump_field(N_FD, LBOX, SEED_FD_BASE, 6, 0.05);
    add_scaled(&mut fbase, &bump_fd, 1.0);
    let _ = fbase.renormalize();
    let mut worst_all = 0.0_f64;
    for (si, scheme) in [Scheme::Central4, Scheme::Corner].into_iter().enumerate() {
        for (ci, (label, c, l, guards)) in fd_cases().into_iter().enumerate() {
            let mut o = Opts { scheme, t: T_FROZEN, c, l, deg_ref: 1.0, anchor: None, wall: None };
            if guards {
                o = o.with_guards(1.0, 10.0);
            }
            let (_out, g) = fs_gum_kern::eval(&fbase, &o, true, THREADS);
            let g = g.expect("grad");
            let mut worst = 0.0_f64;
            for dir in 0..3u64 {
                let seed = 9000 + 100 * (si as u64) + 10 * (ci as u64) + dir;
                let mut u = bump_field(N_FD, LBOX, seed, 6, 1.0);
                tangent_project(&fbase, &mut u);
                let _ = normalize_l2(&mut u);
                let mut fp = clone_field(&fbase);
                add_scaled(&mut fp, &u, EPS_FD);
                let _ = fp.renormalize();
                let mut fm = clone_field(&fbase);
                add_scaled(&mut fm, &u, -EPS_FD);
                let _ = fm.renormalize();
                let (op, _) = fs_gum_kern::eval(&fp, &o, false, THREADS);
                let (om, _) = fs_gum_kern::eval(&fm, &o, false, THREADS);
                let fd = (op.obj - om.obj) / (2.0 * EPS_FD);
                let an = dot_flat(&g.data, &u);
                let r = (fd - an).abs() / fd.abs().max(1.0e-30);
                worst = worst.max(r);
            }
            run.num(&format!("KB_{}_{label}", scheme_name(scheme)), worst);
            worst_all = worst_all.max(worst);
        }
    }
    run.num("KB_worst_all", worst_all);
    run.gate(
        "K-B",
        format!(
            "FD-vs-analytic gradient THROUGH the new layer, N={N_FD}, both schemes, 8 objectives x 3 dirs: worst rel {worst_all:.2e} (tol 1e-5)"
        ),
        worst_all < TOL_GRAD,
    );
    vp(&format!("K-B FD gate done: worst rel {worst_all:.2e} ({:.1}s)", tb.elapsed().as_secs_f64()));

    // ---- K-C: NEW bit goldens for the tiled order --------------------------
    let mut golden = fs_blake3::Blake3::new();
    golden.update(fs_gum_kern::GUM_KERN_BIT_SEMANTICS.as_bytes());
    golden.update(b"gum-kern:goldens:v1");
    for scheme in [Scheme::Central4, Scheme::Corner] {
        let o = Opts::routhian(scheme, L_FD).with_guards(1.0, 10.0);
        let (out, g) = fs_gum_kern::eval(&f48, &o, true, 1);
        let g = g.expect("grad");
        for v in out_fields(&out) {
            golden.update(&v.to_bits().to_le_bytes());
        }
        for v in &g.data {
            golden.update(&v.to_bits().to_le_bytes());
        }
        let nrm = fs_gum_kern::norm_flat(&g.data, N_RUN, 1);
        golden.update(&nrm.to_bits().to_le_bytes());
        run.num(&format!("KC_obj_{}", scheme_name(scheme)), out.obj);
        run.num(&format!("KC_gnorm_{}", scheme_name(scheme)), nrm);
        if verbose {
            println!(
                "  K-C {}: obj {:.17e}  gnorm {:.17e}",
                scheme_name(scheme),
                out.obj,
                nrm
            );
        }
    }
    let golden = golden.finalize().to_hex();
    run.gate(
        "K-C",
        format!(
            "NEW bit goldens (tiled order, both schemes, Out+grad+norm): {golden} == frozen {GOLDEN_EVAL}"
        ),
        golden == GOLDEN_EVAL,
    );
    vp(&format!("  K-C golden hash: {golden}"));

    // ---- K-D: THE HEADLINE GATE — 1/2/3/4-thread hash identity -------------
    let td = Instant::now();
    // mini-descent input at N = 32 (over-spun rotor + guards, halo seed)
    let mut f32h = hedgehog(rp, N_HASH);
    let (out_h32, _) = fs_gum_kern::eval(&f32h, &Opts::estatic(Scheme::Corner), false, 1);
    let deg_ref32 = out_h32.deg;
    let fgap_ref32 = out_h32.e6 + out_h32.e0 - fs_gum_field::bps_floor() - FLOOR_BAND;
    let bump32 = bump_field(N_HASH, LBOX, SEED_PERT, 6, 0.02);
    let halo32 = halo_seed_field(N_HASH, LBOX, 0.05, 3.2, 0.6);
    add_scaled(&mut f32h, &bump32, 1.0);
    add_scaled(&mut f32h, &halo32, 1.0);
    let _ = f32h.renormalize();
    let o32 = Opts::routhian(Scheme::Corner, KAPPA_MAIN * out_h32.i)
        .with_guards(deg_ref32, fgap_ref32);
    let o48 = Opts::routhian(Scheme::Corner, L_FD).with_guards(1.0, 10.0);
    let mut hashes: Vec<String> = Vec::new();
    for threads in 1..=4usize {
        // (a) eval fwd+grad at N = 48, both schemes
        let mut parts: Vec<Vec<f64>> = Vec::new();
        for scheme in [Scheme::Central4, Scheme::Corner] {
            let o = Opts { scheme, ..o48 };
            let (out, g) = fs_gum_kern::eval(&f48, &o, true, threads);
            let g = g.expect("grad");
            parts.push(out_fields(&out).to_vec());
            parts.push(g.data);
        }
        // (b) renormalize a perturbed copy: drift + full field bytes
        let mut fr = clone_field(&f48);
        add_scaled(&mut fr, &bump, 0.5);
        let drift = fs_gum_kern::renormalize(&mut fr, threads);
        parts.push(vec![drift]);
        parts.push(fr.data().to_vec());
        // (c) 80-iteration arrested-Newton descent at N = 32
        let mut fd32 = clone_field(&f32h);
        let prm = AnfParams::new(MAXIT_HASH);
        let res = fs_gum_kern::anf(&mut fd32, &o32, &prm, "hash", threads);
        let (fs, us) = anf_stream(&res);
        parts.push(fs);
        parts.push(fd32.data().to_vec());
        let slices: Vec<&[f64]> = parts.iter().map(Vec::as_slice).collect();
        let h = hash_stream("gum-kern:thread-identity:v1", &slices, &us);
        vp(&format!("  K-D threads={threads}: {}", h.to_hex()));
        hashes.push(h.to_hex());
    }
    let all_equal = hashes.iter().all(|h| h == &hashes[0]);
    run.num("KD_all_equal", f64::from(u8::from(all_equal)));
    run.gate(
        "K-D",
        format!(
            "HEADLINE: serial vs threaded BIT-IDENTICAL at 1/2/3/4 threads (eval fwd+grad both schemes at N={N_RUN}, renormalize field+drift, {MAXIT_HASH}-iter ANF descent at N={N_HASH}): hash {}",
            &hashes[0][..16]
        ),
        all_equal,
    );
    vp(&format!("K-D thread-identity done ({:.1}s)", td.elapsed().as_secs_f64()));
    run
}

// ---------------------------------------------------------------------------
// K-E: the G-C main protocol through the new layer
// ---------------------------------------------------------------------------

struct Descent {
    first: Record,
    last: Record,
    max_rise_r: f64,
    max_rise_obj: f64,
    kappa_cross: Option<usize>,
}

fn descent_stats(res: &AnfResult) -> Descent {
    let first = res.series[0];
    let last = *res.series.last().expect("series");
    let mut max_rise_r = f64::NEG_INFINITY;
    let mut max_rise_obj = f64::NEG_INFINITY;
    for w in res.series.windows(2) {
        max_rise_r = max_rise_r.max(w[1].r - w[0].r);
        max_rise_obj = max_rise_obj.max(w[1].obj - w[0].obj);
    }
    let thr = fs_gum_statics::kappa_threshold();
    let mut kappa_cross = None;
    for w in res.series.windows(2) {
        if w[0].kappa >= thr && w[1].kappa < thr {
            kappa_cross = Some(w[1].it);
            break;
        }
    }
    Descent { first, last, max_rise_r, max_rise_obj, kappa_cross }
}

#[allow(clippy::too_many_lines)]
fn protocol(rp: &fs_gum_field::radial::RadialProfile, run: &mut Run) {
    println!("\n---- K-E: G-C main protocol at N={N_RUN} (new layer, {THREADS} threads; old engine for the endpoint bands) ----");
    let fh = hedgehog(rp, N_RUN);
    let bump2 = bump_field(N_RUN, LBOX, SEED_PERT, 6, 0.02);
    let halo_seed = halo_seed_field(N_RUN, LBOX, 0.05, 3.2, 0.6);

    // protocol closure over an engine: returns (static result, main result)
    let run_protocol = |name: &str, use_new: bool| -> (AnfResult, AnfResult, f64) {
        let t0 = Instant::now();
        let (out_h, _) = if use_new {
            fs_gum_kern::eval(&fh, &Opts::estatic(Scheme::Corner), false, THREADS)
        } else {
            fs_gum_statics::eval(&fh, &Opts::estatic(Scheme::Corner), false)
        };
        let deg_ref = out_h.deg;
        let fgap_ref = out_h.e6 + out_h.e0 - fs_gum_field::bps_floor() - FLOOR_BAND;
        let opts_static = Opts::estatic(Scheme::Corner).with_guards(deg_ref, fgap_ref);
        let mut f = clone_field(&fh);
        add_scaled(&mut f, &bump2, 1.0);
        let _ = f.renormalize();
        let mut prm = AnfParams::new(MAXIT_STATIC);
        prm.print_every = 300;
        let res_s = if use_new {
            fs_gum_kern::anf(&mut f, &opts_static, &prm, &format!("{name}-static"), THREADS)
        } else {
            fs_gum_statics::anf(&mut f, &opts_static, &prm, &format!("{name}-static"))
        };
        let mut q_static = vec![0.0_f64; 4 * N_RUN * N_RUN * N_RUN];
        save_cells(&f, &mut q_static);
        let l_main = KAPPA_MAIN * out_h.i;
        restore_cells(&mut f, &q_static);
        add_scaled(&mut f, &bump2, 1.0);
        add_scaled(&mut f, &halo_seed, 1.0);
        let _ = f.renormalize();
        let opts_main = Opts::routhian(Scheme::Corner, l_main).with_guards(deg_ref, fgap_ref);
        let mut prm = AnfParams::new(MAXIT_MAIN);
        prm.print_every = 200;
        let res_m = if use_new {
            fs_gum_kern::anf(&mut f, &opts_main, &prm, &format!("{name}-main"), THREADS)
        } else {
            fs_gum_statics::anf(&mut f, &opts_main, &prm, &format!("{name}-main"))
        };
        let dt = t0.elapsed().as_secs_f64();
        println!("  [{name}] protocol wall time {dt:.1}s");
        (res_s, res_m, dt)
    };

    let (new_s, new_m, t_new) = run_protocol("new", true);
    let (_old_s, old_m, t_old) = run_protocol("old", false);
    run.num("KE_t_new", t_new);
    run.num("KE_t_old", t_old);

    let dn = descent_stats(&new_m);
    let do_ = descent_stats(&old_m);
    let thr = fs_gum_statics::kappa_threshold();

    // qualitative referee on the NEW run (the G-C1/2/3 thresholds)
    let e1 = dn.last.r < dn.first.r - 0.05 && dn.max_rise_r <= 2.0e-3;
    run.gate(
        "K-E1",
        format!(
            "new-layer main R decreases: {:.6} -> {:.6} (dR {:+.5} <= -0.05), max recorded rise {:+.2e} <= 2e-3",
            dn.first.r,
            dn.last.r,
            dn.last.r - dn.first.r,
            dn.max_rise_r
        ),
        e1,
    );
    let dk = dn.first.kappa - dn.last.kappa;
    let e2 = dk >= 0.05;
    run.gate(
        "K-E2",
        format!(
            "new-layer kappa falls: {:.5} -> {:.5} (drop {dk:+.5} >= 0.05); crossing of threshold {thr:.5}: {} (old engine: {})",
            dn.first.kappa,
            dn.last.kappa,
            dn.kappa_cross.map_or("none".to_string(), |it| format!("by it {it}")),
            do_.kappa_cross.map_or("none".to_string(), |it| format!("by it {it}")),
        ),
        e2,
    );
    let dh = dn.last.halo - dn.first.halo;
    let e3 = dh >= 0.015;
    run.gate(
        "K-E3",
        format!(
            "new-layer halo rises (the 4A signature): {:.4} -> {:.4} (rise {dh:+.4} >= 0.015)",
            dn.first.halo, dn.last.halo
        ),
        e3,
    );
    let e4 = dn.max_rise_obj <= 0.0 && new_s.series.windows(2).all(|w| w[1].obj <= w[0].obj);
    run.gate(
        "K-E4",
        format!(
            "new-layer objectives monotone by construction (static + main): max recorded rise {:+.2e} <= 0",
            dn.max_rise_obj
        ),
        e4,
    );

    // endpoint tolerance bands vs the old engine (documented: need NOT
    // bit-match — the arrest cascade amplifies the order-change ulps)
    let d_r = dn.last.r - do_.last.r;
    let d_kap = dn.last.kappa - do_.last.kappa;
    let d_halo = dn.last.halo - do_.last.halo;
    let d_deg = dn.last.deg - do_.last.deg;
    for (k, v) in [
        ("KE_new_Rf", dn.last.r),
        ("KE_old_Rf", do_.last.r),
        ("KE_new_kapf", dn.last.kappa),
        ("KE_old_kapf", do_.last.kappa),
        ("KE_new_halof", dn.last.halo),
        ("KE_old_halof", do_.last.halo),
        ("KE_new_degf", dn.last.deg),
        ("KE_old_degf", do_.last.deg),
        ("KE_new_arrests", new_m.arrests as f64),
        ("KE_old_arrests", old_m.arrests as f64),
    ] {
        run.num(k, v);
    }
    let e5 = d_r.abs() <= BAND_R
        && d_kap.abs() <= BAND_KAPPA
        && d_halo.abs() <= BAND_HALO
        && d_deg.abs() <= BAND_DEG;
    run.gate(
        "K-E5",
        format!(
            "endpoint vs old engine within tolerance bands: dR {d_r:+.2e} (band {BAND_R}), dkappa {d_kap:+.2e} ({BAND_KAPPA}), dhalo {d_halo:+.2e} ({BAND_HALO}), ddeg {d_deg:+.2e} ({BAND_DEG}); arrests new {} vs old {} — trajectories legitimately diverge at ulp scale, bands are the sanctioned golden class",
            new_m.arrests, old_m.arrests
        ),
        e5,
    );
    println!(
        "  endpoint: R {:.17e} vs {:.17e} | kappa {:.9} vs {:.9} | halo {:.6} vs {:.6} | deg {:.7} vs {:.7}",
        dn.last.r, do_.last.r, dn.last.kappa, do_.last.kappa, dn.last.halo, do_.last.halo,
        dn.last.deg, do_.last.deg
    );
}

// ---------------------------------------------------------------------------
// K-F: measured performance
// ---------------------------------------------------------------------------

fn bench<F: FnMut() -> f64>(reps: usize, mut f: F) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..reps {
        let t0 = Instant::now();
        black_box(f());
        best = best.min(t0.elapsed().as_secs_f64());
    }
    best
}

struct PerfRow {
    label: String,
    pts: f64,
    t_old: f64,
    t_new1: f64,
    t_new4: f64,
}

fn perf(rp: &fs_gum_field::radial::RadialProfile, run: &mut Run) {
    println!("\n---- K-F: measured ns/point, serial-old vs serial-new vs threaded-4 ----");
    let mut rows: Vec<PerfRow> = Vec::new();
    let mut corner_grad_96 = (0.0_f64, 0.0_f64); // (serial-new, threaded-4)
    let mut extra_96 = (0.0_f64, 0.0_f64); // threaded-2 / threaded-3 corner grad
    for &n in &[48usize, 96, 128] {
        let f = hedgehog(rp, n);
        let reps = if n >= 128 { 2 } else { 3 };
        for scheme in [Scheme::Central4, Scheme::Corner] {
            let pts = match scheme {
                Scheme::Central4 => (n * n * n) as f64,
                Scheme::Corner => ((n + 1) * (n + 1) * (n + 1)) as f64,
                Scheme::Central2 => unreachable!(),
            };
            // the ANF-class objective (rotor + guards active)
            let (oh, _) = fs_gum_statics::eval(&f, &Opts::estatic(scheme), false);
            let o = Opts::routhian(scheme, KAPPA_MAIN * oh.i).with_guards(
                oh.deg,
                oh.e6 + oh.e0 - fs_gum_field::bps_floor() - FLOOR_BAND,
            );
            let mut ws = fs_gum_kern::Ws::new();
            for need_grad in [false, true] {
                let t_old = bench(reps, || fs_gum_statics::eval(&f, &o, need_grad).0.obj);
                let t_new1 =
                    bench(reps, || fs_gum_kern::eval_ws(&f, &o, need_grad, 1, &mut ws).0.obj);
                let t_new4 = bench(reps, || {
                    fs_gum_kern::eval_ws(&f, &o, need_grad, THREADS, &mut ws).0.obj
                });
                let label = format!(
                    "{} {:<8} N={n}",
                    if need_grad { "grad" } else { "fwd " },
                    scheme_name(scheme)
                );
                if need_grad && scheme == Scheme::Corner && n == 96 {
                    corner_grad_96 = (t_new1, t_new4);
                    let t2 =
                        bench(2, || fs_gum_kern::eval_ws(&f, &o, need_grad, 2, &mut ws).0.obj);
                    let t3 =
                        bench(2, || fs_gum_kern::eval_ws(&f, &o, need_grad, 3, &mut ws).0.obj);
                    extra_96 = (t2, t3);
                }
                rows.push(PerfRow { label, pts, t_old, t_new1, t_new4 });
            }
        }
    }
    println!(
        "  {:<22} | {:>9} | {:>9} | {:>9} | {:>7} | {:>7}",
        "sweep", "old ns/pt", "new ns/pt", "4T ns/pt", "4T/new", "new/old"
    );
    for r in &rows {
        println!(
            "  {:<22} | {:>9.1} | {:>9.1} | {:>9.1} | {:>6.2}x | {:>6.2}x",
            r.label,
            r.t_old * 1e9 / r.pts,
            r.t_new1 * 1e9 / r.pts,
            r.t_new4 * 1e9 / r.pts,
            r.t_new1 / r.t_new4,
            r.t_old / r.t_new1
        );
        run.num(&format!("KF_{}_old", r.label), r.t_old * 1e9 / r.pts);
        run.num(&format!("KF_{}_new1", r.label), r.t_new1 * 1e9 / r.pts);
        run.num(&format!("KF_{}_new4", r.label), r.t_new4 * 1e9 / r.pts);
    }
    let (t1, t4) = corner_grad_96;
    let speedup = t1 / t4;
    let (t2, t3) = extra_96;
    println!(
        "  corner grad N=96 thread ladder: 1T {:.3}s, 2T {:.3}s ({:.2}x), 3T {:.3}s ({:.2}x), 4T {:.3}s ({:.2}x)",
        t1,
        t2,
        t1 / t2,
        t3,
        t1 / t3,
        t4,
        speedup
    );
    run.num("KF_speedup_corner96", speedup);
    run.gate(
        "K-F",
        format!(
            "threaded-4 speedup on the ANF hot path (corner gradient, N=96): {speedup:.2}x >= 3.0x vs serial-new ({t1:.3}s -> {t4:.3}s)"
        ),
        speedup >= 3.0,
    );
}

// ---------------------------------------------------------------------------

fn main() {
    let t_start = Instant::now();
    println!(
        "fs-gum-kern gates (Phase G1) — TILE_I={}, threads<={THREADS}, N_fd={N_FD}, N_hash={N_HASH}, N_run={N_RUN}, LBOX={LBOX}, t={T_FROZEN:.15}",
        fs_gum_kern::TILE_I
    );
    println!("bit contract: {}", fs_gum_kern::GUM_KERN_BIT_SEMANTICS);
    let rp = radial_solve(T_FROZEN, 4000, 6.0);

    println!("\n================ core run 1 (K-A..K-D) ================");
    let mut run1 = core(&rp, true);
    println!("\n================ core run 2 (replay, silent) ================");
    let run2 = core(&rp, false);

    // K-G: bit-identical replay of the core
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
    let fp = fs_blake3::hash_domain("gum-kern:fingerprint:v1", &bytes);

    protocol(&rp, &mut run1);
    perf(&rp, &mut run1);

    println!("\n================ gate table ================");
    let mut all_pass = true;
    for g in run1.gates.iter() {
        println!("  {:<5} {}  -> {}", g.id, g.desc, if g.pass { "PASS" } else { "FAIL" });
        all_pass &= g.pass;
    }
    println!(
        "  K-G   bit-identical two-run replay of {} core gate numbers  -> {}",
        run2.nums.len(),
        if identical { "PASS" } else { "FAIL" }
    );
    all_pass &= identical;
    println!("\nrun fingerprint (BLAKE3 over every core f64, fixed order): {}", fp.to_hex());
    println!("total runtime {:.1}s", t_start.elapsed().as_secs_f64());
    println!("\nOVERALL: {}", if all_pass { "PASS" } else { "FAIL" });
    if !all_pass {
        std::process::exit(1);
    }
}
