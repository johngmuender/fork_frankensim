//! One campaign run: seed (two-knot product ansatz at separation d = 2 m h
//! in a channel, or the single-knot reference), guarded ANF relax at the
//! frozen eps = 0.05 dial, record energies + diagnostics, write one JSON.
//!
//! Usage: twoknot_run <channel> <m> <maxit> <outdir>
//!   channel  attract | repulse | align | single
//!   m        half-separation in grid cells (d = 2 m h); ignored for single
//!   maxit    ANF iteration cap (uniform across the campaign)
//!   outdir   directory for the per-run JSON
//!
//! Deterministic: no RNG anywhere (analytic seed), plain sequential loops;
//! same args => bit-identical JSON.  Runs are independent processes — the
//! task farm is `farm.sh`.

use fs_gum_field::radial::radial_solve;
use fs_gum_field::{bps_floor, rstar, Scheme, T_FROZEN};
use fs_gum_twoknot::{
    anf_a, effective_separation, eval_a, mu_frozen, seed_single_knot, seed_two_knot, AnfParamsA,
    Channel, FieldA, OptsA, OutA, H_FROZEN, NT, NX, RADIAL_N, RADIAL_RMAX,
};

use std::fmt::Write as _;
use std::time::Instant;

fn je(v: f64) -> String {
    if v.is_finite() {
        format!("{v:e}")
    } else {
        "null".to_string()
    }
}

fn out_json(o: &OutA) -> String {
    format!(
        "{{\"e2\": {}, \"e4\": {}, \"e6\": {}, \"e0\": {}, \"i\": {}, \"deg\": {}, \"estat\": {}, \"floor_gap\": {}, \"epen\": {}, \"efpen\": {}, \"obj\": {}}}",
        je(o.e2),
        je(o.e4),
        je(o.e6),
        je(o.e0),
        je(o.i),
        je(o.deg),
        je(o.estat),
        je(o.floor_gap),
        je(o.epen),
        je(o.efpen),
        je(o.obj)
    )
}

#[allow(clippy::too_many_lines)]
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("usage: twoknot_run <attract|repulse|align|single> <m> <maxit> <outdir>");
        std::process::exit(2);
    }
    let chan_s = args[1].as_str();
    let m: usize = args[2].parse().expect("m");
    let maxit: usize = args[3].parse().expect("maxit");
    let outdir = args[4].clone();
    let single = chan_s == "single";
    let channel = if single { None } else { Some(Channel::parse(chan_s).expect("channel")) };

    let t_all = Instant::now();
    let h = H_FROZEN;
    let s = m as f64 * h; // half-separation; d = 2 s
    let d = 2.0 * s;
    let x_sep = d / rstar();
    let b_charge = if single { 1.0 } else { 2.0 };

    println!(
        "twoknot_run: channel={chan_s} m={m} d={d:.6} x=d/R*={x_sep:.4} grid {NX}x{NT}x{NT} h={h} maxit={maxit}"
    );

    // radial profile (frozen protocol)
    let rp = radial_solve(T_FROZEN, RADIAL_N, RADIAL_RMAX);

    // seed
    let mut f = FieldA::new(NX, NT, NT, h);
    if single {
        seed_single_knot(&mut f, &rp);
    } else {
        seed_two_knot(&mut f, &rp, s, channel.unwrap());
    }
    let norm_drift = f.renormalize();

    // seed measurements (bare)
    let (seed_c, _) = eval_a(&f, &OptsA::estatic(Scheme::Corner), false);
    let (seed_c4, _) = eval_a(&f, &OptsA::estatic(Scheme::Central4), false);
    let seed_deff = if single { 0.0 } else { effective_separation(&f) };
    let seed_tail = f.boundary_tail();

    // guards referenced to the seed (statics protocol, B-charge floor unit)
    let deg_ref = seed_c.deg;
    let deg_unit = deg_ref / b_charge;
    let fgap_seed = seed_c.e6 + seed_c.e0 - bps_floor() * seed_c.deg / deg_unit;
    let fgap_ref = fgap_seed - fs_gum_statics::FLOOR_BAND;
    let opts = OptsA::estatic(Scheme::Corner).with_guards(deg_ref, deg_unit, fgap_ref);

    println!(
        "  seed: Estat(corner)={:.9} deg={:.6} fgap={:+.6} -> wall {:+.6} deff={:.4} tail={:.2e} drift={:.2e}",
        seed_c.estat, seed_c.deg, fgap_seed, fgap_ref, seed_deff, seed_tail, norm_drift
    );

    // guarded relaxation
    let mut prm = AnfParamsA::new(maxit);
    prm.instr_every = 10;
    prm.print_every = 50;
    let t_anf = Instant::now();
    let res = anf_a(&mut f, &opts, &prm, chan_s);
    let anf_secs = t_anf.elapsed().as_secs_f64();

    // endpoint measurements
    let fin_c = res.out;
    let (fin_c4, _) = eval_a(&f, &OptsA::estatic(Scheme::Central4), false);
    let fin_deff = if single { 0.0 } else { effective_separation(&f) };
    let fin_tail = f.boundary_tail();

    println!(
        "  final: Estat(corner)={:.9} Estat(central4)={:.9} deg={:.6} deff={:.4} status={} iters={} arrests={} ({:.1}s, {:.3} s/iter)",
        fin_c.estat,
        fin_c4.estat,
        fin_c.deg,
        fin_deff,
        res.status.as_str(),
        res.iters,
        res.arrests,
        anf_secs,
        anf_secs / res.iters.max(1) as f64
    );

    // ---- JSON ---------------------------------------------------------------
    let mut sj = String::from("[");
    for (i, r) in res.series.iter().enumerate() {
        if i > 0 {
            sj.push_str(", ");
        }
        let _ = write!(
            sj,
            "{{\"it\": {}, \"obj\": {}, \"estat\": {}, \"e2\": {}, \"e4\": {}, \"e6\": {}, \"e0\": {}, \"deg\": {}, \"floor_gap\": {}, \"epen\": {}, \"efpen\": {}, \"gnorm\": {}, \"dt\": {}, \"arrests\": {}}}",
            r.it,
            je(r.obj),
            je(r.estat),
            je(r.e2),
            je(r.e4),
            je(r.e6),
            je(r.e0),
            je(r.deg),
            je(r.floor_gap),
            je(r.epen),
            je(r.efpen),
            je(r.gnorm),
            je(r.dt),
            r.arrests
        );
    }
    sj.push(']');

    let json = format!(
        "{{\n\
         \"protocol\": {{\"channel\": \"{chan_s}\", \"m\": {m}, \"d\": {}, \"x_sep\": {}, \"nx\": {NX}, \"nt\": {NT}, \"h\": {}, \"t\": {}, \"mu\": {}, \"rstar\": {}, \"maxit\": {maxit}, \"scheme\": \"corner\", \"deg_ref\": {}, \"deg_unit\": {}, \"fgap_ref\": {}, \"mu_deg\": {}, \"mu_floor\": {}, \"radial_n\": {RADIAL_N}, \"radial_rmax\": {}}},\n\
         \"seed\": {{\"corner\": {}, \"central4\": {}, \"d_eff\": {}, \"boundary_tail\": {}, \"norm_drift\": {}}},\n\
         \"final\": {{\"corner\": {}, \"central4\": {}, \"d_eff\": {}, \"boundary_tail\": {}}},\n\
         \"anf\": {{\"status\": \"{}\", \"iters\": {}, \"arrests\": {}, \"gn0\": {}, \"gn\": {}, \"seconds\": {}, \"sec_per_iter\": {}}},\n\
         \"series\": {sj}\n\
         }}\n",
        je(d),
        je(x_sep),
        je(h),
        je(T_FROZEN),
        je(mu_frozen()),
        je(rstar()),
        je(deg_ref),
        je(deg_unit),
        je(fgap_ref),
        je(fs_gum_statics::MU_DEG),
        je(fs_gum_statics::MU_FLOOR),
        je(RADIAL_RMAX),
        out_json(&seed_c),
        out_json(&seed_c4),
        je(seed_deff),
        je(seed_tail),
        je(norm_drift),
        out_json(&fin_c),
        out_json(&fin_c4),
        je(fin_deff),
        je(fin_tail),
        res.status.as_str(),
        res.iters,
        res.arrests,
        je(res.gn0),
        je(res.gn),
        je(anf_secs),
        je(anf_secs / res.iters.max(1) as f64),
    );

    let fname = if single {
        format!("{outdir}/run_single.json")
    } else {
        format!("{outdir}/run_{chan_s}_m{m:02}.json")
    };
    std::fs::create_dir_all(&outdir).expect("outdir");
    std::fs::write(&fname, json).expect("write json");
    println!("  wrote {fname}  (total {:.1}s)", t_all.elapsed().as_secs_f64());
}
