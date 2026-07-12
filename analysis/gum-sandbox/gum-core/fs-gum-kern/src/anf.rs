//! Arrested Newton flow re-hosted on the tiled deterministic layer — the
//! exact `fs_gum_statics::anf` protocol (same types, same arrest/accept
//! controller, same instrumented series) with every field sweep routed
//! through `fs-gum-kern`: eval (tiled reductions + gather adjoint),
//! renormalize / velocity update / position step / tangent re-projection
//! as parallel maps, g.norm as a tiled reduction.
//!
//! The controller itself (revert + zero-velocity arrest, dt cascade,
//! accept bookkeeping) stays SERIAL and verbatim: it consumes ~12 scalars
//! per iteration and is inherently sequential.  Because every number it
//! branches on is bit-identical across thread counts by construction,
//! the whole descent — trajectory, arrests, series, endpoint — is
//! bit-identical at 1..=N threads (asserted by gate K-D).
//!
//! Versus the OLD engine the descent is NOT expected to bit-match: the
//! tiled/gather accumulation order differs at machine-eps class and the
//! arrest cascade amplifies ulps into diverging trajectories (survey
//! DIMENSION 1, hazard c) — the endpoint comparison is a tolerance-band
//! gate (K-E), per the sanctioned cross-backend golden policy.

use fs_gum_field::Field3;
use fs_gum_statics::diag::{halo_fraction, restore_cells, save_cells};
use fs_gum_statics::{AnfParams, AnfResult, Opts, Out, Record, Status};

use crate::engine::{axpy_sub, eval_ws, norm_flat, renormalize, step_q, tangent_project, Ws};

fn record(f: &Field3, o: &Opts, out: &Out, it: usize, gn: f64, dt: f64, arrests: usize) -> Record {
    let (halo, _cen) = halo_fraction(f);
    Record {
        it,
        r: out.r,
        obj: out.obj,
        epen: out.epen,
        efpen: out.efpen,
        floor_gap: out.floor_gap,
        estat: out.estat,
        e2: out.e2,
        e4: out.e4,
        e6: out.e6,
        e0: out.e0,
        i: out.i,
        kappa: o.l.map_or(0.0, |l| l / out.i),
        deg: out.deg,
        halo,
        gnorm: gn,
        dt,
        arrests,
    }
}

/// Run arrested Newton flow on the field in place, on `threads` threads.
/// On return the field holds the last accepted configuration.  Protocol
/// and hyper-parameters are `fs_gum_statics::anf` verbatim.
pub fn anf(f: &mut Field3, o: &Opts, p: &AnfParams, label: &str, threads: usize) -> AnfResult {
    let n = f.n();
    let h = f.h();
    let h3 = h * h * h;
    let ncell4 = 4 * n * n * n;

    let mut ws = Ws::new();
    let (mut out, gopt) = eval_ws(f, o, true, threads, &mut ws);
    let mut g = gopt.expect("gradient requested");
    let mut obj = out.obj;
    let gn0 = norm_flat(&g.data, n, threads) / h3;
    let mut gn = gn0;
    let mut v = vec![0.0_f64; ncell4];
    let mut qb = vec![0.0_f64; ncell4];
    save_cells(f, &mut qb);
    let mut dt = p.dt0;
    let mut arrests = 0usize;
    let mut series = Vec::new();
    let mut status = Status::IterCap;

    series.push(record(f, o, &out, 0, gn0, dt, arrests));
    if p.print_every > 0 {
        let r0 = series[0];
        println!(
            "  [{label}] it=0  R={:.7} I={:.4} kappa={:.5} deg={:.5} halo={:.4} gnorm={:.3e}",
            out.r, out.i, r0.kappa, out.deg, r0.halo, gn0
        );
    }
    let mut it = 0usize;
    while it < p.maxit {
        it += 1;
        // v -= (dt/h^3) g;  q += dt v;  renormalize — all tiled maps
        axpy_sub(&mut v, dt / h3, &g.data, threads);
        step_q(f, &v, dt, threads);
        let _drift = renormalize(f, threads);
        let (out_new, gnew) = eval_ws(f, o, true, threads, &mut ws);
        if out_new.obj > obj {
            // arrest: revert, kill the velocity, shrink dt (serial controller)
            restore_cells(f, &qb);
            for vm in v.iter_mut() {
                *vm = 0.0;
            }
            dt *= 0.6;
            arrests += 1;
            if dt < 1.0e-7 {
                status = Status::StallDt;
                break;
            }
            // g of the accepted point is still valid (kept)
        } else {
            out = out_new;
            obj = out.obj;
            g = gnew.expect("gradient requested");
            save_cells(f, &mut qb);
            // re-project velocity onto the new tangent space (tiled map)
            tangent_project(f, &mut v, threads);
            dt = (dt * 1.01).min(p.dt_max);
        }
        gn = norm_flat(&g.data, n, threads) / h3;
        if it % p.instr_every == 0 || it == p.maxit {
            let rec = record(f, o, &out, it, gn, dt, arrests);
            series.push(rec);
            if p.print_every > 0 && (it % p.print_every == 0 || it == p.maxit) {
                println!(
                    "  [{label}] it={it}  R={:.7} pen={:.2e} fpen={:.2e} fgap={:+.4} I={:.4} kappa={:.5} deg={:.5} halo={:.4} gnorm={:.3e} dt={:.2e} arr={arrests}",
                    out.r, out.epen, out.efpen, out.floor_gap, out.i, rec.kappa, out.deg, rec.halo, gn, dt
                );
            }
        }
        if gn < p.gtol_ratio * gn0 {
            status = Status::Gtol;
            series.push(record(f, o, &out, it, gn, dt, arrests));
            break;
        }
    }
    restore_cells(f, &qb);
    if p.print_every > 0 {
        println!(
            "  [{label}] DONE ({}) it={it} R={:.7} pen={:.2e} gnorm/gnorm0={:.3e} arrests={arrests}",
            status.as_str(),
            out.r,
            out.epen,
            gn / gn0
        );
    }
    AnfResult { out, series, status, iters: it, arrests, gn0, gn }
}
