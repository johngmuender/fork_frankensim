//! Arrested Newton flow — the exact field3d `anf()` protocol.
//!
//! Velocity field v on the quaternion components; acceleration is the
//! NEGATIVE h^3-normalised gradient (v -= (dt/h^3) g); the field steps
//! q += dt v and is renormalised; on ANY objective increase the step is
//! reverted, the velocity zeroed and dt shrunk 0.6x (arrest); on acceptance
//! the velocity is re-projected onto the tangent space of the updated field
//! and dt grows 1.01x up to dt_max = 0.05; dt < 1e-7 is a stall.  The
//! objective is monotone non-increasing BY CONSTRUCTION.

use fs_gum_field::Field3;

use crate::diag::{aidx, halo_fraction, restore_cells, save_cells};
use crate::engine::{eval, Opts, Out};

/// ANF hyper-parameters (defaults are the frozen field3d values).
#[derive(Clone, Copy, Debug)]
pub struct AnfParams {
    pub maxit: usize,
    pub dt0: f64,
    pub dt_max: f64,
    /// Stop when |g| < gtol_ratio * |g_initial|.
    pub gtol_ratio: f64,
    /// Record the instrumented series every this many iterations.
    pub instr_every: usize,
    /// Print a progress line every this many iterations (0 = silent).
    pub print_every: usize,
}

impl AnfParams {
    #[must_use]
    pub fn new(maxit: usize) -> Self {
        AnfParams {
            maxit,
            dt0: 0.01,
            dt_max: 0.05,
            gtol_ratio: 1.0e-6,
            instr_every: 10,
            print_every: 0,
        }
    }
}

/// Termination status.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// Iteration budget exhausted.
    IterCap,
    /// dt shrank below 1e-7 (arrest stall).
    StallDt,
    /// Gradient ratio below gtol_ratio.
    Gtol,
}

impl Status {
    #[must_use]
    pub fn as_str(&self) -> &'static str {
        match self {
            Status::IterCap => "iter_cap",
            Status::StallDt => "stall_dt",
            Status::Gtol => "gtol",
        }
    }
}

/// One instrumented sample of the descent (field3d `record`).
#[derive(Clone, Copy, Debug)]
pub struct Record {
    pub it: usize,
    pub r: f64,
    pub obj: f64,
    pub epen: f64,
    pub efpen: f64,
    pub floor_gap: f64,
    pub estat: f64,
    pub e2: f64,
    pub e4: f64,
    pub e6: f64,
    pub e0: f64,
    pub i: f64,
    pub kappa: f64,
    pub deg: f64,
    pub halo: f64,
    pub gnorm: f64,
    pub dt: f64,
    pub arrests: usize,
}

/// Descent result: final Out (last ACCEPTED point; the field is left at that
/// point), the instrumented series, and the run metadata.
pub struct AnfResult {
    pub out: Out,
    pub series: Vec<Record>,
    pub status: Status,
    pub iters: usize,
    pub arrests: usize,
    pub gn0: f64,
    pub gn: f64,
}

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

/// Run arrested Newton flow on the field in place.  On return the field
/// holds the last accepted configuration.
pub fn anf(f: &mut Field3, o: &Opts, p: &AnfParams, label: &str) -> AnfResult {
    let n = f.n();
    let h = f.h();
    let h3 = h * h * h;
    let ncell4 = 4 * n * n * n;

    let (mut out, gopt) = eval(f, o, true);
    let mut g = gopt.expect("gradient requested");
    let mut obj = out.obj;
    let gn0 = g.norm() / h3;
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
        // v -= (dt/h^3) g;  q += dt v;  renormalize
        let a = dt / h3;
        for (vm, &gm) in v.iter_mut().zip(g.data.iter()) {
            *vm -= a * gm;
        }
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let mut q = f.get(i, j, k);
                    for (c, qv) in q.iter_mut().enumerate() {
                        *qv += dt * v[aidx(n, c, i, j, k)];
                    }
                    f.set(i, j, k, q);
                }
            }
        }
        let _drift = f.renormalize();
        let (out_new, gnew) = eval(f, o, true);
        if out_new.obj > obj {
            // arrest: revert, kill the velocity, shrink dt
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
            // re-project velocity onto the new tangent space
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let q = f.get(i, j, k);
                        let mut dot = 0.0_f64;
                        for (c, &qv) in q.iter().enumerate() {
                            dot += v[aidx(n, c, i, j, k)] * qv;
                        }
                        for (c, &qv) in q.iter().enumerate() {
                            v[aidx(n, c, i, j, k)] -= dot * qv;
                        }
                    }
                }
            }
            dt = (dt * 1.01).min(p.dt_max);
        }
        gn = g.norm() / h3;
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
