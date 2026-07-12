//! Arrested Newton flow on the anisotropic field — the exact field3d /
//! fs-gum-statics `anf()` protocol with anisotropic bounds.  The halo
//! diagnostic of the cubic Record is not carried (it is a single-knot
//! instrument); everything that feeds the descent is identical.

use crate::engine_a::{eval_a, GradA, OptsA, OutA};
use crate::field_a::{aidx_a, FieldA};

/// ANF hyper-parameters (frozen field3d defaults).
#[derive(Clone, Copy, Debug)]
pub struct AnfParamsA {
    pub maxit: usize,
    pub dt0: f64,
    pub dt_max: f64,
    pub gtol_ratio: f64,
    pub instr_every: usize,
    pub print_every: usize,
}

impl AnfParamsA {
    #[must_use]
    pub fn new(maxit: usize) -> Self {
        AnfParamsA {
            maxit,
            dt0: 0.01,
            dt_max: 0.05,
            gtol_ratio: 1.0e-6,
            instr_every: 10,
            print_every: 0,
        }
    }
}

/// Termination status (fs-gum-statics `Status`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StatusA {
    IterCap,
    StallDt,
    Gtol,
}

impl StatusA {
    #[must_use]
    pub fn as_str(&self) -> &'static str {
        match self {
            StatusA::IterCap => "iter_cap",
            StatusA::StallDt => "stall_dt",
            StatusA::Gtol => "gtol",
        }
    }
}

/// One instrumented sample.
#[derive(Clone, Copy, Debug)]
pub struct RecordA {
    pub it: usize,
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
    pub deg: f64,
    pub gnorm: f64,
    pub dt: f64,
    pub arrests: usize,
}

/// Descent result.
pub struct AnfResultA {
    pub out: OutA,
    pub series: Vec<RecordA>,
    pub status: StatusA,
    pub iters: usize,
    pub arrests: usize,
    pub gn0: f64,
    pub gn: f64,
}

fn record(out: &OutA, it: usize, gn: f64, dt: f64, arrests: usize) -> RecordA {
    RecordA {
        it,
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
        deg: out.deg,
        gnorm: gn,
        dt,
        arrests,
    }
}

fn save_cells_a(f: &FieldA, buf: &mut [f64]) {
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let q = f.get(i, j, k);
                for (a, &v) in q.iter().enumerate() {
                    buf[aidx_a(nx, ny, nz, a, i, j, k)] = v;
                }
            }
        }
    }
}

fn restore_cells_a(f: &mut FieldA, buf: &[f64]) {
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let mut q = [0.0_f64; 4];
                for (a, qv) in q.iter_mut().enumerate() {
                    *qv = buf[aidx_a(nx, ny, nz, a, i, j, k)];
                }
                f.set(i, j, k, q);
            }
        }
    }
}

/// Run arrested Newton flow in place (fs-gum-statics `anf`, anisotropic).
#[allow(clippy::too_many_lines)]
pub fn anf_a(f: &mut FieldA, o: &OptsA, p: &AnfParamsA, label: &str) -> AnfResultA {
    let (nx, ny, nz) = (f.nx(), f.ny(), f.nz());
    let h = f.h();
    let h3 = h * h * h;
    let ncell4 = 4 * nx * ny * nz;

    let (mut out, gopt) = eval_a(f, o, true);
    let mut g: GradA = gopt.expect("gradient requested");
    let mut obj = out.obj;
    let gn0 = g.norm() / h3;
    let mut gn = gn0;
    let mut v = vec![0.0_f64; ncell4];
    let mut qb = vec![0.0_f64; ncell4];
    save_cells_a(f, &mut qb);
    let mut dt = p.dt0;
    let mut arrests = 0usize;
    let mut series = Vec::new();
    let mut status = StatusA::IterCap;

    series.push(record(&out, 0, gn0, dt, arrests));
    if p.print_every > 0 {
        println!(
            "  [{label}] it=0  obj={:.9} Estat={:.7} deg={:.5} gnorm={:.3e}",
            out.obj, out.estat, out.deg, gn0
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
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let mut q = f.get(i, j, k);
                    for (c, qv) in q.iter_mut().enumerate() {
                        *qv += dt * v[aidx_a(nx, ny, nz, c, i, j, k)];
                    }
                    f.set(i, j, k, q);
                }
            }
        }
        let _drift = f.renormalize();
        let (out_new, gnew) = eval_a(f, o, true);
        if out_new.obj > obj {
            // arrest: revert, kill the velocity, shrink dt
            restore_cells_a(f, &qb);
            for vm in v.iter_mut() {
                *vm = 0.0;
            }
            dt *= 0.6;
            arrests += 1;
            if dt < 1.0e-7 {
                status = StatusA::StallDt;
                break;
            }
        } else {
            out = out_new;
            obj = out.obj;
            g = gnew.expect("gradient requested");
            save_cells_a(f, &mut qb);
            // re-project velocity onto the new tangent space
            for i in 0..nx {
                for j in 0..ny {
                    for k in 0..nz {
                        let q = f.get(i, j, k);
                        let mut dot = 0.0_f64;
                        for (c, &qv) in q.iter().enumerate() {
                            dot += v[aidx_a(nx, ny, nz, c, i, j, k)] * qv;
                        }
                        for (c, &qv) in q.iter().enumerate() {
                            v[aidx_a(nx, ny, nz, c, i, j, k)] -= dot * qv;
                        }
                    }
                }
            }
            dt = (dt * 1.01).min(p.dt_max);
        }
        gn = g.norm() / h3;
        if it % p.instr_every == 0 || it == p.maxit {
            let rec = record(&out, it, gn, dt, arrests);
            series.push(rec);
            if p.print_every > 0 && (it % p.print_every == 0 || it == p.maxit) {
                println!(
                    "  [{label}] it={it}  obj={:.9} Estat={:.7} pen={:.2e} fpen={:.2e} fgap={:+.4} deg={:.5} gnorm={:.3e} dt={:.2e} arr={arrests}",
                    out.obj, out.estat, out.epen, out.efpen, out.floor_gap, out.deg, gn, dt
                );
            }
        }
        if gn < p.gtol_ratio * gn0 {
            status = StatusA::Gtol;
            series.push(record(&out, it, gn, dt, arrests));
            break;
        }
    }
    restore_cells_a(f, &qb);
    if p.print_every > 0 {
        println!(
            "  [{label}] DONE ({}) it={it} obj={:.9} gnorm/gnorm0={:.3e} arrests={arrests}",
            status.as_str(),
            out.obj,
            gn / gn0
        );
    }
    AnfResultA { out, series, status, iters: it, arrests, gn0, gn }
}
