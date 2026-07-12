//! GUM core Phase G2 — `fs-gum-gpu`: the wgpu f64 GPU backend for the gum
//! sweep kernels, exercised HERE on lavapipe (Mesa llvmpipe, software
//! Vulkan) as a correctness rig.  A real GPU inherits this code unchanged —
//! only the adapter that `GpuEngine::new` selects differs.
//!
//! What runs on the device (all f64 WGSL, `Features::SHADER_F64` / naga
//! `Capabilities::FLOAT64`):
//!
//!  * the field, resident as 4 per-component f64 storage buffers in
//!    fs-gum-field's padded SoA layout (P = N + 2*GHOST, offset
//!    (ip*P + jp)*P + kp), ghost rind included;
//!  * `fill_ghosts` (fixed-vacuum rind), `renormalize` (pointwise unit fix);
//!  * the central4 forward sweep producing per-workgroup partial sums for
//!    (s_e2, s_e4, s_det, s_det2, s_e0, s_i) — fixed workgroup size 64, one
//!    thread per interior cell, workgroup-shared reduction in fixed
//!    ascending lane order; the ~KB partials buffer is read back and
//!    combined on the CPU in fixed ascending workgroup order;
//!  * the gradient as pointwise flux pass (P[3][4], gq[4] per point, the
//!    exact `fs_gum_statics::engine::point_flux` math) plus a GATHER
//!    adjoint (transposed central4 stencil, ghost flux dropped), cell
//!    density terms, h^3 weight and tangent projection fused in;
//!  * the ANF inner loop (velocity kick, field step, velocity
//!    re-projection, snapshot save/restore as device-side buffer copies).
//!
//! What stays on the CPU: the fixed-order combine of workgroup partials
//! (~KBs per eval), the guard bookkeeping (dev/devf and the guard-dressed
//! kernel weights w6s/cstw/e0c — ~10 scalars written to a params buffer),
//! and the ANF arrest logic (one objective readback per iteration).
//!
//! Determinism / golden policy (Phase G2 contract, stated per gate in
//! RESULTS.md): same-device replay of any kernel sequence is expected
//! bit-identical and is GATED (dispatch twice, compare buffer bits);
//! CPU-vs-GPU agreement is the tolerance-band/metric class sanctioned by
//! ROADMAP_v4 for backend boundaries — the shader compiler may FMA-contract
//! (observed on llvmpipe), and the gather adjoint plus workgroup reduction
//! deliberately reorder the CPU engine's accumulation — so sector sums and
//! gradient fields are gated within stated bands, never bit-compared.
//!
//! Epistemic notice (binding): within-model numerical engineering on a
//! speculative theory's functional.  A passing gate certifies the port and
//! the backend plumbing, never anything about nature.

pub mod shaders;

use core::f64::consts::PI;

use fs_gum_field::{bps_floor, Field3, GHOST};
use fs_gum_statics::{Opts, Out, SGN6};

const FOURPI: f64 = 4.0 * PI;

/// Number of per-workgroup partial sums produced by the forward kernel.
pub const NSUMS: usize = 6;

fn ceil_div(a: u32, b: u32) -> u32 {
    a.div_ceil(b)
}

fn f64s_to_bytes(v: &[f64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(v.len() * 8);
    for x in v {
        out.extend_from_slice(&x.to_le_bytes());
    }
    out
}

fn bytes_to_f64s(b: &[u8]) -> Vec<f64> {
    b.chunks_exact(8).map(|c| f64::from_le_bytes(c.try_into().expect("8-byte chunk"))).collect()
}

/// The device-resident engine for one grid size.
pub struct GpuEngine {
    n: usize,
    p: usize,
    h: f64,
    adapter_info: wgpu::AdapterInfo,
    device: wgpu::Device,
    queue: wgpu::Queue,
    params: [f64; 16],
    nwg_fwd: u32,
    nwg_cell: u32,
    nwg_cell4: u32,
    nwg_p3: u32,
    // pipelines + bind groups (kernel inventory order)
    pl_fill_ghosts: wgpu::ComputePipeline,
    bg_fill_ghosts: wgpu::BindGroup,
    pl_forward: wgpu::ComputePipeline,
    bg_forward: wgpu::BindGroup,
    pl_flux: wgpu::ComputePipeline,
    bg_flux: wgpu::BindGroup,
    pl_gather: wgpu::ComputePipeline,
    bg_gather: wgpu::BindGroup,
    pl_renorm: wgpu::ComputePipeline,
    bg_renorm: wgpu::BindGroup,
    pl_axpy_v: wgpu::ComputePipeline,
    bg_axpy_v: wgpu::BindGroup,
    pl_step_q: wgpu::ComputePipeline,
    bg_step_q: wgpu::BindGroup,
    pl_project_v: wgpu::ComputePipeline,
    bg_project_v: wgpu::BindGroup,
    // buffers
    buf_f: [wgpu::Buffer; 4],
    buf_snap: [wgpu::Buffer; 4],
    buf_partials: wgpu::Buffer,
    buf_partials_stage: wgpu::Buffer,
    /// Kept for the handle inventory (referenced only via bind groups).
    _buf_flux: wgpu::Buffer,
    buf_grad: wgpu::Buffer,
    buf_grad_stage: wgpu::Buffer,
    buf_v: wgpu::Buffer,
    buf_params: wgpu::Buffer,
    buf_field_stage: wgpu::Buffer,
}

impl GpuEngine {
    /// Initialise the backend for an `n`-cell grid over `[-half, half]^3`:
    /// headless Vulkan instance (lavapipe in this container), SHADER_F64
    /// required, all pipelines compiled, all buffers allocated.
    pub fn new(n: usize, half: f64) -> Result<Self, String> {
        let dm = shaders::Dims::new(n);
        let p = dm.p as usize;
        let h = 2.0 * half / (n as f64);
        let ncell = dm.ncell();
        let p3 = dm.p3();
        let nwg_fwd = ceil_div(ncell, shaders::WG);

        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::VULKAN,
            ..wgpu::InstanceDescriptor::new_without_display_handle()
        });
        let adapter =
            pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                ..Default::default()
            }))
            .map_err(|e| format!("no Vulkan adapter: {e}"))?;
        if !adapter.features().contains(wgpu::Features::SHADER_F64) {
            return Err(format!(
                "adapter {:?} lacks SHADER_F64",
                adapter.get_info().name
            ));
        }
        let adapter_info = adapter.get_info();
        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
                label: Some("fs-gum-gpu"),
                required_features: wgpu::Features::SHADER_F64,
                ..Default::default()
            }))
            .map_err(|e| format!("request_device: {e}"))?;

        let mk_buf = |label: &str, len_f64: u64, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: len_f64 * 8,
                usage,
                mapped_at_creation: false,
            })
        };
        use wgpu::BufferUsages as U;
        let buf_f = core::array::from_fn(|c| {
            mk_buf(&format!("field-q{c}"), p3 as u64, U::STORAGE | U::COPY_DST | U::COPY_SRC)
        });
        let buf_snap = core::array::from_fn(|c| {
            mk_buf(&format!("snap-q{c}"), p3 as u64, U::COPY_DST | U::COPY_SRC)
        });
        let buf_partials =
            mk_buf("partials", u64::from(nwg_fwd) * NSUMS as u64, U::STORAGE | U::COPY_SRC);
        let buf_partials_stage =
            mk_buf("partials-stage", u64::from(nwg_fwd) * NSUMS as u64, U::MAP_READ | U::COPY_DST);
        let buf_flux = mk_buf("flux", u64::from(ncell) * 16, U::STORAGE);
        let buf_grad = mk_buf("grad", u64::from(ncell) * 4, U::STORAGE | U::COPY_SRC);
        let buf_grad_stage = mk_buf("grad-stage", u64::from(ncell) * 4, U::MAP_READ | U::COPY_DST);
        let buf_v = mk_buf("v", u64::from(ncell) * 4, U::STORAGE | U::COPY_DST);
        let buf_params = mk_buf("params", 16, U::STORAGE | U::COPY_DST);
        let buf_field_stage = mk_buf("field-stage", u64::from(p3) * 4, U::MAP_READ | U::COPY_DST);

        let mk_pl = |label: &str, src: String| {
            let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some(label),
                source: wgpu::ShaderSource::Wgsl(src.into()),
            });
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: None,
                module: &module,
                entry_point: Some("main"),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        let mk_bg = |label: &str, pl: &wgpu::ComputePipeline, bufs: &[&wgpu::Buffer]| {
            let entries: Vec<wgpu::BindGroupEntry> = bufs
                .iter()
                .enumerate()
                .map(|(i, b)| wgpu::BindGroupEntry {
                    binding: i as u32,
                    resource: b.as_entire_binding(),
                })
                .collect();
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some(label),
                layout: &pl.get_bind_group_layout(0),
                entries: &entries,
            })
        };

        let pl_fill_ghosts = mk_pl("fill_ghosts", shaders::fill_ghosts(dm));
        let bg_fill_ghosts = mk_bg(
            "fill_ghosts",
            &pl_fill_ghosts,
            &[&buf_f[0], &buf_f[1], &buf_f[2], &buf_f[3]],
        );
        let pl_forward = mk_pl("forward", shaders::forward(dm));
        let bg_forward = mk_bg(
            "forward",
            &pl_forward,
            &[&buf_f[0], &buf_f[1], &buf_f[2], &buf_f[3], &buf_partials, &buf_params],
        );
        let pl_flux = mk_pl("flux", shaders::flux(dm));
        let bg_flux = mk_bg(
            "flux",
            &pl_flux,
            &[&buf_f[0], &buf_f[1], &buf_f[2], &buf_f[3], &buf_flux, &buf_params],
        );
        let pl_gather = mk_pl("gather", shaders::gather(dm));
        let bg_gather = mk_bg(
            "gather",
            &pl_gather,
            &[&buf_f[0], &buf_f[1], &buf_f[2], &buf_f[3], &buf_flux, &buf_grad, &buf_params],
        );
        let pl_renorm = mk_pl("renormalize", shaders::renormalize(dm));
        let bg_renorm = mk_bg(
            "renormalize",
            &pl_renorm,
            &[&buf_f[0], &buf_f[1], &buf_f[2], &buf_f[3]],
        );
        let pl_axpy_v = mk_pl("axpy_v", shaders::axpy_v(dm));
        let bg_axpy_v = mk_bg("axpy_v", &pl_axpy_v, &[&buf_v, &buf_grad, &buf_params]);
        let pl_step_q = mk_pl("step_q", shaders::step_q(dm));
        let bg_step_q = mk_bg(
            "step_q",
            &pl_step_q,
            &[&buf_f[0], &buf_f[1], &buf_f[2], &buf_f[3], &buf_v, &buf_params],
        );
        let pl_project_v = mk_pl("project_v", shaders::project_v(dm));
        let bg_project_v = mk_bg(
            "project_v",
            &pl_project_v,
            &[&buf_f[0], &buf_f[1], &buf_f[2], &buf_f[3], &buf_v],
        );

        let h3 = h * h * h;
        let mut params = [0.0_f64; 16];
        params[6] = h3;
        params[9] = 8.0 / (12.0 * h); // c1
        params[10] = 1.0 / (12.0 * h); // c2s
        params[11] = 2.0 / FOURPI;
        let eng = GpuEngine {
            n,
            p,
            h,
            adapter_info,
            device,
            queue,
            params,
            nwg_fwd,
            nwg_cell: ceil_div(ncell, shaders::WG),
            nwg_cell4: ceil_div(4 * ncell, shaders::WG),
            nwg_p3: ceil_div(p3, shaders::WG),
            pl_fill_ghosts,
            bg_fill_ghosts,
            pl_forward,
            bg_forward,
            pl_flux,
            bg_flux,
            pl_gather,
            bg_gather,
            pl_renorm,
            bg_renorm,
            pl_axpy_v,
            bg_axpy_v,
            pl_step_q,
            bg_step_q,
            pl_project_v,
            bg_project_v,
            buf_f,
            buf_snap,
            buf_partials,
            buf_partials_stage,
            _buf_flux: buf_flux,
            buf_grad,
            buf_grad_stage,
            buf_v,
            buf_params,
            buf_field_stage,
        };
        eng.flush_params();
        Ok(eng)
    }

    /// Adapter identity line (environment proof for RESULTS.md).
    #[must_use]
    pub fn adapter_desc(&self) -> String {
        format!(
            "name={:?} backend={:?} type={:?} driver={:?} driver_info={:?}",
            self.adapter_info.name,
            self.adapter_info.backend,
            self.adapter_info.device_type,
            self.adapter_info.driver,
            self.adapter_info.driver_info
        )
    }

    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }

    #[must_use]
    pub fn h(&self) -> f64 {
        self.h
    }

    fn flush_params(&self) {
        self.queue.write_buffer(&self.buf_params, 0, &f64s_to_bytes(&self.params));
    }

    fn dispatch(&self, pl: &wgpu::ComputePipeline, bg: &wgpu::BindGroup, nwg: u32) {
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(pl);
            pass.set_bind_group(0, bg, &[]);
            pass.dispatch_workgroups(nwg, 1, 1);
        }
        self.queue.submit([enc.finish()]);
    }

    fn read_stage(&self, src: &wgpu::Buffer, stage: &wgpu::Buffer, bytes: u64) -> Vec<f64> {
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.copy_buffer_to_buffer(src, 0, stage, 0, bytes);
        self.queue.submit([enc.finish()]);
        let slice = stage.slice(0..bytes);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("map_async"));
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("poll");
        let out = {
            let view = slice.get_mapped_range().expect("mapped range");
            bytes_to_f64s(&view)
        };
        stage.unmap();
        out
    }

    // ---- field transfer ---------------------------------------------------

    /// Upload a field: the full padded per-component arrays, ghosts included.
    pub fn upload(&self, f: &Field3) {
        assert_eq!(f.n(), self.n, "grid mismatch");
        let p = self.p;
        let mut comp = vec![0.0_f64; p * p * p];
        for (c, buf) in self.buf_f.iter().enumerate() {
            let mut o = 0usize;
            for ip in 0..p {
                for jp in 0..p {
                    for kp in 0..p {
                        comp[o] = f.get_pad(c, ip, jp, kp);
                        o += 1;
                    }
                }
            }
            self.queue.write_buffer(buf, 0, &f64s_to_bytes(&comp));
        }
    }

    /// Upload raw per-component padded data (4 * P^3, component-major) —
    /// used by the fill_ghosts gate to plant corrupted ghosts.
    pub fn upload_raw(&self, data: &[f64]) {
        let p3 = self.p * self.p * self.p;
        assert_eq!(data.len(), 4 * p3, "raw field size");
        for (c, buf) in self.buf_f.iter().enumerate() {
            self.queue.write_buffer(buf, 0, &f64s_to_bytes(&data[c * p3..(c + 1) * p3]));
        }
    }

    /// Download the raw padded field (4 * P^3, component-major).
    #[must_use]
    pub fn download_raw(&self) -> Vec<f64> {
        let p3 = (self.p * self.p * self.p) as u64;
        let mut enc = self.device.create_command_encoder(&Default::default());
        for (c, buf) in self.buf_f.iter().enumerate() {
            enc.copy_buffer_to_buffer(buf, 0, &self.buf_field_stage, c as u64 * p3 * 8, p3 * 8);
        }
        self.queue.submit([enc.finish()]);
        let bytes = 4 * p3 * 8;
        let slice = self.buf_field_stage.slice(0..bytes);
        slice.map_async(wgpu::MapMode::Read, |r| r.expect("map_async"));
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("poll");
        let out = {
            let view = slice.get_mapped_range().expect("mapped range");
            bytes_to_f64s(&view)
        };
        self.buf_field_stage.unmap();
        out
    }

    /// Download the REAL cells into a Field3 (ghosts of `f` untouched).
    pub fn download_field(&self, f: &mut Field3) {
        assert_eq!(f.n(), self.n, "grid mismatch");
        let raw = self.download_raw();
        let (n, p) = (self.n, self.p);
        let p3 = p * p * p;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let o = ((i + GHOST) * p + (j + GHOST)) * p + (k + GHOST);
                    f.set(i, j, k, [raw[o], raw[p3 + o], raw[2 * p3 + o], raw[3 * p3 + o]]);
                }
            }
        }
    }

    // ---- kernels ----------------------------------------------------------

    /// (a) reset the ghost rind to vacuum on-device.
    pub fn fill_ghosts(&self) {
        self.dispatch(&self.pl_fill_ghosts, &self.bg_fill_ghosts, self.nwg_p3);
    }

    /// (d) pointwise unit-norm renormalisation of the real cells on-device.
    pub fn renormalize(&self) {
        self.dispatch(&self.pl_renorm, &self.bg_renorm, self.nwg_cell);
    }

    /// (b) forward sweep; returns the raw per-workgroup partials
    /// (nwg * 6 f64) — the replay gate compares these bitwise.
    #[must_use]
    pub fn forward_partials(&self) -> Vec<f64> {
        self.dispatch(&self.pl_forward, &self.bg_forward, self.nwg_fwd);
        self.read_stage(
            &self.buf_partials,
            &self.buf_partials_stage,
            u64::from(self.nwg_fwd) * NSUMS as u64 * 8,
        )
    }

    /// Fixed-order CPU combine of the workgroup partials: ascending
    /// workgroup index per quantity.  Returns
    /// (s_e2, s_e4, s_det, s_det2, s_e0, s_i).
    #[must_use]
    pub fn combine_partials(partials: &[f64]) -> [f64; NSUMS] {
        let mut sums = [0.0_f64; NSUMS];
        for (s, sum) in sums.iter_mut().enumerate() {
            let mut acc = 0.0_f64;
            for wg in 0..partials.len() / NSUMS {
                acc += partials[wg * NSUMS + s];
            }
            *sum = acc;
        }
        sums
    }

    /// Forward objective evaluation: device sweep + fixed-order CPU combine
    /// + the exact `fs_gum_statics::engine::eval` forward bookkeeping
    /// (sector normalisations, Routhian, one-sided guards, objective).
    #[must_use]
    pub fn forward_out(&self, o: &Opts) -> Out {
        let sums = Self::combine_partials(&self.forward_partials());
        self.out_from_sums(o, &sums)
    }

    /// Objective bookkeeping on already-read-back partials (lets one
    /// forward dispatch serve several objective specifications — the
    /// forward sweep itself is objective-independent).
    #[must_use]
    pub fn out_from_partials(&self, o: &Opts, partials: &[f64]) -> Out {
        self.out_from_sums(o, &Self::combine_partials(partials))
    }

    fn out_from_sums(&self, o: &Opts, s: &[f64; NSUMS]) -> Out {
        let h3 = self.h * self.h * self.h;
        let e2 = s[0] * h3 / FOURPI;
        let e4 = s[1] * h3 / FOURPI;
        let e6 = s[3] * h3 / FOURPI;
        let e0 = s[4] * h3 / FOURPI;
        let i_val = 2.0 * s[5] * h3;
        let deg = SGN6 * s[2] * h3 / (2.0 * PI * PI);
        let estat = o.t * (e2 + e4) + e6 + e0;
        let rot = match o.l {
            Some(l) => l * l / (2.0 * i_val),
            None => 0.0,
        };
        let r = estat + rot;
        let floor_gap = e6 + e0 - bps_floor() * deg / o.deg_ref;
        let epen = match o.anchor {
            Some(a) => {
                let dev = (deg - (o.deg_ref - a.band)).min(0.0);
                0.5 * a.mu * dev * dev
            }
            None => 0.0,
        };
        let efpen = match o.wall {
            Some(w) => {
                let devf = (floor_gap - w.fgap_ref).min(0.0);
                0.5 * w.mu * devf * devf
            }
            None => 0.0,
        };
        let obj = o.c[0] * e2 + o.c[1] * e4 + o.c[2] * e6 + o.c[3] * e0 + rot + epen + efpen;
        Out { e2, e4, e6, e0, i: i_val, deg, estat, r, floor_gap, epen, efpen, obj }
    }

    /// Gradient pass for the objective state `out` (obtained from
    /// [`forward_out`] on the SAME field state): computes the guard-dressed
    /// kernel weights on the CPU (the ~10-scalar params upload), then
    /// dispatches flux + gather.  The gradient stays device-resident in the
    /// Grad-layout buffer.
    pub fn gradient(&mut self, o: &Opts, out: &Out) {
        // guard-dressed weights, verbatim from engine::eval pass 2
        let dev = match o.anchor {
            Some(a) => (out.deg - (o.deg_ref - a.band)).min(0.0),
            None => 0.0,
        };
        let (devf, wall_mu) = match o.wall {
            Some(w) => ((out.floor_gap - w.fgap_ref).min(0.0), w.mu),
            None => (0.0, 0.0),
        };
        let w6s = (o.c[2] + wall_mu * devf) / (2.0 * PI);
        let mut cst = 0.0_f64;
        if dev < 0.0 {
            cst += o.anchor.map_or(0.0, |a| a.mu) * dev;
        }
        if devf < 0.0 {
            cst -= wall_mu * devf * bps_floor() / o.deg_ref;
        }
        let cstw = cst * SGN6 / (2.0 * PI * PI);
        let e0c = -(o.c[3] + wall_mu * devf) / FOURPI;
        let rotfac = match o.l {
            Some(l) => -(l * l) / (2.0 * out.i * out.i) * 4.0,
            None => 0.0,
        };
        self.params[0] = o.c[0];
        self.params[1] = o.c[1];
        self.params[2] = w6s;
        self.params[3] = cstw;
        self.params[4] = e0c;
        self.params[5] = rotfac;
        self.flush_params();
        let mut enc = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.pl_flux);
            pass.set_bind_group(0, &self.bg_flux, &[]);
            pass.dispatch_workgroups(self.nwg_cell, 1, 1);
        }
        {
            let mut pass = enc.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.pl_gather);
            pass.set_bind_group(0, &self.bg_gather, &[]);
            pass.dispatch_workgroups(self.nwg_cell, 1, 1);
        }
        self.queue.submit([enc.finish()]);
    }

    /// Download the gradient (Grad layout, 4 * N^3).
    #[must_use]
    pub fn download_grad(&self) -> Vec<f64> {
        self.read_stage(
            &self.buf_grad,
            &self.buf_grad_stage,
            (4 * self.n * self.n * self.n) as u64 * 8,
        )
    }

    // ---- ANF device ops ---------------------------------------------------

    /// v = 0 (device-side clear).
    pub fn v_zero(&self) {
        let mut enc = self.device.create_command_encoder(&Default::default());
        enc.clear_buffer(&self.buf_v, 0, None);
        self.queue.submit([enc.finish()]);
    }

    /// v -= (dt/h^3) g.
    pub fn v_axpy(&mut self, acoef: f64) {
        self.params[7] = acoef;
        self.flush_params();
        self.dispatch(&self.pl_axpy_v, &self.bg_axpy_v, self.nwg_cell4);
    }

    /// q += dt v on the real cells.
    pub fn q_step(&mut self, dt: f64) {
        self.params[8] = dt;
        self.flush_params();
        self.dispatch(&self.pl_step_q, &self.bg_step_q, self.nwg_cell);
    }

    /// v -= (v.q) q per real cell.
    pub fn v_project(&self) {
        self.dispatch(&self.pl_project_v, &self.bg_project_v, self.nwg_cell);
    }

    /// Snapshot the field into the device-side backup buffers.
    pub fn snapshot_save(&self) {
        let p3 = (self.p * self.p * self.p) as u64 * 8;
        let mut enc = self.device.create_command_encoder(&Default::default());
        for c in 0..4 {
            enc.copy_buffer_to_buffer(&self.buf_f[c], 0, &self.buf_snap[c], 0, p3);
        }
        self.queue.submit([enc.finish()]);
    }

    /// Restore the field from the device-side backup buffers.
    pub fn snapshot_restore(&self) {
        let p3 = (self.p * self.p * self.p) as u64 * 8;
        let mut enc = self.device.create_command_encoder(&Default::default());
        for c in 0..4 {
            enc.copy_buffer_to_buffer(&self.buf_snap[c], 0, &self.buf_f[c], 0, p3);
        }
        self.queue.submit([enc.finish()]);
    }
}

// ---------------------------------------------------------------------------
// device-resident arrested Newton flow (the ANF smoke driver)
// ---------------------------------------------------------------------------

/// Result of a device-resident ANF run.
pub struct GpuAnfResult {
    /// Last ACCEPTED objective state (the field on device holds it).
    pub out: Out,
    /// Accepted objective after every iteration (index 0 = initial state).
    pub objs: Vec<f64>,
    pub iters: usize,
    pub arrests: usize,
    pub status: &'static str,
}

/// Arrested Newton flow with the field, gradient, velocity and snapshot all
/// device-resident; the CPU sees only the per-iteration objective readback
/// (the forward partials, ~KBs) and runs the frozen field3d arrest logic
/// (revert + zero velocity + dt *= 0.6 on any objective increase; dt *=
/// 1.01 capped at dt_max on acceptance; stall below 1e-7).  The field must
/// already be uploaded with ghosts filled.
pub fn anf_gpu(e: &mut GpuEngine, o: &Opts, maxit: usize, dt0: f64, dt_max: f64) -> GpuAnfResult {
    let h = e.h();
    let h3 = h * h * h;
    let mut out = e.forward_out(o);
    e.gradient(o, &out);
    e.snapshot_save();
    e.v_zero();
    let mut dt = dt0;
    let mut arrests = 0usize;
    let mut objs = vec![out.obj];
    let mut status = "iter_cap";
    let mut it = 0usize;
    while it < maxit {
        it += 1;
        e.v_axpy(dt / h3);
        e.q_step(dt);
        e.renormalize();
        let out_new = e.forward_out(o);
        if out_new.obj > out.obj {
            e.snapshot_restore();
            e.v_zero();
            dt *= 0.6;
            arrests += 1;
            if dt < 1.0e-7 {
                objs.push(out.obj);
                status = "stall_dt";
                break;
            }
            // gradient of the accepted point is still resident (kept)
        } else {
            out = out_new;
            e.gradient(o, &out);
            e.snapshot_save();
            e.v_project();
            dt = (dt * 1.01).min(dt_max);
        }
        objs.push(out.obj);
    }
    // leave the field at the last accepted point
    e.snapshot_restore();
    GpuAnfResult { out, objs, iters: it, arrests, status }
}
