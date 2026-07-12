//! WGSL f64 compute kernels for the central4 gum sweep, generated as
//! strings with the grid extents baked in (one shader set per N).
//!
//! f64 in WGSL is the naga extension empirically validated by the Phase G-A
//! survey (scale_survey_v4.json, DIMENSION 3): `f64` type + `lf` literals,
//! gated by `naga::valid::Capabilities::FLOAT64` = wgpu
//! `Features::SHADER_F64`.  All h-derived and guard-dressed constants are
//! passed through a small f64 storage buffer (`params`) rather than as WGSL
//! literals so their bit patterns are exactly the CPU engine's.
//!
//! Determinism contract (device-side): every kernel writes a disjoint
//! output location per invocation; the only cross-invocation reduction is
//! the forward kernel's per-workgroup partial sum, which lane 0 accumulates
//! over the 64 lanes in FIXED ascending lane order after a barrier.  The
//! final combine over workgroup partials happens on the CPU in fixed
//! ascending workgroup order (`lib.rs`).  Same-device replay of any
//! dispatch is therefore bit-identical (gated); CPU-vs-GPU agreement is the
//! tolerance-band class (llvmpipe/real drivers may FMA-contract — observed
//! here at the 1-2 ulp per-operation level).
//!
//! Params buffer layout (f64, 16 slots):
//!   0 c2w   objective weight on E2 (c[0])
//!   1 c4w   objective weight on E4 (c[1])
//!   2 w6s   guard-dressed det weight (c[2] + MU_F devf) / (2 pi)
//!   3 cstw  guard constant (MU dev - MU_F devf FLOOR/deg_ref) SGN6/(2 pi^2)
//!   4 e0c   E0 density coefficient -(c[3] + MU_F devf) / (4 pi)
//!   5 rotfac  Routhian cell-term factor -(L^2/2I^2) 4  (0 for statics)
//!   6 h3    cell volume
//!   7 acoef ANF velocity kick dt/h^3
//!   8 dt    ANF step
//!   9 c1    4th-order stencil coefficient 8/(12h)
//!  10 c2s   4th-order stencil coefficient 1/(12h)
//!  11 two4pi  2/(4 pi)
//!  12-15 reserved (0)

/// Fixed workgroup size for every kernel (part of the bit contract).
pub const WG: u32 = 64;

/// The 6-term det pair expansion (field3d PAIRS; must match
/// `fs_gum_statics::engine`).
const PAIRS: [(usize, usize, usize, usize, f64); 6] = [
    (0, 1, 2, 3, 1.0),
    (0, 2, 1, 3, -1.0),
    (0, 3, 1, 2, 1.0),
    (1, 2, 0, 3, 1.0),
    (1, 3, 0, 2, -1.0),
    (2, 3, 0, 1, 1.0),
];

fn sg_lit(sg: f64) -> &'static str {
    if sg > 0.0 {
        "1.0lf"
    } else {
        "-1.0lf"
    }
}

/// Grid extents baked into every shader.
#[derive(Clone, Copy)]
pub struct Dims {
    pub n: u32,
    pub p: u32,
}

impl Dims {
    pub fn new(n: usize) -> Self {
        let n = u32::try_from(n).expect("n fits u32");
        Dims { n, p: n + 2 * fs_gum_field::GHOST as u32 }
    }
    pub fn ncell(&self) -> u32 {
        self.n * self.n * self.n
    }
    pub fn p3(&self) -> u32 {
        self.p * self.p * self.p
    }
    fn consts(&self) -> String {
        format!(
            "const NN: u32 = {n}u;\nconst PP: u32 = {p}u;\nconst NCELL: u32 = {nc}u;\nconst P3: u32 = {p3}u;\nconst NCELL4: u32 = {nc4}u;\nconst NG: u32 = 2u;\n\nfn offp(ip: u32, jp: u32, kp: u32) -> u32 {{\n    return (ip * PP + jp) * PP + kp;\n}}\n",
            n = self.n,
            p = self.p,
            nc = self.ncell(),
            p3 = self.p3(),
            nc4 = 4 * self.ncell(),
        )
    }
}

/// Field bindings (4 per-component f64 storage buffers, matching
/// fs-gum-field's SoA layout: component c at padded offset
/// `(ip*P + jp)*P + kp` of buffer `f{c}`).
fn field_bindings(read_write: bool) -> String {
    let access = if read_write { "read_write" } else { "read" };
    let mut s = String::new();
    for c in 0..4 {
        s.push_str(&format!(
            "@group(0) @binding({c}) var<storage, {access}> f{c}: array<f64>;\n"
        ));
    }
    s
}

/// Per-component loads of q and the 4th-order central derivatives at real
/// cell (i, j, k) (padded ip = i + NG, ...).  Emits `q[c]` and `d[ax][c]`,
/// same expression order as `fs_gum_field::stencil::derivs4`.
fn loads_q_d() -> String {
    let mut s = String::new();
    s.push_str("    var q: array<f64, 4>;\n    var d: array<array<f64, 4>, 3>;\n");
    for c in 0..4 {
        s.push_str(&format!(
            "    q[{c}] = f{c}[offp(ip, jp, kp)];\n\
             \x20   d[0][{c}] = (f{c}[offp(ip + 1u, jp, kp)] - f{c}[offp(ip - 1u, jp, kp)]) * c1\n\
             \x20       - (f{c}[offp(ip + 2u, jp, kp)] - f{c}[offp(ip - 2u, jp, kp)]) * c2s;\n\
             \x20   d[1][{c}] = (f{c}[offp(ip, jp + 1u, kp)] - f{c}[offp(ip, jp - 1u, kp)]) * c1\n\
             \x20       - (f{c}[offp(ip, jp + 2u, kp)] - f{c}[offp(ip, jp - 2u, kp)]) * c2s;\n\
             \x20   d[2][{c}] = (f{c}[offp(ip, jp, kp + 1u)] - f{c}[offp(ip, jp, kp - 1u)]) * c1\n\
             \x20       - (f{c}[offp(ip, jp, kp + 2u)] - f{c}[offp(ip, jp, kp - 2u)]) * c2s;\n"
        ));
    }
    s
}

/// nn / dd dot products (fs_gum_statics::engine::nn_dd expression order).
fn nn_dd_code() -> String {
    let mut s = String::new();
    for ax in 0..3 {
        s.push_str(&format!(
            "    let nn{ax} = d[{ax}][0] * d[{ax}][0] + d[{ax}][1] * d[{ax}][1] + d[{ax}][2] * d[{ax}][2] + d[{ax}][3] * d[{ax}][3];\n"
        ));
    }
    for (idx, (i, j)) in [(0usize, 1usize), (0, 2), (1, 2)].iter().enumerate() {
        s.push_str(&format!(
            "    let dd{idx} = d[{i}][0] * d[{j}][0] + d[{i}][1] * d[{j}][1] + d[{i}][2] * d[{j}][2] + d[{i}][3] * d[{j}][3];\n"
        ));
    }
    s
}

/// det[q, Dx, Dy, Dz] via the 6-term PAIRS expansion, fixed term order
/// (fs_gum_statics::engine::det_pairs).
fn det_pairs_code() -> String {
    let mut s = String::from("    var det = 0.0lf;\n");
    for &(a, b, c, e, sg) in &PAIRS {
        s.push_str(&format!(
            "    {{\n        let m = q[{a}] * d[0][{b}] - q[{b}] * d[0][{a}];\n        let pc = d[1][{c}] * d[2][{e}] - d[1][{e}] * d[2][{c}];\n        det = det + {sg} * (m * pc);\n    }}\n",
            sg = sg_lit(sg)
        ));
    }
    s
}

/// (a) fill_ghosts: reset every ghost cell of the padded rind to the vacuum
/// (1, 0, 0, 0); real cells untouched.
pub fn fill_ghosts(dm: Dims) -> String {
    format!(
        "{consts}{fields}\n@compute @workgroup_size({WG})\nfn main(@builtin(global_invocation_id) gid: vec3<u32>) {{\n    let t = gid.x;\n    if (t >= P3) {{\n        return;\n    }}\n    let ip = t / (PP * PP);\n    let jp = (t / PP) % PP;\n    let kp = t % PP;\n    let real = ip >= NG && ip < NG + NN && jp >= NG && jp < NG + NN && kp >= NG && kp < NG + NN;\n    if (!real) {{\n        f0[t] = 1.0lf;\n        f1[t] = 0.0lf;\n        f2[t] = 0.0lf;\n        f3[t] = 0.0lf;\n    }}\n}}\n",
        consts = dm.consts(),
        fields = field_bindings(true),
    )
}

/// (b) central4 forward sweep: one thread per interior cell, per-workgroup
/// partial sums for (s_e2, s_e4, s_det, s_det2, s_e0, s_i) via a
/// workgroup-shared reduction in FIXED ascending lane order; partials to
/// `partials[wg*6 + s]`.  Final combine on the CPU in fixed workgroup order.
pub fn forward(dm: Dims) -> String {
    let mut sh_decl = String::new();
    let mut sh_store = String::new();
    let mut acc_decl = String::new();
    let mut acc_loop = String::new();
    let mut acc_out = String::new();
    for s in 0..6 {
        sh_decl.push_str(&format!("var<workgroup> sh{s}: array<f64, {WG}>;\n"));
        sh_store.push_str(&format!("    sh{s}[lid] = s{s};\n"));
        acc_decl.push_str(&format!("        var t{s} = 0.0lf;\n"));
        acc_loop.push_str(&format!("            t{s} = t{s} + sh{s}[t];\n"));
        acc_out.push_str(&format!("        partials[o + {s}u] = t{s};\n"));
    }
    format!(
        "{consts}{fields}@group(0) @binding(4) var<storage, read_write> partials: array<f64>;\n@group(0) @binding(5) var<storage, read> params: array<f64>;\n\n{sh_decl}\n@compute @workgroup_size({WG})\nfn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lv: vec3<u32>) {{\n    let lid = lv.x;\n    let cell = wid.x * {WG}u + lid;\n    var s0 = 0.0lf;\n    var s1 = 0.0lf;\n    var s2 = 0.0lf;\n    var s3 = 0.0lf;\n    var s4 = 0.0lf;\n    var s5 = 0.0lf;\n    if (cell < NCELL) {{\n    let i = cell / (NN * NN);\n    let j = (cell / NN) % NN;\n    let k = cell % NN;\n    let ip = i + NG;\n    let jp = j + NG;\n    let kp = k + NG;\n    let c1 = params[9];\n    let c2s = params[10];\n{loads}{nndd}{det}    s0 = nn0 + nn1 + nn2;\n    s1 = (nn0 * nn1 - dd0 * dd0) + (nn0 * nn2 - dd1 * dd1) + (nn1 * nn2 - dd2 * dd2);\n    s2 = det;\n    s3 = det * det;\n    s4 = 1.0lf - q[0];\n    s5 = q[1] * q[1] + q[2] * q[2];\n    }}\n{sh_store}    workgroupBarrier();\n    if (lid == 0u) {{\n{acc_decl}        for (var t = 0u; t < {WG}u; t = t + 1u) {{\n{acc_loop}        }}\n        let o = wid.x * 6u;\n{acc_out}    }}\n}}\n",
        consts = dm.consts(),
        fields = field_bindings(false),
        loads = loads_q_d(),
        nndd = nn_dd_code(),
        det = det_pairs_code(),
    )
}

/// (c1) pointwise flux pass: per interior cell the (2+4) flux P[3][4] and
/// the sextic/guard q-derivative gq[4] (fs_gum_statics::engine::point_flux,
/// identical expression order), written to the intermediate buffer as
/// 16 f64 per cell: [P[0][0..4], P[1][0..4], P[2][0..4], gq[0..4]].
pub fn flux(dm: Dims) -> String {
    // OTHERS table of engine::point_flux: per axis (o1, o2, dd idx, dd idx)
    const OTHERS: [(usize, usize, usize, usize); 3] = [(1, 2, 0, 1), (0, 2, 0, 2), (0, 1, 1, 2)];
    let mut paxes = String::new();
    for ax in 0..3 {
        let (o1, o2, k1, k2) = OTHERS[ax];
        paxes.push_str(&format!(
            "    {{\n        let t1 = c2w + c4w * (nn{o1} + nn{o2});\n        for (var a = 0u; a < 4u; a = a + 1u) {{\n            p[{ax}][a] = two4pi * (t1 * d[{ax}][a] - c4w * (dd{k1} * d[{o1}][a] + dd{k2} * d[{o2}][a]));\n        }}\n    }}\n"
        ));
    }
    let mut pairs6 = String::new();
    for &(a, b, c, e, sg) in &PAIRS {
        let sg = sg_lit(sg);
        pairs6.push_str(&format!(
            "    {{\n        let m = q[{a}] * d[0][{b}] - q[{b}] * d[0][{a}];\n        let pc = d[1][{c}] * d[2][{e}] - d[1][{e}] * d[2][{c}];\n        let t2 = {sg} * w6 * pc;\n        gq[{a}] = gq[{a}] + t2 * d[0][{b}];\n        gq[{b}] = gq[{b}] - t2 * d[0][{a}];\n        p[0][{b}] = p[0][{b}] + t2 * q[{a}];\n        p[0][{a}] = p[0][{a}] - t2 * q[{b}];\n        let t3 = {sg} * w6 * m;\n        p[1][{c}] = p[1][{c}] + t3 * d[2][{e}];\n        p[1][{e}] = p[1][{e}] - t3 * d[2][{c}];\n        p[2][{e}] = p[2][{e}] + t3 * d[1][{c}];\n        p[2][{c}] = p[2][{c}] - t3 * d[1][{e}];\n    }}\n"
        ));
    }
    format!(
        "{consts}{fields}@group(0) @binding(4) var<storage, read_write> flux: array<f64>;\n@group(0) @binding(5) var<storage, read> params: array<f64>;\n\n@compute @workgroup_size({WG})\nfn main(@builtin(global_invocation_id) gid: vec3<u32>) {{\n    let cell = gid.x;\n    if (cell >= NCELL) {{\n        return;\n    }}\n    let i = cell / (NN * NN);\n    let j = (cell / NN) % NN;\n    let k = cell % NN;\n    let ip = i + NG;\n    let jp = j + NG;\n    let kp = k + NG;\n    let c1 = params[9];\n    let c2s = params[10];\n    let c2w = params[0];\n    let c4w = params[1];\n    let w6s = params[2];\n    let cstw = params[3];\n    let two4pi = params[11];\n{loads}{nndd}    var p: array<array<f64, 4>, 3>;\n{paxes}{det}    let w6 = w6s * det + cstw;\n    var gq: array<f64, 4>;\n    for (var a = 0u; a < 4u; a = a + 1u) {{\n        gq[a] = 0.0lf;\n    }}\n{pairs6}    let base = cell * 16u;\n    for (var ax = 0u; ax < 3u; ax = ax + 1u) {{\n        for (var a = 0u; a < 4u; a = a + 1u) {{\n            flux[base + ax * 4u + a] = p[ax][a];\n        }}\n    }}\n    for (var a = 0u; a < 4u; a = a + 1u) {{\n        flux[base + 12u + a] = gq[a];\n    }}\n}}\n",
        consts = dm.consts(),
        fields = field_bindings(false),
        loads = loads_q_d(),
        nndd = nn_dd_code(),
        det = det_pairs_code(),
    )
}

/// (c2) gather adjoint: per REAL cell, gather the transposed central4
/// stencil contributions from the eval points that read it (ghost flux
/// dropped = the boundary conditions of the scatter form), add the cell
/// density terms, apply the h^3 weight and the pointwise tangent
/// projection, and write the gradient (Grad layout a*NCELL + cell).
///
/// FIXED per-cell add order (the gather contract): gq(self); then per axis
/// x, y, z the neighbour fluxes at offsets (-1, +1, -2, +2); then the cell
/// terms.  This is a deliberate reordering of the CPU scatter's
/// eval-point-ascending add order — CPU-vs-GPU gradient agreement is a
/// tolerance-band gate, never a bit gate.
pub fn gather(dm: Dims) -> String {
    let mut qload = String::new();
    for c in 0..4 {
        qload.push_str(&format!("    qv[{c}] = f{c}[offp(ip, jp, kp)];\n"));
    }
    format!(
        "{consts}{fields}@group(0) @binding(4) var<storage, read> flux: array<f64>;\n@group(0) @binding(5) var<storage, read_write> grad: array<f64>;\n@group(0) @binding(6) var<storage, read> params: array<f64>;\n\nfn cbase(i: u32, j: u32, k: u32) -> u32 {{\n    return ((i * NN + j) * NN + k) * 16u;\n}}\n\n@compute @workgroup_size({WG})\nfn main(@builtin(global_invocation_id) gid: vec3<u32>) {{\n    let cell = gid.x;\n    if (cell >= NCELL) {{\n        return;\n    }}\n    let i = cell / (NN * NN);\n    let j = (cell / NN) % NN;\n    let k = cell % NN;\n    let ip = i + NG;\n    let jp = j + NG;\n    let kp = k + NG;\n    let c1 = params[9];\n    let c2s = params[10];\n    let e0c = params[4];\n    let rotfac = params[5];\n    let h3 = params[6];\n    var qv: array<f64, 4>;\n{qload}    var g: array<f64, 4>;\n    let self16 = cbase(i, j, k);\n    for (var a = 0u; a < 4u; a = a + 1u) {{\n        var gg = flux[self16 + 12u + a];\n        // axis x (P slot 0): transpose of the scatter with ghost flux dropped\n        if (i >= 1u) {{\n            gg = gg + c1 * flux[cbase(i - 1u, j, k) + a];\n        }}\n        if (i + 1u < NN) {{\n            gg = gg - c1 * flux[cbase(i + 1u, j, k) + a];\n        }}\n        if (i >= 2u) {{\n            gg = gg - c2s * flux[cbase(i - 2u, j, k) + a];\n        }}\n        if (i + 2u < NN) {{\n            gg = gg + c2s * flux[cbase(i + 2u, j, k) + a];\n        }}\n        // axis y (P slot 1)\n        if (j >= 1u) {{\n            gg = gg + c1 * flux[cbase(i, j - 1u, k) + 4u + a];\n        }}\n        if (j + 1u < NN) {{\n            gg = gg - c1 * flux[cbase(i, j + 1u, k) + 4u + a];\n        }}\n        if (j >= 2u) {{\n            gg = gg - c2s * flux[cbase(i, j - 2u, k) + 4u + a];\n        }}\n        if (j + 2u < NN) {{\n            gg = gg + c2s * flux[cbase(i, j + 2u, k) + 4u + a];\n        }}\n        // axis z (P slot 2)\n        if (k >= 1u) {{\n            gg = gg + c1 * flux[cbase(i, j, k - 1u) + 8u + a];\n        }}\n        if (k + 1u < NN) {{\n            gg = gg - c1 * flux[cbase(i, j, k + 1u) + 8u + a];\n        }}\n        if (k >= 2u) {{\n            gg = gg - c2s * flux[cbase(i, j, k - 2u) + 8u + a];\n        }}\n        if (k + 2u < NN) {{\n            gg = gg + c2s * flux[cbase(i, j, k + 2u) + 8u + a];\n        }}\n        g[a] = gg;\n    }}\n    // cell density terms (E0 + Routhian), h^3 weight, tangent projection\n    g[0] = g[0] + e0c;\n    g[1] = g[1] + rotfac * qv[1];\n    g[2] = g[2] + rotfac * qv[2];\n    for (var a = 0u; a < 4u; a = a + 1u) {{\n        g[a] = g[a] * h3;\n    }}\n    let dt = g[0] * qv[0] + g[1] * qv[1] + g[2] * qv[2] + g[3] * qv[3];\n    for (var a = 0u; a < 4u; a = a + 1u) {{\n        grad[a * NCELL + cell] = g[a] - dt * qv[a];\n    }}\n}}\n",
        consts = dm.consts(),
        fields = field_bindings(false),
    )
}

/// (d) renormalize: pointwise unit-norm fix of every REAL cell (ghosts
/// untouched), IEEE sqrt, same expression order as Field3::renormalize.
pub fn renormalize(dm: Dims) -> String {
    let mut load = String::new();
    let mut store = String::new();
    for c in 0..4 {
        load.push_str(&format!("    let q{c} = f{c}[o];\n"));
        store.push_str(&format!("    f{c}[o] = q{c} / nrm;\n"));
    }
    format!(
        "{consts}{fields}\n@compute @workgroup_size({WG})\nfn main(@builtin(global_invocation_id) gid: vec3<u32>) {{\n    let cell = gid.x;\n    if (cell >= NCELL) {{\n        return;\n    }}\n    let i = cell / (NN * NN);\n    let j = (cell / NN) % NN;\n    let k = cell % NN;\n    let o = offp(i + NG, j + NG, k + NG);\n{load}    let nrm = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);\n{store}}}\n",
        consts = dm.consts(),
        fields = field_bindings(true),
    )
}

/// ANF: velocity kick v -= (dt/h^3) g over the flat 4*N^3 buffer.
pub fn axpy_v(dm: Dims) -> String {
    format!(
        "{consts}@group(0) @binding(0) var<storage, read_write> v: array<f64>;\n@group(0) @binding(1) var<storage, read> grad: array<f64>;\n@group(0) @binding(2) var<storage, read> params: array<f64>;\n\n@compute @workgroup_size({WG})\nfn main(@builtin(global_invocation_id) gid: vec3<u32>) {{\n    let t = gid.x;\n    if (t >= NCELL4) {{\n        return;\n    }}\n    v[t] = v[t] - params[7] * grad[t];\n}}\n",
        consts = dm.consts(),
    )
}

/// ANF: field step q += dt v on the REAL cells (v in Grad layout).
pub fn step_q(dm: Dims) -> String {
    let mut upd = String::new();
    for c in 0..4 {
        upd.push_str(&format!(
            "    f{c}[o] = f{c}[o] + dt * v[{c}u * NCELL + cell];\n"
        ));
    }
    format!(
        "{consts}{fields}@group(0) @binding(4) var<storage, read> v: array<f64>;\n@group(0) @binding(5) var<storage, read> params: array<f64>;\n\n@compute @workgroup_size({WG})\nfn main(@builtin(global_invocation_id) gid: vec3<u32>) {{\n    let cell = gid.x;\n    if (cell >= NCELL) {{\n        return;\n    }}\n    let i = cell / (NN * NN);\n    let j = (cell / NN) % NN;\n    let k = cell % NN;\n    let o = offp(i + NG, j + NG, k + NG);\n    let dt = params[8];\n{upd}}}\n",
        consts = dm.consts(),
        fields = field_bindings(true),
    )
}

/// ANF: velocity tangent re-projection v -= (v.q) q per REAL cell.
pub fn project_v(dm: Dims) -> String {
    let mut qload = String::new();
    for c in 0..4 {
        qload.push_str(&format!("    let q{c} = f{c}[o];\n"));
    }
    format!(
        "{consts}{fields}@group(0) @binding(4) var<storage, read_write> v: array<f64>;\n\n@compute @workgroup_size({WG})\nfn main(@builtin(global_invocation_id) gid: vec3<u32>) {{\n    let cell = gid.x;\n    if (cell >= NCELL) {{\n        return;\n    }}\n    let i = cell / (NN * NN);\n    let j = (cell / NN) % NN;\n    let k = cell % NN;\n    let o = offp(i + NG, j + NG, k + NG);\n{qload}    let dt = v[0u * NCELL + cell] * q0 + v[1u * NCELL + cell] * q1 + v[2u * NCELL + cell] * q2 + v[3u * NCELL + cell] * q3;\n    v[0u * NCELL + cell] = v[0u * NCELL + cell] - dt * q0;\n    v[1u * NCELL + cell] = v[1u * NCELL + cell] - dt * q1;\n    v[2u * NCELL + cell] = v[2u * NCELL + cell] - dt * q2;\n    v[3u * NCELL + cell] = v[3u * NCELL + cell] - dt * q3;\n}}\n",
        consts = dm.consts(),
        fields = field_bindings(false),
    )
}

/// All kernels for a given N, in inventory order (name, WGSL source).
pub fn all_kernels(n: usize) -> Vec<(&'static str, String)> {
    let dm = Dims::new(n);
    vec![
        ("fill_ghosts", fill_ghosts(dm)),
        ("forward", forward(dm)),
        ("flux", flux(dm)),
        ("gather", gather(dm)),
        ("renormalize", renormalize(dm)),
        ("axpy_v", axpy_v(dm)),
        ("step_q", step_q(dm)),
        ("project_v", project_v(dm)),
    ]
}
