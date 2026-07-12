//! Tier-1 CI gate (no GPU required, milliseconds): every WGSL f64 kernel
//! parses AND validates under naga with `Capabilities::FLOAT64`, and the
//! f64 extension is genuinely capability-gated (rejected with
//! `Capabilities::empty()` — the survey's negative probe).

use naga::valid::{Capabilities, ValidationFlags, Validator};

#[test]
fn wgsl_kernels_validate_with_float64() {
    for n in [8usize, 32, 48] {
        for (name, src) in fs_gum_gpu::shaders::all_kernels(n) {
            let module = naga::front::wgsl::parse_str(&src)
                .unwrap_or_else(|e| panic!("{name} (n={n}): parse error: {e:?}"));
            let mut v = Validator::new(ValidationFlags::all(), Capabilities::FLOAT64);
            v.validate(&module)
                .unwrap_or_else(|e| panic!("{name} (n={n}): validation error: {e:?}"));
        }
    }
}

#[test]
fn f64_rejected_without_capability() {
    let (_, src) = &fs_gum_gpu::shaders::all_kernels(8)[1]; // forward kernel
    let module = naga::front::wgsl::parse_str(src).expect("parse");
    let mut v = Validator::new(ValidationFlags::all(), Capabilities::empty());
    assert!(
        v.validate(&module).is_err(),
        "f64 must be rejected without Capabilities::FLOAT64"
    );
}
