# PORTABILITY — building fs-la and fs-fft from out-of-workspace crates (Phase F3)

Date: 2026-07-12. Scope: why standalone crates in this sandbox could not
path-depend on `crates/fs-la` / `crates/fs-fft`, the minimal fix applied to
the fork's workspace, and the resulting swap status of the two gum-core
crates that had worked around it.

## Diagnosis

The broken edge is a single one, and it is direct in both crates:

```
fs-la  --[dependencies]-->  fs-exec  --[dependencies]-->  asupersync = { path = "../../../asupersync" }
fs-fft --[dependencies]-->  fs-exec  --(same edge)
```

`asupersync` is a SIBLING REPOSITORY expected next to the fork checkout
(`/home/user/asupersync` here). It is absent, so any Cargo resolution that
includes fs-exec fails at manifest load:

```
error: failed to get `asupersync` as a dependency of package `fs-exec v0.0.1`
  ... failed to read `/home/user/asupersync/Cargo.toml`
```

Reproduced before the fix with a scratch crate (`cargo new` +
`fs-la = { path = ... }`): the failure is exactly the above, via both fs-la
and fs-fft. Nothing else in either crate's normal/build dependency closure
(fs-alloc, fs-math, fs-simd, fs-rand, fs-soa, fs-obs, fs-substrate, fs-qty,
fs-soa-derive, fs-blake3) leaves the repository.

Two further couplings sit behind the same edge in fs-la only:

1. `crates/fs-la/build.rs` shells out to `git -C ../asupersync` and reads
   the sibling repo's sources to compute the GEMM build fingerprint
   (`FS_LA_GEMM_BUILD_FINGERPRINT`), and fails closed without
   `FRANKENSIM_DEPGRAPH_RECEIPT`/`FRANKENSIM_DEPGRAPH_SALT` (the salt is
   supplied by the workspace's `.cargo/config.toml`, which out-of-tree
   consumers do not inherit).
2. `fs_exec` types appear in fs-la's public API (the cancellation-aware
   pooled GEMM family) and in fs-fft's pooled `FftNd` API.

Note the fork's WORKSPACE ROOT does not resolve in this checkout at all
(e.g. `cargo check -p fs-la` fails while loading member manifests, via
fs-dd -> fs-geom -> fs-exec -> asupersync). That was true before Phase F3
and is unchanged by it.

## Fix (minimal, feature-based)

Both crates gained an `exec` cargo feature that is ON BY DEFAULT and
carries the entire fs-exec coupling:

- `crates/fs-la/Cargo.toml`, `crates/fs-fft/Cargo.toml`:
  `fs-exec = { path = "../fs-exec", optional = true }`;
  `[features] default = ["exec"]; exec = ["dep:fs-exec"]`.
- `crates/fs-la/src/gemm.rs`: `#[cfg(feature = "exec")]` on the pooled/
  cancellation-aware GEMM surface (`gemm_f64_parallel_with_cancel`,
  `gemm_f64_parallel_with_pool{,_declared,_budgeted}`, `gemm_panel_run_id`,
  `GemmRunReport`/`GemmCancelled`/`GemmRunError`, the band-kernel and
  memory-preflight internals) and on the build-identity constants
  (`GEMM_BUILD_FINGERPRINT`, `GEMM_GRAPH_EVIDENCE*`, `GEMM_DEPGRAPH_*`,
  `gemm_build_identity`, `gemm_graph_evidence`). The serial and
  `std::thread`-scoped kernels (`gemm_f64`, `gemm_f32`, `gemm_mixed`,
  `gemm_f64_op`, `gemm_f64_parallel{,_with}`), `Trans`, the memory-envelope
  types, and ALL other modules (`eigen`, `eigen_complex`, `factor`,
  `batched*`, `mixed`, `rand_nla`) remain unconditional.
- `crates/fs-la/src/lib.rs`: the re-export list is split into an
  unconditional block and an `exec`-gated block.
- `crates/fs-la/build.rs`: early-returns when `CARGO_FEATURE_EXEC` is
  unset — no asupersync checkout, no depgraph evidence env, no fingerprint
  (the constants that consume it are compiled out in that configuration).
  With `exec` on, the script is behaviorally identical to before.
- `crates/fs-fft/src/lib.rs`: `#[cfg(feature = "exec")]` on
  `PencilBlockKernel`, `PencilColumnKernel`, and the pooled `impl FftNd`
  block (`forward_pooled`/`inverse_pooled`/`run_pooled`). The serial
  `Fft`, `RealFft`, DCT, and `FftNd` surface is unconditional.
- Executor-coupled integration tests in both crates are declared with
  `required-features = ["exec"]` (same pattern as fs-fft's existing
  `frontier-sixstep` target).

Default-ON (rather than the default-OFF first suggested) is what keeps the
change minimal: ~20 workspace crates consume the pooled GEMM API and would
each have needed a feature line under default-OFF, whereas default-ON
changes NOTHING for any in-workspace consumer — with `exec` enabled the
compiled token stream is the pre-F3 one. Out-of-tree consumers opt out:

```toml
fs-la  = { path = ".../crates/fs-la",  default-features = false }
fs-fft = { path = ".../crates/fs-fft", default-features = false }
```

Caveats of the no-`exec` configuration: `cargo test` on fs-la/fs-fft
themselves still wants default features (unit-test modules exercise the
pooled paths), and the GEMM build fingerprint does not exist (it is a
tune-key identity for the executor path, meaningless without it).

## Verification

- Scratch crate outside the workspace (`f3-smoke`), path-depping both
  crates with `default-features = false`: resolves, builds with zero
  warnings (clippy clean), and a smoke test passes — `jacobi_eigh` on a
  3x3 symmetric matrix recovers eigenvalues {1, 3, 4} to 1e-12, and an
  `FftNd` 4x8 forward/inverse round-trip is exact to 4.4e-16.
- Workspace root: `cargo check -p fs-la -p fs-fft` fails identically
  before and after F3 (missing sibling repo, pre-existing) — no
  regression is possible to observe in this checkout; in a checkout WITH
  asupersync, `exec` being default-ON makes every in-workspace build
  identical to pre-F3 by construction.

## Swap status

### fs-gum-cosserat — DE-VENDORED (done)

`src/fsla_vendored.rs` (verbatim copies of `eigen::jacobi_eigh` and
`eigen_complex::eig` + helpers) is DELETED. The crate now depends on fs-la
directly (`default-features = false`); `src/eigh.rs` imports
`fs_la::eigen::jacobi_eigh` and `fs_la::eigen_complex::eig`. Before the
swap the vendored bodies were re-diffed against fs-la: byte-identical.
Gate outcome after the swap: **26/26 PASS**, 26/26 claims certified,
fs-checker 3/3, replay-deterministic, golden Merkle root UNCHANGED at
`9ddc0fe40548bf9c97cd505ef30339bbe8b4b8ff537ef6f3ef64377cb5f57b15`
(gate output bit-identical to the pre-swap run modulo the wall-clock
line); crate unit tests 13/13.

### fs-gum-topo — NOT swapped (deliberate)

fs-fft is now importable from the crate (same `default-features = false`
recipe), but the internal radix-2 FFT (`src/fft.rs`) stays: it is
golden-gated, and fs-fft's Stockham kernels are not guaranteed to
reproduce its butterfly order bit-for-bit, so a swap risks a golden
re-freeze for no functional gain. `src/fft.rs`'s module docs and
RESULTS.md record the future-swap recipe; keep the structured `NotPow2`
rejection at the API boundary when swapping (fs-fft's `Fft::new` panics
on non-power-of-two sizes rather than returning a typed error).
