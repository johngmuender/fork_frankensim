# `-C target-cpu=native` re-measured on the post-restart host (2026-07-16)

The G-A survey's finding #3 ("native is a measured 10–25% regression on
the real kernels") was measured on the original host (Cascade Lake-class
Xeon @ 2.80 GHz). The 2026-07-16 container restart moved the session to a
**Sapphire-Rapids-class Xeon @ 2.10 GHz** (avx512_fp16/bf16, AMX,
AVX-VNNI). Re-measured per the user's request.

## Protocol
Standalone A/B crate (scratchpad `natbench`) path-depping fs-gum-field +
fs-gum-statics, timing the real kernels at N=96 (ε=0.05 hedgehog,
best-of-5): forward sector sweep (central4), gradient eval (central4 and
corner — the ANF hot path), plus a synthetic 8-chain mul_add loop.
Baseline vs `RUSTFLAGS="-C target-cpu=native"` (verified: 0 vs 13 zmm
sites in the binaries). Three interleaved A/B pairs, run under the
constant ~3.4-core load of the concurrent N=192 crossing run — the
interleaving makes the RATIOS contention-fair; absolute ns/pt are
mildly inflated.

## Results (ns/pt, best over 3 pairs; ratio = native/baseline)

| kernel | baseline | native | ratio |
|---|---|---|---|
| forward central4 | 52.8 | 44.3 | **0.84 — native 16% FASTER** |
| gradient central4 | 205.6 | 208.6 | 1.01 (parity) |
| gradient corner (ANF hot path) | 391.8 | 392.0 | 1.00 (parity) |
| synthetic FMA chain (GFLOP/s) | 0.7 | 11.3 | 16× (vectorized) |

## Conclusions
1. **The old-host regression does NOT reproduce here**: native is
   neutral on the gradient sweeps and +16% on the forward sweep. The
   regression was host/µarch-specific (Cascade-Lake AVX-512 frequency/
   codegen interplay), not a property of the kernels. GAP_ANALYSIS_v4's
   finding #3 is annotated accordingly.
2. **The strategic conclusion is unchanged**: flags are not where the
   headroom is — the ANF-dominating gradient sweeps sit at parity, so
   explicit k-contiguous SIMD remains the documented unclaimed step.
   But the *policy* softens: on this host, building gum-core binaries
   with `target-cpu=native` is safe-to-mildly-beneficial. Determinism
   note: a native build changes instruction selection and MAY change
   fp contraction — any switch of shipped gate binaries to native
   requires a re-golden pass (same class as the G1 tiled-order bump);
   until someone does that, shipped goldens stay on the baseline ISA.
3. **This host is substantially faster per-core anyway**: baseline
   corner-gradient 392 ns/pt vs 553–622 on the old host — the standing
   perf baselines in scale_survey_v4.json are old-host numbers and
   should be re-baselined if precise budgeting is needed.
