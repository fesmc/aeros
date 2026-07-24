# M0a scaling results (Option A) — first measurement

Run: AWI Albedo, `smp` partition, exclusive 128-core nodes (2×64, 1 thread/core),
`ifx` + SHTns-omp 3.7.5 + FFTW-omp. **Double precision, dynamics-only** (transforms
+ nonlinear gridpoint pass + semi-implicit spectral-solve stand-in; **no physics**).
100 timed steps after 10 warmup, Δt=1800 s. Full Gaussian grid (no reduced/octahedral).
Data: `../A_results.csv`.

## Scaling (speedup vs 1 thread)

**T31L16** (1-thread step = 2.41 ms):

| thr | t_step(ms) | speedup | par_eff | yr/day | core-s/yr | xform% |
|----:|----:|----:|----:|----:|----:|----:|
| 1 | 2.41 | 1.00 | 1.00 | 2048 | 42 | 89% |
| 2 | 1.77 | 1.36 | 0.68 | 2794 | 62 | 90% |
| **4** | **1.22** | **1.97** | **0.49** | **4032** | 86 | 89% |
| 8 | 2.00 | 1.20 | 0.15 | 2463 | 281 | 92% |
| 16 | 2.88 | 0.84 | 0.05 | 1712 | 807 | 95% |
| 32 | 6.18 | 0.39 | 0.01 | 797 | 3467 | 97% |
| 64 | 23.1 | 0.10 | 0.00 | 213 | 25946 | 99% |
| 128 | 21.5 | 0.11 | 0.00 | 230 | 48130 | 95% |

**T42L19** (1-thread step = 5.75 ms):

| thr | t_step(ms) | speedup | par_eff | yr/day | core-s/yr | xform% |
|----:|----:|----:|----:|----:|----:|----:|
| 1 | 5.75 | 1.00 | 1.00 | 858 | 101 | 88% |
| 2 | 4.41 | 1.30 | 0.65 | 1117 | 155 | 90% |
| **4** | **2.46** | **2.34** | **0.58** | **2005** | 172 | 88% |
| 8 | 3.38 | 1.70 | 0.21 | 1459 | 474 | 90% |
| 16 | 4.73 | 1.22 | 0.08 | 1043 | 1325 | 93% |
| 32 | 7.20 | 0.80 | 0.02 | 685 | 4034 | 96% |
| 64 | 29.4 | 0.20 | 0.00 | 168 | 32943 | 99% |
| 128 | — | — | — | — | — | (node-flaky; collapse inferable) |

## Findings

1. **The core saturates at ~4 threads.** Peak throughput at 4 threads for both
   truncations (T31: 4032 yr/day; T42: 2005 yr/day). Beyond 8 threads throughput
   *collapses* — at 64–128 threads it is 10–20× slower than 1 thread. Effective
   parallelism N_eff ≈ 2–2.3. This **contradicts §4.3's "target good scaling to
   ~32 threads"** and **confirms risk 2** (work-per-thread limit) for the
   SHTns-per-field transform.

2. **Transforms dominate and are what fail to scale.** 88–99% of step time is
   SHTns Legendre transforms; the fraction *rises* with thread count because the
   gridpoint/solve OpenMP regions do scale a little while the transforms do not.
   At lmax=31–42 the per-field Legendre step has too little work for SHTns's
   internal threading to help.

3. **Single-core cost is ~2× cheaper than the §3.6 estimate.**
   - T31L16: **42** core-s/yr measured vs **57** estimated.
   - T42L19: **101** core-s/yr measured vs **227** estimated (the §11-flagged
     unverified number was ~2.3× too pessimistic).
   These are dynamics-only, double; Float32 ≈ 2× faster again.

## The escape routes this does NOT rule out

- **Dynamics-only.** Physics (radiation/condensation/surface) is embarrassingly
  parallel over columns (§4.3) and scales to any thread count. It dilutes the
  transform-bound fraction, so full-model N_eff > this. This measured the
  worst-case component in isolation.
- **SHTns-per-field ≠ the design's transform.** §4 specifies a *custom
  batched-DGEMM* Legendre transform (all levels of a wavenumber in one call,
  threaded BLAS). SHTns transforms one 2-D field at a time and does not thread
  the small-lmax step. Whether the batched approach scales where SHTns does not
  is now the key open M0a question.
- **Float32** (~2× and better SIMD/cache, where §3.4 says the speed lives).

## Variant 2: batched-DGEMM transform (custom SHT, §4)

Second harness `aeros_bench_batched.x`: identical gridpoint/solve workload, but
the transform is a custom Gauss-grid SHT with one BLAS GEMM per zonal wavenumber
`m` batching all levels, parallelized **OpenMP-over-m** (sequential MKL inside),
FFTW in longitude. Round-trip validated to 3e-15. Same sweep, compute nodes.

Speedup vs each variant's own 1-thread time; yr/day (double, dynamics-only):

**T42L19**

| thr | SHTns spd | SHTns yr/day | Batched spd | Batched yr/day |
|----:|----:|----:|----:|----:|
| 1 | 1.0 | 863 | 1.0 | 477 |
| 4 | 2.4 | 2053 | 3.5 | 1678 |
| 8 | 1.7 | 1461 | 4.0 | 1897 |
| 16 | 1.1 | 938 | **4.2** | **1997** |
| 32 | 0.8 | 659 | 4.1 | 1931 |
| 64 | 0.2 | 168 | 4.1 | 1940 |
| 128 | 0.05 | 42 | 2.9 | 1400 |

**T31L16**: SHTns peaks 3986 yr/day @4 then collapses; batched holds ~3900–4070
yr/day flat from 4→64 threads (plateau ~3.2×), best 4068 @4.

### What the comparison shows

1. **Batched scales where SHTns collapses.** SHTns peaks at 4 threads then falls
   off a cliff (T42 @128 = 42 yr/day). Batched reaches a **broad plateau**
   (T42 ~4.2× over 8–64 threads; T31 ~3.2× over 4–64) with **no collapse**.
2. **But it does not raise the throughput ceiling.** Batched's *best* ≈ SHTns's
   best (T42: 1997 vs 2053; T31: 4068 vs 3986). Better scaling only claws back a
   worse single-core constant (naive dgemm-per-m is ~1.7–1.8× slower serial than
   tuned SHTns).
3. **And it is less core-efficient at that throughput.** ~2000 yr/day T42 costs
   SHTns 4 cores (168 core-s/yr) but batched 16 cores (693 core-s/yr).
4. **Neither reaches N_eff ≥ ~24.** Effective parallelism: SHTns ≈ 2, batched
   ≈ 4.2 (T42). Both far below the §3.6 threshold for bare T42 at 5000 coupled.
   The batched plateau (~4×) is set by the OpenMP-over-m load imbalance
   (only 32–43 `m`-values); a finer (m,level)-pair split *might* raise it, at the
   cost of GEMM level-batching — untested.

## Bottom line for §10.3

**The "scale the core to ~32 threads → bare T42 fits the 5000 yr/day coupled
target" route is closed by BOTH transforms.** SHTns per-field stalls at N_eff≈2;
the batched-DGEMM transform scales to a broad plateau but only ~4×, does not
raise the single-run ceiling (~2000 yr/day T42 double, ~4000 est. Float32), and
is less core-efficient there. A single coupled T42 run therefore lands at the
design's own fallback of **~2500–3000 yr/day, not 5000** (§9 risk 1 confirmed).

Consequences:
- **The 5000 target requires T31** (dynamics ~4000 yr/day double @4 threads →
  ~8000 Float32; coupled ~4000–5000 plausible), **+ §3.7 for the science**
  (whether T31+Δ reproduces T85 is exactly what Option B tests).
- **Ensemble-parallelism is the right throughput model** (§4.3 vindicated): the
  core wants only ~4 threads, so spend the other ~124 cores on ensemble members.
  For ensembles, **SHTns-per-field @~4 threads is the more core-efficient
  choice** — the batched transform's non-collapse matters only when a *single*
  member must be pushed past ~4 cores, and even then it buys no higher ceiling.
- Remaining single-run headroom lives in **Float32** (~2×) and possibly a
  (m,level)-pair batched transform — not in throwing more cores at SHTns.
