# A_scaling — aeros M0a cost / OpenMP-scaling harness (Option A)

A minimal Fortran benchmark over **SHTns** that reproduces the *transform +
nonlinear + implicit-solve* workload of one semi-implicit spectral
primitive-equation timestep. It exists to answer the three M0a questions from
the design doc (§3.6, §4.3, milestone M0a):

1. real **core-seconds per model-year** at T31L16 and T42L19,
2. the **transform share** of runtime,
3. the **OpenMP scaling curve** to many threads (the number that decides
   T31+§3.7 vs bare T42).

**This is a benchmark, not a dynamical core.** It does no physics and does not
time-integrate a real atmosphere. It runs representative arithmetic on
representative data so the timings are defensible.

## Build

```bash
./build.sh          # sources Intel oneAPI setvars, runs make -> bin/aeros_bench.x
./build.sh clean    # clean rebuild
```

- Compiler: **Intel `ifx`** (LLVM-based). `ifort`/`gfortran` cannot link the
  OpenMP SHTns build. The `Makefile` forces `FC := ifx` because the login
  environment exports `FC=f77`.
- The Intel runtime dir (holding `libiomp5.so`) is **rpath'd into the binary**,
  so it runs on a compute node without sourcing setvars. `build.sh` prints an
  `ldd` check; no library should be *not found*. Verified clean:
  `libimf`, `libiomp5`, `libintlc` all resolve via the rpath.
- External libs (verified OpenMP builds):
  - SHTns: `fesm-utils/SHTns/shtns-omp` (`include/shtns.f03`, `lib/libshtns_omp.a`)
  - FFTW:  `fesm-utils/fftw/fftw-omp`   (`lib/libfftw3_omp.a`, `libfftw3.a`)

## Run / interface

```bash
OMP_NUM_THREADS=<N> ./bin/aeros_bench.x <bench.nml>
```

- `argv(1)` is a Fortran namelist with group `&bench`:
  `lmax, nlat, nphi, nlev, dt, nstep, nwarm, label`.
- **Threads come from `OMP_NUM_THREADS`** (via `omp_get_max_threads()`), *not*
  from the namelist. The value is passed to `shtns_use_threads()` before grid
  init, so it controls both the SHTns internal transform threading and the
  harness's own `!$omp parallel do` regions. The runme sets `OMP_NUM_THREADS`,
  `OMP_PROC_BIND=cores`, `OMP_STACKSIZE=512M`.
- Provided parameter files: `par/bench_T31L16.nml`, `par/bench_T42L19.nml`
  (`nstep=100`, `nwarm=10`, `dt=1800`).

## What each phase represents

Per timestep, with three separate `omp_get_wtime()` accumulators:

- **`t_transform`** — SHTns spherical-harmonic transforms. Inverse: per level a
  vector `SHsphtor_to_spat` (vor,div → u,v) plus scalar `SH_to_spat` for vor,
  div, T, q; once per step for lnps. Forward: the mirror
  (`spat_to_SHsphtor`, `spat_to_SH`). Threading is **internal to SHTns**.
- **`t_gridpoint`** — OpenMP-parallel (`!$omp parallel do` over grid points, each
  point owning all its levels) nonlinear tendency pass: absolute-vorticity flux
  `eta = vor+f`, kinetic energy `0.5(u²+v²)`, heat/moisture advection products
  `u·T, v·T, u·q, v·q`, and the `T·div`, `q·div` compression terms —
  ~30 genuine multiply-adds per point-level (no dummy spin loop).
- **`t_solve`** — OpenMP-parallel over the `nlm` spectral modes: a size-`nlev`
  complex tridiagonal (Thomas) solve standing in for the vertical /
  gravity-wave implicit solve, then an `l`-dependent hyperdiffusion multiply.

Prognostic fields are held **fixed** at their non-trivial initial amplitude
across all timed steps (the benchmark needs representative arithmetic, not a
stable integration). A checksum of the tendencies is accumulated and printed to
defeat dead-code elimination. `nwarm` warmup steps precede the `nstep` timed
steps.

## Caveats (all bias the numbers toward *upper bounds*)

- **Double precision.** The SHTns OpenMP library exposes only
  `C_DOUBLE` / `C_DOUBLE_COMPLEX` interfaces, so the whole harness runs in
  double. The design wants **Float32**, which would be roughly **~2× faster** —
  so these timings are a **conservative upper bound**. Reported as
  `precision = double`.
- **Full Gaussian grid.** The SHTns basic API has no reduced/octahedral grid,
  so a full `SHT_GAUSS` grid is used. A reduced grid removes **~1/3** of the
  gridpoints (design §4), so gridpoint cost here is also an upper bound.
- **Physics excluded** — radiation/condensation/surface are ~1% of budget
  (design §3.6) and are not represented.
- **Coriolis** is a smooth representative field, not geophysically placed; it
  only supplies real O(1e-4) data to the arithmetic and does not affect timing.

## `results.txt` schema

Written to the **current working directory** (and echoed to stdout): a
human-readable block followed by machine-parseable `key = value` lines:

| key | meaning |
|---|---|
| `label` | run label from the namelist |
| `lmax`, `nlev` | truncation and level count |
| `nlat`, `nphi`, `nlm` | grid + spectral sizes, **as reported by SHTns** |
| `nthreads` | threads actually used (`shtns_use_threads` return) |
| `precision` | `single` / `double` (currently `double`, see caveats) |
| `dt` | timestep (s) |
| `nstep` | timed steps |
| `t_step_s` | mean wallclock per step (sum of the three instrumented phases / nstep) |
| `t_transform_frac`, `t_gridpoint_frac`, `t_solve_frac` | phase fractions (sum to 1) |
| `steps_per_year` | `365*86400/dt` |
| `wallclock_s_per_year` | `t_step_s * steps_per_year` |
| `core_s_per_year` | `wallclock_s_per_year * nthreads` |
| `model_years_per_day`, `sypd` | `86400 / wallclock_s_per_year` |

## Variant 2: batched-DGEMM transform (`bin/aeros_bench_batched.x`)

A **second** transform variant. Everything except the spherical-harmonic
transform is **byte-identical** to variant 1 (`src/aeros_bench.f90`,
`variant = spectral_shtns`): the `&bench` namelist, the `OMP_NUM_THREADS`
thread handling, the gridpoint nonlinear pass, and the semi-implicit
spectral-solve stand-in are copied verbatim, so the only thing that differs
between the two binaries is **how the transform is computed and parallelized**.
This is the decisive M0a test: does the design's *custom batched Legendre*
transform (§4) scale where SHTns-per-field did not?

### Algorithm (real SHT on a Gaussian grid, triangular T(lmax), mres=1)

- **Gauss-Legendre grid** computed in-code (`gauleg`): `nlat` nodes
  `mu_j = cos(theta_j)` and weights `w_j` via Newton-Raphson on `P_nlat`;
  `nphi` equally-spaced longitudes. `nlm = (lmax+1)(lmax+2)/2` (= 528 at T31,
  946 at T42 — matches SHTns exactly), `nspat = nlat*nphi`.
- **Fully-normalized associated Legendre functions** `Pbar_l^m(mu_j)`
  (`build_plm`), stable seminormalized recurrence, normalized so
  `sum_j w_j Pbar_l^m^2 = 1`. Precomputed once into m-major blocks
  `plm(nlat, nlm)`; block `m` is `plm(:, moff(m)+1 : moff(m)+ (lmax-m+1))`.
- **Inverse (spectral -> grid)**: per `m`, ONE real `dgemm` batches all `nlev`
  levels (real & imag packed as `2*nlev` columns):
  `Gpack(nlat, 2nl) = Pm(nlat,Km) . Spack(Km,2nl)`, then an inverse real-FFT
  (FFTW c2r) along longitude for every `(lat,level)`.
- **Forward (grid -> spectral)**: forward real-FFT (FFTW r2c) per `(lat,level)`,
  apply Gauss weight `w_j` and the `1/nphi` FFT normalization, then per `m` ONE
  transposed `dgemm`: `Spack(Km,2nl) = Pm^T(Km,nlat) . Ghat(nlat,2nl)`.

The `(u,v) <-> (vor,div)` vector transform is represented, as the design allows,
by the appropriate pair of **scalar** Legendre GEMMs with the derivative/`(l,m)`
factor applied as a cheap `O(nlm*nlev)` spectral multiply (`vfac_mul`). The
per-step transform **count matches variant 1 exactly**: inverse = 6 batched
scalar transforms (`u,v,vor,div,T,q`) + 1 for `lnps`; forward = 4 batched
(`Fu,Fv->dvor,ddiv` + `FT,Fq`) + 1 for `lnps`. Only the parallelization/batching
differs, not the transform workload.

### Parallelization (the thing under measurement)

- **Primary: OpenMP over `m`** for the Legendre GEMMs
  (`!$omp parallel do schedule(dynamic)` over `m=0..lmax`), with **sequential
  MKL BLAS** inside each GEMM (`-qmkl=sequential`) so OUR over-`m` loop owns the
  parallelism, not MKL's internal threading.
- **FFTs**: OpenMP over `(lat,level)` pairs, per-thread scratch, a shared serial
  FFTW plan created with `FFTW_UNALIGNED` (`fftw_execute_dft_{r2c,c2r}` on
  distinct data is thread-safe). Serial `libfftw3` (not the `_omp` build).
- **Expected plateau**: at T31/T42 there are only **32-43** `m`-values with
  strong load imbalance (`m=0` carries the most degrees), so a dynamic schedule
  helps but the curve is expected to **plateau below ~43 threads**. Measuring
  *where* it plateaus is the whole point.

### Validation (`--validate`)

`OMP_NUM_THREADS=N ./bin/aeros_bench_batched.x <bench.nml> --validate` builds a
random band-limited spectral field (real `m=0` coefficients), does
inverse-then-forward, and reports the max round-trip abs error. Measured
(double precision):

| truncation | max round-trip abs err |
|---|---|
| T31L16 | 4.49e-15 |
| T42L19 | 3.13e-15 |

i.e. machine precision — the transform is exact, as required before trusting any
timing. (Exactness relies on `nlat >= lmax+1` so Gauss quadrature is exact, and
on `lmax < nphi/2` so the Nyquist mode is unused.)

### Build / link line

`build.sh` sources the MKL setvars (for `MKLROOT`) then the ifx compilers
setvars, and `make` builds both binaries. Variant 2's link line:

```
ifx -qopenmp -O3 -xHost -fp-model=precise -qmkl=sequential -I$(FFT)/include \
    src/aeros_bench_batched.f90 -qmkl=sequential -L$(FFT)/lib -lfftw3 -lm \
    -Wl,-rpath,$(INTEL_RT) -Wl,-rpath,$(MKL_RT) -o bin/aeros_bench_batched.x
```

`ldd bin/aeros_bench_batched.x` is clean in a stripped env (`env -i`): the MKL
(`libmkl_intel_lp64/sequential/core`) and `libiomp5` runtimes resolve via the
rpath, so the binary runs on a compute node without sourcing setvars. MKL used
here is `intel-oneapi-mkl/2022.1.0-akthm3n` (2024 compilers ship no bundled
MKL).

### `results.txt` gains `variant`

Variant 2 emits all the same keys as variant 1, plus one new key
**`variant = batched`** (variant 1 is understood as `variant = spectral_shtns`;
it was left untouched). `nlat/nphi/nlm` are the in-code grid sizes (identical to
what SHTns reported for variant 1).

### Login-node smoke read (shared node, eyeball only — NOT the real curve)

Unlike SHTns-per-field (which got *slower* with threads), the batched transform
**does drop with threads** on the login node:

| threads | T31 t_step (ms) | T42 t_step (ms) |
|---:|---:|---:|
| 1 | 4.05 | 10.49 |
| 2 | 2.31 | — |
| 4 | 1.30 | — |
| 8 | 1.33 | 2.80 |

Single-core cost is *higher* than SHTns (a naive dgemm-per-m vs a tuned
library), but it **scales** (T31 ~3.1x at 4 threads then plateaus ~4-8 as the
32-`m` limit predicts; T42 ~3.75x at 8). The real scaling curve to high thread
counts must be measured on dedicated pinned compute nodes via the runme.

## Login-node smoke test — what to expect

On the shared login node at T31/T42 the step time **does not drop** with more
threads — it rises. This is expected and is the M0a signal, not a defect:

- the step is **~88–92% SHTns transforms**, and SHTns's internal Legendre
  threading does not pay off at `lmax = 31–42` (too little work per thread —
  exactly the design §4.3 / risk 2 "work-per-thread" limit);
- the login node is shared, and spin-waiting OpenMP threads oversubscribe
  contended cores.

That OpenMP itself scales was confirmed separately: a compute-bound probe on
this toolchain went 16.9 → 8.5 → 4.2 → 2.1 s across 1→2→4→8 threads
(near-linear). The **real** scaling curve for the design decision must be
measured on **dedicated, pinned compute nodes** via the runme — which is what
M0a is for.
