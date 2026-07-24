# B_multires — Option B, stage b1: multi-resolution error DIAGNOSIS

Feasibility probe for design §3.7 ("multi-resolution error correction"), built on
**SpeedyWeather.jl**. b1's narrow job is to establish that there is a coherent,
sizable, spatially-structured **resolution error in the stationary-wave field**
between a high (T85) and a low (T31) truncation run under *identical boundary
conditions* — i.e. that there is a well-defined ΔF worth correcting.

**b1 DIAGNOSES ONLY. It does not apply any correction — that is b2 (spec'd later).**

## Environment

- `julia` 1.12.1 (juliaup).
- **SpeedyWeather v0.21.1** (resolved by `Pkg.add`; API read directly from the
  installed source — this stage was written against 0.21.1, not older examples).
- NCDatasets 0.14.15, CairoMakie 0.15.13.
- Instantiate on the **login node** (`albedo0`, has internet) — this also caches
  the SpeedyWeatherAssets boundary-condition artifacts (orography, land-sea mask,
  SST climatology) into the depot so the internet-less compute nodes can run:

  ```bash
  julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
  julia --project=. -e 'using SpeedyWeather'   # force artifact download + full precompile
  ```

## What it does

`b1_probe.jl <config.toml>`:

1. Builds a **PrimitiveWetModel** (moist primitive equations, real-Earth
   orography + land-sea mask) at both truncations.
2. Integrates `spinup_days` (discarded) then accumulates a time-mean over
   `mean_days` of geopotential, temperature and wind at the sigma level nearest
   σ=0.5 (~500 hPa), via a per-timestep accumulator callback.
3. Coarse-grains the T85 time-mean to the T31 spectral representation
   (`spectral_truncation` to degree `trunc_lo`), and brings both means onto a
   common **full Gaussian grid at the T31 resolution**.
4. Diagnoses the stationary wave Z\*(λ,φ) = time-mean minus zonal mean, and the
   resolution error ΔZ\* = Z\*_T85(→T31) − Z\*_T31 (and the same for T and u).
5. Writes `output/b1_fields.nc`, `output/b1_summary.toml` (area-weighted RMS
   numbers + locations of ΔZ\* maxima), and PNG maps.

### Model construction (exact call)

```julia
spectral_grid = SpectralGrid(; trunc, nlayers=nlev)      # trunc = 85 or 31; nlev identical
model = PrimitiveWetModel(spectral_grid;
    ocean   = SeasonalOceanClimatology(spectral_grid),   # PRESCRIBED seasonal SST (non-interactive)
    sea_ice = PrescribedSeaIce(spectral_grid),           # non-interactive sea ice
)
```

### How identical boundary conditions are enforced

The **only** intended difference between the two runs is the spectral truncation.
Everything else is a shared source, spectrally truncated at each run's own
resolution (that resolution-dependent smoothing of orography *is* the effect
under study):

- **Orography**: default `EarthOrography` (identical real-Earth source at both).
- **Land-sea mask**: default `EarthLandSeaMask` (identical source).
- **SST**: `SeasonalOceanClimatology` — a *prescribed* monthly SST climatology
  read from assets and interpolated in time. This replaces the SpeedyWeather
  default `SlabOcean`, which is an interactive mixed-layer that would let the two
  runs drift apart through ocean state. This is the crux of "identical BCs".
- **Sea ice**: `PrescribedSeaIce` — replaces the default interactive
  `ThermodynamicSeaIce`, so the surface is fully non-interactive over ocean.
- **Vertical**: same `nlev` and same default sigma coordinate at both.
- **Physics**: identical default parameterizations at both truncations (no
  per-resolution parameter tuning is applied here).

### Confounds worth flagging (design §3.7 risk 3)

- **Timestep** differs per resolution (SpeedyWeather picks Δt from a T31
  reference scaled by CFL: ~40 min at T31, shorter at T85). This is intended and
  permitted by the spec.
- **Hyperdiffusion** (`HyperDiffusion`) is resolution-dependent *by
  construction* — it must scale with truncation for numerical stability. It is
  numerics, not a physics tuning, but it is a genuine non-truncation difference
  between the runs and so is part of the diagnosed ΔF.
- **Interactive land surface**: the land model (`LandModel`, prognostic
  soil/bucket temperature) is left interactive and identical at both resolutions.
  It is a deterministic response to the atmosphere, not a divergent slow
  reservoir like the ocean, so it is not prescribed here. Its difference between
  T85 and T31 is legitimately part of the resolution error.

Per §3.7 risk 3, ΔF conflates genuine dynamical resolution error with any
resolution-dependent numerics/physics; that is acceptable for the b1/b2 goal
(T31+ΔF ≈ T85) but means ΔF cannot be read as pure dynamics.

## Diagnostic conventions

- **~500 hPa**: the full sigma level nearest σ=0.5 is used. At `nlev=8` this is
  σ=0.4375 (~440 hPa). Both runs use the same level index.
- **Common grid**: full Gaussian grid at the T31 `nlat_half` — regular in
  longitude, Gaussian in latitude — so zonal means and RMS are straightforward.
- **RMS**: cos(latitude) area-weighted.

## Running

Production (compute node, via runme — do **not** run T85 on the login node):

```bash
# runme info.json exe alias "speedy" -> B_multires/run_speedy.sh, par path as $1
# par/b1.toml: trunc_hi=85, trunc_lo=31, nlev=8, spinup_days=200, mean_days=360
```

`run_speedy.sh` cd's here, sets `JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}`
(SpeedyWeather runs on the KernelAbstractions `CPU()` backend and threads its
kernels/transforms over Julia threads — no special model/architecture setting is
needed for multi-threaded CPU execution), and uses an absolute `--project`.

Smoke test (login node, tiny — validates the whole pipeline):

```bash
julia --project=. b1_probe.jl par/b1_smoke.toml   # T42/T21, 5+10 days -> output_smoke/
```

Plotting is decoupled and can be re-run after the heavy job:

```bash
julia --project=. plot_b1.jl output/b1_fields.nc output
```

## Files

- `b1_probe.jl`   — diagnosis driver.
- `plot_b1.jl`    — CairoMakie maps (Z\*_hi, Z\*_lo, ΔZ\*, ΔT\*, Δu\*), standalone.
- `run_speedy.sh` — runme wrapper (absolute `--project`, threads from SLURM).
- `par/b1.toml`   — production config (T85/T31, 200+360 days).
- `par/b1_smoke.toml` — smoke config (T42/T21, 5+10 days).
- `output/`       — production outputs (gitignored). `output_smoke/` for the smoke test.

## Production config values (par/b1.toml) and rationale

| key | value | why |
|---|---|---|
| `trunc_hi` | 85 | L&L reference resolution that builds LGM ice sheets |
| `trunc_lo` | 31 | cheap transient truncation that does not |
| `nlev` | 8 | identical vertical config; SpeedyWeather default |
| `spinup_days` | 200 | reach statistical equilibrium of the stationary-wave field |
| `mean_days` | 360 | full-year time mean, removes seasonal aliasing |
| `target_sigma` | 0.5 | ~500 hPa (nearest full sigma level) |

## Smoke-test result (T42/T21, NOT scientific)

Pipeline validated end-to-end: both runs complete, diagnostics compute, NetCDF +
summary + 3 PNGs written. RMS Z\* ≈ 427 m (hi), ΔZ\* ≈ 118 m, ratio ≈ 0.28.
Even un-equilibrated, the ΔZ\* maxima already concentrate over the expected
orography — Tibet (~84°E, 41°N), the Rockies (~247°E, 41°N) and Greenland
(~321°E, 69°N) — the design-relevant signature.

---

# Stage b2 — multi-resolution error CORRECTION test

b1 established that a well-defined resolution error ΔF *exists* (ΔZ\* ≈ 27% of the
T85 stationary-wave signal, coherent over orography). **b2 tests the actual §3.7
assumption: does injecting ΔF back into T31 as a constant additive forcing make
T31 reproduce T85's stationary waves?**

b2 is built as NEW files; b1 (`b1_probe.jl`, `run_speedy.sh`, `plot_b1.jl`,
`par/b1*.toml`) is left untouched.

## The §3.7 formulation implemented (approach A — faithful first cut)

§3.7 is explicit: *"diagnose the model-error operator, do not nudge toward a
state."* So the correction is a **constant additive forcing on the prognostic
tendencies**, not a relaxation toward a target state.

1. **Run T85** (identical prescribed BCs as b1: `SeasonalOceanClimatology` +
   `PrescribedSeaIce`, default `EarthOrography`/`EarthLandSeaMask`, same `nlev`),
   spin up, then accumulate the **time-mean of the full prognostic spectral
   state** (vorticity, divergence, temperature, humidity, log-surface-pressure —
   all levels, spectral coefficients). Coarse-grain (spectral-truncate) that mean
   to T31 → **X85** (the target, expressed as a T31-resolution spectral state).
2. **Diagnose ΔF = −RHS_T31(X85):** set a T31 model to X85 and evaluate the
   **dynamical-core tendency once** (`dynamics_tendencies!`). X85 is not T31's own
   equilibrium, so this tendency is nonzero; the constant forcing that makes X85 a
   fixed point of the T31 core is its negative.
3. **Inject ΔF** into a T31 run via a custom `AbstractForcing`
   (`ConstantTendencyForcing`, see `b_common.jl`), which adds the precomputed
   constant ΔF onto the tendencies every timestep. Run **T31+ΔF** with the same
   spinup/mean protocol as bare T31.
4. **Test reproduction.** Compute the stationary waves Z\*(~500 hPa), T\*, U\* for
   **X85 target**, **bare T31**, **T31+ΔF** (all three from their mean spectral
   states, diagnosed identically), and report the **fractional RMS error
   reduction** `(RMS(bare−target) − RMS(corr−target)) / RMS(bare−target)`.
5. **Conservation check** (§3.7 risk 2): area-weighted global mean of each ΔF
   field (should be ~zero-mean).

### SpeedyWeather v0.21.1 API hooks used

- **Reading the tendency at a prescribed state**: `set!` the prognostic spectral
  state on both leapfrog steps, `reset_tendencies!`, `transform!(vars, lf, model)`
  (spectral→grid diagnostics, incl. geopotential inside `dynamics_tendencies!`),
  then `dynamics_tendencies!(vars, lf, model)` and read
  `vars.tendencies.{vorticity,divergence,temperature,humidity,pressure}`.
- **Injecting constant forcing**: a custom `<: SpeedyWeather.AbstractForcing` with
  a `forcing!(vars, forcing, lf, model)` method, passed as `forcing=` to
  `PrimitiveWetModel`. **Injection point matters** (verified against the source):
  the model builds the spectral vor/div/pres tendencies by *accumulating*
  (`add=true`), so those ΔF are added in **spectral** space; but the temperature
  and humidity spectral tendencies are *overwritten* by a grid→spectral transform
  of the grid tendency (`temperature_tendency!`, `humidity_tendency!`), so those
  ΔF must be injected into the **grid-space** tendency
  (`vars.tendencies.grid.{temperature,humidity}`), the same path `HeldSuarez`
  uses. `ConstantTendencyForcing` stores `dtemp_grid/dhumid_grid =
  transform(−RHS)` for exactly this reason.
- **Scaling**: SpeedyWeather radius-scales vorticity/divergence during a run
  (`scale = radius`); the mean state is accumulated, X85 set, and ΔF diagnosed all
  in that same scaled representation (radius is identical at both truncations), so
  the scaling cancels and ΔF is injected in the units the leapfrog integrates.

### Injection-applied verification (guards against a silent null result)

Before trusting the climate result, b2 confirms the forcing is actually applied:
it evaluates the **forced** dynamical-core tendency at X85 and checks it is ≈ 0
(X85 becomes a fixed point of the forced core). Smoke-test residual fractions
`max|forced| / max|baseline|`: **vor 0.0, div 5e-4, temp ~2e-6, humid ~2e-7,
pres 0.0** — i.e. ΔF cancels −RHS to Float32 round-off (the div residual is
round-off through the `∇²` of the large geopotential in `bernoulli_potential!`).
A forcing that silently was not applied would leave `residual_frac ≈ 1`.

## Files (b2)

- `b_common.jl`     — shared machinery: `SpectralMeanAccumulator` (full-state
  time-mean callback), `ConstantTendencyForcing` (+ `forcing!`), `set_state!` /
  `eval_dynamics!`, coarse-graining, and diagnostic helpers. `include`d by b2;
  does **not** touch b1.
- `b2_probe.jl`     — b2 driver (reads config `ARGS[1]`).
- `plot_b2.jl`      — CairoMakie maps (Z\* target/bare/corr triptych; Z\*, T\*, U\*
  residual maps), standalone → `output_b2/`.
- `run_speedy_b2.sh`— runme wrapper (resolves config to ABSOLUTE **before** `cd`,
  absolute `--project`, `JULIA_NUM_THREADS` from SLURM). Analogous to
  `run_speedy.sh`. (To wire into runme, add a `speedy_b2` exe alias pointing here;
  `.runme/` was intentionally left untouched.)
- `par/b2.toml`     — production (T85/T31, nlev=8, 200+360 days).
- `par/b2_smoke.toml` — smoke (T42/T21, nlev=6, 2+3 days) → `output_b2_smoke/`.
- `output_b2/`      — production outputs (`b2_fields.nc`, `b2_summary.toml`, PNGs).

## Running

```bash
# smoke (login node, tiny — validates the whole pipeline incl. injection check):
julia --project=. b2_probe.jl par/b2_smoke.toml

# production (compute node — do NOT run T85 on the login node): three integrations
# (T85 reference, bare T31, T31+ΔF) + ΔF diagnosis, via run_speedy_b2.sh.
```

## Caveats (documented honestly in `b2_summary.toml`)

- **Single-tendency-at-mean-state neglects transient rectification.** The
  instantaneous RHS at the time-mean state ≠ the time-mean RHS, because of the
  nonlinear eddy terms. If T31+ΔF does not fully reproduce the target, this is the
  leading candidate — and a partial result is itself informative: it says the
  constant-forcing-from-a-single-eval approximation is insufficient and the
  iterated / mean-tendency version is the next step.
- **ΔF is the dynamical-core tendency error** (`dynamics_tendencies!`). Physics
  parameterizations respond freely each step and are not part of ΔF, so any
  residual physics tendency at X85 also limits reproduction. This is aligned with
  §3.7's "correct selected terms, not everything," but is a deliberate scope
  choice to record.
- **Full-field correction** (all of vor/div/T/humid/lnps), not §3.7's later
  "selected terms" refinement.
- Same b1 confounds: resolution-dependent hyperdiffusion + timestep, interactive
  land; ΔF conflates genuine dynamical resolution error with resolution numerics.
- **SpeedyWeather ≠ aeros** → suggestive, not conclusive.

## Smoke-test result (T42/T21, NOT scientific)

Pipeline validated end-to-end: T85-analog run, ΔF diagnosis, **injection verified
applied** (see above), bare-T31 and T31+ΔF runs, three-way diagnostics, NetCDF +
`b2_summary.toml` + PNGs all written. Conservation of ΔF is ~zero-mean
(`|mean|/RMS`: temp 0.006, humid 0.013, vor/div/pres 0). The error-**reduction**
numbers at smoke scale are negative (correction "worsens" the field) purely
because a 2-day spinup / 3-day mean at T42/T21 is nowhere near equilibrium — not a
scientific result. The production config (200+360 days, T85/T31) is the real test.
