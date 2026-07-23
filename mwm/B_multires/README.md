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
