# aeros M1 — dry dynamical core: scope

Status: **scope draft** (2026-07-23), to be picked up later. Companion to
[`design.md`](design.md) §7 (milestones) and the `mwm/` feasibility probes. Start
here; read `design.md` §3–§4 for the reasoning behind the choices below.

## Goal

A dry spectral primitive-equation dynamical core, validated against **Held–Suarez**,
with **truncation and vertical levels as namelist parameters from day one** and a
**first-class, pluggable correction framework** baked in — so §3.7/§3.8 (the
multi-resolution error correction, the pivot of the whole plan) can be tested at
M2b without re-architecting. No moist physics, no radiation, no coupling yet.

## What the MWM (mwm/) already settled — and how it shapes M1

The feasibility probes are done; their verdicts are load-bearing here.

1. **Parallelism: the core saturates at ~4 threads** (M0a, `mwm/A_scaling/RESULTS.md`).
   Neither SHTns-per-field nor a custom batched-DGEMM transform reaches N_eff ≥ ~24.
   → **Design for ensemble parallelism, not core scaling.** The driver must make it
   trivial to run many independent members (one per ~4 cores) on a node. Truncation
   affordable when coupled is **T31**, not bare T42 — but keep truncation a namelist
   knob (see below).
2. **Transform: use SHTns per-field** (`aeros_spectral.f90` over `shtns.f03`). It is
   faster single-core and more core-efficient per ensemble member than the batched
   transform, which scales to a broad plateau but **does not raise the throughput
   ceiling** and is less core-efficient. Keep the batched-DGEMM path documented as a
   fallback only if single-run wallclock ever becomes the binding constraint. The
   MWM already wrote and validated a working SHTns Fortran wrapper + build
   (`ifx` + `libshtns_omp` + `libfftw3_omp`, rpath'd) — reuse it as the M1 starting
   point.
3. **Single-core cost is ~2× better than the §3.6 estimate** (T42L19 ≈ 100 core-s/yr
   measured vs 227 estimated; T31L16 ≈ 42 vs 57). Double precision; Float32 ≈ 2×
   faster again → **Float32 throughout the core** (§4).
4. **The correction is the make-or-break, and it is delicate** (b1/b2, `mwm/B_multires/`).
   b1 confirmed a real, sizable (27% of signal), orographically-structured T85−T31
   resolution error exists to correct. b2 showed a *naive constant* correction is
   insufficient and that even the better *mean-tendency* correction is numerically
   fragile (nudging toward an unbalanced high-res mean excites gravity waves). The
   lesson for M1: **the correction framework must support, from the start**, (a)
   per-field selection (§3.7 "selected terms"), (b) per-scale application (large-scale
   / spectral-filtered), (c) switchability and a `correction=0` twin (§3.8), and (d)
   correct per-field tendency scaling when injecting a correction (the MWM worked out
   the SpeedyWeather convention; aeros needs its own, documented, analogue).

## Architecture

- **Core**: semi-implicit leapfrog Eulerian spectral primitive equations,
  hydrostatic, triangular truncation, σ or hybrid σ–p vertical (§4). Float32.
- **Transform**: `aeros_spectral.f90` wrapping `shtns.f03` (per-field), OpenMP.
  Batched-DGEMM transform kept as a documented alternative, not the default.
- **Grid**: Gaussian. Decide reduced/octahedral vs full — reduced saves ~⅓ of the
  gridpoint work (§4) but is **not in the SHTns basic API**, so it needs either a
  custom reduced-grid path or acceptance of the full-grid cost. Flag as an M1
  decision; full grid is the safe first cut.
- **Resolution as namelist parameters**: `trunc` (T31/T42/T63) and `nlev` (L16–20,
  with 3–4 levels below 850 hPa — §4.1) must be namelist changes from day one
  (§9 risk 1, §10.3). Benchmark T42 vs T63 at M1 to keep the scaling curve current.
- **Correction framework (the M1-specific requirement)**: a pluggable additive
  tendency-correction layer applied inside the tendency computation, structured as a
  sum of independent terms
  `ΔF_total = Σ_k ΔF_k`, each term carrying: which prognostic fields it acts on,
  a spectral-scale selector (e.g. apply only for total wavenumber ≤ L_c), a global
  on/off switch, and correct per-field tendency scaling. At M1 the only term is a
  zero/no-op (validated to change nothing), but the plumbing — namelist config,
  conservation diagnostics, the `ΔF=0` twin — is built and tested now so M2b plugs
  §3.7's `ΔF_res(t)` and §3.8's `ΔF_bias` into an existing, exercised interface.
- **Driver / ensembles**: a thin driver that runs one member from a namelist, plus a
  convention (mirroring `mwm/.runme/`) for launching N independent members per node.
  Each member ~4 OpenMP threads; the node's parallelism comes from member count.

## Validation

- **Held–Suarez** (Held & Suarez 1994): the standard dry-core benchmark. Reproduce
  the zonal-mean zonal wind, temperature, and eddy statistics. This is the M1 pass/fail.
- **Benchmark** T42 vs T63 core-s/yr and the OpenMP curve (extends M0a to the real core).
- **Conservation**: mass and energy to expected tolerance; verify the (no-op)
  correction term changes nothing to machine precision.

## Build / infrastructure (from design.md §8)

- Scaffold from `~/models/chion` (single dependency, library + drivers,
  `precision=sp|dp`). Library target `aeros-static`, `nml` params, `ncio` +
  `variable_io` output, one driver. Add `packages/aeros.toml` to the configme registry.
- External libs per variant (`serial`/`omp`): `fesm-utils/{fftw,SHTns}/…`. SHTns is
  the key one; the Fortran wrapper is ours (prototyped in `mwm/A_scaling`).
- Open: **library-first vs standalone** (§10.8) — lean library-first given the eventual
  CLIMBER-X ocean coupling (§6.1), but not yet decided.

## Risks carried into M1

1. **Gibbs oscillations at ice margins** (§4.2) — spectral core's worst pathology,
   under-documented; margin diagnostics from day one. (Coarse core sees smoothed
   topography, so less acute at M1, but build the diagnostics now.)
2. **Vertical CFL over steep orography** may cap Δt before horizontal CFL (§4) —
   implicit/flux-limited vertical advection, smoothed orography.
3. **The correction stays hard** (b2). M1 can't de-risk this — it only builds the
   framework. M2b is where it is settled; keep the Eurasian ice sheet a first-order
   validation target (§9 risk 1), and expect to need large-scale / selected-terms /
   balanced correction (see `mwm/B_multires/B2B_OUTCOME.md`).

## Definition of done

- Dry core runs at T31 and T42 from a namelist, Float32, on the SHTns wrapper.
- Passes Held–Suarez.
- Correction framework present, exercised with a no-op term + `ΔF=0` twin + conservation
  check, ready for M2b to plug real terms in.
- T42-vs-T63 benchmark recorded; ensemble-launch convention in place.
