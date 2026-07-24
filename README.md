# aeros

A fast global atmospheric model for coupled climate–ice-sheet simulation on
10⁴–10⁶ yr timescales, with best performance over polar domains.

Spectral primitive-equation core (T31–T42, L16–20) built as the static library
`libaeros.a`, OpenMP-only, intended to couple to CLIMBER-X's ocean and to Yelmo.

Status: **dry dynamical core integrates (M1.4)**. Hybrid σ–p coordinate,
(ζ,D)↔(u,v), primitive-equation tendencies by the transform method, and a
semi-implicit leapfrog with a RAW time filter and ∇⁶ diffusion. `aeros_update`
advances the model. What is missing is forcing: Held–Suarez is M1.5, moisture
and radiation M2, so a run from the shipped initial condition correctly does
nothing at all.

## Docs

- [docs/design.md](docs/design.md) — the design plan: what aeros is for, what
  is settled, what is still open, and the milestone sequence.
- [docs/m0a_results.md](docs/m0a_results.md) — the measured cost of the
  dynamical core, and why aeros is double precision.
- [docs/m1_results.md](docs/m1_results.md) — what the core does once it runs:
  stability, conservation and the numbers M1.5 will need.

## Install

`aeros` is a configme package. From a machine with
[configme](https://github.com/fesmc/configme) installed:

```bash
configme install aeros
```

That clones aeros and `fesm-utils`, generates the repo-root `Makefile` for your
machine and compiler, and builds the dependencies. In an existing checkout,
`configme -m <machine> -c <compiler>` regenerates the Makefile alone.

`fesm-utils` supplies all three external dependencies: `nml`/`ncio`/`coords`,
the **SHTns** spherical-harmonic transform library, and **FFTW** beneath it.

## Build

```bash
make aeros-static
```

builds `libaeros.a` (OpenMP by default; `make openmp=0` for serial). `make run`
builds the standalone driver, `make tests` the acceptance tests, `make all`
both.

aeros is **double precision throughout**, with no build-time switch. That
reverses the design plan's "Float32 throughout the core", and it was settled by
measurement: SHTns has no single-precision CPU path, so a single-precision core
converts up and back at every transform, and transforms are 86–97% of a
dynamics timestep. Measured, `sp` was ~17% *slower*. See
[docs/m0a_results.md](docs/m0a_results.md) §5. Single precision survives as
`wp_ext`, the kind used at the coupling boundary.

## Run

```bash
runme -r -o output/test -n par/aeros.nml
```

stages a run directory and runs it in the background; add `-s -q <queue>` to
submit to SLURM instead. See [runme](https://github.com/fesmc/runme).
