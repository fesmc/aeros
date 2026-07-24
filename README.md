# aeros

A fast global atmospheric model for coupled climate–ice-sheet simulation on
10⁴–10⁶ yr timescales, with best performance over polar domains.

Spectral primitive-equation core (T31–T42, L16–20) built as the static library
`libaeros.a`, OpenMP-only, intended to couple to CLIMBER-X's ocean and to Yelmo.

Status: **scaffold (M0)**. The build, the SHTns wrapper, the Gaussian grid and
the public interface are in place; the dynamical core is not.

## Docs

- [docs/design.md](docs/design.md) — the design plan: what aeros is for, what
  is settled, what is still open, and the milestone sequence.

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

`make precision=dp` builds a double-precision variant into its own
`libaeros/{include,bin}-dp/`, so it coexists with the default single-precision
build. Both are real configurations, not a production setting plus a debug
toggle — SHTns' interface is double, so `sp` trades a copy-convert at every
transform against half the memory traffic everywhere else. Which wins is what
M0a measures.

## Run

```bash
runme -r -o output/test -n par/aeros.nml
```

stages a run directory and runs it in the background; add `-s -q <queue>` to
submit to SLURM instead. See [runme](https://github.com/fesmc/runme).
