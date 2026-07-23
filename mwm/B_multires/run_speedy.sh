#!/bin/bash
# runme exe alias "speedy" -> this wrapper. Takes the config path as $1.
# The runme rundir differs from this dir, so use an absolute --project path.
set -e

PROJDIR=/albedo/work/projects/p_forclima/robinson/models/aeros/mwm/B_multires

# runme invokes this from the rundir with the config as a RELATIVE basename
# (e.g. "b1.toml"). Resolve it to an absolute path BEFORE we cd away, else the
# driver can't find it once cwd changes to PROJDIR.
CFG="$1"
case "$CFG" in
  /*) : ;;
  *)  CFG="$(readlink -f "$CFG")" ;;
esac

cd "$PROJDIR"

# Honour the SLURM cpus-per-task for Julia threading. SpeedyWeather runs on the
# KernelAbstractions CPU() backend, which parallelises its kernels and spectral
# transforms over Julia threads -- so JULIA_NUM_THREADS is what controls the
# on-node thread count (no special model/architecture setting is required for
# multi-threaded CPU execution).
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

exec julia --project="$PROJDIR" "$PROJDIR/b1_probe.jl" "$CFG"
