#!/bin/bash
# runme wrapper for the Option-B stage-b2 correction probe. Analogous to
# run_speedy.sh but points at b2_probe.jl. Takes the config path as $1.
#
# NB: the runme rundir differs from this dir, so we resolve the config to an
# ABSOLUTE path BEFORE any cd (run_speedy.sh's config handling, done right), then
# cd to the project dir and use an absolute --project.
set -e

PROJDIR=/albedo/work/projects/p_forclima/robinson/models/aeros/mwm/B_multires

# Resolve the (possibly relative) config to an absolute path BEFORE cd'ing away.
CFG="$1"
case "$CFG" in
  /*) : ;;
  *)  CFG="$(readlink -f "$CFG")" ;;
esac

cd "$PROJDIR"

# Honour SLURM cpus-per-task for Julia threading (SpeedyWeather threads its
# KernelAbstractions CPU() kernels/transforms over Julia threads).
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

exec julia --project="$PROJDIR" "$PROJDIR/b2_probe.jl" "$CFG"
