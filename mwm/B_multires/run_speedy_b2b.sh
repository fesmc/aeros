#!/bin/bash
# runme wrapper for the Option-B stage-b2 APPROACH B probe (nudging-derived <N>).
# Analogous to run_speedy_b2.sh but points at b2b_probe.jl. Takes the config path
# as $1.
#
# NB: the runme rundir differs from this dir, so we resolve the config to an
# ABSOLUTE path BEFORE any cd, then cd to the project dir and use an absolute
# --project.
set -e

PROJDIR=/albedo/work/projects/p_forclima/robinson/models/aeros/mwm/B_multires

# Resolve the (possibly relative) config to an absolute path BEFORE cd'ing away.
CFG="$1"
case "$CFG" in
  /*) : ;;
  *)  CFG="$(readlink -f "$CFG")" ;;
esac

cd "$PROJDIR"

# Honour SLURM cpus-per-task for Julia threading.
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

exec julia --project="$PROJDIR" "$PROJDIR/b2b_probe.jl" "$CFG"
