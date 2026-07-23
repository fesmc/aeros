#!/bin/bash
# Submit the M0a OpenMP scaling sweep for BOTH transform variants:
#   shtns   -> aeros_bench.x          (SHTns per-field transform)
#   batched -> aeros_bench_batched.x  (custom batched-DGEMM transform, §4)
# Each point runs on a full EXCLUSIVE node so the scaling curve is clean.
# Run from mwm/ :  ./run_sweep_A.sh [variant ...]   (default: both)
# Results: runs/A_<variant>_<cfg>_omp<N>/results.txt
set -eu
cd "$(dirname "$0")"

declare -A EXE=( [shtns]=bench [batched]=bench_batched )
VARIANTS=("${@:-shtns batched}")
# allow "shtns batched" as one arg or separate args
VARIANTS=(${VARIANTS[@]})

CFGS=(T31L16 T42L19)
THREADS=(1 2 4 8 16 32 64 128)
QUEUE=30min

for v in "${VARIANTS[@]}"; do
  exe="${EXE[$v]:-}"
  [[ -n "$exe" ]] || { echo "unknown variant: $v (use shtns|batched)" >&2; exit 1; }
  for cfg in "${CFGS[@]}"; do
    par="A_scaling/par/bench_${cfg}.nml"
    [[ -f "$par" ]] || { echo "MISSING par: $par" >&2; exit 1; }
    for omp in "${THREADS[@]}"; do
      out="runs/A_${v}_${cfg}_omp${omp}"
      echo ">>> submit $v $cfg omp=$omp -> $out"
      runme -e "$exe" -n "$par" --omp "$omp" -o "$out" -q "$QUEUE" -s -r
    done
  done
done
echo "Submitted. Watch: squeue --me ; collect: ./collect_A.sh"
