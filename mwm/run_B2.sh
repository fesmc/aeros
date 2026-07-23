#!/bin/bash
# Submit the Option-B stage-b2 correction probe (SpeedyWeather: T85 ref, bare T31,
# and T31+ΔF). One heavy job. Run from mwm/ :  ./run_B2.sh [config]
#   config defaults to B_multires/par/b2.toml
set -eu
cd "$(dirname "$0")"

cfg="${1:-B_multires/par/b2.toml}"
[[ -f "$cfg" ]] || { echo "MISSING config: $cfg" >&2; exit 1; }

# b2 runs three integrations (one T85 is the cost driver) → 12h queue, 32 threads
# (avoids the >~8-thread spectral-core collapse measured in A).
runme -e speedy_b2 -n "$cfg" --omp 32 -o runs/B_b2 -q 12h -s -r
echo "B2 job submitted. Watch: squeue --me ; results in B_multires/output_b2/"
