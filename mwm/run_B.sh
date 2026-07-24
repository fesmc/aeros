#!/bin/bash
# Submit the Option-B stage-b1 multi-resolution probe (SpeedyWeather T85 vs T31).
# One heavy job on a full exclusive node. Run from mwm/ :  ./run_B.sh [config]
#   config defaults to B_multires/par/b1.toml
set -eu
cd "$(dirname "$0")"

cfg="${1:-B_multires/par/b1.toml}"
if [[ ! -f "$cfg" ]]; then echo "MISSING config: $cfg" >&2; exit 1; fi

# 12h queue: T85 nlev=8 for a few hundred days is the cost driver. Exclusive node
# regardless; 32 threads avoids the >~8-thread transform collapse measured in A
# (SpeedyWeather T85 has more per-transform work, so it should thread better than
# the A T42 case, but 128 risks oversubscription).
runme -e speedy -n "$cfg" --omp 32 -o runs/B_b1 -q 12h -s -r
echo "B job submitted. Watch: squeue --me ; results in runs/B_b1/ and B_multires/output/"
