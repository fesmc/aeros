#!/bin/bash
# Gather results.txt from all Option-A sweep rundirs into one CSV table.
# Run from mwm/ :  ./collect_A.sh > A_results.csv
set -eu
cd "$(dirname "$0")"

echo "variant,label,lmax,nlev,nlat,nphi,nlm,nthreads,precision,t_step_s,t_transform_frac,t_gridpoint_frac,t_solve_frac,core_s_per_year,model_years_per_day"
for d in runs/A_*/; do
  r="$d/results.txt"
  [[ -f "$r" ]] || continue
  # prefer the machine-parseable block (last match wins); default variant for v1
  get() { grep -E "^[[:space:]]*$1[[:space:]]*=" "$r" | tail -1 | sed 's/.*=[[:space:]]*//' | tr -d ' '; }
  variant="$(get variant)"; [[ -n "$variant" ]] || variant="shtns"
  echo "${variant},$(get label),$(get lmax),$(get nlev),$(get nlat),$(get nphi),$(get nlm),$(get nthreads),$(get precision),$(get t_step_s),$(get t_transform_frac),$(get t_gridpoint_frac),$(get t_solve_frac),$(get core_s_per_year),$(get model_years_per_day)"
done
