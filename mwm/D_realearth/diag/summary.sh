#!/bin/bash
# usage: summary.sh run1 run2 ...  (names under mwm/D_realearth/diag, logs in logs/)
printf "%-6s %6s %6s %6s %6s %6s %6s %5s %5s %6s | tropics T bias sig .36/.56/.73/.86/.96\n" run TOA OLR OLRcs CRESW CRELW cover P W LH
for n in "$@"; do
  l=logs/$n.log; [ -f mwm/D_realearth/diag/$n.nc ] || { echo "$n: not finished"; continue; }
  toa=$(grep "^   TOA: SWin" $l | tail -1); cs=$(grep "^   TOA clear" $l | tail -1); w=$(grep "^   WATER: W" $l | tail -1)
  net=$(echo "$toa" | awk '{print $7}'); olr=$(echo "$toa" | awk '{print $5}'); lh=$(echo "$toa" | awk '{print $14}')
  olrcs=$(echo "$cs" | awk '{print $6}'); sw=$(echo "$cs" | awk '{print $10}'); lw=$(echo "$cs" | awk '{print $12}'); cov=$(echo "$cs" | awk '{print $NF}')
  p=$(echo "$w" | awk '{print $7}'); W=$(echo "$w" | awk '{print $3}')
  tb=$(python3 mwm/D_realearth/diag/tbias_vs_era5.py mwm/D_realearth/diag/$n.nc | sed -n 5,9p | awk '{printf "%s/", $4}')
  printf "%-6s %6s %6s %6s %6s %6s %6s %5s %5s %6s | %s\n" $n $net $olr $olrcs $sw $lw $cov $p $W $lh $tb
done
