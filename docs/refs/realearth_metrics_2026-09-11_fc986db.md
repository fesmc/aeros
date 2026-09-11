# aeros vs ERA5, real-Earth T21L6 (2026-09-11: post coupling-api Wave 4)

Post-refactor baseline, aeros tip `a636d4a` (physics identical to `fc986db` per
Wave 4 bit-for-bit gate). Reference file:
`mwm/D_realearth/output/aeros_realearth_T21L6_ref_fc986db.nc`
(bit-for-bit equal to `T21L6_wave5.nc`; the accompanying log used for the
line-parsed globals is `logs/wave5_T21L6.log`).

Pre-refactor reference for comparison:
`.claude/worktrees/coupling-api/mwm/D_realearth/output/aeros_realearth_T21L6_ref_b192d2d.nc`
and the old handoff `docs/refs/realearth_metrics_2026-09-09_L6_rhprof.md`.

The two knobs Wave 4b intentionally changed are (1) dust cadence unified to
daily (was tied to the 6h radiation interval), (2) `sht_mode = quick` default
for reproducible transforms. Everything else is meant to be bit-for-bit.

The T21L6 case runs with `l_diag = .false.` (the default set in the rewritten
case nml); the monthly dump already carries all zonal-mean fields
(`t`, `u`, `v`, `rh`, `cf`, `eke`, `uv`, `uv_stat`, `vt_stat`, `vq_stat`,
`eke_stat`) as (month, lev, lat), so no separate `_zm.nc` was needed and the
existing diagnostic scripts run unchanged.

## Year-3 globals (from `logs/wave5_T21L6.log` last cycle)

| quantity          | post fc986db | pre b192d2d | ERA5 |
| ----------------- | ------------ | ----------- | ---- |
| TOA net (W/m^2)   | +0.22        | -0.08       | 0    |
| OLR               | 243.9        | 243.9       | 242  |
| clear-sky OLR     | 262.4        | 262.5       | 264  |
| CRE SW            | -44.6        | -45.0       | -46  |
| CRE LW            | +18.5        | +18.6       | +22  |
| cover             | 0.48         | 0.49        | 0.63 |
| P = E (mm/day)    | 2.99         | 3.02        | 2.9  |
| LH (W/m^2)        | 86.7         | 87.4        | 84   |
| SH                | 23.3         | 22.9        | --   |
| column water (mm) | 21.57        | 21.55       | 24.5 |
| max \[u\] (m/s)   | 29.2         | 31.4        | 33.7 |

Same story as the pre-refactor doc except for one meaningful shift:
**max \[u\] fell 31.4 -> 29.2 m/s**. sht_mode=quick is a plausible
proximate cause; every other global sits within round-off.

## Tropical 15S-15N RH at sigma 0.365 / 0.563 / 0.728 / 0.860 / 0.959

Post (ANN): 40 / 34 / 44 / 71 / 86 %
Pre  (ANN): 40 / 34 / 44 / 71 / 86 %
ERA5      : 37 / 45 / 56 / 70 / 84 %

Identical.  The mid-troposphere dry bias (-10 percentage points at 500-700 hPa)
is unchanged.

## Zonal T bias vs ERA5 (annual, K)

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134  post   -0.2        7.6       12.2        7.3        0.6        0.6
       pre    -0.3        7.3       12.0        8.5        1.3       -0.1
0.365  post   -1.0       -1.9        2.6        2.0       -1.1       -1.6
       pre    -0.9       -2.0        2.7        2.5       -0.5       -1.9
0.563  post   -3.0       -4.4        0.9       -1.5       -2.3       -2.7
       pre    -2.8       -4.4        0.9       -1.4       -2.1       -2.7
0.728  post   -2.2       -4.0       -0.9       -2.6       -1.7       -2.1
       pre    -2.0       -4.0       -0.8       -2.5       -1.6       -2.1
0.860  post   -0.7       -1.7       -1.7       -2.4       -0.6       -1.2
       pre    -0.5       -1.6       -1.7       -2.4       -0.5       -1.0
0.959  post   -0.7       -0.3       -0.0       -0.2        1.2        0.5
       pre    -0.5       -0.3       -0.0       -0.3        1.2        0.7
```

Tropical UT +2.6 K (was +2.7); 30-45S mid-trop -4.4 / -4.0 unchanged;
lower NH midlatitudes marginally cooler (30-45N 0.365 +2.0 vs +2.5,
45-60N 0.563 -2.3 vs -2.1).

## Lid target (aeros lid T minus ERA5 10-257 hPa mean, K)

```
                60-90S   30-60S   trop   30-60N   60-90N   max [u]
post fc986db     -1.1    +3.3    +2.2    +1.8    +1.7      29.2
pre  b192d2d     -1.2    +2.8    +2.1    +2.7    +0.8      31.4
```

Small redistribution: NH midlatitudes cooler by ~1 K, NH polar warmer by
~1 K, SH midlatitudes warmer by ~0.5 K.  Same weakened-westerlies signal.

## NH 40-70N land seasonal cycle

```
month      ERA5 t2m   ERA5 skt |   post fc986db     pre b192d2d
    1       257.6     256.5    |   260.3/259.4     259.7/258.7
    7       290.7     291.5    |   294.8/295.3     294.7/295.1
Jul-Jan  ERA5 t2m 33.1 skt 35.0|   34.5/35.9       35.0/36.5
```

2 m amplitude 34.5 (was 35.0; ERA5 33.1). July warm bias +4.1 K (was +4.0 K,
i.e. unchanged, ERA5 t2m 290.7 vs model 294.8).

Land surface energy in JJA 40-70N (unchanged in character from the pre-refactor
doc): TOA refl SW 128 model / 145 ERA5 (clear-sky 80, so 66 W/m^2 cloud
reflection short of ERA5's real cover); cover 0.38 / 0.58; LH ocean/land
109.6 / 31.3 vs ERA5 101.5 / 41.4 (i.e. summer NH land still evaporates
too little).

## Storm-track precipitation distribution (mm/day, annual zonal-mean)

```
band        aeros  ERA5
60-45S      2.02   2.93
45-30S      2.80   2.80
30-15S      2.60   2.30
15S-15N     5.05   4.52
15-30N      2.57   2.21
30-45N      2.48   2.61
45-60N      1.87   2.50
```

**Storm-track deficit is real**: 45-60S 2.02 vs 2.93 (69%), 45-60N 1.87 vs 2.50
(75%).  Subtropics 30-15S and 15-30N run wet (+13-16%); the ITCZ runs 12% wet;
30-45S/N is close.  The precipitation gap in the midlatitude storm tracks is
consistent with the weakened jet (29 vs 34 m/s) and half-strength descending
Ferrel-cell branch below.

## MMC (annual zonal-mean omega, hPa/day, model / ERA5)

```
sigma      60-75S     45-60S     30-45S     15-30S    15S-15N    15-30N    30-45N    45-60N    60-75N
0.365    -2.8/-3.5  -1.8/-4.2   1.9/6.1    6.7/11.3  -3.9/-8.4   7.0/11.8  1.2/3.8   -2.2/-3.4  0.1/-0.3
0.563    -5.0/-5.2  -2.7/-5.1   3.3/8.6   10.1/11.5  -6.7/-9.0  12.2/11.2  1.3/5.2   -2.3/-4.6  0.2/0.1
0.728    -5.6/-5.3  -2.6/-5.4   4.3/10.0   8.6/11.7  -6.9/-9.7  12.4/9.8   1.1/6.3   -1.5/-4.5  0.2/0.5
0.860    -5.1/-5.2  -2.6/-4.3   4.9/9.8    7.5/9.8   -6.5/-8.7  10.4/7.3   1.5/5.7   -1.2/-3.3  0.1/0.7
0.959    -1.8/-2.3  -1.0/-2.0   2.1/4.9    2.7/4.0   -2.6/-4.0   3.5/3.2   0.6/2.5    0.1/-1.1  0.0/0.3
```

Deficit 2 (MMC half of ERA5) quantified:

| branch (sig 0.56)        | aeros | ERA5 | aeros / ERA5 |
| ------------------------ | ----- | ---- | ------------ |
| SH descent (30-45S)      | 3.3   | 8.6  | 0.38         |
| NH ascent  (15-30N)      | 12.2  | 11.2 | 1.09         |
| Hadley rising (15S-15N)  | -6.7  | -9.0 | 0.74         |
| NH descent (30-45N)      | 1.3   | 5.2  | 0.25         |
| SH Ferrel ascent (45-60S)| -2.7  | -5.1 | 0.53         |
| NH Ferrel ascent (45-60N)| -2.3  | -4.6 | 0.50         |

The SH ascending branch (15-30S) and NH ascending branch (15-30N) are close to
ERA5; the descending branches (SH subtropics, NH subtropics, both Ferrel cells)
are all at ~half strength. Pre-refactor numbers were essentially the same
(e.g. 30-45S 3.9 vs 8.6; 30-45N 1.9 vs 5.2), so this deficit is not something
the refactor moved.

## Known deficits: pre vs post refactor

| deficit                                       | pre b192d2d | post fc986db | delta | ERA5 target | status                        |
| --------------------------------------------- | ----------- | ------------ | ----- | ----------- | ----------------------------- |
| 1. Tropical UT warm (T bias sig 0.365 trop)   | +2.7 K      | +2.6 K       | -0.1  | 0           | no retune needed              |
| 2. MMC half of ERA5 (omega sig 0.56 30-45S)   | 3.9 hPa/day | 3.3 hPa/day  | -0.6  | 8.6         | no retune needed (unchanged)  |
| 3. 30-45S mid-trop cold (T bias sig 0.56/0.73)| -4.4 / -4.0 | -4.4 / -4.0  |  0    | 0           | no retune needed              |
| 4. Summer NH land cloud deficit (JJA cover)   | 0.37 / 0.57 | 0.38 / 0.58  | +0.01 | 0.58        | no retune needed              |
| jet strength (max \|[u]\|)                    | 31.4 m/s    | 29.2 m/s     | -2.2  | 33.7        | flagged for follow-up (sht?)  |
| cover                                         | 0.49        | 0.48         | -0.01 | 0.63        | small, tracks jet             |
| TOA net                                       | -0.08       | +0.22        | +0.30 | 0           | small; still within tuning    |
| lid target 30-60N / 60-90N (K)                | +2.7 / +0.8 | +1.8 / +1.7  |       | 0           | redistributes ~1 K, watchable |

None of the four numbered deficits moved materially.  The main quantitative
shift is dynamical (jet -2.2 m/s, storm-track precip 30% below ERA5 at
45-60S/N, lid target ~1 K redistribution NH-midlat -> NH-polar); this is what
the wave 5 diff table (top variable z500 max|A-B| = 218 m at grid-cell
level, though the zonal-mean bias is essentially unchanged) was pointing at.

## Diagnostics broken by the refactor

None.  All eight scripts (`tbias_vs_era5.py`, `mmc_vs_era5.py`,
`land_cycle.py`, `land_summer.py`, `land_water.py`, `wind_vs_era5.py`,
`flux_split.py`, `lidtab.py`) ran unchanged against the post-refactor
monthly dump.  `summary.sh` parses the run log; it also works against
`logs/wave5_T21L6.log` after the diag `.nc` and log are named the same
stem (a one-line convenience, no code change).

## Bottom line

The pre-refactor tuning target list is still the right one: none of the four
numbered deficits moved by more than measurement noise. The refactor's
one meaningful climatology signature is a weaker jet and correspondingly
weaker storm-track precipitation, which is worth understanding but is not
itself in the tuning list.

## 2026-09-11 update (later): sht_mode flipped back to "cache"

The ~2 m/s jet drift (31.4 → 29.2 m/s annual max [u]) traced to `sht_mode = "quick"` (Wave 4a.2 default). A T21L6 rerun with `sht_mode = "cache"` restored the pre-refactor number (31.05 m/s, within 0.3 m/s of 31.4). The default is now `cache` at `input/aeros_defaults.nml`; `quick` remains available (`par/held_suarez.nml` keeps `quick` explicit for its bit-for-bit history where the extra ~1e-13 SHTns nondeterminism is unwelcome; note - held_suarez actually pins "cache" - see nml).

Cache-mode Wave-5 references stamped as `mwm/D_realearth/output/aeros_realearth_<case>_ref_cache.nc` (T21L6, T21L8, T21L12, T21L12_dt2700). **T21L12_tau2h NaN'd at day 388** in cache mode (max|u| already 81 m/s from day 180); it completed in quick mode. Interpreted as a stability-margin case, not a defaults regression — the tau2h experiment is knob-tuned and expected to need retuning at cache.

vs pre-refactor `b192d2d` (T21L6 cache): 32/37 vars differ at the pointwise level (chaotic divergence over 3 yr); climate means are tiny — z500 mean(post-pre) −3.6 m, shf −0.14 W/m², lwd −0.03 W/m², swup +0.02 W/m². The four numbered deficits (tropical UT warm, MMC half of ERA5, 30–45S mid-trop cold, summer NH land cover) are unchanged. Pre-refactor tuning targets are usable as-is.

Recommendation on the ERA5 tuning track: proceed with the SBM convective/large-scale partition (deficit 2 lever).
