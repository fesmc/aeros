# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, 2026-09-08 late: ice-phase condensation and cloud thresholds)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600, 2694 SYPD. On top of 2c7d59d: `cond_rh_crit_ice = 0.80`, `cld_rhc_ice = 0.60` (new keys; thresholds blended by the liquid fraction), `cld_lwc 0.034e-3`.

Year-3 global means: TOA net +0.2 W/m², OLR 235.8 (ERA5 242), clear-sky OLR 258.2 (264), CRE SW −55.0 / LW +22.3 (−46 / +22), cover 0.50 (0.63), P = E = 2.95 mm/day (2.9), LH 85.6 (84; ocean 102 / land 47 vs 101 / 41), SH 16.0 (17), column water 24.9 mm (24.5).

Winter cap DJF RH at 60–75N: 73/75/75/76/86 % at σ 0.36–0.96 (ERA5 70/68/69/74/84; before 84–94). NH 40–70N land: January LW down 219 (ERA5 209; before 233), 2 m Jan 259.4 / Jul 293.9 (257.6 / 290.7), Jul−Jan 2 m 34.5, skin 35.9 (33.1 / 35.0); July SW net 217 (196), LH 93 (73), cover 0.47 (0.57).

Lid T minus the ERA5 10–257 hPa layer mean: 60–90S +1.9, 30–60S +3.0, tropics +1.8, 30–60N +2.0, 60–90N +2.5.

Zonal T bias vs ERA5 at the same pressure (annual, K; 60–75S / 30–45S / 15S–15N / 30–45N / 45–60N / 60–75N):

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134          2.5        7.0       11.8        7.5        1.0        1.4
0.365          0.3       -1.9        3.0        2.5       -0.1       -0.1
0.563         -2.5       -3.7        1.7       -0.5       -1.6       -2.8
0.728         -1.8       -3.4       -0.0       -1.7       -1.2       -2.4
0.860         -0.4       -1.1       -0.7       -1.5       -0.2       -1.4
0.959         -0.7       -0.1        0.6        0.1        1.1       -0.1
```

DJF 69N-band −0.5/−3.5/−3.5/−2.5/−0.4 at σ 0.36–0.96.

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc` · SW: `mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 3.0 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 216.7 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 285.6 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 31.1 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 29.1 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 73.1 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 35.5 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2694.1 | 708.7 | NaN |
