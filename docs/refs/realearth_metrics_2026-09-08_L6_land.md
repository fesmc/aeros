# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, 2026-09-08 night: land winter recalibration + 2 m diagnostic)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600, 2694 SYPD (95 s). On top of 035ec69 (aeros_bd): `surf_ri_crit = 1` (Frierson's stable limit), `beta_frozen = 0.1` (frozen-ground evaporation), and the monthly `t2m` is now a surface-layer 2 m diagnostic (the lowest level is `tlow`). SpeedyWeather's t2m in the table is still its lowest level.

Year-3 global means: TOA net +0.2 W/m², OLR 240.8 (ERA5 242), clear-sky OLR 257.8 (264), CRE SW −50.1 / LW +17.1 (−46 / +22), cover 0.55 (0.63), P = E = 3.08 mm/day (2.9), LH 89 (84; ocean 106 / land 49 vs 101 / 41), SH 17.6 (17; ocean 13 / land 29 vs 12 / 29), column water 24.9 mm (24.5).

NH 40–70N land (`diag/land_cycle.py`, `land_budget.py`): 2 m T Jan 261.4 / Jul 293.6 (ERA5 257.6 / 290.7), skin 260.8 / 293.9 (skt 256.5 / 291.5); Jul−Jan 2 m 32.2, skin 33.1 (ERA5 33.1 / 35.0). January budget: LW down 233 (ERA5 209), SH −12 (−8), LH 8 (4), cover 0.68 (0.70), skin−2 m −0.6 (−1.2). July: LW down 349 (341), SW net 225 (196), SH 43 (40), LH 100 (73), cover 0.47 (0.57).

Lid T minus the ERA5 10–257 hPa layer mean: 60–90S +0.9, 30–60S +2.2, tropics +2.5, 30–60N +2.6, 60–90N +1.2.

Zonal T bias vs ERA5 at the same pressure (annual, K; 60–75S / 30–45S / 15S–15N / 30–45N / 45–60N / 60–75N):

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134          1.5        6.4       12.4        8.4        1.1        0.2
0.365         -1.3       -2.2        2.8        2.6       -0.1       -1.3
0.563         -3.3       -4.1        1.3       -0.8       -1.9       -3.1
0.728         -1.9       -3.7       -0.4       -2.0       -1.3       -2.2
0.860         -0.1       -1.3       -1.0       -1.8       -0.2       -1.2
0.959          0.0       -0.3        0.3       -0.1        1.1        0.2
```

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc` · SW: `mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 3.1 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 216.8 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 285.5 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 33.5 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 31.0 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 71.6 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 32.8 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2694.1 | 708.7 | NaN |
