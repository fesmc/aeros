# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, 2026-09-08: pressure-scaled near-IR distribution + lid cold trap)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600, 2694 SYPD (101 s). Commits bb6ac9c (ecCKD near-IR distributed on a pressure-scaled path) and b25f96e (lid cold trap, `cond_lid_q_max = 1e-5`) on top of the 2026-09-07 evening state (`..._2026-09-07_L6_surface.md`); clouds unchanged (`cld_rhc_top 1.03`, `cld_lwc 0.031e-3`).

Year-3 global means: TOA net −0.5 W/m², OLR 237.8 (ERA5 242), clear-sky OLR 256.9 (264), CRE SW −53.9 / LW +19.1 (−46 / +22), cover 0.57 (0.63), P = E = 3.21 mm/day (2.9), LH 93 (84; ocean 110 / land 52 vs 101 / 41), SH 15.9 (17), column water 24.7 mm (24.5). Lid vapour 0.01–0.02 g/kg (was 0.13–0.14; ERA5 0.00–0.01).

Zonal T bias vs ERA5 at the same pressure (annual, K; bands 60–75S / 30–45S / 15S–15N / 30–45N / 45–60N / 60–75N):

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134        -20.3        7.8       15.2       11.3       -9.1      -21.8
0.365         -3.0        0.4        3.0        4.5        0.2       -2.4
0.563         -1.8       -4.0        1.4       -0.9       -2.4       -1.5
0.728         -0.9       -4.2       -0.3       -2.7       -2.3       -1.0
0.860          0.3       -1.9       -1.0       -2.7       -1.4       -0.6
0.959         -0.3       -0.8        0.3       -1.0       -0.7       -0.2
```

σ 0.365 at 30–45N +4.5 (was +7.7 the evening before; near-IR distribution −1.3, lid trap −2); lid 30–45N +11 (was +17); tropics +3.0 at σ 0.365, within 1.5 K below. DJF 69N-band −3.2/−2.0/−2.1/−2.3/−2.1 at σ 0.36–0.96.

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc` · SW: `mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 3.2 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 213.9 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 285.0 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 52.4 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 44.0 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 84.6 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 32.7 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2694.1 | 708.7 | NaN |
