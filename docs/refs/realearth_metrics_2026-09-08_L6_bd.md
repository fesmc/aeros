# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, 2026-09-08 evening: prescribed lid Brewer-Dobson circulation)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600, 2694 SYPD (104 s). On top of `..._2026-09-08_L6.md` (b25f96e): `aeros_bd` with `input/bd_wstar_era5.nc` (effective w*(lat,month) backed out of the run relaxed toward the ERA5 10–257 hPa layer-mean lid, ×0.85 calibration; `bd_scale 1`, `bd_per_K 0`), `cld_lwc 0.029e-3`.

Year-3 global means: TOA net −0.3 W/m², OLR 240.1 (ERA5 242), clear-sky OLR 257.8 (264), CRE SW −51.3 / LW +17.7 (−46 / +22), cover 0.57 (0.63), P = E = 3.15 mm/day (2.9), LH 91 (84; ocean 107 / land 53 vs 101 / 41), SH 15.3 (17), column water 24.8 mm (24.5).

Lid T minus the ERA5 10–257 hPa layer mean (annual): 60–90S 0.0, 30–60S +2.9, tropics +2.5, 30–60N +2.0, 60–90N +0.5 K (before the term: −21 / −3 / +5 / 0 / −21).

Zonal T bias vs ERA5 at the same pressure (annual, K; bands 60–75S / 30–45S / 15S–15N / 30–45N / 45–60N / 60–75N):

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134          0.6        7.3       12.5        7.5        0.8       -0.6
0.365         -1.4       -1.9        2.7        2.3       -0.5       -1.2
0.563         -2.8       -4.2        1.3       -1.0       -2.3       -3.3
0.728         -1.5       -3.8       -0.4       -2.1       -1.6       -2.5
0.860          0.1       -1.4       -1.1       -2.0       -0.7       -1.6
0.959         -0.5       -0.4        0.2       -0.4        0.2       -0.9
```

σ 0.365 at 30–45N +2.3 (morning +4.5, yesterday +7.7); tropics +2.7. DJF 69N-band −1.7/−3.6/−3.5/−3.0/−2.3 at σ 0.36–0.96.

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc` · SW: `mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 3.1 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 216.7 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 285.2 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 31.5 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 29.4 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 70.7 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 30.7 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2694.1 | 708.7 | NaN |
