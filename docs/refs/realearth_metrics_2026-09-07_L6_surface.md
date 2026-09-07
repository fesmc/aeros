# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, end of 2026-09-07 evening: Monin-Obukhov surface exchange + cloud re-balance)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600, 2694 SYPD (97 s). Commit after 6b8b55a (surface) + this cloud commit. Defaults: `cld_rhc_top = 1.03`, `cld_lwc = 0.031e-3`, `cld_iwc = 0.014e-3`, `vert_vanleer = .true.`, `conv_rh_ref = 0.6`, `conv_rh_ref_bl = 0.85`, bulk coefficients from the log law with the true skin-air static-energy jump and the Businger-Dyer unstable branch.

Year-3 global means: TOA net +0.1 W/m², OLR 236.3 (ERA5 242), clear-sky OLR 257.4 (264), clear-sky absorbed 291, CRE SW −54.8 / LW +21.2 (ERA5 −46 / +22), cloud cover 0.57 (0.63), P = E = 3.11 mm/day (2.9), LH 90 (84; ocean 108 / land 48 vs ERA5 101 / 41), SH 16.5 (17; ocean 12 / land 27 vs 12 / 29), column water 25.3 mm (24.5).

Zonal T bias vs ERA5 at the same pressure (annual, K; bands 60–75S / 30–45S / 15S–15N / 30–45N / 45–60N / 60–75N; `mwm/D_realearth/diag/tbias_vs_era5.py`):

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134        -19.6       12.5       20.8       17.4       -5.9      -21.3
0.365         -2.5        3.0        3.6        7.7        2.4       -1.9
0.563         -1.4       -3.5        1.2       -0.5       -2.2       -1.7
0.728         -0.7       -4.2       -0.3       -3.2       -2.7       -1.2
0.860          0.6       -1.9       -0.8       -3.3       -2.0       -0.7
0.959         -0.2       -0.8        0.5       -1.5       -1.2       -0.1
```

DJF: 69N-band (60–75N) bias −2.3/−1.6/−1.4/−1.3/−0.7 K at σ 0.36–0.96 (start of day −4.1/−8.0/−6.1/−5.1/−3.8); 30–45N −0.9/−4.0/−3.9/−2.9 at σ 0.56–0.96; tropics +3.4 at σ 0.365, within 1 K below. Remaining: warm σ 0.365 outside the tropics (+7.7 at 30–45N annual) over a cold 0.56–0.86 at 30–45S/N (−3 to −4): the static-stability item. Lid unchanged (−20 K polar, +17..+21 subtropics/tropics).

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc` · SW: `mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 3.1 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 218.4 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 285.0 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 56.9 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 47.1 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 101.6 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 33.1 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2694.1 | 708.7 | NaN |
