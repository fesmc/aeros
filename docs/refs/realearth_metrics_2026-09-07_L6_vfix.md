# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, 2026-09-07: humidity-transport sign fix d2f96fc + recalibration: vert_vanleer, cld_rhc_top 1.40, cld_lwc 0.035e-3)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600. Year-3 means: TOA net +1.6 W/m², absorbed SW 241 (ERA5 243), OLR 240 (242), clear-sky OLR 250 (264), CRE SW −47.5 LW +10 (ERA5 −46/+22), cloud cover 0.45 (0.63), P = E = 3.00 mm/day, LH 87 W/m².

DJF northward energy transport into the cap >60N 2.18 PW (ERA5 2.19), latent 0.62 (0.55); eddy latent transport at 58N 0.96 PW (was 0.00 before the fix). DJF T bias at 69N: −4.2/−8.3/−6.2/−5.2/−4.4 K at σ 0.36–0.96 (was −4.4/−5.9/−7.2/−9.9/−12.3 before the fix); the free troposphere is 4–10 K cold at 500–850 hPa at all latitudes and seasons — the dominant remaining bias (see m2_handoff.md 2026-09-07). Winter cap: OLR 160 (177), cover 0.82, RH 82–101% at all levels at 64N (ERA5 66–85).


| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 3.0 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 213.6 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 282.8 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 52.0 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 43.7 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 82.0 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 34.8 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2694.1 | 708.7 | NaN |
