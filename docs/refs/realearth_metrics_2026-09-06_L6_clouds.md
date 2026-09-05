# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, mixed-phase saturation + coupled cloud calibration, 2026-09-06)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600. Year-3 means: TOA net +1.3 W/m², absorbed SW 288 clear / 231 all-sky, OLR 261 clear / 230 all-sky, CRE SW −57 LW +31 (ERA5 −46 / +22), cloud cover 0.52 (ERA5 0.63), P = E = 2.72 mm/day, LH 79 W/m², water residual −0.010 mm/day.

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc`  
SW: `/Users/alrobi001/models/aeros/mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 2.7 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 225.1 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 283.3 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 56.7 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 47.1 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 93.1 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 36.9 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2830.6 | 1037.9 | NaN |
