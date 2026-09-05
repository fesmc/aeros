# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6 default, 2026-09-06)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600, TOA yr1/yr2/yr3 −3.0/−3.7/+3.2 W/m², P = E = 2.2 mm/day, water residual −0.007 mm/day.

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc`  
SW: `/Users/alrobi001/models/aeros/mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 2.2 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 225.8 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 283.5 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 52.0 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 43.5 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 85.1 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 31.9 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2810.4 | 1037.9 | NaN |
