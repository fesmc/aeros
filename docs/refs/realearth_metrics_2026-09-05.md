# aeros vs SpeedyWeather vs ERA5, real-Earth T21

aeros: `/Users/alrobi001/models/aeros/mwm/D_realearth/output/aeros_realearth_T21L12_dt2700.nc`  
SW: `/Users/alrobi001/models/aeros/mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 2.1 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 214.2 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 288.1 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 63.1 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 36.5 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 46.8 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 29.2 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 1605.8 | 1016.4 | NaN |
