# aeros vs SpeedyWeather vs ERA5, real-Earth T21

aeros: `mwm/D_realearth/output/aeros_realearth_T21L8.nc`  
SW: `/Users/alrobi001/models/aeros/mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 1.9 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 222.7 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 285.8 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 55.6 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 41.0 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 78.2 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 32.5 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2208.2 | 1037.9 | NaN |
