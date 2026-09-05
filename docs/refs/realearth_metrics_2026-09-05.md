# aeros vs SpeedyWeather vs ERA5, real-Earth T21

aeros: `mwm/D_realearth/output/aeros_realearth_T21L12.nc`  
SW: `/Users/alrobi001/models/aeros/mwm/D_realearth/output/sw_realearth_T21L8.nc`

| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 1.8 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 217.9 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 287.5 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 60.5 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 37.8 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 67.0 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 30.3 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 1573.9 | 1037.9 | NaN |
