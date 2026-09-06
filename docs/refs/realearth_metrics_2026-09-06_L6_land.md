# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, land recalibration 2026-09-06: seasonal ERA5 albedo, c_soil 5e6, snow bump off)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600. Year-3 means: TOA net +2.1 W/m², CRE SW −58 LW +32, cloud cover 0.52, P = E = 2.70 mm/day, LH 79 W/m².

NH-land lowest-level T vs ERA5 t2m by band (DJF / JJA, K): 40–50N 263.4/293.3 vs 266.7/293.9; 50–60N 253.0/289.7 vs 259.0/288.6; 60–70N 241.3/280.7 vs 248.8/282.8; 70–85N 236.0/269.4 vs 242.4/271.2. ERA5 DJF at 925 hPa over 60–70N land is 255 K: the winter bias at comparable height is ~13 K and is the whole lower troposphere (zonal mean, 700–925 hPa, 7–13 K).

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc`  
SW: `/Users/alrobi001/models/aeros/mwm/D_realearth/output/sw_realearth_T21L8.nc`
| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 2.7 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 224.8 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 283.2 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 55.9 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 46.5 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 79.9 | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 38.6 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2724.1 | 1037.9 | NaN |
