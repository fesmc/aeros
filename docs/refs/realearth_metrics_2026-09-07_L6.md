# aeros vs SpeedyWeather vs ERA5, real-Earth T21 (aeros L6, end of 2026-09-07: transport sign fix, Betts sub-cloud reference, near-IR SW re-anchoring, cloud recalibration)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600, 2700 SYPD (96 s). Commit 9e7b227. Defaults: `vert_vanleer = .true.`, `cld_rhc_top = 1.40`, `cld_lwc = 0.035e-3`, `conv_rh_ref = 0.6`, `conv_rh_ref_bl = 0.85`, ecCKD near-IR SW_A1 0.28 / SW_B2 0.033.

Year-3 global means: TOA net +2.2 W/m², absorbed SW 246.6 (ERA5 243), OLR 244.4 (242), clear-sky OLR 252 (264), CRE SW −43.9 / LW +7.8 (ERA5 −46 / +22), cloud cover 0.52 (0.63), P = E = 2.80 mm/day (2.9), LH 81 (84), SH 21 (17), column water 20.6 mm (24.5).
Atmospheric budget: solar absorbed by the atmosphere 68 W/m² (ERA5 79; clear-sky columns 78 vs 76.5), LW cooling 170 (185).

Zonal T bias vs ERA5 (annual, K): 500 hPa −2.8 (tropics) / −3.7 (47N) / −6.5 (58N); 700–850 hPa −3 to −7 everywhere; 300 hPa +2 (tropics) / +6 to +8 (30–47N); lid (σ 0.13) −20 at 60N, +11 JJA 41N. RH at 500 hPa 62–76 % vs ERA5 47–60.

DJF: transport into the cap >60N 1.45 PW (ERA5 2.19), eddy latent at 58N 1.04 PW (before the transport fix 0.00), eddy sensible 6.1 PW at 58N (strong); 69N T bias −4.1/−8.0/−6.1/−5.1/−3.8 at σ 0.36–0.96 (before the fix −4.4/−5.9/−7.2/−9.9/−12.3); cap OLR 161 (177), cap cover 0.84, RH 81–102 % at 64N (ERA5 66–85). Tracer mass-rescaling term r: rms 0.4 %/step (was 2.5 %), polar row 1.4 %.

Tropical ocean: lowest-level RH 85–88 % (ERA5 85–88), T 292.6 (2 m 298.7, ~295 mixed down), LH 131 (132).

aeros: `mwm/D_realearth/output/aeros_realearth_T21L6.nc`  
SW: `mwm/D_realearth/output/sw_realearth_T21L8.nc` (with the eddy diagnostics, 709 SYPD)
| metric | aeros | SpeedyWeather | ERA5 |
|---|---|---|---|
| global-mean precip [mm/day] | 2.8 | 3.3 | 2.9 |
| global-mean T, top level / 70 hPa [K] | 215.1 | 235.2 | 206.6 |
| global-mean T, lowest level / 2 m [K] | 282.4 | 278.3 | 287.5 |
| max annual zonal-mean u [m/s] | 53.1 | 63.0 | 33.7 |
| max annual [u] at 200 hPa [m/s] | 45.1 | 41.2 | 31.3 |
| DJF NH Z500* RMS 20-80N [m] | 106.1 (78–106 across same-physics runs: 1-yr noise) | 68.3 | 76.1 |
| t2m |Jul-Jan|, NH 40-70N land [K] | 35.2 | 29.2 | 33.6 |
| throughput [SYPD, 10 threads] | 2700 | 709 | NaN |
