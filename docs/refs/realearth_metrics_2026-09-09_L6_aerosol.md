# aeros vs ERA5, real-Earth T21L6 (2026-09-09: CLIMBER-X aerosol suite and land surface resistance)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600 (run `IEl032`). On top of 11ddd15: background aerosol `sw_aer_ot 0.12 / sw_aer_im 0.02` with the ecCKD near-IR H2O fit back at SESAM's (ac75600); sulfate direct effect from the CMIP6 SO4 load, year 2003 (9186c70); `surf_r_s_land 100` s/m (4edcfd0); dust optical thickness maps, `dust_scale 0.3 / dust_im 0.01` (e89a263); sulfate indirect effect on the liquid cloud water (this commit); `cld_lwc 0.032e-3`.

Year-3 global means: TOA net +0.08 W/m², OLR 235.6 (ERA5 242), clear-sky OLR 258.1 (264), CRE SW −53.4 / LW +22.5 (−46 / +22), cover 0.48 (0.63), P = E = 2.85 mm/day (2.9), LH 82.7 (84), column water 24.5 mm (24.5), max [u] 31.1 m/s (33.7).

Offline clear-sky SW (validate_era5, ecCKD, annual, model − ERA5): global surface down +1.1 / TOA net −4.9; 40–60S ocean −0.3 / −0.5; Sahara −7 / +0.6; E Asia +11 / −4; Europe +15 / −2; Siberia +11 / −12 (annual-mean albedo artefact at TOA). Indirect effect: global SW CRE −0.3 (ERA5 cloud), −1.1 (diagnosed cloud).

Lid T minus the ERA5 10–257 hPa layer mean: 60–90S +1.4, 30–60S +2.9, tropics +1.7, 30–60N +2.2, 60–90N +2.4.

Zonal T bias vs ERA5 (annual, K; 60–75S / 30–45S / 15S–15N / 30–45N / 45–60N / 60–75N):

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134          2.0        7.0       11.6        7.6        1.2        1.2
0.365          0.2       -2.0        3.0        2.5        0.0       -0.2
0.563         -3.0       -4.1        1.5       -0.9       -2.0       -2.8
0.728         -2.3       -3.7       -0.3       -2.1       -1.6       -2.4
0.860         -0.7       -1.2       -0.9       -1.9       -0.6       -1.7
0.959         -0.9       -0.2        0.5       -0.1        0.8       -0.2
```

NH 40–70N land (model / ERA5): water budget now on ERA5 in spring (E April–June 0.98/1.68/2.18 vs 0.96/1.65/2.40 mm/day; P 1.37/1.71/2.26 vs 1.52/1.88/2.31), Jul–Aug P low (1.96/1.58 vs 2.49/2.37). Energy: July LW down 348 (341), SW net 228 (196), SH 74 (40), LH 61 (73), cover 0.37 (0.57); JJA TOA reflected 138 (145). 2 m: January 259.4-ish (257.6), July 295.1 (290.7); the summer warm bias is +4 to +5 K, radiative (clear-sky continental aerosol not represented by dust/sulfate + the diagnostic cloud's summer-land deficit).

```
month      ERA5 t2m   ERA5 skt |    IEl032.nc
    1       257.6     256.5 |  258.6/257.5
    2       259.8     258.8 |  258.3/257.3
    3       266.1     265.4 |  264.9/264.4
    4       274.2     273.7 |  275.8/275.7
    5       281.6     281.5 |  284.3/284.6
    6       287.9     288.4 |  291.1/291.5
    7       290.7     291.5 |  295.2/295.6
    8       288.9     289.5 |  294.7/295.0
    9       283.0     283.1 |  288.3/288.6
   10       275.0     274.5 |  277.9/277.9
   11       265.8     264.9 |  268.6/268.1
   12       259.7     258.6 |  262.2/261.4
Jul-Jan  ERA5 t2m 33.1 skt 35.0 |   36.6/38.1 
(model columns: lowest-level T / skin)
```
