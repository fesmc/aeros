# aeros vs ERA5, real-Earth T21L6 (2026-09-09 late: CLIMBER-X RH-profile reference in the SBM)

aeros 3 yr (2 yr spin-up + 1 yr monthly), dt = 3600 (run `FINAL4` = `P6`, defaults only). On top of e56e375: the SBM cloud-layer reference humidity is CLIMBER-X's `rh_prof` (`conv_rh_ref 0.80` = ram to `conv_rh_zpbl` 500 m, exp(−(z−z_pbl)/3000 m) to `conv_rh_ztop` 3000 m, constant 0.35 above); `cld_rhc_ice 0.40`, `cld_lwc 0.022e-3`.

Year-3 global means: TOA net −0.08 W/m², OLR 243.9 (ERA5 242), clear-sky OLR 262.5 (264; before 258), CRE SW −45.0 / LW +18.6 (−46 / +22), cover 0.49 (0.63), P = E = 3.02 mm/day (2.9), LH 87.4 (84), column water 21.6 mm (24.5; before 24.5), max [u] 31.4 m/s (33.7).

Tropics 15S–15N RH at σ 0.365/0.563/0.728/0.860/0.959: 40/34/44/71/86 % (ERA5 37/45/56/70/84; before 67/55/57/61/85). The mid-troposphere is now ~10 % dry (the column water deficit). Tropical T bias at σ 0.365 unchanged, +2.7 K: the budget there is convection +2.45 K/day vs radiation −2.2, the surface-anchored moist-adiabat reference (no entrainment); not a humidity property.

Lid T minus the ERA5 10–257 hPa layer mean: 60–90S −1.2, 30–60S +2.8, tropics +2.1, 30–60N +2.7, 60–90N +0.8.

Zonal T bias vs ERA5 (annual, K; 60–75S / 30–45S / 15S–15N / 30–45N / 45–60N / 60–75N):

```
sigma        60-75S     30-45S    15S-15N     30-45N     45-60N     60-75N
0.134         -0.3        7.3       12.0        8.5        1.3       -0.1
0.365         -0.9       -2.0        2.7        2.5       -0.5       -1.9
0.563         -2.8       -4.4        0.9       -1.4       -2.1       -2.7
0.728         -2.0       -4.0       -0.8       -2.5       -1.6       -2.1
0.860         -0.5       -1.6       -1.7       -2.4       -0.5       -1.0
0.959         -0.5       -0.3       -0.0       -0.3        1.2        0.7
```

The 30–45S mid-troposphere cold (−4.4/−4.0) is the pre-existing deficit, 0.3 K worse; 45–75N 0.56–0.73 −1.6 to −2.7.

NH 40–70N land: January LW down 220 (209), 2 m amplitude 35.0 (33.6); July SW net 238 (196; the thinner liquid cloud), SH 82 (40), LH 62 (73), cover 0.37 (0.57), 2 m July 294.7 (290.7). Unchanged in character from the aerosol note.

```
month      ERA5 t2m   ERA5 skt |    FINAL4.nc
    1       257.6     256.5 |  259.7/258.7
    2       259.8     258.8 |  261.3/260.4
    3       266.1     265.4 |  266.9/266.5
    4       274.2     273.7 |  276.6/276.6
    5       281.6     281.5 |  284.9/285.2
    6       287.9     288.4 |  291.2/291.6
    7       290.7     291.5 |  294.7/295.1
    8       288.9     289.5 |  293.8/294.2
    9       283.0     283.1 |  287.7/288.0
   10       275.0     274.5 |  278.0/278.0
   11       265.8     264.9 |  269.8/269.4
   12       259.7     258.6 |  263.2/262.4
Jul-Jan  ERA5 t2m 33.1 skt 35.0 |   35.0/36.5 
(model columns: lowest-level T / skin)
```
