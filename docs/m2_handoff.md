# M2 handoff — where to pick up

## ►► NEXT SESSION (written 2026-09-06): state, performance, recipe, and what to do

### Where the model is

All on `main` (last commits 59657fa, 90e9408, 63c7c4e), 26/26 tests, built with
`make all openmp=1 && make rce-long openmp=1`. The production configuration is the
real-Earth T21 run: ERA5 orography, ERA5 land-sea mask with the bucket land + snow albedo,
prescribed 2D seasonal ERA5 SST, seasonal insolation, ecCKD radiation with diagnostic clouds,
Frierson SBM convection, SpeedyWeather-style condensation. **Namelist:**
`mwm/D_realearth/aeros_realearth_T21L6.nml` (only physical keys; every numerical/core knob
comes from `input/rce_defaults.nml`, which now IS the production recipe).

Since the 2026-07-30 handoff the integration and the coupling were rebuilt (details in the
log sections below). What is different and load-bearing:

| item | now | was |
|---|---|---|
| time scheme | semi-implicit Lorenz 4-cycle (`time_scheme="ncycle"`), no time filter | RAW-filtered leapfrog |
| physics coupling | increments on the one time chain; surface fluxes + vdiff implicit on the stepped state | forward-split onto n+1 (**half rate**, a 2Δt computational mode) |
| tracer transport | FV tracer mass divided by the spectral Δp(n+1); water conserved to ~1% of P | divided by the FV Δp (**created ~30% of the precipitation**) |
| surface exchange | bulk-Richardson log-law coefficient per point (SW/Frierson 2006); SW surface wind (0.95 slowdown + 5 m/s gust) | fixed 1.5e-3, 1 m/s wind floor |
| speed-limit drag | quadratic above 70 m/s (SW's default is 80) | none |
| dynamical core | ∇⁸/4 h, divergence ∇⁸→∇⁴ above σ=0.2 on 1 h, α=1, no sponge — as defaults | opt-in |
| vertical grid | **L6** default (stretched hybrid; L8/L12 namelist choices) | L12 |
| dt | **3600 s** at T21 | 900 s |
| monthly climatology | grid diagnostics refreshed at each sample; + precip2d, z500, vlow, evap, shf | stale (state of the last report) |

**Every number in this file older than 2026-09-05 (RCE tuning, land amplitudes, seasonal
validation, the SpeedyWeather omega comparison) was produced with half-rate physics, a
water-creating transport and stale monthly sampling. Treat them as history, not baselines.**

### Performance (Apple M5, 10 threads, real-Earth T21, dt = 3600, 2-yr runs with monthly sampling)

| levels | SYPD | comment |
|---|---|---|
| L4 | 3650 | stable; harness metrics as good as L8 |
| L5 | 3250 | |
| L6 | 2810 | **default since 2026-09-06**; the 3-yr climatology below |
| L8 | 2190 | default 2026-09-05/06; `docs/refs/realearth_metrics_2026-09-05.md` |
| L12 | 1870 | previous default |
| SpeedyWeather T21L8 default | 1040 (1360 without the diagnostic callback) | leapfrog, dt 3600 |

Cost is close to linear in levels with a fixed overhead (transforms and physics are both
per level); aeros is ~2.5× faster than SpeedyWeather per step. L20 blew up at day 354 at
dt=3600 (not investigated). SpeedyWeather's 8 levels are SPEEDY's (Molteni 2003).

### Climatology, L6, 3 yr (aeros / SpeedyWeather / ERA5) — `docs/refs/realearth_metrics_2026-09-06_L6_clouds.md`

(after the mixed-phase saturation and the coupled cloud calibration, evening 2026-09-06)
precip = evap **2.7** / 3.3 / 2.9 mm/day · TOA +1.3 W/m² (yr-3 mean) · absorbed SW 231 (ERA5 243),
OLR 230 (242), CRE SW −57 / LW +31 (ERA5 −46 / +22), cloud cover 0.52 (0.63) · lowest-level T
283.3 / 278 / 287.5 K · max [u] 57 / 63 / 34 m/s · DJF NH Z500* RMS **93** / 68 / 76 m · NH-land
t2m Jul−Jan **36.9** / 29 / 33.6 K · 2830 SYPD. Figures in `docs/figures/realearth_*.png`.
Before the cloud retune (mixed phase only): P 2.2, Z500* 68.5, t2m amplitude 29.8 — the
brighter planet (LH 79 vs 59 W/m²) over-drives the stationary wave and the land seasonal
cycle; the land calibration was already flagged void. Earlier L6 (liquid, old clouds):
`docs/refs/realearth_metrics_2026-09-06_L6.md` (P 2.2, Z500* 85, t2m 31.9); L8:
`..._2026-09-05.md`. aeros' zonal temperature structure is closer to ERA5 than
SpeedyWeather's (SW is cold aloft with a 235 K lid).

**Known deficits, in the order to attack them:**
1. **Saturated upper troposphere / lid** (RH ~100% above 200 hPa; ERA5 30–50%, SW dry at its
   top). Establish the mechanism first: the FV transport carries q into the cold lid, and the
   condensation there (rh_crit·q_sat with a tiny q_sat) may be removing next to nothing
   per step, or rain re-evaporating aloft, or the σ-coordinate vertical flux at k=1. Compare
   with SW's `rh` at its top level in the harness output; the per-operator water tally
   (`l_diag`) can be extended per level.
2. **Spurious polar precipitation** (Antarctic ~4 mm/day in DJF, Arctic in JJA; ERA5 < 0.5).
   The prescribed SST over the polar ocean is the 271 K freeze floor with no sea ice, so cold
   air over "open water" evaporates and condenses. Sea ice (`l_seaice`, needs `ocean_mode=1`)
   or an ice-concentration mask on the ocean flux (ERA5 `siconc` is on disk) is the physical
   fix; the lowest-level cold trap may add to it.
3. **Evaporation ~30% low.** Near-surface wind speed 5.2 m/s vs ERA5 ~7 m/s at 10 m: the
   zonal means match ERA5, the transient eddies are weak (T21), and the lowest layer sits at
   RH 80%. Not a surface-scheme problem any more; a T31 harness run (both models) would show
   what resolution buys.
4. Upper jets ~1.6× ERA5 in both models.

### Recipe

```bash
make all openmp=1 && make rce-long openmp=1
# the 26 tests
for t in libaeros/bin/test_*.x; do OMP_NUM_THREADS=5 $t > logs/suite_$(basename $t .x).log 2>&1 || echo FAIL $t; done
# aeros 3-yr real-Earth climatology (2-yr spin-up + 1 yr monthly), ~2 min
OMP_NUM_THREADS=10 libaeros/bin/rce_long.x mwm/D_realearth/aeros_realearth_T21L6.nml
# SpeedyWeather yardstick (~3.5 min), then the comparison (figures + metrics table)
julia --project=mwm/D_realearth -t 10 mwm/D_realearth/sw_realearth.jl 730 1 21 8
julia mwm/D_realearth/compare_realearth.jl mwm/D_realearth/output/aeros_realearth_T21L6.nc
```

The `rce_long` report prints a `TOA:` and a `WATER:` line every `print_every` steps; the
water-budget residual must stay ~0.01 mm/day (it was +0.6–0.9 before the transport fix).
With `l_diag = .true.` it adds a per-operator water tally and per-term heating tables (the
`surf` column is now zero: the surface fluxes live inside `vdiff`). `nlev`, `dt`,
`drag_vmax`, `gust_speed`, `z0_ocean/z0_land`, `c_h/c_e/c_d` (≤ 0 = log-law), `cond_time_scale`,
`conv_tau`, `time_scheme`, `ncycle_variant` are `&rce` keys with defaults in
`input/rce_defaults.nml`. Gotchas: `make all` does not build `rce_long` (separate target);
`rce_long` reports no wall time (use `/usr/bin/time -p`); the monthly file is written only
at a clean end (not on a NaN).

### First actions next session

1. ~~Make **L6** the default~~ — done 2026-09-06 (library + driver defaults, harness namelist
   `aeros_realearth_T21L6.nml`, climatology above).
2. Deficit 1 (upper-troposphere moisture), then deficit 2 (polar precipitation).
3. ~~Re-do the land t2m calibration~~ — done 2026-09-07 (seasonal albedo, c_soil 5e6); the residual is the winter atmospheric cold bias.

---

## Log 2026-09-05/06 — the time integration was the speed AND stability problem (detail)

**Two structural defects in the physics coupling, both gone (branch `fix/physics-coupling`).**

1. **Forward-split physics landed at HALF its rate.** `apply_phys_heating` added the increment
   Δt·P (diagnosed from X^n) to X^(n+1) only; each leapfrog chain advances 2Δt per step, so
   every diabatic term (surface, convection, condensation, radiation) heated at half the
   diagnosed rate while the single-level gridpoint humidity dried at the full rate (latent
   heat leak ~L·P/2, sensible heat into the air half what the soil/ocean lost). Verified with a
   radiation-only run: level-1 cooling 2.4 K/day against a diagnosed ~5; the in-solve rate
   coupling gives 4.7. **Every M2 physics number / tuning before this date sits on top of
   this** (rh_crit, albedo, cloud calibration, the +15 W/m² TOA residual, the ERA5 and
   SpeedyWeather comparisons). Held–Suarez is unaffected (in-solve `dtdt`).
2. **The "dt=900 CFL limit" was a 2Δt computational mode of that split**, not CFL: at dt=1800
   the lowest-layer global-mean T alternates by 7–18 K between consecutive time levels and the
   physics heating flips sign every step (an eastward-moving tropical hot spot). Reproduced in
   4 s on the production land+topo+seasonal config (NaN day 22.5); killed by removing
   reevaporation / convection / dry adjustment, by eps_filter ≥ 0.15, or dt ≤ 1350. Root:
   stiff operators (exact saturation adjustment, reevap-to-threshold) diagnosed on one
   leapfrog chain and applied to the other. The earlier "model-top / spectral-pole / CFL"
   framings in this file and in the memory notes were this mode.

**What landed (26/26 tests):**
- **Physics coupling is in-solve only**: `wrk%dt_phys/dt` is added to `dtdt` before the
  semi-implicit solve (SPEEDY/SpeedyWeather coupling); the forward-split path and
  `couple_diabatic` are removed. Operators softened to SpeedyWeather's values:
  condensation + reevaporation relax over `cond_time_scale = 3` steps (new knob, module +
  `&rce` + defaults), convection/dry-adjust `conv_tau = 4 h` (was 2 h).
- **Two-time-level dynamics: the semi-implicit Lorenz N-cycle** (Hotta et al. 2016),
  `time_scheme = "ncycle"` (N=4, variant A default; ABBA is ~4th order), now the **default**
  in the library and the driver; `"leapfrog"` kept for comparison (the four leapfrog-specific
  tests pin it). One time chain for T, q and the physics → no cross-chain physics, no
  computational mode, no Robert filter (`eps_filter`/`raw_alpha` are leapfrog-only). The
  accumulator G is in the restart (`gacc_*`, `ncycle_present`); `test_ncycle` (26th test):
  rest state exact, stability triple, convergence (A 2nd order, ABBA ~4th), restart bit-exact
  mid-cycle. Module header of `aeros_timestep.f90` has the derivation.
- **SW-faithful core is the driver default** (`input/rce_defaults.nml`): ∇⁸/4 h,
  `tau_diff_div = 1 h`, taper to ∇⁴ above σ=0.2, `si_alpha = 1`, no sponge; those keys plus
  the sponge strengths are now optional in `&rce`. A namelist with only the physical keys
  reproduces the explicit production config bit-for-bit.

**Result.** The production config (land + ERA5 topography + seasonal SST, T21L12) runs
NaN-free for a year at **dt = 3600 s** (dt=1800/2700/3600 all stable; 5400 dies at day 3,
untested whether CFL or physics). Throughput on 10 cores (Apple M5): **~1850 SYPD** vs 834
before, and vs SpeedyWeather 0.21.1 default T21L8 (slab ocean, sea ice, land, orography):
1357 SYPD at dt=3600. aeros is ~2.5× faster per step than SW; per model-year now 1.4× faster.

**Update 2026-09-05 (later, commits 6bf625d / 9eb5868).** (i) The N-cycle now applies the
physics increment to the stepped state instead of through the weighted accumulator G (the
Lorenz weights extrapolated the relaxation terms: −5.9..+0.1 K/day sub-cycle swings and a
dt-dependent OLR; gone: OLR 131.6/131.3/130.2 at dt 900/1800/3600). (ii) The real-Earth 1-1
harness exists: `mwm/D_realearth/` (`sw_realearth.jl`, `compare_realearth.jl`,
`aeros_realearth_T21L12_dt2700.nml`; figures `docs/figures/realearth_*.png`, table
`docs/refs/realearth_metrics_2026-09-05.md`). (iii) **`monthly_out` was sampling stale grid
diagnostics** (refreshed only at `print_every`) — every earlier seasonal/land number (land
Task 1/2 amplitudes, `rce_seasonal_*`) was taken from months that were snapshots of the last
report; fixed. (iv) **dt = 3600 dies in the second winter** (lid jet 91 m/s at level 1 →
past the N-cycle advective limit ~0.7 of leapfrog's); **dt = 2700 is the working production
step** (3 yr NaN-free, ~1600 SYPD with monthly sampling). (v) Ozone is NOT the issue: the
real-Earth lid is 214 K vs ERA5 207 K at 70 hPa (the 188 K was the aquaplanet's).

First real-Earth numbers (T21, 2 yr spin-up + 1 yr): aeros / SpeedyWeather / ERA5 —
precip 2.1 / 3.3 / 2.9 mm/day; lid T 214 / 235 / 207 K; lowest-level T 288 / 278 / 288 K;
max [u] 63 / 63 / 34 m/s; DJF NH Z500* RMS 47 / 68 / 76 m; NH-land t2m |Jul−Jan| 29 / 29 /
34 K; 1606 / 1016 SYPD. aeros' zonal T is closer to ERA5 than SW's (SW is cold aloft with a
warm lid); aeros' deficits are **too little precipitation/evaporation (LH ~35 W/m²)**, a
**saturated upper troposphere/lid (RH ~100% above 200 hPa)**, **spurious polar precipitation**,
and weak stationary waves; both models over-strong upper jets.

**Update 2026-09-05 (evening, commits 9236f67 / <transport>).** (vi) **SpeedyWeather's
speed-limit drag** (quadratic above 80 m/s) is now in aeros; it is what lets SW hold dt = 1 h,
and with it aeros runs the real-Earth config 3 yr at **dt = 3600** (1870 SYPD), climate
unchanged vs 2700. (vii) **The FV tracer transport was creating water**: 0.5–0.9 mm/day, a
third of the precipitation — q was divided by the FV-advanced air mass and the next step
re-diagnosed layer masses from the spectral p_s. Fixed by dividing by the spectral layer
mass at n+1 (constancy of a uniform tracer becomes approximate; cf clamped at 1). The
`rce_long` report now prints a global water budget (W, P, E, residual) and, with `l_diag`, a
per-operator water tally — use them. (viii) Surface fluxes use SpeedyWeather's surface wind
(0.95 slowdown + 5 m/s gust in quadrature) instead of a 1 m/s floor.

Real-Earth metrics after (vi)–(viii): P = E = **1.85** mm/day (SW 3.3, ERA5 2.9), W 30 mm
(ERA5 ~25), lowest-level T 287.5 (ERA5 287.5), max [u] 60 (SW 63, ERA5 34), **DJF NH Z500*
RMS 67 m** (SW 68, ERA5 76; was 40–50 with the leak), NH-land t2m amplitude 30 K (ERA5 34).
Remaining first-order deficit: **evaporation ~35% low** — aeros' lowest-level wind is 3 m/s
vs SW's 5.6 (its lowest level is ~500 m vs aeros' ~120 m), so V0 is 5.8 vs 7.7 m/s; the
exchange coefficient (1.5e-3 fixed vs SW's Richardson ~1e-3) is not the gap.

**Update 2026-09-05 (night, commit 6f55928).** Surface exchange coefficient is now
SpeedyWeather's bulk-Richardson log-law (Frierson 2006), per point, one coefficient for
heat/moisture/momentum; the BL diffusion was already the same scheme. Findings on the
**evaporation deficit** (aeros 1.7 vs SW 3.3 vs ERA5 2.9 mm/day): zonal-mean near-surface
winds match ERA5 and SW (|[u]| 2–3 m/s), the wind *speed* is low (4.6 vs SW 5.6 vs ERA5 ~7)
because aeros' eddies are weaker; and the lowest-level RH is 82% (SW 75%) — aeros' lowest
level sits at ~160 m (L12 stretched grid) vs SW's ~500 m, and a thin, moist lowest layer
gives a small saturation deficit at the same coefficient. With the fixed 1.5e-3 coefficient
E was 1.85; with the log-law (~7e-4 at 160 m over sea) 1.70. Momentum drag is weaker too, so
the winter lid jet reaches ~95 m/s and **dt=3600 dies at day 422; the harness runs dt=2700**
(drag_vmax=70 survives 2 yr at 3600 with 94 m/s peaks — an option). **Vertical levels:** the
same real-Earth run at **L8** is +15% faster (2145 SYPD), NaN-free, and its climate metrics
are as good or better (precip 2.1, Z* 69 m, lid 223 K); **L20 blows up at day 354** at
dt=3600 (untested at smaller dt). The L12 grid's motivation (design §4.1: polar boundary
layer / inversions for the ice-sheet SMB) has not yet been tested by any metric.

**Update 2026-09-06 (commits 59657fa, 90e9408).** (ix) **L8 is the default vertical grid**
(library + driver; L12 stays a namelist choice, Held–Suarez keeps L20). Level sweep on the
real-Earth run, 2 yr at dt=3600 with the speed limit at 70 m/s: L4 3650 / L5 3250 / L6 2860 /
L8 2190 / L12 1870 SYPD, all NaN-free, and the harness metrics are as good or better with
fewer levels (precip 2.3→1.9 mm/day from L4 to L8, Z500* 76–89 m, t2m amplitude 29–31 K).
The design-doc case for L16–20 (polar boundary layer for ice-sheet SMB) has no metric yet;
L20 blew up at day 354 at dt=3600 (not investigated). SpeedyWeather's 8 is SPEEDY's
(Molteni 2003) count — its physics was built on it. (x) **Speed limit 70 m/s at dt=3600** is
the production default (the log-law surface drag lets the winter lid jet reach ~95 m/s
otherwise). (xi) **Surface fluxes are the implicit lower boundary of the vdiff solve**
(sensible heat in dry static energy, evaporation in q, like the momentum stress): the
lowest-layer response is independent of that layer's thickness and of dt; `aeros_surface`
only diagnoses the conductances; ocean/land reservoirs step after the solve; `test_surface`
rewritten. **L8 real-Earth, 3 yr:** P = E = 2.0 mm/day (SW 3.3, ERA5 2.9), TOA −0.2 W/m²,
Z500* 77 m (ERA5 76), t2m amplitude 32.5 K (ERA5 33.6), lowest-level T 284.7 (ERA5 287.5),
2190 SYPD (SW 1040). The evaporation gap is now the wind speed (5.2 m/s lowest level vs
ERA5 ~7 at 10 m; the zonal means match, the eddies are weak) and the moist lowest layer
(RH 80%).

**Update 2026-09-06 (L6 default).** (xii) **L6 is the default vertical grid** (library + driver;
L8/L12 stay namelist choices, Held–Suarez keeps L20). L6 real-Earth 3 yr, dt=3600, NaN-free,
water residual −0.007 mm/day: P = E = 2.2 mm/day (L8 2.0; SW 3.3; ERA5 2.9), TOA −3.0/−3.7/+3.2
W/m² over the three years, lowest-level T 283.5 K (L8 284.7; ERA5 287.5), lid 226 K (L8 222),
max [u] 52 m/s (L8 56), DJF NH Z500* RMS 85 m (L8 77; ERA5 76), NH-land t2m amplitude 31.9 K
(L8 32.5; ERA5 33.6), **2810 SYPD** (L8 2190; SW 1040). The metrics move within the level-sweep
scatter; the stationary wave is the one that worsens (85 vs 76 m).

**Update 2026-09-06 (deficits 1 and 2 diagnosed, no fix yet).** Tools: the L6 harness run with
`l_diag`, the per-term heating tables, `output/rce_rh.nc`, ERA5 pressure-level q/T/r and TOA
fluxes (`tisr/tsr/ttr/tsrc/ttrc/tcc`), 2-yr sensitivity runs (ozone off, clear-sky radiation,
L12). Findings:
- **The lid is not a moisture problem, it is a temperature problem.** aeros' top layer (L6:
  0–250 hPa) is 231 K in the tropics against ERA5's mass-weighted 213 K, and 199–206 K at the
  poles against 213–220 K: a 40 K equator-pole contrast at the lid where ERA5's is reversed.
  That warm layer holds 0.62 g/kg at the 0.95·q_sat(liquid) cap (ERA5 0.04; ice saturation would
  cap it at 0.44, RH_ice is 134–158%). SpeedyWeather's 0–250 hPa is just as warm (232 K) and
  holds 0.25 g/kg; its "dry lid" is an RH artefact of a warm low-pressure top layer.
- **It is not the thick top layer.** L12 with the current code: 245 K at 198 hPa (ERA5 220),
  257 K at 313 hPa (ERA5 240), RH 66–93% there. The tropical warm bias grows from ~0 at the
  surface to +20–25 K at 200 hPa at every level count.
- **Not ozone** (`l_o3 = .false.`, now an `&rce` key: lid 230.3 K, 0.58 g/kg, unchanged),
  **not cloud LW trapping alone** (clear-sky run: lid warmer, 235 K; TOA +42 W/m²).
- **Heat budget of the tropical upper troposphere** (per-term tables, L6 and L12): convection
  is zero above ~450 hPa (the surface parcel's h_b 330–340 kJ/kg never reaches h* 375–410
  kJ/kg of the warm layers); grid-scale condensation heats 200–300 hPa at +0.3–0.9 K/day
  (≈0.5 mm/day of tropical rain forms in the top layer); radiation cools it at only −0.2 to
  −0.5 K/day (Earth: −1 to −1.5). A self-consistent warm state: warm aloft → h* too high →
  no deep convection → the resolved ascent does the convecting with 100% precipitation
  efficiency → warm, saturated, cirrus-capped.
- **The planet runs dim.** Global TOA (aeros L6 / ERA5): OLR 189 / 242, absorbed SW 185 / 243
  (planetary albedo 0.47 / 0.29), clear-sky OLR ~256 / 264, cloud LW effect ~+67 / +22,
  cloud SW effect ~−113 / −46, total cloud cover 0.75 / 0.63. The diagnostic cloud scheme
  (tuned on ERA5 columns) driven by the model's own saturated upper troposphere makes a
  permanent cirrus deck (cf 0.4–0.8 at the lid). The net TOA is balanced by construction,
  but at 185 W/m² throughput the atmosphere's radiative cooling (~0.4 K/day column mean) and
  hence P = E (2.2 vs 2.9 mm/day) and the tropical LH (51 vs ~110 W/m²) are all ~25% low.
  Deficits 1 and 3 (evaporation) are one loop.
- **Deficit 2 (polar precipitation) is not there at L6/L8.** Cap means (aeros / ERA5, mm/day):
  DJF 60–90S 0.91 / 0.75, JJA 60–90N 1.26 / 1.41, DJF 60–90N 0.43 / 0.88, JJA 60–90S 0.49 /
  0.90 — the winter poles are too dry, not too wet; the only excess is ~0.5 mm/day over the
  summer pole itself. The 4 mm/day in the earlier note was an L12/old-physics feature. The
  real precipitation gaps are the ITCZ (4.3 vs 6.4 annual peak) and the storm tracks.

**Update 2026-09-06 (evening): the warm tropical upper troposphere survives every switch.**
New diagnostics: the report prints clear-sky TOA fluxes, CRE and cloud cover (`TOA clear:` line);
`cld_rhc_sfc/top`, `cld_lwc/iwc` are `&rce` keys. Every run below is the L6 real-Earth config, 2 yr,
one change at a time; "lid" = tropical (|lat|<15) top-layer T / q (ERA5 layer mean 213 K / 0.04 g/kg):

| run | lid T [K] | lid q [g/kg] | OLR | CRE SW / LW | P [mm/day] |
|---|---|---|---|---|---|
| default | 230.9 | 0.62 | 189 | −113 / +67 | 2.18 |
| ozone off | 230.3 | 0.58 | 186 | | 2.27 |
| clear-sky radiation | 235.0 | 0.91 | 256 | | 2.63 |
| no high cloud (`cld_rhc_top=1`) | 229.9 | 0.55 | 222 | −59 / +35 | 2.64 |
| convection off | 229.2 | 0.52 | 178 | −126 / +76 | 2.07 |
| van Leer vertical q transport | 231.0 | 0.60 | 199 | −100 / +60 | 2.34 |
| dt = 900 | 231.8 | 0.66 | 195 | −106 / +60 | 2.31 |
| SW water absorption by mass (scratch) | 229.0 | 0.51 | 181 | | 2.11 |
| topography off | 232.5 | 0.73 | 186 | | 2.24 |
| aquaplanet, real SST + seasons | 233.1 | 0.79 | 180 | | 2.47 |
| uniform insolation, real SST/land | 231.3 | 0.65 | 186 | | 2.19 |
| idealized SST, seasonal insolation | 215.8 | 0.11 | 151 (TOA +42) | | 1.08 |
| RCE (uniform insolation, idealized SST) | 214.2 | 0.09 | 146 (TOA +23) | | 1.07 |
| L12 default (198 hPa level) | 245 at 198 hPa | | 166 | | |

- The column physics can make the right tropical profile: the two idealized-SST runs sit on the
  moist adiabat with a 214–216 K lid, deep convection to the top, and 0.1 g/kg aloft — but with a
  half-strength hydrological cycle (LH 31 W/m²) and a dead equatorial column in one of them.
- With the real SST the tropical column is warm and stable above ~450 hPa at any level count;
  the convective band (Frierson) never reaches it (h_b 330–340 kJ/kg vs h* 375–410); grid-scale
  condensation does the deep moist processing there (+0.3–0.9 K/day) and radiation cannot remove
  it. Not caused by clouds, ozone, the SW distribution, the convection scheme, transport, dt, land,
  topography or the insolation pattern individually. It is a circulation-regime property: the
  boundary layer is too dry (12.6 vs 14.8 g/kg at 950 hPa, RH 75% vs 85%), h_b is ~5–10 kJ/kg low,
  the mean meridional circulation is 3–4× too weak (ITCZ ascent −7 hPa/day at 500 hPa vs ERA5
  ~−25), low-level winds 5.2 vs ~7 m/s, evaporation ~35% low — one loop.
- Offline column probe (scratch `probe_lid.f90`): the SESAM/ecCKD SW deposits the near-IR
  water absorption by the drop of a 2-exponential transmission, so the first ~1.5 kg/m² of water
  from the top absorb ~40 W/m² (+1.5 K/day) — 2–3× a reference profile; the coupled test says it
  is secondary (−2 K on the lid). Worth a proper near-IR k-distribution later.
- SW-cloud fix (`cld_rhc_top=1`) brings the CRE to −59/+35 (ERA5 −46/+22) and P to 2.64 but leaves
  TOA at +17 W/m² (albedo was tuned against the old clouds) — a retune is a 4-knob / 3-target job.
- Polar-cap precipitation is NOT a deficit at L6/L8 (see the morning entry).

**Update 2026-09-06 (night): mixed-phase saturation and the coupled cloud calibration.**
- **Mixed-phase saturation** (577b9ca): `aeros_qsat` blends liquid Tetens (≥ 0 °C) and ice Tetens
  (≤ −23 °C) with the IFS quadratic weight; `aeros_latent(T)` returns the matching L_v→L_s and is
  used by condensation, convection (per-layer L held through a step, so the Betts–Miller /
  Frierson / Manabe energy closures stay exact; the Frierson shallow branch now shifts T_ref by
  the latent heat of the redistributed water), the surface latent heat flux (L_s over ice/snow)
  and the diagnostics. RH in the monthly file is now the ERA5 convention. Effect alone: lid water
  −15%, P 2.2, Z500* 68.5 m, t2m amplitude 29.8 K.
- **The energy-budget report is a time mean** (aed113e) — the snapshot at the Dec-31 report
  date had SWin 352 (perihelion) and a ±10 W/m² noise on the net.
- **Cloud calibration** (this commit): the clear-sky budget was already ERA5's (absorbed 288 / 289,
  OLR 262 / 264), so the whole dim-planet excess was the cloud effect. The LW CRE is set by the
  lid layer's cloud FRACTION, not its ice content: the 0–250 hPa layer is ~13 km deep, so any
  ice makes it black at 231 K under 65% cover (iwc/4 → LW +61). Sweep R1–R12 (`logs/cld_R*.log`):
  `cld_rhc_top = 1.07` (threshold ≈ 1 in the lid, 0.87 at 365 hPa, unchanged below), `cld_lwc
  0.020e-3`, `cld_iwc 0.014e-3` → CRE −56/+31, cover 0.52, net +1.6 with no albedo change,
  P 2.70, LH 79 W/m². The LW floor of ~+30 with no lid cloud is the mid-level layers' geometry.
- **Cost:** Z500* 93 m (ERA5 76) and NH-land t2m amplitude 36.9 K (33.6): the stronger
  hydrological cycle over-drives both. Next: the land calibration (already void), and the
  stationary wave (orographic + thermal forcing at T21).
- The tropical lid is unchanged (225 K top-level mean; the mechanism entry above stands).

**Update 2026-09-07: land recalibration (through the harness).** Findings, then what changed:
- The harness metric compares the model's LOWEST LEVEL (~350 m at L6) with ERA5 2 m. Under the
  winter inversion that hides half the bias: at 60–70N land the model's DJF lowest level is
  241 K vs ERA5 2 m 249 K and ERA5 925 hPa 255 K. **The DJF cold is the whole lower troposphere**
  (zonal mean at 64–69N: 700 hPa −7 K, 850 −10 K, 925 −13 K; JJA within 1 K), land and ocean
  alike — an atmospheric (winter heat-transport / boundary) problem, not the land scheme.
- **Seasonal land albedo** (`land_albedo_seasonal = .true.`, `aeros_land_set_albedo`): the 12
  monthly ERA5 `fal` slices interpolated daily (mid-month centres, as the seasonal SST). The
  annual-mean map is snow-contaminated (0.36 at 65N, vs 0.15 in July) and left JJA high-latitude
  land 6 K cold; with the monthly field JJA matches ERA5 within 1–2 K in every NH band, and the
  diagnostic snow bump (`snow_albedo`) is redundant (off in the harness). Bug fixed on the way:
  the per-step albedo-map refresh in `aeros_timestep` was gated on the snow bump, so without it a
  seasonal albedo would have stuck at its January value.
- **`c_soil` IS a lever now** (the 2026-07 "not the lever" note was on stale sampling): 2e6 → 1e7
  cuts the NH-land amplitude by ~5 K (DJF +2, JJA −3.5). Chosen **5e6** (physical: 2.5e6 J m⁻³ K⁻¹
  over the ~2 m annual penetration depth), which keeps JJA right rather than buying winter with a
  summer cold bias. Single-month land means over 60–90 grid points vary by ±2–3 K between 2-yr
  runs; calibrate on DJF/JJA means.
- Harness 3 yr (`docs/refs/realearth_metrics_2026-09-06_L6_land.md`): P 2.70, TOA +2.1, Z500* 80 m
  (ERA5 76; was 93 with the annual-mean albedo), lowest-level |Jul−Jan| NH-land 38.6 K (the
  single-month lowest-level metric; DJF/JJA band table in the ref file), 2720 SYPD.
- **Next:** the winter high-latitude lower-troposphere cold bias (7–13 K) — candidates: winter
  low cloud / LW down over cold land and sea ice (`cld_rhc_sfc`, mixed-phase RH), the polar ocean
  freeze floor without sea ice, eddy heat transport at T21. A 2 m diagnostic from the surface layer
  would make the harness metric honest.

**Update 2026-09-07: the winter cold bias was the humidity transport advecting AGAINST the
meridional wind (d2f96fc).** Chain of evidence, in the order it was found:
- New monthly diagnostics (`skin`, `lwd`, `olr`, `olr_cs`, `cover`, per-level `cf`; e272fbb) showed
  the winter cap atmosphere loses LESS heat locally than ERA5 yet is colder: the implied northward
  energy transport was half of ERA5 at every latitude 30–75N (cap >60N: 52 vs 99 W/m² convergence,
  latent 13 vs 32). Not the surface (prescribed SST = ERA5 skt incl. sea ice; the old `sst` output
  was the land-filled ocean field, not the skin), not radiation. Unchanged at L8, L12 (dt 1800 —
  both NaN at the top at dt 3600 with the current physics), T31, tau_diff 12 h, `vert_vanleer`,
  `cond_rh_crit` 1.0.
- Eddy-flux diagnostics (`vt`, `vq`, `eke` + `*_stat`, in rce_long and the SW harness; 3b25ed2):
  aeros' transient eddies carried heat poleward (16 K m/s at 850 hPa, twice SW) but MOISTURE
  EQUATORWARD at 500–850 hPa; SW's [v'q'] is positive at every level.
- The tracer mass-rescaling term dp_FV/dp_spec − 1 (54b8f7d) had a persistent pattern of 1.6–1.9 %
  per step in mid-latitudes and 13 % at the polar row: the FV meridional mass flux ran opposite
  to the spectral one. Root cause: `sweep_merid` used vf = +½(v_j + v_j+1) with rows north→south,
  i.e. the along-row velocity of a NORTHWARD grid v (aeros_vordiv: v = −V_θ) has the wrong sign.
  A tracer band under v = +5 m/s moved south. Every transport test was direction-blind;
  `test_direction` now pins north and east.
- After the fix (same calibration): r rms 2.5 → 0.35 %/step, eddy latent transport 0 → 0.9 PW at
  58N, 69N DJF lowest-level bias −12.3 → −1.4 K, NH-land amplitude 38.6 → 32.0 K, jet 56 → 50 m/s;
  but the free troposphere is moist everywhere (60–90 % RH at 360–730 hPa vs ERA5 30–65 %) and
  the planet was +11 W/m² (OLR 216). The reversed transport had been compensated by the cloud,
  condensation and land calibrations: all pre-fix moisture/cloud numbers in this file are void,
  and the "subtropical RH 94 %", "saturated lid" and "water-creating transport" histories should
  be re-read with this in mind.
- **Recalibration** (`mwm/D_realearth/cal/`, 15 runs; `docs/refs/realearth_metrics_2026-09-07_L6_vfix.md`):
  `cond_rh_crit` and the new `conv_rh_ref` key move the clear-sky OLR by < 3 W/m²; it is 250 vs
  ERA5 264 in every latitude band because the troposphere is 3–6 K cold at 500–850 hPa
  (prescribed SST, lowest level within 1–2 K). Balance therefore needs CRE ≈ −46/+10: defaults
  now `vert_vanleer = .true.`, `cld_rhc_top = 1.40` (little lid cloud), `cld_lwc = 0.035e-3`
  → TOA +1.6, absorbed 241, OLR 240, P 3.0, cap transport 2.18 PW (ERA5 2.19), Z500* 82,
  amplitude 34.8 (33.6). Cost: with thin high cloud the free-troposphere cold bias is 4–10 K at
  500–850 hPa at all latitudes and seasons (69N DJF −8/−4 K at σ 0.56/0.96) and the winter cap
  is still overcast (RH 82–101 % at 64N, cover 0.82, OLR 160 vs 177).
- **Cold free troposphere, traced (2026-09-07, later):** not vertical resolution (L8/L12 at dt 1800:
  same −3..−8 K at 500–700 hPa), not the boundary layer (Betts sub-cloud reference `conv_rh_ref_bl`
  0.85 puts the tropical BL at ERA5's 85–88 % RH but leaves the column 1 K colder; kept, with
  `conv_rh_ref` 0.6). The atmospheric energy budget (new `swn`, `swup` monthly outputs) showed the
  atmosphere absorbing 55 W/m² of solar against ERA5's 79 (surface 190 vs 163), so the troposphere
  cools 19 W/m² less by LW — by being cold. The clear-sky column validation had this all along
  (63 vs 76.5, surface down SW +14): SESAM's near-IR H2O 2-exponential was anchored on the TOA
  net only. Re-anchored on ERA5 columns (SW_A1 0.28, SW_B2 0.033 in aeros_ecckd; SESAM's copy in
  aeros_radiation untouched): clear-sky atm absorption 78, surface down 241.5, TOA net −2.5,
  bands within 5 (residual tilt tropics +5 / poles −4 is ozone/airmass, not water). Coupled:
  atm SW 68 (ERA5 79; the model's column water 20.6 vs 24.5 mm accounts for most of the rest),
  500 hPa bias −3.7 → −2.8 K (tropics) / −4.6 → −3.7 (47N), TOA +2.2 with the same clouds.
  700–850 hPa is still 3–7 K cold with 300 hPa 5–8 K warm outside the tropics: a lapse-rate /
  subsidence structure, next on the list with the column water.
- **Next:** the cold free troposphere is now the model's first-order bias and the reason the
  clear-sky OLR is low. Hypothesis: the tropical boundary layer is too dry (lowest-level RH 80 %
  vs ERA5 87 %, LH 87 W/m² — the "evaporation 30 % low" item) so the SBM reference moist adiabat
  is ~3 K too cold aloft; then the winter-cap moisture/cloud (a drying process for near-saturated
  cold columns), then the consistent FV mass flux (option 3: build the face air fluxes from the
  spectral divergence, r as the metric; r is 0.35 %/step now, 1.2 % at the polar row).





**►► NEXT (in order):**
1. **Re-validate the climate at full-rate physics.** The 1-yr aquaplanet at dt=3600 has OLR
   ~117 W/m² (was ~155), a 188 K level-1 mean (172 K min) and an 81 m/s lid jet — the
   stratosphere now cools at the true radiative rate with no ozone heating. Expect the
   `cond_rh_crit`, albedo and cloud calibration to move. Rebuild the ERA5 comparisons
   (`scripts/rce_validate_era5.jl`, `land_validate_era5.jl`) from scratch.
2. **1-1 SpeedyWeather comparison at real-Earth BCs** (T21, orography + land-sea mask +
   seasonal SST, SW default physics as the yardstick; `mwm/B_multires` has the SW
   PrimitiveWetModel setup, `mwm/C_omega` the diagnostics pattern). Same-hardware SYPD.
3. Characterize the N-cycle dt-sensitivity (30-d OLR 159 at dt=900 vs 144 at 3600: the
   physics time scales are in steps, as in SW) and the dt=5400 failure; try ABBA as default.
4. Slab + sea-ice validation, then T31/T42.

The scratch configs/logs for all of the above are in this session's scratchpad only; the
repro is trivial: `output/land_cal/land_topo_snow.nml` with `time_scheme="leapfrog"`,
`dt=1800`, 30 days.

---

## ►► NEXT SESSION: the descent-concentration residual is ~closed — pick the next M-item

**RESOLVED (`docs/refs/sw_faithful_dynamics.md`).** The residual moist bias traced to aeros's Hadley
descent not concentrating at 30°, which traced to aeros having **no eddy momentum flux** — its RCE
could not sustain baroclinic eddies (seeded eddies blew up at the model top by ~day 27).
Reproducing SpeedyWeather's core numerics fixed it (see below). **At equilibrium (~200–300 d) the
SW-faithful config reaches descent +5.7 hPa/day @ 25° and subtropical RH 72%** — close to SW's
+6.4/62% and a huge move from the pre-fix 94% RH / no-eddies / descent-smeared-25→69°. (Beware: the
descent needs ~200 d to equilibrate; a 100-day snapshot mis-reports it as +3.15/78%.)

**Understood but NOT worth chasing further:** aeros's sustained eddies are *mis-scaled* — 43% of the
eddy KE sits in an m=8 spike (SW is synoptic m=1–5) that is an instability of aeros's too-strong,
narrow subtropical jet (51 m/s @ 30°). BUT this is **cosmetic**: aeros concentrates its descent
*thermally* via the steep APE-Control SST gradient at 30°, where SW does it via eddies — different
routes, similar result. Swapping to SW's smoother cos²(φ) SST (`sst_shape=1.0`) fixes the eddy scale
(m=8 43%→5%) but **weakens/diffuses the descent** (+5.7→+1.3), so the SST/jet is *not* a productive
lever — the APE profile is load-bearing. Closing the last ~0.7 hPa/day and 10% RH to SW is deep
coupled-dynamics work with low marginal value; better to move to the next M-item (land calibration,
seasonal cycle, slab/sea-ice validation) and revisit only if a real application needs it.

**►► THE FIX THAT LANDED THIS SESSION (all on `main`, opt-in / default off, bit-for-bit, 25 tests
pass) — `docs/refs/sw_faithful_dynamics.md`, `hadley_edge_emf.md`:** aeros could not sustain eddies;
the blow-up was a **leapfrog computational mode** in the balanced field (the "model-top" framing was
a red herring). The unlock was **`eps_filter=0.1` — SpeedyWeather's actual Robert-filter
coefficient** (aeros ran 0.06, too weak → blows up; 0.15 is too strong → kills eddies). Plus two
opt-in SW-faithful core options: **`si_alpha`** (semi-implicit decentering; 1.0 = backward-implicit;
commit `7d7484e`) and **`tau_diff_div` + divergence-only order taper** (SW's gravity-wave absorber).
Full SW config = `si_alpha=1, ndiff=8, tau_diff=4, tau_diff_div=1, diff_taper, diff_ndiff_top=4,
diff_taper_sigma=0.2, l_sponge=.false., eps_filter=0.1, raw_alpha=0.53`. Result (T21 RCE, 300 d
stable, eddy KE ~6): descent 25→69° tightened to **19–30°**, jet 14→**31 m/s @ 30°**, subtropical RH
94→**78%** (SW 62%). New diagnostics (commits `402fa7a`, `b0321f7`, `d3d1748`): the eddy
momentum-flux dump (`uvpr`/`ubar`/`vbar`), the eddy zonal-wavenumber spectrum (`accum_espec`, KE +
flux co-spectrum, and KE by latitude/level), and `scripts/hadley_edge_compare.jl`; the SW extraction
(`mwm/C_omega/omega_structure.jl`) mirrors all of them.

**Superseded framing below** (kept for the ruled-out chain): the descent gap was thought to be the
cell-edge with candidates weak-jet / eddy-termination / T21. The eddy-termination piece is now the
active thread (the jet–eddy loop); T21 is not it (SW works at T21, and aeros T42 was a prior null).

**Landed this session (all on `main`, opt-in / default off, bit-for-bit; 25 tests pass):**
`couple_diabatic` (in-solve diabatic coupling, needs eps~0.15; commit `51bd92a`), `vert_vanleer`
(van Leer vertical humidity transport; commit `3de0b69`), the zonal-mean per-term heating dump,
the `q_force` prescribed-heating hook, and the `t_ref` knob. Comparison artifacts:
`mwm/C_omega/` (SpeedyWeather T21 run), `docs/refs/hadley_core_diff_*.md`,
`docs/refs/speedy_omega_structure.md`.

### ✅ DONE this session: the moist-bias column physics (on `main`, commit `a5395d6`)
Diagnosed against SpeedyWeather (same Simplified Betts–Miller) that aeros had two saturating
sinks and no drying sink. Three fixes, all landed, 25 tests pass:
- **vdiff** now diffuses **dry static energy** over a **Richardson-diagnosed BL** (Frierson
  2006; SpeedyWeather's scheme). Was diffusing T → isothermal/over-stable lower troposphere.
- **condensation** defaults **rh_crit=0.95** + **reevaporation** of falling precip (was pinning
  the grid mean at 100%).
- **convection default → `sbm_frierson`** (surface-to-LZB depth + integrated Pt/Pq trigger —
  what SpeedyWeather runs; stays active where the strictly-buoyant band went dormant).
- Result (rotating RCE, cond_rh_crit=0.95): cover 0.999 → **0.70**, TOA −18.6 → **+3.8** W/m²,
  tropics moist-adiabatic and dried to ~73%. `omega` subsidence diagnostic added (commit
  `c8120b8`). **Superseded framing:** the earlier "SBM goes dormant, fix the trigger" scope
  (`docs/frierson_sbm_scope.md`) was only half the story — the trigger port alone did *not*
  fix it; the vdiff + condensation sinks were the load-bearing column-physics fixes, and the
  true remaining lever is the weak circulation above.

### Landed this session (all on `main`, pushed)
- **Land×cloud blow-up FIXED.** Root cause was *not* the ocean/model-top (my first handoff
  framing was wrong — corrected below): it is a **land × prognostic-cloud interaction** — the
  thin, patchy Sundqvist cover over land's zonal asymmetry drove a near-equatorial jet.
  Fixed by **calibrating** the cloud scheme (cover 0.07 → 0.56 ≈ ERA5 0.63; slower 3 h
  timescale de-patchifies the field). `cf` advection was added too (it's correct physics but
  was stability-neutral here — calibration is what fixed it).
- **`aeros_moisture` → `aeros_transport`** — the FV moisture core generalized into a reusable
  positive-definite tracer engine (qv bit-identical); `cf` now advects through it. Restart is
  bit-identical with land+cloud+sea-ice all active.
- **Adaptive hyperdiffusion** (`diff_adapt`, `diff_taper`; both default OFF, bit-for-bit) — a
  SpeedyWeather-style numerical safety net (vorticity-scaled strength + σ-tapered order). It
  independently holds the land×cloud jet bounded (max|u| ~10 m/s); the **order-taper** is what
  bites (a thermally-driven mode; vorticity-scaling alone can't touch it). Available if the
  Frierson work needs a top-stability net.
- **Suite now 25 tests**, all pass; `make all openmp=1` clean.
- **Reference studies:** `docs/refs/speedy_comparison.md` (SpeedyWeather/SPEEDY vs aeros) and
  `docs/refs/reduced_grid_feasibility.md` (reduced grid = **NO-GO**: SHTns is regular-grid-
  only; and the coupled model is physics-bound, not transform-bound, so the reachable speed
  lever is physics, established by profiling — not the transform).

---

## Superseded: the land×cloud "cloud calibration" task (DONE — see above)

**CORRECTION (diagnosed 2026-07-28, supersedes the "slab model-top instability"
framing that was here):** the coupled blow-up is **not** an ocean/model-top
problem. It is a **land × prognostic-cloud interaction**, and the ocean is a
bystander. Bisection with `rce_long` (`logs/topinvestig/`):
- Neither feature alone blows (land-only, clouds-only: stable ≥19 days).
- **Only `l_land + l_cloud_prog` together blow up** — an explosively-growing
  *near-equatorial upper-tropospheric jet* (equator, lev 3–5: ~10 m/s day 12 → 31
  m/s day 14 → top NaN day 15). Sea ice is irrelevant (land+seaice, clouds+seaice
  both stable); **the ocean is irrelevant** — prescribed SST blows at the same
  step as a 2 m or 10 m slab (bit-identical timing).
- **It is the too-thin, patchy prognostic cover.** With land: diagnostic clouds →
  stable; clear-sky → stable; prognostic **thin** (default knobs) → blows;
  prognostic **forced thick** (`cloud_rhc_sfc=0.40, cloud_rhc_top=0.60`) → stable.
  Without land, prognostic clouds are stable. The prognostic scheme touches the
  model **only through radiation** (`aeros_cloud_prog_apply` is `intent(in)` on
  T/q; driver forces `rad%clouds=.TRUE.` when `l_cloud_prog`, rce_long.f90:205), so
  the mechanism is: land's zonal asymmetry × a thin, *spatially-patchy* cloud field
  → patchy cloud-radiative heating near the equator → a resonant jet nothing damps.
  This is exactly the "default knobs give thin cover; calibration is future work"
  caveat the cloud scheme shipped with — now shown to be a *stability* problem, not
  just a bias, because the scheme was only ever tested on the zonally-symmetric
  aquaplanet.

**So the bare slab (incl. sea ice) is NOT the blocker** — a 2 m slab with hard
damping runs stably ≥19 days, and land+sea-ice is stable. The §29 slab *thermo-
dynamic* runaway (cloud cover 0.66→0.86 over years) is a separate, slower issue.

**Fix direction (get AR's call before implementing — cloud-scheme design):**
1. **Calibrate the prognostic cloud scheme to realistic cover (~0.6)** — the real
   fix; forced-thick is already stable. Tune the Sundqvist RH_crit profile /
   formation-evaporation balance so equilibrium cover matches obs instead of ~0.07.
2. **Add horizontal advection of `cf`** (the agent deferred it) — an unadvected
   cloud field is patchier than a real one; transport would smooth the forcing.
3. **Consider a limiter/smoother on the cloud-radiative heating** if the equatorial
   jet proves partly numerical (the day-12→14 explosive acceleration looks resonant).
Recommended: (1)+(2) together; (3) only if a jet remains after realistic cover.
Repro configs in `logs/topinvestig/` (`L_progcloud` blows, `L_progthick`/`A_progcloud`
stable).

**Wave 2 landed (on `main`, integrated + verified together).** Three features
built in parallel worktrees, merged with an albedo-unification pass, all opt-in
(default off → the 21 prior tests bit-for-bit), suite now **24** tests, all pass:

- **Land surface + land–sea mask (`l_land`, feat/land).** `src/physics/
  aeros_land.f90`: LSM + land albedo from ERA5 (`lsm`/`fal`) via `aeros_bcinput`;
  bucket soil moisture `w` + slab soil temperature `t_soil` (the land skin temp,
  replacing SST on land); β-limited evapotranspiration. Surface/radiation branch
  on the mask; ocean path bit-for-bit. Smoke run: 34% land, deserts hot/dry.
- **Prognostic cloud fraction (`l_cloud_prog`, feat/clouds).** `src/physics/
  aeros_cloud_prog.f90`: Sundqvist cf as a gridpoint prognostic `now%cf_g`
  (source = large-scale + convective detrainment, sink = evaporation). cf-only,
  in-cloud water diagnosed by the existing tuned optics (`aeros_cloud_water`).
  On the *same* atmospheric state the prognostic sink holds cover ~7× below the
  runaway diagnostic — the runaway mechanism is addressed. **Horizontal advection
  of cf deferred** (documented); default knobs give a thin cover — calibration is
  future work.
- **Sea ice (`l_seaice`, feat/seaice).** Semtner 0-layer thermodynamic ice in
  `aeros_ocean` (thickness/fraction/ice-surface-temp), ice-albedo feedback,
  replacing the freeze-floor clamp. Energy-conserving to 1e-16; ice forms poleward
  in a cold slab smoke run.
- **All three extend the restart** (append-only, gated by `*_present` attrs) and
  feed a **single unified per-cell surface-albedo field `rad%alb_map`** — land and
  ice compose into it (disjoint cells; ice never on land). Both land and sea-ice
  had independently added a per-cell albedo to radiation; the merge unified them.
- **Verified together.** 24/24 tests; a combined **land+clouds+sea-ice restart is
  bit-identical** for all meaningful state (spectral, moisture, land soil, ocean,
  radiation cache; every cloud cell with cf>1e-12 exactly 0.0 diff). Residual: 62
  cf_g cells at ~1e-96 (denormal-tail underflow, physically zero cloud) — benign,
  not clamped.

**►► The blocker that Wave 2 exposed and made critical: the slab ocean is
model-top unstable within ~15 days.** Sea ice *requires* a slab ocean
(`ocean_mode=1`), but a coupled slab run blows up at level 1 (the model top,
axisymmetric — the §29 slab runaway / the `rce-instability-diagnosis` note) in
~10–15 days, while ice needs ~75 days to form. **So sea-ice thermodynamics cannot
be exercised in a coupled run, and no slab/coupled climate can be validated, until
this is fixed.** It is no longer deferrable. Likely levers (in order): a proper
model-top treatment (higher top / more levels / gravity-wave drag — the same fix
topography wants), stronger/deeper sponge as an interim, and the free-SST energy
imbalance itself (§29: cloud over-reflection + no realistic cover). Do this next.

**Then:** land calibration + seasonal cycle on; cf horizontal advection +
knob calibration; sea-ice validated once the slab is stable; resolution/eddy work.

Original Wave-1 status and the full leverage analysis follow.

---

## Wave 1 recap — topography / restart (DONE)

**Wave 1 landed (on `main`, integrated + verified together).** Two features were
built in parallel worktrees, merged, and validated as a pair:

- **Topography (orography) forcing — DONE.** `src/aeros_bcinput.f90` (a general,
  reusable lon/lat→Gaussian bilinear regridder — the foundation for *all* future
  changing boundary conditions: LSM, albedo, SST, ice topo) + `src/physics/
  aeros_topography.f90` (ERA5 `z` → phis, no ÷g). Wired into `rce_long` behind
  `l_topography` (default `.false.` → flat aquaplanet, bit-for-bit unchanged),
  `topo_file`, `topo_ramp_days`. The orography is **spectrally truncated to T21**
  before it feeds the dynamics (raw grid-scale phis aliases and rings the model
  top — standard spectral-model practice), and ramped in over `topo_ramp_days` via
  a **pure function of absolute model time** (restart-safe). A 2-yr real-Earth run
  is NaN-free, TOA net ~−2 W/m² with no secular drift; `docs/figures/topo_phis.png`
  shows continents/mountains in the right places. **Gotcha:** topography-on needs
  *stronger top damping* than the aquaplanet (blow-up recurs at ~0.88 of full
  amplitude regardless of ramp length — the lever is top damping, not ramp
  duration). The validation run used `tau_diff=3.0`, `sponge_sigma=0.20`,
  `sponge_kr/kt ×2` — a run-config choice; **code defaults are untouched.** A
  proper gravity-wave-drag / higher-top fix is the durable answer (Wave 2+).
- **Restart / checkpoint — DONE.** `aeros_timestep_write_restart` /
  `aeros_timestep_read_restart` (in `aeros_timestep.f90`) serialize the *complete*
  integrator state to one netCDF file: both leapfrog spectral levels
  (`now`=Xⁿ, `old`=Xⁿ⁻¹), grid humidity `qv_g`, **ocean SST+fnet (spun-up ocean)**,
  the scalars `nstep`/`mass_target`/`lnr_cum`/`lnr_last`/time, and the radiation
  cache. Note the cache is more than `rad%heat`: the slab ocean consumes the cached
  radiative *surface fluxes* every step, so those are saved too. Wired into
  `rce_long` behind `restart_in` / `restart_out` / `restart_interval` (all default
  to cold-start). Metadata (nlm/nlev/nlon/nlat) is validated on read; `mass_target`
  is restored explicitly (not recaptured).
- **Both verified — separately and together.** `test_restart` (bit-exact split-run)
  and `test_topography` (regridder exactness + ramp) are acceptance tests **20 and
  21**; all 21 pass, `make all openmp=1` clean. Integration check: a continuous
  40-step topo-on run vs. a 20-step + checkpoint + resume-to-40 run are
  **bit-identical** across all 29 restart fields with the topography ramp active
  across the checkpoint boundary.

**Wave 2 (recommended order, unchanged from the leverage analysis below):**
1. **Land surface + land–sea mask** — reuse `aeros_bcinput` to read the ERA5 `lsm`
   / `fal` fields; add land albedo/roughness/heat-capacity maps + a soil-moisture/
   temperature scheme (bucket is the pragmatic first cut) + evapotranspiration.
   This is what creates monsoons and the subtropical dry zones.
2. **Prognostic cloud fraction** — replace the diagnostic RH→cover that runs away
   (0.66→0.86) in slab runs. Sundqvist-type prognostic `cf` (source = large-scale
   + convective detrainment, sink = evaporation) fits the L12 column physics.
3. **Sea-ice thermodynamics + albedo** — replace the ocean freeze-floor clamp;
   note this edits `aeros_ocean`, whose state is now serialized by restart — extend
   `write/read_restart` when ice state is added.

Original leverage analysis and full M2 detail follow, unchanged.

---

## Topography (Wave 1) — original brief, now DONE

The ecCKD revalidation and the radiation-cost fix are **done** (§29, on `main`,
pushed). Short version: ecCKD is the fast default, the coupled RCE needed **no
re-tune** (`cond_rh_crit=0.93` still best, TOA net +5.6 W/m², clear-sky OLR −5.1
with no `LW_VAP_OPAC`), and radiation is now **~6× faster / bit-exact** (the
per-column k-table was cached; cadence 3 h→6 h). **~1656 model-yr/day** at T21L12
on 10 cores — a multi-decade coupled run is minutes.

**The model is a validated aquaplanet with no Earth forcing.** Every ERA5
mismatch that remains (near-saturated RH, weak jet, and the slab-ocean cloud
runaway below) traces to **missing boundary conditions**, not to the column
physics or radiation. So stop tuning against ERA5 and start adding the forcings —
in leverage order.

**Do topography first.** Highest leverage, and the hook already exists:
`aeros_timestep_set_phis(ts, phis2)` is called in `drivers/rce_long.f90` but fed
**zeros** (`phis2 = 0.0`). Feeding a real surface-geopotential field unlocks
stationary waves, orographic precip, and regional structure — the first thing that
makes zonal-mean comparisons to ERA5 fair.
- Read an orography field (ERA5 `z`/geopotential at the surface, on the same
  144×73 grid as the Tier-1/2 data) → `phis2`.
- Balance the initial state to the topography (avoid a spin-up shock); check the
  run stays NaN-free and mass/energy close.
- For paleo (the insol/Laskar orbital path, §27), this is also where LGM ice-sheet
  topography enters.

**Then, in order:** land surface + land–sea mask (the real build — soil
moisture/temp, land albedo/roughness, evapotranspiration; this is what creates
monsoons and the subtropical dry zones) → **seasonal cycle *on*** in coupled runs
(infrastructure exists, §27/M3a; RCE currently runs annual-mean) → resolution
bump + eddy deficit (T21 jet ~25% of ERA5; §25–26 showed T42 alone didn't fix it)
→ gravity-wave drag → cloud/sea-ice/aerosol/ozone refinement.

**Two gaps are now on the critical path (not deferrable), exposed by §29 item 4.**
The 5-yr slab-ocean run does **not** reach a physical equilibrium: free SST warms
the tropics, and the **near-saturated moist bias makes cloud cover run away
0.66→0.86** (SW over-reflection → TOA net −13.5 W/m²), while the **freeze-floor
slab (no sea ice)** blocks high-latitude closure. So:
- **A real cloud-fraction scheme** (prognostic, not the diagnostic RH→cover that
  saturates) is needed before any slab/coupled climate is trustworthy.
- **Sea-ice thermodynamics + albedo** (replacing the freeze-floor clamp).
These were "Tier-3 refinements"; the slab run promotes them.

**Tools/gotchas (still current):** `rce_long` runs ecCKD by default; cadence and
scheme are now namelist-optional (`rad_interval`, `rad_scheme`) defaulting to
6 h/ecCKD via `input/rce_defaults.nml` (a keyless namelist inherits them, doesn't
error). For a SESAM A/B set `rad_scheme = 1`. `validate_era5 ... ecckd` selects
ecCKD offline (arg 3 defaults to `sesam`). Worktrees need the **`insol`** and
**`fesm-utils`** symlinks (→ `~/models/`, gitignored) and a copy of the generated
`Makefile` — and **parallel worktree agents must not use `git stash`** (the stash
stack is global across worktrees and collides; build a separate baseline binary
instead). Always **build + run the merged tree yourself** after integrating agent
branches — a clean auto-merge can still be semantically broken or a silent
perf regression (both happened in §29).

---

Supersedes the earlier M2.5a–e handoff. `make
all openmp=1` is clean, and **all 19 acceptance tests pass**. Two blockers that
this doc used to open with are now cleared:

- **The M2 RCE runaway** — diagnosed and resolved (surface momentum drag + dry
  convective adjustment); detail below, unchanged.
- **The Tier-2 ERA5 data** — the cloud/RH/u/v fields that were missing are now
  delivered (`~/data/era5/monthly-{pressure,single}-levels`, 1991–2020 monthly
  clim on the same 144×73×37 grid as Tier-1). This unblocked the cloudy-sky
  radiation work.

**Latest work (M2.5f, this doc's most recent): the cloudy-sky radiation branch
is built and validated offline against ERA5** (m2_results §17–18). All-sky LW
(`aeros_lw_cloudy_column`, max-overlap run-twice blend) and SW
(`aeros_sw_cloudy_column`, cloud reflector + blend) operators, resolved per-layer
grey-cloud optics, opt-in (`cf=0` → clear-sky bit-for-bit, all tests untouched).
Driven on ERA5's `cc`/`clwc`/`ciwc`, the modelled cloud radiative effect matches
ERA5 to a few W/m² (TOA LW +18.2/+21.8, SW −44.5/−46.1, net −26.3/−24.2) with
faithful patterns (`docs/figures/era5_cre_validation.png`). **Not yet wired into
the coupled model** — that (a diagnostic cloud field in `aeros_radiation_apply`)
is the next step, see below.

## The RCE runaway is resolved

The coupled radiative-convective run blew up (~day 63–186 depending on config).
Diagnosed end-to-end with new instrumentation (below); **five suspects were ruled
out by direct test** — a radiation vertical "sawtooth" (turned out to be a bug in
the diagnostic itself), OLR saturation (refuted by a clean `probe_lw` dOLR/dTs ≈
2.75 W/m²/K sweep), the convection closure, global energy balance, and the
model-top sponge. **Two real root causes, both now addressed:**

1. **Missing surface momentum drag → the rotating blow-up (the main one).** The
   RCE baroclinicity is *surface-trapped* (gradient + jet peak at the lowest
   level, decaying upward), so it spins a low-level jet — and nothing braked it
   (`aeros_surface` did heat/moisture only; `aeros_vdiff` was zero-flux at the
   surface; the sponge is top-only). HS stayed bounded only because it kept
   explicit low-level Rayleigh friction, dropped when the RCE replaced HS forcing.
   **Fix (committed):** a bulk surface stress `τ = ρ c_d |u| u`, imposed
   *implicitly* as the momentum bottom BC of `aeros_vdiff` (forward-split momentum
   on the leapfrog is unstable — the vdiff lesson). `surf%c_d` (0 = off/
   bit-unchanged; ~1.5e-3 to enable). → the **rotating RCE now completes 200 days
   NaN-free**, jet at a realistic upper level (~16 m/s, 25–36°), T_low ~292–300 K,
   **without eddies**. So "needs 3D eddies" was wrong; the defect was the momentum
   sink.
2. **No dry convective adjustment → a hot dry surface heat-trap.** Once surface
   fluxes heat + dry the lowest layer (RH→~0), the *moist* SBM buoyancy test finds
   it stable and the convecting band detaches upward, leaving a super-adiabatic
   layer no moist scheme ventilates. **Fix (committed, opt-in):** standalone dry
   convective adjustment before the moist scheme, relaxed over τ (`cnv%dry_adjust`,
   off by default). Correct physics; ventilates the layer; not required to bound
   the run once drag is on.

## Two bounded-RCE vehicles now exist (`drivers/rce_long.f90`)

- **Rotating, realistic (preferred):** `c_d > 0` (+ default rotation/gradient).
  Real equator-pole gradient and upper-level jet. The right vehicle for ERA5.
- **Non-rotating uniform:** `l_nonrotating` + `l_uniform_insol` (uniform IC). No
  jet at all (f=0); isolates column physics. Confirmed the column physics is
  sound (T_low bounded, OLR responds).

Both leave a residual **+15–17 W/m² net TOA** (OLR ~20 W/m² too low — the Tier-1
"too opaque" bias) → slow secondary warming, not a blow-up. Balance it with
albedo ~0.47 or an LW-opacity fix for a true steady equilibrium.

## What landed in the RCE-resolution session (all on `main`, pushed)

- **Surface momentum drag** (`surf%c_d`, implicit via vdiff). The rotating-RCE fix.
- **Dry convective adjustment** (`cnv%dry_adjust`, opt-in, relaxed over τ).
- **Convection single-layer-band fix** (skip `ktop==kb`: not overturning).
- **Diagnostic suite in `rce_long`:** per-term heating split (`enable_diag`),
  buoyancy `h_b` vs `h*_env` + RH, hot-latitude local TOA, eddy (RMS T′/KE/[v′T′]),
  baroclinicity vertical structure, `max|u|` location; a per-term capture-bug fix.
- **`drivers/probe_lw.f90`** (`make probe-lw`): offline clear-sky radiation column
  harness — the tool that exonerated radiation. Reusable.
- **`aeros_timestep_set_sponge`** (sponge strengths tunable after init).

New `rce_long` namelist knobs: `l_diag`, `l_dry_adjust`, `l_uniform_insol`,
`l_nonrotating`, `c_d`, `sponge_kr/kt/sigma`.

## Next steps (recommended order)

1. **Balance TOA in the rotating vehicle** (albedo ~0.47, or start the OLR-opacity
   fix) → a *true* steady RCE equilibrium, the reference state for validation.
2. **Do eddies now grow? — ANSWERED: no, and it is not a resolution problem**
   (§25). Seeded (`seed_asym>0`) baroclinic eddies decay from the seed and
   saturate at eddy KE ~1 m²/s² (~100× below the observed storm track) with
   `[v'T']` at noise level; the jet rises only from ~7 to ~10 m/s (ERA5 27–31).
   Slackening ∇⁶ 4× and going to **T42 both give the same answer** — so the weak
   jet is neither a hyperdiffusion artifact nor a truncation limit. The limiter is
   the RCE's weak, surface-trapped baroclinicity; closing the jet gap needs
   stronger baroclinic forcing (meridional SST gradient / seasonal cycle / land),
   not more resolution. **Closed as a characterized, expected limitation (§26):**
   a thermal-wind diagnosis shows the surface T gradient matches ERA5 but collapses
   aloft (no eddy ventilation of the midlat upper troposphere — the same missing-
   eddy loop). Crucially this is the *expected* behaviour of a perpetual-annual-mean
   aquaplanet, which omits the land–sea contrast, stationary waves and seasonal
   cycle that drive Earth's storm tracks; the fair circulation yardstick is the APE
   aquaplanet intercomparison, not ERA5. The thermal/moisture columns (§22–24),
   where the ERA5 comparison *is* fair, validate well. Not a tuning target.
3. **Model humidity bias in the coupled RCE — RESOLVED for now** (§22–23).
   Localized against ERA5 (`scripts/rce_humidity_vs_era5.jl`, `rce_long` RH dump):
   the RCE was near-saturated at nearly all lat/height (+40% RH vs ERA5), both the
   rotating and the non-rotating-uniform vehicle alike, dried only where deep
   convection fires — so **column physics, not missing subsidence**. Root cause:
   convection shuts itself off at equilibrium (`hb < h*_env` everywhere) and
   large-scale condensation then pins RH at 100% with no drying sink. Fixed by the
   sub-grid-saturation condensation ceiling `cond_rh_crit` (now an `rce_long`
   namelist knob; module default stays 1.0). At **`cond_rh_crit = 0.93`** the
   overcast is gone — rotating-vehicle TOA net −36 → **+8 W/m²**, cover **0.65 vs
   ERA5 0.63**. *Limitation:* it is a ceiling, not subsidence — it buys balanced
   TOA + realistic cover but leaves free-trop RH ~85–90% (ERA5 ~45%); the residual
   RH bias and the cover's flat latitudinal structure need a subsidence-drying
   treatment (option C), deferred.
4. **The clear-sky OLR bias** (§14, §20): addressed with a single documented
   `LW_VAP_OPAC` correction (0.80) — OLR bias −15.9 → −7.1, and it improved the
   cloudy LW CRE too (−3.6 → −2.1). A single knob can't null both OLR and
   surface-LW (they trade off), so a residual −7 OLR remains; the structural fix
   is ecCKD. Not a blocker.
5. **ERA5 moist-line (Tier 2) — DONE** (§24). The coupled RCE (rotating,
   `cond_rh_crit=0.93`) validated against ERA5 zonal-mean T/u/RH via
   `scripts/rce_validate_era5.jl` (`docs/figures/rce_validate_era5.png`). Verdict:
   gross thermal structure right (mid/lower-trop |T bias| < 3 K), jet at the right
   latitude/level (±36°, 198 hPa) but only ~25% of ERA5 strength (7.4 vs 27–31
   m/s), RH still too moist (the §23 ceiling), plus an upper-trop warm bias (+12 K
   @ 198 hPa, coarse L12) and a polar warm bias (prescribed-SST artifact, no cold
   poles). Every deficit is physically attributable. **Natural next: step 2** —
   the weak jet is exactly the missing-eddy signature; seed a T21 asymmetry and
   see if baroclinic eddies grow and spin the jet toward ERA5 strength, else T42.

## Gotchas the next session will want

- **`rce_long` requires every namelist key it reads to be present** — `nml_read`
  (no `defaults_file` in this driver) is *fatal* on a missing key ("parameter not
  found"). When you add a knob, add it to every `logs/rce_*.nml` you run.
- **`probe_lw`** is a diagnostic driver, not a test. `make probe-lw`.
- **`c_d=0` and `cnv%dry_adjust=.FALSE.` are the defaults** — every prior run/test
  is bit-unchanged; the fixes are opt-in.
- **The `Makefile` is generated** from `config/Makefile` by `configme` and is
  gitignored. Add driver/test rules to **`config/Makefile`** (the committed
  template), not just the root `Makefile`.
- **Worktrees** (per the global worktree discipline) need two things the fresh
  worktree lacks: copy the root `Makefile` in, and symlink `fesm-utils` →
  `/Users/alrobi001/models/fesm-utils`. Both are gitignored.
- **macOS:** `script` has no `-qfc` (that's util-linux — logs silently fail); and
  `timeout` is `gtimeout`. For live logs, redirect stdout to a file directly.
- **The RCE driver is `drivers/rce_long.f90`** — namelist-configured, heavily
  instrumented; a diagnostic tool, not an acceptance test.
- **Slab (`ocean_mode=1`) is unstable on the OLD core** — the pre-SW-faithful slab
  recipe (`rce_slab5yr.nml`) now NaNs ~day 89 (git-bisected to `a5395d6`'s default
  physics: Richardson-BL vdiff + precip reevaporation, validated on prescribed-SST
  only, never the slab). Not a bug in the opt-in features (SW-faithful core and
  `sst_shape` are bit-for-bit). Fix: run slab on the **SW-faithful recipe**
  (`si_alpha=1`, `ndiff=8`/`tau_diff=4`, `tau_diff_div=1`, `diff_taper=T`,
  `l_sponge=F`, `eps_filter=0.1`), which is stable with the new physics on; or, on
  the old core, set `vdiff_richardson=.false.` **or** `cond_reevap=0`. No slab
  stability test exists yet — worth adding.

## Carried M2 items (not blocking)

Cloudy-sky radiation (#3 above); CLIMBER-X ocean (slab is the interim); shallow
convection is temperature-only (§12.2); Manabe multi-layer adjustment; vertical
van Leer; per-row polar sub-cycling; `tau_diff`/∇⁶ retune.

**In progress / done (capabilities):**
- **Seasonal cycle + Laskar orbital insolation — DONE (§27).** `aeros_insolation`
  wraps the fesmc/insol package; namelist `time_bp` selects the epoch. Validated
  (`probe_insol`): annual-global mean = S0/4, polar day/night, live Milankovitch
  redistribution at 21 ka. All 18 tests pass. Needs the `insol` symlink (like
  `fesm-utils`) in each worktree.
- **ecCKD correlated-k radiation — DONE (§28), merged.** Complete clear + all-sky,
  LW + SW, opt-in behind `SCHEME_ECCKD`; SESAM stays the bit-identical default.
  Compact Malkmus/Goody band model (analytic inverse-Gaussian k-distribution), 15
  g-points, no external data table. Beats SESAM on OLR (−5.1), CO2 forcing (3.6),
  and LW CRE (−0.1) with `LW_VAP_OPAC` retired; cost 0.84–1.78× SESAM at the 1–3 h
  cadence. New test `test_ecckd` (suite now **19**). **Recommended scheme going
  forward — switching the default from SESAM to ecCKD is a deliberate call still to
  be made** (it changes all results / breaks bit-reproducibility of past runs).

## Land calibration handoff (2026-07-30) — two open tasks for the next session

**State.** Coupled land runs work and are validated against ERA5. The seasonal
cycle (prescribed 2D ERA5 SST + seasonal insolation + real land mask, flat) runs
stably and land fixes most of the aquaplanet's missing seasonality (NH-midlat
amplitude ~9.5 K aquaplanet → ~21 K land; ERA5 ~37 K). All on branch `land-cal`,
merged to main. See [[model-top-instability]], [[slab-instability-a5395d6]].

**THE dt=900 REQUIREMENT (critical).** Coupled/asymmetric runs (land, or seeded
eddies, or slab) blow up at the model top on the SW-faithful core at the usual
`dt=1800` — it is a **CFL/time-stability limit**, not diffusion or spectral-pole
(git-diagnosed: every damping/inertia lever only *delayed* it; halving `dt`→900
*cures* it, 4 yr steady). So every asymmetric run must use **`dt=900`** (≈2× cost).
The aquaplanet-symmetric seasonal run is fine at 1800.

**Infra in place (all opt-in, default-off, 25 tests pass):**
- 2D seasonal SST: `scripts/make_seasonal_sst_2d.jl` → `input/sst_seasonal_2d.nc`
  (ocean-masked, land-filled ERA5 skt); `aeros_bcinput_read_field(...,itime=m)`
  regrids per month; driver dispatches on the `sst` var rank.
- `monthly_out` now also writes `t2m(lon,lat,month)` (2D near-surface air T).
- `scripts/land_validate_era5.jl` — land-masked t2m amplitude vs ERA5.
- Reference run + configs: `output/land_cal/land_csoil2e6.nml` (dt=900, 4 yr,
  accumulate last 2). `c_soil` sweep nmls alongside.

**Task 1 — orography. ✅ DONE (2026-07-30).** `l_topography=.true.` on the
calibrated flat-land dt=900 config (`output/land_cal/land_topo.nml`), 4 yr,
NaN-free. **Stability-neutral:** the model-top jet (`max|u|` lev 1, 61–84 m/s
over yr 3–4) is bounded/stationary and matches the flat config's own top-jet
regime (flat vs topo lev-1 baroclinicity `peak|[u]|` 52 vs 50 m/s at day 90) —
topography added no stress, no shorter-dt/longer-ramp needed. T21 phis loaded OK
(max ~5.0e4 m²/s² ≈ 5100 m; expected Gibbs undershoot at the min). **Seasonal-cycle
validation** (`land_validate_era5.jl`, `docs/figures/land_t2m_amp_{flat,topo}.png`):
topography **helps the NH mid/high lats** — land t2m amplitude at 50°N 18.7→**25.0**
(+6.3 K), 65°N 21.6→**24.6** (+3.0), toward ERA5 (32.2 / 37.5); elevated, isolated
continental interiors cool harder in winter. SH −35° drops 11.4→8.8 (little SH
land, sensitive to the few narrow continents), so the global land-mean edges down
10.1→9.5. The residual 65°N gap (24.6 vs 37.5) is Task 2's snow-albedo lever.

**Task 2 — close the NH land amplitude gap (missing physics, NOT c_soil).**
`c_soil` is NOT the lever — a 10× sweep (2e6→2e5) moved the land-mean amplitude
<0.5 K; the soil already equilibrates. The NH-midlat gap (model 21 vs ERA5 37 K at
65°N; model winters too warm + ~2-month phase lag) is structural:

  (a) **Snow-albedo feedback — ✅ DONE (2026-07-30).** Diagnostic, `t_soil`-keyed:
      land albedo ramps from the bare ERA5 `fal` value toward `snow_albedo` as the
      skin temp falls below T0, over `snow_dt` K (`f_snow` 0 at/above T0, 1 at
      T0−snow_dt). Refreshed each step by `aeros_land_couple_radiation`; the
      timestep albedo block now fires on land-snow too, not only `l_seaice`. New
      `&rce` knobs `snow_albedo` (default −1 = **off, bit-for-bit**) / `snow_dt`
      (5 K); needs `freeze_floor=.FALSE.` (the default) so t_soil can go sub-zero.
      No prognostic snow / no hysteresis (v1). Test: `test_land` §6b (bump on/off,
      ramp midpoint, ocean untouched). 25/25 tests pass. Config
      `output/land_cal/land_topo_snow.nml` (α_snow=0.6), 4 yr, NaN-free, top-jet
      regime unchanged. **Result** (`docs/figures/land_t2m_amp_toposnow.png`):
      amplitude 65°N 24.6→**26.0**, 50°N 25.0→**26.7**, −35° 8.8→**10.8** (SH
      recovered), land-mean 9.5→**10.3** (ERA5 14.9). Correct-signed but **modest**,
      and the reason is diagnosed: ERA5 `fal` is **already 0.36 at 65°N** (the
      annual-mean field is snow-contaminated), so the bump to 0.6 has little
      headroom exactly where the gap is largest (50°N, fal 0.18, has more room and
      moves more). **The residual 65°N gap (26 vs 37.5) is therefore NOT an albedo
      problem — it is (b).** Follow-ups if pursued: raise α_snow to 0.7–0.8 (cheap,
      headroom-limited); or prescribe ERA5 *seasonal* `fal` per-month (it is a
      12-month field) as the bare baseline with the snow feedback on top (a
      boundary condition + feedback, more faithful to the observed albedo cycle).

  (b) **Lowest model level ~170 m, not 2 m → misses the winter surface inversion.**
      Now the leading residual for the NH high-lat gap (albedo is spent). Partly a
      fair-comparison artifact; a proper 2 m surface-layer diagnostic (from the
      surface flux + stability) is the next lever, not more albedo tuning.
