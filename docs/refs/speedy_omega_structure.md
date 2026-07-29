# SpeedyWeather.jl T21 aquaplanet — zonal-mean overturning structure

Native SpeedyWeather.jl circulation at **T21, L8, moist aquaplanet**, extracted
for an apples-to-apples comparison against aeros's own T21 RCE aquaplanet
overturning. Answers the standing ambiguity: **is the quoted "~13 hPa/day"
SpeedyWeather subtropical subsidence the DESCENT or the peak ASCENT?**

**Short answer: it is the subtropical DESCENT.** SpeedyWeather's zonal-mean
subtropical subsidence is ω ≈ **+6 hPa/day (hemispherically symmetrized) to
+10–11 hPa/day (stronger hemisphere)** at ~30°, ~500–700 hPa. Its peak *tropical
ascent* is only ω ≈ **−2 to −3.4 hPa/day** at ±14–19°. The 13 number is descent;
the ascent is nowhere near it.

---

## 1. Setup (what was run)

- Model: `PrimitiveWetModel` (SpeedyWeather.jl **v0.21.1**,
  `~/.julia/packages/SpeedyWeather/i8kOF/`), CPU/Float32.
- Truncation / vertical: **T21** (`trunc=21`), **L8** equally-spaced σ
  (σ_full = 0.0625, 0.1875, 0.3125, 0.4375, 0.5625, 0.6875, 0.8125, 0.9375),
  octahedral Gaussian grid (32-ring, 1600 points). Δt = 1 h.
- Aquaplanet:
  - `ocean = AquaPlanet(spectral_grid)` — SST constant in time & longitude,
    **SST(φ) = (302−273)·cos²φ + 273 K** (Te=302 K equator, Tp=273 K poles). This
    is SpeedyWeather's default aquaplanet SST; noted here as required.
  - `land_sea_mask = AquaPlanetMask(spectral_grid)` — all sea, no land.
- Everything else default and standard: convection = **BettsMillerConvection**
  (SpeedyWeather's Simplified Betts–Miller, Frierson 2007; rhbm≈0.8),
  large-scale condensation = **ImplicitCondensation** (RH threshold 0.95),
  one-band SW/LW radiation, bulk-Richardson surface fluxes + vertical diffusion.
- Spin-up **400 days**, then time-average the next **800 days** (19 200 hourly
  samples). An earlier 300+150-day run gave the same structure; 800-day
  averaging was used to beat down transient-eddy noise and hemispheric wobble.

Run/extraction scripts: `mwm/C_omega/omega_structure.jl` (+ `postproc.jl`,
`Project.toml`); raw output `mwm/C_omega/output/speedy_omega_T21L8.nc` and
`summary.txt`.

## 2. Method (how ω, heating, RH were extracted)

- **ω [hPa/day], >0 = subsidence** (matches aeros's convention). SpeedyWeather's
  σ-velocity lives in `vars.dynamics.w` = **a·σ̇** at half levels (radius-scaled,
  "positive down"). So σ̇ = w / a (a = planet radius 6.371e6 m), and
  ω = pₛ·σ̇ (dominant term; the σ·∂pₛ/∂t term ≈ 0 in the time mean), with
  pₛ = exp(`vars.grid.pressure`) [Pa]. Full-level ω is the half-level average.
  ω accumulated per timestep and time-averaged (Eulerian zonal-mean ω = the
  Hadley-cell vertical velocity).
- **Relative humidity [%]**: RH = q / q_sat(T, σ·pₛ) using SpeedyWeather's own
  `saturation_humidity`, accumulated per step.
- **Precipitation [mm/day]**, convective vs large-scale, from the instantaneous
  `rain_rate_convection` / `rain_rate_large_scale` diagnostics — the exact
  column-integrated latent heating and the cleanest ITCZ locator.
- **Net diabatic (physics) heating [K/day]**: re-run `parameterization_tendencies!`
  each step and read `tendencies.grid.temperature`/a. This is latent + radiative
  + turbulent; in the tropical mid-troposphere it is latent-heating-dominated. It
  is weak and noisy (|<1.3| K/day) because latent heating and radiative cooling
  largely cancel in the zonal+time mean, so precip (above) is the primary
  latent-heating diagnostic. Full field in the NetCDF.
- All native-grid fields regridded (via spectral) to a FullGaussianGrid 64×32 and
  zonal-averaged. Tables below are **hemispherically symmetrized** onto |lat|
  (the symmetric-SST aquaplanet develops a persistent hemispheric asymmetry — the
  stronger-hemisphere numbers are given where they matter).

---

## 3. ω(|lat|, σ) — hemispherically symmetrized, hPa/day (>0 subsidence)

| \|lat\| | σ=.062 | .188 | .312 | .438 | .562 | .688 | .812 | .938 | precip mm/d |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|  2.8 | −0.2 |  0.2 |  0.1 | −0.4 | −0.4 | −0.4 | −0.2 |  0.0 | 2.83 |
|  8.3 | −0.1 | −0.2 | −1.5 | −1.7 | −0.7 | −0.7 | −1.0 | −0.6 | 3.03 |
| 13.8 |  0.0 | −1.0 | **−2.3** | −2.2 | −1.6 | −1.4 | −1.4 | −0.8 | 3.92 |
| 19.4 | −0.1 | −1.5 | −1.6 | −0.6 | −1.0 | −0.9 | −0.5 | −0.1 | **4.41** |
| 24.9 | −0.1 | −0.9 |  0.3 |  2.4 |  2.8 |  2.6 |  2.2 |  1.0 | 3.99 |
| 30.5 | −0.1 |  0.8 |  3.0 |  5.1 | **6.4** |  6.2 |  4.8 |  2.0 | 3.38 |
| 36.0 |  0.1 |  2.3 |  4.7 |  4.9 |  4.8 |  4.4 |  3.5 |  1.5 | 3.01 |
| 41.5 |  0.5 |  2.2 |  2.8 |  1.0 | −0.7 | −1.0 | −0.4 |  0.1 | 3.43 |
| 47.1 |  0.4 |  0.4 | −1.1 | −3.2 | −4.3 | −3.9 | −2.5 | −0.9 | 3.94 |
| 52.6 |  0.0 | −0.9 | −2.6 | −3.8 | −3.9 | −3.3 | −2.4 | −1.0 | 3.18 |

(Full 32-lat, both-hemisphere ω/RH/heat/precip in the NetCDF. Extratropics 45–65°
show the mid-latitude storm-track ascent; the global-max ascent is polar/eddy,
ω≈−8.5 hPa/day at 85.8°, not an ITCZ feature.)

**Stronger-hemisphere (unsymmetrized) subtropical descent peaks at ω = +10.3
hPa/day at 30.5°, σ≈0.56–0.69** (the other hemisphere ~+4.4). Averaged over the
temporal/hemispheric variability the subtropical subsidence spans **~+6 to +11
hPa/day** — i.e. the "~13 hPa/day" of prior notes.

## 4. RH(|lat|, σ) — symmetrized, %

| \|lat\| | .188 | .312 | .438 | .562 | .688 |
|---:|---:|---:|---:|---:|---:|
|  2.8 | 75 | 74 | 76 | 76 | 76 |
| 13.8 | 80 | 75 | 75 | 76 | 76 |
| 19.4 | 73 | 68 | 71 | 75 | 75 |
| 24.9 | 62 | 56 | 63 | 71 | 74 |
| 30.5 | 52 | 48 | 53 | 64 | 73 |
| 36.0 | 46 | **49** | **45** | 56 | 72 |
| 41.5 | 43 | 54 | 48 | 56 | 72 |

Subtropical (25–40°) free-troposphere RH bottoms out at **~45%** (σ≈0.31–0.56) and
averages **~62%** over the 15–35° band at ~500 hPa. The deep-tropical column is
moist (74–80%). Subtropics are **dry**, consistent with the prior "~66% RH" note.

## 5. Headline numbers (symmetrized unless noted)

| Quantity | SpeedyWeather T21 aquaplanet |
|---|---|
| Peak tropical ascent ω | **−2.3 hPa/day at ±13.8°** (−3.4 in the stronger hemisphere at 19.4°), σ≈0.31–0.44 |
| Equatorial ω (mid-trop) | **−0.4 hPa/day** (quiescent) |
| Subtropical descent ω | **+6.4 hPa/day symmetrized at 30.5°; +10.3 stronger-hemisphere**, σ≈0.56–0.69 |
| ITCZ structure | **Split / off-equatorial** — ascent & precip maxima at ±14–19°, quiescent (weakly subsiding) equator, modest equatorial precip minimum (2.8 vs 4.4 mm/day) |
| Subtropical free-trop RH | **~45% (min, 30–40°) to ~62% (15–35° mean)** — dry |
| Precip peak | **4.4 mm/day at ±19°** (convective; large-scale ≈0) |
| **Is "~13 hPa/day" descent or ascent?** | **DESCENT.** Subtropical subsidence ≈ +6→+11 hPa/day; peak ascent only ≈ −2→−3.4 |

## 6. Direct comparison with aeros (T21 RCE aquaplanet)

| Feature | **aeros** | **SpeedyWeather** |
|---|---|---|
| ITCZ | Double, off-eq ±8–14°, quiescent equator | **Same** — split, off-eq ±14–19°, quiescent equator |
| Peak ascent ω | **−8.8** hPa/day (8–14°) | −2.3 to −3.4 hPa/day (14–19°) |
| Subtropical descent ω | **+2 to +3** hPa/day (19–30°) | **+6 to +11** hPa/day (30°) |
| Subtropical free-trop RH | **~94%** (moist) | **~45–62%** (dry) |

**Interpretation for the aeros diagnosis.** Both models produce the *same
qualitative* split/double-ITCZ with a quiescent equator — so aeros's double-ITCZ
is not anomalous. The gap is in the **overturning partition and the subtropics**:

- SpeedyWeather's tropical ascent is actually *weaker* than aeros's, but its
  **subtropical subsidence is ~3–5× stronger** (+6→+11 vs +2→+3 hPa/day) and
  *concentrated* at ~30°.
- That strong, persistent subsidence is exactly what **dries the subtropical free
  troposphere to ~45–62% RH**. aeros, with weak subsidence, cannot dry its
  subtropics and stays stuck near ~94%.

So the actionable difference is not the ITCZ shape but the **strength of the
subtropical descending branch** (and the resulting subtropical drying). aeros's
subsidence is too weak; SpeedyWeather's is strong and concentrated. The "~13
hPa/day" reference figure is the *descent*, and aeros is short of it by 4–5×.

---

*Sources: SpeedyWeather.jl v0.21.1 source (`~/.julia/packages/SpeedyWeather/i8kOF/`);
run scripts `mwm/C_omega/*.jl`; output `mwm/C_omega/output/speedy_omega_T21L8.nc`.
Sign of ω verified against SpeedyWeather's `vertical_velocity!` ("positive down")
and calibrated to physical hPa/day (instantaneous peaks ±90–160; 800-day
zonal+time means an order smaller, as expected for episodic convective ascent vs
steady subsidence).*
