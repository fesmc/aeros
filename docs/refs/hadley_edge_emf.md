# Cell-edge diagnosis: aeros has no eddy momentum flux, SpeedyWeather does

The residual subtropical moist bias (upper troposphere pinned ~94% RH vs
SpeedyWeather/ERA5 ~45–62%) was pinned in `descent_concentration_scope.md` to
aeros's Hadley **descent not concentrating at ~30°**. This note reports the
diagnostic that identifies *why*, at the same T21 aquaplanet footing.

## The diagnostic (landed this branch)

- **aeros** (`drivers/rce_long.f90`): time-mean zonal-mean jet `ubar`, MMC `vbar`,
  and **eddy momentum flux `uvpr` = [u*v*]** (zonal-mean of the deviations from
  the instantaneous zonal mean), accumulated over the equilibrated 2nd half —
  same window as the `omega` time-mean. Written to the `l_diag` dump.
- **SpeedyWeather** (`mwm/C_omega/omega_structure.jl`): the same three fields
  (`uwind`, `vwind`, `uvpr`), with `[u*v*] = [uv] − [u][v]` (stationary+transient
  eddy flux) from per-step point accumulation.
- **Comparison** (`scripts/hadley_edge_compare.jl`): loads both, builds the mass
  streamfunction Ψ from `[vbar]`, the EMF convergence, and a 4-row figure
  (`docs/figures/hadley_edge_compare.png`).

## Result — same T21 aquaplanet

| | aeros (`rce_allfix`, seed 0) | SpeedyWeather |
|---|---|---|
| Eddy momentum flux \|[u*v*]\| | **≡ 0** (axisymmetric) | **up to 3.9 m²/s²** |
| Subtropical jet | 14 m/s, mid-level (σ≈0.31), smeared | deep 25–30 m/s free-trop (44 at σ=0.06), edge at 30° |
| Hadley descent | +1.9 hPa/day, spread **25→69°**, peak 36° | +6.4 hPa/day, tight **25–36°**, peak **30.5°** |
| Subtropical free-trop RH | ~94% | 62% |

SW converges `[u*v*]` into the midlatitude jet (~3–4 m²/s² at 35–58°) and diverges
it out of the subtropics (~−0.29 m/s/day at 15–35°). That eddy-driven termination
is what concentrates the descent at 30° and dries the subtropics. **aeros has none
of it** — its descending branch spreads diffusely poleward, so nothing dries the
subtropical free troposphere below the `rh_crit` ceiling and it sits near 94%.

## Interpretation

1. **This refutes the §25–26 framing** that weak eddies are "the expected behaviour
   of a perpetual-annual-mean aquaplanet, not a tuning target." SpeedyWeather is
   the *same* perpetual T21 aquaplanet and eddies vigorously. aeros's eddy-free
   state is an aeros deficiency, not a fair-yardstick artifact. Supersedes
   `[[weak-hadley-forward-split-coupling]]`.

2. **aeros cannot sustain an eddying state — it blows up at the model top.**
   With `seed_asym=0` the flow stays perfectly axisymmetric (EMF≡0, survives). With
   a seed, baroclinic eddies begin to grow at the subtropical latitudes (eddyKE
   0 → ~0.2 m²/s² by day ~25) **then trigger a NaN at level 1 (model top)**:
   `seed_asym=0.05` dies day 27, `seed_asym=1.0` dies day 3. This is the same
   model-top instability flagged across the handoff (§29 slab runaway, land×cloud
   jet, topography blow-up) and repeatedly deferred to "a proper gravity-wave-drag
   / higher-top fix."

**So the residual moist bias traces to the model-top instability:** it forbids the
eddying state → no eddy momentum flux → no subtropical jet → no concentrated
descent → undried subtropics. The durable fix is model-top stabilization
(gravity-wave drag / higher top / deeper sponge), so aeros can hold eddies.

## Reproduce

```
# aeros (unseeded, eddy-free reference):
./libaeros/bin/rce_long.x logs/rce_allfix.nml logs/rce_emf_dump.nc
# SpeedyWeather (300 d spinup + 400 d mean):
julia --project=mwm/C_omega mwm/C_omega/omega_structure.jl 300 400
# comparison + figure:
julia scripts/hadley_edge_compare.jl logs/rce_emf_dump.nc \
      mwm/C_omega/output/speedy_omega_T21L8.nc docs/figures/hadley_edge_compare.png
```
