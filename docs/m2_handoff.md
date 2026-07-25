# M2 handoff — where to pick up

Written at the end of the session that added radiation, a surface energy budget,
and a model-top sponge (M2.4a–f). `origin/main` is at the merge of
`m2/radiation`. `make all openmp=1` is clean; all **17** acceptance tests pass.
Full detail is in [`m2_results.md`](m2_results.md) §13; this is the orientation
and the next step.

## What landed this session (M2.4a–f)

The radiation **operators** are done and tested — this was the big remaining M2
piece:

- **Clear-sky longwave and shortwave** (`aeros_radiation`), a port of SESAM's
  broadband band kernel onto the resolved column (design.md §5's option; ecCKD
  behind the same `scheme` selector later). Validated offline in `test_radiation`
  against physical targets: OLR 256 W/m², CO₂-doubling forcing +2.8 W/m²,
  planetary albedo 0.117, exact flux-divergence identity. §13.1.
- **Insolation** — present-day daily-mean (declination + hour angle), validated
  by S₀/π and S₀/4. A stopgap for the Laskar `insol` forcing.
- **Surface energy budget** (`aeros_surface`), prescribed-SST bulk turbulent
  fluxes — the source that closes the atmosphere's budget once Held–Suarez is
  removed. Unit-tested (`test_surface`). §13.2.
- **Grid seam + cadence** (`aeros_radiation_apply`), heating cached and
  recomputed every 3 h, forward-split onto `wrk%dt_phys`.
- **Model-top sponge** (C1 Rayleigh + C2 Newtonian, off by default) in
  `aeros_timestep`, and **prescribed ozone** — the stability fixes of M2.4f.

Everything is off by default (`enabled = .FALSE.`, `sponge_on = .FALSE.`), so HS
and every existing test are bit-unchanged.

## The immediate next step: the coupled-RCE lowest-layer instability

**This is the blocker on validating anything against climate.** `test_rce` runs
the full stack (surface + radiation + convection + condensation, HS off, sponge
on) and is stable and active over a few hundred steps — but it is a short-run
check, NOT an equilibrium test, because a long run does not stay up.

The confirming long integration (see the scratch driver approach below) blows up
after **~1–2 model months** via a **grid-scale hot spot in the lowest layer
(level 12), growing to ~375 K**. The sequence of fixes already applied and the
failure they leave is in `m2_results.md` §13.3 — read it first. In short: the
model-top over-cooling is *solved* (forward-split heating + ozone + sponge); what
remains is a **surface-flux / moist-physics interaction in the boundary layer**.

- **Where to start.** Instrument which term drives level 12: surface sensible
  heat (`aeros_surface`), convective heating (`aeros_convection`, forward-split),
  or condensation (`aeros_condensation`, centered). The hot spot is confined to
  the lowest layer and grows slowly — it reads like a positive feedback or a
  grid-scale (spectral) noise mode there, not a single blow-up.
- **Hypotheses worth testing, cheapest first.** (a) dt=900 s + stronger `tau_diff`
  — separates numerical fragility from a physical feedback; if it survives, it is
  CFL/diffusion, if not, it is physical. (b) A surface-flux magnitude limiter or
  a minimum-`|U|` that is too low. (c) Whether the forward-split surface heating
  wants the *updated* (convectively adjusted) lowest-layer T rather than the
  start-of-step snapshot. (d) Boundary-layer vertical diffusion — deferred at
  M2.4c on the bet that convection alone would carry the surface fluxes up; the
  hot spot may be that bet failing (the fluxes pile up in level 12).
- **Do NOT** paper over it with a temperature clamp. It is a real coupling
  problem and the point of the RCE is that the climate emerges.

Once the long run is stable and reaches a bounded RCE, the moist line
(condensation + convection) can finally be tuned against something, and the ERA5
validation below becomes possible.

## Then: ERA5 ingestion + climate validation (blocked on data)

AR is preparing the data; the spec was given (Tier 1 unblocks radiation
validation, Tier 2 the moist stack). Two Copernicus datasets, monthly-mean
climatology, 2.5°, regridded to sigma on ingest:

- **pressure levels**: temperature, specific_humidity, ozone_mass_mixing_ratio,
  geopotential (Tier 1); cloud fraction + cloud water, relative_humidity, u, v
  (Tier 2).
- **single levels**: surface_pressure, skin_temperature, 2m_temperature,
  orography, land_sea_mask, forecast_albedo; TOA/surface radiative fluxes
  **including the clear-sky variants** (ttrc, tsrc, strc, ssrc, strdc, ssrdc) —
  those let us validate the clear-sky operator directly; total_precipitation,
  evaporation, sensible/latent heat flux (Tier 2).

The plan once it lands: regrid → validate clear-sky fluxes zonally against the
operator (this can be done column-by-column even before the RCE is stable) →
then turn on seasonal insolation and the real ozone field and tune.

## Other open M2 items (not blocking)

- **Radiation deferred pieces.** The cloudy-sky LW/SW branch (`lwr_clouds`,
  cloud albedo) waits on a cloud field. Ozone is a crude analytic lognormal —
  swap for the ERA5 zonal-mean field. Insolation is present-day perpetual —
  Laskar orbital forcing later. The slab-ocean surface energy balance (design.md
  §6.1) replaces the prescribed SST later.
- **Shallow convection is temperature-only** (M2.3, §12.2). Restore the
  moisture-redistributing shallow branch when there is something to tune its
  reference against — i.e. with a working RCE.
- **Manabe multi-layer adjustment** — lower priority since SBM is the default.
- **Water closure characterization** — check the ~1–3×10⁻³ FV-vs-spectral gap at
  T42 and over a long run.
- **Transport leftovers** — vertical van Leer, per-row polar sub-cycling.
- **Carried from M1** — `tau_diff`/∇⁶ retune, `SHqst_to_spat` 11→10, the Albedo
  16–128 thread sweep.

## Gotchas the next session will want

- **Physics seam, coupling paths.** HS forcing and condensation add to
  `wrk%dtdt` (grid, centered leapfrog). **Surface sensible heat, convection, and
  radiation** add forward increments to `wrk%dt_phys` (grid), applied to the n+1
  state in `aeros_timestep_step` step 6 — off the centered step, because their
  heating is large and vertically sharp. `dt_phys` is zeroed and applied whenever
  surface, convection OR radiation is on. Humidity (`qv_g`, incl. evaporation) is
  gridpoint and forward.
- **The sponge is a numerical device, not physics.** It lives in `aeros_timestep`
  next to the diffusion, `sponge_on = .FALSE.` by default. Its rates and `tref`
  are class fields; tune them there. Enabling it changes results, so keep it off
  for HS/benchmark reproducibility.
- **Reproducing the long run.** There is no committed long-RCE driver (it is not
  a unit test — it does not stay stable). Copy `tests/test_rce.f90`, raise
  `nstep`, add a trajectory print, and compile against the worktree lib with the
  flags `make -n rce` prints. Track `min/max T` per level and their level index —
  that is what localized the hot spot to level 12.
- **Worktree build needs the `fesm-utils` symlink.** A fresh worktree lacks it:
  `ln -sf /Users/alrobi001/models/fesm-utils <worktree>/fesm-utils`, and copy a
  generated `Makefile` in (gitignored; `config/Makefile` is the configme
  template). Use `git -C <worktree>` for git regardless.
- **Generated `Makefile` is gitignored.** M2.4 added the `radiation`, `surface`
  and `rce` test targets to the `config/Makefile` template *and* the object
  rules to `config/Makefile_aeros.mk`; a fresh checkout regenerates the root
  Makefile from those.
- **`timeout` is not on macOS** (it's `gtimeout`); don't wrap test runs in it.
- **Moist IC trap**: seed humidity from a realistic lapse-rate profile, not a
  fraction of `q_sat` against a warm column (`q_sat` explodes at low `p`).
