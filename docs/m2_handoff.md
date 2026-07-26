# M2 handoff — where to pick up (after M2.5a–d)

Supersedes the earlier M2.4 handoff. `main` is at the M2.5 series; `make all
openmp=1` is clean and **all 18 acceptance tests pass**. Full detail is in
[`m2_results.md`](m2_results.md) §14–16; this is the orientation, the one open
blocker, and the recommended choices for it.

## What landed this session (M2.5a–d)

- **§14 — ERA5 clear-sky radiation validation (M2.5a).** The clear-sky LW/SW
  operators driven on ERA5's native pressure-level columns (no regrid) vs ERA5
  clear-sky fluxes, lat×lon maps. Patterns match; biases modest and
  interpretable — **OLR is ~16 W/m² too low (the scheme is a touch too opaque)**,
  surface SW ~16 too high (no aerosol). `drivers/validate_era5.f90`,
  `scripts/plot_era5_validation.jl`. Tier 1 only; Tier 2 (cloud/RH/u/v) is not in
  the delivered data.
- **§15 — the RCE lowest-layer instability, diagnosed and mitigated (M2.5b–c).**
  It is **perfectly zonally symmetric** — an axisymmetric column instability, not
  a grid-scale hot spot. Two fixes landed: `aeros_vdiff` (boundary-layer vertical
  diffusion, **implicit on the n+1 state** — forward-splitting a diffusion term
  on a leapfrog is unconditionally unstable) and condensation moved to the
  forward-split path (its coupled-RCE heating is no longer "small and smooth").
  Together they push the blow-up from day 34 to day ~200.
- **§16 — slab ocean (M2.5d) + the defaults-file paradigm.** `aeros_ocean` owns
  the SST (prescribed | slab surface-energy-balance), the plug point for a real
  ocean. Parameter files are now overrides over `input/aeros_defaults.nml`.

Everything new is **off by default** (`vdiff`, `surface`, `radiation` enabled
`.FALSE.`; ocean `mode = prescribed`), so HS and every prior run/test are
bit-unchanged.

## The one open blocker: the subtropical axisymmetric column instability

A bounded, equilibrated RCE is the prerequisite for the ERA5 moist-line
validation, and it is **not yet reached**. The fast blow-up is fixed; what
remains is a slower runaway that survives everything tried so far. The
characterization is sharp (all via `drivers/rce_long.f90`, which has the
instrumentation — per-level min/max T with location, zonal-mean-vs-departure
split, per-term heating, TOA/surface energy balance, and namelist knobs for
every physics toggle, `vdiff_k0`, `albedo`, `ocean_mode/depth`, `seed_asym`):

- The runaway is a **zonal-mean lowest-layer temperature at a subtropical
  latitude** climbing to 360–380 K (`max |T − zonalmean| = 0.00 K` — purely
  axisymmetric with the default m=0 start).
- It is **robust** to: vdiff strength, condensation coupling, a cloud-proxy
  albedo that balances TOA to ~0, a slab ocean that makes the surface flux
  self-limiting, and an eddy seed that breaks the m=0 symmetry. Each delays or
  changes it; none removes it (slab+albedo → NaN day 176; slab+albedo+eddies →
  day 88).
- So it is **not** a coupling/numerics bug and **not** the energy imbalance — it
  is a property of the moist axisymmetric column at T21L12 with the current
  convection scheme.

### Recommended choices, cheapest first

1. **Diagnose the mechanism before fixing it (½ day).** Run `rce_long` to the
   onset and read, at the hot latitude, the per-term heating (surface / vdiff /
   convection / condensation / radiation) and the SBM band it selects. The
   question to settle: is the subtropical descending column getting a spurious
   convective/condensational heating, or is the lowest layer simply not being
   ventilated? This tells you which of the below is the real fix rather than
   guessing.
2. **Resolution (1 day).** Rerun the same config at **T42L20**. If the
   instability weakens or vanishes, it is under-resolution of the subtropical
   subsidence/boundary layer and the path is simply to run at higher resolution
   (design.md's production target anyway). If it persists at T42L20, it is the
   scheme — go to 3/4.
3. **Clouds (the real physics, ~1 week) — recommended regardless.** The
   cloudy-sky LW/SW branch (`lwr_clouds`, cloud albedo) is the deferred M2.4
   radiation piece. It is worth doing next on its own merits: it fixes the TOA
   balance *physically* (albedo from clouds, ~0.3, rather than the cloud-proxy
   knob), and cloud radiative effects in the subtropics (low-cloud cooling of the
   descending branch) are exactly the missing feedback that could stabilise the
   hot latitude. This is my pick for the substantive next step — it is needed for
   climate validation anyway, and it plausibly addresses both the global balance
   and the local runaway at once. Needs a cloud field (a diagnostic cloud
   fraction/water from RH is enough to start).
4. **Convection at the descending branch (if 1 points there).** If the diagnosis
   shows SBM misbehaving in the dry subsiding column, revisit its trigger/closure
   there (the shallow branch is temperature-only, §12.2) rather than adding more
   diffusion.

Do **not** reach for a temperature clamp or stronger sponge — the point of the
RCE is that the climate emerges.

## Then: ERA5 climate validation (Tier 1 done, Tier 2 blocked)

- **Clear-sky radiation (Tier 1): done** (§14), and it flagged the OLR −16 W/m²
  bias — worth a look when clouds land, as the two interact.
- **Moist line (Tier 2): blocked twice** — on a bounded RCE (above) and on the
  Tier-2 fields (cloud fraction/water, RH, u, v), which were **not** in the
  delivered ERA5 data. When both clear, turn on seasonal insolation + the real
  ozone field and tune the moist stack against ERA5.

## Other open M2 items (not blocking)

- **Radiation deferred:** cloudy-sky branch (see #3 above), ecCKD behind the
  `scheme` selector, Laskar orbital insolation (present-day perpetual now).
- **Ocean:** the slab is the interim; CLIMBER-X's ocean is the eventual
  component (design.md §6.1) — plugs in where the slab is.
- **Shallow convection is temperature-only** (§12.2); restore the
  moisture-redistributing branch once there is a working RCE to tune against.
- Carried from before: Manabe multi-layer adjustment, water-closure
  characterization at T42 over a long run, vertical van Leer, per-row polar
  sub-cycling, `tau_diff`/∇⁶ retune, the Albedo thread sweep.

## Gotchas the next session will want

- **The RCE driver is `drivers/rce_long.f90`** (committed, unlike the old scratch
  driver). Namelist-configured; `logs/rce_*.nml` have working examples. It is a
  diagnostic tool, not an acceptance test — it does not stay stable yet.
- **Physics coupling paths.** Surface SH, convection, condensation, radiation and
  vdiff's temperature change all ride the **forward-split** `wrk%dt_phys` (applied
  to n+1 after the dynamics step). vdiff is the exception in HOW: it is a
  diffusion operator, so it is applied **implicitly on the n+1 state** at step 3c
  next to the horizontal diffusion/sponge — never forward-split. Humidity
  (incl. evaporation, condensation drying) is gridpoint-forward. The ocean SST
  steps after surface+radiation from their surface fluxes.
- **Defaults paradigm.** Case `.nml` files carry only overrides;
  `input/aeros_defaults.nml` is the full schema. `nml_read` takes an optional
  `defaults_file` (threaded through the load chain, not module state — `nml` is
  stateless/shared). Absent it, legacy "case file must be complete" behaviour
  holds, which is why the `_init`-only tests are unaffected. Drivers require the
  defaults file.
- **A bounded RCE needs a cloud-proxy albedo with the current clear-sky scheme.**
  Ocean albedo 0.06 gives a +90 W/m² TOA imbalance; ~0.4 balances it. This goes
  away once clouds provide the albedo.
- **`timeout` is not on macOS** (it is `gtimeout`); don't wrap runs in it.
- **Moist IC trap:** seed humidity from a realistic lapse-rate profile, not a
  fraction of `q_sat` against a warm column aloft.
