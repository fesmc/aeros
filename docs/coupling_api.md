# Aeros library coupling API — specification (A1)

Status: DRAFT for review. This document specifies the aeros public library API that CLIMBER-X's `src/atm_aeros/` adapter will call. It commits to what phases A–D will deliver and — importantly — what they will not. See [aeros ↔ CLIMBER-X coupling plan](../../.claude/projects/-Users-alrobi001-models-aeros/memory/aeros-climber-x-coupling-plan.md) for the surrounding context.

## 1. Goal

Aeros becomes usable in two shapes without duplicating code:

1. **Standalone**, driven by an executable driver (aeros_run, rce_long, held_suarez, mwm/D_realearth, etc.). Existing drivers keep their behavior bit-for-bit.
2. **Coupled**, called as a Fortran library by a host (initially CLIMBER-X's `atm_aeros` adapter), which hands aeros boundary fields once per coupling cycle (1 day) and reads back grid-mean and per-tile surface fluxes.

The interface below is what the coupled shape needs. The standalone drivers must be refactorable to call the same routines so both shapes exercise the same code paths.

## 2. Non-goals for phase A

- No changes to CLIMBER-X. Any adapter, regridder, or configme install lives on the `aeros-dev` branch of climber-x and is Phase B.
- No 5-tile decomposition. Aeros stays 2-tile (ocean/land) through Phase A. The import/export **types** below are shaped so the 5-tile expansion is a Phase C change to their contents, not a redesign.
- No 4-band surface radiation split, no radiation Jacobians, no external-SEB mode on land tiles, no O3 3D reader. All Phase C.
- No new physics; no calibration change; no algorithmic change. Everything Phase A produces must be bit-for-bit on the existing regression tests.

## 3. Library entry points

Six public routines, all in a new module `aeros_lib` (or extension of the existing `aeros` module — decision in A2). Prefix `aeros_lib_` used below for clarity; final name TBD.

```fortran
subroutine aeros_lib_init(ams, cfg_file, defaults_file, ic_source, ierr)
    type(aeros_class),  intent(out) :: ams
    character(len=*),   intent(in)  :: cfg_file          ! case namelist file
    character(len=*),   intent(in), optional :: defaults_file  ! e.g. input/aeros_defaults.nml
    character(len=*),   intent(in), optional :: ic_source      ! "resting" | "rce_warm" | "restart" | "coupler"
    integer,            intent(out) :: ierr
end subroutine
```

Init reads all namelists, allocates grid/SHT/state/timestep/physics, loads all `boundary_source = "file"` inputs, installs the initial condition per `ic_source` (default `"resting"` matches today's facade behavior), and reports readiness. **Does not** start the time loop.

```fortran
subroutine aeros_lib_step(ams, dt_couple, imports, exports, ierr)
    type(aeros_class),     intent(inout) :: ams
    real(dp),              intent(in)    :: dt_couple    ! seconds (e.g. 86400)
    type(aeros_imports_t), intent(in), optional  :: imports
    type(aeros_exports_t), intent(out), optional :: exports
    integer,               intent(out)   :: ierr
end subroutine
```

Advances the model by `dt_couple` seconds by sub-cycling `aeros_timestep_step` `nint(dt_couple / par%dt)` times. Consumes `imports` at the *start* of the cycle (writes into the internal boundary buffers) and populates `exports` at the *end* by reading grid-mean and tile-averaged accumulators. Both bundles are optional so standalone drivers can skip them entirely (bit-for-bit equivalence with today).

```fortran
subroutine aeros_lib_finalize(ams, ierr)
    type(aeros_class),  intent(inout) :: ams
    integer,            intent(out)   :: ierr
end subroutine

subroutine aeros_lib_restart_write(ams, filename, ierr)
    type(aeros_class),  intent(in)  :: ams
    character(len=*),   intent(in)  :: filename
    integer,            intent(out) :: ierr
end subroutine

subroutine aeros_lib_restart_read(ams, filename, ierr)
    type(aeros_class),  intent(inout) :: ams
    character(len=*),   intent(in)    :: filename
    integer,            intent(out)   :: ierr
end subroutine
```

Finalize deallocates in the same order as the current `aeros_end`. Restart wraps today's `aeros_timestep_write_restart` / `_read_restart` unchanged.

```fortran
subroutine aeros_lib_diagnose(ams, diag_out, ierr)
    type(aeros_class),  intent(in)  :: ams
    type(aeros_diag_t), intent(out) :: diag_out
    integer,            intent(out) :: ierr
end subroutine
```

Optional. Returns the diagnostic accumulator snapshot without stepping. Phase A ships a minimal `aeros_diag_t` — the rce_long-style monthly/annual accumulators that a coupled host cannot compute externally.

## 4. Error handling

Every public routine returns `ierr` (0 = success, non-zero = failure). Internal callers keep today's `error stop 1` semantics as a fallback where recovery is not possible (out-of-memory, missing config), but a coupled host must be able to survive a boundary-input mismatch, a restart-file mismatch, or a step-count anomaly without process death. **Phase A converts every `stop 1` / `error stop 1` on the public path to an `ierr` return**; internal `stop`s in physics kernels stay as-is for now (Phase C hardening).

## 5. State object

The existing `type(aeros_class)` ([src/aeros.f90:51](../src/aeros.f90)) stays as the single opaque handle. Fields today: `par, pool, grd, vgrid, now, ts, time, time_init, nstep`. Phase A additions:

- `bounds : aeros_boundary_t` — a container of every boundary provider (see §7).
- `diag_acc : aeros_diag_acc_t` — library-side accumulators for the diagnostics a coupled host cannot compute externally (grid-mean SLP, 700-hPa winds, column-mean T/RH/gams, monthly T2m). Phase A ships the minimum needed to answer "what does the coupler need each cycle?" — Phase C extends.

`ts%*` (physics state) remains unchanged in shape. Direct member writes to `ts%surf`, `ts%cnv`, `ts%rad`, `ts%cpr`, `ts%ocn`, `ts%land`, `ts%bd`, `ts%vd`, `ts%mst`, `ts%cnd` from rce_long move behind per-subsystem `aeros_<sub>_load(nml_file)` routines that `aeros_lib_init` calls in a fixed order. Namelist groups: `&surface, &convection, &condensation, &radiation, &clouds, &ocean, &land, &brewer_dobson, &vdiff, &moisture_transport, &insolation`. Existing scalar toggles in `&rce` split across those groups; a new `&physics` group holds the on/off master flags (`l_surf, l_cnv, l_cnd, l_rad, l_land, l_topography, l_seaice, l_cloud_prog, l_bd`, …) that today live in `&rce`.

Module-scope `_configure` calls (`aeros_cloud_configure`, `aeros_sw_aerosol_configure`, `aeros_ecckd_sw_configure`) move behind the `_load` routines so their state is set through the same call chain as `_init`, restart-visible, and survives a coupled `init → restart_read` flow.

## 6. The bundle types

Names, units, and physical meanings mirror CLIMBER-X's `cmn_to_atm` / `atm_to_cmn` so the adapter is a translation, not a rethink. Shapes are on the **aeros native grid** (`nlon, nlat`) — regridding to/from the CLIMBER-X 72×36 exchange grid is the adapter's job (fesm-utils conservative remap).

### 6.1 `aeros_imports_t` — what the coupler pushes in each cycle

Cycle-frequency fields (every `aeros_lib_step`):

```
sst           (nlon,nlat)      real(wp)  ! sea surface temperature [K]
sic           (nlon,nlat)      real(wp)  ! sea-ice fraction        [0-1]
t_ice_sfc     (nlon,nlat)      real(wp)  ! sea-ice surface T       [K]
t_skin_land   (nlon,nlat)      real(wp)  ! land skin T             [K]
q_skin_land   (nlon,nlat)      real(wp)  ! land skin q             [kg/kg]
albedo        (nlon,nlat)      real(wp)  ! broadband albedo (Phase A single-band; Phase C: 4-band VIS/NIR × dir/dif)
z0m           (nlon,nlat)      real(wp)  ! momentum roughness      [m]
sh            (nlon,nlat)      real(wp)  ! sensible heat flux from surface   [W m-2]
lh            (nlon,nlat)      real(wp)  ! latent   heat flux from surface   [W m-2]
evp           (nlon,nlat)      real(wp)  ! evaporation             [kg m-2 s-1]
lwu           (nlon,nlat)      real(wp)  ! upward LW from surface  [W m-2]
so4_load      (nlon,nlat)      real(wp)  ! sulfate load            [kg m-2]
```

Year-boundary fields (only re-read when `imports%year_boundary = .true.`):

```
co2_ppm       scalar           real(wp)
ch4_ppb       scalar           real(wp)   ! Phase A: accepted, may not be wired
n2o_ppb       scalar           real(wp)   ! Phase A: accepted, may not be wired
time_bp       scalar           real(wp)   ! years BP, drives insol orbital params
```

The `imports` argument is optional; when absent, aeros uses its own boundary providers (file-cached), preserving standalone behavior. When present, aeros writes the values into its `bounds` buffers via `external_buffer` providers (§7), and any provider still configured as `file` on the aeros side takes precedence over the corresponding import.

Phase A wires the subset that has aeros-side readers today: `sst`, `albedo` (single-band), `so4_load`, `co2_ppm`, `time_bp`. The rest are shaped in the type but ignored — Phase C wires them.

### 6.2 `aeros_exports_t` — what the coupler reads out each cycle

```
! Surface fluxes / near-surface state
ps            (nlon,nlat)      real(wp)  ! surface pressure        [Pa]
slp           (nlon,nlat)      real(wp)  ! sea-level pressure      [Pa]
t2m           (nlon,nlat)      real(wp)  ! 2 m T                   [K]
q2m           (nlon,nlat)      real(wp)  ! 2 m q                   [kg/kg]
u10m          (nlon,nlat)      real(wp)  ! 10 m u                  [m/s]
v10m          (nlon,nlat)      real(wp)  ! 10 m v                  [m/s]
wind10m       (nlon,nlat)      real(wp)  ! 10 m wind speed         [m/s]
taux          (nlon,nlat)      real(wp)  ! zonal      wind stress  [N m-2]
tauy          (nlon,nlat)      real(wp)  ! meridional wind stress  [N m-2]

! Radiation at surface
sw_dw_sur     (nlon,nlat)      real(wp)  ! downward SW    [W m-2]
lw_dw_sur     (nlon,nlat)      real(wp)  ! downward LW    [W m-2]
sw_net_sur    (nlon,nlat)      real(wp)  ! net SW         [W m-2]  (Phase A: convenience)

! Water
rain          (nlon,nlat)      real(wp)  ! liquid precipitation    [kg m-2 s-1]
snow          (nlon,nlat)      real(wp)  ! solid  precipitation    [kg m-2 s-1]
prc_total     (nlon,nlat)      real(wp)  ! rain + snow             [kg m-2 s-1]

! Grid-mean climate diagnostics (SESAM-style)
tam           (nlon,nlat)      real(wp)  ! column-mean T           [K]
ram           (nlon,nlat)      real(wp)  ! column-mean RH          [0-1]
gams          (nlon,nlat)      real(wp)  ! near-surface lapse rate [K m-1]
u700, v700    (nlon,nlat)      real(wp)  ! 700 hPa winds           [m/s]
cld           (nlon,nlat)      real(wp)  ! total cloud cover       [0-1]

! Cycle-mean TOA (nice-to-have, Phase A optional)
olr_toa       (nlon,nlat)      real(wp)
sw_toa_net    (nlon,nlat)      real(wp)
sw_toa_dw     (nlon,nlat)      real(wp)
```

Each field is a **cycle-mean** accumulated during the sub-cycle loop, not an instantaneous snapshot. `aeros_lib_step` reads accumulators from `ams%diag_acc` at the end of the cycle and writes them into `exports`, then zeroes the accumulators.

Phase A wires what aeros can already compute cheaply: `ps, slp, t2m, q2m, u10m, v10m, taux, tauy, sw_dw_sur, lw_dw_sur, rain, snow, tam, ram, gams, cld` and the TOA trio. Phase C adds `u700, v700` (needs vertical interp), any missing tile-decomposed fields, and radiation Jacobians.

### 6.3 `aeros_diag_t` — optional detailed diagnostics for the host

Extended snapshot for a host that wants more than the cycle-mean export bundle — zonal-mean climatology, spectral energy, energy/water budgets. Phase A ships an empty shell (name-reserved). Phase C fills it with what today's `rce_long` accumulators produce.

## 7. Boundary provider abstraction

New module `aeros_boundary` (working name). One provider per boundary field, discriminated union:

```fortran
integer, parameter :: BC_CONSTANT = 0
integer, parameter :: BC_FILE     = 1
integer, parameter :: BC_EXTERNAL = 2

type :: aeros_bc_scalar_t
    integer  :: kind
    real(wp) :: value                     ! CONSTANT
    character(len=256) :: file            ! FILE (year-slice)
    integer  :: year                      !
    real(wp) :: current                   ! last-read value
end type

type :: aeros_bc_2d_t
    integer  :: kind
    real(wp), allocatable :: buffer(:,:)  ! (nlon,nlat), always allocated after init
    character(len=256) :: file            ! FILE
    character(len=64)  :: varname         !
    logical  :: is_monthly                ! FILE: 12-month climatology
    real(wp), allocatable :: clim(:,:,:)  ! FILE monthly: (nlon,nlat,12)
    integer  :: itime                     ! FILE single-slice: 1-based
end type
```

`aeros_boundary_t` bundles one field per known boundary (`sst, sic, t_ice_sfc, t_skin_land, q_skin_land, albedo, z0m, phis, lsm, dust, so4, o3, co2, ch4, n2o, time_bp, tsi`). Each is configured at init from its namelist (§5) as one of:

- `BC_CONSTANT`: scalar; `buffer` is a uniform fill.
- `BC_FILE`: reads via `aeros_bcinput_read_field` at init (one-shot or 12-month clim). Daily interp for `is_monthly` fields.
- `BC_EXTERNAL`: `buffer` is written by the caller (via `aeros_lib_step(imports=…)`); no I/O.

Every physics module that today reaches into a driver-owned field (rce_long's `sst_clim`, aeros_land's `land%albedo_clim`, aeros_radiation's `rad%dust_clim`, aeros_topography's `phis_full`) is refactored to consume from `ams%bounds%<field>%buffer` instead. Bit-for-bit preserved by defaulting every namelist to `BC_FILE` with today's paths.

Timing (each `aeros_lib_step` cycle):

1. If `imports` present, copy provided fields into their `BC_EXTERNAL` buffers.
2. For `BC_FILE` monthly fields, refresh `buffer` from `clim` by day-of-year interp (same routines as today).
3. For `BC_CONSTANT`, no-op.
4. Physics sees a fresh `buffer(:,:)`.

Fields with no aeros reader today (SIC, T_skin, q_skin, z0m map, O3(z), CH4/N2O) get their type declared but the physics consumer is a Phase C wire-up.

## 8. Configuration

The A2 refactor introduces a per-subsystem `_load(nml_file)` for each physics module (see §5) that populates its config struct from a namelist group. `aeros_lib_init` calls them in a fixed order:

```
aeros_par_load          -> par           (&aeros)
aeros_vgrid_load        -> vgrid         (&vgrid)
aeros_correction_load   -> correction    (&correction)
aeros_physics_load      -> physics flags (&physics)
aeros_surface_load      -> ts%surf       (&surface)
aeros_moist_load        -> ts%mst        (&moisture_transport)
aeros_convection_load   -> ts%cnv        (&convection)
aeros_condensation_load -> ts%cnd        (&condensation)
aeros_ocean_load        -> ts%ocn        (&ocean)
aeros_land_load         -> ts%land       (&land)
aeros_vdiff_load        -> ts%vd         (&vdiff)
aeros_bd_load           -> ts%bd         (&brewer_dobson)
aeros_cloud_load        -> ts%cpr        (&clouds)
aeros_radiation_load    -> ts%rad        (&radiation)  ! also calls aeros_sw_aerosol_configure, ecCKD configure
aeros_insol_load        -> insol config  (&insolation)
aeros_boundary_load     -> ams%bounds    (per §7, from &boundaries + subsystem file paths)
```

Standalone drivers become thin wrappers: `aeros_lib_init(ams, "input/aeros_realearth_T21L6.nml")` replaces the ~500-line rce_long config-write block. rce_long's driver-specific bits (`&rce` group's dump paths, spin-up thresholds, uniform-insol probe) stay in rce_long as before, but every physics knob it currently sets moves into the appropriate `&<subsystem>` group.

## 9. Restart contract

Phase A does not change the restart file format. `aeros_lib_restart_{write,read}` wraps today's `aeros_timestep_write_restart` / `_read_restart` verbatim.

Two additions in Phase A:

1. **Boundary buffers with `BC_EXTERNAL` kind** are written to the restart file, because after a coupled restart the coupler has not yet pushed new values and physics needs a valid starting field. `BC_FILE` and `BC_CONSTANT` buffers are not written (they're re-derived).
2. **Module-scope configure state** (cloud, sw_aerosol, ecCKD) becomes restart-visible by being carried inside the corresponding `_load`-populated struct, not in file-scope module state. This closes the pre-existing hazard where a restart could resume with different aerosol/cloud settings than the original run.

## 10. Calendar

Runtime `calendar = "julian" | "360d"` in the `&aeros` namelist group (default `"julian"` — Julian remains bit-for-bit). One `days_year` variable in `aeros_defs` set at `aeros_par_load`; every `DAYS_YR` / `DAY_YEAR` / `year_seconds` `parameter` today (8 files, ~50–150 LOC per the calendar survey) rewires to consume it. `insol_init` already accepts `day_year ∈ {360, 365, 366}` so the insolation library needs only its call-site plumbed. Monthly interpolation is already `DAYS_YR/12`-derived, so 360-day → 30-day months for free. Pre-existing 0.25 d/yr Julian phase drift (365 vs 365.25 mismatch) is documented but not fixed here.

## 11. Grid

Phase A stays on aeros's native T-Gaussian grid. `par%trunc` is a runtime variable that already accepts arbitrary values (validated `trunc >= 1` in [aeros_defs.f90:362](../src/aeros_defs.f90)); T18 will need `nlon, nlat` matched (56, 28 is the natural choice — dealiased Nlon ≥ 3T + 1 = 55, rounded to a friendly FFT size). A5 exercises this in `mwm/D_realearth_T18L6` as a shakedown for hardcoded T21 assumptions.


### 11.1 Grid dimensions after `aeros_lib_init`

`ams%par%nlon` and `ams%par%nlat` retain their **namelist request** value after init, not the realised grid dimensions. The default request is `-1` (auto-infer from `trunc`), so a host that reads `%par%nlon` after init sees `-1` unless the namelist pinned an explicit value. The realised dimensions live on the grid object:

```fortran
call aeros_lib_init(ams, "input/aeros_realearth_T21L6.nml", ierr=ierr)
allocate(imp%sst(ams%grd%nlon, ams%grd%nlat))   ! CORRECT
!                    ^^^^^^^^ not ams%par%nlon
```

Every `aeros_imports_t` / `aeros_exports_t` field on the aeros native grid must be sized against `ams%grd%*`. This was the source of the Wave 4d bug where `aeros_lib_step`'s auto-allocated monthly bin used `%par%nlon = -1` and silently corrupted memory; the library itself now consistently uses `%grd%*`, but hosts must too.

## 12. What Phase A explicitly ships

| Piece | Concrete deliverable |
|-------|----------------------|
| Library API | `aeros_lib_init/step/finalize/restart_{read,write}/diagnose` — public routines added to the existing `aeros` module ([src/aeros.f90](../src/aeros.f90)), not a new module. |
| State handle | `aeros_class` unchanged shape, with `ams%bounds` and `ams%diag_acc` additions |
| Bundle types | `aeros_imports_t`, `aeros_exports_t`, `aeros_diag_t` (fully populated — see below) |
| Provider abstraction | `aeros_boundary_t`, `aeros_bc_2d_t`, `aeros_bc_scalar_t` — `BC_FILE` and `BC_EXTERNAL` wired for the fields listed in §6.1 |
| Per-subsystem loaders | `aeros_<sub>_load(nml_file)` for the subsystems in §8, consumed from one many-group namelist file (as today) |
| Diagnostic accumulators | Full — the rce_long monthly means, second-half time-mean, energy/water budgets, spectra all move library-side into `ams%diag_acc`, populating `aeros_diag_t`. A coupled `aeros_lib_step(1 day)` produces the same monthly climatology and budget bookkeeping that standalone rce_long does today. |
| Error surface | `ierr` returns on all public routines, `stop 1` retained internally |
| Calendar | `calendar = "julian" \| "360d"` runtime option |
| T18 harness | `mwm/D_realearth_T18L6` (not bit-for-bit) |
| Drivers | aeros_run (kept and updated, canonical minimal example), held_suarez, mwm/D_realearth, and rce_long all refactored to call `aeros_lib_*`. Physics knobs in `&rce` migrate to `&surface/&convection/…` **breaking** — `&rce` stops accepting them the day the loaders land; namelist files (`aeros_defaults.nml`, mwm/D_realearth/*.nml, per-subsystem overrides) migrated in the same commit series. |
| Regression | `tests/test_restart.f90`, `test_ncycle.f90`, `test_moist_run.f90`, `test_held_suarez.f90`, `test_timestep.f90`, `test_rce.f90` all bit-for-bit green with Julian calendar |

## 13. What Phase A explicitly does not ship

- 5-tile decomposition (2 tiles today, stays 2 tiles in A)
- 4-band VIS/NIR × dir/dif surface radiation split
- Radiation Jacobians (`dswd/dalb`, `dswd/dz`)
- External-SEB mode on land tiles (aeros closes its own SEB in every A-shipped run)
- O3 3D-in-space reader (analytic profile stays)
- CH4/N2O/CFC wiring (types declared, not consumed)
- Roughness map (z0 stays scalar)
- Any CLIMBER-X code
- Any regridder
- Any new physics or calibration change

## 14. Design decisions (locked)

1. Public routines added to the existing `aeros` module, not a new `aeros_lib`.
2. One namelist file with many groups (as today).
3. Full diagnostic accumulators land in Phase A — no coupled-run monthly-diagnostics gap.
4. Breaking split of `&rce` physics knobs into `&surface/&convection/…` — no compatibility period.
5. rce_long refactored to call `aeros_lib_*` — one code path.
6. aeros_run kept and updated as the canonical minimal example.
