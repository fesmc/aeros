# M2 handoff — where to pick up

## ►► NEXT SESSION: topography — the first real boundary-condition forcing

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
