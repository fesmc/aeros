# M2 handoff — where to pick up (after the RCE resolution)

Supersedes the earlier M2.5a–d handoff. `main` is pushed (`origin/main`), `make
all openmp=1` is clean, and **all 18 acceptance tests pass**. This session
diagnosed and **resolved the M2 RCE runaway** (the long-standing blocker). Detail
below; the earlier handoff's content is folded in where still relevant.

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

## What landed this session (all on `main`, pushed)

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
2. **Do eddies now grow?** With `c_d>0` + `seed_asym>0`, does the drag-stabilized
   realistic jet develop baroclinic eddies (the `eddy_diag` metrics)? Tells us
   whether T21 supports meaningful eddy transport or T42 is needed.
3. **Clouds + the OLR −20 bias** (deferred M2.4 radiation piece): the cloudy-sky
   LW/SW branch fixes TOA balance physically and is needed for Tier-2 anyway.
4. **ERA5 moist-line (Tier 2):** the bounded-RCE blocker is now cleared, but it is
   **still blocked on the missing Tier-2 ERA5 fields** (cloud/RH/u/v were not in
   the delivered data). Tier 1 (clear-sky) is done (§14).

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

Cloudy-sky radiation (#3 above), ecCKD behind the `scheme` selector, Laskar
orbital insolation; CLIMBER-X ocean (slab is the interim); shallow convection is
temperature-only (§12.2); Manabe multi-layer adjustment; vertical van Leer;
per-row polar sub-cycling; `tau_diff`/∇⁶ retune.
