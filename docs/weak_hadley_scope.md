# Weak Hadley overturning — couple diabatic heating to the semi-implicit dynamics

**Status:** diagnosed, not started. For a fresh session.
**Goal:** fix aeros's Hadley/divergent circulation, which is **~10× too weak**, by coupling
the diabatic heating to the dynamics the standard way (into the thermodynamic tendency,
*ahead of* the semi-implicit solve) instead of the current forward-split increment applied
after the whole dynamics step. This is the root of the residual free-tropospheric moist bias
(and very likely the long-standing weak-jet deficit, §25–26).

**Chosen approach (AR, this session): approach 1 — move ALL diabatic heating into the
tendency (full standard coupling). "Do it well."** Not the radiation-only half-measure.

---

## 0. How we got here (the diagnosis chain)

The near-saturated free-troposphere moist bias (overcast, RH 85–100% vs ERA5 ~45%) was
attacked in two stages this session:

1. **Column physics — DONE, on `main`** (commit `a5395d6`). Diagnosed against SpeedyWeather
   (same Simplified Betts–Miller scheme) that aeros had two saturating sinks and no drying
   sink. Fixed:
   - **vertical diffusion** now diffuses **dry static energy** (not T → was driving the BL
     isothermal/over-stable) over a **Richardson-number diagnosed depth** (Frierson 2006
     K-profile, Ri_c=10; SpeedyWeather's `BulkRichardsonDiffusion`). Legacy fixed-depth kept
     (`vdiff_richardson=.false.`). `src/physics/aeros_vdiff.f90`.
   - **condensation** defaults to **rh_crit=0.95** (was 1.0, pinned grid mean at 100%) and
     **reevaporates** falling precip into sub-saturated layers (SpeedyWeather's
     `ImplicitCondensation`). `src/physics/aeros_condensation.f90`.
   - **convection default → `sbm_frierson`** (surface-to-LZB depth, integrated Pt/Pq trigger
     — the scheme SpeedyWeather runs; stays active where the strictly-buoyant band went
     dormant). `src/physics/aeros_convection.f90`.
   - **Result** (rotating RCE T21, `logs/rce_allfix.nml`, cond_rh_crit=0.95): cover
     0.999 → **0.70** (ERA5 0.63), TOA net −18.6 → **+3.8** W/m², tropics now moist-adiabatic
     and dried to **~73%**, convection active, no instability. 25 tests pass (new
     reevaporation + sbm_frierson cases).

2. **The residual, and the real cause — THIS TASK.** Free-trop RH is still **~86%** (global
   mean), pinned near 94–95% in the **non-convecting subtropics/extratropics**, because the
   resolved Hadley subsidence is only **~1 hPa/day**. Established, with evidence:
   - **Not column physics** (fixed above; the tropics dried correctly).
   - **Not resolution.** aeros T42 (`logs/rce_allfix_t42.nml`) is *the same* — 87% RH,
     ~1 hPa/day subsidence. And **SpeedyWeather at T21** (agent-run, matched aquaplanet) is
     **dry: 66% free-trop RH, ~13 hPa/day subtropical subsidence, coherent subtropical
     minimum ~44–50%** — identical to its T42. So at equal T21, SpeedyWeather has a vigorous
     Hadley cell and we do not.
   - **Not damping / not forcing.** Sponge off → no change; weaker hyperdiffusion → no
     change; weaker RAW filter → NaNs (no help). SST gradient is *steeper* than
     SpeedyWeather's in the subtropics (13.5 K vs 7.3 K over 0–30°).
   - **The circulation is genuinely weak** (not a diagnostic artifact): zonal-mean
     upper-branch meridional wind max **|v| ≈ 0.57 m/s** (should be ~2–4), subtropical jet
     |u| ~15 m/s (ERA5 ~30). Confirmed by the `omega` dump (see tools below).

| | aeros T21 | aeros T42 | **SpeedyWeather T21** | SpeedyWeather T42 |
|---|---|---|---|---|
| free-trop RH | 86% | 87% | **66%** | 66% |
| subtropical ω | ~1 | ~1 | **~13 hPa/day** | ~13 |
| cover | 0.70 | 0.71 | ~0.63 | 0.63 |

---

## 1. Root cause (confirmed structurally in the code)

`wrk%dtdt` — the temperature tendency that feeds the semi-implicit solve
(`aeros_tendency.f90:581`, transformed to `tnd%temp` at `:463`) — carries **only the
dynamical heating**: horizontal advection − vertical advection + the adiabatic `κ T ω/p`.
**No diabatic term is in it.**

Every diabatic term — **surface** sensible heat, **convection**, **condensation**, and
**radiation** — is accumulated on the grid in `wrk%dt_phys` [K] (`aeros_timestep.f90:836–933`)
and applied as a **forward-split increment onto the n+1 temperature at step 6**
(`apply_phys_heating`, `aeros_timestep.f90:1021–1033`) — **after** the semi-implicit advance
(step 2), **after** horizontal diffusion + sponge + vdiff (step 3), and **after** the RAW
time filter (step 4).

**Why this starves the Hadley cell.** The Hadley/divergent circulation is the dynamics'
*response* to diabatic heating: heating → warms the column → raises the geopotential →
∇²Φ drives the divergence tendency → continuity gives ω (the overturning). In a standard
semi-implicit spectral model (IFS, SPEEDY, SpeedyWeather) the diabatic heating is added to
the thermodynamic tendency **before** the semi-implicit solve, so the gravity-wave adjustment
generates the divergent response *within the step*. aeros runs the dynamics on a purely
adiabatic tendency and bolts the heating on afterward, so the divergent response is (a)
lagged a full step and (b) fighting the RAW filter, which is turned up high (`eps_filter=0.06`;
0.02 NaNs) precisely to kill the computational mode this split excites. The circulation the
heating should drive is being damped out.

---

## 2. The task

Move the diabatic heating into `tnd%temp` **before** `aeros_semiimp_step` (and remove the
step-6 forward-split), so the divergent circulation is generated by the semi-implicit dynamics
in response to the heating — the standard coupling.

### Key facts for the implementation

- **Ordering is already favourable.** Physics (`dt_phys` accumulation) runs at lines 836–933,
  *before* `aeros_tendency_spectral` at :962 and the semi-implicit step at :981. So `dt_phys`
  is available to inject into the tendency with no reordering.
- **Grid → spectral.** `dt_phys` is a gridpoint [K] field; `tnd%temp` is spectral [K/s].
  Transform `dt_phys` (as a rate) per level via `aeros_sht_analysis` and add to `tnd%temp`.
- **The leapfrog factor is the crux.** The forward-split applies exactly `dt_phys` [K] per
  step to n+1. The tendency is advanced over `h = 2 dt` (leapfrog). So a heating *rate* added
  to `tnd%temp` contributes `h·rate = 2 dt·rate` per step — to reproduce `dt_phys` you'd add
  `dt_phys/(2 dt)`, but that is a **centered** contribution across n−1→n+1, which is exactly
  the treatment that excited the computational mode before. Understand the semi-implicit /
  leapfrog / RAW discretization (`aeros_semiimp.f90`, `raw_filter` in `aeros_timestep.f90`)
  and choose the coupling deliberately.

### The two real hazards (this is the "do it well" part)

1. **Computational-mode instability.** Forward-split was adopted because centered convective
   heating — large and sign-alternating in the vertical — NaN'd (m2_results §12.1). Full
   coupling must stay stable. Levers: the semi-implicit solve + RAW filter may now hold it
   (α=0.53, eps=0.06); if not, consider (a) the Williams/RAW filter tuning, (b) a gentler
   convective relaxation (larger `conv_tau`), (c) whether only convection needs special care
   while surface/condensation/radiation couple cleanly. Adaptive hyperdiffusion (`diff_adapt`,
   `diff_taper`, both available, default off) is a top-stability net if needed.
2. **Condensation heating/drying pairing.** Condensation removes `dqc` from `qv_g` (gridpoint,
   forward) **and** adds `L/cp·dqc` to `dt_phys` — paired from the *same* `dqc`, so the column
   MSE is conserved by construction (`test_condensation` pins this as an equality). Humidity is
   **off-spectral** in aeros (gridpoint, advected by `aeros_transport`), so it cannot go
   through the spectral semi-implicit. If the heating moves to the spectral centered tendency
   while the drying stays gridpoint-forward, the two are on different discretizations and the
   exact MSE pairing breaks. Decide how to keep the energy/water budget closed (accept a small
   discretization mismatch and verify it's bounded? a correction term? keep condensation
   heating forward-split and only move the others?). This is a genuine design decision — get it
   right, don't paper over it.

---

## 3. Acceptance

1. **Hadley cell strengthens toward SpeedyWeather**: subtropical ω from ~1 toward ~10 hPa/day,
   upper-branch |v| from ~0.6 toward ~2–3 m/s (use the `omega`/`v` dump).
2. **Free-trop RH drops** in the non-convecting subtropics/extratropics from ~94% toward
   SpeedyWeather's ~66% global mean; a coherent subtropical minimum appears.
3. **Cover stays ~0.63 and TOA stays near balance** (don't regress the column-physics win).
4. **Stable** (no NaN) over a 100-day RCE, ideally with the land case too.
5. **All 25 acceptance tests still pass**; energy/water budgets stay closed (watch
   `test_condensation`'s MSE equality if you touch the condensation pairing).
6. Bonus check: does the **jet** also strengthen toward ERA5 (~30 m/s)? If so, this fix closes
   two biases at once, and §25–26's "inherent aquaplanet limitation" conclusion should be
   revised (SpeedyWeather, also an aquaplanet, has a strong cell at T21).

---

## 4. Tools & references

- **`omega` diagnostic (landed, commit `c8120b8`).** `wrk%omega` [Pa/s] stored in the
  tendency under `l_diag`; the RCE dump writes the **time-mean** (run's 2nd half) zonal-mean
  `omega` [hPa/day, >0 = subsidence] to the RH NetCDF, alongside `v`, `u`, `rh`, `cf`, `t`.
  This is the primary instrument for the Hadley strength.
- **Configs:** `logs/rce_allfix.nml` (T21, all column fixes, cond_rh_crit=0.95, l_diag),
  `logs/rce_allfix_t42.nml` (T42). Run: `./libaeros/bin/rce_long.x <nml> <out.nc>`.
- **SpeedyWeather reference:** it couples physics into the tendency the standard way; its
  T21/T42 aquaplanet gives 66% free-trop RH, ~13 hPa/day subtropical ω. `docs/refs/
  speedy_comparison.md`. (Two agents this session extracted its `BulkRichardsonDiffusion`,
  `ImplicitCondensation`, and `BettsMillerConvection` in detail.)
- **Code:** `aeros_timestep.f90` (`aeros_timestep_step` — the whole step; step 6 is the
  forward-split to remove; `apply_phys_heating`, `raw_filter`), `aeros_semiimp.f90` (the
  implicit gravity-wave solve), `aeros_tendency.f90` (`wrk%dtdt`/`wrk%dt_phys`,
  `tnd%temp`), m2_results §12.1 (why forward-split was chosen).
- **Caveat on the SpeedyWeather comparison:** its default SST profile (302/273 cos²) and
  gray radiation differ from aeros's; resolution was controlled identically in the T21 test,
  which is what isolates the coupling as the cause.
