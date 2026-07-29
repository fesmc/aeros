# Hadley core difference — heating structure vs. heating→circulation coupling

Read-only comparison of **aeros** (T21 L12 aquaplanet RCE, `logs/rce_allfix.nml`) and
**SpeedyWeather.jl** (`~/.julia/packages/SpeedyWeather/i8kOF/src`) aimed at one question:
why does SpeedyWeather have a vigorous Hadley cell (subtropical ω ~13 hPa/day, free-trop
RH ~66%) at the same T21 aquaplanet where aeros has a weak one (ω ~1–1.4 hPa/day, RH ~95%
at σ 0.2–0.6)?

**Headline finding, stated up front.** The diabatic heating *structure* is **not** the
material difference — aeros's net heating profile is structurally correct and close to
SpeedyWeather's (latent heating on the moist adiabat in ascent, ~1.5–2 K/day radiative
cooling in the free troposphere, stratospheric SW balance). The material difference is
**where the heating is inserted in the time integration**: SpeedyWeather adds every
diabatic tendency to the thermodynamic tendency *before* the semi-implicit correction and
leapfrog, so the divergent gravity-wave adjustment to the heating happens **within the
step**; aeros forward-splits **all** diabatic heating onto the n+1 state *after* the
semi-implicit solve, after diffusion/sponge/vdiff, and after the RAW filter. The
divergence tendency never sees the heating directly. This is the "core difference, more
than damping" — it is a structural coupling difference in the integrator, not a
physics-scheme or damping difference. It is corroborated by an internal aeros control:
Held-Suarez forcing *is* on the centered in-solve path and produces a normal circulation,
while the RCE diabatic heating on the forward-split path does not.

> **Correction to the task framing.** The task states the diabatic heating was "coupled
> into the thermodynamic tendency before the semi-implicit solve." The code does **not**
> yet do this. As of `main` (commit `e40c323`), coupling into the tendency is the *planned
> next task* (`docs/weak_hadley_scope.md`, memory `weak-hadley-forward-split-coupling`).
> The ~10× |v| strengthening + stabilization the task describes came from the **column-
> physics moist-bias fix** (commit `a5395d6`: dry-static-energy vdiff, condensation
> rh_crit=0.95 + reevaporation, `sbm_frierson` default), **not** from tendency coupling.
> The residual weak-ω / saturated-mid-troposphere problem is *because* the heating is still
> forward-split. This document analyzes the code as it actually is.

---

## 1. Heating-structure comparison (vertical profile / magnitude, with code evidence)

Both models are Simplified-Betts–Miller + large-scale-condensation + grey/band radiation +
bulk surface fluxes + Richardson boundary-layer diffusion. Term by term the structures are
close; none of the differences below is large enough to explain a 10× circulation gap.

### 1a. Radiation (the descending-branch cooling that subsidence must balance)

| | SpeedyWeather (default `PrimitiveWetModel`) | aeros (RCE default) |
|---|---|---|
| Scheme | `OneBandLongwave` — Frierson (2006) grey two-stream. `parameterizations/radiation/longwave_radiation.jl:160-284`, transmissivity `longwave_transmissivity.jl:25-70` | ecCKD Malkmus/Goody correlated-k, 10 LW g-points / 4 bands. `src/physics/aeros_ecckd.f90:7-60`; selected `aeros_radiation.f90:203`, `input/rce_defaults.nml:8` |
| Cooling profile | optical depth `τ(σ)=τ₀(0.1σ+0.9σ⁴)`, `τ₀_eq=6`, `τ₀_pole=1.5` (`longwave_transmissivity.jl:60-66`, `:27-34`). σ⁴ dominance ⇒ bottom-weighted; net flux-divergence cooling ~1–2 K/day in free trop, decaying to the tropopause | net-upward LW flux built on interfaces, layer heating = flux divergence `heat(k)=(g/cp)·(fnet(k)−fnet(k-1))/dp` (`aeros_radiation.f90:509-511`). Emergent ~1.5–2 K/day free-trop cooling, decaying to tropopause |
| Model-top balance | SW ozone absorption `σ→50·max(0,0.2−σ)` (`shortwave_radiation.jl:133-149`) | SW ozone deposited in stratosphere, peak ~20 hPa (`aeros_radiation.f90:855-891`, `:188-195`) |
| Prescribed option | `UniformCooling` (Paulius–Garner): flat −1.5 K/day above 207.5 K, stratospheric relaxation (`longwave_radiation.jl:4-50`) | none packaged; cooling is always the correlated-k flux divergence |

**Verdict:** aeros's radiative cooling is neither too weak nor Newtonian — it is a genuine
correlated-k flux-divergence profile, full column, magnitude set by the actual T/q. Its
free-troposphere cooling is *comparable in shape and magnitude* to SpeedyWeather's default
grey scheme. Subtropical subsidence is **not** starved by a missing cooling term; it is
starved because the divergent response to the ascending-branch heating is missing (§3).

### 1b. Convection (ascending-branch latent heating)

Both are Frierson (2007) Simplified Betts–Miller with target RH **0.7**.

- SpeedyWeather `BettsMillerConvection`: `time_scale = 4 h`, `relative_humidity = 0.7`;
  reference profile from surface parcel to **level of zero buoyancy**, relaxes T,q over τ;
  deep (enthalpy-conserving shift) vs shallow branches. `parameterizations/convection.jl:8-14,
  65-66, 89-138`. Heating deposited surface→LZB.
- aeros default `sbm_frierson` (`conv_scheme="sbm_frierson"`, `rce_defaults.nml:16`):
  `conv_tau = 7200 s` (~2 h), `rh_ref = 0.7`; connected layer surface→LZB of a lifted
  **surface** parcel, integral Pt/Pq trigger, implicit relaxation
  `t(k)+=a·(tref+shift−t)`, `q(k)+=a·(qref−q)`. `src/physics/aeros_convection.f90:137-143,
  499-616`. Heating deposited on the moist adiabat surface→LZB; **stays active at
  equilibrium** (the `sbm` variant that went dormant, `:54-64`, is not the default).

**Verdict:** near-identical structure. The only knob difference is τ (aeros 2 h vs
SpeedyWeather 4 h) — aeros relaxes *faster*, if anything a stronger heating, not weaker.
Latent heating is at the right levels in both.

### 1c. Large-scale condensation

- SpeedyWeather `ImplicitCondensation`: RH threshold **0.95**, τ = 3·Δt, implicit
  saturation adjustment `δT=−Lᵥ/cp·δq`, reevaporation into subsaturated layers below.
  `parameterizations/large_scale_condensation.jl:6-24, 93-158`.
- aeros: `rh_crit = 0.95` (`cond_rh_crit=0.95`), top-down saturation adjustment (3-Newton),
  `dt_phys += (Lᵥ/cp)·dqc`, reevaporation `reevap=30`. `src/physics/aeros_condensation.f90:
  191-286, 217, 272`.

**Verdict:** same threshold, same MSE-conserving saturation-adjustment structure, both
reevaporate. No material difference.

### 1d. Surface fluxes + boundary-layer diffusion

- SpeedyWeather: bulk `flux=ρ·C·V·(Ts−Ta)`, ocean drag 0.9e-3 into the lowest layer
  (`surface_fluxes/heat.jl:94,108`); `BulkRichardsonDiffusion`, Ri_c=10, diffuses dry static
  energy + humidity within the diagnosed BL (`vertical_diffusion.jl:4-31,164-240`).
- aeros: `SH=ρ·cp·c_h·V·(SST−T₁)`, `c_h=c_e=1.5e-3`, `u_min=1`, into the lowest layer
  (`src/physics/aeros_surface.f90:203,217`); `aeros_vdiff` Richardson K-profile (Frierson
  2006), Ri_c=10, diffuses **dry static energy** and humidity, **implicit** backward-Euler
  (`aeros_vdiff.f90:262-274,304-380`, applied `aeros_timestep.f90:1010`).

**Verdict:** same design (SpeedyWeather's scheme is literally what aeros copied in commit
`a5395d6`). aeros `c_h` is slightly larger (1.5e-3 vs 0.9e-3) — again a *stronger* surface
source, not weaker.

### 1e. Vertical grid (rules out level placement as the cause)

- SpeedyWeather default: **L8**, equidistant σ, uniform Δσ=0.125
  (`dynamics/spectral_grid.jl:11`, `vertical_coordinates.jl:44-54`, `geometry.jl:85-91`).
- aeros: **L12**, stretched hybrid σ (stretch_a=0.4, stretch_r=2.0, σ_t=0.2), top 10 hPa,
  ~5 full levels below 850 hPa (`src/dynamics/aeros_vertical.f90:126-284`).

**Verdict:** aeros has *more* levels and *finer* low-level resolution than SpeedyWeather's
default, with a higher top. If level placement mattered, aeros would be favored. The scope
doc confirms resolution was directly ruled out (aeros T42 identical; SpeedyWeather T21 = T42
= dry). Grid is **not** the cause.

### 1f. Net heating profile — the structures agree

| Region | aeros | SpeedyWeather | Match? |
|---|---|---|---|
| Ascent mid/upper trop | convective latent heating on moist adiabat, surface→LZB | same (BM, surface→LZB) | ✔ |
| Free trop (subsidence) | correlated-k LW cooling ~1.5–2 K/day, top-decaying | grey LW cooling ~1–2 K/day, top-decaying | ✔ |
| Strat / model top | SW ozone absorption balances LW | SW ozone absorption balances LW | ✔ |
| Boundary layer | surface SH into lowest layer, mixed by Ri-vdiff | same | ✔ |

The heating that *should* drive the cell is present and correctly shaped in aeros. What
differs is whether the dynamics can *respond* to it within the timestep — §3.

---

## 2. Humidity / dynamics coupling assessment

This is where aeros and SpeedyWeather genuinely diverge, and it compounds the integrator
problem below.

**SpeedyWeather — humidity fully coupled to the dynamics (two paths):**
1. Specific humidity is a **spectral prognostic** variable, advected in spectral space
   (`models/primitive_wet.jl:158-162`; `dynamics/tendencies.jl:652-668`).
2. It couples back into the momentum/mass fields through **virtual temperature**
   `Tᵥ=virtual_temperature(T,q)` used in the geopotential and pressure-gradient terms
   (`dynamics/tendencies.jl:509,640`).
3. Its **latent heating** enters `temp_tend`, which is added to the tendency *before* the
   implicit solve (§3).

So in SpeedyWeather, moisture drives circulation through *both* virtual-temperature buoyancy
*and* latent heating, and both are inside the semi-implicit dynamics.

**aeros — humidity almost entirely decoupled from the divergent dynamics:**
1. Humidity is **off-spectral**: a gridpoint field `now%qv_g` transported by a separate
   finite-volume scheme using the stepped winds (`aeros_timestep.f90:1051`, step 8;
   `aeros_transport`). Transport itself is fine, but q never becomes a spectral field.
2. **No virtual-temperature coupling.** The hydrostatic/geopotential operator takes `temp`,
   and its own header notes Tᵥ is a deferred M2 edit (`aeros_vertical.f90:436-439`); a search
   of `aeros_tendency.f90` and `aeros_timestep.f90` finds **no** virtual-temperature term
   (no `0.608`, no `virtual`). So moisture does **not** modify the geopotential/pressure
   gradient. Path (2) above is **severed**.
3. Latent heating (convection + condensation) is deposited into `wrk%dt_phys` and applied
   **forward-split** after the solve (`aeros_condensation.f90:191-286`,
   `aeros_convection.f90:321`; applied `aeros_timestep.f90:1021-1033`). Path (3) is
   **outside** the dynamics.

**Consequence.** In aeros the *only* moisture→dynamics link is the forward-split latent
heating — and that link is off the centered path. The moist-circulation feedback (latent
heating → column warming → geopotential rise → divergence → ascent → moisture convergence →
more latent heating) is broken at **two** places: no virtual-temperature buoyancy, and
latent heating that never enters the divergence tendency. The mid-troposphere stays
saturated because condensation keeps firing there but the ascent that would ventilate it is
never dynamically generated. Note the recent change *did* make condensation's heating and
drying **share one discretization** (both forward-split gridpoint,
`aeros_timestep.f90:1027-1028`) — good for column-MSE conservation, but it consolidated them
on the *wrong side* of the semi-implicit solve for driving circulation.

**Relative weight.** The severed latent-heating→divergence link (§3) is the dominant
starvation; the missing virtual-temperature coupling is a real but secondary amplifier.
Both point the same direction, so tests should target the integrator coupling first.

---

## 3. The integrator coupling difference (the core mechanism)

**SpeedyWeather — diabatic heating is inside the semi-implicit solve.**
`time_stepping/time_integration.jl:129-157` (the `PrimitiveWetModel` timestep):
```
142  parameterization_tendencies!(vars, model)          # radiation, convection,
                                                         # condensation, surface, vdiff
149  dynamics_tendencies!(vars, lf2, model)             # dynamical core adds to same tend
150  implicit_correction!(vars, model.implicit, model)  # semi-implicit gravity-wave solve
157  leapfrog!(vars, dt, lf1, model)                    # advance
```
The physics tendencies are in `temp_tend`/`humid_tend` **before** `implicit_correction!` and
`leapfrog!`. The semi-implicit gravity-wave adjustment therefore acts on the *combined*
(diabatic + dynamical) tendency: heating → divergence → ω is generated within the step.

**aeros — diabatic heating is bolted on after everything.**
`aeros_timestep_step` (`aeros_timestep.f90`):
- `wrk%dtdt` / `tnd%temp`, the tendency that feeds the semi-implicit solve, carries **only**
  advective + adiabatic heating: `dtdt = −(u·∇T) + κTω/p` (`aeros_tendency.f90:581,588`,
  transformed at `:463`). No diabatic term is in it.
- Step 1 accumulates surface + convection + condensation + radiation into `wrk%dt_phys`
  (`aeros_timestep.f90:836-933`).
- Step 2 does the semi-implicit advance (`:980-981`); steps 3/3b/3c apply
  diffusion/sponge/vdiff; step 4 applies the RAW filter (`:1013`).
- **Step 6** finally adds `dt_phys` to the n+1 temperature, forward-split
  (`apply_phys_heating`, `:1021-1033`).

So the divergence tendency is computed and the gravity-wave solve is done **before** the
heating is ever applied. The overturning can only respond one step later, indirectly,
through advective feedback — and that lagged response is further damped by the RAW filter
(`eps_filter=0.06`, kept high specifically to suppress the computational mode the sharp
forward-split heating otherwise excites; `aeros_timestep.f90:815, 1021-1030`). The
ascending-branch heating is present but produces almost no divergence, hence ω ~1 hPa/day.

**Internal control already in the code.** aeros's own Held-Suarez forcing writes its
Newtonian relaxation to `wrk%dtdt` — the **centered, in-solve path**
(`aeros_held_suarez.f90:210`) — not to `dt_phys`. HS therefore produces a normal
circulation, while the RCE diabatic heating on the forward-split path does not. That
HS-vs-RCE contrast, entirely within aeros, is itself strong evidence that the integrator
coupling — not the heating magnitude or the dynamical core — is the cause.

---

## 4. Prioritized isolated-test protocol (bisection)

Each test isolates one layer. Tests are ordered so the cheapest, most decisive come first;
the first two are runnable **inside aeros alone** and can settle the hypothesis before any
cross-model work.

### Test 1 (FIRST — decisive, aeros-only) — prescribed analytic Q, centered vs forward-split
**What it isolates:** the heating→circulation transfer of aeros's own dynamical core,
holding the heating *identical* and flipping only the insertion point. Removes every
physics-scheme and cross-model confound.
**Config (aeros):** all physics off (`l_surf=l_cnv=l_cnd=l_rad=l_vdiff=.false.`, sponge as
usual). Inject a fixed analytic heating `Q(σ,φ) = Q₀·g(σ)·cos φ` (e.g. `g(σ)` a half-sine
peaking at σ≈0.4, `Q₀`≈2 K/day, balanced by uniform cooling so the global mean is ~0). Run
it **two ways**:
- (1a) add Q to `wrk%dtdt` before the solve (centered/in-solve — the standard path, the same
  slot HS uses);
- (1b) add Q to `wrk%dt_phys` (forward-split — the current RCE path).
Compare zonal-mean ω and upper-branch |v| (the `omega` NetCDF diagnostic, commit `c8120b8`).
**Incriminating outcome:** if (1a) gives ω an order of magnitude stronger than (1b) for the
*same* Q, the forward-split coupling is confirmed as the cause, independent of SpeedyWeather
and of all physics schemes. This is the single most decisive experiment and it is essentially
the minimal version of the planned fix.

### Test 2 (aeros-only control) — dry Held-Suarez vs RCE, same core
**What it isolates:** whether aeros's dry dynamical core can build a strong meridional cell
*at all* when the forcing is on the centered path.
**Config (aeros):** run the existing Held-Suarez setup (`aeros_held_suarez.f90`, forcing →
`dtdt`, `:210`) and read the same ω/|v| diagnostics; compare against the RCE run.
**Incriminating outcome:** HS (in-solve relaxation) yields a normal ω/|v| while RCE
(forward-split) is ~10× weaker ⇒ the core is fine and the coupling is the difference. (If HS
were *also* weak, the core/diffusion would be implicated instead — but the existing HS
validation says it is not.)

### Test 3 (cross-model core check) — identical prescribed Q, all physics off, both models
**What it isolates:** whether the two dynamical cores transfer heating→circulation the same
way, once aeros's coupling is put on the standard path.
**Config:** same `Q(σ,φ)` as Test 1. aeros: variant (1a) (Q in `dtdt`). SpeedyWeather:
`PrimitiveWetModel`/`PrimitiveDryModel` with all parameterizations disabled and a custom
`temp_tend += Q` forcing (its natural slot, before `implicit_correction!`). Match truncation
(T21), levels, and Δt as closely as possible.
**Incriminating outcome:** if aeros-(1a) now matches SpeedyWeather's ω for the same Q, the
cores agree and the *entire* gap is the coupling. A residual gap would point to a secondary
core/diffusion difference (spectral truncation of the divergent modes, hyperdiffusion on
divergence, semi-implicit reference state) to chase next.

### Test 4 (radiation-driven overturning) — prescribed SST + radiation only, no moist physics
**What it isolates:** the radiatively-driven component of the cell, free of convection/
condensation and of moisture coupling.
**Config:** aeros `l_rad=.true., l_surf=.true., l_cnv=l_cnd=l_vdiff=.false.` on the fixed
aquaplanet SST; SpeedyWeather with radiation + surface only, convection/condensation off,
`AquaPlanet` SST (`ocean.jl:226-265`, `temp_equator=302, temp_poles=273`).
**Incriminating outcome:** aeros still weak here ⇒ the starvation is not moist-physics-
specific but generic to any forward-split diabatic term (expected, since radiation is also
in `dt_phys`). If aeros were *strong* here, it would shift suspicion onto the moist terms
specifically.

### Test 5 (rule-out cleanups — low priority, several already done)
- **Time filter:** rerun aeros with smaller `eps_filter` once heating is on the centered
  path (the high value exists to tame the forward-split mode; coupling should let it drop).
- **Levels/σ placement:** run SpeedyWeather at L12 with aeros's stretched σ, and/or aeros at
  SpeedyWeather's L8 uniform σ, to confirm grid is neutral (scope doc already indicates it
  is).
- **Virtual temperature:** with Test 1 settled, add Tᵥ to aeros's geopotential and re-measure
  ω to size the secondary moisture-buoyancy path (§2).
These only matter after Tests 1–3; the scope doc already rules out sponge, hyperdiffusion,
RAW strength, SST gradient, and resolution.

---

## 5. Top hypothesis and the single first test

**Top hypothesis.** The core difference is **the time-integration coupling of diabatic
heating, not the heating structure**. aeros's per-term heating profiles (latent heating on
the moist adiabat, ~1.5–2 K/day correlated-k free-troposphere cooling, stratospheric SW
balance, Richardson BL mixing) are structurally correct and close to SpeedyWeather's.
SpeedyWeather inserts all of that heating into the thermodynamic tendency **before** the
semi-implicit correction and leapfrog (`time_integration.jl:142→150→157`), so the divergent
gravity-wave adjustment builds the Hadley overturning within each step. aeros carries **only**
adiabatic/advective heating in `dtdt`/`tnd%temp` and forward-splits **all** diabatic heating
onto the n+1 state after the solve, diffusion, and RAW filter
(`apply_phys_heating`, `aeros_timestep.f90:1021-1033`). The divergence tendency never sees
the heating, so ω stays ~1 hPa/day and condensation keeps the mid-troposphere saturated for
lack of ventilating ascent. A missing virtual-temperature coupling (§2) removes a second
moisture→dynamics path and compounds the effect. This is exactly the "more than damping"
core difference the user senses, and it matches the confirmed diagnosis in
`docs/weak_hadley_scope.md`.

**Single first test to run: Test 1 — aeros-only, prescribed identical analytic `Q(σ,lat)`
with all physics off, run once with `Q` added to `wrk%dtdt` (centered, in-solve) and once
with `Q` added to `wrk%dt_phys` (forward-split), comparing zonal-mean ω and upper-branch
|v|.** It is cheap (no SpeedyWeather run needed), it changes only the insertion point of an
otherwise identical forcing, and a ~10× ω difference between the two variants directly
incriminates the coupling and predicts that moving the diabatic heating into `tnd%temp`
(the planned fix) will restore the cell. Test 2 (HS-vs-RCE, already-built) is the immediate
confirmatory control.

---

### Key file:line evidence
- SpeedyWeather physics-before-implicit: `time_stepping/time_integration.jl:142,149,150,157`
- SpeedyWeather humidity spectral + Tᵥ coupling: `models/primitive_wet.jl:158-162`;
  `dynamics/tendencies.jl:509,640,652-668`
- SpeedyWeather LW cooling: `parameterizations/radiation/longwave_transmissivity.jl:60-66`;
  `longwave_radiation.jl:160-284`; convection `parameterizations/convection.jl:8-14,89-138`
- aeros `dtdt` = adiabatic only: `src/dynamics/aeros_tendency.f90:581,588,463`
- aeros forward-split of all diabatic heating: `src/dynamics/aeros_timestep.f90:836-933,
  1021-1033` (`apply_phys_heating`)
- aeros HS on centered path (internal control): `src/physics/aeros_held_suarez.f90:210`
- aeros no virtual temperature: absence in `aeros_tendency.f90` / `aeros_timestep.f90`;
  deferral noted `src/dynamics/aeros_vertical.f90:436-439`
- aeros LW flux-divergence cooling: `src/physics/aeros_radiation.f90:509-511`; convection
  `sbm_frierson` `src/physics/aeros_convection.f90:499-616`; grid `aeros_vertical.f90:126-284`
