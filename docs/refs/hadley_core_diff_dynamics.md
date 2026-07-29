# Heating → divergence → ω: the core dynamical-core difference

Why aeros's Hadley/divergent **vertical** branch (ω) stays weak even when the
horizontal overturning |v| is strong, compared with SpeedyWeather.jl. A
read-only, code-level comparison of the two semi-implicit spectral cores.

Scope: this is *beyond* the known filter/diffusion damping. It isolates where
the two cores differ in the *forcing* of divergence by diabatic heating and in
the linear operator that converts that forcing to a balanced ω.

File:line evidence throughout. SpeedyWeather paths are under
`~/.julia/packages/SpeedyWeather/i8kOF/src/`; aeros paths under `src/`.

---

## 1. How SpeedyWeather converts heating → divergence → ω

**(a) Physics heating is accumulated into the temperature tendency BEFORE the
implicit solve.** In the primitive-equation `timestep!`
(`time_stepping/time_integration.jl:129-165`) the order is:

```
parameterization_tendencies!(vars, model)   # physics → vars.tendencies.*  (l.142)
dynamics_tendencies!(vars, lf2, model)       # adiabatic dynamics, += into same arrays (l.149)
implicit_correction!(vars, model.implicit, model)   # semi-implicit solve  (l.150)
horizontal_diffusion! ; leapfrog! ; transform!
```

Every dynamics routine accumulates with `+=` onto the tendency arrays that
already hold the physics contributions — e.g. `temperature_tendency!`
(`dynamics/tendencies.jl:569-615`, "`+=` because the tendencies already contain
parameterizations"). So at the moment `implicit_correction!` runs,
`vars.tendencies.temperature` = **physics heating + adiabatic dynamics**.

**(b) The implicit solve turns that temperature tendency into a divergence
tendency.** `implicit_primitive_single_kernel!`
(`time_stepping/implicit.jl:396-462`):

```
geopotential[lm,k] = Σ_r R[k,r] * temp_tend[lm,r]                 # l.416-423
G[lm,k]  = div_tend[lm,k] + ξ*eigenvalue*(U[k]*pres_tend[lm] + geopotential[lm,k])  # l.432-434
div_tend = S⁻¹ * G                                               # l.436-443  (δD = S⁻¹ G)
temp_tend += ξ L δD ;  pres_tend += ξ W δD                       # l.447-461
```

`R` is the hydrostatic geopotential-integration matrix (`implicit.jl:277-281`),
`eigenvalue = -l(l-1)` is ∇², and `ξ = α·dt`. So the diabatic heating enters
`δD = S⁻¹·(D_tend + ξ∇²(R·T_tend + U·pₛ_tend))` **directly and within the same
step**: heat a column → warmer T_tend → larger ∇²Φ → divergence tendency, with
`S⁻¹ = (I − ξ²∇²(RL + UW))⁻¹` (`implicit.jl:318-329`) inverting the gravity-wave
operator so the result is the *balanced* divergent response to the heating.

**(c) The linear operators use the ACTUAL mean temperature profile, including
its stratification.** `initialize!(ImplicitPrimitiveEquation,…)`
(`implicit.jl:210-347`): `temp_profile .= temp_average` (l.255) — the model's
horizontally-averaged profile, **not** isothermal. The temperature operator `L`
(`implicit.jl:283-310`) contains the vertical-advection-of-reference-temperature
terms `L1`, built from `ΔT_above`/`ΔT_below` (l.292-303) — the **reference
static stability**. `L = Diagonal(L0)*L1 + Diagonal(L2)*L3 + Diagonal(L4)`.

**(d) ω / σ̇ is diagnosed from the divergence** by Simmons–Burridge continuity,
`vertical_velocity!` (`dynamics/tendencies.jl:343-377`), and the adiabatic term
κTᵥ·Dlnp/Dt is assembled in `temperature_tendency!` (`tendencies.jl:617-650`)
from `div_sum_above`/`pres_flux`. Standard.

**(e) α = 1 (backward implicit) by default** for the primitive core
(`implicit.jl:130`, `α::NF = 1`), i.e. the strongest in-step gravity-wave /
divergence adjustment.

Net: a persistent heating is converted, every step, into the correct balanced
divergent circulation (a Gill/Matsuno-type response computed implicitly), with a
realistic-stratification operator.

---

## 2. How aeros does it

**(a) Diabatic heating is FORWARD-SPLIT — it never enters the semi-implicit
solve.** In `aeros_timestep_step` (`src/dynamics/aeros_timestep.f90:746-1066`)
surface, convection, condensation and radiation all write to `ts%wrk%dt_phys`,
a temperature **increment [K]** (l.820-933), *not* to the tendency that the
dynamics carries. The adiabatic RHS is assembled into `wrk%dtdt` in
`column_terms` (`src/dynamics/aeros_tendency.f90:581`), and
`aeros_tendency_spectral` analyses **only `wrk%dtdt`** into `tnd%temp`
(`aeros_tendency.f90:463`). `wrk%dtdt` is set *only* at line 581 — advection +
κTω/p — so `tnd%temp` is **pure adiabatic dynamics; it contains no diabatic
heating.**

The semi-implicit solve then runs on that heating-free tendency
(`aeros_timestep.f90:980-981` → `aeros_semiimp_step`). Only afterwards, in
**step 6, after the time-level swap and after the filter**, is the physics
heating applied — forward, onto T alone:

```
apply_phys_heating(s, now%spec, ts%wrk%dt_phys, …)   # l.1031-1033
  spec%temp(lm,k) += analysis(dt_phys(:,:,k))        # l.1085-1090
```

So within a step, a diabatic heating produces **zero** divergence and **zero**
pₛ response. It shows up in the divergence equation only one step later and only
indirectly, through the nonlinear geopotential `Φ = G·T` (`aeros_tendency.f90:562`,
`-∇²(Φ+K)` at `:459-461`) and the *linear* part of the next solve — after
having passed through hyperdiffusion (`diffuse`, divergence damped at
`aeros_timestep.f90:1631`) and the RAW filter (`raw_filter`, divergence damped at
`:1708-1710`, `eps=0.15`).

**(b) The semi-implicit operator itself DOES couple T→D — but for an ISOTHERMAL
reference, and it OMITS the reference stratification.** `aeros_semiimp_step`
(`aeros_semiimp.f90:436-542`) forms `tstar = told + ½h(rt − L_T(D))`
(l.497-499), `ldstar = lin_div(tstar,…) = λ(G·tstar + R_dT_ref·pstar)`
(l.502, `lin_div` :410-434), and solves `(I + ¼h²λW) Dⁿ⁺¹ = … + h(rd + ldstar −
ldnow)` (l.505-514). So *whatever is in* `rt = tnd%temp` **is** converted to
divergence — but per §2(a) `rt` holds no diabatic heating, so this pathway is
fed only the adiabatic tendency.

Crucially, the operator is linearized about an **isothermal** `t_ref`, enforced
by an `error stop` (`aeros_semiimp.f90:178-184`), and the thermodynamic matrix
`tmat` (τ) is built with only the `κT_ref·α` and `κT_ref·dlnp` terms
(`:220-226`) — it has **no reference-stratification term**, because that term is
∝ `T_ref,k+1 − T_ref,k`, which vanishes for isothermal `t_ref` (header
`:57-61`). This is exactly SpeedyWeather's `L1` (`ΔT_above/ΔT_below`) term —
**present there, structurally absent here.**

**(c) ω / geopotential are otherwise the SAME scheme.** `column_terms`
(`aeros_tendency.f90:475-600`): geopotential by upward hydrostatic integral
(`aeros_hydrostatic`, matching `gmat`), `-∇²(Φ+K)` for the divergence tendency,
and ω/p by Simmons–Burridge continuity (`:556-559`, `omga = cfac·G −
(scum·dlnp + α·cmass)/dp`). This is structurally identical to SpeedyWeather's
geopotential (`dynamics/geopotential.jl`) and `vertical_velocity!`. **The
ω-diagnosis is not the difference.**

**(d) α = 0.5 (centered).** The correction form `X_bar = (Xⁿ⁺¹+Xⁿ)/2`
(`aeros_semiimp.f90:69-98`) is centered-implicit, i.e. effective α = 0.5, half
SpeedyWeather's default backward α = 1.

---

## 3. Ranked structural differences that weaken aeros's ω

**#1 — Diabatic heating is forward-split, outside the semi-implicit solve.**
**[confirmed in code]** `aeros_timestep.f90:820-933` (physics → `dt_phys`),
`:1031-1033` + `:1085-1090` (`apply_phys_heating`, forward onto `now%temp` after
the swap), vs `tnd%temp` = adiabatic only (`aeros_tendency.f90:463,581`).
SpeedyWeather accumulates physics into `tendencies.temperature` *before*
`implicit_correction!` (`time_integration.jl:142-150`) where `ξ∇²R·temp_tend`
generates divergence (`implicit.jl:432-443`). *Mechanism:* in aeros a heating
produces no in-step divergence/pₛ; the divergent circulation is only forced
one-step-lagged through the nonlinear ∇²Φ, and that weak forcing is then
attacked by the same-step ∇⁶ diffusion and the RAW filter on divergence. This is
the largest and most certain difference and, on current `main`, is by itself
sufficient to explain a weak ω with a plausible |v| (|v| is set by the explicit
pressure-gradient/geopotential, which uses the real T; ω needs the *divergence*
response that never gets forced in-step).

**#2 — Isothermal reference + omitted reference-stratification in the implicit
temperature operator.** **[confirmed in code; effect plausible, needs test]**
`aeros_semiimp.f90:178-184` (isothermal enforced), `:220-226` (τ has no ΔT_ref
term) vs SpeedyWeather `implicit.jl:255` (real `temp_average`) and `:288-310`
(`L1` built from `ΔT_above/ΔT_below`). *Mechanism:* the implicit operator `S⁻¹ =
(I − ξ²∇²(RL+UW))⁻¹` sets how much *balanced* divergence a given in-step heating
generates and with what vertical-mode structure. A warm isothermal reference
(required to be ≥ column T for stability, header `:63-65`) is **over-stable for
the tropical troposphere** and its vertical normal modes (equivalent depths,
`aeros_semiimp_gwspeed`) project heating onto ascent differently than a realistic
stratified profile. Even once heating is coupled into the tendency (the next
task), this makes aeros under-convert heating to divergence relative to
SpeedyWeather — the residual reason ω can stay weak after |v| recovers. Needs a
test because the exact-steady-state (R=0) argument means its clearest effect is
on the *transient efficiency vs. damping* balance, not the asymptotic state.

**#3 — Centered (α=0.5) vs backward (α=1) implicit, plus RAW/diffusion on
divergence.** **[confirmed in code; known-suspect damping]**
`aeros_semiimp.f90:69-98` (α=0.5) vs `implicit.jl:130` (α=1);
divergence is RAW-filtered (`aeros_timestep.f90:1708-1710`, eps=0.15) and
∇⁶-diffused (`:1631`) every step. Weaker in-step adjustment + explicit damping on
the small divergent residual. Ranked below #1/#2 per the brief: this is the
"damping" the user already suspects; it compounds #1/#2 but is not the core.

**#4 — (refuted as the cause) geopotential build and ω-diagnosis.** **[confirmed
equivalent]** Both build Φ by the same upward hydrostatic T-integral
(`aeros_hydrostatic`/`gmat` vs `R` matrix), both force divergence with
`-∇²(Φ+KE)` (`aeros_tendency.f90:459-461` vs `bernoulli_potential!`
`tendencies.jl:936-965`), both diagnose ω/σ̇ by identical Simmons–Burridge
continuity (`column_terms:551-559` vs `vertical_velocity!:343-377`). A given
divergence yields the *same* ω in both. The difference is entirely upstream, in
how much divergence the heating generates — **#1 then #2.**

---

## 4. The single most likely core difference

**#1: aeros applies diabatic heating forward-split, so it never enters the
semi-implicit divergence solve.** On current `main` this is the whole story: the
divergence equation is forced by heating only through a one-step-lagged nonlinear
∇²Φ that is immediately damped, whereas SpeedyWeather converts heating to
balanced divergence *in-step* via `S⁻¹·ξ∇²R·temp_tend`. Coupling the heating into
`tnd%temp` before `aeros_semiimp_step` (the intended change) is exactly what
restores the in-step T→D pathway and is expected to spin up both |v| and ω.

**If, after that coupling, ω is still weak (the reported |v|→10×, ω-still-weak
result), the residual core difference is #2:** aeros's semi-implicit operator is
linearized about an isothermal, over-stable reference and structurally omits the
reference-stratification term (SpeedyWeather's `L1`). That term is precisely what
makes the implicit T↔D coupling carry the correct static stability, i.e. the
correct amount of ascent per unit heating. Restoring it requires a
non-isothermal reference in `aeros_semiimp` — currently forbidden by design
(`aeros_semiimp.f90:178-184`) — which is the deeper structural fix to test after
the coupling lands.
