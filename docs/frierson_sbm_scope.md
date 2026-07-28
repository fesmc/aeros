# Free-tropospheric drying via a working Simplified Betts–Miller — scope / handoff

**Status:** design memo, for a fresh session. Not started.
**Goal:** cure the near-saturated free-troposphere moist bias (m2_results §23) — RH ~85–100%
everywhere convection is not actively firing, vs ERA5 ~45% — which forces the overcast /
cloud-cover runaway and blocks a physical coupled TOA balance.

---

## 0. Important correction to the framing

Earlier in the orchestration this was called "add subsidence drying," and a comparison
study noted SpeedyWeather uses Simplified Betts–Miller (Frierson 2007) while "aeros is
Betts–Miller-*style*." Reading `src/physics/aeros_convection.f90` corrects that:

**aeros already implements the Frierson (2007) Simplified Betts–Miller scheme.** It is
`sbm_adjust` (`SCHEME_SBM`), the **default** convection scheme, and it is genuinely the
Frierson closure: a whole-column reference moist adiabat, `q_ref = rh_ref·q_sat`, a deep
(precipitating) and a shallow (non-precipitating) branch, enthalpy-conserving, relaxed
implicitly over `tau`. The *pairwise adjacent-layer* scheme (the one that "converges only
geometrically") is **Manabe** (`SCHEME_MANABE`), a separate opt-in alternative — not what
the RCE runs use.

**So the task is NOT to implement Frierson.** It is: **diagnose why aeros's Frierson-SBM
goes dormant at RCE equilibrium and does not hold the free troposphere near `rh_ref`, when
SpeedyWeather's implementation of the same scheme does — and fix that divergence.**

The neither-references-use-a-bolt-on-subsidence-term point still stands: SPEEDY dries via
mass-flux compensating subsidence; SpeedyWeather dries via *this* scheme kept active. We
have the scheme; we need it to stay on.

## 1. The precise failure (from §23 + the code)

- `sbm_adjust` anchors the reference adiabat on the **actual near-surface parcel MSE**,
  `hb = cp·T(nlev) + Φ(nlev) + L·q(nlev)` (actual `q`, not `q_sat`), and only adjusts the
  band where the parcel is buoyant in the MSE sense, `hb > h*_env(k)`.
- §23 (`rce_long` `l_diag` buoyancy readout): once the column equilibrates, `hb − h*_env < 0`
  at **every** level (surface ≈ −1.5 kJ/kg), so the band is empty (`kb == 0`) and the scheme
  returns the column unchanged. Convection is **dormant**; nothing then dries the free
  troposphere, and large-scale condensation only caps supersaturation at 100%.
- The code's own header explains the anchor choice: the *saturated* anchor was rejected
  because it made the reference "so warm that `rh_ref q_sat(T_ref)` exceeds the environmental
  humidity everywhere and no column ever rains." The actual-`q` anchor fixed *that*, but is a
  prime suspect for the opposite failure (a boundary layer that, once dried, can never
  re-trigger).

## 2. Hypotheses to test (diagnosis-first — do not just start editing)

Use the existing `rce_long` `l_diag` buoyancy diagnostic (`hb` vs `h*_env` per level, per
lat) at equilibrium as the primary instrument. Candidate root causes, roughly in likelihood
order:

1. **Anchor / trigger.** The actual-`q` surface anchor + strict `hb > h*_env` band means once
   convection dries the boundary layer, `hb` drops and convection cannot re-fire, even though
   surface fluxes are still adding heat+moisture. Frierson lifts a **boundary-layer-mean /
   mixed parcel** and re-establishes CAPE each step; compare the exact trigger and parcel
   definition against Frierson 2007 §2 and the SpeedyWeather `SimplifiedBettsMiller` source.
2. **Surface-flux ↔ boundary-layer coupling.** In a real RCE, surface evaporation keeps the
   sub-cloud layer marginally buoyant so convection fires continuously. Check whether aeros's
   surface fluxes + vdiff actually maintain sub-cloud MSE, or whether the boundary layer
   starves (interaction with `aeros_surface`, `aeros_vdiff`, and the seam ordering).
3. **Deep vs shallow branch.** If at equilibrium `Pq = ∫(q − q_ref)dp ≤ 0`, the **shallow**
   branch fires and relaxes **temperature only** (q untouched by design). Verify the model
   isn't stuck shallow — which would explain "T adjusts but q never dries."
4. **`rh_ref` value.** aeros uses 0.7; Frierson/SpeedyWeather use **0.8**. Almost certainly
   not the root cause (0.7 is drier), but align it once the mechanism is fixed.

## 3. Fix + acceptance (once the divergence is identified)

- Fix the identified cause (most likely the parcel/trigger so convection stays active under
  continuous surface forcing), keeping the change **opt-in / bit-for-bit off** if it alters
  `sbm_adjust` behaviour, or as a new `conv_scheme` variant if it is a distinct closure.
- **Acceptance:**
  1. All acceptance tests bit-for-bit unless the change is opt-in and enabled (then
     `test_convection` gets a case for the new behaviour).
  2. In the rotating RCE (`logs/rce_revalidate.nml`), free-tropospheric **RH drops from
     ~85–100% toward ERA5 ~45%** (use `scripts/rce_humidity_vs_era5.jl` / the `rce_long` RH
     dump), and **convection is non-dormant at equilibrium** (`hb − h*_env` shows an active
     band).
  3. **Cloud cover** falls from ~1.0 toward ERA5 ~0.63, and **TOA net** moves toward balance
     without the `cond_rh_crit` sub-saturation crutch (which was the §23 stopgap for exactly
     this missing drying — it can then be relaxed back toward 1.0).
  4. No new instability (run it with the land case; adaptive hyperdiffusion is now available
     as a safety net if needed).

## 4. Reference material

- **Frierson (2007)**, *The Dynamics of Idealized Convection Schemes…*, J. Atmos. Sci. — §2
  defines the reference profiles, the deep/shallow split, and the trigger.
- **SpeedyWeather.jl** `SimplifiedBettsMiller` — a clean, readable implementation of the same
  scheme to diff against (`rhbm = 0.8`, CAPE trigger, deep+shallow); see
  `docs/refs/speedy_comparison.md`.
- **Isca** "Simple Betts–Miller" module docs — the same Frierson scheme, well documented.
- aeros: `src/physics/aeros_convection.f90` (`sbm_adjust`, `moist_adiabat_temp`), the `l_diag`
  buoyancy diagnostic in `drivers/rce_long.f90`, and m2_results §22–23.

## 5. Why this is high priority

The moist bias "keeps confounding everything" (user): it drives the overcast cloud runaway,
the OLR bias, and the coupled TOA imbalance, and it is the reason the diagnostic cloud scheme
runs away and the prognostic one needed a sink to be robust. Fixing the drying at its source
(convection staying active and drying to `rh_ref`) is the highest-leverage single fix
remaining — it should improve clouds, OLR, and TOA balance simultaneously, and let the
`cond_rh_crit` crutch retire.
