# Subtropical descent concentration — the residual moist bias is a dynamics problem

**Status (2026-07-29): LARGELY RESOLVED — see [`docs/refs/sw_faithful_dynamics.md`](refs/sw_faithful_dynamics.md).**
The descent didn't concentrate because aeros had **no eddy momentum flux** — its RCE couldn't sustain
baroclinic eddies (they blew up at the model top). Reproducing SpeedyWeather's core numerics (unlock:
`eps_filter=0.1`, SW's Robert coefficient; + `si_alpha` backward-implicit + divergence diffusion) now
sustains eddies: descent tightened 25→69° down to **19–30°**, jet 14→**31 m/s @ 30°**, subtropical RH
94→**78%**. Remaining gap (+3.15 vs SW +6.4): the eddies are **mis-scaled** (m=8 jet-instability spike)
because aeros's subtropical jet is too strong/narrow — a coupled jet–eddy loop, now the active thread
(handoff ►►). The original diagnosis chain below is retained as the ruled-out record.

**The residual bias.** The free-tropospheric moist bias (subtropical/extratropical **upper
troposphere pinned near ~94% RH** vs SpeedyWeather/ERA5 ~45–62%) traces to aeros's **subtropical
subsidence being too weak and diffuse**: the Hadley descending branch spreads poleward
(near-uniform ~1.5–2 hPa/day from 19° to the pole) instead of **concentrating at ~30°** like
SpeedyWeather (+6→+11 hPa/day at 30°). Nothing dries the subtropical free troposphere below the
`rh_crit=0.95` condensation ceiling, so it sits at saturation.

---

## 0. The diagnosis chain — what has been ruled out (don't re-litigate)

1. **Column physics** (commit `a5395d6`): fixed the lower troposphere, cover (→0.70), TOA (→+3.8).
   Residual free-trop RH ~86%.
2. **Diabatic coupling** — `couple_diabatic` (commit `51bd92a`, **opt-in, default off = forward
   split, bit-for-bit**). Couples all diabatic heating into the semi-implicit tendency the
   standard way (rate `dt_phys/dt`, no leapfrog factor); needs `eps_filter~0.15` to hold the
   convective computational mode (`diff_taper`/`diff_adapt`/larger `conv_tau` did NOT; a stronger
   RAW filter did). Stable 100 d. **Strengthened the circulation but did NOT fix the RH bias.**
3. **Two-model comparison vs SpeedyWeather T21 aquaplanet** (`mwm/C_omega/`,
   `docs/refs/hadley_core_diff_{dynamics,forcing}.md`, `docs/refs/speedy_omega_structure.md`):
   - aeros's overturning **mass flux is NOT weak** — max|Ψ| 4.9×10¹⁰ kg/s > SpeedyWeather's
     3.8×10¹⁰.
   - The **dry core converts prescribed heating→ω correctly** (`q_force` hook, physics off) and is
     **filter-insensitive** (ω identical at eps 0.02 vs 0.15).
   - **Both models have a split (double) ITCZ** — not an aeros bug.
   - SpeedyWeather's quoted "~13 hPa/day" is its subtropical **DESCENT**, not peak ascent (its
     ascent is only −2 to −3.4).
   - **The gap:** aeros descent weak/diffuse (~2, spread 19°→pole, no 30° peak); SW **concentrated**
     (+6→+11 at 30°). aeros dries **bottom-up**, SW **top-down**; aeros's subtropical upper
     troposphere is uniformly pinned at the ~94% `rh_crit` ceiling.
4. **Moisture side RULED OUT (this session):**
   - `t_ref` sweep (250/275/300): ω insensitive to the isothermal semi-implicit reference.
   - **`vert_vanleer`** (commit `3de0b69`, opt-in van Leer vertical humidity transport vs the old
     first-order upwind): **no effect on RH** (subtropical UT-RH 94.2→94.1 over 100-day RCE).
     Vertical numerical diffusion is NOT the cause. Kept as a correct numerical improvement.
   - `rh_crit`: SpeedyWeather condenses toward **saturation (100%)** yet reaches ~50% RH via
     subsidence drying — so lowering aeros's `rh_crit=0.95` is a band-aid (would over-dry the ITCZ).

**Net:** core ✓, filter ✓, reference ✓, heating magnitude/scheme ✓, ITCZ shape ✓, humidity
transport ✓, condensation floor ✓ — all ruled out. **The root is the Hadley descending branch not
concentrating in the subtropics.**

---

## 1. The task

Understand and fix why aeros's descending branch spreads poleward instead of concentrating at
~30°. This is what sets the Hadley cell edge, so the candidates are dynamical:

- **Subtropical jet too weak** — aeros |u| ~20 m/s vs ERA5 ~30. A weak jet → ill-defined cell
  edge → diffuse descent. Check the jet structure and its relation to the descent latitude.
- **Baroclinic-eddy termination** — the Hadley edge is set by where baroclinic instability /
  eddy momentum fluxes take over. Are aeros's midlatitude eddies too weak/misplaced at T21?
- **Resolution (T21)** — T21 may be too coarse to resolve the jet + eddies that concentrate the
  descent. Test T42 (`logs/rce_allfix_t42.nml`). **Caveat:** SpeedyWeather's structure was the
  *same* at T21 and T42, so resolution may NOT be it for SW — but aeros must be checked directly
  (aeros T42 was "same" in the earlier weak-Hadley framing, but that was before this diagnosis;
  re-examine the *descent concentration* specifically, not just the RH average).
- **Numerical smearing** — the eps=0.15 filter (for `couple_diabatic`) or hyperdiffusion smearing
  the upper-level jet/eddies. But the bias is present in the **forward-split default (eps=0.06)**
  too, so the filter is not the sole cause.

Start by **not assuming** — build the meridional overturning streamfunction and the
eddy-momentum-flux / jet diagnostics for aeros vs SpeedyWeather and see what differs at the cell
edge.

---

## 2. Tools & configs (all landed)

- **omega/heating/RH dump** (`l_diag`): `rce_long` writes zonal-mean `omega`, `v`, `u`, `rh`,
  `cf`, `t`, `q`, and per-term heating (`q_cnv`/`q_cnd`/`q_rad`/`q_surf`/`q_net`) to the dump nc.
- **`q_force` hook** (`aeros_timestep_set_qforce`, driver `l_qforce`/`qforce_amp`): prescribed
  analytic heating, physics off, to probe the dry-core circulation response.
- **`couple_diabatic`** (opt-in, needs eps~0.15) and **`vert_vanleer`** (opt-in) — both default off.
- **`t_ref`** knob (isothermal semi-implicit reference).
- **SpeedyWeather reference:** `mwm/C_omega/omega_structure.jl` + `postproc.jl` (T21 L8 aquaplanet
  run + extraction); output regenerable, `docs/refs/speedy_omega_structure.md` has the numbers.
- **Configs:** `logs/rce_allfix.nml` (T21), `logs/rce_allfix_t42.nml` (T42). Run:
  `./libaeros/bin/rce_long.x <nml> <out.nc>`.
- **Overturning streamfunction** from ω: Ψ(φ,p) = −(2πa²/g) ∫ cosφ′ ω dφ′ (cumulative from a
  pole). Used ad hoc this session (gave aeros Ψ 4.9e10 vs SW 3.8e10); consider landing it as a
  diagnostic.

## 3. Acceptance

1. Subtropical descent **concentrates toward ~30°** and strengthens toward SpeedyWeather's
   +6→+11 hPa/day.
2. Subtropical/extratropical **upper-trop RH drops** from ~94% toward ~50–62%.
3. Cover/TOA hold; stable; 25 tests pass.
4. Bonus: the **subtropical jet** strengthens toward ERA5 ~30 m/s (the same fix likely closes both).

## 4. Expectation-setting

This is **open-ended research, not a targeted patch.** It may reveal a T21 resolution limit or a
genuine eddy/jet-dynamics deficiency rather than a one-line fix. Scope it realistically and report
early if it looks like a resolution wall.
