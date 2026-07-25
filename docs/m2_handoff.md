# M2 handoff — where to pick up

Written at the end of the session that made convection usable coupled (M2.3e–g).
`origin/main` is at the merge of `m2/sbm`. `make all openmp=1` is clean; all
**14** acceptance tests pass. Full detail for everything below is in
[`m2_results.md`](m2_results.md); this is the orientation and the next step.

## What landed this session (M2.3e–g)

The immediate next step from the previous handoff — *make convection usable
coupled* — is **done**. Both blockers (`m2_results.md` §11) are closed:

- **M2.3e forward-split heating** — moist-physics temperature changes are now a
  forward increment in `wrk%dt_phys`, applied to the `n+1` state after the
  dynamics step (decoupled from the centered leapfrog), the same forward
  treatment gridpoint humidity gets. This is what killed the computational-mode
  instability. Convection uses it; condensation keeps the centered `dtdt` path.
- **M2.3f Simplified Betts–Miller** — Frierson (2007), now the **default**
  scheme (`conv_scheme = "sbm"`); Manabe retained as `"manabe"`. Reference moist
  adiabat on the actual boundary-layer θe, MSE-buoyancy band (folds in the LCL),
  Frierson energy closure, implicit finite-τ relaxation. One Newton solve per
  level, no iteration. Conserves ∫(cp T + L q)dp to machine precision; precip =
  drying exactly. Knobs `conv_tau` (7200 s), `conv_rhref` (0.7).
- **M2.3g coupled** — `test_moist_run` runs convection + condensation together,
  200 steps, stable (no NaN, |u|≈9 m/s), convection actively raining, q≥0, water
  closing to 1.2×10⁻³.

Convection is still **off by default** (`convect = .FALSE.`): it is validated by
*conservation and stability*, not yet by *climate*. That waits on radiation.

## The immediate next step: radiation + ERA5

This is the big remaining M2 piece and the gate on validating everything in the
moist line (condensation, transport, convection) against observations rather than
against conservation laws. AR said ERA5 data is coming.

- **Design** is `design.md` §5: a bespoke ecCKD gas-optics table is preferred; the
  SESAM scheme is the pragmatic fallback. Call radiation every 1–3 h with
  per-step surface-flux rescaling (§6 budget is ample).
- **Why it unblocks convection:** with radiative cooling driving the instability
  and ERA5 to compare against, `convect` can be turned on and SBM's `τ` and
  `rh_ref` tuned (they are namelist knobs precisely so they can be). Until then
  turning convection on is untested against climate — the machinery works, the
  parameters are Frierson's defaults, unvalidated for this model.
- **Suggested order:** radiation operator (grid seam, like the other physics) →
  ERA5 ingestion for boundary/validation → turn on the moist physics stack and
  tune against ERA5.

## Other open M2 items (not blocking)

- **Shallow convection is temperature-only.** SBM's shallow branch (dry columns,
  `Pq ≤ 0`) currently does a non-precipitating heat redistribution and leaves
  humidity untouched — a deliberate simplification, because the additive `q_ref`
  shift of full Betts–Miller can drive a very dry column's humidity negative.
  Restore the moisture-redistributing shallow branch when there is something to
  tune its reference against (i.e. with radiation/ERA5). `m2_results.md` §12.2.
- **Manabe multi-layer adjustment** — now lower priority since SBM is the
  default and has no convergence loop, but if Manabe is ever used seriously its
  pairwise sweep is still O(tens of passes) on a cold start. The multi-layer fix
  is described in the `aeros_convection` header. `m2_results.md` §11 problem 1.
- **Water closure characterization.** The FV-vs-spectral air-mass gap is ~1–3×10⁻³
  at T21 (§8–9, §12.3). Check at T42 to confirm it shrinks with resolution, and
  over a long run to confirm it does not accumulate secularly (the way the mass
  leak did).
- **Transport accuracy leftovers** (§10): vertical van Leer (horizontal is done),
  and per-row polar sub-cycling instead of global (a cost optimization).
- **Carried from M1:** `tau_diff`/∇⁶-vs-∇⁴ retune (jets run ~10% strong);
  `SHqst_to_spat` 11→10 transforms; the Albedo 16–128 thread sweep to close the
  cliff question beyond 10 threads.

## Gotchas the next session will want

- **Physics seam, now two coupling paths.** HS forcing and condensation add to
  `wrk%dtdt` (grid) and ride the centered leapfrog. **Convection** adds a forward
  increment to `wrk%dt_phys` (grid), applied to the `n+1` spectral state in
  `aeros_timestep_step` step 6 — forward, off the centered step. New radiation
  heating is smooth and can use either; `dtdt` is simplest unless it turns out
  large and stiff. Humidity (`qv_g`) is gridpoint and advected in step 8.
- **Worktree build needs the `fesm-utils` symlink.** A fresh worktree lacks it
  (`fesm-utils` is a symlink to `/Users/alrobi001/models/fesm-utils` in main).
  Before `make` in a worktree: `ln -sf /Users/alrobi001/models/fesm-utils
  <worktree>/fesm-utils`, and copy a generated `Makefile` in (it is gitignored;
  `config/Makefile` is the configme template). The `Bash` tool's cwd *persists*
  between calls here, but use `git -C <worktree>` for git regardless.
- **Generated `Makefile` is gitignored.** After merging a branch that adds a
  *new* test target, refresh main's `Makefile` from the template or the new test
  won't build on a fresh checkout. This session only modified existing targets,
  so no refresh was needed — but it has bitten before.
- **`timeout` is not on macOS** (it's `gtimeout`); don't wrap test runs in it.
- **Moist IC trap**: seeding humidity as a fraction of `q_sat` against an
  isothermal-warm column injects absurd values aloft (`q_sat` explodes at low
  `p`); use a realistic lapse-rate profile. Both moist-run initial states here
  follow that rule.
- **SBM test columns are SBM-specific.** A "saturated but cold" column that
  looks conditionally unstable to Manabe reads as *shallow* to SBM (its absolute
  humidity is below 70% of the warm reference adiabat's saturation, so nothing
  rains). Deep convection needs a *warm, humid* unstable column; shallow needs a
  *moist boundary layer under a dry free troposphere*. See `test_convection`.
