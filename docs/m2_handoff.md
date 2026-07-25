# M2 handoff — where to pick up

Written at the end of the session that built M2.1–M2.3d. `origin/main` is at the
merge of `m2/convection`. `make all openmp=1` is clean; all **14** acceptance
tests pass. Full detail for everything below is in [`m2_results.md`](m2_results.md);
this is the orientation and the next step.

## What landed this session

- **M2.1 mass fixer** — global `p_s` rescaling, exact, off by default, on in
  Held–Suarez. Drives the §5.3 leak to machine zero. `p_s`-vs-`ln p_s` left open.
- **M2.2 correction framework** — `aeros_correction`, the §3.7/§3.8 pluggable
  additive tendency correction, at the **spectral** seam (it is linear, so
  unlike physics it need not go on the grid). No-op term, `ΔF=0` twin bit-exact.
- **M2.3a transport** — humidity is now a **gridpoint** prognostic (`spec%qv`
  removed), advected by positive-definite finite-volume transport off the
  spectral core. Conservation/positivity/constancy all machine-precision.
- **M2.3b condensation** — saturation adjustment, grid seam, drying to `qv_g`
  and latent heating to the spectral `T` tendency; column MSE closed by
  construction.
- **M2.3c van Leer limiter** — second-order horizontal transport (91.7% peak
  retention vs ~40% upwind). Corrected a wrong claim: the limiter does **not**
  shrink the moist water closure (that is the FV-vs-spectral air-mass gap, which
  scales with resolution, not tracer accuracy).
- **M2.3d convection** — Manabe operator, correct and unit-tested, but **off by
  default** because it is not yet usable coupled (see the next step).

Also measured on the running core: the thread cliff does **not** reproduce
locally (§1), T63 costs 3.27× T42/yr (§2), precision (dp) and the batched
transform (don't build) are closed (§3–4).

## The immediate next step: make convection usable coupled

`aeros_convection` is correct in isolation (`test_convection` passes: MSE
9.4×10⁻¹⁶, precip = water removed, neutralizes, dry branch exact, q≥0) but
`convect=.FALSE.` ships because enabling it has **two** problems, both design
decisions, detailed in `m2_results.md` §11 and the module header:

1. **Slow convergence for deep convection.** Pairwise Manabe propagates
   instability one layer per sweep (`MAXSWEEP=80` gives ~0.015 K residual;
   machine-neutral needs ~400). A cold-start moist column convects deeply in
   most columns → large cost. **Fix: multi-layer segment adjustment** —
   neutralize a whole connected unstable saturated segment at once. The algebra
   is in the header: `h_ref = (C0 + ΣΦ_i dp_i)/Σdp_i`, then a 1-D solve per
   layer `cp T_i + L q_sat(T_i) = h_ref − Φ_i`. Converges in a few passes.

2. **The coupled leapfrog goes unstable.** The convective heating is large and
   sign-alternating in the vertical; through the centered leapfrog RHS it
   excites the computational mode and the run NaNs in tens of steps (condensation
   doesn't — its heating is small/smooth). **Fix: change how moist physics
   couples** — apply the tendency as a forward adjustment on the `n+1` state
   (decoupled from the centered step), or give convection a finite relaxation
   timescale. The latter is literally the Betts–Miller direction, so this
   decision interacts with the BM-vs-Manabe question.

Recommended order: do (2) first (it blocks *any* coupled convection and also
matters for how strong a condensation heating the model can take), then (1).
When both are done, enable `convect` in `test_moist_run` (the coupled smoke test
already has the hook, currently commented) and add a coupled-convection stability
assertion.

**A design question to settle first (ask AR):** whether to fix Manabe's coupling
as above, or take (2)'s relaxation-timescale fix as the cue to implement
Betts–Miller now (the `scheme` selector and the whole seam are already built for
it — it is a second branch in `adjust_column` plus a reference-profile routine).
AR chose Manabe-first explicitly, so confirm before pivoting.

## Other open M2 items (not blocking)

- **Radiation + ERA5 validation.** The whole moist line has been validated by
  *conservation*, not climate — by design, until there is radiation and ERA5.
  AR said ERA5 data is coming. This is the big remaining M2 piece (design §5:
  bespoke ecCKD gas-optics table preferred; SESAM scheme is the pragmatic
  fallback). Call radiation every 1–3 h with per-step surface-flux rescaling.
- **Water closure characterization.** The FV-vs-spectral air-mass gap is
  ~2–3×10⁻³ at T21L12 (§8–9). It should be checked at T42 to confirm it shrinks
  with resolution as claimed, and checked over a long run to confirm it does not
  accumulate secularly (the way the mass leak did).
- **Transport accuracy leftovers** (§10): vertical van Leer (horizontal is done),
  and per-row polar sub-cycling instead of global (a cost optimization).
- **Carried from M1:** `tau_diff`/∇⁶-vs-∇⁴ retune (jets run ~10% strong);
  `SHqst_to_spat` 11→10 transforms; the Albedo 16–128 thread sweep to close the
  cliff question beyond 10 threads.

## Gotchas the next session will want

- **Worktree discipline** (CLAUDE.md): use `git -C <worktree>` and absolute
  paths; a bare `git status` reports main, not the worktree.
- **Generated `Makefile` is gitignored**; `config/Makefile` is the configme
  template. After merging a branch that added a test target, refresh main's
  `Makefile` from the template (or copy the worktree's) or the new test won't
  build on a fresh checkout. This bit twice already.
- **`timeout` is not on macOS** (it's `gtimeout`); don't wrap test runs in it.
- **Physics seam**: HS forcing, convection, condensation all act at the grid
  seam between `aeros_tendency_grid` and `aeros_tendency_spectral`, adding to
  `wrk%dtdt`; the correction acts at the spectral seam after
  `aeros_tendency_spectral`. Humidity (`qv_g`) is gridpoint and advected in step
  7 of `aeros_timestep_step`, using the time-`n` winds still in `wrk`.
- **Moist IC trap**: seeding humidity as a fraction of `q_sat` against an
  isothermal-warm column injects absurd values aloft (`q_sat` explodes at low
  `p`); use a realistic lapse-rate profile. Cost a debugging cycle already.
