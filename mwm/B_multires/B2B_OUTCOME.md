# Stage b2 approach B (nudging / mean-tendency correction) — outcome

**Status: machinery correct, diagnosis run numerically fragile. Not a working
end-to-end result. Parked deliberately — see below.**

## What b2b does
Approach A (committed, `b2_probe.jl`) diagnosed ΔF from a *single* dynamical-core
tendency evaluation at the T85 time-mean X85 — insufficient (~19% Z* reduction,
T*/U* worse) because it ignores transient (eddy) rectification and physics response.

b2b diagnoses the correction via the CAPT **nudging tendency**: nudge T31 toward
X85 with a relaxation −(x−X85)/τ, average the applied nudging tendency ⟨N⟩ over the
run, then apply ⟨N⟩ as a constant forcing to a *free* T31 and test reproduction of
the T85 stationary waves. Averaging ⟨N⟩ over the fluctuating trajectory (rather than
evaluating at the mean) is meant to capture rectification + physics.

## What was fixed
A real SpeedyWeather **tendency-scaling bug**, resolved from the source (the `drag`
scheme + `scale_tendencies!`), not guessed. The correct per-field radius scaling for
injecting a tendency correction in SpeedyWeather v0.21.1 is:
- **spectral** vor/div: coefficient ×`scale[]` (=R) on the already-R-scaled state ⇒ R²·physical;
- **spectral** pres: ×R on the unscaled state ⇒ R¹·physical;
- **grid** temp/humid: add in **physical** units — `scale_tendencies!` multiplies the
  grid temp/humid tendencies by R *after* `forcing!` runs.

After the fix, the injection check is clean for every field (residual_frac
2e-7 … 0.046), i.e. ⟨N⟩ is applied correctly. This scaling knowledge transfers
directly to how aeros should inject its own correction terms.

## The remaining blocker (not a bug)
The nudging **diagnosis run itself is numerically unstable**: nudging a coarse
(T31/T21) spectral model toward the **unbalanced** high-res time-mean X85 excites
gravity waves and slowly blows up (NaN ~step 31 of the averaging window at the T21
smoke scale). This is a genuine methodological fragility — the concrete form of the
"a constant/strong correction is delicate" risk — not a coding error.

Mitigations exist but each is a **methodology choice best made in the real,
controlled aeros model** (where physics/hyperdiffusion are identical across
resolutions and there is no SpeedyWeather≠aeros caveat):
- nudge only **large scales** (spectral filter ⟨N⟩ / the relaxation) — the
  stationary-wave correction lives at large scales anyway (§3.7 "selected terms");
- nudge only the **slow, balanced variables** (e.g. vorticity + temperature; let
  divergence/pressure adjust freely) — the `par/b2b_vordiv.toml` variant is a start;
- gentler τ, or a balanced-state initialization;
- diagnose ⟨N⟩ from **−⟨RHS⟩** (approach-A's `eval_dynamics` path) accumulated over a
  lightly-constrained trajectory — guaranteed-correct units, decoupled from the nudge.

## Decision (2026-07-23)
Pivot to building aeros correction-ready rather than continue hardening stability in
SpeedyWeather. Even a success here would be only *suggestive* (SpeedyWeather ≠ aeros).
The definitive test is design milestone M2b in the controlled model. See
`../../docs/M1_scope.md`.

To reproduce: `julia --project=. b2b_probe.jl par/b2b_smoke.toml` (tiny, will NaN in
the nudged averaging phase as described).
