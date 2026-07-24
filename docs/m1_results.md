# M1 — the dry dynamical core, measured

What the core does now that it integrates. Everything below is from
`tests/test_semiimp.f90` and `tests/test_timestep.f90` at **T21 L20**, shipped
hybrid coordinate (`sigma_t = 0.2`, `p_top = 1000 Pa`), `dt = 1800 s`,
gfortran -O2 on an M5 MacBook. Reproduce with `make tests`.

---

## 1. The semi-implicit solve buys a factor of 3–6

| | |
|---|---|
| fastest linear gravity wave (external mode) | **335.2 m s⁻¹** |
| analytic `sqrt(γ R T_ref)`, T_ref = 300 K | 347.2 m s⁻¹ |
| explicit stability limit at T21 | **884 s** |
| explicit stability limit at T31 | **604 s** |
| configured step | 1800 s |

The 3.4% deficit against the analytic Lamb speed is the ~1% of atmospheric mass
above a 10 hPa top plus the L20 discretization of the mode. The speed is
**coordinate-independent to 7×10⁻⁴** between pure σ and the hybrid, which is a
sharper check than either number against theory — a coordinate-dependent error
in the linear operator cancels in neither.

Measured stability at `dt = 1800 s`, 200 steps from a moving state:

| scheme | peak growth of max\|ζ\| | diverged |
|---|---|---|
| semi-implicit, dt = 1800 s | 1.19 | no |
| explicit, dt = 1800 s | 2.1×10⁴ | **yes, by step ~24** |
| explicit, dt = 300 s | 1.00 | no |

The third row is what makes the second one evidence rather than a broken code
path.

## 2. The linearization is not the obvious one

The coefficient of `ln p_s` in the implicit divergence equation is
**`R_d T_ref`, not `R_d T_ref c_k`** — even though the nonlinear
pressure-gradient force in `aeros_tendency` carries `c_k` (its `cfac`), so
matching it looks like the correct thing to do.

A surface-pressure perturbation also moves the half-level pressures, hence
`α_k` and `Δln p_j`, so the geopotential responds too:

```
dΦ_k/d(ln p_s) = R_d T_ref (1 − c_k)      (isothermal column)
```

which is exactly the identity Simmons & Burridge's `α_k` is constructed to
satisfy — the one `test_tendency` measures at 3×10⁻¹². Added to the
`−R_d T_ref c_k` from the pressure-gradient term it leaves `R_d T_ref`.

Restoring `c_k` makes the operator wrong by `R_d T_ref (1 − c_k)`, which in the
upper layers of a hybrid coordinate is **100% of the term** (`c_k → 0` there).
`test_semiimp` finite-differences the real nonlinear tendency against the
operator and is what holds this down: agreement is 2×10⁻¹⁰ for a temperature
perturbation and 2×10⁻⁹ for a surface-pressure one.

## 3. Diffusion damps per **two** steps, not per step

Leapfrog carries two interleaved chains, and the implicit ∇⁶ factor is applied
over the full `2 dt` from `X^(n−1)` to `X^(n+1)`. So within one chain the
amplitude falls by `1/(1 + 2 dt/τ)` every **two** steps, and the apparent
per-step decay is the square root of that.

Measured at `l = lmax`, `tau_diff = 6 h`: **0.85206 per 2 dt** against the
analytic 0.85714 — 0.6% low, which is the linear coupling to divergence that a
vorticity perturbation cannot avoid exciting. Comparing the per-*step* ratio to
the per-application factor is the natural mistake and is wrong by a factor of
two in the exponent.

## 4. Conservation

Two states, because they answer different questions.

**Balanced resting atmosphere over topography, 100 steps.** Nothing here can
move except a genuine leak in the integrator:

| | drift over 100 steps | per year |
|---|---|---|
| mass | 2.0×10⁻¹⁶ | −3.5×10⁻¹⁴ |
| energy | 1.9×10⁻¹⁶ | −3.2×10⁻¹⁴ |
| angular momentum | 4.3×10⁻¹⁶ | 7.6×10⁻¹⁴ |

All three at machine precision. The state also holds to **3×10⁻¹¹ m s⁻¹** of
spurious wind and 7×10⁻¹⁰ Pa of surface pressure. The standalone driver
confirms it over 1753 steps at T31L16: surface pressure exactly 10⁵ Pa,
temperature exactly 300 K, winds ~10⁻²⁷ m s⁻¹.

**Moving state, 400 steps.** Mass is a bounded oscillation of order 10⁻⁵ that
crosses zero repeatedly and shows no secular trend — the discrete-time
signature of a gravity-wave adjustment, since what the scheme conserves exactly
is `d/dt ∫p_s dA`, not `∫p_s dA` itself. Energy drifts at 3×10⁻⁵ over the
integration and angular momentum at 3×10⁻⁴.

**Energy and angular momentum are reported, not asserted.** Leapfrog, the time
filter and a deliberately dissipative hyperdiffusion all remove energy; no
honest implementation of this scheme conserves it to machine precision, and
what a paleo integration needs is the drift rate rather than a claim. The rest
state is where the assertion belongs.

## 5. Held–Suarez (M1.5)

`make held-suarez && libaeros/bin/held_suarez.x par/held_suarez.nml`. 1200 days,
the last 1000 averaged. ~2 min at T31L20 and ~6 min at T42L20 on ten cores, so
this is a routine regression run rather than a campaign.

**There are no published target numbers.** Held & Suarez present figures and ask
that yours look like theirs; CESM, ECMWF and MITgcm all restate it that way, and
the AMS page is paywalled. So the numbers below are aeros' own, and their role
is (a) to be compared by eye against the published figures — see
[docs/figures/held_suarez.png](figures/held_suarez.png), reproducible with
`tools/hs_report.py` — and (b) to be the regression reference for later changes.

| | T31L20 | T42L20 |
|---|---|---|
| jet max `[u]`, NH | **34.6 m s⁻¹** at 46.4°, 233 hPa | **32.4 m s⁻¹** at 46.0°, 232 hPa |
| jet max `[u]`, SH | 35.3 m s⁻¹ at −42.7°, 233 hPa | 32.9 m s⁻¹ at −43.3°, 233 hPa |
| surface `[u]` max | 9.7 m s⁻¹ at 46.4° | 8.8 m s⁻¹ at 46.0° |
| surface `[u]` tropics | −8.6 m s⁻¹ | −8.7 m s⁻¹ |
| max `[v'T']` | 18.7 K m s⁻¹ at 39.0°, 883 hPa | 21.1 K m s⁻¹ at 37.7°, 882 hPa |
| max `[u'v']` | 74.4 m² s⁻² at 35.3°, 235 hPa | 75.7 m² s⁻² at 34.9°, 235 hPa |
| max eddy KE | 268 m² s⁻² | **367 m² s⁻²** |
| max `[T'²]` | 31.3 K² | 40.4 K² |
| surface `[T]`, equator−pole | 306.8 − 263.0 = 43.8 K | 305.9 − 263.8 = 42.2 K |
| hemispheric asymmetry of `[u]` | 0.086 | 0.084 |

Structurally this is the published circulation: two westerly jets near ±45° at
~250 hPa, tropical and polar surface easterlies with midlatitude surface
westerlies, an isothermal 200 K equatorial cap aloft where the `T_min` floor
binds, `[v'T']` poleward in both hemispheres with its maximum in the lower
troposphere near 850 hPa, and `[u'v']` converging into each jet from the
upper-troposphere maxima near ±35°. All the signs are right:

```
              [v'T']    [u'v']        (means over 20-70 deg)
    NH        +7.0      +12.2         poleward heat, momentum into the jet
    SH        -7.0      -11.0
```

and the two hemispheres agree to ~1% in the fluxes and ~2% in the jet, from a
forcing that is symmetric by construction and an integration that was seeded
asymmetrically.

### 5.1 T31 is missing a quarter of the eddy activity

The one number that moves substantially between truncations is **eddy kinetic
energy: 268 m² s⁻² at T31 against 367 at T42, i.e. T31 is 27% low**. Eddy
temperature variance moves the same way (31.3 vs 40.4 K²) and the eddy heat flux
with it (18.7 vs 21.1 K m s⁻¹). The jet responds as one would expect from
weaker eddies: **T31's is ~7% STRONGER** (34.6 vs 32.4 m s⁻¹), there being less
eddy momentum flux to decelerate it.

This is the resolution bind of design.md §3.6 and §9 risk 1, now measured on
aeros' own core rather than inferred: at T31 the storm track is under-resolved
in amplitude while the jet is over-strong, and both of those set where
precipitation lands on an ice sheet. It does not decide the truncation — §3.7's
correction is what M2b tests, and the jet LATITUDE is nearly identical at both
(46.4° vs 46.0°) — but it puts a number on what bare T31 costs.

### 5.2 Diffusion is not part of the benchmark, and it matters

Held & Suarez deliberately leave the horizontal diffusion to the modeller, since
it is part of what an intercomparison compares. aeros brings ∇⁶ with a 6 h
e-folding at `l = lmax`; Held & Suarez' own spectral model used ∇⁴, and CESM's
Held–Suarez configuration uses ∇⁴ with a 0.5 day timescale. aeros' jets run
~10% stronger than the ~30 m s⁻¹ usually read off the published figure, and a
more scale-selective, weaker diffusion is the first thing to suspect. Any
comparison against a published figure has to state which diffusion was used.

### 5.3 The mass leak — the real finding

**Global mass drifts linearly at ~6.6×10⁻⁶ per year**, and this is the one
result from M1.5 that should change what happens next.

The adiabatic tests (§4) show mass conserved to 2×10⁻¹⁶ on a resting state and
bounded within ~10⁻⁵ on a moving one. Under Held–Suarez, in a statistically
steady eddying state, it becomes a **secular** loss. T31L20 at `dt = 1800 s`:

```
    day  300     600     900    1200
    dM/M  -5.1e-6 -1.1e-5 -1.8e-5 -2.4e-5      increments -6.3, -6.5, -6.2e-6
```

Perfectly linear. Extrapolated to design.md §1's 10⁵ yr target that is **0.66 of
the atmosphere**, so it is not a rounding detail.

It is not a bug and not the time filter. Measured at day 400, T31:

| configuration | dM/M |
|---|---|
| ν = 0.06, dt = 1800 s | −7.24×10⁻⁶ |
| ν = 0.15, dt = 1800 s | −7.75×10⁻⁶ |
| ν = 0.06, dt = 900 s | −3.30×10⁻⁶ |

**2.5× the filter strength changes the drift by 7%; halving the timestep halves
it.** So the drift rate is ∝ `dt`, i.e. the per-step error is the scheme's
ordinary O(dt²) time-truncation error, and what makes it accumulate rather than
cancel is that the conserved quantity is a NONLINEAR functional of a prognostic:
the model integrates `ln p_s`, and mass is `∫exp(ln p_s) dA`. What the
discretization conserves exactly is `∫ p_s d(ln p_s)/dt dA = 0` — measured at
1.4×10⁻¹⁶ in `test_tendency` — which is the continuous statement, not the
discrete one. In an eddying state those per-step errors have a systematic sign.

Running with the filter disabled entirely is not an option for comparison: the
leapfrog computational mode grows and the run is NaN within 100 days.

Three ways out, none of them implemented, all of them a decision rather than a
fix:

1. **A global mass fixer** — rescale `p_s` each step so the global integral is
   restored. One line, standard practice in spectral models, and exact. It puts
   a small non-local correction into the surface pressure.
2. **Carry `p_s` rather than `ln p_s`.** Removes the nonlinearity at the root,
   and costs the reason `ln p_s` was chosen: the pressure-gradient term is
   linear in `ln p_s`.
3. **Accept it and shorten `dt`.** Since the rate is ∝ `dt`, halving the step
   halves the drift — at twice the cost, and it never reaches zero.

Note that **T42 leaks less than T31** (−1.46×10⁻⁵ against −2.42×10⁻⁵ at day
1050), so the rate depends on the flow rather than on the truncation, and
resolution is not a way out either.

## 6. What is still open

### 6.1 Carried forward from M0a

1. **The OpenMP thread sweep still needs a real machine.** Extend `threads` in
   `par/bench_m0a.nml` to 16/32/64 on the HPC. That measurement is what
   design.md §3.6 and §9 risk 2 actually turn on, and it is still the largest
   unresolved engineering question in the project.
2. **`SHqst_to_spat` could take 11 transforms down to 10** by folding ζ in with
   (u,v). There is now a running core to measure it against.
3. Cost figures: the dry core needs **11 transforms per level, not 8**, so
   T42L20 is **19%** of §3.6's coupled budget, not the 14% first reported. See
   `docs/m0a_results.md` §4.2 for both columns.
