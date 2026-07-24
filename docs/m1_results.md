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

## 5. What M1.5 needs from this

1. **Held–Suarez needs a non-zero `p_top`.** It is nominally defined on
   σ ∈ [0,1], but a zero-pressure top makes the top layer hydrostatically
   inconsistent — Simmons & Burridge impose `α₁ = ln 2` there rather than
   deriving it, and `test_tendency` measures the residual (2.4×10⁻¹⁰ at k = 1
   against 10⁻²² elsewhere). Use `p_top = 100`–`1000 Pa` and say so, or expect
   top-layer noise. `stretch_a = 1.0`, `sigma_t = 0.0`, `p_top = 1000.0` gives
   the 20 evenly spaced levels the benchmark asks for, with a consistent top.
2. **`tau_diff` is a tuning parameter and 6 h is a starting value, not a
   result.** Retune when the truncation changes, against the enstrophy
   spectrum rather than the look of a single field.
3. **Take the target values from Held & Suarez (1994) directly**, not from
   recollection.
4. The initial condition to perturb is `init_isothermal` in `src/aeros.f90`.

## 6. Still open, carried forward from M0a

1. **The OpenMP thread sweep still needs a real machine.** Extend `threads` in
   `par/bench_m0a.nml` to 16/32/64 on the HPC. That measurement is what
   design.md §3.6 and §9 risk 2 actually turn on, and it is still the largest
   unresolved engineering question in the project.
2. **`SHqst_to_spat` could take 11 transforms down to 10** by folding ζ in with
   (u,v). There is now a running core to measure it against.
3. Cost figures: the dry core needs **11 transforms per level, not 8**, so
   T42L20 is **19%** of §3.6's coupled budget, not the 14% first reported. See
   `docs/m0a_results.md` §4.2 for both columns.
