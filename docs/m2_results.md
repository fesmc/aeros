# M2 — measured on the running core

Results produced against the real dynamical core, as opposed to the M0a proxy
harness (`drivers/bench_m0a.F90`) or the separate `mwm/A_scaling` benchmark.
Both of those stood in for a core that did not exist yet; this one does.

## 1. The thread cliff does not reproduce

`docs/m1_results.md` §6.1 carried an open question inherited from the merged
`mwm/` workstream: `mwm/A_scaling/RESULTS.md` measured the SHTns-per-field
transform saturating at 4 threads and then *collapsing* — 10–20× slower than
one thread at 64–128 — and the question was whether aeros' own wrapper
inherits that.

**Within the range that can be tested locally, it does not.**

Method: time `libaeros/bin/held_suarez.x` — the real core, so the real 11
transforms per level, the real semi-implicit solve, and the Held–Suarez forcing
at the grid seam — at each thread count, differencing a 100-day run against a
20-day one so that SHTns table setup, netCDF open and the initial snapshot
cancel out rather than being amortized over an assumed step count. L20
throughout. Raw data in the log; the sweep script is reproducible from it.

| thr | T31 ms/step | ×    | T42 ms/step | ×    | T63 ms/step | ×    |
|----:|------------:|-----:|------------:|-----:|------------:|-----:|
| 1   | 5.393       | 1.00 | 14.164      | 1.00 | 30.917      | 1.00 |
| 2   | 3.628       | 1.49 | 8.203       | 1.73 | 17.212      | 1.80 |
| 4   | 1.763       | 3.06 | 4.279       | 3.31 | 10.375      | 2.98 |
| 8   | 1.792       | 3.01 | 3.854       | 3.67 | 8.047       | 3.84 |
| 10  | 1.776       | 3.04 | 3.669       | 3.86 | 7.203       | 4.29 |

Against `mwm/A_scaling` at the same thread counts, where T31 peaked at 1.97× on
4 threads and had already fallen to 1.20× by 8, and T42 peaked at 2.34× and
fell to 1.70×:

1. **No collapse anywhere in range.** T31 plateaus at ~3× from 4 threads on.
   T42 and T63 are still *improving* at 10 threads. The `mwm/` curves had
   turned over by 8 in both truncations.
2. **Scaling improves with truncation** — best speedup 3.06 (T31), 3.86 (T42),
   4.29 (T63). That is the same diagnosis `mwm/` reached for why SHTns stalls,
   read forwards: the constraint is too little work per Legendre transform, so
   more work per transform scales better. It also means the truncation the
   design most wants to afford is the one that parallelizes worst.
3. **Effective parallelism is 3–4.3, not ~24.** So `mwm/`'s structural
   conclusion survives intact even though its curve does not: the core wants
   ~4 threads, and design.md §3.6's N_eff ≈ 24 for bare T42 is out of reach.
   What changes is that running at 8–10 threads is no longer *harmful*.

**Caveats, and they are not small.** Ten cores cannot see the 16–128 range
where `mwm/`'s collapse was worst, so this rules out the collapse *starting* at
8 threads — where `mwm/` saw it begin — and nothing beyond that. This is
gfortran on macOS on a laptop against `ifx` on an AWI Albedo Xeon node, so the
comparison is directional rather than like-for-like, and the absolute yr/day
figures here are not the numbers for the target machine. The sweep should be
repeated on Albedo to 128 threads before the cliff is called closed.

## 2. T42 vs T63 cost

The benchmark `docs/M1_scope.md` lists in its definition of done.

| | T42L20 | T63L20 |
|---|---|---|
| dt | 1800 s | 1200 s |
| ms/step, 1 thread | 14.164 | 30.917 |
| core-s/yr, 1 thread | 248 | 813 |
| core-s/yr, 4 threads | 300 | 1091 |

T63 costs **3.27×** T42 per simulated year. That splits into 2.18× per step
and 1.5× from the shorter timestep — dt is not a free parameter here, it falls
with truncation on CFL grounds, and a cost comparison quoted per step rather
than per simulated year would understate T63 by a third.

Note 2.18× per step is well below the (63/42)³ = 3.375 a pure T³ scaling would
give. At these truncations the Legendre transform is not yet the asymptotically
dominant term, which is the same reason the thread scaling improves with
truncation.

## 3. Precision: sp is closed for the core, not for the physics

`docs/m1_results.md` §6.0 recorded a disagreement with the `mwm/` line over
precision — this line measuring dp ~17% faster than sp, the `mwm/` line
asserting "Float32 throughout the core" and sp ≈ 2× faster — and left it open
on the grounds that the two might have measured different transform paths.

**It is not a disagreement between two measurements. `mwm/` never measured
sp.** Every table in `mwm/A_scaling/RESULTS.md` is labelled double precision;
the 2× appears twice as a rule of thumb ("Float32 ≈ 2× faster again", "~4000
**est.** Float32") and never as a measured column. This line's number
(`docs/m0a_results.md` §5) is a measurement.

The physical reason stands and is not a matter of tuning: SHTns has no
single-precision CPU path (`SHT_FP32` is GPU-only), so an sp core converts up
and back at every transform, and transforms are 86–97% of a step. **dp for the
core, decided.**

Two consequences worth keeping straight:

- Adopting a different transform would reopen it. §4 below closes that route,
  so in practice it stays shut.
- sp remains genuinely live **inside column physics**, which never touches
  SHTns. `src/aeros_defs.f90` already says a physics module may declare it
  locally. That is a different decision from the model's working precision and
  it should be made with a measurement at M2, when there is physics to measure.

## 4. The batched-DGEMM transform is not worth building

`mwm/A_scaling/RESULTS.md` §"Bottom line" already reaches this conclusion on
its own data, and §1 above does not disturb it:

- batched does not raise the ceiling — T42 best 1997 yr/day against SHTns'
  2053, T31 4068 against 3986;
- it is ~1.7–1.8× slower on one thread;
- at equal throughput it costs 4× the cores (693 core-s/yr against 168).

Its one advantage is not collapsing past 8 threads, and §1 finds aeros does not
collapse there either. Under an ensemble-parallel throughput model —
~4 threads per member, the node's parallelism coming from member count — the
batched transform is strictly worse. **Keep it documented as a fallback for the
case where a single member must be pushed past ~4 cores, and do not build it.**

## 5. What "ensemble parallelism" does and does not buy

Worth stating plainly, because the phrase does real work in `mwm/`'s
conclusions and can be read as solving more than it does.

Ensemble parallelism does **not** make one run faster. For a single 10⁵ yr
transient coupled to an ice sheet — the actual target of design.md §1 — the
core wants ~4 threads and no number of additional cores changes its wallclock.
That is a constraint on the design, and §1's measurement does not lift it; it
only removes the penalty for overshooting.

Where genuine ensembles do appear: the `ΔF = 0` twin that §3.8 mitigation 3
requires of every production run (2×), the reference-dataset spread of
mitigation 4 (3×), and any parameter-calibration campaign. Those are throughput
for a *campaign*, not wallclock for a run, and they are what the remaining
~124 cores of a node are for.

One consequence for the throughput target: `mwm/`'s route to design.md §3.5's
5000 yr/day ran "~4000 double → ~8000 Float32 → coupled 4000–5000 plausible".
§3 above removes the Float32 factor from that chain, and it was load-bearing.
The target needs re-deriving against measured numbers on the target machine
rather than inheriting.

## 6. The mass fixer works, and reports a number that is easy to misread

`docs/m1_results.md` §5.3 left three options for the ~6.6×10⁻⁶/yr mass leak.
The global fixer is implemented (`mass_fixer`, off by default, on in
`par/held_suarez.nml`); carrying `p_s` instead of `ln p_s` remains open and is
deliberately not foreclosed — the switch makes the two comparable, since
`mass_fixer = .FALSE.` recovers the unfixed integrator exactly.

**It works.** Held–Suarez T31L20, 400 days, fixer on: mass drift
**−1.6×10⁻¹⁵ per year**, against **−6.608×10⁻⁶** for the same configuration
unfixed — which itself reproduces §5.3's reported figure exactly. The
correction is closed-form rather than iterative (multiplying `p_s` by a
constant multiplies `∫p_s dA` by exactly that constant), so it touches the
(0,0) coefficient and nothing else. `tests/test_timestep.f90` asserts that
bitwise: vorticity, divergence, temperature and every `lnps` wavenumber but
(0,0) are unchanged to the last bit after a step.

**And it reports a number that must not be quoted as the drift.** The fixer
accumulates every correction, and the temptation is to read that total as
"what the unfixed model would have lost". It is not, and it is consistently
about twice as large:

| | unfixed drift | fixer put back | ratio |
|---|---|---|---|
| adiabatic moving state, 200 steps | 7.45×10⁻⁶ | 1.91×10⁻⁵ | 2.56 |
| Held–Suarez T31L20, 400 days | −6.61×10⁻⁶/yr | −1.32×10⁻⁵/yr | 1.99 |

Two candidate explanations, and the measurements separate them. It is **not**
trajectory divergence: four unfixed Held–Suarez runs seeded with perturbations
differing in the eighth decimal give −6.608, −6.627, −6.664 and −6.668×10⁻⁶/yr,
a spread under 1%, and the adiabatic moving state is not chaotic at all yet
shows the same factor.

What it is: the unfixed run's mass carries a **bounded oscillation** as well as
a secular drift — the leapfrog computational mode in mass, already documented
in `m1_results.md` §4 as the ~10⁻⁵ bounded excursion on a moving state. Over
many steps that oscillation largely cancels in an accumulated total. The fixer
resets it every step and therefore pays for it every step. Per step the two
agree exactly — `test_timestep` pins the reported correction to the step's
actual leak at 10⁻¹⁴, eight decades below the leak itself — but the sum is the
fixer's *workload*, not a twin's drift.

So the cumulative figure is worth watching (a fixer whose workload grows is
covering for something) and worthless as a drift measurement. Measuring the
drift needs an unfixed twin run — the same discipline design.md §3.8 already
requires of the correction terms, arrived at here for an unrelated reason.

## 7. The tendency-correction framework

`docs/M1_scope.md` lists a "first-class, pluggable correction framework" in its
definition of done, and `m1_results.md` §6.0 recorded that it had not been
built. It is now (`src/dynamics/aeros_correction.f90`), carrying a validated
no-op term.

**The seam is spectral, not the grid seam the physics uses.** That is a
deliberate split rather than a convenience. Physics must be applied on the grid
because it is nonlinear in the state — a horizontally varying drag does not
commute with the curl, so `−k_v·v` is not `−k_v·ζ` and `−k_v·D`. A correction
is a different object: additive and linear, so it commutes with everything, and
applying it spectrally is free of that constraint. What that buys:

1. **The scale selector costs nothing.** §3.7 wants the correction applied only
   at large scales, and `B2B_OUTCOME.md` lists spectral filtering as the first
   mitigation for the instability its whole-field version hit. Spectrally that
   is a mask on total wavenumber; in grid space it would be a transform
   round-trip per field per step.
2. Correction and physics stay independently switchable — different claims
   about the model, so they should not share a seam.
3. It matches how `ΔF` is diagnosed: a coarse-grained difference of two models'
   tendencies (§3.7 step 4).

The cost, stated plainly: a correction that genuinely had to be nonlinear in
grid space could not use this path. Neither §3.7 nor §3.8 describes one.

**Units.** Correction terms are stored in the same units as the tendency they
are added to, with no scaling factor: `vor`, `div` [s⁻²], `temp` [K s⁻¹],
`lnps` [s⁻¹]. This is aeros' analogue of the convention `mwm/B_multires` had to
recover from SpeedyWeather's source, where spectral vor/div needed ×R², pres
×R¹, and grid temp/humid physical units — and getting it wrong produced a
silently wrong-magnitude correction. aeros carries unscaled physical tendencies
throughout, so the convention is "no convention", which is worth writing down
precisely because the last model needed a page for its own.

**What the tests establish** (`tests/test_correction.f90`, ten checks):

- the `ΔF = 0` twin is **bit-exact** over a 20-step integration with a term
  that is enabled *and filled* — so the switch is consulted before the addition
  and not merely per-term;
- a registered but unfilled term is **inert to the last bit** over the same
  integration, which is what separates "the plumbing is in" from "a correction
  is being applied";
- a filled term adds exactly what it holds, asserted as `result == ref + corr`
  rather than `result − ref − corr == 0` (the latter is a claim about
  floating-point reversibility, not about the code, and fails at ~10⁻²¹ purely
  from the magnitude ratio);
- the scale selector cuts exactly at `lcut`, checked in **both** directions —
  168 coefficients corrected and 1856 untouched at T21 with `lcut = 5` — since
  a selector that silently passed everything would look like a working
  correction;
- per-field selection leaves the unselected fields exactly as the dynamics left
  them.

**What is deliberately not here.** Diagnosing `ΔF` is M2b — it needs two
truncations under identical boundary conditions and a coarse-graining operator,
none of which exists yet. And full flux-form conservation verification (§3.7
risk 2, "verify to machine precision") is not built against a term that is
identically zero: `aeros_correction_report` measures the channel a correction
would inject mass through — the (0,0) coefficient of the `lnps` correction, the
only mode that moves the global mean of `ln p_s` — but a complete statement
needs the grid-space integral `∫ p_s ∂(ln p_s)/∂t dA`, because mass is a
nonlinear functional of the prognostic. That is the same nonlinearity §6 above
is about, and it belongs with the first real `ΔF`.

## 8. Prognostic humidity: positive-definite transport (M2.3a)

`qv` was allocated and carried since M1 but never evolved. It now is — and the
first decision it forced was structural: **humidity is a gridpoint prognostic,
off the spectral core entirely.** `spec%qv` is gone. Advecting a positive field
spectrally is not positive — truncating the exp-shaped gradients at fronts and
the ITCZ produces Gibbs ringing and negative water (design.md §4.2) — so q is
carried in `now%qv_g` and advected by a finite-volume scheme on the Gaussian
grid that never touches the transform.

**The scheme** (`src/dynamics/aeros_moisture.f90`) is flux-form, first-order
upwind, forward in time, with its own finite-volume air-mass budget — the
"option (ii)" chosen over leapfrogging q with a water mass-fixer. It advances
the tracer mass `q·Δp` and a finite-volume air mass by the *same* face fluxes
and divides (Lin–Rood), which is what makes the three required properties hold
as machine-precision facts rather than approximations. `tests/test_moisture.f90`
measures each:

| property | measured | why it matters |
|---|---|---|
| constancy (q≡const preserved, divergent wind) | 2.6×10⁻¹⁵ | the first thing a tracer scheme breaks; needs air and tracer on identical fluxes |
| conservation (total water, 100 steps) | 3.0×10⁻¹⁵ | flux form + single-valued faces; poles carry zero flux since cos(lat)=0 |
| positivity (sharp blob, 200 steps) | min q ≥ 0 | the property spectral advection cannot deliver |
| polar sub-step fires + still conserves | 11 sub-steps, 1.2×10⁻¹⁵ | the Gaussian grid's zonal Courant blows up at the poles |

**Mass consistency — the price of positivity on a spectral core.** The scheme's
FV air mass differs from the spectral air mass by O(truncation), because one
divergence is a grid finite-difference and the other a spectral derivative. This
is not a leak: the layer masses are re-diagnosed from the true (spectral)
surface pressure every step, so it is a bounded per-step consistency error, and
the max-principle division keeps q bounded regardless. It is diagnosable, and
the unit conservation test isolates it — a conservation check with a *fixed* ps
and a wind that is divergent against it fails by ~10⁻², precisely that gap, so
the machine-precision checks use a uniform ps where FV and spectral air mass
coincide. The interaction with an evolving ps is an integration-level check, for
when condensation and a moist Held–Suarez run exist.

**The polar Courant problem.** A spectral core has no grid CFL; a grid transport
does, and on a Gaussian grid the zonal cell width collapses toward the poles
while the wind does not, so the zonal Courant number reaches ~9 at T21. The fix
is to sub-step the whole transport, `nsub` set from the largest Courant on the
grid — *not* a Fourier/polar filter, which is the classical spectral remedy and
is exactly what must not be used, since filtering a positive field in wavenumber
space makes it negative. Sub-stepping preserves every invariant; it only costs
arithmetic, and transport is grid-local and cheap against the transforms.

**Deliberately deferred, and said so in the module header:**

- *Accuracy.* The fluxes are first-order upwind — the correct-and-provable
  baseline, whose conservation/positivity/constancy are unconditional, but too
  diffusive for a real humidity field. A van Leer limiter is the next commit
  (§10) and keeps the invariants (they depend on monotonicity, not order).
- *The dry path.* A dry run still sub-steps the transport of identically-zero
  humidity every step — correct (0 stays 0) but wasted work. A moisture-active
  switch belongs with condensation, where there is a natural place for it.
- *Per-row sub-cycling.* Global sub-stepping makes the whole grid sub-step at
  the polar rate; a per-row zonal sub-cycle would confine the extra work to the
  rows that need it. An optimization, not a correctness item.

Condensation — supersaturation removal, latent heating, precipitation — is the
next commit.

## 9. Large-scale condensation (M2.3b)

The first moist physics: where the air is supersaturated, the excess vapour
condenses, the latent heat warms the air, and the condensate rains out. No
stored cloud water, no re-evaporation, no convection (next), no ice phase
(condensate is liquid, L_v, everywhere) — the minimum that closes a moist energy
and water budget, built to be exactly that.

**Two seams, kept consistent.** Condensation couples the *gridpoint* humidity to
the *spectral* temperature, so it acts at both at once. The drying is applied
straight to `qv_g`; the heating, `L_v/cp_d · dq_c/dt`, is added to the
temperature tendency `wrk%dtdt` at the same grid seam as the Held–Suarez
forcing, so it rides the same transform and the same leapfrog as the dynamical
heating. Both come from the same condensed amount `dq_c`, so column moist static
energy is conserved by construction — dry static energy gains `L_v dq_c`, latent
energy loses it. `tests/test_condensation.f90` checks this as an equality, per
cell, at 1.5×10⁻¹⁵ (relative), alongside: subsaturated columns untouched
exactly, condensing columns brought to saturation (3×10⁻¹³), vapour removed =
precipitation (1.9×10⁻¹⁶), q ≥ 0. The saturation adjustment is three Newton
iterations of the implicit equation `dq_c = (q − q_sat(T))/(1 + (L_v/cp_d)
dq_sat/dT)`; `q_sat` is Tetens with the full `p − (1−ε)e_s` denominator (14.66
g/kg at 20 °C, 1000 hPa).

**The coupled run holds together** (`tests/test_moist_run.f90`, a 200-step moist
Held–Suarez at T21L12): 200 steps, no NaN, q ≥ 0 throughout, wind bounded, 3.3%
of the vapour rained out. This is the test the two operator unit tests cannot be
— that the latent heating actually reaches the spectral temperature and the
system stays stable — and it caught two real things in the building:

1. *Condensation off, the water closure is 2.4×10⁻⁴* over 200 steps — the
   finite-volume transport's O(truncation) dispersion against the spectral
   surface pressure, measured on a running model for the first time (§8
   predicted it; here it is).
2. *Condensation on, it rises to 2.1×10⁻³.* Not a leak, and — see §10 — not
   the transport's own conservation error, which is machine precision. It is
   the FV-vs-spectral air-mass gap (§8): the diagnostic weights q by the
   spectral layer masses, the transport conserves against its finite-volume
   ones, and the two differ by O(truncation). It scales with the roughness of
   q, which is why condensation (sharpening q) raises it from ~2×10⁻⁴. It
   shrinks with resolution, not with a better tracer scheme — the design does
   not permit machine-precision water conservation on a spectral core with
   positive-definite transport, which is the trade §8 named.

**The two seams are on different time discretizations** — the drying is a
forward step on `qv_g`, the heating goes through the 2 dt leapfrog and its
filter — so over a run they track to the filter's accuracy rather than exactly.
That is the same small inconsistency every leapfrog spectral model carries
between grid physics and spectral dynamics; it is inside the 2.1×10⁻³ above,
not separately asserted.

Convection — the moist adjustment that keeps the tropics from grid-scale
saturating — is the next commit.

## 10. The van Leer limiter (M2.3c)

The first-order upwind transport of §8 was the correct-and-provable baseline;
this replaces it with a second-order van Leer scheme, and the replacement is
structural rather than a one-line swap. Upwind needs only the two cells at a
face; a limited reconstruction needs a three-cell stencil in the flux
direction, and multidim positivity is cleanest via operator-split
one-dimensional sweeps. So `transport_substep` is now a zonal sweep, then a
meridional sweep, then the vertical sweep — each flux-form, each monotone, each
updating the air mass and the tracer mass by the same fluxes.

**The invariants are untouched**, exactly because they rest on monotonicity, not
on order: `test_moisture` still reports constancy 5.6×10⁻¹⁵, conservation
4.1×10⁻¹⁶, positivity, and conservation through the polar sub-step 4.7×10⁻¹⁶.
Van Leer is applied in the two horizontal sweeps, where humidity gradients and
Courant numbers both live; the vertical stays first-order upwind (a handful of
levels, a gentle mass-coordinate flux) — van Leer there is a later refinement.

**What it buys is accuracy, measured.** A broad smooth bell advected one full
solid-body rotation (371 steps) — for which the exact answer is the bell it
started as — retains **91.7% of its peak amplitude**. First-order upwind smears
the same bell well below half over that distance. That is the whole point: a
humidity field with real fronts and an ITCZ edge survives being advected, rather
than diffusing to mush.

**What it does NOT buy — a correction to §9.** I expected the limiter to shrink
the moist water-budget closure (2.1×10⁻³). It does the opposite: 2.1×10⁻³ →
3.0×10⁻³. The closure is not tracer-transport error — the transport conserves
its own water to machine precision — it is the FV-vs-spectral air-mass gap (§8),
and that gap scales with the *roughness* of q. A less-diffused, more accurate q
has larger gradients, so the truncation-level gap acts on more structure and the
closure rises. The number improves with resolution, not with a better tracer
scheme. This is worth stating plainly because it was a wrong prediction in §9,
now corrected there: accuracy of the field and closure of the budget are
different axes, and the limiter is on the first.

## 11. Moist convective adjustment — operator done, coupling not (M2.3d)

Large-scale condensation only removes gridbox-mean supersaturation; without
convection a column heated from below saturates layer by layer and rains
stratiform everywhere. `aeros_convection` adds the Manabe (1965) moist
convective adjustment: sweep the column for unstable adjacent pairs and relax
each to neutrality — dry pairs mix potential temperature (conserving enthalpy),
saturated moist-unstable pairs relax to a moist adiabat (conserving
∫(cp T + L q) dp, both saturated) and precipitate. Manabe was chosen over
Betts–Miller because it conserves by construction — no timescale, no reference
RH, no energy-closure correction — and the module carries a `scheme` selector so
Betts–Miller is a second branch, not a rewrite.

**The operator is correct and unit-tested** (`tests/test_convection.f90`, all
passing): moist static energy conserved 9.4×10⁻¹⁶, precipitation equals the
water removed 5.5×10⁻¹⁵, the column is moist-stable afterward, the dry branch
conserves enthalpy exactly and moves no water, q stays ≥ 0, and a stable column
is left untouched.

**It is OFF by default and not yet usable coupled.** Enabling it in the full
model surfaced two real problems, neither a bug in the operator:

1. *Convergence is slow for deep convection.* Pairwise Manabe propagates
   instability one layer per sweep, so a deep cold-start column converges only
   geometrically — machine-neutral needs ~400 sweeps at L20, and even the
   physically-neutral ~0.015 K takes 80. In a running model each step's
   instability is a small increment (a handful of sweeps), but a cold-start
   moist atmosphere convects deeply in most columns and the cost is large. The
   fix is a **multi-layer adjustment** — neutralize a whole connected unstable
   segment at once (h* = const across it, conserving segment MSE, which reduces
   to one reference value plus a 1-D solve per layer), converging in a few
   passes. Documented in the module header.

2. *The coupled leapfrog goes unstable.* The convective heating is large and
   sign-alternating in the vertical (cool low, warm high), and pushing it
   through the centered leapfrog RHS every step excites the leapfrog
   computational mode — the run NaNs within tens of steps. Condensation does not
   do this because its heating is small and smooth. This is the textbook
   difficulty of stiff physics through leapfrog, and the fix is a physics-
   coupling change: apply the moist-physics tendency to the n+1 state as a
   forward adjustment (decoupled from the centered step), or give convection a
   finite relaxation timescale — which is, not coincidentally, the Betts–Miller
   direction.

So `convect = .FALSE.` ships; the dry benchmark and the condensation-coupled
moist run are unaffected and all fourteen tests pass. The two items above are
the next session's work, and both are design decisions (multi-layer adjustment;
how physics couples to the leapfrog) rather than mechanical fixes.
