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

So at M2.3d `convect = .FALSE.` shipped; the dry benchmark and the
condensation-coupled moist run were unaffected. **Both problems are now solved**
— problem 2 by the forward-split coupling, and problem 1 by switching the default
scheme to Simplified Betts–Miller, which has no convergence-to-neutrality loop at
all. See §12.

## 12. Convection made usable coupled — Simplified Betts–Miller (M2.3e–g)

Two commits of coupling and scheme work, plus the coupled test, close out §11's
two problems. Convection now runs in the full model without blowing up — the
first time. It is still `convect = .FALSE.` by default (validation waits on
radiation and ERA5, §13/design §5), but the machinery is done and tested.

### 12.1 Forward-split physics heating (M2.3e) — the coupling fix

The instability of §11 problem 2 was never in the operator; it was in how the
heating reached the spectral state. Physics heating rode `wrk%dtdt` through
`aeros_tendency_spectral` and into the **centered** leapfrog (`X^{n+1} =
X^{n-1} + 2Δt·RHS`). A term evaluated at time *n* and applied through the
centered step feeds the 2Δt computational mode directly; convection's large,
vertically sign-alternating heating excites it and the run NaNs in tens of
steps. Condensation survives only because its heating is small and smooth.

The fix is a proper coupling change, not a damping band-aid: moist-physics
temperature changes are now accumulated on the grid as a forward **increment**
in `wrk%dt_phys` [K], transformed to spectral, and added to the `n+1` state
*after* the dynamics step (`aeros_timestep_step` step 6) — a forward Euler
time-split, decoupled from the centered core. This is the same forward-in-time
treatment the gridpoint humidity already gets, so T and q are now symmetric.
Convection uses this path; condensation keeps the centered `dtdt` path (its
heating is benign and the path is validated). The cost is one extra spectral
analysis per level per step, only when convection is on.

### 12.2 Simplified Betts–Miller (M2.3f) — the scheme

`conv_scheme = "sbm"` (Frierson 2007) is now the default; Manabe stays as the
parameter-free reference (`"manabe"`). Per column:

- **Reference moist adiabat** anchored on the *actual* boundary-layer moist
  static energy `h_b = cp T_s + Φ_s + L q_s` (the θe of the surface air) — not
  the saturated value, which gives a reference so warm that `rh_ref·q_sat(T_ref)`
  exceeds the environmental humidity everywhere and no column ever rains.
  `T_ref(k)` solves `cp T_ref + Φ + L q_sat(T_ref) = h_b` by Newton — one solve
  per level, no iteration to neutrality (so §11 problem 1 simply does not arise).
- **Convecting band** = the layer where the parcel is buoyant in the MSE sense
  `h_b > h*_env`. Using MSE buoyancy rather than a temperature comparison folds
  in the LCL for free: an unsaturated boundary layer has `h*_env > h_b` at the
  surface, so the sub-cloud layer drops out and the band is the cloud layer
  (level of free convection to level of neutral buoyancy) — no dry-adiabat/LCL
  bookkeeping.
- **Reference humidity** `q_ref = rh_ref·q_sat(T_ref)`, `rh_ref = 0.7`.
- **Energy closure** (Frierson): with `Pq = ∫(q−q_ref)dp` and
  `Pt = ∫cp(T_ref−T)dp` over the band — the *deep* branch (`Pq > 0`) shifts
  `T_ref` by the constant that makes `∫cp(T_ref−T)dp = L·Pq`, so heating balances
  the latent heat of the precipitation, then relaxes T and q. The *shallow*
  branch (`Pq ≤ 0`, too dry to rain) does a zero-net-heating **temperature-only**
  relaxation, leaving humidity untouched — non-precipitating and `q ≥ 0` by
  construction. (This is a deliberate simplification of Betts–Miller's
  moisture-redistributing shallow branch: an additive `q_ref` shift can drive a
  very dry column's humidity negative, and there is nothing to tune the shallow
  reference against until there is radiation and observations. Revisit then.)
- **Implicit relaxation** over one physics step: `a = (Δt/τ)/(1+Δt/τ)`,
  `X ← X + a(X_ref − X)`, unconditionally stable in τ (`τ = 7200 s`).

Both branches conserve `∫(cp T + L q)dp` to machine precision (heating = latent
heat of the water removed) and give precip = column drying exactly. Namelist:
`conv_scheme`, `conv_tau`, `conv_rhref`.

`tests/test_convection` (rewritten) pins all of it: SBM deep (|ΔMSE|/MSE
1.3×10⁻¹⁶, precip = drying 6×10⁻¹⁶, no supersaturation, free troposphere warms),
SBM shallow (non-precipitating, moisture untouched, heat redistributed, q ≥ 0),
60-step convergence (energy conserved every step, precip decays as the column
neutralizes, no NaN), stable-column-untouched, and the Manabe moist and dry
branches retained under their explicit scheme.

### 12.3 Coupled result (M2.3g)

`tests/test_moist_run` now runs dynamics + transport + **convection** +
condensation together, 200 steps at T21L12, from a warm, humid, conditionally-
unstable initial state so the columns actually deep-convect:

| quantity | value |
|---|---|
| NaN over 200 steps | none |
| max \|u\| at the end | 9.2 m s⁻¹ (bounded) |
| convective precipitation | 1.3×10¹⁴ kg (positive — convection is active) |
| q < 0 ever | no |
| water-budget closure \|resid\|/W | 1.2×10⁻³ |

The stability is the whole point: the same convective heating through the
centered leapfrog NaNs in tens of steps; forward-split, it is stable and the
water budget still closes to the FV-vs-spectral gap (§8), unchanged in character.
Convective precip is smaller than the condensation total — the finite-τ
relaxation warms and dries gradually and large-scale condensation mops up the
rest — but it is unambiguously firing.

## 13. Radiation and a surface budget — operators done, RCE not (M2.4)

The big M2 piece: get away from Held–Suarez's artificial thermal forcing and
onto a real energy budget, so the moist line can eventually be validated against
climate rather than conservation. The **operators** landed and are tested; the
**coupled radiative–convective equilibrium** did not close, and that is the open
item carried to the next session.

### 13.1 The shortwave/longwave scheme (M2.4a–b)

Ported from CLIMBER-X/SESAM's broadband transmissivity–emissivity scheme
(`src/atm/{lwr,swr}.f90`; PIK report 81), design.md §5's pragmatic option. What
is taken is the *band physics* — the σT⁴ source, the H₂O/CO₂/O₃ absorber paths,
the empirical `d_vap`/`d_co2`/`d_o3` (LW) and two-exponential near-IR (SW)
transmission fits with the 1.66 diffusivity. What is discarded is SESAM's
statistical-dynamical *driver*, which reconstructs an analytic column from
lapse-rate descriptors and only ever wants boundary fluxes; aeros has a resolved
T/q/p column and needs a heating **profile**, so the kernel is driven with the
actual layer fields. On a resolved grid this is, structurally, design.md §5's
CCM3-style option 2 with SESAM's validated coefficients; ecCKD (option 1) slots
behind the same `scheme` selector later.

Two deliberate deviations: fluxes are staggered to the half levels so the layer
heating is a genuine flux divergence (SESAM co-locates them, wanting only
boundaries), and the water-vapour path is the exact hydrostatic `q dp/g` rather
than the exp-profile integral.

Offline single-column checks (`test_radiation`, a midlatitude column):

| quantity | value | sanity |
|---|---|---|
| OLR | 256 W/m² | clear-sky midlat, < σTs⁴ = 390 (greenhouse) |
| surface down-LW | 287 W/m² | physical (dry column, 12 kg/m²) |
| peak tropospheric cooling | −4.3 K/day | clear-sky range |
| flux-divergence identity | exact | conservative |
| CO₂ 280→560 forcing | +2.8 W/m² | canonical clear-sky TOA ~2.6–3 |
| planetary albedo (dark ocean) | 0.117 | physical |
| column SW absorption | 72 W/m² | ~ observed clear-sky |

The CO₂ forcing landing in the right pocket is the payoff over grey radiation —
the thing a one-band scheme cannot do.

**Insolation** is a present-day, circular-orbit, obliquity-only daily-mean
routine (declination + sunset hour angle → daily-mean TOA flux and
daylight-weighted airmass cosine zenith). Daily-mean, no diurnal cycle: matches
the multi-hour radiation cadence, and there is no surface heat store to lag yet.
Validated by the two textbook identities — equator-equinox = S₀/π = 433,
annual-global = S₀/4 = 340, polar night zero. A stopgap for the Laskar-2004
`insol` forcing (design.md §8).

### 13.2 Grid seam and surface budget (M2.4c–d)

Radiation is applied at the grid seam like the other physics: the full transfer
is recomputed every `interval` seconds (3 h default) and cached, the heating
held fixed between (design.md §5). Skin temperature is the surface module's
prescribed SST; ozone is the analytic profile of §13.3.

The **surface budget** (`aeros_surface`) is what closes the atmosphere's energy
once Held–Suarez is gone: a prescribed SST (APE control profile, 27 °C equator
to freezing poleward of 60°) drives bulk aerodynamic sensible and latent fluxes
into the lowest layer. Prescribed SST, not a slab — an infinite reservoir the
atmosphere equilibrates to, no surface energy balance to solve (design.md §6.1
for the slab later). No boundary-layer diffusion: convection carries the surface
fluxes up the column, the standard idealized-RCE closure. Unit-tested for flux
signs, exact flux/tendency bookkeeping, no-flux-at-equilibrium, and the q ≥ 0
guard under deposition.

### 13.3 Why the coupled run would not stay up (M2.4e–f)

`test_rce` — the full stack with Held–Suarez off — is stable and every part
active over a few hundred steps (OLR ~210 W/m², albedo ~0.11, surface fluxes and
precip firing, winds bounded). **But it is not an equilibrium test, because a
long integration does not stay stable.** The confirming long run (a scratch
driver, ~167 days) exposed a sequence of instabilities, each fix pushing the
failure out:

| configuration | blows up at | mechanism |
|---|---|---|
| radiation on centered path | ~day 9 | model top over-cools; sharp heating excites the computational mode |
| + forward-split radiation + ozone | ~day 20 | top still over-cools |
| + forward-split surface sensible heat | ~day 55 | top cools to 198 K, thermal wind runs away |
| + model-top sponge (C1+C2) | ~day 37 | top controlled; **lowest layer** develops a grid-scale hot spot (→375 K) instead |

Three fixes were correct and necessary and are kept:

1. **Forward-split heating.** Radiation and the surface sensible flux are large
   and vertically sharp (a cooling lid, a single-layer flux); on the centered
   leapfrog they excite the computational mode, exactly as convection does.
   Moved onto the forward-split `wrk%dt_phys` path — the handoff anticipated
   this ("dtdt is simplest unless it turns out large and stiff").
2. **Prescribed ozone**, with its shortwave heating deposited in the
   stratosphere (by ozone amount), so the model top gets a balancing heating.
   Too weak alone (~1 K/day vs the LW cooling), but part of the picture.
3. **Model-top sponge** (C1 Rayleigh drag + C2 Newtonian relaxation toward
   216 K, ramped over the top layers) — the cap Held–Suarez's `T_eq` floor used
   to provide. It controls the top over-cooling.

**The open item** is the surface-layer hot spot: with the top controlled, a
grid-scale hot spot grows in the **lowest** layer (level 12) and the run NaNs
after ~1–2 model months. It is a surface-flux / moist-physics interaction,
distinct from the model-top problem the sponge solves, and it needs the next
session (m2_handoff.md). Reaching a validated RCE — the prerequisite for the
ERA5 climate validation — waits on that.

## 14. ERA5 clear-sky radiation validation (M2.5a)

The clear-sky operators finally meet observations. Because the column kernels
(`aeros_lw_clearsky_column`, `aeros_sw_clearsky_column`) are grid-agnostic in
the vertical — they take an arbitrary `nlev` column of T, q, o3, `dp` and
interface height — they can be driven on **ERA5's own native pressure levels**
rather than on the model sigma grid. That is the deliberate choice here: no
regrid, so a transfer error is never confounded with an interpolation error,
and the radiative transfer is the only thing under test. This does not depend
on a stable RCE (§13.3) and so ran first.

`drivers/validate_era5.f90` reads the ERA5 1991–2020 monthly climatology
(`cdo ymonmean`, 2.5°, 144×73, 37 levels), averages the 12 months to an annual
mean, and for each 2.5° cell builds a column top→surface from the valid
(above-ground, non-fill) pressure levels — interfaces at level midpoints, base
at the ERA5 surface pressure, heights from the ERA5 geopotential — and runs
both operators. Shortwave is driven by ERA5 `tisr` for the TOA insolation (so
the stopgap insolation of §13.1 is not a confounder) and the annual-mean
insolation-weighted airmass; surface albedo is ERA5 `fal`; CO₂ is 380 ppm (the
period mean). ERA5 flux diagnostics are daily accumulations [J m⁻²], divided by
86400, and net fluxes are downward-positive, so `OLR = −ttrc`, etc. The output
is lat×lon maps of model, ERA5 and bias for six clear-sky fluxes
(`docs/figures/era5_rad_validation.png`, plotted by
`scripts/plot_era5_validation.jl`).

Area-weighted global-mean clear-sky biases (model − ERA5), annual mean:

| flux | model | ERA5 | bias |
|---|---|---|---|
| OLR (TOA up) | 248.2 | 264.1 | **−15.9** |
| surface down LW | 321.3 | 314.9 | +6.4 |
| surface net LW | −76.3 | −82.4 | +6.1 |
| TOA net SW | 283.7 | 288.8 | −5.1 |
| surface down SW | 258.2 | 242.3 | +15.9 |
| surface net SW | 222.6 | 212.2 | +10.3 |

The **spatial patterns match ERA5 closely** for all six fields; the biases are
modest and physically interpretable, not random:

- **OLR is ~16 W/m² too low, concentrated in the warm, moist tropics** — the
  broadband water-vapour band is slightly too opaque where the vapour path is
  largest. This is the one real tuning target the maps expose.
- **Downward/net surface LW is ~6 W/m² too high over deserts and high terrain**
  (Sahara, Arabia, Tibet, Antarctica) — too much down-LW in dry columns, the
  same band bias seen from below.
- **Surface SW (down and net) is uniformly too high (+16 / +10 W/m²)** — the
  clean-sky scheme carries no aerosol, so more shortwave reaches the surface
  than in ERA5's clear sky. The sign is exactly as expected; closing it needs an
  aerosol path, not a fit change.
- **TOA net SW bias is small (−5 W/m²)** with terrain/albedo-edge speckle.

The high-terrain speckle in the LW bias maps is the column-truncation and
crude interface-height construction over elevated surfaces (low surface
pressure, few levels), not a transfer error — cosmetic for this clear-sky check.

This is Tier 1 of the ERA5 spec (m2_handoff.md). Tier 2 (cloud fraction/water,
RH, u, v for the moist stack) is absent from the delivered data, so only
clear-sky is validated here; the cloudy-sky branch and the moist-line tuning
still wait on those fields and on a stable RCE.

## 15. The RCE lowest-layer instability, diagnosed and mitigated (M2.5b–c)

Task #7 from the M2 handoff -- the coupled-RCE blow-up after ~1 model month --
was diagnosed with an instrumented long-run driver (`drivers/rce_long.f90`,
per-level min/max T with location, a zonal-mean-vs-departure split, and the
per-level heating decomposition) and two principled fixes landed. It is not yet
a clean equilibrium, but the immediate blocker is gone.

### 15.1 The diagnosis: it is axisymmetric, and two mechanisms compound

The handoff read the blow-up as a *grid-scale* hot spot. It is the opposite: the
departure from the zonal mean is **exactly zero** (`max |T - zonalmean| = 0.00 K`)
until the very end. The run starts m=0 and, with zonally symmetric forcing,
never populates m>0 -- it is trapped in the axisymmetric manifold, and the "hot
spot" is a zonal-mean lowest-layer temperature at a subtropical latitude running
away. The max always landing at longitude index 1 was a tie-break over a zonally
uniform field, not a grid point.

Two compounding causes, isolated by toggling terms in the long run:

1. **Lowest-layer pile-up (physical).** The prescribed-SST surface injects heat
   and moisture into level 12 only; SBM convection mixes the cloud layer and
   excludes the sub-cloud layer, and there was no boundary-layer diffusion. The
   lowest layer built a 15-20 K warm/moist inversion. Diagnostics that ruled out
   the handoff's cheaper hypotheses: stronger ∇⁶ made it blow up *sooner* (not a
   horizontal grid mode), and halving `dt` only delayed it.
2. **Condensation on the centered leapfrog (numerical).** Once the pile-up is
   controlled, the latent release reaches ~100 K/day in the lowest layer at the
   hot latitude; on the centered path it excites the 2Δt computational mode (the
   zonal-mean there oscillates 379→297→357→382 K). `vdiff + condensation off`
   runs 200 days; `vdiff + convection off` does not.

### 15.2 The fixes

- **M2.5b -- boundary-layer vertical diffusion** (`aeros_vdiff`): down-gradient
  diffusion of T, q, u, v with a fixed eddy diffusivity tapering over a
  boundary-layer depth, mixing the surface source up the column. It is applied
  **implicitly on the n+1 state**, next to the horizontal diffusion and the
  sponge -- *not* forward-split. Forward-splitting a diffusion term onto the
  leapfrog is unconditionally unstable (the increment `dt·L·Xⁿ` on
  `Xⁿ⁺¹ = Xⁿ⁻¹ + …` gives a root pair with product −1, and it is scale-selective,
  so it blows up at the grid scale however small K is). The first, forward-split,
  implementation confirmed this by blowing up in a day; the implicit version is
  the one kept. Off by default; all 17 acceptance tests bit-unchanged.
- **M2.5c -- condensation forward-split**: move the latent heating from the
  centered `wrk%dtdt` to the forward-split `wrk%dt_phys`, the path convection,
  radiation and the surface fluxes use. Same fix M2.3e applied to convection,
  now warranted for condensation because its coupled-RCE heating is no longer
  "small and smooth." Drying and heating now share one discretization, so their
  budgets track exactly rather than to the filter's accuracy.

### 15.3 Where it stands

Together the two fixes carry the run from the day-34 blow-up to **200 model days
NaN-free** (vdiff `k0=25`, full stack). What remains is not a blow-up but a slow
secular warming of the subtropical lowest layer that has not equilibrated by day
200. Seeding a zonal asymmetry confirms baroclinic eddies then grow (T departs
from the zonal mean, winds ~20 m/s) and the early evolution is healthier, but the
run still eventually runs away -- so the residual is a **tuning/physics** problem
(surface exchange, convective `τ`, the K profile and depth, resolution at
T21L12), not the axisymmetric constraint alone and not a remaining coupling bug.
Reaching a bounded, validated RCE -- and with it the ERA5 moist-line validation
(§14 did the clear-sky line, independent of this) -- is the next tuning step.

## 16. Slab ocean and the defaults-file paradigm (M2.5d)

### 16.1 A responsive lower boundary — aeros_ocean

The prescribed SST is an infinite reservoir: it pumps heat into a column without
limit, which §15.3 showed a bounded RCE cannot tolerate. `aeros_ocean` makes the
sea surface a first-class, swappable component. It OWNS the SST (moved out of
aeros_surface, whose header had always anticipated this) and offers two modes
behind one interface:

- **prescribed** — a fixed SST (the APE control profile), the M2.4 behaviour and
  the default, so every existing run and test is bit-unchanged.
- **slab** — SST is prognostic, a well-mixed layer of depth `depth` (heat
  capacity `C = ρ_w cp_w depth`) integrated from the net surface energy balance
  `C dSST/dt = SW_net_sfc + LW_down_sfc − σSST⁴ − SH − LH` (forward Euler; freeze
  floor at T0, no sea ice yet). The slab makes the surface flux self-limiting.

This is deliberately the plug point for a real ocean (design.md §6.1): a
dynamical ocean replaces the module under the same contract — the atmosphere
hands it the net surface flux, it returns the SST. aeros_surface now takes the
SST as an argument; radiation reads it as the skin temperature and stores
`sw_net_sur` for the balance; aeros_timestep steps the ocean after surface and
radiation. `test_ocean` pins the APE profile, prescribed-inertness, and the slab
sign / equilibrium / freeze-floor.

**Result — necessary but not sufficient.** The slab is correct and, with a
cloud-proxy albedo to balance TOA, does what it should locally. But it does NOT
by itself bound the RCE: the subtropical **axisymmetric column instability of
§15.3 survives** the slab, a balanced albedo, and an eddy seed (slab+albedo→NaN
day 176; slab+albedo+eddies→day 88; the zonal-mean lowest layer at a subtropical
latitude still runs to 370+ K with runaway surface fluxes). That instability is
the real remaining blocker and is a resolution / convection-scheme / axisymmetric
dynamics problem, not a coupling one — see the handoff.

### 16.2 Defaults-file paradigm

Parameter files were required to be exhaustive (`nml.f90`'s `ERROR_NO_PARAM` is
true), which had quietly broken `aeros_run` and `held_suarez` the moment the M2
physics groups appeared. Now `input/aeros_defaults.nml` holds the complete
schema and case files carry only their overrides. Implemented by threading an
optional `defaults_file` through the load chain (`aeros_init` → `par_load`,
`vgrid_load`, `timestep_init` → each physics `*_load`), each `nml_read` passing
it on; `nml.f90`'s own `defaults_file` overlay (value from defaults, overridden
by the case file when present) does the rest. It is carried in the signature, not
in module state — `nml` stays stateless and shared across models. Absent
`defaults_file`, behaviour is exactly the legacy path, so the `_init`-only test
and driver paths are untouched.

## 17. ERA5 cloud radiative effect — the Tier-2 target (M2.5e)

The Tier-2 ERA5 fields the earlier sessions lacked (cloud fraction/water, RH,
u, v) are now delivered as 1991–2020 monthly climatologies on the same
144×73×37 grid as the Tier-1 clear-sky data (`~/data/era5/monthly-pressure-
levels`, `monthly-single-levels`). This unblocks the cloudy-sky radiation work.
The first step is to pin the **observational target** the model's cloudy-sky
branch (§13, still clear-sky only) must reproduce: the cloud radiative effect.

ERA5 ships paired all-sky and clear-sky flux diagnostics, so the CRE follows
directly with no model in the loop — CRE ≡ (all-sky − clear-sky) net downward
flux, positive = clouds warm the system, at both the TOA (`ttr`/`ttrc` thermal,
`tsr`/`tsrc` solar) and the surface (`str`/`strc`, `ssr`/`ssrc`). Fluxes are
ERA5 daily accumulations [J m⁻²], /86400 → W m⁻², downward-positive (the §14
convention). `scripts/plot_era5_cre.jl` reads the `_clim.nc` files, annual-means
the 12 months, and writes `docs/figures/era5_cre.png` (six CRE maps) and
`docs/figures/era5_moisture_clim.png` (column water vapour, total cloud cover,
and zonal-mean RH and cloud-fraction cross-sections). It is a pure-ERA5
diagnostic — a standalone Julia script, no Fortran driver.

Area-weighted global-mean CRE (annual mean, W m⁻²):

| | LW | SW | net |
|---|---|---|---|
| TOA | +21.8 | −46.1 | −24.2 |
| surface | +24.6 | −48.6 | −24.0 |

These sit right on the literature (CERES ≈ +26 / −47 / −21 at TOA; ERA5's own
clear-sky reference makes the net a touch stronger-negative, as expected). The
patterns are textbook: TOA LW warming peaks over the ITCZ and warm-pool deep
convection; SW cooling is strongest over the marine stratocumulus decks and the
storm tracks; net CRE cools the tropics and midlatitudes and is near-neutral
over the subtropical deserts. TCWV (24.3 kg m⁻²) and total cloud cover (0.63)
likewise match observations, and the zonal-mean RH shows the moist boundary
layer / dry subtropical mid-troposphere / moist ITCZ tower structure. The only
artifact is a cosmetic sub-surface step at the bottom-left of the cross-sections
over the high Antarctic plateau (1000/925 hPa levels below ground) — the same
column-truncation cosmetic noted for the §14 high-terrain speckle.

This is Tier 2's target. The cloudy-sky LW/SW branch that must reproduce it is
built and validated in §18.

## 18. Cloudy-sky radiation — the operators reproduce the CRE target (M2.5f)

The all-sky longwave and shortwave operators are ported the same way the
clear-sky ones were: the *band physics* of SESAM's cloudy branch on aeros's
resolved column, not SESAM's analytic slab. Cloud optics are built per layer
from the condensate paths, so — like the clear-sky kernels — the routines are
intensive and grid-agnostic (they run on ERA5's 37 levels here, or a refined
radiation grid remapped to the transport grid later, unchanged). Both are opt-in
siblings of the clear-sky kernels; `cf=0` recovers clear-sky bit-for-bit, so
every existing run and the 18 acceptance tests are untouched.

**Longwave (`aeros_lw_cloudy_column`).** A grey cloud is pure absorption, so it
folds into the clear-sky band transmission: each layer carries a transmission
`exp(-(k_liq LWP + k_ice IWP))` from its in-cloud water paths (standard
Stephens-type mass absorption coefficients), multiplied into the gas
transmission between interfaces (grey absorbers multiply; the gas band fit is
still evaluated once on the accumulated path). Overlap is SESAM's `lwr_total`
structure: a clear and an overcast column blended by the column cloud fraction
at **maximum overlap**, `F = (1−CF) F_clear + CF F_overcast`, `CF = max_k cf_k`,
the overcast column carrying the in-cloud condensate `grid-mean/CF` so the blend
conserves grid-mean water. The first cut used per-layer *random* overlap and
over-trapped (column cloud → `1−∏(1−cf_k)` → near-overcast; +15 W/m² TOA LW CRE
bias); maximum overlap, consistent with the shortwave, fixed it.

**Shortwave (`aeros_sw_cloudy_column`).** Shortwave is scattering, not grey
absorption, so the cloud enters as an albedo, not a per-layer transmission: the
tuned clear-sky column is run once, an overcast column is built by placing a
cloud reflector above it, and the two are blended by `CF`. The per-layer cloud
optical depth comes from the in-cloud water paths by geometric optics
`τ = 1.5 WP/(ρ r_e)`; a conservative-scattering two-stream on the column depth
gives the reflectance `R = γτ/(1+γτ)`, `γ=(1−g)/2μ`, with a small near-IR
absorptance; the cloud reflector is added above the clear column
(`alb_ov = R + T² alb_clr/(1−R alb_clr)`) so the overcast surface fluxes are the
clear ones times `T/(1−R alb_clr)`. Absorption is the exact TOA-minus-surface
residual, distributed by the clear heating plus the cloud optical depth, so the
column energy identity holds to machine precision (as it does for the longwave
flux divergence).

**Validation.** `drivers/validate_era5.f90` now drives the cloudy operators on
ERA5's own `cc`/`clwc`/`ciwc` columns and compares the resulting cloud radiative
effect (all-sky − clear-sky) against ERA5's, cell by cell
(`docs/figures/era5_cre_validation.png`, three rows model/ERA5/bias for TOA LW,
SW and net CRE). Area-weighted global means [W m⁻²]:

| CRE (TOA) | model | ERA5 | bias |
|---|---|---|---|
| LW  | +19.7 | +21.8 | −2.1 |
| SW  | −44.5 | −46.1 | +1.6 |
| net | −24.7 | −24.2 | −0.5 |

(These are with the §20 water-vapour opacity correction; before it the LW CRE
was +18.2, bias −3.6, and the net −26.3.) All three land within a couple of W/m²
and the **spatial patterns match** — LW warming peaking over the ITCZ and
warm-pool deep convection, SW cooling over the marine stratocumulus decks and
storm tracks, net cooling concentrated in the tropics. The small residual −2.1
LW CRE is consistent with maximum overlap running slightly weak where cloud is
in separated decks (maximum-random overlap is the identified refinement); the SW
and net are essentially on target. This validates the *radiative transfer* of
the cloudy branch against observations, offline — decoupled from coupled-model
cloud generation, which §19 now supplies.

## 19. Diagnostic clouds wired into the coupled model (M2.5g)

The all-sky operators need a cloud column; the model has no prognostic cloud
water or fraction (condensation removes condensate on the spot, §9), so
`aeros_cloud` diagnoses one from the resolved T/q/p: an RH-based cloud fraction
(Sundqvist, `cf = 1 − √(1 − b)`, `b = (RH − rhc)/(1 − rhc)`, critical RH a
constant-by-default but two-endpoint vertical profile) and a temperature-
dependent in-cloud condensate (`q_ic = f · q_sat`, split liquid/ice by a
temperature ramp), returned as the grid-mean `clwc`/`ciwc` the kernels expect.
`aeros_radiation_apply` calls it per column and switches to the cloudy kernels
when the new `radiation%clouds` flag is on (default off → clear-sky, every run
and the 18 tests bit-unchanged). `test_radiation` checks the diagnosis gives
physical `cf`/water and that fed through the all-sky kernel it reduces OLR.

**The coupled TOA-balance check exposes that the first-cut cloud parameters are
far too aggressive for the model's moist RCE state.** Running the bounded
rotating vehicle (`rce_long`, `c_d>0`, prescribed SST, physical ocean albedo
0.06) clouds-off vs clouds-on: the global TOA net flux swings from **+88 W/m²**
(dark ocean, no cloud reflection) to **−100 W/m²** — OLR cut 208→140, planetary
albedo ~0.85, a near-total overcast. The sign is right (clouds cool the planet
and the balance moves the right way) but the magnitude is ~8× the realistic net
CRE. The cause is not the kernels (the ERA5-driven validation, §18, is on target
with *realistic* input cloud) but the diagnosis: with `RH_crit = 0.70` and
`q_ic = 0.04 q_sat`, the model's near-saturated deep convective columns produce
near-overcast, optically thick cloud everywhere. Tuning the diagnosis against
ERA5 is §21; the wiring, opt-in and validated end-to-end, is what lands here.

## 21. Tuning the diagnostic cloud scheme against ERA5 (M2.5i)

Tuning was done against ERA5's *own* columns, not the coupled model: a
diagnosed-cloud path in `validate_era5` runs `aeros_cloud` on ERA5 T/q/p, then
the cloudy kernels, and compares the diagnosed CRE and cloud cover to ERA5's
(`docs/figures/era5_cre_diagnosed.png`). This makes the scheme correct *given
correct inputs* and keeps a model humidity bias from being hidden in the cloud
tuning.

**The first-cut `q_ic = f·q_sat` was structurally wrong** and the ERA5 harness
showed why: on realistic columns it gave the CRE ~2× too weak in the LW and ~2×
too strong in the SW, robust to `RH_crit`. Tying condensate to `q_sat` starves
cold high cloud (cirrus, `q_sat` tiny → negligible ice water → no LW effect) and
over-thickens warm low cloud (large `q_sat` → bright → too much SW). **Fix
(committed): a specified in-cloud water content** [kg m⁻³] converted to a mixing
ratio by the air density, ice much thinner than liquid, so cirrus gets a real
ice-water path (LW) and low cloud is not over-bright (SW). Liquid and ice are
then independent knobs, and because LW is ice-dominated (high cloud) and SW
liquid-dominated (low cloud), they decouple.

Tuned (`RH_crit` 0.52→0.35 in σ, LWC 0.017, IWC 0.018 g m⁻³), diagnosed on ERA5:

| CRE (TOA) | diagnosed | ERA5 | bias |
|---|---|---|---|
| LW  | +19.6 | +21.8 | −2.2 |
| SW  | −43.2 | −46.1 | +2.9 |
| net | −23.6 | −24.2 | **+0.7** |

The net CRE is within 0.7 W/m² and LW/SW each within 3; the patterns broadly
match (LW warming over deep convection, SW cooling over the stratocumulus decks
and storm tracks). Two residuals: total cloud cover is low (0.47 vs 0.63 — the
scheme trades cover for optical depth but gets the CRE right), and the spatial
distribution has compensating errors (a little too much subtropical cloud
effect, too little in the deep tropics — the RH diagnosis makes low/subtropical
cloud more readily than deep-tropical high cloud).

**The coupled check now cleanly separates a model humidity bias.** Re-running
the `rce_long` TOA-balance check (clouds on) with the tuned scheme: the TOA net
flux goes from clouds-off **+88** to clouds-on **−36 W/m²** — much improved from
the pre-tuning −100 (overcast), but still too negative. Since the scheme is
validated correct on ERA5's columns, this residual overcast is not the cloud
scheme: it is the model's RCE being **too moist** (near-saturated over a deep
layer → high `cf`), so it makes more cloud than ERA5 would. That is exactly the
finding the ERA5-column tuning was chosen to expose (§19's decision): the cloud
scheme is now sound, and the coupled overcast is a humidity/RCE problem to chase
separately.

## 20. The clear-sky OLR opacity bias — a vapour-opacity correction (M2.5h)

§14 found the clear-sky OLR ~16 W/m² too low, concentrated in the warm moist
tropics: the SESAM broadband water-vapour band is too opaque where the vapour
path is largest. The SESAM coefficients are one internally-consistent fit not
meant to be re-fit piecemeal (and the structural fix is the eventual ecCKD
scheme), so rather than repartition them this adds a **single documented
correction** `LW_VAP_OPAC` multiplying the whole vapour optical depth in
`d_vap = 1/(1 + vap_opac·τ_vap)`. `vap_opac = 1` is SESAM verbatim (bit-for-bit);
`< 1` reduces the opacity. Because the term is negligible at small path and
dominant at large path, one scale corrects the tropics and barely touches the
dry columns — it targets exactly where the bias lives.

A single knob cannot null the TOA and surface longwave bias at once — they
trade off, so it is set to minimize the combined clear-sky LW error:

| `vap_opac` | OLR bias | sfc-down-LW bias |
|---|---|---|
| 1.00 (SESAM) | −15.9 | +6.4 |
| 0.80 (kept)  | −7.1  | −5.0 |
| 0.65         | +0.9  | −16.7 |

At **0.80** the clear-sky OLR bias more than halves (−15.9 → −7.1) while the
surface-down-LW magnitude stays comparable (+6.4 → −5.0, the one cost — its sign
flips). The shortwave is untouched (a separate band); the CO₂-doubling forcing
stays canonical (2.8 → 3.0 W/m², the vapour band masking the CO₂ bands slightly
less); the greenhouse and all `test_radiation` checks pass. The correction also
helps the cloudy branch: with the clear-sky OLR less biased there is more
headroom for cloud, so the §18 TOA cloud radiative effect improves — LW CRE
−3.6 → −2.1, net CRE −2.0 → −0.5 against ERA5. `vap_opac` is a named parameter
(one place, documented) to be retired when ecCKD lands.

## 22. The coupled-RCE humidity bias, localized against ERA5 (M2.5j)

§21 left the overcast coupled TOA (`rce_long` clouds-on −36 W/m²) attributed —
by inference — to a *model* moisture bias, since the cloud scheme is validated
correct on ERA5 columns. This section measures it directly. `rce_long` now dumps
its end-of-run zonal-mean RH, diagnosed cloud fraction, T, q and cover (arg 2 =
output NetCDF; not a namelist key, so no `rce_*.nml` needs editing), and
`scripts/rce_humidity_vs_era5.jl` lays it against the ERA5 1991–2020 RH
climatology (`r`) — figure `docs/figures/rce_humidity_vs_era5.png`. Two vehicles
share an identical surface boundary (prescribed SST) and differ *only* in the
dynamics, to separate a structural (missing-subsidence) bias from a column-physics
one: **rot** (rotating, realistic insolation) and **uni** (non-rotating, uniform
insolation — each column ≈ a local 1-D RCE). Both 100-day, clouds on.

**The bias is column-wide and enormous, not a subtropical detail.** The model is
near-saturated (RH ~85–100%) at essentially every latitude and height; ERA5 runs
17% aloft, ~45% in the mid-troposphere, ~78% in the boundary layer. Area-weighted
mean RH, model top → surface:

| p [hPa] | rot | uni | ERA5 | rot−ERA5 |
|---|---|---|---|---|
| 74  | 99.8 | 99.6 | 16.9 | **+82.9** |
| 198 | 93.1 | 98.1 | 40.8 | +52.3 |
| 313 | 89.2 | 88.3 | 47.7 | +41.5 |
| 519 | 87.5 | 98.2 | 43.3 | +44.2 |
| 693 | 96.5 | 76.2 | 48.6 | +47.9 |
| 891 | 93.9 | 83.6 | 72.4 | +21.5 |
| 981 | 97.3 | 87.0 | 77.1 | +20.2 |

Cloud cover follows: rot 1.00, uni 0.96, ERA5 0.63 — the model is globally
overcast, which is the −36 W/m² directly.

**The decisive result: it is column physics, not missing subsidence.** Two
signatures both point the same way. (1) The non-rotating uniform vehicle — which
has essentially no circulation, so no dynamical subsidence anywhere — is *equally*
saturated (uni ≈ rot). If the bias were the missing subtropical Hadley drying, it
would appear in rot only where the subtropics ought to be dry, and uni would be
free of it. It is in both. (2) In the rot cross-section the *only* sub-saturated
region is the tropical upper troposphere (200–600 hPa, ±30°) — exactly where
deep convection is actively firing. Everywhere convection is not actively
overturning, RH pins at ~100%. So the model has **no free-tropospheric drying
mechanism except active deep convection**: SBM relaxes q toward `rh_ref`·q_sat
only in buoyant columns, and large-scale condensation only *caps* supersaturation
at 100% and leaves it there — nothing ever brings a layer back below saturation.
The classic 1-D-RCE "everything saturates without subsidence drying" attractor,
here severe enough to overcast the planet. (Aloft, RH ~100% at 74 hPa vs ERA5's
17% is a second, smaller issue: no tropopause cold-trap / freeze-drying and q_sat
is tiny there, so a trace of trapped q saturates it; little condensate results,
so it matters less for the cover than the deep 500–900 hPa moist layer.)

This is the barrier to a physical coupled TOA balance, and it is now localized:
the fix must give the non-convecting free troposphere a way to dry (a subsidence
/ environmental-descent drying, a convective downdraft/detrainment that removes
column water, or a re-evaporation path that can pull q below saturation), not a
cloud-scheme knob.

## 23. Fixing the overcast: a sub-grid-saturation condensation ceiling (M2.5k)

Two side diagnostics narrowed the fix before it was made. **(1) Convection is
dormant at equilibrium.** The `rce_long` buoyancy diagnostic (`l_diag`) shows the
boundary-layer parcel MSE below the saturated environmental MSE at *every* level
(`hb − h*_env < 0`, surface −1.5 kJ/kg) once the column has equilibrated, so
SBM's trigger never fires (`cnv = 0` at all levels) — though at the RH-90% initial
state it *is* buoyant (levels 3–11, +7 to +29 kJ/kg) and fires during spin-up.
Convection runs, saturates the column, and thereby raises `h*_env` (the `L·q_sat`
term) until it shuts *itself* off; large-scale condensation then holds the column
at RH~100% because it removes only the excess above saturation. So a convective
downdraft (the first-choice fix) is moot — there is no active convection at
equilibrium to carry one. **(2) Boundary-layer diffusion is a minor, low-level
contributor.** Confining vdiff (`vdiff_sigma` 0.7→0.9) dries only ~800 hPa
(RH 80→57%) and leaves the mid/upper troposphere saturated and the planet
overcast (cover 0.96–1.00 throughout the sweep); it never reaches the levels that
set the cover.

**The fix is the drying sink the column lacks, and it already existed as a
dormant knob.** `aeros_condensation` relaxes q toward `rh_crit·q_sat` (module
default 1.0 = true saturation adjustment); a value below 1 is the standard
sub-grid-saturation stand-in — condensation removes water before the *grid mean*
saturates, representing the saturated fraction of a partly-cloudy box. `rce_long`
had hard-coded it to 1.0; it is now a namelist knob (`cond_rh_crit`). Lowering it
gives the free troposphere the sink it was missing: q can no longer sit above
`cond_rh_crit`, and — by keeping the environment sub-saturated — it also lets
convection stay active. On the uniform (single-column) vehicle the sweep is a
clean monotone control on RH, cover and energy:

| `cond_rh_crit` | free-trop RH | cover (uni) | TOA net (uni) |
|---|---|---|---|
| 1.0 (was) | ~98% | 0.96 | −29 |
| 0.9  | ~90% | 0.60 | +7.9 |
| 0.8  | ~80% | 0.44 | +25.7 |
| 0.7  | ~70% | 0.31 | +38.9 |

**Chosen value `cond_rh_crit = 0.93`**, verified on the rotating vehicle (the
−36 W/m² case of §22):

| vehicle | TOA net | cover | (ERA5) |
|---|---|---|---|
| rot, was (1.0) | −36 | ~1.0 | — |
| rot, 0.93 | **+8** | **0.65** | cover 0.63 |
| uni, 0.93 | **+0.7** | 0.66 | — |

The overcast is gone: cover 0.65 against ERA5's 0.63, and the TOA net is +8
(rot) / +0.7 (uni), from −36. Figure `docs/figures/rce_humidity_fixed.png` (same
layout as §22's `rce_humidity_vs_era5.png`, the "before").

**The honest limitation — what `cond_rh_crit` does and does not fix.** It is a
*ceiling*, not subsidence drying, and it cannot match cover and the RH profile at
once. At 0.93 the cover and TOA are right, but the free-troposphere RH is still
~85–90% against ERA5's ~45%; driving RH to ERA5 values (`cond_rh_crit ≈ 0.6`)
under-clouds the planet (cover 0.31) and swings the TOA to +30. This trade-off is
the fingerprint of the genuinely-missing process — large-scale subsidence, which
would dry RH while cloud cover is set independently. `cond_rh_crit = 0.93` buys a
physically-balanced coupled TOA and a realistic global cloud cover *now*; closing
the residual RH bias (and the cover's flat latitudinal structure, which the RCE's
weak circulation cannot shape into ERA5's dry-subtropic / moist-storm-track
pattern) is deferred to a subsidence-drying treatment (§22's option C). The knob
default stays 1.0 (bit-reproducible); 0.93 is set per-run in the RCE configs.

## 24. ERA5 moist-line validation of the coupled RCE (M2.5l)

The M2 validation target (`docs/design.md`): the coupled radiative-convective
equilibrium judged against the ERA5 1991–2020 climatology in the zonal mean —
temperature, the zonal-wind jet, and RH. This closes the loop the clear-sky line
(§14) and the CRE (§17–23) left open: the whole moist stack run to equilibrium
and compared to *climate*, not conservation. The `rce_long` dump now carries
zonal-mean `u`/`v` alongside `T`/`q`/`RH`, and `scripts/rce_validate_era5.jl`
regrids ERA5 `t`/`u`/`r` onto the model grid for the comparison. Driven on the
balanced rotating vehicle (`cond_rh_crit = 0.93`, §23); figure
`docs/figures/rce_validate_era5.png`.

**Temperature — the mid/lower troposphere is right, two biases aloft and at the
pole.** Area-mean |bias| < 3 K from 900 to 420 hPa; the moist-adiabatic structure
the convection scheme builds matches ERA5 through the bulk of the troposphere.
Two departures: (a) an **upper-tropospheric warm bias** (+8 K at 313 hPa, +12 K
at 198 hPa), strongest in midlatitudes — the L12 upper levels span a deep layer
and the detraining moist adiabat runs too warm there; (b) a **polar
lower-tropospheric warm bias** (off-scale at the SH pole), a prescribed-SST
artifact — this vehicle has no cold poles or sea ice, so the near-surface polar
air cannot cool to ERA5 values. Near-surface is otherwise slightly cold
(−1.4/−2.8 K at 940/981 hPa).

**Jet — right position, a quarter of the strength.** The model's zonal-mean jet
peaks at **±36°, 198 hPa** — ERA5's jet latitude and level almost exactly — but at
only **7.4/7.5 m/s against ERA5's 27/31**. So the thermal-wind response to the
equator–pole gradient puts the jet in the right place, but it reaches ~25% of the
observed speed: the axisymmetric RCE has no baroclinic-eddy momentum-flux
convergence to spin the jet up to full strength. The `u`-bias is a deep
easterly-side deficit centered on both jet cores. This is the direct motivation
for the eddy question (handoff step 2): whether a seeded T21 jet grows the eddies
that would close this gap, or T42 is needed.

**RH — as §23.** Realistic where deep convection dries (the tropical
upper-troposphere minimum is reproduced), too moist elsewhere: the `cond_rh_crit`
ceiling holds the free troposphere near 85–90% against ERA5's ~45%, the standing
subsidence-drying limitation.

**Verdict.** The coupled RCE reproduces the gross thermal structure and the jet
*position*; every deficit is consistent and physically attributable — weak jet
(no eddies), moist free troposphere (no subsidence drying), warm upper-troposphere
(coarse L12 + detrainment), warm poles (prescribed SST) — none random. This is a
credible first Tier-2 validation and a clean baseline for the next steps.

## 25. Do baroclinic eddies grow? T21 vs T42, and why the jet stays weak (M2.5m)

§24 left the jet at the right latitude/level but ~25% of ERA5 strength, the
signature of absent eddy momentum-flux convergence. This tests directly whether
seeded baroclinic eddies grow and spin the jet up, and whether the truncation
(T21) is the limiter — the standing handoff step-2 question. All runs are the
balanced rotating RCE (`cond_rh_crit = 0.93`, `c_d = 1.5e-3`), 200 days;
`eddy_diag` reports mid-tropospheric eddy KE and `[v'T']` each 10 days, and the
seed (`seed_asym`) perturbs zonal wavenumbers 1–6.

| config | eddy KE (saturated) | zonal-mean jet | jet position |
|---|---|---|---|
| T21 control (no seed) | ~0 (machine zero) | 6.9 / 7.0 | ±36°, 198 hPa |
| T21 seed, ∇⁶ τ=6 h | ~0.7 | 9.6 / 9.1 | ±25°, 313 hPa |
| T21 seed, ∇⁶ τ=24 h | ~1.0 | ~9–10 | (noisy) |
| **T42 control** | ~0 | 7.7 / 7.9 | ±35°, 198 hPa |
| **T42 seed** | ~0.9 | 10.4 / 10.3 | ±24° |
| ERA5 | ~10²–10³ | 27 / 31 | ±30–36° |

**Three findings, one conclusion.** (1) *The axisymmetric trap is real:* with no
seed the run stays exactly on the m=0 manifold — eddy KE at machine zero — because
zonally symmetric forcing never populates m>0. (2) *The seed does not take off:*
from an initial pulse (eddy KE ~3–4 m²/s²) it **decays** over ~40 days and
saturates at ~1 m²/s² — ~100× below the observed storm track — with `[v'T']` at
noise level. (3) *Neither cure works:* slackening the ∇⁶ hyperdiffusion 4×
(τ 6→24 h) barely changes the saturated level, and **T42 gives essentially the
same answer as T21** (eddy KE ~0.9 vs ~0.7; jet 10.4 vs 9.6). So the weak jet is
neither a numerical-diffusion artifact nor a truncation limit — doubling the
resolution does not unlock the eddies.

**What it is.** The limiter is the RCE state itself: its baroclinicity is weak and
surface-trapped (the equator–pole gradient and jet peak at the lowest level,
decaying upward — the same structure that drove the §15 low-level-jet blow-up),
so the baroclinic instability that would build a storm track is marginal at any
resolution here. Seeded eddies do transport *some* momentum — they lift the
zonal-mean jet from ~7 to ~10 m/s and shift it equatorward (±36°→±24°), a real but
small effect — but cannot close the gap to ERA5's 27–31 m/s. Closing it needs
*stronger baroclinic forcing* (a realistic meridional SST gradient, a seasonal
cycle, land–sea contrast), not more spectral resolution — a modeling step beyond
the aquaplanet RCE and deferred. For the M2 validation the takeaway is bounded:
the jet's *position* is right and its *weakness* is now explained and quantified,
not a mystery.

## 26. Why the jet is weak: a thermal-wind diagnosis, and the right yardstick (M2.5n)

The jet weakness (§24–25) has a precise mechanism. The zonal-mean jet is
thermal-wind balanced, so its shear integrates the meridional temperature
gradient over height. Comparing the equator-minus-midlatitude temperature
contrast ΔT by level, model vs ERA5:

| p [hPa] | model ΔT | ERA5 ΔT |
|---|---|---|
| 981 (sfc) | +16.5 | +15.6 |
| 891 | +15.6 | +15.2 |
| 693 | +10.4 | +14.0 |
| 519 | +7.8 | +14.8 |
| 420 | +2.3 | +15.8 |
| 313 | +1.1 | +14.8 |
| 198 | −0.1 | +1.6 |

**At the surface the model's gradient matches ERA5 almost exactly (+16.5 vs
+15.6 K) — the APE SST forcing works — but it collapses to ~0 by 300–420 hPa,
while ERA5 holds ~15 K up to 300 hPa.** ERA5's deep gradient builds a 34 m/s jet
aloft; the model's vanishing gradient leaves 7.5 m/s, surface-trapped. That is the
whole of the jet deficit.

**What sets the collapse — a two-latitude heating budget.** `rce_long`'s per-term
diagnostic (`l_diag`) now reports at the tropical hot latitude *and* at ~45°
(`term_table` takes a latitude index). At 45° at equilibrium, in the free
troposphere: **convection is off** (`cnv = 0` at every level; `hb − h*_env` is
−18 to −45 kJ/kg), **condensation heating is surface-confined** (`cnd ≈ 0` above
the boundary layer), and **radiation cools** (−0.3 to −1.2 K/day), balanced by
weak subsidence/diffusion. So moist physics is *not* actively warming the
midlatitude free troposphere — the earlier "midlats driven onto a moist adiabat"
guess is refuted for the equilibrium state. The warm midlatitude upper
troposphere is a spin-up legacy: convection warmed it early, and once it shut off
there is **no baroclinic-eddy ventilation** to pull it back to a cold,
baroclinic state. Collapsed gradient and missing eddies (§25) are one self-locking
problem: no eddies → warm midlat aloft → weak upper-level gradient → weak
surface-trapped thermal wind → no eddies; §25 showed it does not break at T42.

**The right yardstick.** Much of this mismatch against ERA5 is *expected*, not a
model defect: this is an aquaplanet with perpetual annual-mean forcing, and Earth's
extratropical storm tracks are organized and energized by exactly what the setup
omits — land–sea contrast and stationary waves, and a seasonal cycle. Judging the
model's *circulation* against ERA5 (full reality) therefore overstates the
"error"; the fair reference is the aquaplanet intercomparison (APE, Neale &
Hoskins), where the target jet is itself weaker and more zonally symmetric than
ERA5's. Two things remain genuinely worth noting even against that fairer bar:
the surface gradient and jet *position* are right (the dynamics respond
correctly), and the free-tropospheric gradient collapse is a real characteristic
of this steady moist RCE. The thermal/moisture columns (§22–24) — where the
comparison to reality *is* fair — validate well; the extratropical circulation is
where a steady aquaplanet is simply not expected to reproduce ERA5, and chasing
that with parameterization tuning would be fitting the wrong target. Closing it
belongs to the capabilities a paleoclimate model needs anyway: a seasonal cycle
(orbital insolation) and, eventually, land and a dynamical ocean.

## 27. Seasonal cycle and orbital forcing via the insol package (M3a)

The first of those capabilities. All the RCE work above ran on **steady
annual-mean insolation**; a paleoclimate model needs the seasonal cycle (ice-sheet
SMB is a summer-melt quantity) and orbital (Milankovitch) forcing as a genuine
axis (design.md §6). Rather than grow the earlier present-day, circular,
obliquity-only stopgap, aeros now uses the **fesmc/insol package** (Laskar et al.
2004 orbital elements, `~/models/insol`), which delivers both at once.

- **New module `aeros_insolation`** wraps insol: an `insol_class` built once with
  the grid latitudes, supplying annual-mean and per-day `sw_toa` / `coszen`
  (`coszm_kind="flux"`, the insolation-weighted airmass cosine the shortwave band
  wants). `aeros_radiation` calls it in init (annual mean) and, in `seasonal`
  mode, on the recompute cadence. The old `aeros_insolation_daily` is retained
  only as the `test_radiation` analytic reference.
- **Orbital axis:** a namelist knob `time_bp` (orbital year before present, 1950
  CE; 0 = present-day) threads straight into insol, so the *same* code gives any
  epoch's seasonal cycle — the paleoclimate forcing, for free.
- **Build:** insol is a prebuilt dependency reached through an `insol` symlink
  (gitignored, like `fesm-utils`); it carries no external deps and ships its own
  LA2004 tables, which are read through the symlink so nothing is duplicated into
  the aeros tree. Radiative calendar `DAY_YEAR = 365`.

**Validation** (`drivers/probe_insol.f90` → `make probe-insol`;
`scripts/plot_insol_seasonal.jl` → `docs/figures/insol_seasonal.png`):

| check | result | expected |
|---|---|---|
| annual + global mean | 340.28 W/m² | S0/4 = 340.25 |
| solstice (d172) NH pole | 523.8 W/m² | bright (24-h day) |
| solstice (d172) SH pole | 0.0 W/m² | polar night |
| LGM (21 ka BP) NH-pole solstice | 535.6 W/m² | > present (orbital) |

The annual-global mean lands on S0/4 to 0.03 W/m², the solstice reproduces the
polar day/night bow-tie, and the 21 ka run shifts NH summer-pole insolation up
~12 W/m² at a fixed annual-global mean — the Milankovitch redistribution, live.
All 18 acceptance tests still pass (the insolation swap preserves the physical
relations the radiation/RCE tests assert: equator ≈ S0/π, polar night = 0, global
mean ≈ S0/4). Unlike the extratropical circulation (§26), the seasonal cycle *is*
fair to validate against ERA5, which has seasons — the natural next use is a
seasonal coupled run and the SMB-relevant seasonal temperature cycle.

## 28. ecCKD correlated-k radiation — the structural replacement for SESAM (M3b)

The radiation endpoint design.md §5 named. The ported SESAM broadband scheme
needed a hand-tuned vapour-opacity fudge (`LW_VAP_OPAC`, §20) to half-cancel a
structural clear-sky OLR bias, and a single knob can't null OLR and surface-LW at
once. A correlated-k scheme fixes that at the source. Built opt-in behind the
existing `scheme` selector (`SCHEME_ECCKD`); **SESAM stays the default and is
bit-for-bit unchanged** (`scheme=1`), so every prior result and all previous tests
hold. Developed on a branch in parallel with §27 (a background agent, iterated to
these results); merged clean.

**Gas optics without a data table.** No ecCKD/CKDMIP look-up table exists offline
and the generator needs tens of GB of line-by-line input, so the k-tables are
built from a compact **Malkmus/Goody statistical band model**: each band's
k-distribution is analytically the inverse-Gaussian (Wald), so its g-point
`(k, weight)` pairs follow from Gauss–Legendre nodes through the closed-form
inverse-Gaussian CDF — no fitting, no run-time data. Planck fraction integrated
exactly per band per layer (the dominant T-dependence), a mild Malkmus `S(T)`, and
a cheap `(p/p0)^κ` pressure broadening (κ a-priori Lorentz, **not** ERA5-fit — the
vertical structure must stay predictive at non-present orbits/climates).

- **10 LW g-points** (4 bands: H2O rotation, CO2 15 µm, window + O3 9.6 µm,
  H2O vib-rotation) + **5 SW** (visible with a 3-g-point correlated-k O3, near-IR
  reusing SESAM's validated 2-exponential H2O) = **15 total**, inside the 8–32
  budget.
- **Tuning discipline:** ≤2 ERA5-anchored magnitude scales per band; only a per-
  band H2O scale and the CO2 band strength were used (LW), and **none** on the SW
  (a-priori O3 validated as-is). Band *shape* stayed a-priori throughout.

**Validation vs ERA5 1991–2020, and vs SESAM — better on every headline, fudge
retired (`LW_VAP_OPAC = 1`):**

| metric | ecCKD | SESAM (fudged) | ERA5 target |
|---|---|---|---|
| clear-sky OLR bias | **−5.1** | −7.1 | ~0 |
| 2×CO2 forcing | **3.6** | 3.0 | ≈3.7 |
| clear-sky TOA net SW bias | **−4.5** | −5.1 | small |
| TOA LW cloud effect | +21.7 (**−0.1**) | +19.7 (−2.1) | +21.8 |
| TOA net cloud effect | −23.3 (+0.9) | −24.7 (−0.5) | −24.3 |

The clouds fold SESAM's validated grey per-layer optics into the g-point
transmissions with **no new spectral points** (LW runs the transfer worker twice,
clear + overcast with cloud optical depth added into every g-point at max overlap;
SW places the two-stream cloud reflector on the ecCKD clear base). The LW CRE
gain (−0.1 vs SESAM −2.1) is free: the accurate clear-sky OLR leaves the right
cloud headroom (the §20 mechanism, now structural).

**Cost (the fast-model constraint):** per-column LW 1.78× SESAM, SW 0.84×
(*faster* — drops SESAM's unused paths), cloudy LW 1.77×; radiation runs on the
1–3 h cadence, so the runtime share stays small. Energy closure (heating = flux
divergence) is exact by construction, LW and SW.

**Coverage:** a new acceptance test `tests/test_ecckd.f90` (suite 18→19) drives the
`SCHEME_ECCKD` path and locks in energy closure, CO2 monotonicity, and the
reduce-to-SESAM SW property. `validate_era5` gained a scheme arg (`sesam|ecckd`)
and a CO2 arg. ecCKD is the recommended scheme going forward; SESAM is retained as
the fast, bit-reproducible default until the switch is made deliberately.

## 29. Coupled-RCE revalidation under ecCKD, and the radiation-cost fix (M3c)

ecCKD was made the production default (`scheme` default flipped SESAM→ecCKD, §28
follow-up). Because the radiation that sets the coupled equilibrium changed, every
earlier coupled-RCE result and tuning (§22–26, all under SESAM) needed a
revalidation pass. Vehicle throughout: the balanced rotating slab, `cond_rh_crit
= 0.93`, clouds on (`logs/rce_revalidate.nml`). All runs NaN-free.

**Revalidation — nothing forced a re-tune.**

1. **Coupled TOA + cover sweep.** `cond_rh_crit` ∈ {0.88, 0.90, 0.93, 0.96}, 100 d:
   OLR falls / cover rises monotonically with crit (169.0/0.555, 164.3/0.590,
   154.9/0.658, 143.8/0.735 W/m²·—). **0.93 stays the best joint point** — cover
   0.658 vs ERA5 0.63 (+0.03), TOA net **+5.6 W/m²** (was +8 under SESAM, §23).
   TOA only closes at 0.96 (+0.4) but that overcasts to 0.735 (+0.11 vs ERA5); the
   residual is the Tier-1 opacity/albedo term, not a cloud-tuning target.
2. **Moist-line vs ERA5** (`scripts/rce_validate_era5.jl`,
   `docs/figures/rce_validate_era5_ecckd.png`). Upper-trop warm bias **eased ~2 K**
   by ecCKD's more accurate LW aloft: 198 hPa +9.9 K (was +12 K, §24); mid-trop
   ±1 K. **Jet unchanged** — 8.6 vs ERA5 27 m/s at 198 hPa: a dynamics/resolution
   limit (§25–26), radiation-independent, as predicted.
3. **Humidity/cover** (`rce_humidity_vs_era5.jl`,
   `docs/figures/rce_humidity_vs_era5_ecckd.png`). Story unchanged from SESAM:
   near-saturated free troposphere (RH +40 to +42% vs ERA5's 40–55% mid-trop),
   ~ERA5 cover (0.66 vs 0.63) because the cover diagnostic saturates. This is a
   missing-circulation signature (no eddies, weak Hadley — §25/§26), not a
   closure-only bias; a calibration target only once topography/land/eddies exist.
4. **Slab-ocean equilibrium** (`logs/rce_slab5yr.nml`, `ocean_mode=1`, 10 m, 5 yr;
   264.8 s wall = **53 s/model-yr threaded**, confirming the throughput below).
   From the fixed-SST +5.46 W/m² start the slab warms, OLR climbs 157→174 W/m²,
   and TOA net crosses zero at ~yr 1 — then **overshoots to a quasi-steady net
   ≈ −13.5 W/m² by ~yr 3** (plateau, not drift). This is **not a clean
   energy-balanced equilibrium**, and the cause is diagnostic gold: cloud cover
   runs away **0.66 → 0.86** (heavily overcast) as free SST warms the tropics
   (surface-level T max ~307 K) and evaporation surges (LH 43→66 W/m²) — the
   near-saturated moist bias (item 3) means the cover diagnostic saturates *upward*,
   over-reflecting SW and pulling net negative. The freeze-floor slab (no sea ice)
   compounds it at high latitudes (surface-level T min ~212 K at 80°), preventing
   TOA closure. **The fixed-SST `cond_rh_crit=0.93` tuning does not survive slab
   coupling** — the cloud–SW feedback is too strong without subsidence drying. A
   real cloud-fraction scheme and sea ice (both listed below) are prerequisites for
   a physical slab equilibrium; this is not a re-tune of `cond_rh_crit`.
5. **Clear-sky OLR bias structurally gone** (`validate_era5 ... ecckd`): ecCKD
   **−5.07 W/m²** vs SESAM −7.07, with **no `LW_VAP_OPAC` fudge** (§20 fudge
   retired). Consistent offline (−5.1) ↔ coupled all-sky residual (+5.6, item 1).

**The radiation-cost fix — the real story of the session.** §28's "radiation runs
on the cadence, so the runtime share stays small" was **wrong**. A profile (T21L12,
single-thread, 10 model-days) showed radiation was **~90% of runtime**: dynamics
1.16 s, +non-rad physics 2.54 s, +clear-sky radiation 13.5 s, +cloud all-sky
24.7 s — despite the 3 h cadence. Root cause: the ecCKD reference k-table (an
80-iteration inverse-Gaussian bisection per g-point) was **rebuilt for every one of
the 2048 columns**, though it depends only on compile-time Malkmus parameters.

Fix (committed, bit-exact): build it once into a module cache (`aeros_ecckd_init`
primes it serially before the OpenMP region; `ensure_ktable` is the lazy fallback).
Full-run **~6× single-thread** (24.7 → 4.0 s / 10 d), radiation output
byte-for-byte unchanged, 19/19 green. The cadence default was also raised **3 h →
6 h** (`aeros_rad_class%interval`; near-linear speedup, coupled TOA net drift only
−0.15 W/m², within run-to-run noise); `rce_long` exposes `rad_interval`/`rad_scheme`
as optional overrides (`input/rce_defaults.nml`).

Measured throughput (T21L12, 10-core): **487 → 1656 model-yr/day** at 6 h (≈52 s
/model-yr), a **3.4×** wall-clock win. (Threaded scaling is only ~2.8× on 10 cores
now that radiation — the parallel-heavy term — is cheap; the dynamical core is the
floor. The design-target 5000 yr/day is not reached at this thread efficiency.)

A parallel single-pass all-sky-LW rewrite (clear+overcast sharing one transfer) was
prototyped and **rejected**: the k-table cache already collapses the run-twice cost,
so it gave no measurable full-run gain, and its `exp(a)·exp(b)` overcast is
roundoff-different from `exp(a+b)`, shifting the coupled day-100 net by +0.11 W/m².
The bit-exact run-twice path was kept.

**Bottom line:** ecCKD holds the coupled equilibrium together with no re-tune,
modestly improves both the clear-sky OLR bias (−7.1→−5.1) and the upper-trop warm
bias (−2 K), and is now the fast default. The revalidation is a *stability +
radiation-physics* confirmation under idealized boundary conditions — **not** an
ERA5 climate calibration, which is gated on the missing forcings (below).
