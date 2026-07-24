# M0a — measured cost of the dynamical core

Machine: Apple M5, 10 cores, gfortran 15, macOS. Run with
`make bench && libaeros/bin/bench_m0a.x par/bench_m0a.nml`.
Raw data: `logs/bench_m0a.csv` (45 cases per precision).

**Read this alongside the caveats in §4 before quoting any number.** The
headline result is large enough to change design decisions, and most of that
change survives the caveats — but not all of it.

---

## 1. The number M0a existed to produce

`design.md` §11 records that the central figure in the plan was never measured:

> T42L19 throughput (~380 SYPD/core) is an **extrapolation**, not a
> measurement: SpeedyWeather's T31L8 2,300 SYPD scaled by (42/31)³×(19/8) ≈ 6.
> M0a exists to replace this estimate with a number.

Measured, single core, T42**L20** (one level *more* than the estimate):

| | design.md estimate | measured (sp) | measured (dp) |
|---|---|---|---|
| SYPD / core | ~380 | **1,803** | **2,166** |
| core-s / model yr | 227 | **47.9** | **39.9** |

**The extrapolation was 4.7–5.7× too pessimistic.** The error is in the
scaling, not the anchor: the (T/T)³ law it used is too steep in this range
(§3 below).

Same comparison across the table (dp, 1 thread, vs §3.6's estimates):

| config | §3.6 estimate | measured | ratio |
|---|---|---|---|
| T31L16 | ~57 | **20.8** | 2.7× cheaper |
| T42L19/L20 | 227 | **39.9** | 5.7× cheaper |
| T85L19/L20 | ~3,600 | **245.8** | 14.6× cheaper |

## 2. What this does to the central bind

§3.6 concludes that *"T42L19 does not fit a 5,000 yr/day coupled target"* at
82% of a 276 core-s/yr budget, and §9 risk 1 calls this "the resolution bind":
the truncation we can afford (T31) is the one Löfverström & Liakka show fails
to build the Eurasian ice sheet. §3.7's multi-resolution correction exists
primarily to escape it.

Against the same 276 core-s/yr budget, using measured numbers:

| config | §3.6 said | measured, 8 transforms/level | corrected, 11 (§4.2) |
|---|---|---|---|
| T31L16 | 21% | 7.5% | **10%** |
| T42L19/L20 | 82% ⚠ | 14% | **19%** |
| T85L19/L20 | 1300% ✗ | 89% | **123%** |

**On this evidence the cost argument for T31 over T42 disappears.** T42 costs
~19% of the coupled budget, leaving ~81% for ocean, ice, land and physics
where §3.6 expected ~18%.

Two things follow, and they should not be conflated:

- **§3.7 is no longer *forced by cost*.** The bind it was invented to escape is
  much weaker than the plan assumed.
- **§3.7 may still be *wanted for accuracy*.** L&L's result is that T42 fails
  Eurasia and T85 succeeds; that is a statement about physics, not budget, and
  nothing here touches it. If anything §3.7 gets *more* attractive, because the
  T85 reference is far cheaper than budgeted (§3 below).

## 3. Scaling laws

**Truncation scales well below cubic.** Measured at L20, 1 thread, dp:

| step | measured | cubic (T/T)³ |
|---|---|---|
| T31 → T42 | **1.84×** | 2.48× |
| T42 → T85 | **6.16×** | 8.31× |

This is the source of the §1 discrepancy. At these truncations the O(N³)
Legendre term has not yet overtaken the O(N² log N) FFT, so cubic scaling
overcharges.

Consequence for §3.7: it budgets the T85 reference at *"T85L19 ≈ 16× T42L19 ≈
63× T31L16"*. Measured, **T85L20 = 6.2× T42L20 = 11.8× T31L16** — the
correction's diagnostic runs are ~5× cheaper than planned. §3.7 estimates two
years of T85 per 500 years of transient at ~25% of the T31 atmosphere cost;
measured, that is ~5%.

**Levels scale linearly**, as §4.1 assumed. T42, L8 → L20 costs 2.68× for 2.5×
the levels (§4.1 estimated L19 ≈ 2.4× L8). There is no reason to economize on
vertical resolution, which is the conclusion §4.1 reached on physical grounds.

**Transforms dominate.** Single-threaded, the semi-implicit solve is 3.0–13.5%
of a step and transforms are the rest. Optimization effort belongs in the
transform layer, as designed — though this ratio will shift once column physics
exists (§4).

The solve's share *rises* with thread count, to 24% at T31L8 on 10 threads.
That is not the solve getting slower: it is the transform loop hitting the
threading ceiling (§4 caveat 4) while the solve — which parallelizes over
`nlm` ≈ 528 coefficients rather than over 32–80 transforms — keeps scaling.
Worth remembering when the same sweep runs on 32 cores, where the two parts
may well cross over.

## 4. Caveats — what these numbers are not

1. **Dynamics only.** No radiation, condensation, surface tiles, ocean, ice or
   land. §3.6 budgets ~1% for physics and ~4% for the polar nests, so there is
   ample headroom, but those shares are themselves estimates.
2. **`ntr_per_level = 8` was an assumption, and M1.3 has now measured it: the
   real count is 11.** The tendency evaluation synthesizes (u,v) [2
   scalar-equivalents], ζ [1], D [1], T [1] and ∇T [2], and analyses A [2],
   Φ+K [1] and dT/dt [1]. Cost scales linearly in this, so **every core-s/yr
   figure above should be multiplied by 11/8 (+37%)** and every SYPD figure
   divided by it.

   Redoing §2 with the measured count: T31L16 goes 7.5% → **10%** of the
   coupled budget and T42L20 goes 14% → **19%**. The conclusion is unchanged —
   §3.6 expected 82% for T42 — but the margin is smaller than first reported.

   SHTns' combined `SHqst_to_spat` could deliver ζ alongside (u,v) in one call
   and is documented as "significantly faster" than separate transforms, which
   would bring the count to 10. Worth measuring once the core runs.
3. **`dt = 1800 s` is held fixed across truncations**, which flatters T85: a
   T85 run needs roughly half the timestep of T42 on CFL grounds. Adjusting
   `dt ∝ 1/T` gives T85L20 ≈ **490** core-s/yr rather than 246, i.e. ~178% of
   the 276 budget rather than 89%, and T85/T42 ≈ 12× rather than 6.2×. **The
   T42 conclusion in §2 is unaffected** — T42 is benchmarked at the ~30 min
   step §4 specifies for it — but every T85 figure above should be read as an
   optimistic bound.
4. **Threading saturates at 4 threads on this machine** and gets no better to
   10. That is the M5's performance/efficiency-core split, not a statement
   about a homogeneous 32-core node. **The scaling question §3.6 and §9 risk 2
   turn on is still open and still needs the HPC.** Extend `threads` in the
   namelist to 16, 32, 64 there.
5. Untuned vs tuned matters: all results use SHTns' tuned path
   (`quick=.FALSE.`), which is ~20% faster than the quick-init path the
   acceptance tests use.

## 5. Precision: dp is faster than sp

`design.md` §4 calls for "Float32 throughout the core". Measured, single
thread, step time as a ratio dp/sp:

| | L8 | L16 | L20 |
|---|---|---|---|
| T31 | 1.29 | 1.04 | 0.87 |
| T42 | 0.83 | 0.83 | **0.83** |
| T85 | 0.89 | 0.79 | **0.83** |

**At every production-relevant configuration, dp is ~17% faster than sp.**

The mechanism is the one flagged when the scaffold was built: SHTns has no
single-precision CPU path, so an sp build converts the grid-side array up to
double at every transform and back down afterwards. With transforms at 84–97%
of the step, that copy costs more than single precision saves on the rest.

sp's advantage is real but lives in the *column physics* — radiation,
condensation, snowpack — which never touch SHTns and which §3.6 puts at ~1–5%
of budget. So sp optimizes a twentieth of the work while taxing the rest.

**Recommendation: work in dp internally**, and re-measure at M2 when physics
exists.

**Adopted (AR, 2026-07-24).** `wp = dp`, fixed in `aeros_defs`; the
`precision=` make switch and the sp specifics in `aeros_spectral` are gone, so
the conversion path has left the hot loop entirely. Single precision remains as
`wp_ext` — the coupling-boundary kind, matching yelmo and fesm-utils — with the
conversion to happen at the facade when coupling arrives at M4. A physics module
that wants sp locally may still declare it locally; that is a separate decision,
to be made with a measurement at M2.

## 6. Incidental findings

- **SHTns tuning costs 8.2 s per config at T85** and is not remembered between
  configs, so a 10-thread pool paid it ten times. `SHT_LOAD_SAVE_CFG`
  (`cache=.TRUE.`) fixes it: 8.2 s → 3.7 ms on the second and later configs.
  The cache is per working directory, and runme stages each simulation into a
  fresh one, so a cold production run still pays it once — negligible over a
  paleo integration, dominant over a short test.
- **`OMP_PROC_BIND` / `OMP_PLACES` abort under macOS libgomp** ("Affinity not
  supported on this configuration", exit 144). They are wanted on the HPC and
  must not be set here.
- **The bundled `shtns.f03` is stale**: it declares `shtns_set_batch`, while
  the library exports `shtns_set_many`. Anything else in that header may be
  equally out of date. Batched transforms remain untested and are the one
  plausible route to making SHTns' own threading useful — see
  `aeros_sht_pool_class`.
