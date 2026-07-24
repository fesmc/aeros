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
