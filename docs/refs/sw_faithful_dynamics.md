# SpeedyWeather-faithful dynamical core: sustaining baroclinic eddies at T21

`hadley_edge_emf.md` pinned the residual subtropical moist bias to aeros having
**no eddy momentum flux** (SpeedyWeather ~3.9 m²/s²) — its RCE cannot hold an
eddying state: seed it and the eddies grow then NaN at the model top by ~day 27.
This note records the fix: reproducing SpeedyWeather's dynamical-core numerics.

## What actually blows up (diagnosis)

Seeded eddies grow the *right* structure — an upper-level subtropical jet marching
to 30° — then the **balanced/rotational eddy KE explodes super-exponentially**
(doubling per day at KE~1 m²/s²) and NaNs. It is a **leapfrog computational mode**
in the balanced field:

- The semi-implicit **gravity-wave** treatment does not touch it. α=1
  backward-implicit (`si_alpha`, below) changes the divergent circulation
  (ω differs ~3 hPa/day) but the seeded run still dies ~day 26.
- Reproducing SW's whole diffusion recipe (∇⁸/4h + divergence ∇⁸→∇⁴/1h) does not
  fix it either; `couple_diabatic` makes it *worse* (a convective computational
  mode at low filter).
- Only the **Robert time filter** moves it — and aeros ran it too weak.

## The unlock: SpeedyWeather's Robert filter coefficient

SpeedyWeather's leapfrog uses `robert_filter = 0.1`, `williams_filter = 0.53`
(`SpeedyWeather/src/time_stepping/leapfrog.jl`). aeros's `rce_allfix` ran
**`eps_filter = 0.06`** — too weak, so the computational mode blows up; the only
value tried before was 0.15, which is too strong and **damps the eddies to death**
(EMF→0). **SW's actual 0.1 is the sweet spot**: the run is stable *and* the eddies
survive (KE ~5.9 m²/s²). `raw_alpha = 0.53` already matches SW's Williams α.

## The SpeedyWeather-faithful configuration (all opt-in; defaults bit-for-bit)

| knob | value | SpeedyWeather |
|---|---|---|
| `si_alpha` | 1.0 | backward-implicit (α=1) |
| `ndiff`, `tau_diff` | 8, 4.0 | ∇⁸, 4 h (vorticity+temperature) |
| `tau_diff_div` | 1.0 | 1 h divergence diffusion |
| `diff_taper`, `diff_ndiff_top`, `diff_taper_sigma` | T, 4, 0.2 | divergence ∇⁸→∇⁴ above σ=0.2 |
| `l_sponge` | F | SW has no sponge |
| `eps_filter`, `raw_alpha` | 0.1, 0.53 | Robert 0.1, Williams 0.53 |

## Result (T21 RCE aquaplanet, 100 d, seeded, stable)

| metric | baseline aeros | SW-config, no eddies | SW-config + eddies | SpeedyWeather |
|---|---|---|---|---|
| subtropical descent peak | +1.9 @ 36° | +1.7 @ 25° | **+3.15 @ 25°** | +6.4 @ 30.5° |
| descent span (ω>1, σ~0.5) | 25–69° | 25–30° | **19–30°** | 25–36° |
| jet | 14 @ 25° | 40 @ 30° | **31 @ 30°** | ~30–44 |
| subtropical free-trop RH | 94% | 87% | **78%** | 62% |
| eddy KE (mid-trop) | 0 / blows up | 0 | **5.9 (stable)** | ~few |

The SW config *itself* (even axisymmetric) concentrates the descent to 25–30° and
builds a 40 m/s jet at 30°; the **eddies** then nearly double the descent strength
(+1.7→+3.15), dry the subtropics further (87→78%), and brake the jet (40→31, a real
eddy-momentum feedback). Both push toward SW/ERA5.

## Open (the gap that remains)

Not all the way to SW yet: descent +3.15 vs +6.4, RH 78 vs 62%. **Verified stable
over 300 days** — eddy KE equilibrates at ~6–7 m²/s² (not a pre-blowup transient).

The gap is now characterized: with proper averaging (150-day mean) the eddy
momentum flux is **|[u*v*]| ~0.52 m²/s², coherent but ~7× weaker than SW's 3.9**
(the earlier ~0.04 was transient cancellation in a too-short average). So aeros's
eddies have **SW-like energy (KE~6) but under-transport** — they flux momentum, just
inefficiently. That 7× flux deficit is why the descent reaches +3.15 not +6.4.

**Root of the flux deficit — the eddies are mis-scaled (eddy KE spectrum, both
models T21, `accum_espec` in rce_long / the SW extraction):**

| m | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| SpeedyWeather KE% | 24 | 18 | 12 | 15 | 11 | 11 | 5 | **2** |
| aeros KE% | 10 | 9 | 4 | 22 | 6 | 2 | 2 | **43** |

SpeedyWeather's eddy energy is smooth and synoptic (m=1–5) and it transports (flux
carried by m=2–5). aeros locks **43% of its eddy energy in an m=8 spike that SW
barely has (2%)** — and m=8 carries ~no momentum flux (aeros's flux comes from the
weaker m=1–4). So aeros's eddies are not low on total energy, they are **mis-scaled**:
half the energy sits in a small mode SW lacks, starving the transporting synoptic
scales. T42 is not the fix (prior null result; SW works at T21). The seed only
perturbs m=1–6, so **m=8 is naturally selected by aeros's flow**, not seeded.

Open: *why* does aeros select an m=8 mode SW doesn't? Numerical (a spurious mode)
vs physical (the narrow 31–40 m/s jet making m=8 the most-unstable scale). Killing
the m=8 spike / shifting energy to synoptic scales should close the descent/RH gap.

## Code (this branch)

- **`si_alpha`** — semi-implicit decentering weight in `aeros_semiimp`
  (`X_a = alpha X^new + (1-alpha) X^old`): 0.5 centered (default, bit-for-bit),
  1.0 backward-implicit. Guarded by `test_semiimp` at both 0.5 and 1.0.
- **`tau_diff_div`** + divergence-only order taper in `aeros_timestep`'s `diffuse`
  — vorticity/temperature keep the untapered full order; divergence gets the
  shorter timescale and the σ-tapered order (SW's gravity-wave absorber).
- Both default off → the 25 acceptance tests are bit-for-bit unchanged.
