# SpeedyWeather.jl & SPEEDY vs. aeros — comparative study

A read-only comparison of **SpeedyWeather.jl** (Julia) and **SPEEDY / speedy.f90**
(the ICTP simplified GCM, Molteni/Kucharski) against **aeros**, aimed at two
current aeros pain points: **(a) numerical instability** (model-top over-cooling
→ unbounded thermal wind; a new prognostic-cloud scheme driving an explosive
near-equatorial jet) and **(b) speed**.

Scope note: this compares against aeros as it exists *today* (T21L12, **double
precision**, ∇⁶ hyperdiffusion, explicit Rayleigh+Newtonian sponge, ecCKD/SESAM
radiation, Sundqvist prognostic clouds), not the T42/L19/Float32 target sketched
in `docs/design.md`. Where the two disagree it is flagged.

Sources are collected inline as URLs and listed at the end.

---

## 1. Comparison table

| Axis | **SpeedyWeather.jl** | **SPEEDY / speedy.f90** | **aeros (today)** |
|---|---|---|---|
| Equations | Primitive eqns (±humidity), also shallow-water & barotropic | Hydrostatic primitive eqns | Hydrostatic primitive eqns |
| Spectral method | Spherical-harmonic transform, **grid-flexible** (own `SpeedyTransforms`, libsharp-lineage) | Spherical-harmonic transform (GFDL/Bourke core) | Spherical-harmonic transform via **SHTns** (external C lib) |
| Prognostic vars | **Vorticity–divergence**, T, ln pₛ, q | Vorticity–divergence, T, ln pₛ, q | **Vorticity–divergence**, T, ln pₛ; q off-spectral |
| Truncation (default) | **T31** (regular Gaussian), runs to T340/T511/T1023 | **T30** (96×48 Gaussian) | **T21** (dev); design target T42 |
| Grid | RingGrids: Gaussian, octahedral Gaussian, octahedral Clenshaw–Curtis, HEALPix | Full Gaussian 96×48 | Gaussian (SHTns); reduced/octahedral is a design target |
| Vertical | σ (sorted top→surface), default **L8**, equally spaced | σ, **L8** at 0.025…0.95 (top ≈25 hPa) | σ, **L12**, top ≈10 hPa |
| Time stepping | **Leapfrog + semi-implicit + Robert–Asselin–Williams** | **Leapfrog + semi-implicit + Robert–Asselin(–Williams)** | **Leapfrog + semi-implicit + RAW** (Williams 2009, α=0.53) |
| Hyperdiffusion | **∇⁸** (Laplacian power 4), implicit in spectral space; **adaptive** (boosts where \|vorticity\|>vor_max) and **power tapered to ∇⁴ in the stratosphere** | ∇⁸-class horizontal diffusion | **∇⁶** (power 3), implicit; fixed order at all levels |
| Model-top / stratosphere control | Diffusion-based: adaptive vor-scaling + stratospheric power tapering; **no explicit Rayleigh sponge** | Prescribed **zonal-mean stratospheric temperature relaxation** ("correction") + diffusion | Explicit **Rayleigh drag + Newtonian cooling sponge** above σ≈0.12 (off by default) |
| Precision | **Float32 default** (Float64/Float16 selectable — number-format-flexible) | Float64 (standard Fortran double) | **Float64 throughout, deliberately** (SHTns has no SP CPU path) |
| Parallelism | Julia **multithreading**: over levels (dycore), over gridpoints (physics); **no MPI**; **GPU in progress** (KernelAbstractions + Reactant) | Serial / modest OpenMP in forks | **OpenMP only**, over columns + (level,wavenumber); no MPI/GPU by design |
| Clouds | **Diagnostic** (RH-based), simple | **Diagnostic** (RH→cover), affects only LW "window" band | **Prognostic** cloud fraction (Sundqvist budget) — *the outlier* |
| Convection | Simplified mass-flux (Kucharski/Tiedtke-like) | Simplified mass-flux, SMSE closure | Betts–Miller-style |
| Radiation | Simple SPEEDY-style bands | **2 SW bands, 4 LW bands**, transmissivity fits | **ecCKD correlated-k** (default) + SESAM band LW; 6 h cadence |
| Throughput | ~**440–500 SYPD single-thread** at T31L8 Float32 (~3 min/yr); design.md cites ~2300 SYPD on one M4 core | ~**125–240 SYPD/core** at T30L8 | ~**1656 model-yr/day at T21L12 on 10 cores** (~165/core) |
| Language / niche | Julia; interactive, differentiable, ML-friendly research GCM | Fortran-77/90; cheap climate/DA workhorse | Fortran; **paleo, 10⁴–10⁶ yr coupled ice-sheet** runs, polar-priority |

Headline rows to remember: **both references keep clouds diagnostic and run L8**;
**SpeedyWeather's speed is a Float32 story built on its own transform**; and
**SpeedyWeather folds model-top stability into the diffusion operator rather than
bolting on a Rayleigh sponge.**

---

## 2. Top 5 differences that matter for aeros

1. **Precision is architectural, not a knob.** SpeedyWeather is Float32-by-default
   and its 4–5× per-core speed spread comes from precision + SIMD + cache, not
   from Julia. It can do this only because it *owns its transform*
   (`SpeedyTransforms`, number-format-flexible). aeros outsourced the transform to
   **SHTns, which has no single-precision CPU path** (`SHT_FP32` is GPU-only), and
   measured dp ~17% *faster* than sp precisely because of up/down conversion at
   every transform (`src/aeros_defs.f90`). So the single biggest lever behind
   SpeedyWeather's speed is **locked out for aeros as long as it depends on
   SHTns on CPU.** (JOSS paper; Klöwer et al. 2020/2022; `aeros_defs.f90`.)

2. **Model-top stability by diffusion vs. by sponge.** SpeedyWeather has **no
   classic Rayleigh sponge**. It stabilises the top two ways inside the
   hyperdiffusion: (i) **adaptive** diffusion that *increases where absolute
   vorticity exceeds `vor_max`* — i.e. it damps exactly the layers that are
   spinning up — and (ii) **tapering the Laplacian power from ∇⁸ to ∇⁴ in the
   stratosphere** (`power_stratosphere`, default 2), which makes near-top damping
   stronger and less scale-selective. aeros instead uses an explicit
   **Rayleigh-drag + Newtonian-cooling sponge** (`aeros_timestep.f90`, off by
   default). Both are valid; SpeedyWeather's is physically integrated and
   *state-adaptive*, which is directly relevant to aeros's two blow-ups.

3. **Clouds are diagnostic in both references — deliberately.** SPEEDY uses
   RH-based diagnostic cover that only modulates the LW *window*-band
   transmissivity; SpeedyWeather likewise uses simple diagnostic clouds. **Neither
   carries a prognostic cloud-fraction budget.** aeros's Sundqvist prognostic
   scheme (`aeros_cloud_prog.f90`) is a departure from *both references and from
   aeros's own design.md §6* ("Clouds: diagnostic initially") — and it is the
   source of the near-equatorial jet blow-up. This is a strong **over-engineering
   signal** (see §4a).

4. **Levels: references run L8; aeros runs L12 (design wants L16–20).** SPEEDY and
   SpeedyWeather both live comfortably at 8 σ-levels with tops around 25 hPa.
   aeros already exceeds them at L12/10 hPa. The references show a *few-level, low-top*
   configuration is a viable, stable regime — but they buy that stability with the
   diffusion/temperature-relaxation top treatment in (2), not with more levels.
   aeros's own design.md §4.1 argues the *opposite* case (more levels for the polar
   boundary layer); the references are the cheaper, simpler counter-example.

5. **Parallelism model is the same (threads, no MPI) — but their scaling ceiling is
   lower.** SpeedyWeather multithreads the dycore *over vertical levels*, so at L8
   its dycore parallelism caps near 8; physics threads over gridpoints. aeros
   threads over *columns and (level,wavenumber) pairs*, which is why it already
   sustains 10 cores at T21L12. aeros's OpenMP-over-columns design is arguably
   **better-suited to many cores** than SpeedyWeather's level-threaded dycore —
   this is a place aeros is already ahead, not behind.

---

## 3. What actually makes SpeedyWeather fast (and how much is reachable)

- **Float32.** The dominant lever. Single precision roughly halves memory traffic
  and doubles SIMD width; Klöwer et al. (2022) pushed ShallowWaters.jl to Float16
  for ~**4× on A64FX**. The whole SpeedyWeather number-format story rests on
  owning a format-flexible transform.
- **Own transform, in-place, ring-structured.** `SpeedyTransforms` +
  `RingGrids` + `LowerTriangularMatrices` give in-place, allocation-light
  transforms on iso-latitude rings, and let them choose **reduced/octahedral
  grids** (octahedral Gaussian, Clenshaw–Curtis) that cut gridpoints ~30–50% with
  negligible accuracy loss.
- **GPU (emerging).** As of 2025–26 SpeedyWeather is being made GPU-capable via
  **KernelAbstractions + Reactant** (JuliaCon 2026 talk "Towards a differentiable
  and GPU-capable GCM"); not yet the production path but the clear direction.
- **Multithreading**, no MPI, no distributed memory (JOSS, 2024).

**Reachability for aeros — the crux.** aeros's ~1656 yr/day at T21L12 on 10 cores
is already a *healthy* number (comparable to PlaSim T21 and to SpeedyWeather once
work is normalised for truncation/precision). The **Float32 advantage is the one
big lever aeros cannot pull through SHTns on CPU** — and aeros has already measured
that forcing sp *through* SHTns is a net loss. Therefore:

> SpeedyWeather's headline speed advantage is **largely a precision advantage that
> the SHTns-Float64 constraint locks out.** It is *not* reachable by tuning within
> the SHTns-on-CPU world. It becomes reachable only by **changing the transform
> story** — either (a) move the transform to GPU, where SHTns *does* expose a
> single-precision path (`SHT_FP32`), or (b) replace SHTns with a self-owned,
> format-flexible CPU transform (the exact move SpeedyWeather made). Both are large;
> (a) also fights the OpenMP-only, embedded-polar-module, long-run niche. So the
> honest reading is: **do not chase Float32 on CPU; treat it as the argument for a
> future GPU transform, and meanwhile take SpeedyWeather's *non-precision* speed
> ideas** (reduced/octahedral grids, in-place/allocation-light transforms, batched
> Legendre) which are compatible with SHTns-Float64.

---

## 4. Strategies aeros could adopt

### (a) Stability

**A1 — Adaptive, vorticity-scaled hyperdiffusion (highest value).**
Adopt SpeedyWeather's `vor_max` idea: *increase* hyperdiffusion in layers where
absolute vorticity exceeds a threshold. This targets **exactly** aeros's two
failure modes — the unbounded thermal wind at the top and the prognostic-cloud
near-equatorial jet — by damping the runaway *where and when it happens*, without
globally over-smoothing.
- *Benefit:* directly attacks both blow-ups; state-adaptive so it costs nothing
  in the quiescent regime.
- *Effort:* low–moderate — aeros's ∇⁶ hyperdiffusion is already implicit in
  spectral space (`aeros_timestep.f90`); add a per-level strength multiplier keyed
  to that level's vorticity maximum.
- *Fit:* pure win under SHTns-Float64/OpenMP; no new dependency; complements the
  existing sponge rather than replacing it.

**A2 — Stratospheric hyperdiffusion tapering (pairs with A1).**
Reduce the Laplacian power near the model top (SpeedyWeather: ∇⁸→∇⁴ above a
tapering σ). Lower order = stronger, less scale-selective damping at the top =
an *implicit* sponge that cannot spin an unbounded thermal wind.
- *Benefit:* a physically-integrated alternative/complement to the explicit
  Rayleigh sponge for the over-cooling top.
- *Effort:* low — make `ndiff` (currently a single order) σ-dependent.
- *Fit:* clean under the current architecture; may let the explicit sponge be
  gentler or off.

**A3 — Keep clouds diagnostic; fix the runaway by weakening cloud→radiation
coupling, not by adding a prognostic budget.** Both references avoid diagnostic-cloud
runaway not with a prognostic scheme but by making the **cloud-radiation coupling
weak/bounded** — SPEEDY lets cloud cover modulate *only* the LW window-band
transmissivity. aeros's `aeros_cloud_prog.f90` header shows the prognostic scheme
was introduced to cure a diagnostic RH→cover runaway (0.66→0.86 in coupled RCE);
but the cure introduced the equatorial-jet blow-up.
- *Benefit:* removes the blow-up's source entirely; matches both references and
  aeros's own design.md §6.
- *Effort:* low if reverting; the machinery already exists behind a `scheme`
  selector.
- *Fit / over-engineering flag:* **the references deliberately stayed diagnostic.**
  At T21–T31 a full Sundqvist prognostic cloud-fraction budget is very likely
  over-engineering. Recommended path: keep clouds diagnostic, and if the RH
  runaway returns, damp the *feedback* (SPEEDY-style window-only coupling, or a
  bounded/relaxed cover as CLIMBER-X does with `0.1·new+0.9·old`) before reaching
  for a prognostic budget.

**A4 — Prescribed stratospheric temperature relaxation (SPEEDY's approach).**
SPEEDY nudges the stratosphere toward a prescribed zonal-mean temperature. aeros's
Newtonian-cooling sponge (`sponge_tref=216 K`) is already a cousin of this; SPEEDY's
version is a documented, stable recipe for a low-top model. Note design.md §5's
caveat: a *prescribed* stratospheric correction blocks CO₂ stratospheric cooling —
so keep it as a stability aid, not a radiative substitute, and keep it out of the
ecCKD path.
- *Benefit:* proven cure for low-top over-cooling.
- *Effort:* low (extends existing sponge).
- *Fit:* fine for stability runs; do not let it contaminate paleo CO₂ sensitivity.

### (b) Speed

**S1 — Reduced / octahedral Gaussian grid (best speed win under the constraints).**
SpeedyWeather's grid-flexibility (octahedral Gaussian / Clenshaw–Curtis) is a
non-precision speedup: ~30–50% fewer gridpoints for negligible accuracy loss, and
it removes redundant polar columns *without* touching polar spectral resolution —
ideal for a polar-priority model. design.md §4 already calls for this.
- *Benefit:* ~30% off gridpoint physics + transform input; disproportionately
  removes polar redundancy.
- *Effort:* moderate — SHTns supports reduced grids; needs the grid plumbing.
- *Fit:* fully compatible with SHTns-Float64 and OpenMP; strong niche match.

**S2 — Allocation-light, in-place, batched transforms.** SpeedyWeather's speed
also comes from in-place transforms and batching all levels of a wavenumber into
one BLAS call. design.md §4 already prescribes "all levels of a wavenumber batched
into one DGEMM" — confirm it is implemented and that per-step transform allocations
are minimised.
- *Benefit:* better cache/SIMD use; batching improves with more levels (helps the
  L12→L16–20 direction).
- *Effort:* moderate; an implementation-quality task, not an architecture change.
- *Fit:* native to SHTns-Float64/OpenMP.

**S3 — Radiation cadence + per-step surface-flux rescaling.** aeros already calls
radiation on a 6 h cadence; SPEEDY/IFS confirm this is the right lever and that the
key is to **rescale surface LW/SW every step** from skin T and albedo (design.md
§5). Already largely in place — verify the per-step rescale exists and is cheap.
- *Benefit:* radiation ≈ a few % of runtime while keeping the diurnal/surface
  response; matters most over ice sheets (SMB).
- *Effort:* low — mostly a verification/tuning item.
- *Fit:* native.

**S4 — GPU transform as the *only* route to Float32-class speed (long-term, flagged
not recommended now).** The Float32 speedup is unreachable on CPU through SHTns.
The *sole* way to capture it without abandoning SHTns is the GPU path, where SHTns
exposes `SHT_FP32`. This fights OpenMP-only and the embedded-polar-module design,
and SpeedyWeather's own GPU work is still maturing (2025–26).
- *Benefit:* potentially the 2–4× precision speedup — but only on GPU.
- *Effort:* very high; a strategic pivot.
- *Fit:* **fights** the current niche/architecture. Treat as a future option to
  weigh against a self-owned transform, not a near-term task.

**S5 — Ensemble-parallelism over core-parallelism for long paleo work.** Matches
both design.md §4.3 and SpeedyWeather's threading ceiling: past ~32 effective
threads, one-member-per-run beats one big parallel run. Not a SpeedyWeather import
per se, but their level-threaded dycore ceiling is empirical support for it.

---

## 5. Over-engineering signals (things the references did *simpler*)

- **Prognostic cloud fraction** — both references diagnostic; aeros is the outlier
  and it is the current blow-up. (§4a A3.) **Strongest signal.**
- **L12 / 10 hPa top** — references run L8 / ~25 hPa and are stable; aeros's extra
  levels are justified by the *polar-BL* argument (design.md §4.1), not by dynamical
  necessity. Fine, but it is more than the references need.
- **ecCKD correlated-k radiation** — far more sophisticated than SPEEDY's 2 SW / 4
  LW band scheme. Justified by aeros's paleo-CO₂ niche (design.md §5 rejects
  SPEEDY's constant-α CO₂), so *not* over-engineering — but worth noting aeros is
  deliberately heavier here.
- **Explicit Rayleigh+Newtonian sponge** — SpeedyWeather achieves top stability
  with adaptive diffusion instead. aeros's sponge is not wrong, but the adaptive-
  diffusion route (A1/A2) may be a cleaner primary mechanism.

---

## 6. Sources

- Klöwer et al. (2024), *SpeedyWeather.jl: Reinventing atmospheric GCMs towards
  interactivity and extensibility*, JOSS 9(98):6323 —
  https://doi.org/10.21105/joss.06323 ;
  https://www.navidconstantinou.com/publications/speedyweather.pdf
- SpeedyWeather.jl repository — https://github.com/SpeedyWeather/SpeedyWeather.jl
- SpeedyWeather.jl documentation (functions index / HyperDiffusion, vertical
  coordinates) — https://docs.juliahub.com/General/SpeedyWeather/stable/functions/ ;
  https://speedyweather.github.io/SpeedyWeatherDocumentation/dev/
- SpeedyWeather.jl performance (~500 SYPD single-thread T31L8) —
  https://juliapackages.com/p/speedyweather
- SpeedyWeather.jl GPU direction (KernelAbstractions + Reactant), JuliaCon 2026 —
  https://pretalx.com/juliacon-2026/talk/GB8WXW/
- Klöwer, Düben & Palmer (2020), *Number Formats, Error Mitigation, and Scope for
  16-bit Arithmetics…*, JAMES 12(10) — https://doi.org/10.1029/2020MS002246
- Klöwer, Hatfield et al. (2022), *Fluid Simulations Accelerated With 16 Bits…*,
  JAMES 14(2) — https://doi.org/10.1029/2021MS002684
- Molteni (2003), *Atmospheric simulations using a GCM with simplified physical
  parametrizations (SPEEDY)*, Clim. Dyn. 20 — https://doi.org/10.1007/s00382-002-0268-2
- Kucharski et al. (2013), *On the need of intermediate-complexity GCMs: a SPEEDY
  example*, BAMS 94 — https://doi.org/10.1175/BAMS-D-11-00238.1
- Molteni & Kucharski, *Description of the ICTP AGCM (SPEEDY) v41* —
  http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf ;
  http://users.ictp.it/~kucharsk/speedy-doc.html
- speedy.f90 (modern Fortran rewrite, S. Hatfield) —
  https://github.com/samhatfield/speedy.f90
- SPEEDY-NEMO coupled ICM — https://doi.org/10.1007/s00382-023-07097-8
- aeros internal: `src/aeros_defs.f90` (precision), `src/dynamics/aeros_timestep.f90`
  (∇⁶ hyperdiffusion, sponge, RAW), `src/physics/aeros_cloud_prog.f90` (Sundqvist),
  `src/physics/aeros_radiation.f90` (SESAM/ecCKD), `docs/design.md` (§3.4 anchors,
  §4 numerics, §4.1 levels, §5 radiation, §6 physics).
