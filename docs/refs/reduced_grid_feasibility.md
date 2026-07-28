# Reduced / octahedral Gaussian grid for aeros — feasibility & design memo

A read-only feasibility spike on strategy **S1** in `docs/refs/speedy_comparison.md`
(reduced/octahedral Gaussian grid for speed). Scope: determine whether aeros can
adopt it with the **current transform library (SHTns)**, map the blast radius,
estimate the honest payoff for aeros' transform-dominated cost profile, and give a
GO / NO-GO. **No code was changed.**

**Bottom line up front: NO-GO.** Two independent blockers, either sufficient on its
own: (1) **SHTns cannot represent a reduced grid** — its grid contract is a single
scalar `nphi` shared by every latitude ring; (2) even if it could, the payoff is
**~1–5%, not ~30%**, because aeros spends **84–97% of the step in the spherical-
harmonic transform**, and a reduced grid's transform savings are only realizable
*inside* the transform library, which SHTns will not do. The better near-term speed
lever attacks the 84–97% directly: **cut the number of transform calls per step and
batch levels** (§5).

---

## 1. SHTns verdict — NOT SUPPORTED

**SHTns 3.7.5 (the version aeros links, `fesm-utils/SHTns/shtns-3.7.5/`) is a
regular-grid-only library. A reduced/octahedral grid is impossible through its API,
not merely un-plumbed in aeros.** Evidence, from the shipped source:

1. **The grid is one scalar `nphi`, not a per-ring array.** The public struct
   (`shtns.h:80–93`) exposes:
   ```c
   const unsigned int nlat;   // latitudes
   const unsigned int nphi;   // longitudes  -- ONE number for the whole grid
   const unsigned int nspat;  // real numbers in a spatial field
   ```
   and `nspat == nlat*nphi` (enforced; aeros itself asserts this in
   `aeros_sht_init`, `aeros_spectral.f90:236`). A reduced grid needs `nphi(j)`
   per ring `j`. There is no such array anywhere in the struct or the API.

2. **Every grid-setup entry point takes scalar `(nlat, nphi)`:**
   `shtns_init(..., int nlat, int nphi)`, `shtns_set_grid(cfg, flags, eps, nlat,
   nphi)`, `shtns_set_grid_auto(cfg, ..., int *nlat, int *nphi)` (`shtns.h:155–161`).
   None accepts a longitude-count vector.

3. **The grid-type enum offers only regular grids** (`shtns.h:57–65`):
   `sht_gauss`, `sht_auto` (= gauss), `sht_reg_fast`, `sht_reg_dct`,
   `sht_quick_init`, `sht_reg_poles` (Clenshaw–Curtis), `sht_gauss_fly`. All are
   **uniform-`nphi`** grids — Gauss or equispaced. There is no reduced/octahedral/
   thinned type.

4. **The longitude geometry is hard-wired uniform.** The coordinate macros assume
   a single spacing: `PHI_RAD(shtns, ip) = (2π/(nphi*mres))*ip` (`shtns.h:129`).
   The longitudinal step is an FFT over a fixed `nphi` per ring.

5. **The README lists exactly two grid families** (`shtns-3.7.5/README.md:21`):
   "support for **regular grids** (but they require twice the number of nodes than
   Gauss grid)" — i.e. Gauss + equispaced, both rectangular. No reduced grid.

6. **Zero occurrences** of `reduced`, `octahedral`, `nphi[`, per-ring or
   variable-longitude logic in the entire C source (`grep` over all `*.c`/`*.h`).

7. **aeros already documents this** in `aeros_spectral.f90:411–417`:
   > "This is the FULL Gaussian grid. docs/design.md section 4 wants a
   > reduced/octahedral grid (~30% fewer gridpoints), which SHTns does not
   > support: its longitudinal FFT assumes a constant nlon per latitude."

   This spike confirms that comment against the upstream source and header.

### The critical distinction: SHTns ≠ libsharp

A web search for "SHTns reduced grid" returns text about "variable pixels per
iso-latitude ring" — **that describes libsharp / ducc0** (Reinecke & Seljebotn
2013, *Libsharp*), a *different* library, and HEALPix. **SpeedyWeather.jl does not
use SHTns.** Its grid flexibility comes from its own `SpeedyTransforms` (libsharp
lineage), which is exactly why it can do octahedral Gaussian, Clenshaw–Curtis and
HEALPix. The comparison doc (`speedy_comparison.md` §1, §3) is right that
*SpeedyWeather* owns a grid-flexible transform; the error to avoid is assuming that
capability transfers to aeros through SHTns. **It does not.** The libraries that do
reduced grids (libsharp, ducc0) are the ones SpeedyWeather built on; adopting one
of them *is* "replace the transform," which is out of scope here and is the same
move `speedy_comparison.md` §3 flags as large.

**Conclusion:** S1 cannot be done on CPU without replacing SHTns. That alone is a
NO-GO for the stated scope ("if SHTns cannot do reduced grids on CPU, that likely
kills S1"). Sections 2–4 show that even setting SHTns aside, the payoff would not
justify the work.

---

## 2. Blast-radius map

aeros stores every field as a dense `(nlon, nlat, nlev)` array keyed off a single
`grd%nlon` (`aeros_defs.f90:201–225`). A reduced grid makes longitude ragged
(variable per row). Mapping the assumption, file by file:

| File / area | What assumes uniform `nlon` | Severity | Notes |
|---|---|---|---|
| `aeros_spectral.f90` | The SHTns config itself — `nspat=nlon*nlat`, scalar `nphi`, uniform `lon(i)=2π i/nlon`, `surface_integral` uses one `dlon` | **BREAKING** | The transform *is* the grid contract. No aeros-side change can make SHTns ragged. Root blocker. |
| `aeros_grid.f90` | `area = w_j·dlon·a²` with one `dlon`; `coriolis`, `lon(:)` all `(nlon,...)` | **Breaking→moderate** | Area/Coriolis become per-row `dlon(j)`; conceptually easy *if* the transform allowed it, but it doesn't. |
| `aeros_defs.f90` (`aeros_grid_class`) | `lon(nlon)`, `area(nlon,nlat)`, `coriolis(nlon,nlat)`, all state arrays `(nlon,nlat,nlev)` | **Breaking** | Either go ragged (drop dense arrays) or pad + carry `nlon_of(j)` valid-count. Touches every consumer. |
| `dynamics/aeros_moisture.f90` (FV tracer transport) | Single `dlam=2π/nlon`; **zonal sweep** wraps `ip=i+1; if(ip>nlon)ip=1` per row; **meridional sweep** exchanges N–S fluxes between *aligned* cells of equal-length rows; `nsub` from worst-case polar Courant | **BREAKING (hardest)** | Zonal sweep with per-row `nlon(j)`/`dlam(j)` is tractable. The **meridional sweep between rings of different `nlon` needs a conservative longitude remap at every ring interface** — a new scheme, not a tweak. This is where "positive-definite + conservative + constant-preserving" (the module's three proven properties) must be re-derived. High risk to the mass budget. |
| `physics/*.f90` (radiation, convection, condensation, surface, vdiff, cloud, land, ocean, held_suarez) | Every module loops `do j; do i=1,nlon` | **Moderate (but pervasive)** | Column physics is per-column, so ragged is *semantically* fine — but every loop bound and every `(nlon,nlat)` work array changes to `nlon_of(j)`. Mechanical, ~a dozen files. |
| `aeros_io.f90` | netCDF write is rectangular: `count=[nlon,nlat,nlev,1]` (`:177–193`) | **Moderate** | netCDF/CF is inherently rectangular. Output must either pad-to-full (defeats the point on disk) or write a ragged/CF-unstructured layout that breaks every downstream analysis/plotting tool. |
| `aeros_diagnostics.f90`, `aeros_budget.f90` | Zonal means, global integrals via Gauss quadrature over `(nlon,nlat)` | **Moderate** | Global integral generalizes to per-row `dlon(j)` cleanly; zonal-mean plots and spectra need per-row handling. |
| `aeros_state.f90`, `aeros_tendency.f90`, `aeros_vordiv.f90`, `aeros.f90` (driver) | Allocate/loop `(nlon,nlat,nlev)`; grid↔spectral round-trips | **Breaking** | The tendency core round-trips spectral↔grid every step through SHTns on the *full* grid. A reduced physics grid would require **full↔reduced interpolation on every transform boundary** (see §3). |

**Is "padded rectangular array + per-row valid count" a viable compromise?**
Partly. It keeps arrays dense and lets the *physics and grid-space dynamics* skip
redundant polar points. But it **cannot help the transform**: SHTns must be handed
a full `nlat*nphi` rectangle, so before every `spat_to_SH`/`SH_to_spat` you must
scatter the padded-reduced field back to a full uniform grid (and gather after).
That interpolation runs at every transform boundary — 11 scalar-equivalent
transforms × nlev per step (`m0a_results.md` note 2) — and its cost is charged
against a savings that is already tiny (§3). The compromise reduces the *physics*
bill (already ~1–5%) while **adding** an interpolation bill to the *transform* path
(84–97%). That is the wrong side of the ledger.

---

## 3. Expected speedup — the honest number is ~1–5%, not ~30%

The generic "reduced grids cut 30–50% of gridpoints" figure is **irrelevant to
aeros' cost structure**, for a specific and decisive reason.

**aeros' measured profile:** spherical-harmonic transforms are **84–97% of the
step** (`m0a_results.md:154`; consistent with design.md §3/§4 and the `aeros_sht_*`
being "the performance kernel"). Column physics + grid-space dynamics are the
remaining **~3–16%**.

**Why the reduced-grid win does not reach the transform:** a reduced grid speeds up
a transform only if the transform library **skips the unresolved high-`m` Legendre
work on polar rings** — i.e. makes the (order `m` × latitude) work matrix
triangular instead of rectangular. That is precisely what libsharp/SpeedyWeather
do, and it is *internal to the transform*. SHTns has no such path; its Legendre
kernels are built around the fixed `nlat×nphi` rectangle. So even if aeros fed it
reduced data, SHTns would do the full rectangular transform anyway. **The transform
— 84–97% of the cost — sees no benefit.**

What a reduced grid *can* speed up under a fixed SHTns:

- **Column physics** (~1–5% of step): ~30% fewer columns → **~0.3–1.5% of the
  step** saved.
- **Grid-space dynamics** (FV moisture transport, gradient/energy assembly; part
  of the ~3–16%): a fraction of ~30% → optimistically **another ~1–3%**.
- **Minus interpolation overhead:** full↔reduced scatter/gather on every transform
  boundary, charged against the dominant 84–97% path. Plausibly **eats most or all
  of the above.**

**Net, at T21L12 and T42L19 alike: ~1–5% gross, quite possibly net-negative**,
because the two truncations aeros targets are exactly the regime where transforms
dominate. (The gridpoint-reduction win grows with truncation — at T1279 IFS,
gridpoint work dominates and 30% is real — but aeros lives two-to-three truncation
decades below that crossover.) Note also aeros is **already on the smallest
quadratically-unaliased Gauss grid** (`3T+1`, `aeros_sht_grid_size`,
`aeros_spectral.f90:399–436`): there is no grid-oversizing to reclaim.

**Verdict on payoff: the win is in the noise for aeros' profile, and the blast
radius (§2) is a multi-week refactor of the FV moisture core and every physics
loop. The ratio is indefensible.**

---

## 4. Niche fit — orthogonal, neither helps nor hurts

design.md §2 describes a **global coarse core + embedded polar fine-grids** (real
~20 km polar-stereographic nests, §2.3), and §3.2 motivates the spectral core by
**triangular-truncation isotropy** (effective resolution at 85°N = equator).

- A reduced grid removes redundant polar **longitude gridpoints**; it does **not**
  touch **spectral** resolution. Triangular-truncation isotropy is a property of
  the coefficient set, not of `nphi(j)`, so §3.2's polar argument is unaffected.
- The embedded polar modules are **separate high-resolution nests**, coupled by
  conservative remap (§2.1). They are orthogonal to how the global core samples
  longitude. A reduced global core would neither fight nor help them.

So S1 does not conflict with the polar-priority design — but the design's own
selling point ("disproportionately removes redundant polar columns", design.md §4)
is a hollow benefit here: those polar columns are **physics columns**, which are
~1–5% of the budget. Removing cheap points cheaply saves cheaply.

---

## 5. Recommendation — NO-GO, with a better lever

**NO-GO on S1 (reduced/octahedral grid).** Two independent, each-sufficient
blockers:

1. **SHTns cannot represent the grid** (§1). Capturing S1 requires replacing the
   transform with a libsharp/ducc0-class library — out of scope, and the single
   largest change one could make to aeros.
2. **Even with the grid, the payoff is ~1–5% (likely net-negative)** for aeros'
   84–97%-transform profile (§3), against a breaking refactor of the FV moisture
   core and pervasive changes to every physics loop and the IO layer (§2).

The reduced-grid idea is sound for gridpoint-dominated models (IFS at high
truncation) and for models that own a grid-flexible transform (SpeedyWeather). It
is a **category error for aeros**, whose cost is transform-bound at low truncation
through a rectangular-grid library.

### Cheaper alternatives that actually attack the 84–97%

The right target is the transform count and its per-call efficiency, not the
gridpoint count. In descending ROI:

- **S2a — Cut the number of transforms per step (highest ROI, low effort).**
  aeros currently issues **11 scalar-equivalent transforms per step**
  (`m0a_results.md` note 2), one field at a time (`aeros_tendency.f90:400–464`).
  Each transform removed is **~9% of the dominant cost**. SHTns' combined
  `SHqst_to_spat` can deliver ζ alongside (u,v) in one call and is documented
  "significantly faster" than separate transforms — the module already flags this
  as a TODO (`aeros_tendency.f90:76`, `m0a_results.md:120`). Realistic: 11→10 (and
  possibly fewer) ≈ **~9%+ off the step**, i.e. more than the *entire* optimistic
  reduced-grid win, for a few days of work and **zero** grid/physics/IO churn.

- **S2b — Batch levels into single SHTns calls (design.md §4's "all levels of a
  wavenumber in one DGEMM").** SHTns ships `shtns_set_many(cfg, howmany, spec_dist)`
  (`shtns.h:166–167`), documented as working on CPU, not just GPU. aeros presently
  loops transforms per level (`aeros_tendency.f90:409–411`, the `do k` loop). Level
  batching improves cache/SIMD reuse of the Legendre matrices and **improves with
  more levels** — directly helping the L12→L16–20 direction design.md §4.1 wants.
  Moderate effort; native to SHTns-Float64/OpenMP; needs benchmarking against the
  measured per-level path since the level loop is also aeros' parallel unit
  (`aeros_spectral.f90:87–113`) — so batching trades thread-parallelism granularity
  for per-call efficiency and must be measured, not assumed.

- **On-the-fly vs stored Legendre:** already handled — SHTns' tuner
  (`quick=.FALSE.`) picks the faster of `sht_gauss` (stored) vs `sht_gauss_fly`
  (on-the-fly) per config, and aeros uses the tuned path (`m0a_results.md:135`).
  No lever left here.

- **Float32 on CPU: confirmed dead** (`m0a_results.md §5`: dp is ~17% *faster* than
  sp through SHTns). Not a lever until/unless the transform moves to GPU.

**Recommended next speed action:** implement **S2a** (combine `SHqst_to_spat`,
reduce the 11-transform count) and **benchmark S2b** (level batching via
`shtns_set_many`). Together they target the 84–97% the reduced grid cannot touch,
with none of its blast radius.

---

## 6. Sources

- SHTns 3.7.5 source, shipped in-tree: `fesm-utils/SHTns/shtns-3.7.5/shtns.h`
  (grid API, `shtns_info` struct, `shtns_type` enum, `PHI_RAD` macro,
  `shtns_set_many`), `README.md` (grid families).
- aeros internal: `src/aeros_spectral.f90` (SHTns wrapper, grid-size rule, the
  reduced-grid comment at :411–417), `src/aeros_grid.f90`, `src/aeros_defs.f90`
  (`aeros_grid_class`), `src/dynamics/aeros_moisture.f90` (FV tracer transport),
  `src/dynamics/aeros_tendency.f90` (per-step transform calls), `src/aeros_io.f90`
  (rectangular netCDF), `docs/m0a_results.md` (transforms 84–97% of step; 11
  transforms/step; dp>sp), `docs/design.md` §2–4, `docs/refs/speedy_comparison.md`
  (S1).
- Reinecke & Seljebotn (2013), *Libsharp — spherical harmonic transforms
  revisited*, A&A 554 A112 — the reduced-grid (variable-`nphi`-per-ring) library
  lineage SpeedyWeather builds on; **not** SHTns.
- Schaeffer (2013), *Efficient spherical harmonic transforms…*, arXiv:1202.6522 —
  SHTns' own paper (regular Gauss/equispaced grids).
