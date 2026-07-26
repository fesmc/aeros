# aeros — design plan (draft for review)

A fast global atmospheric model for coupled climate–ice-sheet simulation on
10⁴–10⁶ yr timescales, with best performance over polar domains.

Status: **proposal**, rev. 3. This is a long project; the document is meant to be
picked up cold. Start at §0 for what has changed, §10 for what is still open.

**Settled so far:** primitive-equation spectral core (§3.2); OpenMP-only (§4.3);
real fine grids on both hemispheres, no elevation classes (§2.2); CLIMBER-X's
ocean adopted as-is (§6.1); no acceleration (§1.3); no PDD (§6); throughput
target 5,000 yr/day *fully coupled* (§3.5).

**Still open:** truncation T31+correction vs bare T42 — decided by the OpenMP
scaling measurement at M0a (§10.3); radiation scheme (§10.4); library-first
(§10.8).

**The two ideas that carry the most risk and the most upside** are §3.7
(multi-resolution error correction) and §3.8 (two-term correction). M2b tests
them and is the pivot of the whole plan.

---

## 0. Revision note — what changed and why

### rev. 3 (2026-07-22)

Throughput target re-framed as **fully coupled** (AR: CLIMBER-X's 10 kyr/day
includes ocean, ice, land — and that is aeros' target configuration too). This
tightened the budget enough to create the central resolution bind (§3.6), which
§3.7 exists to resolve. Added §3.7 (multi-resolution error correction, AR's
proposal), §3.8 (constant bias term alongside it, AR), §4.3 (OpenMP-only), §6.1
(adopt CLIMBER-X's ocean), and settled the elevation-class question against
elevation classes (§2.2). Added milestone M2b as the pivotal experiment.

### rev. 2

Rev. 1 recommended a 3-layer quasi-geostrophic core (ECBilt-class) on the
grounds that primitive equations were ~140× too expensive. **That was wrong on
both the cost and the science, and the recommendation is reversed.**

**Cost error.** Rev. 1 anchored PE cost on PLASIM-GENIE's ~70 yr/node-day. That
is a bad anchor: it includes a 3-D GOLDSTEIN ocean, PlaSim is 98% of the
runtime, and PlaSim calls full SW+LW radiation *every timestep* because `nrad`
is an on/off switch rather than a frequency. The correct anchor is
**SpeedyWeather.jl at T31L8: ~2,300 SYPD on one Apple M4 core** — essentially
identical to ECBilt's 2,500 SYPD. Primitive equations at comparable truncation
are not meaningfully more expensive than quasi-geostrophy.

**Arithmetic error.** Rev. 1's budget table was wrong by ~470× through a units
mistake (86400/4000 = 21.6 s per model year, not 46 ms; and 2500 SYPD = 34.5
core-s/yr, not 0.035). The corrected table is §3.6. The conclusion it produced —
"the fine-grid surface scheme is the binding constraint" — was an artefact of
that error and is false. Fine-grid surface physics is ~1% of budget.

**Science.** Löfverström & Liakka (2018) is decisive and was not available in
rev. 1: T21 fails to build *either* NH ice sheet. ECBilt is T21.

---

## 1. Niche

Ganopolski et al. (2010) give the mission statement: a 50 km 3-D
thermomechanical ice-sheet model is **~1% of coupled-model CPU**. The entire
climate-matrix / glacial-index / emulator industry exists only because the
*atmosphere* is the bottleneck.

Those methods' documented failure modes become the requirements list. Each is
something aeros should **produce**, not interpolate:

| missing in matrix/index methods | evidence |
|---|---|
| stationary-wave response to ice topography | Abe-Ouchi 2007 names it as the missing term; 2013 patches it with *one* fixed anomaly pattern ramped linearly in NH ice area |
| orographic precipitation | "inability to smoothly track orographic precipitation" — Pollard et al. 2013's own assessment |
| insolation as a genuine axis | in ANICE it enters only via absorbed insolation; all snapshots share one orbit |
| millennial variability | Berends smooths proxies with a 4 kyr running mean because D–O events "are not present in our model forcing or climate reference runs" |
| nonlinear CO₂ × orbit × ice interaction | measured at up to 15–20% by Abe-Ouchi 2007 — and then set to zero |

Two sensitivity results define how much this matters. In CLIMBER-2, disabling
the **elevation–desertification** effect gives **+40% LGM ice volume**;
disabling the **slope-precipitation** effect makes NH ice **vanish entirely** at
100–80 ka and reach only half of LGM volume.

### 1.1 The scientific differentiator, stated precisely

aeros' competition is CLIMBER-X/SESAM at ~10,000 SYPD. SESAM is more capable
than a caricature EMBM: it has **prognostic EKE** with Eady-growth-rate
production and √EKE-scaled diffusivity, **Charney–Eliassen topographic
stationary waves**, and **in-atmosphere orographic precipitation**
(`C_slope = c√K |∇z_s| ρ₀ q_a`). Any honest justification must be narrower than
"SESAM lacks dynamics." It is:

1. **SESAM's stationary-wave decomposition is additive** (`p_sl = p_sl,T +
   p_sl,O`). Liakka (2012, Tellus A 64, 11088) finds topographic stationary
   waves *in isolation* have "essentially no impact on the equilibrium features
   of the ice sheet"; the effect emerges only through a **strong non-linear
   interaction between thermal and topographical forcing**. An additive scheme
   has no cross-term.
2. **Charney–Eliassen is a 1-D barotropic channel with no meridional
   propagation.** Löfverström et al. (2014, App. A): real stationary waves
   "tend to follow great circles rather than latitude circles." Liakka et al.
   (2016): linear-model results are "a qualitative first-order estimate."
3. **A diffusive closure structurally cannot carry the flux.** Liakka &
   Löfverström (2018, CP 14, 887): for a high Laurentide the **stationary-wave**
   contribution to Arctic energy-flux convergence *dominates over transient
   eddies*. Eddy diffusivity — however scaled — represents the transient flux.
4. **Linearity is regime-dependent and sets ice-sheet geometry.** Liakka et al.
   (2012): a linear atmosphere yields an ice sheet with mass centre shifted far
   eastward; the nonlinear atmosphere yields something "virtually identical to
   uncoupled simulations." A linear parameterization cannot know which regime it
   is in.

Magnitude check: switching CLIMBER-2's azonal SLP off changes remote ice volume
by **30% to 270–400%** (Beghin et al. 2014).

**Honest counterweight.** Abe-Ouchi et al. (2007): the stationary-wave effect is
"secondary compared to the albedo and lapse rate effects" for summer temperature
*over* the ice. Stationary waves are second-order for the SMB **mean** and
first-order for ice-sheet **geometry**. Claim the latter, not the former.

### 1.2 Two arguments to retire

- **Interannual variability.** Changing *interannual* variability moves Greenland
  SMB by **<5%**; using a daily climatology instead of resolved *intra-annual*
  variability moves it by **40%** (TC 18, 4831, 2024). What is needed is a
  correct seasonal cycle and synoptic-scale intra-annual distribution — not
  realistic ENSO/NAO. CLIMBER-X's "no interannual internal variability" caveat is
  therefore a weaker indictment than it appears.
- **Mean-state skill over Greenland.** A coarse global model's Greenland is
  *worse* than REMBO's. Do not claim otherwise; that is what the polar module
  (§2) is for.

### 1.3 Non-goal: acceleration

Ruled out as a cost strategy. Ganopolski et al. (2010), N = 5/10/20 vs a
synchronous 120 kyr baseline: NH ice volume survives to N = 20, but millennial
variability is **"considerably suppressed for a two-fold acceleration"** and
deep-ocean temperature is degraded at every factor tested. Lunt et al. (2006):
1–2 °C Southern Ocean / Antarctic surface error at ×10, with an explicit
recommendation not to accelerate Southern-Ocean-dependent processes. Speed must
be structural.

---

## 2. Architecture: global core + embedded polar modules

Cost splits into two categories that scale differently — *nonlocal* work
(transport, spectral transforms, CFL-limited steps) and *pointwise* work
(radiation, condensation, surface energy balance, snowpack). The second is
embarrassingly parallel and, per §3.6, **cheap** — around 1% of budget.

The reason for a two-tier design is therefore **not** cost. It is that **no
affordable global resolution resolves what matters at the ice margins**:

- East Greenland's precipitation maximum is a ~50 km orographic feature.
- Last-glacial inception occurred in the mountains of Baffin Island; GCMs miss it
  "likely because of their coarse horizontal resolution that smooths topography"
  (Birch, Cronin & Tziperman 2017). A regional WRF configuration was needed.
- The same failure signature recurs across model classes, always the same shape:
  CLIMBER-2/SEMI puts Greenland's precipitation maximum on the *south-west*
  flank instead of the south-east; REMBO gives "too much precipitation on the
  southwest coast and not enough on the southern tip"; iLOVECLIM gives "too
  large a thickness in its central part owing to the overestimation of
  precipitation."

| | **global core** (T42, ~8k columns) | **polar modules** (~20 km, ~40k points) |
|---|---|---|
| transport, dynamics, stationary waves | prognostic | — |
| full radiative transfer | yes, 1–3 h | no — linearised in z, albedo, cloud |
| condensation, orographic precip | coarse closure | **yes**, from resolved slopes |
| surface energy balance, snowpack, SMB | no | **yes** |

The global model's job over Greenland is to get the **synoptic moisture delivery
direction and the 500 hPa flow** right — not the precipitation.

### 2.1 Two-way coupling

Conservative flux aggregation: the fine grid computes column tendencies;
conservative area-remap returns them as the coarse tendency. Conservative by
construction if the remap is. This is a deliberate upgrade over REMBO, whose
coupling is explicitly uni-directional — the ice-sheet→climate feedback (albedo,
surface temperature, topography) is exactly what one-way downscaling cannot
carry.

Second upgrade over REMBO: REMBO's diffusivities `D_T(φ, z_s)`, `D_Q(φ)` are
prescribed functions fitted to reproduce Greenland's *present-day* seasonal
cycle and cannot respond to a changed circulation. With a real global model
above, supply resolved winds directly.

### 2.2 Why a real fine grid, not elevation classes

Sellevold et al. (2019), CESM 10-EC vs RACMO2.3, Greenland:

| gradient | elev. classes | RACMO | error |
|---|---|---|---|
| net solar radiation | −3.5 W m⁻² km⁻¹ | −19.6 | **~6× too weak** |
| albedo | 0.019 km⁻¹ | 0.081 | **~4× too weak** |
| melt | −425 mm yr⁻¹ km⁻¹ | −717 | ~41% too weak |
| net SMB | 439 mm yr⁻¹ km⁻¹ | 369 | *close — by compensating errors* |

In CLM5, precipitation is not downscaled at all (only rain/snow partition), nor
shortwave, nor wind. An SMB gradient that looks right while every underlying
energy term is wrong by 4–6× is exactly the failure you cannot afford in a
deglaciation where the ablation zone migrates vertically.

**Decision (AR question, 2026-07-22): do both hemispheres on real fine grids;
do not use elevation classes for the ice sheets.** The reasoning is now settled
by §3.6 rather than by preference:

- **Cost does not force the compromise.** Fully interactive NH + SH ice sheets,
  with NH domains sized for *LGM* extent (Laurentide ~16×10⁶ km², Eurasia
  ~6×10⁶ km², plus Greenland and Antarctica ≈ **~96,000 points at 20 km**),
  costs ~12 core-s/yr with a BESSI-class multilayer snowpack — **~4% of the
  coupled budget.** Elevation classes would save ~3% of budget in exchange for
  the errors tabulated above.
- **Directionality is the signal.** Elevation classes are a *scalar* height
  binning: they cannot represent upslope/downslope asymmetry, precipitation
  shadows, or katabatic drainage direction. For Greenland and Antarctica that
  structure *is* what distinguishes a good SMB field from a bad one — and it is
  exactly what every coarse model gets wrong in the same characteristic way
  (§2, precipitation maximum on the wrong flank).
- **The remapping machinery already exists** — `fesm-utils` `coupler.f90` +
  SCRIP conservative maps + `coords`. Two-way conservative aggregation (§2.1) is
  the mechanism, and it preserves directionality by construction.

Elevation classes remain useful in two narrower roles: (a) a **cheap third tier**
for non-nested mountainous regions (Himalaya, Andes) where a nest is not
warranted; (b) a **fallback first implementation** if the two-way conservative
remap proves harder than expected, letting M3 proceed while M4 is debugged.

**Nest sizing note:** because the ice sheets grow, each nest domain must be sized
for *maximum* (LGM) extent, not present-day. A nest that the ice sheet can grow
out of is a silent failure mode.

### 2.3 Fine-grid geometry

Polar-stereographic patches matching the ice-sheet model grids directly — no
pole singularity, no double remap through a lat-lon intermediate. `fesm-utils`
already provides `coords`, `coupler.f90` (`coupler_add_grid`, `coupler_prime`,
`remap`) and SCRIP conservative maps.

---

## 3. The core: resolved primitive equations

### 3.1 The decisive constraint

⭐ **Löfverström & Liakka (2018, TC 12, 1499)** — four identical CAM LGM runs at
T85/T42/T31/T21, each used to force an ice-sheet model:

| resolution | result |
|---|---|
| **T85** (1.4°) | reproduces LGM ice sheets to high accuracy |
| **T42** (2.8°) | **fails to build the Eurasian ice sheet**; North America OK |
| **T31** (3.8°) | same failure |
| **T21** (5.6°) | **fails both**; North America splits into two disconnected domes |

Mechanisms, both quantified: weakened stationary planetary waves reduce cold-air
advection; and topographic smoothing (500–1000 m lost over Eurasia) drives cloud
fraction from ~50% at T85 to **nearly 100% at T21**, adding **10–30 W m⁻²**
downwelling LW. Their stated floor: "around T31, and possibly somewhat higher
(nominal T42 or even T85)."

The same paper notes the common denominator among studies showing this failure
pattern — they all used CLIMBER-2 forcing — and raises "a fundamental problem
with low-resolution climate models that transcends model complexity."

**This rules out T21.** ECBilt is T21, which removes quasi-geostrophy as a
serious option: the one place it is cheap enough to be attractive is the one
place it is demonstrably unable to build NH ice sheets.

### 3.2 Recommendation

**Semi-implicit leapfrog Eulerian spectral primitive equations**, hydrostatic,
triangular truncation, on a reduced/octahedral Gaussian grid. **T42, L16–20.**

Triangular truncation is **isotropic on the sphere** — no polar filter, no polar
CFL penalty, effective resolution at 85°N identical to the equator. For a
polar-priority model this is structural, not a tuning choice, and it is the
strongest technical argument for a spectral core. Contrast lat-lon grid-point
(polar Fourier filter damps exactly the region of interest) and cubed-sphere /
icosahedral (quasi-uniform, but worse per-DOF accuracy at coarse resolution;
ICON takes five dynamics sub-steps per physics step to handle acoustic modes we
do not need at 300 km).

### 3.3 Rejected alternatives

- **Quasi-geostrophic (ECBilt-class)** — §3.1. Also: ECBilt's radiation is an
  empirical Green's-function *anomaly* model relative to an NCEP reference
  climatology, which cannot be trusted at 180 ppm; and the authors state "the
  quasi-geostrophic structure of the model limits its ability to simulate
  equatorial variability."
- **Statistical-dynamical (SESAM-class)** — no cost advantage worth the loss
  (§1.1), and CLIMBER-X already occupies that niche at 10,000 SYPD. Building a
  second one is not a project.
- **Semi-Lagrangian advection** — the famous 6× IFS speedup exists at
  T_L511–T_L1279 where the Eulerian advective CFL is brutal. At T42 a
  semi-implicit Eulerian core already runs at 30 min. The cost is
  non-conservation of mass and tracers requiring a mass fixer — a drift risk
  over 10⁵ yr that buys nothing here.
- **Stretched-grid (Schmidt-transform) spectral** — refinement over one pole is
  "compensated by a coarser resolution on the opposite side of the sphere."
  Stretching ×3 on Greenland gives Antarctica ~T14. Fatal, since both poles
  matter. It also destroys the isotropy that motivated the spectral choice.
- **Variable-resolution grid-point (VR-CESM)** — mature and effective (SMB
  "mainly improves with refinement"), and it preserves two-way interaction. But
  it requires a spectral-element or unstructured core, forfeiting §3.2's polar
  isotropy, and is not cheap.

### 3.4 Throughput anchors

| model | throughput | config |
|---|---|---|
| **SpeedyWeather.jl** | **~2,300 SYPD, 1 M4 core** (500–600 on older Intel) | T31L8, Float32 |
| SPEEDY (Fortran) | 125–240 SYPD/core | T30L8 |
| PlaSim / ExoPlaSim | ~1,400–2,900 SYPD | T21L10 |
| ECBilt (iLOVECLIM) | 2,500 SYPD, 1 core | T21L3, coupled |
| CLIMBER-X | 10,000 SYPD, 16 cores | 5°, statistical-dynamical |
| FAMOUS | ~120 SYPD | 7.5°×5°L11, full HadCM3 physics |
| CCSM3 T31 | ~35 SYPD | — |

The 4–5× spread for *identical* SpeedyWeather code between M4 and older Intel
shows **single precision, SIMD and cache behaviour dominate**. A Fortran rewrite
has no intrinsic advantage; it must earn speed the same way.

### 3.5 Target (set by AR, 2026-07-22)

- **1,000 yr/day = floor**; **3,000–5,000 yr/day = target**, coarse+fine,
  coupled. 10,000 (CLIMBER-X) accepted as probably out of reach with prognostic
  dynamics.

### 3.6 Cost budget (corrected twice — see §0)

**Important framing correction (AR, 2026-07-22): CLIMBER-X's 10,000 yr/day is
the FULLY COUPLED model** — atmosphere, ocean, sea ice, land, ice sheets. That
is also the target configuration for aeros, so the budget below is for
*everything*, not for the atmosphere alone.

At **5,000 yr/day on a 16-core node**: 86400/5000 = **17.3 s wallclock per model
year**, × 16 = **~276 core-seconds per model year, total, all components.**

| component | core-s / model yr | share of 276 |
|---|---|---|
| ECBilt-class QG T21L3 (coupled) | 34.5 | 13% |
| SpeedyWeather PE T31L8 | 37.6 | 14% |
| **PE T31L16** (est.) | **~57** | **~21%** |
| **PE T42L19** (est.) | **227** | **~82%** ⚠ |
| PE T85L19 (est.) | ~3,600 | 1300% ✗ |
| dEBM SMB, 41k points | 0.57 | 0.2% |
| BESSI snowpack, 41k points | 5.0 | 1.8% |
| BESSI snowpack, ~96k points (NH+SH, LGM extent) | ~12 | 4.3% |

T31L16 from T42L19 × (31/42)³ × (16/19) × (31/42, fewer steps) ≈ 0.25.

**Three conclusions:**

1. **PE at T31L8 costs the same as QG at T21L3** (37.6 vs 34.5). There is no
   cost argument for quasi-geostrophy. Unchanged from rev. 2 and still the
   single most important number in the plan.
2. **Surface physics is cheap even at full ambition.** A BESSI-class multilayer
   snowpack on *both* hemispheres at 20 km, with NH domains sized for LGM
   extent (~96k points), is **~4% of the coupled budget**. Do not compromise the
   surface scheme to save cost. See §2.2.
3. **T42L19 does not fit a 5,000 yr/day coupled target.** At 82% of budget
   before ocean, ice and land are counted — and assuming near-perfect 16-core
   scaling, which spectral transforms at ~8,200 columns will not deliver —
   realistic coupled throughput is **~2,500–3,000 yr/day**. T31L16 fits
   comfortably at ~21%.

**This creates the central bind of the whole design:** the resolution affordable
in a coupled configuration (T31) is precisely the resolution Löfverström &
Liakka show *fails to build the Eurasian ice sheet* (§3.1). §3.7 is the proposed
way out.

**The budget scales with *effective* threads, so state it that way.** At
5,000 yr/day the wallclock budget is 17.3 s per model year; the core-second
budget is `17.3 × N_eff`:

| N_eff | budget (core-s/yr) | T31L16 (57) | T42L19 (227) |
|---|---|---|---|
| 8 | 138 | 41% | 165% ✗ |
| 16 | 276 | 21% | 82% ⚠ |
| **32** | **553** | **10%** | **41%** ✓ |

**This softens the §3.6 bind considerably.** At 32 *effective* threads T42L19
fits at 41%, leaving ~59% for ocean, ice and land — so bare T42 becomes viable
without §3.7 if scaling is good. The bind is real only if effective parallelism
stalls below ~24. Which makes §4.3 the decisive engineering question, not a
detail.

For long ensembles, prefer **ensemble-parallel** (one member per run) over
domain-parallel — see §4.3.

### 3.7 Multi-resolution error correction (proposed by AR, 2026-07-22)

**The idea.** Run the transient model at an affordable truncation (T31), and
periodically run a short high-resolution (T85) integration to diagnose the
*resolution error*, applying that as a correction until the next refresh. If it
works, it resolves the §3.6 bind: T85-like behaviour at ~T31 cost.

**Formulation — diagnose the model-error operator, do not nudge toward a
state.** This distinction is the whole design.

At refresh *k*:
1. Freeze boundary conditions `B(t_k)` — topography, ice mask/extent, SST, sea
   ice, CO₂, orbit.
2. Integrate **T85** for a short window (months to ~2 yr) under `B(t_k)`.
3. Integrate **T31** for the *same* window under *identical* `B(t_k)`.
4. Diagnose `ΔF_k = (terms)_T85 − (terms)_T31`, coarse-grained to T31.
5. Apply `ΔF_k` as a **fixed additive forcing** in the transient run for
   `t ∈ [t_k, t_{k+1})`.

The transient model continues to respond freely to evolving forcing and to all
its own feedbacks; only the *resolution error* is held fixed between refreshes.

**Why this is not a bias correction.** The documented killer of anomaly/matrix
methods is a correction diagnosed in one climate state and applied in another —
quantified at **2.3 m SLE** of spread for MIS-11c Greenland purely from the
choice of reference climatology (§6). Here the correction is (a) re-diagnosed as
the climate state evolves, and (b) generated by *the same model with identical
physics*, differing only in truncation. Neither property holds for a
conventional bias correction.

**Timescale separation is excellent.** Stationary waves adjust in ~10–30 days;
ice topography evolves on 10²–10³ yr — roughly four orders of magnitude. This is
squarely the regime where heterogeneous-multiscale / equation-free coarse
projective integration is well posed, rather than a hopeful approximation.

**Cost.** T85L19 ≈ 16× T42L19 ≈ 63× T31L16. Two years of T85 per 500 yr of
transient adds **~25%** to the T31 atmosphere cost, i.e. ~5% of the coupled
budget. Refreshing every 1,000 yr halves that. Affordable enough that refresh
frequency can be tuned by accuracy rather than by cost.

**Refresh trigger: state-based, not purely time-based.** Trigger on the slow
variable that actually drives the correction — e.g. RMS ice-surface-elevation
change exceeding a threshold (~100–200 m), or ice-area change exceeding a few
percent — with a time-based backstop. Ice topography is the right trigger
because the stationary-wave response is what degrades at T31.

**Correct selected terms, not everything.** The evidence identifies what
degrades at coarse truncation: (a) stationary-wave forcing / resolved orographic
drag, (b) cloud fraction over smoothed ice topography (50% → ~100% from T85 to
T21, worth 10–30 W m⁻² of downwelling LW), (c) eddy heat/momentum flux
convergence. A blanket tendency correction on all prognostic fields is far more
likely to destabilise and much harder to diagnose.

**Risks:**

1. **Untested bridge.** Löfverström & Liakka show T85 works and T31 does not.
   Nobody has shown T31+Δ reproduces T85 ice sheets. This is a genuine research
   risk and must be an explicit early milestone (M2b), not an assumption. If it
   fails, the fallback is T42 at ~2,500–3,000 yr/day — above the floor, below
   target.
2. **Conservation.** `ΔF` must be constructed in flux form and sum to zero
   globally (or be explicitly accounted), or it injects spurious energy/water
   over 10⁵ yr. Verify to machine precision.
3. **Physics/resolution confound.** Parameterizations are tuned per resolution,
   so `ΔF` conflates genuine dynamical resolution error with parameterization
   mis-tuning. Arguably acceptable — the goal is T31+Δ ≈ T85 regardless of
   source — but `ΔF` cannot then be interpreted as physics.
4. **Feedback suppression.** Holding `ΔF` fixed suppresses the resolution
   error's *own* response to evolution between refreshes. Mitigated by the
   timescale separation above and by state-based triggering, but it argues
   against long refresh intervals during rapid ice change (deglaciation).
5. **Is T85 enough?** L&L validate T85 for LGM. It is the reference here, so any
   T85 error propagates into `ΔF` unexamined.

**Prior art to check before claiming novelty** (recalled, *unverified* — the
web-search budget was exhausted): heterogeneous multiscale methods /
equation-free coarse projective integration (Kevrekidis and co.); "nudging
tendency" model-error correction in the CAPT tradition; the machine-learning
corrected coarse-GCM line (Watt-Meyer et al. and related). Superparameterization
is the spatial analogue of the same idea. **Verify these before writing this up
as new.**

### 3.8 Two-term correction: resolution error + irreducible bias

**Decision (AR, 2026-07-22): reserve a second, constant, present-day bias term
alongside the state-dependent resolution correction.** Even at T85 the model has
real biases against observations, and for ice-sheet work the absolute position
of the ablation zone matters, not just the anomaly.

The total correction decomposes into two terms with **different provenance,
different time-dependence, and different risk profiles**, and they must be kept
separately diagnosable:

```
ΔF_total(t) = ΔF_bias           +  ΔF_res(t)
              ^ constant           ^ state-dependent, refreshed (§3.7)
              T85 − observations   T85 − T31, under identical BCs
```

| | `ΔF_bias` | `ΔF_res(t)` |
|---|---|---|
| source | T85 vs ERA5 / MAR / RACMO | T85 vs T31, same physics |
| time-dependence | **constant** | re-diagnosed per refresh |
| what it fixes | irreducible model error | truncation error |
| risk | **high** — see below | low (§3.7) |

**The honest tension.** `ΔF_bias` is structurally the same object this plan
spends §6 arguing against: a correction diagnosed in the present-day climate and
applied across glacial states. That practice produced **2.3 m SLE of spread for
MIS-11c Greenland from the choice of reference dataset alone** (MAR ~3.2 vs
RACMO ~5.5), and it is the documented weak point of CLIMBER-X's fixed modern
`prc_bias_i`/`t2m_bias_i`. Including `ΔF_bias` is a deliberate, eyes-open
reintroduction of that risk because the alternative — accepting a misplaced
present-day ablation zone — is also unacceptable.

**Five mitigations, all of which should be implemented, not just intended:**

1. **Define `ΔF_bias` as the residual *after* `ΔF_res`, at the highest available
   resolution.** It is `T85 − obs`, never `T31 − obs`. This keeps it as small as
   the model allows and prevents it silently absorbing truncation error that
   `ΔF_res(t)` is supposed to carry — the failure mode that makes conventional
   bias corrections state-dependent in disguise.
2. **Require it to be small, and treat growth as a bug.** Track
   `‖ΔF_bias‖ / ‖ΔF_res‖`. If the constant term dominates, the model is being
   held to reality by force rather than physics, and the paleo results are not
   trustworthy.
3. **Make it switchable and always report both.** Every production run should
   have a `ΔF_bias = 0` twin. If conclusions differ qualitatively, say so.
4. **Bound the reference-dataset sensitivity explicitly.** Diagnose `ΔF_bias`
   against ERA5, MAR *and* RACMO and report the spread — this is precisely the
   MIS-11c test, run on our own model instead of inherited. If the spread is
   large, that is the honest uncertainty on any paleo result.
5. **Prefer terms whose bias is plausibly structural** (e.g. systematically
   deficient orographic drag, surface roughness over ice) over terms likely to
   be state-dependent (cloud fraction, precipitation partition). A constant
   correction to a state-dependent error is the exact mechanism that fails.

**Do not apply `ΔF_bias` to the ice sheet's SMB directly.** Correct the
*atmospheric* fields that feed the surface scheme, so the SMB remains a physical
consequence of a corrected climate rather than a tuned field. Otherwise the
elevation feedback runs on a corrected quantity and the geometry error described
in §6 (co-location of prescribed anomalies with a modelled ablation zone)
returns.

---

## 4. Numerics

- **Reduced / octahedral Gaussian grid** — reduced Gaussian alone saves "in
  excess of one-third" the gridpoints (Hortal & Simmons 1991); octahedral saves
  a further ~22% of compute at ECMWF. ~30% off gridpoint physics for negligible
  accuracy loss, and it disproportionately removes *redundant polar columns*
  without touching polar spectral resolution.
- **Float32 throughout the core.**
- **Legendre transforms: on-the-fly recurrences, all levels of a wavenumber
  batched into one DGEMM, explicit SIMD.** The O(N³) scare is irrelevant here —
  SHTns is already cubic at ℓ=63 and still beats lower-complexity algorithms up
  to ℓ=1023. Batching *improves* with more levels, partly offsetting L8→L19.
  Estimated 25–40% of runtime at T42L19. Do not implement a fast/butterfly
  Legendre transform; ECMWF only needed it at T3999+.
- **Double Fourier series** (O(N² log N)) is the one real alternative, but no
  published crossover wavenumber exists and it forfeits polar isotropy. Not
  worth the novelty risk at T42.
- **Δt ≈ 30 min** (linear from SPEEDY/SpeedyWeather's 40 min at T31).
- **Implicit or flux-limited vertical advection.** In σ-coordinates σ̇ picks up
  `u·∇h_s`; over ice-sheet margins with |∇h| ~ 0.02 and u ~ 20 m/s this gives
  effective w ~ 0.4 m/s and Δt < ~2500 s for Δz = 1 km. **Vertical CFL over ice
  margins, not horizontal CFL, may cap the timestep** once levels are refined
  near the surface. Argues for hybrid σ–p and smoothed orography.

### 4.1 Vertical levels — the highest-value decision

**Recommendation: L16–20, with 3–4 levels below 850 hPa.**

- Roeckner et al. (2006), ECHAM5 across T21L19→T159L31: "**L19 vertical
  resolution is adequate for T31 and T42**"; at L19 there is no convergence
  benefit above T42.
- **Levels are nearly free.** ECHAM6 uses **10 min at both L47 and L95** — level
  count does not shorten the timestep in this regime. L19 costs ~2.4× the work
  of L8 with no timestep penalty.
- **The polar boundary layer demands it.** Lapse-rate feedback is the single
  largest driver of polar amplification; "the climatological near-surface
  inversion over the Arctic strongly suppresses vertical mixing." Byrkjedal et
  al. (2008): 31 layers "cannot be expected to realistically resolve the Arctic
  stable boundary layer"; 90 layers substantially improved agreement with
  ERA-40. A polar winter SBL is 50–200 m deep with 10–25 K inversions — a single
  ~400 m SPEEDY-style bulk layer physically cannot hold one, making near-surface
  polar temperature a diagnostic of the flux scheme rather than of the dynamics.
  **If polar performance is the point of the model, put 3–4 levels below 850 hPa
  before adding anything aloft.**
- **Fix `nlev` before tuning anything.** Wan et al. (2008): higher vertical
  resolution shifts the westerly jets equatorward. Storm-track position — hence
  ice-sheet precipitation — is a function of level count.
- Low-top is acceptable: CMIP5 low-top and high-top have similar *mean*
  stratospheric climate (low-top has weak variability, ~half the observed SSW
  frequency). Give up on the QBO.

### 4.3 Parallelism: OpenMP only (AR, 2026-07-22)

**Decision: OpenMP, shared memory, no MPI.** Target good scaling to ~32 threads;
do not chase 128.

**The key point is that MPI would not help.** The limit is *work per thread*,
not the programming model. T31L16 has ~4,600 columns and ~500 spectral
coefficients per level; T42L19 has ~8,200 columns. Spread over 128 threads that
is tens of columns each, with a global transpose every timestep. Distributed
memory would add communication cost to a problem that is already
communication-bound at that thread count. CLIMBER-X's use of **max 32 threads
for scalability reasons** is the relevant empirical datapoint, and aeros' core
is of comparable size. OpenMP is therefore not a compromise here — it is the
correct choice, and it keeps the code far simpler.

**Where the parallelism actually comes from:**

- **Column physics** (radiation, condensation, surface, snowpack) — trivially
  parallel over columns, scales to whatever thread count exists. This is also
  the part that grows with the fine grids (§2.2), so it scales *better* as
  ambition increases.
- **Spectral transforms** — parallelise over (level, wavenumber) pairs and lean
  on threaded BLAS for the batched DGEMM (§4). Batching all levels of a
  wavenumber into one call is what makes this both fast and threadable, and it
  improves with more levels — another argument for L16–19 over L8.
- **The §3.7 T85 diagnostic runs** — ~8× more columns than T31, so they scale to
  higher thread counts naturally. Use the full node there.

**Use the remaining threads for ensembles, not for the core.** On a 128-thread
node, four 32-thread ensemble members will always beat one 128-thread run. For
paleo work — parameter ensembles, orbital sensitivity, `ΔF_bias = 0` twins
(§3.8) — this is the right shape anyway.

**Consequence for the design:** effective thread count is the single number that
decides T31+§3.7 vs bare T42 (§3.6 table). **M0a must measure the OpenMP scaling
curve to 8/16/32/64 threads before the truncation is fixed.**

### 4.2 Known risk: Gibbs oscillations at ice margins

A spectral core's worst pathology — spurious ripples in surface pressure and
negative specific humidity — bites hardest at steep orography, i.e. exactly the
Greenland and Antarctic margins that drive SMB. Mitigations: spectrally-filtered
orography, positive-definite gridpoint humidity, reduced Gaussian quadrature,
and margin-region diagnostics from day one. The coarse core sees deliberately
smoothed topography and is not asked to resolve the margin; steep-orography
effects live on the fine grid (§2). Budget real effort here — this is
under-documented in the literature and easy to underestimate.

---

## 5. Radiation

**Do not inherit SPEEDY's scheme.** Its LW band-2 transmissivity is
`τ = exp(−α_CO2·Δp′)` with **α_CO2 = 6.0 a constant**; there is no CO₂ mixing
ratio variable anywhere in the model. Its stratosphere is additionally nudged by
a prescribed zonally-symmetric correction, so CO₂ stratospheric cooling cannot
emerge. Disqualifying.

**Do not use grey radiation.** Against RRTMG it produces spurious −5 K/day peak
cooling (vs ~−2), a 322 hPa tropical radiative tropopause (vs 134), ~30% too
weak upper-tropospheric warming amplification, a subtropical jet at 26° (vs
36°), and 4–5× too weak Hadley weakening. Worst: clear-sky λ_LW rises to
~8 W/m²K at 330 K where the spectral value stays flat at ~2 — an artificial
stabilising feedback that does not exist. Mechanistically, CO₂ forcing is "a
swap of surface emission for stratospheric emission" within the 667 cm⁻¹ band;
this is irreducibly spectral and a one-band model has no analogue. **The error
is state-dependent and will not cancel between glacial and hothouse states.**

**Recommendation, in order:**

1. **A bespoke ecCKD gas-optics table at 16–32 g-points**, trained across
   180–2000 ppm from CKDMIP line-by-line data. ecCKD is a tool that generates
   *your own* table at a chosen g-point count. Real spectral structure, real
   stratospheric CO₂ cooling, validated forcing. Evidence that aggressive
   reduction is safe: ICON-A cut 480→240 g-points for ~2× speedup "without
   physically significant effects"; RRTMGP-NN reduced to 112 SW/128 LW with SW
   surface RMSE 0.78→0.80 W/m².
2. **CCM3-style broadband absorptivity–emissivity** with the
   Ramanathan/Kiehl closed-form CO₂ band absorptance
   `A(u,T,P) = 2A₀ ln{1 + u/√[4 + u(1+1/β)]}`. Physically grounded,
   concentration-aware, no lookup tables. Watch the O(nlev²) structure.
3. **SESAM's scheme** (`~/models/climber-x/src/atm/{lwr,swr}.f90`) — 2 SW bands,
   LW on 15 levels. **Confirmed by AR to work to ~2000 ppm** (unpublished),
   which removes the usual objection. It is designed around a
   column-integrated statistical-dynamical atmosphere, so adapting it to a
   resolved L19 PE column is not free — but it is local, validated, and cheap,
   and is the pragmatic fallback if (1) proves slow to build.

**Call frequency — the single most transferable finding.** Call full radiation
every **1–3 h**, but **rescale surface LW and SW fluxes every timestep and
gridpoint** from local skin temperature and albedo, conserving energy through
the flux profile. ECMWF measured this at **~2% of the radiation scheme, ≈0.2% of
total runtime**, and it removes the large surface-temperature errors and
diurnal-cycle lag that infrequent radiation otherwise causes — precisely where
skin temperature and albedo change fastest, which over ice sheets is exactly
where SMB is decided. RRTMG_LW ships an `idrv` option returning dF↑/dTₛ per
level for this purpose. SPEEDY already does a version of it.

For reference, IFS calls radiation hourly on a grid with ~10× fewer points than
the dynamics, and radiation is then only ~3.5% of runtime.

---

## 6. Physics and coupling

- **Moisture:** prognostic q with an explicit condensation scheme. CLIMBER-X has
  *none* — precipitation is a residence-time closure. Positive-definite
  gridpoint transport (§4.2).
- **Clouds:** diagnostic initially. Note both the T21 runaway-cloudiness result
  (§3.1) and CLIMBER-X's need for `0.1*new + 0.9*old` temporal relaxation on
  cloud fields — diagnostic cloud chains are delicate.
- **Surface:** sub-grid tiling per cell (ocean / sea ice / land / ice / lake),
  following CLIMBER-X's `nm=5`.
- **SMB: insolation-explicit, never PDD.** Bauer & Ganopolski (2017): *no* PDD
  parameter set works across a full glacial cycle — inception needs smaller melt
  factors than termination, America differs from Europe, and tuning to LGM
  volume makes Holocene deglaciation impossible. dEBM: the effective degree-day
  factor shifts >10% under mid-Holocene orbit alone (8.7→9.8 mm K⁻¹ d⁻¹ wet
  snow) vs ~0.3 for a century of warming. Robinson & Goelzer (2014): insolation
  is 20–50% of the peak-Eemian melt anomaly. **This is why `insol` is a
  first-class dependency.** Per §3.6 the budget is ample: use SEMIC / dEBM /
  BESSI-class, not a corner-cut.
- **Ice coupling: synchronous, annual.** Published deglacial transients couple at
  1–10 yr (iLOVECLIM 1, CESM2-CISM2 1, UKESM 1, MPI-ESM 10).
- **No anomaly/bias-correction forcing.** MIS-11c Greenland: identical CESM
  anomalies bias-corrected onto MAR give ~3.2 m SLE, onto RACMO ~5.5 m SLE —
  **2.3 m SLE of spread from the reference dataset alone.**

### 6.1 Ocean (AR, 2026-07-22: start from CLIMBER-X's)

**Decision: do not rebuild the ocean. Adopt CLIMBER-X's** and treat it as a
fixed component for the foreseeable milestones. Rationale: it is local,
validated, cheap, and already coupled to the ice-sheet and carbon components
aeros will eventually need. Rebuilding two components at once multiplies risk
for no scientific return on the stated goals, which are atmospheric.

Two consequences worth recording now:

- **Interface work is real but bounded.** CLIMBER-X's ocean expects fluxes from a
  column-integrated statistical-dynamical atmosphere. A resolved L16–19 PE
  column produces the same quantities, but the coupling layer (flux
  aggregation, timestep mismatch, the 1-day vs 30-min step) must be written
  deliberately rather than inherited.
- **The ocean will eventually become the limiting component, and that is
  acceptable for now.** A PE atmosphere generates synoptic variability that a
  coarse frictional-geostrophic ocean will simply average — fine for SMB, which
  is what aeros is for. It becomes limiting only for AMOC–ice interactions and
  millennial variability (D–O), where the atmosphere would then no longer be the
  weak link. **Flag, do not fix.** Revisit only after M5.

**Interim lower boundary (AR, M2.5d).** Until that coupling exists, the sea
surface is `aeros_ocean`: prescribed SST (the aquaplanet control) or a
well-mixed slab that closes the surface energy balance (`C dSST/dt = SW_net +
LW_down − σSST⁴ − SH − LH`). This is exactly the plug point CLIMBER-X's ocean
will attach at — the atmosphere hands the module the net surface flux and reads
back an SST, the same contract the slab already satisfies — so adopting the real
ocean later is a module swap, not a re-plumbing.

### 6.2 10⁵–10⁶ yr is not reachable transiently

At T42L19 (~227 core-s/yr), 10⁶ yr ≈ **7 core-years**. Even at target throughput
on one node that is ~8 months of wallclock for a *single* 10⁶ yr run. 10⁵ yr
(~25 days) is comfortable; 10⁶ yr is not, and ensembles at 10⁶ yr are out of
reach.

**Plan for a matrix or emulator layer on top of aeros, not instead of it.** The
two workable published strategies:

- **Snapshot matrix** (Berends et al. 2018): two snapshots, axes = scalar CO₂
  index + *spatially variable* ice index from local absorbed insolation;
  temperature and orography interpolate **linearly**, precipitation
  **logarithmically** (copy this detail). A 120 kyr run takes ~12 h on one core.
  Scales to 3.6 Myr with 11 snapshots.
- **GP emulator** (CLISEM v1.0): 100 training runs of 40 yr; inputs =
  eccentricity, obliquity, longitude of perihelion, CO₂, plus an ice descriptor
  calibrated on **≥12 predefined ice geometries** (with 8, the ice sheet cannot
  grow realistically); outputs monthly T and P directly on the ice grid.
  **Coupling interval 500–2000 yr** is the stated sweet spot. 3 Myr demonstrated.

Both published emulators fail in the same direction — PALEO-PGEM emulates LGM
cooling at 4.1 ± 0.2 °C vs 5.9 °C simulated, ~30% low at the edge of training.
Since ice growth is a positive feedback, under-predicting the extreme
under-predicts glaciation. Because aeros would be generating its *own* training
snapshots, the training envelope is at least known — a decisive advantage over
an ML atmosphere emulator, whose failure mode is to stay stable and plausible
while being quietly wrong about the forced response (ACE2's own abstract:
sensitivities to changing SST and CO₂ "are not entirely realistic"; NeuralGCM
"does not extrapolate to substantially different future climates"; all have
fixed land–sea mask, fixed orography, fixed insolation geometry).

---

## 7. Milestones

- **M0a — cost sizing and OpenMP scaling.** Before physics: benchmark a bare
  spectral transform + semi-implicit step at T31L16 and T42L19 on the target
  node, and **measure the OpenMP scaling curve at 8/16/32/64 threads**. Per
  §3.6 and §4.3 this single measurement decides the truncation. Cheap to run,
  and everything downstream depends on it.
- **M0 — scaffold.** configme + fesm-utils, `chion` as the template. Library
  target `aeros-static`, nml parameters, ncio + `variable_io` output, one driver.
  Requires adding `packages/aeros.toml` to the configme registry.
- **M1 — dry dynamical core.** Spectral PE via SHTns, reduced Gaussian grid,
  Float32. Validate against **Held–Suarez**. Benchmark T42 vs T63 here — make
  resolution a namelist change from day one (§9 risk 1).
- **M2 — moist + radiation, global.** Prognostic q, condensation, gas-optics
  radiation with per-timestep surface flux rescaling, diagnostic clouds, surface
  tiles. Validate vs ERA5: zonal-mean T, jet latitude, storm-track position,
  P−E, and **polar inversion strength**.
- **M2b — multi-resolution correction (§3.7). The critical experiment.** Run
  T85 and T31 under identical LGM boundary conditions; diagnose `ΔF`; verify
  (a) it conserves to machine precision, (b) T31+`ΔF` reproduces T85 stationary
  waves and cloud fraction over ice, (c) the correction is stable over ≥10 kyr.
  **Then force an ice-sheet model with all three and check T31+`ΔF` builds the
  Eurasian ice sheet where bare T31 does not.** This single test decides whether
  the 5,000 yr/day target is reachable or whether the design falls back to T42
  at ~2,500–3,000. Do it before M3.
- **M3 — polar modules, one-way.** Stereographic nests, linearised radiation,
  slope precipitation, SEMIC/dEBM/BESSI-class SMB. Validate against RACMO/MAR
  **gradients**, not just integrals (§2.2).
- **M4 — two-way coupling.** Conservative flux aggregation; verify global energy
  and water closure to machine precision with nests active.
- **M5 — paleo.** LGM stationary waves vs PMIP; **the Eurasian ice sheet is a
  first-order validation target, not an afterthought** (§9 risk 1); inception at
  115 ka; then coupled to Yelmo.

M0a, M1 and M5-Eurasia are the three places this can fail, and all are reachable
early.

---

## 8. Build / infrastructure

- `~/models/chion` is the cleanest scaffold to copy (recent, single dependency,
  library + drivers, `precision=sp|dp` variant handling).
- `config/` holds `Makefile` (template with `<COMPILER_CONFIGURATION>`),
  `common.mk` (dependency wiring), `Makefile_aeros.mk` (explicit per-object
  rules, no wildcards).
- Convention: `.mod` files and `libX.a` live in the **same** directory, so `-I`
  and `-L` point at the same path.
- fesm-utils provides `nml`, `ncio`, `variable_io`, `varslice`, `series`,
  `timestepping`, `timeout`, `timer`, `coupler`, and the `coords` library
  (projections, SCRIP remapping, kdtree, gaussian filter).
- External libs per variant (`serial`/`omp`):
  `fesm-utils/{fftw,lis,SHTns}/<lib>-{serial,omp}/{include,lib}`. **SHTns is the
  key one** — and there is **no Fortran wrapper in fesm-utils**, so
  `aeros_spectral.f90` over `shtns.f03` is ours to write.
- `insol` API: `insol_init` / `calc_insol_day` / `calc_insol_days` /
  `calc_insol_ave` (point, 1-D, 2-D), Laskar (2004) tables, `time_bp` in years
  before 1950.
- **`rembo`'s build is not configme-migrated** — do not copy it. `rembo1` is,
  and shows the pattern for exporting artifacts to a downstream orchestrator.

---

## 9. Risks

1. **The resolution bind (§3.6): the truncation affordable when coupled (T31)
   is the one shown to fail.** Löfverström & Liakka's T42 run failed Eurasia and
   T31 failed it too; their floor is "T31, possibly T42 or even T85." Meanwhile
   T42L19 is 82% of the coupled budget before ocean/ice/land.
   **Mitigation: §3.7, tested at M2b.** If §3.7 works, this risk closes. If it
   does not, the honest fallback is T42 at ~2,500–3,000 yr/day coupled — above
   the 1,000 floor, below the 5,000 target. Either way, make truncation a
   namelist change from day one and make the Eurasian ice sheet a first-order
   validation target, never an afterthought.
2. **Parallel efficiency — now the decisive engineering risk.** T31 has ~4,600
   columns, T42 ~8,200; spectral transforms need a global transpose every step.
   Per §3.6, effective thread count decides the truncation: N_eff ≥ ~24 makes
   bare T42 viable, below that T31 + §3.7 is required. MPI would not rescue this
   (§4.3) — the limit is work per thread. **Measured at M0a, before anything
   else is built.**
3. **Gibbs oscillations at ice margins** (§4.2) — the spectral core's worst
   pathology occurs exactly where SMB is decided, and it is under-documented.
4. **Duplicating CLIMBER-X at higher cost.** The differentiator is §1.1 and only
   §1.1. Design the validation around nonlinear thermal×orographic stationary
   waves and their contribution to poleward energy flux. Do not over-claim on
   variability (§1.2) or on Greenland mean state (§1.2).
5. **`ΔF_bias` quietly becoming the model** (§3.8). A constant present-day
   correction is the same object that produced 2.3 m SLE of MIS-11c spread. The
   mitigations in §3.8 — define it at T85 not T31, track its magnitude against
   `ΔF_res`, keep a `ΔF_bias = 0` twin for every production run, bound the
   reference-dataset spread — are not optional hygiene; they are what keeps the
   paleo results meaningful. If the constant term grows to dominate the
   state-dependent one, the model is being held to reality by force.

---

## 10. Open decisions

1. ~~Throughput budget~~ — **settled**: 1,000 floor / 3,000–5,000 target.
2. ~~Core formulation~~ — **recommendation: spectral primitive equations, T42,
   L16–20**, reversed from rev. 1. Needs sign-off.
3. **Transient truncation: T31 (with §3.7 correction) or T42 (bare)?** The
   central open question, and **§4.3's scaling measurement decides it**. At
   N_eff ≥ ~24 bare T42 fits the 5,000 yr/day target directly; below that, T31 +
   §3.7 is the route. Recommendation: **make truncation a namelist parameter,
   build §3.7 regardless** (it also buys ensemble headroom and closes risk 1),
   and fix the production truncation at M2b.
4. **Radiation: build an ecCKD table, or port SESAM's scheme?** (1) is better
   physics for a resolved column; (3) is local, validated to 2000 ppm and much
   faster to stand up. Leaning (3) for M2 and (1) later, since §3.7/§3.8 are the
   higher-risk novelty and deserve the development attention.
5. ~~Number and extent of nests~~ — **settled**: fully interactive NH + SH, real
   fine grids at ~20 km, nests sized for LGM extent, no elevation classes
   (§2.2). ~4% of coupled budget.
6. ~~Ocean~~ — **settled**: adopt CLIMBER-X's, do not rebuild (§6.1).
7. ~~Parallelism~~ — **settled**: OpenMP only, target ~32 threads, ensembles
   beyond that (§4.3).
8. **Library-first?** Callable from a yelmox-style driver from day one, or
   standalone first? Given §6.1 (coupling to CLIMBER-X's ocean) the answer is
   probably yes, but it is not yet decided.
9. **Present-day reference dataset(s) for `ΔF_bias`** (§3.8) — ERA5 for the
   atmosphere is obvious; the Greenland/Antarctic surface reference (MAR vs
   RACMO) is not, and the choice carries the 2.3 m SLE sensitivity. Plan to use
   both and report the spread.

---

## 11. Provenance

Built from parallel surveys of: local build infrastructure (`configme`,
`yelmo`, `chion`, `fesm-utils`, `insol`, `rembo`); local physics codebases
(`rembo`, `climber-x`); the coupling/downscaling literature; the climate-matrix,
acceleration and emulator literature; and the intermediate-complexity
atmospheric model literature.

Caveats carried forward:

- The local `climber-x` tree is **modified** relative to `climber-x-orig`
  (implicit diffusion, budget checks, implicit `sam` solve, convergence filter).
  Cite `climber-x-orig` or the papers, not this tree.
- The local `rembo` tree is **mid-refactor**: the moisture/condensation path is
  inside `if (.FALSE.)`, advection is zeroed, the REMBO1 path is disabled.
  `HEAD` is not a validated configuration.
- **T42L19 throughput (~380 SYPD/core) is an extrapolation**, not a measurement:
  SpeedyWeather's T31L8 2,300 SYPD scaled by (42/31)³×(19/8) ≈ 6. M0a exists to
  replace this estimate with a number.
- The SPEEDY single-core figure (125 SYPD) came from a source that could not be
  fully verified; SpeedyWeather's own measured ~240 SYPD for Fortran SPEEDY is
  the better citation.
