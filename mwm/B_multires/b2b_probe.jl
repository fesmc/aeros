# b2b_probe.jl  --  Option B, stage b2, APPROACH B: nudging-derived constant forcing.
#
# Approach A (b2_probe.jl) diagnosed the resolution-error forcing from a SINGLE
# dynamical-core tendency evaluation at the T85 time-mean state X85 (ΔF=−RHS_T31(X85)).
# It was INSUFFICIENT: Z* error reduced only ~19% and T*/U* got WORSE, because
#   (1) transient rectification -- instantaneous RHS at the mean state ≠ time-mean RHS
#       (nonlinear eddy terms), and
#   (2) it corrected only the dynamical core; physics responded freely.
#
# Approach B diagnoses the constant forcing by NUDGING (CAPT / Watt-Meyer "nudging
# tendency"): run T31 with a relaxation −(x−X85)/τ added to the prognostic tendencies,
# and time-average the applied nudging tendency ⟨N⟩ over the averaging window. ⟨N⟩ is
# the constant forcing that holds free T31 at ≈X85; averaged over the eddying
# trajectory it captures transient rectification, and by compensating the model's NET
# tendency each step it captures the physics response too. Then run FREE T31 + ⟨N⟩
# (constant forcing, no nudging) and test stationary-wave reproduction, reporting the
# three-way error reduction DIRECTLY alongside approach A's numbers.
#
# Selected-terms variant (§3.7): `correct_fields` selects which prognostic fields are
# nudged/corrected (par/b2b.toml = all; par/b2b_vordiv.toml = vorticity+divergence
# only, to test whether restricting the correction avoids the T*/U* degradation).
#
# Usage:  julia --project=. b2b_probe.jl par/b2b.toml
#
# NB: production T85 is expensive -- compute node only (run_speedy_b2b.sh). Use
# par/b2b_smoke.toml for the login-node plumbing smoke test.

using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using NCDatasets
using TOML
using Printf
using Statistics

const HERE = @__DIR__
include(joinpath(HERE, "b_common.jl"))     # read-only reuse (approach A machinery)
include(joinpath(HERE, "b2b_common.jl"))   # NEW: NudgingForcing + ⟨N⟩ builder

# approach-A production headline (par/b2.toml, from output_b2/b2_summary.toml) for the
# side-by-side comparison. These are FIXED reference numbers, not recomputed here.
const A_FRAC = (Zstar = 0.1948, Zstar_NH = 0.1438, Tstar = -0.323, Ustar = -0.0974)

# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
length(ARGS) >= 1 || error("usage: julia --project=. b2b_probe.jl <config.toml>")
const CFGPATH = ARGS[1]
isfile(CFGPATH) || error("config file not found: $CFGPATH")
const CFG = TOML.parsefile(CFGPATH)

getcfg(k, default) = get(CFG, k, default)
const TRUNC_HI     = Int(getcfg("trunc_hi", 85))
const TRUNC_LO     = Int(getcfg("trunc_lo", 31))
const NLEV         = Int(getcfg("nlev", 8))
const SPINUP_DAYS  = Int(getcfg("spinup_days", 200))
const MEAN_DAYS    = Int(getcfg("mean_days", 360))
const TARGET_SIGMA = Float64(getcfg("target_sigma", 0.5))
const NUDGE_TAU_H  = Float64(getcfg("nudge_tau_hours", 6.0))   # τ [hours]; see caveats
const NUDGE_TAU_S  = NUDGE_TAU_H * 3600.0                      # τ [s]
const CORRECT_FIELDS = String.(getcfg("correct_fields", ["vor", "div", "temp", "humid", "pres"]))
const STAGE        = getcfg("stage", "b2b")
const OUTDIR_REL   = getcfg("output_dir", "output_b2b")
const OUTDIR       = isabspath(OUTDIR_REL) ? OUTDIR_REL : joinpath(HERE, OUTDIR_REL)
mkpath(OUTDIR)

@info "b2b_probe config" stage = STAGE trunc_hi = TRUNC_HI trunc_lo = TRUNC_LO nlev = NLEV spinup_days = SPINUP_DAYS mean_days = MEAN_DAYS target_sigma = TARGET_SIGMA nudge_tau_hours = NUDGE_TAU_H correct_fields = CORRECT_FIELDS output = OUTDIR threads = Threads.nthreads()

# ----------------------------------------------------------------------------
# model construction: IDENTICAL prescribed BCs as b1/b2, only `trunc` (+ optional
# forcing) differ. Non-interactive lower boundary (prescribed SST + sea ice).
# ----------------------------------------------------------------------------
function build_model(trunc; forcing = nothing)
    sg = SpectralGrid(; trunc, nlayers = NLEV)
    model = PrimitiveWetModel(sg;
        ocean   = SeasonalOceanClimatology(sg),
        sea_ice = PrescribedSeaIce(sg),
        forcing = forcing,
    )
    return sg, model
end

# spin up then accumulate the time-mean prognostic spectral state (no forcing here)
function run_and_mean(model; spinup_days, mean_days, label)
    sim = initialize!(model)
    @info "  spinup" label = label days = spinup_days dt = model.time_stepping.Δt_sec
    run!(sim, period = Day(spinup_days), output = false)
    acc = SpectralMeanAccumulator()
    add!(model, :specmean => acc)
    @info "  averaging" label = label days = mean_days
    run!(sim, period = Day(mean_days), output = false)
    return mean_state(acc), sim
end

# nudged diagnosis run. To avoid the cold-start shock of relaxing the near-rest
# default IC toward the fully-structured X85 (which blows the run up), we (1) run a
# plain, un-nudged spinup to reach the T31 climate, (2) switch nudging ON and run a
# second spinup so the state settles near X85, then (3) accumulate BOTH the time-mean
# prognostic state (to verify mean ≈ X85) and the time-mean nudging tendency ⟨N⟩ over
# the averaging window.
function run_nudged_and_mean(model, forcing; spinup_days, mean_days, label)
    sim = initialize!(model)
    set_active!(forcing, false)
    @info "  plain (un-nudged) spinup" label = label days = spinup_days dt = model.time_stepping.Δt_sec
    run!(sim, period = Day(spinup_days), output = false)
    set_active!(forcing, true)
    @info "  nudged spinup (settle toward X85)" label = label days = spinup_days tau_hours = NUDGE_TAU_H
    run!(sim, period = Day(spinup_days), output = false)
    reset_accumulation!(forcing)
    set_accumulate!(forcing, true)
    acc = SpectralMeanAccumulator()
    add!(model, :specmean => acc)
    @info "  nudged averaging (accumulating state mean + ⟨N⟩)" label = label days = mean_days
    run!(sim, period = Day(mean_days), output = false)
    set_accumulate!(forcing, false)
    return mean_state(acc), mean_nudging(forcing), sim
end

# ----------------------------------------------------------------------------
# 1. T85 reference run + coarse-grain to T31 target state X85  (same as b2)
# ----------------------------------------------------------------------------
@info "=== [1] HI resolution T$TRUNC_HI reference ==="
sg_hi, model_hi = build_model(TRUNC_HI)
mean_hi, _ = run_and_mean(model_hi; spinup_days = SPINUP_DAYS, mean_days = MEAN_DAYS, label = "T$TRUNC_HI")

# ----------------------------------------------------------------------------
# 2. T31 diagnosis sim + coarse-grain T85 mean -> X85 (target state)  (same as b2)
# ----------------------------------------------------------------------------
@info "=== [2] build T31 diagnosis model + coarse-grain T85 mean -> X85 ==="
sg_lo, model_lo = build_model(TRUNC_LO)
sim_lo = initialize!(model_lo)
SpeedyWeather.initialize!(sim_lo; period = Day(1), output = false)   # scale prognostic by radius
vars_lo = sim_lo.variables
const RADIUS = Float64(model_lo.planet.radius)
const G = Float64(model_lo.planet.gravity)

const KLEV  = argmin(abs.(Array(model_lo.geometry.σ_levels_full) .- TARGET_SIGMA))
const SIGMA = Float64(model_lo.geometry.σ_levels_full[KLEV])
@info "  target level" k = KLEV sigma = SIGMA approx_hPa = 1000 * SIGMA

function to_T31(field3d_or_2d, template)
    tr = coarsegrain_spectral(field3d_or_2d, TRUNC_LO)
    dst = zero(SpeedyWeather.get_step(template, 2))
    copyto!(dst, tr)
    return dst
end
X85 = (
    vor   = to_T31(mean_hi.vor,   vars_lo.prognostic.vorticity),
    div   = to_T31(mean_hi.div,   vars_lo.prognostic.divergence),
    temp  = to_T31(mean_hi.temp,  vars_lo.prognostic.temperature),
    humid = to_T31(mean_hi.humid, vars_lo.prognostic.humidity),
    pres  = to_T31(mean_hi.pres,  vars_lo.prognostic.pressure),
)

# ----------------------------------------------------------------------------
# 2b. MAGNITUDE REFERENCE: approach A's −RHS_T31(X85), a single stable dynamical-core
#     tendency eval at X85. ⟨N⟩ and −RHS(X85) are the SAME physical quantity (the
#     forcing that holds X85) diagnosed two ways, so they must be the SAME ORDER of
#     magnitude per field. This is the units/scaling sanity gate: a ~radius (6.4e6)
#     mismatch would flag a dropped/double radius factor in the nudge scaling.
# ----------------------------------------------------------------------------
@info "=== [2b] magnitude reference −RHS_T31(X85) (approach A object) ==="
set_state!(vars_lo, X85)
rhs_ref = eval_dynamics!(vars_lo, model_lo, 2)   # RHS(X85); −RHS = approach A's ΔF
refmax = (
    vor   = maximum(abs, Array(rhs_ref.vor)),
    div   = maximum(abs, Array(rhs_ref.div)),
    temp  = maximum(abs, Array(rhs_ref.temp)),
    humid = maximum(abs, Array(rhs_ref.humid)),
    pres  = maximum(abs, Array(rhs_ref.pres)),
)
@info "  max|−RHS_T31(X85)| (scaled tendency units)" vor = refmax.vor div = refmax.div temp = refmax.temp humid = refmax.humid pres = refmax.pres

# ----------------------------------------------------------------------------
# 3. NUDGING diagnosis: run T31 with −(x−X85)/τ on the selected fields, accumulate ⟨N⟩
# ----------------------------------------------------------------------------
@info "=== [3] nudging diagnosis of ⟨N⟩ ===" tau_hours = NUDGE_TAU_H tau_seconds = NUDGE_TAU_S correct_fields = CORRECT_FIELDS
nudging = NudgingForcing(X85, NUDGE_TAU_S, CORRECT_FIELDS)
_, model_nudge = build_model(TRUNC_LO; forcing = nudging)
mean_nudged, Nmean, _ = run_nudged_and_mean(model_nudge, nudging;
    spinup_days = SPINUP_DAYS, mean_days = MEAN_DAYS, label = "T$TRUNC_LO nudged")
@info "  ⟨N⟩ magnitudes (scaled tendency units)" nsteps = Nmean.nsteps max_vor = maximum(abs, Array(Nmean.vor)) max_div = maximum(abs, Array(Nmean.div)) max_temp = maximum(abs, Array(Nmean.temp)) max_humid = maximum(abs, Array(Nmean.humid)) max_pres = maximum(abs, Array(Nmean.pres))

# MAGNITUDE SANITY GATE: max|⟨N⟩| next to max|−RHS(X85)| per field. Same order = OK.
magcmp = Dict{String, Any}()
for (name, nf, rf) in (
        ("vor",   Nmean.vor,   refmax.vor),
        ("div",   Nmean.div,   refmax.div),
        ("temp",  Nmean.temp,  refmax.temp),
        ("humid", Nmean.humid, refmax.humid),
        ("pres",  Nmean.pres,  refmax.pres),
    )
    nmax = maximum(abs, Array(nf))
    ratio = rf > 0 ? nmax / rf : (nmax == 0 ? 0.0 : Inf)
    magcmp[name] = (Nmax = nmax, RHSmax = rf, ratio = ratio)
    @info "  MAG GATE ⟨N⟩ vs −RHS(X85)" field = name maxN = nmax maxNegRHS = rf ratio_N_over_RHS = ratio
end

# ----------------------------------------------------------------------------
# diagnostics grid + helpers (same conventions as b2)
# ----------------------------------------------------------------------------
const TARGET_GRID = FullGaussianGrid(sg_lo.nlat_half)
const NLON = get_nlon_max(TARGET_GRID)
const NLAT = get_nlat(TARGET_GRID)
const LOND = Array(get_lond(TARGET_GRID))
const LATD = Array(get_latd(TARGET_GRID))
const S_LO = model_lo.spectral_transform

# zero the l=m=0 (global-mean) coefficient per layer, so an RMS reflects the spatial
# STRUCTURE rather than a huge DC offset (temperature/log-pₛ have large global means
# that would otherwise dominate the nudged-mean-vs-X85 relative RMS).
function demean_spectral(field)
    g = copy(field)
    if ndims(g) == 1 || size(g, 2) == 1
        g[1] = 0
    else
        for k in 1:size(g, 2)
            g[1, k] = 0
        end
    end
    return g
end

# area-weighted RMS of a spectral field (per level, averaged over levels), on grid
function spectral_area_rms(field)
    if ndims(field) == 1 || size(field, 2) == 1
        M = tomatrix(regrid_field(transform(field, S_LO), sg_lo, TARGET_GRID), NLON, NLAT)
        return area_rms(M, LATD)
    else
        nk = size(field, 2)
        acc = 0.0
        for k in 1:nk
            M = tomatrix(regrid_field(transform(field[:, k], S_LO), sg_lo, TARGET_GRID), NLON, NLAT)
            acc += area_rms(M, LATD)^2
        end
        return sqrt(acc / nk)
    end
end

# ----------------------------------------------------------------------------
# 4. VERIFY the nudged run's time-mean ≈ X85 (nudging actually constrains the state)
# ----------------------------------------------------------------------------
@info "=== [4] verify nudged-mean ≈ X85 ==="
nudged_match = Dict{String, Any}()
for (name, mfield, xfield, corrected) in (
        ("vor",   mean_nudged.vor,   X85.vor,   nudging.do_vor),
        ("div",   mean_nudged.div,   X85.div,   nudging.do_div),
        ("temp",  mean_nudged.temp,  X85.temp,  nudging.do_temp),
        ("humid", mean_nudged.humid, X85.humid, nudging.do_humid),
        ("pres",  mean_nudged.pres,  X85.pres,  nudging.do_pres),
    )
    # compare spatial STRUCTURE (l=m=0 mean removed) so the huge temperature/log-pₛ
    # DC offsets do not swamp the metric.
    diff = demean_spectral(mfield .- xfield)
    xstruct = demean_spectral(xfield)
    rms_diff = spectral_area_rms(diff)
    rms_tgt  = spectral_area_rms(xstruct)
    ratio = rms_tgt > 0 ? rms_diff / rms_tgt : 0.0
    nudged_match[name] = (rms_diff = rms_diff, rms_target = rms_tgt, ratio = ratio, corrected = corrected)
    @info "  nudged-mean vs X85 (structure, DC removed)" field = name corrected = corrected rms_diff = rms_diff rms_target = rms_tgt rel_rms = ratio
end

# ----------------------------------------------------------------------------
# 5. build the constant forcing ⟨N⟩ and VERIFY it is injected on free T31
# ----------------------------------------------------------------------------
@info "=== [5] build ⟨N⟩ constant forcing + verify injection ==="
forcing_const = build_forcing_from_meanN(Nmean, model_lo)

# unforced tendency at X85, then forced tendency at X85; the forced−unforced delta
# must equal ⟨N⟩ (in spectral space, temp/humid within the grid round-trip). This is
# the state-independent-forcing analogue of approach A's fixed-point check.
set_state!(vars_lo, X85)
rhs0 = eval_dynamics!(vars_lo, model_lo, 2)

_, model_fchk = build_model(TRUNC_LO; forcing = build_forcing_from_meanN(Nmean, model_lo))
sim_fchk = initialize!(model_fchk)
SpeedyWeather.initialize!(sim_fchk; period = Day(1), output = false)
vars_fchk = sim_fchk.variables
set_state!(vars_fchk, X85)
rhs_F = eval_dynamics!(vars_fchk, model_fchk, 2)

inj = Dict{String, Any}()
for (name, d0, dF, nmean) in (
        ("vor",   rhs0.vor,   rhs_F.vor,   Nmean.vor),
        ("div",   rhs0.div,   rhs_F.div,   Nmean.div),
        ("temp",  rhs0.temp,  rhs_F.temp,  Nmean.temp),
        ("humid", rhs0.humid, rhs_F.humid, Nmean.humid),
        ("pres",  rhs0.pres,  rhs_F.pres,  Nmean.pres),
    )
    delta   = Array(dF) .- Array(d0)          # applied forcing effect on the tendency
    target  = Array(nmean)                    # intended ⟨N⟩
    mismatch = maximum(abs, delta .- target)
    scale    = maximum(abs, target)
    frac = scale > 0 ? mismatch / scale : (maximum(abs, delta) == 0 ? 0.0 : 1.0)
    inj[name] = (target_max = scale, applied_max = maximum(abs, delta), mismatch = mismatch, residual_frac = frac)
    @info "  injection check (applied−⟨N⟩)" field = name target_max = scale applied_max = maximum(abs, delta) residual_frac = frac
end

# ----------------------------------------------------------------------------
# 6. bare T31 and FREE T31+⟨N⟩ runs + mean states
# ----------------------------------------------------------------------------
@info "=== [6a] bare T31 run ==="
_, model_bare = build_model(TRUNC_LO)
mean_bare, _ = run_and_mean(model_bare; spinup_days = SPINUP_DAYS, mean_days = MEAN_DAYS, label = "T$TRUNC_LO bare")

@info "=== [6b] FREE T31 + ⟨N⟩ run (constant forcing, NO nudging) ==="
_, model_corr = build_model(TRUNC_LO; forcing = build_forcing_from_meanN(Nmean, model_lo))
mean_corr, _ = run_and_mean(model_corr; spinup_days = SPINUP_DAYS, mean_days = MEAN_DAYS, label = "T$TRUNC_LO+⟨N⟩")

as_state(m) = (vor = m.vor, div = m.div, temp = m.temp, humid = m.humid, pres = m.pres)
X_bare = as_state(mean_bare)
X_corr = as_state(mean_corr)

# ----------------------------------------------------------------------------
# 7. stationary-wave diagnostics for the three cases on a common grid  (same as b2)
# ----------------------------------------------------------------------------
@info "=== [7] stationary-wave diagnostics (Z*, T*, U*) ==="
function level_fields(X)
    set_state!(vars_lo, X)
    eval_dynamics!(vars_lo, model_lo, 2)
    Z = tomatrix(regrid_field(vars_lo.grid.geopotential[:, KLEV], sg_lo, TARGET_GRID), NLON, NLAT) ./ G
    T = tomatrix(regrid_field(vars_lo.grid.temperature[:, KLEV], sg_lo, TARGET_GRID), NLON, NLAT)
    U = tomatrix(regrid_field(vars_lo.grid.u[:, KLEV], sg_lo, TARGET_GRID), NLON, NLAT)
    return (Z = Z, T = T, U = U)
end

f_tgt  = level_fields(X85)
f_bare = level_fields(X_bare)
f_corr = level_fields(X_corr)

Zs_tgt  = eddy(f_tgt.Z);  Zs_bare = eddy(f_bare.Z);  Zs_corr = eddy(f_corr.Z)
Ts_tgt  = eddy(f_tgt.T);  Ts_bare = eddy(f_bare.T);  Ts_corr = eddy(f_corr.T)
Us_tgt  = eddy(f_tgt.U);  Us_bare = eddy(f_bare.U);  Us_corr = eddy(f_corr.U)

resZ_bare = Zs_bare .- Zs_tgt;  resZ_corr = Zs_corr .- Zs_tgt
resT_bare = Ts_bare .- Ts_tgt;  resT_corr = Ts_corr .- Ts_tgt
resU_bare = Us_bare .- Us_tgt;  resU_corr = Us_corr .- Us_tgt

function reduction(res_bare, res_corr)
    rb = area_rms(res_bare, LATD)
    rc = area_rms(res_corr, LATD)
    return (rms_bare = rb, rms_corr = rc, frac_reduction = rb > 0 ? (rb - rc) / rb : 0.0)
end
redZ = reduction(resZ_bare, resZ_corr)
redT = reduction(resT_bare, resT_corr)
redU = reduction(resU_bare, resU_corr)

nh = LATD .> 0
redZ_nh = reduction(resZ_bare[:, nh], resZ_corr[:, nh])

tops_resZ_corr = top_maxima(resZ_corr, LOND, LATD; n = 6)

@info "RESULT Z* error reduction" rms_bare = redZ.rms_bare rms_corr = redZ.rms_corr frac_reduction = redZ.frac_reduction frac_reduction_NH = redZ_nh.frac_reduction approachA = A_FRAC.Zstar
@info "RESULT T* error reduction" rms_bare = redT.rms_bare rms_corr = redT.rms_corr frac_reduction = redT.frac_reduction approachA = A_FRAC.Tstar
@info "RESULT U* error reduction" rms_bare = redU.rms_bare rms_corr = redU.rms_corr frac_reduction = redU.frac_reduction approachA = A_FRAC.Ustar

# ----------------------------------------------------------------------------
# 8. conservation check of ⟨N⟩ (§3.7 risk 2)  (same convention as b2)
# ----------------------------------------------------------------------------
@info "=== [8] conservation check of ⟨N⟩ ==="
cons = Dict{String, Any}()
for (name, dF) in (
        ("vor",   Nmean.vor),
        ("div",   Nmean.div),
        ("temp",  Nmean.temp),
        ("humid", Nmean.humid),
        ("pres",  Nmean.pres),
    )
    gm = spectral_global_mean(dF, S_LO)
    gm_scalar = gm isa AbstractVector ? mean(gm) : gm
    gm_phys = gm_scalar / RADIUS
    rms = spectral_area_rms(dF) / RADIUS
    ratio = rms > 0 ? abs(gm_phys) / rms : 0.0
    cons[name] = (global_mean_phys = gm_phys, rms_phys = rms, ratio = ratio)
    @info "  ⟨N⟩ conservation" field = name global_mean_phys = gm_phys rms_phys = rms ratio_mean_over_rms = ratio
end

# ----------------------------------------------------------------------------
# 9. write NetCDF (three stationary-wave states + residuals)  (same schema as b2)
# ----------------------------------------------------------------------------
const NCPATH = joinpath(OUTDIR, "b2b_fields.nc")
jperm = sortperm(LATD)
lat_asc = LATD[jperm]
permlat(M) = M[:, jperm]

NCDataset(NCPATH, "c") do ds
    defDim(ds, "lon", NLON)
    defDim(ds, "lat", NLAT)
    v = defVar(ds, "lon", Float64, ("lon",)); v[:] = LOND; v.attrib["units"] = "degrees_east"
    v = defVar(ds, "lat", Float64, ("lat",)); v[:] = lat_asc; v.attrib["units"] = "degrees_north"

    function put(name, M, units, longname)
        var = defVar(ds, name, Float32, ("lon", "lat"))
        var[:, :] = Float32.(permlat(M))
        var.attrib["units"] = units
        var.attrib["long_name"] = longname
    end
    put("Zstar_target", Zs_tgt,  "m", "Stationary-wave geop. height, T$TRUNC_HI->T$TRUNC_LO target (X85)")
    put("Zstar_bare",   Zs_bare, "m", "Stationary-wave geop. height, bare T$TRUNC_LO")
    put("Zstar_corr",   Zs_corr, "m", "Stationary-wave geop. height, T$TRUNC_LO + <N> (nudging-derived)")
    put("resZ_bare", resZ_bare, "m", "Z* error bare T$TRUNC_LO minus target")
    put("resZ_corr", resZ_corr, "m", "Z* error T$TRUNC_LO+<N> minus target")
    put("Tstar_target", Ts_tgt,  "K", "Stationary-wave temperature, target (X85)")
    put("Tstar_bare",   Ts_bare, "K", "Stationary-wave temperature, bare T$TRUNC_LO")
    put("Tstar_corr",   Ts_corr, "K", "Stationary-wave temperature, T$TRUNC_LO + <N>")
    put("resT_bare", resT_bare, "K", "T* error bare minus target")
    put("resT_corr", resT_corr, "K", "T* error corrected minus target")
    put("Ustar_target", Us_tgt,  "m/s", "Stationary-wave zonal wind, target (X85)")
    put("Ustar_bare",   Us_bare, "m/s", "Stationary-wave zonal wind, bare T$TRUNC_LO")
    put("Ustar_corr",   Us_corr, "m/s", "Stationary-wave zonal wind, T$TRUNC_LO + <N>")
    put("resU_bare", resU_bare, "m/s", "U* error bare minus target")
    put("resU_corr", resU_corr, "m/s", "U* error corrected minus target")

    ds.attrib["stage"] = STAGE
    ds.attrib["approach"] = "B (nudging-derived constant forcing <N>)"
    ds.attrib["trunc_hi"] = TRUNC_HI
    ds.attrib["trunc_lo"] = TRUNC_LO
    ds.attrib["nlev"] = NLEV
    ds.attrib["sigma_level"] = SIGMA
    ds.attrib["approx_pressure_hPa"] = 1000 * SIGMA
    ds.attrib["spinup_days"] = SPINUP_DAYS
    ds.attrib["mean_days"] = MEAN_DAYS
    ds.attrib["nudge_tau_hours"] = NUDGE_TAU_H
    ds.attrib["correct_fields"] = join(CORRECT_FIELDS, ",")
    ds.attrib["note"] = "b2 approach B; <N> = time-mean nudging tendency (CAPT), injected as constant additive tendency forcing"
end
@info "wrote NetCDF" path = NCPATH

# ----------------------------------------------------------------------------
# 10. summary TOML (with approach-A side-by-side)
# ----------------------------------------------------------------------------
const SUMPATH = joinpath(OUTDIR, "b2b_summary.toml")
r4(x) = round(x; digits = 4)
r3(x) = round(x; digits = 3)
open(SUMPATH, "w") do io
    println(io, "# b2 approach B -- nudging-derived constant forcing <N> -- summary")
    println(io, "stage = \"$STAGE\"")
    println(io, "approach = \"B (nudging tendency / CAPT)\"")
    println(io, "trunc_hi = $TRUNC_HI")
    println(io, "trunc_lo = $TRUNC_LO")
    println(io, "nlev = $NLEV")
    println(io, "spinup_days = $SPINUP_DAYS")
    println(io, "mean_days = $MEAN_DAYS")
    println(io, "nudge_tau_hours = $NUDGE_TAU_H")
    println(io, "correct_fields = [$(join(map(s -> "\"$s\"", CORRECT_FIELDS), ", "))]")
    println(io, "sigma_level = $(r4(SIGMA))")
    println(io, "approx_pressure_hPa = $(round(1000*SIGMA; digits=1))")
    println(io, "nsteps_hi = $(mean_hi.nsteps)")
    println(io, "nsteps_nudged = $(Nmean.nsteps)")
    println(io, "nsteps_bare = $(mean_bare.nsteps)")
    println(io, "nsteps_corr = $(mean_corr.nsteps)")
    println(io, "common_grid = \"FullGaussianGrid nlat_half=$(sg_lo.nlat_half) ($(NLON)x$(NLAT))\"")
    println(io, "radius_m = $RADIUS")
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# HEADLINE: three-way error reduction, approach B vs approach A.")
    println(io, "# frac_reduction = (RMS(bare-target) - RMS(corr-target)) / RMS(bare-target)")
    println(io, "# 1.0 = perfect reproduction, 0.0 = no effect, <0 = made worse.")
    println(io, "# approach A (par/b2.toml): Z* +0.195 (NH +0.144), T* -0.323, U* -0.097.")
    println(io, "# ------------------------------------------------------------------")
    println(io, "[error_reduction.Zstar]  # ~$(round(1000*SIGMA))hPa geop. height, metres")
    println(io, "rms_bare_minus_target = $(r3(redZ.rms_bare))")
    println(io, "rms_corr_minus_target = $(r3(redZ.rms_corr))")
    println(io, "frac_reduction        = $(r4(redZ.frac_reduction))")
    println(io, "frac_reduction_NH     = $(r4(redZ_nh.frac_reduction))")
    println(io, "approachA_frac_reduction    = $(A_FRAC.Zstar)")
    println(io, "approachA_frac_reduction_NH = $(A_FRAC.Zstar_NH)")
    println(io)
    println(io, "[error_reduction.Tstar]  # kelvin")
    println(io, "rms_bare_minus_target = $(r4(redT.rms_bare))")
    println(io, "rms_corr_minus_target = $(r4(redT.rms_corr))")
    println(io, "frac_reduction        = $(r4(redT.frac_reduction))")
    println(io, "approachA_frac_reduction = $(A_FRAC.Tstar)")
    println(io)
    println(io, "[error_reduction.Ustar]  # m/s")
    println(io, "rms_bare_minus_target = $(r4(redU.rms_bare))")
    println(io, "rms_corr_minus_target = $(r4(redU.rms_corr))")
    println(io, "frac_reduction        = $(r4(redU.frac_reduction))")
    println(io, "approachA_frac_reduction = $(A_FRAC.Ustar)")
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# MAGNITUDE SANITY GATE: max|<N>| vs max|-RHS_T31(X85)| (approach A's")
    println(io, "# object). Same physical quantity (forcing that holds X85) two ways, so")
    println(io, "# same order of magnitude expected. ratio ~ O(1) = scaling OK; ratio ~")
    println(io, "# radius (6.4e6) would flag a dropped/double radius factor.")
    println(io, "# ------------------------------------------------------------------")
    for name in ("vor", "div", "temp", "humid", "pres")
        c = magcmp[name]
        println(io, "[magnitude_gate.$name]")
        println(io, "max_N          = $(c.Nmax)")
        println(io, "max_negRHS_X85 = $(c.RHSmax)")
        println(io, "ratio_N_over_RHS = $(r4(c.ratio))")
    end
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# NUDGED-MEAN vs X85 verification: does the nudging actually constrain")
    println(io, "# the state? rel_rms = RMS(mean_nudged - X85) / RMS(X85) on the spatial")
    println(io, "# STRUCTURE (l=m=0 global mean removed, so the large temp/log-pₛ DC does")
    println(io, "# not swamp it), per field. Corrected fields should be small; non-")
    println(io, "# corrected fields evolve freely (rel_rms not meaningful -- corrected=false).")
    println(io, "# ------------------------------------------------------------------")
    for name in ("vor", "div", "temp", "humid", "pres")
        c = nudged_match[name]
        println(io, "[nudged_mean_match.$name]")
        println(io, "corrected  = $(c.corrected)")
        println(io, "rms_diff   = $(c.rms_diff)")
        println(io, "rms_target = $(c.rms_target)")
        println(io, "rel_rms    = $(r4(c.ratio))")
    end
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# INJECTION VERIFICATION: forced-minus-unforced dynamical-core tendency")
    println(io, "# at X85 must equal <N>. residual_frac = max|applied - <N>| / max|<N>|.")
    println(io, "# (Unlike approach A, the forced tendency at X85 is NOT ~0 -- <N> is the")
    println(io, "# time-mean nudging tendency, not -RHS(X85); this check confirms the")
    println(io, "# constant forcing is injected, not that X85 is a fixed point.)")
    println(io, "# ------------------------------------------------------------------")
    for name in ("vor", "div", "temp", "humid", "pres")
        c = inj[name]
        println(io, "[injection_check.$name]")
        println(io, "target_max    = $(c.target_max)")
        println(io, "applied_max   = $(c.applied_max)")
        println(io, "mismatch      = $(c.mismatch)")
        println(io, "residual_frac = $(r4(c.residual_frac))")
    end
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# CONSERVATION of <N> (§3.7 risk 2): global area-weighted mean vs RMS,")
    println(io, "# physical per-second units. ratio = |global_mean| / RMS.")
    println(io, "# ------------------------------------------------------------------")
    for name in ("vor", "div", "temp", "humid", "pres")
        c = cons[name]
        println(io, "[conservation.$name]")
        println(io, "global_mean_phys_per_s = $(c.global_mean_phys)")
        println(io, "rms_phys_per_s         = $(c.rms_phys)")
        println(io, "ratio_mean_over_rms    = $(r4(c.ratio))")
    end
    println(io)
    println(io, "[[residual_Zstar_corr_top_maxima]]  # (lon,lat,value_m) where corrected error remains")
    for (lonp, latp, valp) in tops_resZ_corr
        println(io, "point = [$(round(lonp; digits=1)), $(round(latp; digits=1)), $(round(valp; digits=2))]")
    end
    println(io)
    println(io, "# CAVEATS (see README Stage b2 approach B):")
    println(io, "# - tau sensitivity: tau = $(NUDGE_TAU_H) h. Too strong (small tau) over-")
    println(io, "#   constrains -> <N> absorbs genuine dynamics/suppresses variability; too")
    println(io, "#   weak (large tau) -> state not held, <N> degenerates toward approach A.")
    println(io, "# - <N> captures the FULL model-error tendency (dynamics+physics) at X85,")
    println(io, "#   the intended improvement over approach A's dynamical-core-only DeltaF.")
    println(io, "# - selected-terms: only correct_fields are nudged/corrected; the rest")
    println(io, "#   evolve freely (par/b2b_vordiv.toml corrects vor+div only).")
    println(io, "# - same b1 confounds: resolution-dependent hyperdiffusion + timestep,")
    println(io, "#   interactive land; SpeedyWeather != aeros; present-day not LGM;")
    println(io, "#   single ~$(round(1000*SIGMA))hPa level.")
end
@info "wrote summary" path = SUMPATH

# ----------------------------------------------------------------------------
# 11. plotting (decoupled; failure here must not lose numerical outputs)
# ----------------------------------------------------------------------------
try
    include(joinpath(HERE, "plot_b2b.jl"))
    Base.invokelatest(plot_b2b, NCPATH, OUTDIR)
    @info "wrote plots" dir = OUTDIR
catch err
    @warn "plotting failed (numerical outputs are saved; run plot_b2b.jl separately)" exception = (err, catch_backtrace())
end

@info "b2b_probe done" output = OUTDIR
