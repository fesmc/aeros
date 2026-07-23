# b2_probe.jl  --  Option B, stage b2: multi-resolution error CORRECTION test.
#
# b1 established that a coherent, sizable stationary-wave resolution error ΔZ*
# (~27% of the T85 signal) exists between T85 and T31 under identical BCs. b2 tests
# the actual design §3.7 assumption: does injecting the diagnosed error back into
# T31 as a CONSTANT ADDITIVE FORCING on the prognostic tendencies make T31 reproduce
# T85's stationary waves?  (§3.7: "diagnose the model-error operator, do not nudge
# toward a state.")
#
# Pipeline:
#   1. Run T85 (identical prescribed BCs as b1), spin up, accumulate the time-mean
#      of the full prognostic spectral state. Coarse-grain (spectral-truncate) that
#      mean to T31 -> X85 (target state, expressed at T31 resolution).
#   2. Diagnose ΔF = −RHS_T31(X85): set a T31 model to X85, evaluate the dynamical-
#      core tendency ONCE. X85 is not T31's equilibrium so this tendency is nonzero;
#      the constant forcing that makes X85 a fixed point is its negative.
#   3. Inject ΔF as a constant forcing into a T31 run (custom AbstractForcing), and
#      VERIFY the injection is actually applied (forced tendency at X85 ≈ 0).
#   4. Run bare T31 and T31+ΔF with the same spinup/mean protocol; compute the
#      stationary waves Z*(~500hPa), T*, U* for X85 (target), bare T31, T31+ΔF and
#      report the fractional error reduction.
#   5. Conservation check on ΔF (§3.7 risk 2): global (area-weighted) integral of
#      each ΔF field.
#
# Usage:  julia --project=. b2_probe.jl par/b2.toml
#
# NB: production T85 is expensive -- run on a compute node (run_speedy_b2.sh), NOT
# the login node. Use par/b2_smoke.toml for the login-node plumbing smoke test.

using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using NCDatasets
using TOML
using Printf
using Statistics

const HERE = @__DIR__
include(joinpath(HERE, "b_common.jl"))

# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
length(ARGS) >= 1 || error("usage: julia --project=. b2_probe.jl <config.toml>")
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
const STAGE        = getcfg("stage", "b2")
const OUTDIR_REL   = getcfg("output_dir", "output_b2")
const OUTDIR       = isabspath(OUTDIR_REL) ? OUTDIR_REL : joinpath(HERE, OUTDIR_REL)
mkpath(OUTDIR)

@info "b2_probe config" stage = STAGE trunc_hi = TRUNC_HI trunc_lo = TRUNC_LO nlev = NLEV spinup_days = SPINUP_DAYS mean_days = MEAN_DAYS target_sigma = TARGET_SIGMA output = OUTDIR threads = Threads.nthreads()

# ----------------------------------------------------------------------------
# model construction: IDENTICAL prescribed BCs as b1, only `trunc` (and optional
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

# spin up then accumulate the time-mean prognostic spectral state
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

# ----------------------------------------------------------------------------
# 1. T85 reference run + coarse-grain to T31 target state X85
# ----------------------------------------------------------------------------
@info "=== [1] HI resolution T$TRUNC_HI reference ==="
sg_hi, model_hi = build_model(TRUNC_HI)
mean_hi, _ = run_and_mean(model_hi; spinup_days = SPINUP_DAYS, mean_days = MEAN_DAYS, label = "T$TRUNC_HI")

# ----------------------------------------------------------------------------
# 2. T31 diagnosis simulation: put it in the runtime (radius-scaled) state, then
#    build X85 by coarse-graining the T85 mean into T31-shaped spectral fields.
# ----------------------------------------------------------------------------
@info "=== [2] build T31 diagnosis model + coarse-grain T85 mean -> X85 ==="
sg_lo, model_lo = build_model(TRUNC_LO)
sim_lo = initialize!(model_lo)
SpeedyWeather.initialize!(sim_lo; period = Day(1), output = false)   # scale prognostic by radius
vars_lo = sim_lo.variables
const RADIUS = Float64(model_lo.planet.radius)
const G = Float64(model_lo.planet.gravity)

# target sigma level index (nearest full sigma level), same convention as b1
const KLEV  = argmin(abs.(Array(model_lo.geometry.σ_levels_full) .- TARGET_SIGMA))
const SIGMA = Float64(model_lo.geometry.σ_levels_full[KLEV])
@info "  target level" k = KLEV sigma = SIGMA approx_hPa = 1000 * SIGMA

# coarse-grain each T85 mean field to T31 by spectral truncation + copy into a
# T31-shaped zero (copyto! copies the overlapping (l,m) coefficients = truncation)
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

# diagnose ΔF = −RHS_T31(X85): one dynamical-core evaluation at X85 (no forcing yet)
@info "=== [2b] diagnose ΔF = −RHS_T31(X85) (single tendency eval) ==="
set_state!(vars_lo, X85)
rhs = eval_dynamics!(vars_lo, model_lo, 2)
@info "  RHS(X85) magnitudes (scaled units)" max_vor = maximum(abs, Array(rhs.vor)) max_div = maximum(abs, Array(rhs.div)) max_temp = maximum(abs, Array(rhs.temp)) max_humid = maximum(abs, Array(rhs.humid)) max_pres = maximum(abs, Array(rhs.pres))

# ----------------------------------------------------------------------------
# 3. build the constant forcing ΔF and VERIFY it is actually applied
# ----------------------------------------------------------------------------
@info "=== [3] build ΔF forcing + verify injection ==="
forcing = build_forcing(rhs, model_lo)
sg_f, model_f = build_model(TRUNC_LO; forcing = forcing)
sim_f = initialize!(model_f)
SpeedyWeather.initialize!(sim_f; period = Day(1), output = false)
vars_f = sim_f.variables
set_state!(vars_f, X85)
rhs_forced = eval_dynamics!(vars_f, model_f, 2)   # forcing! applied inside dynamics_tendencies!

# per-field: forced tendency should be ≈ 0 (X85 is a fixed point of the forced core)
inj = Dict{String, Any}()
for (name, base, forced) in (
        ("vor",   rhs.vor,   rhs_forced.vor),
        ("div",   rhs.div,   rhs_forced.div),
        ("temp",  rhs.temp,  rhs_forced.temp),
        ("humid", rhs.humid, rhs_forced.humid),
        ("pres",  rhs.pres,  rhs_forced.pres),
    )
    mb = maximum(abs, Array(base))
    mf = maximum(abs, Array(forced))
    frac = mb > 0 ? mf / mb : 0.0
    inj[name] = (base = mb, forced = mf, residual_frac = frac)
    @info "  injection check" field = name max_baseline = mb max_forced = mf residual_frac = frac
end

# ----------------------------------------------------------------------------
# 4. bare T31 and T31+ΔF runs + mean states
# ----------------------------------------------------------------------------
@info "=== [4a] bare T31 run ==="
_, model_bare = build_model(TRUNC_LO)
mean_bare, _ = run_and_mean(model_bare; spinup_days = SPINUP_DAYS, mean_days = MEAN_DAYS, label = "T$TRUNC_LO bare")

@info "=== [4b] T31+ΔF run ==="
_, model_corr = build_model(TRUNC_LO; forcing = build_forcing(rhs, model_lo))
mean_corr, _ = run_and_mean(model_corr; spinup_days = SPINUP_DAYS, mean_days = MEAN_DAYS, label = "T$TRUNC_LO+ΔF")

# express bare/corr means as T31-shaped states (already T31; wrap for set_state!)
as_state(m) = (vor = m.vor, div = m.div, temp = m.temp, humid = m.humid, pres = m.pres)
X_bare = as_state(mean_bare)
X_corr = as_state(mean_corr)

# ----------------------------------------------------------------------------
# 5. stationary-wave diagnostics for the three cases on a common grid
# ----------------------------------------------------------------------------
@info "=== [5] stationary-wave diagnostics (Z*, T*, U*) ==="
const TARGET_GRID = FullGaussianGrid(sg_lo.nlat_half)
const NLON = get_nlon_max(TARGET_GRID)
const NLAT = get_nlat(TARGET_GRID)
const LOND = Array(get_lond(TARGET_GRID))
const LATD = Array(get_latd(TARGET_GRID))

# put a state into the (baseline) T31 diagnosis sim, evaluate to populate grid
# geopotential/temperature/u, and return the level-KLEV fields on the common grid.
function level_fields(X)
    set_state!(vars_lo, X)
    eval_dynamics!(vars_lo, model_lo, 2)   # populates vars_lo.grid.{geopotential,temperature,u}
    Z = tomatrix(regrid_field(vars_lo.grid.geopotential[:, KLEV], sg_lo, TARGET_GRID), NLON, NLAT) ./ G
    T = tomatrix(regrid_field(vars_lo.grid.temperature[:, KLEV], sg_lo, TARGET_GRID), NLON, NLAT)
    U = tomatrix(regrid_field(vars_lo.grid.u[:, KLEV], sg_lo, TARGET_GRID), NLON, NLAT)
    return (Z = Z, T = T, U = U)
end

f_tgt  = level_fields(X85)
f_bare = level_fields(X_bare)
f_corr = level_fields(X_corr)

# stationary-wave (eddy) fields
Zs_tgt  = eddy(f_tgt.Z);  Zs_bare = eddy(f_bare.Z);  Zs_corr = eddy(f_corr.Z)
Ts_tgt  = eddy(f_tgt.T);  Ts_bare = eddy(f_bare.T);  Ts_corr = eddy(f_corr.T)
Us_tgt  = eddy(f_tgt.U);  Us_bare = eddy(f_bare.U);  Us_corr = eddy(f_corr.U)

# residual (error) fields relative to the target
resZ_bare = Zs_bare .- Zs_tgt;  resZ_corr = Zs_corr .- Zs_tgt
resT_bare = Ts_bare .- Ts_tgt;  resT_corr = Ts_corr .- Ts_tgt
resU_bare = Us_bare .- Us_tgt;  resU_corr = Us_corr .- Us_tgt

# error-reduction metric: RMS(residual) bare vs corrected, and fractional reduction
function reduction(res_bare, res_corr)
    rb = area_rms(res_bare, LATD)
    rc = area_rms(res_corr, LATD)
    return (rms_bare = rb, rms_corr = rc,
            frac_reduction = rb > 0 ? (rb - rc) / rb : 0.0)
end
redZ = reduction(resZ_bare, resZ_corr)
redT = reduction(resT_bare, resT_corr)
redU = reduction(resU_bare, resU_corr)

# NH-only for Z* (stationary-wave/orographic signal lives in the NH)
nh = LATD .> 0
redZ_nh = reduction(resZ_bare[:, nh], resZ_corr[:, nh])

# where residual error remains (top maxima of corrected Z* residual)
tops_resZ_corr = top_maxima(resZ_corr, LOND, LATD; n = 6)

@info "RESULT Z* error reduction" rms_bare = redZ.rms_bare rms_corr = redZ.rms_corr frac_reduction = redZ.frac_reduction frac_reduction_NH = redZ_nh.frac_reduction
@info "RESULT T* error reduction" rms_bare = redT.rms_bare rms_corr = redT.rms_corr frac_reduction = redT.frac_reduction
@info "RESULT U* error reduction" rms_bare = redU.rms_bare rms_corr = redU.rms_corr frac_reduction = redU.frac_reduction

# ----------------------------------------------------------------------------
# 6. conservation check of ΔF (§3.7 risk 2)
# ----------------------------------------------------------------------------
# ΔF = −rhs (spectral). Report the area-weighted global mean (from the l=m=0
# harmonic) and area-weighted RMS of each ΔF field, in PHYSICAL per-second units
# (the model integrates radius-scaled tendencies, so physical = scaled / radius).
# The dimensionless ratio |mean|/RMS is the scale-invariant "zero-mean" metric.
@info "=== [6] conservation check of ΔF ==="
const S_LO = model_lo.spectral_transform

# area-weighted RMS of a spectral field (per level, then averaged over levels)
function spectral_area_rms(field)
    if ndims(field) == 1 || size(field, 2) == 1
        M = tomatrix(regrid_field(transform(field, S_LO), sg_lo, TARGET_GRID), NLON, NLAT)
        return area_rms(M, LATD)
    else
        nk = size(field, 2)
        acc = 0.0
        for k in 1:nk
            fk = field[:, k]
            M = tomatrix(regrid_field(transform(fk, S_LO), sg_lo, TARGET_GRID), NLON, NLAT)
            acc += area_rms(M, LATD)^2
        end
        return sqrt(acc / nk)
    end
end

cons = Dict{String, Any}()
for (name, dF) in (
        ("vor",   -rhs.vor),
        ("div",   -rhs.div),
        ("temp",  -rhs.temp),
        ("humid", -rhs.humid),
        ("pres",  -rhs.pres),
    )
    gm = spectral_global_mean(dF, S_LO)              # scaled units, per level or scalar
    gm_scalar = gm isa AbstractVector ? mean(gm) : gm
    gm_phys = gm_scalar / RADIUS                      # -> physical per-second
    rms = spectral_area_rms(dF) / RADIUS              # physical per-second
    ratio = rms > 0 ? abs(gm_phys) / rms : 0.0
    cons[name] = (global_mean_phys = gm_phys, rms_phys = rms, ratio = ratio)
    @info "  ΔF conservation" field = name global_mean_phys = gm_phys rms_phys = rms ratio_mean_over_rms = ratio
end

# ----------------------------------------------------------------------------
# 7. write NetCDF (three stationary-wave states + residuals)
# ----------------------------------------------------------------------------
const NCPATH = joinpath(OUTDIR, "b2_fields.nc")
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
    put("Zstar_corr",   Zs_corr, "m", "Stationary-wave geop. height, T$TRUNC_LO + ΔF")
    put("resZ_bare", resZ_bare, "m", "Z* error bare T$TRUNC_LO minus target")
    put("resZ_corr", resZ_corr, "m", "Z* error T$TRUNC_LO+ΔF minus target")
    put("Tstar_target", Ts_tgt,  "K", "Stationary-wave temperature, target (X85)")
    put("Tstar_bare",   Ts_bare, "K", "Stationary-wave temperature, bare T$TRUNC_LO")
    put("Tstar_corr",   Ts_corr, "K", "Stationary-wave temperature, T$TRUNC_LO + ΔF")
    put("resT_bare", resT_bare, "K", "T* error bare minus target")
    put("resT_corr", resT_corr, "K", "T* error corrected minus target")
    put("Ustar_target", Us_tgt,  "m/s", "Stationary-wave zonal wind, target (X85)")
    put("Ustar_bare",   Us_bare, "m/s", "Stationary-wave zonal wind, bare T$TRUNC_LO")
    put("Ustar_corr",   Us_corr, "m/s", "Stationary-wave zonal wind, T$TRUNC_LO + ΔF")
    put("resU_bare", resU_bare, "m/s", "U* error bare minus target")
    put("resU_corr", resU_corr, "m/s", "U* error corrected minus target")

    ds.attrib["stage"] = STAGE
    ds.attrib["trunc_hi"] = TRUNC_HI
    ds.attrib["trunc_lo"] = TRUNC_LO
    ds.attrib["nlev"] = NLEV
    ds.attrib["sigma_level"] = SIGMA
    ds.attrib["approx_pressure_hPa"] = 1000 * SIGMA
    ds.attrib["spinup_days"] = SPINUP_DAYS
    ds.attrib["mean_days"] = MEAN_DAYS
    ds.attrib["note"] = "b2 correction test; ΔF = −RHS_T31(X85) injected as constant additive tendency forcing"
end
@info "wrote NetCDF" path = NCPATH

# ----------------------------------------------------------------------------
# 8. summary TOML
# ----------------------------------------------------------------------------
const SUMPATH = joinpath(OUTDIR, "b2_summary.toml")
r4(x) = round(x; digits = 4)
r3(x) = round(x; digits = 3)
open(SUMPATH, "w") do io
    println(io, "# b2 multi-resolution error CORRECTION test -- summary")
    println(io, "stage = \"$STAGE\"")
    println(io, "trunc_hi = $TRUNC_HI")
    println(io, "trunc_lo = $TRUNC_LO")
    println(io, "nlev = $NLEV")
    println(io, "spinup_days = $SPINUP_DAYS")
    println(io, "mean_days = $MEAN_DAYS")
    println(io, "sigma_level = $(r4(SIGMA))")
    println(io, "approx_pressure_hPa = $(round(1000*SIGMA; digits=1))")
    println(io, "nsteps_hi = $(mean_hi.nsteps)")
    println(io, "nsteps_bare = $(mean_bare.nsteps)")
    println(io, "nsteps_corr = $(mean_corr.nsteps)")
    println(io, "common_grid = \"FullGaussianGrid nlat_half=$(sg_lo.nlat_half) ($(NLON)x$(NLAT))\"")
    println(io, "radius_m = $RADIUS")
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# HEADLINE: does ΔF move T31 toward T85? (fractional RMS error reduction)")
    println(io, "# frac_reduction = (RMS(bare-target) - RMS(corr-target)) / RMS(bare-target)")
    println(io, "# 1.0 = perfect reproduction, 0.0 = no effect, <0 = made worse")
    println(io, "# ------------------------------------------------------------------")
    println(io, "[error_reduction.Zstar]  # ~$(round(1000*SIGMA))hPa geop. height, metres")
    println(io, "rms_bare_minus_target = $(r3(redZ.rms_bare))")
    println(io, "rms_corr_minus_target = $(r3(redZ.rms_corr))")
    println(io, "frac_reduction        = $(r4(redZ.frac_reduction))")
    println(io, "frac_reduction_NH     = $(r4(redZ_nh.frac_reduction))")
    println(io)
    println(io, "[error_reduction.Tstar]  # kelvin")
    println(io, "rms_bare_minus_target = $(r4(redT.rms_bare))")
    println(io, "rms_corr_minus_target = $(r4(redT.rms_corr))")
    println(io, "frac_reduction        = $(r4(redT.frac_reduction))")
    println(io)
    println(io, "[error_reduction.Ustar]  # m/s")
    println(io, "rms_bare_minus_target = $(r4(redU.rms_bare))")
    println(io, "rms_corr_minus_target = $(r4(redU.rms_corr))")
    println(io, "frac_reduction        = $(r4(redU.frac_reduction))")
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# INJECTION VERIFICATION: forced dynamical-core tendency at X85 should")
    println(io, "# be ~0 (residual_frac = max|forced| / max|baseline|). vor/div/pres are")
    println(io, "# injected in spectral space; temp/humid in grid space (see b_common.jl).")
    println(io, "# ------------------------------------------------------------------")
    for name in ("vor", "div", "temp", "humid", "pres")
        c = inj[name]
        println(io, "[injection_check.$name]")
        println(io, "max_baseline  = $(c.base)")
        println(io, "max_forced    = $(c.forced)")
        println(io, "residual_frac = $(r4(c.residual_frac))")
    end
    println(io)
    println(io, "# ------------------------------------------------------------------")
    println(io, "# CONSERVATION of ΔF (§3.7 risk 2): global area-weighted mean vs RMS,")
    println(io, "# physical per-second units (scaled tendency / radius). A pure")
    println(io, "# stationary-wave correction should be ~zero global mean.")
    println(io, "# ratio = |global_mean| / RMS (dimensionless, scale-invariant).")
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
    println(io, "# CAVEATS (see README Stage b2):")
    println(io, "# - single-tendency-at-mean-state neglects transient rectification")
    println(io, "#   (instantaneous RHS at time-mean state != time-mean RHS); leading")
    println(io, "#   candidate if reproduction is partial. The iterated / mean-tendency")
    println(io, "#   version would be the next step.")
    println(io, "# - full-field correction (not §3.7's 'selected terms' refinement).")
    println(io, "# - resolution-dependent hyperdiffusion + timestep, interactive land")
    println(io, "#   (same b1 confounds); ΔF conflates dynamics with resolution numerics.")
    println(io, "# - SpeedyWeather != aeros -> suggestive, not conclusive.")
    println(io, "# - ΔF here is the DYNAMICAL-CORE tendency error (dynamics_tendencies!);")
    println(io, "#   physics parameterizations respond freely each step and are not part")
    println(io, "#   of ΔF, so any residual physics tendency at X85 also limits reproduction.")
end
@info "wrote summary" path = SUMPATH

# ----------------------------------------------------------------------------
# 9. plotting (decoupled; failure here must not lose numerical outputs)
# ----------------------------------------------------------------------------
try
    include(joinpath(HERE, "plot_b2.jl"))
    Base.invokelatest(plot_b2, NCPATH, OUTDIR)
    @info "wrote plots" dir = OUTDIR
catch err
    @warn "plotting failed (numerical outputs are saved; run plot_b2.jl separately)" exception = (err, catch_backtrace())
end

@info "b2_probe done" output = OUTDIR
