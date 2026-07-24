# b1_probe.jl  --  Option B, stage b1: multi-resolution error DIAGNOSIS
#
# Runs a moist primitive-equation atmosphere (SpeedyWeather.jl PrimitiveWetModel)
# at two spectral truncations (hi = T85, lo = T31) under IDENTICAL boundary
# conditions (same real-Earth orography + land-sea mask, prescribed seasonal SST,
# prescribed sea ice -- non-interactive lower boundary). The only intended
# difference between the two runs is the spectral truncation.
#
# It accumulates a time-mean of the ~500 hPa geopotential, temperature and wind,
# coarse-grains the hi-res mean to the lo-res spectral representation, and
# diagnoses the resolution error in the stationary-wave (zonal-anomaly) field.
#
# b1 DIAGNOSES ONLY -- it does NOT apply any correction (that is b2).
#
# Usage:  julia --project=. b1_probe.jl par/b1.toml

using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using NCDatasets
using TOML
using Printf
using Statistics

# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
length(ARGS) >= 1 || error("usage: julia --project=. b1_probe.jl <config.toml>")
const CFGPATH = ARGS[1]
isfile(CFGPATH) || error("config file not found: $CFGPATH")
const CFG = TOML.parsefile(CFGPATH)

getcfg(k, default) = get(CFG, k, default)
const TRUNC_HI    = Int(getcfg("trunc_hi", 85))
const TRUNC_LO    = Int(getcfg("trunc_lo", 31))
const NLEV        = Int(getcfg("nlev", 8))
const SPINUP_DAYS = Int(getcfg("spinup_days", 200))
const MEAN_DAYS   = Int(getcfg("mean_days", 360))
const TARGET_SIGMA = Float64(getcfg("target_sigma", 0.5))   # ~500 hPa
const STAGE       = getcfg("stage", "b1")
const OUTDIR_REL  = getcfg("output_dir", "output")
const HERE        = @__DIR__
const OUTDIR      = isabspath(OUTDIR_REL) ? OUTDIR_REL : joinpath(HERE, OUTDIR_REL)
mkpath(OUTDIR)

@info "b1_probe config" stage=STAGE trunc_hi=TRUNC_HI trunc_lo=TRUNC_LO nlev=NLEV spinup_days=SPINUP_DAYS mean_days=MEAN_DAYS target_sigma=TARGET_SIGMA output=OUTDIR threads=Threads.nthreads()

# ----------------------------------------------------------------------------
# time-mean accumulator callback (one sigma level, fields Phi, T, u, v)
# ----------------------------------------------------------------------------
mutable struct LevelMeanAccumulator <: SpeedyWeather.AbstractCallback
    target_sigma::Float64
    k::Int
    n::Int
    Phi::Any
    T::Any
    u::Any
    v::Any
end
LevelMeanAccumulator(; target_sigma=0.5) =
    LevelMeanAccumulator(target_sigma, 0, 0, nothing, nothing, nothing, nothing)

function SpeedyWeather.initialize!(cb::LevelMeanAccumulator, vars, model)
    cb.k = argmin(abs.(Array(model.geometry.σ_levels_full) .- cb.target_sigma))
    cb.n = 0
    cb.Phi = zero(vars.grid.geopotential[:, cb.k])
    cb.T   = zero(vars.grid.temperature[:, cb.k])
    cb.u   = zero(vars.grid.u[:, cb.k])
    cb.v   = zero(vars.grid.v[:, cb.k])
    return nothing
end

function SpeedyWeather.callback!(cb::LevelMeanAccumulator, vars, model)
    cb.n += 1
    cb.Phi .+= vars.grid.geopotential[:, cb.k]
    cb.T   .+= vars.grid.temperature[:, cb.k]
    cb.u   .+= vars.grid.u[:, cb.k]
    cb.v   .+= vars.grid.v[:, cb.k]
    return nothing
end

SpeedyWeather.finalize!(::LevelMeanAccumulator, args...) = nothing

# ----------------------------------------------------------------------------
# build + run one resolution, return time-mean level fields on native grid
# ----------------------------------------------------------------------------
function run_resolution(trunc, nlev, spinup_days, mean_days, target_sigma)
    spectral_grid = SpectralGrid(; trunc, nlayers=nlev)

    # IDENTICAL BOUNDARY CONDITIONS across resolutions, the only difference is `trunc`:
    #   orography      = EarthOrography  (default, real Earth, spectrally truncated at own res)
    #   land_sea_mask  = EarthLandSeaMask (default, real Earth)
    #   ocean          = SeasonalOceanClimatology -> PRESCRIBED seasonal SST from assets (non-interactive)
    #   sea_ice        = PrescribedSeaIce         -> non-interactive sea ice
    # (default SlabOcean/ThermodynamicSeaIce would let the two runs diverge through ocean state.)
    model = PrimitiveWetModel(spectral_grid;
        ocean   = SeasonalOceanClimatology(spectral_grid),
        sea_ice = PrescribedSeaIce(spectral_grid),
    )

    simulation = initialize!(model)

    @info "  spinup" trunc=trunc days=spinup_days dt=model.time_stepping.Δt_sec
    run!(simulation, period=Day(spinup_days), output=false)

    acc = LevelMeanAccumulator(; target_sigma)
    add!(model, :levelmean => acc)
    @info "  averaging" trunc=trunc days=mean_days
    run!(simulation, period=Day(mean_days), output=false)

    acc.n > 0 || error("no timesteps accumulated at trunc=$trunc")
    # divide by the Int step count directly: Float32 field / Int stays Float32,
    # keeping the number format consistent with the (Float32) spectral transforms.
    return (
        spectral_grid = spectral_grid,
        gravity = Float64(model.planet.gravity),
        k = acc.k,
        sigma = Float64(model.geometry.σ_levels_full[acc.k]),
        nsteps = acc.n,
        Phi = acc.Phi ./ acc.n,     # geopotential [m^2/s^2] on native grid
        T   = acc.T   ./ acc.n,     # temperature [K]
        u   = acc.u   ./ acc.n,     # zonal wind [m/s]
        v   = acc.v   ./ acc.n,     # meridional wind [m/s]
    )
end

# ----------------------------------------------------------------------------
# spectral machinery to bring both means onto a common grid (full Gaussian @ lo)
# ----------------------------------------------------------------------------
# lo-res field -> common target grid (no truncation, just regrid via spectral)
function regrid_lo(field, sg_lo, target_grid)
    S_src = SpectralTransform(sg_lo.spectrum, sg_lo.grid; NF=Float32, nlayers=1)
    spec  = transform(field, S_src)
    S_tgt = SpectralTransform(spec.spectrum, target_grid; NF=Float32, nlayers=1)
    return transform(spec, S_tgt)
end

# hi-res field -> spectrally truncate to lo truncation -> common target grid
function coarsegrain_hi(field, sg_hi, trunc_lo, target_grid)
    S_src   = SpectralTransform(sg_hi.spectrum, sg_hi.grid; NF=Float32, nlayers=1)
    spec    = transform(field, S_src)
    spec_tr = spectral_truncation(spec, trunc_lo)   # low-pass to degree trunc_lo (= T{trunc_lo})
    S_tgt   = SpectralTransform(spec_tr.spectrum, target_grid; NF=Float32, nlayers=1)
    return transform(spec_tr, S_tgt)
end

# ----------------------------------------------------------------------------
# matrix (nlon x nlat) helpers for diagnostics on the common full-Gaussian grid
# ----------------------------------------------------------------------------
tomatrix(field, nlon, nlat) = reshape(Array(field), nlon, nlat)   # [lon, lat], lat 90->-90

# stationary-wave (eddy) = deviation from zonal mean; returns [lon,lat] matrix
function eddy(M)
    zm = mean(M; dims=1)          # 1 x nlat zonal mean
    return M .- zm
end

# cos-lat area-weighted RMS over a [lon,nlat] matrix
function area_rms(M, latd)
    cw = cosd.(latd)              # length nlat weights
    num = 0.0; den = 0.0
    @inbounds for j in axes(M, 2)
        w = cw[j]
        s = 0.0
        for i in axes(M, 1)
            s += M[i, j]^2
        end
        num += w * s
        den += w * size(M, 1)
    end
    return sqrt(num / den)
end

# location (lon,lat) and value of the largest |field| point
function argmax_loc(M, lond, latd)
    idx = argmax(abs.(M))
    i, j = Tuple(idx)
    return (lon=lond[i], lat=latd[j], value=M[i, j])
end

# top-N |field| maxima as list of (lon,lat,value), coarsely de-duplicated
function top_maxima(M, lond, latd; n=5, mindeg=15.0)
    A = abs.(M)
    picks = Tuple{Float64,Float64,Float64}[]
    Awork = copy(A)
    for _ in 1:n
        idx = argmax(Awork)
        i, j = Tuple(idx)
        push!(picks, (lond[i], latd[j], M[i, j]))
        # suppress a neighbourhood so we don't return the same feature n times
        for jj in axes(Awork, 2), ii in axes(Awork, 1)
            dlat = abs(latd[jj] - latd[j])
            dlon = abs(lond[ii] - lond[i]); dlon = min(dlon, 360 - dlon)
            if dlat < mindeg && dlon < mindeg
                Awork[ii, jj] = 0.0
            end
        end
    end
    return picks
end

# ----------------------------------------------------------------------------
# run both resolutions
# ----------------------------------------------------------------------------
@info "=== HI resolution T$TRUNC_HI ==="
hi = run_resolution(TRUNC_HI, NLEV, SPINUP_DAYS, MEAN_DAYS, TARGET_SIGMA)
@info "=== LO resolution T$TRUNC_LO ==="
lo = run_resolution(TRUNC_LO, NLEV, SPINUP_DAYS, MEAN_DAYS, TARGET_SIGMA)

# common comparison grid: full Gaussian grid at the LO resolution
const TARGET_GRID = FullGaussianGrid(lo.spectral_grid.nlat_half)
const NLON = get_nlon_max(TARGET_GRID)
const NLAT = get_nlat(TARGET_GRID)
const LOND = Array(get_lond(TARGET_GRID))
const LATD = Array(get_latd(TARGET_GRID))
@assert length(LOND) == NLON
@assert length(LATD) == NLAT
const G = lo.gravity   # identical Earth gravity at both resolutions

@info "common grid" nlon=NLON nlat=NLAT k_hi=hi.k sigma_hi=hi.sigma k_lo=lo.k sigma_lo=lo.sigma nsteps_hi=hi.nsteps nsteps_lo=lo.nsteps

# bring each field onto the common grid, as [lon,lat] matrices
# geopotential height Z = Phi / g  [m]
Z_hi = tomatrix(coarsegrain_hi(hi.Phi, hi.spectral_grid, TRUNC_LO, TARGET_GRID), NLON, NLAT) ./ G
Z_lo = tomatrix(regrid_lo(lo.Phi, lo.spectral_grid, TARGET_GRID),               NLON, NLAT) ./ G
T_hi = tomatrix(coarsegrain_hi(hi.T, hi.spectral_grid, TRUNC_LO, TARGET_GRID),  NLON, NLAT)
T_lo = tomatrix(regrid_lo(lo.T, lo.spectral_grid, TARGET_GRID),                 NLON, NLAT)
U_hi = tomatrix(coarsegrain_hi(hi.u, hi.spectral_grid, TRUNC_LO, TARGET_GRID),  NLON, NLAT)
U_lo = tomatrix(regrid_lo(lo.u, lo.spectral_grid, TARGET_GRID),                 NLON, NLAT)

# stationary-wave (eddy) fields
Zstar_hi = eddy(Z_hi); Zstar_lo = eddy(Z_lo); dZstar = Zstar_hi .- Zstar_lo
Tstar_hi = eddy(T_hi); Tstar_lo = eddy(T_lo); dTstar = Tstar_hi .- Tstar_lo
Ustar_hi = eddy(U_hi); Ustar_lo = eddy(U_lo); dUstar = Ustar_hi .- Ustar_lo

# full-field (zonal mean incl.) differences too, for T and u
dZ = Z_hi .- Z_lo
dT = T_hi .- T_lo
dU = U_hi .- U_lo

# ----------------------------------------------------------------------------
# quantify
# ----------------------------------------------------------------------------
rms_Zstar_hi = area_rms(Zstar_hi, LATD)
rms_Zstar_lo = area_rms(Zstar_lo, LATD)
rms_dZstar   = area_rms(dZstar,   LATD)
ratio_Z      = rms_dZstar / rms_Zstar_hi

rms_Tstar_hi = area_rms(Tstar_hi, LATD)
rms_Tstar_lo = area_rms(Tstar_lo, LATD)
rms_dTstar   = area_rms(dTstar,   LATD)
ratio_T      = rms_dTstar / rms_Tstar_hi

rms_Ustar_hi = area_rms(Ustar_hi, LATD)
rms_Ustar_lo = area_rms(Ustar_lo, LATD)
rms_dUstar   = area_rms(dUstar,   LATD)
ratio_U      = rms_dUstar / rms_Ustar_hi

rms_dZ = area_rms(dZ, LATD)
rms_dT = area_rms(dT, LATD)
rms_dU = area_rms(dU, LATD)

# NH-only ratio for Z* (the stationary-wave / orographic signal lives in the NH)
nh = LATD .> 0
rms_dZstar_nh   = area_rms(dZstar[:, nh],   LATD[nh])
rms_Zstar_hi_nh = area_rms(Zstar_hi[:, nh], LATD[nh])
ratio_Z_nh      = rms_dZstar_nh / rms_Zstar_hi_nh

maxloc_dZstar = argmax_loc(dZstar, LOND, LATD)
tops_dZstar   = top_maxima(dZstar, LOND, LATD; n=6)

@info "RESULT Z* (~500hPa geopotential height, m)" rms_hi=rms_Zstar_hi rms_lo=rms_Zstar_lo rms_dZstar=rms_dZstar ratio=ratio_Z ratio_NH=ratio_Z_nh
@info "RESULT dZ* max" lon=maxloc_dZstar.lon lat=maxloc_dZstar.lat value=maxloc_dZstar.value

# ----------------------------------------------------------------------------
# save NetCDF of the stationary-wave fields + differences on the common grid
# ----------------------------------------------------------------------------
const NCPATH = joinpath(OUTDIR, "b1_fields.nc")
# store latitudes ascending (-90 -> 90) for tidy plotting
jperm = sortperm(LATD)
lat_asc = LATD[jperm]
permlat(M) = permutedims(M[:, jperm], (1, 2))   # keep [lon,lat] with lat ascending

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
    put("Zstar_hi", Zstar_hi, "m",   "Stationary-wave geop. height, T$TRUNC_HI coarse-grained to T$TRUNC_LO")
    put("Zstar_lo", Zstar_lo, "m",   "Stationary-wave geop. height, T$TRUNC_LO")
    put("dZstar",   dZstar,   "m",   "Resolution error in Z* (hi->lo minus lo)")
    put("Tstar_hi", Tstar_hi, "K",   "Stationary-wave temperature, T$TRUNC_HI -> T$TRUNC_LO")
    put("Tstar_lo", Tstar_lo, "K",   "Stationary-wave temperature, T$TRUNC_LO")
    put("dTstar",   dTstar,   "K",   "Resolution error in T*")
    put("Ustar_hi", Ustar_hi, "m/s", "Stationary-wave zonal wind, T$TRUNC_HI -> T$TRUNC_LO")
    put("Ustar_lo", Ustar_lo, "m/s", "Stationary-wave zonal wind, T$TRUNC_LO")
    put("dUstar",   dUstar,   "m/s", "Resolution error in u*")
    put("Z_hi",     Z_hi,     "m",   "Time-mean geop. height, T$TRUNC_HI -> T$TRUNC_LO")
    put("Z_lo",     Z_lo,     "m",   "Time-mean geop. height, T$TRUNC_LO")
    put("dZ",       dZ,       "m",   "Full-field resolution error in Z")

    ds.attrib["stage"] = STAGE
    ds.attrib["trunc_hi"] = TRUNC_HI
    ds.attrib["trunc_lo"] = TRUNC_LO
    ds.attrib["nlev"] = NLEV
    ds.attrib["sigma_level"] = lo.sigma
    ds.attrib["approx_pressure_hPa"] = 1000 * lo.sigma
    ds.attrib["spinup_days"] = SPINUP_DAYS
    ds.attrib["mean_days"] = MEAN_DAYS
    ds.attrib["note"] = "b1 diagnosis only; identical BCs (prescribed SST/sea-ice), only truncation differs"
end
@info "wrote NetCDF" path=NCPATH

# ----------------------------------------------------------------------------
# summary TOML
# ----------------------------------------------------------------------------
const SUMPATH = joinpath(OUTDIR, "b1_summary.toml")
open(SUMPATH, "w") do io
    println(io, "# b1 multi-resolution error diagnosis -- summary")
    println(io, "stage = \"$STAGE\"")
    println(io, "trunc_hi = $TRUNC_HI")
    println(io, "trunc_lo = $TRUNC_LO")
    println(io, "nlev = $NLEV")
    println(io, "spinup_days = $SPINUP_DAYS")
    println(io, "mean_days = $MEAN_DAYS")
    println(io, "sigma_level = $(round(lo.sigma; digits=4))")
    println(io, "approx_pressure_hPa = $(round(1000*lo.sigma; digits=1))")
    println(io, "nsteps_hi = $(hi.nsteps)")
    println(io, "nsteps_lo = $(lo.nsteps)")
    println(io, "common_grid = \"FullGaussianGrid nlat_half=$(lo.spectral_grid.nlat_half) ($(NLON)x$(NLAT))\"")
    println(io)
    println(io, "[geopotential_height_stationary_wave]  # metres")
    println(io, "rms_Zstar_hi = $(round(rms_Zstar_hi; digits=3))")
    println(io, "rms_Zstar_lo = $(round(rms_Zstar_lo; digits=3))")
    println(io, "rms_dZstar   = $(round(rms_dZstar; digits=3))")
    println(io, "ratio_dZstar_over_Zstar_hi = $(round(ratio_Z; digits=4))")
    println(io, "rms_dZstar_NH = $(round(rms_dZstar_nh; digits=3))")
    println(io, "ratio_dZstar_over_Zstar_hi_NH = $(round(ratio_Z_nh; digits=4))")
    println(io, "rms_dZ_fullfield = $(round(rms_dZ; digits=3))")
    println(io)
    println(io, "[temperature_stationary_wave]  # kelvin")
    println(io, "rms_Tstar_hi = $(round(rms_Tstar_hi; digits=4))")
    println(io, "rms_Tstar_lo = $(round(rms_Tstar_lo; digits=4))")
    println(io, "rms_dTstar   = $(round(rms_dTstar; digits=4))")
    println(io, "ratio_dTstar_over_Tstar_hi = $(round(ratio_T; digits=4))")
    println(io, "rms_dT_fullfield = $(round(rms_dT; digits=4))")
    println(io)
    println(io, "[zonal_wind_stationary_wave]  # m/s")
    println(io, "rms_Ustar_hi = $(round(rms_Ustar_hi; digits=4))")
    println(io, "rms_Ustar_lo = $(round(rms_Ustar_lo; digits=4))")
    println(io, "rms_dUstar   = $(round(rms_dUstar; digits=4))")
    println(io, "ratio_dUstar_over_Ustar_hi = $(round(ratio_U; digits=4))")
    println(io, "rms_dU_fullfield = $(round(rms_dU; digits=4))")
    println(io)
    println(io, "[dZstar_max]  # single largest |resolution error| in Z*")
    println(io, "lon = $(round(maxloc_dZstar.lon; digits=2))")
    println(io, "lat = $(round(maxloc_dZstar.lat; digits=2))")
    println(io, "value_m = $(round(maxloc_dZstar.value; digits=2))")
    println(io)
    println(io, "[[dZstar_top_maxima]]  # (lon, lat, value_m), de-duplicated features")
    for (lonp, latp, valp) in tops_dZstar
        println(io, "point = [$(round(lonp; digits=1)), $(round(latp; digits=1)), $(round(valp; digits=2))]")
    end
end
@info "wrote summary" path=SUMPATH

# ----------------------------------------------------------------------------
# plotting (decoupled; failure here must not lose the numerical outputs above)
# ----------------------------------------------------------------------------
try
    include(joinpath(HERE, "plot_b1.jl"))
    Base.invokelatest(plot_b1, NCPATH, OUTDIR)
    @info "wrote plots" dir=OUTDIR
catch err
    @warn "plotting failed (numerical outputs are saved; run plot_b1.jl separately)" exception=(err, catch_backtrace())
end

@info "b1_probe done" output=OUTDIR
