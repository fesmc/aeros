# sw_realearth.jl -- SpeedyWeather.jl PrimitiveWetModel at real-Earth boundary
# conditions, default physics, as the 1-1 yardstick for aeros (docs/m2_handoff.md,
# 2026-09-05 "NEXT" item 2).
#
#   truncation TRUNC (default 21), NLAYERS sigma levels (default 8, SW default)
#   orography      = EarthOrography (default)
#   land_sea_mask  = EarthLandSeaMask (default)
#   ocean          = SeasonalOceanClimatology  (PRESCRIBED seasonal SST, like aeros)
#   sea_ice        = PrescribedSeaIce
#   land           = default LandModel (soil temperature/moisture prognostic)
#   physics        = SpeedyWeather defaults (BettsMillerConvection, ImplicitCondensation,
#                    one-band SW/LW, BulkRichardsonDiffusion, ...)
#
# Runs SPINUP days, then accumulates a MONTHLY climatology over MEANYR years with
# the same layout aeros' rce_long `monthly_out` writes:
#   t,u,v,rh (lat,lev,month)  zonal means;  precip (lat,month);  t2m (lon,lat,month)
#   precip2d (lon,lat,month);  z500 (lon,lat,month) 500 hPa geopotential height
# plus the wall-clock throughput (SYPD) of the accumulation years.
#
# Usage: julia --project=. -t N sw_realearth.jl [spinup_days] [mean_years] [trunc] [nlayers]

using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using NCDatasets
using Statistics
using Dates
using Printf

const SPINUP  = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 730
const MEANYR  = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
const TRUNC   = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 21
const NLAYERS = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 8
const HERE    = @__DIR__
const OUTDIR  = joinpath(HERE, "output")
const NMON    = 12
const P500    = 50000.0
mkpath(OUTDIR)

println("SpeedyWeather ", pkgversion(SpeedyWeather), "  threads ", Threads.nthreads(),
        "  T", TRUNC, " L", NLAYERS, "  spinup ", SPINUP, " d  mean ", MEANYR, " yr")

# ---------------------------------------------------------------------------
# regrid helpers: native grid -> full Gaussian [lon, lat] (lat 90 -> -90)
# ---------------------------------------------------------------------------
function make_regridder(sg)
    target = FullGaussianGrid(sg.nlat_half)
    S_src  = SpectralTransform(sg.spectrum, sg.grid;  NF=Float32, nlayers=1)
    S_tgt  = SpectralTransform(sg.spectrum, target;   NF=Float32, nlayers=1)
    return (target=target, S_src=S_src, S_tgt=S_tgt)
end

function regrid_level(vec_native, sg, R)
    fld  = SpeedyWeather.RingGrids.Field(Float32.(vec_native), sg.grid)
    spec = transform(fld, R.S_src)
    out  = transform(spec, R.S_tgt)
    nlon = get_nlon_max(R.target); nlat = get_nlat(R.target)
    return reshape(Array(out), nlon, nlat)
end

# ---------------------------------------------------------------------------
# monthly accumulator callback (native grid)
# ---------------------------------------------------------------------------
mutable struct MonthlyAccumulator <: SpeedyWeather.AbstractCallback
    σf::Vector{Float64}
    cnt::Vector{Int}
    temp::Any; uwnd::Any; vwnd::Any; rh::Any      # [point, k, month]
    tlow::Any; prec::Any; z500::Any               # [point, month]
end
MonthlyAccumulator() = MonthlyAccumulator(Float64[], zeros(Int, NMON),
    nothing, nothing, nothing, nothing, nothing, nothing, nothing)

function SpeedyWeather.initialize!(cb::MonthlyAccumulator, vars, model)
    cb.σf  = Float64.(Array(model.geometry.σ_levels_full))
    T = vars.grid.temperature
    np, nk = size(T)
    cb.cnt  = zeros(Int, NMON)
    cb.temp = zeros(Float64, np, nk, NMON)
    cb.uwnd = zeros(Float64, np, nk, NMON)
    cb.vwnd = zeros(Float64, np, nk, NMON)
    cb.rh   = zeros(Float64, np, nk, NMON)
    cb.tlow = zeros(Float64, np, NMON)
    cb.prec = zeros(Float64, np, NMON)
    cb.z500 = zeros(Float64, np, NMON)
    return nothing
end

function SpeedyWeather.callback!(cb::MonthlyAccumulator, vars, model)
    m   = Dates.month(vars.prognostic.clock.time)
    atm = model.atmosphere
    σf  = cb.σf
    nk  = length(σf)
    T   = Array(vars.grid.temperature)
    q   = Array(vars.grid.humidity)
    u   = Array(vars.grid.u)
    v   = Array(vars.grid.v)
    Φ   = Array(vars.grid.geopotential)
    ps  = exp.(Array(vars.grid.pressure))
    np  = length(ps)
    cb.cnt[m] += 1
    @inbounds for k in 1:nk, ij in 1:np
        p_k  = σf[k]*ps[ij]
        qsat = SpeedyWeather.saturation_humidity(T[ij, k], p_k, atm)
        cb.temp[ij, k, m] += T[ij, k]
        cb.uwnd[ij, k, m] += u[ij, k]
        cb.vwnd[ij, k, m] += v[ij, k]
        cb.rh[ij, k, m]   += qsat > 0 ? 100.0*q[ij, k]/qsat : 0.0
    end
    @inbounds for ij in 1:np
        cb.tlow[ij, m] += T[ij, nk]
        # Z500 by log-p interpolation between full levels (extrapolate below the
        # lowest full level where the surface is above 500 hPa)
        k = nk - 1
        while k > 1 && σf[k]*ps[ij] > P500
            k -= 1
        end
        lp1 = log(σf[k]*ps[ij]); lp2 = log(σf[k+1]*ps[ij])
        w   = (log(P500) - lp1)/(lp2 - lp1)
        cb.z500[ij, m] += (Φ[ij, k] + w*(Φ[ij, k+1] - Φ[ij, k]))/9.81
    end
    pr = Array(vars.parameterizations.rain_rate) .+ Array(vars.parameterizations.snow_rate)
    cb.prec[:, m] .+= pr .* (1000.0*86400.0)      # m/s -> mm/day
    return nothing
end
SpeedyWeather.finalize!(::MonthlyAccumulator, args...) = nothing

# ---------------------------------------------------------------------------
# build + run
# ---------------------------------------------------------------------------
spectral_grid = SpectralGrid(; trunc=TRUNC, nlayers=NLAYERS)
model = PrimitiveWetModel(spectral_grid;
    ocean   = SeasonalOceanClimatology(spectral_grid),
    sea_ice = PrescribedSeaIce(spectral_grid),
)
simulation = initialize!(model)
Δt = model.time_stepping.Δt_sec
println("components: ", join(["$(f)=$(nameof(typeof(getproperty(model, f))))" for f in
    (:orography, :land_sea_mask, :ocean, :land, :sea_ice, :convection,
     :large_scale_condensation, :shortwave_radiation, :longwave_radiation,
     :vertical_diffusion, :time_stepping)], ", "))
println("dt = ", Δt, " s")

run!(simulation, period=Day(SPINUP), output=false)
acc = MonthlyAccumulator()
add!(model, :monthly => acc)
t0 = time()
run!(simulation, period=Day(365*MEANYR), output=false)
wall = time() - t0
sypd = MEANYR/wall*86400
@printf("accumulation: %d yr in %.1f s -> %.0f SYPD (dt=%.0f s)\n", MEANYR, wall, sypd, Δt)

# ---------------------------------------------------------------------------
# reduce: regrid, zonal means, write NetCDF in the aeros monthly layout
# ---------------------------------------------------------------------------
R    = make_regridder(spectral_grid)
latd = Array(get_latd(R.target))            # 90 -> -90
nlat = length(latd)
nlon = get_nlon_max(R.target)
lond = collect(0:nlon-1) .* (360.0/nlon)
σf   = acc.σf
nk   = NLAYERS

t_o  = fill(NaN, nlat, nk, NMON); u_o = copy(t_o); v_o = copy(t_o); rh_o = copy(t_o)
pr_o = fill(NaN, nlat, NMON)
t2m_o = fill(NaN, nlon, nlat, NMON); pr2d_o = copy(t2m_o); z500_o = copy(t2m_o)
for m in 1:NMON
    n = acc.cnt[m]; n > 0 || continue
    for k in 1:nk
        t_o[:, k, m]  = vec(mean(regrid_level(acc.temp[:, k, m] ./ n, spectral_grid, R); dims=1))
        u_o[:, k, m]  = vec(mean(regrid_level(acc.uwnd[:, k, m] ./ n, spectral_grid, R); dims=1))
        v_o[:, k, m]  = vec(mean(regrid_level(acc.vwnd[:, k, m] ./ n, spectral_grid, R); dims=1))
        rh_o[:, k, m] = vec(mean(regrid_level(acc.rh[:, k, m]   ./ n, spectral_grid, R); dims=1))
    end
    P = regrid_level(acc.prec[:, m] ./ n, spectral_grid, R)
    pr2d_o[:, :, m] = P
    pr_o[:, m]      = vec(mean(P; dims=1))
    t2m_o[:, :, m]  = regrid_level(acc.tlow[:, m] ./ n, spectral_grid, R)
    z500_o[:, :, m] = regrid_level(acc.z500[:, m] ./ n, spectral_grid, R)
end

NCPATH = joinpath(OUTDIR, "sw_realearth_T$(TRUNC)L$(NLAYERS).nc")
NCDataset(NCPATH, "c") do ds
    defDim(ds, "lat", nlat); defDim(ds, "lon", nlon); defDim(ds, "lev", nk); defDim(ds, "month", NMON)
    v = defVar(ds, "lat", Float64, ("lat",)); v[:] = latd; v.attrib["units"] = "degrees_north"
    v = defVar(ds, "lon", Float64, ("lon",)); v[:] = lond; v.attrib["units"] = "degrees_east"
    v = defVar(ds, "lev", Float64, ("lev",)); v[:] = σf;   v.attrib["units"] = "sigma"
    v = defVar(ds, "month", Float64, ("month",)); v[:] = collect(1.0:NMON)
    function put(name, A, dims, units, ln)
        var = defVar(ds, name, Float32, dims); var[:] = Float32.(A)
        var.attrib["units"] = units; var.attrib["long_name"] = ln
    end
    put("t",  t_o,  ("lat","lev","month"), "K",   "monthly zonal-mean temperature")
    put("u",  u_o,  ("lat","lev","month"), "m/s", "monthly zonal-mean zonal wind")
    put("v",  v_o,  ("lat","lev","month"), "m/s", "monthly zonal-mean meridional wind")
    put("rh", rh_o, ("lat","lev","month"), "%",   "monthly zonal-mean relative humidity")
    put("precip",   pr_o,   ("lat","month"),       "mm/day", "monthly zonal-mean precipitation")
    put("t2m",      t2m_o,  ("lon","lat","month"), "K",      "monthly lowest-level air temperature")
    put("precip2d", pr2d_o, ("lon","lat","month"), "mm/day", "monthly precipitation")
    put("z500",     z500_o, ("lon","lat","month"), "m",      "monthly 500 hPa geopotential height")
    ds.attrib["model"] = "SpeedyWeather.jl $(pkgversion(SpeedyWeather)) PrimitiveWetModel, real Earth, prescribed seasonal SST"
    ds.attrib["trunc"] = TRUNC; ds.attrib["nlayers"] = NLAYERS
    ds.attrib["dt_sec"] = Δt; ds.attrib["spinup_days"] = SPINUP; ds.attrib["mean_years"] = MEANYR
    ds.attrib["sypd"] = sypd; ds.attrib["threads"] = Threads.nthreads()
end
println("wrote ", NCPATH)
