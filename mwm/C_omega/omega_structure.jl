# omega_structure.jl -- SpeedyWeather.jl T21 moist aquaplanet: extract the
# zonal-mean, time-mean overturning structure (omega, latent/diabatic heating,
# relative humidity) for an apples-to-apples comparison against aeros.
#
# Standard PrimitiveWetModel aquaplanet at trunc=21, L8:
#   convection = SimplifiedBettsMiller (Frierson 2007, default)
#   large-scale condensation (default ImplicitCondensation)
#   radiation (default one-band SW/LW), surface fluxes, vertical diffusion (default)
#   ocean = AquaPlanet  (SST = (Te-Tp) cos^2(lat) + Tp, Te=302 K, Tp=273 K, zonally symmetric)
#   land_sea_mask = AquaPlanetMask (all sea, no land)
#
# Diagnostics (zonal-mean, time-mean on latitude x sigma):
#   omega [hPa/day], >0 = subsidence  (from SpeedyWeather sigma-velocity w = a*sigmadot)
#   net diabatic (physics) heating [K/day]  (physics re-run each step; latent-dominated in the tropics)
#   relative humidity [%]  (q / q_sat(T, sigma*ps))
#   convective + large-scale precipitation [mm/day] vs latitude (exact column latent heating)
#
# Usage: julia --project=. omega_structure.jl [spinup_days] [mean_days]

using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using NCDatasets
using Statistics
using Printf

const TRUNC   = 21
const NLAYERS = 8
const SPINUP  = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 300
const MEANDAY = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 150
const HERE    = @__DIR__
const OUTDIR  = joinpath(HERE, "output")
mkpath(OUTDIR)

@info "omega_structure config" trunc=TRUNC nlayers=NLAYERS spinup_days=SPINUP mean_days=MEANDAY threads=Threads.nthreads()

# ---------------------------------------------------------------------------
# accumulator callback: native-grid, per-step time means
# ---------------------------------------------------------------------------
mutable struct OverturnAccumulator <: SpeedyWeather.AbstractCallback
    n::Int
    a::Float64                    # planet radius [m]
    σf::Vector{Float64}           # full sigma levels
    omega::Any                    # [point, k]  Pa/s  (per-step, >0 subsidence)
    rh::Any                       # [point, k]  fraction
    temp::Any                     # [point, k]  K
    humid::Any                    # [point, k]  kg/kg
    heat::Any                     # [point, k]  K/s   (net physics heating)
    rain_c::Any                   # [point]     m/s   convective
    rain_l::Any                   # [point]     m/s   large-scale
    uwnd::Any                     # [point, k]  m/s   zonal wind
    vwnd::Any                     # [point, k]  m/s   meridional wind
    uvp::Any                      # [point, k]  m2/s2 u*v product (for eddy momentum flux)
    # eddy zonal-wavenumber spectrum (regrid to full Gaussian, latitude DFT,
    # sampled every 10 steps) -- apples-to-apples with aeros's accum_espec.
    ke_spec::Any                  # (mmax) eddy KE by zonal wavenumber (accumulated)
    mf_spec::Any                  # (mmax) [u*v*] co-spectrum by zonal wavenumber
    nspec::Int                    # spectrum samples
    R::Any                        # regridder to full Gaussian grid
    ctab::Any                     # (mmax, nlon) DFT cos table
    stab::Any                     # (mmax, nlon) DFT -sin table
    mmax::Int
    cslat::Any                    # (nlat) cos(latitude) area weight
end
OverturnAccumulator() = OverturnAccumulator(0, 0.0, Float64[], nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, 0, nothing, nothing, nothing, 0, nothing)

function SpeedyWeather.initialize!(cb::OverturnAccumulator, vars, model)
    cb.n = 0
    cb.a = Float64(model.planet.radius)
    cb.σf = Float64.(Array(model.geometry.σ_levels_full))
    T = vars.grid.temperature
    cb.omega  = zeros(Float64, size(T))
    cb.rh     = zeros(Float64, size(T))
    cb.temp   = zeros(Float64, size(T))
    cb.humid  = zeros(Float64, size(T))
    cb.heat   = zeros(Float64, size(T))
    cb.rain_c = zeros(Float64, size(vars.parameterizations.rain_rate_convection))
    cb.rain_l = zeros(Float64, size(vars.parameterizations.rain_rate_large_scale))
    cb.uwnd = zeros(Float64, size(T))
    cb.vwnd = zeros(Float64, size(T))
    cb.uvp  = zeros(Float64, size(T))
    # eddy wavenumber spectrum setup
    sg = model.spectral_grid
    cb.R = make_regridder(sg)
    cb.mmax = Int(sg.trunc)
    nlon = get_nlon_max(cb.R.target); nlat = get_nlat(cb.R.target)
    latd = Array(get_latd(cb.R.target))
    cb.cslat = cosd.(latd)
    cb.ctab = zeros(Float64, cb.mmax, nlon)
    cb.stab = zeros(Float64, cb.mmax, nlon)
    for i in 1:nlon, m in 1:cb.mmax
        ang = 2π*m*(i-1)/nlon
        cb.ctab[m,i] = cos(ang); cb.stab[m,i] = -sin(ang)
    end
    cb.ke_spec = zeros(Float64, cb.mmax)
    cb.mf_spec = zeros(Float64, cb.mmax)
    cb.nspec = 0
    return nothing
end

function SpeedyWeather.callback!(cb::OverturnAccumulator, vars, model)
    cb.n += 1
    atm = model.atmosphere
    a   = cb.a
    σf  = cb.σf
    nlayers = length(σf)

    T  = Array(vars.grid.temperature)          # [point, k]  K
    q  = Array(vars.grid.humidity)             # [point, k]  kg/kg
    w  = Array(vars.dynamics.w)                # [point, k]  a*sigmadot at half level k+1/2
    ps = exp.(Array(vars.grid.pressure))       # [point]     surface pressure [Pa] (grid.pressure = ln ps)
    npoints = length(ps)

    # --- omega at full levels: omega = ps * sigmadot, sigmadot = w/a ---
    # half-level sigmadot: shalf[:,1]=0 (top), shalf[:,k+1]=w[:,k]/a (k=1..nlayers, w[:,nlayers]=0 surface)
    @inbounds for k in 1:nlayers
        for ij in 1:npoints
            sdot_above = k == 1 ? 0.0 : w[ij, k-1] / a
            sdot_below = w[ij, k] / a                 # w[:,nlayers]=0 gives surface sigmadot=0
            omega = ps[ij] * 0.5 * (sdot_above + sdot_below)   # Pa/s, >0 down = subsidence
            cb.omega[ij, k] += omega
            # relative humidity
            p_k  = σf[k] * ps[ij]
            qsat = SpeedyWeather.saturation_humidity(T[ij, k], p_k, atm)
            cb.rh[ij, k]   += qsat > 0 ? q[ij, k] / qsat : 0.0
            cb.temp[ij, k] += T[ij, k]
            cb.humid[ij, k]+= q[ij, k]
        end
    end

    # --- precipitation rates (instantaneous, m/s) from the real timestep ---
    cb.rain_c .+= Array(vars.parameterizations.rain_rate_convection)
    cb.rain_l .+= Array(vars.parameterizations.rain_rate_large_scale)

    # --- winds + u*v product (for the eddy momentum flux [u*v*] = [uv]-[u][v]) ---
    ug = Array(vars.grid.u)                    # [point, k]  m/s
    vg = Array(vars.grid.v)                    # [point, k]  m/s
    cb.uwnd .+= ug
    cb.vwnd .+= vg
    cb.uvp  .+= ug .* vg

    # --- eddy wavenumber spectrum, sampled every 10 steps (regrid + latitude DFT) ---
    if cb.n % 10 == 0
        sg = model.spectral_grid
        dsig = 1.0/nlayers                      # SW L8 equal-sigma (relative mass weight)
        nlon = size(cb.ctab, 2)
        for k in 1:nlayers
            Uk = regrid_level(ug[:, k], sg, cb.R)    # [lon, lat]
            Vk = regrid_level(vg[:, k], sg, cb.R)
            nlat = size(Uk, 2)
            for j in 1:nlat
                wj  = cb.cslat[j]*dsig
                uzm = sum(@view Uk[:, j])/nlon
                vzm = sum(@view Vk[:, j])/nlon
                for m in 1:cb.mmax
                    ur = 0.0; ui = 0.0; vr = 0.0; vi = 0.0
                    @inbounds for i in 1:nlon
                        up = Uk[i, j] - uzm; vp = Vk[i, j] - vzm
                        c = cb.ctab[m, i]; s = cb.stab[m, i]
                        ur += up*c; ui += up*s; vr += vp*c; vi += vp*s
                    end
                    ur /= nlon; ui /= nlon; vr /= nlon; vi /= nlon
                    cb.ke_spec[m] += wj*(ur*ur + ui*ui + vr*vr + vi*vi)
                    cb.mf_spec[m] += wj*2.0*(ur*vr + ui*vi)
                end
            end
        end
        cb.nspec += 1
    end

    # --- net diabatic (physics) heating profile: re-run physics on current state ---
    # parameterization_tendencies! resets tendencies+diagnostics, fills tend.grid.temperature
    # (radius-scaled K/s). Safe: next real timestep resets everything; only param diagnostics
    # are overwritten with identical recomputed values. Read rain rates BEFORE this (done above).
    try
        SpeedyWeather.reset_tendencies!(vars)
        SpeedyWeather.parameterization_tendencies!(vars, model)
        dT = Array(vars.tendencies.grid.temperature) ./ a    # unscale radius -> K/s
        cb.heat .+= dT
    catch err
        @warn "physics re-run failed at step $(cb.n)" exception=err
    end
    return nothing
end

SpeedyWeather.finalize!(::OverturnAccumulator, args...) = nothing

# ---------------------------------------------------------------------------
# regrid native (octahedral) single-level field -> FullGaussianGrid [nlon,nlat]
# ---------------------------------------------------------------------------
function make_regridder(sg)
    target = FullGaussianGrid(sg.nlat_half)
    S_src  = SpectralTransform(sg.spectrum, sg.grid;  NF=Float32, nlayers=1)
    return (target=target, S_src=S_src)
end

function regrid_level(vec_native, sg, R)
    # vec_native :: Vector on the octahedral grid -> Field -> spectral -> full Gaussian matrix
    fld  = SpeedyWeather.RingGrids.Field(Float32.(vec_native), sg.grid)
    spec = transform(fld, R.S_src)
    S_tgt = SpectralTransform(spec.spectrum, R.target; NF=Float32, nlayers=1)
    out  = transform(spec, S_tgt)
    nlon = get_nlon_max(R.target); nlat = get_nlat(R.target)
    return reshape(Array(out), nlon, nlat)   # [lon, lat], lat 90 -> -90
end

# zonal mean of a native single-level field -> Vector over latitude
function zonal_mean_level(vec_native, sg, R)
    M = regrid_level(vec_native, sg, R)       # [lon,lat]
    return vec(mean(M; dims=1))               # [lat]
end

# ---------------------------------------------------------------------------
# build + run
# ---------------------------------------------------------------------------
spectral_grid = SpectralGrid(; trunc=TRUNC, nlayers=NLAYERS)
model = PrimitiveWetModel(spectral_grid;
    ocean         = AquaPlanet(spectral_grid),
    land_sea_mask = AquaPlanetMask(spectral_grid),
)
simulation = initialize!(model)

@info "convection / condensation" convection=typeof(model.convection) condensation=typeof(model.large_scale_condensation)
@info "SST" ocean=typeof(model.ocean)

@info "spinup" days=SPINUP dt_sec=model.time_stepping.Δt_sec
run!(simulation, period=Day(SPINUP), output=false)

acc = OverturnAccumulator()
add!(model, :overturn => acc)
@info "averaging" days=MEANDAY
run!(simulation, period=Day(MEANDAY), output=false)
acc.n > 0 || error("no steps accumulated")
@info "accumulated" nsteps=acc.n

# ---------------------------------------------------------------------------
# reduce to zonal-mean, time-mean on (lat, sigma)
# ---------------------------------------------------------------------------
n = acc.n
R = make_regridder(spectral_grid)
latd = Array(get_latd(R.target))           # 90 -> -90
nlat = length(latd)
σf   = acc.σf

# per-point time means
omega_m = acc.omega ./ n        # Pa/s
rh_m    = acc.rh    ./ n        # fraction
temp_m  = acc.temp  ./ n        # K
humid_m = acc.humid ./ n        # kg/kg
heat_m  = acc.heat  ./ n        # K/s
rc_m    = acc.rain_c ./ n       # m/s
rl_m    = acc.rain_l ./ n       # m/s
u_m     = acc.uwnd  ./ n        # m/s
v_m     = acc.vwnd  ./ n        # m/s
uv_m    = acc.uvp   ./ n        # m2/s2  time-mean of point product u*v

# zonal means on (lat, sigma)
OMEGA = zeros(nlat, NLAYERS)    # hPa/day
RH    = zeros(nlat, NLAYERS)    # %
TEMP  = zeros(nlat, NLAYERS)    # K
HEAT  = zeros(nlat, NLAYERS)    # K/day
UWND  = zeros(nlat, NLAYERS)    # m/s   [ubar]  (the jet)
VWND  = zeros(nlat, NLAYERS)    # m/s   [vbar]  (the MMC)
UVPR  = zeros(nlat, NLAYERS)    # m2/s2 eddy momentum flux [u*v*] = [uv]-[u][v]
for k in 1:NLAYERS
    OMEGA[:, k] = zonal_mean_level(omega_m[:, k], spectral_grid, R) .* 864.0     # Pa/s -> hPa/day
    RH[:, k]    = zonal_mean_level(rh_m[:, k],    spectral_grid, R) .* 100.0     # -> %
    TEMP[:, k]  = zonal_mean_level(temp_m[:, k],  spectral_grid, R)
    HEAT[:, k]  = zonal_mean_level(heat_m[:, k],  spectral_grid, R) .* 86400.0   # K/s -> K/day
    ubar_k      = zonal_mean_level(u_m[:, k],  spectral_grid, R)                 # [ubar]
    vbar_k      = zonal_mean_level(v_m[:, k],  spectral_grid, R)                 # [vbar]
    uvbar_k     = zonal_mean_level(uv_m[:, k], spectral_grid, R)                 # [ubar*vbar overline] = [<uv>]
    UWND[:, k]  = ubar_k
    VWND[:, k]  = vbar_k
    UVPR[:, k]  = uvbar_k .- ubar_k .* vbar_k                                    # stationary+transient eddy flux
end
PRECIP_C = zonal_mean_level(rc_m, spectral_grid, R) .* (1000.0 * 86400.0)        # mm/day
PRECIP_L = zonal_mean_level(rl_m, spectral_grid, R) .* (1000.0 * 86400.0)        # mm/day

# ---------------------------------------------------------------------------
# save NetCDF
# ---------------------------------------------------------------------------
jperm = sortperm(latd)          # ascending -90 -> 90
lat_asc = latd[jperm]
NCPATH = joinpath(OUTDIR, "speedy_omega_T$(TRUNC)L$(NLAYERS).nc")
NCDataset(NCPATH, "c") do ds
    defDim(ds, "lat", nlat); defDim(ds, "sigma", NLAYERS)
    v = defVar(ds, "lat", Float64, ("lat",));   v[:] = lat_asc; v.attrib["units"]="degrees_north"
    v = defVar(ds, "sigma", Float64, ("sigma",)); v[:] = σf
    function put2(name, M, units, ln)
        var = defVar(ds, name, Float32, ("lat","sigma")); var[:,:] = Float32.(M[jperm, :])
        var.attrib["units"]=units; var.attrib["long_name"]=ln
    end
    put2("omega", OMEGA, "hPa/day", "zonal-mean vertical velocity (>0 subsidence)")
    put2("rh",    RH,    "%",       "zonal-mean relative humidity")
    put2("temp",  TEMP,  "K",       "zonal-mean temperature")
    put2("heat",  HEAT,  "K/day",   "zonal-mean net diabatic (physics) heating")
    put2("uwind", UWND,  "m/s",     "time-mean zonal-mean zonal wind (jet)")
    put2("vwind", VWND,  "m/s",     "time-mean zonal-mean meridional wind (MMC)")
    put2("uvpr",  UVPR,  "m2/s2",   "time-mean zonal-mean eddy momentum flux [u*v*]=[uv]-[u][v]")
    v = defVar(ds, "precip_convective", Float32, ("lat",)); v[:] = Float32.(PRECIP_C[jperm]); v.attrib["units"]="mm/day"
    v = defVar(ds, "precip_largescale", Float32, ("lat",)); v[:] = Float32.(PRECIP_L[jperm]); v.attrib["units"]="mm/day"
    ds.attrib["model"]="SpeedyWeather.jl PrimitiveWetModel aquaplanet"
    ds.attrib["trunc"]=TRUNC; ds.attrib["nlayers"]=NLAYERS
    ds.attrib["spinup_days"]=SPINUP; ds.attrib["mean_days"]=MEANDAY; ds.attrib["nsteps"]=n
    ds.attrib["SST"]="AquaPlanet cos^2 lat, Te=302K Tp=273K"
end
@info "wrote NetCDF" path=NCPATH

# ---------------------------------------------------------------------------
# console summary + write a compact results text used to build the markdown
# ---------------------------------------------------------------------------
abslat = abs.(latd)

# a representative free-troposphere level index (~500 hPa, sigma~0.5)
k500 = argmin(abs.(σf .- 0.5))
kmid = argmin(abs.(σf .- 0.44))   # ~440 hPa

# column-max ascent per latitude (most negative omega over levels 2..nlayers-1)
col_asc = [minimum(OMEGA[j, 2:NLAYERS-1]) for j in 1:nlat]   # most negative
# TROPICAL (ITCZ) peak ascent: restrict to |lat| <= 25 deg
trop = findall(x -> abs(x) <= 25, latd)
jasc = trop[argmin(col_asc[trop])]
peak_asc_val = col_asc[jasc]; peak_asc_lat = latd[jasc]
kasc = argmin(OMEGA[jasc, :])
# global-max ascent (may be extratropical/eddy)
jglob = argmin(col_asc); glob_asc_val = col_asc[jglob]; glob_asc_lat = latd[jglob]

# subtropical descent: max column-mean subsidence in 15-35 band
band = findall(x -> 15 <= abs(x) <= 35, latd)
col_desc = [maximum(OMEGA[j, 2:NLAYERS-1]) for j in 1:nlat]   # most positive per lat
jdesc = band[argmax(col_desc[band])]
sub_desc_val = col_desc[jdesc]; sub_desc_lat = latd[jdesc]

# free-trop subtropical RH (15-35 band, ~500 hPa)
rh_sub = mean(RH[band, k500])

# precip structure
jprc = argmax(PRECIP_C .+ PRECIP_L)
tot_precip = PRECIP_C .+ PRECIP_L

open(joinpath(OUTDIR, "summary.txt"), "w") do io
    @printf(io, "nsteps=%d  k500=%d(sigma=%.3f)\n", n, k500, σf[k500])
    @printf(io, "PEAK_ASCENT(tropical|lat|<=25) omega=%.2f hPa/day at lat=%.1f (sigma=%.3f)\n",
            peak_asc_val, peak_asc_lat, σf[kasc])
    @printf(io, "GLOBAL_MAX_ASCENT omega=%.2f hPa/day at lat=%.1f (may be extratropical)\n",
            glob_asc_val, glob_asc_lat)
    @printf(io, "SUBTROP_DESCENT omega=%.2f hPa/day at lat=%.1f (15-35 band)\n", sub_desc_val, sub_desc_lat)
    @printf(io, "SUBTROP_FREETROP_RH=%.1f %% (15-35 band, sigma=%.3f)\n", rh_sub, σf[k500])
    @printf(io, "PRECIP peak (conv+ls)=%.2f mm/day at lat=%.1f\n", tot_precip[jprc], latd[jprc])
    println(io, "\n# per-latitude column-max ascent and column-max descent (hPa/day), precip (mm/day)")
    println(io, "# lat  col_ascent  col_descent  precip_conv  precip_ls")
    for j in jperm
        @printf(io, "%6.1f  %8.2f  %8.2f  %8.3f  %8.3f\n", latd[j], col_asc[j], col_desc[j], PRECIP_C[j], PRECIP_L[j])
    end
    println(io, "\n# OMEGA(lat,sigma) hPa/day  rows=lat(90..-90) cols=sigma")
    print(io, "  lat\\sig"); for k in 1:NLAYERS; @printf(io, " %7.3f", σf[k]); end; println(io)
    for j in 1:nlat
        @printf(io, "%7.1f", latd[j]); for k in 1:NLAYERS; @printf(io, " %7.2f", OMEGA[j,k]); end; println(io)
    end
    println(io, "\n# RH(lat,sigma) %")
    print(io, "  lat\\sig"); for k in 1:NLAYERS; @printf(io, " %7.3f", σf[k]); end; println(io)
    for j in 1:nlat
        @printf(io, "%7.1f", latd[j]); for k in 1:NLAYERS; @printf(io, " %7.1f", RH[j,k]); end; println(io)
    end
    println(io, "\n# HEAT(lat,sigma) K/day (net physics)")
    print(io, "  lat\\sig"); for k in 1:NLAYERS; @printf(io, " %7.3f", σf[k]); end; println(io)
    for j in 1:nlat
        @printf(io, "%7.1f", latd[j]); for k in 1:NLAYERS; @printf(io, " %7.2f", HEAT[j,k]); end; println(io)
    end
end

println("\n================ SUMMARY ================")
@printf("TROPICAL PEAK ASCENT : omega = %.2f hPa/day at lat %.1f  (sigma=%.3f)\n", peak_asc_val, peak_asc_lat, σf[kasc])
@printf("GLOBAL MAX ASCENT    : omega = %.2f hPa/day at lat %.1f  (extratrop eddy?)\n", glob_asc_val, glob_asc_lat)
@printf("SUBTROP DESC : omega = %.2f hPa/day at lat %.1f  (15-35 band)\n", sub_desc_val, sub_desc_lat)
@printf("SUBTROP RH  : %.1f %% (free trop, sigma=%.3f, 15-35 band)\n", rh_sub, σf[k500])
@printf("PRECIP peak : %.2f mm/day at lat %.1f\n", tot_precip[jprc], latd[jprc])

# --- eddy zonal-wavenumber spectrum (apples-to-apples with aeros) ---
if acc.nspec > 0
    ke = acc.ke_spec ./ acc.nspec
    mf = acc.mf_spec ./ acc.nspec
    ket = sum(ke); mft = sum(mf)
    open(joinpath(OUTDIR, "eddy_spectrum.txt"), "w") do io
        @printf(io, "# SpeedyWeather eddy zonal-wavenumber spectrum (mass+area wtd, %d samples)\n", acc.nspec)
        @printf(io, "# total eddy KE %.3e   total [u*v*] %.3e\n", ket, mft)
        @printf(io, "#  m    KE(m)      KE%%     [u*v*](m)   flux%%\n")
        for m in 1:acc.mmax
            @printf(io, "%4d %11.3e %7.1f %12.3e %7.1f\n", m, ke[m], 100*ke[m]/max(ket,eps()),
                    mf[m], 100*mf[m]/max(abs(mft),eps()))
        end
    end
    println("\n=== SpeedyWeather eddy KE spectrum by zonal wavenumber m ===")
    @printf("total eddy KE %.3e   total [u*v*] %.3e\n", ket, mft)
    @printf("  m    KE%%    flux%%\n")
    for m in 1:min(acc.mmax, 12)
        @printf("%4d %7.1f %7.1f\n", m, 100*ke[m]/max(ket,eps()), 100*mf[m]/max(abs(mft),eps()))
    end
end
println("wrote ", joinpath(OUTDIR,"summary.txt"), " and ", NCPATH)
@info "DONE"
