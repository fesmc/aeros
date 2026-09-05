# compare_realearth.jl -- aeros vs SpeedyWeather vs ERA5 at real-Earth boundary
# conditions (docs/m2_handoff.md 2026-09-05, "NEXT" item 2).
#
# Inputs (monthly climatologies in the rce_long `monthly_out` layout):
#   aeros : output/aeros_realearth_T21L12*.nc   (rce_long, mwm/D_realearth/aeros_*.nml)
#   SW    : output/sw_realearth_T21L8.nc        (sw_realearth.jl)
#   ERA5  : ~/data/era5 1991-2020 monthly climatologies (pressure + single levels)
#
# Figures (output/):
#   realearth_zonal_annual.png   zonal-mean T / u / RH cross-sections, 3 models
#   realearth_precip.png         zonal-mean precipitation, annual / DJF / JJA
#   realearth_z500star_djf.png   DJF 500 hPa stationary wave Z* = Z - [Z], NH
#   realearth_t2m_amp.png        t2m seasonal amplitude (Jul-Jan) maps
# and a metrics table, printed and written to output/realearth_metrics.md.
#
# Usage: julia compare_realearth.jl [aeros.nc] [sw.nc]

using NCDatasets, CairoMakie, Statistics, Printf

const HERE   = @__DIR__
const OUT    = joinpath(HERE, "output")
const ERA5_PL = expanduser("~/data/era5/monthly-pressure-levels")
const ERA5_SL = expanduser("~/data/era5/monthly-single-levels")
plc(v) = joinpath(ERA5_PL, "era5_monthly-pressure-levels_$(v)_1991-2020_clim.nc")
slc(v) = joinpath(ERA5_SL, "era5_monthly-single-levels_$(v)_1991-2020_clim.nc")

aeros_nc = length(ARGS) >= 1 ? ARGS[1] : joinpath(OUT, "aeros_realearth_T21L12_dt2700.nc")
sw_nc    = length(ARGS) >= 2 ? ARGS[2] : joinpath(OUT, "sw_realearth_T21L8.nc")

# --- helpers ----------------------------------------------------------------
function interp1(x, y, xq)
    p = sortperm(x); xs = x[p]; ys = y[p]
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    k = searchsortedlast(xs, xq); k == length(xs) && return ys[end]
    t = (xq - xs[k]) / (xs[k+1] - xs[k])
    return ys[k]*(1-t) + ys[k+1]*t
end
function regrid2d(lonS, latS, F, lonT, latT)   # F(lonS,latS) -> (lonT,latT), lon periodic via wrap
    lonSw = vcat(lonS, lonS[1] + 360.0); Fw = vcat(F, F[1:1, :])
    tmp = [interp1(latS, Fw[i, :], la) for i in eachindex(lonSw), la in latT]
    return [interp1(lonSw, tmp[:, jl], lo) for lo in lonT, jl in eachindex(latT)]
end
regrid_latp(latS, pS, F, latT, pT) = begin       # F(latS,pS) -> (latT,pT)
    tmp = [interp1(latS, F[:, jp], la) for la in latT, jp in eachindex(pS)]
    [interp1(pS, tmp[jl, :], p) for jl in eachindex(latT), p in pT]
end
gmean(F, lat) = begin                            # area-weighted mean of F(lon,lat) or F(lat)
    w = cosd.(lat); w ./= sum(w)
    ndims(F) == 1 ? sum(F .* w) : sum(vec(mean(F; dims=1)) .* w)
end
season(A, ms) = dropdims(mean(A[:, :, ms]; dims=3); dims=3)    # (x,y,month) -> (x,y)
const DJF = [12, 1, 2]; const JJA = [6, 7, 8]; const ANN = 1:12

# --- model readers (rce_long monthly layout; NCDatasets = dim1-first) --------
function read_model(path)
    ds = NCDataset(path)
    m = (lat=Array(ds["lat"][:]), lon=Array(ds["lon"][:]), sig=Array(ds["lev"][:]),
         t=Array(ds["t"][:,:,:]), u=Array(ds["u"][:,:,:]), rh=Array(ds["rh"][:,:,:]),
         precip=Array(ds["precip"][:,:]), t2m=Array(ds["t2m"][:,:,:]),
         precip2d=Array(ds["precip2d"][:,:,:]), z500=Array(ds["z500"][:,:,:]),
         sypd = haskey(ds.attrib, "sypd") ? ds.attrib["sypd"] : NaN,
         name = haskey(ds.attrib, "model") ? split(ds.attrib["model"], " ")[1] : "aeros")
    close(ds); return m
end
A = read_model(aeros_nc); A = merge(A, (name="aeros",))
# aeros throughput from the run's `time -p` log next to its namelist (rce_long does not stamp it)
let base = splitext(basename(aeros_nc))[1], tf = joinpath(HERE, "logs", base * ".time"), nf = joinpath(HERE, base * ".nml")
    if isfile(tf) && isfile(nf)
        real = parse(Float64, split(filter(l -> startswith(l, "real"), readlines(tf))[1])[2])
        getkey(k) = parse(Float64, strip(split(filter(l -> occursin(Regex("^\\s*$k\\s*="), l), readlines(nf))[1], "=")[2]))
        yrs = getkey("nstep") * getkey("dt") / (365.0 * 86400.0)
        global A = merge(A, (sypd = yrs / real * 86400.0,))
    end
end
S = read_model(sw_nc);    S = merge(S, (name="SpeedyWeather",))
# missing-value months (rce_long writes -9999 for empty months)
for M in (A, S)
    for f in (:t, :u, :rh, :precip, :t2m, :precip2d, :z500)
        X = getfield(M, f); X[X .<= -9000] .= NaN
    end
end
lat = A.lat; lon = A.lon
@assert length(lat) == length(S.lat) && length(lon) == length(S.lon) "both models must share the T21 grid"

# --- ERA5 -------------------------------------------------------------------
function era5_pl(var)
    ds = NCDataset(plc(var))
    X = Array(ds[var][:,:,:,:]); latE = Array(ds["latitude"][:]); plev = Array(ds["pressure_level"][:])
    lonE = Array(ds["longitude"][:]); close(ds)
    return lonE, latE, plev, X                  # (lon, lat, plev, month)
end
function era5_sl(var)
    ds = NCDataset(slc(var)); X = Array(ds[var][:,:,:]); lonE = Array(ds["longitude"][:]); latE = Array(ds["latitude"][:]); close(ds)
    return lonE, latE, X                        # (lon, lat, month)
end
lonE, latE, plev, tE = era5_pl("t");  _, _, _, uE = era5_pl("u");  _, _, _, rE = era5_pl("r");  _, _, _, zE = era5_pl("z")
_, _, tpE  = era5_sl("tp");  _, _, t2mE = era5_sl("t2m");  _, _, lsmE = era5_sl("lsm")
tpE = tpE .* 1000.0                             # m/day -> mm/day
k500 = argmin(abs.(plev .- 500.0))
zE500 = zE[:, :, k500, :] ./ 9.80665            # geopotential -> height [m]
zonal(X) = dropdims(mean(X; dims=1); dims=1)    # (lon,lat,...) -> (lat,...)

# model sigma -> pressure with ps ~ 1000 hPa (cross-sections only)
pA = A.sig .* 1000.0; pS = S.sig .* 1000.0
pgrid = [925.0, 850, 775, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70]
latT = sort(lat)                                 # -90..90 for plotting

# --- 1. zonal-mean annual cross-sections ------------------------------------
function field_on(M, f, p, ms)                   # model (lat,lev,month) -> (latT,pgrid)
    X = dropdims(mean(getfield(M, f)[:, :, ms]; dims=3); dims=3)
    regrid_latp(M.lat, p, X, latT, pgrid)
end
function era_on(X, ms)                            # ERA5 (lon,lat,plev,month) -> (latT,pgrid)
    Z = zonal(dropdims(mean(X[:, :, :, ms]; dims=4); dims=4))   # (lat,plev)
    regrid_latp(latE, plev, Z, latT, pgrid)
end
fig = Figure(size=(1500, 1100))
rows = [("T [K]", :t, tE, 180:10:300, :thermal), ("u [m/s]", :u, uE, -30:5:60, :RdBu),
        ("RH [%]", :rh, rE, 0:10:100, :Blues)]
for (r, (lab, f, XE, lv, cm)) in enumerate(rows)
    panels = [(A.name, field_on(A, f, pA, ANN)), (S.name, field_on(S, f, pS, ANN)), ("ERA5", era_on(XE, ANN))]
    for (c, (nm, Z)) in enumerate(panels)
        ax = Axis(fig[r, c], title="$nm  $lab (annual)", xlabel="lat", ylabel="p [hPa]", yreversed=true)
        cf = contourf!(ax, latT, pgrid, Z; levels=lv, colormap=cm, extendlow=:auto, extendhigh=:auto)
        c == 3 && Colorbar(fig[r, 4], cf)
    end
end
save(joinpath(OUT, "realearth_zonal_annual.png"), fig); println("wrote realearth_zonal_annual.png")

# --- 2. zonal-mean precipitation --------------------------------------------
fig = Figure(size=(1300, 420))
for (c, (nm, ms)) in enumerate([("annual", ANN), ("DJF", DJF), ("JJA", JJA)])
    ax = Axis(fig[1, c], title="zonal-mean precipitation, $nm", xlabel="lat", ylabel="mm/day")
    pa = [interp1(A.lat, vec(mean(A.precip[:, ms]; dims=2)), la) for la in latT]
    ps_ = [interp1(S.lat, vec(mean(S.precip[:, ms]; dims=2)), la) for la in latT]
    pe = [interp1(latE, zonal(season(tpE, ms)), la) for la in latT]
    lines!(ax, latT, pa, label=A.name); lines!(ax, latT, ps_, label=S.name); lines!(ax, latT, pe, label="ERA5", color=:black)
    c == 1 && axislegend(ax)
end
save(joinpath(OUT, "realearth_precip.png"), fig); println("wrote realearth_precip.png")

# --- 3. DJF Z500 stationary wave, NH ----------------------------------------
zstar(Z) = Z .- mean(Z; dims=1)                   # (lon,lat) minus zonal mean
zA = zstar(season(A.z500, DJF)); zS = zstar(season(S.z500, DJF))
zEm = zstar(regrid2d(lonE, latE, season(zE500, DJF), lon, lat))
nh = findall(l -> 20 <= l <= 80, lat)
rmsw(Z) = begin w = cosd.(lat[nh]); sqrt(sum(vec(mean(Z[:, nh].^2; dims=1)) .* w)/sum(w)) end
fig = Figure(size=(1500, 420))
for (c, (nm, Z)) in enumerate([(A.name, zA), (S.name, zS), ("ERA5", zEm)])
    ax = Axis(fig[1, c], title=@sprintf("%s  DJF Z500* [m]  (NH 20-80 RMS %.0f m)", nm, rmsw(Z)), xlabel="lon", ylabel="lat")
    ord = sortperm(lat)
    cf = contourf!(ax, lon, lat[ord], Z[:, ord]; levels=-200:25:200, colormap=:RdBu, extendlow=:auto, extendhigh=:auto)
    ylims!(ax, 0, 90); c == 3 && Colorbar(fig[1, 4], cf)
end
save(joinpath(OUT, "realearth_z500star_djf.png"), fig); println("wrote realearth_z500star_djf.png")

# --- 4. t2m seasonal amplitude (Jul - Jan) ----------------------------------
ampA = A.t2m[:, :, 7] .- A.t2m[:, :, 1]; ampS = S.t2m[:, :, 7] .- S.t2m[:, :, 1]
t2mEm = cat([regrid2d(lonE, latE, t2mE[:, :, m], lon, lat) for m in 1:12]...; dims=3)
ampE = t2mEm[:, :, 7] .- t2mEm[:, :, 1]
lsmM = regrid2d(lonE, latE, dropdims(mean(lsmE; dims=3); dims=3), lon, lat); landM = lsmM .>= 0.5
fig = Figure(size=(1500, 420))
for (c, (nm, Z)) in enumerate([(A.name, ampA), (S.name, ampS), ("ERA5", ampE)])
    ax = Axis(fig[1, c], title="$nm  t2m Jul - Jan [K]", xlabel="lon", ylabel="lat")
    ord = sortperm(lat)
    cf = contourf!(ax, lon, lat[ord], Z[:, ord]; levels=-40:5:40, colormap=:RdBu, extendlow=:auto, extendhigh=:auto)
    c == 3 && Colorbar(fig[1, 4], cf)
end
save(joinpath(OUT, "realearth_t2m_amp.png"), fig); println("wrote realearth_t2m_amp.png")

# --- 5. metrics --------------------------------------------------------------
nhland = findall(l -> 40 <= l <= 70, lat)
landamp(amp) = begin                              # mean |Jul-Jan| over NH 40-70 land points
    v = Float64[]; for j in nhland, i in eachindex(lon); landM[i, j] && push!(v, abs(amp[i, j])); end; mean(v) end
tann(M) = dropdims(mean(M.t; dims=3); dims=3)     # (lat,lev)
uann(M) = dropdims(mean(M.u; dims=3); dims=3)
lidT(M) = gmean(tann(M)[:, 1], M.lat)
jetmax(M) = begin U = uann(M); i = argmax(U); (U[i], M.lat[i[1]], M.sig[i[2]]) end
tE_ann = zonal(dropdims(mean(tE; dims=4); dims=4)); uE_ann = zonal(dropdims(mean(uE; dims=4); dims=4))
k70 = argmin(abs.(plev .- 70.0)); k200 = argmin(abs.(plev .- 200.0))
rowsM = [
    ("global-mean precip [mm/day]", gmean(dropdims(mean(A.precip2d; dims=3); dims=3), A.lat), gmean(dropdims(mean(S.precip2d; dims=3); dims=3), S.lat), gmean(dropdims(mean(tpE; dims=3); dims=3), latE)),
    ("global-mean T, top level / 70 hPa [K]", lidT(A), lidT(S), gmean(tE_ann[:, k70], latE)),
    ("global-mean T, lowest level / 2 m [K]", gmean(tann(A)[:, end], A.lat), gmean(tann(S)[:, end], S.lat), gmean(dropdims(mean(t2mE; dims=3); dims=3), latE)),
    ("max annual zonal-mean u [m/s]", jetmax(A)[1], jetmax(S)[1], maximum(uE_ann)),
    ("max annual [u] at 200 hPa [m/s]", maximum(regrid_latp(A.lat, pA, uann(A), latT, [200.0])), maximum(regrid_latp(S.lat, pS, uann(S), latT, [200.0])), maximum(uE_ann[:, k200])),
    ("DJF NH Z500* RMS 20-80N [m]", rmsw(zA), rmsw(zS), rmsw(zEm)),
    ("t2m |Jul-Jan|, NH 40-70N land [K]", landamp(ampA), landamp(ampS), landamp(ampE)),
    ("throughput [SYPD, 10 threads]", A.sypd, S.sypd, NaN),
]
io = IOBuffer()
println(io, "| metric | aeros | SpeedyWeather | ERA5 |"); println(io, "|---|---|---|---|")
for (nm, a, s, e) in rowsM
    @printf(io, "| %s | %.1f | %.1f | %.1f |\n", nm, a, s, e)
end
tbl = String(take!(io)); print(tbl)
open(joinpath(OUT, "realearth_metrics.md"), "w") do f
    println(f, "# aeros vs SpeedyWeather vs ERA5, real-Earth T21\n"); println(f, "aeros: `$aeros_nc`  \nSW: `$sw_nc`\n"); print(f, tbl)
end
println("wrote realearth_metrics.md")
