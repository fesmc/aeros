# Model RCE humidity bias vs ERA5 (M2 handoff step 3).
#
# The coupled RCE overcasts (rce_long TOA clouds-on -36 W/m^2) even though the
# diagnostic cloud scheme is validated correct on ERA5's own columns
# (m2_results §21). That points the finger at a *model* moisture bias: the RCE
# equilibrates too moist (near-saturated over a deep layer -> high cloud
# fraction), so it makes more cloud than ERA5 would. This script lays the model's
# own zonal-mean RH(lat, p) against the ERA5 1991-2020 climatology to localize the
# bias in latitude and height, and compares two vehicles that share an identical
# surface boundary (prescribed SST) and differ only in the dynamics:
#
#   rot  — rotating, realistic insolation (has a meridional circulation)
#   uni  — non-rotating, uniform insolation (each column ~ a local 1D RCE)
#
# If the rot subtropics are moist where ERA5 is dry, that is the missing
# subsidence-drying (structural, T21 has weak Hadley cells). If even the uni
# columns sit near-saturated where ERA5 convective columns are ~70-80%, the
# column physics itself (convection/condensation/evaporation balance) is too
# moist — a model bug to chase.
#
# The model dumps come from drivers/rce_long.f90 (arg 2 = output NetCDF path):
#   rce_long.x logs/rce_hum_rot.nml output/rce_rh_rot.nc
#   rce_long.x logs/rce_hum_uni.nml output/rce_rh_uni.nc
#
#   julia scripts/rce_humidity_vs_era5.jl [rot.nc] [uni.nc]
#
# Writes docs/figures/rce_humidity_vs_era5.png and prints area-weighted summaries.

using NCDatasets
using CairoMakie
using Statistics
using Printf

const ERA5_PL = "/Users/alrobi001/data/era5/monthly-pressure-levels"
const ERA5_SL = "/Users/alrobi001/data/era5/monthly-single-levels"
const OUTDIR  = "docs/figures"

rot_nc  = length(ARGS) >= 1 ? ARGS[1] : "output/rce_rh_rot.nc"
uni_nc  = length(ARGS) >= 2 ? ARGS[2] : "output/rce_rh_uni.nc"
outname = length(ARGS) >= 3 ? ARGS[3] : "rce_humidity_vs_era5.png"

plc(v) = joinpath(ERA5_PL, "era5_monthly-pressure-levels_$(v)_1991-2020_clim.nc")
slc(v) = joinpath(ERA5_SL, "era5_monthly-single-levels_$(v)_1991-2020_clim.nc")

# --- readers -------------------------------------------------------------

# Model dump: return a field oriented (lat, lev) regardless of on-disk order.
function read_model(path)
    ds = NCDataset(path)
    lat = ds["lat"][:]; lev = ds["lev"][:]
    nlat = length(lat); nlev = length(lev)
    grab(v) = begin
        A = Array(ds[v][:, :])
        size(A) == (nlat, nlev) ? A : permutedims(A)   # -> (lat, lev)
    end
    m = (lat = lat, lev = lev,
         rh = grab("rh"), cf = grab("cf"), t = grab("t"),
         q = grab("q"), pfull = grab("pfull"), cover = ds["cover"][:])
    close(ds)
    return m
end

# ERA5 annual + zonal mean pressure-level field -> (lat, plev)
function era5_pl_zonal(var)
    ds = NCDataset(plc(var))
    A = Array(ds[var][:, :, :, :])                    # (lon, lat, plev, month)
    lat = ds["latitude"][:]; plev = ds["pressure_level"][:]
    close(ds)
    zm = dropdims(mean(A; dims = (1, 4)); dims = (1, 4))   # (lat, plev)
    return lat, plev, zm
end

function era5_sl(var)
    ds = NCDataset(slc(var))
    A = Array(ds[var][:, :, :])                       # (lon, lat, month)
    lat = ds["latitude"][:]
    close(ds)
    return lat, dropdims(mean(A; dims = (1, 3)); dims = (1, 3))  # zonal+annual (lat)
end

# --- helpers -------------------------------------------------------------

# 1D linear interp of y(x) at xq (x need not be ascending); clamps to ends.
function interp1(x, y, xq)
    p = sortperm(x); xs = x[p]; ys = y[p]
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    k = searchsortedlast(xs, xq)
    t = (xq - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] * (1 - t) + ys[k+1] * t
end

# Regrid ERA5 field E(latE, plevE) onto the model (lat, p) grid, where the model
# pressure axis pm[lev] is the lat-mean layer pressure. Bilinear via two passes.
function regrid_to_model(latE, plevE, E, latM, pm)
    tmp = Matrix{Float64}(undef, length(latM), length(plevE))  # interp in lat
    for (jl, la) in enumerate(latM), (jp, _) in enumerate(plevE)
        tmp[jl, jp] = interp1(latE, E[:, jp], la)
    end
    out = Matrix{Float64}(undef, length(latM), length(pm))     # interp in p
    for jl in axes(out, 1), (jp, p) in enumerate(pm)
        out[jl, jp] = interp1(plevE, tmp[jl, :], p)
    end
    return out
end

# area-weighted (cos lat) mean profile over latitude -> (lev,)
function latmean(A, lat)
    w = cos.(deg2rad.(lat)); w ./= sum(w)
    return vec(sum(A .* w; dims = 1))
end

# --- load ----------------------------------------------------------------

rot = read_model(rot_nc)
uni = read_model(uni_nc)

latE, plevE, rhE = era5_pl_zonal("r")     # RH %
_,    _,     ccE = era5_pl_zonal("cc")    # cloud fraction 0-1
latT, tccE       = era5_sl("tcc")         # total cloud cover 0-1

# model lat-mean pressure axis (ps ~ uniform in RCE, so this is well-defined)
pm  = latmean(rot.pfull, rot.lat)         # (lev,) hPa, index 1 = model top
ord = sortperm(rot.lat)                   # ascending latitude for plotting
latM = rot.lat[ord]

rhE_m  = regrid_to_model(latE, plevE, rhE, latM, pm)   # ERA5 RH on model grid
ccE_m  = regrid_to_model(latE, plevE, ccE, latM, pm)
rot_rh = rot.rh[ord, :]; rot_cf = rot.cf[ord, :]
uni_rh = uni.rh[ord, :]; uni_cf = uni.cf[ord, :]
drh    = rot_rh .- rhE_m
dcf    = rot_cf .- ccE_m

# --- text summary --------------------------------------------------------

wmean(v, lat) = (w = cos.(deg2rad.(lat)); sum(v .* w) / sum(w))
@printf("\nArea-weighted total cloud cover:  rot %.2f   uni %.2f   ERA5 %.2f\n",
        wmean(rot.cover, rot.lat), wmean(uni.cover, uni.lat), wmean(tccE, latT))
println("\nArea-weighted mean RH (%) by level (model top -> surface):")
println("   p[hPa]   rot    uni    ERA5   rot-ERA5")
rotP = latmean(rot_rh, latM); uniP = latmean(uni_rh, latM); eraP = latmean(rhE_m, latM)
for k in eachindex(pm)
    @printf("  %7.1f  %5.1f  %5.1f  %5.1f  %+6.1f\n",
            pm[k], rotP[k], uniP[k], eraP[k], rotP[k] - eraP[k])
end

# --- figure --------------------------------------------------------------

mkpath(OUTDIR)
fig = Figure(size = (1500, 1150), fontsize = 15)

# Row 1: RH cross-sections
ax = Axis(fig[1, 1]; title = "Model RH — rotating", xlabel = "latitude",
          ylabel = "pressure (hPa)", yreversed = true)
hm = heatmap!(ax, latM, pm, rot_rh; colormap = :viridis, colorrange = (0, 100))
Colorbar(fig[1, 2], hm)
ax = Axis(fig[1, 3]; title = "ERA5 RH (r)", xlabel = "latitude",
          ylabel = "pressure (hPa)", yreversed = true)
hm = heatmap!(ax, latM, pm, rhE_m; colormap = :viridis, colorrange = (0, 100))
Colorbar(fig[1, 4], hm)

# Row 2: RH difference + cloud fraction
ax = Axis(fig[2, 1]; title = "RH bias  rot − ERA5  (%)", xlabel = "latitude",
          ylabel = "pressure (hPa)", yreversed = true)
b = 40
hm = heatmap!(ax, latM, pm, drh; colormap = :balance, colorrange = (-b, b))
Colorbar(fig[2, 2], hm)
ax = Axis(fig[2, 3]; title = "Cloud fraction bias  rot − ERA5", xlabel = "latitude",
          ylabel = "pressure (hPa)", yreversed = true)
hm = heatmap!(ax, latM, pm, dcf; colormap = :balance, colorrange = (-0.5, 0.5))
Colorbar(fig[2, 4], hm)

# Row 3: cloud fraction model vs ERA5
ax = Axis(fig[3, 1]; title = "Model cloud fraction — rotating", xlabel = "latitude",
          ylabel = "pressure (hPa)", yreversed = true)
hm = heatmap!(ax, latM, pm, rot_cf; colormap = :viridis, colorrange = (0, 0.6))
Colorbar(fig[3, 2], hm)
ax = Axis(fig[3, 3]; title = "ERA5 cloud fraction (cc)", xlabel = "latitude",
          ylabel = "pressure (hPa)", yreversed = true)
hm = heatmap!(ax, latM, pm, ccE_m; colormap = :viridis, colorrange = (0, 0.6))
Colorbar(fig[3, 4], hm)

# Row 4: total cloud cover vs lat, and area-mean RH profiles
ax = Axis(fig[4, 1]; title = "Total cloud cover vs latitude",
          xlabel = "latitude", ylabel = "cover")
lines!(ax, rot.lat[ord], rot.cover[ord]; label = "rot", color = :firebrick)
lines!(ax, uni.lat[ord], uni.cover[ord]; label = "uni", color = :orange)
lines!(ax, sort(latT), tccE[sortperm(latT)]; label = "ERA5", color = :black)
axislegend(ax; position = :lb)

ax = Axis(fig[4, 3]; title = "Area-mean RH profile",
          xlabel = "RH (%)", ylabel = "pressure (hPa)", yreversed = true)
lines!(ax, rotP, pm; label = "rot", color = :firebrick)
lines!(ax, uniP, pm; label = "uni", color = :orange)
lines!(ax, eraP, pm; label = "ERA5", color = :black)
axislegend(ax; position = :lt)

Label(fig[0, :], "Coupled RCE humidity & cloud vs ERA5 (zonal/annual mean)";
      fontsize = 20, font = :bold)

save(joinpath(OUTDIR, outname), fig)
println("\nwrote $(joinpath(OUTDIR, outname))")
