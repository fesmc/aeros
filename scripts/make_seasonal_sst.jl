# Build a zonally-symmetric seasonal-SST boundary condition for the aquaplanet
# from the ERA5 1991-2020 monthly climatology, for the prescribed-SST seasonal
# validation (option A: the slab is unstable on the SW-faithful core, so we drive
# the atmosphere with a realistic seasonally-varying surface instead).
#
# ERA5 skin temperature (skt) is ocean-masked with the land-sea mask (lsm) and
# zonal-averaged over ocean points to SST(lat, month), then interpolated to the
# model's Gaussian latitudes. Latitudes with no ocean (e.g. the Antarctic cap)
# fall back to the all-points zonal mean so the aquaplanet has a defined SST
# everywhere. The result is zonally symmetric -- the consistent choice for an
# aquaplanet with no land.
#
#   julia scripts/make_seasonal_sst.jl [model_ref.nc] [out.nc]
#
# model_ref.nc supplies the target model latitudes (any rce_long dump has `lat`).
# Writes input/sst_seasonal_t21.nc: sst(lat, month) [K] on the model lat grid.

using NCDatasets
using Statistics
using Printf

const ERA5_SL = "/Users/alrobi001/data/era5/monthly-single-levels"
slc(v) = joinpath(ERA5_SL, "era5_monthly-single-levels_$(v)_1991-2020_clim.nc")

ref_nc = length(ARGS) >= 1 ? ARGS[1] : "output/seasonal_smoke/rh_swf.nc"
out_nc = length(ARGS) >= 2 ? ARGS[2] : "input/sst_seasonal_t21.nc"

# --- model target latitudes ---------------------------------------------
latM = NCDataset(ref_nc) do ds; Array(ds["lat"][:]); end
nlat = length(latM)

# --- ERA5 skt + lsm, ocean-masked zonal mean -> SST(latE, 12) -----------
dskt = NCDataset(slc("skt")); dlsm = NCDataset(slc("lsm"))
skt = Array(dskt["skt"][:, :, :])          # (lon, lat, month)
lsm = Array(dlsm["lsm"][:, :, :])          # (lon, lat, month)
latE = Array(dskt["latitude"][:])
close(dskt); close(dlsm)

nlonE, nlatE, nmon = size(skt)
@assert nmon == 12
ocean = mean(lsm; dims = 3)[:, :, 1] .< 0.5   # (lon, lat): ocean where lsm < 0.5

sstE = Matrix{Float64}(undef, nlatE, nmon)
for m in 1:nmon, j in 1:nlatE
    oc = @view ocean[:, j]
    col = @view skt[:, j, m]
    sstE[j, m] = any(oc) ? mean(col[oc]) : mean(col)   # ocean mean, else all-points
end

# --- interpolate ERA5 lat -> model lat (per month) ----------------------
function interp1(x, y, xq)
    p = sortperm(x); xs = x[p]; ys = y[p]
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    k = searchsortedlast(xs, xq)
    t = (xq - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] * (1 - t) + ys[k+1] * t
end

sstM = Matrix{Float64}(undef, nlat, nmon)
for m in 1:nmon, j in 1:nlat
    sstM[j, m] = interp1(latE, sstE[:, m], latM[j])
end

# --- write --------------------------------------------------------------
isfile(out_nc) && rm(out_nc)
NCDataset(out_nc, "c") do ds
    defDim(ds, "lat", nlat); defDim(ds, "month", nmon)
    v = defVar(ds, "lat", Float64, ("lat",)); v[:] = latM
    v.attrib["units"] = "degrees_north"
    mo = defVar(ds, "month", Float64, ("month",)); mo[:] = collect(1.0:nmon)
    s = defVar(ds, "sst", Float64, ("lat", "month"))
    s[:, :] = sstM
    s.attrib["units"] = "K"
    s.attrib["long_name"] = "zonally-symmetric ocean-masked ERA5 skt climatology"
end

@printf("wrote %s : sst(%d lat, %d month)\n", out_nc, nlat, nmon)
# quick sanity: seasonal amplitude at a few latitudes
for L in (50.0, 0.0, -50.0)
    j = argmin(abs.(latM .- L))
    amp = maximum(sstM[j, :]) - minimum(sstM[j, :])
    @printf("  lat %+5.1f : SST %.1f-%.1f K  (seasonal amp %.1f K)\n",
            latM[j], minimum(sstM[j, :]), maximum(sstM[j, :]), amp)
end
