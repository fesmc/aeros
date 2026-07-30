# Land seasonal calibration: near-surface-air-T seasonal amplitude over LAND,
# model vs ERA5. The aquaplanet seasonal validation showed the NH-midlat
# amplitude far too weak because there was no land; this turns land on (coupled
# land run, dt=900) and judges the land seasonal amplitude -- the quantity c_soil
# (soil heat capacity) controls -- against ERA5 t2m over land.
#
# Both fields are put on the model grid; ERA5 lsm (regridded, threshold 0.5)
# defines land. For each latitude the land-only zonal-mean t2m gives a seasonal
# cycle; its max-min is the amplitude. A land-area-weighted global land mean is
# printed too.
#
#   julia scripts/land_validate_era5.jl [model_monthly.nc]
# Writes docs/figures/land_t2m_amp.png and a text summary.

using NCDatasets
using CairoMakie
using Statistics
using Printf

const ERA5_SL = "/Users/alrobi001/data/era5/monthly-single-levels"
const OUTDIR  = "docs/figures"
slc(v) = joinpath(ERA5_SL, "era5_monthly-single-levels_$(v)_1991-2020_clim.nc")
model_nc = length(ARGS) >= 1 ? ARGS[1] : "output/land_cal/monthly_csoil2e6.nc"

# --- model near-surface T (lon,lat,month) + grid --------------------------
dm = NCDataset(model_nc)
lonM = Array(dm["lon"][:]); latM = Array(dm["lat"][:])
t2mM = Array(dm["t2m"][:, :, :])              # (lon, lat, month)
close(dm)
nlonM, nlatM = length(lonM), length(latM)

# --- ERA5 t2m + lsm (lon,lat,month)/(lon,lat) ----------------------------
de = NCDataset(slc("t2m")); dl = NCDataset(slc("lsm"))
t2mE = Array(de["t2m"][:, :, :])              # (lon, lat, month)
lsmE = mean(Array(dl["lsm"][:, :, :]); dims = 3)[:, :, 1]
lonE = Array(de["longitude"][:]); latE = Array(de["latitude"][:])
close(de); close(dl)

# --- bilinear regrid ERA5 (lon,lat) -> model grid (lon periodic) ----------
function interp1(x, y, xq)
    p = sortperm(x); xs = x[p]; ys = y[p]
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    k = searchsortedlast(xs, xq); t = (xq - xs[k]) / (xs[k+1] - xs[k])
    return ys[k]*(1-t) + ys[k+1]*t
end
function regrid2d(lonS, latS, F, lonT, latT)   # F(lonS,latS) -> (lonT,latT)
    tmp = [interp1(latS, F[i, :], la) for i in eachindex(lonS), la in latT]  # (lonS, latT)
    return [interp1(lonS, tmp[:, jl], lo) for lo in lonT, jl in eachindex(latT)]
end

lsmM  = regrid2d(lonE, latE, lsmE, lonM, latM)          # ERA5 land frac on model grid
landM = lsmM .>= 0.5
t2mE_m = cat([regrid2d(lonE, latE, t2mE[:, :, m], lonM, latM) for m in 1:12]...; dims = 3)

# --- land-only zonal-mean seasonal amplitude by latitude ------------------
function land_amp(t2m, land)                    # t2m(lon,lat,12), land(lon,lat)
    nlat = size(t2m, 2)
    amp = fill(NaN, nlat); cyc = fill(NaN, nlat, 12)
    for j in 1:nlat
        msk = land[:, j]
        any(msk) || continue
        for m in 1:12
            cyc[j, m] = mean(t2m[msk, j, m])
        end
        amp[j] = maximum(cyc[j, :]) - minimum(cyc[j, :])
    end
    return amp, cyc
end
ampM, cycM = land_amp(t2mM,   landM)
ampE, cycE = land_amp(t2mE_m, landM)            # same land mask -> fair comparison

# --- figure ---------------------------------------------------------------
mkpath(OUTDIR)
fig = Figure(size = (1400, 560), fontsize = 15)
ax1 = Axis(fig[1,1]; title="Land near-surface-T seasonal amplitude (max−min)",
           xlabel="latitude", ylabel="amplitude (K)")
ok = .!isnan.(ampM)
lines!(ax1, latM[ok], ampM[ok], label="model (land)", linewidth=3)
lines!(ax1, latM[ok], ampE[ok], label="ERA5 t2m (land)", linewidth=3, linestyle=:dash)
axislegend(ax1; position=:ct)
# seasonal cycle at a NH-midlat land band
jj = findmin(abs.(latM .- 55))[2]
ax2 = Axis(fig[1,2]; title=@sprintf("Land seasonal cycle @ %.0f°N", latM[jj]),
           xlabel="month", ylabel="t2m (K)", xticks=1:2:12)
lines!(ax2, 1:12, cycM[jj,:], label="model", linewidth=3)
lines!(ax2, 1:12, cycE[jj,:], label="ERA5", linewidth=3, linestyle=:dash)
axislegend(ax2; position=:cb)
save(joinpath(OUTDIR, "land_t2m_amp.png"), fig)

# --- summary --------------------------------------------------------------
@printf("\nLand near-surface-T seasonal amplitude (K), model vs ERA5:\n  lat   model   ERA5\n")
for L in (65.0, 50.0, 35.0, 0.0, -35.0)
    j = findmin(abs.(latM .- L))[2]
    isnan(ampM[j]) && continue
    @printf("  %+4.0f   %5.1f   %5.1f\n", L, ampM[j], ampE[j])
end
w = cosd.(latM)
gm(a) = (m = .!isnan.(a); sum(a[m].*w[m])/sum(w[m]))
@printf("\nland-area-weighted mean amplitude: model %.1f K   ERA5 %.1f K\n", gm(ampM), gm(ampE))
println("wrote $(OUTDIR)/land_t2m_amp.png")
