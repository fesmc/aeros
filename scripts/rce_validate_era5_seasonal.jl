# Seasonal validation of the prescribed-SST RCE against the ERA5 monthly
# climatology (M2 seasonal, option A).
#
# The slab aquaplanet is unstable on the SW-faithful core (docs/m2_handoff.md), so
# the seasonal cycle is driven by prescribed seasonal SST (ERA5 skt climatology,
# scripts/make_seasonal_sst.jl) plus seasonal insolation. This judges the
# *atmospheric* seasonal response -- the Jan/Jul zonal-mean T, u, RH cross-sections
# and the near-surface-air-T seasonal amplitude & phase by latitude -- against
# ERA5 1991-2020. (Eddies are off: seeded eddies + seasonal forcing blow up at the
# model top, so the jet is unbraked and biased strong -- a documented limitation.)
#
#   rce_long.x output/seasonal_ppSST/seasonal_ppSST.nml   (writes monthly.nc)
#   julia scripts/rce_validate_era5_seasonal.jl [monthly.nc]
#
# Writes docs/figures/rce_seasonal_xsec.png, rce_seasonal_amp.png and prints a
# jet-migration / amplitude summary.

using NCDatasets
using CairoMakie
using Statistics
using Printf

const ERA5_PL = "/Users/alrobi001/data/era5/monthly-pressure-levels"
const ERA5_SL = "/Users/alrobi001/data/era5/monthly-single-levels"
const OUTDIR  = "docs/figures"
plc(v) = joinpath(ERA5_PL, "era5_monthly-pressure-levels_$(v)_1991-2020_clim.nc")
slc(v) = joinpath(ERA5_SL, "era5_monthly-single-levels_$(v)_1991-2020_clim.nc")

model_nc = length(ARGS) >= 1 ? ARGS[1] : "output/seasonal_ppSST/monthly.nc"

# --- readers -------------------------------------------------------------
# NCDatasets reads dim1-first, so the model (dim1=lat,dim2=lev,dim3=month) comes
# in as (lat, lev, month) directly.
function read_model(path)
    ds = NCDataset(path)
    lat = Array(ds["lat"][:]); sig = Array(ds["lev"][:])
    m = (lat=lat, p=sig .* 1000.0,                      # sigma -> hPa (aquaplanet ps~1000)
         t=Array(ds["t"][:,:,:]), u=Array(ds["u"][:,:,:]), rh=Array(ds["rh"][:,:,:]),
         sst=Array(ds["sst"][:,:]))
    close(ds); return m
end

# ERA5 pressure-level var -> zonal mean (lat, plev, month)
function era5_pl(var)
    ds = NCDataset(plc(var))
    A = Array(ds[var][:,:,:,:])                          # (lon, lat, plev, month)
    lat = Array(ds["latitude"][:]); plev = Array(ds["pressure_level"][:])
    close(ds)
    return lat, plev, dropdims(mean(A; dims=1); dims=1)  # (lat, plev, month)
end

# ERA5 single-level var -> zonal mean (lat, month)
function era5_sl(var)
    ds = NCDataset(slc(var))
    A = Array(ds[var][:,:,:])                            # (lon, lat, month)
    lat = Array(ds["latitude"][:]); close(ds)
    return lat, dropdims(mean(A; dims=1); dims=1)        # (lat, month)
end

function interp1(x, y, xq)
    p = sortperm(x); xs = x[p]; ys = y[p]
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    k = searchsortedlast(xs, xq); t = (xq - xs[k]) / (xs[k+1] - xs[k])
    return ys[k]*(1-t) + ys[k+1]*t
end

# regrid ERA5 (lat,plev) at one month onto the model (lat,p) grid
function regrid(latE, plevE, E, latM, pM)
    tmp = [interp1(latE, E[:, jp], la) for la in latM, jp in eachindex(plevE)]
    return [interp1(plevE, tmp[jl, :], p) for jl in eachindex(latM), p in pM]
end

# --- load ----------------------------------------------------------------
m = read_model(model_nc)
ord = sortperm(m.lat); latM = m.lat[ord]; pM = m.p
latE, plevE, tE = era5_pl("t")
_,    _,     uE = era5_pl("u")
_,    _,     rE = era5_pl("r")
latEs, t2mE = era5_sl("t2m")                             # near-surface air T

JAN, JUL = 1, 7
mfield(A, mo) = A[ord, :, mo]                            # model (lat,lev) at month mo

# --- Figure 1: Jan/Jul cross-sections, model vs ERA5 --------------------
mkpath(OUTDIR)
fig1 = Figure(size = (1750, 1250), fontsize = 15)
function xsec!(fig, r, c, x, p, F, ttl; crange, cmap)
    ax = Axis(fig[r, c]; title=ttl, xlabel="latitude", ylabel="pressure (hPa)", yreversed=true)
    hm = heatmap!(ax, x, p, F; colormap=cmap, colorrange=crange)
    Colorbar(fig[r, c+1], hm); return ax
end
specs = [(:t,  tE, "T (K)",   (200,305), :thermal),
         (:u,  uE, "u (m/s)", (-45,45),  :balance),
         (:rh, rE, "RH (%)",  (0,100),   :viridis)]
for (ri, (sym, E, lbl, cr, cm)) in enumerate(specs)
    Mjan = mfield(getproperty(m, sym), JAN); Mjul = mfield(getproperty(m, sym), JUL)
    Ejan = regrid(latE, plevE, E[:,:,JAN], latM, pM)
    Ejul = regrid(latE, plevE, E[:,:,JUL], latM, pM)
    xsec!(fig1, ri, 1, latM, pM, Mjan, "model $lbl — Jan"; crange=cr, cmap=cm)
    xsec!(fig1, ri, 3, latM, pM, Ejan, "ERA5 $lbl — Jan";  crange=cr, cmap=cm)
    xsec!(fig1, ri, 5, latM, pM, Mjul, "model $lbl — Jul"; crange=cr, cmap=cm)
    xsec!(fig1, ri, 7, latM, pM, Ejul, "ERA5 $lbl — Jul";  crange=cr, cmap=cm)
end
Label(fig1[0, :], "Seasonal zonal-mean cross-sections: model (prescribed seasonal SST) vs ERA5",
      fontsize=19, font=:bold)
save(joinpath(OUTDIR, "rce_seasonal_xsec.png"), fig1)

# --- Figure 2: near-surface-air-T seasonal amplitude & phase, jet migration ---
# model near-surface air T = lowest sigma level; ERA5 = t2m.
Tsfc_M = m.t[ord, end, :]                                # (lat, month)
ampM = vec(maximum(Tsfc_M; dims=2) .- minimum(Tsfc_M; dims=2))
phM  = [argmax(Tsfc_M[j, :]) for j in eachindex(latM)]
ampE = [maximum(t2mE[j,:]) - minimum(t2mE[j,:]) for j in eachindex(latEs)]
phE  = [argmax(t2mE[j,:]) for j in eachindex(latEs)]

# jet: peak |u| latitude by month (winter-hemisphere migration)
function jetlat(U, lat, mo)   # U model (lat,lev,month)
    col = dropdims(maximum(U[:, :, mo]; dims=2); dims=2)  # peak over levels, per lat
    j = argmax(col); return lat[j], col[j]
end
uM = m.u[ord, :, :]

fig2 = Figure(size = (1500, 520), fontsize = 15)
ax1 = Axis(fig2[1,1]; title="Near-surface-air-T seasonal amplitude (max−min)",
           xlabel="latitude", ylabel="amplitude (K)")
lines!(ax1, latM, ampM, label="model (lowest level)", linewidth=3)
lines!(ax1, latEs, ampE, label="ERA5 t2m", linewidth=3, linestyle=:dash)
axislegend(ax1; position=:ct)
ax2 = Axis(fig2[1,2]; title="Month of warmest near-surface air",
           xlabel="latitude", ylabel="month of max", yticks=1:2:12)
scatter!(ax2, latM, phM, label="model", markersize=9)
scatter!(ax2, latEs, phE, label="ERA5 t2m", markersize=6)
axislegend(ax2; position=:ct)
ax3 = Axis(fig2[1,3]; title="Jet latitude by month (winter-hemisphere migration)",
           xlabel="month", ylabel="peak-u latitude")
jlat = [jetlat(uM, latM, mo)[1] for mo in 1:12]
lines!(ax3, 1:12, jlat, linewidth=3)
scatter!(ax3, 1:12, jlat, markersize=9)
save(joinpath(OUTDIR, "rce_seasonal_amp.png"), fig2)

# --- text summary --------------------------------------------------------
@printf("\nNear-surface-air-T seasonal amplitude (K):\n  lat   model   ERA5\n")
for L in (60.0, 30.0, 0.0, -30.0, -60.0)
    jm = argmin(abs.(latM .- L)); je = argmin(abs.(latEs .- L))
    @printf("  %+4.0f   %5.1f   %5.1f\n", L, ampM[jm], ampE[je])
end
@printf("\nJet peak-u latitude: Jan % .0f°, Jul % .0f°  (winter-hemisphere)\n",
        jetlat(uM, latM, JAN)[1], jetlat(uM, latM, JUL)[1])
println("\nwrote $(OUTDIR)/rce_seasonal_xsec.png and rce_seasonal_amp.png")
