# ERA5 moist-line validation of the coupled RCE (M2 Tier 2).
#
# The M2 validation target (docs/design.md): the coupled radiative-convective
# equilibrium against the ERA5 1991-2020 climatology in the zonal mean —
# temperature structure, the zonal-wind jet, and relative humidity. This is the
# "moist-line": the moist physics stack (radiation + convection + condensation +
# surface + vdiff) run to equilibrium and judged against climate, not against
# conservation. The clear-sky line (§14) and the cloud radiative effect (§17-18,
# §21-23) are done; this closes the loop on the equilibrated dynamical/thermal
# state.
#
# Driven on the rotating vehicle at the tuned cond_rh_crit (§23), whose dump now
# carries zonal-mean u and v as well as T, q, RH:
#   rce_long.x logs/rce_rot_cr_0.93.nml output/rce_rh_rot_cr_0.93.nc
#
#   julia scripts/rce_validate_era5.jl [model.nc]
#
# Writes docs/figures/rce_validate_era5.png and prints a jet / T-bias summary.

using NCDatasets
using CairoMakie
using Statistics
using Printf

const ERA5_PL = "/Users/alrobi001/data/era5/monthly-pressure-levels"
const OUTDIR  = "docs/figures"

model_nc = length(ARGS) >= 1 ? ARGS[1] : "output/rce_rh_rot_cr_0.93.nc"
outname  = length(ARGS) >= 2 ? ARGS[2] : "rce_validate_era5.png"
plc(v) = joinpath(ERA5_PL, "era5_monthly-pressure-levels_$(v)_1991-2020_clim.nc")

# --- readers (shared shape with rce_humidity_vs_era5.jl) -----------------

function read_model(path)
    ds = NCDataset(path)
    lat = ds["lat"][:]; lev = ds["lev"][:]
    nlat = length(lat); nlev = length(lev)
    grab(v) = (A = Array(ds[v][:, :]); size(A) == (nlat, nlev) ? A : permutedims(A))
    m = (lat=lat, lev=lev, t=grab("t"), rh=grab("rh"),
         u=grab("u"), v=grab("v"), pfull=grab("pfull"))
    close(ds); return m
end

function era5_pl_zonal(var)
    ds = NCDataset(plc(var))
    A = Array(ds[var][:, :, :, :])                    # (lon, lat, plev, month)
    lat = ds["latitude"][:]; plev = ds["pressure_level"][:]
    close(ds)
    zm = dropdims(mean(A; dims=(1, 4)); dims=(1, 4))   # (lat, plev)
    return lat, plev, zm
end

# --- regrid helpers ------------------------------------------------------

function interp1(x, y, xq)
    p = sortperm(x); xs = x[p]; ys = y[p]
    xq <= xs[1]   && return ys[1]
    xq >= xs[end] && return ys[end]
    k = searchsortedlast(xs, xq)
    t = (xq - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] * (1 - t) + ys[k+1] * t
end

function regrid_to_model(latE, plevE, E, latM, pm)
    tmp = Matrix{Float64}(undef, length(latM), length(plevE))
    for (jl, la) in enumerate(latM), jp in eachindex(plevE)
        tmp[jl, jp] = interp1(latE, E[:, jp], la)
    end
    out = Matrix{Float64}(undef, length(latM), length(pm))
    for jl in axes(out, 1), (jp, p) in enumerate(pm)
        out[jl, jp] = interp1(plevE, tmp[jl, :], p)
    end
    return out
end

latmean(A, lat) = (w = cos.(deg2rad.(lat)); w ./= sum(w); vec(sum(A .* w; dims=1)))

# peak zonal wind in a hemisphere: (lat, p, value) at the max of u
function jet(u, lat, p; north=true)
    best = (-Inf, 0.0, 0.0)
    for jl in eachindex(lat), jp in eachindex(p)
        (north ? lat[jl] > 0 : lat[jl] < 0) || continue
        u[jl, jp] > best[1] && (best = (u[jl, jp], lat[jl], p[jp]))
    end
    return best  # (umax, lat, p)
end

# --- load ----------------------------------------------------------------

m = read_model(model_nc)
latE, plevE, tE = era5_pl_zonal("t")
_,    _,     uE = era5_pl_zonal("u")
_,    _,     rE = era5_pl_zonal("r")

ord  = sortperm(m.lat)
latM = m.lat[ord]
pm   = latmean(m.pfull, m.lat)          # model lat-mean pressure axis (hPa)

tM, uM, rhM = m.t[ord, :], m.u[ord, :], m.rh[ord, :]
tE_m  = regrid_to_model(latE, plevE, tE, latM, pm)
uE_m  = regrid_to_model(latE, plevE, uE, latM, pm)
rE_m  = regrid_to_model(latE, plevE, rE, latM, pm)

# --- text summary --------------------------------------------------------

println("\nJet (peak zonal-mean u):")
for (nm, u, la, p) in (("model", uM, latM, pm), ("ERA5 ", uE_m, latM, pm))
    jn = jet(u, la, p; north=true); js = jet(u, la, p; north=false)
    @printf("  %s   NH % .1f m/s @ %+.0f°, %.0f hPa   SH % .1f m/s @ %+.0f°, %.0f hPa\n",
            nm, jn[1], jn[2], jn[3], js[1], js[2], js[3])
end

tbM = latmean(tM, latM); tbE = latmean(tE_m, latM)
println("\nArea-mean T (K) by level (model top -> surface):")
println("   p[hPa]   model   ERA5   model-ERA5")
for k in eachindex(pm)
    @printf("  %7.1f  %6.1f  %6.1f  %+6.1f\n", pm[k], tbM[k], tbE[k], tbM[k]-tbE[k])
end

# --- figure --------------------------------------------------------------

mkpath(OUTDIR)
fig = Figure(size = (1500, 1150), fontsize = 15)

function xsec!(row, col, x, p, F, ttl; crange, cmap)
    ax = Axis(fig[row, col]; title = ttl, xlabel = "latitude",
              ylabel = "pressure (hPa)", yreversed = true)
    hm = heatmap!(ax, x, p, F; colormap = cmap, colorrange = crange)
    Colorbar(fig[row, col+1], hm)
end

# Row 1: temperature
xsec!(1, 1, latM, pm, tM,   "Model T (K) — rotating"; crange=(200, 300), cmap=:thermal)
xsec!(1, 3, latM, pm, tE_m, "ERA5 T (K)";              crange=(200, 300), cmap=:thermal)
xsec!(2, 1, latM, pm, tM .- tE_m, "T bias  model − ERA5 (K)"; crange=(-20, 20), cmap=:balance)

# Row 2 (right of T-bias): zonal wind
xsec!(2, 3, latM, pm, uM .- uE_m, "u bias  model − ERA5 (m/s)"; crange=(-30, 30), cmap=:balance)
xsec!(3, 1, latM, pm, uM,   "Model u (m/s) — rotating"; crange=(-40, 40), cmap=:balance)
xsec!(3, 3, latM, pm, uE_m, "ERA5 u (m/s)";             crange=(-40, 40), cmap=:balance)

# Row 4: RH cross-sections + jet-profile line
xsec!(4, 1, latM, pm, rhM,  "Model RH (%) — rotating"; crange=(0, 100), cmap=:viridis)
xsec!(4, 3, latM, pm, rE_m, "ERA5 RH (%)";             crange=(0, 100), cmap=:viridis)

Label(fig[0, :], "Coupled RCE vs ERA5 — moist-line validation (zonal/annual mean)";
      fontsize = 20, font = :bold)

save(joinpath(OUTDIR, outname), fig)
println("\nwrote $(joinpath(OUTDIR, outname))")
