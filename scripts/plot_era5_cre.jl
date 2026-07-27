# ERA5 cloud radiative effect (CRE) and moisture/cloud climatology — the
# observational target for the aeros cloudy-sky radiation branch (M2, Tier 2).
#
# CRE = all-sky − clear-sky net downward flux (positive = clouds warm), derived
# purely from the ERA5 1991–2020 monthly climatology; no model is involved. LW,
# SW and net are computed at both the TOA and the surface from the paired all-
# sky / clear-sky ERA5 flux diagnostics:
#
#   TOA:      ttr/ttrc (net thermal),  tsr/tsrc (net solar)
#   surface:  str/strc (net thermal),  ssr/ssrc (net solar)
#
# ERA5 fluxes are daily accumulations [J m⁻²], /86400 → W m⁻², downward-positive
# (same convention as drivers/validate_era5.f90, §14). A companion figure shows
# the ERA5 moisture/cloud climatology (column water vapour, total cloud cover,
# and zonal-mean RH and cloud-fraction cross-sections).
#
#   julia scripts/plot_era5_cre.jl
#
# Writes docs/figures/era5_cre.png and docs/figures/era5_moisture_clim.png.

using NCDatasets
using CairoMakie
using Statistics
using Printf

const ERA5_SL = "/Users/alrobi001/data/era5/monthly-single-levels"
const ERA5_PL = "/Users/alrobi001/data/era5/monthly-pressure-levels"
const OUTDIR  = "docs/figures"

sl(v) = joinpath(ERA5_SL, "era5_monthly-single-levels_$(v)_1991-2020_clim.nc")
pl(v) = joinpath(ERA5_PL, "era5_monthly-pressure-levels_$(v)_1991-2020_clim.nc")

# --- readers -------------------------------------------------------------

# annual-mean single-level field → (lon, lat), optionally scaled (flux: 1/86400)
function read_sl(var; scale = 1.0)
    ds = NCDataset(sl(var))
    A = Array(ds[var][:, :, :])                     # (lon, lat, month)
    lon = ds["longitude"][:]; lat = ds["latitude"][:]
    close(ds)
    am = dropdims(mean(A; dims = 3); dims = 3) .* scale
    return lon, lat, am
end

# annual + zonal mean pressure-level field → (lat, plev)
function read_pl_zonal(var)
    ds = NCDataset(pl(var))
    A = Array(ds[var][:, :, :, :])                  # (lon, lat, plev, month)
    lat = ds["latitude"][:]; plev = ds["pressure_level"][:]
    close(ds)
    zm = dropdims(mean(A; dims = (1, 4)); dims = (1, 4))  # (lat, plev)
    return lat, plev, zm
end

# area-weighted global mean over a (lon, lat) field, lat in degrees
function gmean(A, lat)
    w = cos.(deg2rad.(lat))
    num = 0.0; den = 0.0
    for j in axes(A, 2), i in axes(A, 1)
        v = A[i, j]
        (ismissing(v) || isnan(v)) && continue
        num += v * w[j]; den += w[j]
    end
    return den > 0 ? num / den : NaN
end

# --- CRE fields ----------------------------------------------------------

lon, lat, ttr  = read_sl("ttr";  scale = 1 / 86400)
_,   _,   ttrc = read_sl("ttrc"; scale = 1 / 86400)
_,   _,   tsr  = read_sl("tsr";  scale = 1 / 86400)
_,   _,   tsrc = read_sl("tsrc"; scale = 1 / 86400)
_,   _,   str  = read_sl("str";  scale = 1 / 86400)
_,   _,   strc = read_sl("strc"; scale = 1 / 86400)
_,   _,   ssr  = read_sl("ssr";  scale = 1 / 86400)
_,   _,   ssrc = read_sl("ssrc"; scale = 1 / 86400)

# CRE = all-sky − clear-sky (downward-positive net flux; + = clouds warm)
cre = [
    ("TOA LW CRE",      ttr  .- ttrc),
    ("TOA SW CRE",      tsr  .- tsrc),
    ("TOA net CRE",     (ttr .+ tsr) .- (ttrc .+ tsrc)),
    ("Surface LW CRE",  str  .- strc),
    ("Surface SW CRE",  ssr  .- ssrc),
    ("Surface net CRE", (str .+ ssr) .- (strc .+ ssrc)),
]

println("Area-weighted global-mean CRE (W m⁻²):")
for (name, A) in cre
    @printf("  %-16s % .2f\n", name, gmean(A, lat))
end

# --- Figure 1: CRE maps --------------------------------------------------

ord  = sortperm(lat)          # ERA5 lat is 90→−90; ascending for a sane axis
lats = lat[ord]

fig = Figure(size = (1500, 720), fontsize = 15)
for (k, (name, A)) in enumerate(cre)
    r = (k - 1) ÷ 3 + 1
    c = (k - 1) % 3 + 1
    F = A[:, ord]
    babs = quantile(collect(skipmissing(vec(F))), 0.98)
    babs = max(babs, abs(quantile(collect(skipmissing(vec(F))), 0.02)))
    gl = GridLayout(fig[r, c])
    ax = Axis(gl[1, 1];
              title = @sprintf("%s   ⟨%+.1f⟩", name, gmean(A, lat)),
              xticksvisible = false, yticksvisible = false,
              xticklabelsvisible = false, yticklabelsvisible = false,
              aspect = DataAspect())
    hm = heatmap!(ax, lon, lats, F; colormap = :balance, colorrange = (-babs, babs))
    Colorbar(gl[2, 1], hm; vertical = false, flipaxis = false, height = 10)
    rowgap!(gl, 4)
end
Label(fig[0, :], "ERA5 cloud radiative effect (all-sky − clear-sky, annual mean, W m⁻²)",
      fontsize = 19, font = :bold)
mkpath(OUTDIR)
save(joinpath(OUTDIR, "era5_cre.png"), fig)
println("wrote $(joinpath(OUTDIR, "era5_cre.png"))")

# --- Figure 2: moisture / cloud climatology ------------------------------

_, _, tcwv = read_sl("tcwv")               # kg m⁻² (column water vapour)
_, _, tcc  = read_sl("tcc")                # 0–1 (total cloud cover)
latp, plev, rh = read_pl_zonal("r")        # % relative humidity
_,    _,    cc = read_pl_zonal("cc")       # 0–1 cloud fraction

ordp  = sortperm(latp)
latps = latp[ordp]
# pressure axis: surface (large p) at the bottom
pord  = sortperm(plev; rev = true)
plevs = plev[pord]

fig2 = Figure(size = (1300, 780), fontsize = 15)

# (a) column water vapour
gl = GridLayout(fig2[1, 1])
ax = Axis(gl[1, 1]; title = @sprintf("Total column water vapour  ⟨%.1f kg m⁻²⟩", gmean(tcwv, lat)),
          xticksvisible = false, yticksvisible = false,
          xticklabelsvisible = false, yticklabelsvisible = false, aspect = DataAspect())
hm = heatmap!(ax, lon, lats, tcwv[:, ord]; colormap = :viridis)
Colorbar(gl[2, 1], hm; vertical = false, flipaxis = false, height = 10); rowgap!(gl, 4)

# (b) total cloud cover
gl = GridLayout(fig2[1, 2])
ax = Axis(gl[1, 1]; title = @sprintf("Total cloud cover  ⟨%.2f⟩", gmean(tcc, lat)),
          xticksvisible = false, yticksvisible = false,
          xticklabelsvisible = false, yticklabelsvisible = false, aspect = DataAspect())
hm = heatmap!(ax, lon, lats, tcc[:, ord]; colormap = :viridis, colorrange = (0, 1))
Colorbar(gl[2, 1], hm; vertical = false, flipaxis = false, height = 10); rowgap!(gl, 4)

# (c) zonal-mean relative humidity (lat × pressure)
gl = GridLayout(fig2[2, 1])
ax = Axis(gl[1, 1]; title = "Zonal-mean relative humidity (%)",
          xlabel = "latitude", ylabel = "pressure (hPa)", yreversed = true)
hm = heatmap!(ax, latps, plevs, rh[ordp, pord]; colormap = :viridis, colorrange = (0, 100))
Colorbar(gl[1, 2], hm)

# (d) zonal-mean cloud fraction (lat × pressure)
gl = GridLayout(fig2[2, 2])
ax = Axis(gl[1, 1]; title = "Zonal-mean cloud fraction",
          xlabel = "latitude", ylabel = "pressure (hPa)", yreversed = true)
hm = heatmap!(ax, latps, plevs, cc[ordp, pord]; colormap = :viridis, colorrange = (0, 0.5))
Colorbar(gl[1, 2], hm)

Label(fig2[0, :], "ERA5 moisture & cloud climatology (1991–2020 annual mean)",
      fontsize = 19, font = :bold)
save(joinpath(OUTDIR, "era5_moisture_clim.png"), fig2)
println("wrote $(joinpath(OUTDIR, "era5_moisture_clim.png"))")
