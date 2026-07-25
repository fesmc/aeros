# Lat x lon maps of the clear-sky radiation validation against ERA5.
#
# For each of six clear-sky fluxes, a three-panel row: the model field driven
# on ERA5 columns, the ERA5 clear-sky diagnostic, and the bias (model - ERA5).
# Reads the NetCDF written by drivers/validate_era5.f90.
#
#   julia scripts/plot_era5_validation.jl [input.nc] [output.png]

using NCDatasets
using CairoMakie
using Statistics
using Printf

infile  = length(ARGS) >= 1 ? ARGS[1] : "output/era5_rad_validation.nc"
outfile = length(ARGS) >= 2 ? ARGS[2] : "docs/figures/era5_rad_validation.png"

const MV = -9990.0
mask(a) = replace(x -> (x <= MV || isnan(x)) ? missing : x, a)

# (title, model var, era var, bias var)
fluxes = [
    ("OLR (TOA up)",       "olr_mod",       "olr_era",       "olr_bias"),
    ("Surface down LW",    "lwdn_sfc_mod",  "lwdn_sfc_era",  "lwdn_sfc_bias"),
    ("Surface net LW",     "lwnet_sfc_mod", "lwnet_sfc_era", "lwnet_sfc_bias"),
    ("TOA net SW",         "swnet_toa_mod", "swnet_toa_era", "swnet_toa_bias"),
    ("Surface down SW",    "swdn_sfc_mod",  "swdn_sfc_era",  "swdn_sfc_bias"),
    ("Surface net SW",     "swnet_sfc_mod", "swnet_sfc_era", "swnet_sfc_bias"),
]

ds  = NCDataset(infile)
lon = ds["lon"][:]
lat = ds["lat"][:]

# ERA5 latitude runs 90 -> -90; flip to ascending for a sane axis.
latflip = issorted(lat) ? false : true
ord = latflip ? sortperm(lat) : collect(eachindex(lat))
lats = lat[ord]

fig = Figure(size = (1500, 340 * length(fluxes)), fontsize = 15)

# area weights for the global-mean annotation
w = cos.(deg2rad.(lats))
function gmean(A)
    num = 0.0; den = 0.0
    for j in eachindex(lats), i in eachindex(lon)
        v = A[i, j]
        ismissing(v) && continue
        num += v * w[j]; den += w[j]
    end
    den > 0 ? num / den : NaN
end

for (r, (title, vm, ve, vb)) in enumerate(fluxes)
    M = mask(ds[vm][:, ord]); E = mask(ds[ve][:, ord]); B = mask(ds[vb][:, ord])

    # shared color range for model & ERA5
    both = skipmissing(vcat(vec(M), vec(E)))
    lo, hi = quantile(collect(both), (0.02, 0.98))
    # symmetric range for the bias
    babs = quantile(collect(skipmissing(vec(B))), 0.98)
    babs = max(babs, abs(quantile(collect(skipmissing(vec(B))), 0.02)))

    for (c, (A, sub, cmap, rng)) in enumerate([
            (M, "model  ⟨$(@sprintf("%.1f", gmean(M)))⟩", :viridis, (lo, hi)),
            (E, "ERA5   ⟨$(@sprintf("%.1f", gmean(E)))⟩", :viridis, (lo, hi)),
            (B, "bias   ⟨$(@sprintf("%+.1f", gmean(B)))⟩", :balance, (-babs, babs)),
        ])
        gl = GridLayout(fig[r, c])
        ax = Axis(gl[1, 1];
                  title = c == 1 ? "$title  —  $sub" : sub,
                  xticksvisible = false, yticksvisible = false,
                  xticklabelsvisible = false, yticklabelsvisible = false,
                  aspect = DataAspect())
        hm = heatmap!(ax, lon, lats, A; colormap = cmap, colorrange = rng)
        Colorbar(gl[2, 1], hm; vertical = false, flipaxis = false, height = 10)
        rowgap!(gl, 4)
    end
end

Label(fig[0, :], "Clear-sky radiation: aeros operators on ERA5 columns vs ERA5 (annual mean, W m⁻²)",
      fontsize = 19, font = :bold)

mkpath(dirname(outfile))
save(outfile, fig)
close(ds)
println("wrote $outfile")
