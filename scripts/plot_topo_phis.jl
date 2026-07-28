# Map of the model-grid surface geopotential (T21 orography) written by
# drivers/rce_long.f90 when l_topography=.true. A visual sanity check that the
# ERA5 orography regridded onto the model Gaussian grid has the right continents
# and mountains and a sane magnitude.
#
#   julia scripts/plot_topo_phis.jl [input.nc] [outdir]

using NCDatasets
using CairoMakie
using Printf

const GRAV = 9.80665

infile = length(ARGS) >= 1 ? ARGS[1] : "output/rce_phis.nc"
outdir = length(ARGS) >= 2 ? ARGS[2] : "docs/figures"

ds   = NCDataset(infile)
lon  = ds["lon"][:]
lat  = ds["lat"][:]
phis = ds["phis"][:, :]          # (lon, lat)
close(ds)

# Latitude runs north -> south on the model grid; flip to ascending for a sane
# axis (and reorder the field to match).
ord   = sortperm(lat)
lats  = lat[ord]
elev  = (phis ./ GRAV)[:, ord]   # surface elevation [m] = phis/g

vmax = maximum(abs, elev)
fig  = Figure(size = (900, 460))
ax   = Axis(fig[1, 1];
            title  = @sprintf("aeros T21 orography (phis/g)  [max %.0f m, min %.0f m]",
                              maximum(elev), minimum(elev)),
            xlabel = "longitude [°E]", ylabel = "latitude [°N]")
hm = heatmap!(ax, lon, lats, elev; colormap = :terrain, colorrange = (-vmax, vmax))
Colorbar(fig[1, 2], hm; label = "surface elevation [m]")

mkpath(outdir)
outfile = joinpath(outdir, "topo_phis.png")
save(outfile, fig)
println("wrote ", outfile)
