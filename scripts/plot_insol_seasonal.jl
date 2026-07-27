# Seasonal-cycle capability check: the model's TOA insolation over a year, from
# aeros_insolation (fesmc/insol, Laskar 2004). Reads output/insol_seasonal.nc
# (written by probe_insol.x) and plots the latitude x day-of-year march plus the
# annual-mean profile.
#
#   ./libaeros/bin/probe_insol.x 42 0
#   julia scripts/plot_insol_seasonal.jl
#
# Writes docs/figures/insol_seasonal.png.

using NCDatasets, CairoMakie, Statistics, Printf

ds = NCDataset("output/insol_seasonal.nc")
lat = ds["lat"][:]; day = ds["day"][:]
sw  = Array(ds["sw_toa"][:, :]); swann = ds["sw_ann"][:]
close(ds)
sw = size(sw) == (length(lat), length(day)) ? sw : permutedims(sw)   # (lat, day)

ord = sortperm(lat); lats = lat[ord]

fig = Figure(size = (1150, 460), fontsize = 15)

ax = Axis(fig[1, 1]; title = "TOA insolation (W m⁻²) — seasonal march",
          xlabel = "day of year", ylabel = "latitude")
hm = heatmap!(ax, day, lats, permutedims(sw[ord, :]); colormap = :solar,
              colorrange = (0, maximum(sw)))
Colorbar(fig[1, 2], hm)
# solstices / equinoxes
vlines!(ax, [80.0, 172.0, 266.0, 355.0]; color = (:white, 0.5), linestyle = :dash)

ax2 = Axis(fig[1, 3]; title = "Annual-mean insolation",
           xlabel = "W m⁻²", ylabel = "latitude")
lines!(ax2, swann[ord], lats; color = :firebrick, linewidth = 2)
vlines!(ax2, [Float64(1361) / 4]; color = :black, linestyle = :dash)  # S0/4

Label(fig[0, :], "aeros seasonal insolation via insol (Laskar 2004), present-day";
      fontsize = 18, font = :bold)
colsize!(fig.layout, 1, Relative(0.6))

mkpath("docs/figures")
save("docs/figures/insol_seasonal.png", fig)
println("wrote docs/figures/insol_seasonal.png")
