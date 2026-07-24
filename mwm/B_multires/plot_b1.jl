# plot_b1.jl  --  maps for the b1 multi-resolution error diagnosis.
#
# Decoupled from the heavy integration: reads the NetCDF written by b1_probe.jl
# and produces PNG maps with CairoMakie (headless). Can be run standalone after
# the compute job:
#
#     julia --project=. plot_b1.jl output/b1_fields.nc output
#
# or imported and called as plot_b1(ncfile, outdir).

using CairoMakie
using NCDatasets
using Printf
CairoMakie.activate!(type = "png")

# symmetric diverging map of one [lon,lat] variable
function _map!(fig, pos, lon, lat, M, title, unit)
    amax = maximum(abs, M)
    amax = amax == 0 ? 1.0 : amax
    ax = Axis(fig[pos...]; title = title, xlabel = "lon", ylabel = "lat",
              xticks = 0:60:360, yticks = -90:30:90)
    hm = heatmap!(ax, lon, lat, M; colormap = :RdBu, colorrange = (-amax, amax))
    Colorbar(fig[pos[1], pos[2]+1], hm; label = unit)
    xlims!(ax, extrema(lon)); ylims!(ax, extrema(lat))
    return ax
end

function plot_b1(ncfile::AbstractString, outdir::AbstractString)
    mkpath(outdir)
    ds = NCDataset(ncfile, "r")
    lon = ds["lon"][:]; lat = ds["lat"][:]
    thi = get(ds.attrib, "trunc_hi", "hi"); tlo = get(ds.attrib, "trunc_lo", "lo")
    php = get(ds.attrib, "approx_pressure_hPa", 500)
    read2(name) = Array(ds[name][:, :])

    # --- Z* triptych: hi->lo, lo, difference ---------------------------------
    Zhi = read2("Zstar_hi"); Zlo = read2("Zstar_lo"); dZ = read2("dZstar")
    fig = Figure(size = (1500, 420))
    _map!(fig, (1, 1), lon, lat, Zhi, "Z* T$thi -> T$tlo (~$(round(php))hPa)", "m")
    _map!(fig, (1, 3), lon, lat, Zlo, "Z* T$tlo", "m")
    _map!(fig, (1, 5), lon, lat, dZ,  "ΔZ* = T$thi(->T$tlo) - T$tlo", "m")
    p1 = joinpath(outdir, "b1_Zstar.png"); save(p1, fig)

    # --- ΔZ* alone (larger) --------------------------------------------------
    fig2 = Figure(size = (720, 460))
    _map!(fig2, (1, 1), lon, lat, dZ, "Resolution error ΔZ* (~$(round(php))hPa)", "m")
    p2 = joinpath(outdir, "b1_dZstar.png"); save(p2, fig2)

    # --- ΔT* and Δu* ---------------------------------------------------------
    dT = read2("dTstar"); dU = read2("dUstar")
    fig3 = Figure(size = (1050, 440))
    _map!(fig3, (1, 1), lon, lat, dT, "Resolution error ΔT*", "K")
    _map!(fig3, (1, 3), lon, lat, dU, "Resolution error Δu*", "m/s")
    p3 = joinpath(outdir, "b1_dTstar_dUstar.png"); save(p3, fig3)

    close(ds)
    return [p1, p2, p3]
end

# standalone entry point
if abspath(PROGRAM_FILE) == @__FILE__
    ncfile = length(ARGS) >= 1 ? ARGS[1] : error("usage: julia plot_b1.jl <b1_fields.nc> [outdir]")
    outdir = length(ARGS) >= 2 ? ARGS[2] : dirname(ncfile)
    paths = plot_b1(ncfile, outdir)
    println("wrote:"); foreach(p -> println("  ", p), paths)
end
