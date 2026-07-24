# plot_b2b.jl  --  maps for the b2 approach-B (nudging-derived <N>) correction test.
#
# Decoupled from the heavy integration: reads the NetCDF written by b2b_probe.jl and
# produces PNG maps with CairoMakie (headless). Standalone:
#
#     julia --project=. plot_b2b.jl output_b2b/b2b_fields.nc output_b2b
#
# or imported and called as plot_b2b(ncfile, outdir).

using CairoMakie
using NCDatasets
using Printf
CairoMakie.activate!(type = "png")

function _map!(fig, pos, lon, lat, M, title, unit; amax = nothing)
    a = amax === nothing ? maximum(abs, M) : amax
    a = a == 0 ? 1.0 : a
    ax = Axis(fig[pos...]; title = title, xlabel = "lon", ylabel = "lat",
              xticks = 0:60:360, yticks = -90:30:90)
    hm = heatmap!(ax, lon, lat, M; colormap = :RdBu, colorrange = (-a, a))
    Colorbar(fig[pos[1], pos[2] + 1], hm; label = unit)
    xlims!(ax, extrema(lon)); ylims!(ax, extrema(lat))
    return ax
end

function plot_b2b(ncfile::AbstractString, outdir::AbstractString)
    mkpath(outdir)
    ds = NCDataset(ncfile, "r")
    lon = ds["lon"][:]; lat = ds["lat"][:]
    thi = get(ds.attrib, "trunc_hi", "hi"); tlo = get(ds.attrib, "trunc_lo", "lo")
    php = get(ds.attrib, "approx_pressure_hPa", 500)
    tau = get(ds.attrib, "nudge_tau_hours", "?")
    cf  = get(ds.attrib, "correct_fields", "")
    read2(name) = Array(ds[name][:, :])

    paths = String[]
    tag = "<N> tau=$(tau)h [$cf]"

    # --- Z* triptych: target / bare / corrected (shared colour scale) ----------
    Ztg = read2("Zstar_target"); Zba = read2("Zstar_bare"); Zco = read2("Zstar_corr")
    amaxZ = maximum(abs, vcat(vec(Ztg), vec(Zba), vec(Zco)))
    fig = Figure(size = (1500, 420))
    _map!(fig, (1, 1), lon, lat, Ztg, "Z* target (T$thi->T$tlo, ~$(round(php))hPa)", "m"; amax = amaxZ)
    _map!(fig, (1, 3), lon, lat, Zba, "Z* bare T$tlo", "m"; amax = amaxZ)
    _map!(fig, (1, 5), lon, lat, Zco, "Z* T$tlo + $tag", "m"; amax = amaxZ)
    p = joinpath(outdir, "b2b_Zstar_three.png"); save(p, fig); push!(paths, p)

    # --- Z* residual maps: bare-target vs corrected-target (shared scale) -------
    rZb = read2("resZ_bare"); rZc = read2("resZ_corr")
    amaxR = maximum(abs, vcat(vec(rZb), vec(rZc)))
    fig2 = Figure(size = (1050, 440))
    _map!(fig2, (1, 1), lon, lat, rZb, "ΔZ* residual: bare T$tlo - target", "m"; amax = amaxR)
    _map!(fig2, (1, 3), lon, lat, rZc, "ΔZ* residual: T$tlo+$tag - target", "m"; amax = amaxR)
    p = joinpath(outdir, "b2b_Zstar_residual.png"); save(p, fig2); push!(paths, p)

    # --- T* and U* residuals: bare vs corrected --------------------------------
    for (fld, unit) in (("T", "K"), ("U", "m/s"))
        rb = read2("res$(fld)_bare"); rc = read2("res$(fld)_corr")
        amax = maximum(abs, vcat(vec(rb), vec(rc)))
        fig3 = Figure(size = (1050, 440))
        _map!(fig3, (1, 1), lon, lat, rb, "Δ$(fld)* residual: bare - target", unit; amax = amax)
        _map!(fig3, (1, 3), lon, lat, rc, "Δ$(fld)* residual: +$tag - target", unit; amax = amax)
        p = joinpath(outdir, "b2b_$(fld)star_residual.png"); save(p, fig3); push!(paths, p)
    end

    close(ds)
    return paths
end

# standalone entry point
if abspath(PROGRAM_FILE) == @__FILE__
    ncfile = length(ARGS) >= 1 ? ARGS[1] : error("usage: julia plot_b2b.jl <b2b_fields.nc> [outdir]")
    outdir = length(ARGS) >= 2 ? ARGS[2] : dirname(ncfile)
    paths = plot_b2b(ncfile, outdir)
    println("wrote:"); foreach(p -> println("  ", p), paths)
end
