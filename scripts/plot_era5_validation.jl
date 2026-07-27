# Lat x lon maps of the radiation validation against ERA5, from the NetCDF
# written by drivers/validate_era5.f90.
#
# Two figures:
#   era5_rad_validation.png  - six clear-sky fluxes, each a row of
#                              model / ERA5 / bias (model - ERA5).
#   era5_cre_validation.png  - the TOA cloud radiative effect (LW, SW, net),
#                              model vs ERA5 vs bias: the §17 target, now with
#                              the aeros cloudy operators reproducing it.
#
#   julia scripts/plot_era5_validation.jl [input.nc] [outdir]

using NCDatasets
using CairoMakie
using Statistics
using Printf

infile = length(ARGS) >= 1 ? ARGS[1] : "output/era5_rad_validation.nc"
outdir = length(ARGS) >= 2 ? ARGS[2] : "docs/figures"

const MV = -9990.0
mask(a) = replace(x -> (x <= MV || isnan(x)) ? missing : x, a)

ds  = NCDataset(infile)
lon = ds["lon"][:]
lat = ds["lat"][:]

# ERA5 latitude runs 90 -> -90; flip to ascending for a sane axis.
ord  = issorted(lat) ? collect(eachindex(lat)) : sortperm(lat)
lats = lat[ord]
w    = cos.(deg2rad.(lats))

function gmean(A)
    num = 0.0; den = 0.0
    for j in eachindex(lats), i in eachindex(lon)
        v = A[i, j]
        ismissing(v) && continue
        num += v * w[j]; den += w[j]
    end
    return den > 0 ? num / den : NaN
end

# Build a figure: one row per (title, model var, era var, bias var), each row
# model / ERA5 / bias. `cmap` is the model/ERA5 colormap; :viridis for one-sign
# fluxes, :balance for the (signed) cloud radiative effect.
function make_figure(fluxes, suptitle, outpng; cmap = :viridis)
    fig = Figure(size = (1500, 340 * length(fluxes)), fontsize = 15)
    for (r, (title, vm, ve, vb)) in enumerate(fluxes)
        M = mask(ds[vm][:, ord]); E = mask(ds[ve][:, ord]); B = mask(ds[vb][:, ord])

        both = collect(skipmissing(vcat(vec(M), vec(E))))
        lo, hi = quantile(both, (0.02, 0.98))
        if cmap === :balance                     # symmetric range about zero
            a = max(abs(lo), abs(hi)); lo, hi = -a, a
        end
        babs = quantile(collect(skipmissing(vec(B))), 0.98)
        babs = max(babs, abs(quantile(collect(skipmissing(vec(B))), 0.02)))

        for (c, (A, sub, cm, rng)) in enumerate([
                (M, "model  ⟨$(@sprintf("%+.1f", gmean(M)))⟩", cmap, (lo, hi)),
                (E, "ERA5   ⟨$(@sprintf("%+.1f", gmean(E)))⟩", cmap, (lo, hi)),
                (B, "bias   ⟨$(@sprintf("%+.1f", gmean(B)))⟩", :balance, (-babs, babs)),
            ])
            gl = GridLayout(fig[r, c])
            ax = Axis(gl[1, 1];
                      title = c == 1 ? "$title  —  $sub" : sub,
                      xticksvisible = false, yticksvisible = false,
                      xticklabelsvisible = false, yticklabelsvisible = false,
                      aspect = DataAspect())
            hm = heatmap!(ax, lon, lats, A; colormap = cm, colorrange = rng)
            Colorbar(gl[2, 1], hm; vertical = false, flipaxis = false, height = 10)
            rowgap!(gl, 4)
        end
    end
    Label(fig[0, :], suptitle, fontsize = 19, font = :bold)
    mkpath(dirname(outpng))
    save(outpng, fig)
    println("wrote $outpng")
    return
end

# --- clear-sky fluxes ----------------------------------------------------
clear = [
    ("OLR (TOA up)",       "olr_mod",       "olr_era",       "olr_bias"),
    ("Surface down LW",    "lwdn_sfc_mod",  "lwdn_sfc_era",  "lwdn_sfc_bias"),
    ("Surface net LW",     "lwnet_sfc_mod", "lwnet_sfc_era", "lwnet_sfc_bias"),
    ("TOA net SW",         "swnet_toa_mod", "swnet_toa_era", "swnet_toa_bias"),
    ("Surface down SW",    "swdn_sfc_mod",  "swdn_sfc_era",  "swdn_sfc_bias"),
    ("Surface net SW",     "swnet_sfc_mod", "swnet_sfc_era", "swnet_sfc_bias"),
]
make_figure(clear,
    "Clear-sky radiation: aeros operators on ERA5 columns vs ERA5 (annual mean, W m⁻²)",
    joinpath(outdir, "era5_rad_validation.png"))

# --- cloud radiative effect (TOA) ----------------------------------------
cre = [
    ("TOA LW cloud effect",  "cre_lw_toa_mod",  "cre_lw_toa_era",  "cre_lw_toa_bias"),
    ("TOA SW cloud effect",  "cre_sw_toa_mod",  "cre_sw_toa_era",  "cre_sw_toa_bias"),
    ("TOA net cloud effect", "cre_net_toa_mod", "cre_net_toa_era", "cre_net_toa_bias"),
]
make_figure(cre,
    "Cloud radiative effect: aeros cloudy operators on ERA5 columns vs ERA5 (annual mean, W m⁻²)",
    joinpath(outdir, "era5_cre_validation.png"); cmap = :balance)

# --- diagnosed-cloud CRE (aeros_cloud diagnosis on ERA5 columns) ----------
dcre = [
    ("TOA LW cloud effect",  "dcre_lw_toa_mod",  "cre_lw_toa_era",  "dcre_lw_toa_bias"),
    ("TOA SW cloud effect",  "dcre_sw_toa_mod",  "cre_sw_toa_era",  "dcre_sw_toa_bias"),
    ("TOA net cloud effect", "dcre_net_toa_mod", "cre_net_toa_era", "dcre_net_toa_bias"),
]
make_figure(dcre,
    "Diagnosed-cloud CRE: aeros_cloud diagnosis on ERA5 columns vs ERA5 (annual mean, W m⁻²)",
    joinpath(outdir, "era5_cre_diagnosed.png"); cmap = :balance)

close(ds)
