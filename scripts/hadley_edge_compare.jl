# hadley_edge_compare.jl -- side-by-side aeros vs SpeedyWeather cell-edge diagnostics.
#
# The residual subtropical moist bias traces to aeros's Hadley DESCENT not
# concentrating at ~30 deg. The candidate mechanism is dynamical: the subtropical
# jet and the eddy momentum flux [u*v*] that terminates the Hadley cell. This
# script loads both models' zonal-mean, time-mean structure and compares, at the
# same T21 aquaplanet footing:
#   omega(phi,sigma)     -- the overturning vertical branch (>0 subsidence)
#   Psi(phi,sigma)       -- mass streamfunction from [vbar] (the closed cell)
#   [ubar](phi,sigma)    -- the jet
#   [u*v*](phi,sigma)    -- eddy momentum flux, and its meridional convergence
#
# Usage: julia --project=mwm/C_omega scripts/hadley_edge_compare.jl \
#            <aeros_dump.nc> <speedy_omega.nc> [out.png]
#
# aeros stores fields (lev,lat); SpeedyWeather stores (lat,sigma). load2d() below
# normalizes both to (nlat, nlev) with sigma ascending top->surface.

using NCDatasets, Statistics, Printf, CairoMakie

const A_EARTH = 6.371e6      # m
const G0      = 9.81         # m/s2
const PS0     = 1.0e5        # Pa (aquaplanet reference surface pressure)

aeros_nc = length(ARGS) >= 1 ? ARGS[1] : error("need aeros dump nc")
speedy_nc = length(ARGS) >= 2 ? ARGS[2] : error("need speedy nc")
outpng   = length(ARGS) >= 3 ? ARGS[3] : joinpath(@__DIR__, "..", "docs", "figures", "hadley_edge_compare.png")

# --- load a 2d field oriented (nlat, nlev), given the lat vector to disambiguate ---
function load2d(ds, name, lat)
    haskey(ds, name) || return nothing
    M = Array(ds[name][:, :])
    nl = length(lat)
    size(M, 1) == nl && return M          # already (lat, lev)
    size(M, 2) == nl && return permutedims(M)  # (lev, lat) -> (lat, lev)
    error("$name shape $(size(M)) matches neither nlat=$nl")
end

# --- hemispheric symmetrization onto |lat| bins ---
function symmetrize(lat, field)
    absl = abs.(lat)
    bins = sort(unique(round.(absl; digits=1)))
    if ndims(field) == 1
        return bins, [mean(field[findall(x -> isapprox(round(abs(x); digits=1), b), lat)]) for b in bins]
    else
        nlev = size(field, 2)
        S = [mean(field[findall(x -> isapprox(round(abs(x); digits=1), b), lat), k]) for b in bins, k in 1:nlev]
        return bins, S
    end
end

# sigma half-level edges (0 at top, 1 at surface) reconstructed from full levels
function sigma_edges(sig)
    n = length(sig)
    e = zeros(n + 1)
    e[1] = 0.0; e[n+1] = 1.0
    for k in 2:n
        e[k] = 0.5 * (sig[k-1] + sig[k])
    end
    return e
end

# mass streamfunction Psi(phi,sigma) = (2*pi*a*cos(phi)*ps/g) * int_0^sigma [v] dsig'
# integrated from the model top downward; units 1e10 kg/s.
function streamfunction(lat, sig, vbar)
    e = sigma_edges(sig)
    dsig = diff(e)                              # thickness per full level
    nlat, nlev = size(vbar)
    Psi = zeros(nlat, nlev)
    for j in 1:nlat
        c = 2pi * A_EARTH * cosd(lat[j]) * PS0 / G0
        acc = 0.0
        for k in 1:nlev
            acc += vbar[j, k] * dsig[k]
            Psi[j, k] = c * acc
        end
    end
    return Psi ./ 1e10                          # -> 1e10 kg/s
end

# eddy momentum flux convergence -1/(a cos^2 phi) d/dphi([u*v*] cos^2 phi), m/s/day
function emf_convergence(lat, uv)
    nlat, nlev = size(uv)
    phi = deg2rad.(lat)
    cs2 = cosd.(lat) .^ 2
    F = uv .* cs2                               # [u*v*] cos^2 phi
    C = zeros(nlat, nlev)
    for k in 1:nlev
        for j in 2:nlat-1
            dFdphi = (F[j+1, k] - F[j-1, k]) / (phi[j+1] - phi[j-1])
            C[j, k] = -dFdphi / (A_EARTH * cs2[j])
        end
    end
    return C .* 86400.0                         # m/s per day
end

# ---------------------------------------------------------------------------
da = NCDataset(aeros_nc)
ds = NCDataset(speedy_nc)

lat_a = Array(da["lat"][:]); sig_a = Array(da["lev"][:])
lat_s = Array(ds["lat"][:]); sig_s = Array(ds["sigma"][:])

omega_a = load2d(da, "omega", lat_a); omega_s = load2d(ds, "omega", lat_s)
ubar_a  = load2d(da, "ubar",  lat_a); ubar_s  = load2d(ds, "uwind", lat_s)
vbar_a  = load2d(da, "vbar",  lat_a); vbar_s  = load2d(ds, "vwind", lat_s)
uv_a    = load2d(da, "uvpr",  lat_a); uv_s    = load2d(ds, "uvpr",  lat_s)
# aeros may fall back to snapshot u/v if the time-mean fields are absent
ubar_a === nothing && (ubar_a = load2d(da, "u", lat_a))
vbar_a === nothing && (vbar_a = load2d(da, "v", lat_a))

Psi_a = streamfunction(lat_a, sig_a, vbar_a)
Psi_s = streamfunction(lat_s, sig_s, vbar_s)
emfc_a = uv_a === nothing ? nothing : emf_convergence(lat_a, uv_a)
emfc_s = uv_s === nothing ? nothing : emf_convergence(lat_s, uv_s)

# --- headline numbers (symmetrized) ---
function report(tag, lat, sig, omega, ubar, uv, emfc)
    bins, O = symmetrize(lat, omega)
    _, U    = symmetrize(lat, ubar)
    println("\n=== $tag ===")
    # subtropical descent: max column subsidence in 15-40 band (levels 2:end-1)
    nlev = size(O, 2)
    cdesc = [maximum(O[i, 2:nlev-1]) for i in 1:length(bins)]
    sb = findall(b -> 12 <= b <= 45, bins)
    id = sb[argmax(cdesc[sb])]
    @printf("subtropical descent peak: omega=%+.2f hPa/day at |lat|=%.1f\n", cdesc[id], bins[id])
    # descent spread: fraction of 15-60 band with omega>1 at mid-sigma
    kmid = argmin(abs.(sig .- 0.5))
    desc_lats = bins[findall(i -> O[i, kmid] > 1.0 && bins[i] >= 12, 1:length(bins))]
    isempty(desc_lats) || @printf("  descent (omega>1 @ sig~0.5) spans |lat| %.0f..%.0f deg\n",
                                   minimum(desc_lats), maximum(desc_lats))
    # jet
    jmax = argmax(U); (jj, kk) = Tuple(jmax)
    @printf("jet [ubar] max: %.1f m/s at |lat|=%.1f, sigma=%.2f\n", U[jj, kk], bins[jj], sig[kk])
    if uv !== nothing
        _, UV = symmetrize(lat, uv)
        _, EC = symmetrize(lat, emfc)
        uvmax = maximum(abs.(UV))
        @printf("eddy mom flux |[u*v*]| max: %.2f m2/s2\n", uvmax)
        # subtropical (15-35) EMF convergence at mid-sigma
        sbj = findall(b -> 15 <= b <= 35, bins)
        @printf("subtropical EMF convergence (15-35, sig~0.5): %+.3f m/s/day\n", mean(EC[sbj, kmid]))
    else
        println("eddy mom flux: (absent)")
    end
end

report("aeros T21 RCE aquaplanet", lat_a, sig_a, omega_a, ubar_a, uv_a, emfc_a)
report("SpeedyWeather T21 aquaplanet", lat_s, sig_s, omega_s, ubar_s, uv_s, emfc_s)

# ---------------------------------------------------------------------------
# figure: rows = omega, Psi, [ubar], [u*v*]; cols = aeros, SpeedyWeather
# ---------------------------------------------------------------------------
mkpath(dirname(outpng))
fig = Figure(size = (1100, 1300))

function panel!(fig, row, col, lat, sig, F, title; colormap=:balance, symrange=nothing)
    ax = Axis(fig[row, col], title=title, xlabel="latitude", ylabel="sigma",
              yreversed=true)
    rng = symrange === nothing ? maximum(abs.(filter(isfinite, F))) : symrange
    rng = rng == 0 ? 1.0 : rng
    hm = heatmap!(ax, lat, sig, F; colormap=colormap, colorrange=(-rng, rng))
    Colorbar(fig[row, col+2], hm)
    return ax
end

# use a common color scale per row so aeros/SW are directly comparable
orng  = max(maximum(abs.(omega_a)), maximum(abs.(omega_s)))
prng  = max(maximum(abs.(Psi_a)),   maximum(abs.(Psi_s)))
urng  = max(maximum(abs.(ubar_a)),  maximum(abs.(ubar_s)))
uvrng = uv_s === nothing ? 1.0 : max(uv_a === nothing ? 0.0 : maximum(abs.(uv_a)), maximum(abs.(uv_s)))

panel!(fig, 1, 1, lat_a, sig_a, omega_a, "aeros  omega (hPa/day)"; symrange=orng)
panel!(fig, 1, 4, lat_s, sig_s, omega_s, "SpeedyWeather  omega (hPa/day)"; symrange=orng)
panel!(fig, 2, 1, lat_a, sig_a, Psi_a, "aeros  Psi (1e10 kg/s)"; symrange=prng)
panel!(fig, 2, 4, lat_s, sig_s, Psi_s, "SpeedyWeather  Psi (1e10 kg/s)"; symrange=prng)
panel!(fig, 3, 1, lat_a, sig_a, ubar_a, "aeros  [ubar] (m/s)"; symrange=urng)
panel!(fig, 3, 4, lat_s, sig_s, ubar_s, "SpeedyWeather  [ubar] (m/s)"; symrange=urng)
if uv_a !== nothing
    panel!(fig, 4, 1, lat_a, sig_a, uv_a, "aeros  [u*v*] (m2/s2)"; symrange=uvrng)
end
if uv_s !== nothing
    panel!(fig, 4, 4, lat_s, sig_s, uv_s, "SpeedyWeather  [u*v*] (m2/s2)"; symrange=uvrng)
end

save(outpng, fig)
println("\nwrote figure -> ", outpng)
close(da); close(ds)
