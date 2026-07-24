# b_common.jl -- shared machinery for the B_multires (Option B) probes.
#
# This file is NEW for stage b2. It is deliberately standalone (it does NOT edit
# or depend on b1_probe.jl, which stays reproducible). b2_probe.jl `include`s it.
#
# It provides:
#   * diagnostic helpers (eddy, cos-lat RMS, regridding a native-grid field onto a
#     common full-Gaussian grid via spectral, spectral coarse-graining) -- the same
#     conventions as b1, duplicated here so b1 is left untouched.
#   * SpectralMeanAccumulator -- a callback that accumulates the time-mean of the
#     FULL prognostic spectral state (vorticity, divergence, temperature, humidity,
#     log-surface-pressure), which b1's grid-level accumulator did not do.
#   * ConstantTendencyForcing -- a custom SpeedyWeather AbstractForcing that injects
#     a precomputed, constant ΔF onto the prognostic tendencies every timestep. This
#     is the §3.7 "additive tendency forcing" (NOT a relaxation toward a state).
#
# SpeedyWeather v0.21.1. API read from the installed source, not old examples.

using SpeedyWeather
using SpeedyWeather.SpeedyTransforms
using SpeedyWeather.RingGrids
using Statistics

# ============================================================================
# 1. diagnostic helpers (same conventions as b1_probe.jl)
# ============================================================================

# native single-level field (on the model's own grid) -> common target grid,
# by a round trip through spectral space (no truncation beyond the source res).
function regrid_field(field, sg, target_grid)
    S_src = SpectralTransform(sg.spectrum, sg.grid; NF = Float32, nlayers = 1)
    spec  = transform(field, S_src)
    S_tgt = SpectralTransform(spec.spectrum, target_grid; NF = Float32, nlayers = 1)
    return transform(spec, S_tgt)
end

tomatrix(field, nlon, nlat) = reshape(Array(field), nlon, nlat)   # [lon, lat], lat 90->-90

# stationary-wave (eddy) = deviation from zonal mean; returns [lon,lat] matrix
function eddy(M)
    zm = mean(M; dims = 1)
    return M .- zm
end

# cos-lat area-weighted RMS over a [lon,nlat] matrix
function area_rms(M, latd)
    cw = cosd.(latd)
    num = 0.0; den = 0.0
    @inbounds for j in axes(M, 2)
        w = cw[j]
        s = 0.0
        for i in axes(M, 1)
            s += M[i, j]^2
        end
        num += w * s
        den += w * size(M, 1)
    end
    return sqrt(num / den)
end

# cos-lat area-weighted mean over a [lon,nlat] matrix
function area_mean(M, latd)
    cw = cosd.(latd)
    num = 0.0; den = 0.0
    @inbounds for j in axes(M, 2)
        w = cw[j]
        s = 0.0
        for i in axes(M, 1)
            s += M[i, j]
        end
        num += w * s
        den += w * size(M, 1)
    end
    return num / den
end

# location (lon,lat) and value of the largest |field| point
function argmax_loc(M, lond, latd)
    idx = argmax(abs.(M))
    i, j = Tuple(idx)
    return (lon = lond[i], lat = latd[j], value = M[i, j])
end

# top-N |field| maxima, coarsely de-duplicated
function top_maxima(M, lond, latd; n = 5, mindeg = 15.0)
    picks = Tuple{Float64, Float64, Float64}[]
    Awork = abs.(M)
    for _ in 1:n
        idx = argmax(Awork)
        i, j = Tuple(idx)
        push!(picks, (lond[i], latd[j], M[i, j]))
        for jj in axes(Awork, 2), ii in axes(Awork, 1)
            dlat = abs(latd[jj] - latd[j])
            dlon = abs(lond[ii] - lond[i]); dlon = min(dlon, 360 - dlon)
            if dlat < mindeg && dlon < mindeg
                Awork[ii, jj] = 0.0
            end
        end
    end
    return picks
end

# ============================================================================
# 2. time-mean accumulator for the FULL prognostic spectral state
# ============================================================================
# Accumulates the leapfrog-index-2 (current) spectral prognostic state each step.
# During a run the prognostic vorticity/divergence are radius-scaled (SpeedyWeather
# scales them for the dynamical core); we accumulate in that scaled representation
# and keep everything scaled consistently (radius is identical at both truncations),
# so no unscaling is needed here.
mutable struct SpectralMeanAccumulator <: SpeedyWeather.AbstractCallback
    n::Int
    vor::Any
    div::Any
    temp::Any
    humid::Any
    pres::Any
end
SpectralMeanAccumulator() = SpectralMeanAccumulator(0, nothing, nothing, nothing, nothing, nothing)

function SpeedyWeather.initialize!(cb::SpectralMeanAccumulator, vars, model)
    cb.n     = 0
    cb.vor   = zero(SpeedyWeather.get_step(vars.prognostic.vorticity,   2))
    cb.div   = zero(SpeedyWeather.get_step(vars.prognostic.divergence,  2))
    cb.temp  = zero(SpeedyWeather.get_step(vars.prognostic.temperature, 2))
    cb.humid = zero(SpeedyWeather.get_step(vars.prognostic.humidity,    2))
    cb.pres  = zero(SpeedyWeather.get_step(vars.prognostic.pressure,    2))
    return nothing
end

function SpeedyWeather.callback!(cb::SpectralMeanAccumulator, vars, model)
    cb.n += 1
    cb.vor   .+= SpeedyWeather.get_step(vars.prognostic.vorticity,   2)
    cb.div   .+= SpeedyWeather.get_step(vars.prognostic.divergence,  2)
    cb.temp  .+= SpeedyWeather.get_step(vars.prognostic.temperature, 2)
    cb.humid .+= SpeedyWeather.get_step(vars.prognostic.humidity,    2)
    cb.pres  .+= SpeedyWeather.get_step(vars.prognostic.pressure,    2)
    return nothing
end

SpeedyWeather.finalize!(::SpectralMeanAccumulator, args...) = nothing

# return the mean state as a NamedTuple of spectral fields (scaled representation)
function mean_state(cb::SpectralMeanAccumulator)
    cb.n > 0 || error("SpectralMeanAccumulator: no steps accumulated")
    n = cb.n
    return (vor = cb.vor ./ n, div = cb.div ./ n, temp = cb.temp ./ n,
            humid = cb.humid ./ n, pres = cb.pres ./ n, nsteps = n)
end

# ============================================================================
# 3. constant additive tendency forcing (§3.7 "diagnose the error operator")
# ============================================================================
# Injects a precomputed, constant ΔF onto the prognostic tendencies every step.
#
# SpeedyWeather assembles the spectral prognostic tendencies two different ways:
#   * vorticity, divergence, pressure : built by ACCUMULATING (add=true) onto the
#     spectral tendency array -> a spectral forcing added in forcing!() survives.
#   * temperature, humidity          : the spectral tendency is OVERWRITTEN by a
#     grid->spectral transform of the GRID tendency (temperature_tendency! line 610,
#     humidity_tendency!), so a spectral forcing would be wiped. These must be
#     injected into the GRID tendency (vars.tendencies.grid.{temperature,humidity}),
#     which vertical_advection! accumulates onto (-=) -- the same path HeldSuarez
#     uses. We store the grid-space forcing = itransform(ΔF_spectral).
#
# All fields are in the model's radius-scaled internal tendency units (the units in
# which the leapfrog integrates), so injecting them reproduces −RHS exactly.
struct ConstantTendencyForcing{SV, SP, GF} <: SpeedyWeather.AbstractForcing
    dvor::SV          # spectral (coeffs × layers)
    ddiv::SV          # spectral
    dpres::SP         # spectral (coeffs)  [l=m=0 mode is zeroed by the model anyway]
    dtemp_grid::GF    # grid field (points × layers)
    dhumid_grid::GF   # grid field
end

SpeedyWeather.initialize!(::ConstantTendencyForcing, ::SpeedyWeather.AbstractModel) = nothing

function SpeedyWeather.forcing!(
        vars,
        forcing::ConstantTendencyForcing,
        lf::Integer,
        model::SpeedyWeather.AbstractModel,
    )
    vars.tendencies.vorticity        .+= forcing.dvor
    vars.tendencies.divergence       .+= forcing.ddiv
    vars.tendencies.pressure         .+= forcing.dpres
    vars.tendencies.grid.temperature .+= forcing.dtemp_grid
    vars.tendencies.grid.humidity    .+= forcing.dhumid_grid
    return nothing
end

# build a ConstantTendencyForcing from a diagnosed spectral RHS (ΔF = −RHS)
function build_forcing(rhs, model)
    S = model.spectral_transform
    dtemp_grid  = transform(-rhs.temp,  S)   # spectral ΔF -> grid tendency
    dhumid_grid = transform(-rhs.humid, S)
    return ConstantTendencyForcing(-rhs.vor, -rhs.div, -rhs.pres, dtemp_grid, dhumid_grid)
end

# ============================================================================
# 4. set a prognostic state and evaluate the model's dynamical-core tendency
# ============================================================================
# Overwrite the prognostic spectral state (both leapfrog steps) with `X`.
function set_state!(vars, X)
    for lf in 1:2
        SpeedyWeather.get_step(vars.prognostic.vorticity,   lf) .= X.vor
        SpeedyWeather.get_step(vars.prognostic.divergence,  lf) .= X.div
        SpeedyWeather.get_step(vars.prognostic.temperature, lf) .= X.temp
        SpeedyWeather.get_step(vars.prognostic.humidity,    lf) .= X.humid
        SpeedyWeather.get_step(vars.prognostic.pressure,    lf) .= X.pres
    end
    return vars
end

# Evaluate the dynamical-core RHS at the current prognostic state (leapfrog index lf).
# Populates the grid diagnostics (geopotential, temperature, u, v from transform!;
# geopotential from geopotential! inside dynamics_tendencies!) and returns a copy of
# the resulting spectral prognostic tendencies. `model.forcing` (if any) is applied
# via forcing!() inside dynamics_tendencies!, exactly as in a real timestep.
function eval_dynamics!(vars, model, lf::Integer)
    SpeedyWeather.reset_tendencies!(vars)
    SpeedyTransforms.transform!(vars, lf, model)        # spectral state -> grid diagnostics
    SpeedyWeather.dynamics_tendencies!(vars, lf, model) # forcing! + dynamical core
    return (vor   = copy(vars.tendencies.vorticity),
            div   = copy(vars.tendencies.divergence),
            temp  = copy(vars.tendencies.temperature),
            humid = copy(vars.tendencies.humidity),
            pres  = copy(vars.tendencies.pressure))
end

# spectral coarse-graining: truncate a spectral field to degree `trunc_lo`
coarsegrain_spectral(field, trunc_lo) = SpeedyTransforms.spectral_truncation(field, trunc_lo)

# global (area-weighted) mean of a spectral field, per layer, from the l=m=0 harmonic.
# For a LowerTriangularArray the (l=0,m=0) coefficient is index 1; dividing by the
# sphere norm gives the area-weighted global mean. Returns a Vector over layers
# (or a scalar for a 2D/surface field).
function spectral_global_mean(field, S)
    nrm = S.norm_sphere
    if ndims(field) == 1 || size(field, 2) == 1
        return real(field[1]) / nrm
    else
        return [real(field[1, k]) / nrm for k in 1:size(field, 2)]
    end
end
