# b2b_common.jl -- NEW machinery for the B_multires stage-b2 APPROACH B probe.
#
# Approach A (b2_probe.jl / b_common.jl, committed, left untouched) diagnosed the
# resolution-error forcing from a SINGLE dynamical-core tendency evaluation at the
# T85 time-mean state X85 (ΔF = −RHS_T31(X85)). That was INSUFFICIENT: it removed
# only ~19% of the Z* error and made T*/U* worse, for two reasons -- (1) transient
# rectification (the instantaneous RHS at the time-mean state ≠ the time-mean RHS
# because of nonlinear eddy terms) and (2) it corrected only the dynamical core, so
# physics responded freely and reintroduced error.
#
# Approach B fixes both by diagnosing the constant forcing via NUDGING (the CAPT /
# Watt-Meyer "nudging tendency" object): run T31 with a relaxation term
# −(x − X85)/τ added to the prognostic tendencies, and accumulate the time-mean of
# that nudging tendency ⟨N⟩ over the averaging window. ⟨N⟩ is the constant forcing
# that holds free T31 at ≈X85, and because it is averaged over the actual fluctuating
# (eddying) trajectory rather than evaluated at the mean state it captures transient
# rectification; because it compensates the model's NET tendency each step it
# captures the physics response too.
#
# This file provides ONLY the new pieces (the nudging forcing + a builder that turns
# ⟨N⟩ into a plain ConstantTendencyForcing). It relies on b_common.jl for the
# ConstantTendencyForcing type, SpectralMeanAccumulator, set_state!/eval_dynamics!
# and the diagnostic helpers, which are `include`d first by b2b_probe.jl.
#
# SpeedyWeather v0.21.1. API read from the installed source.

using SpeedyWeather
using SpeedyWeather.SpeedyTransforms

# ============================================================================
# NudgingForcing -- state-dependent relaxation −(x − X85)/τ on the prognostic
# tendencies, injected at the SAME field-dependent points approach A uses:
#   * vorticity, divergence, pressure : added in SPECTRAL space (the model builds
#     those tendencies by accumulating, so a spectral term survives).
#   * temperature, humidity           : the model OVERWRITES the spectral tendency
#     with a grid→spectral transform of the GRID tendency, so the term must be
#     injected into the GRID tendency. We compute the nudge in spectral space and
#     transform it to grid, so after the model's grid→spectral round trip the
#     intended spectral term is recovered (verified to Float32 round-off, as in A).
#
# Units (CRITICAL -- verified against the SpeedyWeather source):
#   * The dynamical core integrates in a radius-scaled representation. The PROGNOSTIC
#     vorticity/divergence STATE is scaled by scale=R (scale_prognostic!), i.e.
#     ζ_state = R·ζ_phys; temperature/humidity/log-pₛ state are UNSCALED (physical).
#   * The TENDENCY arrays carry an EXTRA power of R on top of the state scaling:
#       - spectral vorticity/divergence tendency  = R²·(physical [1/s²])   (cf.
#         StochasticStirring, which scales its vorticity forcing by scale[]^2),
#       - grid temperature/humidity tendency       = R¹·(physical [K/s])    (cf.
#         HeldSuarez, which scales its relaxation frequency by radius),
#       - spectral log-pₛ tendency                 = R¹·(physical [1/s])    (winds
#         in surface_pressure_tendency! are R-scaled, lnpₛ gradients are not).
#   A proper relaxation −(x − X85)/τ therefore needs, in the tendency array units,
#   an extra factor of scale=R for EVERY field:
#       - vor/div: (X_state − x_state)/τ = R·(phys diff)/τ ; ×R -> R²·(phys)/τ  ✓
#       - temp/humid/pres: (X − x)/τ = (phys diff)/τ ; ×R -> R¹·(phys)/τ        ✓
#   Omitting this factor (an earlier bug) made the nudge ~R≈6.4e6× too weak: it
#   neither constrained the diagnosis run nor survived Float32 re-injection into the
#   large per-step tendencies (the free run came out byte-identical to bare). The
#   factor `progn.scale[]` (=R during a run) is applied uniformly below. X85 must be
#   in the same scaled representation as the runtime prognostic (it is: built from a
#   scaled-run mean), so the vor/div state difference is already correctly R-scaled.
#
# The forcing OPTIONALLY accumulates the applied nudge each step (accumulate flag,
# toggled by the driver after spinup) so ⟨N⟩ can be read off with mean_nudging().
# Only the fields selected by `correct_fields` are nudged (and accumulated); the
# rest evolve freely -- this is the §3.7 "correct selected terms, not everything".
mutable struct NudgingForcing{SV, SP} <: SpeedyWeather.AbstractForcing
    # target state X85 (spectral, radius-scaled representation)
    X_vor::SV; X_div::SV; X_temp::SV; X_humid::SV; X_pres::SP
    # scratch: last-applied spectral nudge per field (reused each step, no alloc)
    tmp_vor::SV; tmp_div::SV; tmp_temp::SV; tmp_humid::SV; tmp_pres::SP
    # accumulators for ⟨N⟩ (spectral)
    acc_vor::SV; acc_div::SV; acc_temp::SV; acc_humid::SV; acc_pres::SP
    inv_tau::Float64                 # 1/τ  [1/s]
    do_vor::Bool; do_div::Bool; do_temp::Bool; do_humid::Bool; do_pres::Bool
    active::Bool                     # gate: nudging applied only when true
    accumulate::Bool
    n::Int
end

# build from the target state X85 (NamedTuple of spectral fields) + timescale +
# the set of field names to correct. Starts INACTIVE (active=false) so the driver
# can do a plain, un-nudged spinup first (reach the T31 climate) before switching
# nudging on -- this avoids the cold-start shock of relaxing a near-rest default
# initial condition toward the fully-structured X85, which blows the run up.
function NudgingForcing(X85, tau_seconds::Real, correct_fields)
    cf = Set(String.(correct_fields))
    z(a) = zero(a)
    return NudgingForcing(
        X85.vor, X85.div, X85.temp, X85.humid, X85.pres,
        z(X85.vor), z(X85.div), z(X85.temp), z(X85.humid), z(X85.pres),
        z(X85.vor), z(X85.div), z(X85.temp), z(X85.humid), z(X85.pres),
        1.0 / Float64(tau_seconds),
        "vor" in cf, "div" in cf, "temp" in cf, "humid" in cf, "pres" in cf,
        false, false, 0,
    )
end

SpeedyWeather.initialize!(::NudgingForcing, ::SpeedyWeather.AbstractModel) = nothing

function SpeedyWeather.forcing!(
        vars,
        f::NudgingForcing,
        lf::Integer,
        model::SpeedyWeather.AbstractModel,
    )
    f.active || return nothing
    prog = vars.prognostic
    # PER-FIELD radius scaling (corrected — an earlier uniform ×R blew up temp/humid).
    # forcing! runs FIRST in dynamics_tendencies!; SpeedyWeather's scale_tendencies!
    # later multiplies the GRID temp/humid tendencies by R, but leaves the SPECTRAL
    # vor/div/pres tendencies alone. Therefore:
    #   * vor/div (spectral): coeff ×R on the already-R-scaled state diff -> R²·phys,
    #     matching the drag scheme (c = drag·scale[] on scaled vorticity).  s_spec.
    #   * pres (spectral, unscaled state): ×R -> R¹·phys.                    s_spec.
    #   * temp/humid: add in PHYSICAL units (no R here) -> scale_tendencies! ×R -> R¹. s_phys.
    s_spec = prog.scale[] * f.inv_tau
    s_phys = f.inv_tau

    if f.do_vor
        x = SpeedyWeather.get_step(prog.vorticity, lf)
        @. f.tmp_vor = (f.X_vor - x) * s_spec
        vars.tendencies.vorticity .+= f.tmp_vor
    end
    if f.do_div
        x = SpeedyWeather.get_step(prog.divergence, lf)
        @. f.tmp_div = (f.X_div - x) * s_spec
        vars.tendencies.divergence .+= f.tmp_div
    end
    if f.do_pres
        x = SpeedyWeather.get_step(prog.pressure, lf)
        @. f.tmp_pres = (f.X_pres - x) * s_spec
        vars.tendencies.pressure .+= f.tmp_pres
    end
    if f.do_temp
        x = SpeedyWeather.get_step(prog.temperature, lf)
        @. f.tmp_temp = (f.X_temp - x) * s_phys
        vars.tendencies.grid.temperature .+= transform(f.tmp_temp, model.spectral_transform)
    end
    if f.do_humid
        x = SpeedyWeather.get_step(prog.humidity, lf)
        @. f.tmp_humid = (f.X_humid - x) * s_phys
        vars.tendencies.grid.humidity .+= transform(f.tmp_humid, model.spectral_transform)
    end

    if f.accumulate
        f.n += 1
        f.do_vor   && (f.acc_vor   .+= f.tmp_vor)
        f.do_div   && (f.acc_div   .+= f.tmp_div)
        f.do_temp  && (f.acc_temp  .+= f.tmp_temp)
        f.do_humid && (f.acc_humid .+= f.tmp_humid)
        f.do_pres  && (f.acc_pres  .+= f.tmp_pres)
    end
    return nothing
end

# accumulation control (driver toggles these around the averaging window)
function reset_accumulation!(f::NudgingForcing)
    f.n = 0
    f.acc_vor   .= 0
    f.acc_div   .= 0
    f.acc_temp  .= 0
    f.acc_humid .= 0
    f.acc_pres  .= 0
    return f
end
set_accumulate!(f::NudgingForcing, b::Bool) = (f.accumulate = b; f)
set_active!(f::NudgingForcing, b::Bool) = (f.active = b; f)

# ⟨N⟩ as a NamedTuple of spectral fields (zero for non-corrected fields).
function mean_nudging(f::NudgingForcing)
    f.n > 0 || error("NudgingForcing: no nudging steps accumulated")
    n = f.n
    return (vor = f.acc_vor ./ n, div = f.acc_div ./ n, temp = f.acc_temp ./ n,
            humid = f.acc_humid ./ n, pres = f.acc_pres ./ n, nsteps = n)
end

# ============================================================================
# turn ⟨N⟩ into a plain ConstantTendencyForcing (from b_common.jl) for the free
# T31 + ⟨N⟩ run. NB: NO sign flip -- ⟨N⟩ is already the additive forcing that holds
# the state at X85 (⟨N⟩ ≈ −⟨model tendency⟩ at statistical steady state), whereas
# approach A's build_forcing takes ΔF = −RHS. temp/humid are transformed to grid
# space, matching the ConstantTendencyForcing injection convention.
function build_forcing_from_meanN(Nmean, model)
    S = model.spectral_transform
    dtemp_grid  = transform(Nmean.temp,  S)
    dhumid_grid = transform(Nmean.humid, S)
    return ConstantTendencyForcing(Nmean.vor, Nmean.div, Nmean.pres, dtemp_grid, dhumid_grid)
end
