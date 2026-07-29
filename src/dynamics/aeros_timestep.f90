module aeros_timestep
    ! The time integrator: semi-implicit leapfrog, filtered, with hyperdiffusion.
    !
    ! This is the module that turns a tendency evaluation into a model. It owns
    ! the three time levels, the scratch the tendency needs, the semi-implicit
    ! solve and the diffusion, and it is the only place in aeros where a
    ! prognostic variable changes value.
    !
    ! === One step ============================================================
    !
    !   1. R      = aeros_tendency_calc(X^n)            the full nonlinear RHS
    !   2. zeta^(n+1) = zeta^(n-1) + h R_zeta           plain leapfrog
    !      (D, T, lnps)^(n+1) from aeros_semiimp_step   implicit gravity waves
    !   3. horizontal diffusion, implicitly, on zeta, D, T
    !   4. RAW time filter on X^n, using X^(n-1), X^n, X^(n+1)
    !   5. shift: (old, now) <- (now, new), by swapping descriptors
    !   6. global mass fixer, optional
    !
    ! h = 2 dt, except on the very first step of a run, where there is no
    ! X^(n-1) and the scheme starts as X^(n-1) = X^n with h = dt -- forward in
    ! the explicit terms, trapezoidal in the implicit ones. The semi-implicit
    ! factorization depends on h, so aeros_semiimp_set_step is called once at
    ! the transition and never again; the filter is skipped on that step,
    ! there being no third time level for it to act on.
    !
    ! === The time filter =====================================================
    !
    ! Leapfrog admits a computational mode that must be damped or it grows. The
    ! classical Robert-Asselin filter does that but is only FIRST-order
    ! accurate and damps the physical solution's mean as well -- which is the
    ! wrong property for a model whose purpose is 10^4-10^6 yr integrations
    ! (docs/design.md section 1).
    !
    ! So aeros uses Williams (2009)'s RAW variant, which splits the same
    ! displacement between the two time levels:
    !
    !     d       = (nu/2) ( X^(n-1) - 2 X^n + X^(n+1) )
    !     X^n     <- X^n     + alpha d
    !     X^(n+1) <- X^(n+1) - (1 - alpha) d
    !
    ! It costs three extra lines, restores third-order accuracy in the
    ! amplitude and conserves the mean exactly. alpha = 1 recovers the
    ! classical filter, which is what makes the change auditable rather than a
    ! reformulation: set raw_alpha = 1.0 in the namelist to get Robert-Asselin
    ! back. Williams' recommended 0.53 is the default.
    !
    ! === Horizontal diffusion ================================================
    !
    ! del^6 (namelist `ndiff`, an even order) on vorticity, divergence and
    ! temperature, applied IMPLICITLY -- an explicit hyperdiffusion at this
    ! order has a stability limit far shorter than the dynamics' own step, so
    ! doing it any other way would give back what the semi-implicit solve just
    ! bought. Diagonal in spectral space, so "implicitly" is one division:
    !
    !     X <- X / ( 1 + (h/tau_diff) [ l(l+1) / lmax(lmax+1) ]^(ndiff/2) )
    !
    ! written so that `tau_diff` is exactly the e-folding time of the LAST
    ! resolved wave, which is the only number in a hyperdiffusion anyone can
    ! reason about. It is a namelist parameter because it is a tuning knob, not
    ! a constant: the value that keeps a spectral core clean depends on
    ! truncation and on how vigorous the flow is.
    !
    ! NOT diffused: ln(p_s), which is a mass and diffusing it would leak one.
    !
    ! ANGULAR MOMENTUM: the operator damps l = 1 vorticity, i.e. solid-body
    ! rotation, which is the model's total angular momentum. The classical
    ! remedy is to diffuse (del^2 + 2/a^2) instead. It is not applied here
    ! because at del^6 the l = 1 factor is (2/(lmax(lmax+1)))^3 -- 8e-9 of the
    ! l = lmax rate at T31, an e-folding time of ~10^5 yr against a 6 h one --
    ! so the correction would change nothing measurable while making the
    ! operator no longer a pure power of the Laplacian. aeros_budget measures
    ! the drift rather than assuming it away.
    !
    ! === Adaptive hyperdiffusion (SpeedyWeather-style) =======================
    !
    ! Two OPT-IN refinements of the operator above, both default OFF and both
    ! bit-for-bit reducing to the fixed del^ndiff scheme when off. They are the
    ! two ideas SpeedyWeather uses to stay stable at the model top with no
    ! Rayleigh sponge (docs/refs/speedy_comparison.md 4a); aeros keeps its
    ! sponge, so these are a complementary, independently-useful numerical
    ! safety net against a runaway thermal wind or jet. Both stay IMPLICIT --
    ! forward-split diffusion on the leapfrog is unconditionally unstable
    ! (aeros_vdiff header, the house lesson) -- so each only reshapes the
    ! per-degree, now per-level, damping factor that the division applies.
    !
    ! (1) VORTICITY-SCALED STRENGTH (`diff_adapt`). Each step, per level, the
    ! diffusion coefficient is multiplied by
    !
    !     mult_k = 1 + diff_adapt_gain * max(0, |zeta|_k/diff_zeta_ref - 1)
    !
    ! capped at diff_adapt_max. |zeta|_k is the level's RMS vorticity [s-1],
    ! computed from the spectral coefficients by Parseval rather than a grid
    ! synthesis -- surface-mean-square = sum_lm w_m |vor_lm|^2 / (4 pi), with
    ! w_m = 1 for m = 0 and 2 for m > 0 (negative orders are not stored for a
    ! real field), and 4 pi = c00^2 read back from the transform. That is one
    ! O(nlm) reduction per level and no transform. RMS (not grid max) is the
    ! robust, cheap choice: a growing thermal wind / jet is a large-scale,
    ! low-degree structure, exactly what the RMS is dominated by, and it cannot
    ! be tripped by a single noisy gridpoint. Below the threshold mult_k = 1 and
    ! the operator is untouched; where and when a level spikes it is damped
    ! harder, in that level and that step only.
    !
    ! (2) SIGMA-TAPERED ORDER (`diff_taper`). The order is lowered toward the
    ! model top -- del^ndiff through the troposphere, ramping to del^diff_ndiff_top
    ! (a broader-scale, less scale-selective damping) above diff_taper_sigma.
    ! A lower order is a stronger, less selective sink at the top: an implicit
    ! sponge that cannot spin an unbounded thermal wind. This makes the
    ! per-degree ratio level-dependent, `dratio_lev(0:lmax, 1:nlev)`, with a
    ! per-level EVEN order (ramped in even steps so the exponent stays an integer
    ! power of the Laplacian -- and so diff_ndiff_top = ndiff is bit-for-bit the
    ! 1-D dratio). Precomputed at init and on any knob change; the per-step cost
    ! is zero.
    !
    ! === The mass fixer ======================================================
    !
    ! Optional (`mass_fixer`), off by default, applied LAST so that the state
    ! the caller receives is the fixed one.
    !
    ! docs/m1_results.md section 5.3 measures global mass drifting linearly at
    ! ~6.6e-6 per year under Held-Suarez -- 0.66 of the atmosphere over the
    ! 10^5 yr this model exists to run -- while the adiabatic tests hold it at
    ! 2e-16. The cause is not a bug and not the filter: halving dt halves the
    ! drift, so the per-step error is the scheme's ordinary O(dt^2) one, and
    ! what makes it accumulate rather than cancel is that the model integrates
    ! ln(p_s) while mass is int exp(ln p_s) dA. What the discretization
    ! conserves exactly is the CONTINUOUS statement.
    !
    ! The fix is exact and closed-form rather than iterative. Multiplying p_s
    ! by a constant r multiplies int p_s dA by exactly r, so
    !
    !     r = M_target / M ,      ln(p_s) <- ln(p_s) + ln(r)
    !
    ! and adding a constant to ln(p_s) changes the (0,0) spectral coefficient
    ! and NOTHING else. The fixer therefore cannot alter the surface-pressure
    ! pattern, only its global level -- which is the property that makes it
    ! safe to apply every step.
    !
    ! WHAT `lnr_cum` IS, AND WHAT IT IS NOT. Every correction is accumulated in
    ! `lnr_cum`, so the fixer's work is always visible rather than silent, and
    ! a fixer whose workload starts growing is a fixer that is covering for
    ! something -- watch it.
    !
    ! But it is NOT the drift an unfixed twin would have shown, and it must not
    ! be quoted as one. Measured, it runs 2-2.5x larger: 1.9e-5 against a
    ! 7.5e-6 unfixed drift on the adiabatic moving state over 200 steps
    ! (tests/test_timestep.f90), and -1.3e-5/yr against -6.6e-6/yr under
    ! Held-Suarez at T31L20 over 400 days.
    !
    ! The gap is not chaos -- the adiabatic state is not chaotic, and unfixed
    ! Held-Suarez trajectories seeded differently agree on their drift to
    ! under a percent. It is that the unfixed run's mass carries a bounded
    ! oscillation as well as a secular drift (the leapfrog computational mode
    ! in mass, docs/m1_results.md section 4), and over many steps that
    ! oscillation largely cancels in the accumulated total. The fixer resets it
    ! every step and therefore pays for it every step. Each individual
    ! `lnr_last` is still exactly the leak of that step from the fixed state --
    ! that is asserted to 1e-14 -- but the SUM is the fixer's workload, not the
    ! twin's drift.
    !
    ! So measuring the unfixed drift needs an unfixed twin run. That is the
    ! same discipline docs/design.md section 3.8 already requires of the
    ! correction terms, and for the same reason.
    !
    ! It is deliberately NOT a substitute for the deeper choice. Carrying p_s
    ! rather than ln(p_s) removes the nonlinearity at its root, at the cost of
    ! the reason ln(p_s) was chosen (the pressure-gradient term is linear in
    ! it). That remains open; this switch is what makes the two comparable,
    ! since `mass_fixer = .FALSE.` recovers the unfixed integrator exactly.
    !
    ! SIGMA-SURFACE DIFFUSION: diffusing T along model levels rather than
    ! pressure surfaces drives spurious heat fluxes over steep orography. It
    ! does not arise at M1 -- the surface is flat -- and the correction belongs
    ! with topography at M2, where docs/design.md section 4.2 already flags the
    ! ice-margin Gibbs problem it is entangled with.

    use aeros_defs,     only : dp, wp, wp_sh, io_unit_err, r_earth, grav, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_spec_class, aeros_state_class
    use aeros_spectral, only : aeros_sht_class, aeros_sht_pool_class, &
                                aeros_sht_pool_get, aeros_sht_synthesis, &
                                aeros_sht_analysis, aeros_sht_surface_integral, &
                                aeros_sht_lm
    use aeros_state,    only : aeros_spec_alloc, aeros_spec_end, aeros_spec_copy, &
                                aeros_spec_swap
    use aeros_vertical, only : aeros_vgrid_class
    use aeros_vordiv,   only : aeros_uv_from_vordiv, aeros_vordiv_from_uv
    use aeros_tendency, only : aeros_tend_class, aeros_work_class, &
                                aeros_tend_alloc, aeros_tend_end, &
                                aeros_work_alloc, aeros_work_alloc_diag, &
                                aeros_work_end, &
                                aeros_tendency_grid, aeros_tendency_spectral
    use aeros_semiimp,  only : aeros_semiimp_class, aeros_semiimp_init, &
                                aeros_semiimp_end, aeros_semiimp_set_step, &
                                aeros_semiimp_step, aeros_semiimp_print
    use aeros_correction, only : aeros_correction_class, aeros_correction_init, &
                                aeros_correction_load, aeros_correction_end, &
                                aeros_correction_apply, aeros_correction_report
    use aeros_transport, only : aeros_transport_class, aeros_transport_init, &
                                aeros_transport_end, aeros_transport_transport, &
                                aeros_transport_report
    use aeros_condensation, only : aeros_cond_class, aeros_condensation_init, &
                                aeros_condensation_load, aeros_condensation_end, &
                                aeros_condensation_apply, aeros_condensation_report
    use aeros_convection, only : aeros_conv_class, aeros_convection_init, &
                                aeros_convection_load, aeros_convection_end, &
                                aeros_convection_apply, aeros_convection_report
    use aeros_surface,    only : aeros_surf_class, aeros_surface_init, &
                                aeros_surface_load, aeros_surface_end, &
                                aeros_surface_apply, aeros_surface_report
    use aeros_radiation,  only : aeros_rad_class, aeros_radiation_init, &
                                aeros_radiation_load, aeros_radiation_end, &
                                aeros_radiation_apply, aeros_radiation_report
    use aeros_cloud_prog, only : aeros_cloud_prog_class, aeros_cloud_prog_init, &
                                aeros_cloud_prog_load, aeros_cloud_prog_end, &
                                aeros_cloud_prog_apply
    use aeros_vdiff,      only : aeros_vdiff_class, aeros_vdiff_init, &
                                aeros_vdiff_load, aeros_vdiff_end, &
                                aeros_vdiff_apply, aeros_vdiff_report
    use aeros_ocean,      only : aeros_ocean_class, aeros_ocean_init, &
                                aeros_ocean_load, aeros_ocean_end, &
                                aeros_ocean_step, aeros_ocean_report, &
                                aeros_ocean_albedo_update
    ! Land surface (feat/land): the second lower-boundary type. Composes the
    ! skin temperature/beta the surface+radiation see and steps the soil.
    use aeros_land,       only : aeros_land_class, aeros_land_init, &
                                aeros_land_load, aeros_land_end, &
                                aeros_land_pre, aeros_land_step, &
                                aeros_land_couple_radiation, aeros_land_report
    use aeros_held_suarez, only : aeros_hs_class, aeros_hs_init, aeros_hs_end, &
                                aeros_hs_apply, aeros_hs_print
    use ncio,             only : nc_create, nc_write, nc_write_dim, &
                                nc_write_attr, nc_read_attr, nc_read

    implicit none

    private

    type aeros_timestep_class
        ! Everything the integrator owns, allocated once and reused.

        real(wp) :: dt        = 0.0_wp    ! dynamics timestep [s]
        real(wp) :: eps_filter = 0.0_wp   ! Robert-Asselin coefficient nu [-]
        real(wp) :: raw_alpha  = 1.0_wp   ! Williams alpha; 1 = classical filter
        real(wp) :: tau_diff   = 0.0_wp   ! diffusion e-folding time at l=lmax [s]
        integer  :: ndiff      = 6        ! diffusion order (6 = del^6)
        logical  :: semi_implicit = .TRUE.

        ! Couple the diabatic heating (surface, convection, condensation,
        ! radiation) into the thermodynamic tendency BEFORE the semi-implicit
        ! solve (step 1b), the standard coupling, so the divergent circulation
        ! responds to the heating within the step. Default OFF = the validated
        ! forward-split increment on n+1 (step 6), bit-for-bit. The in-solve path
        ! needs a stronger time filter (eps_filter ~0.15) to hold the convective
        ! computational mode the centered heating excites; see m2_results 12.1.
        logical  :: couple_diabatic = .FALSE.

        ! Steps taken since init. Zero means the next one is the start-up step.
        integer  :: nstep = 0

        ! The other two time levels. `now` belongs to the caller.
        type(aeros_spec_class) :: old     ! X^(n-1)
        type(aeros_spec_class) :: new     ! X^(n+1), scratch between steps

        type(aeros_tend_class)    :: tnd
        type(aeros_work_class)    :: wrk
        type(aeros_semiimp_class) :: si

        ! Idealized forcing (M1.5). At M2 this becomes a physics driver with
        ! several components; the seam it is applied at does not change.
        type(aeros_hs_class)      :: hs

        ! Diagnostic prescribed heating [K/s], (nlon,nlat,nlev). When allocated,
        ! added to wrk%dtdt at step 1b -- the same in-solve seam the diabatic
        ! heating couples at -- independent of the physics toggles. A controlled
        ! analytic Q for isolating the dry dynamical core's heating->ω response
        ! from the moist physics. Unallocated (the default) is a bit-for-bit no-op.
        real(wp), allocatable :: q_force(:,:,:)

        ! The additive tendency correction (design.md 3.7, 3.8). A DIFFERENT
        ! seam from `hs` above -- spectral, not grid -- and the reason is in
        ! aeros_correction's header. Disabled unless a namelist says otherwise,
        ! so every test and every driver that does not ask for it gets the
        ! uncorrected model bit for bit.
        type(aeros_correction_class) :: cor

        ! Positive-definite grid tracer transport (aeros_transport). One engine
        ! instance, allocated once, advects every prognostic that lives off the
        ! spectral core: now%qv_g always (zero and staying zero for a dry run),
        ! and now%cf_g when the prognostic cloud scheme is on. Shared geometry
        ! and scratch; the fields are transported sequentially.
        type(aeros_transport_class) :: mst

        ! Large-scale condensation (aeros_condensation), at the grid seam. Off
        ! unless a namelist turns it on; a dry run never enters it.
        type(aeros_cond_class) :: cnd

        ! Moist convective adjustment (aeros_convection), at the grid seam and
        ! BEFORE condensation. Off unless a namelist turns it on.
        type(aeros_conv_class) :: cnv

        ! Surface energy/moisture budget (aeros_surface): prescribed-SST bulk
        ! turbulent fluxes into the lowest layer, at the grid seam BEFORE
        ! convection (it is the source convection overturns). Off unless a
        ! namelist turns it on; the source that closes the budget once
        ! Held-Suarez relaxation is removed.
        type(aeros_surf_class) :: surf

        ! Radiation (aeros_radiation): clear-sky LW+SW at the grid seam, cached
        ! and recomputed on a multi-hour cadence, heating added to wrk%dtdt.
        ! Uses the surface module's prescribed SST as the skin temperature. Off
        ! unless a namelist turns it on.
        type(aeros_rad_class) :: rad

        ! Prognostic cloud fraction (aeros_cloud_prog): the opt-in Sundqvist
        ! source/sink budget for now%cf_g, applied at the grid seam AFTER
        ! convection and condensation and BEFORE radiation, so radiation's
        ! all-sky path consumes the just-updated prognostic fraction. Off unless
        ! a namelist turns it on; when off, radiation uses the diagnostic
        ! RH->cover path and every current result is bit-for-bit unchanged.
        type(aeros_cloud_prog_class) :: cpr

        ! Sea surface temperature (aeros_ocean): prescribed or a slab. Owns the
        ! SST that the surface fluxes and radiation are evaluated against; in slab
        ! mode it closes the surface energy balance. Always initialized (radiation
        ! needs a skin temperature even when the surface fluxes are off).
        type(aeros_ocean_class) :: ocn

        ! Boundary-layer vertical diffusion (aeros_vdiff): mixes the surface
        ! source (aeros_surface) up the column, BETWEEN the surface fluxes and
        ! convection. Off unless a namelist turns it on. Without it the lowest
        ! layer builds a warm/moist inversion that convection turns into a
        ! grid-scale hot spot -- the M2 RCE instability (aeros_vdiff header).
        type(aeros_vdiff_class) :: vd

        ! === land state (feat/land) ==========================================
        ! Land surface as a second boundary type under the same flux seam as the
        ! ocean: a land-sea mask, a bucket soil moisture and a slab soil
        ! temperature (the land skin temperature). Disabled by default, so an
        ! all-ocean run is bit-for-bit and its prognostics are never allocated.
        type(aeros_land_class) :: land

        ! [ l(l+1) / lmax(lmax+1) ]^(ndiff/2), precomputed per degree. The
        ! step-dependent part is one multiply, so this never needs rebuilding.
        real(wp), allocatable :: dratio(:)   ! (0:lmax)

        ! === Adaptive hyperdiffusion (see the module header) =================
        ! Both opt-in, both default OFF => bit-for-bit the fixed del^ndiff scheme.
        !
        ! (1) Vorticity-scaled strength: mult_k on the coefficient, per level,
        !     keyed to the level's RMS vorticity exceeding diff_zeta_ref.
        logical  :: diff_adapt      = .FALSE.       ! vorticity-scaled strength
        real(wp) :: diff_zeta_ref   = 1.0e-4_wp     ! reference RMS |zeta| [s-1]
        real(wp) :: diff_adapt_gain = 1.0_wp        ! excess-vorticity gain [-]
        real(wp) :: diff_adapt_max  = 10.0_wp       ! cap on the multiplier [-]

        ! (2) Sigma-tapered order: order ramps from ndiff to diff_ndiff_top above
        !     diff_taper_sigma. dratio_lev is the per-level [l(l+1)/..]^(n_k/2);
        !     ndiff_lev is the per-level even order it was built with. Always
        !     allocated and built; when diff_taper is off every level uses ndiff,
        !     so dratio_lev(:,k) == dratio(:) exactly and the scheme is unchanged.
        logical  :: diff_taper       = .FALSE.      ! sigma-tapered order
        integer  :: diff_ndiff_top   = 4            ! order at the model top
        real(wp) :: diff_taper_sigma = 0.15_wp      ! ramp threshold [sigma]
        real(wp), allocatable :: dratio_lev(:,:)    ! (0:lmax, nlev)
        integer,  allocatable :: ndiff_lev(:)       ! (nlev)

        ! === Upper-level sponge ==============================================
        !
        ! A model-top damping layer, off by default. It removes the runaway a
        ! tropospheric lid (10 hPa, L12) develops once the artificial
        ! Held-Suarez stratospheric temperature cap is gone and clear-sky
        ! radiation cools the top layers without a strong balancing heating:
        ! the top over-cools, its equator-pole gradient drives an unbounded
        ! thermal wind, and the run blows up. Two pieces, both implicit and
        ! ramped over the top layers:
        !   C1  Rayleigh drag  -- vor,div relaxed toward zero (kills the winds)
        !   C2  Newtonian cool -- temp relaxed toward sponge_tref, and its
        !                         horizontal structure toward zero (removes the
        !                         gradient that drives the thermal wind)
        ! sponge_kr_lev / sponge_kt_lev are the per-level rates [s-1], nonzero
        ! only where sigma_full < sponge_sigma, squared-ramped to the top.
        logical  :: sponge_on    = .FALSE.
        real(wp) :: sponge_kr    = 1.0_wp/43200.0_wp   ! max Rayleigh rate [s-1] (0.5 d)
        real(wp) :: sponge_kt    = 1.0_wp/86400.0_wp   ! max Newtonian rate [s-1] (1 d)
        real(wp) :: sponge_sigma = 0.12_wp             ! sponge top of ramp [sigma]
        real(wp) :: sponge_tref  = 216.0_wp            ! reference strat T [K]
        real(wp), allocatable :: sponge_kr_lev(:)      ! (nlev) [s-1]
        real(wp), allocatable :: sponge_kt_lev(:)      ! (nlev) [s-1]

        ! === The mass fixer (see the module header) ==========================

        logical  :: mass_fixer = .FALSE.

        ! The mass the fixer restores to [kg]. Negative means "not yet set":
        ! the first fixed step captures the initial state's mass and corrects
        ! nothing. A RESTART must set it explicitly, with
        ! aeros_timestep_set_mass_target, or the target is recaptured from a
        ! state that has already drifted and the fixer silently ratifies the
        ! drift instead of removing it.
        real(dp) :: mass_target = -1.0_dp

        ! What the fixer has had to do. lnr_last is this step's log correction
        ! factor, and IS exactly that step's leak; lnr_cum is the running
        ! total, which is the fixer's cumulative workload and runs 2-2.5x an
        ! unfixed twin's drift -- see the module header before quoting it.
        real(dp) :: lnr_last = 0.0_dp
        real(dp) :: lnr_cum  = 0.0_dp

        ! The spectral coefficient of the constant field 1, read back from a
        ! transform at init rather than written down as sqrt(4 pi). It is the
        ! only number the fixer needs from SHTns' normalization convention,
        ! and reading it keeps the fixer correct if that convention ever
        ! changes -- the same reason aeros_grid reads its geometry back rather
        ! than recomputing it.
        real(dp) :: c00  = 0.0_dp
        integer  :: lm00 = 0            ! index of the (l=0, m=0) mode

        ! Grid scratch for the surface pressure the fixer integrates.
        real(dp), allocatable :: ps_fix(:,:)
    end type aeros_timestep_class

    public :: aeros_timestep_class
    public :: aeros_timestep_init
    public :: aeros_timestep_end
    public :: aeros_timestep_step
    public :: aeros_timestep_set_phis
    public :: aeros_timestep_set_qforce
    public :: aeros_timestep_set_mass_target
    public :: aeros_timestep_diagnose
    public :: aeros_timestep_enable_diag
    public :: aeros_timestep_set_sponge
    public :: aeros_timestep_set_diff_taper
    public :: aeros_timestep_print
    public :: aeros_timestep_write_restart
    public :: aeros_timestep_read_restart

contains

    subroutine aeros_timestep_init(ts, par, pool, grd, vg, filename, defaults_file)
        ! Allocate the time levels and the scratch, and build the solve.
        !
        ! `filename` is optional and only the correction layer needs it: its
        ! configuration is a namelist group of its own rather than a handful of
        ! scalars in aeros_param_class, because a term carries per-field and
        ! per-scale selectors and there may be several. Omitting it leaves the
        ! correction empty and disabled, which is what the acceptance tests
        ! want -- they measure the bare integrator.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_param_class),    intent(in)    :: par
        type(aeros_sht_pool_class), intent(in)    :: pool
        type(aeros_grid_class),     intent(in)    :: grd
        type(aeros_vgrid_class),    intent(in)    :: vg
        character(len=*), intent(in), optional    :: filename
        character(len=*), intent(in), optional    :: defaults_file

        integer  :: nlm, lmax, l, k
        real(wp) :: lref, ramp
        complex(wp_sh), allocatable :: c_one(:)
        real(dp), allocatable :: scr(:,:)

        nlm  = pool%sht(1)%nlm
        lmax = pool%sht(1)%lmax

        if (par%dt <= 0.0_wp) then
            write(io_unit_err,*) "aeros_timestep_init:: error: dt must be > 0, got ", par%dt
            error stop 1
        end if
        if (par%ndiff < 2 .or. mod(par%ndiff,2) /= 0) then
            write(io_unit_err,*) "aeros_timestep_init:: error: ndiff must be a positive even "// &
                                    "order (2 = del^2, 4 = del^4, 6 = del^6), got ", par%ndiff
            error stop 1
        end if
        if (par%tau_diff <= 0.0_wp) then
            write(io_unit_err,*) "aeros_timestep_init:: error: tau_diff must be > 0 h, got ", &
                                    par%tau_diff
            error stop 1
        end if
        if (par%eps_filter < 0.0_wp .or. par%eps_filter >= 1.0_wp) then
            write(io_unit_err,*) "aeros_timestep_init:: error: eps_filter must be in [0,1), got ", &
                                    par%eps_filter
            error stop 1
        end if
        if (par%raw_alpha <= 0.0_wp .or. par%raw_alpha > 1.0_wp) then
            write(io_unit_err,*) "aeros_timestep_init:: error: raw_alpha must be in (0,1], got ", &
                                    par%raw_alpha
            error stop 1
        end if

        call aeros_timestep_end(ts)

        ts%dt            = par%dt
        ts%eps_filter    = par%eps_filter
        ts%raw_alpha     = par%raw_alpha
        ts%tau_diff      = par%tau_diff*3600.0_wp      ! namelist is in hours
        ts%ndiff         = par%ndiff
        ts%semi_implicit = par%semi_implicit
        ts%mass_fixer    = par%mass_fixer
        ts%nstep         = 0

        ts%mass_target = -1.0_dp
        ts%lnr_last    = 0.0_dp
        ts%lnr_cum     = 0.0_dp

        call aeros_spec_alloc(ts%old, nlm, vg%nlev)
        call aeros_spec_alloc(ts%new, nlm, vg%nlev)
        call aeros_tend_alloc(ts%tnd, nlm, vg%nlev)
        call aeros_work_alloc(ts%wrk, grd%nlon, grd%nlat, vg%nlev)

        ! Factorized for the start-up step; switched to 2 dt after it.
        call aeros_semiimp_init(ts%si, vg, lmax, ts%dt)

        call aeros_hs_init(ts%hs, grd, par%held_suarez)

        if (present(filename)) then
            call aeros_correction_load(ts%cor, filename, pool%sht(1), vg%nlev, &
                                        defaults_file=defaults_file)
        else
            call aeros_correction_init(ts%cor, pool%sht(1), vg%nlev)
        end if

        call aeros_transport_init(ts%mst, grd, vg%nlev)

        if (present(filename)) then
            call aeros_condensation_load(ts%cnd, filename, grd, defaults_file=defaults_file)
            call aeros_convection_load(ts%cnv, filename, grd, defaults_file=defaults_file)
            call aeros_surface_load(ts%surf, filename, grd, defaults_file=defaults_file)
            call aeros_ocean_load(ts%ocn, filename, grd, defaults_file=defaults_file)
            call aeros_radiation_load(ts%rad, filename, grd, defaults_file=defaults_file)
            call aeros_cloud_prog_load(ts%cpr, filename, grd, defaults_file=defaults_file)
            call aeros_vdiff_load(ts%vd, filename, grd, defaults_file=defaults_file)
            call aeros_land_load(ts%land, filename, grd, defaults_file=defaults_file)   ! feat/land
        else
            call aeros_condensation_init(ts%cnd, grd, .FALSE.)
            call aeros_convection_init(ts%cnv, grd, .FALSE.)
            call aeros_surface_init(ts%surf, grd, .FALSE.)
            call aeros_ocean_init(ts%ocn, grd)
            call aeros_radiation_init(ts%rad, grd, .FALSE.)
            call aeros_cloud_prog_init(ts%cpr, grd, .FALSE.)
            call aeros_vdiff_init(ts%vd, grd, .FALSE.)
            call aeros_land_init(ts%land, grd, .FALSE.)   ! feat/land
        end if

        allocate(ts%dratio(0:lmax))
        lref = real(lmax*(lmax+1), wp)
        do l = 0, lmax
            ts%dratio(l) = (real(l*(l+1), wp)/lref)**(ts%ndiff/2)
        end do

        ! Per-level ratios for the sigma-tapered order (module header). Built
        ! from ts%ndiff / diff_taper / diff_ndiff_top / diff_taper_sigma, which
        ! carry their type defaults here (taper OFF); a driver that wants the
        ! taper sets those and calls aeros_timestep_set_diff_taper to rebuild.
        ! With taper off every level uses ndiff, so dratio_lev(:,k) == dratio(:).
        allocate(ts%dratio_lev(0:lmax, vg%nlev), ts%ndiff_lev(vg%nlev))
        call build_diff_ratio(ts, vg)

        ! The (0,0) mode index and the (0,0) coefficient of the constant field 1
        ! (= sqrt(4 pi) in this normalization), read back from a transform
        ! rather than written down, so a normalization change cannot silently
        ! break the mass fixer or the sponge's temperature reference. Always
        ! computed: both features need it, and it is a one-off transform.
        ts%lm00 = aeros_sht_lm(pool%sht(1), 0, 0)
        allocate(c_one(nlm), scr(grd%nlon, grd%nlat))
        scr = 1.0_dp
        call aeros_sht_analysis(pool%sht(1), scr, c_one)
        ts%c00 = real(c_one(ts%lm00), dp)
        deallocate(c_one, scr)
        if (abs(ts%c00) < 1.0e-12_dp) then
            write(io_unit_err,*) "aeros_timestep_init:: error: read a (0,0) coefficient of ", &
                                    ts%c00, " for the constant field 1, which cannot be right."
            error stop 1
        end if

        if (ts%mass_fixer) allocate(ts%ps_fix(grd%nlon, grd%nlat))

        ! Per-level sponge rates: squared ramp from zero at sigma_sponge to the
        ! maximum at the model top. Nonzero only in the top layers.
        allocate(ts%sponge_kr_lev(vg%nlev), ts%sponge_kt_lev(vg%nlev))
        call build_sponge_ramp(ts, vg)

        return

    end subroutine aeros_timestep_init

    subroutine build_sponge_ramp(ts, vg)
        ! (Re)build the per-level sponge rates from ts%sponge_kr/kt/sigma: a
        ! squared ramp from zero at sigma_sponge to the maximum at the model top,
        ! nonzero only where sigma_full < sponge_sigma. Split out of init so the
        ! strengths can be changed after init (aeros_timestep_set_sponge).

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_vgrid_class),    intent(in)    :: vg
        real(wp) :: ramp
        integer  :: k

        ts%sponge_kr_lev = 0.0_wp
        ts%sponge_kt_lev = 0.0_wp
        do k = 1, vg%nlev
            if (vg%sigma_full(k) < ts%sponge_sigma) then
                ramp = ((ts%sponge_sigma - vg%sigma_full(k))/ts%sponge_sigma)**2
                ts%sponge_kr_lev(k) = ts%sponge_kr*ramp
                ts%sponge_kt_lev(k) = ts%sponge_kt*ramp
            end if
        end do

        return

    end subroutine build_sponge_ramp

    subroutine aeros_timestep_set_sponge(ts, vg, kr, kt, sigma)
        ! Set the sponge strengths and rebuild the per-level ramp. For sweeping
        ! the model-top damping after init without recompiling (the top
        ! thermal-wind blow-up is the RCE's terminal event).

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_vgrid_class),    intent(in)    :: vg
        real(wp), intent(in) :: kr, kt, sigma

        ts%sponge_kr    = kr
        ts%sponge_kt    = kt
        ts%sponge_sigma = sigma
        call build_sponge_ramp(ts, vg)

        return

    end subroutine aeros_timestep_set_sponge

    subroutine build_diff_ratio(ts, vg)
        ! (Re)build the per-level hyperdiffusion ratios dratio_lev(0:lmax,nlev)
        ! and the per-level order ndiff_lev. When diff_taper is off, every level
        ! uses ts%ndiff and dratio_lev(:,k) reproduces the 1-D dratio(:) exactly.
        ! When on, the order ramps from ndiff (at/below diff_taper_sigma) to
        ! diff_ndiff_top at the model top, in EVEN steps so the exponent stays an
        ! integer power -- which is what makes diff_ndiff_top == ndiff a bit-for-
        ! bit no-op. sigma_full is (1:nlev), small at the top (k = 1). Split out
        ! of init so the taper can be reconfigured after init.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_vgrid_class),    intent(in)    :: vg

        integer  :: k, l, lmax, n
        real(wp) :: lref, frac, order_r

        lmax = ubound(ts%dratio, 1)
        lref = real(lmax*(lmax+1), wp)

        do k = 1, vg%nlev
            n = ts%ndiff
            if (ts%diff_taper .and. vg%sigma_full(k) < ts%diff_taper_sigma) then
                frac    = (ts%diff_taper_sigma - vg%sigma_full(k))/ts%diff_taper_sigma
                order_r = real(ts%ndiff, wp) &
                            - frac*real(ts%ndiff - ts%diff_ndiff_top, wp)
                n = 2*nint(order_r/2.0_wp)                  ! nearest even order
                n = max(ts%diff_ndiff_top, min(ts%ndiff, n))
            end if
            ts%ndiff_lev(k) = n
            do l = 0, lmax
                ts%dratio_lev(l,k) = (real(l*(l+1), wp)/lref)**(n/2)
            end do
        end do

        return

    end subroutine build_diff_ratio

    subroutine aeros_timestep_set_diff_taper(ts, vg, taper, ndiff_top, taper_sigma)
        ! Turn the sigma-tapered diffusion order on/off and (re)build the
        ! per-level ratios. For configuring the implicit top sponge after init
        ! without recompiling (the top thermal-wind blow-up is model-top
        ! dynamical). ndiff_top and taper_sigma are optional; omitted ones keep
        ! their current value.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_vgrid_class),    intent(in)    :: vg
        logical,  intent(in)           :: taper
        integer,  intent(in), optional :: ndiff_top
        real(wp), intent(in), optional :: taper_sigma

        ts%diff_taper = taper
        if (present(ndiff_top))   ts%diff_ndiff_top   = ndiff_top
        if (present(taper_sigma)) ts%diff_taper_sigma = taper_sigma

        if (ts%diff_ndiff_top < 2 .or. mod(ts%diff_ndiff_top,2) /= 0 &
            .or. ts%diff_ndiff_top > ts%ndiff) then
            write(io_unit_err,*) "aeros_timestep_set_diff_taper:: error: diff_ndiff_top "// &
                "must be a positive even order <= ndiff, got ", ts%diff_ndiff_top
            error stop 1
        end if

        call build_diff_ratio(ts, vg)

        return

    end subroutine aeros_timestep_set_diff_taper

    subroutine aeros_timestep_end(ts)

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts

        ts%nstep = 0

        call aeros_spec_end(ts%old)
        call aeros_spec_end(ts%new)
        call aeros_tend_end(ts%tnd)
        call aeros_work_end(ts%wrk)
        call aeros_semiimp_end(ts%si)
        call aeros_hs_end(ts%hs)
        call aeros_correction_end(ts%cor)
        call aeros_transport_end(ts%mst)
        call aeros_condensation_end(ts%cnd)
        call aeros_convection_end(ts%cnv)
        call aeros_surface_end(ts%surf)
        call aeros_ocean_end(ts%ocn)
        call aeros_radiation_end(ts%rad)
        call aeros_cloud_prog_end(ts%cpr)
        call aeros_vdiff_end(ts%vd)
        call aeros_land_end(ts%land)   ! feat/land

        if (allocated(ts%dratio)) deallocate(ts%dratio)
        if (allocated(ts%dratio_lev)) deallocate(ts%dratio_lev)
        if (allocated(ts%ndiff_lev)) deallocate(ts%ndiff_lev)
        if (allocated(ts%sponge_kr_lev)) deallocate(ts%sponge_kr_lev)
        if (allocated(ts%sponge_kt_lev)) deallocate(ts%sponge_kt_lev)
        if (allocated(ts%ps_fix)) deallocate(ts%ps_fix)

        ts%mass_target = -1.0_dp
        ts%lnr_last    = 0.0_dp
        ts%lnr_cum     = 0.0_dp

        return

    end subroutine aeros_timestep_end

    subroutine aeros_timestep_enable_diag(ts)
        ! Turn on the per-term heating diagnostics. Call AFTER aeros_timestep_init
        ! and after the physics toggles are set (the work arrays are already
        ! sized by init; this only adds the diagnostic split arrays). A no-op for
        ! any run that does not ask, so HS and production stay bit unchanged.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts

        call aeros_work_alloc_diag(ts%wrk)

        return

    end subroutine aeros_timestep_enable_diag

    subroutine aeros_timestep_step(ts, pool, vg, grd, now)
        ! Advance `now` by one dynamics timestep.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_vgrid_class),    intent(in)    :: vg
        type(aeros_grid_class),     intent(in)    :: grd
        type(aeros_state_class),    intent(inout) :: now

        type(aeros_sht_class), pointer :: s
        real(wp) :: h
        integer  :: k, lm

        s => pool%sht(1)

        ! Per-cell surface albedo -> radiation (rad%alb_map). Two schemes can set
        ! it and they COMPOSE (land and sea ice occupy disjoint cells -- ice never
        ! forms on land), so this is the one place both are folded into the single
        ! field radiation reads. Only entered when sea ice is on, because the ice
        ! albedo evolves and must be refreshed each step; the land-only map is
        ! static and built once at init (unchanged, bit-for-bit). When neither is
        ! on, rad%alb_map stays unallocated and radiation uses its scalar albedo.
        if (ts%ocn%l_seaice .and. allocated(ts%ocn%alb)) then
            if (ts%land%enabled) then
                ! Land albedo on land, ocean scalar on sea (rebuilt fresh so a
                ! melted-ice cell reverts to ocean albedo)...
                call aeros_land_couple_radiation(ts%land, ts%rad%alb_map, ts%rad%albedo)
                ! ...then the ice albedo on ice-covered (ocean) cells.
                where (ts%ocn%a_ice > 0.0_wp) ts%rad%alb_map = ts%ocn%alb
            else
                ! No land: ocn%alb already encodes ocean-vs-ice albedo everywhere.
                if (.not. allocated(ts%rad%alb_map)) &
                    allocate(ts%rad%alb_map(ts%ocn%nlon, ts%ocn%nlat))
                ts%rad%alb_map = ts%ocn%alb
            end if
        end if

        ! === 0. Step length, and the start-up ================================
        if (ts%nstep == 0) then
            ! No X^(n-1) yet: start from rest against the current state.
            h = ts%dt
            call aeros_spec_copy(ts%old, now%spec)

            ! Capture the fixer's target BEFORE any dynamics have run. Taking
            ! it after the first step instead would bake that step's error into
            ! the target and make it permanent.
            if (ts%mass_fixer .and. ts%mass_target < 0.0_dp) then
                ts%mass_target = mass_of(ts, s, now%spec)
            end if
        else
            h = 2.0_wp*ts%dt
        end if

        if (ts%si%dt_step /= h) call aeros_semiimp_set_step(ts%si, h)

        ! === 1. The right-hand sides, dynamics then physics ==================
        !
        ! Split at the grid-space seam so the forcing enters BEFORE the
        ! transforms back to spectral space -- which for a horizontally varying
        ! Rayleigh drag is a correctness requirement, not an optimization. See
        ! aeros_tendency_grid.
        call aeros_tendency_grid(pool, vg, grd, now%spec, ts%wrk)

        if (ts%hs%enabled) call aeros_hs_apply(ts%hs, vg, ts%wrk)

        ! Diabatic heating accumulator. Surface fluxes, convection, condensation
        ! and radiation each add a temperature increment [K] here over one dt.
        ! Default: applied forward onto n+1 at step 6. With couple_diabatic on:
        ! step 1b converts the sum to a rate and couples it into the thermodynamic
        ! tendency ahead of the semi-implicit solve, so the divergent circulation
        ! is generated in response to the heating within the step. Zero it once if
        ! any term is on.
        if (ts%surf%enabled .or. ts%cnv%enabled .or. ts%rad%enabled &
            .or. ts%cnd%enabled) &
            ts%wrk%dt_phys = 0.0_wp

        ! vdiff writes its diagnostic only when it runs; zero it here so a step
        ! with vdiff off reports 0 rather than the previous step's value.
        if (ts%wrk%diag) ts%wrk%dt_vdiff = 0.0_wp

        ! Land (feat/land): compose the skin temperature and evaporation
        ! efficiency the surface fluxes and radiation are driven against -- the
        ! ocean SST/beta=1 over sea, the prognostic soil temperature and the
        ! moisture-limited beta on land. A no-op when land is off, so the ocean
        ! path below is bit-for-bit.
        if (ts%land%enabled) call aeros_land_pre(ts%land, ts%ocn%sst)

        ! Surface turbulent fluxes, before convection: sensible heat warms the
        ! lowest layer (a forward increment to wrk%dt_phys) and evaporation
        ! moistens it (a forward increment to now%qv_g). This is the
        ! boundary-layer source that vertical diffusion then mixes up the column.
        ! Off unless enabled. With land on, the fluxes see the composed skin
        ! temperature and the beta evaporation efficiency; the ocean-only path is
        ! the untouched call.
        if (ts%surf%enabled) then
            if (ts%land%enabled) then
                call aeros_surface_apply(ts%surf, vg, ts%land%skin, ts%wrk%t_g, now%qv_g, &
                                            ts%wrk%lnps_g, ts%wrk%u, ts%wrk%v, &
                                            ts%wrk%dt_phys, ts%dt, evap_eff=ts%land%beta)
            else
                call aeros_surface_apply(ts%surf, vg, ts%ocn%sst, ts%wrk%t_g, now%qv_g, &
                                            ts%wrk%lnps_g, ts%wrk%u, ts%wrk%v, &
                                            ts%wrk%dt_phys, ts%dt)
            end if
        end if
        ! Diagnostic split: snapshot the CUMULATIVE dt_phys after each term (each
        ! only adds to it), differenced back into per-term increments after
        ! radiation. A disabled term leaves dt_phys unchanged, so its increment is
        ! zero -- no need to guard on the physics toggle.
        if (ts%wrk%diag) ts%wrk%dt_surf = ts%wrk%dt_phys

        ! Moist convective adjustment, before condensation: it overturns the
        ! column and precipitates, and large-scale condensation then removes
        ! whatever gridbox-mean supersaturation it left. Both act on the
        ! time-level-n gridpoint T and q. Its humidity change is applied to
        ! now%qv_g; its heating to wrk%dt_phys. See aeros_convection.
        if (ts%cnv%enabled) &
            call aeros_convection_apply(ts%cnv, vg, ts%wrk%t_g, now%qv_g, &
                                            ts%wrk%lnps_g, ts%wrk%dt_phys, ts%dt)
        if (ts%wrk%diag) ts%wrk%dt_cnv = ts%wrk%dt_phys

        ! Large-scale condensation, at the same grid seam: it dries the
        ! gridpoint humidity in place and adds its latent heating to wrk%dt_phys.
        ! Default (forward-split) keeps the heating and the forward drying on one
        ! discretization; with couple_diabatic the heating moves to the centered
        ! step-1b tendency while the drying stays forward on the off-spectral
        ! humidity, so column MSE is then conserved only in the mean (a watched
        ! budget diagnostic, not a per-step equality). Returns at once when dry.
        call aeros_condensation_apply(ts%cnd, vg, ts%wrk%t_g, now%qv_g, &
                                        ts%wrk%lnps_g, ts%wrk%dt_phys, ts%dt)
        if (ts%wrk%diag) ts%wrk%dt_cnd = ts%wrk%dt_phys

        ! Prognostic cloud fraction (aeros_cloud_prog), at the same grid seam and
        ! AFTER convection and condensation so it sees the grid-box saturation
        ! state the step ends in, and the convective precip for its detrainment
        ! source. Advances now%cf_g in place (genuine prognostic state, persisted
        ! across steps); radiation below consumes it. Off unless enabled -- when
        ! off, now%cf_g stays zero and radiation uses the diagnostic RH->cover.
        if (ts%cpr%enabled) then
            if (ts%cnv%enabled) then
                call aeros_cloud_prog_apply(ts%cpr, vg, ts%wrk%t_g, now%qv_g, &
                    ts%wrk%lnps_g, now%cf_g, ts%dt, precip=ts%cnv%precip)
            else
                call aeros_cloud_prog_apply(ts%cpr, vg, ts%wrk%t_g, now%qv_g, &
                    ts%wrk%lnps_g, now%cf_g, ts%dt)
            end if
        end if

        ! Radiation, at the same grid seam: clear-sky LW+SW heating, recomputed
        ! on a multi-hour cadence and cached, accumulated into wrk%dt_phys as a
        ! forward increment (NOT the centered path -- its top-of-atmosphere
        ! cooling is sharp and stiff). Uses the ocean's SST as the skin
        ! temperature (always allocated, even when the surface fluxes are off).
        ! When the prognostic cloud scheme is on it consumes now%cf_g; otherwise
        ! it diagnoses the cloud fraction from RH (bit-for-bit unchanged).
        ! Off unless enabled.
        ! Two independent options compose here (feat/land + feat/clouds):
        !   * skin temperature -- ts%land%skin (soil temp on land, SST over sea)
        !     when land is on, else ts%ocn%sst. Radiation also reads the land/ice
        !     albedo through rad%alb_map when it is allocated (set upstream).
        !   * cloud fraction -- the prognostic now%cf_g is passed as cf_prog when
        !     the Sundqvist scheme is on; otherwise radiation diagnoses cf from RH
        !     (bit-for-bit unchanged). Optional arg, so it is a 2x2 on (land,cpr).
        if (ts%rad%enabled) then
            if (ts%cpr%enabled) then
                if (ts%land%enabled) then
                    call aeros_radiation_apply(ts%rad, vg, grd, ts%wrk%t_g, now%qv_g, &
                                                ts%wrk%lnps_g, ts%land%skin, ts%nstep, &
                                                ts%dt, ts%wrk%dt_phys, cf_prog=now%cf_g)
                else
                    call aeros_radiation_apply(ts%rad, vg, grd, ts%wrk%t_g, now%qv_g, &
                                                ts%wrk%lnps_g, ts%ocn%sst, ts%nstep, &
                                                ts%dt, ts%wrk%dt_phys, cf_prog=now%cf_g)
                end if
            else
                if (ts%land%enabled) then
                    call aeros_radiation_apply(ts%rad, vg, grd, ts%wrk%t_g, now%qv_g, &
                                                ts%wrk%lnps_g, ts%land%skin, ts%nstep, &
                                                ts%dt, ts%wrk%dt_phys)
                else
                    call aeros_radiation_apply(ts%rad, vg, grd, ts%wrk%t_g, now%qv_g, &
                                                ts%wrk%lnps_g, ts%ocn%sst, ts%nstep, &
                                                ts%dt, ts%wrk%dt_phys)
                end if
            end if
        end if
        if (ts%wrk%diag) ts%wrk%dt_rad = ts%wrk%dt_phys
        ! Difference the cumulative snapshots into per-term increments [K/step].
        ! Top-down so each subtraction still sees the earlier cumulative value;
        ! dt_surf is already the increment (dt_phys started this step at zero).
        if (ts%wrk%diag) then
            ts%wrk%dt_rad = ts%wrk%dt_rad - ts%wrk%dt_cnd
            ts%wrk%dt_cnd = ts%wrk%dt_cnd - ts%wrk%dt_cnv
            ts%wrk%dt_cnv = ts%wrk%dt_cnv - ts%wrk%dt_surf
        end if

        ! Slab ocean: step the SST from the net surface energy flux just computed
        ! (surface turbulent fluxes + surface radiative fluxes against this SST).
        ! A no-op in prescribed mode, so the aquaplanet control is unchanged. This
        ! is what makes the surface flux self-limiting -- the fix for the local
        ! subtropical runaway a fixed SST cannot bound (aeros_ocean header).
        ! A no-op in prescribed mode (guarded inside).
        call aeros_ocean_step(ts%ocn, ts%rad%sw_net_sur, ts%rad%lw_dw_sur, &
                                ts%surf%shf, ts%surf%lhf, ts%dt)

        ! Land (feat/land): step the soil temperature and the bucket soil
        ! moisture from the fluxes just computed, over land points -- the soil
        ! analogue of the slab-ocean update above. Precipitation reaching the
        ! ground is the convective plus large-scale condensate. Guarded so the
        ! ocean run neither steps the soil nor forms the precip-sum temporary.
        if (ts%land%enabled) &
            call aeros_land_step(ts%land, ts%rad%sw_net_sur, ts%rad%lw_dw_sur, &
                                    ts%surf%shf, ts%surf%lhf, ts%surf%evap, &
                                    ts%cnv%precip + ts%cnd%precip, ts%dt)

        ! === 1b. Couple the diabatic heating into the thermodynamic tendency ==
        ! (couple_diabatic only; default OFF uses the step-6 forward-split.)
        ! Surface, convection, condensation and radiation accumulated a heating
        ! INCREMENT [K] over one dt in wrk%dt_phys. Convert it back to a RATE
        ! [K/s] -- dt_phys/dt -- and add it to the dynamical heating rate dtdt,
        ! so the whole diabatic forcing enters tnd%temp BEFORE the semi-implicit
        ! solve and the gravity-wave adjustment generates the divergent (Hadley)
        ! circulation in response to it, within the step. This is the standard
        ! coupling (IFS/SPEEDY/SpeedyWeather), and the treatment Held-Suarez
        ! Newtonian heating already gets here (aeros_held_suarez adds to dtdt).
        !
        ! The rate is the true physical rate, no leapfrog factor: the semi-
        ! implicit advances temperature over h = 2 dt but each step moves the
        ! clock forward only dt, so a rate Q deposits Q*dt per clock-step in the
        ! mean -- exactly what the forward-split (dt_phys = Q*dt onto n+1)
        ! deposits. Any future attenuation of a term stays an explicit factor on
        ! the rate, out here, never folded into the rate itself.
        !
        ! Guarded on the same physics-enabled condition that zeroed dt_phys, so
        ! the pure-dynamics / Held-Suarez path leaves dtdt untouched, bit-exact.
        if (ts%couple_diabatic .and. (ts%surf%enabled .or. ts%cnv%enabled &
            .or. ts%rad%enabled .or. ts%cnd%enabled)) &
            ts%wrk%dtdt = ts%wrk%dtdt + ts%wrk%dt_phys/ts%dt

        ! Diagnostic prescribed heating, on the same in-solve seam. A controlled
        ! analytic Q [K/s] for isolating the dry core's heating->ω response; a
        ! no-op unless a driver has set it (aeros_timestep_set_qforce).
        if (allocated(ts%q_force)) ts%wrk%dtdt = ts%wrk%dtdt + ts%q_force

        call aeros_tendency_spectral(pool, vg, ts%wrk, ts%tnd)

        ! The additive correction, on the assembled spectral tendency. Not at
        ! the grid seam with the physics above: see aeros_correction's header.
        ! Returns immediately when disabled, so the `DeltaF = 0` twin is
        ! bit-exact rather than approximately equal.
        call aeros_correction_apply(ts%cor, ts%tnd)

        ! === 2. Advance ======================================================
        !
        ! Vorticity has no linearized gravity-wave part, so it leapfrogs
        ! plainly whichever scheme the other three use.
        do k = 1, vg%nlev
            do lm = 1, s%nlm
                ts%new%vor(lm,k) = ts%old%vor(lm,k) + h*ts%tnd%vor(lm,k)
            end do
        end do

        if (ts%semi_implicit) then
            call aeros_semiimp_step(ts%si, s, ts%old, now%spec, ts%tnd, ts%new)
        else
            ! Fully explicit leapfrog. Not a production configuration -- it is
            ! limited by the ~300 m s-1 external gravity wave rather than by
            ! the wind -- but it is what makes the semi-implicit solve's value
            ! measurable rather than asserted (tests/test_timestep.f90).
            do k = 1, vg%nlev
                do lm = 1, s%nlm
                    ts%new%div(lm,k)  = ts%old%div(lm,k)  + h*ts%tnd%div(lm,k)
                    ts%new%temp(lm,k) = ts%old%temp(lm,k) + h*ts%tnd%temp(lm,k)
                end do
            end do
            do lm = 1, s%nlm
                ts%new%lnps(lm) = ts%old%lnps(lm) + h*ts%tnd%lnps(lm)
            end do
        end if

        ! === 3. Horizontal diffusion, implicit ===============================
        call diffuse(ts, s, h)

        ! === 3b. Upper-level sponge, implicit ================================
        if (ts%sponge_on) call sponge(ts, s, h)

        ! === 3c. Boundary-layer vertical diffusion, implicit =================
        ! Mixes the surface source (aeros_surface) up the column. Implicit on the
        ! n+1 state like the two diffusive operators above -- NOT forward-split,
        ! which for a diffusion term is unconditionally unstable on the leapfrog
        ! (aeros_vdiff header). Transforms the stepped T, winds and lnps to the
        ! grid, diffuses T/q/u/v there, and transforms back.
        if (ts%vd%enabled) call apply_vdiff(ts, s, vg, now)

        ! === 4. Time filter ==================================================
        if (ts%nstep > 0) call raw_filter(ts, s, now%spec, vg%nlev)

        ! === 5. Shift the time levels ========================================
        ! old <- now <- new, by exchanging descriptors: the retired X^(n-1)
        ! becomes the scratch that the next step writes X^(n+1) into.
        call aeros_spec_swap(ts%old, now%spec)
        call aeros_spec_swap(now%spec, ts%new)

        ! === 6. Forward-split physics heating (default; not couple_diabatic) ==
        ! Surface, convective, condensational AND radiative heating were
        ! accumulated on the grid in wrk%dt_phys as an increment [K]. Apply it
        ! FORWARD onto the n+1 state now -- transform per level and add to
        ! now%temp -- decoupled from the centered leapfrog, which would otherwise
        ! turn their large, vertically sharp heating into a computational-mode
        ! instability (m2_results 12.1). This mirrors the forward treatment of
        ! gridpoint humidity in step 8, so for condensation the drying and heating
        ! share one discretization. (Boundary-layer vertical diffusion is
        ! separate: a diffusion operator, applied implicitly in step 3c.) When
        ! couple_diabatic is on, step 1b coupled the heating in-solve instead and
        ! this is skipped.
        if (.not. ts%couple_diabatic .and. (ts%surf%enabled .or. ts%cnv%enabled &
            .or. ts%rad%enabled .or. ts%cnd%enabled)) &
            call apply_phys_heating(s, now%spec, ts%wrk%dt_phys, vg%nlev)

        ! === 7. Mass fixer ===================================================
        ! After the swap, so it acts on X^(n+1) -- the state the caller is
        ! about to receive -- and after the filter, so nothing follows it that
        ! could reintroduce what it just removed.
        if (ts%mass_fixer) call mass_fix(ts, s, now%spec)

        ! === 8. Grid tracer transport ========================================
        ! On the grid, forward in time, using the winds and surface pressure at
        ! the time level just stepped from -- which are still in ts%wrk, filled
        ! by aeros_tendency_grid and untouched since. The transported fields
        ! (now%qv_g, now%cf_g) persist across steps; nothing here is spectral.
        ! See aeros_transport.
        !
        ! Humidity first, as always. The engine reports its Courant/sub-step
        ! diagnostics from this (moisture) call, and both calls see identical
        ! winds so the numbers are the same regardless of order.
        call aeros_transport_transport(ts%mst, vg, ts%wrk%u, ts%wrk%v, &
                                        ts%wrk%lnps_g, now%qv_g, ts%dt)

        ! Cloud fraction, through the SAME engine instance (shared geometry and
        ! scratch), when the prognostic cloud scheme is on. cf starts in [0,1]
        ! and the scheme's max-principle keeps it there, so it transports safely
        ! with no clamp -- a clamp would break the conservative flux form.
        if (ts%cpr%enabled) &
            call aeros_transport_transport(ts%mst, vg, ts%wrk%u, ts%wrk%v, &
                                            ts%wrk%lnps_g, now%cf_g, ts%dt)

        ts%nstep = ts%nstep + 1

        return

    end subroutine aeros_timestep_step

    subroutine apply_phys_heating(s, spec, dt_phys, nlev)
        ! Add the forward-split physics temperature increment -- accumulated on
        ! the grid in [K] -- to the spectral temperature of the n+1 state, one
        ! analysis per level. The grid field is consumed (analysis overwrites
        ! its input). See step 6 of aeros_timestep_step for why the convective
        ! heating is applied here, forward, rather than through the centered RHS.
        ! Used only in the default (not couple_diabatic) configuration.

        implicit none

        type(aeros_sht_class),  intent(in)    :: s
        type(aeros_spec_class), intent(inout) :: spec
        real(wp),               intent(inout) :: dt_phys(:,:,:)
        integer,                intent(in)    :: nlev

        complex(wp_sh) :: tlm(s%nlm)
        integer :: k, lm

        do k = 1, nlev
            call aeros_sht_analysis(s, dt_phys(:,:,k), tlm)
            do lm = 1, s%nlm
                spec%temp(lm,k) = spec%temp(lm,k) + tlm(lm)
            end do
        end do

        return

    end subroutine apply_phys_heating

    subroutine apply_vdiff(ts, s, vg, now)
        ! Boundary-layer vertical diffusion, implicit on the n+1 state. The
        ! stepped spectral temperature, winds and surface pressure (ts%new) are
        ! synthesized to the grid, aeros_vdiff diffuses T, q, u and v there in
        ! place (a tridiagonal per column, unconditionally stable), and T and the
        ! winds are transformed back into ts%new. Humidity is a gridpoint field
        ! (now%qv_g) and is diffused directly. Called before the time filter and
        ! the swap, alongside the horizontal diffusion and sponge, because a
        ! diffusion operator must be implicit on the stepped state, not
        ! forward-split (aeros_vdiff header).

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_class),      intent(in)    :: s
        type(aeros_vgrid_class),    intent(in)    :: vg
        type(aeros_state_class),    intent(inout) :: now

        real(wp), allocatable :: t3(:,:,:), u3(:,:,:), v3(:,:,:), lnps2(:,:)
        complex(wp_sh) :: tlm(s%nlm), vor(s%nlm), div(s%nlm)
        integer :: k

        allocate(t3(s%nlon, s%nlat, vg%nlev), u3(s%nlon, s%nlat, vg%nlev), &
                 v3(s%nlon, s%nlat, vg%nlev), lnps2(s%nlon, s%nlat))

        do k = 1, vg%nlev
            call aeros_sht_synthesis(s, ts%new%temp(:,k), t3(:,:,k))
            call aeros_uv_from_vordiv(s, ts%new%vor(:,k), ts%new%div(:,k), &
                                        u3(:,:,k), v3(:,:,k))
        end do
        call aeros_sht_synthesis(s, ts%new%lnps, lnps2)

        ! Diagnostic: capture vdiff's grid-T change as t3(after) - t3(before).
        ! It rides t3 in place through aeros_vdiff_apply, so hold the pre-diffuse
        ! temperature and difference below. [K/step], like the per-term physics
        ! heating diagnostics (dt_surf/dt_cnv/dt_cnd/dt_rad).
        if (ts%wrk%diag) ts%wrk%dt_vdiff = t3

        call aeros_vdiff_apply(ts%vd, vg, t3, now%qv_g, u3, v3, lnps2, ts%dt, &
                                ts%surf%c_d, ts%surf%u_min)

        if (ts%wrk%diag) ts%wrk%dt_vdiff = t3 - ts%wrk%dt_vdiff

        do k = 1, vg%nlev
            call aeros_sht_analysis(s, t3(:,:,k), tlm)
            ts%new%temp(:,k) = tlm
            call aeros_vordiv_from_uv(s, u3(:,:,k), v3(:,:,k), vor, div)
            ts%new%vor(:,k) = vor
            ts%new%div(:,k) = div
        end do

        deallocate(t3, u3, v3, lnps2)

        return

    end subroutine apply_vdiff

    subroutine mass_fix(ts, sht, spec)
        ! Restore the global mass by the one rescaling of p_s that does it
        ! exactly. See the module header for why this exists and what it is
        ! deliberately not.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_class),      intent(in)    :: sht
        type(aeros_spec_class),     intent(inout) :: spec

        real(dp) :: lnr

        lnr = log(ts%mass_target/mass_of(ts, sht, spec))

        ! p_s <- r p_s, which in ln(p_s) is the addition of a constant, which
        ! is a change to the (0,0) coefficient alone. Every other wavenumber is
        ! untouched, so the surface-pressure pattern is exactly what the
        ! dynamics produced.
        spec%lnps(ts%lm00) = spec%lnps(ts%lm00) + cmplx(lnr*ts%c00, 0.0_dp, wp_sh)

        ts%lnr_last = lnr
        ts%lnr_cum  = ts%lnr_cum + lnr

        return

    end subroutine mass_fix

    real(dp) function mass_of(ts, sht, spec) result(mass)
        ! Global dry-air mass [kg] of a spectral state.
        !
        ! Deliberately the SAME quadrature aeros_budget uses: a fixer that
        ! measured mass differently from the diagnostic would drive the
        ! diagnostic somewhere other than zero, and the disagreement would look
        ! like a leak.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_class),      intent(in)    :: sht
        type(aeros_spec_class),     intent(in)    :: spec

        ! ln(p_s) -> p_s on the grid. One single-level synthesis, against the
        ! ~11 per level the step has already done.
        call aeros_sht_synthesis(sht, spec%lnps, ts%ps_fix)
        ts%ps_fix = exp(ts%ps_fix)

        mass = aeros_sht_surface_integral(sht, ts%ps_fix)/grav

        return

    end function mass_of

    subroutine aeros_timestep_set_mass_target(ts, mass)
        ! Set the mass the fixer restores to, explicitly.
        !
        ! Exists for restarts. Without it a restarted run recaptures the target
        ! from its own initial state, which is the already-drifted state the
        ! previous run ended on -- so the fixer would hold the drift rather
        ! than remove it, and would do so silently.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        real(dp),                   intent(in)    :: mass

        if (mass <= 0.0_dp) then
            write(io_unit_err,*) "aeros_timestep_set_mass_target:: error: mass must be > 0, got ", &
                                    mass
            error stop 1
        end if

        ts%mass_target = mass

        return

    end subroutine aeros_timestep_set_mass_target

    subroutine aeros_timestep_write_restart(ts, now, time, filename)
        ! Write the full integrator state to a self-describing netCDF restart
        ! file, from which aeros_timestep_read_restart resumes bit-identically.
        !
        ! WHAT IS SAVED, and why exactly this and no more:
        !   - BOTH persistent spectral time levels -- X^n (now%spec) and
        !     X^(n-1) (ts%old). Leapfrog needs n and n-1 to take its 2 dt step;
        !     ts%new is scratch rewritten every step, so it is NOT saved.
        !   - now%qv_g, the gridpoint humidity prognostic. It is never
        !     transformed to spectral space (aeros_defs), so there is nothing
        !     else to reconstruct it from -- it MUST be saved.
        !   - the ocean's SST and net-flux fields and its mode: the spun-up
        !     lower boundary the restart exists to reuse.
        !   - ts%nstep and the model `time`. nstep drives the leapfrog start-up
        !     branch AND the radiation recompute cadence (mod(nstep,nrad)), so
        !     restoring it is what makes the cached radiative heating correct.
        !   - mass_target / lnr_cum / lnr_last. mass_target MUST be carried:
        !     recapturing it from the restarted state would ratify the drift
        !     (see aeros_timestep_class and aeros_timestep_set_mass_target).
        !   - ts%rad%heat, the cached heating applied every step between
        !     recomputes. The radiation cadence has no internal "last recompute"
        !     accumulator -- it is purely mod(nstep,nrad) -- so nstep + this
        !     cache are together sufficient for a bit-exact first restarted step.
        !   - the cached radiative surface/TOA fluxes (sw_net_sur, lw_dw_sur,
        !     sw_dw_sur, olr, sw_up_toa). These, like rad%heat, are refreshed
        !     ONLY on a recompute step and held between; the slab ocean consumes
        !     sw_net_sur and lw_dw_sur every step (aeros_ocean_step), so a
        !     restart landing on a non-recompute step must carry them or the
        !     first SST update -- and everything downstream of it -- diverges.
        !
        ! Grid diagnostics (now%u/v/temp_g/ps) are recomputed from spectral each
        ! step and are deliberately NOT saved.
        !
        ! Complex spectral arrays are stored as separate real/imag components.
        ! Metadata (nlm, nlev, nlon, nlat) goes in as global attributes and the
        ! read routine errors on any mismatch.

        implicit none

        type(aeros_timestep_class), intent(in) :: ts
        type(aeros_state_class),    intent(in) :: now
        real(wp),                   intent(in) :: time
        character(len=*),           intent(in) :: filename

        integer :: nlm, nlev, nlon, nlat, rad_present, land_present, cf_present, seaice_present

        nlm  = now%spec%nlm
        nlev = now%spec%nlev
        nlon = ts%ocn%nlon
        nlat = ts%ocn%nlat
        rad_present = 0
        if (allocated(ts%rad%heat)) rad_present = 1
        land_present = 0
        if (ts%land%enabled .and. allocated(ts%land%w)) land_present = 1   ! feat/land
        cf_present = 0
        if (ts%cpr%enabled) cf_present = 1                                 ! feat/clouds
        seaice_present = 0
        if (ts%ocn%l_seaice .and. allocated(ts%ocn%h_ice)) seaice_present = 1  ! feat/seaice

        call nc_create(filename)

        ! Plain index dimensions; the restart is not a viewable field file.
        call nc_write_dim(filename,"lon",x=1,dx=1,nx=nlon,units="1")
        call nc_write_dim(filename,"lat",x=1,dx=1,nx=nlat,units="1")
        call nc_write_dim(filename,"lev",x=1,dx=1,nx=nlev,units="1")
        call nc_write_dim(filename,"nlm",x=1,dx=1,nx=nlm,units="1")

        ! Validation metadata and scalars, as global attributes.
        call nc_write_attr(filename,"nlm",              nlm)
        call nc_write_attr(filename,"nlev",             nlev)
        call nc_write_attr(filename,"nlon",             nlon)
        call nc_write_attr(filename,"nlat",             nlat)
        call nc_write_attr(filename,"nstep",            ts%nstep)
        call nc_write_attr(filename,"ocn_mode",         ts%ocn%mode)
        call nc_write_attr(filename,"rad_heat_present", rad_present)
        call nc_write_attr(filename,"land_present",     land_present)    ! feat/land
        call nc_write_attr(filename,"seaice_present",   seaice_present)  ! feat/seaice
        call nc_write_attr(filename,"time",             real(time,dp))
        call nc_write_attr(filename,"mass_target",      ts%mass_target)
        call nc_write_attr(filename,"lnr_cum",          ts%lnr_cum)
        call nc_write_attr(filename,"lnr_last",         ts%lnr_last)

        ! The two persistent spectral time levels.
        call write_spec("now", now%spec)   ! X^n
        call write_spec("old", ts%old)     ! X^(n-1)

        ! Gridpoint humidity prognostic.
        call nc_write(filename,"qv_g", now%qv_g, dim1="lon",dim2="lat",dim3="lev", &
                        units="kg kg-1", long_name="specific humidity (prognostic)")

        ! Ocean lower boundary.
        call nc_write(filename,"ocn_sst", ts%ocn%sst, dim1="lon",dim2="lat", &
                        units="K", long_name="sea surface temperature")
        call nc_write(filename,"ocn_fnet", ts%ocn%fnet, dim1="lon",dim2="lat", &
                        units="W m-2", long_name="net surface energy flux")

        ! === sea-ice state (feat/seaice) ===
        ! Prognostic ice thickness / fraction / surface temperature, saved only
        ! when the thermodynamic sea ice is active (seaice_present). Append-only:
        ! a restart without sea ice carries neither the attribute's arrays nor any
        ! change to the blocks above. The surface albedo field is a diagnostic of
        ! these three (aeros_ocean_albedo_update) and is rebuilt on read, not saved.
        if (seaice_present == 1) then
            call nc_write(filename,"ocn_h_ice", ts%ocn%h_ice, dim1="lon",dim2="lat", &
                            units="m", long_name="sea-ice thickness")
            call nc_write(filename,"ocn_a_ice", ts%ocn%a_ice, dim1="lon",dim2="lat", &
                            units="1", long_name="sea-ice fraction")
            call nc_write(filename,"ocn_t_ice_sfc", ts%ocn%t_ice_sfc, dim1="lon",dim2="lat", &
                            units="K", long_name="sea-ice surface temperature")
        end if

        ! Cached radiative heating (only if it has been computed).
        if (rad_present == 1) &
            call nc_write(filename,"rad_heat", ts%rad%heat, &
                            dim1="lon",dim2="lat",dim3="lev", &
                            units="K s-1", long_name="cached radiative heating rate")

        ! Cached radiative surface/TOA fluxes, refreshed only on a recompute
        ! and consumed between recomputes (the slab ocean reads sw_net_sur and
        ! lw_dw_sur every step). Always allocated once radiation is initialised.
        call nc_write(filename,"rad_sw_net_sur", ts%rad%sw_net_sur, &
                        dim1="lon",dim2="lat",units="W m-2")
        call nc_write(filename,"rad_lw_dw_sur",  ts%rad%lw_dw_sur, &
                        dim1="lon",dim2="lat",units="W m-2")
        call nc_write(filename,"rad_sw_dw_sur",  ts%rad%sw_dw_sur, &
                        dim1="lon",dim2="lat",units="W m-2")
        call nc_write(filename,"rad_olr",        ts%rad%olr, &
                        dim1="lon",dim2="lat",units="W m-2")
        call nc_write(filename,"rad_sw_up_toa",  ts%rad%sw_up_toa, &
                        dim1="lon",dim2="lat",units="W m-2")

        ! === land state (feat/land) ==========================================
        ! The two land prognostics -- bucket soil moisture and slab soil
        ! temperature -- are new state with no spectral counterpart; a restart
        ! must carry them or a resumed land run diverges silently. Written only
        ! when land is on; the mask and albedo maps are static and rebuilt from
        ! file on init, so they are not saved. (land_present is written up top
        ! with the other scalar attributes.)
        if (land_present == 1) then
            call nc_write(filename,"land_w", ts%land%w, dim1="lon",dim2="lat", &
                            units="m", long_name="bucket soil moisture")
            call nc_write(filename,"land_t_soil", ts%land%t_soil, dim1="lon",dim2="lat", &
                            units="K", long_name="slab soil temperature")
        end if
        ! === prognostic cloud state (feat/clouds) ===========================
        ! now%cf_g is new prognostic state (aeros_cloud_prog): a gridpoint field
        ! never reconstructed from anything else, so it MUST be saved when the
        ! scheme is on or the run diverges on restart. Gated by cf_present so a
        ! diagnostic-cloud run writes only the flag (=0) and no field, leaving
        ! its restart otherwise unchanged. Append-only block; see rad_present.
        call nc_write_attr(filename,"cf_present", cf_present)
        if (cf_present == 1) &
            call nc_write(filename,"cf_g", now%cf_g, &
                            dim1="lon",dim2="lat",dim3="lev", &
                            units="1", long_name="prognostic cloud fraction")

        return

    contains

        subroutine write_spec(prefix, spec)
            ! Write one spectral level set as real/imag component variables.
            character(len=*),       intent(in) :: prefix
            type(aeros_spec_class), intent(in) :: spec

            call nc_write(filename, prefix//"_vor_re",  real(spec%vor,dp), &
                            dim1="nlm",dim2="lev")
            call nc_write(filename, prefix//"_vor_im",  aimag(spec%vor), &
                            dim1="nlm",dim2="lev")
            call nc_write(filename, prefix//"_div_re",  real(spec%div,dp), &
                            dim1="nlm",dim2="lev")
            call nc_write(filename, prefix//"_div_im",  aimag(spec%div), &
                            dim1="nlm",dim2="lev")
            call nc_write(filename, prefix//"_temp_re", real(spec%temp,dp), &
                            dim1="nlm",dim2="lev")
            call nc_write(filename, prefix//"_temp_im", aimag(spec%temp), &
                            dim1="nlm",dim2="lev")
            call nc_write(filename, prefix//"_lnps_re", real(spec%lnps,dp), &
                            dim1="nlm")
            call nc_write(filename, prefix//"_lnps_im", aimag(spec%lnps), &
                            dim1="nlm")
            return
        end subroutine write_spec

    end subroutine aeros_timestep_write_restart

    subroutine aeros_timestep_read_restart(ts, now, time, filename)
        ! Overwrite an already-initialised timestep and state with a restart
        ! file written by aeros_timestep_write_restart, resuming bit-identically.
        !
        ! CONTRACT: the caller must have built `ts` and `now` normally first
        ! (same trunc/nlev/grid), so every array is already allocated to the
        ! right shape; this routine only overwrites values. The metadata in the
        ! file is validated against the live shapes and any mismatch is fatal.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_state_class),    intent(inout) :: now
        real(wp),                   intent(out)   :: time
        character(len=*),           intent(in)    :: filename

        integer  :: nlm, nlev, nlon, nlat, rad_present, land_present, cf_present, seaice_present
        integer  :: nlm_f, nlev_f, nlon_f, nlat_f
        real(dp) :: time_f

        nlm  = now%spec%nlm
        nlev = now%spec%nlev
        nlon = ts%ocn%nlon
        nlat = ts%ocn%nlat

        call nc_read_attr(filename,"nlm",  nlm_f)
        call nc_read_attr(filename,"nlev", nlev_f)
        call nc_read_attr(filename,"nlon", nlon_f)
        call nc_read_attr(filename,"nlat", nlat_f)

        if (nlm_f /= nlm .or. nlev_f /= nlev .or. &
            nlon_f /= nlon .or. nlat_f /= nlat) then
            write(io_unit_err,*) "aeros_timestep_read_restart:: error: restart file "// &
                                    "metadata does not match the current model."
            write(io_unit_err,*) "  file  (nlm,nlev,nlon,nlat) = ", &
                                    nlm_f, nlev_f, nlon_f, nlat_f
            write(io_unit_err,*) "  model (nlm,nlev,nlon,nlat) = ", &
                                    nlm, nlev, nlon, nlat
            error stop 1
        end if

        call nc_read_attr(filename,"nstep",            ts%nstep)
        call nc_read_attr(filename,"ocn_mode",         ts%ocn%mode)
        call nc_read_attr(filename,"rad_heat_present", rad_present)
        call nc_read_attr(filename,"time",             time_f)
        call nc_read_attr(filename,"mass_target",      ts%mass_target)
        call nc_read_attr(filename,"lnr_cum",          ts%lnr_cum)
        call nc_read_attr(filename,"lnr_last",         ts%lnr_last)

        time = real(time_f, wp)

        call read_spec("now", now%spec)   ! X^n
        call read_spec("old", ts%old)     ! X^(n-1)

        call nc_read(filename,"qv_g",     now%qv_g)
        call nc_read(filename,"ocn_sst",  ts%ocn%sst)
        call nc_read(filename,"ocn_fnet", ts%ocn%fnet)

        ! === sea-ice state (feat/seaice) ===
        ! Restore the ice prognostics when the file carries them (older restarts,
        ! and any l_seaice=.false. run, have seaice_present=0 -> nothing to read).
        ! The ice fields are allocated lazily by aeros_ocean_init when l_seaice is
        ! on; guard on allocation so a mismatched pairing fails loudly rather than
        ! writing into an unallocated array. Rebuild the diagnostic albedo field
        ! from the restored state so the first restarted radiation step sees it.
        call nc_read_attr(filename,"seaice_present", seaice_present)
        if (seaice_present == 1) then
            if (.not. (allocated(ts%ocn%h_ice) .and. allocated(ts%ocn%a_ice) &
                       .and. allocated(ts%ocn%t_ice_sfc))) then
                write(io_unit_err,*) "aeros_timestep_read_restart:: error: restart "// &
                    "carries sea-ice state but the ocean was not initialised with "// &
                    "l_seaice = .true."
                error stop 1
            end if
            call nc_read(filename,"ocn_h_ice",     ts%ocn%h_ice)
            call nc_read(filename,"ocn_a_ice",     ts%ocn%a_ice)
            call nc_read(filename,"ocn_t_ice_sfc", ts%ocn%t_ice_sfc)
            call aeros_ocean_albedo_update(ts%ocn)
        end if

        if (rad_present == 1) then
            ! ts%rad%heat is allocated lazily on the first apply; a restart must
            ! fill it BEFORE the first restarted step so the cached heating is
            ! the one the continuous run would have applied.
            if (.not. allocated(ts%rad%heat)) &
                allocate(ts%rad%heat(nlon, nlat, nlev))
            call nc_read(filename,"rad_heat", ts%rad%heat)
        end if

        ! Cached radiative surface/TOA fluxes (always allocated by init).
        call nc_read(filename,"rad_sw_net_sur", ts%rad%sw_net_sur)
        call nc_read(filename,"rad_lw_dw_sur",  ts%rad%lw_dw_sur)
        call nc_read(filename,"rad_sw_dw_sur",  ts%rad%sw_dw_sur)
        call nc_read(filename,"rad_olr",        ts%rad%olr)
        call nc_read(filename,"rad_sw_up_toa",  ts%rad%sw_up_toa)

        ! === land state (feat/land) ==========================================
        ! Restore the two land prognostics. The caller must already have built
        ! ts%land the same way (same namelist/files) so the mask, albedo and the
        ! arrays are allocated; this only overwrites the prognostic values.
        call nc_read_attr(filename,"land_present", land_present)
        if (land_present == 1) then
            if (.not. allocated(ts%land%w)) then
                write(io_unit_err,*) "aeros_timestep_read_restart:: error: restart has land "// &
                                        "state but ts%land is not allocated (land disabled?)."
                error stop 1
            end if
            call nc_read(filename,"land_w",      ts%land%w)
            call nc_read(filename,"land_t_soil", ts%land%t_soil)
        end if
        ! === prognostic cloud state (feat/clouds) ===========================
        ! Restore now%cf_g when the file carries it (cf_present == 1). cf_g is
        ! always allocated by aeros_state_alloc, so a diagnostic-cloud restart
        ! simply leaves it at zero. Append-only block; see rad_present.
        call nc_read_attr(filename,"cf_present", cf_present)
        if (cf_present == 1) call nc_read(filename,"cf_g", now%cf_g)

        return

    contains

        subroutine read_spec(prefix, spec)
            ! Read one spectral level set from its real/imag component variables.
            character(len=*),       intent(in)    :: prefix
            type(aeros_spec_class), intent(inout) :: spec

            real(dp) :: re2(spec%nlm, spec%nlev), im2(spec%nlm, spec%nlev)
            real(dp) :: re1(spec%nlm), im1(spec%nlm)

            call nc_read(filename, prefix//"_vor_re",  re2)
            call nc_read(filename, prefix//"_vor_im",  im2)
            spec%vor  = cmplx(re2, im2, wp_sh)
            call nc_read(filename, prefix//"_div_re",  re2)
            call nc_read(filename, prefix//"_div_im",  im2)
            spec%div  = cmplx(re2, im2, wp_sh)
            call nc_read(filename, prefix//"_temp_re", re2)
            call nc_read(filename, prefix//"_temp_im", im2)
            spec%temp = cmplx(re2, im2, wp_sh)
            call nc_read(filename, prefix//"_lnps_re", re1)
            call nc_read(filename, prefix//"_lnps_im", im1)
            spec%lnps = cmplx(re1, im1, wp_sh)
            return
        end subroutine read_spec

    end subroutine aeros_timestep_read_restart

    subroutine sponge(ts, sht, h)
        ! Implicit model-top sponge (C1 Rayleigh drag + C2 Newtonian cooling),
        ! applied to the n+1 state in the top layers only. Same implicit form as
        ! diffuse: X <- (X + h k X_target)/(1 + h k). Vorticity and divergence
        ! relax toward zero; temperature toward sponge_tref in the mean and
        ! toward zero in its horizontal structure -- the latter is what removes
        ! the equator-pole gradient driving the runaway thermal wind.

        implicit none
        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_class),      intent(in)    :: sht
        real(wp), intent(in) :: h

        real(wp) :: fr, ft, tref_c
        integer  :: k, lm

        do k = 1, ts%new%nlev
            if (ts%sponge_kr_lev(k) <= 0.0_wp .and. ts%sponge_kt_lev(k) <= 0.0_wp) cycle

            fr = 1.0_wp/(1.0_wp + h*ts%sponge_kr_lev(k))
            ft = 1.0_wp/(1.0_wp + h*ts%sponge_kt_lev(k))
            tref_c = ts%sponge_tref*real(ts%c00, wp)   ! (0,0) coeff of a constant field

            do lm = 1, sht%nlm
                ts%new%vor(lm,k) = fr*ts%new%vor(lm,k)
                ts%new%div(lm,k) = fr*ts%new%div(lm,k)
                if (lm == ts%lm00) then
                    ts%new%temp(lm,k) = ft*(ts%new%temp(lm,k) &
                                            + cmplx(h*ts%sponge_kt_lev(k)*tref_c, 0.0_wp, wp_sh))
                else
                    ts%new%temp(lm,k) = ft*ts%new%temp(lm,k)
                end if
            end do
        end do

        return
    end subroutine sponge

    subroutine diffuse(ts, sht, h)
        ! Implicit del^ndiff damping of vorticity, divergence and temperature.
        !
        ! The per-level factor is fac(l,k) = 1/(1 + mult_k (h/tau) dratio_lev(l,k)),
        ! where dratio_lev carries the (optionally sigma-tapered) order and mult_k
        ! the (optional) vorticity scaling. With both adaptive mechanisms off,
        ! mult_k = 1 and dratio_lev(:,k) = dratio(:), so fac reduces exactly to the
        ! fixed 1-D scheme -- bit for bit. See the module header.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_class),      intent(in)    :: sht
        real(wp), intent(in) :: h

        real(wp) :: fac(0:sht%lmax, ts%new%nlev), rate, mult
        integer  :: l, k, lm

        rate = h/ts%tau_diff
        do k = 1, ts%new%nlev
            mult = 1.0_wp
            if (ts%diff_adapt) mult = diff_adapt_mult(ts, sht, k)
            do l = 0, sht%lmax
                fac(l,k) = 1.0_wp/(1.0_wp + mult*rate*ts%dratio_lev(l,k))
            end do
        end do

        !$omp parallel do collapse(2) schedule(static) private(k,lm)
        do k = 1, ts%new%nlev
            do lm = 1, sht%nlm
                ts%new%vor(lm,k)  = fac(sht%l_of_lm(lm),k)*ts%new%vor(lm,k)
                ts%new%div(lm,k)  = fac(sht%l_of_lm(lm),k)*ts%new%div(lm,k)
                ts%new%temp(lm,k) = fac(sht%l_of_lm(lm),k)*ts%new%temp(lm,k)
            end do
        end do
        !$omp end parallel do

        return

    end subroutine diffuse

    real(wp) function diff_adapt_mult(ts, sht, k) result(mult)
        ! The vorticity-scaled diffusion multiplier for level k (module header):
        !   mult = min( 1 + gain*max(0, |zeta|_k/zeta_ref - 1), adapt_max ).
        ! |zeta|_k is the level's RMS vorticity [s-1] from Parseval on the
        ! spectral coefficients -- surface-mean-square = sum_lm w_m |vor|^2/(4 pi),
        ! w_m = 1 (m=0) or 2 (m>0), 4 pi = c00^2 -- so no grid synthesis is
        ! needed. l = 0 (which vorticity never populates) is skipped.

        implicit none

        type(aeros_timestep_class), intent(in) :: ts
        type(aeros_sht_class),      intent(in) :: sht
        integer,                    intent(in) :: k

        real(dp) :: s2, w
        real(wp) :: zrms, excess
        integer  :: lm

        s2 = 0.0_dp
        do lm = 1, sht%nlm
            if (sht%l_of_lm(lm) == 0) cycle
            w = merge(1.0_dp, 2.0_dp, sht%m_of_lm(lm) == 0)
            s2 = s2 + w*real(ts%new%vor(lm,k)*conjg(ts%new%vor(lm,k)), dp)
        end do
        zrms = real(sqrt(s2)/ts%c00, wp)            ! c00^2 = 4 pi

        excess = zrms/ts%diff_zeta_ref - 1.0_wp
        mult   = 1.0_wp + ts%diff_adapt_gain*max(0.0_wp, excess)
        mult   = min(mult, ts%diff_adapt_max)

        return

    end function diff_adapt_mult

    subroutine raw_filter(ts, sht, cur, nlev)
        ! Williams (2009) Robert-Asselin-Williams filter, in place on X^n and
        ! X^(n+1). With raw_alpha = 1 this is the classical Robert-Asselin
        ! filter, acting on X^n alone.
        !
        ! Runs BEFORE the shift, so X^n is still the caller's state and is
        ! passed in as `cur`.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_class),      intent(in)    :: sht
        type(aeros_spec_class),     intent(inout) :: cur     ! X^n
        integer, intent(in) :: nlev

        complex(wp_sh) :: d
        real(wp) :: half_nu, alpha, beta
        integer  :: k, lm

        half_nu = 0.5_wp*ts%eps_filter
        alpha   = ts%raw_alpha
        beta    = 1.0_wp - ts%raw_alpha

        if (ts%eps_filter <= 0.0_wp) return

        !$omp parallel do collapse(2) schedule(static) private(k,lm,d)
        do k = 1, nlev
            do lm = 1, sht%nlm

                d = half_nu*(ts%old%vor(lm,k) - 2.0_wp*cur%vor(lm,k) + ts%new%vor(lm,k))
                cur%vor(lm,k)    = cur%vor(lm,k)    + alpha*d
                ts%new%vor(lm,k) = ts%new%vor(lm,k) - beta*d

                d = half_nu*(ts%old%div(lm,k) - 2.0_wp*cur%div(lm,k) + ts%new%div(lm,k))
                cur%div(lm,k)    = cur%div(lm,k)    + alpha*d
                ts%new%div(lm,k) = ts%new%div(lm,k) - beta*d

                d = half_nu*(ts%old%temp(lm,k) - 2.0_wp*cur%temp(lm,k) + ts%new%temp(lm,k))
                cur%temp(lm,k)    = cur%temp(lm,k)    + alpha*d
                ts%new%temp(lm,k) = ts%new%temp(lm,k) - beta*d

            end do
        end do
        !$omp end parallel do

        do lm = 1, sht%nlm
            d = half_nu*(ts%old%lnps(lm) - 2.0_wp*cur%lnps(lm) + ts%new%lnps(lm))
            cur%lnps(lm)    = cur%lnps(lm)    + alpha*d
            ts%new%lnps(lm) = ts%new%lnps(lm) - beta*d
        end do

        return

    end subroutine raw_filter

    subroutine aeros_timestep_set_phis(ts, phis)
        ! Set the surface geopotential [m2 s-2].
        !
        ! A prescribed boundary field, zero until topography arrives at M2 (and
        ! zero for Held-Suarez, which is defined over a flat surface). It lives
        ! in the tendency's work arrays because that is where it is consumed;
        ! this exists so a caller sets it through the interface rather than
        ! reaching into scratch, and so there is one identified place to change
        ! when topography stops being a boundary condition and starts being
        ! part of the model configuration.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        real(wp), intent(in) :: phis(:,:)

        if (size(phis,1) /= ts%wrk%nlon .or. size(phis,2) /= ts%wrk%nlat) then
            write(io_unit_err,*) "aeros_timestep_set_phis:: error: shape mismatch, got ", &
                                    size(phis,1), size(phis,2), " expected ", &
                                    ts%wrk%nlon, ts%wrk%nlat
            error stop 1
        end if

        ts%wrk%phis = phis

        return

    end subroutine aeros_timestep_set_phis

    subroutine aeros_timestep_set_qforce(ts, q)
        ! Set the diagnostic prescribed heating rate [K/s], (nlon,nlat,nlev).
        ! Added to wrk%dtdt at step 1b (the in-solve seam), independent of the
        ! physics toggles -- a controlled analytic Q for isolating the dry
        ! dynamical core's heating->ω response. Not a production forcing.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        real(wp), intent(in) :: q(:,:,:)

        if (size(q,1) /= ts%wrk%nlon .or. size(q,2) /= ts%wrk%nlat &
            .or. size(q,3) /= ts%wrk%nlev) then
            write(io_unit_err,*) "aeros_timestep_set_qforce:: error: shape mismatch, got ", &
                                    shape(q), " expected ", &
                                    ts%wrk%nlon, ts%wrk%nlat, ts%wrk%nlev
            error stop 1
        end if

        if (.not. allocated(ts%q_force)) &
            allocate(ts%q_force(ts%wrk%nlon, ts%wrk%nlat, ts%wrk%nlev))
        ts%q_force = q

        return

    end subroutine aeros_timestep_set_qforce

    subroutine aeros_timestep_diagnose(ts, pool, vg, grd, now)
        ! Synthesize the grid-space fields of `now` from its spectral state.
        !
        ! NOT done every step. The integrator itself never needs u, v, T or p_s
        ! on the grid outside the tendency evaluation, which builds them in its
        ! own scratch; only output, the conservation budgets and (from M2) the
        ! column physics do. Calling this on the output interval rather than on
        ! the dynamics step is ~4 transforms per level saved on 99 steps out of
        ! 100, and it keeps now%u honestly consistent with now%spec whenever it
        ! has been called, rather than lagging it by a step.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_vgrid_class),    intent(in)    :: vg
        type(aeros_grid_class),     intent(in)    :: grd
        type(aeros_state_class),    intent(inout) :: now

        type(aeros_sht_class), pointer :: s
        integer :: i, j, k

        s => pool%sht(1)

        call aeros_sht_synthesis(s, now%spec%lnps, ts%wrk%lnps_g)
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                now%ps(i,j) = exp(ts%wrk%lnps_g(i,j))
            end do
        end do

        !$omp parallel do num_threads(pool%nthreads) private(k,s) schedule(static)
        do k = 1, vg%nlev
            s => aeros_sht_pool_get(pool)
            call aeros_uv_from_vordiv(s, now%spec%vor(:,k), now%spec%div(:,k), &
                                        now%u(:,:,k), now%v(:,:,k))
            call aeros_sht_synthesis(s, now%spec%temp(:,k), now%temp_g(:,:,k))
        end do
        !$omp end parallel do

        ! Humidity is NOT synthesized here. now%qv_g is a prognostic in its own
        ! right, advanced on the grid by aeros_transport, and persists across
        ! steps -- there is no spectral qv to rebuild it from.

        return

    end subroutine aeros_timestep_diagnose

    subroutine aeros_timestep_print(ts, io_unit)
        ! Report the integrator's configuration.

        implicit none

        type(aeros_timestep_class), intent(in) :: ts
        integer, intent(in), optional :: io_unit

        integer :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        write(iou,*) ""
        write(iou,"(a)")        " == time integration =="
        write(iou,"(a,f9.1,a)") "   dt                         ", ts%dt, " s"
        if (ts%semi_implicit) then
            write(iou,"(a)")    "   scheme                      semi-implicit leapfrog"
        else
            write(iou,"(a)")    "   scheme                      EXPLICIT leapfrog (not production)"
        end if
        write(iou,"(a,f9.3)")   "   filter nu                  ", ts%eps_filter
        write(iou,"(a,f9.3,a)") "   filter alpha               ", ts%raw_alpha, &
                                    "  (1 = classical Robert-Asselin)"
        write(iou,"(a,i9)")     "   diffusion order            ", ts%ndiff
        write(iou,"(a,f9.2,a)") "   diffusion e-folding at lmax", ts%tau_diff/3600.0_wp, " h"
        if (ts%diff_adapt) then
            write(iou,"(a,es9.2,a,f6.2,a,f6.2)") &
                "   diff vorticity-scaling   ON  zeta_ref", ts%diff_zeta_ref, &
                " s-1  gain", ts%diff_adapt_gain, "  max", ts%diff_adapt_max
        else
            write(iou,"(a)")    "   diff vorticity-scaling      off"
        end if
        if (ts%diff_taper) then
            write(iou,"(a,i2,a,i2,a,f6.3)") &
                "   diff order taper         ON  del^", ts%ndiff, " -> del^", &
                ts%diff_ndiff_top, " above sigma", ts%diff_taper_sigma
        else
            write(iou,"(a)")    "   diff order taper            off"
        end if
        if (ts%mass_fixer) then
            write(iou,"(a)")    "   mass fixer                  ON  (p_s rescaled each step)"
        else
            write(iou,"(a)")    "   mass fixer                  off"
        end if

        call aeros_semiimp_print(ts%si, iou)
        call aeros_hs_print(ts%hs, iou)
        call aeros_correction_report(ts%cor, iou)
        call aeros_transport_report(ts%mst, iou)
        call aeros_convection_report(ts%cnv, iou)
        call aeros_condensation_report(ts%cnd, iou)
        call aeros_surface_report(ts%surf, iou)
        call aeros_ocean_report(ts%ocn, iou)
        call aeros_radiation_report(ts%rad, iou)
        call aeros_vdiff_report(ts%vd, iou)

        return

    end subroutine aeros_timestep_print

end module aeros_timestep
