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
    use aeros_vordiv,   only : aeros_uv_from_vordiv
    use aeros_tendency, only : aeros_tend_class, aeros_work_class, &
                                aeros_tend_alloc, aeros_tend_end, &
                                aeros_work_alloc, aeros_work_end, &
                                aeros_tendency_grid, aeros_tendency_spectral
    use aeros_semiimp,  only : aeros_semiimp_class, aeros_semiimp_init, &
                                aeros_semiimp_end, aeros_semiimp_set_step, &
                                aeros_semiimp_step, aeros_semiimp_print
    use aeros_correction, only : aeros_correction_class, aeros_correction_init, &
                                aeros_correction_load, aeros_correction_end, &
                                aeros_correction_apply, aeros_correction_report
    use aeros_moisture, only : aeros_moist_class, aeros_moisture_init, &
                                aeros_moisture_end, aeros_moisture_transport, &
                                aeros_moisture_report
    use aeros_condensation, only : aeros_cond_class, aeros_condensation_init, &
                                aeros_condensation_load, aeros_condensation_end, &
                                aeros_condensation_apply, aeros_condensation_report
    use aeros_held_suarez, only : aeros_hs_class, aeros_hs_init, aeros_hs_end, &
                                aeros_hs_apply, aeros_hs_print

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

        ! The additive tendency correction (design.md 3.7, 3.8). A DIFFERENT
        ! seam from `hs` above -- spectral, not grid -- and the reason is in
        ! aeros_correction's header. Disabled unless a namelist says otherwise,
        ! so every test and every driver that does not ask for it gets the
        ! uncorrected model bit for bit.
        type(aeros_correction_class) :: cor

        ! Prognostic humidity transport (aeros_moisture). q is a gridpoint
        ! field carried in now%qv_g and advected here, off the spectral core.
        ! Always allocated; it transports whatever is in qv_g, which for a dry
        ! run is zero and stays zero.
        type(aeros_moist_class) :: mst

        ! Large-scale condensation (aeros_condensation), at the grid seam. Off
        ! unless a namelist turns it on; a dry run never enters it.
        type(aeros_cond_class) :: cnd

        ! [ l(l+1) / lmax(lmax+1) ]^(ndiff/2), precomputed per degree. The
        ! step-dependent part is one multiply, so this never needs rebuilding.
        real(wp), allocatable :: dratio(:)   ! (0:lmax)

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
    public :: aeros_timestep_set_mass_target
    public :: aeros_timestep_diagnose
    public :: aeros_timestep_print

contains

    subroutine aeros_timestep_init(ts, par, pool, grd, vg, filename)
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

        integer  :: nlm, lmax, l
        real(wp) :: lref
        complex(wp_sh), allocatable :: c_one(:)

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
            call aeros_correction_load(ts%cor, filename, pool%sht(1), vg%nlev)
        else
            call aeros_correction_init(ts%cor, pool%sht(1), vg%nlev)
        end if

        call aeros_moisture_init(ts%mst, grd, vg%nlev)

        if (present(filename)) then
            call aeros_condensation_load(ts%cnd, filename, grd)
        else
            call aeros_condensation_init(ts%cnd, grd, .FALSE.)
        end if

        allocate(ts%dratio(0:lmax))
        lref = real(lmax*(lmax+1), wp)
        do l = 0, lmax
            ts%dratio(l) = (real(l*(l+1), wp)/lref)**(ts%ndiff/2)
        end do

        if (ts%mass_fixer) then
            allocate(ts%ps_fix(grd%nlon, grd%nlat))

            ! Read the (0,0) coefficient of the constant field 1 out of a
            ! transform. Analysis overwrites its input, which is why this uses
            ! the scratch it has just allocated rather than a saved field.
            ts%lm00 = aeros_sht_lm(pool%sht(1), 0, 0)

            allocate(c_one(nlm))
            ts%ps_fix = 1.0_dp
            call aeros_sht_analysis(pool%sht(1), ts%ps_fix, c_one)
            ts%c00 = real(c_one(ts%lm00), dp)
            deallocate(c_one)

            ! The constant field 1 must have a non-zero (0,0) coefficient in
            ! any normalization. A zero means the transform or the indexing is
            ! not what the fixer believes, and it is better to say so at init
            ! than to divide by it a thousand steps in.
            if (abs(ts%c00) < 1.0e-12_dp) then
                write(io_unit_err,*) "aeros_timestep_init:: error: the mass fixer read a "// &
                                        "(0,0) coefficient of ", ts%c00, " for the constant "// &
                                        "field 1, which cannot be right."
                error stop 1
            end if
        end if

        return

    end subroutine aeros_timestep_init

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
        call aeros_moisture_end(ts%mst)
        call aeros_condensation_end(ts%cnd)

        if (allocated(ts%dratio)) deallocate(ts%dratio)
        if (allocated(ts%ps_fix)) deallocate(ts%ps_fix)

        ts%mass_target = -1.0_dp
        ts%lnr_last    = 0.0_dp
        ts%lnr_cum     = 0.0_dp

        return

    end subroutine aeros_timestep_end

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

        ! Large-scale condensation, at the same grid seam and for the same
        ! reason: it dries the gridpoint humidity in place and adds its latent
        ! heating to wrk%dtdt, so the heating rides the transform and the
        ! leapfrog with the dynamical temperature tendency. Returns at once when
        ! dry. See aeros_condensation.
        call aeros_condensation_apply(ts%cnd, vg, ts%wrk%t_g, now%qv_g, &
                                        ts%wrk%lnps_g, ts%wrk%dtdt, ts%dt)

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

        ! === 4. Time filter ==================================================
        if (ts%nstep > 0) call raw_filter(ts, s, now%spec, vg%nlev)

        ! === 5. Shift the time levels ========================================
        ! old <- now <- new, by exchanging descriptors: the retired X^(n-1)
        ! becomes the scratch that the next step writes X^(n+1) into.
        call aeros_spec_swap(ts%old, now%spec)
        call aeros_spec_swap(now%spec, ts%new)

        ! === 6. Mass fixer ===================================================
        ! After the swap, so it acts on X^(n+1) -- the state the caller is
        ! about to receive -- and after the filter, so nothing follows it that
        ! could reintroduce what it just removed.
        if (ts%mass_fixer) call mass_fix(ts, s, now%spec)

        ! === 7. Humidity transport ===========================================
        ! On the grid, forward in time, using the winds and surface pressure at
        ! the time level just stepped from -- which are still in ts%wrk, filled
        ! by aeros_tendency_grid and untouched since. q (now%qv_g) persists
        ! across steps; nothing here is spectral. See aeros_moisture.
        call aeros_moisture_transport(ts%mst, vg, ts%wrk%u, ts%wrk%v, &
                                        ts%wrk%lnps_g, now%qv_g, ts%dt)

        ts%nstep = ts%nstep + 1

        return

    end subroutine aeros_timestep_step

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

    subroutine diffuse(ts, sht, h)
        ! Implicit del^ndiff damping of vorticity, divergence and temperature.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_sht_class),      intent(in)    :: sht
        real(wp), intent(in) :: h

        real(wp) :: fac(0:sht%lmax), rate
        integer  :: l, k, lm

        rate = h/ts%tau_diff
        do l = 0, sht%lmax
            fac(l) = 1.0_wp/(1.0_wp + rate*ts%dratio(l))
        end do

        !$omp parallel do collapse(2) schedule(static) private(k,lm)
        do k = 1, ts%new%nlev
            do lm = 1, sht%nlm
                ts%new%vor(lm,k)  = fac(sht%l_of_lm(lm))*ts%new%vor(lm,k)
                ts%new%div(lm,k)  = fac(sht%l_of_lm(lm))*ts%new%div(lm,k)
                ts%new%temp(lm,k) = fac(sht%l_of_lm(lm))*ts%new%temp(lm,k)
            end do
        end do
        !$omp end parallel do

        return

    end subroutine diffuse

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
        ! right, advanced on the grid by aeros_moisture, and persists across
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
        if (ts%mass_fixer) then
            write(iou,"(a)")    "   mass fixer                  ON  (p_s rescaled each step)"
        else
            write(iou,"(a)")    "   mass fixer                  off"
        end if

        call aeros_semiimp_print(ts%si, iou)
        call aeros_hs_print(ts%hs, iou)
        call aeros_correction_report(ts%cor, iou)
        call aeros_moisture_report(ts%mst, iou)
        call aeros_condensation_report(ts%cnd, iou)

        return

    end subroutine aeros_timestep_print

end module aeros_timestep
