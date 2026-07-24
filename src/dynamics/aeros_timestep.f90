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
    ! SIGMA-SURFACE DIFFUSION: diffusing T along model levels rather than
    ! pressure surfaces drives spurious heat fluxes over steep orography. It
    ! does not arise at M1 -- the surface is flat -- and the correction belongs
    ! with topography at M2, where docs/design.md section 4.2 already flags the
    ! ice-margin Gibbs problem it is entangled with.

    use aeros_defs,     only : dp, wp, wp_sh, io_unit_err, r_earth, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_spec_class, aeros_state_class
    use aeros_spectral, only : aeros_sht_class, aeros_sht_pool_class, &
                                aeros_sht_pool_get, aeros_sht_synthesis
    use aeros_state,    only : aeros_spec_alloc, aeros_spec_end, aeros_spec_copy, &
                                aeros_spec_swap
    use aeros_vertical, only : aeros_vgrid_class
    use aeros_vordiv,   only : aeros_uv_from_vordiv
    use aeros_tendency, only : aeros_tend_class, aeros_work_class, &
                                aeros_tend_alloc, aeros_tend_end, &
                                aeros_work_alloc, aeros_work_end, &
                                aeros_tendency_calc
    use aeros_semiimp,  only : aeros_semiimp_class, aeros_semiimp_init, &
                                aeros_semiimp_end, aeros_semiimp_set_step, &
                                aeros_semiimp_step, aeros_semiimp_print

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

        ! [ l(l+1) / lmax(lmax+1) ]^(ndiff/2), precomputed per degree. The
        ! step-dependent part is one multiply, so this never needs rebuilding.
        real(wp), allocatable :: dratio(:)   ! (0:lmax)
    end type aeros_timestep_class

    public :: aeros_timestep_class
    public :: aeros_timestep_init
    public :: aeros_timestep_end
    public :: aeros_timestep_step
    public :: aeros_timestep_set_phis
    public :: aeros_timestep_diagnose
    public :: aeros_timestep_print

contains

    subroutine aeros_timestep_init(ts, par, pool, grd, vg)
        ! Allocate the time levels and the scratch, and build the solve.

        implicit none

        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_param_class),    intent(in)    :: par
        type(aeros_sht_pool_class), intent(in)    :: pool
        type(aeros_grid_class),     intent(in)    :: grd
        type(aeros_vgrid_class),    intent(in)    :: vg

        integer  :: nlm, lmax, l
        real(wp) :: lref

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
        ts%nstep         = 0

        call aeros_spec_alloc(ts%old, nlm, vg%nlev)
        call aeros_spec_alloc(ts%new, nlm, vg%nlev)
        call aeros_tend_alloc(ts%tnd, nlm, vg%nlev)
        call aeros_work_alloc(ts%wrk, grd%nlon, grd%nlat, vg%nlev)

        ! Factorized for the start-up step; switched to 2 dt after it.
        call aeros_semiimp_init(ts%si, vg, lmax, ts%dt)

        allocate(ts%dratio(0:lmax))
        lref = real(lmax*(lmax+1), wp)
        do l = 0, lmax
            ts%dratio(l) = (real(l*(l+1), wp)/lref)**(ts%ndiff/2)
        end do

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

        if (allocated(ts%dratio)) deallocate(ts%dratio)

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
        else
            h = 2.0_wp*ts%dt
        end if

        if (ts%si%dt_step /= h) call aeros_semiimp_set_step(ts%si, h)

        ! === 1. The nonlinear right-hand sides ===============================
        call aeros_tendency_calc(pool, vg, grd, now%spec, ts%wrk, ts%tnd)

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

        ! Humidity is not prognostic until M2; carry it so `new` never holds
        ! stale data from two steps ago.
        ts%new%qv = ts%old%qv

        ! === 3. Horizontal diffusion, implicit ===============================
        call diffuse(ts, s, h)

        ! === 4. Time filter ==================================================
        if (ts%nstep > 0) call raw_filter(ts, s, now%spec, vg%nlev)

        ! === 5. Shift the time levels ========================================
        ! old <- now <- new, by exchanging descriptors: the retired X^(n-1)
        ! becomes the scratch that the next step writes X^(n+1) into.
        call aeros_spec_swap(ts%old, now%spec)
        call aeros_spec_swap(now%spec, ts%new)

        ts%nstep = ts%nstep + 1

        return

    end subroutine aeros_timestep_step

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
            call aeros_sht_synthesis(s, now%spec%qv(:,k),   now%qv_g(:,:,k))
        end do
        !$omp end parallel do

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

        call aeros_semiimp_print(ts%si, iou)

        return

    end subroutine aeros_timestep_print

end module aeros_timestep
