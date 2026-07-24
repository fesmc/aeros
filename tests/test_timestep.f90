program test_timestep
    ! Acceptance test for the time integrator.
    !
    ! Everything before M1.4 was a single tendency evaluation. These are the
    ! first tests in aeros that INTEGRATE, so they are the first that can say
    ! anything about stability or about what the model conserves.
    !
    ! Four things, in order of how much they would hurt to get wrong:
    !
    ! 1. A BALANCED RESTING ATMOSPHERE STAYS AT REST. aeros_tendency already
    !    showed the tendency of that state vanishes; this shows the integrator
    !    wrapped around it does not manufacture motion out of the start-up
    !    step, the filter, the diffusion or the semi-implicit correction. It is
    !    run over topography, on the shipped hybrid coordinate, for long enough
    !    that a slow instability would show.
    !
    ! 2. THE SEMI-IMPLICIT SOLVE EARNS ITS KEEP. The same moving state is
    !    integrated three ways at the same 30-minute step: semi-implicit
    !    (stable), explicit (must diverge), and explicit at a step below the
    !    gravity-wave limit (stable again). The third run is what makes the
    !    second one evidence rather than a broken code path.
    !
    ! 3. THE DIFFUSION DAMPS AT THE RATE IT ADVERTISES. `tau_diff` claims to be
    !    the e-folding time of the l = lmax wave, so a small single-mode
    !    perturbation there must decay by the analytic implicit factor.
    !
    ! 4. CONSERVATION, measured on two states because they answer different
    !    questions. On the BALANCED RESTING atmosphere all three integrals hold
    !    to machine precision, and that is asserted: nothing there can drift
    !    except a genuine leak in the integrator. On a MOVING state mass is a
    !    bounded oscillation of order 1e-5 -- the discrete-time signature of a
    !    gravity-wave adjustment, since what the scheme conserves exactly is
    !    d/dt int p_s dA rather than int p_s dA itself -- and energy and
    !    angular momentum are REPORTED rather than asserted, because leapfrog,
    !    the time filter and a deliberate hyperdiffusion sink all remove
    !    energy. What a paleo integration needs from those two is a drift rate,
    !    not a claim.
    !
    ! 5. THE MASS FIXER DOES WHAT IT CLAIMS AND NOTHING ELSE. It restores the
    !    mass exactly, it moves only the (0,0) surface-pressure coefficient,
    !    and each step's reported correction is exactly that step's leak. The
    !    third is the one that matters most: a fixer nobody can audit is a way
    !    of not knowing what the model is doing.
    !
    !    The CUMULATIVE correction is reported here too, and deliberately not
    !    asserted against the unfixed drift, because it is a larger number and
    !    that is not a bug -- see the comment in test_mass_fixer.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, R_d, p0, grav, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_spec_class, aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_budget
    use aeros_vertical
    use aeros_tendency
    use aeros_semiimp
    use aeros_timestep

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 20

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg

    integer :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)                 ! the shipped hybrid default

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_timestep:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    call test_rest(nfail)
    call test_stability(nfail)
    call test_diffusion(nfail)
    call test_conservation(nfail)
    call test_mass_fixer(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_timestep:: PASS"
    else
        write(*,*) "test_timestep:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine default_par(par, dt, semi_implicit, eps_filter, tau_diff, mass_fixer)
        ! The shipped par/aeros.nml values, with the ones a test varies exposed.

        implicit none

        type(aeros_param_class), intent(out) :: par
        real(wp), intent(in) :: dt
        logical,  intent(in) :: semi_implicit
        real(wp), intent(in), optional :: eps_filter, tau_diff
        logical,  intent(in), optional :: mass_fixer

        par%trunc    = trunc
        par%nlon     = -1
        par%nlat     = -1
        par%nlev     = nlev
        par%nthreads = -1

        par%dt            = dt
        par%semi_implicit = semi_implicit
        par%held_suarez   = .FALSE.
        par%eps_filter    = 0.06_wp
        par%raw_alpha     = 0.53_wp
        par%ndiff         = 6
        par%tau_diff      = 6.0_wp
        par%mass_fixer    = .FALSE.

        if (present(eps_filter)) par%eps_filter = eps_filter
        if (present(tau_diff))   par%tau_diff   = tau_diff
        if (present(mass_fixer)) par%mass_fixer = mass_fixer

        return

    end subroutine default_par

    subroutine rest_state(now, phis, topography)
        ! An isothermal atmosphere at rest, in exact discrete hydrostatic
        ! balance with its own surface geopotential -- the same construction as
        ! tests/test_tendency.f90, and for the same reason: ln p_s is defined
        ! SPECTRALLY and Phi_s computed from the synthesized grid field, so the
        ! test measures the discretization and not a truncation residual
        ! between two independently written analytic fields.

        implicit none

        type(aeros_state_class), intent(inout) :: now
        real(wp), intent(out) :: phis(:,:)
        logical,  intent(in)  :: topography

        type(aeros_sht_class), pointer :: s
        real(wp), parameter :: tiso = 280.0_wp
        real(wp) :: lnps_g(grd%nlon,grd%nlat)
        real(dp) :: y00
        integer  :: k, i, j, lm

        s => pool%sht(1)
        y00 = sqrt(4.0_dp*acos(-1.0_dp))

        call aeros_spec_zero(now%spec)

        do k = 1, nlev
            now%spec%temp(aeros_sht_lm(s,0,0),k) = cmplx(real(tiso,dp)*y00, 0.0_dp, wp_sh)
        end do
        now%spec%lnps(aeros_sht_lm(s,0,0)) = cmplx(log(real(p0,dp))*y00, 0.0_dp, wp_sh)

        if (topography) then
            ! Surface elevations of order a kilometre, so grad(ln p_s) is
            ! genuinely large and every term in the balance is active.
            lm = aeros_sht_lm(s,2,1); now%spec%lnps(lm) = cmplx(-0.06_dp,  0.03_dp, wp_sh)
            lm = aeros_sht_lm(s,3,0); now%spec%lnps(lm) = cmplx( 0.05_dp,  0.0_dp,  wp_sh)
            lm = aeros_sht_lm(s,5,3); now%spec%lnps(lm) = cmplx( 0.02_dp, -0.01_dp, wp_sh)
        end if

        call aeros_sht_synthesis(s, now%spec%lnps, lnps_g)

        do j = 1, grd%nlat
            do i = 1, grd%nlon
                phis(i,j) = R_d*tiso*(log(p0) - lnps_g(i,j))
            end do
        end do

        return

    end subroutine rest_state

    subroutine moving_state(now)
        ! A sheared, rotational flow with realistic magnitudes, over a flat
        ! surface. Same construction as test_tendency's moving state.

        implicit none

        type(aeros_state_class), intent(inout) :: now

        type(aeros_sht_class), pointer :: s
        real(wp), parameter :: tiso = 280.0_wp
        real(dp) :: amp, y00
        integer  :: k, lm, l, m

        s => pool%sht(1)
        y00 = sqrt(4.0_dp*acos(-1.0_dp))

        call aeros_spec_zero(now%spec)

        now%spec%lnps(aeros_sht_lm(s,0,0)) = cmplx(log(real(p0,dp))*y00, 0.0_dp, wp_sh)
        now%spec%lnps(aeros_sht_lm(s,2,1)) = cmplx(-0.03_dp, 0.01_dp, wp_sh)

        do k = 1, nlev
            now%spec%temp(aeros_sht_lm(s,0,0),k) = cmplx(real(tiso,dp)*y00, 0.0_dp, wp_sh)
            now%spec%temp(aeros_sht_lm(s,2,0),k) = &
                    cmplx(-20.0_dp*real(k,dp)/real(nlev,dp), 0.0_dp, wp_sh)
            do lm = 1, s%nlm
                l = s%l_of_lm(lm); m = s%m_of_lm(lm)
                if (l < 1 .or. l > 6) cycle
                amp = 1.0e-5_dp/real(l*l, dp)*(0.5_dp + 0.5_dp*real(k,dp)/real(nlev,dp))
                now%spec%vor(lm,k) = cmplx(amp*cos(real(3*l+m,dp)), amp*sin(real(l+2*m,dp)), wp_sh)
                now%spec%div(lm,k) = cmplx(0.2_dp*amp*sin(real(l+m,dp)), &
                                            0.2_dp*amp*cos(real(2*l+m,dp)), wp_sh)
                if (m == 0) then
                    now%spec%vor(lm,k) = cmplx(real(now%spec%vor(lm,k)), 0.0_wp_sh, wp_sh)
                    now%spec%div(lm,k) = cmplx(real(now%spec%div(lm,k)), 0.0_wp_sh, wp_sh)
                end if
            end do
        end do

        return

    end subroutine moving_state

    ! === 1. The resting atmosphere ===========================================

    subroutine test_rest(nfail)

        implicit none

        integer, intent(inout) :: nfail

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(wp) :: ps0(grd%nlon,grd%nlat)
        real(dp) :: umax, dps, tmax
        integer  :: n
        integer, parameter :: nstep = 100

        write(*,*) ""
        write(*,*) " -- a balanced resting atmosphere, over topography, 100 steps"

        call default_par(par, dt=1800.0_wp, semi_implicit=.TRUE.)
        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        call rest_state(now, phis, topography=.TRUE.)
        call aeros_timestep_set_phis(ts, phis)

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        ps0 = now%ps

        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)
        end do

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)

        umax = max(maxval(abs(real(now%u,dp))), maxval(abs(real(now%v,dp))))
        dps  = maxval(abs(real(now%ps,dp) - real(ps0,dp)))
        tmax = maxval(abs(real(now%temp_g,dp) - 280.0_dp))

        write(*,"(a40,es12.3,a)") "   max |wind| after 100 steps ", umax, " m s-1"
        write(*,"(a40,es12.3,a)") "   max |dp_s|                 ", dps,  " Pa"
        write(*,"(a40,es12.3,a)") "   max |T - 280|              ", tmax, " K"

        ! 1 mm/s of spurious wind after two days is ~5 orders of magnitude
        ! below anything the model would call weather.
        call check(umax < 1.0e-3_dp, "the resting state generates no wind", nfail)
        call check(dps  < 1.0e-2_dp, "surface pressure does not drift", nfail)
        call check(tmax < 1.0e-4_dp, "temperature does not drift", nfail)

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine test_rest

    ! === 2. Stability ========================================================

    subroutine test_stability(nfail)

        implicit none

        integer, intent(inout) :: nfail

        real(dp) :: g_si, g_ex, g_ex_short
        logical  :: d_si, d_ex, d_ex_short

        write(*,*) ""
        write(*,*) " -- stability at and below the gravity-wave limit"
        write(*,*) "    (peak growth of max|zeta| over 200 steps, and whether it diverged)"

        call integrate_growth(1800.0_wp, .TRUE.,  200, g_si,       d_si)
        call integrate_growth(1800.0_wp, .FALSE., 200, g_ex,       d_ex)
        call integrate_growth( 300.0_wp, .FALSE., 200, g_ex_short, d_ex_short)

        write(*,"(a40,es12.3,l4)") "   semi-implicit, dt = 1800 s ", g_si,       d_si
        write(*,"(a40,es12.3,l4)") "   explicit,      dt = 1800 s ", g_ex,       d_ex
        write(*,"(a40,es12.3,l4)") "   explicit,      dt =  300 s ", g_ex_short, d_ex_short

        call check(.not. d_si .and. g_si < 10.0_dp, &
                    "semi-implicit is stable at dt = 1800 s", nfail)
        call check(d_ex, &
                    "explicit diverges at dt = 1800 s", nfail)
        call check(.not. d_ex_short .and. g_ex_short < 10.0_dp, &
                    "explicit is stable below the gravity-wave limit", nfail)

        return

    end subroutine test_stability

    subroutine integrate_growth(dt, semi_implicit, nstep, growth, diverged)
        ! Integrate the moving state and report the PEAK growth of the
        ! vorticity, and whether the run went non-finite.
        !
        ! Peak rather than final, and a separate `diverged` flag rather than a
        ! threshold on the return value, because once a field is all-NaN
        ! `maxval` reports 0 -- so a diverged run looks quieter than a healthy
        ! one if the last value is all that is kept. Stops as soon as either
        ! happens, so an overflowing run does not spend its remaining steps
        ! producing NaNs.

        implicit none

        real(wp), intent(in)  :: dt
        logical,  intent(in)  :: semi_implicit
        integer,  intent(in)  :: nstep
        real(dp), intent(out) :: growth
        logical,  intent(out) :: diverged

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(dp) :: z0, z
        integer  :: n

        call default_par(par, dt=dt, semi_implicit=semi_implicit)
        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        call moving_state(now)
        phis = 0.0_wp
        call aeros_timestep_set_phis(ts, phis)

        z0       = maxval(abs(now%spec%vor))
        growth   = 1.0_dp
        diverged = .FALSE.

        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)

            if (.not. all(now%spec%vor == now%spec%vor)) then
                diverged = .TRUE.
                exit
            end if

            z      = maxval(abs(now%spec%vor))
            growth = max(growth, z/z0)

            if (growth > 1.0e3_dp) then
                diverged = .TRUE.
                exit
            end if
        end do

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine integrate_growth

    ! === 3. Diffusion ========================================================

    subroutine test_diffusion(nfail)
        ! The l = lmax mode must decay at exactly the implicit factor the
        ! namelist asks for.
        !
        ! The perturbation is deliberately TINY (1e-12 s-1 of vorticity, ~1e-6
        ! m s-1 of wind) so the nonlinear terms, which are quadratic in it, are
        ! ~1e-6 of the linear ones. What survives is the diffusion plus a
        ! linear Rossby term that rotates the mode's phase rather than changing
        ! its amplitude -- so the measured amplitude decay is the diffusion's,
        ! and the residual disagreement below is the size of everything else.
        !
        ! The time filter is switched off here: it mixes three time levels and
        ! would smear a pure exponential into something with no closed form.
        !
        ! WHAT IS COMPARED, and it is not the per-step decay. Leapfrog carries
        ! two interleaved chains, and the implicit factor is applied over the
        ! full 2 dt from X^(n-1) to X^(n+1) -- so within one chain the
        ! amplitude falls by 1/(1 + 2 dt/tau) every TWO steps, and the apparent
        ! per-step decay is the square root of that. Comparing the per-step
        ! ratio to the per-application factor is the natural mistake and is
        ! wrong by a factor of two in the exponent. The measurement below
        ! therefore walks one chain (even steps only) and takes the geometric
        ! mean of five applications, which also averages out the residual
        ! dynamical coupling that makes individual pairs scatter by ~2%.

        implicit none

        integer, intent(inout) :: nfail

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts
        type(aeros_sht_class), pointer :: s
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(wp), parameter :: dt = 1800.0_wp, tau_h = 6.0_wp
        real(dp) :: a_mid, a_end, measured, predicted, tau, y00
        integer  :: n, k, lm
        integer, parameter :: nstep = 20, nhalf = 10

        write(*,*) ""
        write(*,*) " -- del^6 damping of the l = lmax wave"

        s => pool%sht(1)
        y00 = sqrt(4.0_dp*acos(-1.0_dp))

        call default_par(par, dt=dt, semi_implicit=.TRUE., eps_filter=0.0_wp, tau_diff=tau_h)
        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        call aeros_spec_zero(now%spec)
        do k = 1, nlev
            now%spec%temp(aeros_sht_lm(s,0,0),k) = cmplx(280.0_dp*y00, 0.0_dp, wp_sh)
        end do
        now%spec%lnps(aeros_sht_lm(s,0,0)) = cmplx(log(real(p0,dp))*y00, 0.0_dp, wp_sh)

        lm = aeros_sht_lm(s, s%lmax, 3)
        do k = 1, nlev
            now%spec%vor(lm,k) = cmplx(1.0e-12_dp, 0.0_dp, wp_sh)
        end do

        phis = 0.0_wp
        call aeros_timestep_set_phis(ts, phis)

        ! Measured between steps nhalf and nstep, i.e. after the start-up step
        ! and its h = dt factor have washed out of the chain.
        a_mid = 0.0_dp
        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)
            if (n == nhalf) a_mid = abs(now%spec%vor(lm,nlev/2))
        end do
        a_end = abs(now%spec%vor(lm,nlev/2))

        ! At l = lmax the ratio [l(l+1)/lmax(lmax+1)]^3 is exactly 1, so one
        ! application of the operator multiplies by 1/(1 + 2 dt/tau). There are
        ! (nstep-nhalf)/2 of them along the chain.
        tau       = real(tau_h, dp)*3600.0_dp
        predicted = 1.0_dp/(1.0_dp + 2.0_dp*real(dt,dp)/tau)
        measured  = (a_end/a_mid)**(2.0_dp/real(nstep - nhalf, dp))

        write(*,"(a40,f12.8)")  "   measured factor per 2 dt ", measured
        write(*,"(a40,f12.8)")  "   predicted 1/(1 + 2dt/tau)", predicted
        write(*,"(a40,es12.3)") "   relative difference      ", &
                                    abs(measured - predicted)/predicted

        ! 2% covers the linear coupling to divergence that a vorticity
        ! perturbation cannot avoid exciting; the diffusion factor itself is
        ! exact arithmetic.
        call check(abs(measured - predicted)/predicted < 2.0e-2_dp, &
                    "the l = lmax wave decays at the advertised rate", nfail)

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine test_diffusion

    ! === 4. Conservation =====================================================

    subroutine test_conservation(nfail)

        implicit none

        integer, intent(inout) :: nfail

        write(*,*) ""
        write(*,*) " -- conservation, balanced resting atmosphere (100 steps)"
        call conserve_run(rest=.TRUE., nstep=100, &
                            tol_mass=1.0e-14_dp, tol_ener=1.0e-14_dp, &
                            tol_amom=1.0e-14_dp, assert_all=.TRUE., nfail=nfail)

        write(*,*) ""
        write(*,*) " -- conservation, moving state (400 steps)"
        call conserve_run(rest=.FALSE., nstep=400, &
                            tol_mass=1.0e-4_dp, tol_ener=1.0e-3_dp, &
                            tol_amom=1.0e-2_dp, assert_all=.FALSE., nfail=nfail)

        return

    end subroutine test_conservation

    subroutine conserve_run(rest, nstep, tol_mass, tol_ener, tol_amom, assert_all, nfail)
        ! Integrate and compare the global budgets against their initial values.
        !
        ! `assert_all` distinguishes the two questions. On the resting state
        ! all three integrals are asserted at machine precision, because
        ! nothing there can move them except a leak. On a moving state only
        ! BOUNDEDNESS is asserted: the numbers oscillate with the flow, and a
        ! tight tolerance would be a claim the scheme does not support.

        implicit none

        logical,  intent(in) :: rest
        integer,  intent(in) :: nstep
        real(dp), intent(in) :: tol_mass, tol_ener, tol_amom
        logical,  intent(in) :: assert_all
        integer, intent(inout) :: nfail

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts
        type(aeros_budget_class)   :: b0, b1
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(wp), parameter :: dt = 1800.0_wp
        real(dp) :: elapsed, d_mass, d_ener, d_amom
        integer  :: n

        call default_par(par, dt=dt, semi_implicit=.TRUE.)
        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        if (rest) then
            call rest_state(now, phis, topography=.TRUE.)
        else
            call moving_state(now)
            phis = 0.0_wp
        end if
        call aeros_timestep_set_phis(ts, phis)

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        call aeros_budget_calc(b0, vg, grd, now, phis)

        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)
        end do

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        call aeros_budget_calc(b1, vg, grd, now, phis)

        elapsed = real(nstep, dp)*real(dt, dp)

        d_mass = abs(b1%mass   - b0%mass)  /b0%mass
        d_ener = abs(b1%energy - b0%energy)/b0%energy
        d_amom = abs(b1%angmom - b0%angmom)/b0%angmom

        write(*,"(a40,es12.3)") "   |dM/M| over the integration ", d_mass
        write(*,"(a40,es12.3)") "   |dE/E|                      ", d_ener
        write(*,"(a40,es12.3)") "   |dA/A|                      ", d_amom
        call aeros_budget_report(b1, b0, elapsed, "drift rates")

        if (assert_all) then
            call check(d_mass < tol_mass, "mass is conserved to machine precision", nfail)
            call check(d_ener < tol_ener, "energy is conserved to machine precision", nfail)
            call check(d_amom < tol_amom, "angular momentum is conserved to machine precision", nfail)
        else
            call check(d_mass < tol_mass, "mass stays bounded", nfail)
            call check(d_ener < tol_ener, "energy stays bounded", nfail)
            call check(d_amom < tol_amom, "angular momentum stays bounded", nfail)
        end if

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine conserve_run

    ! === 5. The mass fixer ===================================================

    subroutine test_mass_fixer(nfail)
        ! The fixer has to do three separable things, and getting one right
        ! while getting another wrong would be worse than not having it:
        !
        !   1. RESTORE THE MASS EXACTLY. Not approximately, not eventually --
        !      the correction is closed-form, so machine precision is the only
        !      acceptable answer, and the unfixed run is measured alongside so
        !      that the check is against a drift that actually exists.
        !
        !   2. TOUCH NOTHING ELSE. Only the (0,0) surface-pressure coefficient
        !      may change. A fixer that moved any other wavenumber would be
        !      quietly smoothing the surface pressure field, and after a few
        !      thousand steps nobody would be able to tell what had done it.
        !      Checked bitwise on a single step, where the two trajectories
        !      have not yet had a chance to separate.
        !
        !   3. REPORT WHAT IT DID. `lnr_last` must be exactly the leak the
        !      unfixed step produced, so that switching the fixer on hides
        !      nothing. Also checked on the single step, for the same reason.
        !
        ! The cumulative correction is printed but NOT asserted against the
        ! unfixed drift, because the two are different quantities and expecting
        ! them to match is the mistake this comment exists to prevent. Here the
        ! fixer puts back ~1.9e-5 while the unfixed run drifts ~7.5e-6. That is
        ! not chaos -- this state is not chaotic -- it is that the unfixed run's
        ! mass oscillates as well as drifts, the oscillation largely cancels in
        ! an accumulated total, and the fixer resets it every step and so pays
        ! for it every step. Per step the two agree exactly, which is what
        ! check 3 pins down; over a run they do not, and only an unfixed twin
        ! measures the drift.

        implicit none

        integer, intent(inout) :: nfail

        integer,  parameter :: nstep = 200
        real(dp) :: d_on, d_off, removed, lnr_expect, lnr_got
        real(dp) :: dvor, ddiv, dtemp, dlnps_other, dlnps_00

        write(*,*) ""
        write(*,*) " -- the global mass fixer"

        ! --- 1. the drift, with and without ---------------------------------
        call fixer_run(.FALSE., nstep, d_off)
        call fixer_run(.TRUE.,  nstep, d_on, removed)

        write(*,"(a40,es12.3)") "   |dM/M| over 200 steps, fixer off", d_off
        write(*,"(a40,es12.3)") "   |dM/M| over 200 steps, fixer on ", d_on
        write(*,"(a40,es12.3)") "   what the fixer put back to do it", removed

        call check(d_on < 1.0e-13_dp, "the fixer restores mass to machine precision", nfail)
        call check(d_off > 1.0e3_dp*d_on, &
                    "and there was a drift to remove (the check is not vacuous)", nfail)

        ! --- 2 and 3. one step, side by side --------------------------------
        call fixer_one_step(lnr_expect, lnr_got, dvor, ddiv, dtemp, dlnps_other, dlnps_00)

        write(*,"(a40,es12.3)") "   max |d zeta| from the fixer     ", dvor
        write(*,"(a40,es12.3)") "   max |d D|                       ", ddiv
        write(*,"(a40,es12.3)") "   max |d T|                       ", dtemp
        write(*,"(a40,es12.3)") "   max |d lnps|, (l,m) /= (0,0)    ", dlnps_other
        write(*,"(a40,es12.3)") "   |d lnps| at (0,0)               ", dlnps_00
        write(*,"(a40,es12.3)") "   reported lnr vs the actual leak ", abs(lnr_got - lnr_expect)

        call check(dvor == 0.0_dp .and. ddiv == 0.0_dp .and. dtemp == 0.0_dp, &
                    "the fixer leaves vorticity, divergence and temperature untouched", nfail)
        call check(dlnps_other == 0.0_dp, &
                    "the fixer leaves every surface-pressure wavenumber but (0,0) untouched", nfail)
        call check(dlnps_00 > 0.0_dp, &
                    "the fixer did change the (0,0) coefficient", nfail)
        ! 1e-14 is not slack, it is the floor: the two masses are dp sums over
        ! ~10^4 grid points accumulated in different orders, so they agree to
        ! ~1e-15 relative and lnr inherits that as an absolute error. The leak
        ! itself is ~5e-6 here, so this still pins the reported correction to
        ! eight decimal places of the thing it claims to be.
        call check(abs(lnr_got - lnr_expect) < 1.0e-14_dp, &
                    "the reported correction is exactly the leak it removed", nfail)

        return

    end subroutine test_mass_fixer

    subroutine fixer_run(fixer, nsteps, d_mass, removed)
        ! Integrate the moving state and return the relative mass change, and
        ! optionally what the fixer had to put back to achieve it.

        implicit none

        logical,  intent(in)  :: fixer
        integer,  intent(in)  :: nsteps
        real(dp), intent(out) :: d_mass
        real(dp), intent(out), optional :: removed

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts
        type(aeros_budget_class)   :: b0, b1
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(wp), parameter :: dt = 1800.0_wp
        integer  :: n

        phis = 0.0_wp

        call default_par(par, dt=dt, semi_implicit=.TRUE., mass_fixer=fixer)
        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        call moving_state(now)
        call aeros_timestep_set_phis(ts, phis)

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        call aeros_budget_calc(b0, vg, grd, now, phis)

        do n = 1, nsteps
            call aeros_timestep_step(ts, pool, vg, grd, now)
        end do

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        call aeros_budget_calc(b1, vg, grd, now, phis)

        d_mass = abs(b1%mass - b0%mass)/b0%mass
        if (present(removed)) removed = exp(-ts%lnr_cum) - 1.0_dp

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine fixer_run

    subroutine fixer_one_step(lnr_expect, lnr_got, dvor, ddiv, dtemp, dlnps_other, dlnps_00)
        ! Take one step from the same state twice, fixer off then on, and
        ! compare the two spectral states coefficient by coefficient.
        !
        ! One step, because that is the only place the comparison is clean: the
        ! two runs share an initial state, and the fixer acts after the state
        ! update, so everything the dynamics produced must be bit-identical.
        ! By step two the corrected run is on a different trajectory and a
        ! difference no longer says anything about what the fixer touched.

        implicit none

        real(dp), intent(out) :: lnr_expect, lnr_got
        real(dp), intent(out) :: dvor, ddiv, dtemp, dlnps_other, dlnps_00

        type(aeros_param_class)    :: par
        type(aeros_timestep_class) :: ts
        type(aeros_budget_class)   :: b0, b1
        type(aeros_state_class)    :: s_off, s_on
        type(aeros_sht_class), pointer :: sh
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(wp), parameter :: dt = 1800.0_wp
        real(dp) :: m0, m1
        integer  :: k, lm, lm00

        sh => pool%sht(1)
        lm00 = aeros_sht_lm(sh, 0, 0)
        phis = 0.0_wp

        ! -- fixer off
        call default_par(par, dt=dt, semi_implicit=.TRUE., mass_fixer=.FALSE.)
        call aeros_state_alloc(s_off, grd, sh%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)
        call moving_state(s_off)
        call aeros_timestep_set_phis(ts, phis)

        call aeros_timestep_diagnose(ts, pool, vg, grd, s_off)
        call aeros_budget_calc(b0, vg, grd, s_off, phis)
        m0 = b0%mass

        call aeros_timestep_step(ts, pool, vg, grd, s_off)

        call aeros_timestep_diagnose(ts, pool, vg, grd, s_off)
        call aeros_budget_calc(b1, vg, grd, s_off, phis)
        m1 = b1%mass

        ! What the fixer should report having removed. Any systematic
        ! difference between this quadrature and the fixer's own cancels in the
        ! ratio, so the comparison below is exact rather than approximate.
        lnr_expect = log(m0/m1)
        call aeros_timestep_end(ts)

        ! -- fixer on, same initial state
        call default_par(par, dt=dt, semi_implicit=.TRUE., mass_fixer=.TRUE.)
        call aeros_state_alloc(s_on, grd, sh%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)
        call moving_state(s_on)
        call aeros_timestep_set_phis(ts, phis)

        call aeros_timestep_step(ts, pool, vg, grd, s_on)

        lnr_got = ts%lnr_last

        ! -- compare
        dvor = 0.0_dp; ddiv = 0.0_dp; dtemp = 0.0_dp
        do k = 1, nlev
            do lm = 1, sh%nlm
                dvor  = max(dvor,  abs(s_on%spec%vor(lm,k)  - s_off%spec%vor(lm,k)))
                ddiv  = max(ddiv,  abs(s_on%spec%div(lm,k)  - s_off%spec%div(lm,k)))
                dtemp = max(dtemp, abs(s_on%spec%temp(lm,k) - s_off%spec%temp(lm,k)))
            end do
        end do

        dlnps_other = 0.0_dp
        do lm = 1, sh%nlm
            if (lm == lm00) cycle
            dlnps_other = max(dlnps_other, abs(s_on%spec%lnps(lm) - s_off%spec%lnps(lm)))
        end do
        dlnps_00 = abs(s_on%spec%lnps(lm00) - s_off%spec%lnps(lm00))

        call aeros_timestep_end(ts)
        call aeros_state_end(s_off)
        call aeros_state_end(s_on)

        return

    end subroutine fixer_one_step

    subroutine check(ok, label, nfail)

        implicit none

        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail

        if (ok) then
            write(*,*) "  ok   : ", trim(label)
        else
            write(*,*) "  FAIL : ", trim(label)
            nfail = nfail + 1
        end if

        return

    end subroutine check

end program test_timestep
