program test_ncycle
    ! Acceptance test for the two-time-level integrator: the semi-implicit
    ! Lorenz N-cycle (aeros_timestep, "The Lorenz N-cycle").
    !
    ! 1. A BALANCED RESTING ATMOSPHERE STAYS AT REST, over topography, 100
    !    steps: the weighted accumulator and the correction solve called with
    !    X^(n-1) = X^n manufacture no motion.
    !
    ! 2. STABILITY. The moving state at a 30-minute step: semi-implicit
    !    N-cycle stable; explicit N-cycle diverges (gravity-wave limit);
    !    explicit at a step below the limit stable again -- the same triple
    !    tests/test_timestep.f90 runs for the leapfrog.
    !
    ! 3. CONVERGENCE. Adiabatic dry dynamics from the moving state over a fixed
    !    interval at dt and dt/2 against a dt/8 reference: the error must drop
    !    by at least a factor ~3.5, i.e. the scheme is at least second order
    !    on the nonlinear problem (the implicit part is Crank-Nicolson).
    !    Reported for every variant; asserted for the default.
    !
    ! 4. RESTART IS BIT-EXACT MID-CYCLE. A continuous run against one checkpointed
    !    inside a cycle (the accumulator G travels in the file).
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, R_d, p0, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_spec_class, aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_tendency
    use aeros_timestep

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg

    integer :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_ncycle:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    call test_rest(nfail)
    call test_stability(nfail)
    call test_convergence(nfail)
    call test_restart_midcycle(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_ncycle:: PASS"
    else
        write(*,*) "test_ncycle:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine default_par(par, dt, semi_implicit, scheme, variant)

        implicit none

        type(aeros_param_class), intent(out) :: par
        real(wp), intent(in) :: dt
        logical,  intent(in) :: semi_implicit
        character(len=*), intent(in), optional :: scheme, variant

        par%trunc    = trunc
        par%nlon     = -1
        par%nlat     = -1
        par%nlev     = nlev
        par%nthreads = -1

        par%dt            = dt
        par%semi_implicit = semi_implicit
        par%held_suarez   = .FALSE.
        par%eps_filter    = 0.1_wp
        par%raw_alpha     = 0.53_wp
        par%ndiff         = 8
        par%tau_diff      = 4.0_wp
        par%mass_fixer    = .FALSE.
        par%time_scheme    = "ncycle"
        par%ncycle         = 4
        par%ncycle_variant = "A"
        if (present(scheme))  par%time_scheme    = scheme
        if (present(variant)) par%ncycle_variant = variant

        return

    end subroutine default_par

    subroutine rest_state(now, phis)
        ! Isothermal atmosphere at rest in exact discrete hydrostatic balance
        ! with its surface geopotential (tests/test_timestep.f90).

        implicit none

        type(aeros_state_class), intent(inout) :: now
        real(wp), intent(out) :: phis(:,:)

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
        lm = aeros_sht_lm(s,2,1); now%spec%lnps(lm) = cmplx(-0.06_dp,  0.03_dp, wp_sh)
        lm = aeros_sht_lm(s,3,0); now%spec%lnps(lm) = cmplx( 0.05_dp,  0.0_dp,  wp_sh)
        lm = aeros_sht_lm(s,5,3); now%spec%lnps(lm) = cmplx( 0.02_dp, -0.01_dp, wp_sh)

        call aeros_sht_synthesis(s, now%spec%lnps, lnps_g)
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                phis(i,j) = R_d*tiso*(log(p0) - lnps_g(i,j))
            end do
        end do

        return

    end subroutine rest_state

    subroutine moving_state(now)
        ! Sheared rotational flow over a flat surface (tests/test_timestep.f90).

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
        real(wp) :: phis(grd%nlon,grd%nlat), ps0(grd%nlon,grd%nlat)
        real(dp) :: umax, dps, tmax
        integer  :: n

        write(*,*) ""
        write(*,*) " -- a balanced resting atmosphere, over topography, 100 substeps"

        call default_par(par, dt=1800.0_wp, semi_implicit=.TRUE.)
        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        call rest_state(now, phis)
        call aeros_timestep_set_phis(ts, phis)
        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        ps0 = now%ps

        do n = 1, 100
            call aeros_timestep_step(ts, pool, vg, grd, now)
        end do
        call aeros_timestep_diagnose(ts, pool, vg, grd, now)

        umax = max(maxval(abs(real(now%u,dp))), maxval(abs(real(now%v,dp))))
        dps  = maxval(abs(real(now%ps,dp) - real(ps0,dp)))
        tmax = maxval(abs(real(now%temp_g,dp) - 280.0_dp))

        write(*,"(a40,es12.3,a)") "   max |wind| after 100 steps ", umax, " m s-1"
        write(*,"(a40,es12.3,a)") "   max |dp_s|                 ", dps,  " Pa"
        write(*,"(a40,es12.3,a)") "   max |T - 280|              ", tmax, " K"

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
        write(*,*) " -- stability at and below the gravity-wave limit (200 substeps)"

        call integrate_growth(1800.0_wp, .TRUE.,  200, g_si, d_si)
        call integrate_growth(1800.0_wp, .FALSE., 200, g_ex, d_ex)
        call integrate_growth(150.0_wp,  .FALSE., 200, g_ex_short, d_ex_short)

        write(*,"(a40,es12.3,l3)") "   semi-implicit, 1800 s: growth ", g_si, d_si
        write(*,"(a40,es12.3,l3)") "   explicit, 1800 s:      growth ", g_ex, d_ex
        write(*,"(a40,es12.3,l3)") "   explicit, 150 s:       growth ", g_ex_short, d_ex_short

        call check(.not. d_si .and. g_si < 10.0_dp, &
                    "semi-implicit N-cycle is stable at a 30 min step", nfail)
        call check(d_ex, "explicit N-cycle diverges at a 30 min step", nfail)
        call check(.not. d_ex_short .and. g_ex_short < 10.0_dp, &
                    "explicit N-cycle is stable below the gravity-wave limit", nfail)

        return

    end subroutine test_stability

    subroutine integrate_growth(dt, semi_implicit, nstep, growth, diverged)

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

        z0 = maxval(abs(now%spec%vor)); growth = 1.0_dp; diverged = .FALSE.
        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)
            if (.not. all(now%spec%vor == now%spec%vor)) then
                diverged = .TRUE.; exit
            end if
            z = maxval(abs(now%spec%vor)); growth = max(growth, z/z0)
            if (growth > 1.0e3_dp) then
                diverged = .TRUE.; exit
            end if
        end do

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine integrate_growth

    ! === 3. Convergence ======================================================
    subroutine test_convergence(nfail)

        implicit none

        integer, intent(inout) :: nfail

        real(wp), parameter :: period = 6.0_wp*3600.0_wp    ! integrate 6 h
        real(wp), parameter :: dt0 = 1200.0_wp
        character(len=4), parameter :: variants(4) = [character(len=4) :: "A", "B", "AB", "ABBA"]
        type(aeros_spec_class) :: ref, x1, x2
        real(dp) :: e1, e2, ratio, ratio_lf
        integer  :: iv

        write(*,*) ""
        write(*,*) " -- convergence, adiabatic moving state, 6 h: error at dt vs dt/2"
        write(*,*) "    against a dt/8 reference (ratio ~4 = second order)"

        call aeros_spec_alloc(ref, pool%sht(1)%nlm, nlev)
        call aeros_spec_alloc(x1,  pool%sht(1)%nlm, nlev)
        call aeros_spec_alloc(x2,  pool%sht(1)%nlm, nlev)

        do iv = 1, size(variants)
            call integrate_to(period, dt0/8.0_wp, "ncycle", trim(variants(iv)), ref)
            call integrate_to(period, dt0,        "ncycle", trim(variants(iv)), x1)
            call integrate_to(period, dt0/2.0_wp, "ncycle", trim(variants(iv)), x2)
            e1 = spec_dist(x1, ref); e2 = spec_dist(x2, ref)
            ratio = e1/max(e2, tiny(1.0_dp))
            write(*,"(a,a4,a,es11.3,a,es11.3,a,f7.2)") "   ncycle ", variants(iv), &
                    ": err(dt) ", e1, "  err(dt/2) ", e2, "  ratio ", ratio
            if (iv == 1) call check(ratio > 3.5_dp .and. e1 == e1, &
                    "N-cycle (variant A) converges at least second order", nfail)
        end do

        ! The filtered leapfrog, for the record (reported, not asserted).
        call integrate_to(period, dt0/8.0_wp, "leapfrog", "A", ref)
        call integrate_to(period, dt0,        "leapfrog", "A", x1)
        call integrate_to(period, dt0/2.0_wp, "leapfrog", "A", x2)
        e1 = spec_dist(x1, ref); e2 = spec_dist(x2, ref)
        ratio_lf = e1/max(e2, tiny(1.0_dp))
        write(*,"(a,es11.3,a,es11.3,a,f7.2)") "   leapfrog  : err(dt) ", e1, &
                "  err(dt/2) ", e2, "  ratio ", ratio_lf

        call aeros_spec_end(ref); call aeros_spec_end(x1); call aeros_spec_end(x2)

        return

    end subroutine test_convergence

    subroutine integrate_to(period, dt, scheme, variant, out)
        ! Adiabatic, undiffused (tau_diff huge) integration of the moving
        ! state over `period` at step dt; the final spectral state in `out`.

        implicit none

        real(wp), intent(in) :: period, dt
        character(len=*), intent(in) :: scheme, variant
        type(aeros_spec_class), intent(inout) :: out

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts
        real(wp) :: phis(grd%nlon,grd%nlat)
        integer  :: n, nstep

        call default_par(par, dt=dt, semi_implicit=.TRUE., scheme=scheme, variant=variant)
        par%tau_diff = 1.0e6_wp          ! effectively no diffusion: a pure time-scheme test
        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        call moving_state(now)
        phis = 0.0_wp
        call aeros_timestep_set_phis(ts, phis)

        nstep = nint(period/dt)
        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)
        end do
        call aeros_spec_copy(out, now%spec)

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine integrate_to

    real(dp) function spec_dist(a, b) result(d)
        ! RMS distance between two spectral states, vorticity + divergence
        ! (scaled by 1e-5 s-1), temperature (K) and ln p_s, all coefficients.

        implicit none

        type(aeros_spec_class), intent(in) :: a, b

        real(dp) :: sv, sd, st, sp
        integer  :: n

        sv = sum(abs(a%vor  - b%vor )**2)/1.0e-10_dp
        sd = sum(abs(a%div  - b%div )**2)/1.0e-10_dp
        st = sum(abs(a%temp - b%temp)**2)
        sp = sum(abs(a%lnps - b%lnps)**2)
        n  = size(a%vor) + size(a%div) + size(a%temp) + size(a%lnps)
        d  = sqrt((sv + sd + st + sp)/real(n, dp))

        return

    end function spec_dist

    ! === 4. Restart mid-cycle ================================================
    subroutine test_restart_midcycle(nfail)

        implicit none

        integer, intent(inout) :: nfail

        character(len=*), parameter :: fnm = "test_ncycle_restart_chk.nc"
        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now_a, now_b
        type(aeros_timestep_class) :: ts_a, ts_b
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(dp) :: dmax, time
        integer  :: n

        write(*,*) ""
        write(*,*) " -- restart inside a cycle: 10 substeps continuous vs 6 + checkpoint + 4"

        call default_par(par, dt=1800.0_wp, semi_implicit=.TRUE.)
        phis = 0.0_wp

        ! A: continuous.
        call aeros_state_alloc(now_a, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts_a, par, pool, grd, vg)
        call moving_state(now_a)
        call aeros_timestep_set_phis(ts_a, phis)
        do n = 1, 10
            call aeros_timestep_step(ts_a, pool, vg, grd, now_a)
        end do

        ! B: 6 substeps, checkpoint (substep 2 of the second cycle), fresh
        ! integrator, resume, 4 more.
        call aeros_state_alloc(now_b, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts_b, par, pool, grd, vg)
        call moving_state(now_b)
        call aeros_timestep_set_phis(ts_b, phis)
        do n = 1, 6
            call aeros_timestep_step(ts_b, pool, vg, grd, now_b)
        end do
        time = 6.0_dp*1800.0_dp
        call aeros_timestep_write_restart(ts_b, now_b, time, fnm)
        call aeros_timestep_end(ts_b)
        call aeros_timestep_init(ts_b, par, pool, grd, vg)
        call aeros_timestep_set_phis(ts_b, phis)
        call aeros_timestep_read_restart(ts_b, now_b, time, fnm)
        do n = 1, 4
            call aeros_timestep_step(ts_b, pool, vg, grd, now_b)
        end do

        dmax = max(maxval(abs(now_a%spec%vor  - now_b%spec%vor)), &
                   maxval(abs(now_a%spec%div  - now_b%spec%div)), &
                   maxval(abs(now_a%spec%temp - now_b%spec%temp)), &
                   maxval(abs(now_a%spec%lnps - now_b%spec%lnps)))
        write(*,"(a40,es12.3)") "   max |continuous - restarted|  ", dmax
        call check(dmax == 0.0_dp, "restart inside a cycle is bit-exact", nfail)

        call aeros_timestep_end(ts_a); call aeros_timestep_end(ts_b)
        call aeros_state_end(now_a);   call aeros_state_end(now_b)

        return

    end subroutine test_restart_midcycle

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

end program test_ncycle
