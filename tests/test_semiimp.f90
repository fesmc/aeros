program test_semiimp
    ! Acceptance test for the semi-implicit gravity-wave solve.
    !
    ! Two things have to be true, and they are independent.
    !
    ! 1. L IS THE LINEARIZATION OF R. The correction form only cancels -- and
    !    the scheme is only reliably stable -- if the operator treated
    !    implicitly is a term-for-term linearization of the discrete nonlinear
    !    right-hand side, not merely of the continuous equations it
    !    approximates. So L is checked by FINITE DIFFERENCE against
    !    aeros_tendency itself: perturb the reference state, evaluate the real
    !    nonlinear tendency, and compare the difference quotient with L applied
    !    to the same perturbation.
    !
    !    This is the test that decides the coefficient of ln(p_s) in L_D. The
    !    nonlinear pressure-gradient force carries c_k (aeros_tendency's cfac),
    !    so R_d T_ref c_k is the natural-looking choice, and it is wrong: a
    !    surface-pressure perturbation also moves the half levels, so the
    !    geopotential responds by R_d T_ref (1 - c_k), and the two sum to
    !    R_d T_ref. With c_k restored the lnps check below fails by ~100% in
    !    the upper layers of a hybrid coordinate, where c_k -> 0.
    !
    !    Only the components with a clean linearization are asserted. A
    !    divergence perturbation also spins up the CORIOLIS term, which is
    !    linear in D but deliberately NOT part of L -- semi-implicit treats
    !    gravity waves, not rotation -- so d(D)/dt is reported for that case
    !    and not asserted.
    !
    ! 2. THE SOLVE IMPLEMENTS THE SCHEME. aeros_semiimp_step eliminates T and
    !    ln(p_s), solves an nlev x nlev system per degree and back-substitutes.
    !    None of that is visible in its output, so the result is checked
    !    against the scheme's DEFINING equation,
    !
    !        X^new = X^old + h R + h [ L(X_bar) - L(X^now) ],  X_bar = (X^new+X^old)/2
    !
    !    re-evaluated with the public operator. That exercises the LU, the
    !    elimination algebra and the back-substitution together, at machine
    !    precision, for arbitrary inputs -- and it duplicates none of the
    !    module's internals, so it is a real check rather than a restatement.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, R_d, cp_d, p0, &
                                aeros_grid_class, aeros_spec_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_tendency
    use aeros_semiimp

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 20

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)     :: grd
    type(aeros_vgrid_class)    :: vg
    type(aeros_semiimp_class)  :: si

    integer  :: nfail
    real(wp) :: c_sigma, c_hybrid, c_lamb

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_semiimp:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    ! === 1. The gravity-wave spectrum ========================================
    !
    ! For an isothermal atmosphere the external (Lamb) mode travels at
    ! sqrt(gamma R T) with gamma = cp/cv, and that is a property of the
    ! ATMOSPHERE, not of the coordinate used to discretize it. So the same
    ! reference state on pure sigma and on the shipped hybrid must give the
    ! same speed -- which is a sharper test of the operator than either number
    ! against theory, since a coordinate-dependent error cancels in neither.
    write(*,*) ""
    write(*,*) " -- gravity-wave spectrum"

    c_lamb = sqrt(cp_d/(cp_d - R_d)*R_d*300.0_wp)

    call aeros_vgrid_init(vg, nlev, stretch_a=1.0_wp, sigma_t=0.0_wp, &
                            p_top=1000.0_wp, t_ref=300.0_wp)
    call aeros_semiimp_init(si, vg, trunc, 1800.0_wp)
    c_sigma = aeros_semiimp_gwspeed(si)
    call aeros_semiimp_print(si)
    call aeros_semiimp_end(si)
    call aeros_vgrid_end(vg)

    call aeros_vgrid_init(vg, nlev, t_ref=300.0_wp)
    call aeros_semiimp_init(si, vg, trunc, 1800.0_wp)
    c_hybrid = aeros_semiimp_gwspeed(si)
    call aeros_semiimp_end(si)
    call aeros_vgrid_end(vg)

    write(*,"(a40,f9.2,a)") "   analytic sqrt(gamma R T) ", c_lamb,   " m s-1"
    write(*,"(a40,f9.2,a)") "   pure sigma               ", c_sigma,  " m s-1"
    write(*,"(a40,f9.2,a)") "   hybrid (shipped default) ", c_hybrid, " m s-1"
    write(*,"(a40,es12.3)") "   hybrid vs sigma, relative", &
                                abs(c_hybrid - c_sigma)/c_sigma

    ! 5% covers the ~1% of the atmosphere's mass above a 10 hPa top plus the
    ! L20 vertical discretization of the mode.
    call check(abs(c_hybrid - c_lamb)/c_lamb < 0.05_wp, &
                "external mode is within 5% of sqrt(gamma R T)", nfail)
    call check(abs(c_hybrid - c_sigma)/c_sigma < 0.01_wp, &
                "external mode does not depend on the coordinate", nfail)

    ! === 2. L against a finite-differenced R ==================================
    call aeros_vgrid_init(vg, nlev)              ! the shipped default
    call test_linearization(pool, grd, vg, nfail)
    call aeros_vgrid_end(vg)

    ! === 3. The solve against the scheme it implements =======================
    call aeros_vgrid_init(vg, nlev)
    call test_solve(pool, grd, vg, nfail)
    call aeros_vgrid_end(vg)

    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_semiimp:: PASS"
    else
        write(*,*) "test_semiimp:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine base_state(sht, vg, spec)
        ! The linearization point itself: isothermal at t_ref, uniform surface
        ! pressure at ps_ref, motionless. It is simultaneously a rest state
        ! over a flat surface, so the nonlinear tendency vanishes there and the
        ! difference quotients below need no base-state subtraction to be
        ! meaningful -- though they subtract it anyway.

        implicit none

        type(aeros_sht_class),  intent(in)    :: sht
        type(aeros_vgrid_class), intent(in)   :: vg
        type(aeros_spec_class), intent(inout) :: spec

        real(dp) :: y00
        integer  :: k

        call aeros_spec_zero(spec)

        ! Orthonormal harmonics: Y_00 = 1/sqrt(4 pi), so a uniform field of
        ! value X has coefficient X*sqrt(4 pi).
        y00 = sqrt(4.0_dp*acos(-1.0_dp))

        do k = 1, vg%nlev
            spec%temp(aeros_sht_lm(sht,0,0),k) = cmplx(vg%t_ref(k)*y00, 0.0_dp, wp_sh)
        end do
        spec%lnps(aeros_sht_lm(sht,0,0)) = cmplx(log(vg%ps_ref)*y00, 0.0_dp, wp_sh)

        return

    end subroutine base_state

    subroutine test_linearization(pool, grd, vg, nfail)

        implicit none

        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_vgrid_class), intent(in) :: vg
        integer, intent(inout) :: nfail

        type(aeros_sht_class), pointer :: s
        type(aeros_spec_class) :: base, pert
        type(aeros_work_class) :: wrk
        type(aeros_tend_class) :: t0, t1
        type(aeros_semiimp_class) :: si

        complex(wp_sh), allocatable :: ldiv0(:,:), ltemp0(:,:), lps0(:)
        complex(wp_sh), allocatable :: ldiv1(:,:), ltemp1(:,:), lps1(:)

        integer  :: nlm, k, lm
        real(dp) :: e_div, e_temp, e_lps

        s   => pool%sht(1)
        nlm =  s%nlm

        call aeros_spec_alloc(base,  nlm, vg%nlev)
        call aeros_spec_alloc(pert,  nlm, vg%nlev)
        call aeros_work_alloc(wrk, grd%nlon, grd%nlat, vg%nlev)
        call aeros_tend_alloc(t0, nlm, vg%nlev)
        call aeros_tend_alloc(t1, nlm, vg%nlev)
        call aeros_semiimp_init(si, vg, s%lmax, 1800.0_wp)

        allocate(ldiv0(nlm,vg%nlev), ltemp0(nlm,vg%nlev), lps0(nlm))
        allocate(ldiv1(nlm,vg%nlev), ltemp1(nlm,vg%nlev), lps1(nlm))

        write(*,*) ""
        write(*,*) " -- L against a finite-differenced nonlinear tendency"

        wrk%phis = 0.0_wp

        call base_state(s, vg, base)
        call aeros_tendency_calc(pool, vg, grd, base, wrk, t0)
        call aeros_semiimp_linear(si, s, base, ldiv0, ltemp0, lps0)

        write(*,"(a40,es12.3,a)") "   base state |d D/dt| ", &
                                    maxval(abs(t0%div)), " s-2"
        call check(maxval(abs(t0%div)) < 1.0e-18_dp, &
                    "the linearization point is a rest state", nfail)

        ! --- 2a. Temperature. Perturb one (l,m) at EVERY level, with a profile,
        ! so the whole upper-triangular G is exercised in one evaluation.
        call aeros_spec_copy(pert, base)
        do k = 1, vg%nlev
            lm = aeros_sht_lm(s,3,2)
            pert%temp(lm,k) = pert%temp(lm,k) &
                                + cmplx(1.0e-2_dp*cos(real(k,dp)), &
                                        0.6e-2_dp*sin(real(2*k,dp)), wp_sh)
        end do

        call aeros_tendency_calc(pool, vg, grd, pert, wrk, t1)
        call aeros_semiimp_linear(si, s, pert, ldiv1, ltemp1, lps1)

        e_div  = reldiff(t1%div  - t0%div,  ldiv1  - ldiv0)
        e_temp = maxval(abs(t1%temp - t0%temp))
        e_lps  = maxval(abs(t1%lnps - t0%lnps))

        write(*,"(a40,es12.3)") "   T pert: d D/dt vs L_D ", e_div
        write(*,"(a40,es12.3)") "   T pert: |d T/dt|      ", e_temp
        write(*,"(a40,es12.3)") "   T pert: |d lnps/dt|   ", e_lps
        call check(e_div  < 1.0e-6_dp,  "L_D reproduces the response to T", nfail)
        call check(e_temp < 1.0e-15_dp, "a T perturbation alone does not move T", nfail)
        call check(e_lps  < 1.0e-18_dp, "a T perturbation alone does not move lnps", nfail)

        ! --- 2b. Surface pressure. THE coefficient test.
        call aeros_spec_copy(pert, base)
        lm = aeros_sht_lm(s,3,2)
        pert%lnps(lm) = pert%lnps(lm) + cmplx(2.0e-5_dp, -1.0e-5_dp, wp_sh)

        call aeros_tendency_calc(pool, vg, grd, pert, wrk, t1)
        call aeros_semiimp_linear(si, s, pert, ldiv1, ltemp1, lps1)

        e_div  = reldiff(t1%div - t0%div, ldiv1 - ldiv0)
        e_temp = maxval(abs(t1%temp - t0%temp))

        write(*,"(a40,es12.3)") "   P pert: d D/dt vs L_D ", e_div
        write(*,"(a40,es12.3)") "   P pert: |d T/dt|      ", e_temp
        call check(e_div  < 1.0e-6_dp,  "L_D reproduces the response to ln(p_s)", nfail)
        call check(e_temp < 1.0e-15_dp, "a lnps perturbation alone does not move T", nfail)

        ! --- 2c. Divergence. Checks tau and nu; d(D)/dt is reported only,
        ! because it also carries the explicitly-treated Coriolis term.
        call aeros_spec_copy(pert, base)
        do k = 1, vg%nlev
            lm = aeros_sht_lm(s,3,2)
            pert%div(lm,k) = pert%div(lm,k) &
                                + cmplx(1.0e-9_dp*cos(real(k,dp)), &
                                        0.7e-9_dp*sin(real(3*k,dp)), wp_sh)
        end do

        call aeros_tendency_calc(pool, vg, grd, pert, wrk, t1)
        call aeros_semiimp_linear(si, s, pert, ldiv1, ltemp1, lps1)

        e_temp = reldiff(t1%temp - t0%temp, ltemp1 - ltemp0)
        e_lps  = reldiff1(t1%lnps - t0%lnps, lps1 - lps0)
        ! L_D is identically zero for a divergence perturbation, so this is an
        ! ABSOLUTE tendency, not a ratio: it is the explicitly-treated Coriolis
        ! response, printed to show it is present and small rather than absent.
        e_div  = maxval(abs(t1%div - t0%div))

        write(*,"(a40,es12.3)")   "   D pert: d T/dt vs L_T ", e_temp
        write(*,"(a40,es12.3)")   "   D pert: d lnps/dt vs L_P ", e_lps
        write(*,"(a40,es12.3,a)") "   D pert: Coriolis d D/dt ", e_div, " s-2 (not asserted)"
        call check(e_temp < 1.0e-6_dp, "L_T reproduces the response to D", nfail)
        call check(e_lps  < 1.0e-9_dp, "L_P reproduces the response to D", nfail)

        deallocate(ldiv0, ltemp0, lps0, ldiv1, ltemp1, lps1)
        call aeros_semiimp_end(si)
        call aeros_tend_end(t0); call aeros_tend_end(t1)
        call aeros_work_end(wrk)
        call aeros_spec_end(base); call aeros_spec_end(pert)

        return

    end subroutine test_linearization

    subroutine test_solve(pool, grd, vg, nfail)
        ! Does aeros_semiimp_step solve the equation it says it solves?

        implicit none

        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_vgrid_class), intent(in) :: vg
        integer, intent(inout) :: nfail

        type(aeros_sht_class), pointer :: s
        type(aeros_spec_class)    :: old, now, new, bar
        type(aeros_tend_class)    :: tnd
        type(aeros_semiimp_class) :: si

        complex(wp_sh), allocatable :: ldb(:,:), ltb(:,:), lpb(:)
        complex(wp_sh), allocatable :: ldn(:,:), ltn(:,:), lpn(:)

        real(wp), parameter :: h = 3600.0_wp    ! 2*dt at dt = 1800 s
        integer  :: nlm, k, lm, l
        real(dp) :: amp, r_div, r_temp, r_lps

        s   => pool%sht(1)
        nlm =  s%nlm

        call aeros_spec_alloc(old, nlm, vg%nlev)
        call aeros_spec_alloc(now, nlm, vg%nlev)
        call aeros_spec_alloc(new, nlm, vg%nlev)
        call aeros_spec_alloc(bar, nlm, vg%nlev)
        call aeros_tend_alloc(tnd, nlm, vg%nlev)
        call aeros_semiimp_init(si, vg, s%lmax, h)

        allocate(ldb(nlm,vg%nlev), ltb(nlm,vg%nlev), lpb(nlm))
        allocate(ldn(nlm,vg%nlev), ltn(nlm,vg%nlev), lpn(nlm))

        write(*,*) ""
        write(*,*) " -- the solve against the scheme's defining equation"

        ! Arbitrary, deterministic, and with realistic magnitudes: old and now
        ! genuinely different, and a non-zero tendency, so nothing in the
        ! elimination can cancel by accident.
        call base_state(s, vg, old)
        call base_state(s, vg, now)
        do k = 1, vg%nlev
            do lm = 1, nlm
                l = s%l_of_lm(lm)
                if (l < 1 .or. l > 8) cycle
                amp = 1.0e-6_dp/real(l*l, dp)
                old%div(lm,k)  = cmplx(amp*cos(real(l+2*k,dp)), amp*sin(real(3*l+k,dp)), wp_sh)
                now%div(lm,k)  = cmplx(amp*cos(real(2*l+k,dp)), amp*sin(real(l+3*k,dp)), wp_sh)
                old%temp(lm,k) = old%temp(lm,k) + cmplx(0.5_dp*cos(real(l*k,dp)), 0.0_dp, wp_sh)
                now%temp(lm,k) = now%temp(lm,k) + cmplx(0.4_dp*sin(real(l+k,dp)), 0.0_dp, wp_sh)
                tnd%div(lm,k)  = cmplx(1.0e-11_dp*sin(real(l+k,dp)), 0.0_dp, wp_sh)
                tnd%temp(lm,k) = cmplx(1.0e-4_dp*cos(real(l-k,dp)), 0.0_dp, wp_sh)
            end do
        end do
        do lm = 1, nlm
            l = s%l_of_lm(lm)
            if (l < 1 .or. l > 8) cycle
            old%lnps(lm) = cmplx(1.0e-3_dp*cos(real(l,dp)), 0.0_dp, wp_sh)
            now%lnps(lm) = cmplx(0.9e-3_dp*sin(real(l,dp)), 0.0_dp, wp_sh)
            tnd%lnps(lm) = cmplx(1.0e-8_dp*cos(real(2*l,dp)), 0.0_dp, wp_sh)
        end do

        call aeros_semiimp_step(si, s, old, now, tnd, new)

        ! X_bar = (X^new + X^old)/2, then re-evaluate the scheme.
        bar%div  = 0.5_wp*(new%div  + old%div)
        bar%temp = 0.5_wp*(new%temp + old%temp)
        bar%lnps = 0.5_wp*(new%lnps + old%lnps)

        call aeros_semiimp_linear(si, s, bar, ldb, ltb, lpb)
        call aeros_semiimp_linear(si, s, now, ldn, ltn, lpn)

        r_div  = maxval(abs(new%div  - old%div  - h*(tnd%div  + ldb - ldn))) &
                    /max(maxval(abs(new%div  - old%div)),  tiny(1.0_dp))
        r_temp = maxval(abs(new%temp - old%temp - h*(tnd%temp + ltb - ltn))) &
                    /max(maxval(abs(new%temp - old%temp)), tiny(1.0_dp))
        r_lps  = maxval(abs(new%lnps - old%lnps - h*(tnd%lnps + lpb - lpn))) &
                    /max(maxval(abs(new%lnps - old%lnps)), tiny(1.0_dp))

        write(*,"(a40,es12.3)") "   divergence equation residual ", r_div
        write(*,"(a40,es12.3)") "   temperature equation residual ", r_temp
        write(*,"(a40,es12.3)") "   ln(p_s) equation residual ", r_lps

        call check(r_div  < 1.0e-10_dp, "the divergence solve satisfies its equation", nfail)
        call check(r_temp < 1.0e-12_dp, "temperature back-substitution is consistent", nfail)
        call check(r_lps  < 1.0e-12_dp, "ln(p_s) back-substitution is consistent", nfail)

        deallocate(ldb, ltb, lpb, ldn, ltn, lpn)
        call aeros_semiimp_end(si)
        call aeros_tend_end(tnd)
        call aeros_spec_end(old); call aeros_spec_end(now)
        call aeros_spec_end(new); call aeros_spec_end(bar)

        return

    end subroutine test_solve

    real(dp) function reldiff(a, b) result(e)
        ! max|a-b| / max|b|, the natural measure for "does a reproduce b".

        implicit none

        complex(wp_sh), intent(in) :: a(:,:), b(:,:)

        real(dp) :: scale

        scale = maxval(abs(b))
        if (scale <= 0.0_dp) then
            e = maxval(abs(a))
        else
            e = maxval(abs(a - b))/scale
        end if

        return

    end function reldiff

    real(dp) function reldiff1(a, b) result(e)

        implicit none

        complex(wp_sh), intent(in) :: a(:), b(:)

        real(dp) :: scale

        scale = maxval(abs(b))
        if (scale <= 0.0_dp) then
            e = maxval(abs(a))
        else
            e = maxval(abs(a - b))/scale
        end if

        return

    end function reldiff1

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

end program test_semiimp
