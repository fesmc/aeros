program test_vertical
    ! Acceptance test for the hybrid sigma-pressure coordinate and the discrete
    ! hydrostatic relation.
    !
    ! Two of these checks are EXACT rather than tolerant, and deliberately so:
    !
    !   sum(dp) = p_s        mass. If this is only approximate, the model leaks
    !                        mass at a rate set by round-off times the number of
    !                        timesteps -- which over 10^5 yr is not round-off.
    !
    !   sigma_t = 0 => A = 0 An exactly-pure-sigma configuration is what
    !                        Held-Suarez (M1.5) is defined on. "Close to sigma"
    !                        would make the benchmark comparison approximate for
    !                        a reason that has nothing to do with the dynamics.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, R_d, p0, MV
    use aeros_vertical

    implicit none

    type(aeros_vgrid_class) :: vg

    integer, parameter :: nlev = 20

    real(wp) :: tol
    integer  :: nfail

    nfail = 0
    tol   = 20.0_wp*epsilon(1.0_wp)

    ! === 1. Pure sigma =======================================================
    ! sigma_t = 0 must give A = 0 and B = sigma identically.
    call aeros_vgrid_init(vg, nlev, sigma_t=0.0_wp, p_top=0.0_wp)

    call check(all(vg%A == 0.0_wp), "sigma_t=0 gives A = 0 exactly", nfail)
    call check(all(vg%B == vg%sigma_half), "sigma_t=0 gives B = sigma exactly", nfail)
    call check(vg%B(nlev) == 1.0_wp, "surface half level has B = 1 exactly", nfail)
    call check(vg%sigma_half(0) == 0.0_wp, "model top at sigma = 0", nfail)

    call test_column(vg, "pure sigma", tol, nfail)
    call aeros_vgrid_end(vg)

    ! === 2. Uniform sigma (the Held-Suarez configuration) ====================
    ! stretch_a = 1 must give evenly spaced sigma. Held & Suarez (1994) specify
    ! 20 equally spaced sigma levels, so M1.5 runs exactly this.
    call aeros_vgrid_init(vg, nlev, stretch_a=1.0_wp, sigma_t=0.0_wp, p_top=0.0_wp)
    call test_uniform(vg, tol, nfail)
    call test_column(vg, "uniform sigma", tol, nfail)
    call aeros_vgrid_end(vg)

    ! === 3. Stretched hybrid, ALL DEFAULTS ===================================
    ! Deliberately takes no arguments: this asserts the configuration aeros
    ! actually ships with, not a set of parameters invented for the test.
    call aeros_vgrid_init(vg, nlev)

    ! docs/design.md section 4.1: 3-4 levels below 850 hPa, because the polar
    ! winter boundary layer is 50-200 m deep with 10-25 K inversions and a
    ! single bulk layer physically cannot hold one. This is the design's
    ! highest-value vertical requirement, so it is asserted, not just intended.
    call check(count(vg%p_ref_full > 85000.0_wp) >= 3, &
                "at least 3 full levels below 850 hPa", nfail)
    write(*,"(a40,i6)") "   levels below 850 hPa ", count(vg%p_ref_full > 85000.0_wp)

    ! Hybrid: pure pressure aloft, terrain-following at the surface.
    call check(vg%B(0) == 0.0_wp, "B = 0 at the model top", nfail)
    call check(vg%B(nlev) == 1.0_wp, "B = 1 at the surface", nfail)
    call check(any(vg%B(1:nlev-1) == 0.0_wp), "hybrid has pure-pressure levels aloft", nfail)

    ! At the reference surface pressure the hybrid must reproduce its own sigma
    ! profile exactly -- that is the property that makes sigma_t a free
    ! parameter rather than a change of level placement.
    call check(maxval(abs(vg%p_ref_half - vg%sigma_half*vg%ps_ref)) &
                    <= tol*vg%ps_ref, &
                "hybrid reduces to its sigma profile at ps_ref", nfail)

    call test_column(vg, "stretched hybrid", tol, nfail)
    call aeros_vgrid_print(vg)
    call aeros_vgrid_end(vg)

    ! === 4. Explicit A/B table ===============================================
    call test_table(nfail)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_vertical:: PASS"
    else
        write(*,*) "test_vertical:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine test_column(vg, label, tol, nfail)
        ! Pressures, layer masses and the hydrostatic integral for a column at
        ! a surface pressure DIFFERENT from the reference, so that a hybrid
        ! coordinate is actually exercised rather than collapsing to sigma.

        implicit none

        type(aeros_vgrid_class), intent(in) :: vg
        character(len=*), intent(in) :: label
        real(wp), intent(in) :: tol
        integer, intent(inout) :: nfail

        real(wp), parameter :: ps   = 97000.0_wp   ! a plausible high-altitude surface
        real(wp), parameter :: tiso = 250.0_wp     ! isothermal test atmosphere
        real(wp), parameter :: phis = 3000.0_wp*9.80665_wp

        real(wp) :: p_half(0:vg%nlev), p_full(vg%nlev), dp_lev(vg%nlev)
        real(wp) :: alpha(vg%nlev), temp(vg%nlev)
        real(wp) :: phi_full(vg%nlev), phi_half(0:vg%nlev)
        real(wp) :: total, err, expect
        integer  :: k
        logical  :: ok

        write(*,"(a)") ""
        write(*,"(a,a)") " -- ", trim(label)

        call aeros_vgrid_pressure(vg, ps, p_half, p_full, dp_lev)

        ! Endpoints.
        call check(p_half(vg%nlev) == ps, "p_half(nlev) = ps exactly", nfail)

        ! Mass. sum(dp) telescopes to ps - p_top, and must do so EXACTLY.
        total  = sum(dp_lev)
        expect = ps - p_half(0)
        call check(total == expect, "sum(dp) = ps - p_top exactly", nfail)

        ! Monotonicity and positive layer masses.
        ok = .TRUE.
        do k = 1, vg%nlev
            if (dp_lev(k) <= 0.0_wp) ok = .FALSE.
            if (p_full(k) <= p_half(k-1) .or. p_full(k) >= p_half(k)) ok = .FALSE.
        end do
        call check(ok, "layers have positive mass, p_full inside its layer", nfail)

        ! alpha. Bounded in (0,1) for every layer -- it is a fractional
        ! position within the layer -- and exactly ln 2 for a top layer whose
        ! upper interface is at zero pressure.
        call aeros_vgrid_alpha(vg, p_half, alpha)
        call check(all(alpha > 0.0_wp .and. alpha < 1.0_wp), "alpha in (0,1)", nfail)
        if (p_half(0) <= 0.0_wp) then
            call check(alpha(1) == log(2.0_wp), "alpha(1) = ln 2 at a zero-pressure top", nfail)
        end if

        ! Hydrostatic integral against the analytic isothermal atmosphere.
        !
        ! For T = const the discrete half-level relation telescopes exactly:
        !   Phi_(k-1/2) = Phi_s + R T ln(ps / p_(k-1/2))
        ! so this is an EXACTNESS test of the integration, not an accuracy test
        ! of a discretization. Any drift here is a bug in the recursion.
        temp = tiso
        call aeros_hydrostatic(vg, phis, temp, p_half, phi_full, phi_half)

        err = 0.0_wp
        do k = vg%nlev, 1, -1
            if (p_half(k-1) <= 0.0_wp) cycle
            expect = phis + R_d*tiso*log(ps/p_half(k-1))
            err = max(err, abs(phi_half(k-1) - expect)/max(abs(expect), 1.0_wp))
        end do
        call check(err < tol, "isothermal half-level geopotential is exact", nfail)
        write(*,"(a40,es12.3)") "   hydrostatic relative error ", err

        ! Full-level geopotential must sit inside its own layer.
        ok = .TRUE.
        do k = 1, vg%nlev
            if (phi_full(k) <= phi_half(k)) ok = .FALSE.
            if (phi_half(k-1) /= MV) then
                if (phi_full(k) >= phi_half(k-1)) ok = .FALSE.
            end if
        end do
        call check(ok, "phi_full lies between its bounding half levels", nfail)

        ! Geopotential increases upward, and the surface value is honoured.
        call check(phi_half(vg%nlev) == phis, "phi_half(nlev) = phis exactly", nfail)
        ok = .TRUE.
        do k = vg%nlev, 2, -1
            if (phi_full(k-1) <= phi_full(k)) ok = .FALSE.
        end do
        call check(ok, "geopotential increases upward", nfail)

        return

    end subroutine test_column

    subroutine test_uniform(vg, tol, nfail)

        implicit none

        type(aeros_vgrid_class), intent(in) :: vg
        real(wp), intent(in) :: tol
        integer, intent(inout) :: nfail

        real(wp) :: dsig, err
        integer  :: j

        dsig = 1.0_wp/real(vg%nlev, wp)
        err  = 0.0_wp
        do j = 0, vg%nlev
            err = max(err, abs(vg%sigma_half(j) - real(j, wp)*dsig))
        end do

        call check(err < tol, "stretch_a=1 gives evenly spaced sigma (Held-Suarez)", nfail)

        return

    end subroutine test_uniform

    subroutine test_table(nfail)
        ! An explicit A/B table must be taken verbatim. This is the path a
        ! published level set (ECHAM L19, IFS L31) would arrive by, so it has
        ! to bypass the generator entirely rather than being fitted by it.

        implicit none

        integer, intent(inout) :: nfail

        integer, parameter :: n = 4
        type(aeros_vgrid_class) :: vgt
        real(wp) :: a_tab(0:n), b_tab(0:n)
        real(wp) :: p_half(0:n), p_full(n), dp_lev(n)

        a_tab = [    0.0_wp, 5000.0_wp, 8000.0_wp, 4000.0_wp, 0.0_wp]
        b_tab = [    0.0_wp,    0.0_wp,    0.1_wp,    0.5_wp, 1.0_wp]

        call aeros_vgrid_init(vgt, n, a_half=a_tab, b_half=b_tab)

        call check(all(vgt%A == a_tab), "explicit A table used verbatim", nfail)
        call check(all(vgt%B == b_tab), "explicit B table used verbatim", nfail)

        call aeros_vgrid_pressure(vgt, p0, p_half, p_full, dp_lev)
        call check(sum(dp_lev) == p0 - p_half(0), "table column conserves mass", nfail)
        call check(all(dp_lev > 0.0_wp), "table column has positive layer masses", nfail)

        call aeros_vgrid_end(vgt)

        return

    end subroutine test_table

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

end program test_vertical
