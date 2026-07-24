program test_vordiv
    ! Acceptance test for the vorticity/divergence <-> (u,v) mapping.
    !
    ! This mapping has no free parameters, so it is either exactly right or
    ! wrong -- there is nothing to tune and no tolerance to negotiate. Every
    ! check below is therefore against a CLOSED-FORM answer, not against
    ! another run of the same code:
    !
    !   solid-body rotation   u = U0 cos(lat), v = 0 has zeta = (2 U0/a) sin(lat)
    !                         and D = 0 identically. This is the check that
    !                         catches a sign error, a factor of a, or a swapped
    !                         spheroidal/toroidal channel -- all of which a
    !                         round-trip test would happily pass.
    !
    !   pure divergence       the mirror case, built from a divergent potential,
    !                         so that a swapped channel cannot hide.
    !
    !   round trip            identity on l >= 1.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, pi, r_earth, radians_to_degrees
    use aeros_spectral
    use aeros_vordiv

    implicit none

    type(aeros_sht_class) :: sht

    integer, parameter :: trunc = 21

    real(wp) :: tol
    integer  :: nfail

    nfail = 0
    tol   = 1.0e3_wp*epsilon(1.0_wp)   ! transforms are ~10 ulp; this is not tight-fitting

    call aeros_sht_init(sht, trunc, quick=.TRUE.)

    write(*,"(a,i0,a,i0,a,i0)") " test_vordiv:: T", sht%trunc, "  grid ", sht%nlon, "x", sht%nlat

    call test_solid_body(sht, tol, nfail)
    call test_pure_divergence(sht, tol, nfail)
    call test_roundtrip(sht, tol, nfail)
    call test_l0_discarded(sht, nfail)

    call aeros_sht_end(sht)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_vordiv:: PASS"
    else
        write(*,*) "test_vordiv:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine test_solid_body(sht, tol, nfail)
        ! Solid-body rotation: u = U0 cos(lat), v = 0.
        !
        ! Analytically zeta = (2 U0/a) sin(lat) and D = 0. In spectral space
        ! the vorticity is a single coefficient, (l,m) = (1,0), with the value
        ! (2 U0/a) / N where N = sqrt(3/4pi) is the orthonormal Y_10
        ! normalization -- so this pins the normalization too, not just the
        ! signs.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(wp), intent(in) :: tol
        integer, intent(inout) :: nfail

        real(wp), parameter :: u0 = 40.0_wp    ! [m s-1]

        real(wp)       :: u(sht%nlon,sht%nlat), v(sht%nlon,sht%nlat)
        complex(wp_sh) :: vor(sht%nlm), div(sht%nlm)
        real(wp)       :: zeta_g(sht%nlon,sht%nlat)
        real(dp)       :: expect, err, scale, nrm
        integer        :: i, j, lm10

        write(*,*) ""
        write(*,*) " -- solid-body rotation"

        ! sinlat is cos(colatitude); cos(lat) = sin(colatitude).
        do j = 1, sht%nlat
            do i = 1, sht%nlon
                u(i,j) = u0*real(sin(sht%colat(j)), wp)
                v(i,j) = 0.0_wp
            end do
        end do

        call aeros_vordiv_from_uv(sht, u, v, vor, div)

        ! Divergence must vanish identically.
        scale = maxval(abs(vor))
        err   = maxval(abs(div))/scale
        call check(err < tol, "solid-body rotation is non-divergent", nfail)
        write(*,"(a44,es12.3)") "   |D|/|zeta| ", err

        ! Vorticity must be the single (1,0) coefficient, with the right value.
        lm10   = aeros_sht_lm(sht, 1, 0)
        nrm    = sqrt(3.0_dp/(4.0_dp*pi))
        expect = (2.0_dp*real(u0, dp)/r_earth)/nrm

        err = abs(real(vor(lm10), dp) - expect)/abs(expect)
        call check(err < tol, "zeta_(1,0) has the analytic value", nfail)
        write(*,"(a44,es12.3,a,es12.3)") "   zeta_10 ", real(vor(lm10), dp), &
                                            "  expected ", expect
        write(*,"(a44,es12.3)") "   relative error ", err

        ! Every other coefficient must be zero.
        vor(lm10) = (0.0_wp_sh, 0.0_wp_sh)
        err = maxval(abs(vor))/abs(expect)
        call check(err < tol, "no other vorticity coefficient is excited", nfail)

        ! And in grid space, zeta = (2 U0/a) sin(lat).
        call aeros_vordiv_from_uv(sht, u, v, vor, div)
        call aeros_sht_synthesis(sht, vor, zeta_g)

        err   = 0.0_dp
        scale = 2.0_dp*real(u0, dp)/r_earth
        do j = 1, sht%nlat
            do i = 1, sht%nlon
                err = max(err, abs(real(zeta_g(i,j), dp) - scale*sht%sinlat(j)))
            end do
        end do
        err = err/scale
        call check(err < tol, "grid-space zeta = (2 U0/a) sin(lat)", nfail)
        write(*,"(a44,es12.3)") "   grid-space relative error ", err

        ! Sign, stated as a physical fact rather than a formula: westerly flow
        ! has CYCLONIC (positive) vorticity in the northern hemisphere. A
        ! flipped toroidal sign passes every magnitude check above and fails
        ! only here.
        call check(real(zeta_g(1,1), dp) > 0.0_dp, &
                    "westerly flow gives positive vorticity in the north", nfail)
        call check(real(zeta_g(1,sht%nlat), dp) < 0.0_dp, &
                    "westerly flow gives negative vorticity in the south", nfail)

        return

    end subroutine test_solid_body

    subroutine test_pure_divergence(sht, tol, nfail)
        ! The mirror of the rotational case: a wind field built from divergence
        ! alone must come back with zero vorticity.
        !
        ! Built by synthesis from a single divergence coefficient rather than
        ! by writing down a grid field, because the point is to check that the
        ! two channels do not leak into one another -- and that is exactly what
        ! a hand-written grid field with a small rotational contamination would
        ! obscure.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(wp), intent(in) :: tol
        integer, intent(inout) :: nfail

        real(wp)       :: u(sht%nlon,sht%nlat), v(sht%nlon,sht%nlat)
        complex(wp_sh) :: vor(sht%nlm), div(sht%nlm)
        complex(wp_sh) :: vor2(sht%nlm), div2(sht%nlm)
        real(dp)       :: err
        integer        :: lm21

        write(*,*) ""
        write(*,*) " -- pure divergence"

        vor = (0.0_wp_sh, 0.0_wp_sh)
        div = (0.0_wp_sh, 0.0_wp_sh)
        lm21 = aeros_sht_lm(sht, 2, 1)
        div(lm21) = cmplx(1.0e-6_dp, -4.0e-7_dp, wp_sh)

        call aeros_uv_from_vordiv(sht, vor, div, u, v)
        call aeros_vordiv_from_uv(sht, u, v, vor2, div2)

        err = maxval(abs(vor2))/maxval(abs(div))
        call check(err < tol, "a divergent field carries no vorticity", nfail)
        write(*,"(a44,es12.3)") "   |zeta|/|D| ", err

        err = maxval(abs(div2 - div))/maxval(abs(div))
        call check(err < tol, "divergence survives the round trip", nfail)
        write(*,"(a44,es12.3)") "   divergence relative error ", err

        return

    end subroutine test_pure_divergence

    subroutine test_roundtrip(sht, tol, nfail)
        ! (zeta,D) -> (u,v) -> (zeta,D) on a full spectrum, identity for l >= 1.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(wp), intent(in) :: tol
        integer, intent(inout) :: nfail

        real(wp)       :: u(sht%nlon,sht%nlat), v(sht%nlon,sht%nlat)
        complex(wp_sh) :: vor(sht%nlm), div(sht%nlm)
        complex(wp_sh) :: vor2(sht%nlm), div2(sht%nlm)
        real(dp)       :: amp, err, scale
        integer        :: lm, l, m

        write(*,*) ""
        write(*,*) " -- round trip"

        ! A spectrum with realistic amplitudes: ~1e-5 s-1 vorticity at the
        ! largest scales, decaying with degree.
        do lm = 1, sht%nlm
            l = sht%l_of_lm(lm)
            m = sht%m_of_lm(lm)
            if (l == 0) then
                vor(lm) = (0.0_wp_sh, 0.0_wp_sh)
                div(lm) = (0.0_wp_sh, 0.0_wp_sh)
                cycle
            end if
            amp = 1.0e-5_dp/real(l*l, dp)
            vor(lm) = cmplx(amp*cos(real(3*l+m,dp)), amp*sin(real(l+2*m,dp)), wp_sh)
            div(lm) = cmplx(0.3_dp*amp*sin(real(l+m,dp)), 0.3_dp*amp*cos(real(2*l+m,dp)), wp_sh)
            if (m == 0) then
                vor(lm) = cmplx(real(vor(lm)), 0.0_wp_sh, wp_sh)
                div(lm) = cmplx(real(div(lm)), 0.0_wp_sh, wp_sh)
            end if
        end do

        call aeros_uv_from_vordiv(sht, vor, div, u, v)
        call aeros_vordiv_from_uv(sht, u, v, vor2, div2)

        scale = max(maxval(abs(vor)), maxval(abs(div)))
        err   = max(maxval(abs(vor2 - vor)), maxval(abs(div2 - div)))/scale
        call check(err < tol, "(zeta,D) -> (u,v) -> (zeta,D) is the identity", nfail)
        write(*,"(a44,es12.3,a,es12.3)") "   relative error ", err, "  tol ", real(tol, dp)

        ! Winds should be physically plausible, not merely self-consistent: a
        ! 1e-5 s-1 vorticity at planetary scale is a few tens of m/s.
        call check(maxval(abs(u)) > 1.0_wp .and. maxval(abs(u)) < 200.0_wp, &
                    "wind speeds are physically plausible", nfail)
        write(*,"(a44,f9.2,a)") "   max |u| ", maxval(abs(u)), " m s-1"

        return

    end subroutine test_roundtrip

    subroutine test_l0_discarded(sht, nfail)
        ! l = 0 carries no horizontal vector field, so analysis must return
        ! exactly zero there and synthesis must ignore whatever it is given.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        integer, intent(inout) :: nfail

        real(wp)       :: u(sht%nlon,sht%nlat), v(sht%nlon,sht%nlat)
        real(wp)       :: u2(sht%nlon,sht%nlat), v2(sht%nlon,sht%nlat)
        complex(wp_sh) :: vor(sht%nlm), div(sht%nlm)
        complex(wp_sh) :: vor_rt(sht%nlm), div_rt(sht%nlm)
        integer        :: lm00

        write(*,*) ""
        write(*,*) " -- degree 0"

        lm00 = aeros_sht_lm(sht, 0, 0)

        vor = (0.0_wp_sh, 0.0_wp_sh)
        div = (0.0_wp_sh, 0.0_wp_sh)
        vor(aeros_sht_lm(sht,3,2)) = cmplx(1.0e-5_dp, 2.0e-6_dp, wp_sh)

        ! Analysis must report exactly zero at l = 0. Uses its own copies so
        ! that vor/div are left untouched for the synthesis check below -- the
        ! two must be driven from the SAME spectrum, or the comparison is
        ! between two different fields rather than between two l=0 choices.
        call aeros_uv_from_vordiv(sht, vor, div, u, v)
        call aeros_vordiv_from_uv(sht, u, v, vor_rt, div_rt)

        call check(vor_rt(lm00) == (0.0_wp_sh, 0.0_wp_sh), "analysis returns zeta_00 = 0", nfail)
        call check(div_rt(lm00) == (0.0_wp_sh, 0.0_wp_sh), "analysis returns D_00 = 0", nfail)

        ! Synthesis must ignore a nonzero l = 0 input entirely.
        vor(lm00) = cmplx(1.0e-3_dp, 0.0_wp_sh, wp_sh)
        div(lm00) = cmplx(1.0e-3_dp, 0.0_wp_sh, wp_sh)
        call aeros_uv_from_vordiv(sht, vor, div, u2, v2)

        call check(all(u2 == u) .and. all(v2 == v), &
                    "a nonzero degree-0 input is discarded, bit-for-bit", nfail)

        return

    end subroutine test_l0_discarded

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

end program test_vordiv
