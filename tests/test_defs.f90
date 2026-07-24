program test_defs
    ! Acceptance test for aeros_defs: precision policy, constants, and the
    ! namelist round trip.
    !
    ! The precision checks are not decoration. docs/design.md section 4 calls
    ! for a single-precision core, while SHTns' interface is double; the whole
    ! design rests on those two facts being reconciled deliberately rather than
    ! by accident. If wp_sh ever stops being C double, the transforms break
    ! silently, so it is asserted here as well as at compile time.
    !
    ! Exits non-zero on failure.

    use, intrinsic :: iso_c_binding, only : c_double
    use aeros_defs

    implicit none

    type(aeros_param_class) :: par

    integer :: nfail

    nfail = 0

    ! === Precision policy ====================================================
    call check(wp == sp .or. wp == dp, "wp is sp or dp", nfail)
    call check(wp_sh == dp,            "wp_sh is dp",    nfail)
    call check(wp_sh == c_double,      "wp_sh is C double (SHTns interface)", nfail)

    write(*,"(a40,i4,a)") "   wp    ", storage_size(1.0_wp)/8,    " bytes"
    write(*,"(a40,i4,a)") "   wp_sh ", storage_size(1.0_wp_sh)/8, " bytes"

    ! === Constants ===========================================================
    ! Sanity, not accuracy: these catch a mistyped exponent, which is the error
    ! that survives review and then quietly rescales a whole model.
    call check(abs(kappa - R_d/cp_d) < 1.0e-12_dp, "kappa = R_d/cp_d", nfail)
    call check(abs(kappa - 0.2857_dp) < 1.0e-3_dp, "kappa ~ 0.2857",   nfail)
    call check(abs(L_s - (L_v + L_f)) < 1.0e-6_dp, "L_s = L_v + L_f",  nfail)
    call check(r_earth > 6.3e6_dp .and. r_earth < 6.4e6_dp, "r_earth ~ 6371 km", nfail)
    call check(omega > 7.0e-5_dp .and. omega < 7.5e-5_dp,   "omega ~ 7.29e-5 s-1", nfail)
    call check(abs(pi - acos(-1.0_dp)) < 1.0e-15_dp, "pi", nfail)
    call check(abs(degrees_to_radians*radians_to_degrees - 1.0_dp) < 1.0e-15_dp, &
                "degree/radian conversions are inverse", nfail)

    ! === Namelist ============================================================
    call aeros_par_load(par, "par/aeros.nml", "aeros")

    call check(par%trunc > 0,  "trunc read", nfail)
    call check(par%nlev  > 0,  "nlev read",  nfail)
    call check(par%dt    > 0.0_wp, "dt read", nfail)

    write(*,"(a40,i6)")   "   trunc ", par%trunc
    write(*,"(a40,i6)")   "   nlev  ", par%nlev
    write(*,"(a40,f8.1)") "   dt    ", par%dt

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_defs:: PASS"
    else
        write(*,*) "test_defs:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

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

end program test_defs
