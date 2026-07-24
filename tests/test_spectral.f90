program test_spectral
    ! Acceptance test for the SHTns wrapper.
    !
    ! This is the test that proves the external dependency is wired correctly
    ! and that aeros' conventions (orthonormal, no Condon-Shortley phase,
    ! phi-contiguous, 1-based indexing) are what the code below assumes. Every
    ! later result rests on it, so it runs on every build.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : sp, dp, wp, wp_sh, pi, r_earth
    use aeros_spectral

    implicit none

    type(aeros_sht_class) :: sht

    integer, parameter :: trunc = 21

    real(wp),       allocatable :: field(:,:), field2(:,:)
    complex(wp_sh), allocatable :: coeffs(:), coeffs2(:)

    real(dp) :: tol, err, total, expect
    integer  :: lm, l, m, nfail

    nfail = 0

    ! Round-trip accuracy is limited by the GRID-side kind, since that is
    ! where an sp build converts. Spectral space is double in both builds.
    tol = 20.0_dp*real(epsilon(1.0_wp), dp)

    call aeros_sht_init(sht, trunc, quick=.TRUE.)

    write(*,*) "test_spectral:: T", sht%trunc, " grid ", sht%nlon, "x", sht%nlat, &
                " nlm ", sht%nlm

    allocate(field(sht%nlon,sht%nlat))
    allocate(field2(sht%nlon,sht%nlat))
    allocate(coeffs(sht%nlm))
    allocate(coeffs2(sht%nlm))

    ! === 1. Grid sizing ======================================================
    ! The quadratic rule (nlon >= 3T+1) is what keeps advection unaliased. A
    ! silently smaller grid would still run, and would still look plausible.
    call check(sht%nlon >= 3*trunc + 1, "grid is quadratically unaliased", nfail)
    call check(sht%nlat == sht%nlon/2,  "nlat = nlon/2", nfail)
    call check(sht%nlm == (trunc+1)*(trunc+2)/2, "nlm = (T+1)(T+2)/2", nfail)

    ! === 2. Coefficient indexing =============================================
    ! aeros_sht_lm must agree with the (l,m) table SHTns itself reports. These
    ! are two independent statements of the same layout; if they ever disagree
    ! every spectral operator is silently permuted.
    do m = 0, sht%mmax
        do l = m, sht%lmax
            lm = aeros_sht_lm(sht, l, m)
            if (lm < 1 .or. lm > sht%nlm) then
                nfail = nfail + 1
                write(*,*) "  FAIL: lm out of range for (l,m) = ", l, m, " -> ", lm
                cycle
            end if
            if (sht%l_of_lm(lm) /= l .or. sht%m_of_lm(lm) /= m) then
                nfail = nfail + 1
                write(*,*) "  FAIL: index mismatch (l,m) = ", l, m, &
                            " -> lm ", lm, " reports ", sht%l_of_lm(lm), sht%m_of_lm(lm)
            end if
        end do
    end do
    call check(aeros_sht_lm(sht, trunc+1, 0) == -1, "out-of-truncation index returns -1", nfail)

    ! === 3. Synthesis/analysis round trip ====================================
    ! Start in spectral space so the field is exactly representable at this
    ! truncation: a grid field built by any other means carries components
    ! above T that analysis must discard, which would make the comparison a
    ! statement about aliasing rather than about the transform.
    call fill_test_spectrum(coeffs, sht)

    call aeros_sht_synthesis(sht, coeffs, field)
    field2 = field                       ! analysis may destroy its input
    call aeros_sht_analysis(sht, field2, coeffs2)

    err = maxval(abs(coeffs2 - coeffs))/maxval(abs(coeffs))
    call check(err < tol, "spectral -> grid -> spectral round trip", nfail)
    write(*,"(a40,es12.3,a,es12.3)") "   round-trip relative error ", err, "  tol ", tol

    ! === 4. Quadrature =======================================================
    ! The Gauss weights must integrate a constant exactly. Every conservation
    ! check the design requires "to machine precision" (docs/design.md section
    ! 3.7 risk 2) is this identity with a different integrand.
    field  = 1.0_wp
    total  = aeros_sht_surface_integral(sht, field)
    expect = 4.0_dp*pi*r_earth*r_earth
    err    = abs(total - expect)/expect
    call check(err < tol, "surface integral of 1 = 4*pi*a^2", nfail)
    write(*,"(a40,es12.3)") "   quadrature relative error ", err

    ! A non-constant field with zero mean must integrate to zero. This catches
    ! a weight array that is right in total but mirrored wrongly about the
    ! equator -- which the constant test above cannot see.
    coeffs = (0.0_wp_sh, 0.0_wp_sh)
    coeffs(aeros_sht_lm(sht,1,0)) = (1.0_wp_sh, 0.0_wp_sh)
    call aeros_sht_synthesis(sht, coeffs, field)
    total = aeros_sht_surface_integral(sht, field)
    err   = abs(total)/expect
    call check(err < tol, "surface integral of Y_10 = 0", nfail)

    ! === 5. Laplacian ========================================================
    ! Diagonal in spectral space, eigenvalue -l(l+1)/a^2. Checked against the
    ! degree table rather than a recomputed index, and then round-tripped
    ! through the inverse, which must be the identity on all l > 0.
    call fill_test_spectrum(coeffs, sht)
    coeffs2 = coeffs
    call aeros_sht_laplacian(sht, coeffs2)

    err = 0.0_dp
    do lm = 1, sht%nlm
        l = sht%l_of_lm(lm)
        err = max(err, abs(coeffs2(lm) - coeffs(lm)*(-real(l*(l+1),dp)/(r_earth*r_earth))))
    end do
    call check(err == 0.0_dp, "Laplacian eigenvalue -l(l+1)/a^2 exact", nfail)

    call aeros_sht_laplacian(sht, coeffs2, inverse=.TRUE.)
    coeffs(aeros_sht_lm(sht,0,0)) = (0.0_wp_sh, 0.0_wp_sh)   ! l=0 is not recoverable
    err = maxval(abs(coeffs2 - coeffs))/maxval(abs(coeffs))
    call check(err < 1.0e-12_dp, "inverse Laplacian recovers l > 0", nfail)

    ! === 6. Vector transforms ================================================
    ! Only that the spheroidal/toroidal pair round-trips. The mapping to
    ! vorticity, divergence and (u,v) belongs to the dynamical core (M1) and is
    ! deliberately not asserted here.
    call test_vector_roundtrip(sht, tol, nfail)

    call aeros_sht_end(sht)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_spectral:: PASS"
    else
        write(*,*) "test_spectral:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine fill_test_spectrum(coeffs, sht)
        ! A reproducible, smooth-ish test spectrum: amplitude decaying with
        ! degree, so the field is dominated by resolved scales but every
        ! coefficient is non-zero.

        implicit none

        complex(wp_sh), intent(out) :: coeffs(:)
        type(aeros_sht_class), intent(in) :: sht

        integer  :: lm, l, m
        real(dp) :: amp

        do lm = 1, sht%nlm
            l = sht%l_of_lm(lm)
            m = sht%m_of_lm(lm)
            amp = 1.0_dp/real((l+1)*(l+1), dp)
            coeffs(lm) = cmplx(amp*cos(real(3*l + m, dp)), &
                                amp*sin(real(l + 2*m, dp)), wp_sh)
            ! m = 0 coefficients of a real field are real.
            if (m == 0) coeffs(lm) = cmplx(real(coeffs(lm)), 0.0_wp_sh, wp_sh)
        end do

        return

    end subroutine fill_test_spectrum

    subroutine test_vector_roundtrip(sht, tol, nfail)

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(dp), intent(in) :: tol
        integer, intent(inout) :: nfail

        real(wp),       allocatable :: vth(:,:), vph(:,:)
        complex(wp_sh), allocatable :: slm(:), tlm(:), slm2(:), tlm2(:)
        real(dp) :: err

        allocate(vth(sht%nlon,sht%nlat), vph(sht%nlon,sht%nlat))
        allocate(slm(sht%nlm), tlm(sht%nlm), slm2(sht%nlm), tlm2(sht%nlm))

        call fill_test_spectrum(slm, sht)
        call fill_test_spectrum(tlm, sht)

        ! l = 0 carries no horizontal vector field; SHTns returns zero there,
        ! so seeding it would guarantee a spurious mismatch.
        slm(aeros_sht_lm(sht,0,0)) = (0.0_wp_sh, 0.0_wp_sh)
        tlm(aeros_sht_lm(sht,0,0)) = (0.0_wp_sh, 0.0_wp_sh)

        call aeros_sht_synthesis_vec(sht, slm, tlm, vth, vph)
        call aeros_sht_analysis_vec(sht, vth, vph, slm2, tlm2)

        err = max(maxval(abs(slm2 - slm)), maxval(abs(tlm2 - tlm))) &
                / max(maxval(abs(slm)), maxval(abs(tlm)))
        call check(err < tol, "spheroidal/toroidal round trip", nfail)
        write(*,"(a40,es12.3)") "   vector round-trip relative error ", err

        return

    end subroutine test_vector_roundtrip

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

end program test_spectral
