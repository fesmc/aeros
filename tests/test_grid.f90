program test_grid
    ! Acceptance test for the Gaussian grid.
    !
    ! The grid is derived from the transform's own quadrature, so these checks
    ! are really about that derivation staying faithful: total area, latitude
    ! ordering, hemispheric symmetry and the Coriolis parameter.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : sp, dp, wp, pi, r_earth, omega, aeros_grid_class
    use aeros_spectral, only : aeros_sht_class, aeros_sht_init, aeros_sht_end, &
                                aeros_sht_grid_size
    use aeros_grid,     only : aeros_grid_init, aeros_grid_end

    implicit none

    type(aeros_sht_class)  :: sht
    type(aeros_grid_class) :: grd

    real(dp) :: total, expect, err, tol
    integer  :: nlon, nlat, j, nfail

    nfail = 0
    tol   = 20.0_dp*real(epsilon(1.0_wp), dp)

    ! === Grid sizing rule, independent of SHTns ==============================
    ! The truncations named in docs/design.md section 3.1 must map to the
    ! standard grids the literature uses, or every comparison against
    ! Loefverstroem & Liakka is against a different model.
    call aeros_sht_grid_size(21, nlon, nlat)
    call check(nlon ==  64 .and. nlat ==  32, "T21 -> 64x32",   nfail)
    call aeros_sht_grid_size(31, nlon, nlat)
    call check(nlon ==  96 .and. nlat ==  48, "T31 -> 96x48",   nfail)
    call aeros_sht_grid_size(42, nlon, nlat)
    call check(nlon == 128 .and. nlat ==  64, "T42 -> 128x64",  nfail)
    call aeros_sht_grid_size(85, nlon, nlat)
    call check(nlon == 256 .and. nlat == 128, "T85 -> 256x128", nfail)

    ! === Realized grid =======================================================
    call aeros_sht_init(sht, 31, quick=.TRUE.)
    call aeros_grid_init(grd, sht)

    call check(grd%ncol == grd%nlon*grd%nlat, "ncol = nlon*nlat", nfail)

    ! Total area. Uses the cell areas the model itself will weight every
    ! conservation diagnostic with, accumulated in dp.
    total = 0.0_dp
    do j = 1, grd%nlat
        total = total + real(grd%nlon, dp)*real(grd%area(1,j), dp)
    end do
    expect = 4.0_dp*pi*r_earth*r_earth
    err    = abs(total - expect)/expect
    call check(err < tol, "sum(area) = 4*pi*a^2", nfail)
    write(*,"(a40,es12.3)") "   area relative error ", err

    ! Latitude ordering: north to south, no exact poles (Gauss nodes are
    ! interior). A reversed array here would transpose every field written.
    call check(grd%lat(1) > grd%lat(grd%nlat), "latitude runs north to south", nfail)
    call check(grd%lat(1) < 90.0_wp .and. grd%lat(grd%nlat) > -90.0_wp, &
                "Gauss latitudes exclude the poles", nfail)

    ! Hemispheric symmetry of the nodes and weights. Latitudes are compared
    ! relative to their own scale (180 degrees), not absolutely -- an absolute
    ! tolerance here is a units error that happens to pass in single precision.
    err = 0.0_dp
    do j = 1, grd%nlat
        err = max(err, abs(real(grd%lat(j), dp) + real(grd%lat(grd%nlat-j+1), dp))/180.0_dp)
    end do
    call check(err < tol, "latitudes are antisymmetric about the equator", nfail)
    write(*,"(a40,es12.3)") "   latitude asymmetry ", err

    err = 0.0_dp
    do j = 1, grd%nlat
        err = max(err, abs(grd%gauss_w(j) - grd%gauss_w(grd%nlat-j+1)))
    end do
    call check(err == 0.0_dp, "Gauss weights are exactly symmetric", nfail)

    ! Coriolis: positive in the north, negative in the south, magnitude 2*omega
    ! at the poles.
    call check(grd%coriolis(1,1) > 0.0_wp, "f > 0 in the northern hemisphere", nfail)
    call check(grd%coriolis(1,grd%nlat) < 0.0_wp, "f < 0 in the southern hemisphere", nfail)
    call check(abs(real(grd%coriolis(1,1), dp)) < 2.0_dp*omega, "|f| < 2*omega", nfail)

    call aeros_grid_end(grd)
    call aeros_sht_end(sht)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_grid:: PASS"
    else
        write(*,*) "test_grid:: FAIL, ", nfail, " check(s) failed"
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

end program test_grid
