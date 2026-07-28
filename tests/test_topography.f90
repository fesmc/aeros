program test_topography
    ! Acceptance test for the boundary-field regridder (aeros_bcinput) and the
    ! topography scaling/ramp (aeros_topography).
    !
    ! The regridder is the reusable foundation for every prescribed lower
    ! boundary (topography today; land-sea mask, albedo, SST later), so its two
    ! defining properties are pinned here directly on the bilinear kernel, with
    ! no netCDF file required:
    !
    !   (a) exactness on a constant  -- a flat source field must reproduce
    !       everywhere to round-off, or a "zero over ocean" field would leak a
    !       spurious signal;
    !   (b) accuracy on a smooth analytic lon/lat field -- the interpolant must
    !       track a known function within the bilinear truncation error, and hit
    !       source nodes exactly.
    !
    ! The topography feed's third property -- that l_topography = .false. gives
    ! EXACTLY zero surface geopotential (bit-for-bit aquaplanet) -- is pinned on
    ! aeros_topography_scale, together with the restart-safe time ramp.
    !
    ! Exits non-zero on failure.

    use aeros_defs,       only : dp, wp, pi
    use aeros_bcinput,    only : aeros_bcinput_bilinear
    use aeros_topography, only : aeros_topography_scale

    implicit none

    integer  :: nfail, i, j
    real(wp) :: tol, val, expect, err, maxerr
    real(wp), allocatable :: lon_src(:), lat_src(:), f(:,:)
    integer  :: nx, ny
    real(wp) :: tlon, tlat, day

    nfail  = 0
    tol    = 1.0e-5_wp

    ! === Build a synthetic regular lon/lat source grid =====================
    ! Longitude 0..360 (periodic, last node at 355 so the wrap gap to 360 is
    ! real); latitude 87.5 -> -87.5 stored NORTH-TO-SOUTH like ERA5, to exercise
    ! the descending-latitude branch.
    nx = 72; ny = 36
    allocate(lon_src(nx), lat_src(ny), f(nx,ny))
    do i = 1, nx
        lon_src(i) = 360.0_wp*real(i-1, wp)/real(nx, wp)      ! 0,5,...,355
    end do
    do j = 1, ny
        lat_src(j) = 90.0_wp - 180.0_wp*(real(j, wp) - 0.5_wp)/real(ny, wp)  ! 87.5 .. -87.5
    end do

    ! === (a) Exactness on a constant field =================================
    f = 42.0_wp
    maxerr = 0.0_wp
    do j = 1, 40
        tlon = 360.0_wp*real(j-1, wp)/40.0_wp + 1.234_wp   ! includes the wrap gap > 355
        do i = 1, 40
            tlat = 89.0_wp - 178.0_wp*real(i-1, wp)/39.0_wp
            val  = aeros_bcinput_bilinear(lon_src, lat_src, f, tlon, tlat)
            maxerr = max(maxerr, abs(val - 42.0_wp))
        end do
    end do
    ! Round-off scales with the field value (42), not with epsilon(1): the
    ! bilinear weights sum to 1 exactly, so any error is pure rounding.
    call check(maxerr < 42.0_wp*100.0_wp*epsilon(1.0_wp), &
               "constant field reproduced exactly", nfail)
    write(*,"(a40,es12.3)") "   max |interp - const| ", maxerr

    ! === (b) Accuracy on a smooth analytic field ===========================
    ! f(lon,lat) = analytic_field(lon,lat): smooth, periodic in longitude, so the
    ! bilinear interpolant should track it to O(h^2) and hit nodes exactly.
    do j = 1, ny
        do i = 1, nx
            f(i,j) = analytic_field(lon_src(i), lat_src(j))
        end do
    end do

    ! b1: exact at a source node (bilinear reproduces grid values).
    val    = aeros_bcinput_bilinear(lon_src, lat_src, f, lon_src(10), lat_src(7))
    expect = analytic_field(lon_src(10), lat_src(7))
    call check(abs(val - expect) < tol, "bilinear hits source nodes exactly", nfail)

    ! b2: interior accuracy over a scatter of off-node points, including the
    ! longitude wrap gap and both latitude edges (clamped, so looser there).
    maxerr = 0.0_wp
    do j = 1, 30
        tlon = 360.0_wp*real(j-1, wp)/30.0_wp + 2.5_wp
        do i = 1, 30
            tlat = 80.0_wp - 160.0_wp*real(i-1, wp)/29.0_wp   ! within [-87.5,87.5]
            val    = aeros_bcinput_bilinear(lon_src, lat_src, f, tlon, tlat)
            expect = analytic_field(tlon, tlat)
            err    = abs(val - expect)
            maxerr = max(maxerr, err)
        end do
    end do
    ! Bilinear truncation on a 5-deg source grid for this smooth field.
    call check(maxerr < 5.0e-3_wp, "analytic field within bilinear tolerance", nfail)
    write(*,"(a40,es12.3)") "   max |interp - analytic| ", maxerr

    ! === (c) l_topography = .false. => exactly zero ========================
    ! Every combination of time and ramp must give bit-zero when topography is
    ! off -- this is what preserves the flat-aquaplanet, bit-for-bit default.
    call check(aeros_topography_scale(.false., 0.0_wp,      20.0_wp) == 0.0_wp, &
               "topo off: scale = 0 at t=0",              nfail)
    call check(aeros_topography_scale(.false., 1.0e7_wp,    20.0_wp) == 0.0_wp, &
               "topo off: scale = 0 at large t",          nfail)
    call check(aeros_topography_scale(.false., 5.0e5_wp,     0.0_wp) == 0.0_wp, &
               "topo off: scale = 0 with zero ramp",      nfail)

    ! === Restart-safe time ramp (pure function of elapsed time) ============
    day = 86400.0_wp
    call check(aeros_topography_scale(.true., 0.0_wp,       20.0_wp) == 0.0_wp, &
               "ramp: 0 at t=0",                          nfail)
    call check(abs(aeros_topography_scale(.true., 10.0_wp*day, 20.0_wp) - 0.5_wp) < tol, &
               "ramp: 0.5 at half ramp",                  nfail)
    call check(aeros_topography_scale(.true., 20.0_wp*day, 20.0_wp) == 1.0_wp, &
               "ramp: 1 at end of ramp",                  nfail)
    call check(aeros_topography_scale(.true., 40.0_wp*day, 20.0_wp) == 1.0_wp, &
               "ramp: held at 1 past ramp",               nfail)
    call check(aeros_topography_scale(.true., 0.0_wp,       0.0_wp) == 1.0_wp, &
               "ramp: full topography from t=0 when ramp_days=0", nfail)

    deallocate(lon_src, lat_src, f)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_topography:: PASS"
    else
        write(*,*) "test_topography:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    real(wp) function analytic_field(lon, lat) result(v)
        ! A smooth, longitude-periodic test field.
        implicit none
        real(wp), intent(in) :: lon, lat
        real(wp) :: rlon, rlat
        rlon = lon*pi/180.0_wp
        rlat = lat*pi/180.0_wp
        v = 1.0_wp + 0.3_wp*sin(rlon) + 0.2_wp*cos(2.0_wp*rlon) &
                   + 0.5_wp*sin(rlat)
        return
    end function analytic_field

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

end program test_topography
