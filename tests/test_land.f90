program test_land
    ! Acceptance test for the land surface (aeros_land): the mask threshold, the
    ! bucket water budget, the slab soil-temperature response, the skin/beta
    ! composition, and the disabled (bit-for-bit) path.
    !
    ! The soil physics is pinned directly on aeros_land_step / _pre with a small
    ! hand-built land state, no netCDF file required. The mask read is pinned
    ! end-to-end through the bcinput regridder against the ERA5 land-sea mask
    ! when that file is present (skipped, not failed, when it is absent).
    !
    ! Exits non-zero on failure.

    use aeros_defs,    only : dp, wp, sigma_sb, T0
    use aeros_bcinput, only : aeros_bcinput_read_field
    use aeros_land,    only : aeros_land_class, aeros_land_is_land, &
                              aeros_land_pre, aeros_land_step, &
                              aeros_land_couple_radiation

    implicit none

    ! Fresh-water density used by the bucket (kept in step with aeros_land's
    ! private RHO_W); the water-budget arithmetic below depends on it.
    real(wp), parameter :: RHO_W = 1000.0_wp

    integer  :: nfail
    real(wp) :: tol
    type(aeros_land_class) :: land
    real(wp) :: sst(2,1)
    real(wp), allocatable :: alb_map(:,:)
    real(wp) :: dt, w0, precip, evap, dw_exp, runoff_depth
    real(wp) :: sw, lw, shf, lhf, f_exp, t0soil

    nfail = 0
    tol   = 1.0e-9_wp

    ! === (1) mask threshold predicate ======================================
    call check(aeros_land_is_land(1.0_wp,  0.5_wp),      "mask: lsm=1 is land",   nfail)
    call check(aeros_land_is_land(0.5_wp,  0.5_wp),      "mask: lsm=threshold is land", nfail)
    call check(.not. aeros_land_is_land(0.49_wp, 0.5_wp),"mask: lsm<threshold is ocean", nfail)
    call check(.not. aeros_land_is_land(0.0_wp,  0.5_wp),"mask: lsm=0 is ocean",   nfail)

    ! === Build a tiny 2x1 land state (cell 1 land, cell 2 ocean) ===========
    call make_land(land)
    dt = 1800.0_wp

    ! === (2a) bucket water budget closes -- no runoff ======================
    ! Small net input keeps the bucket below capacity: Dw = (precip-evap) dt/rho.
    land%w(1,1) = 0.05_wp
    w0     = land%w(1,1)
    precip = 5.0e-5_wp        ! kg m-2 s-1
    evap   = 1.0e-5_wp
    call aeros_land_step(land, flx(2,1,0.0_wp), flx(2,1,0.0_wp), &
                         flx(2,1,0.0_wp), flx(2,1,0.0_wp), &
                         evap2(evap), prc2(precip), dt)
    dw_exp       = (precip - evap)*dt/RHO_W
    runoff_depth = land%runoff(1,1)*dt/RHO_W
    call check(abs((land%w(1,1) - w0) + runoff_depth - dw_exp) < tol, &
               "bucket budget closes (Dw + runoff = precip - evap), no runoff", nfail)
    call check(land%runoff(1,1) == 0.0_wp, "no runoff below capacity", nfail)
    write(*,"(a40,es12.3)") "   |budget residual| ", &
        abs((land%w(1,1) - w0) + runoff_depth - dw_exp)

    ! === (2b) bucket water budget closes -- with runoff ====================
    ! Start near capacity and pour water in: excess must run off, budget still closes.
    land%w(1,1)      = land%w_field_capacity - 1.0e-4_wp
    land%runoff(1,1) = 0.0_wp
    w0     = land%w(1,1)
    precip = 1.0e-3_wp
    evap   = 0.0_wp
    call aeros_land_step(land, flx(2,1,0.0_wp), flx(2,1,0.0_wp), &
                         flx(2,1,0.0_wp), flx(2,1,0.0_wp), &
                         evap2(evap), prc2(precip), dt)
    dw_exp       = (precip - evap)*dt/RHO_W
    runoff_depth = land%runoff(1,1)*dt/RHO_W
    call check(land%w(1,1) == land%w_field_capacity, "runoff caps w at field capacity", nfail)
    call check(land%runoff(1,1) > 0.0_wp, "runoff positive above capacity", nfail)
    call check(abs((land%w(1,1) - w0) + runoff_depth - dw_exp) < tol, &
               "bucket budget closes (Dw + runoff = precip - evap), with runoff", nfail)

    ! === (3) slab soil temperature responds to the net flux ================
    land%t_soil(1,1) = 290.0_wp
    t0soil = land%t_soil(1,1)
    sw = 300.0_wp; lw = 350.0_wp; shf = 20.0_wp; lhf = 60.0_wp
    call aeros_land_step(land, flx(2,1,sw), flx(2,1,lw), &
                         flx(2,1,shf), flx(2,1,lhf), &
                         evap2(0.0_wp), prc2(0.0_wp), dt)
    f_exp = sw + lw - sigma_sb*t0soil**4 - shf - lhf
    call check(abs((land%t_soil(1,1) - t0soil) - f_exp*dt/land%c_soil) < 1.0e-6_wp, &
               "soil temperature responds to net flux (c_soil dT = F dt)", nfail)
    call check(abs(land%fnet(1,1) - f_exp) < 1.0e-6_wp, "fnet diagnostic equals net flux", nfail)

    ! === (4) ocean cell is untouched by the land step ======================
    land%t_soil(2,1) = 111.0_wp
    land%w(2,1)      = 0.222_wp
    call aeros_land_step(land, flx(2,1,sw), flx(2,1,lw), &
                         flx(2,1,shf), flx(2,1,lhf), &
                         evap2(1.0e-4_wp), prc2(1.0e-3_wp), dt)
    call check(land%t_soil(2,1) == 111.0_wp .and. land%w(2,1) == 0.222_wp, &
               "ocean cell untouched by land step", nfail)

    ! === (5) skin/beta composition (aeros_land_pre) ========================
    land%t_soil(1,1) = 295.0_wp
    land%w(1,1)      = 0.5_wp*land%w_crit      ! half the beta knee -> beta = 0.5
    sst(1,1) = 280.0_wp; sst(2,1) = 301.0_wp
    call aeros_land_pre(land, sst)
    call check(land%skin(1,1) == 295.0_wp, "skin = soil temperature on land", nfail)
    call check(land%skin(2,1) == 301.0_wp, "skin = SST over ocean", nfail)
    call check(abs(land%beta(1,1) - 0.5_wp) < tol, "beta = w/w_crit on unsaturated land", nfail)
    call check(land%beta(2,1) == 1.0_wp, "beta = 1 over ocean", nfail)
    ! saturated soil -> beta capped at 1
    land%w(1,1) = 2.0_wp*land%w_crit
    call aeros_land_pre(land, sst)
    call check(land%beta(1,1) == 1.0_wp, "beta capped at 1 on saturated land", nfail)

    ! === (6) radiation albedo map composition ==============================
    land%albedo(1,1) = 0.35_wp
    call aeros_land_couple_radiation(land, alb_map, 0.06_wp)
    call check(allocated(alb_map), "alb_map allocated by couple_radiation", nfail)
    call check(abs(alb_map(1,1) - 0.35_wp) < tol, "alb_map = land albedo on land", nfail)
    call check(abs(alb_map(2,1) - 0.06_wp) < tol, "alb_map = ocean albedo over ocean", nfail)

    ! === (7) disabled land is inert (bit-for-bit ocean-only) ===============
    call check_disabled(nfail)

    ! === (8) mask read end-to-end through bcinput (file-gated) =============
    call check_mask_from_era5(nfail)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_land:: PASS"
    else
        write(*,*) "test_land:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine make_land(land)
        ! A 2x1 grid: cell (1,1) land, cell (2,1) ocean, all arrays allocated.
        type(aeros_land_class), intent(inout) :: land
        land%enabled          = .TRUE.
        land%nlon             = 2
        land%nlat             = 1
        land%w_field_capacity = 0.15_wp
        land%w_crit           = 0.1125_wp        ! 0.75 * capacity
        land%c_soil           = 2.0e6_wp
        land%freeze_floor     = .FALSE.
        allocate(land%mask(2,1), land%albedo(2,1), land%w(2,1), land%t_soil(2,1))
        allocate(land%skin(2,1), land%beta(2,1), land%runoff(2,1), land%fnet(2,1))
        land%mask(1,1)   = .TRUE.
        land%mask(2,1)   = .FALSE.
        land%albedo      = 0.2_wp
        land%w           = 0.05_wp
        land%t_soil      = 288.0_wp
        land%skin        = 288.0_wp
        land%beta        = 1.0_wp
        land%runoff      = 0.0_wp
        land%fnet        = 0.0_wp
        return
    end subroutine make_land

    function flx(nx, ny, val) result(a)
        ! A constant (nx,ny) flux field.
        integer,  intent(in) :: nx, ny
        real(wp), intent(in) :: val
        real(wp) :: a(nx,ny)
        a = val
        return
    end function flx

    function evap2(val) result(a)
        real(wp), intent(in) :: val
        real(wp) :: a(2,1)
        a = val
        return
    end function evap2

    function prc2(val) result(a)
        real(wp), intent(in) :: val
        real(wp) :: a(2,1)
        a = val
        return
    end function prc2

    subroutine check_disabled(nfail)
        ! A disabled land must touch nothing: the seam calls return at once and
        ! nothing is allocated. This is the module-level guarantee behind the
        ! model-level bit-for-bit "l_land = .false." invariant.
        integer, intent(inout) :: nfail
        type(aeros_land_class) :: dland
        real(wp) :: sstd(3,2)
        dland%enabled = .FALSE.
        dland%nlon = 3; dland%nlat = 2
        sstd = 290.0_wp
        call aeros_land_pre(dland, sstd)
        call aeros_land_step(dland, flx(3,2,1.0_wp), flx(3,2,1.0_wp), &
                             flx(3,2,1.0_wp), flx(3,2,1.0_wp), &
                             flx(3,2,1.0_wp), flx(3,2,1.0_wp), 1800.0_wp)
        call check(.not. allocated(dland%skin) .and. .not. allocated(dland%w) &
                   .and. .not. allocated(dland%t_soil), &
                   "disabled land allocates nothing and its seam calls are no-ops", nfail)
        return
    end subroutine check_disabled

    subroutine check_mask_from_era5(nfail)
        ! Read the ERA5 land-sea mask through the bcinput regridder onto a coarse
        ! regular lon/lat target and threshold it. Pins that the mask is read the
        ! right way up (a known land point is land, a known ocean point ocean)
        ! and that the global land fraction is physically plausible.
        integer, intent(inout) :: nfail
        character(len=*), parameter :: lsm_file = &
            "/Users/alrobi001/data/era5/monthly-single-levels/"// &
            "era5_monthly-single-levels_lsm_1991-2020_clim.nc"
        integer, parameter :: nx = 144, ny = 72
        real(wp) :: tlon(nx), tlat(ny)
        real(wp), allocatable :: lsm(:,:)
        logical,  allocatable :: mask(:,:)
        logical  :: ex
        real(wp) :: frac
        integer  :: i, j, il, jl, io, jo

        inquire(file=lsm_file, exist=ex)
        if (.not. ex) then
            write(*,*) "  skip : ERA5 lsm file absent, mask-read check skipped"
            return
        end if

        do i = 1, nx
            tlon(i) = 360.0_wp*real(i-1, wp)/real(nx, wp)         ! 0 .. ~357.5 E
        end do
        do j = 1, ny
            tlat(j) = 89.0_wp - 178.0_wp*real(j-1, wp)/real(ny-1, wp)  ! 89 .. -89 N
        end do

        allocate(lsm(nx,ny), mask(nx,ny))
        call aeros_bcinput_read_field(lsm_file, "lsm", tlon, tlat, lsm)
        mask = aeros_land_is_land(lsm, 0.5_wp)

        frac = real(count(mask), wp)/real(nx*ny, wp)
        call check(frac > 0.15_wp .and. frac < 0.45_wp, &
                   "ERA5 land fraction is physically plausible", nfail)
        write(*,"(a40,f8.3)") "   ERA5 land fraction ", frac

        ! Central Africa (~20 E, 5 N) is land; central Pacific (~200 E, 0 N) is ocean.
        il = nearest(tlon, nx, 20.0_wp); jl = nearest(tlat, ny, 5.0_wp)
        io = nearest(tlon, nx, 200.0_wp); jo = nearest(tlat, ny, 0.0_wp)
        call check(mask(il,jl),        "central Africa classified as land", nfail)
        call check(.not. mask(io,jo),  "central Pacific classified as ocean", nfail)

        deallocate(lsm, mask)
        return
    end subroutine check_mask_from_era5

    integer function nearest(x, n, target) result(kbest)
        real(wp), intent(in) :: x(:)
        integer,  intent(in) :: n
        real(wp), intent(in) :: target
        real(wp) :: d, dbest
        integer  :: k
        dbest = 1.0e30_wp; kbest = 1
        do k = 1, n
            d = abs(x(k) - target)
            if (d < dbest) then; dbest = d; kbest = k; end if
        end do
        return
    end function nearest

    subroutine check(ok, label, nfail)
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

end program test_land
