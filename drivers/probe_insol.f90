program probe_insol
    ! Seasonal-cycle validation harness. Dumps the daily-mean TOA insolation and
    ! flux-weighted mean cosine zenith versus latitude over a year, from
    ! aeros_insolation (the fesmc/insol wrapper, Laskar 2004 orbit). Prints
    ! physical sanity checks and writes output/insol_seasonal.nc for
    ! scripts/plot_insol_seasonal.jl.
    !
    !   make probe-insol && ./libaeros/bin/probe_insol.x [trunc] [time_bp_yr]
    !
    ! trunc defaults to 42, time_bp (orbital year before present, 1950 CE) to 0.

    use aeros_defs,     only : wp, dp, S0, pi, aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_insolation
    use ncio,           only : nc_create, nc_write_dim, nc_write

    implicit none

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_insol_class) :: ins

    integer  :: trunc, d, j, nlat, jn, js, jsol
    real(dp) :: time_bp
    real(wp) :: gmean, wsum, w
    character(len=64) :: arg
    real(wp), allocatable :: sw(:,:), cz(:,:), swd(:), czd(:)
    real(wp), allocatable :: lat(:), dayx(:), swann(:), czann(:)

    trunc = 42; time_bp = 0.0_dp
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg); read(arg,*) trunc
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg); read(arg,*) time_bp
    end if

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    nlat = grd%nlat
    call aeros_insol_init(ins, grd, time_bp)

    allocate(sw(nlat,DAY_YEAR), cz(nlat,DAY_YEAR), swd(nlat), czd(nlat), &
             lat(nlat), dayx(DAY_YEAR), swann(nlat), czann(nlat))
    lat = grd%lat(1:nlat)

    do d = 1, DAY_YEAR
        call aeros_insol_day(ins, swd, czd, real(d, wp))
        sw(:,d) = swd; cz(:,d) = czd; dayx(d) = real(d, wp)
    end do
    call aeros_insol_annual(ins, swann, czann)

    ! --- sanity checks -----------------------------------------------------
    write(*,"(a,i0,a,f8.1,a)") " probe_insol:: trunc=", trunc, "  time_bp=", time_bp, " yr"

    ! annual-and-global mean insolation ~ S0/4
    gmean = 0.0_wp; wsum = 0.0_wp
    do j = 1, nlat
        w = cos(lat(j)*pi/180.0_wp)
        gmean = gmean + swann(j)*w; wsum = wsum + w
    end do
    gmean = gmean/wsum
    write(*,"(a,f7.2,a,f7.2,a)") "   annual-global mean insolation ", gmean, &
        " W/m2   (S0/4 = ", real(S0,wp)/4.0_wp, ")"

    ! NH summer solstice (~day 172): NH pole in 24-h daylight, SH pole in night
    jsol = 172
    jn = maxloc(lat, dim=1)          ! northernmost row
    js = minloc(lat, dim=1)          ! southernmost row
    write(*,"(a,f6.1,a,f7.2,a,f6.1,a,f7.2,a)") &
        "   solstice d172:  lat", lat(jn), " insol ", sw(jn,jsol), &
        " W/m2 ;  lat", lat(js), " insol ", sw(js,jsol), " W/m2"
    write(*,"(a)") "   (NH pole should be bright, SH pole dark)"

    ! --- dump --------------------------------------------------------------
    call nc_create("output/insol_seasonal.nc")
    call nc_write_dim("output/insol_seasonal.nc", "lat", x=lat, units="degrees_north")
    call nc_write_dim("output/insol_seasonal.nc", "day", x=dayx, units="day of year")
    call nc_write("output/insol_seasonal.nc", "sw_toa", sw, dim1="lat", dim2="day", &
        units="W/m2", long_name="daily-mean TOA downward SW")
    call nc_write("output/insol_seasonal.nc", "coszen", cz, dim1="lat", dim2="day", &
        units="1", long_name="flux-weighted mean cosine zenith")
    call nc_write("output/insol_seasonal.nc", "sw_ann", swann, dim1="lat", &
        units="W/m2", long_name="annual-mean insolation")
    write(*,"(a)") " probe_insol:: wrote output/insol_seasonal.nc"

    call aeros_insol_end(ins)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

end program probe_insol
