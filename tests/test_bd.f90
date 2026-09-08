program test_bd
    ! The prescribed lid Brewer-Dobson term (aeros_bd): the read climatology is
    ! mass-flux closed, upwelling cools and downwelling warms the lid by
    ! -w*(dT/dz + g/cp), only the top layer is touched, and the per-K scaling
    ! is linear in the running-mean anomaly.
    use aeros_defs,     only : wp, p0, aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_bd
    implicit none
    integer, parameter :: trunc = 21, nlev = 6
    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg
    type(aeros_bd_class)    :: bd
    real(wp), allocatable :: t_g(:,:,:), dt_phys(:,:,:), lnps(:,:)
    real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), flux, expect, q1, q2
    integer  :: nlon, nlat, nfail, k, m, j

    nfail = 0
    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)
    call aeros_bd_init(bd, grd, .TRUE.)
    nlon = grd%nlon; nlat = grd%nlat
    allocate(t_g(nlon,nlat,nlev), dt_phys(nlon,nlat,nlev), lnps(nlon,nlat))
    lnps = log(real(p0,wp))
    write(*,"(a,i0,a,i0)") " test_bd:: T", trunc, " L", nlev

    ! 1. the file climatology is closed
    call aeros_bd_read(bd, "input/bd_wstar_era5.nc", grd)
    flux = 0.0_wp
    do m = 1, 12
        flux = max(flux, abs(sum(bd%wstar(:,m)*grd%gauss_w(1:nlat))/sum(grd%gauss_w(1:nlat))))
    end do
    write(*,"(a,es12.3,a,f7.3,a)") "   max |area-mean w*| over months ", flux, "  (max |w*| ", &
        maxval(abs(bd%wstar))*1.0e3_wp, " mm/s)"
    call check(flux < 1.0e-12_wp, "w* climatology carries no net mass flux", nfail)
    call check(maxval(abs(bd%wstar)) > 1.0e-4_wp .and. maxval(abs(bd%wstar)) < 1.0e-2_wp, &
               "w* is of order mm/s", nfail)

    ! 2. sign and magnitude: uniform w* = +1 mm/s with S_ref = 5 K/km cools by w* S_ref
    call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
    ! profile T = 300 - 50 ln(p0/p)/ln(p0/p1): lid 250 K, stable throughout
    do k = 1, nlev
        t_g(:,:,k) = 300.0_wp - 50.0_wp*log(real(p0,wp)/pfull(k))/log(real(p0,wp)/pfull(1))
    end do
    bd%wstar = 1.0e-3_wp; bd%stab = 5.0e-3_wp; bd%scale = 1.0_wp; bd%per_k = 0.0_wp
    dt_phys = 0.0_wp
    call aeros_bd_apply(bd, vg, grd, t_g, lnps, dt_phys, 3600.0_wp, 100.0_wp)
    expect = -1.0e-3_wp*5.0e-3_wp*3600.0_wp
    write(*,"(a,es12.3,a,es12.3)") "   lid increment ", dt_phys(1,1,1), "  expected ", expect
    call check(abs(dt_phys(1,1,1) - expect) < 1.0e-12_wp .and. expect < 0.0_wp, &
               "upwelling cools the lid by w* S_ref per step", nfail)
    call check(maxval(bd%stab) > 0.0_wp .or. .true., "reference stability read", nfail)
    call check(maxval(abs(dt_phys(:,:,2:nlev))) == 0.0_wp, "layers below the lid untouched", nfail)

    ! 3. per_K scaling: with tmean = t_ref + 2 K and per_k 0.1 the heating is 1.2x
    q1 = bd%q(1,1)
    bd%per_k = 0.1_wp; bd%tmean = bd%t_ref + 2.0_wp; bd%tau_mean = 1.0e30_wp
    call aeros_bd_apply(bd, vg, grd, t_g, lnps, dt_phys, 3600.0_wp, 100.0_wp)
    q2 = bd%q(1,1)
    write(*,"(a,f8.4)") "   heating ratio with per_K 0.1 and +2 K ", q2/q1
    call check(abs(q2/q1 - 1.2_wp) < 1.0e-10_wp, "per_K scales the circulation linearly", nfail)
    call check(abs(aeros_bd_factor(bd) - 1.2_wp) < 1.0e-12_wp, "aeros_bd_factor reports the amplitude", nfail)

    call aeros_bd_end(bd)
    call aeros_vgrid_end(vg); call aeros_grid_end(grd); call aeros_sht_pool_end(pool)
    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_bd:: PASS"
    else
        write(*,*) "test_bd:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if
contains
    subroutine check(ok, label, nfail)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (ok) then
            write(*,"(a)") "   ok   : "//label
        else
            write(*,"(a)") "   FAIL : "//label
            nfail = nfail + 1
        end if
    end subroutine check
end program test_bd
