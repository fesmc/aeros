program test_ocean
    ! Acceptance test for the sea surface (aeros_ocean): the prescribed SST
    ! profile, and the slab's surface energy balance.
    !
    ! 1. SST PROFILE. The initial/prescribed SST is the APE control profile:
    !    warmest at the equator, frozen poleward of 60 deg, zonally symmetric.
    !
    ! 2. PRESCRIBED MODE IS INERT. A step in prescribed mode leaves the SST
    !    exactly unchanged, whatever the fluxes -- so every prescribed-SST run is
    !    bit-identical to before the slab existed.
    !
    ! 3. SLAB SIGN. A net heat gain (SW in > everything out) warms the slab; a net
    !    loss cools it. The change is exactly F_net dt / C.
    !
    ! 4. SLAB EQUILIBRIUM. With the fluxes set so F_net = 0 the SST does not move.
    !
    ! 5. FREEZE FLOOR. A strong net loss cannot drive the SST below T0 when the
    !    floor is on.
    !
    ! Exits non-zero on failure.

    use aeros_defs,  only : dp, wp, sigma_sb, T0, aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_ocean

    implicit none

    integer, parameter :: trunc = 21

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_ocean_class) :: ocn

    integer :: nlon, nlat, nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    nlon = grd%nlon; nlat = grd%nlat

    write(*,"(a,i0,a,i0,a,i0)") " test_ocean:: T", trunc, "  grid ", nlon, "x", nlat

    call test_profile(nfail)
    call test_prescribed_inert(nfail)
    call test_slab_sign(nfail)
    call test_slab_equilibrium(nfail)
    call test_freeze_floor(nfail)

    call aeros_ocean_end(ocn)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_ocean:: PASS"
    else
        write(*,*) "test_ocean:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine test_profile(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: t_eq, t_pole, zonal_spread
        integer  :: jeq, jpole, j

        write(*,*) ""
        write(*,*) " -- initial SST is the APE control profile"

        call aeros_ocean_init(ocn, grd)            ! default: prescribed, APE
        jeq   = minloc(abs(grd%lat), 1)
        jpole = maxloc(abs(grd%lat), 1)
        t_eq   = ocn%sst(1, jeq)
        t_pole = ocn%sst(1, jpole)
        zonal_spread = 0.0_wp
        do j = 1, nlat
            zonal_spread = max(zonal_spread, maxval(ocn%sst(:,j)) - minval(ocn%sst(:,j)))
        end do

        write(*,"(a,f7.2,a)") "   SST nearest equator            ", t_eq, " K"
        write(*,"(a,f7.2,a)") "   SST nearest pole               ", t_pole, " K"

        call check(t_eq > real(T0,wp) + 20.0_wp, "equatorial SST is tropical (> 20 C)", nfail)
        call check(t_pole <= real(T0,wp) + 1.0_wp, "polar SST is at freezing", nfail)
        call check(t_eq > t_pole, "SST decreases from equator to pole", nfail)
        call check(zonal_spread == 0.0_wp, "SST is zonally symmetric", nfail)
        return
    end subroutine test_profile

    subroutine test_prescribed_inert(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: sst0(:,:), swn(:,:), lwd(:,:), sh(:,:), lh(:,:)

        write(*,*) ""
        write(*,*) " -- prescribed mode: a step never moves the SST"

        call aeros_ocean_init(ocn, grd)            ! prescribed
        allocate(sst0(nlon,nlat), swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat))
        sst0 = ocn%sst
        swn = 300.0_wp; lwd = 300.0_wp; sh = 20.0_wp; lh = 80.0_wp    ! big imbalance
        call aeros_ocean_step(ocn, swn, lwd, sh, lh, 1800.0_wp)

        call check(maxval(abs(ocn%sst - sst0)) == 0.0_wp, &
                   "prescribed SST is bit-unchanged by a step", nfail)
        deallocate(sst0, swn, lwd, sh, lh)
        return
    end subroutine test_prescribed_inert

    subroutine test_slab_sign(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: swn(:,:), lwd(:,:), sh(:,:), lh(:,:)
        real(wp) :: sst0, fnet, expect, dt
        integer  :: ie, je

        write(*,*) ""
        write(*,*) " -- slab: net gain warms, and by exactly F_net dt / C"

        ocn%mode = OCEAN_SLAB
        ocn%depth = 10.0_wp
        call aeros_ocean_init(ocn, grd)
        allocate(swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat))
        dt = 1800.0_wp
        ie = 1; je = minloc(abs(grd%lat), 1)
        sst0 = ocn%sst(ie,je)

        ! a clear net heat gain into the slab
        swn = 400.0_wp; lwd = 400.0_wp; sh = 0.0_wp; lh = 0.0_wp
        fnet = swn(ie,je) + lwd(ie,je) - sigma_sb*sst0**4 - sh(ie,je) - lh(ie,je)
        expect = sst0 + fnet*dt/ocn%heat_cap

        call aeros_ocean_step(ocn, swn, lwd, sh, lh, dt)

        write(*,"(a,es13.5)") "   F_net at equator [W/m2]        ", fnet
        write(*,"(a,es13.5)") "   |SST - (SST0 + Fnet dt/C)|     ", abs(ocn%sst(ie,je) - expect)

        call check(fnet > 0.0_wp .and. ocn%sst(ie,je) > sst0, "net gain warms the slab", nfail)
        call check(abs(ocn%sst(ie,je) - expect) < 1.0e-10_wp, &
                   "SST change is exactly F_net dt / C", nfail)
        deallocate(swn, lwd, sh, lh)
        return
    end subroutine test_slab_sign

    subroutine test_slab_equilibrium(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: swn(:,:), lwd(:,:), sh(:,:), lh(:,:), sst0(:,:)
        integer  :: i, j

        write(*,*) ""
        write(*,*) " -- slab: F_net = 0 leaves the SST put"

        ocn%mode = OCEAN_SLAB
        call aeros_ocean_init(ocn, grd)
        allocate(swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat), sst0(nlon,nlat))
        sst0 = ocn%sst
        ! choose SH so that SW + LW - sigma T^4 - SH - LH = 0 exactly, per cell
        swn = 0.0_wp; lwd = 0.0_wp; lh = 0.0_wp
        do j = 1, nlat
            do i = 1, nlon
                sh(i,j) = -sigma_sb*sst0(i,j)**4     ! F_net = 0
            end do
        end do
        call aeros_ocean_step(ocn, swn, lwd, sh, lh, 1800.0_wp)

        write(*,"(a,es13.5)") "   max |SST - SST0|               ", maxval(abs(ocn%sst - sst0))
        call check(maxval(abs(ocn%sst - sst0)) < 1.0e-9_wp, &
                   "zero net flux leaves the SST unchanged", nfail)
        deallocate(swn, lwd, sh, lh, sst0)
        return
    end subroutine test_slab_equilibrium

    subroutine test_freeze_floor(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: swn(:,:), lwd(:,:), sh(:,:), lh(:,:)

        write(*,*) ""
        write(*,*) " -- slab: the freeze floor holds SST >= T0"

        ocn%mode = OCEAN_SLAB
        ocn%depth = 1.0_wp                 ! thin slab so one huge step overshoots
        ocn%freeze_floor = .TRUE.
        call aeros_ocean_init(ocn, grd)
        allocate(swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat))
        ! enormous net loss: no SW/LW in, huge sensible loss
        swn = 0.0_wp; lwd = 0.0_wp; sh = 1.0e4_wp; lh = 0.0_wp
        call aeros_ocean_step(ocn, swn, lwd, sh, lh, 1.0e5_wp)

        write(*,"(a,f8.2,a)") "   min SST after huge loss        ", minval(ocn%sst), " K"
        call check(minval(ocn%sst) >= real(T0,wp), "freeze floor holds SST >= T0", nfail)
        call check(minval(ocn%sst) == real(T0,wp), "and the floor actually engaged", nfail)
        deallocate(swn, lwd, sh, lh)
        return
    end subroutine test_freeze_floor

    subroutine check(ok, label, nfail)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (ok) then
            write(*,"(a,a)") "   ok   : ", label
        else
            write(*,"(a,a)") "   FAIL : ", label
            nfail = nfail + 1
        end if
        return
    end subroutine check

end program test_ocean
