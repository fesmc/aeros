program test_seaice
    ! Acceptance test for thermodynamic sea ice (aeros_ocean, l_seaice on): the
    ! Semtner (1976) 0-layer ice + mixed-layer slab that replaces the freeze-floor
    ! clamp.
    !
    ! 1. OFF IS THE FREEZE FLOOR. With l_seaice = .FALSE. the ocean is the
    !    freeze-floor slab, bit-for-bit: a huge net loss clamps SST at T0 exactly
    !    and no ice state is allocated. This is the invariant that keeps every
    !    existing result unchanged.
    !
    ! 2. FREEZE-UP AND GROWTH. Under sustained cold forcing an open, near-freezing
    !    slab cools to the freezing point, forms ice, and the ice thickens step on
    !    step. The ice surface runs colder than the freezing base (conduction), the
    !    ice fraction is 1, and the surface albedo rises to the ice albedo.
    !
    ! 3. MELT-OUT. Under sustained warm forcing the ice thins to zero, the cell
    !    returns to open water, the mixed layer then warms above freezing, and the
    !    albedo falls back to the open-water value.
    !
    ! 4. ENERGY CONSERVATION. Across a freeze-then-melt cycle the change in system
    !    energy E = C (T_mix - t_frz) - rho_ice L_f h equals the time-integral of
    !    the net surface flux the ocean reports (ocn%fnet) -- latent heat included.
    !
    ! Exits non-zero on failure.

    use aeros_defs,  only : dp, wp, sigma_sb, T0, L_f, aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_ocean

    implicit none

    integer, parameter :: trunc = 21

    ! Must match the private RHO_ICE parameter in aeros_ocean (used only to
    ! reconstruct the system energy independently of the model's own bookkeeping).
    real(wp), parameter :: RHO_ICE = 917.0_wp

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_ocean_class) :: ocn

    integer :: nlon, nlat, nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    nlon = grd%nlon; nlat = grd%nlat

    write(*,"(a,i0,a,i0,a,i0)") " test_seaice:: T", trunc, "  grid ", nlon, "x", nlat

    call test_off_is_freeze_floor(nfail)
    call test_freeze_and_grow(nfail)
    call test_melt_out(nfail)
    call test_energy_conservation(nfail)

    call aeros_ocean_end(ocn)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_seaice:: PASS"
    else
        write(*,*) "test_seaice:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine setup_ocean(depth, seaice)
        ! A slab ocean with (optionally) thermodynamic sea ice on.
        real(wp), intent(in) :: depth
        logical,  intent(in) :: seaice
        call aeros_ocean_end(ocn)
        ocn%mode     = OCEAN_SLAB
        ocn%depth    = depth
        ocn%l_seaice = seaice
        call aeros_ocean_init(ocn, grd)
    end subroutine setup_ocean

    subroutine test_off_is_freeze_floor(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: swn(:,:), lwd(:,:), sh(:,:), lh(:,:)

        write(*,*) ""
        write(*,*) " -- l_seaice = .FALSE. is the freeze-floor slab, unchanged"

        call setup_ocean(1.0_wp, .FALSE.)
        allocate(swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat))
        swn = 0.0_wp; lwd = 0.0_wp; sh = 1.0e4_wp; lh = 0.0_wp   ! enormous loss
        call aeros_ocean_step(ocn, swn, lwd, sh, lh, 1.0e5_wp)

        call check(minval(ocn%sst) == real(T0,wp), &
                   "freeze floor still clamps SST at T0 exactly", nfail)
        call check(.not. allocated(ocn%h_ice), &
                   "no ice state allocated when l_seaice is off", nfail)
        deallocate(swn, lwd, sh, lh)
    end subroutine test_off_is_freeze_floor

    subroutine test_freeze_and_grow(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: swn(:,:), lwd(:,:), sh(:,:), lh(:,:)
        real(wp) :: h1, h2, dt
        integer  :: n

        write(*,*) ""
        write(*,*) " -- sustained cold forcing: freeze-up, then the ice grows"

        call setup_ocean(2.0_wp, .TRUE.)
        ocn%sst = ocn%t_frz + 0.5_wp        ! start just above freezing, open water
        allocate(swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat))
        dt  = 86400.0_wp                    ! 1 day
        swn = 0.0_wp; lwd = 150.0_wp; sh = 0.0_wp; lh = 0.0_wp   ! net cooling

        do n = 1, 5
            call aeros_ocean_step(ocn, swn, lwd, sh, lh, dt)
        end do
        h1 = ocn%h_ice(1,1)
        do n = 1, 15
            call aeros_ocean_step(ocn, swn, lwd, sh, lh, dt)
        end do
        h2 = ocn%h_ice(1,1)

        write(*,"(a,f8.4,a,f8.4,a)") "   h_ice after 5d ", h1, " m -> after 20d ", h2, " m"
        write(*,"(a,f8.2,a,f8.2,a)") "   ice surface ", ocn%t_ice_sfc(1,1), &
            " K vs freezing base ", ocn%t_frz, " K"

        call check(h1 > 0.0_wp, "ice has formed under cold forcing", nfail)
        call check(h2 > h1, "the ice thickens step on step", nfail)
        call check(ocn%a_ice(1,1) == 1.0_wp, "ice fraction is 1 where ice is present", nfail)
        call check(ocn%t_ice_sfc(1,1) < ocn%t_frz, &
                   "ice surface is colder than the freezing base (conduction)", nfail)
        call check(ocn%alb(1,1) == ocn%ice_albedo, &
                   "surface albedo rises to the ice albedo", nfail)
        deallocate(swn, lwd, sh, lh)
    end subroutine test_freeze_and_grow

    subroutine test_melt_out(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: swn(:,:), lwd(:,:), sh(:,:), lh(:,:)
        real(wp) :: dt
        integer  :: n

        write(*,*) ""
        write(*,*) " -- sustained warm forcing: the ice melts out, then water warms"

        call setup_ocean(2.0_wp, .TRUE.)
        ! Start fully ice-covered.
        ocn%h_ice     = 0.8_wp
        ocn%a_ice     = 1.0_wp
        ocn%t_ice_sfc = ocn%t_frz - 5.0_wp
        ocn%sst       = ocn%t_frz - 5.0_wp
        call aeros_ocean_albedo_update(ocn)
        call check(ocn%alb(1,1) == ocn%ice_albedo, "starts at the ice albedo", nfail)

        allocate(swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat))
        dt  = 86400.0_wp
        swn = 250.0_wp; lwd = 300.0_wp; sh = 0.0_wp; lh = 0.0_wp  ! strong warming

        do n = 1, 30
            call aeros_ocean_step(ocn, swn, lwd, sh, lh, dt)
        end do

        write(*,"(a,f8.4,a,f8.2,a)") "   h_ice ", ocn%h_ice(1,1), " m   SST ", &
            ocn%sst(1,1), " K"

        call check(ocn%h_ice(1,1) == 0.0_wp, "the ice has melted out", nfail)
        call check(ocn%a_ice(1,1) == 0.0_wp, "ice fraction back to 0", nfail)
        call check(ocn%sst(1,1) > ocn%t_frz, "the open mixed layer warms above freezing", nfail)
        call check(ocn%alb(1,1) == ocn%ocn_albedo, &
                   "albedo falls back to the open-water value", nfail)
        deallocate(swn, lwd, sh, lh)
    end subroutine test_melt_out

    subroutine test_energy_conservation(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: swn(:,:), lwd(:,:), sh(:,:), lh(:,:)
        real(wp) :: dt, rlf, cheat, e0, e1, eflux, relerr
        integer  :: n

        write(*,*) ""
        write(*,*) " -- energy conservation across a freeze-then-melt cycle"

        call setup_ocean(2.0_wp, .TRUE.)
        ocn%sst = ocn%t_frz + 3.0_wp         ! open, warm
        allocate(swn(nlon,nlat), lwd(nlon,nlat), sh(nlon,nlat), lh(nlon,nlat))
        dt    = 86400.0_wp
        rlf   = RHO_ICE*real(L_f, wp)
        cheat = ocn%heat_cap

        e0    = system_energy(1, 1, cheat, rlf)
        eflux = 0.0_wp

        ! Freeze phase: cool hard for 20 days.
        swn = 0.0_wp; lwd = 140.0_wp; sh = 0.0_wp; lh = 0.0_wp
        do n = 1, 20
            call aeros_ocean_step(ocn, swn, lwd, sh, lh, dt)
            eflux = eflux + ocn%fnet(1,1)*dt
        end do
        ! Melt phase: warm hard for 25 days (back to open water and beyond).
        swn = 250.0_wp; lwd = 300.0_wp; sh = 0.0_wp; lh = 0.0_wp
        do n = 1, 25
            call aeros_ocean_step(ocn, swn, lwd, sh, lh, dt)
            eflux = eflux + ocn%fnet(1,1)*dt
        end do

        e1 = system_energy(1, 1, cheat, rlf)
        relerr = abs((e1 - e0) - eflux)/max(abs(eflux), 1.0_wp)

        write(*,"(a,es13.5,a)") "   dE (state)      ", e1 - e0, " J/m2"
        write(*,"(a,es13.5,a)") "   integral F dt   ", eflux,   " J/m2"
        write(*,"(a,es13.5)")   "   relative error  ", relerr

        call check(relerr < 1.0e-8_wp, &
                   "system energy change equals the integrated surface flux", nfail)
        deallocate(swn, lwd, sh, lh)
    end subroutine test_energy_conservation

    real(wp) function system_energy(i, j, cheat, rlf) result(e)
        ! Energy of the ice + mixed-layer column at cell (i,j), relative to an
        ! ice-free mixed layer at the freezing point. Ice is a latent deficit; the
        ! mixed layer is pinned at t_frz whenever ice is present.
        integer,  intent(in) :: i, j
        real(wp), intent(in) :: cheat, rlf
        if (ocn%h_ice(i,j) > 0.0_wp) then
            e = -rlf*ocn%h_ice(i,j)
        else
            e = cheat*(ocn%sst(i,j) - ocn%t_frz)
        end if
    end function system_energy

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
    end subroutine check

end program test_seaice
