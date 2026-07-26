program test_surface
    ! Acceptance test for the prescribed-SST surface fluxes.
    !
    ! The surface module is the source that closes the atmosphere's budget once
    ! Held-Suarez is removed, so the tests are the sign and the bookkeeping of
    ! that source, not a tolerance. The SST it runs against comes from aeros_ocean
    ! (the profile itself is checked in test_ocean):
    !
    ! 1. WARM OCEAN, COLD DRY AIR -> FLUXES UP. Over a surface warmer and moister
    !    than the air above it, sensible heat and evaporation are both positive
    !    (into the atmosphere); the lowest layer warms and moistens.
    !
    ! 2. BULK-FLUX BOOKKEEPING. The heating added to the lowest layer is exactly
    !    the sensible flux over the layer heat capacity, and the moistening is
    !    exactly the evaporation over the layer mass -- the reported SHF and
    !    evap ARE the tendencies applied.
    !
    ! 3. EQUILIBRIUM -> NO FLUX. When the lowest-layer air is at the SST and
    !    saturated, both fluxes vanish: the surface does nothing where there is
    !    no disequilibrium.
    !
    ! 4. q STAYS >= 0. A downward moisture flux (dew) cannot remove more vapour
    !    than the layer holds.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, p0, cp_d, grav, L_v, R_d, R_v, T0, &
                                aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_surface
    use aeros_ocean

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg
    type(aeros_surf_class)  :: surf
    type(aeros_ocean_class) :: ocn

    real(wp), allocatable :: t_g(:,:,:), qv(:,:,:), u(:,:,:), v(:,:,:), dt_phys(:,:,:), lnps(:,:)
    integer :: nlon, nlat, nfail, ks

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)
    call aeros_surface_init(surf, grd, .TRUE.)
    call aeros_ocean_init(ocn, grd)               ! SST source (prescribed APE)

    nlon = grd%nlon; nlat = grd%nlat
    ks = nlev

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_surface:: T", trunc, " L", nlev, &
                                       "  grid ", nlon, "x", nlat

    allocate(t_g(nlon,nlat,nlev), qv(nlon,nlat,nlev), u(nlon,nlat,nlev), &
             v(nlon,nlat,nlev), dt_phys(nlon,nlat,nlev), lnps(nlon,nlat))
    lnps = log(real(p0,wp))

    call test_fluxes_up(nfail)
    call test_equilibrium(nfail)
    call test_no_negative_q(nfail)

    call aeros_ocean_end(ocn)
    call aeros_surface_end(surf)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_surface:: PASS"
    else
        write(*,*) "test_surface:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    ! === 1-2. Fluxes up, and the bookkeeping ==================================

    subroutine test_fluxes_up(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev)
        real(wp) :: qv0, dsh, dqm, worst_h, worst_q
        integer  :: i, j

        write(*,*) ""
        write(*,*) " -- warm ocean under cold dry air: fluxes go up, books balance"

        ! air 10 K colder than the SST and half-saturated, calm-ish wind
        do j = 1, nlat
            t_g(:,j,:) = spread(ocn%sst(:,j) - 10.0_wp, 2, nlev)
        end do
        qv0 = 1.0e-3_wp
        qv  = qv0
        u   = 3.0_wp; v = 0.0_wp
        dt_phys = 0.0_wp

        call aeros_surface_apply(surf, vg, ocn%sst, t_g, qv, lnps, u, v, dt_phys, 1800.0_wp)

        write(*,"(a,es12.3,a)") "   min SHF over grid              ", minval(surf%shf), " W/m2"
        write(*,"(a,es12.3,a)") "   min evap over grid             ", minval(surf%evap), " kg/m2/s"

        call check(minval(surf%shf) > 0.0_wp, "sensible heat flux is upward everywhere", nfail)
        call check(minval(surf%evap) > 0.0_wp, "evaporation is upward everywhere", nfail)
        call check(minval(surf%lhf) > 0.0_wp, "latent heat flux is upward everywhere", nfail)

        ! bookkeeping: the applied tendencies equal flux / layer-capacity, exactly
        worst_h = 0.0_wp; worst_q = 0.0_wp
        do j = 1, nlat
            do i = 1, nlon
                call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
                ! forward temperature increment [K] applied to the lowest layer:
                ! the sensible flux over the layer heat capacity, times the step
                dsh = surf%shf(i,j)*real(grav,wp)/(real(cp_d,wp)*dpc(ks))*1800.0_wp
                worst_h = max(worst_h, abs(dt_phys(i,j,ks) - dsh))
                ! moistening [kg/kg] applied over the step
                dqm = surf%evap(i,j)*real(grav,wp)/dpc(ks)*1800.0_wp
                worst_q = max(worst_q, abs((qv(i,j,ks) - qv0) - dqm))
            end do
        end do
        write(*,"(a,es12.3)") "   worst |dt_phys - SHF/cap|         ", worst_h
        write(*,"(a,es12.3)") "   worst |dq - evap/mass|         ", worst_q

        call check(worst_h < 1.0e-15_wp, "lowest-layer heating equals SHF over heat capacity", nfail)
        call check(worst_q < 1.0e-15_wp, "lowest-layer moistening equals evap over layer mass", nfail)

        ! only the lowest layer is touched
        call check(maxval(abs(dt_phys(:,:,1:nlev-1))) == 0.0_wp, &
                   "only the lowest layer is forced", nfail)
        return
    end subroutine test_fluxes_up

    ! === 4. Equilibrium: no disequilibrium, no flux ===========================

    subroutine test_equilibrium(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: qs, dqs, ps
        integer  :: i, j

        write(*,*) ""
        write(*,*) " -- lowest-layer air at the SST and saturated: no flux"

        ps = real(p0,wp)
        do j = 1, nlat
            do i = 1, nlon
                t_g(i,j,:) = ocn%sst(i,j)
                call aeros_qsat(ocn%sst(i,j), ps, qs, dqs)
                qv(i,j,:) = qs
            end do
        end do
        u = 5.0_wp; v = 5.0_wp
        dt_phys = 0.0_wp

        call aeros_surface_apply(surf, vg, ocn%sst, t_g, qv, lnps, u, v, dt_phys, 1800.0_wp)

        write(*,"(a,es12.3,a)") "   max |SHF|                      ", maxval(abs(surf%shf)), " W/m2"
        write(*,"(a,es12.3,a)") "   max |evap|                     ", maxval(abs(surf%evap)), " kg/m2/s"

        call check(maxval(abs(surf%shf)) < 1.0e-10_wp .and. maxval(abs(surf%evap)) < 1.0e-10_wp, &
                   "no flux when air is at the SST and saturated", nfail)
        return
    end subroutine test_equilibrium

    ! === 5. q stays non-negative under dew ====================================

    subroutine test_no_negative_q(nfail)
        integer, intent(inout) :: nfail

        write(*,*) ""
        write(*,*) " -- a strong downward moisture flux cannot drive q < 0"

        ! Surface far colder than the air, so q_sat(SST) < q and the moisture
        ! flux is downward (deposition), with a large exchange coefficient and
        ! step so it would remove more vapour than the layer holds. The guard
        ! must clip the removal exactly at q = 0.
        call chill_surface()
        t_g = real(T0,wp) + 5.0_wp        ! mild air, above the cold surface ...
        qv  = 1.0e-3_wp                    ! ... and moist relative to q_sat(SST)
        u = 30.0_wp; v = 0.0_wp
        dt_phys = 0.0_wp

        call aeros_surface_apply(surf, vg, ocn%sst, t_g, qv, lnps, u, v, dt_phys, 1.0e4_wp)

        write(*,"(a,es12.3)") "   min q after                    ", minval(qv)
        call check(minval(qv) >= 0.0_wp, "q never goes negative under a downward flux", nfail)
        ! non-vacuous: the deposition was strong enough to reach the floor
        call check(minval(qv) == 0.0_wp, "and the guard actually engaged (q hit 0)", nfail)
        return
    end subroutine test_no_negative_q

    subroutine chill_surface()
        ! Force the SST far below the air (so q_sat(SST) < q, deposition) and a
        ! large moisture exchange (so the removal would overshoot q = 0).
        ocn%sst = real(T0,wp) - 30.0_wp
        surf%c_e = 0.1_wp
        return
    end subroutine chill_surface

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

end program test_surface
