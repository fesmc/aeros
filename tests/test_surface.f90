program test_surface
    ! Acceptance test for the surface exchange (aeros_surface diagnoses the
    ! conductances; aeros_vdiff applies the fluxes as the implicit lower
    ! boundary of its column solve). The vertical diffusion itself is OFF here,
    ! so the solve is diagonal above the lowest layer and only the surface
    ! boundary acts: what is measured is the surface flux alone.
    !
    ! 1. A WARM, WET SURFACE UNDER COOLER, DRIER AIR FLUXES UPWARD everywhere,
    !    and the reported fluxes ARE what the column received: the lowest
    !    layer's dry-static-energy gain equals SHF*dt*g/dp and its moisture gain
    !    equals E*dt*g/dp, to round-off; nothing above the lowest layer moves.
    !
    ! 2. AN AIR COLUMN IN EQUILIBRIUM WITH ITS SURFACE (equal dry static energy,
    !    saturated at the skin) exchanges nothing.
    !
    ! 3. A DOWNWARD MOISTURE FLUX (cold skin under moist air) never drives q
    !    negative -- the implicit boundary cannot -- while the flux is finite.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, p0, cp_d, grav, L_v, R_d, R_v, T0, &
                                aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_surface
    use aeros_vdiff
    use aeros_ocean

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 8
    real(wp), parameter :: DT = 1800.0_wp

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg
    type(aeros_surf_class)  :: surf
    type(aeros_vdiff_class) :: vd
    type(aeros_ocean_class) :: ocn

    real(wp), allocatable :: t_g(:,:,:), qv(:,:,:), u(:,:,:), v(:,:,:), lnps(:,:)
    integer :: nlon, nlat, nfail, ks

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)
    call aeros_surface_init(surf, grd, .TRUE.)
    call aeros_vdiff_init(vd, grd, .FALSE.)       ! interior diffusion off: surface boundary only
    call aeros_ocean_init(ocn, grd)               ! SST source (prescribed APE)

    nlon = grd%nlon; nlat = grd%nlat
    ks = nlev

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_surface:: T", trunc, " L", nlev, &
                                       "  grid ", nlon, "x", nlat

    allocate(t_g(nlon,nlat,nlev), qv(nlon,nlat,nlev), u(nlon,nlat,nlev), &
             v(nlon,nlat,nlev), lnps(nlon,nlat))
    lnps = log(real(p0,wp))

    call test_fluxes_up(nfail)
    call test_equilibrium(nfail)
    call test_no_negative_q(nfail)

    call aeros_ocean_end(ocn)
    call aeros_vdiff_end(vd)
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

    subroutine column_state(tsfc_offset, rh)
        ! Isothermal-ish column: T = SST - tsfc_offset at every level, wind 5 m/s,
        ! humidity at relative humidity rh of the level-1 temperature.
        real(wp), intent(in) :: tsfc_offset, rh
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), qs, dqs
        integer  :: i, j, k
        do j = 1, nlat
            do i = 1, nlon
                call aeros_vgrid_pressure(vg, exp(lnps(i,j)), phalf, pfull, dpc)
                do k = 1, nlev
                    t_g(i,j,k) = ocn%sst(i,j) - tsfc_offset
                    call aeros_qsat(t_g(i,j,k), pfull(k), qs, dqs)
                    qv(i,j,k)  = rh*qs
                end do
            end do
        end do
        u = 5.0_wp; v = 0.0_wp
    end subroutine column_state

    subroutine apply()
        call aeros_surface_apply(surf, vg, ocn%sst, t_g, qv, lnps, u, v)
        call aeros_vdiff_apply(vd, vg, t_g, qv, u, v, lnps, DT, surf)
    end subroutine apply

    ! === 1. Fluxes upward, and consistent with the column change ===========
    subroutine test_fluxes_up(nfail)
        integer, intent(inout) :: nfail
        real(wp), allocatable :: t0(:,:,:), q0(:,:,:)
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), dsh, dqm, worst_h, worst_q, upper
        integer  :: i, j

        write(*,*) ""
        write(*,*) " -- a warm, wet surface under cooler, drier air (5 K, 70% RH)"

        call column_state(5.0_wp, 0.7_wp)
        allocate(t0(nlon,nlat,nlev), q0(nlon,nlat,nlev)); t0 = t_g; q0 = qv
        call apply()

        write(*,"(a,es12.3,a)") "   min SHF over grid              ", minval(surf%shf), " W/m2"
        write(*,"(a,es12.3,a)") "   min evap over grid             ", minval(surf%evap), " kg/m2/s"
        call check(minval(surf%shf) > 0.0_wp, "sensible heat flux is upward everywhere", nfail)
        call check(minval(surf%evap) > 0.0_wp, "evaporation is upward everywhere", nfail)
        call check(minval(surf%lhf) > 0.0_wp, "latent heat flux is upward everywhere", nfail)

        ! The column received exactly the reported flux (no geopotential change
        ! across the solve: the DSE gain is cp dT).
        worst_h = 0.0_wp; worst_q = 0.0_wp
        do j = 1, nlat
            do i = 1, nlon
                call aeros_vgrid_pressure(vg, exp(lnps(i,j)), phalf, pfull, dpc)
                dsh = surf%shf(i,j)*real(grav,wp)/(real(cp_d,wp)*dpc(ks))*DT
                dqm = surf%evap(i,j)*real(grav,wp)/dpc(ks)*DT
                worst_h = max(worst_h, abs((t_g(i,j,ks) - t0(i,j,ks)) - dsh)/max(abs(dsh), 1.0e-30_wp))
                worst_q = max(worst_q, abs((qv(i,j,ks)  - q0(i,j,ks)) - dqm)/max(abs(dqm), 1.0e-30_wp))
            end do
        end do
        upper = max(maxval(abs(t_g(:,:,1:nlev-1) - t0(:,:,1:nlev-1))), &
                    maxval(abs(qv(:,:,1:nlev-1)  - q0(:,:,1:nlev-1))))
        write(*,"(a,es12.3)") "   worst |dT - SHF/(cp m)| / |.|   ", worst_h
        write(*,"(a,es12.3)") "   worst |dq - E/m| / |.|          ", worst_q
        call check(worst_h < 1.0e-10_wp, "lowest-layer heating equals SHF over heat capacity", nfail)
        call check(worst_q < 1.0e-10_wp, "lowest-layer moistening equals evap over layer mass", nfail)
        call check(upper < 1.0e-9_wp, "layers above the lowest are untouched (diffusion off)", nfail)
        deallocate(t0, q0)
    end subroutine test_fluxes_up

    ! === 2. Equilibrium ======================================================
    subroutine test_equilibrium(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), phi_full(nlev), phi_half(0:nlev), qs, dqs
        integer  :: i, j, k, it

        write(*,*) ""
        write(*,*) " -- an air column in equilibrium with its surface exchanges nothing"

        ! Equal dry static energy: T_1 = SST - (Phi_1 - Phi_s)/cp; saturated at the skin.
        u = 5.0_wp; v = 0.0_wp
        do j = 1, nlat
            do i = 1, nlon
                call aeros_vgrid_pressure(vg, exp(lnps(i,j)), phalf, pfull, dpc)
                t_g(i,j,:) = ocn%sst(i,j)
                call aeros_hydrostatic(vg, 0.0_wp, t_g(i,j,:), phalf, phi_full, phi_half)
                do k = 1, nlev
                    t_g(i,j,k) = ocn%sst(i,j) - (phi_full(k) - phi_half(nlev))/real(cp_d,wp)
                end do
                ! iterate so Phi matches the corrected T (fixed point)
                do it = 1, 6
                    call aeros_hydrostatic(vg, 0.0_wp, t_g(i,j,:), phalf, phi_full, phi_half)
                    do k = 1, nlev
                        t_g(i,j,k) = ocn%sst(i,j) - (phi_full(k) - phi_half(nlev))/real(cp_d,wp)
                    end do
                end do
                call aeros_qsat(ocn%sst(i,j), exp(lnps(i,j)), qs, dqs)
                qv(i,j,:) = qs
            end do
        end do
        call apply()
        write(*,"(a,es12.3,a)") "   max |SHF|                      ", maxval(abs(surf%shf)), " W/m2"
        write(*,"(a,es12.3,a)") "   max |evap|                     ", maxval(abs(surf%evap)), " kg/m2/s"
        call check(maxval(abs(surf%shf)) < 1.0e-3_wp .and. maxval(abs(surf%evap)) < 1.0e-12_wp, &
                    "no flux at equilibrium (DSE-balanced, saturated column)", nfail)
    end subroutine test_equilibrium

    ! === 3. Downward moisture flux never makes q negative ===================
    subroutine test_no_negative_q(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: qmin0

        write(*,*) ""
        write(*,*) " -- a cold skin under moist air: dew, and q stays >= 0"

        call column_state(-30.0_wp, 0.9_wp)      ! air 30 K WARMER than the skin
        qmin0 = minval(qv)
        call apply()
        write(*,"(a,es12.3,a)") "   min evap (negative = dew)      ", minval(surf%evap), " kg/m2/s"
        write(*,"(a,es12.3)")   "   min q after the step           ", minval(qv)
        call check(minval(surf%evap) < 0.0_wp, "the flux is downward (dew)", nfail)
        call check(minval(qv) >= 0.0_wp, "q never goes negative under a downward flux", nfail)
        call check(all(qv == qv) .and. minval(qv) < qmin0, "q decreased and stayed finite", nfail)
    end subroutine test_no_negative_q

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
    end subroutine check

end program test_surface
