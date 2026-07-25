program test_radiation
    ! Offline single-column acceptance test for the clear-sky longwave operator
    ! (aeros_lw_clearsky_column, the ported SESAM band kernel).
    !
    ! There is no analytic answer to check against here the way condensation has
    ! a closed budget: a broadband band model is an empirical fit, so the tests
    ! are that the emergent radiative quantities land where a midlatitude column
    ! puts them, that the greenhouse and CO2-forcing SIGNS are right, and that
    ! the flux divergence and the layer heating are the same energy (the one
    ! exact identity the discretization must satisfy).
    !
    ! 1. GREENHOUSE. The outgoing longwave at TOA is less than the surface
    !    blackbody emission -- the atmosphere is opaque enough to hold some in.
    !    And it is a physical number for this column: ~200-290 W/m2.
    !
    ! 2. DOWNWELLING. The surface receives a physical downward longwave,
    !    150-380 W/m2, and less than it emits (the air above is colder).
    !
    ! 3. THE TROPOSPHERE COOLS. Clear-sky longwave cools the troposphere at a
    !    few K/day -- the cooling convection and latent heating balance. The
    !    layer heating rate is negative through the troposphere and O(1) K/day.
    !
    ! 4. FLUX DIVERGENCE = HEATING. The column-integrated heating equals the net
    !    longwave flux convergence, surface-minus-TOA, to machine precision.
    !    This is the telescoping identity of the flux-divergence heating and is
    !    what makes the operator conservative.
    !
    ! 5. CO2 FORCING SIGN AND SCALE. Doubling CO2 reduces the OLR (positive
    !    forcing) by a few W/m2 -- the whole reason a concentration-aware scheme
    !    is worth porting over grey radiation.
    !
    ! Exits non-zero on failure.

    use aeros_defs,      only : dp, wp, p0, grav, cp_d, sigma_sb, aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_radiation,    only : aeros_lw_clearsky_column

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 20

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg

    real(wp) :: phalf(0:nlev), pfull(nlev), dp_lev(nlev)
    real(wp) :: phi_full(nlev), phi_half(0:nlev), z_half(0:nlev)
    real(wp) :: t(nlev), q(nlev), o3(nlev)
    real(wp) :: ts, ps
    integer  :: nfail, k

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)            ! production hybrid, 10 hPa top

    write(*,"(a,i0,a,i0)") " test_radiation:: clear-sky LW, single column, L", nlev

    ! --- one midlatitude-ish column ------------------------------------------
    ps = real(p0, wp)
    call aeros_vgrid_pressure(vg, ps, phalf, pfull, dp_lev)

    ! Temperature: 288 K surface, lapsing linearly in sigma to a 216.5 K
    ! isothermal stratosphere above 200 hPa.
    do k = 1, nlev
        t(k) = strat_troposphere(pfull(k)/ps)
    end do
    ts = 288.0_wp

    ! Humidity: RH 0.75 at the surface tapering with height, times q_sat, with a
    ! stratospheric floor. Seeded against the actual layer temperature so it is a
    ! realistic profile, not q_sat evaluated on a warm column aloft.
    do k = 1, nlev
        q(k) = humidity(t(k), pfull(k), pfull(k)/ps)
    end do

    o3 = 0.0_wp                                 ! ozone off in this first pass

    ! Interface heights from the hydrostatic integral (surface geopotential 0).
    call aeros_hydrostatic(vg, 0.0_wp, t, phalf, phi_full, phi_half)
    z_half = phi_half/real(grav, wp)

    call report_profile()

    call test_fluxes(nfail)
    call test_energy_identity(nfail)
    call test_co2_forcing(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_radiation:: PASS"
    else
        write(*,*) "test_radiation:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    pure real(wp) function strat_troposphere(sig) result(tt)
        real(wp), intent(in) :: sig            ! p/ps
        real(wp), parameter  :: t_sfc = 288.0_wp, t_strat = 216.5_wp, sig_tp = 0.2_wp
        if (sig > sig_tp) then
            tt = t_strat + (t_sfc - t_strat)*(sig - sig_tp)/(1.0_wp - sig_tp)
        else
            tt = t_strat
        end if
        return
    end function strat_troposphere

    real(wp) function humidity(tk, pk, sig) result(qq)
        real(wp), intent(in) :: tk, pk, sig
        real(wp) :: qs, dqsdt, rh
        call aeros_qsat(tk, pk, qs, dqsdt)
        rh = 0.75_wp*max(0.0_wp, (sig - 0.15_wp)/(1.0_wp - 0.15_wp))
        qq = max(3.0e-6_wp, rh*qs)             ! stratospheric floor 3 ppmv-ish
        return
    end function humidity

    subroutine report_profile()
        real(wp) :: fnet(0:nlev), heat(nlev), olr, fdw_sur
        real(wp) :: tcwv
        integer  :: kk
        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      co2_kgkg(280.0_wp), .FALSE., &
                                      fnet, heat, olr, fdw_sur)
        tcwv = 0.0_wp
        do kk = 1, nlev
            tcwv = tcwv + q(kk)*dp_lev(kk)/real(grav, wp)
        end do
        write(*,"(a,f7.2,a)")  "   surface skin T                 ", ts, " K"
        write(*,"(a,f7.2,a)")  "   column water vapour            ", tcwv, " kg/m2"
        write(*,"(a,f7.2,a)")  "   surface blackbody sigma Ts^4   ", sigma_sb*ts**4, " W/m2"
        write(*,"(a,f7.2,a)")  "   OLR (TOA, up)                  ", olr, " W/m2"
        write(*,"(a,f7.2,a)")  "   downwelling LW at surface      ", fdw_sur, " W/m2"
        return
    end subroutine report_profile

    pure real(wp) function co2_kgkg(ppm) result(qc)
        real(wp), intent(in) :: ppm
        qc = ppm*1.0e-6_wp * 44.0095_wp/28.97_wp
        return
    end function co2_kgkg

    ! === 1-3. Boundary fluxes and tropospheric cooling =======================

    subroutine test_fluxes(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: fnet(0:nlev), heat(nlev), olr, fdw_sur
        real(wp) :: bb, cool_kday, worst_trop
        integer  :: kk

        write(*,*) ""
        write(*,*) " -- boundary fluxes and tropospheric cooling"

        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      co2_kgkg(280.0_wp), .FALSE., &
                                      fnet, heat, olr, fdw_sur)
        bb = sigma_sb*ts**4

        call check(olr > 200.0_wp .and. olr < 290.0_wp, &
                   "OLR is a physical midlatitude value (200-290 W/m2)", nfail)
        call check(olr < bb, &
                   "greenhouse: OLR is below the surface blackbody emission", nfail)
        call check(fdw_sur > 150.0_wp .and. fdw_sur < bb, &
                   "surface downwelling LW is physical and below sigma Ts^4", nfail)

        ! tropospheric layers (sigma > 0.25), mean cooling rate
        cool_kday  = 0.0_wp
        worst_trop = 0.0_wp
        do kk = 1, nlev
            if (pfull(kk)/ps > 0.25_wp) then
                cool_kday = cool_kday + heat(kk)*86400.0_wp
                worst_trop = min(worst_trop, heat(kk)*86400.0_wp)
            end if
        end do
        write(*,"(a,f7.3,a)") "   most-cooling tropo layer       ", worst_trop, " K/day"

        call check(worst_trop < 0.0_wp .and. worst_trop > -5.0_wp, &
                   "the troposphere cools at a few K/day (sign and scale)", nfail)
        return
    end subroutine test_fluxes

    ! === 4. Column heating equals net flux convergence =======================

    subroutine test_energy_identity(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: fnet(0:nlev), heat(nlev), olr, fdw_sur
        real(wp) :: col_heat, flux_conv, rel
        integer  :: kk

        write(*,*) ""
        write(*,*) " -- column heating equals net LW flux convergence"

        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      co2_kgkg(280.0_wp), .FALSE., &
                                      fnet, heat, olr, fdw_sur)

        ! integral of cp/g dT/dt dp over the column
        col_heat = 0.0_wp
        do kk = 1, nlev
            col_heat = col_heat + (real(cp_d, wp)/real(grav, wp))*heat(kk)*dp_lev(kk)
        end do
        ! net upward flux: surface interface (k=nlev) minus TOA (k=0)
        flux_conv = fnet(nlev) - fnet(0)

        rel = abs(col_heat - flux_conv)/max(1.0_wp, abs(flux_conv))
        write(*,"(a,f9.2,a)") "   column LW heating (integral)   ", col_heat, " W/m2"
        write(*,"(a,f9.2,a)") "   net flux convergence sfc-TOA   ", flux_conv, " W/m2"
        write(*,"(a,es12.3)") "   relative mismatch              ", rel

        call check(col_heat < 0.0_wp, "the column cools radiatively (net)", nfail)
        call check(rel < 1.0e-12_wp, &
                   "layer heating is exactly the flux divergence", nfail)
        return
    end subroutine test_energy_identity

    ! === 5. CO2 forcing ======================================================

    subroutine test_co2_forcing(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: fnet(0:nlev), heat(nlev)
        real(wp) :: olr1, olr2, fdw1, fdw2, forcing

        write(*,*) ""
        write(*,*) " -- CO2 doubling reduces the OLR"

        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      co2_kgkg(280.0_wp), .FALSE., &
                                      fnet, heat, olr1, fdw1)
        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      co2_kgkg(560.0_wp), .FALSE., &
                                      fnet, heat, olr2, fdw2)
        forcing = olr1 - olr2

        write(*,"(a,f7.2,a)") "   OLR at 280 ppm                 ", olr1, " W/m2"
        write(*,"(a,f7.2,a)") "   OLR at 560 ppm                 ", olr2, " W/m2"
        write(*,"(a,f7.3,a)") "   instantaneous TOA forcing      ", forcing, " W/m2"

        call check(forcing > 0.5_wp .and. forcing < 8.0_wp, &
                   "doubling CO2 cuts OLR by a few W/m2 (positive forcing)", nfail)
        return
    end subroutine test_co2_forcing

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

end program test_radiation
