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
    use aeros_defs,         only : S0, pi
    use aeros_radiation,    only : aeros_lw_clearsky_column, &
                                   aeros_lw_cloudy_column, &
                                   aeros_sw_clearsky_column, aeros_sw_cloudy_column, &
                                   aeros_insolation_daily
    use aeros_cloud,        only : aeros_cloud_diagnose

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
    call test_insolation(nfail)
    call test_shortwave(nfail)
    call test_cloudy_lw(nfail)
    call test_cloudy_sw(nfail)
    call test_cloud_diagnose(nfail)

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

    ! === 6. Insolation ========================================================

    subroutine test_insolation(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: cz, sw, cz2, sw2, annual, lat, doy
        integer  :: d

        write(*,*) ""
        write(*,*) " -- daily-mean insolation"

        ! equator at equinox: daily mean ~ S0/pi ~ 433 W/m2
        call aeros_insolation_daily(0.0_wp, 80.0_wp, real(S0,wp), cz, sw)
        write(*,"(a,f7.2,a)") "   equator equinox insolation     ", sw, " W/m2"
        call check(abs(sw - real(S0,wp)/real(pi,wp)) < 15.0_wp, &
                   "equator equinox daily insolation ~ S0/pi", nfail)
        call check(cz > 0.0_wp .and. cz <= 1.0_wp, &
                   "airmass cosine zenith is in (0,1]", nfail)

        ! north pole in NH winter (doy ~ 355): polar night, zero insolation
        call aeros_insolation_daily(80.0_wp*real(pi,wp)/180.0_wp, 355.0_wp, real(S0,wp), cz2, sw2)
        write(*,"(a,f7.2,a)") "   80N midwinter insolation       ", sw2, " W/m2"
        call check(sw2 == 0.0_wp .and. cz2 == 0.0_wp, &
                   "polar night gives zero insolation", nfail)

        ! annual-and-global mean insolation ~ S0/4 ~ 340 W/m2
        annual = 0.0_wp
        do d = 1, 360, 5
            do lat = -85.0_wp, 85.0_wp, 10.0_wp
                doy = real(d, wp)
                call aeros_insolation_daily(lat*real(pi,wp)/180.0_wp, doy, real(S0,wp), cz, sw)
                annual = annual + sw*cos(lat*real(pi,wp)/180.0_wp)
            end do
        end do
        annual = annual/count_area()
        write(*,"(a,f7.2,a)") "   annual-global mean insolation  ", annual, " W/m2"
        call check(abs(annual - real(S0,wp)/4.0_wp) < 12.0_wp, &
                   "annual-global mean insolation ~ S0/4", nfail)
        return
    end subroutine test_insolation

    real(wp) function count_area() result(a)
        real(wp) :: lat
        integer  :: d
        a = 0.0_wp
        do d = 1, 360, 5
            do lat = -85.0_wp, 85.0_wp, 10.0_wp
                a = a + cos(lat*real(pi,wp)/180.0_wp)
            end do
        end do
        return
    end function count_area

    ! === 7. Shortwave =========================================================

    subroutine test_shortwave(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: cz, swdn, heat(nlev)
        real(wp) :: sw_up, sw_dw, sw_net, col_heat, a_atm, peak_kday
        integer  :: kk

        write(*,*) ""
        write(*,*) " -- clear-sky shortwave (ocean surface, 30N midsummer)"

        call aeros_insolation_daily(30.0_wp*real(pi,wp)/180.0_wp, 172.0_wp, real(S0,wp), cz, swdn)
        call aeros_sw_clearsky_column(nlev, q, o3, .FALSE., dp_lev, swdn, cz, &
                                      0.06_wp, 0.06_wp, heat, sw_up, sw_dw, sw_net)

        a_atm = swdn - sw_up - sw_net
        write(*,"(a,f7.2,a)") "   TOA down SW                    ", swdn, " W/m2"
        write(*,"(a,f7.3,a)") "   planetary albedo               ", sw_up/swdn, " "
        write(*,"(a,f7.2,a)") "   surface down SW                ", sw_dw, " W/m2"
        write(*,"(a,f7.2,a)") "   atmospheric SW absorption      ", a_atm, " W/m2"

        call check(sw_up/swdn > 0.05_wp .and. sw_up/swdn < 0.4_wp, &
                   "planetary albedo is physical for a dark ocean", nfail)
        call check(sw_dw > 0.0_wp .and. sw_dw < swdn, &
                   "surface down SW is positive and below TOA", nfail)
        call check(a_atm > 0.0_wp, "the atmosphere absorbs shortwave", nfail)

        ! column heating positive, and it equals the absorbed flux exactly
        col_heat = 0.0_wp
        peak_kday = 0.0_wp
        do kk = 1, nlev
            col_heat = col_heat + (real(cp_d, wp)/real(grav, wp))*heat(kk)*dp_lev(kk)
            peak_kday = max(peak_kday, heat(kk)*86400.0_wp)
        end do
        write(*,"(a,f7.3,a)") "   peak SW heating                ", peak_kday, " K/day"
        write(*,"(a,es12.3)") "   |col heating - absorption|/abs ", abs(col_heat - a_atm)/a_atm

        call check(peak_kday > 0.1_wp .and. peak_kday < 4.0_wp, &
                   "shortwave warms the column at a fraction of a K/day", nfail)
        call check(abs(col_heat - a_atm)/a_atm < 1.0e-12_wp, &
                   "column SW heating equals the absorbed flux", nfail)
        return
    end subroutine test_shortwave

    ! === 8. All-sky (cloudy) longwave ========================================

    subroutine test_cloudy_lw(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: fnet(0:nlev), heat(nlev), olr_cs, fdw_cs
        real(wp) :: olr0, fdw0, olr_cld, fdw_cld, olr_thick, fdw_thick
        real(wp) :: cf(nlev), clwc(nlev), ciwc(nlev)
        real(wp) :: col_heat, flux_conv, rel, sig
        integer  :: kk

        write(*,*) ""
        write(*,*) " -- all-sky longwave (mid-tropospheric cloud)"

        ! clear-sky reference
        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      co2_kgkg(280.0_wp), .FALSE., &
                                      fnet, heat, olr_cs, fdw_cs)

        ! zero cloud must recover the clear-sky kernel bit-for-bit
        cf = 0.0_wp; clwc = 0.0_wp; ciwc = 0.0_wp
        call aeros_lw_cloudy_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                    co2_kgkg(280.0_wp), .FALSE., cf, clwc, ciwc, &
                                    fnet, heat, olr0, fdw0)
        call check(olr0 == olr_cs .and. fdw0 == fdw_cs, &
                   "zero cloud recovers the clear-sky OLR/down-LW exactly", nfail)

        ! a liquid cloud deck in the mid troposphere (0.4 < sigma < 0.6)
        do kk = 1, nlev
            sig = pfull(kk)/ps
            if (sig > 0.4_wp .and. sig < 0.6_wp) then
                cf(kk)   = 0.9_wp
                clwc(kk) = 1.0e-4_wp        ! ~40 g/m2 layer water path
            end if
        end do
        call aeros_lw_cloudy_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                    co2_kgkg(280.0_wp), .FALSE., cf, clwc, ciwc, &
                                    fnet, heat, olr_cld, fdw_cld)

        write(*,"(a,f7.2,a)") "   clear-sky OLR                  ", olr_cs,  " W/m2"
        write(*,"(a,f7.2,a)") "   all-sky OLR (cloud)            ", olr_cld, " W/m2"
        write(*,"(a,f7.2,a)") "   LW cloud radiative effect TOA  ", olr_cs - olr_cld, " W/m2"
        write(*,"(a,f7.2,a)") "   clear-sky surface down LW      ", fdw_cs,  " W/m2"
        write(*,"(a,f7.2,a)") "   all-sky surface down LW        ", fdw_cld, " W/m2"

        call check(olr_cld < olr_cs, &
                   "cloud reduces OLR (positive LW cloud radiative effect)", nfail)
        call check(fdw_cld > fdw_cs, &
                   "cloud increases the downward LW at the surface", nfail)

        ! energy identity carries over to the all-sky path
        col_heat = 0.0_wp
        do kk = 1, nlev
            col_heat = col_heat + (real(cp_d, wp)/real(grav, wp))*heat(kk)*dp_lev(kk)
        end do
        flux_conv = fnet(nlev) - fnet(0)
        rel = abs(col_heat - flux_conv)/max(1.0_wp, abs(flux_conv))
        call check(rel < 1.0e-12_wp, &
                   "all-sky layer heating is exactly the flux divergence", nfail)

        ! a thicker/more opaque deck cools OLR further (monotone in optical depth)
        do kk = 1, nlev
            sig = pfull(kk)/ps
            if (sig > 0.4_wp .and. sig < 0.6_wp) then
                cf(kk)   = 1.0_wp
                clwc(kk) = 5.0e-4_wp
            end if
        end do
        call aeros_lw_cloudy_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                    co2_kgkg(280.0_wp), .FALSE., cf, clwc, ciwc, &
                                    fnet, heat, olr_thick, fdw_thick)
        call check(olr_thick < olr_cld, &
                   "a thicker cloud reduces OLR further (monotone in opacity)", nfail)
        return
    end subroutine test_cloudy_lw

    ! === 9. All-sky (cloudy) shortwave =======================================

    subroutine test_cloudy_sw(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: cz, swdn
        real(wp) :: heat_cs(nlev), up_cs, dw_cs, net_cs
        real(wp) :: heat0(nlev), up0, dw0, net0
        real(wp) :: heat1(nlev), up1, dw1, net1
        real(wp) :: heat2(nlev), up2, dw2, net2
        real(wp) :: cf(nlev), clwc(nlev), ciwc(nlev)
        real(wp) :: a_atm, col_heat, rel, sig
        integer  :: kk

        write(*,*) ""
        write(*,*) " -- all-sky shortwave (mid-tropospheric cloud, 30N midsummer)"

        call aeros_insolation_daily(30.0_wp*real(pi,wp)/180.0_wp, 172.0_wp, real(S0,wp), cz, swdn)

        ! clear-sky reference
        call aeros_sw_clearsky_column(nlev, q, o3, .FALSE., dp_lev, swdn, cz, &
                                      0.06_wp, 0.06_wp, heat_cs, up_cs, dw_cs, net_cs)

        ! zero cloud must recover the clear-sky column exactly
        cf = 0.0_wp; clwc = 0.0_wp; ciwc = 0.0_wp
        call aeros_sw_cloudy_column(nlev, q, o3, .FALSE., dp_lev, swdn, cz, &
                                    0.06_wp, 0.06_wp, cf, clwc, ciwc, &
                                    heat0, up0, dw0, net0)
        call check(up0 == up_cs .and. dw0 == dw_cs .and. net0 == net_cs, &
                   "zero cloud recovers the clear-sky SW fluxes exactly", nfail)

        ! a liquid cloud deck in the mid troposphere
        do kk = 1, nlev
            sig = pfull(kk)/ps
            if (sig > 0.4_wp .and. sig < 0.6_wp) then
                cf(kk)   = 0.9_wp
                clwc(kk) = 1.0e-4_wp
            end if
        end do
        call aeros_sw_cloudy_column(nlev, q, o3, .FALSE., dp_lev, swdn, cz, &
                                    0.06_wp, 0.06_wp, cf, clwc, ciwc, &
                                    heat1, up1, dw1, net1)

        write(*,"(a,f7.3,a)") "   clear-sky planetary albedo     ", up_cs/swdn, " "
        write(*,"(a,f7.3,a)") "   all-sky planetary albedo       ", up1/swdn, " "
        write(*,"(a,f7.2,a)") "   SW cloud radiative effect TOA  ", -(up1-up_cs), " W/m2"
        write(*,"(a,f7.2,a)") "   clear-sky surface down SW      ", dw_cs, " W/m2"
        write(*,"(a,f7.2,a)") "   all-sky surface down SW        ", dw1, " W/m2"

        call check(up1 > up_cs, &
                   "cloud raises the planetary albedo (negative SW CRE)", nfail)
        call check(dw1 < dw_cs .and. net1 < net_cs, &
                   "cloud reduces the surface shortwave", nfail)

        ! energy: column absorption equals the TOA-minus-surface residual
        a_atm = swdn - up1 - net1
        col_heat = 0.0_wp
        do kk = 1, nlev
            col_heat = col_heat + (real(cp_d, wp)/real(grav, wp))*heat1(kk)*dp_lev(kk)
        end do
        rel = abs(col_heat - a_atm)/max(1.0_wp, abs(a_atm))
        write(*,"(a,es12.3)") "   |col SW heating - absorption|  ", rel
        call check(a_atm > 0.0_wp .and. rel < 1.0e-10_wp, &
                   "all-sky column SW heating equals the absorbed flux", nfail)

        ! a thicker cloud reflects more (monotone albedo)
        do kk = 1, nlev
            sig = pfull(kk)/ps
            if (sig > 0.4_wp .and. sig < 0.6_wp) then
                cf(kk)   = 1.0_wp
                clwc(kk) = 5.0e-4_wp
            end if
        end do
        call aeros_sw_cloudy_column(nlev, q, o3, .FALSE., dp_lev, swdn, cz, &
                                    0.06_wp, 0.06_wp, cf, clwc, ciwc, &
                                    heat2, up2, dw2, net2)
        call check(up2 > up1, &
                   "a thicker cloud raises the albedo further (monotone)", nfail)
        return
    end subroutine test_cloudy_sw

    ! === 10. Diagnostic cloud scheme =========================================

    subroutine test_cloud_diagnose(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: cf(nlev), clwc(nlev), ciwc(nlev)
        real(wp) :: fnet(0:nlev), heat(nlev), olr_cs, fdw_cs, olr_cld, fdw_cld
        real(wp) :: cfmax, clwsum, valid
        integer  :: kk

        write(*,*) ""
        write(*,*) " -- diagnostic cloud scheme"

        call aeros_cloud_diagnose(nlev, t, q, pfull, ps, cf, clwc, ciwc)

        cfmax = 0.0_wp; clwsum = 0.0_wp; valid = 1.0_wp
        do kk = 1, nlev
            cfmax  = max(cfmax, cf(kk))
            clwsum = clwsum + clwc(kk) + ciwc(kk)
            if (cf(kk) < 0.0_wp .or. cf(kk) > 1.0_wp) valid = 0.0_wp
            if (clwc(kk) < 0.0_wp .or. ciwc(kk) < 0.0_wp) valid = 0.0_wp
        end do
        write(*,"(a,f7.3)")   "   max diagnosed cloud fraction   ", cfmax
        write(*,"(a,es12.3)") "   column cloud water (l+i)       ", clwsum

        call check(valid > 0.5_wp, "cloud fraction in [0,1] and water non-negative", nfail)
        call check(cfmax > 0.0_wp .and. cfmax <= 1.0_wp, &
                   "the moist column diagnoses some cloud", nfail)
        call check(clwsum > 0.0_wp, "cloudy layers carry cloud water", nfail)

        ! the diagnosed cloud, fed to the all-sky kernel, warms the LW (reduces OLR)
        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      co2_kgkg(280.0_wp), .FALSE., &
                                      fnet, heat, olr_cs, fdw_cs)
        call aeros_lw_cloudy_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                    co2_kgkg(280.0_wp), .FALSE., cf, clwc, ciwc, &
                                    fnet, heat, olr_cld, fdw_cld)
        call check(olr_cld < olr_cs, &
                   "diagnosed cloud reduces OLR through the all-sky kernel", nfail)
        return
    end subroutine test_cloud_diagnose

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
