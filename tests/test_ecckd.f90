program test_ecckd
    ! Acceptance test for the ecCKD correlated-k radiation (SCHEME_ECCKD): the
    ! opt-in longwave (aeros_ecckd_lw_clearsky_column) and shortwave
    ! (aeros_ecckd_sw_clearsky_column) column kernels. Locks in the Phase-1/2
    ! properties before the cloudy branch (Phase 3) changes anything.
    !
    ! Checks (exit non-zero on any failure):
    !   1. LW ENERGY CLOSURE. Column-integrated heating equals the net LW flux
    !      convergence (surface - TOA) to machine precision -- the telescoping
    !      identity that makes the flux-divergence heating conservative.
    !   2. SW ENERGY CLOSURE. Column-integrated heating equals the atmospheric
    !      absorption (TOA-down minus TOA-up minus surface-net), exactly.
    !   3. PHYSICAL OLR + CO2 MONOTONICITY. OLR sits in a sane midlatitude band,
    !      and doubling CO2 raises the greenhouse (OLR drops) by ~3-4 W/m2.
    !   4. SW REDUCES TO SESAM. The ecCKD SW is SESAM's validated clear-sky SW
    !      with ozone upgraded from a constant band transmission to a correlated-k;
    !      the reused Rayleigh + near-IR H2O structure must keep the two within the
    !      small (a few percent) intended ozone difference -- a regression guard on
    !      the reuse. Asserted to 3% relative on TOA-up, surface-down, surface-net.
    !
    ! This drives the SCHEME_ECCKD kernels explicitly; the SESAM default path is
    ! unchanged (its own tests cover it).

    use aeros_defs,      only : dp, wp, p0, grav, cp_d, sigma_sb, S0, pi, &
                                aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_radiation,    only : aeros_sw_clearsky_column, aeros_insolation_daily
    use aeros_ecckd,        only : aeros_ecckd_lw_clearsky_column, &
                                   aeros_ecckd_sw_clearsky_column, &
                                   aeros_ecckd_lw_cloudy_column, &
                                   aeros_ecckd_sw_cloudy_column

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
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0)") " test_ecckd:: correlated-k LW+SW, single column, L", nlev

    ! --- a midlatitude column, with a stratospheric ozone layer ---------------
    ps = real(p0, wp)
    call aeros_vgrid_pressure(vg, ps, phalf, pfull, dp_lev)
    do k = 1, nlev
        t(k)  = strat_troposphere(pfull(k)/ps)
        q(k)  = humidity(t(k), pfull(k), pfull(k)/ps)
        o3(k) = 1.5e-5_wp*exp(-0.5_wp*(log(pfull(k)/2000.0_wp)/1.2_wp)**2)
    end do
    ts = 288.0_wp
    call aeros_hydrostatic(vg, 0.0_wp, t, phalf, phi_full, phi_half)
    z_half = phi_half/real(grav, wp)

    call test_lw_energy(nfail)
    call test_sw_energy(nfail)
    call test_olr_and_co2(nfail)
    call test_sw_reduces_to_sesam(nfail)
    call test_clouds(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_ecckd:: PASS"
    else
        write(*,*) "test_ecckd:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    pure real(wp) function strat_troposphere(sig) result(tt)
        real(wp), intent(in) :: sig
        real(wp), parameter  :: t_sfc = 288.0_wp, t_strat = 216.5_wp, sig_tp = 0.2_wp
        if (sig > sig_tp) then
            tt = t_strat + (t_sfc - t_strat)*(sig - sig_tp)/(1.0_wp - sig_tp)
        else
            tt = t_strat
        end if
    end function strat_troposphere

    real(wp) function humidity(tk, pk, sig) result(qq)
        real(wp), intent(in) :: tk, pk, sig
        real(wp) :: qs, dqsdt, rh
        call aeros_qsat(tk, pk, qs, dqsdt)
        rh = 0.75_wp*max(0.0_wp, (sig - 0.15_wp)/(1.0_wp - 0.15_wp))
        qq = max(3.0e-6_wp, rh*qs)
    end function humidity

    pure real(wp) function co2_kgkg(ppm) result(qc)
        real(wp), intent(in) :: ppm
        qc = ppm*1.0e-6_wp * 44.0095_wp/28.97_wp
    end function co2_kgkg

    ! === 1. LW energy closure =================================================
    subroutine test_lw_energy(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: fnet(0:nlev), heat(nlev), olr, fdw_sur, col_heat, flux_conv, rel
        integer  :: kk
        write(*,*) ""
        write(*,*) " -- LW column heating equals net flux convergence"
        call aeros_ecckd_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                            co2_kgkg(280.0_wp), .TRUE., &
                                            fnet, heat, olr, fdw_sur)
        col_heat = 0.0_wp
        do kk = 1, nlev
            col_heat = col_heat + (real(cp_d, wp)/real(grav, wp))*heat(kk)*dp_lev(kk)
        end do
        flux_conv = fnet(nlev) - fnet(0)
        rel = abs(col_heat - flux_conv)/max(1.0_wp, abs(flux_conv))
        write(*,"(a,es12.3)") "   relative mismatch              ", rel
        call check(rel < 1.0e-12_wp, "LW layer heating is exactly the flux divergence", nfail)
        call check(col_heat < 0.0_wp, "the column cools radiatively (net LW)", nfail)
    end subroutine test_lw_energy

    ! === 2. SW energy closure =================================================
    subroutine test_sw_energy(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: heat(nlev), up, dw, net, cz, swdn, a_atm, col_heat, rel
        integer  :: kk
        write(*,*) ""
        write(*,*) " -- SW column heating equals the absorbed flux"
        call aeros_insolation_daily(30.0_wp*real(pi,wp)/180.0_wp, 172.0_wp, real(S0,wp), cz, swdn)
        call aeros_ecckd_sw_clearsky_column(nlev, q, o3, .TRUE., dp_lev, swdn, cz, &
                                            0.06_wp, 0.06_wp, heat, up, dw, net)
        a_atm = swdn - up - net
        col_heat = 0.0_wp
        do kk = 1, nlev
            col_heat = col_heat + (real(cp_d, wp)/real(grav, wp))*heat(kk)*dp_lev(kk)
        end do
        rel = abs(col_heat - a_atm)/max(1.0_wp, abs(a_atm))
        write(*,"(a,f8.2,a)") "   atmospheric SW absorption      ", a_atm, " W/m2"
        write(*,"(a,es12.3)") "   relative mismatch              ", rel
        call check(a_atm > 0.0_wp, "the atmosphere absorbs shortwave", nfail)
        call check(rel < 1.0e-10_wp, "SW column heating equals the absorbed flux", nfail)
    end subroutine test_sw_energy

    ! === 3. Physical OLR and CO2 forcing ======================================
    subroutine test_olr_and_co2(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: fnet(0:nlev), heat(nlev), olr1, olr2, fdw1, fdw2, forcing
        write(*,*) ""
        write(*,*) " -- physical OLR and CO2-doubling forcing"
        call aeros_ecckd_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                            co2_kgkg(280.0_wp), .TRUE., fnet, heat, olr1, fdw1)
        call aeros_ecckd_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                            co2_kgkg(560.0_wp), .TRUE., fnet, heat, olr2, fdw2)
        forcing = olr1 - olr2
        write(*,"(a,f7.2,a)") "   OLR at 280 ppm                 ", olr1, " W/m2"
        write(*,"(a,f7.2,a)") "   OLR at 560 ppm                 ", olr2, " W/m2"
        write(*,"(a,f7.3,a)") "   instantaneous TOA forcing      ", forcing, " W/m2"
        call check(olr1 > 220.0_wp .and. olr1 < 290.0_wp, &
                   "OLR is a physical midlatitude value (220-290 W/m2)", nfail)
        call check(olr1 < sigma_sb*ts**4, &
                   "greenhouse: OLR is below the surface blackbody emission", nfail)
        call check(olr2 < olr1, "doubling CO2 reduces OLR (positive forcing)", nfail)
        call check(forcing > 3.0_wp .and. forcing < 5.0_wp, &
                   "CO2-doubling forcing is ~3-4 W/m2 on the reference column", nfail)
    end subroutine test_olr_and_co2

    ! === 4. SW reduces to SESAM (regression guard on the reuse) ================
    subroutine test_sw_reduces_to_sesam(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: he(nlev), hs(nlev), cz, swdn
        real(wp) :: eup, edw, enet, sup, sdw, snet
        write(*,*) ""
        write(*,*) " -- ecCKD SW agrees with SESAM SW (ozone is the only deviation)"
        call aeros_insolation_daily(30.0_wp*real(pi,wp)/180.0_wp, 172.0_wp, real(S0,wp), cz, swdn)
        call aeros_ecckd_sw_clearsky_column(nlev, q, o3, .TRUE., dp_lev, swdn, cz, &
                                            0.06_wp, 0.06_wp, he, eup, edw, enet)
        call aeros_sw_clearsky_column(nlev, q, o3, .TRUE., dp_lev, swdn, cz, &
                                      0.06_wp, 0.06_wp, hs, sup, sdw, snet)
        write(*,"(a,3f9.3)") "   ecCKD up/dw/net               ", eup, edw, enet
        write(*,"(a,3f9.3)") "   SESAM up/dw/net               ", sup, sdw, snet
        call check(abs(eup-sup)/max(1.0_wp,sup)   < 0.03_wp, "TOA up SW within 3% of SESAM", nfail)
        call check(abs(edw-sdw)/max(1.0_wp,sdw)   < 0.03_wp, "surface down SW within 3% of SESAM", nfail)
        call check(abs(enet-snet)/max(1.0_wp,snet) < 0.03_wp, "surface net SW within 3% of SESAM", nfail)
    end subroutine test_sw_reduces_to_sesam

    ! === 5. All-sky (grey cloud folded into the g-points) =====================
    subroutine test_clouds(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: fnet(0:nlev), heat(nlev), olr_cs, fdw_cs, olr0, fdw0, olr_cld, fdw_cld
        real(wp) :: cf(nlev), clwc(nlev), ciwc(nlev)
        real(wp) :: cz, swdn, hcs(nlev), up_cs, dw_cs, net_cs
        real(wp) :: h0(nlev), up0, dw0, net0, h1(nlev), up1, dw1, net1
        real(wp) :: col_heat, flux_conv, rel, sig
        integer  :: kk

        write(*,*) ""
        write(*,*) " -- all-sky: cf=0 reduces to clear; a cloud deck warms LW / brightens SW"

        ! clear-sky references
        call aeros_ecckd_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                            co2_kgkg(280.0_wp), .TRUE., fnet, heat, olr_cs, fdw_cs)
        call aeros_insolation_daily(30.0_wp*real(pi,wp)/180.0_wp, 172.0_wp, real(S0,wp), cz, swdn)
        call aeros_ecckd_sw_clearsky_column(nlev, q, o3, .TRUE., dp_lev, swdn, cz, &
                                            0.06_wp, 0.06_wp, hcs, up_cs, dw_cs, net_cs)

        ! cf=0 must recover the clear-sky columns bit-for-bit
        cf = 0.0_wp; clwc = 0.0_wp; ciwc = 0.0_wp
        call aeros_ecckd_lw_cloudy_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                          co2_kgkg(280.0_wp), .TRUE., cf, clwc, ciwc, &
                                          fnet, heat, olr0, fdw0)
        call aeros_ecckd_sw_cloudy_column(nlev, q, o3, .TRUE., dp_lev, swdn, cz, &
                                          0.06_wp, 0.06_wp, cf, clwc, ciwc, h0, up0, dw0, net0)
        call check(olr0 == olr_cs .and. fdw0 == fdw_cs, &
                   "cf=0 recovers the clear-sky LW exactly", nfail)
        call check(up0 == up_cs .and. dw0 == dw_cs .and. net0 == net_cs, &
                   "cf=0 recovers the clear-sky SW exactly", nfail)

        ! a mid-tropospheric liquid cloud deck
        do kk = 1, nlev
            sig = pfull(kk)/ps
            if (sig > 0.4_wp .and. sig < 0.6_wp) then
                cf(kk) = 0.9_wp; clwc(kk) = 1.0e-4_wp
            end if
        end do
        call aeros_ecckd_lw_cloudy_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                          co2_kgkg(280.0_wp), .TRUE., cf, clwc, ciwc, &
                                          fnet, heat, olr_cld, fdw_cld)
        call aeros_ecckd_sw_cloudy_column(nlev, q, o3, .TRUE., dp_lev, swdn, cz, &
                                          0.06_wp, 0.06_wp, cf, clwc, ciwc, h1, up1, dw1, net1)
        write(*,"(a,f7.2,a)") "   LW cloud radiative effect      ", olr_cs - olr_cld, " W/m2"
        write(*,"(a,f7.2,a)") "   SW cloud radiative effect      ", -(up1 - up_cs), " W/m2"
        call check(olr_cld < olr_cs, "cloud reduces OLR (positive LW cloud effect)", nfail)
        call check(fdw_cld > fdw_cs, "cloud raises the surface downward LW", nfail)
        call check(up1 > up_cs .and. dw1 < dw_cs, "cloud raises albedo / dims surface SW", nfail)

        ! energy closure carries to the all-sky LW
        col_heat = 0.0_wp
        do kk = 1, nlev
            col_heat = col_heat + (real(cp_d, wp)/real(grav, wp))*heat(kk)*dp_lev(kk)
        end do
        flux_conv = fnet(nlev) - fnet(0)
        rel = abs(col_heat - flux_conv)/max(1.0_wp, abs(flux_conv))
        call check(rel < 1.0e-12_wp, "all-sky LW heating is exactly the flux divergence", nfail)
    end subroutine test_clouds

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

end program test_ecckd
