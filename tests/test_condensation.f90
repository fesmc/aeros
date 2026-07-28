program test_condensation
    ! Acceptance test for large-scale condensation.
    !
    ! Condensation is the point where the gridpoint humidity and the spectral
    ! temperature are coupled, and the whole design of the operator is to keep
    ! that coupling a closed budget rather than an approximate one. So the tests
    ! are budget equalities, not tolerances:
    !
    ! 1. A SUBSATURATED COLUMN IS UNTOUCHED. No supersaturation, nothing
    !    condenses: q unchanged, no heating added, no precipitation. The
    !    operator must do exactly nothing where there is nothing to do -- a
    !    scheme that drizzled from unsaturated air would be wrong in a way that
    !    only shows up as a slow global drying.
    !
    ! 2. ENERGY AND WATER CLOSE, PER CELL. Where it does condense, the heating
    !    added to the temperature tendency, times the step, is exactly
    !    L_v/cp_d times the vapour removed -- the latent heat of that water and
    !    no other. This is column moist static energy conservation written per
    !    gridbox, and it is what makes the coupling honest: the heat the
    !    spectral temperature will gain is the heat the gridpoint humidity lost.
    !
    ! 3. THE COLUMN REACHES SATURATION. After the adjustment the air is at (its
    !    critical) saturation, not still supersaturated and not overshot dry --
    !    the two Newton iterations have actually solved the implicit equation.
    !
    ! 4. THE WATER REMOVED IS THE PRECIPITATION. The column vapour loss, as a
    !    mass per area per time, equals the reported precipitation rate. Nothing
    !    is stored and nothing leaks.
    !
    ! 5. q STAYS >= 0. Condensation cannot remove more vapour than is present.
    !
    ! Plus a spot-check that q_sat itself is right, since everything rests on it.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, p0, cp_d, grav, L_v, &
                                aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_condensation

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg
    type(aeros_cond_class)  :: cnd

    real(wp), allocatable :: t_g(:,:,:), qv(:,:,:), qv0(:,:,:), dt_phys(:,:,:), lnps(:,:)
    integer :: nlon, nlat, nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)
    call aeros_condensation_init(cnd, grd, .TRUE.)     ! enabled, rh_crit = 1

    nlon = grd%nlon; nlat = grd%nlat

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_condensation:: T", trunc, " L", nlev, &
                                        "  grid ", nlon, "x", nlat

    allocate(t_g(nlon,nlat,nlev), qv(nlon,nlat,nlev), qv0(nlon,nlat,nlev), &
                dt_phys(nlon,nlat,nlev), lnps(nlon,nlat))

    lnps = log(real(p0,wp))

    call test_qsat(nfail)
    call test_subsaturated_untouched(nfail)
    call test_budgets(nfail)
    call test_reevaporation(nfail)

    call aeros_condensation_end(cnd)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_condensation:: PASS"
    else
        write(*,*) "test_condensation:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine set_temperature(t_g)
        ! A plausible tropospheric profile: warm surface, lapsing upward, so the
        ! saturation humidity falls sharply with height and the same vapour is
        ! subsaturated low down and supersaturated aloft.
        implicit none
        real(wp), intent(out) :: t_g(:,:,:)
        real(wp) :: sigma
        integer  :: k
        do k = 1, nlev
            sigma = (real(k,wp) - 0.5_wp)/real(nlev,wp)      ! ~ p/ps
            t_g(:,:,k) = 300.0_wp - 45.0_wp*(1.0_wp - sigma)  ! 255 K top, 300 K base
        end do
        return
    end subroutine set_temperature

    ! === q_sat spot check ====================================================

    subroutine test_qsat(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: qs1, qs2, d1, d2

        write(*,*) ""
        write(*,*) " -- q_sat: value and monotonicity"

        call aeros_qsat(293.15_wp, real(p0,wp), qs1, d1)   ! 20 C, 1000 hPa
        call aeros_qsat(303.15_wp, real(p0,wp), qs2, d2)   ! 30 C

        write(*,"(a40,f10.4,a)") "   q_sat(20 C, 1000 hPa)          ", qs1*1000.0_wp, " g/kg"
        write(*,"(a40,f10.4,a)") "   q_sat(30 C, 1000 hPa)          ", qs2*1000.0_wp, " g/kg"

        ! ~14.5 g/kg at 20 C is textbook; allow a little slack for the formula.
        call check(abs(qs1 - 0.01467_wp) < 5.0e-4_wp, "q_sat(20 C) is about 14.7 g/kg", nfail)
        call check(qs2 > qs1 .and. d1 > 0.0_wp, "q_sat and dq_sat/dT increase with T", nfail)
        return
    end subroutine test_qsat

    ! === 1. Subsaturated column untouched ====================================

    subroutine test_subsaturated_untouched(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: dq, dh, dpr

        write(*,*) ""
        write(*,*) " -- a subsaturated column is left exactly alone"

        call set_temperature(t_g)
        qv  = 1.0e-4_wp        ! 0.1 g/kg: well below saturation everywhere warm
        qv0 = qv
        dt_phys = 0.0_wp

        call aeros_condensation_apply(cnd, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)

        dq  = maxval(abs(qv - qv0))
        dh  = maxval(abs(dt_phys))
        dpr = maxval(cnd%precip)

        write(*,"(a40,es12.3)") "   max |dq|                       ", dq
        write(*,"(a40,es12.3)") "   max |heating added|           ", dh
        write(*,"(a40,es12.3)") "   max precip                    ", dpr

        call check(dq == 0.0_wp .and. dh == 0.0_wp .and. dpr == 0.0_wp, &
                    "no condensation, no heating, no precip below saturation", nfail)
        return
    end subroutine test_subsaturated_untouched

    ! === 2-5. Budgets on a supersaturated column =============================

    subroutine test_budgets(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev)
        real(wp) :: qs, dqsdt, resid, worst_sat, worst_en, dqc
        real(wp) :: water_removed, precip_water, hcp
        real(dp) :: colw
        integer  :: i, j, k

        write(*,*) ""
        write(*,*) " -- a supersaturated column: budgets close"

        hcp = real(L_v,wp)/real(cp_d,wp)

        call set_temperature(t_g)
        ! Load every level 15% above its own saturation -- a physical
        ! supersaturation the adjustment must remove, not the absurd 500% a
        ! constant q over a 45 K column would impose (that no single-step
        ! Newton solve, or real atmosphere, ever sees).
        call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
        do k = 1, nlev
            call aeros_qsat(t_g(1,1,k), pfull(k), qs, dqsdt)
            qv(:,:,k) = 1.15_wp*qs
        end do
        qv0  = qv
        dt_phys = 0.0_wp

        call aeros_condensation_apply(cnd, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)

        ! Per-cell energy-water closure and saturation, worst case over the grid.
        worst_sat = 0.0_wp
        worst_en  = 0.0_wp
        do j = 1, nlat
            do i = 1, nlon
                call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
                do k = 1, nlev
                    dqc = qv0(i,j,k) - qv(i,j,k)

                    ! energy-water: cp * heating_increment == L * dqc, relative
                    ! to the latent heat itself (which is ~1e4, so an absolute
                    ! bound would just be measuring its roundoff). dt_phys is the
                    ! forward [K] increment, so no dt factor here.
                    if (dqc > 0.0_wp) &
                        worst_en = max(worst_en, &
                                abs(real(cp_d,wp)*dt_phys(i,j,k) - real(L_v,wp)*dqc) &
                                    /(real(L_v,wp)*dqc))

                    ! saturation: where it condensed, q sits at rh_crit*q_sat(T + heating)
                    if (dqc > 0.0_wp) then
                        call aeros_qsat(t_g(i,j,k) + hcp*dqc, pfull(k), qs, dqsdt)
                        resid = abs(qv(i,j,k) - cnd%rh_crit*qs)
                        worst_sat = max(worst_sat, resid)
                    end if
                end do
            end do
        end do

        ! Column water budget at one point: vapour lost == precip * g * dt / (dp sum)
        ! Compare the removed vapour mass with the reported precip mass directly.
        i = 1; j = 1
        call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
        water_removed = 0.0_wp
        do k = 1, nlev
            water_removed = water_removed + (qv0(i,j,k) - qv(i,j,k))*dpc(k)/real(grav,wp)
        end do
        precip_water = cnd%precip(i,j)*1800.0_wp

        write(*,"(a40,es12.3)")   "   worst |cp dT - L dq| per cell  ", worst_en
        write(*,"(a40,es12.3,a)") "   worst supersaturation residual ", worst_sat, " kg/kg"
        write(*,"(a40,es12.3,a)") "   column precip                  ", &
                                    cnd%precip(i,j)*86400.0_wp, " mm/day"
        write(*,"(a40,es12.3)")   "   |water removed - precip|/water ", &
                                    abs(water_removed - precip_water)/water_removed
        write(*,"(a40,es12.3)")   "   min q after                    ", minval(qv)

        call check(worst_en < 1.0e-13_wp, &
                    "heating equals L times the vapour condensed, per cell", nfail)
        call check(worst_sat < 1.0e-6_wp, &
                    "condensing cells are brought to rh_crit * q_sat", nfail)
        call check(abs(water_removed - precip_water)/water_removed < 1.0e-12_dp, &
                    "the vapour removed equals the precipitation", nfail)
        call check(minval(qv) >= 0.0_wp, "q never goes negative", nfail)

        ! Not vacuous: something actually condensed.
        colw = real(sum(qv0 - qv), dp)
        call check(colw > 0.0_dp, "and something did condense (the test is not vacuous)", nfail)
        return
    end subroutine test_budgets

    ! === 6. Reevaporation of falling precip ==================================

    subroutine test_reevaporation(nfail)
        ! Rain condensed aloft falls into dry layers below and reevaporates:
        ! it moistens them, reduces the surface precip, and closes the column
        ! water and energy budgets exactly (the evaporated water leaves the flux
        ! and enters the vapour; its latent heat cools the air).
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev)
        real(wp) :: qs, dqsdt, w0, w1, precip_m, en_heat, dry0, dry1
        real(wp) :: precip_off, precip_on
        integer  :: k

        write(*,*) ""
        write(*,*) " -- reevaporation: rain moistens dry layers below and conserves"

        call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
        call set_temperature(t_g)

        ! supersaturate the top third (the rain source); dry the layers below it
        ! (the reevaporation sink)
        do k = 1, nlev
            call aeros_qsat(t_g(1,1,k), pfull(k), qs, dqsdt)
            if (k <= nlev/3) then
                qv(:,:,k) = 1.6_wp*qs
            else
                qv(:,:,k) = 0.2_wp*qs
            end if
        end do
        qv0  = qv
        dry0 = qv0(1,1,nlev)

        ! baseline: reevaporation off -> all condensate reaches the surface
        cnd%reevap = 0.0_wp
        dt_phys = 0.0_wp; qv = qv0
        call aeros_condensation_apply(cnd, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)
        precip_off = cnd%precip(1,1)

        ! reevaporation on
        cnd%reevap = 30.0_wp
        dt_phys = 0.0_wp; qv = qv0
        call aeros_condensation_apply(cnd, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)
        precip_on = cnd%precip(1,1)
        dry1 = qv(1,1,nlev)

        ! column water and energy budgets (reevaporation on)
        w0 = 0.0_wp; w1 = 0.0_wp; en_heat = 0.0_wp
        do k = 1, nlev
            w0 = w0 + qv0(1,1,k)*dpc(k)/real(grav,wp)
            w1 = w1 + qv(1,1,k)*dpc(k)/real(grav,wp)
            en_heat = en_heat + real(cp_d,wp)*dt_phys(1,1,k)*dpc(k)/real(grav,wp)
        end do
        precip_m = precip_on*1800.0_wp    ! kg m-2 over the step

        write(*,"(a40,es12.3,a)") "   precip, reevap off             ", precip_off*86400.0_wp, " mm/day"
        write(*,"(a40,es12.3,a)") "   precip, reevap on              ", precip_on*86400.0_wp, " mm/day"
        write(*,"(a40,es12.3)")   "   dry lowest layer, dq (moisten) ", dry1 - dry0
        write(*,"(a40,es12.3)")   "   |water in - (vapour+precip)|/in", abs(w0 - (w1 + precip_m))/w0
        write(*,"(a40,es12.3)")   "   |heating - L*precip|/(L*precip)", &
                                    abs(en_heat - real(L_v,wp)*precip_m)/abs(real(L_v,wp)*precip_m)

        call check(abs(w0 - (w1 + precip_m))/w0 < 1.0e-12_wp, &
                    "column water conserved (vapour + surface precip)", nfail)
        call check(dry1 > dry0, "reevaporation moistens a dry layer below the rain", nfail)
        call check(precip_on < precip_off, "reevaporation reduces the surface precip", nfail)
        call check(abs(en_heat - real(L_v,wp)*precip_m)/abs(real(L_v,wp)*precip_m) < 1.0e-12_wp, &
                    "net latent heating equals L * surface precip", nfail)
        call check(minval(qv) >= 0.0_wp, "q stays non-negative", nfail)
        return
    end subroutine test_reevaporation

    subroutine check(ok, label, nfail)
        implicit none
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

end program test_condensation
