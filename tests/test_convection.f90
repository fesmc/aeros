program test_convection
    ! Acceptance test for the moist convective adjustment, both schemes.
    !
    ! The default scheme is Simplified Betts-Miller (Frierson 2007); Manabe is
    ! retained as the parameter-free reference. The properties that matter are
    ! the column budgets -- and they are equalities, because a convection scheme
    ! that conserves only approximately heats or dries the whole atmosphere over
    ! a long run.
    !
    ! Simplified Betts-Miller (sbm_adjust):
    !   1. A STABLE COLUMN (no buoyant layer) is left untouched, no precip.
    !   2. DEEP convection (a warm, humid, conditionally-unstable column):
    !      conserves int (c_p T + L q) dp to machine precision -- the heating
    !      balances the latent heat of the water removed -- precipitates a
    !      positive amount equal to the column drying exactly, creates no
    !      supersaturation, and warms the free troposphere. q stays >= 0.
    !   3. SHALLOW convection (a conditionally-unstable but DRY column): does not
    !      rain, conserves both column moisture and int (c_p T + L q) dp, and
    !      still redistributes (non-vacuous). q stays >= 0.
    !   4. RELAXATION CONVERGES: applied step after step to a deep column, the
    !      run never NaNs, conserves energy every step, drives the precipitation
    !      toward zero (the column neutralizes), and keeps q >= 0.
    !
    ! Manabe (manabe_adjust), kept for the reference scheme:
    !   5. MOIST budgets: conserves int (c_p T + L q) dp, precipitates what it
    !      removes, leaves the column moist-stable, q >= 0.
    !   6. DRY branch: a dry-unstable column mixes potential temperature,
    !      conserving int c_p T dp, with no precipitation.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, p0, cp_d, grav, L_v, R_d, kappa, &
                                aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_convection

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 20

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg
    type(aeros_conv_class)  :: cnv

    real(wp), allocatable :: t_g(:,:,:), qv(:,:,:), qv0(:,:,:), dt_phys(:,:,:), lnps(:,:)
    real(wp) :: pf(nlev), ph(0:nlev), dpc(nlev)
    integer :: nlon, nlat, nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)
    call aeros_convection_init(cnv, grd, .TRUE.)

    nlon = grd%nlon; nlat = grd%nlat
    call aeros_vgrid_pressure(vg, real(p0,wp), ph, pf, dpc)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_convection:: T", trunc, " L", nlev, &
                                        "  grid ", nlon, "x", nlat

    allocate(t_g(nlon,nlat,nlev), qv(nlon,nlat,nlev), qv0(nlon,nlat,nlev), &
                dt_phys(nlon,nlat,nlev), lnps(nlon,nlat))
    lnps = log(real(p0,wp))

    ! --- Simplified Betts-Miller (the default) ---
    cnv%scheme = SCHEME_SBM
    cnv%tau    = 7200.0_wp
    cnv%rh_ref = 0.7_wp
    call test_sbm_stable(nfail)
    call test_sbm_deep(nfail)
    call test_sbm_shallow(nfail)
    call test_sbm_converge(nfail)

    ! --- Manabe (the reference scheme) ---
    cnv%scheme = SCHEME_MANABE
    call test_manabe_moist(nfail)
    call test_manabe_dry(nfail)

    call aeros_convection_end(cnv)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_convection:: PASS"
    else
        write(*,*) "test_convection:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    real(dp) function col_mse(t, q, dpc) result(h)
        ! int (c_p T + L q) dp, the conserved moist static energy (minus the
        ! geopotential term, which the adjustment holds fixed).
        implicit none
        real(wp), intent(in) :: t(:), q(:), dpc(:)
        integer :: k
        h = 0.0_dp
        do k = 1, nlev
            h = h + (real(cp_d,dp)*real(t(k),dp) + real(L_v,dp)*real(q(k),dp))*real(dpc(k),dp)
        end do
    end function col_mse

    real(dp) function col_enthalpy(t, dpc) result(h)
        implicit none
        real(wp), intent(in) :: t(:), dpc(:)
        integer :: k
        h = 0.0_dp
        do k = 1, nlev
            h = h + real(cp_d,dp)*real(t(k),dp)*real(dpc(k),dp)
        end do
    end function col_enthalpy

    real(dp) function col_water(q, dpc) result(w)
        ! int q dp / g, the column vapour mass per area.
        implicit none
        real(wp), intent(in) :: q(:), dpc(:)
        integer :: k
        w = 0.0_dp
        do k = 1, nlev
            w = w + real(q(k),dp)*real(dpc(k),dp)
        end do
        w = w/real(grav,dp)
    end function col_water

    logical function moist_stable(t, q, pf, dpc) result(ok)
        ! True if every saturated adjacent pair is moist-neutral or stable
        ! (saturated MSE non-decreasing upward).
        implicit none
        real(wp), intent(in) :: t(:), q(:), pf(:), dpc(:)
        real(wp) :: phi(nlev), qsk, qsk1, d, hk, hk1
        integer  :: k
        ok = .TRUE.
        phi(nlev) = 0.0_wp
        do k = nlev-1, 1, -1
            phi(k) = phi(k+1) + R_d*0.5_wp*(t(k)+t(k+1))*log(pf(k+1)/pf(k))
        end do
        do k = 1, nlev-1
            call aeros_qsat(t(k),   pf(k),   qsk,  d)
            call aeros_qsat(t(k+1), pf(k+1), qsk1, d)
            if (q(k) >= 0.999_wp*qsk .and. q(k+1) >= 0.999_wp*qsk1) then
                hk  = cp_d*t(k)   + phi(k)   + L_v*qsk
                hk1 = cp_d*t(k+1) + phi(k+1) + L_v*qsk1
                if (hk < hk1 - 100.0_wp) ok = .FALSE.
            end if
        end do
    end function moist_stable

    ! === Simplified Betts-Miller ============================================

    ! 1. A stable column is left untouched.
    subroutine test_sbm_stable(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: qs, d, dtmax, dqmax
        integer  :: k

        write(*,*) ""
        write(*,*) " -- SBM: a column with no buoyant layer is left alone"

        ! Near-isothermal, strongly stable: the surface moist adiabat cools with
        ! height far faster than this environment, so no level is buoyant.
        do k = 1, nlev
            t_g(:,:,k) = 300.0_wp - 5.0_wp*(1.0_wp - pf(k)/real(p0,wp))
            call aeros_qsat(t_g(1,1,k), pf(k), qs, d)
            qv(:,:,k) = 0.4_wp*qs
        end do
        dt_phys = 0.0_wp
        qv0 = qv

        call aeros_convection_apply(cnv, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)

        dtmax = maxval(abs(dt_phys))
        dqmax = maxval(abs(qv - qv0))
        write(*,"(a40,es12.3)") "   max |dT|                       ", dtmax
        write(*,"(a40,es12.3)") "   max |dq|                       ", dqmax
        write(*,"(a40,es12.3)") "   max precip                     ", maxval(cnv%precip)
        call check(dtmax == 0.0_wp .and. dqmax == 0.0_wp .and. maxval(cnv%precip) == 0.0_wp, &
                    "a stable column is unchanged, no precipitation", nfail)
        return
    end subroutine test_sbm_stable

    ! 2. Deep convection: conserves, precipitates what it removes, no supersat.
    subroutine test_sbm_deep(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: qs, d, tcol(nlev), qcol(nlev), rhmax, midwarm
        real(dp) :: h0, h1, dmse, water_removed, precip_water
        integer  :: k

        write(*,*) ""
        write(*,*) " -- SBM: deep convection over a warm, humid, unstable column"

        call deep_column()
        dt_phys = 0.0_wp
        qv0 = qv
        do k = 1, nlev
            tcol(k) = t_g(1,1,k); qcol(k) = qv(1,1,k)
        end do
        h0 = col_mse(tcol, qcol, dpc)

        call aeros_convection_apply(cnv, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)

        do k = 1, nlev
            tcol(k) = t_g(1,1,k) + dt_phys(1,1,k)
            qcol(k) = qv(1,1,k)
        end do
        h1   = col_mse(tcol, qcol, dpc)
        dmse = abs(h1 - h0)/abs(h0)

        water_removed = col_water(qv0(1,1,:), dpc) - col_water(qv(1,1,:), dpc)
        precip_water  = real(cnv%precip(1,1),dp)*1800.0_dp

        ! No supersaturation created, and the mid-troposphere warmed.
        rhmax = 0.0_wp
        do k = 1, nlev
            call aeros_qsat(tcol(k), pf(k), qs, d)
            rhmax = max(rhmax, qcol(k)/qs)
        end do
        midwarm = dt_phys(1,1,nlev/2)

        write(*,"(a40,es12.3)")   "   |dMSE|/MSE                    ", dmse
        write(*,"(a40,es12.3,a)") "   precip                        ", &
                                    cnv%precip(1,1)*86400.0_wp, " mm/day"
        write(*,"(a40,es12.3)")   "   |water removed - precip|/water", &
                                    abs(water_removed - precip_water)/abs(water_removed)
        write(*,"(a40,f12.4)")    "   max RH after                  ", rhmax
        write(*,"(a40,f12.4)")    "   mid-column warming dT          ", midwarm
        write(*,"(a40,es12.3)")   "   min q                         ", minval(qv)

        call check(dmse < 1.0e-13_dp, "moist static energy is conserved", nfail)
        call check(precip_water > 0.0_dp, "deep convection precipitates (not vacuous)", nfail)
        call check(abs(water_removed - precip_water)/abs(water_removed) < 1.0e-12_dp, &
                    "precip equals the water removed", nfail)
        call check(rhmax <= 1.0_wp + 1.0e-6_wp, "no supersaturation is created", nfail)
        call check(midwarm > 0.0_wp, "the free troposphere warms", nfail)
        call check(minval(qv) >= 0.0_wp, "q stays non-negative", nfail)
        return
    end subroutine test_sbm_deep

    ! 3. Shallow convection: non-precipitating, conserves moisture and energy.
    subroutine test_sbm_shallow(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: qs, d, tcol(nlev), qcol(nlev), dtmax
        real(dp) :: h0, h1, dmse, w0, w1
        integer  :: k

        write(*,*) ""
        write(*,*) " -- SBM: shallow convection, moist boundary layer under dry air"

        ! Same unstable temperature and a moist boundary layer (so the parcel is
        ! buoyant and a band exists), but a dry free troposphere: q < q_ref
        ! through the cloud layer, so Pq <= 0 and the shallow branch is taken.
        call deep_column()
        do k = 1, nlev
            call aeros_qsat(t_g(1,1,k), pf(k), qs, d)
            qv(:,:,k) = (0.2_wp + 0.7_wp*pf(k)/real(p0,wp))*qs
        end do
        dt_phys = 0.0_wp
        qv0 = qv
        do k = 1, nlev
            tcol(k) = t_g(1,1,k); qcol(k) = qv(1,1,k)
        end do
        h0 = col_mse(tcol, qcol, dpc)
        w0 = col_water(qv0(1,1,:), dpc)

        call aeros_convection_apply(cnv, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)

        do k = 1, nlev
            tcol(k) = t_g(1,1,k) + dt_phys(1,1,k)
            qcol(k) = qv(1,1,k)
        end do
        h1    = col_mse(tcol, qcol, dpc)
        dmse  = abs(h1 - h0)/abs(h0)
        w1    = col_water(qv(1,1,:), dpc)
        dtmax = maxval(abs(dt_phys))

        write(*,"(a40,es12.3)")   "   |dMSE|/MSE                    ", dmse
        write(*,"(a40,es12.3,a)") "   precip                        ", &
                                    cnv%precip(1,1)*86400.0_wp, " mm/day"
        write(*,"(a40,es12.3)")   "   |d column water|/water        ", abs(w1-w0)/abs(w0)
        write(*,"(a40,f12.4)")    "   max |dT| (heat redistribution)", dtmax
        write(*,"(a40,es12.3)")   "   min q                         ", minval(qv)

        call check(dmse < 1.0e-13_dp, "moist static energy is conserved", nfail)
        call check(cnv%precip(1,1) == 0.0_wp, "shallow convection does not rain", nfail)
        call check(abs(w1-w0) == 0.0_dp, "column moisture is untouched", nfail)
        call check(dtmax > 0.0_wp, "it still redistributes heat (not vacuous)", nfail)
        call check(minval(qv) >= 0.0_wp, "q stays non-negative", nfail)
        return
    end subroutine test_sbm_shallow

    ! 4. The relaxation converges: energy conserved every step, precip -> 0.
    subroutine test_sbm_converge(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(dp) :: h0, h1, worst_mse, p_first, p_last
        real(wp) :: qmin
        logical  :: nan_seen
        integer  :: k, it
        integer, parameter :: nit = 60

        write(*,*) ""
        write(*,*) " -- SBM: repeated application neutralizes and stays conservative"

        call deep_column()
        worst_mse = 0.0_dp; nan_seen = .FALSE.; qmin = huge(1.0_wp)
        p_first = 0.0_dp; p_last = 0.0_dp

        do it = 1, nit
            h0 = col_mse(t_g(1,1,:), qv(1,1,:), dpc)
            dt_phys = 0.0_wp
            call aeros_convection_apply(cnv, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)
            do k = 1, nlev
                t_g(:,:,k) = t_g(:,:,k) + dt_phys(:,:,k)
            end do
            h1 = col_mse(t_g(1,1,:), qv(1,1,:), dpc)
            worst_mse = max(worst_mse, abs(h1 - h0)/abs(h0))
            qmin = min(qmin, minval(qv))
            if (any(t_g /= t_g) .or. any(qv /= qv)) nan_seen = .TRUE.
            if (it == 1)   p_first = real(cnv%precip(1,1),dp)
            if (it == nit) p_last  = real(cnv%precip(1,1),dp)
        end do

        write(*,"(a40,l4)")     "   NaN seen                       ", nan_seen
        write(*,"(a40,es12.3)") "   worst |dMSE|/MSE over steps     ", worst_mse
        write(*,"(a40,es12.3)") "   precip step 1  [kg m-2 s-1]     ", p_first
        write(*,"(a40,es12.3)") "   precip step 60 [kg m-2 s-1]     ", p_last
        write(*,"(a40,es12.3)") "   min q over the run             ", qmin

        call check(.not. nan_seen, "no NaN over 60 applications", nfail)
        call check(worst_mse < 1.0e-13_dp, "energy is conserved every step", nfail)
        call check(p_last < 0.1_dp*p_first, "precipitation decays as the column neutralizes", nfail)
        call check(qmin >= 0.0_wp, "q stays non-negative throughout", nfail)
        return
    end subroutine test_sbm_converge

    subroutine deep_column()
        ! A warm, humid, conditionally-unstable column: surface ~300 K, a lapse
        ! steeper than moist-adiabatic aloft, 85% relative humidity. The classic
        ! deep-convecting sounding.
        implicit none
        real(wp) :: qs, d
        integer  :: k
        do k = 1, nlev
            t_g(:,:,k) = 300.0_wp - 90.0_wp*(1.0_wp - (pf(k)/real(p0,wp))**0.6_wp)
            call aeros_qsat(t_g(1,1,k), pf(k), qs, d)
            qv(:,:,k) = 0.95_wp*qs
        end do
        return
    end subroutine deep_column

    ! === Manabe (reference scheme) ==========================================

    ! 5. Moist budgets.
    subroutine test_manabe_moist(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: qs, d, tcol(nlev), qcol(nlev)
        real(dp) :: h0, h1, worst_mse, water_removed, precip_water
        real(wp) :: qminv
        logical  :: all_stable
        integer  :: k

        write(*,*) ""
        write(*,*) " -- Manabe: a conditionally unstable, saturated column convects"

        do k = 1, nlev
            t_g(:,:,k) = 245.0_wp + 55.0_wp*(pf(k)/real(p0,wp))**1.4_wp
            call aeros_qsat(t_g(1,1,k), pf(k), qs, d)
            qv(:,:,k) = qs
        end do
        dt_phys = 0.0_wp
        qv0 = qv
        do k = 1, nlev
            tcol(k) = t_g(1,1,k); qcol(k) = qv(1,1,k)
        end do
        h0 = col_mse(tcol, qcol, dpc)

        call aeros_convection_apply(cnv, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)

        do k = 1, nlev
            tcol(k) = t_g(1,1,k) + dt_phys(1,1,k)
            qcol(k) = qv(1,1,k)
        end do
        h1 = col_mse(tcol, qcol, dpc)
        worst_mse = abs(h1 - h0)/abs(h0)

        water_removed = col_water(qv0(1,1,:), dpc) - col_water(qv(1,1,:), dpc)
        precip_water  = real(cnv%precip(1,1),dp)*1800.0_dp

        all_stable = moist_stable(tcol, qcol, pf, dpc)
        qminv = minval(qv)

        write(*,"(a40,es12.3)")   "   |dMSE|/MSE                    ", worst_mse
        write(*,"(a40,es12.3,a)") "   precip                        ", &
                                    cnv%precip(1,1)*86400.0_wp, " mm/day"
        write(*,"(a40,es12.3)")   "   |water removed - precip|/water", &
                                    abs(water_removed - precip_water)/water_removed
        write(*,"(a40,l4)")       "   column moist-stable after     ", all_stable
        write(*,"(a40,es12.3)")   "   min q                         ", qminv

        call check(worst_mse < 1.0e-13_dp, "Manabe conserves moist static energy", nfail)
        call check(precip_water > 0.0_dp, "Manabe precipitates (not vacuous)", nfail)
        call check(abs(water_removed - precip_water)/water_removed < 1.0e-12_dp, &
                    "Manabe precip equals the water removed", nfail)
        call check(all_stable, "Manabe leaves the column moist-stable", nfail)
        call check(qminv >= 0.0_wp, "q stays non-negative", nfail)
        return
    end subroutine test_manabe_moist

    ! 6. Dry branch.
    subroutine test_manabe_dry(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: qs, d, tcol(nlev)
        real(dp) :: e0, e1, worst_e
        real(wp) :: th(nlev), dqmax
        logical  :: theta_ok
        integer  :: k

        write(*,*) ""
        write(*,*) " -- Manabe: a dry-unstable, unsaturated column mixes potential temperature"

        do k = 1, nlev
            t_g(:,:,k) = 300.0_wp - 80.0_wp*(1.0_wp - pf(k)/real(p0,wp))
            call aeros_qsat(t_g(1,1,k), pf(k), qs, d)
            qv(:,:,k) = 0.05_wp*qs
        end do
        dt_phys = 0.0_wp
        qv0 = qv
        do k = 1, nlev
            tcol(k) = t_g(1,1,k)
        end do
        e0 = col_enthalpy(tcol, dpc)

        call aeros_convection_apply(cnv, vg, t_g, qv, lnps, dt_phys, 1800.0_wp)

        do k = 1, nlev
            tcol(k) = t_g(1,1,k) + dt_phys(1,1,k)
            th(k)   = tcol(k)*(real(p0,wp)/pf(k))**real(kappa,wp)
        end do
        e1 = col_enthalpy(tcol, dpc)
        worst_e = abs(e1 - e0)/abs(e0)
        dqmax   = maxval(abs(qv - qv0))

        theta_ok = .TRUE.
        do k = 1, nlev-1
            if (th(k) < th(k+1) - 1.0e-3_wp) theta_ok = .FALSE.
        end do

        write(*,"(a40,es12.3)") "   |d enthalpy|/enthalpy          ", worst_e
        write(*,"(a40,es12.3)") "   max |dq| (should be ~0)        ", dqmax
        write(*,"(a40,es12.3)") "   precip (should be 0)           ", maxval(cnv%precip)
        write(*,"(a40,l4)")     "   theta stable after            ", theta_ok

        call check(worst_e < 1.0e-13_dp, "the dry adjustment conserves enthalpy", nfail)
        call check(dqmax == 0.0_wp .and. maxval(cnv%precip) == 0.0_wp, &
                    "the dry adjustment moves no water and does not rain", nfail)
        call check(theta_ok, "the column is dry-stable afterward", nfail)
        return
    end subroutine test_manabe_dry

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

end program test_convection
