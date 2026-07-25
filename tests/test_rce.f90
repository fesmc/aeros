program test_rce
    ! Coupled radiative-convective run: the full moist-physics stack with
    ! Held-Suarez removed and replaced by a real energy budget.
    !
    ! This is the payoff of the radiation and surface work. test_moist_run keeps
    ! Held-Suarez as an artificial thermal forcing; here that is gone. The energy
    ! source is the prescribed-SST surface (sensible + latent fluxes), the sink
    ! is radiation (net longwave cooling), and convection carries the surface
    ! fluxes up the column -- the radiative-convective equilibrium loop. Nothing
    ! relaxes the temperature toward a prescribed target; the climate is whatever
    ! the physics balances to.
    !
    ! A short run cannot reach equilibrium -- the radiative timescale is weeks --
    ! so what is asserted is that the closed loop is STABLE and every part of it
    ! is ACTIVE and BOUNDED:
    !
    !   STABILITY. No NaN in T or q, bounded wind, q >= 0 -- across a stack where
    !   radiation and surface fluxes both feed the centered temperature path and
    !   convection the forward-split one.
    !
    !   EVERY COMPONENT ACTIVE. The surface fluxes move heat and moisture, the
    !   radiation cools (physical global-mean OLR), and it rains. A dead loop
    !   would pass "stable" vacuously; these check the loop is closed.
    !
    !   BOUNDED, NOT RUNNING AWAY. Without a thermal relaxation an unstable
    !   coupling would drift the global-mean temperature without limit; it stays
    !   in a physical band, and the TOA imbalance is finite.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, p0, R_d, grav, S0, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_timestep

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12
    integer, parameter :: nstep = 400

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)     :: grd
    type(aeros_vgrid_class)    :: vg
    type(aeros_param_class)    :: par
    type(aeros_state_class)    :: now
    type(aeros_timestep_class) :: ts

    integer :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_rce:: T", trunc, " L", nlev, &
                                       "  grid ", grd%nlon, "x", grd%nlat

    call run(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_rce:: PASS"
    else
        write(*,*) "test_rce:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine run(nfail)
        implicit none
        integer, intent(inout) :: nfail

        real(wp), allocatable :: phis2(:,:)
        real(wp) :: qs, dqsdt, tval
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev)
        real(dp) :: precip_tot, umax, t1_0, t1_1
        real(dp) :: olr_gm, swup_gm, swtoa_gm, shf_gm, lhf_gm, toa_net
        logical  :: nan_seen, neg_seen
        integer  :: i, j, k, n

        par%trunc = trunc; par%nlon = -1; par%nlat = -1; par%nlev = nlev
        par%nthreads = -1
        par%dt = 1800.0_wp
        par%semi_implicit = .TRUE.
        par%held_suarez   = .FALSE.      ! <-- no artificial thermal forcing
        par%eps_filter = 0.06_wp
        par%raw_alpha  = 0.53_wp
        par%ndiff = 6; par%tau_diff = 6.0_wp
        par%mass_fixer = .FALSE.

        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        ! The full moist-physics stack, on directly (no namelist in this test).
        ts%cnd%enabled = .TRUE.
        ts%cnd%rh_crit = 1.0_wp
        ts%cnv%enabled = .TRUE.
        ts%surf%enabled = .TRUE.         ! prescribed-SST surface fluxes: the source
        ts%rad%enabled  = .TRUE.         ! clear-sky LW+SW: the sink
        ! annual-mean insolation, 280 ppm, ocean albedo -- init defaults; the
        ! insolation was precomputed in aeros_radiation_init.

        allocate(phis2(grd%nlon, grd%nlat))
        phis2 = 0.0_wp
        call aeros_timestep_set_phis(ts, phis2)

        ! Same warm, humid, conditionally-unstable start as test_moist_run, so
        ! the columns convect from the outset. It is NOT the equilibrium; the
        ! run relaxes away from it under the real fluxes.
        call aeros_spec_zero(now%spec)
        do k = 1, nlev
            tval = 300.0_wp - 90.0_wp*(1.0_wp - ((real(k,wp) - 0.5_wp)/real(nlev,wp))**0.6_wp)
            now%spec%temp(aeros_sht_lm(pool%sht(1),0,0),k) = &
                    cmplx(real(tval,dp)*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
        end do
        now%spec%lnps(aeros_sht_lm(pool%sht(1),0,0)) = &
                cmplx(log(real(p0,dp))*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
        now%spec%temp(aeros_sht_lm(pool%sht(1),2,0),:) = cmplx(-5.0_dp, 0.0_dp, wp_sh)

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)

        do j = 1, grd%nlat
            do i = 1, grd%nlon
                call aeros_vgrid_pressure(vg, now%ps(i,j), phalf, pfull, dpc)
                do k = 1, nlev
                    call aeros_qsat(now%temp_g(i,j,k), pfull(k), qs, dqsdt)
                    now%qv_g(i,j,k) = 0.9_wp*qs
                end do
            end do
        end do

        t1_0 = gmean_lowest_t(now)
        precip_tot = 0.0_dp
        nan_seen = .FALSE.; neg_seen = .FALSE.

        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)
            do j = 1, grd%nlat
                do i = 1, grd%nlon
                    precip_tot = precip_tot &
                        + real(ts%cnv%precip(i,j) + ts%cnd%precip(i,j),dp) &
                                    *real(grd%area(i,j),dp)*real(par%dt,dp)
                end do
            end do
            if (any(now%qv_g < 0.0_wp)) neg_seen = .TRUE.
            if (any(now%qv_g /= now%qv_g)) nan_seen = .TRUE.
        end do

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        if (any(now%temp_g /= now%temp_g)) nan_seen = .TRUE.
        umax = real(maxval(abs(now%u)), dp)
        t1_1 = gmean_lowest_t(now)

        ! Global-mean radiation and surface-flux diagnostics from the last step.
        olr_gm   = gmean(ts%rad%olr)
        swup_gm  = gmean(ts%rad%sw_up_toa)
        shf_gm   = gmean(ts%surf%shf)
        lhf_gm   = gmean(ts%surf%lhf)
        swtoa_gm = 0.0_dp
        do j = 1, grd%nlat
            swtoa_gm = swtoa_gm + real(ts%rad%sw_toa(j),dp)*sum(real(grd%area(:,j),dp))
        end do
        swtoa_gm = swtoa_gm/real(sum(grd%area),dp)
        toa_net  = swtoa_gm - swup_gm - olr_gm

        write(*,*) ""
        write(*,"(a40,f10.2,a)")  "   global-mean insolation (TOA)   ", swtoa_gm, " W/m2"
        write(*,"(a40,f10.2,a)")  "   global-mean OLR                ", olr_gm, " W/m2"
        write(*,"(a40,f10.2,a)")  "   global-mean reflected SW       ", swup_gm, " W/m2"
        write(*,"(a40,f10.2,a)")  "   TOA net imbalance              ", toa_net, " W/m2"
        write(*,"(a40,f10.2,a)")  "   global-mean sensible flux      ", shf_gm, " W/m2"
        write(*,"(a40,f10.2,a)")  "   global-mean latent flux        ", lhf_gm, " W/m2"
        write(*,"(a40,f10.2,a)")  "   global-mean lowest-layer T t=0 ", t1_0, " K"
        write(*,"(a40,f10.2,a)")  "   global-mean lowest-layer T end ", t1_1, " K"
        write(*,"(a40,es12.3,a)") "   total precipitated             ", precip_tot, " kg"
        write(*,"(a40,f10.2,a)")  "   max |u| at the end             ", umax, " m/s"

        call check(.not. nan_seen, "the RCE stack runs with no NaN in T or q", nfail)
        call check(.not. neg_seen, "q stays non-negative through the run", nfail)
        call check(umax < 200.0_dp, "the wind stays bounded (no blow-up)", nfail)

        call check(olr_gm > 150.0_dp .and. olr_gm < 320.0_dp, &
                   "global-mean OLR is physical (radiation is active)", nfail)
        call check(swtoa_gm > 320.0_dp .and. swtoa_gm < 350.0_dp, &
                   "global-mean insolation is ~ S0/4", nfail)
        call check(abs(shf_gm) > 0.0_dp .and. abs(lhf_gm) > 0.0_dp, &
                   "surface fluxes are active", nfail)
        call check(precip_tot > 0.0_dp, "it rains (the moist loop is closed)", nfail)

        call check(t1_1 > 200.0_dp .and. t1_1 < 330.0_dp, &
                   "lowest-layer temperature stays in a physical band", nfail)
        call check(abs(toa_net) < 250.0_dp, &
                   "the TOA imbalance is finite (not running away)", nfail)

        call aeros_timestep_end(ts)
        call aeros_state_end(now)
        return
    end subroutine run

    real(dp) function gmean(f) result(m)
        ! area-weighted global mean of a (nlon,nlat) field
        implicit none
        real(wp), intent(in) :: f(:,:)
        real(dp) :: s
        integer  :: i, j
        s = 0.0_dp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                s = s + real(f(i,j),dp)*real(grd%area(i,j),dp)
            end do
        end do
        m = s/real(sum(grd%area),dp)
    end function gmean

    real(dp) function gmean_lowest_t(st) result(m)
        implicit none
        type(aeros_state_class), intent(in) :: st
        m = gmean(st%temp_g(:,:,nlev))
    end function gmean_lowest_t

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

end program test_rce
