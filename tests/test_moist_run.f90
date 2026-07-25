program test_moist_run
    ! Coupled integration smoke test: dynamics + humidity transport +
    ! condensation, stepping together.
    !
    ! The two operator unit tests (test_moisture, test_condensation) check the
    ! transport and the condensation in isolation and to machine precision. What
    ! they cannot check is that the pieces work TOGETHER through a real timestep:
    ! that the latent heating computed on the grid actually reaches the spectral
    ! temperature and the run stays stable, that transport moving humidity into
    ! a cold column makes it rain, and that the whole thing conserves water at
    ! the level the design permits. That is what this does -- a short moist
    ! Held-Suarez integration with everything on.
    !
    ! What is asserted:
    !
    !   STABILITY. 200 steps with no NaN and a bounded wind. The latent heating
    !   is the largest new term in the temperature equation; if its magnitude or
    !   its sign or its path into the spectral state were wrong, the run would
    !   blow up here.
    !
    !   POSITIVITY SURVIVES COUPLING. q >= 0 at every step of the full model, not
    !   just under prescribed winds -- the dynamical winds, the transport and the
    !   condensation sink all act on q, and it must stay non-negative through all
    !   of them.
    !
    !   IT ACTUALLY RAINS. Cumulative precipitation is positive: the coupling is
    !   doing something, so the conservation check below is not vacuous.
    !
    !   WATER CLOSES, to O(truncation). Initial vapour = final vapour + total
    !   precipitation, to the truncation-level consistency the finite-volume
    !   transport carries against the spectral surface pressure (m2_results.md
    !   §8) -- NOT machine precision, and the point is to measure how large that
    !   gap actually is on a running model rather than to assert it away.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, p0, R_d, grav, &
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
    integer, parameter :: nstep = 200

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)     :: grd
    type(aeros_vgrid_class)    :: vg
    type(aeros_param_class)    :: par
    type(aeros_state_class)    :: now
    type(aeros_timestep_class) :: ts

    integer  :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_moist_run:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    call run(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_moist_run:: PASS"
    else
        write(*,*) "test_moist_run:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine run(nfail)
        implicit none
        integer, intent(inout) :: nfail

        real(wp), allocatable :: phis2(:,:)
        real(wp) :: qs, dqsdt, tval
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev)
        real(dp) :: w0, w1, precip_tot, closure, umax
        logical  :: nan_seen, neg_seen
        integer  :: i, j, k, n

        par%trunc = trunc; par%nlon = -1; par%nlat = -1; par%nlev = nlev
        par%nthreads = -1
        par%dt = 1800.0_wp
        par%semi_implicit = .TRUE.
        par%held_suarez   = .TRUE.       ! motion, hence cooling, hence rain
        par%eps_filter = 0.06_wp
        par%raw_alpha  = 0.53_wp
        par%ndiff = 6; par%tau_diff = 6.0_wp
        par%mass_fixer = .FALSE.

        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        ! Turn condensation on directly (no namelist in this test), rh_crit = 1.
        ts%cnd%enabled = .TRUE.
        ts%cnd%rh_crit = 1.0_wp

        allocate(phis2(grd%nlon, grd%nlat))
        phis2 = 0.0_wp
        call aeros_timestep_set_phis(ts, phis2)

        ! A resting atmosphere with a REALISTIC lapse-rate profile -- warm near
        ! the surface, cold aloft -- lightly perturbed so Held-Suarez has an
        ! eddy to grow. The profile matters: an isothermal-warm column would put
        ! a 288 K stratosphere at a few hundred pascals, where q_sat is enormous
        ! (the denominator p - (1-eps)e_s collapses), and seeding against it
        ! would inject absurd humidity aloft. T(k) from a linear sigma profile,
        ! ~219 K at the top layer to ~285 K at the surface.
        call aeros_spec_zero(now%spec)
        do k = 1, nlev
            tval = 216.0_wp + 72.0_wp*(real(k,wp) - 0.5_wp)/real(nlev,wp)
            now%spec%temp(aeros_sht_lm(pool%sht(1),0,0),k) = &
                    cmplx(real(tval,dp)*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
        end do
        now%spec%lnps(aeros_sht_lm(pool%sht(1),0,0)) = &
                cmplx(log(real(p0,dp))*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
        now%spec%temp(aeros_sht_lm(pool%sht(1),2,0),:) = cmplx(-5.0_dp, 0.0_dp, wp_sh)

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)

        ! Seed humidity at 80% RH against the diagnosed T and p. With a realistic
        ! profile the upper-level q_sat is small, so this is a physical moisture
        ! field (surface ~10 g/kg, tending to zero aloft), not the runaway a
        ! warm stratosphere would give.
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                call aeros_vgrid_pressure(vg, now%ps(i,j), phalf, pfull, dpc)
                do k = 1, nlev
                    call aeros_qsat(now%temp_g(i,j,k), pfull(k), qs, dqsdt)
                    now%qv_g(i,j,k) = 0.8_wp*qs
                end do
            end do
        end do

        w0 = water(now)
        precip_tot = 0.0_dp
        nan_seen = .FALSE.; neg_seen = .FALSE.; umax = 0.0_dp

        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)

            ! Accumulate precipitated mass: sum rate*area*dt over the grid.
            do j = 1, grd%nlat
                do i = 1, grd%nlon
                    precip_tot = precip_tot + real(ts%cnd%precip(i,j),dp) &
                                    *real(grd%area(i,j),dp)*real(par%dt,dp)
                end do
            end do

            if (any(now%qv_g < 0.0_wp)) neg_seen = .TRUE.
            if (any(now%qv_g /= now%qv_g)) nan_seen = .TRUE.
        end do

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        if (any(now%temp_g /= now%temp_g)) nan_seen = .TRUE.
        umax = real(maxval(abs(now%u)), dp)
        w1 = water(now)

        ! Water closure: what left the vapour reservoir should be what rained.
        closure = abs((w0 - w1) - precip_tot)/w0

        write(*,*) ""
        write(*,"(a40,es12.3,a)") "   initial column water           ", w0, " kg"
        write(*,"(a40,es12.3,a)") "   total precipitated             ", precip_tot, " kg"
        write(*,"(a40,es12.3)")   "   fraction of vapour rained out  ", (w0-w1)/w0
        write(*,"(a40,es12.3)")   "   water budget closure |resid|/W ", closure
        write(*,"(a40,f10.2,a)")  "   max |u| at the end             ", umax, " m/s"

        call check(.not. nan_seen, "200 coupled steps with no NaN in T or q", nfail)
        call check(.not. neg_seen, "q stays non-negative through the coupled run", nfail)
        call check(umax < 200.0_dp, "the wind stays bounded (no blow-up)", nfail)
        call check(precip_tot > 0.0_dp, "it actually rains (the closure test is not vacuous)", nfail)
        ! ~2e-3 at T21L12 over 200 steps, dominated by the finite-volume
        ! transport's O(truncation) dispersion -- which condensation enlarges by
        ! sharpening the humidity field (it was ~2e-4 on the smooth advected
        ! field, condensation off). Not machine precision, not claimed to be:
        ! the bound catches gross non-conservation, and the van Leer limiter
        ! (the accuracy commit) is what shrinks the number.
        call check(closure < 5.0e-3_dp, &
                    "water closes to the transport's truncation-level consistency", nfail)

        call aeros_timestep_end(ts)
        call aeros_state_end(now)
        return
    end subroutine run

    real(dp) function water(st) result(w)
        implicit none
        type(aeros_state_class), intent(in) :: st
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev)
        real(dp) :: col, tot
        integer  :: i, j, k
        tot = 0.0_dp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                call aeros_vgrid_pressure(vg, st%ps(i,j), phalf, pfull, dpc)
                col = 0.0_dp
                do k = 1, nlev
                    col = col + real(st%qv_g(i,j,k),dp)*real(dpc(k),dp)
                end do
                tot = tot + col*real(grd%area(i,j),dp)/real(grav,dp)
            end do
        end do
        w = tot
    end function water

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

end program test_moist_run
