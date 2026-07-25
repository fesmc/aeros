program test_moist_run
    ! Coupled integration smoke test: dynamics + humidity transport +
    ! convection + condensation, stepping together.
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
    !   STABILITY. 200 steps with no NaN and a bounded wind. The convective
    !   heating is the largest, most sharply vertically structured new term in
    !   the temperature equation; routed through the centered leapfrog it would
    !   excite the computational mode and NaN within tens of steps. It is applied
    !   forward-split instead (aeros_timestep step 6), and staying stable here is
    !   the assertion that that path is right. Condensation's smaller heating
    !   still rides the centered dtdt path.
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
        real(dp) :: w0, w1, precip_tot, conv_precip_tot, closure, umax
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

        ! Turn moist convection on too: Simplified Betts-Miller with its default
        ! tau and rh_ref (set by aeros_convection_init). This is the coupled test
        ! of the forward-split heating path -- convective heating is the largest,
        ! most sharply structured new term, and if it rode the centered leapfrog
        ! the run would go computational-mode unstable and NaN within tens of
        ! steps. Staying stable for 200 steps is the assertion that it does not.
        ts%cnv%enabled = .TRUE.

        allocate(phis2(grd%nlon, grd%nlat))
        phis2 = 0.0_wp
        call aeros_timestep_set_phis(ts, phis2)

        ! A resting atmosphere with a warm, humid, CONDITIONALLY UNSTABLE profile
        ! -- surface ~298 K, a lapse steeper than moist-adiabatic aloft -- so the
        ! columns actually convect and the coupled convective path is exercised,
        ! not just present. Lightly perturbed so Held-Suarez has an eddy to grow.
        ! The shape matters: a warm stratosphere would put huge q_sat at a few
        ! hundred pascals (the denominator p - (1-eps)e_s collapses) and seeding
        ! humidity against it would be absurd; this profile stays cold aloft
        ! (~223 K at the top layer). T(k) in sigma ~ (k-0.5)/nlev.
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

        ! Seed humidity at 80% RH against the diagnosed T and p. With a realistic
        ! profile the upper-level q_sat is small, so this is a physical moisture
        ! field (surface ~10 g/kg, tending to zero aloft), not the runaway a
        ! warm stratosphere would give.
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                call aeros_vgrid_pressure(vg, now%ps(i,j), phalf, pfull, dpc)
                do k = 1, nlev
                    call aeros_qsat(now%temp_g(i,j,k), pfull(k), qs, dqsdt)
                    now%qv_g(i,j,k) = 0.9_wp*qs
                end do
            end do
        end do

        w0 = water(now)
        precip_tot = 0.0_dp; conv_precip_tot = 0.0_dp
        nan_seen = .FALSE.; neg_seen = .FALSE.; umax = 0.0_dp

        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)

            ! Accumulate precipitated mass: sum rate*area*dt over the grid, from
            ! BOTH sinks -- convection rains first, condensation mops up the rest.
            do j = 1, grd%nlat
                do i = 1, grd%nlon
                    precip_tot = precip_tot &
                        + real(ts%cnv%precip(i,j) + ts%cnd%precip(i,j),dp) &
                                    *real(grd%area(i,j),dp)*real(par%dt,dp)
                    conv_precip_tot = conv_precip_tot + real(ts%cnv%precip(i,j),dp) &
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
        write(*,"(a40,es12.3,a)") "     of which convective          ", conv_precip_tot, " kg"
        write(*,"(a40,es12.3)")   "   fraction of vapour rained out  ", (w0-w1)/w0
        write(*,"(a40,es12.3)")   "   water budget closure |resid|/W ", closure
        write(*,"(a40,f10.2,a)")  "   max |u| at the end             ", umax, " m/s"

        call check(.not. nan_seen, "200 coupled steps with no NaN in T or q", nfail)
        call check(.not. neg_seen, "q stays non-negative through the coupled run", nfail)
        call check(umax < 200.0_dp, "the wind stays bounded (no blow-up)", nfail)
        call check(precip_tot > 0.0_dp, "it actually rains (the closure test is not vacuous)", nfail)
        call check(conv_precip_tot > 0.0_dp, "convection is active in the coupled run", nfail)
        ! ~1e-3 at T21L12 over 200 steps with convection and condensation both
        ! on. This is the FV-vs-spectral AIR-mass gap (m2_results.md §8), not
        ! tracer-transport error: the transport conserves its own water to
        ! machine precision, but the diagnostic weights q by the spectral layer
        ! masses, which differ from the finite-volume ones by O(truncation). It
        ! scales with the roughness of q -- the moist sinks sharpen q and raise
        ! it (from ~2e-4 on the smooth advected field), and so, counter-
        ! intuitively, does the van Leer limiter, because a less-diffused q has
        ! larger gradients for the gap to act on. It shrinks with resolution, not
        ! with the limiter. The bound catches gross non-conservation; the number
        ! itself is a measurement.
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
