program test_held_suarez
    ! Acceptance test for the Held & Suarez (1994) forcing.
    !
    ! The benchmark's own validation is a 1200-day integration compared against
    ! published FIGURES -- there are no target numbers in the literature, only
    ! "check that yours looks like ours" -- so that lives in
    ! drivers/held_suarez.f90 and docs/m1_results.md. What can be asserted in
    ! seconds, and is asserted here, is that the forcing FUNCTIONS are the ones
    ! the paper specifies and that they are wired into the model correctly.
    !
    ! 1. THE FORMULAS, against values computed by hand from the paper. If T_eq
    !    is wrong the 1200-day run still produces a plausible-looking
    !    circulation, just not Held-Suarez' -- which is the failure mode a
    !    figure comparison is worst at catching, so it is worth pinning
    !    pointwise.
    !
    ! 2. THE WIRING. The forcing has to reach the model's tendencies with the
    !    right sign and the right magnitude, and -- the part that is easy to get
    !    wrong -- the Rayleigh drag has to be applied on the GRID before the
    !    vector analysis, because k_v varies horizontally and does not commute
    !    with the curl. So the drag is checked against a state whose vorticity
    !    is confined to the free troposphere: k_v is zero there, and a spectral
    !    application of a horizontally averaged k_v would damp it anyway.
    !
    ! 3. RELAXATION TOWARD EQUILIBRIUM. Started AT T_eq with no wind, the
    !    thermal forcing must vanish. Started away from it, the temperature has
    !    to move toward it at the advertised rate.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, kappa, p0, R_d, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_spec_class, aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_tendency
    use aeros_timestep
    use aeros_held_suarez

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 20

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg

    integer :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))

    ! The benchmark's own vertical grid: 20 evenly spaced sigma levels. p_top
    ! is 1000 Pa rather than 0 -- Held & Suarez define sigma on [0,1], but a
    ! zero-pressure top makes the top layer hydrostatically inconsistent
    ! (Simmons & Burridge's alpha_1 = ln 2 convention; see tests/test_tendency).
    call aeros_vgrid_init(vg, nlev, stretch_a=1.0_wp, sigma_t=0.0_wp, &
                            p_top=1000.0_wp, ps_ref=1.0e5_wp, t_ref=300.0_wp)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_held_suarez:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    call test_formulas(nfail)
    call test_wiring(nfail)
    call test_relaxation(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_held_suarez:: PASS"
    else
        write(*,*) "test_held_suarez:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    ! === 1. The specification ================================================

    subroutine test_formulas(nfail)

        implicit none

        integer, intent(inout) :: nfail

        real(wp) :: teq, expect, s2, c2, c4

        write(*,*) ""
        write(*,*) " -- Held & Suarez (1994) eqs. (1)-(3), pointwise"

        ! --- T_eq at the EQUATOR, surface (p = p0). sin^2 = 0, cos^2 = 1,
        ! log(p/p0) = 0, (p/p0)^kappa = 1, so T_eq = 315 K exactly.
        teq = aeros_hs_teq(1.0e5_wp, 0.0_wp, 1.0_wp)
        write(*,"(a40,f10.4,a)") "   T_eq(equator, 1000 hPa) ", teq, " K"
        call check(abs(teq - 315.0_wp) < 1.0e-12_wp, &
                    "T_eq is 315 K at the equatorial surface", nfail)

        ! --- T_eq at the POLE, surface. sin^2 = 1, cos^2 = 0, so
        ! T_eq = 315 - 60 = 255 K, above the 200 K floor.
        teq = aeros_hs_teq(1.0e5_wp, 1.0_wp, 0.0_wp)
        write(*,"(a40,f10.4,a)") "   T_eq(pole, 1000 hPa)    ", teq, " K"
        call check(abs(teq - 255.0_wp) < 1.0e-12_wp, &
                    "T_eq is 315 - dT_y at the polar surface", nfail)

        ! --- T_eq at 45 N, 500 hPa. Every term active; computed here from the
        ! paper's expression written out independently of the module.
        s2 = 0.5_wp; c2 = 0.5_wp
        expect = (315.0_wp - 60.0_wp*s2 - 10.0_wp*log(0.5_wp)*c2)*0.5_wp**kappa
        teq    = aeros_hs_teq(0.5e5_wp, s2, c2)
        write(*,"(a40,f10.4,a)") "   T_eq(45N, 500 hPa)      ", teq, " K"
        write(*,"(a40,f10.4,a)") "   hand-computed           ", expect, " K"
        call check(abs(teq - expect) < 1.0e-12_wp, &
                    "T_eq matches the written-out expression", nfail)

        ! --- The 200 K floor must bind in the stratosphere. At 10 hPa the
        ! bracket is ~330 K but (p/p0)^kappa is ~0.27, so the product is well
        ! under 200 and the floor takes over.
        teq = aeros_hs_teq(1.0e3_wp, 0.0_wp, 1.0_wp)
        write(*,"(a40,f10.4,a)") "   T_eq(equator, 10 hPa)   ", teq, " K"
        call check(abs(teq - 200.0_wp) < 1.0e-12_wp, &
                    "the 200 K stratospheric floor binds aloft", nfail)

        ! --- k_T: k_a above sigma_b everywhere, k_s at the tropical surface.
        write(*,"(a40,f10.4,a)") "   1/k_T(equator, sigma=1) ", &
                                    1.0_wp/(aeros_hs_kt(1.0_wp, 1.0_wp)*86400.0_wp), " day"
        write(*,"(a40,f10.4,a)") "   1/k_T(any, sigma=0.5)   ", &
                                    1.0_wp/(aeros_hs_kt(0.5_wp, 1.0_wp)*86400.0_wp), " day"
        write(*,"(a40,f10.4,a)") "   1/k_T(pole, sigma=1)    ", &
                                    1.0_wp/(aeros_hs_kt(1.0_wp, 0.0_wp)*86400.0_wp), " day"

        call check(abs(aeros_hs_kt(1.0_wp, 1.0_wp) - hs_k_s) < 1.0e-18_wp, &
                    "k_T is k_s at the tropical surface", nfail)
        call check(abs(aeros_hs_kt(0.5_wp, 1.0_wp) - hs_k_a) < 1.0e-18_wp, &
                    "k_T is k_a above sigma_b", nfail)
        call check(abs(aeros_hs_kt(1.0_wp, 0.0_wp) - hs_k_a) < 1.0e-18_wp, &
                    "k_T is k_a at the poles even at the surface", nfail)

        ! --- k_v: zero above sigma_b, k_f at the surface, linear between. At
        ! sigma = 0.85, halfway from 0.7 to 1, it must be exactly k_f/2.
        write(*,"(a40,f10.4,a)") "   1/k_v(sigma=1)          ", &
                                    1.0_wp/(aeros_hs_kv(1.0_wp)*86400.0_wp), " day"
        call check(aeros_hs_kv(0.5_wp) == 0.0_wp, &
                    "k_v vanishes above sigma_b", nfail)
        call check(aeros_hs_kv(0.7_wp) == 0.0_wp, &
                    "k_v vanishes exactly at sigma_b", nfail)
        call check(abs(aeros_hs_kv(1.0_wp) - hs_k_f) < 1.0e-18_wp, &
                    "k_v is k_f at the surface", nfail)
        call check(abs(aeros_hs_kv(0.85_wp) - 0.5_wp*hs_k_f) < 1.0e-18_wp, &
                    "k_v is linear in sigma below sigma_b", nfail)

        ! --- The published constants, so a typo in a parameter is caught here
        ! rather than 1200 days later.
        call check(hs_t_equator == 315.0_wp .and. hs_dt_y == 60.0_wp &
                    .and. hs_dtheta_z == 10.0_wp .and. hs_t_min == 200.0_wp &
                    .and. hs_sigma_b == 0.7_wp, &
                    "the published constants are unchanged", nfail)
        call check(abs(1.0_wp/(hs_k_a*86400.0_wp) - 40.0_wp) < 1.0e-12_wp &
                    .and. abs(1.0_wp/(hs_k_s*86400.0_wp) - 4.0_wp) < 1.0e-12_wp &
                    .and. abs(1.0_wp/(hs_k_f*86400.0_wp) - 1.0_wp) < 1.0e-12_wp, &
                    "the published damping timescales are 40, 4 and 1 days", nfail)

        return

    end subroutine test_formulas

    ! === 2. The wiring =======================================================

    subroutine test_wiring(nfail)
        ! Does the forcing reach wrk%dtdt, wrk%ae and wrk%an, with the right
        ! sign and size?

        implicit none

        integer, intent(inout) :: nfail

        type(aeros_sht_class), pointer :: s
        type(aeros_spec_class) :: spec
        type(aeros_work_class) :: wrk
        type(aeros_hs_class)   :: hs

        real(wp), allocatable :: dtdt0(:,:,:), ae0(:,:,:), an0(:,:,:)
        real(wp) :: ps, sig, teq, kt, kv, expect
        real(wp) :: p_half(0:nlev), p_full(nlev)
        real(dp) :: y00, err_t, err_u
        integer  :: i, j, k

        s => pool%sht(1)
        y00 = sqrt(4.0_dp*acos(-1.0_dp))

        write(*,*) ""
        write(*,*) " -- the forcing reaches the grid-space right-hand sides"

        call aeros_spec_alloc(spec, s%nlm, nlev)
        call aeros_work_alloc(wrk, grd%nlon, grd%nlat, nlev)
        call aeros_hs_init(hs, grd, .TRUE.)

        allocate(dtdt0(grd%nlon,grd%nlat,nlev))
        allocate(ae0(grd%nlon,grd%nlat,nlev), an0(grd%nlon,grd%nlat,nlev))

        ! An isothermal atmosphere at 300 K with a solid-body-ish rotation, so
        ! T - T_eq and the winds are both non-zero everywhere the forcing acts.
        call aeros_spec_zero(spec)
        do k = 1, nlev
            spec%temp(aeros_sht_lm(s,0,0),k) = cmplx(300.0_dp*y00, 0.0_dp, wp_sh)
            spec%vor(aeros_sht_lm(s,1,0),k)  = cmplx(2.0e-5_dp, 0.0_dp, wp_sh)
        end do
        spec%lnps(aeros_sht_lm(s,0,0)) = cmplx(log(1.0e5_dp)*y00, 0.0_dp, wp_sh)

        wrk%phis = 0.0_wp
        call aeros_tendency_grid(pool, vg, grd, spec, wrk)

        dtdt0 = wrk%dtdt; ae0 = wrk%ae; an0 = wrk%an

        call aeros_hs_apply(hs, vg, wrk)

        ! Recompute the expected increment independently, column by column.
        err_t = 0.0_dp; err_u = 0.0_dp
        do j = 1, grd%nlat
            do i = 1, grd%nlon

                ps = exp(wrk%lnps_g(i,j))
                call aeros_vgrid_pressure(vg, ps, p_half, p_full)

                do k = 1, nlev
                    sig = p_full(k)/ps
                    teq = aeros_hs_teq(p_full(k), &
                                real(grd%sinlat(j)**2, wp), 1.0_wp - real(grd%sinlat(j)**2, wp))
                    kt  = aeros_hs_kt(sig, (1.0_wp - real(grd%sinlat(j)**2, wp))**2)
                    kv  = aeros_hs_kv(sig)

                    expect = -kt*(wrk%t_g(i,j,k) - teq)
                    err_t  = max(err_t, abs(real(wrk%dtdt(i,j,k) - dtdt0(i,j,k) - expect, dp)))

                    expect = -kv*wrk%u(i,j,k)
                    err_u  = max(err_u, abs(real(wrk%ae(i,j,k) - ae0(i,j,k) - expect, dp)))
                end do

            end do
        end do

        write(*,"(a40,es12.3,a)") "   dT/dt increment error ", err_t, " K s-1"
        write(*,"(a40,es12.3,a)") "   drag increment error  ", err_u, " m s-2"
        call check(err_t < 1.0e-18_dp, "the thermal forcing is applied as specified", nfail)
        call check(err_u < 1.0e-18_dp, "the Rayleigh drag is applied as specified", nfail)

        ! The forcing must actually DO something -- an increment of zero would
        ! pass the two checks above.
        write(*,"(a40,es12.3,a)") "   max |dT/dt| from forcing ", &
                                    maxval(abs(wrk%dtdt - dtdt0)), " K s-1"
        write(*,"(a40,es12.3,a)") "   max |drag|               ", &
                                    maxval(abs(wrk%ae - ae0)), " m s-2"
        call check(maxval(abs(wrk%dtdt - dtdt0)) > 1.0e-7_wp, &
                    "the thermal forcing is non-trivial", nfail)
        call check(maxval(abs(wrk%ae - ae0)) > 1.0e-7_wp, &
                    "the drag is non-trivial", nfail)

        ! The drag must be CONFINED below sigma_b. With 20 evenly spaced sigma
        ! levels and p_top = 1000 Pa, sigma_full(k) < 0.7 for k <= 13, so the
        ! top 13 levels must be untouched. A spectral application of the drag
        ! -- damping zeta and D by a horizontally averaged k_v -- would leak
        ! into every level, which is what this catches.
        call check(maxval(abs(wrk%ae(:,:,1:13) - ae0(:,:,1:13))) == 0.0_wp .and. &
                    maxval(abs(wrk%an(:,:,1:13) - an0(:,:,1:13))) == 0.0_wp, &
                    "the drag is exactly zero above sigma_b", nfail)

        deallocate(dtdt0, ae0, an0)
        call aeros_hs_end(hs)
        call aeros_work_end(wrk)
        call aeros_spec_end(spec)

        return

    end subroutine test_wiring

    ! === 3. Relaxation =======================================================

    subroutine test_relaxation(nfail)
        ! Integrate, and check the temperature moves toward T_eq.

        implicit none

        integer, intent(inout) :: nfail

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts
        type(aeros_sht_class), pointer :: s
        real(wp) :: phis(grd%nlon,grd%nlat)
        real(dp) :: y00, d0, d1, umax
        integer  :: n, k
        integer, parameter :: nstep = 200      ! ~4 days at dt = 1800 s

        s => pool%sht(1)
        y00 = sqrt(4.0_dp*acos(-1.0_dp))

        write(*,*) ""
        write(*,*) " -- relaxation toward equilibrium, 200 steps"

        par%trunc = trunc; par%nlon = -1; par%nlat = -1; par%nlev = nlev
        par%nthreads = -1
        par%dt = 1800.0_wp
        par%semi_implicit = .TRUE.
        par%held_suarez   = .TRUE.
        par%eps_filter = 0.06_wp
        par%raw_alpha  = 0.53_wp
        par%ndiff      = 6
        par%tau_diff   = 6.0_wp

        call aeros_state_alloc(now, grd, s%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        ! Isothermal at 250 K, at rest -- a long way from T_eq at both the
        ! equator (315 K) and the pole (255 K), so the relaxation has to work
        ! in both directions.
        call aeros_spec_zero(now%spec)
        do k = 1, nlev
            now%spec%temp(aeros_sht_lm(s,0,0),k) = cmplx(250.0_dp*y00, 0.0_dp, wp_sh)
        end do
        now%spec%lnps(aeros_sht_lm(s,0,0)) = cmplx(log(1.0e5_dp)*y00, 0.0_dp, wp_sh)

        phis = 0.0_wp
        call aeros_timestep_set_phis(ts, phis)

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        d0 = teq_distance(now)

        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, now)
        end do

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        d1 = teq_distance(now)

        umax = maxval(abs(real(now%u, dp)))

        write(*,"(a40,f10.4,a)")  "   rms(T - T_eq), initial ", d0, " K"
        write(*,"(a40,f10.4,a)")  "   rms(T - T_eq), 4 days  ", d1, " K"
        write(*,"(a40,f10.4)")    "   fraction remaining     ", d1/d0
        write(*,"(a40,f10.4,a)")  "   max |u| after 4 days   ", umax, " m s-1"

        call check(d1 < d0, "the temperature moves toward T_eq", nfail)

        ! Roughly: the free troposphere relaxes at k_a = 1/40 day, so over 4
        ! days about 10% of the departure should be gone. Bounded on both sides
        ! -- too little means the forcing is not reaching the model, too much
        ! means it is being applied with the wrong rate.
        call check(d1/d0 > 0.7_dp .and. d1/d0 < 0.97_dp, &
                    "it relaxes at roughly the k_a rate", nfail)

        ! The forcing generates a circulation, which is the whole point.
        call check(umax > 1.0_dp .and. umax < 150.0_dp, &
                    "a plausible circulation develops", nfail)

        call aeros_timestep_end(ts)
        call aeros_state_end(now)

        return

    end subroutine test_relaxation

    real(dp) function teq_distance(now) result(d)
        ! Area-weighted rms departure of T from T_eq [K].

        implicit none

        type(aeros_state_class), intent(in) :: now

        real(wp) :: p_half(0:nlev), p_full(nlev)
        real(wp) :: ps, teq, s2
        real(dp) :: acc, wsum, w
        integer  :: i, j, k

        acc = 0.0_dp; wsum = 0.0_dp

        do j = 1, grd%nlat
            s2 = real(grd%sinlat(j)**2, wp)
            do i = 1, grd%nlon
                ps = now%ps(i,j)
                call aeros_vgrid_pressure(vg, ps, p_half, p_full)
                do k = 1, nlev
                    teq  = aeros_hs_teq(p_full(k), s2, 1.0_wp - s2)
                    w    = real(grd%area(i,j), dp)*real(p_half(k) - p_half(k-1), dp)
                    acc  = acc + w*real(now%temp_g(i,j,k) - teq, dp)**2
                    wsum = wsum + w
                end do
            end do
        end do

        d = sqrt(acc/wsum)

        return

    end function teq_distance

    subroutine check(ok, label, nfail)

        implicit none

        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail

        if (ok) then
            write(*,*) "  ok   : ", trim(label)
        else
            write(*,*) "  FAIL : ", trim(label)
            nfail = nfail + 1
        end if

        return

    end subroutine check

end program test_held_suarez
