program rce_long
    ! Long coupled radiative-convective integration, instrumented, for chasing
    ! the lowest-layer hot-spot instability (m2_handoff.md task #7). This is the
    ! scratch driver test_rce.f90 refers to but does not commit: same physics
    ! stack, but run for model months with a per-level T trajectory so the blow
    ! up can be localized in level, space and time, and with every relevant knob
    ! exposed through a namelist so the hypothesis matrix can be run without
    ! recompiling.
    !
    ! Namelist group "rce" (file = arg 1, default rce.nml; all optional):
    !   trunc, nlev, nstep, dt, tau_diff, ndiff, eps_filter, raw_alpha
    !   l_surf, l_cnv, l_cnd, l_rad, l_sponge      physics toggles
    !   conv_tau, c_h, c_e, u_min                  scheme knobs
    !   print_every                                trajectory cadence [steps]
    !
    ! Each report prints, per level: global-mean T, the global max T with its
    ! (i,j) location, and the global min T -- the max-location wandering and a
    ! single level running away is the hot-spot signature. It also prints the
    ! per-level global-mean forward-split heating (surface+convection+radiation,
    ! wrk%dt_phys) and centered heating (condensation, wrk%dtdt) as K/day, so the
    ! term driving a level is visible directly. Stops at the first NaN, naming
    ! the level and step.

    use aeros_defs,     only : dp, wp, wp_sh, p0, R_d, grav, S0, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_timestep
    use nml,            only : nml_read

    implicit none

    ! --- namelist-configurable, with test_rce defaults -----------------------
    character(len=512) :: nmlfile
    integer  :: trunc = 21, nlev = 12, nstep = 4800, ndiff = 6
    integer  :: print_every = 96          ! ~2 model days at dt=1800
    real(wp) :: dt = 1800.0_wp, tau_diff = 6.0_wp
    real(wp) :: eps_filter = 0.06_wp, raw_alpha = 0.53_wp
    real(wp) :: conv_tau = 7200.0_wp, c_h = 1.5e-3_wp, c_e = 1.5e-3_wp, u_min = 1.0_wp
    logical  :: l_surf = .TRUE., l_cnv = .TRUE., l_cnd = .TRUE.
    logical  :: l_rad = .TRUE., l_sponge = .TRUE., l_vdiff = .FALSE.
    real(wp) :: vdiff_k0 = 10.0_wp, vdiff_sigma = 0.7_wp
    real(wp) :: seed_asym = 0.0_wp    ! zonal-asymmetry seed amplitude [K-ish]

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)     :: grd
    type(aeros_vgrid_class)    :: vg
    type(aeros_param_class)    :: par
    type(aeros_state_class)    :: now
    type(aeros_timestep_class) :: ts

    real(wp), allocatable :: phis2(:,:)
    real(wp) :: qs, dqsdt, tval
    real(wp) :: phalf(0:64), pfull(64), dpc(64)
    integer  :: i, j, k, n, m
    logical  :: blew_up

    nmlfile = "rce.nml"
    if (command_argument_count() >= 1) call get_command_argument(1, nmlfile)
    call read_config()

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " rce_long:: T", trunc, " L", nlev, &
                                       "  grid ", grd%nlon, "x", grd%nlat
    write(*,"(a,i0,a,f6.1,a,f5.1,a,i0)") "   nstep=", nstep, " dt=", dt, &
        " tau_diff=", tau_diff, " ndiff=", ndiff
    write(*,"(a,5(l1,1x))") "   surf/cnv/cnd/rad/sponge = ", &
        l_surf, l_cnv, l_cnd, l_rad, l_sponge

    par%trunc = trunc; par%nlon = -1; par%nlat = -1; par%nlev = nlev
    par%nthreads = -1
    par%dt = dt
    par%semi_implicit = .TRUE.
    par%held_suarez   = .FALSE.
    par%eps_filter = eps_filter
    par%raw_alpha  = raw_alpha
    par%ndiff = ndiff; par%tau_diff = tau_diff
    par%mass_fixer = .FALSE.

    call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
    call aeros_timestep_init(ts, par, pool, grd, vg)

    ts%cnd%enabled = l_cnd
    ts%cnd%rh_crit = 1.0_wp
    ts%cnv%enabled = l_cnv
    ts%cnv%tau = conv_tau
    ts%surf%enabled = l_surf
    ts%surf%c_h = c_h; ts%surf%c_e = c_e; ts%surf%u_min = u_min
    ts%rad%enabled  = l_rad
    ts%sponge_on    = l_sponge
    ts%vd%enabled   = l_vdiff
    ts%vd%k0        = vdiff_k0
    ts%vd%sigma     = vdiff_sigma

    allocate(phis2(grd%nlon, grd%nlat)); phis2 = 0.0_wp
    call aeros_timestep_set_phis(ts, phis2)

    ! Warm, humid, conditionally-unstable start (as test_rce/test_moist_run).
    call aeros_spec_zero(now%spec)
    do k = 1, nlev
        tval = 300.0_wp - 90.0_wp*(1.0_wp - ((real(k,wp) - 0.5_wp)/real(nlev,wp))**0.6_wp)
        now%spec%temp(aeros_sht_lm(pool%sht(1),0,0),k) = &
                cmplx(real(tval,dp)*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
    end do
    now%spec%lnps(aeros_sht_lm(pool%sht(1),0,0)) = &
            cmplx(log(real(p0,dp))*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
    now%spec%temp(aeros_sht_lm(pool%sht(1),2,0),:) = cmplx(-5.0_dp, 0.0_dp, wp_sh)

    ! Optional zonal-asymmetry seed: a small perturbation in a few m>0 modes so
    ! baroclinic instability can grow eddies. Without it the run is trapped in
    ! the axisymmetric (m=0) manifold -- with zonally symmetric forcing and an
    ! m=0-only start, the nonlinear terms never populate m>0, so there are no
    ! eddies to flux heat meridionally. Tests whether the residual subtropical
    ! warming is that artifact.
    if (seed_asym > 0.0_wp) then
        do k = 1, nlev
            do m = 1, 6
                now%spec%temp(aeros_sht_lm(pool%sht(1), m+2, m), k) = &
                    cmplx(real(seed_asym,dp), real(seed_asym,dp), wp_sh)
            end do
        end do
    end if

    call aeros_timestep_diagnose(ts, pool, vg, grd, now)
    do j = 1, grd%nlat
        do i = 1, grd%nlon
            call aeros_vgrid_pressure(vg, now%ps(i,j), phalf(0:nlev), pfull(1:nlev), dpc(1:nlev))
            do k = 1, nlev
                call aeros_qsat(now%temp_g(i,j,k), pfull(k), qs, dqsdt)
                now%qv_g(i,j,k) = 0.9_wp*qs
            end do
        end do
    end do

    call report(0)

    blew_up = .FALSE.
    do n = 1, nstep
        call aeros_timestep_step(ts, pool, vg, grd, now)
        if (any(ts%wrk%t_g /= ts%wrk%t_g) .or. any(now%qv_g /= now%qv_g)) then
            write(*,"(a,i0,a,f7.2,a)") " *** NaN at step ", n, "  (day ", &
                real(n,wp)*dt/86400.0_wp, ")"
            call locate_nan()
            blew_up = .TRUE.
            exit
        end if
        if (mod(n, print_every) == 0) call report(n)
    end do

    if (.not. blew_up) write(*,"(a)") " rce_long:: completed without NaN"

    call aeros_timestep_end(ts)
    call aeros_state_end(now)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

contains

    subroutine read_config()
        logical :: ex
        inquire(file=trim(nmlfile), exist=ex)
        if (.not. ex) then
            write(*,"(a)") " rce_long:: no "//trim(nmlfile)//", using defaults"
            return
        end if
        call nml_read(nmlfile, "rce", "trunc", trunc)
        call nml_read(nmlfile, "rce", "nlev", nlev)
        call nml_read(nmlfile, "rce", "nstep", nstep)
        call nml_read(nmlfile, "rce", "dt", dt)
        call nml_read(nmlfile, "rce", "tau_diff", tau_diff)
        call nml_read(nmlfile, "rce", "ndiff", ndiff)
        call nml_read(nmlfile, "rce", "eps_filter", eps_filter)
        call nml_read(nmlfile, "rce", "raw_alpha", raw_alpha)
        call nml_read(nmlfile, "rce", "conv_tau", conv_tau)
        call nml_read(nmlfile, "rce", "c_h", c_h)
        call nml_read(nmlfile, "rce", "c_e", c_e)
        call nml_read(nmlfile, "rce", "u_min", u_min)
        call nml_read(nmlfile, "rce", "print_every", print_every)
        call nml_read(nmlfile, "rce", "l_surf", l_surf)
        call nml_read(nmlfile, "rce", "l_cnv", l_cnv)
        call nml_read(nmlfile, "rce", "l_cnd", l_cnd)
        call nml_read(nmlfile, "rce", "l_rad", l_rad)
        call nml_read(nmlfile, "rce", "l_sponge", l_sponge)
        call nml_read(nmlfile, "rce", "l_vdiff", l_vdiff)
        call nml_read(nmlfile, "rce", "vdiff_k0", vdiff_k0)
        call nml_read(nmlfile, "rce", "vdiff_sigma", vdiff_sigma)
        call nml_read(nmlfile, "rce", "seed_asym", seed_asym)
        return
    end subroutine read_config

    subroutine report(n)
        integer, intent(in) :: n
        real(wp) :: tmx, tmn, tmean, umax
        real(wp) :: hphys, hcnd
        integer  :: imx, jmx, imn, jmn
        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        umax = maxval(abs(now%u))
        write(*,"(a)") ""
        write(*,"(a,i0,a,f8.2,a,f7.1,a)") " -- step ", n, "  day ", &
            real(n,wp)*dt/86400.0_wp, "   max|u| ", umax, " m/s"
        write(*,"(a)") "   lev    Tmean     Tmax  @(i,j)        Tmin    " // &
                       "Hfwd[K/d] Hcnd[K/d]"
        do k = 1, nlev
            call level_stats(k, tmean, tmx, imx, jmx, tmn, imn, jmn)
            ! forward-split (surf+cnv+rad) and centered (cnd) heating, K/day
            hphys = gmean_lev(ts%wrk%dt_phys, k)/dt*86400.0_wp
            hcnd  = gmean_lev(ts%wrk%dtdt, k)*86400.0_wp
            write(*,"(i6,f9.2,f9.2,a,i0,a,i0,a,f9.2,f10.3,f10.3)") &
                k, tmean, tmx, "  (", imx, ",", jmx, ")", tmn, hphys, hcnd
        end do
        call hotspot_split(nlev)
        return
    end subroutine report

    subroutine hotspot_split(k)
        ! Distinguish an axisymmetric (zonal-mean) hot spot from a grid-point
        ! spike at level k: report the warmest zonal mean and, separately, the
        ! largest departure of any point from ITS OWN latitude's zonal mean.
        integer, intent(in) :: k
        real(wp) :: zm(grd%nlat), zmax, devmax, d
        integer  :: i, j, jzmax, idev, jdev
        do j = 1, grd%nlat
            zm(j) = sum(now%temp_g(:,j,k))/real(grd%nlon, wp)
        end do
        zmax = -1.0e30_wp; jzmax = 1
        do j = 1, grd%nlat
            if (zm(j) > zmax) then; zmax = zm(j); jzmax = j; end if
        end do
        devmax = 0.0_wp; idev = 1; jdev = 1
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                d = abs(now%temp_g(i,j,k) - zm(j))
                if (d > devmax) then; devmax = d; idev = i; jdev = j; end if
            end do
        end do
        write(*,"(a,i0,a,f8.2,a,i0,a,f8.2,a,i0,a,i0,a)") &
            "   L", k, " zonal-mean max ", zmax, " K (j=", jzmax, &
            ");  max |T - zonalmean| ", devmax, " K @(", idev, ",", jdev, ")"
        return
    end subroutine hotspot_split

    subroutine level_stats(k, tmean, tmx, imx, jmx, tmn, imn, jmn)
        integer, intent(in)  :: k
        real(wp), intent(out) :: tmean, tmx, tmn
        integer, intent(out) :: imx, jmx, imn, jmn
        real(wp) :: s, wsum, w
        integer  :: i, j
        tmx = -1.0e30_wp; tmn = 1.0e30_wp; imx=1; jmx=1; imn=1; jmn=1
        s = 0.0_wp; wsum = 0.0_wp
        do j = 1, grd%nlat
            w = grd%area(1,j)
            do i = 1, grd%nlon
                if (now%temp_g(i,j,k) > tmx) then
                    tmx = now%temp_g(i,j,k); imx = i; jmx = j
                end if
                if (now%temp_g(i,j,k) < tmn) then
                    tmn = now%temp_g(i,j,k); imn = i; jmn = j
                end if
                s = s + now%temp_g(i,j,k)*grd%area(i,j); wsum = wsum + grd%area(i,j)
            end do
        end do
        tmean = s/wsum
        return
    end subroutine level_stats

    real(wp) function gmean_lev(f, k) result(m)
        real(wp), intent(in) :: f(:,:,:)
        integer,  intent(in) :: k
        real(wp) :: s, wsum
        integer  :: i, j
        s = 0.0_wp; wsum = 0.0_wp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                s = s + f(i,j,k)*grd%area(i,j); wsum = wsum + grd%area(i,j)
            end do
        end do
        m = s/wsum
        return
    end function gmean_lev

    subroutine locate_nan()
        integer :: i, j, k
        do k = 1, nlev
            do j = 1, grd%nlat
                do i = 1, grd%nlon
                    if (ts%wrk%t_g(i,j,k) /= ts%wrk%t_g(i,j,k)) then
                        write(*,"(a,i0,a,i0,a,i0,a)") "     first T NaN at level ", &
                            k, "  (i,j)=(", i, ",", j, ")"
                        return
                    end if
                end do
            end do
        end do
        return
    end subroutine locate_nan

end program rce_long
