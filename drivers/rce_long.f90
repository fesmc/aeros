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

    use aeros_defs,     only : dp, wp, wp_sh, p0, R_d, grav, S0, cp_d, L_v, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_timestep
    use aeros_ocean,    only : aeros_ocean_init
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
    logical  :: l_diag = .TRUE.       ! per-term heating split at the hot latitude
    logical  :: l_dry_adjust = .TRUE. ! dry convective adjustment before the moist scheme
    logical  :: l_uniform_insol = .FALSE. ! flatten insolation to its global mean (no meridional gradient)
    logical  :: l_nonrotating = .FALSE.   ! zero Coriolis (no jet organization; RCE vehicle)
    real(wp) :: vdiff_k0 = 10.0_wp, vdiff_sigma = 0.7_wp
    ! Model-top sponge knobs (defaults match aeros_timestep_class). Exposed to
    ! test the top thermal-wind blow-up: a stronger/deeper sponge that delays or
    ! removes it confirms the terminal event is model-top dynamical.
    real(wp) :: sponge_kr = 1.0_wp/43200.0_wp   ! max Rayleigh rate [s-1]
    real(wp) :: sponge_kt = 1.0_wp/86400.0_wp   ! max Newtonian rate [s-1]
    real(wp) :: sponge_sigma = 0.12_wp          ! sponge ramp top [sigma]
    real(wp) :: seed_asym = 0.0_wp    ! zonal-asymmetry seed amplitude [K-ish]
    real(wp) :: albedo = 0.06_wp      ! surface broadband albedo (cloud proxy knob)
    real(wp) :: co2_ppm = 280.0_wp
    integer  :: ocean_mode = 0       ! 0 prescribed SST, 1 slab
    real(wp) :: ocean_depth = 10.0_wp

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

    ! Non-rotating vehicle: with no Coriolis there is nothing to organise a
    ! meridional asymmetry into a geostrophic/angular-momentum jet, so with
    ! uniform insolation and a uniform IC every column runs the same RCE -- the
    ! clean single-column validation vehicle, reusing the full physics stack.
    if (l_nonrotating) grd%coriolis = 0.0_wp

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
    ts%cnv%dry_adjust = l_dry_adjust
    ts%surf%enabled = l_surf
    ts%surf%c_h = c_h; ts%surf%c_e = c_e; ts%surf%u_min = u_min
    ts%rad%enabled  = l_rad
    ts%sponge_on    = l_sponge
    call aeros_timestep_set_sponge(ts, vg, sponge_kr, sponge_kt, sponge_sigma)
    ts%vd%enabled   = l_vdiff
    ts%vd%k0        = vdiff_k0
    ts%vd%sigma     = vdiff_sigma
    ts%rad%albedo   = albedo
    ts%rad%co2_ppm  = co2_ppm
    ts%ocn%mode     = ocean_mode
    ts%ocn%depth    = ocean_depth
    call aeros_ocean_init(ts%ocn, grd)   ! recompute C for the chosen depth/mode

    ! Per-term heating diagnostics: split the forward-split physics back into
    ! surface/convection/condensation/radiation, capture vdiff's implicit change,
    ! and separate the vertical-advective (ventilation) part of the dynamical
    ! heating -- reported zonal-mean at the hot latitude by term_table below.
    if (l_diag) call aeros_timestep_enable_diag(ts)

    ! Uniform-insolation test: flatten sw_toa and coszen to their area-weighted
    ! global means, removing the equator-pole gradient (hence the Hadley jet and
    ! thermal wind) while keeping the global energy input (~S0/4). If the RCE
    ! bounds only here, the runaway is the axisymmetric meridional dynamics, not
    ! the column physics.
    if (l_uniform_insol) then
        block
            real(wp) :: sw_m, cz_m, wsum, w
            integer  :: jj
            sw_m = 0.0_wp; cz_m = 0.0_wp; wsum = 0.0_wp
            do jj = 1, grd%nlat
                w = sum(grd%area(:,jj))
                sw_m = sw_m + ts%rad%sw_toa(jj)*w
                cz_m = cz_m + ts%rad%coszen(jj)*w
                wsum = wsum + w
            end do
            ts%rad%sw_toa(:)  = sw_m/wsum
            ts%rad%coszen(:)  = cz_m/wsum
            write(*,"(a,f7.2,a,f6.3)") " rce_long:: uniform insolation SWin ", &
                sw_m/wsum, " W/m2  coszen ", cz_m/wsum
        end block
    end if

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
    ! Initial equator-pole gradient (l=2). Skipped under uniform insolation: with
    ! no insolation gradient the intended state is horizontally uniform, so every
    ! latitude runs the identical column RCE with nothing to spin up a jet -- the
    ! clean single-column validation vehicle.
    if (.not. l_uniform_insol) &
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
        call nml_read(nmlfile, "rce", "l_dry_adjust", l_dry_adjust)
        call nml_read(nmlfile, "rce", "l_uniform_insol", l_uniform_insol)
        call nml_read(nmlfile, "rce", "l_nonrotating", l_nonrotating)
        call nml_read(nmlfile, "rce", "sponge_kr", sponge_kr)
        call nml_read(nmlfile, "rce", "sponge_kt", sponge_kt)
        call nml_read(nmlfile, "rce", "sponge_sigma", sponge_sigma)
        call nml_read(nmlfile, "rce", "l_cnd", l_cnd)
        call nml_read(nmlfile, "rce", "l_rad", l_rad)
        call nml_read(nmlfile, "rce", "l_sponge", l_sponge)
        call nml_read(nmlfile, "rce", "l_vdiff", l_vdiff)
        call nml_read(nmlfile, "rce", "l_diag", l_diag)
        call nml_read(nmlfile, "rce", "vdiff_k0", vdiff_k0)
        call nml_read(nmlfile, "rce", "vdiff_sigma", vdiff_sigma)
        call nml_read(nmlfile, "rce", "seed_asym", seed_asym)
        call nml_read(nmlfile, "rce", "albedo", albedo)
        call nml_read(nmlfile, "rce", "co2_ppm", co2_ppm)
        call nml_read(nmlfile, "rce", "ocean_mode", ocean_mode)
        call nml_read(nmlfile, "rce", "ocean_depth", ocean_depth)
        return
    end subroutine read_config

    subroutine report(n)
        integer, intent(in) :: n
        real(wp) :: tmx, tmn, tmean, umax
        real(wp) :: hphys, hcnd
        integer  :: imx, jmx, imn, jmn, kum, jum, ium
        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        ! max|u| with its location: if it lives at the top level and migrates
        ! poleward as it grows, the terminal blow-up is a top thermal-wind jet.
        call locate_umax(umax, ium, jum, kum)
        write(*,"(a)") ""
        write(*,"(a,i0,a,f8.2,a,f7.1,a,i0,a,f6.1,a)") " -- step ", n, "  day ", &
            real(n,wp)*dt/86400.0_wp, "   max|u| ", umax, " m/s @ lev ", kum, &
            " lat ", grd%lat(jum), " N"
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
        call eddy_diag()
        call baroclinicity_diag()
        if (l_diag) call term_table()
        call energy_balance()
        return
    end subroutine report

    subroutine term_table()
        ! Per-term heating split, zonal-mean at the hot latitude j* (the warmest
        ! lowest-layer zonal mean -- the runaway column). Settles the handoff
        ! question: is the descending column being spuriously heated (large +cnv
        ! or +cnd at depth) or simply not ventilated (surf pumping in, vadv/vdiff
        ! not carrying it away)? All columns in K/day; hadv = total dynamical
        ! heating minus its vertical/adiabatic part (vadv).
        real(wp) :: zm, zmax
        real(wp) :: hs, hcv, hcd, hr, hvd, hva, hdyn
        real(wp) :: olrj, swupj, swinj
        integer  :: k, j, jstar
        zmax = -1.0e30_wp; jstar = 1
        do j = 1, grd%nlat
            zm = sum(now%temp_g(:,j,nlev))/real(grd%nlon, wp)
            if (zm > zmax) then; zmax = zm; jstar = j; end if
        end do
        write(*,"(a,i0,a,f6.1,a,f7.2,a)") &
            "   per-term heating [K/day] at hot latitude j=", jstar, &
            " (lat ", grd%lat(jstar), " N, T_low ", zmax, " K):"
        write(*,"(a)") "   lev    Tzm    qzm[g/kg]    surf     cnv     cnd" // &
                       "     rad   vdiff    vadv    hadv"
        do k = 1, nlev
            hs  = zmean_at(ts%wrk%dt_surf,  k, jstar)/dt*86400.0_wp
            hcv = zmean_at(ts%wrk%dt_cnv,   k, jstar)/dt*86400.0_wp
            hcd = zmean_at(ts%wrk%dt_cnd,   k, jstar)/dt*86400.0_wp
            hr  = zmean_at(ts%wrk%dt_rad,   k, jstar)/dt*86400.0_wp
            hvd = zmean_at(ts%wrk%dt_vdiff, k, jstar)/dt*86400.0_wp
            hva = zmean_at(ts%wrk%dt_vadv,  k, jstar)*86400.0_wp
            hdyn= zmean_at(ts%wrk%dtdt,     k, jstar)*86400.0_wp
            ! Tzm/qzm: the zonal-mean T and q at j* that RADIATION sees -- a
            ! vertical sawtooth here is the garbage-in behind the huge rad column.
            write(*,"(i6,f8.2,f10.4,3x,7f8.2)") k, &
                sum(now%temp_g(:,jstar,k))/real(grd%nlon,wp), &
                sum(now%qv_g(:,jstar,k))/real(grd%nlon,wp)*1000.0_wp, &
                hs, hcv, hcd, hr, hvd, hva, hdyn - hva
        end do
        ! LOCAL TOA at the hot latitude: does OLR rise as T_low climbs, or does it
        ! saturate (runaway-greenhouse signature -- absorbed SW can't be shed)?
        if (l_rad) then
            olrj  = sum(ts%rad%olr(:,jstar))/real(grd%nlon,wp)
            swupj = sum(ts%rad%sw_up_toa(:,jstar))/real(grd%nlon,wp)
            swinj = ts%rad%sw_toa(jstar)
            write(*,"(a,f7.2,a,f7.2,a,f7.2,a,f7.2,a)") &
                "   hot-lat TOA: SWin ", swinj, "  SWabs ", swinj-swupj, &
                "  OLR ", olrj, "  net ", swinj-swupj-olrj, " W/m2"
        end if
        call buoyancy_profile(jstar)
        return
    end subroutine term_table

    subroutine buoyancy_profile(jstar)
        ! Why the surface heat is (or is not) ventilated by convection: the
        ! boundary-layer parcel MSE h_b vs the environment saturated MSE h*_env(k)
        ! at the hot latitude. Deep convection needs h_b > h*_env over a deep band
        ! aloft; if the buoyancy is positive only near the surface (or nowhere),
        ! the surface heat cannot be carried up to the radiating levels and it
        ! traps. RH shows whether the boundary layer is drying out of buoyancy.
        integer, intent(in) :: jstar
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), phi(nlev)
        real(wp) :: tzm(nlev), qzm(nlev), qs, dqs, hb, hstar, rh, psj, cp, lv
        integer  :: k
        cp = real(cp_d, wp); lv = real(L_v, wp)
        psj = sum(now%ps(:,jstar))/real(grd%nlon, wp)
        call aeros_vgrid_pressure(vg, psj, phalf, pfull, dpc)
        do k = 1, nlev
            tzm(k) = sum(now%temp_g(:,jstar,k))/real(grd%nlon, wp)
            qzm(k) = sum(now%qv_g(:,jstar,k))/real(grd%nlon, wp)
        end do
        phi(nlev) = 0.0_wp
        do k = nlev-1, 1, -1
            phi(k) = phi(k+1) + R_d*0.5_wp*(tzm(k)+tzm(k+1))*log(pfull(k+1)/pfull(k))
        end do
        hb = cp*tzm(nlev) + phi(nlev) + lv*qzm(nlev)
        write(*,"(a,f8.1,a)") "   BL parcel h_b = ", hb/1000.0_wp, &
            " kJ/kg;  buoyancy hb-h*_env (>0 = convecting band):"
        write(*,"(a)") "   lev   RH[%]   h*_env[kJ/kg]   hb-h*[kJ/kg]"
        do k = 1, nlev
            call aeros_qsat(tzm(k), pfull(k), qs, dqs)
            rh = qzm(k)/max(qs, 1.0e-12_wp)*100.0_wp
            hstar = cp*tzm(k) + phi(k) + lv*qs
            write(*,"(i6,f8.1,f15.1,f14.2)") k, rh, hstar/1000.0_wp, (hb-hstar)/1000.0_wp
        end do
        return
    end subroutine buoyancy_profile

    real(wp) function zmean_at(f, k, j) result(m)
        ! Zonal mean of grid field f at level k, latitude row j.
        real(wp), intent(in) :: f(:,:,:)
        integer,  intent(in) :: k, j
        m = sum(f(:,j,k))/real(grd%nlon, wp)
        return
    end function zmean_at

    subroutine baroclinicity_diag()
        ! Vertical structure of the baroclinicity: per level, the meridional
        ! temperature contrast (max-min of the zonal mean over latitude) and the
        ! peak zonal-mean zonal wind |[u]|. If both pile up at the lowest level
        ! and fall off aloft, the baroclinicity is SURFACE-TRAPPED -- shallow,
        ! small-scale unstable modes that T21 resolves poorly and that feed the
        ! low-level jet; then higher resolution is the fix. If they grow upward
        ! (a deep, upper-level jet as on Earth), the eddies should organise and
        ! the issue is the seed/forcing, not resolution.
        integer  :: j, k
        real(wp) :: tzm, uzm, tmax, tmin, upk
        write(*,"(a)") "   baroclinicity  lev   dT_merid[K]   peak|[u]|[m/s]"
        do k = 1, nlev
            tmax = -1.0e30_wp; tmin = 1.0e30_wp; upk = 0.0_wp
            do j = 1, grd%nlat
                tzm = sum(now%temp_g(:,j,k))/real(grd%nlon, wp)
                uzm = sum(now%u(:,j,k))/real(grd%nlon, wp)
                if (tzm > tmax) tmax = tzm
                if (tzm < tmin) tmin = tzm
                if (abs(uzm) > upk) upk = abs(uzm)
            end do
            write(*,"(a,i6,f13.2,f16.2)") "   ", k, tmax - tmin, upk
        end do
        return
    end subroutine baroclinicity_diag

    subroutine eddy_diag()
        ! Do baroclinic eddies grow (m>0), i.e. is there anything to flux heat
        ! and momentum meridionally? Report, at a mid-troposphere level, the RMS
        ! eddy temperature T' = T - zonalmean and eddy KE 0.5(u'^2+v'^2), plus the
        ! zonal-mean meridional eddy heat flux [v'T'] (the term that relaxes the
        ! equator-pole gradient). Growing with time = baroclinic instability is
        ! working; decaying = the seed is damped and the run stays axisymmetric.
        integer  :: i, j, k
        real(wp) :: tzm, uzm, vzm, tp, up, vp
        real(wp) :: t2, ke, vt, w, wsum
        k = max(1, nlev/2)                 ! ~mid-troposphere
        t2 = 0.0_wp; ke = 0.0_wp; vt = 0.0_wp; wsum = 0.0_wp
        do j = 1, grd%nlat
            tzm = sum(now%temp_g(:,j,k))/real(grd%nlon, wp)
            uzm = sum(now%u(:,j,k))/real(grd%nlon, wp)
            vzm = sum(now%v(:,j,k))/real(grd%nlon, wp)
            w   = grd%area(1,j)
            do i = 1, grd%nlon
                tp = now%temp_g(i,j,k) - tzm
                up = now%u(i,j,k) - uzm
                vp = now%v(i,j,k) - vzm
                t2 = t2 + w*tp*tp
                ke = ke + w*0.5_wp*(up*up + vp*vp)
                vt = vt + w*vp*tp
            end do
            wsum = wsum + w*real(grd%nlon, wp)
        end do
        write(*,"(a,i0,a,es9.2,a,es9.2,a,es9.2)") "   eddy(lev ", k, "): RMS T'' ", &
            sqrt(t2/wsum), " K  eddyKE ", ke/wsum, " m2/s2  [v''T''] ", vt/wsum
        return
    end subroutine eddy_diag

    subroutine locate_umax(umax, ium, jum, kum)
        ! max|u| and its (i,j,k) -- to see whether the growing jet sits at the
        ! model top (the thermal-wind blow-up) and where in latitude.
        real(wp), intent(out) :: umax
        integer,  intent(out) :: ium, jum, kum
        real(wp) :: a
        integer  :: i, j, k
        umax = -1.0_wp; ium = 1; jum = 1; kum = 1
        do k = 1, nlev
            do j = 1, grd%nlat
                do i = 1, grd%nlon
                    a = abs(now%u(i,j,k))
                    if (a > umax) then; umax = a; ium = i; jum = j; kum = k; end if
                end do
            end do
        end do
        return
    end subroutine locate_umax

    subroutine energy_balance()
        ! Global-mean TOA and surface energy budget. A secular drift of the
        ! column temperature is a TOA imbalance (absorbed SW - OLR /= 0); this
        ! separates an energy-balance problem from a mixing/numerics one.
        real(wp) :: olr, swup, swin, shf, lhf, wsum
        integer  :: j
        if (.not. (l_rad .or. l_surf)) return
        olr = 0.0_wp; swup = 0.0_wp; swin = 0.0_wp; shf = 0.0_wp; lhf = 0.0_wp
        if (l_rad) then
            olr  = gmean2(ts%rad%olr)
            swup = gmean2(ts%rad%sw_up_toa)
            wsum = 0.0_wp
            do j = 1, grd%nlat
                swin = swin + ts%rad%sw_toa(j)*sum(grd%area(:,j))
                wsum = wsum + sum(grd%area(:,j))
            end do
            swin = swin/wsum
        end if
        if (l_surf) then
            shf = gmean2(ts%surf%shf); lhf = gmean2(ts%surf%lhf)
        end if
        write(*,"(a,f7.1,a,f7.1,a,f7.2,a,f7.1,a,f7.1,a)") &
            "   TOA: SWin ", swin, "  OLR ", olr, "  net ", swin-swup-olr, &
            " W/m2 | sfc SH ", shf, "  LH ", lhf, " W/m2"
        return
    end subroutine energy_balance

    real(wp) function gmean2(f) result(m)
        real(wp), intent(in) :: f(:,:)
        real(wp) :: s, wsum
        integer  :: i, j
        s = 0.0_wp; wsum = 0.0_wp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                s = s + f(i,j)*grd%area(i,j); wsum = wsum + grd%area(i,j)
            end do
        end do
        m = s/wsum
        return
    end function gmean2

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
