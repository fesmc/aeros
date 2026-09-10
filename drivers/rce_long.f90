program rce_long
    ! Long coupled radiative-convective integration, instrumented, for chasing
    ! the lowest-layer hot-spot instability (m2_handoff.md task #7) and, since
    ! feat/seasonal, the real-Earth climatology harness (mwm/D_realearth): run
    ! for model years with a per-level T trajectory so a blow-up can be
    ! localized in level, space and time, second-half diabatic/eddy
    ! accumulators, and a monthly climatology dump for the ERA5 comparison.
    !
    ! Wave 4a.2 of the coupling API (docs/coupling_api.md sec. 5, 8, 14): the
    ! model is initialised through aeros_lib_init from ONE namelist file --
    ! every dynamics / physics knob in its subsystem group (&aeros, &physics,
    ! &aeros_moisture, &surface, &radiation, &clouds, &ocean, &land, &vdiff,
    ! &bd, ...; schema and calibrated defaults in input/aeros_defaults.nml) --
    ! and this driver owns only its `&rce` group: run length, report cadence,
    ! the monthly-climatology and restart paths and two diagnostic probes
    ! (uniform insolation, analytic heating). Everything it reads back (dt,
    ! level count, the enabled flags, doy0, ...) comes from `ams%par`,
    ! `ams%phys` and `ams%ts%<sub>` -- there are no shadow copies.
    !
    ! Wave 4b: every prescribed field -- the seasonal SST table, the
    ! land-sea mask and bare albedo, the dust and sulfate maps, the
    ! orography -- is a boundary provider (&boundary, ams%bounds) read by
    ! the library; this driver only refreshes the providers once a day
    ! (aeros_boundary_refresh + aeros_boundary_apply, the cadence the old
    ! set_seasonal_sst had) and ramps the band-limited orography in over
    ! &boundary phis_ramp_days. The time loop, the restart calls and the
    ! dumps stay here (Waves 4c-4d move them behind the library).
    !
    ! Each report prints, per level: global-mean T, the global max T with its
    ! (i,j) location, and the global min T -- the max-location wandering and a
    ! single level running away is the hot-spot signature. It also prints the
    ! per-level global-mean forward-split heating (surface+convection+radiation,
    ! wrk%dt_phys) and centered heating (condensation, wrk%dtdt) as K/day, so the
    ! term driving a level is visible directly. Stops at the first NaN, naming
    ! the level and step.

    use aeros_defs,     only : dp, wp, wp_sh, R_d, grav, cp_d, year_seconds
    use aeros_spectral, only : aeros_sht_analysis, aeros_sht_synthesis
    use aeros_vertical, only : aeros_vgrid_pressure
    use aeros_condensation, only : aeros_qsat, aeros_latent
    use aeros_cloud,    only : aeros_cloud_diagnose
    use aeros_timestep, only : aeros_timestep_step, aeros_timestep_diagnose, &
                                aeros_timestep_enable_diag, aeros_timestep_set_phis, &
                                aeros_timestep_set_qforce, aeros_timestep_write_restart
    use aeros_topography, only : aeros_topography_scale
    use aeros_boundary, only : aeros_boundary_refresh, aeros_boundary_kind_name, BC_NONE
    use aeros_diag_acc, only : aeros_diag_acc_init_monthly, AEROS_DIAG_NMON, &
                                aeros_diag_acc_energy, aeros_diag_acc_water, &
                                aeros_diag_acc_second_half, aeros_diag_acc_monthly, &
                                aeros_diag_acc_imon
    use aeros,          only : aeros_class, aeros_lib_init, aeros_lib_finalize, &
                                aeros_lib_restart_read, aeros_set_rce_warm, &
                                aeros_boundary_apply
    use nml,            only : nml_read
    use ncio,           only : nc_create, nc_write_dim, nc_write

    implicit none

    ! The model: parameters, physics flags, transforms, grids, state, integrator
    ! and the diagnostic accumulators, all owned by the library (aeros_class).
    type(aeros_class) :: ams

    character(len=512) :: nmlfile
    ! The parameter schema every read falls back to (nml_read defaults_file):
    ! the calibrated real-Earth state, so a case file lists only its changes.
    character(len=*), parameter :: defaults_file = "input/aeros_defaults.nml"

    ! End-of-run zonal-mean RH/cf dump (arg 2; "" = off). NOT a namelist key.
    ! The comparison against ERA5 is done offline by scripts/rce_humidity_vs_era5.jl.
    character(len=512) :: rh_out = "output/rce_rh.nc"

    ! --- &rce: the driver's own keys (defaults in input/aeros_defaults.nml) ---
    ! The character keys are initialised to "": nml_read leaves the variable
    ! untouched for an empty-string value, so "" in the schema means "off".
    integer  :: nstep                       ! run length [steps]
    integer  :: print_every                 ! trajectory-report cadence [steps]
    logical  :: l_diag                      ! per-term heating split + 2nd-half accumulators
    ! Monthly-mean climatology output ("" = off): zonal-mean T/u/RH/v
    ! (lat,lev,month), the 2D surface fields and the eddy covariances binned by
    ! calendar month (from the seasonal day-of-year) for the ERA5 comparison;
    ! months are accumulated only once the elapsed day exceeds mon_spinup_days.
    character(len=512) :: monthly_out = ""
    real(wp) :: mon_spinup_days
    ! Checkpoint / restart: resume from restart_in ("" = cold start), write
    ! checkpoints to restart_out ("" = none) every restart_interval steps (0 =
    ! only at the end).
    character(len=512) :: restart_in = "", restart_out = ""
    integer  :: restart_interval
    ! Diagnostic probes: flatten the insolation to its global mean (with the
    ! gradient-free IC), and the analytic Q(lat,sigma) heating [K/day] coupled
    ! on the in-solve seam (run with all physics off).
    logical  :: l_uniform_insol, l_qforce
    real(wp) :: qforce_amp

    ! --- driver state ---------------------------------------------------------
    ! The orography ramp (&physics l_topography): the band-limited full field
    ! is the library's phis provider buffer (ams%bounds%phis%buffer, set by
    ! aeros_boundary_apply); it is ramped in linearly over &boundary
    ! phis_ramp_days (a pure function of elapsed time, aeros_topography_scale,
    ! hence restart-safe) and fed to the dynamics as phis2.
    real(wp), allocatable :: phis2(:,:)       ! surface geopotential fed to the dynamics
    real(wp) :: tscale, tscale_prev           ! current / previous ramp factor
    integer  :: nstep_day                     ! steps per day: the boundary refresh cadence
    real(wp) :: restart_time = 0.0_wp
    integer  :: n0 = 0                        ! absolute step offset (nonzero on restart)
    integer, parameter :: NMON = AEROS_DIAG_NMON
    integer  :: n, k, ierr
    logical  :: blew_up

    nmlfile = "rce.nml"
    if (command_argument_count() >= 1) call get_command_argument(1, nmlfile)
    if (command_argument_count() >= 2) call get_command_argument(2, rh_out)
    call read_config()

    ! === Library init =======================================================
    ! Dynamics parameters, the &physics flags, the transform pool, grids,
    ! state, every physics subsystem from its own namelist group (in
    ! aeros_timestep_init's fixed order, incl. the land-albedo composition),
    ! the boundary providers and the diagnostic accumulators. The cold start
    ! is the warm, humid, conditionally unstable RCE state (aeros_set_rce_warm,
    ! seed from &aeros seed_asym); a restart installs the resting shape and
    ! reads the checkpoint below, once the driver-side boundary state is set.
    if (len_trim(restart_in) > 0) then
        call aeros_lib_init(ams, nmlfile, defaults_file=defaults_file, &
                            ic_source="restart", ierr=ierr)
    else
        call aeros_lib_init(ams, nmlfile, defaults_file=defaults_file, &
                            ic_source="rce_warm", ierr=ierr)
    end if
    if (ierr /= 0) then
        write(*,"(a,i0)") " rce_long:: error: aeros_lib_init failed, ierr = ", ierr
        error stop 1
    end if

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " rce_long:: T", ams%par%trunc, " L", ams%par%nlev, &
                                       "  grid ", ams%grd%nlon, "x", ams%grd%nlat
    write(*,"(a,i0,a,f6.1,a,f5.1,a,i0)") "   nstep=", nstep, " dt=", ams%par%dt, &
        " tau_diff=", ams%par%tau_diff, " ndiff=", ams%par%ndiff
    write(*,"(a,5(l1,1x))") "   surf/cnv/cnd/rad/sponge = ", &
        ams%phys%l_surf, ams%phys%l_cnv, ams%phys%l_cnd, ams%phys%l_rad, ams%par%l_sponge
    if (ams%ts%bd%enabled) write(*,"(a,f5.2,a,f6.3,a,f6.1)") &
        " rce_long:: lid Brewer-Dobson w* ON  scale ", ams%ts%bd%scale, &
        "  per_K ", ams%ts%bd%per_k, "  t_ref ", ams%ts%bd%t_ref
    if (ams%phys%l_land) then
        write(*,"(a,a,a,a)") " rce_long:: land-sea mask provider ", &
            trim(aeros_boundary_kind_name(ams%bounds%lsm%kind)), " ", trim(ams%bounds%lsm%file)
        write(*,"(a,i0,a,i0,a,f4.2,a)") " rce_long:: land cells ", count(ams%ts%land%mask), &
            " of ", ams%grd%nlon*ams%grd%nlat, " (mask threshold ", ams%ts%land%lsm_threshold, ")"
    end if

    ! Prescribed SST (&ocean mode = 0): the sst provider (&boundary; the
    ! seasonal ERA5 table by default, refreshed daily below) has already
    ! been applied to ams%ts%ocn%sst by aeros_lib_init; BC_NONE leaves the
    ! ocean module's analytic profile. In slab mode the SST is prognostic and
    ! aeros_ocean_step owns it, whatever the provider says.
    if (ams%ts%ocn%mode == 0 .and. ams%bounds%sst%kind /= BC_NONE) then
        write(*,"(a,a,a,a)") " rce_long:: prescribed SST provider ", &
            trim(aeros_boundary_kind_name(ams%bounds%sst%kind)), " ", trim(ams%bounds%sst%file)
    end if

    ! Per-term heating diagnostics: split the forward-split physics back into
    ! surface/convection/condensation/radiation, capture vdiff's implicit change,
    ! and separate the vertical-advective (ventilation) part of the dynamical
    ! heating -- reported zonal-mean at the hot latitude by term_table below.
    if (l_diag) call aeros_timestep_enable_diag(ams%ts)

    ! Uniform-insolation test: flatten sw_toa and coszen to their area-weighted
    ! global means, removing the equator-pole gradient (hence the Hadley jet and
    ! thermal wind) while keeping the global energy input (~S0/4). If the RCE
    ! bounds only here, the runaway is the axisymmetric meridional dynamics, not
    ! the column physics. The cold start is then re-installed WITHOUT the l=2
    ! equator-pole gradient: with no insolation gradient the intended state is
    ! horizontally uniform, so every latitude runs the identical column RCE.
    if (l_uniform_insol) then
        block
            real(wp) :: sw_m, cz_m, wsum, w
            integer  :: jj
            sw_m = 0.0_wp; cz_m = 0.0_wp; wsum = 0.0_wp
            do jj = 1, ams%grd%nlat
                w = sum(ams%grd%area(:,jj))
                sw_m = sw_m + ams%ts%rad%sw_toa(jj)*w
                cz_m = cz_m + ams%ts%rad%coszen(jj)*w
                wsum = wsum + w
            end do
            ams%ts%rad%sw_toa(:)  = sw_m/wsum
            ams%ts%rad%coszen(:)  = cz_m/wsum
            write(*,"(a,f7.2,a,f6.3)") " rce_long:: uniform insolation SWin ", &
                sw_m/wsum, " W/m2  coszen ", cz_m/wsum
        end block
        if (len_trim(restart_in) == 0) call aeros_set_rce_warm(ams, gradient=.FALSE.)
    end if

    ! Surface geopotential (lower boundary). Aquaplanet phis = 0 by default; with
    ! l_topography the phis provider (ERA5 orography by default) has been read
    ! and band-limited to the model truncation by the library
    ! (aeros_boundary_apply; ams%bounds%phis%buffer holds the band-limited
    ! field) and is ramped in here over phis_ramp_days.
    allocate(phis2(ams%grd%nlon, ams%grd%nlat)); phis2 = 0.0_wp
    if (ams%phys%l_topography) then
        write(*,"(a,a,a,a)") " rce_long:: topography provider ", &
            trim(aeros_boundary_kind_name(ams%bounds%phis%kind)), " ", trim(ams%bounds%phis%file)
        write(*,"(a,i0,a,es11.3,a,es11.3,a,f6.1,a)") " rce_long:: T", ams%par%trunc, &
            "   phis min ", minval(ams%bounds%phis%buffer), "  max ", maxval(ams%bounds%phis%buffer), &
            " m2/s2   ramp ", ams%bounds%phis_ramp_days, " days"
        ! Dump the band-limited orography on the model grid for a visual sanity
        ! check (correct continents/mountains, sane magnitude).
        call nc_create("output/rce_phis.nc")
        call nc_write_dim("output/rce_phis.nc", "lon", x=ams%grd%lon, units="degrees_east")
        call nc_write_dim("output/rce_phis.nc", "lat", x=ams%grd%lat, units="degrees_north")
        call nc_write("output/rce_phis.nc", "phis", ams%bounds%phis%buffer, dim1="lon", dim2="lat", &
                      units="m2 s-2", long_name="surface geopotential (band-limited)")
    end if
    ! Initial feed at t = 0 (scale = 0 while ramping, so the run starts flat).
    tscale      = aeros_topography_scale(ams%phys%l_topography, 0.0_wp, ams%bounds%phis_ramp_days)
    tscale_prev = tscale
    phis2 = tscale*ams%bounds%phis%buffer
    call aeros_timestep_set_phis(ams%ts, phis2)

    ! Diagnostic prescribed heating: an analytic Q(lat,sigma), equatorial and
    ! area-weighted zero-mean in latitude (no net global heating -> no drift),
    ! sin(pi*sigma) in the vertical (mid-troposphere peak, zero at top/surface).
    ! Coupled on the in-solve seam to isolate the dry core's heating->omega
    ! response from moist physics. Run with all physics off.
    if (l_qforce) then
        block
            real(wp), allocatable :: q3(:,:,:), slat(:)
            real(wp) :: amp, mlat, wsum, latr, pi
            integer  :: iq, jq, kq
            pi = 4.0_wp*atan(1.0_wp)
            allocate(q3(ams%grd%nlon, ams%grd%nlat, ams%par%nlev), slat(ams%grd%nlat))
            amp = qforce_amp/86400.0_wp                     ! K/day -> K/s
            wsum = 0.0_wp; mlat = 0.0_wp
            do jq = 1, ams%grd%nlat
                latr    = ams%grd%lat(jq)
                slat(jq) = exp(-(latr/15.0_wp)**2)          ! equatorial heating
                mlat    = mlat + slat(jq)*cos(latr*pi/180.0_wp)
                wsum    = wsum + cos(latr*pi/180.0_wp)
            end do
            mlat = mlat/wsum                                ! area-weighted mean
            do kq = 1, ams%par%nlev
                do jq = 1, ams%grd%nlat
                    do iq = 1, ams%grd%nlon
                        q3(iq,jq,kq) = amp*(slat(jq) - mlat) &
                                        *sin(pi*ams%vgrid%sigma_full(kq))
                    end do
                end do
            end do
            call aeros_timestep_set_qforce(ams%ts, q3)
            write(*,'(a,f6.2,a)') " rce_long:: prescribed Q ON, amp=", qforce_amp, &
                                    " K/day (physics MUST be off for a clean dry test)"
            deallocate(q3, slat)
        end block
    end if

    ! === Restart: load the checkpoint over the initialised shapes ============
    ! Overwrites the time levels, humidity, ocean and radiation cache with the
    ! saved state and resumes from the saved nstep. The cold start (rce_warm)
    ! was installed by aeros_lib_init above.
    if (len_trim(restart_in) > 0) then
        call aeros_lib_restart_read(ams, restart_in, ierr)
        if (ierr /= 0) then
            write(*,"(a,i0)") " rce_long:: error: aeros_lib_restart_read failed, ierr = ", ierr
            error stop 1
        end if
        n0 = ams%ts%nstep
        restart_time = ams%time*real(year_seconds, wp)
        write(*,"(a,a,a,i0,a,f10.1,a)") " rce_long:: restarted from ", &
            trim(restart_in), "  (nstep=", n0, ", time=", restart_time, " s)"
        call report(n0)
    else
        call report(0)
    end if

    blew_up = .FALSE.

    ! The twelve-month climatology bins are an output choice, so they are
    ! allocated here (on the library's accumulator container, which
    ! aeros_lib_init already sized for the always-on buffers) only when a
    ! monthly dump was asked for.
    if (len_trim(monthly_out) > 0) then
        call aeros_diag_acc_init_monthly(ams%diag_acc, ams%grd%nlon, ams%grd%nlat, ams%par%nlev)
        if (.not. ams%ts%rad%seasonal) write(*,"(a)") &
            " rce_long:: WARNING monthly_out set but &radiation seasonal=.FALSE. -- "// &
            "insolation is annual-mean, so all months bin the same forcing."
    end if
    nstep_day = max(1, nint(86400.0_wp/ams%par%dt))
    do n = n0+1, n0+nstep
        ! Boundary providers: once a day, interpolate the monthly fields (SST,
        ! bare land albedo, dust) to the running day-of-year and hand every
        ! buffer to its consumer, before the step reads them in the surface
        ! fluxes and radiation (the cadence and the doy0 + n dt argument of the
        ! old set_seasonal_sst / aeros_land_set_albedo).
        if (mod(n, nstep_day) == 0) then
            call aeros_boundary_refresh(ams%bounds, ams%ts%rad%doy0 + real(n,wp)*ams%par%dt/86400.0_wp)
            call aeros_boundary_apply(ams)
        end if
        ! Advance the topography ramp: a pure function of absolute elapsed time
        ! n*dt, so it is restart-safe -- a run resumed at n0 continues the ramp at
        ! the right amplitude. Only re-set phis while the factor is still changing
        ! (during the ramp); once it reaches 1 the surface is constant.
        if (ams%phys%l_topography .and. ams%bounds%phis_ramp_days > 0.0_wp .and. tscale_prev < 1.0_wp) then
            tscale = aeros_topography_scale(.true., real(n,wp)*ams%par%dt, ams%bounds%phis_ramp_days)
            if (tscale /= tscale_prev) then
                phis2 = tscale*ams%bounds%phis%buffer
                call aeros_timestep_set_phis(ams%ts, phis2)
                tscale_prev = tscale
            end if
        end if
        call aeros_timestep_step(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)
        if (any(ams%ts%wrk%t_g /= ams%ts%wrk%t_g) .or. any(ams%now%qv_g /= ams%now%qv_g)) then
            write(*,"(a,i0,a,f7.2,a)") " *** NaN at step ", n, "  (day ", &
                real(n,wp)*ams%par%dt/86400.0_wp, ")"
            call locate_nan()
            blew_up = .TRUE.
            exit
        end if
        ! Per-step accumulators (library owns the arithmetic; the driver
        ! decides WHEN). Water/energy every step; the 2nd-half tavg and
        ! monthly bin only after the spin-up windows below.
        call aeros_diag_acc_water(ams%diag_acc, ams%ts, ams%phys, ams%grd)
        call aeros_diag_acc_energy(ams%diag_acc, ams%ts, ams%phys, ams%grd)
        if (mod(n, print_every) == 0) call report(n)
        ! Time-mean omega over the run's second half (equilibrium) for a robust
        ! zonal-mean subsidence diagnostic -- a single snapshot is too noisy.
        if (l_diag .and. n > n0 + nstep/2) then
            call aeros_diag_acc_second_half(ams%diag_acc, ams%ts, ams%now, &
                                            ams%vgrid, ams%grd, ams%par)
        end if
        ! Monthly-mean climatology sampling (feat/seasonal): bin the zonal-mean
        ! state by calendar month on a ~daily cadence, once past the spin-up
        ! window (mon_spinup_days) so the slab drift does not pollute the mean.
        ! The sampled grid fields (now%temp_g/u/v/ps) are diagnostics that only
        ! aeros_timestep_diagnose refreshes -- refresh them here, or the sample
        ! is the state at the last report (which is what silently happened up to
        ! 2026-09: months identical between prints).
        if (len_trim(monthly_out) > 0 .and. mod(n, ams%diag_acc%nacc_mon) == 0 .and. &
            real(n,wp)*ams%par%dt/86400.0_wp >= mon_spinup_days) then
            call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)
            call aeros_diag_acc_monthly(ams%diag_acc, ams%ts, ams%now, ams%vgrid, &
                                        ams%grd, ams%par, ams%phys, &
                                        aeros_diag_acc_imon(ams%ts%rad%doy0 + real(n,wp)*ams%par%dt/86400.0_wp))
        end if
        ! Periodic checkpoint (restart_interval > 0). The model time carried into
        ! the file is the absolute elapsed time n*dt.
        if (len_trim(restart_out) > 0 .and. restart_interval > 0) then
            if (mod(n, restart_interval) == 0) &
                call aeros_timestep_write_restart(ams%ts, ams%now, real(n,wp)*ams%par%dt, restart_out)
        end if
    end do

    if (.not. blew_up) write(*,"(a)") " rce_long:: completed without NaN"

    ! End-of-run checkpoint: always written when restart_out is set (this is the
    ! only checkpoint when restart_interval = 0). Skipped on a blow-up, which
    ! would only persist a NaN state.
    if (len_trim(restart_out) > 0 .and. .not. blew_up) then
        call aeros_timestep_write_restart(ams%ts, ams%now, real(n0+nstep,wp)*ams%par%dt, restart_out)
        write(*,"(a,a)") " rce_long:: wrote restart -> ", trim(restart_out)
    end if

    if (len_trim(rh_out) > 0) call dump_rh(trim(rh_out))
    if (len_trim(monthly_out) > 0 .and. .not. blew_up) call dump_monthly(trim(monthly_out))

    deallocate(phis2)
    call aeros_lib_finalize(ams, ierr)

contains

    subroutine dump_rh(fname)
        ! End-of-run zonal-mean state for the humidity-bias diagnosis (handoff
        ! step 3): the model RCE overcasts (rce_long TOA clouds-on -36 W/m2) even
        ! though the cloud scheme is validated correct on ERA5 columns (m2_results
        ! s21), so the residual is a *model* moisture bias. This writes the model's
        ! own zonal-mean RH(lat,sigma) and diagnosed cloud fraction so it can be
        ! laid against the ERA5 RH climatology (scripts/rce_humidity_vs_era5.jl) to
        ! localize where and by how much the column is too moist. Pointwise RH and
        ! cf are zonal-averaged (both are nonlinear in T,q, so we average the
        ! diagnostic, not the inputs). `cover` is the max-overlap total cloud cover
        ! per latitude -- the quantity the overcast TOA reflects directly.
        character(len=*), intent(in) :: fname
        real(wp) :: phc(0:ams%par%nlev), pfc(ams%par%nlev), dpcc(ams%par%nlev)
        real(wp) :: cf(ams%par%nlev), clwc(ams%par%nlev), ciwc(ams%par%nlev)
        real(wp) :: qsl, dqsl, colcov, rnlon
        real(wp), allocatable :: rh_zm(:,:), cf_zm(:,:), t_zm(:,:), q_zm(:,:), p_zm(:,:)
        real(wp), allocatable :: u_zm(:,:), v_zm(:,:), w_zm(:,:)
        real(wp), allocatable :: cnv_zm(:,:), cnd_zm(:,:), rad_zm(:,:), &
                                 surf_zm(:,:), qnet_zm(:,:), vadv_zm(:,:), dyn_zm(:,:), vdf_zm(:,:), &
                                 bd_zm(:,:)
        real(wp), allocatable :: cover(:), latout(:), levout(:)
        real(wp) :: kday
        integer :: ii, jj, kk

        rnlon = real(ams%grd%nlon, wp)
        allocate(rh_zm(ams%grd%nlat,ams%par%nlev), cf_zm(ams%grd%nlat,ams%par%nlev), t_zm(ams%grd%nlat,ams%par%nlev), &
                 q_zm(ams%grd%nlat,ams%par%nlev), p_zm(ams%grd%nlat,ams%par%nlev), &
                 u_zm(ams%grd%nlat,ams%par%nlev), v_zm(ams%grd%nlat,ams%par%nlev), w_zm(ams%grd%nlat,ams%par%nlev), &
                 cover(ams%grd%nlat), latout(ams%grd%nlat), levout(ams%par%nlev))
        rh_zm = 0.0_wp; cf_zm = 0.0_wp; t_zm = 0.0_wp; q_zm = 0.0_wp
        p_zm = 0.0_wp; u_zm = 0.0_wp; v_zm = 0.0_wp; w_zm = 0.0_wp; cover = 0.0_wp
        allocate(cnv_zm(ams%grd%nlat,ams%par%nlev), cnd_zm(ams%grd%nlat,ams%par%nlev), rad_zm(ams%grd%nlat,ams%par%nlev), &
                 surf_zm(ams%grd%nlat,ams%par%nlev), qnet_zm(ams%grd%nlat,ams%par%nlev), vadv_zm(ams%grd%nlat,ams%par%nlev), &
                 dyn_zm(ams%grd%nlat,ams%par%nlev), vdf_zm(ams%grd%nlat,ams%par%nlev), bd_zm(ams%grd%nlat,ams%par%nlev))
        cnv_zm = 0.0_wp; cnd_zm = 0.0_wp; rad_zm = 0.0_wp; surf_zm = 0.0_wp
        vadv_zm = 0.0_wp; dyn_zm = 0.0_wp; vdf_zm = 0.0_wp; bd_zm = 0.0_wp
        latout = ams%grd%lat(1:ams%grd%nlat)
        levout = ams%vgrid%sigma_full(1:ams%par%nlev)
        ! [K/step] time-mean -> [K/day]: /dacc%n_wacc (time), /nlon (zonal), *86400/dt.
        kday = 86400.0_wp/ams%par%dt

        do jj = 1, ams%grd%nlat
            do ii = 1, ams%grd%nlon
                call aeros_vgrid_pressure(ams%vgrid, ams%now%ps(ii,jj), phc, pfc, dpcc)
                call aeros_cloud_diagnose(ams%par%nlev, ams%now%temp_g(ii,jj,:), &
                    ams%now%qv_g(ii,jj,:), pfc, ams%now%ps(ii,jj), cf, clwc, ciwc)
                colcov = 0.0_wp
                do kk = 1, ams%par%nlev
                    call aeros_qsat(ams%now%temp_g(ii,jj,kk), pfc(kk), qsl, dqsl)
                    rh_zm(jj,kk) = rh_zm(jj,kk) + &
                        ams%now%qv_g(ii,jj,kk)/max(qsl,1.0e-12_wp)*100.0_wp
                    cf_zm(jj,kk) = cf_zm(jj,kk) + cf(kk)
                    t_zm(jj,kk)  = t_zm(jj,kk)  + ams%now%temp_g(ii,jj,kk)
                    q_zm(jj,kk)  = q_zm(jj,kk)  + ams%now%qv_g(ii,jj,kk)*1000.0_wp
                    p_zm(jj,kk)  = p_zm(jj,kk)  + pfc(kk)/100.0_wp
                    u_zm(jj,kk)  = u_zm(jj,kk)  + ams%now%u(ii,jj,kk)
                    v_zm(jj,kk)  = v_zm(jj,kk)  + ams%now%v(ii,jj,kk)
                    ! time-mean vertical pressure velocity [hPa/day], >0 = subsidence
                    if (ams%diag_acc%n_wacc > 0) w_zm(jj,kk) = w_zm(jj,kk) &
                        + ams%diag_acc%w_tavg(ii,jj,kk)/real(ams%diag_acc%n_wacc,wp)*864.0_wp
                    ! per-term diabatic heating, time-mean [K/day]
                    if (ams%diag_acc%n_wacc > 0) then
                        cnv_zm(jj,kk)  = cnv_zm(jj,kk)  + ams%diag_acc%cnv_tavg(ii,jj,kk) /real(ams%diag_acc%n_wacc,wp)*kday
                        cnd_zm(jj,kk)  = cnd_zm(jj,kk)  + ams%diag_acc%cnd_tavg(ii,jj,kk) /real(ams%diag_acc%n_wacc,wp)*kday
                        rad_zm(jj,kk)  = rad_zm(jj,kk)  + ams%diag_acc%rad_tavg(ii,jj,kk) /real(ams%diag_acc%n_wacc,wp)*kday
                        surf_zm(jj,kk) = surf_zm(jj,kk) + ams%diag_acc%surf_tavg(ii,jj,kk)/real(ams%diag_acc%n_wacc,wp)*kday
                        vadv_zm(jj,kk) = vadv_zm(jj,kk) + ams%diag_acc%vadv_tavg(ii,jj,kk)/real(ams%diag_acc%n_wacc,wp)*kday
                        dyn_zm(jj,kk)  = dyn_zm(jj,kk)  + ams%diag_acc%dyn_tavg(ii,jj,kk) /real(ams%diag_acc%n_wacc,wp)*kday
                        vdf_zm(jj,kk)  = vdf_zm(jj,kk)  + ams%diag_acc%vdf_tavg(ii,jj,kk) /real(ams%diag_acc%n_wacc,wp)*kday
                        bd_zm(jj,kk)   = bd_zm(jj,kk)   + ams%diag_acc%bd_tavg(ii,jj,kk)  /real(ams%diag_acc%n_wacc,wp)*kday
                    end if
                    colcov = max(colcov, cf(kk))     ! max-overlap column cover
                end do
                cover(jj) = cover(jj) + colcov
            end do
            rh_zm(jj,:) = rh_zm(jj,:)/rnlon
            cf_zm(jj,:) = cf_zm(jj,:)/rnlon
            t_zm(jj,:)  = t_zm(jj,:)/rnlon
            q_zm(jj,:)  = q_zm(jj,:)/rnlon
            p_zm(jj,:)  = p_zm(jj,:)/rnlon
            u_zm(jj,:)  = u_zm(jj,:)/rnlon
            v_zm(jj,:)  = v_zm(jj,:)/rnlon
            w_zm(jj,:)  = w_zm(jj,:)/rnlon
            cnv_zm(jj,:)  = cnv_zm(jj,:)/rnlon
            cnd_zm(jj,:)  = cnd_zm(jj,:)/rnlon
            rad_zm(jj,:)  = rad_zm(jj,:)/rnlon
            surf_zm(jj,:) = surf_zm(jj,:)/rnlon
            vadv_zm(jj,:) = vadv_zm(jj,:)/rnlon
            dyn_zm(jj,:)  = dyn_zm(jj,:)/rnlon
            vdf_zm(jj,:)  = vdf_zm(jj,:)/rnlon
            bd_zm(jj,:)   = bd_zm(jj,:)/rnlon
            qnet_zm(jj,:) = cnv_zm(jj,:) + cnd_zm(jj,:) + rad_zm(jj,:) + surf_zm(jj,:) + vdf_zm(jj,:) &
                            + bd_zm(jj,:)
            cover(jj)   = cover(jj)/rnlon
        end do

        call nc_create(fname)
        call nc_write_dim(fname, "lat", x=latout, units="degrees_north")
        call nc_write_dim(fname, "lev", x=levout, units="sigma")
        call nc_write(fname, "rh",    rh_zm, dim1="lat", dim2="lev", units="%", &
            long_name="zonal-mean relative humidity")
        call nc_write(fname, "cf",    cf_zm, dim1="lat", dim2="lev", units="1", &
            long_name="zonal-mean diagnosed cloud fraction")
        call nc_write(fname, "t",     t_zm,  dim1="lat", dim2="lev", units="K", &
            long_name="zonal-mean temperature")
        call nc_write(fname, "q",     q_zm,  dim1="lat", dim2="lev", units="g/kg", &
            long_name="zonal-mean specific humidity")
        call nc_write(fname, "pfull", p_zm,  dim1="lat", dim2="lev", units="hPa", &
            long_name="zonal-mean layer pressure")
        call nc_write(fname, "omega", w_zm,  dim1="lat", dim2="lev", units="hPa/day", &
            long_name="zonal-mean vertical pressure velocity (>0 subsidence)")
        call nc_write(fname, "u",     u_zm,  dim1="lat", dim2="lev", units="m/s", &
            long_name="zonal-mean zonal wind")
        call nc_write(fname, "v",     v_zm,  dim1="lat", dim2="lev", units="m/s", &
            long_name="zonal-mean meridional wind")
        if (ams%diag_acc%n_wacc > 0) then
            call nc_write(fname, "ubar", ams%diag_acc%ubar_tavg/real(ams%diag_acc%n_wacc,wp), dim1="lat", dim2="lev", &
                units="m/s", long_name="time-mean zonal-mean zonal wind (2nd half)")
            call nc_write(fname, "vbar", ams%diag_acc%vbar_tavg/real(ams%diag_acc%n_wacc,wp), dim1="lat", dim2="lev", &
                units="m/s", long_name="time-mean zonal-mean meridional wind (2nd half)")
            call nc_write(fname, "uvpr", ams%diag_acc%uv_tavg/real(ams%diag_acc%n_wacc,wp), dim1="lat", dim2="lev", &
                units="m2/s2", long_name="time-mean zonal-mean eddy momentum flux [u*v*]")
        end if
        call nc_write(fname, "cover", cover, dim1="lat", units="1", &
            long_name="max-overlap total cloud cover")
        call nc_write(fname, "q_cnv",  cnv_zm,  dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean convective heating")
        call nc_write(fname, "q_cnd",  cnd_zm,  dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean condensation heating")
        call nc_write(fname, "q_rad",  rad_zm,  dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean radiative heating")
        call nc_write(fname, "q_surf", surf_zm, dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean surface-flux heating")
        call nc_write(fname, "q_vdiff", vdf_zm, dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean vertical-diffusion heating (incl. the surface flux)")
        call nc_write(fname, "q_bd",   bd_zm,   dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean prescribed Brewer-Dobson heating (lid)")
        call nc_write(fname, "q_net",  qnet_zm, dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean net diabatic heating (cnv+cnd+rad+surf+vdiff+bd)")
        call nc_write(fname, "q_vadv", vadv_zm, dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean vertical advection + adiabatic heating")
        call nc_write(fname, "q_dyn",  dyn_zm,  dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean total dynamical heating (horizontal advection = q_dyn - q_vadv)")
        write(*,"(a)") " rce_long:: wrote zonal-mean RH/cf dump -> "//trim(fname)

        if (ams%diag_acc%n_wacc > 0 .and. allocated(ams%diag_acc%ke_spec)) then
            block
                integer  :: m
                real(wp) :: ket, mft
                ket = sum(ams%diag_acc%ke_spec)/real(ams%diag_acc%n_wacc,wp); mft = sum(ams%diag_acc%mf_spec)/real(ams%diag_acc%n_wacc,wp)
                write(*,"(a)") " eddy zonal-wavenumber spectrum (mass+area wtd, 2nd-half mean):"
                write(*,"(a,es10.3,a,es10.3)") "   total eddy KE ", ket, "   total [u*v*] ", mft
                write(*,"(a)") "     m    KE(m)      KE%     [u*v*](m)   flux%"
                do m = 1, ams%diag_acc%mmax_spec
                    write(*,"(i6,es11.3,f8.1,es12.3,f8.1)") m, ams%diag_acc%ke_spec(m)/real(ams%diag_acc%n_wacc,wp), &
                        100.0_wp*ams%diag_acc%ke_spec(m)/max(sum(ams%diag_acc%ke_spec),tiny(1.0_wp)), &
                        ams%diag_acc%mf_spec(m)/real(ams%diag_acc%n_wacc,wp), &
                        100.0_wp*ams%diag_acc%mf_spec(m)/max(abs(sum(ams%diag_acc%mf_spec)),tiny(1.0_wp))
                end do
            end block
            ! Where the m=8 spike lives vs a transporting mode (m=4): KE by
            ! latitude and by level, % of that wavenumber's total.
            block
                integer  :: j, k
                integer, parameter :: ms(2) = [4, 8]
                integer  :: mi, mm
                real(wp) :: tot
                do mi = 1, 2
                    mm = ms(mi)
                    if (mm > ams%diag_acc%mmax_spec) cycle
                    write(*,"(a,i0,a)") " m=", mm, " eddy KE by latitude (%, N+S folded):"
                    tot = max(sum(ams%diag_acc%ke_latm(:,mm)), tiny(1.0_wp))
                    do j = 1, ams%grd%nlat
                        if (100.0_wp*ams%diag_acc%ke_latm(j,mm)/tot > 3.0_wp) &
                            write(*,"(a,f6.1,a,f6.1)") "    lat ", ams%grd%lat(j), " : ", &
                                100.0_wp*ams%diag_acc%ke_latm(j,mm)/tot
                    end do
                    write(*,"(a,i0,a)") " m=", mm, " eddy KE by level (%, top=1):"
                    tot = max(sum(ams%diag_acc%ke_levm(:,mm)), tiny(1.0_wp))
                    do k = 1, ams%par%nlev
                        write(*,"(a,i3,a,f6.1)") "    lev ", k, " : ", 100.0_wp*ams%diag_acc%ke_levm(k,mm)/tot
                    end do
                end do
            end block
        end if

        deallocate(rh_zm, cf_zm, t_zm, q_zm, p_zm, u_zm, v_zm, cover, latout, levout)
        return
    end subroutine dump_rh

    subroutine dump_monthly(fname)
        ! Write the 12-month zonal-mean climatology (feat/seasonal): T/u/v/RH
        ! (lat,lev,month) and SST/precip (lat,month). Empty months (dacc%mon_cnt=0,
        ! e.g. a sub-annual run) are written as the -9999 missing value.
        character(len=*), intent(in) :: fname
        real(wp), parameter :: MISS = -9999.0_wp
        real(wp), allocatable :: latout(:), lonout(:), levout(:), monout(:)
        real(wp), allocatable :: t_o(:,:,:), u_o(:,:,:), v_o(:,:,:), rh_o(:,:,:)
        real(wp), allocatable :: sst_o(:,:), pr_o(:,:), t2m_o(:,:,:), tlow_o(:,:,:)
        real(wp), allocatable :: pr2d_o(:,:,:), z500_o(:,:,:), vlow_o(:,:,:), evap_o(:,:,:), shf_o(:,:,:)
        real(wp), allocatable :: cf_o(:,:,:), skin_o(:,:,:), lwd_o(:,:,:), olr_o(:,:,:), &
                                 olrcs_o(:,:,:), cover_o(:,:,:)
        real(wp), allocatable :: vt_o(:,:,:), vq_o(:,:,:), eke_o(:,:,:), &
                                 vts_o(:,:,:), vqs_o(:,:,:), ekes_o(:,:,:), uv_o(:,:,:), uvs_o(:,:,:)
        real(wp), allocatable :: mr_o(:,:,:), mrt_o(:,:,:), mrq_o(:,:,:), mr2_o(:,:,:)
        real(wp), allocatable :: swn_o(:,:,:), swup_o(:,:,:)
        real(wp) :: ub(ams%grd%nlon), vb(ams%grd%nlon), tb(ams%grd%nlon), qb(ams%grd%nlon), rn
        integer :: m, jj, kk
        allocate(latout(ams%grd%nlat), lonout(ams%grd%nlon), levout(ams%par%nlev), monout(NMON))
        allocate(t_o(ams%grd%nlat,ams%par%nlev,NMON), u_o(ams%grd%nlat,ams%par%nlev,NMON), &
                 v_o(ams%grd%nlat,ams%par%nlev,NMON), rh_o(ams%grd%nlat,ams%par%nlev,NMON), &
                 sst_o(ams%grd%nlat,NMON), pr_o(ams%grd%nlat,NMON), &
                 t2m_o(ams%grd%nlon,ams%grd%nlat,NMON), tlow_o(ams%grd%nlon,ams%grd%nlat,NMON), pr2d_o(ams%grd%nlon,ams%grd%nlat,NMON), &
                 z500_o(ams%grd%nlon,ams%grd%nlat,NMON), vlow_o(ams%grd%nlon,ams%grd%nlat,NMON), &
                 evap_o(ams%grd%nlon,ams%grd%nlat,NMON), shf_o(ams%grd%nlon,ams%grd%nlat,NMON))
        allocate(mr_o(ams%grd%nlat,ams%par%nlev,NMON), mrt_o(ams%grd%nlat,ams%par%nlev,NMON), mrq_o(ams%grd%nlat,ams%par%nlev,NMON))
        allocate(mr2_o(ams%grd%nlon,ams%grd%nlat,NMON))
        allocate(swn_o(ams%grd%nlon,ams%grd%nlat,NMON), swup_o(ams%grd%nlon,ams%grd%nlat,NMON))
        allocate(vt_o(ams%grd%nlat,ams%par%nlev,NMON), vq_o(ams%grd%nlat,ams%par%nlev,NMON), eke_o(ams%grd%nlat,ams%par%nlev,NMON), &
                 vts_o(ams%grd%nlat,ams%par%nlev,NMON), vqs_o(ams%grd%nlat,ams%par%nlev,NMON), ekes_o(ams%grd%nlat,ams%par%nlev,NMON), &
                 uv_o(ams%grd%nlat,ams%par%nlev,NMON), uvs_o(ams%grd%nlat,ams%par%nlev,NMON))
        allocate(cf_o(ams%grd%nlat,ams%par%nlev,NMON), skin_o(ams%grd%nlon,ams%grd%nlat,NMON), &
                 lwd_o(ams%grd%nlon,ams%grd%nlat,NMON), olr_o(ams%grd%nlon,ams%grd%nlat,NMON), &
                 olrcs_o(ams%grd%nlon,ams%grd%nlat,NMON), cover_o(ams%grd%nlon,ams%grd%nlat,NMON))
        latout = ams%grd%lat(1:ams%grd%nlat)
        lonout = ams%grd%lon(1:ams%grd%nlon)
        levout = ams%vgrid%sigma_full(1:ams%par%nlev)
        do m = 1, NMON
            monout(m) = real(m, wp)
            if (ams%diag_acc%mon_cnt(m) > 0) then
                t_o(:,:,m)  = ams%diag_acc%mon_t(:,:,m) /real(ams%diag_acc%mon_cnt(m),wp)
                u_o(:,:,m)  = ams%diag_acc%mon_u(:,:,m) /real(ams%diag_acc%mon_cnt(m),wp)
                v_o(:,:,m)  = ams%diag_acc%mon_v(:,:,m) /real(ams%diag_acc%mon_cnt(m),wp)
                rh_o(:,:,m) = ams%diag_acc%mon_rh(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                sst_o(:,m)  = ams%diag_acc%mon_sst(:,m) /real(ams%diag_acc%mon_cnt(m),wp)
                pr_o(:,m)   = ams%diag_acc%mon_pr(:,m)  /real(ams%diag_acc%mon_cnt(m),wp)
                t2m_o(:,:,m) = ams%diag_acc%mon_t2m(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                tlow_o(:,:,m) = ams%diag_acc%mon_tlow(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                pr2d_o(:,:,m) = ams%diag_acc%mon_pr2d(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                z500_o(:,:,m) = ams%diag_acc%mon_z500(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                vlow_o(:,:,m) = ams%diag_acc%mon_vlow(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                evap_o(:,:,m) = ams%diag_acc%mon_evap(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                shf_o(:,:,m)  = ams%diag_acc%mon_shf(:,:,m) /real(ams%diag_acc%mon_cnt(m),wp)
                cf_o(:,:,m)   = ams%diag_acc%mon_cf(:,:,m)  /real(ams%diag_acc%mon_cnt(m),wp)
                skin_o(:,:,m) = ams%diag_acc%mon_skin(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                lwd_o(:,:,m)  = ams%diag_acc%mon_lwd(:,:,m) /real(ams%diag_acc%mon_cnt(m),wp)
                olr_o(:,:,m)  = ams%diag_acc%mon_olr(:,:,m) /real(ams%diag_acc%mon_cnt(m),wp)
                olrcs_o(:,:,m) = ams%diag_acc%mon_olrcs(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                cover_o(:,:,m) = ams%diag_acc%mon_cover(:,:,m)/real(ams%diag_acc%mon_cnt(m),wp)
                rn = real(ams%diag_acc%mon_cnt(m),wp)
                vt_o(:,:,m)  = ams%diag_acc%mon_vt(:,:,m) /rn
                uv_o(:,:,m)  = ams%diag_acc%mon_uv(:,:,m) /rn
                vq_o(:,:,m)  = ams%diag_acc%mon_vq(:,:,m) /rn*1000.0_wp
                eke_o(:,:,m) = ams%diag_acc%mon_eke(:,:,m)/rn
                mr_o(:,:,m)  = ams%diag_acc%mon_mr(:,:,m) /rn
                mrt_o(:,:,m) = ams%diag_acc%mon_mrt(:,:,m)/rn
                mrq_o(:,:,m) = ams%diag_acc%mon_mrq(:,:,m)/rn*1000.0_wp
                mr2_o(:,:,m) = ams%diag_acc%mon_mr2(:,:,m)/rn
                swn_o(:,:,m) = ams%diag_acc%mon_swn(:,:,m)/rn; swup_o(:,:,m) = ams%diag_acc%mon_swup(:,:,m)/rn
                ! stationary part from the monthly-mean fields: [v*T*] etc.
                do kk = 1, ams%par%nlev
                    do jj = 1, ams%grd%nlat
                        ub = ams%diag_acc%mon_u3(:,jj,kk,m)/rn; vb = ams%diag_acc%mon_v3(:,jj,kk,m)/rn
                        tb = ams%diag_acc%mon_t3(:,jj,kk,m)/rn; qb = ams%diag_acc%mon_q3(:,jj,kk,m)/rn
                        ub = ub - sum(ub)/real(ams%grd%nlon,wp); vb = vb - sum(vb)/real(ams%grd%nlon,wp)
                        tb = tb - sum(tb)/real(ams%grd%nlon,wp); qb = qb - sum(qb)/real(ams%grd%nlon,wp)
                        vts_o(jj,kk,m)  = sum(vb*tb)/real(ams%grd%nlon,wp)
                        uvs_o(jj,kk,m)  = sum(ub*vb)/real(ams%grd%nlon,wp)
                        vqs_o(jj,kk,m)  = sum(vb*qb)/real(ams%grd%nlon,wp)*1000.0_wp
                        ekes_o(jj,kk,m) = 0.5_wp*sum(ub*ub + vb*vb)/real(ams%grd%nlon,wp)
                    end do
                end do
            else
                vt_o(:,:,m) = MISS; vq_o(:,:,m) = MISS; eke_o(:,:,m) = MISS
                mr_o(:,:,m) = MISS; mrt_o(:,:,m) = MISS; mrq_o(:,:,m) = MISS; mr2_o(:,:,m) = MISS
                swn_o(:,:,m) = MISS; swup_o(:,:,m) = MISS
                vts_o(:,:,m) = MISS; vqs_o(:,:,m) = MISS; ekes_o(:,:,m) = MISS
                uv_o(:,:,m) = MISS; uvs_o(:,:,m) = MISS
                vlow_o(:,:,m) = MISS; evap_o(:,:,m) = MISS; shf_o(:,:,m) = MISS
                cf_o(:,:,m) = MISS; skin_o(:,:,m) = MISS; lwd_o(:,:,m) = MISS
                olr_o(:,:,m) = MISS; olrcs_o(:,:,m) = MISS; cover_o(:,:,m) = MISS
                t_o(:,:,m) = MISS; u_o(:,:,m) = MISS; v_o(:,:,m) = MISS
                rh_o(:,:,m) = MISS; sst_o(:,m) = MISS; pr_o(:,m) = MISS
                t2m_o(:,:,m) = MISS; tlow_o(:,:,m) = MISS; pr2d_o(:,:,m) = MISS; z500_o(:,:,m) = MISS
            end if
        end do
        call nc_create(fname)
        call nc_write_dim(fname, "lat",   x=latout, units="degrees_north")
        call nc_write_dim(fname, "lon",   x=lonout, units="degrees_east")
        call nc_write_dim(fname, "lev",   x=levout, units="sigma")
        call nc_write_dim(fname, "month", x=monout, units="1")
        call nc_write(fname, "t",   t_o,  dim1="lat", dim2="lev", dim3="month", &
            units="K",  long_name="monthly zonal-mean temperature")
        call nc_write(fname, "u",   u_o,  dim1="lat", dim2="lev", dim3="month", &
            units="m/s", long_name="monthly zonal-mean zonal wind")
        call nc_write(fname, "v",   v_o,  dim1="lat", dim2="lev", dim3="month", &
            units="m/s", long_name="monthly zonal-mean meridional wind")
        call nc_write(fname, "rh",  rh_o, dim1="lat", dim2="lev", dim3="month", &
            units="%",  long_name="monthly zonal-mean relative humidity")
        call nc_write(fname, "sst", sst_o, dim1="lat", dim2="month", &
            units="K",  long_name="monthly zonal-mean surface temperature")
        call nc_write(fname, "precip", pr_o, dim1="lat", dim2="month", &
            units="mm/day", long_name="monthly zonal-mean precipitation")
        call nc_write(fname, "t2m", t2m_o, dim1="lon", dim2="lat", dim3="month", &
            units="K", long_name="monthly near-surface (lowest level) air temperature")
        call nc_write(fname, "tlow", tlow_o, dim1="lon", dim2="lat", dim3="month", &
            units="K", long_name="monthly lowest-level air temperature (~350 m)")
        call nc_write(fname, "precip2d", pr2d_o, dim1="lon", dim2="lat", dim3="month", &
            units="mm/day", long_name="monthly precipitation")
        call nc_write(fname, "z500", z500_o, dim1="lon", dim2="lat", dim3="month", &
            units="m", long_name="monthly 500 hPa geopotential height (hydrostatic, log-p interpolated)")
        call nc_write(fname, "vlow", vlow_o, dim1="lon", dim2="lat", dim3="month", &
            units="m/s", long_name="monthly lowest-level wind speed")
        call nc_write(fname, "evap", evap_o, dim1="lon", dim2="lat", dim3="month", &
            units="mm/day", long_name="monthly surface evaporation")
        call nc_write(fname, "shf", shf_o, dim1="lon", dim2="lat", dim3="month", &
            units="W/m2", long_name="monthly surface sensible heat flux (up)")
        call nc_write(fname, "cf", cf_o, dim1="lat", dim2="lev", dim3="month", &
            units="1", long_name="monthly zonal-mean layer cloud fraction (as seen by radiation)")
        call nc_write(fname, "skin", skin_o, dim1="lon", dim2="lat", dim3="month", &
            units="K", long_name="monthly skin temperature (land skin over land, SST over ocean)")
        call nc_write(fname, "lwd", lwd_o, dim1="lon", dim2="lat", dim3="month", &
            units="W/m2", long_name="monthly surface downward longwave")
        call nc_write(fname, "olr", olr_o, dim1="lon", dim2="lat", dim3="month", &
            units="W/m2", long_name="monthly outgoing longwave at TOA")
        call nc_write(fname, "olr_cs", olrcs_o, dim1="lon", dim2="lat", dim3="month", &
            units="W/m2", long_name="monthly clear-sky outgoing longwave at TOA")
        call nc_write(fname, "cover", cover_o, dim1="lon", dim2="lat", dim3="month", &
            units="1", long_name="monthly column cloud cover (max overlap)")
        call nc_write(fname, "vt", vt_o, dim1="lat", dim2="lev", dim3="month", &
            units="K m/s", long_name="monthly zonal-mean eddy heat flux [v'T'] (stationary + transient)")
        call nc_write(fname, "vq", vq_o, dim1="lat", dim2="lev", dim3="month", &
            units="g/kg m/s", long_name="monthly zonal-mean eddy moisture flux [v'q'] (stationary + transient)")
        call nc_write(fname, "eke", eke_o, dim1="lat", dim2="lev", dim3="month", &
            units="m2/s2", long_name="monthly zonal-mean eddy kinetic energy (stationary + transient)")
        call nc_write(fname, "uv", uv_o, dim1="lat", dim2="lev", dim3="month", &
                      units="m2 s-2", long_name="eddy momentum flux [u'v'] (transient+stationary), zonal mean")
        call nc_write(fname, "uv_stat", uvs_o, dim1="lat", dim2="lev", dim3="month", &
                      units="m2 s-2", long_name="stationary eddy momentum flux [u*v*] of the monthly means")
        call nc_write(fname, "vt_stat", vts_o, dim1="lat", dim2="lev", dim3="month", &
            units="K m/s", long_name="monthly stationary eddy heat flux [v*T*] of the monthly-mean fields")
        call nc_write(fname, "vq_stat", vqs_o, dim1="lat", dim2="lev", dim3="month", &
            units="g/kg m/s", long_name="monthly stationary eddy moisture flux [v*q*] of the monthly-mean fields")
        call nc_write(fname, "eke_stat", ekes_o, dim1="lat", dim2="lev", dim3="month", &
            units="m2/s2", long_name="monthly stationary eddy kinetic energy of the monthly-mean fields")
        call nc_write(fname, "mr", mr_o, dim1="lat", dim2="lev", dim3="month", &
            units="1", long_name="monthly zonal-mean tracer mass-rescaling term dp_FV/dp_spec - 1 (per step)")
        call nc_write(fname, "mrt", mrt_o, dim1="lat", dim2="lev", dim3="month", &
            units="K", long_name="monthly zonal-mean eddy covariance [r' T'] of the rescaling term with temperature")
        call nc_write(fname, "mrq", mrq_o, dim1="lat", dim2="lev", dim3="month", &
            units="g/kg", long_name="monthly zonal-mean eddy covariance [r' q'] of the rescaling term with humidity")
        call nc_write(fname, "mr2d", mr2_o, dim1="lon", dim2="lat", dim3="month", &
            units="1", long_name="monthly-mean tracer mass-rescaling term dp_FV/dp_spec - 1 at the lowest level")
        call nc_write(fname, "swn", swn_o, dim1="lon", dim2="lat", dim3="month", &
            units="W/m2", long_name="monthly surface net absorbed shortwave")
        call nc_write(fname, "swup", swup_o, dim1="lon", dim2="lat", dim3="month", &
            units="W/m2", long_name="monthly reflected shortwave at TOA")
        write(*,"(a,i0,a,a)") " rce_long:: wrote monthly climatology (", &
            count(ams%diag_acc%mon_cnt>0), " months populated) -> ", trim(fname)
        deallocate(latout, lonout, levout, monout, t_o, u_o, v_o, rh_o, sst_o, pr_o, t2m_o, tlow_o, &
                   pr2d_o, z500_o, vlow_o, evap_o, shf_o, cf_o, skin_o, lwd_o, olr_o, &
                   olrcs_o, cover_o, vt_o, vq_o, eke_o, vts_o, vqs_o, ekes_o, mr_o, mrt_o, mrq_o, mr2_o, swn_o, swup_o, &
                   uv_o, uvs_o)
    end subroutine dump_monthly
    subroutine read_config()
        ! The driver's own `&rce` keys. Every one is optional in the case file
        ! and falls back to input/aeros_defaults.nml (&rce there is the schema).
        ! The model's own configuration is read by aeros_lib_init.
        logical :: ex
        inquire(file=trim(nmlfile), exist=ex)
        if (.not. ex) then
            write(*,"(a)") " rce_long:: error: namelist file not found: "//trim(nmlfile)
            error stop 1
        end if
        call nml_read(nmlfile, "rce", "nstep",             nstep,             defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "print_every",       print_every,       defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "l_diag",            l_diag,            defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "monthly_out",       monthly_out,       defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "mon_spinup_days",   mon_spinup_days,   defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "restart_in",        restart_in,        defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "restart_out",       restart_out,       defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "restart_interval",  restart_interval,  defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "l_uniform_insol",   l_uniform_insol,   defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "l_qforce",          l_qforce,          defaults_file=defaults_file)
        call nml_read(nmlfile, "rce", "qforce_amp",        qforce_amp,        defaults_file=defaults_file)
        return
    end subroutine read_config

    subroutine report(n)
        integer, intent(in) :: n
        real(wp) :: tmx, tmn, tmean, umax
        real(wp) :: hphys, hcnd
        integer  :: imx, jmx, imn, jmn, kum, jum, ium
        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)
        ! max|u| with its location: if it lives at the top level and migrates
        ! poleward as it grows, the terminal blow-up is a top thermal-wind jet.
        call locate_umax(umax, ium, jum, kum)
        write(*,"(a)") ""
        write(*,"(a,i0,a,f8.2,a,f7.1,a,i0,a,f6.1,a)") " -- step ", n, "  day ", &
            real(n,wp)*ams%par%dt/86400.0_wp, "   max|u| ", umax, " m/s @ lev ", kum, &
            " lat ", ams%grd%lat(jum), " N"
        write(*,"(a)") "   lev    Tmean     Tmax  @(i,j)        Tmin    " // &
                       "Hfwd[K/d] Hcnd[K/d]"
        do k = 1, ams%par%nlev
            call level_stats(k, tmean, tmx, imx, jmx, tmn, imn, jmn)
            ! forward-split (surf+cnv+rad) and centered (cnd) heating, K/day
            hphys = gmean_lev(ams%ts%wrk%dt_phys, k)/ams%par%dt*86400.0_wp
            hcnd  = gmean_lev(ams%ts%wrk%dtdt, k)*86400.0_wp
            write(*,"(i6,f9.2,f9.2,a,i0,a,i0,a,f9.2,f10.3,f10.3)") &
                k, tmean, tmx, "  (", imx, ",", jmx, ")", tmn, hphys, hcnd
        end do
        call hotspot_split(ams%par%nlev)
        call eddy_diag()
        call baroclinicity_diag()
        if (l_diag) then
            call term_table(hot_lat_index())          ! tropical (warmest low level)
            call term_table(lat_index_near(45.0_wp))  ! midlatitude
        end if
        call energy_balance()
        if (ams%ts%cpr%enabled) call cloud_prog_cover()
        return
    end subroutine report

    subroutine cloud_prog_cover()
        ! Area-weighted global-mean max-overlap total cloud cover from the
        ! PROGNOSTIC cloud fraction (now%cf_g), plus its column max. This is the
        ! quantity the diagnostic scheme runs away on (0.66 -> 0.86); the point
        ! of the prognostic budget is that it does not.
        real(wp) :: colcov, wsum, cov, w, cfmax
        integer  :: ii, jj, kk
        cov = 0.0_wp; wsum = 0.0_wp; cfmax = 0.0_wp
        do jj = 1, ams%grd%nlat
            do ii = 1, ams%grd%nlon
                colcov = 0.0_wp
                do kk = 1, ams%par%nlev
                    colcov = max(colcov, ams%now%cf_g(ii,jj,kk))
                end do
                w = ams%grd%area(ii,jj)
                cov  = cov + colcov*w
                wsum = wsum + w
                cfmax = max(cfmax, colcov)
            end do
        end do
        write(*,"(a,f7.3,a,f7.3)") "   prognostic cloud cover  mean ", &
            cov/wsum, "   max-column ", cfmax
        return
    end subroutine cloud_prog_cover

    integer function hot_lat_index() result(jstar)
        ! Latitude row with the warmest zonal-mean lowest layer (the tropics).
        real(wp) :: zm, zmax
        integer  :: j
        zmax = -1.0e30_wp; jstar = 1
        do j = 1, ams%grd%nlat
            zm = sum(ams%now%temp_g(:,j,ams%par%nlev))/real(ams%grd%nlon, wp)
            if (zm > zmax) then; zmax = zm; jstar = j; end if
        end do
        return
    end function hot_lat_index

    integer function lat_index_near(target) result(jbest)
        ! Latitude row nearest a target latitude [deg].
        real(wp), intent(in) :: target
        real(wp) :: d, dbest
        integer  :: j
        dbest = 1.0e30_wp; jbest = 1
        do j = 1, ams%grd%nlat
            d = abs(ams%grd%lat(j) - target)
            if (d < dbest) then; dbest = d; jbest = j; end if
        end do
        return
    end function lat_index_near

    subroutine term_table(jstar)
        ! Per-term heating split, zonal-mean at latitude index jstar. Called for
        ! the tropical hot latitude (warmest lowest layer) AND a midlatitude
        ! (~45 deg): the equator-to-midlat contrast in the cnv/cnd heating aloft
        ! reveals whether moist physics (convection, condensation latent heat) is
        ! warming the midlatitude free troposphere onto a moist adiabat and so
        ! collapsing the upper-level meridional temperature gradient the
        ! thermal-wind jet needs (m2_results §26). All columns in K/day; hadv =
        ! total dynamical heating minus its vertical/adiabatic part (vadv).
        integer, intent(in) :: jstar
        real(wp) :: zmax
        real(wp) :: hs, hcv, hcd, hr, hvd, hva, hdyn
        real(wp) :: olrj, swupj, swinj
        integer  :: k
        zmax = sum(ams%now%temp_g(:,jstar,ams%par%nlev))/real(ams%grd%nlon, wp)
        write(*,"(a,i0,a,f6.1,a,f7.2,a)") &
            "   per-term heating [K/day] at j=", jstar, &
            " (lat ", ams%grd%lat(jstar), " N, T_low ", zmax, " K):"
        write(*,"(a)") "   lev    Tzm    qzm[g/kg]    surf     cnv     cnd" // &
                       "     rad   vdiff    vadv    hadv"
        do k = 1, ams%par%nlev
            hs  = zmean_at(ams%ts%wrk%dt_surf,  k, jstar)/ams%par%dt*86400.0_wp
            hcv = zmean_at(ams%ts%wrk%dt_cnv,   k, jstar)/ams%par%dt*86400.0_wp
            hcd = zmean_at(ams%ts%wrk%dt_cnd,   k, jstar)/ams%par%dt*86400.0_wp
            hr  = zmean_at(ams%ts%wrk%dt_rad,   k, jstar)/ams%par%dt*86400.0_wp
            hvd = zmean_at(ams%ts%wrk%dt_vdiff, k, jstar)/ams%par%dt*86400.0_wp
            hva = zmean_at(ams%ts%wrk%dt_vadv,  k, jstar)*86400.0_wp
            hdyn= zmean_at(ams%ts%wrk%dtdt,     k, jstar)*86400.0_wp
            ! Tzm/qzm: the zonal-mean T and q at j* that RADIATION sees -- a
            ! vertical sawtooth here is the garbage-in behind the huge rad column.
            write(*,"(i6,f8.2,f10.4,3x,7f8.2)") k, &
                sum(ams%now%temp_g(:,jstar,k))/real(ams%grd%nlon,wp), &
                sum(ams%now%qv_g(:,jstar,k))/real(ams%grd%nlon,wp)*1000.0_wp, &
                hs, hcv, hcd, hr, hvd, hva, hdyn - hva
        end do
        ! LOCAL TOA at the hot latitude: does OLR rise as T_low climbs, or does it
        ! saturate (runaway-greenhouse signature -- absorbed SW can't be shed)?
        if (ams%phys%l_rad) then
            olrj  = sum(ams%ts%rad%olr(:,jstar))/real(ams%grd%nlon,wp)
            swupj = sum(ams%ts%rad%sw_up_toa(:,jstar))/real(ams%grd%nlon,wp)
            swinj = ams%ts%rad%sw_toa(jstar)
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
        real(wp) :: phalf(0:ams%par%nlev), pfull(ams%par%nlev), dpc(ams%par%nlev), phi(ams%par%nlev)
        real(wp) :: tzm(ams%par%nlev), qzm(ams%par%nlev), qs, dqs, hb, hstar, rh, psj, cp
        integer  :: k
        cp = real(cp_d, wp)
        psj = sum(ams%now%ps(:,jstar))/real(ams%grd%nlon, wp)
        call aeros_vgrid_pressure(ams%vgrid, psj, phalf, pfull, dpc)
        do k = 1, ams%par%nlev
            tzm(k) = sum(ams%now%temp_g(:,jstar,k))/real(ams%grd%nlon, wp)
            qzm(k) = sum(ams%now%qv_g(:,jstar,k))/real(ams%grd%nlon, wp)
        end do
        phi(ams%par%nlev) = 0.0_wp
        do k = ams%par%nlev-1, 1, -1
            phi(k) = phi(k+1) + R_d*0.5_wp*(tzm(k)+tzm(k+1))*log(pfull(k+1)/pfull(k))
        end do
        hb = cp*tzm(ams%par%nlev) + phi(ams%par%nlev) + aeros_latent(tzm(ams%par%nlev))*qzm(ams%par%nlev)
        write(*,"(a,f8.1,a)") "   BL parcel h_b = ", hb/1000.0_wp, &
            " kJ/kg;  buoyancy hb-h*_env (>0 = convecting band):"
        write(*,"(a)") "   lev   RH[%]   h*_env[kJ/kg]   hb-h*[kJ/kg]"
        do k = 1, ams%par%nlev
            call aeros_qsat(tzm(k), pfull(k), qs, dqs)
            rh = qzm(k)/max(qs, 1.0e-12_wp)*100.0_wp
            hstar = cp*tzm(k) + phi(k) + aeros_latent(tzm(k))*qs
            write(*,"(i6,f8.1,f15.1,f14.2)") k, rh, hstar/1000.0_wp, (hb-hstar)/1000.0_wp
        end do
        return
    end subroutine buoyancy_profile

    real(wp) function zmean_at(f, k, j) result(m)
        ! Zonal mean of grid field f at level k, latitude row j.
        real(wp), intent(in) :: f(:,:,:)
        integer,  intent(in) :: k, j
        m = sum(f(:,j,k))/real(ams%grd%nlon, wp)
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
        do k = 1, ams%par%nlev
            tmax = -1.0e30_wp; tmin = 1.0e30_wp; upk = 0.0_wp
            do j = 1, ams%grd%nlat
                tzm = sum(ams%now%temp_g(:,j,k))/real(ams%grd%nlon, wp)
                uzm = sum(ams%now%u(:,j,k))/real(ams%grd%nlon, wp)
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
        k = max(1, ams%par%nlev/2)                 ! ~mid-troposphere
        t2 = 0.0_wp; ke = 0.0_wp; vt = 0.0_wp; wsum = 0.0_wp
        do j = 1, ams%grd%nlat
            tzm = sum(ams%now%temp_g(:,j,k))/real(ams%grd%nlon, wp)
            uzm = sum(ams%now%u(:,j,k))/real(ams%grd%nlon, wp)
            vzm = sum(ams%now%v(:,j,k))/real(ams%grd%nlon, wp)
            w   = ams%grd%area(1,j)
            do i = 1, ams%grd%nlon
                tp = ams%now%temp_g(i,j,k) - tzm
                up = ams%now%u(i,j,k) - uzm
                vp = ams%now%v(i,j,k) - vzm
                t2 = t2 + w*tp*tp
                ke = ke + w*0.5_wp*(up*up + vp*vp)
                vt = vt + w*vp*tp
            end do
            wsum = wsum + w*real(ams%grd%nlon, wp)
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
        do k = 1, ams%par%nlev
            do j = 1, ams%grd%nlat
                do i = 1, ams%grd%nlon
                    a = abs(ams%now%u(i,j,k))
                    if (a > umax) then; umax = a; ium = i; jum = j; kum = k; end if
                end do
            end do
        end do
        return
    end subroutine locate_umax

    subroutine energy_balance()
        ! Global-mean TOA and surface energy budget, time-mean since the last
        ! report. A secular drift of the column temperature is a TOA imbalance
        ! (absorbed SW - OLR /= 0); this separates an energy-balance problem from
        ! a mixing/numerics one.
        real(wp) :: olr, swup, swin, shf, lhf, olrc, swupc, cov, rn
        if (.not. (ams%phys%l_rad .or. ams%phys%l_surf)) return
        if (ams%diag_acc%tb_n == 0) return
        rn   = 1.0_wp/real(ams%diag_acc%tb_n, wp)
        olr  = ams%diag_acc%tb_olr*rn;  swup  = ams%diag_acc%tb_swup*rn;  swin = ams%diag_acc%tb_swin*rn
        olrc = ams%diag_acc%tb_olrc*rn; swupc = ams%diag_acc%tb_swupc*rn; cov  = ams%diag_acc%tb_cov*rn
        shf  = ams%diag_acc%tb_shf*rn;  lhf   = ams%diag_acc%tb_lhf*rn
        ams%diag_acc%tb_olr = 0.0_wp; ams%diag_acc%tb_swup = 0.0_wp; ams%diag_acc%tb_swin = 0.0_wp; ams%diag_acc%tb_olrc = 0.0_wp
        ams%diag_acc%tb_swupc = 0.0_wp; ams%diag_acc%tb_cov = 0.0_wp; ams%diag_acc%tb_shf = 0.0_wp; ams%diag_acc%tb_lhf = 0.0_wp; ams%diag_acc%tb_n = 0
        write(*,"(a,f7.1,a,f7.1,a,f7.2,a,f7.1,a,f7.1,a)") &
            "   TOA: SWin ", swin, "  OLR ", olr, "  net ", swin-swup-olr, &
            " W/m2 | sfc SH ", shf, "  LH ", lhf, " W/m2  (mean since last report)"
        ! Cloud radiative effect = all-sky minus clear-sky (ERA5 1991-2020: absorbed SW
        ! 243, OLR 242, clear 289 / 264, CRE SW -46 LW +22, total cover 0.63).
        if (ams%phys%l_rad .and. ams%ts%rad%clouds) write(*,"(a,f7.1,a,f7.1,a,f6.1,a,f6.1,a,f5.2)") &
            "   TOA clear: SWabs ", swin-swupc, "  OLR ", olrc, " | CRE SW ", swupc-swup, &
            "  LW ", olrc-olr, " W/m2 | cloud cover ", cov
        call water_balance()
        return
    end subroutine energy_balance

    subroutine water_balance()
        ! Global water budget since the last report: precipitable water W and
        ! its change rate against the mean E - P; the residual is water the
        ! model created or destroyed (should be ~0).
        real(wp) :: wcol(ams%grd%nlon, ams%grd%nlat), phw(0:ams%par%nlev), pfw(ams%par%nlev), dpw(ams%par%nlev)
        real(wp) :: w, pm, em, dwdt
        integer  :: i, j, k
        if (ams%diag_acc%wb_n_acc == 0) return
        do j = 1, ams%grd%nlat
            do i = 1, ams%grd%nlon
                call aeros_vgrid_pressure(ams%vgrid, ams%now%ps(i,j), phw, pfw, dpw)
                wcol(i,j) = 0.0_wp
                do k = 1, ams%par%nlev
                    wcol(i,j) = wcol(i,j) + ams%now%qv_g(i,j,k)*dpw(k)/real(grav,wp)
                end do
            end do
        end do
        w  = gmean2(wcol)
        pm = ams%diag_acc%wb_p_acc/real(ams%diag_acc%wb_n_acc,wp); em = ams%diag_acc%wb_e_acc/real(ams%diag_acc%wb_n_acc,wp)
        if (ams%diag_acc%wb_w_last >= 0.0_wp) then
            dwdt = (w - ams%diag_acc%wb_w_last)/(real(ams%diag_acc%wb_n_acc,wp)*ams%par%dt)*86400.0_wp     ! kg/m2/day = mm/day
            write(*,"(a,f7.2,a,f6.2,a,f6.2,a,f7.3,a,f7.3,a)") &
                "   WATER: W ", w, " mm | P ", pm, "  E ", em, "  dW/dt ", dwdt, &
                "  resid (dW/dt - E + P) ", dwdt - em + pm, " mm/day"
        else
            write(*,"(a,f7.2,a,f6.2,a,f6.2,a)") "   WATER: W ", w, " mm | P ", pm, "  E ", em, " mm/day"
        end if
        if (l_diag) then
            write(*,"(a,6f8.3,a)") "   WATER by operator (surf cnv cnd cld vdiff transp): ", &
                ams%ts%wb*86400.0_dp/(real(ams%diag_acc%wb_n_acc,dp)*real(ams%par%dt,dp)), " mm/day"
            ams%ts%wb = 0.0_dp
        end if
        ams%diag_acc%wb_w_last = w; ams%diag_acc%wb_p_acc = 0.0_wp; ams%diag_acc%wb_e_acc = 0.0_wp; ams%diag_acc%wb_n_acc = 0
        return
    end subroutine water_balance

    real(wp) function gmean2(f) result(m)
        real(wp), intent(in) :: f(:,:)
        real(wp) :: s, wsum
        integer  :: i, j
        s = 0.0_wp; wsum = 0.0_wp
        do j = 1, ams%grd%nlat
            do i = 1, ams%grd%nlon
                s = s + f(i,j)*ams%grd%area(i,j); wsum = wsum + ams%grd%area(i,j)
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
        real(wp) :: zm(ams%grd%nlat), zmax, devmax, d
        integer  :: i, j, jzmax, idev, jdev
        do j = 1, ams%grd%nlat
            zm(j) = sum(ams%now%temp_g(:,j,k))/real(ams%grd%nlon, wp)
        end do
        zmax = -1.0e30_wp; jzmax = 1
        do j = 1, ams%grd%nlat
            if (zm(j) > zmax) then; zmax = zm(j); jzmax = j; end if
        end do
        devmax = 0.0_wp; idev = 1; jdev = 1
        do j = 1, ams%grd%nlat
            do i = 1, ams%grd%nlon
                d = abs(ams%now%temp_g(i,j,k) - zm(j))
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
        do j = 1, ams%grd%nlat
            w = ams%grd%area(1,j)
            do i = 1, ams%grd%nlon
                if (ams%now%temp_g(i,j,k) > tmx) then
                    tmx = ams%now%temp_g(i,j,k); imx = i; jmx = j
                end if
                if (ams%now%temp_g(i,j,k) < tmn) then
                    tmn = ams%now%temp_g(i,j,k); imn = i; jmn = j
                end if
                s = s + ams%now%temp_g(i,j,k)*ams%grd%area(i,j); wsum = wsum + ams%grd%area(i,j)
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
        do j = 1, ams%grd%nlat
            do i = 1, ams%grd%nlon
                s = s + f(i,j,k)*ams%grd%area(i,j); wsum = wsum + ams%grd%area(i,j)
            end do
        end do
        m = s/wsum
        return
    end function gmean_lev

    subroutine locate_nan()
        integer :: i, j, k
        do k = 1, ams%par%nlev
            do j = 1, ams%grd%nlat
                do i = 1, ams%grd%nlon
                    if (ams%ts%wrk%t_g(i,j,k) /= ams%ts%wrk%t_g(i,j,k)) then
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
