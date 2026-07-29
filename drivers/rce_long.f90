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
    use aeros_cloud,    only : aeros_cloud_diagnose
    use aeros_convection, only : SCHEME_SBM, SCHEME_SBM_FRIERSON, SCHEME_MANABE
    use aeros_timestep
    use aeros_radiation, only : SCHEME_ECCKD
    use aeros_ocean,    only : aeros_ocean_init
    use aeros_topography, only : aeros_topography_load, aeros_topography_scale
    use aeros_land,     only : aeros_land_init, aeros_land_couple_radiation
    use nml,            only : nml_read
    use ncio,           only : nc_create, nc_write_dim, nc_write

    implicit none

    ! --- namelist-configurable, with test_rce defaults -----------------------
    character(len=512) :: nmlfile
    ! End-of-run zonal-mean RH/cf dump (arg 2; "" = off). NOT a namelist key --
    ! kept off the namelist so it needs no entry in every rce_*.nml (nml_read is
    ! fatal on a missing key in this driver). The comparison against ERA5 is done
    ! offline by scripts/rce_humidity_vs_era5.jl.
    character(len=512) :: rh_out = "output/rce_rh.nc"
    integer  :: trunc = 21, nlev = 12, nstep = 4800, ndiff = 6
    integer  :: print_every = 96          ! ~2 model days at dt=1800
    real(wp) :: dt = 1800.0_wp, tau_diff = 6.0_wp
    real(wp) :: tau_diff_div = 0.0_wp     ! divergence-only diffusion timescale [h]; 0 => = tau_diff
    real(wp) :: t_ref = 300.0_wp          ! isothermal semi-implicit reference [K]
    logical  :: vert_vanleer = .FALSE.    ! van Leer (not donor-cell) vertical q transport
    ! Couple diabatic heating in-solve (step 1b) vs the default forward-split.
    ! Needs eps_filter ~0.15 to hold the convective computational mode.
    logical  :: couple_diabatic = .FALSE.
    real(wp) :: eps_filter = 0.06_wp, raw_alpha = 0.53_wp
    real(wp) :: si_alpha = 0.5_wp         ! semi-implicit decentering: 0.5 centered, 1.0 backward (SW)
    real(wp) :: conv_tau = 7200.0_wp, c_h = 1.5e-3_wp, c_e = 1.5e-3_wp, u_min = 1.0_wp
    character(len=32) :: conv_scheme = "sbm_frierson" ! "sbm" | "sbm_frierson" | "manabe"
    real(wp) :: c_d = 1.5e-3_wp       ! surface momentum drag (0 = off; brakes the low-level jet)
    logical  :: l_surf = .TRUE., l_cnv = .TRUE., l_cnd = .TRUE.
    logical  :: l_rad = .TRUE., l_sponge = .TRUE., l_vdiff = .FALSE.
    logical  :: l_diag = .TRUE.       ! per-term heating split at the hot latitude
    logical  :: l_dry_adjust = .TRUE. ! dry convective adjustment before the moist scheme
    logical  :: l_rad_clouds = .FALSE. ! diagnostic all-sky clouds in the radiation
    logical  :: l_uniform_insol = .FALSE. ! flatten insolation to its global mean (no meridional gradient)
    logical  :: l_nonrotating = .FALSE.   ! zero Coriolis (no jet organization; RCE vehicle)
    ! Diagnostic prescribed heating (dry-core heating->ω isolation). When on, set
    ! all physics off; an analytic Q(lat,sigma) = amp * [equatorial, zero-mean] *
    ! sin(pi*sigma) is coupled on the in-solve seam. amp in K/day.
    logical  :: l_qforce = .FALSE.
    real(wp) :: qforce_amp = 2.0_wp       ! prescribed heating amplitude [K/day]
    real(wp) :: vdiff_k0 = 10.0_wp, vdiff_sigma = 0.7_wp
    logical  :: vdiff_richardson = .TRUE.    ! Ri-diagnosed BL depth (else fixed)
    real(wp) :: vdiff_ri_crit = 10.0_wp      ! critical bulk Richardson number
    ! Model-top sponge knobs (defaults match aeros_timestep_class). Exposed to
    ! test the top thermal-wind blow-up: a stronger/deeper sponge that delays or
    ! removes it confirms the terminal event is model-top dynamical.
    real(wp) :: sponge_kr = 1.0_wp/43200.0_wp   ! max Rayleigh rate [s-1]
    real(wp) :: sponge_kt = 1.0_wp/86400.0_wp   ! max Newtonian rate [s-1]
    real(wp) :: sponge_sigma = 0.12_wp          ! sponge ramp top [sigma]
    ! Adaptive hyperdiffusion knobs (defaults match aeros_timestep_class; both
    ! OFF => the fixed del^ndiff scheme, bit for bit). See the aeros_timestep
    ! header. (1) vorticity-scaled strength, (2) sigma-tapered order.
    logical  :: diff_adapt       = .FALSE.
    real(wp) :: diff_zeta_ref    = 1.0e-4_wp    ! reference RMS |zeta| [s-1]
    real(wp) :: diff_adapt_gain  = 1.0_wp       ! excess-vorticity gain [-]
    real(wp) :: diff_adapt_max   = 10.0_wp      ! cap on the multiplier [-]
    logical  :: diff_taper       = .FALSE.
    integer  :: diff_ndiff_top   = 4            ! order at the model top
    real(wp) :: diff_taper_sigma = 0.15_wp      ! taper ramp top [sigma]
    real(wp) :: seed_asym = 0.0_wp    ! zonal-asymmetry seed amplitude [K-ish]
    real(wp) :: albedo = 0.06_wp      ! surface broadband albedo (cloud proxy knob)
    real(wp) :: co2_ppm = 280.0_wp
    ! Critical RH for large-scale condensation: q relaxes to cond_rh_crit*q_sat,
    ! so the column cannot equilibrate wetter than this where condensation is the
    ! moisture sink. 1.0 = true saturation adjustment (the model default); a lower
    ! value is the sub-grid-saturation stand-in for the missing subsidence drying
    ! that otherwise leaves the free troposphere pinned at RH~100% (m2_results §23).
    real(wp) :: cond_rh_crit = 0.95_wp
    real(wp) :: cond_reevap  = 30.0_wp   ! reevaporation efficiency of falling precip
    integer  :: ocean_mode = 0       ! 0 prescribed SST, 1 slab
    real(wp) :: ocean_depth = 10.0_wp
    ! --- thermodynamic sea ice (feat/seaice) --------------------------------
    ! Default .FALSE. keeps the freeze-floor slab, bit-for-bit unchanged. When on
    ! (with a slab ocean) the mixed-layer deficit forms ice, ice grows/melts from
    ! the Semtner 0-layer surface balance, and the ice albedo feeds back into the
    ! radiation. Optional overrides like the rad/topo/restart knobs below.
    logical  :: l_seaice   = .FALSE.
    real(wp) :: ice_albedo = 0.60_wp   ! surface albedo over ice [-]
    real(wp) :: k_ice      = 2.0_wp    ! ice thermal conductivity [W m-1 K-1]
    real(wp) :: t_frz      = 271.35_wp ! seawater freezing point [K]
    real(wp) :: rad_interval = 21600.0_wp   ! radiation recompute cadence [s]; default 6 h
    integer  :: rad_scheme = SCHEME_ECCKD   ! LW/SW scheme (2=ecCKD default, 1=SESAM)
    ! --- real surface topography (orography) --------------------------------
    ! Default .FALSE. keeps the flat aquaplanet (phis = 0), bit-for-bit unchanged.
    ! When on, phis is read from topo_file (ERA5 surface geopotential, m2 s-2) and
    ! ramped in over topo_ramp_days to spread the spin-up shock. The ramp is a
    ! pure function of elapsed time (see aeros_topography_scale), hence restart-safe.
    logical  :: l_topography = .FALSE.
    character(len=512) :: topo_file = "/Users/alrobi001/data/era5/era5_orography.nc"
    real(wp) :: topo_ramp_days = 20.0_wp    ! linear 0->1 spin-up ramp [days]; 0 = full at t=0

    ! --- land surface + land-sea mask (feat/land) ---------------------------
    ! Default .FALSE. keeps the all-ocean aquaplanet, bit-for-bit unchanged. When
    ! on, the land-sea mask and land albedo are read from ERA5 climatologies and
    ! land points carry a prognostic bucket soil moisture and slab soil
    ! temperature (the land skin temperature) instead of the ocean SST.
    logical  :: l_land = .FALSE.
    character(len=512) :: lsm_file = &
        "/Users/alrobi001/data/era5/monthly-single-levels/era5_monthly-single-levels_lsm_1991-2020_clim.nc"
    character(len=512) :: land_albedo_file = &
        "/Users/alrobi001/data/era5/monthly-single-levels/era5_monthly-single-levels_fal_1991-2020_clim.nc"
    real(wp) :: w_field_capacity = 0.15_wp  ! bucket capacity [m]
    real(wp) :: w_crit           = -1.0_wp  ! beta knee [m]; <0 => 0.75*capacity
    real(wp) :: c_soil           = 2.0e6_wp ! soil areal heat capacity [J m-2 K-1]
    real(wp) :: land_albedo      = 0.20_wp  ! fallback land albedo (no albedo file)

    ! --- prognostic cloud fraction (feat/clouds) ----------------------------
    ! Default off = the diagnostic RH->cover scheme, bit-for-bit unchanged. When
    ! on, radiation consumes the Sundqvist prognostic cf (aeros_cloud_prog).
    logical  :: l_cloud_prog   = .FALSE.
    real(wp) :: cloud_rhc_sfc  = 0.40_wp    ! critical RH at the surface (calibrated)
    real(wp) :: cloud_rhc_top  = 0.60_wp    ! critical RH at the model top (higher aloft)
    real(wp) :: cloud_tau_form = 10800.0_wp ! formation timescale [s] (~3 h; de-patchifies)
    real(wp) :: cloud_tau_evap = 10800.0_wp ! evaporation timescale [s] (~3 h)
    real(wp) :: cloud_c_detr   = 0.5_wp     ! convective detrainment anvil ceiling [-]

    ! --- checkpoint / restart ------------------------------------------------
    ! restart_in  : path to a restart file to resume from ("" = cold start).
    ! restart_out : path to write checkpoints to ("" = none).
    ! restart_interval : write a checkpoint every N steps (0 = only at end, if
    !                    restart_out is set). Defaults preserve current behaviour.
    character(len=512) :: restart_in  = ""
    character(len=512) :: restart_out = ""
    integer            :: restart_interval = 0
    real(wp)           :: restart_time = 0.0_wp
    integer            :: n0 = 0          ! absolute step offset (nonzero on restart)

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)     :: grd
    type(aeros_vgrid_class)    :: vg
    type(aeros_param_class)    :: par
    type(aeros_state_class)    :: now
    type(aeros_timestep_class) :: ts

    real(wp), allocatable :: phis2(:,:)       ! surface geopotential fed to the dynamics
    real(wp), allocatable :: phis_full(:,:)   ! full (unramped) topography [m2 s-2]
    real(wp), allocatable :: w_tavg(:,:,:)    ! time-mean omega accumulator [Pa/s]
    integer            :: n_wacc = 0          ! samples in w_tavg
    ! Time-mean per-term diabatic heating accumulators [K/step], same 2nd-half
    ! sampling as w_tavg -- to see WHERE (lat,sigma) the latent heating fires.
    real(wp), allocatable :: cnv_tavg(:,:,:), cnd_tavg(:,:,:), &
                             rad_tavg(:,:,:), surf_tavg(:,:,:)
    ! Time-mean zonal-mean jet/MMC/eddy-momentum-flux accumulators (lat,lev),
    ! same 2nd-half sampling as w_tavg. uv_tavg = [u*v*] (zonal-mean of the
    ! deviations from the instantaneous zonal mean) -- the northward eddy flux
    ! of zonal momentum whose meridional convergence sets the Hadley-cell edge.
    real(wp), allocatable :: uv_tavg(:,:), ubar_tavg(:,:), vbar_tavg(:,:)
    ! Eddy zonal-wavenumber spectrum (mass+area-weighted, time-mean over the 2nd
    ! half): ke_spec(m) = eddy KE at zonal wavenumber m, mf_spec(m) = the [u*v*]
    ! momentum-flux co-spectrum. Answers whether the eddies carry SW-like KE but
    ! at the wrong scale / with no coherent flux. ctab/stab are the precomputed
    ! DFT trig table so the per-step accumulation is multiply-add only.
    integer               :: mmax_spec = 0
    real(wp), allocatable :: ke_spec(:), mf_spec(:)     ! (1:mmax_spec)
    real(wp), allocatable :: ctab(:,:), stab(:,:)       ! (mmax_spec, nlon)
    ! Where each wavenumber's eddy KE lives -- to place the m=8 spike (on the jet
    ! => physical jet instability; spread/grid-tied => numerical).
    real(wp), allocatable :: ke_latm(:,:)               ! (nlat, mmax) mass-wtd over levels
    real(wp), allocatable :: ke_levm(:,:)               ! (nlev, mmax) area-wtd over lats
    real(wp) :: tscale, tscale_prev           ! current / previous ramp factor
    real(wp) :: qs, dqsdt, tval
    real(wp) :: phalf(0:64), pfull(64), dpc(64)
    integer  :: i, j, k, n, m
    logical  :: blew_up

    nmlfile = "rce.nml"
    if (command_argument_count() >= 1) call get_command_argument(1, nmlfile)
    if (command_argument_count() >= 2) call get_command_argument(2, rh_out)
    call read_config()

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev, t_ref=t_ref)

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
    par%si_alpha   = si_alpha
    par%ndiff = ndiff; par%tau_diff = tau_diff
    par%mass_fixer = .FALSE.

    call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
    call aeros_timestep_init(ts, par, pool, grd, vg)

    ts%couple_diabatic = couple_diabatic
    ts%mst%vert_vanleer = vert_vanleer

    ts%cnd%enabled = l_cnd
    ts%cnd%rh_crit = cond_rh_crit
    ts%cnd%reevap  = cond_reevap
    ts%cnv%enabled = l_cnv
    ts%cnv%tau = conv_tau
    ts%cnv%dry_adjust = l_dry_adjust
    select case (trim(conv_scheme))
    case ("sbm");          ts%cnv%scheme = SCHEME_SBM
    case ("sbm_frierson"); ts%cnv%scheme = SCHEME_SBM_FRIERSON
    case ("manabe");       ts%cnv%scheme = SCHEME_MANABE
    case default
        write(*,*) "rce_long:: error: unknown conv_scheme '"//trim(conv_scheme)// &
                   "' (expected 'sbm', 'sbm_frierson' or 'manabe')"
        error stop 1
    end select
    ts%surf%enabled = l_surf
    ts%surf%c_h = c_h; ts%surf%c_e = c_e; ts%surf%u_min = u_min
    ts%surf%c_d = c_d
    ts%rad%enabled  = l_rad
    ts%rad%clouds   = l_rad_clouds
    ! Prognostic cloud fraction: when on, force radiation's all-sky path on too
    ! (clouds=.TRUE.) so it consumes the prognostic cf rather than clear-sky.
    ts%cpr%enabled  = l_cloud_prog
    ts%cpr%rhc_sfc  = cloud_rhc_sfc
    ts%cpr%rhc_top  = cloud_rhc_top
    ts%cpr%tau_form = cloud_tau_form
    ts%cpr%tau_evap = cloud_tau_evap
    ts%cpr%c_detr   = cloud_c_detr
    if (l_cloud_prog) ts%rad%clouds = .TRUE.
    ts%sponge_on    = l_sponge
    call aeros_timestep_set_sponge(ts, vg, sponge_kr, sponge_kt, sponge_sigma)
    ! Adaptive hyperdiffusion (numerical safety net for the top thermal-wind /
    ! jet blow-up). (1) vorticity-scaled strength: scalar fields set directly.
    ts%diff_adapt      = diff_adapt
    ts%diff_zeta_ref   = diff_zeta_ref
    ts%diff_adapt_gain = diff_adapt_gain
    ts%diff_adapt_max  = diff_adapt_max
    ! (2) sigma-tapered order: the setter rebuilds the per-level ratios.
    if (diff_taper) call aeros_timestep_set_diff_taper(ts, vg, .TRUE., &
                            ndiff_top=diff_ndiff_top, taper_sigma=diff_taper_sigma)
    ! Divergence-only diffusion timescale (SpeedyWeather-style). Namelist is in
    ! hours; 0 => keep tau_diff (bit-for-bit single-timescale scheme).
    ts%tau_diff_div = tau_diff_div*3600.0_wp
    ts%vd%enabled    = l_vdiff
    ts%vd%k0         = vdiff_k0
    ts%vd%sigma      = vdiff_sigma
    ts%vd%richardson = vdiff_richardson
    ts%vd%ri_crit    = vdiff_ri_crit
    ts%rad%albedo   = albedo
    ts%rad%co2_ppm  = co2_ppm
    ts%rad%interval = rad_interval
    ts%rad%scheme   = rad_scheme
    ts%ocn%mode     = ocean_mode
    ts%ocn%depth    = ocean_depth
    ts%ocn%l_seaice   = l_seaice
    ts%ocn%ice_albedo = ice_albedo
    ts%ocn%ocn_albedo = albedo        ! open-water albedo = the `albedo` knob (one source of truth)
    ts%ocn%k_ice      = k_ice
    ts%ocn%t_frz      = t_frz
    call aeros_ocean_init(ts%ocn, grd)   ! recompute C, (re)allocate ice state

    ! Land surface (feat/land). Configure the land state from the namelist and
    ! init it (reads the land-sea mask and albedo maps when on). Then compose the
    ! radiation surface-albedo map: the ocean scalar albedo (set just above) over
    ! sea, the land albedo on land. Off by default -> all-ocean, bit-for-bit.
    ts%land%enabled          = l_land
    ts%land%lsm_file         = lsm_file
    ts%land%albedo_file      = land_albedo_file
    ts%land%w_field_capacity = w_field_capacity
    ts%land%w_crit           = w_crit
    ts%land%c_soil           = c_soil
    ts%land%land_albedo      = land_albedo
    call aeros_land_init(ts%land, grd, l_land)
    if (l_land) then
        call aeros_land_couple_radiation(ts%land, ts%rad%alb_map, ts%rad%albedo)
        write(*,"(a,a)") " rce_long:: land-sea mask from ", trim(lsm_file)
        write(*,"(a,i0,a,i0,a)") " rce_long:: land cells ", count(ts%land%mask), &
            " of ", grd%nlon*grd%nlat, " (mask threshold 0.5)"
    end if

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

    ! Surface geopotential (lower boundary). Aquaplanet phis = 0 by default; with
    ! l_topography the real orography is read once and ramped in over the run.
    allocate(phis2(grd%nlon, grd%nlat));     phis2 = 0.0_wp
    allocate(phis_full(grd%nlon, grd%nlat)); phis_full = 0.0_wp
    if (l_topography) then
        call aeros_topography_load(topo_file, grd%lon, grd%lat, phis_full)
        write(*,"(a,a)") " rce_long:: topography from ", trim(topo_file)
        write(*,"(a,es11.3,a,es11.3)") " rce_long:: raw   phis min ", &
            minval(phis_full), "  max ", maxval(phis_full)
        ! Spectrally truncate to the model resolution (T21). A spectral core can
        ! only carry the surface geopotential at its resolved scales: the raw
        ! interpolated field still has grid-scale structure at coastlines and
        ! mountain flanks, which aliases in the transforms and rings (Gibbs),
        ! destabilizing the lowest layer. Band-limiting via one analysis ->
        ! synthesis round-trip (SHTns truncates at T21) gives the smooth
        ! orography the dynamics is consistent with. Standard spectral-model
        ! practice.
        block
            complex(wp_sh), allocatable :: phis_lm(:)
            real(dp),       allocatable :: gwork(:,:)
            allocate(phis_lm(pool%sht(1)%nlm))
            allocate(gwork(grd%nlon, grd%nlat))
            gwork = real(phis_full, dp)
            call aeros_sht_analysis(pool%sht(1), gwork, phis_lm)   ! grid -> T21 spectral (overwrites gwork)
            call aeros_sht_synthesis(pool%sht(1), phis_lm, gwork)  ! T21 spectral -> band-limited grid
            phis_full = real(gwork, wp)
            deallocate(phis_lm, gwork)
        end block
        write(*,"(a,es11.3,a,es11.3,a,f6.1,a)") " rce_long:: T21   phis min ", &
            minval(phis_full), "  max ", maxval(phis_full), " m2/s2   ramp ", &
            topo_ramp_days, " days"
        ! Dump the band-limited orography on the model grid for a visual sanity
        ! check (correct continents/mountains, sane magnitude).
        call nc_create("output/rce_phis.nc")
        call nc_write_dim("output/rce_phis.nc", "lon", x=grd%lon, units="degrees_east")
        call nc_write_dim("output/rce_phis.nc", "lat", x=grd%lat, units="degrees_north")
        call nc_write("output/rce_phis.nc", "phis", phis_full, dim1="lon", dim2="lat", &
                      units="m2 s-2", long_name="surface geopotential (T21)")
    end if
    ! Initial feed at t = 0 (scale = 0 while ramping, so the run starts flat).
    tscale      = aeros_topography_scale(l_topography, 0.0_wp, topo_ramp_days)
    tscale_prev = tscale
    phis2 = tscale*phis_full
    call aeros_timestep_set_phis(ts, phis2)

    ! Diagnostic prescribed heating: an analytic Q(lat,sigma), equatorial and
    ! area-weighted zero-mean in latitude (no net global heating -> no drift),
    ! sin(pi*sigma) in the vertical (mid-troposphere peak, zero at top/surface).
    ! Coupled on the in-solve seam to isolate the dry core's heating->ω response
    ! from moist physics. Run with all physics off.
    if (l_qforce) then
        block
            real(wp), allocatable :: q3(:,:,:), slat(:)
            real(wp) :: amp, mlat, wsum, latr, pi
            integer  :: iq, jq, kq
            pi = 4.0_wp*atan(1.0_wp)
            allocate(q3(grd%nlon, grd%nlat, nlev), slat(grd%nlat))
            amp = qforce_amp/86400.0_wp                     ! K/day -> K/s
            wsum = 0.0_wp; mlat = 0.0_wp
            do jq = 1, grd%nlat
                latr    = grd%lat(jq)
                slat(jq) = exp(-(latr/15.0_wp)**2)          ! equatorial heating
                mlat    = mlat + slat(jq)*cos(latr*pi/180.0_wp)
                wsum    = wsum + cos(latr*pi/180.0_wp)
            end do
            mlat = mlat/wsum                                ! area-weighted mean
            do kq = 1, nlev
                do jq = 1, grd%nlat
                    do iq = 1, grd%nlon
                        q3(iq,jq,kq) = amp*(slat(jq) - mlat) &
                                        *sin(pi*vg%sigma_full(kq))
                    end do
                end do
            end do
            call aeros_timestep_set_qforce(ts, q3)
            write(*,'(a,f6.2,a)') " rce_long:: prescribed Q ON, amp=", qforce_amp, &
                                    " K/day (physics MUST be off for a clean dry test)"
            deallocate(q3, slat)
        end block
    end if

    ! === Restart branch: load state instead of constructing an IC ===========
    ! When restart_in is set, the timestep and state are already allocated by
    ! the init/config above (same trunc/nlev/grid); read_restart overwrites them
    ! with the saved leapfrog levels, humidity, ocean and radiation cache, and
    ! the run resumes from the saved nstep. The whole IC-construction block is
    ! skipped in that case.
    if (len_trim(restart_in) > 0) then
        call aeros_timestep_read_restart(ts, now, restart_time, restart_in)
        n0 = ts%nstep
        write(*,"(a,a,a,i0,a,f10.1,a)") " rce_long:: restarted from ", &
            trim(restart_in), "  (nstep=", n0, ", time=", restart_time, " s)"
        call report(n0)
    else

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

    end if   ! restart_in vs cold-start IC

    blew_up = .FALSE.
    allocate(w_tavg(grd%nlon,grd%nlat,nlev)); w_tavg = 0.0_wp; n_wacc = 0
    allocate(cnv_tavg(grd%nlon,grd%nlat,nlev), cnd_tavg(grd%nlon,grd%nlat,nlev), &
             rad_tavg(grd%nlon,grd%nlat,nlev), surf_tavg(grd%nlon,grd%nlat,nlev))
    cnv_tavg = 0.0_wp; cnd_tavg = 0.0_wp; rad_tavg = 0.0_wp; surf_tavg = 0.0_wp
    allocate(uv_tavg(grd%nlat,nlev), ubar_tavg(grd%nlat,nlev), vbar_tavg(grd%nlat,nlev))
    uv_tavg = 0.0_wp; ubar_tavg = 0.0_wp; vbar_tavg = 0.0_wp
    ! Eddy wavenumber spectrum: resolve m = 1..trunc; build the DFT trig table.
    mmax_spec = trunc
    allocate(ke_spec(mmax_spec), mf_spec(mmax_spec))
    allocate(ctab(mmax_spec, grd%nlon), stab(mmax_spec, grd%nlon))
    allocate(ke_latm(grd%nlat, mmax_spec), ke_levm(nlev, mmax_spec))
    ke_spec = 0.0_wp; mf_spec = 0.0_wp; ke_latm = 0.0_wp; ke_levm = 0.0_wp
    block
        integer :: mm, ii
        real(wp) :: ang, twopi
        twopi = 8.0_wp*atan(1.0_wp)
        do ii = 1, grd%nlon
            do mm = 1, mmax_spec
                ang = twopi*real(mm*(ii-1), wp)/real(grd%nlon, wp)
                ctab(mm,ii) = cos(ang); stab(mm,ii) = -sin(ang)
            end do
        end do
    end block
    do n = n0+1, n0+nstep
        ! Advance the topography ramp: a pure function of absolute elapsed time
        ! n*dt, so it is restart-safe -- a run resumed at n0 continues the ramp at
        ! the right amplitude. Only re-set phis while the factor is still changing
        ! (during the ramp); once it reaches 1 the surface is constant.
        if (l_topography .and. topo_ramp_days > 0.0_wp .and. tscale_prev < 1.0_wp) then
            tscale = aeros_topography_scale(.true., real(n,wp)*dt, topo_ramp_days)
            if (tscale /= tscale_prev) then
                phis2 = tscale*phis_full
                call aeros_timestep_set_phis(ts, phis2)
                tscale_prev = tscale
            end if
        end if
        call aeros_timestep_step(ts, pool, vg, grd, now)
        if (any(ts%wrk%t_g /= ts%wrk%t_g) .or. any(now%qv_g /= now%qv_g)) then
            write(*,"(a,i0,a,f7.2,a)") " *** NaN at step ", n, "  (day ", &
                real(n,wp)*dt/86400.0_wp, ")"
            call locate_nan()
            blew_up = .TRUE.
            exit
        end if
        if (mod(n, print_every) == 0) call report(n)
        ! Time-mean omega over the run's second half (equilibrium) for a robust
        ! zonal-mean subsidence diagnostic -- a single snapshot is too noisy.
        if (l_diag .and. n > n0 + nstep/2) then
            w_tavg = w_tavg + ts%wrk%omega
            cnv_tavg  = cnv_tavg  + ts%wrk%dt_cnv
            cnd_tavg  = cnd_tavg  + ts%wrk%dt_cnd
            rad_tavg  = rad_tavg  + ts%wrk%dt_rad
            surf_tavg = surf_tavg + ts%wrk%dt_surf
            call accum_zmflux()
            call accum_espec()
            n_wacc = n_wacc + 1
        end if
        ! Periodic checkpoint (restart_interval > 0). The model time carried into
        ! the file is the absolute elapsed time n*dt.
        if (len_trim(restart_out) > 0 .and. restart_interval > 0) then
            if (mod(n, restart_interval) == 0) &
                call aeros_timestep_write_restart(ts, now, real(n,wp)*dt, restart_out)
        end if
    end do

    if (.not. blew_up) write(*,"(a)") " rce_long:: completed without NaN"

    ! End-of-run checkpoint: always written when restart_out is set (this is the
    ! only checkpoint when restart_interval = 0). Skipped on a blow-up, which
    ! would only persist a NaN state.
    if (len_trim(restart_out) > 0 .and. .not. blew_up) then
        call aeros_timestep_write_restart(ts, now, real(n0+nstep,wp)*dt, restart_out)
        write(*,"(a,a)") " rce_long:: wrote restart -> ", trim(restart_out)
    end if

    if (len_trim(rh_out) > 0) call dump_rh(trim(rh_out))

    call aeros_timestep_end(ts)
    call aeros_state_end(now)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

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
        real(wp) :: phc(0:nlev), pfc(nlev), dpcc(nlev)
        real(wp) :: cf(nlev), clwc(nlev), ciwc(nlev)
        real(wp) :: qsl, dqsl, colcov, rnlon
        real(wp), allocatable :: rh_zm(:,:), cf_zm(:,:), t_zm(:,:), q_zm(:,:), p_zm(:,:)
        real(wp), allocatable :: u_zm(:,:), v_zm(:,:), w_zm(:,:)
        real(wp), allocatable :: cnv_zm(:,:), cnd_zm(:,:), rad_zm(:,:), &
                                 surf_zm(:,:), qnet_zm(:,:)
        real(wp), allocatable :: cover(:), latout(:), levout(:)
        real(wp) :: kday
        integer :: ii, jj, kk

        rnlon = real(grd%nlon, wp)
        allocate(rh_zm(grd%nlat,nlev), cf_zm(grd%nlat,nlev), t_zm(grd%nlat,nlev), &
                 q_zm(grd%nlat,nlev), p_zm(grd%nlat,nlev), &
                 u_zm(grd%nlat,nlev), v_zm(grd%nlat,nlev), w_zm(grd%nlat,nlev), &
                 cover(grd%nlat), latout(grd%nlat), levout(nlev))
        rh_zm = 0.0_wp; cf_zm = 0.0_wp; t_zm = 0.0_wp; q_zm = 0.0_wp
        p_zm = 0.0_wp; u_zm = 0.0_wp; v_zm = 0.0_wp; w_zm = 0.0_wp; cover = 0.0_wp
        allocate(cnv_zm(grd%nlat,nlev), cnd_zm(grd%nlat,nlev), rad_zm(grd%nlat,nlev), &
                 surf_zm(grd%nlat,nlev), qnet_zm(grd%nlat,nlev))
        cnv_zm = 0.0_wp; cnd_zm = 0.0_wp; rad_zm = 0.0_wp; surf_zm = 0.0_wp
        latout = grd%lat(1:grd%nlat)
        levout = vg%sigma_full(1:nlev)
        ! [K/step] time-mean -> [K/day]: /n_wacc (time), /nlon (zonal), *86400/dt.
        kday = 86400.0_wp/dt

        do jj = 1, grd%nlat
            do ii = 1, grd%nlon
                call aeros_vgrid_pressure(vg, now%ps(ii,jj), phc, pfc, dpcc)
                call aeros_cloud_diagnose(nlev, now%temp_g(ii,jj,:), &
                    now%qv_g(ii,jj,:), pfc, now%ps(ii,jj), cf, clwc, ciwc)
                colcov = 0.0_wp
                do kk = 1, nlev
                    call aeros_qsat(now%temp_g(ii,jj,kk), pfc(kk), qsl, dqsl)
                    rh_zm(jj,kk) = rh_zm(jj,kk) + &
                        now%qv_g(ii,jj,kk)/max(qsl,1.0e-12_wp)*100.0_wp
                    cf_zm(jj,kk) = cf_zm(jj,kk) + cf(kk)
                    t_zm(jj,kk)  = t_zm(jj,kk)  + now%temp_g(ii,jj,kk)
                    q_zm(jj,kk)  = q_zm(jj,kk)  + now%qv_g(ii,jj,kk)*1000.0_wp
                    p_zm(jj,kk)  = p_zm(jj,kk)  + pfc(kk)/100.0_wp
                    u_zm(jj,kk)  = u_zm(jj,kk)  + now%u(ii,jj,kk)
                    v_zm(jj,kk)  = v_zm(jj,kk)  + now%v(ii,jj,kk)
                    ! time-mean vertical pressure velocity [hPa/day], >0 = subsidence
                    if (n_wacc > 0) w_zm(jj,kk) = w_zm(jj,kk) &
                        + w_tavg(ii,jj,kk)/real(n_wacc,wp)*864.0_wp
                    ! per-term diabatic heating, time-mean [K/day]
                    if (n_wacc > 0) then
                        cnv_zm(jj,kk)  = cnv_zm(jj,kk)  + cnv_tavg(ii,jj,kk) /real(n_wacc,wp)*kday
                        cnd_zm(jj,kk)  = cnd_zm(jj,kk)  + cnd_tavg(ii,jj,kk) /real(n_wacc,wp)*kday
                        rad_zm(jj,kk)  = rad_zm(jj,kk)  + rad_tavg(ii,jj,kk) /real(n_wacc,wp)*kday
                        surf_zm(jj,kk) = surf_zm(jj,kk) + surf_tavg(ii,jj,kk)/real(n_wacc,wp)*kday
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
            qnet_zm(jj,:) = cnv_zm(jj,:) + cnd_zm(jj,:) + rad_zm(jj,:) + surf_zm(jj,:)
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
        if (n_wacc > 0) then
            call nc_write(fname, "ubar", ubar_tavg/real(n_wacc,wp), dim1="lat", dim2="lev", &
                units="m/s", long_name="time-mean zonal-mean zonal wind (2nd half)")
            call nc_write(fname, "vbar", vbar_tavg/real(n_wacc,wp), dim1="lat", dim2="lev", &
                units="m/s", long_name="time-mean zonal-mean meridional wind (2nd half)")
            call nc_write(fname, "uvpr", uv_tavg/real(n_wacc,wp), dim1="lat", dim2="lev", &
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
        call nc_write(fname, "q_net",  qnet_zm, dim1="lat", dim2="lev", units="K/day", &
            long_name="zonal-mean net diabatic heating")
        write(*,"(a)") " rce_long:: wrote zonal-mean RH/cf dump -> "//trim(fname)

        if (n_wacc > 0 .and. allocated(ke_spec)) then
            block
                integer  :: m
                real(wp) :: ket, mft
                ket = sum(ke_spec)/real(n_wacc,wp); mft = sum(mf_spec)/real(n_wacc,wp)
                write(*,"(a)") " eddy zonal-wavenumber spectrum (mass+area wtd, 2nd-half mean):"
                write(*,"(a,es10.3,a,es10.3)") "   total eddy KE ", ket, "   total [u*v*] ", mft
                write(*,"(a)") "     m    KE(m)      KE%     [u*v*](m)   flux%"
                do m = 1, mmax_spec
                    write(*,"(i6,es11.3,f8.1,es12.3,f8.1)") m, ke_spec(m)/real(n_wacc,wp), &
                        100.0_wp*ke_spec(m)/max(sum(ke_spec),tiny(1.0_wp)), &
                        mf_spec(m)/real(n_wacc,wp), &
                        100.0_wp*mf_spec(m)/max(abs(sum(mf_spec)),tiny(1.0_wp))
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
                    if (mm > mmax_spec) cycle
                    write(*,"(a,i0,a)") " m=", mm, " eddy KE by latitude (%, N+S folded):"
                    tot = max(sum(ke_latm(:,mm)), tiny(1.0_wp))
                    do j = 1, grd%nlat
                        if (100.0_wp*ke_latm(j,mm)/tot > 3.0_wp) &
                            write(*,"(a,f6.1,a,f6.1)") "    lat ", grd%lat(j), " : ", &
                                100.0_wp*ke_latm(j,mm)/tot
                    end do
                    write(*,"(a,i0,a)") " m=", mm, " eddy KE by level (%, top=1):"
                    tot = max(sum(ke_levm(:,mm)), tiny(1.0_wp))
                    do k = 1, nlev
                        write(*,"(a,i3,a,f6.1)") "    lev ", k, " : ", 100.0_wp*ke_levm(k,mm)/tot
                    end do
                end do
            end block
        end if

        deallocate(rh_zm, cf_zm, t_zm, q_zm, p_zm, u_zm, v_zm, cover, latout, levout)
        return
    end subroutine dump_rh

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
        call nml_read(nmlfile, "rce", "t_ref", t_ref, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "couple_diabatic", couple_diabatic, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "vert_vanleer", vert_vanleer, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "tau_diff", tau_diff)
        call nml_read(nmlfile, "rce", "ndiff", ndiff)
        call nml_read(nmlfile, "rce", "eps_filter", eps_filter)
        call nml_read(nmlfile, "rce", "raw_alpha", raw_alpha)
        call nml_read(nmlfile, "rce", "conv_tau", conv_tau)
        call nml_read(nmlfile, "rce", "c_h", c_h)
        call nml_read(nmlfile, "rce", "c_e", c_e)
        call nml_read(nmlfile, "rce", "c_d", c_d)
        call nml_read(nmlfile, "rce", "u_min", u_min)
        call nml_read(nmlfile, "rce", "print_every", print_every)
        call nml_read(nmlfile, "rce", "l_surf", l_surf)
        call nml_read(nmlfile, "rce", "l_cnv", l_cnv)
        call nml_read(nmlfile, "rce", "l_dry_adjust", l_dry_adjust)
        call nml_read(nmlfile, "rce", "l_uniform_insol", l_uniform_insol)
        call nml_read(nmlfile, "rce", "l_nonrotating", l_nonrotating)
        call nml_read(nmlfile, "rce", "l_qforce", l_qforce, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "qforce_amp", qforce_amp, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "sponge_kr", sponge_kr)
        call nml_read(nmlfile, "rce", "sponge_kt", sponge_kt)
        call nml_read(nmlfile, "rce", "sponge_sigma", sponge_sigma)
        call nml_read(nmlfile, "rce", "l_cnd", l_cnd)
        call nml_read(nmlfile, "rce", "l_rad", l_rad)
        call nml_read(nmlfile, "rce", "l_rad_clouds", l_rad_clouds)
        call nml_read(nmlfile, "rce", "l_sponge", l_sponge)
        call nml_read(nmlfile, "rce", "l_vdiff", l_vdiff)
        call nml_read(nmlfile, "rce", "l_diag", l_diag)
        call nml_read(nmlfile, "rce", "vdiff_k0", vdiff_k0)
        call nml_read(nmlfile, "rce", "vdiff_sigma", vdiff_sigma)
        ! Boundary-layer scheme (optional; inherits input/rce_defaults.nml).
        call nml_read(nmlfile, "rce", "vdiff_richardson", vdiff_richardson, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "vdiff_ri_crit", vdiff_ri_crit, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "seed_asym", seed_asym)
        call nml_read(nmlfile, "rce", "albedo", albedo)
        call nml_read(nmlfile, "rce", "co2_ppm", co2_ppm)
        call nml_read(nmlfile, "rce", "cond_rh_crit", cond_rh_crit)
        call nml_read(nmlfile, "rce", "cond_reevap", cond_reevap, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "ocean_mode", ocean_mode)
        call nml_read(nmlfile, "rce", "ocean_depth", ocean_depth)
        ! Sea-ice knobs (feat/seaice): optional (inherit input/rce_defaults.nml)
        ! so existing rce_*.nml files that omit them keep the freeze-floor slab.
        call nml_read(nmlfile, "rce", "l_seaice", l_seaice, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "ice_albedo", ice_albedo, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "k_ice", k_ice, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "t_frz", t_frz, &
                      defaults_file="input/rce_defaults.nml")
        ! Optional overrides: a namelist that omits these inherits input/rce_defaults.nml
        ! (6 h / ecCKD) instead of erroring. Every other key above stays required.
        call nml_read(nmlfile, "rce", "rad_interval", rad_interval, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "rad_scheme", rad_scheme, &
                      defaults_file="input/rce_defaults.nml")
        ! Topography (feat/topography). Optional overrides like the rad knobs
        ! above: a namelist that omits them inherits input/rce_defaults.nml
        ! (l_topography off), so existing run namelists are unaffected.
        ! Convection scheme selector. Optional (inherits input/rce_defaults.nml =
        ! "sbm"), so existing namelists that omit it keep the default scheme.
        call nml_read(nmlfile, "rce", "conv_scheme", conv_scheme, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "l_topography", l_topography, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "topo_file", topo_file, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "topo_ramp_days", topo_ramp_days, &
                      defaults_file="input/rce_defaults.nml")
        ! Land (feat/land). Optional overrides like the topo knobs above: a
        ! namelist that omits them inherits input/rce_defaults.nml (l_land off),
        ! so existing run namelists are unaffected.
        call nml_read(nmlfile, "rce", "l_land", l_land, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "lsm_file", lsm_file, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "land_albedo_file", land_albedo_file, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "w_field_capacity", w_field_capacity, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "w_crit", w_crit, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "c_soil", c_soil, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "land_albedo", land_albedo, &
                      defaults_file="input/rce_defaults.nml")
        ! Restart knobs (feat/restart): optional (inherit input/rce_defaults.nml)
        ! so existing rce_*.nml files that omit them keep cold-starting, not erroring.
        call nml_read(nmlfile, "rce", "restart_in", restart_in, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "restart_out", restart_out, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "restart_interval", restart_interval, &
                      defaults_file="input/rce_defaults.nml")
        ! Prognostic cloud knobs (feat/clouds): optional (inherit rce_defaults.nml)
        ! so existing rce_*.nml files that omit them keep the diagnostic scheme.
        call nml_read(nmlfile, "rce", "l_cloud_prog", l_cloud_prog, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "cloud_rhc_sfc", cloud_rhc_sfc, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "cloud_rhc_top", cloud_rhc_top, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "cloud_tau_form", cloud_tau_form, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "cloud_tau_evap", cloud_tau_evap, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "cloud_c_detr", cloud_c_detr, &
                      defaults_file="input/rce_defaults.nml")
        ! Adaptive hyperdiffusion knobs (feat/adaptive-diff): optional (inherit
        ! rce_defaults.nml, both OFF) so existing rce_*.nml files that omit them
        ! keep the fixed del^ndiff scheme, bit for bit.
        call nml_read(nmlfile, "rce", "diff_adapt", diff_adapt, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "diff_zeta_ref", diff_zeta_ref, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "diff_adapt_gain", diff_adapt_gain, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "diff_adapt_max", diff_adapt_max, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "diff_taper", diff_taper, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "diff_ndiff_top", diff_ndiff_top, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "diff_taper_sigma", diff_taper_sigma, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "tau_diff_div", tau_diff_div, &
                      defaults_file="input/rce_defaults.nml")
        call nml_read(nmlfile, "rce", "si_alpha", si_alpha, &
                      defaults_file="input/rce_defaults.nml")
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
        if (l_diag) then
            call term_table(hot_lat_index())          ! tropical (warmest low level)
            call term_table(lat_index_near(45.0_wp))  ! midlatitude
        end if
        call energy_balance()
        if (l_cloud_prog) call cloud_prog_cover()
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
        do jj = 1, grd%nlat
            do ii = 1, grd%nlon
                colcov = 0.0_wp
                do kk = 1, nlev
                    colcov = max(colcov, now%cf_g(ii,jj,kk))
                end do
                w = grd%area(ii,jj)
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
        do j = 1, grd%nlat
            zm = sum(now%temp_g(:,j,nlev))/real(grd%nlon, wp)
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
        do j = 1, grd%nlat
            d = abs(grd%lat(j) - target)
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
        zmax = sum(now%temp_g(:,jstar,nlev))/real(grd%nlon, wp)
        write(*,"(a,i0,a,f6.1,a,f7.2,a)") &
            "   per-term heating [K/day] at j=", jstar, &
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

    subroutine accum_zmflux()
        ! Accumulate, over the equilibrated 2nd half (same window as w_tavg), the
        ! time-mean zonal-mean jet [u], meridional wind [v], and eddy momentum
        ! flux [u*v*] = zonal mean of (u-[u])(v-[v]) at every (lat,lev). The
        ! meridional convergence of [u*v*] is the eddy forcing that terminates
        ! the Hadley cell -- the cell-edge diagnostic this branch is after.
        integer  :: i, j, k
        real(wp) :: uzm, vzm, up, vp, cov, rn
        rn = real(grd%nlon, wp)
        do k = 1, nlev
            do j = 1, grd%nlat
                uzm = sum(now%u(:,j,k))/rn
                vzm = sum(now%v(:,j,k))/rn
                cov = 0.0_wp
                do i = 1, grd%nlon
                    up = now%u(i,j,k) - uzm
                    vp = now%v(i,j,k) - vzm
                    cov = cov + up*vp
                end do
                ubar_tavg(j,k) = ubar_tavg(j,k) + uzm
                vbar_tavg(j,k) = vbar_tavg(j,k) + vzm
                uv_tavg(j,k)   = uv_tavg(j,k)   + cov/rn
            end do
        end do
        return
    end subroutine accum_zmflux

    subroutine accum_espec()
        ! Eddy KE and [u*v*] momentum-flux co-spectrum by zonal wavenumber m,
        ! mass (dsigma) + area weighted and summed over the globe, accumulated
        ! over the 2nd half. One-sided (factor 2 for +-m): sum_m ke_spec(m) is the
        ! total eddy KE, sum_m mf_spec(m) the total momentum flux.
        integer  :: i, j, k, m
        real(wp) :: uzm, vzm, up, vp, w, dsig, rn, ur, ui, vr, vi, c, s, kem
        rn = real(grd%nlon, wp)
        do k = 1, nlev
            dsig = vg%sigma_half(k) - vg%sigma_half(k-1)
            do j = 1, grd%nlat
                w   = grd%area(1,j)*dsig
                uzm = sum(now%u(:,j,k))/rn
                vzm = sum(now%v(:,j,k))/rn
                do m = 1, mmax_spec
                    ur = 0.0_wp; ui = 0.0_wp; vr = 0.0_wp; vi = 0.0_wp
                    do i = 1, grd%nlon
                        up = now%u(i,j,k) - uzm
                        vp = now%v(i,j,k) - vzm
                        c  = ctab(m,i); s = stab(m,i)
                        ur = ur + up*c; ui = ui + up*s
                        vr = vr + vp*c; vi = vi + vp*s
                    end do
                    ur = ur/rn; ui = ui/rn; vr = vr/rn; vi = vi/rn
                    kem = ur*ur + ui*ui + vr*vr + vi*vi
                    ke_spec(m)   = ke_spec(m)   + w*kem
                    ke_latm(j,m) = ke_latm(j,m) + dsig*kem
                    ke_levm(k,m) = ke_levm(k,m) + grd%area(1,j)*kem
                    mf_spec(m) = mf_spec(m) + w*2.0_wp*(ur*vr + ui*vi)
                end do
            end do
        end do
        return
    end subroutine accum_espec

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
