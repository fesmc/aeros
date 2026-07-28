module aeros_radiation
    ! Radiative transfer at the grid seam. Longwave first (this file's initial
    ! state), shortwave and the grid-apply/cadence coupling to follow.
    !
    ! === Provenance: the SESAM band kernel, on a resolved column ============
    !
    ! The longwave scheme is a port of CLIMBER-X/SESAM's broadband
    ! transmissivity-emissivity longwave (climber-x src/atm/lwr.f90, Ganopolski
    ! & Willeit; PIK report 81). What is taken is the *band physics*:
    !
    !   - the blackbody source B = emis sigma T^4 at each layer,
    !   - the absorber paths for water vapour, CO2 and ozone,
    !   - the empirical broadband transmission fits d_vap, d_co2, d_o3 as
    !     functions of the cumulative absorber path (with the 1.66 diffusivity),
    !   - the emissivity-method flux integral.
    !
    ! What is NOT taken is SESAM's *driver*. SESAM is a statistical-dynamical
    ! model: its lw_radiation reconstructs an analytic column (a four-segment
    ! grid built from lapse-rate descriptors gams/gamb/gamt, tam, htrop, ram,
    ! hrm via t_prof/rh_prof) and only ever uses boundary fluxes (surface,
    ! tropopause, TOA) to close a column energy budget. aeros already carries a
    ! resolved T/q/p column on the sigma grid and needs a resolved heating
    ! *profile*, so the reconstruction is discarded and the kernel is driven
    ! with the actual layer temperature, humidity, ozone and pressure.
    !
    ! On a resolved column this is, structurally, design.md section 5's option 2
    ! (CCM3-style broadband absorptivity-emissivity) carrying SESAM's tuned,
    ! validated-to-~2000-ppm coefficients. ecCKD (option 1) is the eventual
    ! endpoint and slots in behind the same `scheme` selector.
    !
    ! === Two deviations from SESAM, both deliberate ========================
    !
    !   1. Fluxes are staggered to the half levels (interfaces), with the
    !      blackbody source and the absorber amounts defined per full layer.
    !      SESAM co-locates B and the flux at its radiation levels because it
    !      only wants boundary fluxes; a GCM wants a layer heating rate, which
    !      is a flux divergence across a layer, so the fluxes belong on the
    !      interfaces that bound the layers. The transmission physics is
    !      unchanged: transmission between two interfaces is the d_* fit of the
    !      total absorber path of the layers between them (a band model is
    !      nonlinear in path, so paths accumulate and the fit is evaluated once
    !      on the total -- exactly as SESAM does, never a product of per-layer
    !      transmissions).
    !
    !   2. The water-vapour path is the exact hydrostatic layer path
    !      0.1 * q * dp/g [g cm-2] rather than SESAM's exp-profile integral
    !      am_wv = rhos q exp(...)/kappa (...). On a resolved grid the layer is
    !      thin and the analytic within-layer exponential (whose scale height
    !      comes from log(q_k/q_{k+1}), fragile when q is non-monotone) buys
    !      nothing; the hydrostatic path is the same physical column mass and is
    !      robust. CO2 and O3 keep SESAM's pressure-broadening z-weighting
    !      verbatim, evaluated at the actual layer heights.
    !
    ! === Scope of the current state ========================================
    !
    ! Clear-sky longwave and shortwave, plus an all-sky (cloudy) longwave that
    ! folds a resolved per-layer grey cloud into the clear-sky band kernel
    ! (aeros_lw_cloudy_column). SESAM's own cloudy branch (lwr_clouds) blended a
    ! single analytic slab; on the resolved sigma column aeros instead carries
    ! the cloud as a per-layer emissivity built from the condensate paths, the
    ! same "resolve it, don't reconstruct it" deviation the clear-sky port made.
    ! The all-sky path validates against ERA5's all-sky fluxes (ttr, str) and the
    ! §17 cloud radiative effect. Still to come: the all-sky shortwave sibling,
    ! and the coupled grid apply -- the column kernels are exposed and exercised
    ! offline by test_radiation and drivers/validate_era5 for now. `enabled`
    ! defaults .FALSE. as for the other physics.
    !
    ! The column kernels take an arbitrary nlev column and are agnostic to the
    ! grid they run on (that is why validate_era5 can drive them on ERA5's 37
    ! levels, not the model's). Radiation on a refined vertical grid remapped
    ! back to the transport grid is therefore a caller-side wrapper, not a kernel
    ! change: interpolate T/q/cloud up, run the kernel, remap the heating down.

    use aeros_defs,     only : dp, wp, io_unit_err, R_d, cp_d, grav, T0, p0, &
                               sigma_sb, S0, pi, aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure, &
                               aeros_hydrostatic
    use aeros_cloud,    only : aeros_cloud_diagnose
    use aeros_insolation, only : aeros_insol_class, aeros_insol_init, &
                               aeros_insol_end, aeros_insol_annual, aeros_insol_day
    use aeros_ecckd,    only : aeros_ecckd_init, &
                               aeros_ecckd_lw_clearsky_column, &
                               aeros_ecckd_sw_clearsky_column, &
                               aeros_ecckd_lw_cloudy_column, &
                               aeros_ecckd_sw_cloudy_column
    use nml,            only : nml_read

    implicit none

    private

    ! === Scheme selector ===================================================
    integer, parameter, public :: SCHEME_SESAM = 1   ! ported SESAM band kernel
    integer, parameter, public :: SCHEME_ECCKD = 2   ! ecCKD table (future)

    ! === SESAM longwave band coefficients ==================================
    ! Values are CLIMBER-X's tuned defaults (nml/atm_par.nml) and the hardcoded
    ! constants in lwr.f90. Kept as named parameters, not namelist knobs: they
    ! are one internally-consistent fit and are not meant to be tuned piecemeal.

    real(wp), parameter :: LW_EMIS   = 1.0_wp      ! atmospheric emissivity
    real(wp), parameter :: LW_BETA0  = 1.66_wp     ! diffusivity factor

    ! water vapour: d_vap = 1/(1 + vap_opac (a_vap x^b + a2 x^b2 + a3 x^3)),
    !   x = beta0*u
    ! LW_VAP_OPAC is a single documented correction to the vapour opacity, NOT a
    ! re-fit of the SESAM coefficients: it multiplies the whole vapour optical
    ! depth, so 1.0 is SESAM verbatim and < 1 makes the band less opaque. Against
    ! ERA5 clear-sky the SESAM fit runs ~16 W/m2 too opaque in the warm moist
    ! tropics (OLR too low, §14); because the term is negligible at small path
    ! and dominant at large path, one scale corrects the tropics without touching
    ! the dry columns. Tuned to null the global-mean clear-sky OLR bias; retire
    ! it when the correlated-k (ecCKD) scheme lands. See m2_results §20.
    real(wp), parameter :: LW_VAP_OPAC = 0.80_wp
    real(wp), parameter :: LW_A_VAP    = 1.6_wp
    real(wp), parameter :: LW_BETA_VAP = 0.45_wp
    real(wp), parameter :: LW_A2_VAP   = 0.1_wp
    real(wp), parameter :: LW_BETA2_VAP= 1.5_wp
    real(wp), parameter :: LW_A3_VAP   = 0.01_wp
    real(wp), parameter :: LW_AK_WV    = 1.0_wp    ! (kept for provenance; see note)

    ! CO2: d_co2 = (1 - min(0.2,0.1(u/1000)^2)) (1 + a0 a1 x^b)/(1 + a0 x^b)
    real(wp), parameter :: LW_AK_CO2   = 0.8_wp    ! pressure-broadening exponent
    real(wp), parameter :: LW_BETA_CO2 = 0.45_wp
    real(wp), parameter :: LW_A0_CO2   = 0.247_wp
    real(wp), parameter :: LW_A1_CO2   = 0.755_wp

    ! ozone: d_o3 = 1 - a_o3 u^b
    real(wp), parameter :: LW_AK_O3    = 0.6_wp    ! pressure-broadening exponent
    real(wp), parameter :: LW_A_O3     = 8.246_wp
    real(wp), parameter :: LW_BETA_O3  = 0.539_wp

    ! === Grey-cloud longwave optics ========================================
    ! Per-layer cloud emissivity eps = 1 - exp(-(k_liq LWP + k_ice IWP)), with
    ! the condensate water paths in [g m-2]. The coefficients are standard
    ! broadband LW flux (diffusivity-included) mass absorption values -- Stephens
    ! (1978)-type for liquid, smaller for the larger ice crystals -- not tuned to
    ! SESAM's slab. Named parameters, tunable in one place.
    real(wp), parameter :: LW_KABS_LIQ = 0.10_wp   ! LW mass absorption, liquid [m2 g-1]
    real(wp), parameter :: LW_KABS_ICE = 0.06_wp   ! LW mass absorption, ice    [m2 g-1]

    ! Below this cloud fraction a cloudy cell is treated as clear (guards the
    ! grid-mean -> in-cloud division by cf). Shared by the LW and SW cloud paths.
    real(wp), parameter :: CLD_CF_FLOOR = 1.0e-3_wp

    ! molar masses for the CO2 volume->mass mixing ratio conversion
    real(wp), parameter :: M_CO2 = 44.0095_wp      ! g/mol
    real(wp), parameter :: M_AIR = 28.97_wp        ! g/mol

    ! === SESAM shortwave band coefficients =================================
    ! Two bands: visible+UV (fraction SW_FRAC_VU) and near-infrared. Water
    ! vapour absorbs only in the near-IR, via a two-exponential band; Rayleigh
    ! scattering (SW_RSCAT) and ozone (SW_C_ITF_O) act in the visible. Clean
    ! sky only in this first pass -- aerosol and cloud branches dropped
    ! (aerosol_ot = 0 collapses SESAM's scattering albedo to SW_RSCAT).

    real(wp), parameter :: SW_FRAC_VU = 0.45_wp    ! visible+UV fraction of TSI
    real(wp), parameter :: SW_RSCAT   = 0.17_wp    ! Rayleigh visible reflectance
    real(wp), parameter :: SW_C_ITF_O = 0.98_wp    ! ozone visible transmission
    real(wp), parameter :: SW_A1_W    = 0.21_wp    ! near-IR H2O band, weights
    real(wp), parameter :: SW_A2_W    = 1.0_wp - SW_A1_W
    real(wp), parameter :: SW_B1_W    = 6.27_wp    ! near-IR H2O band, exponents
    real(wp), parameter :: SW_B2_W    = 0.0267_wp
    real(wp), parameter :: SW_COSZ_O  = 1.0_wp/1.66_wp  ! diffuse-beam cosine

    ! === Grey-cloud shortwave optics =======================================
    ! Cloud SW optical depth per layer from the in-cloud water paths by
    ! geometric optics tau = 1.5 WP/(rho r_e); a conservative-scattering
    ! two-stream then gives the cloud reflectance R = gamma tau/(1+gamma tau),
    ! gamma = (1-g)/(2 mu), plus a small (near-IR) cloud absorptance. Standard
    ! droplet/crystal sizes and densities; tunable in one place. Unlike the LW
    ! grey absorption, SW is scattering, so the cloud enters as an albedo and the
    ! column is run clear + overcast and blended by cloud fraction (max overlap).
    real(wp), parameter :: SW_R_LIQ   = 10.0e-6_wp  ! liquid effective radius [m]
    real(wp), parameter :: SW_R_ICE   = 30.0e-6_wp  ! ice effective radius [m]
    real(wp), parameter :: SW_RHO_LIQ = 1000.0_wp   ! liquid water density [kg m-3]
    real(wp), parameter :: SW_RHO_ICE = 917.0_wp    ! ice density [kg m-3]
    real(wp), parameter :: SW_CLD_G   = 0.85_wp     ! cloud asymmetry factor [-]
    real(wp), parameter :: SW_CLD_ABS = 0.08_wp     ! max cloud SW absorptance [-]
    real(wp), parameter :: SW_CLD_TAU_A = 8.0_wp    ! absorptance e-folding optical depth
    real(wp), parameter :: SW_CLD_MAX = 0.999_wp    ! max column cloud fraction [-]

    ! === Orbit (present-day, for the stopgap insolation) ===================
    real(wp), parameter :: OBLIQUITY  = 23.44_wp*pi/180.0_wp  ! [rad]
    real(wp), parameter :: DAYS_YEAR  = 365.25_wp
    real(wp), parameter :: DOY_VE     = 80.0_wp    ! vernal equinox day-of-year

    ! === Prescribed ozone (analytic, zonally uniform) ======================
    ! A lognormal-in-pressure profile peaking in the stratosphere. Crude but
    ! enough to give the model top an ozone shortwave heating that balances its
    ! longwave cooling; without it a clear-sky top over-cools without bound.
    ! Zonal/seasonal ozone structure waits for an ERA5 ozone field.
    real(wp), parameter :: O3_MAX   = 1.5e-5_wp    ! peak mass mixing ratio [kg kg-1]
    real(wp), parameter :: O3_PPEAK = 2000.0_wp    ! peak pressure [Pa] (~20 hPa)
    real(wp), parameter :: O3_WIDTH = 1.2_wp       ! lognormal width in ln(p)

    ! === Configuration / state =============================================
    type, public :: aeros_rad_class
        logical :: enabled = .FALSE.        ! off by default, like convection
        ! ecCKD (correlated-k, §28) is the production default: better clear-sky OLR,
        ! CO2 forcing and LW cloud effect than SESAM, with no vapour-opacity fudge.
        ! SESAM (SCHEME_SESAM) is retained as the fast, bit-reproducible fallback.
        integer :: scheme  = SCHEME_ECCKD

        integer :: nlon = 0
        integer :: nlat = 0

        real(wp) :: co2_ppm  = 280.0_wp     ! CO2 volume mixing ratio [ppmv]
        logical  :: l_o3     = .TRUE.        ! prescribed ozone (LW + stratospheric SW)
        logical  :: clouds   = .FALSE.       ! diagnostic all-sky clouds (off => clear-sky)
        real(wp) :: q_co2    = 0.0_wp        ! CO2 mass mixing ratio [kg kg-1], derived

        ! Insolation / shortwave.
        real(wp) :: tsi      = S0            ! total solar irradiance [W m-2]
        real(wp) :: albedo   = 0.06_wp       ! surface broadband albedo (ocean)
        logical  :: seasonal = .FALSE.       ! .FALSE. = annual-mean insolation
        real(wp) :: doy0     = 0.0_wp        ! start day-of-year (seasonal mode)
        real(dp) :: time_bp  = 0.0_dp        ! orbital year before present (Laskar); 0 = present-day
        type(aeros_insol_class) :: ins       ! insol (Laskar 2004) insolation state

        ! Call cadence: recompute the full transfer every `interval` seconds and
        ! hold the heating rate fixed between, per design.md section 5.
        real(wp) :: interval = 10800.0_wp    ! 3 h

        ! Per-latitude insolation (nlat): annual-mean, or refreshed per call in
        ! seasonal mode.
        real(wp), allocatable :: sw_toa(:)   ! daily-mean TOA down SW [W m-2]
        real(wp), allocatable :: coszen(:)   ! airmass cosine zenith [-]

        ! Cached heating rate [K s-1], (nlon,nlat,nlev), applied every step.
        real(wp), allocatable :: heat(:,:,:)

        ! Diagnostics from the last recompute, (nlon,nlat) [W m-2].
        real(wp), allocatable :: olr(:,:)        ! outgoing LW at TOA
        real(wp), allocatable :: lw_dw_sur(:,:)  ! surface downward LW
        real(wp), allocatable :: sw_dw_sur(:,:)  ! surface downward SW
        real(wp), allocatable :: sw_net_sur(:,:) ! surface net absorbed SW
        real(wp), allocatable :: sw_up_toa(:,:)  ! reflected SW at TOA
    end type aeros_rad_class

    public :: aeros_radiation_init
    public :: aeros_radiation_load
    public :: aeros_radiation_end
    public :: aeros_radiation_apply
    public :: aeros_radiation_report
    public :: aeros_lw_clearsky_column
    public :: aeros_lw_cloudy_column
    public :: aeros_sw_clearsky_column
    public :: aeros_sw_cloudy_column
    public :: aeros_insolation_daily

contains

    subroutine aeros_radiation_init(rad, grd, enabled)
        ! Geometry, the derived CO2 mass ratio, and the per-latitude insolation.
        ! The 3D heating cache is allocated lazily on the first apply, where the
        ! level count is known.

        implicit none
        type(aeros_rad_class),   intent(inout) :: rad
        type(aeros_grid_class),  intent(in)    :: grd
        logical,                 intent(in)    :: enabled

        call aeros_radiation_end(rad)

        rad%enabled = enabled
        rad%nlon    = grd%nlon
        rad%nlat    = grd%nlat
        rad%q_co2   = co2_mass_ratio(rad%co2_ppm)

        allocate(rad%sw_toa(grd%nlat), rad%coszen(grd%nlat))
        allocate(rad%olr(grd%nlon, grd%nlat), rad%lw_dw_sur(grd%nlon, grd%nlat))
        allocate(rad%sw_dw_sur(grd%nlon, grd%nlat), rad%sw_up_toa(grd%nlon, grd%nlat))
        allocate(rad%sw_net_sur(grd%nlon, grd%nlat))
        rad%olr = 0.0_wp; rad%lw_dw_sur = 0.0_wp
        rad%sw_dw_sur = 0.0_wp; rad%sw_up_toa = 0.0_wp; rad%sw_net_sur = 0.0_wp

        ! Insolation from the insol package (Laskar 2004 orbit at rad%time_bp).
        ! The annual-mean per-latitude sw_toa/coszen seed the annual-mean mode;
        ! in seasonal mode they are overwritten each recompute from the current
        ! day (aeros_radiation_apply). coszen is the insolation-weighted mean
        ! airmass cosine, as the shortwave band scheme expects.
        call aeros_insol_init(rad%ins, grd, rad%time_bp)
        call aeros_insol_annual(rad%ins, rad%sw_toa, rad%coszen)

        ! Build the ecCKD reference k-table once, serially, before any OpenMP
        ! region (the per-column recompute then only reads the cached table).
        call aeros_ecckd_init()

        return
    end subroutine aeros_radiation_init

    subroutine aeros_radiation_load(rad, filename, grd, defaults_file)
        ! Read the `radiation` namelist group, then finish init.

        implicit none
        type(aeros_rad_class),   intent(inout) :: rad
        character(len=*),        intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_grid_class),  intent(in)    :: grd

        logical  :: enabled, seasonal
        integer  :: scheme
        real(wp) :: co2_ppm, albedo, tsi, interval, doy0
        real(dp) :: time_bp
        logical  :: l_o3, clouds

        enabled  = rad%enabled
        scheme   = rad%scheme
        co2_ppm  = rad%co2_ppm
        l_o3     = rad%l_o3
        clouds   = rad%clouds
        albedo   = rad%albedo
        tsi      = rad%tsi
        seasonal = rad%seasonal
        interval = rad%interval
        doy0     = rad%doy0
        time_bp  = rad%time_bp

        call nml_read(filename, "radiation", "enabled",  enabled, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "scheme",   scheme, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "co2_ppm",  co2_ppm, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "l_o3",     l_o3, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "clouds",   clouds, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "albedo",   albedo, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "tsi",      tsi, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "seasonal", seasonal, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "interval", interval, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "doy0",     doy0, defaults_file=defaults_file)
        call nml_read(filename, "radiation", "time_bp",  time_bp, defaults_file=defaults_file)

        rad%scheme   = scheme
        rad%co2_ppm  = co2_ppm
        rad%l_o3     = l_o3
        rad%clouds   = clouds
        rad%albedo   = albedo
        rad%tsi      = tsi
        rad%seasonal = seasonal
        rad%interval = interval
        rad%doy0     = doy0
        rad%time_bp  = time_bp

        call aeros_radiation_init(rad, grd, enabled)

        return
    end subroutine aeros_radiation_load

    subroutine aeros_radiation_end(rad)
        implicit none
        type(aeros_rad_class), intent(inout) :: rad
        if (allocated(rad%sw_toa))    deallocate(rad%sw_toa)
        if (allocated(rad%coszen))    deallocate(rad%coszen)
        if (allocated(rad%heat))      deallocate(rad%heat)
        if (allocated(rad%olr))       deallocate(rad%olr)
        if (allocated(rad%lw_dw_sur)) deallocate(rad%lw_dw_sur)
        if (allocated(rad%sw_dw_sur)) deallocate(rad%sw_dw_sur)
        if (allocated(rad%sw_net_sur)) deallocate(rad%sw_net_sur)
        if (allocated(rad%sw_up_toa)) deallocate(rad%sw_up_toa)
        call aeros_insol_end(rad%ins)
        rad%enabled = .FALSE.
        rad%nlon = 0; rad%nlat = 0
        return
    end subroutine aeros_radiation_end

    pure function co2_mass_ratio(co2_ppm) result(q_co2)
        ! ppmv -> kg kg-1
        implicit none
        real(wp), intent(in) :: co2_ppm
        real(wp) :: q_co2
        q_co2 = co2_ppm*1.0e-6_wp * M_CO2/M_AIR
        return
    end function co2_mass_ratio

    pure subroutine aeros_ozone_profile(pfull, o3)
        ! Analytic ozone mass mixing ratio [kg kg-1] on the full levels, from
        ! the lognormal-in-pressure profile above. Zonally uniform.
        implicit none
        real(wp), intent(in)  :: pfull(:)    ! (nlev) [Pa]
        real(wp), intent(out) :: o3(:)       ! (nlev) [kg kg-1]
        real(wp) :: x
        integer  :: k
        do k = 1, size(pfull)
            x = log(pfull(k)/O3_PPEAK)/O3_WIDTH
            o3(k) = O3_MAX*exp(-0.5_wp*x*x)
        end do
        return
    end subroutine aeros_ozone_profile

    subroutine aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                        q_co2, l_o3, fnet, heat, olr, fdw_sur)
        ! Clear-sky longwave for one column, model ordering (k=1 top ..
        ! k=nlev surface). Returns the net-upward LW flux at every interface,
        ! the layer heating rate, and the two boundary fluxes of interest.
        !
        ! Ported band physics (SESAM): blackbody source per layer, cumulative
        ! water-vapour / CO2 / ozone paths, the d_vap/d_co2/d_o3 transmission
        ! fits, and the emissivity-method flux integral. Fluxes live on the
        ! nlev+1 interfaces; layer absorber amounts and B on the nlev layers.
        !
        ! Interface convention (this routine, local): i = 0 top of atmosphere,
        ! i = nlev surface. Layer k lies between interface k-1 (above) and
        ! interface k (below). z_half and dp_lev come from the caller's vgrid.

        implicit none
        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: t(:)        ! (nlev) layer temperature [K]
        real(wp), intent(in)  :: q(:)        ! (nlev) specific humidity [kg kg-1]
        real(wp), intent(in)  :: o3(:)       ! (nlev) ozone mass mixing ratio [kg kg-1]
        real(wp), intent(in)  :: dp_lev(:)   ! (nlev) layer thickness [Pa], > 0
        real(wp), intent(in)  :: z_half(0:)  ! (0:nlev) interface height [m], z_half(nlev)=surface
        real(wp), intent(in)  :: ts          ! surface skin temperature [K]
        real(wp), intent(in)  :: q_co2       ! CO2 mass mixing ratio [kg kg-1]
        logical,  intent(in)  :: l_o3        ! include ozone

        real(wp), intent(out) :: fnet(0:)    ! (0:nlev) net UPWARD LW flux [W m-2]
        real(wp), intent(out) :: heat(:)     ! (nlev) LW heating rate [K s-1]
        real(wp), intent(out) :: olr         ! outgoing LW at TOA [W m-2]
        real(wp), intent(out) :: fdw_sur     ! downward LW at surface [W m-2]

        ! layer quantities, in local surface->top order (l = 1 surface .. nlev top)
        real(wp) :: b(nlev)                  ! blackbody source sigma T^4
        real(wp) :: uwv(nlev), uco2(nlev), uo3(nlev)  ! layer absorber amounts
        real(wp) :: zc(0:nlev)               ! interface heights [cm], zc(0)=surface
        real(wp) :: bsfc, fup(0:nlev), fdw(0:nlev)
        real(wp) :: sw, sc, so, tau
        real(wp) :: h0, kco2, expc(0:nlev), zmid, dz
        integer  :: k, l, i, kk

        h0   = (R_d*T0/grav)*100.0_wp        ! atmospheric scale height [cm]
        kco2 = (LW_AK_CO2 + 1.0_wp)/h0       ! [1/cm]

        ! Interface heights in cm, local surface->top order. Model k=nlev is the
        ! surface layer, whose lower interface z_half(nlev) is the ground.
        do i = 0, nlev
            zc(i) = z_half(nlev - i)*100.0_wp
            expc(i) = exp(-kco2*zc(i))
        end do

        ! Layer quantities in local order: local layer l corresponds to model
        ! layer k = nlev - l + 1, bounded below by interface l-1, above by l.
        do l = 1, nlev
            k = nlev - l + 1
            b(l) = LW_EMIS*sigma_sb*t(k)**4

            ! water vapour: exact hydrostatic path 0.1 q dp/g  [g cm-2]
            uwv(l) = 0.1_wp * q(k) * dp_lev(k)/grav

            ! CO2: SESAM pressure-broadened path (dimensionless mixing ratio x
            ! length scale), evaluated between the layer's two interfaces.
            uco2(l) = q_co2/kco2 * (expc(l-1) - expc(l))

            ! ozone: SESAM pressure-broadened path [g cm-2]
            if (l_o3) then
                zmid = 0.5_wp*(zc(l-1) + zc(l))
                dz   = zc(l) - zc(l-1)
                uo3(l) = (p0/(R_d*T0))*1.0e-3_wp * exp(-zmid*(LW_AK_O3+1.0_wp)/h0) &
                         * o3(k) * dz
            else
                uo3(l) = 0.0_wp
            end if
        end do

        bsfc = LW_EMIS*sigma_sb*ts**4

        ! --- Upward flux at each interface i (positive up) -----------------
        ! F_up(i) = Bsfc T(sfc->i) + sum over layers l at/below i of
        !           B_l ( T(top_l -> i) - T(bot_l -> i) ),
        ! T built from the cumulative path of the layers strictly between the
        ! two interfaces. Surface is interface 0.
        fup(0) = bsfc
        do i = 1, nlev
            fup(i) = bsfc*trans(0, i)
            do l = 1, i
                fup(i) = fup(i) + b(l)*(trans(l, i) - trans(l-1, i))
            end do
        end do

        ! --- Downward flux at each interface i (positive down) -------------
        ! No downward flux incident at TOA. F_dw(i) = sum over layers above i.
        fdw(nlev) = 0.0_wp
        do i = nlev-1, 0, -1
            fdw(i) = 0.0_wp
            do l = i+1, nlev
                fdw(i) = fdw(i) + b(l)*(trans(l-1, i) - trans(l, i))
            end do
        end do

        ! --- Net upward flux on interfaces, back to model ordering ---------
        ! Local interface i (surface->top) is model interface nlev-i.
        do i = 0, nlev
            fnet(nlev - i) = fup(i) - fdw(i)
        end do

        olr     = fup(nlev)     ! TOA, positive up
        fdw_sur = fdw(0)        ! surface downwelling

        ! --- Layer heating: divergence of net-upward flux ------------------
        ! Model layer k bounded by interfaces k-1 (above) and k (below).
        ! Absorbed energy = F_net_up(below) - F_net_up(above); dp>0.
        do k = 1, nlev
            heat(k) = (grav/cp_d) * (fnet(k) - fnet(k-1))/dp_lev(k)
        end do

        return

    contains

        pure real(wp) function trans(ia, ib) result(tr)
            ! Broadband transmission between local interfaces ia and ib, from
            ! the cumulative absorber path of the layers strictly between them.
            integer, intent(in) :: ia, ib
            integer :: lo, hi, m
            real(wp) :: aw, ac, ao, dv, dc, dobn

            lo = min(ia, ib); hi = max(ia, ib)
            if (hi - lo <= 0) then
                tr = 1.0_wp
                return
            end if

            aw = 0.0_wp; ac = 0.0_wp; ao = 0.0_wp
            do m = lo+1, hi
                aw = aw + uwv(m)
                ac = ac + uco2(m)
                ao = ao + uo3(m)
            end do

            ! water vapour, PIK report 81 eq. 6.5
            dv = 1.0_wp/(1.0_wp + LW_VAP_OPAC*( &
                                  LW_A_VAP *(LW_BETA0*aw)**LW_BETA_VAP  &
                                + LW_A2_VAP*(LW_BETA0*aw)**LW_BETA2_VAP &
                                + LW_A3_VAP*(LW_BETA0*aw)**3))
            ! CO2, PIK report 81 eq. 6.6 with the high-CO2 correction factor
            dc = (1.0_wp - min(0.2_wp, 0.1_wp*(ac/1000.0_wp)**2)) &
                 * (1.0_wp + LW_A0_CO2*LW_A1_CO2*(LW_BETA0*ac)**LW_BETA_CO2) &
                 / (1.0_wp + LW_A0_CO2*(LW_BETA0*ac)**LW_BETA_CO2)
            ! ozone, PIK report 81 eq. 6.7
            dobn = 1.0_wp - LW_A_O3*(ao**LW_BETA_O3)

            tr = dv*dc*dobn
            return
        end function trans

    end subroutine aeros_lw_clearsky_column

    subroutine aeros_lw_cloudy_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      q_co2, l_o3, cf, clwc, ciwc, &
                                      fnet, heat, olr, fdw_sur)
        ! All-sky longwave for one column by SESAM's run-twice-and-blend
        ! (lwr_total): a clear column and an overcast column are computed with the
        ! same band kernel and blended by the column cloud fraction at maximum
        ! overlap,
        !
        !   F = (1 - CF) F_clear + CF F_overcast,   CF = max_k cf_k.
        !
        ! The overcast column carries the cloud as a per-layer grey transmission
        ! tcl_k = exp(-(k_liq LWP_k + k_ice IWP_k)) built from the IN-CLOUD water
        ! paths LWP_k = 1e3 clwc_k dp_k/g / CF (the grid-mean condensate spread
        ! over the cloud fraction, so the blend conserves the grid-mean water),
        ! multiplied into the clear-sky gas transmission between interfaces (grey
        ! absorbers multiply; the gas band fit is still evaluated once on the
        ! accumulated path). Maximum overlap matches the shortwave kernel and
        ! avoids the overcast bias of per-layer random overlap (which drives the
        ! column toward 1 - prod(1-cf_k) and over-traps). cf=0 (or zero
        ! condensate) returns the clear-sky column bit-for-bit.
        !
        ! Cloud optics are intensive, so the routine is grid-agnostic: it runs on
        ! ERA5's pressure levels (the all-sky validation) or a refined radiation
        ! grid remapped to the transport grid, unchanged -- the grid choice lives
        ! in the caller.

        implicit none
        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: t(:)        ! (nlev) layer temperature [K]
        real(wp), intent(in)  :: q(:)        ! (nlev) specific humidity [kg kg-1]
        real(wp), intent(in)  :: o3(:)       ! (nlev) ozone mass mixing ratio [kg kg-1]
        real(wp), intent(in)  :: dp_lev(:)   ! (nlev) layer thickness [Pa], > 0
        real(wp), intent(in)  :: z_half(0:)  ! (0:nlev) interface height [m]
        real(wp), intent(in)  :: ts          ! surface skin temperature [K]
        real(wp), intent(in)  :: q_co2       ! CO2 mass mixing ratio [kg kg-1]
        logical,  intent(in)  :: l_o3        ! include ozone
        real(wp), intent(in)  :: cf(:)       ! (nlev) cloud fraction [0-1]
        real(wp), intent(in)  :: clwc(:)     ! (nlev) cloud liquid water [kg kg-1]
        real(wp), intent(in)  :: ciwc(:)     ! (nlev) cloud ice water    [kg kg-1]

        real(wp), intent(out) :: fnet(0:)    ! (0:nlev) net UPWARD LW flux [W m-2]
        real(wp), intent(out) :: heat(:)     ! (nlev) LW heating rate [K s-1]
        real(wp), intent(out) :: olr         ! outgoing LW at TOA [W m-2]
        real(wp), intent(out) :: fdw_sur     ! downward LW at surface [W m-2]

        ! clear-sky reference
        real(wp) :: fnet_cs(0:nlev), heat_cs(nlev), olr_cs, fdw_cs
        ! overcast column, local surface->top order (l = 1 surface .. nlev top)
        real(wp) :: b(nlev), uwv(nlev), uco2(nlev), uo3(nlev), tcl(nlev)
        real(wp) :: zc(0:nlev), expc(0:nlev)
        real(wp) :: bsfc, fup(0:nlev), fdw(0:nlev)
        real(wp) :: h0, kco2, zmid, dz, lwp, iwp, cff, olr_ov, fdw_ov
        integer  :: k, l, i

        ! clear-sky column (the reference; also the no-cloud short-circuit)
        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
            q_co2, l_o3, fnet_cs, heat_cs, olr_cs, fdw_cs)

        ! column cloud fraction, maximum overlap
        cff = 0.0_wp
        do k = 1, nlev
            cff = max(cff, min(1.0_wp, max(0.0_wp, cf(k))))
        end do
        if (cff <= CLD_CF_FLOOR) then
            fnet = fnet_cs; heat = heat_cs; olr = olr_cs; fdw_sur = fdw_cs
            return
        end if

        ! --- overcast column: gas kernel x grey cloud (in-cloud water/CF) ------
        h0   = (R_d*T0/grav)*100.0_wp
        kco2 = (LW_AK_CO2 + 1.0_wp)/h0
        do i = 0, nlev
            zc(i) = z_half(nlev - i)*100.0_wp
            expc(i) = exp(-kco2*zc(i))
        end do
        do l = 1, nlev
            k = nlev - l + 1
            b(l) = LW_EMIS*sigma_sb*t(k)**4
            uwv(l) = 0.1_wp * q(k) * dp_lev(k)/grav
            uco2(l) = q_co2/kco2 * (expc(l-1) - expc(l))
            if (l_o3) then
                zmid = 0.5_wp*(zc(l-1) + zc(l))
                dz   = zc(l) - zc(l-1)
                uo3(l) = (p0/(R_d*T0))*1.0e-3_wp * exp(-zmid*(LW_AK_O3+1.0_wp)/h0) &
                         * o3(k) * dz
            else
                uo3(l) = 0.0_wp
            end if
            ! in-cloud condensate = grid-mean/CF; overcast grey transmission
            lwp = 1.0e3_wp * max(0.0_wp, clwc(k)) * dp_lev(k)/grav / cff
            iwp = 1.0e3_wp * max(0.0_wp, ciwc(k)) * dp_lev(k)/grav / cff
            tcl(l) = exp(-(LW_KABS_LIQ*lwp + LW_KABS_ICE*iwp))
        end do

        bsfc = LW_EMIS*sigma_sb*ts**4
        fup(0) = bsfc
        do i = 1, nlev
            fup(i) = bsfc*trans(0, i)
            do l = 1, i
                fup(i) = fup(i) + b(l)*(trans(l, i) - trans(l-1, i))
            end do
        end do
        fdw(nlev) = 0.0_wp
        do i = nlev-1, 0, -1
            fdw(i) = 0.0_wp
            do l = i+1, nlev
                fdw(i) = fdw(i) + b(l)*(trans(l-1, i) - trans(l, i))
            end do
        end do
        olr_ov = fup(nlev)
        fdw_ov = fdw(0)

        ! --- blend clear and overcast by CF (maximum overlap) -----------------
        ! overcast net-upward on model interface i is fup(nlev-i) - fdw(nlev-i)
        do i = 0, nlev
            fnet(i) = (1.0_wp-cff)*fnet_cs(i) + cff*(fup(nlev-i) - fdw(nlev-i))
        end do
        olr     = (1.0_wp-cff)*olr_cs + cff*olr_ov
        fdw_sur = (1.0_wp-cff)*fdw_cs + cff*fdw_ov
        do k = 1, nlev
            heat(k) = (grav/cp_d) * (fnet(k) - fnet(k-1))/dp_lev(k)
        end do

        return

    contains

        pure real(wp) function trans(ia, ib) result(tr)
            ! Broadband gas transmission between local interfaces ia and ib
            ! (clear-sky d_vap d_co2 d_o3 fits on the accumulated path) times the
            ! product of the grey cloud transmissions of the layers between them.
            integer, intent(in) :: ia, ib
            integer :: lo, hi, m
            real(wp) :: aw, ac, ao, dv, dc, dobn, tcld

            lo = min(ia, ib); hi = max(ia, ib)
            if (hi - lo <= 0) then
                tr = 1.0_wp
                return
            end if

            aw = 0.0_wp; ac = 0.0_wp; ao = 0.0_wp; tcld = 1.0_wp
            do m = lo+1, hi
                aw = aw + uwv(m)
                ac = ac + uco2(m)
                ao = ao + uo3(m)
                tcld = tcld * tcl(m)
            end do

            dv = 1.0_wp/(1.0_wp + LW_VAP_OPAC*( &
                                  LW_A_VAP *(LW_BETA0*aw)**LW_BETA_VAP  &
                                + LW_A2_VAP*(LW_BETA0*aw)**LW_BETA2_VAP &
                                + LW_A3_VAP*(LW_BETA0*aw)**3))
            dc = (1.0_wp - min(0.2_wp, 0.1_wp*(ac/1000.0_wp)**2)) &
                 * (1.0_wp + LW_A0_CO2*LW_A1_CO2*(LW_BETA0*ac)**LW_BETA_CO2) &
                 / (1.0_wp + LW_A0_CO2*(LW_BETA0*ac)**LW_BETA_CO2)
            dobn = 1.0_wp - LW_A_O3*(ao**LW_BETA_O3)

            tr = dv*dc*dobn*tcld
            return
        end function trans

    end subroutine aeros_lw_cloudy_column

    subroutine aeros_insolation_daily(lat, doy, tsi, coszen, swdown)
        ! Daily-mean insolation for one latitude and day-of-year. Present-day
        ! orbit, circular, obliquity-only -- a stopgap for the eventual Laskar
        ! (2004) `insol` forcing (design.md section 8). Returns the daily-mean
        ! TOA downward SW on a horizontal surface, and the daylight-weighted
        ! mean cosine zenith the shortwave uses as its airmass. No diurnal
        ! cycle: radiation is called on a multi-hour cadence and there is no
        ! surface heat store yet to lag, so the daily mean is the honest input.

        implicit none
        real(wp), intent(in)  :: lat        ! latitude [rad]
        real(wp), intent(in)  :: doy        ! day of year [1..365.25]
        real(wp), intent(in)  :: tsi        ! solar constant [W m-2]
        real(wp), intent(out) :: coszen     ! daylight-mean cos(zenith) [-]
        real(wp), intent(out) :: swdown     ! daily-mean TOA down SW [W m-2]

        real(wp) :: decl, lam, cosh0, h0, sinterm

        ! solar declination, circular-orbit approximation
        lam  = 2.0_wp*pi*(doy - DOY_VE)/DAYS_YEAR
        decl = asin(sin(OBLIQUITY)*sin(lam))

        ! sunset hour angle, clipped for polar day / night
        cosh0 = -tan(lat)*tan(decl)
        if (cosh0 >= 1.0_wp) then
            h0 = 0.0_wp                       ! polar night
        else if (cosh0 <= -1.0_wp) then
            h0 = pi                           ! polar day
        else
            h0 = acos(cosh0)
        end if

        ! sinterm is proportional to both the daily-mean insolation and, divided
        ! by h0, the daylight-mean cosine zenith.
        sinterm = h0*sin(lat)*sin(decl) + cos(lat)*cos(decl)*sin(h0)
        sinterm = max(0.0_wp, sinterm)

        swdown = (tsi/pi)*sinterm
        if (h0 > 0.0_wp .and. swdown > 0.0_wp) then
            coszen = sinterm/h0
        else
            coszen = 0.0_wp
        end if

        return
    end subroutine aeros_insolation_daily

    subroutine aeros_sw_clearsky_column(nlev, q, o3, l_o3, dp_lev, swdown_toa, &
                                        coszen, alb_vis, alb_ir, heat, &
                                        sw_up_toa, sw_dw_sur, sw_net_sur)
        ! Clear-sky shortwave for one column, model ordering (k=1 top ..
        ! k=nlev surface). Two bands; near-IR water-vapour absorption and
        ! visible ozone absorption are the atmospheric heating, Rayleigh acts
        ! in the visible.
        !
        ! The boundary quantities -- planetary albedo (TOA up) and the surface
        ! net/downward flux -- are SESAM's clear-sky formulas verbatim (its
        ! tuned two-path surface transmission). The column atmospheric
        ! absorption is then the TOA-minus-surface residual, split into a
        ! near-IR water part (distributed by the water direct-beam absorption
        ! weight) and a visible ozone part (distributed by the ozone amount).
        ! Depositing the ozone absorption where the ozone is -- the stratosphere
        ! -- is what balances the longwave cooling of the model top; lumping it
        ! near the surface (with the water) leaves the top with no heating and
        ! it over-cools. The upward-beam re-absorption is folded into the
        ! boundary residual, not resolved per layer -- a few W/m2, second order.

        implicit none
        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: q(:)        ! (nlev) specific humidity [kg kg-1]
        real(wp), intent(in)  :: o3(:)       ! (nlev) ozone mass mixing ratio [kg kg-1]
        logical,  intent(in)  :: l_o3        ! resolve ozone SW heating
        real(wp), intent(in)  :: dp_lev(:)   ! (nlev) layer thickness [Pa], > 0
        real(wp), intent(in)  :: swdown_toa  ! daily-mean TOA down SW [W m-2]
        real(wp), intent(in)  :: coszen      ! airmass cosine zenith [-]
        real(wp), intent(in)  :: alb_vis     ! surface albedo, visible [-]
        real(wp), intent(in)  :: alb_ir      ! surface albedo, near-IR [-]

        real(wp), intent(out) :: heat(:)     ! (nlev) SW heating rate [K s-1]
        real(wp), intent(out) :: sw_up_toa   ! TOA upward (reflected) SW [W m-2]
        real(wp), intent(out) :: sw_dw_sur   ! surface downward SW [W m-2]
        real(wp), intent(out) :: sw_net_sur  ! surface net absorbed SW [W m-2]

        real(wp) :: cosz, icos
        real(wp) :: uwv_cum(0:nlev), fdw_ir(0:nlev), w(nlev), wo3(nlev)
        real(wp) :: m_col, itf_ir_full, alb_vu, alb_ir_p, alb_p
        real(wp) :: m_d1, m_d2, itf_d1, itf_d2, itf_atm_vu, itf_atm_ir
        real(wp) :: a_atm, a_o3, a_w, wsum, wo3sum, m
        integer  :: k, i

        heat = 0.0_wp
        if (swdown_toa <= 0.0_wp) then          ! polar night
            sw_up_toa = 0.0_wp; sw_dw_sur = 0.0_wp; sw_net_sur = 0.0_wp
            return
        end if

        cosz = max(coszen, 0.1_wp)
        icos = 1.0_wp/cosz + 1.0_wp/SW_COSZ_O   ! down direct + diffuse up airmass

        ! cumulative water-vapour path from TOA to each interface [g cm-2]
        uwv_cum(0) = 0.0_wp
        do k = 1, nlev
            uwv_cum(k) = uwv_cum(k-1) + 0.1_wp*q(k)*dp_lev(k)/grav
        end do

        ! --- planetary albedo (SESAM clear-sky, aerosol off) ---------------
        m_col = uwv_cum(nlev)*icos
        itf_ir_full = itf_w_ir(m_col)
        alb_vu   = (SW_RSCAT + (1.0_wp-SW_RSCAT)**2*alb_vis &
                    /(1.0_wp - SW_RSCAT*alb_vis))*SW_C_ITF_O
        alb_ir_p = alb_ir*itf_ir_full
        alb_p    = SW_FRAC_VU*alb_vu + (1.0_wp-SW_FRAC_VU)*alb_ir_p
        sw_up_toa = swdown_toa*alb_p

        ! --- surface net and downward (SESAM two-path clear-sky) -----------
        m_d1 = uwv_cum(nlev)/cosz
        m_d2 = m_d1 + uwv_cum(nlev)*(1.0_wp - exp(-0.25_wp))*2.0_wp/SW_COSZ_O
        itf_d1 = itf_w_ir(m_d1)
        itf_d2 = itf_w_ir(m_d2)
        ! visible: direct + one surface<->Rayleigh reflection, times ozone
        itf_atm_vu = (1.0_wp-SW_RSCAT)*(1.0_wp-alb_vis)*SW_C_ITF_O &
                   + (1.0_wp-SW_RSCAT)*alb_vis*SW_RSCAT*(1.0_wp-alb_vis) &
                     /(1.0_wp - SW_RSCAT*alb_vis)*SW_C_ITF_O
        ! near-IR: no atmospheric scattering, so only the direct absorbed beam
        itf_atm_ir = (1.0_wp-alb_ir)*itf_d1
        sw_net_sur = swdown_toa*(SW_FRAC_VU*itf_atm_vu + (1.0_wp-SW_FRAC_VU)*itf_atm_ir)
        sw_dw_sur  = swdown_toa*(SW_FRAC_VU*itf_atm_vu/(1.0_wp-alb_vis) &
                               + (1.0_wp-SW_FRAC_VU)*itf_atm_ir/(1.0_wp-alb_ir))

        ! --- distribute column absorption over layers ----------------------
        a_atm = max(0.0_wp, swdown_toa - sw_up_toa - sw_net_sur)

        ! Split off the visible ozone absorption (the post-Rayleigh visible beam
        ! that ozone takes) so it can be deposited where the ozone is. The
        ! remainder is the near-IR water absorption. When ozone is not resolved,
        ! a_o3 = 0 and everything rides the water weight, as before.
        if (l_o3) then
            a_o3 = swdown_toa*SW_FRAC_VU*(1.0_wp - SW_RSCAT)*(1.0_wp - SW_C_ITF_O)
            a_o3 = min(a_o3, a_atm)
        else
            a_o3 = 0.0_wp
        end if
        a_w = a_atm - a_o3

        ! near-IR water weight: the drop in direct-beam transmission per layer
        do i = 0, nlev
            m = uwv_cum(i)/cosz
            fdw_ir(i) = itf_w_ir(m)
        end do
        wsum = 0.0_wp
        do k = 1, nlev
            w(k) = max(0.0_wp, fdw_ir(k-1) - fdw_ir(k))
            wsum = wsum + w(k)
        end do
        if (wsum <= 0.0_wp) then                 ! dry column: weight by mass
            do k = 1, nlev
                w(k) = dp_lev(k); wsum = wsum + dp_lev(k)
            end do
        end if

        ! ozone weight: the ozone mass path per layer (stratosphere-peaked)
        wo3sum = 0.0_wp
        do k = 1, nlev
            wo3(k) = max(0.0_wp, o3(k))*dp_lev(k)
            wo3sum = wo3sum + wo3(k)
        end do
        if (wo3sum <= 0.0_wp) then               ! no ozone: fall back to water
            wo3 = w; wo3sum = wsum
        end if

        do k = 1, nlev
            heat(k) = (grav/cp_d)*(a_w*(w(k)/wsum) + a_o3*(wo3(k)/wo3sum))/dp_lev(k)
        end do

        return

    contains

        pure real(wp) function itf_w_ir(m) result(tr)
            ! near-IR water-vapour band transmission, path m [g cm-2]
            real(wp), intent(in) :: m
            tr = SW_A1_W*exp(-SW_B1_W*m) + SW_A2_W*exp(-SW_B2_W*m)
            return
        end function itf_w_ir

    end subroutine aeros_sw_clearsky_column

    subroutine aeros_sw_cloudy_column(nlev, q, o3, l_o3, dp_lev, swdown_toa, &
                                      coszen, alb_vis, alb_ir, cf, clwc, ciwc, &
                                      heat, sw_up_toa, sw_dw_sur, sw_net_sur)
        ! All-sky shortwave for one column. Shortwave is scattering, not grey
        ! absorption, so -- unlike the longwave -- the cloud cannot fold into a
        ! per-layer transmission. Instead the tuned clear-sky column
        ! (aeros_sw_clearsky_column) is run once, an overcast column is built by
        ! placing a cloud reflector on top of it, and the two are blended by the
        ! column cloud fraction (maximum overlap), exactly SESAM's run-twice
        ! structure:  F = (1-CF) F_clear + CF F_overcast.
        !
        ! The cloud reflector is resolved from the in-cloud water paths: each
        ! layer contributes tau_k = 1.5 (LWP_k/(rho_l r_l) + IWP_k/(rho_i r_i)),
        ! with the in-cloud path = grid-mean/CF so the blend conserves the
        ! grid-mean condensate. A conservative-scattering two-stream on the
        ! column optical depth tau = sum tau_k gives the cloud reflectance
        ! R = gamma tau/(1+gamma tau), gamma=(1-g)/(2 mu); a small saturating
        ! absorptance A takes the near-IR cloud absorption; T = 1-R-A. Adding the
        ! cloud reflector above the clear column (planetary albedo alb_clr):
        !
        !   alb_ov = R + T^2 alb_clr/(1 - R alb_clr),   Phi = T/(1 - R alb_clr)
        !
        ! so the overcast surface fluxes are the clear ones times Phi. The
        ! overcast atmospheric absorption is the exact TOA-minus-surface residual,
        ! distributed by the clear-sky heating profile (gas) plus the per-layer
        ! cloud optical depth (cloud), conserving energy. Optics are intensive,
        ! so the routine is grid-agnostic, as the clear-sky one.

        implicit none
        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: q(:)        ! (nlev) specific humidity [kg kg-1]
        real(wp), intent(in)  :: o3(:)       ! (nlev) ozone mass mixing ratio [kg kg-1]
        logical,  intent(in)  :: l_o3
        real(wp), intent(in)  :: dp_lev(:)   ! (nlev) layer thickness [Pa], > 0
        real(wp), intent(in)  :: swdown_toa  ! daily-mean TOA down SW [W m-2]
        real(wp), intent(in)  :: coszen      ! airmass cosine zenith [-]
        real(wp), intent(in)  :: alb_vis     ! surface albedo, visible [-]
        real(wp), intent(in)  :: alb_ir      ! surface albedo, near-IR [-]
        real(wp), intent(in)  :: cf(:)       ! (nlev) cloud fraction [0-1]
        real(wp), intent(in)  :: clwc(:)     ! (nlev) cloud liquid water [kg kg-1]
        real(wp), intent(in)  :: ciwc(:)     ! (nlev) cloud ice water    [kg kg-1]

        real(wp), intent(out) :: heat(:)     ! (nlev) SW heating rate [K s-1]
        real(wp), intent(out) :: sw_up_toa   ! TOA upward (reflected) SW [W m-2]
        real(wp), intent(out) :: sw_dw_sur   ! surface downward SW [W m-2]
        real(wp), intent(out) :: sw_net_sur  ! surface net absorbed SW [W m-2]

        real(wp) :: heat_cs(nlev), up_cs, dw_cs, net_cs
        real(wp) :: tau_c(nlev), cldw(nlev)
        real(wp) :: cff, cosz, gamma, tau, rcld, acld, tcld
        real(wp) :: alb_clr, alb_ov, phi, up_ov, dw_ov, net_ov, a_atm_ov
        real(wp) :: lwp, iwp, wsum, gsum
        integer  :: k

        ! clear-sky column (the tuned reference; also the polar-night guard)
        call aeros_sw_clearsky_column(nlev, q, o3, l_o3, dp_lev, swdown_toa, &
            coszen, alb_vis, alb_ir, heat_cs, up_cs, dw_cs, net_cs)

        ! column cloud fraction: maximum overlap
        cff = 0.0_wp
        do k = 1, nlev
            cff = max(cff, min(1.0_wp, max(0.0_wp, cf(k))))
        end do
        cff = min(cff, SW_CLD_MAX)

        ! no cloud or no sun: the clear-sky column is the answer
        if (cff <= CLD_CF_FLOOR .or. swdown_toa <= 0.0_wp) then
            heat = heat_cs
            sw_up_toa = up_cs; sw_dw_sur = dw_cs; sw_net_sur = net_cs
            return
        end if

        ! per-layer cloud optical depth from the in-cloud water paths
        ! (grid-mean/CF), geometric optics tau = 1.5 WP/(rho r_e)
        tau = 0.0_wp
        do k = 1, nlev
            lwp = max(0.0_wp, clwc(k))*dp_lev(k)/grav / cff      ! [kg m-2] in-cloud
            iwp = max(0.0_wp, ciwc(k))*dp_lev(k)/grav / cff
            tau_c(k) = 1.5_wp*(lwp/(SW_RHO_LIQ*SW_R_LIQ) + iwp/(SW_RHO_ICE*SW_R_ICE))
            tau = tau + tau_c(k)
        end do

        ! conservative-scattering two-stream cloud reflectance + small absorptance
        cosz  = max(coszen, 0.1_wp)
        gamma = 0.5_wp*(1.0_wp - SW_CLD_G)/cosz
        rcld  = gamma*tau/(1.0_wp + gamma*tau)
        acld  = SW_CLD_ABS*(1.0_wp - exp(-tau/SW_CLD_TAU_A))
        rcld  = min(rcld, 1.0_wp - acld)
        tcld  = max(0.0_wp, 1.0_wp - rcld - acld)

        ! add the cloud reflector above the clear column
        alb_clr = up_cs/swdown_toa
        alb_ov  = rcld + tcld*tcld*alb_clr/(1.0_wp - rcld*alb_clr)
        phi     = tcld/(1.0_wp - rcld*alb_clr)
        up_ov   = swdown_toa*alb_ov
        dw_ov   = dw_cs*phi
        net_ov  = net_cs*phi

        ! overcast atmospheric absorption = exact boundary residual, distributed
        ! by the clear-sky heating (gas) plus the cloud optical depth (cloud)
        a_atm_ov = max(0.0_wp, swdown_toa - up_ov - net_ov)
        wsum = 0.0_wp; gsum = 0.0_wp
        do k = 1, nlev
            cldw(k) = tau_c(k) + max(0.0_wp, heat_cs(k))*dp_lev(k)
            wsum = wsum + cldw(k)
            gsum = gsum + max(0.0_wp, heat_cs(k))*dp_lev(k)
        end do
        if (wsum <= 0.0_wp) then                 ! no weight: fall back to mass
            do k = 1, nlev
                cldw(k) = dp_lev(k); wsum = wsum + dp_lev(k)
            end do
        end if

        ! blend clear and overcast by the column cloud fraction
        sw_up_toa  = (1.0_wp-cff)*up_cs  + cff*up_ov
        sw_dw_sur  = (1.0_wp-cff)*dw_cs  + cff*dw_ov
        sw_net_sur = (1.0_wp-cff)*net_cs + cff*net_ov
        do k = 1, nlev
            heat(k) = (1.0_wp-cff)*heat_cs(k) &
                    + cff*(grav/cp_d)*(a_atm_ov*(cldw(k)/wsum))/dp_lev(k)
        end do

        return
    end subroutine aeros_sw_cloudy_column

    subroutine aeros_radiation_apply(rad, vg, grd, t_g, qv_g, lnps_g, t_s, &
                                     nstep, dt, dt_phys)
        ! Radiative heating at the grid seam. The full column transfer (LW + SW)
        ! is recomputed every `interval` seconds; between recomputes the cached
        ! heating rate is held fixed (design.md section 5).
        !
        ! The heating is applied FORWARD-SPLIT, as an increment [K] accumulated
        ! in wrk%dt_phys and added to the n+1 state after the dynamics step --
        ! the same path convection uses, NOT the centered leapfrog. Radiation is
        ! smooth in time but its vertical structure is sharp and, once the model
        ! top starts cooling, large; routed through the centered leapfrog it
        ! excites the computational mode and the run blows up (the handoff's
        ! "large and stiff" case). Forward-split is stable for a damping process
        ! like radiation at this timestep.
        !
        ! Surface (skin) temperature t_s is the prescribed SST; it is the LW
        ! lower boundary and, with the albedo, the SW surface. Ozone is the
        ! prescribed analytic profile when l_o3 -- its shortwave heating in the
        ! stratosphere is what keeps the model top from over-cooling.

        implicit none
        type(aeros_rad_class),   intent(inout) :: rad
        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_grid_class),  intent(in)    :: grd
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(in)    :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat) ln[Pa]
        real(wp), intent(in)    :: t_s(:,:)       ! (nlon,nlat) skin temp [K]
        integer,  intent(in)    :: nstep          ! step counter
        real(wp), intent(in)    :: dt             ! [s]
        real(wp), intent(inout) :: dt_phys(:,:,:) ! (nlon,nlat,nlev) [K] increment

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: phi_full(vg%nlev), phi_half(0:vg%nlev), z_half(0:vg%nlev)
        real(wp) :: o3col(vg%nlev)
        real(wp) :: cf(vg%nlev), clwc(vg%nlev), ciwc(vg%nlev)
        real(wp) :: fnet(0:vg%nlev), heat_lw(vg%nlev), heat_sw(vg%nlev)
        real(wp) :: olr, fdw_lw, sw_up, sw_dw, sw_net, doy, cz, sw
        integer  :: i, j, k, nlev, nrad

        if (.not. rad%enabled) return

        nlev = vg%nlev
        if (.not. allocated(rad%heat)) then
            allocate(rad%heat(rad%nlon, rad%nlat, nlev))
            rad%heat = 0.0_wp
        end if

        nrad = max(1, nint(rad%interval/dt))

        ! --- recompute the transfer on the cadence -------------------------
        if (mod(nstep, nrad) == 0) then

            if (rad%seasonal) then
                ! Refresh insolation for the current day from insol (Laskar orbit
                ! at rad%time_bp). aeros_insol_day wraps the day into the radiative
                ! year (DAY_YEAR), so pass the raw elapsed day-of-year.
                doy = rad%doy0 + real(nstep, wp)*dt/86400.0_wp
                call aeros_insol_day(rad%ins, rad%sw_toa, rad%coszen, doy)
            end if

            !$omp parallel do collapse(2) schedule(static) &
            !$omp   private(i,j,phalf,pfull,dpc,phi_full,phi_half,z_half,o3col, &
            !$omp           cf,clwc,ciwc,fnet,heat_lw,heat_sw,olr,fdw_lw,sw_up,sw_dw,sw_net)
            do j = 1, rad%nlat
                do i = 1, rad%nlon
                    call aeros_vgrid_pressure(vg, exp(lnps_g(i,j)), phalf, pfull, dpc)

                    if (rad%l_o3) then
                        call aeros_ozone_profile(pfull, o3col)
                    else
                        o3col = 0.0_wp
                    end if

                    ! interface heights (aquaplanet: surface geopotential 0)
                    call aeros_hydrostatic(vg, 0.0_wp, t_g(i,j,:), phalf, &
                                           phi_full, phi_half)
                    z_half = phi_half/grav

                    if (rad%clouds) then
                        ! diagnose the cloud column, then the all-sky operators
                        call aeros_cloud_diagnose(nlev, t_g(i,j,:), qv_g(i,j,:), &
                            pfull, exp(lnps_g(i,j)), cf, clwc, ciwc)

                        ! All-sky: dispatch on the scheme selector (opt-in ecCKD
                        ! folds the same grey cloud optics into its g-points).
                        if (rad%scheme == SCHEME_ECCKD) then
                            call aeros_ecckd_lw_cloudy_column(nlev, t_g(i,j,:), &
                                qv_g(i,j,:), o3col, dpc, z_half, t_s(i,j), rad%q_co2, &
                                rad%l_o3, cf, clwc, ciwc, fnet, heat_lw, olr, fdw_lw)
                            call aeros_ecckd_sw_cloudy_column(nlev, qv_g(i,j,:), o3col, &
                                rad%l_o3, dpc, rad%sw_toa(j), rad%coszen(j), rad%albedo, &
                                rad%albedo, cf, clwc, ciwc, heat_sw, sw_up, sw_dw, sw_net)
                        else
                            call aeros_lw_cloudy_column(nlev, t_g(i,j,:), qv_g(i,j,:), &
                                o3col, dpc, z_half, t_s(i,j), rad%q_co2, rad%l_o3, &
                                cf, clwc, ciwc, fnet, heat_lw, olr, fdw_lw)
                            call aeros_sw_cloudy_column(nlev, qv_g(i,j,:), o3col, &
                                rad%l_o3, dpc, rad%sw_toa(j), rad%coszen(j), rad%albedo, &
                                rad%albedo, cf, clwc, ciwc, heat_sw, sw_up, sw_dw, sw_net)
                        end if
                    else
                        ! Clear-sky longwave: dispatch on the scheme selector.
                        ! SESAM is the default (bit-for-bit unchanged); the opt-in
                        ! correlated-k kernels (LW and SW) are reached only via
                        ! SCHEME_ECCKD.
                        if (rad%scheme == SCHEME_ECCKD) then
                            call aeros_ecckd_lw_clearsky_column(nlev, t_g(i,j,:), &
                                qv_g(i,j,:), o3col, dpc, z_half, t_s(i,j), &
                                rad%q_co2, rad%l_o3, fnet, heat_lw, olr, fdw_lw)
                            call aeros_ecckd_sw_clearsky_column(nlev, qv_g(i,j,:), &
                                o3col, rad%l_o3, dpc, rad%sw_toa(j), rad%coszen(j), &
                                rad%albedo, rad%albedo, heat_sw, sw_up, sw_dw, sw_net)
                        else
                            call aeros_lw_clearsky_column(nlev, t_g(i,j,:), qv_g(i,j,:), &
                                o3col, dpc, z_half, t_s(i,j), rad%q_co2, rad%l_o3, &
                                fnet, heat_lw, olr, fdw_lw)
                            call aeros_sw_clearsky_column(nlev, qv_g(i,j,:), o3col, &
                                rad%l_o3, dpc, rad%sw_toa(j), rad%coszen(j), rad%albedo, &
                                rad%albedo, heat_sw, sw_up, sw_dw, sw_net)
                        end if
                    end if

                    rad%heat(i,j,:) = heat_lw + heat_sw

                    rad%olr(i,j)        = olr
                    rad%lw_dw_sur(i,j)  = fdw_lw
                    rad%sw_dw_sur(i,j)  = sw_dw
                    rad%sw_net_sur(i,j) = sw_net
                    rad%sw_up_toa(i,j)  = sw_up
                end do
            end do
            !$omp end parallel do
        end if

        ! --- accumulate the cached heating as a forward increment ----------
        do k = 1, nlev
            do j = 1, rad%nlat
                do i = 1, rad%nlon
                    dt_phys(i,j,k) = dt_phys(i,j,k) + rad%heat(i,j,k)*dt
                end do
            end do
        end do

        return
    end subroutine aeros_radiation_apply

    subroutine aeros_radiation_report(rad, io_unit)
        implicit none
        type(aeros_rad_class), intent(in) :: rad
        integer,               intent(in) :: io_unit

        write(io_unit, '(a)')        "  radiation:"
        write(io_unit, '(a,l1)')     "    enabled  = ", rad%enabled
        write(io_unit, '(a,i0)')     "    scheme   = ", rad%scheme
        write(io_unit, '(a,f8.2)')   "    co2_ppm  = ", rad%co2_ppm
        write(io_unit, '(a,l1)')     "    l_o3     = ", rad%l_o3
        write(io_unit, '(a,l1)')     "    clouds   = ", rad%clouds
        return
    end subroutine aeros_radiation_report

end module aeros_radiation
