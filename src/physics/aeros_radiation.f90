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
    ! Clear-sky longwave only. lwr_clouds (SESAM's cloudy branch) waits until
    ! clouds exist as a field; the clear-sky path is what validates against
    ! ERA5 clear-sky fluxes (ttrc, strdc). No grid apply and no timestep wiring
    ! yet: aeros_lw_clearsky_column is exposed and exercised offline by
    ! test_radiation on a single column. `enabled` defaults .FALSE. as for the
    ! other physics.

    use aeros_defs,     only : dp, wp, io_unit_err, R_d, cp_d, grav, T0, p0, &
                               sigma_sb, S0, aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure, &
                               aeros_hydrostatic
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

    ! water vapour: d_vap = 1/(1 + a_vap x^b + a2 x^b2 + a3 x^3),  x = beta0*u
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

    ! molar masses for the CO2 volume->mass mixing ratio conversion
    real(wp), parameter :: M_CO2 = 44.0095_wp      ! g/mol
    real(wp), parameter :: M_AIR = 28.97_wp        ! g/mol

    ! === Configuration / state =============================================
    type, public :: aeros_rad_class
        logical :: enabled = .FALSE.        ! off by default, like convection
        integer :: scheme  = SCHEME_SESAM

        integer :: nlon = 0
        integer :: nlat = 0

        real(wp) :: co2_ppm  = 280.0_wp     ! CO2 volume mixing ratio [ppmv]
        logical  :: l_o3     = .FALSE.       ! include ozone absorption
        real(wp) :: q_co2    = 0.0_wp        ! CO2 mass mixing ratio [kg kg-1], derived
    end type aeros_rad_class

    public :: aeros_radiation_init
    public :: aeros_radiation_load
    public :: aeros_radiation_end
    public :: aeros_radiation_report
    public :: aeros_lw_clearsky_column

contains

    subroutine aeros_radiation_init(rad, grd, enabled)
        ! Minimal init: geometry, defaults, and the derived CO2 mass ratio.

        implicit none
        type(aeros_rad_class),   intent(inout) :: rad
        type(aeros_grid_class),  intent(in)    :: grd
        logical,                 intent(in)    :: enabled

        rad%enabled = enabled
        rad%nlon    = grd%nlon
        rad%nlat    = grd%nlat
        rad%q_co2   = co2_mass_ratio(rad%co2_ppm)

        return
    end subroutine aeros_radiation_init

    subroutine aeros_radiation_load(rad, filename, grd)
        ! Read the `radiation` namelist group, then finish init.

        implicit none
        type(aeros_rad_class),   intent(inout) :: rad
        character(len=*),        intent(in)    :: filename
        type(aeros_grid_class),  intent(in)    :: grd

        logical  :: enabled
        integer  :: scheme
        real(wp) :: co2_ppm
        logical  :: l_o3

        enabled = rad%enabled
        scheme  = rad%scheme
        co2_ppm = rad%co2_ppm
        l_o3    = rad%l_o3

        call nml_read(filename, "radiation", "enabled", enabled)
        call nml_read(filename, "radiation", "scheme",  scheme)
        call nml_read(filename, "radiation", "co2_ppm", co2_ppm)
        call nml_read(filename, "radiation", "l_o3",    l_o3)

        rad%scheme  = scheme
        rad%co2_ppm = co2_ppm
        rad%l_o3    = l_o3

        call aeros_radiation_init(rad, grd, enabled)

        return
    end subroutine aeros_radiation_load

    subroutine aeros_radiation_end(rad)
        implicit none
        type(aeros_rad_class), intent(inout) :: rad
        rad%enabled = .FALSE.
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
            dv = 1.0_wp/(1.0_wp + LW_A_VAP *(LW_BETA0*aw)**LW_BETA_VAP  &
                                + LW_A2_VAP*(LW_BETA0*aw)**LW_BETA2_VAP &
                                + LW_A3_VAP*(LW_BETA0*aw)**3)
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

    subroutine aeros_radiation_report(rad, io_unit)
        implicit none
        type(aeros_rad_class), intent(in) :: rad
        integer,               intent(in) :: io_unit

        write(io_unit, '(a)')        "  radiation:"
        write(io_unit, '(a,l1)')     "    enabled  = ", rad%enabled
        write(io_unit, '(a,i0)')     "    scheme   = ", rad%scheme
        write(io_unit, '(a,f8.2)')   "    co2_ppm  = ", rad%co2_ppm
        write(io_unit, '(a,l1)')     "    l_o3     = ", rad%l_o3
        return
    end subroutine aeros_radiation_report

end module aeros_radiation
