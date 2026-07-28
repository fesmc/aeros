module aeros_ocean
    ! The lower boundary's temperature: a sea surface the atmosphere exchanges
    ! heat and moisture with. This module OWNS the SST field; aeros_surface reads
    ! it for the turbulent fluxes and aeros_radiation for the skin temperature.
    !
    ! === Two modes, one interface =========================================
    !
    ! OCEAN_PRESCRIBED -- SST is fixed (an infinite reservoir), the canonical
    !   aquaplanet control. This is the M2.4 behaviour, kept as the default so
    !   every existing run and test is unchanged.
    !
    ! OCEAN_SLAB -- SST is prognostic: a well-mixed slab of water of depth `depth`
    !   with heat capacity C = rho_w cp_w depth, integrated from the net surface
    !   energy flux,
    !
    !       C dSST/dt = SW_net_sfc + LW_down_sfc - sigma SST^4 - SH - LH .
    !
    !   The slab is what makes the surface flux SELF-LIMITING: a prescribed SST
    !   pumps heat into a column without bound, but a slab warms in response and
    !   the flux gradient closes, which is what a bounded radiative-convective
    !   equilibrium needs (docs/m2_results.md section 15; design.md section 6.1).
    !   The slab does NOT by itself fix a global TOA imbalance -- that is set by
    !   the albedo/greenhouse and is a separate knob.
    !
    ! This is the plug point for a real ocean: a dynamical ocean model replaces
    ! this module wholesale -- the atmosphere hands it the net surface flux and
    ! it returns the SST, the same contract the slab satisfies.
    !
    ! === Initial / prescribed SST =========================================
    !
    ! The APE "control" profile (Neale & Hoskins 2001), moved here from
    ! aeros_surface: T_s(phi) = T0 + sst_eq (1 - sin^2(1.5 phi)) for |phi| <
    ! sst_lat, T0 poleward. In prescribed mode it is the fixed SST; in slab mode
    ! it is the initial condition the slab evolves away from.
    !
    ! No sea ice: with a freeze floor (`freeze_floor`, default on) a slab SST is
    ! clamped at T0 rather than modelling ice, the honest minimal guard until a
    ! sea-ice model exists.
    !
    ! === Thermodynamic sea ice (opt-in, `l_seaice`) ========================
    !
    ! When `l_seaice` is on, the freeze-floor clamp is REPLACED by a Semtner
    ! (1976) 0-layer thermodynamic ice model (aeros_seaice_step_cell below). The
    ! mixed-layer heat deficit at the freezing point forms ice (latent heat of
    ! fusion) instead of being discarded by the clamp; ice grows/melts from the
    ! surface energy balance and the conductive flux through the slab, and the
    ! atmosphere sees the ICE-SURFACE skin temperature (carried in `sst`) rather
    ! than the mixed-layer temperature. Where ice is present the surface albedo
    ! rises to `ice_albedo` -- the feedback that closes the high-latitude balance,
    ! wired into radiation through the `alb` field. Off by default: `l_seaice =
    ! .FALSE.` is the exact freeze-floor behaviour, bit-for-bit.

    use aeros_defs,     only : dp, wp, io_unit_err, sigma_sb, T0, pi, L_f, &
                               aeros_grid_class
    use nml,            only : nml_read

    implicit none

    private

    integer, parameter, public :: OCEAN_PRESCRIBED = 0
    integer, parameter, public :: OCEAN_SLAB       = 1

    ! Water properties for the slab heat capacity (sea water).
    real(wp), parameter :: RHO_W = 1025.0_wp    ! density [kg m-3]
    real(wp), parameter :: CP_W  = 3990.0_wp    ! specific heat [J kg-1 K-1]

    ! Sea-ice properties (Semtner 0-layer). Named parameters; the tunable knobs
    ! (albedo, conductivity, freezing point) are namelist-configurable fields on
    ! the class below.
    real(wp), parameter :: RHO_ICE = 917.0_wp   ! sea-ice density [kg m-3]
    integer,  parameter :: ICE_NEWTON = 10      ! surface-balance Newton iterations

    type, public :: aeros_ocean_class
        integer :: mode = OCEAN_PRESCRIBED

        integer :: nlon = 0, nlat = 0

        ! Slab parameters.
        real(wp) :: depth        = 10.0_wp      ! mixed-layer depth [m]
        real(wp) :: heat_cap     = 0.0_wp       ! C = rho_w cp_w depth [J m-2 K-1]
        logical  :: freeze_floor = .TRUE.       ! clamp SST >= T0 (no sea ice)

        ! === Thermodynamic sea ice (opt-in) ================================
        ! When l_seaice is on the freeze floor is bypassed in favour of the
        ! Semtner 0-layer ice model. Params (namelist-configurable via the driver):
        logical  :: l_seaice   = .FALSE.        ! master switch (off = freeze_floor)
        real(wp) :: ice_albedo = 0.60_wp        ! surface albedo over ice [-]
        real(wp) :: ocn_albedo = 0.06_wp        ! surface albedo over open water [-]
        real(wp) :: k_ice      = 2.0_wp         ! ice thermal conductivity [W m-1 K-1]
        real(wp) :: t_frz      = 271.35_wp      ! seawater freezing point [K]

        ! Prescribed / initial SST profile (APE control).
        real(wp) :: sst_eq  = 27.0_wp           ! equator-minus-freezing SST [K]
        real(wp) :: sst_lat = 60.0_wp           ! frozen poleward of this [deg]

        ! Sea surface temperature [K], (nlon,nlat). Read by surface and radiation.
        ! With sea ice on this carries the ICE-SURFACE skin temperature where ice
        ! is present, and the mixed-layer temperature over open water.
        real(wp), allocatable :: sst(:,:)

        ! Net surface energy flux from the last step [W m-2], (nlon,nlat), a
        ! diagnostic (positive INTO the ocean). With sea ice on this is the true
        ! flux whose *dt equals the change in system energy (latent + sensible).
        real(wp), allocatable :: fnet(:,:)

        ! === Sea-ice prognostic state (allocated only when l_seaice) =======
        real(wp), allocatable :: h_ice(:,:)     ! ice thickness [m]
        real(wp), allocatable :: a_ice(:,:)     ! ice fraction [0-1] (0/1 in this cut)
        real(wp), allocatable :: t_ice_sfc(:,:) ! ice surface temperature [K]

        ! Surface broadband albedo field [-], (nlon,nlat): ocn_albedo over open
        ! water, ice_albedo over ice. Radiation reads this (via rad%alb_map) so
        ! the ice-albedo feedback is active. Allocated only when l_seaice.
        real(wp), allocatable :: alb(:,:)
    end type aeros_ocean_class

    public :: aeros_ocean_init
    public :: aeros_ocean_load
    public :: aeros_ocean_end
    public :: aeros_ocean_step
    public :: aeros_ocean_report
    public :: aeros_ocean_albedo_update

contains

    subroutine aeros_ocean_init(ocn, grd)
        ! Allocate and build the APE SST profile (the fixed SST, or the slab's
        ! initial condition).

        implicit none
        type(aeros_ocean_class), intent(inout) :: ocn
        type(aeros_grid_class),  intent(in)    :: grd

        integer  :: i, j
        real(wp) :: latr, dsst

        call aeros_ocean_end(ocn)

        ocn%nlon = grd%nlon
        ocn%nlat = grd%nlat
        ocn%heat_cap = RHO_W*CP_W*ocn%depth

        allocate(ocn%sst(grd%nlon, grd%nlat))
        allocate(ocn%fnet(grd%nlon, grd%nlat))

        do j = 1, grd%nlat
            latr = grd%lat(j)*pi/180.0_wp
            if (abs(grd%lat(j)) < ocn%sst_lat) then
                dsst = ocn%sst_eq*(1.0_wp - sin(1.5_wp*latr)**2)
            else
                dsst = 0.0_wp
            end if
            do i = 1, grd%nlon
                ocn%sst(i,j) = T0 + dsst
            end do
        end do

        ocn%fnet = 0.0_wp

        ! Sea-ice state: ice-free start (the APE polar SST is T0 > t_frz, so no
        ! ice is present initially), open-water albedo everywhere.
        if (ocn%l_seaice) then
            allocate(ocn%h_ice(grd%nlon, grd%nlat))
            allocate(ocn%a_ice(grd%nlon, grd%nlat))
            allocate(ocn%t_ice_sfc(grd%nlon, grd%nlat))
            allocate(ocn%alb(grd%nlon, grd%nlat))
            ocn%h_ice     = 0.0_wp
            ocn%a_ice     = 0.0_wp
            ocn%t_ice_sfc = ocn%t_frz
            ocn%alb       = ocn%ocn_albedo
        end if

        return
    end subroutine aeros_ocean_init

    subroutine aeros_ocean_load(ocn, filename, grd, defaults_file)
        ! Read the `ocean` namelist group, then init.

        implicit none
        type(aeros_ocean_class), intent(inout) :: ocn
        character(len=*),        intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_grid_class),  intent(in)    :: grd

        integer  :: mode
        real(wp) :: depth, sst_eq, sst_lat
        logical  :: freeze_floor
        logical  :: l_seaice
        real(wp) :: ice_albedo, ocn_albedo, k_ice, t_frz

        mode = ocn%mode; depth = ocn%depth; freeze_floor = ocn%freeze_floor
        sst_eq = ocn%sst_eq; sst_lat = ocn%sst_lat
        l_seaice = ocn%l_seaice; ice_albedo = ocn%ice_albedo
        ocn_albedo = ocn%ocn_albedo; k_ice = ocn%k_ice; t_frz = ocn%t_frz

        call nml_read(filename, "ocean", "mode",         mode, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "depth",        depth, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "freeze_floor", freeze_floor, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "sst_eq",       sst_eq, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "sst_lat",      sst_lat, defaults_file=defaults_file)
        ! Sea-ice knobs (opt-in). Off by default -> exact freeze-floor behaviour.
        call nml_read(filename, "ocean", "l_seaice",     l_seaice, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "ice_albedo",   ice_albedo, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "ocn_albedo",   ocn_albedo, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "k_ice",        k_ice, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "t_frz",        t_frz, defaults_file=defaults_file)

        ocn%mode = mode; ocn%depth = depth; ocn%freeze_floor = freeze_floor
        ocn%sst_eq = sst_eq; ocn%sst_lat = sst_lat
        ocn%l_seaice = l_seaice; ocn%ice_albedo = ice_albedo
        ocn%ocn_albedo = ocn_albedo; ocn%k_ice = k_ice; ocn%t_frz = t_frz

        call aeros_ocean_init(ocn, grd)

        return
    end subroutine aeros_ocean_load

    subroutine aeros_ocean_end(ocn)
        implicit none
        type(aeros_ocean_class), intent(inout) :: ocn
        if (allocated(ocn%sst))  deallocate(ocn%sst)
        if (allocated(ocn%fnet)) deallocate(ocn%fnet)
        if (allocated(ocn%h_ice))     deallocate(ocn%h_ice)
        if (allocated(ocn%a_ice))     deallocate(ocn%a_ice)
        if (allocated(ocn%t_ice_sfc)) deallocate(ocn%t_ice_sfc)
        if (allocated(ocn%alb))       deallocate(ocn%alb)
        ocn%nlon = 0; ocn%nlat = 0
        return
    end subroutine aeros_ocean_end

    subroutine aeros_ocean_step(ocn, sw_net_sfc, lw_dw_sfc, shf, lhf, dt)
        ! Advance the SST by one step from the net surface energy flux. A no-op in
        ! prescribed mode. Forward Euler: with a 10 m slab the per-step change is
        ! O(1e-3 K), far inside stability, and the slab's own response time (~C /
        ! 4 sigma T^3, months) is the physically relevant scale.
        !
        !   F_net = SW_net_sfc + LW_down_sfc - sigma SST^4 - SH - LH   (into ocean)

        implicit none
        type(aeros_ocean_class), intent(inout) :: ocn
        real(wp), intent(in) :: sw_net_sfc(:,:)   ! absorbed SW at surface [W m-2]
        real(wp), intent(in) :: lw_dw_sfc(:,:)    ! downward LW at surface [W m-2]
        real(wp), intent(in) :: shf(:,:)          ! sensible heat flux (up) [W m-2]
        real(wp), intent(in) :: lhf(:,:)          ! latent heat flux (up) [W m-2]
        real(wp), intent(in) :: dt                ! [s]

        integer  :: i, j
        real(wp) :: f

        if (ocn%mode /= OCEAN_SLAB) return

        if (ocn%l_seaice) then
            ! Thermodynamic sea ice replaces the freeze-floor clamp: the
            ! mixed-layer deficit forms ice, ice grows/melts from the surface
            ! energy balance, and `sst` carries the ice-surface skin temperature.
            !$omp parallel do collapse(2) schedule(static) private(i,j)
            do j = 1, ocn%nlat
                do i = 1, ocn%nlon
                    call aeros_seaice_step_cell(sw_net_sfc(i,j), lw_dw_sfc(i,j), &
                            shf(i,j), lhf(i,j), dt, ocn%heat_cap, ocn%k_ice, ocn%t_frz, &
                            ocn%sst(i,j), ocn%h_ice(i,j), ocn%a_ice(i,j), &
                            ocn%t_ice_sfc(i,j), ocn%fnet(i,j))
                end do
            end do
            !$omp end parallel do
            call aeros_ocean_albedo_update(ocn)
            return
        end if

        !$omp parallel do collapse(2) schedule(static) private(i,j,f)
        do j = 1, ocn%nlat
            do i = 1, ocn%nlon
                f = sw_net_sfc(i,j) + lw_dw_sfc(i,j) &
                    - sigma_sb*ocn%sst(i,j)**4 - shf(i,j) - lhf(i,j)
                ocn%fnet(i,j) = f
                ocn%sst(i,j)  = ocn%sst(i,j) + f*dt/ocn%heat_cap
                if (ocn%freeze_floor .and. ocn%sst(i,j) < T0) ocn%sst(i,j) = T0
            end do
        end do
        !$omp end parallel do

        return
    end subroutine aeros_ocean_step

    subroutine aeros_seaice_step_cell(sw_net, lw_dw, shf, lhf, dt, cheat, kice, tfrz, &
                                      sst, h, aice, tice, fnet)
        ! One column of the Semtner (1976) 0-layer thermodynamic ice + mixed-layer
        ! slab. Energy-conserving by construction: `fnet` is set to the flux whose
        ! *dt equals the change in system energy E = C(T_mix - tfrz) - rho_ice L_f h.
        !
        ! Three regimes, on the state at entry:
        !   OPEN WATER (h = 0): the ordinary slab update. If the step would cool
        !     the mixed layer below tfrz, the heat deficit freezes ice instead of
        !     being clamped away (the mixed layer is pinned at tfrz while ice
        !     exists) and `sst` switches to the ice surface (initially tfrz).
        !   ICE PRESENT (h > 0): solve the ice-surface energy balance
        !       sw_net + lw_dw - sigma T^4 - shf - lhf + kice (tfrz - T)/h = 0
        !     for the surface temperature T (Newton), capped at the melt point tfrz.
        !     The net atmospheric flux at that (capped) temperature grows ice at the
        !     base (deficit) or melts it from the top (surplus). `sst` = T.
        !   MELT-THROUGH: if the melt removes all the ice, the leftover melt energy
        !     warms the mixed layer above tfrz and the cell returns to open water.
        !
        ! sw_net/lw_dw/shf/lhf are the atmosphere's fluxes for this step, evaluated
        ! against the entry skin temperature; they are held constant across the
        ! surface-temperature solve (a 0-layer forward step).

        implicit none
        real(wp), intent(in)    :: sw_net, lw_dw, shf, lhf, dt
        real(wp), intent(in)    :: cheat        ! mixed-layer heat capacity C [J m-2 K-1]
        real(wp), intent(in)    :: kice, tfrz
        real(wp), intent(inout) :: sst, h, aice, tice
        real(wp), intent(out)   :: fnet

        real(wp) :: rlf, fatm, tnew, qdef, q, dh, eleft, g, gp, t
        integer  :: it

        rlf = RHO_ICE*real(L_f, wp)              ! volumetric latent heat [J m-3]

        if (h <= 0.0_wp) then
            ! --- open water: ordinary slab update, freeze-up on overshoot ------
            fatm = sw_net + lw_dw - sigma_sb*sst**4 - shf - lhf
            tnew = sst + fatm*dt/cheat
            fnet = fatm
            if (tnew >= tfrz) then
                sst  = tnew                      ! stays open water
                h    = 0.0_wp
                aice = 0.0_wp
                tice = tfrz
            else
                ! deficit below freezing forms ice; mixed layer pinned at tfrz
                qdef = cheat*(tfrz - tnew)       ! [J m-2] > 0
                h    = qdef/rlf
                aice = 1.0_wp
                tice = tfrz
                sst  = tfrz                       ! skin = new (thin) ice surface
            end if
            return
        end if

        ! --- ice present: solve the surface energy balance for T --------------
        t = tice
        do it = 1, ICE_NEWTON
            g  = sw_net + lw_dw - sigma_sb*t**4 - shf - lhf + kice*(tfrz - t)/h
            gp = -4.0_wp*sigma_sb*t**3 - kice/h
            t  = t - g/gp
            if (t < 150.0_wp) t = 150.0_wp
        end do
        if (t > tfrz) t = tfrz                   ! cap at the melt point
        tice = t

        ! Net atmospheric flux into the ice top at the (capped) surface
        ! temperature. Negative => deficit balanced by basal freezing (grow);
        ! positive => surplus at the melt cap melts the top. dh has the opposite
        ! sign so that dE = -rlf*dh = q*dt (exact energy conservation).
        q  = sw_net + lw_dw - sigma_sb*tice**4 - shf - lhf
        dh = -q*dt/rlf
        h  = h + dh
        fnet = q

        if (h <= 0.0_wp) then
            ! melted through: leftover melt energy warms the mixed layer
            eleft = (-h)*rlf                     ! [J m-2] >= 0
            h     = 0.0_wp
            aice  = 0.0_wp
            tice  = tfrz
            sst   = tfrz + eleft/cheat
        else
            aice = 1.0_wp
            sst  = tice
        end if

        return
    end subroutine aeros_seaice_step_cell

    subroutine aeros_ocean_albedo_update(ocn)
        ! Rebuild the surface albedo field from the ice state: open-water albedo
        ! blended to ice albedo by ice fraction. A no-op unless the ice fields are
        ! allocated (l_seaice). Radiation consumes this field for the ice-albedo
        ! feedback.

        implicit none
        type(aeros_ocean_class), intent(inout) :: ocn
        integer :: i, j

        if (.not. allocated(ocn%alb)) return

        do j = 1, ocn%nlat
            do i = 1, ocn%nlon
                ocn%alb(i,j) = (1.0_wp - ocn%a_ice(i,j))*ocn%ocn_albedo &
                             + ocn%a_ice(i,j)*ocn%ice_albedo
            end do
        end do

        return
    end subroutine aeros_ocean_albedo_update

    subroutine aeros_ocean_report(ocn, io_unit)
        implicit none
        type(aeros_ocean_class), intent(in) :: ocn
        integer,                 intent(in) :: io_unit
        write(io_unit, '(a)')    "  ocean:"
        if (ocn%mode == OCEAN_SLAB) then
            write(io_unit, '(a,f8.2,a)') "    mode    = slab, depth ", ocn%depth, " m"
            if (ocn%l_seaice) then
                write(io_unit, '(a,f5.2,a,f6.1,a)') "    sea ice = on, albedo ", &
                    ocn%ice_albedo, ", t_frz ", ocn%t_frz, " K"
                if (allocated(ocn%h_ice)) &
                    write(io_unit, '(a,f8.3,a,f8.3,a)') "    h_ice   = ", &
                        minval(ocn%h_ice), " to ", maxval(ocn%h_ice), " m"
            else
                write(io_unit, '(a,l1)') "    freeze_floor = ", ocn%freeze_floor
            end if
        else
            write(io_unit, '(a)')        "    mode    = prescribed SST"
        end if
        if (allocated(ocn%sst)) &
            write(io_unit, '(a,f8.2,a,f8.2,a)') "    SST     = ", &
                minval(ocn%sst), " to ", maxval(ocn%sst), " K"
        return
    end subroutine aeros_ocean_report

end module aeros_ocean
