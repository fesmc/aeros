module aeros_land
    ! The land surface as the model's second lower-boundary type. Where
    ! aeros_ocean owns the sea surface (a prescribed or slab SST), this module
    ! owns the continents: a land-sea mask, a single-layer bucket soil moisture,
    ! and a slab soil temperature that IS the land skin temperature -- the land
    ! analogue of the slab ocean's SST.
    !
    ! === One flux seam, two surface types =================================
    !
    ! Land and ocean are two surface types under the SAME turbulent-flux seam
    ! (aeros_surface). The coupling is deliberately symmetric with the ocean:
    !
    !   * aeros_land_pre  composes, once per step, the skin temperature the
    !     surface fluxes and radiation are driven against -- the ocean SST,
    !     overwritten by the prognostic soil temperature on land points -- and
    !     the evaporation efficiency beta (1 over ocean, moisture-limited over
    !     land). aeros_surface then runs UNCHANGED against these fields; it never
    !     learns there is such a thing as land.
    !
    !   * aeros_land_step, called after radiation like aeros_ocean_step, advances
    !     the two prognostics from the fluxes just computed:
    !       c_soil dT_soil/dt = SW_net + LW_down - sigma T_soil^4 - SH - LH
    !       dw/dt             = (precip - evap)/rho_w,  capped at field capacity
    !     runoff is whatever would carry w above field capacity.
    !
    ! === Beta-limited evapotranspiration ==================================
    !
    ! Actual evaporation over land is beta * potential, beta = min(1, w/w_crit).
    ! A wet soil (w >= w_crit) evaporates at the potential rate; a dry soil
    ! throttles evaporation toward zero, so a desert stays a desert instead of
    ! being pinned at the ocean's saturation flux. This is the Manabe (1969)
    ! bucket, the honest minimal hydrology: one layer, one capacity, one
    ! efficiency curve. A multi-layer soil, vegetation and a real runoff routing
    ! all replace pieces here without touching the flux seam.
    !
    ! === Albedo and roughness =============================================
    !
    ! Land albedo is read per-point from ERA5 forecast albedo (fal) via the
    ! bcinput regridder; radiation reads it through rad%alb_map (composed by
    ! aeros_land_couple_radiation: the ocean scalar albedo everywhere, the land
    ! map on land points). Roughness enters as land exchange coefficients
    ! (c_h_land / c_e_land), a land constant feeding the same bulk formulae;
    ! ocean points keep aeros_surface's coefficients.
    !
    ! === Opt-in, bit-reproducible =========================================
    !
    ! Disabled by default (enabled = .FALSE.): nothing is allocated, every seam
    ! call returns at once, and the model is the all-ocean aquaplanet bit for
    ! bit. The two prognostics (w, t_soil) are new state and are serialized by
    ! aeros_timestep_write_restart / _read_restart; the mask and albedo maps are
    ! static and rebuilt from file on init, so they are not saved.

    use aeros_defs,    only : dp, wp, io_unit_err, sigma_sb, T0, grav, pi, &
                              aeros_grid_class
    use aeros_bcinput, only : aeros_bcinput_read_field
    use nml,           only : nml_read

    implicit none

    private

    ! Fresh water density for the bucket (soil water, not sea water).
    real(wp), parameter :: RHO_W = 1000.0_wp        ! [kg m-3]

    type, public :: aeros_land_class
        logical :: enabled = .FALSE.

        integer :: nlon = 0, nlat = 0

        ! Source files (ERA5 climatologies, regridded via aeros_bcinput).
        character(len=512) :: lsm_file    = ""      ! land-sea mask (var lsm)
        character(len=512) :: albedo_file = ""      ! forecast albedo (var fal)

        ! Parameters.
        real(wp) :: lsm_threshold    = 0.5_wp       ! land where lsm >= threshold
        real(wp) :: w_field_capacity = 0.15_wp      ! bucket capacity [m]
        real(wp) :: w_crit           = -1.0_wp      ! beta knee [m]; <0 => 0.75*fc
        real(wp) :: c_soil           = 2.0e6_wp     ! soil heat capacity [J m-2 K-1]
        real(wp) :: land_albedo      = 0.20_wp      ! fallback land albedo [-]
        real(wp) :: c_h_land         = 1.5e-3_wp    ! land sensible exchange coeff
        real(wp) :: c_e_land         = 1.5e-3_wp    ! land moisture exchange coeff
        logical  :: freeze_floor     = .FALSE.      ! clamp t_soil >= T0 (no snow)

        ! Initial conditions for the prognostics on a cold start (overwritten on
        ! restart). Kept off the namelist: they equilibrate quickly under the
        ! small soil heat capacity and are not a tuning knob.
        real(wp) :: t_soil_init = 288.0_wp          ! initial soil temperature [K]

        ! Static maps (rebuilt from file on init).
        logical,  allocatable :: mask(:,:)          ! .TRUE. on land
        real(wp), allocatable :: albedo(:,:)        ! land albedo (used on land)

        ! Prognostic state (serialized in the restart).
        real(wp), allocatable :: w(:,:)             ! soil moisture [m]
        real(wp), allocatable :: t_soil(:,:)        ! soil skin temperature [K]

        ! Per-step scratch, composed by aeros_land_pre.
        real(wp), allocatable :: skin(:,:)          ! skin T for surface+radiation [K]
        real(wp), allocatable :: beta(:,:)          ! evaporation efficiency [-]

        ! Diagnostics from the last step.
        real(wp), allocatable :: runoff(:,:)        ! runoff [kg m-2 s-1]
        real(wp), allocatable :: fnet(:,:)          ! net flux into soil [W m-2]
    end type aeros_land_class

    public :: aeros_land_init
    public :: aeros_land_load
    public :: aeros_land_end
    public :: aeros_land_pre
    public :: aeros_land_step
    public :: aeros_land_couple_radiation
    public :: aeros_land_report
    public :: aeros_land_is_land

contains

    elemental logical function aeros_land_is_land(lsm, threshold) result(is_land)
        ! The mask predicate, isolated so the threshold convention lives in one
        ! place and can be pinned by the acceptance test with no file. A cell is
        ! land where the (0-1) land fraction reaches the threshold.
        implicit none
        real(wp), intent(in) :: lsm, threshold
        is_land = (lsm >= threshold)
        return
    end function aeros_land_is_land

    subroutine aeros_land_init(land, grd, enabled)
        ! Allocate the state and, when enabled, read the mask and albedo maps and
        ! set the prognostic initial conditions. Disabled leaves everything
        ! deallocated -- the all-ocean, bit-for-bit default.

        implicit none
        type(aeros_land_class), intent(inout) :: land
        type(aeros_grid_class), intent(in)    :: grd
        logical,                intent(in)    :: enabled

        real(wp), allocatable :: tmp(:,:)
        integer :: i, j

        call aeros_land_end(land)

        land%enabled = enabled
        land%nlon    = grd%nlon
        land%nlat    = grd%nlat

        if (.not. enabled) return

        if (len_trim(land%lsm_file) == 0) then
            write(io_unit_err,*) "aeros_land_init:: error: l_land is on but lsm_file is empty."
            error stop 1
        end if

        ! Resolve the beta knee sentinel: default to 3/4 of field capacity.
        if (land%w_crit < 0.0_wp) land%w_crit = 0.75_wp*land%w_field_capacity

        allocate(land%mask(grd%nlon, grd%nlat))
        allocate(land%albedo(grd%nlon, grd%nlat))
        allocate(land%w(grd%nlon, grd%nlat))
        allocate(land%t_soil(grd%nlon, grd%nlat))
        allocate(land%skin(grd%nlon, grd%nlat))
        allocate(land%beta(grd%nlon, grd%nlat))
        allocate(land%runoff(grd%nlon, grd%nlat))
        allocate(land%fnet(grd%nlon, grd%nlat))
        allocate(tmp(grd%nlon, grd%nlat))

        ! --- land-sea mask: regrid the (0-1) fraction, then threshold ---------
        call aeros_bcinput_read_field(trim(land%lsm_file), "lsm", &
                                      grd%lon, grd%lat, tmp)
        land%mask = aeros_land_is_land(tmp, land%lsm_threshold)

        ! --- land albedo: the ERA5 forecast albedo map, or a constant ---------
        if (len_trim(land%albedo_file) > 0) then
            call aeros_bcinput_read_field(trim(land%albedo_file), "fal", &
                                          grd%lon, grd%lat, land%albedo)
        else
            land%albedo = land%land_albedo
        end if

        ! --- prognostic initial conditions (cold start) -----------------------
        ! Soil starts at a uniform temperature and a half-full bucket; both
        ! equilibrate quickly. Ocean points carry harmless values (the seam
        ! overwrites them with the SST) but are set for cleanliness.
        land%w      = 0.5_wp*land%w_field_capacity
        land%t_soil = land%t_soil_init
        land%skin   = land%t_soil_init
        land%beta   = 1.0_wp
        land%runoff = 0.0_wp
        land%fnet   = 0.0_wp

        deallocate(tmp)

        return
    end subroutine aeros_land_init

    subroutine aeros_land_load(land, filename, grd, defaults_file)
        ! Read the `land` namelist group, then init. Mirrors aeros_ocean_load.

        implicit none
        type(aeros_land_class), intent(inout) :: land
        character(len=*),       intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_grid_class), intent(in)    :: grd

        logical            :: enabled, freeze_floor
        character(len=512) :: lsm_file, albedo_file
        real(wp)           :: lsm_threshold, w_field_capacity, w_crit, c_soil
        real(wp)           :: land_albedo, c_h_land, c_e_land

        enabled          = land%enabled
        lsm_file         = land%lsm_file
        albedo_file      = land%albedo_file
        lsm_threshold    = land%lsm_threshold
        w_field_capacity = land%w_field_capacity
        w_crit           = land%w_crit
        c_soil           = land%c_soil
        land_albedo      = land%land_albedo
        c_h_land         = land%c_h_land
        c_e_land         = land%c_e_land
        freeze_floor     = land%freeze_floor

        call nml_read(filename, "land", "enabled",          enabled, defaults_file=defaults_file)
        call nml_read(filename, "land", "lsm_file",         lsm_file, defaults_file=defaults_file)
        call nml_read(filename, "land", "albedo_file",      albedo_file, defaults_file=defaults_file)
        call nml_read(filename, "land", "lsm_threshold",    lsm_threshold, defaults_file=defaults_file)
        call nml_read(filename, "land", "w_field_capacity", w_field_capacity, defaults_file=defaults_file)
        call nml_read(filename, "land", "w_crit",           w_crit, defaults_file=defaults_file)
        call nml_read(filename, "land", "c_soil",           c_soil, defaults_file=defaults_file)
        call nml_read(filename, "land", "land_albedo",      land_albedo, defaults_file=defaults_file)
        call nml_read(filename, "land", "c_h_land",         c_h_land, defaults_file=defaults_file)
        call nml_read(filename, "land", "c_e_land",         c_e_land, defaults_file=defaults_file)
        call nml_read(filename, "land", "freeze_floor",     freeze_floor, defaults_file=defaults_file)

        land%lsm_file         = lsm_file
        land%albedo_file      = albedo_file
        land%lsm_threshold    = lsm_threshold
        land%w_field_capacity = w_field_capacity
        land%w_crit           = w_crit
        land%c_soil           = c_soil
        land%land_albedo      = land_albedo
        land%c_h_land         = c_h_land
        land%c_e_land         = c_e_land
        land%freeze_floor     = freeze_floor

        call aeros_land_init(land, grd, enabled)

        return
    end subroutine aeros_land_load

    subroutine aeros_land_end(land)
        implicit none
        type(aeros_land_class), intent(inout) :: land
        if (allocated(land%mask))   deallocate(land%mask)
        if (allocated(land%albedo)) deallocate(land%albedo)
        if (allocated(land%w))      deallocate(land%w)
        if (allocated(land%t_soil)) deallocate(land%t_soil)
        if (allocated(land%skin))   deallocate(land%skin)
        if (allocated(land%beta))   deallocate(land%beta)
        if (allocated(land%runoff)) deallocate(land%runoff)
        if (allocated(land%fnet))   deallocate(land%fnet)
        land%enabled = .FALSE.
        land%nlon = 0; land%nlat = 0
        return
    end subroutine aeros_land_end

    subroutine aeros_land_pre(land, sst)
        ! Compose, once per step, the skin temperature and evaporation efficiency
        ! the surface fluxes and radiation see: the ocean SST and beta = 1 over
        ! ocean, overwritten by the prognostic soil temperature and the
        ! moisture-limited beta on land. A no-op (and no allocation touched) when
        ! disabled, so the ocean path is bit-for-bit.

        implicit none
        type(aeros_land_class), intent(inout) :: land
        real(wp), intent(in) :: sst(:,:)         ! (nlon,nlat) ocean SST [K]

        integer :: i, j

        if (.not. land%enabled) return

        !$omp parallel do collapse(2) schedule(static) private(i,j)
        do j = 1, land%nlat
            do i = 1, land%nlon
                if (land%mask(i,j)) then
                    land%skin(i,j) = land%t_soil(i,j)
                    land%beta(i,j) = min(1.0_wp, land%w(i,j)/land%w_crit)
                else
                    land%skin(i,j) = sst(i,j)
                    land%beta(i,j) = 1.0_wp
                end if
            end do
        end do
        !$omp end parallel do

        return
    end subroutine aeros_land_pre

    subroutine aeros_land_step(land, sw_net_sfc, lw_dw_sfc, shf, lhf, evap, precip, dt)
        ! Advance the two land prognostics by one step from the fluxes just
        ! computed, over land points only. The soil analogue of aeros_ocean_step.
        !
        !   soil temperature (slab):
        !     c_soil dT_soil/dt = SW_net + LW_down - sigma T_soil^4 - SH - LH
        !   bucket soil moisture:
        !     dw = (precip - evap) dt / rho_w,  excess above field capacity runs off
        !
        ! A no-op when disabled, so a land-off run never touches this.

        implicit none
        type(aeros_land_class), intent(inout) :: land
        real(wp), intent(in) :: sw_net_sfc(:,:)   ! absorbed SW at surface [W m-2]
        real(wp), intent(in) :: lw_dw_sfc(:,:)    ! downward LW at surface [W m-2]
        real(wp), intent(in) :: shf(:,:)          ! sensible heat flux (up) [W m-2]
        real(wp), intent(in) :: lhf(:,:)          ! latent heat flux (up) [W m-2]
        real(wp), intent(in) :: evap(:,:)         ! evaporation [kg m-2 s-1]
        real(wp), intent(in) :: precip(:,:)       ! precipitation [kg m-2 s-1]
        real(wp), intent(in) :: dt                ! [s]

        integer  :: i, j
        real(wp) :: f, dw, wcap, excess

        if (.not. land%enabled) return

        wcap = land%w_field_capacity

        !$omp parallel do collapse(2) schedule(static) private(i,j,f,dw,excess)
        do j = 1, land%nlat
            do i = 1, land%nlon
                if (.not. land%mask(i,j)) cycle

                ! --- slab soil temperature ---------------------------------
                f = sw_net_sfc(i,j) + lw_dw_sfc(i,j) &
                    - sigma_sb*land%t_soil(i,j)**4 - shf(i,j) - lhf(i,j)
                land%fnet(i,j)   = f
                land%t_soil(i,j) = land%t_soil(i,j) + f*dt/land%c_soil
                if (land%freeze_floor .and. land%t_soil(i,j) < T0) &
                    land%t_soil(i,j) = T0

                ! --- bucket soil moisture ----------------------------------
                ! Source: precipitation reaching the ground. Sink: the
                ! evaporation the surface fluxes already removed to the
                ! atmosphere (beta-limited, so consistent with the LH used
                ! above). Water depth in metres; rho_w converts the mass flux.
                dw = (precip(i,j) - evap(i,j))*dt/RHO_W
                land%w(i,j) = land%w(i,j) + dw

                land%runoff(i,j) = 0.0_wp
                if (land%w(i,j) > wcap) then
                    excess = land%w(i,j) - wcap
                    land%w(i,j)      = wcap
                    land%runoff(i,j) = excess*RHO_W/dt      ! [kg m-2 s-1]
                else if (land%w(i,j) < 0.0_wp) then
                    ! Evaporation drew the bucket dry (beta should make this
                    ! rare); clamp and let the residual be an implicit source.
                    land%w(i,j) = 0.0_wp
                end if
            end do
        end do
        !$omp end parallel do

        return
    end subroutine aeros_land_step

    subroutine aeros_land_couple_radiation(land, alb_map, ocean_albedo)
        ! Build the per-point surface albedo map radiation reads: the ocean
        ! scalar albedo everywhere, overwritten by the land albedo on land
        ! points. Allocated on first call; a no-op when disabled, so radiation
        ! keeps its scalar albedo and the ocean run is bit-for-bit.

        implicit none
        type(aeros_land_class),            intent(in)    :: land
        real(wp), allocatable,             intent(inout) :: alb_map(:,:)
        real(wp),                          intent(in)    :: ocean_albedo

        integer :: i, j

        if (.not. land%enabled) return

        if (.not. allocated(alb_map)) allocate(alb_map(land%nlon, land%nlat))

        do j = 1, land%nlat
            do i = 1, land%nlon
                if (land%mask(i,j)) then
                    alb_map(i,j) = land%albedo(i,j)
                else
                    alb_map(i,j) = ocean_albedo
                end if
            end do
        end do

        return
    end subroutine aeros_land_couple_radiation

    subroutine aeros_land_report(land, io_unit)
        implicit none
        type(aeros_land_class), intent(in) :: land
        integer,                intent(in) :: io_unit
        integer  :: nland
        real(wp) :: frac

        write(io_unit, '(a)')    "  land:"
        write(io_unit, '(a,l1)') "    enabled = ", land%enabled
        if (.not. land%enabled) return
        nland = count(land%mask)
        frac  = real(nland, wp)/real(land%nlon*land%nlat, wp)
        write(io_unit, '(a,f6.3,a,i0,a)') "    land fraction = ", frac, &
            " (", nland, " cells)"
        write(io_unit, '(a,f6.3,a,f6.3,a)') "    w_field_capacity/w_crit = ", &
            land%w_field_capacity, " / ", land%w_crit, " m"
        write(io_unit, '(a,es9.2,a)') "    c_soil = ", land%c_soil, " J m-2 K-1"
        if (allocated(land%t_soil)) &
            write(io_unit, '(a,f7.2,a,f7.2,a)') "    t_soil = ", &
                minval(land%t_soil, mask=land%mask), " to ", &
                maxval(land%t_soil, mask=land%mask), " K"
        return
    end subroutine aeros_land_report

end module aeros_land
