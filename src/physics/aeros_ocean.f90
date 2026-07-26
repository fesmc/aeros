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

    use aeros_defs,     only : dp, wp, io_unit_err, sigma_sb, T0, pi, &
                               aeros_grid_class
    use nml,            only : nml_read

    implicit none

    private

    integer, parameter, public :: OCEAN_PRESCRIBED = 0
    integer, parameter, public :: OCEAN_SLAB       = 1

    ! Water properties for the slab heat capacity (sea water).
    real(wp), parameter :: RHO_W = 1025.0_wp    ! density [kg m-3]
    real(wp), parameter :: CP_W  = 3990.0_wp    ! specific heat [J kg-1 K-1]

    type, public :: aeros_ocean_class
        integer :: mode = OCEAN_PRESCRIBED

        integer :: nlon = 0, nlat = 0

        ! Slab parameters.
        real(wp) :: depth        = 10.0_wp      ! mixed-layer depth [m]
        real(wp) :: heat_cap     = 0.0_wp       ! C = rho_w cp_w depth [J m-2 K-1]
        logical  :: freeze_floor = .TRUE.       ! clamp SST >= T0 (no sea ice)

        ! Prescribed / initial SST profile (APE control).
        real(wp) :: sst_eq  = 27.0_wp           ! equator-minus-freezing SST [K]
        real(wp) :: sst_lat = 60.0_wp           ! frozen poleward of this [deg]

        ! Sea surface temperature [K], (nlon,nlat). Read by surface and radiation.
        real(wp), allocatable :: sst(:,:)

        ! Net surface energy flux from the last step [W m-2], (nlon,nlat), a
        ! diagnostic (positive INTO the ocean).
        real(wp), allocatable :: fnet(:,:)
    end type aeros_ocean_class

    public :: aeros_ocean_init
    public :: aeros_ocean_load
    public :: aeros_ocean_end
    public :: aeros_ocean_step
    public :: aeros_ocean_report

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

        mode = ocn%mode; depth = ocn%depth; freeze_floor = ocn%freeze_floor
        sst_eq = ocn%sst_eq; sst_lat = ocn%sst_lat

        call nml_read(filename, "ocean", "mode",         mode, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "depth",        depth, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "freeze_floor", freeze_floor, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "sst_eq",       sst_eq, defaults_file=defaults_file)
        call nml_read(filename, "ocean", "sst_lat",      sst_lat, defaults_file=defaults_file)

        ocn%mode = mode; ocn%depth = depth; ocn%freeze_floor = freeze_floor
        ocn%sst_eq = sst_eq; ocn%sst_lat = sst_lat

        call aeros_ocean_init(ocn, grd)

        return
    end subroutine aeros_ocean_load

    subroutine aeros_ocean_end(ocn)
        implicit none
        type(aeros_ocean_class), intent(inout) :: ocn
        if (allocated(ocn%sst))  deallocate(ocn%sst)
        if (allocated(ocn%fnet)) deallocate(ocn%fnet)
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

    subroutine aeros_ocean_report(ocn, io_unit)
        implicit none
        type(aeros_ocean_class), intent(in) :: ocn
        integer,                 intent(in) :: io_unit
        write(io_unit, '(a)')    "  ocean:"
        if (ocn%mode == OCEAN_SLAB) then
            write(io_unit, '(a,f8.2,a)') "    mode    = slab, depth ", ocn%depth, " m"
            write(io_unit, '(a,l1)')     "    freeze_floor = ", ocn%freeze_floor
        else
            write(io_unit, '(a)')        "    mode    = prescribed SST"
        end if
        if (allocated(ocn%sst)) &
            write(io_unit, '(a,f8.2,a,f8.2,a)') "    SST     = ", &
                minval(ocn%sst), " to ", maxval(ocn%sst), " K"
        return
    end subroutine aeros_ocean_report

end module aeros_ocean
