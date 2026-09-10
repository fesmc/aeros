module aeros_boundary
    ! Boundary-provider abstraction for the library API refactor.
    !
    ! Every lower-boundary or forcing field that today lives in a physics
    ! module or a driver as a private (nlon,nlat) array (rce_long's sst_clim,
    ! aeros_land's albedo_clim, aeros_radiation's dust_clim / so4_map,
    ! aeros_topography's phis_full, ...) gets ONE provider here. A provider
    ! has a discriminated `kind`:
    !
    !   BC_CONSTANT  -- a scalar `value` (a uniform fill for 2D buffers).
    !   BC_FILE      -- read from `file` (via aeros_bcinput_read_field): a
    !                   single-slice `itime` or a 12-month climatology
    !                   (`is_monthly = .TRUE.`). Monthly fields are held in
    !                   `clim(:,:,1..12)`; the daily working `buffer(:,:)`
    !                   is refreshed by day-of-year interpolation.
    !   BC_EXTERNAL  -- the buffer is written by the caller (a coupled host
    !                   through aeros_lib_step(imports=...)); no I/O here.
    !
    ! `aeros_boundary_t` bundles one provider per known boundary field
    ! (spec docs/coupling_api.md sec 7); each is configured at init from the
    ! `&boundary` namelist group.
    !
    ! === Wave-2 scope (this commit) ========================================
    !
    ! This module is ADDITIVE: the types, the loader, init and refresh are
    ! defined and the boundary bundle is a new component on aeros_class, but
    ! no existing physics module or driver consumes from it yet. That wiring
    ! (routing aeros_land's albedo_clim, aeros_radiation's dust_clim /
    ! so4_map, aeros_topography's phis_full and rce_long's sst_clim through
    ! ams%bounds%<field>%buffer) is Wave 3, once aeros_lib_step is the entry
    ! point that calls aeros_boundary_refresh at the start of each cycle.
    ! Keeping the two waves separate preserves bit-for-bit trivially: the
    ! existing code paths are untouched, and every aeros_defaults.nml default
    ! written below matches today's paths and monthly conventions so a Wave-3
    ! rewire also lands bit-for-bit.
    !
    ! Bit-for-bit note: the daily interpolation kernel (mid-month centres,
    ! cyclic, x = mod(day, 365)/365*12 - 0.5) matches aeros_land_set_albedo
    ! (aeros_land.f90:430-435), aeros_radiation_set_dust (aeros_radiation.f90:
    ! 496-504) and rce_long's set_seasonal_sst -- one shared kernel here so
    ! Wave 3 is a one-line consumer rather than a re-derivation.

    use aeros_defs,    only : dp, wp, io_unit_err, aeros_grid_class
    use aeros_bcinput, only : aeros_bcinput_read_field
    use nml,           only : nml_read

    implicit none

    private

    integer, parameter, public :: BC_CONSTANT = 0
    integer, parameter, public :: BC_FILE     = 1
    integer, parameter, public :: BC_EXTERNAL = 2

    ! A scalar (year-frequency) boundary: CO2, CH4, N2O, TSI, time_bp.
    ! `year` selects the record in FILE mode; `current` holds the value the
    ! physics reads. FILE-mode readers for these are Phase C -- Wave 2 ships
    ! the type only.
    type, public :: aeros_bc_scalar_t
        integer            :: kind    = BC_CONSTANT
        real(wp)           :: value   = 0.0_wp     ! CONSTANT
        character(len=256) :: file    = ""         ! FILE (year-slice)
        integer            :: year    = 2000       ! FILE: year to select
        real(wp)           :: current = 0.0_wp     ! last-read value (physics-visible)
    end type aeros_bc_scalar_t

    ! A 2D boundary: (nlon,nlat) grid-consistent buffer. `buffer` is ALWAYS
    ! allocated after init (a physics reader never has to test allocated).
    ! `clim` is allocated only for monthly FILE fields.
    type, public :: aeros_bc_2d_t
        integer               :: kind       = BC_CONSTANT
        real(wp), allocatable :: buffer(:,:)           ! (nlon,nlat), physics reads
        character(len=256)    :: file       = ""       ! FILE
        character(len=64)     :: varname    = ""       ! FILE variable name
        character(len=64)     :: lon_name   = "longitude"
        character(len=64)     :: lat_name   = "latitude"
        logical               :: is_monthly = .FALSE.  ! FILE: 12-month climatology
        real(wp), allocatable :: clim(:,:,:)           ! (nlon,nlat,12) FILE monthly
        integer               :: itime      = 1        ! FILE single-slice: 1-based
        real(wp)              :: value      = 0.0_wp   ! CONSTANT (uniform fill)
        real(wp)              :: scale      = 1.0_wp   ! post-read multiplier (dust_scale et al.)
        real(wp)              :: floor      = -huge(1.0_wp) ! post-read floor (max(0, ...) et al.)
    end type aeros_bc_2d_t

    ! One provider per known boundary. Field names and units mirror
    ! docs/coupling_api.md sec 6.1 / 7 and CLIMBER-X's cmn_to_atm bundle.
    type, public :: aeros_boundary_t
        integer :: nlon = 0
        integer :: nlat = 0

        ! === 2D fields ====================================================
        type(aeros_bc_2d_t) :: sst           ! sea surface temperature [K]
        type(aeros_bc_2d_t) :: sic           ! sea-ice fraction [0-1]
        type(aeros_bc_2d_t) :: t_ice_sfc     ! sea-ice surface T [K]
        type(aeros_bc_2d_t) :: t_skin_land   ! land skin T [K]
        type(aeros_bc_2d_t) :: q_skin_land   ! land skin q [kg/kg]
        type(aeros_bc_2d_t) :: albedo        ! broadband albedo [-]
        type(aeros_bc_2d_t) :: z0m           ! momentum roughness [m]
        type(aeros_bc_2d_t) :: phis          ! surface geopotential [m2 s-2]
        type(aeros_bc_2d_t) :: lsm           ! land-sea mask (0..1)
        type(aeros_bc_2d_t) :: dust          ! dust optical thickness
        type(aeros_bc_2d_t) :: so4           ! sulfate load [kg m-2]
        type(aeros_bc_2d_t) :: o3            ! ozone column (2D placeholder; 3D is Phase C)

        ! Surface-geopotential spin-up ramp [days] (rce_long's topo_ramp_days:
        ! phis is scaled linearly from 0 to full over this many days). Wave 4a
        ! reads the key only; the ramp consumer moves into the library with
        ! the phis provider in Wave 4b. 0 = no ramp.
        real(wp) :: phis_ramp_days = 0.0_wp

        ! === Scalar (year-boundary) fields ================================
        type(aeros_bc_scalar_t) :: co2       ! CO2 [ppmv]
        type(aeros_bc_scalar_t) :: ch4       ! CH4 [ppbv] (Phase C reader)
        type(aeros_bc_scalar_t) :: n2o       ! N2O [ppbv] (Phase C reader)
        type(aeros_bc_scalar_t) :: time_bp   ! orbital year before present [yr]
        type(aeros_bc_scalar_t) :: tsi       ! total solar irradiance [W m-2]
    end type aeros_boundary_t

    public :: aeros_boundary_load
    public :: aeros_boundary_init
    public :: aeros_boundary_refresh
    public :: aeros_boundary_end
    public :: aeros_boundary_interp_monthly

contains

    subroutine aeros_boundary_load(bounds, filename, grd, defaults_file)
        ! Read the `&boundary` namelist group and populate each provider's
        ! kind / file / varname / value. Then init (allocate buffers and read
        ! BC_FILE fields).
        !
        ! Namelist convention: one key per provider field. For a 2D provider
        ! `X` the keys are X_kind, X_file, X_varname, X_is_monthly, X_itime,
        ! X_value, X_scale, X_floor, X_lon_name, X_lat_name; for a scalar,
        ! X_kind, X_value, X_file, X_year. Missing keys leave the compile-time
        ! defaults (BC_CONSTANT with value = 0), which is bit-for-bit "no
        ! provider" -- physics that hasn't been rewired to consume from
        ! ams%bounds ignores the empty buffers.

        implicit none

        type(aeros_boundary_t), intent(inout) :: bounds
        character(len=*),       intent(in)    :: filename
        type(aeros_grid_class), intent(in)    :: grd
        character(len=*),       intent(in), optional :: defaults_file

        call read_2d(bounds%sst,         "sst",         filename, defaults_file)
        call read_2d(bounds%sic,         "sic",         filename, defaults_file)
        call read_2d(bounds%t_ice_sfc,   "t_ice_sfc",   filename, defaults_file)
        call read_2d(bounds%t_skin_land, "t_skin_land", filename, defaults_file)
        call read_2d(bounds%q_skin_land, "q_skin_land", filename, defaults_file)
        call read_2d(bounds%albedo,      "albedo",      filename, defaults_file)
        call read_2d(bounds%z0m,         "z0m",         filename, defaults_file)
        call read_2d(bounds%phis,        "phis",        filename, defaults_file)
        call nml_read(filename, "boundary", "phis_ramp_days", bounds%phis_ramp_days, &
                      defaults_file=defaults_file)
        call read_2d(bounds%lsm,         "lsm",         filename, defaults_file)
        call read_2d(bounds%dust,        "dust",        filename, defaults_file)
        call read_2d(bounds%so4,         "so4",         filename, defaults_file)
        call read_2d(bounds%o3,          "o3",          filename, defaults_file)

        call read_scalar(bounds%co2,     "co2",     filename, defaults_file)
        call read_scalar(bounds%ch4,     "ch4",     filename, defaults_file)
        call read_scalar(bounds%n2o,     "n2o",     filename, defaults_file)
        call read_scalar(bounds%time_bp, "time_bp", filename, defaults_file)
        call read_scalar(bounds%tsi,     "tsi",     filename, defaults_file)

        call aeros_boundary_init(bounds, grd)

        return
    end subroutine aeros_boundary_load

    subroutine read_2d(fld, name, filename, defaults_file)
        implicit none
        type(aeros_bc_2d_t), intent(inout) :: fld
        character(len=*),    intent(in)    :: name, filename
        character(len=*),    intent(in), optional :: defaults_file

        integer            :: kind_, itime
        character(len=256) :: file
        character(len=64)  :: varname, lon_name, lat_name
        logical            :: is_monthly
        real(wp)           :: value, scale, floor_

        kind_      = fld%kind
        file       = fld%file
        varname    = fld%varname
        lon_name   = fld%lon_name
        lat_name   = fld%lat_name
        is_monthly = fld%is_monthly
        itime      = fld%itime
        value      = fld%value
        scale      = fld%scale
        floor_     = fld%floor

        call nml_read(filename, "boundary", trim(name)//"_kind",       kind_,      defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_file",       file,       defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_varname",    varname,    defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_lon_name",   lon_name,   defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_lat_name",   lat_name,   defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_is_monthly", is_monthly, defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_itime",      itime,      defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_value",      value,      defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_scale",      scale,      defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_floor",      floor_,     defaults_file=defaults_file)

        fld%kind       = kind_
        fld%file       = file
        fld%varname    = varname
        fld%lon_name   = lon_name
        fld%lat_name   = lat_name
        fld%is_monthly = is_monthly
        fld%itime      = itime
        fld%value      = value
        fld%scale      = scale
        fld%floor      = floor_

        return
    end subroutine read_2d

    subroutine read_scalar(fld, name, filename, defaults_file)
        implicit none
        type(aeros_bc_scalar_t), intent(inout) :: fld
        character(len=*),        intent(in)    :: name, filename
        character(len=*),        intent(in), optional :: defaults_file

        integer            :: kind_, year
        real(wp)           :: value
        character(len=256) :: file

        kind_ = fld%kind
        value = fld%value
        file  = fld%file
        year  = fld%year

        call nml_read(filename, "boundary", trim(name)//"_kind",  kind_, defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_value", value, defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_file",  file,  defaults_file=defaults_file)
        call nml_read(filename, "boundary", trim(name)//"_year",  year,  defaults_file=defaults_file)

        fld%kind    = kind_
        fld%value   = value
        fld%file    = file
        fld%year    = year
        fld%current = value

        return
    end subroutine read_scalar

    subroutine aeros_boundary_init(bounds, grd)
        ! Allocate `buffer` for every 2D provider (a physics reader can
        ! assume `buffer` is allocated once init has run). For BC_FILE,
        ! either read the single time slice or read the 12-month clim into
        ! `clim` and seed `buffer` with the annual mean. For BC_CONSTANT,
        ! fill the buffer with `value`. BC_EXTERNAL leaves the buffer at
        ! zero -- a coupled host writes it in aeros_lib_step.

        implicit none
        type(aeros_boundary_t), intent(inout) :: bounds
        type(aeros_grid_class), intent(in)    :: grd

        call aeros_boundary_end(bounds)

        bounds%nlon = grd%nlon
        bounds%nlat = grd%nlat

        call init_2d(bounds%sst,         grd)
        call init_2d(bounds%sic,         grd)
        call init_2d(bounds%t_ice_sfc,   grd)
        call init_2d(bounds%t_skin_land, grd)
        call init_2d(bounds%q_skin_land, grd)
        call init_2d(bounds%albedo,      grd)
        call init_2d(bounds%z0m,         grd)
        call init_2d(bounds%phis,        grd)
        call init_2d(bounds%lsm,         grd)
        call init_2d(bounds%dust,        grd)
        call init_2d(bounds%so4,         grd)
        call init_2d(bounds%o3,          grd)

        ! Scalars: `current` reflects the CONSTANT value; FILE readers for
        ! CO2/CH4/N2O/TSI/time_bp are Phase C, so nothing to read at init.
        bounds%co2%current     = bounds%co2%value
        bounds%ch4%current     = bounds%ch4%value
        bounds%n2o%current     = bounds%n2o%value
        bounds%time_bp%current = bounds%time_bp%value
        bounds%tsi%current     = bounds%tsi%value

        return
    end subroutine aeros_boundary_init

    subroutine init_2d(fld, grd)
        implicit none
        type(aeros_bc_2d_t),    intent(inout) :: fld
        type(aeros_grid_class), intent(in)    :: grd

        integer :: mm

        allocate(fld%buffer(grd%nlon, grd%nlat))
        fld%buffer = 0.0_wp

        select case (fld%kind)
        case (BC_CONSTANT)
            fld%buffer = fld%value

        case (BC_FILE)
            if (len_trim(fld%file) == 0 .or. len_trim(fld%varname) == 0) then
                write(io_unit_err,*) "aeros_boundary_init:: BC_FILE with empty file/varname; leaving buffer zero"
                return
            end if
            if (fld%is_monthly) then
                allocate(fld%clim(grd%nlon, grd%nlat, 12))
                do mm = 1, 12
                    call aeros_bcinput_read_field(trim(fld%file), trim(fld%varname), &
                                                  grd%lon, grd%lat, fld%clim(:,:,mm), &
                                                  lon_name=trim(fld%lon_name), &
                                                  lat_name=trim(fld%lat_name), itime=mm)
                end do
                fld%clim   = max(fld%floor, fld%clim*fld%scale)
                fld%buffer = sum(fld%clim, dim=3)/12.0_wp
            else
                call aeros_bcinput_read_field(trim(fld%file), trim(fld%varname), &
                                              grd%lon, grd%lat, fld%buffer, &
                                              lon_name=trim(fld%lon_name), &
                                              lat_name=trim(fld%lat_name), itime=fld%itime)
                fld%buffer = max(fld%floor, fld%buffer*fld%scale)
            end if

        case (BC_EXTERNAL)
            ! Coupled host writes on the first aeros_lib_step; zero for now.
            fld%buffer = 0.0_wp

        case default
            write(io_unit_err,*) "aeros_boundary_init:: unknown BC kind ", fld%kind
            error stop 1
        end select

        return
    end subroutine init_2d

    subroutine aeros_boundary_refresh(bounds, doy)
        ! Refresh every 2D provider for the running day-of-year:
        !   BC_CONSTANT  -- no-op (buffer holds `value`).
        !   BC_FILE monthly -- daily interp of `clim` into `buffer`.
        !   BC_FILE single -- no-op (buffer was read at init).
        !   BC_EXTERNAL  -- no-op (host writes at cycle start; done outside
        !                   this routine, in aeros_lib_step's imports copy).
        implicit none
        type(aeros_boundary_t), intent(inout) :: bounds
        real(wp),               intent(in)    :: doy

        call refresh_2d(bounds%sst,         doy)
        call refresh_2d(bounds%sic,         doy)
        call refresh_2d(bounds%t_ice_sfc,   doy)
        call refresh_2d(bounds%t_skin_land, doy)
        call refresh_2d(bounds%q_skin_land, doy)
        call refresh_2d(bounds%albedo,      doy)
        call refresh_2d(bounds%z0m,         doy)
        call refresh_2d(bounds%phis,        doy)
        call refresh_2d(bounds%lsm,         doy)
        call refresh_2d(bounds%dust,        doy)
        call refresh_2d(bounds%so4,         doy)
        call refresh_2d(bounds%o3,          doy)

        return
    end subroutine aeros_boundary_refresh

    subroutine refresh_2d(fld, doy)
        implicit none
        type(aeros_bc_2d_t), intent(inout) :: fld
        real(wp),            intent(in)    :: doy

        if (fld%kind /= BC_FILE) return
        if (.not. fld%is_monthly) return
        if (.not. allocated(fld%clim)) return

        call aeros_boundary_interp_monthly(fld%clim, doy, fld%buffer)

        return
    end subroutine refresh_2d

    pure subroutine aeros_boundary_interp_monthly(clim, doy, out)
        ! Mid-month-centred linear cyclic interpolation of a monthly
        ! (nlon,nlat,12) climatology to day-of-year `doy`. Shared kernel
        ! matching aeros_land_set_albedo (aeros_land.f90:430-435),
        ! aeros_radiation_set_dust (aeros_radiation.f90:496-504) and
        ! rce_long's set_seasonal_sst.
        implicit none
        real(wp), intent(in)  :: clim(:,:,:)
        real(wp), intent(in)  :: doy
        real(wp), intent(out) :: out(:,:)

        real(wp), parameter :: DAYS_YR = 365.0_wp
        real(wp) :: x, w
        integer  :: im0, m0, m1

        x   = modulo(doy, DAYS_YR)/DAYS_YR*12.0_wp - 0.5_wp
        im0 = floor(x)
        w   = x - real(im0, wp)
        m0  = modulo(im0,     12) + 1
        m1  = modulo(im0 + 1, 12) + 1
        out = (1.0_wp - w)*clim(:,:,m0) + w*clim(:,:,m1)

        return
    end subroutine aeros_boundary_interp_monthly

    subroutine aeros_boundary_end(bounds)
        implicit none
        type(aeros_boundary_t), intent(inout) :: bounds

        call end_2d(bounds%sst)
        call end_2d(bounds%sic)
        call end_2d(bounds%t_ice_sfc)
        call end_2d(bounds%t_skin_land)
        call end_2d(bounds%q_skin_land)
        call end_2d(bounds%albedo)
        call end_2d(bounds%z0m)
        call end_2d(bounds%phis)
        call end_2d(bounds%lsm)
        call end_2d(bounds%dust)
        call end_2d(bounds%so4)
        call end_2d(bounds%o3)

        bounds%nlon = 0; bounds%nlat = 0
        return
    end subroutine aeros_boundary_end

    subroutine end_2d(fld)
        implicit none
        type(aeros_bc_2d_t), intent(inout) :: fld
        if (allocated(fld%buffer)) deallocate(fld%buffer)
        if (allocated(fld%clim))   deallocate(fld%clim)
        return
    end subroutine end_2d

end module aeros_boundary
