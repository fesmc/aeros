module aeros_boundary
    ! Boundary providers: the ONE source of every externally prescribed field.
    !
    ! Every lower-boundary or forcing field the physics consumes but does not
    ! compute -- the prescribed SST, the land-sea mask, the bare land albedo,
    ! the dust and sulfate maps, the surface geopotential -- has one provider
    ! here (spec docs/coupling_api.md sec 7). A provider has a discriminated
    ! `kind`:
    !
    !   BC_NONE      -- no provider: the consumer keeps its own internal
    !                   default (the ocean module's analytic SST profile, the
    !                   land module's scalar land_albedo, no dust, so4 = 0,
    !                   phis = 0). The buffer is allocated and zero.
    !   BC_CONSTANT  -- a scalar `value` (a uniform fill for 2D buffers).
    !   BC_FILE      -- read from `file` (via aeros_bcinput_read_field): a
    !                   single slice `itime` (0 = the time-mean of the file),
    !                   or a 12-month climatology (`is_monthly = .TRUE.`) held
    !                   in `clim(:,:,1..12)` and interpolated to the day into
    !                   `buffer(:,:)` by aeros_boundary_refresh.
    !   BC_EXTERNAL  -- the buffer is written by the caller (a coupled host
    !                   through aeros_lib_step(imports=...)); no I/O here.
    !
    ! The physics modules never read a file and never own a climatology: the
    ! library (aeros_boundary_apply in src/aeros.f90) copies each provider's
    ! `buffer` into the working field the physics kernel reads (ocn%sst,
    ! land%mask / land%albedo, rad%dust_ot / rad%so4_map, the dynamics' phis)
    ! after every refresh. A coupler therefore pushes any of these fields as
    ! BC_EXTERNAL with no physics change (Wave 4b).
    !
    ! `active`: a provider is materialised (file read, buffer applied) only
    ! when a consumer of it is on -- the library sets it from the &physics
    ! flags before init (l_land -> lsm/albedo, l_rad -> dust/so4,
    ! l_topography -> phis, &ocean mode = 0 -> sst), so a dry-core run does
    ! not depend on the ERA5 files its namelist defaults name.
    !
    ! `version`: bumped whenever the buffer is (re)written -- init, a monthly
    ! refresh, an import copy, a restart read -- so a consumer that does
    ! expensive work on a change (the phis band-limit) can tell a change from
    ! a repeat.
    !
    ! Zonal tables: a BC_FILE variable with no longitude dimension is taken
    ! as a zonal (lat[,month]) table already on the model latitudes and is
    ! broadcast along longitude -- the legacy aquaplanet sst(lat,month) from
    ! scripts/make_seasonal_sst.jl, available to any 2D field.
    !
    ! Bit-for-bit note: the daily interpolation kernel (mid-month centres,
    ! cyclic, x = mod(day, 365)/365*12 - 0.5) is the one rce_long's
    ! set_seasonal_sst, aeros_land_set_albedo and aeros_radiation_set_dust
    ! used before Wave 4b; the per-month regrid + max(floor, x*scale) matches
    ! the module readers they replace.

    use aeros_defs,    only : dp, wp, io_unit_err, aeros_grid_class
    use aeros_bcinput, only : aeros_bcinput_read_field
    use nml,           only : nml_read
    use ncio,          only : nc_read, nc_dims

    implicit none

    private

    integer, parameter, public :: BC_NONE     = -1
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
        integer               :: kind       = BC_NONE
        logical               :: active     = .TRUE.   ! a consumer is on (set by the library)
        integer               :: version    = 0        ! bumped on every buffer write
        real(wp), allocatable :: buffer(:,:)           ! (nlon,nlat), physics reads
        character(len=256)    :: file       = ""       ! FILE
        character(len=64)     :: varname    = ""       ! FILE variable name
        character(len=64)     :: lon_name   = "longitude"
        character(len=64)     :: lat_name   = "latitude"
        logical               :: is_monthly = .FALSE.  ! FILE: 12-month climatology
        real(wp), allocatable :: clim(:,:,:)           ! (nlon,nlat,12) FILE monthly
        integer               :: itime      = 1        ! FILE single-slice: 1-based; 0 = time-mean
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
        type(aeros_bc_2d_t) :: albedo        ! bare land broadband albedo [-]
        type(aeros_bc_2d_t) :: z0m           ! momentum roughness [m]
        type(aeros_bc_2d_t) :: phis          ! surface geopotential [m2 s-2]; band-limited IN PLACE by aeros_boundary_apply
        type(aeros_bc_2d_t) :: lsm           ! land-sea mask (0..1)
        type(aeros_bc_2d_t) :: dust          ! dust optical thickness
        type(aeros_bc_2d_t) :: so4           ! sulfate load [kg m-2]
        type(aeros_bc_2d_t) :: o3            ! ozone column (2D placeholder; 3D is Phase C)

        ! Surface-geopotential spin-up ramp [days]: linearly from 0 to full
        ! over this many days (aeros_topography_scale). 0 = no ramp. The
        ! library applies the ramped field in aeros_boundary_apply -- Wave 4d
        ! (was rce_long's driver-side scale) -- so a coupled host inherits it.
        real(wp) :: phis_ramp_days = 0.0_wp

        ! Wave 4d: refresh cadence de-duplication. aeros_lib_step calls
        ! refresh once per integer day rather than once per step; this holds
        ! the day of the last refresh so a repeated refresh on the same day
        ! is a no-op. -huge = "no refresh yet".
        integer :: last_refresh_int_doy = -huge(0)

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
    public :: aeros_boundary_kind_name

contains

    subroutine aeros_boundary_load(bounds, filename, grd, defaults_file)
        ! Read the `&boundary` namelist group and populate each provider's
        ! kind / file / varname / value. Then init (allocate buffers and read
        ! BC_FILE fields). The `active` flags are left as the caller set them.
        !
        ! Namelist convention: one key per provider field. For a 2D provider
        ! `X` the keys are X_kind, X_file, X_varname, X_is_monthly, X_itime,
        ! X_value, X_scale, X_floor, X_lon_name, X_lat_name; for a scalar,
        ! X_kind, X_value, X_file, X_year. Missing keys leave the compile-time
        ! defaults (BC_NONE for 2D fields: the consumer's own default).

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
        ! fill the buffer with `value`. BC_NONE and BC_EXTERNAL leave the
        ! buffer at zero -- a coupled host writes the latter in aeros_lib_step.
        ! An inactive provider (no consumer on) is not read.

        implicit none
        type(aeros_boundary_t), intent(inout) :: bounds
        type(aeros_grid_class), intent(in)    :: grd

        call aeros_boundary_end(bounds)

        bounds%nlon = grd%nlon
        bounds%nlat = grd%nlat

        call init_2d(bounds%sst,         "sst",         grd)
        call init_2d(bounds%sic,         "sic",         grd)
        call init_2d(bounds%t_ice_sfc,   "t_ice_sfc",   grd)
        call init_2d(bounds%t_skin_land, "t_skin_land", grd)
        call init_2d(bounds%q_skin_land, "q_skin_land", grd)
        call init_2d(bounds%albedo,      "albedo",      grd)
        call init_2d(bounds%z0m,         "z0m",         grd)
        call init_2d(bounds%phis,        "phis",        grd)
        call init_2d(bounds%lsm,         "lsm",         grd)
        call init_2d(bounds%dust,        "dust",        grd)
        call init_2d(bounds%so4,         "so4",         grd)
        call init_2d(bounds%o3,          "o3",          grd)

        ! Scalars: `current` reflects the CONSTANT value; FILE readers for
        ! CO2/CH4/N2O/TSI/time_bp are Phase C, so nothing to read at init.
        bounds%co2%current     = bounds%co2%value
        bounds%ch4%current     = bounds%ch4%value
        bounds%n2o%current     = bounds%n2o%value
        bounds%time_bp%current = bounds%time_bp%value
        bounds%tsi%current     = bounds%tsi%value

        return
    end subroutine aeros_boundary_init

    subroutine init_2d(fld, name, grd)
        implicit none
        type(aeros_bc_2d_t),    intent(inout) :: fld
        character(len=*),       intent(in)    :: name
        type(aeros_grid_class), intent(in)    :: grd

        integer :: mm

        allocate(fld%buffer(grd%nlon, grd%nlat))
        fld%buffer  = 0.0_wp
        fld%version = fld%version + 1

        if (.not. fld%active) return

        select case (fld%kind)
        case (BC_NONE)
            ! The consumer keeps its own default.

        case (BC_CONSTANT)
            fld%buffer = fld%value

        case (BC_FILE)
            if (len_trim(fld%file) == 0 .or. len_trim(fld%varname) == 0) then
                write(io_unit_err,*) "aeros_boundary_init:: error: "//trim(name)// &
                    "_kind = 1 (file) with an empty "//trim(name)//"_file / _varname."
                error stop 1
            end if
            if (fld%is_monthly) then
                allocate(fld%clim(grd%nlon, grd%nlat, 12))
                if (is_zonal_table(fld)) then
                    call read_zonal(fld, name, grd, fld%clim)
                else
                    do mm = 1, 12
                        call aeros_bcinput_read_field(trim(fld%file), trim(fld%varname), &
                                                      grd%lon, grd%lat, fld%clim(:,:,mm), &
                                                      lon_name=trim(fld%lon_name), &
                                                      lat_name=trim(fld%lat_name), itime=mm)
                    end do
                end if
                fld%clim   = max(fld%floor, fld%clim*fld%scale)
                fld%buffer = sum(fld%clim, dim=3)/12.0_wp
            else
                if (is_zonal_table(fld)) then
                    call read_zonal(fld, name, grd, fld%buffer)
                else if (fld%itime > 0) then
                    call aeros_bcinput_read_field(trim(fld%file), trim(fld%varname), &
                                                  grd%lon, grd%lat, fld%buffer, &
                                                  lon_name=trim(fld%lon_name), &
                                                  lat_name=trim(fld%lat_name), itime=fld%itime)
                else
                    ! itime = 0: the time-mean of the file (aeros_bcinput's
                    ! default when no slice is named).
                    call aeros_bcinput_read_field(trim(fld%file), trim(fld%varname), &
                                                  grd%lon, grd%lat, fld%buffer, &
                                                  lon_name=trim(fld%lon_name), &
                                                  lat_name=trim(fld%lat_name))
                end if
                fld%buffer = max(fld%floor, fld%buffer*fld%scale)
            end if

        case (BC_EXTERNAL)
            ! Coupled host writes on the first aeros_lib_step; zero for now.
            fld%buffer = 0.0_wp

        case default
            write(io_unit_err,*) "aeros_boundary_init:: unknown "//trim(name)//"_kind ", fld%kind
            error stop 1
        end select

        return
    end subroutine init_2d

    logical function is_zonal_table(fld) result(zonal)
        ! A BC_FILE variable with no longitude dimension is a zonal table
        ! (lat[,month]) already on the model latitudes.
        implicit none
        type(aeros_bc_2d_t), intent(in) :: fld
        character(len=64), allocatable :: names(:)
        integer :: i
        call nc_dims(trim(fld%file), trim(fld%varname), names=names)
        zonal = .TRUE.
        do i = 1, size(names)
            if (trim(names(i)) == trim(fld%lon_name)) zonal = .FALSE.
        end do
        deallocate(names)
        return
    end function is_zonal_table

    subroutine read_zonal(fld, name, grd, out)
        ! Read a zonal table -- (lat) for a single slice into out(:,:), or
        ! (lat,month) into out(:,:,1:12) -- and broadcast it along longitude.
        ! The table must be on the model latitudes: nlat rows (and 12 months).
        implicit none
        type(aeros_bc_2d_t),    intent(in)  :: fld
        character(len=*),       intent(in)  :: name
        type(aeros_grid_class), intent(in)  :: grd
        real(wp),               intent(out) :: out(..)

        integer, allocatable  :: dims(:)
        real(wp), allocatable :: zm1(:), zm2(:,:)
        integer :: ii, mm

        call nc_dims(trim(fld%file), trim(fld%varname), dims=dims)

        select rank (out)
        rank (2)
            if (size(dims) /= 1 .or. dims(1) /= grd%nlat) then
                write(io_unit_err,*) "aeros_boundary_init:: error: zonal "//trim(name)// &
                    " table must be (nlat) on the model grid, got dims ", dims
                error stop 1
            end if
            allocate(zm1(grd%nlat))
            call nc_read(trim(fld%file), trim(fld%varname), zm1)
            do ii = 1, grd%nlon
                out(ii, :) = zm1(:)
            end do
            deallocate(zm1)
        rank (3)
            if (size(dims) /= 2 .or. dims(1) /= grd%nlat .or. dims(2) /= 12) then
                write(io_unit_err,*) "aeros_boundary_init:: error: zonal monthly "//trim(name)// &
                    " table must be (nlat, 12) on the model grid, got dims ", dims
                error stop 1
            end if
            allocate(zm2(grd%nlat, 12))
            call nc_read(trim(fld%file), trim(fld%varname), zm2)
            do mm = 1, 12
                do ii = 1, grd%nlon
                    out(ii, :, mm) = zm2(:, mm)
                end do
            end do
            deallocate(zm2)
        rank default
            write(io_unit_err,*) "aeros_boundary_init:: read_zonal: unsupported rank"
            error stop 1
        end select

        deallocate(dims)
        return
    end subroutine read_zonal

    subroutine aeros_boundary_refresh(bounds, doy)
        ! Refresh every 2D provider for the running day-of-year:
        !   BC_FILE monthly -- daily interp of `clim` into `buffer` (version bump).
        !   everything else -- no-op (CONSTANT holds `value`; FILE single-slice
        !                   was read at init; EXTERNAL is written by the host
        !                   in aeros_lib_step's imports copy; NONE is empty).
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

        if (.not. fld%active) return
        if (fld%kind /= BC_FILE) return
        if (.not. fld%is_monthly) return
        if (.not. allocated(fld%clim)) return

        call aeros_boundary_interp_monthly(fld%clim, doy, fld%buffer)
        fld%version = fld%version + 1

        return
    end subroutine refresh_2d

    pure subroutine aeros_boundary_interp_monthly(clim, doy, out)
        ! Mid-month-centred linear cyclic interpolation of a monthly
        ! (nlon,nlat,12) climatology to day-of-year `doy`.
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

    pure function aeros_boundary_kind_name(kind) result(s)
        implicit none
        integer, intent(in) :: kind
        character(len=8) :: s
        select case (kind)
        case (BC_NONE);     s = "none"
        case (BC_CONSTANT); s = "constant"
        case (BC_FILE);     s = "file"
        case (BC_EXTERNAL); s = "external"
        case default;       s = "?"
        end select
        return
    end function aeros_boundary_kind_name

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
