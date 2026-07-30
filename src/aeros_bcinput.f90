module aeros_bcinput
    ! Boundary-condition input: read a prescribed 2D field from a netCDF file on
    ! a regular lon/lat grid and interpolate it onto the model's Gaussian grid.
    !
    ! This is deliberately GENERAL and topography-agnostic. It is the foundation
    ! for every changing lower-boundary field the model will eventually read from
    ! disk -- surface geopotential (topography), land-sea mask, surface albedo,
    ! prescribed SST, an ice-sheet surface. Each of those is "a named 2D field on
    ! a regular lon/lat source grid, wanted on the model grid", which is exactly
    ! the one contract this module provides. Keeping the netCDF read and the
    ! regridding here means a new boundary field is a one-line caller, not a new
    ! copy of the I/O and interpolation.
    !
    ! === What it assumes about the source ==================================
    !
    ! The source is a regular (separable) lon/lat grid: a 1D longitude coordinate
    ! and a 1D latitude coordinate, with the field varying on their tensor
    ! product. This covers ERA5, CMIP, and essentially every reanalysis/GCM
    ! boundary product. Longitude is treated as PERIODIC (the wrap between the
    ! last and first source columns is interpolated across, so a target point in
    ! the gap near 360 deg is handled correctly). Latitude may be stored either
    ! north-to-south (ERA5) or south-to-north (many models) -- both are detected
    ! and handled; the target latitude is clamped to the source range rather than
    ! extrapolated past the outermost source rows.
    !
    ! A leading time dimension (CF order time,lat,lon) is tolerated and collapsed
    ! by averaging over time, so a single-time orography file and a multi-month
    ! climatology both read the same way. ncio unpacks scale_factor/add_offset on
    ! read, so short-packed ERA5 fields come back in physical units.
    !
    ! === Interpolation =====================================================
    !
    ! Bilinear, self-contained (no external interpolation dependency): robust,
    ! exact on a constant field, and linear-order accurate -- adequate for a
    ! smooth lower-boundary field remapped from a finer source grid to T21/T42.
    ! A conservative or higher-order remap can replace the kernel here without
    ! touching any caller.

    use aeros_defs, only : dp, wp, io_unit_err
    use ncio,       only : nc_read, nc_size, nc_ndims, nc_dims

    implicit none

    private

    public :: aeros_bcinput_read_field
    public :: aeros_bcinput_bilinear

contains

    subroutine aeros_bcinput_read_field(filename, varname, tgt_lon, tgt_lat, field, &
                                        lon_name, lat_name, itime)
        ! Read `varname` from `filename` and bilinearly interpolate it onto the
        ! target grid defined by (tgt_lon [deg east], tgt_lat [deg north]),
        ! returning field(size(tgt_lon), size(tgt_lat)).
        !
        ! lon_name/lat_name name the source coordinate variables; they default to
        ! the CF-standard "longitude"/"latitude" (ERA5). The field variable may
        ! be 2D (lat,lon) or 3D (time,lat,lon). For a 3D field, `itime` selects a
        ! single time slice (1-based); when absent the time dimension is averaged
        ! out (bit-for-bit unchanged from before this argument existed).

        implicit none

        character(len=*), intent(in)  :: filename, varname
        real(wp),         intent(in)  :: tgt_lon(:), tgt_lat(:)
        real(wp),         intent(out) :: field(:,:)
        character(len=*), intent(in), optional :: lon_name, lat_name
        integer,          intent(in), optional :: itime

        character(len=64) :: lonnm, latnm
        integer  :: nsx, nsy, ndims, nt, it
        real(wp), allocatable :: lon_src(:), lat_src(:), src(:,:)
        real(wp), allocatable :: raw2(:,:), raw3(:,:,:)
        integer  :: i, j, nlon, nlat

        lonnm = "longitude"; latnm = "latitude"
        if (present(lon_name)) lonnm = lon_name
        if (present(lat_name)) latnm = lat_name

        nlon = size(tgt_lon); nlat = size(tgt_lat)
        if (size(field,1) /= nlon .or. size(field,2) /= nlat) then
            write(io_unit_err,*) "aeros_bcinput_read_field:: error: field shape ", &
                shape(field), " does not match target grid ", nlon, nlat
            error stop 1
        end if

        ! --- source coordinates -------------------------------------------
        nsx = nc_size(trim(filename), trim(lonnm))
        nsy = nc_size(trim(filename), trim(latnm))
        allocate(lon_src(nsx), lat_src(nsy))
        call nc_read(trim(filename), trim(lonnm), lon_src)
        call nc_read(trim(filename), trim(latnm), lat_src)

        ! --- source field, collapsing any leading time dimension ----------
        ndims = nc_ndims(trim(filename), trim(varname))
        allocate(src(nsx, nsy))
        select case (ndims)
        case (2)
            allocate(raw2(nsx, nsy))
            call nc_read(trim(filename), trim(varname), raw2)
            src = raw2
            deallocate(raw2)
        case (3)
            ! CF order (time,lat,lon) -> Fortran (lon,lat,time). ncio/nc_dims
            ! reports dims in Fortran (array) order, so the time dimension is
            ! LAST. Read (nlon,nlat,nt) and average over time.
            block
                integer, allocatable :: ddims(:)
                call nc_dims(trim(filename), trim(varname), dims=ddims)
                nt = ddims(size(ddims))     ! time == last (slowest) array dim
                deallocate(ddims)
            end block
            allocate(raw3(nsx, nsy, nt))
            call nc_read(trim(filename), trim(varname), raw3)
            if (present(itime)) then
                ! Single time slice (e.g. one month of a seasonal climatology).
                if (itime < 1 .or. itime > nt) then
                    write(io_unit_err,*) "aeros_bcinput_read_field:: error: itime ", &
                        itime, " out of range 1..", nt, " for '", trim(varname), "'"
                    error stop 1
                end if
                src = raw3(:,:,itime)
            else
                ! Time-average (default, bit-for-bit unchanged).
                src = raw3(:,:,1)
                do it = 2, nt
                    src = src + raw3(:,:,it)
                end do
                src = src/real(nt, wp)
            end if
            deallocate(raw3)
        case default
            write(io_unit_err,*) "aeros_bcinput_read_field:: error: variable '", &
                trim(varname), "' has unsupported rank ", ndims, &
                " (expected 2 or 3)"
            error stop 1
        end select

        ! --- regrid --------------------------------------------------------
        do j = 1, nlat
            do i = 1, nlon
                field(i,j) = aeros_bcinput_bilinear(lon_src, lat_src, src, &
                                                    tgt_lon(i), tgt_lat(j))
            end do
        end do

        deallocate(lon_src, lat_src, src)

        return

    end subroutine aeros_bcinput_read_field

    real(wp) function aeros_bcinput_bilinear(lon_src, lat_src, f, tlon, tlat) result(val)
        ! Bilinear interpolation of f(nx,ny) -- indexed f(ilon,ilat) on the
        ! separable source grid (lon_src, lat_src) -- to the point (tlon, tlat).
        !
        ! Longitude is periodic (period 360 deg); latitude is clamped to the
        ! source range. Both source coordinates must be strictly monotonic;
        ! latitude may run in either direction.

        implicit none

        real(wp), intent(in) :: lon_src(:), lat_src(:)
        real(wp), intent(in) :: f(:,:)
        real(wp), intent(in) :: tlon, tlat

        integer  :: nx, ny, i0, i1, j0, j1
        real(wp) :: wlon, wlat, xt, f00, f10, f01, f11

        nx = size(lon_src); ny = size(lat_src)

        ! --- longitude: periodic bracket ----------------------------------
        ! Wrap the target into [lon_src(1), lon_src(1)+360). Then either it sits
        ! between two interior columns, or in the wrap gap between the last and
        ! the (shifted) first column.
        xt = tlon
        do while (xt <  lon_src(1))          ; xt = xt + 360.0_wp; end do
        do while (xt >= lon_src(1) + 360.0_wp); xt = xt - 360.0_wp; end do

        if (xt >= lon_src(nx)) then
            ! wrap gap between column nx and column 1 (+360)
            i0 = nx; i1 = 1
            wlon = (xt - lon_src(nx)) / (lon_src(1) + 360.0_wp - lon_src(nx))
        else
            i0 = 1
            do while (i0 < nx .and. lon_src(i0+1) <= xt)
                i0 = i0 + 1
            end do
            i1 = i0 + 1
            wlon = (xt - lon_src(i0)) / (lon_src(i1) - lon_src(i0))
        end if

        ! --- latitude: clamp + bracket (either direction) -----------------
        call bracket_monotonic(lat_src, tlat, j0, j1, wlat)

        ! --- bilinear blend -----------------------------------------------
        f00 = f(i0,j0); f10 = f(i1,j0)
        f01 = f(i0,j1); f11 = f(i1,j1)
        val = (1.0_wp - wlon)*(1.0_wp - wlat)*f00 &
            +          wlon *(1.0_wp - wlat)*f10 &
            + (1.0_wp - wlon)*         wlat *f01 &
            +          wlon *         wlat *f11

        return

    end function aeros_bcinput_bilinear

    subroutine bracket_monotonic(x, xt, i0, i1, w)
        ! Locate xt within the strictly monotonic array x, returning the
        ! bracketing indices (i0,i1) and the fractional weight w in [0,1] such
        ! that the interpolant is (1-w)*f(i0) + w*f(i1). Values outside the range
        ! are clamped to the nearest edge (w pinned to 0 or 1), never
        ! extrapolated. Handles ascending and descending x.

        implicit none

        real(wp), intent(in)  :: x(:)
        real(wp), intent(in)  :: xt
        integer,  intent(out) :: i0, i1
        real(wp), intent(out) :: w

        integer :: n, k
        logical :: ascending

        n = size(x)
        if (n == 1) then
            i0 = 1; i1 = 1; w = 0.0_wp
            return
        end if

        ascending = (x(n) > x(1))

        if (ascending) then
            if (xt <= x(1)) then
                i0 = 1; i1 = 2; w = 0.0_wp; return
            else if (xt >= x(n)) then
                i0 = n-1; i1 = n; w = 1.0_wp; return
            end if
            k = 1
            do while (k < n-1 .and. x(k+1) <= xt)
                k = k + 1
            end do
        else
            if (xt >= x(1)) then
                i0 = 1; i1 = 2; w = 0.0_wp; return
            else if (xt <= x(n)) then
                i0 = n-1; i1 = n; w = 1.0_wp; return
            end if
            k = 1
            do while (k < n-1 .and. x(k+1) >= xt)
                k = k + 1
            end do
        end if

        i0 = k; i1 = k + 1
        w  = (xt - x(i0)) / (x(i1) - x(i0))

        return

    end subroutine bracket_monotonic

end module aeros_bcinput
