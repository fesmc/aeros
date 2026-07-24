module aeros_grid
    ! The global Gaussian grid the spectral core lives on.
    !
    ! Everything here is derived from the SHTns configuration rather than
    ! recomputed independently, so there is exactly one source of truth for the
    ! Gauss nodes and weights. A grid that disagrees with the transform's own
    ! quadrature by even a rounding error breaks global conservation, which is
    ! the property docs/design.md section 3.7 (risk 2) and section 7 (M4) both
    ! require to hold to machine precision.

    use aeros_defs,     only : sp, dp, wp, aeros_grid_class, &
                                pi, r_earth, omega, radians_to_degrees
    use aeros_spectral, only : aeros_sht_class

    implicit none

    private

    public :: aeros_grid_init
    public :: aeros_grid_end

contains

    subroutine aeros_grid_init(grd, sht)
        ! Build the grid from a configured spherical-harmonic transform.

        implicit none

        type(aeros_grid_class), intent(inout) :: grd
        type(aeros_sht_class),  intent(in)    :: sht

        real(dp) :: dlon, cell
        integer  :: i, j

        grd%nlon = sht%nlon
        grd%nlat = sht%nlat
        grd%ncol = sht%nlon*sht%nlat

        allocate(grd%lon(grd%nlon))
        allocate(grd%lat(grd%nlat))
        allocate(grd%colat(grd%nlat))
        allocate(grd%sinlat(grd%nlat))
        allocate(grd%gauss_w(grd%nlat))
        allocate(grd%area(grd%nlon,grd%nlat))
        allocate(grd%coriolis(grd%nlon,grd%nlat))

        grd%colat   = sht%colat
        grd%sinlat  = sht%sinlat
        grd%gauss_w = sht%gauss_w

        do i = 1, grd%nlon
            grd%lon(i) = real(sht%lon(i)*radians_to_degrees, wp)
        end do

        ! Latitude runs NORTH TO SOUTH: SHTns' first colatitude row is theta
        ! near 0, i.e. the north pole. Reordering here would silently
        ! transpose every field the transforms write.
        !
        ! From sinlat (SHTns' Gauss nodes) rather than from colat: asin is odd,
        ! so lat comes out exactly antisymmetric about the equator. Going via
        ! `90 - colat` instead leaves a ~1e-14 degree asymmetry -- harmless in
        ! itself, but it puts a spurious hemispheric asymmetry into a quantity
        ! that later shows up in insolation and in every zonal-mean diagnostic.
        do j = 1, grd%nlat
            grd%lat(j) = real(asin(sht%sinlat(j))*radians_to_degrees, wp)
        end do

        ! Cell area: the Gauss weight IS the latitude-band measure, so the cell
        ! area is exact by construction and sum(area) = 4*pi*a^2 to round-off.
        ! Do not substitute a cos(lat)*dlat*dlon formula, which is only
        ! second-order accurate and would break that identity.
        dlon = 2.0_dp*pi/real(grd%nlon, dp)
        do j = 1, grd%nlat
            cell = grd%gauss_w(j)*dlon*r_earth*r_earth
            do i = 1, grd%nlon
                grd%area(i,j) = real(cell, wp)
            end do
        end do

        ! f = 2*omega*sin(lat), and sinlat IS sin(lat) -- no trig call, and
        ! exactly antisymmetric, so f sums to zero over the globe.
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                grd%coriolis(i,j) = real(2.0_dp*omega*sht%sinlat(j), wp)
            end do
        end do

        return

    end subroutine aeros_grid_init

    subroutine aeros_grid_end(grd)

        implicit none

        type(aeros_grid_class), intent(inout) :: grd

        grd%nlon = 0; grd%nlat = 0; grd%ncol = 0

        if (allocated(grd%lon))      deallocate(grd%lon)
        if (allocated(grd%lat))      deallocate(grd%lat)
        if (allocated(grd%colat))    deallocate(grd%colat)
        if (allocated(grd%sinlat))   deallocate(grd%sinlat)
        if (allocated(grd%gauss_w))  deallocate(grd%gauss_w)
        if (allocated(grd%area))     deallocate(grd%area)
        if (allocated(grd%coriolis)) deallocate(grd%coriolis)

        return

    end subroutine aeros_grid_end

end module aeros_grid
