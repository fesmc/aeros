module aeros_topography
    ! Real surface topography (orography) as the model's lower boundary.
    !
    ! The dynamics enters the ground through the surface geopotential phis =
    ! g*z_s [m2 s-2], the geopotential at the bottom of the hydrostatic
    ! integration (aeros_hydrostatic sets ph(nlev) = phis). An aquaplanet uses
    ! phis = 0 everywhere; real orography replaces that flat surface with the
    ! Earth's, lifting the pressure surfaces over the continents and mountains
    ! and letting the flow feel them.
    !
    ! This module is a thin, named front door: given a file it returns phis on
    ! the model grid. It leans entirely on aeros_bcinput for the netCDF read and
    ! the regridding -- the topography-specific knowledge here is only that the
    ! ERA5 surface-geopotential variable IS phis directly (units m2 s-2, already
    ! a geopotential -- NOT an elevation, so there is no multiply by g), stored
    ! under the short name "z" in the ERA5 single-level convention.
    !
    ! A separate module (rather than an inline read in the driver) so that every
    ! driver that wants orography -- rce_long today, the coupled model later --
    ! shares one definition of "load the Earth's topography", and so the ERA5
    ! variable-name convention lives in exactly one place.

    use aeros_defs,    only : wp
    use aeros_bcinput, only : aeros_bcinput_read_field

    implicit none

    private

    real(wp), parameter :: SECONDS_PER_DAY = 86400.0_wp

    public :: aeros_topography_load
    public :: aeros_topography_scale

contains

    pure real(wp) function aeros_topography_scale(l_topo, time_seconds, ramp_days) result(s)
        ! Scaling applied to the full topography before it is handed to the
        ! dynamics: phis(t) = aeros_topography_scale(...) * phis_full.
        !
        ! It is a PURE FUNCTION OF ELAPSED MODEL TIME -- no internal accumulator
        ! -- so it is restart-safe: a checkpoint that restores (nstep, time)
        ! resumes the ramp at exactly the right value with no extra state to
        ! save. The caller passes elapsed time = nstep*dt.
        !
        !   l_topo = .false.  -> 0 always (flat aquaplanet, bit-for-bit unchanged)
        !   ramp_days <= 0     -> 1 always (full topography from t = 0)
        !   otherwise          -> linear 0 -> 1 over ramp_days, then held at 1
        !
        ! The ramp exists to spread the spin-up shock of switching on mountains:
        ! raising the Earth's orography instantaneously under a balanced flow
        ! launches large gravity waves, so the surface is lifted gradually.

        implicit none

        logical,  intent(in) :: l_topo
        real(wp), intent(in) :: time_seconds
        real(wp), intent(in) :: ramp_days

        real(wp) :: r

        if (.not. l_topo) then
            s = 0.0_wp
            return
        end if

        if (ramp_days <= 0.0_wp) then
            s = 1.0_wp
            return
        end if

        r = time_seconds / (ramp_days*SECONDS_PER_DAY)
        s = min(max(r, 0.0_wp), 1.0_wp)

        return

    end function aeros_topography_scale

    subroutine aeros_topography_load(filename, tgt_lon, tgt_lat, phis, varname)
        ! Return phis(size(tgt_lon), size(tgt_lat)) [m2 s-2] read from `filename`
        ! and interpolated onto the model grid (tgt_lon [deg east], tgt_lat [deg
        ! north]). `varname` defaults to "z" (ERA5 surface geopotential); the
        ! value is used as phis with no unit conversion.

        implicit none

        character(len=*), intent(in)  :: filename
        real(wp),         intent(in)  :: tgt_lon(:), tgt_lat(:)
        real(wp),         intent(out) :: phis(:,:)
        character(len=*), intent(in), optional :: varname

        character(len=64) :: vnm

        vnm = "z"
        if (present(varname)) vnm = varname

        call aeros_bcinput_read_field(trim(filename), trim(vnm), &
                                      tgt_lon, tgt_lat, phis)

        return

    end subroutine aeros_topography_load

end module aeros_topography
