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
    ! The field itself comes from the phis boundary provider (aeros_boundary:
    ! ERA5 surface geopotential `z`, already m2 s-2 -- NOT an elevation, so no
    ! multiply by g -- or a coupler push). This module holds the two things
    ! every source of orography is put through before the dynamics sees it:
    !
    !   aeros_topography_bandlimit -- the spectral truncation to the model
    !                                 resolution (applied by the library
    !                                 whenever the provider buffer changes);
    !   aeros_topography_scale     -- the restart-safe spin-up ramp (applied
    !                                 by a driver that wants one).

    use aeros_defs,     only : dp, wp, wp_sh
    use aeros_spectral, only : aeros_sht_class, aeros_sht_analysis, aeros_sht_synthesis

    implicit none

    private

    real(wp), parameter :: SECONDS_PER_DAY = 86400.0_wp

    public :: aeros_topography_scale
    public :: aeros_topography_bandlimit

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

    subroutine aeros_topography_bandlimit(sht, phis)
        ! Spectrally truncate phis (in place, on the model grid) to the model
        ! resolution. A spectral core can only carry the surface geopotential
        ! at its resolved scales: the raw interpolated field still has grid-
        ! scale structure at coastlines and mountain flanks, which aliases in
        ! the transforms and rings (Gibbs), destabilizing the lowest layer.
        ! One analysis -> synthesis round-trip (SHTns truncates at the model
        ! truncation) gives the smooth orography the dynamics is consistent
        ! with. Standard spectral-model practice. Applied to every source of
        ! orography -- file or coupler -- alike.

        implicit none

        type(aeros_sht_class), intent(in)    :: sht
        real(wp),              intent(inout) :: phis(:,:)

        complex(wp_sh), allocatable :: phis_lm(:)
        real(dp),       allocatable :: gwork(:,:)

        allocate(phis_lm(sht%nlm))
        allocate(gwork(size(phis,1), size(phis,2)))
        gwork = real(phis, dp)
        call aeros_sht_analysis(sht, gwork, phis_lm)    ! grid -> spectral (overwrites gwork)
        call aeros_sht_synthesis(sht, phis_lm, gwork)   ! spectral -> band-limited grid
        phis = real(gwork, wp)
        deallocate(phis_lm, gwork)

        return

    end subroutine aeros_topography_bandlimit

end module aeros_topography
