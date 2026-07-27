module aeros_insolation
    ! Top-of-atmosphere insolation for aeros via the fesmc/insol package (Laskar
    ! et al. 2004 orbital elements, ~/models/insol). Replaces the earlier
    ! present-day, circular, obliquity-only stopgap (aeros_radiation's
    ! aeros_insolation_daily, retained there only as a unit-test reference):
    ! insol gives the true seasonal cycle AND Milankovitch/orbital forcing through
    ! time_bp -- the paleoclimate axis this model exists to resolve.
    !
    ! Two fields per latitude are handed to the shortwave scheme:
    !   sw_toa  daily-mean TOA downward SW [W m-2]
    !   coszen  flux-weighted mean cosine of the solar zenith angle [-]
    ! (coszm_kind = "flux", i.e. the insolation-weighted airmass cosine the band
    ! shortwave wants). The insol state holds the orbital tables and the grid
    ! geometry; it is built once at init and read-only after, so it is safe to
    ! share between threads.

    use aeros_defs,  only : wp, dp, S0, aeros_grid_class
    use insolation,  only : insol_class, insol_init, insol_end, &
                            calc_insol_day, calc_insol_ave

    implicit none

    private

    ! Radiative calendar: the seasonal cycle repeats every DAY_YEAR days, and the
    ! model's day-of-year is taken modulo this. 365 keeps the season lengths ~real
    ! (insol accepts 360, 365 or 366).
    integer, parameter, public :: DAY_YEAR = 365

    type, public :: aeros_insol_class
        type(insol_class) :: is
        logical  :: ready   = .FALSE.
        real(dp) :: time_bp = 0.0_dp        ! orbital year before present (1950 CE)
        integer  :: nlat    = 0
    end type aeros_insol_class

    public :: aeros_insol_init
    public :: aeros_insol_end
    public :: aeros_insol_annual
    public :: aeros_insol_day

contains

    subroutine aeros_insol_init(ins, grd, time_bp, fldr)
        ! Load the orbital tables and precompute the grid geometry once. `fldr` is
        ! where the LA2004 tables live (default insol/input (the insol symlink), relative to the run
        ! directory).
        type(aeros_insol_class), intent(inout) :: ins
        type(aeros_grid_class),  intent(in)    :: grd
        real(dp),         intent(in)           :: time_bp
        character(len=*), intent(in), optional :: fldr

        real(dp) :: lats(grd%nlat)
        character(len=256) :: dir

        call aeros_insol_end(ins)

        ! LA2004 orbital tables ship with the insol package; read them through the
        ! `insol` symlink (the same one that provides libinsol), so nothing is
        ! duplicated into the aeros tree.
        dir = "insol/input"
        if (present(fldr)) dir = fldr

        lats = real(grd%lat(1:grd%nlat), dp)
        call insol_init(ins%is, lats, fldr=trim(dir), S0=real(S0, dp), &
                        day_year=DAY_YEAR, coszm_kind="flux")

        ins%time_bp = time_bp
        ins%nlat    = grd%nlat
        ins%ready   = .TRUE.

        return
    end subroutine aeros_insol_init

    subroutine aeros_insol_end(ins)
        type(aeros_insol_class), intent(inout) :: ins
        if (ins%ready) call insol_end(ins%is)
        ins%ready = .FALSE.
        ins%nlat  = 0
        return
    end subroutine aeros_insol_end

    subroutine aeros_insol_annual(ins, sw_toa, coszen)
        ! Annual-mean daily insolation and the year's flux-weighted mean cosine
        ! zenith per latitude -- the annual-mean insolation mode. The mean cosine
        ! is insolation-weighted (calc_insol_ave returns only the flux, so the
        ! weighting is done here, day by day).
        type(aeros_insol_class), intent(in)  :: ins
        real(wp), intent(out) :: sw_toa(:), coszen(:)

        real(wp) :: sw(ins%nlat), cz(ins%nlat), wsum(ins%nlat)
        integer  :: d

        call calc_insol_ave(ins%is, sw_toa, [1, DAY_YEAR], ins%time_bp)

        coszen = 0.0_wp; wsum = 0.0_wp
        do d = 1, DAY_YEAR
            call calc_insol_day(ins%is, sw, d, ins%time_bp, coszm=cz)
            coszen = coszen + sw*cz
            wsum   = wsum + sw
        end do
        where (wsum > 0.0_wp)
            coszen = coszen/wsum
        elsewhere
            coszen = 0.0_wp
        end where

        return
    end subroutine aeros_insol_annual

    subroutine aeros_insol_day(ins, sw_toa, coszen, doy)
        ! Daily-mean insolation and flux-weighted mean cosine zenith for the given
        ! day of year, per latitude -- the seasonal mode, called on the radiation
        ! recompute cadence.
        type(aeros_insol_class), intent(in)  :: ins
        real(wp), intent(out) :: sw_toa(:), coszen(:)
        real(wp), intent(in)  :: doy

        integer :: day
        day = modulo(nint(doy) - 1, DAY_YEAR) + 1     ! -> [1, DAY_YEAR]
        call calc_insol_day(ins%is, sw_toa, day, ins%time_bp, coszm=coszen)

        return
    end subroutine aeros_insol_day

end module aeros_insolation
