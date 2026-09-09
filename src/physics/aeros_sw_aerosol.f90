module aeros_sw_aerosol
    ! Background aerosol optics for the clear-sky shortwave, ported from
    ! CLIMBER-X (src/atm/swr.f90, the SESAM aerosol branch that the first
    ! aeros port dropped). One optical thickness `ot` and one imaginary
    ! refractive index `im`, uniform and constant -- CLIMBER-X's default
    ! (ars_ot 0.1, ars_im 0.01, `l_dust_rad = F`): a climate-independent
    ! background, equally valid for present day, climate change and paleo.
    ! A dust/AOD field (CLIMBER-X's dust_ot, its fake_dust_* files) replaces
    ! the constant behind the same two numbers when that is wanted.
    !
    ! Two effects, both from the SESAM fits (Petoukhov et al. 2000):
    !
    !   scattering  alb_sct_vu = 1 - (1 - r_scat) exp(-f1 f2 f3)
    !               alb_sct_ir = 1 - exp(-f1 f2 f3)
    !               f1 = cosz^p1, f2 = b_ar^p2, f3 = alf1 - alf2 ln(1 + alf3 im),
    !               b_ar = 0.55 ot  (the visible Rayleigh albedo r_scat is the
    !               ot = 0 limit; the near-IR has no Rayleigh term)
    !   absorption  itf_a(m) = exp(-gam1 b_ar m im^gam2)  on a path of airmass m
    !               (1/cosz down, 1.66 diffuse up), both bands.
    !
    ! With ot = 0 both reduce exactly to the aerosol-free expressions (f2 = 0,
    ! itf_a = 1), so the default is bit-for-bit the pre-aerosol scheme.
    !
    ! Why a constant and not an ERA5-fitted map: the offline validation
    ! (validate_era5, 2026-09-08) shows the scheme's clear-sky surface SW
    ! excess to be airmass-scaled and water-vapour-independent, uniform over
    ! the oceans and larger over the aerosol-laden continents. The uniform
    ! part is what this module carries; the continental excess is present-day
    ! dust and pollution and is left to a dust field, as in CLIMBER-X.
    use aeros_defs, only : wp
    implicit none
    private

    ! SESAM aerosol fit constants (CLIMBER-X swr.f90 module parameters)
    real(wp), parameter :: AER_P1    = -1.97_wp
    real(wp), parameter :: AER_P2    = 0.82_wp
    real(wp), parameter :: AER_ALF1  = 7.73e-2_wp
    real(wp), parameter :: AER_ALF2  = 2.39e-2_wp
    real(wp), parameter :: AER_ALF3  = 1.51e2_wp
    real(wp), parameter :: AER_GAM1  = 2.75_wp
    real(wp), parameter :: AER_GAM2  = 0.636_wp
    real(wp), parameter :: AER_B_OT  = 0.55_wp       ! b_ar = 0.55 ot

    ! Module state: the background aerosol. 0 = off (bit-for-bit clean sky);
    ! the calibrated production values live in the namelist defaults.
    real(wp), protected, public :: SW_AER_OT = 0.0_wp   ! optical thickness [-]
    real(wp), protected, public :: SW_AER_IM = 0.01_wp  ! imaginary refractive index [-]

    public :: aeros_sw_aerosol_configure
    public :: aeros_sw_aerosol_scat
    public :: aeros_sw_aerosol_itf

contains

    subroutine aeros_sw_aerosol_configure(ot, im)
        ! Set the background aerosol (any subset).
        implicit none
        real(wp), intent(in), optional :: ot, im
        if (present(ot)) SW_AER_OT = max(0.0_wp, ot)
        if (present(im)) SW_AER_IM = max(0.0_wp, im)
        return
    end subroutine aeros_sw_aerosol_configure

    pure subroutine aeros_sw_aerosol_scat(cosz, r_scat, alb_sct_vu, alb_sct_ir)
        ! Scattering (Rayleigh + aerosol) albedo of the atmosphere, visible and
        ! near-IR, for the direct beam at cosine zenith `cosz` (already floored
        ! by the caller). r_scat is the visible Rayleigh reflectance.
        implicit none
        real(wp), intent(in)  :: cosz, r_scat
        real(wp), intent(out) :: alb_sct_vu, alb_sct_ir
        real(wp) :: b_ar, f1, f2, f3, e
        b_ar = AER_B_OT*SW_AER_OT
        if (b_ar <= 0.0_wp) then
            alb_sct_vu = r_scat
            alb_sct_ir = 0.0_wp
            return
        end if
        f1 = cosz**AER_P1
        f2 = b_ar**AER_P2
        f3 = AER_ALF1 - AER_ALF2*log(1.0_wp + AER_ALF3*SW_AER_IM)
        e  = exp(-f1*f2*f3)
        alb_sct_vu = 1.0_wp - (1.0_wp - r_scat)*e
        alb_sct_ir = 1.0_wp - e
        return
    end subroutine aeros_sw_aerosol_scat

    pure real(wp) function aeros_sw_aerosol_itf(m) result(itf)
        ! Aerosol absorption transmission on a path of airmass factor `m`
        ! (1/cosz for the direct beam, 1/cosz + 1.66 for the TOA round trip,
        ! the two-path surface factors for the reflected beam).
        implicit none
        real(wp), intent(in) :: m
        real(wp) :: b_ar
        b_ar = AER_B_OT*SW_AER_OT
        if (b_ar <= 0.0_wp) then
            itf = 1.0_wp
            return
        end if
        itf = exp(-AER_GAM1*b_ar*m*SW_AER_IM**AER_GAM2)
        return
    end function aeros_sw_aerosol_itf

end module aeros_sw_aerosol
