module aeros_sw_aerosol
    ! Background aerosol optics for the clear-sky shortwave, ported from
    ! CLIMBER-X (src/atm/swr.f90, the SESAM aerosol branch that the first
    ! aeros port dropped). One optical thickness `ot` and one imaginary
    ! refractive index `im`, uniform and constant -- CLIMBER-X's default
    ! (ars_ot 0.1, ars_im 0.01, `l_dust_rad = F`): a climate-independent
    ! background, equally valid for present day, climate change and paleo.
    ! A dust optical-thickness field (CLIMBER-X's dust_ot; here the monthly
    ! maps of scripts/dust_ot_from_deposition.py from its fake_dust_* files)
    ! adds to the background per column (aeros_radiation dust_file).
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

    ! Sulfate direct effect (CLIMBER-X `l_so4_de`, Bauer et al. 2008): a
    ! non-absorbing scattering layer of load so4 [kg m-2] just above the
    ! surface, upscatter fraction beta_so4 and scattering efficiency sigma_so4
    ! [m2 kg-1]; its albedo beta sigma so4 / cosz combines with the surface
    ! albedo as alb_sa = alb_a + (1-alb_a)^2 alb_s / (1 - alb_a alb_s). The
    ! load is a field (CMIP6 historical/SSP ensemble mean in CLIMBER-X's
    ! so4_*.nc, time in years since 2000) or a constant: 0 for paleo.
    real(wp), protected, public :: SW_BETA_SO4  = 0.23_wp
    real(wp), protected, public :: SW_SIGMA_SO4 = 3.5e3_wp

    public :: aeros_sw_aerosol_configure
    public :: aeros_sw_aerosol_mix
    public :: aeros_sw_so4_albedo
    public :: aeros_sw_so4_cloud_factor
    public :: aeros_sw_aerosol_scat
    public :: aeros_sw_aerosol_itf

contains

    subroutine aeros_sw_aerosol_configure(ot, im, beta_so4, sigma_so4)
        ! Set the background aerosol (any subset).
        implicit none
        real(wp), intent(in), optional :: ot, im, beta_so4, sigma_so4
        if (present(ot)) SW_AER_OT = max(0.0_wp, ot)
        if (present(im)) SW_AER_IM = max(0.0_wp, im)
        if (present(beta_so4))  SW_BETA_SO4  = beta_so4
        if (present(sigma_so4)) SW_SIGMA_SO4 = sigma_so4
        return
    end subroutine aeros_sw_aerosol_configure

    pure subroutine aeros_sw_aerosol_mix(ot1, im1, ot2, im2, ot, im)
        ! Two aerosol species (background + dust) as one (ot, im) pair for the
        ! SESAM fits: optical thicknesses add, and the effective imaginary
        ! index conserves the absorption exponent, ot im^gam2 = sum ot_k im_k^gam2.
        implicit none
        real(wp), intent(in)  :: ot1, im1, ot2, im2
        real(wp), intent(out) :: ot, im
        real(wp) :: a
        ot = ot1 + ot2
        if (ot <= 0.0_wp) then
            im = im1
            return
        end if
        a  = (ot1*im1**AER_GAM2 + ot2*im2**AER_GAM2)/ot
        im = a**(1.0_wp/AER_GAM2)
        return
    end subroutine aeros_sw_aerosol_mix

    pure real(wp) function aeros_sw_so4_cloud_factor(so4) result(f)
        ! Sulfate indirect effect (CLIMBER-X `l_so4_ie`, clouds.f90): the
        ! cloud optical thickness scales with the droplet number as
        ! f_mod^(1/3), f_mod = (1 - exp(-alpha N))/(1 - exp(-alpha N_nat)),
        ! N the sulfate number at cloud base = N_ant,0 e^-1 + N_nat, where
        ! N_ant,0 = so4/(rho 4/3 pi r^3)/H is the anthropogenic surface number
        ! from the column load (r 0.17 um, rho 1769 kg/m3, H 1.5 km) and
        ! N_nat 75e6 m-3 the natural CCN. Applied to the liquid water content
        ! (optical thickness at fixed droplet size), so SW and LW both see it,
        ! as CLIMBER-X's clot does. so4 = 0 gives exactly 1.
        implicit none
        real(wp), intent(in) :: so4
        real(wp), parameter :: PI_ = 3.14159265358979323846_wp
        real(wp), parameter :: R_SO4 = 0.17e-6_wp, RHO_SO4 = 1.769e3_wp, H_SO4 = 1500.0_wp
        real(wp), parameter :: N_SO4_NAT = 75.0e6_wp, ALPHA_C = 2.5e-9_wp
        real(wp) :: n_ant0, n_so4, f_mod
        if (so4 <= 0.0_wp) then
            f = 1.0_wp
            return
        end if
        n_ant0 = so4/(RHO_SO4*4.0_wp/3.0_wp*PI_*R_SO4**3)/H_SO4
        n_so4  = n_ant0*exp(-1.0_wp) + N_SO4_NAT
        f_mod  = (1.0_wp - exp(-ALPHA_C*n_so4))/(1.0_wp - exp(-ALPHA_C*N_SO4_NAT))
        f = f_mod**0.33_wp
        return
    end function aeros_sw_so4_cloud_factor

    pure real(wp) function aeros_sw_so4_albedo(so4, cosz, alb_s) result(alb_sa)
        ! Combined surface + sulfate-layer albedo (Bauer et al. 2008 eq. 5-6).
        ! so4 <= 0 returns the surface albedo unchanged (bit-for-bit).
        implicit none
        real(wp), intent(in) :: so4, cosz, alb_s
        real(wp) :: alb_a
        if (so4 <= 0.0_wp) then
            alb_sa = alb_s
            return
        end if
        alb_a  = min(0.9_wp, SW_BETA_SO4*SW_SIGMA_SO4*so4/cosz)
        alb_sa = alb_a + (1.0_wp - alb_a)**2*alb_s/(1.0_wp - alb_a*alb_s)
        return
    end function aeros_sw_so4_albedo

    pure subroutine aeros_sw_aerosol_scat(ot, im, cosz, r_scat, alb_sct_vu, alb_sct_ir)
        ! Scattering (Rayleigh + aerosol) albedo of the atmosphere, visible and
        ! near-IR, for the direct beam at cosine zenith `cosz` (already floored
        ! by the caller), aerosol optical thickness `ot` (the background plus
        ! any dust). r_scat is the visible Rayleigh reflectance.
        implicit none
        real(wp), intent(in)  :: ot, im, cosz, r_scat
        real(wp), intent(out) :: alb_sct_vu, alb_sct_ir
        real(wp) :: b_ar, f1, f2, f3, e
        b_ar = AER_B_OT*ot
        if (b_ar <= 0.0_wp) then
            alb_sct_vu = r_scat
            alb_sct_ir = 0.0_wp
            return
        end if
        f1 = cosz**AER_P1
        f2 = b_ar**AER_P2
        f3 = AER_ALF1 - AER_ALF2*log(1.0_wp + AER_ALF3*im)
        e  = exp(-f1*f2*f3)
        alb_sct_vu = 1.0_wp - (1.0_wp - r_scat)*e
        alb_sct_ir = 1.0_wp - e
        return
    end subroutine aeros_sw_aerosol_scat

    pure real(wp) function aeros_sw_aerosol_itf(ot, im, m) result(itf)
        ! Aerosol absorption transmission of optical thickness `ot` on a path
        ! of airmass factor `m` (1/cosz for the direct beam, 1/cosz + 1.66 for
        ! the TOA round trip, the two-path surface factors for the reflected beam).
        implicit none
        real(wp), intent(in) :: ot, im, m
        real(wp) :: b_ar
        b_ar = AER_B_OT*ot
        if (b_ar <= 0.0_wp) then
            itf = 1.0_wp
            return
        end if
        itf = exp(-AER_GAM1*b_ar*m*im**AER_GAM2)
        return
    end function aeros_sw_aerosol_itf

end module aeros_sw_aerosol
