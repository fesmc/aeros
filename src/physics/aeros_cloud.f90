module aeros_cloud
    ! Diagnostic cloud scheme for the radiation. The model carries no prognostic
    ! cloud water or cloud fraction -- large-scale condensation removes the
    ! condensate immediately (aeros_condensation, section 9) -- so radiation's
    ! all-sky path (aeros_lw_cloudy_column, aeros_sw_cloudy_column) needs the
    ! cloud diagnosed from the resolved T/q/p column. This module is that
    ! diagnosis and nothing else consumes it yet.
    !
    ! Per layer, from temperature, specific humidity and pressure:
    !
    !   - Cloud fraction, RH-based (Sundqvist 1989):
    !         cf = 1 - sqrt(1 - b),   b = (RH - rhc)/(1 - rhc)  clamped to [0,1],
    !     with a critical relative humidity rhc(sigma). rhc is a constant by
    !     default but written as a two-endpoint vertical profile (surface, top)
    !     so it can be given height structure without a signature change.
    !
    !   - In-cloud condensate, a specified water CONTENT [kg m-3] converted to a
    !     mixing ratio by the air density rho = p/(R_d T):
    !         q_liq = (1-f_ice) LWC / rho,   q_ice = f_ice IWC / rho,
    !     with the liquid/ice split by a temperature ramp f_ice(T) (all liquid
    !     above CLD_T_LIQ, all ice below CLD_T_ICE, linear between). A specified
    !     content -- not q_ic proportional to q_sat -- is deliberate: tying the
    !     condensate to q_sat starves cold high cloud (q_sat tiny => no ice water
    !     => no longwave effect) while over-thickening warm low cloud, so the
    !     cloud radiative effect comes out ~2x too weak in the LW and ~2x too
    !     strong in the SW (validated on ERA5 columns). A content that is thin for
    !     ice (cirrus) and thicker for liquid gives the observed CRE balance.
    !
    !   - Outputs are the GRID-MEAN mixing ratios the radiation kernels expect
    !     (the kernels recover the in-cloud value by dividing by the column
    !     cloud fraction):
    !         clwc = cf (1-f_ice) LWC/rho,   ciwc = cf f_ice IWC/rho.
    !
    ! Grid-agnostic column routine, like the radiation kernels: an arbitrary nlev
    ! column, invariant to the grid it runs on. Parameters are named constants in
    ! one place (as the SESAM band coefficients are in aeros_radiation); the only
    ! user switch is radiation's `clouds` on/off flag. The RH_crit endpoints and
    ! water contents are tuned so the scheme, driven on ERA5 columns, reproduces
    ! ERA5's cloud radiative effect (m2_results §21): net CRE −23.6 vs −24.2, LW
    ! and SW each within 3 W/m². The one weak point is total cloud cover (0.47 vs
    ! 0.63) -- the scheme trades cover for optical depth but gets the CRE right.

    use aeros_defs,         only : wp, R_d
    use aeros_condensation, only : aeros_qsat

    implicit none
    private

    ! critical relative humidity profile endpoints (sigma=1 surface, sigma=0
    ! top); equal => constant with height, the default.
    real(wp), parameter :: CLD_RHC_SFC = 0.52_wp
    real(wp), parameter :: CLD_RHC_TOP = 0.35_wp

    ! in-cloud water content [kg m-3]; ice much thinner than liquid (cirrus)
    real(wp), parameter :: CLD_LWC = 0.017e-3_wp        ! liquid water content
    real(wp), parameter :: CLD_IWC = 0.018e-3_wp       ! ice water content

    ! liquid/ice partition temperature ramp [K]
    real(wp), parameter :: CLD_T_LIQ   = 273.15_wp     ! all liquid above
    real(wp), parameter :: CLD_T_ICE   = 250.15_wp     ! all ice below (-23 C)

    real(wp), parameter :: CLD_CF_MAX  = 0.999_wp      ! max diagnosed fraction

    public :: aeros_cloud_diagnose

contains

    subroutine aeros_cloud_diagnose(nlev, t, q, pfull, psfc, cf, clwc, ciwc)
        ! Diagnose grid-mean cloud fraction and cloud liquid/ice water for one
        ! column (any ordering; the routine is pointwise per layer).
        implicit none
        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: t(:)       ! (nlev) layer temperature [K]
        real(wp), intent(in)  :: q(:)       ! (nlev) specific humidity [kg kg-1]
        real(wp), intent(in)  :: pfull(:)   ! (nlev) layer pressure [Pa]
        real(wp), intent(in)  :: psfc       ! surface pressure [Pa]
        real(wp), intent(out) :: cf(:)      ! (nlev) cloud fraction [0-1]
        real(wp), intent(out) :: clwc(:)    ! (nlev) grid-mean cloud liquid [kg kg-1]
        real(wp), intent(out) :: ciwc(:)    ! (nlev) grid-mean cloud ice    [kg kg-1]

        integer  :: k
        real(wp) :: qs, dqsdt, rh, rhc, b, fice, rho

        do k = 1, nlev
            call aeros_qsat(t(k), pfull(k), qs, dqsdt)
            rh  = q(k)/max(qs, 1.0e-12_wp)
            rhc = rh_crit(pfull(k)/psfc)

            ! Sundqvist cloud fraction
            if (rhc < 1.0_wp) then
                b = (rh - rhc)/(1.0_wp - rhc)
                b = min(1.0_wp, max(0.0_wp, b))
                cf(k) = min(CLD_CF_MAX, 1.0_wp - sqrt(1.0_wp - b))
            else
                cf(k) = 0.0_wp
            end if

            ! in-cloud water content [kg m-3] -> mixing ratio via air density,
            ! then grid-mean water paths (x cloud fraction)
            if (cf(k) > 0.0_wp) then
                fice = ice_fraction(t(k))
                rho  = pfull(k)/(R_d*t(k))
                clwc(k) = cf(k)*(1.0_wp - fice)*CLD_LWC/rho
                ciwc(k) = cf(k)*fice*CLD_IWC/rho
            else
                clwc(k) = 0.0_wp
                ciwc(k) = 0.0_wp
            end if
        end do

        return

    contains

        pure real(wp) function rh_crit(sig) result(rhc)
            ! Critical RH, linear in sigma between the two endpoints. Constant by
            ! default (endpoints equal); the seam for a height-dependent profile.
            real(wp), intent(in) :: sig
            rhc = CLD_RHC_TOP + (CLD_RHC_SFC - CLD_RHC_TOP)*min(1.0_wp, max(0.0_wp, sig))
            return
        end function rh_crit

        pure real(wp) function ice_fraction(tk) result(fice)
            real(wp), intent(in) :: tk
            if (tk >= CLD_T_LIQ) then
                fice = 0.0_wp
            else if (tk <= CLD_T_ICE) then
                fice = 1.0_wp
            else
                fice = (CLD_T_LIQ - tk)/(CLD_T_LIQ - CLD_T_ICE)
            end if
            return
        end function ice_fraction

    end subroutine aeros_cloud_diagnose

end module aeros_cloud
