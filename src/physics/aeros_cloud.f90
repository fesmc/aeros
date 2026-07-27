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
    !   - In-cloud condensate, a temperature-dependent specified content:
    !         q_ic = CLD_QCFRAC * q_sat(T,p),
    !     so the condensate scales with the available vapour and falls off with
    !     height and cold as q_sat does. It is split into liquid and ice by a
    !     temperature ramp f_ice(T) (all liquid above CLD_T_LIQ, all ice below
    !     CLD_T_ICE, linear between).
    !
    !   - Outputs are the GRID-MEAN water paths the radiation kernels expect
    !     (the kernels recover the in-cloud value by dividing by the column
    !     cloud fraction):
    !         clwc = cf (1 - f_ice) q_ic,   ciwc = cf f_ice q_ic.
    !
    ! Grid-agnostic column routine, like the radiation kernels: an arbitrary nlev
    ! column, invariant to the grid it runs on. Parameters are named constants in
    ! one place (as the SESAM band coefficients are in aeros_radiation); the only
    ! user switch is radiation's `clouds` on/off flag. The values here are
    ! standard first-cut choices and are expected to need tuning so the
    ! diagnosed cloud climatology matches ERA5 cc / the cloud radiative effect.

    use aeros_defs,         only : wp
    use aeros_condensation, only : aeros_qsat

    implicit none
    private

    ! critical relative humidity profile endpoints (sigma=1 surface, sigma=0
    ! top); equal => constant with height, the default.
    real(wp), parameter :: CLD_RHC_SFC = 0.70_wp
    real(wp), parameter :: CLD_RHC_TOP = 0.70_wp

    ! in-cloud condensate as a fraction of the saturation specific humidity
    real(wp), parameter :: CLD_QCFRAC  = 0.04_wp

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
        real(wp) :: qs, dqsdt, rh, rhc, b, qic, fice

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

            ! in-cloud condensate, grid-mean water paths
            if (cf(k) > 0.0_wp) then
                qic  = CLD_QCFRAC * qs
                fice = ice_fraction(t(k))
                clwc(k) = cf(k)*(1.0_wp - fice)*qic
                ciwc(k) = cf(k)*fice*qic
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
