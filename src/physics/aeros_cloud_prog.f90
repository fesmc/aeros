module aeros_cloud_prog
    ! Prognostic cloud fraction (Sundqvist 1989), the opt-in replacement for the
    ! diagnostic RH->cover scheme (aeros_cloud). The diagnostic scheme slaves the
    ! cloud fraction to the instantaneous grid-box relative humidity, so once the
    ! model's free troposphere runs moist the cover saturates and runs away
    ! (0.66 -> 0.86 in the coupled slab RCE): a pure function of RH cannot lag,
    ! and there is no sink to balance a positive cloud-RH-moisture feedback.
    !
    ! Here the cloud fraction is a genuine PROGNOSTIC: a gridpoint field
    ! (aeros_state_class%cf_g, off the spectral core for the same reason
    ! humidity is -- truncating a bounded [0,1] field overshoots out of range)
    ! carrying its own budget in time,
    !
    !     dC/dt = formation - evaporation + convective detrainment ,
    !
    ! integrated once per dynamics step at the grid seam, right after convection
    ! and condensation have set the grid-box T and q. The budget has a sink, so
    ! the equilibrium cover sits BELOW the diagnostic Sundqvist value wherever
    ! the air is subsaturated, and the finite formation timescale breaks the
    ! instantaneous feedback -- which is the whole point.
    !
    ! === The three terms ====================================================
    !
    !   FORMATION relaxes C up toward the Sundqvist diagnostic target
    !
    !       C_eq = 1 - sqrt( (1 - RH)/(1 - RH_crit) ) ,   RH >= RH_crit ,
    !
    !     over a formation timescale tau_form. RH_crit(sigma) is the classic
    !     critical-RH profile, higher aloft than near the surface. Only the
    !     upward pull (C rising toward C_eq) is formation; C_eq is the ceiling
    !     the resolved humidity can support, reached in the limit tau_form -> 0
    !     (which recovers the diagnostic scheme with no sink).
    !
    !   EVAPORATION relaxes C down toward zero at a rate scaled by the grid-box
    !     SUBSATURATION (1 - RH): cloud in unsaturated air evaporates, and the
    !     drier the air the faster. This is the sink the diagnostic scheme
    !     lacks. At exact saturation (RH = 1) it vanishes, so a saturated column
    !     is pulled cleanly to overcast by the formation term alone.
    !
    !   CONVECTIVE DETRAINMENT adds cloud in convecting columns, tied to the
    !     convection scheme's precipitation (its measure of convective
    !     intensity): an anvil source in the free troposphere that the resolved
    !     RH would otherwise miss. Off when convection is off (no precip).
    !
    ! The update is written as an exponential relaxation over the step, so each
    ! term moves C by a bounded fraction of its gap and C provably stays in
    ! [0, CLD_CF_MAX] from any C in that range -- no clipping needed for
    ! stability, only the final guard against round-off.
    !
    ! === Scope: cloud fraction only, in-cloud water diagnosed ===============
    !
    ! This carries the cloud FRACTION prognostically but NOT the in-cloud
    ! condensate. Radiation's cloud optics (aeros_cloud_water) already convert a
    ! grid-mean cloud fraction into liquid/ice water paths through a specified
    ! in-cloud water content tuned to reproduce the ERA5 cloud radiative effect
    ! (m2_results s21); carrying a prognostic q_c as well would mean re-deriving
    ! and re-tuning those optics and doubling the new restart state. So the
    ! prognostic scheme replaces only the FRACTION -- the runaway quantity --
    ! and drives the identical optics as the diagnostic path. Prognostic q_c is
    ! a later refinement, at which point Sundqvist's C-q_c coupling is added.

    use aeros_defs,         only : wp, io_unit_err, aeros_grid_class
    use aeros_vertical,     only : aeros_vgrid_class, aeros_vgrid_pressure
    use aeros_condensation, only : aeros_qsat
    use nml,                only : nml_read

    implicit none
    private

    ! Reference saturation RH the Sundqvist target is measured against.
    real(wp), parameter :: CLD_RH_SAT   = 1.0_wp

    ! Maximum cloud fraction (matches the diagnostic scheme's ceiling).
    real(wp), parameter :: CLD_CF_MAX   = 0.999_wp

    ! Sigma above which convective detrainment deposits cloud (free troposphere;
    ! below this is the sub-cloud boundary layer, no anvil source).
    real(wp), parameter :: CLD_DETR_SIG = 0.8_wp

    ! Convective precip [kg m-2 s-1] that saturates the detrainment source
    ! (~8.6 mm/day, a solidly convecting column).
    real(wp), parameter :: CLD_PRECIP_REF = 1.0e-4_wp

    type aeros_cloud_prog_class
        ! Prognostic cloud-fraction scheme configuration. Allocation-free: the
        ! state it evolves (cf_g) lives in aeros_state_class, so this holds only
        ! scalars and never needs an allocate/deallocate cycle.
        logical  :: enabled = .FALSE.
        integer  :: nlon = 0, nlat = 0

        ! Critical-RH profile endpoints (sigma = 1 surface, sigma = 0 top). The
        ! classic profile is higher aloft than near the surface.
        real(wp) :: rhc_sfc = 0.70_wp
        real(wp) :: rhc_top = 0.90_wp

        real(wp) :: tau_form = 3600.0_wp   ! formation relaxation timescale [s]
        real(wp) :: tau_evap = 3600.0_wp   ! evaporation relaxation timescale [s]
        real(wp) :: c_detr   = 0.5_wp      ! convective detrainment anvil ceiling [-]
    end type aeros_cloud_prog_class

    public :: aeros_cloud_prog_class
    public :: aeros_cloud_prog_init
    public :: aeros_cloud_prog_load
    public :: aeros_cloud_prog_end
    public :: aeros_cloud_prog_apply
    public :: aeros_cloud_prog_report

contains

    subroutine aeros_cloud_prog_init(cpr, grd, enabled)
        ! Set the grid sizes and the (disabled-by-default) configuration. The
        ! driver overrides the individual knobs after this, exactly as it does
        ! for the other grid-seam physics.
        implicit none
        type(aeros_cloud_prog_class), intent(inout) :: cpr
        type(aeros_grid_class),       intent(in)    :: grd
        logical,                      intent(in)    :: enabled

        cpr%enabled = enabled
        cpr%nlon    = grd%nlon
        cpr%nlat    = grd%nlat

        return
    end subroutine aeros_cloud_prog_init

    subroutine aeros_cloud_prog_load(cpr, filename, grd, defaults_file)
        ! Read the scheme configuration from the `aeros_cloud` namelist group,
        ! mirroring the other physics `_load` routines.
        implicit none
        type(aeros_cloud_prog_class), intent(inout) :: cpr
        character(len=*),             intent(in)    :: filename
        type(aeros_grid_class),       intent(in)    :: grd
        character(len=*), intent(in), optional      :: defaults_file

        logical  :: enabled
        real(wp) :: rhc_sfc, rhc_top, tau_form, tau_evap, c_detr

        enabled  = .FALSE.
        rhc_sfc  = 0.70_wp
        rhc_top  = 0.90_wp
        tau_form = 3600.0_wp
        tau_evap = 3600.0_wp
        c_detr   = 0.5_wp

        call nml_read(filename, "aeros_cloud", "l_cloud_prog",  enabled,  defaults_file=defaults_file)
        call nml_read(filename, "aeros_cloud", "cloud_rhc_sfc", rhc_sfc,  defaults_file=defaults_file)
        call nml_read(filename, "aeros_cloud", "cloud_rhc_top", rhc_top,  defaults_file=defaults_file)
        call nml_read(filename, "aeros_cloud", "cloud_tau_form", tau_form, defaults_file=defaults_file)
        call nml_read(filename, "aeros_cloud", "cloud_tau_evap", tau_evap, defaults_file=defaults_file)
        call nml_read(filename, "aeros_cloud", "cloud_c_detr",  c_detr,   defaults_file=defaults_file)

        call aeros_cloud_prog_init(cpr, grd, enabled)
        cpr%rhc_sfc  = rhc_sfc
        cpr%rhc_top  = rhc_top
        cpr%tau_form = tau_form
        cpr%tau_evap = tau_evap
        cpr%c_detr   = c_detr

        call validate(cpr)

        return
    end subroutine aeros_cloud_prog_load

    subroutine validate(cpr)
        implicit none
        type(aeros_cloud_prog_class), intent(in) :: cpr
        if (cpr%rhc_sfc <= 0.0_wp .or. cpr%rhc_sfc >= 1.0_wp .or. &
            cpr%rhc_top <= 0.0_wp .or. cpr%rhc_top >= 1.0_wp) then
            write(io_unit_err,*) "aeros_cloud_prog:: error: cloud_rhc_* must be in (0,1), got ", &
                                    cpr%rhc_sfc, cpr%rhc_top
            error stop 1
        end if
        if (cpr%tau_form <= 0.0_wp .or. cpr%tau_evap <= 0.0_wp) then
            write(io_unit_err,*) "aeros_cloud_prog:: error: cloud_tau_* must be > 0 s, got ", &
                                    cpr%tau_form, cpr%tau_evap
            error stop 1
        end if
        if (cpr%c_detr < 0.0_wp) then
            write(io_unit_err,*) "aeros_cloud_prog:: error: cloud_c_detr must be >= 0, got ", cpr%c_detr
            error stop 1
        end if
        return
    end subroutine validate

    subroutine aeros_cloud_prog_end(cpr)
        implicit none
        type(aeros_cloud_prog_class), intent(inout) :: cpr
        cpr%enabled = .FALSE.
        cpr%nlon = 0; cpr%nlat = 0
        return
    end subroutine aeros_cloud_prog_end

    subroutine aeros_cloud_prog_apply(cpr, vg, t_g, qv_g, lnps_g, cf_g, dt, precip)
        ! Advance the prognostic cloud fraction one step, in place, at the grid
        ! seam. t_g/qv_g are the post-convection, post-condensation grid-box
        ! temperature and humidity (so RH is the resolved saturation state the
        ! step ends in); precip (optional) is the convection scheme's column
        ! precipitation, the convective-intensity proxy for detrainment.
        implicit none
        type(aeros_cloud_prog_class), intent(in)    :: cpr
        type(aeros_vgrid_class),      intent(in)    :: vg
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(in)    :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat)      ln[Pa]
        real(wp), intent(inout) :: cf_g(:,:,:)    ! (nlon,nlat,nlev) [0-1]
        real(wp), intent(in)    :: dt             ! [s]
        real(wp), intent(in), optional :: precip(:,:)  ! (nlon,nlat) [kg m-2 s-1]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: psfc, qs, dqsdt, rh, rhc, sig, cf_eq, b
        real(wp) :: aform, aevap, sub, cf, conv_act, adetr
        logical  :: have_precip
        integer  :: i, j, k, nlev

        if (.not. cpr%enabled) return

        nlev = vg%nlev
        have_precip = present(precip)

        ! Per-step relaxation fractions from the timescales (bounded in [0,1)).
        aform = 1.0_wp - exp(-dt/cpr%tau_form)
        aevap = 1.0_wp - exp(-dt/cpr%tau_evap)

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,phalf,pfull,dpc,psfc,qs,dqsdt,rh,rhc,sig,cf_eq,b, &
        !$omp           sub,cf,conv_act,adetr)
        do j = 1, cpr%nlat
            do i = 1, cpr%nlon
                psfc = exp(lnps_g(i,j))
                call aeros_vgrid_pressure(vg, psfc, phalf, pfull, dpc)

                conv_act = 0.0_wp
                if (have_precip) &
                    conv_act = min(1.0_wp, max(0.0_wp, precip(i,j)/CLD_PRECIP_REF))

                do k = 1, nlev
                    call aeros_qsat(t_g(i,j,k), pfull(k), qs, dqsdt)
                    rh  = qv_g(i,j,k)/max(qs, 1.0e-12_wp)
                    sig = pfull(k)/psfc
                    rhc = rh_crit(cpr, sig)

                    ! Sundqvist diagnostic target (the formation ceiling).
                    if (rh >= rhc .and. rhc < CLD_RH_SAT) then
                        b     = (CLD_RH_SAT - rh)/(CLD_RH_SAT - rhc)
                        cf_eq = 1.0_wp - sqrt(max(0.0_wp, b))
                        cf_eq = min(CLD_CF_MAX, max(0.0_wp, cf_eq))
                    else
                        cf_eq = 0.0_wp
                    end if

                    cf = cf_g(i,j,k)

                    ! Formation: relax up toward the target (upward pull only).
                    if (cf_eq > cf) cf = cf + (cf_eq - cf)*aform

                    ! Evaporation: relax down toward zero, scaled by the
                    ! subsaturation (zero at saturation, full below RH_crit).
                    sub = (CLD_RH_SAT - rh)/max(CLD_RH_SAT - rhc, 1.0e-3_wp)
                    sub = min(1.0_wp, max(0.0_wp, sub))
                    cf  = cf - cf*sub*aevap

                    ! Convective detrainment: an anvil source in the convecting
                    ! free troposphere, tied to the column convective intensity.
                    ! It relaxes cf up toward a bounded anvil target
                    ! c_detr*conv_act (a source only -- never a sink), so a
                    ! strongly convecting column reaches an anvil cover of at most
                    ! c_detr, NOT overcast: relaxing every convecting layer toward
                    ! 1 would blanket the tropics and defeat the whole point.
                    if (conv_act > 0.0_wp .and. sig < CLD_DETR_SIG) then
                        adetr = min(CLD_CF_MAX, cpr%c_detr*conv_act)
                        if (adetr > cf) cf = cf + (adetr - cf)*aform
                    end if

                    cf_g(i,j,k) = min(CLD_CF_MAX, max(0.0_wp, cf))
                end do
            end do
        end do
        !$omp end parallel do

        return
    end subroutine aeros_cloud_prog_apply

    pure real(wp) function rh_crit(cpr, sig) result(rhc)
        ! Critical RH, linear in sigma between the two endpoints (higher aloft).
        implicit none
        type(aeros_cloud_prog_class), intent(in) :: cpr
        real(wp),                     intent(in) :: sig
        real(wp) :: s
        s = min(1.0_wp, max(0.0_wp, sig))
        rhc = cpr%rhc_top + (cpr%rhc_sfc - cpr%rhc_top)*s
        return
    end function rh_crit

    subroutine aeros_cloud_prog_report(cpr, io_unit)
        implicit none
        type(aeros_cloud_prog_class), intent(in) :: cpr
        integer,                      intent(in) :: io_unit
        write(io_unit, '(a)')      "  cloud (prognostic):"
        write(io_unit, '(a,l1)')   "    enabled  = ", cpr%enabled
        write(io_unit, '(a,f6.3,a,f6.3)') "    rhc_sfc/top = ", cpr%rhc_sfc, " / ", cpr%rhc_top
        write(io_unit, '(a,f8.1,a,f8.1)') "    tau_form/evap [s] = ", cpr%tau_form, " / ", cpr%tau_evap
        return
    end subroutine aeros_cloud_prog_report

end module aeros_cloud_prog
