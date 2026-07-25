module aeros_condensation
    ! Large-scale (stratiform) condensation: the first moist physics.
    !
    ! Where the air is supersaturated, the excess vapour condenses, the latent
    ! heat warms the air, and the condensate falls out as precipitation. That is
    ! the whole of it at M2.3b -- no cloud water is stored, nothing re-evaporates,
    ! there is no convection (the next commit) and no ice phase (condensate is
    ! liquid, latent heat L_v, everywhere). It is the minimum that closes a moist
    ! energy and water budget, and it is built to be exactly that budget rather
    ! than more.
    !
    ! === The two seams it has to touch ======================================
    !
    ! Condensation couples humidity, which is a GRIDPOINT prognostic, to
    ! temperature, which is a SPECTRAL one. So it acts at two places at once, and
    ! keeping them consistent is the point:
    !
    !   DRYING is applied straight to the gridpoint humidity (qv_g): the vapour
    !   that condenses is removed from the cell it condensed in. No transform.
    !
    !   HEATING is added to the temperature TENDENCY on the grid (wrk%dtdt),
    !   exactly where the Held-Suarez forcing is added, so it rides the same
    !   transform back to spectral space and the same semi-implicit leapfrog as
    !   the dynamical heating. It is a rate, L_v/cp_d * dq_c / dt, so that one
    !   step's worth of it warms the air by L_v/cp_d * dq_c -- the latent heat of
    !   exactly the water that left.
    !
    ! Because both come from the SAME condensed amount dq_c, the column moist
    ! static energy is conserved by construction: the dry static energy gains
    ! L_v dq_c and the latent energy loses it. The water budget closes the same
    ! way: the vapour lost equals the precipitation produced. tests/
    ! test_condensation checks both as equalities, not tolerances.
    !
    ! The two seams are on different time discretizations -- the drying is a
    ! forward step on qv_g, the heating goes through the 2 dt leapfrog and its
    ! filter -- so over a run they track to the filter's accuracy rather than
    ! exactly. That is the same small inconsistency every leapfrog spectral model
    ! carries between its grid physics and its spectral dynamics; it is measured
    ! at the integration level, not asserted away here.
    !
    ! === The saturation adjustment ==========================================
    !
    ! Condensing changes the temperature, which changes the saturation humidity
    ! it is condensing toward, so the amount is the root of an implicit
    ! equation. One Newton step linearizes it:
    !
    !   dq_c = (q - q_sat(T)) / (1 + (L_v/cp_d) dq_sat/dT)
    !
    ! and a second iteration removes essentially all of the linearization error.
    ! Only supersaturated cells condense (dq_c >= 0); nothing evaporates, because
    ! there is no condensate to evaporate yet.

    use aeros_defs,     only : dp, wp, io_unit_err, R_d, R_v, cp_d, grav, L_v, &
                                aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure

    use nml,            only : nml_read

    implicit none

    private

    ! Tetens (1930) saturation vapour pressure over liquid water [Pa]:
    !   e_s(T) = ES0 exp( ATET (T - T_TET) / (T - BTET) )
    ! with the derivative de_s/dT = e_s ATET (T_TET - BTET) / (T - BTET)^2.
    real(wp), parameter :: ES0   = 610.78_wp
    real(wp), parameter :: ATET  = 17.27_wp
    real(wp), parameter :: T_TET = 273.16_wp
    real(wp), parameter :: BTET  = 35.86_wp
    real(wp), parameter :: EPS   = R_d/R_v          ! ~0.622

    type aeros_cond_class
        ! Configuration and the last step's precipitation field.

        logical :: enabled = .FALSE.

        ! Critical relative humidity for condensation. 1.0 is true saturation
        ! adjustment; a lower value (large-scale condensation before gridbox-mean
        ! saturation, standing in for sub-grid variability) is left as a knob but
        ! defaults to 1.0 -- the honest, parameter-free choice for a first cut.
        real(wp) :: rh_crit = 1.0_wp

        integer :: nlon = 0, nlat = 0

        ! Precipitation rate from the last apply, [kg m-2 s-1] = mm s-1 of
        ! water. A diagnostic, not a prognostic: it leaves the atmosphere the
        ! instant it forms.
        real(wp), allocatable :: precip(:,:)
    end type aeros_cond_class

    public :: aeros_cond_class
    public :: aeros_condensation_init
    public :: aeros_condensation_load
    public :: aeros_condensation_end
    public :: aeros_condensation_apply
    public :: aeros_condensation_report
    public :: aeros_qsat            ! exposed for the test

contains

    subroutine aeros_condensation_init(cnd, grd, enabled)

        implicit none

        type(aeros_cond_class), intent(inout) :: cnd
        type(aeros_grid_class), intent(in)    :: grd
        logical, intent(in) :: enabled

        call aeros_condensation_end(cnd)

        cnd%enabled = enabled
        cnd%rh_crit = 1.0_wp
        cnd%nlon    = grd%nlon
        cnd%nlat    = grd%nlat

        allocate(cnd%precip(grd%nlon, grd%nlat))
        cnd%precip = 0.0_wp

        return

    end subroutine aeros_condensation_init

    subroutine aeros_condensation_load(cnd, filename, grd)
        ! Configure from the `aeros_moisture` namelist group.

        implicit none

        type(aeros_cond_class), intent(inout) :: cnd
        character(len=*),       intent(in)    :: filename
        type(aeros_grid_class), intent(in)    :: grd

        logical  :: moist
        real(wp) :: rh_crit

        moist   = .FALSE.
        rh_crit = 1.0_wp
        call nml_read(filename, "aeros_moisture", "moist",   moist)
        call nml_read(filename, "aeros_moisture", "rh_crit", rh_crit)

        call aeros_condensation_init(cnd, grd, moist)

        if (rh_crit <= 0.0_wp .or. rh_crit > 1.0_wp) then
            write(io_unit_err,*) "aeros_condensation_load:: error: rh_crit must be in (0,1], got ", &
                                    rh_crit
            error stop 1
        end if
        cnd%rh_crit = rh_crit

        return

    end subroutine aeros_condensation_load

    subroutine aeros_condensation_end(cnd)

        implicit none

        type(aeros_cond_class), intent(inout) :: cnd

        if (allocated(cnd%precip)) deallocate(cnd%precip)
        cnd%enabled = .FALSE.
        cnd%nlon = 0; cnd%nlat = 0

        return

    end subroutine aeros_condensation_end

    subroutine aeros_condensation_apply(cnd, vg, t_g, qv_g, lnps_g, dtdt, dt)
        ! Condense the supersaturation: dry qv_g, add the latent heating to the
        ! temperature tendency dtdt, and record the column precipitation.
        !
        ! t_g and lnps_g are the current gridpoint temperature and log surface
        ! pressure (aeros_tendency's wrk); qv_g is the gridpoint humidity, dried
        ! in place; dtdt is the grid-space temperature tendency the heating is
        ! added to, before it is transformed with the dynamics.

        implicit none

        type(aeros_cond_class),  intent(inout) :: cnd
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(inout) :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat)      ln[Pa]
        real(wp), intent(inout) :: dtdt(:,:,:)    ! (nlon,nlat,nlev) [K s-1]
        real(wp), intent(in)    :: dt             ! [s]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: qs, dqsdt, gam, dqc, pcol, hcp
        integer  :: i, j, k, it

        if (.not. cnd%enabled) return

        hcp = real(L_v, wp)/real(cp_d, wp)

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,it,phalf,pfull,dpc,qs,dqsdt,gam,dqc,pcol)
        do j = 1, cnd%nlat
            do i = 1, cnd%nlon
                call aeros_vgrid_pressure(vg, exp(lnps_g(i,j)), phalf, pfull, dpc)
                pcol = 0.0_wp

                do k = 1, vg%nlev
                    ! Newton iterations of the saturation adjustment, each
                    ! against the temperature the accumulated heating implies.
                    ! Three, not two: per step the supersaturation is tiny
                    ! (radiative/adiabatic cooling is a fraction of a kelvin) and
                    ! one iteration would do, but the extra pair is cheap
                    ! insurance for the larger excursions of a spin-up, where the
                    ! quadratic convergence still lands it at saturation.
                    dqc = 0.0_wp
                    do it = 1, 3
                        call aeros_qsat(t_g(i,j,k) + hcp*dqc, pfull(k), qs, dqsdt)
                        gam = hcp*dqsdt
                        dqc = dqc + (qv_g(i,j,k) - dqc - cnd%rh_crit*qs)/(1.0_wp + gam)
                    end do

                    ! Only condensation, and never more vapour than is present.
                    dqc = max(0.0_wp, min(dqc, qv_g(i,j,k)))

                    if (dqc > 0.0_wp) then
                        qv_g(i,j,k) = qv_g(i,j,k) - dqc
                        dtdt(i,j,k) = dtdt(i,j,k) + hcp*dqc/dt
                        pcol        = pcol + dqc*dpc(k)
                    end if
                end do

                ! Column condensate falls out immediately: mass per area per
                ! time. dp/g is the layer mass per area.
                cnd%precip(i,j) = pcol/(real(grav, wp)*dt)
            end do
        end do
        !$omp end parallel do

        return

    end subroutine aeros_condensation_apply

    subroutine aeros_qsat(t, p, qs, dqsdt)
        ! Saturation specific humidity [kg kg-1] and its temperature derivative
        ! [kg kg-1 K-1], over liquid water (Tetens).
        !
        !   q_sat = eps e_s / (p - (1 - eps) e_s)
        !
        ! The (1-eps) e_s term in the denominator matters at low pressure and is
        ! kept; dropping it (q_sat ~ eps e_s / p) is the common shortcut and is
        ! wrong by several percent in the upper troposphere.

        implicit none

        real(wp), intent(in)  :: t, p
        real(wp), intent(out) :: qs, dqsdt

        real(wp) :: es, desdt, den

        es    = ES0*exp(ATET*(t - T_TET)/(t - BTET))
        desdt = es*ATET*(T_TET - BTET)/(t - BTET)**2

        den   = p - (1.0_wp - EPS)*es
        ! Guard the (very high, very warm) case where saturation vapour pressure
        ! approaches the total pressure and the denominator collapses; there the
        ! air is effectively all vapour and q_sat is meaningless, so cap it.
        if (den <= 0.1_wp*p) then
            qs    = 1.0_wp
            dqsdt = 0.0_wp
            return
        end if

        qs    = EPS*es/den
        dqsdt = EPS*p*desdt/den**2

        return

    end subroutine aeros_qsat

    subroutine aeros_condensation_report(cnd, io_unit)

        implicit none

        type(aeros_cond_class), intent(in) :: cnd
        integer, intent(in), optional :: io_unit

        integer :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        write(iou,*) ""
        write(iou,"(a)") " == condensation =="
        if (cnd%enabled) then
            write(iou,"(a)")        "   large-scale condensation    ON"
            write(iou,"(a,f9.3)")   "   critical relative humidity ", cnd%rh_crit
        else
            write(iou,"(a)")        "   large-scale condensation    off (dry)"
        end if

        return

    end subroutine aeros_condensation_report

end module aeros_condensation
