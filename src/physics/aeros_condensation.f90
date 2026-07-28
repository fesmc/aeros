module aeros_condensation
    ! Large-scale (stratiform) condensation: the first moist physics.
    !
    ! Where the air is supersaturated, the excess vapour condenses, the latent
    ! heat warms the air, and the condensate falls out as precipitation. No cloud
    ! water is stored and there is no ice phase (condensate is liquid, latent heat
    ! L_v, everywhere). Condensation removes vapour down to rh_crit*q_sat (default
    ! 0.95, not full saturation), and the falling precipitation RE-EVAPORATES into
    ! sub-saturated layers below (moistening + latent cooling) before what remains
    ! reaches the surface -- both following SpeedyWeather's ImplicitCondensation,
    ! and both drying sinks the column would otherwise lack (see the rh_crit and
    ! reevap notes below).
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
    !   HEATING is added to the forward-split increment on the grid (wrk%dt_phys,
    !   [K]), the same path convection, radiation and the surface fluxes use: it
    !   is applied to the n+1 state after the dynamics step, forward in time,
    !   decoupled from the centered leapfrog. It is an increment L_v/cp_d * dq_c
    !   -- one step's warming from the latent heat of exactly the water that left.
    !
    !   This was originally on the centered dtdt path, on the argument that
    !   condensation's heating is small and smooth. In the coupled RCE that
    !   assumption breaks: at the hot subtropical latitude the heating reaches
    !   ~100 K/day in the lowest layer and, on the centered leapfrog, excites the
    !   computational mode -- the same failure convection had before it was moved
    !   forward-split (M2.3e/section 12.1). So the heating now rides the same
    !   forward path, symmetric with the forward drying of qv_g. See
    !   docs/m2_results.md and m2_handoff.md (RCE level-12 instability).
    !
    ! Because both come from the SAME condensed amount dq_c, the column moist
    ! static energy is conserved by construction: the dry static energy gains
    ! L_v dq_c and the latent energy loses it. The water budget closes the same
    ! way: the vapour lost equals the precipitation produced. tests/
    ! test_condensation checks both as equalities, not tolerances. Drying and
    ! heating now share the SAME forward discretization, so they track exactly
    ! rather than to the filter's accuracy as when the heating was centered.
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

        ! Critical relative humidity for condensation: condensation removes vapour
        ! down to rh_crit*q_sat, not full saturation, standing in for the saturated
        ! fraction of a partly-cloudy gridbox. Default 0.95 (SpeedyWeather's
        ! ImplicitCondensation threshold); 1.0 recovers exact saturation
        ! adjustment. Pinning the grid mean at 100% (rh_crit=1) is a major cause of
        ! the overcast moist bias, so the sub-saturating default is deliberate.
        real(wp) :: rh_crit = 0.95_wp

        ! Reevaporation efficiency [1/(kg/kg)] of falling precipitation into
        ! sub-saturated layers below (SpeedyWeather's `reevaporation`): the
        ! fraction of the falling flux evaporated in a layer is
        ! min(reevap*(q_sat - q), 1), moistening the layer toward rh_crit*q_sat and
        ! cooling it by the latent heat. 0 disables it (all condensate falls
        ! straight to the surface).
        real(wp) :: reevap = 30.0_wp

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
        cnd%rh_crit = 0.95_wp
        cnd%reevap  = 30.0_wp
        cnd%nlon    = grd%nlon
        cnd%nlat    = grd%nlat

        allocate(cnd%precip(grd%nlon, grd%nlat))
        cnd%precip = 0.0_wp

        return

    end subroutine aeros_condensation_init

    subroutine aeros_condensation_load(cnd, filename, grd, defaults_file)
        ! Configure from the `aeros_moisture` namelist group.

        implicit none

        type(aeros_cond_class), intent(inout) :: cnd
        character(len=*),       intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_grid_class), intent(in)    :: grd

        logical  :: moist
        real(wp) :: rh_crit, reevap

        moist   = .FALSE.
        rh_crit = 0.95_wp
        reevap  = 30.0_wp
        call nml_read(filename, "aeros_moisture", "moist",   moist, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "rh_crit", rh_crit, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "reevap",  reevap, defaults_file=defaults_file)

        call aeros_condensation_init(cnd, grd, moist)

        if (rh_crit <= 0.0_wp .or. rh_crit > 1.0_wp) then
            write(io_unit_err,*) "aeros_condensation_load:: error: rh_crit must be in (0,1], got ", &
                                    rh_crit
            error stop 1
        end if
        if (reevap < 0.0_wp) then
            write(io_unit_err,*) "aeros_condensation_load:: error: reevap must be >= 0, got ", reevap
            error stop 1
        end if
        cnd%rh_crit = rh_crit
        cnd%reevap  = reevap

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

    subroutine aeros_condensation_apply(cnd, vg, t_g, qv_g, lnps_g, dt_phys, dt)
        ! Condense the supersaturation: dry qv_g, add the latent heating to the
        ! forward-split increment dt_phys, and record the column precipitation.
        !
        ! t_g and lnps_g are the current gridpoint temperature and log surface
        ! pressure (aeros_tendency's wrk); qv_g is the gridpoint humidity, dried
        ! in place; dt_phys is the forward-split grid temperature increment [K]
        ! the heating is added to, applied to the n+1 state after the dynamics
        ! step (NOT the centered leapfrog -- see the module header).

        implicit none

        type(aeros_cond_class),  intent(inout) :: cnd
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(inout) :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat)      ln[Pa]
        real(wp), intent(inout) :: dt_phys(:,:,:) ! (nlon,nlat,nlev) [K] increment
        real(wp), intent(in)    :: dt             ! [s]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: qs, dqsdt, gam, dqc, pfall, hcp, qthr, frac, dqe
        integer  :: i, j, k, it

        if (.not. cnd%enabled) return

        hcp = real(L_v, wp)/real(cp_d, wp)

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,it,phalf,pfull,dpc,qs,dqsdt,gam,dqc,pfall,qthr,frac,dqe)
        do j = 1, cnd%nlat
            do i = 1, cnd%nlon
                call aeros_vgrid_pressure(vg, exp(lnps_g(i,j)), phalf, pfull, dpc)

                ! Sweep top-to-bottom carrying the falling precipitation flux
                ! pfall [kg kg-1 Pa] (= g x mass/area): each layer first
                ! reevaporates some of the rain falling from above, then condenses
                ! its own supersaturation into it.
                pfall = 0.0_wp

                do k = 1, vg%nlev
                    call aeros_qsat(t_g(i,j,k), pfull(k), qs, dqsdt)
                    qthr = cnd%rh_crit*qs

                    ! Reevaporation: rain falling into a layer below the
                    ! condensation threshold evaporates, moistening it toward the
                    ! threshold and cooling it by the latent heat. The evaporated
                    ! fraction of the falling flux scales with the dryness
                    ! (SpeedyWeather's min(reevap*(q_sat - q), 1)), capped so it
                    ! never overshoots the threshold (which would just re-condense).
                    if (pfall > 0.0_wp .and. cnd%reevap > 0.0_wp .and. &
                        qv_g(i,j,k) < qthr) then
                        frac = min(cnd%reevap*(qs - qv_g(i,j,k)), 1.0_wp)
                        dqe  = frac*pfall/dpc(k)                 ! [kg kg-1]
                        dqe  = min(dqe, qthr - qv_g(i,j,k))
                        if (dqe > 0.0_wp) then
                            qv_g(i,j,k)    = qv_g(i,j,k) + dqe
                            dt_phys(i,j,k) = dt_phys(i,j,k) - hcp*dqe  ! latent cooling
                            pfall          = pfall - dqe*dpc(k)
                        end if
                    end if

                    ! Newton iterations of the saturation adjustment, each
                    ! against the temperature the accumulated heating implies.
                    ! Three, not two: per step the supersaturation is tiny
                    ! (radiative/adiabatic cooling is a fraction of a kelvin) and
                    ! one iteration would do, but the extra pair is cheap
                    ! insurance for the larger excursions of a spin-up, where the
                    ! quadratic convergence still lands it at the threshold.
                    dqc = 0.0_wp
                    do it = 1, 3
                        call aeros_qsat(t_g(i,j,k) + hcp*dqc, pfull(k), qs, dqsdt)
                        gam = hcp*dqsdt
                        dqc = dqc + (qv_g(i,j,k) - dqc - cnd%rh_crit*qs)/(1.0_wp + gam)
                    end do

                    ! Only condensation, and never more vapour than is present.
                    dqc = max(0.0_wp, min(dqc, qv_g(i,j,k)))

                    if (dqc > 0.0_wp) then
                        qv_g(i,j,k)    = qv_g(i,j,k) - dqc
                        dt_phys(i,j,k) = dt_phys(i,j,k) + hcp*dqc   ! [K] increment
                        pfall          = pfall + dqc*dpc(k)
                    end if
                end do

                ! Whatever rain survives to the surface is the precipitation:
                ! mass per area per time. dp/g is the layer mass per area.
                cnd%precip(i,j) = pfall/(real(grav, wp)*dt)
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
            write(iou,"(a,f9.3)")   "   reevaporation efficiency   ", cnd%reevap
        else
            write(iou,"(a)")        "   large-scale condensation    off (dry)"
        end if

        return

    end subroutine aeros_condensation_report

end module aeros_condensation
