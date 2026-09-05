module aeros_condensation
    ! Large-scale (stratiform) condensation: the first moist physics.
    !
    ! Where the air is supersaturated, the excess vapour condenses, the latent
    ! heat warms the air, and the condensate falls out as precipitation. No cloud
    ! water is stored. Saturation is MIXED-PHASE (aeros_qsat: over liquid at and
    ! above 0 C, over ice at and below -23 C, IFS/ERA5's quadratic blend between)
    ! and the latent heat follows the same blend (aeros_latent: L_v -> L_s), so the
    ! cold upper troposphere condenses toward ice saturation and releases the heat
    ! of sublimation. Condensation removes vapour down to rh_crit*q_sat (default
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
    !   HEATING is added to the grid increment wrk%dt_phys [K], the same path
    !   convection, radiation and the surface fluxes use; aeros_timestep turns
    !   it into a rate and couples it into the thermodynamic tendency ahead of
    !   the semi-implicit solve (aeros_timestep header, "Physics coupling"). It
    !   is an increment L_v/cp_d * dq_c -- the latent heat of exactly the water
    !   that left.
    !
    ! Because both come from the SAME condensed amount dq_c, the column moist
    ! static energy is conserved by construction: the dry static energy gains
    ! L_v dq_c and the latent energy loses it. The water budget closes the same
    ! way: the vapour lost equals the precipitation produced. tests/
    ! test_condensation checks both as equalities, not tolerances.
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
    !
    ! === The time scale =====================================================
    !
    ! Only 1/time_scale of the adjustment (condensation AND reevaporation) is
    ! applied per step -- SpeedyWeather's ImplicitCondensation `time_scale`,
    ! default 3 steps. The heating is coupled into the leapfrog from time level
    ! n and lands on n+1: an operator that removes ALL of the supersaturation
    ! in one step is a stiff restoring term diagnosed on one leapfrog chain and
    ! applied to the other, and it drives the 2 dt computational mode
    ! (aeros_timestep header). Relaxing over a few steps keeps the per-step
    ! kick small enough for the Robert filter to hold; the equilibrium state is
    ! the same. time_scale = 1 recovers the exact one-step adjustment (what
    ! tests/test_condensation checks the budgets on).

    use aeros_defs,     only : dp, wp, io_unit_err, R_d, R_v, cp_d, grav, L_v, L_s, &
                                aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure

    use nml,            only : nml_read

    implicit none

    private

    ! Tetens (1930) saturation vapour pressure over liquid water [Pa]:
    !   e_s(T) = ES0 exp( ATET (T - T_TET) / (T - BTET) )
    ! with the derivative de_s/dT = e_s ATET (T_TET - BTET) / (T - BTET)^2;
    ! over ice the IFS Tetens form e_i = ESI0 exp( AICE (T - T_TET) / (T - BICE) ).
    ! Phase blend (IFS/ERA5 convention): liquid fraction alpha = 1 at T >= T_TET,
    ! 0 at T <= T_ICEP, alpha = ((T - T_ICEP)/(T_TET - T_ICEP))^2 between.
    real(wp), parameter :: ES0   = 610.78_wp
    real(wp), parameter :: ATET  = 17.27_wp
    real(wp), parameter :: T_TET = 273.16_wp
    real(wp), parameter :: BTET  = 35.86_wp
    real(wp), parameter :: ESI0  = 611.21_wp
    real(wp), parameter :: AICE  = 22.587_wp
    real(wp), parameter :: BICE  = -0.7_wp
    real(wp), parameter :: T_ICEP = 250.16_wp       ! all ice at and below (-23 C)
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

        ! Fraction 1/time_scale of the saturation adjustment (and of the
        ! reevaporation) applied per step, in steps (SpeedyWeather's
        ! ImplicitCondensation time_scale, default 3). 1 = exact one-step
        ! adjustment. See the module header, "The time scale".
        real(wp) :: time_scale = 3.0_wp

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
    public :: aeros_qsat            ! mixed-phase saturation (shared by all physics)
    public :: aeros_latent          ! matching latent heat L(T) = alpha L_v + (1-alpha) L_s
    public :: aeros_liquid_fraction ! the phase blend alpha(T)

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
        cnd%time_scale = 3.0_wp
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
        real(wp) :: rh_crit, reevap, time_scale

        moist      = .FALSE.
        rh_crit    = 0.95_wp
        reevap     = 30.0_wp
        time_scale = 3.0_wp
        call nml_read(filename, "aeros_moisture", "moist",   moist, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "rh_crit", rh_crit, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "reevap",  reevap, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "cond_time_scale", time_scale, &
                      defaults_file=defaults_file)

        call aeros_condensation_init(cnd, grd, moist)

        if (time_scale < 1.0_wp) then
            write(io_unit_err,*) "aeros_condensation_load:: error: cond_time_scale must be >= 1 "// &
                                    "step, got ", time_scale
            error stop 1
        end if
        cnd%time_scale = time_scale

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

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,it,phalf,pfull,dpc,qs,dqsdt,gam,dqc,pfall,hcp,qthr,frac,dqe)
        do j = 1, cnd%nlat
            do i = 1, cnd%nlon
                call aeros_vgrid_pressure(vg, exp(lnps_g(i,j)), phalf, pfull, dpc)

                ! Sweep top-to-bottom carrying the falling precipitation flux
                ! pfall [kg kg-1 Pa] (= g x mass/area): each layer first
                ! reevaporates some of the rain falling from above, then condenses
                ! its own supersaturation into it.
                pfall = 0.0_wp

                do k = 1, vg%nlev
                    ! latent heat of the layer's phase (L_v warm, L_s cold), per K of
                    ! humidity change; held at the pre-adjustment temperature so the
                    ! heat equals L times the water that leaves, exactly
                    hcp = aeros_latent(t_g(i,j,k))/real(cp_d, wp)
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
                        dqe  = dqe/cnd%time_scale                ! relaxed, "The time scale"
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

                    ! Only condensation, and never more vapour than is present;
                    ! 1/time_scale of the adjustment per step ("The time scale").
                    dqc = max(0.0_wp, min(dqc, qv_g(i,j,k)))/cnd%time_scale

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

    pure real(wp) function aeros_liquid_fraction(t) result(alpha)
        ! Liquid fraction of the saturation blend: 1 at/above 0 C, 0 at/below
        ! -23 C, quadratic between (IFS/ERA5).
        implicit none
        real(wp), intent(in) :: t
        real(wp) :: x
        if (t >= T_TET) then
            alpha = 1.0_wp
        else if (t <= T_ICEP) then
            alpha = 0.0_wp
        else
            x = (t - T_ICEP)/(T_TET - T_ICEP)
            alpha = x*x
        end if
        return
    end function aeros_liquid_fraction

    pure real(wp) function aeros_latent(t) result(lat)
        ! Latent heat [J kg-1] matching aeros_qsat's phase blend: vaporization
        ! for the liquid part, sublimation for the ice part.
        implicit none
        real(wp), intent(in) :: t
        real(wp) :: alpha
        alpha = aeros_liquid_fraction(t)
        lat = alpha*real(L_v, wp) + (1.0_wp - alpha)*real(L_s, wp)
        return
    end function aeros_latent

    subroutine aeros_qsat(t, p, qs, dqsdt)
        ! Saturation specific humidity [kg kg-1] and its temperature derivative
        ! [kg kg-1 K-1], mixed-phase: e_s = alpha e_liq + (1 - alpha) e_ice with
        ! the blend alpha(T) of aeros_liquid_fraction (pure liquid Tetens at and
        ! above 0 C, so warm results are unchanged).
        !
        !   q_sat = eps e_s / (p - (1 - eps) e_s)
        !
        ! The (1-eps) e_s term in the denominator matters at low pressure and is
        ! kept; dropping it (q_sat ~ eps e_s / p) is the common shortcut and is
        ! wrong by several percent in the upper troposphere.

        implicit none

        real(wp), intent(in)  :: t, p
        real(wp), intent(out) :: qs, dqsdt

        real(wp) :: es, desdt, den, el, deldt, ei, deidt, alpha, daldt

        el    = ES0*exp(ATET*(t - T_TET)/(t - BTET))
        deldt = el*ATET*(T_TET - BTET)/(t - BTET)**2
        alpha = aeros_liquid_fraction(t)
        if (alpha >= 1.0_wp) then
            es = el; desdt = deldt
        else
            ei    = ESI0*exp(AICE*(t - T_TET)/(t - BICE))
            deidt = ei*AICE*(T_TET - BICE)/(t - BICE)**2
            daldt = 0.0_wp
            if (alpha > 0.0_wp) daldt = 2.0_wp*(t - T_ICEP)/(T_TET - T_ICEP)**2
            es    = alpha*el + (1.0_wp - alpha)*ei
            desdt = alpha*deldt + (1.0_wp - alpha)*deidt + daldt*(el - ei)
        end if

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
            write(iou,"(a,f9.3)")   "   adjustment time scale [steps]", cnd%time_scale
        else
            write(iou,"(a)")        "   large-scale condensation    off (dry)"
        end if

        return

    end subroutine aeros_condensation_report

end module aeros_condensation
