module aeros_convection
    ! Moist convective adjustment: the physics that keeps the tropics from
    ! saturating at the grid scale.
    !
    ! Large-scale condensation (aeros_condensation) only removes supersaturation
    ! from the gridbox mean. Left alone, a moist column heated from below never
    ! overturns -- it just saturates layer by layer and rains stratiform
    ! everywhere, which is not how the tropics work. Convection is the fast
    ! vertical redistribution that a resolved column cannot do on its own at this
    ! resolution: it moves heat and moisture up, condenses, precipitates, and
    ! leaves the column neutrally stratified.
    !
    ! === Manabe (1965) moist convective adjustment ==========================
    !
    ! The scheme sweeps the column for unstable adjacent layer pairs and relaxes
    ! each to neutrality, iterating until the whole column is stable:
    !
    !   DRY pairs (unsaturated, potential temperature decreasing upward) mix to
    !   a common potential temperature, conserving enthalpy int c_p T dp. No
    !   water, no precipitation -- just heat moved to where buoyancy wants it.
    !
    !   MOIST pairs (both saturated, saturated moist static energy decreasing
    !   upward) relax to a moist adiabat -- equal saturated MSE
    !   h* = c_p T + Phi + L q_sat(T) between the two -- while conserving the
    !   column moist static energy int (c_p T + L q) dp and holding both at
    !   saturation. The water that leaves as the layers re-saturate at their new
    !   temperatures is the precipitation.
    !
    ! Conservation is the property to hold onto and the one the test pins: the
    ! dry adjustment conserves int c_p T dp exactly, the moist adjustment
    ! conserves int (c_p T + L q) dp exactly, and precipitation equals the water
    ! removed. The moist relaxation of a pair is a 2x2 nonlinear system (the two
    ! new temperatures, both layers saturated, MSE conserved, moist-neutral),
    ! solved by Newton against dq_sat/dT.
    !
    ! === Manabe vs Simplified Betts-Miller ==================================
    !
    ! Manabe conserves by construction -- it redistributes rather than relaxing
    ! toward an external reference, so there is no timescale and no reference
    ! humidity -- but its pairwise sweep propagates instability one layer per
    ! pass, so a deep column converges only over tens of passes. It is kept as
    ! the parameter-free reference scheme (conv_scheme = "manabe").
    !
    ! The default, and the scheme built for a stable-but-fast running model, is
    ! Simplified Betts-Miller (Frierson 2007; conv_scheme = "sbm"): see
    ! sbm_adjust below. One Newton solve per level, no iteration to neutrality,
    ! and it relaxes toward the reference over a finite timescale rather than
    ! adjusting instantaneously -- gentler on the integrator and tunable once
    ! there is radiation and observations. The seam carries both: aeros_conv_class
    ! has a `scheme` selector, aeros_convection_apply dispatches per column to
    ! adjust_column, and adjust_column branches on the scheme. The apply loop,
    ! the tendency plumbing and the precipitation accounting are shared.
    !
    ! === Coupling: forward-split, NOT through the centered leapfrog =========
    !
    ! Convection changes temperature (up AND down -- it is a redistribution) and
    ! humidity. It acts at the grid seam, but unlike condensation it does NOT
    ! feed its heating through wrk%dtdt and the centered leapfrog: convective
    ! heating is large and sign-alternating in the vertical, and a term
    ! evaluated at time n and applied through the centered 2 dt step excites the
    ! leapfrog computational mode (the run NaNs in tens of steps). Instead the
    ! temperature change is written as a forward INCREMENT into wrk%dt_phys and
    ! applied to the n+1 state after the dynamics step (aeros_timestep_step),
    ! decoupled from the centered core -- the same forward-in-time treatment the
    ! gridpoint humidity already gets. The humidity change is applied to qv_g,
    ! and both come from the same adjustment so the column budget closes. It runs
    ! BEFORE large-scale condensation, which then mops up whatever
    ! gridbox-mean supersaturation the convection did not remove.

    use aeros_defs,       only : dp, wp, io_unit_err, R_d, cp_d, grav, L_v, &
                                    kappa, p0, aeros_grid_class
    use aeros_vertical,   only : aeros_vgrid_class, aeros_vgrid_pressure
    use aeros_condensation, only : aeros_qsat

    use nml,              only : nml_read

    implicit none

    private

    integer, parameter :: SCHEME_MANABE = 1
    integer, parameter :: SCHEME_SBM    = 2

    ! Column sweeps for the pairwise relaxation to converge.
    !
    ! Manabe's pairwise adjustment propagates instability one layer per sweep,
    ! so a deep, strongly unstable column converges only geometrically -- at
    ! T21L20 the worst-case cold-start residual falls from ~50 J/kg at 60 sweeps
    ! to ~15 J/kg (0.015 K) at 80 to machine-neutral near 400. 80 leaves a column
    ! physically neutral (~0.01 K) at a bounded cost. It is not a concern for a
    ! running model, where each step's instability is a small increment on an
    ! already-neutralized column and a handful of sweeps suffice; only a cold
    ! start pays the full price, once, and even then a slightly-residual column
    ! simply finishes overturning over the next few steps. Only unstable columns
    ! iterate -- a stable one exits after the first sweep.
    !
    ! The proper fix for deep convection is a MULTI-LAYER adjustment: neutralize
    ! a whole connected unstable saturated segment at once (h* = const across it,
    ! conserving the segment MSE, which reduces to a single reference value plus
    ! a 1-D solve per layer) rather than one pair at a time. It converges in a
    ! few passes instead of tens. It is the natural next refinement; the pairwise
    ! form here is the canonical 1965 algorithm and is correct, just slower to
    ! converge on the cold-start extreme.
    integer, parameter :: MAXSWEEP = 80

    type aeros_conv_class
        logical :: enabled = .FALSE.
        integer :: scheme  = SCHEME_SBM

        ! Simplified Betts-Miller knobs (unused by Manabe). tau is the
        ! convective relaxation timescale; rh_ref the reference relative
        ! humidity the moist adiabat is dried toward. Frierson (2007) values.
        real(wp) :: tau    = 7200.0_wp   ! [s] ~ 2 h
        real(wp) :: rh_ref = 0.7_wp      ! [-]

        integer :: nlon = 0, nlat = 0

        ! Precipitation rate from the last apply, [kg m-2 s-1].
        real(wp), allocatable :: precip(:,:)
    end type aeros_conv_class

    public :: aeros_conv_class
    public :: SCHEME_MANABE, SCHEME_SBM     ! exposed for the test to select
    public :: aeros_convection_init
    public :: aeros_convection_load
    public :: aeros_convection_end
    public :: aeros_convection_apply
    public :: aeros_convection_report

contains

    subroutine aeros_convection_init(cnv, grd, enabled)

        implicit none

        type(aeros_conv_class), intent(inout) :: cnv
        type(aeros_grid_class), intent(in)    :: grd
        logical, intent(in) :: enabled

        call aeros_convection_end(cnv)

        cnv%enabled = enabled
        cnv%scheme  = SCHEME_SBM
        cnv%tau     = 7200.0_wp
        cnv%rh_ref  = 0.7_wp
        cnv%nlon    = grd%nlon
        cnv%nlat    = grd%nlat

        allocate(cnv%precip(grd%nlon, grd%nlat))
        cnv%precip = 0.0_wp

        return

    end subroutine aeros_convection_init

    subroutine aeros_convection_load(cnv, filename, grd, defaults_file)
        ! Configure from the `aeros_moisture` namelist group. `convect` turns it
        ! on; `conv_scheme` selects ("sbm" or "manabe"); `conv_tau` and
        ! `conv_rhref` are the Simplified Betts-Miller knobs. An unknown scheme
        ! name is an error, not a silent fallback, so a typo does not quietly
        ! run a scheme the user did not ask for.

        implicit none

        type(aeros_conv_class), intent(inout) :: cnv
        character(len=*),       intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_grid_class), intent(in)    :: grd

        logical            :: convect
        character(len=32)  :: conv_scheme
        real(wp)           :: conv_tau, conv_rhref

        convect     = .FALSE.
        conv_scheme = "sbm"
        conv_tau    = 7200.0_wp
        conv_rhref  = 0.7_wp
        call nml_read(filename, "aeros_moisture", "convect",     convect, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "conv_scheme", conv_scheme, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "conv_tau",    conv_tau, defaults_file=defaults_file)
        call nml_read(filename, "aeros_moisture", "conv_rhref",  conv_rhref, defaults_file=defaults_file)

        call aeros_convection_init(cnv, grd, convect)

        select case (trim(conv_scheme))
        case ("sbm")
            cnv%scheme = SCHEME_SBM
        case ("manabe")
            cnv%scheme = SCHEME_MANABE
        case default
            write(io_unit_err,*) "aeros_convection_load:: error: unknown conv_scheme '"// &
                                    trim(conv_scheme)//"' (expected 'sbm' or 'manabe')"
            error stop 1
        end select

        if (conv_tau <= 0.0_wp) then
            write(io_unit_err,*) "aeros_convection_load:: error: conv_tau must be > 0, got ", conv_tau
            error stop 1
        end if
        if (conv_rhref <= 0.0_wp .or. conv_rhref > 1.0_wp) then
            write(io_unit_err,*) "aeros_convection_load:: error: conv_rhref must be in (0,1], got ", conv_rhref
            error stop 1
        end if
        cnv%tau    = conv_tau
        cnv%rh_ref = conv_rhref

        return

    end subroutine aeros_convection_load

    subroutine aeros_convection_end(cnv)

        implicit none

        type(aeros_conv_class), intent(inout) :: cnv

        if (allocated(cnv%precip)) deallocate(cnv%precip)
        cnv%enabled = .FALSE.
        cnv%nlon = 0; cnv%nlat = 0

        return

    end subroutine aeros_convection_end

    subroutine aeros_convection_apply(cnv, vg, t_g, qv_g, lnps_g, dt_phys, dt)
        ! Adjust each column for moist (and dry) convective instability: dry
        ! qv_g, add the temperature change to the forward-split physics increment
        ! dt_phys (NOT the centered-leapfrog rate dtdt), record precip.

        implicit none

        type(aeros_conv_class),  intent(inout) :: cnv
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(inout) :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat)      ln[Pa]
        real(wp), intent(inout) :: dt_phys(:,:,:) ! (nlon,nlat,nlev) [K]
        real(wp), intent(in)    :: dt             ! [s]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: tcol(vg%nlev), qcol(vg%nlev), t0(vg%nlev), q0(vg%nlev)
        real(wp) :: pcol
        integer  :: i, j, k, nlev

        if (.not. cnv%enabled) return

        nlev = vg%nlev

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,phalf,pfull,dpc,tcol,qcol,t0,q0,pcol)
        do j = 1, cnv%nlat
            do i = 1, cnv%nlon
                call aeros_vgrid_pressure(vg, exp(lnps_g(i,j)), phalf, pfull, dpc)

                do k = 1, nlev
                    t0(k)   = t_g(i,j,k);  tcol(k) = t_g(i,j,k)
                    q0(k)   = qv_g(i,j,k); qcol(k) = qv_g(i,j,k)
                end do

                call adjust_column(cnv%scheme, tcol, qcol, pfull, dpc, nlev, &
                                        cnv%tau, cnv%rh_ref, dt)

                ! Temperature change as a forward increment [K]; humidity set to
                ! the adjusted profile; precip is the column water removed.
                pcol = 0.0_wp
                do k = 1, nlev
                    dt_phys(i,j,k) = dt_phys(i,j,k) + (tcol(k) - t0(k))
                    qv_g(i,j,k)    = qcol(k)
                    pcol           = pcol + (q0(k) - qcol(k))*dpc(k)
                end do
                cnv%precip(i,j) = pcol/(real(grav, wp)*dt)
            end do
        end do
        !$omp end parallel do

        return

    end subroutine aeros_convection_apply

    subroutine adjust_column(scheme, t, q, pfull, dp, nlev, tau, rh_ref, dt)
        ! Dispatch to the convective scheme. tau/rh_ref/dt are the Simplified
        ! Betts-Miller parameters and are ignored by Manabe.

        implicit none

        integer,  intent(in)    :: scheme
        real(wp), intent(inout) :: t(:), q(:)
        real(wp), intent(in)    :: pfull(:), dp(:)
        integer,  intent(in)    :: nlev
        real(wp), intent(in)    :: tau, rh_ref, dt

        select case (scheme)
        case (SCHEME_SBM)
            call sbm_adjust(t, q, pfull, dp, nlev, tau, rh_ref, dt)
        case (SCHEME_MANABE)
            call manabe_adjust(t, q, pfull, dp, nlev)
        case default
            call sbm_adjust(t, q, pfull, dp, nlev, tau, rh_ref, dt)
        end select

        return

    end subroutine adjust_column

    subroutine sbm_adjust(t, q, pfull, dp, nlev, tau, rh_ref, dt)
        ! Simplified Betts-Miller (Frierson 2007), one column.
        !
        ! Reference moist adiabat anchored on the ACTUAL near-surface moist static
        ! energy (the boundary-layer theta_e),
        !
        !   h_b = cp T_s + Phi_s + L q_s      (q_s the actual surface humidity),
        !
        ! not the saturated value: the saturated anchor gives a reference so warm
        ! that rh_ref q_sat(T_ref) exceeds the environmental humidity everywhere
        ! and no column ever rains. The reference temperature is the saturated
        ! parcel carrying that MSE at each level,
        !
        !   cp T_ref(k) + Phi(k) + L q_sat(T_ref(k), p_k) = h_b      (Newton),
        !
        ! and the convecting band is the layer where the parcel is buoyant, in the
        ! MSE sense h_b > h*_env(k) = cp T + Phi + L q_sat(T). Using MSE buoyancy
        ! rather than a raw T comparison folds in the LCL for free: an unsaturated
        ! boundary layer has h*_env > h_b at the surface, so the sub-cloud layer
        ! falls out of the band and no dry-adiabat / LCL bookkeeping is needed --
        ! the band is the cloud layer (level of free convection to level of
        ! neutral buoyancy). Reference humidity q_ref = rh_ref q_sat(T_ref).
        !
        ! Energy is closed as in Frierson. With column integrals over the band,
        !   Pq = int (q - q_ref) dp   (water available to remove),
        !   Pt = int cp (T_ref - T) dp   (candidate heating):
        !
        !   DEEP (Pq > 0): shift T_ref by the constant that makes
        !     int cp (T_ref - T) dp = L Pq,  so the heating balances L * precip,
        !     then relax T and q toward the reference. Precipitates.
        !
        !   SHALLOW (Pq <= 0, too dry to rain): a non-precipitating, humidity-
        !     conserving heat adjustment. Shift T_ref so int cp (T_ref - T) dp = 0
        !     and relax TEMPERATURE ONLY; q is left untouched. (This is a
        !     deliberate simplification of Betts-Miller's moisture-redistributing
        !     shallow branch -- a moisture relaxation toward a shifted reference
        !     can drive a very dry column negative, and there is nothing to tune
        !     its reference against until there is radiation and observations.)
        !
        ! Both branches relax IMPLICITLY over one physics step,
        !   a = (dt/tau)/(1 + dt/tau),   X <- X + a (X_ref - X),
        ! which is unconditionally stable in tau and conserves int (cp T + L q) dp
        ! to machine precision (heating = latent heat of the water removed), with
        ! precip equal to the column drying exactly and q >= 0 by construction.
        ! Only the band is touched; a column with no buoyant layer is returned
        ! unchanged. Dry (unsaturated) convective adjustment is out of scope --
        ! Betts-Miller is a moist scheme.

        implicit none

        real(wp), intent(inout) :: t(:), q(:)
        real(wp), intent(in)    :: pfull(:), dp(:)
        integer,  intent(in)    :: nlev
        real(wp), intent(in)    :: tau, rh_ref, dt

        real(wp) :: phi(nlev), tref(nlev), qref(nlev), hstar(nlev)
        real(wp) :: hb, qs, dqs, cp, lv, a
        real(wp) :: sumdp, pq, pt, shift, ct
        integer  :: k, ktop, kb

        cp = real(cp_d, wp); lv = real(L_v, wp)

        ! Geopotential relative to the surface, from the environment profile.
        phi(nlev) = 0.0_wp
        do k = nlev - 1, 1, -1
            phi(k) = phi(k+1) + R_d*0.5_wp*(t(k) + t(k+1))*log(pfull(k+1)/pfull(k))
        end do

        ! Boundary-layer MSE (actual humidity) anchors the reference adiabat, and
        ! the saturated environmental MSE sets the buoyancy of a lifted parcel.
        hb = cp*t(nlev) + phi(nlev) + lv*q(nlev)
        do k = 1, nlev
            call aeros_qsat(t(k), pfull(k), qs, dqs)
            hstar(k) = cp*t(k) + phi(k) + lv*qs
        end do

        ! Convecting band: the first contiguous layer, scanning up from the
        ! surface, where the parcel is buoyant (h_b > h*_env). kb is its base
        ! (level of free convection), ktop its top (level of neutral buoyancy).
        kb = 0; ktop = 0
        do k = nlev, 1, -1
            if (hb > hstar(k)) then
                if (kb == 0) kb = k
                ktop = k
            else if (kb /= 0) then
                exit
            end if
        end do
        if (kb == 0) return         ! no buoyant layer -> no convection

        ! A single-layer band (ktop == kb) is not convection: there is no layer
        ! above to overturn into, and the deep closure then anchors the reference
        ! adiabat on the one layer itself and dumps L*precip back into it as local
        ! heating -- a spurious positive feedback that fires when surface fluxes
        ! push the lowest layer to saturation (h_b > h*_env there). That is
        ! large-scale condensation's job, which runs at the same seam right after;
        ! convection needs a genuine cloud depth. Return the column unchanged.
        if (ktop == kb) return

        ! Reference profiles over the band: the moist adiabat carrying h_b, dried
        ! to rh_ref of its own saturation.
        do k = ktop, kb
            tref(k) = moist_adiabat_temp(hb - phi(k), pfull(k), t(k))
            call aeros_qsat(tref(k), pfull(k), qs, dqs)
            qref(k) = rh_ref*qs
        end do

        ! Column integrals over the band.
        sumdp = 0.0_wp; pq = 0.0_wp; pt = 0.0_wp
        do k = ktop, kb
            sumdp = sumdp + dp(k)
            pq    = pq + (q(k)    - qref(k))*dp(k)
            pt    = pt + cp*(tref(k) - t(k))*dp(k)
        end do

        a = (dt/tau)/(1.0_wp + dt/tau)

        if (pq > 0.0_wp) then
            ! Deep: shift T_ref so the heating balances L * precip exactly, then
            ! relax both temperature and humidity toward the reference.
            shift = (lv*pq - pt)/(cp*sumdp)
            do k = ktop, kb
                t(k) = t(k) + a*(tref(k) + shift - t(k))
                q(k) = q(k) + a*(qref(k) - q(k))
            end do
        else
            ! Shallow: shift T_ref to zero net heating and relax TEMPERATURE ONLY.
            ! No precipitation, column moisture untouched, q >= 0 by construction.
            ct = pt/(cp*sumdp)
            do k = ktop, kb
                t(k) = t(k) + a*(tref(k) - ct - t(k))
            end do
        end if

        return

    end subroutine sbm_adjust

    function moist_adiabat_temp(rhs, p, tguess) result(tsol)
        ! Solve cp T + L q_sat(T, p) = rhs for T -- the temperature of a
        ! saturated parcel with the given cp T + L q_sat -- by Newton from
        ! tguess against dq_sat/dT.

        implicit none

        real(wp), intent(in) :: rhs, p, tguess
        real(wp) :: tsol

        real(wp) :: qs, dqs, f, fp, cp, lv
        integer  :: it

        cp = real(cp_d, wp); lv = real(L_v, wp)

        tsol = tguess
        do it = 1, 8
            call aeros_qsat(tsol, p, qs, dqs)
            f    = cp*tsol + lv*qs - rhs
            fp   = cp + lv*dqs
            tsol = tsol - f/fp
        end do

        return

    end function moist_adiabat_temp

    subroutine manabe_adjust(t, q, pfull, dp, nlev)
        ! One column, adjusted to stability by repeated pairwise relaxation.

        implicit none

        real(wp), intent(inout) :: t(:), q(:)
        real(wp), intent(in)    :: pfull(:), dp(:)
        integer,  intent(in)    :: nlev

        real(wp) :: phi(nlev), theta(nlev)
        real(wp) :: qsk, qsk1, dqs, hstar_k, hstar_k1
        logical  :: changed, sat_k, sat_k1
        integer  :: k, sweep

        do sweep = 1, MAXSWEEP
            changed = .FALSE.

            ! Geopotential relative to the surface, from the current profile.
            ! Recomputed each sweep because the adjustment changes T. Only the
            ! layer-to-layer differences enter the moist-neutral test, so the
            ! surface reference is arbitrary.
            phi(nlev) = 0.0_wp
            do k = nlev - 1, 1, -1
                phi(k) = phi(k+1) + R_d*0.5_wp*(t(k) + t(k+1))*log(pfull(k+1)/pfull(k))
            end do

            do k = 1, nlev
                theta(k) = t(k)*(real(p0,wp)/pfull(k))**real(kappa,wp)
            end do

            do k = 1, nlev - 1
                ! k is the upper layer (smaller pressure), k+1 the lower.
                call aeros_qsat(t(k),   pfull(k),   qsk,  dqs)
                call aeros_qsat(t(k+1), pfull(k+1), qsk1, dqs)

                ! Saturated to within a hair counts as saturated: convection
                ! acts on cloud, and the transport/condensation leave layers
                ! sitting essentially at q_sat.
                sat_k  = q(k)   >= 0.999_wp*qsk
                sat_k1 = q(k+1) >= 0.999_wp*qsk1

                if (sat_k .and. sat_k1) then
                    ! Moist instability: saturated MSE decreasing upward.
                    hstar_k  = cp_d*t(k)   + phi(k)   + L_v*qsk
                    hstar_k1 = cp_d*t(k+1) + phi(k+1) + L_v*qsk1
                    if (hstar_k < hstar_k1 - 1.0e-6_wp) then
                        call moist_pair(t, q, pfull, dp, phi(k) - phi(k+1), k)
                        changed = .TRUE.
                    end if
                else
                    ! Dry instability: potential temperature decreasing upward.
                    if (theta(k) < theta(k+1) - 1.0e-6_wp) then
                        call dry_pair(t, pfull, dp, k)
                        changed = .TRUE.
                    end if
                end if
            end do

            if (.not. changed) exit
        end do

        return

    end subroutine manabe_adjust

    subroutine dry_pair(t, pfull, dp, k)
        ! Mix a dry-unstable pair to a common potential temperature, conserving
        ! enthalpy int c_p T dp. c_p cancels, so it drops out of the algebra.

        implicit none

        real(wp), intent(inout) :: t(:)
        real(wp), intent(in)    :: pfull(:), dp(:)
        integer,  intent(in)    :: k

        real(wp) :: ex1, ex2, theta_new

        ex1 = (pfull(k)  /real(p0,wp))**real(kappa,wp)
        ex2 = (pfull(k+1)/real(p0,wp))**real(kappa,wp)

        ! int c_p T dp conserved => sum T dp = theta_new sum (p/p0)^kappa dp.
        theta_new = (t(k)*dp(k) + t(k+1)*dp(k+1)) / (ex1*dp(k) + ex2*dp(k+1))

        t(k)   = theta_new*ex1
        t(k+1) = theta_new*ex2

        return

    end subroutine dry_pair

    subroutine moist_pair(t, q, pfull, dp, dphi, k)
        ! Relax a moist-unstable saturated pair to a moist adiabat, conserving
        ! the moist static energy int (c_p T + L q) dp and holding both layers
        ! at saturation. Two unknowns (the new temperatures), two constraints,
        ! Newton against dq_sat/dT.
        !
        !   MSE:   c_p(T1 dp1 + T2 dp2) + L(q_s(T1) dp1 + q_s(T2) dp2) = H0
        !   neutral: (c_p T1 + L q_s(T1)) - (c_p T2 + L q_s(T2)) + dphi = 0
        !
        ! dphi = Phi(k) - Phi(k+1) > 0 is held fixed across the Newton solve and
        ! refreshed by the outer sweep. On convergence q_i = q_s(T_i), and the
        ! water that left is the precipitation the caller accounts.

        implicit none

        real(wp), intent(inout) :: t(:), q(:)
        real(wp), intent(in)    :: pfull(:), dp(:), dphi
        integer,  intent(in)    :: k

        real(wp) :: t1, t2, qs1, qs2, d1, d2
        real(wp) :: h0, rh, rn, j11, j12, j21, j22, det, det_t1, det_t2
        real(wp) :: cp, lv, dp1, dp2
        integer  :: it

        cp = real(cp_d, wp); lv = real(L_v, wp)
        dp1 = dp(k); dp2 = dp(k+1)

        t1 = t(k); t2 = t(k+1)
        call aeros_qsat(t1, pfull(k),   qs1, d1)
        call aeros_qsat(t2, pfull(k+1), qs2, d2)

        ! Conserved moist static energy, from the pre-adjustment saturated state.
        h0 = cp*(t1*dp1 + t2*dp2) + lv*(qs1*dp1 + qs2*dp2)

        do it = 1, 6
            call aeros_qsat(t1, pfull(k),   qs1, d1)
            call aeros_qsat(t2, pfull(k+1), qs2, d2)

            rh = cp*(t1*dp1 + t2*dp2) + lv*(qs1*dp1 + qs2*dp2) - h0
            rn = (cp*t1 + lv*qs1) - (cp*t2 + lv*qs2) + dphi

            j11 = (cp + lv*d1)*dp1;   j12 = (cp + lv*d2)*dp2
            j21 =  cp + lv*d1;        j22 = -(cp + lv*d2)

            det = j11*j22 - j12*j21
            if (abs(det) < 1.0e-30_wp) exit

            det_t1 = rh*j22 - j12*rn
            det_t2 = j11*rn - rh*j21

            t1 = t1 - det_t1/det
            t2 = t2 - det_t2/det
        end do

        call aeros_qsat(t1, pfull(k),   qs1, d1)
        call aeros_qsat(t2, pfull(k+1), qs2, d2)

        t(k)   = t1;  q(k)   = qs1
        t(k+1) = t2;  q(k+1) = qs2

        return

    end subroutine moist_pair

    subroutine aeros_convection_report(cnv, io_unit)

        implicit none

        type(aeros_conv_class), intent(in) :: cnv
        integer, intent(in), optional :: io_unit

        integer :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        write(iou,*) ""
        write(iou,"(a)") " == convection =="
        if (cnv%enabled) then
            select case (cnv%scheme)
            case (SCHEME_SBM)
                write(iou,"(a)")           "   moist convective adjustment ON  (Simplified Betts-Miller)"
                write(iou,"(a,f9.1,a)")    "   relaxation timescale       ", cnv%tau, " s"
                write(iou,"(a,f9.3)")      "   reference relative humidity", cnv%rh_ref
            case (SCHEME_MANABE)
                write(iou,"(a)")           "   moist convective adjustment ON  (Manabe)"
            end select
        else
            write(iou,"(a)") "   moist convective adjustment off"
        end if

        return

    end subroutine aeros_convection_report

end module aeros_convection
