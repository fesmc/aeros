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
    ! === Why this shape, with Betts-Miller in mind ==========================
    !
    ! Manabe was chosen over Betts-Miller for the first cut because it conserves
    ! by construction -- it redistributes rather than relaxing toward an external
    ! reference, so there is no timescale and no reference humidity to tune and
    ! no energy-closure correction to get right. Betts-Miller is the better
    ! *tunable* scheme once there is radiation and observations to tune against.
    !
    ! So the seam is built for both: aeros_conv_class carries a `scheme`
    ! selector, aeros_convection_apply dispatches per column to adjust_column,
    ! and adjust_column branches on the scheme. Adding Betts-Miller is a second
    ! branch and a reference-profile routine -- the apply loop, the tendency
    ! plumbing, the precipitation accounting and the tests do not change.
    !
    ! === Coupling: same two seams as condensation ===========================
    !
    ! Convection changes temperature (up AND down -- it is a redistribution) and
    ! humidity. Like condensation it acts at the grid seam: the temperature
    ! change enters wrk%dtdt as a rate so it rides the transform and the
    ! leapfrog, the humidity change is applied to the gridpoint qv_g, and both
    ! come from the same adjustment so the column budget closes. It runs BEFORE
    ! large-scale condensation, which then mops up whatever supersaturation the
    ! convection did not remove.

    use aeros_defs,       only : dp, wp, io_unit_err, R_d, cp_d, grav, L_v, &
                                    kappa, p0, aeros_grid_class
    use aeros_vertical,   only : aeros_vgrid_class, aeros_vgrid_pressure
    use aeros_condensation, only : aeros_qsat

    use nml,              only : nml_read

    implicit none

    private

    integer, parameter :: SCHEME_MANABE = 1

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
        integer :: scheme  = SCHEME_MANABE

        integer :: nlon = 0, nlat = 0

        ! Precipitation rate from the last apply, [kg m-2 s-1].
        real(wp), allocatable :: precip(:,:)
    end type aeros_conv_class

    public :: aeros_conv_class
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
        cnv%scheme  = SCHEME_MANABE
        cnv%nlon    = grd%nlon
        cnv%nlat    = grd%nlat

        allocate(cnv%precip(grd%nlon, grd%nlat))
        cnv%precip = 0.0_wp

        return

    end subroutine aeros_convection_init

    subroutine aeros_convection_load(cnv, filename, grd)
        ! Configure from the `aeros_moisture` namelist group. `convect` turns it
        ! on; `conv_scheme` selects, with only "manabe" implemented -- an
        ! unknown name is an error, not a silent fallback, so a typo does not
        ! quietly run a scheme the user did not ask for.

        implicit none

        type(aeros_conv_class), intent(inout) :: cnv
        character(len=*),       intent(in)    :: filename
        type(aeros_grid_class), intent(in)    :: grd

        logical            :: convect
        character(len=32)  :: conv_scheme

        convect     = .FALSE.
        conv_scheme = "manabe"
        call nml_read(filename, "aeros_moisture", "convect",     convect)
        call nml_read(filename, "aeros_moisture", "conv_scheme", conv_scheme)

        call aeros_convection_init(cnv, grd, convect)

        select case (trim(conv_scheme))
        case ("manabe")
            cnv%scheme = SCHEME_MANABE
        case default
            write(io_unit_err,*) "aeros_convection_load:: error: unknown conv_scheme '"// &
                                    trim(conv_scheme)//"' (only 'manabe' is implemented)"
            error stop 1
        end select

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

    subroutine aeros_convection_apply(cnv, vg, t_g, qv_g, lnps_g, dtdt, dt)
        ! Adjust each column for moist (and dry) convective instability: dry
        ! qv_g, add the temperature change to dtdt as a rate, record precip.

        implicit none

        type(aeros_conv_class),  intent(inout) :: cnv
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(inout) :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat)      ln[Pa]
        real(wp), intent(inout) :: dtdt(:,:,:)    ! (nlon,nlat,nlev) [K s-1]
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

                call adjust_column(cnv%scheme, tcol, qcol, pfull, dpc, nlev)

                ! Heating as a rate; humidity set to the adjusted profile;
                ! precip is the column water removed.
                pcol = 0.0_wp
                do k = 1, nlev
                    dtdt(i,j,k) = dtdt(i,j,k) + (tcol(k) - t0(k))/dt
                    qv_g(i,j,k) = qcol(k)
                    pcol        = pcol + (q0(k) - qcol(k))*dpc(k)
                end do
                cnv%precip(i,j) = pcol/(real(grav, wp)*dt)
            end do
        end do
        !$omp end parallel do

        return

    end subroutine aeros_convection_apply

    subroutine adjust_column(scheme, t, q, pfull, dp, nlev)
        ! Dispatch to the convective scheme. Betts-Miller becomes a second case.

        implicit none

        integer,  intent(in)    :: scheme
        real(wp), intent(inout) :: t(:), q(:)
        real(wp), intent(in)    :: pfull(:), dp(:)
        integer,  intent(in)    :: nlev

        select case (scheme)
        case (SCHEME_MANABE)
            call manabe_adjust(t, q, pfull, dp, nlev)
        case default
            call manabe_adjust(t, q, pfull, dp, nlev)
        end select

        return

    end subroutine adjust_column

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
            write(iou,"(a)") "   moist convective adjustment ON  (Manabe)"
        else
            write(iou,"(a)") "   moist convective adjustment off"
        end if

        return

    end subroutine aeros_convection_report

end module aeros_convection
