module aeros_budget
    ! Global integrals of the model state: mass, energy, angular momentum.
    !
    ! These are the quantities docs/design.md requires to close -- section 3.7
    ! risk 2 and section 7 (M4) both say "to machine precision", and section
    ! 1's whole premise is integrations of 10^4-10^6 yr, over which a drift of
    ! one part in 10^8 per step is not a rounding error but a climate.
    !
    ! Until M1.4 only mass was checkable: a single tendency evaluation says
    ! nothing about what an integration conserves. All three are checkable now,
    ! and they behave differently, which is the point of measuring rather than
    ! asserting:
    !
    !   MASS is conserved to machine precision and structurally so. The
    !   surface-pressure tendency is minus the divergence of a vertically
    !   integrated flux, the Gaussian quadrature integrates that exactly, and
    !   the vertical mass flux vanishes at both boundaries by construction
    !   (aeros_tendency). Nothing in the timestep can break it except the time
    !   filter, which is why the filter is applied to ln(p_s) too and the
    !   diffusion is not.
    !
    !   ENERGY and ANGULAR MOMENTUM are NOT conserved to machine precision, and
    !   no honest implementation of this scheme would claim otherwise. The
    !   Simmons-Burridge vertical discretization conserves energy in the
    !   continuous-time limit, but leapfrog, the time filter and the
    !   hyperdiffusion each remove some; diffusion is a deliberate sink. What
    !   matters is the RATE, and that it is dominated by the terms one meant to
    !   include.
    !
    ! Everything accumulates in dp regardless of wp, for the same reason
    ! aeros_sht_surface_integral does: a sum over nlon*nlat*nlev terms in
    ! single precision would cap the check far above the precision it is
    ! supposed to measure.
    !
    ! The state's GRID fields are the input, so aeros_timestep_diagnose must
    ! have been called for the time level being measured. Deliberately not the
    ! spectral state: a budget computed from spectral coefficients would be a
    ! statement about the transform, and what has to close is the model's
    ! physical state on the grid the physics will run on.

    use aeros_defs,     only : dp, wp, io_unit_err, grav, cp_d, omega, r_earth, &
                                year_seconds, aeros_grid_class, aeros_state_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure

    implicit none

    private

    type aeros_budget_class
        real(dp) :: mass   = 0.0_dp    ! total dry air mass [kg]
        real(dp) :: energy = 0.0_dp    ! total energy [J]
        real(dp) :: angmom = 0.0_dp    ! axial angular momentum [kg m2 s-1]
    end type aeros_budget_class

    public :: aeros_budget_class
    public :: aeros_budget_calc
    public :: aeros_budget_report

contains

    subroutine aeros_budget_calc(bud, vg, grd, now, phis)
        ! Global integrals from the grid-space state.
        !
        !   mass   = int  p_s / g  dA
        !   energy = int  [ sum_k (cp T_k + K_k) dp_k / g  +  Phi_s p_s / g ] dA
        !   angmom = int  sum_k (u_k + Omega a cos(lat)) a cos(lat) dp_k / g  dA
        !
        ! The energy is the hydrostatic total energy -- enthalpy rather than
        ! internal energy, which is what a pressure-coordinate model conserves;
        ! the surface term is what carries the work done against orography and
        ! is zero only while the surface is flat.
        !
        ! The angular momentum's second term is the OMEGA contribution, the
        ! angular momentum of the atmosphere co-rotating with the Earth. It
        ! dwarfs the wind term by ~2 orders of magnitude and is not a constant
        ! of the discrete motion either, since it moves with the mass -- which
        ! is exactly why it must be included rather than diagnosing the winds
        ! alone.

        implicit none

        type(aeros_budget_class), intent(out) :: bud
        type(aeros_vgrid_class),  intent(in)  :: vg
        type(aeros_grid_class),   intent(in)  :: grd
        type(aeros_state_class),  intent(in)  :: now
        real(wp), intent(in), optional :: phis(:,:)   ! surface geopotential [m2 s-2]

        real(wp) :: p_half(0:vg%nlev), p_full(vg%nlev), dp_lev(vg%nlev)
        real(dp) :: colmass, colener, colamom, coslat, area, ps, ginv
        real(dp) :: m_tot, e_tot, a_tot
        integer  :: i, j, k

        ginv  = 1.0_dp/grav
        m_tot = 0.0_dp; e_tot = 0.0_dp; a_tot = 0.0_dp

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,p_half,p_full,dp_lev,colmass,colener,colamom,coslat,area,ps) &
        !$omp   reduction(+:m_tot,e_tot,a_tot)
        do j = 1, grd%nlat
            do i = 1, grd%nlon

                ps   = real(now%ps(i,j), dp)
                area = real(grd%area(i,j), dp)

                ! cos(latitude) = sin(colatitude); grd%colat comes from SHTns'
                ! own Gauss nodes, so this is the same geometry the transforms
                ! and the quadrature weights use.
                coslat = sin(real(grd%colat(j), dp))

                call aeros_vgrid_pressure(vg, real(ps, wp), p_half, p_full, dp_lev)

                colener = 0.0_dp
                colamom = 0.0_dp
                do k = 1, vg%nlev
                    colener = colener + real(dp_lev(k), dp) &
                                *( cp_d*real(now%temp_g(i,j,k), dp) &
                                    + 0.5_dp*( real(now%u(i,j,k), dp)**2 &
                                                + real(now%v(i,j,k), dp)**2 ) )
                    colamom = colamom + real(dp_lev(k), dp) &
                                *( real(now%u(i,j,k), dp) + omega*r_earth*coslat ) &
                                *r_earth*coslat
                end do

                colmass = ps

                if (present(phis)) colener = colener + real(phis(i,j), dp)*ps

                m_tot = m_tot + area*colmass*ginv
                e_tot = e_tot + area*colener*ginv
                a_tot = a_tot + area*colamom*ginv

            end do
        end do
        !$omp end parallel do

        bud%mass   = m_tot
        bud%energy = e_tot
        bud%angmom = a_tot

        return

    end subroutine aeros_budget_calc

    subroutine aeros_budget_report(bud, ref, elapsed, label, io_unit)
        ! Print a budget against a reference one, as a drift RATE.
        !
        ! Per unit time rather than per step, and relative rather than
        ! absolute, because the question a paleo integration asks is "how much
        ! of this is left after 10^5 yr", and that is only answerable from a
        ! rate. `elapsed` is the model time between the two budgets [s].

        implicit none

        type(aeros_budget_class), intent(in) :: bud
        type(aeros_budget_class), intent(in) :: ref
        real(dp), intent(in) :: elapsed
        character(len=*), intent(in) :: label
        integer, intent(in), optional :: io_unit

        real(dp) :: yr
        integer  :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        if (elapsed <= 0.0_dp) then
            write(io_unit_err,*) "aeros_budget_report:: error: elapsed must be > 0, got ", elapsed
            error stop 1
        end if

        ! Drift per YEAR, on the same calendar the facade converts with.
        yr = year_seconds

        write(iou,*) ""
        write(iou,"(a,a)") " == conservation: ", trim(label)
        write(iou,"(a,es14.6,a)")  "   mass            ", bud%mass,   " kg"
        write(iou,"(a,es14.6,a)")  "   energy          ", bud%energy, " J"
        write(iou,"(a,es14.6,a)")  "   angular momentum", bud%angmom, " kg m2 s-1"
        write(iou,"(a,es12.3)")    "   mass   drift [/yr] ", reldrift(bud%mass,   ref%mass,   elapsed, yr)
        write(iou,"(a,es12.3)")    "   energy drift [/yr] ", reldrift(bud%energy, ref%energy, elapsed, yr)
        write(iou,"(a,es12.3)")    "   angmom drift [/yr] ", reldrift(bud%angmom, ref%angmom, elapsed, yr)
        write(iou,*) ""

        return

    end subroutine aeros_budget_report

    real(dp) function reldrift(x, x0, elapsed, yr) result(r)

        implicit none

        real(dp), intent(in) :: x, x0, elapsed, yr

        if (x0 == 0.0_dp) then
            r = 0.0_dp
        else
            r = (x - x0)/x0*(yr/elapsed)
        end if

        return

    end function reldrift

end module aeros_budget
