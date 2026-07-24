module aeros_held_suarez
    ! The Held & Suarez (1994) idealized forcing.
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the intercomparison
    !   of the dynamical cores of atmospheric general circulation models.
    !   Bull. Amer. Meteor. Soc., 75, 1825-1830.
    !
    ! Newtonian relaxation of temperature toward a zonally symmetric
    ! radiative-equilibrium profile, plus Rayleigh drag on the low-level winds.
    ! No moisture, no radiation, no topography, no seasonal cycle. It exists to
    ! test a dynamical core against other dynamical cores, and it is
    ! docs/design.md's stated validation for M1.
    !
    ! === The specification ===================================================
    !
    !   T_eq(phi,p) = max{ 200 K,
    !                      [ 315 - dT_y sin^2(phi) - dtheta_z log(p/p0) cos^2(phi) ]
    !                      (p/p0)^kappa }
    !
    !   k_T(phi,sigma) = k_a + (k_s - k_a) max(0, (sigma - sigma_b)/(1 - sigma_b))
    !                                          cos^4(phi)
    !
    !   k_v(sigma)     = k_f max(0, (sigma - sigma_b)/(1 - sigma_b))
    !
    !   dT/dt |_forcing = -k_T (T - T_eq)
    !   dv/dt |_forcing = -k_v v
    !
    !   dT_y = 60 K, dtheta_z = 10 K, sigma_b = 0.7, p0 = 10^5 Pa,
    !   k_a = 1/40 day^-1, k_s = 1/4 day^-1, k_f = 1 day^-1.
    !
    ! Two details the restatements differ on, both settled here deliberately:
    !
    !   SIGMA is p/p_s, the ORIGINAL definition, not p/p0. Some restatements
    !   (MITgcm's, for one) write p/p0, which is the same thing over a flat
    !   surface at exactly 1000 hPa and not otherwise. p/p_s is what the paper
    !   says and what a sigma-coordinate model means by sigma, and it is the
    !   one that keeps the boundary layer at a fixed fraction of the column
    !   when the surface pressure moves.
    !
    !   KAPPA is the model's own R_d/c_p (aeros_defs, 0.285716), not the 2/7
    !   the paper quotes as a round number. They differ by 6 parts in a
    !   million. Using the model's value keeps T_eq an exact potential
    !   temperature in the model's own thermodynamics, which is what the
    !   specification means.
    !
    ! === The constants are not parameters ====================================
    !
    ! dT_y, dtheta_z, sigma_b and the three rates are `parameter`s here rather
    ! than namelist entries. That is the point of a benchmark: a run with
    ! dT_y = 50 K is not Held-Suarez and cannot be compared with anyone else's
    ! Held-Suarez. Only whether the forcing is ON is configurable
    ! (`held_suarez` in the `aeros` namelist group). If M-something later wants
    ! a modified-HS sensitivity study, that is a new module with its own name,
    ! not a loosened version of this one.
    !
    ! === Where this is applied, and why it matters ===========================
    !
    ! Between aeros_tendency_grid and aeros_tendency_spectral, on the GRID.
    !
    ! For the thermal relaxation that is merely convenient. For the drag it is
    ! required: k_v depends on sigma and therefore on the local surface
    ! pressure, so -k_v*v is not -k_v*zeta and -k_v*D in spectral space, k_v
    ! not commuting with the curl. Applied to the grid-space momentum RHS
    ! before the vector analysis, it is exact.
    !
    ! === Time discretization =================================================
    !
    ! The forcing is EXPLICIT, evaluated at the current time level along with
    ! the rest of the right-hand side.
    !
    ! A relaxation term in a leapfrog destabilizes the computational mode --
    ! -2 dt k X^n feeds the alternating solution rather than damping it -- and
    ! the textbook remedy is to lag the damping to X^(n-1). That is not done
    ! here: the fastest rate is k_f = 1 day^-1, so k_f * 2 dt = 0.042 at
    ! dt = 1800 s, an order of magnitude below the RAW filter's own damping of
    ! that mode, and lagging would cost a full extra synthesis of the previous
    ! state every step. If dt ever rises far enough that k_f*2dt approaches
    ! the filter strength, revisit this before blaming the dynamics.

    use aeros_defs,     only : dp, wp, io_unit_err, kappa, &
                                aeros_grid_class, aeros_spec_class
    use aeros_spectral, only : aeros_sht_class, aeros_sht_lm
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure
    use aeros_tendency, only : aeros_work_class

    implicit none

    private

    ! === Held & Suarez (1994), section 2 =====================================

    real(wp), parameter, public :: hs_t_equator = 315.0_wp   ! [K]
    real(wp), parameter, public :: hs_t_min     = 200.0_wp   ! [K], stratospheric floor
    real(wp), parameter, public :: hs_dt_y      =  60.0_wp   ! equator-pole contrast [K]
    real(wp), parameter, public :: hs_dtheta_z  =  10.0_wp   ! vertical stability [K]
    real(wp), parameter, public :: hs_sigma_b   =   0.7_wp   ! top of the damped layer
    real(wp), parameter, public :: hs_p0        = 1.0e5_wp   ! reference pressure [Pa]

    real(wp), parameter, public :: hs_k_a = 1.0_wp/(40.0_wp*86400.0_wp)  ! [s-1]
    real(wp), parameter, public :: hs_k_s = 1.0_wp/( 4.0_wp*86400.0_wp)  ! [s-1]
    real(wp), parameter, public :: hs_k_f = 1.0_wp/(        86400.0_wp)  ! [s-1]

    type aeros_hs_class
        ! The latitude factors, precomputed once.
        !
        ! From grd%sinlat -- SHTns' Gauss nodes verbatim -- rather than from
        ! grd%lat, so sin^2 and cos^4 are EXACTLY equal in the two hemispheres.
        ! The forcing is latitudinally symmetric by construction, and the
        ! headline diagnostic of a Held-Suarez run is how symmetric the
        ! resulting circulation is; a forcing asymmetric at 1e-16 is harmless
        ! physically but there is no reason to introduce it.

        logical :: enabled = .FALSE.
        integer :: nlat = 0

        real(wp), allocatable :: sinlat2(:)   ! sin^2(lat), (nlat)
        real(wp), allocatable :: coslat2(:)   ! cos^2(lat)
        real(wp), allocatable :: coslat4(:)   ! cos^4(lat)
    end type aeros_hs_class

    public :: aeros_hs_class
    public :: aeros_hs_init
    public :: aeros_hs_end
    public :: aeros_hs_apply
    public :: aeros_hs_perturb
    public :: aeros_hs_print
    public :: aeros_hs_teq
    public :: aeros_hs_kt
    public :: aeros_hs_kv

contains

    subroutine aeros_hs_init(hs, grd, enabled)

        implicit none

        type(aeros_hs_class),   intent(inout) :: hs
        type(aeros_grid_class), intent(in)    :: grd
        logical, intent(in) :: enabled

        integer :: j

        call aeros_hs_end(hs)

        hs%enabled = enabled
        hs%nlat    = grd%nlat

        allocate(hs%sinlat2(grd%nlat), hs%coslat2(grd%nlat), hs%coslat4(grd%nlat))

        do j = 1, grd%nlat
            hs%sinlat2(j) = real(grd%sinlat(j)*grd%sinlat(j), wp)
            hs%coslat2(j) = 1.0_wp - hs%sinlat2(j)
            hs%coslat4(j) = hs%coslat2(j)*hs%coslat2(j)
        end do

        return

    end subroutine aeros_hs_init

    subroutine aeros_hs_end(hs)

        implicit none

        type(aeros_hs_class), intent(inout) :: hs

        hs%enabled = .FALSE.
        hs%nlat    = 0

        if (allocated(hs%sinlat2)) deallocate(hs%sinlat2)
        if (allocated(hs%coslat2)) deallocate(hs%coslat2)
        if (allocated(hs%coslat4)) deallocate(hs%coslat4)

        return

    end subroutine aeros_hs_end

    subroutine aeros_hs_apply(hs, vg, wrk)
        ! Add the forcing to the grid-space right-hand sides.
        !
        ! Threaded over columns, the same loop shape and for the same reason as
        ! aeros_tendency's column_terms: this is the embarrassingly parallel
        ! part, and it is the part that grows when real physics arrives at M2.

        implicit none

        type(aeros_hs_class),    intent(in)    :: hs
        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_work_class),  intent(inout) :: wrk

        real(wp) :: p_half(0:vg%nlev), p_full(vg%nlev)
        real(wp) :: ps, sig, teq, kt, kv
        integer  :: i, j, k

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,p_half,p_full,ps,sig,teq,kt,kv)
        do j = 1, wrk%nlat
            do i = 1, wrk%nlon

                ps = exp(wrk%lnps_g(i,j))
                call aeros_vgrid_pressure(vg, ps, p_half, p_full)

                do k = 1, vg%nlev

                    sig = p_full(k)/ps

                    teq = aeros_hs_teq(p_full(k), hs%sinlat2(j), hs%coslat2(j))
                    kt  = aeros_hs_kt(sig, hs%coslat4(j))
                    kv  = aeros_hs_kv(sig)

                    wrk%dtdt(i,j,k) = wrk%dtdt(i,j,k) - kt*(wrk%t_g(i,j,k) - teq)

                    wrk%ae(i,j,k) = wrk%ae(i,j,k) - kv*wrk%u(i,j,k)
                    wrk%an(i,j,k) = wrk%an(i,j,k) - kv*wrk%v(i,j,k)

                end do

            end do
        end do
        !$omp end parallel do

        return

    end subroutine aeros_hs_apply

    real(wp) function aeros_hs_teq(p, sinlat2, coslat2) result(teq)
        ! Radiative-equilibrium temperature [K], Held & Suarez eq. (1).

        implicit none

        real(wp), intent(in) :: p         ! pressure [Pa]
        real(wp), intent(in) :: sinlat2   ! sin^2(latitude)
        real(wp), intent(in) :: coslat2   ! cos^2(latitude)

        real(wp) :: pfac

        pfac = p/hs_p0

        teq = (hs_t_equator - hs_dt_y*sinlat2 - hs_dtheta_z*log(pfac)*coslat2) &
                *pfac**kappa

        ! The floor is an isothermal stratosphere, and it is what stops T_eq
        ! from falling without bound as p -> 0.
        teq = max(hs_t_min, teq)

        return

    end function aeros_hs_teq

    real(wp) function aeros_hs_kt(sigma, coslat4) result(kt)
        ! Thermal relaxation rate [s-1], Held & Suarez eq. (2).
        !
        ! k_a everywhere, increasing to k_s (ten times faster) at the surface
        ! in the tropics -- a crude stand-in for the fact that the tropical
        ! boundary layer is convectively coupled to the surface.

        implicit none

        real(wp), intent(in) :: sigma     ! p/p_s
        real(wp), intent(in) :: coslat4   ! cos^4(latitude)

        kt = hs_k_a + (hs_k_s - hs_k_a)*hfac(sigma)*coslat4

        return

    end function aeros_hs_kt

    real(wp) function aeros_hs_kv(sigma) result(kv)
        ! Rayleigh drag rate [s-1], Held & Suarez eq. (3).

        implicit none

        real(wp), intent(in) :: sigma     ! p/p_s

        kv = hs_k_f*hfac(sigma)

        return

    end function aeros_hs_kv

    real(wp) function hfac(sigma) result(f)
        ! max(0, (sigma - sigma_b)/(1 - sigma_b)): zero above sigma_b, rising
        ! linearly to 1 at the surface. The shared boundary-layer profile of
        ! both damping rates.

        implicit none

        real(wp), intent(in) :: sigma

        f = max(0.0_wp, (sigma - hs_sigma_b)/(1.0_wp - hs_sigma_b))

        return

    end function hfac

    subroutine aeros_hs_perturb(spec, sht, amp)
        ! Seed baroclinic instability with a small temperature perturbation.
        !
        ! The forcing and the initial state are both zonally symmetric, so
        ! nothing breaks the symmetry except round-off, and a run started
        ! without a seed spends an unpredictable stretch of time waiting for
        ! round-off to grow. Since the benchmark is about the STATISTICS of the
        ! equilibrated flow, the seed only has to be small and to project onto
        ! the unstable wavenumbers.
        !
        ! DETERMINISTIC rather than random: two runs of the same configuration
        ! should produce the same numbers, or a regression check has nothing to
        ! compare against. Held & Suarez do not specify a seed -- the
        ! equilibrium statistics do not depend on it -- so the choice is ours
        ! and being reproducible is the only thing that matters.
        !
        ! Zonal wavenumbers 1-8 at mid-troposphere, which is where the fastest
        ! growing baroclinic modes live.

        implicit none

        type(aeros_spec_class), intent(inout) :: spec
        type(aeros_sht_class),  intent(in)    :: sht
        real(wp), intent(in) :: amp           ! perturbation amplitude [K]

        integer :: k, l, m, lm

        do k = 1, spec%nlev
            do m = 1, min(8, sht%mmax)
                do l = m, min(m+4, sht%lmax)
                    lm = aeros_sht_lm(sht, l, m)
                    if (lm < 1) cycle
                    spec%temp(lm,k) = spec%temp(lm,k) &
                            + cmplx(real(amp,dp)*cos(real(3*l + 5*m + k, dp)), &
                                    real(amp,dp)*sin(real(7*l + 2*m + 3*k, dp)), dp)
                end do
            end do
        end do

        return

    end subroutine aeros_hs_perturb

    subroutine aeros_hs_print(hs, io_unit)

        implicit none

        type(aeros_hs_class), intent(in) :: hs
        integer, intent(in), optional :: io_unit

        integer :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        if (.not. hs%enabled) return

        write(iou,*) ""
        write(iou,"(a)")        " == Held-Suarez (1994) forcing =="
        write(iou,"(a,f9.2,a)") "   T at the equatorial surface", hs_t_equator, " K"
        write(iou,"(a,f9.2,a)") "   stratospheric floor        ", hs_t_min, " K"
        write(iou,"(a,f9.2,a)") "   equator-pole contrast      ", hs_dt_y, " K"
        write(iou,"(a,f9.2,a)") "   vertical stability         ", hs_dtheta_z, " K"
        write(iou,"(a,f9.2)")   "   sigma_b                    ", hs_sigma_b
        write(iou,"(a,f9.2,a)") "   1/k_a                      ", 1.0_wp/(hs_k_a*86400.0_wp), " day"
        write(iou,"(a,f9.2,a)") "   1/k_s                      ", 1.0_wp/(hs_k_s*86400.0_wp), " day"
        write(iou,"(a,f9.2,a)") "   1/k_f                      ", 1.0_wp/(hs_k_f*86400.0_wp), " day"
        write(iou,*) ""

        return

    end subroutine aeros_hs_print

end module aeros_held_suarez
