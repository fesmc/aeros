module aeros_vdiff
    ! Boundary-layer vertical diffusion at the grid seam.
    !
    ! === Why this exists ====================================================
    !
    ! The prescribed-SST surface (aeros_surface) injects sensible heat and
    ! moisture into the LOWEST model layer only. The design bet at M2.4c was
    ! that convection would carry that source up the column -- the standard
    ! idealized-RCE closure, so no boundary-layer diffusion was added. In the
    ! coupled RCE that bet fails: the Simplified Betts-Miller scheme
    ! (aeros_convection) mixes only the cloud layer (level of free convection to
    ! level of neutral buoyancy) and by construction EXCLUDES the sub-cloud
    ! layer, so the surface fluxes pile up in the lowest layer with no vertical
    ! sink. The lowest layer builds a 15-20 K warm/moist inversion over ~1 month,
    ! convection then converts the unstable sub-cloud layer into runaway heating,
    ! and the run blows up in a grid-scale lowest-layer hot spot (m2_handoff.md
    ! task #7; the diagnosis is in m2_results.md). Making horizontal diffusion
    ! stronger does not help (it is not a horizontal grid mode) and a smaller
    ! step only delays it: the missing physics is the vertical mixing.
    !
    ! This module supplies it: a simple down-gradient vertical diffusion of
    ! temperature, humidity and momentum through the boundary layer, which
    ! spreads the surface source up the column so the lowest layer equilibrates
    ! toward the SST instead of running away. It runs BETWEEN the surface fluxes
    ! and convection, so convection acts on the already-mixed column.
    !
    ! === The scheme =========================================================
    !
    ! Down-gradient diffusion with a fixed eddy diffusivity K within a
    ! boundary-layer depth, tapering to zero at its top:
    !
    !   dX/dt = (g/dp_k) [ F_{k-1/2} - F_{k+1/2} ],   F_{k+1/2} = C_{k+1/2}(X_k - X_{k+1})
    !
    ! with the interface conductance C_{k+1/2} = rho K / dz [kg m-2 s-1] and
    ! zero-flux boundaries at the model top AND the surface -- the surface flux
    ! is aeros_surface's job, so vdiff only REDISTRIBUTES what is already in the
    ! column, conserving each field's mass-weighted column integral. Solved
    ! IMPLICITLY (backward Euler, one tridiagonal per column): vertical diffusion
    ! is stiff at these layer thicknesses, so an explicit solve would be
    ! CFL-limited far below the dynamics step. The matrix is the same for all
    ! four fields; only the right-hand side differs.
    !
    ! K is fixed (one namelist knob) with a linear taper in sigma from its full
    ! value at the surface to zero at vdiff_sigma. Not a Richardson-number or TKE
    ! closure: for an idealized aquaplanet the fixed-K boundary layer is the
    ! smallest thing that removes the pile-up, and a stability-dependent K is a
    ! later refinement if the boundary-layer structure needs it.
    !
    ! === Coupling ===========================================================
    !
    ! Unlike the source terms (surface, convection, radiation), vertical
    ! diffusion is NOT forward-split. A diffusion operator applied explicitly on
    ! a leapfrog is unconditionally unstable -- the increment dt*L*X^n added to
    ! X^{n+1} = X^{n-1} + ... gives a root pair whose product is -1, so one root
    ! always has |r| > 1, and it is scale-selective, so it blows up at the grid
    ! scale however small K is. (That is exactly why the horizontal
    ! hyperdiffusion is done implicitly on the stepped state, not through the
    ! tendency.) So vdiff is applied IMPLICITLY on the n+1 state: aeros_timestep
    ! synthesizes the just-stepped temperature, winds and surface pressure to the
    ! grid, calls this to diffuse T, q, u and v in place, and transforms back --
    ! next to the horizontal diffusion and the sponge, all three implicit on the
    ! same state. Humidity is a gridpoint field and is diffused in place directly.

    use aeros_defs,     only : dp, wp, io_unit_err, R_d, grav, aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure, &
                               aeros_hydrostatic
    use nml,            only : nml_read

    implicit none

    private

    type, public :: aeros_vdiff_class
        logical :: enabled = .FALSE.

        integer :: nlon = 0, nlat = 0

        real(wp) :: k0    = 10.0_wp     ! eddy diffusivity [m2 s-1]
        real(wp) :: sigma = 0.7_wp      ! boundary-layer top [sigma]; K=0 above
    end type aeros_vdiff_class

    public :: aeros_vdiff_init
    public :: aeros_vdiff_load
    public :: aeros_vdiff_end
    public :: aeros_vdiff_apply
    public :: aeros_vdiff_report

contains

    subroutine aeros_vdiff_init(vd, grd, enabled)
        implicit none
        type(aeros_vdiff_class), intent(inout) :: vd
        type(aeros_grid_class),  intent(in)    :: grd
        logical,                 intent(in)    :: enabled

        call aeros_vdiff_end(vd)
        vd%enabled = enabled
        vd%nlon    = grd%nlon
        vd%nlat    = grd%nlat
        return
    end subroutine aeros_vdiff_init

    subroutine aeros_vdiff_load(vd, filename, grd)
        implicit none
        type(aeros_vdiff_class), intent(inout) :: vd
        character(len=*),        intent(in)    :: filename
        type(aeros_grid_class),  intent(in)    :: grd

        logical  :: enabled
        real(wp) :: k0, sigma

        enabled = vd%enabled
        k0 = vd%k0; sigma = vd%sigma

        call nml_read(filename, "vdiff", "enabled", enabled)
        call nml_read(filename, "vdiff", "k0",      k0)
        call nml_read(filename, "vdiff", "sigma",   sigma)

        vd%k0 = k0; vd%sigma = sigma

        call aeros_vdiff_init(vd, grd, enabled)
        return
    end subroutine aeros_vdiff_load

    subroutine aeros_vdiff_end(vd)
        implicit none
        type(aeros_vdiff_class), intent(inout) :: vd
        vd%enabled = .FALSE.
        vd%nlon = 0; vd%nlat = 0
        return
    end subroutine aeros_vdiff_end

    subroutine aeros_vdiff_apply(vd, vg, t, qv, u, v, lnps_g, dt)
        ! One implicit vertical-diffusion sweep over every column, diffusing
        ! temperature, humidity and both wind components IN PLACE.
        !
        ! Applied to the n+1 state (see the coupling note in the header): the
        ! caller synthesizes the just-stepped spectral state to the grid, calls
        ! this, and transforms back. Vertical diffusion is a diffusion operator
        ! and MUST be implicit on the stepped state -- forward-splitting it onto
        ! the leapfrog is unconditionally unstable at the grid scale (the same
        ! reason the horizontal hyperdiffusion is implicit), no matter how small
        ! the diffusivity. The backward-Euler tridiagonal here is unconditionally
        ! stable and conserves each field's mass-weighted column integral (zero
        ! flux at the top and the surface).

        implicit none
        type(aeros_vdiff_class), intent(inout) :: vd
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(inout) :: t(:,:,:)       ! (nlon,nlat,nlev) [K]
        real(wp), intent(inout) :: qv(:,:,:)      ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(inout) :: u(:,:,:)       ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(inout) :: v(:,:,:)
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat) ln[Pa]
        real(wp), intent(in)    :: dt             ! [s]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: phi_full(vg%nlev), phi_half(0:vg%nlev), zf(vg%nlev)
        real(wp) :: cface(vg%nlev)              ! C_{k+1/2}, k=1..nlev-1
        real(wp) :: sub(vg%nlev), dia(vg%nlev), sup(vg%nlev)
        real(wp) :: col(vg%nlev), out(vg%nlev)
        real(wp) :: ps, rmk, rhof, thalf, dz, ramp, kface
        integer  :: i, j, k, nlev

        if (.not. vd%enabled) return

        nlev = vg%nlev

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,phalf,pfull,dpc,phi_full,phi_half,zf,cface, &
        !$omp           sub,dia,sup,col,out,ps,rmk,rhof,thalf,dz,ramp,kface)
        do j = 1, vd%nlat
            do i = 1, vd%nlon
                ps = exp(lnps_g(i,j))
                call aeros_vgrid_pressure(vg, ps, phalf, pfull, dpc)

                ! layer-midpoint heights (aquaplanet: surface geopotential 0);
                ! only differences enter, so the reference cancels
                call aeros_hydrostatic(vg, 0.0_wp, t(i,j,:), phalf, phi_full, phi_half)
                do k = 1, nlev
                    zf(k) = phi_full(k)/grav
                end do

                ! interface conductances C_{k+1/2} = rho K / dz, k=1..nlev-1;
                ! zero-flux top (no C above layer 1) and surface (none below nlev)
                do k = 1, nlev-1
                    ramp = (vg%sigma_half(k) - vd%sigma)/(1.0_wp - vd%sigma)
                    ramp = max(0.0_wp, min(1.0_wp, ramp))
                    kface = vd%k0*ramp
                    if (kface <= 0.0_wp) then
                        cface(k) = 0.0_wp
                        cycle
                    end if
                    thalf = 0.5_wp*(t(i,j,k) + t(i,j,k+1))
                    rhof  = phalf(k)/(R_d*thalf)
                    dz    = zf(k) - zf(k+1)
                    cface(k) = rhof*kface/max(dz, 1.0_wp)
                end do

                ! tridiagonal (backward Euler); same matrix for all four fields
                do k = 1, nlev
                    rmk = dt/(dpc(k)/grav)
                    if (k > 1) then
                        sub(k) = -rmk*cface(k-1)
                    else
                        sub(k) = 0.0_wp
                    end if
                    if (k < nlev) then
                        sup(k) = -rmk*cface(k)
                    else
                        sup(k) = 0.0_wp
                    end if
                    dia(k) = 1.0_wp - sub(k) - sup(k)
                end do

                col = t(i,j,:);  call tridiag(sub, dia, sup, col, out, nlev)
                t(i,j,:)  = out
                col = qv(i,j,:); call tridiag(sub, dia, sup, col, out, nlev)
                qv(i,j,:) = out          ! positivity-preserving for this M-matrix
                col = u(i,j,:);  call tridiag(sub, dia, sup, col, out, nlev)
                u(i,j,:)  = out
                col = v(i,j,:);  call tridiag(sub, dia, sup, col, out, nlev)
                v(i,j,:)  = out
            end do
        end do
        !$omp end parallel do

        return
    end subroutine aeros_vdiff_apply

    pure subroutine tridiag(a, b, c, d, x, n)
        ! Thomas algorithm for a tridiagonal system with sub-diagonal a, diagonal
        ! b, super-diagonal c and right-hand side d. Non-destructive in the
        ! inputs (works on local copies), so the same matrix serves every field.
        implicit none
        real(wp), intent(in)  :: a(:), b(:), c(:), d(:)
        real(wp), intent(out) :: x(:)
        integer,  intent(in)  :: n
        real(wp) :: cp(n), dp_(n), m
        integer  :: k

        cp(1)  = c(1)/b(1)
        dp_(1) = d(1)/b(1)
        do k = 2, n
            m      = b(k) - a(k)*cp(k-1)
            cp(k)  = c(k)/m
            dp_(k) = (d(k) - a(k)*dp_(k-1))/m
        end do
        x(n) = dp_(n)
        do k = n-1, 1, -1
            x(k) = dp_(k) - cp(k)*x(k+1)
        end do
        return
    end subroutine tridiag

    subroutine aeros_vdiff_report(vd, io_unit)
        implicit none
        type(aeros_vdiff_class), intent(in) :: vd
        integer,                 intent(in) :: io_unit
        write(io_unit, '(a)')    "  vdiff:"
        write(io_unit, '(a,l1)') "    enabled = ", vd%enabled
        if (.not. vd%enabled) return
        write(io_unit, '(a,f8.2,a)') "    k0      = ", vd%k0, " m2/s"
        write(io_unit, '(a,f8.3)')   "    sigma   = ", vd%sigma
        return
    end subroutine aeros_vdiff_report

end module aeros_vdiff
