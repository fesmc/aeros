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

    use aeros_defs,     only : dp, wp, io_unit_err, R_d, cp_d, grav, aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure, &
                               aeros_hydrostatic
    use nml,            only : nml_read

    implicit none

    private

    type, public :: aeros_vdiff_class
        logical :: enabled = .FALSE.

        integer :: nlon = 0, nlat = 0

        ! Boundary-layer scheme. .TRUE. (default) = the Richardson-number
        ! diagnosed depth of Frierson (2006) -- the mixed layer is the contiguous
        ! stack of layers up from the surface whose bulk Richardson number is
        ! below `ri_crit`, and K follows the Frierson K-profile that vanishes at
        ! the diagnosed top (this is SpeedyWeather's BulkRichardsonDiffusion).
        ! .FALSE. = the legacy fixed-depth scheme: a single K0 tapered linearly in
        ! sigma to zero at `sigma`. Retained so the value of the diagnosis can be
        ! measured against a fixed depth. BOTH diffuse dry static energy
        ! s = cp T + Phi (not T), so the mixing relaxes toward the dry adiabat
        ! (neutral), never toward isothermal (spuriously stable).
        logical  :: richardson = .TRUE.

        ! Fixed-depth knobs (used only when richardson = .FALSE.).
        real(wp) :: k0    = 10.0_wp     ! eddy diffusivity [m2 s-1]
        real(wp) :: sigma = 0.7_wp      ! boundary-layer top [sigma]; K=0 above

        ! Richardson-scheme knobs (Frierson 2006 eq. 12-20; SpeedyWeather defaults).
        real(wp) :: ri_crit    = 10.0_wp     ! critical bulk Richardson number
        real(wp) :: z0         = 3.21e-5_wp  ! roughness length [m]
        real(wp) :: von_karman = 0.4_wp      ! von Karman constant
        real(wp) :: surf_frac  = 0.1_wp      ! surface-layer fraction f_b
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

    subroutine aeros_vdiff_load(vd, filename, grd, defaults_file)
        implicit none
        type(aeros_vdiff_class), intent(inout) :: vd
        character(len=*),        intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_grid_class),  intent(in)    :: grd

        logical  :: enabled, richardson
        real(wp) :: k0, sigma, ri_crit

        enabled = vd%enabled
        k0 = vd%k0; sigma = vd%sigma
        richardson = vd%richardson; ri_crit = vd%ri_crit

        call nml_read(filename, "vdiff", "enabled",    enabled, defaults_file=defaults_file)
        call nml_read(filename, "vdiff", "k0",         k0, defaults_file=defaults_file)
        call nml_read(filename, "vdiff", "sigma",      sigma, defaults_file=defaults_file)
        call nml_read(filename, "vdiff", "richardson", richardson, defaults_file=defaults_file)
        call nml_read(filename, "vdiff", "ri_crit",    ri_crit, defaults_file=defaults_file)

        vd%k0 = k0; vd%sigma = sigma
        vd%richardson = richardson; vd%ri_crit = ri_crit

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

    subroutine aeros_vdiff_apply(vd, vg, t, qv, u, v, lnps_g, dt, cd_surf, u_min)
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
        real(wp), intent(in)    :: cd_surf        ! surface drag coeff (0 = off)
        real(wp), intent(in)    :: u_min          ! wind floor for the drag [m/s]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: phi_full(vg%nlev), phi_half(0:vg%nlev), zf(vg%nlev)
        real(wp) :: cface(vg%nlev)              ! C_{k+1/2}, k=1..nlev-1
        real(wp) :: kfull(vg%nlev)              ! K at full levels [m2 s-1]
        real(wp) :: sub(vg%nlev), dia(vg%nlev), sup(vg%nlev)
        real(wp) :: diam(vg%nlev), col(vg%nlev), out(vg%nlev)
        real(wp) :: ps, rmk, rhof, thalf, dz, ramp, kface, sig
        real(wp) :: rho_s, wind, drag, cp
        integer  :: i, j, k, nlev, kh

        if (.not. vd%enabled) return

        nlev = vg%nlev
        cp   = real(cp_d, wp)

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,phalf,pfull,dpc,phi_full,phi_half,zf,cface,kfull,kh, &
        !$omp           sub,dia,sup,diam,col,out,ps,rmk,rhof,thalf,dz,ramp,kface,sig, &
        !$omp           rho_s,wind,drag)
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

                ! Full-level diffusivity K(k) [m2 s-1]. Richardson: the Frierson
                ! (2006) K-profile over the bulk-Ri-diagnosed mixed layer. Fixed:
                ! K0 tapered linearly in sigma to zero at vd%sigma.
                if (vd%richardson) then
                    call richardson_kprofile(nlev, t(i,j,:), qv(i,j,:), &
                            u(i,j,:), v(i,j,:), phi_full, zf, vd%ri_crit, &
                            vd%z0, vd%von_karman, vd%surf_frac, cp, kfull, kh)
                else
                    do k = 1, nlev
                        sig  = pfull(k)/ps
                        ramp = (sig - vd%sigma)/(1.0_wp - vd%sigma)
                        ramp = max(0.0_wp, min(1.0_wp, ramp))
                        kfull(k) = vd%k0*ramp
                    end do
                end if

                ! interface conductances C_{k+1/2} = rho K_face / dz, k=1..nlev-1;
                ! zero-flux top (no C above layer 1) and surface (none below nlev).
                ! K_face is the mean of the two adjacent full-level K, so an
                ! interface at the mixed-layer top (K=0 just above) carries no flux.
                do k = 1, nlev-1
                    kface = 0.5_wp*(kfull(k) + kfull(k+1))
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

                ! Heat: diffuse DRY STATIC ENERGY s = cp T + Phi, not T. Mixing s
                ! relaxes the column toward the dry adiabat (neutral); mixing T
                ! would relax toward isothermal (spuriously stable) and build a
                ! warm, convection-choking lower troposphere. Phi is held fixed
                ! across the solve, so conserving int s dp conserves the enthalpy
                ! int cp T dp exactly. Convert back with the same Phi.
                do k = 1, nlev
                    col(k) = cp*t(i,j,k) + phi_full(k)
                end do
                call tridiag(sub, dia, sup, col, out, nlev)
                do k = 1, nlev
                    t(i,j,k) = (out(k) - phi_full(k))/cp
                end do

                ! Moisture: zero-flux surface, the shared matrix.
                col = qv(i,j,:); call tridiag(sub, dia, sup, col, out, nlev)
                qv(i,j,:) = out          ! positivity-preserving for this M-matrix

                ! Momentum: same matrix, but the surface is NOT zero-flux -- a
                ! bulk stress tau = rho c_d |u| u drains the lowest layer. As an
                ! implicit (backward-Euler) sink, tau = rho c_d |u^n| u^{n+1}, it
                ! adds rmk_surf rho_s c_d |u^n| to the k=nlev diagonal (|u| lagged,
                ! isotropic so u and v share it). Unconditionally stable, and it is
                ! the momentum boundary flux the heat/moisture fluxes already have.
                diam = dia
                if (cd_surf > 0.0_wp) then
                    wind  = max(u_min, sqrt(u(i,j,nlev)**2 + v(i,j,nlev)**2))
                    rho_s = phalf(nlev)/(R_d*t(i,j,nlev))
                    drag  = (dt*grav/dpc(nlev))*rho_s*cd_surf*wind
                    diam(nlev) = dia(nlev) + drag
                end if
                col = u(i,j,:);  call tridiag(sub, diam, sup, col, out, nlev)
                u(i,j,:)  = out
                col = v(i,j,:);  call tridiag(sub, diam, sup, col, out, nlev)
                v(i,j,:)  = out
            end do
        end do
        !$omp end parallel do

        return
    end subroutine aeros_vdiff_apply

    subroutine richardson_kprofile(nlev, tcol, qcol, ucol, vcol, phi, zf, &
                                   ri_crit, z0, kappa, fb, cp, kfull, kh)
        ! Frierson (2006) boundary-layer K-profile with the bulk Richardson-number
        ! mixed-layer depth of Frierson (2007) -- SpeedyWeather's
        ! BulkRichardsonDiffusion. The mixed layer is the contiguous stack of
        ! layers up from the surface with bulk Ri < ri_crit; K is zero above it.
        !
        !   Ri(k) = Phi_k (s_v,k - s_v,sfc) / (s_v,sfc |U_k|^2),  s_v = cp T_v + Phi
        !   K0    = kappa |U_sfc| (kappa/ln(Z/z0))(1 - Ri_sfc/ri_crit)   (~ kappa u*)
        !   K(z)  = K0 min(z, fb h) * [z<fb h ? 1 : zfac] * [Ri(kh)>0 ? Rifac : 1]
        ! with h the mixed-layer-top height, fb the surface-layer fraction, and
        ! zfac (eq. 18) / Rifac (eq. 20) the Frierson (2006) shape factors. Ri uses
        ! virtual dry static energy s_v = cp T_v + Phi referenced to the surface.

        implicit none
        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: tcol(:), qcol(:), ucol(:), vcol(:)
        real(wp), intent(in)  :: phi(:), zf(:)
        real(wp), intent(in)  :: ri_crit, z0, kappa, fb, cp
        real(wp), intent(out) :: kfull(:)
        integer,  intent(out) :: kh

        real(wp), parameter :: v2min = 1.0e-4_wp   ! |U|^2 floor [m2 s-2]
        real(wp) :: ri(nlev), tv, sv, v2, theta0, theta1, ri_n, sqrtc
        real(wp) :: logzz0, usfc, k0, h, z, zm, kk, rr
        integer  :: k

        ! bulk Richardson number of each level relative to the surface (k=nlev)
        tv       = tcol(nlev)*(1.0_wp + 0.608_wp*qcol(nlev))
        theta0   = cp*tv
        theta1   = theta0 + phi(nlev)
        v2       = max(ucol(nlev)**2 + vcol(nlev)**2, v2min)
        ri(nlev) = phi(nlev)*(theta1 - theta0)/(theta0*v2)
        do k = 1, nlev-1
            tv    = tcol(k)*(1.0_wp + 0.608_wp*qcol(k))
            sv    = cp*tv + phi(k)
            v2    = max(ucol(k)**2 + vcol(k)**2, v2min)
            ri(k) = phi(k)*(sv - theta1)/(theta1*v2)
        end do

        ! mixed-layer top: uppermost layer reached by a contiguous run of
        ! Ri < ri_crit up from the surface
        kh = nlev
        do while (kh > 0)
            if (ri(kh) < ri_crit) then
                kh = kh - 1
            else
                exit
            end if
        end do
        kh = kh + 1

        kfull = 0.0_wp
        if (kh > nlev) return              ! surface layer already stable: no BL

        h      = max(zf(kh), z0)
        logzz0 = log(max(zf(nlev), z0)/z0)
        ri_n   = min(max(ri(nlev), 0.0_wp), ri_crit)
        sqrtc  = (kappa/logzz0)*(1.0_wp - ri_n/ri_crit)
        usfc   = sqrt(max(ucol(nlev)**2 + vcol(nlev)**2, v2min))
        k0     = kappa*usfc*sqrtc

        do k = kh, nlev
            z  = max(zf(k), z0)
            zm = min(z, fb*h)
            kk = k0*zm
            if (z >= fb*h) &
                kk = kk*(z/(fb*h))*(1.0_wp - (z - fb*h)/((1.0_wp - fb)*h))**2
            if (ri(kh) > 0.0_wp) then
                rr = ri(kh)/ri_crit
                kk = kk/(1.0_wp + rr*logzz0/(1.0_wp - rr))
            end if
            kfull(k) = max(kk, 0.0_wp)
        end do

        return
    end subroutine richardson_kprofile

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
        write(io_unit, '(a,l1)') "    diffuses dry static energy s = cp T + Phi"
        if (vd%richardson) then
            write(io_unit, '(a)')      "    BL depth = Richardson-diagnosed (Frierson 2006)"
            write(io_unit, '(a,f8.2)') "    ri_crit = ", vd%ri_crit
        else
            write(io_unit, '(a)')      "    BL depth = fixed"
            write(io_unit, '(a,f8.2,a)') "    k0      = ", vd%k0, " m2/s"
            write(io_unit, '(a,f8.3)')   "    sigma   = ", vd%sigma
        end if
        return
    end subroutine aeros_vdiff_report

end module aeros_vdiff
