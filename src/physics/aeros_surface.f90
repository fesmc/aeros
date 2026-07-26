module aeros_surface
    ! The surface energy and moisture budget at the grid seam. This is what
    ! closes the atmosphere's budget once Held-Suarez relaxation is removed:
    ! radiation only cools, so without a surface source the column has nothing
    ! to balance it against and just cools secularly. The surface supplies that
    ! source as turbulent sensible and latent fluxes out of the sea surface
    ! (whose temperature is aeros_ocean's, prescribed or a slab).
    !
    ! === The sea surface temperature comes from aeros_ocean =================
    !
    ! The SST the fluxes are driven against is owned by aeros_ocean and passed in
    ! -- either prescribed (a fixed reservoir, the canonical aquaplanet control,
    ! Neale & Hoskins 2001) or a slab whose temperature responds to the net
    ! surface flux. This module does not know or care which; it only reads the
    ! SST field. (Through M2.4 the SST lived here as a prescribed field; M2.5d
    ! moved it to aeros_ocean when the slab was added.)
    !
    ! === The fluxes =========================================================
    !
    ! Bulk aerodynamic formulae with fixed exchange coefficients and a wind
    ! floor (gustiness), evaluated against the lowest model layer:
    !   SH = rho1 cp C_H |U| (SST - T_1)                    [W m-2]
    !   E  = rho1    C_E |U| (q_sat(SST,p_s) - q_1)         [kg m-2 s-1]
    ! Sensible heat warms the lowest layer through wrk%dt_phys, the FORWARD-split
    ! path convection and radiation use: the flux goes into one layer only, so
    ! its sharp vertical structure would excite the computational mode on the
    ! centered leapfrog (the Held-Suarez forcing it replaces was distributed over
    ! the column, and was safe centered; a single-layer flux is not).
    ! Evaporation moistens the lowest layer as a forward increment to qv_g, the
    ! same treatment convection and transport give the gridpoint humidity.
    !
    ! The surface fluxes go into the lowest layer; aeros_vdiff (M2.5b) then mixes
    ! them up the column. Radiation reads the same SST as the skin temperature,
    ! and in slab mode aeros_ocean closes the surface energy balance from these
    ! fluxes plus the surface radiative fluxes.

    use aeros_defs,     only : dp, wp, io_unit_err, R_d, R_v, cp_d, grav, L_v, &
                               T0, pi, aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure
    use aeros_condensation, only : aeros_qsat
    use nml,            only : nml_read

    implicit none

    private

    real(wp), parameter :: RV_OVER_RD_M1 = R_v/R_d - 1.0_wp   ! ~0.608, virtual T

    type, public :: aeros_surf_class
        logical :: enabled = .FALSE.

        integer :: nlon = 0, nlat = 0

        ! Bulk-flux parameters.
        real(wp) :: c_h   = 1.5e-3_wp      ! sensible-heat exchange coefficient
        real(wp) :: c_e   = 1.5e-3_wp      ! moisture exchange coefficient
        real(wp) :: u_min = 1.0_wp         ! wind floor / gustiness [m s-1]
        ! Surface momentum drag coefficient. The bulk stress tau = rho C_d |u| u
        ! is the boundary-layer momentum sink; without it a surface-trapped jet
        ! (the RCE's low-level jet) spins up unbraked. Applied IMPLICITLY as the
        ! momentum bottom boundary condition of the vertical diffusion (forward-
        ! splitting momentum on the leapfrog is unstable), so the coefficient is
        ! carried here but the stress is imposed in aeros_vdiff. Zero = no drag
        ! (bit-unchanged default); ~1.5e-3 to enable, like c_h/c_e.
        real(wp) :: c_d   = 0.0_wp         ! momentum drag coefficient

        ! The sea surface temperature the fluxes are evaluated against is owned by
        ! aeros_ocean and passed to aeros_surface_apply -- prescribed or a slab.

        ! Diagnostics from the last apply, [W m-2] and [kg m-2 s-1].
        real(wp), allocatable :: shf(:,:)     ! sensible heat flux
        real(wp), allocatable :: lhf(:,:)     ! latent heat flux
        real(wp), allocatable :: evap(:,:)    ! evaporation
    end type aeros_surf_class

    public :: aeros_surface_init
    public :: aeros_surface_load
    public :: aeros_surface_end
    public :: aeros_surface_apply
    public :: aeros_surface_report

contains

    subroutine aeros_surface_init(surf, grd, enabled)
        ! Allocate the fields and build the prescribed SST from latitude.

        implicit none
        type(aeros_surf_class), intent(inout) :: surf
        type(aeros_grid_class), intent(in)    :: grd
        logical,                intent(in)    :: enabled

        call aeros_surface_end(surf)

        surf%enabled = enabled
        surf%nlon    = grd%nlon
        surf%nlat    = grd%nlat

        allocate(surf%shf(grd%nlon, grd%nlat))
        allocate(surf%lhf(grd%nlon, grd%nlat))
        allocate(surf%evap(grd%nlon, grd%nlat))

        surf%shf = 0.0_wp; surf%lhf = 0.0_wp; surf%evap = 0.0_wp

        return
    end subroutine aeros_surface_init

    subroutine aeros_surface_load(surf, filename, grd, defaults_file)
        ! Configure from the `surface` namelist group, then init.

        implicit none
        type(aeros_surf_class), intent(inout) :: surf
        character(len=*),       intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_grid_class), intent(in)    :: grd

        logical  :: enabled
        real(wp) :: c_h, c_e, u_min, c_d

        enabled = surf%enabled
        c_h = surf%c_h; c_e = surf%c_e; u_min = surf%u_min; c_d = surf%c_d

        call nml_read(filename, "surface", "enabled", enabled, defaults_file=defaults_file)
        call nml_read(filename, "surface", "c_h",     c_h, defaults_file=defaults_file)
        call nml_read(filename, "surface", "c_e",     c_e, defaults_file=defaults_file)
        call nml_read(filename, "surface", "u_min",   u_min, defaults_file=defaults_file)
        call nml_read(filename, "surface", "c_d",     c_d, defaults_file=defaults_file)

        surf%c_h = c_h; surf%c_e = c_e; surf%u_min = u_min; surf%c_d = c_d

        call aeros_surface_init(surf, grd, enabled)

        return
    end subroutine aeros_surface_load

    subroutine aeros_surface_end(surf)
        implicit none
        type(aeros_surf_class), intent(inout) :: surf
        if (allocated(surf%shf))  deallocate(surf%shf)
        if (allocated(surf%lhf))  deallocate(surf%lhf)
        if (allocated(surf%evap)) deallocate(surf%evap)
        surf%enabled = .FALSE.
        surf%nlon = 0; surf%nlat = 0
        return
    end subroutine aeros_surface_end

    subroutine aeros_surface_apply(surf, vg, sst, t_g, qv_g, lnps_g, u, v, dt_phys, dt)
        ! Sensible and latent surface fluxes into the lowest model layer.
        !
        ! sst is the sea surface temperature (aeros_ocean) the fluxes are driven
        ! against. t_g, lnps_g, u, v are the current gridpoint fields
        ! (aeros_tendency's wrk); qv_g is the gridpoint humidity, moistened in
        ! place by evaporation; dt_phys is the forward-split temperature
        ! increment [K] the sensible heating is added to. The lowest layer is
        ! k = vg%nlev.

        implicit none
        type(aeros_surf_class),  intent(inout) :: surf
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: sst(:,:)       ! (nlon,nlat) sea surface T [K]
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(inout) :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat) ln[Pa]
        real(wp), intent(in)    :: u(:,:,:)       ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(in)    :: v(:,:,:)       ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(inout) :: dt_phys(:,:,:) ! (nlon,nlat,nlev) [K] increment
        real(wp), intent(in)    :: dt             ! [s]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: ps, t1, q1, wind, rho1, qs_s, dqs, sh, ev, dq
        integer  :: i, j, ks

        if (.not. surf%enabled) return

        ks = vg%nlev                              ! lowest (surface) layer

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,phalf,pfull,dpc,ps,t1,q1,wind,rho1,qs_s,dqs,sh,ev,dq)
        do j = 1, surf%nlat
            do i = 1, surf%nlon
                ps = exp(lnps_g(i,j))
                call aeros_vgrid_pressure(vg, ps, phalf, pfull, dpc)

                t1 = t_g(i,j,ks)
                q1 = qv_g(i,j,ks)

                ! lowest-layer air density (virtual temperature) and wind speed
                rho1 = pfull(ks)/(R_d*t1*(1.0_wp + RV_OVER_RD_M1*q1))
                wind = max(surf%u_min, sqrt(u(i,j,ks)**2 + v(i,j,ks)**2))

                ! sensible heat flux [W m-2], positive up (into the atmosphere)
                sh = rho1*cp_d*surf%c_h*wind*(sst(i,j) - t1)

                ! evaporation [kg m-2 s-1] from the saturation deficit at the SST
                call aeros_qsat(sst(i,j), ps, qs_s, dqs)
                ev = rho1*surf%c_e*wind*(qs_s - q1)

                ! sensible heating of the lowest layer, forward-split increment
                ! [K]: it is confined to one layer, so like convection and
                ! radiation its sharp vertical structure would go
                ! computational-mode unstable on the centered leapfrog.
                dt_phys(i,j,ks) = dt_phys(i,j,ks) + sh*grav/(cp_d*dpc(ks))*dt

                ! evaporative moistening, forward increment; never drive q < 0
                ! when the flux is downward (dew onto a colder surface).
                dq = ev*grav/dpc(ks)*dt
                if (qv_g(i,j,ks) + dq < 0.0_wp) then
                    dq = -qv_g(i,j,ks)
                    ev = dq*dpc(ks)/(grav*dt)
                end if
                qv_g(i,j,ks) = qv_g(i,j,ks) + dq

                surf%shf(i,j)  = sh
                surf%lhf(i,j)  = L_v*ev
                surf%evap(i,j) = ev
            end do
        end do
        !$omp end parallel do

        return
    end subroutine aeros_surface_apply

    subroutine aeros_surface_report(surf, io_unit)
        implicit none
        type(aeros_surf_class), intent(in) :: surf
        integer,                intent(in) :: io_unit

        write(io_unit, '(a)')      "  surface:"
        write(io_unit, '(a,l1)')   "    enabled = ", surf%enabled
        if (.not. surf%enabled) return
        write(io_unit, '(a,es10.2,a,es10.2)') "    C_H/C_E = ", surf%c_h, " / ", surf%c_e
        return
    end subroutine aeros_surface_report

end module aeros_surface
