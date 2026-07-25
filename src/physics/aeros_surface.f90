module aeros_surface
    ! The surface energy and moisture budget at the grid seam. This is what
    ! closes the atmosphere's budget once Held-Suarez relaxation is removed:
    ! radiation only cools, so without a surface source the column has nothing
    ! to balance it against and just cools secularly. The surface supplies that
    ! source as turbulent sensible and latent fluxes out of a fixed sea surface.
    !
    ! === Minimal by design: prescribed-SST aquaplanet =======================
    !
    ! The surface temperature is PRESCRIBED, not solved. A fixed SST is an
    ! infinite heat reservoir: the atmosphere equilibrates to it through the
    ! surface fluxes, and there is no surface energy balance to close. This is
    ! the canonical aquaplanet setup (Neale & Hoskins 2001) and the smallest
    ! thing that lets radiative-convective equilibrium exist. A slab ocean with
    ! a real surface energy balance -- where SST responds to the net surface
    ! flux -- is a later milestone (design.md section 6.1).
    !
    ! Default SST is the APE "control" profile:
    !   T_s(phi) = 273.15 + 27 (1 - sin^2(1.5 phi))   for |phi| < 60 deg
    !   T_s(phi) = 273.15                             poleward of 60 deg
    ! 27 C at the equator to freezing at 60 deg. Swap for an ERA5 zonal-mean SST
    ! once that lands (t_s is a public field; only the init needs changing).
    !
    ! === The fluxes =========================================================
    !
    ! Bulk aerodynamic formulae with fixed exchange coefficients and a wind
    ! floor (gustiness), evaluated against the lowest model layer:
    !   SH = rho1 cp C_H |U| (T_s - T_1)                    [W m-2]
    !   E  = rho1    C_E |U| (q_sat(T_s,p_s) - q_1)         [kg m-2 s-1]
    ! Sensible heat warms the lowest layer through wrk%dtdt (centered path, like
    ! the Held-Suarez forcing it replaces: single-signed, smooth in time).
    ! Evaporation moistens the lowest layer as a forward increment to qv_g, the
    ! same treatment convection and transport give the gridpoint humidity.
    !
    ! === What is deferred ===================================================
    !
    ! No boundary-layer diffusion. The surface fluxes go into the lowest layer
    ! and CONVECTION carries them up the column -- the standard idealized-RCE
    ! closure. A dedicated vertical diffusion is a later addition if the
    ! boundary layer needs it; it is noted, not smuggled in.

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

        ! Prescribed SST amplitude and cap latitude (APE control defaults).
        real(wp) :: sst_eq  = 27.0_wp      ! equator-minus-freezing SST [K]
        real(wp) :: sst_lat = 60.0_wp      ! poleward of this SST is frozen [deg]

        ! Prescribed sea surface temperature [K], (nlon,nlat). Public so
        ! radiation can read it as the skin temperature.
        real(wp), allocatable :: t_s(:,:)

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

        integer  :: i, j
        real(wp) :: latr, dsst

        call aeros_surface_end(surf)

        surf%enabled = enabled
        surf%nlon    = grd%nlon
        surf%nlat    = grd%nlat

        allocate(surf%t_s(grd%nlon, grd%nlat))
        allocate(surf%shf(grd%nlon, grd%nlat))
        allocate(surf%lhf(grd%nlon, grd%nlat))
        allocate(surf%evap(grd%nlon, grd%nlat))

        ! APE control profile, zonally symmetric.
        do j = 1, grd%nlat
            latr = grd%lat(j)*pi/180.0_wp
            if (abs(grd%lat(j)) < surf%sst_lat) then
                dsst = surf%sst_eq*(1.0_wp - sin(1.5_wp*latr)**2)
            else
                dsst = 0.0_wp
            end if
            do i = 1, grd%nlon
                surf%t_s(i,j) = T0 + dsst
            end do
        end do

        surf%shf = 0.0_wp; surf%lhf = 0.0_wp; surf%evap = 0.0_wp

        return
    end subroutine aeros_surface_init

    subroutine aeros_surface_load(surf, filename, grd)
        ! Configure from the `surface` namelist group, then init.

        implicit none
        type(aeros_surf_class), intent(inout) :: surf
        character(len=*),       intent(in)    :: filename
        type(aeros_grid_class), intent(in)    :: grd

        logical  :: enabled
        real(wp) :: c_h, c_e, u_min, sst_eq, sst_lat

        enabled = surf%enabled
        c_h = surf%c_h; c_e = surf%c_e; u_min = surf%u_min
        sst_eq = surf%sst_eq; sst_lat = surf%sst_lat

        call nml_read(filename, "surface", "enabled", enabled)
        call nml_read(filename, "surface", "c_h",     c_h)
        call nml_read(filename, "surface", "c_e",     c_e)
        call nml_read(filename, "surface", "u_min",   u_min)
        call nml_read(filename, "surface", "sst_eq",  sst_eq)
        call nml_read(filename, "surface", "sst_lat", sst_lat)

        surf%c_h = c_h; surf%c_e = c_e; surf%u_min = u_min
        surf%sst_eq = sst_eq; surf%sst_lat = sst_lat

        call aeros_surface_init(surf, grd, enabled)

        return
    end subroutine aeros_surface_load

    subroutine aeros_surface_end(surf)
        implicit none
        type(aeros_surf_class), intent(inout) :: surf
        if (allocated(surf%t_s))  deallocate(surf%t_s)
        if (allocated(surf%shf))  deallocate(surf%shf)
        if (allocated(surf%lhf))  deallocate(surf%lhf)
        if (allocated(surf%evap)) deallocate(surf%evap)
        surf%enabled = .FALSE.
        surf%nlon = 0; surf%nlat = 0
        return
    end subroutine aeros_surface_end

    subroutine aeros_surface_apply(surf, vg, t_g, qv_g, lnps_g, u, v, dtdt, dt)
        ! Sensible and latent surface fluxes into the lowest model layer.
        !
        ! t_g, lnps_g, u, v are the current gridpoint fields (aeros_tendency's
        ! wrk); qv_g is the gridpoint humidity, moistened in place by
        ! evaporation; dtdt is the grid temperature tendency the sensible
        ! heating is added to. The lowest layer is k = vg%nlev.

        implicit none
        type(aeros_surf_class),  intent(inout) :: surf
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(inout) :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat) ln[Pa]
        real(wp), intent(in)    :: u(:,:,:)       ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(in)    :: v(:,:,:)       ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(inout) :: dtdt(:,:,:)    ! (nlon,nlat,nlev) [K s-1]
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
                sh = rho1*cp_d*surf%c_h*wind*(surf%t_s(i,j) - t1)

                ! evaporation [kg m-2 s-1] from the saturation deficit at the SST
                call aeros_qsat(surf%t_s(i,j), ps, qs_s, dqs)
                ev = rho1*surf%c_e*wind*(qs_s - q1)

                ! sensible heating of the lowest layer, centered path
                dtdt(i,j,ks) = dtdt(i,j,ks) + sh*grav/(cp_d*dpc(ks))

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
        write(io_unit, '(a,f8.2,a,f8.2,a)') "    SST     = ", &
            minval(surf%t_s), " to ", maxval(surf%t_s), " K"
        write(io_unit, '(a,es10.2,a,es10.2)') "    C_H/C_E = ", surf%c_h, " / ", surf%c_e
        return
    end subroutine aeros_surface_report

end module aeros_surface
