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
    ! Both are applied by aeros_vdiff as the IMPLICIT lower boundary of its
    ! column solve (sensible heat in dry static energy, so the flux is driven
    ! by the potential-temperature difference; evaporation in q), exactly like
    ! the momentum stress: this module only diagnoses the conductances. The
    ! implicit form is what makes the lowest layer's response independent of
    ! its thickness and of the time step -- the surface-layer treatment that
    ! keeps the fluxes resolution-insensitive -- and the fluxes reported (shf,
    ! lhf, evap) are exactly what the column received.
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
    real(wp), parameter :: VON_KARMAN    = 0.4_wp

    type, public :: aeros_surf_class
        logical :: enabled = .FALSE.

        integer :: nlon = 0, nlat = 0

        ! Bulk-flux parameters.
        ! Exchange coefficients. <= 0 (the default) = the bulk-Richardson
        ! log-law coefficient diagnosed per point (below, SpeedyWeather's
        ! BulkRichardsonDrag / Frierson 2006): C = kappa^2 / ln^2(z_1/z_0) *
        ! (1 - Ri/Ri_c)^2, one coefficient for heat, moisture and momentum. A
        ! positive value pins a fixed coefficient instead (intercomparison,
        ! tests). z_1 is the lowest-level height, z_0 the roughness length.
        real(wp) :: c_h   = -1.0_wp        ! sensible-heat exchange coefficient
        real(wp) :: c_e   = -1.0_wp        ! moisture exchange coefficient
        real(wp) :: z0_ocean = 1.0e-4_wp   ! roughness length over sea [m]
        real(wp) :: z0_land  = 0.5_wp      ! roughness length over land [m]
        real(wp) :: ri_crit  = 10.0_wp     ! critical bulk Richardson number
        real(wp) :: cd_min   = 1.0e-5_wp   ! floor of the diagnosed coefficient
        ! Surface wind the fluxes see (SpeedyWeather's SurfaceCondition):
        !   |V_0| = sqrt( (wind_slowdown |V_1|)^2 + gust_speed^2 )
        ! the lowest-level wind slowed toward the surface plus sub-grid gusts in
        ! quadrature. The gust term is what keeps the fluxes alive under weak
        ! resolved winds (a 1 m/s floor gave a third of the observed evaporation).
        real(wp) :: gust_speed    = 5.0_wp     ! [m s-1]
        real(wp) :: wind_slowdown = 0.95_wp    ! [-]
        ! Surface momentum drag coefficient. The bulk stress tau = rho C_d |u| u
        ! is the boundary-layer momentum sink; without it a surface-trapped jet
        ! (the RCE's low-level jet) spins up unbraked. Applied IMPLICITLY as the
        ! momentum bottom boundary condition of the vertical diffusion (forward-
        ! splitting momentum on the leapfrog is unstable), so the coefficient is
        ! carried here but the stress is imposed in aeros_vdiff. Zero = no drag
        ! (bit-unchanged default); ~1.5e-3 to enable, like c_h/c_e.
        real(wp) :: c_d   = -1.0_wp        ! momentum drag coefficient (0 = no drag)
        ! The diagnosed (or fixed) momentum coefficient per point, refreshed by
        ! aeros_surface_apply; aeros_vdiff imposes the stress with it.
        real(wp), allocatable :: cdrag(:,:)
        ! The exchange conductances the implicit boundary-layer solve applies
        ! (aeros_vdiff): ch_v = rho C_H |V_0| and ce_v = beta rho C_E |V_0|
        ! [kg m-2 s-1], with the skin temperature and its saturation humidity
        ! the fluxes are driven against. Refreshed by aeros_surface_apply.
        real(wp), allocatable :: ch_v(:,:), ce_v(:,:), skin(:,:), qs_skin(:,:)

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
    public :: aeros_surface_wind
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


        allocate(surf%cdrag(grd%nlon, grd%nlat)); surf%cdrag = 0.0_wp



        allocate(surf%ch_v(grd%nlon, grd%nlat), surf%ce_v(grd%nlon, grd%nlat), &



                 surf%skin(grd%nlon, grd%nlat), surf%qs_skin(grd%nlon, grd%nlat))



        surf%ch_v = 0.0_wp; surf%ce_v = 0.0_wp; surf%skin = 0.0_wp; surf%qs_skin = 0.0_wp
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
        real(wp) :: c_h, c_e, gust, slow, c_d, z0o, z0l, ric

        enabled = surf%enabled
        c_h = surf%c_h; c_e = surf%c_e; gust = surf%gust_speed; slow = surf%wind_slowdown
        c_d = surf%c_d; z0o = surf%z0_ocean; z0l = surf%z0_land; ric = surf%ri_crit

        call nml_read(filename, "surface", "enabled", enabled, defaults_file=defaults_file)
        call nml_read(filename, "surface", "c_h",     c_h, defaults_file=defaults_file)
        call nml_read(filename, "surface", "c_e",     c_e, defaults_file=defaults_file)
        call nml_read(filename, "surface", "gust_speed",    gust, defaults_file=defaults_file)
        call nml_read(filename, "surface", "wind_slowdown", slow, defaults_file=defaults_file)
        call nml_read(filename, "surface", "c_d",     c_d, defaults_file=defaults_file)
        call nml_read(filename, "surface", "z0_ocean", z0o, defaults_file=defaults_file)
        call nml_read(filename, "surface", "z0_land",  z0l, defaults_file=defaults_file)
        call nml_read(filename, "surface", "ri_crit",  ric, defaults_file=defaults_file)

        surf%c_h = c_h; surf%c_e = c_e; surf%gust_speed = gust; surf%wind_slowdown = slow
        surf%c_d = c_d; surf%z0_ocean = z0o; surf%z0_land = z0l; surf%ri_crit = ric

        call aeros_surface_init(surf, grd, enabled)

        return
    end subroutine aeros_surface_load

    subroutine aeros_surface_end(surf)
        implicit none
        type(aeros_surf_class), intent(inout) :: surf
        if (allocated(surf%shf))  deallocate(surf%shf)
        if (allocated(surf%cdrag)) deallocate(surf%cdrag)

        if (allocated(surf%ch_v))  deallocate(surf%ch_v)

        if (allocated(surf%ce_v))  deallocate(surf%ce_v)

        if (allocated(surf%skin))  deallocate(surf%skin)

        if (allocated(surf%qs_skin)) deallocate(surf%qs_skin)
        if (allocated(surf%lhf))  deallocate(surf%lhf)
        if (allocated(surf%evap)) deallocate(surf%evap)
        surf%enabled = .FALSE.
        surf%nlon = 0; surf%nlat = 0
        return
    end subroutine aeros_surface_end

    subroutine aeros_surface_apply(surf, vg, sst, t_g, qv_g, lnps_g, u, v, evap_eff, is_land)
        ! Diagnose the surface exchange for this step: the bulk transfer
        ! coefficient per point (class header), the surface wind, the air
        ! density, and from them the conductances ch_v / ce_v [kg m-2 s-1] plus
        ! the skin temperature and its saturation humidity. Nothing is applied
        ! here: aeros_vdiff imposes the sensible and latent fluxes as the
        ! IMPLICIT lower boundary of its column solve (the momentum stress was
        ! always there), so the lowest layer's response does not depend on that
        ! layer's thickness or on the time step, and the fluxes it reports are
        ! exactly what the column received.
        !
        ! sst is the skin temperature the fluxes are driven against -- the ocean
        ! SST over sea (aeros_ocean), the prognostic soil temperature over land
        ! (aeros_land composes the combined field). t_g, lnps_g, u, v are the
        ! current gridpoint fields; qv_g the gridpoint humidity. The lowest layer
        ! is k = vg%nlev.
        !
        ! evap_eff is an OPTIONAL per-point evaporation efficiency (beta):
        ! potential evaporation is multiplied by it (a drying soil evaporates
        ! below potential). is_land selects the land roughness length.
        implicit none
        type(aeros_surf_class),  intent(inout) :: surf
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: sst(:,:)       ! (nlon,nlat) skin T [K]
        real(wp), intent(in)    :: t_g(:,:,:)     ! (nlon,nlat,nlev) [K]
        real(wp), intent(in)    :: qv_g(:,:,:)    ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: lnps_g(:,:)    ! (nlon,nlat) ln[Pa]
        real(wp), intent(in)    :: u(:,:,:)       ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(in)    :: v(:,:,:)       ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(in), optional :: evap_eff(:,:) ! (nlon,nlat) beta [-]
        logical,  intent(in), optional :: is_land(:,:)  ! (nlon,nlat) roughness selector
        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: ps, t1, q1, wind, rho1, qs_s, dqs, beta
        real(wp) :: z1, z0, cbulk, ri, ch, ce
        logical  :: have_beta, have_land
        integer  :: i, j, ks
        if (.not. surf%enabled) return
        have_beta = present(evap_eff)
        have_land = present(is_land)
        ks = vg%nlev                              ! lowest (surface) layer
        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,phalf,pfull,dpc,ps,t1,q1,wind,rho1,qs_s,dqs,beta,z1,z0,cbulk,ri,ch,ce)
        do j = 1, surf%nlat
            do i = 1, surf%nlon
                ps = exp(lnps_g(i,j))
                call aeros_vgrid_pressure(vg, ps, phalf, pfull, dpc)
                t1 = t_g(i,j,ks)
                q1 = qv_g(i,j,ks)
                ! lowest-layer air density (virtual temperature) and the surface wind
                rho1 = pfull(ks)/(R_d*t1*(1.0_wp + RV_OVER_RD_M1*q1))
                wind = aeros_surface_wind(surf, u(i,j,ks), v(i,j,ks))
                ! Bulk transfer coefficient (class header): the neutral log-law
                ! value for the lowest-level height z_1 (hypsometric, from the
                ! layer's own T) over the roughness z_0, reduced by the bulk
                ! Richardson number of that layer -- SpeedyWeather's form of
                ! Frierson (2006) eqs 12-15, where the surface air is the level-1
                ! air brought down dry-adiabatically, so the static-energy jump
                ! over the layer is g z_1 and Ri = (g z_1)^2 / (cp T_v |V_0|^2).
                z1 = real(R_d,wp)*t1/real(grav,wp)*log(ps/pfull(ks))
                z0 = surf%z0_ocean
                if (have_land) then
                    if (is_land(i,j)) z0 = surf%z0_land
                end if
                z1 = max(z1, z0)
                cbulk = (VON_KARMAN/log(z1/z0))**2
                ri = (real(grav,wp)*z1)**2/(real(cp_d,wp)*t1*(1.0_wp + RV_OVER_RD_M1*q1)*wind*wind)
                ri = min(max(ri, 0.0_wp), surf%ri_crit)
                cbulk = max(surf%cd_min, cbulk*(1.0_wp - ri/surf%ri_crit)**2)
                ch = surf%c_h; if (ch <= 0.0_wp) ch = cbulk
                ce = surf%c_e; if (ce <= 0.0_wp) ce = cbulk
                surf%cdrag(i,j) = surf%c_d
                if (surf%c_d < 0.0_wp) surf%cdrag(i,j) = cbulk
                beta = 1.0_wp
                if (have_beta) beta = evap_eff(i,j)
                call aeros_qsat(sst(i,j), ps, qs_s, dqs)
                surf%ch_v(i,j)    = rho1*ch*wind
                surf%ce_v(i,j)    = beta*rho1*ce*wind
                surf%skin(i,j)    = sst(i,j)
                surf%qs_skin(i,j) = qs_s
            end do
        end do
        !$omp end parallel do
        return
    end subroutine aeros_surface_apply

    pure real(wp) function aeros_surface_wind(surf, u1, v1) result(wind)
        ! The surface wind speed the bulk fluxes (and the momentum drag in
        ! aeros_vdiff) see: lowest-level wind slowed toward the surface plus
        ! sub-grid gusts in quadrature (class header).
        implicit none
        type(aeros_surf_class), intent(in) :: surf
        real(wp), intent(in) :: u1, v1
        wind = sqrt((surf%wind_slowdown**2)*(u1*u1 + v1*v1) + surf%gust_speed**2)
    end function aeros_surface_wind

    subroutine aeros_surface_report(surf, io_unit)
        implicit none
        type(aeros_surf_class), intent(in) :: surf
        integer,                intent(in) :: io_unit

        write(io_unit, '(a)')      "  surface:"
        write(io_unit, '(a,l1)')   "    enabled = ", surf%enabled
        if (.not. surf%enabled) return
        write(io_unit, '(a,es10.2,a,es10.2,a)') "    C_H/C_E = ", surf%c_h, " / ", surf%c_e, &
            "  (<= 0: bulk-Richardson log-law)"
        write(io_unit, '(a,es10.2,a,es10.2,a,f5.1)') "    z0 ocean/land = ", surf%z0_ocean, " / ", &
            surf%z0_land, "  Ri_c = ", surf%ri_crit
        write(io_unit, '(a,f6.2,a,f5.2)') "    gust_speed / wind_slowdown = ", surf%gust_speed, " / ", surf%wind_slowdown
        return
    end subroutine aeros_surface_report

end module aeros_surface
