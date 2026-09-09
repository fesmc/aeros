module aeros_bd
    ! Prescribed Brewer-Dobson residual circulation in the lid layer.
    !
    ! The top layer (10-257 hPa at L6) holds the model's whole stratosphere. A
    ! one-layer stratosphere with a 10 hPa top cannot carry the wave-driven
    ! residual circulation that cools the tropical lower stratosphere by
    ! upwelling and warms the polar one by downwelling; without it the lid sits
    ! 5 K warm in the tropics and 20-30 K cold at the poles (ERA5 layer mean),
    ! and that contrast drives the jet to 52 m/s by thermal wind (ERA5 34) and
    ! warms 350 hPa by 2.5 K through the lid's downward longwave (2026-09-08).
    !
    ! This module imposes the missing adiabatic term only:
    !     dT/dt|_BD = -w* S_ref                                  [K/s]
    ! with a prescribed residual vertical velocity w*(lat, month) [m/s, > 0
    ! upward] and the reference lid stability S_ref = dT/dz + g/cp (lid vs the
    ! layer below) of the run the climatology was derived from, both read from
    ! the file. The lid temperature and its radiation stay interactive -- CO2
    ! cooling of the lid can still emerge; only the circulation is fixed. (The
    ! model's own lid stability is NOT used: measured against 365 hPa it is
    ! dominated by the lid temperature itself, so a downwelling that warms the
    ! lid raises its own heating -- the polar lid overshot by 7 K and the jet
    ! fell to 26 m/s in that form, 2026-09-08.)
    ! The climatology is closed (area-weighted mean w* = 0 per month) so it
    ! carries no net mass flux. Two knobs give the prescribed circulation some
    ! climate dependence without inventing physics:
    !   scale   -- a global amplitude factor (1 = the file), and
    !   per_k   -- a fractional change of w* per kelvin of the running annual
    !              global-mean lowest-level temperature anomaly from t_ref
    !              (stratosphere-resolving models: +2..3 %/K under warming;
    !              0 = off, the default).
    ! The seasonal cycle is tied to the calendar month, like the seasonal SST.
    !
    ! Applied at the grid seam as a forward-split increment on the top layer
    ! only, like the other column physics. The file w* was backed out of the
    ! L6 run relaxed toward the ERA5 layer-mean lid (input/bd_wstar_era5.nc,
    ! 2026-09-08): an EFFECTIVE residual velocity for this grid, 1-3 mm/s,
    ! larger than the 70-hPa textbook 0.3 mm/s because the lid is thick and its
    ! stability is measured against 365 hPa. Off unless a file is given.

    use aeros_defs,     only : wp, io_unit_err, aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class
    use ncio,           only : nc_read, nc_size
    use nml,            only : nml_read

    implicit none

    private

    type, public :: aeros_bd_class
        logical :: enabled = .FALSE.
        integer :: nlon = 0, nlat = 0
        real(wp), allocatable :: wstar(:,:)   ! (nlat, 12) residual velocity [m/s], > 0 up
        real(wp), allocatable :: stab(:,:)    ! (nlat, 12) reference lid stability [K/m]
        real(wp) :: scale    = 1.0_wp         ! amplitude factor on w*
        real(wp) :: per_k    = 0.0_wp         ! fractional change of w* per K of tmean - t_ref
        real(wp) :: t_ref    = 285.0_wp       ! reference global-mean lowest-level T [K]
        real(wp) :: tau_mean = 365.0_wp*86400.0_wp  ! running-mean time scale [s]
        real(wp) :: tmean    = 285.0_wp       ! running global-mean lowest-level T [K]
        real(wp), allocatable :: q(:,:)       ! last heating rate of the lid [K/s]
    end type aeros_bd_class

    public :: aeros_bd_init
    public :: aeros_bd_load
    public :: aeros_bd_read
    public :: aeros_bd_end
    public :: aeros_bd_apply
    public :: aeros_bd_factor

contains

    subroutine aeros_bd_init(bd, grd, enabled)
        implicit none
        type(aeros_bd_class),   intent(inout) :: bd
        type(aeros_grid_class), intent(in)    :: grd
        logical,                intent(in)    :: enabled
        call aeros_bd_end(bd)
        bd%enabled = enabled
        bd%nlon = grd%nlon; bd%nlat = grd%nlat
        allocate(bd%wstar(grd%nlat, 12), bd%stab(grd%nlat, 12), bd%q(grd%nlon, grd%nlat))
        bd%wstar = 0.0_wp; bd%stab = 0.0_wp; bd%q = 0.0_wp
        bd%tmean = bd%t_ref
        return
    end subroutine aeros_bd_init

    subroutine aeros_bd_load(bd, filename, grd, defaults_file)
        ! Read the `bd` namelist group, then init. Mirrors aeros_land_load.
        ! If `file` is non-empty and `enabled`, also reads the w*(lat,month)
        ! climatology from that netCDF file.
        implicit none
        type(aeros_bd_class),   intent(inout) :: bd
        character(len=*),       intent(in)    :: filename
        type(aeros_grid_class), intent(in)    :: grd
        character(len=*), intent(in), optional :: defaults_file

        logical            :: enabled
        character(len=512) :: file
        real(wp)           :: scale, per_k, t_ref, tau_mean_days

        enabled       = bd%enabled
        file          = ""
        scale         = bd%scale
        per_k         = bd%per_k
        t_ref         = bd%t_ref
        tau_mean_days = bd%tau_mean/86400.0_wp

        call nml_read(filename, "bd", "enabled",       enabled,       defaults_file=defaults_file)
        call nml_read(filename, "bd", "file",          file,          defaults_file=defaults_file)
        call nml_read(filename, "bd", "scale",         scale,         defaults_file=defaults_file)
        call nml_read(filename, "bd", "per_k",         per_k,         defaults_file=defaults_file)
        call nml_read(filename, "bd", "t_ref",         t_ref,         defaults_file=defaults_file)
        call nml_read(filename, "bd", "tau_mean_days", tau_mean_days, defaults_file=defaults_file)

        ! An empty file forces the module off (no climatology to apply).
        if (len_trim(file) == 0) enabled = .FALSE.

        call aeros_bd_init(bd, grd, enabled)
        bd%scale    = scale
        bd%per_k    = per_k
        bd%t_ref    = t_ref
        bd%tau_mean = tau_mean_days*86400.0_wp
        bd%tmean    = t_ref

        if (enabled) call aeros_bd_read(bd, trim(file), grd)

        return
    end subroutine aeros_bd_load

    subroutine aeros_bd_read(bd, filename, grd)
        ! Read wstar and stab (month, lat) from a netCDF file with a "lat" axis
        ! [deg N], interpolate linearly onto the model latitudes and close the
        ! mass flux of w* (remove the Gauss-weighted mean per month).
        implicit none
        type(aeros_bd_class),   intent(inout) :: bd
        character(len=*),       intent(in)    :: filename
        type(aeros_grid_class), intent(in)    :: grd
        real(wp), allocatable :: latf(:), wf(:,:), sf(:,:)
        real(wp) :: wsum, wmean
        integer  :: nlatf, nmon, m
        nlatf = nc_size(filename, "lat"); nmon = nc_size(filename, "month")
        if (nmon /= 12) then
            write(io_unit_err,*) "aeros_bd_read:: error: expected 12 months in ", trim(filename)
            error stop 1
        end if
        allocate(latf(nlatf), wf(nlatf, 12), sf(nlatf, 12))
        call nc_read(filename, "lat", latf)
        call nc_read(filename, "wstar", wf)
        call nc_read(filename, "stab", sf)
        if (latf(1) > latf(nlatf)) then          ! ascending latitudes
            latf = latf(nlatf:1:-1); wf = wf(nlatf:1:-1, :); sf = sf(nlatf:1:-1, :)
        end if
        wsum = sum(grd%gauss_w(1:grd%nlat))
        do m = 1, 12
            call interp_lat(latf, wf(:, m), grd%lat(1:grd%nlat), bd%wstar(:, m))
            call interp_lat(latf, sf(:, m), grd%lat(1:grd%nlat), bd%stab(:, m))
            wmean = sum(bd%wstar(:, m)*grd%gauss_w(1:grd%nlat))/wsum
            bd%wstar(:, m) = bd%wstar(:, m) - wmean
        end do
        deallocate(latf, wf, sf)
        return
    contains
        pure subroutine interp_lat(xs, ys, xt, yt)
            real(wp), intent(in)  :: xs(:), ys(:), xt(:)
            real(wp), intent(out) :: yt(:)
            integer :: j, jf, n
            n = size(xs)
            do j = 1, size(xt)
                if (xt(j) <= xs(1)) then
                    yt(j) = ys(1)
                else if (xt(j) >= xs(n)) then
                    yt(j) = ys(n)
                else
                    jf = 1
                    do while (xs(jf+1) < xt(j)); jf = jf + 1; end do
                    yt(j) = ys(jf) + (ys(jf+1) - ys(jf))*(xt(j) - xs(jf))/(xs(jf+1) - xs(jf))
                end if
            end do
        end subroutine interp_lat
    end subroutine aeros_bd_read

    subroutine aeros_bd_end(bd)
        implicit none
        type(aeros_bd_class), intent(inout) :: bd
        if (allocated(bd%wstar)) deallocate(bd%wstar)
        if (allocated(bd%stab))  deallocate(bd%stab)
        if (allocated(bd%q))     deallocate(bd%q)
        bd%enabled = .FALSE.; bd%nlon = 0; bd%nlat = 0
        return
    end subroutine aeros_bd_end

    pure real(wp) function aeros_bd_factor(bd) result(f)
        ! The amplitude the circulation runs at: scale x (1 + per_k (tmean - t_ref)),
        ! floored at zero.
        implicit none
        type(aeros_bd_class), intent(in) :: bd
        f = max(0.0_wp, bd%scale*(1.0_wp + bd%per_k*(bd%tmean - bd%t_ref)))
        return
    end function aeros_bd_factor

    subroutine aeros_bd_apply(bd, vg, grd, t_g, lnps_g, dt_phys, dt, day)
        ! Add the lid's Brewer-Dobson increment -w* S_ref dt to dt_phys(:,:,1)
        ! at day-of-year `day` (mid-month-centred linear interpolation of the
        ! monthly w* and S_ref, cyclic, as the seasonal SST). Updates the
        ! running global-mean lowest-level temperature first.
        implicit none
        type(aeros_bd_class),    intent(inout) :: bd
        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_grid_class),  intent(in)    :: grd
        real(wp), intent(in)    :: t_g(:,:,:)      ! (nlon,nlat,nlev) [K]
        real(wp), intent(in)    :: lnps_g(:,:)     ! (nlon,nlat) ln[Pa]
        real(wp), intent(inout) :: dt_phys(:,:,:)  ! (nlon,nlat,nlev) [K] increment
        real(wp), intent(in)    :: dt, day
        real(wp), parameter :: DAYS_YR = 365.0_wp
        real(wp) :: x, w, tg, ws, sr, a, fac
        integer  :: i, j, im0, m0, m1, nlev, ks
        if (.not. bd%enabled) return
        nlev = vg%nlev; ks = nlev
        ! running annual mean of the global-mean lowest-level temperature
        tg = 0.0_wp
        do j = 1, grd%nlat
            tg = tg + sum(t_g(:, j, ks))/real(grd%nlon, wp)*grd%gauss_w(j)
        end do
        tg = tg/sum(grd%gauss_w(1:grd%nlat))
        a  = dt/bd%tau_mean
        bd%tmean = bd%tmean + a*(tg - bd%tmean)
        fac = aeros_bd_factor(bd)
        ! month interpolation
        x   = modulo(day, DAYS_YR)/DAYS_YR*12.0_wp - 0.5_wp
        im0 = floor(x); w = x - real(im0, wp)
        m0  = modulo(im0, 12) + 1; m1 = modulo(im0 + 1, 12) + 1
        do j = 1, grd%nlat
            ws = fac*((1.0_wp - w)*bd%wstar(j, m0) + w*bd%wstar(j, m1))
            sr = (1.0_wp - w)*bd%stab(j, m0) + w*bd%stab(j, m1)
            do i = 1, grd%nlon
                bd%q(i,j) = -ws*sr
                dt_phys(i,j,1) = dt_phys(i,j,1) + bd%q(i,j)*dt
            end do
        end do
        return
    end subroutine aeros_bd_apply

end module aeros_bd
