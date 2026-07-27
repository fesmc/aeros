program validate_era5
    ! Offline validation of the clear-sky radiation operators against ERA5.
    !
    ! The clear-sky longwave and shortwave column kernels (aeros_radiation) are
    ! grid-agnostic in the vertical: they take an arbitrary nlev column of
    ! T, q, o3, layer thickness and interface height. That is the whole trick
    ! here -- rather than regrid ERA5 onto the model sigma grid (and confound a
    ! transfer error with an interpolation error), we drive the *exact*
    ! operators on ERA5's native 37 pressure levels, one 2.5-degree column at a
    ! time, and compare the emergent TOA and surface clear-sky fluxes against
    ! ERA5's own clear-sky diagnostics on the same grid. No preprocessing, no
    ! regridding: the radiative transfer is the only thing under test.
    !
    ! Data: ERA5 1991-2020 monthly climatology (cdo ymonmean), 2.5 deg
    ! (144 x 73), 37 pressure levels. All fields are averaged over the 12
    ! climatological months to an annual mean before the comparison. ERA5 flux
    ! diagnostics are stored as daily accumulations [J m-2]; divide by 86400 to
    ! get [W m-2]. Net fluxes are downward-positive in ERA5's convention, so
    ! OLR = -ttrc, etc. (see the assembly of each *_era below).
    !
    ! Inputs (arg 1 = ERA5 base dir, default below):
    !   pressure-levels/  t, q, o3, z   (K, kg/kg, kg/kg, m2 s-2)  clim
    !   single-levels/    sp, skt, z, fal                          clim
    !                     ttrc, tsrc, strc, ssrc, strdc, ssrdc, tisr  clim
    !
    ! Output (arg 2, default output/era5_rad_validation.nc): for each of six
    ! clear-sky fluxes, the modelled field, the ERA5 field and their bias
    ! (model - ERA5), all (lon,lat) [W m-2], plus the annual-mean insolation
    ! actually fed to the shortwave. A Julia script maps these.

    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use aeros_defs,      only : dp, wp, grav, cp_d, sigma_sb, R_d, p0, S0, pi
    use aeros_radiation, only : aeros_lw_clearsky_column, aeros_lw_cloudy_column, &
                                aeros_sw_clearsky_column, aeros_sw_cloudy_column, &
                                aeros_insolation_daily
    use aeros_cloud,     only : aeros_cloud_diagnose
    use ncio,            only : nc_create, nc_write_dim, nc_write, nc_read, nc_size

    implicit none

    ! --- CO2 for the 1991-2020 climatology period. Global-mean CO2 rose from
    !     ~355 ppm (1991) to ~412 ppm (2020); ~380 ppm is the period mean.
    real(wp), parameter :: CO2_PPM = 380.0_wp
    real(wp), parameter :: M_CO2 = 44.0095_wp, M_AIR = 28.97_wp
    real(wp), parameter :: ACC   = 86400.0_wp     ! J m-2 (daily accum) -> W m-2

    character(len=512) :: base, fout, dpl, dsl
    integer  :: nlon, nlat, nplev, ntime
    integer  :: i, j

    ! ERA5 fields, annual mean. Pressure-level arrays are (lon,lat,plev),
    ! stored surface->top in the file (1000..1 hPa).
    real(wp), allocatable :: lon(:), lat(:), plev(:)
    real(wp), allocatable :: t3(:,:,:), q3(:,:,:), o33(:,:,:), z3(:,:,:)
    real(wp), allocatable :: cc3(:,:,:), clwc3(:,:,:), ciwc3(:,:,:)  ! cloud
    real(wp), allocatable :: sp2(:,:), skt2(:,:), zsfc2(:,:), fal2(:,:)
    real(wp), allocatable :: ttrc(:,:), tsrc(:,:), strc(:,:), ssrc(:,:)
    real(wp), allocatable :: strdc(:,:), ssrdc(:,:), tisr(:,:)
    real(wp), allocatable :: ttr(:,:), tsr(:,:)     ! all-sky TOA thermal/solar
    real(wp), allocatable :: tcc(:,:)               ! total cloud cover [0-1]
    real(wp), allocatable :: coszen(:)              ! annual-mean, per latitude

    ! Outputs (lon,lat): model, ERA5, bias for six fluxes, plus insolation.
    real(wp), allocatable :: olr_m(:,:),   olr_e(:,:)
    real(wp), allocatable :: lwd_m(:,:),   lwd_e(:,:)     ! surface downward LW
    real(wp), allocatable :: lwn_m(:,:),   lwn_e(:,:)     ! surface net LW
    real(wp), allocatable :: swt_m(:,:),   swt_e(:,:)     ! TOA net SW
    real(wp), allocatable :: swd_m(:,:),   swd_e(:,:)     ! surface downward SW
    real(wp), allocatable :: swn_m(:,:),   swn_e(:,:)     ! surface net SW
    real(wp), allocatable :: insol(:,:)                   ! TOA down SW fed in

    ! All-sky TOA (model, ERA5) and the cloud radiative effect (model, ERA5).
    real(wp), allocatable :: olra_m(:,:),  olra_e(:,:)    ! all-sky OLR
    real(wp), allocatable :: swta_m(:,:),  swta_e(:,:)    ! all-sky TOA net SW
    real(wp), allocatable :: crelw_m(:,:), crelw_e(:,:)   ! TOA LW cloud effect
    real(wp), allocatable :: cresw_m(:,:), cresw_e(:,:)   ! TOA SW cloud effect
    real(wp), allocatable :: crenet_m(:,:),crenet_e(:,:)  ! TOA net cloud effect

    ! Diagnosed-cloud path: the aeros_cloud scheme run on the ERA5 columns, then
    ! the cloudy kernels -- the harness for tuning the diagnosis against ERA5.
    real(wp), allocatable :: dcrelw_m(:,:), dcresw_m(:,:), dcrenet_m(:,:)
    real(wp), allocatable :: dcc_m(:,:)                   ! diagnosed cloud cover

    real(wp) :: q_co2

    ! --- paths ---------------------------------------------------------------
    base = "/Users/alrobi001/data/era5"
    fout = "output/era5_rad_validation.nc"
    if (command_argument_count() >= 1) call get_command_argument(1, base)
    if (command_argument_count() >= 2) call get_command_argument(2, fout)
    dpl = trim(base)//"/monthly-pressure-levels/"
    dsl = trim(base)//"/monthly-single-levels/"

    q_co2 = CO2_PPM*1.0e-6_wp * M_CO2/M_AIR

    ! --- dimensions ----------------------------------------------------------
    nlon  = nc_size(trim(pl("t")), "longitude")
    nlat  = nc_size(trim(pl("t")), "latitude")
    nplev = nc_size(trim(pl("t")), "pressure_level")
    ntime = nc_size(trim(pl("t")), "valid_time")
    write(*,"(a,4(i0,1x))") " validate_era5:: nlon,nlat,nplev,ntime = ", &
                            nlon, nlat, nplev, ntime

    allocate(lon(nlon), lat(nlat), plev(nplev), coszen(nlat))
    call nc_read(trim(pl("t")), "longitude", lon)
    call nc_read(trim(pl("t")), "latitude",  lat)
    call nc_read(trim(pl("t")), "pressure_level", plev)   ! hPa

    ! --- read and time-average all fields ------------------------------------
    call read_pl("t",  t3)
    call read_pl("q",  q3)
    call read_pl("o3", o33)
    call read_pl("z",  z3)
    call read_pl("cc",   cc3)      ! cloud fraction [0-1]
    call read_pl("clwc", clwc3)    ! cloud liquid water [kg kg-1]
    call read_pl("ciwc", ciwc3)    ! cloud ice water    [kg kg-1]
    call read_sl("sp",  sp2)
    call read_sl("skt", skt2)
    call read_sl("z",   zsfc2)         ! surface geopotential
    call read_sl("fal", fal2)
    call read_sl("ttrc",  ttrc)
    call read_sl("tsrc",  tsrc)
    call read_sl("strc",  strc)
    call read_sl("ssrc",  ssrc)
    call read_sl("strdc", strdc)
    call read_sl("ssrdc", ssrdc)
    call read_sl("tisr",  tisr)
    call read_sl("ttr",   ttr)         ! all-sky TOA net thermal
    call read_sl("tsr",   tsr)         ! all-sky TOA net solar
    call read_sl("tcc",   tcc)         ! total cloud cover

    ! --- annual-mean insolation-weighted cosine zenith per latitude ----------
    ! Same construction aeros_radiation_init uses; the shortwave needs an
    ! airmass. TOA down SW itself is taken from ERA5 tisr so the insolation
    ! scheme is not a confounder in the transfer comparison.
    call annual_coszen()

    ! --- allocate outputs ----------------------------------------------------
    allocate(olr_m(nlon,nlat), olr_e(nlon,nlat))
    allocate(lwd_m(nlon,nlat), lwd_e(nlon,nlat))
    allocate(lwn_m(nlon,nlat), lwn_e(nlon,nlat))
    allocate(swt_m(nlon,nlat), swt_e(nlon,nlat))
    allocate(swd_m(nlon,nlat), swd_e(nlon,nlat))
    allocate(swn_m(nlon,nlat), swn_e(nlon,nlat))
    allocate(insol(nlon,nlat))
    allocate(olra_m(nlon,nlat), olra_e(nlon,nlat))
    allocate(swta_m(nlon,nlat), swta_e(nlon,nlat))
    allocate(crelw_m(nlon,nlat), crelw_e(nlon,nlat))
    allocate(cresw_m(nlon,nlat), cresw_e(nlon,nlat))
    allocate(crenet_m(nlon,nlat),crenet_e(nlon,nlat))
    allocate(dcrelw_m(nlon,nlat), dcresw_m(nlon,nlat), dcrenet_m(nlon,nlat))
    allocate(dcc_m(nlon,nlat))

    ! --- ERA5 target fields in the operator's conventions [W m-2] ------------
    olr_e = -ttrc/ACC                    ! OLR (up) = -(top net thermal)
    lwd_e =  strdc/ACC                   ! surface downward LW
    lwn_e =  strc/ACC                    ! surface net LW (down - up)
    swt_e =  tsrc/ACC                    ! TOA net SW (down - up)
    swd_e =  ssrdc/ACC                   ! surface downward SW
    swn_e =  ssrc/ACC                    ! surface net SW (absorbed)

    olra_e = -ttr/ACC                    ! all-sky OLR
    swta_e =  tsr/ACC                    ! all-sky TOA net SW
    ! Cloud radiative effect = all-sky minus clear-sky net-downward flux
    ! (positive = clouds warm), the §17 convention. At TOA net-down LW = -OLR.
    crelw_e  = olr_e  - olra_e           ! = (ttr - ttrc)/ACC  > 0
    cresw_e  = swta_e - swt_e            ! = (tsr - tsrc)/ACC  < 0
    crenet_e = crelw_e + cresw_e

    ! --- column-by-column transfer -------------------------------------------
    !$omp parallel do collapse(2) schedule(dynamic) private(i,j)
    do j = 1, nlat
        do i = 1, nlon
            call one_column(i, j)
        end do
    end do
    !$omp end parallel do

    call report_global_means()
    call write_output()

    write(*,"(a)") " validate_era5:: wrote "//trim(fout)

contains

    ! ------------------------------------------------------------------------
    ! Column driver: build the ERA5 column top->surface, run both operators,
    ! store modelled fluxes.
    ! ------------------------------------------------------------------------
    subroutine one_column(i, j)
        integer, intent(in) :: i, j
        integer  :: n, k, kk
        real(wp) :: ps, ts, alb, swin
        real(wp) :: pcol(nplev), tcol(nplev), qcol(nplev), ocol(nplev), hcol(nplev)
        real(wp) :: ccol(nplev), lcol(nplev), icol(nplev)
        real(wp) :: cfd(nplev), lcd(nplev), icd(nplev)   ! diagnosed cloud
        real(wp) :: dpc(nplev), zhalf(0:nplev), phalf(0:nplev)
        real(wp) :: fnet(0:nplev), heat(nplev)
        real(wp) :: olr, fdw_lw, sw_up, sw_dw, sw_net
        real(wp) :: olr_a, fdw_a, swup_a, swdw_a, swnet_a
        real(wp) :: olr_d, fdw_d, swup_d, swdw_d, swnet_d, ccmax

        ps  = sp2(i,j)                    ! Pa
        ts  = skt2(i,j)
        alb = min(max(fal2(i,j), 0.0_wp), 0.9_wp)

        ! Collect valid pressure levels (above ground, not fill), reversed to
        ! top(k=1) -> surface(k=n). File order is surface->top (1000..1 hPa).
        n = 0
        do kk = nplev, 1, -1
            if (plev(kk)*100.0_wp >= ps) cycle          ! at/below ground
            if (ieee_is_nan(t3(i,j,kk))) cycle
            n = n + 1
            pcol(n) = plev(kk)*100.0_wp                 ! Pa
            tcol(n) = t3(i,j,kk)
            qcol(n) = max(0.0_wp, q3(i,j,kk))
            ocol(n) = max(0.0_wp, o33(i,j,kk))
            hcol(n) = z3(i,j,kk)/grav                   ! geopotential height [m]
            ccol(n) = min(1.0_wp, max(0.0_wp, cc3(i,j,kk)))
            lcol(n) = max(0.0_wp, clwc3(i,j,kk))
            icol(n) = max(0.0_wp, ciwc3(i,j,kk))
        end do
        if (n < 3) then                                 ! degenerate column
            call set_missing(i, j); return
        end if

        ! Interfaces: top at half the top level's pressure, interior at level
        ! midpoints, base at the surface pressure. Layer k = [phalf(k-1),phalf(k)].
        phalf(0) = 0.5_wp*pcol(1)
        do k = 1, n-1
            phalf(k) = 0.5_wp*(pcol(k) + pcol(k+1))
        end do
        phalf(n) = ps
        do k = 1, n
            dpc(k) = phalf(k) - phalf(k-1)
        end do

        ! Interface heights: interior midpoints, surface from orography, top by
        ! upward extrapolation. Used only for the CO2/O3 pressure-broadening
        ! path, so robustness matters more than precision.
        zhalf(n) = zsfc2(i,j)/grav
        do k = 1, n-1
            zhalf(k) = 0.5_wp*(hcol(k) + hcol(k+1))
        end do
        zhalf(0) = hcol(1) + (hcol(1) - zhalf(1))
        ! guard monotonicity (heights decrease with index toward the surface)
        do k = 1, n
            if (zhalf(k-1) <= zhalf(k)) zhalf(k-1) = zhalf(k) + 10.0_wp
        end do

        call aeros_lw_clearsky_column(n, tcol(1:n), qcol(1:n), ocol(1:n), &
            dpc(1:n), zhalf(0:n), ts, q_co2, .TRUE., &
            fnet(0:n), heat(1:n), olr, fdw_lw)

        swin = tisr(i,j)/ACC
        call aeros_sw_clearsky_column(n, qcol(1:n), ocol(1:n), .TRUE., &
            dpc(1:n), swin, coszen(j), alb, alb, &
            heat(1:n), sw_up, sw_dw, sw_net)

        ! all-sky: the same operators with ERA5's cloud column
        call aeros_lw_cloudy_column(n, tcol(1:n), qcol(1:n), ocol(1:n), &
            dpc(1:n), zhalf(0:n), ts, q_co2, .TRUE., &
            ccol(1:n), lcol(1:n), icol(1:n), &
            fnet(0:n), heat(1:n), olr_a, fdw_a)
        call aeros_sw_cloudy_column(n, qcol(1:n), ocol(1:n), .TRUE., &
            dpc(1:n), swin, coszen(j), alb, alb, &
            ccol(1:n), lcol(1:n), icol(1:n), &
            heat(1:n), swup_a, swdw_a, swnet_a)

        olr_m(i,j) = olr
        lwd_m(i,j) = fdw_lw
        lwn_m(i,j) = fdw_lw - sigma_sb*ts**4          ! net = down - up
        swt_m(i,j) = swin - sw_up                      ! TOA net SW
        swd_m(i,j) = sw_dw
        swn_m(i,j) = sw_net
        insol(i,j) = swin

        olra_m(i,j) = olr_a                            ! all-sky OLR
        swta_m(i,j) = swin - swup_a                    ! all-sky TOA net SW
        crelw_m(i,j)  = olr - olr_a                    ! clear - all-sky OLR
        cresw_m(i,j)  = (swin - swup_a) - (swin - sw_up)
        crenet_m(i,j) = crelw_m(i,j) + cresw_m(i,j)

        ! diagnosed-cloud path: the aeros_cloud scheme on this ERA5 column, then
        ! the same cloudy operators -- the tuning target is that this reproduce
        ! the observed-cloud CRE above and ERA5's total cloud cover.
        call aeros_cloud_diagnose(n, tcol(1:n), qcol(1:n), pcol(1:n), ps, &
            cfd(1:n), lcd(1:n), icd(1:n))
        call aeros_lw_cloudy_column(n, tcol(1:n), qcol(1:n), ocol(1:n), &
            dpc(1:n), zhalf(0:n), ts, q_co2, .TRUE., &
            cfd(1:n), lcd(1:n), icd(1:n), fnet(0:n), heat(1:n), olr_d, fdw_d)
        call aeros_sw_cloudy_column(n, qcol(1:n), ocol(1:n), .TRUE., &
            dpc(1:n), swin, coszen(j), alb, alb, &
            cfd(1:n), lcd(1:n), icd(1:n), heat(1:n), swup_d, swdw_d, swnet_d)

        ccmax = 0.0_wp
        do kk = 1, n
            ccmax = max(ccmax, cfd(kk))                ! max overlap, as the kernels use
        end do
        dcrelw_m(i,j)  = olr - olr_d
        dcresw_m(i,j)  = (swin - swup_d) - (swin - sw_up)
        dcrenet_m(i,j) = dcrelw_m(i,j) + dcresw_m(i,j)
        dcc_m(i,j)     = ccmax
        return
    end subroutine one_column

    subroutine set_missing(i, j)
        integer, intent(in) :: i, j
        real(wp), parameter :: mv = -9999.0_wp
        olr_m(i,j) = mv; lwd_m(i,j) = mv; lwn_m(i,j) = mv
        swt_m(i,j) = mv; swd_m(i,j) = mv; swn_m(i,j) = mv; insol(i,j) = mv
        olra_m(i,j) = mv; swta_m(i,j) = mv
        crelw_m(i,j) = mv; cresw_m(i,j) = mv; crenet_m(i,j) = mv
        dcrelw_m(i,j) = mv; dcresw_m(i,j) = mv; dcrenet_m(i,j) = mv; dcc_m(i,j) = mv
        return
    end subroutine set_missing

    ! ------------------------------------------------------------------------
    ! I/O helpers
    ! ------------------------------------------------------------------------
    function pl(v) result(f)
        character(len=*), intent(in) :: v
        character(len=512) :: f
        f = trim(dpl)//"era5_monthly-pressure-levels_"//trim(v)//"_1991-2020_clim.nc"
        return
    end function pl

    function sl(v) result(f)
        character(len=*), intent(in) :: v
        character(len=512) :: f
        f = trim(dsl)//"era5_monthly-single-levels_"//trim(v)//"_1991-2020_clim.nc"
        return
    end function sl

    subroutine read_pl(v, out)
        character(len=*), intent(in) :: v
        real(wp), allocatable, intent(out) :: out(:,:,:)
        real(wp), allocatable :: tmp(:,:,:,:)
        integer :: it
        allocate(tmp(nlon,nlat,nplev,ntime))
        call nc_read(trim(pl(v)), v, tmp)
        allocate(out(nlon,nlat,nplev))
        out = tmp(:,:,:,1)
        do it = 2, ntime
            out = out + tmp(:,:,:,it)
        end do
        out = out/real(ntime, wp)
        deallocate(tmp)
        return
    end subroutine read_pl

    subroutine read_sl(v, out)
        character(len=*), intent(in) :: v
        real(wp), allocatable, intent(out) :: out(:,:)
        real(wp), allocatable :: tmp(:,:,:)
        integer :: it
        allocate(tmp(nlon,nlat,ntime))
        call nc_read(trim(sl(v)), v, tmp)
        allocate(out(nlon,nlat))
        out = tmp(:,:,1)
        do it = 2, ntime
            out = out + tmp(:,:,it)
        end do
        out = out/real(ntime, wp)
        deallocate(tmp)
        return
    end subroutine read_sl

    subroutine annual_coszen()
        integer  :: jj, d
        real(wp) :: cz, sw, czsum, wsum
        do jj = 1, nlat
            czsum = 0.0_wp; wsum = 0.0_wp
            do d = 1, 360
                call aeros_insolation_daily(lat(jj)*pi/180.0_wp, real(d,wp), &
                                            real(S0,wp), cz, sw)
                czsum = czsum + cz*sw
                wsum  = wsum + sw
            end do
            if (wsum > 0.0_wp) then
                coszen(jj) = czsum/wsum
            else
                coszen(jj) = 0.0_wp
            end if
        end do
        return
    end subroutine annual_coszen

    ! Area-weighted global mean over valid cells only.
    real(wp) function gmean(a) result(m)
        real(wp), intent(in) :: a(:,:)
        real(wp) :: num, den, w
        integer  :: ii, jj
        num = 0.0_wp; den = 0.0_wp
        do jj = 1, nlat
            w = cos(lat(jj)*pi/180.0_wp)
            do ii = 1, nlon
                if (a(ii,jj) <= -9990.0_wp) cycle
                num = num + a(ii,jj)*w
                den = den + w
            end do
        end do
        m = num/max(den, tiny(1.0_wp))
        return
    end function gmean

    subroutine report_global_means()
        write(*,"(a)") ""
        write(*,"(a)") "  clear-sky global means [W m-2]      model     ERA5     bias"
        call line("OLR (TOA up)        ", olr_m, olr_e)
        call line("surface down LW     ", lwd_m, lwd_e)
        call line("surface net LW      ", lwn_m, lwn_e)
        call line("TOA net SW          ", swt_m, swt_e)
        call line("surface down SW     ", swd_m, swd_e)
        call line("surface net SW      ", swn_m, swn_e)
        write(*,"(a)") ""
        write(*,"(a)") "  all-sky / cloud effect [W m-2]      model     ERA5     bias"
        call line("all-sky OLR         ", olra_m, olra_e)
        call line("all-sky TOA net SW  ", swta_m, swta_e)
        call line("TOA LW cloud effect ", crelw_m,  crelw_e)
        call line("TOA SW cloud effect ", cresw_m,  cresw_e)
        call line("TOA net cloud effect", crenet_m, crenet_e)
        write(*,"(a)") ""
        write(*,"(a)") "  diagnosed cloud (aeros_cloud on ERA5)  model     ERA5     bias"
        call line("TOA LW cloud effect ", dcrelw_m,  crelw_e)
        call line("TOA SW cloud effect ", dcresw_m,  cresw_e)
        call line("TOA net cloud effect", dcrenet_m, crenet_e)
        call line("total cloud cover   ", dcc_m,     tcc)
        return
    end subroutine report_global_means

    subroutine line(lbl, m, e)
        character(len=*), intent(in) :: lbl
        real(wp), intent(in) :: m(:,:), e(:,:)
        real(wp) :: gm, ge
        gm = gmean(m); ge = gmean(e)
        write(*,"(a,3f9.2)") "    "//lbl, gm, ge, gm-ge
        return
    end subroutine line

    subroutine write_output()
        call nc_create(trim(fout))
        call nc_write_dim(trim(fout), "lon", x=lon, units="degrees_east")
        call nc_write_dim(trim(fout), "lat", x=lat, units="degrees_north")
        call wmap("olr_mod",  olr_m); call wmap("olr_era",  olr_e)
        call wbias("olr", olr_m, olr_e)
        call wmap("lwdn_sfc_mod", lwd_m); call wmap("lwdn_sfc_era", lwd_e)
        call wbias("lwdn_sfc", lwd_m, lwd_e)
        call wmap("lwnet_sfc_mod", lwn_m); call wmap("lwnet_sfc_era", lwn_e)
        call wbias("lwnet_sfc", lwn_m, lwn_e)
        call wmap("swnet_toa_mod", swt_m); call wmap("swnet_toa_era", swt_e)
        call wbias("swnet_toa", swt_m, swt_e)
        call wmap("swdn_sfc_mod", swd_m); call wmap("swdn_sfc_era", swd_e)
        call wbias("swdn_sfc", swd_m, swd_e)
        call wmap("swnet_sfc_mod", swn_m); call wmap("swnet_sfc_era", swn_e)
        call wbias("swnet_sfc", swn_m, swn_e)
        call wmap("insol", insol)
        ! all-sky and cloud radiative effect
        call wmap("olr_allsky_mod", olra_m); call wmap("olr_allsky_era", olra_e)
        call wbias("olr_allsky", olra_m, olra_e)
        call wmap("swnet_toa_allsky_mod", swta_m); call wmap("swnet_toa_allsky_era", swta_e)
        call wbias("swnet_toa_allsky", swta_m, swta_e)
        call wmap("cre_lw_toa_mod", crelw_m); call wmap("cre_lw_toa_era", crelw_e)
        call wbias("cre_lw_toa", crelw_m, crelw_e)
        call wmap("cre_sw_toa_mod", cresw_m); call wmap("cre_sw_toa_era", cresw_e)
        call wbias("cre_sw_toa", cresw_m, cresw_e)
        call wmap("cre_net_toa_mod", crenet_m); call wmap("cre_net_toa_era", crenet_e)
        call wbias("cre_net_toa", crenet_m, crenet_e)
        ! diagnosed-cloud path (aeros_cloud on ERA5 columns) vs ERA5
        call wmap("dcre_lw_toa_mod", dcrelw_m); call wbias("dcre_lw_toa", dcrelw_m, crelw_e)
        call wmap("dcre_sw_toa_mod", dcresw_m); call wbias("dcre_sw_toa", dcresw_m, cresw_e)
        call wmap("dcre_net_toa_mod", dcrenet_m); call wbias("dcre_net_toa", dcrenet_m, crenet_e)
        call nc_write(trim(fout), "tcc_mod", dcc_m, dim1="lon", dim2="lat", &
                      units="1", missing_value=-9999.0_wp)
        call nc_write(trim(fout), "tcc_era", tcc, dim1="lon", dim2="lat", units="1")
        return
    end subroutine write_output

    subroutine wmap(name, a)
        character(len=*), intent(in) :: name
        real(wp), intent(in) :: a(:,:)
        call nc_write(trim(fout), name, a, dim1="lon", dim2="lat", &
                      units="W m-2", missing_value=-9999.0_wp)
        return
    end subroutine wmap

    subroutine wbias(name, m, e)
        character(len=*), intent(in) :: name
        real(wp), intent(in) :: m(:,:), e(:,:)
        real(wp) :: b(nlon,nlat)
        integer  :: ii, jj
        do jj = 1, nlat
            do ii = 1, nlon
                if (m(ii,jj) <= -9990.0_wp) then
                    b(ii,jj) = -9999.0_wp
                else
                    b(ii,jj) = m(ii,jj) - e(ii,jj)
                end if
            end do
        end do
        call nc_write(trim(fout), trim(name)//"_bias", b, dim1="lon", dim2="lat", &
                      units="W m-2", missing_value=-9999.0_wp)
        return
    end subroutine wbias

end program validate_era5
