module aeros_diagnostics
    ! Zonal-mean and eddy statistics, accumulated online.
    !
    ! The Held-Suarez benchmark is a statement about the CLIMATE of a
    ! dynamical core, not about any of its states: Held & Suarez ask for a
    ! 1200-day integration and compare zonal-mean time-mean fields over the
    ! last 1000 days. Nothing here is meaningful for a single snapshot.
    !
    ! === Why online =========================================================
    !
    ! Because the alternative does not fit on disk. A 1000-day average from
    ! daily T42L20 snapshots is 1000 x 128 x 64 x 20 x 8 bytes per 3-D field,
    ! i.e. 1.3 GB for u alone and ~5 GB for the four fields the eddy statistics
    ! need -- to produce 64 x 20 arrays of output. Accumulating as the model
    ! runs turns that into ~40 kB, and it is also the only form in which the
    ! averaging window is a property of the run rather than of whatever
    ! post-processing script was used.
    !
    ! === The decomposition ==================================================
    !
    ! Primes are departures from the INSTANTANEOUS zonal mean, and the time
    ! average is taken afterwards:
    !
    !     [X](j,k,t) = (1/nlon) sum_i X(i,j,k,t)
    !     X'         = X - [X]
    !     <[u'v']>   = < [uv] - [u][v] >
    !
    ! so the reported eddy fluxes are TOTAL eddy fluxes -- transient plus
    ! stationary -- which is what Held & Suarez plot. Accumulating [uv] and
    ! [u][v] separately and differencing at the end would instead give the
    ! transient part only, and would silently be a different diagnostic.
    !
    ! === Sampling ===========================================================
    !
    ! The caller decides when to sample, by calling aeros_diag_accumulate. It
    ! must be after aeros_timestep_diagnose (or aeros_update, which ends with
    ! one), since the grid-space fields are the input. Sampling every dynamics
    ! step would double the model's transform cost for no statistical benefit:
    ! synoptic eddies evolve over days, and a daily sample over 1000 days is
    ! already 1000 effectively independent-ish samples. Every sample carries
    ! equal weight, so the sampling interval must be constant across the
    ! averaging window.

    use aeros_defs, only : dp, wp, io_unit_err, grav, &
                            aeros_grid_class, aeros_state_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure

    implicit none

    private

    type aeros_diag_class
        ! Running sums over the averaging window, divided through by
        ! aeros_diag_finalize. `nsample` is the only normalization, so a
        ! partially accumulated object is still valid -- just noisier.

        integer  :: nlat = 0, nlev = 0
        integer  :: nsample = 0
        ! Window bounds, in whatever unit the caller passes to
        ! aeros_diag_accumulate. Recorded rather than interpreted: the
        ! benchmark driver works in days and the facade in years, and the
        ! averaging window is only ever compared with itself.
        real(dp) :: time_start = 0.0_dp
        real(dp) :: time_end   = 0.0_dp
        logical  :: finalized  = .FALSE.

        ! Zonal-mean, time-mean fields, (nlat,nlev).
        real(dp), allocatable :: u(:,:)        ! [m s-1]
        real(dp), allocatable :: v(:,:)        ! [m s-1]
        real(dp), allocatable :: temp(:,:)     ! [K]

        ! Eddy statistics, (nlat,nlev).
        real(dp), allocatable :: uv(:,:)       ! <[u'v']> [m2 s-2]
        real(dp), allocatable :: vt(:,:)       ! <[v'T']> [K m s-1]
        real(dp), allocatable :: eke(:,:)      ! <[(u'^2+v'^2)/2]> [m2 s-2]
        real(dp), allocatable :: tvar(:,:)     ! <[T'^2]> [K2]

        ! Single-level, (nlat).
        real(dp), allocatable :: ps(:)         ! [Pa]
        real(dp), allocatable :: u_sfc(:)      ! lowest-level [u] [m s-1]

        ! Pressure of each full level at the time-mean zonal-mean surface
        ! pressure, (nlat,nlev). Not a constant in a sigma coordinate, and the
        ! axis every Held-Suarez figure is plotted against.
        real(dp), allocatable :: pressure(:,:) ! [Pa]
    end type aeros_diag_class

    public :: aeros_diag_class
    public :: aeros_diag_init
    public :: aeros_diag_end
    public :: aeros_diag_accumulate
    public :: aeros_diag_finalize
    public :: aeros_diag_summary

contains

    subroutine aeros_diag_init(dgn, grd, nlev)

        implicit none

        type(aeros_diag_class), intent(inout) :: dgn
        type(aeros_grid_class), intent(in)    :: grd
        integer, intent(in) :: nlev

        call aeros_diag_end(dgn)

        dgn%nlat = grd%nlat
        dgn%nlev = nlev

        allocate(dgn%u(grd%nlat,nlev), dgn%v(grd%nlat,nlev), dgn%temp(grd%nlat,nlev))
        allocate(dgn%uv(grd%nlat,nlev), dgn%vt(grd%nlat,nlev))
        allocate(dgn%eke(grd%nlat,nlev), dgn%tvar(grd%nlat,nlev))
        allocate(dgn%ps(grd%nlat), dgn%u_sfc(grd%nlat))
        allocate(dgn%pressure(grd%nlat,nlev))

        dgn%u = 0.0_dp; dgn%v = 0.0_dp; dgn%temp = 0.0_dp
        dgn%uv = 0.0_dp; dgn%vt = 0.0_dp
        dgn%eke = 0.0_dp; dgn%tvar = 0.0_dp
        dgn%ps = 0.0_dp; dgn%u_sfc = 0.0_dp
        dgn%pressure = 0.0_dp

        dgn%nsample   = 0
        dgn%finalized = .FALSE.

        return

    end subroutine aeros_diag_init

    subroutine aeros_diag_end(dgn)

        implicit none

        type(aeros_diag_class), intent(inout) :: dgn

        dgn%nlat = 0; dgn%nlev = 0; dgn%nsample = 0
        dgn%finalized = .FALSE.

        if (allocated(dgn%u))        deallocate(dgn%u)
        if (allocated(dgn%v))        deallocate(dgn%v)
        if (allocated(dgn%temp))     deallocate(dgn%temp)
        if (allocated(dgn%uv))       deallocate(dgn%uv)
        if (allocated(dgn%vt))       deallocate(dgn%vt)
        if (allocated(dgn%eke))      deallocate(dgn%eke)
        if (allocated(dgn%tvar))     deallocate(dgn%tvar)
        if (allocated(dgn%ps))       deallocate(dgn%ps)
        if (allocated(dgn%u_sfc))    deallocate(dgn%u_sfc)
        if (allocated(dgn%pressure)) deallocate(dgn%pressure)

        return

    end subroutine aeros_diag_end

    subroutine aeros_diag_accumulate(dgn, grd, now, time)
        ! Add one sample. `now` must have current grid-space fields.

        implicit none

        type(aeros_diag_class),  intent(inout) :: dgn
        type(aeros_grid_class),  intent(in)    :: grd
        type(aeros_state_class), intent(in)    :: now
        real(wp), intent(in) :: time    ! caller's units, for the record only

        real(dp) :: rn, um, vm, tm, uvm, vtm, ekem, tvm, up, vp, tp
        integer  :: i, j, k

        if (dgn%finalized) then
            write(io_unit_err,*) "aeros_diag_accumulate:: error: cannot add a sample to a "// &
                                    "finalized average. Call aeros_diag_init to start a new one."
            error stop 1
        end if

        rn = 1.0_dp/real(grd%nlon, dp)

        if (dgn%nsample == 0) dgn%time_start = real(time, dp)
        dgn%time_end = real(time, dp)

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,um,vm,tm,uvm,vtm,ekem,tvm,up,vp,tp)
        do k = 1, dgn%nlev
            do j = 1, dgn%nlat

                ! Instantaneous zonal means.
                um = 0.0_dp; vm = 0.0_dp; tm = 0.0_dp
                do i = 1, grd%nlon
                    um = um + real(now%u(i,j,k), dp)
                    vm = vm + real(now%v(i,j,k), dp)
                    tm = tm + real(now%temp_g(i,j,k), dp)
                end do
                um = um*rn; vm = vm*rn; tm = tm*rn

                ! Instantaneous eddy covariances about them.
                uvm = 0.0_dp; vtm = 0.0_dp; ekem = 0.0_dp; tvm = 0.0_dp
                do i = 1, grd%nlon
                    up = real(now%u(i,j,k), dp)      - um
                    vp = real(now%v(i,j,k), dp)      - vm
                    tp = real(now%temp_g(i,j,k), dp) - tm
                    uvm  = uvm  + up*vp
                    vtm  = vtm  + vp*tp
                    ekem = ekem + 0.5_dp*(up*up + vp*vp)
                    tvm  = tvm  + tp*tp
                end do

                dgn%u(j,k)    = dgn%u(j,k)    + um
                dgn%v(j,k)    = dgn%v(j,k)    + vm
                dgn%temp(j,k) = dgn%temp(j,k) + tm
                dgn%uv(j,k)   = dgn%uv(j,k)   + uvm*rn
                dgn%vt(j,k)   = dgn%vt(j,k)   + vtm*rn
                dgn%eke(j,k)  = dgn%eke(j,k)  + ekem*rn
                dgn%tvar(j,k) = dgn%tvar(j,k) + tvm*rn

            end do
        end do
        !$omp end parallel do

        do j = 1, dgn%nlat
            um = 0.0_dp
            do i = 1, grd%nlon
                um = um + real(now%ps(i,j), dp)
            end do
            dgn%ps(j) = dgn%ps(j) + um*rn

            um = 0.0_dp
            do i = 1, grd%nlon
                um = um + real(now%u(i,j,dgn%nlev), dp)
            end do
            dgn%u_sfc(j) = dgn%u_sfc(j) + um*rn
        end do

        dgn%nsample = dgn%nsample + 1

        return

    end subroutine aeros_diag_accumulate

    subroutine aeros_diag_finalize(dgn, vg)
        ! Divide the sums through, and build the pressure axis.

        implicit none

        type(aeros_diag_class),  intent(inout) :: dgn
        type(aeros_vgrid_class), intent(in)    :: vg

        real(wp) :: p_half(0:vg%nlev), p_full(vg%nlev)
        real(dp) :: rn
        integer  :: j, k

        if (dgn%finalized) return

        if (dgn%nsample < 1) then
            write(io_unit_err,*) "aeros_diag_finalize:: error: no samples were accumulated. "// &
                                    "Check that the averaging window opens before the run ends."
            error stop 1
        end if

        rn = 1.0_dp/real(dgn%nsample, dp)

        dgn%u = dgn%u*rn; dgn%v = dgn%v*rn; dgn%temp = dgn%temp*rn
        dgn%uv = dgn%uv*rn; dgn%vt = dgn%vt*rn
        dgn%eke = dgn%eke*rn; dgn%tvar = dgn%tvar*rn
        dgn%ps = dgn%ps*rn; dgn%u_sfc = dgn%u_sfc*rn

        ! Full-level pressures at the mean zonal-mean surface pressure. Not a
        ! constant in a sigma coordinate, which is why it is stored per
        ! latitude rather than left for the reader to reconstruct.
        do j = 1, dgn%nlat
            call aeros_vgrid_pressure(vg, real(dgn%ps(j), wp), p_half, p_full)
            do k = 1, dgn%nlev
                dgn%pressure(j,k) = real(p_full(k), dp)
            end do
        end do

        dgn%finalized = .TRUE.

        return

    end subroutine aeros_diag_finalize

    subroutine aeros_diag_summary(dgn, grd, io_unit)
        ! The numbers a Held-Suarez result is judged on, printed.
        !
        ! Held & Suarez publish figures rather than numbers, so a comparison is
        ! ultimately visual. These are the quantities one reads off those
        ! figures, reduced to scalars so that a run can be regression-checked
        ! against a previous one without a plotting step.

        implicit none

        type(aeros_diag_class), intent(in) :: dgn
        type(aeros_grid_class), intent(in) :: grd
        integer, intent(in), optional :: io_unit

        real(dp) :: umax, ulat, upres, usfc_max, usfc_lat
        real(dp) :: vt_max, vt_lat, uv_max, eke_max, tsfc_eq, tsfc_pole
        real(dp) :: asym_u, num, den
        integer  :: iou, j, k, jmax, kmax, jm2

        iou = 6
        if (present(io_unit)) iou = io_unit

        ! --- Jet: the maximum of the zonal-mean zonal wind, and where it is.
        umax = -1.0e30_dp; jmax = 1; kmax = 1
        do k = 1, dgn%nlev
            do j = 1, dgn%nlat
                if (dgn%u(j,k) > umax) then
                    umax = dgn%u(j,k); jmax = j; kmax = k
                end if
            end do
        end do
        ulat  = real(grd%lat(jmax), dp)
        upres = dgn%pressure(jmax,kmax)

        ! --- Surface westerlies.
        usfc_max = -1.0e30_dp; jm2 = 1
        do j = 1, dgn%nlat
            if (dgn%u_sfc(j) > usfc_max) then
                usfc_max = dgn%u_sfc(j); jm2 = j
            end if
        end do
        usfc_lat = real(grd%lat(jm2), dp)

        ! --- Eddy fluxes and variances.
        vt_max = -1.0e30_dp; vt_lat = 0.0_dp
        uv_max = -1.0e30_dp; eke_max = -1.0e30_dp
        do k = 1, dgn%nlev
            do j = 1, dgn%nlat
                if (dgn%vt(j,k) > vt_max) then
                    vt_max = dgn%vt(j,k); vt_lat = real(grd%lat(j), dp)
                end if
                uv_max  = max(uv_max,  dgn%uv(j,k))
                eke_max = max(eke_max, dgn%eke(j,k))
            end do
        end do

        ! --- Surface temperature at the equator and the pole. `grd%lat` runs
        ! north to south, so the equator is the middle and the pole is the end.
        tsfc_eq   = dgn%temp(dgn%nlat/2, dgn%nlev)
        tsfc_pole = dgn%temp(1, dgn%nlev)

        ! --- Hemispheric asymmetry of [u]. The forcing is symmetric, so this
        ! measures how far the run is from a converged climate: it falls as the
        ! averaging window lengthens and never reaches zero.
        num = 0.0_dp; den = 0.0_dp
        do k = 1, dgn%nlev
            do j = 1, dgn%nlat
                num = num + (dgn%u(j,k) - dgn%u(dgn%nlat+1-j,k))**2
                den = den + dgn%u(j,k)**2
            end do
        end do
        asym_u = sqrt(num/max(den, tiny(1.0_dp)))

        write(iou,*) ""
        write(iou,"(a)") " == Held-Suarez climate =="
        write(iou,"(a,i0,a,f10.2,a,f10.2,a)") "   samples ", dgn%nsample, &
                        " over [", dgn%time_start, ",", dgn%time_end, "]"
        write(iou,*) ""
        write(iou,"(a,f9.2,a)")  "   max zonal-mean [u]          ", umax, " m s-1"
        write(iou,"(a,f9.2,a)")  "     at latitude               ", ulat, " deg"
        write(iou,"(a,f9.2,a)")  "     at pressure               ", upres/100.0_dp, " hPa"
        write(iou,"(a,f9.2,a)")  "   max surface [u]             ", usfc_max, " m s-1"
        write(iou,"(a,f9.2,a)")  "     at latitude               ", usfc_lat, " deg"
        write(iou,"(a,f9.2,a)")  "   max eddy heat flux [v'T']   ", vt_max, " K m s-1"
        write(iou,"(a,f9.2,a)")  "     at latitude               ", vt_lat, " deg"
        write(iou,"(a,f9.2,a)")  "   max eddy momentum flux [u'v']", uv_max, " m2 s-2"
        write(iou,"(a,f9.2,a)")  "   max eddy kinetic energy     ", eke_max, " m2 s-2"
        write(iou,"(a,f9.2,a)")  "   surface [T] at the equator  ", tsfc_eq, " K"
        write(iou,"(a,f9.2,a)")  "   surface [T] at the pole     ", tsfc_pole, " K"
        write(iou,"(a,es12.3)")  "   hemispheric asymmetry of [u]", asym_u
        write(iou,*) ""

        return

    end subroutine aeros_diag_summary

end module aeros_diagnostics
