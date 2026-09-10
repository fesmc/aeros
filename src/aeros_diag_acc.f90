module aeros_diag_acc
    ! Library-side diagnostic accumulators (Wave 2 of the coupling API,
    ! docs/coupling_api.md sec. 5, 6.2, 6.3, 12).
    !
    ! Owns the ~30 accumulator buffers that today's rce_long driver keeps at
    ! program scope: the equilibrated second-half time-mean fields (omega, the
    ! per-term diabatic heating split, the jet / eddy momentum flux, the eddy
    ! zonal-wavenumber spectrum), the twelve-month climatology bins (T/u/v/RH
    ! zonal-mean, 2D surface climatology, radiation / cloud, eddy fluxes,
    ! tracer-mass-rescaling covariances, mean 3D fields for stationary eddy
    ! diagnostics), and the running water and energy budget scalars.
    !
    ! Wave 2 (this module) only moves the STORAGE library-side; the accumulate
    ! arithmetic still lives at the rce_long call sites, referencing the fields
    ! through the type. Wave 3 migrates the accumulate/finalize methods here.
    ! The layout, allocation timing and zeroing preserve today's rce_long
    ! semantics exactly, so the existing regression climatologies stay
    ! bit-for-bit.
    !
    ! Init contract:
    !   aeros_diag_acc_init(acc, nlon, nlat, nlev, mmax_spec, lon, dt)
    !       allocates the always-on tavg buffers (nlon,nlat,nlev), the
    !       (nlat,nlev) zonal jet/eddy accumulators, the (mmax_spec) spectra
    !       and their (mmax_spec,nlon) DFT trig tables, zeroes every field
    !       (arrays and scalars), and sets `nacc_mon` from dt.
    !   aeros_diag_acc_init_monthly(acc, nlon, nlat, nlev)
    !       allocates the twelve-month climatology bins and zeroes them, plus
    !       the mon_cnt(NMON) sample counter. Called only when the driver
    !       requests a monthly dump (rce_long: monthly_out /= "").
    !   aeros_diag_acc_end(acc)
    !       deallocates all buffers, restores the type to its uninitialized
    !       state (allocation status = FALSE, scalars = 0).
    !
    ! The type carries every scalar counter (n_wacc, mon_cnt, nacc_mon,
    ! wb_*, tb_*, mmax_spec) that the accumulators need to be interpreted --
    ! so ownership of the accumulator lifecycle is a single struct.

    use aeros_defs,         only : wp, dp, R_d, grav, &
                                    aeros_param_class, aeros_grid_class, &
                                    aeros_state_class
    use aeros_timestep,     only : aeros_timestep_class
    use aeros_vertical,     only : aeros_vgrid_class, aeros_vgrid_pressure
    use aeros_physics,      only : aeros_physics_class
    use aeros_condensation, only : aeros_qsat
    use aeros_cloud,        only : aeros_cloud_diagnose

    implicit none

    private

    integer, parameter, public :: AEROS_DIAG_NMON = 12

    type, public :: aeros_diag_acc_t

        ! -- second-half time-mean, per-column (nlon,nlat,nlev) ---------------
        real(wp), allocatable :: w_tavg(:,:,:)      ! omega
        real(wp), allocatable :: cnv_tavg(:,:,:)    ! forward-split convection heating
        real(wp), allocatable :: cnd_tavg(:,:,:)    ! forward-split condensation heating
        real(wp), allocatable :: rad_tavg(:,:,:)    ! forward-split radiation heating
        real(wp), allocatable :: surf_tavg(:,:,:)   ! forward-split surface heating
        real(wp), allocatable :: vadv_tavg(:,:,:)   ! vertical advection + adiabatic
        real(wp), allocatable :: dyn_tavg(:,:,:)    ! total dynamical rate
        real(wp), allocatable :: vdf_tavg(:,:,:)    ! vertical diffusion incl. surface flux
        real(wp), allocatable :: bd_tavg(:,:,:)     ! prescribed Brewer-Dobson lid heating
        integer               :: n_wacc = 0         ! samples in the second-half tavg

        ! -- second-half zonal-mean jet / eddy momentum (nlat,nlev) -----------
        real(wp), allocatable :: uv_tavg(:,:)       ! [u*v*] zonal-mean eddy momentum flux
        real(wp), allocatable :: ubar_tavg(:,:)     ! zonal-mean [u]
        real(wp), allocatable :: vbar_tavg(:,:)     ! zonal-mean [v]

        ! -- eddy zonal-wavenumber spectrum -----------------------------------
        ! Mass + area weighted, second-half time-mean.
        integer               :: mmax_spec = 0
        real(wp), allocatable :: ke_spec(:)         ! (mmax_spec)  eddy KE per wavenumber
        real(wp), allocatable :: mf_spec(:)         ! (mmax_spec)  [u*v*] co-spectrum
        real(wp), allocatable :: ctab(:,:)          ! (mmax_spec, nlon)  cos DFT table
        real(wp), allocatable :: stab(:,:)          ! (mmax_spec, nlon)  sin DFT table
        real(wp), allocatable :: ke_latm(:,:)       ! (nlat, mmax_spec)  where each m's KE lives
        real(wp), allocatable :: ke_levm(:,:)       ! (nlev, mmax_spec)  same, by level

        ! -- monthly climatology bins (feat/seasonal) --------------------------
        ! (nlat,nlev,NMON) zonal-mean fields
        real(wp), allocatable :: mon_t(:,:,:), mon_u(:,:,:), mon_v(:,:,:), mon_rh(:,:,:)
        real(wp), allocatable :: mon_cf(:,:,:)
        real(wp), allocatable :: mon_vt(:,:,:), mon_vq(:,:,:), mon_eke(:,:,:), mon_uv(:,:,:)
        real(wp), allocatable :: mon_mr(:,:,:), mon_mrt(:,:,:), mon_mrq(:,:,:)
        ! (nlat,NMON) zonal-mean surface fields
        real(wp), allocatable :: mon_sst(:,:), mon_pr(:,:)
        ! (nlon,nlat,NMON) 2D surface / near-surface fields
        real(wp), allocatable :: mon_t2m(:,:,:), mon_tlow(:,:,:)
        real(wp), allocatable :: mon_pr2d(:,:,:), mon_z500(:,:,:)
        real(wp), allocatable :: mon_vlow(:,:,:), mon_evap(:,:,:), mon_shf(:,:,:)
        real(wp), allocatable :: mon_skin(:,:,:), mon_lwd(:,:,:), mon_olr(:,:,:), &
                                 mon_olrcs(:,:,:), mon_cover(:,:,:)
        real(wp), allocatable :: mon_swn(:,:,:), mon_swup(:,:,:)
        real(wp), allocatable :: mon_mr2(:,:,:)
        ! (nlon,nlat,nlev,NMON) full 3D means for stationary-eddy diagnostics
        real(wp), allocatable :: mon_u3(:,:,:,:), mon_v3(:,:,:,:), &
                                 mon_t3(:,:,:,:), mon_q3(:,:,:,:)
        ! samples per calendar-month bin
        integer,  allocatable :: mon_cnt(:)
        ! monthly-sampling cadence [dynamics steps]; ~1 sample/day set in _init
        integer               :: nacc_mon = 1

        ! Spin-up gates (Wave 4d): the driver's WHEN policy, applied inside
        ! the library method so aeros_lib_step can invoke it unconditionally.
        !   spinup_second_half_nstep -- skip _second_half until ts%nstep is
        !     strictly greater than this. rce_long: n0 + nstep/2. Coupled
        !     host: 0 (accumulate from step 1).
        !   spinup_monthly_days -- skip _monthly until ts%nstep*par%dt/86400
        !     >= spinup_monthly_days. rce_long: &rce mon_spinup_days. Coupled
        !     host: 0 (bin from day 0).
        integer  :: spinup_second_half_nstep = 0
        real(wp) :: spinup_monthly_days      = 0.0_wp

        ! -- water budget (report-cadence global means) ------------------------
        real(wp) :: wb_w_last = -1.0_wp
        real(wp) :: wb_p_acc  = 0.0_wp
        real(wp) :: wb_e_acc  = 0.0_wp
        integer  :: wb_n_acc  = 0

        ! -- energy budget (report-cadence global means) -----------------------
        real(wp) :: tb_swin  = 0.0_wp
        real(wp) :: tb_swup  = 0.0_wp
        real(wp) :: tb_olr   = 0.0_wp
        real(wp) :: tb_swupc = 0.0_wp
        real(wp) :: tb_olrc  = 0.0_wp
        real(wp) :: tb_cov   = 0.0_wp
        real(wp) :: tb_shf   = 0.0_wp
        real(wp) :: tb_lhf   = 0.0_wp
        integer  :: tb_n     = 0

    end type aeros_diag_acc_t

    public :: aeros_diag_acc_init
    public :: aeros_diag_acc_init_monthly
    public :: aeros_diag_acc_reset
    public :: aeros_diag_acc_end

    ! Wave 4c: the library owns the accumulator arithmetic. Each of these
    ! methods is a per-step / per-daily-sample increment on `acc`, called by
    ! the driver (which owns the WHEN: spin-up threshold, month bin from the
    ! calendar, report cadence) or by aeros_lib_step on the coupled path.
    public :: aeros_diag_acc_energy         ! per-step TOA/surface energy sums
    public :: aeros_diag_acc_water          ! per-step precip/evap sums
    public :: aeros_diag_acc_second_half    ! per-step 2nd-half tavg + eddy diag
    public :: aeros_diag_acc_monthly        ! per-daily-sample monthly bin
    public :: aeros_diag_acc_imon           ! calendar-day -> month bin index

contains

    subroutine aeros_diag_acc_init(acc, nlon, nlat, nlev, mmax_spec, dt)
        ! Allocate the always-on second-half tavg buffers and the eddy
        ! wavenumber spectrum, and initialise every scalar counter to zero.
        ! Mirrors the block starting at drivers/rce_long.f90:717 line for line,
        ! including the DFT trig table build.

        implicit none

        type(aeros_diag_acc_t), intent(inout) :: acc
        integer,  intent(in)                  :: nlon, nlat, nlev
        integer,  intent(in)                  :: mmax_spec
        real(wp), intent(in)                  :: dt

        integer  :: mm, ii
        real(wp) :: ang, twopi

        allocate(acc%w_tavg(nlon,nlat,nlev));                             acc%w_tavg    = 0.0_wp
        allocate(acc%cnv_tavg(nlon,nlat,nlev));                           acc%cnv_tavg  = 0.0_wp
        allocate(acc%cnd_tavg(nlon,nlat,nlev));                           acc%cnd_tavg  = 0.0_wp
        allocate(acc%rad_tavg(nlon,nlat,nlev));                           acc%rad_tavg  = 0.0_wp
        allocate(acc%surf_tavg(nlon,nlat,nlev));                          acc%surf_tavg = 0.0_wp
        allocate(acc%vadv_tavg(nlon,nlat,nlev));                          acc%vadv_tavg = 0.0_wp
        allocate(acc%dyn_tavg(nlon,nlat,nlev));                           acc%dyn_tavg  = 0.0_wp
        allocate(acc%vdf_tavg(nlon,nlat,nlev));                           acc%vdf_tavg  = 0.0_wp
        allocate(acc%bd_tavg(nlon,nlat,nlev));                            acc%bd_tavg   = 0.0_wp
        acc%n_wacc = 0

        allocate(acc%uv_tavg(nlat,nlev));                                 acc%uv_tavg   = 0.0_wp
        allocate(acc%ubar_tavg(nlat,nlev));                               acc%ubar_tavg = 0.0_wp
        allocate(acc%vbar_tavg(nlat,nlev));                               acc%vbar_tavg = 0.0_wp

        ! Eddy wavenumber spectrum. rce_long today sets mmax_spec = trunc; the
        ! table is pre-built once so the per-step accumulation is multiply-add.
        acc%mmax_spec = mmax_spec
        allocate(acc%ke_spec(mmax_spec));                                 acc%ke_spec = 0.0_wp
        allocate(acc%mf_spec(mmax_spec));                                 acc%mf_spec = 0.0_wp
        allocate(acc%ctab(mmax_spec, nlon));                              acc%ctab    = 0.0_wp
        allocate(acc%stab(mmax_spec, nlon));                              acc%stab    = 0.0_wp
        allocate(acc%ke_latm(nlat, mmax_spec));                           acc%ke_latm = 0.0_wp
        allocate(acc%ke_levm(nlev, mmax_spec));                           acc%ke_levm = 0.0_wp
        twopi = 8.0_wp*atan(1.0_wp)
        do ii = 1, nlon
            do mm = 1, mmax_spec
                ang = twopi*real(mm*(ii-1), wp)/real(nlon, wp)
                acc%ctab(mm,ii) = cos(ang)
                acc%stab(mm,ii) = -sin(ang)
            end do
        end do

        ! Cadence: ~one monthly-climatology sample per day. Also the accumulate
        ! interval for the (much cheaper) monthly-bin fields.
        acc%nacc_mon = max(1, nint(86400.0_wp/dt))

        acc%wb_w_last = -1.0_wp
        acc%wb_p_acc  = 0.0_wp
        acc%wb_e_acc  = 0.0_wp
        acc%wb_n_acc  = 0

        acc%tb_swin  = 0.0_wp
        acc%tb_swup  = 0.0_wp
        acc%tb_olr   = 0.0_wp
        acc%tb_swupc = 0.0_wp
        acc%tb_olrc  = 0.0_wp
        acc%tb_cov   = 0.0_wp
        acc%tb_shf   = 0.0_wp
        acc%tb_lhf   = 0.0_wp
        acc%tb_n     = 0

        return

    end subroutine aeros_diag_acc_init

    subroutine aeros_diag_acc_init_monthly(acc, nlon, nlat, nlev)
        ! Allocate the twelve-month climatology bins and zero every field.
        ! Called only when the driver requests a monthly dump (rce_long:
        ! monthly_out /= "").

        implicit none

        type(aeros_diag_acc_t), intent(inout) :: acc
        integer, intent(in)                   :: nlon, nlat, nlev
        integer, parameter :: NM = AEROS_DIAG_NMON

        allocate(acc%mon_t (nlat,nlev,NM));     acc%mon_t  = 0.0_wp
        allocate(acc%mon_u (nlat,nlev,NM));     acc%mon_u  = 0.0_wp
        allocate(acc%mon_v (nlat,nlev,NM));     acc%mon_v  = 0.0_wp
        allocate(acc%mon_rh(nlat,nlev,NM));     acc%mon_rh = 0.0_wp

        allocate(acc%mon_sst(nlat,NM));         acc%mon_sst = 0.0_wp
        allocate(acc%mon_pr (nlat,NM));         acc%mon_pr  = 0.0_wp
        allocate(acc%mon_cnt(NM));              acc%mon_cnt = 0

        allocate(acc%mon_t2m (nlon,nlat,NM));   acc%mon_t2m  = 0.0_wp
        allocate(acc%mon_tlow(nlon,nlat,NM));   acc%mon_tlow = 0.0_wp
        allocate(acc%mon_pr2d(nlon,nlat,NM));   acc%mon_pr2d = 0.0_wp
        allocate(acc%mon_z500(nlon,nlat,NM));   acc%mon_z500 = 0.0_wp
        allocate(acc%mon_vlow(nlon,nlat,NM));   acc%mon_vlow = 0.0_wp
        allocate(acc%mon_evap(nlon,nlat,NM));   acc%mon_evap = 0.0_wp
        allocate(acc%mon_shf (nlon,nlat,NM));   acc%mon_shf  = 0.0_wp

        allocate(acc%mon_cf(nlat,nlev,NM));     acc%mon_cf = 0.0_wp
        allocate(acc%mon_skin (nlon,nlat,NM));  acc%mon_skin  = 0.0_wp
        allocate(acc%mon_lwd  (nlon,nlat,NM));  acc%mon_lwd   = 0.0_wp
        allocate(acc%mon_olr  (nlon,nlat,NM));  acc%mon_olr   = 0.0_wp
        allocate(acc%mon_olrcs(nlon,nlat,NM));  acc%mon_olrcs = 0.0_wp
        allocate(acc%mon_cover(nlon,nlat,NM));  acc%mon_cover = 0.0_wp

        allocate(acc%mon_vt (nlat,nlev,NM));    acc%mon_vt  = 0.0_wp
        allocate(acc%mon_vq (nlat,nlev,NM));    acc%mon_vq  = 0.0_wp
        allocate(acc%mon_eke(nlat,nlev,NM));    acc%mon_eke = 0.0_wp
        allocate(acc%mon_uv (nlat,nlev,NM));    acc%mon_uv  = 0.0_wp

        allocate(acc%mon_u3(nlon,nlat,nlev,NM));  acc%mon_u3 = 0.0_wp
        allocate(acc%mon_v3(nlon,nlat,nlev,NM));  acc%mon_v3 = 0.0_wp
        allocate(acc%mon_t3(nlon,nlat,nlev,NM));  acc%mon_t3 = 0.0_wp
        allocate(acc%mon_q3(nlon,nlat,nlev,NM));  acc%mon_q3 = 0.0_wp

        allocate(acc%mon_mr (nlat,nlev,NM));    acc%mon_mr  = 0.0_wp
        allocate(acc%mon_mrt(nlat,nlev,NM));    acc%mon_mrt = 0.0_wp
        allocate(acc%mon_mrq(nlat,nlev,NM));    acc%mon_mrq = 0.0_wp
        allocate(acc%mon_mr2(nlon,nlat,NM));    acc%mon_mr2 = 0.0_wp

        allocate(acc%mon_swn (nlon,nlat,NM));   acc%mon_swn  = 0.0_wp
        allocate(acc%mon_swup(nlon,nlat,NM));   acc%mon_swup = 0.0_wp

        return

    end subroutine aeros_diag_acc_init_monthly

    subroutine aeros_diag_acc_reset(acc)
        ! Zero every allocated accumulator and its counter, leaving allocation
        ! status untouched -- Wave 3 will call this at the end of a coupled
        ! cycle so the next cycle starts clean without a re-allocate.

        implicit none

        type(aeros_diag_acc_t), intent(inout) :: acc

        if (allocated(acc%w_tavg))    acc%w_tavg    = 0.0_wp
        if (allocated(acc%cnv_tavg))  acc%cnv_tavg  = 0.0_wp
        if (allocated(acc%cnd_tavg))  acc%cnd_tavg  = 0.0_wp
        if (allocated(acc%rad_tavg))  acc%rad_tavg  = 0.0_wp
        if (allocated(acc%surf_tavg)) acc%surf_tavg = 0.0_wp
        if (allocated(acc%vadv_tavg)) acc%vadv_tavg = 0.0_wp
        if (allocated(acc%dyn_tavg))  acc%dyn_tavg  = 0.0_wp
        if (allocated(acc%vdf_tavg))  acc%vdf_tavg  = 0.0_wp
        if (allocated(acc%bd_tavg))   acc%bd_tavg   = 0.0_wp
        acc%n_wacc = 0

        if (allocated(acc%uv_tavg))   acc%uv_tavg   = 0.0_wp
        if (allocated(acc%ubar_tavg)) acc%ubar_tavg = 0.0_wp
        if (allocated(acc%vbar_tavg)) acc%vbar_tavg = 0.0_wp

        if (allocated(acc%ke_spec))   acc%ke_spec = 0.0_wp
        if (allocated(acc%mf_spec))   acc%mf_spec = 0.0_wp
        if (allocated(acc%ke_latm))   acc%ke_latm = 0.0_wp
        if (allocated(acc%ke_levm))   acc%ke_levm = 0.0_wp
        ! ctab/stab are precomputed constants -- do NOT reset.

        if (allocated(acc%mon_t))    acc%mon_t    = 0.0_wp
        if (allocated(acc%mon_u))    acc%mon_u    = 0.0_wp
        if (allocated(acc%mon_v))    acc%mon_v    = 0.0_wp
        if (allocated(acc%mon_rh))   acc%mon_rh   = 0.0_wp
        if (allocated(acc%mon_sst))  acc%mon_sst  = 0.0_wp
        if (allocated(acc%mon_pr))   acc%mon_pr   = 0.0_wp
        if (allocated(acc%mon_cnt))  acc%mon_cnt  = 0
        if (allocated(acc%mon_t2m))  acc%mon_t2m  = 0.0_wp
        if (allocated(acc%mon_tlow)) acc%mon_tlow = 0.0_wp
        if (allocated(acc%mon_pr2d)) acc%mon_pr2d = 0.0_wp
        if (allocated(acc%mon_z500)) acc%mon_z500 = 0.0_wp
        if (allocated(acc%mon_vlow)) acc%mon_vlow = 0.0_wp
        if (allocated(acc%mon_evap)) acc%mon_evap = 0.0_wp
        if (allocated(acc%mon_shf))  acc%mon_shf  = 0.0_wp
        if (allocated(acc%mon_cf))   acc%mon_cf   = 0.0_wp
        if (allocated(acc%mon_skin)) acc%mon_skin = 0.0_wp
        if (allocated(acc%mon_lwd))  acc%mon_lwd  = 0.0_wp
        if (allocated(acc%mon_olr))  acc%mon_olr  = 0.0_wp
        if (allocated(acc%mon_olrcs))acc%mon_olrcs= 0.0_wp
        if (allocated(acc%mon_cover))acc%mon_cover= 0.0_wp
        if (allocated(acc%mon_vt))   acc%mon_vt   = 0.0_wp
        if (allocated(acc%mon_vq))   acc%mon_vq   = 0.0_wp
        if (allocated(acc%mon_eke))  acc%mon_eke  = 0.0_wp
        if (allocated(acc%mon_uv))   acc%mon_uv   = 0.0_wp
        if (allocated(acc%mon_u3))   acc%mon_u3   = 0.0_wp
        if (allocated(acc%mon_v3))   acc%mon_v3   = 0.0_wp
        if (allocated(acc%mon_t3))   acc%mon_t3   = 0.0_wp
        if (allocated(acc%mon_q3))   acc%mon_q3   = 0.0_wp
        if (allocated(acc%mon_mr))   acc%mon_mr   = 0.0_wp
        if (allocated(acc%mon_mrt))  acc%mon_mrt  = 0.0_wp
        if (allocated(acc%mon_mrq))  acc%mon_mrq  = 0.0_wp
        if (allocated(acc%mon_mr2))  acc%mon_mr2  = 0.0_wp
        if (allocated(acc%mon_swn))  acc%mon_swn  = 0.0_wp
        if (allocated(acc%mon_swup)) acc%mon_swup = 0.0_wp

        acc%wb_w_last = -1.0_wp
        acc%wb_p_acc  = 0.0_wp
        acc%wb_e_acc  = 0.0_wp
        acc%wb_n_acc  = 0

        acc%tb_swin  = 0.0_wp
        acc%tb_swup  = 0.0_wp
        acc%tb_olr   = 0.0_wp
        acc%tb_swupc = 0.0_wp
        acc%tb_olrc  = 0.0_wp
        acc%tb_cov   = 0.0_wp
        acc%tb_shf   = 0.0_wp
        acc%tb_lhf   = 0.0_wp
        acc%tb_n     = 0

        return

    end subroutine aeros_diag_acc_reset

    subroutine aeros_diag_acc_end(acc)
        ! Deallocate every buffer and clear the counters. Idempotent.

        implicit none

        type(aeros_diag_acc_t), intent(inout) :: acc

        if (allocated(acc%w_tavg))    deallocate(acc%w_tavg)
        if (allocated(acc%cnv_tavg))  deallocate(acc%cnv_tavg)
        if (allocated(acc%cnd_tavg))  deallocate(acc%cnd_tavg)
        if (allocated(acc%rad_tavg))  deallocate(acc%rad_tavg)
        if (allocated(acc%surf_tavg)) deallocate(acc%surf_tavg)
        if (allocated(acc%vadv_tavg)) deallocate(acc%vadv_tavg)
        if (allocated(acc%dyn_tavg))  deallocate(acc%dyn_tavg)
        if (allocated(acc%vdf_tavg))  deallocate(acc%vdf_tavg)
        if (allocated(acc%bd_tavg))   deallocate(acc%bd_tavg)
        if (allocated(acc%uv_tavg))   deallocate(acc%uv_tavg)
        if (allocated(acc%ubar_tavg)) deallocate(acc%ubar_tavg)
        if (allocated(acc%vbar_tavg)) deallocate(acc%vbar_tavg)
        if (allocated(acc%ke_spec))   deallocate(acc%ke_spec)
        if (allocated(acc%mf_spec))   deallocate(acc%mf_spec)
        if (allocated(acc%ctab))      deallocate(acc%ctab)
        if (allocated(acc%stab))      deallocate(acc%stab)
        if (allocated(acc%ke_latm))   deallocate(acc%ke_latm)
        if (allocated(acc%ke_levm))   deallocate(acc%ke_levm)

        if (allocated(acc%mon_t))    deallocate(acc%mon_t)
        if (allocated(acc%mon_u))    deallocate(acc%mon_u)
        if (allocated(acc%mon_v))    deallocate(acc%mon_v)
        if (allocated(acc%mon_rh))   deallocate(acc%mon_rh)
        if (allocated(acc%mon_sst))  deallocate(acc%mon_sst)
        if (allocated(acc%mon_pr))   deallocate(acc%mon_pr)
        if (allocated(acc%mon_cnt))  deallocate(acc%mon_cnt)
        if (allocated(acc%mon_t2m))  deallocate(acc%mon_t2m)
        if (allocated(acc%mon_tlow)) deallocate(acc%mon_tlow)
        if (allocated(acc%mon_pr2d)) deallocate(acc%mon_pr2d)
        if (allocated(acc%mon_z500)) deallocate(acc%mon_z500)
        if (allocated(acc%mon_vlow)) deallocate(acc%mon_vlow)
        if (allocated(acc%mon_evap)) deallocate(acc%mon_evap)
        if (allocated(acc%mon_shf))  deallocate(acc%mon_shf)
        if (allocated(acc%mon_cf))   deallocate(acc%mon_cf)
        if (allocated(acc%mon_skin)) deallocate(acc%mon_skin)
        if (allocated(acc%mon_lwd))  deallocate(acc%mon_lwd)
        if (allocated(acc%mon_olr))  deallocate(acc%mon_olr)
        if (allocated(acc%mon_olrcs))deallocate(acc%mon_olrcs)
        if (allocated(acc%mon_cover))deallocate(acc%mon_cover)
        if (allocated(acc%mon_vt))   deallocate(acc%mon_vt)
        if (allocated(acc%mon_vq))   deallocate(acc%mon_vq)
        if (allocated(acc%mon_eke))  deallocate(acc%mon_eke)
        if (allocated(acc%mon_uv))   deallocate(acc%mon_uv)
        if (allocated(acc%mon_u3))   deallocate(acc%mon_u3)
        if (allocated(acc%mon_v3))   deallocate(acc%mon_v3)
        if (allocated(acc%mon_t3))   deallocate(acc%mon_t3)
        if (allocated(acc%mon_q3))   deallocate(acc%mon_q3)
        if (allocated(acc%mon_mr))   deallocate(acc%mon_mr)
        if (allocated(acc%mon_mrt))  deallocate(acc%mon_mrt)
        if (allocated(acc%mon_mrq))  deallocate(acc%mon_mrq)
        if (allocated(acc%mon_mr2))  deallocate(acc%mon_mr2)
        if (allocated(acc%mon_swn))  deallocate(acc%mon_swn)
        if (allocated(acc%mon_swup)) deallocate(acc%mon_swup)

        acc%n_wacc    = 0
        acc%mmax_spec = 0
        acc%nacc_mon  = 1

        acc%wb_w_last = -1.0_wp
        acc%wb_p_acc  = 0.0_wp
        acc%wb_e_acc  = 0.0_wp
        acc%wb_n_acc  = 0

        acc%tb_swin  = 0.0_wp
        acc%tb_swup  = 0.0_wp
        acc%tb_olr   = 0.0_wp
        acc%tb_swupc = 0.0_wp
        acc%tb_olrc  = 0.0_wp
        acc%tb_cov   = 0.0_wp
        acc%tb_shf   = 0.0_wp
        acc%tb_lhf   = 0.0_wp
        acc%tb_n     = 0

        acc%spinup_second_half_nstep = 0
        acc%spinup_monthly_days      = 0.0_wp

        return

    end subroutine aeros_diag_acc_end

    ! === Wave 4c accumulator arithmetic ======================================
    !
    ! Every one of these methods mirrors, line-for-line, the block that used to
    ! live at rce_long program scope, and reads state directly off `ts`, `now`,
    ! `grd`, `vgrid`, `par`, `phys`. Bit-for-bit on the rce_long regression:
    ! the port is a code move, not a rewrite.
    !
    ! The driver still owns the WHEN (spin-up threshold, report cadence, which
    ! calendar month bin, whether to accumulate at all). The library owns the
    ! WHAT (the sums). aeros_lib_step on the coupled path uses the same
    ! methods to fill ams%diag_acc for the export bundle.

    integer function aeros_diag_acc_imon(doy) result(im)
        ! Bin a day-of-year into 1..NMON using the 365-day radiative calendar
        ! rce_long uses (aeros_insol_day wrap). frac is the 0-based day within
        ! the year; im = floor(frac / (365/12)) + 1, clipped to NMON.

        implicit none

        real(wp), intent(in) :: doy
        real(wp), parameter :: DAYS_YR = 365.0_wp
        real(wp) :: frac

        frac = modulo(doy - 1.0_wp, DAYS_YR)
        im   = min(AEROS_DIAG_NMON, int(frac/(DAYS_YR/real(AEROS_DIAG_NMON,wp))) + 1)

        return

    end function aeros_diag_acc_imon

    real(wp) function gmean2_local(f, grd) result(m)
        ! Area-weighted global mean of a 2D field, verbatim the rce_long
        ! local gmean2. Kept in this module (private) so callers do not need
        ! to re-implement it.
        implicit none
        real(wp),               intent(in) :: f(:,:)
        type(aeros_grid_class), intent(in) :: grd
        real(wp) :: s, wsum
        integer  :: i, j
        s = 0.0_wp; wsum = 0.0_wp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                s    = s    + f(i,j)*grd%area(i,j)
                wsum = wsum + grd%area(i,j)
            end do
        end do
        m = s/wsum
    end function gmean2_local

    subroutine aeros_diag_acc_water(acc, ts, phys, grd)
        ! Per-step accumulator for the global-mean precip / evap water budget
        ! (rce_long inline block above the main step). Only accumulates when
        ! l_surf and both precip fields are allocated -- same guard rce_long
        ! used.

        implicit none

        type(aeros_diag_acc_t),      intent(inout) :: acc
        type(aeros_timestep_class),  intent(in)    :: ts
        type(aeros_physics_class),   intent(in)    :: phys
        type(aeros_grid_class),      intent(in)    :: grd

        real(wp), allocatable :: tmp(:,:)

        if (.not. phys%l_surf) return
        if (.not. (allocated(ts%cnv%precip) .and. allocated(ts%cnd%precip))) return

        allocate(tmp(grd%nlon,grd%nlat))
        tmp = (ts%cnv%precip + ts%cnd%precip)*86400.0_wp
        acc%wb_p_acc = acc%wb_p_acc + gmean2_local(tmp, grd)
        tmp = ts%surf%evap*86400.0_wp
        acc%wb_e_acc = acc%wb_e_acc + gmean2_local(tmp, grd)
        deallocate(tmp)
        acc%wb_n_acc = acc%wb_n_acc + 1

        return

    end subroutine aeros_diag_acc_water

    subroutine aeros_diag_acc_energy(acc, ts, phys, grd)
        ! Per-step accumulator for the global-mean TOA + surface energy
        ! budget (rce_long energy_accumulate, verbatim).

        implicit none

        type(aeros_diag_acc_t),      intent(inout) :: acc
        type(aeros_timestep_class),  intent(in)    :: ts
        type(aeros_physics_class),   intent(in)    :: phys
        type(aeros_grid_class),      intent(in)    :: grd

        real(wp) :: swin, wsum
        integer  :: j

        if (.not. (phys%l_rad .or. phys%l_surf)) return
        if (phys%l_rad) then
            acc%tb_olr   = acc%tb_olr   + gmean2_local(ts%rad%olr,          grd)
            acc%tb_swup  = acc%tb_swup  + gmean2_local(ts%rad%sw_up_toa,    grd)
            acc%tb_olrc  = acc%tb_olrc  + gmean2_local(ts%rad%olr_cs,       grd)
            acc%tb_swupc = acc%tb_swupc + gmean2_local(ts%rad%sw_up_toa_cs, grd)
            acc%tb_cov   = acc%tb_cov   + gmean2_local(ts%rad%cover,        grd)
            swin = 0.0_wp; wsum = 0.0_wp
            do j = 1, grd%nlat
                swin = swin + ts%rad%sw_toa(j)*sum(grd%area(:,j))
                wsum = wsum + sum(grd%area(:,j))
            end do
            acc%tb_swin = acc%tb_swin + swin/wsum
        end if
        if (phys%l_surf) then
            acc%tb_shf = acc%tb_shf + gmean2_local(ts%surf%shf, grd)
            acc%tb_lhf = acc%tb_lhf + gmean2_local(ts%surf%lhf, grd)
        end if
        acc%tb_n = acc%tb_n + 1

        return

    end subroutine aeros_diag_acc_energy

    subroutine aeros_diag_acc_second_half(acc, ts, now, vgrid, grd, par)
        ! Per-step 2nd-half tavg accumulator: omega, the forward-split heating
        ! split, the eddy zonal-mean momentum flux, the eddy wavenumber
        ! spectrum. Bumps n_wacc. Called by the driver only when its
        ! spin-up policy says "we are in the equilibrated 2nd half".

        implicit none

        type(aeros_diag_acc_t),      intent(inout) :: acc
        type(aeros_timestep_class),  intent(in)    :: ts
        type(aeros_state_class),     intent(in)    :: now
        type(aeros_vgrid_class),     intent(in)    :: vgrid
        type(aeros_grid_class),      intent(in)    :: grd
        type(aeros_param_class),     intent(in)    :: par

        integer  :: i, j, k, m
        real(wp) :: uzm, vzm, up, vp, cov, rn
        real(wp) :: w, dsig, ur, ui, vr, vi, c, s, kem

        ! The per-term heating split arrays are only allocated when the
        ! driver has called aeros_timestep_enable_diag. Without them the
        ! 2nd-half tavg is meaningless -- skip it.
        if (.not. (allocated(ts%wrk%dt_cnv) .and. allocated(ts%wrk%omega))) return
        ! Wave 4d driver spin-up gate: rce_long's `n > n0 + nstep/2` moved
        ! here so aeros_lib_step can call the method every sub-step and the
        ! WHEN stays a driver choice (default 0 = coupled host, accumulate
        ! from step 1).
        if (ts%nstep <= acc%spinup_second_half_nstep) return

        acc%w_tavg    = acc%w_tavg    + ts%wrk%omega
        acc%cnv_tavg  = acc%cnv_tavg  + ts%wrk%dt_cnv
        acc%cnd_tavg  = acc%cnd_tavg  + ts%wrk%dt_cnd
        acc%rad_tavg  = acc%rad_tavg  + ts%wrk%dt_rad
        acc%surf_tavg = acc%surf_tavg + ts%wrk%dt_surf
        acc%vadv_tavg = acc%vadv_tavg + ts%wrk%dt_vadv*par%dt
        acc%dyn_tavg  = acc%dyn_tavg  + ts%wrk%dtdt*par%dt
        acc%vdf_tavg  = acc%vdf_tavg  + ts%wrk%dt_vdiff
        acc%bd_tavg   = acc%bd_tavg   + ts%wrk%dt_bd

        ! Zonal-mean jet / eddy momentum flux (accum_zmflux).
        rn = real(grd%nlon, wp)
        do k = 1, par%nlev
            do j = 1, grd%nlat
                uzm = sum(now%u(:,j,k))/rn
                vzm = sum(now%v(:,j,k))/rn
                cov = 0.0_wp
                do i = 1, grd%nlon
                    up  = now%u(i,j,k) - uzm
                    vp  = now%v(i,j,k) - vzm
                    cov = cov + up*vp
                end do
                acc%ubar_tavg(j,k) = acc%ubar_tavg(j,k) + uzm
                acc%vbar_tavg(j,k) = acc%vbar_tavg(j,k) + vzm
                acc%uv_tavg(j,k)   = acc%uv_tavg(j,k)   + cov/rn
            end do
        end do

        ! Eddy wavenumber spectrum (accum_espec).
        if (acc%mmax_spec > 0) then
            do k = 1, par%nlev
                dsig = vgrid%sigma_half(k) - vgrid%sigma_half(k-1)
                do j = 1, grd%nlat
                    w   = grd%area(1,j)*dsig
                    uzm = sum(now%u(:,j,k))/rn
                    vzm = sum(now%v(:,j,k))/rn
                    do m = 1, acc%mmax_spec
                        ur = 0.0_wp; ui = 0.0_wp; vr = 0.0_wp; vi = 0.0_wp
                        do i = 1, grd%nlon
                            up = now%u(i,j,k) - uzm
                            vp = now%v(i,j,k) - vzm
                            c  = acc%ctab(m,i); s = acc%stab(m,i)
                            ur = ur + up*c; ui = ui + up*s
                            vr = vr + vp*c; vi = vi + vp*s
                        end do
                        ur = ur/rn; ui = ui/rn; vr = vr/rn; vi = vi/rn
                        kem = ur*ur + ui*ui + vr*vr + vi*vi
                        acc%ke_spec(m)   = acc%ke_spec(m)   + w*kem
                        acc%ke_latm(j,m) = acc%ke_latm(j,m) + dsig*kem
                        acc%ke_levm(k,m) = acc%ke_levm(k,m) + grd%area(1,j)*kem
                        acc%mf_spec(m)   = acc%mf_spec(m)   + w*2.0_wp*(ur*vr + ui*vi)
                    end do
                end do
            end do
        end if

        acc%n_wacc = acc%n_wacc + 1

        return

    end subroutine aeros_diag_acc_second_half

    subroutine aeros_diag_acc_monthly(acc, ts, now, vgrid, grd, par, phys, imon)
        ! Per-daily-sample monthly-climatology accumulator (rce_long
        ! accum_monthly, verbatim). Reads the current phis off ts%wrk%phis
        ! (the field aeros_timestep_set_phis last wrote), so the driver's
        ! topography ramp is reflected without a separate argument.

        implicit none

        type(aeros_diag_acc_t),      intent(inout) :: acc
        type(aeros_timestep_class),  intent(in)    :: ts
        type(aeros_state_class),     intent(in)    :: now
        type(aeros_vgrid_class),     intent(in)    :: vgrid
        type(aeros_grid_class),      intent(in)    :: grd
        type(aeros_param_class),     intent(in)    :: par
        type(aeros_physics_class),   intent(in)    :: phys
        integer,                     intent(in)    :: imon

        real(wp), parameter :: P500 = 50000.0_wp
        real(wp) :: phc(0:par%nlev), pfc(par%nlev), dpcc(par%nlev)
        real(wp) :: qsl, dqsl, rnlon, pr_ij, wgt
        integer  :: ii, jj, kk
        real(wp) :: tzm(grd%nlat,par%nlev), uzm(grd%nlat,par%nlev), vzm(grd%nlat,par%nlev), rzm(grd%nlat,par%nlev)
        real(wp) :: sstz(grd%nlat), prz(grd%nlat)
        real(wp) :: phif(par%nlev), phih, alph
        real(wp) :: cfz(grd%nlat,par%nlev), cfc(par%nlev), clwc(par%nlev), ciwc(par%nlev)
        real(wp) :: qzm(grd%nlat,par%nlev), vtzm(grd%nlat,par%nlev), vqzm(grd%nlat,par%nlev), uvzm(grd%nlat,par%nlev)
        real(wp) :: u2zm(grd%nlat,par%nlev), v2zm(grd%nlat,par%nlev)
        real(wp) :: mrzm(grd%nlat,par%nlev), mrtzm(grd%nlat,par%nlev), mrqzm(grd%nlat,par%nlev), mr
        logical  :: have_pr

        ! Guard against subsystems that never allocated their diagnostic
        ! fields (dry Held-Suarez has no moisture, no radiation cache);
        ! the monthly bin is not meaningful without them.
        if (.not. allocated(ts%mst%mass_ratio)) return
        if (.not. allocated(ts%ocn%sst))        return
        ! Wave 4d driver spin-up gate: rce_long's `elapsed_days >= mon_spinup_days`
        ! moved here (default 0 = coupled host, no spin-up).
        if (real(ts%nstep, wp)*par%dt/86400.0_wp < acc%spinup_monthly_days) return

        have_pr = allocated(ts%cnv%precip) .and. allocated(ts%cnd%precip)
        rnlon = real(grd%nlon, wp)
        tzm = 0.0_wp; uzm = 0.0_wp; vzm = 0.0_wp; rzm = 0.0_wp
        sstz = 0.0_wp; prz = 0.0_wp; cfz = 0.0_wp
        qzm = 0.0_wp; vtzm = 0.0_wp; vqzm = 0.0_wp; u2zm = 0.0_wp; v2zm = 0.0_wp; uvzm = 0.0_wp
        mrzm = 0.0_wp; mrtzm = 0.0_wp; mrqzm = 0.0_wp
        do jj = 1, grd%nlat
            do ii = 1, grd%nlon
                call aeros_vgrid_pressure(vgrid, now%ps(ii,jj), phc, pfc, dpcc)
                if (ts%cpr%enabled) then
                    cfc = now%cf_g(ii,jj,:)
                else
                    call aeros_cloud_diagnose(par%nlev, now%temp_g(ii,jj,:), &
                        now%qv_g(ii,jj,:), pfc, now%ps(ii,jj), cfc, clwc, ciwc)
                end if
                cfz(jj,:) = cfz(jj,:) + cfc
                do kk = 1, par%nlev
                    call aeros_qsat(now%temp_g(ii,jj,kk), pfc(kk), qsl, dqsl)
                    tzm(jj,kk) = tzm(jj,kk) + now%temp_g(ii,jj,kk)
                    uzm(jj,kk) = uzm(jj,kk) + now%u(ii,jj,kk)
                    vzm(jj,kk) = vzm(jj,kk) + now%v(ii,jj,kk)
                    rzm(jj,kk) = rzm(jj,kk) + &
                        now%qv_g(ii,jj,kk)/max(qsl,1.0e-12_wp)*100.0_wp
                    qzm(jj,kk)  = qzm(jj,kk)  + now%qv_g(ii,jj,kk)
                    vtzm(jj,kk) = vtzm(jj,kk) + now%v(ii,jj,kk)*now%temp_g(ii,jj,kk)
                    uvzm(jj,kk) = uvzm(jj,kk) + now%u(ii,jj,kk)*now%v(ii,jj,kk)
                    vqzm(jj,kk) = vqzm(jj,kk) + now%v(ii,jj,kk)*now%qv_g(ii,jj,kk)
                    u2zm(jj,kk) = u2zm(jj,kk) + now%u(ii,jj,kk)**2
                    v2zm(jj,kk) = v2zm(jj,kk) + now%v(ii,jj,kk)**2
                    mr = ts%mst%mass_ratio(ii,jj,kk) - 1.0_wp
                    mrzm(jj,kk)  = mrzm(jj,kk)  + mr
                    mrtzm(jj,kk) = mrtzm(jj,kk) + mr*now%temp_g(ii,jj,kk)
                    mrqzm(jj,kk) = mrqzm(jj,kk) + mr*now%qv_g(ii,jj,kk)
                end do
                sstz(jj) = sstz(jj) + ts%ocn%sst(ii,jj)
                if (have_pr) then
                    pr_ij = (ts%cnv%precip(ii,jj) + ts%cnd%precip(ii,jj))*86400.0_wp
                    prz(jj) = prz(jj) + pr_ij
                    acc%mon_pr2d(ii,jj,imon) = acc%mon_pr2d(ii,jj,imon) + pr_ij
                end if
                phih = ts%wrk%phis(ii,jj)
                do kk = par%nlev, 1, -1
                    if (kk == 1) then
                        alph = log(2.0_wp)
                    else
                        alph = 1.0_wp - phc(kk-1)/dpcc(kk)*log(phc(kk)/phc(kk-1))
                    end if
                    phif(kk) = phih + real(R_d,wp)*now%temp_g(ii,jj,kk)*alph
                    if (kk > 1) phih = phih + real(R_d,wp)*now%temp_g(ii,jj,kk)*log(phc(kk)/phc(kk-1))
                end do
                kk = par%nlev - 1
                do while (kk > 1 .and. pfc(kk) > P500)
                    kk = kk - 1
                end do
                wgt = (log(P500) - log(pfc(kk)))/(log(pfc(kk+1)) - log(pfc(kk)))
                acc%mon_z500(ii,jj,imon) = acc%mon_z500(ii,jj,imon) + &
                    (phif(kk) + wgt*(phif(kk+1) - phif(kk)))/real(grav,wp)
            end do
        end do
        acc%mon_t(:,:,imon)  = acc%mon_t(:,:,imon)  + tzm/rnlon
        acc%mon_u(:,:,imon)  = acc%mon_u(:,:,imon)  + uzm/rnlon
        acc%mon_v(:,:,imon)  = acc%mon_v(:,:,imon)  + vzm/rnlon
        acc%mon_rh(:,:,imon) = acc%mon_rh(:,:,imon) + rzm/rnlon
        acc%mon_sst(:,imon)  = acc%mon_sst(:,imon)  + sstz/rnlon
        acc%mon_pr(:,imon)   = acc%mon_pr(:,imon)   + prz/rnlon
        acc%mon_cf(:,:,imon) = acc%mon_cf(:,:,imon) + cfz/rnlon
        acc%mon_vt(:,:,imon)  = acc%mon_vt(:,:,imon)  + (vtzm - vzm*tzm/rnlon)/rnlon
        acc%mon_uv(:,:,imon)  = acc%mon_uv(:,:,imon)  + (uvzm - uzm*vzm/rnlon)/rnlon
        acc%mon_vq(:,:,imon)  = acc%mon_vq(:,:,imon)  + (vqzm - vzm*qzm/rnlon)/rnlon
        acc%mon_eke(:,:,imon) = acc%mon_eke(:,:,imon) + 0.5_wp*(u2zm + v2zm - (uzm**2 + vzm**2)/rnlon)/rnlon
        acc%mon_mr(:,:,imon)  = acc%mon_mr(:,:,imon)  + mrzm/rnlon
        acc%mon_mrt(:,:,imon) = acc%mon_mrt(:,:,imon) + (mrtzm - mrzm*tzm/rnlon)/rnlon
        acc%mon_mrq(:,:,imon) = acc%mon_mrq(:,:,imon) + (mrqzm - mrzm*qzm/rnlon)/rnlon
        acc%mon_mr2(:,:,imon) = acc%mon_mr2(:,:,imon) + ts%mst%mass_ratio(:,:,par%nlev) - 1.0_wp
        acc%mon_u3(:,:,:,imon) = acc%mon_u3(:,:,:,imon) + now%u(:,:,1:par%nlev)
        acc%mon_v3(:,:,:,imon) = acc%mon_v3(:,:,:,imon) + now%v(:,:,1:par%nlev)
        acc%mon_t3(:,:,:,imon) = acc%mon_t3(:,:,:,imon) + now%temp_g(:,:,1:par%nlev)
        acc%mon_q3(:,:,:,imon) = acc%mon_q3(:,:,:,imon) + now%qv_g(:,:,1:par%nlev)
        if (ts%land%enabled .and. allocated(ts%land%skin)) then
            acc%mon_skin(:,:,imon) = acc%mon_skin(:,:,imon) + ts%land%skin
        else
            acc%mon_skin(:,:,imon) = acc%mon_skin(:,:,imon) + ts%ocn%sst
        end if
        if (allocated(ts%rad%olr)) then
            acc%mon_lwd(:,:,imon)   = acc%mon_lwd(:,:,imon)   + ts%rad%lw_dw_sur
            acc%mon_olr(:,:,imon)   = acc%mon_olr(:,:,imon)   + ts%rad%olr
            acc%mon_olrcs(:,:,imon) = acc%mon_olrcs(:,:,imon) + ts%rad%olr_cs
            acc%mon_cover(:,:,imon) = acc%mon_cover(:,:,imon) + ts%rad%cover
            acc%mon_swn(:,:,imon)   = acc%mon_swn(:,:,imon)   + ts%rad%sw_net_sur
            acc%mon_swup(:,:,imon)  = acc%mon_swup(:,:,imon)  + ts%rad%sw_up_toa
        end if
        acc%mon_tlow(:,:,imon) = acc%mon_tlow(:,:,imon) + now%temp_g(:,:,par%nlev)
        if (ts%surf%enabled) then
            acc%mon_t2m(:,:,imon) = acc%mon_t2m(:,:,imon) + ts%surf%t2m
        else
            acc%mon_t2m(:,:,imon) = acc%mon_t2m(:,:,imon) + now%temp_g(:,:,par%nlev)
        end if
        acc%mon_vlow(:,:,imon) = acc%mon_vlow(:,:,imon) + sqrt(now%u(:,:,par%nlev)**2 + now%v(:,:,par%nlev)**2)
        if (allocated(ts%surf%evap)) then
            acc%mon_evap(:,:,imon) = acc%mon_evap(:,:,imon) + ts%surf%evap*86400.0_wp
            acc%mon_shf(:,:,imon)  = acc%mon_shf(:,:,imon)  + ts%surf%shf
        end if
        acc%mon_cnt(imon) = acc%mon_cnt(imon) + 1

        ! phys is passed for future gating (unused today; keep signature stable).
        if (phys%l_rad .neqv. phys%l_rad) return

        return

    end subroutine aeros_diag_acc_monthly

end module aeros_diag_acc
