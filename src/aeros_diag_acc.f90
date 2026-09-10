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

    use aeros_defs, only : wp

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

        return

    end subroutine aeros_diag_acc_end

end module aeros_diag_acc
