!==============================================================================
! aeros_bench  —  M0a cost / OpenMP-scaling harness for the aeros core (Option A)
!------------------------------------------------------------------------------
! This is a BENCHMARK, not a dynamical core. It reproduces the *transform +
! nonlinear + implicit-solve* workload of one semi-implicit spectral
! primitive-equation timestep so we can measure, on the target node:
!
!   (1) real core-seconds per model-year at T31L16 and T42L19,
!   (2) the share of runtime spent in spherical-harmonic transforms,
!   (3) the OpenMP scaling curve to many threads.
!
! Physics (radiation, condensation, surface) is deliberately excluded — per the
! design doc §3.6 it is ~1% of the budget. Prognostic fields are held fixed at
! their (non-trivial) initial amplitude across all timed steps: the benchmark
! needs *representative arithmetic on representative data* each step, not a
! stable time integration. A checksum of the tendencies is accumulated and
! printed to defeat dead-code elimination.
!
! Precision: the SHTns OpenMP library exposes only C_DOUBLE / C_DOUBLE_COMPLEX
! interfaces (see shtns.f03), so the transforms — and hence this whole harness —
! run in DOUBLE precision. The design wants Float32; a single-precision core
! would be roughly ~2x faster, so the numbers here are a conservative UPPER
! BOUND on cost.
!
! Grid: full Gaussian grid (SHT_GAUSS). The SHTns basic API has no reduced /
! octahedral grid, so gridpoint work here is an UPPER BOUND — a reduced grid
! removes ~1/3 of the gridpoints (design §4).
!==============================================================================
program aeros_bench
   use, intrinsic :: iso_c_binding
   use omp_lib
   implicit none

   include "shtns.f03"      ! SHTns C-bound interfaces, enums, and shtns_info type

   integer,  parameter :: dp  = c_double            ! spatial reals
   integer,  parameter :: dpc = c_double_complex    ! spectral coeffs

   ! ---- namelist inputs -----------------------------------------------------
   integer            :: lmax, nlat, nphi, nlev, nstep, nwarm
   real(dp)           :: dt
   character(len=64)  :: label
   namelist /bench/ lmax, nlat, nphi, nlev, dt, nstep, nwarm, label

   ! ---- SHTns handle and grid metadata --------------------------------------
   type(c_ptr)                 :: sht
   type(shtns_info), pointer   :: si
   integer                     :: mmax, nlm, nspat, nlat_r, nphi_r
   integer                     :: nthreads, threads_used

   ! ---- prognostic spectral fields (held fixed across timed steps) ----------
   complex(dpc), allocatable :: vor(:,:), div(:,:), tmp(:,:), qhu(:,:)  ! (nlm,nlev)
   complex(dpc), allocatable :: lnps(:)                                 ! (nlm)
   ! ---- tendency spectral fields (recomputed every step) --------------------
   complex(dpc), allocatable :: dvor(:,:), ddiv(:,:), dtmp(:,:), dqhu(:,:)
   complex(dpc), allocatable :: dlnps(:)

   ! ---- spatial (gridpoint) work arrays, all (nspat,nlev) unless noted -------
   real(dp), allocatable :: ug(:,:), vg(:,:)          ! wind components
   real(dp), allocatable :: vorg(:,:), divg(:,:)      ! rel. vorticity, divergence
   real(dp), allocatable :: tg(:,:),  qg(:,:)         ! temperature, humidity
   real(dp), allocatable :: Fu(:,:),  Fv(:,:)         ! momentum tendency (spheroidal/toroidal)
   real(dp), allocatable :: FT(:,:),  Fq(:,:)         ! thermodynamic / moisture tendency
   real(dp), allocatable :: lnpsg(:), Flnps(:)        ! (nspat) log-surface-pressure
   real(dp), allocatable :: fcor(:)                   ! (nspat) Coriolis (representative)

   ! ---- implicit vertical structure + hyperdiffusion ------------------------
   real(dp), allocatable :: asub(:), bdiag(:), csup(:)   ! (nlev) tridiag coeffs
   real(dp), allocatable :: damp(:)                      ! (nlm) l-dependent damping
   integer,  allocatable :: lval(:)                      ! (nlm) degree l of each mode

   ! ---- timing --------------------------------------------------------------
   real(dp) :: t_transform, t_gridpoint, t_solve
   real(dp) :: t0, t_wall0, t_wall, t_total
   real(dp) :: t_step_s, phase_sum
   real(dp) :: fr_tr, fr_gp, fr_sv
   real(dp) :: checksum

   ! ---- reporting -----------------------------------------------------------
   real(dp) :: steps_per_year, wall_per_year, core_per_year, sypd
   character(len=512) :: nmlpath
   character(len=8)   :: precision_str = "double"
   integer :: k, istep, ios, u

   ! ==========================================================================
   ! 1. Read namelist
   ! ==========================================================================
   ! defaults
   lmax = 31; nlat = 48; nphi = 96; nlev = 16
   dt = 1800.0_dp; nstep = 100; nwarm = 10; label = "unnamed"

   call get_command_argument(1, nmlpath, status=ios)
   if (ios /= 0 .or. len_trim(nmlpath) == 0) then
      write(*,*) "ERROR: usage: aeros_bench.x <bench.nml>"
      stop 1
   end if
   open(newunit=u, file=trim(nmlpath), status="old", action="read", iostat=ios)
   if (ios /= 0) then
      write(*,*) "ERROR: cannot open namelist: ", trim(nmlpath)
      stop 1
   end if
   read(u, nml=bench, iostat=ios)
   if (ios /= 0) then
      write(*,*) "ERROR: failed to read &bench group from ", trim(nmlpath)
      stop 1
   end if
   close(u)
   mmax = lmax   ! triangular truncation, mres = 1

   ! ==========================================================================
   ! 2. SHTns setup.  Thread count comes from OMP_NUM_THREADS (set by runme),
   !    NOT from the namelist.  shtns_use_threads must precede grid init.
   ! ==========================================================================
   nthreads = omp_get_max_threads()
   threads_used = shtns_use_threads(nthreads)
   sht = shtns_create(lmax, mmax, 1, SHT_ORTHONORMAL)
   call shtns_set_grid(sht, SHT_GAUSS, 1.0e-10_dp, nlat, nphi)

   ! Read grid metadata straight from the SHTns struct to confirm the grid.
   call c_f_pointer(sht, si)
   nlm    = si%nlm
   nlat_r = si%nlat
   nphi_r = si%nphi
   nspat  = si%nspat

   ! ==========================================================================
   ! 3. Allocate
   ! ==========================================================================
   allocate(vor(nlm,nlev), div(nlm,nlev), tmp(nlm,nlev), qhu(nlm,nlev), lnps(nlm))
   allocate(dvor(nlm,nlev), ddiv(nlm,nlev), dtmp(nlm,nlev), dqhu(nlm,nlev), dlnps(nlm))
   allocate(ug(nspat,nlev), vg(nspat,nlev), vorg(nspat,nlev), divg(nspat,nlev))
   allocate(tg(nspat,nlev), qg(nspat,nlev))
   allocate(Fu(nspat,nlev), Fv(nspat,nlev), FT(nspat,nlev), Fq(nspat,nlev))
   allocate(lnpsg(nspat), Flnps(nspat), fcor(nspat))
   allocate(asub(nlev), bdiag(nlev), csup(nlev), damp(nlm), lval(nlm))

   ! ==========================================================================
   ! 4. Initialize
   ! ==========================================================================
   call build_lval(lmax, mmax, nlm, lval)      ! degree l of each spectral mode
   call init_spectral(nlm, nlev, lmax, lval, vor, div, tmp, qhu, lnps)

   ! representative Coriolis field: smooth, nonzero. Its geophysical placement
   ! is irrelevant to timing; only that it is real O(1e-4) data feeding the
   ! gridpoint arithmetic.
   block
      integer :: i, ilat
      real(dp), parameter :: omega = 7.292e-5_dp, pi = 3.141592653589793_dp
      real(dp) :: theta
      do i = 1, nspat
         ilat  = mod(i-1, nlat) + 1
         theta = pi * (real(ilat,dp) - 0.5_dp) / real(nlat,dp)
         fcor(i) = 2.0_dp * omega * cos(theta)
      end do
   end block

   ! implicit vertical-structure tridiagonal (diagonally dominant, constant
   ! across modes) — stands in for the gravity-wave / vertical-mode solve.
   asub  = -1.0_dp
   csup  = -1.0_dp
   bdiag =  2.5_dp

   ! spectral hyperdiffusion: multiply mode (l,m) by exp(-nu*(l(l+1))^2).
   block
      integer :: lm
      real(dp) :: nu, lnorm
      nu = 1.0e-2_dp
      lnorm = real(lmax*(lmax+1), dp)
      do lm = 1, nlm
         damp(lm) = exp( -nu * ( real(lval(lm)*(lval(lm)+1),dp)/lnorm )**2 )
      end do
   end block

   t_transform = 0.0_dp; t_gridpoint = 0.0_dp; t_solve = 0.0_dp
   checksum = 0.0_dp

   ! ==========================================================================
   ! 5. Warmup (untimed) then timed loop
   ! ==========================================================================
   do istep = 1, nwarm
      call do_step()
   end do

   t_transform = 0.0_dp; t_gridpoint = 0.0_dp; t_solve = 0.0_dp
   checksum = 0.0_dp
   t_wall0 = omp_get_wtime()
   do istep = 1, nstep
      call do_step()
   end do
   t_wall = omp_get_wtime() - t_wall0

   ! ==========================================================================
   ! 6. Derived metrics
   ! ==========================================================================
   phase_sum = t_transform + t_gridpoint + t_solve
   t_step_s  = phase_sum / real(nstep, dp)         ! per-step, sum of instrumented phases
   fr_tr = t_transform / phase_sum
   fr_gp = t_gridpoint / phase_sum
   fr_sv = t_solve     / phase_sum

   steps_per_year = 365.0_dp * 86400.0_dp / dt
   wall_per_year  = t_step_s * steps_per_year
   core_per_year  = wall_per_year * real(threads_used, dp)
   sypd           = 86400.0_dp / wall_per_year

   ! ==========================================================================
   ! 7. Report — human block + machine-parseable key = value lines, to both
   !    stdout and results.txt in the current working directory.
   ! ==========================================================================
   call report(6)                                   ! stdout
   open(newunit=u, file="results.txt", status="replace", action="write")
   call report(u)
   close(u)

   call shtns_destroy(sht)

contains

   !---------------------------------------------------------------------------
   ! One timestep: inverse transforms -> gridpoint nonlinear -> forward
   ! transforms -> implicit solve.  Prognostics stay fixed; tendencies are
   ! (re)written each call.
   !---------------------------------------------------------------------------
   subroutine do_step()
      integer :: kk

      ! --- (1) inverse transforms (t_transform) ------------------------------
      t0 = omp_get_wtime()
      do kk = 1, nlev
         ! vector transform: (spheroidal=div, toroidal=vor) -> (Vt=v, Vp=u)
         call SHsphtor_to_spat(sht, div(:,kk), vor(:,kk), vg(:,kk), ug(:,kk))
         call SH_to_spat(sht, vor(:,kk), vorg(:,kk))
         call SH_to_spat(sht, div(:,kk), divg(:,kk))
         call SH_to_spat(sht, tmp(:,kk), tg(:,kk))
         call SH_to_spat(sht, qhu(:,kk), qg(:,kk))
      end do
      call SH_to_spat(sht, lnps, lnpsg)
      t_transform = t_transform + (omp_get_wtime() - t0)

      ! --- (2) gridpoint nonlinear pass (t_gridpoint), OpenMP over points -----
      call gridpoint_pass()

      ! --- (3) forward transforms (t_transform) ------------------------------
      t0 = omp_get_wtime()
      do kk = 1, nlev
         ! (Vt=Fv, Vp=Fu) -> (spheroidal=ddiv, toroidal=dvor)
         call spat_to_SHsphtor(sht, Fv(:,kk), Fu(:,kk), ddiv(:,kk), dvor(:,kk))
         call spat_to_SH(sht, FT(:,kk), dtmp(:,kk))
         call spat_to_SH(sht, Fq(:,kk), dqhu(:,kk))
      end do
      call spat_to_SH(sht, Flnps, dlnps)
      t_transform = t_transform + (omp_get_wtime() - t0)

      ! --- (4) semi-implicit spectral solve (t_solve), OpenMP over modes ------
      call solve_pass()
   end subroutine do_step

   !---------------------------------------------------------------------------
   ! Gridpoint nonlinear tendencies.  ~30 genuine multiply-adds per point-level
   ! of primitive-equation nonlinear terms: absolute-vorticity flux, kinetic
   ! energy, thermodynamic and moisture advection, and div-compression terms.
   ! Parallelized over points (i owns all its levels) so Flnps has no race.
   !---------------------------------------------------------------------------
   subroutine gridpoint_pass()
      integer  :: i, kk
      real(dp) :: eta, ke, uu, vv, tt, qq, dvg, uT, vT, uq, vq

      t0 = omp_get_wtime()
      !$omp parallel do default(shared) private(i,kk,eta,ke,uu,vv,tt,qq,dvg,uT,vT,uq,vq) schedule(static)
      do i = 1, nspat
         Flnps(i) = 0.0_dp
         do kk = 1, nlev
            uu  = ug(i,kk);  vv = vg(i,kk)
            tt  = tg(i,kk);  qq = qg(i,kk)
            dvg = divg(i,kk)
            eta = vorg(i,kk) + fcor(i)                 ! absolute vorticity
            ke  = 0.5_dp * (uu*uu + vv*vv)             ! kinetic energy
            uT  = uu*tt;  vT = vv*tt                    ! heat advection products
            uq  = uu*qq;  vq = vv*qq                    ! moisture advection products
            ! momentum tendency (abs-vorticity flux + KE gradient stand-in)
            Fu(i,kk) =  eta*vv + ke
            Fv(i,kk) = -eta*uu - ke
            ! thermodynamic: -(u.grad)T stand-in - kappa*T*div compression
            FT(i,kk) = -(uT + vT) - 0.2856_dp * tt * dvg
            ! moisture: -(u.grad)q stand-in + q*div compression
            Fq(i,kk) = -(uq + vq) + qq * dvg
            ! log-surface-pressure tendency accrues -div over the column
            Flnps(i) = Flnps(i) - dvg
         end do
      end do
      !$omp end parallel do
      t_gridpoint = t_gridpoint + (omp_get_wtime() - t0)
   end subroutine gridpoint_pass

   !---------------------------------------------------------------------------
   ! Semi-implicit solve stand-in.  OpenMP over spectral modes; each mode does a
   ! size-nlev complex tridiagonal (Thomas) solve on the divergence tendency
   ! (vertical gravity-wave structure) then an l-dependent hyperdiffusion
   ! multiply on all four tendency fields.  A checksum is accumulated to prevent
   ! the compiler eliminating the (otherwise unused) results.
   !---------------------------------------------------------------------------
   subroutine solve_pass()
      integer  :: lm
      real(dp) :: partial

      t0 = omp_get_wtime()
      partial = 0.0_dp
      !$omp parallel do default(shared) private(lm) schedule(static) reduction(+:partial)
      do lm = 1, nlm
         call solve_mode(nlev, nlm, lm, damp(lm), asub, bdiag, csup, &
                         dvor, ddiv, dtmp, dqhu)
         partial = partial + real(ddiv(lm,nlev)) + real(dtmp(lm,1))
      end do
      !$omp end parallel do
      checksum = checksum + partial
      t_solve = t_solve + (omp_get_wtime() - t0)
   end subroutine solve_pass

   !---------------------------------------------------------------------------
   ! Per-mode vertical solve + hyperdiffusion.  Automatic (stack) arrays give
   ! each OpenMP thread its own scratch — hence OMP_STACKSIZE in the runme.
   !---------------------------------------------------------------------------
   subroutine solve_mode(nl, nm, lm, d, a, b, c, xvor, xdiv, xtmp, xqhu)
      integer,      intent(in)    :: nl, nm, lm
      real(dp),     intent(in)    :: d, a(nl), b(nl), c(nl)
      complex(dpc), intent(inout) :: xvor(nm,nl), xdiv(nm,nl), xtmp(nm,nl), xqhu(nm,nl)
      complex(dpc) :: rhs(nl), x(nl), dp1(nl)
      real(dp)     :: cp(nl), m
      integer      :: kk

      ! gather the divergence-tendency column across levels
      do kk = 1, nl
         rhs(kk) = xdiv(lm,kk)
      end do
      ! Thomas algorithm (real tridiagonal, complex rhs)
      cp(1)  = c(1) / b(1)
      dp1(1) = rhs(1) / b(1)
      do kk = 2, nl
         m      = b(kk) - a(kk)*cp(kk-1)
         cp(kk) = c(kk) / m
         dp1(kk)= (rhs(kk) - a(kk)*dp1(kk-1)) / m
      end do
      x(nl) = dp1(nl)
      do kk = nl-1, 1, -1
         x(kk) = dp1(kk) - cp(kk)*x(kk+1)
      end do
      ! write solved divergence back, and apply hyperdiffusion to all fields
      do kk = 1, nl
         xdiv(lm,kk) = x(kk)      * d
         xvor(lm,kk) = xvor(lm,kk)* d
         xtmp(lm,kk) = xtmp(lm,kk)* d
         xqhu(lm,kk) = xqhu(lm,kk)* d
      end do
   end subroutine solve_mode

   !---------------------------------------------------------------------------
   ! Degree l of each spectral coefficient, in SHTns ORTHONORMAL m-major order:
   ! m = 0..mmax (step mres=1), and for each m, l = m..lmax.
   !---------------------------------------------------------------------------
   subroutine build_lval(lmx, mmx, nm, lv)
      integer, intent(in)  :: lmx, mmx, nm
      integer, intent(out) :: lv(nm)
      integer :: mm, ll, cnt
      cnt = 0
      do mm = 0, mmx
         do ll = mm, lmx
            cnt = cnt + 1
            if (cnt <= nm) lv(cnt) = ll
         end do
      end do
      if (cnt /= nm) then
         write(*,*) "WARNING: reconstructed nlm=", cnt, " != SHTns nlm=", nm
      end if
   end subroutine build_lval

   !---------------------------------------------------------------------------
   ! Non-trivial smooth initial spectral fields: a handful of O(1) low modes
   ! decaying with l plus a small deterministic perturbation, so transforms and
   ! arithmetic operate on real (non-zero) data.
   !---------------------------------------------------------------------------
   subroutine init_spectral(nm, nl, lmx, lv, w_vor, w_div, w_tmp, w_qhu, w_lnps)
      integer,      intent(in)  :: nm, nl, lmx
      integer,      intent(in)  :: lv(nm)
      complex(dpc), intent(out) :: w_vor(nm,nl), w_div(nm,nl), w_tmp(nm,nl), w_qhu(nm,nl)
      complex(dpc), intent(out) :: w_lnps(nm)
      integer  :: lm, kk
      real(dp) :: amp, ph, lr
      do kk = 1, nl
         do lm = 1, nm
            lr  = real(lv(lm)+1, dp)
            amp = 1.0_dp / lr                          ! red spectrum
            ph  = 0.1_dp * real(lm, dp) + 0.05_dp * real(kk, dp)
            w_vor(lm,kk) = cmplx(amp*cos(ph),        0.6_dp*amp*sin(ph),      dpc)
            w_div(lm,kk) = cmplx(0.5_dp*amp*sin(ph), 0.4_dp*amp*cos(2*ph),    dpc)
            w_tmp(lm,kk) = cmplx(amp*cos(0.5_dp*ph), 0.3_dp*amp*sin(0.5_dp*ph),dpc)
            w_qhu(lm,kk) = cmplx(0.3_dp*amp*sin(ph), 0.2_dp*amp*cos(ph),      dpc)
         end do
      end do
      do lm = 1, nm
         lr = real(lv(lm)+1, dp)
         w_lnps(lm) = cmplx((1.0_dp/lr)*cos(0.2_dp*real(lm,dp)), &
                            0.3_dp*(1.0_dp/lr)*sin(0.2_dp*real(lm,dp)), dpc)
      end do
      ! keep the l=0,m=0 (global mean) modes real / well-defined
      w_vor(1,:) = cmplx(0.0_dp,0.0_dp,dpc)
      w_div(1,:) = cmplx(0.0_dp,0.0_dp,dpc)
   end subroutine init_spectral

   !---------------------------------------------------------------------------
   ! Emit the results block to unit `iu` (stdout or results.txt).
   !---------------------------------------------------------------------------
   subroutine report(iu)
      integer, intent(in) :: iu
      write(iu,'(a)') "=============================================================="
      write(iu,'(a)') " aeros M0a scaling benchmark (Option A) — results"
      write(iu,'(a)') "=============================================================="
      write(iu,'(a)') " NOTE: double precision (SHTns lib is double); Float32 core"
      write(iu,'(a)') "       would be ~2x faster -> these are a conservative upper"
      write(iu,'(a)') "       bound. Full Gaussian grid (no reduced/octahedral in"
      write(iu,'(a)') "       SHTns basic API): gridpoint cost is also an upper"
      write(iu,'(a)') "       bound (~1/3 fewer points on a reduced grid). Physics"
      write(iu,'(a)') "       excluded (~1% of budget, design §3.6)."
      write(iu,'(a)') "--------------------------------------------------------------"
      write(iu,'(a,a)')      " label                = ", trim(label)
      write(iu,'(a,i0)')     " lmax                 = ", lmax
      write(iu,'(a,i0)')     " nlev                 = ", nlev
      write(iu,'(a,i0,a,i0)')" grid (SHTns) nlat    = ", nlat_r, "   nphi = ", nphi_r
      write(iu,'(a,i0)')     " nlm (SHTns)          = ", nlm
      write(iu,'(a,i0)')     " nspat (SHTns)        = ", nspat
      write(iu,'(a,i0)')     " nthreads             = ", threads_used
      write(iu,'(a,a)')      " precision            = ", trim(precision_str)
      write(iu,'(a,f10.2)')  " dt (s)               = ", dt
      write(iu,'(a,i0)')     " nstep                = ", nstep
      write(iu,'(a,es14.6)') " t_step_s             = ", t_step_s
      write(iu,'(a,f8.4)')   " t_transform_frac     = ", fr_tr
      write(iu,'(a,f8.4)')   " t_gridpoint_frac     = ", fr_gp
      write(iu,'(a,f8.4)')   " t_solve_frac         = ", fr_sv
      write(iu,'(a,es14.6)') " steps_per_year       = ", steps_per_year
      write(iu,'(a,es14.6)') " wallclock_s_per_year = ", wall_per_year
      write(iu,'(a,es14.6)') " core_s_per_year      = ", core_per_year
      write(iu,'(a,es14.6)') " model_years_per_day  = ", sypd
      write(iu,'(a,es14.6)') " sypd                 = ", sypd
      write(iu,'(a,es14.6)') " (checksum)           = ", checksum
      write(iu,'(a,es14.6)') " (wallclock_timed_s)  = ", t_wall
      write(iu,'(a)') "--------------------------------------------------------------"
      write(iu,'(a)') "# machine-parseable key = value"
      write(iu,'(a,a)')      "label = ", trim(label)
      write(iu,'(a,i0)')     "lmax = ", lmax
      write(iu,'(a,i0)')     "nlev = ", nlev
      write(iu,'(a,i0)')     "nlat = ", nlat_r
      write(iu,'(a,i0)')     "nphi = ", nphi_r
      write(iu,'(a,i0)')     "nlm = ", nlm
      write(iu,'(a,i0)')     "nthreads = ", threads_used
      write(iu,'(a,a)')      "precision = ", trim(precision_str)
      write(iu,'(a,f0.2)')   "dt = ", dt
      write(iu,'(a,i0)')     "nstep = ", nstep
      write(iu,'(a,es16.8)') "t_step_s = ", t_step_s
      write(iu,'(a,f0.6)')   "t_transform_frac = ", fr_tr
      write(iu,'(a,f0.6)')   "t_gridpoint_frac = ", fr_gp
      write(iu,'(a,f0.6)')   "t_solve_frac = ", fr_sv
      write(iu,'(a,es16.8)') "steps_per_year = ", steps_per_year
      write(iu,'(a,es16.8)') "wallclock_s_per_year = ", wall_per_year
      write(iu,'(a,es16.8)') "core_s_per_year = ", core_per_year
      write(iu,'(a,es16.8)') "model_years_per_day = ", sypd
      write(iu,'(a,es16.8)') "sypd = ", sypd
      write(iu,'(a)') "=============================================================="
   end subroutine report

end program aeros_bench
