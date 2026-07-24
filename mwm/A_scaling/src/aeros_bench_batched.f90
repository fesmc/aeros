!==============================================================================
! aeros_bench_batched  —  M0a cost / OpenMP-scaling harness, VARIANT 2
!------------------------------------------------------------------------------
! Same benchmark as src/aeros_bench.f90 (variant 1, "spectral_shtns"), with ONE
! difference: the spherical-harmonic transforms are done by a CUSTOM batched
! Legendre transform instead of SHTns-one-2D-field-at-a-time.
!
!   * Legendre step: for each zonal wavenumber m, ALL nlev levels are batched
!     into ONE real GEMM (MKL sequential dgemm).  Parallelism is OpenMP over m
!     (schedule(dynamic)); BLAS itself runs single-threaded inside each GEMM.
!   * Longitude step: FFTW real<->complex, one length-nphi FFT per (lat,level),
!     OpenMP-parallel over (lat,level) with per-thread scratch and a shared
!     serial FFTW plan (fftw_execute_dft_* on distinct data is thread-safe).
!
! This variant exists to answer the decisive M0a question: does the batched /
! threaded-BLAS transform scale past ~4 threads where SHTns-per-field did not?
!
! Everything OTHER than the transform is byte-identical to variant 1: the
! namelist interface, the thread handling (OMP_NUM_THREADS via
! omp_get_max_threads), the gridpoint nonlinear pass, the semi-implicit
! spectral-solve stand-in, and the results.txt schema.  The comparison is only
! clean if the non-transform work is identical, so it is copied verbatim.
!
! Grid: full Gaussian grid (Gauss-Legendre latitudes computed here), triangular
! truncation T(lmax), mres = 1.  Double precision (matches variant 1; Float32
! would be ~2x faster -> these are a conservative upper bound).
!==============================================================================
program aeros_bench_batched
   use, intrinsic :: iso_c_binding
   use omp_lib
   implicit none

   include "fftw3.f03"      ! FFTW C-bound interfaces + FFTW_ESTIMATE etc.

   integer,  parameter :: dp  = c_double            ! spatial reals
   integer,  parameter :: dpc = c_double_complex    ! spectral coeffs

   ! ---- namelist inputs (IDENTICAL to variant 1) ----------------------------
   integer            :: lmax, nlat, nphi, nlev, nstep, nwarm
   real(dp)           :: dt
   character(len=64)  :: label
   namelist /bench/ lmax, nlat, nphi, nlev, dt, nstep, nwarm, label

   ! ---- grid / spectral metadata (computed here, not from a library) --------
   integer :: mmax, nlm, nspat, nhalf
   integer :: nthreads, threads_used
   logical :: do_validate

   ! ---- custom-transform precomputed tables ---------------------------------
   real(dp), allocatable :: mu(:), gw(:)         ! (nlat) Gauss nodes, weights
   real(dp), allocatable :: plm(:,:)             ! (nlat,nlm) normalized P_l^m
   real(dp), allocatable :: vfac(:)              ! (nlm) cheap vector-transform factor
   integer,  allocatable :: moff(:), kcnt(:)     ! (0:lmax) block offset / degree count
   integer,  allocatable :: lval(:), mval(:)     ! (nlm) degree l, order m of each mode
   complex(dpc), allocatable :: Fc(:,:,:)        ! (0:nphi/2,nlat,nlev) Fourier scratch
   type(c_ptr) :: r2c_plan, c2r_plan

   ! ---- prognostic spectral fields (held fixed across timed steps) ----------
   complex(dpc), allocatable :: vor(:,:), div(:,:), tmp(:,:), qhu(:,:)  ! (nlm,nlev)
   complex(dpc), allocatable :: lnps(:)                                 ! (nlm)
   ! ---- tendency spectral fields (recomputed every step) --------------------
   complex(dpc), allocatable :: dvor(:,:), ddiv(:,:), dtmp(:,:), dqhu(:,:)
   complex(dpc), allocatable :: dlnps(:)
   ! ---- vector-transform spectral scratch (u,v <-> vor,div stand-ins) -------
   complex(dpc), allocatable :: uspec(:,:), vspec(:,:)

   ! ---- spatial (gridpoint) work arrays, all (nspat,nlev) unless noted -------
   real(dp), allocatable :: ug(:,:), vg(:,:)          ! wind components
   real(dp), allocatable :: vorg(:,:), divg(:,:)      ! rel. vorticity, divergence
   real(dp), allocatable :: tg(:,:),  qg(:,:)         ! temperature, humidity
   real(dp), allocatable :: Fu(:,:),  Fv(:,:)         ! momentum tendency
   real(dp), allocatable :: FT(:,:),  Fq(:,:)         ! thermodynamic / moisture tendency
   real(dp), allocatable :: lnpsg(:), Flnps(:)        ! (nspat) log-surface-pressure
   real(dp), allocatable :: fcor(:)                   ! (nspat) Coriolis (representative)

   ! ---- implicit vertical structure + hyperdiffusion ------------------------
   real(dp), allocatable :: asub(:), bdiag(:), csup(:)   ! (nlev) tridiag coeffs
   real(dp), allocatable :: damp(:)                      ! (nlm) l-dependent damping

   ! ---- timing --------------------------------------------------------------
   real(dp) :: t_transform, t_gridpoint, t_solve
   real(dp) :: t0, t_wall0, t_wall
   real(dp) :: t_step_s, phase_sum
   real(dp) :: fr_tr, fr_gp, fr_sv
   real(dp) :: checksum

   ! ---- reporting -----------------------------------------------------------
   real(dp) :: steps_per_year, wall_per_year, core_per_year, sypd
   character(len=512) :: nmlpath, arg2
   character(len=8)   :: precision_str = "double"
   character(len=16)  :: variant = "batched"
   integer :: istep, ios, u

   ! ==========================================================================
   ! 1. Read namelist (same &bench group / same argv(1) as variant 1).
   !    Optional argv(2) == "--validate" runs the round-trip check and exits.
   ! ==========================================================================
   lmax = 31; nlat = 48; nphi = 96; nlev = 16
   dt = 1800.0_dp; nstep = 100; nwarm = 10; label = "unnamed"

   call get_command_argument(1, nmlpath, status=ios)
   if (ios /= 0 .or. len_trim(nmlpath) == 0) then
      write(*,*) "ERROR: usage: aeros_bench_batched.x <bench.nml> [--validate]"
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

   do_validate = .false.
   call get_command_argument(2, arg2, status=ios)
   if (ios == 0 .and. trim(arg2) == "--validate") do_validate = .true.

   mmax  = lmax                              ! triangular truncation, mres = 1
   nlm   = (lmax+1)*(lmax+2)/2               ! spectral coeffs, m-major l=m..lmax
   nspat = nlat*nphi                         ! full Gaussian grid points
   nhalf = nphi/2 + 1                        ! real-FFT half-spectrum length

   ! ==========================================================================
   ! 2. Threads: from OMP_NUM_THREADS (via omp_get_max_threads), NOT namelist.
   !    This is what our !$omp parallel do over m and over (lat,level) uses.
   ! ==========================================================================
   nthreads     = omp_get_max_threads()
   threads_used = nthreads

   ! ==========================================================================
   ! 3. Allocate
   ! ==========================================================================
   allocate(mu(nlat), gw(nlat), plm(nlat,nlm), vfac(nlm))
   allocate(moff(0:lmax), kcnt(0:lmax), lval(nlm), mval(nlm))
   allocate(Fc(0:nphi/2, nlat, nlev))
   allocate(vor(nlm,nlev), div(nlm,nlev), tmp(nlm,nlev), qhu(nlm,nlev), lnps(nlm))
   allocate(dvor(nlm,nlev), ddiv(nlm,nlev), dtmp(nlm,nlev), dqhu(nlm,nlev), dlnps(nlm))
   allocate(uspec(nlm,nlev), vspec(nlm,nlev))
   allocate(ug(nspat,nlev), vg(nspat,nlev), vorg(nspat,nlev), divg(nspat,nlev))
   allocate(tg(nspat,nlev), qg(nspat,nlev))
   allocate(Fu(nspat,nlev), Fv(nspat,nlev), FT(nspat,nlev), Fq(nspat,nlev))
   allocate(lnpsg(nspat), Flnps(nspat), fcor(nspat))
   allocate(asub(nlev), bdiag(nlev), csup(nlev), damp(nlm))

   ! ==========================================================================
   ! 4. Build the transform tables + FFTW plans
   ! ==========================================================================
   call build_index(lmax, nlm, moff, kcnt, lval, mval)
   call gauleg(nlat, mu, gw)                 ! Gauss-Legendre nodes + weights
   call build_plm(lmax, nlat, nlm, mu, moff, plm)
   Fc = (0.0_dp, 0.0_dp)                      ! high-m rows stay zero across inverse calls

   ! serial FFTW plans, shared across threads; UNALIGNED so new-array execute is
   ! safe with any per-thread scratch alignment.
   block
      real(dp)     :: planr(nphi)
      complex(dpc) :: planc(0:nphi/2)
      r2c_plan = fftw_plan_dft_r2c_1d(nphi, planr, planc, ior(FFTW_ESTIMATE, FFTW_UNALIGNED))
      c2r_plan = fftw_plan_dft_c2r_1d(nphi, planc, planr, ior(FFTW_ESTIMATE, FFTW_UNALIGNED))
   end block

   ! cheap (l,m) factor for the vector-transform stand-in (u,v <- vor,div)
   block
      integer :: lm
      do lm = 1, nlm
         vfac(lm) = 1.0_dp / real(lval(lm)+1, dp)
      end do
   end block

   ! ==========================================================================
   ! 5. Validation mode: random band-limited spectral field, inverse then
   !    forward, report max round-trip error.  A fast WRONG transform is useless.
   ! ==========================================================================
   if (do_validate) then
      call run_validation()
      call fftw_destroy_plan(r2c_plan)
      call fftw_destroy_plan(c2r_plan)
      stop 0
   end if

   ! ==========================================================================
   ! 6. Initialize benchmark fields (identical to variant 1)
   ! ==========================================================================
   call init_spectral(nlm, nlev, lval, vor, div, tmp, qhu, lnps)

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

   asub  = -1.0_dp
   csup  = -1.0_dp
   bdiag =  2.5_dp

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
   ! 7. Warmup (untimed) then timed loop
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
   ! 8. Derived metrics (identical to variant 1)
   ! ==========================================================================
   phase_sum = t_transform + t_gridpoint + t_solve
   t_step_s  = phase_sum / real(nstep, dp)
   fr_tr = t_transform / phase_sum
   fr_gp = t_gridpoint / phase_sum
   fr_sv = t_solve     / phase_sum

   steps_per_year = 365.0_dp * 86400.0_dp / dt
   wall_per_year  = t_step_s * steps_per_year
   core_per_year  = wall_per_year * real(threads_used, dp)
   sypd           = 86400.0_dp / wall_per_year

   ! ==========================================================================
   ! 9. Report
   ! ==========================================================================
   call report(6)
   open(newunit=u, file="results.txt", status="replace", action="write")
   call report(u)
   close(u)

   call fftw_destroy_plan(r2c_plan)
   call fftw_destroy_plan(c2r_plan)

contains

   !---------------------------------------------------------------------------
   ! One timestep: inverse transforms -> gridpoint nonlinear -> forward
   ! transforms -> implicit solve.  Same field workload / transform COUNT as
   ! variant 1 (see notes below); only the transform IMPLEMENTATION differs.
   !---------------------------------------------------------------------------
   subroutine do_step()

      ! --- (1) inverse transforms (t_transform) ------------------------------
      ! Variant-1 count reproduced: per level a vector inverse (vor,div -> u,v)
      ! = 2 scalar Legendre transforms, plus scalar inverses of vor,div,T,q, and
      ! lnps once.  Total 6 batched(nlev) scalar transforms + 1 batched(1).
      t0 = omp_get_wtime()
      ! vector stand-in: build u,v spectral from vor,div by cheap (l,m) multiply
      call vfac_mul(vor, uspec)
      call vfac_mul(div, vspec)
      call sh_to_grid(uspec, ug,   nlev)
      call sh_to_grid(vspec, vg,   nlev)
      call sh_to_grid(vor,   vorg, nlev)
      call sh_to_grid(div,   divg, nlev)
      call sh_to_grid(tmp,   tg,   nlev)
      call sh_to_grid(qhu,   qg,   nlev)
      call sh_to_grid1(lnps, lnpsg)
      t_transform = t_transform + (omp_get_wtime() - t0)

      ! --- (2) gridpoint nonlinear pass (t_gridpoint) ------------------------
      call gridpoint_pass()

      ! --- (3) forward transforms (t_transform) ------------------------------
      ! Variant-1 count reproduced: per level a vector forward (Fv,Fu ->
      ! div,vor tendencies) = 2 scalar transforms, plus scalar forwards of
      ! FT,Fq, and lnps once.  Total 4 batched(nlev) scalar transforms + 1.
      t0 = omp_get_wtime()
      call grid_to_sh(Fu, uspec, nlev)
      call grid_to_sh(Fv, vspec, nlev)
      call vfac_mul(uspec, dvor)
      call vfac_mul(vspec, ddiv)
      call grid_to_sh(FT, dtmp, nlev)
      call grid_to_sh(Fq, dqhu, nlev)
      call grid_to_sh1(Flnps, dlnps)
      t_transform = t_transform + (omp_get_wtime() - t0)

      ! --- (4) semi-implicit spectral solve (t_solve) ------------------------
      call solve_pass()
   end subroutine do_step

   !---------------------------------------------------------------------------
   ! Cheap (l,m) spectral multiply — the vector-transform stand-in.  O(nlm*nlev)
   ! multiplies (counted inside t_transform, as in variant 1's vector call).
   !---------------------------------------------------------------------------
   subroutine vfac_mul(S, R)
      complex(dpc), intent(in)  :: S(nlm,nlev)
      complex(dpc), intent(out) :: R(nlm,nlev)
      integer :: lm, lev
      do lev = 1, nlev
         do lm = 1, nlm
            R(lm,lev) = vfac(lm) * S(lm,lev)
         end do
      end do
   end subroutine vfac_mul

   !===========================================================================
   ! CUSTOM BATCHED INVERSE TRANSFORM  (spectral -> grid)
   !   Step A: per m, ONE dgemm batching all nl levels (real & imag packed as
   !           2*nl columns): Gpack(nlat,2nl) = Pm(nlat,Km) . Spack(Km,2nl).
   !           OpenMP parallel over m, sequential BLAS inside.
   !   Step B: per (lat,level) inverse real-FFT along longitude (OpenMP over the
   !           lat*level pairs, per-thread scratch, shared serial plan).
   !===========================================================================
   subroutine sh_to_grid(S, g, nl)
      integer,      intent(in)  :: nl
      complex(dpc), intent(in)  :: S(nlm,nl)
      real(dp),     intent(out) :: g(nspat,nl)
      integer      :: m, k, j, lev, kk, off, p
      real(dp)     :: Spack(lmax+1, 2*nl)
      real(dp)     :: Gpack(nlat,   2*nl)
      real(dp)     :: rscr(nphi)
      complex(dpc) :: cscr(0:nphi/2)

      ! keep the unused high-m rows zero (a previous forward filled them)
      if (lmax < nphi/2) Fc(lmax+1:nphi/2, :, :) = (0.0_dp, 0.0_dp)

      !$omp parallel do default(shared) private(m,k,j,lev,kk,off,Spack,Gpack) schedule(dynamic)
      do m = 0, lmax
         kk  = kcnt(m)
         off = moff(m)
         do lev = 1, nl
            do k = 1, kk
               Spack(k,     lev) = real (S(off+k,lev))
               Spack(k, nl+ lev) = aimag(S(off+k,lev))
            end do
         end do
         call dgemm('N','N', nlat, 2*nl, kk, 1.0_dp, &
                    plm(1,off+1), nlat, Spack, lmax+1, 0.0_dp, Gpack, nlat)
         do lev = 1, nl
            do j = 1, nlat
               Fc(m,j,lev) = cmplx(Gpack(j,lev), Gpack(j,nl+lev), dpc)
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(j,lev,p,cscr,rscr) collapse(2) schedule(static)
      do lev = 1, nl
         do j = 1, nlat
            do p = 0, nphi/2
               cscr(p) = Fc(p,j,lev)
            end do
            call fftw_execute_dft_c2r(c2r_plan, cscr, rscr)
            do p = 1, nphi
               g(j + (p-1)*nlat, lev) = rscr(p)
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine sh_to_grid

   !===========================================================================
   ! CUSTOM BATCHED FORWARD TRANSFORM  (grid -> spectral)
   !   Step A: per (lat,level) forward real-FFT along longitude; apply Gauss
   !           weight w_j and the 1/nphi FFT normalization.
   !   Step B: per m, ONE dgemm (transpose): Spack(Km,2nl) = Pm^T(Km,nlat) .
   !           Ghat(nlat,2nl).  OpenMP parallel over m, sequential BLAS inside.
   !===========================================================================
   subroutine grid_to_sh(g, S, nl)
      integer,      intent(in)  :: nl
      real(dp),     intent(in)  :: g(nspat,nl)
      complex(dpc), intent(out) :: S(nlm,nl)
      integer      :: m, k, j, lev, kk, off, p
      real(dp)     :: Spack(lmax+1, 2*nl)
      real(dp)     :: Gpack(nlat,   2*nl)
      real(dp)     :: rscr(nphi), wj
      complex(dpc) :: cscr(0:nphi/2)

      !$omp parallel do default(shared) private(j,lev,p,cscr,rscr,wj) collapse(2) schedule(static)
      do lev = 1, nl
         do j = 1, nlat
            do p = 1, nphi
               rscr(p) = g(j + (p-1)*nlat, lev)
            end do
            call fftw_execute_dft_r2c(r2c_plan, rscr, cscr)
            wj = gw(j) / real(nphi, dp)
            do p = 0, nphi/2
               Fc(p,j,lev) = cscr(p) * wj
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(m,k,j,lev,kk,off,Spack,Gpack) schedule(dynamic)
      do m = 0, lmax
         kk  = kcnt(m)
         off = moff(m)
         do lev = 1, nl
            do j = 1, nlat
               Gpack(j,     lev) = real (Fc(m,j,lev))
               Gpack(j, nl+ lev) = aimag(Fc(m,j,lev))
            end do
         end do
         call dgemm('T','N', kk, 2*nl, nlat, 1.0_dp, &
                    plm(1,off+1), nlat, Gpack, nlat, 0.0_dp, Spack, lmax+1)
         do lev = 1, nl
            do k = 1, kk
               S(off+k,lev) = cmplx(Spack(k,lev), Spack(k,nl+lev), dpc)
            end do
         end do
      end do
      !$omp end parallel do
   end subroutine grid_to_sh

   ! single-level wrappers (lnps) — reuse the batched routines with nl=1
   subroutine sh_to_grid1(S, g)
      complex(dpc), intent(in)  :: S(nlm)
      real(dp),     intent(out) :: g(nspat)
      complex(dpc) :: S2(nlm,1)
      real(dp)     :: g2(nspat,1)
      S2(:,1) = S
      call sh_to_grid(S2, g2, 1)
      g = g2(:,1)
   end subroutine sh_to_grid1

   subroutine grid_to_sh1(g, S)
      real(dp),     intent(in)  :: g(nspat)
      complex(dpc), intent(out) :: S(nlm)
      real(dp)     :: g2(nspat,1)
      complex(dpc) :: S2(nlm,1)
      g2(:,1) = g
      call grid_to_sh(g2, S2, 1)
      S = S2(:,1)
   end subroutine grid_to_sh1

   !---------------------------------------------------------------------------
   ! Round-trip validation: for a random band-limited spectral field with real
   ! m=0 coefficients, forward(inverse(S)) must equal S to ~machine precision.
   !---------------------------------------------------------------------------
   subroutine run_validation()
      complex(dpc), allocatable :: Sref(:,:), Srt(:,:)
      real(dp),     allocatable :: grid(:,:)
      integer  :: lm, lev
      real(dp) :: amp, ph, err
      allocate(Sref(nlm,nlev), Srt(nlm,nlev), grid(nspat,nlev))
      do lev = 1, nlev
         do lm = 1, nlm
            amp = 1.0_dp / real(lval(lm)+1, dp)
            ph  = 0.7_dp*real(lm,dp) + 0.31_dp*real(lev,dp)
            if (mval(lm) == 0) then
               Sref(lm,lev) = cmplx(amp*cos(ph), 0.0_dp, dpc)   ! zonal mean is real
            else
               Sref(lm,lev) = cmplx(amp*cos(ph), amp*sin(1.3_dp*ph), dpc)
            end if
         end do
      end do
      call sh_to_grid(Sref, grid, nlev)
      call grid_to_sh(grid, Srt,  nlev)
      err = maxval(abs(Srt - Sref))
      write(*,'(a)')        "=============================================================="
      write(*,'(a)')        " aeros_bench_batched — transform round-trip validation"
      write(*,'(a)')        "=============================================================="
      write(*,'(a,a)')      " label              = ", trim(label)
      write(*,'(a,i0,a,i0,a,i0)') " lmax = ", lmax, "  nlat = ", nlat, "  nphi = ", nphi
      write(*,'(a,i0,a,i0)')" nlm  = ", nlm, "  nlev = ", nlev
      write(*,'(a,es14.6)') " max_roundtrip_abs_err = ", err
      if (err <= 1.0e-10_dp) then
         write(*,'(a)')     " VALIDATION: PASS (<= 1e-10)"
      else
         write(*,'(a)')     " VALIDATION: FAIL (> 1e-10) — transform is wrong"
      end if
      write(*,'(a)')        "=============================================================="
      deallocate(Sref, Srt, grid)
   end subroutine run_validation

   !===========================================================================
   ! The following are COPIED VERBATIM from variant 1 (src/aeros_bench.f90):
   ! gridpoint_pass, solve_pass, solve_mode, init_spectral.  Keeping them
   ! byte-identical is what makes the transform-only comparison clean.
   !===========================================================================

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
            Fu(i,kk) =  eta*vv + ke
            Fv(i,kk) = -eta*uu - ke
            FT(i,kk) = -(uT + vT) - 0.2856_dp * tt * dvg
            Fq(i,kk) = -(uq + vq) + qq * dvg
            Flnps(i) = Flnps(i) - dvg
         end do
      end do
      !$omp end parallel do
      t_gridpoint = t_gridpoint + (omp_get_wtime() - t0)
   end subroutine gridpoint_pass

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

   subroutine solve_mode(nl, nm, lm, d, a, b, c, xvor, xdiv, xtmp, xqhu)
      integer,      intent(in)    :: nl, nm, lm
      real(dp),     intent(in)    :: d, a(nl), b(nl), c(nl)
      complex(dpc), intent(inout) :: xvor(nm,nl), xdiv(nm,nl), xtmp(nm,nl), xqhu(nm,nl)
      complex(dpc) :: rhs(nl), x(nl), dp1(nl)
      real(dp)     :: cp(nl), m
      integer      :: kk

      do kk = 1, nl
         rhs(kk) = xdiv(lm,kk)
      end do
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
      do kk = 1, nl
         xdiv(lm,kk) = x(kk)      * d
         xvor(lm,kk) = xvor(lm,kk)* d
         xtmp(lm,kk) = xtmp(lm,kk)* d
         xqhu(lm,kk) = xqhu(lm,kk)* d
      end do
   end subroutine solve_mode

   subroutine init_spectral(nm, nl, lv, w_vor, w_div, w_tmp, w_qhu, w_lnps)
      integer,      intent(in)  :: nm, nl
      integer,      intent(in)  :: lv(nm)
      complex(dpc), intent(out) :: w_vor(nm,nl), w_div(nm,nl), w_tmp(nm,nl), w_qhu(nm,nl)
      complex(dpc), intent(out) :: w_lnps(nm)
      integer  :: lm, kk
      real(dp) :: amp, ph, lr
      do kk = 1, nl
         do lm = 1, nm
            lr  = real(lv(lm)+1, dp)
            amp = 1.0_dp / lr
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
      w_vor(1,:) = cmplx(0.0_dp,0.0_dp,dpc)
      w_div(1,:) = cmplx(0.0_dp,0.0_dp,dpc)
   end subroutine init_spectral

   !===========================================================================
   ! Custom-transform helpers
   !===========================================================================

   !---------------------------------------------------------------------------
   ! Spectral index tables: m-major layout, m=0..lmax and l=m..lmax within each
   ! m — the SAME ordering as SHTns ORTHONORMAL, so the copied solve/gridpoint
   ! code sees identically-shaped data.  moff(m) = 0-based start of block m,
   ! kcnt(m) = number of degrees in block m (= lmax-m+1).
   !---------------------------------------------------------------------------
   subroutine build_index(lmx, nm, off, cnt, lv, mv)
      integer, intent(in)  :: lmx, nm
      integer, intent(out) :: off(0:lmx), cnt(0:lmx), lv(nm), mv(nm)
      integer :: mm, ll, c
      c = 0
      do mm = 0, lmx
         off(mm) = c
         cnt(mm) = lmx - mm + 1
         do ll = mm, lmx
            c = c + 1
            lv(c) = ll
            mv(c) = mm
         end do
      end do
      if (c /= nm) write(*,*) "WARNING: index count", c, " != nlm ", nm
   end subroutine build_index

   !---------------------------------------------------------------------------
   ! Gauss-Legendre nodes (mu = cos theta) and weights on (-1,1), computed by
   ! Newton-Raphson on the Legendre polynomial P_n (standard algorithm).
   !---------------------------------------------------------------------------
   subroutine gauleg(n, x, w)
      integer,  intent(in)  :: n
      real(dp), intent(out) :: x(n), w(n)
      integer,  parameter :: maxit = 100
      real(dp), parameter :: pi = 3.141592653589793238462643_dp
      real(dp), parameter :: eps = 3.0e-15_dp
      integer  :: i, j, it, mroots
      real(dp) :: z, z1, p1, p2, p3, pp

      mroots = (n + 1) / 2
      do i = 1, mroots
         z = cos(pi * (real(i,dp) - 0.25_dp) / (real(n,dp) + 0.5_dp))
         do it = 1, maxit
            p1 = 1.0_dp
            p2 = 0.0_dp
            do j = 1, n
               p3 = p2
               p2 = p1
               p1 = ((2.0_dp*real(j,dp) - 1.0_dp)*z*p2 - (real(j,dp) - 1.0_dp)*p3) / real(j,dp)
            end do
            pp = real(n,dp) * (z*p1 - p2) / (z*z - 1.0_dp)   ! P_n'
            z1 = z
            z  = z1 - p1/pp                                   ! Newton step
            if (abs(z - z1) <= eps) exit
         end do
         x(i)         =  z
         x(n + 1 - i) = -z
         w(i)         = 2.0_dp / ((1.0_dp - z*z) * pp * pp)
         w(n + 1 - i) = w(i)
      end do
   end subroutine gauleg

   !---------------------------------------------------------------------------
   ! Fully-normalized associated Legendre functions P_l^m(mu_j), normalized so
   ! that  sum_j w_j P_l^m(mu_j)^2 = 1  (Gauss quadrature exact for nlat>=lmax+1).
   ! With this normalization + the 1/nphi FFT factor + Gauss weights, the
   ! forward transform is the exact inverse of the synthesis (round-trip = id).
   ! Stored m-major: plm(:, moff(m)+ (l-m+1)) = P_l^m.
   !---------------------------------------------------------------------------
   subroutine build_plm(lmx, nl, nm, muj, off, p)
      integer,  intent(in)  :: lmx, nl, nm
      real(dp), intent(in)  :: muj(nl)
      integer,  intent(in)  :: off(0:lmx)
      real(dp), intent(out) :: p(nl,nm)
      integer  :: mm, ll, j, col
      real(dp) :: sinj(nl), pmm(nl), a, b
      real(dp), parameter :: invsqrt2 = 0.7071067811865475244_dp

      do j = 1, nl
         sinj(j) = sqrt(max(0.0_dp, 1.0_dp - muj(j)*muj(j)))
      end do

      ! P_0^0 = 1/sqrt(2)
      pmm(:) = invsqrt2
      do mm = 0, lmx
         if (mm > 0) then
            ! sectoral recurrence: P_m^m = sqrt((2m+1)/(2m)) sin * P_{m-1}^{m-1}
            a = sqrt(real(2*mm+1, dp) / real(2*mm, dp))
            do j = 1, nl
               pmm(j) = a * sinj(j) * pmm(j)
            end do
         end if
         col = off(mm) + 1                 ! l = mm
         do j = 1, nl
            p(j, col) = pmm(j)
         end do
         if (mm < lmx) then
            ! l = mm+1:  a_{m+1}^m = sqrt(2m+3)
            a = sqrt(real(2*mm+3, dp))
            do j = 1, nl
               p(j, col+1) = a * muj(j) * pmm(j)
            end do
         end if
         do ll = mm+2, lmx
            a = sqrt( real((2*ll+1)*(2*ll-1), dp) / real((ll-mm)*(ll+mm), dp) )
            b = sqrt( real(2*ll+1, dp) * real((ll-1-mm)*(ll-1+mm), dp) &
                      / ( real(2*ll-3, dp) * real((ll-mm)*(ll+mm), dp) ) )
            col = off(mm) + (ll - mm + 1)
            do j = 1, nl
               p(j, col) = a*muj(j)*p(j, col-1) - b*p(j, col-2)
            end do
         end do
      end do
   end subroutine build_plm

   !---------------------------------------------------------------------------
   ! Emit results (same schema as variant 1) + one new key: variant = batched.
   !---------------------------------------------------------------------------
   subroutine report(iu)
      integer, intent(in) :: iu
      write(iu,'(a)') "=============================================================="
      write(iu,'(a)') " aeros M0a scaling benchmark — VARIANT 2 (batched DGEMM)"
      write(iu,'(a)') "=============================================================="
      write(iu,'(a)') " Transform: custom batched Legendre (all levels/m in one"
      write(iu,'(a)') "       dgemm, OpenMP over m, sequential MKL BLAS) + FFTW."
      write(iu,'(a)') " Non-transform work (gridpoint + solve) identical to"
      write(iu,'(a)') "       variant 1 (spectral_shtns). Double precision -> a"
      write(iu,'(a)') "       conservative upper bound (Float32 ~2x faster)."
      write(iu,'(a)') "--------------------------------------------------------------"
      write(iu,'(a,a)')      " label                = ", trim(label)
      write(iu,'(a,a)')      " variant              = ", trim(variant)
      write(iu,'(a,i0)')     " lmax                 = ", lmax
      write(iu,'(a,i0)')     " nlev                 = ", nlev
      write(iu,'(a,i0,a,i0)')" grid nlat            = ", nlat, "   nphi = ", nphi
      write(iu,'(a,i0)')     " nlm                  = ", nlm
      write(iu,'(a,i0)')     " nspat                = ", nspat
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
      write(iu,'(a,a)')      "variant = ", trim(variant)
      write(iu,'(a,i0)')     "lmax = ", lmax
      write(iu,'(a,i0)')     "nlev = ", nlev
      write(iu,'(a,i0)')     "nlat = ", nlat
      write(iu,'(a,i0)')     "nphi = ", nphi
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

end program aeros_bench_batched
