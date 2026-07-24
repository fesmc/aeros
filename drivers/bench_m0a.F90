program bench_m0a
    ! M0a: cost sizing for the spectral core.
    !
    ! docs/design.md section 11 admits the central number in the plan is an
    ! extrapolation: "T42L19 throughput (~380 SYPD/core) is an extrapolation,
    ! not a measurement: SpeedyWeather's T31L8 2,300 SYPD scaled by
    ! (42/31)^3 x (19/8) ~ 6. M0a exists to replace this estimate with a
    ! number." The whole section 3.6 budget table -- and therefore the
    ! resolution bind, and therefore section 3.7's reason to exist -- rests on
    ! it. This program produces that number.
    !
    ! WHAT IS MEASURED
    !
    !   transforms  One timestep's worth of spherical-harmonic transforms,
    !               timed as a single threaded loop over (field, level) pairs
    !               -- the real parallel pool, not one transform scaled up, so
    !               threading efficiency at the actual pool size is included.
    !
    !   solve       One timestep's semi-implicit solve. This is a COST PROXY,
    !               not physics: the true vertical-structure matrix arrives
    !               with the dynamical core at M1. What is faithful is its
    !               SHAPE -- a pre-factorized nlev x nlev system per spectral
    !               coefficient, with the matrix depending on l but not m, so
    !               factorization is amortized and the per-step cost is a pair
    !               of triangular solves, O(nlev^2) per coefficient. That is
    !               the standard semi-implicit spectral arrangement and it is
    !               what sets the flop count.
    !
    ! WHAT IS ASSUMED, and the one number to argue with:
    !
    !   ntr_per_level = 8 scalar-equivalent transforms per level per timestep:
    !     synthesis  (u,v) [vector = 2] + T + q                 = 4
    !     analysis   (vor,div) tendency [vector = 2] + T' + q'   = 4
    !   A vector transform counts as two scalar-equivalents. Single-level work
    !   (ln ps and its gradient) is neglected -- it is O(1/nlev) of the total.
    !   The headline result scales LINEARLY in this number, so it is a namelist
    !   parameter and is printed with every result.
    !
    ! WHAT IS NOT MEASURED: column physics (radiation, condensation, surface),
    ! which section 3.6 puts at ~1% of budget, and the polar nests (~4%). This
    ! is the dynamical core only.

    use aeros_defs,     only : sp, dp, wp, wp_sh, io_unit_err
    use aeros_spectral
    use nml,            only : nml_read

    implicit none

#ifndef AEROS_MACHINE
#define AEROS_MACHINE "unknown"
#endif
#ifndef AEROS_COMPILER
#define AEROS_COMPILER "unknown"
#endif

    integer, parameter :: nmax = 8    ! max entries in a sweep list

    character(len=512) :: path_par, csv_out
    character(len=256) :: label

    integer  :: truncs(nmax), nlevs(nmax), threads(nmax)
    integer  :: ntr_per_level, nrep
    real(wp) :: dt

    integer :: it, il, ih

    path_par = "par/bench_m0a.nml"
    if (command_argument_count() >= 1) call get_command_argument(1, path_par)

    truncs  = 0
    nlevs   = 0
    threads = 0

    call nml_read(path_par, "bench", "label",         label)
    call nml_read(path_par, "bench", "truncs",        truncs)
    call nml_read(path_par, "bench", "nlevs",         nlevs)
    call nml_read(path_par, "bench", "threads",       threads)
    call nml_read(path_par, "bench", "ntr_per_level", ntr_per_level)
    call nml_read(path_par, "bench", "dt",            dt)
    call nml_read(path_par, "bench", "nrep",          nrep)
    call nml_read(path_par, "bench", "csv_out",       csv_out)

    ! Optional label override, so the same namelist can produce distinguishable
    ! rows for the sp and dp builds without editing a file between runs.
    if (command_argument_count() >= 2) call get_command_argument(2, label)

    write(*,*) ""
    write(*,"(a)")         " === aeros M0a: dynamical-core cost sizing ==="
    write(*,"(a,a)")       "   label            ", trim(label)
    write(*,"(a,a)")       "   machine          ", AEROS_MACHINE
    write(*,"(a,a)")       "   compiler         ", AEROS_COMPILER
    write(*,"(a,i0,a)")    "   precision        ", storage_size(1.0_wp)/8, " byte reals (grid side)"
    write(*,"(a,i0)")      "   transforms/level ", ntr_per_level
    write(*,"(a,f0.1,a)")  "   dt               ", dt, " s"
    write(*,"(a,i0)")      "   repeats          ", nrep
    write(*,*) ""

    call csv_open(csv_out)

    ! Loop order is trunc -> threads -> nlev, because an SHTns pool depends on
    ! (trunc, threads) but NOT on nlev, and building one is expensive: tuning
    ! costs 8.2 s per config at T85 (measured). Building it inside the nlev
    ! loop would re-tune identical configs and make setup dominate the run.
    do it = 1, nmax
        if (truncs(it) <= 0) cycle

        write(*,"(a)") ""
        write(*,"(a,i0)") " --- T", truncs(it)
        write(*,"(a)") "   nlev  threads    transf_ms     solve_ms      step_ms  solve%       SYPD    core-s/yr"

        do ih = 1, nmax
            if (threads(ih) <= 0) cycle
            call run_threadcase(truncs(it), nlevs, threads(ih), ntr_per_level, &
                                    nrep, dt, label, csv_out)
        end do
    end do

    write(*,*) ""
    write(*,"(a,a)") " bench_m0a:: results appended to ", trim(csv_out)

contains

    subroutine run_threadcase(trunc, nlevs, nthreads, ntr_per_level, nrep, dt, &
                                    label, csv_out)
        ! One SHTns pool, every nlev in the sweep measured against it.

        implicit none

        integer,  intent(in) :: trunc, nthreads, ntr_per_level, nrep
        integer,  intent(in) :: nlevs(:)
        real(wp), intent(in) :: dt
        character(len=*), intent(in) :: label, csv_out

        type(aeros_sht_pool_class) :: pool
        real(dp) :: t_tr, t_sv, t_step, sypd, core_s, steps_per_yr
        integer  :: nlm, nlon, nlat, il, nlev

        ! quick=.FALSE. -- let SHTns auto-tune. Measuring the untuned path
        ! would understate the model by ~20%, i.e. would measure the wrong
        ! thing. cache=.TRUE. so the pool tunes once instead of once per
        ! thread; without it a 10-config T85 pool costs 82 s of setup.
        call aeros_sht_pool_init(pool, trunc, nthreads=nthreads, &
                                    quick=.FALSE., cache=.TRUE.)

        nlm  = pool%sht(1)%nlm
        nlon = pool%sht(1)%nlon
        nlat = pool%sht(1)%nlat

        ! A "model year" of dynamics: 365.2422 d, the tropical year the orbital
        ! forcing will be on, not 365.
        steps_per_yr = 365.2422_dp*86400.0_dp/real(dt, dp)

        do il = 1, size(nlevs)
            if (nlevs(il) <= 0) cycle
            nlev = nlevs(il)

            t_tr = time_transforms(pool, nlev, ntr_per_level, nrep, nthreads)
            t_sv = time_solve(nlm, nlev, pool%sht(1)%l_of_lm, pool%sht(1)%lmax, &
                                nrep, nthreads)

            t_step = t_tr + t_sv
            sypd   = 86400.0_dp/(steps_per_yr*t_step)
            core_s = steps_per_yr*t_step*real(nthreads, dp)

            write(*,"(i7,i9,f13.4,f13.4,f13.4,f8.1,f11.1,f13.1)") &
                    nlev, nthreads, t_tr*1.0e3_dp, t_sv*1.0e3_dp, t_step*1.0e3_dp, &
                    100.0_dp*t_sv/t_step, sypd, core_s

            call csv_row(csv_out, label, trunc, nlon, nlat, nlm, nlev, nthreads, &
                            ntr_per_level, dt, t_tr, t_sv, t_step, sypd, core_s)
        end do

        call aeros_sht_pool_end(pool)

        return

    end subroutine run_threadcase

    real(dp) function time_transforms(pool, nlev, ntr_per_level, nrep, nthreads) result(tstep)
        ! One timestep's transforms, as ONE threaded loop over (field, level).
        !
        ! Each unit is a synthesis + an analysis, i.e. 2 scalar-equivalents, so
        ! the loop runs ntr_per_level/2 * nlev units. Fields are per-unit, not
        ! shared, because threads write them concurrently.

        implicit none

        type(aeros_sht_pool_class), intent(in) :: pool
        integer, intent(in) :: nlev, ntr_per_level, nrep, nthreads

        type(aeros_sht_class), pointer :: s
        real(wp),       allocatable :: field(:,:,:)
        complex(wp_sh), allocatable :: coeffs(:,:)
        real(dp) :: t0, t1
        integer  :: nunit, u, r, nlon, nlat, nlm

        nunit = (ntr_per_level/2)*nlev
        nlon  = pool%sht(1)%nlon
        nlat  = pool%sht(1)%nlat
        nlm   = pool%sht(1)%nlm

        allocate(field(nlon,nlat,nunit))
        allocate(coeffs(nlm,nunit))

        do u = 1, nunit
            call seed_spectrum(coeffs(:,u), pool%sht(1), u)
        end do
        field = 0.0_wp

        ! Warm up: first touch of every page, and SHTns' internal buffers.
        !$omp parallel do num_threads(nthreads) private(u,s) schedule(static)
        do u = 1, nunit
            s => aeros_sht_pool_get(pool)
            call aeros_sht_synthesis(s, coeffs(:,u), field(:,:,u))
        end do
        !$omp end parallel do

        t0 = wall_time()
        do r = 1, nrep
            !$omp parallel do num_threads(nthreads) private(u,s) schedule(static)
            do u = 1, nunit
                s => aeros_sht_pool_get(pool)
                call aeros_sht_synthesis(s, coeffs(:,u), field(:,:,u))
                call aeros_sht_analysis(s, field(:,:,u), coeffs(:,u))
            end do
            !$omp end parallel do
        end do
        t1 = wall_time()

        tstep = (t1 - t0)/real(nrep, dp)

        deallocate(field, coeffs)

        return

    end function time_transforms

    real(dp) function time_solve(nlm, nlev, l_of_lm, lmax, nrep, nthreads) result(tstep)
        ! One timestep's semi-implicit solve, as a cost proxy.
        !
        ! Pre-factorized nlev x nlev system per spectral coefficient, matrix
        ! indexed by l alone. The factorization is done once (as it would be in
        ! the model, where it changes only when the reference state or dt does)
        ! and is NOT timed; the per-step cost is the triangular solve pair.

        implicit none

        integer, intent(in) :: nlm, nlev, lmax, nrep, nthreads
        integer, intent(in) :: l_of_lm(:)

        real(dp),       allocatable :: lu(:,:,:)     ! (nlev,nlev,0:lmax)
        complex(wp_sh), allocatable :: rhs(:,:)      ! (nlev,nlm)
        real(dp) :: t0, t1
        integer  :: l, lm, k, r

        allocate(lu(nlev,nlev,0:lmax))
        allocate(rhs(nlev,nlm))

        do l = 0, lmax
            call build_matrix(lu(:,:,l), nlev, l)
            call lu_factor(lu(:,:,l), nlev)
        end do

        do lm = 1, nlm
            do k = 1, nlev
                rhs(k,lm) = cmplx(real(k,dp)/real(nlev,dp), real(lm,dp)/real(nlm,dp), wp_sh)
            end do
        end do

        t0 = wall_time()
        do r = 1, nrep
            !$omp parallel do num_threads(nthreads) private(lm,l) schedule(static)
            do lm = 1, nlm
                l = l_of_lm(lm)
                call lu_solve(lu(:,:,l), nlev, rhs(:,lm))
            end do
            !$omp end parallel do
        end do
        t1 = wall_time()

        tstep = (t1 - t0)/real(nrep, dp)

        deallocate(lu, rhs)

        return

    end function time_solve

    subroutine build_matrix(a, n, l)
        ! A diagonally dominant nlev x nlev matrix standing in for the
        ! semi-implicit vertical structure operator at degree l. Its NUMBERS
        ! are arbitrary; its size, density and l-dependence are not, and those
        ! are what the timing depends on.

        implicit none

        real(dp), intent(out) :: a(:,:)
        integer,  intent(in)  :: n, l

        integer :: i, j

        do j = 1, n
            do i = 1, n
                a(i,j) = 1.0_dp/real(1 + abs(i-j), dp)
            end do
            a(j,j) = a(j,j) + real(n, dp) + real(l*(l+1), dp)/1.0e3_dp
        end do

        return

    end subroutine build_matrix

    subroutine lu_factor(a, n)
        ! Doolittle LU, no pivoting -- the real operator is diagonally
        ! dominant, and pivoting would change the cost model, not just the
        ! numerics. Done once, never timed.

        implicit none

        real(dp), intent(inout) :: a(:,:)
        integer,  intent(in) :: n

        integer  :: i, j, k
        real(dp) :: s

        do k = 1, n
            do j = k, n
                s = a(k,j)
                do i = 1, k-1
                    s = s - a(k,i)*a(i,j)
                end do
                a(k,j) = s
            end do
            do i = k+1, n
                s = a(i,k)
                do j = 1, k-1
                    s = s - a(i,j)*a(j,k)
                end do
                a(i,k) = s/a(k,k)
            end do
        end do

        return

    end subroutine lu_factor

    subroutine lu_solve(a, n, b)
        ! Forward + back substitution against a pre-factorized matrix.
        ! 2 x n^2/2 complex-by-real operations: the per-step semi-implicit cost.

        implicit none

        real(dp),       intent(in)    :: a(:,:)
        integer,        intent(in)    :: n
        complex(wp_sh), intent(inout) :: b(:)

        integer        :: i, j
        complex(wp_sh) :: s

        do i = 2, n
            s = b(i)
            do j = 1, i-1
                s = s - a(i,j)*b(j)
            end do
            b(i) = s
        end do

        do i = n, 1, -1
            s = b(i)
            do j = i+1, n
                s = s - a(i,j)*b(j)
            end do
            b(i) = s/a(i,i)
        end do

        return

    end subroutine lu_solve

    subroutine seed_spectrum(coeffs, sht, iseed)

        implicit none

        complex(wp_sh), intent(out) :: coeffs(:)
        type(aeros_sht_class), intent(in) :: sht
        integer, intent(in) :: iseed

        integer  :: lm, l
        real(dp) :: amp

        do lm = 1, sht%nlm
            l = sht%l_of_lm(lm)
            amp = 1.0_dp/real((l+1)*(l+1), dp)
            coeffs(lm) = cmplx(amp*cos(real(l+iseed,dp)), amp*sin(real(l+2*iseed,dp)), wp_sh)
        end do

        return

    end subroutine seed_spectrum

    real(dp) function wall_time() result(t)
        ! system_clock rather than omp_get_wtime, so the benchmark still builds
        ! and runs under openmp=0.

        implicit none

        integer(kind=8) :: count, count_rate

        call system_clock(count, count_rate)
        t = real(count, dp)/real(count_rate, dp)

        return

    end function wall_time

    subroutine csv_open(csv_out)
        ! Append-only, with a header written just once. The laptop run and the
        ! eventual HPC run are meant to land in the same file and be compared
        ! directly, so nothing here ever truncates it.

        implicit none

        character(len=*), intent(in) :: csv_out

        logical :: exists
        integer :: iou

        inquire(file=csv_out, exist=exists)
        if (exists) return

        open(newunit=iou, file=csv_out, status="new", action="write")
        write(iou,"(a)") "label,machine,compiler,wp_bytes,trunc,nlon,nlat,nlm,nlev,"// &
                            "threads,ntr_per_level,dt,t_transf_s,t_solve_s,t_step_s,sypd,core_s_per_yr"
        close(iou)

        return

    end subroutine csv_open

    subroutine csv_row(csv_out, label, trunc, nlon, nlat, nlm, nlev, nthreads, &
                            ntr_per_level, dt, t_tr, t_sv, t_step, sypd, core_s)

        implicit none

        character(len=*), intent(in) :: csv_out, label
        integer,  intent(in) :: trunc, nlon, nlat, nlm, nlev, nthreads, ntr_per_level
        real(wp), intent(in) :: dt
        real(dp), intent(in) :: t_tr, t_sv, t_step, sypd, core_s

        integer :: iou

        open(newunit=iou, file=csv_out, status="old", position="append", action="write")
        write(iou,"(a,a,a,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,f0.1,5(a,es14.6))") &
                trim(label), ",", AEROS_MACHINE, ",", AEROS_COMPILER, ",", &
                storage_size(1.0_wp)/8, ",", trunc, ",", nlon, ",", nlat, ",", nlm, ",", &
                nlev, ",", nthreads, ",", ntr_per_level, ",", dt, &
                ",", t_tr, ",", t_sv, ",", t_step, ",", sypd, ",", core_s
        close(iou)

        return

    end subroutine csv_row

end program bench_m0a
