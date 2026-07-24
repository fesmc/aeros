program bench_sht
    ! Where does the parallelism in the spectral transforms actually come from?
    !
    ! docs/design.md section 4.3 assumes OpenMP scaling to ~32 threads and
    ! names it the decisive engineering risk (section 9 risk 2), because the
    ! effective thread count is what picks T31+correction over bare T42. But
    ! "OpenMP scaling" is ambiguous for a spectral core, and the two readings
    ! have completely different prospects at aeros' truncations:
    !
    !   MODE A -- SHTns threads ONE transform over many cores, and aeros loops
    !             over levels serially. This is the mode a single call to
    !             shtns_use_threads(N) buys. At T31 a single scalar transform
    !             is ~500 coefficients over a 96x48 grid: a few hundred
    !             microseconds of work. Threading THAT is a synchronization
    !             problem, not a throughput one.
    !
    !   MODE B -- SHTns runs single-threaded, and aeros threads its own loop
    !             over levels, one SHTns config per thread. The parallel unit
    !             is then a whole transform, which is ~L16-19 independent
    !             units per field per step -- embarrassingly parallel, no
    !             synchronization inside.
    !
    ! Mode B is what section 4.3's "parallelise over (level, wavenumber) pairs"
    ! actually describes, and it is what FastEarth3D already does for the same
    ! reason (fe_sht's clone_cfg). This benchmark measures both so the choice
    ! is made on numbers.
    !
    ! Reports time per scalar transform (one synthesis + one analysis counts as
    ! two), so the truncations are directly comparable.

    use omp_lib
    use aeros_defs,     only : dp, wp, wp_sh
    use aeros_spectral

    implicit none

    integer, parameter :: truncs(3) = [31, 42, 85]
    integer, parameter :: nlev      = 16
    integer, parameter :: nrep      = 20

    integer :: threads(5)
    integer :: it, ith, nt, maxth

    maxth = omp_get_max_threads()
    threads = [1, 2, 4, 8, maxth]

    write(*,*) ""
    write(*,*) "bench_sht:: nlev = ", nlev, " repeats = ", nrep
    write(*,*) "            max threads available = ", maxth
    write(*,*) "            wp = ", storage_size(1.0_wp)/8, " bytes"
    write(*,*) ""

    do it = 1, size(truncs)

        write(*,"(a)") ""
        write(*,"(a,i0)") " === T", truncs(it)
        write(*,"(a)") "  threads   modeA us/transf   modeB us/transf   A speedup   B speedup"

        call run_truncation(truncs(it), threads, maxth)

    end do

contains

    subroutine run_truncation(trunc, threads, maxth)

        implicit none

        integer, intent(in) :: trunc
        integer, intent(in) :: threads(:)
        integer, intent(in) :: maxth

        real(dp) :: ta, tb, ta1, tb1
        integer  :: ith, nt

        ta1 = -1.0_dp
        tb1 = -1.0_dp

        do ith = 1, size(threads)
            nt = threads(ith)
            if (nt > maxth) cycle
            if (ith > 1) then
                if (nt == threads(ith-1)) cycle    ! maxth may repeat 8
            end if

            ta = time_mode_a(trunc, nt)
            tb = time_mode_b(trunc, nt)

            if (ta1 < 0.0_dp) ta1 = ta
            if (tb1 < 0.0_dp) tb1 = tb

            write(*,"(i9,f18.2,f18.2,f12.2,f12.2)") nt, ta*1.0e6_dp, tb*1.0e6_dp, ta1/ta, tb1/tb
        end do

        return

    end subroutine run_truncation

    real(dp) function time_mode_a(trunc, nt) result(tper)
        ! SHTns threads one transform; the level loop is serial.

        implicit none

        integer, intent(in) :: trunc, nt

        type(aeros_sht_class) :: sht
        real(wp),       allocatable :: field(:,:,:)
        complex(wp_sh), allocatable :: coeffs(:,:)
        real(dp) :: t0, t1
        integer  :: k, r

        call aeros_sht_init(sht, trunc, nthreads=nt, quick=.TRUE.)

        allocate(field(sht%nlon,sht%nlat,nlev))
        allocate(coeffs(sht%nlm,nlev))
        call seed(field, coeffs, sht, nlev)

        ! warm up
        do k = 1, nlev
            call aeros_sht_synthesis(sht, coeffs(:,k), field(:,:,k))
        end do

        t0 = omp_get_wtime()
        do r = 1, nrep
            do k = 1, nlev
                call aeros_sht_synthesis(sht, coeffs(:,k), field(:,:,k))
                call aeros_sht_analysis(sht, field(:,:,k), coeffs(:,k))
            end do
        end do
        t1 = omp_get_wtime()

        tper = (t1 - t0)/real(2*nrep*nlev, dp)

        deallocate(field, coeffs)
        call aeros_sht_end(sht)

        return

    end function time_mode_a

    real(dp) function time_mode_b(trunc, nt) result(tper)
        ! SHTns single-threaded; aeros threads the level loop, one config per
        ! thread. Configs are created serially -- SHTns config creation is not
        ! thread-safe (it plans FFTs).

        implicit none

        integer, intent(in) :: trunc, nt

        type(aeros_sht_class), allocatable :: sht(:)
        real(wp),       allocatable :: field(:,:,:)
        complex(wp_sh), allocatable :: coeffs(:,:)
        real(dp) :: t0, t1
        integer  :: k, r, i, tid

        allocate(sht(nt))
        do i = 1, nt
            call aeros_sht_init(sht(i), trunc, nthreads=1, quick=.TRUE.)
        end do

        allocate(field(sht(1)%nlon,sht(1)%nlat,nlev))
        allocate(coeffs(sht(1)%nlm,nlev))
        call seed(field, coeffs, sht(1), nlev)

        do k = 1, nlev
            call aeros_sht_synthesis(sht(1), coeffs(:,k), field(:,:,k))
        end do

        t0 = omp_get_wtime()
        do r = 1, nrep
            !$omp parallel do num_threads(nt) private(k,tid) schedule(static)
            do k = 1, nlev
                tid = omp_get_thread_num() + 1
                call aeros_sht_synthesis(sht(tid), coeffs(:,k), field(:,:,k))
                call aeros_sht_analysis(sht(tid), field(:,:,k), coeffs(:,k))
            end do
            !$omp end parallel do
        end do
        t1 = omp_get_wtime()

        tper = (t1 - t0)/real(2*nrep*nlev, dp)

        deallocate(field, coeffs)
        do i = 1, nt
            call aeros_sht_end(sht(i))
        end do
        deallocate(sht)

        return

    end function time_mode_b

    subroutine seed(field, coeffs, sht, nlev)

        implicit none

        real(wp),       intent(out) :: field(:,:,:)
        complex(wp_sh), intent(out) :: coeffs(:,:)
        type(aeros_sht_class), intent(in) :: sht
        integer, intent(in) :: nlev

        integer  :: lm, k, l
        real(dp) :: amp

        field = 0.0_wp
        do k = 1, nlev
            do lm = 1, sht%nlm
                l = sht%l_of_lm(lm)
                amp = 1.0_dp/real((l+1)*(l+1), dp)
                coeffs(lm,k) = cmplx(amp*cos(real(l+k,dp)), amp*sin(real(l+2*k,dp)), wp_sh)
            end do
        end do

        return

    end subroutine seed

end program bench_sht
