module aeros
    ! Public facade: the one module a host program or coupled driver uses.
    !
    ! docs/design.md section 10.8 settles this as library-first -- aeros is
    ! built as libaeros.a and called, rather than being a program that other
    ! things are bolted onto. The interface below is therefore the contract
    ! that matters, and it is deliberately the same three calls a CLIMBER-X- or
    ! yelmox-style driver already expects:
    !
    !     aeros_init(ams, filename)      allocate and configure
    !     aeros_update(ams, time)        advance to `time`
    !     aeros_end(ams)                 release
    !
    ! The signature survived M1 unchanged, which was the point of choosing it
    ! at M0: a host driver says "catch up to this time" and never learns the
    ! atmosphere's internal timestep.
    !
    ! === Time ================================================================
    !
    ! The interface is in YEARS and the dynamics run in SECONDS, and the two
    ! cannot in general be made to divide. aeros resolves that by deriving the
    ! step count from the ABSOLUTE time since init rather than from each
    ! interval:
    !
    !     n_target = nint( (time - time_init) * year_seconds / dt )
    !
    ! and taking n_target - n_done steps. Rounding then never accumulates: a
    ! sequence of calls on an interval that is 1.4 steps long takes 1, 2, 1, 3,
    ! ... steps and stays within half a step of the requested time forever,
    ! where rounding each interval separately would drift by 0.4 steps per
    ! call. `ams%time` is set to the requested time exactly, so the clock a
    ! coupled driver sees is the one it asked for.

    use aeros_defs,     only : sp, dp, wp, wp_sh, io_unit_err, year_seconds, p0, &
                                aeros_param_class, aeros_grid_class, aeros_state_class, &
                                aeros_par_load
    use aeros_spectral, only : aeros_sht_class, aeros_sht_pool_class, &
                                aeros_sht_pool_init, aeros_sht_pool_end, aeros_sht_lm
    use aeros_grid,     only : aeros_grid_init, aeros_grid_end
    use aeros_state,    only : aeros_state_alloc, aeros_state_end, aeros_spec_zero
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_load, aeros_vgrid_end, &
                                aeros_vgrid_print
    use aeros_timestep, only : aeros_timestep_class, aeros_timestep_init, &
                                aeros_timestep_end, aeros_timestep_step, &
                                aeros_timestep_diagnose, aeros_timestep_print

    implicit none

    private

    type aeros_class
        type(aeros_param_class)    :: par
        type(aeros_sht_pool_class) :: pool
        type(aeros_grid_class)     :: grd
        type(aeros_vgrid_class)    :: vgrid
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts

        real(wp) :: time                  ! current model time [yr]
        real(wp) :: time_init             ! time at aeros_init [yr]
        integer  :: nstep = 0             ! dynamics steps taken since init
    end type aeros_class

    public :: aeros_class
    public :: aeros_init
    public :: aeros_update
    public :: aeros_end
    public :: aeros_print_config

contains

    subroutine aeros_init(ams, filename, group, time_init)
        ! Configure the transform pool, the grids, the state and the integrator.

        implicit none

        type(aeros_class), intent(inout) :: ams
        character(len=*),  intent(in)    :: filename
        character(len=*),  intent(in), optional :: group
        real(wp),          intent(in), optional :: time_init

        call aeros_par_load(ams%par, filename, group)

        ! One SHTns config per thread; see aeros_sht_pool_class for why aeros
        ! threads the level loop rather than letting SHTns thread a transform.
        ! cache=.TRUE. because the pool builds nthreads IDENTICAL configs and
        ! SHTns does not remember its auto-tuning between them -- without the
        ! cache a 10-thread pool pays the tuning ten times.
        call aeros_sht_pool_init(ams%pool, ams%par%trunc, &
                                    nlon=ams%par%nlon, nlat=ams%par%nlat, &
                                    nthreads=ams%par%nthreads, cache=.TRUE.)

        call aeros_grid_init(ams%grd, ams%pool%sht(1))

        call aeros_vgrid_load(ams%vgrid, ams%par%nlev, filename)

        call aeros_state_alloc(ams%now, ams%grd, ams%pool%sht(1)%nlm, ams%par%nlev)

        call aeros_timestep_init(ams%ts, ams%par, ams%pool, ams%grd, ams%vgrid)

        ! Initial condition. A resting isothermal atmosphere over a flat
        ! surface -- the only state the dry core can be started from at M1
        ! without inventing a climatology, and the one whose correct behaviour
        ! (staying exactly where it is) the acceptance tests measure. It is
        ! also the state Held-Suarez is spun up from at M1.5, perturbed.
        !
        ! NOT the zeroed state aeros_state_alloc leaves: ln(p_s) = 0 there
        ! means a surface pressure of 1 Pa, which is not a small error but a
        ! different planet.
        call init_isothermal(ams)

        ams%time_init = 0.0_wp
        if (present(time_init)) ams%time_init = time_init
        ams%time  = ams%time_init
        ams%nstep = 0

        return

    end subroutine aeros_init

    subroutine init_isothermal(ams)
        ! A motionless isothermal atmosphere at the reference state.

        implicit none

        ! `target` so the pool's configs can be pointed at: gfortran refuses a
        ! pointer assignment to a component of a non-target dummy. Harmless
        ! here -- the pointer does not outlive the call.
        type(aeros_class), intent(inout), target :: ams

        type(aeros_sht_class), pointer :: s
        real(dp) :: y00
        integer  :: k, lm00

        s => ams%pool%sht(1)

        call aeros_spec_zero(ams%now%spec)

        ! Orthonormal harmonics with no Condon-Shortley phase (aeros_spectral):
        ! Y_00 = 1/sqrt(4 pi), so a uniform field of value X is the single
        ! coefficient X*sqrt(4 pi).
        y00  = sqrt(4.0_dp*acos(-1.0_dp))
        lm00 = aeros_sht_lm(s, 0, 0)

        do k = 1, ams%vgrid%nlev
            ams%now%spec%temp(lm00,k) = cmplx(real(ams%vgrid%t_ref(k),dp)*y00, 0.0_dp, wp_sh)
        end do
        ams%now%spec%lnps(lm00) = cmplx(log(real(ams%vgrid%ps_ref,dp))*y00, 0.0_dp, wp_sh)

        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)

        return

    end subroutine init_isothermal

    subroutine aeros_update(ams, time)
        ! Advance the model to `time` [yr].

        implicit none

        type(aeros_class), intent(inout) :: ams
        real(wp), intent(in) :: time

        integer :: n, nstep_target, nstep_do

        if (time < ams%time) then
            write(io_unit_err,*) "aeros_update:: error: cannot step backwards, from ", &
                                    ams%time, " to ", time
            stop 1
        end if

        ! Step count from the absolute elapsed time, so the rounding error is
        ! bounded rather than accumulated. See the module header.
        nstep_target = nint(real(time - ams%time_init, dp)*year_seconds &
                                /real(ams%par%dt, dp))
        nstep_do     = nstep_target - ams%nstep

        do n = 1, nstep_do
            call aeros_timestep_step(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)
        end do

        ams%nstep = max(nstep_target, ams%nstep)
        ams%time  = time

        ! Leave the grid-space fields consistent with the spectral state the
        ! caller now holds. The integrator does not do this every step -- see
        ! aeros_timestep_diagnose -- so a caller that reads ams%now%u without
        ! going through aeros_update would be reading whatever the last update
        ! left. Doing it here makes the facade's contract simply "after
        ! aeros_update, ams%now is a complete, consistent state".
        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)

        return

    end subroutine aeros_update

    subroutine aeros_end(ams)

        implicit none

        type(aeros_class), intent(inout) :: ams

        call aeros_timestep_end(ams%ts)
        call aeros_state_end(ams%now)
        call aeros_vgrid_end(ams%vgrid)
        call aeros_grid_end(ams%grd)
        call aeros_sht_pool_end(ams%pool)

        ams%time = 0.0_wp; ams%time_init = 0.0_wp; ams%nstep = 0

        return

    end subroutine aeros_end

    subroutine aeros_print_config(ams, io_unit)
        ! Summarize the resolved configuration.
        !
        ! Prints the REALIZED grid, not the requested one: nlon/nlat default to
        ! -1 in the namelist and are chosen by the truncation, and SHTns is
        ! free to adjust what it is given. A run log that reports the request
        ! rather than the result is how a silently different resolution gets
        ! into a 10^5 yr integration.

        implicit none

        type(aeros_class), intent(in) :: ams
        integer, intent(in), optional :: io_unit

        integer :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        write(iou,*) ""
        write(iou,*) "== aeros configuration =="
        write(iou,"(a24,i8)")    "truncation  T", ams%pool%sht(1)%trunc
        write(iou,"(a24,i8)")    "nlon         ", ams%pool%sht(1)%nlon
        write(iou,"(a24,i8)")    "nlat         ", ams%pool%sht(1)%nlat
        write(iou,"(a24,i8)")    "columns      ", ams%grd%ncol
        write(iou,"(a24,i8)")    "nlev         ", ams%par%nlev
        write(iou,"(a24,i8)")    "coeffs/level ", ams%pool%sht(1)%nlm
        write(iou,"(a24,i8)")    "threads      ", ams%pool%nthreads
        write(iou,"(a24,f8.1)")  "dt [s]       ", ams%par%dt
        write(iou,"(a24,i8)")    "wp bytes     ", storage_size(1.0_wp)/8
        write(iou,"(a24,i8)")    "wp_sh bytes  ", storage_size(1.0_wp_sh)/8
        write(iou,*) ""

        call aeros_vgrid_print(ams%vgrid, iou)
        call aeros_timestep_print(ams%ts, iou)

        return

    end subroutine aeros_print_config

end module aeros
