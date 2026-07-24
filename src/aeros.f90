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
    ! At M0 `aeros_update` advances the clock and does nothing else. The
    ! dynamical core lands at M1.

    use aeros_defs,     only : sp, dp, wp, wp_sh, io_unit_err, &
                                aeros_param_class, aeros_grid_class, aeros_state_class, &
                                aeros_par_load
    use aeros_spectral, only : aeros_sht_class, aeros_sht_init, aeros_sht_end
    use aeros_grid,     only : aeros_grid_init, aeros_grid_end
    use aeros_state,    only : aeros_state_alloc, aeros_state_end

    implicit none

    private

    type aeros_class
        type(aeros_param_class) :: par
        type(aeros_sht_class)   :: sht
        type(aeros_grid_class)  :: grd
        type(aeros_state_class) :: now

        real(wp) :: time      ! current model time [yr]
    end type aeros_class

    public :: aeros_class
    public :: aeros_init
    public :: aeros_update
    public :: aeros_end
    public :: aeros_print_config

contains

    subroutine aeros_init(ams, filename, group, time_init)
        ! Configure the transform, the grid and the state from a namelist.

        implicit none

        type(aeros_class), intent(inout) :: ams
        character(len=*),  intent(in)    :: filename
        character(len=*),  intent(in), optional :: group
        real(wp),          intent(in), optional :: time_init

        call aeros_par_load(ams%par, filename, group)

        call aeros_sht_init(ams%sht, ams%par%trunc, nlon=ams%par%nlon, nlat=ams%par%nlat, &
                                nthreads=ams%par%nthreads)

        call aeros_grid_init(ams%grd, ams%sht)

        call aeros_state_alloc(ams%now, ams%grd, ams%sht%nlm, ams%par%nlev)

        ams%time = 0.0_wp
        if (present(time_init)) ams%time = time_init

        return

    end subroutine aeros_init

    subroutine aeros_update(ams, time)
        ! Advance the model to `time` [yr].
        !
        ! M0 STUB: advances the clock only. The semi-implicit leapfrog loop
        ! (docs/design.md section 3.2, 4) replaces the body at M1; the
        ! signature is expected to survive that, since a host driver only ever
        ! needs to say "catch up to this time".

        implicit none

        type(aeros_class), intent(inout) :: ams
        real(wp), intent(in) :: time

        if (time < ams%time) then
            write(io_unit_err,*) "aeros_update:: error: cannot step backwards, from ", &
                                    ams%time, " to ", time
            stop 1
        end if

        ams%time = time

        return

    end subroutine aeros_update

    subroutine aeros_end(ams)

        implicit none

        type(aeros_class), intent(inout) :: ams

        call aeros_state_end(ams%now)
        call aeros_grid_end(ams%grd)
        call aeros_sht_end(ams%sht)

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
        write(iou,"(a24,i8)")    "truncation  T", ams%sht%trunc
        write(iou,"(a24,i8)")    "nlon         ", ams%sht%nlon
        write(iou,"(a24,i8)")    "nlat         ", ams%sht%nlat
        write(iou,"(a24,i8)")    "columns      ", ams%grd%ncol
        write(iou,"(a24,i8)")    "nlev         ", ams%par%nlev
        write(iou,"(a24,i8)")    "coeffs/level ", ams%sht%nlm
        write(iou,"(a24,f8.1)")  "dt [s]       ", ams%par%dt
        write(iou,"(a24,i8)")    "wp bytes     ", storage_size(1.0_wp)/8
        write(iou,"(a24,i8)")    "wp_sh bytes  ", storage_size(1.0_wp_sh)/8
        write(iou,*) ""

        return

    end subroutine aeros_print_config

end module aeros
