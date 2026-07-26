program aeros_run
    ! Standalone driver for the global core.
    !
    ! This is the executable `runme` stages and runs (.runme/info.json). It
    ! takes the parameter file as its first argument -- runme's
    ! `par_path_as_argument = true` contract -- so a staged run directory is
    ! self-describing: the namelist sitting next to the executable is the one
    ! that produced the output.
    !
    ! M0: initialize, step the clock, write output, finish. The time loop is
    ! already shaped the way M1 needs it, so the dynamical core lands inside
    ! aeros_update without the driver changing.

    use aeros_defs, only : wp, io_unit_err
    use aeros,      only : aeros_class, aeros_init, aeros_update, aeros_end, &
                            aeros_print_config
    use aeros_io,   only : aeros_io_class, aeros_io_init, aeros_write_init, &
                            aeros_write_state
    use nml,        only : nml_read

    implicit none

    type(aeros_class)    :: ams
    type(aeros_io_class) :: io

    character(len=512) :: path_par
    character(len=512) :: file_out
    character(len=*), parameter :: defaults_file = "input/aeros_defaults.nml"

    real(wp) :: time_init, time_end, dt_out, time
    integer  :: n, nstep
    logical  :: have_defaults

    ! -- Parameter file: argv(1), or the in-tree default when run by hand.
    if (command_argument_count() >= 1) then
        call get_command_argument(1, path_par)
    else
        path_par = "par/aeros.nml"
    end if

    ! -- The defaults schema: every parameter's fallback, so the case file
    !    (path_par) need only list its overrides. Staged next to the executable
    !    (input/ is symlinked into a run dir); required, since the case file is
    !    no longer expected to be complete.
    inquire(file=defaults_file, exist=have_defaults)
    if (.not. have_defaults) then
        write(io_unit_err,*) "aeros_run:: error: defaults file not found: ", defaults_file
        stop 1
    end if

    call nml_read(path_par, "ctrl", "time_init", time_init, defaults_file=defaults_file)
    call nml_read(path_par, "ctrl", "time_end",  time_end,  defaults_file=defaults_file)
    call nml_read(path_par, "ctrl", "dt_out",    dt_out,    defaults_file=defaults_file)
    call nml_read(path_par, "ctrl", "file_out",  file_out,  defaults_file=defaults_file)

    if (dt_out <= 0.0_wp) then
        write(io_unit_err,*) "aeros_run:: error: dt_out must be > 0, got ", dt_out
        stop 1
    end if

    ! -- Initialize.
    call aeros_init(ams, path_par, group="aeros", time_init=time_init, &
                    defaults_file=defaults_file)
    call aeros_print_config(ams)

    call aeros_io_init(io)
    call aeros_write_init(trim(file_out), ams%grd, ams%vgrid, time_init)
    call aeros_write_state(io, trim(file_out), ams%grd, ams%now, time_init, 1)

    ! -- Time loop, on the output interval. aeros_update is responsible for
    ! subdividing into its own dt: a host driver should never have to know the
    ! atmosphere's internal timestep, which is what makes the same call usable
    ! from a coupled driver at M4.
    nstep = nint((time_end - time_init)/dt_out)

    do n = 1, nstep
        time = time_init + real(n, wp)*dt_out
        call aeros_update(ams, time)
        call aeros_write_state(io, trim(file_out), ams%grd, ams%now, time, n+1)
    end do

    call aeros_end(ams)

    write(*,*) "aeros_run:: done. Output written to ", trim(file_out)

end program aeros_run
