program aeros_run
    ! Standalone driver for the global core -- the canonical minimal example
    ! (docs/coupling_api.md sec. 14, decision 6): initialise through
    ! aeros_lib_init, drive the time loop with aeros_lib_step on the output
    ! interval, finalise through aeros_lib_finalize. This is the same three
    ! calls a coupled host makes; the driver's only extra job is to write the
    ! state after each cycle.
    !
    ! runme stages this executable (.runme/info.json) with the case file as
    ! its first argument -- runme's `par_path_as_argument = true` contract --
    ! so a run directory is self-describing: the namelist next to the
    ! executable is the one that produced the output.

    use aeros_defs, only : dp, wp, io_unit_err, year_seconds
    use aeros,      only : aeros_class, aeros_lib_init, aeros_lib_step, &
                            aeros_lib_finalize, aeros_print_config
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
    real(dp) :: dt_out_sec
    integer  :: n, nstep, ierr
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
    call aeros_lib_init(ams, path_par, defaults_file=defaults_file, &
                        ic_source="resting", ierr=ierr)
    if (ierr /= 0) then
        write(io_unit_err,*) "aeros_run:: aeros_lib_init failed, ierr = ", ierr
        stop 1
    end if
    ! Wall-clock origin. aeros_lib_init leaves ams%time = 0; if the case wants
    ! to label its output from a non-zero time_init [yr] the driver sets it
    ! here. lib_step advances ams%time by nsub*dt/year_seconds each cycle.
    ams%time_init = time_init
    ams%time      = time_init

    call aeros_print_config(ams)

    call aeros_io_init(io)
    call aeros_write_init(trim(file_out), ams%grd, ams%vgrid, time_init)
    call aeros_write_state(io, trim(file_out), ams%grd, ams%now, time_init, 1)

    ! -- Time loop, on the output interval. aeros_lib_step subdivides
    ! dt_couple = dt_out into its own dt (nint(dt_couple / par%dt) sub-steps),
    ! so this driver never has to know the atmosphere's internal timestep --
    ! which is what makes the same call usable from a coupled driver.
    nstep      = nint((time_end - time_init)/dt_out)
    dt_out_sec = real(dt_out, dp)*real(year_seconds, dp)

    do n = 1, nstep
        time = time_init + real(n, wp)*dt_out
        call aeros_lib_step(ams, dt_couple=dt_out_sec, ierr=ierr)
        if (ierr /= 0) then
            write(io_unit_err,*) "aeros_run:: aeros_lib_step failed, ierr = ", ierr, &
                                    " at step ", n
            stop 1
        end if
        call aeros_write_state(io, trim(file_out), ams%grd, ams%now, time, n+1)
    end do

    call aeros_lib_finalize(ams, ierr=ierr)
    if (ierr /= 0) then
        write(io_unit_err,*) "aeros_run:: aeros_lib_finalize failed, ierr = ", ierr
        stop 1
    end if

    write(*,*) "aeros_run:: done. Output written to ", trim(file_out)

end program aeros_run
