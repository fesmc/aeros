program held_suarez
    ! The Held & Suarez (1994) benchmark run.
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the intercomparison
    !   of the dynamical cores of atmospheric general circulation models.
    !   Bull. Amer. Meteor. Soc., 75, 1825-1830.
    !
    ! docs/design.md section 7 makes this the validation for M1. It is a
    ! statement about the CLIMATE of a dynamical core, so the run is long and
    ! the output is statistics: 1200 days, with the first 200 discarded as
    ! spin-up and the last 1000 averaged.
    !
    ! A separate driver rather than a switch in aeros_run.x because the
    ! benchmark needs things a general run does not -- a fixed vertical grid, a
    ! seeded initial condition, an averaging window, zonal-mean output -- and
    ! because a benchmark whose configuration can be half-applied is not a
    ! benchmark. Times here are in DAYS, which is the unit the paper and every
    ! comparison are written in; aeros' own interface is in years and the
    ! conversion happens once, below.
    !
    ! Times in the `hs` namelist group are days; everything else comes from the
    ! usual `aeros` and `aeros_vert` groups, so the truncation and the vertical
    ! grid are configured exactly as any other run.

    use aeros_defs,        only : dp, wp, io_unit_err, year_seconds
    use aeros,             only : aeros_class, aeros_init, aeros_update, aeros_end, &
                                    aeros_set_isothermal, aeros_diagnose, &
                                    aeros_print_config
    use aeros_budget,      only : aeros_budget_class, aeros_budget_calc, aeros_budget_report
    use aeros_diagnostics, only : aeros_diag_class, aeros_diag_init, aeros_diag_end, &
                                    aeros_diag_accumulate, aeros_diag_finalize, &
                                    aeros_diag_summary
    use aeros_held_suarez, only : aeros_hs_perturb
    use aeros_io,          only : aeros_io_class, aeros_io_init, aeros_write_init, &
                                    aeros_write_state, aeros_write_diag
    use nml,               only : nml_read

    implicit none

    type(aeros_class)        :: ams
    type(aeros_io_class)     :: io
    type(aeros_diag_class)   :: dgn
    type(aeros_budget_class) :: bud0, bud1

    character(len=512) :: path_par, file_out, file_diag

    real(wp) :: day_end, day_avg_start, dt_sample, dt_out, t_init, perturb
    real(wp) :: day, time, day_out_next
    real(wp), allocatable :: phis(:,:)
    integer  :: n, nsample, nout

    ! -- Parameter file: argv(1), or the in-tree default when run by hand.
    if (command_argument_count() >= 1) then
        call get_command_argument(1, path_par)
    else
        path_par = "par/held_suarez.nml"
    end if

    call nml_read(path_par, "hs", "day_end",       day_end)
    call nml_read(path_par, "hs", "day_avg_start", day_avg_start)
    call nml_read(path_par, "hs", "dt_sample",     dt_sample)
    call nml_read(path_par, "hs", "dt_out",        dt_out)
    call nml_read(path_par, "hs", "t_init",        t_init)
    call nml_read(path_par, "hs", "perturb",       perturb)
    call nml_read(path_par, "hs", "file_out",      file_out)
    call nml_read(path_par, "hs", "file_diag",     file_diag)

    if (dt_sample <= 0.0_wp) then
        write(io_unit_err,*) "held_suarez:: error: dt_sample must be > 0 days, got ", dt_sample
        stop 1
    end if
    if (day_avg_start >= day_end) then
        write(io_unit_err,*) "held_suarez:: error: the averaging window is empty; ", &
                                "day_avg_start = ", day_avg_start, " >= day_end = ", day_end
        stop 1
    end if

    ! === Initialize ==========================================================

    call aeros_init(ams, path_par, group="aeros", time_init=0.0_wp)

    if (.not. ams%par%held_suarez) then
        write(io_unit_err,*) "held_suarez:: error: `held_suarez` is .FALSE. in the `aeros` "// &
                                "namelist group. This driver exists to run that forcing; "// &
                                "running it without would silently produce an adiabatic "// &
                                "result labelled as a benchmark."
        stop 1
    end if

    call aeros_print_config(ams)

    ! Initial condition: motionless and isothermal, well below the equilibrium
    ! the forcing relaxes toward. Held & Suarez do not specify one -- the point
    ! of the benchmark is the equilibrated statistics -- but starting cold and
    ! at rest is conventional and makes the 200-day spin-up interpretable: the
    ! free troposphere relaxes at 1/40 day, so 200 days is five e-folds.
    call aeros_set_isothermal(ams, temp=t_init)

    ! Break the zonal symmetry. Without a seed the only thing that can start
    ! the baroclinic instability is round-off, which works but on an
    ! unpredictable schedule. The seed is deterministic so that two runs of the
    ! same configuration give the same numbers.
    call aeros_hs_perturb(ams%now%spec, ams%pool%sht(1), perturb)
    call aeros_diagnose(ams)

    ! No topography: the benchmark is explicit that the lower boundary is flat.
    ! Set here rather than relied upon, so the run states its own assumption.
    allocate(phis(ams%grd%nlon, ams%grd%nlat))
    phis = 0.0_wp

    call aeros_io_init(io)
    call aeros_diag_init(dgn, ams%grd, ams%par%nlev)

    call aeros_write_init(trim(file_out), ams%grd, ams%vgrid, 0.0_wp, units_time="days")
    call aeros_write_state(io, trim(file_out), ams%grd, ams%now, 0.0_wp, 1)

    call aeros_budget_calc(bud0, ams%vgrid, ams%grd, ams%now, phis)

    ! === Integrate ===========================================================

    nsample      = nint(day_end/dt_sample)
    nout         = 1
    day_out_next = dt_out

    write(*,*) ""
    write(*,"(a,f9.1,a,f9.1,a)") " held_suarez:: integrating ", day_end, &
                                    " days, averaging from day ", day_avg_start, ""
    write(*,"(a,f9.2,a,i0,a)")   "               sampling every ", dt_sample, &
                                    " days (", nsample, " samples)"
    write(*,*) ""
    ! The last column is the mass the fixer has put back so far. With the fixer
    ! off it is identically zero; with it on, the "mass drift" column goes to
    ! machine zero and this one carries the fixer's running workload. It is NOT
    ! the drift an unfixed run would show -- see the note after the budgets.
    write(*,"(a)") "        day      max|u|     min[T]     max[T]        mass drift"// &
                    "       fixer removed"

    do n = 1, nsample

        day  = real(n, wp)*dt_sample
        time = day*86400.0_wp/real(year_seconds, wp)

        call aeros_update(ams, time)

        ! aeros_update leaves the grid fields consistent, so the sample costs
        ! no extra transforms -- which is why the sampling interval and the
        ! update interval are the same thing here.
        if (day >= day_avg_start) call aeros_diag_accumulate(dgn, ams%grd, ams%now, day)

        if (day >= day_out_next - 0.5_wp*dt_sample) then
            nout = nout + 1
            call aeros_write_state(io, trim(file_out), ams%grd, ams%now, day, nout)
            day_out_next = day_out_next + dt_out
        end if

        ! A progress line often enough to see a run go wrong, rarely enough not
        ! to bury the result. The mass drift is the cheap canary: it is the one
        ! quantity that should not move at all.
        if (mod(n, max(1, nsample/40)) == 0) then
            call aeros_budget_calc(bud1, ams%vgrid, ams%grd, ams%now, phis)
            write(*,"(f11.1,3f11.2,2es18.3)") day, &
                    maxval(abs(ams%now%u)), minval(ams%now%temp_g), maxval(ams%now%temp_g), &
                    (bud1%mass - bud0%mass)/bud0%mass, &
                    exp(-ams%ts%lnr_cum) - 1.0_dp
        end if

    end do

    ! === Report ==============================================================

    ! Only MASS is expected to be conserved here. The Newtonian relaxation is
    ! an energy source and the Rayleigh drag an angular-momentum sink, both by
    ! design, so their drift rates measure the forcing rather than the
    ! integrator -- and over the spin-up they are large, because the atmosphere
    ! is being heated from t_init toward equilibrium. The adiabatic
    ! conservation test is tests/test_timestep.f90.
    call aeros_budget_calc(bud1, ams%vgrid, ams%grd, ams%now, phis)
    write(*,*) ""
    write(*,"(a)") " Budgets. The forcing is an energy source and a momentum sink by"
    write(*,"(a)") " construction, so only the mass drift says anything about the integrator."
    call aeros_budget_report(bud1, bud0, real(day_end, dp)*86400.0_dp, &
                                "Held-Suarez run")

    if (ams%ts%mass_fixer) then
        write(*,*) ""
        write(*,"(a)") " The mass fixer was ON, so the mass drift above is the fixer's"
        write(*,"(a)") " residual and not the integrator's. What it had to put back:"
        write(*,"(a,es12.3)")   "   relative mass restored, total    ", &
                                    exp(-ams%ts%lnr_cum) - 1.0_dp
        write(*,"(a,es12.3,a)") "   the same as a rate               ", &
                                    (exp(-ams%ts%lnr_cum) - 1.0_dp) &
                                        /(real(day_end, dp)/365.25_dp), " per year"
        write(*,"(a)") " This is the fixer's WORKLOAD, not the drift an unfixed run would"
        write(*,"(a)") " show -- measured, it runs ~2x larger, because the unfixed run's"
        write(*,"(a)") " mass oscillation partly cancels over time while the fixer pays"
        write(*,"(a)") " for it every step. For the drift itself, run an unfixed twin."
    end if

    call aeros_diag_finalize(dgn, ams%vgrid)
    call aeros_diag_summary(dgn, ams%grd)

    call aeros_write_diag(io, trim(file_diag), ams%grd, ams%vgrid, dgn)

    call aeros_diag_end(dgn)
    call aeros_end(ams)

    write(*,*) "held_suarez:: done."
    write(*,*) "  state       -> ", trim(file_out)
    write(*,*) "  diagnostics -> ", trim(file_diag)

end program held_suarez
