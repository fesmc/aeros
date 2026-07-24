module aeros_io
    ! netCDF output, driven by variable tables.
    !
    ! What a run writes is decided by two markdown tables under `input/`, read
    ! at init through fesm-utils' `variable_io` -- the same mechanism yelmo
    ! uses, and for the same reason: the list of output variables, their units
    ! and their long names belong next to each other in one editable place
    ! rather than scattered through call sites, and changing an output should
    ! not be a recompile.
    !
    !   input/aeros-variables-state.md   instantaneous state, on the grid
    !   input/aeros-variables-diag.md    zonal-mean time-mean statistics
    !
    ! Two tables rather than one because they have different shapes, different
    ! sources and different lifetimes: the state is (lon,lat,lev,time) taken
    ! from aeros_state_class every output interval, the diagnostics are
    ! (lat,lev) taken from aeros_diag_class once at the end of a run.
    !
    ! Each table has a dispatcher below, which ERRORS on a variable name it
    ! does not recognize. A name in the table with no case in the dispatcher is
    ! a typo rather than a request, and silently skipping it would produce a
    ! file that is quietly missing a field.
    !
    ! Deliberately takes the grid, the state and the diagnostics directly
    ! rather than the top-level aeros_class: output has no business knowing
    ! about the SHTns handle, and keeping this below the facade in the
    ! dependency order means the facade can use it rather than the other way
    ! round.

    use aeros_defs,        only : sp, dp, wp, io_unit_err, &
                                    aeros_grid_class, aeros_state_class
    use aeros_vertical,    only : aeros_vgrid_class
    use aeros_diagnostics, only : aeros_diag_class
    use variable_io
    use ncio

    implicit none

    private

    type aeros_io_class
        ! The loaded tables, held for the life of a run so that a staged run
        ! directory is not re-read on every output step.
        type(var_io_type), allocatable :: state(:)
        type(var_io_type), allocatable :: diag(:)
    end type aeros_io_class

    public :: aeros_io_class
    public :: aeros_io_init
    public :: aeros_write_init
    public :: aeros_write_state
    public :: aeros_write_diag

contains

    subroutine aeros_io_init(io, path)
        ! Load the variable tables.
        !
        ! `path` is the directory holding them, `input` by default -- which is
        ! symlinked into a staged run directory by runme (.runme/info.json
        ! `links`), so the same relative path works in-tree and in a run.

        implicit none

        type(aeros_io_class), intent(inout) :: io
        character(len=*), intent(in), optional :: path

        character(len=512) :: dir

        dir = "input"
        if (present(path)) dir = trim(path)

        call load_var_io_table(io%state, trim(dir)//"/aeros-variables-state.md")
        call load_var_io_table(io%diag,  trim(dir)//"/aeros-variables-diag.md")

        return

    end subroutine aeros_io_init

    subroutine aeros_write_init(filename, grd, vg, time_init, units_time)
        ! Create an output file and write its dimensions and static fields.

        implicit none

        character(len=*), intent(in) :: filename
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_vgrid_class), intent(in) :: vg
        real(wp), intent(in) :: time_init
        character(len=*), intent(in), optional :: units_time

        character(len=56) :: tunits

        tunits = "years"
        if (present(units_time)) tunits = trim(units_time)

        call nc_create(filename)

        call nc_write_dim(filename,"lon",x=grd%lon,units="degrees_east")
        call nc_write_dim(filename,"lat",x=grd%lat,units="degrees_north")
        call nc_write_dim(filename,"lev",x=1,dx=1,nx=vg%nlev,units="1")
        call nc_write_dim(filename,"time",x=time_init,dx=1.0_wp,nx=1, &
                            units=trim(tunits),unlimited=.TRUE.)

        ! Static grid fields. `area` is written because every conservation
        ! diagnostic downstream is a weighted sum with it, and a reader that
        ! recomputes the weights from lat is the classic way to get a
        ! conservation check that passes against the wrong quadrature.
        call nc_write(filename,"area",grd%area,dim1="lon",dim2="lat", &
                        units="m2",long_name="Grid cell area")
        call nc_write(filename,"coriolis",grd%coriolis,dim1="lon",dim2="lat", &
                        units="s-1",long_name="Coriolis parameter")

        ! The vertical coordinate itself. Without A and B a level index cannot
        ! be turned back into a pressure, which makes every 3-D field in the
        ! file uninterpretable on its own.
        call nc_write(filename,"sigma",vg%sigma_full,dim1="lev", &
                        units="1",long_name="Sigma at full levels, at ps_ref")
        call nc_write(filename,"p_ref",vg%p_ref_full,dim1="lev", &
                        units="Pa",long_name="Full-level pressure at ps_ref")

        return

    end subroutine aeros_write_init

    subroutine aeros_write_state(io, filename, grd, now, time, n)
        ! Append one time slice of the instantaneous state.

        implicit none

        type(aeros_io_class),    intent(in) :: io
        character(len=*),        intent(in) :: filename
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_state_class), intent(in) :: now
        real(wp), intent(in) :: time
        integer,  intent(in) :: n      ! time index, 1-based

        integer :: ncid, q

        call nc_open(filename,ncid,writable=.TRUE.)

        call nc_write(filename,"time",time,dim1="time", &
                        start=[n],count=[1],grid_mapping="",ncid=ncid)

        do q = 1, size(io%state)
            call aeros_write_var_state(filename, io%state(q), grd, now, n, ncid)
        end do

        call nc_close(ncid)

        return

    end subroutine aeros_write_state

    subroutine aeros_write_var_state(filename, v, grd, now, n, ncid)

        implicit none

        character(len=*),        intent(in) :: filename
        type(var_io_type),       intent(in) :: v
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_state_class), intent(in) :: now
        integer, intent(in) :: n
        integer, intent(in), optional :: ncid

        character(len=32), allocatable :: dims(:)

        ! Local copy of the dims with "time" appended, so the table itself does
        ! not have to carry a dimension every state variable shares.
        allocate(dims(v%ndims+1))
        dims(1:v%ndims) = v%dims
        dims(v%ndims+1) = "time"

        select case(trim(v%varname))

            case("ps")
                call nc_write(filename,trim(v%varname),now%ps, &
                            start=[1,1,n],count=[grd%nlon,grd%nlat,1], &
                            units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("u")
                call nc_write(filename,trim(v%varname),now%u, &
                            start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                            units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("v")
                call nc_write(filename,trim(v%varname),now%v, &
                            start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                            units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("temp")
                call nc_write(filename,trim(v%varname),now%temp_g, &
                            start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                            units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)
            case("qv")
                call nc_write(filename,trim(v%varname),now%qv_g, &
                            start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                            units=v%units,long_name=v%long_name,dims=dims,ncid=ncid)

            case DEFAULT
                write(io_unit_err,*) "aeros_write_var_state:: error: variable is in the table "// &
                                        "but has no case here."
                write(io_unit_err,*) "variable = ", trim(v%varname)
                write(io_unit_err,*) "table    = input/aeros-variables-state.md"
                error stop 1

        end select

        return

    end subroutine aeros_write_var_state

    subroutine aeros_write_diag(io, filename, grd, vg, dgn)
        ! Write the accumulated zonal-mean statistics to their own file.
        !
        ! One file per averaging window, with no time dimension: the window is
        ! the file. Its bounds go in as global attributes so a result can never
        ! be separated from the period it averages.

        implicit none

        type(aeros_io_class),    intent(in) :: io
        character(len=*),        intent(in) :: filename
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_vgrid_class), intent(in) :: vg
        type(aeros_diag_class),  intent(in) :: dgn

        integer :: q

        if (.not. dgn%finalized) then
            write(io_unit_err,*) "aeros_write_diag:: error: the average has not been finalized. "// &
                                    "Call aeros_diag_finalize before writing."
            error stop 1
        end if

        call nc_create(filename)

        call nc_write_dim(filename,"lat",x=grd%lat,units="degrees_north")
        call nc_write_dim(filename,"lev",x=1,dx=1,nx=vg%nlev,units="1")

        call nc_write_attr(filename,"time_start",real(dgn%time_start,wp))
        call nc_write_attr(filename,"time_end",  real(dgn%time_end,wp))
        call nc_write_attr(filename,"nsample",   dgn%nsample)

        call nc_write(filename,"sigma",vg%sigma_full,dim1="lev", &
                        units="1",long_name="Sigma at full levels, at ps_ref")

        do q = 1, size(io%diag)
            call aeros_write_var_diag(filename, io%diag(q), dgn)
        end do

        return

    end subroutine aeros_write_diag

    subroutine aeros_write_var_diag(filename, v, dgn)

        implicit none

        character(len=*),       intent(in) :: filename
        type(var_io_type),      intent(in) :: v
        type(aeros_diag_class), intent(in) :: dgn

        select case(trim(v%varname))

            case("u_zm")
                call nc_write(filename,trim(v%varname),dgn%u, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("v_zm")
                call nc_write(filename,trim(v%varname),dgn%v, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("temp_zm")
                call nc_write(filename,trim(v%varname),dgn%temp, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("uv_eddy")
                call nc_write(filename,trim(v%varname),dgn%uv, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("vt_eddy")
                call nc_write(filename,trim(v%varname),dgn%vt, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("eke")
                call nc_write(filename,trim(v%varname),dgn%eke, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("temp_var")
                call nc_write(filename,trim(v%varname),dgn%tvar, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("pressure")
                call nc_write(filename,trim(v%varname),dgn%pressure, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("ps_zm")
                call nc_write(filename,trim(v%varname),dgn%ps, &
                            units=v%units,long_name=v%long_name,dims=v%dims)
            case("u_sfc")
                call nc_write(filename,trim(v%varname),dgn%u_sfc, &
                            units=v%units,long_name=v%long_name,dims=v%dims)

            case DEFAULT
                write(io_unit_err,*) "aeros_write_var_diag:: error: variable is in the table "// &
                                        "but has no case here."
                write(io_unit_err,*) "variable = ", trim(v%varname)
                write(io_unit_err,*) "table    = input/aeros-variables-diag.md"
                error stop 1

        end select

        return

    end subroutine aeros_write_var_diag

end module aeros_io
