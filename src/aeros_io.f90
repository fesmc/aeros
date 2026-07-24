module aeros_io
    ! netCDF output for the global core, via fesm-utils' ncio.
    !
    ! Deliberately takes the grid and state directly rather than the top-level
    ! aeros_class: output has no business knowing about the SHTns handle, and
    ! keeping it below the facade in the dependency order means the facade can
    ! use it rather than the other way round.
    !
    ! Variable selection is hard-coded here at M0. It moves to fesm-utils'
    ! variable_io table (docs/design.md section 7, M0) once there are enough
    ! fields for a table to be worth more than a list.

    use aeros_defs, only : sp, dp, wp, aeros_grid_class, aeros_state_class
    use ncio

    implicit none

    private

    public :: aeros_write_init
    public :: aeros_write_step

contains

    subroutine aeros_write_init(filename, grd, nlev, time_init, units_time)
        ! Create the output file and write its dimensions.

        implicit none

        character(len=*), intent(in) :: filename
        type(aeros_grid_class), intent(in) :: grd
        integer,  intent(in) :: nlev
        real(wp), intent(in) :: time_init
        character(len=*), intent(in), optional :: units_time

        character(len=56) :: tunits

        tunits = "years"
        if (present(units_time)) tunits = trim(units_time)

        call nc_create(filename)

        call nc_write_dim(filename,"lon",x=grd%lon,units="degrees_east")
        call nc_write_dim(filename,"lat",x=grd%lat,units="degrees_north")
        call nc_write_dim(filename,"lev",x=1,dx=1,nx=nlev,units="1")
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

        return

    end subroutine aeros_write_init

    subroutine aeros_write_step(filename, grd, now, time, n)
        ! Append one time slice.

        implicit none

        character(len=*), intent(in) :: filename
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_state_class), intent(in) :: now
        real(wp), intent(in) :: time
        integer,  intent(in) :: n      ! time index, 1-based

        integer :: ncid

        call nc_open(filename,ncid,writable=.TRUE.)

        call nc_write(filename,"time",time,dim1="time", &
                        start=[n],count=[1],grid_mapping="",ncid=ncid)

        call nc_write(filename,"ps",now%ps,dim1="lon",dim2="lat",dim3="time", &
                        start=[1,1,n],count=[grd%nlon,grd%nlat,1], &
                        units="Pa",long_name="Surface pressure",ncid=ncid)

        call nc_write(filename,"u",now%u,dim1="lon",dim2="lat",dim3="lev",dim4="time", &
                        start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                        units="m s-1",long_name="Zonal wind",ncid=ncid)

        call nc_write(filename,"v",now%v,dim1="lon",dim2="lat",dim3="lev",dim4="time", &
                        start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                        units="m s-1",long_name="Meridional wind",ncid=ncid)

        call nc_write(filename,"temp",now%temp_g,dim1="lon",dim2="lat",dim3="lev",dim4="time", &
                        start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                        units="K",long_name="Temperature",ncid=ncid)

        call nc_write(filename,"qv",now%qv_g,dim1="lon",dim2="lat",dim3="lev",dim4="time", &
                        start=[1,1,1,n],count=[grd%nlon,grd%nlat,now%nlev,1], &
                        units="kg kg-1",long_name="Specific humidity",ncid=ncid)

        call nc_close(ncid)

        return

    end subroutine aeros_write_step

end module aeros_io
