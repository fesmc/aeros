module aeros_state
    ! Allocation and initialization of the prognostic state.
    !
    ! M0 allocates and zeroes; nothing here evolves anything. The dynamical
    ! core arrives at M1 (docs/design.md section 7).

    use aeros_defs, only : wp, wp_sh, aeros_state_class, aeros_grid_class

    implicit none

    private

    public :: aeros_state_alloc
    public :: aeros_state_end

contains

    subroutine aeros_state_alloc(now, grd, nlm, nlev)
        ! Allocate the spectral and grid-space state and zero it.

        implicit none

        type(aeros_state_class), intent(inout) :: now
        type(aeros_grid_class),  intent(in)    :: grd
        integer, intent(in) :: nlm    ! spectral coefficients per level
        integer, intent(in) :: nlev   ! vertical levels

        call aeros_state_end(now)

        now%nlm  = nlm
        now%nlev = nlev

        allocate(now%vor(nlm,nlev))
        allocate(now%div(nlm,nlev))
        allocate(now%temp(nlm,nlev))
        allocate(now%qv(nlm,nlev))
        allocate(now%lnps(nlm))

        allocate(now%u(grd%nlon,grd%nlat,nlev))
        allocate(now%v(grd%nlon,grd%nlat,nlev))
        allocate(now%temp_g(grd%nlon,grd%nlat,nlev))
        allocate(now%qv_g(grd%nlon,grd%nlat,nlev))
        allocate(now%ps(grd%nlon,grd%nlat))

        now%vor  = (0.0_wp_sh, 0.0_wp_sh)
        now%div  = (0.0_wp_sh, 0.0_wp_sh)
        now%temp = (0.0_wp_sh, 0.0_wp_sh)
        now%qv   = (0.0_wp_sh, 0.0_wp_sh)
        now%lnps = (0.0_wp_sh, 0.0_wp_sh)

        now%u      = 0.0_wp
        now%v      = 0.0_wp
        now%temp_g = 0.0_wp
        now%qv_g   = 0.0_wp
        now%ps     = 0.0_wp

        return

    end subroutine aeros_state_alloc

    subroutine aeros_state_end(now)

        implicit none

        type(aeros_state_class), intent(inout) :: now

        now%nlm = 0; now%nlev = 0

        if (allocated(now%vor))    deallocate(now%vor)
        if (allocated(now%div))    deallocate(now%div)
        if (allocated(now%temp))   deallocate(now%temp)
        if (allocated(now%qv))     deallocate(now%qv)
        if (allocated(now%lnps))   deallocate(now%lnps)

        if (allocated(now%u))      deallocate(now%u)
        if (allocated(now%v))      deallocate(now%v)
        if (allocated(now%temp_g)) deallocate(now%temp_g)
        if (allocated(now%qv_g))   deallocate(now%qv_g)
        if (allocated(now%ps))     deallocate(now%ps)

        return

    end subroutine aeros_state_end

end module aeros_state
