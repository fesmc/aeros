module aeros_state
    ! Allocation and lifetime of the prognostic state.
    !
    ! Two levels, matching the two types in aeros_defs: `aeros_spec_class` is
    ! the spectral prognostics alone, which is what the time integrator carries
    ! three copies of; `aeros_state_class` is one time level complete with its
    ! grid-space diagnostics, of which the model needs exactly one.

    use aeros_defs, only : wp, wp_sh, io_unit_err, &
                            aeros_spec_class, aeros_state_class, aeros_grid_class

    implicit none

    private

    public :: aeros_spec_alloc
    public :: aeros_spec_end
    public :: aeros_spec_zero
    public :: aeros_spec_copy
    public :: aeros_spec_swap
    public :: aeros_state_alloc
    public :: aeros_state_end

contains

    subroutine aeros_spec_alloc(spec, nlm, nlev)
        ! Allocate the spectral prognostics and zero them.

        implicit none

        type(aeros_spec_class), intent(inout) :: spec
        integer, intent(in) :: nlm    ! spectral coefficients per level
        integer, intent(in) :: nlev   ! vertical levels

        call aeros_spec_end(spec)

        spec%nlm  = nlm
        spec%nlev = nlev

        allocate(spec%vor(nlm,nlev))
        allocate(spec%div(nlm,nlev))
        allocate(spec%temp(nlm,nlev))
        allocate(spec%lnps(nlm))

        call aeros_spec_zero(spec)

        return

    end subroutine aeros_spec_alloc

    subroutine aeros_spec_end(spec)

        implicit none

        type(aeros_spec_class), intent(inout) :: spec

        spec%nlm = 0; spec%nlev = 0

        if (allocated(spec%vor))  deallocate(spec%vor)
        if (allocated(spec%div))  deallocate(spec%div)
        if (allocated(spec%temp)) deallocate(spec%temp)
        if (allocated(spec%lnps)) deallocate(spec%lnps)

        return

    end subroutine aeros_spec_end

    subroutine aeros_spec_zero(spec)

        implicit none

        type(aeros_spec_class), intent(inout) :: spec

        spec%vor  = (0.0_wp_sh, 0.0_wp_sh)
        spec%div  = (0.0_wp_sh, 0.0_wp_sh)
        spec%temp = (0.0_wp_sh, 0.0_wp_sh)
        spec%lnps = (0.0_wp_sh, 0.0_wp_sh)

        return

    end subroutine aeros_spec_zero

    subroutine aeros_spec_copy(dst, src)
        ! Deep copy, dst <- src. Both must already be allocated to the same
        ! shape: this is the assignment used inside a timestep, where the
        ! allocation is done once at init, so an allocate-on-assignment here
        ! would hide a shape error rather than report one.

        implicit none

        type(aeros_spec_class), intent(inout) :: dst
        type(aeros_spec_class), intent(in)    :: src

        if (dst%nlm /= src%nlm .or. dst%nlev /= src%nlev) then
            write(io_unit_err,*) "aeros_spec_copy:: error: shape mismatch, ", &
                                    dst%nlm, dst%nlev, " vs ", src%nlm, src%nlev
            error stop 1
        end if

        dst%vor  = src%vor
        dst%div  = src%div
        dst%temp = src%temp
        dst%lnps = src%lnps

        return

    end subroutine aeros_spec_copy

    subroutine aeros_spec_swap(a, b)
        ! Exchange the CONTENTS of two spectral states without copying them.
        !
        ! This is how the leapfrog shift (old <- now <- new) is done: two swaps
        ! rotate the three time levels and leave the retired one as scratch for
        ! the next step. move_alloc moves the descriptor, so the cost is
        ! independent of the array size and no temporary is ever materialized.

        implicit none

        type(aeros_spec_class), intent(inout) :: a, b

        complex(wp_sh), allocatable :: tmp2(:,:), tmp1(:)
        integer :: itmp

        itmp = a%nlm;  a%nlm  = b%nlm;  b%nlm  = itmp
        itmp = a%nlev; a%nlev = b%nlev; b%nlev = itmp

        call move_alloc(a%vor,  tmp2); call move_alloc(b%vor,  a%vor);  call move_alloc(tmp2, b%vor)
        call move_alloc(a%div,  tmp2); call move_alloc(b%div,  a%div);  call move_alloc(tmp2, b%div)
        call move_alloc(a%temp, tmp2); call move_alloc(b%temp, a%temp); call move_alloc(tmp2, b%temp)
        call move_alloc(a%lnps, tmp1); call move_alloc(b%lnps, a%lnps); call move_alloc(tmp1, b%lnps)

        return

    end subroutine aeros_spec_swap

    subroutine aeros_state_alloc(now, grd, nlm, nlev)
        ! Allocate one complete time level -- spectral plus grid -- and zero it.

        implicit none

        type(aeros_state_class), intent(inout) :: now
        type(aeros_grid_class),  intent(in)    :: grd
        integer, intent(in) :: nlm    ! spectral coefficients per level
        integer, intent(in) :: nlev   ! vertical levels

        call aeros_state_end(now)

        now%nlev = nlev

        call aeros_spec_alloc(now%spec, nlm, nlev)

        allocate(now%u(grd%nlon,grd%nlat,nlev))
        allocate(now%v(grd%nlon,grd%nlat,nlev))
        allocate(now%temp_g(grd%nlon,grd%nlat,nlev))
        allocate(now%qv_g(grd%nlon,grd%nlat,nlev))
        allocate(now%cf_g(grd%nlon,grd%nlat,nlev))
        allocate(now%ps(grd%nlon,grd%nlat))

        now%u      = 0.0_wp
        now%v      = 0.0_wp
        now%temp_g = 0.0_wp
        now%qv_g   = 0.0_wp
        now%cf_g   = 0.0_wp
        now%ps     = 0.0_wp

        return

    end subroutine aeros_state_alloc

    subroutine aeros_state_end(now)

        implicit none

        type(aeros_state_class), intent(inout) :: now

        now%nlev = 0

        call aeros_spec_end(now%spec)

        if (allocated(now%u))      deallocate(now%u)
        if (allocated(now%v))      deallocate(now%v)
        if (allocated(now%temp_g)) deallocate(now%temp_g)
        if (allocated(now%qv_g))   deallocate(now%qv_g)
        if (allocated(now%cf_g))   deallocate(now%cf_g)
        if (allocated(now%ps))     deallocate(now%ps)

        return

    end subroutine aeros_state_end

end module aeros_state
