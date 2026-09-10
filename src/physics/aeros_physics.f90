module aeros_physics
    ! Master on/off flags for the physics subsystems.
    !
    ! Wave 1 of the coupling API (docs/coupling_api.md sec. 5, 8) reserved the
    ! single `&physics` namelist group; Wave 4a.1 makes it the ONE source of
    ! truth for the on/off state of every subsystem on the library path:
    ! aeros_init loads it right after the dynamics parameters and
    ! aeros_timestep_init hands each flag to the matching subsystem `_load`
    ! (which no longer reads an `enabled` key of its own). rce_long still
    ! sets the `ts%<sub>%enabled` fields by hand from `&rce`; Wave 4a.2 moves
    ! it onto this group.
    !
    ! Type defaults here are the real-Earth production configuration; the
    ! shipped input/aeros_defaults.nml overrides them to the bare dynamical
    ! core (every switch off) so aeros_run.x / held_suarez.x keep their
    ! pre-Wave-4a behaviour -- a case file turns things on. `l_topography` has
    ! no library consumer yet (the phis provider is Wave 4b) and `l_moist` is
    ! unassigned (condensation is `l_cnd`; the transport is always on).

    use nml, only : nml_read

    implicit none

    private

    type, public :: aeros_physics_class
        logical :: l_surf        = .true.
        logical :: l_cnv         = .true.
        logical :: l_cnd         = .true.
        logical :: l_rad         = .true.
        logical :: l_land        = .true.
        logical :: l_topography  = .true.
        logical :: l_seaice      = .false.
        logical :: l_cloud_prog  = .false.
        logical :: l_bd          = .false.
        logical :: l_vdiff       = .true.
        logical :: l_moist       = .true.
    end type aeros_physics_class

    public :: aeros_physics_load

contains

    subroutine aeros_physics_load(phys, filename, defaults_file)
        ! Populate the master physics flags from the `&physics` namelist group.
        ! Every field of the class is namelist-readable; anything omitted keeps
        ! the current (initial) value on the passed struct.

        implicit none

        type(aeros_physics_class), intent(inout) :: phys
        character(len=*),          intent(in)    :: filename
        character(len=*), intent(in), optional   :: defaults_file

        call nml_read(filename, "physics", "l_surf",       phys%l_surf,       defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_cnv",        phys%l_cnv,        defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_cnd",        phys%l_cnd,        defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_rad",        phys%l_rad,        defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_land",       phys%l_land,       defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_topography", phys%l_topography, defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_seaice",     phys%l_seaice,     defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_cloud_prog", phys%l_cloud_prog, defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_bd",         phys%l_bd,         defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_vdiff",      phys%l_vdiff,      defaults_file=defaults_file)
        call nml_read(filename, "physics", "l_moist",      phys%l_moist,      defaults_file=defaults_file)

        return

    end subroutine aeros_physics_load

end module aeros_physics
