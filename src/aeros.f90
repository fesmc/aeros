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
                                aeros_vgrid_print, aeros_vgrid_pressure
    use aeros_physics,  only : aeros_physics_class, aeros_physics_load
    use aeros_condensation, only : aeros_qsat
    use aeros_timestep, only : aeros_timestep_class, aeros_timestep_init, &
                                aeros_timestep_end, aeros_timestep_step, &
                                aeros_timestep_diagnose, aeros_timestep_print, &
                                aeros_timestep_write_restart, &
                                aeros_timestep_read_restart, aeros_timestep_set_phis
    use aeros_boundary, only : aeros_boundary_t, aeros_bc_2d_t, aeros_boundary_load, &
                                aeros_boundary_end, aeros_boundary_refresh, &
                                BC_NONE, BC_EXTERNAL
    use aeros_ocean,    only : OCEAN_PRESCRIBED
    use aeros_land,     only : aeros_land_is_land, aeros_land_couple_radiation
    use aeros_topography, only : aeros_topography_bandlimit
    use ncio,           only : nc_write, nc_read, nc_exists_var
    use aeros_diag_acc, only : aeros_diag_acc_t, aeros_diag_acc_init, &
                                aeros_diag_acc_reset, aeros_diag_acc_end
    use aeros_bundles,  only : aeros_imports_t, aeros_exports_t, aeros_diag_t

    implicit none

    private

    type aeros_class
        type(aeros_param_class)    :: par
        ! Master physics on/off flags (&physics). The one source of every
        ! subsystem's enabled state on the library path (Wave 4a): loaded
        ! right after `par` and handed to the subsystem loaders by
        ! aeros_timestep_init.
        type(aeros_physics_class)  :: phys
        type(aeros_sht_pool_class) :: pool
        type(aeros_grid_class)     :: grd
        type(aeros_vgrid_class)    :: vgrid
        type(aeros_state_class)    :: now
        type(aeros_timestep_class) :: ts

        ! Boundary providers (docs/coupling_api.md sec 7): the one source of
        ! every externally prescribed field. Populated by aeros_boundary_load
        ! in aeros_init; every 2D provider's `buffer` is allocated on the
        ! model grid. The physics never reads a provider directly:
        ! aeros_boundary_apply copies each refreshed buffer into the working
        ! field its kernel reads (ts%ocn%sst, ts%land%mask / %albedo,
        ! ts%rad%dust_ot / %so4_map, the dynamics' phis) -- at the end of
        ! aeros_init, after every refresh in aeros_lib_step, and wherever a
        ! driver refreshes the providers itself (rce_long, daily).
        type(aeros_boundary_t) :: bounds
        ! Provider version of the phis buffer last band-limited and handed to
        ! the dynamics (aeros_boundary_apply): the analysis/synthesis
        ! round-trip is done once per change, never repeated on an unchanged
        ! buffer (a second round-trip is not bit-identical).
        integer :: bc_phis_applied = 0

        ! Library-side diagnostic accumulators (docs/coupling_api.md sec. 5,
        ! 6.2, 6.3). Wave 4a: aeros_init allocates and zeroes the always-on
        ! buffers (aeros_diag_acc_init), so the container is usable after init;
        ! the monthly climatology bins stay a driver choice
        ! (aeros_diag_acc_init_monthly). rce_long accumulates into this
        ! container (Wave 4a.2); Wave 4c moves the accumulation itself behind
        ! library methods.
        type(aeros_diag_acc_t) :: diag_acc

        real(wp) :: time                  ! current model time [yr]
        real(wp) :: time_init             ! time at aeros_init [yr]
        integer  :: nstep = 0             ! dynamics steps taken since init
    end type aeros_class

    public :: aeros_class
    public :: aeros_init
    public :: aeros_update
    public :: aeros_end
    public :: aeros_set_isothermal
    public :: aeros_set_rce_warm
    public :: aeros_diagnose
    public :: aeros_print_config
    public :: aeros_boundary_apply

    ! Coupled library entry points (docs/coupling_api.md sec. 3, 12, 14).
    ! Wave-3.a addition: public interface added, but no in-tree driver calls
    ! them yet -- Wave 4 migrates rce_long onto them. Because these are
    ! unused-code paths on the Phase A regression set, bit-for-bit is
    ! trivially preserved.
    public :: aeros_lib_init
    public :: aeros_lib_step
    public :: aeros_lib_finalize
    public :: aeros_lib_restart_write
    public :: aeros_lib_restart_read
    public :: aeros_lib_diagnose

contains

    subroutine aeros_init(ams, filename, group, time_init, defaults_file)
        ! Configure the transform pool, the grids, the state and the integrator.
        !
        ! `defaults_file`, when given, is a namelist holding the full parameter
        ! schema: every read takes its value from there and `filename` need only
        ! carry the overrides for this case (see nml_read's defaults_file). Absent,
        ! the legacy behaviour holds and `filename` must be complete.

        implicit none

        type(aeros_class), intent(inout) :: ams
        character(len=*),  intent(in)    :: filename
        character(len=*),  intent(in), optional :: group
        real(wp),          intent(in), optional :: time_init
        character(len=*),  intent(in), optional :: defaults_file

        call aeros_par_load(ams%par, filename, group, defaults_file=defaults_file)

        ! Master physics flags (&physics; Wave 4a), loaded right after the
        ! dynamics parameters: aeros_timestep_init hands them to the subsystem
        ! loaders, which take their on/off state from nowhere else.
        call aeros_physics_load(ams%phys, filename, defaults_file=defaults_file)

        ! One SHTns config per thread; see aeros_sht_pool_class for why aeros
        ! threads the level loop rather than letting SHTns thread a transform.
        ! sht_mode = "quick" (the schema default since Wave 4a.2) skips the
        ! auto-tuning altogether and is the reproducible choice; "cache" tunes
        ! once per working directory and saves/loads the tuning -- the pool
        ! builds nthreads IDENTICAL configs and SHTns does not remember its
        ! tuning between them, so without the cache a 10-thread pool would pay
        ! it ten times. The two are not bit-identical (see
        ! aeros_param_class%sht_mode).
        select case (trim(ams%par%sht_mode))
        case ("cache")
            call aeros_sht_pool_init(ams%pool, ams%par%trunc, &
                                        nlon=ams%par%nlon, nlat=ams%par%nlat, &
                                        nthreads=ams%par%nthreads, cache=.TRUE.)
        case ("quick")
            call aeros_sht_pool_init(ams%pool, ams%par%trunc, &
                                        nlon=ams%par%nlon, nlat=ams%par%nlat, &
                                        nthreads=ams%par%nthreads, quick=.TRUE.)
        case default
            write(io_unit_err,*) "aeros_init:: error: sht_mode must be 'cache' or 'quick', got '"// &
                                    trim(ams%par%sht_mode)//"'"
            stop 1
        end select

        call aeros_grid_init(ams%grd, ams%pool%sht(1))

        ! Non-rotating planet (l_nonrotating): with no Coriolis there is
        ! nothing to organise a meridional asymmetry into a jet, so with
        ! uniform forcing every column runs the same RCE -- the single-column
        ! validation vehicle on the full physics stack. Zeroed here, between
        ! the grid and the integrator init, which is where everything that
        ! reads grd%coriolis is built.
        if (ams%par%l_nonrotating) ams%grd%coriolis = 0.0_wp

        call aeros_vgrid_load(ams%vgrid, ams%par%nlev, filename, defaults_file=defaults_file)

        call aeros_state_alloc(ams%now, ams%grd, ams%pool%sht(1)%nlm, ams%par%nlev)

        call aeros_timestep_init(ams%ts, ams%par, ams%pool, ams%grd, ams%vgrid, filename, &
                                    defaults_file=defaults_file, phys=ams%phys)

        ! Boundary providers: load the &boundary group and read every active
        ! BC_FILE field onto the model grid. A provider is active only when a
        ! consumer of it is on -- the &physics flags are the one source of
        ! that, so a run with land off never opens the mask file, a dry core
        ! never opens the SST table, and the schema can name the calibrated
        ! real-Earth files for every field without making every run depend
        ! on them. (l_topography off = flat, whatever the phis provider says;
        ! a slab ocean owns its SST, so the sst provider is only for &ocean
        ! mode = 0 and only when something reads the surface temperature.)
        ams%bounds%lsm%active    = ams%phys%l_land
        ams%bounds%albedo%active = ams%phys%l_land
        ams%bounds%dust%active   = ams%phys%l_rad
        ams%bounds%so4%active    = ams%phys%l_rad
        ams%bounds%phis%active   = ams%phys%l_topography
        ams%bounds%sst%active    = ams%ts%ocn%mode == OCEAN_PRESCRIBED .and. &
                                   (ams%phys%l_surf .or. ams%phys%l_rad)
        call aeros_boundary_load(ams%bounds, filename, ams%grd, defaults_file=defaults_file)
        ! Seed the monthly fields at the start day and hand every buffer to
        ! its consumer: the model is complete on return from init, with the
        ! SST / albedo / dust of doy0 (the values rce_long's set_seasonal_sst,
        ! aeros_land_set_albedo and aeros_radiation_set_dust installed at
        ! the start before Wave 4b) and the band-limited orography.
        call aeros_boundary_refresh(ams%bounds, ams%ts%rad%doy0)
        call aeros_boundary_apply(ams)

        ! Library-side diagnostic accumulators (Wave 4a): allocate and zero the
        ! always-on buffers so ams%diag_acc is usable by a driver and by
        ! aeros_lib_step, as rce_long's local dacc is after its own init call.
        ! The twelve-month climatology bins are NOT allocated here -- they are
        ! an output choice (rce_long's monthly_out), so a driver that wants
        ! them calls aeros_diag_acc_init_monthly(ams%diag_acc, ...) after init
        ! and nothing else pays for their storage.
        call aeros_diag_acc_end(ams%diag_acc)
        call aeros_diag_acc_init(ams%diag_acc, ams%grd%nlon, ams%grd%nlat, ams%par%nlev, &
                                    mmax_spec=ams%par%trunc, dt=ams%par%dt)

        ! Initial condition. A resting isothermal atmosphere over a flat
        ! surface -- the only state the dry core can be started from at M1
        ! without inventing a climatology, and the one whose correct behaviour
        ! (staying exactly where it is) the acceptance tests measure. It is
        ! also the state Held-Suarez is spun up from at M1.5, perturbed.
        !
        ! NOT the zeroed state aeros_state_alloc leaves: ln(p_s) = 0 there
        ! means a surface pressure of 1 Pa, which is not a small error but a
        ! different planet.
        call aeros_set_isothermal(ams)

        ams%time_init = 0.0_wp
        if (present(time_init)) ams%time_init = time_init
        ams%time  = ams%time_init
        ams%nstep = 0

        return

    end subroutine aeros_init

    subroutine aeros_set_isothermal(ams, temp, ps)
        ! Reset the state to a motionless isothermal atmosphere.
        !
        ! `temp` and `ps` default to the vertical grid's reference state, which
        ! is the initial condition aeros_init installs. A benchmark that wants
        ! to start somewhere else -- Held-Suarez conventionally starts
        ! isothermal well below its own equilibrium, so the 40-day relaxation
        ! has several e-folds to work over the 200-day spin-up -- passes its
        ! own.
        !
        ! Ends by refreshing the grid-space fields, so the state is complete on
        ! return and a driver can perturb it and integrate without knowing that
        ! the two halves exist.

        implicit none

        ! `target` so the pool's configs can be pointed at: gfortran refuses a
        ! pointer assignment to a component of a non-target dummy. Harmless
        ! here -- the pointer does not outlive the call.
        type(aeros_class), intent(inout), target :: ams
        real(wp), intent(in), optional :: temp   ! [K]
        real(wp), intent(in), optional :: ps     ! [Pa]

        type(aeros_sht_class), pointer :: s
        real(wp) :: t_, ps_
        real(dp) :: y00
        integer  :: k, lm00

        s => ams%pool%sht(1)

        ps_ = ams%vgrid%ps_ref

        call aeros_spec_zero(ams%now%spec)

        ! Orthonormal harmonics with no Condon-Shortley phase (aeros_spectral):
        ! Y_00 = 1/sqrt(4 pi), so a uniform field of value X is the single
        ! coefficient X*sqrt(4 pi).
        y00  = sqrt(4.0_dp*acos(-1.0_dp))
        lm00 = aeros_sht_lm(s, 0, 0)

        if (present(ps)) ps_ = ps

        do k = 1, ams%vgrid%nlev
            t_ = ams%vgrid%t_ref(k)
            if (present(temp)) t_ = temp
            ams%now%spec%temp(lm00,k) = cmplx(real(t_,dp)*y00, 0.0_dp, wp_sh)
        end do
        ams%now%spec%lnps(lm00) = cmplx(log(real(ps_,dp))*y00, 0.0_dp, wp_sh)

        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)

        return

    end subroutine aeros_set_isothermal

    subroutine aeros_set_rce_warm(ams, seed_asym, gradient)
        ! Reset the state to the warm, humid, conditionally unstable cold start
        ! of the RCE runs (drivers/rce_long.f90's IC block, moved here so the
        ! library path reproduces it bit for bit; also test_rce/test_moist_run):
        !
        !     T(k)   = 300 - 90 (1 - ((k - 1/2)/nlev)^0.6)     [K], resting
        !     ln p_s = ln(p0)
        !     T(l=2, m=0) = -5 (spectral units): the initial equator-pole
        !                   gradient -- skipped when `gradient` is .FALSE.,
        !                   the uniform-insolation single-column case, where
        !                   the intended state is horizontally uniform
        !     T(l=m+2, m=1..6) = seed_asym (1+i) at every level: an optional
        !                   zonal-asymmetry seed so baroclinic instability can
        !                   grow eddies. Without it the run is trapped in the
        !                   axisymmetric manifold -- with zonally symmetric
        !                   forcing and an m=0-only start the nonlinear terms
        !                   never populate m>0.
        !     qv     = 0.9 qsat(T, p) on the grid
        !
        ! `seed_asym` defaults to the &aeros key (par%seed_asym); `gradient`
        ! to .TRUE. Ends with the grid-space fields refreshed (they are what
        ! the humidity is built from), so the state is complete on return.

        implicit none

        type(aeros_class), intent(inout), target :: ams
        real(wp), intent(in), optional :: seed_asym
        logical,  intent(in), optional :: gradient

        type(aeros_sht_class), pointer :: s
        real(wp) :: tval, seed, qs, dqsdt
        real(wp), allocatable :: phalf(:), pfull(:), dpc(:)
        integer  :: i, j, k, m, nlev, lm00
        logical  :: grad

        s    => ams%pool%sht(1)
        nlev = ams%vgrid%nlev
        lm00 = aeros_sht_lm(s, 0, 0)

        seed = ams%par%seed_asym
        if (present(seed_asym)) seed = seed_asym
        grad = .TRUE.
        if (present(gradient)) grad = gradient

        call aeros_spec_zero(ams%now%spec)
        do k = 1, nlev
            tval = 300.0_wp - 90.0_wp*(1.0_wp - ((real(k,wp) - 0.5_wp)/real(nlev,wp))**0.6_wp)
            ams%now%spec%temp(lm00,k) = &
                    cmplx(real(tval,dp)*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
        end do
        ams%now%spec%lnps(lm00) = &
                cmplx(log(real(p0,dp))*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)

        if (grad) ams%now%spec%temp(aeros_sht_lm(s,2,0),:) = cmplx(-5.0_dp, 0.0_dp, wp_sh)

        if (seed > 0.0_wp) then
            do k = 1, nlev
                do m = 1, 6
                    ams%now%spec%temp(aeros_sht_lm(s, m+2, m), k) = &
                        cmplx(real(seed,dp), real(seed,dp), wp_sh)
                end do
            end do
        end if

        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)

        allocate(phalf(0:nlev), pfull(nlev), dpc(nlev))
        do j = 1, ams%grd%nlat
            do i = 1, ams%grd%nlon
                call aeros_vgrid_pressure(ams%vgrid, ams%now%ps(i,j), phalf, pfull, dpc)
                do k = 1, nlev
                    call aeros_qsat(ams%now%temp_g(i,j,k), pfull(k), qs, dqsdt)
                    ams%now%qv_g(i,j,k) = 0.9_wp*qs
                end do
            end do
        end do
        deallocate(phalf, pfull, dpc)

        return

    end subroutine aeros_set_rce_warm

    subroutine aeros_boundary_apply(ams)
        ! Hand every provider buffer to the physics field that reads it.
        !
        ! The kernels keep their working fields (ts%ocn%sst, ts%land%mask and
        ! %albedo, ts%rad%dust_ot and %so4_map, the tendency's phis) and are
        ! untouched by the providers; this is the one routine that connects
        ! the two, so a coupler that pushes any of these fields as
        ! BC_EXTERNAL through aeros_lib_step(imports=...) needs no physics
        ! change. Called at the end of aeros_init and after every
        ! aeros_boundary_refresh (aeros_lib_step; rce_long's daily refresh).
        !
        ! A provider is applied only when it is active (its consumer is on)
        ! and has a kind other than BC_NONE -- otherwise the consumer keeps
        ! its own default (the ocean's analytic SST, the land's scalar
        ! land_albedo, an all-ocean mask, no dust map, no sulfate map, a
        ! flat surface).
        !
        !   sst    -> ts%ocn%sst, prescribed-SST mode only (a slab owns it).
        !   lsm    -> ts%land%mask, thresholded at &land lsm_threshold.
        !   albedo -> ts%land%albedo (the bare land albedo); the per-point
        !             surface albedo map radiation reads (rad%alb_map) is
        !             recomposed whenever land is on, so a mask or albedo
        !             change reaches radiation here.
        !   dust   -> ts%rad%dust_ot (allocated on first use).
        !   so4    -> ts%rad%so4_map (allocated on first use).
        !   phis   -> spectrally band-limited to the model truncation IN
        !             PLACE (aeros_topography_bandlimit; the buffer then
        !             holds the orography the dynamics is consistent with,
        !             which is what a driver that ramps it reads) and set on
        !             the dynamics (aeros_timestep_set_phis) -- once per
        !             buffer version, since a repeated round-trip is not
        !             bit-identical. Only with &physics l_topography; off,
        !             the surface stays flat whatever the provider holds.
        !             A driver may still scale what the dynamics sees
        !             (rce_long's spin-up ramp); the buffer is the full field.

        implicit none

        type(aeros_class), intent(inout) :: ams

        if (in_use(ams%bounds%sst) .and. ams%ts%ocn%mode == OCEAN_PRESCRIBED) &
            ams%ts%ocn%sst = ams%bounds%sst%buffer

        if (ams%ts%land%enabled) then
            if (in_use(ams%bounds%lsm)) &
                ams%ts%land%mask = aeros_land_is_land(ams%bounds%lsm%buffer, ams%ts%land%lsm_threshold)
            if (in_use(ams%bounds%albedo)) ams%ts%land%albedo = ams%bounds%albedo%buffer
            call aeros_land_couple_radiation(ams%ts%land, ams%ts%rad%alb_map, ams%ts%rad%albedo)
        end if

        if (ams%ts%rad%enabled) then
            if (in_use(ams%bounds%dust)) then
                if (.not. allocated(ams%ts%rad%dust_ot)) &
                    allocate(ams%ts%rad%dust_ot(ams%grd%nlon, ams%grd%nlat))
                ams%ts%rad%dust_ot = ams%bounds%dust%buffer
            end if
            if (in_use(ams%bounds%so4)) then
                if (.not. allocated(ams%ts%rad%so4_map)) &
                    allocate(ams%ts%rad%so4_map(ams%grd%nlon, ams%grd%nlat))
                ams%ts%rad%so4_map = ams%bounds%so4%buffer
            end if
        end if

        if (ams%phys%l_topography .and. in_use(ams%bounds%phis)) then
            if (ams%bounds%phis%version /= ams%bc_phis_applied) then
                call aeros_topography_bandlimit(ams%pool%sht(1), ams%bounds%phis%buffer)
                call aeros_timestep_set_phis(ams%ts, ams%bounds%phis%buffer)
                ams%bc_phis_applied = ams%bounds%phis%version
            end if
        end if

        return

    contains

        pure logical function in_use(fld)
            type(aeros_bc_2d_t), intent(in) :: fld
            in_use = fld%active .and. fld%kind /= BC_NONE
        end function in_use

    end subroutine aeros_boundary_apply

    subroutine aeros_diagnose(ams)
        ! Refresh the grid-space fields from the spectral state.
        !
        ! aeros_update already ends with this, so a driver only needs it after
        ! reaching into the spectral state directly -- installing an initial
        ! condition, or seeding a perturbation.

        implicit none

        type(aeros_class), intent(inout) :: ams

        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)

        return

    end subroutine aeros_diagnose

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

        call aeros_diag_acc_end(ams%diag_acc)
        call aeros_boundary_end(ams%bounds)
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

    ! =========================================================================
    ! Coupled library API (docs/coupling_api.md sec. 3).
    !
    ! Six public routines a coupled host calls in place of the standalone
    ! aeros_init / aeros_update / aeros_end trio, with `ierr` returns instead
    ! of `stop 1` on the public path. The signatures track spec sec. 3 and
    ! decision 1 of sec. 14 (public routines live on the existing `aeros`
    ! module, not a new `aeros_lib`).
    !
    ! Wave-3.a scope is the interface itself: entry points added, bundle types
    ! shaped in aeros_bundles, and the internal routing (aeros_boundary_load /
    ! _refresh, aeros_diag_acc_reset, aeros_timestep_step,
    ! aeros_timestep_{write,read}_restart) wired. In-tree drivers do NOT call
    ! these yet, so the code below is unused on the Phase A regression set
    ! and bit-for-bit is trivial. Wave 4 migrates rce_long onto aeros_lib_*.
    !
    ! Note on step-count knob: spec sec. 3 gives aeros_lib_step a dt_couple
    ! argument in seconds; the sub-cycle count is `nint(dt_couple / par%dt)`.
    ! No new `nstep_per_cycle` namelist key is added -- the caller decides
    ! cadence at each call.
    ! =========================================================================

    subroutine aeros_lib_init(ams, cfg_file, defaults_file, ic_source, ierr)
        ! Configure the transform pool, grids, state, integrator and boundary
        ! providers. Does NOT start the time loop.
        !
        ! `ic_source` selects the initial condition (spec sec. 3):
        !   "resting" (default) -- motionless isothermal (today's facade).
        !   "restart"           -- reserved; call aeros_lib_restart_read after
        !                          init to load a checkpoint. init installs
        !                          "resting" so the shapes are ready.
        !   "rce_warm"          -- rce_long's warm, humid cold start
        !                          (aeros_set_rce_warm, seed from &aeros
        !                          seed_asym). Wave 4a; rce_long's cold start
        !                          since 4a.2.
        !   "coupler"           -- reserved for Phase C. Not accepted here.

        implicit none

        type(aeros_class), intent(inout)     :: ams
        character(len=*),  intent(in)        :: cfg_file
        character(len=*),  intent(in), optional :: defaults_file
        character(len=*),  intent(in), optional :: ic_source
        integer,           intent(out)       :: ierr

        character(len=32) :: ic

        ierr = 0

        ic = "resting"
        if (present(ic_source)) ic = ic_source

        select case (trim(ic))
        case ("resting", "restart", "rce_warm")
            ! "restart" defers the state overwrite to aeros_lib_restart_read;
            ! init still installs the resting shape so allocations exist.
            ! "rce_warm" is installed below, once the grids and state exist.
        case default
            write(io_unit_err,*) "aeros_lib_init:: unsupported ic_source='", &
                                    trim(ic), "' (Phase A: 'resting', 'restart' or 'rce_warm')"
            ierr = 1
            return
        end select

        call aeros_init(ams, cfg_file, defaults_file=defaults_file)

        if (trim(ic) == "rce_warm") call aeros_set_rce_warm(ams)

        return

    end subroutine aeros_lib_init

    subroutine aeros_lib_step(ams, dt_couple, imports, exports, ierr)
        ! Advance the model by dt_couple seconds by sub-cycling
        ! aeros_timestep_step nint(dt_couple / par%dt) times.
        !
        ! Consumes `imports` at the start of the cycle (spec sec. 7 timing
        ! steps 1-4: copy into BC_EXTERNAL buffers, refresh BC_FILE monthly
        ! buffers by day-of-year, then physics sees a fresh buffer set) and
        ! populates `exports` at the end. Both bundles are optional -- an
        ! absent one is treated as "not exchanging that side this cycle",
        ! preserving today's standalone behaviour.
        !
        ! The imports copy feeds the BC_EXTERNAL buffers of every provider
        ! (sst, lsm, albedo, dust, so4, phis have consumers; the rest are
        ! Phase C) and aeros_boundary_apply hands the refreshed set to the
        ! physics. exports is still a shell (Wave 4c).

        implicit none

        type(aeros_class),     intent(inout)          :: ams
        real(dp),              intent(in)             :: dt_couple
        type(aeros_imports_t), intent(in),  optional  :: imports
        type(aeros_exports_t), intent(out), optional  :: exports
        integer,               intent(out)            :: ierr

        integer  :: nsub, n
        real(wp) :: doy

        ierr = 0

        if (dt_couple <= 0.0_dp) then
            write(io_unit_err,*) "aeros_lib_step:: dt_couple must be > 0, got ", dt_couple
            ierr = 1
            return
        end if

        nsub = nint(dt_couple/real(ams%par%dt, dp))
        if (nsub < 1) nsub = 1

        ! Day-of-year for the boundary refresh: the model's one calendar,
        ! &radiation doy0 plus the elapsed steps -- the same expression the
        ! physics uses inside the step (insolation, Brewer-Dobson) and the
        ! one rce_long passes to its daily refresh, evaluated for the step
        ! about to run (rce_long refreshes before step n with n*dt). The
        ! 365-day wrap is aeros_boundary_interp_monthly's; a 360-day
        ! calendar (spec sec. 10) is wired there when the calendar lands.
        doy = ams%ts%rad%doy0 + real(ams%nstep + 1, wp)*ams%par%dt/86400.0_wp

        ! === Step 1-4: imports, boundary refresh, hand-over to the physics ===
        if (present(imports)) call import_copy(ams, imports)
        call aeros_boundary_refresh(ams%bounds, doy)
        call aeros_boundary_apply(ams)

        ! === Sub-cycle ========================================================
        do n = 1, nsub
            call aeros_timestep_step(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)
        end do

        ams%nstep = ams%nstep + nsub
        ams%time  = ams%time + real(nsub, wp)*ams%par%dt/real(year_seconds, wp)

        ! Leave the grid-space fields consistent with the spectral state.
        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)

        ! === Cycle-end exports and accumulator reset =========================
        if (present(exports)) call export_populate(ams, exports)

        ! Zero the library-side accumulators for the next cycle.
        !
        ! Wave 4a: ams%diag_acc IS allocated now (aeros_init), so this reset
        ! is no longer a no-op -- it wipes the buffers every coupling cycle.
        ! That is the per-cycle export contract, not rce_long's run-long
        ! second-half accumulation: a driver that needs run-long means must
        ! NOT accumulate into ams%diag_acc through aeros_lib_step yet. Wave 4c
        ! replaces this blanket reset with method-based accumulation
        ! (accumulate / finalize / reset as explicit calls on the container).
        call aeros_diag_acc_reset(ams%diag_acc)

        return

    end subroutine aeros_lib_step

    subroutine import_copy(ams, imports)
        ! Copy each allocated field of `imports` into its BC_EXTERNAL buffer
        ! on ams%bounds (bumping the buffer version). A field the host did
        ! not allocate is silently skipped -- the existing buffer (BC_FILE /
        ! BC_CONSTANT / previous external value) is used. Fields with no
        ! aeros-side reader today are copied anyway so the buffer is
        ! populated before Phase C wires them.

        implicit none

        type(aeros_class),     intent(inout) :: ams
        type(aeros_imports_t), intent(in)    :: imports

        if (allocated(imports%sst))         call push(ams%bounds%sst,         imports%sst)
        if (allocated(imports%sic))         call push(ams%bounds%sic,         imports%sic)
        if (allocated(imports%t_ice_sfc))   call push(ams%bounds%t_ice_sfc,   imports%t_ice_sfc)
        if (allocated(imports%t_skin_land)) call push(ams%bounds%t_skin_land, imports%t_skin_land)
        if (allocated(imports%q_skin_land)) call push(ams%bounds%q_skin_land, imports%q_skin_land)
        if (allocated(imports%albedo))      call push(ams%bounds%albedo,      imports%albedo)
        if (allocated(imports%z0m))         call push(ams%bounds%z0m,         imports%z0m)
        if (allocated(imports%so4_load))    call push(ams%bounds%so4,         imports%so4_load)
        if (allocated(imports%dust))        call push(ams%bounds%dust,        imports%dust)
        if (allocated(imports%lsm))         call push(ams%bounds%lsm,         imports%lsm)
        if (allocated(imports%phis))        call push(ams%bounds%phis,        imports%phis)

        ! Year-boundary scalars: `current` is what physics reads.
        if (imports%year_boundary) then
            if (ams%bounds%co2%kind     == BC_EXTERNAL) ams%bounds%co2%current     = imports%co2_ppm
            if (ams%bounds%ch4%kind     == BC_EXTERNAL) ams%bounds%ch4%current     = imports%ch4_ppb
            if (ams%bounds%n2o%kind     == BC_EXTERNAL) ams%bounds%n2o%current     = imports%n2o_ppb
            if (ams%bounds%time_bp%kind == BC_EXTERNAL) ams%bounds%time_bp%current = imports%time_bp
        end if

        return

    contains

        subroutine push(fld, values)
            type(aeros_bc_2d_t), intent(inout) :: fld
            real(wp),            intent(in)    :: values(:,:)
            if (fld%kind /= BC_EXTERNAL) return
            fld%buffer  = values
            fld%version = fld%version + 1
            return
        end subroutine push

    end subroutine import_copy

    subroutine export_populate(ams, exports)
        ! Populate the allocated fields of `exports` with cycle-mean values.
        !
        ! Wave-3.a: this is a shell. The cycle-mean accumulators live in
        ! ams%diag_acc, whose accumulate/finalize arithmetic ports from
        ! rce_long in a later wave. For now the routine leaves allocated
        ! export arrays at zero -- the caller sees a well-shaped bundle with
        ! zeros, which is enough to compile and link the API. Phase A
        ! calibration lands the field-by-field wire-up when a driver first
        ! consumes it.

        implicit none

        type(aeros_class),     intent(inout) :: ams
        type(aeros_exports_t), intent(out)   :: exports

        ! intent(out) leaves scalar fields at their type defaults and every
        ! allocatable component un-allocated. If the caller pre-allocated
        ! fields on the exports bundle they'd be dropped here -- so the
        ! contract for the shell is: caller does not pre-allocate exports;
        ! the wave that fills these fields will switch to intent(inout).
        !
        ! Keep the reference to `ams` alive so gfortran does not warn about
        ! unused-argument on this stub.
        if (ams%par%nlon < 0) return

        return
    end subroutine export_populate

    subroutine aeros_lib_finalize(ams, ierr)
        ! Delegate to aeros_end. Idempotent: safe to call on an ams that was
        ! never lib_init'd.

        implicit none

        type(aeros_class), intent(inout) :: ams
        integer,           intent(out)   :: ierr

        ierr = 0
        call aeros_end(ams)

        return
    end subroutine aeros_lib_finalize

    subroutine aeros_lib_restart_write(ams, filename, ierr)
        ! Write the checkpoint: aeros_timestep_write_restart unchanged (spec
        ! sec. 9: Phase A does not change the restart file format), plus one
        ! variable `bc_<field>` per BC_EXTERNAL 2D provider -- after a
        ! coupled restart the host has not pushed anything yet and the
        ! physics needs a valid starting field. BC_FILE / BC_CONSTANT
        ! buffers are re-derived at init and are not written, so a run
        ! without external providers writes exactly today's file.
        !
        ! diag_acc state is not persisted (the accumulator arithmetic is
        ! still rce_long's; Wave 4c).

        implicit none

        type(aeros_class), intent(inout) :: ams   ! inout (not in) so the timestep call is legal
        character(len=*),  intent(in)    :: filename
        integer,           intent(out)   :: ierr

        real(wp) :: time_s

        ierr = 0

        time_s = ams%time*real(year_seconds, wp)
        call aeros_timestep_write_restart(ams%ts, ams%now, time_s, filename)

        call write_bc(ams%bounds%sst,         "sst")
        call write_bc(ams%bounds%sic,         "sic")
        call write_bc(ams%bounds%t_ice_sfc,   "t_ice_sfc")
        call write_bc(ams%bounds%t_skin_land, "t_skin_land")
        call write_bc(ams%bounds%q_skin_land, "q_skin_land")
        call write_bc(ams%bounds%albedo,      "albedo")
        call write_bc(ams%bounds%z0m,         "z0m")
        call write_bc(ams%bounds%phis,        "phis")
        call write_bc(ams%bounds%lsm,         "lsm")
        call write_bc(ams%bounds%dust,        "dust")
        call write_bc(ams%bounds%so4,         "so4")
        call write_bc(ams%bounds%o3,          "o3")

        return

    contains

        subroutine write_bc(fld, name)
            type(aeros_bc_2d_t), intent(in) :: fld
            character(len=*),    intent(in) :: name
            if (fld%kind /= BC_EXTERNAL) return
            call nc_write(filename, "bc_"//name, fld%buffer, dim1="lon", dim2="lat", &
                          long_name="boundary provider buffer (BC_EXTERNAL): "//name)
            return
        end subroutine write_bc

    end subroutine aeros_lib_restart_write

    subroutine aeros_lib_restart_read(ams, filename, ierr)
        ! Read a checkpoint into an already-initialised `ams`. Wraps
        ! aeros_timestep_read_restart, then reconciles ams%time / ams%nstep
        ! from the file so subsequent aeros_lib_step calls advance from the
        ! restored point.

        implicit none

        type(aeros_class), intent(inout) :: ams
        character(len=*),  intent(in)    :: filename
        integer,           intent(out)   :: ierr

        real(wp) :: time_s

        ierr = 0

        call aeros_timestep_read_restart(ams%ts, ams%now, time_s, filename)

        ams%time  = time_s/real(year_seconds, wp)
        ams%nstep = ams%ts%nstep

        ! BC_EXTERNAL provider buffers (written by aeros_lib_restart_write).
        ! The physics working fields they feed were themselves restored by
        ! the timestep restart, so nothing is re-applied here; the next
        ! aeros_lib_step applies whatever the host pushes -- or, absent a
        ! push, these. A file without them (an older checkpoint, or a run
        ! that had no external provider) leaves the zero buffers of init.
        call read_bc(ams%bounds%sst,         "sst")
        call read_bc(ams%bounds%sic,         "sic")
        call read_bc(ams%bounds%t_ice_sfc,   "t_ice_sfc")
        call read_bc(ams%bounds%t_skin_land, "t_skin_land")
        call read_bc(ams%bounds%q_skin_land, "q_skin_land")
        call read_bc(ams%bounds%albedo,      "albedo")
        call read_bc(ams%bounds%z0m,         "z0m")
        call read_bc(ams%bounds%phis,        "phis")
        call read_bc(ams%bounds%lsm,         "lsm")
        call read_bc(ams%bounds%dust,        "dust")
        call read_bc(ams%bounds%so4,         "so4")
        call read_bc(ams%bounds%o3,          "o3")

        ! phis is not part of the timestep restart, and the stored buffer is
        ! the band-limited field aeros_boundary_apply left there: hand it to
        ! the dynamics as is and mark that version applied, so the next
        ! apply does not put it through a second (non-identical) round-trip.
        if (ams%bounds%phis%kind == BC_EXTERNAL .and. ams%bc_phis_applied /= ams%bounds%phis%version) then
            if (ams%phys%l_topography) call aeros_timestep_set_phis(ams%ts, ams%bounds%phis%buffer)
            ams%bc_phis_applied = ams%bounds%phis%version
        end if

        call aeros_timestep_diagnose(ams%ts, ams%pool, ams%vgrid, ams%grd, ams%now)

        return

    contains

        subroutine read_bc(fld, name)
            type(aeros_bc_2d_t), intent(inout) :: fld
            character(len=*),    intent(in)    :: name
            if (fld%kind /= BC_EXTERNAL) return
            if (.not. nc_exists_var(filename, "bc_"//name)) return
            call nc_read(filename, "bc_"//name, fld%buffer)
            fld%version = fld%version + 1
            return
        end subroutine read_bc

    end subroutine aeros_lib_restart_read

    subroutine aeros_lib_diagnose(ams, diag_out, ierr)
        ! Return a diagnostic snapshot without stepping.
        !
        ! Phase A ships an empty shell (docs/coupling_api.md sec. 6.3);
        ! Phase C fills aeros_diag_t with the rce_long accumulator contents.
        ! For now the routine just clears the shell and reports success.

        implicit none

        type(aeros_class), intent(in)  :: ams
        type(aeros_diag_t), intent(out) :: diag_out
        integer,            intent(out) :: ierr

        ierr = 0
        diag_out%reserved = 0

        ! Keep `ams` referenced so the shell does not warn on unused arg.
        if (ams%nstep < 0) ierr = 1

        return
    end subroutine aeros_lib_diagnose

end module aeros
