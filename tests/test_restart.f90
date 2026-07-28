program test_restart
    ! Acceptance test for checkpoint/restart (aeros_timestep_write_restart /
    ! aeros_timestep_read_restart): a restart must resume BIT-IDENTICALLY.
    !
    ! The correctness bar is a split-run equivalence. With everything on that
    ! the restart file carries -- moisture, a slab ocean, ecCKD radiation -- we:
    !
    !   1. run K steps continuously                      (the reference), and
    !   2. SEPARATELY run K/2 steps, write a restart, read it into a FRESHLY
    !      initialised timestep, and run the remaining K/2.
    !
    ! The two final states must be bit-for-bit equal: the spectral prognostics
    ! vor/div/temp/lnps at BOTH leapfrog time levels (X^n = now%spec and
    ! X^(n-1) = ts%old), the gridpoint humidity qv_g, the ocean SST, and the
    ! scalars (nstep, mass_target, lnr_cum). Any field the restart drops or
    ! rounds shows up here as a non-zero difference.
    !
    ! The cadence is chosen so the FIRST restarted step APPLIES the cached
    ! radiative heating rather than recomputing it (rad_interval = 2 dt, restart
    ! at an odd nstep): that is the step that would differ if ts%rad%heat were
    ! not saved, so this configuration actively guards the radiation-cache save.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, p0, R_d, grav, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_state_class, aeros_spec_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_ocean,    only : aeros_ocean_init, OCEAN_SLAB
    use aeros_timestep

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12
    integer, parameter :: nstepK = 10          ! total reference steps
    integer, parameter :: nhalf  = nstepK/2    ! restart written here (odd -> cache applied)
    real(wp), parameter :: dt = 1800.0_wp

    character(len=*), parameter :: rstfile = "test_restart_chk.nc"

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)     :: grd
    type(aeros_vgrid_class)    :: vg

    ! Reference (continuous) run, and the restart run.
    type(aeros_timestep_class) :: ts_ref, ts_a, ts_b
    type(aeros_state_class)    :: now_ref, now_a, now_b

    ! Saved reference final state.
    type(aeros_spec_class) :: ref_now, ref_old
    real(wp), allocatable  :: ref_qv(:,:,:), ref_sst(:,:)
    integer                :: ref_nstep
    real(dp)               :: ref_mass_target, ref_lnr_cum

    real(wp) :: time_a, time_b
    integer  :: nfail, n

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_restart:: T", trunc, " L", nlev, &
                                       "  grid ", grd%nlon, "x", grd%nlat
    write(*,"(a,i0,a,i0)") "   K=", nstepK, "  restart after ", nhalf

    ! --- 1. Reference: K continuous steps ---------------------------------
    call setup(ts_ref, now_ref)
    do n = 1, nstepK
        call aeros_timestep_step(ts_ref, pool, vg, grd, now_ref)
    end do
    ! Snapshot the reference final state.
    call aeros_spec_alloc(ref_now, pool%sht(1)%nlm, nlev)
    call aeros_spec_alloc(ref_old, pool%sht(1)%nlm, nlev)
    call aeros_spec_copy(ref_now, now_ref%spec)
    call aeros_spec_copy(ref_old, ts_ref%old)
    allocate(ref_qv(grd%nlon,grd%nlat,nlev), ref_sst(grd%nlon,grd%nlat))
    ref_qv  = now_ref%qv_g
    ref_sst = ts_ref%ocn%sst
    ref_nstep       = ts_ref%nstep
    ref_mass_target = ts_ref%mass_target
    ref_lnr_cum     = ts_ref%lnr_cum

    ! --- 2a. Restart producer: K/2 steps, then checkpoint -----------------
    call setup(ts_a, now_a)
    do n = 1, nhalf
        call aeros_timestep_step(ts_a, pool, vg, grd, now_a)
    end do
    call aeros_timestep_write_restart(ts_a, now_a, real(nhalf,wp)*dt, rstfile)
    write(*,"(a,a)") "   wrote restart -> ", rstfile

    ! --- 2b. Fresh timestep, read the checkpoint, finish the run ----------
    ! ts_b/now_b are built by the SAME cold-start path, then completely
    ! overwritten from the file: this is the "read into a freshly initialised
    ! timestep" contract.
    call setup(ts_b, now_b)
    call aeros_timestep_read_restart(ts_b, now_b, time_b, rstfile)

    write(*,"(a,i0,a,f9.1)") "   restored nstep=", ts_b%nstep, "  time=", time_b
    call check(ts_b%nstep == nhalf, "restored nstep matches the checkpoint", nfail)
    call check(time_b == real(nhalf,wp)*dt, "restored model time matches", nfail)

    do n = nhalf+1, nstepK
        call aeros_timestep_step(ts_b, pool, vg, grd, now_b)
    end do

    ! --- 3. The two final states must be bit-identical --------------------
    write(*,*) ""
    write(*,*) " -- restarted final state vs continuous reference"

    call check_spec("now (X^n)  ", now_b%spec, ref_now, nfail)
    call check_spec("old (X^n-1)", ts_b%old,   ref_old, nfail)

    call check(maxval(abs(now_b%qv_g - ref_qv)) == 0.0_wp, &
               "qv_g bit-identical", nfail)
    call check(maxval(abs(ts_b%ocn%sst - ref_sst)) == 0.0_wp, &
               "ocn%sst bit-identical", nfail)

    call check(ts_b%nstep == ref_nstep, "final nstep matches reference", nfail)
    call check(ts_b%mass_target == ref_mass_target, "mass_target matches", nfail)
    call check(ts_b%lnr_cum == ref_lnr_cum, "lnr_cum matches", nfail)

    ! --- cleanup ----------------------------------------------------------
    call aeros_timestep_end(ts_ref); call aeros_state_end(now_ref)
    call aeros_timestep_end(ts_a);   call aeros_state_end(now_a)
    call aeros_timestep_end(ts_b);   call aeros_state_end(now_b)
    call aeros_spec_end(ref_now);    call aeros_spec_end(ref_old)
    deallocate(ref_qv, ref_sst)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_restart:: PASS"
    else
        write(*,*) "test_restart:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine setup(ts, now)
        ! Cold-start the timestep and state identically for every run: the full
        ! moist + slab-ocean + ecCKD-radiation stack, and the warm, humid,
        ! conditionally-unstable initial condition of the RCE driver. The
        ! radiation cadence is 2 dt so the restart (at an odd nstep) lands on a
        ! cache-apply step.
        type(aeros_timestep_class), intent(inout) :: ts
        type(aeros_state_class),    intent(inout) :: now

        type(aeros_param_class) :: par
        real(wp), allocatable :: phis2(:,:)
        real(wp) :: tval, qs, dqsdt
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev)
        integer  :: i, j, k

        par%trunc = trunc; par%nlon = -1; par%nlat = -1; par%nlev = nlev
        par%nthreads = -1
        par%dt = dt
        par%semi_implicit = .TRUE.
        par%held_suarez   = .FALSE.
        par%eps_filter = 0.06_wp
        par%raw_alpha  = 0.53_wp
        par%ndiff = 6; par%tau_diff = 6.0_wp
        par%mass_fixer = .FALSE.

        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)

        ! Physics: surface fluxes, convection, condensation, radiation on.
        ts%cnd%enabled = .TRUE.
        ts%cnv%enabled = .TRUE.
        ts%cnv%dry_adjust = .TRUE.
        ts%surf%enabled = .TRUE.
        ts%surf%c_d = 1.5e-3_wp
        ts%rad%enabled = .TRUE.
        ts%rad%interval = 2.0_wp*dt        ! nrad = 2 -> odd-nstep step applies cache
        ts%sponge_on = .TRUE.
        call aeros_timestep_set_sponge(ts, vg, ts%sponge_kr, ts%sponge_kt, ts%sponge_sigma)

        ! Slab ocean so the SST is a genuinely prognostic restart field.
        ts%ocn%mode  = OCEAN_SLAB
        ts%ocn%depth = 10.0_wp
        call aeros_ocean_init(ts%ocn, grd)

        allocate(phis2(grd%nlon, grd%nlat)); phis2 = 0.0_wp
        call aeros_timestep_set_phis(ts, phis2)
        deallocate(phis2)

        ! Warm, humid, conditionally-unstable start (as rce_long / test_rce).
        call aeros_spec_zero(now%spec)
        do k = 1, nlev
            tval = 300.0_wp - 90.0_wp*(1.0_wp - ((real(k,wp) - 0.5_wp)/real(nlev,wp))**0.6_wp)
            now%spec%temp(aeros_sht_lm(pool%sht(1),0,0),k) = &
                    cmplx(real(tval,dp)*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
        end do
        now%spec%lnps(aeros_sht_lm(pool%sht(1),0,0)) = &
                cmplx(log(real(p0,dp))*sqrt(16.0_dp*atan(1.0_dp)), 0.0_dp, wp_sh)
        now%spec%temp(aeros_sht_lm(pool%sht(1),2,0),:) = cmplx(-5.0_dp, 0.0_dp, wp_sh)

        call aeros_timestep_diagnose(ts, pool, vg, grd, now)
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                call aeros_vgrid_pressure(vg, now%ps(i,j), phalf(0:nlev), &
                                          pfull(1:nlev), dpc(1:nlev))
                do k = 1, nlev
                    call aeros_qsat(now%temp_g(i,j,k), pfull(k), qs, dqsdt)
                    now%qv_g(i,j,k) = 0.9_wp*qs
                end do
            end do
        end do

        return
    end subroutine setup

    subroutine check_spec(label, a, b, nfail)
        ! Bit-identity of every spectral prognostic between two level sets.
        character(len=*),       intent(in)    :: label
        type(aeros_spec_class), intent(in)    :: a, b
        integer,                intent(inout) :: nfail
        logical :: ok
        ok = maxval(abs(a%vor  - b%vor))  == 0.0_wp .and. &
             maxval(abs(a%div  - b%div))  == 0.0_wp .and. &
             maxval(abs(a%temp - b%temp)) == 0.0_wp .and. &
             maxval(abs(a%lnps - b%lnps)) == 0.0_wp
        call check(ok, "spectral "//label//" bit-identical (vor/div/temp/lnps)", nfail)
        return
    end subroutine check_spec

    subroutine check(ok, label, nfail)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (ok) then
            write(*,"(a,a)") "   ok   : ", label
        else
            write(*,"(a,a)") "   FAIL : ", label
            nfail = nfail + 1
        end if
        return
    end subroutine check

end program test_restart
