program test_adaptive_diff
    ! Acceptance test for the two opt-in adaptive-hyperdiffusion mechanisms
    ! (aeros_timestep header, "Adaptive hyperdiffusion"). Both default OFF and
    ! both must reduce, bit for bit, to the fixed del^ndiff scheme -- that is
    ! the property the whole feature is built around, so it is what is asserted
    ! first and most tightly.
    !
    ! Three things:
    !
    ! (a) BOTH OFF => THE FIXED del^ndiff SCHEME, BIT FOR BIT. With the taper
    !     off, the per-level ratio dratio_lev(:,k) must equal the 1-D dratio(:)
    !     to the last bit at every level, and the vorticity multiplier must be
    !     exactly 1 -- so a full integration with both flags off is identical to
    !     one that never knew about them.
    !
    ! (b) THE TAPER IS A NO-OP WHEN diff_ndiff_top == ndiff. Turning the taper
    !     on with the top order equal to the interior order must leave
    !     dratio_lev == dratio bit for bit; dropping the top order to del^4 must
    !     then change the top levels (a lower order = a larger interior-degree
    !     ratio = broader-scale damping) and NOT the tropospheric ones.
    !
    ! (c) THE VORTICITY SCALING RAISES THE EFFECTIVE COEFFICIENT ON A SPIKING
    !     LEVEL, AND ONLY THERE. Two runs from one initial state differ only in
    !     whether diff_adapt is on. After a single step the level whose RMS
    !     vorticity exceeds diff_zeta_ref is damped strictly harder; the level
    !     below the threshold is bit-for-bit unchanged (its multiplier is 1).
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, p0, &
                                aeros_param_class, aeros_grid_class, &
                                aeros_spec_class, aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_timestep

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12    ! the L12 top the taper is designed for

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg

    integer :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_adaptive_diff:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    call test_off_is_plain(nfail)
    call test_taper_reduces(nfail)
    call test_vorticity_scaling(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    if (nfail == 0) then
        write(*,*) ""
        write(*,*) " test_adaptive_diff:: ALL PASSED"
    else
        write(*,"(a,i0,a)") "  test_adaptive_diff:: ", nfail, " FAILED"
        error stop 1
    end if

contains

    subroutine build_par(par)
        implicit none
        type(aeros_param_class), intent(out) :: par
        par%trunc = trunc; par%nlon = -1; par%nlat = -1; par%nlev = nlev
        par%nthreads = -1
        par%dt = 1800.0_wp
        par%semi_implicit = .TRUE.
        par%held_suarez   = .FALSE.
        par%eps_filter = 0.0_wp
        par%raw_alpha  = 0.53_wp
        par%ndiff = 6
        par%tau_diff = 6.0_wp
        par%mass_fixer = .FALSE.
        return
    end subroutine build_par

    subroutine seed_rest(now, s, lm_pert, vor_lev)
        ! A resting isothermal atmosphere (280 K, p0) plus one vorticity mode at
        ! lm_pert with a per-level amplitude vor_lev(:).
        implicit none
        type(aeros_state_class), intent(inout) :: now
        type(aeros_sht_class),   intent(in)    :: s
        integer,  intent(in) :: lm_pert
        real(wp), intent(in) :: vor_lev(nlev)
        real(dp) :: y00
        integer  :: k
        y00 = sqrt(4.0_dp*acos(-1.0_dp))
        call aeros_spec_zero(now%spec)
        do k = 1, nlev
            now%spec%temp(aeros_sht_lm(s,0,0),k) = cmplx(280.0_dp*y00, 0.0_dp, wp_sh)
            now%spec%vor(lm_pert,k) = cmplx(real(vor_lev(k),dp), 0.0_dp, wp_sh)
        end do
        now%spec%lnps(aeros_sht_lm(s,0,0)) = cmplx(log(real(p0,dp))*y00, 0.0_dp, wp_sh)
        return
    end subroutine seed_rest

    ! === (a) both off == the fixed del^ndiff scheme ==========================
    subroutine test_off_is_plain(nfail)
        implicit none
        integer, intent(inout) :: nfail
        type(aeros_param_class)    :: par
        type(aeros_timestep_class) :: ts
        integer :: k, l, lmax
        logical :: exact

        write(*,*) ""
        write(*,*) " -- (a) both flags off reduce to the 1-D dratio"

        call build_par(par)
        call aeros_timestep_init(ts, par, pool, grd, vg)
        lmax = ubound(ts%dratio, 1)

        call check(.not. ts%diff_adapt .and. .not. ts%diff_taper, &
                    "adaptive mechanisms default OFF", nfail)

        exact = .TRUE.
        do k = 1, nlev
            do l = 0, lmax
                if (ts%dratio_lev(l,k) /= ts%dratio(l)) exact = .FALSE.
            end do
        end do
        call check(exact, "dratio_lev(:,k) == dratio(:) bit for bit at every level", nfail)

        exact = .TRUE.
        do k = 1, nlev
            if (ts%ndiff_lev(k) /= ts%ndiff) exact = .FALSE.
        end do
        call check(exact, "every level's order is ndiff when the taper is off", nfail)

        call aeros_timestep_end(ts)
        return
    end subroutine test_off_is_plain

    ! === (b) taper is a no-op at diff_ndiff_top == ndiff =====================
    subroutine test_taper_reduces(nfail)
        implicit none
        integer, intent(inout) :: nfail
        type(aeros_param_class)    :: par
        type(aeros_timestep_class) :: ts
        integer :: k, l, lmax
        logical :: exact, top_changed, tropo_same

        write(*,*) ""
        write(*,*) " -- (b) taper reduces to the 1-D dratio when ndiff_top == ndiff"

        call build_par(par)
        call aeros_timestep_init(ts, par, pool, grd, vg)
        lmax = ubound(ts%dratio, 1)

        ! Taper ON but top order == interior order: must be a bit-for-bit no-op.
        call aeros_timestep_set_diff_taper(ts, vg, .TRUE., ndiff_top=ts%ndiff, &
                                            taper_sigma=0.15_wp)
        exact = .TRUE.
        do k = 1, nlev
            do l = 0, lmax
                if (ts%dratio_lev(l,k) /= ts%dratio(l)) exact = .FALSE.
            end do
        end do
        call check(exact, "diff_ndiff_top == ndiff leaves dratio_lev == dratio", nfail)

        ! Drop the top order to del^4: the top level(s) must change (a lower
        ! order raises the interior-degree ratio), a tropospheric level must not.
        call aeros_timestep_set_diff_taper(ts, vg, .TRUE., ndiff_top=4, &
                                            taper_sigma=0.15_wp)
        top_changed = (ts%ndiff_lev(1) < ts%ndiff) .and. &
                        (ts%dratio_lev(1,1) > ts%dratio(1))
        call check(top_changed, "del^4 top: order lowered and interior-l ratio raised", nfail)

        tropo_same = .TRUE.
        do k = 1, nlev
            if (vg%sigma_full(k) >= 0.15_wp) then
                if (ts%ndiff_lev(k) /= ts%ndiff) tropo_same = .FALSE.
                do l = 0, lmax
                    if (ts%dratio_lev(l,k) /= ts%dratio(l)) tropo_same = .FALSE.
                end do
            end if
        end do
        call check(tropo_same, "levels below the taper sigma keep the ndiff ratio", nfail)

        call aeros_timestep_end(ts)
        return
    end subroutine test_taper_reduces

    ! === (c) vorticity scaling raises the coefficient where zeta spikes ======
    subroutine test_vorticity_scaling(nfail)
        implicit none
        integer, intent(inout) :: nfail
        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: now_off, now_on
        type(aeros_timestep_class) :: ts_off, ts_on
        type(aeros_sht_class), pointer :: s
        real(wp) :: phis(grd%nlon,grd%nlat), vor_lev(nlev)
        real(dp) :: a_off_hi, a_on_hi, a_off_lo, a_on_lo
        integer  :: lm_pert, k_hi, k_lo

        write(*,*) ""
        write(*,*) " -- (c) vorticity scaling damps a spiking level harder, only there"

        s => pool%sht(1)
        k_hi = 1                 ! model top: give it a large vorticity
        k_lo = nlev              ! near surface: keep it small
        lm_pert = aeros_sht_lm(s, s%lmax, 3)

        ! RMS|zeta| ~ sqrt(2)|vor|/sqrt(4pi): 3e-3 -> ~1.2e-3 (>> ref), 1e-6 -> tiny.
        vor_lev = 0.0_wp
        vor_lev(k_hi) = 3.0e-3_wp
        vor_lev(k_lo) = 1.0e-6_wp
        phis = 0.0_wp

        ! Run A: adaptive OFF.
        call build_par(par)
        call aeros_state_alloc(now_off, grd, s%nlm, nlev)
        call aeros_timestep_init(ts_off, par, pool, grd, vg)
        call seed_rest(now_off, s, lm_pert, vor_lev)
        call aeros_timestep_set_phis(ts_off, phis)
        call aeros_timestep_step(ts_off, pool, vg, grd, now_off)
        a_off_hi = abs(now_off%spec%vor(lm_pert,k_hi))
        a_off_lo = abs(now_off%spec%vor(lm_pert,k_lo))

        ! Run B: adaptive ON, same initial state and dynamics.
        call aeros_state_alloc(now_on, grd, s%nlm, nlev)
        call aeros_timestep_init(ts_on, par, pool, grd, vg)
        ts_on%diff_adapt      = .TRUE.
        ts_on%diff_zeta_ref   = 1.0e-4_wp
        ts_on%diff_adapt_gain = 1.0_wp
        ts_on%diff_adapt_max  = 10.0_wp
        call seed_rest(now_on, s, lm_pert, vor_lev)
        call aeros_timestep_set_phis(ts_on, phis)
        call aeros_timestep_step(ts_on, pool, vg, grd, now_on)
        a_on_hi = abs(now_on%spec%vor(lm_pert,k_hi))
        a_on_lo = abs(now_on%spec%vor(lm_pert,k_lo))

        write(*,"(a40,2es14.6)") "   high-vor level |vor| off / on ", a_off_hi, a_on_hi
        write(*,"(a40,2es14.6)") "   low-vor  level |vor| off / on ", a_off_lo, a_on_lo

        call check(a_on_hi < 0.999_dp*a_off_hi, &
                    "spiking level is damped strictly harder with diff_adapt on", nfail)
        call check(a_on_lo == a_off_lo, &
                    "sub-threshold level is bit-for-bit unchanged (multiplier = 1)", nfail)

        call aeros_timestep_end(ts_off)
        call aeros_timestep_end(ts_on)
        call aeros_state_end(now_off)
        call aeros_state_end(now_on)
        return
    end subroutine test_vorticity_scaling

    subroutine check(ok, label, nfail)
        implicit none
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        integer, intent(inout) :: nfail
        if (ok) then
            write(*,*) "  ok   : ", trim(label)
        else
            write(*,*) "  FAIL : ", trim(label)
            nfail = nfail + 1
        end if
        return
    end subroutine check

end program test_adaptive_diff
