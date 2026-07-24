program test_correction
    ! Acceptance test for the additive tendency-correction framework.
    !
    ! At M2 there is no real correction to test -- the only term is a no-op --
    ! so what has to be established is that the FRAMEWORK is correct and inert,
    ! in that order. Both halves matter and they fail differently: a framework
    ! that is wired in but wrong will produce a plausible, wrong DeltaF at M2b;
    ! a framework that is not actually inert will have been silently changing
    ! every result taken between now and then.
    !
    ! 1. THE `DeltaF = 0` TWIN IS BIT-EXACT. design.md section 3.8 mitigation 3
    !    requires every production run to have a twin with the correction off,
    !    and requires the comparison to mean something. So `correction=.FALSE.`
    !    with a FILLED, enabled term must reproduce the uncorrected integration
    !    to the last bit -- not to a tolerance. Tested over an integration, not
    !    a single call, because that is how it will be used.
    !
    ! 2. A REGISTERED BUT UNFILLED TERM IS INERT. The M2 configuration is a
    !    no-op term with the master switch ON. It exercises every line of the
    !    apply path -- the loop runs, the mask is read, the addition happens --
    !    and must still change nothing at all. This is the check that separates
    !    "the plumbing is in" from "a correction is being applied".
    !
    ! 3. A FILLED TERM ADDS EXACTLY WHAT IT HOLDS. Sum, not scaled sum: aeros
    !    stores corrections in the tendency's own units, and the whole point of
    !    writing that convention down is that mwm/B_multires lost time to
    !    SpeedyWeather's (B2B_OUTCOME.md). Checked as an exact equality.
    !
    ! 4. THE SPECTRAL-SCALE SELECTOR CUTS WHERE IT SAYS. `lcut` is section
    !    3.7's "apply only at large scales" and B2B_OUTCOME.md's first
    !    mitigation for the instability the whole-field version hit. Every
    !    coefficient with l <= lcut must be corrected and every one with
    !    l > lcut must be untouched -- both directions, since a selector that
    !    silently passed everything would look like a working correction.
    !
    ! 5. PER-FIELD SELECTION IS RESPECTED. A term that acts on vorticity alone
    !    must leave divergence, temperature and ln(p_s) exactly as the dynamics
    !    left them.
    !
    ! Exits non-zero on failure.

    use aeros_defs,       only : dp, wp, wp_sh, p0, &
                                    aeros_param_class, aeros_grid_class, &
                                    aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_state
    use aeros_vertical
    use aeros_tendency
    use aeros_correction
    use aeros_timestep

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 8

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg

    integer :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_correction:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    call test_apply_exact(nfail)
    call test_scale_selector(nfail)
    call test_field_selector(nfail)
    call test_noop_is_inert(nfail)
    call test_twin_is_bit_exact(nfail)

    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_correction:: PASS"
    else
        write(*,*) "test_correction:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine fill_tend(tnd, seed)
        ! A deterministic, non-trivial tendency to correct.

        implicit none

        type(aeros_tend_class), intent(inout) :: tnd
        real(dp), intent(in) :: seed

        integer :: k, lm

        do k = 1, nlev
            do lm = 1, pool%sht(1)%nlm
                tnd%vor(lm,k)  = cmplx(seed*real(lm+k,dp), seed*real(lm-k,dp), wp_sh)
                tnd%div(lm,k)  = cmplx(seed*real(2*lm,dp), seed*real(k,dp),    wp_sh)
                tnd%temp(lm,k) = cmplx(seed*real(lm*k,dp), 0.0_dp,             wp_sh)
            end do
        end do
        do lm = 1, pool%sht(1)%nlm
            tnd%lnps(lm) = cmplx(seed*real(lm,dp), seed, wp_sh)
        end do

        return

    end subroutine fill_tend

    subroutine fill_term_vor(cor, it, val)
        ! Put a known constant into a term's vorticity correction.

        implicit none

        type(aeros_correction_class), intent(inout) :: cor
        integer,  intent(in) :: it
        real(dp), intent(in) :: val

        cor%term(it)%vor = cmplx(val, 0.5_dp*val, wp_sh)

        return

    end subroutine fill_term_vor

    ! === 3. A filled term adds exactly what it holds =========================

    subroutine test_apply_exact(nfail)

        implicit none

        integer, intent(inout) :: nfail

        type(aeros_correction_class) :: cor
        type(aeros_tend_class) :: tnd, ref
        real(dp), parameter :: val = 3.25e-9_dp
        integer  :: it, k, lm, nlm, nwrong

        write(*,*) ""
        write(*,*) " -- a filled term adds exactly what it holds"

        nlm = pool%sht(1)%nlm

        call aeros_tend_alloc(tnd, nlm, nlev)
        call aeros_tend_alloc(ref, nlm, nlev)
        call fill_tend(tnd, 1.0e-7_dp)
        call fill_tend(ref, 1.0e-7_dp)

        call aeros_correction_init(cor, pool%sht(1), nlev)
        cor%enabled = .TRUE.
        it = aeros_correction_add_term(cor, pool%sht(1), "constant", &
                                        on_vor=.TRUE., on_div=.FALSE., &
                                        on_temp=.FALSE., on_lnps=.FALSE., lcut=-1)
        call fill_term_vor(cor, it, val)

        call aeros_correction_apply(cor, tnd)

        ! Compared as `result == ref + correction`, NOT as
        ! `result - ref - correction == 0`. The second is not the same test: it
        ! is a statement about floating-point reversibility, and it fails at
        ! ~1e-21 here purely because the tendency is ~1e-4 and the correction
        ! ~1e-9. What has to be exact is that the code performs this addition
        ! and no other operation -- no scaling, no accumulation order of its
        ! own -- and that is what this form asserts.
        nwrong = 0
        do k = 1, nlev
            do lm = 1, nlm
                if (tnd%vor(lm,k) /= ref%vor(lm,k) + cmplx(val, 0.5_dp*val, wp_sh)) &
                    nwrong = nwrong + 1
            end do
        end do

        write(*,"(a,i0,a,i0)") "    coefficients not exactly ref + corr: ", nwrong, &
                                " of ", nlev*nlm
        call check(nwrong == 0, "the correction is added exactly, with no scaling", nfail)

        call aeros_correction_end(cor)
        call aeros_tend_end(tnd)
        call aeros_tend_end(ref)

        return

    end subroutine test_apply_exact

    ! === 4. The spectral-scale selector ======================================

    subroutine test_scale_selector(nfail)

        implicit none

        integer, intent(inout) :: nfail

        type(aeros_correction_class) :: cor
        type(aeros_tend_class) :: tnd, ref
        type(aeros_sht_class), pointer :: s
        real(dp), parameter :: val = 1.5e-9_dp
        integer, parameter  :: lcut = 5
        integer  :: it, k, lm, nlm, n_in, n_out, bad_in, bad_out

        write(*,*) ""
        write(*,"(a,i0)") "  -- the scale selector, lcut = ", lcut

        s => pool%sht(1)
        nlm = s%nlm

        call aeros_tend_alloc(tnd, nlm, nlev)
        call aeros_tend_alloc(ref, nlm, nlev)
        call fill_tend(tnd, 1.0e-7_dp)
        call fill_tend(ref, 1.0e-7_dp)

        call aeros_correction_init(cor, s, nlev)
        cor%enabled = .TRUE.
        it = aeros_correction_add_term(cor, s, "large-scale", &
                                        on_vor=.TRUE., on_div=.FALSE., &
                                        on_temp=.FALSE., on_lnps=.FALSE., lcut=lcut)
        call fill_term_vor(cor, it, val)

        call aeros_correction_apply(cor, tnd)

        ! Both directions. A selector that passed everything would satisfy the
        ! first count and fail the second; one that passed nothing, the reverse.
        n_in = 0; n_out = 0; bad_in = 0; bad_out = 0
        do k = 1, nlev
            do lm = 1, nlm
                if (s%l_of_lm(lm) <= lcut) then
                    n_in = n_in + 1
                    if (tnd%vor(lm,k) /= ref%vor(lm,k) + cmplx(val, 0.5_dp*val, wp_sh)) &
                        bad_in = bad_in + 1
                else
                    n_out = n_out + 1
                    if (tnd%vor(lm,k) /= ref%vor(lm,k)) bad_out = bad_out + 1
                end if
            end do
        end do

        write(*,"(a,i0,a,i0,a)") "    coefficients at l <= lcut     ", n_in,  " (", bad_in,  " wrong)"
        write(*,"(a,i0,a,i0,a)") "    coefficients at l >  lcut     ", n_out, " (", bad_out, " wrong)"

        call check(n_in > 0 .and. n_out > 0, &
                    "the truncation actually straddles lcut (the test is not vacuous)", nfail)
        call check(bad_in  == 0, "every coefficient at l <= lcut is corrected", nfail)
        call check(bad_out == 0, "every coefficient at l >  lcut is untouched", nfail)

        call aeros_correction_end(cor)
        call aeros_tend_end(tnd)
        call aeros_tend_end(ref)

        return

    end subroutine test_scale_selector

    ! === 5. Per-field selection ==============================================

    subroutine test_field_selector(nfail)

        implicit none

        integer, intent(inout) :: nfail

        type(aeros_correction_class) :: cor
        type(aeros_tend_class) :: tnd, ref
        real(dp) :: d_div, d_temp, d_lnps, d_vor
        integer  :: it, k, lm, nlm

        write(*,*) ""
        write(*,*) " -- per-field selection: vorticity only"

        nlm = pool%sht(1)%nlm

        call aeros_tend_alloc(tnd, nlm, nlev)
        call aeros_tend_alloc(ref, nlm, nlev)
        call fill_tend(tnd, 1.0e-7_dp)
        call fill_tend(ref, 1.0e-7_dp)

        call aeros_correction_init(cor, pool%sht(1), nlev)
        cor%enabled = .TRUE.
        it = aeros_correction_add_term(cor, pool%sht(1), "vor-only", &
                                        on_vor=.TRUE., on_div=.FALSE., &
                                        on_temp=.FALSE., on_lnps=.FALSE., lcut=-1)
        call fill_term_vor(cor, it, 2.0e-9_dp)

        call aeros_correction_apply(cor, tnd)

        d_vor = 0.0_dp; d_div = 0.0_dp; d_temp = 0.0_dp; d_lnps = 0.0_dp
        do k = 1, nlev
            do lm = 1, nlm
                d_vor  = max(d_vor,  abs(tnd%vor(lm,k)  - ref%vor(lm,k)))
                d_div  = max(d_div,  abs(tnd%div(lm,k)  - ref%div(lm,k)))
                d_temp = max(d_temp, abs(tnd%temp(lm,k) - ref%temp(lm,k)))
            end do
        end do
        do lm = 1, nlm
            d_lnps = max(d_lnps, abs(tnd%lnps(lm) - ref%lnps(lm)))
        end do

        write(*,"(a40,es12.3)") "   max |d vor|  (should be > 0)    ", d_vor
        write(*,"(a40,es12.3)") "   max |d div|                     ", d_div
        write(*,"(a40,es12.3)") "   max |d temp|                    ", d_temp
        write(*,"(a40,es12.3)") "   max |d lnps|                    ", d_lnps

        call check(d_vor > 0.0_dp, "the selected field is corrected", nfail)
        call check(d_div == 0.0_dp .and. d_temp == 0.0_dp .and. d_lnps == 0.0_dp, &
                    "the unselected fields are untouched", nfail)

        call aeros_correction_end(cor)
        call aeros_tend_end(tnd)
        call aeros_tend_end(ref)

        return

    end subroutine test_field_selector

    ! === 2. The no-op term is inert ==========================================

    subroutine test_noop_is_inert(nfail)
        ! The M2 shipped configuration: master switch ON, one registered term
        ! that is never filled. Integrated, not just applied once.

        implicit none

        integer, intent(inout) :: nfail

        real(dp) :: worst

        write(*,*) ""
        write(*,*) " -- a registered but unfilled term, master switch ON, 20 steps"

        call integrate_diff(corr_on=.TRUE., fill=.FALSE., worst=worst)

        write(*,"(a40,es12.3)") "   max |state - uncorrected state| ", worst
        call check(worst == 0.0_dp, &
                    "a no-op term changes nothing, bit for bit", nfail)

        return

    end subroutine test_noop_is_inert

    ! === 1. The twin is bit-exact ============================================

    subroutine test_twin_is_bit_exact(nfail)
        ! The master switch off, with a term that is enabled and FILLED. This
        ! is the configuration section 3.8 mitigation 3 asks for, and the
        ! filling is what makes it a real test: if the switch were only
        ! consulted per-term, or consulted after the addition, this would fail.

        implicit none

        integer, intent(inout) :: nfail

        real(dp) :: worst_off, worst_on

        write(*,*) ""
        write(*,*) " -- a FILLED term with the master switch off, 20 steps"

        call integrate_diff(corr_on=.FALSE., fill=.TRUE., worst=worst_off)
        call integrate_diff(corr_on=.TRUE.,  fill=.TRUE., worst=worst_on)

        write(*,"(a40,es12.3)") "   switch off: |state - reference| ", worst_off
        write(*,"(a40,es12.3)") "   switch on:  |state - reference| ", worst_on

        call check(worst_off == 0.0_dp, &
                    "correction = .FALSE. reproduces the uncorrected run bit for bit", nfail)
        call check(worst_on > 0.0_dp, &
                    "and with the switch on the same term does change the run", nfail)

        return

    end subroutine test_twin_is_bit_exact

    subroutine integrate_diff(corr_on, fill, worst)
        ! Integrate 20 steps twice -- once with no correction layer at all,
        ! once with one configured as asked -- and return the largest spectral
        ! difference between the two final states.

        implicit none

        logical,  intent(in)  :: corr_on
        logical,  intent(in)  :: fill
        real(dp), intent(out) :: worst

        type(aeros_param_class)    :: par
        type(aeros_state_class)    :: a, b
        type(aeros_timestep_class) :: ts
        real(wp) :: phis(grd%nlon,grd%nlat)
        integer, parameter :: nstep = 20
        integer :: n, k, lm, it, nlm

        nlm  = pool%sht(1)%nlm
        phis = 0.0_wp

        call default_par(par)

        ! -- reference: no correction layer configured at all
        call aeros_state_alloc(a, grd, nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)
        call start_state(a)
        call aeros_timestep_set_phis(ts, phis)
        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, a)
        end do
        call aeros_timestep_end(ts)

        ! -- the configured one
        call aeros_state_alloc(b, grd, nlm, nlev)
        call aeros_timestep_init(ts, par, pool, grd, vg)
        ts%cor%enabled = corr_on
        it = aeros_correction_add_term(ts%cor, pool%sht(1), "test", &
                                        on_vor=.TRUE., on_div=.FALSE., &
                                        on_temp=.TRUE., on_lnps=.FALSE., lcut=-1)
        if (fill) then
            ts%cor%term(it)%vor  = cmplx(1.0e-12_dp, 0.0_dp, wp_sh)
            ts%cor%term(it)%temp = cmplx(1.0e-7_dp,  0.0_dp, wp_sh)
        end if
        call start_state(b)
        call aeros_timestep_set_phis(ts, phis)
        do n = 1, nstep
            call aeros_timestep_step(ts, pool, vg, grd, b)
        end do
        call aeros_timestep_end(ts)

        worst = 0.0_dp
        do k = 1, nlev
            do lm = 1, nlm
                worst = max(worst, abs(a%spec%vor(lm,k)  - b%spec%vor(lm,k)))
                worst = max(worst, abs(a%spec%div(lm,k)  - b%spec%div(lm,k)))
                worst = max(worst, abs(a%spec%temp(lm,k) - b%spec%temp(lm,k)))
            end do
        end do
        do lm = 1, nlm
            worst = max(worst, abs(a%spec%lnps(lm) - b%spec%lnps(lm)))
        end do

        call aeros_state_end(a)
        call aeros_state_end(b)

        return

    end subroutine integrate_diff

    subroutine default_par(par)

        implicit none

        type(aeros_param_class), intent(out) :: par

        par%trunc    = trunc
        par%nlon     = -1
        par%nlat     = -1
        par%nlev     = nlev
        par%nthreads = -1

        par%dt            = 1800.0_wp
        par%semi_implicit = .TRUE.
        par%held_suarez   = .FALSE.
        par%eps_filter    = 0.06_wp
        par%raw_alpha     = 0.53_wp
        par%ndiff         = 6
        par%tau_diff      = 6.0_wp
        par%mass_fixer    = .FALSE.

        return

    end subroutine default_par

    subroutine start_state(now)
        ! A mildly moving isothermal atmosphere -- enough that the dynamics are
        ! doing something for the correction to be distinguishable from.

        implicit none

        type(aeros_state_class), intent(inout) :: now

        type(aeros_sht_class), pointer :: s
        real(wp), parameter :: tiso = 280.0_wp
        real(dp) :: amp, y00
        integer  :: k, lm, l, m

        s => pool%sht(1)
        y00 = sqrt(4.0_dp*acos(-1.0_dp))

        call aeros_spec_zero(now%spec)

        now%spec%lnps(aeros_sht_lm(s,0,0)) = cmplx(log(real(p0,dp))*y00, 0.0_dp, wp_sh)
        now%spec%lnps(aeros_sht_lm(s,2,1)) = cmplx(-0.03_dp, 0.01_dp, wp_sh)

        do k = 1, nlev
            now%spec%temp(aeros_sht_lm(s,0,0),k) = cmplx(real(tiso,dp)*y00, 0.0_dp, wp_sh)
            do lm = 1, s%nlm
                l = s%l_of_lm(lm); m = s%m_of_lm(lm)
                if (l < 1 .or. l > 6) cycle
                amp = 1.0e-5_dp/real(l*l, dp)
                now%spec%vor(lm,k) = cmplx(amp*cos(real(3*l+m,dp)), amp*sin(real(l+2*m,dp)), wp_sh)
                if (m == 0) now%spec%vor(lm,k) = cmplx(real(now%spec%vor(lm,k)), 0.0_wp_sh, wp_sh)
            end do
        end do

        return

    end subroutine start_state

    subroutine check(ok, label, nfail)

        implicit none

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

end program test_correction
