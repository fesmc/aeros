program test_cloud_prog
    ! Acceptance test for the prognostic cloud-fraction scheme (aeros_cloud_prog,
    ! Sundqvist 1989). The scheme carries cloud fraction as a genuine gridpoint
    ! prognostic with a source/sink budget, so the checks are on the qualitative
    ! behaviour the budget must guarantee, not on tolerances:
    !
    ! 1. BOUNDED. From any cloud fraction in [0,1] and any humidity, the updated
    !    fraction stays in [0,1]. The relaxation form is a convex move over the
    !    step, so this holds without clipping -- the guard is only round-off.
    !
    ! 2. SATURATION DRIVES OVERCAST. A column held at (or above) saturation is
    !    pulled to cloud fraction ~1: formation with no evaporation sink.
    !
    ! 3. A DRY COLUMN CLEARS. A subsaturated column relaxes any pre-existing
    !    cloud to ~0: the evaporation sink the diagnostic scheme lacks.
    !
    ! 4. FORMS WHEN RH > RH_crit, DECAYS WHEN SUBSATURATED. The two signs of
    !    dC/dt, isolated.
    !
    ! 5. THE SINK KEEPS COVER BELOW THE DIAGNOSTIC. At a fixed sub-saturated RH
    !    above RH_crit the steady fraction sits strictly below the diagnostic
    !    Sundqvist value -- the property that stops the runaway.
    !
    ! 6. DISABLED IS A NO-OP. With the scheme off, cf is left exactly untouched
    !    -- the mechanism that makes l_cloud_prog=.FALSE. bit-for-bit identical
    !    to the diagnostic model (the full-suite bit-repro covers the rest).
    !
    ! Exits non-zero on failure.

    use aeros_defs,       only : dp, wp, p0, aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_condensation, only : aeros_qsat
    use aeros_cloud_prog

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)       :: grd
    type(aeros_vgrid_class)      :: vg
    type(aeros_cloud_prog_class) :: cpr

    real(wp), allocatable :: t_g(:,:,:), qv(:,:,:), cf(:,:,:), cf0(:,:,:), lnps(:,:)
    real(wp) :: dt
    integer  :: nlon, nlat, nfail

    nfail = 0
    dt = 1800.0_wp

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)
    call aeros_cloud_prog_init(cpr, grd, .TRUE.)     ! enabled, default knobs

    nlon = grd%nlon; nlat = grd%nlat

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_cloud_prog:: T", trunc, " L", nlev, &
                                        "  grid ", nlon, "x", nlat

    allocate(t_g(nlon,nlat,nlev), qv(nlon,nlat,nlev), cf(nlon,nlat,nlev), &
             cf0(nlon,nlat,nlev), lnps(nlon,nlat))
    lnps = log(real(p0,wp))

    call test_disabled_noop(nfail)
    call test_bounded(nfail)
    call test_saturated_to_one(nfail)
    call test_dry_to_zero(nfail)
    call test_form_and_decay(nfail)
    call test_sink_below_diagnostic(nfail)

    call aeros_cloud_prog_end(cpr)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_cloud_prog:: PASS"
    else
        write(*,*) "test_cloud_prog:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine set_temperature(t)
        ! A plausible tropospheric profile (as in test_condensation).
        real(wp), intent(out) :: t(:,:,:)
        real(wp) :: sigma
        integer  :: k
        do k = 1, nlev
            sigma = (real(k,wp) - 0.5_wp)/real(nlev,wp)
            t(:,:,k) = 300.0_wp - 45.0_wp*(1.0_wp - sigma)
        end do
        return
    end subroutine set_temperature

    subroutine set_rh(rh_target)
        ! Load qv so every level sits at the given relative humidity.
        real(wp), intent(in) :: rh_target
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), qs, dqsdt
        integer  :: i, j, k
        do j = 1, nlat
            do i = 1, nlon
                call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
                do k = 1, nlev
                    call aeros_qsat(t_g(i,j,k), pfull(k), qs, dqsdt)
                    qv(i,j,k) = rh_target*qs
                end do
            end do
        end do
        return
    end subroutine set_rh

    subroutine test_disabled_noop(nfail)
        integer, intent(inout) :: nfail
        type(aeros_cloud_prog_class) :: off
        write(*,*) ""
        write(*,*) " -- disabled scheme leaves cf exactly untouched"
        call aeros_cloud_prog_init(off, grd, .FALSE.)
        call set_temperature(t_g)
        call set_rh(1.2_wp)              ! supersaturated: would form if enabled
        cf  = 0.37_wp
        cf0 = cf
        call aeros_cloud_prog_apply(off, vg, t_g, qv, lnps, cf, dt)
        call check(all(cf == cf0), "cf is bit-for-bit unchanged when disabled", nfail)
        call aeros_cloud_prog_end(off)
        return
    end subroutine test_disabled_noop

    subroutine test_bounded(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), qs, dqsdt, rr
        logical  :: ok
        integer  :: i, j, k, n
        write(*,*) ""
        write(*,*) " -- cf stays in [0,1] over many steps of varied humidity"
        call set_temperature(t_g)
        ! A spatially varied humidity from sub- to super-saturated, and a varied
        ! initial cloud, to exercise both source and sink everywhere.
        do j = 1, nlat
            do i = 1, nlon
                call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
                do k = 1, nlev
                    call aeros_qsat(t_g(i,j,k), pfull(k), qs, dqsdt)
                    rr = 0.2_wp + 1.1_wp*real(mod(i*7 + j*3 + k*5, 13), wp)/12.0_wp
                    qv(i,j,k) = rr*qs
                    cf(i,j,k) = real(mod(i*3 + j*2 + k, 11), wp)/10.0_wp
                end do
            end do
        end do
        ok = .TRUE.
        do n = 1, 60
            call aeros_cloud_prog_apply(cpr, vg, t_g, qv, lnps, cf, dt)
            if (minval(cf) < 0.0_wp .or. maxval(cf) > 1.0_wp) ok = .FALSE.
            if (any(cf /= cf)) ok = .FALSE.
        end do
        write(*,"(a40,es12.4,a,es12.4)") "   min / max cf after 60 steps  ", &
            minval(cf), " / ", maxval(cf)
        call check(ok, "cf remains in [0,1] and finite at every step", nfail)
        return
    end subroutine test_bounded

    subroutine test_saturated_to_one(nfail)
        integer, intent(inout) :: nfail
        integer :: n
        write(*,*) ""
        write(*,*) " -- a saturated column is driven to overcast"
        call set_temperature(t_g)
        call set_rh(1.0_wp)              ! exactly saturated
        cf = 0.0_wp
        do n = 1, 200
            call aeros_cloud_prog_apply(cpr, vg, t_g, qv, lnps, cf, dt)
        end do
        write(*,"(a40,f10.5)") "   min cf after 200 steps        ", minval(cf)
        call check(minval(cf) > 0.99_wp, "saturated column: cf -> ~1", nfail)
        return
    end subroutine test_saturated_to_one

    subroutine test_dry_to_zero(nfail)
        integer, intent(inout) :: nfail
        integer :: n
        write(*,*) ""
        write(*,*) " -- a dry column clears any pre-existing cloud"
        call set_temperature(t_g)
        call set_rh(0.2_wp)              ! well below RH_crit
        cf = 0.85_wp
        do n = 1, 200
            call aeros_cloud_prog_apply(cpr, vg, t_g, qv, lnps, cf, dt)
        end do
        write(*,"(a40,es12.4)") "   max cf after 200 steps        ", maxval(cf)
        call check(maxval(cf) < 0.01_wp, "dry column: cf -> ~0", nfail)
        return
    end subroutine test_dry_to_zero

    subroutine test_form_and_decay(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: cf_after_form
        write(*,*) ""
        write(*,*) " -- forms above RH_crit, decays when subsaturated"
        ! Formation: RH above RH_crit everywhere, from zero cloud.
        call set_temperature(t_g)
        call set_rh(0.95_wp)
        cf = 0.0_wp
        call aeros_cloud_prog_apply(cpr, vg, t_g, qv, lnps, cf, dt)
        cf_after_form = minval(cf)
        write(*,"(a40,f10.5)") "   min cf after one forming step ", cf_after_form
        call check(minval(cf) > 0.0_wp, "cf forms where RH > RH_crit", nfail)

        ! Decay: RH below RH_crit, from a filled cloud; one step must reduce it.
        call set_rh(0.4_wp)
        cf  = 0.7_wp
        cf0 = cf
        call aeros_cloud_prog_apply(cpr, vg, t_g, qv, lnps, cf, dt)
        write(*,"(a40,f10.5)") "   cf after one decaying step    ", maxval(cf)
        call check(maxval(cf) < maxval(cf0), "cf decays where subsaturated", nfail)
        return
    end subroutine test_form_and_decay

    subroutine test_sink_below_diagnostic(nfail)
        integer, intent(inout) :: nfail
        real(wp) :: phalf(0:nlev), pfull(nlev), dpc(nlev), qs, dqsdt
        real(wp) :: rh, rhc, sig, cf_eq, b, cf_steady
        integer  :: n, k
        write(*,*) ""
        write(*,*) " -- the evaporation sink holds cover below the diagnostic value"
        ! Fixed RH above RH_crit but subsaturated: run to steady state.
        call set_temperature(t_g)
        call set_rh(0.95_wp)
        cf = 0.0_wp
        do n = 1, 400
            call aeros_cloud_prog_apply(cpr, vg, t_g, qv, lnps, cf, dt)
        end do
        ! Diagnostic Sundqvist value at a mid-level reference column (1,1,k).
        k = nlev/2
        call aeros_vgrid_pressure(vg, real(p0,wp), phalf, pfull, dpc)
        call aeros_qsat(t_g(1,1,k), pfull(k), qs, dqsdt)
        rh  = qv(1,1,k)/qs
        sig = pfull(k)/real(p0,wp)
        rhc = cpr%rhc_top + (cpr%rhc_sfc - cpr%rhc_top)*min(1.0_wp,max(0.0_wp,sig))
        b   = (1.0_wp - rh)/(1.0_wp - rhc)
        cf_eq = 1.0_wp - sqrt(max(0.0_wp, b))
        cf_steady = cf(1,1,k)
        write(*,"(a40,f10.5)") "   steady prognostic cf          ", cf_steady
        write(*,"(a40,f10.5)") "   diagnostic Sundqvist cf       ", cf_eq
        call check(cf_steady < cf_eq .and. cf_steady > 0.0_wp, &
                    "steady cover is positive but below the diagnostic value", nfail)
        return
    end subroutine test_sink_below_diagnostic

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

end program test_cloud_prog
