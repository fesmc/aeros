program test_tendency
    ! Acceptance test for the primitive-equation right-hand sides.
    !
    ! There is no closed-form solution to compare a nonlinear tendency against,
    ! so the tests are BALANCE tests: states whose exact tendency is known to be
    ! zero, and integral identities that must hold whatever the state.
    !
    ! The central one is the ISOTHERMAL ATMOSPHERE AT REST OVER TOPOGRAPHY.
    ! Simmons & Burridge's discretization is constructed so that this state has
    ! IDENTICALLY zero tendencies -- the discrete pressure-gradient force
    ! cancels the discrete geopotential gradient term for term, not
    ! approximately. One test exercises the hydrostatic integral, the shared
    ! coefficient c_k, the omega term, the vertical mass flux and the spectral
    ! Laplacian at once, and it fails for any of the plausible ways of getting
    ! the vertical discretization subtly wrong.
    !
    ! For an isothermal resting atmosphere in hydrostatic balance,
    !
    !     ln p_s = ln p_00 - Phi_s/(R T)
    !
    ! and the balance the discretization must reproduce is
    !
    !     c_k = 1 - p_s d/dp_s [ sum_(j>k) dlnp_j + alpha_k ]
    !
    ! which is what alpha_k is FOR. In pure sigma both sides are 1 because
    ! dlnp and alpha depend on sigma alone. In a hybrid coordinate both sides
    ! vary, and they still agree -- that is Simmons & Burridge's result, and it
    ! holds here to 3e-12.
    !
    ! THE ONE EXCEPTION, and it is measured rather than assumed: a model top at
    ! ZERO pressure. There Simmons & Burridge impose alpha_1 = ln 2 rather than
    ! deriving it -- the layer-mean geopotential of a layer reaching p = 0
    ! diverges, so a finite convention is needed -- and the identity above
    ! fails for that layer alone. Measured, the residual is confined to k = 1
    ! (2.4e-10 against 1e-22 in every other layer). It is a property of the
    ! zero-pressure top, NOT of pure sigma: aeros' default p_top = 10 hPa
    ! balances in every layer.
    !
    ! Practical consequence for M1.5: Held-Suarez is nominally defined on
    ! sigma from 0 to 1, so a literal implementation inherits a top layer that
    ! is not hydrostatically consistent. Use a small non-zero p_top there, or
    ! expect noise in the top layer.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, wp_sh, R_d, p0, grav, aeros_grid_class, &
                                aeros_state_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_vordiv
    use aeros_state
    use aeros_tendency

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 12

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)     :: grd
    type(aeros_vgrid_class)    :: vg
    type(aeros_state_class)    :: now
    type(aeros_work_class)     :: wrk
    type(aeros_tend_class)     :: tnd

    integer :: nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_tendency:: T", trunc, " L", nlev, &
                                        "  grid ", grd%nlon, "x", grd%nlat

    ! === 1. Rest state, flat surface, pure sigma =============================
    call aeros_vgrid_init(vg, nlev, sigma_t=0.0_wp, p_top=0.0_wp)
    call setup(pool, grd, vg, now, wrk, tnd)
    call rest_state(pool, grd, vg, now, wrk, topography=.FALSE.)
    call check_rest(pool, grd, vg, now, wrk, tnd, &
                        "flat, pure sigma", 1.0e-12_dp, 1, nfail)
    call teardown(now, wrk, tnd)
    call aeros_vgrid_end(vg)

    ! === 2. Rest state, WITH topography, hybrid (the shipped default) ========
    ! The demanding one: grad(ln p_s) is now large, so every term is active,
    ! and with a non-zero model top EVERY layer must cancel exactly.
    call aeros_vgrid_init(vg, nlev)
    call setup(pool, grd, vg, now, wrk, tnd)
    call rest_state(pool, grd, vg, now, wrk, topography=.TRUE.)
    call check_rest(pool, grd, vg, now, wrk, tnd, &
                        "topography, hybrid (default)", 1.0e-9_dp, 1, nfail)
    call teardown(now, wrk, tnd)
    call aeros_vgrid_end(vg)

    ! === 3. Rest state, WITH topography, pure sigma to a ZERO-pressure top ===
    ! Exact below the top layer; the top layer carries the alpha_1 = ln 2
    ! convention and does not balance. Asserted from k = 2, with k = 1 reported
    ! so that a regression there is still visible.
    call aeros_vgrid_init(vg, nlev, sigma_t=0.0_wp, p_top=0.0_wp)
    call setup(pool, grd, vg, now, wrk, tnd)
    call rest_state(pool, grd, vg, now, wrk, topography=.TRUE.)
    call check_rest(pool, grd, vg, now, wrk, tnd, &
                        "topography, pure sigma, p_top = 0", 1.0e-9_dp, 2, nfail)
    call teardown(now, wrk, tnd)
    call aeros_vgrid_end(vg)

    ! === 4. A moving state: integral identities ==============================
    call aeros_vgrid_init(vg, nlev, sigma_t=0.0_wp, p_top=0.0_wp)
    call setup(pool, grd, vg, now, wrk, tnd)
    call test_moving(pool, grd, vg, now, wrk, tnd, nfail)
    call teardown(now, wrk, tnd)
    call aeros_vgrid_end(vg)

    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_tendency:: PASS"
    else
        write(*,*) "test_tendency:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine setup(pool, grd, vg, now, wrk, tnd)

        implicit none

        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_grid_class),  intent(in)    :: grd
        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_state_class), intent(inout) :: now
        type(aeros_work_class),  intent(inout) :: wrk
        type(aeros_tend_class),  intent(inout) :: tnd

        call aeros_state_alloc(now, grd, pool%sht(1)%nlm, vg%nlev)
        call aeros_work_alloc(wrk, grd%nlon, grd%nlat, vg%nlev)
        call aeros_tend_alloc(tnd, pool%sht(1)%nlm, vg%nlev)

        return

    end subroutine setup

    subroutine teardown(now, wrk, tnd)

        implicit none

        type(aeros_state_class), intent(inout) :: now
        type(aeros_work_class),  intent(inout) :: wrk
        type(aeros_tend_class),  intent(inout) :: tnd

        call aeros_state_end(now)
        call aeros_work_end(wrk)
        call aeros_tend_end(tnd)

        return

    end subroutine teardown

    subroutine rest_state(pool, grd, vg, now, wrk, topography)
        ! Build an isothermal atmosphere at rest, in exact discrete hydrostatic
        ! balance with its own surface geopotential.
        !
        ! Order matters. ln p_s is defined SPECTRALLY and synthesized to the
        ! grid; Phi_s is then computed FROM that grid field. Doing it the other
        ! way round -- writing an analytic Phi_s and analyzing the implied
        ! ln p_s -- would leave a truncation residual between the two, and the
        ! test would measure that residual instead of the discretization.

        implicit none

        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_grid_class),  intent(in)    :: grd
        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_state_class), intent(inout) :: now
        type(aeros_work_class),  intent(inout) :: wrk
        logical, intent(in) :: topography

        type(aeros_sht_class), pointer :: s
        real(wp), parameter :: tiso = 280.0_wp
        integer  :: k, i, j, lm

        s => pool%sht(1)

        now%vor  = (0.0_wp_sh, 0.0_wp_sh)
        now%div  = (0.0_wp_sh, 0.0_wp_sh)
        now%lnps = (0.0_wp_sh, 0.0_wp_sh)

        ! Isothermal: only the (0,0) coefficient. With orthonormal harmonics
        ! Y_00 = 1/sqrt(4 pi), so the coefficient is T*sqrt(4 pi).
        now%temp = (0.0_wp_sh, 0.0_wp_sh)
        do k = 1, vg%nlev
            now%temp(aeros_sht_lm(s,0,0),k) = &
                    cmplx(real(tiso,dp)*sqrt(4.0_dp*acos(-1.0_dp)), 0.0_dp, wp_sh)
        end do

        ! Surface pressure: a uniform part plus, optionally, structure.
        now%lnps(aeros_sht_lm(s,0,0)) = &
                cmplx(log(real(p0,dp))*sqrt(4.0_dp*acos(-1.0_dp)), 0.0_dp, wp_sh)

        if (topography) then
            ! Enough structure to make grad(ln p_s) genuinely large: this
            ! corresponds to surface elevations of order a kilometre.
            lm = aeros_sht_lm(s,2,1); now%lnps(lm) = cmplx(-0.06_dp,  0.03_dp, wp_sh)
            lm = aeros_sht_lm(s,3,0); now%lnps(lm) = cmplx( 0.05_dp,  0.0_dp,  wp_sh)
            lm = aeros_sht_lm(s,5,3); now%lnps(lm) = cmplx( 0.02_dp, -0.01_dp, wp_sh)
        end if

        ! Grid surface pressure, then the surface geopotential that balances it.
        call aeros_sht_synthesis(s, now%lnps, wrk%lnps_g)

        do j = 1, grd%nlat
            do i = 1, grd%nlon
                wrk%phis(i,j) = R_d*tiso*(log(p0) - wrk%lnps_g(i,j))
            end do
        end do

        return

    end subroutine rest_state

    subroutine check_rest(pool, grd, vg, now, wrk, tnd, label, tol, kmin, nfail)
        ! All four tendencies of a resting balanced state, measured against the
        ! natural scale of the terms that cancelled to produce them.

        implicit none

        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_grid_class),  intent(in)    :: grd
        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_state_class), intent(inout) :: now
        type(aeros_work_class),  intent(inout) :: wrk
        type(aeros_tend_class),  intent(inout) :: tnd
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: tol
        integer,  intent(in) :: kmin      ! first level included in the assertion
        integer,  intent(inout) :: nfail

        real(dp) :: scale, e_vor, e_div, e_t, e_ps, e_top, dscale
        real(dp) :: elev_max
        integer  :: i, j, k

        call aeros_tendency_calc(pool, vg, grd, now, wrk, tnd)

        ! Scale: the magnitude of the pressure-gradient term that had to cancel.
        ! Normalising by the tendency itself would be meaningless (it is ~0);
        ! normalising by the terms that produced it is the honest measure.
        scale = 0.0_dp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                scale = max(scale, abs(real(wrk%ek(i,j,vg%nlev), dp)))
            end do
        end do
        scale = max(scale, 1.0_dp)

        elev_max = 0.0_dp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                elev_max = max(elev_max, abs(real(wrk%phis(i,j), dp))/grav)
            end do
        end do

        ! d/dt has units of [X]/s; compare against [scale]/a^2 for the momentum
        ! tendencies, which is what laplacian(Phi) contributes.
        dscale = scale/(6.371e6_dp**2)

        e_vor = 0.0_dp; e_div = 0.0_dp; e_t = 0.0_dp
        do k = kmin, vg%nlev
            e_vor = max(e_vor, maxval(abs(tnd%vor(:,k))))
            e_div = max(e_div, maxval(abs(tnd%div(:,k))))
            e_t   = max(e_t,   maxval(abs(tnd%temp(:,k))))
        end do
        e_vor = e_vor/dscale
        e_div = e_div/dscale
        e_t   = e_t/280.0_dp
        e_ps  = maxval(abs(tnd%lnps))

        e_top = maxval(abs(tnd%div(:,1)))/dscale

        write(*,*) ""
        write(*,"(a,a)") " -- rest state: ", trim(label)
        write(*,"(a40,f10.1,a)") "   max surface elevation ", elev_max, " m"
        write(*,"(a40,es12.3)")  "   |d zeta/dt| (scaled) ", e_vor
        write(*,"(a40,es12.3)")  "   |d D/dt|    (scaled) ", e_div
        write(*,"(a40,es12.3)")  "   |d T/dt|    (scaled) ", e_t
        write(*,"(a40,es12.3)")  "   |d lnps/dt|          ", e_ps
        write(*,"(a40,es12.3)")  "   top layer |d D/dt|   ", e_top
        if (kmin > 1) write(*,"(a,i0,a)") "   (asserted from level ", kmin, ")"

        call check(e_vor < tol, "vorticity tendency vanishes",  nfail)
        call check(e_div < tol, "divergence tendency vanishes", nfail)
        call check(e_t   < tol, "temperature tendency vanishes", nfail)
        call check(e_ps  < tol, "surface pressure tendency vanishes", nfail)

        return

    end subroutine check_rest

    subroutine test_moving(pool, grd, vg, now, wrk, tnd, nfail)
        ! Identities that must hold for ANY state, checked on a moving one.

        implicit none

        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_grid_class),  intent(in)    :: grd
        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_state_class), intent(inout) :: now
        type(aeros_work_class),  intent(inout) :: wrk
        type(aeros_tend_class),  intent(inout) :: tnd
        integer, intent(inout) :: nfail

        type(aeros_sht_class), pointer :: s
        real(wp), parameter :: tiso = 280.0_wp
        real(dp) :: amp, mass_tend, mass_scale
        integer  :: k, lm, l, m, i, j
        logical  :: ok

        s => pool%sht(1)

        write(*,*) ""
        write(*,*) " -- moving state"

        ! A sheared, rotational flow with realistic magnitudes.
        now%vor  = (0.0_wp_sh, 0.0_wp_sh)
        now%div  = (0.0_wp_sh, 0.0_wp_sh)
        now%temp = (0.0_wp_sh, 0.0_wp_sh)
        now%lnps = (0.0_wp_sh, 0.0_wp_sh)

        now%lnps(aeros_sht_lm(s,0,0)) = &
                cmplx(log(real(p0,dp))*sqrt(4.0_dp*acos(-1.0_dp)), 0.0_dp, wp_sh)
        now%lnps(aeros_sht_lm(s,2,1)) = cmplx(-0.03_dp, 0.01_dp, wp_sh)

        do k = 1, vg%nlev
            now%temp(aeros_sht_lm(s,0,0),k) = &
                    cmplx(real(tiso,dp)*sqrt(4.0_dp*acos(-1.0_dp)), 0.0_dp, wp_sh)
            now%temp(aeros_sht_lm(s,2,0),k) = cmplx(-20.0_dp*real(k,dp)/real(vg%nlev,dp), &
                                                        0.0_dp, wp_sh)
            do lm = 1, s%nlm
                l = s%l_of_lm(lm); m = s%m_of_lm(lm)
                if (l < 1 .or. l > 6) cycle
                amp = 1.0e-5_dp/real(l*l, dp)*(0.5_dp + 0.5_dp*real(k,dp)/real(vg%nlev,dp))
                now%vor(lm,k) = cmplx(amp*cos(real(3*l+m,dp)), amp*sin(real(l+2*m,dp)), wp_sh)
                now%div(lm,k) = cmplx(0.2_dp*amp*sin(real(l+m,dp)), &
                                        0.2_dp*amp*cos(real(2*l+m,dp)), wp_sh)
                if (m == 0) then
                    now%vor(lm,k) = cmplx(real(now%vor(lm,k)), 0.0_wp_sh, wp_sh)
                    now%div(lm,k) = cmplx(real(now%div(lm,k)), 0.0_wp_sh, wp_sh)
                end if
            end do
        end do

        wrk%phis = 0.0_wp

        call aeros_tendency_calc(pool, vg, grd, now, wrk, tnd)

        ! Winds should be plausible -- if they are not, nothing below means
        ! anything.
        write(*,"(a40,f9.2,a)") "   max |u| ", maxval(abs(wrk%u)), " m s-1"
        call check(maxval(abs(wrk%u)) > 1.0_wp .and. maxval(abs(wrk%u)) < 150.0_wp, &
                    "winds are physically plausible", nfail)

        ! Tendencies must be finite and of a sane magnitude. A vorticity
        ! tendency of 1e-9 s-2 changes zeta by 1e-5 s-1 in ~3 hours, which is
        ! the right order for synoptic flow.
        ok = all(tnd%vor == tnd%vor) .and. all(tnd%div == tnd%div) &
                .and. all(tnd%temp == tnd%temp) .and. all(tnd%lnps == tnd%lnps)
        call check(ok, "no NaNs in the tendencies", nfail)

        write(*,"(a40,es12.3,a)") "   max |d zeta/dt| ", maxval(abs(tnd%vor)), " s-2"
        write(*,"(a40,es12.3,a)") "   max |d T/dt|    ", maxval(abs(tnd%temp)), " K s-1"

        ! GLOBAL MASS. The surface-pressure tendency is minus the divergence of
        ! a vertically integrated mass flux, so its global integral must
        ! vanish: the model may move air around but may not create it. This is
        ! the one conservation statement available before there is a time
        ! integrator, and it is the one that matters most over 10^5 yr.
        mass_tend  = 0.0_dp
        mass_scale = 0.0_dp
        do j = 1, grd%nlat
            do i = 1, grd%nlon
                mass_tend  = mass_tend &
                        + real(wrk%dlnpsdt(i,j),dp)*exp(real(wrk%lnps_g(i,j),dp)) &
                            *real(grd%area(i,j),dp)
                mass_scale = mass_scale &
                        + abs(real(wrk%dlnpsdt(i,j),dp))*exp(real(wrk%lnps_g(i,j),dp)) &
                            *real(grd%area(i,j),dp)
            end do
        end do

        write(*,"(a40,es12.3)") "   global mass tendency / scale ", abs(mass_tend)/mass_scale
        call check(abs(mass_tend)/mass_scale < 1.0e-12_dp, &
                    "global mass tendency vanishes", nfail)

        return

    end subroutine test_moving

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

end program test_tendency
