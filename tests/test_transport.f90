program test_transport
    ! Acceptance test for the positive-definite humidity transport.
    !
    ! There is no physics here yet -- condensation is the next commit -- so what
    ! has to be nailed down is that the TRANSPORT has the three properties the
    ! whole finite-volume choice was made for, and that it has them as
    ! machine-precision facts rather than approximately. A tracer scheme that is
    ! merely close on each of these is a tracer scheme that will, over a 10^5 yr
    ! integration, invent or destroy water.
    !
    ! 1. CONSTANCY. q = const stays const, to machine precision, under a wind
    !    with real horizontal divergence -- so the vertical mass flux is
    !    exercised, not just the horizontal sweep. This is the property that
    !    fails first if the tracer and the air mass are advanced by even
    !    slightly different fluxes, and it is the reason the module runs its own
    !    finite-volume air budget rather than borrowing the spectral one.
    !
    ! 2. CONSERVATION. Total water is conserved to machine precision under
    !    solid-body rotation over many steps. Flux form with single-valued
    !    faces makes this exact; the test is that it is exact and not just
    !    small.
    !
    ! 3. POSITIVITY. A positive blob with sharp edges stays >= 0 everywhere, at
    !    every step, as it is advected -- including where the leading edge
    !    steepens. This is the property spectral advection cannot deliver and
    !    the entire reason humidity is off the transform.
    !
    ! 4. THE POLAR COURANT SUB-STEP FIRES, AND CHANGES NOTHING IT SHOULDN'T. A
    !    wind that does not vanish at the poles drives the zonal Courant number
    !    above one there; the transport must sub-step to stay stable, and must
    !    still conserve and stay positive when it does. A scheme that needed
    !    sub-stepping and did not have it would simply go unstable near the pole
    !    -- silently, since the rest of the grid looks fine.
    !
    ! Exits non-zero on failure.

    use aeros_defs,     only : dp, wp, p0, pi, degrees_to_radians, &
                                aeros_grid_class
    use aeros_spectral
    use aeros_grid
    use aeros_vertical
    use aeros_transport

    implicit none

    integer, parameter :: trunc = 21
    integer, parameter :: nlev  = 8

    type(aeros_sht_pool_class), target :: pool
    type(aeros_grid_class)  :: grd
    type(aeros_vgrid_class) :: vg
    type(aeros_transport_class) :: mst

    real(wp), allocatable :: lnps(:,:), u(:,:,:), v(:,:,:), q(:,:,:)
    integer :: nlon, nlat, nfail

    nfail = 0

    call aeros_sht_pool_init(pool, trunc, quick=.TRUE.)
    call aeros_grid_init(grd, pool%sht(1))
    call aeros_vgrid_init(vg, nlev)
    call aeros_transport_init(mst, grd, nlev)

    nlon = grd%nlon; nlat = grd%nlat

    write(*,"(a,i0,a,i0,a,i0,a,i0)") " test_transport:: T", trunc, " L", nlev, &
                                        "  grid ", nlon, "x", nlat

    allocate(lnps(nlon,nlat), u(nlon,nlat,nlev), v(nlon,nlat,nlev), q(nlon,nlat,nlev))

    ! A surface pressure that varies over the grid, so the layer masses the
    ! transport diagnoses are not uniform and the mass weighting is genuinely
    ! tested.
    call set_lnps(grd, lnps)

    call test_constancy(nfail)
    call test_conservation(nfail)
    call test_positivity(nfail)
    call test_polar_substep(nfail)
    call test_accuracy(nfail)

    call aeros_transport_end(mst)
    call aeros_vgrid_end(vg)
    call aeros_grid_end(grd)
    call aeros_sht_pool_end(pool)

    write(*,*) ""
    if (nfail == 0) then
        write(*,*) "test_transport:: PASS"
    else
        write(*,*) "test_transport:: FAIL, ", nfail, " check(s) failed"
        stop 1
    end if

contains

    subroutine set_lnps(grd, lnps)
        implicit none
        type(aeros_grid_class), intent(in)  :: grd
        real(wp),               intent(out) :: lnps(:,:)
        real(wp) :: lon, lat
        integer  :: i, j
        do j = 1, grd%nlat
            lat = grd%lat(j)*degrees_to_radians
            do i = 1, grd%nlon
                lon = grd%lon(i)*degrees_to_radians
                lnps(i,j) = log(real(p0,wp)*(1.0_wp + 0.05_wp*sin(lon)*cos(lat)))
            end do
        end do
        return
    end subroutine set_lnps

    subroutine solid_body(grd, u, v, u0)
        ! u = u0 cos(lat), v = 0: rigid rotation. Non-divergent, and the zonal
        ! Courant number u/(a cos dlambda) is uniform in latitude -- cos cancels
        ! -- so this case does NOT trigger the polar sub-step. That is what
        ! makes it the clean case for conservation.
        implicit none
        type(aeros_grid_class), intent(in)  :: grd
        real(wp),               intent(out) :: u(:,:,:), v(:,:,:)
        real(wp),               intent(in)  :: u0
        real(wp) :: lat
        integer  :: i, j, k
        do k = 1, nlev
            do j = 1, grd%nlat
                lat = grd%lat(j)*degrees_to_radians
                do i = 1, grd%nlon
                    u(i,j,k) = u0*cos(lat)
                    v(i,j,k) = 0.0_wp
                end do
            end do
        end do
        return
    end subroutine solid_body

    subroutine divergent_wind(grd, u, v, u0)
        ! u = u0 sin(lon) cos(lat), v = u0 cos(lon): a wind with genuine
        ! horizontal divergence, so the column mass convergence -- and hence the
        ! vertical mass flux -- is non-zero. The constancy test needs this: a
        ! non-divergent wind would leave the vertical flux path untested.
        implicit none
        type(aeros_grid_class), intent(in)  :: grd
        real(wp),               intent(out) :: u(:,:,:), v(:,:,:)
        real(wp),               intent(in)  :: u0
        real(wp) :: lon, lat
        integer  :: i, j, k
        do k = 1, nlev
            do j = 1, grd%nlat
                lat = grd%lat(j)*degrees_to_radians
                do i = 1, grd%nlon
                    lon = grd%lon(i)*degrees_to_radians
                    u(i,j,k) = u0*sin(lon)*cos(lat)*real(k,wp)/real(nlev,wp)
                    v(i,j,k) = u0*cos(lon)*0.5_wp
                end do
            end do
        end do
        return
    end subroutine divergent_wind

    subroutine blob(grd, q, qbg, qamp)
        ! A positive background plus a localized bump near the equator. The bump
        ! has sharp flanks, which is where a non-monotone scheme would ring
        ! negative.
        implicit none
        type(aeros_grid_class), intent(in)  :: grd
        real(wp),               intent(out) :: q(:,:,:)
        real(wp),               intent(in)  :: qbg, qamp
        real(wp) :: lon, lat, r2
        integer  :: i, j, k
        do k = 1, nlev
            do j = 1, grd%nlat
                lat = grd%lat(j)*degrees_to_radians
                do i = 1, grd%nlon
                    lon = grd%lon(i)*degrees_to_radians
                    r2 = (lon - pi)**2 + (lat**2)*4.0_wp
                    q(i,j,k) = qbg + qamp*exp(-r2*6.0_wp)
                end do
            end do
        end do
        return
    end subroutine blob

    ! === 1. Constancy ========================================================

    subroutine test_constancy(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp), parameter :: qc = 0.7_wp
        real(wp) :: worst
        integer  :: n

        write(*,*) ""
        write(*,*) " -- constancy: q = 0.7 under a divergent wind, 50 steps"

        call divergent_wind(grd, u, v, 40.0_wp)
        q = qc

        do n = 1, 50
            call aeros_transport_transport(mst, vg, u, v, lnps, q, 1800.0_wp)
        end do

        worst = maxval(abs(q - qc))
        write(*,"(a40,es12.3)") "   max |q - 0.7| after 50 steps  ", worst
        write(*,"(a40,i8)")     "   sub-steps on the last call    ", mst%last_nsub

        call check(worst < 1.0e-13_wp, &
                    "a uniform q stays uniform under a divergent wind", nfail)
        return
    end subroutine test_constancy

    ! === 2. Conservation =====================================================

    subroutine test_conservation(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(dp) :: w0, w1, rel
        integer  :: n

        write(*,*) ""
        write(*,*) " -- conservation: a blob under solid-body rotation, 100 steps"

        ! Uniform surface pressure here: the scheme conserves the tracer mass
        ! against its OWN finite-volume air budget, so a clean conservation test
        ! needs the diagnosed air mass to equal that budget -- which it does
        ! only when ps is uniform (a varying ps with the wind held fixed makes
        ! div(dp V) non-zero, i.e. an air-mass source the fixed ps denies, and
        ! that inconsistency, not the scheme, would break the check). With ps
        ! co-evolving in the real model the two agree to O(truncation); that is
        ! measured at the integration level, not here.
        lnps = log(real(p0,wp))
        call solid_body(grd, u, v, 60.0_wp)
        call blob(grd, q, 0.001_wp, 0.02_wp)

        w0 = aeros_transport_water(mst, grd, vg, exp(lnps), q)
        do n = 1, 100
            call aeros_transport_transport(mst, vg, u, v, lnps, q, 1800.0_wp)
        end do
        w1 = aeros_transport_water(mst, grd, vg, exp(lnps), q)

        rel = abs(w1 - w0)/w0
        write(*,"(a40,es12.3,a)") "   initial water                 ", w0, " kg"
        write(*,"(a40,es12.3)")   "   |dW/W| over 100 steps         ", rel

        call check(rel < 1.0e-13_dp, "total water is conserved to machine precision", nfail)
        return
    end subroutine test_conservation

    ! === 3. Positivity =======================================================

    subroutine test_positivity(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(wp) :: qmin_ever, qmin
        integer  :: n

        write(*,*) ""
        write(*,*) " -- positivity: a blob with sharp flanks, 200 steps"

        lnps = log(real(p0,wp))
        call solid_body(grd, u, v, 70.0_wp)
        call blob(grd, q, 0.0_wp, 0.02_wp)      ! zero background: the flanks ARE the test

        qmin_ever = minval(q)
        do n = 1, 200
            call aeros_transport_transport(mst, vg, u, v, lnps, q, 1800.0_wp)
            qmin = minval(q)
            qmin_ever = min(qmin_ever, qmin)
        end do

        write(*,"(a40,es12.3)") "   min q over 200 steps          ", qmin_ever
        call check(qmin_ever >= 0.0_wp, "q never goes negative", nfail)
        return
    end subroutine test_positivity

    ! === 4. Polar sub-step ===================================================

    subroutine test_polar_substep(nfail)
        implicit none
        integer, intent(inout) :: nfail
        real(dp) :: w0, w1, rel
        real(wp) :: qmin_ever
        integer  :: n

        write(*,*) ""
        write(*,*) " -- polar Courant: a rigid u = 80 m/s that does not vanish at the pole"

        ! u constant in latitude (NOT u0 cos lat): the zonal cell shrinks toward
        ! the pole while u does not, so the Courant number there is large and
        ! the sub-step must engage.
        lnps = log(real(p0,wp))
        u = 80.0_wp
        v = 0.0_wp
        call blob(grd, q, 0.0_wp, 0.02_wp)

        w0 = aeros_transport_water(mst, grd, vg, exp(lnps), q)
        qmin_ever = minval(q)
        do n = 1, 50
            call aeros_transport_transport(mst, vg, u, v, lnps, q, 1800.0_wp)
            qmin_ever = min(qmin_ever, minval(q))
        end do
        w1 = aeros_transport_water(mst, grd, vg, exp(lnps), q)
        rel = abs(w1 - w0)/w0

        write(*,"(a40,i8)")       "   sub-steps forced              ", mst%last_nsub
        write(*,"(a40,es12.3)")   "   max Courant seen              ", mst%last_cfl
        write(*,"(a40,es12.3)")   "   |dW/W|                        ", rel
        write(*,"(a40,es12.3)")   "   min q                         ", qmin_ever

        call check(mst%last_nsub > 1, "the polar Courant number forced sub-stepping", nfail)
        call check(rel < 1.0e-13_dp, "and water is still conserved through it", nfail)
        call check(qmin_ever >= 0.0_wp, "and q is still positive through it", nfail)
        return
    end subroutine test_polar_substep

    ! === 5. Accuracy: the limiter is second order, not upwind ================

    subroutine test_accuracy(nfail)
        ! Advect a broad, smooth, well-resolved bell one full solid-body
        ! rotation and measure how much of its peak amplitude survives. A full
        ! rotation returns the exact field to its start, so any loss is
        ! numerical diffusion. First-order upwind smears a bell like this well
        ! below half amplitude over a rotation; the van Leer limiter must keep
        ! most of it -- that difference is the whole point of the limiter, and
        ! the reason the moist water budget (m2_results.md §9) improves with it.

        implicit none
        integer, intent(inout) :: nfail

        real(wp) :: u0, peak0, peak1, retain
        integer  :: nrot, n

        write(*,*) ""
        write(*,*) " -- accuracy: peak retention over one full rotation of a smooth bell"

        lnps = log(real(p0,wp))
        u0 = 60.0_wp
        call solid_body(grd, u, v, u0)
        call smooth_bell(grd, q)

        peak0 = maxval(q)

        ! One full rotation: angular velocity u0/a, period 2*pi*a/u0.
        nrot = nint(2.0_wp*real(pi,wp)*real(6.371e6_wp,wp)/u0/1800.0_wp)
        do n = 1, nrot
            call aeros_transport_transport(mst, vg, u, v, lnps, q, 1800.0_wp)
        end do

        peak1  = maxval(q)
        retain = peak1/peak0

        write(*,"(a40,i8)")     "   steps per rotation             ", nrot
        write(*,"(a40,f10.3)")  "   peak amplitude retained        ", retain

        call check(retain > 0.70_wp, &
                    "van Leer keeps most of the peak (upwind would not)", nfail)
        call check(minval(q) >= 0.0_wp, "and stays positive", nfail)
        return
    end subroutine test_accuracy

    subroutine smooth_bell(grd, q)
        ! A broad cosine-squared bell centred on the equator, wide enough to be
        ! well resolved at T21 -- so the loss measured is the scheme's, not the
        ! grid's inability to represent the shape.
        implicit none
        type(aeros_grid_class), intent(in)  :: grd
        real(wp),               intent(out) :: q(:,:,:)
        real(wp) :: lon, lat, r
        integer  :: i, j, k
        do k = 1, nlev
            do j = 1, grd%nlat
                lat = grd%lat(j)*degrees_to_radians
                do i = 1, grd%nlon
                    lon = grd%lon(i)*degrees_to_radians
                    r = sqrt((lon - pi)**2 + lat**2)
                    if (r < 1.2_wp) then
                        q(i,j,k) = 0.01_wp*cos(r/1.2_wp*0.5_wp*pi)**2
                    else
                        q(i,j,k) = 0.0_wp
                    end if
                end do
            end do
        end do
        return
    end subroutine smooth_bell

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

end program test_transport
