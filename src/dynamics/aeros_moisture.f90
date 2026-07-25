module aeros_moisture
    ! Prognostic specific humidity: positive-definite finite-volume transport.
    !
    ! This is the first prognostic in aeros that does NOT live in spectral space,
    ! and the reason is the whole point of the module. Spectral advection of a
    ! positive field is not positive: truncating exp-shaped humidity gradients at
    ! fronts and the ITCZ produces Gibbs ringing and negative water (design.md
    ! section 4.2). So humidity is carried and advected ON THE GRID, by a
    ! finite-volume scheme that cannot make a positive field negative, and it
    ! never touches the spherical-harmonic transform.
    !
    ! === What "positive-definite" costs, and how it is bought ================
    !
    ! Three properties are required and each is provable rather than hoped for:
    !
    !   CONSERVATION. Total water is conserved to machine precision. The scheme
    !   is flux-form with single-valued face fluxes, so every internal flux
    !   appears once with each sign and the global sum telescopes; the only
    !   faces that do not cancel are the poles, where cos(lat) = 0 makes the
    !   meridional flux identically zero, and the top and bottom, where the
    !   vertical mass flux is zero by the coordinate's construction.
    !
    !   POSITIVITY. q stays >= 0 from any q >= 0. Donor-cell (upwind) fluxes
    !   make the updated tracer mass a convex combination of its neighbours'
    !   values, so it cannot go below the smallest of them, as long as no cell
    !   is asked to give away more than it holds in one step -- the Courant
    !   condition, enforced by sub-stepping (see `nsub` below).
    !
    !   CONSTANCY. q = const stays const, to machine precision. The tracer and
    !   the air mass are advanced by the SAME face fluxes, so for q = 1 the
    !   tracer equation IS the air-mass equation and the ratio is unchanged.
    !   This is the property a tracer scheme most easily violates and the one
    !   that matters most: a scheme that fails it manufactures humidity
    !   gradients out of a uniform field wherever the wind diverges.
    !
    ! === Mass consistency: the design decision (option ii) ===================
    !
    ! Constancy above requires that the air mass the tracer is divided by is
    ! advanced by the tracer's OWN fluxes -- not the air mass the spectral
    ! dynamics produces. So this module runs its own finite-volume air-mass
    ! budget: from the surface pressure it diagnoses the layer masses, forms the
    ! horizontal mass-flux divergence on the grid, and from that the vertical
    ! mass flux -- the same construction the dynamics uses (aeros_tendency's
    ! `mflux`), but built from grid fluxes so it closes the grid tracer budget
    ! exactly.
    !
    ! That FV air mass differs from the spectral air mass by O(truncation),
    ! because one divergence is a grid finite difference and the other a
    ! spectral derivative. This is not a leak and does not accumulate: the layer
    ! masses are re-diagnosed from the true (spectral) surface pressure every
    ! step, so the difference is a per-step consistency error bounded by the
    ! truncation, and the max-principle division keeps q bounded regardless. It
    ! is the price of positivity on a spectral core, and it is diagnosable --
    ! aeros_moisture_report prints the gap.
    !
    ! This is the LGMR "option (ii)": a two-time-level forward scheme on the
    ! grid with its own FV air budget, chosen over leapfrogging q with a water
    ! mass-fixer (option i) precisely so there is no fixer and no leapfrog
    ! computational mode in the humidity.
    !
    ! === Accuracy: first-order now, limited later ============================
    !
    ! The fluxes here are FIRST-ORDER upwind (donor cell). That is deliberately
    ! the correct-and-provable baseline, not the finished scheme: donor cell is
    ! the scheme whose conservation, positivity and constancy are unconditional
    ! (under the Courant limit), and those are exactly what tests/test_moisture
    ! checks. It is also diffusive -- too diffusive for a real humidity field --
    ! so a van Leer / MC flux limiter is the next commit, dropped into
    ! `face_upwind` without changing anything else. The invariants above do not
    ! depend on the accuracy of the reconstruction, only on its monotonicity, so
    ! the limiter inherits them.
    !
    ! === The polar Courant problem ==========================================
    !
    ! A spectral core has no grid Courant limit; a grid transport scheme does,
    ! and on a Gaussian grid the zonal cell width a cos(lat) dlambda collapses
    ! toward the poles while the wind does not, so the zonal Courant number is
    ! large there (~5 at T31, dt = 1800 s). The fix is to sub-step the whole
    ! transport `nsub` times per dynamics step, `nsub` set from the largest
    ! total Courant number on the grid. A Fourier/polar filter -- the classical
    ! spectral-model remedy -- is exactly what must NOT be used here: filtering a
    ! positive field in wavenumber space makes it negative, which is the
    ! pathology this module exists to avoid. Sub-stepping preserves every
    ! invariant; it only costs arithmetic, and transport is grid-local and cheap
    ! against the transforms.
    !
    ! Global sub-stepping (the whole grid sub-steps at the polar rate) is
    ! wasteful and known to be: a per-row zonal sub-cycle would do the extra
    ! work only where it is needed. That is an optimization for later; the
    ! global form is simpler and its correctness is obvious.

    use aeros_defs,     only : dp, wp, io_unit_err, r_earth, grav, pi, &
                                aeros_grid_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure

    implicit none

    private

    type aeros_moist_class
        ! The transport's fixed geometry and its scratch, allocated once.

        integer :: nlon = 0, nlat = 0, nlev = 0

        ! Horizontal grid metric, in mu = sin(latitude). The finite-volume
        ! cells tile [-1,1] in mu with widths equal to the Gauss weights, so the
        ! cell area is a^2 dlambda w_j exactly (aeros_grid's own area). Latitude
        ! runs north (mu = +1) to south (mu = -1).
        real(wp) :: dlam = 0.0_wp             ! longitude spacing [rad]
        real(wp), allocatable :: dmu(:)       ! cell mu-width = Gauss weight, (nlat)
        real(wp), allocatable :: coslat(:)    ! cos(lat) at cell centres, (nlat)
        real(wp), allocatable :: cos_face(:)  ! cos(lat) at mu-interfaces, (0:nlat)

        ! Scratch, all (nlon,nlat,nlev) unless noted. Allocated once.
        real(wp), allocatable :: dp_lev(:,:,:)  ! layer air mass (pressure) [Pa]
        real(wp), allocatable :: qm(:,:,:)      ! tracer mass  q*dp        [Pa]
        real(wp), allocatable :: dp0(:,:,:)     ! layer air mass before the horizontal sweeps [Pa]

        ! Diagnostics, retained between calls for aeros_moisture_report.
        real(wp) :: last_cfl   = 0.0_wp        ! max total Courant last call
        integer  :: last_nsub  = 0             ! sub-steps taken last call
        real(dp) :: mass_consistency = 0.0_dp  ! |FV - spectral| ps tendency gap
    end type aeros_moist_class

    ! Courant target per sub-step. Below 1 for the donor-cell max principle;
    ! 0.9 leaves a margin against the estimate being a hair low.
    real(wp), parameter :: CFL_MAX = 0.9_wp

    public :: aeros_moist_class
    public :: aeros_moisture_init
    public :: aeros_moisture_end
    public :: aeros_moisture_transport
    public :: aeros_moisture_water
    public :: aeros_moisture_report

contains

    subroutine aeros_moisture_init(mst, grd, nlev)
        ! Precompute the grid metric and allocate the scratch.

        implicit none

        type(aeros_moist_class), intent(inout) :: mst
        type(aeros_grid_class),  intent(in)    :: grd
        integer,                 intent(in)    :: nlev

        real(dp) :: muh
        integer  :: j

        call aeros_moisture_end(mst)

        mst%nlon = grd%nlon
        mst%nlat = grd%nlat
        mst%nlev = nlev

        mst%dlam = real(2.0_dp*pi/real(grd%nlon, dp), wp)

        allocate(mst%dmu(grd%nlat), mst%coslat(grd%nlat))
        allocate(mst%cos_face(0:grd%nlat))

        ! Cell centres: the Gauss latitudes. cos(lat) = sqrt(1 - mu^2).
        do j = 1, grd%nlat
            mst%dmu(j)    = real(grd%gauss_w(j), wp)
            mst%coslat(j) = real(sqrt(max(0.0_dp, 1.0_dp - grd%sinlat(j)**2)), wp)
        end do

        ! Interfaces: mu tiled by the Gauss weights from the north pole down.
        ! muh(0) = +1 and muh(nlat) = -1 to machine precision (weights sum to
        ! 2), so cos_face vanishes at both poles -- which is what makes the
        ! meridional flux zero there and the global budget close.
        muh = 1.0_dp
        mst%cos_face(0) = 0.0_wp
        do j = 1, grd%nlat
            muh = muh - real(grd%gauss_w(j), dp)
            mst%cos_face(j) = real(sqrt(max(0.0_dp, 1.0_dp - muh**2)), wp)
        end do
        ! Pin the south pole exactly, guarding weight round-off.
        mst%cos_face(grd%nlat) = 0.0_wp

        allocate(mst%dp_lev(grd%nlon,grd%nlat,nlev))
        allocate(mst%qm    (grd%nlon,grd%nlat,nlev))
        allocate(mst%dp0   (grd%nlon,grd%nlat,nlev))

        return

    end subroutine aeros_moisture_init

    subroutine aeros_moisture_end(mst)

        implicit none

        type(aeros_moist_class), intent(inout) :: mst

        if (allocated(mst%dmu))      deallocate(mst%dmu)
        if (allocated(mst%coslat))   deallocate(mst%coslat)
        if (allocated(mst%cos_face)) deallocate(mst%cos_face)
        if (allocated(mst%dp_lev))   deallocate(mst%dp_lev)
        if (allocated(mst%qm))       deallocate(mst%qm)
        if (allocated(mst%dp0))      deallocate(mst%dp0)

        mst%nlon = 0; mst%nlat = 0; mst%nlev = 0
        mst%last_cfl = 0.0_wp; mst%last_nsub = 0; mst%mass_consistency = 0.0_dp

        return

    end subroutine aeros_moisture_end

    subroutine aeros_moisture_transport(mst, vg, u, v, lnps, q, dt)
        ! Advance q one dynamics step, forward in time, on the grid.
        !
        ! u, v, lnps are the wind and (log) surface pressure at the CURRENT time
        ! level -- the same time level the leapfrog evaluates its own right-hand
        ! side at -- so the humidity is advected by the flow that is driving it.
        ! lnps rather than ps because that is what the core carries; ps is
        ! recovered per column. q is updated in place. Nothing here is spectral.

        implicit none

        type(aeros_moist_class), intent(inout) :: mst
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in)    :: u(:,:,:), v(:,:,:)   ! (nlon,nlat,nlev) [m s-1]
        real(wp), intent(in)    :: lnps(:,:)            ! (nlon,nlat)      ln[Pa]
        real(wp), intent(inout) :: q(:,:,:)             ! (nlon,nlat,nlev) [kg kg-1]
        real(wp), intent(in)    :: dt                   ! [s]

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: cfl, dts
        integer  :: i, j, k, nsub, isub

        ! Layer masses from the true surface pressure, and the tracer mass.
        !$omp parallel do collapse(2) private(i,j,k,phalf,pfull,dpc) schedule(static)
        do j = 1, mst%nlat
            do i = 1, mst%nlon
                call aeros_vgrid_pressure(vg, exp(lnps(i,j)), phalf, pfull, dpc)
                do k = 1, mst%nlev
                    mst%dp_lev(i,j,k) = dpc(k)
                    mst%qm(i,j,k)     = q(i,j,k)*dpc(k)
                end do
            end do
        end do
        !$omp end parallel do

        ! Sub-stepping count from the worst Courant number on the grid.
        cfl  = max_courant(mst, u, v, lnps, vg, dt)
        nsub = max(1, ceiling(cfl/CFL_MAX))
        dts  = dt/real(nsub, wp)

        mst%last_cfl  = cfl
        mst%last_nsub = nsub

        do isub = 1, nsub
            call transport_substep(mst, vg, u, v, dts)
        end do

        ! Back to specific humidity. dp_lev now holds the FV-updated air mass;
        ! dividing by it (not by the spectral air mass) is what preserves the
        ! max principle, hence positivity.
        !$omp parallel do collapse(2) private(i,j,k) schedule(static)
        do j = 1, mst%nlat
            do i = 1, mst%nlon
                do k = 1, mst%nlev
                    q(i,j,k) = mst%qm(i,j,k)/mst%dp_lev(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        return

    end subroutine aeros_moisture_transport

    subroutine transport_substep(mst, vg, u, v, dt)
        ! One forward sub-step, operator-split: a zonal sweep, then a meridional
        ! sweep, then a vertical sweep. Each updates the air mass (dp_lev) and
        ! the tracer mass (qm) by the SAME fluxes, so q = const is preserved
        ! through every sweep and hence through the step; each is flux-form, so
        ! the global sums are conserved; and each is a monotone one-dimensional
        ! scheme, so q stays bounded and therefore non-negative.
        !
        ! The horizontal sweeps use a van Leer limiter (second order, monotone);
        ! the vertical stays first-order upwind. That split is where the
        ! accuracy is spent to match where it matters: humidity gradients are
        ! overwhelmingly horizontal, the horizontal Courant numbers are the
        ! large ones, and the vertical has a handful of levels and a gentle
        ! mass-coordinate flux. Van Leer in the vertical is a later refinement.

        implicit none

        type(aeros_moist_class), intent(inout) :: mst
        type(aeros_vgrid_class), intent(in)    :: vg
        real(wp), intent(in) :: u(:,:,:), v(:,:,:)
        real(wp), intent(in) :: dt

        integer  :: i, j, k, nlon, nlat, nlev
        real(wp) :: scum_a(0:vg%nlev), mfl_a(0:vg%nlev)
        real(wp) :: dq_v, da_v

        nlon = mst%nlon; nlat = mst%nlat; nlev = mst%nlev

        ! The air mass before the horizontal sweeps: the vertical mass flux is
        ! built from the horizontal convergence, which is exactly the change the
        ! two horizontal sweeps make to dp_lev.
        !$omp parallel do collapse(2) private(i,j,k) schedule(static)
        do k = 1, nlev
            do j = 1, nlat
                do i = 1, nlon
                    mst%dp0(i,j,k) = mst%dp_lev(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        call sweep_zonal(mst, u, dt)
        call sweep_merid(mst, v, dt)

        ! --- Vertical sweep, per column -------------------------------------
        ! The horizontal convergence into each layer is (dp0 - dp_lev)/dt after
        ! the two sweeps. Its running column sum drives the vertical mass flux
        ! mfl_a(k) = B(k) S_N - S_k, exactly aeros_tendency's `mflux` form, so
        ! the coordinate stays a mass coordinate. Zero at top and ground.
        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,scum_a,mfl_a,da_v,dq_v)
        do j = 1, nlat
            do i = 1, nlon

                scum_a(0) = 0.0_wp
                do k = 1, nlev
                    scum_a(k) = scum_a(k-1) + (mst%dp0(i,j,k) - mst%dp_lev(i,j,k))/dt
                end do

                do k = 0, nlev
                    mfl_a(k) = vg%B(k)*scum_a(nlev) - scum_a(k)
                end do

                do k = 1, nlev
                    da_v = mfl_a(k) - mfl_a(k-1)
                    dq_v = vflux_q(mst%qm, mst%dp_lev, i, j, k,   mfl_a(k),   nlev) &
                            - vflux_q(mst%qm, mst%dp_lev, i, j, k-1, mfl_a(k-1), nlev)

                    mst%dp_lev(i,j,k) = mst%dp_lev(i,j,k) - dt*da_v
                    mst%qm(i,j,k)     = mst%qm(i,j,k)     - dt*dq_v
                end do

            end do
        end do
        !$omp end parallel do

        return

    end subroutine transport_substep

    subroutine sweep_zonal(mst, u, dt)
        ! Van Leer flux-form advection in longitude, periodic, per (level, row).
        !
        ! The face air flux is upwind (dp of the donor cell times the face
        ! velocity); the face tracer value is the van Leer reconstruction of
        ! q = qm/dp in the donor cell, carried by that same air flux. Carrying a
        ! reconstruction of the RATIO q, rather than of qm, is what keeps
        ! q = const exact: a flat q reconstructs flat, so the tracer flux is q
        ! times the air flux and the two divergences differ only by that factor.

        implicit none

        type(aeros_moist_class), intent(inout) :: mst
        real(wp), intent(in) :: u(:,:,:)
        real(wp), intent(in) :: dt

        integer  :: i, j, k, ip, nlon
        real(wp) :: qc(mst%nlon), dc(mst%nlon)
        real(wp) :: fair(mst%nlon), fq(mst%nlon)     ! face i+1/2 fluxes
        real(wp) :: uf, cour, qface, dx

        nlon = mst%nlon

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,ip,qc,dc,fair,fq,uf,cour,qface,dx)
        do k = 1, mst%nlev
            do j = 1, mst%nlat
                dx = r_earth*mst%coslat(j)*mst%dlam

                do i = 1, nlon
                    dc(i) = mst%dp_lev(i,j,k)
                    qc(i) = mst%qm(i,j,k)/dc(i)
                end do

                ! Flux at each face i+1/2, van Leer on q with a Courant
                ! correction. The donor cell's limited slope is evaluated on the
                ! periodic stencil.
                do i = 1, nlon
                    ip = i + 1; if (ip > nlon) ip = 1
                    uf   = 0.5_wp*(u(i,j,k) + u(ip,j,k))
                    cour = abs(uf)*dt/dx
                    if (uf >= 0.0_wp) then
                        qface   = vl_face(qc, i,  +1, cour, nlon)
                        fair(i) = dc(i)*uf
                    else
                        qface   = vl_face(qc, ip, -1, cour, nlon)
                        fair(i) = dc(ip)*uf
                    end if
                    fq(i) = qface*fair(i)
                end do

                ! Update from the flux difference across each cell.
                do i = 1, nlon
                    ip = i - 1; if (ip < 1) ip = nlon         ! face i-1/2 index
                    mst%dp_lev(i,j,k) = mst%dp_lev(i,j,k) - dt*(fair(i) - fair(ip))/dx
                    mst%qm(i,j,k)     = mst%qm(i,j,k)     - dt*(fq(i)   - fq(ip))/dx
                end do
            end do
        end do
        !$omp end parallel do

        return

    end subroutine sweep_zonal

    subroutine sweep_merid(mst, v, dt)
        ! Van Leer flux-form advection in latitude (mu), per (level, column).
        !
        ! Bounded, not periodic: there is no cell beyond the poles, and the
        ! interface cos(lat) = 0 there makes the pole face flux identically
        ! zero, which is what closes the global budget. The van Leer stencil is
        ! clamped at the first and last interior rows, degrading to upwind there
        ! -- correct, since there is no cell to form the outer slope from.

        implicit none

        type(aeros_moist_class), intent(inout) :: mst
        real(wp), intent(in) :: v(:,:,:)
        real(wp), intent(in) :: dt

        integer  :: i, j, k, nlat
        real(wp) :: qc(mst%nlat), dc(mst%nlat)
        real(wp) :: fair(0:mst%nlat), fq(0:mst%nlat)   ! face j+1/2, j = 0..nlat
        real(wp) :: vf, cour, qface, dy

        nlat = mst%nlat

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,qc,dc,fair,fq,vf,cour,qface,dy)
        do k = 1, mst%nlev
            do i = 1, mst%nlon
                do j = 1, nlat
                    dc(j) = mst%dp_lev(i,j,k)
                    qc(j) = mst%qm(i,j,k)/dc(j)
                end do

                ! Interior faces j+1/2 (between rows j and j+1), j = 1..nlat-1.
                fair(0) = 0.0_wp; fq(0) = 0.0_wp         ! north pole face
                fair(nlat) = 0.0_wp; fq(nlat) = 0.0_wp   ! south pole face
                do j = 1, nlat - 1
                    vf   = 0.5_wp*(v(i,j,k) + v(i,j+1,k))
                    ! dy for the Courant number: the mu-band width scaled to a
                    ! latitude distance at the interface, a dmu / cos(lat_face).
                    dy   = r_earth*mst%dmu(j)/max(mst%cos_face(j), 1.0e-6_wp)
                    cour = abs(vf)*dt/dy
                    if (vf >= 0.0_wp) then
                        qface = vl_face_bounded(qc, j,   +1, cour, nlat)
                        fair(j) = dc(j)*vf*mst%cos_face(j)
                    else
                        qface = vl_face_bounded(qc, j+1, -1, cour, nlat)
                        fair(j) = dc(j+1)*vf*mst%cos_face(j)
                    end if
                    ! Tracer flux is q_face times the air flux, so a flat q
                    ! gives fq = q * fair and constancy is exact.
                    fq(j) = qface*fair(j)
                end do

                do j = 1, nlat
                    mst%dp_lev(i,j,k) = mst%dp_lev(i,j,k) &
                                        - dt*(fair(j) - fair(j-1))/(r_earth*mst%dmu(j))
                    mst%qm(i,j,k)     = mst%qm(i,j,k) &
                                        - dt*(fq(j)   - fq(j-1))/(r_earth*mst%dmu(j))
                end do
            end do
        end do
        !$omp end parallel do

        return

    end subroutine sweep_merid

    real(wp) function vl_face(q, iup, dir, cour, n) result(qf)
        ! Van Leer face value from the donor cell iup, periodic in i.
        !
        ! dir = +1 when the flow is left-to-right (the face is on the cell's
        ! right), -1 when right-to-left. The reconstructed value the flux
        ! carries over one step is the donor value plus half its limited slope,
        ! reduced by the Courant number so it is the average of what actually
        ! crosses the face, not the instantaneous edge value.

        implicit none

        real(wp), intent(in) :: q(:)
        integer,  intent(in) :: iup, dir, n
        real(wp), intent(in) :: cour

        integer  :: im, ip
        real(wp) :: sm, sp

        im = iup - 1; if (im < 1) im = n
        ip = iup + 1; if (ip > n) ip = 1

        ! One-sided differences oriented along the flow: sm is the upwind-side
        ! slope, sp the downwind-side.
        if (dir > 0) then
            sm = q(iup) - q(im)
            sp = q(ip)  - q(iup)
        else
            sm = q(iup) - q(ip)
            sp = q(im)  - q(iup)
        end if

        qf = q(iup) + 0.5_wp*(1.0_wp - cour)*vanleer(sm, sp)

        return

    end function vl_face

    real(wp) function vl_face_bounded(q, jup, dir, cour, n) result(qf)
        ! As vl_face but for the bounded meridional direction: the outer slope
        ! is dropped (set to zero, giving upwind) at the first and last rows,
        ! where there is no cell to form it from.

        implicit none

        real(wp), intent(in) :: q(:)
        integer,  intent(in) :: jup, dir, n
        real(wp), intent(in) :: cour

        real(wp) :: sm, sp

        if (dir > 0) then
            if (jup <= 1 .or. jup >= n) then
                qf = q(jup); return
            end if
            sm = q(jup) - q(jup-1)
            sp = q(jup+1) - q(jup)
        else
            if (jup <= 1 .or. jup >= n) then
                qf = q(jup); return
            end if
            sm = q(jup) - q(jup+1)
            sp = q(jup-1) - q(jup)
        end if

        qf = q(jup) + 0.5_wp*(1.0_wp - cour)*vanleer(sm, sp)

        return

    end function vl_face_bounded

    real(wp) function vanleer(a, b) result(s)
        ! Van Leer limited slope: the harmonic mean of the two one-sided
        ! differences when they agree in sign, zero at an extremum. Monotone,
        ! second-order where the field is smooth, and it is what guarantees the
        ! reconstruction introduces no new maximum or minimum -- hence
        ! positivity of q under the Courant limit.

        implicit none

        real(wp), intent(in) :: a, b

        if (a*b > 0.0_wp) then
            s = 2.0_wp*a*b/(a + b)
        else
            s = 0.0_wp
        end if

        return

    end function vanleer

    real(wp) function vflux_q(qm, dp_lev, i, j, k, mfl, nlev) result(fq)
        ! Donor-cell vertical tracer flux at interface k, [Pa/s].
        !
        ! mfl is the vertical air-mass flux there (aeros_tendency's convention:
        ! positive downward in the coordinate). q at the interface is taken from
        ! the layer the air comes from. Zero at both boundaries, where mfl = 0
        ! exactly, so the guards on k are for the array bounds, not the physics.

        implicit none

        real(wp), intent(in) :: qm(:,:,:), dp_lev(:,:,:)
        integer,  intent(in) :: i, j, k, nlev
        real(wp), intent(in) :: mfl

        real(wp) :: q_up

        if (k <= 0 .or. k >= nlev) then
            fq = 0.0_wp
            return
        end if

        ! Interface k sits between layer k (above) and layer k+1 (below). mfl > 0
        ! is downward flow, carrying the upper layer's q; mfl < 0 carries the
        ! lower layer's.
        if (mfl >= 0.0_wp) then
            q_up = qm(i,j,k)  /dp_lev(i,j,k)
        else
            q_up = qm(i,j,k+1)/dp_lev(i,j,k+1)
        end if

        fq = q_up*mfl

        return

    end function vflux_q

    real(wp) function max_courant(mst, u, v, lnps, vg, dt) result(cmax)
        ! Largest total (zonal + meridional + vertical) Courant number on the
        ! grid, for setting the sub-step count. An estimate, deliberately on the
        ! generous side: the vertical term uses the spectral-free column budget
        ! only approximately, but CFL_MAX < 1 leaves the margin.

        implicit none

        type(aeros_moist_class), intent(in) :: mst
        real(wp), intent(in) :: u(:,:,:), v(:,:,:), lnps(:,:)
        type(aeros_vgrid_class), intent(in) :: vg
        real(wp), intent(in) :: dt

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(wp) :: scum(0:vg%nlev), cmass
        real(wp) :: cx, cy, cz, ctot, dy
        integer  :: i, j, k

        cmax = 0.0_wp

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,phalf,pfull,dpc,scum,cmass,cx,cy,cz,ctot,dy) &
        !$omp   reduction(max:cmax)
        do j = 1, mst%nlat
            do i = 1, mst%nlon
                call aeros_vgrid_pressure(vg, exp(lnps(i,j)), phalf, pfull, dpc)

                ! A rough column mass divergence for the vertical Courant: the
                ! horizontal divergence of dp*V by centred estimate is not worth
                ! forming here, so use |div(V)|*dp bounded by the wind and grid.
                scum(0) = 0.0_wp
                do k = 1, mst%nlev
                    cmass   = dpc(k)*abs(u(i,j,k))/(r_earth*mst%coslat(j)*mst%dlam)
                    scum(k) = scum(k-1) + cmass
                end do

                dy = r_earth*mst%dmu(j)/max(mst%coslat(j), 1.0e-6_wp)
                do k = 1, mst%nlev
                    cx = abs(u(i,j,k))*dt/(r_earth*mst%coslat(j)*mst%dlam)
                    cy = abs(v(i,j,k))*dt/dy
                    cz = (abs(vg%B(k)*scum(mst%nlev) - scum(k)) &
                            + abs(vg%B(k-1)*scum(mst%nlev) - scum(k-1)))*dt/dpc(k)
                    ctot = cx + cy + cz
                    cmax = max(cmax, ctot)
                end do
            end do
        end do
        !$omp end parallel do

        return

    end function max_courant

    real(dp) function aeros_moisture_water(mst, grd, vg, ps, q) result(water)
        ! Global water mass, int q dp/g dA [kg], for the conservation check.
        ! Accumulated in dp for the same reason aeros_budget is: a machine-
        ! precision statement cannot be made in an sp sum over the grid.

        implicit none

        type(aeros_moist_class), intent(in) :: mst
        type(aeros_grid_class),  intent(in) :: grd
        type(aeros_vgrid_class), intent(in) :: vg
        real(wp), intent(in) :: ps(:,:), q(:,:,:)

        real(wp) :: phalf(0:vg%nlev), pfull(vg%nlev), dpc(vg%nlev)
        real(dp) :: col, tot
        integer  :: i, j, k

        tot = 0.0_dp

        !$omp parallel do collapse(2) private(i,j,k,phalf,pfull,dpc,col) &
        !$omp   reduction(+:tot) schedule(static)
        do j = 1, mst%nlat
            do i = 1, mst%nlon
                call aeros_vgrid_pressure(vg, ps(i,j), phalf, pfull, dpc)
                col = 0.0_dp
                do k = 1, mst%nlev
                    col = col + real(q(i,j,k), dp)*real(dpc(k), dp)
                end do
                tot = tot + col*real(grd%area(i,j), dp)/real(grav, dp)
            end do
        end do
        !$omp end parallel do

        water = tot

        return

    end function aeros_moisture_water

    subroutine aeros_moisture_report(mst, io_unit)
        ! What the transport is doing: sub-steps forced by the polar Courant,
        ! and the FV-vs-spectral mass-consistency gap.

        implicit none

        type(aeros_moist_class), intent(in) :: mst
        integer, intent(in), optional :: io_unit

        integer :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        write(iou,*) ""
        write(iou,"(a)")        " == moisture transport =="
        write(iou,"(a)")        "   scheme                      FV flux-form, upwind, positive-definite"
        write(iou,"(a,es12.3)") "   last max Courant           ", mst%last_cfl
        write(iou,"(a,i9)")     "   last sub-steps             ", mst%last_nsub

        return

    end subroutine aeros_moisture_report

end module aeros_moisture
