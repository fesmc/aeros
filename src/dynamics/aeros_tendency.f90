module aeros_tendency
    ! Right-hand sides of the primitive equations, by the transform method.
    !
    ! Given the spectral state (zeta, D, T, ln p_s) this returns the spectral
    ! tendencies of the same four. Nothing here integrates in time -- the
    ! semi-implicit leapfrog, the Asselin filter and the horizontal diffusion
    ! are M1.4, and the Held-Suarez forcing is M1.5.
    !
    ! === Discretization ======================================================
    !
    ! Simmons & Burridge (1981), the same vertical scheme as the geopotential
    ! in aeros_vertical, and for the same reason: the discrete hydrostatic
    ! operator, the pressure-gradient term and the omega term must SHARE their
    ! coefficients, or the discretization does not conserve energy. Two
    ! coefficients carry that pairing, and both are built once per column:
    !
    !   alpha_k   the Simmons-Burridge layer coefficient (aeros_vertical)
    !
    !   c_k       the layer-mean of B p_s / p, which is the factor multiplying
    !             grad(ln p_s) in BOTH the pressure-gradient force and the
    !             first term of omega/p:
    !
    !                 c_k = (p_s/dp_k) [ B_(k-1/2) dlnp_k + dB_k alpha_k ]
    !
    !             It is 1 in a pure sigma coordinate -- verified to machine
    !             precision for every layer when p_top > 0, and for every layer
    !             but the top when p_top = 0, where Simmons & Burridge's
    !             alpha_1 = ln 2 convention makes c_1 = ln 2. That top-layer
    !             value is deliberate: it is more important that the three
    !             terms share a coefficient than that any one of them is
    !             individually exact in a layer whose upper interface is at
    !             zero pressure.
    !
    ! Column quantities, top-down, with G_k = v_k . grad(ln p_s):
    !
    !   C_k   = dp_k D_k + dB_k p_s G_k        = div(v_k dp_k), exactly
    !   S_k   = sum_(j<=k) C_j                   cumulative mass divergence
    !   M_(k+1/2) = B_(k+1/2) S_N - S_k          vertical mass flux; 0 at both ends
    !
    !   d(ln p_s)/dt = -S_N / p_s
    !   (omega/p)_k  = c_k G_k - (1/dp_k)[ S_(k-1) dlnp_k + alpha_k C_k ]
    !
    ! Momentum, in vector-invariant form, as an EAST/NORTH grid field A:
    !
    !   A = -(zeta + f) k x v  -  vertical advection  -  R_d T c grad(ln p_s)
    !
    !   d(zeta)/dt = curl(A)
    !   d(D)/dt    = div(A) - laplacian(Phi + K),      K = (u^2 + v^2)/2
    !
    ! The gradient of (Phi + K) is deliberately NOT formed on the grid and
    ! re-analyzed. Taking its Laplacian spectrally instead is cheaper (one
    ! scalar transform rather than a vector pair) and, more importantly, it
    ! avoids computing a gradient and then a divergence of the same field on a
    ! quadratic grid -- which aliases. The curl of A is unaffected either way,
    ! since curl(grad) = 0.
    !
    ! Thermodynamics, advective form:
    !
    !   dT/dt = -v . grad(T) - vertical advection + kappa T (omega/p)
    !
    ! === Transform count =====================================================
    !
    ! Per level, synthesize (u,v) [2 scalar-equivalents], zeta [1], D [1],
    ! T [1], grad(T) [2]; analyze A [2], Phi+K [1], dT/dt [1]. That is ELEVEN
    ! scalar-equivalent transforms per level per call, against the EIGHT
    ! assumed by drivers/bench_m0a.F90.
    !
    ! So the M0a cost estimate is ~37% optimistic, and every core-s/yr figure
    ! in docs/m0a_results.md should be read as multiplied by 11/8. The headline
    ! conclusion survives comfortably -- T42L20 goes from ~14% of the coupled
    ! budget to ~19%, still nothing like the 82% design.md section 3.6
    ! assumed -- but the number itself has moved and the doc says so.
    !
    ! zeta and D are synthesized rather than reconstructed because the column
    ! loop genuinely needs both on the grid: zeta for the absolute-vorticity
    ! flux, D for the layer mass budget. SHTns' combined SHqst_to_spat could
    ! deliver zeta alongside (u,v) in one call and is documented as
    ! "significantly faster" than separate transforms -- worth trying once
    ! there is a working core to measure it against.
    !
    ! Plus four single-level transforms for ln p_s, which are O(1/nlev).
    !
    ! === Robert form =========================================================
    !
    ! Still not enabled (see aeros_vordiv). The formulation above never divides
    ! by cos(latitude): the (u,v) it works with come from SHTns' vector
    ! transforms, which handle the metric internally, and the only place a
    ! 1/cos could appear -- the horizontal gradients -- is likewise handled by
    ! SHsph_to_spat. So the problem Robert form solves does not arise here, and
    ! turning it on would change every vector transform in the model to buy
    ! nothing. Revisit only if a polar-noise problem actually shows up.

    use aeros_defs,     only : dp, wp, wp_sh, io_unit_err, MV, R_d, kappa, &
                                aeros_grid_class, aeros_spec_class
    use aeros_spectral, only : aeros_sht_class, aeros_sht_pool_class, &
                                aeros_sht_pool_get, aeros_sht_analysis, &
                                aeros_sht_synthesis, aeros_sht_laplacian
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure, &
                                aeros_vgrid_alpha, aeros_hydrostatic
    use aeros_vordiv,   only : aeros_uv_from_vordiv, aeros_vordiv_from_uv, &
                                aeros_gradient

    implicit none

    private

    type aeros_tend_class
        ! Spectral tendencies, same shape as the spectral half of the state.
        integer :: nlm = 0, nlev = 0
        complex(wp_sh), allocatable :: vor(:,:)   ! [s-2]
        complex(wp_sh), allocatable :: div(:,:)   ! [s-2]
        complex(wp_sh), allocatable :: temp(:,:)  ! [K s-1]
        complex(wp_sh), allocatable :: lnps(:)    ! [s-1]
    end type aeros_tend_class

    type aeros_work_class
        ! Grid-space scratch, allocated once and reused every timestep.
        !
        ! Shared across threads: each thread writes its own level slice in the
        ! level loops, and its own column in the column loop, so there is no
        ! overlap. Allocated rather than automatic because these are megabytes
        ! -- at T42L20 each 3-D field is 13 MB -- which does not belong on an
        ! OpenMP thread stack.
        !
        ! Layout is (nlon, nlat, nlev): a LEVEL is contiguous, because that is
        ! what the transforms require. A column is therefore strided, which the
        ! column loop pays for by copying into local nlev buffers. The opposite
        ! layout would suit the column physics and penalise every transform;
        ! this is the cheaper side of that trade, but it is a trade.

        integer :: nlon = 0, nlat = 0, nlev = 0

        real(wp), allocatable :: u(:,:,:), v(:,:,:)        ! winds [m s-1]
        real(wp), allocatable :: vor_g(:,:,:)              ! relative vorticity [s-1]
        real(wp), allocatable :: div_g(:,:,:)              ! divergence [s-1]
        real(wp), allocatable :: t_g(:,:,:)                ! temperature [K]
        real(wp), allocatable :: dtdx(:,:,:), dtdy(:,:,:)  ! grad T [K m-1]

        real(wp), allocatable :: ae(:,:,:), an(:,:,:)      ! momentum RHS, east/north
        real(wp), allocatable :: ek(:,:,:)                 ! Phi + K [m2 s-2]
        real(wp), allocatable :: dtdt(:,:,:)               ! dT/dt [K s-1]

        real(wp), allocatable :: lnps_g(:,:)               ! ln(p_s)
        real(wp), allocatable :: dlnpsdx(:,:), dlnpsdy(:,:)! grad ln p_s [m-1]
        real(wp), allocatable :: dlnpsdt(:,:)              ! [s-1]

        ! Surface geopotential [m2 s-2]. A prescribed boundary field, zero
        ! until topography arrives at M2 (and zero for Held-Suarez, which is
        ! defined over a flat surface).
        real(wp), allocatable :: phis(:,:)
    end type aeros_work_class

    public :: aeros_tend_class
    public :: aeros_work_class
    public :: aeros_tend_alloc
    public :: aeros_tend_end
    public :: aeros_work_alloc
    public :: aeros_work_end
    public :: aeros_tendency_calc

contains

    subroutine aeros_tend_alloc(tnd, nlm, nlev)

        implicit none

        type(aeros_tend_class), intent(inout) :: tnd
        integer, intent(in) :: nlm, nlev

        call aeros_tend_end(tnd)

        tnd%nlm  = nlm
        tnd%nlev = nlev

        allocate(tnd%vor(nlm,nlev), tnd%div(nlm,nlev), tnd%temp(nlm,nlev))
        allocate(tnd%lnps(nlm))

        tnd%vor  = (0.0_wp_sh, 0.0_wp_sh)
        tnd%div  = (0.0_wp_sh, 0.0_wp_sh)
        tnd%temp = (0.0_wp_sh, 0.0_wp_sh)
        tnd%lnps = (0.0_wp_sh, 0.0_wp_sh)

        return

    end subroutine aeros_tend_alloc

    subroutine aeros_tend_end(tnd)

        implicit none

        type(aeros_tend_class), intent(inout) :: tnd

        tnd%nlm = 0; tnd%nlev = 0
        if (allocated(tnd%vor))  deallocate(tnd%vor)
        if (allocated(tnd%div))  deallocate(tnd%div)
        if (allocated(tnd%temp)) deallocate(tnd%temp)
        if (allocated(tnd%lnps)) deallocate(tnd%lnps)

        return

    end subroutine aeros_tend_end

    subroutine aeros_work_alloc(wrk, nlon, nlat, nlev)

        implicit none

        type(aeros_work_class), intent(inout) :: wrk
        integer, intent(in) :: nlon, nlat, nlev

        call aeros_work_end(wrk)

        wrk%nlon = nlon; wrk%nlat = nlat; wrk%nlev = nlev

        allocate(wrk%u(nlon,nlat,nlev), wrk%v(nlon,nlat,nlev))
        allocate(wrk%vor_g(nlon,nlat,nlev), wrk%div_g(nlon,nlat,nlev))
        allocate(wrk%t_g(nlon,nlat,nlev))
        allocate(wrk%dtdx(nlon,nlat,nlev), wrk%dtdy(nlon,nlat,nlev))
        allocate(wrk%ae(nlon,nlat,nlev), wrk%an(nlon,nlat,nlev))
        allocate(wrk%ek(nlon,nlat,nlev), wrk%dtdt(nlon,nlat,nlev))
        allocate(wrk%lnps_g(nlon,nlat))
        allocate(wrk%dlnpsdx(nlon,nlat), wrk%dlnpsdy(nlon,nlat))
        allocate(wrk%dlnpsdt(nlon,nlat))
        allocate(wrk%phis(nlon,nlat))

        wrk%phis = 0.0_wp

        return

    end subroutine aeros_work_alloc

    subroutine aeros_work_end(wrk)

        implicit none

        type(aeros_work_class), intent(inout) :: wrk

        wrk%nlon = 0; wrk%nlat = 0; wrk%nlev = 0

        if (allocated(wrk%u))       deallocate(wrk%u)
        if (allocated(wrk%v))       deallocate(wrk%v)
        if (allocated(wrk%vor_g))   deallocate(wrk%vor_g)
        if (allocated(wrk%div_g))   deallocate(wrk%div_g)
        if (allocated(wrk%t_g))     deallocate(wrk%t_g)
        if (allocated(wrk%dtdx))    deallocate(wrk%dtdx)
        if (allocated(wrk%dtdy))    deallocate(wrk%dtdy)
        if (allocated(wrk%ae))      deallocate(wrk%ae)
        if (allocated(wrk%an))      deallocate(wrk%an)
        if (allocated(wrk%ek))      deallocate(wrk%ek)
        if (allocated(wrk%dtdt))    deallocate(wrk%dtdt)
        if (allocated(wrk%lnps_g))  deallocate(wrk%lnps_g)
        if (allocated(wrk%dlnpsdx)) deallocate(wrk%dlnpsdx)
        if (allocated(wrk%dlnpsdy)) deallocate(wrk%dlnpsdy)
        if (allocated(wrk%dlnpsdt)) deallocate(wrk%dlnpsdt)
        if (allocated(wrk%phis))    deallocate(wrk%phis)

        return

    end subroutine aeros_work_end

    subroutine aeros_tendency_calc(pool, vg, grd, now, wrk, tnd)
        ! One evaluation of the primitive-equation right-hand sides.

        implicit none

        type(aeros_sht_pool_class), intent(in), target :: pool
        type(aeros_vgrid_class),    intent(in)    :: vg
        type(aeros_grid_class),     intent(in)    :: grd
        type(aeros_spec_class),     intent(in)    :: now
        type(aeros_work_class),     intent(inout) :: wrk
        type(aeros_tend_class),     intent(inout) :: tnd

        type(aeros_sht_class), pointer :: s
        complex(wp_sh) :: ek_lm(pool%sht(1)%nlm)
        integer :: k, nlev

        nlev = vg%nlev

        ! === 1. Single-level: surface pressure and its gradient ==============
        s => pool%sht(1)
        call aeros_sht_synthesis(s, now%lnps, wrk%lnps_g)
        call aeros_gradient(s, now%lnps, wrk%dlnpsdx, wrk%dlnpsdy)

        ! === 2. Spectral -> grid, one thread per level =======================
        !$omp parallel do num_threads(pool%nthreads) private(k,s) schedule(static)
        do k = 1, nlev
            s => aeros_sht_pool_get(pool)
            call aeros_uv_from_vordiv(s, now%vor(:,k), now%div(:,k), &
                                        wrk%u(:,:,k), wrk%v(:,:,k))
            call aeros_sht_synthesis(s, now%vor(:,k), wrk%vor_g(:,:,k))
            call aeros_sht_synthesis(s, now%div(:,k), wrk%div_g(:,:,k))
            call aeros_sht_synthesis(s, now%temp(:,k), wrk%t_g(:,:,k))
            call aeros_gradient(s, now%temp(:,k), wrk%dtdx(:,:,k), wrk%dtdy(:,:,k))
        end do
        !$omp end parallel do

        ! === 3. Column loop: everything vertical ============================
        call column_terms(vg, grd, wrk)

        ! === 4. Grid -> spectral, one thread per level ======================
        !$omp parallel do num_threads(pool%nthreads) private(k,s,ek_lm) schedule(static)
        do k = 1, nlev
            s => aeros_sht_pool_get(pool)

            ! curl and divergence of the momentum RHS, in one vector analysis.
            call aeros_vordiv_from_uv(s, wrk%ae(:,:,k), wrk%an(:,:,k), &
                                        tnd%vor(:,k), tnd%div(:,k))

            ! -laplacian(Phi + K). NB analysis destroys wrk%ek, which is
            ! rebuilt from scratch on every call.
            call aeros_sht_analysis(s, wrk%ek(:,:,k), ek_lm)
            call aeros_sht_laplacian(s, ek_lm)
            tnd%div(:,k) = tnd%div(:,k) - ek_lm

            call aeros_sht_analysis(s, wrk%dtdt(:,:,k), tnd%temp(:,k))
        end do
        !$omp end parallel do

        ! === 5. Single-level: surface pressure tendency ======================
        s => pool%sht(1)
        call aeros_sht_analysis(s, wrk%dlnpsdt, tnd%lnps)

        return

    end subroutine aeros_tendency_calc

    subroutine column_terms(vg, grd, wrk)
        ! The vertical half of the right-hand sides, column by column.
        !
        ! Threaded over columns, which is the embarrassingly-parallel part of
        ! the model (docs/design.md section 4.3) and the part that grows when
        ! physics arrives. Local nlev buffers rather than strided access into
        ! the (nlon,nlat,nlev) work arrays -- see aeros_work_class.

        implicit none

        type(aeros_vgrid_class), intent(in)    :: vg
        type(aeros_grid_class),  intent(in)    :: grd
        type(aeros_work_class),  intent(inout) :: wrk

        integer :: i, j, k, nlev

        real(wp) :: p_half(0:vg%nlev), p_full(vg%nlev), dp_lev(vg%nlev)
        real(wp) :: alpha(vg%nlev), dlnp(vg%nlev), cfac(vg%nlev)
        real(wp) :: uc(vg%nlev), vc(vg%nlev), tc(vg%nlev)
        real(wp) :: dc(vg%nlev), gc(vg%nlev)
        real(wp) :: cmass(vg%nlev), scum(0:vg%nlev), mflux(0:vg%nlev)
        real(wp) :: omga(vg%nlev), phi(vg%nlev)
        real(wp) :: ps, zeta_f, pgf, vadv_u, vadv_v, vadv_t, kin
        logical  :: ztop

        nlev = vg%nlev
        ztop = (vg%A(0) <= 0.0_wp .and. vg%B(0) <= 0.0_wp)

        !$omp parallel do collapse(2) schedule(static) &
        !$omp   private(i,j,k,p_half,p_full,dp_lev,alpha,dlnp,cfac,uc,vc,tc,dc,gc, &
        !$omp           cmass,scum,mflux,omga,phi,ps,zeta_f,pgf,vadv_u,vadv_v,vadv_t,kin)
        do j = 1, wrk%nlat
            do i = 1, wrk%nlon

                ps = exp(wrk%lnps_g(i,j))

                call aeros_vgrid_pressure(vg, ps, p_half, p_full, dp_lev)
                call aeros_vgrid_alpha(vg, p_half, alpha, dlnp)

                ! dlnp is MV in a layer whose upper interface is at zero
                ! pressure. Every term that uses it there is multiplied by
                ! B_(k-1/2) = 0 or S_(k-1) = 0, so zero is the correct value to
                ! substitute -- but substituting it explicitly is what keeps
                ! MV out of the arithmetic.
                if (ztop) dlnp(1) = 0.0_wp

                ! The shared coefficient: layer mean of B p_s / p.
                do k = 1, nlev
                    cfac(k) = (ps/dp_lev(k)) &
                                *(vg%B(k-1)*dlnp(k) + (vg%B(k) - vg%B(k-1))*alpha(k))
                end do

                ! Gather the column.
                do k = 1, nlev
                    uc(k) = wrk%u(i,j,k)
                    vc(k) = wrk%v(i,j,k)
                    tc(k) = wrk%t_g(i,j,k)
                    dc(k) = wrk%div_g(i,j,k)
                    gc(k) = uc(k)*wrk%dlnpsdx(i,j) + vc(k)*wrk%dlnpsdy(i,j)
                end do

                ! Layer mass divergence and its cumulative sum.
                scum(0) = 0.0_wp
                do k = 1, nlev
                    cmass(k) = dp_lev(k)*dc(k) &
                                + (vg%B(k) - vg%B(k-1))*ps*gc(k)
                    scum(k)  = scum(k-1) + cmass(k)
                end do

                ! Surface pressure tendency.
                wrk%dlnpsdt(i,j) = -scum(nlev)/ps

                ! Vertical mass flux at half levels. Zero at both ends by
                ! construction: B(0) = 0 gives M(0) = 0, and B(nlev) = 1 gives
                ! M(nlev) = S_N - S_N = 0. No air leaves through the top or the
                ! ground, exactly, whatever the winds do.
                do k = 0, nlev
                    mflux(k) = vg%B(k)*scum(nlev) - scum(k)
                end do

                ! omega/p.
                do k = 1, nlev
                    omga(k) = cfac(k)*gc(k) &
                                - (scum(k-1)*dlnp(k) + alpha(k)*cmass(k))/dp_lev(k)
                end do

                ! Geopotential.
                call aeros_hydrostatic(vg, wrk%phis(i,j), tc, p_half, phi)

                ! Assemble the right-hand sides.
                do k = 1, nlev

                    zeta_f = wrk%vor_g(i,j,k) + grd%coriolis(i,j)

                    vadv_u = vert_adv(uc, mflux, dp_lev, k, nlev)
                    vadv_v = vert_adv(vc, mflux, dp_lev, k, nlev)
                    vadv_t = vert_adv(tc, mflux, dp_lev, k, nlev)

                    pgf = R_d*tc(k)*cfac(k)

                    wrk%ae(i,j,k) =  zeta_f*vc(k) - vadv_u - pgf*wrk%dlnpsdx(i,j)
                    wrk%an(i,j,k) = -zeta_f*uc(k) - vadv_v - pgf*wrk%dlnpsdy(i,j)

                    kin = 0.5_wp*(uc(k)*uc(k) + vc(k)*vc(k))
                    wrk%ek(i,j,k) = phi(k) + kin

                    wrk%dtdt(i,j,k) = -(uc(k)*wrk%dtdx(i,j,k) + vc(k)*wrk%dtdy(i,j,k)) &
                                        - vadv_t + kappa*tc(k)*omga(k)

                end do

            end do
        end do
        !$omp end parallel do

        return

    end subroutine column_terms

    real(wp) function vert_adv(x, mflux, dp_lev, k, nlev) result(a)
        ! Vertical advection of a column variable, Simmons & Burridge form:
        !
        !   (eta_dot dX/deta)_k = (1/(2 dp_k))
        !         [ M_(k+1/2) (X_(k+1) - X_k) + M_(k-1/2) (X_k - X_(k-1)) ]
        !
        ! The end conditions need no special case: M is exactly zero at both
        ! boundaries, so the out-of-range differences are multiplied by zero.
        ! They are still guarded, because forming X(nlev+1) would be an
        ! out-of-bounds read regardless of the factor it is multiplied by.

        implicit none

        real(wp), intent(in) :: x(:), mflux(0:), dp_lev(:)
        integer,  intent(in) :: k, nlev

        real(wp) :: upper, lower

        lower = 0.0_wp
        if (k < nlev) lower = mflux(k)*(x(k+1) - x(k))

        upper = 0.0_wp
        if (k > 1) upper = mflux(k-1)*(x(k) - x(k-1))

        a = 0.5_wp*(lower + upper)/dp_lev(k)

        return

    end function vert_adv

end module aeros_tendency
