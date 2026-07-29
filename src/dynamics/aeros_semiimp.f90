module aeros_semiimp
    ! The semi-implicit gravity-wave solve.
    !
    ! Explicit leapfrog on the primitive equations is limited by the fastest
    ! GRAVITY wave, not by the wind. The external mode travels at ~300 m s-1,
    ! against a jet of ~50, so an explicit step at T31 would be ~600 s where the
    ! advective CFL permits ~2700. docs/design.md section 3.2 therefore
    ! specifies semi-implicit: treat the linearized gravity-wave terms
    ! implicitly and leave everything else explicit, which buys back the factor
    ! of ~4 and lets section 4's "dt ~ 30 min" hold.
    !
    ! === What is linearized ==================================================
    !
    ! About a MOTIONLESS, horizontally uniform, ISOTHERMAL basic state --
    ! vg%t_ref and vg%ps_ref, which aeros_vertical defines for exactly this
    ! purpose. Three operators, all built from A, B and t_ref together:
    !
    !   L_D(T',P')_k = -laplacian( sum_j G_kj T'_j + R_d T_ref,k P' )
    !   L_T(D)       = -tau D
    !   L_P(D)       = -nu . D                                 (P = ln p_s)
    !
    !   G_kj    = R_d alpha_k          (j = k)      the hydrostatic matrix,
    !           = R_d dlnp_j           (j > k)      i.e. Phi' = G T'
    !           = 0                    (j < k)
    !
    !   tau_kj  = kappa T_ref,k alpha_k             (j = k)
    !           = kappa T_ref,k dlnp_k dp_j / dp_k  (j < k)
    !           = 0                                 (j > k)
    !
    !   nu_k    = dp_ref,k / ps_ref
    !
    ! Every coefficient is evaluated at the REFERENCE pressures by calling the
    ! same aeros_vertical routines the nonlinear tendency calls. That is not
    ! tidiness: the correction form below only cancels if L is a term-for-term
    ! linearization of the discrete R, so the two must share their
    ! coefficients or the scheme is neither exact on the reference state nor
    ! reliably stable.
    !
    ! THE COEFFICIENT OF P' IS R_d T_ref, NOT R_d T_ref c_k -- even though the
    ! nonlinear pressure-gradient force in aeros_tendency carries c_k, and
    ! putting it here is the obvious-looking way to make the two agree. It is
    ! wrong. A surface-pressure perturbation also moves the half-level
    ! pressures, hence alpha_k and dlnp_j, so the GEOPOTENTIAL responds to P'
    ! as well:
    !
    !     d Phi_k / d P = R_d T_ref (1 - c_k)      (isothermal column)
    !
    ! which is precisely the identity Simmons & Burridge's alpha_k is
    ! constructed to satisfy, and the one tests/test_tendency.f90 measures at
    ! 3e-12. Adding it to the -R_d T_ref c_k of the pressure-gradient term
    ! leaves R_d T_ref exactly. Restoring c_k here would make L differ from the
    ! true linearization by R_d T_ref (1 - c_k) -- which in the upper layers of
    ! a hybrid coordinate is 100% of the term, c_k being ~0 there. The
    ! finite-difference check in tests/test_semiimp.f90 is what holds this
    ! down.
    !
    ! ISOTHERMAL IS LOAD-BEARING. The vertical advection of the basic-state
    ! temperature is proportional to T_ref,k+1 - T_ref,k, so it vanishes
    ! identically for a uniform t_ref and no term for it appears in tau. A
    ! non-isothermal reference needs that term, and aeros_semiimp_init refuses
    ! one rather than silently dropping it.
    !
    ! Warm rather than cold, per aeros_vertical: the scheme is stable when the
    ! reference temperature is at or above the true column temperature and
    ! unstable well below it, so the error is one-sided.
    !
    ! === The scheme ==========================================================
    !
    ! CORRECTION form, over a step h (h = 2 dt leapfrog, h = dt on the startup
    ! step). The linear terms are evaluated at a DECENTERED average
    ! X_a = alpha X^new + (1 - alpha) X^old (`si%alpha`; 0.5 = centered / neutral
    ! gravity waves, the default; 1.0 = backward-implicit / SpeedyWeather, which
    ! damps the fast gravity + computational modes in the solve):
    !
    !     X^new = X^old + h R(X^now) + h [ L(X_a) - L(X^now) ]
    !
    ! R is the FULL nonlinear right-hand side from aeros_tendency -- nothing is
    ! split out of it. L appears twice with opposite signs, so the correction
    ! vanishes identically whenever the state is steady, and in particular the
    ! resting balanced atmosphere that aeros_tendency's acceptance test builds
    ! stays at rest to the same precision with the integrator wrapped around it.
    ! That is the property the alternative ("split R into linear + residual")
    ! does not have unless the split is exact.
    !
    ! Eliminating T^new and P^new from the D equation leaves, per coefficient,
    !
    !   ( I + a^2 h^2 lambda_l W ) D^new = ( I - a(1-a) h^2 lambda_l W ) D^old
    !                                    + h [ R_D + L_D(T*,P*) - L_D(T^now,P^now) ]
    !
    !     a = alpha,  W  = G tau + R_d T_ref (x) nu           [nlev x nlev]
    !     lambda_l = l(l+1)/a^2
    !     T* = T^old + alpha h ( R_T - L_T(D^now) ),  P* likewise
    !
    ! (alpha = 1/2 gives the centered ( I +/- h^2/4 lambda W ) form, bit for bit.)
    ! then T^new and P^new follow by back-substitution with D_a = alpha D^new +
    ! (1-alpha) D^old. The matrix
    ! depends on l only through the scalar lambda_l, which is what makes the
    ! solve nlev x nlev per DEGREE rather than global -- the reason the core
    ! carries vorticity and divergence instead of u and v (section 3.2).
    !
    ! W's eigenvalues are g H_n, the equivalent depths of the vertical normal
    ! modes; sqrt of the largest is the external gravity-wave speed, which
    ! aeros_semiimp_gwspeed reports.
    !
    ! === The solve ===========================================================
    !
    ! (lmax+1) dense nlev x nlev LU factorizations, rebuilt whenever h changes
    ! (aeros_semiimp_set_step) and reused for every coefficient of that degree.
    ! At T42L20 that is 138 kB of factors and ~0.9 Mflop of triangular solves
    ! per timestep, against ~11 x nlev spectral transforms -- unmeasurable.
    !
    ! LU is in-tree rather than LAPACK: aeros links no linear-algebra library
    ! (config/common.mk), and a pivoted nlev x nlev factorization is forty
    ! lines. W is NOT symmetric (G tau is a product of two triangular matrices
    ! and the rank-one term is not symmetric either), so Cholesky is not
    ! available and partial pivoting is not optional.

    use aeros_defs,     only : dp, wp, wp_sh, io_unit_err, R_d, kappa, r_earth, &
                                aeros_spec_class
    use aeros_spectral, only : aeros_sht_class
    use aeros_vertical, only : aeros_vgrid_class, aeros_vgrid_pressure, &
                                aeros_vgrid_alpha
    use aeros_tendency, only : aeros_tend_class

    implicit none

    private

    type aeros_semiimp_class
        ! The linearized operators and their factorized Helmholtz matrices.

        integer  :: nlev = 0
        integer  :: lmax = 0
        real(wp) :: dt_step = 0.0_wp        ! the h the factors were built for [s]
        ! Decentering weight on X^(n+1): 0.5 = centered (neutral gravity waves),
        ! 1.0 = backward-implicit (SpeedyWeather; damps them). The factors depend
        ! on it (I + alpha^2 h^2 lambda W), so it is fixed at init.
        real(wp) :: alpha = 0.5_wp

        ! Reference-state column coefficients, at ps_ref.
        real(wp), allocatable :: dp_ref(:)     ! (nlev) layer thickness [Pa]
        real(wp), allocatable :: alpha_ref(:)  ! (nlev) [-]
        real(wp), allocatable :: dlnp_ref(:)   ! (nlev) [-]

        ! The linear operators.
        real(wp), allocatable :: gmat(:,:)     ! (nlev,nlev) Phi' = G T'
        real(wp), allocatable :: tmat(:,:)     ! (nlev,nlev) dT'/dt = -tau D
        real(wp), allocatable :: nu(:)         ! (nlev)      dP/dt  = -nu . D
        real(wp), allocatable :: wmat(:,:)     ! (nlev,nlev) G tau + R_d T_ref (x) nu

        ! Response of the momentum equation to ln(p_s): pcoef_k = R_d T_ref,k.
        ! See the module header for why there is no c_k here.
        real(wp), allocatable :: pcoef(:)      ! (nlev) [m2 s-2]

        ! LU factors of I + (h^2/4) lambda_l W, one per degree.
        real(wp), allocatable :: lu(:,:,:)     ! (nlev,nlev,0:lmax)
        integer,  allocatable :: piv(:,:)      ! (nlev,0:lmax)
    end type aeros_semiimp_class

    public :: aeros_semiimp_class
    public :: aeros_semiimp_init
    public :: aeros_semiimp_end
    public :: aeros_semiimp_set_step
    public :: aeros_semiimp_linear
    public :: aeros_semiimp_step
    public :: aeros_semiimp_gwspeed
    public :: aeros_semiimp_print

contains

    subroutine aeros_semiimp_init(si, vg, lmax, dt_step, alpha)
        ! Build the linear operators and factorize for a step of `dt_step`.
        !
        ! alpha is the decentering weight on X^(n+1) (optional, default 0.5 =
        ! centered). 1.0 = backward-implicit (SpeedyWeather). Stored on si and
        ! used by aeros_semiimp_set_step / _step; 0.5 reproduces the centered
        ! scheme bit for bit.

        implicit none

        type(aeros_semiimp_class), intent(inout) :: si
        type(aeros_vgrid_class),   intent(in)    :: vg
        integer,  intent(in) :: lmax
        real(wp), intent(in) :: dt_step        ! [s]
        real(wp), intent(in), optional :: alpha

        real(wp) :: p_half(0:vg%nlev), p_full(vg%nlev)
        real(wp) :: alpha_use
        integer  :: nlev, k, j, m
        logical  :: ztop

        nlev = vg%nlev

        alpha_use = 0.5_wp
        if (present(alpha)) alpha_use = alpha
        if (alpha_use < 0.5_wp .or. alpha_use > 1.0_wp) then
            write(io_unit_err,*) "aeros_semiimp_init:: error: alpha (decentering) must be in "// &
                                    "[0.5, 1.0], got ", alpha_use
            error stop 1
        end if

        if (any(vg%t_ref /= vg%t_ref(1))) then
            write(io_unit_err,*) "aeros_semiimp_init:: error: the reference temperature is not "// &
                                    "isothermal. The linearized thermodynamic operator omits the "// &
                                    "vertical advection of T_ref, which is exact only for a uniform "// &
                                    "profile. Add that term before using one."
            error stop 1
        end if

        call aeros_semiimp_end(si)

        si%nlev = nlev
        si%lmax = lmax
        si%alpha = alpha_use

        allocate(si%dp_ref(nlev), si%alpha_ref(nlev), si%dlnp_ref(nlev))
        allocate(si%gmat(nlev,nlev), si%tmat(nlev,nlev), si%wmat(nlev,nlev))
        allocate(si%nu(nlev), si%pcoef(nlev))
        allocate(si%lu(nlev,nlev,0:lmax), si%piv(nlev,0:lmax))

        ! === Reference column ================================================
        call aeros_vgrid_pressure(vg, vg%ps_ref, p_half, p_full, si%dp_ref)
        call aeros_vgrid_alpha(vg, p_half, si%alpha_ref, si%dlnp_ref)

        ! Same substitution, for the same reason, as aeros_tendency's column
        ! loop: with a zero-pressure top dlnp(1) is MV, and every term that
        ! would use it there is multiplied by B(0) = 0 or by an empty sum.
        ztop = (vg%A(0) <= 0.0_wp .and. vg%B(0) <= 0.0_wp)
        if (ztop) si%dlnp_ref(1) = 0.0_wp

        ! === G: the hydrostatic matrix =======================================
        ! Linearization of aeros_hydrostatic, whose upward integral gives
        !   Phi_k = Phi_s + alpha_k R T_k + sum_(j>k) R T_j dlnp_j
        si%gmat = 0.0_wp
        do k = 1, nlev
            si%gmat(k,k) = R_d*si%alpha_ref(k)
            do j = k+1, nlev
                si%gmat(k,j) = R_d*si%dlnp_ref(j)
            end do
        end do

        ! === tau: the thermodynamic matrix ===================================
        ! From kappa T (omega/p) with the winds and grad(ln p_s) set to zero,
        ! so (omega/p)_k = -(1/dp_k)[ dlnp_k sum_(j<k) dp_j D_j + alpha_k dp_k D_k ].
        si%tmat = 0.0_wp
        do k = 1, nlev
            si%tmat(k,k) = kappa*vg%t_ref(k)*si%alpha_ref(k)
            do j = 1, k-1
                si%tmat(k,j) = kappa*vg%t_ref(k)*si%dlnp_ref(k)*si%dp_ref(j)/si%dp_ref(k)
            end do
        end do

        ! === nu: the surface-pressure vector =================================
        do k = 1, nlev
            si%nu(k) = si%dp_ref(k)/vg%ps_ref
        end do

        ! === The pressure-gradient coupling and W ============================
        do k = 1, nlev
            si%pcoef(k) = R_d*vg%t_ref(k)
        end do

        do j = 1, nlev
            do k = 1, nlev
                si%wmat(k,j) = si%pcoef(k)*si%nu(j)
                do m = 1, nlev
                    si%wmat(k,j) = si%wmat(k,j) + si%gmat(k,m)*si%tmat(m,j)
                end do
            end do
        end do

        call aeros_semiimp_set_step(si, dt_step)

        return

    end subroutine aeros_semiimp_init

    subroutine aeros_semiimp_end(si)

        implicit none

        type(aeros_semiimp_class), intent(inout) :: si

        si%nlev = 0; si%lmax = 0; si%dt_step = 0.0_wp

        if (allocated(si%dp_ref))    deallocate(si%dp_ref)
        if (allocated(si%alpha_ref)) deallocate(si%alpha_ref)
        if (allocated(si%dlnp_ref))  deallocate(si%dlnp_ref)
        if (allocated(si%gmat))      deallocate(si%gmat)
        if (allocated(si%tmat))      deallocate(si%tmat)
        if (allocated(si%nu))        deallocate(si%nu)
        if (allocated(si%wmat))      deallocate(si%wmat)
        if (allocated(si%pcoef))     deallocate(si%pcoef)
        if (allocated(si%lu))        deallocate(si%lu)
        if (allocated(si%piv))       deallocate(si%piv)

        return

    end subroutine aeros_semiimp_end

    subroutine aeros_semiimp_set_step(si, dt_step)
        ! Factorize I + (h^2/4) lambda_l W for every degree.
        !
        ! Separate from init because h is not constant over a run: the leapfrog
        ! start-up step is dt and every step after it is 2 dt, so the factors
        ! are rebuilt at least once. Rebuilding is (lmax+1) nlev^3/3 flops --
        ! ~0.1 Mflop at T42L20, i.e. cheaper than one spectral transform --
        ! which is why h is stored on the type and never passed to the solve.
        ! A caller cannot then apply factors built for a different step.

        implicit none

        type(aeros_semiimp_class), intent(inout) :: si
        real(wp), intent(in) :: dt_step        ! [s]

        real(wp) :: hh, lam, a2
        integer  :: l, k, j

        if (dt_step <= 0.0_wp) then
            write(io_unit_err,*) "aeros_semiimp_set_step:: error: step must be > 0, got ", dt_step
            error stop 1
        end if

        si%dt_step = dt_step

        ! Helmholtz coefficient (alpha h)^2; alpha = 0.5 recovers the centered h^2/4.
        hh = si%alpha*si%alpha*dt_step*dt_step
        a2 = r_earth*r_earth

        do l = 0, si%lmax
            lam = real(l*(l+1), dp)/a2
            do j = 1, si%nlev
                do k = 1, si%nlev
                    si%lu(k,j,l) = hh*lam*si%wmat(k,j)
                end do
                si%lu(j,j,l) = si%lu(j,j,l) + 1.0_wp
            end do
            call lu_factor(si%lu(:,:,l), si%piv(:,l), si%nlev)
        end do

        return

    end subroutine aeros_semiimp_set_step

    subroutine aeros_semiimp_linear(si, sht, spec, ldiv, ltemp, lps)
        ! Apply L to a spectral state: the linearized tendencies of divergence,
        ! temperature and ln(p_s) that the scheme treats implicitly.
        !
        ! The integrator does not call this -- aeros_semiimp_step folds it into
        ! the column loop to avoid materializing three extra spectral arrays
        ! per step. It is public because the operator has to be checkable
        ! against the nonlinear tendency it claims to linearize, which is what
        ! tests/test_semiimp.f90 does, and because normal-mode initialization
        ! (should it ever be wanted) is written in terms of it.

        implicit none

        type(aeros_semiimp_class), intent(in)  :: si
        type(aeros_sht_class),     intent(in)  :: sht
        type(aeros_spec_class),    intent(in)  :: spec
        complex(wp_sh), intent(out) :: ldiv(:,:)   ! (nlm,nlev) [s-2]
        complex(wp_sh), intent(out) :: ltemp(:,:)  ! (nlm,nlev) [K s-1]
        complex(wp_sh), intent(out) :: lps(:)      ! (nlm) [s-1]

        real(wp) :: lam, a2
        integer  :: lm, l, k, nlev

        complex(wp_sh) :: dcol(si%nlev), tcol(si%nlev)
        complex(wp_sh) :: ldcol(si%nlev), ltcol(si%nlev), lpc

        nlev = si%nlev
        a2   = r_earth*r_earth

        !$omp parallel do schedule(static) &
        !$omp   private(lm,l,k,lam,dcol,tcol,ldcol,ltcol,lpc)
        do lm = 1, sht%nlm

            l   = sht%l_of_lm(lm)
            lam = real(l*(l+1), dp)/a2

            do k = 1, nlev
                dcol(k) = spec%div(lm,k)
                tcol(k) = spec%temp(lm,k)
            end do

            call lin_column(si, lam, dcol, tcol, spec%lnps(lm), ldcol, ltcol, lpc)

            do k = 1, nlev
                ldiv(lm,k)  = ldcol(k)
                ltemp(lm,k) = ltcol(k)
            end do
            lps(lm) = lpc

        end do
        !$omp end parallel do

        return

    end subroutine aeros_semiimp_linear

    subroutine lin_column(si, lam, dcol, tcol, p, ldcol, ltcol, lp)
        ! L for one spectral coefficient's column.

        implicit none

        type(aeros_semiimp_class), intent(in) :: si
        real(wp), intent(in) :: lam                    ! l(l+1)/a^2
        complex(wp_sh), intent(in)  :: dcol(:), tcol(:), p
        complex(wp_sh), intent(out) :: ldcol(:), ltcol(:), lp

        complex(wp_sh) :: acc
        integer :: k, j

        call lin_div(si, lam, tcol, p, ldcol)

        ! L_T(D) = -tau D. tau is lower triangular.
        do k = 1, si%nlev
            acc = (0.0_wp_sh, 0.0_wp_sh)
            do j = 1, k
                acc = acc + si%tmat(k,j)*dcol(j)
            end do
            ltcol(k) = -acc
        end do

        ! L_P(D) = -nu . D.
        acc = (0.0_wp_sh, 0.0_wp_sh)
        do k = 1, si%nlev
            acc = acc + si%nu(k)*dcol(k)
        end do
        lp = -acc

        return

    end subroutine lin_column

    subroutine lin_div(si, lam, tcol, p, ldcol)
        ! L_D alone: -laplacian( G T' + R_d T_ref P' ), which in spectral space
        ! is +lambda_l times the bracket. G is upper triangular.

        implicit none

        type(aeros_semiimp_class), intent(in) :: si
        real(wp), intent(in) :: lam
        complex(wp_sh), intent(in)  :: tcol(:), p
        complex(wp_sh), intent(out) :: ldcol(:)

        complex(wp_sh) :: acc
        integer :: k, j

        do k = 1, si%nlev
            acc = si%pcoef(k)*p
            do j = k, si%nlev
                acc = acc + si%gmat(k,j)*tcol(j)
            end do
            ldcol(k) = lam*acc
        end do

        return

    end subroutine lin_div

    subroutine aeros_semiimp_step(si, sht, old, now, tnd, new)
        ! Advance divergence, temperature and ln(p_s) over one step.
        !
        ! Vorticity is NOT touched: it has no linearized gravity-wave part, so
        ! it is a plain leapfrog and belongs with the rest of the explicit
        ! update in aeros_timestep. This routine is the gravity-wave solve and
        ! nothing else.
        !
        ! The step length is si%dt_step, set by aeros_semiimp_set_step.

        implicit none

        type(aeros_semiimp_class), intent(in) :: si
        type(aeros_sht_class),     intent(in) :: sht
        type(aeros_spec_class),    intent(in) :: old    ! X^(n-1)
        type(aeros_spec_class),    intent(in) :: now    ! X^n
        type(aeros_tend_class),    intent(in) :: tnd    ! R(X^n)
        type(aeros_spec_class), intent(inout) :: new    ! X^(n+1), out

        real(wp) :: h, a1, astep, hh_old, lam, a2
        integer  :: lm, l, k, j, nlev

        complex(wp_sh) :: dold(si%nlev), dnow(si%nlev), told(si%nlev), tnow(si%nlev)
        complex(wp_sh) :: rd(si%nlev), rt(si%nlev)
        complex(wp_sh) :: ltnow(si%nlev), ldnow(si%nlev), ldstar(si%nlev)
        complex(wp_sh) :: tstar(si%nlev), rhs(si%nlev), dbar(si%nlev)
        complex(wp_sh) :: pold, pnow, rp, lpnow, pstar, lpbar
        complex(wp_sh) :: acc

        nlev = si%nlev
        h    = si%dt_step
        ! Decentering coefficients (alpha = 0.5 -> the centered scheme, bit for bit).
        a1     = si%alpha                     ! weight on X^(n+1) in the implicit average
        astep  = a1*h                         ! starred explicit half-step (was 0.5 h)
        hh_old = a1*(1.0_wp - a1)*h*h          ! rhs D^(n-1) damping (was h^2/4)
        a2   = r_earth*r_earth

        !$omp parallel do schedule(static) &
        !$omp   private(lm,l,k,j,lam,dold,dnow,told,tnow,rd,rt,ltnow,ldnow,ldstar, &
        !$omp           tstar,rhs,dbar,pold,pnow,rp,lpnow,pstar,lpbar,acc)
        do lm = 1, sht%nlm

            l   = sht%l_of_lm(lm)
            lam = real(l*(l+1), dp)/a2

            ! Gather the column. The spectral arrays are (nlm, nlev), so a
            ! single coefficient's column is strided -- the same trade
            ! aeros_tendency's column loop makes, and for the same reason.
            do k = 1, nlev
                dold(k) = old%div(lm,k)
                dnow(k) = now%div(lm,k)
                told(k) = old%temp(lm,k)
                tnow(k) = now%temp(lm,k)
                rd(k)   = tnd%div(lm,k)
                rt(k)   = tnd%temp(lm,k)
            end do
            pold = old%lnps(lm)
            pnow = now%lnps(lm)
            rp   = tnd%lnps(lm)

            ! --- Linear terms at the current time level, L(X^n).
            call lin_column(si, lam, dnow, tnow, pnow, ldnow, ltnow, lpnow)

            ! --- The starred state: T and P advanced alpha of a step explicitly.
            do k = 1, nlev
                tstar(k) = told(k) + astep*(rt(k) - ltnow(k))
            end do
            pstar = pold + astep*(rp - lpnow)

            call lin_div(si, lam, tstar, pstar, ldstar)

            ! --- Helmholtz right-hand side, then the solve.
            do k = 1, nlev
                acc = (0.0_wp_sh, 0.0_wp_sh)
                do j = 1, nlev
                    acc = acc + si%wmat(k,j)*dold(j)
                end do
                rhs(k) = dold(k) - hh_old*lam*acc &
                            + h*(rd(k) + ldstar(k) - ldnow(k))
            end do

            call lu_solve(si%lu(:,:,l), si%piv(:,l), rhs, nlev)

            do k = 1, nlev
                new%div(lm,k) = rhs(k)
                dbar(k)       = a1*rhs(k) + (1.0_wp - a1)*dold(k)
            end do

            ! --- Back-substitute for T and P.
            do k = 1, nlev
                acc = (0.0_wp_sh, 0.0_wp_sh)
                do j = 1, k
                    acc = acc + si%tmat(k,j)*dbar(j)
                end do
                new%temp(lm,k) = told(k) + h*(rt(k) - acc - ltnow(k))
            end do

            acc = (0.0_wp_sh, 0.0_wp_sh)
            do k = 1, nlev
                acc = acc + si%nu(k)*dbar(k)
            end do
            lpbar = -acc
            new%lnps(lm) = pold + h*(rp + lpbar - lpnow)

        end do
        !$omp end parallel do

        return

    end subroutine aeros_semiimp_step

    real(wp) function aeros_semiimp_gwspeed(si) result(cmax)
        ! Fastest linear gravity-wave phase speed [m s-1], sqrt of the largest
        ! eigenvalue of W.
        !
        ! W's eigenvalues are g H_n, the equivalent depths of the vertical
        ! normal modes; the largest belongs to the external (Lamb) mode and is
        ! the one that sets the explicit stability limit,
        !
        !     dt_explicit < a / ( cmax sqrt(lmax(lmax+1)) )
        !
        ! which is the number that says whether the semi-implicit solve is
        ! earning its keep. Power iteration: W is not symmetric, but its
        ! spectrum is real and positive and the dominant eigenvalue is well
        ! separated (the external mode is ~10 km of equivalent depth against
        ! ~100 m for the first internal one), so a few dozen iterations
        ! converge to well past the precision this diagnostic needs.

        implicit none

        type(aeros_semiimp_class), intent(in) :: si

        real(wp) :: x(si%nlev), y(si%nlev)
        real(wp) :: lam, nrm
        integer  :: it, k, j

        x   = 1.0_wp/sqrt(real(si%nlev, wp))
        lam = 0.0_wp

        do it = 1, 200
            do k = 1, si%nlev
                y(k) = 0.0_wp
                do j = 1, si%nlev
                    y(k) = y(k) + si%wmat(k,j)*x(j)
                end do
            end do
            nrm = sqrt(sum(y*y))
            if (nrm <= 0.0_wp) exit
            x = y/nrm
            ! Rayleigh-style estimate from the un-normalized image, which is
            ! the eigenvalue itself once x has converged to the eigenvector.
            lam = dot_product(x, y)
        end do

        cmax = sqrt(max(lam, 0.0_wp))

        return

    end function aeros_semiimp_gwspeed

    subroutine aeros_semiimp_print(si, io_unit)
        ! Report what the solve was built for.
        !
        ! The explicit limit is printed next to the configured step because
        ! their ratio is the whole justification for this module, and because a
        ! run log that does not record it cannot answer "was dt safe?" after
        ! the fact.

        implicit none

        type(aeros_semiimp_class), intent(in) :: si
        integer, intent(in), optional :: io_unit

        real(wp) :: cmax, dt_exp
        integer  :: iou

        iou = 6
        if (present(io_unit)) iou = io_unit

        cmax   = aeros_semiimp_gwspeed(si)
        dt_exp = r_earth/(cmax*sqrt(real(si%lmax*(si%lmax+1), wp)))

        write(iou,*) ""
        write(iou,"(a)")         " == semi-implicit solve =="
        write(iou,"(a,f9.1,a)")  "   fastest gravity wave       ", cmax, " m s-1"
        write(iou,"(a,f9.1,a)")  "   explicit stability limit   ", dt_exp, " s"
        write(iou,"(a,f9.1,a)")  "   step being factorized      ", si%dt_step, " s"
        write(iou,"(a,f9.1)")    "   ratio (what this buys)     ", si%dt_step/dt_exp
        write(iou,*) ""

        return

    end subroutine aeros_semiimp_print

    ! === Dense LU, in-tree ===================================================

    subroutine lu_factor(a, piv, n)
        ! Doolittle LU with partial pivoting, in place. `piv(k)` is the row
        ! swapped into position k.
        !
        ! Partial pivoting is not optional here: W is not symmetric and, more
        ! to the point, its rows differ by orders of magnitude (dp_k spans a
        ! factor of ~50 between the top layer and the bottom one), so the
        ! unpivoted factorization loses digits in exactly the configurations
        ! the model is meant to run.

        implicit none

        real(wp), intent(inout) :: a(:,:)
        integer,  intent(out)   :: piv(:)
        integer,  intent(in)    :: n

        real(wp) :: amax, tmp
        integer  :: k, i, j, ip

        do k = 1, n

            ip   = k
            amax = abs(a(k,k))
            do i = k+1, n
                if (abs(a(i,k)) > amax) then
                    amax = abs(a(i,k))
                    ip   = i
                end if
            end do

            if (amax <= 0.0_wp) then
                write(io_unit_err,*) "aeros_semiimp:: error: singular Helmholtz matrix at column ", k
                error stop 1
            end if

            piv(k) = ip
            if (ip /= k) then
                do j = 1, n
                    tmp     = a(k,j)
                    a(k,j)  = a(ip,j)
                    a(ip,j) = tmp
                end do
            end if

            do i = k+1, n
                a(i,k) = a(i,k)/a(k,k)
                do j = k+1, n
                    a(i,j) = a(i,j) - a(i,k)*a(k,j)
                end do
            end do

        end do

        return

    end subroutine lu_factor

    subroutine lu_solve(a, piv, b, n)
        ! Solve L U x = P b in place, for a COMPLEX right-hand side against the
        ! real factors -- the spectral coefficients are complex, the operators
        ! are not.

        implicit none

        real(wp),       intent(in)    :: a(:,:)
        integer,        intent(in)    :: piv(:)
        complex(wp_sh), intent(inout) :: b(:)
        integer,        intent(in)    :: n

        complex(wp_sh) :: tmp
        integer :: k, j

        do k = 1, n
            if (piv(k) /= k) then
                tmp     = b(k)
                b(k)    = b(piv(k))
                b(piv(k)) = tmp
            end if
            do j = k+1, n
                b(j) = b(j) - a(j,k)*b(k)
            end do
        end do

        do k = n, 1, -1
            do j = k+1, n
                b(k) = b(k) - a(k,j)*b(j)
            end do
            b(k) = b(k)/a(k,k)
        end do

        return

    end subroutine lu_solve

end module aeros_semiimp
