module aeros_vertical
    ! The vertical coordinate, and the hydrostatic relation on it.
    !
    ! HYBRID SIGMA-PRESSURE, following Simmons & Burridge (1981, MWR 109, 758),
    ! which is what ECMWF, ECHAM and most spectral models use. Half-level
    ! pressures are
    !
    !     p_(k+1/2) = A_(k+1/2) + B_(k+1/2) * p_s
    !
    ! so the coordinate is terrain-following near the ground (B -> 1) and
    ! becomes pure pressure aloft (B -> 0). Setting sigma_t = 0 gives A = 0 and
    ! B = sigma everywhere, i.e. **pure sigma**, exactly -- which is the
    ! configuration Held-Suarez is defined on (M1.5), so the same code runs the
    ! validation benchmark and the production model.
    !
    ! WHY HYBRID FROM THE START, when M1 only needs sigma: docs/design.md
    ! section 4 wants it for the vertical CFL problem over ice-sheet margins,
    ! where sigma-coordinate vertical velocity picks up u.grad(h_s) -- with
    ! |grad h| ~ 0.02 and u ~ 20 m/s that is an effective w ~ 0.4 m/s, and it
    ! is a property of the coordinate, not of the physics. Retrofitting a
    ! vertical coordinate after the dynamics are written is miserable; adding
    ! the two coefficient arrays now costs nothing.
    !
    ! INDEX CONVENTION, fixed here and assumed by every dynamics module:
    !
    !     full levels  k = 1 .. nlev,     k = 1 at the MODEL TOP
    !     half levels  j = 0 .. nlev,     j = 0 at the top, j = nlev at surface
    !     p_half(0) = p_top (may be 0),   p_half(nlev) = p_s
    !     dp(k) = p_half(k) - p_half(k-1) > 0
    !
    ! Top-down, matching the latitude ordering (north first) and ECMWF's
    ! convention. Every loop that walks the column downward runs k = 1, nlev.
    !
    ! Routines here are COLUMN routines. The model never materializes p_half
    ! for the whole globe: pressures are recomputed inside the threaded column
    ! loop where they are used, which is both cheaper in memory traffic and the
    ! natural shape for the physics (M2).

    use aeros_defs, only : dp, wp, io_unit_err, MV, R_d, grav, p0
    use nml,        only : nml_read

    implicit none

    private

    type aeros_vgrid_class
        ! The vertical coordinate, plus the reference state the semi-implicit
        ! solve linearizes about.

        integer :: nlev = 0

        ! Hybrid coefficients at half levels, (0:nlev). These are the vertical
        ! coordinate: everything else here is derived from them, and an
        ! explicit A/B table from another model can be dropped straight in
        ! (see aeros_vgrid_init's `a_half`/`b_half` arguments).
        real(wp), allocatable :: A(:)          ! [Pa]
        real(wp), allocatable :: B(:)          ! [-]

        ! The sigma profile the coefficients were generated from, i.e. the
        ! level placement at p_s = ps_ref. Diagnostic, but worth keeping: it is
        ! what one actually looks at when judging a level distribution.
        real(wp), allocatable :: sigma_half(:) ! (0:nlev)
        real(wp), allocatable :: sigma_full(:) ! (1:nlev)

        ! === Reference state =================================================
        !
        ! The semi-implicit scheme (M1.4) splits the gravity-wave terms about a
        ! motionless, horizontally uniform basic state and treats that part
        ! implicitly. It is defined here rather than in the timestepping module
        ! because it is a property of the vertical discretization -- the
        ! structure matrix it produces is built from A, B and t_ref together --
        ! and because the same reference profile is wanted for diagnostics.
        !
        ! An ISOTHERMAL reference is the standard choice and is deliberate: the
        ! scheme is stable when the reference temperature is at or above the
        ! true column temperature, and unstable when it is far below, so a warm
        ! isothermal profile is the safe side of an asymmetric error.
        real(wp) :: ps_ref = 0.0_wp            ! reference surface pressure [Pa]
        real(wp), allocatable :: t_ref(:)      ! reference temperature [K], (nlev)
        real(wp), allocatable :: p_ref_half(:) ! (0:nlev) [Pa]
        real(wp), allocatable :: p_ref_full(:) ! (1:nlev) [Pa]
    end type aeros_vgrid_class

    public :: aeros_vgrid_class
    public :: aeros_vgrid_init
    public :: aeros_vgrid_end
    public :: aeros_vgrid_load
    public :: aeros_vgrid_pressure
    public :: aeros_vgrid_alpha
    public :: aeros_hydrostatic
    public :: aeros_vgrid_print

contains

    subroutine aeros_vgrid_load(vg, nlev, filename, group)
        ! Build the vertical grid from a namelist group (default `aeros_vert`).

        implicit none

        type(aeros_vgrid_class), intent(inout) :: vg
        integer, intent(in) :: nlev
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: group

        character(len=56) :: nml_group
        real(wp) :: stretch_a, stretch_r, sigma_t, p_top, ps_ref, t_ref

        nml_group = "aeros_vert"
        if (present(group)) nml_group = trim(group)

        call nml_read(filename, nml_group, "stretch_a", stretch_a)
        call nml_read(filename, nml_group, "stretch_r", stretch_r)
        call nml_read(filename, nml_group, "sigma_t",   sigma_t)
        call nml_read(filename, nml_group, "p_top",     p_top)
        call nml_read(filename, nml_group, "ps_ref",    ps_ref)
        call nml_read(filename, nml_group, "t_ref",     t_ref)

        call aeros_vgrid_init(vg, nlev, stretch_a=stretch_a, stretch_r=stretch_r, &
                                sigma_t=sigma_t, p_top=p_top, ps_ref=ps_ref, t_ref=t_ref)

        return

    end subroutine aeros_vgrid_load

    subroutine aeros_vgrid_init(vg, nlev, stretch_a, stretch_r, sigma_t, p_top, &
                                    ps_ref, t_ref, a_half, b_half)
        ! Construct the hybrid coefficients.
        !
        ! Either generated from the stretching parameters, or taken verbatim
        ! from `a_half`/`b_half` when a published table is wanted (ECHAM L19,
        ! IFS L31, ...). The generator exists so that nlev is a free parameter
        ! during development; a production configuration may well end up
        ! pinning a table instead, and nothing downstream can tell the
        ! difference.

        implicit none

        type(aeros_vgrid_class), intent(inout) :: vg
        integer,  intent(in) :: nlev
        real(wp), intent(in), optional :: stretch_a  ! 1 = uniform in sigma
        real(wp), intent(in), optional :: stretch_r  ! stretching exponent
        real(wp), intent(in), optional :: sigma_t    ! hybrid transition; 0 = pure sigma
        real(wp), intent(in), optional :: p_top      ! model top pressure [Pa]
        real(wp), intent(in), optional :: ps_ref     ! reference surface pressure [Pa]
        real(wp), intent(in), optional :: t_ref      ! isothermal reference temperature [K]
        real(wp), intent(in), optional :: a_half(0:) ! explicit table [Pa]
        real(wp), intent(in), optional :: b_half(0:) ! explicit table [-]

        real(wp) :: a_, r_, st_, ptop_, psref_, tref_
        real(wp) :: x, sig, sig_top
        integer  :: j, k

        ! Defaults are the PRODUCTION configuration: a stretched hybrid with a
        ! 10 hPa top. Held-Suarez (M1.5) overrides them to uniform pure sigma.
        !
        ! stretch_a = 0.4, stretch_r = 2.0 chosen by scanning the level
        ! distribution at L20: it gives 5 full levels below 850 hPa, which
        ! meets docs/design.md section 4.1's "3-4 levels below 850 hPa" with
        ! margin, while keeping 3 levels above 200 hPa. Pushing further down
        ! (a = 0.2) buys 8 levels below 850 but leaves only 2 above 200 hPa,
        ! and section 4.1 also warns that level count moves the jet latitude
        ! and hence the storm track -- so starving the upper troposphere to
        ! feed the boundary layer trades one first-order error for another.
        a_     = 0.4_wp
        r_     = 2.0_wp
        st_    = 0.2_wp
        ptop_  = 1000.0_wp
        psref_ = p0
        tref_  = 300.0_wp

        if (present(stretch_a)) a_     = stretch_a
        if (present(stretch_r)) r_     = stretch_r
        if (present(sigma_t))   st_    = sigma_t
        if (present(p_top))     ptop_  = p_top
        if (present(ps_ref))    psref_ = ps_ref
        if (present(t_ref))     tref_  = t_ref

        if (nlev < 2) then
            write(io_unit_err,*) "aeros_vgrid_init:: error: nlev must be >= 2, got ", nlev
            error stop 1
        end if
        if (a_ <= 0.0_wp .or. a_ > 1.0_wp) then
            write(io_unit_err,*) "aeros_vgrid_init:: error: stretch_a must be in (0,1], got ", a_
            error stop 1
        end if
        if (st_ < 0.0_wp .or. st_ >= 1.0_wp) then
            write(io_unit_err,*) "aeros_vgrid_init:: error: sigma_t must be in [0,1), got ", st_
            error stop 1
        end if

        call aeros_vgrid_end(vg)

        vg%nlev   = nlev
        vg%ps_ref = psref_

        allocate(vg%A(0:nlev), vg%B(0:nlev))
        allocate(vg%sigma_half(0:nlev), vg%sigma_full(nlev))
        allocate(vg%t_ref(nlev))
        allocate(vg%p_ref_half(0:nlev), vg%p_ref_full(nlev))

        if (present(a_half) .and. present(b_half)) then

            if (size(a_half) /= nlev+1 .or. size(b_half) /= nlev+1) then
                write(io_unit_err,*) "aeros_vgrid_init:: error: a_half/b_half must have nlev+1 = ", &
                                        nlev+1, " entries, got ", size(a_half), size(b_half)
                error stop 1
            end if
            ! Verbatim. A supplied table is the authority on the coordinate,
            ! so it is VALIDATED and never adjusted -- silently rewriting a
            ! published level set would make the model's coordinate differ from
            ! the one its provenance claims.
            if (b_half(nlev) /= 1.0_wp .or. a_half(nlev) /= 0.0_wp) then
                write(io_unit_err,*) "aeros_vgrid_init:: error: table must end at the surface, "// &
                                        "A(nlev)=0 and B(nlev)=1; got ", a_half(nlev), b_half(nlev)
                error stop 1
            end if
            if (b_half(0) /= 0.0_wp) then
                write(io_unit_err,*) "aeros_vgrid_init:: error: table must start pure-pressure, "// &
                                        "B(0)=0; got ", b_half(0)
                error stop 1
            end if

            vg%A = a_half
            vg%B = b_half
            do j = 0, nlev
                vg%sigma_half(j) = (vg%A(j)/psref_) + vg%B(j)
            end do

        else

            ! Level placement in sigma. The spacing is
            !
            !     dsigma/dx = a + r(1-a)(1-x)^(r-1),    x = j/nlev
            !
            ! so it is proportional to `a` at the surface and to a + r(1-a) at
            ! the top: SMALL `a` means fine resolution near the ground, which
            ! is what docs/design.md section 4.1 asks for ("3-4 levels below
            ! 850 hPa ... if polar performance is the point of the model, put
            ! them there before adding anything aloft"). The polar winter
            ! boundary layer is 50-200 m deep with 10-25 K inversions, and a
            ! single bulk layer cannot hold one.
            !
            ! a = 1 collapses this to sigma = x, UNIFORM spacing -- the
            ! Held-Suarez configuration (20 evenly spaced sigma levels).
            sig_top = ptop_/psref_
            do j = 0, nlev
                x   = real(j, wp)/real(nlev, wp)
                sig = 1.0_wp - (a_*(1.0_wp - x) + (1.0_wp - a_)*(1.0_wp - x)**r_)
                ! Squeeze onto [sigma_top, 1] so a non-zero model top is honoured.
                vg%sigma_half(j) = sig_top + (1.0_wp - sig_top)*sig
            end do

            ! Hybrid ramp. B is zero above the transition and rises linearly to
            ! 1 at the surface; A takes up the remainder, so that at
            ! p_s = ps_ref the half-level pressures are exactly sigma*ps_ref
            ! whatever sigma_t is.
            !
            ! sigma_t = 0 gives B = sigma and A = 0 identically: PURE SIGMA,
            ! bit-exact, not merely close. tests/test_vertical.f90 asserts it.
            do j = 0, nlev
                if (st_ <= 0.0_wp) then
                    vg%B(j) = vg%sigma_half(j)
                else if (vg%sigma_half(j) <= st_) then
                    vg%B(j) = 0.0_wp
                else
                    vg%B(j) = (vg%sigma_half(j) - st_)/(1.0_wp - st_)
                end if
                vg%A(j) = (vg%sigma_half(j) - vg%B(j))*psref_
            end do

            ! Pin the generated endpoints exactly rather than trusting the
            ! arithmetic: the surface half level MUST be p_s (B = 1, A = 0) or
            ! mass is not conserved, and sum(dp) = p_s - p_top is the identity
            ! everything downstream leans on. Inside the generator branch only
            ! -- a supplied table was validated above, not corrected.
            vg%B(nlev) = 1.0_wp
            vg%A(nlev) = 0.0_wp
            vg%B(0)    = 0.0_wp
            vg%A(0)    = ptop_
            vg%sigma_half(nlev) = 1.0_wp
            vg%sigma_half(0)    = ptop_/psref_

        end if

        do k = 1, nlev
            vg%sigma_full(k) = 0.5_wp*(vg%sigma_half(k) + vg%sigma_half(k-1))
        end do

        ! Reference state.
        vg%t_ref = tref_
        call aeros_vgrid_pressure(vg, psref_, vg%p_ref_half, vg%p_ref_full)

        call check_monotonic(vg)

        return

    end subroutine aeros_vgrid_init

    subroutine aeros_vgrid_end(vg)

        implicit none

        type(aeros_vgrid_class), intent(inout) :: vg

        vg%nlev   = 0
        vg%ps_ref = 0.0_wp

        if (allocated(vg%A))          deallocate(vg%A)
        if (allocated(vg%B))          deallocate(vg%B)
        if (allocated(vg%sigma_half)) deallocate(vg%sigma_half)
        if (allocated(vg%sigma_full)) deallocate(vg%sigma_full)
        if (allocated(vg%t_ref))      deallocate(vg%t_ref)
        if (allocated(vg%p_ref_half)) deallocate(vg%p_ref_half)
        if (allocated(vg%p_ref_full)) deallocate(vg%p_ref_full)

        return

    end subroutine aeros_vgrid_end

    subroutine check_monotonic(vg)
        ! A non-monotonic coordinate produces negative layer masses, which show
        ! up much later as a mysterious instability. Fail at init instead.

        implicit none

        type(aeros_vgrid_class), intent(in) :: vg

        integer :: j

        do j = 1, vg%nlev
            if (vg%A(j) + vg%B(j)*vg%ps_ref <= vg%A(j-1) + vg%B(j-1)*vg%ps_ref) then
                write(io_unit_err,*) "aeros_vgrid_init:: error: non-monotonic vertical coordinate"
                write(io_unit_err,*) "  half level ", j-1, " p = ", vg%p_ref_half(j-1)
                write(io_unit_err,*) "  half level ", j,   " p = ", vg%p_ref_half(j)
                error stop 1
            end if
        end do

        return

    end subroutine check_monotonic

    subroutine aeros_vgrid_pressure(vg, ps, p_half, p_full, dp_lev)
        ! Half- and full-level pressures for one column, given surface pressure.
        !
        ! The full-level pressure is the arithmetic mean of the bounding half
        ! levels. That is a definition, not an approximation to something else:
        ! the hydrostatic integral in aeros_hydrostatic does NOT use it (it
        ! uses alpha, below), so p_full is only ever a coordinate for physics
        ! and output. Simmons & Burridge's implicit alternative buys nothing
        ! here and costs a log per level.

        implicit none

        type(aeros_vgrid_class), intent(in)  :: vg
        real(wp), intent(in)  :: ps                ! surface pressure [Pa]
        real(wp), intent(out) :: p_half(0:)        ! (0:nlev) [Pa]
        real(wp), intent(out) :: p_full(:)         ! (nlev) [Pa]
        real(wp), intent(out), optional :: dp_lev(:) ! (nlev) layer thickness [Pa]

        integer :: j, k

        do j = 0, vg%nlev
            p_half(j) = vg%A(j) + vg%B(j)*ps
        end do

        do k = 1, vg%nlev
            p_full(k) = 0.5_wp*(p_half(k) + p_half(k-1))
        end do

        if (present(dp_lev)) then
            do k = 1, vg%nlev
                dp_lev(k) = p_half(k) - p_half(k-1)
            end do
        end if

        return

    end subroutine aeros_vgrid_pressure

    subroutine aeros_vgrid_alpha(vg, p_half, alpha, dlnp)
        ! Simmons & Burridge's layer coefficients.
        !
        !   dlnp(k)  = ln( p_(k+1/2) / p_(k-1/2) )       [their delta ln p]
        !   alpha(k) = 1 - (p_(k-1/2)/dp_k) * dlnp(k)
        !
        ! alpha is the fraction of the layer's geopotential thickness that lies
        ! between the lower interface and the full level, and it is what makes
        ! the discrete hydrostatic operator the exact adjoint of the discrete
        ! pressure-gradient term -- which is what conserves energy in the
        ! discretization. Do NOT replace it with 0.5 or with a p_full-based
        ! expression, however similar the numbers look.
        !
        ! TOP LAYER: when p_top = 0 both dlnp(1) and the alpha expression are
        ! singular. Simmons & Burridge set alpha(1) = ln 2, the limit of the
        ! layer-mean geopotential for a layer whose top is at zero pressure.
        ! dlnp(1) is then never used by the hydrostatic integral above the
        ! first full level, and is returned as MV so that a caller which does
        ! use it fails loudly rather than silently propagating an infinity.

        implicit none

        type(aeros_vgrid_class), intent(in)  :: vg
        real(wp), intent(in)  :: p_half(0:)     ! (0:nlev) [Pa]
        real(wp), intent(out) :: alpha(:)       ! (nlev) [-]
        real(wp), intent(out), optional :: dlnp(:)  ! (nlev) [-]

        real(wp) :: dpk, lr
        integer  :: k

        do k = 1, vg%nlev
            dpk = p_half(k) - p_half(k-1)
            if (p_half(k-1) <= 0.0_wp) then
                alpha(k) = log(2.0_wp)
                if (present(dlnp)) dlnp(k) = MV
            else
                lr       = log(p_half(k)/p_half(k-1))
                alpha(k) = 1.0_wp - (p_half(k-1)/dpk)*lr
                if (present(dlnp)) dlnp(k) = lr
            end if
        end do

        return

    end subroutine aeros_vgrid_alpha

    subroutine aeros_hydrostatic(vg, phis, temp, p_half, phi_full, phi_half)
        ! Geopotential from the discrete hydrostatic relation, one column.
        !
        !   Phi_(k-1/2) = Phi_(k+1/2) + R T_k ln( p_(k+1/2) / p_(k-1/2) )
        !   Phi_k       = Phi_(k+1/2) + alpha_k R T_k
        !
        ! integrated UPWARD from the surface, where Phi_(nlev+1/2) = Phi_s.
        !
        ! `temp` should be the VIRTUAL temperature once moisture exists (M2);
        ! at M1 the core is dry and the two coincide. The argument is named
        ! `temp` rather than `tv` so that the moist change is a visible edit at
        ! every call site rather than a silent reinterpretation.
        !
        ! Phi_half(0) is left as MV when the model top is at zero pressure: the
        ! geopotential there is genuinely infinite, and nothing needs it.

        implicit none

        type(aeros_vgrid_class), intent(in)  :: vg
        real(wp), intent(in)  :: phis           ! surface geopotential [m2 s-2]
        real(wp), intent(in)  :: temp(:)        ! (nlev) temperature [K]
        real(wp), intent(in)  :: p_half(0:)     ! (0:nlev) [Pa]
        real(wp), intent(out) :: phi_full(:)    ! (nlev) [m2 s-2]
        real(wp), intent(out), optional :: phi_half(0:)

        real(wp) :: alpha(vg%nlev)
        real(wp) :: ph(0:vg%nlev)
        integer  :: k

        call aeros_vgrid_alpha(vg, p_half, alpha)

        ph(vg%nlev) = phis

        do k = vg%nlev, 1, -1
            phi_full(k) = ph(k) + alpha(k)*R_d*temp(k)
            if (p_half(k-1) > 0.0_wp) then
                ph(k-1) = ph(k) + R_d*temp(k)*log(p_half(k)/p_half(k-1))
            else
                ph(k-1) = MV
            end if
        end do

        if (present(phi_half)) phi_half = ph

        return

    end subroutine aeros_hydrostatic

    subroutine aeros_vgrid_print(vg, io_unit)
        ! Print the level distribution at the reference surface pressure.
        !
        ! Worth printing at every run start: docs/design.md section 4.1 makes
        ! level count and placement a first-order control on the jet latitude
        ! and hence the storm track, so a run log that does not record them
        ! cannot be reproduced from its namelist alone once the generator
        ! changes.

        implicit none

        type(aeros_vgrid_class), intent(in) :: vg
        integer, intent(in), optional :: io_unit

        integer :: iou, k, n_below_850

        iou = 6
        if (present(io_unit)) iou = io_unit

        n_below_850 = count(vg%p_ref_full > 85000.0_wp)

        write(iou,*) ""
        write(iou,"(a,i0,a)") " == vertical grid (nlev = ", vg%nlev, ") =="
        write(iou,"(a,f9.2,a)") "   reference surface pressure ", vg%ps_ref/100.0_wp, " hPa"
        write(iou,"(a,f9.2,a)") "   model top                  ", vg%p_ref_half(0)/100.0_wp, " hPa"
        write(iou,"(a,i0)")     "   full levels below 850 hPa   ", n_below_850
        write(iou,*) ""
        write(iou,"(a)") "      k    p_full [hPa]   dp [hPa]     sigma_full         A [Pa]        B"
        do k = 1, vg%nlev
            write(iou,"(i7,f15.3,f11.3,f15.5,f15.2,f10.5)") &
                    k, vg%p_ref_full(k)/100.0_wp, &
                    (vg%p_ref_half(k) - vg%p_ref_half(k-1))/100.0_wp, &
                    vg%sigma_full(k), vg%A(k), vg%B(k)
        end do
        write(iou,*) ""

        return

    end subroutine aeros_vgrid_print

end module aeros_vertical
