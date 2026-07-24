module aeros_spectral
    ! Fortran wrapper around SHTns, the spherical-harmonic transform library.
    !
    ! This is the performance kernel of the model. docs/design.md section 4
    ! estimates the Legendre transforms alone at 25-40% of runtime at T42L19,
    ! and section 3.2 chooses a spectral core precisely because triangular
    ! truncation is isotropic on the sphere -- no polar filter, no polar CFL
    ! penalty, effective resolution at 85 degrees N identical to the equator.
    ! For a polar-priority model that property is structural, and it lives here.
    !
    ! fesm-utils builds SHTns but supplies no Fortran module for it (only the
    ! bare ISO_C_BINDING header, `shtns.f03`), so this wrapper is aeros' own
    ! (docs/design.md section 8). It exists to keep iso_c_binding, SHTns'
    ! 0-based coefficient indexing and its double-precision interface out of
    ! the rest of the model.
    !
    ! CONVENTIONS, fixed here and assumed everywhere downstream:
    !
    !   normalization  orthonormal, no Condon-Shortley phase
    !                  (SHT_ORTHONORMAL + SHT_NO_CS_PHASE), matching
    !                  FastEarth3D so the two can be compared directly.
    !   layout         SHT_PHI_CONTIGUOUS: spatial fields are (nlon, nlat),
    !                  longitude contiguous, latitude running north to south.
    !   truncation     triangular, lmax = mmax = trunc, mres = 1.
    !   spectral kind  always double (wp_sh), in both sp and dp builds --
    !                  SHTns has no single-precision CPU path.
    !
    ! NOT YET HERE, deliberately: the vorticity/divergence <-> (u,v) mapping.
    ! The raw spheroidal/toroidal transforms below are its building blocks, but
    ! the factors of a and of -l(l+1) that connect them to the primitive
    ! equations belong with the dynamical core (M1), where they can be
    ! validated against Held-Suarez rather than asserted here.

    use, intrinsic :: iso_c_binding
    use aeros_defs, only : sp, dp, wp, wp_sh, pi, r_earth, io_unit_err

    implicit none

    private

    ! SHTns' Fortran 2003 interface, installed into the SHTns prefix include
    ! dir by fesm-utils' build and found via -I$(SHTNSROOT)/include. It is a
    ! bare include, not a module: it declares parameters, the shtns_info bind(C)
    ! type and one large interface block.
    include 'shtns.f03'

    type aeros_sht_class
        ! One configured transform: a spectral truncation plus its grid.
        !
        ! A single SHTns config is NOT safe for concurrent transform calls. The
        ! threading inside SHTns (the omp build) handles one transform on many
        ! threads, which is the parallelism docs/design.md section 4.3 assumes.
        ! If a level loop is ever threaded around whole transforms instead, each
        ! thread needs its own config -- see shtns_create/shtns_set_grid, and
        ! note that creation is not thread-safe (FFTW planning).

        type(c_ptr) :: cfg = c_null_ptr    ! opaque SHTns configuration handle

        integer :: trunc = 0               ! triangular truncation T
        integer :: lmax  = 0               ! max degree  (= trunc)
        integer :: mmax  = 0               ! max order   (= trunc)
        integer :: mres  = 1               ! order step
        integer :: nlon  = 0               ! longitudes (SHTns nphi)
        integer :: nlat  = 0               ! Gaussian latitudes
        integer :: nlm   = 0               ! spectral coefficients, m >= 0

        ! Degree and order of each coefficient, (nlm), copied from SHTns at
        ! init so that callers never dereference a C pointer. Used by every
        ! spectral operator: Laplacian, diffusion, truncation, energy spectra.
        integer, allocatable :: l_of_lm(:)
        integer, allocatable :: m_of_lm(:)

        ! Grid geometry, read back from SHTns rather than recomputed, so there
        ! is exactly one source of truth for the Gauss nodes.
        !
        ! `sinlat` is SHTns' own cos(theta) array verbatim -- the Gauss nodes
        ! themselves, and the ONLY exactly-antisymmetric representation of
        ! them. Anything derived through acos and back (a colatitude, then a
        ! cosine of it) is symmetric only to a rounding error, which is enough
        ! to break a hemispheric-symmetry check at double precision. Prefer
        ! sinlat wherever the quantity wanted is sin(lat) or cos(colat).
        real(dp), allocatable :: sinlat(:)  ! sin(latitude) = cos(colatitude), (nlat)
        real(dp), allocatable :: colat(:)   ! colatitude [rad], (nlat), ascending
        real(dp), allocatable :: lon(:)     ! longitude [rad], (nlon)
        real(dp), allocatable :: gauss_w(:) ! quadrature weights, (nlat), sum = 2
    end type aeros_sht_class

    ! Transforms are generic over the grid-side kind so that both an sp and a
    ! dp build call the same name. The sp specifics convert through a
    ! double-precision automatic array; the dp specifics pass straight through
    ! with no copy. Automatic (not module) buffers, so the routines stay
    ! callable from inside a parallel region.
    !
    ! CONTRACT, uniform across both kinds: analysis MAY DESTROY its input
    ! field, hence intent(inout). SHTns does this for real in the dp path; the
    ! sp path does not, but relying on that would make the same call behave
    ! differently in the two builds.

    interface aeros_sht_analysis
        module procedure aeros_sht_analysis_sp
        module procedure aeros_sht_analysis_dp
    end interface aeros_sht_analysis

    interface aeros_sht_synthesis
        module procedure aeros_sht_synthesis_sp
        module procedure aeros_sht_synthesis_dp
    end interface aeros_sht_synthesis

    interface aeros_sht_analysis_vec
        module procedure aeros_sht_analysis_vec_sp
        module procedure aeros_sht_analysis_vec_dp
    end interface aeros_sht_analysis_vec

    interface aeros_sht_synthesis_vec
        module procedure aeros_sht_synthesis_vec_sp
        module procedure aeros_sht_synthesis_vec_dp
    end interface aeros_sht_synthesis_vec

    interface aeros_sht_surface_integral
        module procedure aeros_sht_surface_integral_sp
        module procedure aeros_sht_surface_integral_dp
    end interface aeros_sht_surface_integral

    public :: aeros_sht_class
    public :: aeros_sht_init
    public :: aeros_sht_end
    public :: aeros_sht_grid_size
    public :: aeros_sht_analysis
    public :: aeros_sht_synthesis
    public :: aeros_sht_analysis_vec
    public :: aeros_sht_synthesis_vec
    public :: aeros_sht_laplacian
    public :: aeros_sht_surface_integral
    public :: aeros_sht_lm

contains

    subroutine aeros_sht_init(sht, trunc, nlon, nlat, nthreads, quick)
        ! Create the SHTns configuration and its Gaussian grid.

        implicit none

        type(aeros_sht_class), intent(inout) :: sht
        integer, intent(in) :: trunc                 ! triangular truncation T
        integer, intent(in), optional :: nlon, nlat  ! grid size override
        integer, intent(in), optional :: nthreads    ! -1 or absent: leave to OMP_NUM_THREADS
        logical, intent(in), optional :: quick       ! skip SHTns' algorithm auto-tuning

        type(shtns_info), pointer :: info
        integer(c_short), pointer :: li(:), mi(:)
        real(c_double),   pointer :: cos_theta(:)
        real(c_double), allocatable :: wts_half(:)
        real(c_double) :: eps
        integer :: nlon_, nlat_, norm, layout, nh, i, ierr
        logical :: quick_

        quick_ = .FALSE.
        if (present(quick)) quick_ = quick

        ! Grid: the smallest quadratically-unaliased Gaussian grid for this
        ! truncation, unless the caller pins it.
        call aeros_sht_grid_size(trunc, nlon_, nlat_)
        if (present(nlon)) then
            if (nlon > 0) nlon_ = nlon
        end if
        if (present(nlat)) then
            if (nlat > 0) nlat_ = nlat
        end if

        sht%trunc = trunc
        sht%lmax  = trunc
        sht%mmax  = trunc
        sht%mres  = 1

        ! Thread count must be set before any config is created.
        if (present(nthreads)) then
            if (nthreads > 0) ierr = shtns_use_threads(nthreads)
        end if

        ! eps is the polar-optimization threshold: SHTns skips Legendre
        ! contributions below it near the poles. 1e-10 is SHTns' own
        ! recommendation for accuracy-critical use; 0 disables the
        ! optimization entirely. Do not raise it without re-running the
        ! round-trip test in tests/test_spectral.f90 -- polar accuracy is the
        ! whole reason for the spectral core (docs/design.md section 3.2).
        eps  = 1.0e-10_c_double
        norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

        ! SHT_GAUSS auto-tunes the transform algorithm at init, which costs
        ! seconds but is amortized instantly over a paleo integration.
        ! SHT_QUICK_INIT skips it (transforms then run ~20% slower) and exists
        ! for tests and short benchmarks that would otherwise be dominated by
        ! setup. If init cost ever matters in production, add
        ! SHT_LOAD_SAVE_CFG to cache the tuning to disk.
        layout = SHT_PHI_CONTIGUOUS
        if (quick_) then
            layout = layout + SHT_QUICK_INIT
        else
            layout = layout + SHT_GAUSS
        end if

        sht%cfg = shtns_create(sht%lmax, sht%mmax, sht%mres, norm)
        ! NB argument order: SHTns takes (nlat, nphi), latitude first.
        call shtns_set_grid(sht%cfg, layout, eps, nlat_, nlon_)

        ! Read the realized configuration back out of SHTns rather than
        ! assuming it honoured the request.
        call c_f_pointer(sht%cfg, info)
        sht%nlat = info%nlat
        sht%nlon = info%nphi
        sht%nlm  = info%nlm

        if (info%nspat /= sht%nlon*sht%nlat) then
            ! Would mean SHTns padded the spatial array, which breaks the
            ! assumption that a (nlon,nlat) Fortran array can be passed
            ! straight through. It cannot happen without SHT_ALLOW_PADDING,
            ! but the cost of checking is one comparison at init.
            write(io_unit_err,*) "aeros_sht_init:: error: padded spatial layout, nspat = ", &
                                    info%nspat, " /= nlon*nlat = ", sht%nlon*sht%nlat
            stop 1
        end if

        ! Degree/order of each coefficient (SHTns stores them as unsigned
        ! short; every value is <= lmax, so a signed short holds them).
        call c_f_pointer(info%li, li, [sht%nlm])
        call c_f_pointer(info%mi, mi, [sht%nlm])
        allocate(sht%l_of_lm(sht%nlm))
        allocate(sht%m_of_lm(sht%nlm))
        sht%l_of_lm = int(li)
        sht%m_of_lm = int(mi)

        ! Gauss nodes, from SHTns' own cos(theta) array.
        call c_f_pointer(info%ct, cos_theta, [sht%nlat])
        allocate(sht%sinlat(sht%nlat))
        allocate(sht%colat(sht%nlat))
        sht%sinlat = cos_theta
        sht%colat  = acos(cos_theta)

        ! Uniform longitudes.
        allocate(sht%lon(sht%nlon))
        do i = 1, sht%nlon
            sht%lon(i) = 2.0_dp*pi*real(i-1, dp)/real(sht%nlon, dp)
        end do

        ! Gauss weights: SHTns returns one hemisphere; mirror to the full grid.
        nh = info%nlat_2
        allocate(wts_half(nh))
        call shtns_gauss_wts(sht%cfg, wts_half)
        allocate(sht%gauss_w(sht%nlat))
        sht%gauss_w(1:nh) = wts_half
        sht%gauss_w(sht%nlat:sht%nlat-nh+1:-1) = wts_half
        deallocate(wts_half)

        return

    end subroutine aeros_sht_init

    subroutine aeros_sht_end(sht)
        ! Release the SHTns configuration and the cached geometry.

        implicit none

        type(aeros_sht_class), intent(inout) :: sht

        if (c_associated(sht%cfg)) call shtns_destroy(sht%cfg)
        sht%cfg = c_null_ptr

        sht%trunc = 0; sht%lmax = 0; sht%mmax = 0
        sht%nlon  = 0; sht%nlat = 0; sht%nlm  = 0

        if (allocated(sht%l_of_lm)) deallocate(sht%l_of_lm)
        if (allocated(sht%m_of_lm)) deallocate(sht%m_of_lm)
        if (allocated(sht%sinlat))  deallocate(sht%sinlat)
        if (allocated(sht%colat))   deallocate(sht%colat)
        if (allocated(sht%lon))     deallocate(sht%lon)
        if (allocated(sht%gauss_w)) deallocate(sht%gauss_w)

        return

    end subroutine aeros_sht_end

    subroutine aeros_sht_grid_size(trunc, nlon, nlat)
        ! Smallest quadratically-unaliased Gaussian grid for truncation T.
        !
        ! The quadratic ("linear-in-name-only") rule nlon >= 3T+1 is what makes
        ! the transform of a product of two truncated fields exact -- i.e. it
        ! is what keeps advection free of aliasing. T31 -> 96x48, T42 -> 128x64,
        ! T85 -> 256x128, the standard grids in the literature aeros is
        ! calibrated against (docs/design.md section 3.1).
        !
        ! nlon is rounded up to a 2/3/5-smooth even number so FFTW gets a
        ! length it has a fast codelet for.
        !
        ! This is the FULL Gaussian grid. docs/design.md section 4 wants a
        ! reduced/octahedral grid (~30% fewer gridpoints), which SHTns does not
        ! support: its longitudinal FFT assumes a constant nlon per latitude.
        ! Reducing the grid therefore means carrying a second, reduced grid for
        ! the column physics and interpolating to it -- a real saving on the
        ! physics but not on the transforms, and not free. Left for after M0a
        ! decides how much of the budget the transforms actually take.

        implicit none

        integer, intent(in)  :: trunc
        integer, intent(out) :: nlon, nlat

        integer :: n

        n = 3*trunc + 1
        do while (.not. is_smooth_even(n))
            n = n + 1
        end do

        nlon = n
        nlat = n/2

        return

    end subroutine aeros_sht_grid_size

    logical function is_smooth_even(n) result(ok)
        ! True if n is even and factors entirely into 2, 3 and 5.

        implicit none

        integer, intent(in) :: n
        integer :: m

        ok = .FALSE.
        if (mod(n,2) /= 0) return

        m = n
        do while (mod(m,2) == 0); m = m/2; end do
        do while (mod(m,3) == 0); m = m/3; end do
        do while (mod(m,5) == 0); m = m/5; end do

        ok = (m == 1)

        return

    end function is_smooth_even

    integer function aeros_sht_lm(sht, l, m) result(lm)
        ! 1-based index of coefficient (l,m) in a spectral array.
        !
        ! Mirrors SHTns' LM() macro, shifted by one: SHTns indexes from 0, the
        ! rest of aeros from 1. Returns -1 for an (l,m) outside the truncation
        ! rather than an out-of-range index.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        integer, intent(in) :: l, m

        if (m < 0 .or. m > sht%mmax*sht%mres .or. l < m .or. l > sht%lmax &
                .or. mod(m, sht%mres) /= 0) then
            lm = -1
            return
        end if

        lm = ((m/sht%mres)*(2*sht%lmax + 2 - (m + sht%mres)))/2 + l + 1

        return

    end function aeros_sht_lm

    ! === Scalar transforms ===================================================

    subroutine aeros_sht_analysis_dp(sht, field, coeffs)
        ! Spatial field (nlon,nlat) -> spectral coefficients (nlm).
        ! NB: SHTns overwrites `field`.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(dp), intent(inout), contiguous :: field(:,:)
        complex(wp_sh), intent(out), contiguous :: coeffs(:)

        call spat_to_SH(sht%cfg, field, coeffs)

        return

    end subroutine aeros_sht_analysis_dp

    subroutine aeros_sht_analysis_sp(sht, field, coeffs)
        ! Single-precision grid side: convert up, transform, discard.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(sp), intent(inout), contiguous :: field(:,:)
        complex(wp_sh), intent(out), contiguous :: coeffs(:)

        real(dp) :: buf(sht%nlon, sht%nlat)

        buf = real(field, dp)
        call spat_to_SH(sht%cfg, buf, coeffs)

        return

    end subroutine aeros_sht_analysis_sp

    subroutine aeros_sht_synthesis_dp(sht, coeffs, field)
        ! Spectral coefficients (nlm) -> spatial field (nlon,nlat).

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        complex(wp_sh), intent(in), contiguous :: coeffs(:)
        real(dp), intent(out), contiguous :: field(:,:)

        call SH_to_spat(sht%cfg, coeffs, field)

        return

    end subroutine aeros_sht_synthesis_dp

    subroutine aeros_sht_synthesis_sp(sht, coeffs, field)

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        complex(wp_sh), intent(in), contiguous :: coeffs(:)
        real(sp), intent(out), contiguous :: field(:,:)

        real(dp) :: buf(sht%nlon, sht%nlat)

        call SH_to_spat(sht%cfg, coeffs, buf)
        field = real(buf, sp)

        return

    end subroutine aeros_sht_synthesis_sp

    ! === Vector transforms ===================================================
    !
    ! Helmholtz decomposition of a horizontal field into a spheroidal
    ! (divergent) and a toroidal (rotational) potential:
    !
    !     (v_theta, v_phi) = grad_h(S) + curl_h(T)
    !
    ! with v_theta positive SOUTHWARD (colatitude increases southward), which
    ! is SHTns' convention, not the meteorological one. Converting to (u,v) and
    ! relating S and T to divergence and vorticity is the dynamical core's job
    ! (M1) -- see the module header.

    subroutine aeros_sht_analysis_vec_dp(sht, vth, vph, slm, tlm)
        ! (v_theta, v_phi) -> spheroidal + toroidal coefficients.
        ! NB: SHTns overwrites both input fields.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(dp), intent(inout), contiguous :: vth(:,:), vph(:,:)
        complex(wp_sh), intent(out), contiguous :: slm(:), tlm(:)

        call spat_to_SHsphtor(sht%cfg, vth, vph, slm, tlm)

        return

    end subroutine aeros_sht_analysis_vec_dp

    subroutine aeros_sht_analysis_vec_sp(sht, vth, vph, slm, tlm)

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(sp), intent(inout), contiguous :: vth(:,:), vph(:,:)
        complex(wp_sh), intent(out), contiguous :: slm(:), tlm(:)

        real(dp) :: bth(sht%nlon, sht%nlat)
        real(dp) :: bph(sht%nlon, sht%nlat)

        bth = real(vth, dp)
        bph = real(vph, dp)
        call spat_to_SHsphtor(sht%cfg, bth, bph, slm, tlm)

        return

    end subroutine aeros_sht_analysis_vec_sp

    subroutine aeros_sht_synthesis_vec_dp(sht, slm, tlm, vth, vph)
        ! Spheroidal + toroidal coefficients -> (v_theta, v_phi).

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        complex(wp_sh), intent(in), contiguous :: slm(:), tlm(:)
        real(dp), intent(out), contiguous :: vth(:,:), vph(:,:)

        call SHsphtor_to_spat(sht%cfg, slm, tlm, vth, vph)

        return

    end subroutine aeros_sht_synthesis_vec_dp

    subroutine aeros_sht_synthesis_vec_sp(sht, slm, tlm, vth, vph)

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        complex(wp_sh), intent(in), contiguous :: slm(:), tlm(:)
        real(sp), intent(out), contiguous :: vth(:,:), vph(:,:)

        real(dp) :: bth(sht%nlon, sht%nlat)
        real(dp) :: bph(sht%nlon, sht%nlat)

        call SHsphtor_to_spat(sht%cfg, slm, tlm, bth, bph)
        vth = real(bth, sp)
        vph = real(bph, sp)

        return

    end subroutine aeros_sht_synthesis_vec_sp

    ! === Spectral operators ==================================================

    subroutine aeros_sht_laplacian(sht, coeffs, inverse)
        ! Horizontal Laplacian on the sphere, in place.
        !
        !     lap(Y_lm) = -l(l+1)/a^2 * Y_lm
        !
        ! Diagonal in spectral space -- the reason the semi-implicit
        ! gravity-wave solve decouples by wavenumber (docs/design.md section
        ! 3.2). With inverse=.TRUE. the l=0 mode is set to zero, which is the
        ! only choice that makes the inverse well-posed (the global mean of a
        ! Laplacian is zero, so its inverse is defined up to a constant).

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        complex(wp_sh), intent(inout) :: coeffs(:)
        logical, intent(in), optional :: inverse

        real(dp) :: fac, a2
        integer  :: lm, l
        logical  :: inv

        inv = .FALSE.
        if (present(inverse)) inv = inverse

        a2 = r_earth*r_earth

        do lm = 1, sht%nlm
            l = sht%l_of_lm(lm)
            if (inv) then
                if (l == 0) then
                    coeffs(lm) = (0.0_wp_sh, 0.0_wp_sh)
                    cycle
                end if
                fac = -a2/real(l*(l+1), dp)
            else
                fac = -real(l*(l+1), dp)/a2
            end if
            coeffs(lm) = coeffs(lm)*fac
        end do

        return

    end subroutine aeros_sht_laplacian

    real(dp) function aeros_sht_surface_integral_dp(sht, field) result(total)
        ! Global integral of a spatial field over the sphere, by Gauss-Legendre
        ! quadrature in latitude and the uniform rule in longitude. Returns
        ! 4*pi*a^2 for field = 1.
        !
        ! Always accumulated in dp, in both builds. Every conservation check
        ! the design demands "to machine precision" (docs/design.md section 3.7
        ! risk 2, section 7 M4) runs through here, and an sp accumulator over
        ! nlon*nlat terms would cap that check several orders of magnitude
        ! above machine precision.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(dp), intent(in) :: field(:,:)

        real(dp) :: dlon, row
        integer  :: i, j

        dlon  = 2.0_dp*pi/real(sht%nlon, dp)
        total = 0.0_dp

        do j = 1, sht%nlat
            row = 0.0_dp
            do i = 1, sht%nlon
                row = row + field(i,j)
            end do
            total = total + sht%gauss_w(j)*row*dlon
        end do

        total = total*r_earth*r_earth

        return

    end function aeros_sht_surface_integral_dp

    real(dp) function aeros_sht_surface_integral_sp(sht, field) result(total)

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(sp), intent(in) :: field(:,:)

        real(dp) :: dlon, row
        integer  :: i, j

        dlon  = 2.0_dp*pi/real(sht%nlon, dp)
        total = 0.0_dp

        do j = 1, sht%nlat
            row = 0.0_dp
            do i = 1, sht%nlon
                row = row + real(field(i,j), dp)
            end do
            total = total + sht%gauss_w(j)*row*dlon
        end do

        total = total*r_earth*r_earth

        return

    end function aeros_sht_surface_integral_sp

end module aeros_spectral
