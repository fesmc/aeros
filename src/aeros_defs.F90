module aeros_defs
    ! Core definitions for aeros: precision, physical constants and the
    ! parameter / grid / state derived types.
    !
    ! This module contains no physics and no algorithms. See docs/design.md.

    use, intrinsic :: iso_fortran_env, only : error_unit
    use, intrinsic :: iso_c_binding,   only : c_double
    use precision, only : sp, dp
    use nml,       only : nml_read

    implicit none

    private

    ! === Precision ===========================================================
    !
    ! wp is selectable at COMPILE TIME: `make ... precision=dp` defines
    ! AEROS_DP. See config/Makefile for why both settings are real
    ! configurations rather than a production setting plus a validation toggle.
    !
    ! The short version: SHTns' interface is double precision, so an sp build
    ! pays a copy-convert on every spectral transform that a dp build does not.
    ! Against that, sp halves the memory traffic through the physics and the
    ! semi-implicit solve. M0a measures which wins; until then neither is
    ! privileged, and no code outside this module may assume `wp == sp`.
#ifdef AEROS_DP
    integer, parameter, public :: wp = dp
#else
    integer, parameter, public :: wp = sp
#endif

    ! Spectral coefficients are ALWAYS double precision, in both builds.
    !
    ! Two reasons, and the first alone is decisive: SHTns takes and returns
    ! complex(c_double), with no single-precision CPU path (SHT_FP32 is
    ! GPU-only), so a single-precision spectral array could not be passed
    ! across the interface at all. Second, spectral coefficients are the
    ! accumulators of the semi-implicit solve -- summed over every timestep of
    ! a 10^5 yr integration -- which is exactly where sp round-off compounds.
    !
    ! The cost is negligible: at T42L19 the full spectral state is under 1 MB.
    integer, parameter, public :: wp_sh = dp

    public :: sp, dp

    ! Fail the build if wp_sh ever stops matching C double, which would
    ! silently break the array passing in aeros_spectral. When the kinds match
    ! this declares a plain real(wp_sh); otherwise the kind is -1 (invalid) and
    ! the compiler rejects it. Borrowed from FastEarth3D's fe_precision.
    real(kind=merge(wp_sh, -1, wp_sh == c_double)), private :: enforce_wp_sh_eq_c_double

    ! === IO units and sentinels ==============================================

    integer,  parameter, public :: io_unit_err = error_unit
    real(wp), parameter, public :: MV          = -9999.0_wp

    ! === Mathematical constants ==============================================

    real(dp), parameter, public :: pi = 3.14159265358979323846_dp
    real(dp), parameter, public :: degrees_to_radians = pi/180.0_dp
    real(dp), parameter, public :: radians_to_degrees = 180.0_dp/pi

    ! === Physical constants ==================================================
    !
    ! Values follow ECMWF IFS Part IV (Physical Processes) Chapter 12, which is
    ! also what SPEEDY, PlaSim and SpeedyWeather use to within round-off.
    ! Declared dp regardless of wp: these are compile-time constants, so a
    ! double-precision value costs nothing and prevents a constant from
    ! becoming the least accurate number in an expression.

    real(dp), parameter, public :: omega   = 7.292e-5_dp     ! Earth angular velocity [rad s-1]
    real(dp), parameter, public :: r_earth = 6.371e6_dp      ! Earth mean radius [m]
    real(dp), parameter, public :: grav    = 9.80665_dp      ! gravitational acceleration [m s-2]

    real(dp), parameter, public :: R_gas   = 8.31446_dp      ! universal gas constant [J mol-1 K-1]
    real(dp), parameter, public :: R_d     = 287.05_dp       ! specific gas constant, dry air [J kg-1 K-1]
    real(dp), parameter, public :: R_v     = 461.51_dp       ! specific gas constant, water vapour [J kg-1 K-1]
    real(dp), parameter, public :: cp_d    = 1004.64_dp      ! specific heat at const. pressure, dry air [J kg-1 K-1]
    real(dp), parameter, public :: cp_v    = 1846.1_dp       ! specific heat at const. pressure, vapour [J kg-1 K-1]
    real(dp), parameter, public :: kappa   = R_d/cp_d        ! R_d/cp_d [-] (~0.2857)

    real(dp), parameter, public :: T0      = 273.15_dp       ! freezing point [K]
    real(dp), parameter, public :: p0      = 1.0e5_dp        ! reference pressure [Pa]

    real(dp), parameter, public :: L_v     = 2.5008e6_dp     ! latent heat of vaporization [J kg-1]
    real(dp), parameter, public :: L_f     = 3.337e5_dp      ! latent heat of fusion [J kg-1]
    real(dp), parameter, public :: L_s     = L_v + L_f       ! latent heat of sublimation [J kg-1]

    real(dp), parameter, public :: sigma_sb = 5.670374e-8_dp ! Stefan-Boltzmann [W m-2 K-4]
    real(dp), parameter, public :: S0       = 1361.0_dp      ! solar constant [W m-2]

    ! === Parameters ==========================================================

    type aeros_param_class
        ! Model configuration, read from a namelist group (default `aeros`).

        ! -- Spectral truncation and grid.
        !
        ! `trunc` is the triangular truncation T: lmax = mmax = trunc. Per
        ! docs/design.md section 10.3 this is THE open question of the design
        ! (T31 + resolution correction vs bare T42), which is why it is a
        ! namelist parameter from day one and never a compile-time constant.
        !
        ! nlon/nlat default to the smallest quadratically-unaliased Gaussian
        ! grid for `trunc` when left at -1 (see aeros_grid). Set them
        ! explicitly only to study aliasing or to match another model's grid.
        integer :: trunc
        integer :: nlon
        integer :: nlat

        ! -- Vertical levels. docs/design.md section 4.1: L16-20 with 3-4
        ! levels below 850 hPa. Fixed before any tuning, since level count
        ! moves the jet latitude and therefore the storm track.
        integer :: nlev

        ! -- Timestepping [s]. ~30 min at T42 (section 4).
        real(wp) :: dt

        ! -- Number of OpenMP threads for the spectral transforms. -1 leaves
        ! SHTns to use whatever OMP_NUM_THREADS says. Explicit values exist for
        ! the M0a scaling measurement (section 7), which sweeps 8/16/32/64.
        integer :: nthreads
    end type aeros_param_class

    ! === Grid ================================================================

    type aeros_grid_class
        ! The global Gaussian grid the spectral core lives on.
        !
        ! Latitude ordering is NORTH TO SOUTH, inherited from SHTns' colatitude
        ! convention (theta = 0 at the north pole). Do not reorder: the grid
        ! arrays index the same memory the transforms write.
        !
        ! Spatial fields are shaped (nlon, nlat) -- longitude contiguous. This
        ! matches SHTns' SHT_PHI_CONTIGUOUS layout, so no transpose is needed
        ! at the transform boundary, and it makes a latitude row contiguous,
        ! which is the loop the column physics parallelizes over.

        integer :: nlon                       ! longitudes
        integer :: nlat                       ! Gaussian latitudes
        integer :: ncol                       ! nlon*nlat, total columns

        real(wp), allocatable :: lon(:)       ! longitude [degrees east], (nlon)
        real(wp), allocatable :: lat(:)       ! latitude [degrees north], (nlat), descending
        real(dp), allocatable :: colat(:)     ! colatitude [rad], (nlat), ascending

        ! Gauss-Legendre quadrature weights, (nlat), summing to 2. Kept in dp
        ! in both builds: every global conservation check the design demands
        ! (section 3.7 risk 2, "verify to machine precision") runs through
        ! these, and an sp weight caps that check at ~1e-7.
        real(dp), allocatable :: gauss_w(:)

        real(wp), allocatable :: area(:,:)    ! cell area [m2], (nlon,nlat)
        real(wp), allocatable :: coriolis(:,:)! Coriolis parameter f [s-1], (nlon,nlat)
    end type aeros_grid_class

    ! === Prognostic state ====================================================

    type aeros_state_class
        ! The primitive-equation prognostic state (docs/design.md section 3.2),
        ! held in spectral space, with its grid-space counterpart alongside.
        !
        ! M0 allocates and zeroes these; nothing evolves them yet. The layout
        ! is declared now because it fixes the transform interface that M0a
        ! benchmarks and M1 fills in.
        !
        ! Spectral arrays are (nlm, nlev) -- coefficient index contiguous, so a
        ! single level is one contiguous SHTns call, and the level loop (which
        ! section 4.3 parallelizes) strides over whole transforms.

        integer :: nlm                        ! spectral coefficients per level
        integer :: nlev                       ! vertical levels

        ! -- Spectral prognostics. Vorticity/divergence rather than u,v: it is
        ! what makes the semi-implicit gravity-wave solve diagonal in l.
        complex(wp_sh), allocatable :: vor(:,:)   ! relative vorticity [s-1]
        complex(wp_sh), allocatable :: div(:,:)   ! divergence [s-1]
        complex(wp_sh), allocatable :: temp(:,:)  ! temperature [K]
        complex(wp_sh), allocatable :: qv(:,:)    ! specific humidity [kg kg-1]
        complex(wp_sh), allocatable :: lnps(:)    ! ln(surface pressure) [-], (nlm)

        ! -- Grid-space diagnostics, recomputed from the spectral state each
        ! step and consumed by the column physics.
        real(wp), allocatable :: u(:,:,:)         ! zonal wind [m s-1], (nlon,nlat,nlev)
        real(wp), allocatable :: v(:,:,:)         ! meridional wind [m s-1]
        real(wp), allocatable :: temp_g(:,:,:)    ! temperature [K]
        real(wp), allocatable :: qv_g(:,:,:)      ! specific humidity [kg kg-1]
        real(wp), allocatable :: ps(:,:)          ! surface pressure [Pa], (nlon,nlat)
    end type aeros_state_class

    public :: aeros_param_class
    public :: aeros_grid_class
    public :: aeros_state_class

    public :: aeros_par_load

contains

    subroutine aeros_par_load(par, filename, group)
        ! Read the model parameters from a Fortran namelist.

        implicit none

        type(aeros_param_class), intent(out) :: par
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: group

        character(len=56) :: nml_group

        nml_group = "aeros"
        if (present(group)) nml_group = trim(group)

        call nml_read(filename, nml_group, "trunc",    par%trunc)
        call nml_read(filename, nml_group, "nlon",     par%nlon)
        call nml_read(filename, nml_group, "nlat",     par%nlat)
        call nml_read(filename, nml_group, "nlev",     par%nlev)
        call nml_read(filename, nml_group, "dt",       par%dt)
        call nml_read(filename, nml_group, "nthreads", par%nthreads)

        if (par%trunc < 1) then
            write(io_unit_err,*) "aeros_par_load:: error: trunc must be >= 1, got ", par%trunc
            stop 1
        end if

        if (par%nlev < 1) then
            write(io_unit_err,*) "aeros_par_load:: error: nlev must be >= 1, got ", par%nlev
            stop 1
        end if

        return

    end subroutine aeros_par_load

end module aeros_defs
