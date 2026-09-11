module aeros_bundles
    ! Coupling-cycle exchange bundles (docs/coupling_api.md sec. 6).
    !
    ! Three types:
    !   aeros_imports_t -- what a coupled host pushes into aeros each cycle
    !                      (SST, albedo, so4, ...; year-boundary scalars CO2,
    !                      time_bp, ...). Shape (nlon,nlat) on the aeros native
    !                      grid. Fields declared but not allocated by default
    !                      -- the caller allocates only what it wants to
    !                      supply; aeros_lib_step tests `allocated()` before
    !                      copying into the corresponding BC_EXTERNAL buffer.
    !
    !   aeros_exports_t -- cycle-mean grid-mean fields aeros returns each cycle
    !                      (surface fluxes, near-surface state, radiation at
    !                      surface, precipitation, SESAM-style column diags).
    !                      Same allocation convention: the host allocates the
    !                      fields it wants.
    !
    !   aeros_diag_t    -- optional detailed diagnostic snapshot (sec. 6.3).
    !                      Phase A ships an empty shell (name-reserved);
    !                      Phase C fills it with the rce_long accumulator
    !                      contents.
    !
    ! Wave-3.a scope: types declared, aeros_lib_step's copy logic is a shell.
    ! Phase A wires the subset with aeros-side readers today:
    !   imports: sst, albedo (single-band), so4, co2, time_bp
    !   exports: ps, slp, t2m, q2m, u10m, v10m, taux, tauy, sw_dw_sur,
    !            lw_dw_sur, rain, snow, tam, ram, gams, cld + TOA trio
    ! The unwired fields are shaped in the type but ignored -- Phase C wires
    ! them.
    !
    ! Leaf module (no aeros-internal dependencies beyond aeros_defs) so the
    ! bundle types can be used by both aeros and a coupled adapter without
    ! either pulling in the other.

    use aeros_defs, only : wp

    implicit none

    private

    ! === Imports: what the coupler pushes in each cycle ======================
    type, public :: aeros_imports_t
        ! Cycle-frequency 2D fields (all shape (nlon,nlat)).
        real(wp), allocatable :: sst(:,:)          ! sea surface temperature [K]
        real(wp), allocatable :: sic(:,:)          ! sea-ice fraction [0-1]
        real(wp), allocatable :: t_ice_sfc(:,:)    ! sea-ice surface T [K]
        real(wp), allocatable :: t_skin_land(:,:)  ! land skin T [K]
        real(wp), allocatable :: q_skin_land(:,:)  ! land skin q [kg/kg]
        real(wp), allocatable :: albedo(:,:)       ! broadband albedo [-]
        real(wp), allocatable :: z0m(:,:)          ! momentum roughness [m]
        real(wp), allocatable :: sh(:,:)           ! sensible heat flux [W m-2]
        real(wp), allocatable :: lh(:,:)           ! latent   heat flux [W m-2]
        real(wp), allocatable :: evp(:,:)          ! evaporation [kg m-2 s-1]
        real(wp), allocatable :: lwu(:,:)          ! upward LW from surface [W m-2]
        real(wp), allocatable :: so4_load(:,:)     ! sulfate load [kg m-2]
        real(wp), allocatable :: dust(:,:)         ! dust optical thickness [-]
        real(wp), allocatable :: lsm(:,:)          ! land-sea mask / land fraction [0-1]
        real(wp), allocatable :: phis(:,:)         ! surface geopotential [m2 s-2] (raw; band-limited by aeros)

        ! Year-boundary scalar fields. Re-read only when year_boundary = T.
        logical  :: year_boundary = .false.
        real(wp) :: co2_ppm  = 0.0_wp
        real(wp) :: ch4_ppb  = 0.0_wp
        real(wp) :: n2o_ppb  = 0.0_wp
        real(wp) :: time_bp  = 0.0_wp
    end type aeros_imports_t

    ! === Exports: cycle-mean fields the coupler reads out each cycle =========
    type, public :: aeros_exports_t
        ! Surface fluxes / near-surface state.
        real(wp), allocatable :: ps(:,:)           ! surface pressure [Pa]
        real(wp), allocatable :: slp(:,:)          ! sea-level pressure [Pa]
        real(wp), allocatable :: t2m(:,:)          ! 2 m T [K]
        real(wp), allocatable :: q2m(:,:)          ! 2 m q [kg/kg]
        real(wp), allocatable :: u10m(:,:)         ! 10 m u [m/s]
        real(wp), allocatable :: v10m(:,:)         ! 10 m v [m/s]
        real(wp), allocatable :: wind10m(:,:)      ! 10 m wind speed [m/s]
        real(wp), allocatable :: taux(:,:)         ! zonal wind stress [N m-2]
        real(wp), allocatable :: tauy(:,:)         ! meridional wind stress [N m-2]

        ! Radiation at surface.
        real(wp), allocatable :: sw_dw_sur(:,:)    ! downward SW [W m-2]
        real(wp), allocatable :: lw_dw_sur(:,:)    ! downward LW [W m-2]
        real(wp), allocatable :: sw_net_sur(:,:)   ! net SW [W m-2]

        ! Water.
        real(wp), allocatable :: rain(:,:)         ! liquid precip [kg m-2 s-1]
        real(wp), allocatable :: snow(:,:)         ! solid  precip [kg m-2 s-1]
        real(wp), allocatable :: prc_total(:,:)    ! rain + snow  [kg m-2 s-1]
        ! Convective vs large-scale precip breakdown; sum equals prc_total.
        real(wp), allocatable :: prc_cnv(:,:)      ! SBM convective precip [kg m-2 s-1]
        real(wp), allocatable :: prc_cnd(:,:)      ! large-scale precip    [kg m-2 s-1]

        ! Grid-mean climate diagnostics (SESAM-style).
        real(wp), allocatable :: tam(:,:)          ! column-mean T [K]
        real(wp), allocatable :: ram(:,:)          ! column-mean RH [0-1]
        real(wp), allocatable :: gams(:,:)         ! near-surface lapse rate [K m-1]
        real(wp), allocatable :: u700(:,:)         ! 700 hPa u [m/s]
        real(wp), allocatable :: v700(:,:)         ! 700 hPa v [m/s]
        real(wp), allocatable :: cld(:,:)          ! total cloud cover [0-1]

        ! Cycle-mean TOA (nice-to-have, Phase A).
        real(wp), allocatable :: olr_toa(:,:)      ! [W m-2]
        real(wp), allocatable :: sw_toa_net(:,:)   ! [W m-2]
        real(wp), allocatable :: sw_toa_dw(:,:)    ! [W m-2]
    end type aeros_exports_t

    ! === Detailed diagnostic snapshot (empty shell in Phase A) ===============
    !
    ! docs/coupling_api.md sec. 6.3: "Phase A ships an empty shell (name-
    ! reserved). Phase C fills it with what today's rce_long accumulators
    ! produce."
    type, public :: aeros_diag_t
        integer :: reserved = 0  ! placeholder to keep the type non-empty
    end type aeros_diag_t

end module aeros_bundles
