module aeros_ecckd
    ! ecCKD-style correlated-k radiation for aeros -- the eventual endpoint named
    ! in design.md section 5 (option 1), reached via the `scheme` selector in
    ! aeros_radiation (SCHEME_ECCKD). SESAM (SCHEME_SESAM) stays the default and
    ! is untouched; this module is the opt-in path.
    !
    ! === What this is, and the gas-optics-data decision =====================
    !
    ! A "few-g-point" correlated-k scheme: a handful of spectral bands, each with
    ! a small k-distribution (a few g-points) per absorbing gas (H2O, CO2, O3).
    ! The flux transfer per g-point is the same absorptivity-emissivity integral
    ! the SESAM port uses (no scattering in the clear-sky longwave), but the
    ! transmission is now a genuine exponential in the g-point's absorption
    ! coefficient, T_g = exp(-k_g u), evaluated per gas and summed over g-points
    ! weighted by the band's Planck fraction. That is the structural fix over
    ! SESAM's single broadband fit: a real window band that stays transparent, a
    ! real 15um CO2 band that saturates logarithmically (so CO2 doubling comes out
    ! near the canonical ~3.7 W/m2, not SESAM's ~3.0), and gas overlap resolved
    ! band-by-band rather than as a product of three broadband fits.
    !
    ! GAS-OPTICS DATA (design.md section 5, hard-constraint 2). There is no ecCKD
    ! look-up table on this machine and the ecmwf/ecckd generator + CKDMIP
    ! line-by-line archive are not available offline (verified). So the k-tables
    ! are NOT a bundled binary LUT (option 2a); they are built from a compact,
    ! published statistical band model (option 2b) -- a Malkmus/Goody random-band
    ! model whose few band parameters are hard-coded, from which the per-g-point
    ! k-distribution follows analytically by Gaussian quadrature in g-space. No
    ! run-time data dependency, fully reproducible, and physically grounded.
    !
    ! === Current state: PHASE 0 (scaffolding) ==============================
    !
    ! This first commit wires the seam and stands up the correlated-k *machinery*
    ! (the band / g-point loop and the per-g-point emissivity-method flux solver)
    ! with a deliberately minimal one-band, one-g-point PLACEHOLDER table so the
    ! path compiles, runs, and returns physical, SESAM-ballpark clear-sky longwave
    ! fluxes end to end. It is a grey band dressed in the g-point structure, NOT
    ! the delivered scheme: the multi-band Malkmus-derived k-tables and Planck
    ! band weights replace the placeholder table in Phase 1, with zero change to
    ! the solver below. Everything marked PHASE-0 is the throwaway placeholder.
    !
    ! Not yet here (Phase 1+): the real multi-band k-tables; temperature/pressure
    ! scaling of k; the shortwave sibling; the all-sky (grey-cloud) branch. Until
    ! then aeros_radiation routes the ecCKD shortwave and all cloudy columns
    ! through the SESAM kernels, so grey clouds keep working (hard-constraint).

    use aeros_defs, only : wp, sigma_sb, grav, cp_d, R_d, T0

    implicit none

    private

    ! === Correlated-k table structure ======================================
    ! A band carries a Planck weight (its fraction of the blackbody flux; for a
    ! true multi-band table this is a function of temperature, tabulated -- here,
    ! Phase 0, a single band spanning the whole spectrum, weight 1) and a small
    ! set of g-points. Each g-point has a quadrature weight (summing to 1 over the
    ! band) and an absorption coefficient per gas [cm2 g-1] acting on that gas's
    ! mass path [g cm-2].
    integer, parameter :: NG_MAX = 8       ! max g-points per band (table cap)

    type :: ecckd_band
        real(wp) :: planck_wt = 1.0_wp     ! Planck fraction in this band [-]
        integer  :: ng        = 1          ! number of g-points
        real(wp) :: gw(NG_MAX)  = 0.0_wp   ! g-point weights (sum to 1)
        real(wp) :: kh2o(NG_MAX)= 0.0_wp   ! H2O absorption coeff [cm2 g-1]
        real(wp) :: kco2(NG_MAX)= 0.0_wp   ! CO2 absorption coeff [cm2 g-1]
        real(wp) :: ko3(NG_MAX) = 0.0_wp   ! O3  absorption coeff [cm2 g-1]
    end type ecckd_band

    ! === PHASE-0 PLACEHOLDER longwave table ================================
    ! One band, one g-point: a grey broadband absorber. The coefficients are
    ! round, hand-set values chosen only to land the clear-sky OLR in the
    ! SESAM/ERA5 ballpark for a midlatitude column -- they are NOT a spectroscopic
    ! fit and carry no provenance. They exist so the machinery runs end to end and
    ! are replaced wholesale by the Malkmus-derived multi-band table in Phase 1.
    !
    !   kh2o : broadband LW water-vapour mass absorption, flux-diffusivity folded
    !          out (the 1.66 is applied in the solver). ~0.55 cm2 g-1 gives a
    !          column emissivity ~0.8 at a few g cm-2 of vapour.
    !   kco2 : a token CO2 absorption so the seam demonstrably responds to the CO2
    !          mixing ratio (a single exponential cannot reproduce the log
    !          saturation of the real 15um band -- that is Phase 1's job).
    !   ko3  : token ozone absorption.
    integer, parameter :: N_LW_BANDS = 1

    real(wp), parameter :: LW_BETA0 = 1.66_wp        ! diffusivity factor

    public :: aeros_ecckd_lw_clearsky_column

contains

    subroutine ecckd_lw_table(band)
        ! Fill the Phase-0 placeholder longwave table. Phase 1 replaces this with
        ! the Malkmus-derived multi-band k-distribution (and its T-dependent
        ! Planck weights).
        type(ecckd_band), intent(out) :: band(N_LW_BANDS)

        ! --- PHASE-0: single grey band, single g-point ---------------------
        band(1)%planck_wt = 1.0_wp
        band(1)%ng        = 1
        band(1)%gw(1)     = 1.0_wp
        band(1)%kh2o(1)   = 0.55_wp
        band(1)%kco2(1)   = 0.30_wp
        band(1)%ko3(1)    = 40.0_wp
        return
    end subroutine ecckd_lw_table

    subroutine aeros_ecckd_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                              q_co2, l_o3, fnet, heat, olr, fdw_sur)
        ! Clear-sky longwave for one column via the correlated-k machinery.
        ! Signature is identical to aeros_lw_clearsky_column (the SESAM kernel) so
        ! the two are interchangeable behind the `scheme` selector.
        !
        ! Model ordering on input (k=1 top .. k=nlev surface). Fluxes on the
        ! nlev+1 interfaces, absorber amounts and the Planck source per layer.
        ! Local interface convention (as SESAM): i=0 surface, i=nlev TOA; layer
        ! l (surface->top) is model layer k = nlev-l+1.

        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: t(:)        ! (nlev) layer temperature [K]
        real(wp), intent(in)  :: q(:)        ! (nlev) specific humidity [kg kg-1]
        real(wp), intent(in)  :: o3(:)       ! (nlev) ozone mass mixing ratio [kg kg-1]
        real(wp), intent(in)  :: dp_lev(:)   ! (nlev) layer thickness [Pa], > 0
        real(wp), intent(in)  :: z_half(0:)  ! (0:nlev) interface height [m]
        real(wp), intent(in)  :: ts          ! surface skin temperature [K]
        real(wp), intent(in)  :: q_co2       ! CO2 mass mixing ratio [kg kg-1]
        logical,  intent(in)  :: l_o3        ! include ozone

        real(wp), intent(out) :: fnet(0:)    ! (0:nlev) net UPWARD LW flux [W m-2]
        real(wp), intent(out) :: heat(:)     ! (nlev) LW heating rate [K s-1]
        real(wp), intent(out) :: olr         ! outgoing LW at TOA [W m-2]
        real(wp), intent(out) :: fdw_sur     ! downward LW at surface [W m-2]

        type(ecckd_band) :: band(N_LW_BANDS)
        ! layer absorber mass paths, local surface->top order [g cm-2]
        real(wp) :: uwv(nlev), uco2(nlev), uo3(nlev)
        ! per-band accumulators of the net-upward flux on interfaces (model order)
        real(wp) :: fnet_b(0:nlev), olr_b, fdw_b
        integer  :: k, l, ib, ig

        call ecckd_lw_table(band)

        ! --- layer absorber mass paths (hydrostatic, g cm-2), local order ---
        ! Correlated-k acts on the absorber amount directly, so all three gases
        ! use the same simple mass path (0.1 q dp/g) -- no band-model z-weighting.
        do l = 1, nlev
            k = nlev - l + 1
            uwv(l)  = 0.1_wp * q(k)    * dp_lev(k)/grav
            uco2(l) = 0.1_wp * q_co2   * dp_lev(k)/grav
            if (l_o3) then
                uo3(l) = 0.1_wp * max(0.0_wp, o3(k)) * dp_lev(k)/grav
            else
                uo3(l) = 0.0_wp
            end if
        end do

        ! --- band / g-point loop: accumulate net-upward flux and boundaries ---
        fnet = 0.0_wp; olr = 0.0_wp; fdw_sur = 0.0_wp
        do ib = 1, N_LW_BANDS
            do ig = 1, band(ib)%ng
                call gpoint_flux(nlev, t, ts, uwv, uco2, uo3, &
                    band(ib)%planck_wt*band(ib)%gw(ig), &
                    band(ib)%kh2o(ig), band(ib)%kco2(ig), band(ib)%ko3(ig), &
                    fnet_b, olr_b, fdw_b)
                fnet    = fnet    + fnet_b
                olr     = olr     + olr_b
                fdw_sur = fdw_sur + fdw_b
            end do
        end do

        ! --- layer heating: divergence of net-upward flux -------------------
        do k = 1, nlev
            heat(k) = (grav/cp_d) * (fnet(k) - fnet(k-1))/dp_lev(k)
        end do
        return
    end subroutine aeros_ecckd_lw_clearsky_column

    subroutine gpoint_flux(nlev, t, ts, uwv, uco2, uo3, wsrc, kh2o, kco2, ko3, &
                           fnet, olr, fdw_sur)
        ! Emissivity-method net-upward LW flux on the interfaces for ONE g-point.
        ! The Planck source in this g-point is wsrc*sigma T^4 (wsrc = band Planck
        ! fraction times the g-point weight); transmission between two interfaces
        ! is exp(-beta0 * sum_k (kh2o u_h2o + kco2 u_co2 + ko3 u_o3)) over the
        ! layers strictly between them (correlated-k: k is constant across the
        ! g-point, so paths add and the exponential is evaluated on the total).
        !
        ! Returns net-upward flux in MODEL interface order (i=0 TOA .. i=nlev
        ! surface), plus the two boundary fluxes.

        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: t(:)          ! (nlev) model order (k=1 top)
        real(wp), intent(in)  :: ts
        real(wp), intent(in)  :: uwv(:), uco2(:), uo3(:)   ! local surface->top
        real(wp), intent(in)  :: wsrc, kh2o, kco2, ko3
        real(wp), intent(out) :: fnet(0:)      ! (0:nlev) net up, MODEL order
        real(wp), intent(out) :: olr, fdw_sur

        real(wp) :: b(nlev)                    ! Planck source per layer, local order
        real(wp) :: bsfc, fup(0:nlev), fdw(0:nlev)
        integer  :: k, l, i

        ! Planck source per layer (this g-point's share), local surface->top
        do l = 1, nlev
            k = nlev - l + 1
            b(l) = wsrc * sigma_sb * t(k)**4
        end do
        bsfc = wsrc * sigma_sb * ts**4

        ! upward flux at each local interface i (surface = 0)
        fup(0) = bsfc
        do i = 1, nlev
            fup(i) = bsfc*trans(0, i)
            do l = 1, i
                fup(i) = fup(i) + b(l)*(trans(l, i) - trans(l-1, i))
            end do
        end do
        ! downward flux (no incident LW at TOA)
        fdw(nlev) = 0.0_wp
        do i = nlev-1, 0, -1
            fdw(i) = 0.0_wp
            do l = i+1, nlev
                fdw(i) = fdw(i) + b(l)*(trans(l-1, i) - trans(l, i))
            end do
        end do

        ! net upward, back to model interface order (local i -> model nlev-i)
        do i = 0, nlev
            fnet(nlev - i) = fup(i) - fdw(i)
        end do
        olr     = fup(nlev)
        fdw_sur = fdw(0)
        return

    contains

        pure real(wp) function trans(ia, ib) result(tr)
            ! exp(-beta0 * accumulated absorber optical depth) between interfaces
            integer, intent(in) :: ia, ib
            integer  :: lo, hi, m
            real(wp) :: tau
            lo = min(ia, ib); hi = max(ia, ib)
            if (hi - lo <= 0) then
                tr = 1.0_wp; return
            end if
            tau = 0.0_wp
            do m = lo+1, hi
                tau = tau + kh2o*uwv(m) + kco2*uco2(m) + ko3*uo3(m)
            end do
            tr = exp(-LW_BETA0*tau)
            return
        end function trans

    end subroutine gpoint_flux

end module aeros_ecckd
