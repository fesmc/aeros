module aeros_ecckd
    ! ecCKD-style correlated-k longwave for aeros -- the endpoint named in
    ! design.md section 5 (option 1), reached via the `scheme` selector in
    ! aeros_radiation (SCHEME_ECCKD). SESAM (SCHEME_SESAM) stays the default and
    ! is untouched; this module is the opt-in path.
    !
    ! === Gas optics: a compact Malkmus band-model correlated-k (option 2b) ===
    !
    ! There is no ecCKD look-up table on this machine and the ecmwf/ecckd
    ! generator + CKDMIP line-by-line archive are not available offline
    ! (verified), so the gas optics are NOT a bundled binary LUT (option 2a).
    ! Instead each spectral band carries a Malkmus/Goody random-band model with a
    ! few hard-coded band parameters, from which a few-g-point k-distribution is
    ! derived analytically. No run-time data dependency; fully reproducible.
    !
    ! THE KEY IDENTITY. The Malkmus band transmission
    !     T(u) = exp{ -a [ sqrt(1 + (2 kbar/a) u) - 1 ] }
    ! (weak-line limit T -> exp(-kbar u); strong-line limit T -> exp(-sqrt(2 a kbar u)),
    ! the square-root growth that makes a saturated band's forcing logarithmic) is
    ! the Laplace transform of an INVERSE-GAUSSIAN (Wald) distribution in k, with
    ! mean mu = kbar and shape lambda = a*kbar. So the correlated-k g-points are
    ! got exactly, with no fitting: pick Gauss-Legendre nodes g_j in [0,1] and
    ! invert the analytic inverse-Gaussian CDF to k(g_j). That is the whole gas
    ! optics -- a handful of band parameters in, g-point (k,weight) tables out.
    !
    ! This is the structural fix over SESAM's single broadband transmission fit:
    !   - a real 8-12 um WINDOW band that stays transparent (so OLR escapes and
    !     the warm-moist-tropics OLR bias of SESAM, §14/§20, closes WITHOUT the
    !     LW_VAP_OPAC fudge -- vap-opac is 1.0 here, the whole point);
    !   - a real 15 um CO2 band whose strong-line sqrt-growth gives a logarithmic
    !     CO2 forcing near the canonical ~3.7 W/m2 (SESAM's broadband fit gives ~3.0);
    !   - gas overlap resolved band by band, not as a product of three broadband fits.
    !
    ! Bands and gases (LW, terrestrial 0..3000 cm-1), g-points chosen lean
    ! (2 where a band tolerates it, 3 where the k-range is wide):
    !   band 1   0- 540 cm-1  H2O rotation                     3 g
    !   band 2 540- 800 cm-1  CO2 15um (major) + H2O overlap    3 g
    !   band 3 800-1250 cm-1  window: H2O continuum + O3 9.6um  2 g
    !   band 4 1250-inf cm-1  H2O vib-rotation + CO2 wing        2 g
    ! => 10 longwave g-points total.
    !
    ! Temperature: the dominant T-dependence is the Planck fraction in each band
    ! (integrated exactly per layer, essentially free); on top of that each gas
    ! carries a mild analytic band-strength S(T) scaling of kbar. Pressure: cheap
    ! -- a power-law broadening of the absorber path (p/p0)^kappa per gas, so the
    ! k-distribution shape is built once at reference and only the path is scaled.
    !
    ! ERA5 anchoring: the a-priori Malkmus parameters carry the band SHAPE
    ! (relative g-point k, band strengths, and the Lorentz pressure exponents
    ! KAPPA -- all a-priori, NOT ERA5-fit); a single documented per-band scale on
    ! the major gas's kbar (KSCALE) sets the MAGNITUDE, plus the CO2 band strength
    ! that anchors the doubling. <= 2 anchored quantities per band.
    !
    ! Phase-1 result, validated on ERA5 1991-2020 clear-sky columns
    ! (drivers/validate_era5.f90 with arg3=ecckd), area-weighted global means,
    ! with LW_VAP_OPAC retired (the SESAM fudge is NOT used here):
    !   OLR bias        -5.1 W/m2   (SESAM -7.1 with the fudge, -15.9 without; §20)
    !   surface-down LW -7.8 W/m2   (SESAM -5.0 with the fudge)
    !   2xCO2 forcing    3.6 W/m2   (SESAM 3.0; canonical ~3.7)
    ! and ~1.8x the SESAM LW cost per column (10 g-points, O(nlev^2) solver).
    !
    ! Not yet here (later phases): the shortwave sibling and the all-sky
    ! (grey-cloud) branch -- aeros_radiation routes the ecCKD shortwave and every
    ! cloudy column through the SESAM kernels, so grey clouds keep working.

    use aeros_defs, only : wp, sigma_sb, grav, cp_d, p0, pi

    implicit none

    private

    ! === Spectral bands ====================================================
    integer,  parameter :: NB    = 4                 ! number of LW bands
    integer,  parameter :: NGAS  = 3                 ! 1=H2O, 2=CO2, 3=O3
    integer,  parameter :: NGT   = 3                 ! max g-points per band
    integer,  parameter :: NGB(NB) = [3, 3, 2, 2]    ! g-points in each band

    ! upper wavenumber edge of each band [cm-1]; lower edge of band 1 is 0, upper
    ! edge of band NB is +infinity (both handled in the Planck fraction).
    real(wp), parameter :: WN_HI(NB) = [540.0_wp, 800.0_wp, 1250.0_wp, 3000.0_wp]

    real(wp), parameter :: C2    = 1.438776_wp       ! hc/k [cm K]
    real(wp), parameter :: BETA0 = 1.66_wp           ! LW flux diffusivity factor
    real(wp), parameter :: TREF  = 260.0_wp          ! reference T for S(T) [K]

    ! === Malkmus reference band parameters =================================
    ! Per band (row) and gas (col: H2O, CO2, O3). A zero KBAR means the gas does
    ! not absorb in that band. These are compact band-averaged values (the a-priori
    ! SHAPE); KSCALE (below) is the ERA5-anchored magnitude knob on the major gas.
    !
    !   KBAR   band-mean mass absorption at reference [cm2 g-1], path in g cm-2
    !   ALS    Malkmus line-structure parameter a (small = strong-line/saturating,
    !          large = weak-line/grey)
    !   KAPPA  pressure-broadening exponent: path scaled by (p/p0)^kappa
    !   NST    band-strength temperature exponent: kbar scaled by (TREF/T)^nst
    !
    ! reshape fills column-major: each bracketed row below is one GAS across the
    ! four bands (H2O bands1-4, then CO2 bands1-4, then O3 bands1-4).
    !                                       band1     band2     band3     band4
    real(wp), parameter :: KBAR(NB,NGAS) = reshape([ &
        12.0_wp,   1.5_wp,   0.10_wp,   6.0_wp,   &   ! H2O
         0.0_wp, 150.0_wp,   0.0_wp,    2.0_wp,   &   ! CO2
         0.0_wp,   0.0_wp, 120.0_wp,    0.0_wp    &   ! O3
        ], [NB,NGAS])
    real(wp), parameter :: ALS(NB,NGAS) = reshape([ &
         0.40_wp,  1.00_wp,  6.00_wp,   0.50_wp,  &   ! H2O
         1.00_wp,  0.35_wp,  1.00_wp,   1.00_wp,  &   ! CO2
         1.00_wp,  1.00_wp,  0.60_wp,   1.00_wp   &   ! O3
        ], [NB,NGAS])
    ! a-priori (textbook Lorentz p-broadening, NOT ERA5-tuned): H2O ~1.0, CO2 0.7,
    ! O3 0.5. Structure comes from the band physics; magnitude from KSCALE below.
    real(wp), parameter :: KAPPA(NB,NGAS) = reshape([ &
         1.00_wp,  1.00_wp,  1.00_wp,   1.00_wp,  &   ! H2O
         0.00_wp,  0.70_wp,  0.00_wp,   0.70_wp,  &   ! CO2
         0.00_wp,  0.00_wp,  0.50_wp,   0.00_wp   &   ! O3
        ], [NB,NGAS])
    real(wp), parameter :: NST(NB,NGAS) = reshape([ &
         1.00_wp,  0.50_wp,  0.00_wp,   0.50_wp,  &   ! H2O rotation strengthens as T falls
         0.00_wp,  0.00_wp,  0.00_wp,   0.00_wp,  &   ! CO2
         0.00_wp,  0.00_wp,  0.00_wp,   0.00_wp   &   ! O3
        ], [NB,NGAS])

    ! ERA5-anchored per-band magnitude scale on each gas's kbar. The major gas of
    ! each band is the tuning knob (<=2 per band); others stay 1.0 (a-priori).
    ! Tuned to null the §14 global-mean clear-sky OLR bias (ERA5 264.1 W/m2) with
    ! the vapour-opacity fudge OFF. CO2 band-2 scale also anchors the doubling.
    real(wp), parameter :: KSCALE(NB,NGAS) = reshape([ &
         1.70_wp,  1.20_wp,  1.20_wp,   1.70_wp,  &   ! H2O
         1.00_wp,  1.00_wp,  1.00_wp,   1.00_wp,  &   ! CO2
         1.00_wp,  1.00_wp,  1.00_wp,   1.00_wp   &   ! O3
        ], [NB,NGAS])

    public :: aeros_ecckd_lw_clearsky_column

contains

    subroutine aeros_ecckd_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                              q_co2, l_o3, fnet, heat, olr, fdw_sur)
        ! Clear-sky longwave for one column via the Malkmus correlated-k. Signature
        ! is identical to aeros_lw_clearsky_column (the SESAM kernel) so the two are
        ! interchangeable behind the `scheme` selector. Model ordering on input
        ! (k=1 top .. k=nlev surface); fluxes on the nlev+1 interfaces.
        !
        ! Local convention (as SESAM): interface i=0 surface .. i=nlev TOA; local
        ! layer l (surface->top) is model layer k = nlev-l+1.

        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: t(:)        ! (nlev) layer temperature [K]
        real(wp), intent(in)  :: q(:)        ! (nlev) specific humidity [kg kg-1]
        real(wp), intent(in)  :: o3(:)       ! (nlev) ozone mass mixing ratio [kg kg-1]
        real(wp), intent(in)  :: dp_lev(:)   ! (nlev) layer thickness [Pa], > 0
        real(wp), intent(in)  :: z_half(0:)  ! (0:nlev) interface height [m] (kept for signature parity)
        real(wp), intent(in)  :: ts          ! surface skin temperature [K]
        real(wp), intent(in)  :: q_co2       ! CO2 mass mixing ratio [kg kg-1]
        logical,  intent(in)  :: l_o3        ! include ozone

        real(wp), intent(out) :: fnet(0:)    ! (0:nlev) net UPWARD LW flux [W m-2]
        real(wp), intent(out) :: heat(:)     ! (nlev) LW heating rate [K s-1]
        real(wp), intent(out) :: olr         ! outgoing LW at TOA [W m-2]
        real(wp), intent(out) :: fdw_sur     ! downward LW at surface [W m-2]

        ! g-point k-table (reference), built once per call: k [cm2 g-1] and weight
        real(wp) :: kj(NB,NGAS,NGT), gw(NB,NGT)
        ! per-layer absorber mass paths [g cm-2], local order (l=1 surface)
        real(wp) :: umass(nlev,NGAS)
        real(wp) :: pratio(nlev), tlay(nlev)            ! layer p/p0 and T, local order
        real(wp) :: pf_layer(nlev,NB)                   ! band Planck fraction per layer
        real(wp) :: bsrc(nlev)                          ! grey layer source sigma T^4, local
        real(wp) :: bsurf_b(NB)                         ! surface Planck fraction * sigma Ts^4
        real(wp) :: ctau(0:nlev)                        ! cumulative optical depth, this g-point
        real(wp) :: fup(0:nlev), fdw(0:nlev)
        real(wp) :: p_half, pfl, tau_l
        integer  :: k, l, i, ib, ig, igas

        call build_ktable(kj, gw)

        ! --- per-layer absorber paths, pressure ratio, temperature, Planck ------
        ! Reconstruct layer pressure from dp (top interface ~ 0): cheap-in-p.
        p_half = 0.0_wp
        do k = 1, nlev                       ! model order, k=1 top
            l   = nlev - k + 1               ! local order, l=1 surface
            pfl = p_half + 0.5_wp*dp_lev(k)  ! layer-mid pressure
            p_half = p_half + dp_lev(k)

            umass(l,1) = 0.1_wp * max(0.0_wp, q(k))  * dp_lev(k)/grav   ! H2O
            umass(l,2) = 0.1_wp * max(0.0_wp, q_co2) * dp_lev(k)/grav   ! CO2
            if (l_o3) then
                umass(l,3) = 0.1_wp * max(0.0_wp, o3(k)) * dp_lev(k)/grav ! O3
            else
                umass(l,3) = 0.0_wp
            end if

            pratio(l) = pfl/p0
            tlay(l)   = t(k)
            bsrc(l)   = sigma_sb*t(k)**4
            do ib = 1, NB
                pf_layer(l,ib) = planck_frac(ib, t(k))
            end do
        end do
        do ib = 1, NB
            bsurf_b(ib) = planck_frac(ib, ts)*sigma_sb*ts**4
        end do

        ! --- accumulate net-upward flux over bands and g-points ----------------
        fnet = 0.0_wp; olr = 0.0_wp; fdw_sur = 0.0_wp
        do ib = 1, NB
            do ig = 1, NGB(ib)
                ! cumulative optical depth from surface (local i=0) upward
                ctau(0) = 0.0_wp
                do l = 1, nlev
                    tau_l = 0.0_wp
                    do igas = 1, NGAS
                        if (kj(ib,igas,ig) <= 0.0_wp) cycle
                        tau_l = tau_l + kj(ib,igas,ig) * umass(l,igas) &
                              * pratio(l)**KAPPA(ib,igas) &
                              * (TREF/tlay(l))**NST(ib,igas)
                    end do
                    ctau(l) = ctau(l-1) + tau_l
                end do

                call band_flux(nlev, gw(ib,ig), pf_layer(:,ib), bsrc, &
                               bsurf_b(ib), ctau, fup, fdw)

                do i = 0, nlev
                    fnet(nlev - i) = fnet(nlev - i) + (fup(i) - fdw(i))
                end do
                olr     = olr     + fup(nlev)
                fdw_sur = fdw_sur + fdw(0)
            end do
        end do

        ! --- layer heating: divergence of net-upward flux ----------------------
        do k = 1, nlev
            heat(k) = (grav/cp_d) * (fnet(k) - fnet(k-1))/dp_lev(k)
        end do
        return
    end subroutine aeros_ecckd_lw_clearsky_column

    subroutine band_flux(nlev, wg, pf_layer, bsrc, bsurf, ctau, fup, fdw)
        ! Emissivity-method net LW flux on the interfaces for one band+g-point.
        ! Source in this g-point for local layer l is wg*pf_layer(l)*bsrc(l);
        ! surface is wg*bsurf. Transmission between local interfaces from the
        ! cumulative optical depth (correlated-k: k constant across the g-point,
        ! so optical depths add and the exponential is on the total path).
        integer,  intent(in)  :: nlev
        real(wp), intent(in)  :: wg
        real(wp), intent(in)  :: pf_layer(:)    ! (nlev) band Planck fraction, local order
        real(wp), intent(in)  :: bsrc(:)        ! (nlev) grey layer source sigma T^4, local order
        real(wp), intent(in)  :: bsurf          ! band Planck fraction * sigma Ts^4
        real(wp), intent(in)  :: ctau(0:)       ! (0:nlev) cumulative optical depth (surface=0)
        real(wp), intent(out) :: fup(0:), fdw(0:)

        real(wp) :: b(nlev), bs
        integer  :: l, i

        do l = 1, nlev
            b(l) = wg*pf_layer(l)*bsrc(l)
        end do
        bs = wg*bsurf

        fup(0) = bs
        do i = 1, nlev
            fup(i) = bs*tr(0, i)
            do l = 1, i
                fup(i) = fup(i) + b(l)*(tr(l, i) - tr(l-1, i))
            end do
        end do
        fdw(nlev) = 0.0_wp
        do i = nlev-1, 0, -1
            fdw(i) = 0.0_wp
            do l = i+1, nlev
                fdw(i) = fdw(i) + b(l)*(tr(l-1, i) - tr(l, i))
            end do
        end do
        return
    contains
        pure real(wp) function tr(ia, ib) result(x)
            integer, intent(in) :: ia, ib
            x = exp(-BETA0*abs(ctau(max(ia,ib)) - ctau(min(ia,ib))))
        end function tr
    end subroutine band_flux

    ! === Correlated-k table construction ===================================

    subroutine build_ktable(kj, gw)
        ! Build reference g-point absorption coefficients and weights for every
        ! band and gas from the Malkmus parameters, via the inverse-Gaussian CDF.
        real(wp), intent(out) :: kj(NB,NGAS,NGT), gw(NB,NGT)
        real(wp) :: g(NGT), w(NGT), mu, lam
        integer  :: ib, ig, ng, igas

        kj = 0.0_wp; gw = 0.0_wp
        do ib = 1, NB
            ng = NGB(ib)
            call gauss_legendre01(ng, g, w)
            do ig = 1, ng
                gw(ib,ig) = w(ig)
            end do
            do igas = 1, NGAS
                mu = KSCALE(ib,igas)*KBAR(ib,igas)     ! mean k = kbar (scaled)
                if (mu <= 0.0_wp) cycle
                lam = ALS(ib,igas)*mu                  ! shape lambda = a*kbar
                do ig = 1, ng
                    kj(ib,igas,ig) = invgauss_icdf(g(ig), mu, lam)
                end do
            end do
        end do
        return
    end subroutine build_ktable

    pure subroutine gauss_legendre01(n, g, w)
        ! Gauss-Legendre nodes/weights mapped to [0,1] (n = 2 or 3).
        integer,  intent(in)  :: n
        real(wp), intent(out) :: g(:), w(:)
        g = 0.0_wp; w = 0.0_wp
        select case (n)
        case (2)
            g(1) = 0.5_wp*(1.0_wp - 0.5773502691896257_wp)
            g(2) = 0.5_wp*(1.0_wp + 0.5773502691896257_wp)
            w(1) = 0.5_wp; w(2) = 0.5_wp
        case (3)
            g(1) = 0.5_wp*(1.0_wp - 0.7745966692414834_wp)
            g(2) = 0.5_wp
            g(3) = 0.5_wp*(1.0_wp + 0.7745966692414834_wp)
            w(1) = 0.2777777777777778_wp
            w(2) = 0.4444444444444444_wp
            w(3) = 0.2777777777777778_wp
        case default
            g(1) = 0.5_wp; w(1) = 1.0_wp
        end select
    end subroutine gauss_legendre01

    pure real(wp) function invgauss_icdf(pin, mu, lam) result(k)
        ! Inverse CDF of the inverse-Gaussian (Wald) distribution with mean mu and
        ! shape lam, by bisection in log(k). This is the Malkmus k-distribution:
        ! g-point absorption coefficient at cumulative probability p.
        real(wp), intent(in) :: pin, mu, lam
        real(wp) :: p, lo, hi, mid, fmid
        integer  :: it
        p  = min(max(pin, 1.0e-6_wp), 1.0_wp - 1.0e-6_wp)
        lo = log(mu) - 20.0_wp          ! k in [mu*e^-20, mu*e^+20]
        hi = log(mu) + 20.0_wp
        do it = 1, 80
            mid  = 0.5_wp*(lo + hi)
            fmid = invgauss_cdf(exp(mid), mu, lam) - p
            if (fmid > 0.0_wp) then
                hi = mid
            else
                lo = mid
            end if
        end do
        k = exp(0.5_wp*(lo + hi))
    end function invgauss_icdf

    pure real(wp) function invgauss_cdf(k, mu, lam) result(f)
        ! CDF of the inverse-Gaussian: F(k)=Phi(sqrt(lam/k)(k/mu-1))
        !                                   + exp(2 lam/mu) Phi(-sqrt(lam/k)(k/mu+1))
        real(wp), intent(in) :: k, mu, lam
        real(wp) :: s, a1, a2, ex
        if (k <= 0.0_wp) then
            f = 0.0_wp; return
        end if
        s  = sqrt(lam/k)
        a1 = s*(k/mu - 1.0_wp)
        a2 = s*(k/mu + 1.0_wp)
        ex = min(2.0_wp*lam/mu, 60.0_wp)          ! guard overflow
        f  = phi(a1) + exp(ex)*phi(-a2)
        f  = min(max(f, 0.0_wp), 1.0_wp)
    end function invgauss_cdf

    pure real(wp) function phi(z) result(p)
        ! standard normal CDF via erfc
        real(wp), intent(in) :: z
        p = 0.5_wp*erfc(-z/sqrt(2.0_wp))
    end function phi

    ! === Planck fraction in a band =========================================

    pure real(wp) function planck_frac(ib, t) result(pf)
        ! Fraction of the blackbody flux at temperature t that falls in band ib.
        ! Band 1 lower edge is 0 (fraction 1 below it); band NB upper edge is inf
        ! (fraction 0 above it). Uses the exact fractional-emissive-power series.
        integer,  intent(in) :: ib
        real(wp), intent(in) :: t
        real(wp) :: flo, fhi
        if (ib == 1) then
            flo = 1.0_wp                          ! from 0 cm-1: whole blackbody
        else
            flo = bb_frac_above(C2*WN_HI(ib-1)/t)
        end if
        if (ib == NB) then
            fhi = 0.0_wp                          ! to +inf
        else
            fhi = bb_frac_above(C2*WN_HI(ib)/t)
        end if
        pf = max(0.0_wp, flo - fhi)
    end function planck_frac

    pure real(wp) function bb_frac_above(x) result(f)
        ! Fraction of blackbody emissive power at wavenumbers ABOVE x = C2*nu/T,
        ! i.e. integral_x^inf of the normalized Planck function. Standard series
        ! f(x) = (15/pi^4) sum_{n>=1} e^{-nx}(x^3/n + 3x^2/n^2 + 6x/n^3 + 6/n^4).
        real(wp), intent(in) :: x
        real(wp) :: s, en, xn
        integer  :: n
        if (x <= 0.0_wp) then
            f = 1.0_wp; return
        end if
        s = 0.0_wp
        do n = 1, 8
            xn = real(n, wp)
            en = exp(-xn*x)
            s  = s + en*( x*x*x/xn + 3.0_wp*x*x/xn**2 + 6.0_wp*x/xn**3 + 6.0_wp/xn**4 )
        end do
        f = (15.0_wp/pi**4)*s
        f = min(max(f, 0.0_wp), 1.0_wp)
    end function bb_frac_above

end module aeros_ecckd
