program probe_lw
    ! Offline single-column checks on the clear-sky radiation kernel, isolated
    ! from the coupled model -- the harness that settled the M2 RCE runaway.
    ! Four probes:
    !   run(nlev)      smooth analytic sounding -> per-layer heat and interface
    !                  fnet, at L12 and L20. Smooth heat here rules out a scheme
    !                  sawtooth (the coupled 'sawtooth' was a diagnostic bug).
    !   run_coupled    the coupled hot column (hot lowest layer, cold aloft) ->
    !                  LW heat and OLR, showing a trapped-surface profile emits
    !                  little (OLR set by the cold radiating levels aloft).
    !   run_sweep      uniform warming at fixed RH -> OLR(Ts): the scheme's true
    !                  dOLR/dTs (~2.75 W/m2/K, no saturation) -- so a flat coupled
    !                  OLR is a transport problem, not radiation.
    ! Not an acceptance test; a diagnostic, alongside drivers/rce_long.f90.
    use aeros_defs,         only : dp, wp, p0, grav, R_d, cp_d, sigma_sb, &
                                   aeros_grid_class
    use aeros_vertical
    use aeros_spectral
    use aeros_grid
    use aeros_condensation, only : aeros_qsat
    use aeros_radiation,    only : aeros_lw_clearsky_column
    implicit none

    call run(12)
    call run(20)
    call run_coupled()
    call run_sweep()

contains

    subroutine run_sweep()
        ! dOLR/dTs of the LW scheme: warm the whole sounding uniformly, hold RH
        ! fixed (so water vapour rises with T -- the real feedback), print OLR.
        ! A healthy scheme gives ~2 W/m2/K; saturation near a low ceiling is the
        ! runaway-greenhouse pathology (absorbed SW ~240-280 can never be shed).
        integer, parameter :: nlev = 20
        type(aeros_sht_pool_class), target :: pool
        type(aeros_grid_class)  :: grd
        type(aeros_vgrid_class) :: vg
        real(wp) :: phalf(0:nlev), pfull(nlev), dp_lev(nlev)
        real(wp) :: phi_full(nlev), phi_half(0:nlev), z_half(0:nlev)
        real(wp) :: t(nlev), q(nlev), o3(nlev), fnet(0:nlev), heat(nlev)
        real(wp) :: olr, fdw_sur, ts, ps, qco2, tcwv, off
        integer  :: k, io

        call aeros_sht_pool_init(pool, 21, quick=.TRUE.)
        call aeros_grid_init(grd, pool%sht(1))
        call aeros_vgrid_init(vg, nlev)
        ps = real(p0, wp)
        call aeros_vgrid_pressure(vg, ps, phalf, pfull, dp_lev)
        o3 = 0.0_wp
        qco2 = 280.0e-6_wp * 44.0095_wp/28.97_wp

        write(*,"(a)") "=== OLR sweep (uniform warming, RH fixed), L20"
        write(*,"(a)") "  dTs[K]   Ts[K]   TCWV[kg/m2]   OLR[W/m2]"
        do io = 0, 8
            off = real(io,wp)*10.0_wp
            do k = 1, nlev
                t(k) = strat_trop(pfull(k)/ps) + off
                q(k) = humid(t(k), pfull(k), pfull(k)/ps)
            end do
            ts = 288.0_wp + off
            call aeros_hydrostatic(vg, 0.0_wp, t, phalf, phi_full, phi_half)
            z_half = phi_half/real(grav, wp)
            call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                          qco2, .FALSE., fnet, heat, olr, fdw_sur)
            tcwv = 0.0_wp
            do k = 1, nlev
                tcwv = tcwv + q(k)*dp_lev(k)/real(grav,wp)
            end do
            write(*,"(f8.1,f8.1,f13.2,f13.2)") off, ts, tcwv, olr
        end do
        write(*,*) ""
        call aeros_vgrid_end(vg); call aeros_grid_end(grd)
        call aeros_sht_pool_end(pool)
    end subroutine run_sweep

    subroutine run_coupled()
        ! Feed the EXACT coupled day-5 hot-column zonal-mean T/q (from rce_long)
        ! through the same column routine, on the L12 vgrid. If this sawtooths,
        ! the profile shape causes it; if smooth, an input outside T/q does.
        integer, parameter :: nlev = 12
        type(aeros_sht_pool_class), target :: pool
        type(aeros_grid_class)  :: grd
        type(aeros_vgrid_class) :: vg
        real(wp) :: phalf(0:nlev), pfull(nlev), dp_lev(nlev)
        real(wp) :: phi_full(nlev), phi_half(0:nlev), z_half(0:nlev)
        real(wp) :: o3(nlev), fnet(0:nlev), heat(nlev), olr, fdw_sur
        real(wp) :: ts, ps, qco2
        integer  :: k
        real(wp) :: t(nlev), q(nlev)
        t = [218.29_wp,236.52_wp,244.72_wp,253.84_wp,262.38_wp,270.26_wp, &
             278.20_wp,283.48_wp,286.30_wp,288.17_wp,289.43_wp,290.41_wp]
        q = [0.2883_wp,0.7576_wp,1.0638_wp,1.7249_wp,2.6856_wp,3.8177_wp, &
             7.9297_wp,10.1355_wp,10.9115_wp,11.3498_wp,11.6442_wp,11.8594_wp]*1.0e-3_wp

        call aeros_sht_pool_init(pool, 21, quick=.TRUE.)
        call aeros_grid_init(grd, pool%sht(1))
        call aeros_vgrid_init(vg, nlev)
        ps = real(p0, wp)
        call aeros_vgrid_pressure(vg, ps, phalf, pfull, dp_lev)
        ts = 290.41_wp; o3 = 0.0_wp
        qco2 = 280.0e-6_wp * 44.0095_wp/28.97_wp
        call aeros_hydrostatic(vg, 0.0_wp, t, phalf, phi_full, phi_half)
        z_half = phi_half/real(grav, wp)
        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      qco2, .FALSE., fnet, heat, olr, fdw_sur)
        write(*,"(a,f7.2)") "=== COUPLED day-5 profile, L12   OLR ", olr
        write(*,"(a)") "  k     T[K]   q[g/kg]    fnet_below   heat[K/day]"
        do k = 1, nlev
            write(*,"(i3,f9.2,f9.3,f14.3,f13.3)") k, t(k), q(k)*1000.0_wp, &
                fnet(k), heat(k)*86400.0_wp
        end do
        write(*,*) ""
        call aeros_vgrid_end(vg); call aeros_grid_end(grd)
        call aeros_sht_pool_end(pool)
    end subroutine run_coupled

    subroutine run(nlev)
        integer, intent(in) :: nlev
        type(aeros_sht_pool_class), target :: pool
        type(aeros_grid_class)  :: grd
        type(aeros_vgrid_class) :: vg
        real(wp) :: phalf(0:nlev), pfull(nlev), dp_lev(nlev)
        real(wp) :: phi_full(nlev), phi_half(0:nlev), z_half(0:nlev)
        real(wp) :: t(nlev), q(nlev), o3(nlev)
        real(wp) :: fnet(0:nlev), heat(nlev), olr, fdw_sur
        real(wp) :: ts, ps, qco2
        integer  :: k

        call aeros_sht_pool_init(pool, 21, quick=.TRUE.)
        call aeros_grid_init(grd, pool%sht(1))
        call aeros_vgrid_init(vg, nlev)

        ps = real(p0, wp)
        call aeros_vgrid_pressure(vg, ps, phalf, pfull, dp_lev)
        do k = 1, nlev
            t(k) = strat_trop(pfull(k)/ps)
            q(k) = humid(t(k), pfull(k), pfull(k)/ps)
        end do
        ts = 288.0_wp
        o3 = 0.0_wp
        qco2 = 280.0e-6_wp * 44.0095_wp/28.97_wp

        call aeros_hydrostatic(vg, 0.0_wp, t, phalf, phi_full, phi_half)
        z_half = phi_half/real(grav, wp)

        call aeros_lw_clearsky_column(nlev, t, q, o3, dp_lev, z_half, ts, &
                                      qco2, .FALSE., fnet, heat, olr, fdw_sur)

        write(*,"(a,i0,a,f7.2,a,f7.2)") "=== L", nlev, "  OLR ", olr, &
            "  fdw_sur ", fdw_sur
        write(*,"(a)") "  k   pfull[hPa]   T[K]    q[g/kg]   dp[hPa]   " // &
                       "fnet_below   heat[K/day]"
        do k = 1, nlev
            write(*,"(i3,f11.2,f8.2,f9.3,f10.2,f12.3,f13.3)") k, pfull(k)/100.0_wp, &
                t(k), q(k)*1000.0_wp, dp_lev(k)/100.0_wp, fnet(k), heat(k)*86400.0_wp
        end do
        write(*,"(a,f12.3)") "  fnet(0)=TOA up = ", fnet(0)
        write(*,*) ""

        call aeros_vgrid_end(vg); call aeros_grid_end(grd)
        call aeros_sht_pool_end(pool)
    end subroutine run

    pure real(wp) function strat_trop(sig) result(tt)
        real(wp), intent(in) :: sig
        real(wp), parameter :: t_sfc=288.0_wp, t_str=216.5_wp, sig_tp=0.2_wp
        if (sig > sig_tp) then
            tt = t_str + (t_sfc - t_str)*(sig - sig_tp)/(1.0_wp - sig_tp)
        else
            tt = t_str
        end if
    end function strat_trop

    real(wp) function humid(tk, pk, sig) result(qq)
        real(wp), intent(in) :: tk, pk, sig
        real(wp) :: qs, dqsdt, rh
        call aeros_qsat(tk, pk, qs, dqsdt)
        rh = 0.75_wp*max(0.0_wp, (sig - 0.15_wp)/(1.0_wp - 0.15_wp))
        qq = max(3.0e-6_wp, rh*qs)
    end function humid

end program probe_lw
