module aeros_vordiv
    ! Vorticity/divergence <-> (u,v).
    !
    ! The primitive-equation core carries relative vorticity and divergence as
    ! its momentum prognostics rather than the winds themselves, because that
    ! is what makes the semi-implicit gravity-wave solve diagonal in degree l
    ! (docs/design.md section 3.2). The physics, the advection and every
    ! diagnostic want (u, v). This module is the translation, and it is the
    ! only place in aeros where the two representations meet.
    !
    ! === The mapping =========================================================
    !
    ! Helmholtz decomposition of a horizontal wind field on a sphere of radius
    ! a, with streamfunction psi and velocity potential chi:
    !
    !     v = grad_h(chi) + k x grad_h(psi)
    !     D = laplacian_h(chi),   zeta = laplacian_h(psi)
    !
    ! SHTns works in colatitude theta with its own spheroidal/toroidal
    ! potentials S and T, on the UNIT sphere:
    !
    !     V_theta = dS/dtheta + (1/sin theta) dT/dphi
    !     V_phi   = (1/sin theta) dS/dphi - dT/dtheta
    !
    ! Note the SIGN of the toroidal term, and that theta increases SOUTHWARD,
    ! so V_theta is a southward velocity. Converting:
    !
    !     u = V_phi        v = -V_theta
    !     S = chi/a        T = -psi/a
    !
    ! and, using laplacian(Y_lm) = -l(l+1)/a^2 Y_lm,
    !
    !     S_lm = -a/(l(l+1)) * D_lm        D_lm = -l(l+1)/a * S_lm
    !     T_lm = +a/(l(l+1)) * zeta_lm     zeta_lm = +l(l+1)/a * T_lm
    !
    ! SHTns' documentation does not state the sign convention in a form that
    ! could be relied on, so it was MEASURED rather than assumed: synthesizing
    ! T_10 = 1 alone gives V_phi = +N sin(theta) with N = sqrt(3/4pi), and
    ! S_10 = 1 alone gives V_theta = -N sin(theta), both to the digit. The
    ! solid-body rotation test in tests/test_vordiv.f90 is what keeps that
    ! true; if SHTns ever flips a sign, that test fails rather than the model
    ! quietly running backwards.
    !
    ! === l = 0 ===============================================================
    !
    ! There is no degree-0 horizontal vector field on a sphere: a globally
    ! uniform wind cannot exist. So S_00 = T_00 = 0 always, and any zeta_00 or
    ! D_00 handed in is DISCARDED rather than being an error. That is physical,
    ! not a convenience -- the global means of vorticity and divergence are
    ! identically zero -- but it does mean the round trip is the identity on
    ! l >= 1 only.
    !
    ! === Robert form =========================================================
    !
    ! SHTns exposes shtns_robert_form(), which makes vector transforms work
    ! with u*sin(theta) instead of u -- the U = u cos(latitude) formulation
    ! classical spectral models use, which avoids dividing by cos(latitude)
    ! near the poles. The symbol IS exported (unlike shtns_set_batch, which the
    ! bundled header declares but the library does not have).
    !
    ! It is deliberately NOT enabled here. Turning it on is a property of the
    ! CONFIG, so it would silently change what every vector transform in the
    ! model returns, including the raw ones in aeros_spectral. The place to
    ! decide it is M1.3, where the nonlinear tendency terms fix which form is
    ! actually wanted; enabling it before then would be choosing a formulation
    ! by accident.

    use aeros_defs,     only : dp, wp, wp_sh, r_earth, io_unit_err
    use aeros_spectral, only : aeros_sht_class, aeros_sht_analysis_vec, &
                                aeros_sht_synthesis_vec

    implicit none

    private

    public :: aeros_uv_from_vordiv
    public :: aeros_vordiv_from_uv

contains

    subroutine aeros_uv_from_vordiv(sht, vor, div, u, v)
        ! Spectral vorticity and divergence -> grid winds, one level.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        complex(wp_sh), intent(in), contiguous :: vor(:)  ! relative vorticity [s-1], (nlm)
        complex(wp_sh), intent(in), contiguous :: div(:)  ! divergence [s-1], (nlm)
        real(wp), intent(out), contiguous :: u(:,:)       ! zonal wind [m s-1], (nlon,nlat)
        real(wp), intent(out), contiguous :: v(:,:)       ! meridional wind [m s-1]

        complex(wp_sh) :: slm(sht%nlm), tlm(sht%nlm)
        real(wp)       :: vth(sht%nlon, sht%nlat)
        real(dp)       :: fac
        integer        :: lm, l, i, j

        do lm = 1, sht%nlm
            l = sht%l_of_lm(lm)
            if (l == 0) then
                slm(lm) = (0.0_wp_sh, 0.0_wp_sh)
                tlm(lm) = (0.0_wp_sh, 0.0_wp_sh)
            else
                fac = r_earth/real(l*(l+1), dp)
                slm(lm) = -fac*div(lm)
                tlm(lm) =  fac*vor(lm)
            end if
        end do

        ! v_phi IS u, so it is written straight into the output array; v_theta
        ! is southward and needs the sign flip.
        call aeros_sht_synthesis_vec(sht, slm, tlm, vth, u)

        do j = 1, sht%nlat
            do i = 1, sht%nlon
                v(i,j) = -vth(i,j)
            end do
        end do

        return

    end subroutine aeros_uv_from_vordiv

    subroutine aeros_vordiv_from_uv(sht, u, v, vor, div)
        ! Grid winds -> spectral vorticity and divergence, one level.
        !
        ! Unlike the raw transforms in aeros_spectral, this does NOT destroy
        ! its inputs: the (-v, u) -> (V_theta, V_phi) relabelling needs a copy
        ! anyway, so intent(in) is free and the safer contract.

        implicit none

        type(aeros_sht_class), intent(in) :: sht
        real(wp), intent(in), contiguous :: u(:,:)         ! (nlon,nlat)
        real(wp), intent(in), contiguous :: v(:,:)
        complex(wp_sh), intent(out), contiguous :: vor(:)  ! (nlm)
        complex(wp_sh), intent(out), contiguous :: div(:)

        complex(wp_sh) :: slm(sht%nlm), tlm(sht%nlm)
        real(wp)       :: vth(sht%nlon, sht%nlat)
        real(wp)       :: vph(sht%nlon, sht%nlat)
        real(dp)       :: fac
        integer        :: lm, l, i, j

        do j = 1, sht%nlat
            do i = 1, sht%nlon
                vth(i,j) = -v(i,j)
                vph(i,j) =  u(i,j)
            end do
        end do

        call aeros_sht_analysis_vec(sht, vth, vph, slm, tlm)

        do lm = 1, sht%nlm
            l = sht%l_of_lm(lm)
            if (l == 0) then
                vor(lm) = (0.0_wp_sh, 0.0_wp_sh)
                div(lm) = (0.0_wp_sh, 0.0_wp_sh)
            else
                fac = real(l*(l+1), dp)/r_earth
                div(lm) = -fac*slm(lm)
                vor(lm) =  fac*tlm(lm)
            end if
        end do

        return

    end subroutine aeros_vordiv_from_uv

end module aeros_vordiv
