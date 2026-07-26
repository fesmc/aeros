module aeros_correction
    ! The pluggable additive tendency-correction layer.
    !
    ! This is the framework docs/M1_scope.md names as its M1-specific
    ! requirement and docs/design.md sections 3.7 and 3.8 are the reason for.
    ! At M2 it carries no real correction -- the only term is a validated
    ! no-op -- and that is the point: the plumbing, the namelist, the
    ! diagnostics and the `correction = .FALSE.` twin are built and exercised
    ! now, so that M2b plugs section 3.7's DeltaF_res(t) and section 3.8's
    ! DeltaF_bias into an interface that already works.
    !
    !     DeltaF_total = sum_k DeltaF_k
    !
    ! Each term carries which prognostic fields it acts on, a spectral-scale
    ! selector, and its own on/off switch, because section 3.7 is explicit that
    ! a blanket correction on all fields "is far more likely to destabilise and
    ! much harder to diagnose", and because mwm/B_multires/B2B_OUTCOME.md found
    ! the naive whole-field version numerically fragile in practice.
    !
    ! === Where it is applied, and why there ==================================
    !
    ! IN SPECTRAL SPACE, on the assembled tendency, immediately after
    ! aeros_tendency_spectral. This is deliberately NOT the grid-space seam
    ! where the physics goes.
    !
    ! Physics has to be applied on the grid because it is nonlinear in the
    ! state: a horizontally varying drag does not commute with the curl, so
    ! -k_v v is not -k_v zeta and -k_v D (see aeros_tendency_grid). A
    ! correction is a different object. It is ADDITIVE AND LINEAR -- a fixed
    ! field added to a tendency -- so it commutes with everything, and applying
    ! it spectrally costs nothing in correctness.
    !
    ! What it buys:
    !
    !   1. The scale selector is free. Section 3.7 wants the correction applied
    !      only at large scales (the stationary-wave error it targets lives
    !      there), and B2B_OUTCOME.md lists spectral filtering as the first
    !      mitigation for the instability it hit. In spectral space that is a
    !      mask on total wavenumber. In grid space it would be a transform
    !      round-trip per field per step.
    !   2. Correction and physics stay independently switchable. They are
    !      different claims about the model and should not share a seam.
    !   3. It matches how DeltaF is diagnosed -- a coarse-grained difference of
    !      two models' tendencies, section 3.7 step 4.
    !
    ! The cost, stated plainly: a correction that genuinely had to be nonlinear
    ! in grid space could not use this path. Neither section 3.7 nor section
    ! 3.8 describes one.
    !
    ! === Units: aeros' tendency-scaling convention ===========================
    !
    ! A correction term is stored in THE SAME UNITS AS THE TENDENCY IT IS ADDED
    ! TO, with no scaling factor of any kind:
    !
    !     vor   [s-2]      div  [s-2]
    !     temp  [K s-1]    lnps [s-1]
    !
    ! This is aeros' analogue of the convention mwm/B_multires had to recover
    ! from SpeedyWeather's source (B2B_OUTCOME.md: spectral vor/div scaled by
    ! R^2, pres by R^1, grid temp/humid in physical units, and getting it wrong
    ! silently produced a wrong-magnitude correction). aeros carries unscaled
    ! physical tendencies throughout, so the convention is "no convention" --
    ! which is worth writing down precisely because the last model needed a
    ! page to explain its own.
    !
    ! === What is NOT here ====================================================
    !
    ! Diagnosing DeltaF. That is M2b: it needs two truncations run under
    ! identical boundary conditions and a coarse-graining operator, none of
    ! which exists yet. This module is the consumer of that result.
    !
    ! Flux-form conservation verification (section 3.7 risk 2, "verify to
    ! machine precision"). aeros_correction_report measures the channel through
    ! which a correction injects mass -- the (0,0) coefficient of the lnps
    ! correction, which is the only spectral mode that moves the global mean of
    ! ln(p_s) -- and that is the right diagnostic to watch. But a full
    ! statement about mass requires the grid-space integral int p_s dlnps/dt dA
    ! (mass is a nonlinear functional of the prognostic; see the mass fixer in
    ! aeros_timestep), and it is not worth building against a term that is
    ! identically zero. It belongs with the first real DeltaF.

    use aeros_defs,     only : dp, wp, wp_sh, io_unit_err, &
                                aeros_spec_class
    use aeros_spectral, only : aeros_sht_class
    use aeros_tendency, only : aeros_tend_class

    use nml,            only : nml_read

    implicit none

    private

    ! The number of terms the namelist can carry. Section 3.8 needs two
    ! (DeltaF_res and DeltaF_bias); the rest is headroom for diagnosing them in
    ! pieces, which section 3.8 mitigation 1 effectively requires.
    integer, parameter, public :: NTERM_MAX = 8

    type aeros_corr_term_class
        character(len=32) :: name = ""

        logical :: enabled = .FALSE.

        ! Which prognostic fields this term acts on (section 3.7, "correct
        ! selected terms, not everything").
        logical :: on_vor  = .FALSE.
        logical :: on_div  = .FALSE.
        logical :: on_temp = .FALSE.
        logical :: on_lnps = .FALSE.

        ! Spectral-scale selector: apply only where total wavenumber l <= lcut.
        ! Negative means every resolved scale. Resolved into `mask` at init so
        ! the inner loop is a lookup rather than a comparison.
        integer :: lcut = -1
        logical, allocatable :: mask(:)          ! (nlm)

        ! The correction itself, in tendency units. Allocated only for the
        ! fields the term acts on, so an unused field costs nothing and a bug
        ! that writes to one is a crash rather than a silent contribution.
        complex(wp_sh), allocatable :: vor(:,:)  ! (nlm,nlev)
        complex(wp_sh), allocatable :: div(:,:)
        complex(wp_sh), allocatable :: temp(:,:)
        complex(wp_sh), allocatable :: lnps(:)   ! (nlm)
    end type aeros_corr_term_class

    type aeros_correction_class
        ! The master switch IS section 3.8 mitigation 3's `DeltaF = 0` twin:
        ! .FALSE. must reproduce the uncorrected model bit for bit, and
        ! tests/test_correction.f90 asserts exactly that.
        logical :: enabled = .FALSE.

        integer :: nterm = 0
        type(aeros_corr_term_class) :: term(NTERM_MAX)

        integer :: nlm  = 0
        integer :: nlev = 0
    end type aeros_correction_class

    public :: aeros_corr_term_class
    public :: aeros_correction_class
    public :: aeros_correction_init
    public :: aeros_correction_load
    public :: aeros_correction_end
    public :: aeros_correction_add_term
    public :: aeros_correction_apply
    public :: aeros_correction_report

contains

    subroutine aeros_correction_init(cor, sht, nlev)
        ! An empty, disabled correction layer sized to the model.

        implicit none

        type(aeros_correction_class), intent(inout) :: cor
        type(aeros_sht_class),        intent(in)    :: sht
        integer,                      intent(in)    :: nlev

        call aeros_correction_end(cor)

        cor%enabled = .FALSE.
        cor%nterm   = 0
        cor%nlm     = sht%nlm
        cor%nlev    = nlev

        return

    end subroutine aeros_correction_init

    subroutine aeros_correction_load(cor, filename, sht, nlev, defaults_file)
        ! Configure from the `aeros_correction` namelist group.
        !
        ! The namelist decides which terms exist, what they act on and at what
        ! scales. It does NOT carry the correction VALUES -- those are fields,
        ! not parameters, and at M2b they arrive from a diagnosis run through
        ! the term's allocated arrays. A term configured here and never filled
        ! is a no-op, which is the M2 configuration.

        implicit none

        type(aeros_correction_class), intent(inout) :: cor
        character(len=*),             intent(in)    :: filename
        character(len=*), intent(in), optional :: defaults_file
        type(aeros_sht_class),        intent(in)    :: sht
        integer,                      intent(in)    :: nlev

        character(len=32) :: names(NTERM_MAX)
        logical :: enab(NTERM_MAX), ovor(NTERM_MAX), odiv(NTERM_MAX)
        logical :: otemp(NTERM_MAX), olnps(NTERM_MAX)
        integer :: lcut(NTERM_MAX)
        integer :: nterm, k, it

        call aeros_correction_init(cor, sht, nlev)

        call nml_read(filename, "aeros_correction", "correction",   cor%enabled, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "nterm",        nterm, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "term_name",    names, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "term_enabled", enab, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "term_on_vor",  ovor, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "term_on_div",  odiv, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "term_on_temp", otemp, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "term_on_lnps", olnps, defaults_file=defaults_file)
        call nml_read(filename, "aeros_correction", "term_lcut",    lcut, defaults_file=defaults_file)

        if (nterm < 0 .or. nterm > NTERM_MAX) then
            write(io_unit_err,*) "aeros_correction_load:: error: nterm must be in [0,", &
                                    NTERM_MAX, "], got ", nterm
            error stop 1
        end if

        do k = 1, nterm
            it = aeros_correction_add_term(cor, sht, trim(names(k)), &
                                            on_vor  = ovor(k),  on_div  = odiv(k), &
                                            on_temp = otemp(k), on_lnps = olnps(k), &
                                            lcut = lcut(k), enabled = enab(k))
        end do

        return

    end subroutine aeros_correction_load

    integer function aeros_correction_add_term(cor, sht, name, on_vor, on_div, &
                                                on_temp, on_lnps, lcut, enabled) &
                                                result(it)
        ! Register a term and allocate the fields it acts on, zeroed.
        !
        ! Returns the term's index, which is how a caller then fills it. Zeroed
        ! is the correct initial value and not merely a safe one: a registered
        ! but unfilled term must be a no-op, so that "the framework is wired
        ! in" and "a correction is being applied" stay separate statements.

        implicit none

        type(aeros_correction_class), intent(inout) :: cor
        type(aeros_sht_class),        intent(in)    :: sht
        character(len=*),             intent(in)    :: name
        logical, intent(in) :: on_vor, on_div, on_temp, on_lnps
        integer, intent(in) :: lcut
        logical, intent(in), optional :: enabled

        integer :: lm

        if (cor%nterm >= NTERM_MAX) then
            write(io_unit_err,*) "aeros_correction_add_term:: error: more than NTERM_MAX = ", &
                                    NTERM_MAX, " terms"
            error stop 1
        end if

        cor%nterm = cor%nterm + 1
        it = cor%nterm

        associate (t => cor%term(it))
            t%name    = name
            t%enabled = .TRUE.
            if (present(enabled)) t%enabled = enabled

            t%on_vor  = on_vor
            t%on_div  = on_div
            t%on_temp = on_temp
            t%on_lnps = on_lnps
            t%lcut    = lcut

            allocate(t%mask(cor%nlm))
            do lm = 1, cor%nlm
                t%mask(lm) = (lcut < 0) .or. (sht%l_of_lm(lm) <= lcut)
            end do

            if (on_vor)  then
                allocate(t%vor(cor%nlm, cor%nlev));  t%vor  = (0.0_wp_sh, 0.0_wp_sh)
            end if
            if (on_div)  then
                allocate(t%div(cor%nlm, cor%nlev));  t%div  = (0.0_wp_sh, 0.0_wp_sh)
            end if
            if (on_temp) then
                allocate(t%temp(cor%nlm, cor%nlev)); t%temp = (0.0_wp_sh, 0.0_wp_sh)
            end if
            if (on_lnps) then
                allocate(t%lnps(cor%nlm));           t%lnps = (0.0_wp_sh, 0.0_wp_sh)
            end if
        end associate

        return

    end function aeros_correction_add_term

    subroutine aeros_correction_apply(cor, tnd)
        ! DeltaF_total = sum_k DeltaF_k, added to the assembled tendency.
        !
        ! Called from aeros_timestep_step immediately after
        ! aeros_tendency_spectral. When `enabled` is .FALSE. this returns
        ! without touching anything, which is what makes the twin bit-exact
        ! rather than merely close.

        implicit none

        type(aeros_correction_class), intent(in)    :: cor
        type(aeros_tend_class),       intent(inout) :: tnd

        integer :: it, k, lm

        if (.not. cor%enabled) return

        do it = 1, cor%nterm
            associate (t => cor%term(it))
                if (.not. t%enabled) cycle

                if (t%on_vor .or. t%on_div .or. t%on_temp) then
                    do k = 1, cor%nlev
                        do lm = 1, cor%nlm
                            if (.not. t%mask(lm)) cycle
                            if (t%on_vor)  tnd%vor(lm,k)  = tnd%vor(lm,k)  + t%vor(lm,k)
                            if (t%on_div)  tnd%div(lm,k)  = tnd%div(lm,k)  + t%div(lm,k)
                            if (t%on_temp) tnd%temp(lm,k) = tnd%temp(lm,k) + t%temp(lm,k)
                        end do
                    end do
                end if

                if (t%on_lnps) then
                    do lm = 1, cor%nlm
                        if (.not. t%mask(lm)) cycle
                        tnd%lnps(lm) = tnd%lnps(lm) + t%lnps(lm)
                    end do
                end if
            end associate
        end do

        return

    end subroutine aeros_correction_apply

    subroutine aeros_correction_report(cor, io_unit)
        ! What each term is and how hard it is pulling.
        !
        ! The lnps (0,0) column is the one to watch for conservation: it is the
        ! only spectral mode that moves the global mean of ln(p_s), so it is
        ! the channel through which a correction injects mass. Section 3.7 risk
        ! 2 requires that channel to be zero (or explicitly accounted) before a
        ! correction is trusted over 10^5 yr. See the module header for what
        ! this diagnostic does and does not establish.

        implicit none

        type(aeros_correction_class), intent(in) :: cor
        integer, intent(in), optional :: io_unit

        integer  :: iou, it, nact
        real(dp) :: rms_vor, rms_div, rms_temp, rms_lnps, m00

        iou = 6
        if (present(io_unit)) iou = io_unit

        write(iou,*) ""
        write(iou,"(a)") " == tendency correction (design.md 3.7, 3.8) =="

        if (.not. cor%enabled) then
            write(iou,"(a)") "   correction                  OFF  (the DeltaF = 0 twin)"
        else
            write(iou,"(a)") "   correction                  ON"
        end if
        write(iou,"(a,i9)") "   terms registered           ", cor%nterm

        do it = 1, cor%nterm
            associate (t => cor%term(it))
                write(iou,*) ""
                write(iou,"(a,i0,a,a)")  "   term ", it, ": ", trim(t%name)
                if (t%enabled) then
                    write(iou,"(a)")     "     state                     enabled"
                else
                    write(iou,"(a)")     "     state                     disabled"
                end if
                write(iou,"(a,4l3)")     "     acts on (vor div T lnps) ", &
                                            t%on_vor, t%on_div, t%on_temp, t%on_lnps
                if (t%lcut < 0) then
                    write(iou,"(a)")     "     scales                    all resolved"
                else
                    write(iou,"(a,i0)")  "     scales                    l <= ", t%lcut
                    write(iou,"(a,i0,a,i0)") "     coefficients passing      ", &
                                            count(t%mask), " of ", cor%nlm
                end if

                call term_norms(t, cor%nlm, cor%nlev, rms_vor, rms_div, rms_temp, &
                                    rms_lnps, m00, nact)

                write(iou,"(a,i0)")      "     coefficients non-zero     ", nact
                if (t%on_vor)  write(iou,"(a,es12.3,a)") "     rms d(zeta)/dt           ", &
                                                            rms_vor,  " s-2"
                if (t%on_div)  write(iou,"(a,es12.3,a)") "     rms d(D)/dt              ", &
                                                            rms_div,  " s-2"
                if (t%on_temp) write(iou,"(a,es12.3,a)") "     rms d(T)/dt              ", &
                                                            rms_temp, " K s-1"
                if (t%on_lnps) then
                    write(iou,"(a,es12.3,a)") "     rms d(lnps)/dt           ", rms_lnps, " s-1"
                    write(iou,"(a,es12.3,a)") "     (0,0) part -- MASS INJECTION", m00, " s-1"
                end if
            end associate
        end do

        return

    end subroutine aeros_correction_report

    subroutine term_norms(t, nlm, nlev, rms_vor, rms_div, rms_temp, rms_lnps, m00, nact)
        ! RMS of each active field over the coefficients the mask passes, the
        ! (0,0) part of the lnps correction, and how many coefficients are
        ! actually non-zero.
        !
        ! `nact` is the honest answer to "is anything being applied at all",
        ! which for a no-op term must be zero and is the cheapest possible
        ! statement that the framework is wired in but idle.

        implicit none

        type(aeros_corr_term_class), intent(in)  :: t
        integer,                     intent(in)  :: nlm, nlev
        real(dp),                    intent(out) :: rms_vor, rms_div, rms_temp, rms_lnps, m00
        integer,                     intent(out) :: nact

        integer  :: k, lm, n
        real(dp) :: s_vor, s_div, s_temp, s_lnps

        s_vor = 0.0_dp; s_div = 0.0_dp; s_temp = 0.0_dp; s_lnps = 0.0_dp
        nact  = 0
        n     = 0
        m00   = 0.0_dp

        do k = 1, nlev
            do lm = 1, nlm
                if (.not. t%mask(lm)) cycle
                n = n + 1
                if (t%on_vor) then
                    s_vor = s_vor + real(abs(t%vor(lm,k)), dp)**2
                    if (t%vor(lm,k) /= (0.0_wp_sh, 0.0_wp_sh)) nact = nact + 1
                end if
                if (t%on_div) then
                    s_div = s_div + real(abs(t%div(lm,k)), dp)**2
                    if (t%div(lm,k) /= (0.0_wp_sh, 0.0_wp_sh)) nact = nact + 1
                end if
                if (t%on_temp) then
                    s_temp = s_temp + real(abs(t%temp(lm,k)), dp)**2
                    if (t%temp(lm,k) /= (0.0_wp_sh, 0.0_wp_sh)) nact = nact + 1
                end if
            end do
        end do

        if (t%on_lnps) then
            do lm = 1, nlm
                if (.not. t%mask(lm)) cycle
                s_lnps = s_lnps + real(abs(t%lnps(lm)), dp)**2
                if (t%lnps(lm) /= (0.0_wp_sh, 0.0_wp_sh)) nact = nact + 1
            end do
            ! Coefficient 1 is (l=0,m=0) -- the global mean of the correction.
            m00 = real(t%lnps(1), dp)
        end if

        rms_vor  = 0.0_dp; rms_div = 0.0_dp; rms_temp = 0.0_dp; rms_lnps = 0.0_dp
        if (n > 0) then
            rms_vor  = sqrt(s_vor /real(n, dp))
            rms_div  = sqrt(s_div /real(n, dp))
            rms_temp = sqrt(s_temp/real(n, dp))
        end if
        if (count(t%mask) > 0) rms_lnps = sqrt(s_lnps/real(count(t%mask), dp))

        return

    end subroutine term_norms

    subroutine aeros_correction_end(cor)

        implicit none

        type(aeros_correction_class), intent(inout) :: cor

        integer :: it

        do it = 1, NTERM_MAX
            associate (t => cor%term(it))
                if (allocated(t%mask)) deallocate(t%mask)
                if (allocated(t%vor))  deallocate(t%vor)
                if (allocated(t%div))  deallocate(t%div)
                if (allocated(t%temp)) deallocate(t%temp)
                if (allocated(t%lnps)) deallocate(t%lnps)
                t%name    = ""
                t%enabled = .FALSE.
                t%on_vor  = .FALSE.; t%on_div  = .FALSE.
                t%on_temp = .FALSE.; t%on_lnps = .FALSE.
                t%lcut    = -1
            end associate
        end do

        cor%enabled = .FALSE.
        cor%nterm   = 0
        cor%nlm     = 0
        cor%nlev    = 0

        return

    end subroutine aeros_correction_end

end module aeros_correction
