! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Scientific calculations for the coastal-constraint process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>#### Salinity and tidal-inundation control factors on mangrove productivity
!>
!> The routines below implement the coastal constraints of Zhao et al. (2025),
!> ORCHIDEE-MAN, J. Adv. Model. Earth Syst., doi:10.1029/2025MS005185, Sections 2.2.2
!> and 2.2.3. Equation numbers refer to that paper. The two control factors are combined
!> into one scalar factor that multiplies gross assimilation, as there.
!>
!> Zhao et al. assume that mangroves grow at mean sea level, so that the tide height
!> relative to mean sea level is also the height of the water above the ground.
!>
!> Two parts differ from the ORCHIDEE-MAN implementation:
!>  - the tide height is reconstructed inside the model from the harmonic constants of
!>    eight constituents instead of being predicted offline with OTIS; nodal corrections
!>    and astronomical arguments are not applied, see calc_tide_arguments;
!>  - the inundated fraction of belowground root (Equation 2) uses the JSBACH root
!>    distribution, which is uniform over the rooting depth, instead of the exponential
!>    profile of ORCHIDEE.
!>
MODULE mo_coast_process
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_jsb_math_constants, ONLY: pi, deg2rad
  ! cs_salinity_optimum imported for the one-sided salinity option.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-09-23
  USE mo_coast_constants,    ONLY: cs_a2, cs_a1, cs_a0, cs_quadratic_maximum, cs_salinity_optimum, &
    &                              ntide, tide_speed,                                                &
    &                              h_agrs_max, h_agrp_max, h_agrp_per_biomass

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_salinity_control, calc_tide_arguments, calc_tide_height, calc_bgr_inundated
  PUBLIC :: calc_agr_height_stilt, calc_agr_height_pneum
  PUBLIC :: calc_bgr_ventilated, calc_inundation_control, calc_coastal_factor

  CHARACTER(len=*), PARAMETER :: modname = 'mo_coast_process'

CONTAINS

  ! ======================================================================================================== !
  !>
  !>#### Compute the control factor of salinity on productivity
  !>
  !> Zhao et al. (2025), Equation 1:
  !>
  !>     Cs = max( -0.0214 * S**2 + 1.28 * S - 18.3 , Cs_min )
  !>
  !> The relationship is an inverted parabola: both high and low salinity reduce mangrove
  !> productivity. S is the long-term mean surface salinity in PSS, the predictor the
  !> coefficients were fitted against. With Cs_min = 0.5 and without normalisation the
  !> factor exceeds its minimum only between 25.92 and 33.89 PSS.
  !>
  !> With l_one_sided true, a salinity below the optimum of the quadratic is replaced by
  !> the optimum, so only salinity above it reduces the factor. This is a sensitivity
  !> option and not part of the publication.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ! Argument l_one_sided added for the sensitivity switch l_cs_one_sided.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-09-23
  ELEMENTAL PURE FUNCTION calc_salinity_control(salinity, cs_min, l_normalise, l_one_sided) RESULT(cs)

    REAL(wp), INTENT(in) :: &
      & salinity,           & !< Long-term mean sea surface salinity [PSS]
      & cs_min                !< Minimum value of the control factor []
    LOGICAL,  INTENT(in) :: &
      & l_normalise,        & !< Scale so that the maximum of the factor is one
      & l_one_sided           !< No reduction below the optimum salinity

    REAL(wp) :: cs            !< RESULT: control factor of salinity []
    REAL(wp) :: s_eff         !< Salinity entering the quadratic [PSS]

    s_eff = salinity
    IF (l_one_sided) s_eff = MAX(salinity, cs_salinity_optimum)
    cs = cs_a2 * s_eff * s_eff + cs_a1 * s_eff + cs_a0
    IF (l_normalise) cs = cs / cs_quadratic_maximum
    cs = MIN(1._wp, MAX(cs, cs_min))

  END FUNCTION calc_salinity_control

  ! ======================================================================================================== !
  !>
  !>#### Compute the phase terms of the tidal constituents at the model time
  !>
  !> We here evaluate cos(omega_c * t) and sin(omega_c * t) for the eight constituents.
  !> They depend on time only, not on the grid cell, so they are computed once per call
  !> of the update routine and reused for every cell.
  !>
  !> Nodal corrections and astronomical arguments are not applied, so t is counted from
  !> an arbitrary reference. The frequencies of the constituents are not commensurate, so
  !> this leaves the long-term distribution of the tide height, and therefore the
  !> frequency and duration of inundation, unchanged, but shifts the tide in calendar
  !> time. Over the 18.6-year nodal cycle the nodal modulation would change the M2
  !> amplitude by up to about 3.7 per cent and the K1 amplitude by up to about 11 per cent.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  PURE SUBROUTINE calc_tide_arguments(time_hours, cos_arg, sin_arg)

    REAL(wp), INTENT(in)  :: time_hours          !< Time since the reference [h]
    REAL(wp), INTENT(out) :: cos_arg(ntide)      !< cos(omega_c * t) []
    REAL(wp), INTENT(out) :: sin_arg(ntide)      !< sin(omega_c * t) []

    INTEGER  :: itide
    REAL(wp) :: arg

    DO itide = 1, ntide
      arg = MODULO(tide_speed(itide) * time_hours, 360._wp) * deg2rad
      cos_arg(itide) = COS(arg)
      sin_arg(itide) = SIN(arg)
    END DO

  END SUBROUTINE calc_tide_arguments

  ! ======================================================================================================== !
  !>
  !>#### Compute the tide height relative to mean sea level
  !>
  !> We here sum the constituents of the TPXO9-atlas solution:
  !>
  !>     eta = sum_c [ hre_c * cos(omega_c * t) - him_c * sin(omega_c * t) ]
  !>
  !> which follows from the TPXO convention amp = |hre + i him|, G = atan2(-him, hre) and
  !> eta = sum amp * cos(omega * t - G).
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  PURE FUNCTION calc_tide_height(hre, him, cos_arg, sin_arg) RESULT(eta)

    REAL(wp), INTENT(in) :: &
      & hre(ntide),         & !< Real part of the complex amplitude of each constituent [m]
      & him(ntide),         & !< Imaginary part of the complex amplitude [m]
      & cos_arg(ntide),     & !< cos(omega_c * t) from calc_tide_arguments []
      & sin_arg(ntide)        !< sin(omega_c * t) from calc_tide_arguments []

    REAL(wp) :: eta           !< RESULT: tide height above mean sea level [m]

    INTEGER :: itide

    eta = 0._wp
    DO itide = 1, ntide
      eta = eta + hre(itide) * cos_arg(itide) - him(itide) * sin_arg(itide)
    END DO

  END FUNCTION calc_tide_height

  ! ======================================================================================================== !
  !>
  !>#### Compute the inundated fraction of belowground root
  !>
  !> Zhao et al. (2025), Equation 2, with the ground at mean sea level. The fraction is the
  !> part of the root distribution that lies below the water level. Three cases:
  !>
  !>   tide below the rooting depth       no root inundated
  !>   tide between rooting depth and 0   the roots below the water level are inundated
  !>   tide above the ground              all roots inundated
  !>
  !> Roots are distributed uniformly over the rooting depth in JSBACH, so the sum over
  !> soil layers reduces to (root_depth + eta) / root_depth in the middle case.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_bgr_inundated(tide_height, root_depth) RESULT(f_inun)

    REAL(wp), INTENT(in) :: &
      & tide_height,        & !< Tide height above mean sea level, i.e. above ground [m]
      & root_depth            !< Rooting depth [m]

    REAL(wp) :: f_inun        !< RESULT: inundated fraction of belowground root []

    IF (root_depth <= 0._wp) THEN
      f_inun = 0._wp
    ELSE IF (tide_height >= 0._wp) THEN
      f_inun = 1._wp
    ELSE IF (tide_height <= -root_depth) THEN
      f_inun = 0._wp
    ELSE
      f_inun = (root_depth + tide_height) / root_depth
    END IF

  END FUNCTION calc_bgr_inundated

  ! ======================================================================================================== !
  !>
  !>#### Compute the height of stilt roots from their biomass
  !>
  !> Zhao et al. (2025), Equation 5, after Ohira et al. (2013):
  !>
  !>     H = 10**((log10(B * 1000 / D) + 2.945) / 2.546) * 7.56 / (pi * 100) + 0.5
  !>
  !> The base of the logarithm is not stated in the publication. The decadic logarithm is
  !> used because the expression is a power law of the form 10**(...) and because it gives
  !> 1.095 m for 4 Mg ha-1, a plausible stilt-root height; the natural logarithm would give
  !> 1.71 m.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_agr_height_stilt(b_agrs, tree_density) RESULT(h_agrs)

    REAL(wp), INTENT(in) :: &
      & b_agrs,             & !< Stilt-root biomass [Mg ha-1]
      & tree_density          !< Stand density [trees ha-1]

    REAL(wp) :: h_agrs        !< RESULT: stilt-root height [m]
    REAL(wp) :: biomass_per_tree

    IF (b_agrs <= 0._wp .OR. tree_density <= 0._wp) THEN
      h_agrs = 0._wp
      RETURN
    END IF

    biomass_per_tree = b_agrs * 1000._wp / tree_density
    h_agrs = 10._wp**((LOG10(biomass_per_tree) + 2.945_wp) / 2.546_wp) &
      &      * 7.56_wp / (pi * 100._wp) + 0.5_wp
    h_agrs = MIN(h_agrs, h_agrs_max)

  END FUNCTION calc_agr_height_stilt

  ! ======================================================================================================== !
  !>
  !>#### Compute the height of pneumatophores from their biomass
  !>
  !> Zhao et al. (2025), Equation 6.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_agr_height_pneum(b_agrp) RESULT(h_agrp)

    REAL(wp), INTENT(in) :: b_agrp   !< Pneumatophore biomass [Mg ha-1]
    REAL(wp) :: h_agrp               !< RESULT: pneumatophore height [m]

    h_agrp = MAX(0._wp, b_agrp) * h_agrp_per_biomass
    h_agrp = MIN(h_agrp, h_agrp_max)

  END FUNCTION calc_agr_height_pneum

  ! ======================================================================================================== !
  !>
  !>#### Compute the fraction of belowground root ventilated by aboveground roots
  !>
  !> Zhao et al. (2025), Equation 3:
  !>
  !>     F_ven = a * (B_s / B_s_lim) * (H_s - H_tide) / H_s
  !>           + b * (B_p / B_p_lim) * (H_p - H_tide) / H_p
  !>
  !> Lenticels are assumed to be distributed uniformly along the aerial roots, so the
  !> submerged part of an aerial root does not ventilate. The biomass ratio limits the
  !> ventilating capacity of young or sparse stands. Each term is limited to the range
  !> from zero to its weight, which covers tide heights below the ground as well as above
  !> the tip of the aerial roots.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_bgr_ventilated( &
    & b_agrs, b_agrp, h_tide, b_agrs_lim, b_agrp_lim, weight_s, weight_p, tree_density) &
    & RESULT(f_ven)

    REAL(wp), INTENT(in) :: &
      & b_agrs,             & !< Stilt-root biomass [Mg ha-1]
      & b_agrp,             & !< Pneumatophore biomass [Mg ha-1]
      & h_tide,             & !< Tide height above the ground [m]
      & b_agrs_lim,         & !< Stilt-root biomass giving full ventilation [Mg ha-1]
      & b_agrp_lim,         & !< Pneumatophore biomass giving full ventilation [Mg ha-1]
      & weight_s,           & !< Weight of the stilt-root term []
      & weight_p,           & !< Weight of the pneumatophore term []
      & tree_density          !< Stand density [trees ha-1]

    REAL(wp) :: f_ven         !< RESULT: ventilated fraction of belowground root []

    ! Local variables
    REAL(wp) ::             &
      & h_s,                & !< Stilt-root height [m]
      & h_p,                & !< Pneumatophore height [m]
      & term_s,             & !< Contribution of the stilt roots []
      & term_p                !< Contribution of the pneumatophores []
    ! ------------------------------------------------------------------------------------------------------ !

    h_s = calc_agr_height_stilt(b_agrs, tree_density)
    h_p = calc_agr_height_pneum(b_agrp)

    term_s = 0._wp
    IF (h_s > 0._wp .AND. b_agrs_lim > 0._wp) &
      & term_s = weight_s * MIN(1._wp, b_agrs / b_agrs_lim) &
      &          * MIN(1._wp, MAX(0._wp, (h_s - h_tide) / h_s))

    term_p = 0._wp
    IF (h_p > 0._wp .AND. b_agrp_lim > 0._wp) &
      & term_p = weight_p * MIN(1._wp, b_agrp / b_agrp_lim) &
      &          * MIN(1._wp, MAX(0._wp, (h_p - h_tide) / h_p))

    f_ven = MIN(1._wp, MAX(0._wp, term_s + term_p))

  END FUNCTION calc_bgr_ventilated

  ! ======================================================================================================== !
  !>
  !>#### Compute the control factor of tidal inundation on productivity
  !>
  !> Zhao et al. (2025), Equations 4 and 7:
  !>
  !>     F_ano = F_inun - F_ven
  !>     Cin   = max( 1 - (1 - Cin_min) * F_ano , Cin_min )
  !>
  !> The anoxic fraction of belowground root is the inundated fraction reduced by the
  !> fraction the aerial roots can still ventilate, limited to the range from 0 to 1.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_inundation_control(f_inun, f_ven, cin_min) RESULT(cin)

    REAL(wp), INTENT(in) :: &
      & f_inun,             & !< Inundated fraction of belowground root []
      & f_ven,              & !< Ventilated fraction of belowground root []
      & cin_min               !< Minimum value of the control factor []

    REAL(wp) :: cin           !< RESULT: control factor of tidal inundation []
    REAL(wp) :: f_ano         !< Anoxic fraction of belowground root []

    f_ano = MIN(1._wp, MAX(0._wp, f_inun - f_ven))
    cin   = MAX(1._wp - (1._wp - cin_min) * f_ano, cin_min)

  END FUNCTION calc_inundation_control

  ! ======================================================================================================== !
  !>
  !>#### Combine the two control factors into one scalar coastal stress factor
  !>
  !>     f_coast = Cs * Cin
  !>
  !> Both factors multiply gross assimilation, as in Zhao et al. (2025). A switch that is
  !> false replaces its factor by one, so with both switches false the result is exactly
  !> one. The two factors are treated as independent; with both at their minimum of 0.5
  !> the combined factor is 0.25.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_coastal_factor(cs, cin, l_apply_salinity, l_apply_tide) RESULT(f_coast)

    REAL(wp), INTENT(in) :: &
      & cs,                 & !< Control factor of salinity []
      & cin                   !< Control factor of tidal inundation []
    LOGICAL,  INTENT(in) :: &
      & l_apply_salinity,   & !< Apply the salinity factor
      & l_apply_tide          !< Apply the tidal-inundation factor

    REAL(wp) :: f_coast       !< RESULT: combined factor on gross assimilation []

    f_coast = 1._wp
    IF (l_apply_salinity) f_coast = f_coast * cs
    IF (l_apply_tide)     f_coast = f_coast * cin

  END FUNCTION calc_coastal_factor

#endif
END MODULE mo_coast_process
