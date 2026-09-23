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
!> ORCHIDEE-MAN, J. Adv. Model. Earth Syst., doi:10.1029/2025MS005185, Sections
!> 2.2.2 and 2.2.3. Equation numbers refer to that paper. Both control factors
!> multiply gross assimilation, which is how they are applied there.
!>
!> Zhao et al. assume that mangroves grow at mean sea level, so that the tide height
!> is also the inundation depth. That assumption removes the need for a mangrove
!> surface elevation, for which no global data set exists, at the cost of not being
!> able to distinguish high-intertidal from low-intertidal stands.
!>
MODULE mo_coast_process
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_salinity_control, calc_inundation_control
  PUBLIC :: calc_agr_height_stilt, calc_agr_height_pneum
  PUBLIC :: calc_bgr_ventilated, calc_coastal_factor
  PUBLIC :: cs_quadratic_maximum

  ! Coefficients of Zhao et al. (2025), Equation 1, fitted to 26 observations of
  ! mangrove gross primary productivity against the long-term mean (2000 to 2014)
  ! surface salinity of Bio-ORACLE v2.0.
  REAL(wp), PARAMETER :: cs_a2 = -0.0214_wp
  REAL(wp), PARAMETER :: cs_a1 =  1.28_wp
  REAL(wp), PARAMETER :: cs_a0 = -18.3_wp

  !> Maximum of the published quadratic, reached at 29.907 PSS.
  !> The paper describes the maximum of the factor as one. The published
  !> coefficients give 0.8402; a maximum of one would require a constant of -18.14
  !> rather than -18.3. The value is provided so that the normalisation offered by
  !> the namelist switch l_normalise_cs is exact.
  REAL(wp), PARAMETER :: cs_quadratic_maximum = cs_a0 - cs_a1 * cs_a1 / (4._wp * cs_a2)

  ! Maximum heights of the two types of aboveground root, Zhao et al. (2025), Table 1
  REAL(wp), PARAMETER :: h_agrs_max = 2.0_wp    !< Stilt roots [m]
  REAL(wp), PARAMETER :: h_agrp_max = 0.4_wp    !< Pneumatophores [m]

  REAL(wp), PARAMETER :: pi_local = 3.14159265358979_wp

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
  !> The relationship is an inverted parabola: both high and low salinity reduce
  !> mangrove productivity. S is the long-term mean surface salinity in PSS, which is
  !> the predictor the coefficients were fitted against.
  !>
  !> Note that the factor exceeds its minimum only for salinities between 25.92 and
  !> 33.89 PSS. Outside that interval it is held at Cs_min and carries no spatial
  !> information, so the result depends strongly on how the salinity field is sampled.
  !> Sampling open ocean rather than coastal water leaves nearly all mangrove cells at
  !> the minimum.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_salinity_control(salinity, cs_min, l_normalise) RESULT(cs)

    REAL(wp), INTENT(in) :: &
      & salinity,           & !< Long-term mean sea surface salinity [PSS]
      & cs_min                !< Minimum value of the control factor []
    LOGICAL,  INTENT(in) :: &
      & l_normalise           !< Scale so that the maximum of the factor is one

    REAL(wp) :: cs            !< RESULT: control factor of salinity []

    cs = cs_a2 * salinity * salinity + cs_a1 * salinity + cs_a0
    IF (l_normalise) cs = cs / cs_quadratic_maximum
    cs = MIN(1._wp, MAX(cs, cs_min))

  END FUNCTION calc_salinity_control

  ! ======================================================================================================== !
  !>
  !>#### Compute the height of stilt roots from their biomass
  !>
  !> Zhao et al. (2025), Equation 5, after Ohira et al. (2013):
  !>
  !>     H = 10**((log10(B * 1000 / D) + 2.945) / 2.546) * 7.56 / (pi * 100) + 0.5
  !>
  !> The base of the logarithm is not stated in the publication. The decadic logarithm
  !> is used here because the expression is a power law of the form 10**(...), and
  !> because it gives 1.095 m for a biomass of 4 Mg ha-1, which is a plausible stilt
  !> root height; the natural logarithm would give 1.71 m.
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
      &      * 7.56_wp / (pi_local * 100._wp) + 0.5_wp
    h_agrs = MIN(h_agrs, h_agrs_max)

  END FUNCTION calc_agr_height_stilt

  ! ======================================================================================================== !
  !>
  !>#### Compute the height of pneumatophores from their biomass
  !>
  !> Zhao et al. (2025), Equation 6. The conversion factor was derived from a mean
  !> pneumatophore height of about 0.3 m at a mean biomass of 4.55 Mg ha-1.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_agr_height_pneum(b_agrp) RESULT(h_agrp)

    REAL(wp), INTENT(in) :: b_agrp   !< Pneumatophore biomass [Mg ha-1]
    REAL(wp) :: h_agrp               !< RESULT: pneumatophore height [m]

    h_agrp = MAX(0._wp, b_agrp) * 0.066_wp
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
  !> submerged part of an aerial root loses its ventilating function. The biomass ratio
  !> limits the ventilating capacity of young or sparse stands.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION calc_bgr_ventilated( &
    & b_agrs, b_agrp, h_tide, b_agrs_lim, b_agrp_lim, weight_s, weight_p, tree_density) &
    & RESULT(f_ven)

    REAL(wp), INTENT(in) :: &
      & b_agrs,             & !< Stilt-root biomass [Mg ha-1]
      & b_agrp,             & !< Pneumatophore biomass [Mg ha-1]
      & h_tide,             & !< Tide height above mean sea level [m]
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
      & term_s = weight_s * MIN(1._wp, b_agrs / b_agrs_lim) * MAX(0._wp, (h_s - h_tide) / h_s)

    term_p = 0._wp
    IF (h_p > 0._wp .AND. b_agrp_lim > 0._wp) &
      & term_p = weight_p * MIN(1._wp, b_agrp / b_agrp_lim) * MAX(0._wp, (h_p - h_tide) / h_p)

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
  !> fraction that the aerial roots are still able to ventilate. All three fractions
  !> are limited to the interval from zero to one.
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
  !>#### Combine the two control factors
  !>
  !> Both factors multiply gross assimilation. With both application switches false the
  !> result is exactly one, so a simulation without the coastal constraints reproduces
  !> the previous behaviour.
  !>
  !> The two factors are applied independently. Tidal flushing removes accumulated
  !> salt, so a dependence of the salinity factor on the inundation regime would be
  !> physically reasonable; no such coupling is described in Zhao et al. (2025) and
  !> none is applied here. Note also that with both factors at their minimum of 0.5
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
