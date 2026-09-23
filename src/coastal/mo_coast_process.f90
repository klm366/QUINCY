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
!> multiply gross assimilation, as there.
!>
!> Zhao et al. assume that mangroves grow at mean sea level, so that the tide height
!> relative to mean sea level is also the height of the water above the ground. That
!> removes the need for a mangrove surface elevation, for which no global data set
!> exists, at the cost of not distinguishing high-intertidal from low-intertidal stands.
!>
!> Two parts differ from the ORCHIDEE-MAN implementation and are stated here:
!>  - the tide height is reconstructed inside the model from the harmonic constants of
!>    eight constituents instead of being predicted offline with OTIS. Nodal
!>    corrections and astronomical arguments are not applied, see calc_tide_height;
!>  - the inundated fraction of belowground root (Equation 2) uses the JSBACH root
!>    distribution, which is uniform over the rooting depth, instead of the exponential
!>    profile of ORCHIDEE.
!>
MODULE mo_coast_process
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_salinity_control, calc_tide_height, calc_bgr_inundated
  PUBLIC :: calc_agr_height_stilt, calc_agr_height_pneum
  PUBLIC :: calc_bgr_ventilated, calc_inundation_control, calc_coastal_factor
  PUBLIC :: cs_quadratic_maximum, ntide, tide_speed, days_since_2000, mangrove_lct_id

  !> lctlib id of the mangrove land cover type. In the npft=13 usecase the tiles pft01 to
  !> pft13 are created with the ids 1,2,3,4,9,10,11,12,14,15,19,20,21
  !> (mo_jsb_model_usecases), so the mangrove clone is 21. The coastal constraints are
  !> specific to mangroves and are evaluated on this land cover type only.
  INTEGER, PARAMETER :: mangrove_lct_id = 21

  ! Coefficients of Zhao et al. (2025), Equation 1, fitted to 26 observations of
  ! mangrove gross primary productivity against the long-term mean (2000 to 2014)
  ! surface salinity of Bio-ORACLE v2.0.
  REAL(wp), PARAMETER :: cs_a2 = -0.0214_wp
  REAL(wp), PARAMETER :: cs_a1 =  1.28_wp
  REAL(wp), PARAMETER :: cs_a0 = -18.3_wp

  !> Maximum of the published quadratic, reached at 29.907 PSS.
  !> The paper describes the maximum of the factor as one. The published coefficients
  !> give 0.8402; a maximum of one would need a constant of -18.14 instead of -18.3.
  REAL(wp), PARAMETER :: cs_quadratic_maximum = cs_a0 - cs_a1 * cs_a1 / (4._wp * cs_a2)

  !> Number of tidal constituents, in the order M2, S2, N2, K2, K1, O1, P1, Q1
  INTEGER,  PARAMETER :: ntide = 8

  !> Angular speeds of the constituents [degree per hour]
  REAL(wp), PARAMETER :: tide_speed(ntide) = [ &
    & 28.9841042_wp, 30.0000000_wp, 28.4397295_wp, 30.0821373_wp, &
    & 15.0410686_wp, 13.9430356_wp, 14.9589314_wp, 13.3986609_wp ]

  ! Maximum heights of the two types of aboveground root, Zhao et al. (2025), Table 1
  REAL(wp), PARAMETER :: h_agrs_max = 2.0_wp    !< Stilt roots [m]
  REAL(wp), PARAMETER :: h_agrp_max = 0.4_wp    !< Pneumatophores [m]

  REAL(wp), PARAMETER :: pi_local = 3.14159265358979323846_wp
  REAL(wp), PARAMETER :: deg2rad  = pi_local / 180._wp

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
  !> mangrove productivity. S is the long-term mean surface salinity in PSS, the
  !> predictor the coefficients were fitted against.
  !>
  !> The factor exceeds its minimum only between 25.92 and 33.89 PSS. Outside that
  !> interval it stays at Cs_min, so the result depends strongly on how the salinity
  !> field is sampled; open-ocean values leave most mangrove cells at the minimum.
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
  !>#### Compute the tide height relative to mean sea level
  !>
  !> We here reconstruct the sea surface elevation from the complex harmonic constants of
  !> the TPXO9-atlas solution:
  !>
  !>     eta(t) = sum_c [ hre_c * cos(omega_c * t) - him_c * sin(omega_c * t) ]
  !>
  !> which follows from the TPXO convention amp = |hre + i him|, G = atan2(-him, hre)
  !> and eta = sum amp * cos(omega * t - G).
  !>
  !> Nodal corrections and astronomical arguments are not applied, so the time t is
  !> counted from an arbitrary reference. For the constituents, whose frequencies are
  !> not commensurate, this leaves the long-term distribution of the tide height and
  !> therefore the frequency and duration of inundation unchanged, but shifts the tide
  !> in calendar time. Nodal modulation would change the M2 amplitude by up to about
  !> 3.7 per cent and the K1 amplitude by up to about 11 per cent over 18.6 years.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  PURE FUNCTION calc_tide_height(hre, him, time_hours) RESULT(eta)

    REAL(wp), INTENT(in) :: &
      & hre(ntide),         & !< Real part of the complex amplitude of each constituent [m]
      & him(ntide),         & !< Imaginary part of the complex amplitude [m]
      & time_hours            !< Time since the reference [h]

    REAL(wp) :: eta           !< RESULT: tide height above mean sea level [m]

    INTEGER  :: ic
    REAL(wp) :: arg

    eta = 0._wp
    DO ic = 1, ntide
      arg = MODULO(tide_speed(ic) * time_hours, 360._wp) * deg2rad
      eta = eta + hre(ic) * COS(arg) - him(ic) * SIN(arg)
    END DO

  END FUNCTION calc_tide_height

  ! ======================================================================================================== !
  !>
  !>#### Compute the inundated fraction of belowground root
  !>
  !> Zhao et al. (2025), Equation 2, with the ground at mean sea level. The fraction is
  !> the part of the root distribution that lies below the water level. Three cases:
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
  !> The base of the logarithm is not stated in the publication. The decadic logarithm
  !> is used because the expression is a power law of the form 10**(...) and because it
  !> gives 1.095 m for 4 Mg ha-1, a plausible stilt root height; the natural logarithm
  !> would give 1.71 m.
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
  !> Zhao et al. (2025), Equation 6. The factor was derived from a mean pneumatophore
  !> height of about 0.3 m at a mean biomass of 4.55 Mg ha-1.
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
  !> submerged part of an aerial root does not ventilate. The biomass ratio limits the
  !> ventilating capacity of young or sparse stands. Each term is limited to the range
  !> from zero to its weight, which covers tide heights below the ground as well as
  !> above the tip of the aerial roots.
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
  !>#### Combine the two control factors
  !>
  !> Both factors multiply gross assimilation. With both application switches false the
  !> result is exactly one.
  !>
  !> The factors are applied independently. Tidal flushing removes accumulated salt, so
  !> a dependence of the salinity factor on the inundation regime would be physically
  !> reasonable; none is described in Zhao et al. (2025) and none is applied here. With
  !> both factors at their minimum of 0.5 the combined factor is 0.25.
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

  ! ======================================================================================================== !
  !>
  !>#### Number of days from 1 January 2000 to a date of the proleptic Gregorian calendar
  !>
  !> We here count days with the Julian day number of Fliegel and Van Flandern (1968), so
  !> that the tide can be evaluated at the model time. With a 365-day calendar the count
  !> would slip by one day at each leap day; the tide statistics are not affected.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION days_since_2000(year, month, day) RESULT(ndays)

    INTEGER, INTENT(in) :: year, month, day
    INTEGER :: ndays
    INTEGER :: a, jdn

    a   = (month - 14) / 12
    jdn = (1461 * (year + 4800 + a)) / 4 + (367 * (month - 2 - 12 * a)) / 12 &
      &   - (3 * ((year + 4900 + a) / 100)) / 4 + day - 32075
    ndays = jdn - 2451545

  END FUNCTION days_since_2000

#endif
END MODULE mo_coast_process
