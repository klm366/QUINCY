! Author: Dr. Sri Nandini-Weiß
! Date: 2026-07-15
! Generalised from PFT13 to all plant functional types
! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Scientific calculations for the trait-variability process
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
!>#### Calculates diagnostic traits from completed previous-year climate
!>
!> The trait equations were previously scalar parameters valid for PFT13 only. They are
!> now arrays indexed by the lctlib id of the tile's primary land cover type, and one
!> equation form is evaluated for every plant functional type that has a fitted model:
!>
!>     trait_raw = a0 + a1 * MAT + a2 * MAP + a3 * Tmax + a4 * MAP * Tmax
!>
!> Coefficients are zero where a PFT does not use a term. PFT13 is the only PFT using
!> the interaction term a4. Its coefficients are unchanged from the validated PFT13
!> implementation and are evaluated in the same term order, so PFT13 results are
!> identical to the previous version.
!>
!> A PFT without a fitted model is not evaluated at all; it is skipped by the
!> has_dynamic_* mask and retains its static lctlib parameters. Evaluating an equation
!> with zero coefficients is not equivalent to a static fallback.
!>
MODULE mo_trait_process
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_physical_constants, ONLY: tmelt, molar_mass_C

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_trait_diagnostics
  PUBLIC :: has_dynamic_vcmax25, has_dynamic_jmax25, has_dynamic_sla
  PUBLIC :: trait_vcmax25_to_carboxrate
  PUBLIC :: trait_jmax25_to_etransport
  PUBLIC :: trait_sla_to_specificleafarea_c
  PUBLIC :: max_lct_id

  ! Coefficient tables are indexed by lctlib id rather than held as scalars.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-09-23
  !> Largest lctlib id that can carry a fitted trait model. The mangrove lctlib file has
  !> 22 columns and Read_lctlib discards the glacier column for the jsbach usecases, so
  !> model%lctlib runs from 1 to 21 and the mangrove clone is id 21.
  INTEGER, PARAMETER :: max_lct_id = 24

  ! ======================================================================================================== !
  ! Unit conversions for active trait use
  ! ======================================================================================================== !
  ! Vcmax25 and Jmax25 diagnostics are held in umol m-2 s-1, whereas CarboxRate and
  ! ETransport are mol m-2 s-1. init_lctlib applies the factor 1.e-6 once during model
  ! initialisation, so the consuming routines already work in mol m-2 s-1.
  REAL(wp), PARAMETER :: trait_vcmax25_to_carboxrate = 1.e-6_wp
  REAL(wp), PARAMETER :: trait_jmax25_to_etransport  = 1.e-6_wp

  ! The SLA conversion is now per plant functional type, because the fitted
  !   equations do not share one input unit.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-09-23
  ! lctlib%specificLeafArea_C is a specific leaf area in m2(leaf) mol(C)-1. The type
  ! declaration comment in mo_jsb_lctlib_class states the reciprocal unit and is
  ! incorrect; the lctlib data file header gives the correct one, and all three
  ! consuming routines multiply a carbon pool by it to obtain an area. A larger
  ! predicted SLA therefore has to give a larger specificLeafArea_C.
  !
  ! The fitted equations do not share one input unit. The PFT13 equation is in
  ! m2 gC-1 and needs only molar_mass_C. The PFT01 and PFT04 equations are in
  ! mm2 mg-1 leaf dry mass and need the additional dry-mass-to-carbon conversion
  ! below. Applying molar_mass_C to a mm2 mg-1 equation overestimates
  ! specificLeafArea_C by a factor of 480, so the conversion is carried per PFT
  ! rather than as a single module constant.
  !
  !   SLA[mm2 mg-1] * 1.e-3                  = SLA[m2 gDW-1]
  !                 / leaf_carbon_fraction   = SLA[m2 gC-1]
  !                 * molar_mass_C           = specificLeafArea_C[m2 molC-1]
  !
  ! The same expression is already applied to the QUINCY lctlib sla field during
  ! initialisation. carbon_per_dryweight_leaf in mo_veg_constants is assigned only from
  ! init_veg_constants, which is called for MODEL_QUINCY and not for MODEL_JSBACH, so
  ! the value is repeated here as a local parameter instead of being imported.
  REAL(wp), PARAMETER :: leaf_carbon_fraction = 0.48_wp
  REAL(wp), PARAMETER :: k_sla = 1.e-3_wp * molar_mass_C / leaf_carbon_fraction

  !> Factor converting the bounded SLA prediction to specificLeafArea_C [m2 molC-1]
  REAL(wp), PARAMETER :: trait_sla_to_specificleafarea_c(max_lct_id) = [ &
    & k_sla,        &  !  1 tropical evergreen trees, fitted in mm2 mg-1
    & 0._wp,        &  !  2
    & 0._wp,        &  !  3
    & k_sla,        &  !  4 extra-tropical deciduous trees, fitted in mm2 mg-1
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, & !  5 - 12
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, & ! 13 - 20
    & molar_mass_C, &  ! 21 mangrove clone, fitted in m2 gC-1
    & 0._wp, 0._wp, 0._wp ]                                     ! 22 - 24

  ! ======================================================================================================== !
  ! Fitted trait models
  ! ======================================================================================================== !
  ! Predictors are taken from the previous calendar year provided by CLIMBUF. MAT and
  ! the warmest monthly mean temperature are converted from K to degC below; annual
  ! precipitation in kg m-2 is numerically identical to mm yr-1 and needs no conversion.
  ! The offline warmest-month predictor is the maximum of the twelve monthly mean
  ! temperatures, which is what CLIMBUF accumulates in max_mmtemp_of_yr.
  !
  ! Coefficients for lctlib ids 1, 2 and 4 are taken from the frozen offline table
  ! step_d_frozen_coefficients_20260917.csv, rows with Activate = True. Coefficients for
  ! id 21 are unchanged from the validated PFT13 implementation.
  !
  !   id   trait     model     LOSO R2   fitted unit      state
  !    1   Vcmax25   Tmax      0.250     umol m-2 s-1     present, not activated
  !    2   Vcmax25   MAP       0.387     umol m-2 s-1     present, not activated
  !    1   SLA       Tmax      0.145     mm2 mg-1         present, not activated
  !    4   SLA       MAT_MAP   0.126     mm2 mg-1         activated
  !   21   Vcmax25   MAT+MAP   frozen    umol m-2 s-1     activated
  !   21   Jmax25    MAT+MAP   frozen    umol m-2 s-1     activated
  !   21   SLA       MAP,Tmax  frozen    m2 gC-1          activated
  !
  ! The Vcmax25 models for ids 1 and 2 are not activated because the frozen table
  ! reports a bound crossing inside their calibration range, so the bounds rather than
  ! the regression determine the value over part of the fitted domain. The SLA model for
  ! id 1 is not activated because its calibration window spans about 10 degC and the
  ! equation saturates outside it. The SLA model for id 4 has the widest calibration
  ! range and no bound crossing on the fitted units, but note that its precipitation
  ! slope is significant only at p = 0.089 and that the residuals are heteroscedastic
  ! (Breusch-Pagan p = 0.0027).
  !
  ! Evaluating the id 4 equation over the full rectangle spanned by its calibration
  ! ranges, rather than over the fitted units alone, gives values below the lower bound
  ! in about 5 % of cases, at the warm and wet corner. No observed site combines those
  ! two extremes, but a global grid will. Clamping frequency should be diagnosed per
  ! grid cell rather than assumed from the offline table.

  ! Coefficients for lctlib ids 1 and 2 added from the frozen offline table.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-09-23
  !--------------------------------------- Vcmax25 [umol m-2 s-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_vcmax25(max_lct_id) = [ &
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., & !  1 -  8
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., & !  9 - 16
    & .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE., .FALSE. ]  ! 17 - 24

  REAL(wp), PARAMETER :: vcmax25_a0(max_lct_id) = [ &
    & -123.55107783642907_wp, 83.56109481863965_wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 49.5690779208497_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: vcmax25_a1(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, -0.391997574087549_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: vcmax25_a2(max_lct_id) = [ &
    & 0._wp, -0.0186094309088575_wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0.001594603433354_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: vcmax25_a3(max_lct_id) = [ &
    & 6.289335485105433_wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: vcmax25_a4(max_lct_id) = 0._wp
  REAL(wp), PARAMETER :: vcmax25_lower(max_lct_id) = [ &
    & 22.963071432422563_wp, 20.62060675208818_wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 16.5700_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: vcmax25_upper(max_lct_id) = [ &
    & 92.9706204665844_wp, 83.72033981692316_wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 66.2500_wp, 0._wp, 0._wp, 0._wp ]

  !---------------------------------------- Jmax25 [umol m-2 s-1] ----------------------------------------!
  ! No PFT other than the mangrove clone has an activatable Jmax25 model. Every Jmax25
  ! row of the frozen table is classified as a static empirical value, as diagnostic
  ! only, or as insufficient data; the best leave-one-site-out R2 is 0.099, below the
  ! usable threshold of 0.1. Whether Jmax25 should instead follow Vcmax25 through a
  ! fixed or temperature-dependent Jmax to Vcmax ratio is not decided here and no such
  ! pathway is implemented.
  LOGICAL,  PARAMETER :: has_dynamic_jmax25(max_lct_id) = [ &
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
    & .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE., .FALSE. ]

  REAL(wp), PARAMETER :: jmax25_a0(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 180.411473978400_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: jmax25_a1(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, -3.09078405961564_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: jmax25_a2(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, -0.00485175158785102_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: jmax25_a3(max_lct_id) = 0._wp
  REAL(wp), PARAMETER :: jmax25_a4(max_lct_id) = 0._wp
  REAL(wp), PARAMETER :: jmax25_lower(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 46.3775_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: jmax25_upper(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 131.4725_wp, 0._wp, 0._wp, 0._wp ]

  ! Coefficients for lctlib ids 1 and 4 added from the frozen offline table.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-09-23
  !------------------------------------------- SLA [fitted unit] -----------------------------------------!
  ! The raw SLA diagnostic is held in the unit of each fitted equation, which is
  ! mm2 mg-1 for lctlib ids 1 and 4 and m2 gC-1 for id 21. Only the converted
  ! specificLeafArea_C is comparable between plant functional types.
  LOGICAL,  PARAMETER :: has_dynamic_sla(max_lct_id) = [ &
    & .FALSE., .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE., .FALSE., .FALSE., &
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
    & .FALSE., .FALSE., .FALSE., .FALSE., .TRUE.,  .FALSE., .FALSE., .FALSE. ]

  REAL(wp), PARAMETER :: sla_a0(max_lct_id) = [ &
    & 45.5576887921345_wp, 0._wp, 0._wp, 24.818121977359752_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0.060960327_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: sla_a1(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, -0.519495959904962_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: sla_a2(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, -0.0015950749148317_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & -0.000045172134_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: sla_a3(max_lct_id) = [ &
    & -1.270761413409342_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & -0.0015723334_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: sla_a4(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0.0000015332445_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: sla_lower(max_lct_id) = [ &
    & 5.325_wp, 0._wp, 0._wp, 7.193295071270618_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0.006295_wp, 0._wp, 0._wp, 0._wp ]
  REAL(wp), PARAMETER :: sla_upper(max_lct_id) = [ &
    & 21.20910493387759_wp, 0._wp, 0._wp, 35.342860388508186_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0.0323334_wp, 0._wp, 0._wp, 0._wp ]

  CHARACTER(len=*), PARAMETER :: modname = 'mo_trait_process'

CONTAINS

  ! ======================================================================================================== !
  !>
  !>#### Compute the annual trait diagnostics of one plant functional type
  !>
  !> We here evaluate the trait equation for each of the three traits that has a fitted
  !> model for the given lctlib id, and limit the result to the fitted bounds. Traits
  !> without a model are returned as zero and have to be identified by the caller
  !> through has_dynamic_vcmax25, has_dynamic_jmax25 and has_dynamic_sla rather than
  !> from the returned value.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ! Renamed from calc_pft13_trait_diagnostics; takes the lctlib id and
  !   evaluates only the traits that have a fitted model for that type.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-09-23
  ELEMENTAL PURE SUBROUTINE calc_trait_diagnostics( &
    & lib_id, mat_prev_year, map_prev_year, max_mmtemp_prev_year, &
    & vcmax25_raw, vcmax25_pred,                                  &
    & jmax25_raw, jmax25_pred,                                    &
    & sla_raw, sla_pred, jmax_vcmax_ratio)

    INTEGER,  INTENT(in)  ::    &
      & lib_id                    !< Id of the tile's primary land cover type in lctlib

    REAL(wp), INTENT(in)  ::    &
      & mat_prev_year,          & !< Previous-year mean air temperature [K]
      & map_prev_year,          & !< Previous-year total precipitation [kg m-2]
      & max_mmtemp_prev_year      !< Previous-year warmest monthly mean temperature [K]

    REAL(wp), INTENT(out) ::    &
      & vcmax25_raw,            & !< Unbounded Vcmax25 diagnostic [umol m-2 s-1]
      & vcmax25_pred,           & !< Bounded Vcmax25 diagnostic [umol m-2 s-1]
      & jmax25_raw,             & !< Unbounded Jmax25 diagnostic [umol m-2 s-1]
      & jmax25_pred,            & !< Bounded Jmax25 diagnostic [umol m-2 s-1]
      & sla_raw,                & !< Unbounded SLA diagnostic [fitted unit]
      & sla_pred,               & !< Bounded SLA diagnostic [fitted unit]
      & jmax_vcmax_ratio          !< Unbounded Jmax25 to Vcmax25 ratio [-]

    ! Local variables
    REAL(wp) ::                 &
      & mat_prev_c,             & !< Previous-year mean air temperature [degC]
      & tmax_prev_c               !< Previous-year warmest monthly mean temperature [degC]
    ! ------------------------------------------------------------------------------------------------------ !

    vcmax25_raw      = 0._wp
    vcmax25_pred     = 0._wp
    jmax25_raw       = 0._wp
    jmax25_pred      = 0._wp
    sla_raw          = 0._wp
    sla_pred         = 0._wp
    jmax_vcmax_ratio = 0._wp

    IF (lib_id < 1 .OR. lib_id > max_lct_id) RETURN

    mat_prev_c  = mat_prev_year - tmelt
    tmax_prev_c = max_mmtemp_prev_year - tmelt

    IF (has_dynamic_vcmax25(lib_id)) THEN
      vcmax25_raw = vcmax25_a0(lib_id)                        &
        & + vcmax25_a1(lib_id) * mat_prev_c                   &
        & + vcmax25_a2(lib_id) * map_prev_year                &
        & + vcmax25_a3(lib_id) * tmax_prev_c                  &
        & + vcmax25_a4(lib_id) * map_prev_year * tmax_prev_c
      vcmax25_pred = MIN(vcmax25_upper(lib_id), MAX(vcmax25_lower(lib_id), vcmax25_raw))
    END IF

    IF (has_dynamic_jmax25(lib_id)) THEN
      jmax25_raw = jmax25_a0(lib_id)                          &
        & + jmax25_a1(lib_id) * mat_prev_c                    &
        & + jmax25_a2(lib_id) * map_prev_year                 &
        & + jmax25_a3(lib_id) * tmax_prev_c                   &
        & + jmax25_a4(lib_id) * map_prev_year * tmax_prev_c
      jmax25_pred = MIN(jmax25_upper(lib_id), MAX(jmax25_lower(lib_id), jmax25_raw))
    END IF

    IF (has_dynamic_sla(lib_id)) THEN
      sla_raw = sla_a0(lib_id)                                &
        & + sla_a1(lib_id) * mat_prev_c                       &
        & + sla_a2(lib_id) * map_prev_year                    &
        & + sla_a3(lib_id) * tmax_prev_c                      &
        & + sla_a4(lib_id) * map_prev_year * tmax_prev_c
      sla_pred = MIN(sla_upper(lib_id), MAX(sla_lower(lib_id), sla_raw))
    END IF

    ! The ratio is a diagnostic and is formed from the unbounded values, as in the
    ! previous PFT13 implementation. It is defined only where both traits have a model.
    IF (has_dynamic_vcmax25(lib_id) .AND. has_dynamic_jmax25(lib_id)) THEN
      IF (ABS(vcmax25_raw) > EPSILON(1._wp)) THEN
        jmax_vcmax_ratio = jmax25_raw / vcmax25_raw
      ELSE
        jmax_vcmax_ratio = 0._wp
      END IF
    END IF

  END SUBROUTINE calc_trait_diagnostics

#endif
END MODULE mo_trait_process
