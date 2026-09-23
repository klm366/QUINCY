! Author: Dr. Sri Nandini-Weiß
! Date: 2026-07-15
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
!>#### Calculates diagnostic PFT13 traits from completed previous-year climate
!>
MODULE mo_trait_process
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_physical_constants, ONLY: tmelt, molar_mass_C

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_pft13_trait_diagnostics
  PUBLIC :: pft13_vcmax25_to_carboxrate
  PUBLIC :: pft13_jmax25_to_etransport
  PUBLIC :: pft13_sla_to_specificleafarea_c

  ! Step 6 unit conversions for active JSBACH use.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-07-27
  ! Vcmax25 and Jmax25 diagnostics are stored in umol m-2 s-1, whereas the
  ! assimilation routines use mol m-2 s-1. The SLA regression is stored in
  ! m2 gC-1, whereas specificLeafArea_C is m2 molC-1.
  ! Convert bounded annual trait predictions to the units required by
  ! the existing JSBACH assimilation and SLA process interfaces.
  ! Convert bounded annual PFT13 predictions to the units required by
  ! the existing JSBACH assimilation and SLA process interfaces.
  REAL(wp), PARAMETER :: pft13_vcmax25_to_carboxrate      = 1.e-6_wp
  REAL(wp), PARAMETER :: pft13_jmax25_to_etransport      = 1.e-6_wp
  REAL(wp), PARAMETER :: pft13_sla_to_specificleafarea_c = molar_mass_C

  REAL(wp), PARAMETER :: pft13_vcmax25_intercept = 49.5690779208497_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_mat_coef  = -0.391997574087549_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_map_coef  = 0.001594603433354_wp
  REAL(wp), PARAMETER :: pft13_jmax25_intercept  = 180.411473978400_wp
  REAL(wp), PARAMETER :: pft13_jmax25_mat_coef   = -3.09078405961564_wp
  REAL(wp), PARAMETER :: pft13_jmax25_map_coef   = -0.00485175158785102_wp
  REAL(wp), PARAMETER :: pft13_sla_intercept     = 0.060960327_wp
  REAL(wp), PARAMETER :: pft13_sla_map_coef      = -0.000045172134_wp
  REAL(wp), PARAMETER :: pft13_sla_tmax_coef     = -0.0015723334_wp
  REAL(wp), PARAMETER :: pft13_sla_map_tmax_coef = 0.0000015332445_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_lower     = 16.5700_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_upper     = 66.2500_wp
  REAL(wp), PARAMETER :: pft13_jmax25_lower      = 46.3775_wp
  REAL(wp), PARAMETER :: pft13_jmax25_upper      = 131.4725_wp
  REAL(wp), PARAMETER :: pft13_sla_lower         = 0.006295_wp
  REAL(wp), PARAMETER :: pft13_sla_upper         = 0.0323334_wp

CONTAINS

  ELEMENTAL PURE SUBROUTINE calc_pft13_trait_diagnostics( &
    & mat_prev_year, map_prev_year, max_mmtemp_prev_year, &
    & pft13_vcmax25_raw, pft13_vcmax25_pred,             &
    & pft13_jmax25_raw, pft13_jmax25_pred,               &
    & pft13_sla_raw, pft13_sla_pred, pft13_jmax_vcmax_ratio)

    REAL(wp), INTENT(in) :: &
      & mat_prev_year,          & !< Previous-year mean air temperature [K]
      & map_prev_year,          & !< Previous-year total precipitation [kg m-2]
      & max_mmtemp_prev_year      !< Previous-year warmest monthly mean temperature [K]

    REAL(wp), INTENT(out) :: &
      & pft13_vcmax25_raw,      &
      & pft13_vcmax25_pred,     &
      & pft13_jmax25_raw,       &
      & pft13_jmax25_pred,      &
      & pft13_sla_raw,          &
      & pft13_sla_pred,         &
      & pft13_jmax_vcmax_ratio

    REAL(wp) :: mat_prev_c, tmax_prev_c

    mat_prev_c  = mat_prev_year - tmelt
    tmax_prev_c = max_mmtemp_prev_year - tmelt

    pft13_vcmax25_raw = pft13_vcmax25_intercept &
      & + pft13_vcmax25_mat_coef * mat_prev_c    &
      & + pft13_vcmax25_map_coef * map_prev_year

    pft13_jmax25_raw = pft13_jmax25_intercept &
      & + pft13_jmax25_mat_coef * mat_prev_c  &
      & + pft13_jmax25_map_coef * map_prev_year

    pft13_sla_raw = pft13_sla_intercept                         &
      & + pft13_sla_map_coef * map_prev_year                    &
      & + pft13_sla_tmax_coef * tmax_prev_c                     &
      & + pft13_sla_map_tmax_coef * map_prev_year * tmax_prev_c

    pft13_vcmax25_pred = MIN(pft13_vcmax25_upper, &
      & MAX(pft13_vcmax25_lower, pft13_vcmax25_raw))
    pft13_jmax25_pred = MIN(pft13_jmax25_upper, &
      & MAX(pft13_jmax25_lower, pft13_jmax25_raw))
    pft13_sla_pred = MIN(pft13_sla_upper, &
      & MAX(pft13_sla_lower, pft13_sla_raw))

    IF (ABS(pft13_vcmax25_raw) > EPSILON(1._wp)) THEN
      pft13_jmax_vcmax_ratio = pft13_jmax25_raw / pft13_vcmax25_raw
    ELSE
      pft13_jmax_vcmax_ratio = 0._wp
    END IF

  END SUBROUTINE calc_pft13_trait_diagnostics

#endif
END MODULE mo_trait_process
