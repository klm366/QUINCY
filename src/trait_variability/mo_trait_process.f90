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
!> Rewritten for the trait-climate relationships of all JSBACH plant functional types
!>   (manuscript 1, Table 3). Author: Dr. Sri Nandini-Weiß. Date: 2026-10-01.
!>
!> One equation form is evaluated for every lctlib id that has a fitted model:
!>
!>     trait_raw = c0 + sum_k c_k * x_k + c_int * x_i * x_j
!>
!> with the seven predictors of Verheijen et al. (2013), all of the previous year and
!> provided by CLIMBUF:  x1 MAT [degC], x2 MAP [mm], x3 Tmin = coldest monthly mean
!> temperature [degC], x4 Tmax = warmest monthly mean temperature [degC], x5 Reh = air
!> relative humidity [%], x6 SoilMoist = volumetric soil water 0-1 m [m3 m-3], x7 NSWR =
!> net surface shortwave radiation [W m-2]. At most one interaction (x_i * x_j) per model.
!> The prediction is limited to the 2.5 and 97.5 percentiles of the individual
!> observations of the fitted group (the pooled group for pooled relationships), and the
!> pair (Vcmax25, Jmax25) is kept inside the 95 % confidence interval of the Jmax25-Vcmax25
!> regression of the PFT, moved to its border by the shortest distance (Verheijen et al. 2013).
!> C4 grass: Vcmax25 and the PEPcase CO2 specificity (held in jmax25, mmol m-2 s-1) come from
!> one leaf-nitrogen relationship (Simioni et al. 2004) and need no constraint.
!>
!> A PFT without a fitted model is skipped through the has_dynamic_* masks and keeps
!> its static lctlib parameters. The coefficient tables are generated, never typed:
!> Manuscript1/scripts/trait_variability_12PFT/make_trait_fortran_table_20261001.py
!>
MODULE mo_trait_process
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
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

  !--- BEGIN GENERATED TRAIT TABLE (make_trait_fortran_table_20261001.py) ---!
  ! Generated from Table 3 of manuscript 1 by
  !   Manuscript1/scripts/trait_variability_12PFT/make_trait_fortran_table_20261001.py
  ! Do not edit by hand; the README of that script gives the evidence for every entry.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-10-01
  !
  ! Switched-on models (lctlib id: source):
  !   sla      id  1: Table 3 SLA PFT01, 84 community means, R2adj 0.155
  !   sla      id  2: Table 3 SLA PFT02, 49 community means, R2adj 0.289
  !   sla      id  3: Table 3 SLA PFT03, 306 community means, R2adj 0.092
  !   sla      id  4: Table 3 SLA PFT04, 298 community means, R2adj 0.301
  !   sla      id  9: Table 3 SLA PFT05, 5 community means, R2adj 0.922
  !   sla      id 10: Table 3 SLA PFT06, 222 community means, R2adj 0.157
  !   sla      id 11: Table 3 SLA PFT07, 617 community means, R2adj 0.055
  !   sla      id 12: Table 3 SLA PFT08, 168 community means, R2adj 0.270
  !   vcmax25  id  1: Table 3 Vcmax25 TROP_TREE, 41 community means, R2adj 0.829
  !   vcmax25  id  2: Table 3 Vcmax25 TROP_TREE, 41 community means, R2adj 0.829
  !   vcmax25  id  3: Table 3 Vcmax25 PFT03, 67 community means, R2adj 0.135
  !   vcmax25  id  4: Table 3 Vcmax25 PFT04, 40 community means, R2adj 0.163
  !   vcmax25  id  9: Table 3 Vcmax25 SHRUB, 5 community means, R2adj 0.825, not significant
  !   vcmax25  id 10: Table 3 Vcmax25 SHRUB, 5 community means, R2adj 0.825, not significant
  !   vcmax25  id 11: Table 3 Vcmax25 PFT07, 14 community means, R2adj 0.887
  !   vcmax25  id 12: Table 3 Vcmax25 PFT08, 15 community means, R2adj 0.766
  !   jmax25   id  1: Table 3 Jmax25 TROP_TREE, 25 community means, R2adj 0.858
  !   jmax25   id  2: Table 3 Jmax25 TROP_TREE, 25 community means, R2adj 0.858
  !   jmax25   id  3: Table 3 Jmax25 PFT03, 30 community means, R2adj 0.212
  !   jmax25   id  4: Table 3 Jmax25 PFT04, 22 community means, R2adj 0.391
  !   jmax25   id  9: Table 3 Jmax25 SHRUB, 6 community means, R2adj 0.987
  !   jmax25   id 10: Table 3 Jmax25 SHRUB, 6 community means, R2adj 0.987
  !   jmax25   id 11: Table 3 Jmax25 PFT07, 12 community means, R2adj 0.615
  !   jmax25   id 12: Table 3 PEPcase from leaf N PFT08, 15 community means, R2adj 0.766

  INTEGER, PARAMETER :: n_pred = 7   !< MAT, MAP, Tmin, Tmax, Reh, SoilMoist, NSWR

  !> specificLeafArea_C [m2 molC-1] = SLA [m2 kgC-1] * 1e-3 * molar_mass_C [g molC-1]
  REAL(wp), PARAMETER :: trait_sla_to_specificleafarea_c(max_lct_id) = &
    & MERGE(1.e-3_wp * molar_mass_C, 0._wp, [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ] )

  !---------------------------------------- vcmax25 [umol m-2 s-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_vcmax25(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: vcmax25_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -268.2952351659212_wp, 0._wp, 0._wp, -4.862990762997263_wp, &
    & 12.57237362356458_wp, 0._wp, 0._wp, 0.5116157495673235_wp, &   ! id 1
    & -268.2952351659212_wp, 0._wp, 0._wp, -4.862990762997263_wp, &
    & 12.57237362356458_wp, 0._wp, 0._wp, 0.5116157495673235_wp, &   ! id 2
    & 120.6554017963595_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, -124.86975723384_wp, -0.2044553300440425_wp, &   ! id 3
    & -102.1426150227984_wp, 0._wp, 0._wp, 0._wp, &
    & 7.094408678047357_wp, 0._wp, 547.1730980883472_wp, 0._wp, &   ! id 4
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 5
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 6
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 7
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 8
    & -59.5774183852557_wp, 0._wp, 0._wp, 0._wp, &
    & -3.242551573283955_wp, 2.24038378759219_wp, 0._wp, 0._wp, &   ! id 9
    & -59.5774183852557_wp, 0._wp, 0._wp, 0._wp, &
    & -3.242551573283955_wp, 2.24038378759219_wp, 0._wp, 0._wp, &   ! id 10
    & 2.374977236694917_wp, 0._wp, 0._wp, 7.686198996894808_wp, &
    & 0._wp, 0.9922995842311706_wp, 0._wp, 0._wp, &   ! id 11
    & 175.1408103061371_wp, 0._wp, -0.197104038103008_wp, 0._wp, &
    & 0._wp, -2.076646529285447_wp, 0._wp, 0._wp, &   ! id 12
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 13
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 14
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 15
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 16
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 17
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 18
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 19
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 20
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 21
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 22
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 23
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp  &   ! id 24
    & ], [n_pred + 1, max_lct_id] )
  INTEGER,  PARAMETER :: vcmax25_int_i(max_lct_id) = [ &
    & 0, 0, 0, 4, 0, 0, &   ! 1-6
    & 0, 0, 0, 0, 3, 2, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: vcmax25_int_j(max_lct_id) = [ &
    & 0, 0, 0, 6, 0, 0, &   ! 1-6
    & 0, 0, 0, 0, 5, 5, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  REAL(wp), PARAMETER :: vcmax25_int_c(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, &   ! 1-3
    & -24.74431500787971_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, -0.1050158360569416_wp, 0.0028202557733853_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: vcmax25_lower(max_lct_id) = [ &
    & 11.09456459260558_wp, 11.09456459260558_wp, 13.77350876449804_wp, &   ! 1-3
    & 14.40481642809861_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 10.29968882335603_wp, &   ! 7-9
    & 10.29968882335603_wp, 37.92741993263682_wp, 23.74083909673173_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: vcmax25_upper(max_lct_id) = [ &
    & 133.454178116485_wp, 133.454178116485_wp, 156.4429685544303_wp, &   ! 1-3
    & 136.273_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 167.2304899549145_wp, &   ! 7-9
    & 167.2304899549145_wp, 168.1426608546983_wp, 52.64179268817207_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24

  !---------------------------------------- jmax25 [umol m-2 s-1; C4: PEP mmol m-2 s-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_jmax25(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: jmax25_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -476.6850550569813_wp, 0._wp, 0._wp, 0._wp, &
    & 19.82287274392588_wp, 0._wp, 0._wp, 0._wp, &   ! id 1
    & -476.6850550569813_wp, 0._wp, 0._wp, 0._wp, &
    & 19.82287274392588_wp, 0._wp, 0._wp, 0._wp, &   ! id 2
    & 294.9168012466707_wp, 0._wp, -0.0274658976587964_wp, 3.539684104192398_wp, &
    & -6.829757743673446_wp, 0._wp, 0._wp, 0._wp, &   ! id 3
    & 494.4991500309876_wp, 0._wp, -0.4313493271800985_wp, 0._wp, &
    & 0._wp, -4.704245824306836_wp, 0._wp, 0._wp, &   ! id 4
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 5
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 6
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 7
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 8
    & -1770.989345629837_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 25.61412288588961_wp, 0._wp, 0._wp, &   ! id 9
    & -1770.989345629837_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 25.61412288588961_wp, 0._wp, 0._wp, &   ! id 10
    & 273.7641064035606_wp, 0._wp, -0.1310088294983611_wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 11
    & 2522.134421282032_wp, 0._wp, -2.987729131787312_wp, 0._wp, &
    & 0._wp, -31.47808330912353_wp, 0._wp, 0._wp, &   ! id 12
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 13
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 14
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 15
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 16
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 17
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 18
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 19
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 20
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 21
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 22
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 23
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp  &   ! id 24
    & ], [n_pred + 1, max_lct_id] )
  INTEGER,  PARAMETER :: jmax25_int_i(max_lct_id) = [ &
    & 0, 0, 0, 2, 0, 0, &   ! 1-6
    & 0, 0, 0, 0, 0, 2, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: jmax25_int_j(max_lct_id) = [ &
    & 0, 0, 0, 5, 0, 0, &   ! 1-6
    & 0, 0, 0, 0, 0, 5, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  REAL(wp), PARAMETER :: jmax25_int_c(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, &   ! 1-3
    & 0.0051982327221508_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 0._wp, 0.0427498107818121_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jmax25_lower(max_lct_id) = [ &
    & 20.73588319665885_wp, 20.73588319665885_wp, 33.579190081284_wp, &   ! 1-3
    & 36.53395095735517_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 29.8617672157933_wp, &   ! 7-9
    & 29.8617672157933_wp, 83.5069053092781_wp, 227.1935926334411_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jmax25_upper(max_lct_id) = [ &
    & 204.3555049219492_wp, 204.3555049219492_wp, 295.8259711717423_wp, &   ! 1-3
    & 236.4235308120971_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 371.7935210624482_wp, &   ! 7-9
    & 371.7935210624482_wp, 417.891511317337_wp, 665.2780774193552_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24

  !---------------------------------------- sla [m2 kgC-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_sla(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: sla_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -279.7208907751893_wp, 0._wp, 0._wp, 11.00277273719332_wp, &
    & 0._wp, 4.424054756857863_wp, 0._wp, 0._wp, &   ! id 1
    & -406.0998584188323_wp, 0._wp, 0._wp, 15.45502369445305_wp, &
    & 0._wp, 6.420833363075849_wp, 0._wp, 0._wp, &   ! id 2
    & 75.19276203563669_wp, 0._wp, 0._wp, 0._wp, &
    & -2.772030730255589_wp, 0._wp, -158.1843977319195_wp, 0._wp, &   ! id 3
    & 73.16215807976064_wp, 0._wp, 0.009399337683634001_wp, -1.187593454477916_wp, &
    & 0._wp, -0.3691255590465796_wp, -36.32915146457854_wp, 0._wp, &   ! id 4
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 5
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 6
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 7
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 8
    & -98.97344850178872_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 461.4794287684795_wp, 0._wp, &   ! id 9
    & 27.48125129593543_wp, 0._wp, 0._wp, 1.604995601772025_wp, &
    & 0._wp, 0._wp, 0._wp, 0.1025912983876993_wp, &   ! id 10
    & -10.56585456353008_wp, 0._wp, 0._wp, 0._wp, &
    & 2.91497453424556_wp, 0._wp, 0._wp, 0.5604027026654583_wp, &   ! id 11
    & 30.05443835739599_wp, 2.009185168544826_wp, 0._wp, 0._wp, &
    & -2.534103368832049_wp, 0.5997004777083121_wp, 0._wp, 0._wp, &   ! id 12
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 13
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 14
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 15
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 16
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 17
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 18
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 19
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 20
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 21
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 22
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 23
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp  &   ! id 24
    & ], [n_pred + 1, max_lct_id] )
  INTEGER,  PARAMETER :: sla_int_i(max_lct_id) = [ &
    & 3, 3, 4, 0, 0, 0, &   ! 1-6
    & 0, 0, 0, 3, 4, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: sla_int_j(max_lct_id) = [ &
    & 5, 5, 6, 0, 0, 0, &   ! 1-6
    & 0, 0, 0, 7, 7, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  REAL(wp), PARAMETER :: sla_int_c(max_lct_id) = [ &
    & -0.160340199470438_wp, -0.231027021415858_wp, 8.127477467234918_wp, &   ! 1-3
    & 0._wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & -0.0176790628763099_wp, -0.0278107710366178_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: sla_lower(max_lct_id) = [ &
    & 12.76212526736464_wp, 12.50232484868917_wp, 5.731088334816514_wp, &   ! 1-3
    & 18.59976356252941_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 14.77604166666667_wp, &   ! 7-9
    & 18.45375_wp, 14.90500308125_wp, 16.06348834269771_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: sla_upper(max_lct_id) = [ &
    & 67.93377784952843_wp, 71.85600272791854_wp, 63.79884431401091_wp, &   ! 1-3
    & 106.4120408257519_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 64.69364743392471_wp, &   ! 7-9
    & 91.47745416666666_wp, 120.4988486111614_wp, 82.72740632793503_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24

  !---------------- Jmax25 = a + b Vcmax25, 95 % confidence interval: t s sqrt(1/W + (V - Vw)^2 / Sxx) ----------------!
  LOGICAL,  PARAMETER :: has_jv_constraint(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .FALSE., .FALSE., .TRUE. , .FALSE., &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: jv_intercept(max_lct_id) = [ &
    & 17.87304671918965_wp, 18.70992415341316_wp, 32.48282779576203_wp, &   ! 1-3
    & 38.91858321372095_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, -76.55931157954257_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jv_slope(max_lct_id) = [ &
    & 1.125757390041126_wp, 1.234139038529975_wp, 1.69288851212829_wp, &   ! 1-3
    & 1.252978849096706_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 3.014455281628512_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jv_t_s(max_lct_id) = [ &
    & 30.8224813587535_wp, 59.43005447764674_wp, 117.6510884281517_wp, &   ! 1-3
    & 132.6648417933358_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 208.0170310190944_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jv_inv_w(max_lct_id) = [ &
    & 0.0212123907555919_wp, 0.02939835701329015_wp, 0.01216904666114307_wp, &   ! 1-3
    & 0.01072158827303688_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 0.0216791477999037_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jv_vmean(max_lct_id) = [ &
    & 52.30823593890563_wp, 81.6552952753792_wp, 53.32215087809054_wp, &   ! 1-3
    & 63.73299165052104_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 93.86806651222535_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jv_inv_sxx(max_lct_id) = [ &
    & 2.937787992122675e-05_wp, 3.735029833044135e-05_wp, 1.966795432646996e-05_wp, &   ! 1-3
    & 3.316494896695333e-05_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 5.347846451836339e-05_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  !--- END GENERATED TRAIT TABLE ---!

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
  ! Seven predictors, one optional interaction, per-PFT bounds and the Jmax25-Vcmax25
  !   constraint of Verheijen et al. (2013).
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-10-01
  ELEMENTAL PURE SUBROUTINE calc_trait_diagnostics( &
    & lib_id, mat_prev_year, map_prev_year, min_mmtemp_prev_year, max_mmtemp_prev_year, &
    & reh_prev_year, soil_moist_prev_year, nswr_prev_year,       &
    & vcmax25_raw, vcmax25_pred,                                  &
    & jmax25_raw, jmax25_pred,                                    &
    & sla_raw, sla_pred, jmax_vcmax_ratio)

    INTEGER,  INTENT(in)  ::    &
      & lib_id                    !< Id of the tile's primary land cover type in lctlib

    REAL(wp), INTENT(in)  ::    &
      & mat_prev_year,          & !< Previous-year mean air temperature [K]
      & map_prev_year,          & !< Previous-year total precipitation [kg m-2] = [mm]
      & min_mmtemp_prev_year,   & !< Previous-year coldest monthly mean temperature [K]
      & max_mmtemp_prev_year,   & !< Previous-year warmest monthly mean temperature [K]
      & reh_prev_year,          & !< Previous-year mean air relative humidity [%]
      & soil_moist_prev_year,   & !< Previous-year mean volumetric soil water 0-1 m [m3 m-3]
      & nswr_prev_year            !< Previous-year mean net surface shortwave radiation [W m-2]

    REAL(wp), INTENT(out) ::    &
      & vcmax25_raw,            & !< Unbounded Vcmax25 diagnostic [umol m-2 s-1]
      & vcmax25_pred,           & !< Bounded Vcmax25 diagnostic [umol m-2 s-1]
      & jmax25_raw,             & !< Unbounded Jmax25 diagnostic [umol m-2 s-1]
      & jmax25_pred,            & !< Bounded and Vcmax25-constrained Jmax25 diagnostic [umol m-2 s-1]
      & sla_raw,                & !< Unbounded SLA diagnostic [m2 kg-1 C]
      & sla_pred,               & !< Bounded SLA diagnostic [m2 kg-1 C]
      & jmax_vcmax_ratio          !< Bounded Jmax25 to Vcmax25 ratio [-]

    REAL(wp) :: x(n_pred)         !< Predictors in the units of the fitted equations
    ! ------------------------------------------------------------------------------------------------------ !

    vcmax25_raw      = 0._wp
    vcmax25_pred     = 0._wp
    jmax25_raw       = 0._wp
    jmax25_pred      = 0._wp
    sla_raw          = 0._wp
    sla_pred         = 0._wp
    jmax_vcmax_ratio = 0._wp

    IF (lib_id < 1 .OR. lib_id > max_lct_id) RETURN
    ! The added predictors exist only after one complete year (or a restart written by this
    ! version); relative humidity is positive whenever they exist. Until then all traits
    ! stay 0 and the consuming processes keep the static lctlib values.
    IF (reh_prev_year <= 0._wp) RETURN

    x(1) = mat_prev_year - tmelt
    x(2) = map_prev_year
    x(3) = min_mmtemp_prev_year - tmelt
    x(4) = max_mmtemp_prev_year - tmelt
    x(5) = reh_prev_year
    x(6) = soil_moist_prev_year
    x(7) = nswr_prev_year

    IF (has_dynamic_vcmax25(lib_id)) THEN
      vcmax25_raw  = eval_model(vcmax25_coef(:,lib_id), vcmax25_int_i(lib_id), vcmax25_int_j(lib_id), &
        &                       vcmax25_int_c(lib_id), x)
      vcmax25_pred = MIN(vcmax25_upper(lib_id), MAX(vcmax25_lower(lib_id), vcmax25_raw))
    END IF

    IF (has_dynamic_jmax25(lib_id)) THEN
      jmax25_raw  = eval_model(jmax25_coef(:,lib_id), jmax25_int_i(lib_id), jmax25_int_j(lib_id), &
        &                      jmax25_int_c(lib_id), x)
      jmax25_pred = MIN(jmax25_upper(lib_id), MAX(jmax25_lower(lib_id), jmax25_raw))
    END IF

    ! Second constraint of Verheijen et al. (2013, Sect. 2.3): the pair (Vcmax25, Jmax25) stays
    ! inside the 95 % confidence interval of the regression Jmax25 = a + b Vcmax25 of the PFT;
    ! a pair outside is moved to the nearest point of the interval border ("adjusted to the
    ! confidence interval border, based on the shortest distance needed to reach the border").
    ! Both traits move. The bounds are applied once more afterwards and take precedence.
    IF (has_dynamic_vcmax25(lib_id) .AND. has_dynamic_jmax25(lib_id) .AND. has_jv_constraint(lib_id)) THEN
      CALL project_on_jv_band(lib_id, vcmax25_pred, jmax25_pred)
      vcmax25_pred = MIN(vcmax25_upper(lib_id), MAX(vcmax25_lower(lib_id), vcmax25_pred))
      jmax25_pred  = MIN(jmax25_upper(lib_id),  MAX(jmax25_lower(lib_id),  jmax25_pred))
    END IF

    IF (has_dynamic_sla(lib_id)) THEN
      sla_raw  = eval_model(sla_coef(:,lib_id), sla_int_i(lib_id), sla_int_j(lib_id), &
        &                   sla_int_c(lib_id), x)
      sla_pred = MIN(sla_upper(lib_id), MAX(sla_lower(lib_id), sla_raw))
    END IF

    IF (has_dynamic_vcmax25(lib_id) .AND. has_dynamic_jmax25(lib_id) .AND. vcmax25_pred > 0._wp) THEN
      jmax_vcmax_ratio = jmax25_pred / vcmax25_pred
    END IF

  END SUBROUTINE calc_trait_diagnostics

  !> Half-width of the 95 % confidence interval of the fitted line Jmax25 = a + b Vcmax25 at v:
  !>   t(0.975, n-2) s sqrt(1/W + (v - Vw)^2 / Sxx)   (weighted least squares on community means)
  PURE REAL(wp) FUNCTION jv_halfwidth(lib_id, v) RESULT(h)
    INTEGER,  INTENT(in) :: lib_id
    REAL(wp), INTENT(in) :: v
    h = jv_t_s(lib_id) * SQRT(jv_inv_w(lib_id) + (v - jv_vmean(lib_id))**2 * jv_inv_sxx(lib_id))
  END FUNCTION jv_halfwidth

  !> Move (v, j) to the nearest point of the border of the confidence band when it lies outside.
  !> The nearest border point lies within the vertical distance d0 of the point, so the search
  !> covers v - d0 ... v + d0: 200 evenly spaced trial points, then 60 golden-section steps
  !> around the best one. The squared distance is (v' - v)^2 + (a + b v' +- h(v') - j)^2.
  PURE SUBROUTINE project_on_jv_band(lib_id, v, j)
    INTEGER,  INTENT(in)    :: lib_id
    REAL(wp), INTENT(inout) :: v, j
    INTEGER,  PARAMETER :: n_trial = 200, n_golden = 60
    REAL(wp), PARAMETER :: gr = 0.6180339887498949_wp
    REAL(wp) :: dev, sgn, d0, step, vt, dbest, vbest, lo, hi, x1, x2, f1, f2
    INTEGER  :: k
    dev = j - (jv_intercept(lib_id) + jv_slope(lib_id) * v)
    IF (ABS(dev) <= jv_halfwidth(lib_id, v)) RETURN
    sgn = SIGN(1._wp, dev)
    d0  = ABS(dev) - jv_halfwidth(lib_id, v)
    step = 2._wp * d0 / REAL(n_trial, wp)
    vbest = v
    dbest = dist2(v)
    DO k = 0, n_trial
      vt = v - d0 + REAL(k, wp) * step
      IF (dist2(vt) < dbest) THEN
        dbest = dist2(vt)
        vbest = vt
      END IF
    END DO
    lo = vbest - step
    hi = vbest + step
    DO k = 1, n_golden
      x1 = hi - gr * (hi - lo)
      x2 = lo + gr * (hi - lo)
      f1 = dist2(x1)
      f2 = dist2(x2)
      IF (f1 < f2) THEN
        hi = x2
      ELSE
        lo = x1
      END IF
    END DO
    v = 0.5_wp * (lo + hi)
    j = jv_intercept(lib_id) + jv_slope(lib_id) * v + sgn * jv_halfwidth(lib_id, v)
  CONTAINS
    PURE REAL(wp) FUNCTION dist2(vv)
      REAL(wp), INTENT(in) :: vv
      dist2 = (vv - v)**2 + (jv_intercept(lib_id) + jv_slope(lib_id) * vv + sgn * jv_halfwidth(lib_id, vv) - j)**2
    END FUNCTION dist2
  END SUBROUTINE project_on_jv_band

  !> Evaluate c0 + sum_k c_k x_k + c_int x_i x_j (no interaction when int_i = 0)
  PURE REAL(wp) FUNCTION eval_model(c, int_i, int_j, int_c, x) RESULT(y)
    REAL(wp), INTENT(in) :: c(0:n_pred), int_c, x(n_pred)
    INTEGER,  INTENT(in) :: int_i, int_j
    y = c(0) + SUM(c(1:n_pred) * x(1:n_pred))
    IF (int_i > 0 .AND. int_j > 0) y = y + int_c * x(int_i) * x(int_j)
  END FUNCTION eval_model

#endif
END MODULE mo_trait_process
