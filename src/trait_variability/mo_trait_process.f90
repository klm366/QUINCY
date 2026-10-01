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
!> observations of the PFT, and Jmax25 is further kept inside the 95 % prediction
!> interval of the Jmax25-Vcmax25 relationship of the PFT (Verheijen et al. 2013).
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
  !   sla      id  1: Table 3 SLA PFT01, 84 community means, R2adj 0.158
  !   sla      id  2: Table 3 SLA PFT02, 49 community means, R2adj 0.261
  !   sla      id  3: Table 3 SLA PFT03, 306 community means, R2adj 0.129
  !   sla      id  4: Table 3 SLA PFT04, 298 community means, R2adj 0.385
  !   sla      id  9: Table 3 SLA PFT05, 5 community means, R2adj 0.922
  !   sla      id 10: Table 3 SLA PFT06, 222 community means, R2adj 0.159
  !   sla      id 11: Table 3 SLA PFT07, 617 community means, R2adj 0.051
  !   sla      id 12: Table 3 SLA PFT08, 168 community means, R2adj 0.281
  !   sla      id 21: validated PFT13 equation (Zhao et al. 2025 data), x 1000 to m2 kgC-1; refit before E4-E6
  !   vcmax25  id  1: Table 3 Vcmax25 TROP_TREE, 43 community means, R2adj 0.722
  !   vcmax25  id  2: Table 3 Vcmax25 TROP_TREE, 43 community means, R2adj 0.722
  !   vcmax25  id  3: Table 3 Vcmax25 PFT03, 67 community means, R2adj 0.148
  !   vcmax25  id  4: Table 3 Vcmax25 PFT04, 40 community means, R2adj 0.103
  !   vcmax25  id  9: Table 3 Vcmax25 SHRUB, 9 community means, R2adj 0.453
  !   vcmax25  id 10: Table 3 Vcmax25 SHRUB, 9 community means, R2adj 0.453
  !   vcmax25  id 11: Table 3 Vcmax25 PFT07, 16 community means, R2adj 0.660
  !   jmax25   id  1: Table 3 Jmax25 TROP_TREE, 25 community means, R2adj 0.793
  !   jmax25   id  2: Table 3 Jmax25 TROP_TREE, 25 community means, R2adj 0.793
  !   jmax25   id  3: Table 3 Jmax25 PFT03, 30 community means, R2adj 0.181
  !   jmax25   id  4: Table 3 Jmax25 PFT04, 22 community means, R2adj 0.334
  !   jmax25   id  9: Table 3 Jmax25 SHRUB, 6 community means, R2adj 0.987
  !   jmax25   id 10: Table 3 Jmax25 SHRUB, 6 community means, R2adj 0.987
  !   jmax25   id 11: Table 3 Jmax25 PFT07, 12 community means, R2adj 0.607

  INTEGER, PARAMETER :: n_pred = 7   !< MAT, MAP, Tmin, Tmax, Reh, SoilMoist, NSWR

  !> specificLeafArea_C [m2 molC-1] = SLA [m2 kgC-1] * 1e-3 * molar_mass_C [g molC-1]
  REAL(wp), PARAMETER :: trait_sla_to_specificleafarea_c(max_lct_id) = &
    & MERGE(1.e-3_wp * molar_mass_C, 0._wp, [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &
    & .FALSE., .FALSE., .TRUE. , .FALSE., .FALSE., .FALSE. ] )

  !---------------------------------------- vcmax25 [umol m-2 s-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_vcmax25(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .FALSE., &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: vcmax25_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -253.0663040852799_wp, 15.00394165720491_wp, 0._wp, -7.556128984938753_wp, &
    & 0._wp, 0._wp, 0._wp, 0.5145896398886425_wp, &   ! id 1
    & -253.0663040852799_wp, 15.00394165720491_wp, 0._wp, -7.556128984938753_wp, &
    & 0._wp, 0._wp, 0._wp, 0.5145896398886425_wp, &   ! id 2
    & 117.9024145761717_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, -122.5822804387149_wp, -0.1945306920803238_wp, &   ! id 3
    & -48.65089210040083_wp, 0._wp, 0._wp, 0._wp, &
    & 4.633025926592828_wp, 0._wp, 397.7900858942706_wp, 0._wp, &   ! id 4
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 5
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 6
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 7
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 8
    & -171.6651414517909_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 3.154184615531433_wp, 0._wp, 0._wp, &   ! id 9
    & -171.6651414517909_wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 3.154184615531433_wp, 0._wp, 0._wp, &   ! id 10
    & -223.0456797436773_wp, 0._wp, -0.1015575924786065_wp, 0._wp, &
    & 0._wp, 3.325336493258569_wp, 0._wp, 1.068648468466633_wp, &   ! id 11
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 12
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
    & 0, 0, 0, 0, 0, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: vcmax25_int_j(max_lct_id) = [ &
    & 0, 0, 0, 6, 0, 0, &   ! 1-6
    & 0, 0, 0, 0, 0, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  REAL(wp), PARAMETER :: vcmax25_int_c(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, &   ! 1-3
    & -18.29635245698318_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 0._wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: vcmax25_lower(max_lct_id) = [ &
    & 11.527_wp, 17.81604085742264_wp, 13.9_wp, &   ! 1-3
    & 14.34127924501023_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 33.4456130959505_wp, &   ! 7-9
    & 9.84634988002226_wp, 25.80500000000001_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: vcmax25_upper(max_lct_id) = [ &
    & 100.1885653526678_wp, 136.6696587743677_wp, 150._wp, &   ! 1-3
    & 137.6355633903134_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 98.85785366991574_wp, &   ! 7-9
    & 162.2853080659865_wp, 159.8355688734351_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24

  !---------------------------------------- jmax25 [umol m-2 s-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_jmax25(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .FALSE., &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: jmax25_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -448.9012318641051_wp, 30.94002184522258_wp, 0._wp, -11.73518063906627_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 1
    & -448.9012318641051_wp, 30.94002184522258_wp, 0._wp, -11.73518063906627_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 2
    & 275.8241080163108_wp, 0._wp, -0.0267167765190061_wp, 3.092970649851722_wp, &
    & -6.206379507256733_wp, 0._wp, 0._wp, 0._wp, &   ! id 3
    & 419.1128890998157_wp, 0._wp, -0.4148109305940413_wp, 0._wp, &
    & 0._wp, -3.705326375937275_wp, 0._wp, 0._wp, &   ! id 4
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 5
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 6
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 7
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 8
    & -1376.685831060944_wp, -1.509603958452587_wp, 0._wp, 0._wp, &
    & 0._wp, 20.31211961559538_wp, 0._wp, 0._wp, &   ! id 9
    & -1376.685831060944_wp, -1.509603958452587_wp, 0._wp, 0._wp, &
    & 0._wp, 20.31211961559538_wp, 0._wp, 0._wp, &   ! id 10
    & 266.6722368246149_wp, 0._wp, -0.1309068619392587_wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 11
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 12
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
    & 0, 0, 0, 0, 0, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: jmax25_int_j(max_lct_id) = [ &
    & 0, 0, 0, 5, 0, 0, &   ! 1-6
    & 0, 0, 0, 0, 0, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 0, 0, 0, 0 ]   ! 19-24
  REAL(wp), PARAMETER :: jmax25_int_c(max_lct_id) = [ &
    & 0._wp, 0._wp, 0._wp, &   ! 1-3
    & 0.0049779151120457_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & 0._wp, 0._wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jmax25_lower(max_lct_id) = [ &
    & 19.97766083435633_wp, 39.05392636413955_wp, 30.178_wp, &   ! 1-3
    & 38.01053582465638_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 55.26207246250373_wp, &   ! 7-9
    & 29.8617672157933_wp, 87.60923430593337_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jmax25_upper(max_lct_id) = [ &
    & 157.7534218077332_wp, 203.9390792626152_wp, 288.8185598211251_wp, &   ! 1-3
    & 234.5847499999999_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 139.0621272011965_wp, &   ! 7-9
    & 363.8726624883944_wp, 405.0749528077252_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24

  !---------------------------------------- sla [m2 kgC-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_sla(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .TRUE. , .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: sla_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -264.6117368399383_wp, 0._wp, 0._wp, 10.33250222803353_wp, &
    & 0._wp, 4.180136475906664_wp, 0._wp, 0._wp, &   ! id 1
    & -365.3639897185811_wp, 0._wp, 0._wp, 13.82997170773674_wp, &
    & 0._wp, 5.793254405902545_wp, 0._wp, 0._wp, &   ! id 2
    & 13.56126062360548_wp, 0._wp, 0._wp, 2.030778774025066_wp, &
    & 0.4057461812019727_wp, 0._wp, 0._wp, 0._wp, &   ! id 3
    & 34.21442679658357_wp, 0._wp, 0.0125931316422182_wp, 0.577795280554252_wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 4
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
    & 28.02674937543264_wp, 0._wp, 0._wp, 1.655310094473117_wp, &
    & 0._wp, 0._wp, 0._wp, 0.0992435244330258_wp, &   ! id 10
    & -6.276717779767665_wp, 0._wp, 0._wp, 0._wp, &
    & 2.434364223211217_wp, 0._wp, 0._wp, 0.5388930344167117_wp, &   ! id 11
    & 27.15646333733383_wp, 2.041478408581396_wp, 0._wp, 0._wp, &
    & -2.502033923273533_wp, 0.6263801514712164_wp, 0._wp, 0._wp, &   ! id 12
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
    & 60.960327_wp, 0._wp, -0.045172134_wp, 0._wp, &
    & -1.5723334_wp, 0._wp, 0._wp, 0._wp, &   ! id 21
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 22
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp, &   ! id 23
    & 0._wp, 0._wp, 0._wp, 0._wp, &
    & 0._wp, 0._wp, 0._wp, 0._wp  &   ! id 24
    & ], [n_pred + 1, max_lct_id] )
  INTEGER,  PARAMETER :: sla_int_i(max_lct_id) = [ &
    & 3, 3, 3, 2, 0, 0, &   ! 1-6
    & 0, 0, 0, 3, 4, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 2, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: sla_int_j(max_lct_id) = [ &
    & 5, 5, 4, 3, 0, 0, &   ! 1-6
    & 0, 0, 0, 7, 7, 0, &   ! 7-12
    & 0, 0, 0, 0, 0, 0, &   ! 13-18
    & 0, 0, 4, 0, 0, 0 ]   ! 19-24
  REAL(wp), PARAMETER :: sla_int_c(max_lct_id) = [ &
    & -0.1494397974485725_wp, -0.2056242448377997_wp, -0.1009197492173917_wp, &   ! 1-3
    & -0.0019184209736277_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 0._wp, &   ! 7-9
    & -0.0179888020229775_wp, -0.0247122616635386_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0.0015332445_wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: sla_lower(max_lct_id) = [ &
    & 13.20334757916667_wp, 12.50438243829334_wp, 5.389355316558443_wp, &   ! 1-3
    & 20.21256818030039_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 14.77604166666667_wp, &   ! 7-9
    & 18.12314799788327_wp, 14.90500308125_wp, 16.39583333333334_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 6.295_wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: sla_upper(max_lct_id) = [ &
    & 71.28134515611995_wp, 79.15244789062507_wp, 60.95868055555551_wp, &   ! 1-3
    & 108.3470271114079_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 64.69364743392471_wp, &   ! 7-9
    & 91.32900833333336_wp, 121.4789133602352_wp, 88.21763080454095_wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 32.3334_wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24

  !------------------------------ Jmax25 ~ Vcmax25 constraint (95 % prediction interval) ------------------------------!
  LOGICAL,  PARAMETER :: has_jv_constraint(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .FALSE., &   ! 7-12
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(wp), PARAMETER :: jv_intercept(max_lct_id) = [ &
    & 45.25889879043241_wp, 45.25889879043241_wp, 20.65816067404732_wp, &   ! 1-3
    & 32.85235719187544_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 7.623259478946235_wp, &   ! 7-9
    & 7.623259478946235_wp, 56.81274011540388_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jv_slope(max_lct_id) = [ &
    & 0.8646143127609888_wp, 0.8646143127609888_wp, 1.789236573781073_wp, &   ! 1-3
    & 1.269745209871273_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 2.062374327615337_wp, &   ! 7-9
    & 2.062374327615337_wp, 1.855105060605188_wp, 0._wp, &   ! 10-12
    & 0._wp, 0._wp, 0._wp, &   ! 13-15
    & 0._wp, 0._wp, 0._wp, &   ! 16-18
    & 0._wp, 0._wp, 0._wp, &   ! 19-21
    & 0._wp, 0._wp, 0._wp ]   ! 22-24
  REAL(wp), PARAMETER :: jv_halfwidth(max_lct_id) = [ &
    & 65.39367679802368_wp, 65.39367679802368_wp, 85.76928970066835_wp, &   ! 1-3
    & 61.29195258447776_wp, 0._wp, 0._wp, &   ! 4-6
    & 0._wp, 0._wp, 100.1374661337129_wp, &   ! 7-9
    & 100.1374661337129_wp, 85.87038700008657_wp, 0._wp, &   ! 10-12
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
    REAL(wp) :: jv_centre         !< Jmax25 expected from Vcmax25 [umol m-2 s-1]
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
      ! Second constraint of Verheijen et al. (2013): Jmax25 stays inside the 95 % prediction
      ! interval of the Jmax25-Vcmax25 relationship of paired observations of the PFT. The
      ! PFT bounds are applied once more afterwards and take precedence.
      IF (has_dynamic_vcmax25(lib_id) .AND. has_jv_constraint(lib_id)) THEN
        jv_centre   = jv_intercept(lib_id) + jv_slope(lib_id) * vcmax25_pred
        jmax25_pred = MIN(jv_centre + jv_halfwidth(lib_id), MAX(jv_centre - jv_halfwidth(lib_id), jmax25_pred))
        jmax25_pred = MIN(jmax25_upper(lib_id), MAX(jmax25_lower(lib_id), jmax25_pred))
      END IF
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

  !> Evaluate c0 + sum_k c_k x_k + c_int x_i x_j (no interaction when int_i = 0)
  PURE REAL(wp) FUNCTION eval_model(c, int_i, int_j, int_c, x) RESULT(y)
    REAL(wp), INTENT(in) :: c(0:n_pred), int_c, x(n_pred)
    INTEGER,  INTENT(in) :: int_i, int_j
    y = c(0) + SUM(c(1:n_pred) * x(1:n_pred))
    IF (int_i > 0 .AND. int_j > 0) y = y + int_c * x(int_i) * x(int_j)
  END FUNCTION eval_model

#endif
END MODULE mo_trait_process
