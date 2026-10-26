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

  USE mo_kind,                   ONLY: lp, wp
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
  REAL(lp), PARAMETER :: trait_vcmax25_to_carboxrate = 1.e-6_lp
  REAL(lp), PARAMETER :: trait_jmax25_to_etransport  = 1.e-6_lp

  !--- BEGIN GENERATED TRAIT TABLE (make_trait_fortran_table_20261001.py) ---!
  ! Generated from Table 3 of manuscript 1 by
  !   Manuscript1/scripts/trait_variability_12PFT/make_trait_fortran_table_20261001.py
  ! Do not edit by hand; the README of that script gives the evidence for every entry.
  ! Author: Dr. Sri Nandini-Weiß
  ! Date: 2026-10-01; id 21 (PFT13) added 2026-10-26 by add_pft13_to_trait_table_20261026.py (Zhao regressions of July 2026)
  !
  ! Switched-on models (lctlib id: source):
  !   sla      id  1: Table 3 SLA PFT01, 96 community means, R2adj 0.267
  !   sla      id  2: Table 3 SLA PFT02, 54 community means, R2adj 0.256
  !   sla      id  3: Table 3 SLA PFT03, 295 community means, R2adj 0.105
  !   sla      id  4: Table 3 SLA PFT04, 214 community means, R2adj 0.139
  !   sla      id  9: Table 3 SLA PFT05, 100 community means, R2adj 0.383
  !   sla      id 10: Table 3 SLA PFT06, 192 community means, R2adj 0.150
  !   sla      id 11: Table 3 SLA PFT07, 410 community means, R2adj 0.043
  !   sla      id 12: Table 3 SLA PFT08, 144 community means, R2adj 0.252
  !   sla      id 14: PFT09 = Table 3 SLA PFT07, 410 community means, R2adj 0.043
  !   sla      id 15: PFT10 = Table 3 SLA PFT08, 144 community means, R2adj 0.252
  !   sla      id 19: PFT11 = Table 3 SLA PFT07, 410 community means, R2adj 0.043
  !   sla      id 20: PFT12 = Table 3 SLA PFT08, 144 community means, R2adj 0.252
  !   sla      id 21: PFT13 = Zhao 2025 SLA weighted MAP + Tmax + MAP x Tmax (July 2026 equation), 7 site means, bounds 2.5-97.5 %
  !   vcmax25  id  1: Table 3 Vcmax25 TROP_TREE, 48 community means, R2adj 0.336
  !   vcmax25  id  2: Table 3 Vcmax25 TROP_TREE, 48 community means, R2adj 0.336
  !   vcmax25  id  3: Table 3 Vcmax25 PFT03, 92 community means, R2adj 0.142
  !   vcmax25  id  4: Table 3 Vcmax25 PFT04, 52 community means, R2adj 0.067
  !   vcmax25  id  9: Table 3 Vcmax25 SHRUB, 24 community means, R2adj 0.586
  !   vcmax25  id 10: Table 3 Vcmax25 SHRUB, 24 community means, R2adj 0.586
  !   vcmax25  id 11: Table 3 Vcmax25 PFT07, 10 community means, R2adj 0.525
  !   vcmax25  id 12: Table 3 Vcmax25 PFT08, 13 community means, R2adj 0.917
  !   vcmax25  id 14: PFT09 = Table 3 Vcmax25 PFT07, 10 community means, R2adj 0.525
  !   vcmax25  id 15: PFT10 = Table 3 Vcmax25 PFT08, 13 community means, R2adj 0.917
  !   vcmax25  id 19: PFT11 = Table 3 Vcmax25 PFT07, 10 community means, R2adj 0.525
  !   vcmax25  id 20: PFT12 = Table 3 Vcmax25 PFT08, 13 community means, R2adj 0.917
  !   vcmax25  id 21: PFT13 = Zhao 2025 Vcmax25 weighted ridge MAT + MAP (alpha 1.58489), 4 site means, bounds 2.5-97.5 %
  !   jmax25   id  1: Table 3 Jmax25 TROP_TREE, 30 community means, R2adj 0.335
  !   jmax25   id  2: Table 3 Jmax25 TROP_TREE, 30 community means, R2adj 0.335
  !   jmax25   id  3: Table 3 Jmax25 PFT03, 53 community means, R2adj 0.178
  !   jmax25   id  4: Table 3 Jmax25 PFT04, 35 community means, R2adj 0.044, not significant
  !   jmax25   id  9: Table 3 Jmax25 SHRUB, 11 community means, R2adj 0.802
  !   jmax25   id 10: Table 3 Jmax25 SHRUB, 11 community means, R2adj 0.802
  !   jmax25   id 11: Table 3 Jmax25 PFT07, 8 community means, R2adj 0.972
  !   jmax25   id 12: Table 3 PEPcase from leaf N PFT08, 13 community means, R2adj 0.917
  !   jmax25   id 14: PFT09 = Table 3 Jmax25 PFT07, 8 community means, R2adj 0.972
  !   jmax25   id 15: PFT10 = Table 3 PEPcase from leaf N PFT08, 13 community means, R2adj 0.917
  !   jmax25   id 19: PFT11 = Table 3 Jmax25 PFT07, 8 community means, R2adj 0.972
  !   jmax25   id 20: PFT12 = Table 3 PEPcase from leaf N PFT08, 13 community means, R2adj 0.917
  !   jmax25   id 21: PFT13 = Zhao 2025 Jmax25 weighted ridge MAT + MAP (alpha 0.001), 4 site means, bounds 2.5-97.5 %

  INTEGER, PARAMETER :: n_pred = 7   !< MAT, MAP, Tmin, Tmax, Reh, SoilMoist, NSWR

  !> specificLeafArea_C [m2 molC-1] = SLA [m2 kgC-1] * 1e-3 * molar_mass_C [g molC-1]
  REAL(lp), PARAMETER :: trait_sla_to_specificleafarea_c(max_lct_id) = &
    & MERGE(1.e-3_lp * molar_mass_C, 0._lp, [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &
    & .FALSE., .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE., &
    & .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE. ] )

  !---------------------------------------- vcmax25 [umol m-2 s-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_vcmax25(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &   ! 7-12
    & .FALSE., .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(lp), PARAMETER :: vcmax25_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -94.47924558931044_lp, 0._lp, 0._lp, 0._lp, &
    & 2.313088537543604_lp, 0._lp, 0._lp, 0.4037436605526389_lp, &   ! id 1
    & -94.47924558931044_lp, 0._lp, 0._lp, 0._lp, &
    & 2.313088537543604_lp, 0._lp, 0._lp, 0.4037436605526389_lp, &   ! id 2
    & -162.9367826194435_lp, 0._lp, 0._lp, 0._lp, &
    & 15.74362074986928_lp, 3.31801949627197_lp, 0._lp, 0._lp, &   ! id 3
    & 90.27981423009771_lp, 0._lp, 0._lp, 0._lp, &
    & -1.757819593693203_lp, 0._lp, 0._lp, 0._lp, &   ! id 4
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 5
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 6
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 7
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 8
    & -472.6884567938937_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 7.270123890870744_lp, 0._lp, 3.249830928976943_lp, &   ! id 9
    & -472.6884567938937_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 7.270123890870744_lp, 0._lp, 3.249830928976943_lp, &   ! id 10
    & -153.9286292158123_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 2.961578531630056_lp, 0._lp, 0._lp, &   ! id 11
    & 174.624404565725_lp, 0._lp, -0.1905342753341689_lp, 0._lp, &
    & 0._lp, -2.112830766726948_lp, 0._lp, 0._lp, &   ! id 12
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 13
    & -153.9286292158123_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 2.961578531630056_lp, 0._lp, 0._lp, &   ! id 14
    & 174.624404565725_lp, 0._lp, -0.1905342753341689_lp, 0._lp, &
    & 0._lp, -2.112830766726948_lp, 0._lp, 0._lp, &   ! id 15
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 16
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 17
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 18
    & -153.9286292158123_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 2.961578531630056_lp, 0._lp, 0._lp, &   ! id 19
    & 174.624404565725_lp, 0._lp, -0.1905342753341689_lp, 0._lp, &
    & 0._lp, -2.112830766726948_lp, 0._lp, 0._lp, &   ! id 20
    & 49.5690779208497_lp, -0.3919975740875487_lp, 0.0015946034333538_lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 21
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 22
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 23
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp  &   ! id 24
    & ], [n_pred + 1, max_lct_id] )
  INTEGER,  PARAMETER :: vcmax25_int_i(max_lct_id) = [ &
    & 0, 0, 4, 0, 0, 0, &   ! 1-6
    & 0, 0, 5, 5, 0, 2, &   ! 7-12
    & 0, 0, 2, 0, 0, 0, &   ! 13-18
    & 0, 2, 0, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: vcmax25_int_j(max_lct_id) = [ &
    & 0, 0, 5, 0, 0, 0, &   ! 1-6
    & 0, 0, 7, 7, 0, 5, &   ! 7-12
    & 0, 0, 5, 0, 0, 0, &   ! 13-18
    & 0, 5, 0, 0, 0, 0 ]   ! 19-24
  REAL(lp), PARAMETER :: vcmax25_int_c(max_lct_id) = [ &
    & 0._lp, 0._lp, -0.2345030224960338_lp, &   ! 1-3
    & 0._lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, -0.0465182543301034_lp, &   ! 7-9
    & -0.0465182543301034_lp, 0._lp, 0.002752778875191_lp, &   ! 10-12
    & 0._lp, 0._lp, 0.002752778875191_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 0._lp, 0.002752778875191_lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: vcmax25_lower(max_lct_id) = [ &
    & 11.06368295050782_lp, 11.06368295050782_lp, 13.28885986192926_lp, &   ! 1-3
    & 14.22689957189845_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 10.6689198554594_lp, &   ! 7-9
    & 10.6689198554594_lp, 30.45208988313193_lp, 23.70323369489726_lp, &   ! 10-12
    & 0._lp, 30.45208988313193_lp, 23.70323369489726_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 30.45208988313193_lp, 23.70323369489726_lp, 16.57_lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: vcmax25_upper(max_lct_id) = [ &
    & 89.00283004858638_lp, 89.00283004858638_lp, 153.6415546244458_lp, &   ! 1-3
    & 123.035_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 139.4915661265199_lp, &   ! 7-9
    & 139.4915661265199_lp, 153.7795019722433_lp, 47.73483935427322_lp, &   ! 10-12
    & 0._lp, 153.7795019722433_lp, 47.73483935427322_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 153.7795019722433_lp, 47.73483935427322_lp, 66.25_lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24

  !---------------------------------------- jmax25 [umol m-2 s-1; C4: PEP mmol m-2 s-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_jmax25(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &   ! 7-12
    & .FALSE., .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(lp), PARAMETER :: jmax25_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & 139.0292395613747_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, -181.5146046148709_lp, 0._lp, &   ! id 1
    & 139.0292395613747_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, -181.5146046148709_lp, 0._lp, &   ! id 2
    & 145.9372851383584_lp, 0._lp, -0.0215855093358501_lp, 2.129086224364318_lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 3
    & 312.8161562535175_lp, 0._lp, -0.2748774426378292_lp, 0._lp, &
    & 0._lp, -2.396253107759837_lp, 0._lp, 0._lp, &   ! id 4
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 5
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 6
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 7
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 8
    & -2022.434634262861_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 30.01856408213668_lp, 0._lp, 11.77794491682234_lp, &   ! id 9
    & -2022.434634262861_lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 30.01856408213668_lp, 0._lp, 11.77794491682234_lp, &   ! id 10
    & 338.3193349218297_lp, 0._lp, -0.179806056159694_lp, 0._lp, &
    & -16.5039898049527_lp, 0._lp, 0._lp, 0._lp, &   ! id 11
    & 2514.30667462955_lp, 0._lp, -2.888143797045799_lp, 0._lp, &
    & 0._lp, -32.02656877576316_lp, 0._lp, 0._lp, &   ! id 12
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 13
    & 338.3193349218297_lp, 0._lp, -0.179806056159694_lp, 0._lp, &
    & -16.5039898049527_lp, 0._lp, 0._lp, 0._lp, &   ! id 14
    & 2514.30667462955_lp, 0._lp, -2.888143797045799_lp, 0._lp, &
    & 0._lp, -32.02656877576316_lp, 0._lp, 0._lp, &   ! id 15
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 16
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 17
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 18
    & 338.3193349218297_lp, 0._lp, -0.179806056159694_lp, 0._lp, &
    & -16.5039898049527_lp, 0._lp, 0._lp, 0._lp, &   ! id 19
    & 2514.30667462955_lp, 0._lp, -2.888143797045799_lp, 0._lp, &
    & 0._lp, -32.02656877576316_lp, 0._lp, 0._lp, &   ! id 20
    & 180.4114739783996_lp, -3.090784059615645_lp, -0.004851751587851_lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 21
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 22
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 23
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp  &   ! id 24
    & ], [n_pred + 1, max_lct_id] )
  INTEGER,  PARAMETER :: jmax25_int_i(max_lct_id) = [ &
    & 0, 0, 0, 2, 0, 0, &   ! 1-6
    & 0, 0, 5, 5, 2, 2, &   ! 7-12
    & 0, 2, 2, 0, 0, 0, &   ! 13-18
    & 2, 2, 0, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: jmax25_int_j(max_lct_id) = [ &
    & 0, 0, 0, 5, 0, 0, &   ! 1-6
    & 0, 0, 7, 7, 4, 5, &   ! 7-12
    & 0, 4, 5, 0, 0, 0, &   ! 13-18
    & 4, 5, 0, 0, 0, 0 ]   ! 19-24
  REAL(lp), PARAMETER :: jmax25_int_c(max_lct_id) = [ &
    & 0._lp, 0._lp, 0._lp, &   ! 1-3
    & 0.0033312767514771_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, -0.1699359375444796_lp, &   ! 7-9
    & -0.1699359375444796_lp, 0.0159034147012709_lp, 0.0417269870162614_lp, &   ! 10-12
    & 0._lp, 0.0159034147012709_lp, 0.0417269870162614_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 0.0159034147012709_lp, 0.0417269870162614_lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: jmax25_lower(max_lct_id) = [ &
    & 24.16951975828193_lp, 24.16951975828193_lp, 28.06337962877644_lp, &   ! 1-3
    & 30.19890222510368_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 23.0215_lp, &   ! 7-9
    & 23.0215_lp, 79.39470606026514_lp, 226.6235649685856_lp, &   ! 10-12
    & 0._lp, 79.39470606026514_lp, 226.6235649685856_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 79.39470606026514_lp, 226.6235649685856_lp, 46.3775_lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: jmax25_upper(max_lct_id) = [ &
    & 149.0475028332668_lp, 149.0475028332668_lp, 278.6637523163492_lp, &   ! 1-3
    & 231.8420771897988_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 336.0337052686671_lp, &   ! 7-9
    & 336.0337052686671_lp, 411.5204761081155_lp, 590.897828465_lp, &   ! 10-12
    & 0._lp, 411.5204761081155_lp, 590.897828465_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 411.5204761081155_lp, 590.897828465_lp, 131.4725_lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24

  !---------------------------------------- sla [m2 kgC-1] ----------------------------------------!
  LOGICAL,  PARAMETER :: has_dynamic_sla(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .TRUE. , &   ! 7-12
    & .FALSE., .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(lp), PARAMETER :: sla_coef(0:n_pred, max_lct_id) = RESHAPE( [ &
    & -105.7553888743224_lp, 3.931508715408319_lp, 0._lp, 0._lp, &
    & 0._lp, 2.124152174332068_lp, 0._lp, 0._lp, &   ! id 1
    & -236.6227638032601_lp, 0._lp, 0._lp, 13.11899982630946_lp, &
    & 9.973174987822761_lp, 0._lp, 0._lp, 0._lp, &   ! id 2
    & 9.162280110357637_lp, 0._lp, 0._lp, 0.324948068332669_lp, &
    & -0.3107325110523031_lp, 0.2050825773850658_lp, 0._lp, 0._lp, &   ! id 3
    & -109.4308979608486_lp, 0._lp, 0._lp, 0._lp, &
    & 8.452387298901678_lp, 0._lp, 0._lp, 1.177019798629968_lp, &   ! id 4
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 5
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 6
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 7
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 8
    & -323.6808203085866_lp, 0._lp, 0._lp, 0._lp, &
    & 11.31704730874332_lp, 4.65938485811432_lp, 0._lp, 0._lp, &   ! id 9
    & 24.10008598074331_lp, 0._lp, 0._lp, 1.064444536088275_lp, &
    & 0._lp, 0._lp, 0._lp, 0.1296978504878302_lp, &   ! id 10
    & 169.4538427659104_lp, 0._lp, 0._lp, 0._lp, &
    & -7.505459505117702_lp, -1.808116256037515_lp, 0._lp, 0._lp, &   ! id 11
    & 46.69062199091803_lp, 1.914970476755935_lp, 0._lp, 0._lp, &
    & -3.041226452355346_lp, 0.538940696629775_lp, 0._lp, 0._lp, &   ! id 12
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 13
    & 169.4538427659104_lp, 0._lp, 0._lp, 0._lp, &
    & -7.505459505117702_lp, -1.808116256037515_lp, 0._lp, 0._lp, &   ! id 14
    & 46.69062199091803_lp, 1.914970476755935_lp, 0._lp, 0._lp, &
    & -3.041226452355346_lp, 0.538940696629775_lp, 0._lp, 0._lp, &   ! id 15
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 16
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 17
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 18
    & 169.4538427659104_lp, 0._lp, 0._lp, 0._lp, &
    & -7.505459505117702_lp, -1.808116256037515_lp, 0._lp, 0._lp, &   ! id 19
    & 46.69062199091803_lp, 1.914970476755935_lp, 0._lp, 0._lp, &
    & -3.041226452355346_lp, 0.538940696629775_lp, 0._lp, 0._lp, &   ! id 20
    & 60.960327_lp, 0._lp, -0.045172134_lp, 0._lp, &
    & -1.5723334_lp, 0._lp, 0._lp, 0._lp, &   ! id 21
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 22
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp, &   ! id 23
    & 0._lp, 0._lp, 0._lp, 0._lp, &
    & 0._lp, 0._lp, 0._lp, 0._lp  &   ! id 24
    & ], [n_pred + 1, max_lct_id] )
  INTEGER,  PARAMETER :: sla_int_i(max_lct_id) = [ &
    & 1, 3, 0, 4, 0, 0, &   ! 1-6
    & 0, 0, 4, 3, 4, 0, &   ! 7-12
    & 0, 4, 0, 0, 0, 0, &   ! 13-18
    & 4, 0, 2, 0, 0, 0 ]   ! 19-24
  INTEGER,  PARAMETER :: sla_int_j(max_lct_id) = [ &
    & 5, 4, 0, 7, 0, 0, &   ! 1-6
    & 0, 0, 5, 7, 5, 0, &   ! 7-12
    & 0, 5, 0, 0, 0, 0, &   ! 13-18
    & 5, 0, 4, 0, 0, 0 ]   ! 19-24
  REAL(lp), PARAMETER :: sla_int_c(max_lct_id) = [ &
    & -0.0667063210536905_lp, -0.4859260510270168_lp, 0._lp, &   ! 1-3
    & -0.0631846740477614_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, -0.1447157009411945_lp, &   ! 7-9
    & -0.0135616523030853_lp, 0.107263811987881_lp, 0._lp, &   ! 10-12
    & 0._lp, 0.107263811987881_lp, 0._lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 0.107263811987881_lp, 0._lp, 0.0015332445_lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: sla_lower(max_lct_id) = [ &
    & 12.15038759689921_lp, 11.40891704138565_lp, 5.87167866757069_lp, &   ! 1-3
    & 18.45562500000001_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 5.757479417163754_lp, &   ! 7-9
    & 19.51249253817709_lp, 9.532108886355051_lp, 10.4242947048611_lp, &   ! 10-12
    & 0._lp, 9.532108886355051_lp, 10.4242947048611_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 9.532108886355051_lp, 10.4242947048611_lp, 6.2949999999999_lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: sla_upper(max_lct_id) = [ &
    & 67.40245284375003_lp, 79.05765961041668_lp, 53.34975961538458_lp, &   ! 1-3
    & 108.3196476874077_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 78.14227609930269_lp, &   ! 7-9
    & 92.7491145833333_lp, 102.3525588768114_lp, 82.83968941998566_lp, &   ! 10-12
    & 0._lp, 102.3525588768114_lp, 82.83968941998566_lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 102.3525588768114_lp, 82.83968941998566_lp, 32.3334_lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24

  !---------------- Jmax25 = a + b Vcmax25, 95 % confidence interval: t s sqrt(1/W + (V - Vw)^2 / Sxx) ----------------!
  LOGICAL,  PARAMETER :: has_jv_constraint(max_lct_id) = [ &
    & .TRUE. , .TRUE. , .TRUE. , .TRUE. , .FALSE., .FALSE., &   ! 1-6
    & .FALSE., .FALSE., .TRUE. , .TRUE. , .TRUE. , .FALSE., &   ! 7-12
    & .FALSE., .TRUE. , .FALSE., .FALSE., .FALSE., .FALSE., &   ! 13-18
    & .TRUE. , .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. ]   ! 19-24
  REAL(lp), PARAMETER :: jv_intercept(max_lct_id) = [ &
    & 17.19790667215345_lp, 1.643979536640135_lp, 16.78972173533669_lp, &   ! 1-3
    & 0.5811002567258718_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 4.009492492908969_lp, &   ! 7-9
    & -8.649271910209256_lp, 16.43213782181125_lp, 0._lp, &   ! 10-12
    & 0._lp, 16.43213782181125_lp, 0._lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 16.43213782181125_lp, 0._lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: jv_slope(max_lct_id) = [ &
    & 1.282613394719342_lp, 1.84312985359049_lp, 1.599309541848464_lp, &   ! 1-3
    & 1.866319668908151_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 1.659710002039026_lp, &   ! 7-9
    & 2.106980986827142_lp, 2.165989675174372_lp, 0._lp, &   ! 10-12
    & 0._lp, 2.165989675174372_lp, 0._lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 2.165989675174372_lp, 0._lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: jv_t_s(max_lct_id) = [ &
    & 49.85370269695549_lp, 143.620899488636_lp, 94.44820931660082_lp, &   ! 1-3
    & 71.32467963167075_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 65.97576784040486_lp, &   ! 7-9
    & 72.63719730660719_lp, 136.4903588893588_lp, 0._lp, &   ! 10-12
    & 0._lp, 136.4903588893588_lp, 0._lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 136.4903588893588_lp, 0._lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: jv_inv_w(max_lct_id) = [ &
    & 0.01050648906193656_lp, 0.02384677584107125_lp, 0.004841312604420015_lp, &   ! 1-3
    & 0.006124271873458078_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 0.06657893853842818_lp, &   ! 7-9
    & 0.04665489445441422_lp, 0.03798376558125448_lp, 0._lp, &   ! 10-12
    & 0._lp, 0.03798376558125448_lp, 0._lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 0.03798376558125448_lp, 0._lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: jv_vmean(max_lct_id) = [ &
    & 46.38518665881921_lp, 49.67744428313673_lp, 67.70495908616182_lp, &   ! 1-3
    & 59.63214065278283_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 46.85444939518062_lp, &   ! 7-9
    & 55.70827223120193_lp, 77.71715018090839_lp, 0._lp, &   ! 10-12
    & 0._lp, 77.71715018090839_lp, 0._lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 77.71715018090839_lp, 0._lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
  REAL(lp), PARAMETER :: jv_inv_sxx(max_lct_id) = [ &
    & 4.098147591517516e-05_lp, 0.0001213869902037006_lp, 5.961480209478039e-06_lp, &   ! 1-3
    & 1.861174616234492e-05_lp, 0._lp, 0._lp, &   ! 4-6
    & 0._lp, 0._lp, 0.0001074514206772897_lp, &   ! 7-9
    & 3.23173236756354e-05_lp, 5.524934557270716e-05_lp, 0._lp, &   ! 10-12
    & 0._lp, 5.524934557270716e-05_lp, 0._lp, &   ! 13-15
    & 0._lp, 0._lp, 0._lp, &   ! 16-18
    & 5.524934557270716e-05_lp, 0._lp, 0._lp, &   ! 19-21
    & 0._lp, 0._lp, 0._lp ]   ! 22-24
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

    REAL(lp), INTENT(in)  ::    &
      & mat_prev_year,          & !< Previous-year mean air temperature [K]
      & map_prev_year,          & !< Previous-year total precipitation [kg m-2] = [mm]
      & min_mmtemp_prev_year,   & !< Previous-year coldest monthly mean temperature [K]
      & max_mmtemp_prev_year,   & !< Previous-year warmest monthly mean temperature [K]
      & reh_prev_year,          & !< Previous-year mean air relative humidity [%]
      & soil_moist_prev_year,   & !< Previous-year mean volumetric soil water 0-1 m [m3 m-3]
      & nswr_prev_year            !< Previous-year mean net surface shortwave radiation [W m-2]

    REAL(lp), INTENT(out) ::    &
      & vcmax25_raw,            & !< Unbounded Vcmax25 diagnostic [umol m-2 s-1]
      & vcmax25_pred,           & !< Bounded Vcmax25 diagnostic [umol m-2 s-1]
      & jmax25_raw,             & !< Unbounded Jmax25 diagnostic [umol m-2 s-1]
      & jmax25_pred,            & !< Bounded and Vcmax25-constrained Jmax25 diagnostic [umol m-2 s-1]
      & sla_raw,                & !< Unbounded SLA diagnostic [m2 kg-1 C]
      & sla_pred,               & !< Bounded SLA diagnostic [m2 kg-1 C]
      & jmax_vcmax_ratio          !< Bounded Jmax25 to Vcmax25 ratio [-]

    REAL(lp) :: x(n_pred)         !< Predictors in the units of the fitted equations
    ! ------------------------------------------------------------------------------------------------------ !

    vcmax25_raw      = 0._lp
    vcmax25_pred     = 0._lp
    jmax25_raw       = 0._lp
    jmax25_pred      = 0._lp
    sla_raw          = 0._lp
    sla_pred         = 0._lp
    jmax_vcmax_ratio = 0._lp

    IF (lib_id < 1 .OR. lib_id > max_lct_id) RETURN
    ! The added predictors exist only after one complete year (or a restart written by this
    ! version); relative humidity is positive whenever they exist. Until then all traits
    ! stay 0 and the consuming processes keep the static lctlib values.
    IF (reh_prev_year <= 0._lp) RETURN

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

    IF (has_dynamic_vcmax25(lib_id) .AND. has_dynamic_jmax25(lib_id) .AND. vcmax25_pred > 0._lp) THEN
      jmax_vcmax_ratio = jmax25_pred / vcmax25_pred
    END IF

  END SUBROUTINE calc_trait_diagnostics

  !> Half-width of the 95 % confidence interval of the fitted line Jmax25 = a + b Vcmax25 at v:
  !>   t(0.975, n-2) s sqrt(1/W + (v - Vw)^2 / Sxx)   (weighted least squares on community means)
  PURE REAL(lp) FUNCTION jv_halfwidth(lib_id, v) RESULT(h)
    INTEGER,  INTENT(in) :: lib_id
    REAL(lp), INTENT(in) :: v
    h = jv_t_s(lib_id) * SQRT(jv_inv_w(lib_id) + (v - jv_vmean(lib_id))**2 * jv_inv_sxx(lib_id))
  END FUNCTION jv_halfwidth

  !> Move (v, j) to the nearest point of the border of the confidence band when it lies outside.
  !> The nearest border point lies within the vertical distance d0 of the point, so the search
  !> covers v - d0 ... v + d0: 200 evenly spaced trial points, then 60 golden-section steps
  !> around the best one. The squared distance is (v' - v)^2 + (a + b v' +- h(v') - j)^2.
  PURE SUBROUTINE project_on_jv_band(lib_id, v, j)
    INTEGER,  INTENT(in)    :: lib_id
    REAL(lp), INTENT(inout) :: v, j
    INTEGER,  PARAMETER :: n_trial = 200, n_golden = 60
    REAL(lp), PARAMETER :: gr = 0.6180339887498949_lp
    REAL(lp) :: dev, sgn, d0, step, vt, dbest, vbest, lo, hi, x1, x2, f1, f2
    INTEGER  :: k
    dev = j - (jv_intercept(lib_id) + jv_slope(lib_id) * v)
    IF (ABS(dev) <= jv_halfwidth(lib_id, v)) RETURN
    sgn = SIGN(1._lp, dev)
    d0  = ABS(dev) - jv_halfwidth(lib_id, v)
    step = 2._lp * d0 / REAL(n_trial, lp)
    vbest = v
    dbest = dist2(v)
    DO k = 0, n_trial
      vt = v - d0 + REAL(k, lp) * step
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
    v = 0.5_lp * (lo + hi)
    j = jv_intercept(lib_id) + jv_slope(lib_id) * v + sgn * jv_halfwidth(lib_id, v)
  CONTAINS
    PURE REAL(lp) FUNCTION dist2(vv)
      REAL(lp), INTENT(in) :: vv
      dist2 = (vv - v)**2 + (jv_intercept(lib_id) + jv_slope(lib_id) * vv + sgn * jv_halfwidth(lib_id, vv) - j)**2
    END FUNCTION dist2
  END SUBROUTINE project_on_jv_band

  !> Evaluate c0 + sum_k c_k x_k + c_int x_i x_j (no interaction when int_i = 0)
  PURE REAL(lp) FUNCTION eval_model(c, int_i, int_j, int_c, x) RESULT(y)
    REAL(lp), INTENT(in) :: c(0:n_pred), int_c, x(n_pred)
    INTEGER,  INTENT(in) :: int_i, int_j
    y = c(0) + SUM(c(1:n_pred) * x(1:n_pred))
    IF (int_i > 0 .AND. int_j > 0) y = y + int_c * x(int_i) * x(int_j)
  END FUNCTION eval_model

#endif
END MODULE mo_trait_process
