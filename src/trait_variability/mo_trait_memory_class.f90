! Author: Dr. Sri Nandini-Weiß
! Date: 2026-07-15 (generic names for all PFTs: 2026-10-01)
!> Memory class for the trait-variability process
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
!>#### Contains memory definitions for climate-dependent trait diagnostics
!>
MODULE mo_trait_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: lp, wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,    ONLY: VEG_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_trait_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 7

  ! Trait diagnostics of the PFT of each vegetation tile. Every PFT tile has its own copy, and
  ! the output name carries the tile name in front of the variable:
  !   trait_<tile>_<variable>, e.g. trait_pft04_sla_pred, trait_pft13_vcmax25_pred, trait_veg_sla_pred
  ! (the box and veg tiles hold the fraction-weighted mean of their child tiles). The short name
  ! (e.g. sla_pred) is the same on every tile, as the aggregation to the parent tile requires
  ! (mo_jsb_tile.f90 finds the child variable by its short name).
  ! A value of 0 means: no climate relationship for this PFT; the consumers then keep lctlib.
  TYPE, EXTENDS(t_jsb_memory) :: t_trait_memory
    TYPE(t_jsb_var_real2d) :: &
      & vcmax25_raw,        & !< Vcmax25 before the bounds [umol m-2 s-1]
      & vcmax25_pred,       & !< Vcmax25 after the bounds and the Jmax25-Vcmax25 constraint [umol m-2 s-1]
      & jmax25_raw,         & !< Jmax25 before the bounds [umol m-2 s-1]; C4: PEPcase CO2 specificity [mmol m-2 s-1]
      & jmax25_pred,        & !< Jmax25 after the bounds and the constraint [umol m-2 s-1]; C4: PEP [mmol m-2 s-1]
      & sla_raw,            & !< SLA before the bounds [m2 kg-1 leaf carbon]
      & sla_pred,           & !< SLA after the bounds [m2 kg-1 leaf carbon]
      & jmax_vcmax_ratio      !< jmax25_raw / vcmax25_raw [-]
  CONTAINS
    PROCEDURE :: Init => Init_trait_memory
  END TYPE t_trait_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_trait_memory_class'

CONTAINS

  SUBROUTINE Init_trait_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,    ONLY: MEDIUM
    USE mo_jsb_io,         ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,       ONLY: Get_grid, Get_vgrid

    CLASS(t_trait_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),      INTENT(in)            :: prefix     !< process name, 'trait'
    CHARACTER(len=*),      INTENT(in)            :: suffix     !< tile name, e.g. 'pft04'
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: lib_id
    INTEGER,               INTENT(in)            :: model_id

    TYPE(t_jsb_grid),  POINTER :: hgrid
    TYPE(t_jsb_vgrid), POINTER :: surface
    INTEGER :: table
    CHARACTER(len=LEN(prefix)+LEN(suffix)+1) :: pre   !< 'trait_pft04': the tile name goes before the variable

    table = tables(1)
    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    pre = TRIM(prefix)
    IF (suffix /= '') pre = TRIM(prefix)//'_'//TRIM(suffix)

    IF (ANY(lct_ids == VEG_TYPE)) THEN

      CALL mem%Add_var('vcmax25_raw', mem%vcmax25_raw,                                       &
        & hgrid, surface,                                                                  &
        & t_cf('vcmax25_raw', 'umol m-2 s-1', 'Vcmax25 of the PFT before the bounds'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & TRIM(pre), '', output_level=MEDIUM,                                              &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('vcmax25_pred', mem%vcmax25_pred,                                     &
        & hgrid, surface,                                                                  &
        & t_cf('vcmax25_pred', 'umol m-2 s-1', 'Vcmax25 of the PFT after bounds and Jmax25-Vcmax25 constraint'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & TRIM(pre), '', output_level=MEDIUM,                                              &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('jmax25_raw', mem%jmax25_raw,                                         &
        & hgrid, surface,                                                                  &
        & t_cf('jmax25_raw', 'umol m-2 s-1 (C4: mmol m-2 s-1)',                            &
        &      'Jmax25 of the PFT before the bounds (C4: PEPcase CO2 specificity)'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & TRIM(pre), '', output_level=MEDIUM,                                              &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('jmax25_pred', mem%jmax25_pred,                                       &
        & hgrid, surface,                                                                  &
        & t_cf('jmax25_pred', 'umol m-2 s-1 (C4: mmol m-2 s-1)',                           &
        &      'Jmax25 of the PFT after bounds and constraint (C4: PEPcase CO2 specificity)'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & TRIM(pre), '', output_level=MEDIUM,                                              &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('sla_raw', mem%sla_raw,                                               &
        & hgrid, surface,                                                                  &
        & t_cf('sla_raw', 'm2 kg-1', 'SLA of the PFT before the bounds, per kg leaf carbon'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & TRIM(pre), '', output_level=MEDIUM,                                              &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('sla_pred', mem%sla_pred,                                             &
        & hgrid, surface,                                                                  &
        & t_cf('sla_pred', 'm2 kg-1', 'SLA of the PFT after the bounds, per kg leaf carbon'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & TRIM(pre), '', output_level=MEDIUM,                                              &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('jmax_vcmax_ratio', mem%jmax_vcmax_ratio,                             &
        & hgrid, surface,                                                                  &
        & t_cf('jmax_vcmax_ratio', '1', 'Ratio jmax25_raw / vcmax25_raw of the PFT'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & TRIM(pre), '', output_level=MEDIUM,                                              &
        & lrestart=.TRUE., initval_r=0.0_wp)

    END IF

  END SUBROUTINE Init_trait_memory

#endif
END MODULE mo_trait_memory_class
