! Author: Dr. Sri Nandini-Weiß
! Date: 2026-07-15
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

  USE mo_kind,             ONLY: wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,    ONLY: VEG_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_trait_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 7

  TYPE, EXTENDS(t_jsb_memory) :: t_trait_memory
    TYPE(t_jsb_var_real2d) :: &
      & pft13_vcmax25_raw,        & !< Unbounded PFT13 Vcmax25 diagnostic [umol m-2 s-1]
      & pft13_vcmax25_pred,       & !< Bounded PFT13 Vcmax25 diagnostic [umol m-2 s-1]
      & pft13_jmax25_raw,         & !< Unbounded PFT13 Jmax25 diagnostic [umol m-2 s-1]
      & pft13_jmax25_pred,        & !< Bounded PFT13 Jmax25 diagnostic [umol m-2 s-1]
      & pft13_sla_raw,            & !< Unbounded PFT13 SLA diagnostic [m2 gC-1]
      & pft13_sla_pred,           & !< Bounded PFT13 SLA diagnostic [m2 gC-1]
      & pft13_jmax_vcmax_ratio      !< Unbounded PFT13 Jmax25 to Vcmax25 ratio [-]
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
    CHARACTER(len=*),      INTENT(in)            :: prefix
    CHARACTER(len=*),      INTENT(in)            :: suffix
    INTEGER,               INTENT(in)            :: lct_ids(:)
    INTEGER,               INTENT(in)            :: lib_id
    INTEGER,               INTENT(in)            :: model_id

    TYPE(t_jsb_grid),  POINTER :: hgrid
    TYPE(t_jsb_vgrid), POINTER :: surface
    INTEGER :: table

    table = tables(1)
    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    IF (ANY(lct_ids == VEG_TYPE)) THEN

      CALL mem%Add_var('pft13_vcmax25_raw', mem%pft13_vcmax25_raw,                         &
        & hgrid, surface,                                                                  &
        & t_cf('pft13_vcmax25_raw', 'umol m-2 s-1', 'Unbounded PFT13 Vcmax25 diagnostic'),&
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix, output_level=MEDIUM,                                             &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('pft13_vcmax25_pred', mem%pft13_vcmax25_pred,                       &
        & hgrid, surface,                                                                  &
        & t_cf('pft13_vcmax25_pred', 'umol m-2 s-1', 'Bounded PFT13 Vcmax25 diagnostic'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix, output_level=MEDIUM,                                             &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('pft13_jmax25_raw', mem%pft13_jmax25_raw,                           &
        & hgrid, surface,                                                                  &
        & t_cf('pft13_jmax25_raw', 'umol m-2 s-1', 'Unbounded PFT13 Jmax25 diagnostic'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix, output_level=MEDIUM,                                             &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('pft13_jmax25_pred', mem%pft13_jmax25_pred,                         &
        & hgrid, surface,                                                                  &
        & t_cf('pft13_jmax25_pred', 'umol m-2 s-1', 'Bounded PFT13 Jmax25 diagnostic'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix, output_level=MEDIUM,                                             &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('pft13_sla_raw', mem%pft13_sla_raw,                                 &
        & hgrid, surface,                                                                  &
        & t_cf('pft13_sla_raw', 'm2 gC-1', 'Unbounded PFT13 SLA diagnostic'),              &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix, output_level=MEDIUM,                                             &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('pft13_sla_pred', mem%pft13_sla_pred,                               &
        & hgrid, surface,                                                                  &
        & t_cf('pft13_sla_pred', 'm2 gC-1', 'Bounded PFT13 SLA diagnostic'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix, output_level=MEDIUM,                                             &
        & lrestart=.TRUE., initval_r=0.0_wp)

      CALL mem%Add_var('pft13_jmax_vcmax_ratio', mem%pft13_jmax_vcmax_ratio,               &
        & hgrid, surface,                                                                  &
        & t_cf('pft13_jmax_vcmax_ratio', '1', 'Unbounded PFT13 Jmax25 to Vcmax25 ratio'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix, output_level=MEDIUM,                                             &
        & lrestart=.TRUE., initval_r=0.0_wp)

    END IF

  END SUBROUTINE Init_trait_memory

#endif
END MODULE mo_trait_memory_class
