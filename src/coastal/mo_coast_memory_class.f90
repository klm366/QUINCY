! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Memory class for the coastal-constraint process
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
!>#### Contains memory definitions for the coastal control factors
!>
!> The boundary conditions are read from bc_land_coastal.nc by mo_coast_init: the mask of
!> the cells that carry mangrove, the mangrove-area-weighted long-term mean sea surface
!> salinity from Bio-ORACLE v2.0, and the complex harmonic
!> constants of eight tidal constituents from the TPXO9-atlas solution. The remaining
!> variables are recomputed every time step. All variables are allocated on vegetation
!> tiles only.
!>
MODULE mo_coast_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,             ONLY: wp
  USE mo_jsb_memory_class, ONLY: t_jsb_memory
  USE mo_jsb_var_class,    ONLY: t_jsb_var_real2d
  USE mo_jsb_lct_class,    ONLY: VEG_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_coast_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 25

  TYPE, EXTENDS(t_jsb_memory) :: t_coast_memory
    TYPE(t_jsb_var_real2d) :: &
      & coast_mask,           & !< Cells on which the coastal process is evaluated (1) or not (0) [-]
      & salinity,             & !< Long-term mean sea surface salinity [PSS]
      & tide_hre_m2,          & !< Real part of the complex M2 elevation amplitude [m]
      & tide_him_m2,          & !< Imaginary part of the complex M2 elevation amplitude [m]
      & tide_hre_s2,          & !< Real part of the complex S2 elevation amplitude [m]
      & tide_him_s2,          & !< Imaginary part of the complex S2 elevation amplitude [m]
      & tide_hre_n2,          & !< Real part of the complex N2 elevation amplitude [m]
      & tide_him_n2,          & !< Imaginary part of the complex N2 elevation amplitude [m]
      & tide_hre_k2,          & !< Real part of the complex K2 elevation amplitude [m]
      & tide_him_k2,          & !< Imaginary part of the complex K2 elevation amplitude [m]
      & tide_hre_k1,          & !< Real part of the complex K1 elevation amplitude [m]
      & tide_him_k1,          & !< Imaginary part of the complex K1 elevation amplitude [m]
      & tide_hre_o1,          & !< Real part of the complex O1 elevation amplitude [m]
      & tide_him_o1,          & !< Imaginary part of the complex O1 elevation amplitude [m]
      & tide_hre_p1,          & !< Real part of the complex P1 elevation amplitude [m]
      & tide_him_p1,          & !< Imaginary part of the complex P1 elevation amplitude [m]
      & tide_hre_q1,          & !< Real part of the complex Q1 elevation amplitude [m]
      & tide_him_q1,          & !< Imaginary part of the complex Q1 elevation amplitude [m]
      & tide_height,          & !< Tide height above mean sea level [m]
      & bgr_fract_inundated,  & !< Inundated fraction of belowground root [-]
      & bgr_fract_ventilated, & !< Ventilated fraction of belowground root [-]
      & bgr_fract_anoxic,     & !< Anoxic fraction of belowground root [-]
      & control_salinity,     & !< Control factor of salinity on assimilation [-]
      & control_inundation,   & !< Control factor of tidal inundation on assimilation [-]
      & control_coastal         !< Combined coastal factor on assimilation [-]
  CONTAINS
    PROCEDURE :: Init => Init_coast_memory
  END TYPE t_coast_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_coast_memory_class'

CONTAINS

  ! ======================================================================================================== !
  !>
  !>#### Initialisation of the coastal-constraint process memory
  !>
  !> The subroutine is a type-bound procedure of [[t_coast_memory]] and is called during
  !> initialisation for all tiles on which the coastal process runs.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE Init_coast_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,    ONLY: MEDIUM
    USE mo_jsb_io,         ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class, ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,       ONLY: Get_grid, Get_vgrid

    CLASS(t_coast_memory), INTENT(inout), TARGET :: mem

    CHARACTER(len=*), INTENT(in) :: &
      & prefix,                     & !< Prefix for variable names, i.e. process name
      & suffix                        !< Suffix for variable names, i.e. tile name

    INTEGER, INTENT(in) :: &
      & lct_ids(:),        & !< Land cover type ids
      & lib_id,            & !< Id of the tile's primary land cover type
      & model_id             !< Id of the current model instance

    TYPE(t_jsb_grid),  POINTER :: hgrid
    TYPE(t_jsb_vgrid), POINTER :: surface
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_coast_memory'
    ! ------------------------------------------------------------------------------------------------------ !

    table = tables(1)
    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    IF (ANY(lct_ids == VEG_TYPE)) THEN

      ! Boundary conditions. They are read from bc_land_coastal.nc at the start of every
      ! simulation and are therefore not written to the restart file.
      CALL mem%Add_var('coast_mask', mem%coast_mask,                                        &
        & hgrid, surface,                                                                    &
        & t_cf('coast_mask', '1', 'Cells on which the coastal process is evaluated'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('salinity', mem%salinity,                                            &
        & hgrid, surface,                                                                    &
        & t_cf('salinity', '1e-3', 'Long-term mean sea surface salinity'),                  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_m2', mem%tide_hre_m2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_m2', 'm', 'Real part of the M2 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_m2', mem%tide_him_m2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_m2', 'm', 'Imaginary part of the M2 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_s2', mem%tide_hre_s2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_s2', 'm', 'Real part of the S2 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_s2', mem%tide_him_s2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_s2', 'm', 'Imaginary part of the S2 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_n2', mem%tide_hre_n2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_n2', 'm', 'Real part of the N2 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_n2', mem%tide_him_n2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_n2', 'm', 'Imaginary part of the N2 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_k2', mem%tide_hre_k2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_k2', 'm', 'Real part of the K2 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_k2', mem%tide_him_k2,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_k2', 'm', 'Imaginary part of the K2 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_k1', mem%tide_hre_k1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_k1', 'm', 'Real part of the K1 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_k1', mem%tide_him_k1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_k1', 'm', 'Imaginary part of the K1 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_o1', mem%tide_hre_o1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_o1', 'm', 'Real part of the O1 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_o1', mem%tide_him_o1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_o1', 'm', 'Imaginary part of the O1 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_p1', mem%tide_hre_p1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_p1', 'm', 'Real part of the P1 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_p1', mem%tide_him_p1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_p1', 'm', 'Imaginary part of the P1 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_hre_q1', mem%tide_hre_q1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_hre_q1', 'm', 'Real part of the Q1 elevation amplitude'),                &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('tide_him_q1', mem%tide_him_q1,                                        &
        & hgrid, surface,                                                                     &
        & t_cf('tide_him_q1', 'm', 'Imaginary part of the Q1 elevation amplitude'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                  &
        & prefix, suffix, output_level=MEDIUM,                                                &
        & lrestart=.FALSE., initval_r=0.0_wp)


      ! Quantities recomputed every time step
      CALL mem%Add_var('tide_height', mem%tide_height,                                      &
        & hgrid, surface,                                                                    &
        & t_cf('tide_height', 'm', 'Tide height above mean sea level'),                     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('bgr_fract_inundated', mem%bgr_fract_inundated,                      &
        & hgrid, surface,                                                                    &
        & t_cf('bgr_fract_inundated', '-', 'Inundated fraction of belowground root'),       &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('bgr_fract_ventilated', mem%bgr_fract_ventilated,                    &
        & hgrid, surface,                                                                    &
        & t_cf('bgr_fract_ventilated', '-', 'Ventilated fraction of belowground root'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('bgr_fract_anoxic', mem%bgr_fract_anoxic,                            &
        & hgrid, surface,                                                                    &
        & t_cf('bgr_fract_anoxic', '-', 'Anoxic fraction of belowground root'),             &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=0.0_wp)

      CALL mem%Add_var('control_salinity', mem%control_salinity,                            &
        & hgrid, surface,                                                                    &
        & t_cf('control_salinity', '-', 'Control factor of salinity on assimilation'),      &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=1.0_wp)

      CALL mem%Add_var('control_inundation', mem%control_inundation,                        &
        & hgrid, surface,                                                                    &
        & t_cf('control_inundation', '-', 'Control factor of tidal inundation'),            &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=1.0_wp)

      CALL mem%Add_var('control_coastal', mem%control_coastal,                              &
        & hgrid, surface,                                                                    &
        & t_cf('control_coastal', '-', 'Combined coastal factor on assimilation'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
        & prefix, suffix, output_level=MEDIUM,                                               &
        & lrestart=.FALSE., initval_r=1.0_wp)

    END IF


  END SUBROUTINE Init_coast_memory

#endif
END MODULE mo_coast_memory_class
