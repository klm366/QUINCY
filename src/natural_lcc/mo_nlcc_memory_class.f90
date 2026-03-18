!> Contains the memory class for the nlcc process.
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
!>
MODULE mo_nlcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: !X e.g: LAND_TYPE, LAKE_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: NONE, BASIC, MEDIUM, FULL
  USE mo_jsb_physical_constants, ONLY: !X e.g: tmelt

  ! Use of prcesses in this module
  dsl4jsb_Use_processes NLCC_

  ! Use process configurations
  dsl4jsb_Use_config(NLCC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_nlcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_nlcc_memory

    ! Common variables independent of lct_type
    TYPE(t_jsb_var_real2d) :: &
      & temp_sum_day, &
      & temp_sum_month, &
      & min_mmtemp_of_yr, &
      & max_mmtemp_of_yr, &
      & min_mmtemp20, &
      & max_mmtemp20, &
      & gdd_sum_year, &
      & gdd_prev_year, &
      & max_snow_year, &
      & max_snow_prev_year, &
      & max_thaw_year, &
      & max_thaw_prev_year, &
      & bare_fpc, &
      & sum_green_bio_memory, &
      & desert_fpc
    TYPE(t_jsb_var_real3d) :: &
      & act_fpc, &
      & pot_fpc, &
      & cover_fract_pot, &
      & cover_fract, &
      & npp_pft, &
      & bio_exist
    ! Additional variables for LAND lct_type
!!$  TR  TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & s_star,               & !< Surface dry static energy (s^star, see manual)    [m2 s-2]


    ! Additional variables for LAKE lct_type
!!$ TR    TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & t_wtr,                & !< Lake surface temperature (water)     [K]

  CONTAINS
    PROCEDURE :: Init => Init_nlcc_memory
  END TYPE t_nlcc_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_nlcc_memory_class'

CONTAINS

  SUBROUTINE Init_nlcc_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_nlcc_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: lib_id     !< id of primary lct in lctlib
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, pft_vgrid           ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_nlcc_memory'

    dsl4jsb_Def_config(NLCC_)

    model => Get_model(model_id)

    dsl4jsb_Get_config(NLCC_)

    table = tables(1)

    hgrid     => Get_grid(mem%grid_id)
    surface   => Get_vgrid('surface')
    pft_vgrid => Get_vgrid('pfts')

    CALL mem%Add_var('temp_sum_day', mem%temp_sum_day,                        &
      & hgrid, surface,                                                       &
      & t_cf('temp_sum_day', 'K', 'sum of air temperature in day'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('temp_sum_month', mem%temp_sum_month,                    &
      & hgrid, surface,                                                       &
      & t_cf('temp_sum_month', 'K', 'sum of air temperature in month'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('min_mmtemp_of_yr', mem%min_mmtemp_of_yr,                &
      & hgrid, surface,                                                       &
      & t_cf('min_mmtemp_of_yr', 'K', 'min monthly temperature of year'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=1000.0_wp )

    CALL mem%Add_var('max_mmtemp_of_yr', mem%max_mmtemp_of_yr,                &
      & hgrid, surface,                                                       &
      & t_cf('max_mmtemp_of_yr', 'K', 'max temperature of year'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=-1000.0_wp )

    CALL mem%Add_var('min_mmtemp20', mem%min_mmtemp20,                        &
      & hgrid, surface,                                                       &
      & t_cf('min_mmtemp20', 'K', 'min monthly temperature climatology'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                   &
      & initval_r=1000.0_wp )

    CALL mem%Add_var('max_mmtemp20', mem%max_mmtemp20,                        &
      & hgrid, surface,                                                       &
      & t_cf('max_mmtemp20', 'K', 'max temperature climatology'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                   &
      & initval_r=-1000.0_wp )

    CALL mem%Add_var('gdd_sum_year', mem%gdd_sum_year,                        &
      & hgrid, surface,                                                       &
      & t_cf('gdd_sum_year', 'degC days', 'gdd current year'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('gdd_prev_year', mem%gdd_prev_year,                      &
      & hgrid, surface,                                                       &
      & t_cf('gdd_prev_year', 'degC days', 'gdd previous year'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('max_snow_year', mem%max_snow_year,                      &
      & hgrid, surface,                                                       &
      & t_cf('max_snow_year', 'm', 'max snow depth current year'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('max_snow_prev_year', mem%max_snow_prev_year,            &
      & hgrid, surface,                                                       &
      & t_cf('max_snow_prev_year', 'm', 'max snow depth previous year'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('max_thaw_year', mem%max_thaw_year,                &
      & hgrid, surface,                                                       &
      & t_cf('max_thaw_year', 'm', 'max thaw depth current year'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('max_thaw_prev_year', mem%max_thaw_prev_year,            &
      & hgrid, surface,                                                       &
      & t_cf('max_thaw_prev_year', 'm', 'max thaw depth previous year'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                    &
      & initval_r=99.9_wp )

    IF (dsl4jsb_config(NLCC_)%climbuf_only .OR. &
       .NOT. dsl4jsb_config(NLCC_)%active) RETURN

    CALL mem%Add_var('act_fpc', mem%act_fpc,                                  &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('act_fpc', '-', 'actual foliage projective cover'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=MEDIUM,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('bare_fpc', mem%bare_fpc,                                &
      & hgrid, surface,                                                       &
      & t_cf('bare_fpc', '-', 'bare foliage projective cover'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('pot_fpc', mem%pot_fpc,                                  &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('pot_fpc', '-', 'potential foliage projective cover'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=MEDIUM,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var('cover_fract_pot', mem%cover_fract_pot,                  &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('cover_fract_pot', '-', 'potential cover fractions'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('bio_exist', mem%bio_exist,                              &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('bio_exist', '-', 'bio-climatic limits'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=MEDIUM,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('sum_green_bio_memory', mem%sum_green_bio_memory,        &
      & hgrid, surface,                                                       &
      & t_cf('sum_green_bio_memory', '-', 'sum_green_bio_memory'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('desert_fpc', mem%desert_fpc,                            &
      & hgrid, surface,                                                       &
      & t_cf('desert_fpc', '-', 'desert foliage projective cover'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    ! diagnostic output variables
    CALL mem%Add_var('cover_fract', mem%cover_fract,                          &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('cover_fract', '-', 'PFT cover fractions'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('npp_pft', mem%npp_pft,                                  &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('npp_pft', 'mol/m2/s', 'NPP per PFT'),                           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                    &
      & lrestart_cont = .TRUE.,                                               &
      & initval_r=0.0_wp )

  END SUBROUTINE Init_nlcc_memory

#endif
END MODULE mo_nlcc_memory_class
