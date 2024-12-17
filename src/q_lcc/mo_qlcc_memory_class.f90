!> Contains the memory class for the qlcc process.
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_qlcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: !X e.g: LAND_TYPE, LAKE_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: NONE, BASIC, MEDIUM, FULL
  USE mo_jsb_physical_constants, ONLY: !X e.g: tmelt
  USE mo_exception,       ONLY: message

  ! Use of prcesses in this module
  dsl4jsb_Use_processes QLCC_

  ! Use process configurations
  dsl4jsb_Use_config(QLCC_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_qlcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_qlcc_memory

    ! Common variables independent of lct_type
    TYPE(t_jsb_var_real2d) :: &
      seconds_day,            &
      seconds_month,          &
      temp_sum_day,           &
      temp_sum_month,         &
      min_mmtemp_of_yr,       &
      max_mmtemp_of_yr,       &
      min_mmtemp20,           &
      max_mmtemp20,           &
      gdd_sum_year,           &
      gdd_prev_year,          &
      bare_fpc,               &
      sum_green_bio_memory,   &
      desert_fpc,             &
      tree_fract,             & ! land fraction covered with tree PFTs
      shrub_fract,            & ! land fraction covered with shrub PFTs
      grass_fract,            & ! land fraction covered with grass PFTs
      crop_fract,             & ! land fraction covered with crop PFTs
      pasture_fract,          & ! land fraction covered with pasture PFTs
      C3pft_fract,            & ! land fraction covered with C3 PFTs
      C4pft_fract,            & ! land fraction covered with C4 PFTs
      baresoil_fract            ! bare land fraction

    TYPE(t_jsb_var_real3d) :: &
      act_fpc,                &
      pot_fpc,                &
      cover_fract_pot,        &
      cover_fract,            &
      bio_exist

    ! Diagnostic global sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) :: &
      tree_area_gsum,         & ! global land area covered with tree PFTs
      shrub_area_gsum,        & ! global land area covered with shrub PFTs
      grass_area_gsum,        & ! global land area covered with grass PFTs
      crop_area_gsum,         & ! global land area covered with crop PFTs
      pasture_area_gsum,      & ! global land area covered with pasture PFTs
      C3pft_area_gsum,        & ! global land area covered with C3 PFTs
      C4pft_area_gsum,        & ! global land area covered with C4 PFTs
      baresoil_area_gsum        ! global bare land area
    ! Additional variables for LAND lct_type
!!$  TR  TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & s_star,               & !< Surface dry static energy (s^star, see manual)    [m2 s-2]


    ! Additional variables for LAKE lct_type
!!$ TR    TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & t_wtr,                & !< Lake surface temperature (water)     [K]

  CONTAINS
    PROCEDURE :: Init => Init_qlcc_memory
  END TYPE t_qlcc_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qlcc_memory_class'

CONTAINS

  SUBROUTINE Init_qlcc_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_qlcc_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, pft_vgrid           ! Vertical grid
    INTEGER :: table
    INTEGER :: output_level

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_qlcc_memory'

    CALL message(TRIM(routine), 'Starting up ...')

    model => Get_model(model_id)

    table = tables(1)

    hgrid     => Get_grid(mem%grid_id)
    surface   => Get_vgrid('surface')
    pft_vgrid => Get_vgrid('pfts')

    CALL mem%Add_var('seconds_day', mem%seconds_day,                          &
      & hgrid, surface,                                                       &
      & t_cf('seconds_day', 'sum of seconds in day', 'time sum day'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('seconds_month', mem%seconds_month,                      &
      & hgrid, surface,                                                       &
      & t_cf('seconds_month', 'sum of seconds in month', 'time sum month'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('temp_sum_day', mem%temp_sum_day,                        &
      & hgrid, surface,                                                       &
      & t_cf('temp_sum_day', 'sum of air temperature in day', 'temperature sum day'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('temp_sum_month', mem%temp_sum_month,                    &
      & hgrid, surface,                                                       &
      & t_cf('temp_sum_month', 'sum of air temperature in month', 'temperature sum month'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('min_mmtemp_of_yr', mem%min_mmtemp_of_yr,                &
      & hgrid, surface,                                                       &
      & t_cf('min_mmtemp_of_yr', 'min temperature of year', 'temperature min'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=1000.0_wp )

    CALL mem%Add_var('max_mmtemp_of_yr', mem%max_mmtemp_of_yr,                &
      & hgrid, surface,                                                       &
      & t_cf('max_mmtemp_of_yr', 'max temperature of year', 'temperature max'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=-1000.0_wp )

    CALL mem%Add_var('min_mmtemp20', mem%min_mmtemp20,                        &
      & hgrid, surface,                                                       &
      & t_cf('min_mmtemp20', 'min temperature climatology', 'temp min climate'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                   &
      & initval_r=1000.0_wp )

    CALL mem%Add_var('max_mmtemp20', mem%max_mmtemp20,                        &
      & hgrid, surface,                                                       &
      & t_cf('max_mmtemp20', 'max temperature climatology', 'temp max climate'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                   &
      & initval_r=-1000.0_wp )

    CALL mem%Add_var('gdd_sum_year', mem%gdd_sum_year,                        &
      & hgrid, surface,                                                       &
      & t_cf('gdd_sum_year', 'gdd current year', 'gdd sum year'),             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                    &
      & initval_r=0.0_wp )

    CALL mem%Add_var('gdd_prev_year', mem%gdd_prev_year,                      &
      & hgrid, surface,                                                       &
      & t_cf('gdd_prev_year', 'gdd previous year', 'gdd prev year'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=FULL,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('act_fpc', mem%act_fpc,                                  &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('act_fpc', 'actual foliage projective cover', 'act fpc'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=MEDIUM,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('bare_fpc', mem%bare_fpc,                                &
      & hgrid, surface,                                                       &
      & t_cf('bare_fpc', 'bare foliage projective cover', 'bare fpc'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('pot_fpc', mem%pot_fpc,                                  &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('pot_fpc', 'potential foliage projective cover', 'pot fpc'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=MEDIUM,                                                  &
      & initval_r=0.0_wp )

    CALL mem%Add_var('cover_fract_pot', mem%cover_fract_pot,                  &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('cover_fract_pot', 'potential cover fractions', 'cover fract pot'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('bio_exist', mem%bio_exist,                              &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('bio_exist', 'bio-climatic limits', 'bio exist'),                &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=MEDIUM,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('sum_green_bio_memory', mem%sum_green_bio_memory,        &
      & hgrid, surface,                                                       &
      & t_cf('sum_green_bio_memory', 'desert memory', 'sum_green_bio_memory'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=NONE,                                                   &
      & initval_r=0.0_wp )

    CALL mem%Add_var('desert_fpc', mem%desert_fpc,                            &
      & hgrid, surface,                                                       &
      & t_cf('desert_fpc', 'desert foliage projective cover', 'desert fpc'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    ! diagnostic output variables
    CALL mem%Add_var('cover_fract', mem%cover_fract,                          &
      & hgrid, pft_vgrid,                                                     &
      & t_cf('cover_fract', 'PFT cover fractions', 'cover fract'),            &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
      & prefix, suffix,                                                       &
      & output_level=BASIC,                                                   &
      & initval_r=0.0_wp )

    ! disagnostics from pplcc
!     IF ((model%Is_process_enabled(QLCC_) .OR. model%Is_process_enabled(ALCC_)) &
!         .AND. TRIM(suffix) == 'box') THEN
      output_level = BASIC   ! Assign variable to basic output group
!     ELSE
!       output_level = NONE    ! Do not assigne variable to an output group
!     END IF

    ! Area fractions of selected groups of PFTs - only calculated on box tile
    IF (ASSOCIATED(mem%parent)) RETURN

    CALL mem%Add_var('tree_fract', mem%tree_fract,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('tree_fract', '1', 'Tree fraction relative to the grid cell area'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    CALL mem%Add_var('shrub_fract', mem%shrub_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('shrub_fract', '1', 'Shrub fraction relative to the grid cell area'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    CALL mem%Add_var('grass_fract', mem%grass_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('grass_fract', '1', 'Grass fraction relative to the grid cell area'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    CALL mem%Add_var('crop_fract', mem%crop_fract,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('crop_fract', '1', 'Crop fraction relative to the grid cell area'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    CALL mem%Add_var('pasture_fract', mem%pasture_fract,                                  &
      & hgrid, surface,                                                                   &
      & t_cf('pasture_fract', '1', 'Pasture fraction relative to the grid cell area'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    CALL mem%Add_var('C3pft_fract', mem%C3pft_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('C3pft_fract', '1', 'Fraction of C3-PFTs relative to the grid cell area'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    CALL mem%Add_var('C4pft_fract', mem%C4pft_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('C4pft_fract', '1', 'Fraction of C4-PFTs relative to the grid cell area'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    CALL mem%Add_var('baresoil_fract', mem%baresoil_fract,                                &
      & hgrid, surface,                                                                   &
      & t_cf('baresoil_fract', '1', 'Bare soil fraction relative to the grid cell area'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=output_level, initval_r=0.0_wp)

    ! diagnostics from pplcc
    ! Diagnostic global areas for experiment monitoring
    CALL mem%Add_var('tree_area_gsum', mem%tree_area_gsum,                     &
      & hgrid, surface,                                                        &
      & t_cf('tree_area_gsum', 'Mio km2', 'Global tree area'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('shrub_area_gsum', mem%shrub_area_gsum,                   &
      & hgrid, surface,                                                        &
      & t_cf('shrub_area_gsum', 'Mio km2', 'Global shrub area'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('grass_area_gsum', mem%grass_area_gsum,                   &
      & hgrid, surface,                                                        &
      & t_cf('grass_area_gsum', 'Mio km2', 'Global grass area'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('crop_area_gsum', mem%crop_area_gsum,                     &
      & hgrid, surface,                                                        &
      & t_cf('crop_area_gsum', 'Mio km2', 'Global crop area'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('pasture_area_gsum', mem%pasture_area_gsum,               &
      & hgrid, surface,                                                        &
      & t_cf('pasture_area_gsum', 'Mio km2', 'Global pasture area'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('baresoil_area_gsum', mem%baresoil_area_gsum,             &
      & hgrid, surface,                                                        &
      & t_cf('baresoil_area_gsum', 'Mio km2', 'Global baresoil area'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('C3pft_area_gsum', mem%C3pft_area_gsum,                   &
      & hgrid, surface,                                                        &
      & t_cf('C3pft_area_gsum', 'Mio km2', 'Global area of C3 PFTs'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)
    CALL mem%Add_var('C4pft_area_gsum', mem%C4pft_area_gsum,                   &
      & hgrid, surface,                                                        &
      & t_cf('C4pft_area_gsum', 'Mio km2', 'Global area of C4 PFTs'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
      & prefix, suffix,                                                        &
      & lrestart=.FALSE., initval_r=0.0_wp)

  END SUBROUTINE Init_qlcc_memory

#endif
END MODULE mo_qlcc_memory_class
