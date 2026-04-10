!> Contains the memory class for the climbuf process.
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
MODULE mo_climbuf_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d

  dsl4jsb_Use_processes WLCC_, FLCC_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_climbuf_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_climbuf_memory

  TYPE(t_jsb_var_real2d) ::     &
    & q_rel_air_climbuf,        & ! Relative humidity of lowest air layer smoothed in time.
    & q_rel_air_climbuf_yday,   & ! Relative humidity of lowest air layer smoothed in time at the end of the last day.
                                  ! This is only needed for cbalone forcing!
    & prev_day_max_wind_10m,    & ! for windbreak
    & max_wind_10m,             & ! for windbreak
    & max_wind_10m_act,         & ! for windbreak
    & temp_sum_day, &
    & temp_sum_month, &
    & min_mmtemp_of_yr, &
    & max_mmtemp_of_yr, &
    & min_mmtemp20, &
    & max_mmtemp20, &
    & gdd_sum_year, &
    & gdd_prev_year


  CONTAINS
    PROCEDURE :: Init => Init_climbuf_memory
  END TYPE t_climbuf_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_memory_class'

CONTAINS

  SUBROUTINE Init_climbuf_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL, NONE
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_climbuf_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: lib_id     !< id of primary lct in lctlib
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid           ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface         ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_climbuf_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)
    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL mem%Add_var( 'q_air_clim', mem%q_rel_air_climbuf,                                            &
      & hgrid, surface,                                                                               &
      & t_cf('q_rel_air_climbuf', '-', 'Relative humiditiy of lowest air layer smoothed in time.'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                            &
      & prefix, suffix,                                                                               &
      & output_level=MEDIUM,                                                                          &
      & lrestart=.TRUE., initval_r=50.0_wp ) ! initval taken from JS3; lrestart=true only necessary if we restart within a day

    CALL mem%Add_var( 'q_air_clim_yday', mem%q_rel_air_climbuf_yday,                                    &
      & hgrid, surface,                                                                                 &
      & t_cf('q_rel_air_climbuf_yday', '-', 'Relative humiditiy of lowest air layer smoothed in time    &
                                               & at the end of the last day.'),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & lrestart=.FALSE., initval_r=50.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'prev_day_max_wind_10m', mem%prev_day_max_wind_10m,         &
      & hgrid, surface,                                                           &
      & t_cf('prev_day_max_wind_10m', '-', 'Previous day maximum wind speed.'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
      & prefix, suffix,                                                           &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var( 'max_wind_10m', mem%max_wind_10m,                                                &
      & hgrid, surface,                                                                                &
      & t_cf('max_wind_10m', '-', 'Maximum wind speed smoothed in time at the end of the last day..'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & output_level=MEDIUM,                                                                           &
      & lrestart=.TRUE., initval_r=0.0_wp )

    CALL mem%Add_var( 'max_wind_10m_act', mem%max_wind_10m_act,             &
      & hgrid, surface,                                                     &
      & t_cf('max_wind_10m_act', '-', 'This days maximum wind speed.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & output_level=MEDIUM,                                                &
      & lrestart=.TRUE., initval_r=0.0_wp )

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

  END SUBROUTINE Init_climbuf_memory

#endif
END MODULE mo_climbuf_memory_class
