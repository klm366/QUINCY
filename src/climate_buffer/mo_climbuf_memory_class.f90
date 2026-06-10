!> Contains the memory class for the climate buffer process.
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
!>### Memory definition of the climate buffer process
!>
!> The module defines the data structure [[t_climbuf_memory]] containing the climate buffer memory
!> Each tile of the CLIMBUF process has its own version of this memory. Variables added to the
!> structure are easily accessible in the process routines of CLIMBUF and other processes, e.g.
!> with 'dsl4jsb_Get_var2D_onChunk(CLIMBUF_, max_wind_10m)' for using variable 'max_wind_10m'.
!>
!> With the [[t_jsb_memory:Add_var]] calls, memory is allocated to the global variable arrays
!> on CPUs. It is also ensured that the correct variable patches are available on different MPI
!> processes, and the variable is made available on GPUs, in case the model is compiled to be
!> run on GPUs. Besides, calling [[t_jsb_memory:Add_var]] enables writing selected variables to
!> output files and adds the necessary variables to the list of restart file variables.
!>
!> Subroutine [[Init_climbuf_memory]] is a type-bound procedure of [[t_climbuf_memory]]. It is
!> called within the initialization phase for all tiles of the CLIMBUF process.
!>

MODULE mo_climbuf_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_climbuf_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 30   !< Maximum number of variables of this process

  !> Memory of the CLIMBUF process
  !> _Note for the ford documentation_: The table lists the general jsbach memory variables
  !>   ([[t_jsb_memory]]) as well as its extension by the CLIMBUF process variables.
  TYPE, EXTENDS(t_jsb_memory) :: t_climbuf_memory

  TYPE(t_jsb_var_real2d) ::     &
    & q_rel_air_14d,            & !< Relative humidity of lowest air layer smoothed in time
    & q_rel_air_14d_yday,       & !< Relative humidity of lowest air layer of the previous day
    & max_wind_yDay,            & !< Maximum 10m wind speed of the previous day [m/s]
    & max_wind_4yr,             & !< Maximum wind speed smoothed in time (~4 years) [m/s]
    & max_wind_of_day,          & !< Current maximum 10m wind speed of the day [m/s]
    & temp_sum_day,             & !< Temperature sum of the current day
    & temp_sum_month,           & !< Temperature sum of the current month
    & min_mmtemp_of_yr,         & !< Minimum monthly mean temperature of the year [K]
    & max_mmtemp_of_yr,         & !< Maximum monthly mean temperature of the year [K]
    & min_mmtemp_20yr,          & !< Smoothed minimum monthly mean temperature within 20 years [K]
    & max_mmtemp_20yr,          & !< Smoothed maximum monthly mean temperature within 20 years [K]
    & temp_15d,                 & !< Near surface temperature with 15-day smoothing [K]
    & temp_15d_yDay,            & !< 15-day temperature at the end of the previous day [K]
    & precip_15d,               & !< Precipitation with 15-day smoothing [kg m-2 s-1]
    & precip_15d_yDay,          & !< 15-day precipitation at the end of previous day [kg m-2 s-1]
    & gdd_sum_year,             & !< Current number of growing degree days in this year
    & gdd_prev_year,            & !< Growing degree days of the previous year
    & max_snow_year,            & !< Maximum snow depth in current year [m]
    & max_snow_prev_year,       & !< Maximum snow depth in previous year [m]
    & max_thaw_year,            & !< Maximum soil column thaw depth in current year [m]
    & max_thaw_10year             !< Smoothed maximum soil column thaw depth (~10 years) [m]

  TYPE(t_jsb_var_real3d) ::     &
    & t_soil_sl_15d,            & !< 15-day smoothed soil layer temperature [K]
    & t_soil_sl_15d_yDay,       & !< 15-day smoothed soil layer temperature at the end of the previous day [K]
    & wtr_soil_rel_sl_15d,      & !< 15-day smoothed relative soil layer moisture [-]
    & wtr_soil_rel_sl_15d_yDay    !< 15-day smoothed relative soil layer moisture at the end of the previous day [-]

  CONTAINS
    PROCEDURE :: Init => Init_climbuf_memory
  END TYPE t_climbuf_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_memory_class'

CONTAINS

  ! ============================================================================================== !
  !>
  !>#### Initialization of the climate buffer process memory
  !>
  !> In this subroutine we initialize the memory of the climate buffer process.
  !> The subroutine is a type-bound procedure of [[t_climbuf_memory]]. It is called within the
  !> initialization phase for all tiles needed by the climbuf process ([[mo_jsb_tile:init_tile]]:
  !> 'CALL this%mem(iproc)%p%Init(TRIM(process_name), ...)').
  !>
  SUBROUTINE Init_climbuf_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)

    USE mo_jsb_varlist,       ONLY: MEDIUM !, FULL, NONE, BASIC
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_climbuf_memory), INTENT(inout), TARGET :: mem !< Memory of the CLIMBUF process

    CHARACTER(len=*), INTENT(in) :: &
      & prefix,       & !< Prefix for variable names in output/restart files, i.e. process name
      & suffix          !< Suffix for variable names in output/restart files, i.e. corresponding tile

    INTEGER, INTENT(in) :: &
      & lct_ids(:),   & !< Land cover type IDs
      & lib_id,       & !< ID of the tile's primary land cover type (e.g. LAND_TYPE)
      & model_id        !< ID of the current model instance

    TYPE(t_jsb_grid),  POINTER :: hgrid     !< Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface   !< Vertical grid for surface variables (vertical dim. = 1)
    TYPE(t_jsb_vgrid), POINTER :: soil_e    !< Vertical grid for soil energy variables (vertical dim. = nsoil)
    TYPE(t_jsb_vgrid), POINTER :: soil_w    !< Vertical grid for soil hydrology variables (vertical dim. = nsoil)
    INTEGER :: table

    INTEGER  :: nsoil                        !< Number of soil layers
    REAL(wp) :: max_depth                    !< lower boundary of deepest soil layer

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_climbuf_memory'


    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    soil_e  => Get_vgrid('soil_depth_energy')
    soil_w  => Get_vgrid('soil_depth_water')

    nsoil = soil_e%n_levels
    max_depth = soil_e%ubounds(nsoil)

    CALL mem%Add_var( 'q_rel_air_14d', mem%q_rel_air_14d,                                 &
      & hgrid, surface,                                                                   &
      & t_cf('q_rel_air_14d', '%', 'Relative humiditiy of lowest air layer smoothed in time.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=50.0_wp )

    CALL mem%Add_var( 'q_rel_air_14d_yday', mem%q_rel_air_14d_yday,                       &
      & hgrid, surface,                                                                   &
      & t_cf('q_rel_air_14d_yday', '%',                                                   &
      & "Previous day's relative humiditiy of lowest air layer smoothed in time."),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.FALSE., initval_r=50.0_wp )

    CALL mem%Add_var( 'max_wind_yDay', mem%max_wind_yDay,                                 &
      & hgrid, surface,                                                                   &
      & t_cf('max_wind_yDay', 'm/s', 'Previous day maximum wind speed.'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var( 'max_wind_4yr', mem%max_wind_4yr,                                   &
      & hgrid, surface,                                                                   &
      & t_cf('max_wind_4yr', 'm/s', 'Maximum wind speed smoothed in time.'),              &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var( 'max_wind_of_day', mem%max_wind_of_day,                             &
      & hgrid, surface,                                                                   &
      & t_cf('max_wind_of_day', 'm/s', 'Current maximum wind speed of the day.'),         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var('temp_sum_day', mem%temp_sum_day,                                    &
      & hgrid, surface,                                                                   &
      & t_cf('temp_sum_day', 'Ks', 'Accumulated air temperature (to calculate average)'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var('temp_sum_month', mem%temp_sum_month,                                &
      & hgrid, surface,                                                                   &
      & t_cf('temp_sum_month', 'Ks', 'Sum of air temperature (to calculate average'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var('min_mmtemp_of_yr', mem%min_mmtemp_of_yr,                            &
      & hgrid, surface,                                                                   &
      & t_cf('min_mmtemp_of_yr', 'K', 'Current minimum monthly mean temperature of the year'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=1000.0_wp)

    CALL mem%Add_var('max_mmtemp_of_yr', mem%max_mmtemp_of_yr,                            &
      & hgrid, surface,                                                                   &
      & t_cf('max_mmtemp_of_yr', 'K', 'Current maximum monthly mean temperature of the year'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=-1000.0_wp)

    CALL mem%Add_var('min_mmtemp_20yr', mem%min_mmtemp_20yr,                              &
      & hgrid, surface,                                                                   &
      & t_cf('min_mmtemp_20yr', 'K', 'Climatological minimum monthly mean temperature'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=1000.0_wp)

    CALL mem%Add_var('max_mmtemp_20yr', mem%max_mmtemp_20yr,                              &
      & hgrid, surface,                                                                   &
      & t_cf('max_mmtemp_20yr', 'K', 'Climatological maximum monthly mean temperature'),  &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=-1000.0_wp)

    CALL mem%Add_var('temp_15d', mem%temp_15d,                                            &
      & hgrid, surface,                                                                   &
      & t_cf('temp_15d', 'K', 'Near surface Temperature with 15-day smoothing'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=285.0_wp)

    CALL mem%Add_var('temp_15d_yDay', mem%temp_15d_yDay,                                  &
      & hgrid, surface,                                                                   &
      & t_cf('temp_15d_yDay', 'K', '15-day smoothed temperature at end of previous day'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=285.0_wp)

    CALL mem%Add_var('precip_15d', mem%precip_15d,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('precip_15d', 'K', 'Precipitation (Rain+Snow) with 15-day smoothing'),       &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=3.0E-5_wp)

    CALL mem%Add_var('precip_15d_yDay', mem%precip_15d_yDay,                              &
      & hgrid, surface,                                                                   &
      & t_cf('precip_15d_yDay', 'K', '15-day smoothed temperature at end of previous day'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=3.0E-5_wp)

    ! Only calculate GDD on the vegetation tiles
    IF (lct_ids(1) == VEG_TYPE ) THEN
      CALL mem%Add_var('gdd_sum_year', mem%gdd_sum_year,                                  &
        & hgrid, surface,                                                                 &
        & t_cf('gdd_sum_year', 'degree', 'Current growing degree days of this year'),     &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                 &
        & lrestart=.TRUE., initval_r=3.0E-5_wp)

      CALL mem%Add_var('gdd_prev_year', mem%gdd_prev_year,                                &
        & hgrid, surface,                                                                 &
        & t_cf('gdd_prev_year', 'degree', 'Growing degree days of previous year'),        &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                 &
        & output_level=MEDIUM,                                                            &
        & lrestart=.TRUE., initval_r=0.0_wp)
    END IF

    CALL mem%Add_var('max_snow_year', mem%max_snow_year,                                  &
      & hgrid, surface,                                                                   &
      & t_cf('max_snow_year', 'm', 'Maximum snow depth in current year'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var('max_snow_prev_year', mem%max_snow_prev_year,                        &
      & hgrid, surface,                                                                   &
      & t_cf('max_snow_prev_year', 'm', 'Maximum snow depth in previous year'),           &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var('max_thaw_year', mem%max_thaw_year,                                  &
      & hgrid, surface,                                                                   &
      & t_cf('max_thaw_year', 'm', 'Maximum thaw depth in current year'),                 &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=max_depth)

    CALL mem%Add_var('max_thaw_10year', mem%max_thaw_10year,                              &
      & hgrid, surface,                                                                   &
      & t_cf('max_thaw_10year', 'm', 'Maximum annual thaw depth with 10-years smoothing'),&
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=max_depth)

    CALL mem%Add_var('t_soil_sl_15d', mem%t_soil_sl_15d,                                  &
      & hgrid, soil_e,                                                                    &
      & t_cf('soil_temperature_15d', 'K',                                                 &
      &      '15-day smoothed soil layer temperature'),                                   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=280.0_wp)

    CALL mem%Add_var('t_soil_sl_15d_yDay', mem%t_soil_sl_15d_yDay,                        &
      & hgrid, soil_e,                                                                    &
      & t_cf('soil_temperature_15d_yDay', 'K',                                            &
      &      '15-day smoothed soil layer temperature at the end of the previous day'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=280.0_wp)

    CALL mem%Add_var('wtr_soil_rel_sl_15d', mem%wtr_soil_rel_sl_15d,                      &
      & hgrid, soil_w,                                                                    &
      & t_cf('rel_soil_water_sl_15d', '-',                                                &
      &      '15-day smoothed rel soil layer water content'),                             &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & lrestart=.TRUE., initval_r=0.0_wp)

    CALL mem%Add_var('wtr_soil_rel_sl_15d_yDay', mem%wtr_soil_rel_sl_15d_yDay,            &
      & hgrid, soil_w,                                                                    &
      & t_cf('rel_soil_water_sl_15d_yDay', '-',                                           &
      &      '15-day smoothed rel soil layer water content at the end of the previous day'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & output_level=MEDIUM,                                                              &
      & lrestart=.TRUE., initval_r=0.0_wp)

  END SUBROUTINE Init_climbuf_memory

#endif
END MODULE mo_climbuf_memory_class
