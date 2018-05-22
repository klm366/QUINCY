!>
!! @brief Contains the memory class for the disturb process.
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2017-04-6)
!! 2018-04-26, StW: Added Spitfire support
!!
MODULE mo_disturb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  ! Use of processes in this module
  dsl4jsb_Use_processes DISTURB_

  ! Use process configurations
  dsl4jsb_Use_config(DISTURB_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_disturb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 50

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_disturb_memory

  TYPE(t_jsb_var_real2d) ::     &
    burned_fract,               & ! burned fraction of tile
    damaged_fract,              & ! windbreak fraction of tile
    ! Climbuf variables
    q_rel_air_climbuf,          & ! Relative humiditiy of lowest air layer smoothed in time.
    q_rel_air_climbuf_yday,     & ! Relative humiditiy of lowest air layer smoothed in time at the end of the last day.
                                  ! This is only needed for cbalone forcing!
    prev_day_max_wind_10m,      & ! for windbreak
    max_wind_10m,               & ! for windbreak
    max_wind_10m_act,           & ! for windbreak
    ! Spitfire variables (of which many should be on the veg tile and others both on veg and pft tiles)
    burned_fract_diag,          & ! burned fraction of tile excluding burned areas with surviving vegetation
    fract_1hr_wood,             &
    fract_10hr_wood,            &
    fract_100hr_wood,           &
    fract_1000hr_wood,          &
    fuel,                       &
    fuel_1hr,                   & ! Veg + pft?
    fuel_10hr,                  & ! Veg + pft?
    fuel_100hr,                 & ! Veg + pft?
    fuel_1000hr,                & ! Veg + pft?
    FDI,                        & ! Veg
    IR,                         & ! Veg
    ROS,                        & ! Veg
    fire_area,                  & ! Veg
    fire_duration,              & ! Veg
    fire_density,               & ! Veg?
    CombustCompleteness_lg,     &
    CombustCompleteness_1hr,    &
    CombustCompleteness_10hr,   &
    CombustCompleteness_100hr,  &
    CombustCompleteness_1000hr, &
    vegetation_height,          &
    fire_intensity,             &
    postfire_mortality,         &
    burned_c_litter_green_ag,   &
    burned_c_litter_wood_ag,    &
    burned_c_green,             &
    burned_c_reserve,           &
    burned_c_wood,              &
    burned_acid_ag1,            &
    burned_water_ag1,           &
    burned_ethanol_ag1,         &
    burned_nonsoluble_ag1,      &
    burned_acid_ag2,            &
    burned_water_ag2,           &
    burned_ethanol_ag2,         &
    burned_nonsoluble_ag2
  CONTAINS
    PROCEDURE :: Init => Init_disturb_memory
  END TYPE t_disturb_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_memory_class'

CONTAINS

  SUBROUTINE Init_disturb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,           ONLY: jsb_add_var, BASIC, MEDIUM, FULL
    USE mo_jsb_io,                ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_grid, Get_vgrid
    USE mo_disturbance_constants, ONLY: FIRE_NONE, FIRE_JSB3, FIRE_SPITFIRE

    CLASS(t_disturb_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    dsl4jsb_Def_config(DISTURB_)

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface, canopy_layer        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_disturb_memory'

    model => Get_model(model_id)

    dsl4jsb_Get_Config(DISTURB_)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')
    !canopy_layer => Get_vgrid('canopy_layer')

    IF (.NOT. ASSOCIATED(mem%burned_fract%ptr)) THEN

      CALL jsb_add_var( mem, 'burned_fract', mem%burned_fract,                            &
        & hgrid, surface,                                                                 &
        & t_cf('burned_fract', 'frac', 'Fraction of tile that is burned.'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
        & prefix, suffix,                                                                 &
        & output_level=BASIC,                                                             &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'damaged_fract', mem%damaged_fract,                           &
        & hgrid, surface,                                                                  &
        & t_cf('damaged_fract', 'frac', 'Fraction of tile that is damaged by windbreak.'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
        & prefix, suffix,                                                                  &
        & output_level=BASIC,                                                              &
        & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'q_air_clim', mem%q_rel_air_climbuf,                                       &
        & hgrid, surface,                                                                               &
        & t_cf('q_rel_air_climbuf', '-', 'Relative humiditiy of lowest air layer smoothed in time.'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                            &
        & prefix, suffix,                                                                               &
        & lrestart=.TRUE., initval_r=50.0_wp ) ! initval taken from JSBACH3; lrestart=true only necessary if we restart within a day

      CALL jsb_add_var( mem, 'q_air_clim_yday', mem%q_rel_air_climbuf_yday,                               &
        & hgrid, surface,                                                                                 &
        & t_cf('q_rel_air_climbuf_yday', '-', 'Relative humiditiy of lowest air layer smoothed in time    &
                                                 & at the end of the last day.'),                         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
        & prefix, suffix,                                                                                 &
        & output_level=BASIC,                                                                             &
        & lrestart=.FALSE., initval_r=50.0_wp ) ! initval taken from JSBACH3

      CALL jsb_add_var( mem, 'prev_day_max_wind_10m', mem%prev_day_max_wind_10m,    &
        & hgrid, surface,                                                           &
        & t_cf('prev_day_max_wind_10m', '-', 'Previous day maximum wind speed.'),   &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp )

      CALL jsb_add_var( mem, 'max_wind_10m', mem%max_wind_10m,                                           &
        & hgrid, surface,                                                                                &
        & t_cf('max_wind_10m', '-', 'Maximum wind speed smoothed in time at the end of the last day..'), &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
        & prefix, suffix,                                                                                &
        & lrestart=.TRUE., initval_r=0.0_wp ) 

      CALL jsb_add_var( mem, 'max_wind_10m_act', mem%max_wind_10m_act,              &
        & hgrid, surface,                                                           &
        & t_cf('max_wind_10m_act', '-', 'This days maximum wind speed.'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
        & prefix, suffix,                                                           &
        & lrestart=.FALSE., initval_r=0.0_wp )

      SELECT CASE (dsl4jsb_Config(DISTURB_)%fire_algorithm)
        CASE (FIRE_NONE)
          ! Per definition nothing to do
        CASE (FIRE_JSB3)
          ! JSB3-JSBACH-fire does not need any further variables
        CASE (FIRE_SPITFIRE)
          CALL jsb_add_var( mem, 'burned_fract_diag', mem%burned_fract_diag,     &
            & hgrid, surface,                                                    &
            & t_cf('burned_fract_diag', '-', 'fraction burned excluding area of surviving vegetation'),&
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fract_1hr_wood', mem%fract_1hr_wood,           &
            & hgrid, surface,                                                    &
            & t_cf('fract_1hr_wood', '-', 'fraction of wood of fuel class 1hr.'),&
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fract_10hr_wood', mem%fract_10hr_wood,            &
            & hgrid, surface,                                                       &
            & t_cf('fract_10hr_wood', '-', 'fraction of wood of fuel class 10hr.'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fract_100hr_wood', mem%fract_100hr_wood,            &
            & hgrid, surface,                                                         &
            & t_cf('fract_100hr_wood', '-', 'fraction of wood of fuel class 100hr.'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),      &
            & prefix, suffix,                                                         &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fract_1000hr_wood', mem%fract_1000hr_wood,            &
            & hgrid, surface,                                                           &
            & t_cf('fract_1000hr_wood', '-', 'fraction of wood of fuel class 1000hr.'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
            & prefix, suffix,                                                           &
            & lrestart=.FALSE., initval_r=0.0_wp)
          CALL jsb_add_var( mem, 'fuel', mem%fuel,                               &
            & hgrid, surface,                                                    &
            & t_cf('fuel', 'mol m^-2', 'carbon in pft tile fuel.'),              &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fuel_1hr', mem%fuel_1hr,                       &
            & hgrid, surface,                                                    &
            & t_cf('fuel_1hr', 'mol m^-2', 'carbon in fuel class 1hr.'),         &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fuel_10hr', mem%fuel_10hr,                     &
            & hgrid, surface,                                                    &
            & t_cf('fuel_10hr', 'mol m^-2', 'carbon in fuel class 10hr.'),       &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fuel_100hr', mem%fuel_100hr,                   &
            & hgrid, surface,                                                    &
            & t_cf('fuel_100hr', 'mol m^-2', 'carbon in fuel class 100hr.'),     &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fuel_1000hr', mem%fuel_1000hr,                 &
            & hgrid, surface,                                                    &
            & t_cf('fuel_1000hr', 'mol m^-2', 'carbon in fuel class 1000hr.'),   &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'FDI', mem%FDI,                                 &
            & hgrid, surface,                                                    &
            & t_cf('FDI', '-', 'Fire danger index.'),                            &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'IR', mem%IR,                                   &
            & hgrid, surface,                                                    &
            & t_cf('IR', 'kJ m^-2 min^-1', 'Fire reaction intensity.'),          &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'ROS', mem%ROS,                                 &
            & hgrid, surface,                                                    &
            & t_cf('ROS', 'm min^-1','Fire rate of spread.'),                    &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fire_area', mem%fire_area,                     &
            & hgrid, surface,                                                    &
            & t_cf('fire_area', 'm^2','Size of a single fire.'),                 &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fire_duration', mem%fire_duration,             &
            & hgrid, surface,                                                    &
            & t_cf('fire_duration', 'min','Duration of a fire.'),                &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fire_density', mem%fire_density,               &
            & hgrid, surface,                                                    &
            & t_cf('fire_density', 'm^-2','Number of fires per area.'),          &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'CombustCompleteness_lg', mem%CombustCompleteness_lg,           &
            & hgrid, surface,                                                                    &
            & t_cf('CombustCompleteness_lg', '-','Combustion completeness of living grass'),     &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
            & prefix, suffix,                                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'CombustCompleteness_1hr', mem%CombustCompleteness_1hr,         &
            & hgrid, surface,                                                                    &
            & t_cf('CombustCompleteness_1hr', '-','Combustion completeness of 1hr fuels'),       &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
            & prefix, suffix,                                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'CombustCompleteness_10hr', mem%CombustCompleteness_10hr,       &
            & hgrid, surface,                                                                    &
            & t_cf('CombustCompleteness_10hr', '-','Combustion completeness of 10hr fuels'),     &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
            & prefix, suffix,                                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'CombustCompleteness_100hr', mem%CombustCompleteness_100hr,     &
            & hgrid, surface,                                                                    &
            & t_cf('CombustCompleteness_100hr', '-','Combustion completeness of 100hr fuels'),   &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
            & prefix, suffix,                                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'CombustCompleteness_1000hr', mem%CombustCompleteness_1000hr,   &
            & hgrid, surface,                                                                    &
            & t_cf('CombustCompleteness_1000hr', '-','Combustion completeness of 1000hr fuels'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
            & prefix, suffix,                                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'vegetation_height', mem%vegetation_height,        &
            & hgrid, surface,                                                       &
            & t_cf('vegetation_height', 'm','Characteristic height of vegetation'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'fire_intensity', mem%fire_intensity,           &
            & hgrid, surface,                                                    &
            & t_cf('fire_intensity', 'kW m^-1','Energy release by fires'),       &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'postfire_mortality', mem%postfire_mortality,   &
            & hgrid, surface,                                                    &
            & t_cf('postfire_mortality', '','Fraction of BA vegetation killed'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_c_litter_green_ag', mem%burned_c_litter_green_ag, &
            & hgrid, surface,                                                              &
            & t_cf('burned_c_litter_green_ag', 'mol m^-2','burned litter green carbon'),   &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
            & prefix, suffix,                                                              &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_c_litter_wood_ag', mem%burned_c_litter_wood_ag, &
            & hgrid, surface,                                                            &
            & t_cf('burned_c_litter_wood_ag', 'mol m^-2','burned litter wood carbon'),   &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),         &
            & prefix, suffix,                                                            &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_c_green', mem%burned_c_green,              &
            & hgrid, surface,                                                       &
            & t_cf('burned_c_green', 'mol m^-2','burned green carbon'),             &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_c_reserve', mem%burned_c_reserve,       &
            & hgrid, surface,                                                    &
            & t_cf('burned_c_reserve', 'mol m^-2','burned reserve carbon'),      &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_c_wood', mem%burned_c_wood,             &
            & hgrid, surface,                                                    &
            & t_cf('burned_c_wood', 'mol m^-2','burned_wood carbon'),            &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits), &
            & prefix, suffix,                                                    &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_acid_ag1', mem%burned_acid_ag1,            &
            & hgrid, surface,                                                       &
            & t_cf('burned_acid_ag1', 'mol m^-2','burned acid_ag1 carbon'),         &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_water_ag1', mem%burned_water_ag1,          &
            & hgrid, surface,                                                       &
            & t_cf('burned_water_ag1', 'mol m^-2','burned water_ag1 carbon'),       &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_ethanol_ag1', mem%burned_ethanol_ag1,      &
            & hgrid, surface,                                                       &
            & t_cf('burned_ethanol_ag1', 'mol m^-2','burned ethanol_ag1 carbon'),   &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_nonsoluble_ag1', mem%burned_nonsoluble_ag1,    &
            & hgrid, surface,                                                           &
            & t_cf('burned_nonsoluble_ag1', 'mol m^-2','burned nonsoluble_ag1 carbon'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
            & prefix, suffix,                                                           &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_acid_ag2', mem%burned_acid_ag2,            &
            & hgrid, surface,                                                       &
            & t_cf('burned_acid_ag2', 'mol m^-2','burned acid_ag2 carbon'),         &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_water_ag2', mem%burned_water_ag2,          &
            & hgrid, surface,                                                       &
            & t_cf('burned_water_ag2', 'mol m^-2','burned water_ag2 carbon'),       &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_ethanol_ag2', mem%burned_ethanol_ag2,      &
            & hgrid, surface,                                                       &
            & t_cf('burned_ethanol_ag2', 'mol m^-2','burned ethanol_ag2 carbon'),   &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),    &
            & prefix, suffix,                                                       &
            & lrestart=.FALSE., initval_r=0.0_wp )
          CALL jsb_add_var( mem, 'burned_nonsoluble_ag2', mem%burned_nonsoluble_ag2,    &
            & hgrid, surface,                                                           &
            & t_cf('burned_nonsoluble_ag2', 'mol m^-2','burned nonsoluble_ag2 carbon'), &
            & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
            & prefix, suffix,                                                           &
            & lrestart=.FALSE., initval_r=0.0_wp )
      END SELECT

    END IF

  END SUBROUTINE Init_disturb_memory

#endif
END MODULE mo_disturb_memory_class
