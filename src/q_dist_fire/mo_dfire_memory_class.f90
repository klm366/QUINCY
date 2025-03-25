!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> 
!> @brief Contains the memory class for the dist_fire process.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!>
MODULE mo_dfire_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  USE mo_dfire_constants,       ONLY: df_sf_vgrid_dead_fuel_pools_name, &
                                      df_sf_vgrid_live_fuel_pools_name, &
                                      df_sf_vgrid_fuel_classes_name

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_dfire_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 80

  ! ======================================================================================================= !
  !>
  !> Type definition for dfire memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_dfire_memory
    TYPE(t_jsb_var_real2d)               :: t_air_max           !< Aproximation to daily max temp (from prev day) [k]
    TYPE(t_jsb_var_real2d)               :: t_air_min           !< Aproximation to daily min temp (from prev day) [k]

    TYPE(t_jsb_var_real2d)               :: pot_burn_area       !< Potential burnable area (forest + grass area only) [m2?]

    TYPE(t_jsb_var_real2d)               :: mask_is_fire_intense_enough !< Mask (LOGICAL): TRUE=current fire is intense enough to exist [-]
                                                                !< Conversion required from REAL(wp) <-> LOGICAL
                                                                !< True = 1._wp, False = 0._wp

    ! dead fuel
    TYPE(t_jsb_var_real3d)               :: dead_fuel_consumption_1hr !< [g DM m-2]
    TYPE(t_jsb_var_real3d)               :: dead_fuel_consumption_10hr !< [g DM m-2]
    TYPE(t_jsb_var_real3d)               :: dead_fuel_consumption_100hr !< [g DM m-2]
    TYPE(t_jsb_var_real3d)               :: dead_fuel_consumption_1000hr !< [g DM m-2]
    ! live fuel
    TYPE(t_jsb_var_real3d)               :: live_fuel_consumption_1hr !< [g DM m-2]
    !TYPE(t_jsb_var_real3d)               :: live_fuel_consumption_10hr !< [g DM m-2]
    !TYPE(t_jsb_var_real3d)               :: live_fuel_consumption_100hr !< [g DM m-2]
    !TYPE(t_jsb_var_real3d)               :: live_fuel_consumption_1000hr !< [g DM m-2]

    TYPE(t_jsb_var_real3d)               :: dead_frac_consumption !< [-]

    TYPE(t_jsb_var_real2d)               :: area_burned_fract   !< Fraction of the area burned [%]
    TYPE(t_jsb_var_real2d)               :: area_burned         !< Area burned in a gridcell [ha ts-1]m-2]
    TYPE(t_jsb_var_real2d)               :: i_surface           !< Fire intensity [kW m-1]

    !< Precipitation
    !< Multiple methods used to aproximate daily precipitation. Not directly available into IQ.
    TYPE(t_jsb_var_real2d)               :: precip_qc_mavg      !< Precipitation moving average (quincy moving average) [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_daytime      !< Aproximation to daily precipitation (previous daily data, aproximatted by sun movement) [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_mavg         !< Precipitation moving average (Pure implementation of moving average) [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_daily        !< Daily precipitation, it links to precip_daytime, precip_mavg or precip_qc_mavg [kg m-2 s-1]
    TYPE(t_jsb_var_real2d)               :: precip_daily_mmd    !< Daily precipitation, same as precip_daily but in spitfire expected units [mm day-1]

    !< Diagnostics
    TYPE(t_jsb_var_real2d)               :: nesterov_index      !< Cumulative function of tmax and tdew [C2]
    TYPE(t_jsb_var_real2d)               :: fdi                 !< fire danger index [%]
    TYPE(t_jsb_var_real2d)               :: pm_fire             !< Probability of mortality
    TYPE(t_jsb_var_real2d)               :: pm_ck               !< Probability of mortality due to crown damage
    TYPE(t_jsb_var_real2d)               :: pm_tau              !< Probability of mortality due to cambial damage
!    TYPE(t_jsb_var_real2d)               :: tau_l               !< Residence time of the fire [min]
!    TYPE(t_jsb_var_real2d)               :: tau_c               !< Critical time for cambial damage [min]
!    TYPE(t_jsb_var_real2d)               :: bark_thickness      !< Bark thickness [cm]
!    TYPE(t_jsb_var_real2d)               :: DBH                 !< Diameter at breast height [cm]
!    TYPE(t_jsb_var_real2d)               :: par1                !< Parameters for woody PFTs used in bark thickness calculation [-]
!    TYPE(t_jsb_var_real2d)               :: par2                !< Parameters for woody PFTs used in bark thickness calculation [-]
    TYPE(t_jsb_var_real2d)               :: ros_f               !< Forward rate of spread of a surface fire [m min-1]
!    TYPE(t_jsb_var_real2d)               :: ros_b               !< Backwards rate of spread of a surface fire [m min-1]
    TYPE(t_jsb_var_real2d)               :: mean_fire_area      !< Mean fire area [ha]
    TYPE(t_jsb_var_real2d)               :: fire_duration       !< Fire duration [min]
!    TYPE(t_jsb_var_real2d)               :: distance_travelled  !< Total distance travelled by fire [m]
!    TYPE(t_jsb_var_real2d)               :: l2b_ratio           !< Length to breadth ratio [-]

    TYPE(t_jsb_var_real2d)               :: e_n_ignitions                 !< (eq 2) Sum of independent estimates of n_l,ig and n_h,ig
    !<        ignitions events [ha-1 ts-1]
    TYPE(t_jsb_var_real2d)               :: e_n_humans                    !< Rate of human caused ignitions [ha-1 ts-1]
    TYPE(t_jsb_var_real2d)               :: e_n_lightning                 !< Rate of lightning-caused ignitions [ha-1 ts-1]
  CONTAINS
    PROCEDURE :: Init => Init_dfire_memory
  END TYPE t_dfire_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> initialize memory (variables) for the process: dist_fire
  !>
  SUBROUTINE Init_dfire_memory(mem, prefix, suffix, lct_ids, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_dfire_memory), INTENT(inout), TARGET  :: mem           !> dfire memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                            !> Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                    !> Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_dead_fuel_pools            !> Dead fuel pools grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_live_fuel_pools            !> Live fuel pools grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_fuel_classes               !> Fuel classes grid
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_dfire_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(model_id)
    table         =  tables(1)
    hgrid         => Get_grid(model%grid_id)
    vgrid_surface => Get_vgrid('surface')
    vgrid_dead_fuel_pools => Get_vgrid(df_sf_vgrid_dead_fuel_pools_name)
    vgrid_live_fuel_pools => Get_vgrid(df_sf_vgrid_live_fuel_pools_name)
    vgrid_fuel_classes => Get_vgrid(df_sf_vgrid_fuel_classes_name)
    ! ----------------------------------------------------------------------------------------------------- !
    CALL mem%Add_var('nesterov_index', mem%nesterov_index, &
        & hgrid, vgrid_surface, &
        & t_cf('nesterov_index', 'Celsius^2', 'Cumulative function of daily tmax and tdew '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fdi', mem%fdi, &
        & hgrid, vgrid_surface, &
        & t_cf('fdi', '%', 'Fire danger index'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('t_air_max', mem%t_air_max, &
        & hgrid, vgrid_surface, &
        & t_cf('t_air_max', 'K', 'Max temperature of the previous day (daytime_counter)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('t_air_min', mem%t_air_min, &
        & hgrid, vgrid_surface, &
        & t_cf('t_air_min', 'K', 'Min temperature of the previous day (daytime_counter)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('mask_is_fire_intense_enough', mem%mask_is_fire_intense_enough, &
        & hgrid, vgrid_surface, &
        & t_cf('mask_is_fire_instense_enough', '-', 'True when fire is intense enough to exist'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    !< Precipitation related variables
    CALL mem%Add_var('precip_daily', mem%precip_daily, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daily', 'kg m-2 s-1', 'Daily precipitation (precip_mavg, precip_qc_mavg, precip_daytime'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_mavg', mem%precip_mavg, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_mavg', 'kg m-2 s-1', 'Precipitation moving average using pure method'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_qc_mavg', mem%precip_qc_mavg, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_qc_mavg', 'kg m-2 s-1', 'Precipitation moving average using quincy method'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_daily_mmd', mem%precip_daily_mmd, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daily_mmd', 'mm/day', 'Daily precipitation in Spitfire units'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_daytime', mem%precip_daytime, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daytime', 'kg m-2 s-1', 'Daily precipitation from previous day, related to sun movement'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pot_burn_area', mem%pot_burn_area, &
        & hgrid, vgrid_surface, &
        & t_cf('pot_burn_area', 'm2', 'Potential burnable area (forest and grass area)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_1hr', mem%dead_fuel_consumption_1hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_1hr', 'g DM m-2', '1 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_10hr', mem%dead_fuel_consumption_10hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_10hr', 'g DM m-2', '10 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_100hr', mem%dead_fuel_consumption_100hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_100hr', 'g DM m-2', '100 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_1000hr', mem%dead_fuel_consumption_1000hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_1000hr', 'g DM m-2', '1000 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('live_fuel_consumption_1hr', mem%live_fuel_consumption_1hr, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_fuel_consumption_1hr', 'g DM m-2', '1 hour live fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_frac_consumption', mem%dead_frac_consumption, &
        & hgrid, vgrid_fuel_classes, &
        & t_cf('dead_frac_consumption', '-', 'Dead fuel fraction of consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('area_burned_fract', mem%area_burned_fract, &
        & hgrid, vgrid_surface, &
        & t_cf('area_burned_fract', '%', 'Fraction of the area burned'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('area_burned', mem%area_burned, &
        & hgrid, vgrid_surface, &
        & t_cf('area_burned', 'ha ts-1', 'Area burned in a gridcell'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('i_surface', mem%i_surface, &
        & hgrid, vgrid_surface, &
        & t_cf('i_surface', 'kW m-1', 'Fire intensity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASI, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pm_fire', mem%pm_fire, &
        & hgrid, vgrid_surface, &
        & t_cf('pm_fire', '%', 'Probability of post-fire mortality'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pm_tau', mem%pm_tau, &
        & hgrid, vgrid_surface, &
        & t_cf('pm_tau', '%', 'Probability of post-fire mortality due to cambial damage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pm_ck', mem%pm_ck, &
        & hgrid, vgrid_surface, &
        & t_cf('pm_ck', '%', 'Probability of post-fire mortality due to crown damage'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ros_f', mem%ros_f, &
        & hgrid, vgrid_surface, &
        & t_cf('ros_f', 'm min-1', 'Forward rate of spread of a surface fire'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('mean_fire_area', mem%mean_fire_area, &
        & hgrid, vgrid_surface, &
        & t_cf('mean_fire_area', 'ha', 'Mean fire area'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fire_duration', mem%fire_duration, &
        & hgrid, vgrid_surface, &
        & t_cf('fire_duration', 'min', 'Fire duration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('e_n_ignitions', mem%e_n_ignitions, &
        & hgrid, vgrid_surface, &
        & t_cf('e_n_ignitions', 'ha-1 ts-1', 'Sum of independent estimates of n_l,ig and n_h,ig ignitions events'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('e_n_humans', mem%e_n_humans, &
        & hgrid, vgrid_surface, &
        & t_cf('e_n_humans', 'ha-1 ts-1', 'Rate of human-caused ignitions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('e_n_lightning', mem%e_n_lightning, &
        & hgrid, vgrid_surface, &
        & t_cf('e_n_lightning', 'ha-1 ts-1', 'Rate of lightning-caused ignitions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
  END SUBROUTINE Init_dfire_memory

#endif
END MODULE mo_dfire_memory_class
