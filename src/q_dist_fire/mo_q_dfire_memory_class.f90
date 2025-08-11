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
MODULE mo_q_dfire_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  USE mo_q_dfire_constants,      ONLY: df_sf_vgrid_dead_fuel_pools_name, &
                                      df_sf_vgrid_live_fuel_pools_name, &
                                      df_sf_vgrid_fuel_classes_name, &
                                      df_sf_vgrid_fuel_moisture_name, &
                                      df_sf_vgrid_months_in_year_name

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: t_q_dfire_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 80

  ! ======================================================================================================= !
  !>
  !> Type definition for dfire memory
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_dfire_memory
    TYPE(t_jsb_var_real3d)               :: fuel_moisture       !< moisture in fuel [0 to 1]
    ! diagnostic vars related to fuel moisture
    TYPE(t_jsb_var_real2d)               :: ab_moisture_diag    !< Soil moisture of the top layers (average of the 1st and 2nd)

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

    ! diagnostic for ros of spread calculation
    TYPE(t_jsb_var_real2d)               :: reaction_intensity             !< I_R: the energy release rate per unit area of fire front [kJ m−2 min−1]
    TYPE(t_jsb_var_real2d)               :: prop_flux_ratio                !< the proportion of IR that heats adjacent fuel particles to ignition [0 to 1]
    TYPE(t_jsb_var_real2d)               :: wind_factor                    !< wind_factor is a multiplier that accounts for the effect of wind in increasing the effective value of prop_flux_ratio[-]
    TYPE(t_jsb_var_real2d)               :: eff_heating_number             !< the proportion of a fuel particle that is heated to ...
                                                                           !< ignition temperature at the time flaming combustion starts [-]
    TYPE(t_jsb_var_real2d)               :: preig_heat                     !< the amount of heat required to ignite a given mass of fuel [kJ kg−1]

    !< diagnostics: for reaction_intensity
    TYPE(t_jsb_var_real2d)               :: opt_react_velo                 !< Optimum reaction velocity [min-1]
    TYPE(t_jsb_var_real2d)               :: net_fuel_load                  !< Amount of fuel after substraction of its mineral content [kg m-2]
    TYPE(t_jsb_var_real2d)               :: moist_damp_coeff               !< Moisture-dampening coefficient [-]
    TYPE(t_jsb_var_real2d)               :: mineral_damp_coeff             !< Mineral dampening coefficient [-]

    !< Diagnostics
    TYPE(t_jsb_var_real2d)               :: nesterov_index      !< Cumulative function of tmax and tdew [C2]
    TYPE(t_jsb_var_real2d)               :: fdi                 !< fire danger index [%]
    TYPE(t_jsb_var_real2d)               :: ros_f               !< Forward rate of spread of a surface fire [m min-1]
    TYPE(t_jsb_var_real2d)               :: ros_b               !< Backwards rate of spread of a surface fire [m min-1]
    TYPE(t_jsb_var_real2d)               :: mean_fire_area      !< Mean fire area [ha]
    TYPE(t_jsb_var_real2d)               :: fire_duration       !< Fire duration [min]
    TYPE(t_jsb_var_real2d)               :: distance_travelled  !< Total distance travelled by fire [m]
    TYPE(t_jsb_var_real2d)               :: l2b_ratio           !< Length to breadth ratio [-]

    TYPE(t_jsb_var_real2d)               :: e_n_ignitions       !< (eq 2) Sum of independent estimates of n_l,ig and n_h,ig
    !<        ignitions events [ha-1 ts-1]
    TYPE(t_jsb_var_real2d)               :: e_n_humans          !< Rate of human caused ignitions [ha-1 ts-1]
    TYPE(t_jsb_var_real2d)               :: e_n_lightning       !< Rate of lightning-caused ignitions [ha-1 ts-1]

    TYPE(t_jsb_var_real3d)               :: frac_consump_fcl    !< Ratio calculated in eqs. B1, B2, B3 and B4 [0 to 1]
                                                                !< Intermediate variable used in dead fuel consumption

    TYPE(t_jsb_var_real2d)               :: ic_pop_dens         !< Read data from [#events ha-1 day-1]
    TYPE(t_jsb_var_real3d)               :: ic_lightning_events_monthly !< Monthly lightning events data [#events ha-1 day-1]
    TYPE(t_jsb_var_real2d)               :: lightning_events    !< Current monthly timestep data [#events ha-1 day-1]
  CONTAINS
    PROCEDURE :: Init => Init_q_dfire_memory
  END TYPE t_q_dfire_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> initialize memory (variables) for the process: dist_fire
  !>
  SUBROUTINE Init_q_dfire_memory(mem, prefix, suffix, lct_ids, lib_id, model_id)
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_dfire_memory), INTENT(inout), TARGET  :: mem         !> dfire memory
    CHARACTER(len=*),     INTENT(in)            :: prefix          !> process name
    CHARACTER(len=*),     INTENT(in)            :: suffix          !> tile name
    INTEGER,              INTENT(in)            :: lct_ids(:)      !< Primary lct (1) and lcts of descendant tiles
    INTEGER,              INTENT(in)            :: lib_id          !< id of primary lct in lctlib
    INTEGER,              INTENT(in)            :: model_id        !> model ID model\%id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                            !> Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_surface                    !> Vertical grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_dead_fuel_pools            !> Dead fuel pools grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_live_fuel_pools            !> Live fuel pools grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_fuel_classes               !> Fuel classes grid
    TYPE(t_jsb_vgrid), POINTER :: vgrid_fuel_moisture              !> Fuel fuel moisture
    TYPE(t_jsb_vgrid), POINTER :: vgrid_months_in_year             !> Number of months in a year
    INTEGER                    :: table
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_q_dfire_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(model_id)
    table                 =  tables(1)
    hgrid                 => Get_grid(model%grid_id)
    vgrid_surface         => Get_vgrid('surface')
    vgrid_dead_fuel_pools => Get_vgrid(df_sf_vgrid_dead_fuel_pools_name)
    vgrid_live_fuel_pools => Get_vgrid(df_sf_vgrid_live_fuel_pools_name)
    vgrid_fuel_classes    => Get_vgrid(df_sf_vgrid_fuel_classes_name)
    vgrid_fuel_moisture   => Get_vgrid(df_sf_vgrid_fuel_moisture_name)
    vgrid_months_in_year => Get_vgrid(df_sf_vgrid_months_in_year_name)
    ! ----------------------------------------------------------------------------------------------------- !

    CALL mem%Add_var('nesterov_index', mem%nesterov_index, &
        & hgrid, vgrid_surface, &
        & t_cf('nesterov_index', 'Celsius_2', 'Cumulative function of daily tmax and tdew '), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fdi', mem%fdi, &
        & hgrid, vgrid_surface, &
        & t_cf('fdi', '%', 'Fire danger index'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('t_air_max', mem%t_air_max, &
        & hgrid, vgrid_surface, &
        & t_cf('t_air_max', 'K', 'Max temperature of the previous day (daytime_counter)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('t_air_min', mem%t_air_min, &
        & hgrid, vgrid_surface, &
        & t_cf('t_air_min', 'K', 'Min temperature of the previous day (daytime_counter)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 9999.0_wp)

    CALL mem%Add_var('mask_is_fire_intense_enough', mem%mask_is_fire_intense_enough, &
        & hgrid, vgrid_surface, &
        & t_cf('mask_is_fire_instense_enough', '-', 'True when fire is intense enough to exist'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
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
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_mavg', mem%precip_mavg, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_mavg', 'kg m-2 s-1', 'Precipitation moving average using pure method'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_qc_mavg', mem%precip_qc_mavg, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_qc_mavg', 'kg m-2 s-1', 'Precipitation moving average using quincy method'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_daily_mmd', mem%precip_daily_mmd, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daily_mmd', 'mm/day', 'Daily precipitation in Spitfire units'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('precip_daytime', mem%precip_daytime, &
        & hgrid, vgrid_surface, &
        & t_cf('precip_daytime', 'kg m-2 s-1', 'Daily precipitation from previous day, related to sun movement'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('pot_burn_area', mem%pot_burn_area, &
        & hgrid, vgrid_surface, &
        & t_cf('pot_burn_area', 'm2', 'Potential burnable area (forest and grass area)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_1hr', mem%dead_fuel_consumption_1hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_1hr', 'g DM m-2', '1 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_10hr', mem%dead_fuel_consumption_10hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_10hr', 'g DM m-2', '10 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_100hr', mem%dead_fuel_consumption_100hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_100hr', 'g DM m-2', '100 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_fuel_consumption_1000hr', mem%dead_fuel_consumption_1000hr, &
        & hgrid, vgrid_dead_fuel_pools, &
        & t_cf('dead_fuel_consumption_1000hr', 'g DM m-2', '1000 hour dead fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('live_fuel_consumption_1hr', mem%live_fuel_consumption_1hr, &
        & hgrid, vgrid_live_fuel_pools, &
        & t_cf('live_fuel_consumption_1hr', 'g DM m-2', '1 hour live fuel consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('dead_frac_consumption', mem%dead_frac_consumption, &
        & hgrid, vgrid_fuel_classes, &
        & t_cf('dead_frac_consumption', '-', 'Dead fuel fraction of consumption'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ab_moisture_diag', mem%ab_moisture_diag, &
        & hgrid, vgrid_surface, &
        & t_cf('ab_moisture_diag', '-', '(Average of the 1st and 2nd) top soil moisture layers'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('area_burned_fract', mem%area_burned_fract, &
        & hgrid, vgrid_surface, &
        & t_cf('area_burned_fract', '%', 'Fraction of the area burned'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('area_burned', mem%area_burned, &
        & hgrid, vgrid_surface, &
        & t_cf('area_burned', 'ha ts-1', 'Area burned in a gridcell'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('i_surface', mem%i_surface, &
        & hgrid, vgrid_surface, &
        & t_cf('i_surface', 'kW m-1', 'Fire intensity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ros_f', mem%ros_f, &
        & hgrid, vgrid_surface, &
        & t_cf('ros_f', 'm min-1', 'Forward rate of spread of a surface fire'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ros_b', mem%ros_b, &
        & hgrid, vgrid_surface, &
        & t_cf('ros_b', 'm min-1', 'Backward rate of spread of a surface fire'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('mean_fire_area', mem%mean_fire_area, &
        & hgrid, vgrid_surface, &
        & t_cf('mean_fire_area', 'ha', 'Mean fire area'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fire_duration', mem%fire_duration, &
        & hgrid, vgrid_surface, &
        & t_cf('fire_duration', 'min', 'Fire duration'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('e_n_ignitions', mem%e_n_ignitions, &
        & hgrid, vgrid_surface, &
        & t_cf('e_n_ignitions', 'ha-1 ts-1', 'Sum of independent estimates of n_l,ig and n_h,ig ignitions events'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('e_n_humans', mem%e_n_humans, &
        & hgrid, vgrid_surface, &
        & t_cf('e_n_humans', 'ha-1 ts-1', 'Rate of human-caused ignitions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('e_n_lightning', mem%e_n_lightning, &
        & hgrid, vgrid_surface, &
        & t_cf('e_n_lightning', 'ha-1 ts-1', 'Rate of lightning-caused ignitions'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('frac_consump_fcl', mem%frac_consump_fcl, &
        & hgrid, vgrid_fuel_classes, &
        & t_cf('frac_consump_fcl', '-', 'Fraction of consumption (dead)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('distance_travelled', mem%distance_travelled, &
        & hgrid, vgrid_surface, &
        & t_cf('distance_travelled', 'm', 'Total distance travelled by fire'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('l2b_ratio', mem%l2b_ratio, &
        & hgrid, vgrid_surface, &
        & t_cf('l2b_ratio', '-', 'Length to breath ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('reaction_intensity', mem%reaction_intensity, &
        & hgrid, vgrid_surface, &
        & t_cf('reaction_intensity', 'kJ m−2 min−1', 'The energy release rate per unit area of fire front'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('prop_flux_ratio', mem%prop_flux_ratio, &
        & hgrid, vgrid_surface, &
        & t_cf('prop_flux_ratio', '-', 'the proportion of Reaction Intensity that heats adjacent fuel ' // &
            'particles to ignition'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('wind_factor', mem%wind_factor, &
        & hgrid, vgrid_surface, &
        & t_cf('wind_factor', '-', 'Multiplier that accounts for the effect of wind in increasing ' // &
            'the effective value of prop_flux_ratio'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('eff_heating_number', mem%eff_heating_number, &
        & hgrid, vgrid_surface, &
        & t_cf('eff_heating_number', '-', 'the proportion of a fuel particle that is heated to the ' // &
            'proportion of a fuel particle that is heated to'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('preig_heat', mem%preig_heat, &
        & hgrid, vgrid_surface, &
        & t_cf('preig_heat', 'kJ min−1', 'The amount of heat required to ignite a given mass of fuel'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('Optimum reaction velocity', mem%opt_react_velo, &
        & hgrid, vgrid_surface, &
        & t_cf('opt_react_velo', 'min−1', 'Optimum reaction velocity'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('Net fuel load', mem%net_fuel_load, &
        & hgrid, vgrid_surface, &
        & t_cf('net_fuel_load', 'kg m-2', 'Amount of fuel after after substraction of its mineral content'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('Moisture dampening coefficient', mem%moist_damp_coeff, &
        & hgrid, vgrid_surface, &
        & t_cf('moist_damp_coeff', '-', 'Moisture-dampening coefficient'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('Mineral-dampening coefficient', mem%mineral_damp_coeff, &
        & hgrid, vgrid_surface, &
        & t_cf('mineral_damp_coeff', '-', 'Mineral-dampening coefficient'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('fuel_moisture', mem%fuel_moisture, &
        & hgrid, vgrid_fuel_moisture, &
        & t_cf('fuel_moisture', '-', &
            'Relative daily fuel moisture per age class (1hr, 10hr, 100hr, 1000hr) for omega_o and omega_nl'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ic_pop_dens', mem%ic_pop_dens, &
        & hgrid, vgrid_surface, &
        & t_cf('ic_pop_dens', '#humans km-2', &
            'Population density'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('ic_lightning_events_monthly', mem%ic_lightning_events_monthly, &
        & hgrid, vgrid_months_in_year, &
        & t_cf('ic_lightning_events_monthly', '#lightnings h-1 d-1', &
            '# lightning events'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

    CALL mem%Add_var('lightning_events', mem%lightning_events, &
        & hgrid, vgrid_surface, &
        & t_cf('lightning_events', '#lightnings h-1 d-1', &
            '# lightning events'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, & ! DIAGNOSTICS
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

  END SUBROUTINE Init_q_dfire_memory

#endif
END MODULE mo_q_dfire_memory_class
