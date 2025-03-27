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
!> @brief Initialization of the the dfire memory.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
!> @par Revision History
!>
MODULE mo_dfire_init
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_exception,               ONLY: message
  USE mo_jsb_control,             ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dfire_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_init'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Intialize disturbance fire process (after memory has been set up)
  !>
  SUBROUTINE dfire_init(tile)
    USE mo_jsb_class,               ONLY: get_model
    USE mo_jsb_process_class,       ONLY: DFIRE_
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,         ONLY: t_jsb_model

    USE mo_dfire_constants,         ONLY: DFIRE_PRECIP_METHOD_DAYTIME, & !< Method to compute precipitation
                                          DFIRE_PRECIP_METHOD_APROX_MAVG, &
                                          DFIRE_PRECIP_METHOD_MAVG

    dsl4jsb_Use_config(DFIRE_)
    dsl4jsb_Use_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER              :: model
    INTEGER                                      :: iblk, ics, ice, nc, npft
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':dfire_init'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(DFIRE_)
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: dfire_pm_id                                  !< Selected precipitation method
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Real2D_onDomain      :: mask_is_fire_intense_enough 
    dsl4jsb_Real2D_onDomain      :: precip_daytime
    dsl4jsb_Real2D_onDomain      :: precip_qc_mavg
    dsl4jsb_Real2D_onDomain      :: precip_daily_mmd
    dsl4jsb_Real2D_onDomain      :: precip_mavg

    dsl4jsb_Real2D_onDomain      :: nesterov_index
    dsl4jsb_Real2D_onDomain      :: fdi
    dsl4jsb_Real2D_onDomain      :: ros_f
    dsl4jsb_Real2D_onDomain      :: ros_b
    dsl4jsb_Real2D_onDomain      :: mean_fire_area
    dsl4jsb_Real2D_onDomain      :: fire_duration 

    dsl4jsb_Real2D_onDomain      :: area_burned
    dsl4jsb_Real2D_onDomain      :: area_burned_fract
    dsl4jsb_Real3D_onDomain      :: dead_frac_consumption
    dsl4jsb_Real2D_onDomain      :: i_surface

    dsl4jsb_Real3D_onDomain      :: fuel_moisture

    dsl4jsb_Real2D_onDomain      :: e_n_ignitions
    dsl4jsb_Real2D_onDomain      :: e_n_humans
    dsl4jsb_Real2D_onDomain      :: e_n_lightning
    dsl4jsb_Real2D_onDomain      :: ab_moisture_diag

    dsl4jsb_Real2D_onDomain      :: reaction_intensity
    dsl4jsb_Real2D_onDomain      :: prop_flux_ratio
    dsl4jsb_Real2D_onDomain      :: wind_factor
    dsl4jsb_Real2D_onDomain      :: eff_heating_number
    dsl4jsb_Real2D_onDomain      :: preig_heat

    dsl4jsb_Real2D_onDomain      :: opt_react_velo
    dsl4jsb_Real2D_onDomain      :: net_fuel_load
    dsl4jsb_Real2D_onDomain      :: moist_damp_coeff
    dsl4jsb_Real2D_onDomain      :: mineral_damp_coeff

    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_1hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_10hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_100hr
    dsl4jsb_Real3D_onDomain      :: dead_fuel_consumption_1000hr
    dsl4jsb_Real3D_onDomain      :: live_fuel_consumption_1hr

    dsl4jsb_Real2D_onDomain      :: distance_travelled
    dsl4jsb_Real2D_onDomain      :: l2b_ratio
    dsl4jsb_Real2D_onDomain      :: ck
    dsl4jsb_Real2D_onDomain      :: scorch_height

    dsl4jsb_Real3D_onDomain      :: frac_consump_fcl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_active(DFIRE_)) RETURN
    !IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model  => get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(DFIRE_)
    dsl4jsb_Get_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onDomain(DFIRE_, nesterov_index)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, fdi)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, ros_f)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, ros_b)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, mean_fire_area)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, fire_duration)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, e_n_ignitions)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, e_n_humans)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, e_n_lightning)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, ab_moisture_diag)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, reaction_intensity)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, prop_flux_ratio)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, wind_factor)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, eff_heating_number)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, preig_heat)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, opt_react_velo)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, net_fuel_load)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, moist_damp_coeff)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, mineral_damp_coeff)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, mask_is_fire_intense_enough)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_qc_mavg)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_daily_mmd)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_mavg)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, precip_daytime)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, area_burned)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, area_burned_fract)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_frac_consumption)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, i_surface)

    dsl4jsb_Get_var3D_onDomain(DFIRE_, fuel_moisture)

    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_1hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_10hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_100hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, dead_fuel_consumption_1000hr)
    dsl4jsb_Get_var3D_onDomain(DFIRE_, live_fuel_consumption_1hr)

    dsl4jsb_Get_var2D_onDomain(DFIRE_, distance_travelled)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, l2b_ratio)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, ck)
    dsl4jsb_Get_var2D_onDomain(DFIRE_, scorch_height)

    dsl4jsb_Get_var3D_onDomain(DFIRE_, frac_consump_fcl)
    ! ----------------------------------------------------------------------------------------------------- !
    nesterov_index = 0._wp
    fdi = 0._wp
    ros_f = 0._wp
    ros_b = 0._wp
    mean_fire_area = 0._wp
    fire_duration = 0._wp

    e_n_ignitions = 0._wp
    e_n_humans = 0._wp
    e_n_lightning = 0._wp
    ab_moisture_diag = 0._wp

    mask_is_fire_intense_enough = 0._wp !< TRUE=1._wp, FALSE=0._wp
    precip_daytime = 0._wp
    precip_qc_mavg = 0._wp
    precip_daily_mmd = 0._wp
    precip_mavg = 0._wp
    area_burned = 0._wp
    area_burned_fract = 0._wp
    dead_frac_consumption = 0._wp
    i_surface = 0._wp
    dead_fuel_consumption_1hr = 0._wp
    dead_fuel_consumption_10hr = 0._wp
    dead_fuel_consumption_100hr = 0._wp
    dead_fuel_consumption_1000hr = 0._wp
    live_fuel_consumption_1hr = 0._wp

    fuel_moisture = 0._wp

    reaction_intensity = 0._wp
    prop_flux_ratio = 0._wp
    wind_factor = 0._wp
    eff_heating_number = 0._wp
    preig_heat = 0._wp

    distance_travelled = 0._wp
    l2b_ratio = 0._wp
    ck = 0._wp
    scorch_height = 0._wp
    frac_consump_fcl = 0._wp

    !< Chose daily precipitation. Those are calculated in 3 different ways.
    dfire_pm_id = dsl4jsb_Config(DFIRE_)%dfire_precip_method_id
    IF (dfire_pm_id == DFIRE_PRECIP_METHOD_DAYTIME) THEN
      !< daytime
      dsl4jsb_var_ptr(DFIRE_,precip_daily) => dsl4jsb_var_ptr(DFIRE_,precip_daytime)
    ELSE IF (dfire_pm_id == DFIRE_PRECIP_METHOD_MAVG) THEN
      !< moving average
      dsl4jsb_var_ptr(DFIRE_,precip_daily) => dsl4jsb_var_ptr(DFIRE_,precip_mavg)
    ELSE IF (dfire_pm_id == DFIRE_PRECIP_METHOD_APROX_MAVG) THEN
      !< aproximation to moving average (quincy default method)
      dsl4jsb_var_ptr(DFIRE_,precip_daily) => dsl4jsb_var_ptr(DFIRE_,precip_qc_mavg)
    ENDIF

  END SUBROUTINE dfire_init
#endif
END MODULE mo_dfire_init
