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
!> @brief Contains structures and methods for config of the dist_fire process
!>
!>  -flag_ros_f_prescribed
!>  -flag_mean_fire_area_prescribed: only calculated if ros_f is also calculated. Otherwise an error is raised.
!>      However, if area_burned_fract is prescribed. Then it is ignore anyway.
!>  -flag_area_burned_fract_prescribed: only calculated if ros_f is also calculated. Otherwise an error is raised.
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_dfire_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: tolower
  USE mo_jsb_config_class,  ONLY: t_jsb_config
  USE mo_dfire_constants,   ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr
  USE mo_dfire_constants,   ONLY: df_sf_vgrid_dead_fuel_pools_name, df_sf_vgrid_live_fuel_pools_name, &
                                  df_sf_vgrid_fuel_classes_name, df_sf_vgrid_fuel_moisture_name
  USE mo_dfire_constants,   ONLY: df_sf_omegas_nr

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_dfire_config
  PUBLIC :: DFIRE_MODEL_SF

  !>
  !> Fire models
  !>
  !>   spitfire, ..
  !>
  ENUM, BIND(C)
    ENUMERATOR :: DFIRE_MODEL_SF=1
  END ENUM

  !< Precipitation method names for namelist
  CHARACTER(len=*), PARAMETER :: dfire_cfg_precip_method_daytime = "daytime"
  CHARACTER(len=*), PARAMETER :: dfire_cfg_precip_method_mavg = "mavg" !< Default
  CHARACTER(len=*), PARAMETER :: dfire_cfg_precip_method_aprox_mavg = "mapvg"
  CHARACTER(len=*), PARAMETER :: dfire_cfg_methods = "'daytime', 'mavg', 'mapvg'"

  ! ======================================================================================================= !
  !>
  !> config TYPE of the dist_fire process, derived from t_jsb_config
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_dfire_config
    LOGICAL           :: flag_dfire                  !< on/off disturbance fire process [boolean]
    CHARACTER(len=20) :: dfire_modelname             !< name of fire model to run
    INTEGER           :: dfire_model_id              !< fire model ID

    CHARACTER(len=20) :: dfire_precip_method_name    !< precipitation method name, from namelist
    INTEGER           :: dfire_precip_method_id      !< precipitation method id (1=,2=,3=), internal

    REAL(wp)          :: density_human_population    !< Ignition events
    REAL(wp)          :: lightning_frequency         !< avrg number of all lightnings [ha-1 d-1]
    REAL(wp)          :: a_nd                        !< Propensity of people to produce ignition events (Eq. 3, Thonicke et al. (2010))
                                                     !< [ignitions indivual-1 d-1]      
    REAL(wp)          :: fract_crown_scorch(4)       !< Fraction of live fuel consumed in crown scorch

    !< flags
    REAL(wp)          :: ratio_polymeric_litter_to_dead_fuel !< Ratio to include polymeric litter as dead fuel

    !< Dev Flags
    LOGICAL           :: flag_has_dead_feedback      !< Flag to enable dead fuel feedback from Spitfire into Quincy
    LOGICAL           :: flag_has_live_feedback      !< Flag to enable live fuel feedback from Spitfire into Quincy
    LOGICAL           :: flag_has_prob_mort_feedback !< Flag to enable live probability of mortality from Spitfire into Quincy (q_veg_dynamics.f90)
    REAL(wp)          :: tree_height_factor          !< Tree height factor for spitfire
    REAL(wp)          :: tree_diameter_factor        !< Tree diameter factor for spitfire

  CONTAINS
    PROCEDURE :: Init => Init_dfire_config
  END type t_dfire_config
  
  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_config_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> configuration routine of t_dfire_config
  !>
  SUBROUTINE Init_dfire_config(config)
#ifdef __QUINCY_STANDALONE__
    USE mo_namelist,            ONLY: open_nml, POSITIONED, position_nml, close_nml
#else
    USE mo_jsb_namelist_iface,  ONLY: open_nml, POSITIONED, position_nml, close_nml
#endif
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_io_iface,       ONLY: zaxis_generic
    USE mo_jsb_grid,           ONLY: register_vgrid

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_dfire_config), INTENT(inout) :: config    !< config type for dfire
    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                       :: active
    CHARACTER(len=filename_max)   :: ic_filename, &
                                     bc_filename
    LOGICAL                       :: flag_dfire               !< on/off disturbance fire process [boolean]
    CHARACTER(len=20)             :: dfire_modelname          !< name of fire model to run
    INTEGER                       :: dfire_model_id           !< fire model ID, derived from 'dfire_modelname'
    CHARACTER(len=20)             :: dfire_precip_method_name !< precipitation method name, from namelist
    INTEGER                       :: dfire_precip_method_id   !< precipitation method id (1=,2=,3=), internal
    REAL(wp)                      :: density_human_population !< Ignition events
    REAL(wp)                      :: lightning_frequency      !<
    REAL(wp)                      :: a_nd                     !< 
    REAL(wp)                      :: fract_crown_scorch(4)    !< Fraction of live fuel consumption from crown scorch

    TYPE(t_jsb_vgrid), POINTER    :: vgrid_dead_fuel_pools
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_live_fuel_pools
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_fuel_classes
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_fuel_moisture
    
    REAL(wp)                      :: ratio_polymeric_litter_to_dead_fuel !< Ratio polymeric litter to include as dead fuels
    LOGICAL                       :: flag_has_dead_feedback      !< Flag to enable dead fuel feedback from Spitfire into Quincy
    LOGICAL                       :: flag_has_live_feedback      !< Flag to enable live fuel feedback from Spitfire into Quincy
    LOGICAL                       :: flag_has_prob_mort_feedback !< Flag to enable live probability of mortality from Spitfire into Quincy (q_veg_dynamics.f90)
    REAL(wp)                      :: tree_height_factor          !< Tree height factor for spitfire
    REAL(wp)                      :: tree_diameter_factor        !< Tree diameter factor for spitfire

    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_dfire_config'
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
    NAMELIST /dist_fire_ctl/        &
#else
    NAMELIST /dist_fire_nml/        &
#endif
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & flag_dfire,                 &
      & dfire_modelname,            &
      & dfire_precip_method_name,   &
      & density_human_population,   &         ! ignition events 
      & lightning_frequency,        &         ! 
      & a_nd,                       &         ! 
      &  fract_crown_scorch,        &         ! 
      & ratio_polymeric_litter_to_dead_fuel, &
      & flag_has_dead_feedback,     &
      & flag_has_live_feedback,     &
      & flag_has_prob_mort_feedback, &
      & tree_height_factor,         &
      & tree_diameter_factor


    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting dist_fire configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults, see namelist.slm for more details
    active                   = .FALSE.
    bc_filename              = 'bc_land_dfire.nc'
    ic_filename              = 'ic_land_dfire.nc'
    flag_dfire               = .FALSE.
    dfire_modelname          = 'spitfire'         ! [spitfire, ...]
    dfire_precip_method_name = TRIM(dfire_cfg_precip_method_mavg)  ! 'mavg'
    density_human_population = 55.0               ! Ignition events 
    lightning_frequency      = 0.000012           ! 
    a_nd                     = 0.000048           ! 
    fract_crown_scorch       = [1.0, 0.05, 0., 0.]
    ratio_polymeric_litter_to_dead_fuel = 0._wp
    flag_has_dead_feedback   = .TRUE.
    flag_has_live_feedback   = .TRUE.
    flag_has_prob_mort_feedback = .TRUE.
    tree_height_factor       = 1._wp
    tree_diameter_factor     = 1._wp

    ! ----------------------------------------------------------------------------------------------------- !
    ! read the namelist
    ! &DIST_FIRE_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
#ifdef __QUINCY_STANDALONE__
    nml_unit = position_nml('dist_fire_ctl', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, dist_fire_ctl)
#else
    nml_unit = position_nml('dist_fire_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, dist_fire_nml)
#endif
    ! ----------------------------------------------------------------------------------------------------- !
    ! pass values as read from file
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%flag_dfire                = flag_dfire
    config%dfire_modelname           = TRIM(dfire_modelname)
    config%dfire_precip_method_name  = TRIM(dfire_precip_method_name)
    config%density_human_population  = density_human_population       ! Ignition events 
    config%lightning_frequency       = lightning_frequency            ! 
    config%a_nd                      = a_nd                           ! 
    config%fract_crown_scorch        = fract_crown_scorch             ! 
    config%ratio_polymeric_litter_to_dead_fuel = ratio_polymeric_litter_to_dead_fuel ! spitfire
    config%flag_has_dead_feedback    = flag_has_dead_feedback         !
    config%flag_has_live_feedback    = flag_has_live_feedback         !
    config%flag_has_prob_mort_feedback = flag_has_prob_mort_feedback  !
    config%tree_height_factor        = tree_height_factor             !
    config%tree_diameter_factor      = tree_diameter_factor           !
    ! ----------------------------------------------------------------------------------------------------- !
    CALL close_nml(nml_handler)
    ! ----------------------------------------------------------------------------------------------------- !
    ! set the fire model ID from 'dfire_modelname'
    dfire_model_id        = Get_dfire_model_id(TRIM(config%dfire_modelname))
    config%dfire_model_id = dfire_model_id
    ! ----------------------------------------------------------------------------------------------------- !
    ! set the fire precipitation method ID from 'dfire_precip_method_name'
    dfire_precip_method_id        = Get_Precip_Method_id(TRIM(config%dfire_precip_method_name))
    config%dfire_precip_method_id = dfire_precip_method_id
    ! ----------------------------------------------------------------------------------------------------- !
    vgrid_dead_fuel_pools  => new_vgrid(df_sf_vgrid_dead_fuel_pools_name, zaxis_generic, sf_fuel_dead_fuel_nr, &
                                        longname='Dead fuel pools', &
                                        units='-')

    vgrid_live_fuel_pools  => new_vgrid(df_sf_vgrid_live_fuel_pools_name, zaxis_generic, sf_fuel_live_fuel_nr, &
                                        longname='Live fuel pools', &
                                        units='-')

    vgrid_fuel_classes     => new_vgrid(df_sf_vgrid_fuel_classes_name, zaxis_generic, df_sf_fc_full_nr, &
                                        longname='Fuel classes', &
                                        units='-')
    vgrid_fuel_moisture    => new_vgrid(df_sf_vgrid_fuel_moisture_name, zaxis_generic, df_sf_omegas_nr, &
                                        longname='fuel moisture types (omega_o, omega_nl)', &
                                        & units='0 to 1')

    CALL register_vgrid(vgrid_dead_fuel_pools)
    CALL register_vgrid(vgrid_live_fuel_pools)
    CALL register_vgrid(vgrid_fuel_classes)
    CALL register_vgrid(vgrid_fuel_moisture)
  END SUBROUTINE Init_dfire_config


  ! ====================================================================================================== !
  !>
  !> Returns the id for this fire model name
  !>
  FUNCTION Get_dfire_model_id(name) RESULT(return_value)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER                      :: return_value
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_dfire_model_id'
    ! ----------------------------------------------------------------------------------------------------- !
    SELECT CASE(tolower(TRIM(name)))
    CASE ('spitfire')
      return_value = DFIRE_MODEL_SF
    CASE DEFAULT
      CALL finish(routine,'Invalid fire model name. Use spitfire.')
    END SELECT
  END FUNCTION Get_dfire_model_id


  ! ====================================================================================================== !
  !>
  !> Returns the id for this fire precipitation method name
  !>
  FUNCTION Get_Precip_Method_id(name) RESULT(return_value)
    USE mo_dfire_constants,   ONLY: DFIRE_PRECIP_METHOD_DAYTIME, &
                                    DFIRE_PRECIP_METHOD_MAVG, &
                                    DFIRE_PRECIP_METHOD_APROX_MAVG
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER                      :: return_value
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_Precip_Method_id'
    ! ----------------------------------------------------------------------------------------------------- !
    SELECT CASE(tolower(TRIM(name)))
    CASE (dfire_cfg_precip_method_daytime)
      return_value = DFIRE_PRECIP_METHOD_DAYTIME
    CASE (dfire_cfg_precip_method_mavg)
      return_value = DFIRE_PRECIP_METHOD_MAVG
    CASE (dfire_cfg_precip_method_aprox_mavg)
      return_value = DFIRE_PRECIP_METHOD_APROX_MAVG
    CASE DEFAULT
      CALL finish(routine,"Invalid precipitation method name. "// &
          " Only allowed: "// dfire_cfg_methods)
    END SELECT
  END FUNCTION Get_Precip_Method_id
#endif
END MODULE mo_dfire_config_class
