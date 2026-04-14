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
MODULE mo_q_dfire_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: tolower
  USE mo_jsb_config_class,  ONLY: t_jsb_config
  USE mo_q_dfire_constants,   ONLY: sf_fuel_dead_fuel_nr, sf_fuel_live_fuel_nr, df_sf_fc_full_nr
  USE mo_q_dfire_constants,   ONLY: df_sf_vgrid_dead_fuel_pools_name, df_sf_vgrid_live_fuel_pools_name, &
                                  df_sf_vgrid_fuel_classes_name, df_sf_vgrid_fuel_moisture_name, &
                                  df_sf_vgrid_months_in_year_name
  USE mo_q_dfire_constants,   ONLY: df_sf_omegas_nr, & !< fuel moisture (omega_o and omega_n_l)
                                  sf_fuel_classes_nr   !< total number of fuel classes

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_dfire_config
  PUBLIC :: Q_DFIRE_MODEL_SF, Q_DFIRE_ROS_THONICKE, Q_DFIRE_ROS_OBERHAGEMANN

  !>
  !> Fire models
  !>
  !>   spitfire, ..
  !>
  ENUM, BIND(C)
    ENUMERATOR :: Q_DFIRE_MODEL_SF=1
  END ENUM

  !> Rate of spread fire models
  ENUM, BIND(C)
    ENUMERATOR :: Q_DFIRE_ROS_THONICKE = 1 !< Thonicke et al. (2010)
    ENUMERATOR :: Q_DFIRE_ROS_OBERHAGEMANN !< Oberhagemann et al. (2025)
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
  TYPE, EXTENDS(t_jsb_config) :: t_q_dfire_config
    LOGICAL           :: flag_dfire                  !< on/off disturbance fire process [boolean]
    CHARACTER(len=20) :: dfire_modelname             !< name of fire model to run
    INTEGER           :: dfire_model_id              !< fire model ID

    CHARACTER(len=20) :: ros_modelname               !< name of rate of spread implementation
    INTEGER           :: ros_model_id                !< Rate of spread implementation ID

    CHARACTER(len=20) :: dfire_precip_method_name    !< precipitation method name, from namelist
    INTEGER           :: dfire_precip_method_id      !< precipitation method id (1=,2=,3=), internal

    REAL(wp)          :: density_human_population    !< average human population density in the grid cell [individuals km-2]
    REAL(wp)          :: lightning_frequency         !< avrg number of all lightnings [ha-1 d-1]
    REAL(wp)          :: a_nd                        !< Propensity of people to produce ignition events (Eq. 3, Thonicke et al. (2010))
                                                     !< [ignitions indivual-1 d-1]
    REAL(wp)          :: fract_crown_scorch(sf_fuel_classes_nr) !< Fraction of live fuel consumed in crown scorch. 4=fuel classes

    ! input data
    CHARACTER(len=filename_max) :: ic_popdens_filename, ic_lightning_events_filename

    !< Fuel ratio flags (how much fuel allocate to the dead fuel classes)
    REAL(wp)          :: ratio_polymeric_litter_to_dead_fuel !< Ratio to include polymeric litter as dead fuel
    REAL(wp)          :: ratio_dom_assoc_to_dead_fuel        !< Ratio to include dom_assoc as dead fuel
    REAL(wp)          :: ratio_dom_to_dead_fuel              !< Ratio to include dom as dead fuel
    REAL(wp)          :: ratio_soluable_litter_to_dead_fuel  !< Ratio to include soluable as dead fuel
    REAL(wp)          :: ratio_residue_to_dead_fuel          !< Ratio to include residue as dead fuel
    REAL(wp)          :: ratio_residue_assoc_to_dead_fuel    !< Ratio to include residue_assoc as dead fuel

    !< Dev Flags
    LOGICAL           :: flag_has_dead_feedback      !< Flag to enable dead fuel feedback from Spitfire into Quincy
    LOGICAL           :: flag_has_live_feedback      !< Flag to enable live fuel feedback from Spitfire into Quincy
    LOGICAL           :: flag_has_prob_mort_feedback !< Flag to enable live probability of mortality from Spitfire into Quincy (q_veg_dynamics.f90)
    REAL(wp)          :: tree_height_factor          !< Tree height factor for spitfire
    REAL(wp)          :: tree_diameter_factor        !< Tree diameter factor for spitfire
    REAL(wp)          :: dead_fuel_factor            !< Dead fuel factor for spitfire
    REAL(wp)          :: live_fuel_factor            !< Live fuel factor for spitfire
    LOGICAL           :: has_diagnostics             !< Enable diagnostics for spitfire, in case of any inconsistencies it stops the simulation

  CONTAINS
    PROCEDURE :: Init => Init_q_dfire_config
  END type t_q_dfire_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_config_class'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> configuration routine of t_dfire_config
  !>
  SUBROUTINE Init_q_dfire_config(config)
    USE mo_jsb_namelist_iface,  ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_GENERIC
    USE mo_jsb_grid,           ONLY: register_vgrid

    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_dfire_config), INTENT(inout) :: config    !< config type for dfire
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
    CHARACTER(len=20)             :: ros_modelname            !< name of rate of spread model
    INTEGER                       :: ros_model_id             !< rate of spread ID, derived from 'ros_modelname'
    REAL(wp)                      :: density_human_population !< Ignition events
    REAL(wp)                      :: lightning_frequency      !<
    REAL(wp)                      :: a_nd                     !<
    REAL(wp)                      :: fract_crown_scorch(sf_fuel_classes_nr) !< Fraction of live fuel consumption from crown scorch
    CHARACTER(len=filename_max)   :: ic_popdens_filename
    CHARACTER(len=filename_max)   :: ic_lightning_events_filename

    TYPE(t_jsb_vgrid), POINTER    :: vgrid_dead_fuel_pools
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_live_fuel_pools
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_fuel_classes
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_fuel_moisture
    TYPE(t_jsb_vgrid), POINTER    :: vgrid_months_in_year

    REAL(wp)                      :: ratio_polymeric_litter_to_dead_fuel !< Ratio polymeric litter to include as dead fuels
    REAL(wp)                      :: ratio_dom_assoc_to_dead_fuel       !< Ratio to include dom_assoc as dead fuel
    REAL(wp)                      :: ratio_dom_to_dead_fuel             !< Ratio to include dom as dead fuel
    REAL(wp)                      :: ratio_soluable_litter_to_dead_fuel !< Ratio to include soluable as dead fuel
    REAL(wp)                      :: ratio_residue_to_dead_fuel         !< Ratio to include residue as dead fuel
    REAL(wp)                      :: ratio_residue_assoc_to_dead_fuel   !< Ratio to include residue_assoc as dead fuel
    LOGICAL                       :: flag_has_dead_feedback      !< Flag to enable dead fuel feedback from Spitfire into Quincy
    LOGICAL                       :: flag_has_live_feedback      !< Flag to enable live fuel feedback from Spitfire into Quincy
    LOGICAL                       :: flag_has_prob_mort_feedback !< Flag to enable live probability of mortality from Spitfire into Quincy (q_veg_dynamics.f90)
    REAL(wp)                      :: tree_height_factor          !< Tree height factor for spitfire
    REAL(wp)                      :: tree_diameter_factor        !< Tree diameter factor for spitfire
    REAL(wp)                      :: dead_fuel_factor            !< Dead fuel factor for spitfire
    REAL(wp)                      :: live_fuel_factor            !< Live fuel factor for spitfire
    LOGICAL                       :: has_diagnostics             !< Enable diagnostics for spitfire, in case of any inconsistencies it stops the simulation

    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_q_dfire_config'
    ! ----------------------------------------------------------------------------------------------------- !
    NAMELIST /lnd_q_dist_fire_nml/            &
      & active,                               &
      & ic_filename,                          &
      & bc_filename,                          &
      & flag_dfire,                           &
      & dfire_modelname,                      &
      & dfire_precip_method_name,             &
      & ros_modelname,                        &
      & density_human_population,             &
      & lightning_frequency,                  &
      & a_nd,                                 &
      & ic_popdens_filename,                  &
      & ic_lightning_events_filename,         &
      & fract_crown_scorch,                   &
      & ratio_polymeric_litter_to_dead_fuel,  &
      & ratio_dom_assoc_to_dead_fuel,         &
      & ratio_dom_to_dead_fuel,               &
      & ratio_soluable_litter_to_dead_fuel,   &
      & ratio_residue_to_dead_fuel,           &
      & ratio_residue_assoc_to_dead_fuel,     &
      & flag_has_dead_feedback,               &
      & flag_has_live_feedback,               &
      & flag_has_prob_mort_feedback,          &
      & tree_height_factor,                   &
      & tree_diameter_factor,                 &
      & dead_fuel_factor,                     &
      & live_fuel_factor,                     &
      & has_diagnostics

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
    ros_modelname            = 'thonicke'         ! thonicke or oberhagemann
    density_human_population = 55.0               ! Ignition events
    lightning_frequency      = 0.000012           !
    a_nd                     = 0.000048           !
    ic_popdens_filename      = "ic_popdens.nc"
    ic_lightning_events_filename = "ic_lightning_events.nc"
    fract_crown_scorch       = [1.0, 0.05, 0., 0.]
    ratio_polymeric_litter_to_dead_fuel = 0._wp
    ratio_dom_assoc_to_dead_fuel        = 0._wp
    ratio_dom_to_dead_fuel              = 0._wp
    ratio_soluable_litter_to_dead_fuel  = 0._wp
    ratio_residue_to_dead_fuel          = 0._wp
    ratio_residue_assoc_to_dead_fuel    = 0._wp
    flag_has_dead_feedback   = .TRUE.
    flag_has_live_feedback   = .TRUE.
    flag_has_prob_mort_feedback = .TRUE.
    tree_height_factor       = 1._wp
    tree_diameter_factor     = 1._wp
    dead_fuel_factor         = 1.0_wp
    live_fuel_factor         = 1.0_wp
    has_diagnostics          = .FALSE.
    ! ----------------------------------------------------------------------------------------------------- !
    ! read the namelist
    ! &Q_DIST_FIRE_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
    nml_unit = position_nml('lnd_q_dist_fire_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_q_dist_fire_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! pass values as read from file
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%flag_dfire                = flag_dfire
    config%dfire_modelname           = TRIM(dfire_modelname)
    config%dfire_precip_method_name  = TRIM(dfire_precip_method_name)
    config%ros_modelname             = TRIM(ros_modelname)
    config%density_human_population  = density_human_population
    config%lightning_frequency       = lightning_frequency
    config%a_nd                      = a_nd
    config%ic_popdens_filename       = TRIM(ic_popdens_filename)
    config%ic_lightning_events_filename = TRIM(ic_lightning_events_filename)
    config%fract_crown_scorch        = fract_crown_scorch
    config%ratio_polymeric_litter_to_dead_fuel = ratio_polymeric_litter_to_dead_fuel
    config%ratio_dom_assoc_to_dead_fuel        = ratio_dom_assoc_to_dead_fuel
    config%ratio_dom_to_dead_fuel              = ratio_dom_to_dead_fuel
    config%ratio_soluable_litter_to_dead_fuel  = ratio_soluable_litter_to_dead_fuel
    config%ratio_residue_to_dead_fuel          = ratio_residue_to_dead_fuel
    config%ratio_residue_assoc_to_dead_fuel    = ratio_residue_assoc_to_dead_fuel
    config%flag_has_dead_feedback    = flag_has_dead_feedback
    config%flag_has_live_feedback    = flag_has_live_feedback
    config%flag_has_prob_mort_feedback = flag_has_prob_mort_feedback
    config%tree_height_factor        = tree_height_factor
    config%tree_diameter_factor      = tree_diameter_factor
    config%dead_fuel_factor          = dead_fuel_factor
    config%live_fuel_factor          = live_fuel_factor
    config%has_diagnostics           = has_diagnostics

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
    ! set the fire rate of spread method ID from 'ros_modelname'
    ros_model_id = Get_dfire_ros_id(TRIM(config%ros_modelname))
    config%ros_model_id = ros_model_id
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
    vgrid_months_in_year   => new_vgrid(df_sf_vgrid_months_in_year_name, zaxis_generic, 12, & ! 12= months in a year
        longname='Months in a year', &
        & units='#')

    CALL register_vgrid(vgrid_dead_fuel_pools)
    CALL register_vgrid(vgrid_live_fuel_pools)
    CALL register_vgrid(vgrid_fuel_classes)
    CALL register_vgrid(vgrid_fuel_moisture)
    CALL register_vgrid(vgrid_months_in_year)
  END SUBROUTINE Init_q_dfire_config


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
      return_value = Q_DFIRE_MODEL_SF
    CASE DEFAULT
      CALL finish(routine,'Invalid fire model name. Use spitfire.')
    END SELECT
  END FUNCTION Get_dfire_model_id


  ! ====================================================================================================== !
  !>
  !> Returns the id for this Rate of spread model name
  !>
  FUNCTION Get_dfire_ros_id(name) RESULT(return_value)
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER                      :: return_value
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_dfire_ros_id'
    ! ----------------------------------------------------------------------------------------------------- !
    SELECT CASE(tolower(TRIM(name)))
    CASE ('thonicke')
      return_value = Q_DFIRE_ROS_THONICKE
    CASE ('oberhagemann')
      return_value = Q_DFIRE_ROS_OBERHAGEMANN
    CASE DEFAULT
      CALL finish(routine,'Invalid fire Rate of Spread')
    END SELECT
  END FUNCTION Get_dfire_ros_id


  ! ====================================================================================================== !
  !>
  !> Returns the id for this fire precipitation method name
  !>
  FUNCTION Get_Precip_Method_id(name) RESULT(return_value)
    USE mo_q_dfire_constants,   ONLY: Q_DFIRE_PRECIP_METHOD_DAYTIME, &
                                      Q_DFIRE_PRECIP_METHOD_MAVG, &
                                      Q_DFIRE_PRECIP_METHOD_APROX_MAVG
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER                      :: return_value
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_Precip_Method_id'
    ! ----------------------------------------------------------------------------------------------------- !
    SELECT CASE(tolower(TRIM(name)))
    CASE (dfire_cfg_precip_method_daytime)
      return_value = Q_DFIRE_PRECIP_METHOD_DAYTIME
    CASE (dfire_cfg_precip_method_mavg)
      return_value = Q_DFIRE_PRECIP_METHOD_MAVG
    CASE (dfire_cfg_precip_method_aprox_mavg)
      return_value = Q_DFIRE_PRECIP_METHOD_APROX_MAVG
    CASE DEFAULT
      CALL finish(routine,"Invalid precipitation method name. "// &
          " Only allowed: "// dfire_cfg_methods)
    END SELECT
  END FUNCTION Get_Precip_Method_id
#endif
END MODULE mo_q_dfire_config_class
