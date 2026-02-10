!> QUINCY fire-disturbance process config
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### define fire-disturbance config structure, read fire-disturbance namelist and init configuration parameters
!>
MODULE mo_dfire_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish, message_text
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_util,              ONLY: tolower
  USE mo_jsb_config_class,  ONLY: t_jsb_config

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

  ! ======================================================================================================= !
  !>
  !> config TYPE of the dist_fire process, derived from t_jsb_config
  !>
  TYPE, EXTENDS(t_jsb_config) :: t_dfire_config
    LOGICAL           :: flag_dfire                  !< on/off disturbance fire process [boolean]
    CHARACTER(len=20) :: dfire_modelname             !< name of fire model to run
    INTEGER           :: dfire_model_id              !< fire model ID
    ! spitfire
    REAL(wp)          :: density_human_population    !< average human population density in the grid cell [individuals km-2]
    REAL(wp)          :: lightning_frequency         !< avrg number of all lightnings [ha-1 d-1]
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
    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml
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
    REAL(wp)                      :: density_human_population
    REAL(wp)                      :: lightning_frequency
    CHARACTER(len=*), PARAMETER   :: routine = TRIM(modname)//':Init_dfire_config'
    ! ----------------------------------------------------------------------------------------------------- !
    NAMELIST /lnd_dist_fire_nml/        &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & flag_dfire,                 &
      & dfire_modelname,            &
      & density_human_population,   &         ! spitfire
      & lightning_frequency                   ! spitfire
    ! ----------------------------------------------------------------------------------------------------- !
    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat
    ! ----------------------------------------------------------------------------------------------------- !
    CALL message(TRIM(routine), 'Starting dist_fire configuration')
    ! ----------------------------------------------------------------------------------------------------- !
    ! Set defaults
    active                   = .FALSE.
    bc_filename              = 'bc_land_dfire.nc'
    ic_filename              = 'ic_land_dfire.nc'
    flag_dfire               = .FALSE.
    dfire_modelname          = 'spitfire'         ! [spitfire, ...]
    density_human_population = 10.0               ! default for spitfire
    lightning_frequency      = 5.0                ! default for spitfire
    ! ----------------------------------------------------------------------------------------------------- !
    ! read the namelist
    ! &lnd_dist_fire_nml
    nml_handler = open_nml(TRIM(config%namelist_filename))
    ! ----------------------------------------------------------------------------------------------------- !
    nml_unit = position_nml('lnd_dist_fire_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_dist_fire_nml)
    ! ----------------------------------------------------------------------------------------------------- !
    ! pass values as read from file
    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%flag_dfire                = flag_dfire
    config%dfire_modelname           = TRIM(dfire_modelname)
    config%density_human_population  = density_human_population   ! spitfire
    config%lightning_frequency       = lightning_frequency        ! spitfire
    ! ----------------------------------------------------------------------------------------------------- !
    CALL close_nml(nml_handler)
    ! ----------------------------------------------------------------------------------------------------- !
    ! set the fire model ID from 'dfire_modelname'
    dfire_model_id        = Get_dfire_model_id(TRIM(config%dfire_modelname))
    config%dfire_model_id = dfire_model_id
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

#endif
END MODULE mo_dfire_config_class
