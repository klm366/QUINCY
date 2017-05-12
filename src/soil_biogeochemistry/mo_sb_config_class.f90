!>
!! @brief Contains structures and methods for soil_biogeochemistry config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! 2nd version by Jan Engel, MPI BGC (20170104) -> integration into SLM
!!
MODULE mo_sb_config_class
!#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sb_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the soil_biogeochemistry process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_sb_config
    CHARACTER(15)         :: sb_model_scheme
    LOGICAL               :: flag_sb_prescribe_nh4
    LOGICAL               :: flag_sb_prescribe_no3
    LOGICAL               :: flag_sb_prescribe_po4
    LOGICAL               :: flag_n_deposition

    CONTAINS
      PROCEDURE :: Init => Init_sb_config
  END Type t_sb_config

  TYPE(t_sb_config), SAVE :: sb_config


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_sb_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_sb_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    USE mo_base_config_class,  ONLY: base_config
    ! 0.1 InOut
    CLASS(t_sb_config), INTENT(inout) :: config    !< config type for sb
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL       :: active
    CHARACTER(15) :: sb_model_scheme
    LOGICAL               :: flag_sb_prescribe_nh4
    LOGICAL               :: flag_sb_prescribe_no3
    LOGICAL               :: flag_sb_prescribe_po4
    LOGICAL               :: flag_n_deposition
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_config'

    NAMELIST /soil_biogeochemistry_ctl/       & 
                      sb_model_scheme,        &
                      flag_sb_prescribe_nh4,  &
                      flag_sb_prescribe_no3,  &
                      flag_sb_prescribe_po4,  &
                      flag_n_deposition
     

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting soil_biogeochemistry configuration')

    ! Set defaults
    sb_model_scheme             = "simple"  ! simple jsm 
    flag_sb_prescribe_nh4       = .FALSE.
    flag_sb_prescribe_no3       = .FALSE.
    flag_sb_prescribe_po4       = .FALSE.
    flag_n_deposition           = .FALSE.

    ! read the namelist
    ! 
    ! &SOIL_BIOGEOCHEMISTRY_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('SOIL_BIOGEOCHEMISTRY_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, soil_biogeochemistry_ctl)
    CALL close_nml(nml_handler)

    SELECT CASE (TRIM(sb_model_scheme))
      CASE ("simple")
      CASE ("jsm")
         IF(TRIM(base_config%slm_model_name) /= "soil") THEN
            WRITE(message_text,'(3a)') 'setting sb_model_scheme to ',TRIM(sb_model_scheme), & 
                                       'is currently only possible configuration for the standalone soil model'
            CALL finish("mo_sb_config_class:init_sb_config",message_text) 
         ENDIF 
      CASE DEFAULT
        WRITE(message_text,'(2a)') 'invalid sb_model_scheme',sb_model_scheme
        CALL finish("mo_sb_config_class:init_sb_config",message_text) 
    ENDSELECT

    ! pass values as read from file
    config%sb_model_scheme                 = sb_model_scheme
    config%flag_sb_prescribe_nh4           = flag_sb_prescribe_nh4
    config%flag_sb_prescribe_no3           = flag_sb_prescribe_no3
    config%flag_sb_prescribe_po4           = flag_sb_prescribe_po4
    config%flag_n_deposition               = flag_n_deposition


  END SUBROUTINE Init_sb_config

!#endif
END MODULE mo_sb_config_class
