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
    REAL(wp)              :: sb_test         !< 2 delete/ replace

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
    ! 0.1 InOut
    CLASS(t_sb_config), INTENT(inout) :: config    !< config type for sb
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL       :: active
    REAL(wp)      :: sb_test
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_sb_config'

    NAMELIST /soil_biogeochemistry_ctl/       & 
                          sb_test
     

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting soil_biogeochemistry configuration')

    ! Set defaults
    sb_test                     = 5.1_wp


    ! read the namelist
    ! 
    ! &SOIL_BIOGEOCHEMISTRY_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('SOIL_BIOGEOCHEMISTRY_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, soil_biogeochemistry_ctl)
    CALL close_nml(nml_handler)


    ! pass values as read from file
    config%sb_test                  = sb_test


  END SUBROUTINE Init_sb_config

!#endif
END MODULE mo_sb_config_class
