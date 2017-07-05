!>
!! @brief Contains structures and methods for config of the phenology process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Jan Engel, MPI BGC (20161214)
!!
MODULE mo_pheno_config_class
!#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp, eps4
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_pheno_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the phenology process, derived from t_jsb_config
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_pheno_config
    REAL(wp)        :: test_pheno

  CONTAINS
    PROCEDURE :: Init => Init_pheno_config
  END Type t_pheno_config

  
  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_pheno_config
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_pheno_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    IMPLICIT NONE
    ! 0.1 InOut
    CLASS(t_pheno_config), INTENT(inout) :: config    !< config type for pheno
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    !LOGICAL         :: active
    REAL(wp)        :: test_pheno
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_pheno_config'


    NAMELIST /phenology_ctl/       & 
                test_pheno      !    , &


    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting pheno configuration')

    ! Set defaults
    test_pheno                = 0.1_wp

    ! read the namelist
    ! 
    ! &PHENOLOGY_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('PHENOLOGY_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, phenology_ctl)
    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%test_pheno                   = test_pheno


  END SUBROUTINE Init_pheno_config

!#endif
END MODULE mo_pheno_config_class
