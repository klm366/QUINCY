!>
!! @brief Contains structures and methods for config of the assimilation process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04);
!! integration into SLM by Jan Engel, MPI BGC (20160722)
!!
MODULE mo_assimi_config_class
!#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

!  PUBLIC :: t_assimi_config
  PUBLIC :: assimi_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the assimilation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_assimi_config
    LOGICAL          :: flag_optimal_Nfraction     !< optimise leaf internal N allocation
    LOGICAL          :: flag_t_jmax_acclimation    !< acclimation of optimum temperature for Jmax 
    CONTAINS
      PROCEDURE :: Init => Init_assimi_config
  END Type t_assimi_config

  TYPE(t_assimi_config), SAVE :: assimi_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_assimi_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_assimi_config(config)   ! 'config' == t_assimi_config, see ABSTRACT INTERFACE -> SUBROUTINE Init_config(config)
!    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
!    using the old JSBACH3 functions instead of wrappers implemented in JSBACH4 -> mo_jsb_namelist_iface
    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file
    ! 0.1 InOut
    CLASS(t_assimi_config), INTENT(inout) :: config    !< config type for assimi
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    !LOGICAL      :: active
    LOGICAL      :: flag_optimal_Nfraction, flag_t_jmax_acclimation

    NAMELIST /assimilation_ctl/     & 
      flag_optimal_Nfraction,       &
      flag_t_jmax_acclimation

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_assimi_config'

    CALL message(TRIM(routine), 'Starting assimi configuration')

    ! Set defaults
    flag_optimal_Nfraction         = .FALSE.  ! .TRUE. .FALSE.
    flag_t_jmax_acclimation        = .FALSE.  ! .TRUE. .FALSE.

    ! read the namelist
    ! 
    ! &ASSIMILATION_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('ASSIMILATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, assimilation_ctl)
    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%flag_optimal_Nfraction        = flag_optimal_Nfraction
    config%flag_t_jmax_acclimation       = flag_t_jmax_acclimation

  END SUBROUTINE Init_assimi_config

!#endif
END MODULE mo_assimi_config_class
