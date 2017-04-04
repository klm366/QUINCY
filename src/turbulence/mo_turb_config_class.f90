!>
!! @brief Contains structures and methods for turbulence config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2916-04-01)
!!
MODULE mo_turb_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_turb_config

  TYPE, EXTENDS(t_jsb_config) :: t_turb_config
   CONTAINS
     PROCEDURE :: Init => Init_turb_config
  END type t_turb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize turbulence process
  !!
  !! @param[inout]     config     Configuration type of process (t_turb_config)
  !!
  SUBROUTINE Init_turb_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_turb_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_turb_nml/      &
         active,                 &
         ic_filename,            &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_turb_config'

    CALL message(TRIM(routine), 'Starting turb configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_phys.nc'
    ic_filename         = 'ic_land_phys.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_turb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_turb_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename

  END SUBROUTINE Init_turb_config

#endif
END MODULE mo_turb_config_class
