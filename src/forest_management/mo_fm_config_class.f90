!>
!! @brief Contains structures and methods for forest management config
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2018-03-19)
!!
MODULE mo_fm_config_class
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
  PUBLIC :: t_fm_config

  TYPE, EXTENDS(t_jsb_config) :: t_fm_config
    CHARACTER(len=25)         :: fm_scheme
  CONTAINS
     PROCEDURE :: Init => Init_fm_config
  END type t_fm_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fm_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize forest management process
  !!
  !! @param[inout]     config     Configuration type of process (t_fm_config)
  !!
  SUBROUTINE Init_fm_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_fm_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
 
    NAMELIST /jsb_fm_nml/     &
         active,                &
         ic_filename,           &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_fm_config'

    CALL message(TRIM(routine), 'Starting fm configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_fm.nc'
    ic_filename            = 'ic_land_fm.nc'

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_fm_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_fm_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename

    config%fm_scheme       = TRIM("annually_from_file")

  END SUBROUTINE Init_fm_config

#endif
END MODULE mo_fm_config_class
