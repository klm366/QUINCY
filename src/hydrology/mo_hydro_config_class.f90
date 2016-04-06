!>
!! @brief Contains structures and methods for hydrology config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-10-30)
!!
MODULE mo_hydro_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hydro_config

  TYPE, EXTENDS(t_jsb_config) :: t_hydro_config
  CONTAINS
    PROCEDURE :: Init => Init_hydro_config
  END type t_hydro_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_config_class'

CONTAINS

  SUBROUTINE Init_hydro_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_hydro_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_hydro_nml/         &
         active,                         &
         ic_filename,                    &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_hydro_config'

    CALL message(TRIM(routine), 'Starting hydro configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_hydro.nc'  !R: gibts im Moment noch nicht
    ic_filename         = 'ic_land_hydro.nc'  !R: gibts im Moment noch nicht

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_hydro_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_hydro_nml)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename

    CALL close_nml(nml_handler)

  END SUBROUTINE Init_hydro_config

#endif
END MODULE mo_hydro_config_class
