!>
!! @brief Contains structures and methods for lake config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_lake_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lake_config

  TYPE, EXTENDS(t_jsb_config) :: t_lake_config
     REAL(wp)   :: mixed_layer_depth  !< Mixed-layer depth of lakes
   CONTAINS
     PROCEDURE :: Init => Init_lake_config
  END type t_lake_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lake_config_class'

CONTAINS

  SUBROUTINE Init_lake_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_lake_config), INTENT(inout) :: config

    LOGICAL  :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    REAL(wp)                    :: mixed_layer_depth

    NAMELIST /jsb_lake_nml/        &
      active,                      &
      ic_filename, bc_filename,    &
      mixed_layer_depth

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_lake_config'

    CALL message(TRIM(routine), 'Starting lake configuration')

    ! Set defaults
    active            = .TRUE.
    ic_filename       = 'ic_land_srf.nc'
    bc_filename       = 'bc_land_srf.nc'
    mixed_layer_depth = 10._wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_lake_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_lake_nml)

    CALL close_nml(nml_handler)

    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%mixed_layer_depth         = mixed_layer_depth

  END SUBROUTINE Init_lake_config

#endif
END MODULE mo_lake_config_class
