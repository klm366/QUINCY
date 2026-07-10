!>
!! @brief Contains structures and methods for flake config
!!
!! @author
!!  Zoé Rehder, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Zoé Rehder (2023-06-01)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_flake_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USe mo_util_string,       ONLY: real2string
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_flake_config

  TYPE, EXTENDS(t_jsb_config) :: t_flake_config
    REAL(wp) :: depth ! Lake depth [m]
  CONTAINS
      PROCEDURE :: Init => Init_flake_config
  END type t_flake_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_flake_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize flake process
  !!
  !! @param[inout]     config     Configuration type of process (t_flake_config)
  !!  
  SUBROUTINE Init_flake_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_flake_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    REAL(wp)                    :: depth

    NAMELIST /jsb_flake_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename, &
         depth               

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_flake_config'

    CALL message(TRIM(routine), 'Starting flake configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_flake.nc'
    ic_filename         = 'ic_land_flake.nc'
    depth               = 10.0_wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_flake_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_flake_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%depth               = depth

    CALL message(routine, 'Using lake depth = '//TRIM(real2string(depth))//' m')
    
  END SUBROUTINE Init_flake_config

#endif
END MODULE mo_flake_config_class
