!>
!! @brief Contains structures and methods for jedi config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (20.05.2017)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jedi_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_jedi_config

  TYPE, EXTENDS(t_jsb_config) :: t_jedi_config
    LOGICAL  :: l_spinup_jedi
  CONTAINS
     PROCEDURE :: Init => Init_jedi_config
  END type t_jedi_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jedi_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize carbon process
  !!
  !! @param[inout]     config     Configuration type of process (t_jedi_config)
  !!
  SUBROUTINE Init_jedi_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_jedi_config), INTENT(inout) :: config

    LOGICAL                     :: active
    LOGICAL                     :: single_tile_output
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: l_spinup_jedi

    NAMELIST /jsb_jedi_nml/      &
         active,                   &
         single_tile_output,       &
         ic_filename,              &
         bc_filename,              &
         l_spinup_jedi

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_jedi_config'

    CALL message(TRIM(routine), 'Starting carbon configuration')

    ! Set defaults
    active                 = .FALSE.
    single_tile_output     = .TRUE.
    bc_filename            = 'bc_land_carbon.nc'
    ic_filename            = 'ic_land_carbon.nc'
    l_spinup_jedi          = .FALSE.

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_jedi_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_jedi_nml)

    CALL close_nml(nml_handler)

    ! Write resulting values in config memory
    config%active                 = active
    config%single_tile_output     = single_tile_output
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%l_spinup_jedi          = l_spinup_jedi

  END SUBROUTINE Init_jedi_config

#endif
END MODULE mo_jedi_config_class
