!>
!! @brief Contains structures and methods for pplcc config 
!!  -- currently of no use beyond being required by the infrastructure; might change in the future
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2021-02-02)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_pplcc_config_class
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
  PUBLIC :: t_pplcc_config

  TYPE, EXTENDS(t_jsb_config) :: t_pplcc_config

   CONTAINS
     PROCEDURE :: Init => Init_pplcc_config
  END type t_pplcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pplcc_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize pplcc process
  !!
  !! @param[inout]     config     Configuration type of process (t_pplcc_config)
  !!  
  SUBROUTINE Init_pplcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_pplcc_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_pplcc_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pplcc_config'

    CALL message(TRIM(routine), 'Starting pplcc configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_pplcc.nc'
    ic_filename         = 'ic_land_pplcc.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_pplcc_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_pplcc_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    
  END SUBROUTINE Init_pplcc_config

#endif
END MODULE mo_pplcc_config_class
