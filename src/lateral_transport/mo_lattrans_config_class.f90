!>
!! @brief Contains structures and methods for the lateral transport config
!!
!! @author
!! Tobias Stacke, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Tobias Stacke (2022-07-27)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_lattrans_config_class
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
  PUBLIC :: t_lattrans_config

  TYPE, EXTENDS(t_jsb_config) :: t_lattrans_config
     LOGICAL                  :: &
       & l_lattrans                !< enable lateral transport of water between grid cell tiles

   CONTAINS
     PROCEDURE :: Init => Init_lattrans_config
  END type t_lattrans_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lattrans_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize lattrans process
  !!  
  SUBROUTINE Init_lattrans_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_lattrans_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: l_lattrans

    NAMELIST /jsb_lattrans_nml/       &
         active,                      &
         l_lattrans

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lattrans_config'

    CALL message(TRIM(routine), 'Starting lateral transport configuration')

    ! Set defaults
    active              = .TRUE.
    l_lattrans          = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_lattrans_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_lattrans_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%l_lattrans          = l_lattrans
    IF (l_lattrans) THEN 
      CALL message(TRIM(routine), 'Lateral water flow between tiles is enabled')
    ELSE
      CALL message(TRIM(routine), 'Lateral water flow between tiles is disabled')
    END IF
    
  END SUBROUTINE Init_lattrans_config

#endif
END MODULE mo_lattrans_config_class
