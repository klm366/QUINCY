!> Summary: tcq config (currently of no use beyond being required by the infrastructure)
!>
!>### Could contain namelist info for the test conserved quantities, but currently of no use beyond infrastructural requirements
!>
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2021-03-09)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
!>
MODULE mo_tcq_config_class
#ifndef __NO_JSBACH__
  
  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_tcq_config

  TYPE, EXTENDS(t_jsb_config) :: t_tcq_config

   CONTAINS
     PROCEDURE :: Init => Init_tcq_config
  END type t_tcq_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tcq_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize cq process
  !
  SUBROUTINE Init_tcq_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_tcq_config), INTENT(inout) :: config !<Configuration type of process (t_tcq_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    ! -------------------------------------------------------------------------------------------------- !

    NAMELIST /jsb_tcq_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_tcq_config'

    CALL message(TRIM(routine), 'Starting cq configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_tcq.nc'
    ic_filename         = 'ic_land_tcq.nc'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_tcq_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_tcq_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    
  END SUBROUTINE Init_tcq_config

#endif
END MODULE mo_tcq_config_class
