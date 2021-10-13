!> Summary: pplcc (pre- and post lcc) config (currently of no use beyond being required by the infrastructure)
!>
!>### Could contain namelist info for pplcc, but currently of no use beyond infrastructural requirements
!>
!>
!> Copyright (c) 2021 Julia Nabel
!>
!>#### Revision History
!> First version and ongoing changes   by Julia Nabel (2021-02-02)
!>
!>#### License
!>
!> This code is subject to the MPI-M-ICON - License - Agreement in it's most recent form.
!> Please see [this link](https://code.mpimet.mpg.de/attachments/download/20888/MPI-M-ICONLizenzvertragV2.6.pdf) for details.
!> Where software is supplied by third parties, it is indicated in the headers of the routines.
!>
MODULE mo_pplcc_config_class
#ifndef __NO_JSBACH__
  
  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_pplcc_config

  TYPE, EXTENDS(t_jsb_config) :: t_pplcc_config

   CONTAINS
     PROCEDURE :: Init => Init_pplcc_config
  END type t_pplcc_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pplcc_config_class'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize pplcc process
  !
  SUBROUTINE Init_pplcc_config(config)

    USE mo_jsb_namelist_iface, ONLY: process_nml, use_namelist_broadcast
    USE mo_jsb_parallel_iface, ONLY: p_bcast, process_mpi_root_id, mpi_comm

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_pplcc_config), INTENT(inout) :: config !<Configuration type of process (t_pplcc_config)
    ! -------------------------------------------------------------------------------------------------- !
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    ! -------------------------------------------------------------------------------------------------- !
    
    NAMELIST /jsb_pplcc_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pplcc_config'

    CALL message(TRIM(routine), 'Starting pplcc configuration')

    ! Set defaults
    active              = .FALSE.
    bc_filename         = 'bc_land_pplcc.nc'
    ic_filename         = 'ic_land_pplcc.nc'

    CALL process_nml(TRIM(config%namelist_filename), 'jsb_pplcc_nml', read_nml)

    IF (use_namelist_broadcast()) THEN
      CALL p_bcast(active,            process_mpi_root_id, mpi_comm)
      CALL p_bcast(bc_filename,       process_mpi_root_id, mpi_comm)
      CALL p_bcast(ic_filename,       process_mpi_root_id, mpi_comm)
    END IF

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    
  CONTAINS

    SUBROUTINE read_nml(nml_unit)
      INTEGER, INTENT(in) :: nml_unit
      READ(nml_unit, NML=jsb_pplcc_nml)
    END SUBROUTINE read_nml

  END SUBROUTINE Init_pplcc_config

#endif
END MODULE mo_pplcc_config_class
