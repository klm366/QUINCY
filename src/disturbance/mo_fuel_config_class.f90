!>
!! @brief Contains structures and methods for fuel config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (06.04.2017)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_fuel_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_fuel_config

  TYPE, EXTENDS(t_jsb_config) :: t_fuel_config
    INTEGER :: fuel_algorithm  ! Select fuel algorithm: 1 : jsbach (default), 2 : Arora and Boer (2005),
                               !                        3 : Thonicke et al.(2001), 4 : Read/GFED
  CONTAINS
     PROCEDURE :: Init => Init_fuel_config
  END type t_fuel_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize fuel process
  !!
  !! @param[inout]     config     Configuration type of process (t_fuel_config)
  !!
  SUBROUTINE Init_fuel_config(config)

    USE mo_jsb_namelist_iface, ONLY: process_nml, use_namelist_broadcast
    USE mo_jsb_parallel_iface, ONLY: p_bcast, process_mpi_root_id, mpi_comm

    CLASS(t_fuel_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    INTEGER                     :: fuel_algorithm

    NAMELIST /jsb_fuel_nml/     &
        & active,               &
        & ic_filename,          &
        & bc_filename,          &
        & fuel_algorithm

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_fuel_config'

    CALL message(TRIM(routine), 'Starting fuel configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_fuel.nc'
    ic_filename            = 'ic_land_fuel.nc'
    fuel_algorithm         =  1

    ! Read namelist
    CALL process_nml(TRIM(config%namelist_filename), 'jsb_fuel_nml', read_nml)

    IF (use_namelist_broadcast()) THEN
      CALL p_bcast(active,         process_mpi_root_id, mpi_comm)
      CALL p_bcast(bc_filename,    process_mpi_root_id, mpi_comm)
      CALL p_bcast(ic_filename,    process_mpi_root_id, mpi_comm)
      CALL p_bcast(fuel_algorithm, process_mpi_root_id, mpi_comm)
    END IF

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%fuel_algorithm         = fuel_algorithm

  CONTAINS

    SUBROUTINE read_nml(nml_unit)
      INTEGER, INTENT(in) :: nml_unit
      READ(nml_unit, NML=jsb_fuel_nml)
    END SUBROUTINE read_nml

  END SUBROUTINE Init_fuel_config

#endif
END MODULE mo_fuel_config_class
