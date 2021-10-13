!>
!! @brief Contains structures and methods for phenology config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-09-09)
!! Julia Nabel (2018-01-10): added l_forestRegrowth flag for optional functionality: forest regrowth
!!                          (derive maxLAI from biomass via allometric relationships)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_pheno_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pheno_config, Init_pheno_config

  TYPE, EXTENDS(t_jsb_config) :: t_pheno_config
    LOGICAL :: l_forestRegrowth ! If TRUE, the maxLAI is derived from biomass via allometric relationships
    CHARACTER(len=15) :: scheme   !< Phenology scheme to use ('logrop', 'climatology')
     !LOGICAL :: first_call_phenology_this_day  ! TRUE marks if a new day has just started and update_phenology is called for the top tile
    CONTAINS
      PROCEDURE :: Init => Init_pheno_config
  END type t_pheno_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_config_class'

CONTAINS

  SUBROUTINE Init_pheno_config(config)

    USE mo_jsb_namelist_iface, ONLY: process_nml, use_namelist_broadcast
    USE mo_jsb_parallel_iface, ONLY: p_bcast, process_mpi_root_id, mpi_comm

    CLASS(t_pheno_config), INTENT(inout) :: config

    LOGICAL :: active, l_forestRegrowth
    CHARACTER(len=15) :: scheme
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_pheno_nml/  &
      & active,               &
      & l_forestRegrowth,     &
      & scheme,               &
      & ic_filename,          &
      & bc_filename

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pheno_config'

    CALL message(TRIM(routine), 'Starting phenology configuration')

    ! Set defaults
    active              = .TRUE.
    l_forestRegrowth    = .FALSE.
    scheme              = 'logrop'
    bc_filename         = 'bc_land_phenology.nc'
    ic_filename         = 'ic_land_phenology.nc'
    !first_call_phenology_this_day = .FALSE.

    ! Read namelist
    CALL process_nml(TRIM(config%namelist_filename), 'jsb_pheno_nml', read_nml)

    IF (use_namelist_broadcast()) THEN
      CALL p_bcast(active,           process_mpi_root_id, mpi_comm)
      CALL p_bcast(bc_filename,      process_mpi_root_id, mpi_comm)
      CALL p_bcast(ic_filename,      process_mpi_root_id, mpi_comm)
      CALL p_bcast(scheme,           process_mpi_root_id, mpi_comm)
      CALL p_bcast(l_forestRegrowth, process_mpi_root_id, mpi_comm)
    END IF

    ! Write resulting values in config memory
    config%active              = active
    config%l_forestRegrowth    = l_forestRegrowth
    config%scheme              = TRIM(scheme)
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename

  CONTAINS

    SUBROUTINE read_nml(nml_unit)
      INTEGER, INTENT(in) :: nml_unit
      READ(nml_unit, NML=jsb_pheno_nml)
    END SUBROUTINE read_nml

  END SUBROUTINE Init_pheno_config

#endif
END MODULE mo_pheno_config_class
