!>
!! @brief Contains structures and methods for carbon config
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

MODULE mo_carbon_config_class
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
  PUBLIC :: t_carbon_config

  TYPE, EXTENDS(t_jsb_config) :: t_carbon_config
    LOGICAL  :: read_cpools ! Read Cpools from external file when experiment is newly started (cpools are
                            ! initialized with values from file) or restarted (cpools from restart file are
                            ! overwritten by values from file).
                            ! ATTENTION: The run script has to take care that read_cpools is set to .FALSE.
                            !            in the experiment jobs following the first one! Or re-set it manually!
    LOGICAL  :: diag_humus_fluxes ! Activate diagnostic carbon fluxes to and from the yasso humus pools
                                  ! required for analytical humus pool equilibration (compare jsb3.2 svn rev 9934)
    REAL(wp) :: fire_fract_wood_2_atmos

   CONTAINS
     PROCEDURE :: Init => Init_carbon_config
  END type t_carbon_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize carbon process
  !!
  !! @param[inout]     config     Configuration type of process (t_carbon_config)
  !!
  SUBROUTINE Init_carbon_config(config)

    USE mo_jsb_namelist_iface, ONLY: process_nml, use_namelist_broadcast
    USE mo_jsb_parallel_iface, ONLY: p_bcast, process_mpi_root_id, mpi_comm
    USE mo_carbon_constants,   ONLY: FireFracWood2Atmos

    CLASS(t_carbon_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    LOGICAL                     :: read_cpools, diag_humus_fluxes
    REAL(wp)                    :: fire_fract_wood_2_atmos
    ! LOGICAL           :: lemissions

    NAMELIST /jsb_carbon_nml/      &
         active,                   &
         ic_filename,              &
         bc_filename,              &
         read_cpools,              &
         diag_humus_fluxes,        &
         fire_fract_wood_2_atmos

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_carbon_config'

    CALL message(TRIM(routine), 'Starting carbon configuration')

    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_carbon.nc'
    ic_filename            = 'ic_land_carbon.nc'
    read_cpools            = .FALSE.
    diag_humus_fluxes      = .FALSE.
    fire_fract_wood_2_atmos = FireFracWood2Atmos

    ! Read namelist
    CALL process_nml(TRIM(config%namelist_filename), 'jsb_carbon_nml', read_nml)

    IF (use_namelist_broadcast()) THEN
      CALL p_bcast(active,                  process_mpi_root_id, mpi_comm)
      CALL p_bcast(bc_filename,             process_mpi_root_id, mpi_comm)
      CALL p_bcast(ic_filename,             process_mpi_root_id, mpi_comm)
      CALL p_bcast(read_cpools,             process_mpi_root_id, mpi_comm)
      CALL p_bcast(diag_humus_fluxes,       process_mpi_root_id, mpi_comm)
      CALL p_bcast(fire_fract_wood_2_atmos, process_mpi_root_id, mpi_comm)
    ENDIF

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%read_cpools            = read_cpools
    config%diag_humus_fluxes      = diag_humus_fluxes
    config%fire_fract_wood_2_atmos = fire_fract_wood_2_atmos

  CONTAINS

    SUBROUTINE read_nml(nml_unit)
      INTEGER, INTENT(in) :: nml_unit
      READ(nml_unit, NML=jsb_carbon_nml)
    END SUBROUTINE read_nml

  END SUBROUTINE Init_carbon_config

#endif
END MODULE mo_carbon_config_class
