!>
!! @brief Contains structures and methods for seb config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_seb_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_seb_config

  TYPE, EXTENDS(t_jsb_config) :: t_seb_config
     CHARACTER(len=10)           :: scheme_albedo !< Scheme to use for computation of albedo (echam5 only for now)
     REAL(wp)                    :: lake_mixed_layer_depth !< Depth of lake mixed layer                             [m]
     REAL(wp)                    :: lake_min_ice_depth     !< Minimum ice thickness to start ice formation on lakes [m]
     LOGICAL                     :: l_ice_on_lakes         !< Whether to use ice on lakes
   CONTAINS
     PROCEDURE :: Init => Init_seb_config
  END type t_seb_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_config_class'

CONTAINS

  SUBROUTINE Init_seb_config(config)

    USE mo_jsb_namelist_iface, ONLY: process_nml, use_namelist_broadcast
    USE mo_jsb_parallel_iface, ONLY: p_bcast, process_mpi_root_id, mpi_comm

    CLASS(t_seb_config), INTENT(inout) :: config

    LOGICAL  :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(len=10)           :: scheme_albedo
    REAL(wp)                    :: lake_mixed_layer_depth, lake_min_ice_depth
    LOGICAL                     :: l_ice_on_lakes

    NAMELIST /jsb_seb_nml/         &
      active,                      &
      ic_filename, bc_filename,    &
      scheme_albedo,               &
      lake_mixed_layer_depth,      &
      lake_min_ice_depth,          &
      l_ice_on_lakes

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_seb_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting seb configuration')

    ! Set defaults
    active          = .TRUE.
    ic_filename     = 'ic_land_srf.nc'
    bc_filename     = 'bc_land_srf.nc'
    scheme_albedo   = 'echam5'
    lake_mixed_layer_depth = 10._wp
    lake_min_ice_depth     = 0.05_wp
    l_ice_on_lakes  = .TRUE.

    CALL process_nml(TRIM(config%namelist_filename), 'jsb_seb_nml', read_nml)

    IF (use_namelist_broadcast()) THEN
      CALL p_bcast(active,                 process_mpi_root_id, mpi_comm)
      CALL p_bcast(bc_filename,            process_mpi_root_id, mpi_comm)
      CALL p_bcast(ic_filename,            process_mpi_root_id, mpi_comm)
      CALL p_bcast(scheme_albedo,          process_mpi_root_id, mpi_comm)
      CALL p_bcast(lake_mixed_layer_depth, process_mpi_root_id, mpi_comm)
      CALL p_bcast(lake_min_ice_depth,     process_mpi_root_id, mpi_comm)
      CALL p_bcast(l_ice_on_lakes,         process_mpi_root_id, mpi_comm)
    END IF

    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%scheme_albedo             = scheme_albedo
    config%lake_mixed_layer_depth    = lake_mixed_layer_depth
    config%lake_min_ice_depth        = lake_min_ice_depth
    config%l_ice_on_lakes            = l_ice_on_lakes

  CONTAINS

    SUBROUTINE read_nml(nml_unit)
      INTEGER, INTENT(in) :: nml_unit
      READ(nml_unit, NML=jsb_seb_nml)
    END SUBROUTINE read_nml

  END SUBROUTINE Init_seb_config

#endif
END MODULE mo_seb_config_class
