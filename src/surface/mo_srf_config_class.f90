!>
!! @brief Contains structures and methods for srf config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_srf_config, new_srf_config

  TYPE t_srf_config
     LOGICAL :: active
     CHARACTER(len=filename_max) :: ic_filename
     CHARACTER(len=filename_max) :: bc_filename
     CHARACTER(len=filename_max) :: bc_sso_filename
     CHARACTER(len=10)           :: scheme_albedo !< Scheme to use for computation of albedo (echam5 only for now)
     REAL(wp)                    :: wsn_critical  !< Critical snow depth for correction of new surface temperature [m water equivalent]
  END type t_srf_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_config_class'

CONTAINS

  FUNCTION new_srf_config(nml_handler) RESULT(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    INTEGER,          INTENT(in) :: nml_handler
    TYPE(t_srf_config), POINTER :: config

    TYPE(t_srf_config), POINTER :: new

    LOGICAL  :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename, bc_sso_filename
    CHARACTER(len=10)           :: scheme_albedo
    REAL(wp)                    :: wsn_critical

    NAMELIST /jsb_srf_nml/         &
      active,                      &
      ic_filename, bc_filename,    &
      bc_sso_filename,             &
      scheme_albedo,               &
      wsn_critical

    INTEGER :: nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_srf_config'

    CALL message(TRIM(routine), 'Starting srf configuration')

    ALLOCATE(new)

    ! Set defaults
    active          = .TRUE.
    ic_filename     = 'ic_land_srf.nc'
    bc_filename     = 'bc_land_srf.nc'
    bc_sso_filename = 'bc_land_sso.nc'
    scheme_albedo   = 'echam5'
    wsn_critical    = 5.85036E-3_wp

    nml_unit = position_nml('jsb_srf_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_srf_nml)

    new%active                    = active
    new%ic_filename               = ic_filename
    new%bc_filename               = bc_filename
    new%bc_sso_filename           = bc_sso_filename
    new%scheme_albedo             = scheme_albedo
    new%wsn_critical              = wsn_critical

    config => new

  END FUNCTION new_srf_config

#endif
END MODULE mo_srf_config_class
