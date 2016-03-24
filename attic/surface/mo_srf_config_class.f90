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
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_srf_config

  TYPE, EXTENDS(t_jsb_config) :: t_srf_config
     CHARACTER(len=filename_max) :: bc_sso_filename
     CHARACTER(len=10)           :: scheme_albedo !< Scheme to use for computation of albedo (echam5 only for now)
     REAL(wp)                    :: wsn_critical  !< Critical snow depth for correction of new surface temperature [m water equivalent]
   CONTAINS
     PROCEDURE :: Init => Init_srf_config
  END type t_srf_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_srf_config_class'

CONTAINS

  SUBROUTINE Init_srf_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_srf_config), INTENT(inout) :: config

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

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_srf_config'

    CALL message(TRIM(routine), 'Starting srf configuration')

    ! Set defaults
    active          = .TRUE.
    ic_filename     = 'ic_land_srf.nc'
    bc_filename     = 'bc_land_srf.nc'
    bc_sso_filename = 'bc_land_sso.nc'
    scheme_albedo   = 'echam5'
    wsn_critical    = 5.85036E-3_wp

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_srf_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_srf_nml)

    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%bc_sso_filename           = bc_sso_filename
    config%scheme_albedo             = scheme_albedo
    config%wsn_critical              = wsn_critical

    CALL close_nml(nml_handler)

  END SUBROUTINE Init_srf_config

#endif
END MODULE mo_srf_config_class
