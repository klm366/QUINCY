!>
!! @brief Contains structures and methods for carbon config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-29)
!!
MODULE mo_carbon_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_carbon_config, new_carbon_config

  TYPE t_carbon_config
     LOGICAL :: active
     CHARACTER(len=filename_max) :: ic_filename       ! input initial filename
     CHARACTER(len=filename_max) :: bc_filename       ! input boundary-conditions filenname
     LOGICAL :: active_assimilation
!     LOGICAL :: active_Credist ! sehe keinen wissenschaftlichen Sinn darin diesen Schalter unabhängig von den "vegetation cover changes Schaltern" an- und ausschalten zu können.
     LOGICAL :: active_soilresp
     LOGICAL :: active_fireloss
  END type t_carbon_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_carbon_config_class'

CONTAINS

  FUNCTION new_carbon_config(nml_handler) RESULT(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    INTEGER, INTENT(in) :: nml_handler
    TYPE(t_carbon_config), POINTER :: config
    TYPE(t_carbon_config), POINTER :: new

    LOGICAL :: active, active_assimilation, active_soilresp, active_fireloss
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_carbon_nml/                 &
         active,                              &
         ic_filename, bc_filename,            &
         active_assimilation,                 &
         active_soilresp,              & 
         active_fireloss                      

    INTEGER :: nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_carbon_config'

    CALL message(TRIM(routine), 'Starting carbon configuration')

    ALLOCATE(new)

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_carbon.nc'    
    ic_filename         = 'ic_land_carbon.nc'
    active_assimilation = .TRUE.
    active_soilresp     = .TRUE.
    active_fireloss     = .TRUE.

    nml_unit = position_nml('jsb_carbon_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_carbon_nml)

    new%active              = active
    new%ic_filename         = ic_filename
    new%bc_filename         = bc_filename
    new%active_assimilation = active_assimilation
    new%active_soilresp     = active_soilresp
    new%active_fireloss     = active_fireloss

    config => new

  END FUNCTION new_carbon_config

#endif
END MODULE mo_carbon_config_class
