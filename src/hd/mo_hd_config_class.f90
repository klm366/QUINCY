!>
!! @brief Contains structures and methods for hd config
!!
!! @author
!!  Veronika Gayler, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Veroniks Gayler (2013-04-20)
!!
MODULE mo_hd_config_class
!#ifndef __NO_JSBACH__
#if !defined(__NO_JSBACH__) && !defined(__NO_JSBACH_HD__)

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hd_config, new_hd_config

  TYPE t_hd_config
     LOGICAL                     :: active                !< hd module is active
     CHARACTER(len=filename_max) :: ic_filename           !< file with initial conditions
     CHARACTER(len=filename_max) :: bc_filename           !< file with boundary conditions       
     LOGICAL                     :: debug_hd              !< additional printing for debugging
     LOGICAL                     :: diag_water_budget     !< prints to diagnose global water budget
     CHARACTER(len=20)           :: routing_scheme        !< river_routing_scheme: full / equally_to_coast / equally_to_ocean
  END type t_hd_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_config_class'

CONTAINS

  FUNCTION new_hd_config(nml_handler) RESULT(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    INTEGER,          INTENT(in) :: nml_handler
    TYPE(t_hd_config), POINTER :: config

    TYPE(t_hd_config), POINTER :: new

    LOGICAL                       :: active                !< hd module is active
    CHARACTER(len=filename_max)   :: ic_filename           !< file with initial conditions
    CHARACTER(len=filename_max)   :: bc_filename           !< file with boundary conditions       
    LOGICAL                       :: debug_hd              !< additional printing for debugging
    LOGICAL                       :: diag_water_budget     !< prints to diagnose global water budget
    CHARACTER(len=:), ALLOCATABLE :: routing_scheme        !< river_routing_scheme: full / equally_to_coast / equally_to_ocean

    NAMELIST /jsb_hd_nml/          &
      active,                      &
      ic_filename, bc_filename,    &
      debug_hd,                    &
      diag_water_budget,           &
      routing_scheme

    INTEGER :: nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_hd_config'

    CALL message(TRIM(routine), 'Starting hd configuration')

    ALLOCATE(new)

    ! Set defaults
    active          = .FALSE.
    ic_filename     = 'hdstart.nc'
    bc_filename     = 'hdpara.nc'
    debug_hd        = .FALSE.
    diag_water_budget = .TRUE.
    routing_scheme  = 'full'

    nml_unit = position_nml('jsb_hd_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_hd_nml)

    new%active                    = active
    new%ic_filename               = ic_filename
    new%bc_filename               = bc_filename
    new%debug_hd                  = debug_hd
    new%diag_water_budget         = diag_water_budget
    new%routing_scheme            = routing_scheme

    config => new

  END FUNCTION new_hd_config

#endif
END MODULE mo_hd_config_class
