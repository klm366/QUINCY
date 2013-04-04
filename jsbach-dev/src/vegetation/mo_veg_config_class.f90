!>
!! @brief Contains structures and methods for veg config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-24)
!!
MODULE mo_veg_config_class

  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_veg_config, new_veg_config

  TYPE t_veg_config
     LOGICAL :: active
     CHARACTER(len=filename_max) :: ic_filename
     CHARACTER(len=filename_max) :: bc_filename
     CHARACTER(len=10)           :: scheme_stom_cond  !< Scheme to use for computation of stomatal conductance (echam5 only for now)
     REAL(wp)                    :: wsr_capacity      !< Water holding capacity of ground skin reservoir [m water equivalent]
     REAL(wp)                    :: wsn_capacity      !< Water holding capacity of ground snow reservoir (max snow depth) [m water equivalent]
  END type t_veg_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_config_class'

CONTAINS

  FUNCTION new_veg_config(nml_unit) RESULT(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    INTEGER,          INTENT(in) :: nml_unit
    TYPE(t_veg_config), POINTER :: config

    TYPE(t_veg_config), POINTER :: new

    LOGICAL :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    CHARACTER(len=10)           :: scheme_stom_cond
    REAL(wp)                    :: wsr_capacity, wsn_capacity

    NAMELIST /jsb_veg_nml/                    &
         active,                              &
         ic_filename, bc_filename,            &
         scheme_stom_cond,                    &
         wsr_capacity, wsn_capacity

    INTEGER :: istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_veg_config'

    CALL message(TRIM(routine), 'Starting veg configuration')

    ALLOCATE(new)

    ! Set defaults
    active           = .FALSE.
    ic_filename      = 'ic_veg.nc'
    bc_filename      = 'bc_veg.nc'
    scheme_stom_cond = 'echam5'
    wsr_capacity     = 2.E-4_wp
    wsn_capacity     = 2.E-4_wp

    CALL position_nml('jsb_veg_nml', STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_veg_nml)

    new%active             = active
    new%ic_filename        = ic_filename
    new%bc_filename        = bc_filename
    new%scheme_stom_cond   = scheme_stom_cond
    new%wsr_capacity              = wsr_capacity
    new%wsn_capacity              = wsn_capacity

    config => new

  END FUNCTION new_veg_config

END MODULE mo_veg_config_class
