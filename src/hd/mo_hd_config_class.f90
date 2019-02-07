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

  USE mo_exception,         ONLY: message
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_hd_config

  TYPE, EXTENDS(t_jsb_config) :: t_hd_config
     LOGICAL                     :: debug_hd              !< additional printing for debugging
     LOGICAL                     :: diag_water_budget     !< prints to diagnose global water budget
     LOGICAL                     :: enforce_water_budget  !< finish if water conservation problem
     CHARACTER(len=32)           :: routing_scheme        !< river_routing_scheme: full / weighted_to_coast / zero
     CHARACTER(len=filename_max) :: fract_filename        !< file name with fractions (land, ocean, etc)
   CONTAINS
     PROCEDURE :: Init => Init_hd_config
  END type t_hd_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hd_config_class'

CONTAINS

  SUBROUTINE Init_hd_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND

    CLASS(t_hd_config),INTENT(inout) :: config

    LOGICAL                       :: active                !< hd module is active
    CHARACTER(len=filename_max)   :: ic_filename           !< file with initial conditions
    CHARACTER(len=filename_max)   :: bc_filename           !< file with boundary conditions
    CHARACTER(len=filename_max)   :: fract_filename        !< file with fractions
    LOGICAL                       :: debug_hd              !< additional printing for debugging
    LOGICAL                       :: diag_water_budget     !< prints to diagnose global water budget
    LOGICAL                       :: enforce_water_budget  !< finish if water conservation problem
    CHARACTER(len=32)             :: routing_scheme        !< river_routing_scheme: full / weighted_to_coast / zero

    NAMELIST /jsb_hd_nml/          &
      active,                      &
      ic_filename, bc_filename,    &
      fract_filename,              &
      debug_hd,                    &
      diag_water_budget,           &
      enforce_water_budget,        &
      routing_scheme

    INTEGER :: nml_handler, nml_unit, istat

    TYPE(t_jsb_vgrid), POINTER :: hd_o, hd_b, hd_r  ! Vertical grids

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_hd_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting hd configuration')

    ! Set defaults
    active               = .FALSE.
    ic_filename          = 'hdstart.nc'
    bc_filename          = 'bc_land_hd.nc'
    fract_filename       = 'bc_land_frac.nc'
    debug_hd             = .FALSE.
    diag_water_budget    = .TRUE.
    enforce_water_budget = .TRUE.
    routing_scheme       = 'full'

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_hd_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_hd_nml)

    CALL close_nml(nml_handler)

    config%active                    = active
    config%ic_filename               = ic_filename
    config%bc_filename               = bc_filename
    config%fract_filename            = fract_filename
    config%debug_hd                  = debug_hd
    config%diag_water_budget         = diag_water_budget
    config%enforce_water_budget      = enforce_water_budget
    config%routing_scheme            = routing_scheme

    IF (.NOT. active) RETURN

    hd_o  => new_vgrid('hd_nres_overlflow', ZAXIS_DEPTH_BELOW_LAND, 1, levels=(/1._wp/), units='')
    CALL Register_vgrid(hd_o)

    hd_b  => new_vgrid('hd_nres_baseflow', ZAXIS_DEPTH_BELOW_LAND, 1, levels=(/1._wp/), units='')
    CALL Register_vgrid(hd_b)

    hd_r  => new_vgrid('hd_nres_riverflow', ZAXIS_DEPTH_BELOW_LAND, 5, levels=(/1._wp,2._wp,3._wp,4._wp,5._wp/), units='')
    CALL Register_vgrid(hd_r)

  END SUBROUTINE Init_hd_config

#endif
END MODULE mo_hd_config_class
