!>
!! @brief Contains structures and methods for forest-age config
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-20)
!!
MODULE mo_fage_config_class
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message_text, message
  USE mo_util,              ONLY: int2string
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_fage_config

  TYPE, EXTENDS(t_jsb_config) :: t_fage_config
    INTEGER               :: &
      & max_age              !Maximum age tracked  
  CONTAINS
    PROCEDURE             :: Init => Init_fage_config
  END type t_fage_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_config_class'

CONTAINS

  SUBROUTINE Init_fage_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid, new_vgrid
    USE mo_jsb_grid,           ONLY: Register_vgrid
    USE mo_jsb_io,             ONLY: ZAXIS_DEPTH_BELOW_LAND !TODO-JN: ask Reiner why not ZAXIS_SURFACE or other

    CLASS(t_fage_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    INTEGER                     :: max_age

    NAMELIST /jsb_fage_nml/         &
      & active,                     &
      & ic_filename,                &
      & bc_filename,                &
      & max_age

    INTEGER :: nml_handler, nml_unit, istat
    INTEGER :: i
    REAL(wp), ALLOCATABLE :: age_fract_init(:)

    TYPE(t_jsb_vgrid), POINTER :: forest_age ! grid cell fraction per age

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_fage_config'

    IF (debug_on()) CALL message(TRIM(routine), 'Starting fage configuration')

    ! Set defaults
    active            = .TRUE.
    bc_filename       = 'bc_land_fage.nc'  !JN: TODO do I need it at all?
    ic_filename       = 'ic_land_fage.nc'  !JN: TODO
    max_age           = 1 !In case of pfts always 1, only in case of cohorts >> 1 (e.g. 250)

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_fage_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_fage_nml)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    config%max_age             = max_age

    ALLOCATE(age_fract_init(max_age))

!TODO-JN: check max_age > 0 < 1000

    CALL message(TRIM(routine), 'Maximum tracked age: '//TRIM(int2string(max_age)))

    CALL close_nml(nml_handler)

    IF (.NOT. active) RETURN

    DO i=1,max_age
      age_fract_init = REAL(i)
    END DO

    forest_age  => new_vgrid('forest_age_fract', ZAXIS_DEPTH_BELOW_LAND, max_age, &
      & levels  = age_fract_init, &
!      & levels  = REAL((/ i, 1:max_age /),wp),                         &
      & units='y')
    CALL register_vgrid(forest_age)

    WRITE(message_text, *) 'Tracked forest ages [y]: ', max_age
    CALL message(TRIM(routine), message_text)

    DEALLOCATE(age_fract_init)

  END SUBROUTINE Init_fage_config

#endif
END MODULE mo_fage_config_class
