!>
!! @brief Contains structures and methods for JSBACH model config
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_config_class

  USE mo_srf_config_class,  ONLY: t_srf_config
  USE mo_soil_config_class, ONLY: t_soil_config
  USE mo_veg_config_class,  ONLY: t_veg_config
  USE mo_io_units,          ONLY: filename_max
  USE mo_exception,         ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_model_config, new_model_config

  TYPE t_jsb_model_config
    CHARACTER(len=filename_max)  :: grid_filename
    CHARACTER(len=30)            :: grid_name
    CHARACTER(len=20)            :: usecase
    TYPE(t_srf_config),  POINTER :: srf
    TYPE(t_soil_config), POINTER :: soil
    TYPE(t_veg_config),  POINTER :: veg
  END TYPE t_jsb_model_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_config_class'

CONTAINS

  FUNCTION new_model_config(namelist_filename) RESULT(config)

    USE mo_srf_config_class,   ONLY: new_srf_config
    USE mo_soil_config_class,  ONLY: new_soil_config
    USE mo_veg_config_class,   ONLY: new_veg_config
    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CHARACTER(len=*), INTENT(in) :: namelist_filename
    TYPE(t_jsb_model_config), POINTER :: config

    TYPE(t_jsb_model_config), POINTER :: new

    CHARACTER(len=filename_max) :: grid_filename
    CHARACTER(len=30)           :: grid_name

    NAMELIST /jsb_grid_nml/   &
      grid_filename,          &
      grid_name

    CHARACTER(len=20)           :: usecase

    NAMELIST /jsb_model_nml/  &
      usecase

    INTEGER :: nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':new_model_config'

    CALL message(TRIM(routine), 'Starting model configuration from '//TRIM(namelist_filename))

    ALLOCATE(new)

    ! Set defaults
    grid_filename = ''
    grid_name     = ''
    usecase       = ''

    ! Read namelist
    nml_unit = open_nml(TRIM(namelist_filename))

    CALL position_nml('jsb_grid_nml', STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_grid_nml)

    new%grid_filename = grid_filename
    new%grid_name     = grid_name

    CALL position_nml('jsb_model_nml', STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_model_nml)

    new%usecase       = usecase

    new%srf  => new_srf_config(nml_unit)

    new%soil => new_soil_config(nml_unit)

    new%veg  => new_veg_config(nml_unit)

    CALL close_nml(nml_unit)

    config => new

  END FUNCTION new_model_config

END MODULE mo_jsb_config_class
