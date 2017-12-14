MODULE mo_rad_config_class

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_rad_config !, Get_rad_config

  TYPE, EXTENDS(t_jsb_config) :: t_rad_config
  CONTAINS
    PROCEDURE :: Init => Init_rad_config
  END type t_rad_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------

  SUBROUTINE Init_rad_config(config)

    USE mo_namelist, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_rad_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /jsb_rad_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_rad_config'


    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_rad.nc'  !R: gibts im Moment noch nicht
    ic_filename         = 'ic_land_rad.nc'  !R: gibts im Moment noch nicht

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_rad_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_rad_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename

  END SUBROUTINE Init_rad_config

END MODULE mo_rad_config_class
