!>
!! @brief Contains structures and methods for radiation config
!! 
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!! 
!! @par Revision History
!! First version by Rainer Schneck (2015-04-27)
!! implementation to SLM Steffen Richter (2017-12-14)
!! 
MODULE mo_rad_config_class

  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_rad_config !, Get_rad_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the radiation process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_rad_config
    !
    ! no additional variables yet (additional to the ones defined in t_jsb_config)
    !
  CONTAINS
    PROCEDURE :: Init => Init_rad_config
  END type t_rad_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_rad_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_rad_config(config)

    USE mo_namelist,          ONLY: open_nml, POSITIONED, position_nml, close_nml

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_rad_config), INTENT(inout) :: config
    ! ---------------------------
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    NAMELIST /radiation_ctl/      &
         active,                      &
         ic_filename,                 &
         bc_filename

    ! variables for reading model-options from namelist
    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_rad_config'

    CALL message(TRIM(routine), 'Starting radiation configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_rad.nc'  !R: gibts im Moment noch nicht
    ic_filename         = 'ic_land_rad.nc'  !R: gibts im Moment noch nicht

    ! read the namelist
    ! &RADIATION_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('RADIATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, radiation_ctl )

    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename

  END SUBROUTINE Init_rad_config

END MODULE mo_rad_config_class
