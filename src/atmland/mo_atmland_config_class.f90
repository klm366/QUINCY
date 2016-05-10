!>
!! @brief Contains structures and methods for atmland config
!!
!! @author
!!  Jan Engel, MPI-BGC
!!
!! @par Revision History
!! First version by Jan Engel (2016-04-07)
!! based on the mo_veg_config_class, which was initially developed by Reiner Schnur MPI-M Hamburg
!!
MODULE mo_atmland_config_class

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: atmland_config

  !-----------------------------------------------------------------------------------------------------
  !> configuration of the atmland process, derived from t_jsb_config
  !! 
  !! currently it does mainly: reading parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_config) :: t_atmland_config
    CHARACTER(8)     :: option_forcing
    REAL(wp)         :: lat, lon   
    CONTAINS
      PROCEDURE :: Init => Init_atmland_config
  END Type t_atmland_config

  TYPE(t_atmland_config), SAVE :: atmland_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_config_class'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> configuration routine of t_atmland_config
  !! 
  !! currently it does only: read parameters from namelist
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Init_atmland_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml
    ! 0.1 InOut
    CLASS(t_atmland_config), INTENT(inout) :: config    !< config type for atmland
    ! 0.2 Local
    ! variables for reading from namlist, identical to variable-name in namelist
    !LOGICAL      :: active
    CHARACTER(8) :: option_forcing
    REAL(wp)     :: latitude, longitude

    NAMELIST /forcing_ctl/  &  
        option_forcing

    NAMELIST /grid_ctl/     &  
        latitude,           &
        longitude

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_atmland_config'

    CALL message(TRIM(routine), 'Starting atmland configuration')

    ! Set defaults
    option_forcing        = 'generate'
    latitude              = 50.91   ! about the MPI-BGC Jena
    longitude             = 11.56   ! about the MPI-BGC Jena


    ! read the namelist
    ! 
    ! &FORCING_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('FORCING_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, forcing_ctl)
    CALL close_nml(nml_handler)
    ! &GRID_CTL
    nml_handler = open_nml(TRIM(config%namelist_filename))
    nml_unit = position_nml('GRID_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, grid_ctl)
    CALL close_nml(nml_handler)

    ! pass values as read from file
    config%option_forcing = option_forcing
    config%lat            = latitude
    config%lon            = longitude

  END SUBROUTINE Init_atmland_config

!#endif
END MODULE mo_atmland_config_class
