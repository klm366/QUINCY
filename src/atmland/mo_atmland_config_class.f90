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

  TYPE, EXTENDS(t_jsb_config) :: t_atmland_config
    CHARACTER(8)     :: option_forcing
    REAL(wp)         :: lat, lon   
    CONTAINS
      PROCEDURE :: Init => Init_atmland_config
  END Type t_atmland_config

  TYPE(t_atmland_config), SAVE :: atmland_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atmland_config_class'

CONTAINS

  SUBROUTINE Init_atmland_config(config)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_atmland_config), INTENT(inout) :: config

    ! variables for reading from namlist, identical to variable-name in namelist
    !LOGICAL      :: active
    CHARACTER(8) :: create_forcing
    REAL(wp)     :: read_latitude, read_longitude

    NAMELIST /forcing_ctl/  &  
        create_forcing

    NAMELIST /grid_ctl/     &  
        read_latitude,      &
        read_longitude

    ! variables for reading model-options from namelist
    INTEGER  :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_atmland_config'

    CALL message(TRIM(routine), 'Starting atmland configuration')

    ! Set defaults
    create_forcing             = 'generate'
    read_latitude              = 50.91   ! about the MPI-BGC Jena
    read_longitude             = 11.56   ! about the MPI-BGC Jena


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
    config%option_forcing = create_forcing
    config%lat            = read_latitude
    config%lon            = read_longitude

  END SUBROUTINE Init_atmland_config

!#endif
END MODULE mo_atmland_config_class
