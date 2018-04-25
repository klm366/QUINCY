!>
!! @brief Contains structures and methods for fuel config
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (06.04.2017)
!! 2018-04-25, StW: Separating namelist in one general and one for each disturbance algorithm
!!
MODULE mo_fuel_config_class
#ifndef __NO_JSBACH__
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_io_units,          ONLY: filename_max
  USE mo_kind,              ONLY: wp
  USE mo_jsb_config_class,  ONLY: t_jsb_config

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_fuel_config

  TYPE, EXTENDS(t_jsb_config) :: t_fuel_config
    INTEGER :: fuel_algorithm  ! Select fuel algorithm: 0 : none, 1 : jsbach (default), 2 : Spitfire (Thonicke et al.(2001))
  CONTAINS
     PROCEDURE :: Init => Init_fuel_config
  END type t_fuel_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize fuel process
  !!
  !! @param[inout]     config     Configuration type of process (t_fuel_config)
  !!
  SUBROUTINE Init_fuel_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml
    USE mo_fuel_constants,     ONLY: fuel_names, FUEL_NONE, FUEL_JSB3, FUEL_SPITFIRE
    USE mo_input_strings,      ONLY: GetStringIndex

    CLASS(t_fuel_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename

    CHARACTER(len=16)           :: fuel_algorithm
 
    NAMELIST /jsb_fuel_nml/     &
         active,                &
         ic_filename,           &
         bc_filename,           &
         fuel_algorithm

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Init_fuel_config'


    ! Set defaults
    active                 = .FALSE.
    bc_filename            = 'bc_land_fuel.nc'
    ic_filename            = 'ic_land_fuel.nc'
    fuel_algorithm         = 'jsb3'

    ! Read namelist
    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_fuel_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) THEN
      READ(nml_unit, jsb_fuel_nml)
      CALL message(TRIM(routine), 'Obtained fuel configuration from: '//TRIM(config%namelist_filename))
    ENDIF

    ! Write resulting values in config memory
    config%active                 = active
    config%ic_filename            = ic_filename
    config%bc_filename            = bc_filename
    config%fuel_algorithm         = GetStringIndex(TRIM(fuel_algorithm),fuel_names,ErrMsg='fuel algorithm name') - 1

    ! Read algorithm specfic namelists here
    SELECT CASE (config%fuel_algorithm)
      CASE (FUEL_NONE)
        ! By definition nothing to do
      CASE (FUEL_JSB3)
        CALL Init_fuel_jsb3_config(config,nml_handler)
      CASE (FUEL_SPITFIRE)
        CALL Init_fuel_spitfire_config(config,nml_handler)
    END SELECT

    CALL close_nml(nml_handler)

  END SUBROUTINE Init_fuel_config

  SUBROUTINE Init_fuel_jsb3_config(config,nml_handler)
    CLASS(t_fuel_config), INTENT(inout) :: config
    INTEGER,              INTENT(in)    :: nml_handler
  END SUBROUTINE Init_fuel_jsb3_config

  SUBROUTINE Init_fuel_spitfire_config(config,nml_handler)
    CLASS(t_fuel_config), INTENT(inout) :: config
    INTEGER,              INTENT(in)    :: nml_handler
  END SUBROUTINE Init_fuel_spitfire_config

#endif
END MODULE mo_fuel_config_class
