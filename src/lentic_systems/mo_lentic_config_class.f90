! ===============================================================================================================================
! === THIS IS A TEMPLATE. PLEASE REPLACE ...                                                                                  ===
! === ... <FIRST NAME LAST NAME>    with your name e.g. "Martin Mustermann"                                                   ===
! === ... <DATE>                    with the present date in YYYY-MM-DD format e.g. "2016-04-30"                              ===
! === ... lentic with your process name e.g. "seb" for surface energy balance                              ===
! ===                                                                                                                         ===
! === Then go through the code line by line to adapt it to your needs:                                                        ===
! === !X marks lines with examples (e.g:), templates (if necessary:) and implementation points (Implementation:) that you     ===
! === have to adapt to your process.                                                                                          ===
! ===                                                                                                                         ===
! === Finally delete this header.                                                                                             ===
! ===============================================================================================================================


!>
!! @brief Contains structures and methods for lentic config
!!
!! @author
!!  <FIRST NAME LAST NAME>, MPI-M Hamburg
!!
!! @par Revision History
!! First version by <FIRST NAME LAST NAME> (<DATE>)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_lentic_config_class
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
  !PUBLIC :: t_lentic_config

  !TYPE, EXTENDS(t_jsb_config) :: t_lentic_config
     !X e.g: LOGICAL                     :: use_alb_veg_simple

  !  CONTAINS
  !    PROCEDURE :: Init => Init_lentic_config
  ! END type t_lentic_config

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lentic_config_class'

CONTAINS
  ! -------------------------------------------------------------------------------------------------------
  !> Initialize lentic process
  !!
  !! @param[inout]     config     Configuration type of process (t_lentic_config)
  !!  
  SUBROUTINE Init_lentic_config(config)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CLASS(t_lentic_config), INTENT(inout) :: config

    LOGICAL                     :: active
    CHARACTER(len=filename_max) :: ic_filename, bc_filename
    !X e.g: LOGICAL                     :: use_alb_veg_simple

    NAMELIST /jsb_lentic_nml/      &
         active,                      &
         ic_filename,                 &
         bc_filename,                 &
         !X e.g: use_lentic_veg_simple

    INTEGER :: nml_handler, nml_unit, istat

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_lentic_config'

    CALL message(TRIM(routine), 'Starting lentic configuration')

    ! Set defaults
    active              = .TRUE.
    bc_filename         = 'bc_land_lentic.nc'
    ic_filename         = 'ic_land_lentic.nc'
    !X e.g: use_alb_veg_simple        = .FALSE.

    nml_handler = open_nml(TRIM(config%namelist_filename))

    nml_unit = position_nml('jsb_lentic_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_lentic_nml)

    CALL close_nml(nml_handler)

    config%active              = active
    config%ic_filename         = ic_filename
    config%bc_filename         = bc_filename
    !X e.g: config%use_alb_veg_simple        = use_alb_veg_simple
    
  END SUBROUTINE Init_lentic_config

#endif
END MODULE mo_lentic_config_class
