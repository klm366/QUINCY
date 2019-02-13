!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief Contains methods for JSBACH master control.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-02-22)
!!
MODULE mo_jsb_control
#ifndef __NO_JSBACH__

  USE mo_io_units,  ONLY: filename_max

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_runs_standalone
  PUBLIC :: debug_on
  PUBLIC :: get_no_of_models
  PUBLIC :: timer_aggregate, timer_integrate, l_timer
  PUBLIC :: max_no_of_models

  LOGICAL :: is_standalone
  INTEGER :: debug = 1
  LOGICAL :: l_timer

  INTEGER, PARAMETER    :: max_no_of_models = 10
  INTEGER :: no_of_models = 1

  INTEGER, DIMENSION(max_no_of_models) :: timer_aggregate, timer_integrate

  CHARACTER(len=filename_max) :: master_namelist_filename = 'jsbach_master.namelist'

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_control'

CONTAINS

  ! SUBROUTINE init_jsb_master_control(namelist_filename)
  ! 
  !   CHARACTER(len=*), OPTIONAL, INTENT(in) :: namelist_filename
  ! 
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_jsb_master_control'
  ! 
  !   IF (PRESENT(namelist_filename)) master_namelist_filename = TRIM(namelist_filename)
  ! 
  !   CALL read_jsb_control_namelist(TRIM(master_namelist_filename))
  ! 
  ! END SUBROUTINE init_jsb_master_control
  ! 
  ! SUBROUTINE read_jsb_control_namelist(namelist_filename)
  ! 
  !   USE mo_namelist, ONLY: open_nml, POSITIONED, position_nml, close_nml
  ! 
  !   CHARACTER(len=*), INTENT(in) :: namelist_filename
  ! 
  !   INTEGER :: nml_handler, nml_unit, istat, i
  !   LOGICAL :: rwnd
  ! 
  !   CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_jsb_control_namelist'
  ! 
  !   NAMELIST /jsb_control_nml/ &
  !     is_standalone,           &
  !     debug
  ! 
  !   ! Namelist defaults
  !   is_standalone  = .TRUE.
  !   debug          = 0
  ! 
  !   nml_handler = open_nml(TRIM(namelist_filename))
  ! 
  !   nml_unit = position_nml('jsb_control_nml', nml_handler, STATUS=istat)
  !   IF (istat == POSITIONED) READ(nml_unit, jsb_control_nml)
  ! 
  !   CALL close_nml(nml_handler)
  ! 
  ! END SUBROUTINE read_jsb_control_namelist


  LOGICAL FUNCTION jsbach_runs_standalone()

    jsbach_runs_standalone = is_standalone

  END FUNCTION jsbach_runs_standalone

  INTEGER FUNCTION get_no_of_models()

    get_no_of_models = no_of_models

  END FUNCTION get_no_of_models

  LOGICAL FUNCTION debug_on(level)

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: level

    IF (PRESENT(level)) THEN
      SELECT CASE (TRIM(level))
      CASE ('basic')
        debug_on = debug >= 1
      CASE ('detail')
        debug_on = debug >= 2
      CASE ('hsm')
        debug_on = debug >= 3
      END SELECT
    ELSE
      debug_on = debug > 0
    END IF

  END FUNCTION debug_on

#endif
END MODULE mo_jsb_control
