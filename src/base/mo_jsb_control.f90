!>
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
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_jsb_master_control
  PUBLIC :: jsbach_runs_standalone, jsbach_is_restarted
  PUBLIC :: debug_on
  PUBLIC :: model_base_dir
  ! t_jsb_model_nml needs to be public for the NEC compiler
  PUBLIC :: t_jsb_model_nml, jsb_models_nml, get_no_of_models

  LOGICAL :: is_standalone
  LOGICAL :: restart_jsbach
  CHARACTER(len=99) :: model_base_dir
  INTEGER :: debug

  TYPE t_jsb_model_nml
    CHARACTER(len=30)  :: model_name
    CHARACTER(len=10)  :: model_shortname
    CHARACTER(len=132) :: model_description
    CHARACTER(len=filename_max) :: model_namelist_filename
  END TYPE t_jsb_model_nml

  INTEGER, PARAMETER    :: max_no_of_models = 2
  INTEGER :: no_of_models = 0
  TYPE(t_jsb_model_nml) :: jsb_models_nml(max_no_of_models)

  CHARACTER(len=filename_max) :: master_namelist_filename = 'jsbach_master.namelist'

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_control'

CONTAINS

  SUBROUTINE init_jsb_master_control(namelist_filename)

    USE mo_jsb_time, ONLY: init_time

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: namelist_filename

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_jsb_master_control'

    IF (PRESENT(namelist_filename)) master_namelist_filename = TRIM(namelist_filename)

    CALL read_jsb_control_namelist(TRIM(master_namelist_filename))

    IF (is_standalone) THEN
       CALL init_time(TRIM(master_namelist_filename))
    ELSE
       CALL init_time()
    END IF

  END SUBROUTINE init_jsb_master_control

  !>
  !! @brief Reads JSBACH control parameters.
  !!
  !! Reads namelist 'jsb_control_nml' and namelists 'jsb_model_nml' from namelist file.
  !!
  !! @par History
  !! Ideas taken from ICON mo_master_nml.f90 by Rene Redler.
  !!
  SUBROUTINE read_jsb_control_namelist(namelist_filename)

    USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, position_nml, close_nml

    CHARACTER(len=*), INTENT(in) :: namelist_filename

    CHARACTER(len=30)  :: model_name
    CHARACTER(len=10)  :: model_shortname
    CHARACTER(len=132) :: model_description
    CHARACTER(len=filename_max) :: model_namelist_filename

    INTEGER :: nml_handler, nml_unit, istat
    LOGICAL :: rwnd

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_jsb_control_namelist'

    NAMELIST /jsb_control_nml/ &
      is_standalone,           &
      restart_jsbach,          &
      model_base_dir,          &
      debug

    NAMELIST /jsb_model_nml/   &
      model_name,              &
      model_shortname,         &
      model_description,       &
      model_namelist_filename

    ! Namelist defaults
    is_standalone  = .TRUE.
    restart_jsbach = .FALSE.
    model_base_dir = '.'
    debug          = 0

    nml_handler = open_nml(TRIM(namelist_filename))

    nml_unit = position_nml('jsb_control_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_control_nml)

    rwnd = .TRUE.
    no_of_models = 0
    DO
      nml_unit = position_nml('jsb_model_nml', nml_handler, lrewind=rwnd, STATUS=istat)
      IF (istat /= POSITIONED) EXIT
      IF (no_of_models >= max_no_of_models) THEN
        CALL finish(routine, 'no_of_models >= max_no_of_models')
      ENDIF
      rwnd=.false.

      model_name = ''
      model_shortname = ''
      model_description = ''
      model_namelist_filename = ''
      READ(nml_unit, jsb_model_nml)

      no_of_models = no_of_models + 1
      jsb_models_nml(no_of_models)%model_name              = model_name
      jsb_models_nml(no_of_models)%model_shortname         = model_shortname
      jsb_models_nml(no_of_models)%model_description       = model_description
      jsb_models_nml(no_of_models)%model_namelist_filename = model_namelist_filename
    END DO

    CALL close_nml(nml_handler)

  END SUBROUTINE read_jsb_control_namelist

  LOGICAL FUNCTION jsbach_runs_standalone()

    jsbach_runs_standalone = is_standalone

  END FUNCTION jsbach_runs_standalone

  LOGICAL FUNCTION jsbach_is_restarted()

    ! TBD: get this from parent model if not standalone
    jsbach_is_restarted = restart_jsbach

  END FUNCTION jsbach_is_restarted

  INTEGER FUNCTION get_no_of_models()

    get_no_of_models = no_of_models

  END FUNCTION get_no_of_models

  LOGICAL FUNCTION debug_on()

    debug_on = debug > 0

  END FUNCTION debug_on

#endif
END MODULE mo_jsb_control
