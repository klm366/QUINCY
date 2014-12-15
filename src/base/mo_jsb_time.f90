!>
!! @brief Contains time control for JSBACH
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-02-18)
!!
MODULE mo_jsb_time
#ifndef __NO_JSBACH__

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: message, finish
  USE mo_jsb_time_iface, ONLY: t_datetime, get_date_components, get_time_current, get_time_start, &
                               get_year_length, get_month_length, get_day_length, get_year_day,   &
                               get_time_interpolation_weights

  IMPLICIT NONE
  PRIVATE

!!$  PUBLIC :: init_jsb_time
  PUBLIC :: init_time
  PUBLIC :: update_time_at_beginloop, update_time_at_endloop
  PUBLIC :: t_datetime
  PUBLIC :: get_date_components, get_time_current, get_time_start, &
            get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: is_newday, get_year_at_model_start, &
            get_year, get_month, &
            get_day, get_secs_of_day
  PUBLIC :: get_time_interpolation_weights

  INTEGER                :: number_timesteps
  REAL(wp)               :: timestep_length
  TYPE(t_datetime), SAVE :: start_date
  TYPE(t_datetime), SAVE :: current_date
  TYPE(t_datetime), SAVE :: final_date

!!$  !>
!!$  !! Initialize JSBACH time control
!!$  !!
!!$  INTERFACE init_jsb_time
!!$    MODULE PROCEDURE init_time_standalone
!!$    MODULE PROCEDURE init_time
!!$  END INTERFACE init_jsb_time

  LOGICAL, SAVE :: time_initialized = .FALSE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time'

CONTAINS

  !>
  !! @brief Initialize time control for JSBACH standalone
  !!
!!$  SUBROUTINE init_time_standalone(namelist_filename)
  SUBROUTINE init_time(namelist_filename)

    USE mo_jsb_time_iface, ONLY: read_time_namelist, & !, configure_time
                                 configure_time_and_events, &
                                 get_time_nsteps, get_time_dt, &
                                 get_time_start, get_time_stop

    CHARACTER(len=*), OPTIONAL, INTENT(in) :: namelist_filename

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_time'

    CALL message(TRIM(routine), 'Initializing JSBACH time')

    IF (PRESENT(namelist_filename)) CALL read_time_namelist(TRIM(namelist_filename))

    CALL configure_time_and_events()

    number_timesteps  = get_time_nsteps()
    timestep_length   = get_time_dt()
    start_date        = get_time_start()
    final_date        = get_time_stop()

    current_date     = start_date

    time_initialized = .TRUE.

  END SUBROUTINE init_time

!!$  END SUBROUTINE init_time_standalone

!!$  !>
!!$  !! @brief Initialize JSBACH time control from other model
!!$  !!
!!$  SUBROUTINE init_time(nsteps, ini_datetime, end_datetime, dtime)
!!$
!!$    INTEGER,          INTENT(in) :: nsteps
!!$    TYPE(t_datetime), INTENT(in) :: ini_datetime
!!$    TYPE(t_datetime), INTENT(in) :: end_datetime
!!$    REAL(wp),         INTENT(in) :: dtime
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_time'
!!$
!!$    number_timesteps = nsteps
!!$    timestep_length  = dtime
!!$    start_date       = ini_datetime
!!$    final_date       = end_datetime
!!$
!!$    current_date     = start_date
!!$
!!$  END SUBROUTINE init_time

  SUBROUTINE update_time_at_beginloop(istep, is_standalone)

    USE mo_jsb_time_iface, ONLY: start_timestep

    INTEGER, INTENT(in) :: istep
    LOGICAL, INTENT(in) :: is_standalone

    IF (is_standalone) CALL start_timestep(istep)

  END SUBROUTINE update_time_at_beginloop

  SUBROUTINE update_time_at_endloop(istep, is_standalone)

    USE mo_jsb_time_iface, ONLY: finish_timestep

    INTEGER, INTENT(in) :: istep
    LOGICAL, INTENT(in) :: is_standalone

    IF (is_standalone) CALL finish_timestep(istep)

  END SUBROUTINE update_time_at_endloop

  LOGICAL FUNCTION is_newday()

    USE mo_jsb_time_iface, ONLY: get_time_previous

    TYPE(t_datetime) :: current, previous

    INTEGER :: current_day, previous_day

    current  = get_time_current()
    previous = get_time_previous()

    CALL get_date_components(current, day=current_day)
    CALL get_date_components(previous, day=previous_day)

    is_newday = current_day /= previous_day

  END FUNCTION is_newday

  INTEGER FUNCTION get_year(datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: year

    IF (PRESENT(datetime)) THEN
      datetime_tmp = datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF
    CALL get_date_components(datetime_tmp, year=year)

    get_year = year

  END FUNCTION get_year

  INTEGER FUNCTION get_month(datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: month

    IF (PRESENT(datetime)) THEN
      datetime_tmp = datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF

    CALL get_date_components(datetime_tmp, month=month)

    get_month = month

  END FUNCTION get_month

  INTEGER FUNCTION get_day(datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: day

    IF (PRESENT(datetime)) THEN
      datetime_tmp = datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF

    CALL get_date_components(datetime_tmp, day=day)

    get_day = day

  END FUNCTION get_day

  INTEGER FUNCTION get_secs_of_day(datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: hour, minute, second

    IF (PRESENT(datetime)) THEN
      datetime_tmp = datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF

    CALL get_date_components(datetime_tmp, hour=hour, minute=minute, second=second)

    get_secs_of_day = NINT(3600._wp * REAL(hour,wp) + 60._wp * REAL(minute,wp) + REAL(second,wp))

  END FUNCTION get_secs_of_day

  INTEGER FUNCTION get_year_at_model_start()

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_year_at_model_start'

    IF (.NOT. time_initialized) &
      CALL finish(TRIM(routine), 'Time not initialized yet.')

    get_year_at_model_start = get_year(start_date)

  END FUNCTION get_year_at_model_start

#endif
END MODULE mo_jsb_time
