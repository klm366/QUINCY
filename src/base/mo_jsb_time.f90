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
  USE mo_exception,      ONLY: message, message_text, finish
  USE mo_jsb_time_iface, ONLY: t_datetime, get_date_components, get_time_current, get_time_start, &
                               get_time_dt, is_time_experiment_start, is_time_restart, &
                               get_year_length, get_month_length, get_day_length, get_year_day,   &
                               is_time_ltrig_rad_m1, get_time_interpolation_weights, get_asselin_coef

  IMPLICIT NONE
  PRIVATE

!!$  PUBLIC :: init_jsb_time
  PUBLIC :: init_time
  PUBLIC :: update_time_at_beginloop, update_time_at_endloop
  PUBLIC :: t_datetime
  PUBLIC :: get_date_components, get_time_current, get_time_start, get_time_dt, &
            is_time_experiment_start, is_time_restart, &
            get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: is_newday, is_newyear, get_year_at_model_start, is_time_ltrig_rad_m1, &
            get_year, get_month, &
            get_day, get_secs_of_day
  PUBLIC :: get_time_interpolation_weights
  PUBLIC :: get_asselin_coef
  PUBLIC :: timeStep_in_days, time_steps_per_day

  TYPE(t_datetime), SAVE :: start_date
  TYPE(t_datetime), SAVE :: current_date
  TYPE(t_datetime), SAVE :: final_date

  LOGICAL, SAVE :: time_initialized = .FALSE.

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time'

CONTAINS

  SUBROUTINE init_time(namelist_filename, model_shortname, restart_jsbach)

    USE mo_jsb_time_iface, ONLY: read_time_namelist, & !, configure_time
                                 configure_time_and_events, &
                                 get_time_start, get_time_stop

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: namelist_filename
    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname
    LOGICAL, OPTIONAL, INTENT(IN) :: restart_jsbach

    REAL(wp) :: eps

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_time'

    !TODO? currently the time is initalised from the namelist of the first processed model
    !      but then used globally for all models!
    !      If different models should have different time information this would need to change
    !      -> Attention: it is not possible to specify this global information before knowing the names
    !         of the models: knowing the models is a precondition for restart with echam, because echam restart files
    !        contain the short name of a model. The time-information of a restart file in turn needs
    !        to be processed before times and events can be initialised.
    IF (.NOT. time_initialized) THEN
      CALL message(TRIM(routine), 'Initializing JSBACH time')

      ! TBD: remove #ifdef
#ifdef __ICON__
    IF (PRESENT(restart_jsbach)) CONTINUE ! avoid compiler warnings about dummy arguments not being used
    IF (PRESENT(namelist_filename)) THEN
      CALL read_time_namelist(TRIM(namelist_filename))
#else
      IF (PRESENT(namelist_filename) .AND. PRESENT(restart_jsbach)) THEN
        CALL read_time_namelist(TRIM(namelist_filename), restart_jsbach)
#endif
      ENDIF

      IF (PRESENT(model_shortname)) THEN
        CALL configure_time_and_events(model_shortname)
      ELSE
        CALL configure_time_and_events()
      ENDIF

      start_date       = get_time_start()
      final_date       = get_time_stop()
      current_date     = start_date

      ! just to avoid compiler warnings about variable set but never referenced
      IF (.FALSE.) THEN
        PRINT*, start_date%date%year
        PRINT*, final_date%date%year
        PRINT*, current_date%date%year
      END IF

      ! Check if Asselin time filter is used and print message
      eps = get_asselin_coef()
      IF (eps > 0._wp) THEN
        WRITE(message_text,*) 'Using Asselin time filter for surface temperature: eps =', eps
      ELSE
        WRITE(message_text,*) 'Not using Asselin time filter for surface temperature'
      END IF
      CALL message(TRIM(routine), TRIM(message_text))

      time_initialized = .TRUE.
    ENDIF

  END SUBROUTINE init_time


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

    ! CALL get_date_components(current, day=current_day)
    ! CALL get_date_components(previous, day=previous_day)
    current_day = INT(get_year_day(current))
    previous_day = INT(get_year_day(previous))

    is_newday = current_day /= previous_day

  END FUNCTION is_newday

  LOGICAL FUNCTION is_newyear()

    USE mo_jsb_time_iface, ONLY: get_time_previous

    TYPE(t_datetime) :: current, previous

    INTEGER :: current_day, previous_day

    current  = get_time_current()
    previous = get_time_previous()

    ! CALL get_date_components(current, year=current_year)
    ! CALL get_date_components(previous, year=previous_year)
    current_day = INT(get_year_day(current))
    previous_day = INT(get_year_day(previous))

    is_newyear = current_day == 1 .AND. previous_day > 300

  END FUNCTION is_newyear


  INTEGER FUNCTION get_year(this_datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: this_datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: year

    IF (PRESENT(this_datetime)) THEN
      datetime_tmp = this_datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF
    CALL get_date_components(datetime_tmp, year=year)

    get_year = year

  END FUNCTION get_year

  INTEGER FUNCTION get_month(this_datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: this_datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: month

    IF (PRESENT(this_datetime)) THEN
      datetime_tmp = this_datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF

    CALL get_date_components(datetime_tmp, month=month)

    get_month = month

  END FUNCTION get_month

  INTEGER FUNCTION get_day(this_datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: this_datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: day

    IF (PRESENT(this_datetime)) THEN
      datetime_tmp = this_datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF

    CALL get_date_components(datetime_tmp, day=day)

    get_day = day

  END FUNCTION get_day

  INTEGER FUNCTION get_secs_of_day(this_datetime)

    TYPE(t_datetime), OPTIONAL, INTENT(in) :: this_datetime

    TYPE(t_datetime) :: datetime_tmp
    INTEGER :: hour, minute, second

    IF (PRESENT(this_datetime)) THEN
      datetime_tmp = this_datetime
    ELSE
      datetime_tmp  = get_time_current()
    END IF

    CALL get_date_components(datetime_tmp, hour=hour, minute=minute, second=second)

    get_secs_of_day = NINT(3600._wp * REAL(hour,wp) + 60._wp * REAL(minute,wp) + REAL(second,wp))

  END FUNCTION get_secs_of_day

  INTEGER FUNCTION get_year_at_model_start()

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_year_at_model_start'

    IF (.NOT. time_initialized) &
      CALL finish(TRIM(routine), 'Time not initialized yet.')

    get_year_at_model_start = get_year(start_date)

  END FUNCTION get_year_at_model_start


  REAL(wp) FUNCTION timeStep_in_days() ! EXPRESS the time step in days (as fraction of day)
    USE mo_jsb_time_iface,   ONLY: get_time_dt

    REAL(wp)                 :: dtime

    dtime = get_time_dt()
    timeStep_in_days = dtime / 86400._wp

  END FUNCTION timeStep_in_days


  INTEGER FUNCTION time_steps_per_day() ! NUMBER of time steps per day
    USE mo_jsb_time_iface,   ONLY: get_time_dt

    REAL(wp)                 :: dtime

    dtime = get_time_dt()

    ! Set time_steps_per_day
    IF( mod(86400,INT(dtime)) /= 0) THEN ! For computing the mean day temperature it is assumed that each day
                                         ! is computed with a fixed number of time steps. Therefore the program is
                                         ! stopped wheen day length is not an integer multiple of the time step.
       CALL finish("function time_steps_per_day","ERROR: Day length is not an integer multiple of the time step!")
    ELSE
       time_steps_per_day = 86400/int(dtime)
    END IF

  END FUNCTION time_steps_per_day


#endif
END MODULE mo_jsb_time
