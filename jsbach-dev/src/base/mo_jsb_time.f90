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

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish, print_value
  USE mo_jsb_time_iface, ONLY: t_datetime, &
                               get_time_current

  IMPLICIT NONE
  PRIVATE

!!$  PUBLIC :: init_jsb_time
  PUBLIC :: init_time
  PUBLIC :: update_time_at_beginloop, update_time_at_endloop
  PUBLIC :: is_newday, get_year_at_model_start, &
            get_year_at_current_timestep, get_month_at_current_timestep, &
            get_day_at_current_timestep, get_secs_of_day_at_current_timestep

  INTEGER          :: number_timesteps
  REAL(wp)         :: timestep_length
  TYPE(t_datetime) :: start_date
  TYPE(t_datetime) :: current_date
  TYPE(t_datetime) :: final_date

!!$  !>
!!$  !! Initialize JSBACH time control
!!$  !!
!!$  INTERFACE init_jsb_time
!!$    MODULE PROCEDURE init_time_standalone
!!$    MODULE PROCEDURE init_time
!!$  END INTERFACE init_jsb_time

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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_time_standalone'

    CALL message(TRIM(routine), 'Initializing JSBACH time')

    IF (PRESENT(namelist_filename)) CALL read_time_namelist(TRIM(namelist_filename))

    CALL configure_time_and_events()

    number_timesteps  = get_time_nsteps()
    timestep_length   = get_time_dt()
    start_date        = get_time_start()
    final_date        = get_time_stop()

    current_date     = start_date

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

    USE mo_jsb_time_iface, ONLY: get_time_current, get_time_previous, &
                                 get_date_components

    TYPE(t_datetime) :: current, previous

    INTEGER :: current_day, previous_day

    current  = get_time_current()
    previous = get_time_previous()

    CALL get_date_components(current, day=current_day)
    CALL get_date_components(previous, day=previous_day)

    is_newday = current_day /= previous_day

  END FUNCTION is_newday

  INTEGER FUNCTION get_year_at_model_start()

    USE mo_jsb_time_iface, ONLY: get_date_components

    INTEGER :: year

    CALL get_date_components(start_date, year=year)

    get_year_at_model_start = year

  END FUNCTION get_year_at_model_start

  INTEGER FUNCTION get_year_at_current_timestep()

    USE mo_jsb_time_iface, ONLY: get_date_components

    TYPE(t_datetime) :: current
    INTEGER :: year

    CALL get_date_components(current, year=year)

    get_year_at_current_timestep = year

  END FUNCTION get_year_at_current_timestep

  INTEGER FUNCTION get_month_at_current_timestep()

    USE mo_jsb_time_iface, ONLY: get_date_components

    TYPE(t_datetime) :: current
    INTEGER :: month

    CALL get_date_components(current, month=month)

    get_month_at_current_timestep = month

  END FUNCTION get_month_at_current_timestep

  INTEGER FUNCTION get_day_at_current_timestep()

    USE mo_jsb_time_iface, ONLY: get_date_components

    TYPE(t_datetime) :: current
    INTEGER :: day

    CALL get_date_components(current, day=day)

    get_day_at_current_timestep = day

  END FUNCTION get_day_at_current_timestep

  INTEGER FUNCTION get_secs_of_day_at_current_timestep()

    USE mo_jsb_time_iface, ONLY: get_date_components

    TYPE(t_datetime) :: current
    INTEGER :: hour, minute, second

    CALL get_date_components(current, hour=hour, minute=minute, second=second)

    get_secs_of_day_at_current_timestep = NINT(3600._wp * REAL(hour,wp) + 60._wp * REAL(minute,wp) + REAL(second,wp))

  END FUNCTION get_secs_of_day_at_current_timestep

END MODULE mo_jsb_time
