!>
!! @brief Contains interface to time control
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-02-06)
!!
MODULE mo_jsb_time_iface
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish, message, message_text

#ifdef __ICON__
  USE mo_datetime,    ONLY: t_datetime, print_datetime, add_time, rdaylen
  USE mo_time_config, ONLY: time_config !configure_time
  USE mo_time_nml,    ONLY: read_time_namelist
  USE mo_amip_bc,     ONLY: wgt1, wgt2, nmw1,nmw2  ! TBD: this should be made more general
#else
  USE mo_time_conversion, ONLY: t_datetime => time_days, &
                                get_day_length => day_len
  USE mo_time_control,    ONLY: get_year_length => get_year_len, get_month_length => get_month_len, get_year_day, &
                                get_date_components
  USE mo_time_control,    ONLY: start_date, stop_date, current_date, previous_date
  USE mo_interpo,         ONLY: wgt1, wgt2, nmw1, nmw2
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_datetime
  PUBLIC :: get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: get_time_dt, get_time_nsteps, &
            get_time_start, get_time_stop, &
            get_time_current, get_time_previous
  PUBLIC :: read_time_namelist !, configure_time
  PUBLIC :: start_timestep, finish_timestep
  PUBLIC :: configure_time_and_events
  PUBLIC :: get_date_components
  PUBLIC :: wgt1, wgt2, nmw1, nmw2

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time_iface'

CONTAINS

#ifdef __ICON__
  SUBROUTINE configure_time_and_events()

    ! CALL configure_time
    RETURN

  END SUBROUTINE configure_time_and_events
#else
  SUBROUTINE read_time_namelist(filename)

    CHARACTER(LEN=*), INTENT(in) :: filename

  END SUBROUTINE read_time_namelist

  SUBROUTINE configure_time_and_events

    USE mo_time_control, ONLY: lfirst_cycle, ec_manager_init, &
                               init_manager, init_times, init_events, time_set

    INTEGER :: IO_timestep, istep

    ! Only called for standalone jsbach ... commented out for the time being ! TBD
!!$    CALL get_dates(IO_timestep, istep)           ! get timestep
!!$    CALL ec_manager_init(IO_timestep, istep)     ! time manager initialization
!!$    IF (lfirst_cycle) CALL init_manager
!!$    CALL init_times                              ! initialize all dates
!!$
!!$    CALL init_events

    CALL time_set

  END SUBROUTINE configure_time_and_events
#endif

  !>
  !! Get time step length
  !!
  !! Function returns the length of the integration time step
  !!
  REAL(wp) FUNCTION get_time_dt()

#ifdef __ICON__
    USE mo_run_config, ONLY: dtime
#else
    USE mo_time_control, ONLY: delta_time
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_time_dt'

    REAL(wp) :: ztime

#ifdef __ICON__
    ztime = dtime
#else
    ztime = delta_time
#endif

    IF (ztime <= 0.0_wp) &
      CALL finish(routine, 'time step not configured yet.')

    get_time_dt = ztime

  END FUNCTION get_time_dt

  !>
  !! Get number of time steps
  !!
  !! Function returns the number of time steps for current model integration
  !!
  INTEGER FUNCTION get_time_nsteps()

#ifdef __ICON__
    USE mo_run_config,   ONLY: nsteps
#else
    USE mo_time_control, ONLY: start_date, stop_date
    USE mo_time_base,    ONLY: IDAYLEN
    USE mo_time_conversion, ONLY: TC_get
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_time_nsteps'

    INTEGER :: start_day, start_second, stop_day, stop_second
    INTEGER :: isteps

#ifdef __ICON__
    isteps = nsteps
#else
    CALL TC_get(start_date, start_day, start_second)
    CALL TC_get(stop_date, stop_day, stop_second)

    isteps =                                   &
      & ( IDAYLEN * (stop_day - start_day - 1) & ! Number of seconds in full days
      &   + stop_second                        & ! Number of seconds in last day
      &   + (IDAYLEN - start_second - 1)       & ! Number of seconds in first day
      & ) / INT(get_time_dt())
#endif

    IF (isteps <= 0) &
      CALL finish(routine, 'number of time steps not configured yet.')

    get_time_nsteps = isteps

  END FUNCTION get_time_nsteps

  TYPE(t_datetime) FUNCTION get_time_start()

#ifdef __ICON__
    get_time_start = time_config%ini_datetime
#else
    get_time_start = start_date
#endif

  END FUNCTION get_time_start

  TYPE(t_datetime) FUNCTION get_time_stop()

#ifdef __ICON__
    get_time_stop = time_config%end_datetime
#else
    get_time_stop = stop_date
#endif

  END FUNCTION get_time_stop

  !>
  !! @brief Get current time
  !!
  TYPE(t_datetime) FUNCTION get_time_current()

#ifdef __ICON__
    get_time_current = time_config%cur_datetime
#else
    get_time_current = current_date
#endif

  END FUNCTION get_time_current

  !>
  !! @brief Get time at previous time step
  !!
  TYPE(t_datetime) FUNCTION get_time_previous()

#ifdef __ICON__
    USE mo_run_config, ONLY: dtime

    TYPE(t_datetime) :: datetime

    datetime = time_config%cur_datetime
    CALL add_time(-dtime, 0, 0, 0, datetime)
    get_time_previous = datetime
#else
    get_time_previous = previous_date
#endif

  END FUNCTION get_time_previous

  SUBROUTINE start_timestep(istep)

#ifdef __ICON__
    USE mo_run_config, ONLY: dtime
#else
    USE mo_time_control, ONLY: time_set
#endif

    INTEGER, INTENT(in) :: istep

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':start_timestep'

#ifdef __ICON__
    CALL add_time(dtime, 0, 0, 0, time_config%cur_datetime)

    WRITE(message_text,'(a,i10)') 'TIME STEP n: ', istep
    CALL message(TRIM(routine),message_text)
#else
    CALL time_set
#endif

  END SUBROUTINE start_timestep

  SUBROUTINE finish_timestep(istep)

#ifndef __ICON__
    USE mo_time_control, ONLY: time_reset
#endif

    INTEGER, INTENT(in) :: istep

#ifndef __ICON__
    CALL time_reset
#endif

  END SUBROUTINE finish_timestep

#ifdef __ICON__
  SUBROUTINE get_date_components(datetime, year, month, day, hour, minute, second)

    USE mo_datetime, ONLY: time_to_date, check_date

    TYPE(t_datetime),  INTENT(in)  :: datetime
    INTEGER, OPTIONAL, INTENT(out) :: year, month, day, hour, minute, second

    TYPE(t_datetime) :: datetime_tmp

    datetime_tmp = datetime
    CALL time_to_date(datetime_tmp)
    CALL check_date(datetime_tmp)

    IF (PRESENT(year))   year   = datetime_tmp%year
    IF (PRESENT(month))  month  = datetime_tmp%month
    IF (PRESENT(day))    day    = datetime_tmp%day
    IF (PRESENT(hour))   hour   = datetime_tmp%hour
    IF (PRESENT(minute)) minute = datetime_tmp%minute
    IF (PRESENT(second)) second = INT(datetime_tmp%second)

  END SUBROUTINE get_date_components
#endif

  !>
  !! Get length of year
  !!
  !! Function returns the length of a year in days
  !!
#ifdef __ICON__
  INTEGER FUNCTION get_year_length(yr)

    USE mo_datetime, ONLY: aux_datetime

    INTEGER, INTENT(in) :: yr

    TYPE(t_datetime) :: date

    date%year   = yr
    date%month  = 1
    date%day    = 1
    date%hour   = 0
    date%minute = 0
    date%second = 0

    CALL aux_datetime(date)

    get_year_length = date%yealen

  END FUNCTION get_year_length
#endif

  !>
  !! Get length of month
  !!
  !! Function returns the length of a month in days
  !!
#ifdef __ICON__
  INTEGER FUNCTION get_month_length(yr, mo)

    USE mo_datetime, ONLY: aux_datetime

    INTEGER, INTENT(in) :: yr, mo

    TYPE(t_datetime) :: date

    date%year   = yr
    date%month  = mo
    date%day    = 1
    date%hour   = 0
    date%minute = 0
    date%second = 0

    CALL aux_datetime(date)

    get_month_length = date%monlen

  END FUNCTION get_month_length
#endif

  !>
  !! Get length of day
  !!
  !! Function returns the length of a day in seconds
  !!
#ifdef __ICON__
  REAL(wp) FUNCTION get_day_length()

    get_day_length = rdaylen

  END FUNCTION get_day_length
#endif

  !>
  !! Get day in year
  !!
  !! Function returns the day of the year for a given date.
  !! The seconds in the day are given as fraction of the day. 
  !!
#ifdef __ICON__
  REAL(wp) FUNCTION get_year_day(date)

    USE mo_datetime, ONLY: aux_datetime

    TYPE(t_datetime),  INTENT(in)  :: date

    TYPE(t_datetime) :: datetime_tmp

    datetime_tmp = date

    CALL aux_datetime(datetime_tmp)

    get_year_day = REAL(datetime_tmp%yeaday, wp) + datetime_tmp%daytim

  END FUNCTION get_year_day
#endif

#endif
END MODULE mo_jsb_time_iface
