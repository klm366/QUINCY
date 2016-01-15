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

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish, message, message_text

#ifdef __ICON__
  USE mtime,                     ONLY: t_datetime => datetime, newDatetime, deallocateDatetime, &
    &                                  OPERATOR(+), OPERATOR(*),                                &
    &                                  divisionquotienttimespan, getDayOfYearFromDateTime,      &   
    &                                  getNoOfDaysInMonthDateTime, getNoOfDaysInYearDateTime,   & 
    &                                  getTotalMilliSecondsTimeDelta, no_of_sec_in_a_day,       &
    &                                  getNoOfSecondsElapsedInDayDateTime,                      &
    &                                  divideDatetimeDifferenceInSeconds
  USE mo_time_config,            ONLY: time_config
  USE mo_time_nml,               ONLY: read_time_namelist
  USE mo_dynamics_config,        ONLY: iequations
  USE mo_impl_constants,         ONLY: ihs_atm_temp, ihs_atm_theta, ishallow_water, &
    &                                  leapfrog_expl, leapfrog_si
  USE mo_ha_dyn_config,          ONLY: ha_dyn_config
  USE mo_bcs_time_interpolation, ONLY: tiw => current_time_interpolation_weights, &
    &                                  t_time_interpolation_weights,              &
    &                                  calculate_time_interpolation_weights
#else
  USE mo_time_conversion,        ONLY: t_datetime => time_days, &
                                       get_day_length => day_len
  USE mo_time_control,           ONLY: get_year_length => get_year_len, get_month_length => get_month_len, get_year_day, &
                                       get_date_components
  USE mo_time_control,           ONLY: start_date, stop_date, current_date, previous_date
  USE mo_interpo,                ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_semi_impl,              ONLY: eps
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
  PUBLIC :: get_time_interpolation_weights
  PUBLIC :: get_asselin_coef

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

!!$    USE mo_time_control, ONLY: lfirst_cycle, ec_manager_init, &
!!$                               init_manager, init_times, init_events, time_set
    USE mo_time_control, ONLY: time_weights

    INTEGER :: IO_timestep, istep

    ! Only called for standalone jsbach ... commented out for the time being ! TBD
!!$    CALL get_dates(IO_timestep, istep)           ! get timestep
!!$    CALL ec_manager_init(IO_timestep, istep)     ! time manager initialization
!!$    IF (lfirst_cycle) CALL init_manager
!!$    CALL init_times                              ! initialize all dates
!!$
!!$    CALL init_events
!!$    CALL time_set

    CALL time_weights(l_compute_rad=.FALSE.) ! radiation_time is not initialized yet at the time of JSBACH initialization

  END SUBROUTINE configure_time_and_events
#endif

  !>
  !! Get time step length
  !!
  !! Function returns the length of the integration time step
  !!
  REAL(wp) FUNCTION get_time_dt()

#ifndef __ICON__
    USE mo_time_control, ONLY: delta_time
#endif
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_time_dt'
    REAL(wp) :: ztime

#ifdef __ICON__
    TYPE(t_datetime), POINTER :: reference_dt
    reference_dt => newDatetime("1980-06-01T00:00:00.000")
    ztime = 0.001_wp*getTotalMilliSecondsTimeDelta(time_config%tc_dt_model, reference_dt)
    CALL deallocateDatetime(reference_dt)
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

#ifndef __ICON__
    USE mo_time_control, ONLY: start_date, stop_date
    USE mo_time_base,    ONLY: IDAYLEN
    USE mo_time_conversion, ONLY: TC_get
#endif

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_time_nsteps'

#ifndef __ICON__
    INTEGER :: start_day, start_second, stop_day, stop_second
#endif
    INTEGER :: isteps

#ifdef __ICON__
    TYPE(divisionquotienttimespan) :: tq
    
    CALL divideDatetimeDifferenceInSeconds(time_config%tc_stopdate,    &
         &                                 time_config%tc_startdate,   &
         &                                 time_config%tc_dt_model, tq)
    isteps = INT(tq%quotient)
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
    get_time_start = time_config%tc_startdate
#else
    get_time_start = start_date
#endif

  END FUNCTION get_time_start

  TYPE(t_datetime) FUNCTION get_time_stop()

#ifdef __ICON__
    get_time_stop = time_config%tc_stopdate
#else
    get_time_stop = stop_date
#endif

  END FUNCTION get_time_stop

  !>
  !! @brief Get current time
  !!
  TYPE(t_datetime) FUNCTION get_time_current()

#ifdef __ICON__
    get_time_current = time_config%tc_current_date
#else
    get_time_current = current_date
#endif

  END FUNCTION get_time_current

  !>
  !! @brief Get time at previous time step
  !!
  FUNCTION get_time_previous() RESULT(previous)
    TYPE(t_datetime), POINTER :: previous 
#ifdef __ICON__
    IF (.NOT. ASSOCIATED(previous)) previous => newDatetime("1980-06-01T00:00:00.000") 
    previous = time_config%tc_current_date + (-1) * time_config%tc_dt_model 
#else
    previous = previous_date
#endif

  END FUNCTION get_time_previous

  SUBROUTINE start_timestep(istep)

#ifndef __ICON__
    USE mo_time_control, ONLY: time_set
#endif

    INTEGER, INTENT(in) :: istep

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':start_timestep'

#ifdef __ICON__
    !
    ! Luis:
    !
    ! The following commented region is not allowed as time_config is
    ! a protected varaible. Reiner told me that this is meant for
    ! later use in JSBACH standalone.
    ! time_config%tc_current_date = time_config%tc_current_date &
    !      &                       +time_config%tc_dt_model
    !
    WRITE(message_text,'(a,i10)') 'JSBACH: Start TIME STEP n: ', istep
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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':finish_timestep'
    
#ifdef __ICON__
    WRITE(message_text,'(a,i10)') 'JSBACH: Finish TIME STEP n: ', istep
    CALL message(TRIM(routine),message_text)
#else
    CALL time_reset
#endif

  END SUBROUTINE finish_timestep

#ifdef __ICON__
  SUBROUTINE get_date_components(this_datetime, year, month, day, hour, minute, second)

    TYPE(t_datetime),  INTENT(in)  :: this_datetime
    INTEGER, OPTIONAL, INTENT(out) :: year, month, day, hour, minute, second

    IF (PRESENT(year))   year   = this_datetime%date%year
    IF (PRESENT(month))  month  = this_datetime%date%month
    IF (PRESENT(day))    day    = this_datetime%date%day
    IF (PRESENT(hour))   hour   = this_datetime%time%hour
    IF (PRESENT(minute)) minute = this_datetime%time%minute
    IF (PRESENT(second)) second = this_datetime%time%second

  END SUBROUTINE get_date_components
#endif

  !>
  !! Get length of year
  !!
  !! Function returns the length of a year in days
  !!
#ifdef __ICON__
  INTEGER FUNCTION get_year_length(yr)
    INTEGER, INTENT(in) :: yr

    TYPE(t_datetime), POINTER :: datetime_tmp

    datetime_tmp => newDatetime(yr, 1, 1, 0, 0, 0, 0)

    get_year_length = getNoOfDaysInYearDateTime(datetime_tmp)

    CALL deallocateDatetime(datetime_tmp)

  END FUNCTION get_year_length
#endif

  !>
  !! Get length of month
  !!
  !! Function returns the length of a month in days
  !!
#ifdef __ICON__
  INTEGER FUNCTION get_month_length(yr, mo)
    INTEGER, INTENT(in) :: yr, mo

    TYPE(t_datetime), POINTER :: datetime_tmp

    datetime_tmp => newDatetime(yr, mo, 1, 0, 0, 0, 0)

    get_month_length = getNoOfDaysInMonthDateTime(datetime_tmp)

    CALL deallocateDatetime(datetime_tmp)

  END FUNCTION get_month_length
#endif

  !>
  !! Get length of day
  !!
  !! Function returns the length of a day in seconds
  !!
#ifdef __ICON__
  REAL(wp) FUNCTION get_day_length()

    get_day_length = REAL(no_of_sec_in_a_day,wp)

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

    ! t_datetime points to mtimes datetime
    TYPE(t_datetime),  INTENT(in)  :: date

    get_year_day = REAL(getDayOfYearFromDateTime(date), wp) &
         &        +REAL(getNoOfSecondsElapsedInDayDateTime(date),wp) &
         &        /REAL(no_of_sec_in_a_day,wp)

  END FUNCTION get_year_day
#endif

  SUBROUTINE  get_time_interpolation_weights(w1, w2, n1, n2)

    REAL(wp), INTENT(out) :: w1, w2
    INTEGER, INTENT(out)  :: n1, n2

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':get_time_interpolation_weights'
    
#ifdef __ICON__
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

    ! LK: I think, you need to check for the 'reference date' of the interpolation weights
    IF (tiw%initialized) THEN    
      w1 = tiw%weight1
      w2 = tiw%weight2
      n1 = tiw%month1
      n2 = tiw%month2
    ELSE
      CALL message('', '')
      CALL message(routine, 'time interpolation weights not configured yet - generate new once.')
      CALL message('', '')
      current_time_interpolation_weights = calculate_time_interpolation_weights(time_config%tc_current_date)
      w1 = current_time_interpolation_weights%weight1
      w2 = current_time_interpolation_weights%weight2
      n1 = current_time_interpolation_weights%month1
      n2 = current_time_interpolation_weights%month2
    ENDIF
#else
    w1 = wgt1
    w2 = wgt2
    n1 = nmw1
    n2 = nmw2
#endif

  END SUBROUTINE get_time_interpolation_weights

  REAL(wp) FUNCTION get_asselin_coef()

#ifdef __ICON__
    SELECT CASE(iequations)
    CASE(ishallow_water,ihs_atm_temp,ihs_atm_theta)
      IF (ha_dyn_config%itime_scheme == leapfrog_expl .OR. ha_dyn_config%itime_scheme == leapfrog_si) THEN
        get_asselin_coef = ha_dyn_config%asselin_coeff
      END IF
    CASE DEFAULT
      get_asselin_coef = 0._wp
    END SELECT
#else
    get_asselin_coef = eps
#endif

  END FUNCTION get_asselin_coef

#endif
END MODULE mo_jsb_time_iface
