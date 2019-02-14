!>
!! @brief Contains interface to time control
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-02-06)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!! Added ltimer although this doesn't really belong here      by Reiner Schnur (2018-06-04)
!!
MODULE mo_jsb_time_iface
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel, mpi_comm

#ifdef __ICON__
  USE mtime,                     ONLY: t_datetime => datetime, newDatetime, deallocateDatetime,      &
    &                                  OPERATOR(+), OPERATOR(*), OPERATOR(==),                       &
    &                                  OPERATOR(<=), OPERATOR(>), OPERATOR(-),                       &
    &                                  divisionquotienttimespan, getDayOfYearFromDateTime,           &
    &                                  getNoOfDaysInMonthDateTime, getNoOfDaysInYearDateTime,        &
    &                                  getTotalMilliSecondsTimeDelta, no_of_sec_in_a_day,            &
    &                                  getNoOfSecondsElapsedInDayDateTime, getTotalSecondsTimeDelta, &
    &                                  divideDatetimeDifferenceInSeconds, isCurrentEventActive
  USE mo_time_config,            ONLY: time_config !configure_time
  USE mo_time_nml,               ONLY: read_time_namelist
  USE mo_dynamics_config,        ONLY: iequations
  USE mo_impl_constants,         ONLY: ihs_atm_temp, ihs_atm_theta, ishallow_water, &
    &                                  leapfrog_expl, leapfrog_si
  USE mo_ha_dyn_config,          ONLY: ha_dyn_config
  USE mo_bcs_time_interpolation, ONLY: t_time_interpolation_weights,              &
    &                                  calculate_time_interpolation_weights
  USE mo_echam_phy_config,       ONLY: echam_phy_tc, dt_zero
  USE mo_master_config,          ONLY: isRestart
  USE mo_run_config,             ONLY: l_timer_host => ltimer
#else
  USE mo_time_conversion, ONLY: t_datetime => time_days, &
                                get_day_length => day_len
  USE mo_time_control,    ONLY: get_year_length => get_year_len, get_month_length => get_month_len, get_year_day, &
                                get_date_components
  USE mo_time_control,    ONLY: start_date, stop_date, current_date, previous_date, next_date, lstart, lresume, l_trigradm1
  USE mo_interpo,         ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_semi_impl,       ONLY: eps
  USE mo_control,         ONLY: l_timer_host => ltimer

#ifdef STANDALONE
  USE mo_time_event,      ONLY: io_time_event, TIME_INC_DAYS, TIME_INC_MONTHS, TRIG_FIRST, TRIG_LAST
  USE mo_time_control,    ONLY: lresume, dt_start, dt_resume, dt_stop, delta_time, no_cycles, &
                                no_days, no_steps, putdata, putrerun, trigfiles, p_bcast_event

  USE mo_jsb_namelist_iface, ONLY: open_nml, POSITIONED, close_nml, position_nml
#endif

#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_datetime
  PUBLIC :: get_year_length, get_month_length, get_day_length, get_year_day
  PUBLIC :: get_time_dt, get_time_nsteps, &
            get_time_start, get_time_stop, &
            get_time_current, get_time_previous, is_time_experiment_start, is_time_restart, &
            is_time_ltrig_rad_m1
  PUBLIC :: read_time_namelist !, configure_time
  PUBLIC :: start_timestep, finish_timestep
  PUBLIC :: configure_time_and_events
  PUBLIC :: get_date_components
  PUBLIC :: get_time_interpolation_weights
  PUBLIC :: get_asselin_coef
  PUBLIC :: l_timer_host

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time_iface'

CONTAINS

#ifdef __ICON__
  SUBROUTINE configure_time_and_events(model_shortname)

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname

    ! CALL configure_time
    RETURN

  END SUBROUTINE configure_time_and_events
#else
  SUBROUTINE read_time_namelist(filename, restart_jsbach)

    CHARACTER(LEN=*), INTENT(in) :: filename
    LOGICAL, INTENT(IN) :: restart_jsbach

    CHARACTER(len=*), PARAMETER :: routine = modname//':read_time_namelist'

#ifdef STANDALONE
    INTEGER :: nml_handler, nml_unit, istat

    NAMELIST /jsb_time_nml/ &
      no_cycles,            &
      dt_start,             &
      dt_resume,            &
      dt_stop,              &
      no_days,              &
      no_steps,             &
      delta_time,           &
      putdata,              &
      putrerun,             &
      trigfiles

    !-----------------------------------------------------------------------------------------------

    !NOTE: lresume is not read here, but was already read as 'restart_jsbach' from the jsb_control_nml
    ! in read_jsb_control_namelist. However, it needs to be set here to properly work with echam
    IF(restart_jsbach) THEN
      lresume = .TRUE.
    ELSE
      lresume = .FALSE.
    ENDIF

    IF(my_process_is_stdio()) THEN
      ! default values
      no_cycles = 1
      dt_start(:) = 0
      dt_resume(:) = 0
      dt_stop(:) = 0
      no_days = -1
      no_steps = -1
      delta_time = 0
      putdata  = io_time_event(1, TIME_INC_DAYS, TRIG_FIRST, 0)
      putrerun = io_time_event(1, TIME_INC_MONTHS, TRIG_LAST, 0)
      trigfiles = io_time_event(1, TIME_INC_MONTHS, TRIG_FIRST, 0)

      nml_handler = open_nml(TRIM(filename))

      nml_unit =  position_nml('jsb_time_nml', nml_handler, STATUS=istat)
      IF (istat == POSITIONED) READ(nml_unit, jsb_time_nml)

      CALL close_nml(nml_handler)

      !Assertion: in case of lresume = .TRUE. dt_start should not be specified, because this
      !   can lead to problemes (diverging results) if it does not match what echam calculates from
      !   fdate/vdate i.e. do not specifiy dt_start if running with lresume but let echam calculate it!
      IF ((lresume) .AND. ANY(dt_start(:) .NE. 0)) CALL finish(routine, 'If lresume=true dt_start should not be specified' &
        &//' (it seems to lead to inconsistencies in the echam event management).')

    ENDIF

    IF (my_process_is_mpi_parallel()) THEN
    ! jsb3 text: The following are used from ECHAM5 modules and ... need
    ! to be broadcast if running in standalone mode so that the modules can use them
      CALL p_bcast(no_cycles, p_io, comm=mpi_comm)
      CALL p_bcast(dt_start, p_io, comm=mpi_comm)
      CALL p_bcast(dt_resume, p_io, comm=mpi_comm)
      CALL p_bcast(dt_stop, p_io, comm=mpi_comm)
      CALL p_bcast(no_days, p_io, comm=mpi_comm)
      CALL p_bcast(no_steps, p_io, comm=mpi_comm)
      CALL p_bcast(delta_time, p_io, comm=mpi_comm)
      CALL p_bcast_event(putdata, p_io, comm=mpi_comm)
      CALL p_bcast_event(trigfiles, p_io, comm=mpi_comm)
      CALL p_bcast_event(putrerun, p_io, comm=mpi_comm)
    ENDIF
#endif

  END SUBROUTINE read_time_namelist

  SUBROUTINE configure_time_and_events(model_shortname)

#ifdef STANDALONE
    USE mo_time_control, ONLY: lfirst_cycle, ec_manager_init, &
                               init_manager, init_times, init_events, time_set
#endif
    USE mo_time_control, ONLY: time_weights

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname

    CHARACTER(len=*), PARAMETER :: routine = modname//':configure_time_and_events'
    INTEGER :: IO_timestep, istep

#ifdef STANDALONE
    !Assertion:
    IF (.NOT. PRESENT(model_shortname)) CALL finish(routine, &
      & 'Short name missing! In standalone model, a model_short name needs to be passed, please check!')

    CALL get_dates(model_shortname, IO_timestep, istep)           ! get timestep
    CALL ec_manager_init(IO_timestep, istep)     ! time manager initialization
    IF (lfirst_cycle) CALL init_manager
    CALL init_times                              ! initialize all dates
    CALL init_events
#endif
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
    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_dt'
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

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_nsteps'

#ifndef __ICON__
    INTEGER :: start_day, start_second, stop_day, stop_second, idt
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

    idt = INT(get_time_dt())
    isteps =                                   &
      & (( IDAYLEN/idt * (stop_day - start_day - 1))   & ! Number of seconds in full days
      &    + (stop_second                              & ! Number of seconds in last day
      &       + (IDAYLEN - start_second - 1)           & ! Number of seconds in first day
      &      ) / idt)
#endif

!@todo: error message is misleading! should be: error in istep calculation or comparable!?
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
    TYPE(t_datetime) :: previous
#ifdef __ICON__
    previous = time_config%tc_current_date + (-1) * time_config%tc_dt_model
#else
    previous = previous_date
#endif

  END FUNCTION get_time_previous

  !>
  !! @brief Get time at next time step
  !!
  FUNCTION get_time_next() RESULT(next)
    TYPE(t_datetime) :: next
#ifdef __ICON__
    next = time_config%tc_current_date + (1) * time_config%tc_dt_model
#else
    next = next_date
#endif

  END FUNCTION get_time_next

  !
  ! Check whether we are one time step before radiation is calculated in atmosphere
  !
  FUNCTION is_time_ltrig_rad_m1(model_id) RESULT(ltrig_rad_m1)

    INTEGER, INTENT(in) :: model_id
    LOGICAL :: ltrig_rad_m1

#ifdef __ICON__
    !
    ! In the ICON MPI atmosphere, radiation is triggered based on the radiation time delta and the old  time step.
    ! We therefore use the current time step to check whether we're one time step before the radiation is triggered
    ! in the atmosphere
    TYPE(t_datetime) :: current
    INTEGER          :: dt_rad

    current = get_time_current()
    IF (echam_phy_tc(model_id)%dt_rad > dt_zero) THEN
      dt_rad = getTotalSecondsTimeDelta(echam_phy_tc(model_id)%dt_rad, current)
      ltrig_rad_m1 = MOD(getNoOfSecondsElapsedInDayDateTime(current), dt_rad) == 0
    ELSE
      ltrig_rad_m1 = .FALSE.
    END IF

#else
    ltrig_rad_m1 = l_trigradm1
#endif

  END FUNCTION is_time_ltrig_rad_m1

  LOGICAL FUNCTION is_time_experiment_start()

#ifdef __ICON__
    ! In ICON, the first computed time step is one time step after the start date
    is_time_experiment_start = get_time_previous() == time_config%tc_exp_startdate
#else
    ! In ECHAM, the first computed time step is the start date
    is_time_experiment_start = lstart
#endif

  END FUNCTION is_time_experiment_start

  LOGICAL FUNCTION is_time_restart()

#ifdef __ICON__
    ! In ICON, the first computed time step is one time step after the restart date
    is_time_restart = .NOT. is_time_experiment_start() .AND. (get_time_previous() == time_config%tc_startdate)
#else
    is_time_restart = lresume
#endif

  END FUNCTION is_time_restart

  SUBROUTINE start_timestep(istep)

#ifndef __ICON__
    USE mo_time_control, ONLY: time_set
#endif

    INTEGER, INTENT(in) :: istep

    CHARACTER(len=*), PARAMETER :: routine = modname//':start_timestep'

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':finish_timestep'

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

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_time_interpolation_weights'

#ifdef __ICON__
    TYPE(t_time_interpolation_weights) :: current_time_interpolation_weights

!    CALL message('', '')
!    CALL message(routine, 'generate time interpolation weights.')
!    CALL message('', '')
    current_time_interpolation_weights = calculate_time_interpolation_weights(time_config%tc_current_date)
    w1 = current_time_interpolation_weights%weight1
    w2 = current_time_interpolation_weights%weight2
    n1 = current_time_interpolation_weights%month1_index
    n2 = current_time_interpolation_weights%month2_index
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


#if ! defined(__ICON__) && defined(STANDALONE)
  !>
  !! @brief Get dates for jsbach standalone
  !!
  !! Currently only working for ECHAM
  !!
  !! @par History
  !! JN: Initial version copied from jsbach3 mo_jsbach_interface.f90
  !!
  SUBROUTINE get_dates(model_shortname, IO_timestep, istep)

    ! For standalone mode only.
    ! Corresponds to ECHAM5 subroutine IO_init in mo_io.f90
    ! This routine doesn't need to be called from ECHAM in coupled mode.

    USE mo_filename, ONLY: out_expname
    USE mo_io, ONLY: IO_READ, IO_open, IO_close
    USE mo_netCDF, ONLY: FILE_INFO, NF_GLOBAL, io_get_att_int
    USE mo_time_control, ONLY: INIT_STEP, delta_time, dt_start, &
                               lresume, resume_date, start_date, inp_convert_date, write_date
    USE mo_time_conversion, ONLY : time_native, TC_set, TC_convert

    USE mo_mpi,              ONLY : p_start, p_bcast, p_parallel_io, p_parallel, p_io

    CHARACTER(len=*), INTENT(IN) :: model_shortname
    INTEGER, INTENT(out) :: IO_timestep, istep

    TYPE(FILE_INFO) :: IO_file
    INTEGER :: IO_file_id
    INTEGER :: forecast_date, verification_date   ! YYYYMMDD
    INTEGER :: forecast_time, verification_time   ! HHMMSS
    TYPE(time_native) :: date_nat

    ! Initialize time
    IF (p_parallel_io) THEN
       IF (lresume) THEN
         ! copied from jsb3 get_dates:
         ! Get time information from restart file
          IO_file%opened = .FALSE.
          CALL IO_open('restart_'//TRIM(out_expname)//'_'//TRIM(model_shortname)//'_01.nc', IO_file, IO_READ)

          IO_file_id = IO_file%file_id
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'nstep', istep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'timestep', IO_timestep)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'fdate', forecast_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'ftime', forecast_time)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vdate', verification_date)
          CALL io_get_att_int(IO_file_id, NF_GLOBAL, 'vtime', verification_time)
          CALL IO_close(IO_file)
       ELSE
          istep = INIT_STEP
          IO_timestep = delta_time
       ENDIF
    ENDIF

    IF (p_parallel) THEN
       CALL p_bcast(istep, p_io, comm=mpi_comm)
       CALL p_bcast(IO_timestep, p_io, comm=mpi_comm)
       IF (lresume) THEN
          CALL p_bcast(forecast_date, p_io, comm=mpi_comm)
          CALL p_bcast(forecast_time, p_io, comm=mpi_comm)
          CALL p_bcast(verification_date, p_io, comm=mpi_comm)
          CALL p_bcast(verification_time, p_io, comm=mpi_comm)
       ENDIF
    ENDIF

    IF (lresume) THEN
       CALL inp_convert_date(forecast_date, forecast_time, start_date)
       CALL inp_convert_date(verification_date, verification_time, resume_date)
    ELSE
       IF (SUM(dt_start(:)) /= 0) THEN
          CALL TC_set(dt_start(1),dt_start(2),dt_start(3), &
               dt_start(4),dt_start(5),dt_start(6), date_nat)
          CALL TC_convert(date_nat, start_date)
          resume_date = start_date
       ELSE
          CALL finish('mo_jsb_time_iface', 'Start date not set')
       ENDIF
    ENDIF

    IF (p_parallel_io) THEN
       CALL write_date(start_date, 'Start date (initial/restart): ')
       CALL write_date(resume_date, 'Resume date (initial/restart): ')
    ENDIF

  END SUBROUTINE get_dates
#endif

#endif
END MODULE mo_jsb_time_iface
