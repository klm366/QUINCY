!>
!! @brief Contains interface to time control
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version                                              by Reiner Schnur (2013-02-06)
!! New subroutines and some adaptations for jsb4 standalone   by Julia Nabel (2016-02-20)
!!
MODULE mo_jsb_time_iface
#ifndef __NO_JSBACH__

  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, message, message_text
  USE mo_jsb_parallel_iface, ONLY: my_process_is_stdio, p_io, p_bcast, my_process_is_mpi_parallel

#ifdef __ICON__
  USE mo_datetime,        ONLY: t_datetime, print_datetime, add_time, rdaylen
  USE mo_time_config,     ONLY: time_config !configure_time
  USE mo_time_nml,        ONLY: read_time_namelist
  USE mo_time_interpolation_weights, ONLY: wi_limm
  USE mo_dynamics_config, ONLY: iequations
  USE mo_impl_constants,  ONLY: ihs_atm_temp, ihs_atm_theta, ishallow_water, leapfrog_expl, leapfrog_si
  USE mo_ha_dyn_config,   ONLY: ha_dyn_config
#else
  USE mo_time_conversion, ONLY: t_datetime => time_days, &
                                get_day_length => day_len
  USE mo_time_control,    ONLY: get_year_length => get_year_len, get_month_length => get_month_len, get_year_day, &
                                get_date_components
  USE mo_time_control,    ONLY: start_date, stop_date, current_date, previous_date
  USE mo_interpo,         ONLY: wgt1, wgt2, nmw1, nmw2
  USE mo_semi_impl,       ONLY: eps

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
  SUBROUTINE configure_time_and_events(model_shortname)

    CHARACTER(len=*), OPTIONAL, INTENT(IN) :: model_shortname

    ! CALL configure_time
    RETURN

  END SUBROUTINE configure_time_and_events
#else
  SUBROUTINE read_time_namelist(filename, restart_jsbach)

    CHARACTER(LEN=*), INTENT(in) :: filename
    LOGICAL, INTENT(IN) :: restart_jsbach

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_time_namelist'

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
      CALL p_bcast(no_cycles, p_io)
      CALL p_bcast(dt_start, p_io)
      CALL p_bcast(dt_resume, p_io)
      CALL p_bcast(dt_stop, p_io)
      CALL p_bcast(no_days, p_io)
      CALL p_bcast(no_steps, p_io)
      CALL p_bcast(delta_time, p_io)
      CALL p_bcast_event(putdata, p_io)
      CALL p_bcast_event(trigfiles, p_io)
      CALL p_bcast_event(putrerun, p_io)
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

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':configure_time_and_events'
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

    USE mo_datetime, ONLY: date_to_time

    INTEGER, INTENT(in) :: yr

    TYPE(t_datetime) :: date, inidate

    inidate = get_time_start()

    date%calendar = inidate%calendar
    date%year     = yr
    date%month    = 1
    date%day      = 1
    date%hour     = 0
    date%minute   = 0
    date%second   = 0._wp

    CALL date_to_time(date)

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

    USE mo_datetime, ONLY: date_to_time

    INTEGER, INTENT(in) :: yr, mo

    TYPE(t_datetime) :: date, inidate

    inidate = get_time_start()

    date%calendar = inidate%calendar
    date%year     = yr
    date%month    = mo
    date%day      = 1
    date%hour     = 0
    date%minute   = 0
    date%second   = 0

    CALL date_to_time(date)

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

  SUBROUTINE  get_time_interpolation_weights(w1, w2, n1, n2)

    REAL(wp), INTENT(out) :: w1, w2
    INTEGER, INTENT(out)  :: n1, n2

#ifdef __ICON__
    w1 = wi_limm%wgt1
    w2 = wi_limm%wgt2
    n1 = wi_limm%inm1
    n2 = wi_limm%inm2
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
          CALL IO_open('restart_'//TRIM(out_expname)//'_'//TRIM(model_shortname)//'.nc', IO_file, IO_READ)

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
       CALL p_bcast(istep, p_io)
       CALL p_bcast(IO_timestep, p_io)
       IF (lresume) THEN
          CALL p_bcast(forecast_date, p_io)
          CALL p_bcast(forecast_time, p_io)
          CALL p_bcast(verification_date, p_io)
          CALL p_bcast(verification_time, p_io)
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
