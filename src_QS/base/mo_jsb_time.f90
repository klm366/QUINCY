!> Contains time control for JSBACH
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_time
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_jsb_math_constants,  ONLY: one_day
  USE mo_jsb_time_iface,      ONLY: is_time_experiment_start, is_time_ltrig_rad_m1, &
    &                               t_datetime, deallocateDatetime, &
    &                               get_date_components, &
    &                               get_year_day, &
    &                               get_year_length, get_month_length, get_day_length, &
    &                               get_time_start, &
    &                               get_asselin_coef

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: is_time_experiment_start, is_time_ltrig_rad_m1, is_newday, is_newyear
  PUBLIC :: t_datetime, deallocateDatetime
  PUBLIC :: get_year, get_year_day
  PUBLIC :: get_year_length, get_month_length, get_day_length
  PUBLIC :: get_time_start
  PUBLIC :: get_asselin_coef
  PUBLIC :: timesteps_per_day
  PUBLIC :: timestep_in_days

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_time'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> returns true if the current and previous time step are from different days
  LOGICAL FUNCTION is_newday(current, dt)

    USE mo_jsb_time_iface, ONLY: get_time_previous

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_day, previous_day

    previous => get_time_previous(current, dt)

    ! CALL get_date_components(current, day=current_day)
    ! CALL get_date_components(previous, day=previous_day)
    current_day = INT(get_year_day(current))
    previous_day = INT(get_year_day(previous))

    is_newday = current_day /= previous_day

    CALL deallocateDatetime(previous)

  END FUNCTION is_newday

  ! ======================================================================================================= !
  !>
  !> returns true if the current and previous time step are from different years
  LOGICAL FUNCTION is_newyear(current, dt)

    USE mo_jsb_time_iface, ONLY: get_time_previous

    TYPE(t_datetime), POINTER, INTENT(in) :: current
    REAL(wp),                  INTENT(in) :: dt

    TYPE(t_datetime), POINTER :: previous

    INTEGER :: current_year, previous_year

    previous => get_time_previous(current, dt)

    ! CALL get_date_components(current, year=current_year)
    ! CALL get_date_components(previous, year=previous_year)
    current_year  = get_year(current)
    previous_year = get_year(previous)

    is_newyear = current_year /= previous_year

    CALL deallocateDatetime(previous)

  END FUNCTION is_newyear

  ! ======================================================================================================= !
  !>
  !> returns the year of the provided datetime variable
  INTEGER FUNCTION get_year(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: year

    CALL get_date_components(this_datetime, year=year)

    get_year = year

  END FUNCTION get_year

  ! ======================================================================================================= !
  !>
  !> returns the day of the provided datetime variable
  INTEGER FUNCTION get_day(this_datetime)

    TYPE(t_datetime), POINTER, INTENT(in) :: this_datetime

    INTEGER :: day

    CALL get_date_components(this_datetime, day=day)

    get_day = day

  END FUNCTION get_day

  ! ======================================================================================================= !
  !>
  !> returns number of time steps per day
  INTEGER FUNCTION timesteps_per_day(dt)

    REAL(wp), INTENT(in) :: dt

    CHARACTER(len=*), PARAMETER :: routine = modname//':timesteps_per_day'

    timesteps_per_day = NINT(one_day / dt)

    IF (ABS(timesteps_per_day * dt - one_day) > EPSILON(1._wp)) &
      & CALL finish(routine, 'Time step must fit into day length')

  END FUNCTION timesteps_per_day

  ! ======================================================================================================= !
  !>
  !> returns time step in days, i.e., as fraction of day
  REAL(wp) FUNCTION timestep_in_days(dt, model_id) ! EXPRESS the time step in days (as fraction of day)

    REAL(wp), OPTIONAL, INTENT(in) :: dt
    INTEGER,  OPTIONAL, INTENT(in) :: model_id

    REAL(wp) :: dt_loc

    CHARACTER(len=*), PARAMETER :: routine = modname//':timestep_in_days'

    IF (PRESENT(dt)) THEN
      dt_loc = dt
    ! ELSE IF (PRESENT(model_id)) THEN
    !   dt_loc = get_time_dt(model_id)
    ELSE
      CALL finish(routine, 'Must provide either dt oder model_id')
    END IF

    timestep_in_days = dt_loc / 86400._wp

  END FUNCTION timestep_in_days

#endif
END MODULE mo_jsb_time
