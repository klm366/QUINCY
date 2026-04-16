!> Contains routines of the climbuf processes
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

MODULE mo_climbuf_process
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_climbuf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_process'
CONTAINS
  !
  ! ======================================================================================================== !
  !>
  !>#### Compute running averages of climate variables
  !>
  !> We here calculate daily, annual or climatological mean, maximum or minimum values from
  !> atmospheric forcing data.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE SUBROUTINE calc_climbuf (lstart, init_running_means, new_day, new_month, new_year, &
    &                                     dtime, t_air, rain, snow, &
    &                                     seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
    &                                     min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp_20y, max_mmtemp_20y, &
    &                                     temp_15d, temp_15d_yDay, precip_15d, precip_15d_yDay, &
    &                                     gdd5_sum_year, gdd5_prev_year)

    LOGICAL, INTENT(in) :: &
      & lstart,                 & !< This is the first time step in a simulation
      & init_running_means,     & !< Initialize the min/max temperature climatology
      & new_day,                & !< This is the first time step of the day
      & new_month,              & !< This is the first time step of the month
      & new_year                  !< This is the first time step of the year

    REAL(wp), INTENT(in) :: &
      & dtime,                  & !< Time step length [s]
      & seconds_day,            & !< Length of a day [s]
      & seconds_month,          & !< Length of the previous month [s]
      & t_air,                  & !< Near surface air temperature [K]
      & rain,                   & !< Liquid precipitation [kg m-2 s-1]
      & snow                      !< Solid precipitation [kg m-2 s-1]

    REAL(wp), INTENT(inout) :: &
      & temp_sum_day,           & !< Current temperature sum of the day [K]
      & temp_sum_month,         & !< Current temperature sum of the month [K]
      & min_mmtemp_of_yr,       & !< Current minimum monthly mean temperature of the year [K]
      & max_mmtemp_of_yr,       & !< Current maximum monthly mean temperature of the year [K]
      & min_mmtemp_20y,         & !< Climatological mean temperature of the coldest month [K]
      & max_mmtemp_20y,         & !< Climatological mean temperature of the warmest month [K]
      & temp_15d,               & !< 15-day smoothed temperature [K]
      & temp_15d_yDay,          & !< 15-day smoothed temperature at the end of the previous day [K]
      & precip_15d,             & !< 15-day smoothed precipitation [K]
      & precip_15d_yDay,        & !< 15-day smoothed precipitation at the end of the previous day [K]
      & gdd5_sum_year,          & !< Current number of growing degree days of this year
      & gdd5_prev_year            !< Growing degree days of previous year

    ! Local variables
    REAL(wp) :: &
      & prev_day_mean_temp,     & !< Air temperature of previous day [K]
      & prev_month_mean_temp,   & !< Air temperature of previous month [K]
      & F_15d,                  & !< Factor for 15-day smoothing
      & N_15d                     !< Factor for 15-day smoothing

    ! ------------------------------------------------------------------------------------------------------ !


    IF (.NOT. new_day .OR. lstart) THEN
      temp_sum_day = temp_sum_day + t_air * dtime
    ELSE
      prev_day_mean_temp = temp_sum_day / seconds_day
      temp_sum_day = t_air * dtime
    ENDIF

    IF (.NOT. new_month .OR. lstart) THEN
      temp_sum_month = temp_sum_month + t_air * dtime
    ELSE
      prev_month_mean_temp = temp_sum_month / seconds_month
      temp_sum_month = t_air * dtime
    END IF ! new_month

    !-- Growing degree days (base 5 degC)
    !  ---------------------
    IF (new_year) THEN
      gdd5_prev_year = gdd5_sum_year
      gdd5_sum_year = 0._wp
    END IF
    IF (new_day .AND. .NOT. lstart) THEN
      gdd5_sum_year = gdd5_sum_year + MAX(0._wp, prev_day_mean_temp - tmelt - 5._wp)
    END IF

    !-- Multi-year minimum and maximum monthly mean temperatures
    !  ----------------------
    IF (new_month .AND. .NOT. lstart) THEN
      ! Calculate coldest/warmest month of the year so far
      min_mmtemp_of_yr = MIN(min_mmtemp_of_yr, prev_month_mean_temp)
      max_mmtemp_of_yr = MAX(max_mmtemp_of_yr, prev_month_mean_temp)

      IF (new_year) THEN
        ! --- build 20yr running mean of coldest and warmest monthly temperature
        IF (init_running_means) THEN
          min_mmtemp_20y = min_mmtemp_of_yr
          max_mmtemp_20y = max_mmtemp_of_yr
        ELSE
          min_mmtemp_20y = (min_mmtemp_20y * 19._wp + min_mmtemp_of_yr) / 20._wp
          max_mmtemp_20y = (max_mmtemp_20y * 19._wp + max_mmtemp_of_yr) / 20._wp
        ENDIF
        ! --- reset for current year
        min_mmtemp_of_yr =  1000._wp
        max_mmtemp_of_yr = -1000._wp
      END IF ! new_year
    END IF

    !-- 15-day smoothing of temperature and precipitation (needed with yasso)
    !  ---------------------
    F_15d = EXP(-dtime / (15._wp * seconds_day))
    N_15d = 1._wp - F_15d

    temp_15d   = t_air         * N_15d  +  F_15d * temp_15d
    precip_15d = (rain + snow) * N_15d  +  F_15d * precip_15d
    ! TODO: Check, if these yDay variables are actually needed
    IF (new_day .AND. .NOT. lstart) THEN
      ! Save previous day variables (for yasso)
      temp_15d_yDay   = temp_15d
      precip_15d_yDay = precip_15d
    END IF


  END SUBROUTINE calc_climbuf

#endif
END MODULE mo_climbuf_process
