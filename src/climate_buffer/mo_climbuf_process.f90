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

  USE mo_jsb_control, ONLY: debug_on
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
  ELEMENTAL PURE SUBROUTINE calc_climbuf (lstart, init_running_means, new_day, new_month, new_year, delta_time, t_air, &
                                          seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
                                          min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20, &
                                          gdd_sum_year, gdd_prev_year)

    LOGICAL, INTENT(in) :: &
      & lstart,                 & !< This is the first time step in a simulation
      & init_running_means,     & !< Initialize the min/max temperature climatology
      & new_day,                & !< This is the first time step of the day
      & new_month,              & !< This is the first time step of the month
      & new_year                  !< This is the first time step of the year

    REAL(wp), INTENT(in) :: &
      & delta_time,             & !< Time step length [s]
      & seconds_day,            & !< Length of a day [s]
      & seconds_month,          & !< Length of the previous month [s]
      & t_air                     !< Near surface air temperature [K]
    REAL(wp), INTENT(inout) :: &
      & temp_sum_day,           & !< Current temperature sum of the day [K]
      & temp_sum_month,         & !< Current temperature sum of the month [K]
      & min_mmtemp_of_yr,       & !< Current minimum monthly mean temperature of the year [K]
      & max_mmtemp_of_yr,       & !< Current maximum monthly mean temperature of the year [K]
      & min_mmtemp20,           & !< Climatological mean temperature of the coldest month [K]
      & max_mmtemp20,           & !< Climatological mean temperature of the warmest month [K]
      & gdd_sum_year,           & !< Current number of growing degree days of this year
      & gdd_prev_year             !< Growing degree days of previous year

    ! Local variables
    REAL(wp) :: &
      & prev_day_mean_temp,     & !< Air temperature of previous day [K]
      & prev_month_mean_temp      !< Air temperature of previous month [K]
    ! ------------------------------------------------------------------------------------------------------ !

    IF (new_year) THEN
      gdd_prev_year = gdd_sum_year
      gdd_sum_year = 0._wp
    END IF

    IF (.NOT. new_day .OR. lstart) THEN
       temp_sum_day = temp_sum_day + t_air * delta_time
    ELSE
       prev_day_mean_temp = temp_sum_day / seconds_day
       temp_sum_day = t_air * delta_time
       gdd_sum_year = gdd_sum_year + MAX(0._wp,prev_day_mean_temp - tmelt - 5._wp)
    ENDIF

    IF (.NOT. new_month .OR. lstart) THEN
       temp_sum_month = temp_sum_month + t_air * delta_time
    ELSE
       prev_month_mean_temp = temp_sum_month / seconds_month
       temp_sum_month = t_air * delta_time
       ! --- minimum/maximum monthly mean temperature (20 year climatology)
       ! --- calculate coldest/warmest month of the year so far
       min_mmtemp_of_yr = min(min_mmtemp_of_yr, prev_month_mean_temp)
       max_mmtemp_of_yr = max(max_mmtemp_of_yr, prev_month_mean_temp)
       IF (new_year) THEN
          ! --- build 20yr running mean of coldest and warmest monthly temperature
          IF (init_running_means) THEN
             min_mmtemp20 = min_mmtemp_of_yr
             max_mmtemp20 = max_mmtemp_of_yr
          ELSE
             min_mmtemp20 = (min_mmtemp20 * 19._wp + min_mmtemp_of_yr) / 20._wp
             max_mmtemp20 = (max_mmtemp20 * 19._wp + max_mmtemp_of_yr) / 20._wp
          ENDIF
          ! --- reset for current year
          min_mmtemp_of_yr =  1000._wp
          max_mmtemp_of_yr = -1000._wp
       END IF ! new_year
    END IF ! new_month

  END SUBROUTINE calc_climbuf

#endif
END MODULE mo_climbuf_process
