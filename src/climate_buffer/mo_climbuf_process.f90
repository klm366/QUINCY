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
  USE mo_exception,   ONLY: message, finish, message_text
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_climbuf

  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_process'
CONTAINS
  !
  ! ===============================================================================================================================
  !>
  !> Compute running averages of climate variables for nlcc
  !>

  ELEMENTAL PURE SUBROUTINE calc_climbuf (lstart, init_running_means, new_day, new_month, new_year, delta_time, t_air, &
                                          seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
                                          min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20, &
                                          gdd_sum_year, gdd_prev_year)

    LOGICAL, INTENT(in) :: lstart ! first time step in a simulation
    LOGICAL, INTENT(in) :: init_running_means ! initialize the min/max temperature climatology
    LOGICAL, INTENT(in) :: new_day ! first time step in a day
    LOGICAL, INTENT(in) :: new_month ! firest time step in a month
    LOGICAL, INTENT(in) :: new_year ! first time step in a year
    REAL(wp), INTENT(in) :: delta_time ! time step length [s]
    REAL(wp), INTENT(in) :: seconds_day ! lenght of a day [s]
    REAL(wp), INTENT(in) :: seconds_month ! length of the previous month [s]
    REAL(wp), INTENT(in) :: t_air ! air temperature of lowest atmosphere level [K]
    REAL(wp), INTENT(inout) :: temp_sum_day
    REAL(wp), INTENT(inout) :: temp_sum_month
    REAL(wp), INTENT(inout) :: min_mmtemp_of_yr ! temperature of coldest month of previous year [K]
    REAL(wp), INTENT(inout) :: max_mmtemp_of_yr ! temperature of warmest month of previous year [K]
    REAL(wp), INTENT(inout) :: min_mmtemp20 ! temperature of coldest month climatology [K]
    REAL(wp), INTENT(inout) :: max_mmtemp20 ! temperature of warmest month climatology [K]
    REAL(wp), INTENT(inout) :: gdd_sum_year ! growing degree days of current year [C]
    REAL(wp), INTENT(inout) :: gdd_prev_year ! growing degree days of previous year [C]

!   local variables
    REAL(wp) :: prev_day_mean_temp   ! air temperature of previous day [K]
    REAL(wp) :: prev_month_mean_temp ! air temperature of previous month [K]

    IF (new_year) THEN
       gdd_prev_year = gdd_sum_year
       gdd_sum_year = 0._wp
    ENDIF

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
