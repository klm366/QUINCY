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
!>### The CLIMBUF process provides long term mean climatic variables from atmospheric forcing.
!>
!> This module contains the climate buffer process routines. Core routine is [[calc_climbuf]], where
!> medium and long term climate variables are calculated from atmospheric forcing variables.
!>

MODULE mo_climbuf_process
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: lp
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_climbuf_2d, calc_climbuf_3d, get_relative_humidity_air

  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_process'
CONTAINS
  !
  ! ======================================================================================================== !
  !>
  !>#### Compute running averages of 2D climate variables
  !>
  !> We here calculate daily, annual or climatological mean, maximum or minimum values from
  !> atmospheric forcing data. This routine is used for 2D (surface) variables.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE SUBROUTINE calc_climbuf_2d ( &
    & lstart, init_running_means, new_day, new_month, new_year, &
    & dtime, t_air, rain, snow, wind_10m, q_rel_air, snow_depth, thaw_depth, &
    & seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
    & min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp_20yr, max_mmtemp_20yr, &
    & temp_15d, temp_15d_yDay, precip_15d, precip_15d_yDay, &
    & gdd5_sum_year, gdd5_prev_year, &
    & max_snow_of_yr, max_snow_prev_year, &
    & max_thaw_of_yr, max_thaw_prev_year, &
    & max_wind_of_day, max_wind_yDay, max_wind_4yr, &
    & q_rel_air_14d, q_rel_air_14d_yDay)

    LOGICAL, INTENT(in) :: &
      & lstart,                 & !< This is the first time step in a simulation
      & init_running_means,     & !< Initialize the min/max temperature climatology
      & new_day,                & !< This is the first time step of the day
      & new_month,              & !< This is the first time step of the month
      & new_year                  !< This is the first time step of the year

    REAL(lp), INTENT(in) :: &
      & dtime,                  & !< Time step length [s]
      & seconds_day,            & !< Length of a day [s]
      & seconds_month,          & !< Length of the previous month [s]
      & t_air,                  & !< Near surface air temperature [K]
      & rain,                   & !< Liquid precipitation [kg m-2 s-1]
      & snow,                   & !< Solid precipitation [kg m-2 s-1]
      & wind_10m,               & !< 10m wind speed [m/s]
      & q_rel_air,              & !< Relative humidity [%]
      & snow_depth,             & !< Snow depth [m]
      & thaw_depth                !< Thaw depth within soil column (== active layer) [m]

    REAL(lp), INTENT(inout) :: &
      & temp_sum_day,           & !< Current temperature sum of the day [K]
      & temp_sum_month,         & !< Current temperature sum of the month [K]
      & min_mmtemp_of_yr,       & !< Current minimum monthly mean temperature of the year [K]
      & max_mmtemp_of_yr,       & !< Current maximum monthly mean temperature of the year [K]
      & min_mmtemp_20yr,        & !< Climatological mean temperature of the coldest month [K]
      & max_mmtemp_20yr,        & !< Climatological mean temperature of the warmest month [K]
      & temp_15d,               & !< 15-day smoothed temperature [K]
      & temp_15d_yDay,          & !< 15-day smoothed temperature at the end of the previous day [K]
      & precip_15d,             & !< 15-day smoothed precipitation [kg m-2 s-1]
      & precip_15d_yDay,        & !< 15-day smoothed precipitation at the end of the previous day [kg m-2 s-1]
      & gdd5_sum_year,          & !< Current number of growing degree days of this year
      & gdd5_prev_year,         & !< Growing degree days of previous year
      & max_snow_of_yr,         & !< Maximum snow depth in current year [m]
      & max_snow_prev_year,     & !< Maximum snow depth in previous year [m]
      & max_thaw_of_yr,         & !< Maximum soil column thaw depth in current year [m]
      & max_thaw_prev_year,     & !< Maximum soil column thaw depth in previous year [m]
      & max_wind_of_day,        & !< Current maximum wind speed of the day [m s-1]
      & max_wind_yDay,          & !< Maximum wind speed of the previous day [m s-1]
      & max_wind_4yr,           & !< Smoothed maximum wind speed (~4 years) [m s-1]
      & q_rel_air_14d,          & !< Relative humidity smoothed in time (~14 days) [%]
      & q_rel_air_14d_yDay        !< Previous day value of smoothed relative humidity [%]

    ! Local variables
    REAL(lp) :: &
      & prev_day_mean_temp,     & !< Air temperature of previous day [K]
      & prev_month_mean_temp,   & !< Air temperature of previous month [K]
      & persist,                & !< Persistence factor adapted to model time step
      & F_15d,                  & !< Factor for 15-day smoothing
      & N_15d                     !< Factor for 15-day smoothing

    REAL(lp), PARAMETER :: persist_wind     = 0.9995_lp   ! ~4 years
    REAL(lp), PARAMETER :: persist_rel_hum  = 0.95_lp     ! ~14 days
    ! ------------------------------------------------------------------------------------------------------ !

    IF (.NOT. new_day .OR. lstart) THEN
      temp_sum_day = temp_sum_day + t_air * dtime
    ELSE
      prev_day_mean_temp = temp_sum_day / seconds_day
      temp_sum_day = t_air * dtime
    END IF

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
      gdd5_sum_year = 0._lp
    END IF
    IF (new_day .AND. .NOT. lstart) THEN
      gdd5_sum_year = gdd5_sum_year + MAX(0._lp, prev_day_mean_temp - tmelt - 5._lp)
    END IF

    !-- Snow depth
    !  ---------------------
    IF (.NOT. new_year .OR. lstart) THEN
       max_snow_of_yr = MAX(max_snow_of_yr, snow_depth)
    ELSE
       max_snow_prev_year = max_snow_of_yr
       max_snow_of_yr = 0._wp
       max_snow_of_yr = MAX(max_snow_of_yr, snow_depth)
    END IF

    !-- Soil thaw depth
    !  ---------------------
    IF (.NOT. new_year .OR. lstart) THEN
      max_thaw_of_yr = MAX(max_thaw_of_yr, thaw_depth)
    ELSE
      max_thaw_prev_year = max_thaw_of_yr
      max_thaw_of_yr = 0._wp
      max_thaw_of_yr = MAX(max_thaw_of_yr, thaw_depth)
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
          min_mmtemp_20yr = min_mmtemp_of_yr
          max_mmtemp_20yr = max_mmtemp_of_yr
        ELSE
          min_mmtemp_20yr = (min_mmtemp_20yr * 19._lp + min_mmtemp_of_yr) / 20._lp
          max_mmtemp_20yr = (max_mmtemp_20yr * 19._lp + max_mmtemp_of_yr) / 20._lp
        END IF
        ! --- reset for current year
        min_mmtemp_of_yr =  1000._lp
        max_mmtemp_of_yr = -1000._lp
      END IF ! new_year
    END IF

    !-- 15-day smoothing of temperature and precipitation (needed with yasso)
    !  ---------------------
    F_15d = EXP(-dtime / (15._lp * seconds_day))
    N_15d = 1._lp - F_15d

    temp_15d   = t_air         * N_15d  +  F_15d * temp_15d
    precip_15d = (rain + snow) * N_15d  +  F_15d * precip_15d
    IF (new_day .AND. .NOT. lstart) THEN
      ! Save previous day variables (for yasso)
      temp_15d_yDay   = temp_15d
      precip_15d_yDay = precip_15d
    END IF

    !-- 10m wind
    !  ---------------------
    IF (.NOT. new_day .OR. lstart) THEN
      IF (wind_10m > max_wind_of_day) THEN
        max_wind_of_day = wind_10m
      END IF
    ELSE
      max_wind_yDay = max_wind_of_day
      ! Maximum wind speed smoothed in time
      max_wind_4yr = max_wind_4yr * persist_wind + max_wind_of_day * (1._lp - persist_wind)
      max_wind_of_day = wind_10m
    END IF

    !-- Relative humidity
    !  ---------------------
    persist = persist_rel_hum ** (dtime/seconds_day)
    q_rel_air_14d = q_rel_air_14d * persist + q_rel_air * (1._lp - persist)
    ! Save previous day variables
    IF (new_day .AND. .NOT. lstart) THEN
      q_rel_air_14d_yDay = q_rel_air_14d
    END IF

  END SUBROUTINE calc_climbuf_2d
  !
  ! ======================================================================================================== !
  !>
  !>#### Compute running averages of 3D climate variables
  !>
  !> We here calculate daily, annual or climatological mean, maximum or minimum values from
  !> 3D (e.g. layered) data.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE SUBROUTINE calc_climbuf_3d ( &
    & lstart, new_day, new_year, dtime, t_soil, wtr_soil_rel, seconds_day, &
    & t_soil_15d, t_soil_15d_yDay, wtr_soil_rel_15d, wtr_soil_rel_15d_yDay, &
    & max_t_soil_year, permafrost_period)
    USE mo_jsb_physical_constants, ONLY: tmelt

    LOGICAL, INTENT(in) :: &
      & lstart,                 & !< This is the first time step in a simulation
      & new_day,                & !< This is the first time step of the day
      & new_year                  !< This is the first time step of the year

    REAL(wp), INTENT(in) :: &
      & dtime,                  & !< Time step length [s]
      & seconds_day,            & !< Length of a day [s]
      & t_soil,                 & !< Soil layer temperature [K]
      & wtr_soil_rel              !< Relative soil water content on layers [-]

    REAL(wp), INTENT(inout) ::  &
      & t_soil_15d,             & !< 15-day smoothed soil layer temperature [K]
      & t_soil_15d_yDay,        & !< 15-day smoothed soil layer temperature at the end of the previous day [K]
      & wtr_soil_rel_15d,       & !< 15-day smoothed relative soil layer moisture [-]
      & wtr_soil_rel_15d_yDay,  & !< 15-day smoothed relative soil layer moisture at the end of the previous day [-]
      & max_t_soil_year,        & !< Maximum soil layer temperature in previous year [K]
      & permafrost_period         !< Number of years of consecutively frozen soil layer [yr]

    ! Local variables
    REAL(wp) :: &
      & F_15d,                  & !< Factor for 15-day smoothing
      & N_15d                     !< Factor for 15-day smoothing

    ! ------------------------------------------------------------------------------------------------------ !

    !-- 15-day smoothing of soil layer temperature and moisture (needed with vertical yasso)
    !  ---------------------
    F_15d = EXP(-dtime / (15._wp * seconds_day))
    N_15d = 1._wp - F_15d

    t_soil_15d = t_soil * N_15d  +  F_15d * t_soil_15d
    wtr_soil_rel_15d = wtr_soil_rel * N_15d  +  F_15d * wtr_soil_rel_15d
    ! TODO: Check, if these yDay variables are actually needed
    IF (new_day .AND. .NOT. lstart) THEN
      ! Save previous day variables (for yasso)
      t_soil_15d_yDay   = t_soil_15d
      wtr_soil_rel_15d_yDay = wtr_soil_rel_15d
    END IF

    !-- Soil thermal states (temperature and permafrost)
    !  ---------------------
    IF (new_year) THEN
      IF (max_t_soil_year < tmelt) THEN
        permafrost_period = permafrost_period + 1._wp
      ELSE
        permafrost_period = 0._wp
      END IF
      max_t_soil_year = t_soil
    ELSE
      max_t_soil_year = MAX(t_soil, max_t_soil_year)
    END IF

  END SUBROUTINE calc_climbuf_3d

  ! ======================================================================================================== !
  !>
  !>#### Compute relative humidity
  !>
  !> We here calculate relative humidity of the lowest atmospheric level
  !>
  ! -------------------------------------------------------------------------------------------------------- !
#ifndef _OPENACC
  ELEMENTAL PURE &
#endif
  FUNCTION get_relative_humidity_air(air_moisture, air_temperature, air_press, use_tmx) RESULT(relative_humidity_air)

    USE mo_jsb_physical_constants, ONLY: rd  ! gas constant for dry air [J/(K*kg)]
    USE mo_jsb_physical_constants, ONLY: rv  ! gas constant for water vapor [J/(K*kg)]
    USE mo_phy_schemes,            ONLY: qsat_water, qsat_mixed

    !$ACC ROUTINE SEQ
    REAL(lp), INTENT(in) :: &
      & air_moisture,       & !< Specific humidity of air [kg kg-1]
      & air_temperature,    & !< Temperature of air [K]
      & air_press             !< Pressure of air [Pa]
    LOGICAL, INTENT(in) :: &
      & use_tmx               !< Turbulent mixing scheme 'tmx' is used in coupling to the atmosphere:
                              !< if true, we use qsat_water, else qsat_mixed

    REAL(lp) :: &
      & air_qsat,           & !< Saturation specific humidity []
      & relative_humidity_air !< Relative humidity [%]

    IF (use_tmx) THEN
      air_qsat = qsat_water(air_temperature, air_press)
    ELSE
      air_qsat = qsat_mixed(air_temperature, air_press)
    END IF
    relative_humidity_air = 100._lp * air_moisture                                           &
                           * (  (1._lp - rd / rv) * air_qsat + (rd / rv) )                   &
                           / ( ((1._lp - rd / rv) * air_moisture + (rd / rv)) * air_qsat )
    relative_humidity_air = MIN(100._lp, relative_humidity_air)

  END FUNCTION get_relative_humidity_air
#endif
END MODULE mo_climbuf_process
