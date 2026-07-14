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

  USE mo_kind,        ONLY: wp
  USE mo_jsb_physical_constants, ONLY: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_climbuf, get_relative_humidity_air

  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_process'

  REAL(wp), PARAMETER :: pft13_vcmax25_intercept = 49.5690779208497_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_mat_coef  = -0.391997574087549_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_map_coef  = 0.001594603433354_wp
  REAL(wp), PARAMETER :: pft13_jmax25_intercept  = 180.411473978400_wp
  REAL(wp), PARAMETER :: pft13_jmax25_mat_coef   = -3.09078405961564_wp
  REAL(wp), PARAMETER :: pft13_jmax25_map_coef   = -0.00485175158785102_wp
  REAL(wp), PARAMETER :: pft13_sla_intercept     = 0.060960327_wp
  REAL(wp), PARAMETER :: pft13_sla_map_coef      = -0.000045172134_wp
  REAL(wp), PARAMETER :: pft13_sla_tmax_coef     = -0.0015723334_wp
  REAL(wp), PARAMETER :: pft13_sla_map_tmax_coef = 0.0000015332445_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_lower     = 16.5700_wp
  REAL(wp), PARAMETER :: pft13_vcmax25_upper     = 66.2500_wp
  REAL(wp), PARAMETER :: pft13_jmax25_lower      = 46.3775_wp
  REAL(wp), PARAMETER :: pft13_jmax25_upper      = 131.4725_wp
  REAL(wp), PARAMETER :: pft13_sla_lower         = 0.006295_wp
  REAL(wp), PARAMETER :: pft13_sla_upper         = 0.0323334_wp

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
  ELEMENTAL PURE SUBROUTINE calc_climbuf (lstart, init_running_means, new_day, new_month, new_year, is_pft13, &
    &                                     dtime, t_air, rain, snow, wind_10m, q_rel_air, &
    &                                     seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
    &                                     min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp_20yr, max_mmtemp_20yr, &
    &                                     temp_15d, temp_15d_yDay, precip_15d, precip_15d_yDay, &
    &                                     gdd5_sum_year, gdd5_prev_year, &
    &                                     max_wind_of_day, max_wind_yDay, max_wind_4yr, &
    &                                     q_rel_air_14d, q_rel_air_14d_yDay, &
    &                                     temp_sum_year, time_sum_year, precip_sum_year, &
    &                                     mat_prev_year, map_prev_year, max_mmtemp_prev_year, &
    &                                     pft13_vcmax25_raw, pft13_vcmax25_pred, &
    &                                     pft13_jmax25_raw, pft13_jmax25_pred, &
    &                                     pft13_sla_raw, pft13_sla_pred, pft13_jmax_vcmax_ratio)

    LOGICAL, INTENT(in) :: &
      & lstart,                 & !< This is the first time step in a simulation
      & init_running_means,     & !< Initialize the min/max temperature climatology
      & new_day,                & !< This is the first time step of the day
      & new_month,              & !< This is the first time step of the month
      & new_year,               & !< This is the first time step of the year
      & is_pft13                  !< True only for the explicit PFT13 child tile

    REAL(wp), INTENT(in) :: &
      & dtime,                  & !< Time step length [s]
      & seconds_day,            & !< Length of a day [s]
      & seconds_month,          & !< Length of the previous month [s]
      & t_air,                  & !< Near surface air temperature [K]
      & rain,                   & !< Liquid precipitation [kg m-2 s-1]
      & snow,                   & !< Solid precipitation [kg m-2 s-1]
      & wind_10m,               & !< 10m wind speed [m/s]
      & q_rel_air                 !< Relative humidity [%]

    REAL(wp), INTENT(inout) :: &
      & temp_sum_day,           & !< Current temperature sum of the day [K]
      & temp_sum_month,         & !< Current temperature sum of the month [K]
      & min_mmtemp_of_yr,       & !< Current minimum monthly mean temperature of the year [K]
      & max_mmtemp_of_yr,       & !< Current maximum monthly mean temperature of the year [K]
      & min_mmtemp_20yr,        & !< Climatological mean temperature of the coldest month [K]
      & max_mmtemp_20yr,        & !< Climatological mean temperature of the warmest month [K]
      & temp_15d,               & !< 15-day smoothed temperature [K]
      & temp_15d_yDay,          & !< 15-day smoothed temperature at the end of the previous day [K]
      & precip_15d,             & !< 15-day smoothed precipitation [K]
      & precip_15d_yDay,        & !< 15-day smoothed precipitation at the end of the previous day [K]
      & gdd5_sum_year,          & !< Current number of growing degree days of this year
      & gdd5_prev_year,         & !< Growing degree days of previous year
      & max_wind_of_day,        & !< Current maximum wind speed of the day [m s-1]
      & max_wind_yDay,          & !< Maximum wind speed of the previous day [m s-1]
      & max_wind_4yr,           & !< Smoothed maximum wind speed (~4 years) [m s-1]
      & q_rel_air_14d,          & !< Relative humidity smoothed in time (~14 days) [%]
      & q_rel_air_14d_yDay,     & !< Previous day value of smoothed relative humidity [%]
      & temp_sum_year,          & !< Current annual air-temperature integral [K s]
      & time_sum_year,          & !< Current accumulated time within the year [s]
      & precip_sum_year,        & !< Current annual precipitation accumulation [kg m-2]
      & mat_prev_year,          & !< Previous-year mean air temperature [K]
      & map_prev_year,          & !< Previous-year total precipitation [kg m-2]
      & max_mmtemp_prev_year,   & !< Previous-year warmest monthly mean temperature [K]
      & pft13_vcmax25_raw,      & !< Unbounded PFT13 Vcmax25 diagnostic [umol m-2 s-1]
      & pft13_vcmax25_pred,     & !< Bounded PFT13 Vcmax25 diagnostic [umol m-2 s-1]
      & pft13_jmax25_raw,       & !< Unbounded PFT13 Jmax25 diagnostic [umol m-2 s-1]
      & pft13_jmax25_pred,      & !< Bounded PFT13 Jmax25 diagnostic [umol m-2 s-1]
      & pft13_sla_raw,          & !< Unbounded PFT13 SLA diagnostic [m2 gC-1]
      & pft13_sla_pred,         & !< Bounded PFT13 SLA diagnostic [m2 gC-1]
      & pft13_jmax_vcmax_ratio    !< Unbounded PFT13 Jmax25 to Vcmax25 ratio [-]

    ! Local variables
    REAL(wp) :: &
      & prev_day_mean_temp,     & !< Air temperature of previous day [K]
      & prev_month_mean_temp,   & !< Air temperature of previous month [K]
      & persist,                & !< Persistence factor adapted to model time step
      & F_15d,                  & !< Factor for 15-day smoothing
      & N_15d,                  & !< Factor for 15-day smoothing
      & mat_prev_c,             & !< Previous-year mean air temperature [degC]
      & tmax_prev_c               !< Previous-year warmest monthly mean temperature [degC]

    REAL(wp), PARAMETER :: persist_wind     = 0.9995_wp   ! ~4 years
    REAL(wp), PARAMETER :: persist_rel_hum  = 0.95_wp     ! ~14 days
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
          min_mmtemp_20yr = min_mmtemp_of_yr
          max_mmtemp_20yr = max_mmtemp_of_yr
        ELSE
          min_mmtemp_20yr = (min_mmtemp_20yr * 19._wp + min_mmtemp_of_yr) / 20._wp
          max_mmtemp_20yr = (max_mmtemp_20yr * 19._wp + max_mmtemp_of_yr) / 20._wp
        ENDIF
      END IF ! new_year
    END IF

    !-- Previous-year climate predictors and PFT13 trait diagnostics
    !  ----------------------
    IF (new_year .AND. .NOT. lstart) THEN
      IF (time_sum_year > 0._wp) THEN
        mat_prev_year = temp_sum_year / time_sum_year
        map_prev_year = precip_sum_year
        max_mmtemp_prev_year = max_mmtemp_of_yr

        mat_prev_c  = mat_prev_year - tmelt
        tmax_prev_c = max_mmtemp_prev_year - tmelt

        ! Previous behaviour, retained in Git history:
        ! the PFT13-named diagnostics were evaluated on every CLIMBUF tile.
        ! The explicit tile guard below restricts their calculation to PFT13.
        IF (is_pft13) THEN
          pft13_vcmax25_raw = pft13_vcmax25_intercept &
            & + pft13_vcmax25_mat_coef * mat_prev_c &
            & + pft13_vcmax25_map_coef * map_prev_year

          pft13_jmax25_raw = pft13_jmax25_intercept &
            & + pft13_jmax25_mat_coef * mat_prev_c &
            & + pft13_jmax25_map_coef * map_prev_year

          pft13_sla_raw = pft13_sla_intercept &
            & + pft13_sla_map_coef * map_prev_year &
            & + pft13_sla_tmax_coef * tmax_prev_c &
            & + pft13_sla_map_tmax_coef * map_prev_year * tmax_prev_c

          pft13_vcmax25_pred = MIN(pft13_vcmax25_upper, &
            & MAX(pft13_vcmax25_lower, pft13_vcmax25_raw))
          pft13_jmax25_pred = MIN(pft13_jmax25_upper, &
            & MAX(pft13_jmax25_lower, pft13_jmax25_raw))
          pft13_sla_pred = MIN(pft13_sla_upper, &
            & MAX(pft13_sla_lower, pft13_sla_raw))

          IF (ABS(pft13_vcmax25_raw) > EPSILON(1._wp)) THEN
            pft13_jmax_vcmax_ratio = pft13_jmax25_raw / pft13_vcmax25_raw
          ELSE
            pft13_jmax_vcmax_ratio = 0._wp
          END IF
        ELSE
          pft13_vcmax25_raw = 0._wp
          pft13_vcmax25_pred = 0._wp
          pft13_jmax25_raw = 0._wp
          pft13_jmax25_pred = 0._wp
          pft13_sla_raw = 0._wp
          pft13_sla_pred = 0._wp
          pft13_jmax_vcmax_ratio = 0._wp
        END IF
      END IF
    END IF

    IF (new_year) THEN
      temp_sum_year = 0._wp
      time_sum_year = 0._wp
      precip_sum_year = 0._wp
      min_mmtemp_of_yr = 1000._wp
      max_mmtemp_of_yr = -1000._wp
    END IF

    ! Add the current time step after the new-year reset so that it belongs to the current year.
    temp_sum_year = temp_sum_year + t_air * dtime
    time_sum_year = time_sum_year + dtime
    precip_sum_year = precip_sum_year + (rain + snow) * dtime

    !-- 15-day smoothing of temperature and precipitation (needed with yasso)
    !  ---------------------
    F_15d = EXP(-dtime / (15._wp * seconds_day))
    N_15d = 1._wp - F_15d

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
      max_wind_4yr = max_wind_4yr * persist_wind + max_wind_of_day * (1._wp - persist_wind)
      max_wind_of_day = wind_10m
    END IF

    !-- Relative humidity
    !  ---------------------
    persist = persist_rel_hum ** (dtime/seconds_day)
    q_rel_air_14d = q_rel_air_14d * persist + q_rel_air * (1._wp - persist)
    ! Save previous day variables
    IF (new_day .AND. .NOT. lstart) THEN
      q_rel_air_14d_yDay = q_rel_air_14d
    END IF

  END SUBROUTINE calc_climbuf

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
    REAL(wp), INTENT(in) :: &
      & air_moisture,       & !< Specific humidity of air [kg kg-1]
      & air_temperature,    & !< Temperature of air [K]
      & air_press             !< Pressure of air [Pa]
    LOGICAL, INTENT(in) :: &
      & use_tmx               !< Turbulent mixing scheme 'tmx' is used in coupling to the atmosphere:
                              !< if true, we use qsat_water, else qsat_mixed

    REAL(wp) :: &
      & air_qsat,           & !< Saturation specific humidity []
      & relative_humidity_air !< Relative humidity [%]

    IF (use_tmx) THEN
      air_qsat = qsat_water(air_temperature, air_press)
    ELSE
      air_qsat = qsat_mixed(air_temperature, air_press)
    END IF
    relative_humidity_air = 100._wp * air_moisture                                           &
                           * (  (1._wp - rd / rv) * air_qsat + (rd / rv) )                   &
                           / ( ((1._wp - rd / rv) * air_moisture + (rd / rv)) * air_qsat )
    relative_humidity_air = MIN(100._wp, relative_humidity_air)

  END FUNCTION get_relative_humidity_air
#endif
END MODULE mo_climbuf_process
