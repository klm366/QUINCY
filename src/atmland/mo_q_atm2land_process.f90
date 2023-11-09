
!> Summary: quincy atm2land routines
!>
!>### Contains routines for the a2l (atm2land) proc
!>
!>
!> Copyright (c) 2023 Julia Nabel
!>
!>#### Revision History
!> - First version    by Julia Nabel (2023-11-06)
!>
!>#### License
!>
!> QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!> Copyright (C) 2023 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany
!>
!> This file is part of the QUINCY model.
!> The QUINCY model is free software: you can redistribute it and/or modify
!> it under the terms of the GNU General Public License version 3,
!> as published by the Free Software Foundation.
!>
!> The QUINCY model is distributed in the hope that it will be useful,
!> but WITHOUT ANY WARRANTY; without even the implied warranty of
!> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!> GNU General Public License for more details.
!>
!> You should have received a copy of the GNU General Public License
!> along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> ======================================================================================================
MODULE mo_q_atm2land_process
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_math_constants,    ONLY: eps8

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_local_time_and_daytime_counter

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_atm2land_process'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> updates local time and updates or resets quincys daytime counter depending on the local time
  !>
  PURE ELEMENTAL SUBROUTINE update_local_time_and_daytime_counter( &
      &                       global_seconds_day, lon, swpar_srf_down, daytime_counter, local_time_day_seconds)

    USE mo_run_config,         ONLY: dtime
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER ,  INTENT(in)        :: global_seconds_day     !< seconds passed on this day (global time)
    REAL(wp),  INTENT(in)        :: lon                    !< longitude
    REAL(wp),  INTENT(in)        :: swpar_srf_down         !< current shortwave surface downward flux
    REAL(wp),  INTENT(inout)     :: daytime_counter        !< nr of timesteps of current day with sunlight
    REAL(wp),  INTENT(inout)     :: local_time_day_seconds !< seconds passed on this day (local cell time)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: time_shift_in_hours
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_local_time_and_daytime_counter'
    ! ----------------------------------------------------------------------------------------------------- !

    ! reset the daytime_counter if last timestep was midnight
    IF (ABS(local_time_day_seconds - dtime) < eps8) THEN
      daytime_counter = 0.0_wp
    ENDIF

    IF (swpar_srf_down > eps8) THEN
      ! count number of daytime values
      daytime_counter = daytime_counter + 1.0_wp
    ENDIF

    ! calculate the local time shift in whole hours (NINT!) relative to the global time (15 degrees per hour)
    time_shift_in_hours = NINT(MODULO(lon,360._wp)/15._wp)
    ! convert the time shift from whole hours to whole seconds and calculate the seconds passed on this day in local time
    ! (MODULO 86400 is applied to detemine the seconds of the local day, as this can be the day before or after the global day)
    local_time_day_seconds  = REAL(MODULO(global_seconds_day + (time_shift_in_hours * 3600), 86400),wp)

  END SUBROUTINE update_local_time_and_daytime_counter

#endif
END MODULE mo_q_atm2land_process

