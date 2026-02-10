!> Time control and time helpers
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!>#### Provides routines related to time and time control for the quincy model
!>
MODULE mo_qs_time_util

  USE mo_kind,               ONLY: wp
  USE mo_util,               ONLY: int2string
  USE mo_jsb_math_constants, ONLY: eps8
  USE mo_qs_output_txt,      ONLY: lun_sim_runtime_estimate


  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: open_quincy_time_output_files, close_quincy_time_output_files
  PUBLIC :: calc_and_report_estimated_remaining_simulation_time

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_qs_time_util'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> open output files that are used in this module
  !>
  SUBROUTINE open_quincy_time_output_files()
    INTEGER                     :: lun   !< logical unit number
    CHARACTER(100)              :: fname !< filename
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':open_quincy_time_output_files'
    ! ----------------------------------------------------------------------------------------------------- !

    lun   = lun_sim_runtime_estimate
    fname = 'sim_runtime_estimate.txt'
    OPEN (UNIT=lun, FILE=TRIM(fname), ACTION="write", STATUS="unknown")
  END SUBROUTINE open_quincy_time_output_files
  
  ! ======================================================================================================= !
  !>
  !> open output files that are used in this module
  !>
  SUBROUTINE close_quincy_time_output_files()
    INTEGER                     :: lun   !< logical unit number
    LOGICAL                     :: file_open
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':close_quincy_time_output_files'
    ! ----------------------------------------------------------------------------------------------------- !

    ! close 'sim_runtime_estimate.txt'
    lun   = lun_sim_runtime_estimate
    INQUIRE(lun, opened = file_open)
    IF(file_open) CLOSE(lun) 
  END SUBROUTINE close_quincy_time_output_files

  ! ======================================================================================================= !
  !>
  !> calculate and report progress of simulation and estimate remaining time
  !>
  SUBROUTINE calc_and_report_estimated_remaining_simulation_time( &
    &         first_timestep, &
    &         itimestep, &
    &         ntimestep)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,      INTENT(in)  :: first_timestep
    INTEGER,      INTENT(in)  :: itimestep
    INTEGER,      INTENT(in)  :: ntimestep
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), SAVE              :: cpu_time_previous = 0.0_wp !< cpu runtime at previous timestep [s]
    REAL(wp)                    :: cpu_time_now               !< cpu runtime at current timestep [s]
    REAL(wp)                    :: frac_simulation_completed  !< fraction of completion of the estimated simulation time
    INTEGER                     :: progress_interval_time     !< interval between progress messages [s]
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_and_report_estimated_remaining_simulation_time'
    ! ----------------------------------------------------------------------------------------------------- !

    ! init local var
    progress_interval_time = 35

    ! get the time the CPU has spend for the simulation until now [s] (TYPE: REAL with decimal)
    CALL CPU_TIME(cpu_time_now)
    ! report new simulation estimations if a specific amount of time has passed since previous timestep
    IF (cpu_time_now - cpu_time_previous > progress_interval_time) THEN
      frac_simulation_completed = REAL(itimestep - first_timestep + 1, wp) / REAL(ntimestep - first_timestep + 1, wp)
      ! Avoid division by zero when fraction completed is very low
      IF (frac_simulation_completed > 0.0_wp) THEN
        WRITE(lun_sim_runtime_estimate, *) &
          & " Elapsed: ",  TRIM(convert_seconds_to_time_string(cpu_time_now)),  &
          & ", remaining: ", TRIM(convert_seconds_to_time_string((cpu_time_now / frac_simulation_completed) &
          & * (1.0_wp - frac_simulation_completed))), &
          & ", ", TRIM(int2string(NINT(frac_simulation_completed * 100.0_wp))), "% finished."
      ENDIF
      ! Save time of current timestep
      cpu_time_previous = cpu_time_now
    ENDIF
  END SUBROUTINE calc_and_report_estimated_remaining_simulation_time

  ! ======================================================================================================= !
  !>
  !> returns a date string hh:mm:ss based on real number value of seconds
  !>
  FUNCTION convert_seconds_to_time_string(sim_run_cputime) result(time_string)
    REAL(wp), INTENT(in) :: sim_run_cputime   !< time the CPU has spend for the simulation [s]
    CHARACTER(len=26)    :: time_string       !< hh:mm:ss
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER           :: loc_sim_run_cputime
    INTEGER           :: hlp1
    INTEGER           :: second
    INTEGER           :: minute
    INTEGER           :: hour
    CHARACTER(len=3)  :: sec_string
    CHARACTER(len=3)  :: min_string
    CHARACTER(len=20) :: hour_string
    ! ----------------------------------------------------------------------------------------------------- !

    ! Rounding to the nearest "full second", i.e., cutting milliseconds; and typecasting REAL -> INTEGER
    loc_sim_run_cputime = NINT(sim_run_cputime)
    ! Calc time information from cpu runtime, which is provided in seconds
    hour    = loc_sim_run_cputime / 3600
    hlp1    = loc_sim_run_cputime - hour * 3600
    minute  = hlp1 / 60 
    second  = hlp1 - minute * 60
    ! Convert an integer eg. 4 to two string e.g '04' with a leading 0 if smaller then 2
    hour_string = convert_int_to_time_digit_string(hour)
    min_string  = convert_int_to_time_digit_string(minute)
    sec_string  = convert_int_to_time_digit_string(second)
    ! Build time unit string
    time_string = TRIM(hour_string) // ":" // TRIM(min_string) // ":"// TRIM(sec_string)
  END FUNCTION convert_seconds_to_time_string

  ! ======================================================================================================= !
  !>
  !> returns a character converted from an integer with a leading zero if only one digit, e.g., '3 -> 03'
  !>
  FUNCTION convert_int_to_time_digit_string(time) result(time_digit_string)
    INTEGER, INTENT(in) :: time                 !< time value as integer
    CHARACTER(len=12)   :: time_digit_string    !< time value as character
    ! ----------------------------------------------------------------------------------------------------- !

    IF(time < 10) THEN
      time_digit_string  = "0"// TRIM(int2string(time))
    ELSE
      time_digit_string  = TRIM(int2string(time))
    END IF
  END FUNCTION convert_int_to_time_digit_string

END MODULE mo_qs_time_util
