!> Contains mathematical constants
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
MODULE mo_jsb_math_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: jp
  USE mo_math_constants,        ONLY: pi, deg2rad

  IMPLICIT NONE
  PUBLIC


  ! QUINCY simplified time control
  REAL(jp), PARAMETER   :: one_year  = 365._jp                !< length of one year in days
  REAL(jp), PARAMETER   :: one_month = 30._jp                 !< length of one month in days (rounded value of '365/12')
  REAL(jp), PARAMETER   :: one_day   = 86400._jp              !< length of one day in seconds (identical with 'seconds_per_day')
  ! conversion of degree to radians
  REAL(jp), PARAMETER   :: degcircle = 360._jp                !< degrees per circle

  REAL(jp), PARAMETER   :: hourday         = 24.0_jp   !< hours per day
  REAL(jp), PARAMETER   :: seconds_per_day = 86400._jp !< seconds per day

  ! QUINCY math helpers
  REAL(jp), PARAMETER   :: eps12 = 0.000000000001_jp
  REAL(jp), PARAMETER   :: eps8  = 0.00000001_jp
  REAL(jp), PARAMETER   :: eps4  = 0.0001_jp
  REAL(jp), PARAMETER   :: eps1  = 0.1_jp
  REAL(jp), PARAMETER   :: zero  = 0.0_jp           !< just zero


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_math_constants'

!CONTAINS

#endif
END MODULE mo_jsb_math_constants
