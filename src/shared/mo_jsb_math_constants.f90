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

  USE mo_kind,                  ONLY: lp
  USE mo_math_constants,        ONLY: pi, deg2rad

  IMPLICIT NONE
  PUBLIC


  ! QUINCY simplified time control
  REAL(lp), PARAMETER   :: one_year  = 365._lp                !< length of one year in days
  REAL(lp), PARAMETER   :: one_month = 30._lp                 !< length of one month in days (rounded value of '365/12')
  REAL(lp), PARAMETER   :: one_day   = 86400._lp              !< length of one day in seconds (identical with 'seconds_per_day')
  ! conversion of degree to radians
  REAL(lp), PARAMETER   :: degcircle = 360._lp                !< degrees per circle

  REAL(lp), PARAMETER   :: hourday         = 24.0_lp   !< hours per day
  REAL(lp), PARAMETER   :: seconds_per_day = 86400._lp !< seconds per day

  ! QUINCY math helpers
  REAL(lp), PARAMETER   :: eps12 = 0.000000000001_lp
  REAL(lp), PARAMETER   :: eps8  = 0.00000001_lp
  REAL(lp), PARAMETER   :: eps4  = 0.0001_lp
  REAL(lp), PARAMETER   :: eps1  = 0.1_lp
  REAL(lp), PARAMETER   :: zero  = 0.0_lp           !< just zero


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_math_constants'

!CONTAINS

#endif
END MODULE mo_jsb_math_constants
