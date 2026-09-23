! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Utility routines of the coastal-constraint process
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
!>#### Non-scientific helper routines of the coastal process
!>
MODULE mo_coast_util
#ifndef __NO_JSBACH__

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: days_since_2000

  CHARACTER(len=*), PARAMETER :: modname = 'mo_coast_util'

CONTAINS

  ! ======================================================================================================== !
  !>
  !>#### Number of days from 1 January 2000 to a date of the proleptic Gregorian calendar
  !>
  !> We here count days with the Julian day number of Fliegel and Van Flandern (1968),
  !> so that the tide can be evaluated at the model time. The tide is periodic, so the
  !> choice of 1 January 2000 as reference only shifts it in time. With a 365-day
  !> calendar the count would slip by one day at each leap day; the statistics of the
  !> tide height are not affected.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  ELEMENTAL PURE FUNCTION days_since_2000(year, month, day) RESULT(ndays)

    INTEGER, INTENT(in) :: &
      & year,              & !< Calendar year
      & month,             & !< Month of the year (1 to 12)
      & day                  !< Day of the month

    INTEGER :: ndays         !< RESULT: days since 1 January 2000 (negative before)

    INTEGER :: a, jdn

    a   = (month - 14) / 12
    jdn = (1461 * (year + 4800 + a)) / 4 + (367 * (month - 2 - 12 * a)) / 12 &
      &   - (3 * ((year + 4900 + a) / 100)) / 4 + day - 32075
    ndays = jdn - 2451545     ! 2451545 is the Julian day number of 1 January 2000

  END FUNCTION days_since_2000

#endif
END MODULE mo_coast_util
