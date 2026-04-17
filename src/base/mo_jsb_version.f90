!> Contains methods to define and print JSBACH version number and label a simulation.
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

#ifdef __xlC__
@PROCESS STRICT
#endif

MODULE mo_jsb_version
#ifndef __NO_JSBACH__

  USE mo_exception,         ONLY: message, message_text

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: get_icon_land_version

  CHARACTER(len=*), PARAMETER :: icon_land_version = '2025.10-2' !< ICON-Land version string

CONTAINS

  !>
  !! Define JSBACH version and labels
  !!
  !! Define JSBACH version number and get some information about a simulation
  !! at model initialization.
  !!
  FUNCTION get_icon_land_version() RESULT(version)

    CHARACTER(len=:), ALLOCATABLE :: version

    version = icon_land_version

  END FUNCTION get_icon_land_version

#endif
END MODULE mo_jsb_version
