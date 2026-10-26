!> Contains utility functions for the soil and snow energy process.
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
!>#### Provide utility functions for the soil and snow energy process
!>
!> The module is intended to host small auxiliary routines that support the SSE
!> scheme and are shared by several process-related modules.
!>
MODULE mo_sse_util
#ifndef __NO_JSBACH__

  !USE mo_kind,      ONLY: lp
  !USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_util'

  !> Utility helpers for the SSE process can be added here as needed.
CONTAINS

#endif
END MODULE mo_sse_util
