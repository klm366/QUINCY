!> QUINCY fire-disturbance constants
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### declare and define fire-disturbance constants
!>
MODULE mo_dfire_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                 ONLY: wp

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_dfire_constants'

CONTAINS
  
  ! ======================================================================================================= !  
  !>
  !> Initialize constants used in the dfire process
  !>
  SUBROUTINE init_dfire_constants
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_dfire_constants'
    ! ----------------------------------------------------------------------------------------------------- !
    ! ...
  END SUBROUTINE init_dfire_constants
  
#endif
END MODULE mo_dfire_constants
