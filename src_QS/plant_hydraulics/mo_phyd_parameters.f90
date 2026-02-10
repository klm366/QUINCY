!> QUINCY plant hydraulic parameters
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
!>#### declare and init phyd parameters
!>
MODULE mo_phyd_parameters
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: def_parameters

  IMPLICIT NONE
  PUBLIC

  REAL(wp), SAVE :: &
    &         hydraulic_solver_precision = def_parameters

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_phyd_parameters'

CONTAINS

  ! ======================================================================================================= !
  !> initialize parameters for the process: plant hydraulics
  !>
  !>   routine is called in mo_jsb_base
  !>
  SUBROUTINE init_phyd_parameters

    CHARACTER(len=*), PARAMETER :: routine = modname//':init_phyd_parameters'

    hydraulic_solver_precision       = 1e-06_wp
  END SUBROUTINE init_phyd_parameters

#endif
END MODULE mo_phyd_parameters
