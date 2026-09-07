!> Various implementation constants used in ICON-Land
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
!>#### Declaration and definition of various constants used in ICON-Land
!>
MODULE mo_jsb_impl_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: jp

  IMPLICIT NONE
  PUBLIC

  ! general parameters
  INTEGER, PARAMETER :: SHORT_NAME_LEN = 10
  INTEGER, PARAMETER :: MEDIUM_NAME_LEN = 20

  ! water balance check descriptors
  INTEGER, PARAMETER :: WB_IGNORE  = 0
  INTEGER, PARAMETER :: WB_LOGGING = 1
  INTEGER, PARAMETER :: WB_ERROR   = 2

  ! using the functionality of LOGICAL var with REAL var
  REAL(jp), PARAMETER :: false            = 0.0_jp    !< zero
  REAL(jp), PARAMETER :: true             = 1.0_jp    !< one
  REAL(jp), PARAMETER :: test_false_true  = 0.5_jp    !< 0.5, inbetween false and true; 'IF(var < test_false_true) -> FALSE'

  INTEGER,  PARAMETER :: ifs_nsoil         = 4                                  !< Number of soil levels in IFS initial file
  REAL(jp), PARAMETER :: ifs_soil_depth(4) = (/ 0.07_jp,0.28_jp,1._jp,2.89_jp/) !< Soil layer depths in IFS initial file [m]

  REAL(jp), PARAMETER :: def_parameters = -9999.0_jp  !< the default value of model parameters (not constants) prior to their init

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_impl_constants'

! CONTAINS


#endif
END MODULE mo_jsb_impl_constants
