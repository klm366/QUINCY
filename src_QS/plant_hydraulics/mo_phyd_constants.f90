!> QUINCY plant hydraulic constants
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
!>#### declare and define plant hydraulics constants
!>
MODULE mo_phyd_constants
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! parameters for changes in canopy N fractions
  REAL(wp), SAVE ::    &
  kg_h2o_to_mol,       &     !< Kilgoram water mol equivalent
  mpa_to_pa,           &     !< MPa to Pa conversion factor
  pa_to_mpa                  !< Pa to MPA conversion factor    

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_phyd_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from phyd 
  !>
  SUBROUTINE init_phyd_constants

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_pyhdro_constants'

    kg_h2o_to_mol     = 1000.0_wp / 18.0_wp         ! Universal chemical constant
    mpa_to_pa         = 1000000.0_wp                ! Mathematical constant
    pa_to_mpa         = 1.0_wp / mpa_to_pa          ! Mathematical constant

  END SUBROUTINE init_phyd_constants

#endif
END MODULE mo_phyd_constants
