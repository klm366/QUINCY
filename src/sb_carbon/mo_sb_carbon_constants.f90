!> Contains constants for the carbon processes
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_sb_carbon_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_jsb_varlist,            ONLY: VARNAME_LEN
  USE mo_jsb_physical_constants, ONLY: amco2  !< Molar weight of CO2 [g/mol]

  IMPLICIT NONE
  PUBLIC

  ! Yasso requires tile dependant lctlib constants for carbon relocation, accessed via indices

  ! Further constants
  REAL(wp), PARAMETER :: molarMassC_kg           = 12.01_wp * 1.e-3_wp  !< 12.01e-3_wp   -  Mass of 1 mol C   in kg
  !$ACC DECLARE COPYIN(molarMassC_kg)
  REAL(wp), PARAMETER :: molarMassCO2_kg         = amco2 * 1.e-3_wp     !< 44.0095e-3_wp -  Mass of 1 mol CO2 in kg
  !$ACC DECLARE COPYIN(molarMassCO2_kg)

  REAL(wp), PARAMETER :: FireFracWood2Atmos      =  0.2_wp              !< Fraction of wood carbon emitted to the atm by fire.
                                                                        !< In JS3 called def_fire_fract_wood_2_atmos.
  !$ACC DECLARE COPYIN(FireFracWood2Atmos)
  REAL(wp), PARAMETER :: fract_wood_aboveGround  = 0.7_wp               !< Fraction of C above ground in wood pool (for separation
                                                                        !  of woody litter into above and below ground litter pools)
  !$ACC DECLARE COPYIN(fract_wood_aboveGround)
  REAL(wp), PARAMETER :: fract_green_aboveGround = 0.5_wp               !< Fraction of C above ground in green pool (for separation
                                                                        !  of green litter into above and below ground litter pools)
  !$ACC DECLARE COPYIN(fract_green_aboveGround)
  REAL(wp), PARAMETER :: sec_per_day             = 86400._wp            !< seconds per day
  !$ACC DECLARE COPYIN(sec_per_day)
  REAL(wp), PARAMETER :: days_per_year           = 365.25_wp
  !$ACC DECLARE COPYIN(days_per_year)
  REAL(wp), PARAMETER :: sec_per_year            = days_per_year * sec_per_day
  !$ACC DECLARE COPYIN(sec_per_year)
  ! Definitions for lcc

#endif
END MODULE mo_sb_carbon_constants
