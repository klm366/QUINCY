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
MODULE mo_carbon_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: wp
  USE mo_jsb_varlist,            ONLY: VARNAME_LEN
  USE mo_jsb_physical_constants, ONLY: amco2  !< Molar weight of CO2 [g/mol]

  IMPLICIT NONE
  PUBLIC

  ! Yasso requires tile dependant lctlib constants for carbon relocation, accessed via indices
  ENUM, BIND(C)
    ENUMERATOR :: i_lctlib_acid=1, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble, i_lctlib_humus
  END ENUM

  ! Indices of the yasso(_out) transfer array
  ENUM, BIND(C)
    ENUMERATOR :: i_pool_acid_ag=1, i_pool_water_ag, i_pool_ethanol_ag, i_pool_nonsoluble_ag,             &
                & i_pool_acid_bg, i_pool_water_bg, i_pool_ethanol_bg, i_pool_nonsoluble_bg, i_pool_humus, &
                & i_flux_resp_total, i_flux_2_humus_ag, i_flux_2_humus_bg, i_flux_from_humus_diag,        &
                & i_flux_resp_litter_ag, i_flux_resp_litter_bg, i_rel_flux_AWEN, i_flux_from_humus,       &
                & i_flux_litter_2_humus
  END ENUM

  ! Further constants
  REAL(wp), PARAMETER :: molarMassC_kg           = 12.01_wp * 1.e-3_wp  !< 12.01e-3_wp   -  Mass of 1 mol C   in kg
  REAL(wp), PARAMETER :: molarMassCO2_kg         = amco2 * 1.e-3_wp     !< 44.0095e-3_wp -  Mass of 1 mol CO2 in kg

  REAL(wp), PARAMETER :: FireFracWood2Atmos      =     0.2_wp           !< Fraction of wood carbon emitted to the atm by fire.
                                                                        !< In JS3 called def_fire_fract_wood_2_atmos.
  REAL(wp), PARAMETER :: fract_wood_aboveGround  =     0.7_wp           !< Fraction of C above ground in wood pool (for separation
                                                                        !  of woody litter into above and below ground litter pools)
  REAL(wp), PARAMETER :: fract_green_aboveGround =     0.5_wp           !< Fraction of C above ground in green pool (for separation
                                                                        !  of green litter into above and below ground litter pools)
  REAL(wp), PARAMETER :: sec_per_day             = 86400._wp            !< seconds per day
  REAL(wp), PARAMETER :: days_per_year           =   365.25_wp
  REAL(wp), PARAMETER :: sec_per_year            = days_per_year * sec_per_day

  ! Parameters for nitrogen cycle, which are currently not used but may be used in the future
  REAL(wp), PARAMETER :: N2O_rate_nitrification = 0.001_wp        ! JSBACH assumes lower bound of estimation in the literature as:
                                                                  !   (0.1%) due to nitrification
  REAL(wp), PARAMETER :: N2O_rate_denitrification = 0.00125_wp    ! 0.0025_wp  ! & (0.2%) due to denitrification
                                                                  ! actually 0.2% of denitrified N, not of N inputs
                                                                  ! assumption: 50% loss of NO3 by leaching
  REAL(wp), PARAMETER :: N2O_rate               = 0.005_wp        ! JSBACH assumes 0.5-1%
                                                                  !    N2O emissions rate per day (<0.1%-0.2% nitrification );
                                                                  !    0.2-4.7% due to denit in Xu Ri, 2008 * therein
  REAL(wp), PARAMETER :: sminn_NH4_fraction     = 0.4_wp          ! soil mineral NH4 is 40% of SMINN pool and rest in the form of
                                                                  !  NO3. Ref: Xu Ri, 2008

  ! Definitions for lcc

  !> list of potential active variables for the carbon process: For these the carbon interface provides active relocation
  CHARACTER(len=VARNAME_LEN) :: carbon_potential_active_vars(11)  = [character(len=VARNAME_LEN) :: &
     & 'c_green', 'c_woods', 'c_reserve',                                                          &
     & 'c_acid_ag1', 'c_water_ag1', 'c_ethanol_ag1', 'c_nonsoluble_ag1',                           &
     & 'c_acid_ag2', 'c_water_ag2', 'c_ethanol_ag2', 'c_nonsoluble_ag2' ]

  !> list of required passive vars if vars are active (because these are used to relocate the active variables)
  INTEGER, PARAMETER :: nr_of_yasso_pools = 18

  CHARACTER(len=VARNAME_LEN) :: carbon_required_passive_vars(nr_of_yasso_pools) = [character(len=VARNAME_LEN) :: &
    & 'c_acid_ag1', 'c_water_ag1', 'c_ethanol_ag1', 'c_nonsoluble_ag1',               &
    & 'c_acid_bg1', 'c_water_bg1', 'c_ethanol_bg1', 'c_nonsoluble_bg1', 'c_humus_1',  &
    & 'c_acid_ag2', 'c_water_ag2', 'c_ethanol_ag2', 'c_nonsoluble_ag2',               &
    & 'c_acid_bg2', 'c_water_bg2', 'c_ethanol_bg2', 'c_nonsoluble_bg2', 'c_humus_2' ]

  !> helper for automatization (order needs to match var name list in carbon_required_passive_vars)
  INTEGER :: coefficient_ind_of_passive_vars(nr_of_yasso_pools) = (/ &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble,                 &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble, i_lctlib_humus, &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble,                 &
    & i_lctlib_acid, i_lctlib_water, i_lctlib_ethanol, i_lctlib_nonsoluble, i_lctlib_humus  /)

#endif
END MODULE mo_carbon_constants
