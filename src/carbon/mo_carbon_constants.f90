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
  !$ACC DECLARE COPYIN(molarMassC_kg)
  REAL(wp), PARAMETER :: molarMassCO2_kg         = amco2 * 1.e-3_wp     !< 44.0095e-3_wp -  Mass of 1 mol CO2 in kg
  !$ACC DECLARE COPYIN(molarMassCO2_kg)
  REAL(wp), PARAMETER :: FireFracWood2Atmos      =     0.2_wp           !< Fraction of wood carbon emitted to the atm by fire.
                                                                        !< In JS3 called def_fire_fract_wood_2_atmos.
  !$ACC DECLARE COPYIN(FireFracWood2Atmos)
  REAL(wp), PARAMETER :: fract_wood_aboveGround  =     0.7_wp           !< Fraction of C above ground in wood pool
  !$ACC DECLARE COPYIN(fract_wood_aboveGround)
  REAL(wp), PARAMETER :: fract_green_aboveGround =     0.5_wp           !< Fraction of C above ground in green pool
  !$ACC DECLARE COPYIN(fract_green_aboveGround)
  REAL(wp), PARAMETER :: sec_per_day             = 86400._wp            !< seconds per day
  !$ACC DECLARE COPYIN(sec_per_day)
  REAL(wp), PARAMETER :: days_per_year           =   365.25_wp
  !$ACC DECLARE COPYIN(days_per_year)
  REAL(wp), PARAMETER :: sec_per_year            = days_per_year * sec_per_day
  !$ACC DECLARE COPYIN(sec_per_year)

  !   parameters for Yasso peat
  REAL(wp), parameter :: rho_YCpool_acro          = 3.5e4_wp           ! Bulk density of acrotelm peat: 0.035 g/cm3, based on
                                                                       ! Granberg (1999): rho upper 20cm: 0.02 g/cm^3  */
                                                                       !  rho 20-40cm: 0.05 g/cm3, > 40 cm: 0.1-0.12 g/cm3
  REAL(wp), parameter :: rho_YCpool_cato          = 9.1e4_wp           ! Bulk density of catotelm peat: 0.091 g/cm3, based on
                                                                       ! Clymo (1983); Turunen (2002)
  REAL(wp), parameter :: Cfrac_YCpool_acro        = 0.50_wp           ! C fraction in acrotelm peat Malmer and Wallen (2004)
  REAL(wp), parameter :: Cfrac_YCpool_cato        = 0.52_wp           ! Malmer and Wallen (2004)

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
