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

!  ! Yasso requires tile dependant lctlib constants for carbon relocation, accessed via indices
!  INTEGER, PARAMETER  :: i_lctlib_acid           = 1 !< Index of constants for acid yasso pools in the lctlib
!  INTEGER, PARAMETER  :: i_lctlib_water          = 2 !< Index of constants for water yasso pools in the lctlib
!  INTEGER, PARAMETER  :: i_lctlib_ethanol        = 3 !< Index of constants for ethanol yasso pools in the lctlib
!  INTEGER, PARAMETER  :: i_lctlib_nonsoluble     = 4 !< Index of constants for nonsoluble yasso pools in the lctlib
!  INTEGER, PARAMETER  :: i_lctlib_humus          = 5 !< Index of constants for humus yasso pools in the lctlib

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

  ! Yasso pool indicators for c pool in/out, rate constant in
  ENUM, BIND(C)
    ENUMERATOR :: i_y_acid=1, i_y_water, i_y_ethanol, i_y_nonsoluble, i_y_humus
  END ENUM

  ! Yasso flux indicators for flux out
  ENUM, BIND(C)
    ENUMERATOR :: i_f_resp=1, i_f_2_humus, i_f_f_humus, i_f_litresp
  END ENUM

  ! Yasso rate modifiers for litter/humus
  ENUM, BIND(C)
    ENUMERATOR :: i_m_litter=1, i_m_humus
  END ENUM

  ! Misc parameters in interface
  ENUM, BIND(C)
    ENUMERATOR :: i_WoodLitSize=1, i_RootDistBeta, i_FracInunMin
  END ENUM


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

  ! ------------------------------------------------------------------------------------------
  ! parameters for vertical YASSO
  !

  ! decomposition modifier parameters for temperature and depth modifier
  REAL(wp), PARAMETER :: Q10     = 1.5_wp         ! Q10 for deep soil C decomposition, based on CLM
!   REAL(wp), PARAMETER :: ref_dep = 0.5_wp         ! reference depth for depth dependence of decomposition, based on CLM 4.5
  REAL(wp), PARAMETER :: ref_dep = 10._wp         ! reference depth for depth dependence of decomposition, based on CLM 5.0

  ! decomposition modifiers for oxic / anoxic soils
  REAL(wp), PARAMETER :: d_anox = 0.35_wp         ! decomposition modifier for anaerobic decomposition
  REAL(wp), PARAMETER :: d_oxic = 1.0_wp          ! decomposition modifier for oxic decomposition

  ! Yasso parameters (Tuomi '11)
  ! parameters for temperature dependence
  REAL(wp), PARAMETER   :: temp_p1                    = 0.095_wp    ! first  order temperature dependence [1/C]
  REAL(wp), PARAMETER   :: temp_p2                    = -0.0014_wp  ! second order temperature dependence [1/C]
  ! parameter for precipitation dependence
  REAL(wp), PARAMETER   :: precip_p1                  = -1.21_wp    ! first  order precipitation dependence [a/m]
  ! parameters for size dependence
  REAL(wp), PARAMETER   :: size_p1                    = -1.71_wp    ! first  order size dependence [1/cm]
  REAL(wp), PARAMETER   :: size_p2                    = 0.86_wp     ! second order size dependence [1/cm]
  REAL(wp), PARAMETER   :: size_p3                    = -0.306_wp   ! size dependence power [ ]
  REAL(wp), PARAMETER   :: d_size_leaf                = 1.0_wp      ! leaf litter is not size dependent

  ! parameters for moisture dependence, DAYCENT formulation
  REAL(wp), PARAMETER   :: moist_p_a                 = 0.575_wp
  REAL(wp), PARAMETER   :: moist_p_b                 = 1.5_wp
  REAL(wp), PARAMETER   :: moist_p_d                 = 3.0_wp

  ! ------------------------------------------------------------------------------------------
  ! Ratio of C mass to biomass (Spitfire and emission of chemical species by fire)
  REAL(wp),  PARAMETER :: C_2_biomass_ratio = 0.5_wp

 !< organic matter [litter/duff?] bulk densities for below above and ground [kg/m3], see Chojnacky (2009), 0'Donnell (2009), (https://philbusey.com/sports_turf/organic_matter.html#table), Ahn (2009)
  REAL(wp), PARAMETER ::  rho_org         = 100._wp
  REAL(wp), PARAMETER ::  rho_org_ag      =  75._wp
  REAL(wp), PARAMETER ::  rho_org_bg      = 150._wp

#endif
END MODULE mo_carbon_constants
