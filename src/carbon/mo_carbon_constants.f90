!> Contains constants for the carbon processes
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
MODULE mo_carbon_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                   ONLY: lp
  USE mo_jsb_varlist,            ONLY: VARNAME_LEN
  USE mo_jsb_physical_constants, ONLY: amco2  !< Molar weight of CO2 [g/mol]

  IMPLICIT NONE
  PUBLIC

  ! Yasso requires tile dependant lctlib constants for carbon relocation, accessed via indices
  ENUM, BIND(C)
    ENUMERATOR :: LCT_ACID=1, LCT_WATER, LCT_ETHANOL, LCT_NONSOLUBLE, LCT_HUMUS
  END ENUM

  ! Indices of the yasso(_out) transfer array
  ENUM, BIND(C)
    ENUMERATOR :: POOL_ACID_AG=1, POOL_WATER_AG, POOL_ETHANOL_AG, POOL_NONSOLUBLE_AG,       &
          & POOL_ACID_BG, POOL_WATER_BG, POOL_ETHANOL_BG, POOL_NONSOLUBLE_BG, POOL_HUMUS,   &
          & FLUX_RESP_TOTAL, FLUX_2_HUMUS_AG, FLUX_2_HUMUS_BG, FLUX_FROM_HUMUS_DIAG,        &
          & FLUX_RESP_LITTER_AG, FLUX_RESP_LITTER_BG, REL_FLUX_AWEN, FLUX_FROM_HUMUS,       &
          & FLUX_LITTER_2_HUMUS
  END ENUM

  ! Yasso pool indicators for c pool in/out, rate constant in
  ENUM, BIND(C)
    ENUMERATOR :: Y_ACID=1, Y_WATER, Y_ETHANOL, Y_NONSOLUBLE, Y_HUMUS
  END ENUM

  ! Yasso flux indicators for flux out
  ENUM, BIND(C)
    ENUMERATOR :: F_RESP=1, F_2_HUMUS, F_F_HUMUS, F_LITRESP
  END ENUM

  ! Yasso rate modifiers for litter/humus
  ENUM, BIND(C)
    ENUMERATOR :: M_LITTER=1, M_HUMUS
  END ENUM

  ! Misc parameters in interface
  ENUM, BIND(C)
    ENUMERATOR :: WOODLITSIZE=1, ROOT_DIST_BETA, FRAC_INUN_MIN
  END ENUM


  ! Further constants
  REAL(lp), PARAMETER :: molarMassC_kg           = 12.01_lp * 1.e-3_lp  !< 12.01e-3_lp   -  Mass of 1 mol C   in kg
  !$ACC DECLARE COPYIN(molarMassC_kg)
  REAL(lp), PARAMETER :: molarMassCO2_kg         = amco2 * 1.e-3_lp     !< 44.0095e-3_lp -  Mass of 1 mol CO2 in kg
  !$ACC DECLARE COPYIN(molarMassCO2_kg)

  REAL(lp), PARAMETER :: FireFracWood2Atmos      =  0.2_lp              !< Fraction of wood carbon emitted to the atm by fire.
                                                                        !< In JS3 called def_fire_fract_wood_2_atmos.
  !$ACC DECLARE COPYIN(FireFracWood2Atmos)
  REAL(lp), PARAMETER :: fract_wood_aboveGround  = 0.7_lp               !< Fraction of C above ground in wood pool (for separation
                                                                        !  of woody litter into above and below ground litter pools)
  !$ACC DECLARE COPYIN(fract_wood_aboveGround)
  REAL(lp), PARAMETER :: fract_green_aboveGround = 0.5_lp               !< Fraction of C above ground in green pool (for separation
                                                                        !  of green litter into above and below ground litter pools)
  !$ACC DECLARE COPYIN(fract_green_aboveGround)
  REAL(lp), PARAMETER :: sec_per_day             = 86400._lp            !< seconds per day
  !$ACC DECLARE COPYIN(sec_per_day)
  REAL(lp), PARAMETER :: days_per_year           = 365.25_lp
  !$ACC DECLARE COPYIN(days_per_year)
  REAL(lp), PARAMETER :: sec_per_year            = days_per_year * sec_per_day
  !$ACC DECLARE COPYIN(sec_per_year)

  ! Parameters for Yasso peat
  REAL(lp), parameter :: rho_YCpool_acro          = 3.5e4_lp !< Bulk density of acrotelm peat: 0.035 g/cm3, based on
                                                             !< Granberg (1999): rho upper 20cm: 0.02 g/cm^3  */
                                                             !<  rho 20-40cm: 0.05 g/cm3, > 40 cm: 0.1-0.12 g/cm3
  REAL(lp), parameter :: rho_YCpool_cato          = 9.1e4_lp !< Bulk density of catotelm peat: 0.091 g/cm3, based on
                                                             !< Clymo (1983); Turunen (2002)
  REAL(lp), parameter :: Cfrac_YCpool_acro        = 0.50_lp  !< C fraction in acrotelm peat Malmer and Wallen (2004)
  REAL(lp), parameter :: Cfrac_YCpool_cato        = 0.52_lp  !< Malmer and Wallen (2004)

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
    & LCT_ACID, LCT_WATER, LCT_ETHANOL, LCT_NONSOLUBLE,              &
    & LCT_ACID, LCT_WATER, LCT_ETHANOL, LCT_NONSOLUBLE, LCT_HUMUS,   &
    & LCT_ACID, LCT_WATER, LCT_ETHANOL, LCT_NONSOLUBLE,              &
    & LCT_ACID, LCT_WATER, LCT_ETHANOL, LCT_NONSOLUBLE, LCT_HUMUS  /)

  ! Yasso parameters (Tuomi 2011)
  !   reference decompositation rates for T = 0C and unlimited water availability in [yr-1]...
  REAL(lp), PARAMETER   :: ref_decomp_rate_acid       = -0.72_lp    !<  ... acid-soluble fraction
  REAL(lp), PARAMETER   :: ref_decomp_rate_water      = -5.9_lp     !<  ... water-soluble fraction
  REAL(lp), PARAMETER   :: ref_decomp_rate_ethanol    = -0.28_lp    !<  ... ethanol-soluble fraction
  REAL(lp), PARAMETER   :: ref_decomp_rate_nonsoluble = -0.031_lp   !<  ... nonsoluble-soluble fraction
  REAL(lp), PARAMETER   :: ref_decomp_rate_humus      = -0.0016_lp  !<  ... humus

  ! ------------------------------------------------------------------------------------------
  ! Parameters for vertical YASSO
  !

  ! decomposition modifier parameters for temperature and depth modifier
  REAL(lp), PARAMETER :: Q10     = 1.5_lp           !< Q10 for deep soil C decomposition, based on CLM
  REAL(lp), PARAMETER :: ref_dep = 1000._lp         !< reference depth to remove depth dependence of decomposition

  ! decomposition modifiers for oxic / anoxic soils
  REAL(lp), PARAMETER :: d_anox = 0.35_lp           !< decomposition modifier for anaerobic decomposition
  REAL(lp), PARAMETER :: d_oxic = 1.0_lp            !< decomposition modifier for oxic decomposition

  ! Yasso parameters (Tuomi 2011)
  ! Parameters for temperature dependence
  REAL(lp), PARAMETER   :: temp_p1     = 0.095_lp   !< first  order temperature dependence [1/C]
  REAL(lp), PARAMETER   :: temp_p2     = -0.0014_lp !< second order temperature dependence [1/C]
  ! Parameter for precipitation dependence
  REAL(lp), PARAMETER   :: precip_p1   = -1.21_lp   !< first  order precipitation dependence [yr/m]
  ! Parameters for size dependence
  REAL(lp), PARAMETER   :: size_p1     = -1.71_lp   !< first  order size dependence [1/cm]
  REAL(lp), PARAMETER   :: size_p2     = 0.86_lp    !< second order size dependence [1/cm]
  REAL(lp), PARAMETER   :: size_p3     = -0.306_lp  !< size dependence power [ ]
  REAL(lp), PARAMETER   :: d_size_leaf = 1.0_lp     !< leaf litter is not size dependent

  ! Parameters for moisture dependence, DAYCENT formulation
  REAL(lp), PARAMETER   :: moist_p_a = 0.575_lp
  REAL(lp), PARAMETER   :: moist_p_b = 1.5_lp
  REAL(lp), PARAMETER   :: moist_p_d = 3.0_lp

  REAL(lp), PARAMETER   :: minimum_layer_thickness = 0.05_lp !< minimum unfrozen soil layer thickness threshold [m]
  REAL(lp), PARAMETER   :: permafrost_threshold    = 1._lp   !< length threshold to consider a layer permafrost [yr]

  ! ------------------------------------------------------------------------------------------
  ! Ratio of C mass to biomass (Spitfire and emission of chemical species by fire)
  REAL(lp),  PARAMETER :: C_2_biomass_ratio = 0.5_lp !< Conversion factor between SOC and SOM (Pribyl, 2010)

  ! organic matter [litter/duff?] bulk densities for below above and ground [kg/m3],
  ! see Chojnacky (2009), 0'Donnell (2009), (https://philbusey.com/sports_turf/organic_matter.html#table), Ahn (2009)
  REAL(lp), PARAMETER ::  rho_org    = 100._lp
  REAL(lp), PARAMETER ::  rho_org_ag =  75._lp
  REAL(lp), PARAMETER ::  rho_org_bg = 150._lp

#endif
END MODULE mo_carbon_constants
