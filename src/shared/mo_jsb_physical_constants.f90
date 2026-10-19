!> Contains physical constants
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
!> The module provides physical constants, either directly defined here or from the host model via the constants iface adapter.
!>
MODULE mo_jsb_physical_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: lp

  ! Most constants are used from host model (ICON) via the adapter file mo_physical_constants_iface in mo_util_jsbach.f90:
  !   src/shared/mo_util_jsbach.f90 -> mo_physical_constants_iface -> mo_physical_constants.f90
  USE mo_physical_constants_iface ! TODO: add "ONLY" with an explicit list

  IMPLICIT NONE
  PUBLIC

  REAL(lp), PARAMETER :: &
    & dens_soil          = 2700._lp,   & !< Density of soil particulate (Hillel 1980, Johansen 1977) [kg/m3]
    & dens_org           = 1300._lp,   & !< Density of organic part of soil                          [kg/m3]
    & dens_ice           = rhoi,       & !< Density of frozen water in soil (ice)                    [kg/m3]
    & dens_snow          = rhos,       & !< Density of snow                                          [kg/m3] (300.0_lp)
    & dens_snow_min      = 50._lp        !< Minimum density of snow (for fresh snow)                 [kg/m3]

  REAL(lp), PARAMETER :: &
    molarMassDryAir = amd   * 1.e-3_lp, & !< Mass [kg] per amount of substance [mol] of dry air (28.970e-3_lp)
    molarMassCO2    = amco2 * 1.e-3_lp    !< Mass [kg] per amount of substance [mol] of CO2 (44.0095e-3_lp)
  !$ACC DECLARE COPYIN(molarMassDryAir, molarMassCO2)

  REAL(lp), PARAMETER :: von_karman = 0.4_lp    !< von Karmann constant for evaporation

  !! quincy - physical constants
  REAL(lp), PARAMETER :: &
    & r_gas        = argas            ,  &      !< Universal gas constant (J / mol / K)
    & r_gas_dryair = rd               ,  &      !< Specific gas constant for dry air (J / kg / K)
    & r_gas_vapour = rv               ,  &      !< Specific gas constant for water vapour (J / kg / K)
    & Dwv          = 0.0000242_lp     ,  &      !< Diffusion coefficient for water vapour in air (m2 s-1)
    & Dco2         = 0.0000151_lp     ,  &      !< Diffusion coefficient for CO2 in air (m2 s-1)
    & Tzero        = 273.15_lp        ,  &      !< The value of Tzero (QUINCY) is identical with tmelt (ECHAM/ICON/jsbach4)
                                                !<   but in reality tmelt is not always at Tzero [K]
    & rain_as_snow_temp_thresh = Tzero + 2.5_lp !< Threshold below which rain falls as snow - source: OCN-weather generator

  REAL(lp), PARAMETER :: &
    & LatentHeatVaporization = alv , &          !< [J/kg]   latent heat for vaporisation
    & LatentHeatSublimation  = als              !< [J/kg]   latent heat for sublimation
  ! molar masses
  REAL(lp), PARAMETER :: &
    & molar_mass_C    = 12.0112_lp    ,  &      !< g C / mol C
    & molar_mass_C12  = 12.0_lp       ,  &      !< g C / mol C
    & molar_mass_C13  = 13.0_lp       ,  &      !< g C / mol C
    & molar_mass_C14  = 14.0_lp       ,  &      !< g C / mol C
    & molar_mass_N    = 14.0072_lp    ,  &      !< g N / mol N
    & molar_mass_N14  = 14.0_lp       ,  &      !< g N / mol N
    & molar_mass_N15  = 15.0_lp       ,  &      !< g N / mol N
    & molar_mass_P    = 30.997_lp               !< g P / mol P
  ! isotope standards
  REAL(lp), PARAMETER :: &
    & PDB_STANDARD_C13           = 0.0112372_lp                    , &  !< 13C/12C PDB standard
    & AAB_STANDARD_C14           = 0.2261_lp                       , &  !< 0.95*0.238 Bq/gC is the absolute specific activity of
                                                                        !<   the radiocarbon standard (i.e. 95% of the activity of the OxA-I standard)
    & ATM_STANDARD_N15           = 0.0036765_lp                    , &  !< 15N/14N atmospheric standard, Robinson 2001, TREE
    & AvogadroC                  = 6.022_lp * 1.e23_lp             , &  !< Avogadro Number
    & lambda_C14                 = 3.8332_lp * 1.e-12_lp           , &  !< Decay constant of C14 (1/s)
    & AvogadroC_times_lambda_C14 = 6.022_lp * 3.8332_lp * 1.e11_lp , &  !< Avogadro Number * decay constant of C14 to
                                                                        !<   reduce the use of the very large/small AvogadroC/lambda_C14
    & scale_correction_c14       = 1.e12_lp                             !< Conversion from mol to pico-mol for C14 (for numerical precision)

#endif
END MODULE mo_jsb_physical_constants
