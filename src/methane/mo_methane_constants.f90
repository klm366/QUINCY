!>
!! @brief Constants for the methane production and transport module.
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen, 2023-03-21
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_methane_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_physical_constants,  ONLY: molarMassCO2, amch4, amo2, molar_mass_C

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: emission_factors_for_methane
  PUBLIC :: max_em_vegtypes, em_vegtype_names, max_specietypes, specietype_names
  PUBLIC :: LT_UNDEFINED, LT_TROPICAL_FOREST, LT_EXTRATROPICAL_FOREST, LT_SAVANNAH, LT_AGRICULTURE, LT_PEAT
  PUBLIC :: SPECIE_CO2, SPECIE_CH4, SPECIE_NOx, SPECIE_CO, SPECIE_ISOPRENE
  PUBLIC :: T_ref
  PUBLIC :: t_specie_const, CO2_consts, CH4_consts, O2_consts
  PUBLIC :: minimum_porevol, C_2_biomass_ratio
  PUBLIC :: molarMassCO2_kg, molarMassCH4_kg, molarMassC_kg

  CHARACTER(len=*), PARAMETER :: modname = 'mo_methane_constants'

  ! Types and constants for specific chemical species
  ! Diffusivity constants expect temperature in degC
  TYPE t_specie_const
    CHARACTER (LEN=10)         :: specie_name  ! Name of chemical specie
    REAL(wp)                   :: molarMass_kg ! Molar mass of the speice [kg mol-1]
    REAL(wp)                   :: henry_k_inv  ! inverted Henry constant at 298K [L atm mol-1]
    REAL(wp)                   :: henry_const  ! Coefficient in Henry's law (Sander 1999) [K-1]
    REAL(wp)                   :: D_air_0      ! Constant part of diffusivity in air (Lerman, 1979) [cm2 s-1]
    REAL(wp)                   :: D_air_1      ! Linear temperature dependency of diffusivity in air (Lerman, 1979) [cm2 s-1 K-1]
    REAL(wp)                   :: D_water_0    ! Constant part of diffusivity in water (Broecker and Peng, 1974) [cm2 s-1]
    REAL(wp)                   :: D_water_1    ! Linear temperature dependency of diffusivity in water (B&P, 1974) [cm2 s-1 K-1]
    REAL(wp)                   :: D_water_2    ! Square temperature dependency of diffusivity in water (B&P, 1974) [cm2 s-1 K-2]
    REAL(wp)                   :: D_water_scale! Scaling of water diffusivity  (Broecker and Peng, 1974) []
  END TYPE t_specie_const

  ! Name  molarMass_kg   henry_k_inv henry_const   D_air_0     D_air_1   D_water_0   D_water_1     D_water_2 D_water_scale
  TYPE (t_specie_const), PARAMETER :: CO2_consts = t_specie_const( &
    'CO2',molarMassCO2,      29.4_wp,  2400.0_wp, 0.1325_wp, 0.00090_wp, 0.9390_wp, 0.02671_wp, 0.0004095_wp,         1e-5 )
!    'CO2',molarMassCO2,      29.4_wp,  2400.0_wp, 0.1105_wp, 0.00150_wp, 0.9550_wp, 0.02550_wp, 0.00050714_wp,        1e-5 )
! air^2                                            0.1381_wp, 0.00100_wp, 1.0859e-6_wp
  TYPE (t_specie_const), PARAMETER :: CH4_consts = t_specie_const( &
    'CH4',amch4/1000.0_wp, 714.29_wp,  1600.0_wp, 0.1875_wp, 0.00130_wp, 0.9798_wp, 0.02986_wp, 0.0004381_wp,         1e-5 )
!    'CH4',amch4/1000.0_wp, 714.29_wp,  1600.0_wp, 0.1489_wp, 0.00180_wp, 0.9531_wp, 0.02480_wp, 0.00042857_wp,        1e-5 )
! air^2                                            0.1851_wp, 0.00120_wp, 1.4233e-6_wp
  TYPE (t_specie_const), PARAMETER :: O2_consts  = t_specie_const( &
    'O2', amo2 /1000.0_wp, 769.23_wp,  1500.0_wp, 0.1759_wp, 0.00117_wp, 1.1720_wp, 0.03443_wp, 0.0005048_wp,         1e-5 )
!    'O2', amo2 /1000.0_wp, 769.23_wp,  1500.0_wp, 0.1759_wp, 0.00117_wp, 1.0700_wp, 0.01900_wp, 0.0014000_wp,         1e-5 )
! https://www.engineeringtoolbox.com/air-diffusion-coefficient-gas-mixture-temperature-d_2010.html (air: CO2, CH4)
! https://www.engineeringtoolbox.com/diffusion-coefficients-d_1404.html (water)

  ! All emission factors from Kaiser et al., 2012, doi:10.5194/bg-9-527-2012
  ! All in [g(specie)/kg(dry matter)]

  REAL(wp), PARAMETER     :: T_ref = 288.0_wp      ! reference temperature for methane processes [K]

  ! Control parameters for emission factors
  INTEGER,      PARAMETER :: max_em_vegtypes = 5
  CHARACTER(*), PARAMETER :: em_vegtype_names = 'tropical forest,extratropical forest,savannah,agriculture,peat'
  ! LT = Landscape type (not to be confused with vegetation type)
  ENUM, BIND(C)
    ENUMERATOR :: LT_UNDEFINED = 0, LT_TROPICAL_FOREST, LT_EXTRATROPICAL_FOREST, LT_SAVANNAH, LT_AGRICULTURE, LT_PEAT
  END ENUM
  INTEGER,      PARAMETER :: max_specietypes = 5
  CHARACTER(*), PARAMETER :: specietype_names = 'CO2,CH4,NOx,CO,isoprene'
  ENUM, BIND(C)
    ENUMERATOR :: SPECIE_CO2 = 1, SPECIE_CH4, SPECIE_NOx, SPECIE_CO, SPECIE_ISOPRENE
  END ENUM

  ! CO2 emission factors
  REAL(wp), PARAMETER :: efactor_CO2_SA = 1646.0_wp   ! emission factor for savannah
  REAL(wp), PARAMETER :: efactor_CO2_TF = 1626.0_wp   ! emission factor for tropical forest
  REAL(wp), PARAMETER :: efactor_CO2_EF = 1572.0_wp   ! emission factor for extratropical forest
  REAL(wp), PARAMETER :: efactor_CO2_AG = 1308.0_wp   ! emission factor for agriculture
  REAL(wp), PARAMETER :: efactor_CO2_PT = 1703.0_wp   ! emission factor for peat

  ! CH4 emission factors
  REAL(wp), PARAMETER :: efactor_CH4_SA =    2.2_wp   ! emission factor for savannah
  REAL(wp), PARAMETER :: efactor_CH4_TF =    6.6_wp   ! emission factor for tropical forest
  REAL(wp), PARAMETER :: efactor_CH4_EF =    4.8_wp   ! emission factor for extratropical forest
  REAL(wp), PARAMETER :: efactor_CH4_AG =    8.4_wp   ! emission factor for agriculture
  REAL(wp), PARAMETER :: efactor_CH4_PT =   20.8_wp   ! emission factor for peat

  ! NOx emission factors
  REAL(wp), PARAMETER :: efactor_NO_SA =     2.1_wp   ! emission factor for savannah
  REAL(wp), PARAMETER :: efactor_NO_TF =     2.3_wp   ! emission factor for tropical forest
  REAL(wp), PARAMETER :: efactor_NO_EF =     3.4_wp   ! emission factor for extratropical forest
  REAL(wp), PARAMETER :: efactor_NO_AG =     2.3_wp   ! emission factor for agriculture
  REAL(wp), PARAMETER :: efactor_NO_PT =     1.0_wp   ! emission factor for peat

  ! CO emission factors
  REAL(wp), PARAMETER :: efactor_CO_SA =    61.0_wp   ! emission factor for savannah
  REAL(wp), PARAMETER :: efactor_CO_TF =   101.0_wp   ! emission factor for tropical forest
  REAL(wp), PARAMETER :: efactor_CO_EF =   106.0_wp   ! emission factor for extratropical forest
  REAL(wp), PARAMETER :: efactor_CO_AG =    92.0_wp   ! emission factor for agriculture
  REAL(wp), PARAMETER :: efactor_CO_PT =   210.0_wp   ! emission factor for peat

  ! isoprene emission
  REAL(wp), PARAMETER :: efactor_iso_SA =    0.026_wp ! emission factor for savannah
  REAL(wp), PARAMETER :: efactor_iso_TF =    0.22_wp  ! emission factor for tropical forest
  REAL(wp), PARAMETER :: efactor_iso_EF =    0.11_wp  ! emission factor for extratropical forest
  REAL(wp), PARAMETER :: efactor_iso_AG =    0.40_wp  ! emission factor for agriculture
  REAL(wp), PARAMETER :: efactor_iso_PT =    1.4_wp   ! emission factor for peat

  REAL(wp), PARAMETER :: emission_factors_for_methane(max_specietypes,max_em_vegtypes) = RESHAPE( (/ &
   & efactor_CO2_TF, efactor_CH4_TF, efactor_NO_TF, efactor_CO_TF, efactor_iso_TF, &
   & efactor_CO2_EF, efactor_CH4_EF, efactor_NO_EF, efactor_CO_EF, efactor_iso_EF, &
   & efactor_CO2_SA, efactor_CH4_SA, efactor_NO_SA, efactor_CO_SA, efactor_iso_SA, &
   & efactor_CO2_AG, efactor_CH4_AG, efactor_NO_AG, efactor_CO_AG, efactor_iso_AG, &
   & efactor_CO2_PT, efactor_CH4_PT, efactor_NO_PT, efactor_CO_PT, efactor_iso_PT  &
   & /), SHAPE(emission_factors_for_methane) )

  REAL(wp), PARAMETER :: minimum_porevol = 0.01_wp    ! minimum pore volume to consider

  ! Ratio of C mass to biomass (Spitfire and emission of chemical species by fire)
  REAL(wp),  PARAMETER :: C_2_biomass_ratio = 0.5_wp

  REAL(wp), PARAMETER :: molarMassCO2_kg = molarMassCO2
  REAL(wp), PARAMETER :: molarMassCH4_kg = amch4        * 1.e-3_wp
  REAL(wp), PARAMETER :: molarMassC_kg   = molar_mass_C * 1.e-3_wp

#endif
END MODULE mo_methane_constants
