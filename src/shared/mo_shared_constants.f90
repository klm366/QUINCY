!>
!! @brief Contains general parameters required for land calculations 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_shared_constants
 
  USE mo_kind,   ONLY: wp

  IMPLICIT NONE

  PUBLIC

  ! physical constants
  REAL(wp), SAVE :: &
    r_gas        ,  &           !< universal gas constant (J / mol / K) 
    r_gas_dryair ,  &           !< specific gas constant for dry air (J / kg / K) 
    Dwv          ,  &           !< diffusion coefficient for water vapour in air (m2 s-1) 
    Dco2         ,  &           !< diffusion coefficient for CO2 in air (m2 s-1) 
    Tzero                       !< melting point in K -> conversion factor for °C
  ! molar masses
  REAL(wp), SAVE :: &
    molar_mass_C    ,  &        !< g C / mol C
    molar_mass_C12  ,  &        !< g C / mol C
    molar_mass_C13  ,  &        !< g C / mol C
    molar_mass_C14  ,  &        !< g C / mol C
    molar_mass_N    ,  &        !< g N / mol N
    molar_mass_N14  ,  &        !< g N / mol N
    molar_mass_N15  ,  &        !< g N / mol N
    molar_mass_P                !< g P / mol P
  ! isotope standards
  REAL(wp), SAVE :: &
    PDB_STANDARD_C13   ,  &     !< 13C/12C PDB standard
    AAB_STANDARD_C14   ,  &     !< 0.95·0.238 Bq/gC is the absolute specific activity of 
                                 !! the radiocarbon standard (i.e. 95% of the activity of the OxA-I standard)
    AvogadroC          ,  &     !< Avogadro Number
    lambda_C14         ,  &     !< decay constant of C14 (1/s)
    scale_correction_c14        !< conversion from mol to pico-mol for C14 (for numerical precision) 


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from shared
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_shared_constants

    ! physical constants
    r_gas        = 8.3144_wp    
    r_gas_dryair = 287.058_wp   
    Dwv          = 0.0000242_wp 
    Dco2         = 0.0000151_wp 
    Tzero        = 273.15_wp 
    ! molar masses
    molar_mass_C   = 12.0112_wp 
    molar_mass_C12 = 12.0_wp    
    molar_mass_C13 = 13.0_wp    
    molar_mass_C14 = 14.0_wp    
    molar_mass_N   = 14.0072_wp 
    molar_mass_N14 = 14.0_wp    
    molar_mass_N15 = 15.0_wp    
    molar_mass_P   = 30.997_wp  
    ! isotope standards
    PDB_STANDARD_C13     = 0.0112372_wp   
    AAB_STANDARD_C14     = 0.2261_wp      
    AvogadroC            = 6.022_wp*1.e23     
    lambda_C14           = 3.8332_wp*1.e-12     
    scale_correction_c14 = 1000000000000._wp   

  END SUBROUTINE init_shared_constants

END MODULE mo_shared_constants
