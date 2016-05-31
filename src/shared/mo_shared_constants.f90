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

  ! physical constants
  REAL(wp), PARAMETER :: r_gas        = 8.3144_wp        !< universal gas constant (J / mol / K) 
  REAL(wp), PARAMETER :: r_gas_dryair = 287.058_wp       !< specific gas constant for dry air (J / kg / K) 
  REAL(wp), PARAMETER :: Dwv          = 0.0000242_wp     !< diffusion coefficient for water vapour in air (m2 s-1) 
  REAL(wp), PARAMETER :: Dco2         = 0.0000151_wp     !< diffusion coefficient for CO2 in air (m2 s-1) 
  ! molar masses
  REAL(wp), SAVE :: molar_mass_C   = 12.0112_wp          !< g C / mol C
  REAL(wp), SAVE :: molar_mass_C12 = 12.0_wp             !< g C / mol C
  REAL(wp), SAVE :: molar_mass_C13 = 13.0_wp             !< g C / mol C
  REAL(wp), SAVE :: molar_mass_C14 = 14.0_wp             !< g C / mol C
  REAL(wp), SAVE :: molar_mass_N   = 14.0072_wp          !< g N / mol N
  REAL(wp), SAVE :: molar_mass_N14 = 14.0_wp             !< g N / mol N
  REAL(wp), SAVE :: molar_mass_N15 = 15.0_wp             !< g N / mol N
  REAL(wp), SAVE :: molar_mass_P   = 30.997_wp           !< g P / mol P
  ! isotope standards
  REAL(wp), SAVE :: PDB_STANDARD_C13 = 0.0112372_wp      !< 13C/12C PDB standard
  REAL(wp), SAVE :: AAB_STANDARD_C14 = 0.2261_wp         !< 0.95·0.238 Bq/gC is the absolute specific activity of 
                                                         !< the radiocarbon standard (i.e. 95% of the activity of the OxA-I standard)
  REAL(wp), SAVE :: AvogadroC = 6.022_wp*1.e23           !< Avogadro Number
  REAL(wp), SAVE :: lambda_C14 = 3.8332_wp*1.e-12        !< decay constant of C14 (1/s)
  REAL(wp), SAVE :: scale_correction_c14 = 1000000000000._wp   !< conversion from mol to pico-mol for C14 (for numerical precision) 

END MODULE mo_shared_constants
