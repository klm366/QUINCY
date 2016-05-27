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
  REAL(wp), PARAMETER :: r_gas     = 8.3144_wp       ! universal gas constant (J / mol / K) 
  REAL(wp), PARAMETER :: r_gas_dryair = 287.058_wp   ! specific gas constant for dry air (J / kg / K) 
  ! molar masses
  REAL(wp), SAVE :: molar_mass_C   = 12.0112_wp  ! g C / mol C
  REAL(wp), SAVE :: molar_mass_C12 = 12.0_wp     ! g C / mol C
  REAL(wp), SAVE :: molar_mass_C13 = 13.0_wp     ! g C / mol C
  REAL(wp), SAVE :: molar_mass_C14 = 14.0_wp     ! g C / mol C
  REAL(wp), SAVE :: molar_mass_N   = 14.0072_wp  ! g N / mol N
  REAL(wp), SAVE :: molar_mass_N14 = 14.0_wp     ! g N / mol N
  REAL(wp), SAVE :: molar_mass_N15 = 15.0_wp     ! g N / mol N
  REAL(wp), SAVE :: molar_mass_P   = 30.997_wp   ! g P / mol P


END MODULE mo_shared_constants
