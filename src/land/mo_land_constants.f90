!!
!! @brief Contains general parameters required for land calculations 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_land_constants
 
  USE mo_kind
  ! physical constants
  REAL(wp), PARAMETER :: r_gas = 8.3144_wp ! universal gas constant (J / mol / K) 
  ! molar masses
  REAL(wp), SAVE :: molar_mass_C = 12.0112_wp   ! g C / mol C
  REAL(wp), SAVE :: molar_mass_C12 = 12.0_wp   ! g C / mol C
  REAL(wp), SAVE :: molar_mass_C13 = 13.0_wp   ! g C / mol C
  REAL(wp), SAVE :: molar_mass_C14 = 14.0_wp   ! g C / mol C
  REAL(wp), SAVE :: molar_mass_N = 14.0072_wp   ! g N / mol N
  REAL(wp), SAVE :: molar_mass_N14 = 14.0_wp   ! g N / mol N
  REAL(wp), SAVE :: molar_mass_N15 = 15.0_wp   ! g N / mol N
  REAL(wp), SAVE :: molar_mass_P = 30.997_wp ! g P / mol P
  ! conversion factors
  REAL(wp), SAVE :: sw2par = 2.3_wp            ! W / µmol / s
  REAL(wp), SAVE :: Tzero = 273.15_wp          ! melting point in K
  REAL(wp), SAVE :: standard_press_srf = 101300.0_wp ! standard pressure in Pa
  REAL(wp), SAVE :: eps_vpd = 0.622_wp         ! eps_vpd used for unit conversion Pa -> g g-1

  PUBLIC :: calc_spec_humidity_sat

CONTAINS

  FUNCTION calc_spec_humidity_sat(temp,press_srf) RESULT(q_sat)
        ! constants for the calculation of atmospheric water content at saturation
        ! from Campell and Norman's Introduction to Environmental Biophysics
     USE mo_kind, ONLY: wp
     IMPLICIT NONE 
     ! 0.1 Input
     REAL(wp) :: temp            ! temperature (K)
     REAL(wp) :: press_srf       ! atmospheric pressure (Pa) 
     ! 0.3 Output 
     REAL(wp) :: q_sat           ! saturating specific humidity (g/g) 
     ! 0.4 Local
     REAL(wp) :: a = 611.0_wp    ! Pa
     REAL(wp) :: b = 17.502_wp   ! einheitenlos
     REAL(wp) :: c = 240.97      ! °C
     REAL(wp) :: esat            ! saturating vapour pressure (Pa)

     ! calculating saturating vapour pressure [Pa] (Campbell & Norman)
     esat = a * exp(b * (temp - TZero)/(temp - TZero + c))
     ! conversion to saturating specific humidity [kg kg-1]
     q_sat = eps_vpd * esat/(press_srf - (1.0_wp - eps_vpd) * esat)

  END FUNCTION calc_spec_humidity_sat

END MODULE mo_land_constants
