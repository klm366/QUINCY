!>
!! @brief Contains constants for the carbon processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-15)
!!
MODULE mo_carbon_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_physical_constants,   ONLY: amco2

  IMPLICIT NONE
   PUBLIC
   REAL(wp),  PARAMETER :: molarMassC_kg      = 12.01_wp * 1.e-3_wp  ! 12.01e-3_wp   -  Mass of 1 mod C   in kg
   REAL(wp),  PARAMETER :: molarMassCO2_kg    = amco2 * 1.e-3_wp     ! 44.0095e-3_wp -  Mass of 1 mol CO2 in kg

   REAL(wp),  PARAMETER :: FireFracWood2Atmos =  0.2_wp              ! Fraction of wood carbon emitted to the atm by fire.
                                                                     ! In JS3 called def_fire_frac_wood_2_atmos.
   REAL(wp), PARAMETER  :: frac_wood_aboveGround  = 0.7_wp           ! Fraction of Carbon above ground in wood pool (for separation
                                                                     !    of woody litter into above and below ground litter pools)
   REAL(wp), PARAMETER  :: frac_green_aboveGround = 0.5_wp           ! Fraction of Carbon above ground in green pool (for separation
                                                                     !    of green litter into above and below ground litter pools)

   REAL(wp), parameter :: sec_per_day             = 86400._wp        ! seconds per day
   REAL(wp), parameter :: days_per_year           = 365.25_wp
   REAL(wp), parameter :: sec_per_year            = days_per_year*sec_per_day
#endif
END MODULE mo_carbon_constants
