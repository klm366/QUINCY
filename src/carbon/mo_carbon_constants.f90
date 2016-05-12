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
  USE mo_physical_constants,       ONLY: amd, amco2

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::                                     &
                         molarMassDryAir = amd   * 1.e-3_wp, &  ! 28.970e-3_wp. Mass [kg] per amount of substance [mol] of dry air
                         molarMassCO2    = amco2 * 1.e-3_wp    ! 44.0095e-3_wp. Mass [kg] per amount of substance [mol] of CO2

!     LaiMax = 8._wp,                      & ! Maximum LAI (used for nitrogen scaling) 
!     ZenithMin = 0.0174524_wp,            & ! Check for solar zenith angle > 89 degrees
!     SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region


    ! Factors that relates leaf internal CO2-concentration to CO2-concentration of ambient air: 
!    & FCI1C3        = 0.87_wp,              & ! For C3 plants
!    & FCI1C4        = 0.67_wp,              & ! For C4 plants
                      ! [g/mol] dry air
  
#endif
END MODULE mo_carbon_constants
