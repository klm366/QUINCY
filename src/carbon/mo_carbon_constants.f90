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

  IMPLICIT NONE
  PUBLIC

!  REAL(wp), PARAMETER ::                                     &

!     LaiMax = 8._wp,                      & ! Maximum LAI (used for nitrogen scaling) 
!     ZenithMin = 0.0174524_wp,            & ! Check for solar zenith angle > 89 degrees
!     SoilReflectivityParMin = 0.0_wp,     & ! Minimum soil reflectivity in PAR region

#endif
END MODULE mo_carbon_constants
