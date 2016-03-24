!>
!! @brief Contains constants for the srf processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_srf_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER :: &
    !
    ! Parameters for roughness length
    blending_height = 100._wp,  & !< Blending height [m]
    roughness_snow  = 0.001_wp, & !< Roughness length of snow-covered surfaces [m]
    roughness_bare  = 0.005_wp    !< Roughness length of bare land surfaces [m]

#endif
END MODULE mo_srf_constants
