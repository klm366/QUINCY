!>
!! @brief Contains constants for the glacier processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_glacier_constants
#ifndef __NO_JSBACH__

  USE mo_physical_constants, ONLY: &
    rd,              &  !< Gas constant [J/K/kg]
    cpd,             &  !< Specific heat of dry air at constant pressure [J/K/kg]
    cpv,             &  !< Specific heat of water vapor at constant pressure [J/K/kg]
    clw,             &  !< Specific heat of liquid water [J/K/kg]
    cpvd1 => vtmpc2, &  !< cpv/cpd-1
    rhoh2o              !< Density of liquid water [kg/m3]
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
END MODULE mo_glacier_constants
