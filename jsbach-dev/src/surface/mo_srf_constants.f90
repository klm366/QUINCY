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
    & eps = 0.1_wp    ! Time filtering coefficient for surface temperature

END MODULE mo_srf_constants
