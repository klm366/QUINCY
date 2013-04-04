!>
!! @brief Contains the routines for the soil processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2012-09-27)
!!
MODULE mo_soil_process

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: snow_fract

  CHARACTER(len=*), PARAMETER :: modname = 'mo_soil_process'

CONTAINS

  !< Roesch et al. 2002, Climate Dynamics
  !
  REAL(wp) ELEMENTAL FUNCTION snow_fract(snow, orodev)

    USE mo_soil_constants, ONLY: wsn2frac_const, wsn2frac_eps, wsn2frac_sigfac

    REAL(wp), INTENT(in) :: &
      & snow,   &
      & orodev

    snow_fract = wsn2frac_const * TANH(snow * 100._wp) &
                               & * SQRT(snow * 1000._wp / (snow * 1000._wp + wsn2frac_eps + wsn2frac_sigfac * orodev))

  END FUNCTION snow_fract

END MODULE mo_soil_process
