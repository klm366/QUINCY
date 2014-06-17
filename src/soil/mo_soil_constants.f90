!>
!! @brief Contains constants for the soil processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_soil_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER :: &
                                       !
    & ws_crit_fract  = 0.75_wp, &      !< Fraction of field capacity at which soil starts to dry,
                                       !  i.e. at which water stress starts to be above zero
    & ws_wilt_fract  = 0.35_wp, &      !< Fraction of field capacity at which soil is at wilting point,
                                       !  i.e. when no transpiration occurs anymore

    ! Parameters used for snow cover fraction
    & wsn2frac_eps    = 1.E-12_wp,  &
    & wsn2frac_sigfac = 0.15_wp,    &
    & wsn2frac_const  = 0.95_wp        !< Inverse of snow height when snow is considered to completely cover the ground

#endif
END MODULE mo_soil_constants
