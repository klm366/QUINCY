!>
!! @brief Contains constants for the natural disturbances processes
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-15)
!! 2018-04-25, StW: Added used constants and names - removed unused
!!
MODULE mo_disturbance_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! Constants for controlling windbreak algorithms
  CHARACTER(len=*), PARAMETER :: windbreak_names = 'none jsb3'
  INTEGER,          PARAMETER :: WINDBREAK_NONE  = 0
  INTEGER,          PARAMETER :: WINDBREAK_JSB3  = 1

  ! Constants for controlling fire algorithms
  CHARACTER(len=*), PARAMETER :: fire_names      = 'none jsb3 spitfire'
  INTEGER,          PARAMETER :: FIRE_NONE       = 0
  INTEGER,          PARAMETER :: FIRE_JSB3       = 1
  INTEGER,          PARAMETER :: FIRE_SPITFIRE   = 2
  INTEGER,          PARAMETER :: FUEL_JSB3       = FIRE_JSB3
  INTEGER,          PARAMETER :: FUEL_SPITFIRE   = FIRE_SPITFIRE

  ! Constants regarding temporal smooting of disturbance relevant forcing field
  ! StW: would it make sense to move these to jsb_disturb_fire_jsb3/jsb_disturb_windbreak_jsb3 namelists respectively?
  REAL(wp), PARAMETER         :: persist_rel_hum  = 0.95_wp   ! factor concerning the smoothing of relative air humidity (~14 days)
  REAL(wp), PARAMETER         :: persist_wind_10m = 0.9995_wp ! factor concerning the smoothing of daily maximum wind speed

#endif
END MODULE mo_disturbance_constants
