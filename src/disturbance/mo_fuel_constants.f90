!>
!! @brief Contains constants for the fuel process
!!
!! @author
!!  Rainer Schneck, MPI-M Hamburg
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Rainer Schneck (2015-04-15)
!! 2018-04-25, StW: Added constants for controlling different fuel algorithms
!!
MODULE mo_fuel_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  CHARACTER(len=*), PARAMETER :: fuel_names    = 'none jsb3 spitfire'
  INTEGER,          PARAMETER :: FUEL_NONE     = 0
  INTEGER,          PARAMETER :: FUEL_JSB3     = 1
  INTEGER,          PARAMETER :: FUEL_SPITFIRE = 2

#endif
END MODULE mo_fuel_constants
