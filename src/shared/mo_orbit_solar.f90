!>
!! @brief Contains some functions for orbital and solar parameters
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2022-08-26)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_orbit_solar
#ifndef __NO_JSBACH__

  USE mo_jsb_orbit_solar_iface, ONLY: inquire_declination, compute_cos_zenith_angle

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: inquire_declination, compute_cos_zenith_angle

  CHARACTER(len=*), PARAMETER :: modname = 'mo_orbit_solar'

! CONTAINS

#endif
END MODULE mo_orbit_solar
