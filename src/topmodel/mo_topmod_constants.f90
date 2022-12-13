!>
!! @brief Contains constants for the soil and snow energy process
!!
!! @author
!!  Stiig Wilkenskjeld, MPI-M Hamburg
!!
!! @par Revision History
!! Adopted from JSBACH3 code by Thomas Kleinen 2022-12-06
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_top_constants
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
!  USE mo_jsb_physical_constants,  ONLY: ks, ki, clw, cs, ci, rhoh2o, dens_snow, dens_ice

  IMPLICIT NONE
  PUBLIC

<  REAL(wp), PARAMETER ::                    &

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_top_constants'

#endif
END MODULE mo_top_constants
