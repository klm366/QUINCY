!
! This file was automatically generated from:
!    /scratch/snx3000/rschnur/icon-rc-testing-plus-mc8.test/externals/jsbach/src/turbulence/mo_turb_constants.f90
!    with md5 hash: fbe8a4d8f4a31e58e56b145abd328cef
!
!>
!! @brief Contains constants for the turbulence processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_turb_constants
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PARAMETER ::          &
    !
    ! Parameters for roughness length
    & blending_height = 100._wp,  & !< Blending height [m]
    & roughness_snow  = 0.001_wp, & !< Roughness length of snow-covered surfaces [m]
    & roughness_bare  = 0.005_wp    !< Roughness length of bare land surfaces [m]

#endif
END MODULE mo_turb_constants
