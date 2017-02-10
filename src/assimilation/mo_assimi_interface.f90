!>
!! @brief Contains the interface to the assimilation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-19)
!!
MODULE mo_assimi_interface

#include "dsl_jsbach.inc"

  USE mo_kind,                          ONLY: wp
  USE mo_exception,                     ONLY: message, finish

  USE mo_assimi_update_canopy_fluxes,   ONLY: update_canopy_fluxes

  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: Register_assimi_tasks


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_assimi_interface'

!CONTAINS




END MODULE mo_assimi_interface
