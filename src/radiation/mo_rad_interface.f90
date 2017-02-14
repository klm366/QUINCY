!>
!! @brief Contains the interface to the radiation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-19)
!!
MODULE mo_rad_interface

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: message, finish

  USE mo_rad_update_radiation,  ONLY: update_radiation 



  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: Register_rad_tasks

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_rad_interface'

!CONTAINS





END MODULE mo_rad_interface
