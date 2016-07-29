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

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish

  USE mo_assimi_update_canopy_fluxes

  IMPLICIT NONE

!  PRIVATE
!  PUBLIC :: Register_assimi_tasks


CONTAINS

  ! Task
!  CALL update_canopy_fluxes(tile, model%lctlib(tile%lcts%lib_id))





END MODULE mo_assimi_interface
