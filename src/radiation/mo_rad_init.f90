!>
!! @brief Initialization of the the radiation memory.
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161213)
!!
MODULE mo_rad_init

#include "dsl_jsbach.inc"

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rad_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_rad_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Intialize radiation process (after memory has been set up)
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE rad_init(tile)

    USE mo_base_memory,             ONLY: t_slm_tile

    IMPLICIT NONE

    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)              :: tile       !< data structure for dynamic vegetation and soil variables
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':rad_init'


    CALL message(TRIM(routine), 'Currently nothing to do for '//routine)


  END SUBROUTINE rad_init

END MODULE mo_rad_init
