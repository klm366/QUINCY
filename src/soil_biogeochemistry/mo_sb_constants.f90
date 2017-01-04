!>
!! @brief Contains constants for the soil_biogeochemistry process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!!
MODULE mo_sb_constants

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC


!  REAL(wp), SAVE :: &
!    parameter-here

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil_biogeochemistry
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_sb_constants

    IMPLICIT NONE
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sb_constants'



  END SUBROUTINE init_sb_constants

END MODULE mo_sb_constants
