!! @brief
!!  Contains parameter for implementation of some specific model code
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!!  First version by Jan Engel (20190114)
!!
MODULE mo_jsb_impl_constants
 
  USE mo_kind,                  ONLY: wp

  IMPLICIT NONE
  PUBLIC


  ! using the functionality of LOGICAL var with REAL var
  REAL(wp), PARAMETER :: false            = 0.0_wp    !< zero
  REAL(wp), PARAMETER :: true             = 1.0_wp    !< one
  REAL(wp), PARAMETER :: test_false_true  = 0.5_wp    !< 0.5, inbetween false and true; 'IF(var < test_false_true) -> FALSE'


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_jsb_impl_constants'

! CONTAINS



END MODULE mo_jsb_impl_constants
