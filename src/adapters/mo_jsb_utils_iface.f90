!>
!! @brief Contains adapter to different utilities
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-02-29)
!!
MODULE mo_jsb_utils_iface
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text

#ifdef __ICON__
  USE mo_fortran_tools, ONLY: assign_if_present
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: assign_if_present

#ifndef __ICON__

  INTERFACE assign_if_present
    MODULE PROCEDURE assign_if_present_character
    MODULE PROCEDURE assign_if_present_logical
    MODULE PROCEDURE assign_if_present_integer
    MODULE PROCEDURE assign_if_present_integers
    MODULE PROCEDURE assign_if_present_real
  END INTERFACE

CONTAINS

  SUBROUTINE assign_if_present_character (y,x)
    CHARACTER(len=*) ,INTENT(inout)        :: y
    CHARACTER(len=*) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == ' ' )       RETURN
    y = x
  END SUBROUTINE assign_if_present_character
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_logical (y,x)
    LOGICAL ,INTENT(inout)        :: y
    LOGICAL ,INTENT(in) ,OPTIONAL :: x
    IF (PRESENT(x)) y = x
  END SUBROUTINE assign_if_present_logical
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integer (y,x)
    INTEGER ,INTENT(inout)        :: y
    INTEGER ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT. PRESENT(x)) RETURN
    IF ( x == -HUGE(x)  ) RETURN
    y = x
  END SUBROUTINE assign_if_present_integer
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_integers (y,x)
    INTEGER ,INTENT(inout)        :: y (:)
    INTEGER ,INTENT(in) ,OPTIONAL :: x (:)
    INTEGER :: n
    IF (PRESENT(x)) THEN
      n = MIN(SIZE(x), SIZE(y))
      y(1:n) = x(1:n)
    ENDIF
  END SUBROUTINE assign_if_present_integers
!------------------------------------------------------------------------------
  SUBROUTINE assign_if_present_real (y,x)
    REAL(wp), INTENT(inout)        :: y
    REAL(wp) ,INTENT(in) ,OPTIONAL :: x
    IF (.NOT.PRESENT(x)) RETURN
    IF ( x == -HUGE(x) ) RETURN
    y = x
  END SUBROUTINE assign_if_present_real

#endif

#endif
END MODULE mo_jsb_utils_iface
