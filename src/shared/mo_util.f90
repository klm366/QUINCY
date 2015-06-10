!>
!! @brief Contains utility functions
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2014-03-17)
!!
MODULE mo_util
#ifndef __NO_JSBACH__

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: one_of, toupper, int2string

  CHARACTER(len=*), PARAMETER :: modname = 'mo_util'

CONTAINS

  !
  !> Conversion: Lowercase -> Uppercase
  !
  !  Copied from ICON module shared/mo_util_string
  !
  FUNCTION toupper (lowercase)
    CHARACTER(len=*), INTENT(in) :: lowercase
    CHARACTER(len=LEN_TRIM(lowercase)) :: toupper
    !
    INTEGER, PARAMETER :: idel = ICHAR('A')-ICHAR('a')
    INTEGER :: i
    !
    DO i = 1, LEN_TRIM(lowercase)
      IF (ICHAR(lowercase(i:i)) >= ICHAR('a') .AND. ICHAR(lowercase(i:i)) <= ICHAR('z')) THEN
        toupper(i:i) = CHAR( ICHAR(lowercase(i:i)) + idel )
      ELSE
        toupper(i:i) = lowercase(i:i)
      ENDIF
    ENDDO
    !
  END FUNCTION toupper

  !> Function for convenience 
  !
  !  If "in_str" is matching one of the arguments "arg(i)" return the
  !  index "i". Returns "-1" if none of the strings matches.
  !
  !  Copied from ICON module shared/mo_util_string
  !  
  FUNCTION one_of(in_str, arg)
    INTEGER :: one_of
    CHARACTER(len=*), INTENT(IN)           :: in_str    ! input string
    CHARACTER(len=*), INTENT(IN)           :: arg(:)
    ! local variables:
    INTEGER :: i
    
    one_of = -1
    DO i=1,SIZE(arg)
      IF (toupper(TRIM(in_str)) == toupper(TRIM(arg(i))))  one_of=i
    END DO
  END FUNCTION one_of

  !
  ! returns integer n as a string (often needed in printing messages)
  !
  ! Copied from ICON module shared/mo_util_string
  !  
  FUNCTION int2string(n, opt_fmt)
    CHARACTER(len=10) :: int2string ! result
    INTEGER, INTENT(in) :: n
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=10) :: fmt

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(I10)'
    END IF
    WRITE(int2string,fmt) n
    int2string = ADJUSTL(int2string)
    !
  END FUNCTION int2string

#endif
END MODULE mo_util
