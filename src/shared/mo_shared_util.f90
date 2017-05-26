!>
!! @brief diverse subroutines and functions for tile & forcing variables
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-23)
!! 2nd version by Jan Engel (2016-07-15); all content removed, see SVN r198 & trac ticket #63
!!
MODULE mo_shared_util

  USE mo_kind,         ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tolower, int2string

  CHARACTER(len=*), PARAMETER :: modname = 'mo_shared_util'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> conversion: uppercase -> lowercase
  !! 
  !! implemented based on mo_util:toupper in jsbach4.1
  !! 
  !! in jsbach4.1 the mo_util:toupper is: Copied from ICON module shared/mo\_util\_string
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION tolower(string) RESULT(stringlower)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut    
    CHARACTER(len=*), INTENT(in) :: string
    CHARACTER(len=LEN_TRIM(string)) :: stringlower
    ! ---------------------------
    ! 0.2 Local
    INTEGER, PARAMETER :: idel = ICHAR('a')-ICHAR('A')  ! = 32
    INTEGER :: i


    DO i = 1, LEN_TRIM(string)
      IF (ICHAR(string(i:i)) >= ICHAR('A') .AND. ICHAR(string(i:i)) <= ICHAR('Z')) THEN
        stringlower(i:i) = CHAR( ICHAR(string(i:i)) + idel )
      ELSE
        stringlower(i:i) = string(i:i)
      ENDIF
    ENDDO
  END FUNCTION tolower


  !-----------------------------------------------------------------------------------------------------
  !> returns integer n as a string (often needed in printing messages)
  !!
  !! Copied from jsbach4.1 module shared/mo\_util  \n
  !! In jsbach4.1 this function was: Copied from ICON module shared/mo\_util\_string
  !!  
  !-----------------------------------------------------------------------------------------------------
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



END MODULE mo_shared_util
