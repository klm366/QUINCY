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

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_depth_to_layers_2d, one_of, toupper, int2string, real2string, logical2string

  INTERFACE one_of
    MODULE PROCEDURE one_of_str
    MODULE PROCEDURE one_of_int
  END INTERFACE one_of

  CHARACTER(len=*), PARAMETER :: modname = 'mo_util'

CONTAINS

  !!
  ! @brief      Compute soil layer depths from (fixed) layer thicknesses and total depth
  !
  ! @param      depth    Total soil depth
  ! @param      dz       Fixed layer thicknesses
  !
  ! @return     Layer thicknesses until soil depth
  !!
  FUNCTION soil_depth_to_layers_2d(depth, dz) RESULT(depth_l)

    REAL(wp), INTENT(in) :: depth(:,:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: depth_l(SIZE(depth,1), SIZE(dz), SIZE(depth,2))

    REAL(wp) :: bottom, bottom_m1
    INTEGER  :: nsoil, i

    CHARACTER(len=*), PARAMETER :: routine = modname//':soil_depth_to_layer_2d'

    nsoil = SIZE(dz)

    depth_l(:,:,:) = 0._wp
    bottom = 0._wp
    bottom_m1 = 0._wp
    DO i=1,nsoil
      bottom = bottom + dz(i)
      WHERE (depth(:,:) >= bottom)
        depth_l(:,i,:) = dz(i)
      ELSEWHERE (depth(:,:) > bottom_m1)
        depth_l(:,i,:) = depth(:,:) - bottom_m1
      END WHERE
      bottom_m1 = bottom_m1 + dz(i)
    END DO

  END FUNCTION soil_depth_to_layers_2d
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
  FUNCTION one_of_str(in_str, arg)
    INTEGER :: one_of_str
    CHARACTER(len=*), INTENT(IN)           :: in_str    ! input string
    CHARACTER(len=*), INTENT(IN)           :: arg(:)
    ! local variables:
    INTEGER :: i

    one_of_str = -1
    DO i=1,SIZE(arg)
      IF (toupper(TRIM(in_str)) == toupper(TRIM(arg(i))))  one_of_str=i
    END DO
  END FUNCTION one_of_str

  FUNCTION one_of_int(in_int, arg)
    INTEGER :: one_of_int
    INTEGER, INTENT(IN)           :: in_int    ! input integer
    INTEGER, INTENT(IN)           :: arg(:)
    ! local variables:
    INTEGER :: i

    one_of_int = -1
    DO i=1,SIZE(arg)
      IF (in_int == arg(i))  one_of_int=i
    END DO
  END FUNCTION one_of_int

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

  FUNCTION real2string(x, opt_fmt)
    CHARACTER(len=10) :: real2string ! result
    REAL(wp), INTENT(in) :: x
    CHARACTER(len=*), INTENT(in), OPTIONAL :: opt_fmt
    !
    CHARACTER(len=10) :: fmt

    IF (PRESENT(opt_fmt)) THEN
      fmt = opt_fmt
    ELSE
      fmt = '(E10.4)'
    END IF
    WRITE(real2string,fmt) x
    real2string = ADJUSTL(real2string)
    !
  END FUNCTION real2string

  FUNCTION logical2string(value)
    CHARACTER(len=7) :: logical2string ! result
    LOGICAL, INTENT(in) :: value
    !
    IF (value) THEN
      logical2string = '.TRUE.'
    ELSE
      logical2string = '.FALSE.'
    END IF
    !
  END FUNCTION logical2string

#endif
END MODULE mo_util
