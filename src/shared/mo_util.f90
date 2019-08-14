!>
!! @par Copyright
!!  This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!!  Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!
!! @brief diverse subroutines and functions for tile & forcing variables
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!!   First version, copied from JSBACH4 git repo, by Jan Engel (2016-03-23)
!!
MODULE mo_util
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: soil_depth_to_layers_2d
  PUBLIC :: one_of, toupper, int2string, tolower

  INTERFACE one_of
    MODULE PROCEDURE one_of_str
    MODULE PROCEDURE one_of_int
  END INTERFACE one_of

  CHARACTER(len=*), PARAMETER :: modname = 'mo_util'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> FUNCTION return the site-specific soil profile as vector of soil-layer thicknesses
  !! 
  !! in: soil_depth (until bedrock) and soil-layer thicknesses
  !! 
  !! out: site-specific soil-layer thicknesses, soil_depth_sl as 3D variable "dsl4jsb_var3D_onDomain()"
  !-----------------------------------------------------------------------------------------------------
  FUNCTION soil_depth_to_layers_2d(depth, dz) RESULT(depth_l)
  
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: depth(:,:)
    REAL(wp), INTENT(in) :: dz(:)
    REAL(wp)             :: depth_l(SIZE(depth,1), SIZE(dz), SIZE(depth,2))
    ! ---------------------------
    ! 0.2 Local
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
  !> conversion: lowercase -> uppercase
  !! 
  !! Copied from ICON module shared/mo_util_string
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION toupper (lowercase)
    IMPLICIT NONE
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

#endif
END MODULE mo_util
