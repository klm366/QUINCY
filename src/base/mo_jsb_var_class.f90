MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_p

  ENUM, BIND(C)
    ENUMERATOR :: REAL2D=1, READ3D
  END ENUM

  TYPE, ABSTRACT :: t_jsb_var
!!$    INTEGER           :: type = REAL2D
  END TYPE t_jsb_var

  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), POINTER :: ptr(:,:) => NULL()
  END TYPE t_jsb_var_real2d

  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(wp), POINTER :: ptr(:,:,:) => NULL()
  END TYPE t_jsb_var_real3d

  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p
  END TYPE t_jsb_var_p

!!$  INTERFACE t_jsb_var
!!$    PROCEDURE :: Create_jsb_var
!!$  END INTERFACE t_jsb_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'

CONTAINS

!!$  FUNCTION Create_jsb_var(i_shape, i_type) RESULT(return_ptr)
!!$
!!$    INTEGER, INTENT(in) :: i_shape(2)
!!$    INTEGER, INTENT(in), OPTIONAL :: i_type
!!$    TYPE(t_jsb_var), POINTER     :: return_ptr
!!$
!!$    ALLOCATE(return_ptr)
!!$
!!$    IF (SIZE(i_shape) == 2) THEN
!!$      return_value%type = REAL2D
!!$      ELSE IF (SIZE(i_shape) == 3) THEN
!!$        return_value%type = REAL3D
!!$      END IF
!!$    END IF
!!$
!!$    SELECT CASE (return_value%type)
!!$    CASE (REAL2D)
!!$      ALLOCATE(return_value%real2d_ptr(i_shape(1), i_shape(2)))
!!$    CASE (REAL3D)
!!$      ALLOCATE(return_value%real3d_ptr(i_shape(1), i_shape(2), i_shape(3)))
!!$    END SELECT
!!$
!!$  END FUNCTION Create_jsb_var

#endif
END MODULE mo_jsb_var_class
