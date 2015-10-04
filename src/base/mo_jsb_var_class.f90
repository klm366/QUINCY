MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var_p

  ENUM, BIND(C)
    ENUMERATOR :: REAL2D=1, READ3D
  END ENUM

  TYPE, ABSTRACT :: t_jsb_var
!!$    INTEGER           :: type = REAL2D
    CHARACTER(len=5) :: inout = 'IN'
  CONTAINS
    PROCEDURE :: Set_to_OUT => Set_var_to_OUT
    PROCEDURE :: Set_to_IN  => Set_var_to_IN
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

  SUBROUTINE Set_var_to_OUT(this)

    CLASS(t_jsb_var), INTENT(inout) :: this

    this%inout = "OUT"

  END SUBROUTINE Set_var_to_OUT

  SUBROUTINE Set_var_to_IN(this)

    CLASS(t_jsb_var), INTENT(inout) :: this

    this%inout = "IN"

  END SUBROUTINE Set_var_to_IN

#endif
END MODULE mo_jsb_var_class
