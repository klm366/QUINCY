MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp, dp
  USE mo_exception, ONLY: finish, message

  USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var, t_jsb_var_p, REAL2D, REAL3D

  ENUM, BIND(C)
    ENUMERATOR :: REAL2D=1, REAL3D
  END ENUM

  TYPE, ABSTRACT :: t_jsb_var
    REAL(wp), ALLOCATABLE      :: field(:,:,:,:,:)
    CHARACTER(len=VARNAME_LEN) :: name            ! Name of variable
    CHARACTER(len=VARNAME_LEN) :: full_name       ! Full name of variable in varlist, including prefix and suffix
    INTEGER                    :: type            ! REAL2D / REAL3D
    INTEGER                    :: owner_proc_id = -1
    INTEGER, ALLOCATABLE       :: owner_tile_path(:)
    INTEGER, ALLOCATABLE       :: child_idx(:)    ! Position-numbers of same variable for for the child tiles (according to mem%vars(:)).
    REAL(dp)                   :: missval
    LOGICAL                    :: l_aggregate_all = .FALSE. ! When aggregating, consider fractions of all child tiles even
                                                  ! of those where variable doesn't exist
  END TYPE t_jsb_var

!#pragma GCC diagnostic push
!#pragma GCC diagnostic ignored "-Wsurprising"

  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), POINTER :: ptr(:,:) => NULL()
  CONTAINS
    ! PROCEDURE :: Add_real2d                       ! Add two variables of type t_jsb_var_real2d
    ! PROCEDURE :: Add_real2d_and_scalar            ! Add varialbe of type t_jsb_var_real2d with REAL scalar
    ! PROCEDURE :: Add_real2d_and_field             ! Add varialbe of type t_jsb_var_real2d with REAL 2-d field
    ! GENERIC   :: OPERATOR(+) => Add_real2d, Add_real2d_and_scalar, Add_real2d_and_field
    ! PROCEDURE :: Multiply_real2d_with_field       ! Multiply varialbe of type t_jsb_var_real2d with REAL 2-d field
    ! GENERIC   :: OPERATOR(*) => Multiply_real2d_with_field
    FINAL     :: finalize_jsb_var_real2d
  END TYPE t_jsb_var_real2d

  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(wp), POINTER :: ptr(:,:,:) => NULL()
  CONTAINS
    ! PROCEDURE :: Add_real3d
    ! PROCEDURE :: Add_real3d_and_scalar            ! Add varialbe of type t_jsb_var_real3d with REAL scalar
    ! PROCEDURE :: Add_real3d_and_field             ! Add varialbe of type t_jsb_var_real3d with REAL 3d-field
    ! GENERIC   :: OPERATOR(+) => Add_real3d, Add_real3d_and_scalar, Add_real3d_and_field
    ! PROCEDURE :: Multiply_real3d_with_field       ! Multiply varialbe of type t_jsb_var_real3d with REAL 3-d field
    ! GENERIC   :: OPERATOR(*) => Multiply_real3d_with_field
    FINAL     :: finalize_jsb_var_real3d
  END TYPE t_jsb_var_real3d

!#pragma GCC diagnostic pop

  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p
  END TYPE t_jsb_var_p

!!$  INTERFACE t_jsb_var
!!$    PROCEDURE :: Create_jsb_var
!!$  END INTERFACE t_jsb_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'

CONTAINS

  SUBROUTINE finalize_jsb_var_real2d(this)

    TYPE(t_jsb_var_real2d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real2d'

    !CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr)) NULLIFY(this%ptr)
    IF (ALLOCATED(this%field)) DEALLOCATE(this%field)

  END SUBROUTINE finalize_jsb_var_real2d

  SUBROUTINE finalize_jsb_var_real3d(this)

    TYPE(t_jsb_var_real3d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real3d'

    !CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr)) NULLIFY(this%ptr)
    IF (ALLOCATED(this%field)) DEALLOCATE(this%field)

  END SUBROUTINE finalize_jsb_var_real3d

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

  ! FUNCTION Add_real2d(var1, var2) RESULT(var3)

  !   CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
  !   TYPE(t_jsb_var_real2d),  TARGET     :: var3

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field + var2%field
  !   var3%ptr => var3%field(:,:,1,1,1)

  ! END FUNCTION Add_real2d

  ! FUNCTION Add_real2d_and_scalar(var1, scalar) RESULT(var3)

  !   CLASS(t_jsb_var_real2d), INTENT(in) :: var1
  !   REAL(wp),                INTENT(in) :: scalar
  !   TYPE(t_jsb_var_real2d), TARGET              :: var3

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field + scalar
  !   var3%ptr => var3%field(:,:,1,1,1)

  ! END FUNCTION Add_real2d_and_scalar

  ! FUNCTION Add_real2d_and_field(var1, field) RESULT(var3)

  !   CLASS(t_jsb_var_real2d), INTENT(in) :: var1
  !   REAL(wp),                INTENT(in) :: field(:,:)
  !   TYPE(t_jsb_var_real2d),  TARGET     :: var3

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real2d_and_field'

  !   IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
  !     & CALL finish(TRIM(routine), 'Dimension mismatch')

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field + RESHAPE(field, (/SIZE(var1%ptr,1),SIZE(var1%ptr,2),1,1,1/))
  !   var3%ptr => var3%field(:,:,1,1,1)

  ! END FUNCTION Add_real2d_and_field

  ! FUNCTION Add_real3d(var1, var2) RESULT(var3)

  !   CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
  !   TYPE(t_jsb_var_real3d),  TARGET     :: var3

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field + var2%field
  !   var3%ptr => var3%field(:,:,:,1,1)

  ! END FUNCTION Add_real3d

  ! FUNCTION Add_real3d_and_scalar(var1, scalar) RESULT(var3)

  !   CLASS(t_jsb_var_real3d), INTENT(in) :: var1
  !   REAL(wp),                INTENT(in) :: scalar
  !   TYPE(t_jsb_var_real3d),  TARGET     :: var3

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field + scalar
  !   var3%ptr => var3%field(:,:,:,1,1)

  ! END FUNCTION Add_real3d_and_scalar

  ! FUNCTION Add_real3d_and_field(var1, field) RESULT(var3)

  !   CLASS(t_jsb_var_real3d), INTENT(in) :: var1
  !   REAL(wp),                INTENT(in) :: field(:,:,:)
  !   TYPE(t_jsb_var_real3d),  TARGET     :: var3

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real3d_and_field'

  !   IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
  !     & CALL finish(TRIM(routine), 'Dimension mismatch')

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field + RESHAPE(field, (/SIZE(var1%ptr,1),SIZE(var1%ptr,2),SIZE(var1%ptr,3),1,1/))
  !   var3%ptr => var3%field(:,:,:,1,1)

  ! END FUNCTION Add_real3d_and_field

  ! FUNCTION Multiply_real2d_with_field(var1, field) RESULT(var3)

  !   CLASS(t_jsb_var_real2d), INTENT(in) :: var1
  !   REAL(wp),                INTENT(in) :: field(:,:)
  !   TYPE(t_jsb_var_real2d),  TARGET     :: var3

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real2d_with_field'

  !   IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
  !     & CALL finish(TRIM(routine), 'Dimension mismatch')

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field * RESHAPE(field, (/SIZE(var1%ptr,1),SIZE(var1%ptr,2),1,1,1/))
  !   var3%ptr => var3%field(:,:,1,1,1)

  ! END FUNCTION Multiply_real2d_with_field

  ! FUNCTION Multiply_real3d_with_field(var1, field) RESULT(var3)

  !   CLASS(t_jsb_var_real3d), INTENT(in) :: var1
  !   REAL(wp),                INTENT(in) :: field(:,:,:)
  !   TYPE(t_jsb_var_real3d),  TARGET     :: var3

  !   CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real3d_with_field'

  !   IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
  !     & CALL finish(TRIM(routine), 'Dimension mismatch')

  !   var3%name = var1%name
  !   var3%full_name = var1%full_name
  !   var3%owner_proc_id = var1%owner_proc_id
  !   var3%owner_tile_path = var1%owner_tile_path
  !   var3%missval = var1%missval
  !   var3%field = var1%field + RESHAPE(field, (/SIZE(var1%ptr,1),SIZE(var1%ptr,2),SIZE(var1%ptr,3),1,1/))
  !   var3%ptr => var3%field(:,:,:,1,1)

  ! END FUNCTION Multiply_real3d_with_field

#endif
END MODULE mo_jsb_var_class
