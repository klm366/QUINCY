!>
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_jsb_var_class
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp, dp
  USE mo_exception,         ONLY: finish, message
  USE mo_hermetic,          ONLY: t_hermetic
  USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN
  USE mo_jsb_subset,        ONLY: t_subset, ON_DOMAIN, ON_CHUNK

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var, t_jsb_var_p, REAL2D, REAL3D

  ENUM, BIND(C)
    ENUMERATOR :: REAL2D=1, REAL3D
  END ENUM

  !================================================================================================================================
  ! t_jsb_var
  !================================================================================================================================
  TYPE, ABSTRACT :: t_jsb_var
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr2d(:,:) => NULL()
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
      , CONTIGUOUS &
#endif
      :: ptr3d(:,:,:) => NULL()
    CHARACTER(len=VARNAME_LEN)    :: name                      ! Name of variable
    CHARACTER(len=VARNAME_LEN)    :: full_name                 ! Full name of variable in varlist, including prefix and suffix
    INTEGER                       :: type                      ! REAL2D / REAL3D
    INTEGER                       :: subset_type               ! ON_DOMAIN / ON_CHUNK
    CHARACTER(len=:), ALLOCATABLE :: element_name
    CHARACTER(len=:), ALLOCATABLE :: element_shortname         ! Short name of element
    INTEGER                       :: element_id
    INTEGER                       :: owner_model_id
    INTEGER                       :: owner_proc_id = -1
    INTEGER, ALLOCATABLE          :: owner_tile_path(:)
    INTEGER, ALLOCATABLE          :: child_idx(:)              ! Position-numbers of same variable for for the child tiles 
                                                               ! (according to mem%vars(:)).
    ! TYPE(t_subset), POINTER     :: sub
    REAL(dp)                      :: missval
    LOGICAL                       :: l_aggregate_all = .FALSE. ! When aggregating, consider fractions of all child tiles even
                                                               ! of those where variable doesn't exist
    LOGICAL                       :: is_in_output  = .FALSE.   ! Whether variable is requested to be written to output
    LOGICAL                       :: is_in_restart = .FALSE.   ! Whether variable is requested to be written to restart file
  CONTAINS
    PROCEDURE :: Get_subset         => t_jsb_var_get_subset
    PROCEDURE :: force_finalization => t_jsb_var_force_finalization
  END TYPE t_jsb_var

!#pragma GCC diagnostic push
!#pragma GCC diagnostic ignored "-Wsurprising"

  !================================================================================================================================
  ! t_jsb_var_real2d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), POINTER &
#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
     , CONTIGUOUS &
#endif
      :: ptr(:,:) => NULL(), &
      &  ptrm1(:) => NULL()
  CONTAINS
    PROCEDURE :: Associate_pointers => t_jsb_var_real2d_associate_pointers
    PROCEDURE :: Add_real2d                       ! Add two variables of type t_jsb_var_real2d
    PROCEDURE :: Add_real2d_and_scalar            ! Add variable of type t_jsb_var_real2d with REAL scalar
    PROCEDURE :: Add_real2d_and_field             ! Add variable of type t_jsb_var_real2d with REAL 2-d field
    GENERIC   :: OPERATOR(+) => Add_real2d, Add_real2d_and_scalar, Add_real2d_and_field
    PROCEDURE :: Multiply_real2d_with_field       ! Multiply varialbe of type t_jsb_var_real2d with REAL 2-d field
    GENERIC   :: OPERATOR(*) => Multiply_real2d_with_field
    PROCEDURE :: Assign_var_to_var_real2d
    PROCEDURE :: Assign_field1d_to_var_real2d
    PROCEDURE :: Assign_field2d_to_var_real2d
    PROCEDURE :: Assign_scalar_to_var_real2d
    GENERIC   :: ASSIGNMENT(=) => Assign_var_to_var_real2d, &
                                  Assign_field1d_to_var_real2d, Assign_field2d_to_var_real2d, &
      &                           Assign_scalar_to_var_real2d
    FINAL     :: finalize_jsb_var_real2d
  END TYPE t_jsb_var_real2d

  !================================================================================================================================
  ! t_jsb_var_real3d
  !================================================================================================================================
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(wp), POINTER :: ptr(:,:,:) => NULL()
    REAL(wp), POINTER :: ptrm1(:,:) => NULL()
  CONTAINS
    PROCEDURE :: Associate_pointers => t_jsb_var_real3d_associate_pointers
    PROCEDURE :: Add_real3d                       ! Add two variables of type t_jsb_var_real3d
    PROCEDURE :: Add_real3d_and_scalar            ! Add variable of type t_jsb_var_real3d with REAL scalar
    PROCEDURE :: Add_real3d_and_field             ! Add variable of type t_jsb_var_real3d with REAL 3d-field
    GENERIC   :: OPERATOR(+) => Add_real3d, Add_real3d_and_scalar, Add_real3d_and_field
    PROCEDURE :: Multiply_real3d_with_field       ! Multiply varialbe of type t_jsb_var_real3d with REAL 3-d field
    GENERIC   :: OPERATOR(*) => Multiply_real3d_with_field
    PROCEDURE :: Assign_var_to_var_real3d         ! %ptr3d = %ptr3d
    PROCEDURE :: Assign_field2d_to_var_real3d     ! %ptr3d = real3d%ptrm1   ! ON_CHUNK
    PROCEDURE :: Assign_field3d_to_var_real3d     ! %ptr3d = (:,:,:)        ! ON_DOMAIN
    PROCEDURE :: Assign_scalar_to_var_real3d      ! %ptr3d = scalar
    GENERIC   :: ASSIGNMENT(=) => Assign_var_to_var_real3d, &
                                  Assign_field2d_to_var_real3d, Assign_field3d_to_var_real3d, &
      &                           Assign_scalar_to_var_real3d
    FINAL     :: finalize_jsb_var_real3d
  END TYPE t_jsb_var_real3d

!#pragma GCC diagnostic pop

  !================================================================================================================================
  ! t_jsb_var_p
  !================================================================================================================================
  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p => NULL()
  END TYPE t_jsb_var_p

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'

CONTAINS

  FUNCTION t_jsb_var_get_subset(this) RESULT(subset)

    USE mo_jsb_subset,   ONLY: jsbach_subsets
    USE mo_jsb_parallel, ONLY: Get_omp_thread

    CLASS(t_jsb_var), INTENT(in) :: this
    TYPE(t_subset)               :: subset

    ! print*,'FFF ',TRIM(this%full_name), this%owner_model_id, Get_omp_thread()
    subset = jsbach_subsets(this%owner_model_id)%sub(Get_omp_thread())
    ! print*,'GGG ',subset%type,subset%ics,subset%ice,subset%iblk

  END FUNCTION t_jsb_var_get_subset

  SUBROUTINE t_jsb_var_real2d_associate_pointers(this,ic_start,ic_end,iblk_start,iblk_end)

    CLASS(t_jsb_var_real2d), INTENT(inout), TARGET :: this
    INTEGER,                 INTENT(in)    :: ic_start
    INTEGER,                 INTENT(in)    :: ic_end
    INTEGER,                 INTENT(in)    :: iblk_start
    INTEGER,                 INTENT(in)    :: iblk_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_var_real2d_associate_pointers'

    IF (ASSOCIATED(this%ptrm1)) NULLIFY(this%ptrm1)
    IF (iblk_start == iblk_end) THEN
      this%ptrm1 => this%ptr2d(ic_start:ic_end, iblk_start)
    END IF

  END SUBROUTINE t_jsb_var_real2d_associate_pointers

  SUBROUTINE t_jsb_var_real3d_associate_pointers(this,ic_start,ic_end,iblk_start,iblk_end)

    CLASS(t_jsb_var_real3d), INTENT(inout), TARGET :: this
    INTEGER,                 INTENT(in)    :: ic_start
    INTEGER,                 INTENT(in)    :: ic_end
    INTEGER,                 INTENT(in)    :: iblk_start
    INTEGER,                 INTENT(in)    :: iblk_end

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_var_real3d_associate_pointers'

    IF (ASSOCIATED(this%ptrm1)) NULLIFY(this%ptrm1)
    IF (iblk_start == iblk_end) THEN
      this%ptrm1 => this%ptr3d(ic_start:ic_end, :, iblk_start)
    END IF

  END SUBROUTINE t_jsb_var_real3d_associate_pointers

  !================================================================================================================================
  ! force finalization
  !================================================================================================================================
  SUBROUTINE t_jsb_var_force_finalization(this)

    CLASS(t_jsb_var), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':t_jsb_var_force_finalization'

  END SUBROUTINE t_jsb_var_force_finalization

  SUBROUTINE finalize_jsb_var_real2d(this)

    TYPE(t_jsb_var_real2d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real2d'

    CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr)) DEALLOCATE(this%ptr)

  END SUBROUTINE finalize_jsb_var_real2d

  SUBROUTINE finalize_jsb_var_real3d(this)

    TYPE(t_jsb_var_real3d) :: this

    CHARACTER(len=*), PARAMETER :: routine = modname//':finalize_jsb_var_real3d'

    !CALL message(routine, 'Finalizing '//TRIM(this%name))

    IF (ALLOCATED(this%owner_tile_path)) DEALLOCATE(this%owner_tile_path)
    IF (ALLOCATED(this%child_idx)) DEALLOCATE(this%child_idx)

    IF (ASSOCIATED(this%ptr)) DEALLOCATE(this%ptr)

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

  !================================================================================================================================
  ! assign
  !================================================================================================================================
  SUBROUTINE Assign_var_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    CLASS(t_jsb_var_real2d), INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_var_to_var_real2d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      lhs%ptr2d(ics:ice,iblk) = rhs%ptr2d(ics:ice,iblk)
    ELSE
      lhs%ptr2d(:,:) = rhs%ptr2d(:,:)
    END IF

  END SUBROUTINE Assign_var_to_var_real2d

  SUBROUTINE Assign_var_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    CLASS(t_jsb_var_real3d), INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_var_to_var_real3d'

    sub = lhs%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      lhs%ptr3d(ics:ice,:,iblk) = rhs%ptr3d(ics:ice,:,iblk)
    ELSE
      lhs%ptr3d(:,:,:) = rhs%ptr3d(:,:,:)
    END IF

  END SUBROUTINE Assign_var_to_var_real3d

  SUBROUTINE Assign_field1d_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs(:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field1d_to_var_real2d'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc
      lhs%ptr2d(ics:ice,iblk) = rhs(:)
    ELSE
      CALL finish(routine, 'Cannot assign 1d-field to jsb_var on domain')
    END IF

  END SUBROUTINE Assign_field1d_to_var_real2d

  SUBROUTINE Assign_field2d_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    REAL(wp), INTENT(in)    :: rhs(:,:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_var_real2d'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,2) > 1) CALL finish(routine, 'Dimension mismatch for assignment to jsb_var on chunk')
      ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc
      lhs%ptr2d(ics:ice,iblk) = rhs(:,1)
    ELSE
      lhs%ptr2d(:,:) = rhs(:,:)
    END IF

  END SUBROUTINE Assign_field2d_to_var_real2d

  SUBROUTINE Assign_scalar_to_var_real2d(lhs, rhs)

    CLASS(t_jsb_var_real2d), INTENT(inout) :: lhs
    REAL(wp), INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_var_real2d'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc
      lhs%ptr2d(ics:ice,iblk) = rhs
    ELSE
      lhs%ptr2d(:,:) = rhs
    END IF

  END SUBROUTINE Assign_scalar_to_var_real2d

  SUBROUTINE Assign_field2d_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    REAL(wp),                INTENT(in)    :: rhs(:,:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field2d_to_var_real3d'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc
      lhs%ptr3d(ics:ice,:,iblk) = rhs(:,:)
    ELSE
      CALL finish(routine, 'Cannot assign 2d-field to jsb_var on domain')
    END IF

  END SUBROUTINE Assign_field2d_to_var_real3d

  SUBROUTINE Assign_field3d_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    REAL(wp), INTENT(in)    :: rhs(:,:,:)

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_field3d_to_var_real3d'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      IF (SIZE(rhs,3) > 1) CALL finish(routine, 'Dimension mismatch for assignment to jsb_var on chunk')
      ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc
      lhs%ptr3d(ics:ice,:,iblk) = rhs(:,:,1)
    ELSE
      lhs%ptr3d(:,:,:) = rhs(:,:,:)
    END IF

  END SUBROUTINE Assign_field3d_to_var_real3d

  SUBROUTINE Assign_scalar_to_var_real3d(lhs, rhs)

    CLASS(t_jsb_var_real3d), INTENT(inout) :: lhs
    REAL(wp), INTENT(in)    :: rhs

    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Assign_scalar_to_var_real3d'

    sub = lhs%Get_subset()

    IF (sub%type == ON_CHUNK) THEN
      ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc
      lhs%ptr3d(ics:ice,:,iblk) = rhs
    ELSE
      lhs%ptr3d(:,:,:) = rhs
    END IF

  END SUBROUTINE Assign_scalar_to_var_real3d

  !================================================================================================================================
  ! add
  !================================================================================================================================
  FUNCTION Add_real2d(var1, var2) RESULT(res)
  
    CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
    REAL(wp), ALLOCATABLE :: res(:,:)
  
    TYPE(t_subset) :: sub
    INTEGER        :: ics, ice, iblk, nb, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real2d'

    sub = var1%Get_subset()
    ics = sub%ics ; ice = sub%ice ; iblk = sub%iblk ; nb = sub%nb ; nc = sub%nc

    IF (sub%type == ON_CHUNK) THEN
      ALLOCATE(res(nc,iblk:iblk))
      res(:,iblk) = var1%ptr(ics:ice,iblk) + var2%ptr(ics:ice,iblk)
    ELSE
      ALLOCATE(res(nc,nb))
      res(:,:) = var1%ptr(:,:) + var2%ptr(:,:)
    END IF
  
  END FUNCTION Add_real2d

  ! FUNCTION Add_real2d(var1, var2) RESULT(res)
  ! 
  !   CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
  !   REAL(wp) :: res(SIZE(var1%ptr,1),SIZE(var1%ptr,2))
  ! 
  !   res(:,:) = var1%ptr(:,:) + var2%ptr(:,:)
  ! 
  ! END FUNCTION Add_real2d

  FUNCTION Add_real3d(var1, var2) RESULT(res)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
    REAL(wp) :: res(SIZE(var1%ptr,1),SIZE(var1%ptr,2),SIZE(var1%ptr,3))

    res(:,:,:) = var1%ptr(:,:,:) + var2%ptr(:,:,:)

  END FUNCTION Add_real3d

  FUNCTION Add_real2d_and_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), TARGET              :: var3

    var3%name = var1%name
    var3%full_name = var1%full_name
    var3%owner_proc_id = var1%owner_proc_id
    var3%owner_tile_path = var1%owner_tile_path
    var3%missval = var1%missval
    var3%ptr2d = var1%ptr2d + scalar
    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_real2d_and_scalar

  FUNCTION Add_real2d_and_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real2d_and_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%name = var1%name
    var3%full_name = var1%full_name
    var3%owner_proc_id = var1%owner_proc_id
    var3%owner_tile_path = var1%owner_tile_path
    var3%missval = var1%missval
    var3%ptr2d = var1%ptr2d + field
    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Add_real2d_and_field

  FUNCTION Add_real3d_and_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    var3%name = var1%name
    var3%full_name = var1%full_name
    var3%owner_proc_id = var1%owner_proc_id
    var3%owner_tile_path = var1%owner_tile_path
    var3%missval = var1%missval
    var3%ptr3d = var1%ptr3d + scalar
    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_real3d_and_scalar

  FUNCTION Add_real3d_and_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_real3d_and_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%name = var1%name
    var3%full_name = var1%full_name
    var3%owner_proc_id = var1%owner_proc_id
    var3%owner_tile_path = var1%owner_tile_path
    var3%missval = var1%missval
    var3%ptr3d = var1%ptr3d + field
    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Add_real3d_and_field

  !================================================================================================================================
  ! multiply
  !================================================================================================================================
  FUNCTION Multiply_real2d_with_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d),  TARGET     :: var3

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real2d_with_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%name = var1%name
    var3%full_name = var1%full_name
    var3%owner_proc_id = var1%owner_proc_id
    var3%owner_tile_path = var1%owner_tile_path
    var3%missval = var1%missval
    var3%ptr2d = var1%ptr2d * field
    var3%ptr => var3%ptr2d(:,:)

  END FUNCTION Multiply_real2d_with_field

  FUNCTION Multiply_real3d_with_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d),  TARGET     :: var3

    CHARACTER(len=*), PARAMETER :: routine = modname//':Multiply_real3d_with_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%name = var1%name
    var3%full_name = var1%full_name
    var3%owner_proc_id = var1%owner_proc_id
    var3%owner_tile_path = var1%owner_tile_path
    var3%missval = var1%missval
    var3%ptr3d = var1%ptr3d + field
    var3%ptr => var3%ptr3d(:,:,:)

  END FUNCTION Multiply_real3d_with_field

#endif
END MODULE mo_jsb_var_class
