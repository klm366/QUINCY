!>
!! @brief Contains the definition of variable types used across JSBACH4
!!
!! @author
!!  Jan Engel, MPI-BGC Jena (copied from JSBACH4 git at 20160815)
!!
!! @par Revision History
!! First version by Jan Engel (20160815); no changes compared to JSBACH4 except for not USEing mo_jsb_varlist_iface
!!
MODULE mo_jsb_var_class
!#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp, dp
  USE mo_exception, ONLY: finish

!>  @todo currently commented out and defined locally VARNAME_LEN \n
!!  USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN   ! *slmtemp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_var_real2d, t_jsb_var_real3d, t_jsb_var, t_jsb_var_p
  PUBLIC :: t_jsb_struct, t_jsb_struct_p, t_jsb_n_struct, t_jsb_n_struct_p    ! *slmjsb4
  PUBLIC :: t_jsb_var_integer2d, t_jsb_var_logical2d                          ! *slmjsb4
  PUBLIC :: t_jsb_var_real1d, t_jsb_var_integer1d, t_jsb_var_logical1d  !! for SLM compatbility -> used in veg_mem_class  *slmtemp

!> @todo in JSBACH4 this parameter is defined in mo_jsb_varlist_iface (USE mo_jsb_varlist_iface, ONLY: VARNAME_LEN), 
!! thus, delete here and insert the USE statement instead later
!! *slmtemp
  INTEGER, PARAMETER :: VARNAME_LEN = 64

  ENUM, BIND(C)
    ENUMERATOR :: REAL2D=1, REAL3D
  END ENUM

  !-----------------------------------------------------------------------------------------------------
  !> Abstract Variable TYPE - t_jsb_var
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_var
!!$    INTEGER           :: type = REAL2D
    CHARACTER(len=VARNAME_LEN) :: name            ! Name of variable
    CHARACTER(len=VARNAME_LEN) :: full_name       ! Full name of variable in varlist, including prefix and suffix
    INTEGER                    :: owner_proc_id = -1
    INTEGER, ALLOCATABLE       :: owner_tile_path(:)
    REAL(dp)                   :: missval
  CONTAINS
  END TYPE t_jsb_var

  !-----------------------------------------------------------------------------------------------------!
  ! 1D  *slmtemp                                                                                        !
  !-----------------------------------------------------------------------------------------------------!

! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> 1D LOGICAL Variable TYPE - t_jsb_var_logical1d
  !! 
  !! added for temporary implementation of 1D REAL into SLM -> they will be changed to 2D later
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_logical1d
    LOGICAL, POINTER :: ptr => NULL()
  CONTAINS
  END TYPE t_jsb_var_logical1d

! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> 1D INTEGER Variable TYPE - t_jsb_var_integer1d
  !! 
  !! added for temporary implementation of 1D REAL into SLM -> they will be changed to 2D later
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_integer1d
    INTEGER, POINTER :: ptr => NULL()
  CONTAINS
  END TYPE t_jsb_var_integer1d

! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> 1D REAL Variable TYPE - t_jsb_var_real1d
  !! 
  !! added for temporary implementation of 1D REAL into SLM -> they will be changed to 2D later
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real1d
    REAL(wp), POINTER :: ptr => NULL()
  CONTAINS
!    PROCEDURE :: Add_real1d                       ! Add two variables of type t_jsb_var_real2d
!    PROCEDURE :: Add_real1d_and_scalar            ! Add varialbe of type t_jsb_var_real2d with REAL scalar
!    GENERIC   :: OPERATOR(+) => Add_real1d, Add_real1d_and_scalar
  END TYPE t_jsb_var_real1d

  !-----------------------------------------------------------------------------------------------------!
  ! 2D                                                                                                  !
  !-----------------------------------------------------------------------------------------------------!

! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> 2D LOGICAL Variable TYPE - t_jsb_var_logical2d
  !! 
  !! added for temporary implementation of 1D REAL into SLM -> they will be changed to 2D later
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_logical2d
    LOGICAL, POINTER :: ptr(:,:) => NULL()
  CONTAINS
  END TYPE t_jsb_var_logical2d

! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> 2D INTEGER Variable TYPE - t_jsb_var_integer2d
  !! 
  !! added for temporary implementation of 1D REAL into SLM -> they will be changed to 2D later
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_integer2d
    INTEGER, POINTER :: ptr(:,:) => NULL()
  CONTAINS
  END TYPE t_jsb_var_integer2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D REAL Variable TYPE
  !! 
  !! ics:ice,iblk
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), POINTER :: ptr(:,:) => NULL()
  CONTAINS
    PROCEDURE :: Add_real2d                       ! Add two variables of type t_jsb_var_real2d
    PROCEDURE :: Add_real2d_and_scalar            ! Add varialbe of type t_jsb_var_real2d with REAL scalar
    PROCEDURE :: Add_real2d_and_field             ! Add varialbe of type t_jsb_var_real2d with REAL 2-d field
    GENERIC   :: OPERATOR(+) => Add_real2d, Add_real2d_and_scalar, Add_real2d_and_field
    PROCEDURE :: Multiply_real2d_with_field       ! Multiply variable of type t_jsb_var_real2d with REAL 2-d field
    GENERIC   :: OPERATOR(*) => Multiply_real2d_with_field
  END TYPE t_jsb_var_real2d

  !-----------------------------------------------------------------------------------------------------!
  ! 3D                                                                                                  !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> 3D REAL Variable TYPE
  !! 
  !! ics:ice,:,iblk
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real3d
    REAL(wp), POINTER :: ptr(:,:,:) => NULL()
  CONTAINS
    PROCEDURE :: Add_real3d
    PROCEDURE :: Add_real3d_and_scalar            ! Add varialbe of type t_jsb_var_real3d with REAL scalar
    PROCEDURE :: Add_real3d_and_field             ! Add varialbe of type t_jsb_var_real3d with REAL 3d-field
    GENERIC   :: OPERATOR(+) => Add_real3d, Add_real3d_and_scalar, Add_real3d_and_field
    PROCEDURE :: Multiply_real3d_with_field       ! Multiply variable of type t_jsb_var_real3d with REAL 3-d field
    GENERIC   :: OPERATOR(*) => Multiply_real3d_with_field
  END TYPE t_jsb_var_real3d

  !-----------------------------------------------------------------------------------------------------
  !> Pointer to the abstract Variable TYPE
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_var_p
    CLASS(t_jsb_var), POINTER :: p
  END TYPE t_jsb_var_p

! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> TYPE that contains multiple variables
  !! 
  !! this type is aimed to be either included in: \n 
  !!    (1) the process_memory directly (mem\%TYPE) or \n
  !!    (2) a TYPE t_jsb_n_struct
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_struct
    TYPE(t_jsb_var_p), ALLOCATABLE :: vars(:)
    CHARACTER(len=VARNAME_LEN) :: name                !< Name of variable
    CHARACTER(len=VARNAME_LEN) :: full_name           !< Full name of variable in varlist, including prefix and suffix
    INTEGER                    :: owner_proc_id = -1  !< owner process id
    INTEGER, ALLOCATABLE       :: owner_tile_path(:)  !< owner tile path
    INTEGER                    :: no_of_vars = 0      !< number of variables the TYPE contains in vars(:)
    INTEGER                    :: max_no_of_vars = 20 !< ! hard coded limit of max var number - may change later !
    CHARACTER(LEN=50)          :: varlist_name        !< name of variable list
  CONTAINS
  END TYPE t_jsb_struct

! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> POINTER to the t_jsb_struct
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_struct_p
    CLASS(t_jsb_struct), POINTER :: p
  END TYPE t_jsb_struct_p

! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> TYPE that contains multiple TYPEs of t_jsb_struct
  !! 
  !! this type is aimed to be included in the process_memory directly (mem\%TYPE)
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, ABSTRACT :: t_jsb_n_struct
    TYPE(t_jsb_struct_p), ALLOCATABLE :: l2_structs(:)
    CHARACTER(len=VARNAME_LEN) :: name            ! Name of variable
    CHARACTER(len=VARNAME_LEN) :: full_name       ! Full name of variable in varlist, including prefix and suffix
    INTEGER                    :: owner_proc_id = -1
    INTEGER, ALLOCATABLE       :: owner_tile_path(:)
    INTEGER                    :: no_of_l2_structs = 0
    INTEGER                    :: max_no_of_l2_structs = 20 !!! hard coded limit - may change later !!!
    CHARACTER(LEN=50)          :: l2_structlist_name
  CONTAINS
  END TYPE t_jsb_n_struct

! *slmjsb4
  !-----------------------------------------------------------------------------------------------------
  !> POINTER to the t_jsb_n_struct
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_jsb_n_struct_p
    CLASS(t_jsb_n_struct), POINTER :: p
  END TYPE t_jsb_n_struct_p

!!$  INTERFACE t_jsb_var
!!$    PROCEDURE :: Create_jsb_var
!!$  END INTERFACE t_jsb_var

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_var_class'

CONTAINS

!! *slmtemp
!  !-----------------------------------------------------------------------------------------------------
!  !> operator overload: addition of two 1D variables 1D + 1D
!  !! 
!  !! added for temporary implementation of 1D REAL into SLM
!  !!
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_real1d(this, var2) RESULT(var3)

!    CLASS(t_jsb_var_real1d), INTENT(in) :: this, var2
!    TYPE(t_jsb_var_real1d)              :: var3

!    ALLOCATE(var3%ptr)    ! *slmdev needed with current Intel compiler too? - what about deallocation?
!    var3%ptr = this%ptr + var2%ptr
!  END FUNCTION Add_real1d

!! *slmtemp
!  !-----------------------------------------------------------------------------------------------------
!  !> operator overload: addition of two 1D variables 1D + scalar
!  !! 
!  !! added for temporary implementation of 1D REAL into SLM
!  !! 
!  !-----------------------------------------------------------------------------------------------------
!  FUNCTION Add_real1d_and_scalar(this, scalar) RESULT(var3)

!    CLASS(t_jsb_var_real1d), INTENT(in) :: this
!    REAL(wp),                INTENT(in) :: scalar
!    TYPE(t_jsb_var_real1d)              :: var3

!    ALLOCATE(var3%ptr)    ! *slmdev needed with current Intel compiler too? - what about deallocation?
!    var3%ptr = this%ptr + scalar
!  END FUNCTION Add_real1d_and_scalar


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

  FUNCTION Add_real2d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real2d)              :: var3

    var3%ptr(:,:) = var1%ptr(:,:) + var2%ptr(:,:)

  END FUNCTION Add_real2d

  FUNCTION Add_real2d_and_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d)              :: var3

    var3%ptr(:,:) = var1%ptr(:,:) + scalar

  END FUNCTION Add_real2d_and_scalar

  FUNCTION Add_real2d_and_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d)              :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_real2d_and_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%ptr(:,:) = var1%ptr(:,:) + field(:,:)

  END FUNCTION Add_real2d_and_field

  FUNCTION Add_real3d(var1, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
    TYPE(t_jsb_var_real3d)              :: var3

    var3%ptr(:,:,:) = var1%ptr(:,:,:) + var2%ptr(:,:,:)

  END FUNCTION Add_real3d

  FUNCTION Add_real3d_and_scalar(var1, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d)              :: var3

    var3%ptr(:,:,:) = var1%ptr(:,:,:) + scalar

  END FUNCTION Add_real3d_and_scalar

  FUNCTION Add_real3d_and_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d)              :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_real3d_and_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%ptr(:,:,:) = var1%ptr(:,:,:) + field(:,:,:)

  END FUNCTION Add_real3d_and_field

  FUNCTION Multiply_real2d_with_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d)              :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_real2d_with_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%ptr(:,:) = var1%ptr(:,:) * field(:,:)

  END FUNCTION Multiply_real2d_with_field

  FUNCTION Multiply_real3d_with_field(var1, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d)              :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_real3d_with_field'

    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
      & CALL finish(TRIM(routine), 'Dimension mismatch')

    var3%ptr(:,:,:) = var1%ptr(:,:,:) + field(:,:,:)

  END FUNCTION Multiply_real3d_with_field

!#endif
END MODULE mo_jsb_var_class
