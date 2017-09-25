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
  PUBLIC :: t_jsb_var_logical2d                                               ! *slmjsb4
  PUBLIC :: t_jsb_var_real1d, t_jsb_var_logical1d                             !! for SLM compatbility -> used in veg_mem_class  *slmtemp
  PUBLIC :: REAL2D, REAL3D

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
    CHARACTER(len=VARNAME_LEN) :: name            ! Name of variable
    CHARACTER(len=VARNAME_LEN) :: full_name       ! Full name of variable in varlist, including prefix and suffix
    INTEGER                    :: type            ! REAL2D / REAL3D
    INTEGER                    :: owner_proc_id = -1
    INTEGER, ALLOCATABLE       :: owner_tile_path(:)
    INTEGER, ALLOCATABLE       :: child_idx(:)
    REAL(dp)                   :: missval
    LOGICAL                    :: l_aggregate_all = .FALSE. ! When aggregating, consider fractions of all child tiles even
                                                  ! of those where variable doesn't exist
  END TYPE t_jsb_var

  !-----------------------------------------------------------------------------------------------------!
  ! 1D/Scalar  *slmtemp                                                                                        !
  !-----------------------------------------------------------------------------------------------------!

! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> 1D LOGICAL Variable TYPE - t_jsb_var_logical1d
  !! 
  !! added for temporary implementation of 1D REAL into SLM -> they will be changed to 2D later
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_logical1d
    LOGICAL, ALLOCATABLE :: ptr
  !CONTAINS
  END TYPE t_jsb_var_logical1d

! *slmtemp
  !-----------------------------------------------------------------------------------------------------
  !> 1D REAL Variable TYPE - t_jsb_var_real1d
  !! 
  !! added for temporary implementation of 1D REAL into SLM -> they will be changed to 2D later
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real1d
    REAL(wp), ALLOCATABLE :: ptr
  CONTAINS
    PRIVATE
    PROCEDURE               :: Add_real1d                       ! Add two variables of type t_jsb_var_real1d
    PROCEDURE               :: Add_real1d_and_scalar            ! Add variable of type t_jsb_var_real1d with REAL scalar
    PROCEDURE, PASS(this2)  :: Add_scalar_and_real1d            ! Add variable of type t_jsb_var_real1d with REAL scalar
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_real1d, Add_real1d_and_scalar, Add_scalar_and_real1d
    PROCEDURE               :: Subtract_real1d                  ! Subtract two variables of type t_jsb_var_real1d
    PROCEDURE               :: Subtract_real1d_and_scalar       ! Subtract variable of type t_jsb_var_real1d with REAL scalar
    PROCEDURE, PASS(this2)  :: Subtract_scalar_and_real1d       ! Subtract variable of type t_jsb_var_real1d with REAL scalar
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_real1d, Subtract_real1d_and_scalar, Subtract_scalar_and_real1d
    PROCEDURE               :: Multiply_real1d                  ! Multiply two variables of type t_jsb_var_real1d
    PROCEDURE               :: Multiply_real1d_and_scalar       ! Multiply variable of type t_jsb_var_real1d with REAL scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_real1d       ! Multiply variable of type t_jsb_var_real1d with REAL scalar
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_real1d, Multiply_real1d_and_scalar, Multiply_scalar_and_real1d
    PROCEDURE               :: Divide_real1d                    ! Divide two variables of type t_jsb_var_real1d
    PROCEDURE               :: Divide_real1d_and_scalar         ! Divide variable of type t_jsb_var_real1d by REAL scalar
    PROCEDURE, PASS(this2)  :: Divide_scalar_and_real1d         ! Divide variable of type REAL scalar by t_jsb_var_real1d
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_real1d, Divide_real1d_and_scalar, Divide_scalar_and_real1d
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
    LOGICAL, ALLOCATABLE :: ptr(:,:)
  CONTAINS
  END TYPE t_jsb_var_logical2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D REAL Variable TYPE
  !! 
  !! ics:ice,iblk
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_var) :: t_jsb_var_real2d
    REAL(wp), ALLOCATABLE :: ptr(:,:)
  CONTAINS
    PRIVATE
    PROCEDURE               :: Add_real2d                       ! Add two variables of type t_jsb_var_real2d
    PROCEDURE               :: Add_real2d_and_scalar            ! Add variable of type t_jsb_var_real2d with REAL scalar
    PROCEDURE, PASS(this2)  :: Add_scalar_and_real2d            ! Add variable of type t_jsb_var_real2d with REAL scalar
    PROCEDURE               :: Add_real2d_and_field             ! Add variable of type t_jsb_var_real2d with REAL 2-d field
    PROCEDURE, PASS(this2)  :: Add_field_and_real2d             ! Add variable of REAL 2-d field with type t_jsb_var_real2d
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_real2d, Add_real2d_and_scalar, Add_scalar_and_real2d, &
                                              Add_real2d_and_field, Add_field_and_real2d
    PROCEDURE               :: Subtract_real2d                  ! Subtract two variables of type t_jsb_var_real2d
    PROCEDURE               :: Subtract_real2d_and_scalar       ! Subtract variable of type t_jsb_var_real2d with REAL scalar
    PROCEDURE, PASS(this2)  :: Subtract_scalar_and_real2d       ! Subtract variable of type t_jsb_var_real2d with REAL scalar
    PROCEDURE               :: Subtract_real2d_and_field        ! Subtract variable of type t_jsb_var_real2d with REAL 2-d field
    PROCEDURE, PASS(this2)  :: Subtract_field_and_real2d        ! Subtract variable of REAL 2-d field with type t_jsb_var_real2d
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_real2d, Subtract_real2d_and_scalar, Subtract_scalar_and_real2d, &
                                              Subtract_real2d_and_field, Subtract_field_and_real2d
    PROCEDURE               :: Multiply_real2d                  ! Multiply two variables of type t_jsb_var_real2d
    PROCEDURE               :: Multiply_real2d_and_scalar       ! Multiply variable of type t_jsb_var_real2d with REAL scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_real2d       ! Multiply variable of type t_jsb_var_real2d with REAL scalar
    PROCEDURE               :: Multiply_real2d_and_field        ! Multiply variable of type t_jsb_var_real2d with REAL 2-d field
    PROCEDURE, PASS(this2)  :: Multiply_field_and_real2d        ! Multiply variable of REAL 2-d field with type t_jsb_var_real2d
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_real2d, Multiply_real2d_and_scalar, Multiply_scalar_and_real2d, &
                                              Multiply_real2d_and_field, Multiply_field_and_real2d
    PROCEDURE               :: Divide_real2d                    ! Divide two variables of type t_jsb_var_real2d
    PROCEDURE               :: Divide_real2d_and_scalar         ! Divide variable of type t_jsb_var_real2d by REAL scalar
    PROCEDURE               :: Divide_real2d_and_field          ! Divide variable of type t_jsb_var_real2d by REAL 2-d field
    PROCEDURE, PASS(this2)  :: Divide_field_and_real2d          ! Divide variable of REAL 2-d field by type t_jsb_var_real2d
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_real2d, Divide_real2d_and_scalar, &
                                              Divide_real2d_and_field, Divide_field_and_real2d
    PROCEDURE, PUBLIC       :: Init => init_t_jsb_var_real2d
    !FINAL ...
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
    REAL(wp), ALLOCATABLE :: ptr(:,:,:)
  CONTAINS
    PRIVATE
    PROCEDURE               :: Add_real3d                       ! Add two variables of type t_jsb_var_real3d
    PROCEDURE               :: Add_real3d_and_scalar            ! Add variable of type t_jsb_var_real3d with REAL scalar
    PROCEDURE, PASS(this2)  :: Add_scalar_and_real3d            ! Add variable of type t_jsb_var_real3d with REAL scalar
    PROCEDURE               :: Add_real3d_and_field             ! Add variable of type t_jsb_var_real3d with REAL 2-d field
    PROCEDURE, PASS(this2)  :: Add_field_and_real3d             ! Add variable of REAL 2-d field with type t_jsb_var_real3d
    GENERIC, PUBLIC         :: OPERATOR(+) => Add_real3d, Add_real3d_and_scalar, Add_scalar_and_real3d, &
                                              Add_real3d_and_field, Add_field_and_real3d
    PROCEDURE               :: Subtract_real3d                  ! Subtract two variables of type t_jsb_var_real3d
    PROCEDURE               :: Subtract_real3d_and_scalar       ! Subtract variable of type t_jsb_var_real3d with REAL scalar
    PROCEDURE, PASS(this2)  :: Subtract_scalar_and_real3d       ! Subtract variable of type t_jsb_var_real3d with REAL scalar
    PROCEDURE               :: Subtract_real3d_and_field        ! Subtract variable of type t_jsb_var_real3d with REAL 2-d field
    PROCEDURE, PASS(this2)  :: Subtract_field_and_real3d        ! Subtract variable of REAL 2-d field with type t_jsb_var_real3d
    GENERIC, PUBLIC         :: OPERATOR(-) => Subtract_real3d, Subtract_real3d_and_scalar, Subtract_scalar_and_real3d, &
                                              Subtract_real3d_and_field, Subtract_field_and_real3d
    PROCEDURE               :: Multiply_real3d                  ! Multiply two variables of type t_jsb_var_real3d
    PROCEDURE               :: Multiply_real3d_and_scalar       ! Multiply variable of type t_jsb_var_real3d with REAL scalar
    PROCEDURE, PASS(this2)  :: Multiply_scalar_and_real3d       ! Multiply variable of type t_jsb_var_real3d with REAL scalar
    PROCEDURE               :: Multiply_real3d_and_field        ! Multiply variable of type t_jsb_var_real3d with REAL 2-d field
    PROCEDURE, PASS(this2)  :: Multiply_field_and_real3d        ! Multiply variable of REAL 2-d field with type t_jsb_var_real3d
    GENERIC, PUBLIC         :: OPERATOR(*) => Multiply_real3d, Multiply_real3d_and_scalar, Multiply_scalar_and_real3d, &
                                              Multiply_real3d_and_field, Multiply_field_and_real3d
    PROCEDURE               :: Divide_real3d                    ! Divide two variables of type t_jsb_var_real3d
    PROCEDURE               :: Divide_real3d_and_scalar         ! Divide variable of type t_jsb_var_real3d by REAL scalar
    PROCEDURE               :: Divide_real3d_and_field          ! Divide variable of type t_jsb_var_real3d by REAL 2-d field
    PROCEDURE, PASS(this2)  :: Divide_field_and_real3d          ! Divide variable of REAL 2-d field by type t_jsb_var_real3d
    GENERIC, PUBLIC         :: OPERATOR(/) => Divide_real3d, Divide_real3d_and_scalar, &
                                              Divide_real3d_and_field, Divide_field_and_real3d
    PROCEDURE, PUBLIC       :: Init => init_t_jsb_var_real3d
    !FINAL ...
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
  !CONTAINS
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
  !CONTAINS
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

  !-----------------------------------------------------------------------------------------------------!
  ! 1D : operator overload Functions      *slmtemp                                                      !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> 1D + 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real1d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr + var2%ptr
  END FUNCTION Add_real1d

  !-----------------------------------------------------------------------------------------------------
  !> 1D + scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real1d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr + scalar
  END FUNCTION Add_real1d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar + 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_real1d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)

    var3%ptr = scalar + this2%ptr
  END FUNCTION Add_scalar_and_real1d

  !-----------------------------------------------------------------------------------------------------
  !> 1D - 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real1d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr - var2%ptr
  END FUNCTION Subtract_real1d

  !-----------------------------------------------------------------------------------------------------
  !> 1D - scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real1d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr - scalar
  END FUNCTION Subtract_real1d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar - 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_scalar_and_real1d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)

    var3%ptr = scalar - this2%ptr
  END FUNCTION Subtract_scalar_and_real1d

  !-----------------------------------------------------------------------------------------------------
  !> 1D * 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real1d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr * var2%ptr
  END FUNCTION Multiply_real1d

  !-----------------------------------------------------------------------------------------------------
  !> 1D * scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real1d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr * scalar
  END FUNCTION Multiply_real1d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar * 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_and_real1d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)

    var3%ptr = scalar * this2%ptr
  END FUNCTION Multiply_scalar_and_real1d

  !-----------------------------------------------------------------------------------------------------
  !> 1D / 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real1d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr / var2%ptr
  END FUNCTION Divide_real1d

  !-----------------------------------------------------------------------------------------------------
  !> 1D / scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real1d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr = this%ptr / scalar
  END FUNCTION Divide_real1d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar / 1D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_scalar_and_real1d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real1d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real1d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)

    var3%ptr = scalar / this2%ptr
  END FUNCTION Divide_scalar_and_real1d


  !-----------------------------------------------------------------------------------------------------!
  ! 2D : operator overload Functions      *slmjsb4                                                      !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> 2D + 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real2d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) + var2%ptr(:,:)
  END FUNCTION Add_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D + scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real2d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)

    var3%ptr(:,:) = this%ptr(:,:) + scalar
  END FUNCTION Add_real2d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar + 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_real2d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)

    var3%ptr(:,:) = scalar + this2%ptr(:,:)
  END FUNCTION Add_scalar_and_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D + field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real2d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_real2d_and_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) + field(:,:)
  END FUNCTION Add_real2d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field + 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_field_and_real2d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_field_and_real2d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2)))

    var3%ptr(:,:) = field(:,:) + this2%ptr(:,:)
  END FUNCTION Add_field_and_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D - 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real2d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) - var2%ptr(:,:)
  END FUNCTION Subtract_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D - scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real2d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) - scalar
  END FUNCTION Subtract_real2d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar - 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_scalar_and_real2d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2)))

    var3%ptr(:,:) = scalar - this2%ptr(:,:)
  END FUNCTION Subtract_scalar_and_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D - field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real2d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Subtract_real2d_and_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) - field(:,:)
  END FUNCTION Subtract_real2d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field - 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_field_and_real2d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Subtract_field_and_real2d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2)))

    var3%ptr(:,:) = field(:,:) - this2%ptr(:,:)
  END FUNCTION Subtract_field_and_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D * 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real2d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) * var2%ptr(:,:)
  END FUNCTION Multiply_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D * scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real2d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) * scalar
  END FUNCTION Multiply_real2d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar * 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_and_real2d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2)))

    var3%ptr(:,:) = scalar * this2%ptr(:,:)
  END FUNCTION Multiply_scalar_and_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D * field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real2d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_real2d_and_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) * field(:,:)
  END FUNCTION Multiply_real2d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field * 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field_and_real2d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_field_and_real2d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2)))

    var3%ptr(:,:) = field(:,:) * this2%ptr(:,:)
  END FUNCTION Multiply_field_and_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D / 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real2d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) / var2%ptr(:,:)
  END FUNCTION Divide_real2d

  !-----------------------------------------------------------------------------------------------------
  !> 2D / scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real2d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) / scalar
  END FUNCTION Divide_real2d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> 2D / field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real2d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Divide_real2d_by_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2)))

    var3%ptr(:,:) = this%ptr(:,:) / field(:,:)
  END FUNCTION Divide_real2d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field / 2D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_field_and_real2d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real2d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:)
    TYPE(t_jsb_var_real2d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Divide_field_by_real2d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2)))

    var3%ptr(:,:) = field(:,:) / this2%ptr(:,:)
  END FUNCTION Divide_field_and_real2d

  !-----------------------------------------------------------------------------------------------------
  !> Initialisation of var_real2d with a scalar 
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_jsb_var_real2d(this,s)
    IMPLICIT NONE
    CLASS(t_jsb_var_real2d), INTENT(inout) :: this
    REAL(wp),                INTENT(in)    :: s

    this%ptr(:,:)     = s 
  END SUBROUTINE init_t_jsb_var_real2d 


  !-----------------------------------------------------------------------------------------------------!
  ! 3D : operator overload Functions      *slmjsb4                                                      !
  !-----------------------------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------------------------------------
  !> 3D + 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real3d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) + var2%ptr(:,:,:)
  END FUNCTION Add_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D + scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real3d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) + scalar
  END FUNCTION Add_real3d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar + 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_scalar_and_real3d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2), &
!             lbound(this2%ptr,3):ubound(this2%ptr,3)))

    var3%ptr(:,:,:) = scalar + this2%ptr(:,:,:)
  END FUNCTION Add_scalar_and_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D + field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_real3d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_real3d_and_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2) .OR. SIZE(field,3) /= SIZE(this%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) + field(:,:,:)
  END FUNCTION Add_real3d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field + 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Add_field_and_real3d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_field_and_real3d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2) .OR. SIZE(field,3) /= SIZE(this2%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2), &
!             lbound(this2%ptr,3):ubound(this2%ptr,3)))

    var3%ptr(:,:,:) = field(:,:,:) + this2%ptr(:,:,:)
  END FUNCTION Add_field_and_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D - 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real3d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) - var2%ptr(:,:,:)
  END FUNCTION Subtract_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D - scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real3d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) - scalar
  END FUNCTION Subtract_real3d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar - 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_scalar_and_real3d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2), &
!             lbound(this2%ptr,3):ubound(this2%ptr,3)))

    var3%ptr(:,:,:) = scalar - this2%ptr(:,:,:)
  END FUNCTION Subtract_scalar_and_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D - field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_real3d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Subtract_real3d_and_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2) .OR. SIZE(field,3) /= SIZE(this%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) - field(:,:,:)
  END FUNCTION Subtract_real3d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field - 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Subtract_field_and_real3d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Subtract_field_and_real3d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2) .OR. SIZE(field,3) /= SIZE(this2%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2), &
!             lbound(this2%ptr,3):ubound(this2%ptr,3)))

    var3%ptr(:,:,:) = field(:,:,:) - this2%ptr(:,:,:)
  END FUNCTION Subtract_field_and_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D * 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real3d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) * var2%ptr(:,:,:)
  END FUNCTION Multiply_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D * scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real3d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) * scalar
  END FUNCTION Multiply_real3d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> scalar * 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_scalar_and_real3d(scalar,this2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2), &
!             lbound(this2%ptr,3):ubound(this2%ptr,3)))

    var3%ptr(:,:,:) = scalar * this2%ptr(:,:,:)
  END FUNCTION Multiply_scalar_and_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D * field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_real3d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_real3d_and_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2) .OR. SIZE(field,3) /= SIZE(this%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) * field(:,:,:)
  END FUNCTION Multiply_real3d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field * 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Multiply_field_and_real3d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_field_and_real3d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2) .OR. SIZE(field,3) /= SIZE(this2%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2), &
!             lbound(this2%ptr,3):ubound(this2%ptr,3)))

    var3%ptr(:,:,:) = field(:,:,:) * this2%ptr(:,:,:)
  END FUNCTION Multiply_field_and_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D / 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real3d(this, var2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this, var2
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) / var2%ptr(:,:,:)
  END FUNCTION Divide_real3d

  !-----------------------------------------------------------------------------------------------------
  !> 3D / scalar operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real3d_and_scalar(this, scalar) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: scalar
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) / scalar
  END FUNCTION Divide_real3d_and_scalar

  !-----------------------------------------------------------------------------------------------------
  !> 3D / field operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_real3d_and_field(this, field) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Divide_real3d_by_field'

!    IF (SIZE(field,1) /= SIZE(this%ptr,1) .OR. SIZE(field,2) /= SIZE(this%ptr,2) .OR. SIZE(field,3) /= SIZE(this%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this%ptr,1):ubound(this%ptr,1), lbound(this%ptr,2):ubound(this%ptr,2), &
!             lbound(this%ptr,3):ubound(this%ptr,3)))

    var3%ptr(:,:,:) = this%ptr(:,:,:) / field(:,:,:)
  END FUNCTION Divide_real3d_and_field

  !-----------------------------------------------------------------------------------------------------
  !> field / 3D operator overload
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Divide_field_and_real3d(field,this2) RESULT(var3)

    CLASS(t_jsb_var_real3d), INTENT(in) :: this2
    REAL(wp),                INTENT(in) :: field(:,:,:)
    TYPE(t_jsb_var_real3d), ALLOCATABLE :: var3

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Divide_field_by_real3d'

!    IF (SIZE(field,1) /= SIZE(this2%ptr,1) .OR. SIZE(field,2) /= SIZE(this2%ptr,2) .OR. SIZE(field,3) /= SIZE(this2%ptr,3)) &
!      & CALL finish(TRIM(routine), 'Dimension mismatch')

    ALLOCATE(var3, source=this2)
!    NULLIFY(var3%ptr)
!    ALLOCATE(var3%ptr(lbound(this2%ptr,1):ubound(this2%ptr,1), lbound(this2%ptr,2):ubound(this2%ptr,2), &
!             lbound(this2%ptr,3):ubound(this2%ptr,3)))

    var3%ptr(:,:,:) = field(:,:,:) / this2%ptr(:,:,:)
  END FUNCTION Divide_field_and_real3d


  !-----------------------------------------------------------------------------------------------------
  !> Initialisation of var_real3d with a scalar 
  !! 
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE init_t_jsb_var_real3d(this,s)
    IMPLICIT NONE
    CLASS(t_jsb_var_real3d), INTENT(inout) :: this
    REAL(wp),                INTENT(in)    :: s

    this%ptr(:,:,:)     = s 
  END SUBROUTINE init_t_jsb_var_real3d 



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

!!  FUNCTION Add_real2d(var1, var2) RESULT(var3)

!!    CLASS(t_jsb_var_real2d), INTENT(in) :: var1, var2
!!    TYPE(t_jsb_var_real2d)              :: var3

!!    var3%ptr(:,:) = var1%ptr(:,:) + var2%ptr(:,:)

!!  END FUNCTION Add_real2d

!!  FUNCTION Add_real2d_and_scalar(var1, scalar) RESULT(var3)

!!    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
!!    REAL(wp),                INTENT(in) :: scalar
!!    TYPE(t_jsb_var_real2d)              :: var3

!!    var3%ptr(:,:) = var1%ptr(:,:) + scalar

!!  END FUNCTION Add_real2d_and_scalar

!!  FUNCTION Add_real2d_and_field(var1, field) RESULT(var3)

!!    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
!!    REAL(wp),                INTENT(in) :: field(:,:)
!!    TYPE(t_jsb_var_real2d)              :: var3

!!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_real2d_and_field'

!!    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
!!      & CALL finish(TRIM(routine), 'Dimension mismatch')

!!    var3%ptr(:,:) = var1%ptr(:,:) + field(:,:)

!!  END FUNCTION Add_real2d_and_field

!!  FUNCTION Add_real3d(var1, var2) RESULT(var3)

!!    CLASS(t_jsb_var_real3d), INTENT(in) :: var1, var2
!!    TYPE(t_jsb_var_real3d)              :: var3

!!    var3%ptr(:,:,:) = var1%ptr(:,:,:) + var2%ptr(:,:,:)

!!  END FUNCTION Add_real3d

!!  FUNCTION Add_real3d_and_scalar(var1, scalar) RESULT(var3)

!!    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
!!    REAL(wp),                INTENT(in) :: scalar
!!    TYPE(t_jsb_var_real3d)              :: var3

!!    var3%ptr(:,:,:) = var1%ptr(:,:,:) + scalar

!!  END FUNCTION Add_real3d_and_scalar

!!  FUNCTION Add_real3d_and_field(var1, field) RESULT(var3)

!!    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
!!    REAL(wp),                INTENT(in) :: field(:,:,:)
!!    TYPE(t_jsb_var_real3d)              :: var3

!!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Add_real3d_and_field'

!!    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
!!      & CALL finish(TRIM(routine), 'Dimension mismatch')

!!    var3%ptr(:,:,:) = var1%ptr(:,:,:) + field(:,:,:)

!!  END FUNCTION Add_real3d_and_field

!!  FUNCTION Multiply_real2d_with_field(var1, field) RESULT(var3)

!!    CLASS(t_jsb_var_real2d), INTENT(in) :: var1
!!    REAL(wp),                INTENT(in) :: field(:,:)
!!    TYPE(t_jsb_var_real2d)              :: var3

!!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_real2d_with_field'

!!    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2)) &
!!      & CALL finish(TRIM(routine), 'Dimension mismatch')

!!    var3%ptr(:,:) = var1%ptr(:,:) * field(:,:)

!!  END FUNCTION Multiply_real2d_with_field

!!  FUNCTION Multiply_real3d_with_field(var1, field) RESULT(var3)

!!    CLASS(t_jsb_var_real3d), INTENT(in) :: var1
!!    REAL(wp),                INTENT(in) :: field(:,:,:)
!!    TYPE(t_jsb_var_real3d)              :: var3

!!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Multiply_real3d_with_field'

!!    IF (SIZE(field,1) /= SIZE(var1%ptr,1) .OR. SIZE(field,2) /= SIZE(var1%ptr,2) .OR. SIZE(field,3) /= SIZE(var1%ptr,3)) &
!!      & CALL finish(TRIM(routine), 'Dimension mismatch')

!!    var3%ptr(:,:,:) = var1%ptr(:,:,:) + field(:,:,:)

!!  END FUNCTION Multiply_real3d_with_field

!#endif
END MODULE mo_jsb_var_class
