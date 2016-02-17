!>
!! @brief Definition of Hierarchical State Machine.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2015-06-01)
!!
MODULE mo_hsm_class

  USE mo_exception, ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_State, t_Hsm

  INTEGER, PARAMETER :: MAX_STATE_NESTING = 10
  !>
  !! Type for States of the State Machine
  !!
  !! This derived type defines the state class
  !!
  TYPE t_State
    CHARACTER(len=10)           :: name
    CHARACTER(len=50)           :: description
!!$    INTEGER                     :: id
    INTEGER                     :: path(MAX_STATE_NESTING) = -1
    !
    CLASS(t_State), POINTER     :: parent       => NULL()
    CLASS(t_State), POINTER     :: first_child  => NULL()
    CLASS(t_State), POINTER     :: last_child   => NULL()
    !
    CLASS(t_State), POINTER     :: next_sibling => NULL()
    CLASS(t_State), POINTER     :: prev_sibling => NULL()
  CONTAINS
    PROCEDURE :: Init          => Init_state
    PROCEDURE :: Get_parent    => Get_parent_state
    PROCEDURE :: Get_parent_name => Get_parent_name_state
    PROCEDURE :: Has_children  => Has_children_state
    PROCEDURE :: Has_siblings  => Has_siblings_state
    PROCEDURE :: Is_last_child => Is_last_child_state
    PROCEDURE :: Print         => Print_state
    PROCEDURE :: Is_equal      => Is_equal_state
    GENERIC   :: OPERATOR(==)  => Is_equal
  END TYPE t_State

  ! Derived-type constructor
  INTERFACE t_State
    PROCEDURE Construct_state
  END INTERFACE t_State

  !>
  !! Type for Hierarchical State Machine
  !!
  !! This derived type defines the structure of a basic Hierarchical State Machine
  !!
  TYPE t_Hsm
    CHARACTER(len=10)           :: name
    CLASS(t_State), POINTER     :: top     => NULL()  !< Topmost state
    !
    CLASS(t_State), POINTER     :: current => NULL()  !< Current state
    CLASS(t_State), POINTER     :: next    => NULL()  !< Next state (non_null if transition taken)
    CLASS(t_State), POINTER     :: source  => NULL()  !< Source state during last transition
  CONTAINS
    PROCEDURE :: Get_current       => Get_current_hsm
    PROCEDURE :: Get_levels_to_LCA => Get_levels_to_LCA_hsm
  END TYPE t_Hsm

  ! Derived-type constructor
  INTERFACE t_Hsm
    PROCEDURE Construct_hsm
  END INTERFACE t_Hsm

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hsm_class'

CONTAINS

  ! ===============================================================================
  ! Methods for state class
  ! ===============================================================================
  FUNCTION Construct_state(name, description, parent) RESULT(return_ptr)

    CHARACTER(len=*), INTENT(in)      :: name
    CHARACTER(len=*), INTENT(in)      :: description
    CLASS(t_State), POINTER, OPTIONAL :: parent
    CLASS(t_State), POINTER           :: return_ptr

    INTEGER :: depth

    ALLOCATE(return_ptr)

    return_ptr%name = TRIM(name)
    return_ptr%description = TRIM(description)

    IF (.NOT. PRESENT(parent)) THEN
      return_ptr%path(:) = -1
      return_ptr%path(1) = 1
      CALL return_ptr%Print()
      RETURN
    END IF

    return_ptr%path(:) = parent%path(:)
    depth = MINLOC(parent%path, DIM=1)  ! Index of first -1 in path of parent

    IF (.NOT. ASSOCIATED(parent%first_child)) THEN
      parent%first_child => return_ptr
      return_ptr%path(depth) = 1
    END IF

    IF (ASSOCIATED(parent%last_child)) THEN
      parent%last_child%next_sibling => return_ptr
      return_ptr%prev_sibling        => parent%last_child
      return_ptr%path(depth) = return_ptr%prev_sibling%path(depth) + 1
    END IF

    parent%last_child  => return_ptr

    return_ptr%parent => parent

    CALL return_ptr%Print()

  END FUNCTION Construct_state

  SUBROUTINE Init_state(this, name, description, parent)

    CLASS(t_State),   INTENT(inout), TARGET :: this
    CHARACTER(len=*), INTENT(in)            :: name
    CHARACTER(len=*), INTENT(in)            :: description
    CLASS(t_State), POINTER, OPTIONAL       :: parent

    INTEGER :: depth

    this%name = TRIM(name)
    this%description = TRIM(description)

    IF (.NOT. PRESENT(parent)) THEN
      this%path(:) = -1
      this%path(1) = 1
      CALL this%Print()
      RETURN
    END IF

    SELECT TYPE(parent)
    CLASS IS (t_State)

      this%path(:) = parent%path(:)
      depth = MINLOC(parent%path, DIM=1)  ! Index of first -1 in path of parent

      IF (.NOT. ASSOCIATED(parent%first_child)) THEN
        parent%first_child => this
        this%path(depth) = 1
      END IF

      IF (ASSOCIATED(parent%last_child)) THEN
        parent%last_child%next_sibling => this
        this%prev_sibling        => parent%last_child
        this%path(depth) = this%prev_sibling%path(depth) + 1
      END IF

      parent%last_child  => this

      this%parent => parent

    END SELECT

    CALL this%Print()

  END SUBROUTINE Init_state

  FUNCTION Get_parent_state(this) RESULT(return_ptr)

    CLASS(t_State), INTENT(in) :: this
    CLASS(*), POINTER    :: return_ptr

    CLASS(*), POINTER :: p

!!$    return_ptr => this%parent
    IF (ASSOCIATED(this%parent)) THEN
      p => this%parent
      SELECT TYPE(p)
      CLASS IS (t_State)
        return_ptr => p
      END SELECT
    ELSE
      return_ptr => NULL()
    END IF

  END FUNCTION Get_parent_state

  FUNCTION Get_parent_name_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    CHARACTER(LEN=10)          :: return_value

    CLASS(*), POINTER :: p

    p => this%parent
    IF (ASSOCIATED(p)) THEN
      SELECT TYPE(p)
      CLASS IS (t_State)
        return_value = p%name
      END SELECT
    ELSE
      return_value = ''
    END IF

  END FUNCTION Get_parent_name_state

  FUNCTION Has_children_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                   :: return_value

    return_value = ASSOCIATED(this%first_child)

  END FUNCTION Has_children_state

  FUNCTION Has_siblings_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                    :: return_value

    return_value = ASSOCIATED(this%next_sibling) .OR. ASSOCIATED(this%prev_sibling)

  END FUNCTION Has_siblings_state

  FUNCTION Is_last_child_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                   :: return_value

    return_value = .NOT. ASSOCIATED(this%next_sibling)

  END FUNCTION Is_last_child_state

  SUBROUTINE Print_state(this)

    CLASS(t_State), INTENT(in) :: this

    CHARACTER(LEN=2*MAX_STATE_NESTING) :: c_path
    INTEGER :: i

    DO i=1,MAX_STATE_NESTING
      IF (this%path(i) /= -1) THEN
        WRITE(c_path(2*i-1:2*i),'(I1,A1)') this%path(i), '-'
      ELSE
        c_path(2*i-2:2*MAX_STATE_NESTING) = ' '
        EXIT
      END IF
    END DO

    CALL message('State', this%name)
    CALL message('  Description', this%description)
    CALL message('  Tree path', TRIM(c_path))

  END SUBROUTINE Print_state

  ! Only works like this for two states in the same state machine!
  FUNCTION Is_equal_state(this, that) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = ALL(this%path == that%path)

  END FUNCTION Is_equal_state

  ! ===============================================================================
  ! Methods for hsm class
  ! ===============================================================================
  FUNCTION Construct_hsm(name, top) RESULT(return_ptr)

    CHARACTER(len=*), INTENT(in)         :: name
    CLASS(t_State),   INTENT(in), POINTER, OPTIONAL :: top
    TYPE(t_Hsm), POINTER                 :: return_ptr

    ALLOCATE(return_ptr)

    return_ptr%name = TRIM(name)
    IF (PRESENT(top)) THEN
      return_ptr%top => top
    ELSE
      return_ptr%top => t_State("top", "Top")
    END IF

  END FUNCTION Construct_hsm

  FUNCTION Get_current_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm),   INTENT(in) :: this
    CLASS(t_State), POINTER    :: return_ptr

    return_ptr => this%current

  END FUNCTION Get_current_hsm

  !>
  !! Get number of levels from state machines current state to the least common
  !! ancestor with some other state.
  !!
  FUNCTION Get_levels_to_LCA_hsm(this, that) RESULT(return_value)

    CLASS(t_hsm),   INTENT(in) :: this
    CLASS(t_State), INTENT(in) :: that
    INTEGER                    :: return_value

    return_value = 0

  END FUNCTION Get_levels_to_LCA_hsm

END MODULE mo_hsm_class
