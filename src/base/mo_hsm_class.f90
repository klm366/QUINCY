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

  USE mo_exception, ONLY: message, message_text, finish

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_Message, t_Action, t_Queue, t_State, t_Hsm, Handler_interface

  INTEGER, PARAMETER :: MAX_STATE_NESTING = 10
  INTEGER, PARAMETER :: MAX_ACTIONS = 5

  TYPE t_Action
    CHARACTER(len=10) :: name = ""
  CONTAINS
    PROCEDURE :: Is_equal        => Is_equal_action
    GENERIC   :: OPERATOR(==)    => Is_equal
  END TYPE t_Action

  INTERFACE t_Action
    PROCEDURE Construct_action
  END INTERFACE t_Action

  TYPE, ABSTRACT :: t_Message
    CHARACTER(len=50)       :: name = ''
    TYPE(t_Action), POINTER :: action => NULL()
  END TYPE t_Message

  TYPE t_Queue
    CLASS(t_Message), POINTER :: message
    TYPE(t_Queue),    POINTER :: next => NULL()
  END TYPE t_Queue
  !>
  !! Type for States of the State Machine
  !!
  !! This derived type defines the state class
  !!
  TYPE, ABSTRACT :: t_State
    CHARACTER(len=10)           :: name
    CHARACTER(len=50)           :: description
!!$    INTEGER                     :: id
    INTEGER                     :: path(MAX_STATE_NESTING) = -1
    !
    CLASS(t_Queue), POINTER     :: queue => NULL()
    !
    LOGICAL                     :: visited = .FALSE.
    CLASS(t_State), POINTER     :: parent       => NULL()
    CLASS(t_State), POINTER     :: first_child  => NULL()
    CLASS(t_State), POINTER     :: last_child   => NULL()
    !
    CLASS(t_State), POINTER     :: next_sibling => NULL()
    CLASS(t_State), POINTER     :: prev_sibling => NULL()
!!$    PROCEDURE(Handler_interface), POINTER :: Handler           ! State's handler function
  CONTAINS
    PROCEDURE :: Init_state
    PROCEDURE :: Get_parent      => Get_parent_state
    PROCEDURE :: Get_parent_name => Get_parent_name_state
    PROCEDURE :: Has_children    => Has_children_state
    PROCEDURE :: Has_siblings    => Has_siblings_state
    PROCEDURE :: Is_first_child  => Is_first_child_state
    PROCEDURE :: Is_last_child   => Is_last_child_state
    PROCEDURE :: Print_state
    PROCEDURE :: Is_equal        => Is_equal_state
    GENERIC   :: OPERATOR(==)    => Is_equal
    PROCEDURE(Handler_interface), DEFERRED :: Handler
    PROCEDURE :: Process_message => Process_message_state
  END TYPE t_State

  ! Derived-type constructor
!!$  INTERFACE t_State
!!$    PROCEDURE Construct_state
!!$  END INTERFACE t_State

  !>
  !! Type for Hierarchical State Machine
  !!
  !! This derived type defines the structure of a basic Hierarchical State Machine
  !!
  TYPE, ABSTRACT :: t_Hsm
    CHARACTER(len=30)           :: name
    CLASS(t_State), POINTER     :: top     => NULL()    !< Topmost state
    !
    CLASS(t_State), POINTER     :: current => NULL()    !< Current state
    CLASS(t_State), POINTER     :: next    => NULL()    !< Next state (non_null if transition taken)
    CLASS(t_State), POINTER     :: source  => NULL()    !< Source state during last transition
    TYPE(t_Action)              :: actions(MAX_ACTIONS) !< List of possible actions (events)
  CONTAINS
    PROCEDURE :: Init_hsm
    PROCEDURE :: Get_current       => Get_current_hsm
    PROCEDURE :: Get_levels_to_LCA => Get_levels_to_LCA_hsm
    PROCEDURE :: Goto_next         => Goto_next_state_hsm
    PROCEDURE :: Get_state         => Get_state_hsm
    PROCEDURE :: Reset_state       => Reset_state_hsm
    PROCEDURE :: Register_action   => Register_action_hsm
    PROCEDURE :: Print_actions     => Print_actions_hsm
    PROCEDURE :: Process_message   => Process_message_hsm        ! State machine engine
  END TYPE t_Hsm

  ! Derived-type constructor
!!$  INTERFACE t_Hsm
!!$    PROCEDURE Construct_hsm
!!$  END INTERFACE t_Hsm

  ABSTRACT INTERFACE
    FUNCTION Handler_interface(this, hsm, msg_in) RESULT(msg_out)
      IMPORT :: t_State, t_Hsm, t_Message
      CLASS(t_State),   INTENT(inout) :: this
      CLASS(t_Hsm),     POINTER       :: hsm
      CLASS(t_Message), POINTER       :: msg_in
      CLASS(t_Message), POINTER       :: msg_out
    END FUNCTION Handler_interface
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hsm_class'

CONTAINS

  ! ===============================================================================
  ! Methods for state class
  ! ===============================================================================
!!$  FUNCTION Construct_state(name, description, parent, handler) RESULT(return_ptr)
!!$
!!$    CHARACTER(len=*), INTENT(in)      :: name
!!$    CHARACTER(len=*), INTENT(in)      :: description
!!$    CLASS(t_State), POINTER, OPTIONAL :: parent
!!$    PROCEDURE(Handler_interface), POINTER :: handler
!!$    PROCEDURE(), POINTER :: handler
!!$    CLASS(t_State), POINTER           :: return_ptr
!!$
!!$    INTEGER :: depth
!!$
!!$    ALLOCATE(return_ptr)
!!$
!!$    return_ptr%name = TRIM(name)
!!$    return_ptr%description = TRIM(description)
!!$    return_ptr%Handler => handler
!!$
!!$    IF (.NOT. PRESENT(parent)) THEN
!!$      return_ptr%path(:) = -1
!!$      return_ptr%path(1) = 1
!!$      CALL return_ptr%Print()
!!$      RETURN
!!$    END IF
!!$
!!$    return_ptr%path(:) = parent%path(:)
!!$    depth = MINLOC(parent%path, DIM=1)  ! Index of first -1 in path of parent
!!$
!!$    IF (.NOT. ASSOCIATED(parent%first_child)) THEN
!!$      parent%first_child => return_ptr
!!$      return_ptr%path(depth) = 1
!!$    END IF
!!$
!!$    IF (ASSOCIATED(parent%last_child)) THEN
!!$      parent%last_child%next_sibling => return_ptr
!!$      return_ptr%prev_sibling        => parent%last_child
!!$      return_ptr%path(depth) = return_ptr%prev_sibling%path(depth) + 1
!!$    END IF
!!$
!!$    parent%last_child  => return_ptr
!!$
!!$    return_ptr%parent => parent
!!$
!!$    CALL return_ptr%Print()
!!$
!!$  END FUNCTION Construct_state

!!$  SUBROUTINE Init_state(this, name, description, parent, handler)
  SUBROUTINE Init_state(this, name, description, parent)
    CLASS(t_State),   INTENT(inout), TARGET :: this
    CHARACTER(len=*), INTENT(in)            :: name
    CHARACTER(len=*), INTENT(in)            :: description
    CLASS(t_State), POINTER, OPTIONAL       :: parent
!!$    PROCEDURE(Handler_interface), POINTER   :: handler

    INTEGER :: depth

    this%name = TRIM(name)
    this%description = TRIM(description)
!!$    this%Handler => handler

    IF (.NOT. PRESENT(parent)) THEN
      this%path(:) = -1
      this%path(1) = 1
      CALL this%Print_state()
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

    ! CALL this%Print()

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

  FUNCTION Is_first_child_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    LOGICAL                   :: return_value

    return_value = .NOT. ASSOCIATED(this%prev_sibling)

  END FUNCTION Is_first_child_state

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

  FUNCTION Process_message_state(this, hsm, msg) RESULT(return_ptr)

    CLASS(t_State),   INTENT(inout) :: this
    CLASS(t_Hsm),     POINTER       :: hsm
    CLASS(t_Message), POINTER       :: msg
    CLASS(t_Message), POINTER       :: return_ptr

    CALL this%Print_state()
    return_ptr => this%Handler(hsm, msg)

  END FUNCTION Process_message_state
  ! ===============================================================================
  ! Methods for hsm class
  ! ===============================================================================
  SUBROUTINE Init_hsm(self, name)

    CLASS(t_Hsm), INTENT(inout) :: self
    CHARACTER(len=*), INTENT(in)         :: name

    INTEGER :: i

    DO i=1,MAX_ACTIONS
      self%actions(i) = t_Action("")
    END DO
    CALL self%Register_action(t_Action("ENTER"))
    CALL self%Register_action(t_Action("EXIT"))

    self%name = TRIM(name)

  END SUBROUTINE Init_hsm

  FUNCTION Get_current_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm),   INTENT(in) :: this
    CLASS(t_State), POINTER    :: return_ptr

    return_ptr => this%current

  END FUNCTION Get_current_hsm

  SUBROUTINE Reset_state_hsm(this)

    CLASS(t_hsm), INTENT(inout) :: this

    this%current => this%top
    this%top%visited = .FALSE.

  END SUBROUTINE Reset_state_hsm

  ! Set Hsm to next state in pre-order traversal
  FUNCTION Goto_next_state_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(inout) :: this
    CLASS(t_State), POINTER     :: return_ptr

    IF (this%current%visited) THEN   ! State has been visited before on traversal
      this%current%visited = .FALSE.
      IF (this%current == this%top) THEN
        return_ptr => NULL()
        RETURN
      ELSE
        IF (ASSOCIATED(this%current%next_sibling)) THEN
          this%current => this%current%next_sibling
        ELSE
          this%current => this%current%parent
        END IF
      END IF
    ELSE                             ! First visit of this state during traversal
      IF (this%current%Has_children()) THEN
        this%current%visited = .TRUE.
        this%current => this%current%first_child
      ELSE IF (ASSOCIATED(this%current%next_sibling)) THEN
        this%current => this%current%next_sibling
      ELSE
        this%current => this%current%parent
      END IF
    END IF
    return_ptr => this%current

  END FUNCTION Goto_next_state_hsm

  FUNCTION Get_state_hsm(this, path) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(in) :: this
    INTEGER,      INTENT(in) :: path(MAX_STATE_NESTING)
    CLASS(t_State), POINTER  :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_state_hsm'

    INTEGER :: depth, i, j

    return_ptr => this%top

    depth = MINLOC(path, DIM=1) - 1
    DO i=2,depth
      return_ptr => return_ptr%first_child
      DO j=2,path(i)
        return_ptr => return_ptr%next_sibling
      END DO
    END DO

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Path to invalid state')

  END FUNCTION Get_state_hsm
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

  !
  ! State machine engine
  !
  SUBROUTINE Process_message_hsm(this, msg)

    CLASS(t_Hsm),     TARGET, INTENT(inout) :: this
    CLASS(t_Message), POINTER       :: msg

    CLASS(t_State),   POINTER :: state
    CLASS(t_Message), POINTER :: new_msg
    CLASS(t_Hsm),     POINTER :: hsm

    hsm => this

    state => this%current
    this%source => state

!!$    new_msg => state%Process_message(hsm, msg)
    new_msg => state%Handler(hsm, msg)

  END SUBROUTINE Process_message_hsm

  SUBROUTINE Register_action_hsm(this, action)

    CLASS(t_Hsm),   INTENT(inout) :: this
    CLASS(t_Action), INTENT(in)    :: action

    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Register_action_hsm'

    i = 1
    DO WHILE (TRIM(this%actions(i)%name) /= "")
      i = i + 1
    END DO
    IF (i > MAX_ACTIONS) THEN
      CALL finish(TRIM(routine), 'Maximum number of actions exceeded!')
    ELSE
      this%actions(i) = action
      CALL message(TRIM(routine), 'Registering state machine action '//TRIM(this%actions(i)%name))
    END IF

  END SUBROUTINE Register_action_hsm

  SUBROUTINE Print_actions_hsm(this)

    CLASS(t_Hsm), INTENT(in) :: this

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Print_actions_hsm'

    message_text = "State machine actions:"
    DO i=1,MAX_ACTIONS
      IF (TRIM(this%actions(i)%name) /= "") &
        & message_text = TRIM(message_text)//" "//TRIM(this%actions(i)%name)
    END DO
    CALL message(TRIM(routine), TRIM(message_text)//'.')

  END SUBROUTINE Print_actions_hsm

  ! ===============================================================================
  ! Methods for t_Action class
  ! ===============================================================================
  FUNCTION Construct_action(name) RESULT(return_ptr)

    CHARACTER(len=*), INTENT(in) :: name
    CLASS(t_Action), POINTER :: return_ptr

    ALLOCATE(t_Action::return_ptr)
    return_ptr%name = TRIM(name)

  END FUNCTION Construct_action

  FUNCTION Is_equal_action(this, that) RESULT(return_value)

    CLASS(t_Action), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = TRIM(this%name) == TRIM(that%name)

  END FUNCTION Is_equal_action

END MODULE mo_hsm_class
