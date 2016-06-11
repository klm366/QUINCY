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

  USE mo_exception, ONLY: message, message_text, warning, finish
  USE mo_util,      ONLY: int2string

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_Message, t_Action, t_State, t_Hsm

  INTEGER, PARAMETER :: MAX_STATE_NESTING = 10
  INTEGER, PARAMETER :: MAX_ACTIONS = 5

  TYPE t_Action
    CHARACTER(len=10) :: name = ""
  CONTAINS
    PROCEDURE :: Is_equal        => Is_equal_action
    GENERIC   :: OPERATOR(==)    => Is_equal
    PROCEDURE :: Is_not_equal    => Is_not_equal_action
    GENERIC   :: OPERATOR(/=)    => Is_not_equal
  END TYPE t_Action

  INTERFACE t_Action
    PROCEDURE Construct_action
  END INTERFACE t_Action

  TYPE t_Message
    CHARACTER(len=50)       :: name = ''
    TYPE(t_Action), POINTER :: action => NULL()
  END TYPE t_Message

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
    INTEGER, PRIVATE            :: no_of_children = -1
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
    PROCEDURE :: Get_parent         => Get_parent_state
    PROCEDURE :: Get_parent_name    => Get_parent_name_state
    PROCEDURE :: Get_no_of_children => Get_no_of_children_state
    PROCEDURE :: Get_pos_of_child   => Get_pos_of_child_state
    PROCEDURE :: Has_children       => Has_children_state
    PROCEDURE :: Has_siblings       => Has_siblings_state
    PROCEDURE :: Is_first_child     => Is_first_child_state
    PROCEDURE :: Is_last_child      => Is_last_child_state
    PROCEDURE :: Print_state
    PROCEDURE :: Is_equal           => Is_equal_state
    GENERIC   :: OPERATOR(==)       => Is_equal
    PROCEDURE :: Is_not_equal       => Is_not_equal_state
    GENERIC   :: OPERATOR(/=)       => Is_not_equal
    PROCEDURE :: Process_message    => Process_message_state
    PROCEDURE(Handler_interface), DEFERRED :: Handler             ! State's handler function
  END TYPE t_State

  ! Derived-type constructor
!!$  INTERFACE t_State
!!$    PROCEDURE Construct_state
!!$  END INTERFACE t_State

  TYPE t_State_p
    CLASS(t_State), POINTER :: p
  END TYPE t_State_p
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
    TYPE(t_Action), POINTER     :: actions(:) => NULL() !< List of possible actions (events)
    INTEGER                     :: levels_to_lca = -1   !< Temporary storage
  CONTAINS
    PROCEDURE :: Init_hsm
    PROCEDURE :: Get_current       => Get_current_hsm
    PROCEDURE :: Goto              => Goto_state_hsm
    PROCEDURE :: Goto_next         => Goto_next_state_hsm
    PROCEDURE :: Get_next          => Get_next_state_hsm
    PROCEDURE :: Get_state         => Get_state_hsm
    PROCEDURE :: Start_state       => Start_state_hsm
    PROCEDURE :: Reset             => Reset_hsm
    PROCEDURE :: Register_action   => Register_action_hsm
    PROCEDURE :: Get_action        => Get_action_hsm
    PROCEDURE :: Print_actions     => Print_actions_hsm
    PROCEDURE :: Process_message   => Process_message_hsm        ! State machine engine
    PROCEDURE :: Process_entry_to_next
    PROCEDURE :: Get_levels_to_LCA => Get_levels_to_LCA_hsm
    PROCEDURE :: Process_exit_to_LCA
    PROCEDURE :: Start             => Start_hsm                  ! Start state machine
  END TYPE t_Hsm

  ! Derived-type constructor
!!$  INTERFACE t_Hsm
!!$    PROCEDURE Construct_hsm
!!$  END INTERFACE t_Hsm

!!$  ABSTRACT INTERFACE
!!$    FUNCTION Handler_interface(this, hsm, msg_in) RESULT(msg_out)
!!$      IMPORT :: t_State, t_Hsm, t_Message
!!$      CLASS(t_State),   INTENT(inout) :: this
!!$      CLASS(t_Hsm),     POINTER       :: hsm
!!$      CLASS(t_Message), POINTER       :: msg_in
!!$      CLASS(t_Message), POINTER       :: msg_out
!!$    END FUNCTION Handler_interface
!!$  END INTERFACE
  ABSTRACT INTERFACE
    FUNCTION Handler_interface(this, hsm, msg_in) RESULT(msg_out)
      IMPORT :: t_State, t_Hsm, t_Message
      CLASS(t_State),   INTENT(inout) :: this
      CLASS(t_Hsm),     INTENT(inout) :: hsm
      CLASS(t_Message), INTENT(in)    :: msg_in
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
      !CALL this%Print_state()
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

  FUNCTION Get_no_of_children_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(inout) :: this
    INTEGER                       :: return_value

    CLASS(t_State), POINTER :: current

    IF (this%no_of_children >= 0) THEN
      return_value = this%no_of_children
      RETURN
    END IF

    IF (this%Has_children()) THEN
      return_value = 1
      current => this%first_child
      DO WHILE (.NOT. current%Is_last_child())
        current => current%next_sibling
        return_value = return_value + 1
      END DO
    ELSE
      return_value = 0
    END IF
    this%no_of_children = return_value

  END FUNCTION Get_no_of_children_state

  FUNCTION Get_pos_of_child_state(this) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this
    INTEGER                    :: return_value

    INTEGER :: depth

    IF (ALL(this%path > 0)) THEN
      depth = SIZE(this%path)
    ELSE
      depth = MINLOC(this%path, DIM=1) - 1
    END IF

    return_value = this%path(depth)

  END FUNCTION

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
    IF (this%visited) THEN
      CALL message('  Visited', '.TRUE.')
    ELSE
      CALL message('  Visited', '.FALSE.')
    END IF

  END SUBROUTINE Print_state

  ! Only works like this for two states in the same state machine!
  FUNCTION Is_equal_state(this, that) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = ALL(this%path == that%path)

  END FUNCTION Is_equal_state

  FUNCTION Is_not_equal_state(this, that) RESULT(return_value)

    CLASS(t_State), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = ANY(this%path /= that%path)

  END FUNCTION Is_not_equal_state

!!$  FUNCTION Process_message_state(this, hsm, msg) RESULT(return_ptr)
!!$
!!$    CLASS(t_State),   INTENT(inout) :: this
!!$    CLASS(t_Hsm),     POINTER       :: hsm
!!$    CLASS(t_Message), POINTER       :: msg
!!$    CLASS(t_Message), POINTER       :: return_ptr
!!$
!!$    CALL this%Print_state()
!!$    return_ptr => this%Handler(hsm, msg)
!!$
!!$  END FUNCTION Process_message_state
  FUNCTION Process_message_state(this, hsm, msg) RESULT(return_ptr)

    CLASS(t_State),   INTENT(inout) :: this
    CLASS(t_Hsm),     INTENT(inout) :: hsm
    CLASS(t_Message), INTENT(in)    :: msg
!!$    CLASS(t_Message), POINTER       :: msg
    CLASS(t_Message), POINTER       :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_message_state'

!!$    this%msg => msg
!!$    this%hsm => hsm
    !CALL this%Print_state()
    return_ptr => this%Handler(hsm, msg)

    SELECT CASE (TRIM(msg%action%name))
    CASE ('START')
      this%visited = .TRUE.
    CASE ('EXIT')
      this%visited = .FALSE.
    END SELECT

  END FUNCTION Process_message_state

  ! ===============================================================================
  ! Methods for hsm class
  ! ===============================================================================
  SUBROUTINE Init_hsm(self, name)

    CLASS(t_Hsm), INTENT(inout) :: self
    CHARACTER(len=*), INTENT(in)         :: name

    TYPE(t_Action), POINTER :: action
    INTEGER :: i

    ALLOCATE(self%actions(MAX_ACTIONS))
    NULLIFY(action)
    DO i=1,MAX_ACTIONS
      !self%actions(i) = t_Action("")     ! This generates an internal compiler error with gfortran
      action => t_Action("")
      self%actions(i) = action
      DEALLOCATE(action)
    END DO
    CALL self%Register_action(t_Action("ENTER"))
    CALL self%Register_action(t_Action("EXIT"))
    CALL self%Register_action(t_Action("START"))

    self%name = TRIM(name)

  END SUBROUTINE Init_hsm

  FUNCTION Get_current_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm),   INTENT(in) :: this
    CLASS(t_State), POINTER    :: return_ptr

    return_ptr => this%current

  END FUNCTION Get_current_hsm

  SUBROUTINE Reset_hsm(this)

    CLASS(t_hsm), INTENT(inout) :: this

    this%current => this%top
    this%top%visited = .FALSE.

  END SUBROUTINE Reset_hsm

  ! Get next state in pre-order traversal
  FUNCTION Get_next_state_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(inout) :: this
    CLASS(t_State), POINTER     :: return_ptr

    IF (this%current%visited) THEN   ! State has been visited before on traversal
      this%current%visited = .FALSE.
      IF (this%current == this%top) THEN
        return_ptr => NULL()
        RETURN
      ELSE
        IF (ASSOCIATED(this%current%next_sibling)) THEN
          return_ptr => this%current%next_sibling
        ELSE
          return_ptr => this%current%parent
        END IF
      END IF
    ELSE                             ! First visit of this state during traversal
      IF (this%current%Has_children()) THEN
        this%current%visited = .TRUE.
        return_ptr => this%current%first_child
      ELSE IF (ASSOCIATED(this%current%next_sibling)) THEN
        return_ptr => this%current%next_sibling
      ELSE
        return_ptr => this%current%parent
      END IF
    END IF

  END FUNCTION Get_next_state_hsm

  ! Set Hsm to next state in pre-order traversal
  FUNCTION Goto_next_state_hsm(this) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(inout) :: this
    CLASS(t_State), POINTER     :: return_ptr

    !this%current%visited = .NOT. this%current%visited
    this%current => this%Get_next()
    return_ptr => this%current

  END FUNCTION Goto_next_state_hsm

  FUNCTION Get_state_hsm(this, path) RESULT(return_ptr)

    CLASS(t_hsm), INTENT(in) :: this
    INTEGER,      INTENT(in) :: path(:)
    CLASS(t_State), POINTER  :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_state_hsm'

    INTEGER :: depth, i, j

    IF (SIZE(path) > MAX_STATE_NESTING) &
      & CALL finish(TRIM(routine), 'Path too long')

    return_ptr => this%top

    IF (ALL(path > 0)) THEN
      depth = SIZE(path)
    ELSE
      depth = MINLOC(path, DIM=1) - 1
    END IF
    DO i=2,depth
      return_ptr => return_ptr%first_child
      DO j=2,path(i)
        return_ptr => return_ptr%next_sibling
      END DO
    END DO

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Path to invalid state')

  END FUNCTION Get_state_hsm

  ! State transition
  SUBROUTINE Goto_state_hsm(this, that)

    CLASS(t_hsm), INTENT(inout) :: this
    CLASS(t_State), INTENT(in), TARGET :: that

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Goto_state_hsm'

    IF (ASSOCIATED(this%next)) &
      & CALL finish(TRIM(routine), 'Transition already in progress')

    !CALL this%Process_exit_to_LCA(that)
    this%levels_to_lca = this%Get_levels_to_LCA(that)
    this%next => that

  END SUBROUTINE Goto_state_hsm
  !>
  !! Get number of levels from state machines source state to the least common
  !! ancestor with some other target state.
  !!
  FUNCTION Get_levels_to_LCA_hsm(this, that) RESULT(return_value)

    CLASS(t_hsm),   INTENT(in)         :: this
    CLASS(t_State), INTENT(in), TARGET :: that
    INTEGER                            :: return_value

    CLASS(t_State), POINTER :: sstate, tstate

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_levels_to_LCA_hsm'

    return_value = 0

    sstate => this%source
    IF (sstate == that) THEN
      return_value = 1
      RETURN
    END IF

    DO WHILE (ASSOCIATED(sstate))
      tstate => that
      DO WHILE (ASSOCIATED(tstate))
        IF (sstate == tstate) THEN
!!$          CALL message(TRIM(routine), 'Levels to LCA from '//TRIM(this%source%name)//' to '//&
!!$            &                                                TRIM(that%name)//': '//int2string(return_value))
          RETURN
        END IF
        tstate => tstate%parent
      END DO
      return_value = return_value + 1
      sstate => sstate%parent
    END DO

  END FUNCTION Get_levels_to_LCA_hsm

  ! Take start transition (no states need to be exited)
  SUBROUTINE Start_state_hsm(this, that)

    CLASS(t_Hsm),   INTENT(inout)      :: this
    CLASS(t_State), INTENT(in), TARGET :: that

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Start_state_hsm'

!!$    CALL message(TRIM(routine), 'Taking start transition from '//TRIM(this%current%name)//&
!!$      &                                                  ' to '//TRIM(that%name))

    IF (ASSOCIATED(this%next)) &
      & CALL warning(TRIM(routine), 'Transition to '//TRIM(this%next%name)//' already in progress')

    this%next => that

  END SUBROUTINE Start_state_hsm

  SUBROUTINE Start_hsm(this)

    CLASS(t_Hsm), INTENT(inout) :: this

    TYPE(t_Message)           :: start_msg, entry_msg
    CLASS(t_Message), POINTER :: tmp_msg

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Start_hsm'

!!$    CALL message(TRIM(routine), 'Starting hsm ...')

    start_msg%action => this%Get_action('START')
    entry_msg%action => this%Get_action('ENTER')

    this%current => this%top
    !this%next    => NULL()
    tmp_msg => this%current%Process_message(this, entry_msg)

    tmp_msg => this%current%Process_message(this, start_msg)
    DO WHILE (ASSOCIATED(this%next))
      CALL this%Process_entry_to_next()
      this%current => this%next
      NULLIFY(this%next)
      tmp_msg => this%current%Process_message(this, start_msg)
    END DO

!!$    CALL message(TRIM(routine), 'Hsm started.')

  END SUBROUTINE Start_hsm
  !
  ! State machine engine
  !
  SUBROUTINE Process_message_hsm(this, msg)

    CLASS(t_Hsm),     TARGET, INTENT(inout) :: this
!!$    CLASS(t_Message), POINTER       :: msg
    CLASS(t_Message), INTENT(in) :: msg

    CLASS(t_State),   POINTER :: state
    CLASS(t_Message), POINTER :: new_msg, tmp_msg
    TYPE(t_Message)           :: start_msg
    !CLASS(t_Hsm),     POINTER :: hsm

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_message_hsm'

!!$    CALL message(TRIM(routine), 'Processing message '//TRIM(msg%name)//' on hsm ...')

    !hsm => this
    start_msg%action => this%Get_action('START')

!!$    CALL this%Process_exit_to_LCA()

    ALLOCATE(new_msg, source=msg)
    state => this%current
    DO WHILE (ASSOCIATED(state))
      this%source => state
      new_msg => state%Process_message(this, new_msg)
      CALL this%Process_exit_to_LCA()
      IF (ASSOCIATED(this%next)) THEN      ! State transition taken (i.e. already exited to LCA)
        ! Execute entry and start messages to target
        CALL this%Process_entry_to_next()
        this%current => this%next
        NULLIFY(this%next)
        tmp_msg => this%current%Process_message(this, start_msg)
        DO WHILE (ASSOCIATED(this%next))
          CALL this%Process_entry_to_next()
          this%current => this%next
          NULLIFY(this%next)
          tmp_msg => this%current%Process_message(this, start_msg)
        END DO
        IF (.NOT. ASSOCIATED(new_msg)) THEN ! Message processed
          EXIT
        ELSE
          state => this%current
        END IF
      ELSE                                ! No state transition taken
        IF (ASSOCIATED(new_msg)) THEN
          state => state%parent
        ELSE
          EXIT
        END IF
      END IF
    END DO

!!$    CALL message(TRIM(routine), 'Processing finished.')

  END SUBROUTINE Process_message_hsm

  ! Process entry messages from LCA to target
  SUBROUTINE Process_entry_to_next(this)

    CLASS(t_Hsm), INTENT(inout) :: this

    TYPE(t_Message) :: entry_msg
    CLASS(t_Message), POINTER :: tmp_msg
    CLASS(t_State), POINTER :: state
    TYPE(t_State_p) :: entrypath(MAX_STATE_NESTING)
    INTEGER         :: i, j

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_entry_to_next'

!!$    CALL message(TRIM(routine), 'Processing entry messages to next ...')

    entry_msg%action => this%Get_action('ENTER')

    i = 0
    state => this%next
    DO WHILE (state /= this%current)
      i = i + 1
      entrypath(i)%p => state
      state => state%parent
      IF (.NOT. ASSOCIATED(state)) EXIT
    END DO

    DO j=i,1,-1
      state => entrypath(j)%p
      tmp_msg => state%Process_message(this, entry_msg)
    END DO

!!$    CALL message(TRIM(routine), 'Processing of entry messages finished.')

  END SUBROUTINE Process_entry_to_next

  ! Process exit messages from current state up to LCA with target state
  ! (set during a previous state transition)
!  SUBROUTINE Process_exit_to_lca(this, that)
  SUBROUTINE Process_exit_to_lca(this)

    CLASS(t_Hsm),   INTENT(inout) :: this
!    CLASS(t_State), INTENT(in)    :: that

    TYPE(t_Message)           :: exit_msg
    CLASS(t_Message), POINTER :: tmp_msg
    CLASS(t_State),   POINTER :: state
    INTEGER                   :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Process_exit_to_lca'

    IF (this%levels_to_lca == -1) RETURN

!!$    CALL message(TRIM(routine), 'Processing exit messages to lca ...'//int2string(this%levels_to_lca))

    exit_msg%action => this%Get_action('EXIT')

    state => this%current

    IF (ASSOCIATED(this%source)) THEN
      DO WHILE (state /= this%source)
        tmp_msg => state%Process_message(this, exit_msg)
        state => state%parent
        IF (.NOT. ASSOCIATED(state)) EXIT
      END DO
    END IF

!    DO i=this%Get_levels_to_LCA(that),1,-1
    DO i=this%levels_to_lca,1,-1
      tmp_msg => state%Process_message(this, exit_msg)
      state => state%parent
      IF (.NOT. ASSOCIATED(state)) EXIT
    END DO

    this%current => state
    this%levels_to_lca = -1

!!$    CALL message(TRIM(routine), 'Processing of exit messages finished.')

  END SUBROUTINE Process_exit_to_lca

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

  FUNCTION Get_action_hsm(this, name) RESULT(return_ptr)

    CLASS(t_Hsm),     TARGET, INTENT(in) :: this
    CHARACTER(len=*),         INTENT(in) :: name
    CLASS(t_Action),  POINTER            :: return_ptr

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_action_hsm'

    return_ptr => NULL()

    DO i=1,MAX_ACTIONS
      IF (TRIM(this%actions(i)%name) == TRIM(name)) THEN
        ALLOCATE(return_ptr, source=this%actions(i))
        EXIT
      END IF
    END DO

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Action '//TRIM(name)//' not found')

  END FUNCTION Get_action_hsm

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

  FUNCTION Is_not_equal_action(this, that) RESULT(return_value)

    CLASS(t_Action), INTENT(in) :: this, that
    LOGICAL                    :: return_value

    return_value = TRIM(this%name) /= TRIM(that%name)

  END FUNCTION Is_not_equal_action

END MODULE mo_hsm_class
