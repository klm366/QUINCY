!>
!! @brief Contains main JSBACH model class.
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_jsb_model_class
#ifndef __NO_JSBACH__

  USE mo_hsm_class,        ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_config_class, ONLY: t_jsb_model_config
  USE mo_jsb_process_class, ONLY: t_jsb_process_task, t_jsb_process_p

  USE mo_exception,        ONLY: message_text, message, finish
  USE mo_io_units,         ONLY: filename_max

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_model_m, t_jsb_model
  PUBLIC :: new_model, print_model, delete_model

  !>
  !! Type for model structure
  !!
  !! This derived type defines one instance of the JSBACH model
  !!
  TYPE, EXTENDS(t_Hsm) :: t_jsb_model
    ! From base type:
    ! CHARACTER(len=30)           :: name
    ! CLASS(t_State), POINTER     :: top     => NULL()    !< Topmost state
    ! CLASS(t_State), POINTER     :: current => NULL()    !< Current state
    ! CLASS(t_State), POINTER     :: next    => NULL()    !< Next state (non_null if transition taken)
    ! CLASS(t_State), POINTER     :: source  => NULL()    !< Source state during last transition
    ! TYPE(t_Action)              :: actions(MAX_ACTIONS) !< List of possible actions (events)
    ! Extensions:
    INTEGER                           :: id                !< Model id
    CHARACTER(len=10)                 :: shortname         !< Model short name
    CHARACTER(len=132)                :: description       !< Model description
    CHARACTER(len=filename_max)       :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model_config), POINTER :: config            !< Configuration
    INTEGER                           :: grid_id
    TYPE(t_jsb_process_p), POINTER    :: processes(:)
!!$    TYPE(t_hsm),              POINTER :: hsm
  CONTAINS
    ! From base type:
    ! PROCEDURE :: Get_current
    ! PROCEDURE :: Get_levels_to_LCA
    ! PROCEDURE :: Next_state
    ! PROCEDURE :: Get_state
    ! PROCEDURE :: Reset_state
    ! PROCEDURE :: Register_action
    ! PROCEDURE :: Print_actions
    ! PROCEDURE :: Process_message
    ! Extensions:
    PROCEDURE :: Configure_processes => Configure_model_processes
    PROCEDURE :: Get_shortname       => Get_model_shortname
    PROCEDURE :: Get_top_tile
    PROCEDURE :: Goto_next_tile
    PROCEDURE :: Set_top_tile
!!$    PROCEDURE :: Process_task
  END TYPE t_jsb_model

  !>
  !! Container CLASS for JSBACH model
  !!
  !! Contains pointer to an instance of type t_jsb_model in order to allow arrays
  !! of pointers to models.
  !!
  TYPE t_jsb_model_m
    TYPE(t_jsb_model), POINTER :: m
  END TYPE t_jsb_model_m

CONTAINS

  !>
  !! Constructor method for model instance.
  !!
  !! Subroutine to construct new instance of basic model structure.
  !!
  FUNCTION new_model(name, shortname, description, namelist_filename) RESULT(return_ptr)

    USE mo_jsb_process_factory, ONLY: max_no_of_processes

    CHARACTER(len=*),           INTENT(in) :: name
    CHARACTER(len=*),           INTENT(in) :: shortname
    CHARACTER(len=*),           INTENT(in) :: description
    CHARACTER(len=*),           INTENT(in) :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model), POINTER             :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:new_model'

    CALL message(TRIM(routine), 'starting construction of new JSBACH model instance.')

    ALLOCATE(return_ptr)

    CALL return_ptr%Init_hsm(name)

    return_ptr%shortname         = TRIM(shortname)
    return_ptr%description       = TRIM(description)
    return_ptr%namelist_filename = TRIM(namelist_filename)
    return_ptr%config => NULL()

!!$    return_ptr%Hsm => return_ptr%t_Hsm

    ALLOCATE(return_ptr%processes(max_no_of_processes))

    CALL message(TRIM(routine), 'construction of JSBACH model completed.')

  END FUNCTION new_model

  !
  ! State machine engine
  !
!!$  SUBROUTINE Process_task(this, task)
!!$
!!$    CLASS(t_jsb_model), INTENT(inout) :: this
!!$    CLASS(t_Message), POINTER       :: task
!!$
!!$    CALL this%Process_message(task)
!!$
!!$  END SUBROUTINE Process_task

  SUBROUTINE Configure_model_processes(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    INTEGER :: iproc

    DO iproc=1,SIZE(this%processes)
      ALLOCATE(this%processes(iproc)%p)
      CALL this%processes(iproc)%p%Configure(iproc, this%namelist_filename)
    END DO

  END SUBROUTINE Configure_model_processes

  FUNCTION Get_model_name(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%name

  END FUNCTION Get_model_name

  FUNCTION Get_model_shortname(this) RESULT(return_value)

    CLASS(t_jsb_model), INTENT(in) :: this
    CHARACTER(LEN=:), ALLOCATABLE :: return_value

    return_value = this%shortname

  END FUNCTION Get_model_shortname

  FUNCTION Get_top_tile(this) RESULT(return_ptr)

    USE mo_jsb_tile_class,   ONLY: t_jsb_tile

    CLASS(t_jsb_model), INTENT(in) :: this
    TYPE(t_jsb_tile), POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Get_model_top_tile'

    return_ptr => NULL()

    SELECT TYPE (top => this%top)
    TYPE IS (t_jsb_tile)
      return_ptr => top
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'tile is of type t_State, should be t_jsb_tile')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

  END FUNCTION Get_top_tile

  SUBROUTINE Set_top_tile(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Set_model_top_tile'

    CALL this%Reset_state()

  END SUBROUTINE Set_top_tile

  FUNCTION Goto_next_tile(this) RESULT(return_ptr)

    USE mo_jsb_tile_class,   ONLY: t_jsb_tile

    CLASS(t_jsb_model), INTENT(inout) :: this
    TYPE(t_jsb_tile), POINTER     :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Get_model_next_tile'

    next => this%Goto_next()
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    TYPE IS (t_jsb_tile)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

  END FUNCTION Goto_next_tile

  !>
  !! Destructor method for model instance.
  !!
  !! Subroutine to delete model instance from memory.
  !!
  SUBROUTINE delete_model(self)

    TYPE(t_jsb_model), POINTER, INTENT(inout) :: self

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:delete_model'

    CALL message(TRIM(routine), 'starting destruction of JSBACH model instance')

    DEALLOCATE(self)

    CALL message(TRIM(routine), 'destruction of JSBACH model instance completed.')

  END SUBROUTINE delete_model

  !>
  !! Print method for model instance.
  !!
  !! Subroutine to print description of model instance.
  !!
  SUBROUTINE print_model(self)

    USE mo_jsb_control,    ONLY: model_base_dir
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    TYPE(t_jsb_model), INTENT(in) :: self

    TYPE(t_jsb_grid), POINTER :: grid
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:print_model'

    CALL message(TRIM(routine), 'print JSBACH model description')

    WRITE(message_text,'(A,I2,A,A,A,A,A)') 'Model id: ', self%id, ', model name: ', TRIM(self%name), &
                                     ' (', TRIM(self%shortname), ')'
    CALL message('', message_text)
    CALL message('', '... '//TRIM(self%description))

    CALL message('', 'Namelist file: '//TRIM(model_base_dir)//'/'//TRIM(self%namelist_filename))

    grid => Get_grid(self%grid_id)
    IF (ASSOCIATED(grid)) CALL grid%Print

  END SUBROUTINE print_model

#endif
END MODULE mo_jsb_model_class
