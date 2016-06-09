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

  USE mo_exception,         ONLY: message_text, message, finish
  USE mo_kind,              ONLY: wp
  USE mo_io_units,          ONLY: filename_max

  USE mo_hsm_class,         ONLY: t_Hsm, t_State, t_Message, t_Action
  USE mo_jsb_config_class,  ONLY: t_jsb_model_config, t_jsb_config_p
  USE mo_jsb_process_class, ONLY: t_jsb_process, t_jsb_process_p
  USE mo_jsb_task_class,    ONLY: t_jsb_process_task, t_Task_queue, t_jsb_task_msg, t_jsb_task_options
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract
  USE mo_jsb_lctlib_class,  ONLY: t_lctlib

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
    INTEGER                            :: id                !< Model id
    CHARACTER(len=10)                  :: shortname         !< Model short name
    CHARACTER(len=132)                 :: description       !< Model description
    CHARACTER(len=filename_max)        :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model_config),  POINTER :: config            !< Configuration
    TYPE(t_jsb_task_options)           :: options           !< Options needed by tasks
    INTEGER                            :: grid_id
    TYPE(t_jsb_process_p),     POINTER :: processes(:)
    TYPE(t_Task_queue),        POINTER :: queue        => NULL()
    CLASS(t_jsb_process_task), POINTER :: current_task => NULL()
    TYPE(t_lctlib)                     :: lctlib
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
    PROCEDURE :: Get_next_tile
    PROCEDURE :: Reset_tiles
    PROCEDURE :: Get_tile
    PROCEDURE :: Get_process
    PROCEDURE :: Queue_task
!!$    PROCEDURE :: Get_process_configs
!!$    PROCEDURE :: Process_task
    PROCEDURE :: Run_tasks
    PROCEDURE :: Set_options
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

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_class'

CONTAINS

  !>
  !! Constructor method for model instance.
  !!
  !! Subroutine to construct new instance of basic model structure.
  !!
  FUNCTION new_model(id, name, shortname, description, namelist_filename) RESULT(return_ptr)

!!$    USE mo_jsb_process_factory, ONLY: max_no_of_processes, Create_process

    INTEGER,                    INTENT(in) :: id
    CHARACTER(len=*),           INTENT(in) :: name
    CHARACTER(len=*),           INTENT(in) :: shortname
    CHARACTER(len=*),           INTENT(in) :: description
    CHARACTER(len=*),           INTENT(in) :: namelist_filename !< Model namelist file
    TYPE(t_jsb_model), POINTER             :: return_ptr

    INTEGER :: iproc

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:new_model'

    CALL message(TRIM(routine), 'starting construction of new JSBACH model instance.')

    ALLOCATE(return_ptr)

    CALL return_ptr%Init_hsm(name)

    CALL return_ptr%Register_action(t_Action("INTEGRATE"))
    CALL return_ptr%Register_action(t_Action("AGGREGATE"))

    return_ptr%id                = id
    return_ptr%shortname         = TRIM(shortname)
    return_ptr%description       = TRIM(description)
    return_ptr%namelist_filename = TRIM(namelist_filename)
    return_ptr%config => NULL()

!!$    return_ptr%Hsm => return_ptr%t_Hsm

!!$    ALLOCATE(return_ptr%processes(max_no_of_processes))
!!$    DO iproc=1,max_no_of_processes
!!$      return_ptr%processes(iproc)%p => Create_process(iproc, return_ptr%id, TRIM(namelist_filename))
!!$    END DO

    ALLOCATE(return_ptr%queue)

    CALL message(TRIM(routine), 'construction of JSBACH model completed.')

  END FUNCTION new_model

  FUNCTION Get_process(this, iproc) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(in) :: this
    INTEGER, INTENT(in) :: iproc
    CLASS(t_jsb_process), POINTER :: return_ptr

    return_ptr => this%processes(iproc)%p

  END FUNCTION Get_process

  SUBROUTINE Queue_task(this, iproc, name)

    CLASS(t_jsb_model),        INTENT(inout) :: this
    INTEGER,                   INTENT(in)    :: iproc
    CHARACTER(LEN=*),          INTENT(in)    :: name

    CALL this%queue%Append(this%processes(iproc)%p%Get_task(TRIM(name)))

  END SUBROUTINE Queue_task

  SUBROUTINE Configure_model_processes(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    INTEGER :: iproc

    DO iproc=1,SIZE(this%processes)
      IF (ASSOCIATED(this%processes(iproc)%p)) &
        & CALL this%processes(iproc)%p%Configure(this%namelist_filename)
    END DO

  END SUBROUTINE Configure_model_processes

!!$  FUNCTION Get_process_configs(this) RESULT(return_ptr)
!!$
!!$    CLASS(t_jsb_model), INTENT(in) :: this
!!$    CLASS(t_jsb_config_p), POINTER :: return_ptr(:)
!!$
!!$    INTEGER :: iproc
!!$
!!$    ALLOCATE(return_ptr(SIZE(this%processes)))
!!$    DO iproc=1,SIZE(this%processes)
!!$      IF (ASSOCIATED(this%processes(iproc)%p)) THEN
!!$        return_ptr(iproc)%p => this%processes(iproc)%p%config
!!$      ELSE
!!$        return_ptr(iproc)%p => NULL()
!!$      END IF
!!$    END DO
!!$
!!$  END FUNCTION Get_process_configs
  !
  ! State machine engine
  !
!!$  ! Just a stub to the real engine hsm%Process_message
!!$  SUBROUTINE Process_task(this, task)
!!$
!!$    CLASS(t_jsb_model),    INTENT(inout) :: this
!!$    CLASS(t_jsb_task_msg), INTENT(in)    :: msg
!!$
!!$    CALL this%Process_message(msg)
!!$
!!$  END SUBROUTINE Process_task
!!$

  SUBROUTINE Run_tasks(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    CLASS(t_Task_queue), POINTER :: current_queue_element
    TYPE(t_jsb_task_msg) :: msg

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Run_tasks'

    current_queue_element => this%queue
    DO WHILE (ASSOCIATED(current_queue_element))

      ASSOCIATE(task => current_queue_element%task)
        IF (this%processes(task%process_id)%p%config%active) THEN  ! Only if process this task belongs to is active
          this%current_task => task
          CALL this%Start() ! Start state machine
          msg%action => this%Get_action('INTEGRATE')
          CALL this%Process_message(msg)
        END IF
      END ASSOCIATE

      current_queue_element => current_queue_element%next

    END DO

  END SUBROUTINE Run_tasks

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

    CLASS(t_jsb_model), INTENT(in) :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Get_model_top_tile'

    return_ptr => NULL()

    SELECT TYPE (top => this%top)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => top
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'tile is of type t_State, should be t_jsb_tile')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

  END FUNCTION Get_top_tile

  SUBROUTINE Reset_tiles(this)

    CLASS(t_jsb_model), INTENT(inout) :: this

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Set_model_top_tile'

    CALL this%Reset()

  END SUBROUTINE Reset_tiles

  FUNCTION Goto_next_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(inout)  :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Goto_next_tile'

    next => this%Goto_next()
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

  END FUNCTION Goto_next_tile

  FUNCTION Get_next_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(inout)  :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_model_class:Get_next_tile'

    next => this%Get_next()
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
      CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
      CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile')
      CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find next tile')

  END FUNCTION Get_next_tile

  FUNCTION Get_tile(this, tile_id) RESULT(return_ptr)

    CLASS(t_jsb_model), INTENT(in) :: this
    INTEGER,            INTENT(in) :: tile_id(:)
    CLASS(t_jsb_tile_abstract),  POINTER    :: return_ptr

    CLASS(t_State), POINTER :: state

    state => this%Get_state(tile_id)

    SELECT TYPE (state)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => state
    END SELECT

    NULLIFY(state)

  END FUNCTION Get_tile

  SUBROUTINE Set_options(this, dtime, steplen, iblk, ics, ice, nc)

    USE mo_jsb_utils_iface, ONLY: assign_if_present

    CLASS(t_jsb_model), INTENT(inout) :: this
    REAL(wp), OPTIONAL, INTENT(in)    :: dtime
    REAL(wp), OPTIONAL, INTENT(in)    :: steplen
    INTEGER,  OPTIONAL, INTENT(in)    :: iblk
    INTEGER,  OPTIONAL, INTENT(in)    :: ics
    INTEGER,  OPTIONAL, INTENT(in)    :: ice
    INTEGER,  OPTIONAL, INTENT(in)    :: nc

    CALL assign_if_present(this%options%dtime,   dtime)
    CALL assign_if_present(this%options%steplen, steplen)
    CALL assign_if_present(this%options%iblk,    iblk)
    CALL assign_if_present(this%options%ics,     ics)
    CALL assign_if_present(this%options%ice,     ice)
    CALL assign_if_present(this%options%nc,      nc)

  END SUBROUTINE Set_options

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
