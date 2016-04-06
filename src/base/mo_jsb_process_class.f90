!>
!! @brief Contains process class
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2015-09-01)
!!
MODULE mo_jsb_process_class
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message, finish

  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_process_task_p
  USE mo_jsb_config_class,    ONLY: t_jsb_config
  USE mo_jsb_process_factory, ONLY: Get_process_name

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process, t_jsb_process_p
  PUBLIC :: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_, Get_action_name

  TYPE t_jsb_process
    INTEGER                      :: id = 0
    CLASS(t_jsb_config), POINTER :: config => NULL()
    TYPE(t_jsb_process_task_p), ALLOCATABLE :: task_list(:)
    INTEGER                      :: owner_model_id
  CONTAINS
    PROCEDURE :: Configure => Configure_process
    PROCEDURE :: Register_task
    PROCEDURE :: Register_tasks
    PROCEDURE :: Get_task
  END TYPE t_jsb_process

!!$  INTERFACE t_jsb_process
!!$    PROCEDURE Create_jsb_process
!!$  END INTERFACE t_jsb_process

  TYPE t_jsb_process_p
    CLASS(t_jsb_process), POINTER :: p
  END TYPE t_jsb_process_p

  ! Process actions
  ! It is important that AGGREGATE comes before the other actions
  ENUM, BIND(C)
    ENUMERATOR :: SKIP_=0, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_class'

CONTAINS

!!$  FUNCTION Create_jsb_process(iproc, model_id) RESULT(return_ptr)
!!$
!!$    !USE mo_jsb_process_factory, ONLY: Create_process
!!$
!!$    CLASS(t_jsb_process), POINTER :: return_ptr
!!$    INTEGER,           INTENT(in) :: iproc
!!$    INTEGER,           INTENT(in) :: model_id
!!$
!!$    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Create_process'
!!$
!!$    ALLOCATE(return_ptr)
!!$    return_ptr%id             = iproc
!!$    return_ptr%owner_model_id = model_id
!!$
!!$    CALL Register_process_tasks(iproc, model_id)
!!$    CALL return_ptr%Register_tasks(model_id)
!!$
!!$  END FUNCTION Create_jsb_process

  SUBROUTINE Register_task(this, task)

    CLASS(t_jsb_process), INTENT(inout) :: this
    CLASS(t_jsb_process_task), TARGET, INTENT(in) :: task

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Register_task'

    IF (this%id /= task%process_id) &
      CALL finish(TRIM(routine), 'Task is for another process')

    IF (.NOT. ALLOCATED(this%task_list)) THEN
      !ALLOCATE(this%task_list(1))
      this%task_list = (/ t_jsb_process_task_p(task) /)
    ELSE
      this%task_list = (/ this%task_list, t_jsb_process_task_p(task) /)
    END IF

    CALL message(TRIM(routine), 'Registered task "'//TRIM(task%name)// &
      &                         '" for process "'//TRIM(Get_process_name(this%id))//'".')

  END SUBROUTINE Register_task

  ! Dummy subroutine to be overwritten
  SUBROUTINE Register_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

  END SUBROUTINE Register_tasks

  FUNCTION Get_task(this, name) RESULT(return_ptr)

    CLASS(t_jsb_process), INTENT(inout) :: this
    CHARACTER(LEN=*),     INTENT(in)    :: name
    CLASS(t_jsb_process_task), POINTER :: return_ptr

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_task'

    return_ptr => NULL()
    DO i=1,SIZE(this%task_list)
      !IF (.NOT. ASSOCIATED())
      IF (TRIM(this%task_list(i)%p%name) == TRIM(name) ) THEN
        return_ptr => this%task_list(i)%p
        EXIT
      END IF
    END DO

    IF (.NOT.(ASSOCIATED(return_ptr))) &
      & CALL finish(TRIM(routine), 'Task "'//TRIM(name)//'" not found for process "'// &
      &                            TRIM(Get_process_name(this%id))//'".')

  END FUNCTION Get_task

  SUBROUTINE Configure_process(this, namelist_filename)

    !USE mo_jsb_process_factory, ONLY: Create_process_config

    CLASS(t_jsb_process), INTENT(inout) :: this
    CHARACTER(LEN=*),     INTENT(in)    :: namelist_filename

    CLASS(t_jsb_config), POINTER :: config

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Configure_process'

    !config  => Create_process_config(this%id, TRIM(namelist_filename))

    IF (.NOT. ASSOCIATED(this%config)) RETURN

    !SELECT TYPE (config)
    !CLASS IS (t_jsb_config)
    !  this%config => config
      CALL this%config%Init()
    !END SELECT

    !NULLIFY(config)

  END SUBROUTINE Configure_process

  FUNCTION Get_action_name(iact) RESULT(return_value)

    INTEGER, INTENT(in) :: iact
    CHARACTER(len=:), ALLOCATABLE :: return_value

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Get_action_name'

    SELECT CASE(iact)
    CASE (SKIP_)
      return_value = 'SKIP'
    CASE (AGGREGATE_)
      return_value = 'AGGREGATE'
    CASE (ON_LEAFS_)
      return_value = 'ON_LEAFS'
    CASE(ON_TILE_)
      return_value = 'ON_TILE'
    CASE (ON_SUBTREE_)
      return_value = 'ON_SUBTREE'
    CASE (INHERIT_)
      return_value = 'INHERIT'
    END SELECT

  END FUNCTION Get_action_name

#endif
END MODULE mo_jsb_process_class
