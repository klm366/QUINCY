!>
!! @brief Contains process task class
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2015-09-01)
!!
MODULE mo_jsb_task_class
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: message, finish

  USE mo_hsm_class,       ONLY: t_Message
  USE mo_jsb_tile_class,  ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class, ONLY: t_jsb_config_p

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process_task, t_jsb_process_task_p, t_Task_queue, t_jsb_task_options

  TYPE t_jsb_task_options
    REAL(wp) :: dtime   = -1._wp
    REAL(wp) :: steplen = -1._wp
    INTEGER  :: grid_id = -1
    INTEGER  :: iblk    = -1  ! Number of current block (chunk)
    INTEGER  :: ics     = -1  ! Index of chunk start
    INTEGER  :: ice     = -1  ! Index of chunk end
    INTEGER  :: nc      = -1  ! Length of chunk
  END TYPE t_jsb_task_options

  TYPE, ABSTRACT, EXTENDS(t_Message) :: t_jsb_process_task
    ! From base type:
    ! CHARACTER(len=50) :: name
    ! TYPE(t_Action)    :: action
    ! Extensions:
    INTEGER :: process_id     = 0
    INTEGER :: owner_model_id = 0
!!$    PROCEDURE(Integrate_task), POINTER, NOPASS :: Integrate
!!$    PROCEDURE(Aggregate_task), POINTER, NOPASS :: Aggregate
  CONTAINS
    PROCEDURE(Integrate_task), DEFERRED, NOPASS :: Integrate
    PROCEDURE(Aggregate_task), DEFERRED, NOPASS :: Aggregate
    PROCEDURE :: Do_it
  END TYPE t_jsb_process_task

  ABSTRACT INTERFACE
    SUBROUTINE Integrate_task(tile, config, options)
      IMPORT :: t_jsb_tile_abstract, t_jsb_config_p, t_jsb_task_options
      CLASS(t_jsb_tile_abstract),         INTENT(inout) :: tile
      CLASS(t_jsb_config_p),     INTENT(in) :: config(:)
      TYPE(t_jsb_task_options),  INTENT(in) :: options
    END SUBROUTINE Integrate_task
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE Aggregate_task(tile, config, options)
      IMPORT :: t_jsb_tile_abstract, t_jsb_config_p, t_jsb_task_options
      CLASS(t_jsb_tile_abstract),         INTENT(inout) :: tile
      CLASS(t_jsb_config_p),     INTENT(in) :: config(:)
      TYPE(t_jsb_task_options),  INTENT(in) :: options
    END SUBROUTINE Aggregate_task
  END INTERFACE

  TYPE t_jsb_process_task_p
    CLASS(t_jsb_process_task), POINTER :: p
  END TYPE t_jsb_process_task_p

  TYPE t_Task_queue
    CLASS(t_jsb_process_task), POINTER :: task => NULL()
    TYPE(t_Task_queue),    POINTER :: next => NULL()
  CONTAINS
    PROCEDURE :: Append => Append_task
  END TYPE t_Task_queue

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_task_class'

CONTAINS

!!$  SUBROUTINE Do_it(task, mem, config)
  SUBROUTINE Do_it(this, tile, config, options)

    CLASS(t_jsb_process_task), INTENT(in) :: this
    CLASS(t_jsb_tile_abstract),     INTENT(inout) :: tile
    CLASS(t_jsb_config_p),     INTENT(in) :: config(:)
    TYPE(t_jsb_task_options),  INTENT(in) :: options

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Do_it'

    ! If process the task belongs to is not active on this tile, do nothing.
    IF (.NOT. tile%Is_process_active(this%process_id)) RETURN

    SELECT TYPE (this)
    CLASS IS (t_jsb_process_task)

      IF (.NOT. ASSOCIATED(this%action)) &
        & CALL finish(TRIM(routine), 'Action not set for task '//TRIM(this%name))

!!$      IF (.NOT. ASSOCIATED(task%Integrate)) &
!!$        & CALL finish(TRIM(routine), 'Procedure pointer for Integrate not set')
!!$
!!$      IF (.NOT. ASSOCIATED(task%Aggregate)) &
!!$        & CALL finish(TRIM(routine), 'Procedure pointer for Aggregate not set')
!!$
      SELECT CASE(TRIM(this%action%name))
      CASE ('INTEGRATE')
        CALL this%Integrate(tile, config, options)
      CASE ('AGGREGATE')
        CALL this%Aggregate(tile, config, options)
      CASE DEFAULT
        CALL finish(TRIM(routine), 'Unknown action for task'//TRIM(this%name))
      END SELECT

    END SELECT

  END SUBROUTINE Do_it

  SUBROUTINE Append_task(this, task)

    CLASS(t_Task_queue), TARGET, INTENT(inout)  :: this
    CLASS(t_jsb_process_task), TARGET :: task

    CLASS(t_Task_queue), POINTER :: current, next

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Append_task'

    IF (.NOT. ASSOCIATED(this%task)) THEN
      ! If queue is empty just add task as first element
      this%task => task
    ELSE
      ! Else find last task in queue
      current => this
      DO WHILE (ASSOCIATED(current%next))
        current => current%next
      END DO

      ! Append new element
      ALLOCATE(next)
      next%task => task
      current%next => next
    END IF

    CALL message(TRIM(routine), 'Appended task "'//TRIM(task%name)//'" to model queue')

  END SUBROUTINE Append_task

#endif
END MODULE mo_jsb_task_class
