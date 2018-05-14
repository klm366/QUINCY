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
  USE mo_jsb_control,         ONLY: debug_on

  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_process_task_p
  USE mo_jsb_config_class,    ONLY: t_jsb_config

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process, t_jsb_process_p, Get_process_name
  PUBLIC :: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, HD_, RAD_, ASSIMI_, PHENO_, CARBON_ , DISTURB_, FUEL_
  PUBLIC :: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_, Get_action_name

  TYPE t_jsb_process
    INTEGER                      :: id = 0
    CHARACTER(LEN=10)            :: name = ''
    CLASS(t_jsb_config), POINTER :: config => NULL()
    TYPE(t_jsb_process_task_p), ALLOCATABLE :: task_list(:)
    INTEGER                      :: owner_model_id
  CONTAINS
    PROCEDURE :: Configure => Configure_process
    PROCEDURE :: Register_task
    PROCEDURE :: Get_task
  END TYPE t_jsb_process

  ! The constructor for t_jsb_process is Create_process in mo_jsb_process_factory.

  TYPE t_jsb_process_p
    CLASS(t_jsb_process), POINTER :: p
  END TYPE t_jsb_process_p

  ENUM, BIND(C)
    ENUMERATOR :: A2L_=1, L2A_, SEB_, TURB_, SSE_, HYDRO_, HD_, RAD_, ASSIMI_, PHENO_,  CARBON_ , DISTURB_, FUEL_
  END ENUM

  ! Process actions
  ! It is important that AGGREGATE comes before the other actions
  ENUM, BIND(C)
    ENUMERATOR :: SKIP_=0, INHERIT_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_class'

CONTAINS

  FUNCTION Get_process_name(iproc) RESULT(return_value)

    INTEGER, INTENT(in) :: iproc
    CHARACTER(len=:), ALLOCATABLE :: return_value

    SELECT CASE(iproc)
    CASE (A2L_)
      return_value = 'a2l'
    CASE (L2A_)
      return_value = 'l2a'
    CASE (SEB_)
      return_value = 'seb'
    CASE (TURB_)
      return_value = 'turb'
    CASE (SSE_)
      return_value = 'sse'
    CASE (HYDRO_)
      return_value = 'hydro'
    CASE (RAD_)
      return_value = 'rad'
    CASE (ASSIMI_)
      return_value = 'assimi'
    CASE (PHENO_)
      return_value = 'pheno'
#ifndef __NO_JSBACH_HD__
    CASE (HD_)
      return_value = 'hd'
#endif
    CASE (CARBON_)
      return_value = 'carbon'
    CASE (DISTURB_)
      return_value = 'disturb'
    CASE (FUEL_)
      return_value = 'fuel'
    CASE DEFAULT
      return_value = ''
    END SELECT

  END FUNCTION Get_process_name

  SUBROUTINE Register_task(this, task)

    CLASS(t_jsb_process),     INTENT(inout) :: this
    CLASS(t_jsb_process_task), TARGET, INTENT(in)    :: task

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Register_task'

    IF (this%id /= task%process_id) &
      CALL finish(TRIM(routine), 'Task is for another process')

    IF (.NOT. ALLOCATED(this%task_list)) THEN
      this%task_list = (/ t_jsb_process_task_p(task) /)
    ELSE
      this%task_list = (/ this%task_list, t_jsb_process_task_p(task) /)
    END IF

    IF (debug_on()) CALL message(TRIM(routine), 'Registered task "'//TRIM(task%name)// &
      &                                         '" for process "'//TRIM(this%name)//'".')

  END SUBROUTINE Register_task

  FUNCTION Get_task(this, name) RESULT(return_ptr)

    CLASS(t_jsb_process), INTENT(inout) :: this
    CHARACTER(LEN=*),              INTENT(in)    :: name
    CLASS(t_jsb_process_task), POINTER           :: return_ptr

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
      &                            TRIM(this%name)//'".')

  END FUNCTION Get_task

  SUBROUTINE Configure_process(this, namelist_filename)

    CLASS(t_jsb_process), INTENT(inout) :: this
    CHARACTER(LEN=*),     INTENT(in)    :: namelist_filename

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Configure_process'

    IF (.NOT. ASSOCIATED(this%config)) RETURN

    CALL this%config%Init()

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
