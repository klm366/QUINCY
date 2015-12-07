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

  USE mo_hsm_class, ONLY: t_Message, t_hsm
  USE mo_jsb_config_class, ONLY: t_jsb_config
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_process, t_jsb_process_p, t_jsb_process_task
  PUBLIC :: SKIP_, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_, Get_action_name

  TYPE t_jsb_process
    INTEGER                      :: id = 0
    CLASS(t_jsb_config), POINTER :: config => NULL()
  CONTAINS
    PROCEDURE :: Configure => Configure_process
!!$    PROCEDURE :: Register_task
  END TYPE t_jsb_process

  TYPE t_jsb_process_p
    TYPE(t_jsb_process), POINTER :: p
  END TYPE t_jsb_process_p

  TYPE, ABSTRACT, EXTENDS(t_Message) :: t_jsb_process_task
    ! From base type:
    ! CHARACTER(len=50) :: name
    ! TYPE(t_Action)    :: action
!!$    CLASS(t_Message), POINTER :: message
    ! Extensions:
    INTEGER :: process_id
  CONTAINS
    PROCEDURE(Integrate_task), DEFERRED, NOPASS :: Integrate
    PROCEDURE(Aggregate_task), DEFERRED, NOPASS :: Aggregate
  END TYPE t_jsb_process_task

  ABSTRACT INTERFACE
    SUBROUTINE Integrate_task(hsm, iblk, ics, ice, nc)
      IMPORT :: t_hsm
      CLASS(t_hsm),              INTENT(inout) :: hsm
      INTEGER,                   INTENT(in)    :: iblk
      INTEGER,                   INTENT(in)    :: ics
      INTEGER,                   INTENT(in)    :: ice
      INTEGER, OPTIONAL,         INTENT(in)    :: nc
    END SUBROUTINE Integrate_task
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE Aggregate_task(hsm, iblk, ics, ice, nc)
      IMPORT :: t_hsm
      CLASS(t_hsm),              INTENT(inout) :: hsm
      INTEGER,                   INTENT(in)    :: iblk
      INTEGER,                   INTENT(in)    :: ics
      INTEGER,                   INTENT(in)    :: ice
      INTEGER, OPTIONAL,         INTENT(in)    :: nc
    END SUBROUTINE Aggregate_task
  END INTERFACE

  ! Process actions
  ! It is important that AGGREGATE comes before the other actions
  ENUM, BIND(C)
    ENUMERATOR :: SKIP_=0, AGGREGATE_, ON_LEAFS_, ON_TILE_, ON_SUBTREE_, INHERIT_
  END ENUM

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_process_class'

CONTAINS

  SUBROUTINE Configure_process(this, iproc, namelist_filename)

    USE mo_jsb_process_factory, ONLY: Create_process_config

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: iproc
    CHARACTER(LEN=*),     INTENT(in)    :: namelist_filename

    CLASS(t_jsb_config), POINTER :: config

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Configure_process'

    this%id       = iproc

    config  => Create_process_config(iproc, TRIM(namelist_filename))

    IF (.NOT. ASSOCIATED(config)) RETURN

    SELECT TYPE (config)
    CLASS IS (t_jsb_config)
      this%config => config
      CALL this%config%Init()
    END SELECT

    NULLIFY(config)

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
