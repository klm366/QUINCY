!> Initialize insect disturbance variables (QUINCY)
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### Interface to the insect disturbance processes
!>
MODULE mo_q_dinsect_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, Q_DINSECT_
  USE mo_q_dinsect_util,              ONLY: reset_dinsect_fluxes
  USE mo_q_dinsect_main,              ONLY: update_dinsect_effects, update_dinsect_barkbeetle, update_dinsect_defoliator



  ! Use of process configurations (t_PROC_config)
  dsl4jsb_Use_config(Q_DINSECT_)

  ! Use of process memories (t_PROC_memory)
  dsl4jsb_Use_memory(Q_DINSECT_)


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_dinsect_tasks
  PUBLIC :: reset_dinsect_fluxes, update_dinsect_effects, update_dinsect_barkbeetle, update_dinsect_defoliator

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_dinsect_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_dinsect_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => reset_dinsect_fluxes      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_reset_dinsect_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_dinsect_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_dinsect_effects task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_dinsect_effects
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dinsect_effects      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dinsect_effects        !< Aggregates computed task variables
  END TYPE tsk_update_dinsect_effects

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_dinsect_barkbeetle task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_dinsect_barkbeetle
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dinsect_barkbeetle      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dinsect_barkbeetle        !< Aggregates computed task variables
  END TYPE tsk_update_dinsect_barkbeetle


  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_dinsect_defoliator task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_dinsect_defoliator
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dinsect_defoliator      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dinsect_defoliator        !< Aggregates computed task variables
  END TYPE tsk_update_dinsect_defoliator


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_dinsect_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_dinsect_fluxes
    PROCEDURE Create_task_reset_dinsect_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_dinsect_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_dinsect_effects task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  INTERFACE tsk_update_dinsect_effects
    PROCEDURE Create_task_update_dinsect_effects         !< Constructor function for task
  END INTERFACE tsk_update_dinsect_effects


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_dinsect_barkbeetle task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  INTERFACE tsk_update_dinsect_barkbeetle
    PROCEDURE Create_task_update_dinsect_barkbeetle         !< Constructor function for task
  END INTERFACE tsk_update_dinsect_barkbeetle


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_dinsect_defoliator task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  INTERFACE tsk_update_dinsect_defoliator
    PROCEDURE Create_task_update_dinsect_defoliator         !< Constructor function for task
  END INTERFACE tsk_update_dinsect_defoliator

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks:_dinsect_ process
  !!
  !-------------------------------------------------------------------------------------------------------- !
  SUBROUTINE Register_dinsect_tasks(this, model_id)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_dinsect_fluxes(model_id))
    CALL this%Register_task(tsk_update_dinsect_effects(model_id))
    CALL this%Register_task(tsk_update_dinsect_barkbeetle(model_id))
    CALL this%Register_task(tsk_update_dinsect_defoliator(model_id))

  END SUBROUTINE Register_dinsect_tasks

 !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_dinsect_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_dinsect_fluxes(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_reset_dinsect_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_dinsect_fluxes', process_id=Q_DINSECT_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_dinsect_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_dinsect_effects task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  FUNCTION Create_task_update_dinsect_effects(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_update_dinsect_effects::return_ptr)
    CALL return_ptr%Construct(name='update_dinsect_effects', process_id=Q_DINSECT_, owner_model_id=model_id)

  END FUNCTION Create_task_update_dinsect_effects


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_dinsect_barkbeetle task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  FUNCTION Create_task_update_dinsect_barkbeetle(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_update_dinsect_barkbeetle::return_ptr)
    CALL return_ptr%Construct(name='update_dinsect_barkbeetle', process_id=Q_DINSECT_, owner_model_id=model_id)

  END FUNCTION Create_task_update_dinsect_barkbeetle


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_dinsect_defoliator task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  FUNCTION Create_task_update_dinsect_defoliator(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_update_dinsect_defoliator::return_ptr)
    CALL return_ptr%Construct(name='update_dinsect_defoliator', process_id=Q_DINSECT_, owner_model_id=model_id)

  END FUNCTION Create_task_update_dinsect_defoliator

 !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": reset_dinsect_fluxes task
  !!  actually NOT using this routine, i.e., no var is aggregated here
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_reset_dinsect_fluxes(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_reset_dinsect_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_DINSECT_)
    dsl4jsb_Def_memory(Q_DINSECT_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(Q_DINSECT_)
    ! Get process memories
    dsl4jsb_Get_memory(Q_DINSECT_)
    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen pheno Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_reset_dinsect_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_dinsect_effects task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  SUBROUTINE aggr_dinsect_effects(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dinsect_effects'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_DINSECT_)
    dsl4jsb_Def_memory(Q_DINSECT_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(Q_DINSECT_)
    ! Get process memories
    dsl4jsb_Get_memory(Q_DINSECT_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen pheno Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_dinsect_effects

  !-----------------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_dinsect_effects task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  SUBROUTINE aggr_dinsect_barkbeetle(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dinsect_barkbeetle'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_DINSECT_)
    dsl4jsb_Def_memory(Q_DINSECT_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(Q_DINSECT_)
    ! Get process memories
    dsl4jsb_Get_memory(Q_DINSECT_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen pheno Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_dinsect_barkbeetle


  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_dinsect_effects task
  !!
  !-------------------------------------------------------------------------------------------------------- !
  SUBROUTINE aggr_dinsect_defoliatior(tile, options)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dinsect_defoliator'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(Q_DINSECT_)
    dsl4jsb_Def_memory(Q_DINSECT_)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model => Get_model(tile%owner_model_id)
    ! Get process config
    dsl4jsb_Get_config(Q_DINSECT_)
    ! Get process memories
    dsl4jsb_Get_memory(Q_DINSECT_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
      ! do not aggregate after this routine
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen pheno Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_dinsect_defoliator
#endif
END MODULE mo_q_dinsect_interface
