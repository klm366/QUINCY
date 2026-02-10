!> QUINCY fire-disturbance process interface
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
!>#### definition and init of tasks for the fire-disturbance process incl. update and aggregate routines
!>
MODULE mo_dfire_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, &
    &                                       DFIRE_
  USE mo_dfire_sf_fuel,               ONLY: update_dfire_sf_fuel
  USE mo_dfire_sf_ignition,           ONLY: update_dfire_sf_ignition
  USE mo_dfire_sf_spread,             ONLY: update_dfire_sf_spread
  USE mo_dfire_sf_effects,            ONLY: update_dfire_sf_effects
  USE mo_dfire_sf_feedback,           ONLY: update_dfire_sf_feedback

  dsl4jsb_Use_memory(DFIRE_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_dfire_tasks
  PUBLIC :: update_dfire_sf_fuel, update_dfire_sf_ignition, &
    &       update_dfire_sf_spread, update_dfire_sf_effects, &
    &       update_dfire_sf_feedback

  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_fuel task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_dfire_sf_fuel
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dfire_sf_fuel      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dfire_sf_fuel        !< Aggregates computed task variables
  END TYPE tsk_dfire_sf_fuel
  
  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_fuel task
  !>
  INTERFACE tsk_dfire_sf_fuel
    PROCEDURE Create_task_dfire_sf_fuel         !< Constructor function for task
  END INTERFACE tsk_dfire_sf_fuel

  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_ignition task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_dfire_sf_ignition
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dfire_sf_ignition  !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dfire_sf_ignition    !< Aggregates computed task variables
  END TYPE tsk_dfire_sf_ignition

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_ignition task
  !>
  INTERFACE tsk_dfire_sf_ignition
    PROCEDURE Create_task_dfire_sf_ignition         !< Constructor function for task
  END INTERFACE tsk_dfire_sf_ignition

  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_spread task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_dfire_sf_spread
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dfire_sf_spread    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dfire_sf_spread      !< Aggregates computed task variables
  END TYPE tsk_dfire_sf_spread

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_spread task
  !>
  INTERFACE tsk_dfire_sf_spread
    PROCEDURE Create_task_dfire_sf_spread         !< Constructor function for task
  END INTERFACE tsk_dfire_sf_spread

  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_effects task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_dfire_sf_effects
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dfire_sf_effects     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dfire_sf_effects       !< Aggregates computed task variables
  END TYPE tsk_dfire_sf_effects

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_effects task
  !>
  INTERFACE tsk_dfire_sf_effects
    PROCEDURE Create_task_dfire_sf_effects         !< Constructor function for task
  END INTERFACE tsk_dfire_sf_effects

  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_feedback task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_dfire_sf_feedback
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_dfire_sf_feedback    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_dfire_sf_feedback      !< Aggregates computed task variables
  END TYPE tsk_dfire_sf_feedback

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_feedback task
  !>
  INTERFACE tsk_dfire_sf_feedback
    PROCEDURE Create_task_dfire_sf_feedback         !< Constructor function for task
  END INTERFACE tsk_dfire_sf_feedback

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_interface'

CONTAINS
  
  ! ======================================================================================================= !
  !>
  !> Register tasks: DFIRE_ process
  !>
  SUBROUTINE Register_dfire_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_dfire_sf_fuel(model_id))
    CALL this%Register_task(tsk_dfire_sf_ignition(model_id))
    CALL this%Register_task(tsk_dfire_sf_spread(model_id))
    CALL this%Register_task(tsk_dfire_sf_effects(model_id))
    CALL this%Register_task(tsk_dfire_sf_feedback(model_id))
  END SUBROUTINE Register_dfire_tasks

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_fuel task
  !>
  FUNCTION Create_task_dfire_sf_fuel(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dfire_sf_ignition::return_ptr)
    CALL return_ptr%Construct(name='dfire_sf_fuel', process_id=DFIRE_, owner_model_id=model_id)
  END FUNCTION Create_task_dfire_sf_fuel

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_ignition task
  !>
  FUNCTION Create_task_dfire_sf_ignition(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dfire_sf_ignition::return_ptr)
    CALL return_ptr%Construct(name='dfire_sf_ignition', process_id=DFIRE_, owner_model_id=model_id)
  END FUNCTION Create_task_dfire_sf_ignition

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_spread task
  !>
  FUNCTION Create_task_dfire_sf_spread(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dfire_sf_spread::return_ptr)
    CALL return_ptr%Construct(name='dfire_sf_spread', process_id=DFIRE_, owner_model_id=model_id)
  END FUNCTION Create_task_dfire_sf_spread

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_effects task
  !>
  FUNCTION Create_task_dfire_sf_effects(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dfire_sf_effects::return_ptr)
    CALL return_ptr%Construct(name='dfire_sf_effects', process_id=DFIRE_, owner_model_id=model_id)
  END FUNCTION Create_task_dfire_sf_effects

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_feedback task
  !>
  FUNCTION Create_task_dfire_sf_feedback(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_dfire_sf_feedback::return_ptr)
    CALL return_ptr%Construct(name='dfire_sf_feedback', process_id=DFIRE_, owner_model_id=model_id)
  END FUNCTION Create_task_dfire_sf_feedback

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_sf_fuel task
  !>
  SUBROUTINE aggr_dfire_sf_fuel(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dfire_sf_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_dfire_sf_fuel

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_sf_ignition task
  !>
  SUBROUTINE aggr_dfire_sf_ignition(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dfire_sf_ignition'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_dfire_sf_ignition

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_sf_spread task
  !>
  SUBROUTINE aggr_dfire_sf_spread(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dfire_sf_spread'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_dfire_sf_spread

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_sf_effects task
  !>
  SUBROUTINE aggr_dfire_sf_effects(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dfire_sf_effects'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_dfire_sf_effects

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_sf_feedback task
  !>
  SUBROUTINE aggr_dfire_sf_feedback(tile, options)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_dfire_sf_feedback'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(DFIRE_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_dfire_sf_feedback

#endif
END MODULE mo_dfire_interface
