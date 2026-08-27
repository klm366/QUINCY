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
MODULE mo_q_dfire_fuel_sf_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, Q_DFIRE_FUEL_SF_, Q_DFIRE_SF_
  USE mo_q_dfire_sf_fuel,             ONLY: update_q_dfire_sf_fuel
  USE mo_q_dfire_sf_impacts,          ONLY: update_q_dfire_sf_impacts
  USE mo_q_dfire_sf_update,           ONLY: update_q_dfire_sf_update

  dsl4jsb_Use_memory(Q_DFIRE_FUEL_SF_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_q_dfire_fuel_sf_tasks
  PUBLIC :: update_q_dfire_sf_fuel, update_q_dfire_sf_impacts, update_q_dfire_sf_update


  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_fuel_sf task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_dfire_fuel_sf
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_dfire_sf_fuel      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_q_dfire_fuel_sf        !< Aggregates computed task variables
  END TYPE tsk_q_dfire_fuel_sf

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_fuel task
  !>
  INTERFACE tsk_q_dfire_fuel_sf
    PROCEDURE Create_task_q_dfire_fuel_sf         !< Constructor function for task
  END INTERFACE tsk_q_dfire_fuel_sf

  ! ======================================================================================================= !
  !>
  !> Type definition: tsk_dfire_sf_impacts task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_dfire_sf_impacts
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_dfire_sf_impacts      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_q_dfire_sf_impacts        !< Aggregates computed task variables
  END TYPE tsk_q_dfire_sf_impacts

  ! ======================================================================================================= !
  !>
  !> Constructor interface: tsk_dfire_sf_impacts task
  !>
  INTERFACE tsk_q_dfire_sf_impacts
    PROCEDURE Create_task_q_dfire_sf_impacts         !< Constructor function for task
  END INTERFACE tsk_q_dfire_sf_impacts

  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_fuel task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_dfire_sf_update
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_dfire_sf_update      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_q_dfire_sf_update        !< Aggregates computed task variables
  END TYPE tsk_q_dfire_sf_update

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_update task
  !>
  INTERFACE tsk_q_dfire_sf_update
    PROCEDURE Create_task_q_dfire_sf_update         !< Constructor function for task
  END INTERFACE tsk_q_dfire_sf_update

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_fuel_sf_interface'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Register tasks: Q_DFIRE_SF_ process
  !>
  SUBROUTINE Register_q_dfire_fuel_sf_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_q_dfire_fuel_sf(model_id))
    CALL this%Register_task(tsk_q_dfire_sf_impacts(model_id))
    CALL this%Register_task(tsk_q_dfire_sf_update(model_id))
  END SUBROUTINE Register_q_dfire_fuel_sf_tasks

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_fuel task
  !>
  FUNCTION Create_task_q_dfire_fuel_sf(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_dfire_fuel_sf::return_ptr)
    CALL return_ptr%Construct(name='q_dfire_sf_fuel', process_id=Q_DFIRE_FUEL_SF_, owner_model_id=model_id)
  END FUNCTION Create_task_q_dfire_fuel_sf

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_fuel task
  !>
  SUBROUTINE aggr_q_dfire_fuel_sf(tile, options)
    dsl4jsb_Use_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    CLASS(t_jsb_aggregator),  POINTER         :: sum_agg
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_q_dfire_fuel_sf'
    LOGICAL                                   :: is_mask_enabled
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    sum_agg => tile%Get_aggregator("sum")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_SF_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, moisture_extinction  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, fwd_wind_factor      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, fuel_bulk_density    , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_1hr_dfp    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_10hr_dfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_100hr_dfp  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, dead_fuel_1000hr_dfp , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, live_grass_fuel_1hr_lfp    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, cured_live_grass_fuel_1hr_dfp    , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_1hr_lfp     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_10hr_lfp    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_100hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, live_tree_fuel_1000hr_lfp  , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, lvg_moisture  , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, forest_fract     , sum_agg)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, grass_fract      , sum_agg)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_q_dfire_fuel_sf

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_fuel task
  !>
  FUNCTION Create_task_q_dfire_sf_impacts(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_dfire_sf_impacts::return_ptr)
    CALL return_ptr%Construct(name='q_dfire_sf_impacts', process_id=Q_DFIRE_FUEL_SF_, owner_model_id=model_id)
  END FUNCTION Create_task_q_dfire_sf_impacts

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_fuel task
  !>
  SUBROUTINE aggr_q_dfire_sf_impacts(tile, options)
    dsl4jsb_Use_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_q_dfire_sf_impacts'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_SF_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, veg_fire_litterfall_total_els , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, veg_fire_emissions_total_els , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, sb_fire_emissions_total_els , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, sb_fire_litter_ash_total_els , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, pm_fire      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, pm_fire_output , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, pm_prop_ck   , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, pm_tau, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, tau_l, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, tau_c, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, bark_thickness, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, DBH_diag, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, par1_diag, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, par2_diag, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, ck, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, r_ck, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, ck_p, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, scorch_height, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, F_scorch_height, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, crown_length, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, tree_height, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, fuel_consumption_els, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, emissions_spc, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, emissions_els, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, residue_els, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, overflow_els, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_q_dfire_sf_impacts

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_update task
  !>
  FUNCTION Create_task_q_dfire_sf_update(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_dfire_sf_update::return_ptr)
    CALL return_ptr%Construct(name='q_dfire_sf_update', process_id=Q_DFIRE_FUEL_SF_, owner_model_id=model_id)
  END FUNCTION Create_task_q_dfire_sf_update

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_sf_update task
  !>
  SUBROUTINE aggr_q_dfire_sf_update(tile, options)
    dsl4jsb_Use_config(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_q_dfire_sf_update'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_DFIRE_SF_)
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_DFIRE_SF_)
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_SF_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

!    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_SF_, veg_fire_litterfall_total_els , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_q_dfire_sf_update


#endif
END MODULE mo_q_dfire_fuel_sf_interface
