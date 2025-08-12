!>
!>  @par Copyright
!>    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!>    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!>    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!>
!>    This file is part of the QUINCY model.
!>    The QUINCY model is free software: you can redistribute it and/or modify
!>    it under the terms of the GNU General Public License version 3,
!>    as published by the Free Software Foundation.
!>
!>    The QUINCY model is distributed in the hope that it will be useful,
!>    but WITHOUT ANY WARRANTY; without even the implied warranty of
!>    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!>    GNU General Public License for more details.
!>
!>    You should have received a copy of the GNU General Public License
!>    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!> 
!> @brief Contains the interface to the dfire processes
!>
!> @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!>
MODULE mo_q_dfire_fuel_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, Q_DFIRE_FUEL_
  USE mo_q_dfire_fuel,                ONLY: update_q_dfire_fuel
  USE mo_q_dfire_post,                ONLY: update_q_dfire_post

  USE mo_q_dfire_constants,           ONLY: q_dfire_pft_mask ! pfts to compute average

  dsl4jsb_Use_memory(Q_DFIRE_FUEL_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_q_dfire_fuel_tasks
  PUBLIC :: update_q_dfire_fuel, update_q_dfire_post


  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_fuel task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_dfire_fuel
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_dfire_fuel      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_q_dfire_fuel        !< Aggregates computed task variables
  END TYPE tsk_q_dfire_fuel

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_fuel task
  !>
  INTERFACE tsk_q_dfire_fuel
    PROCEDURE Create_task_q_dfire_fuel         !< Constructor function for task
  END INTERFACE tsk_q_dfire_fuel

  ! ======================================================================================================= !
  !>
  !> Type definition: tsk_dfire_post task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_dfire_post
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_dfire_post      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_q_dfire_post        !< Aggregates computed task variables
  END TYPE tsk_q_dfire_post

  ! ======================================================================================================= !
  !>
  !> Constructor interface: tsk_dfire_post task
  !>
  INTERFACE tsk_q_dfire_post
    PROCEDURE Create_task_q_dfire_post         !< Constructor function for task
  END INTERFACE tsk_q_dfire_post

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_fuel_interface'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Register tasks: Q_DFIRE_ process
  !>
  SUBROUTINE Register_q_dfire_fuel_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_q_dfire_fuel(model_id))
    CALL this%Register_task(tsk_q_dfire_post(model_id))
  END SUBROUTINE Register_q_dfire_fuel_tasks

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_fuel task
  !>
  FUNCTION Create_task_q_dfire_fuel(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_dfire_fuel::return_ptr)
    CALL return_ptr%Construct(name='q_dfire_fuel', process_id=Q_DFIRE_FUEL_, owner_model_id=model_id)
  END FUNCTION Create_task_q_dfire_fuel

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_fuel task
  !>
  SUBROUTINE aggr_q_dfire_fuel(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    CLASS(t_jsb_aggregator),  POINTER         :: sum_agg, average_agg
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_q_dfire_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_)
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
    average_agg => tile%Get_aggregator("average")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !

    ! Lctlib related
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, moisture_extinction      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, fwd_wind_factor          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk_masked(Q_DFIRE_FUEL_, fuel_bulk_density , average_agg, q_dfire_pft_mask)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_10hr_dfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_100hr_dfp  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_1000hr_dfp , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_1hr_lfp    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_10hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_100hr_lfp  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_1000hr_lfp , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_grass_fuel_1hr_lfp    , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_10hr_lfp    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_100hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1000hr_lfp  , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, forest_fract     , sum_agg)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, grass_fract      , sum_agg)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_q_dfire_fuel

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_fuel task
  !>
  FUNCTION Create_task_q_dfire_post(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_dfire_post::return_ptr)
    CALL return_ptr%Construct(name='q_dfire_post', process_id=Q_DFIRE_FUEL_, owner_model_id=model_id)
  END FUNCTION Create_task_q_dfire_post

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_fuel task
  !>
  SUBROUTINE aggr_q_dfire_post(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_q_dfire_fuel'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_FUEL_)
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
    dsl4jsb_Get_memory(Q_DFIRE_FUEL_)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Go aggregate
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_10hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_100hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_fuel_1000hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, forest_fract   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, grass_fract   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, fuel_bulk_density   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, moisture_extinction   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, fwd_wind_factor   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, dead_fuel_1hr_dfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_grass_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_grass_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_grass_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_grass_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, live_tree_fuel_1hr_lfp   , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, pm_fire      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, pm_prop_ck   , weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, pm_tau, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, tau_l, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, tau_c, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, bark_thickness, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, DBH_diag, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, par1_diag, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, par2_diag, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, ck, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, r_ck, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, ck_p, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, scorch_height, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, F_scorch_height, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, crown_length, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, tree_height, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, fuel_consumption_els, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, emissions_spc, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, emissions_els, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, residue_els, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_FUEL_, overflow_els, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_q_dfire_post


#endif
END MODULE mo_q_dfire_fuel_interface
