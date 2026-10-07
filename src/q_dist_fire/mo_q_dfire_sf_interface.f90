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
MODULE mo_q_dfire_sf_interface
#ifndef __NO_JSBACH__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, Q_DFIRE_SF_
  USE mo_q_dfire_sf_dynamics,         ONLY: update_q_dfire_sf_dynamics

  dsl4jsb_Use_memory(Q_DFIRE_SF_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_q_dfire_sf_tasks
  PUBLIC :: update_time_average_q_dist_fire, update_q_dfire_sf_dynamics

  ! ======================================================================================================= !
  !>
  !> Type definition: dfire_sf_fuel task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_dfire_sf_dynamics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_dfire_sf_dynamics      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_q_dfire_sf_dynamics        !< Aggregates computed task variables
  END TYPE tsk_q_dfire_sf_dynamics

  ! ======================================================================================================= !
  !>
  !> Constructor interface: dfire_sf_fuel task
  !>
  INTERFACE tsk_q_dfire_sf_dynamics
    PROCEDURE Create_task_q_dfire_sf_dynamics         !< Constructor function for task
  END INTERFACE tsk_q_dfire_sf_dynamics

  ! ======================================================================================================= !
  !>
  !> Constructor interface: time_average_dfire task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_q_dist_fire
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_q_dist_fire    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_q_dist_fire !< Aggregates computed task variables
  END TYPE tsk_time_average_q_dist_fire

  ! ======================================================================================================= !
  !>
  !> Constructor interface: update_time__average_dfire task
  !>
  INTERFACE tsk_time_average_q_dist_fire
    PROCEDURE Create_task_update_time_average_q_dist_fire         !< Constructor function for task
  END INTERFACE tsk_time_average_q_dist_fire

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dfire_sf_interface'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Register tasks: Q_DFIRE_SF_ process
  !>
  SUBROUTINE Register_q_dfire_sf_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_q_dfire_sf_dynamics(model_id))
    CALL this%Register_task(tsk_time_average_q_dist_fire(model_id))
  END SUBROUTINE Register_q_dfire_sf_tasks

  ! ======================================================================================================= !
  !>
  !> Constructor: dfire_sf_dynamics task
  !>
  FUNCTION Create_task_q_dfire_sf_dynamics(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_dfire_sf_dynamics::return_ptr)
    CALL return_ptr%Construct(name='q_dfire_sf_dynamics', process_id=Q_DFIRE_SF_, owner_model_id=model_id)
  END FUNCTION Create_task_q_dfire_sf_dynamics

  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": dfire_sf_dynamics task
  !>
  SUBROUTINE aggr_q_dfire_sf_dynamics(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER         :: model
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice, nc
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_q_dfire_sf_dynamics'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_SF_)
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
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Get_memory(Q_DFIRE_SF_)

    ! Q_DFIRE_SF_ 2D
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, fuel_moisture , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, ab_moisture_diag, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, t_air_max, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, t_air_min, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_1hr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_10hr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_100hr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, dead_fuel_consumption_1000hr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, area_burned_fract, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, area_burned, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, area_burned_fract_output, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, area_burned_output, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, i_surface, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, precip_qc_mavg, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, precip_daytime, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, precip_mavg, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, precip_daily_mmd, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, precip_daily, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, reaction_intensity, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, prop_flux_ratio, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, wind_factor, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, eff_heating_number, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, preig_heat, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, opt_react_velo, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, net_fuel_load, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, moist_damp_coeff, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, mineral_damp_coeff, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, nesterov_index, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, fdi, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, ros_f, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, ros_b, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, mean_fire_area, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, fire_duration, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, distance_travelled, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, l2b_ratio, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, e_n_ignitions, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, e_n_humans, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, e_n_lightning, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, frac_consump_fcl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_DFIRE_SF_, pot_burn_area, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggr_q_dfire_sf_dynamics


  ! ======================================================================================================= !
  !>
  !> calculate time moving averages and daytime averages for Q_DFIRE_SF_
  !>
  FUNCTION Create_task_update_time_average_q_dist_fire(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_q_dist_fire::return_ptr)
    CALL return_ptr%Construct(name='tavrg_q_dist_fire', process_id=Q_DFIRE_SF_, owner_model_id=model_id)

  END FUNCTION Create_task_update_time_average_q_dist_fire


  ! ======================================================================================================= !
  !>
  !> calculate time moving averages and daytime averages for Q_DFIRE_SF_
  !>
  SUBROUTINE update_time_average_q_dist_fire(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_q_dfire_sf_util,     ONLY: calculate_time_average_q_dist_fire,  &
                                      calculate_daily_weather, &
                                      calculate_pure_moving_averages

    dsl4jsb_Use_processes Q_DFIRE_SF_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_q_dist_fire'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_DFIRE_SF_)) RETURN
    !IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
            & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL calculate_time_average_q_dist_fire(tile, options)
    CALL calculate_daily_weather(tile, options)
    CALL calculate_pure_moving_averages(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_q_dist_fire


  ! ======================================================================================================= !
  !>
  !> Implementation of "aggregate": time_average_dfire task
  !>
  SUBROUTINE aggregate_time_average_q_dist_fire(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_DFIRE_SF_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_q_dist_fire'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_DFIRE_SF_)

    ! Q_DFIRE_SF_ 2D

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_q_dist_fire

#endif
END MODULE mo_q_dfire_sf_interface
