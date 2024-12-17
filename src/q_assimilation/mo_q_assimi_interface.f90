!> QUINCY assimilation process interface
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of tasks for the assimilation process (photosynthesis) incl. update and aggregate routines
!>
MODULE mo_q_assimi_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message
  USE mo_jsb_grid_class,      ONLY: t_jsb_grid
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes Q_ASSIMI_, VEG_
  dsl4jsb_Use_config(Q_ASSIMI_)
  dsl4jsb_Use_memory(Q_ASSIMI_)
  dsl4jsb_Use_memory(VEG_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_q_assimi_tasks
  PUBLIC ::  update_canopy_fluxes, update_reset_q_assimi_fluxes, update_time_average_q_assimilation, update_q_npp_buffer

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_assimi_interface'

  ! ======================================================================================================= !
  !> Type definition: reset_q_assimi_fluxes task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_q_assimi_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_q_assimi_fluxes     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_q_assimi_fluxes  !< Aggregates computed task variables
  END TYPE tsk_reset_q_assimi_fluxes
  !> Constructor interface: reset_q_assimi_fluxes task
  !>
  INTERFACE tsk_reset_q_assimi_fluxes
    PROCEDURE Create_task_reset_q_assimi_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_q_assimi_fluxes

  ! ======================================================================================================= !
  !> Type definition: canopy_fluxes task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_fluxes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_fluxes !< Aggregates computed task variables
  END TYPE tsk_canopy_fluxes
  !> Constructor interface: update_canopy_fluxes task
  !>
  INTERFACE tsk_canopy_fluxes
    PROCEDURE Create_task_canopy_fluxes         !< Constructor function for task
  END INTERFACE tsk_canopy_fluxes

  !> Type definition for NPP buffer task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_npp_buffer
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_npp_buffer    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_q_npp_buffer !< Aggregates computed task variables
  END TYPE tsk_q_npp_buffer

  !> Constructor interface for assimilation task
  INTERFACE tsk_q_npp_buffer
    PROCEDURE Create_task_q_npp_buffer     !< Constructor function for task
  END INTERFACE tsk_q_npp_buffer

  ! ======================================================================================================= !
  !> Type definition: time_average_assimilation task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_assimilation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_q_assimilation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_q_assimilation !< Aggregates computed task variables
  END TYPE tsk_time_average_assimilation
  !> Constructor interface: time_average_assimilation task
  !>
  INTERFACE tsk_time_average_assimilation
    PROCEDURE Create_task_time_average_assimilation         !< Constructor function for task
  END INTERFACE tsk_time_average_assimilation

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_q_assimi_fluxes task
  !>
  FUNCTION Create_task_reset_q_assimi_fluxes(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_q_assimi_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_q_assimi_fluxes', process_id=Q_ASSIMI_, owner_model_id=model_id)
  END FUNCTION Create_task_reset_q_assimi_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_canopy_fluxes task
  !>
  FUNCTION Create_task_canopy_fluxes(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_fluxes::return_ptr)
    CALL return_ptr%Construct(name='canopy_fluxes', process_id=Q_ASSIMI_, owner_model_id=model_id)
  END FUNCTION Create_task_canopy_fluxes

  ! ======================================================================================================= !
  !> Constructor: tsk_time_average_assimilation task
  !>
  FUNCTION Create_task_time_average_assimilation(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_assimilation::return_ptr)
    CALL return_ptr%Construct(name='tavrg_assimilation', process_id=Q_ASSIMI_, owner_model_id=model_id)
  END FUNCTION Create_task_time_average_assimilation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for npp_buffer task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "npp_buffer"
  !!
  FUNCTION Create_task_q_npp_buffer(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_npp_buffer::return_ptr)
    CALL return_ptr%Construct(name='q_npp_buffer', process_id=Q_ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_q_npp_buffer


  ! ======================================================================================================= !
  !> Register tasks for assimi process
  !>
  SUBROUTINE Register_q_assimi_tasks(this, model_id)
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    TYPE(t_jsb_model), POINTER :: model

    model => Get_model(model_id)

    CALL this%Register_task(tsk_reset_q_assimi_fluxes(model_id))
    CALL this%Register_task(tsk_canopy_fluxes(model_id))
    CALL this%Register_task(tsk_time_average_assimilation(model_id))
    CALL this%Register_task(tsk_q_npp_buffer(model_id))
  END SUBROUTINE Register_q_assimi_tasks

  ! ======================================================================================================= !
  !> Implementation of "update": reset_q_assimi_fluxes task
  !>
  SUBROUTINE update_reset_q_assimi_fluxes(tile, options)
    USE mo_q_assimi_util, ONLY: reset_q_assimi_fluxes
    dsl4jsb_Use_processes Q_ASSIMI_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_q_assimi_fluxes'

    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN

    iblk = options%iblk

    IF (debug_on()) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL reset_q_assimi_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_reset_q_assimi_fluxes

  ! ======================================================================================================= !
  !> Implementation of "aggregate": reset_q_assimi_fluxes task
  !>
  SUBROUTINE aggregate_reset_q_assimi_fluxes(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_q_assimi_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)

    ! Q_ASSIMI_ 2D
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C13                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C14                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, net_assimilation                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, net_assimilation_boc                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, maint_respiration_leaf                , weighted_by_fract)
    ! Q_ASSIMI_ 3D
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_cl                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, net_assimilation_cl                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, maint_respiration_leaf_cl             , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_reset_q_assimi_fluxes

  ! ======================================================================================================= !
  !> Implementation of "update": canopy_fluxes task
  !>
  SUBROUTINE update_canopy_fluxes(tile, options)
    USE mo_q_assimi_process, ONLY: real_update_canopy_fluxes => update_canopy_fluxes

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggr_canopy_fluxes'

    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN

    iblk = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_update_canopy_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_canopy_fluxes

  ! ======================================================================================================= !
  !> Implementation of "aggregate": canopy_fluxes task
  !>
  SUBROUTINE aggregate_canopy_fluxes(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)
    ! dsl4jsb_Get_memory(A2L_)

    ! Q_ASSIMI_ 2D
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C13      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_C14      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, net_assimilation            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, net_assimilation_boc        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, maint_respiration_leaf      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, aerodyn_cond                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, canopy_cond                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, co2_conc_leaf               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_air                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soa                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_gs                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_ps                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, t_jmax_opt                  , weighted_by_fract)
    ! Q_ASSIMI_ 3D
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, ftranspiration_sl           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, net_assimilation_cl         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, gross_assimilation_cl       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, maint_respiration_leaf_cl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, canopy_cond_cl              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, co2_conc_leaf_cl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, jmax_cl                     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, vcmax_cl                    , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE aggregate_canopy_fluxes

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for Q_ASSIMI_
  !>
  SUBROUTINE update_time_average_q_assimilation(tile, options)
    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_q_assimi_util,         ONLY: calculate_time_average_q_assimilation
    dsl4jsb_Use_processes Q_ASSIMI_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_q_assimilation'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) THEN
      CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    END IF

    CALL calculate_time_average_q_assimilation(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_time_average_q_assimilation

  ! ======================================================================================================= !
  !> Implementation of "aggregate": time_average_assimilation task
  !>
  SUBROUTINE aggregate_time_average_q_assimilation(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_q_assimilation'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_ASSIMI_)

    ! Q_ASSIMI_
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, soa_tsoa_mavg               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soa_tphen_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_gs_tphen_mavg     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_ps_tfrac_mavg     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_gs_tfrac_mavg     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_ps_tcnl_mavg      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_gs_tcnl_mavg      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_air_tfrac_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_air_tcnl_mavg          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_air_daytime            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_ps_daytime        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_gs_daytime        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_air_daytime_dacc       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_ps_daytime_dacc   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, beta_soil_gs_daytime_dacc   , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_time_average_q_assimilation

  ! ================================================================================================================================
  !>
  !> Implementation of NPP buffer
  !! Task "npp_buffer" calculates multi-annual average of potential  NPP for use in the NLCC process
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_q_npp_buffer(tile, options)

    USE mo_q_assimi_process,    ONLY: calc_q_NPPbuf
    USE mo_jsb_time,            ONLY: is_newyear, is_time_experiment_start

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER                    :: iblk, ics, ice !, nc
    REAL(wp)                   :: dtime
    LOGICAL                    :: lstart, init_running_means, init_npp, new_year

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_q_npp_buffer'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
!     dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)

    ! Declare pointers to variables in memory

    dsl4jsb_Real2D_onChunk :: npp
!     dsl4jsb_Real2D_onChunk :: net_assimilation
    dsl4jsb_Real2D_onChunk :: seconds_year
!     dsl4jsb_Real2D_onChunk :: NPP_pot_rate_ca
    dsl4jsb_Real2D_onChunk :: NPP_sum_year
    dsl4jsb_Real2D_onChunk :: NPP_mean_5year
    dsl4jsb_Real2D_onChunk :: land_cover_class
    dsl4jsb_Real2D_onChunk :: bclimit_min_cold_mmtemp
    dsl4jsb_Real2D_onChunk :: bclimit_max_cold_mmtemp
    dsl4jsb_Real2D_onChunk :: bclimit_max_warm_mmtemp
    dsl4jsb_Real2D_onChunk :: bclimit_min_temprange
    dsl4jsb_Real2D_onChunk :: bclimit_min_gdd
    dsl4jsb_Real2D_onChunk :: tau_c_woods
    dsl4jsb_Real2D_onChunk :: sla

    ! Set local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    !nc      = options%nc
    dtime   = options%dtime

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
!     IF (iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set process config and memories
    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(Q_ASSIMI_)
!     dsl4jsb_Get_config(VEG_)

    ! Get process memories
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(VEG_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)                                  ! out

!     dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   net_assimilation)         ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   seconds_year)            ! inout
!     dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   NPP_pot_rate_ca)            ! inout
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   NPP_sum_year)            ! inout
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   NPP_mean_5year)          ! inout
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   land_cover_class)        ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   bclimit_min_cold_mmtemp) ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   bclimit_max_cold_mmtemp) ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   bclimit_max_warm_mmtemp) ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   bclimit_min_temprange)   ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   bclimit_min_gdd)         ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   tau_c_woods)             ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_,   sla)                     ! out

    ! Accumulate NPP_sum_year and update NPP_mean_5year for NLCC
    lstart = is_time_experiment_start(options%current_datetime)
    new_year = is_newyear(options%current_datetime, dtime)

    init_running_means = dsl4jsb_Config(Q_ASSIMI_)%init_running_means
    init_NPP = .TRUE.
!     NPP_pot_rate_ca(:) = net_assimilation(:)
!     NPP_pot_rate_ca(:) = npp(:)
    IF (MINVAL(NPP_mean_5year(:)) < -999._wp) init_NPP = .FALSE.
!     CALL calc_q_NPPbuf(lstart, init_running_means, new_year, dtime, net_assimilation, &
    CALL calc_q_NPPbuf(lstart, init_running_means, new_year, dtime, npp, &
                     seconds_year, NPP_sum_year, NPP_mean_5year)

    IF (.NOT. init_NPP .AND. debug_on() .AND. iblk==1) THEN
       IF (MINVAL(NPP_mean_5year(:)) >= -999._wp) CALL message(TRIM(routine), 'initialisation of NPP_mean_5year')
    ENDIF

    ! write lct information on variables to make them available on pft-tiles via function collect_var for NLCC
    land_cover_class(:) = REAL(dsl4jsb_Lctlib_param(LandcoverClass), wp)
    bclimit_min_cold_mmtemp(:) = dsl4jsb_Lctlib_param(bclimit_min_cold_mmtemp)
    bclimit_max_cold_mmtemp(:) = dsl4jsb_Lctlib_param(bclimit_max_cold_mmtemp)
    bclimit_max_warm_mmtemp(:) = dsl4jsb_Lctlib_param(bclimit_max_warm_mmtemp)
    bclimit_min_temprange(:) = dsl4jsb_Lctlib_param(bclimit_min_temprange)
    bclimit_min_gdd(:) = dsl4jsb_Lctlib_param(bclimit_min_gdd)
    tau_c_woods(:) = dsl4jsb_Lctlib_param(tau_c_woods)
    sla(:) = dsl4jsb_Lctlib_param(sla)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_q_npp_buffer

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "npp_buffer"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_q_npp_buffer(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(Q_ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice !, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_npp_buffer'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(Q_ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, seconds_year,             weighted_by_fract)
!     dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, NPP_pot_rate_ca,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, NPP_sum_year,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, NPP_mean_5year,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, land_cover_class,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, bclimit_min_cold_mmtemp,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, bclimit_max_cold_mmtemp,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, bclimit_max_warm_mmtemp,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, bclimit_min_temprange,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, bclimit_min_gdd,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, tau_c_woods,              weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_ASSIMI_, sla,                      weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_q_npp_buffer

#endif
END MODULE mo_q_assimi_interface
