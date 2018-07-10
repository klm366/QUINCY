!>
!! @brief Contains the interface to the radiation processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-19)
!!
MODULE mo_rad_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, &
                                            A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
  USE mo_rad_update_radiation,        ONLY: update_radiation 

  ! Use of process configurations (t_PROC_config)
  dsl4jsb_Use_config(RAD_)

  ! Use of process memories (t_PROC_memory)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(L2A_)
  dsl4jsb_Use_memory(VEG_)


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_rad_tasks, update_radiation

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_radiation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_radiation      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_radiation !< Aggregates computed task variables
  END TYPE tsk_update_radiation

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_radiation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_radiation
    PROCEDURE Create_task_update_radiation         !< Constructor function for task
  END INTERFACE tsk_update_radiation


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_rad_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: RAD_ process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_rad_tasks(this, model_id)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id


    CALL this%Register_task(tsk_update_radiation(model_id))

  END SUBROUTINE Register_rad_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_radiation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_radiation(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_radiation::return_ptr)
    CALL return_ptr%Construct(name='update_radiation', process_id=RAD_, owner_model_id=model_id)
    
  END FUNCTION Create_task_update_radiation


  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_radiation task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_radiation(tile, options)

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
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_radiation'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(RAD_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(L2A_)
    dsl4jsb_Def_memory(VEG_)
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
    dsl4jsb_Get_config(RAD_)
    ! Get process memories
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(L2A_)
    dsl4jsb_Get_memory(VEG_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      dsl4jsb_Aggregate_onChunk(RAD_, sw_srf_net              , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, swvis_srf_net           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, swnir_srf_net           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, rad_srf_net             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis                 , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir                 , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_soil            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_soil            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_vis_can             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, alb_nir_can             , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_vis_soil           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_nir_soil           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_vis_can            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_nir_can            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, appfd                   , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, rfr_ratio_boc           , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, ppfd_sunlit_cl          , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, ppfd_shaded_cl          , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_sunlit_vis_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_shaded_vis_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, fleaf_sunlit_vis_cl     , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_sunlit_nir_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, arad_shaded_nir_cl      , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(RAD_, fleaf_sunlit_nir_cl     , weighted_by_fract)
      ! ---------------------------
      dsl4jsb_Aggregate_onChunk(L2A_, albedo                  , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(L2A_, net_radiation           , weighted_by_fract)
      ! ---------------------------
      dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_cl         , weighted_by_fract)
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_radiation

#endif
END MODULE mo_rad_interface
