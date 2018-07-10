!>
!! @brief Contains the interface to the phenology processes
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20161214)
!! Second version by Steffen Richter (June 2018)
!!
MODULE mo_pheno_interface
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
  USE mo_pheno_update_phenology,      ONLY: update_phenology

  ! Use of process configurations (t_PROC_config)
  dsl4jsb_Use_config(PHENO_)

  ! Use of process memories (t_PROC_memory)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(VEG_)
  

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pheno_tasks, update_phenology

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_phenology task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_phenology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_phenology      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggr_update_phenology !< Aggregates computed task variables
  END TYPE tsk_update_phenology

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_phenology task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_update_phenology
    PROCEDURE Create_task_update_phenology         !< Constructor function for task
  END INTERFACE tsk_update_phenology


  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_interface'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: PHENO_ process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_pheno_tasks(this, model_id)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id


    CALL this%Register_task(tsk_update_phenology(model_id))

  END SUBROUTINE Register_pheno_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_phenology task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_phenology(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_update_phenology::return_ptr)
    CALL return_ptr%Construct(name='update_phenology', process_id=PHENO_, owner_model_id=model_id)

  END FUNCTION Create_task_update_phenology


  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": update_phenology task
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggr_update_phenology(tile, options)

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
    CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_phenology'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)
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
    dsl4jsb_Get_config(PHENO_)
    ! Get process memories
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(VEG_)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go aggregate
    ! ------------------------------------------------------------------------------------------------------------
    

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      dsl4jsb_Aggregate_onChunk(PHENO_, l_growing_season            , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, gdd                         , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, nd_dormance                 , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(PHENO_, lai_max                     , weighted_by_fract)    ! *slmdev note: currently this var is not changed after init
      ! ---------------------------
      dsl4jsb_Aggregate_onChunk(VEG_, growth_req_n_tlabile_mavg     , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, growth_req_p_tlabile_mavg     , weighted_by_fract)
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen pheno Tasks aggregiert werden sollen.

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggr_update_phenology

#endif
END MODULE mo_pheno_interface
