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

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_pheno_update_phenology,  ONLY: update_phenology
  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_class,               ONLY: Get_model
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,       ONLY: t_jsb_process
  USE mo_jsb_task_class,          ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  USE mo_jsb_process_class, ONLY: SEB_, A2L_, RAD_, HYDRO_, ASSIMI_, PHENO_

  ! Use of process configurations (t_<process>_config)
  dsl4jsb_Use_config(PHENO_)

  ! Use of process memories (t_<process>_memory)
  dsl4jsb_Use_memory(PHENO_)
  
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_pheno_tasks, update_phenology

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_interface'

  !> Type definition for update_phenology task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_phenology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_phenology0    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_phenology !< Aggregates computed task variables
  END TYPE   tsk_phenology

  !> Constructor interface for update_phenology task
  INTERFACE   tsk_phenology
    PROCEDURE Create_task_phenology         !< Constructor function for task
  END INTERFACE   tsk_phenology

CONTAINS
  ! ================================================================================================================================
  !! Constructors for tasks
  
  !> Constructor for update_phenology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "update_phenology"
  !!
  FUNCTION Create_task_phenology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_phenology::return_ptr)

    return_ptr%name = 'pheno'
    return_ptr%process_id = PHENO_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_phenology
  
  ! =======================================================================================================
  !> Register tasks for pheno process
  !!
  !! @param[in,out] this      Instance of pheno process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_pheno_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_phenology(model_id))

  END SUBROUTINE Register_pheno_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the phenology state
  !! Task "phenology" calculates the phenology state of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_phenology0(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_phenology(tile, options)

  END SUBROUTINE update_phenology0

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "phenology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_phenology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER                :: model

    dsl4jsb_Def_config(PHENO_)
    dsl4jsb_Def_memory(PHENO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_phenology'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(PHENO_)

    dsl4jsb_Get_memory(PHENO_)

    IF (tile%Has_children()) THEN
      weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

      dsl4jsb_Aggregate_onChunk(PHENO_, testvar, weighted_by_fract)
    END IF

    !X Implementation: Start your process scheme here...
    !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen pheno Tasks aggregiert werden sollen.

    ! IF (TRIM(dsl4jsb_Config(PHENO_)%scheme) == 'logrop') THEN
    ! 
    ! END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_phenology

#endif
END MODULE mo_pheno_interface
