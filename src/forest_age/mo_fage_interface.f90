!>
!! @brief Contains the interfaces to the forest-age processes
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2017-11-21)
!!
MODULE mo_fage_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish, message, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  dsl4jsb_Use_processes FAGE_

  ! Use of process configurations (t_<process>_config)
  dsl4jsb_Use_config(FAGE_)

  ! Use of process memories (t_<process>_memory)
  dsl4jsb_Use_memory(FAGE_)  
  
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_fage_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_interface'

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_ageing
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_fage    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_fage !< Aggregates computed task variables
  END TYPE   tsk_ageing

  INTERFACE   tsk_ageing
    PROCEDURE Create_task_ageing 
  END INTERFACE   tsk_ageing

CONTAINS
  ! ================================================================================================================================
  !! Constructors for tasks
  
  !> Constructor for ageing task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "ageing"
  !!
  FUNCTION Create_task_ageing(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_ageing::return_ptr)

    return_ptr%name = 'ageing'
    return_ptr%process_id = FAGE_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_ageing

  ! =======================================================================================================
  !> Register tasks for fage process
  !!
  !! @param[in,out] this      Instance of fage process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_fage_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_ageing(model_id))

  END SUBROUTINE Register_fage_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the fage state
  !! Task "ageing" calculates the age  - TODO-JN and more comming
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_fage(tile, options)

    USE mo_fage_process, ONLY: ageing
    USE mo_jsb_time,   ONLY: is_newyear

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age

    INTEGER  :: iblk, ics, ice, nc
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_fage'

    ! ---------------------------
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (.NOT. tile%Is_process_active(FAGE_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    dsl4jsb_Get_var2D_onChunk(FAGE_, mean_age)
    dsl4jsb_Get_var3D_onChunk(FAGE_, fract_per_age)

    CALL ageing(nc, dsl4jsb_Config(FAGE_)%max_age, is_newyear(), fract_per_age(:,:), mean_age(:))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_fage

  ! TODO-JN-COHORTS: does not really make sense to aggregate mean_age?! 
  SUBROUTINE aggregate_fage(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(FAGE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_fpc'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(FAGE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_fage

#endif
END MODULE mo_fage_interface
