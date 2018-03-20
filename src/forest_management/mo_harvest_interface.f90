!>
!! @brief Contains the interfaces to the harv process 
!!
!! @author
!!  Julia Nabel, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Julia Nabel (2018-03-19)
!!
MODULE mo_harv_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes HARV_, FM_

  ! Use of process configurations
  dsl4jsb_Use_config(HARV_)
  dsl4jsb_Use_config(FM_)

  ! Use of process memories
  dsl4jsb_Use_memory(HARV_)
  dsl4jsb_Use_memory(FM_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_harv_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_harv_interface'

  !> Type definition for harv task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_forest_harvest
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_forest_harvest    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_forest_harvest !< Aggregates computed task variables
  END TYPE tsk_forest_harvest

  !> Constructor interface for harv task
  INTERFACE tsk_forest_harvest
    PROCEDURE Create_task_forest_harvest     !< Constructor function for task
  END INTERFACE tsk_forest_harvest

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for forest_harvest task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "forest_harvest"
  !!
  FUNCTION Create_task_forest_harvest(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_forest_harvest::return_ptr)

    return_ptr%name = 'forest_harvest'
    return_ptr%process_id = HARV_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_forest_harvest

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for harv process
  !!
  !! @param[in,out] this      Instance of harv process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_harv_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_forest_harvest(model_id))

  END SUBROUTINE Register_harv_tasks


  ! ================================================================================================================================
  !>
  !! @brief Implementation of task forest_harvest
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_forest_harvest(tile, options)

    USE mo_jsb_time,       ONLY: is_newyear
    USE mo_jsb_tile_class, ONLY: t_jsb_tile_abstract

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model

    INTEGER :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_forest_harvest'

    ! Declare process memories
    dsl4jsb_Def_memory(HARV_)
    dsl4jsb_Def_memory(FM_)

    dsl4jsb_Def_config(HARV_)
    dsl4jsb_Def_config(FM_)

    dsl4jsb_Real2D_onChunk :: harvested_area
    dsl4jsb_Real2D_onChunk :: area_to_harvest

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(HARV_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get model
    model   => Get_model(tile%owner_model_id)

    ! Get process configurations
    dsl4jsb_Get_config(HARV_)
    dsl4jsb_Get_config(FM_)

    ! Get process memories
    dsl4jsb_Get_memory(HARV_)
    dsl4jsb_Get_memory(FM_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(HARV_,  harvested_area)
    dsl4jsb_Get_var2D_onChunk(FM_,  area_to_harvest)

    !
    ! ---------------------------
    ! Go

    IF (is_newyear()) THEN ! only once per year

      print *, area_to_harvest

      IF (TRIM(dsl4jsb_Config(HARV_)%harv_scheme) == 'max_available_age') THEN
          print *, 'harvest forest from max age'
      ELSE
          CALL finish(TRIM(routine),'Unknown harv scheme')
      ENDIF

    ENDIF ! is_newyear

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_forest_harvest

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "forest_harvest"
  !!
  !! @param[in,out] tile    Tile for which aggregation of child tiles is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_forest_harvest(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HARV_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    TYPE(t_jsb_model), POINTER                :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_forest_harvest'

    INTEGER  :: iblk, ics, ice, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(HARV_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HARV_, harvested_area, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_forest_harvest

#endif
END MODULE mo_harv_interface
