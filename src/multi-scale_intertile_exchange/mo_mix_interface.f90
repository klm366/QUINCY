!>
!! @brief Contains the interfaces to the multi-scale intertile exchange process
!!
!! @author
!! Tobias Stacke, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Tobias Stacke (2022-07-22)
!!
!! @par Copyright
!! This code is subject to the MPI-M-Software - License - Agreement in it's most recent form.
!! Please see URL http://www.mpimet.mpg.de/en/science/models/license .
!!

MODULE mo_mix_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes MIX_

  ! Use of process configurations
  dsl4jsb_Use_config(MIX_)

  ! Use of process memories
  dsl4jsb_Use_memory(MIX_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_mix_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_interface'
  
  !> Type definition for hydro_mix task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_hydro_mix
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_hydro_mix     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_hydro_mix  !< Aggregates computed task variables
  END TYPE tsk_hydro_mix
  
  !> Constructor interface for hydro_mix task
  INTERFACE tsk_hydro_mix
    PROCEDURE Create_task_hydro_mix                        !< Constructor function for task
  END INTERFACE tsk_hydro_mix

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for hydro_mix task
  !!
  FUNCTION Create_task_hydro_mix(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_hydro_mix::return_ptr)
    CALL return_ptr%Construct(name='hydro_mix', process_id=MIX_, owner_model_id=model_id)

  END FUNCTION Create_task_hydro_mix
  
  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for multi-scale intertile exchange process
  !!
  SUBROUTINE Register_mix_tasks(process, model_id)

    CLASS(t_jsb_process),  INTENT(inout) :: process
    INTEGER,               INTENT(in)    :: model_id

    CALL process%Register_task(tsk_hydro_mix(model_id))
    
  END SUBROUTINE Register_mix_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "hydro_mix"
  !!
  SUBROUTINE update_hydro_mix(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_hydro_mix'

    ! Declare process memories
    dsl4jsb_Def_config(MIX_)
    dsl4jsb_Def_memory(MIX_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_storage

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(MIX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    
    dsl4jsb_Get_config(MIX_)

    ! Get process memories
    dsl4jsb_Get_memory(MIX_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(MIX_, wtr_storage)     ! inout

    ! Just some dummy stuff to see that the general process implementation is working
    IF (dsl4jsb_Config(MIX_)%l_hydro_mix) THEN
      CALL message(TRIM(routine), 'Transport magic will happen here')
    END IF
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_hydro_mix
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "hydro_mix"
  !!
  SUBROUTINE aggregate_hydro_mix(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(MIX_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_hydro_mix'
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(MIX_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(MIX_, wtr_storage,  weighted_by_fract)
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_hydro_mix
  
#endif
END MODULE mo_mix_interface
