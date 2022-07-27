!>
!! @brief Contains the interfaces to the lateral transport process
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

MODULE mo_lattrans_interface
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
  dsl4jsb_Use_processes LATTRANS_

  ! Use of process configurations
  dsl4jsb_Use_config(LATTRANS_)

  ! Use of process memories
  dsl4jsb_Use_memory(LATTRANS_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_lattrans_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lattrans_interface'
  
  !> Type definition for collect_water task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_collect_water
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_collect_water     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_collect_water  !< Aggregates computed task variables
  END TYPE tsk_collect_water
  
  !> Constructor interface for collect_water task
  INTERFACE tsk_collect_water
    PROCEDURE Create_task_collect_water                        !< Constructor function for task
  END INTERFACE tsk_collect_water

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for collect_water task
  !!
  FUNCTION Create_task_collect_water(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_collect_water::return_ptr)
    CALL return_ptr%Construct(name='collect_water', process_id=LATTRANS_, owner_model_id=model_id)

  END FUNCTION Create_task_collect_water
  
  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for lateral transport process
  !!
  SUBROUTINE Register_lattrans_tasks(process, model_id)

    CLASS(t_jsb_process),  INTENT(inout) :: process
    INTEGER,               INTENT(in)    :: model_id

    CALL process%Register_task(tsk_collect_water(model_id))
    
  END SUBROUTINE Register_lattrans_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "collect_water"
  !!
  SUBROUTINE update_collect_water(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_collect_water'

    ! Declare process memories
    dsl4jsb_Def_config(LATTRANS_)
    dsl4jsb_Def_memory(LATTRANS_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_transport
    dsl4jsb_Real2D_onChunk :: wtr_transstor

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(LATTRANS_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    
    dsl4jsb_Get_config(LATTRANS_)

    ! Get process memories
    dsl4jsb_Get_memory(LATTRANS_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(LATTRANS_, wtr_transstor)     ! inout
    dsl4jsb_Get_var2D_onChunk(LATTRANS_, wtr_transport)     ! out

    ! Just some dummy stuff to see that the general process implementation is working
    IF (dsl4jsb_Config(LATTRANS_)%l_lattrans) THEN
      wtr_transport(:) = 2.0_wp
      wtr_transstor(:) = 2.0_wp
    ELSE
      wtr_transport(:) = 1.0_wp
      wtr_transstor(:) = 1.0_wp
    END IF
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_collect_water
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "collect_water"
  !!
  SUBROUTINE aggregate_collect_water(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(LATTRANS_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_collect_water'
    
    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(LATTRANS_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(LATTRANS_, wtr_transport, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(LATTRANS_, wtr_transstor, weighted_by_fract)
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_collect_water
  
#endif
END MODULE mo_lattrans_interface
