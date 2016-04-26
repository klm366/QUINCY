! ===============================================================================================================================
! === THIS IS A TEMPLATE. PLEASE REPLACE ...                                                                                  ===
! === ... <FIRST NAME LAST NAME>    with your name e.g. "Martin Mustermann"                                                   ===
! === ... <DATE>                    with the present date in YYYY-MM-DD format e.g. "2016-04-30"                              ===
! === ... <PROCESS_NAME_UPPER_CASE> with your process name e.g. "ALBEDO"                                                      ===
! === ... <PROCESS_NAME_LOWER_CASE> with your process name e.g. "albedo"                                                      ===
! === ... <EXPLANATIION_OF_PROCESS> with a short explanation line e.g.                                                        ===
! ===     "computes the surface temperature from latent and sensible heat flux, net radiation and below ground sensible heat" ===
! ===                                                                                                                         ===
! === Then go through the code line by line to adapt it to your needs:                                                        ===
! === !X marks lines with examples (e.g:), templates (if necessary:) and implementations points (Implementation:) that you    ===
! === have to adapt to your process.                                                                                          ===
! ===                                                                                                                         ===
! === Finally delete this header.                                                                                             ===
! ===============================================================================================================================



!>
!! @brief Contains the interfaces to the <PROCESS_NAME_LOWER_CASE> process
!!
!! @author
!!  <FIRST NAME LAST NAME>, MPI-M Hamburg
!!
!! @par Revision History
!! First version by <FIRST NAME LAST NAME> (<DATE>)
!!
MODULE mo_<PROCESS_NAME_LOWER_CASE>_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes <PROCESS_NAME_UPPER_CASE>

  ! Use of process configurations
  dsl4jsb_Use_config(<PROCESS_NAME_UPPER_CASE>)

  ! Use of process memories
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(<PROCESS_NAME_UPPER_CASE>)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_<PROCESS_NAME_LOWER_CASE>_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_<PROCESS_NAME_LOWER_CASE>_interface'
  
  !> Type definition for <PROCESS_NAME_LOWER_CASE> process
  TYPE, EXTENDS(t_jsb_process) :: t_<PROCESS_NAME_LOWER_CASE>_process
  CONTAINS
    PROCEDURE :: Register_tasks => Register_<PROCESS_NAME_LOWER_CASE>_tasks        !< Registers tasks for <PROCESS_NAME_LOWER_CASE> process
  END TYPE t_<PROCESS_NAME_LOWER_CASE>_process
  
  !> Type definition for <PROCESS_NAME_LOWER_CASE> task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_<PROCESS_NAME_LOWER_CASE>
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_<PROCESS_NAME_LOWER_CASE>     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_<PROCESS_NAME_LOWER_CASE>  !< Aggregates computed task variables
  END TYPE tsk_<PROCESS_NAME_LOWER_CASE>
  
  !> Constructor interface for <PROCESS_NAME_LOWER_CASE> task
  INTERFACE tsk_<PROCESS_NAME_LOWER_CASE>
    PROCEDURE Create_task_<PROCESS_NAME_LOWER_CASE>                        !< Constructor function for task
  END INTERFACE tsk_<PROCESS_NAME_LOWER_CASE>

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for <PROCESS_NAME_LOWER_CASE> task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "<PROCESS_NAME_LOWER_CASE>"
  !!
  FUNCTION Create_task_<PROCESS_NAME_LOWER_CASE>(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_<PROCESS_NAME_LOWER_CASE>::return_ptr)

    return_ptr%name = '<PROCESS_NAME_LOWER_CASE>'
    return_ptr%process_id = <PROCESS_NAME_UPPER_CASE>
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_<PROCESS_NAME_LOWER_CASE>
  
  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for <PROCESS_NAME_LOWER_CASE> process
  !!
  !! @param[in,out] this      Instance of <PROCESS_NAME_LOWER_CASE> process class  
  !! @param[in]     model_id  Model id  
  !!
  SUBROUTINE Register_<PROCESS_NAME_LOWER_CASE>_tasks(this, model_id)

    CLASS(t_<PROCESS_NAME_LOWER_CASE>_process), INTENT(inout) :: this
    INTEGER,               INTENT(in)    :: model_id

    CALL this%Register_task(tsk_<PROCESS_NAME_LOWER_CASE>(model_id))

  END SUBROUTINE Register_<PROCESS_NAME_LOWER_CASE>_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "<PROCESS_NAME_LOWER_CASE>"
  !! Task "<PROCESS_NAME_LOWER_CASE>" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_<PROCESS_NAME_LOWER_CASE>(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare process memories
    dsl4jsb_Def_memory(<PROCESS_NAME_UPPER_CASE>)

    ! Declare pointers to variables in memory
    !X e.g: dsl4jsb_Real2D_onChunk :: t_srf

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    !X if necessary: REAL(wp) :: dtime, steplen

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_<PROCESS_NAME_LOWER_CASE>'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    !X if necessary: dtime   = options%dtime
    !X if necessary: steplen = options%steplen
    
    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(<PROCESS_NAME_UPPER_CASE>)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    
    ! Get process memories
    dsl4jsb_Get_memory(<PROCESS_NAME_UPPER_CASE>)

    ! Get process variables
    !X e.g: dsl4jsb_Get_var2D_onChunk(<PROCESS_NAME_UPPER_CASE>, t_srf)     ! IN

    !X Implementation: Start your process scheme here...
    
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
    
  END SUBROUTINE update_<PROCESS_NAME_LOWER_CASE>
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "<PROCESS_NAME_LOWER_CASE>"
  !! 
  !! @param[in,out] tile    Tile for which routine is executed.  
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.  
  !!
  SUBROUTINE aggregate_<PROCESS_NAME_LOWER_CASE>(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_<PROCESS_NAME_LOWER_CASE>'

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
    !X Implementation: Start your process scheme here...
    
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_<PROCESS_NAME_LOWER_CASE>

  
  ! ================================================================================================================================
  !>
  !! @brief Another task for "<PROCESS_NAME_LOWER_CASE>"
  !!
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.  
  !! @param[in,out] tile    Tile for which routine is executed.
  !  ....
  
#endif
END MODULE mo_<PROCESS_NAME_LOWER_CASE>_interface
