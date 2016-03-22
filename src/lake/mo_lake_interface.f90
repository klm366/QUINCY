!>
!! @brief Contains the interfaces to the lake processes
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2013-04-04)
!!
MODULE mo_lake_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module
  
  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes LAKE

  ! Use process configurations
  dsl4jsb_Use_config(LAKE)

  ! Use process memories
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
!!$  dsl4jsb_Use_memory(ATM2LAND)
  dsl4jsb_Use_memory(LAKE)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_lake_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lake_interface'
  
  !> Type definition for lake process
  TYPE, EXTENDS(t_jsb_process) :: t_lake_process
  CONTAINS
    PROCEDURE :: Register_tasks => Register_lake_tasks        !< Registers tasks for lake process
  END TYPE t_lake_process
  
  !> Type definition for lake_surface task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_lake_surface
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_lake_surface     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_lake_surface  !< Aggregates computed task variables
  END TYPE tsk_lake_surface
  
  !> Constructor interface for lake_surface task
  INTERFACE tsk_lake_surface
    PROCEDURE Create_task_lake_surface                        !< Constructor function for task
  END INTERFACE tsk_lake_surface

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks
  
  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for lake_surface task
  !!
  !! @param[in]     model_id   Model id
  !! @return        return_ptr Instance of process task "lake_surface"
  !!
  FUNCTION Create_task_lake_surface(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_lake_surface::return_ptr)

    return_ptr%name = 'lake_surface'
    return_ptr%process_id = LAKE
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_lake_surface
  
  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for lake process
  !!
  !! @param[in]      model_id Model id  
  !! @param[in,out]  this     Instance of lake process class
  !!
  SUBROUTINE Register_lake_tasks(this, model_id)

    CLASS(t_lake_process), INTENT(inout) :: this
    INTEGER,               INTENT(in)    :: model_id

    CALL this%Register_task(tsk_lake_surface(model_id))

  END SUBROUTINE Register_lake_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "lake_surface"
  !! Task "lake_surface" computes the surface temperature, ice cover and ice depth on lakes.
  !!
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.  
  !! @param[in,out] tile    Tile for which routine is executed.
  !!
  SUBROUTINE update_lake_surface(tile, config, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare process memories
    dsl4jsb_Def_memory(LAKE)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_srf
    dsl4jsb_Real2D_onChunk :: t_wtr
    dsl4jsb_Real2D_onChunk :: t_ice
    dsl4jsb_Real2D_onChunk :: hflx_res
    dsl4jsb_Real2D_onChunk :: ice_fract

    ! Local variables
    INTEGER :: iblk, ics, ice, nc
    REAL(wp) :: dtime

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_lake_surface'

    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(LAKE)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(LAKE)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(LAKE, t_srf)       ! IN
    dsl4jsb_Get_var2D_onChunk(LAKE, t_wtr)       ! IN
    dsl4jsb_Get_var2D_onChunk(LAKE, t_ice)       ! IN
    dsl4jsb_Get_var2D_onChunk(LAKE, hflx_res)    ! IN
    dsl4jsb_Get_var2D_onChunk(LAKE, ice_fract)   ! IN

    ! Aggregate temperatures for open water and ice.
    ! Note: Currently, lake is either ice-free or completely ice-covered.
    t_srf(:) = (1._wp - ice_fract(:)) * t_wtr(:) + ice_fract(:) * t_ice(:)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_lake_surface
  
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "lake_surface"
  !!  
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.  
  !! @param[in,out] tile    Tile for which routine is executed.
  !!
  SUBROUTINE aggregate_lake_surface(tile, config, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    CLASS(t_jsb_config_p),      INTENT(in)    :: config(:)
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_lake_surface'

    ! Local variables
    INTEGER :: iblk, ics, ice, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_lake_surface

#endif
END MODULE mo_lake_interface
