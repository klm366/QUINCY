!>
!! @brief Contains the interfaces to the turbulence process
!!
!! @author
!!  Reiner Schnur, MPI-M Hamburg
!!
!! @par Revision History
!! First version by Reiner Schnur (2016-04-01)
!!
MODULE mo_turb_interface
#ifndef __NO_JSBACH__

#include "dsl_jsbach.inc"
  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract
  USE mo_jsb_config_class,   ONLY: t_jsb_config
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes TURB, HYDRO, SOIL

  ! Use process configurations
  dsl4jsb_Use_config(TURB)

  ! Use process memories
  ! dsl4jsb_Use_memory(ATM2LAND)
  USE mo_atmland_memory_class, ONLY: t_atm2land_memory
  dsl4jsb_Use_memory(TURB)
  dsl4jsb_Use_memory(HYDRO)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_turb_tasks !,t_turb_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_turb_interface'

!!$  !> Type definition for turbulence process
!!$  TYPE, EXTENDS(t_jsb_process) :: t_turb_process
!!$  CONTAINS
!!$    PROCEDURE :: Register_tasks => Register_turb_tasks        !< Registers tasks for turbulence process
!!$  END TYPE t_turb_process

  !> Type definition for surface_roughness task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_roughness
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_roughness     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_roughness  !< Aggregates computed task variables
  END TYPE tsk_surface_roughness

  !> Constructor interface for surface_roughness task
  INTERFACE tsk_surface_roughness
    PROCEDURE Create_task_surface_roughness                        !< Constructor function for task
  END INTERFACE tsk_surface_roughness

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_roughness task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_surface_roughness"
  !!
  FUNCTION Create_task_surface_roughness(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_roughness::return_ptr)

    return_ptr%name = 'surface_roughness'
    return_ptr%process_id = TURB
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_surface_roughness

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for turbulence process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of TURB process class
  !!
  SUBROUTINE Register_turb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_surface_roughness(model_id))

  END SUBROUTINE Register_Turb_Tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "surface_roughness"
  !! Task "<PROCESS_NAME_LOWER_CASE>" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_roughness(tile, options)

    USE mo_turb_constants, ONLY: blending_height, roughness_snow
    USE mo_jsb_lct_class,  ONLY: LAKE_TYPE, Contains_lct

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Declare process memories
    dsl4jsb_Def_memory(TURB)
    dsl4jsb_Def_memory(HYDRO)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: rough_m
    dsl4jsb_Real2D_onChunk :: rough_h
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_snow_ice

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_roughness'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(TURB)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process memories
    dsl4jsb_Get_memory(TURB)
    dsl4jsb_Get_memory(HYDRO)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(TURB,  rough_m)         ! IN
    dsl4jsb_Get_var2D_onChunk(TURB,  rough_h)         ! OUT

    ! RMS TBD: This has to be revised!
    IF (tile%Is_process_active(SOIL)) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO, fract_snow_soil) ! IN
      rough_h(:) =          fract_snow_soil(:)  / (LOG(blending_height /        roughness_snow)  ** 2)  &
        &        + (1._wp - fract_snow_soil(:)) / (LOG(blending_height / MIN(1._wp, rough_m(:))) ** 2)
    ELSE IF (Contains_lct(tile%lcts, LAKE_TYPE)) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO, fract_snow_ice) ! IN
      rough_h(:) =          fract_snow_ice(:)  / (LOG(blending_height /        roughness_snow)  ** 2)  &
        &        + (1._wp - fract_snow_ice(:)) / (LOG(blending_height / MIN(1._wp, rough_m(:))) ** 2)
    ELSE
      CALL finish(TRIM(routine), 'Don"t know what to do.')
    END IF

    ! Limit roughness length for heat by 1 m (see below Eq.5.6 in ECHAM5 manual)
    rough_h(:) = MIN(1._wp, blending_height * EXP(-1._wp / SQRT(rough_h(:))))

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_roughness

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_roughness"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_roughness(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_roughness'

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    !X Implementation: Start your process scheme here...

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_roughness

  ! ================================================================================================================================
  !>
  !! @brief Another task for "<PROCESS_NAME_LOWER_CASE>"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !  ....

#endif
END MODULE mo_turb_interface
