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
  USE mo_exception,       ONLY: finish, message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: ON_TILE_
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes MIX_, HYDRO_

  ! Use of process configurations
  dsl4jsb_Use_config(MIX_)

  ! Use of process memories
  dsl4jsb_Use_memory(MIX_)
  dsl4jsb_Use_memory(HYDRO_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables 

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_mix_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_mix_interface'

  !> Type definition for reset_mix task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_mix
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_mix     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_mix  !< Aggregates computed task variables
  END TYPE tsk_reset_mix

  !> Constructor interface for reset_mix task
  INTERFACE tsk_reset_mix
    PROCEDURE Create_task_reset_mix                        !< Constructor function for task
  END INTERFACE tsk_reset_mix

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
  !!
  !> Constructor for reset_mix task
  FUNCTION Create_task_reset_mix(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_mix::return_ptr)
    CALL return_ptr%Construct(name='reset_mix', process_id=MIX_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_mix

  !> Constructor for hydro_mix task
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

    CALL process%Register_task(tsk_reset_mix(model_id))
    CALL process%Register_task(tsk_hydro_mix(model_id))

  END SUBROUTINE Register_mix_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "reset_mix"
  !!
  SUBROUTINE update_reset_mix(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_mix'

    ! Declare process memories
    dsl4jsb_Def_memory(MIX_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_latflow

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(MIX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(MIX_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(MIX_,   wtr_latflow)   ! inout

    ! Reset lateral intertile fluxes prior to update_*_mix
    wtr_latflow(:) = 0._wp

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_reset_mix

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "reset_mix"
  !!
  SUBROUTINE aggregate_reset_mix(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(MIX_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_mix'

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(MIX_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(MIX_,   wtr_latflow,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_reset_mix

  ! ================================================================================================================================
  !>
  !! @brief Implementation of "update" for task "hydro_mix"
  !!
  SUBROUTINE update_hydro_mix(tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    INTEGER  :: nsibl, ic              !< Counter variables for sibling tiles
    REAL(wp) :: dtime

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: sibl_tile !< Sibling tile with MIX_ process running on them

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_hydro_mix'

    ! Declare process memories
    dsl4jsb_Def_config(MIX_)
    dsl4jsb_Def_memory(MIX_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Def_memory_tile(MIX_,   sibl_tile)
    dsl4jsb_Def_memory_tile(HYDRO_, sibl_tile)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_latflow
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: weq_pond
    dsl4jsb_Real2D_onChunk :: weq_pond_sibl_tile
    dsl4jsb_Real2D_onChunk :: wtr_latflow_sibl_tile

    REAL(wp), DIMENSION(options%nc) :: source_fract    !< Fractions of water source tile
    REAL(wp), DIMENSION(options%nc) :: target_fract    !< Fractions of water target tile

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(MIX_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(MIX_)

    ! Get process memories
    dsl4jsb_Get_memory(MIX_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(MIX_,   wtr_latflow)   ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)        ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)      ! inout

    ! Transport water from between tiles
    IF (dsl4jsb_Config(MIX_)%l_hydro_mix) THEN
      ! This check should be moved to init routine as soon as enough stuff is there to justify one
      IF (tile%process_action(HYDRO_) /= ON_TILE_) THEN
        CALL finish(routine, 'MIX_ processes must be activate on the same tiles as HYDRO_')
      END IF

      ! Downhill routing of surface runoff from ridges to depressions as well as from
      !   upland depressions to lowland depressions
      IF (tile%contains_upland .OR. tile%contains_ridge) THEN

        ! Get area fraction of upland tile
        CALL tile%Get_fraction(ics, ice, iblk, fract=source_fract(:))

        ! Move runoff into water storage of the target tile
        nsibl = tile%parent_tile%Get_no_of_children()
        DO ic = 1, nsibl
          ! Loop over all tiles in this hierarchy level
          IF (ic == 1) THEN
            sibl_tile => tile%parent_tile%Get_first_child_tile()
          ELSE
            sibl_tile => sibl_tile%Get_next_sibling_tile()
          END IF
          ! Identify correct target sibling
          IF ( (tile%contains_upland .AND. tile%contains_ridge .AND.                 &
            &   sibl_tile%contains_upland .AND. sibl_tile%contains_depression) .OR.  &
            &  (tile%contains_upland .AND. tile%contains_depression .AND.            &
            &   sibl_tile%contains_lowland .AND. sibl_tile%contains_depression) .OR. &
            &  (tile%contains_lowland .AND. tile%contains_ridge .AND.                &
            &   sibl_tile%contains_lowland .AND. sibl_tile%contains_depression) ) THEN

            CALL sibl_tile%Get_fraction(ics, ice, iblk, fract=target_fract(:))
            ! Get memory for water storage at MIX_ target tile
            dsl4jsb_Get_memory_tile(MIX_,   sibl_tile)
            dsl4jsb_Get_memory_tile(HYDRO_, sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(MIX_,   wtr_latflow, sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, weq_pond,    sibl_tile)
            ! Move surface runoff from source to target tile
            ! @todo this would be the place to add a lag factor later. Then also
            !       the lateral flow needs to be updated as difference between new and old storage
            WHERE(target_fract(:) > 1.0e-10_wp)
              weq_pond_sibl_tile(:)    = weq_pond_sibl_tile(:) + &
                & runoff(:) * dtime / rhoh2o * (source_fract(:) / target_fract(:))
              wtr_latflow_sibl_tile(:) = wtr_latflow_sibl_tile(:) + runoff(:) * (source_fract(:) / target_fract(:))
              wtr_latflow(:)           = wtr_latflow(:) - runoff(:)
              runoff(:)                = 0._wp
            ENDWHERE
          END IF
        END DO
      END IF

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
    dsl4jsb_Def_memory(HYDRO_)

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
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(MIX_,   wtr_latflow,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,     weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_hydro_mix

#endif
END MODULE mo_mix_interface
