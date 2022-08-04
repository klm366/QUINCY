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

    USE mo_mix_util,          ONLY: find_first_descendant

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    INTEGER  :: iblk, ics, ice, nc
    INTEGER  :: ndesc, nsibl, ic       !< Counter variables for descendent and sibling tiles

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), POINTER :: desc_tile !< Descendant tiles with HYDRO_ process running on them
    CLASS(t_jsb_tile_abstract), POINTER :: sibl_tile !< Sibling tile with MIX_ process running on them

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_hydro_mix'

    ! Declare process memories
    dsl4jsb_Def_config(MIX_)
    dsl4jsb_Def_memory(MIX_)

    dsl4jsb_Def_memory_tile(HYDRO_, desc_tile)
    dsl4jsb_Def_memory_tile(MIX_, sibl_tile)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_storage
    dsl4jsb_Real2D_onChunk :: wtr_inflow_desc_tile
    dsl4jsb_Real2D_onChunk :: wtr_outflow_desc_tile
    dsl4jsb_Real2D_onChunk :: wtr_storage_sibl_tile

    REAL(wp), DIMENSION(options%nc) :: upland_fract     !< Fractions of upland tile
    REAL(wp), DIMENSION(options%nc) :: lowland_fract    !< Fractions of lowland tile
    REAL(wp), DIMENSION(options%nc) :: desc_fract       !< Fractions of descendent tiles with HYDRO_ process
    REAL(wp), DIMENSION(options%nc) :: ref_storage      !< reference storage to distribute water to lowland areas

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

    ! Transport water from uphill to downhill tiles
    IF (dsl4jsb_Config(MIX_)%l_hydro_mix) THEN
      ! This check should be moved to init routine as soon as enough stuff is there to justify one
      IF (tile%process_action(HYDRO_) == ON_TILE_) THEN
        CALL finish(routine, 'HYDRO_ process must be lower in the tile hierarchy than MIX_')
      END IF

      ! Get first child tile with active HYDRO_ processes
      CALL find_first_descendant(HYDRO_, tile, desc_tile)
      ndesc = desc_tile%parent_tile%Get_no_of_children()

      IF (tile%is_upland) THEN

        ! Get area fraction of upland tile
        CALL tile%Get_fraction(ics, ice, iblk, fract=upland_fract(:))

        ! Collect water from HYDRO_ descendant tiles in intermediate mix storage of uphill mix tile
        DO ic = 1, ndesc
          IF (ic > 1) THEN
            ! Get next HYDRO_ tile
            desc_tile => desc_tile%Get_next_sibling_tile()
          END IF
          ! Get area fraction of active HYDRO_ tile
          CALL desc_tile%Get_fraction(ics, ice, iblk, fract=desc_fract(:))
          ! Get memory for water inflow flux at HYDRO_ tile
          dsl4jsb_Get_memory_tile(HYDRO_, desc_tile)
          dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, wtr_outflow, desc_tile)
          ! Move water flowing out of upland tiles into mix water storage
          WHERE(upland_fract(:) > 0._wp)
            wtr_storage(:) = wtr_storage(:) + wtr_outflow_desc_tile(:) * (desc_fract(:) / upland_fract(:))
          END WHERE
        ENDDO

        ! Move water to downhill storage of the same type
        nsibl = tile%parent_tile%Get_no_of_children()
        DO ic = 1, nsibl
          ! Loop over siblings
          IF (ic == 1) THEN
            sibl_tile => tile%parent_tile%Get_first_child_tile()
          ELSE
            sibl_tile => sibl_tile%Get_next_sibling_tile()
          END IF
          ! Identify correct sibling
          ! @todo currently this can be done with the tile type, if more tiles are necessary
          !       we either have to do it via name matching (NAME_up --> NAME_low)
          !       or more explicitly via case statements if 1:1 matches are not sufficient anymore
          IF (sibl_tile%is_lowland) THEN
            CALL sibl_tile%Get_fraction(ics, ice, iblk, fract=lowland_fract(:))
            ! Get memory for water storage at MIX_ downhill tile
            dsl4jsb_Get_memory_tile(MIX_, sibl_tile)
            dsl4jsb_Get_var2d_onChunk_tile_name(MIX_, wtr_storage, sibl_tile)
            ! Move water from uphill to downhill tile
            ! @todo this would be the place to add a lag factor later
            WHERE(lowland_fract(:) > 1.0e-10_wp)
              wtr_storage_sibl_tile(:) = wtr_storage(:) * (upland_fract(:) / lowland_fract(:))
              wtr_storage(:) = 0._wp
            ENDWHERE
          END IF
        END DO

      ELSE IF (tile%is_lowland) THEN
        ! Residual storage to track how much water leaves the mix storage
        ref_storage(:) = wtr_storage(:)
        ! Get area fraction of lowland tile
        CALL tile%Get_fraction(ics, ice, iblk, fract=lowland_fract(:))
        ! Distribute water from intermediate mix storage
        DO ic = 1, ndesc
          IF (ic > 1) THEN
            ! Get next HYDRO_ tile
            desc_tile => desc_tile%Get_next_sibling_tile()
          END IF
          ! Get area fraction of active HYDRO_ tile
          CALL desc_tile%Get_fraction(ics, ice, iblk, fract=desc_fract(:))
          ! Get memory for water inflow flux at HYDRO_ tile
          dsl4jsb_Get_memory_tile(HYDRO_, desc_tile)
          dsl4jsb_Get_var2d_onChunk_tile_name(HYDRO_, wtr_inflow, desc_tile)
          ! Move water from mix water storage into water inflow flux towards lowlands
          WHERE(desc_fract(:) > 1.0e-10_wp)
            wtr_inflow_desc_tile(:) = MAX(0._wp, ref_storage(:) * (lowland_fract(:) / desc_fract(:)))
            wtr_storage(:) = wtr_storage(:) - wtr_inflow_desc_tile(:) * (desc_fract(:) / lowland_fract(:))
          END WHERE
        ENDDO

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
