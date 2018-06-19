!>
!! @brief harvest (harv) tasks 
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
  dsl4jsb_Use_processes HARV_, FM_, FAGE_, CARBON_

  ! Use of process configurations
  dsl4jsb_Use_config(HARV_)
  dsl4jsb_Use_config(FM_)
  dsl4jsb_Use_config(FAGE_)

  ! Use of process memories
  dsl4jsb_Use_memory(HARV_)
  dsl4jsb_Use_memory(FM_)
  dsl4jsb_Use_memory(CARBON_)

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
    USE mo_carbon_process,        ONLY: relocate_carbon_damage
    USE mo_carbon_interface,      ONLY: recalc_per_tile_vars, &
      & calculate_current_c_state_sum, check_carbon_conservation_for_operation

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model
    REAL(wp), DIMENSION(options%nc)           :: old_c_state_sum, dummy1, dummy2, this_fract

    INTEGER :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_forest_harvest'

    ! Declare process memories
    dsl4jsb_Def_memory(HARV_)
    dsl4jsb_Def_memory(FM_)
    dsl4jsb_Def_memory(CARBON_)

    dsl4jsb_Def_config(HARV_)
    dsl4jsb_Def_config(FM_)
    dsl4jsb_Def_config(FAGE_)

    dsl4jsb_Real2D_onChunk :: area_to_harvest
    dsl4jsb_Real2D_onChunk :: harvested_area
    dsl4jsb_Real2D_onChunk :: cconservation_harvest
    dsl4jsb_Real3D_onChunk :: to_harvest

    dsl4jsb_Real2D_onChunk :: c_green
    dsl4jsb_Real2D_onChunk :: c_woods
    dsl4jsb_Real2D_onChunk :: c_reserve
    !
    dsl4jsb_Real2D_onChunk :: c_acid_ag1
    dsl4jsb_Real2D_onChunk :: c_water_ag1
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag1
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag1
    dsl4jsb_Real2D_onChunk :: c_acid_ag2
    dsl4jsb_Real2D_onChunk :: c_water_ag2
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag2
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag2
    ! below ground
    dsl4jsb_Real2D_onChunk :: c_acid_bg1
    dsl4jsb_Real2D_onChunk :: c_water_bg1
    dsl4jsb_Real2D_onChunk :: c_ethanol_bg1
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_bg1
    dsl4jsb_Real2D_onChunk :: c_acid_bg2
    dsl4jsb_Real2D_onChunk :: c_water_bg2
    dsl4jsb_Real2D_onChunk :: c_ethanol_bg2
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_bg2

    dsl4jsb_Real2D_onChunk :: c_humus_1
    dsl4jsb_Real2D_onChunk :: c_humus_2


    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! Get model
    model   => Get_model(tile%owner_model_id)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(HARV_)) RETURN
    IF (.NOT. dsl4jsb_Lctlib_param(ForestFlag)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process configurations
    dsl4jsb_Get_config(HARV_)
    dsl4jsb_Get_config(FM_)
    dsl4jsb_Get_config(FAGE_)

    ! Get process memories
    dsl4jsb_Get_memory(HARV_)
    dsl4jsb_Get_memory(FM_)
    dsl4jsb_Get_memory(CARBON_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(HARV_, area_to_harvest)
    dsl4jsb_Get_var2D_onChunk(HARV_, harvested_area)
    dsl4jsb_Get_var2D_onChunk(HARV_, cconservation_harvest)
    
    dsl4jsb_Get_var3D_onChunk(FM_,   to_harvest)


    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_green)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_woods)              ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_reserve)            ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag1)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag1)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag1)    ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_ag2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_ag2)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_ag2)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_ag2)    ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_bg1)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_bg1)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_bg1)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_bg1)    ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_acid_bg2)          ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_water_bg2)         ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_ethanol_bg2)       ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_nonsoluble_bg2)    ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_humus_1)           ! inout
    dsl4jsb_Get_var2D_onChunk(CARBON_,    c_humus_2)           ! inout

    !
    ! ---------------------------
    ! Go
    IF (is_newyear()) THEN ! only once per year

      IF (TRIM(dsl4jsb_Config(HARV_)%harv_scheme) == 'harvest_area') THEN

        harvested_area(:) = 0._wp
        area_to_harvest(:) = MAX(to_harvest(:,tile%lcts(1)%lib_id), 0._wp)

        IF ( tile%Is_process_active(FAGE_) ) THEN
          ! In this case the carbon cannot be relocated immediately
          ! because harvest now adds area to the first age cohort (controbutes to a 'new' cohort)

          !@todo: maybe for a more general concept it could be considered to always do it like this
          !       i.e. only calculate the damage here and then apply it later maybe in a carbon task?
          IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), &
            & TRIM(tile%name)//' - Run with AGE COHORTS: Delayed relocation of harvested carbon.')

        ELSE
          CALL calculate_current_c_state_sum(tile, options, old_c_state_sum(:))

          ! -> for now I use the relocate_carbon_damage

          ! TODO-JN: implement own function using harvest pools and area 
          !    relative to parent/ or absolute area, not relative to this_fract?!
          !          CALL tile%Get_fraction(ics, ice, iblk, fract=this_fract(:))
          !          harvested_area(:)=MIN(area_to_harvest(:), this_fract(:)) 
          !-> in relocate_carbon_damage it is the area relative to this pft 
          !   and here area_to_harvest is already relative to the pft!  
          harvested_area(:)=area_to_harvest(:)

          dummy1(:) = 0._wp
          dummy2(:) = 0._wp

          CALL relocate_carbon_damage(       &
            & harvested_area(:),             & ! in
            & dsl4jsb_Lctlib_param(LeafLit_coef(1:5)), & ! in
            & dsl4jsb_Lctlib_param(WoodLit_coef(1:5)), & ! in
            & c_green(:),                & ! inout
            & c_woods(:),                & ! inout
            & c_reserve(:),              & ! inout
            & dummy1(:),                     & ! inout
            & dummy2(:),                     & ! inout
            & c_acid_ag1(:),            & ! inout
            & c_water_ag1(:),           & ! inout
            & c_ethanol_ag1(:),         & ! inout
            & c_nonsoluble_ag1(:),      & ! inout
            & c_acid_ag2(:),            & ! inout
            & c_water_ag2(:),           & ! inout
            & c_ethanol_ag2(:),         & ! inout
            & c_nonsoluble_ag2(:),      & ! inout
            & c_acid_bg1(:),            & ! inout
            & c_water_bg1(:),           & ! inout
            & c_ethanol_bg1(:),         & ! inout
            & c_nonsoluble_bg1(:),      & ! inout
            & c_acid_bg2(:),            & ! inout
            & c_water_bg2(:),           & ! inout
            & c_ethanol_bg2(:),         & ! inout
            & c_nonsoluble_bg2(:),      & ! inout
            & c_humus_1(:),             & ! inout
            & c_humus_2(:)              & ! inout
            & )

          CALL recalc_per_tile_vars(tile, options, &
            & c_green= c_green(:), c_woods = c_woods(:), c_reserve = c_reserve(:), &
            & c_acid_ag1 = c_acid_ag1(:), c_water_ag1 = c_water_ag1(:), &
            & c_ethanol_ag1 = c_ethanol_ag1(:), c_nonsoluble_ag1 = c_nonsoluble_ag1(:), &
            & c_acid_ag2 = c_acid_ag2(:), c_water_ag2 = c_water_ag2(:), &
            & c_ethanol_ag2 = c_ethanol_ag2(:), c_nonsoluble_ag2 = c_nonsoluble_ag2(:), &
            & c_acid_bg1 = c_acid_bg1(:), c_water_bg1 = c_water_bg1(:), &
            & c_ethanol_bg1 = c_ethanol_bg1(:), c_nonsoluble_bg1 = c_nonsoluble_bg1(:), &
            & c_acid_bg2 = c_acid_bg2(:), c_water_bg2 = c_water_bg2(:), &
            & c_ethanol_bg2 = c_ethanol_bg2(:), c_nonsoluble_bg2 = c_nonsoluble_bg2(:), &
            & c_humus_1 = c_humus_1(:), c_humus_2 = c_humus_2(:))

          ! wind has no extra flux to atmos, only redistribution among cpools - therefore cflux = 0._wp 
          dummy1 = 0._wp
          CALL check_carbon_conservation_for_operation(tile, options, old_c_state_sum(:), &
            & dummy1(:), cconservation_harvest(:))

        ENDIF

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

    dsl4jsb_Aggregate_onChunk(HARV_, area_to_harvest, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_forest_harvest

#endif
END MODULE mo_harv_interface
