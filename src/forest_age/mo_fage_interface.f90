!>
!! @brief Forest-age (fage) tasks -- tasks dealing with age cohorts 
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
  USE mo_exception,       ONLY: finish, message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  dsl4jsb_Use_processes FAGE_, CARBON_, DISTURB_, HARV_, PHENO_

  ! Use of process configurations (t_<process>_config)
  dsl4jsb_Use_config(FAGE_)
  dsl4jsb_Use_config(HARV_)

  ! Use of process memories (t_<process>_memory)
  dsl4jsb_Use_memory(FAGE_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(DISTURB_)
  dsl4jsb_Use_memory(HARV_)
  dsl4jsb_Use_memory(PHENO_)
  
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_fage_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_interface'

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_calculate_ageing
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_calculate_ageing    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_calculate_ageing !< Aggregates computed task variables
  END TYPE   tsk_calculate_ageing

  INTERFACE   tsk_calculate_ageing
    PROCEDURE Create_task_calculate_ageing 
  END INTERFACE   tsk_calculate_ageing

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_nat_dist_age_cohorts
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_nat_dist_age_cohorts    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_nat_dist_age_cohorts !< Aggregates computed task variables
  END TYPE   tsk_nat_dist_age_cohorts

  INTERFACE   tsk_nat_dist_age_cohorts
    PROCEDURE Create_task_nat_dist_age_cohorts 
  END INTERFACE   tsk_nat_dist_age_cohorts


  TYPE, EXTENDS(t_jsb_process_task) :: tsk_harvest_age_cohorts
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_harvest_age_cohorts    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_harvest_age_cohorts !< Aggregates computed task variables
  END TYPE   tsk_harvest_age_cohorts

  INTERFACE   tsk_harvest_age_cohorts
    PROCEDURE Create_task_harvest_age_cohorts 
  END INTERFACE   tsk_harvest_age_cohorts

CONTAINS
  ! ================================================================================================================================
  !! Constructors for tasks
  
  !> Constructor for calculate_ageing task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "calculate_ageing"
  !!
  FUNCTION Create_task_calculate_ageing(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_calculate_ageing::return_ptr)
    CALL return_ptr%Construct(name='calculate_ageing', process_id=FAGE_, owner_model_id=model_id)

  END FUNCTION Create_task_calculate_ageing

  !> Constructor for nat_dist_age_cohorts task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "nat_dist_age_cohorts"
  !!
  FUNCTION Create_task_nat_dist_age_cohorts(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_nat_dist_age_cohorts::return_ptr)
    CALL return_ptr%Construct(name='nat_dist_age_cohorts', process_id=FAGE_, owner_model_id=model_id)

  END FUNCTION Create_task_nat_dist_age_cohorts
  
  !> Constructor for harvest_age_cohorts task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "harvest_age_cohorts"
  !!
  FUNCTION Create_task_harvest_age_cohorts(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_harvest_age_cohorts::return_ptr)
    CALL return_ptr%Construct(name='harvest_age_cohorts', process_id=FAGE_, owner_model_id=model_id)

  END FUNCTION Create_task_harvest_age_cohorts


  ! =======================================================================================================
  !> Register tasks for fage process
  !!
  !! @param[in,out] this      Instance of fage process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_fage_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_calculate_ageing(model_id))
    CALL this%Register_task(tsk_nat_dist_age_cohorts(model_id))
    CALL this%Register_task(tsk_harvest_age_cohorts(model_id))

  END SUBROUTINE Register_fage_tasks

  ! ================================================================================================================================
  !>
  !! @brief Task "calculate_ageing" calculates the required changes upon ageing 
  !!        changes in delta fract, delta cpools, delta veg_fract_correction 
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_calculate_ageing(tile, options)

    USE mo_fage_process,         ONLY: update_fract_per_age_upon_ageing, correct_and_synchronise_fracts, &
      & weighted_avg_var_upon_area_change
    USE mo_jsb_time,             ONLY: is_newyear
    USE mo_carbon_interface,     ONLY: recalc_delta_carbon_rp_of_source_and_sink_upon_area_shift, &
      & rescale_carbon_upon_reference_area_change

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    dsl4jsb_Def_memory_tile(PHENO_, prevChild)
    dsl4jsb_Def_memory_tile(PHENO_, currentChild)
    dsl4jsb_Def_memory_tile(CARBON_, prevChild)
    dsl4jsb_Def_memory_tile(CARBON_, currentChild)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: veg_fract_correction_prevChild
    dsl4jsb_Real2D_onChunk :: veg_fract_correction_currentChild
    dsl4jsb_Real2D_onChunk :: lai_prevChild
    dsl4jsb_Real2D_onChunk :: lai_currentChild
    dsl4jsb_Real2D_onChunk :: LAI_yDayMean_prevChild
    dsl4jsb_Real2D_onChunk :: LAI_yDayMean_currentChild
    dsl4jsb_Real2D_onChunk :: LAI_sum_prevChild
    dsl4jsb_Real2D_onChunk :: LAI_sum_currentChild

    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age, cohort_max_age


    CLASS(t_jsb_tile_abstract), POINTER    :: currentChild, prevChild
    REAL(wp), DIMENSION(options%nc)        :: fract_source_max_age, fract_sum, remainder, &
                                            & source_fract, sink_fract, altered_sink_fract, old_veg_fract_correction_sink

    REAL(wp) :: fractSum
    INTEGER  :: i,j, ncohorts, max_age, iblk, ics, ice, nc
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_calculate_ageing'

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
    dsl4jsb_Get_var3D_onChunk(FAGE_, cohort_max_age)
    dsl4jsb_Get_var3D_onChunk(FAGE_, fract_per_age)

    ncohorts = dsl4jsb_Config(FAGE_)%ncohorts
    max_age = dsl4jsb_Config(FAGE_)%max_age

    IF(is_newyear()) THEN  ! year has changed - i.e. age - i.e. shift fractions!
      tile%l_delta_area_change(ics:ice, iblk) = .TRUE.

      ! Note that this update is done globally at the same time step (like the recalculation of maxLAI in PHENO), 
      !  i.e. for different longitudes the updating happens at different seasons etc 

      ! synchronise fract on child tiles and the fract_per_age vector
      currentChild => tile%Get_first_child_tile()
      DO i=1,ncohorts
        ! Area changes upon ageing can have carbon consequences for all children
        currentChild%l_update_c(ics:ice, iblk) = .TRUE.

        CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))

        CALL correct_and_synchronise_fracts(nc, max_age, ncohorts, &
          & i, tile%lcts(1)%lib_id, cohort_max_age(:,:), fract_per_age(:,:), source_fract(:))

        CALL currentChild%Set_fraction(ics, ice, iblk, fract=source_fract(:))

        IF (i .LT. ncohorts) currentChild => currentChild%Get_next_sibling_tile()
      ENDDO

      !Assert that currentChild is last child
      IF(.NOT. currentChild%Is_last_child()) &
        & CALL finish(TRIM(routine), 'Violation of assertion: expected last child after iterating over ncohorts')

      !----- care for shifts among age cohort boundaries
      ! only iterate until 2 - youngest age cohort can not get area from below
      DO i=ncohorts,2,-1
        !TODO-JN-COHORTS maybe work with a getter?
        prevChild => currentChild%Get_prev_sibling_tile()
        CALL prevChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=sink_fract(:))

        ! a) identify delta fraction for shifts from one to the next age cohort
        !TODO-JN-COHORTS needs to be adapted when working with absolute areas

        !TODO-JN-COHORTS: cohort_max_age would only be required with this dimension
        !                 if there would be dynamic changes individual per grid cell, 
        !                 if only static (as currently) it could be implemented easier (maybe in the config of fage?!)
        DO j=1,nc
          fract_source_max_age(j) = fract_per_age(j,cohort_max_age(j,i-1))
        ENDDO

        ! reduce prev tiles delta fraction  
        prevChild%delta_fract(ics:ice, iblk) = prevChild%delta_fract(ics:ice, iblk) - fract_source_max_age(:)

        ! and add them to the current childs delta 
        currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) + fract_source_max_age(:)

        ! b) accordingly move carbon from one age cohort to the next age cohorts delta
        CALL recalc_delta_carbon_rp_of_source_and_sink_upon_area_shift( &
          & prevChild, currentChild, options, fract_source_max_age(:), source_fract(:))

        ! c) correct veg_fract_correction and lai

        ! ------ veg_fract_correction needs to be averaged because it is required for correct recalculation of the carbon pools  
        !     -> particularly of max_green_veg_c for subsequent maxLAI calculation!
        dsl4jsb_Get_memory_tile(PHENO_, prevChild)
        dsl4jsb_Get_memory_tile(PHENO_, currentChild)

        dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, prevChild)
        dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, veg_fract_correction, currentChild)

        old_veg_fract_correction_sink(:) = veg_fract_correction_currentChild(:)
        CALL weighted_avg_var_upon_area_change( sink_fract(:), fract_source_max_age(:), &
          & veg_fract_correction_prevChild(:), veg_fract_correction_currentChild(:)) 

        ! ++ when the veg_fract_correction changes cpools always need to be rescaled 
        IF(ANY(old_veg_fract_correction_sink(:) .NE. veg_fract_correction_currentChild(:))) THEN
          CALL rescale_carbon_upon_reference_area_change(currentChild, options, &
            & old_veg_fract_correction_sink(:), veg_fract_correction_currentChild(:))
        ENDIF

        ! ------ also the lai needs to be weighted upon area change involved with ageing 
        !        this is not done as for the cpools, because currently (28.03.18) KN and JN interprete LAI as per tile
        !        --> NOTE: if not per tile but per canopy it also would need to be rescaled 
        !                  each time the veg_fract_correction changes and could thus be treated the same as the cpools
        dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, lai, prevChild)
        dsl4jsb_Get_var2d_onChunk_tile_name(PHENO_, lai, currentChild)
        CALL weighted_avg_var_upon_area_change( sink_fract(:), fract_source_max_age(:), lai_prevChild(:), lai_currentChild(:)) 

        dsl4jsb_Get_memory_tile(CARBON_, prevChild)
        dsl4jsb_Get_memory_tile(CARBON_, currentChild)

        dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_yDayMean, prevChild)
        dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_yDayMean, currentChild)
        CALL weighted_avg_var_upon_area_change( sink_fract(:), fract_source_max_age(:), &
          & LAI_yDayMean_prevChild(:), LAI_yDayMean_currentChild(:)) 

        dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_sum, prevChild)
        dsl4jsb_Get_var2d_onChunk_tile_name(CARBON_, LAI_sum, currentChild)
        CALL weighted_avg_var_upon_area_change( sink_fract(:), fract_source_max_age(:), &
          & LAI_sum_prevChild(:), LAI_sum_currentChild(:)) 

        !@todo: with flux aggregation - i.e. when more processes are calculated on the leaves, 
        !                                    more variables will need to be considered here
        ! d) TODO-JN-COHORTS what else to shift? 

        ! move further in the child list
        currentChild => prevChild
      ENDDO

      ! update the area per age tracking vector - for ageing this can be done for all children at once
      CALL update_fract_per_age_upon_ageing(nc, dsl4jsb_Config(FAGE_)%max_age, fract_per_age(:,:), mean_age(:))

    ELSE
      tile%l_delta_area_change(ics:ice, iblk) = .FALSE.
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_calculate_ageing


  SUBROUTINE aggregate_calculate_ageing(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_calculate_ageing'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    IF(ANY(tile%l_delta_area_change(ics:ice, iblk))) THEN 
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
    ENDIF
    !No other aggregations required in the calculation of ageing 
    ! (aggregation is done in the CARBON_ process after applying the changes calculated here to carbon variables)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_calculate_ageing


  ! ================================================================================================================================
  !>
  !! @brief Task "nat_dist_age_cohorts" calculates the required changes upon natural disturbances
  !! 
  !! Note: it is only possible to conduct fire and wind damage in one task, because in jsb upon distrbances 
  !!       only the c pool distribution delivered to the target age cohort (first age cohort) is influenced, 
  !!       but not the densities of the source age cohort 
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_nat_dist_age_cohorts(tile, options)

    USE mo_jsb_time,             ONLY: is_newday
    USE mo_fage_util,            ONLY: get_mean_age
    USE mo_fage_process,         ONLY: update_fract_per_age_equally_upon_disturbance
    USE mo_carbon_interface,     ONLY: recalc_delta_carbon_rp_upon_disturbance
    USE mo_carbon_constants,     ONLY: WIND, FIRE

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Def_memory_tile(DISTURB_, currentChild)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: burned_fract
    dsl4jsb_Real2D_onChunk :: damaged_fract

    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age
    dsl4jsb_Real3D_onChunk :: cohort_max_age
    dsl4jsb_Real3D_onChunk :: disturbed_area

    CLASS(t_jsb_tile_abstract), POINTER    :: firstChild, currentChild
    REAL(wp), DIMENSION(options%nc)        :: source_fract, damagedArea, burnedArea

    INTEGER  :: i, ncohorts, max_age, iblk, ics, ice, nc
    LOGICAL  :: l_sourceisfirstChild
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_nat_dist_age_cohorts'

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
    dsl4jsb_Get_var3D_onChunk(FAGE_, cohort_max_age)
    dsl4jsb_Get_var3D_onChunk(FAGE_, disturbed_area)

    ncohorts = dsl4jsb_Config(FAGE_)%ncohorts
    max_age = dsl4jsb_Config(FAGE_)%max_age

    IF(is_newday()) THEN
      tile%l_delta_area_change(ics:ice, iblk) = .TRUE.

      ! step one: determine changes due to disturbances
      firstChild => tile%Get_first_child_tile()
      currentChild => tile%Get_first_child_tile()

      ! in the first iteration of the loop the first child is source and sink 
      l_sourceisfirstChild = .TRUE.
      DO i=1,ncohorts

        ! Area changes upon disturbances only have carbon consequences for the first child
        !  (All other children only have area changes)
        IF(l_sourceisfirstChild) THEN
          currentChild%l_update_c(ics:ice, iblk) = .TRUE.
        ELSE
          currentChild%l_update_c(ics:ice, iblk) = .FALSE.
        ENDIF

        ! -> damaged_fract and burned_fract were determined by the disturbance process and are a fraction of the tile
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))

        !get current child tiles disturbed areas!
        dsl4jsb_Get_memory_tile(DISTURB_, currentChild)
        dsl4jsb_Get_var2d_onChunk_tile(DISTURB_, burned_fract, currentChild)
        dsl4jsb_Get_var2d_onChunk_tile(DISTURB_, damaged_fract, currentChild)

        !***********************************************
        !*** care for fire damage
        burnedArea(:) = burned_fract(:) * source_fract(:)

        ! a) changes in area
        ! reduce current tiles delta fraction  
        currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - burnedArea(:)
        ! and add it to the first childs delta
        firstChild%delta_fract(ics:ice, iblk) = firstChild%delta_fract(ics:ice, iblk) + burnedArea(:)

        ! b) carbon consequences 
        CALL recalc_delta_carbon_rp_upon_disturbance( &
          currentChild, firstChild, l_sourceisfirstChild, options, burnedArea(:), source_fract(:), FIRE)

        !@todo: with flux aggregation - i.e. when more processes are calculated on the leaves, 
        !                                    probably more variables will need to be considered here
        ! c) TODO-JN-COHORTS what else to shift? 


        !***********************************************
        !*** care for wind damage 
        damagedArea(:) = damaged_fract(:) * (source_fract(:) - burnedArea(:))

        ! wind damage only depends on the forcing and there is a chance that no area is damaged
        IF(ANY(damagedArea(:) .GT. 0._wp)) THEN
          ! a) changes in area
          ! reduce current tiles delta fraction  
          currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - damagedArea(:)
          ! and add it to the first childs delta
          firstChild%delta_fract(ics:ice, iblk) = firstChild%delta_fract(ics:ice, iblk) + damagedArea(:)

          ! b) carbon consequences 
          CALL recalc_delta_carbon_rp_upon_disturbance( &
            currentChild, firstChild, l_sourceisfirstChild, options, damagedArea(:), source_fract(:), WIND)

        !@todo: with flux aggregation - i.e. when more processes are calculated on the leaves, 
        !                                    probably more variables will need to be considered here
          ! c) TODO-JN-COHORTS what else to shift? 
        ENDIF


        ! update the area per age tracking vector regarding this child and the first child
        IF (i .GT. 1) THEN ! area moved within the first age cohort requires no update (it stays between year 0 and 1)
          disturbed_area(:,i) = damagedArea(:) + burnedArea(:)
          IF(ANY(disturbed_area(:,i) .GT. 0._wp)) THEN
            IF (dsl4jsb_Config(FAGE_)%disturbance_distribution_scheme == 'equal' ) THEN
              CALL update_fract_per_age_equally_upon_disturbance(nc, max_age, ncohorts, &
                & i, disturbed_area(:,i), cohort_max_age(:,:), fract_per_age(:,:))
            ELSE
              CALL finish(TRIM(routine), 'Currently only "static" implemented as forest ageing scheme')
            ENDIF
          ENDIF
        ENDIF

        ! Move further in the hierarchy 
        currentChild => currentChild%Get_next_sibling_tile()
        l_sourceisfirstChild = .FALSE.
      ENDDO

      !recalculate mean age
      mean_age (:) = get_mean_age(fract_per_age)

    ELSE
      tile%l_delta_area_change(ics:ice, iblk) = .FALSE.
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_nat_dist_age_cohorts

 
  SUBROUTINE aggregate_nat_dist_age_cohorts(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_nat_dist_age_cohorts'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    IF(ANY(tile%l_delta_area_change(ics:ice, iblk))) THEN 
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
    ENDIF
    !No other aggregations required in the calculation of ageing 
    ! (aggregation is done in the CARBON_ process after applying the changes calculated here to carbon variables)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_nat_dist_age_cohorts


  ! ================================================================================================================================
  !>
  !! @brief Task "harvest_age_cohorts" calculates the required changes upon clearcut harvest 
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  ! NOTE: Here its ASSUMED that the task runs after ageing, therefore youngest age cohort is empty 
  SUBROUTINE update_harvest_age_cohorts(tile, options)

    USE mo_jsb_time,             ONLY: is_newyear
    USE mo_fage_util,            ONLY: get_mean_age
    USE mo_fage_process,         ONLY: update_fract_per_age_and_determine_area_to_harvest_from_cohort
    USE mo_carbon_interface,     ONLY: recalc_delta_carbon_rp_upon_disturbance
    USE mo_carbon_constants,     ONLY: CLEARCUT
    USE mo_harvest_constants,    ONLY: HARVEST_OLDEST_TREES

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_config(HARV_)

    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Def_memory(HARV_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: area_to_harvest
    dsl4jsb_Real2D_onChunk :: harvested_area

    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age
    dsl4jsb_Real3D_onChunk :: cohort_max_age

    CLASS(t_jsb_tile_abstract), POINTER    :: firstChild, currentChild
    REAL(wp), DIMENSION(options%nc)        :: source_fract, stillToHarvest, harvestFromThisAgeCohort

    INTEGER  :: i, oldestAgeCohort, ncohorts, max_age, iblk, ics, ice, nc
    LOGICAL  :: l_sourceisfirstChild
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_harvest_age_cohorts'

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

    dsl4jsb_Get_config(HARV_)
    dsl4jsb_Get_memory(HARV_)
    
    dsl4jsb_Get_var2D_onChunk(FAGE_, mean_age)
    dsl4jsb_Get_var3D_onChunk(FAGE_, fract_per_age)
    dsl4jsb_Get_var3D_onChunk(FAGE_, cohort_max_age)

    dsl4jsb_Get_var2D_onChunk(HARV_, area_to_harvest)
    dsl4jsb_Get_var2D_onChunk(HARV_, harvested_area)

    ncohorts = dsl4jsb_Config(FAGE_)%ncohorts
    max_age = dsl4jsb_Config(FAGE_)%max_age

    tile%l_delta_area_change(ics:ice, iblk) = .FALSE.
    IF(is_newyear()) THEN
      harvested_area(:) = 0._wp

      IF(ANY(area_to_harvest(:) .GT. 0._wp)) THEN
        tile%l_delta_area_change(ics:ice, iblk) = .TRUE.

        firstChild => tile%Get_first_child_tile()
        currentChild => tile%Get_first_child_tile()

        !NOTE: Here its ASSUMED that the task runs after ageing, therefore youngest age cohort is empty 
        l_sourceisfirstChild = .FALSE.
        !Assert that youngest cohort is empty 
        CALL firstChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))
        IF (ANY(source_fract(:) .NE. 0._wp)) &
              & CALL finish(TRIM(routine), 'Violation of assertion: area of youngest age cohort expected to be zero.')

        ! Area changes upon clearcut only have carbon density consequences for the first child
        !  (All other children will only have area changes, but keep their densities)
        firstChild%l_update_c(ics:ice, iblk) = .TRUE.

        IF (dsl4jsb_Config(HARV_)%harv_kind .EQ. HARVEST_OLDEST_TREES) THEN
          ! as opposed to damaged_fract in update_nat_dist_age_cohorts area_to_harvest is already with respect
          ! to this tiles area - i.e. not per age cohort, thus it does not have to be multiplied with age cohort fractions
          stillToHarvest(:) = area_to_harvest(:)

          ! TODO-JN-COHORTS: Identify age cohort with oldest available trees for any calculated gridcell
          oldestAgeCohort = ncohorts

          ! Iterate over the children until currentChild points to that age cohort
          DO i=1,oldestAgeCohort-1
            !Assertion that current its not already the last
            IF(currentChild%Is_last_child()) &
              & CALL finish(TRIM(routine), 'Violation of assertion: unexpected last child.')

            currentChild => currentChild%Get_next_sibling_tile()
            currentChild%l_update_c(ics:ice, iblk) = .FALSE.
          ENDDO

          ! then iterate backwards until either all demanded area is fulfilled or the second lowest age cohort is reached
          DO i=oldestAgeCohort,2,-1
            ! Area changes upon clearcut only have carbon density consequences for the first child
            !  (All other children will only have area changes, but keep their densities)
            currentChild%l_update_c(ics:ice, iblk) = .FALSE.

            ! Determine area to be harvested from this age cohort
            harvestFromThisAgeCohort = 0._wp
            CALL update_fract_per_age_and_determine_area_to_harvest_from_cohort(dsl4jsb_Config(HARV_)%harv_kind, nc, &
              & max_age, cohort_max_age(:,i), cohort_max_age(:,i-1), fract_per_age(:,:), &
              & harvestFromThisAgeCohort(:), stillToHarvest(:))
             
            harvested_area(:) = harvested_area(:) + harvestFromThisAgeCohort(:)

            ! reduce current tiles delta fraction  
            currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - harvestFromThisAgeCohort(:)
            ! and add it to the first childs delta
            firstChild%delta_fract(ics:ice, iblk) = firstChild%delta_fract(ics:ice, iblk) + harvestFromThisAgeCohort(:)

            ! b) carbon consequences 
            CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))
            CALL recalc_delta_carbon_rp_upon_disturbance( &
              currentChild, firstChild, l_sourceisfirstChild, options, harvestFromThisAgeCohort(:), source_fract(:), CLEARCUT)

            ! c) TODO-JN-COHORTS what else to shift? 

            ! Go to next child
            currentChild => currentChild%Get_prev_sibling_tile()
          ENDDO

        ELSE
            CALL finish(TRIM(routine),'Not implemented harvest kind')
        ENDIF

        ! recalculate mean age
        mean_age (:) = get_mean_age(fract_per_age)
      ENDIF
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_harvest_age_cohorts

 
  SUBROUTINE aggregate_harvest_age_cohorts(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Def_memory(HARV_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_harvest_age_cohorts'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(FAGE_)
    dsl4jsb_Get_memory(FAGE_)
    dsl4jsb_Get_memory(HARV_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    IF(ANY(tile%l_delta_area_change(ics:ice, iblk))) THEN 
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(HARV_, harvested_area, weighted_by_fract)
    ENDIF
    !No other aggregations required in the calculation of ageing 
    ! (aggregation is done in the CARBON_ process after applying the changes calculated here to carbon variables)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_harvest_age_cohorts

#endif
END MODULE mo_fage_interface
