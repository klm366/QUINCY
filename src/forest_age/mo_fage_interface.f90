!>
!! @brief Contains the interfaces to the forest-age processes
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

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_nat_disturb_cohort_changes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_nat_disturb_cohort_changes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_nat_disturb_cohort_changes !< Aggregates computed task variables
  END TYPE   tsk_nat_disturb_cohort_changes

  INTERFACE   tsk_nat_disturb_cohort_changes
    PROCEDURE Create_task_nat_disturb_cohort_changes 
  END INTERFACE   tsk_nat_disturb_cohort_changes


  TYPE, EXTENDS(t_jsb_process_task) :: tsk_harvest_cohort_changes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_harvest_cohort_changes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_harvest_cohort_changes !< Aggregates computed task variables
  END TYPE   tsk_harvest_cohort_changes

  INTERFACE   tsk_harvest_cohort_changes
    PROCEDURE Create_task_harvest_cohort_changes 
  END INTERFACE   tsk_harvest_cohort_changes

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

    return_ptr%name = 'calculate_ageing'
    return_ptr%process_id = FAGE_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_calculate_ageing

  !> Constructor for nat_disturb_cohort_changes task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "nat_disturb_cohort_changes"
  !!
  FUNCTION Create_task_nat_disturb_cohort_changes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_nat_disturb_cohort_changes::return_ptr)

    return_ptr%name = 'nat_disturb_cohort_changes'
    return_ptr%process_id = FAGE_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_nat_disturb_cohort_changes
  
  !> Constructor for harvest_cohort_changes task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "harvest_cohort_changes"
  !!
  FUNCTION Create_task_harvest_cohort_changes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_harvest_cohort_changes::return_ptr)

    return_ptr%name = 'harvest_cohort_changes'
    return_ptr%process_id = FAGE_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_harvest_cohort_changes


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
    CALL this%Register_task(tsk_nat_disturb_cohort_changes(model_id))
    CALL this%Register_task(tsk_harvest_cohort_changes(model_id))

  END SUBROUTINE Register_fage_tasks

  ! ================================================================================================================================
  !>
  !! @brief Task "calculate_ageing" calculates the required changes upon ageing 
  !!        changes in delta fract, delta cpools, delta veg_frac_correction 
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_calculate_ageing(tile, options)

    USE mo_fage_process,         ONLY: update_frac_per_age_upon_ageing, correct_and_synchronise_fracts, &
      & weighted_avg_var_upon_area_change
    USE mo_jsb_time,             ONLY: is_newyear
    USE mo_carbon_interface,     ONLY: update_delta_carbon_rp_of_source_and_sink_tile_upon_area_shift, &
      & rescale_carbon_upon_reference_area_change

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    TYPE(t_pheno_memory), POINTER :: pheno__mem__source
    TYPE(t_pheno_memory), POINTER :: pheno__mem__sink

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: veg_frac_correction_source
    dsl4jsb_Real2D_onChunk :: veg_frac_correction_sink
    dsl4jsb_Real2D_onChunk :: lai_source
    dsl4jsb_Real2D_onChunk :: lai_sink


    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age, cohort_max_age


    CLASS(t_jsb_tile_abstract), POINTER    :: currentChild, prevChild
    REAL(wp), DIMENSION(options%nc)        :: fract_source_max_age, fract_sum, remainder, &
                                            & source_fract, sink_fract, altered_sink_fract, old_veg_frac_correction_sink

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
      dsl4jsb_Config(FAGE_)%l_has_area_changes = .TRUE.

      ! Note that this update is done globally at the same time step (like the recalculation of maxLAI in PHENO), 
      !  i.e. for different longitudes the updating happens at different local times 

      ! Area changes upon ageing can have carbon consequences for all children
      dsl4jsb_Config(FAGE_)%l_only_update_c_of_first_child = .FALSE.

      ! synchronise fract on child tiles and the fract_per_age vector
      currentChild => tile%Get_first_child_tile()
      DO i=1,ncohorts
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))

        CALL correct_and_synchronise_fracts(nc, max_age, ncohorts, &
          & i, tile%lcts(1)%lib_id, cohort_max_age(:,:), fract_per_age(:,:), source_fract(:))

        CALL currentChild%Set_fraction(ics, ice, iblk, fract=source_fract(:))

        IF (i .LT. ncohorts) currentChild => currentChild%Get_next_sibling_tile()
      ENDDO

      !Assert that currentChild is last child
      IF(.NOT. currentChild%Is_last_child()) &
        & CALL finish(TRIM(routine), 'Violation of assertion: expected last child after iterating over ncohorts')

      !----- care for shifts among cohort boundaries
      ! only iterate until 2 - youngest cohort can not get area from below
      DO i=ncohorts,2,-1
        !TODO-JN-COHORTS maybe work with a getter?
        prevChild => currentChild%Get_prev_sibling_tile()
        CALL prevChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=sink_fract(:))

        ! a) identify delta fraction for shifts from one to the next cohort
        !TODO-JN needs to be adapted when working with absolute areas

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

        ! b) accordingly move carbon from one cohort to the next cohorts delta
        CALL update_delta_carbon_rp_of_source_and_sink_tile_upon_area_shift( &
          & prevChild, currentChild, options, fract_source_max_age(:), source_fract(:))

        ! c) correct veg_frac_correction - required for correct recalculation the carbon pools  
        !     -> particularly max_green_veg_c for subsequent maxLAI calculation!
        SELECT TYPE (mm => prevChild%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem__source => mm ; END SELECT
        SELECT TYPE (mm => currentChild%mem(pheno_)%p); TYPE IS (t_pheno_memory); pheno__mem__sink => mm ; END SELECT
        veg_frac_correction_source => pheno__mem__source%veg_frac_correction%ptr(ics:ice,iblk)  ! in
        veg_frac_correction_sink => pheno__mem__sink%veg_frac_correction%ptr(ics:ice,iblk)      ! inout

        old_veg_frac_correction_sink(:) = veg_frac_correction_sink(:)
        CALL weighted_avg_var_upon_area_change( sink_fract(:), fract_source_max_age(:), &
          & veg_frac_correction_source(:), veg_frac_correction_sink(:)) 

        IF(ANY(old_veg_frac_correction_sink(:) .NE. veg_frac_correction_sink(:))) THEN
          CALL rescale_carbon_upon_reference_area_change(currentChild, options, &
            & old_veg_frac_correction_sink(:), veg_frac_correction_sink(:))
        ENDIF

        ! TODO-JN: hm, do this also has to be weighted with the new veg_frac_corr? If so it also would need to be
        ! rescaled each time the veg_frac_correction changes!?! - Kim and I currently (28.03.18) interprete LAI as per tile
        lai_source => pheno__mem__source%lai%ptr(ics:ice,iblk)  ! in
        lai_sink => pheno__mem__sink%lai%ptr(ics:ice,iblk)      ! inout
        CALL weighted_avg_var_upon_area_change( sink_fract(:), fract_source_max_age(:), &
          & lai_source(:), lai_sink(:)) 

        ! d) TODO-JN-COHORTS what else to shift? 

        ! move further in the child list
        currentChild => prevChild
      ENDDO

      ! update the area per age tracking vector - for ageing this can be done for all children at once
      CALL update_frac_per_age_upon_ageing(nc, dsl4jsb_Config(FAGE_)%max_age, fract_per_age(:,:), mean_age(:))

    ELSE
      dsl4jsb_Config(FAGE_)%l_has_area_changes = .FALSE.
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

    IF(dsl4jsb_Config(FAGE_)%l_has_area_changes) THEN 
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
    ENDIF
    !No other aggregations required in the calculation of ageing 
    ! (aggregation is done in the CARBON_ process after applying the changes calculated here to carbon variables)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_calculate_ageing


  ! ================================================================================================================================
  !>
  !! @brief Task "nat_disturb_cohort_changes" calculates the required changes upon natural disturbances
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_nat_disturb_cohort_changes(tile, options)

    USE mo_jsb_time,             ONLY: is_newday
    USE mo_fage_util,            ONLY: get_mean_age
    USE mo_fage_process,         ONLY: update_frac_per_age_equally_upon_disturbance
    USE mo_carbon_interface,     ONLY: update_delta_carbon_rp_upon_disturbance
    USE mo_carbon_constants,     ONLY: WIND, FIRE

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Def_memory(DISTURB_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: burned_frac
    dsl4jsb_Real2D_onChunk :: damaged_frac

    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age
    dsl4jsb_Real3D_onChunk :: cohort_max_age

    CLASS(t_jsb_tile_abstract), POINTER    :: firstChild, currentChild
    REAL(wp), DIMENSION(options%nc)        :: source_fract, damagedArea, burnedArea, disturbedArea, helper1, helper2

    INTEGER  :: i, ncohorts, max_age, iblk, ics, ice, nc
    LOGICAL  :: l_sourceisfirstChild
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_nat_disturb_cohort_changes'

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

    ncohorts = dsl4jsb_Config(FAGE_)%ncohorts
    max_age = dsl4jsb_Config(FAGE_)%max_age

    IF(is_newday()) THEN
      dsl4jsb_Config(FAGE_)%l_has_area_changes = .TRUE.

      ! Area changes upon disturbances only have carbon consequences for the first child
      !  (All other children only have area changes)
      dsl4jsb_Config(FAGE_)%l_only_update_c_of_first_child = .TRUE.

      ! step one: determine changes due to disturbances
      firstChild => tile%Get_first_child_tile()
      currentChild => tile%Get_first_child_tile()

      l_sourceisfirstChild = .TRUE.

      DO i=1,ncohorts

        ! -> damaged_frac and burned_frac were determined by the disturbance process and are a fraction of the tile
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))

        !get current child tiles disturbed areas!
        SELECT TYPE (mm => currentChild%mem(disturb_)%p); TYPE IS (t_disturb_memory); disturb__mem => mm ; END SELECT
        burned_frac => disturb__mem%burned_frac%ptr(ics:ice,iblk)
        damaged_frac => disturb__mem%damaged_frac%ptr(ics:ice,iblk)

!print *, i, burned_frac(:)

        !TODO-JN-COHORTS - think about absolute vs relative fracts

!TODO: for fire -- ah, vielleicht zeitlich erst feuer, weil so in dist? obwohl is arbitrary, oder?!
! ++ auch da über ta, ra und densities denken  
! ++ generalisieren?
! ++ fire calculation is switched off currently in disturbances! 

        !TODO-JN: ich glaube nicht, dass ich das brauche -- wieder rausnehmen, denn sonnst müsste ich:
                      !TODO: if only > something than feed back into damaged frac!  mit aggregate und allem!
        !*** care for wind damage -     TODO-JN: check: work with max area to prevent numerical issues!?
        helper2(:) = 0._wp
        helper1(:) = damaged_frac(:) * source_fract(:)
        damagedArea(:) = MERGE(helper1, helper2, helper1>=1.e-10_wp) !TODO-JN make it a constant
!        damagedArea(:) = damaged_frac(:) * source_fract(:)

        ! wind damage only depends on the forcing and there is a chance that no area is damaged
        IF(ANY(damagedArea(:) .GT. 0._wp)) THEN
          ! a) changes in area
          ! reduce current tiles delta fraction  
          currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - damagedArea(:)
          ! and add it to the first childs delta
          firstChild%delta_fract(ics:ice, iblk) = firstChild%delta_fract(ics:ice, iblk) + damagedArea(:)

          ! b) carbon consequences 
          CALL update_delta_carbon_rp_upon_disturbance( &
            currentChild, firstChild, l_sourceisfirstChild, options, damagedArea(:), source_fract(:), WIND)

          ! c) TODO-JN-COHORTS what else to shift? 
        ENDIF

        !*** care for fire damage
        !TODO-JN-COHORTS - ? reduce the available area by what has already been damaged? 
        burnedArea(:) = burned_frac(:) * (source_fract(:) - damagedArea(:))

        ! a) changes in area
        ! reduce current tiles delta fraction  
        currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - burnedArea(:)
        ! and add it to the first childs delta
        firstChild%delta_fract(ics:ice, iblk) = firstChild%delta_fract(ics:ice, iblk) + burnedArea(:)

        ! b) carbon consequences 
        !TODO-JN

        ! c) TODO-JN-COHORTS what else to shift? 

        ! update the area per age tracking vector regarding this child and the first child
        IF (i .GT. 1) THEN ! area moved within the first cohort requires no update (it stays between year 0 and 1)
          disturbedArea(:) = damagedArea(:) + burnedArea(:)
          IF(ANY(disturbedArea(:) .GT. 0._wp)) THEN
            IF (dsl4jsb_Config(FAGE_)%disturbance_distribution_scheme == 'equal' ) THEN
              CALL update_frac_per_age_equally_upon_disturbance(nc, max_age, ncohorts, &
                & i, disturbedArea(:), cohort_max_age(:,:), fract_per_age(:,:))
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
      dsl4jsb_Config(FAGE_)%l_has_area_changes = .FALSE.
    ENDIF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_nat_disturb_cohort_changes

 
  SUBROUTINE aggregate_nat_disturb_cohort_changes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_nat_disturb_cohort_changes'

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

    IF(dsl4jsb_Config(FAGE_)%l_has_area_changes) THEN 
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
    ENDIF
    !No other aggregations required in the calculation of ageing 
    ! (aggregation is done in the CARBON_ process after applying the changes calculated here to carbon variables)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_nat_disturb_cohort_changes


  ! ================================================================================================================================
  !>
  !! @brief Task "harvest_cohort_changes" calculates the required changes upon clearcut harvest 
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  ! NOTE: Here its ASSUMED that the task runs after ageing, therefore youngest cohort is empty 
  SUBROUTINE update_harvest_cohort_changes(tile, options)

    USE mo_jsb_time,             ONLY: is_newyear
    USE mo_fage_util,            ONLY: get_mean_age
    USE mo_fage_process,         ONLY: update_frac_per_age_and_determine_area_to_harvest_from_cohort
    USE mo_carbon_interface,     ONLY: update_delta_carbon_rp_upon_disturbance
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
    REAL(wp), DIMENSION(options%nc)        :: source_fract, stillToHarvest, harvestFromThisCohort

    INTEGER  :: i, oldestCohort, ncohorts, max_age, iblk, ics, ice, nc
    LOGICAL  :: l_sourceisfirstChild
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_harvest_cohort_changes'

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

    dsl4jsb_Config(FAGE_)%l_has_area_changes = .FALSE.
    IF(is_newyear()) THEN
      harvested_area(:) = 0._wp

      IF(ANY(area_to_harvest(:) .GT. 0._wp)) THEN
        dsl4jsb_Config(FAGE_)%l_has_area_changes = .TRUE.

        ! Area changes upon clearcut only have carbon consequences for the first child
        !  (All other children will only have area changes)
        dsl4jsb_Config(FAGE_)%l_only_update_c_of_first_child = .TRUE.

        !NOTE: Here its ASSUMED that the task runs after ageing, therefore youngest cohort is empty 
        !Assertion: TODO - fehler wenn youngest not empty!
        l_sourceisfirstChild = .FALSE.

        firstChild => tile%Get_first_child_tile()
        currentChild => tile%Get_first_child_tile()

        IF (dsl4jsb_Config(HARV_)%harv_kind .EQ. HARVEST_OLDEST_TREES) THEN
          ! as opposed to damaged_frac in update_nat_disturb_cohort_changes area_to_harvest is already with respect
          ! to this tiles area - i.e. not per cohort, thus it does not have to be multiplied with cohort fractions
          stillToHarvest(:) = area_to_harvest(:)

          ! TODO: Identify cohort with oldest available trees for any calculated gridcell?
          oldestCohort = ncohorts

          ! Iterate over the children until currentChild points to that cohort
          DO i=1,oldestCohort-1
            !Assertion that current its not already the last
            IF(currentChild%Is_last_child()) &
              & CALL finish(TRIM(routine), 'Violation of assertion: unexpected last child.')

            currentChild => currentChild%Get_next_sibling_tile()
          ENDDO

          ! then iterate backwards until either all demanded area is fulfilled or the lowest cohort is reached
          DO i=oldestCohort,2,-1
            ! Determine area to be harvested from this cohort
            harvestFromThisCohort = 0._wp
            CALL update_frac_per_age_and_determine_area_to_harvest_from_cohort(dsl4jsb_Config(HARV_)%harv_kind, nc, &
              & max_age, ncohorts, i, cohort_max_age(:,:), fract_per_age(:,:), harvestFromThisCohort(:), stillToHarvest(:))
             
            harvested_area(:) = harvested_area(:) + harvestFromThisCohort(:)

            ! reduce current tiles delta fraction  
            currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - harvestFromThisCohort(:)
            ! and add it to the first childs delta
            firstChild%delta_fract(ics:ice, iblk) = firstChild%delta_fract(ics:ice, iblk) + harvestFromThisCohort(:)

            ! b) carbon consequences 
            CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))
            CALL update_delta_carbon_rp_upon_disturbance( &
              currentChild, firstChild, l_sourceisfirstChild, options, harvestFromThisCohort(:), source_fract(:), CLEARCUT)

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

  END SUBROUTINE update_harvest_cohort_changes

 
  SUBROUTINE aggregate_harvest_cohort_changes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)
    dsl4jsb_Def_memory(HARV_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_harvest_cohort_changes'

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

    IF(dsl4jsb_Config(FAGE_)%l_has_area_changes) THEN 
      dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(HARV_, harvested_area, weighted_by_fract)
    ENDIF
    !No other aggregations required in the calculation of ageing 
    ! (aggregation is done in the CARBON_ process after applying the changes calculated here to carbon variables)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_harvest_cohort_changes

#endif
END MODULE mo_fage_interface
