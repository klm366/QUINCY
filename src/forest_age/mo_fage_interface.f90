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
  USE mo_exception,       ONLY: finish, message, message_text

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  dsl4jsb_Use_processes FAGE_, CARBON_, DISTURB_

  ! Use of process configurations (t_<process>_config)
  dsl4jsb_Use_config(FAGE_)

  ! Use of process memories (t_<process>_memory)
  dsl4jsb_Use_memory(FAGE_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(DISTURB_)
  
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

  ! ================================================================================================================================
  !! Constructors for tasks
  
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

  END SUBROUTINE Register_fage_tasks

  ! ================================================================================================================================
  !>
  !! @brief Task "calculate_ageing" calculates the regquired changes upon ageing  - TODO-JN and more comming
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_calculate_ageing(tile, options)

    USE mo_fage_process,         ONLY: update_frac_per_age_upon_ageing
    USE mo_jsb_time,             ONLY: is_newyear
    USE mo_carbon_interface,     ONLY: update_delta_carbon_of_source_and_sink_tile_upon_area_shift

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age, cohort_max_age

    CLASS(t_jsb_tile_abstract), POINTER    :: currentChild, nextChild
    REAL(wp), DIMENSION(options%nc)        :: fract_this_max_age, source_fract

    REAL(wp) :: fractSum
    INTEGER  :: i,j, ncohorts, iblk, ics, ice, nc
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

    IF(is_newyear()) THEN  ! year has changed - i.e. age - i.e. shift fractions!
      ! Note that this update is done globally at the same time step, 
      !  i.e. for different longitudes the updating happens at different local times.
      dsl4jsb_Config(FAGE_)%l_has_area_changes = .TRUE.

      ! Area changes upon ageing can have carbon consequences for all children
      dsl4jsb_Config(FAGE_)%l_only_update_c_of_first_child = .FALSE.

      ! step one: care for shifts among cohort boundaries
      currentChild => tile%Get_first_child_tile()

      ! only iterate until -1 -- oldest cohort only accumulates in update_calculate_ageing
      DO i=1,ncohorts-1
        !TODO-JN-COHORTS maybe work with a getter?
        nextChild => currentChild%Get_next_sibling_tile()
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))

        ! a) identify delta fraction for shifts from one to the next cohort
        !TODO-JN needs to be adapted when working with absolute areas

        !TODO-JN-COHORTS: cohort_max_age would only be required with this dimension
        !                 if there would be dynamic changes individual per gird cell, 
        !                 if only static (as currently) it could be implemented easier (maybe in the config of fage?!)
        DO j=1,nc
          fract_this_max_age(j) = fract_per_age(j,cohort_max_age(j,i))
        ENDDO

        !TODO-JN-COHORTS maybe generally only work with >= 1.e.10?
        IF(ANY(source_fract(:) - fract_this_max_age(:) .LT. 0._wp)) THEN
          IF (ANY(source_fract(:) - fract_this_max_age(:) + 1.e-10_wp .LT. 0_wp )) THEN
            print *, source_fract(:)
            print *, fract_this_max_age(:)
            CALL finish(TRIM(routine), &
              & 'Violation of assertion: area to be changed needs to be .LE. existing area!')
          ENDIF 
          !TODO-JN-COHORTS two different calculations i.e. numerical errors? -- synchronise fract_per_age
          DO j=1,nc

            IF (fract_this_max_age(j) .GT. source_fract(j)) THEN
              print *, '***************************', i, j
              print *, source_fract(j)
              print *, fract_this_max_age(j)
              !TODO-JN: quick fix, do better!
              fract_this_max_age(j) = source_fract(j)
              fract_per_age(j,cohort_max_age(j,i)) = source_fract(j)
            ENDIF

            IF(i .NE. 1) THEN
              fractSum = SUM(fract_per_age(j,cohort_max_age(j,i-1)+1:cohort_max_age(j,i)))
              IF (fractSum - fract_this_max_age(j) .LT. 0_wp ) THEN
                print *, fractSum
                print *, fract_this_max_age(j)
                CALL finish(TRIM(routine), &
                  & 'Violation of assertion: area to be changed needs to be .LE. existing area!')
              ENDIF
            ENDIF
          ENDDO
        ENDIF

        ! reduce current tiles delta fraction  
        currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - fract_this_max_age(:)

        ! and add them to the next childs delta 
        nextChild%delta_fract(ics:ice, iblk) = nextChild%delta_fract(ics:ice, iblk) + fract_this_max_age(:)

        ! b) accordingly move carbon from one cohort to the next cohorts delta
        CALL update_delta_carbon_of_source_and_sink_tile_upon_area_shift( &
          & currentChild, nextChild, options, fract_this_max_age(:), source_fract(:))

        ! c) TODO-JN-COHORTS what else to shift? 

        ! move further in the child list
        currentChild => nextChild
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
    USE mo_fage_process,         ONLY: update_frac_per_age_equally_upon_disturbance
    USE mo_carbon_interface,     ONLY: update_delta_carbon_upon_disturbance

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

    CLASS(t_jsb_tile_abstract), POINTER    :: firstChild, currentChild, nextChild
    REAL(wp), DIMENSION(options%nc)        :: source_fract, damagedArea, burnedArea, disturbedArea

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

! TODO: wind damage ganz zu ende machen (use statements etc), feuer damage durch auskommentieren ausstellen
! und aus probieren!

!TODO: dann auch feuer -- ah, vielleicht zeitlich erst feuer, weil so in dist? obwohl is arbitrary, oder?!
! ++ auch da über ta, ra und densities denken  
! ++ generalisieren?

        !*** care for wind damage
        damagedArea(:) = damaged_frac(:) * source_fract(:)

        ! wind damage only depends on the forcing and there is a good chance that no area is damaged
        IF(ANY(damagedArea(:) .GT. 0._wp)) THEN

!print *, 'damaged area, i, source_fract', i, damagedArea(:), source_fract(:)

          ! a) changes in area
          ! reduce current tiles delta fraction  
          currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - damagedArea(:)
          ! and add it to the first childs delta
          firstChild%delta_fract(ics:ice, iblk) = firstChild%delta_fract(ics:ice, iblk) + damagedArea(:)

!print *, currentChild%delta_fract(ics:ice, iblk)
!print *, firstChild%delta_fract(ics:ice, iblk)

          ! b) carbon consequences 
          CALL update_delta_carbon_upon_disturbance( &
            currentChild, firstChild, l_sourceisfirstChild, options, damagedArea(:), source_fract(:) )

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


        ! c) TODO-JN-COHORTS what else to shift? 

        ! update the area per age tracking vector regarding this child and the first child
        IF (i .GT. 1) THEN ! area moved within the first cohort requires no update (it stays between year 0 and 1)
          disturbedArea(:) = damagedArea(:) + burnedArea(:)
          IF (dsl4jsb_Config(FAGE_)%disturbance_distribution_scheme == 'equal' ) THEN
            CALL update_frac_per_age_equally_upon_disturbance(nc, ncohorts, max_age, &
              & i, disturbedArea(:), cohort_max_age(:,:), fract_per_age(:,:), mean_age(:))
          ELSE
            CALL finish(TRIM(routine), 'Currently only "static" implemented as forest ageing scheme')
          ENDIF
        ENDIF

        ! Move further in the hierarchy  
        currentChild => currentChild%Get_next_sibling_tile()
        l_sourceisfirstChild = .FALSE.
      ENDDO


! TODO: remove later, only for now to not have any changes for now
!firstChild => tile%Get_first_child_tile()
!currentChild => tile%Get_first_child_tile()
!DO i=1,ncohorts
!currentChild%delta_fract(ics:ice, iblk) = 0._wp
!firstChild%delta_fract(ics:ice, iblk) = 0._wp
!currentChild => currentChild%Get_next_sibling_tile()
!ENDDO

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

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_nat_disturb_cohort_changes


#endif
END MODULE mo_fage_interface
