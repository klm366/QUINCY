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
  dsl4jsb_Use_processes FAGE_

  ! Use of process configurations (t_<process>_config)
  dsl4jsb_Use_config(FAGE_)

  ! Use of process memories (t_<process>_memory)
  dsl4jsb_Use_memory(FAGE_)  
  
  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_fage_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_interface'

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_ageing
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_fage    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_fage !< Aggregates computed task variables
  END TYPE   tsk_ageing

  INTERFACE   tsk_ageing
    PROCEDURE Create_task_ageing 
  END INTERFACE   tsk_ageing

CONTAINS
  ! ================================================================================================================================
  !! Constructors for tasks
  
  !> Constructor for ageing task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "ageing"
  !!
  FUNCTION Create_task_ageing(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_ageing::return_ptr)

    return_ptr%name = 'ageing'
    return_ptr%process_id = FAGE_
    return_ptr%owner_model_id = model_id

  END FUNCTION Create_task_ageing

  ! =======================================================================================================
  !> Register tasks for fage process
  !!
  !! @param[in,out] this      Instance of fage process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_fage_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_ageing(model_id))

  END SUBROUTINE Register_fage_tasks

  ! ================================================================================================================================
  !>
  !! @brief Implementation to update the fage state
  !! Task "ageing" calculates the age  - TODO-JN and more comming
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  ! 
  SUBROUTINE update_fage(tile, options)

    USE mo_fage_process,         ONLY: update_frac_per_age_tracking
    USE mo_jsb_time,             ONLY: is_newyear
    USE mo_carbon_interface,     ONLY: &
                                       & update_delta_carbon_of_source_and_sink_tile_upon_area_shift,            &
                                       & update_carbon_states_with_delta_carbon

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(FAGE_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: mean_age
    dsl4jsb_Real3D_onChunk :: fract_per_age, cohort_max_age

    CLASS(t_jsb_tile_abstract), POINTER    :: currentChild, nextChild
    REAL(wp), DIMENSION(options%nc)        :: fract_this_max_age, source_fract

    INTEGER  :: i,j, ncohorts, iblk, ics, ice, nc
    TYPE(t_jsb_model), POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_fage'

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

      ! step one: care for shifts among cohort boundaries
      currentChild => tile%Get_first_child_tile()

      ! only iterate until -1 -- oldest cohort only accumulates in update_fage
      DO i=1,ncohorts-1
        !TODO-JN-COHORTS maybe work with a getter?
        nextChild => currentChild%Get_next_sibling_tile()

        ! a) identify delta fraction for shifts from one to the next cohort
        !TODO-JN needs to be adapted when working with absolute areas

        !TODO-JN-COHORTS: cohort_max_age would only be required with this dimension
        !                 if there would be dynamic changes individual per gird cell, 
        !                 if only static (as currently) it could be implemented easier (maybe in the config of fage?!)
        DO j=1,nc
          fract_this_max_age(j) = fract_per_age(j,cohort_max_age(j,i))
        ENDDO

        ! reduce current tiles delta fraction  
        currentChild%delta_fract(ics:ice, iblk) = currentChild%delta_fract(ics:ice, iblk) - fract_this_max_age(:)

        ! and add them to the next childs delta 
        nextChild%delta_fract(ics:ice, iblk) = nextChild%delta_fract(ics:ice, iblk) + fract_this_max_age(:)

        ! b) accordingly move carbon from one cohort to the next cohorts delta
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))

        CALL update_delta_carbon_of_source_and_sink_tile_upon_area_shift( &
          & currentChild, nextChild, options, fract_this_max_age(:), source_fract(:))

        ! c) TODO-JN-COHORTS what else to shift? 

        ! move further in the child list
        currentChild => nextChild
      ENDDO

      ! step two: apply changes 
      currentChild => tile%Get_first_child_tile()
      DO i=1,ncohorts

        ! a) TODO-JN-COHORTS what else to update? 

        ! b) update C states with deltas (+ reset C deltas to 0)
        CALL update_carbon_states_with_delta_carbon(currentChild, options)

        ! c) update area fractions 
        ! TODO-JN: maybe make it a function in mo_jsb_tile_class: update_fraction?
        CALL currentChild%Get_fraction(ics, ice, iblk, fract=source_fract(:))
        CALL currentChild%Set_fraction(ics, ice, iblk, &
          fract=source_fract(:)+currentChild%delta_fract(ics:ice, iblk))
        ! + reset deltas to 0
        currentChild%delta_fract(ics:ice, iblk) = 0._wp
        ! move further in the child list
        currentChild => currentChild%Get_next_sibling_tile()
      ENDDO

      ! update the area per age tracking vector
      CALL update_frac_per_age_tracking(nc, dsl4jsb_Config(FAGE_)%max_age, fract_per_age(:,:), mean_age(:))

    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_fage

  ! TODO-JN-COHORTS: does not really make sense to aggregate mean_age?! 
  SUBROUTINE aggregate_fage(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_memory(FAGE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':aggregate_fpc'

    INTEGER  :: iblk, ics, ice, nc

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(FAGE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(FAGE_, mean_age, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_fage

#endif
END MODULE mo_fage_interface
