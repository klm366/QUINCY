!> Contains the interfaces to the climbuf process
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_climbuf_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_jsb_time,        ONLY: is_time_experiment_start
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes A2L_, PHENO_, ASSIMI_, CARBON_, CLIMBUF_, DISTURB_, FAGE_

  ! Use of process configurations
  dsl4jsb_Use_config(CLIMBUF_)
  dsl4jsb_Use_config(FAGE_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(CLIMBUF_)
  dsl4jsb_Use_memory(DISTURB_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_climbuf_tasks !, t_climbuf_process

  CHARACTER(len=*), PARAMETER :: modname = 'mo_climbuf_interface'

  !> Type definition for climbuf task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_climbuf
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_climbuf     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_climbuf  !< Aggregates computed task variables
  END TYPE tsk_climbuf

  !> Constructor interface for climbuf task
  INTERFACE tsk_climbuf
    PROCEDURE Create_task_climbuf                        !< Constructor function for task
  END INTERFACE tsk_climbuf

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for climbuf task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "climbuf"
  !!
  FUNCTION Create_task_climbuf(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_climbuf::return_ptr)
    CALL return_ptr%Construct(name='climbuf', process_id=climbuf_, owner_model_id=model_id)

  END FUNCTION Create_task_climbuf

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for climbuf process
  !!
  !! @param[in,out] this      Instance of climbuf process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_climbuf_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_climbuf(model_id))

  END SUBROUTINE Register_climbuf_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "climbuf"
  !! Task "climbuf" <EXPLANATIION_OF_PROCESS>.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_climbuf(tile, options)

    USE mo_climbuf_process,   ONLY: calc_climbuf
    USE mo_jsb_time,       ONLY: is_newyear, is_newmonth, is_newday, get_day_length, get_previous_month_length
    USE mo_jsb_grid,       ONLY: Get_grid
    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_tile_class, ONLY: t_jsb_tile_abstract

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER  ::  model
    TYPE(t_jsb_grid),  POINTER  ::  grid
    CLASS(t_jsb_tile_abstract),  POINTER  ::  ptr_tile
    LOGICAL  :: lstart, init_running_means, init_mmtemp20
    INTEGER  :: iblk, ics, ice, nc, no_children
    INTEGER  :: i, itile
    REAL(wp) :: seconds_day   ! Number of seconds per day
    REAL(wp) :: seconds_month ! Number of seconds in previous month
    REAL(wp) :: dtime
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_climbuf'
    REAL(wp), POINTER  ::  NPP_mean_5year(:,:)

    ! Declare process memories
    dsl4jsb_Def_config(CLIMBUF_)
    dsl4jsb_Def_config(FAGE_)
    dsl4jsb_Def_memory(CLIMBUF_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(DISTURB_)

    ! Declare pointers to variables in climbuf memory
    dsl4jsb_Real2D_onChunk :: temp_sum_day
    dsl4jsb_Real2D_onChunk :: temp_sum_month
    dsl4jsb_Real2D_onChunk :: min_mmtemp_of_yr
    dsl4jsb_Real2D_onChunk :: max_mmtemp_of_yr
    dsl4jsb_Real2D_onChunk :: min_mmtemp20
    dsl4jsb_Real2D_onChunk :: max_mmtemp20
    dsl4jsb_Real2D_onChunk :: gdd_sum_year
    dsl4jsb_Real2D_onChunk :: gdd_prev_year
    ! Declare pointers to variables in foreign memory
    dsl4jsb_Real2D_onChunk :: t_air

    ! Local variables
    LOGICAL :: new_year, new_month, new_day

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(CLIMBUF_)
    dsl4jsb_Get_config(FAGE_)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(CLIMBUF_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(CLIMBUF_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(DISTURB_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_sum_day)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_sum_month)   ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,min_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,min_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,gdd_sum_year)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,gdd_prev_year)    ! inout
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)             ! in

    IF (is_time_experiment_start(options%current_datetime)) THEN
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF
    new_year = is_newyear(options%current_datetime, options%dtime)
    new_month = is_newmonth(options%current_datetime, options%dtime)
    new_day = is_newday(options%current_datetime, options%dtime)

    seconds_day   = get_day_length()
    seconds_month = get_previous_month_length(options%current_datetime)

    init_running_means = dsl4jsb_Config(CLIMBUF_)%init_running_means
    init_mmtemp20 = .TRUE.
    IF (MAXVAL(min_mmtemp20(:)) > 999._wp .OR. MINVAL(max_mmtemp20(:)) < -999._wp) init_mmtemp20 = .FALSE.
    CALL calc_climbuf(lstart, init_running_means, new_day, new_month, new_year, dtime, t_air, &
                      seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
                      min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20, &
                      gdd_sum_year, gdd_prev_year)

    IF (.NOT. init_mmtemp20 .AND. debug_on() .AND. iblk==1) THEN
       IF (MAXVAL(min_mmtemp20(:)) <= 999._wp .OR. MINVAL(max_mmtemp20(:)) >= -999._wp) &
          CALL message(TRIM(routine), 'initialisation of mmtemp20')
    ENDIF

  END SUBROUTINE update_climbuf

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "climbuf"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_climbuf(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(CLIMBUF_)

    ! Local variables
    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract
    INTEGER  :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_climbuf'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CLIMBUF_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_climbuf

#endif
END MODULE mo_climbuf_interface
