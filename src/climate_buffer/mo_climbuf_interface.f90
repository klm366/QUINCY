!> Contains the interfaces to the climate buffer process
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
!> The CLIMBUF process provides long term mean climatic variables from atmospheric forcing.
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
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes A2L_, CLIMBUF_

  ! Use of process configurations
  dsl4jsb_Use_config(CLIMBUF_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(CLIMBUF_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_climbuf_tasks

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
  !> Constructor for the climate buffer task
  !>
  FUNCTION Create_task_climbuf(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_climbuf::return_ptr)
    CALL return_ptr%Construct(name='climbuf', process_id=climbuf_, owner_model_id=model_id)

  END FUNCTION Create_task_climbuf

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for the climate buffer process
  !>
  SUBROUTINE Register_climbuf_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_climbuf(model_id))

  END SUBROUTINE Register_climbuf_tasks

  ! ======================================================================================================== !
  !>
  !>#### Update the climate buffer
  !>
  !> In this subroutine we update long term mean climatic variables from atmospheric data.
  !>
  SUBROUTINE update_climbuf(tile, options)

    USE mo_climbuf_process,   ONLY: calc_climbuf
    USE mo_jsb_time,          ONLY: is_newyear, is_newmonth, is_newday, get_day_length, get_previous_month_length
    USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Options of the current block

    ! Local variables
    TYPE(t_jsb_model), POINTER  ::  model                 !< Current instance of the ICON-Land

    LOGICAL :: &
      & lstart,              & !< This is the first time step of a simulation
      & new_year,            & !< This is the first time step of a year
      & new_month,           & !< This is the first time step of a month
      & new_day,             & !< This is the first time step of a day
      & init_running_means     !< Initialize (annual) running means, only possible
                               !< after at least one simulated year

    INTEGER :: &
      & iblk,                & !< Currently processed block
      & ics,                 & !< First grid cell of the block
      & ice                    !< Last grid cell of the block

    REAL(wp) :: &
      & seconds_day,         & !< Number of seconds per day
      & seconds_month,       & !< Number of seconds in previous month
      & dtime                  !< Time step lenght [s]

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_climbuf'

    ! Declare process memories
    dsl4jsb_Def_config(CLIMBUF_)   !< Configuration of the climate buffer process
    dsl4jsb_Def_memory(CLIMBUF_)   !< Memory of the climate buffer process
    dsl4jsb_Def_memory(A2L_)       !< Memory of containing atmospheric variables for the land

    ! Declare pointers to variables in climbuf memory
    dsl4jsb_Real2D_onChunk :: temp_sum_day       !< Current temperature sum of the day [K]
    dsl4jsb_Real2D_onChunk :: temp_sum_month     !< Current temperature sum of the month [K]
    dsl4jsb_Real2D_onChunk :: min_mmtemp_of_yr   !< Current minimum monthly mean temperature of the year [K]
    dsl4jsb_Real2D_onChunk :: max_mmtemp_of_yr   !< Current maximum monthly mean temperature of the year [K]
    dsl4jsb_Real2D_onChunk :: min_mmtemp20       !< Climatological minimum monthly mean temperature [K]
    dsl4jsb_Real2D_onChunk :: max_mmtemp20       !< Climatological maximum monthly mean temperature [K]
    dsl4jsb_Real2D_onChunk :: temp_15d           !< 15-day smoothed temperature [K]
    dsl4jsb_Real2D_onChunk :: temp_15d_yDay      !< 15-day smoothed temperature at the end of previous day [K]
    dsl4jsb_Real2D_onChunk :: precip_15d         !< 15-day smoothed precipitation [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: precip_15d_yDay    !< Smoothed precip at the end of previous day [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: gdd_sum_year       !< Current number of growing degree days in the year
    dsl4jsb_Real2D_onChunk :: gdd_prev_year      !< Growing degree days of the previous year
    ! Declare pointers to variables in foreign memory
    dsl4jsb_Real2D_onChunk :: t_air              !< Near surface air temperature [K]
    dsl4jsb_Real2D_onChunk :: rain               !< Liquid precipitation [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: snow               !< Solid precipitation [kg m-2 s-1]

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(CLIMBUF_)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(CLIMBUF_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(CLIMBUF_)
    dsl4jsb_Get_memory(A2L_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   rain)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   snow)              ! in
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_sum_day)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_sum_month)   ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,min_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_mmtemp_of_yr) ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,min_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_mmtemp20)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_15d)         ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_15d_yDay)    ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,precip_15d)       ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,precip_15d_yDay)  ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,gdd_sum_year)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,gdd_prev_year)    ! inout

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
    CALL calc_climbuf(lstart, init_running_means, new_day, new_month, new_year, dtime, t_air, &
      &               rain, snow, seconds_day, seconds_month, temp_sum_day, temp_sum_month,   &
      &               min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp20, max_mmtemp20,         &
      &               temp_15d, temp_15d_yDay, precip_15d, precip_15d_yDay,                   &
      &               gdd_sum_year, gdd_prev_year)

  END SUBROUTINE update_climbuf

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of climate buffer variables
  !>
  !> The climate buffer variables are aggregated in this routine. (Currently not needed)
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE aggregate_climbuf(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    dsl4jsb_Def_memory(CLIMBUF_)                           !< Memory of the climate buffer process

    ! Local variables
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregation method: area weighted fractions

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_climbuf'
    ! ------------------------------------------------------------------------------------------------------ !

    IF (debug_on() .AND. options%iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CLIMBUF_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    IF (debug_on() .AND. options%iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_climbuf

#endif
END MODULE mo_climbuf_interface
