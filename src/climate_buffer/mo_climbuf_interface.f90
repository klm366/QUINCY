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
!>### The CLIMBUF process provides long term mean climatic variables from atmospheric forcing.
!>
!> This module contains the climate buffer interface routines. Core routine is [[update_climbuf]], where
!> the climate variables are calculated from atmospheric forcing variables.
!>
MODULE mo_climbuf_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_kind,               ONLY: wp
  USE mo_jsb_control,        ONLY: debug_on, acc_stream
  USE mo_jsb_time,           ONLY: is_time_experiment_start, is_time_restart, &
    &                              is_newyear, is_newmonth, is_newday, &
    &                              get_day_length, get_previous_month_length
  USE mo_exception,          ONLY: message
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes A2L_, SEB_, SSE_, HYDRO_, CLIMBUF_

  ! Use of process configurations
  dsl4jsb_Use_config(CLIMBUF_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(HYDRO_)
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
    PROCEDURE, NOPASS :: Integrate => update_climbuf     !< Advances task computation for one time step
    PROCEDURE, NOPASS :: Aggregate => aggregate_climbuf  !< Aggregates computed task variables
  END TYPE tsk_climbuf

  !> Constructor interface for climbuf task
  INTERFACE tsk_climbuf
    PROCEDURE Create_task_climbuf                        !< Constructor function for task
  END INTERFACE tsk_climbuf

CONTAINS

  ! ======================================================================================================
  !! Constructors for tasks

  ! ------------------------------------------------------------------------------------------------------
  !> Constructor for the climate buffer task
  !>
  FUNCTION Create_task_climbuf(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id    !< ID of current instance of the model
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "climbuf"

    ALLOCATE(tsk_climbuf::return_ptr)
    CALL return_ptr%Construct(name='climbuf', process_id=climbuf_, owner_model_id=model_id)

  END FUNCTION Create_task_climbuf

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for the climate buffer process
  !>
  SUBROUTINE Register_climbuf_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this          !< Instance of climbuf process task
    INTEGER,              INTENT(in)    :: model_id      !< Model ID

    CALL this%Register_task(tsk_climbuf(model_id))

  END SUBROUTINE Register_climbuf_tasks

  ! ======================================================================================================== !
  !>
  !>#### Update the climate buffer
  !>
  !> In this subroutine we update long term mean climatic variables from atmospheric forcing data.
  !>
  SUBROUTINE update_climbuf(tile, options)

    USE mo_climbuf_process,        ONLY: calc_climbuf_2d, calc_climbuf_3d, get_relative_humidity_air
    USE mo_jsb_tile_class,         ONLY: t_jsb_tile_abstract
    USE mo_jsb_physical_constants, ONLY: rhoh2o

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
      & init_climbuf,        & !< Initialize climate buffer (needed if climbuf process was not active in the
                               !< simulation from which the restart file originates)
      & init_running_means     !< Initialize (annual) running means, only possible after at least one
                               !< full year of simulation after climate buffer initialization

    INTEGER :: &
      & iblk,                & !< Currently processed block
      & ics,                 & !< First grid cell of the block
      & ice,                 & !< Last grid cell of the block
      & ic,                  & !< Looping index for grid cells
      & nc                     !< Number of grid cells in block

    REAL(wp) :: &
      & seconds_day,         & !< Number of seconds per day
      & seconds_month,       & !< Number of seconds in previous month
      & dtime                  !< Time step length [s]

    REAL(wp), ALLOCATABLE :: &
      & snow_depth(:),       & !< Snow depth on soil [m]
      & q_rel_air(:)           !< Relative humidity of the lowest level of the atmosphere [%]

    LOGICAL :: &
      & use_tmx                !< Use tmx scheme; needed as local variable on GPUs

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_climbuf'

    ! Declare process memories
    dsl4jsb_Def_config(CLIMBUF_)   !< Configuration of the climate buffer process
    dsl4jsb_Def_memory(CLIMBUF_)   !< Memory of the climate buffer process
    dsl4jsb_Def_memory(A2L_)       !< Memory for atmospheric variables for the land
    dsl4jsb_Def_memory(SEB_)       !< Memory of surface energy balance process
    dsl4jsb_Def_memory(SSE_)       !< Memory of soil and snow energy process
    dsl4jsb_Def_memory(HYDRO_)     !< Memory of hydrology process

    ! Declare pointers to variables in climbuf memory
    dsl4jsb_Real2D_onChunk :: temp_sum_day       !< Current temperature sum of the day [K]
    dsl4jsb_Real2D_onChunk :: temp_sum_month     !< Current temperature sum of the month [K]
    dsl4jsb_Real2D_onChunk :: min_mmtemp_of_yr   !< Current minimum monthly mean temperature of the year [K]
    dsl4jsb_Real2D_onChunk :: max_mmtemp_of_yr   !< Current maximum monthly mean temperature of the year [K]
    dsl4jsb_Real2D_onChunk :: min_mmtemp_20yr    !< Climatological minimum monthly mean temperature [K]
    dsl4jsb_Real2D_onChunk :: max_mmtemp_20yr    !< Climatological maximum monthly mean temperature [K]
    dsl4jsb_Real2D_onChunk :: temp_15d           !< 15-day smoothed temperature [K]
    dsl4jsb_Real2D_onChunk :: temp_15d_yDay      !< 15-day smoothed temperature at the end of previous day [K]
    dsl4jsb_Real2D_onChunk :: precip_15d         !< 15-day smoothed precipitation [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: precip_15d_yDay    !< 15-day smoothed precipitation at the end of previous day [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: gdd_sum_year       !< Current number of growing degree days in the year
    dsl4jsb_Real2D_onChunk :: gdd_prev_year      !< Growing degree days of the previous year
    dsl4jsb_Real2D_onChunk :: max_snow_year      !< Maximum snow depth in current year [m]
    dsl4jsb_Real2D_onChunk :: max_snow_prev_year !< Maximum snow depth in previous year [m]
    dsl4jsb_Real2D_onChunk :: max_thaw_year      !< Maximum soil column thaw depth in current year [m]
    dsl4jsb_Real2D_onChunk :: max_thaw_prev_year !< Maximum soil column thaw depth in previous year [m]
    dsl4jsb_Real2D_onChunk :: max_wind_of_day    !< Current maximum 10m wind speed of the day
    dsl4jsb_Real2D_onChunk :: max_wind_yDay      !< Maximum 10m wind speed of the previous day
    dsl4jsb_Real2D_onChunk :: max_wind_4yr       !< Maximum 10m wind speed smoothed in time
    dsl4jsb_Real2D_onChunk :: q_rel_air_14d      !< Relative humidity smoothed in time (~14 days) [%]
    dsl4jsb_Real2D_onChunk :: q_rel_air_14d_yDay !< Previous day value of smoothed relative humidity [%]
    dsl4jsb_Real3D_onChunk :: t_soil_sl_15d      !< 15-day smoothed soil layer temperature [K]
    dsl4jsb_Real3D_onChunk :: t_soil_sl_15d_yDay !< 15-day smoothed soil layer temperature at the end of the previous day [K]
    dsl4jsb_Real3D_onChunk :: wtr_soil_rel_sl_15d      !< 15-day smoothed relative soil layer moisture [-]
    dsl4jsb_Real3D_onChunk :: wtr_soil_rel_sl_15d_yDay !< 15-day smoothed relative soil layer moisture at end of previous day [-]
    dsl4jsb_Real3D_onChunk :: max_t_soil_sl_year !< Maximum soil layer temperature in previous year [K]
    dsl4jsb_Real3D_onChunk :: permafrost_period_sl     !< Number of years of consecutively frozen soil layer [year]

    ! Declare pointers to variables in foreign memory
    dsl4jsb_Real2D_onChunk :: t_air              !< Near surface air temperature [K]
    dsl4jsb_Real2D_onChunk :: rain               !< Liquid precipitation [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: snow               !< Solid precipitation [kg m-2 s-1]
    dsl4jsb_Real2D_onChunk :: wind_10m           !< 10m wind speed [m/s]
    dsl4jsb_Real2D_onChunk :: t_unfilt           !< Surface temperature [K]
    dsl4jsb_Real2D_onChunk :: q_air              !< Specific humidity of lowest level of the atmosphere []
    dsl4jsb_Real2D_onChunk :: press_srf          !< Surface pressure [N m-2]
    dsl4jsb_Real2D_onChunk :: weq_snow_soil      !< Snow water equivalent on the soils surface [m]
    dsl4jsb_Real2D_onChunk :: snow_soil_dens     !< Density of snow on soil [kg m-3]
    dsl4jsb_Real2D_onChunk :: thaw_depth         !< Active layer (thaw) depth in soil column [m]
    dsl4jsb_Real3D_onChunk :: t_soil_sl          !< Soil layer temperature [K]
    dsl4jsb_Real3D_onChunk :: wtr_soil_rel_sl    !< Relative soil layer water content [m]

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(CLIMBUF_)

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(CLIMBUF_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(CLIMBUF_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   rain)                    ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   snow)                    ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   wind_10m)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   q_air)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   press_srf)               ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)                ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)           ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)          ! in
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth)              ! in
    dsl4jsb_Get_var3D_onChunk(SSE_,   t_soil_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_rel_sl)         ! in

    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_sum_day)           ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_sum_month)         ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,min_mmtemp_of_yr)       ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_mmtemp_of_yr)       ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,min_mmtemp_20yr)        ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_mmtemp_20yr)        ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_15d)               ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,temp_15d_yDay)          ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,precip_15d)             ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,precip_15d_yDay)        ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,gdd_sum_year)           ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,gdd_prev_year)          ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_snow_year)          ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_snow_prev_year)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_thaw_year)          ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_thaw_prev_year)     ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_wind_of_day)        ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_wind_yDay)          ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,max_wind_4yr)           ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,q_rel_air_14d)          ! inout
    dsl4jsb_Get_var2D_onChunk(CLIMBUF_,q_rel_air_14d_yDay)     ! inout
    dsl4jsb_Get_var3D_onChunk(CLIMBUF_,t_soil_sl_15d)          ! inout
    dsl4jsb_Get_var3D_onChunk(CLIMBUF_,t_soil_sl_15d_yDay)     ! inout
    dsl4jsb_Get_var3D_onChunk(CLIMBUF_,wtr_soil_rel_sl_15d)    ! inout
    dsl4jsb_Get_var3D_onChunk(CLIMBUF_,wtr_soil_rel_sl_15d_yDay)!inout
    dsl4jsb_Get_var3D_onChunk(CLIMBUF_,max_t_soil_sl_year)     ! inout
    dsl4jsb_Get_var3D_onChunk(CLIMBUF_,permafrost_period_sl)   !inout

    ! Generally, the climate buffer is initialized at model initialization (is_time_experiment_start).
    ! There is however the possibility to initialize it in restarted simulations by namelist parameter
    ! 'init_climbuf'. We here make sure, the buffer is initialized only at the first time step of the run.
    IF ((is_time_experiment_start(options%current_datetime) .OR. is_time_restart(options%current_datetime)) &
        & .AND. dsl4jsb_Config(CLIMBUF_)%init_climbuf) THEN
      init_climbuf = .TRUE.
    ELSE
      init_climbuf = .FALSE.
    END IF
    IF (is_time_experiment_start(options%current_datetime) .OR. init_climbuf) THEN
      lstart = .TRUE.
    ELSE
      lstart = .FALSE.
    END IF

    new_year = is_newyear(options%current_datetime, options%dtime)
    new_month = is_newmonth(options%current_datetime, options%dtime)
    new_day = is_newday(options%current_datetime, options%dtime)

    seconds_day   = get_day_length()
    seconds_month = get_previous_month_length(options%current_datetime)

    ALLOCATE (snow_depth(nc), q_rel_air(nc))
    !$ACC DATA ASYNC(acc_stream) CREATE(snow_depth, q_rel_air)

    use_tmx=model%config%use_tmx
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(acc_stream)
    DO ic = 1, nc
      q_rel_air(ic)  = get_relative_humidity_air(q_air(ic), t_unfilt(ic), press_srf(ic), use_tmx)
      snow_depth(ic) = weq_snow_soil(ic) * rhoh2o / snow_soil_dens(ic)
    END DO
    !$ACC END PARALLEL LOOP

    init_running_means = dsl4jsb_Config(CLIMBUF_)%init_running_means
    CALL calc_climbuf_2d(lstart, init_running_means, new_day, new_month, new_year, dtime, t_air, &
      &                  rain, snow, wind_10m, q_rel_air, snow_depth, thaw_depth, &
      &                  seconds_day, seconds_month, temp_sum_day, temp_sum_month, &
      &                  min_mmtemp_of_yr, max_mmtemp_of_yr, min_mmtemp_20yr, max_mmtemp_20yr, &
      &                  temp_15d, temp_15d_yDay, precip_15d, precip_15d_yDay, &
      &                  gdd_sum_year, gdd_prev_year, &
      &                  max_snow_year, max_snow_prev_year, max_thaw_year, max_thaw_prev_year, &
      &                  max_wind_of_day, max_wind_yDay, max_wind_4yr, &
      &                  q_rel_air_14d, q_rel_air_14d_yDay)
    CALL calc_climbuf_3d(lstart, new_day, new_year, dtime, t_soil_sl, wtr_soil_rel_sl, seconds_day, &
      &                  t_soil_sl_15d, t_soil_sl_15d_yDay, wtr_soil_rel_sl_15d, wtr_soil_rel_sl_15d_yDay, &
      &                  max_t_soil_sl_year, permafrost_period_sl)

    !$ACC END DATA
    DEALLOCATE (snow_depth, q_rel_air)

  END SUBROUTINE update_climbuf

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of climate buffer variables
  !>
  !> The climate buffer variables are aggregated in this routine.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE aggregate_climbuf(tile, options)

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile      !< Current tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options   !< Options of the current block

    dsl4jsb_Def_memory(CLIMBUF_)                           !< Memory of the climate buffer process

    ! Local variables
    INTEGER :: iblk, ics, ice

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract  !< Aggregation method: area weighted fractions

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_climbuf'
    ! ------------------------------------------------------------------------------------------------------ !

    ! Get local variables from options argument (needed in dsl4jsb_Aggregate_onChunk)
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. options%iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(CLIMBUF_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(CLIMBUF_, q_rel_air_14d, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, q_rel_air_14d_yday, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_wind_yDay, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_wind_4yr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_wind_of_day, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, temp_sum_day, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, temp_sum_month, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, min_mmtemp_of_yr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_mmtemp_of_yr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, min_mmtemp_20yr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_mmtemp_20yr, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, temp_15d, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, temp_15d_yDay, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, precip_15d, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, precip_15d_yDay, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, gdd_sum_year, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, gdd_prev_year, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_snow_year, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_snow_prev_year, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_thaw_year, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_thaw_prev_year, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, t_soil_sl_15d, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, t_soil_sl_15d_yDay, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, wtr_soil_rel_sl_15d, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, wtr_soil_rel_sl_15d_yDay, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, max_t_soil_sl_year, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(CLIMBUF_, permafrost_period_sl, weighted_by_fract)

    IF (debug_on() .AND. options%iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_climbuf

#endif
END MODULE mo_climbuf_interface
