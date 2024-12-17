!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!!
!! @brief diverse subroutines and functions for the plant-model
!!
!! @author
!>  Albert Jornet Puig, University of Leipzig / MPI-BGC Jena
!>  Ana Bastos, University of Leipzig / MPI-BGC Jena
!!
!! @par Revision History
!!
MODULE mo_dfire_util
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: eps8

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_time_average_dist_fire, calculate_daily_weather, calculate_pure_moving_averages

  CHARACTER(len=*), PARAMETER :: modname = 'mo_dfire_util'

CONTAINS

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for DFIRE_
  !>
  SUBROUTINE calculate_time_average_dist_fire(tile, options)

    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    !USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_lnd_time_averages,       ONLY: mavg_period_dfire, calc_time_mavg
    !------------------------------------------------------------------------------------------------------ !
    dsl4jsb_Use_processes A2L_, DFIRE_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER         :: model                !< instance of the model
    !TYPE(t_lctlib_element), POINTER         :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_canopy         !< Vertical grid
    INTEGER                                 :: ncanopy              !< number of canopy layers, from vgrid
    INTEGER                                 :: iblk, ics, ice, nc   !< dimensions
    REAL(wp)                                :: dtime                !< timestep
    CHARACTER(len=*), PARAMETER             :: routine = modname//':calculate_time_average_dist_fire'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: rain
    dsl4jsb_Real2D_onChunk      :: snow
    ! DFIRE_
    dsl4jsb_Real2D_onChunk      :: precip_qc_mavg
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    write(*,*) "dfire_util.f90:: calculate_time_avg_dist_fire, maybe?"

    IF (.NOT. tile%Is_process_calculated(DFIRE_)) RETURN
    !IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    !lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy  => Get_vgrid('q_canopy_layer')
    ncanopy       =  vgrid_canopy%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(DFIRE_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    ! DFIRE_
    dsl4jsb_Get_var2D_onChunk(DFIRE_, precip_qc_mavg)
    ! ----------------------------------------------------------------------------------------------------- !
    write(*,*) "dfire_util.f90:: calculate_time_avg_dist_fire..."

    precip_qc_mavg(:) = calc_time_mavg(dtime, precip_qc_mavg(:), (rain(:)+snow(:)), mavg_period_dfire)

  END SUBROUTINE calculate_time_average_dist_fire



  ! ======================================================================================================= !
  !>
  !> Subtask to estimate daily weather variables: precipitation, t_air_max and t_air_min
  !>
  !> A day is defined from the lowest point of the sun of the previous day to the current one.
  !> The lowest point of the sun of the current day is estimated from the zenith.
  !>
  SUBROUTINE calculate_daily_weather(tile, options)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    !USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element

    ! Use of process memories
    dsl4jsb_Use_processes A2L_, DFIRE_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(DFIRE_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile               !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options            !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER         :: model              !< instance of the model
    !TYPE(t_lctlib_element), POINTER         :: lctlib             !< land-cover-type library - parameter across pft's
    INTEGER                                 :: iblk, ics, ice, nc
    REAL(wp)                                :: dtime              !< timestep
    LOGICAL, ALLOCATABLE                    :: is_midnight(:)     !< is peak of the day
    LOGICAL, ALLOCATABLE                    :: is_sun_upwards_in_curr_ts(:)
    LOGICAL, ALLOCATABLE                    :: is_sun_downwards_in_prev_ts(:)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calculate_daily_weather'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(DFIRE_)
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: rain
    dsl4jsb_Real2D_onChunk      :: snow
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: cos_zenith_angle
    ! DFIRE_ 2D
    dsl4jsb_Real2D_onChunk      :: precip_daytime         !<
    dsl4jsb_Real2D_onChunk      :: t_air_max              !<
    dsl4jsb_Real2D_onChunk      :: t_air_min              !<
    REAL(wp), ALLOCATABLE, SAVE :: prev_one_cos_zenith_angle(:), prev_two_cos_zenith_angle(:)
    REAL(wp), ALLOCATABLE, SAVE :: precip_daytime_old(:)
    REAL(wp), ALLOCATABLE, SAVE :: t_air_max_old(:)
    REAL(wp), ALLOCATABLE, SAVE :: t_air_min_old(:)
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ---------------------------
    write(*,*) "dfire_util.f90:: calculate_daily_weather, maybe?"
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(DFIRE_)) RETURN   ! ask for all processes considered here ?
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    !lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    ! Get process config
    !..
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(DFIRE_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)
    dsl4jsb_Get_var2D_onChunk(A2L_, cos_zenith_angle)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(DFIRE_, precip_daytime)
    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_max)
    dsl4jsb_Get_var2D_onChunk(DFIRE_, t_air_min)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    write(*,*) "dfire_util.f90:: calculate_daily_weather..., t_air_min(1), t_air_max =", &
            t_air_min(1), t_air_max(1)

    !> 1.0 init
    !!
    IF (.NOT. ALLOCATED(precip_daytime_old)) THEN
      ALLOCATE(precip_daytime_old(nc))
      precip_daytime_old = 0._wp
    ENDIF
    IF (.NOT. ALLOCATED(t_air_max_old)) THEN
      ALLOCATE(t_air_max_old(nc))
      t_air_max_old = 0._wp
    ENDIF
    IF (.NOT. ALLOCATED(t_air_min_old)) THEN
      ALLOCATE(t_air_min_old(nc))
      t_air_min_old = 0._wp
    ENDIF
    IF (.NOT. ALLOCATED(prev_one_cos_zenith_angle)) THEN
      ALLOCATE(prev_one_cos_zenith_angle(nc))
      prev_one_cos_zenith_angle = 0._wp
    ENDIF
    IF (.NOT. ALLOCATED(prev_two_cos_zenith_angle)) THEN
      ALLOCATE(prev_two_cos_zenith_angle(nc))
      prev_two_cos_zenith_angle = 0._wp
    ENDIF

    ALLOCATE(is_midnight(nc))
    ALLOCATE(is_sun_upwards_in_curr_ts(nc))
    ALLOCATE(is_sun_downwards_in_prev_ts(nc))
    is_midnight(:) = .FALSE.
    is_sun_upwards_in_curr_ts(:) = .FALSE.
    is_sun_downwards_in_prev_ts(:) = .FALSE.

    ! 1.1 update min/max values
    t_air_max_old(:) = max(t_air_max_old(:), t_air(:))
    t_air_min_old(:) = min(t_air_min_old(:), t_air(:))
    ! 1.2 update accumulation values
    precip_daytime_old(:) = precip_daytime_old(:) + ((rain(:) + snow(:)) / (one_day / dtime))

    !> 2.0 at the timestep after local midnight (1st timestep of the new day) calc daytime averages
    !!
    !! calculate the min/max of the previous day and pass this value to the according variable \n
    !! the timestep right after midnight is identified based on changes in cos_zenith_angle (turning from negative to positive)

    !> 2.1 1D variables
    !!
    is_sun_upwards_in_curr_ts(:) = ( (cos_zenith_angle(:) - prev_one_cos_zenith_angle(:)) > eps8)
    is_sun_downwards_in_prev_ts(:) = ( (prev_one_cos_zenith_angle(:) - prev_two_cos_zenith_angle(:)) <= eps8 )
    is_midnight(:) = (is_sun_upwards_in_curr_ts(:) .AND. is_sun_downwards_in_prev_ts(:))

    WHERE (is_midnight(:))
      ! check if at least one value had been assigned at the current day, avoiding division by zero
      ! calculate the average over the previous day and pass it to the according *_daytime variable;
      t_air_max(:)            = t_air_max_old(:)
      t_air_min(:)            = t_air_min_old(:)
      precip_daytime(:)       = precip_daytime_old(:)
      ! init variables after min/max daily value
      t_air_max_old(:)        = 0.0_wp     !< zero kelvin
      t_air_min_old(:)        = 9999.0_wp  !< above max
      precip_daytime_old(:)   = 0._wp      !< no precipitation
    END WHERE

    !    IF (debug_on() .AND. iblk == 1) THEN
    IF (is_midnight(1)) THEN
      ! check if at least one value had been assigned at the current day, avoiding division by zero
      ! calculate the average over the previous day and pass it to the according *_daytime variable;
      IF (t_air_max(1) < t_air_min(1))  THEN
        write(*,*) "error: calc_daily_weather:: t_air_min, t_air_max =", t_air_min, t_air_max
        CALL finish("calc_daily_weather", "t_air_min cant be higher than t_air_max")
      ENDIF
    ENDIF
    write(*,*) "dfire_util.f90:: (end) calculate_daily_weather, t_air_min(1), t_air_max(1) =", t_air_min(1), t_air_max(1)

    !> 3.0 update cos_zenith_angle of previous time steps
    !!
    prev_two_cos_zenith_angle(:) = prev_one_cos_zenith_angle(:)
    prev_one_cos_zenith_angle(:) = cos_zenith_angle(:)
  END SUBROUTINE calculate_daily_weather


  ! ======================================================================================================= !
  !>
  !> Subtask to calculate pure moving average
  !>
  !>
  SUBROUTINE calculate_pure_moving_averages(tile, options)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    !USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_quincy_mavg_class,       ONLY: t_mavg

    ! Use of process memories
    dsl4jsb_Use_processes A2L_, DFIRE_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(DFIRE_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile               !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options            !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER         :: model                !< instance of the model
    !TYPE(t_lctlib_element), POINTER         :: lctlib               !< land-cover-type library - parameter across pft's
    INTEGER                                 :: iblk, ics, ice, nc
    REAL(wp)                                :: dtime                !< timestep
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calculate_pure_moving_averages'
    TYPE(t_mavg), SAVE, ALLOCATABLE :: precip_mavg_type(:)
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(DFIRE_)
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: rain
    dsl4jsb_Real2D_onChunk      :: snow
    ! DFIRE_ 2D
    dsl4jsb_Real2D_onChunk      :: precip_mavg
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    write(*,*) "dfire_util.f90:: calculate_pure_moving_averages, maybe?"
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(DFIRE_)) RETURN   ! ask for all processes considered here ?
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    !lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    ! Get process config
    !..
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(DFIRE_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_, rain)
    dsl4jsb_Get_var2D_onChunk(A2L_, snow)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(DFIRE_, precip_mavg)

    write(*,*) "dfire_util.f90:: calculate_pure_moving_averages..."

    IF (.NOT. ALLOCATED(precip_mavg_type)) THEN
      ALLOCATE(precip_mavg_type(nc))
      precip_mavg_type = t_mavg(dtime, 1) !< daily
    ENDIF

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    DO ic = 1, nc
      CALL precip_mavg_type(ic)%add((snow(ic) + rain(ic)))
      precip_mavg(ic) = precip_mavg_type(ic)%get()
    END DO

  END SUBROUTINE calculate_pure_moving_averages

#endif
END MODULE mo_dfire_util
