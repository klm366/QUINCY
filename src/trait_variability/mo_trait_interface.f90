! Author: Dr. Sri Nandini-Weiß
! Date: 2026-07-15
!> Interface to the trait-variability process
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
!>#### Calculates annual diagnostic PFT13 traits on the PFT13 tile
!>
MODULE mo_trait_interface
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_exception,         ONLY: message, finish
  USE mo_jsb_time,          ONLY: is_time_experiment_start, is_newyear
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class, ONLY: t_jsb_process
  USE mo_jsb_task_class,    ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes CLIMBUF_, TRAIT_
  dsl4jsb_Use_memory(CLIMBUF_)
  dsl4jsb_Use_memory(TRAIT_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_trait_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_trait_interface'

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_trait
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_trait
    PROCEDURE, NOPASS :: Aggregate => aggregate_trait
  END TYPE tsk_trait

  INTERFACE tsk_trait
    PROCEDURE Create_task_trait
  END INTERFACE tsk_trait

CONTAINS

  FUNCTION Create_task_trait(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_trait::return_ptr)
    CALL return_ptr%Construct(name='trait', process_id=TRAIT_, owner_model_id=model_id)

  END FUNCTION Create_task_trait

  SUBROUTINE Register_trait_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_trait(model_id))

  END SUBROUTINE Register_trait_tasks

  SUBROUTINE update_trait(tile, options)

    ! Trait generalisation: import the per-type model mask.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-09-23
    USE mo_trait_process, ONLY: calc_trait_diagnostics, max_lct_id, &
      &                      has_dynamic_vcmax25, has_dynamic_jmax25, has_dynamic_sla

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk, ics, ice, nc
    INTEGER :: lib_id
    LOGICAL :: new_year, lstart
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_trait'

    CLASS(t_jsb_tile_abstract), POINTER :: climbuf
    dsl4jsb_Def_memory_tile(CLIMBUF_, climbuf)
    dsl4jsb_Def_memory(TRAIT_)

    dsl4jsb_Real2D_onChunk :: mat_prev_year_climbuf
    dsl4jsb_Real2D_onChunk :: map_prev_year_climbuf
    dsl4jsb_Real2D_onChunk :: max_mmtemp_prev_year_climbuf
    ! The four added previous-year predictors (Author: Dr. Sri Nandini-Weiß, 2026-10-01)
    dsl4jsb_Real2D_onChunk :: min_mmtemp_prev_year_climbuf
    dsl4jsb_Real2D_onChunk :: reh_prev_year_climbuf
    dsl4jsb_Real2D_onChunk :: soil_moist_prev_year_climbuf
    dsl4jsb_Real2D_onChunk :: nswr_prev_year_climbuf

    dsl4jsb_Real2D_onChunk :: vcmax25_raw
    dsl4jsb_Real2D_onChunk :: vcmax25_pred
    dsl4jsb_Real2D_onChunk :: jmax25_raw
    dsl4jsb_Real2D_onChunk :: jmax25_pred
    dsl4jsb_Real2D_onChunk :: sla_raw
    dsl4jsb_Real2D_onChunk :: sla_pred
    dsl4jsb_Real2D_onChunk :: jmax_vcmax_ratio

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(TRAIT_)) RETURN

    ! Generalisation: the process is no longer restricted to the PFT13 tile.
    !   Selection is by lctlib id and by the per-type model mask.
    ! Author: Dr. Sri Nandini-Weiß
    ! Date: 2026-09-23
    ! The trait equations are selected by the lctlib id of the tile's primary land
    ! cover type. lcts(1) is always the primary lct; an id outside the coefficient
    ! table, or a tile carrying no lct, is not a plant functional type tile.
    IF (.NOT. ALLOCATED(tile%lcts)) RETURN
    lib_id = tile%lcts(1)%lib_id
    IF (lib_id < 1 .OR. lib_id > max_lct_id) RETURN

    ! Plant functional types without a fitted model keep their static lctlib
    ! parameters and are not evaluated.
    IF (.NOT. (has_dynamic_vcmax25(lib_id) .OR. has_dynamic_jmax25(lib_id) .OR. &
      &        has_dynamic_sla(lib_id))) RETURN

    new_year = is_newyear(options%current_datetime, options%dtime)
    lstart   = is_time_experiment_start(options%current_datetime)
    IF (.NOT. new_year .OR. lstart) RETURN

    IF (debug_on() .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(TRAIT_)
    dsl4jsb_Get_var2D_onChunk(TRAIT_, vcmax25_raw)
    dsl4jsb_Get_var2D_onChunk(TRAIT_, vcmax25_pred)
    dsl4jsb_Get_var2D_onChunk(TRAIT_, jmax25_raw)
    dsl4jsb_Get_var2D_onChunk(TRAIT_, jmax25_pred)
    dsl4jsb_Get_var2D_onChunk(TRAIT_, sla_raw)
    dsl4jsb_Get_var2D_onChunk(TRAIT_, sla_pred)
    dsl4jsb_Get_var2D_onChunk(TRAIT_, jmax_vcmax_ratio)

    climbuf => tile%parent_tile
    IF (TRIM(climbuf%name) /= 'veg') THEN
      CALL finish(routine, 'Usecase setup ERROR: CLIMBUF needs to run on the veg tile.')
    END IF
    dsl4jsb_Get_memory_tile(CLIMBUF_, climbuf)

    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, mat_prev_year, climbuf)
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, map_prev_year, climbuf)
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, max_mmtemp_prev_year, climbuf)
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, min_mmtemp_prev_year, climbuf)
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, reh_prev_year, climbuf)
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, soil_moist_prev_year, climbuf)
    dsl4jsb_Get_var2D_onChunk_tile_name(CLIMBUF_, nswr_prev_year, climbuf)

    CALL calc_trait_diagnostics(                        &
      & lib_id,                                         &
      & mat_prev_year_climbuf,                          &
      & map_prev_year_climbuf,                          &
      & min_mmtemp_prev_year_climbuf,                   &
      & max_mmtemp_prev_year_climbuf,                   &
      & reh_prev_year_climbuf,                          &
      & soil_moist_prev_year_climbuf,                   &
      & nswr_prev_year_climbuf,                         &
      & vcmax25_raw, vcmax25_pred,         &
      & jmax25_raw, jmax25_pred,           &
      & sla_raw, sla_pred,                 &
      & jmax_vcmax_ratio)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_trait

  SUBROUTINE aggregate_trait(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(TRAIT_)

    INTEGER :: iblk, ics, ice
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_trait'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(TRAIT_)
    weighted_by_fract => tile%Get_aggregator('weighted_by_fract')

    dsl4jsb_Aggregate_onChunk(TRAIT_, vcmax25_raw, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TRAIT_, vcmax25_pred, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TRAIT_, jmax25_raw, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TRAIT_, jmax25_pred, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TRAIT_, sla_raw, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TRAIT_, sla_pred, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TRAIT_, jmax_vcmax_ratio, weighted_by_fract)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_trait

#endif
END MODULE mo_trait_interface
