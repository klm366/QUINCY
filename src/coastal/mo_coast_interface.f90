! Author: Dr. Sri Nandini-Weiß
! Date: 2026-09-23
!> Interface to the coastal-constraint process
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
!>#### Calculates the coastal control factors on the mangrove tile
!>
!> The process runs on the land cover type that represents mangroves, identified by
!> its lctlib id, and is skipped on every other tile. The control factors are written
!> to the process memory; whether they are applied to gross assimilation is decided in
!> mo_assimi_interface from the namelist switches of this process.
!>
MODULE mo_coast_interface
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp
  USE mo_jsb_control,       ONLY: debug_on
  USE mo_exception,         ONLY: message
  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_class,         ONLY: Get_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class, ONLY: t_jsb_process
  USE mo_jsb_task_class,    ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes COAST_
  dsl4jsb_Use_config(COAST_)
  dsl4jsb_Use_memory(COAST_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_coast_tasks

  !> lctlib id of the mangrove land cover type. In the npft=13 usecase the tiles
  !> pft01 to pft13 are created with the ids 1,2,3,4,9,10,11,12,14,15,19,20,21, so the
  !> mangrove clone is 21. The coastal constraints are specific to mangroves and are
  !> not applied to any other type.
  INTEGER, PARAMETER :: mangrove_lct_id = 21

  CHARACTER(len=*), PARAMETER :: modname = 'mo_coast_interface'

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_coast
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_coast
    PROCEDURE, NOPASS :: Aggregate => aggregate_coast
  END TYPE tsk_coast

  INTERFACE tsk_coast
    PROCEDURE Create_task_coast
  END INTERFACE tsk_coast

CONTAINS

  ! ------------------------------------------------------------------------------------------------------ !
  !> Constructor for the coastal-constraint task
  !>
  FUNCTION Create_task_coast(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_coast::return_ptr)
    CALL return_ptr%Construct(name='coast', process_id=COAST_, owner_model_id=model_id)

  END FUNCTION Create_task_coast

  ! ------------------------------------------------------------------------------------------------------ !
  !> Register tasks for the coastal-constraint process
  !>
  SUBROUTINE Register_coast_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_coast(model_id))

  END SUBROUTINE Register_coast_tasks

  ! ======================================================================================================== !
  !>
  !>#### Update the coastal control factors
  !>
  !> We here evaluate the salinity and tidal-inundation control factors of Zhao et al.
  !> (2025) from the boundary conditions held in the process memory. The aboveground
  !> root biomasses that enter the ventilation term are not yet available as prognostic
  !> variables; until they are, the ventilation thresholds of the configuration are
  !> used, which corresponds to a stand with a fully developed aerial root system.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_coast(tile, options)

    USE mo_coast_process, ONLY: calc_salinity_control, calc_inundation_control, &
      &                         calc_bgr_ventilated, calc_coastal_factor

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    INTEGER :: iblk, ics, ice, nc, ic
    INTEGER :: lib_id
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_coast'

    dsl4jsb_Def_config(COAST_)
    dsl4jsb_Def_memory(COAST_)

    dsl4jsb_Real2D_onChunk :: salinity
    dsl4jsb_Real2D_onChunk :: tide_height
    dsl4jsb_Real2D_onChunk :: bgr_fract_inundated
    dsl4jsb_Real2D_onChunk :: bgr_fract_ventilated
    dsl4jsb_Real2D_onChunk :: bgr_fract_anoxic
    dsl4jsb_Real2D_onChunk :: control_salinity
    dsl4jsb_Real2D_onChunk :: control_inundation
    dsl4jsb_Real2D_onChunk :: control_coastal
    ! ------------------------------------------------------------------------------------------------------ !

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(COAST_)) RETURN

    ! The coastal constraints apply to the mangrove land cover type only.
    IF (.NOT. ALLOCATED(tile%lcts)) RETURN
    lib_id = tile%lcts(1)%lib_id
    IF (lib_id /= mangrove_lct_id) RETURN

    IF (debug_on() .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(COAST_)
    dsl4jsb_Get_memory(COAST_)

    dsl4jsb_Get_var2D_onChunk(COAST_, salinity)              ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_height)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, bgr_fract_inundated)   ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, bgr_fract_ventilated)  ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, bgr_fract_anoxic)      ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, control_salinity)      ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, control_inundation)    ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, control_coastal)       ! out

    DO ic = 1, nc

      control_salinity(ic) = calc_salinity_control(          &
        & salinity(ic),                                      &
        & dsl4jsb_Config(COAST_)%cs_min,                     &
        & dsl4jsb_Config(COAST_)%l_normalise_cs)

      bgr_fract_ventilated(ic) = calc_bgr_ventilated(        &
        & dsl4jsb_Config(COAST_)%b_agrs_lim,                 &
        & dsl4jsb_Config(COAST_)%b_agrp_lim,                 &
        & tide_height(ic),                                   &
        & dsl4jsb_Config(COAST_)%b_agrs_lim,                 &
        & dsl4jsb_Config(COAST_)%b_agrp_lim,                 &
        & dsl4jsb_Config(COAST_)%vent_weight_s,              &
        & dsl4jsb_Config(COAST_)%vent_weight_p,              &
        & dsl4jsb_Config(COAST_)%tree_density)

      bgr_fract_anoxic(ic) = MIN(1._wp, MAX(0._wp,           &
        & bgr_fract_inundated(ic) - bgr_fract_ventilated(ic)))

      control_inundation(ic) = calc_inundation_control(      &
        & bgr_fract_inundated(ic),                           &
        & bgr_fract_ventilated(ic),                          &
        & dsl4jsb_Config(COAST_)%cin_min)

      control_coastal(ic) = calc_coastal_factor(             &
        & control_salinity(ic),                              &
        & control_inundation(ic),                            &
        & dsl4jsb_Config(COAST_)%l_apply_salinity_stress,    &
        & dsl4jsb_Config(COAST_)%l_apply_tidal_stress)

    END DO

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_coast

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of the coastal control factors
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE aggregate_coast(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(COAST_)

    INTEGER :: iblk, ics, ice
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_coast'
    ! ------------------------------------------------------------------------------------------------------ !

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(COAST_)
    weighted_by_fract => tile%Get_aggregator('weighted_by_fract')

    dsl4jsb_Aggregate_onChunk(COAST_, salinity,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, tide_height,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, bgr_fract_inundated,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, bgr_fract_ventilated, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, bgr_fract_anoxic,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, control_salinity,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, control_inundation,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, control_coastal,      weighted_by_fract)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_coast

#endif
END MODULE mo_coast_interface
