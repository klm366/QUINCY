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
!> to the process memory every time step; they are multiplied into gross assimilation
!> in mo_assimi_interface when one of the application switches of this process is set.
!>
!> The routines contain no OpenACC directives and run on the host, as the trait process
!> does. A GPU build would need data regions here and ROUTINE SEQ in mo_coast_process.
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
  USE mo_coast_constants,   ONLY: mangrove_lct_id, ntide

  dsl4jsb_Use_processes COAST_, HYDRO_
  dsl4jsb_Use_config(COAST_)
  dsl4jsb_Use_memory(COAST_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_coast_tasks


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
  !> We here evaluate, for the current time step, the tide height from the harmonic
  !> constants, the inundated, ventilated and anoxic fractions of belowground root, the
  !> salinity and tidal-inundation control factors of Zhao et al. (2025), and their
  !> product, the scalar coastal stress factor control_coastal. Only cells inside
  !> coast_mask are evaluated; about 600 of the 20480 cells of R02B04 carry mangrove.
  !>
  !> The aboveground root biomasses in the ventilation term are not prognostic in JSBACH.
  !> The ventilation thresholds of the configuration are passed in their place, which
  !> corresponds to a stand whose aerial root system is fully developed; Zhao et al.
  !> reach this state after about five years.
  !>
  ! -------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_coast(tile, options)

    USE mo_jsb_time,      ONLY: get_year, get_month, get_day, get_secs_of_day
    USE mo_coast_util,    ONLY: days_since_2000
    USE mo_coast_process, ONLY: calc_salinity_control, calc_tide_arguments, calc_tide_height,     &
      &                         calc_bgr_inundated, calc_bgr_ventilated, calc_inundation_control, &
      &                         calc_coastal_factor

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    INTEGER  :: iblk, ics, ice, nc, ic
    INTEGER  :: lib_id
    REAL(wp) :: time_hours            !< Model time since 1 January 2000 [h]
    REAL(wp) :: hre(ntide), him(ntide)          !< Complex amplitudes of the constituents in one cell [m]
    REAL(wp) :: cos_arg(ntide), sin_arg(ntide)  !< Phase terms of the constituents at the model time []
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_coast'

    dsl4jsb_Def_config(COAST_)
    dsl4jsb_Def_memory(COAST_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Real2D_onChunk :: coast_mask
    dsl4jsb_Real2D_onChunk :: salinity
    dsl4jsb_Real2D_onChunk :: tide_hre_m2
    dsl4jsb_Real2D_onChunk :: tide_him_m2
    dsl4jsb_Real2D_onChunk :: tide_hre_s2
    dsl4jsb_Real2D_onChunk :: tide_him_s2
    dsl4jsb_Real2D_onChunk :: tide_hre_n2
    dsl4jsb_Real2D_onChunk :: tide_him_n2
    dsl4jsb_Real2D_onChunk :: tide_hre_k2
    dsl4jsb_Real2D_onChunk :: tide_him_k2
    dsl4jsb_Real2D_onChunk :: tide_hre_k1
    dsl4jsb_Real2D_onChunk :: tide_him_k1
    dsl4jsb_Real2D_onChunk :: tide_hre_o1
    dsl4jsb_Real2D_onChunk :: tide_him_o1
    dsl4jsb_Real2D_onChunk :: tide_hre_p1
    dsl4jsb_Real2D_onChunk :: tide_him_p1
    dsl4jsb_Real2D_onChunk :: tide_hre_q1
    dsl4jsb_Real2D_onChunk :: tide_him_q1
    dsl4jsb_Real2D_onChunk :: tide_height
    dsl4jsb_Real2D_onChunk :: bgr_fract_inundated
    dsl4jsb_Real2D_onChunk :: bgr_fract_ventilated
    dsl4jsb_Real2D_onChunk :: bgr_fract_anoxic
    dsl4jsb_Real2D_onChunk :: control_salinity
    dsl4jsb_Real2D_onChunk :: control_inundation
    dsl4jsb_Real2D_onChunk :: control_coastal
    dsl4jsb_Real2D_onChunk :: root_depth
    ! ------------------------------------------------------------------------------------------------------ !

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(COAST_)) RETURN

    ! The coastal constraints apply to the mangrove land cover type only
    IF (.NOT. ALLOCATED(tile%lcts)) RETURN
    lib_id = tile%lcts(1)%lib_id
    IF (lib_id /= mangrove_lct_id) RETURN

    IF (debug_on() .AND. iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(COAST_)
    dsl4jsb_Get_memory(COAST_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(COAST_, coast_mask)            ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, salinity)              ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_m2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_m2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_s2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_s2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_n2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_n2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_k2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_k2)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_k1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_k1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_o1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_o1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_p1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_p1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_hre_q1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_him_q1)           ! in
    dsl4jsb_Get_var2D_onChunk(COAST_, tide_height)           ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, bgr_fract_inundated)   ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, bgr_fract_ventilated)  ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, bgr_fract_anoxic)      ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, control_salinity)      ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, control_inundation)    ! out
    dsl4jsb_Get_var2D_onChunk(COAST_, control_coastal)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, root_depth)            ! in

    ! Model time, counted in hours from 1 January 2000
    time_hours = 24._wp * REAL(days_since_2000(get_year(options%current_datetime),            &
      &                                        get_month(options%current_datetime),           &
      &                                        get_day(options%current_datetime)), wp)        &
      &        + get_secs_of_day(options%current_datetime) / 3600._wp

    ! The phase terms depend on the time only, so they are evaluated once for all cells
    CALL calc_tide_arguments(time_hours, cos_arg, sin_arg)

    DO ic = 1, nc

      ! Cells without mangrove: neutral values, the equations are not evaluated
      IF (coast_mask(ic) < 0.5_wp) THEN
        tide_height(ic)          = 0._wp
        bgr_fract_inundated(ic)  = 0._wp
        bgr_fract_ventilated(ic) = 0._wp
        bgr_fract_anoxic(ic)     = 0._wp
        control_salinity(ic)     = 1._wp
        control_inundation(ic)   = 1._wp
        control_coastal(ic)      = 1._wp
        CYCLE
      END IF

      hre = [ tide_hre_m2(ic), tide_hre_s2(ic), tide_hre_n2(ic), tide_hre_k2(ic), &
        &     tide_hre_k1(ic), tide_hre_o1(ic), tide_hre_p1(ic), tide_hre_q1(ic) ]
      him = [ tide_him_m2(ic), tide_him_s2(ic), tide_him_n2(ic), tide_him_k2(ic), &
        &     tide_him_k1(ic), tide_him_o1(ic), tide_him_p1(ic), tide_him_q1(ic) ]

      tide_height(ic) = calc_tide_height(hre, him, cos_arg, sin_arg)

      bgr_fract_inundated(ic) = calc_bgr_inundated(tide_height(ic), root_depth(ic))

      bgr_fract_ventilated(ic) = calc_bgr_ventilated(         &
        & dsl4jsb_Config(COAST_)%b_agrs_lim,                  &
        & dsl4jsb_Config(COAST_)%b_agrp_lim,                  &
        & tide_height(ic),                                    &
        & dsl4jsb_Config(COAST_)%b_agrs_lim,                  &
        & dsl4jsb_Config(COAST_)%b_agrp_lim,                  &
        & dsl4jsb_Config(COAST_)%vent_weight_s,               &
        & dsl4jsb_Config(COAST_)%vent_weight_p,               &
        & dsl4jsb_Config(COAST_)%tree_density)

      bgr_fract_anoxic(ic) = MIN(1._wp, MAX(0._wp, bgr_fract_inundated(ic) - bgr_fract_ventilated(ic)))

      control_inundation(ic) = calc_inundation_control(       &
        & bgr_fract_inundated(ic),                            &
        & bgr_fract_ventilated(ic),                           &
        & dsl4jsb_Config(COAST_)%cin_min)

      control_salinity(ic) = calc_salinity_control(           &
        & salinity(ic),                                       &
        & dsl4jsb_Config(COAST_)%cs_min,                      &
        & dsl4jsb_Config(COAST_)%l_normalise_cs)

      control_coastal(ic) = calc_coastal_factor(              &
        & control_salinity(ic),                               &
        & control_inundation(ic),                             &
        & dsl4jsb_Config(COAST_)%l_apply_salinity_stress,     &
        & dsl4jsb_Config(COAST_)%l_apply_tidal_stress)

    END DO

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_coast

  ! ======================================================================================================== !
  !>
  !>#### Aggregation of the coastal control factors
  !>
  !> Only the three control factors are aggregated. They are one on every tile other
  !> than the mangrove tile, so their area-weighted mean is the effective factor of the
  !> parent tile. The boundary conditions, the tide height and the root fractions are
  !> zero on the other tiles; their area-weighted mean would be diluted by the mangrove
  !> fraction and has no physical meaning, so they are read on the mangrove tile only.
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

    dsl4jsb_Aggregate_onChunk(COAST_, control_salinity,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, control_inundation,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(COAST_, control_coastal,      weighted_by_fract)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_coast

#endif
END MODULE mo_coast_interface
