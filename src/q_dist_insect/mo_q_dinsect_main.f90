!> Update dinsect module
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
MODULE mo_q_dinsect_main
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message
  USE mo_jsb_math_constants,      ONLY: one_day, one_year, eps8, eps4
  USE mo_jsb_time,                ONLY: is_newyear
  USE mo_lnd_bgcm_idx
  USE mo_jsb_impl_constants,      ONLY: test_false_true
  USE mo_veg_config_class,        ONLY: get_number_of_veg_compartments
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_SEED_BED_POOL_ID, &
   & DINSECT_BGCM_DEADSAPWOOD_ID, DINSECT_BGCM_DEADHEARTWOOD_ID, DINSECT_BGCM_DEADLEAF_ID, &
   & DINSECT_BGCM_VEG_EXTRACT_ID, DINSECT_BGCM_FAST_LITTER_ID, &
   & DINSECT_BGCM_SLOW_LITTER_ID ,&
   & DINSECT_BGCM_LARVA_ID, DINSECT_BGCM_DROPPING_ID, &
   & DINSECT_BGCM_SEEDB_FASTL_ID

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: update_dinsect_barkbeetle, update_dinsect_defoliator

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_main'

CONTAINS

! ======================================================================================================= !
!>
!> Update insect outbreak effects.
!> for defoliator:   one pool: larva
!>                   two fluxes: 1. veg_extrac (veg_leaf->larva), 2. dropping (larva->soil)
! ----------------------------------------------------------------------------------------------------- !
SUBROUTINE update_dinsect_defoliator(tile, options)
  USE mo_jsb_class,               ONLY: Get_model
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,          ONLY: t_jsb_task_options
  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
  USE mo_jsb_process_class,       ONLY: VEG_, Q_DINSECT_, Q_PHENO_
  USE mo_q_dinsect_constants
#ifdef __QUINCY_STANDALONE__
  USE mo_qs_forcing,              ONLY: forcing_options
#else
  USE mo_jsb_time,                ONLY: get_year, get_year_day, get_secs_of_day
#endif
  USE mo_q_dinsect_damage
  USE mo_q_dinsect_defoliator_effect
  USE mo_q_dinsect_bark_beetle_effect
  USE mo_q_dinsect_util

  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Use_config(Q_DINSECT_)
  dsl4jsb_Use_config(VEG_)
  dsl4jsb_Use_memory(Q_DINSECT_)
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(Q_PHENO_)
  ! ----------------------------------------------------------------------------------------------------- !
  CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile       !< one tile with data structure for one lct
  TYPE(t_jsb_task_options),   INTENT(in)        :: options    !< model options
  TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
  ! ----------------------------------------------------------------------------------------------------- !
  TYPE(t_jsb_model),      POINTER       :: model
  LOGICAL                               :: is_growing_season(options%nc)
  TYPE(t_lctlib_element), POINTER       :: lctlib             !< land-cover-type library - parameter across pft's

  INTEGER                               :: iblk, ics, ice, ic, nc, iy

  INTEGER                               :: year, month, doy,  days_of_mon

  REAL(wp) :: hour
  LOGICAL::  start_of_year
  CHARACTER(LEN=200) :: filename

  LOGICAL                     :: is_bark_beetle, is_defoliator
  REAL(wp)                    :: dtime

  REAL(wp) :: yearly_damage_frac
  REAL(wp) :: yearly_damage_frac_base                        ! the last year's damage fraction, only used for bark beetle
  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_dinsect_defoliator'
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Def_mt2L2D :: veg_pool_mt
  dsl4jsb_Def_mt1L2D :: seed_bed_pool_mt
  dsl4jsb_Def_mt1L2D :: dinsect_deadsapwood_mt
  dsl4jsb_Def_mt1L2D :: dinsect_deadheartwood_mt
  dsl4jsb_Def_mt1L2D :: dinsect_deadleaf_mt
  dsl4jsb_Def_mt1L2D :: dinsect_larva_mt
  dsl4jsb_Def_mt1L2D :: dinsect_dropping_mt
  dsl4jsb_Def_mt2L2D :: dinsect_veg_extract_mt
  dsl4jsb_Def_mt2L2D :: dinsect_fast_litter_mt
  dsl4jsb_Def_mt1L2D :: dinsect_seedb_fastl_mt
  dsl4jsb_Def_mt2L2D :: dinsect_slow_litter_mt
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Def_config(Q_DINSECT_)
  dsl4jsb_Def_config(VEG_)
  dsl4jsb_Def_memory(Q_DINSECT_)
  dsl4jsb_Def_memory(VEG_)
  dsl4jsb_Def_memory(Q_PHENO_)
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Real2D_onChunk      :: mortality_rate                ! global used mortality rate
  dsl4jsb_Real2D_onChunk      :: defoliator_effects_triggered   ! counter for accumulating the sum of timesteps
  dsl4jsb_Real2D_onChunk      :: barkbeetle_effects_triggered
  dsl4jsb_Real2D_onChunk      :: defoliated_state              ! defoliation on or off
  dsl4jsb_Real2D_onChunk      :: is_insect_year                ! is insect disturbance year or not (1: true, 0: false)
  dsl4jsb_Real2D_onChunk      :: dinsect_year_since_start      ! the number of years since the start of insect disturbance
  dsl4jsb_Real2D_onChunk      :: dinsect_df_mort_dt            ! disturbed fraction caused by defoliator [fraction/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_total_mort_dt         ! total disturbed fraction caused by insects [fraction/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c     ! defoliator respiration C flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c13   ! defoliator respiration C13 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c14   ! defoliator respiration C14 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c   ! the defence used labile C from bark beetle attack [mol C/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c13 ! the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c14 ! the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_pool_total_c          ! total C from insect pools [mol C/m2]
  ! Q_PHENO_
  dsl4jsb_Real2D_onChunk      :: growing_season               !< Growing season indicator (1 = growing, 0 = not growing)
  ! ----------------------------------------------------------------------------------------------------- !
  iblk    = options%iblk
  ics     = options%ics
  ice     = options%ice
  nc      = options%nc
  dtime   = options%dtime
  ! ----------------------------------------------------------------------------------------------------- !
  IF (.NOT. tile%Is_process_active(Q_DINSECT_)) RETURN
  IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  ! ----------------------------------------------------------------------------------------------------- !
  model  => Get_model(tile%owner_model_id)
  lctlib => model%lctlib(tile%lcts(1)%lib_id)
  ! ----------------------------------------------------------------------------------------------------- !
  IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
  IF (lctlib%PastureFlag) RETURN !< do not run this routine at tiles like "pasture pft"
  IF (lctlib%CropFlag) RETURN !< do not run this routine at tiles like "crops pft"
  IF (lctlib%GrassFlag) RETURN !< do not run this routine at tiles like "grass pft"
  !< Needle-leaved evergreen trees (5) and breadleaved deciduous trees (3, 4)
  IF (lctlib%LctNumber < 3 .OR. lctlib%LctNumber > 5) RETURN !< Only run on landcover types 3, 4, 5
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Get_config(Q_DINSECT_)
  dsl4jsb_Get_config(VEG_)
  dsl4jsb_Get_memory(Q_DINSECT_)
  dsl4jsb_Get_memory(VEG_)
  dsl4jsb_Get_memory(Q_PHENO_)
  ! ----------------------------------------------------------------------------------------------------- !
  bgcm_store => tile%bgcm_store
  dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
  dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_POOL_ID, seed_bed_pool_mt)
  dsl4jsb_Get_mt2L2D(DINSECT_BGCM_VEG_EXTRACT_ID, dinsect_veg_extract_mt)
  dsl4jsb_Get_mt2L2D(DINSECT_BGCM_SLOW_LITTER_ID, dinsect_slow_litter_mt)
  dsl4jsb_Get_mt2L2D(DINSECT_BGCM_FAST_LITTER_ID, dinsect_fast_litter_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_SEEDB_FASTL_ID, dinsect_seedb_fastl_mt)

  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DEADSAPWOOD_ID, dinsect_deadsapwood_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DEADHEARTWOOD_ID, dinsect_deadheartwood_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DEADLEAF_ID, dinsect_deadleaf_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_LARVA_ID, dinsect_larva_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DROPPING_ID, dinsect_dropping_mt)

  ! ----------------------------------------------------------------------------------------------------- !
  ! Q_DINSECT_ 2D
  dsl4jsb_Get_var2D_onChunk(VEG_, mortality_rate)  !inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, defoliator_effects_triggered) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, barkbeetle_effects_triggered) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_df_mort_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_total_mort_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, defoliated_state) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, is_insect_year) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_year_since_start) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c14) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c14) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_pool_total_c) ! inout
  ! Q_PHENO_
  dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season) ! in

  ! transform REAL growing_season into LOGICAL is_growing_season
  DO ic = 1, nc
    IF (growing_season(ic) > test_false_true) THEN
      is_growing_season(ic) = .TRUE.
    ELSE
      is_growing_season(ic) = .FALSE.
    END IF
  END DO

#ifdef __QUINCY_STANDALONE__
        year = forcing_options%forcing_year
        doy = forcing_options%forcing_doy
        hour = forcing_options%forcing_hour
#else
        year = get_year(options%current_datetime)
        doy = INT(get_year_day(options%current_datetime))
        hour = get_secs_of_day(options%current_datetime) / 3600._wp ! seconds to hours
#endif


  ! ======================================================================================================= !
  !>
  !> Go science
  !>

  !<  ------------------  2.1 get defoliation_dt ----------------------------------
  dinsect_df_mort_dt = dsl4jsb_Config(Q_DINSECT_)%dinsect_dev_df_forcing_dist_frac(doy) / (one_day / dtime) ! distribute DF outbreak accross the whole day
  !<  ------------------  2.2 determine the defoliation  ----------------------------------

  ! Defoliated_state is reset in every timestep
  defoliated_state = 0.0_wp
  IF (lctlib%phenology_type .eq. 1) THEN  ! is evergreen needleleaf
    DO ic = 1, nc
        IF (is_growing_season(ic)) THEN !< is the tree in the growing season?
            IF ( dinsect_df_mort_dt(ic) > eps8) THEN ! is there an outbreak going on?
                defoliated_state(ic) = 1.0_wp ! the flag for controling leaf reflush is turned on
            ENDIF
        ENDIF
    END DO
  ENDIF

  CALL calc_defoliator_effects(&
     & nc                                     , & ! in
     & dtime                                  , & ! in
     & doy                                    , & ! in
     & dinsect_df_mort_dt(:)                   , & ! in
     & veg_pool_mt(:,:,:)                     , & ! in
     & dsl4jsb_Config(Q_DINSECT_)%dinsect_online, & ! in
     & dinsect_larva_mt(:,:)                  , & ! inout (pool)
     & dinsect_veg_extract_mt(:,:,:)          , & ! inout
     & dinsect_dropping_mt(:,:)               , & ! inout
     & dinsect_slow_litter_mt(:,:,:)          , & ! inout
     & dinsect_defoliator_resp_c(:)           , & ! inout
     & dinsect_defoliator_resp_c13(:)         , & ! inout
     & dinsect_defoliator_resp_c14(:)         , & ! inout
     & mask=(dinsect_df_mort_dt(:) > 0.0_wp)) ! optional(in)

  !< Used in mass conservation test
  dinsect_pool_total_c(:) = dinsect_pool_total_c(:) + dinsect_larva_mt(ixC,:)
  ! include defoliator mortality rate into total mortality rate
  dinsect_total_mort_dt = dinsect_total_mort_dt + dinsect_df_mort_dt
END SUBROUTINE update_dinsect_defoliator


! ======================================================================================================= !
!>
!> Update insect outbreak effects.
!> for bark beetle:  three dead pools: 1. dead_sapwood, 2.dead_heartwood, 3.dead_leaf;
!>                   three fluxes: 1. veg_extrac (veg->dead_pools), 2. fast_litter (veg->soil), 3. slow_litter (dead pool->soil)
! ----------------------------------------------------------------------------------------------------- !
SUBROUTINE update_dinsect_barkbeetle(tile, options)
  USE mo_jsb_class,               ONLY: Get_model
  USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
  USE mo_jsb_task_class,          ONLY: t_jsb_task_options
  USE mo_jsb_model_class,         ONLY: t_jsb_model
  USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
  USE mo_jsb_process_class,       ONLY: VEG_, Q_DINSECT_
  USE mo_q_dinsect_constants
#ifdef __QUINCY_STANDALONE__
  USE mo_qs_forcing,              ONLY: forcing_options
#else
  USE mo_jsb_time,                ONLY: get_year, get_year_day, get_secs_of_day
#endif
  USE mo_q_dinsect_damage
  USE mo_q_dinsect_defoliator_effect
  USE mo_q_dinsect_bark_beetle_effect
  USE mo_q_dinsect_util

  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Use_config(Q_DINSECT_)
  dsl4jsb_Use_config(VEG_)
  dsl4jsb_Use_memory(Q_DINSECT_)
  dsl4jsb_Use_memory(VEG_)
  ! ----------------------------------------------------------------------------------------------------- !
  CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile       !< one tile with data structure for one lct
  TYPE(t_jsb_task_options),   INTENT(in)        :: options    !< model options
  TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
  ! ----------------------------------------------------------------------------------------------------- !
  TYPE(t_jsb_model),      POINTER       :: model
  TYPE(t_lctlib_element), POINTER       :: lctlib             !< land-cover-type library - parameter across pft's

  INTEGER                               :: iblk, ics, ice, ic, nc, iy

  INTEGER                               :: year, month, doy,  days_of_mon

  REAL(wp) :: hour
  LOGICAL::  start_of_year
  CHARACTER(LEN=200) :: filename

  LOGICAL                     :: is_bark_beetle, is_defoliator
  REAL(wp)                    :: dtime

  REAL(wp) :: yearly_damage_frac
  REAL(wp) :: yearly_damage_frac_base                        ! the last year's damage fraction, only used for bark beetle
  CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_dinsect_barkbeetle'
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Def_mt2L2D :: veg_pool_mt
  dsl4jsb_Def_mt1L2D :: seed_bed_pool_mt
  dsl4jsb_Def_mt1L2D :: dinsect_deadsapwood_mt
  dsl4jsb_Def_mt1L2D :: dinsect_deadheartwood_mt
  dsl4jsb_Def_mt1L2D :: dinsect_deadleaf_mt
  dsl4jsb_Def_mt1L2D :: dinsect_larva_mt
  dsl4jsb_Def_mt1L2D :: dinsect_dropping_mt
  dsl4jsb_Def_mt2L2D :: dinsect_veg_extract_mt
  dsl4jsb_Def_mt2L2D :: dinsect_fast_litter_mt
  dsl4jsb_Def_mt1L2D :: dinsect_seedb_fastl_mt
  dsl4jsb_Def_mt2L2D :: dinsect_slow_litter_mt
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Def_config(Q_DINSECT_)
  dsl4jsb_Def_config(VEG_)
  dsl4jsb_Def_memory(Q_DINSECT_)
  dsl4jsb_Def_memory(VEG_)
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Real2D_onChunk      :: mortality_rate                ! global used mortality rate
  dsl4jsb_Real2D_onChunk      :: defoliator_effects_triggered   ! counter for accumulating the sum of timesteps
  dsl4jsb_Real2D_onChunk      :: barkbeetle_effects_triggered
  dsl4jsb_Real2D_onChunk      :: defoliated_state              ! defoliation on or off
  dsl4jsb_Real2D_onChunk      :: is_insect_year                ! is insect disturbance year or not (1: true, 0: false)
  dsl4jsb_Real2D_onChunk      :: dinsect_year_since_start      ! the number of years since the start of insect disturbance
  dsl4jsb_Real2D_onChunk      :: dinsect_bb_mort_dt            ! mortality rate caused by bark beetle [fraction/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_total_mort_dt         ! total mortality rate caused by insect outbreaks [fraction/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c     ! defoliator respiration C flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c13   ! defoliator respiration C13 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c14   ! defoliator respiration C14 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c   ! the defence used labile C from bark beetle attack [mol C/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c13 ! the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c14 ! the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_pool_total_c          ! total C from insect pools [mol C/m2]
  dsl4jsb_Real2D_onChunk      :: dinsect_bb_dist_fraction_read ! bark beetle disturbance fraction [fraction/day]
  ! ----------------------------------------------------------------------------------------------------- !
  iblk    = options%iblk
  ics     = options%ics
  ice     = options%ice
  nc      = options%nc
  dtime   = options%dtime
  ! ----------------------------------------------------------------------------------------------------- !
  IF (.NOT. tile%Is_process_active(Q_DINSECT_)) RETURN
  IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  ! ----------------------------------------------------------------------------------------------------- !
  model  => Get_model(tile%owner_model_id)
  lctlib => model%lctlib(tile%lcts(1)%lib_id)
  ! ----------------------------------------------------------------------------------------------------- !
  IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
  IF (lctlib%PastureFlag) RETURN !< do not run this routine at tiles like "pasture pft"
  IF (lctlib%CropFlag) RETURN !< do not run this routine at tiles like "crops pft"
  IF (lctlib%GrassFlag) RETURN !< do not run this routine at tiles like "grass pft"
  IF (lctlib%LctNumber /= 5) RETURN !< Only run on landcover type 5 (Needle-leaved evergreen trees)
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Get_config(Q_DINSECT_)
  dsl4jsb_Get_config(VEG_)
  dsl4jsb_Get_memory(Q_DINSECT_)
  dsl4jsb_Get_memory(VEG_)
  ! ----------------------------------------------------------------------------------------------------- !
  bgcm_store => tile%bgcm_store
  dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
  dsl4jsb_Get_mt1L2D(VEG_BGCM_SEED_BED_POOL_ID, seed_bed_pool_mt)
  dsl4jsb_Get_mt2L2D(DINSECT_BGCM_VEG_EXTRACT_ID, dinsect_veg_extract_mt)
  dsl4jsb_Get_mt2L2D(DINSECT_BGCM_SLOW_LITTER_ID, dinsect_slow_litter_mt)
  dsl4jsb_Get_mt2L2D(DINSECT_BGCM_FAST_LITTER_ID, dinsect_fast_litter_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_SEEDB_FASTL_ID, dinsect_seedb_fastl_mt)

  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DEADSAPWOOD_ID, dinsect_deadsapwood_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DEADHEARTWOOD_ID, dinsect_deadheartwood_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DEADLEAF_ID, dinsect_deadleaf_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_LARVA_ID, dinsect_larva_mt)
  dsl4jsb_Get_mt1L2D(DINSECT_BGCM_DROPPING_ID, dinsect_dropping_mt)

  ! ----------------------------------------------------------------------------------------------------- !
  ! Q_DINSECT_ 2D
  dsl4jsb_Get_var2D_onChunk(VEG_, mortality_rate)  !inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, defoliator_effects_triggered) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, barkbeetle_effects_triggered) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_bb_mort_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_total_mort_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, defoliated_state) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, is_insect_year) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_year_since_start) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c14) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c14) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_pool_total_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_bb_dist_fraction_read) ! in


#ifdef __QUINCY_STANDALONE__
      year = forcing_options%forcing_year
      doy = forcing_options%forcing_doy
      hour = forcing_options%forcing_hour
#else
      year = get_year(options%current_datetime)
      doy = INT(get_year_day(options%current_datetime))
      hour = get_secs_of_day(options%current_datetime) / 3600._wp ! seconds to hours
#endif

  IF (dsl4jsb_Config(Q_DINSECT_)%dinsect_flag_bb_from_namelist) THEN
    ! (fake) read input data (daily), import from namelist, apply same data to every year
    dinsect_bb_mort_dt(:) = dsl4jsb_Config(Q_DINSECT_)%dinsect_dev_bb_forcing_dist_frac(doy) / (one_day / dtime) ! distribute bb outbreak accross the whole day
  ELSE
    dinsect_bb_mort_dt(:) = dinsect_bb_dist_fraction_read(:) / (one_day / dtime)
  ENDIF

  ! ======================================================================================================= !
  !> Go science
  !>

  CALL calc_beetle_effects(&
         & nc                                     , & ! in
         & dtime                                  , & ! in
         & dinsect_bb_mort_dt(:)                  , & ! in
         & veg_pool_mt(:,:,:)                     , & ! in
         & seed_bed_pool_mt(:,:)                  , & ! in
         & dsl4jsb_Config(Q_DINSECT_)%dinsect_online, & ! in
         & dinsect_deadsapwood_mt(:,:)            , & ! inout (pool)
         & dinsect_deadheartwood_mt(:,:)          , & ! inout (pool)
         & dinsect_deadleaf_mt(:,:)               , & ! inout (pool)
         & dinsect_veg_extract_mt(:,:,:)          , & ! inout (flux)
         & dinsect_slow_litter_mt(:,:,:)          , & ! inout (flux)
         & dinsect_fast_litter_mt(:,:,:)          , & ! inout (flux)
         & dinsect_seedb_fastl_mt(:,:)            , & ! inout (flux)
         & dinsect_barkbeetle_labile_c(:)         , & ! inout
         & dinsect_barkbeetle_labile_c13(:)       , & ! inout
         & dinsect_barkbeetle_labile_c14(:)       , & ! inout
         & mask=(dinsect_bb_mort_dt(:) > 0.0_wp)) ! optional(in)

    ! Used in mass conservation test
    dinsect_pool_total_c(:) = dinsect_pool_total_c(:) + &
                              dinsect_deadsapwood_mt(ixC,:) + &
                              dinsect_deadheartwood_mt(ixC,:) + &
                              dinsect_deadleaf_mt(ixC,:)
    ! include bark beetle mortality rate into total mortality rate
    dinsect_total_mort_dt = dinsect_total_mort_dt + dinsect_bb_mort_dt
END SUBROUTINE update_dinsect_barkbeetle


#endif
END MODULE mo_q_dinsect_main
