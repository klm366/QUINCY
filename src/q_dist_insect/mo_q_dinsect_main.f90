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

  PUBLIC :: update_dinsect_effects, update_dinsect_barkbeetle, update_dinsect_defoliator

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_main'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Update insect outbreak effects.
  !> for bark beetle:  three dead pools: 1. dead_sapwood, 2.dead_heartwood, 3.dead_leaf;
  !>                   three fluxes: 1. veg_extrac (veg->dead_pools), 2. fast_litter (veg->soil), 3. slow_litter (dead pool->soil)
  !> for defoliator:   one pool: larva
  !>                   two fluxes: 1. veg_extrac (veg_leaf->larva), 2. dropping (larva->soil)
  ! ----------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_dinsect_effects(tile, options)
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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_dinsect_effects'
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
    dsl4jsb_Real2D_onChunk      :: dinsect_occupy_frac           ! the fraction of snag occupyed area, which cause the waste of seed.
    dsl4jsb_Real2D_onChunk      :: defoliated_state              ! defoliation on or off
    dsl4jsb_Real2D_onChunk      :: is_insect_year                ! is insect disturbance year or not (1: true, 0: false)
    dsl4jsb_Real2D_onChunk      :: dinsect_year_since_start      ! the number of years since the start of insect disturbance
    dsl4jsb_Real2D_onChunk      :: dinsect_mort_dt               ! addtional mortality rate caused by bark beetle
    dsl4jsb_Real2D_onChunk      :: dinsect_defoli_dt             ! defoliation rate caused by defoliator  !dsl4jsb_Real2D_onChunk      :: heartwood_C_zero              ! the amount of heartwood C at the first timestep of bark beetle event
    dsl4jsb_Real2D_onChunk      :: dinsect_flux_frac             ! the fraction of disturbed fraction of target variables.
    dsl4jsb_Real2D_onChunk      :: dinsect_disturbable_frac      ! the fraction of remaining undisturbable tree resources.
    dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c     ! defoliator respiration C flux [mol/m2/timestep]
    dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c13   ! defoliator respiration C13 flux [mol/m2/timestep]
    dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c14   ! defoliator respiration C14 flux [mol/m2/timestep]
    dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c   ! the defence used labile C from bark beetle attack [mol C/m2/timestep]
    dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c13 ! the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
    dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c14 ! the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
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
    !< DF: Needle-leaved evergreen trees (5) and breadleaved deciduous trees (3, 4)
    !< BB: Needle-leaved evergreen trees (5)
    IF (lctlib%LctNumber < 3 .OR. lctlib%LctNumber > 5) RETURN !< Only run on landcover types 3, 4, 5
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
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_occupy_frac) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_mort_dt) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_disturbable_frac) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoli_dt) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_flux_frac) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, defoliated_state) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, is_insect_year) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_year_since_start) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c13) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c14) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c13) ! inout
    dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c14) ! inout

    ! ======================================================================================================= !
    !> Go science
    !>
    !> Get tree mortality for bark beelte OR leaf defoliation rate for defoliator from prescribed data.
    !> In each timestep, additional tree mortality OR defoliation rate is used to calculate effects for two insects, respectively.
    !> Insect-related pools and fluxes are updated.
    !> Vegetation and soil dynamics can be changed by insect effect module if dinsect_online is .True.
    ! ----------------------------------------------------------------------------------------------------- !

#ifdef __QUINCY_STANDALONE__
    year = forcing_options%forcing_year
    doy = forcing_options%forcing_doy
    hour = forcing_options%forcing_hour
#else
    year = get_year(options%current_datetime)
    doy = INT(get_year_day(options%current_datetime))
    hour = get_secs_of_day(options%current_datetime) / 3600._wp ! seconds to hours
#endif

    start_of_year =  is_newyear(options%current_datetime, dtime)

    WRITE(*,*) "dinsect_main:: running... (type=), lctnumber=", &
    dsl4jsb_Config(Q_DINSECT_)%dinsect_type, lctlib%LctNumber
    ! if prescribe is true, the source of flux is either from namelist or from default values
    ! In Phase 2, it can be predicted by insect population model
    DO ic = 1, nc

        ! Determine whether is insect disturbance year or not
        ! After checking year boundaries
        IF (start_of_year) THEN

          is_insect_year(:)=0.0_wp

          IF (dsl4jsb_Config(Q_DINSECT_)%dinsect_prescribe .eqv. .True.) THEN

            ! For now, whether all grids have the the same insect disturbance year. Later, can read in file for each grid.
            DO iy = 1, dsl4jsb_Config(Q_DINSECT_)%dinsect_year_span
              ! Only generate file for prescribed outbreak years

              IF (year == (dsl4jsb_Config(Q_DINSECT_)%dinsect_start_year + (iy - 1))) THEN
                ! Print information only once at the first timestep of each year
                IF (iy == 1) THEN
                  CALL print_dinsect_info(&
                      & dsl4jsb_Config(Q_DINSECT_)%dinsect_type, &
                      & dsl4jsb_Config(Q_DINSECT_)%dinsect_prescribe, &
                      & dsl4jsb_Config(Q_DINSECT_)%dinsect_from_namelist, &
                      & dsl4jsb_Config(Q_DINSECT_)%dinsect_start_year, &
                      & dsl4jsb_Config(Q_DINSECT_)%dinsect_year_span)
                END IF

                is_insect_year(ic)=1.0_wp    ! update the flag for insect disturbance year
                dinsect_year_since_start(ic)=dinsect_year_since_start(ic)+1.0_wp

              END IF

            END DO

          END IF

        END IF

        ! For now, the damage fraction of each grid is the same. Later, can read in file for each grid.

        ! If the disturbance year is determined, derive the flux_frac for current timestep
        IF (is_insect_year(ic) > eps4) THEN

          IF (dsl4jsb_Config(Q_DINSECT_)%dinsect_prescribe .eqv. .True.) THEN

            IF (dsl4jsb_Config(Q_DINSECT_)%dinsect_from_namelist) THEN

              IF (ABS(dinsect_year_since_start(ic) - 1.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year1
                yearly_damage_frac_base = 0.0_wp
              ELSE IF (ABS(dinsect_year_since_start(ic) - 2.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year2
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year1
              ELSE IF (ABS(dinsect_year_since_start(ic) - 3.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year3
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year2
              ELSE IF (ABS(dinsect_year_since_start(ic) - 4.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year4
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year3
              ELSE IF (ABS(dinsect_year_since_start(ic) - 5.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year5
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year4
              ELSE IF (ABS(dinsect_year_since_start(ic) - 6.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year6
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year5
              ELSE IF (ABS(dinsect_year_since_start(ic) - 7.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year7
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year6
              ELSE IF (ABS(dinsect_year_since_start(ic) - 8.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year8
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year7
              ELSE IF (ABS(dinsect_year_since_start(ic) - 9.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year9
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year8
              ELSE IF (ABS(dinsect_year_since_start(ic) - 10.0_wp) < eps8) THEN
                yearly_damage_frac = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year10
                yearly_damage_frac_base = dsl4jsb_Config(Q_DINSECT_)%dinsect_damage_year9
              ELSE
                yearly_damage_frac = 0.0_wp
                yearly_damage_frac_base = 0.0_wp
              END IF

              ! If the insect is bark beetle, the damage fraction is calculated as the relative damage fraction based on the last year's remaining pool
              IF (dsl4jsb_Config(Q_DINSECT_)%dinsect_type == 1) THEN
                yearly_damage_frac = (yearly_damage_frac - yearly_damage_frac_base) / (1.0_wp - yearly_damage_frac_base)
              END IF

              CALL derive_flux_frac_with_sigmoid_from_nml(&
                & dsl4jsb_Config(Q_DINSECT_)%dinsect_start_doy, &
                & dsl4jsb_Config(Q_DINSECT_)%dinsect_duration, &
                & year, &
                & doy, &
                & hour, &
                & yearly_damage_frac, &
                & dinsect_flux_frac(ic))
            ELSE
              CALL message("update_dinsects_effects", &
                "dinsect_from_namelist is OFF! Default values are used!")
              dinsect_flux_frac(ic) = dinsect_const_damage_frac
            END IF

          END IF

        ELSE

          dinsect_flux_frac(ic) = 0.0_wp

        END IF

    END DO


    ! Determine Insect Type
    ! TOOD:  only on evergreen needle leaf (PFT 5)
    is_bark_beetle = (dsl4jsb_Config(Q_DINSECT_)%dinsect_type == 1)
    ! TODO:  evergreen needle leaf(restricted regrowth during the defoliated year) (PFT 5) ;
    !        broad deciduous leaf (regrowth anytime during proper season),
    !        PFT 3 (rain) and 4 (summer)
    is_defoliator = (dsl4jsb_Config(Q_DINSECT_)%dinsect_type == 2)

    ! Insect Specific Calculations
    IF (is_bark_beetle) THEN
        !< -----------------------   1.1 get mortality_dt ---------------------------------
        dinsect_mort_dt(:) =  dinsect_flux_frac(:)
        !< ------------------ 1.2 Calculate effects (after bark beetle event happens) --------------------
        DO ic = 1, nc

          IF (dinsect_mort_dt(ic) > eps8) THEN
              barkbeetle_effects_triggered(ic) = barkbeetle_effects_triggered(ic) + 1.0_wp
          END IF

          IF ((ABS(barkbeetle_effects_triggered(ic)) - 1.0_wp) < eps4) THEN
              dinsect_disturbable_frac(ic) = 1.0_wp
          END IF

        END DO

        CALL calc_beetle_effects(&
        & nc                                     , & ! in
        & dtime                                  , & ! in
        & dinsect_mort_dt(:)                     , & ! in
        & veg_pool_mt(:,:,:)                     , & ! in
        & seed_bed_pool_mt(:,:)                  , & ! in
        & dsl4jsb_Config(Q_DINSECT_)%dinsect_online, & ! in
        !& heartwood_C_zero(:)                    , & ! in
        & dinsect_deadsapwood_mt(:,:)            , & ! inout
        & dinsect_deadheartwood_mt(:,:)          , & ! inout
        & dinsect_deadleaf_mt(:,:)               , & ! inout
        & dinsect_veg_extract_mt(:,:,:)          , & ! inout
        & dinsect_slow_litter_mt(:,:,:)          , & ! inout
        & dinsect_fast_litter_mt(:,:,:)          , & ! inout
        & dinsect_seedb_fastl_mt(:,:)   , & ! inout
        & dinsect_barkbeetle_labile_c(:)         , & ! inout
        & dinsect_barkbeetle_labile_c13(:)       , & ! inout
        & dinsect_barkbeetle_labile_c14(:)       , & ! inout
        & dinsect_occupy_frac(:)                 , & ! inout
        & dinsect_disturbable_frac(:)            , & ! inout
        & mask=(barkbeetle_effects_triggered(:) > 0.0_wp)) ! optional(in)

        DO ic = 1, nc
          IF (barkbeetle_effects_triggered(ic) > 0.0_wp ) THEN
            WRITE(*,*) "dinsect_main.f90:: barkbeetle_effects_triggered, dinsect_disturbable_frac=", &
                      barkbeetle_effects_triggered(ic), dinsect_disturbable_frac
          END IF
        ENDDO

    ELSE IF (is_defoliator) THEN
        !<  ------------------  2.1 get defoliation_dt ----------------------------------
        dinsect_defoli_dt(:) = dinsect_flux_frac(:)

        !<  ------------------  2.2 determine the defoliation  ----------------------------------
        !CALL get_time_info(year, month, days_of_mon, doy)

        ! Initialize or reset flags
        ! defoliated_state is reset in every timestep
        defoliated_state(:) = 0.0_wp

        DO ic = 1, nc

          IF (dinsect_defoli_dt(ic) > eps8) THEN  ! Once defoliation happened, this flag will be turned on
            defoliated_state(ic) = 1.0_wp  ! Update at each timestep in mo_q_veg_growth.f
          END IF

          IF (is_insect_year(ic) > eps8) THEN  ! Outbreak happend in this year
            IF (lctlib%phenology_type .eq. 1) THEN  ! is evergreen needleleaf
                ! (1) the regrowth restriction should be adjusted, if required, according to the research
                ! domain and plant leaf growth strategy; (2) the distribution of grid-based damage fraction
                ! among PFTs should be pre-calculated before inputing (this point will be more easier
                ! and clearer in coding and flexible in scientific using, but need to ask Ana firstly)
                !
                ! TODO: "prevent_late_summer_leaf_regrowth" (set the northern hemisphere as > =July 1 and
                ! southern hemisphere as <= June 30) and make the evergreen PFT got this feature ?
                ! Something like this
                IF (int(month) > 6) THEN  ! any months after defoliation (here we generally say June)
                  defoliated_state(ic) = 1.0_wp  ! the flag for controling leaf reflush will be on, until next year
                END IF
            END IF
          END IF

          !< ---------------  2.3 Calculate effects (after defoliation event happens)-------------------------------
          IF (dinsect_defoli_dt(ic) > eps8) THEN
              defoliator_effects_triggered(ic) = defoliator_effects_triggered(ic) + 1.0_wp

          END IF

        END DO

        CALL calc_defoliator_effects(&
        & nc                                     , & ! in
        & dtime                                  , & ! in
        & doy                                    , & ! in
        & dinsect_defoli_dt(:)                   , & ! in
        & veg_pool_mt(:,:,:)                     , & ! in
        & dsl4jsb_Config(Q_DINSECT_)%dinsect_online, & ! in
        & dinsect_larva_mt(:,:)                  , & ! inout
        & dinsect_veg_extract_mt(:,:,:)          , & ! inout
        & dinsect_dropping_mt(:,:)               , & ! inout
        & dinsect_slow_litter_mt(:,:,:)          , & ! inout
        & dinsect_defoliator_resp_c(:)           , & ! inout
        & dinsect_defoliator_resp_c13(:)         , & ! inout
        & dinsect_defoliator_resp_c14(:)         , & ! inout
        & mask=(defoliator_effects_triggered(:) > 0.0_wp)) ! optional(in)
    END IF

  END SUBROUTINE update_dinsect_effects


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
  ! ----------------------------------------------------------------------------------------------------- !
  dsl4jsb_Real2D_onChunk      :: mortality_rate                ! global used mortality rate
  dsl4jsb_Real2D_onChunk      :: defoliator_effects_triggered   ! counter for accumulating the sum of timesteps
  dsl4jsb_Real2D_onChunk      :: barkbeetle_effects_triggered
  dsl4jsb_Real2D_onChunk      :: dinsect_occupy_frac           ! the fraction of snag occupyed area, which cause the waste of seed.
  dsl4jsb_Real2D_onChunk      :: defoliated_state              ! defoliation on or off
  dsl4jsb_Real2D_onChunk      :: is_insect_year                ! is insect disturbance year or not (1: true, 0: false)
  dsl4jsb_Real2D_onChunk      :: dinsect_year_since_start      ! the number of years since the start of insect disturbance
  dsl4jsb_Real2D_onChunk      :: dinsect_mort_dt               ! addtional mortality rate caused by bark beetle
  dsl4jsb_Real2D_onChunk      :: dinsect_defoli_dt             ! defoliation rate caused by defoliator  !dsl4jsb_Real2D_onChunk      :: heartwood_C_zero              ! the amount of heartwood C at the first timestep of bark beetle event
  dsl4jsb_Real2D_onChunk      :: dinsect_flux_frac             ! the fraction of disturbed fraction of target variables.
  dsl4jsb_Real2D_onChunk      :: dinsect_disturbable_frac      ! the fraction of remaining undisturbable tree resources.
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c     ! defoliator respiration C flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c13   ! defoliator respiration C13 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c14   ! defoliator respiration C14 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c   ! the defence used labile C from bark beetle attack [mol C/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c13 ! the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c14 ! the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
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
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_occupy_frac) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_mort_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_disturbable_frac) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoli_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_flux_frac) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, defoliated_state) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, is_insect_year) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_year_since_start) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c14) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c14) ! inout

  ! ======================================================================================================= !
  !>
  !> Go science
  !>

  CALL derive_flux_frac_with_sigmoid_from_nml(&
    & dsl4jsb_Config(Q_DINSECT_)%dinsect_start_doy, &
    & dsl4jsb_Config(Q_DINSECT_)%dinsect_duration, &
    & year, &
    & doy, &
    & hour, &
    & yearly_damage_frac, &
    & dinsect_flux_frac(ic))

  !<  ------------------  2.1 get defoliation_dt ----------------------------------
  dinsect_defoli_dt(:) = dinsect_flux_frac(:)

  !<  ------------------  2.2 determine the defoliation  ----------------------------------
  !CALL get_time_info(year, month, days_of_mon, doy)

  ! Initialize or reset flags
  ! defoliated_state is reset in every timestep
  defoliated_state(:) = 0.0_wp

  DO ic = 1, nc

    IF (dinsect_defoli_dt(ic) > eps8) THEN  ! Once defoliation happened, this flag will be turned on
      defoliated_state(ic) = 1.0_wp  ! Update at each timestep in mo_q_veg_growth.f
    END IF

    IF (is_insect_year(ic) > eps8) THEN  ! Outbreak happend in this year
      IF (lctlib%phenology_type .eq. 1) THEN  ! is evergreen needleleaf
        ! (1) the regrowth restriction should be adjusted, if required, according to the research
        ! domain and plant leaf growth strategy; (2) the distribution of grid-based damage fraction
        ! among PFTs should be pre-calculated before inputing (this point will be more easier
        ! and clearer in coding and flexible in scientific using, but need to ask Ana firstly)
        !
        ! TODO: "prevent_late_summer_leaf_regrowth" (set the northern hemisphere as > =July 1 and
        ! southern hemisphere as <= June 30) and make the evergreen PFT got this feature ?
        ! Something like this
        IF (int(month) > 6) THEN  ! any months after defoliation (here we generally say June)
          defoliated_state(ic) = 1.0_wp  ! the flag for controling leaf reflush will be on, until next year
        END IF
      END IF
    END IF

    !< ---------------  2.3 Calculate effects (after defoliation event happens)-------------------------------
    IF (dinsect_defoli_dt(ic) > eps8) THEN
        defoliator_effects_triggered(ic) = defoliator_effects_triggered(ic) + 1.0_wp
    END IF
  END DO

  CALL calc_defoliator_effects(&
    & nc                                     , & ! in
    & dtime                                  , & ! in
    & doy                                    , & ! in
    & dinsect_defoli_dt(:)                   , & ! in
    & veg_pool_mt(:,:,:)                     , & ! in
    & dsl4jsb_Config(Q_DINSECT_)%dinsect_online, & ! in
    & dinsect_larva_mt(:,:)                  , & ! inout
    & dinsect_veg_extract_mt(:,:,:)          , & ! inout
    & dinsect_dropping_mt(:,:)               , & ! inout
    & dinsect_slow_litter_mt(:,:,:)          , & ! inout
    & dinsect_defoliator_resp_c(:)           , & ! inout
    & dinsect_defoliator_resp_c13(:)         , & ! inout
    & dinsect_defoliator_resp_c14(:)         , & ! inout
    & mask=(defoliator_effects_triggered(:) > 0.0_wp)) ! optional(in)

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
  dsl4jsb_Real2D_onChunk      :: dinsect_occupy_frac           ! the fraction of snag occupyed area, which cause the waste of seed.
  dsl4jsb_Real2D_onChunk      :: defoliated_state              ! defoliation on or off
  dsl4jsb_Real2D_onChunk      :: is_insect_year                ! is insect disturbance year or not (1: true, 0: false)
  dsl4jsb_Real2D_onChunk      :: dinsect_year_since_start      ! the number of years since the start of insect disturbance
  dsl4jsb_Real2D_onChunk      :: dinsect_mort_dt               ! addtional mortality rate caused by bark beetle
  dsl4jsb_Real2D_onChunk      :: dinsect_defoli_dt             ! defoliation rate caused by defoliator  !dsl4jsb_Real2D_onChunk      :: heartwood_C_zero              ! the amount of heartwood C at the first timestep of bark beetle event
  dsl4jsb_Real2D_onChunk      :: dinsect_flux_frac             ! the fraction of disturbed fraction of target variables.
  dsl4jsb_Real2D_onChunk      :: dinsect_disturbable_frac      ! the fraction of remaining undisturbable tree resources.
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c     ! defoliator respiration C flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c13   ! defoliator respiration C13 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_defoliator_resp_c14   ! defoliator respiration C14 flux [mol/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c   ! the defence used labile C from bark beetle attack [mol C/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c13 ! the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
  dsl4jsb_Real2D_onChunk      :: dinsect_barkbeetle_labile_c14 ! the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
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
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_occupy_frac) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_mort_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_disturbable_frac) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoli_dt) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_flux_frac) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, defoliated_state) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, is_insect_year) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_year_since_start) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_defoliator_resp_c14) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c13) ! inout
  dsl4jsb_Get_var2D_onChunk(Q_DINSECT_, dinsect_barkbeetle_labile_c14) ! inout

#ifdef __QUINCY_STANDALONE__
      year = forcing_options%forcing_year
      doy = forcing_options%forcing_doy
      hour = forcing_options%forcing_hour
#else
      year = get_year(options%current_datetime)
      doy = INT(get_year_day(options%current_datetime))
      hour = get_secs_of_day(options%current_datetime) / 3600._wp ! seconds to hours
#endif

  ! (fake) read input data (daily), import from namelist, apply same data to every year
  dinsect_flux_frac(:) = dsl4jsb_Config(Q_DINSECT_)%dinsect_dev_bb_forcing_dist_frac(doy) / (one_day / dtime) ! distribute bb outbreak accross the whole day

  ! ======================================================================================================= !
  !> Go science
  !>

  ! ! If the insect is bark beetle, the damage fraction is calculated as the relative damage fraction based on the last year's remaining pool
  ! yearly_damage_frac = (yearly_damage_frac - yearly_damage_frac_base) / (1.0_wp - yearly_damage_frac_base)

  ! CALL derive_flux_frac_with_sigmoid_from_nml(&
  !   & dsl4jsb_Config(Q_DINSECT_)%dinsect_start_doy, &
  !   & dsl4jsb_Config(Q_DINSECT_)%dinsect_duration, &
  !   & year, &
  !   & doy, &
  !   & hour, &
  !   & yearly_damage_frac, &
  !   & dinsect_flux_frac(ic))

  !< -----------------------   1.1 get mortality_dt ---------------------------------
  dinsect_mort_dt(:) =  dinsect_flux_frac(:)
  !< ------------------ 1.2 Calculate effects (after bark beetle event happens) --------------------
  ! DO ic = 1, nc

  !   ! IF (dinsect_mort_dt(ic) > eps8) THEN
  !   !   barkbeetle_effects_triggered(ic) = barkbeetle_effects_triggered(ic) + 1.0_wp
  !   ! END IF

  !   IF ((ABS(barkbeetle_effects_triggered(ic)) - 1.0_wp) < eps4) THEN
  !     dinsect_disturbable_frac(ic) = 1.0_wp
  !   END IF

  ! END DO
  dinsect_disturbable_frac = 0.0_wp

  CALL calc_beetle_effects(&
         & nc                                     , & ! in
         & dtime                                  , & ! in
         & dinsect_mort_dt(:)                     , & ! in
         & veg_pool_mt(:,:,:)                     , & ! in
         & seed_bed_pool_mt(:,:)                  , & ! in
         & dsl4jsb_Config(Q_DINSECT_)%dinsect_online, & ! in
          !& heartwood_C_zero(:)                    , & ! in
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
         & dinsect_occupy_frac(:)                 , & ! inout, to remove
         & dinsect_disturbable_frac(:)            , & ! inout, to remove
         & mask=(dinsect_mort_dt(:) > 0.0_wp)) ! optional(in)

    WRITE(*,*) "dinsect_main:: update_dinsect_BB:: lctnumber, dtime(s), doy(d), one_day(?), ", &
            "dist_fract, oneday/dtime, dist_fract(doy)=", &
                lctlib%LctNumber, dtime, doy, one_day, dinsect_flux_frac(1), &
                one_day / dtime, &
                dsl4jsb_Config(Q_DINSECT_)%dinsect_dev_bb_forcing_dist_frac(doy), &
                
         
END SUBROUTINE update_dinsect_barkbeetle


#endif
END MODULE mo_q_dinsect_main
