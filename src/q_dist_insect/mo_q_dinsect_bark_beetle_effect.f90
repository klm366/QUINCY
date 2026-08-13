!> Calculate bark beetle effects on vegetation
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
!>#### Calculates bark beetle effects and associated fluxes
MODULE mo_q_dinsect_bark_beetle_effect
#ifndef __NO_JSBACH__

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_control,             ONLY: debug_on
    USE mo_exception,               ONLY: message
    USE mo_jsb_math_constants,      ONLY: one_day, one_year, eps8, eps4
    USE mo_lnd_bgcm_idx
    USE mo_veg_config_class,        ONLY: get_number_of_veg_compartments
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, &
     & DINSECT_BGCM_DEADSAPWOOD_ID, DINSECT_BGCM_DEADHEARTWOOD_ID, DINSECT_BGCM_DEADLEAF_ID, &
     & DINSECT_BGCM_VEG_EXTRACT_ID, DINSECT_BGCM_FAST_LITTER_ID, &
     & DINSECT_BGCM_SLOW_LITTER_ID ,&
     & DINSECT_BGCM_LARVA_ID, DINSECT_BGCM_DROPPING_ID
    USE mo_q_dinsect_constants,       ONLY: tau_dead_heartwood, tau_dead_sapwood, tau_dead_leaf

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: calc_beetle_effects

    CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_bark_beetle_effect'

  CONTAINS

    !> Calculate bark beetle effects for each timestep
    !>
    !> @brief Handles three main processes:
    !>   1. Vegetation mortality (veg -> dead pools)
    !>   2. Fast litter production (veg -> soil)
    !>   3. Slow litter production (dead pools -> soil)
    SUBROUTINE calc_beetle_effects(     &
      & nc                            , &  !< number of cells in chunk
      & dtime                         , &  !< timestep [s]
      & dinsect_mort_dt               , &  !< mortality rate [fraction per timestep * 1000]
      & veg_pool_mt                   , &  !< vegetation pool [mol/m2]
      & seed_bed_pool_mt              , &  !< bgcm: seed_bed pool [mol/m2]
      & online_flag                   , &  !< online calculation flag
      & dinsect_deadsapwood_mt        , &  !< dead sapwood pool [mol/m2]
      & dinsect_deadheartwood_mt      , &  !< dead heartwood pool [mol/m2]
      & dinsect_deadleaf_mt           , &  !< dead leaf pool [mol/m2]
      & dinsect_veg_extract_mt        , &  !< vegetation extraction flux [mol/m2/timestep]
      & dinsect_slow_litter_mt        , &  !< slow litter production flux [mol/m2/timestep]
      & dinsect_fast_litter_mt        , &  !< fast litter production flux [mol/m2/timestep]
      & dinsect_seedb_fastl_mt , &  !< fast litter seed bed production flux [mol/m2/timestep]
      & dinsect_barkbeetle_labile_c   , &  !< the defence used labile C from bark beetle attack [mol C/m2/timestep]
      & dinsect_barkbeetle_labile_c13 , &  !< the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
      & dinsect_barkbeetle_labile_c14 , &  !< the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
      & mask)                              !< Gridcells to compute (default: all gridcells)
      ! ----------------------------------------------------------------------------------------------------- !
      ! Input/Output Variables
      INTEGER,  INTENT(in)    :: nc                            !< number of cells in chunk
      REAL(wp), INTENT(in)    :: dtime                         !< timestep [s]
      REAL(wp), INTENT(in)    :: dinsect_mort_dt(:)            !< mortality rate [fraction per timestep * 1000]
      !REAL(wp), INTENT(in)    :: heartwood_C_zero(:)           !< heartwood C at the beginning of the bark beetle outbreak [mol/m2]
      REAL(wp), INTENT(in)    :: veg_pool_mt(:,:,:)            !< vegetation pool [mol/m2]
      REAL(wp), INTENT(in)    :: seed_bed_pool_mt(:,:)         !< bgcm: seed_bed pool [mol/m2]
      LOGICAL,  INTENT(in)    :: online_flag                   !< online calculation flag
      REAL(wp), INTENT(inout) :: dinsect_deadsapwood_mt(:,:)   !< dead sapwood pool [mol/m2]
      REAL(wp), INTENT(inout) :: dinsect_deadheartwood_mt(:,:) !< dead heartwood pool [mol/m2]
      REAL(wp), INTENT(inout) :: dinsect_deadleaf_mt(:,:)      !< dead leaf pool [mol/m2]
      REAL(wp), INTENT(inout) :: dinsect_veg_extract_mt(:,:,:) !< extraction flux [mol/m2/timestep]
      REAL(wp), INTENT(inout) :: dinsect_fast_litter_mt(:,:,:) !< fast litter flux [mol/m2/timestep]
      REAL(wp), INTENT(inout) :: dinsect_seedb_fastl_mt(:,:) !< fast litter seed bed flux [mol/m2/timestep]
      REAL(wp), INTENT(inout) :: dinsect_slow_litter_mt(:,:,:) !< slow litter flux [mol/m2/timestep]
      REAL(wp), INTENT(inout) :: dinsect_barkbeetle_labile_c(:)   !< the defence used labile C from bark beetle attack [mol C/m2/timestep]
      REAL(wp), INTENT(inout) :: dinsect_barkbeetle_labile_c13(:) !< the defence used labile C13 from bark beetle attack [mol C13/m2/timestep]
      REAL(wp), INTENT(inout) :: dinsect_barkbeetle_labile_c14(:) !< the defence used labile C14 from bark beetle attack [mol C14/m2/timestep]
      LOGICAL, intent(in), optional :: mask(:)                    !< Gridcells to compute (default: all gridcells)
      ! ----------------------------------------------------------------------------------------------------- !
      ! Local Variables
      INTEGER  :: ic                                  !< loop counter for cells
      REAL(wp), DIMENSION(nc) :: mortality_dt         !< mortality rate [fraction/timestep]
      ! REAL(wp) :: decay_factor_for_BB              ! the decay factor for calculating bark beetle snag effects
      CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_beetle_effects'
      LOGICAL :: mask_all(nc)

      ! ----------------------------------------------------------------------------------------------------- !

      mask_all = .TRUE.
      IF (PRESENT(mask)) THEN
        mask_all = mask
      END IF

      ! >
      !> 1.0 Calculate mortality effects and dead pool increase
      !>
      DO ic = 1, nc
        IF (.NOT. mask_all(ic)) CYCLE

        mortality_dt(ic) = dinsect_mort_dt(ic) ! todo: no need for that, use original variable

        !> 1.1 Update dead pools with killed vegetation
        dinsect_deadleaf_mt(:, ic) = dinsect_deadleaf_mt(:, ic) + &
          & veg_pool_mt(ix_leaf, :, ic) * mortality_dt(ic)
        dinsect_deadsapwood_mt(:, ic) = dinsect_deadsapwood_mt(:, ic) + &
          & veg_pool_mt(ix_sap_wood, :, ic) * mortality_dt(ic)
        dinsect_deadheartwood_mt(:, ic) = dinsect_deadheartwood_mt(:, ic) + &
          & veg_pool_mt(ix_heart_wood, :, ic) * mortality_dt(ic)

        IF (online_flag) THEN

          !> 1.2 Update extraction fluxes for standing dead material
          dinsect_veg_extract_mt(ix_leaf, :, ic) = dinsect_veg_extract_mt(ix_leaf, :, ic) + &
          & veg_pool_mt(ix_leaf, :, ic) * mortality_dt(ic)
          dinsect_veg_extract_mt(ix_sap_wood, :, ic) = dinsect_veg_extract_mt(ix_sap_wood, :, ic) + &
          & veg_pool_mt(ix_sap_wood, :, ic) * mortality_dt(ic)
          dinsect_veg_extract_mt(ix_heart_wood, :, ic) = dinsect_veg_extract_mt(ix_heart_wood, :, ic) + &
          & veg_pool_mt(ix_heart_wood, :, ic) * mortality_dt(ic)

          !> 1.3 Update labile pool losses (defense compounds)
          dinsect_veg_extract_mt(ix_labile, ixC, ic) = dinsect_veg_extract_mt(ix_labile, ixC, ic) + &
          & veg_pool_mt(ix_labile, ixC, ic) * mortality_dt(ic)
          dinsect_veg_extract_mt(ix_labile, ixC13, ic) = dinsect_veg_extract_mt(ix_labile, ixC13, ic) + &
          & veg_pool_mt(ix_labile, ixC13, ic) * mortality_dt(ic)
          dinsect_veg_extract_mt(ix_labile, ixC14, ic) = dinsect_veg_extract_mt(ix_labile, ixC14, ic) + &
          & veg_pool_mt(ix_labile, ixC14, ic) * mortality_dt(ic)

          !> 1.4 Update fast litter production (nutrients from labile pool)
          dinsect_fast_litter_mt(ix_labile, ixN, ic) = dinsect_fast_litter_mt(ix_labile, ixN, ic) + &
          & veg_pool_mt(ix_labile, ixN, ic) * mortality_dt(ic)
          dinsect_fast_litter_mt(ix_labile, ixN15, ic) = dinsect_fast_litter_mt(ix_labile, ixN15, ic) + &
          & veg_pool_mt(ix_labile, ixN15, ic) * mortality_dt(ic)
          dinsect_fast_litter_mt(ix_labile, ixP, ic) = dinsect_fast_litter_mt(ix_labile, ixP, ic) + &
          & veg_pool_mt(ix_labile, ixP, ic) * mortality_dt(ic)

          !> 1.5 Update other vegetation pool transfers to fast litter
          dinsect_fast_litter_mt(ix_fine_root, :, ic) = dinsect_fast_litter_mt(ix_fine_root, :, ic) + &
          & veg_pool_mt(ix_fine_root, :, ic) * mortality_dt(ic)
          dinsect_fast_litter_mt(ix_coarse_root, :, ic) = dinsect_fast_litter_mt(ix_coarse_root, :, ic) + &
          & veg_pool_mt(ix_coarse_root, :, ic) * mortality_dt(ic)
          dinsect_fast_litter_mt(ix_reserve, :, ic) = dinsect_fast_litter_mt(ix_reserve, :, ic) + &
          & veg_pool_mt(ix_reserve, :, ic) * mortality_dt(ic)
          dinsect_fast_litter_mt(ix_fruit, :, ic) = dinsect_fast_litter_mt(ix_fruit, :, ic) + &
          & veg_pool_mt(ix_fruit, :, ic) * mortality_dt(ic)
          dinsect_seedb_fastl_mt(:, ic) = dinsect_seedb_fastl_mt(:, ic) + &
          & seed_bed_pool_mt(:, ic) * mortality_dt(ic)

          !> 1.6 Update respiration fluxes. This part of labile C will be accounted for in veg_update_pools
          dinsect_barkbeetle_labile_c(ic) = veg_pool_mt(ix_labile, ixC, ic) * mortality_dt(ic)
          dinsect_barkbeetle_labile_c13(ic) = veg_pool_mt(ix_labile, ixC13, ic) * mortality_dt(ic)
          dinsect_barkbeetle_labile_c14(ic) = veg_pool_mt(ix_labile, ixC14, ic) * mortality_dt(ic)

          END IF


        !> 2.0 Calculate decomposition of dead pools
        !>
        !> 2.1 Dead leaf decomposition
        IF (MAXVAL(dinsect_deadleaf_mt(:,ic)) > eps8) THEN
          IF (online_flag) THEN
            dinsect_slow_litter_mt(ix_leaf, :, ic) = dinsect_slow_litter_mt(ix_leaf, :, ic) + &
              & dinsect_deadleaf_mt(:, ic) * &
              & (1.0_wp - EXP(-1.0_wp / (tau_dead_leaf * one_year * one_day / dtime)))
          END IF
          dinsect_deadleaf_mt(:, ic) = dinsect_deadleaf_mt(:, ic) - dinsect_deadleaf_mt(:, ic)  * &
          & (1.0_wp - EXP(-1.0_wp / (tau_dead_leaf * one_year * one_day / dtime)))
        END IF

        !> 2.2 Dead sapwood decomposition
        IF (MAXVAL(dinsect_deadsapwood_mt(:,ic)) > eps8) THEN
          IF (online_flag) THEN
            dinsect_slow_litter_mt(ix_sap_wood, :, ic) = dinsect_slow_litter_mt(ix_sap_wood, :, ic) + &
              & dinsect_deadsapwood_mt(:, ic) * &
              & (1.0_wp - EXP(-1.0_wp / (tau_dead_sapwood * one_year * one_day / dtime)))
          END IF
          dinsect_deadsapwood_mt(:, ic) = dinsect_deadsapwood_mt(:, ic) - dinsect_deadsapwood_mt(:, ic) * &
          & (1.0_wp - EXP(-1.0_wp / (tau_dead_sapwood * one_year * one_day / dtime)))
        END IF

        !> 2.3 Dead heartwood decomposition
        IF (MAXVAL(dinsect_deadheartwood_mt(:,ic)) > eps8) THEN
          IF (online_flag) THEN
            dinsect_slow_litter_mt(ix_heart_wood, :, ic) = dinsect_slow_litter_mt(ix_heart_wood, :, ic) + &
              & dinsect_deadheartwood_mt(:, ic) * &
              & (1.0_wp - EXP(-1.0_wp / (tau_dead_heartwood * one_year * one_day / dtime)))
          END IF
          dinsect_deadheartwood_mt(:, ic) = dinsect_deadheartwood_mt(:, ic) - dinsect_deadheartwood_mt(:, ic)&
          &  * (1.0_wp - EXP(-1.0_wp / (tau_dead_heartwood * one_year * one_day / dtime)))
        END IF
      END DO

    END SUBROUTINE calc_beetle_effects

#endif
END MODULE mo_q_dinsect_bark_beetle_effect
