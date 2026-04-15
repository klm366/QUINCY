!> Calculate defoliator effects on vegetation
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
!>#### Calculates insect defoliation effects and associated fluxes
MODULE mo_q_dinsect_defoliator_effect
#ifndef __NO_JSBACH__

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message
  USE mo_jsb_math_constants,      ONLY: one_day, eps12, eps4, eps12
  USE mo_lnd_bgcm_idx
  USE mo_veg_config_class,        ONLY: get_number_of_veg_compartments
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_GROWTH_ID, &
   & DINSECT_BGCM_DEADSAPWOOD_ID, DINSECT_BGCM_DEADHEARTWOOD_ID, DINSECT_BGCM_DEADLEAF_ID, &
   & DINSECT_BGCM_VEG_EXTRACT_ID, DINSECT_BGCM_FAST_LITTER_ID, &
   & DINSECT_BGCM_SLOW_LITTER_ID ,&
   & DINSECT_BGCM_LARVA_ID, DINSECT_BGCM_DROPPING_ID
  USE mo_q_dinsect_constants

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: calc_defoliator_effects

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_dinsect_defoliator_effect'

CONTAINS

!> Calculate the effects of defoliator for each timestep
!>
!> @brief Handles three main processes:
!>   1. Vegetation extraction (veg_leaf -> larva)
!>   2. Dropping production (larva -> soil)
!>   3. Litter production (dead pool -> soil)
SUBROUTINE calc_defoliator_effects( &
  & nc                          , &  !< number of cells in chunk
  & dtime                       , &  !< timestep [s]
  & doy                         , &  !< day of year
  & dinsect_defoli_dt           , &  !< defoliation rate [fraction per timestep  * 1000]
  & veg_pool_mt                 , &  !< vegetation pool [mol/m2]
  & flag_online                 , &  !< online calculation flag
  & dinsect_larva_mt            , &  !< larva pool [mol/m2]
  & dinsect_veg_extract_mt      , &  !< vegetation extraction flux [mol/m2/timestep]
  & dinsect_dropping_mt         , &  !< litter production flux [mol/m2/timestep]
  & dinsect_slow_litter_mt      , &  !< litter production flux [mol/m2/timestep]
  & dinsect_defoliator_resp_c   , &  !< defoliator respiration C flux [mol/m2/timestep]
  & dinsect_defoliator_resp_c13 , &  !< defoliator respiration C13 flux [mol/m2/timestep]
  & dinsect_defoliator_resp_c14 , &  !< defoliator respiration C14 flux [mol/m2/timestep]
  & mask                        )    !< Enabled gridcells mask

    ! ----------------------------------------------------------------------------------------------------- !
    ! Input/Output Variables
    INTEGER,  INTENT(in)    :: nc                             !< number of cells in chunk
    REAL(wp), INTENT(in)    :: dtime                          !< timestep [s]
    INTEGER,  INTENT(in)    :: doy                            !< day of year
    REAL(wp), INTENT(in)    :: dinsect_defoli_dt(:)           !< defoliation rate [fraction per timestep * 1000]
    REAL(wp), INTENT(in)    :: veg_pool_mt(:,:,:)             !< bgcm veg_pool [mol/m2]
    LOGICAL,  INTENT(in)    :: flag_online                    !< online calculation flag
    REAL(wp), INTENT(inout) :: dinsect_defoliator_resp_c(:)   !< defoliator respiration C flux [mol/m2/timestep]
    REAL(wp), INTENT(inout) :: dinsect_defoliator_resp_c13(:) !< defoliator respiration C13 flux [mol/m2/timestep]
    REAL(wp), INTENT(inout) :: dinsect_defoliator_resp_c14(:) !< defoliator respiration C14 flux [mol/m2/timestep]
    REAL(wp), INTENT(inout) :: dinsect_larva_mt(:,:)          !< larva pool [mol/m2]
    REAL(wp), INTENT(inout) :: dinsect_veg_extract_mt(:,:,:)  !< extraction flux [mol/m2/timestep]
    REAL(wp), INTENT(inout) :: dinsect_dropping_mt(:,:)       !< dropping flux [mol/m2/timestep]
    REAL(wp), INTENT(inout) :: dinsect_slow_litter_mt(:,:,:)  !< litter flux [mol/m2/timestep]
    LOGICAL,  INTENT(in), OPTIONAL :: mask(:)                 !< Enabled gridcells mask (default: all enabled)

    ! ----------------------------------------------------------------------------------------------------- !
    ! Local Variables - Loop Indices
    INTEGER  :: ic                                           !< loop counter for cells
    REAL(wp), DIMENSION(nc) :: defoliation_dt                !< defoliation rate [fraction per timestep]

    ! ----------------------------------------------------------------------------------------------------- !
    ! Local Variables - Defoliation Process
    REAL(wp), DIMENSION(nc) :: defoli_c                      !< defoliated carbon [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: defoli_n                      !< defoliated nitrogen [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: defoli_p                      !< defoliated phosphorus [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: defoli_c13                    !< defoliated carbon 13 [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: defoli_c14                    !< defoliated carbon 14 [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: defoli_n15                    !< defoliated nitrogen 15 [mol/m2/timestep]

    ! ----------------------------------------------------------------------------------------------------- !
    ! Local Variables - Larva Consumption
    REAL(wp), DIMENSION(nc) :: eaten_c                       !< C consumed by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: eaten_n                       !< N consumed by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: eaten_p                       !< P consumed by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: eaten_c13                    !< C13 consumed by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: eaten_c14                    !< C14 consumed by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: eaten_n15                    !< N15 consumed by larva [mol/m2/timestep]

    ! ----------------------------------------------------------------------------------------------------- !
    ! Local Variables - Larva Growth
    REAL(wp), DIMENSION(nc) :: assim_c                       !< C assimilated by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: assim_c13                    !< C13 assimilated by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: assim_c14                    !< C14 assimilated by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: assim_n                      !< N assimilated by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: assim_n15                    !< N15 assimilated by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: assim_p                      !< P assimilated by larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: larva_gain_c                  !< net C gain in larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: larva_gain_n                  !< net N gain in larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: larva_gain_p                  !< net P gain in larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: larva_gain_c13                !< net C13 gain in larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: larva_gain_c14                !< net C14 gain in larva [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: larva_gain_n15                !< net N15 gain in larva [mol/m2/timestep]

    ! ----------------------------------------------------------------------------------------------------- !
    ! Local Variables - Respiration
    REAL(wp), DIMENSION(nc) :: growth_resp_c                 !< growth respiration [mol C/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_resp_c13               !< growth respiration C13 [mol C13/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_resp_c14               !< growth respiration C14 [mol C14/m2/timestep]
    REAL(wp), DIMENSION(nc) :: maint_resp_c                  !< maintenance respiration [mol C/m2/timestep]
    REAL(wp), DIMENSION(nc) :: maint_resp_c13                !< maintenance respiration C13 [mol C13/m2/timestep]
    REAL(wp), DIMENSION(nc) :: maint_resp_c14                !< maintenance respiration C14 [mol C14/m2/timestep]

    ! ----------------------------------------------------------------------------------------------------- !
    ! Local Variables - Mortality and Excretion
    REAL(wp), DIMENSION(nc) :: die_c                         !< C loss from mortality [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: die_c13                      !< C13 loss from mortality [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: die_c14                      !< C14 loss from mortality [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: die_n                         !< N loss from mortality [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: die_n15                      !< N15 loss from mortality [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: die_p                         !< P loss from mortality [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: excrete_c                     !< C excretion [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: excrete_c13                  !< C13 excretion [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: excrete_c14                  !< C14 excretion [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: excrete_n                     !< N excretion [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: excrete_p                     !< P excretion [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: excrete_n15                   !< N15 excretion [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: feces_c                      !< C in feces [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: feces_c13                   !< C13 in feces [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: feces_c14                   !< C14 in feces [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: feces_n                     !< N in feces [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: feces_p                     !< P in feces [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: feces_n15                   !< N15 in feces [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: c_avail                     !< C available for growth [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: p_from_c                    !< P from C [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: p_from_n                    !< P from N [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: p_from_p                    !< P from P [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_c                    !< C from growth [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_n                    !< N from growth [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_p                    !< P from growth [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_c13                  !< C13 from growth [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_c14                  !< C14 from growth [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: growth_n15                  !< N15 from growth [mol/m2/timestep]


    ! ----------------------------------------------------------------------------------------------------- !
    ! Local Variables - Litter Production
    REAL(wp), DIMENSION(nc) :: leaf_fragment_c               !< C in leaf fragments [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: leaf_fragment_c13            !< C13 in leaf fragments [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: leaf_fragment_c14            !< C14 in leaf fragments [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: leaf_fragment_n               !< N in leaf fragments [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: leaf_fragment_n15            !< N15 in leaf fragments [mol/m2/timestep]
    REAL(wp), DIMENSION(nc) :: leaf_fragment_p               !< P in leaf fragments [mol/m2/timestep]

    ! ----------------------------------------------------------------------------------------------------- !
    ! Helper variables
    REAL(wp), DIMENSION(nc) :: target_p                      !< target P [mol/m2/timestep]

    ! ----------------------------------------------------------------------------------------------------- !
    LOGICAL, DIMENSION(nc) :: mask_all                       !< Enabled gridcells mask
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_defoliator_effects'

    ! ----------------------------------------------------------------------------------------------------- !
    !< Go science
    !<
    ! ----------------------------------------------------------------------------------------------------- !
    mask_all = .TRUE.
    IF (PRESENT(mask)) THEN
      mask_all = mask
    END IF

    defoliation_dt(:) = dinsect_defoli_dt(:)

    !> 1.0 Calculate defoliation fluxes
    !>
    DO ic = 1, nc
      IF (.NOT. mask_all(ic)) CYCLE !< skip gridcell

      !>
      !> 1.1 Calculate defoliation amounts
      !>
      defoli_c(ic) = veg_pool_mt(ix_leaf, ixC, ic) * defoliation_dt(ic)
      defoli_n(ic) = veg_pool_mt(ix_leaf, ixN, ic) * defoliation_dt(ic)
      defoli_p(ic) = veg_pool_mt(ix_leaf, ixP, ic) * defoliation_dt(ic)
      defoli_c13(ic) = veg_pool_mt(ix_leaf, ixC13, ic) * defoliation_dt(ic)
      defoli_c14(ic) = veg_pool_mt(ix_leaf, ixC14, ic) * defoliation_dt(ic)
      defoli_n15(ic) = veg_pool_mt(ix_leaf, ixN15, ic) * defoliation_dt(ic)

      !> 1.2 Calculate the amount extraction fluxes from vegetation
      !>
      dinsect_veg_extract_mt(ix_leaf, ixC, ic) = dinsect_veg_extract_mt(ix_leaf, ixC, ic) + defoli_c(ic)
      dinsect_veg_extract_mt(ix_leaf, ixN, ic) = dinsect_veg_extract_mt(ix_leaf, ixN, ic) + defoli_n(ic)
      dinsect_veg_extract_mt(ix_leaf, ixP, ic) = dinsect_veg_extract_mt(ix_leaf, ixP, ic) + defoli_p(ic)
      dinsect_veg_extract_mt(ix_leaf, ixC13, ic) = dinsect_veg_extract_mt(ix_leaf, ixC13, ic) + defoli_c13(ic)
      dinsect_veg_extract_mt(ix_leaf, ixC14, ic) = dinsect_veg_extract_mt(ix_leaf, ixC14, ic) + defoli_c14(ic)
      dinsect_veg_extract_mt(ix_leaf, ixN15, ic) = dinsect_veg_extract_mt(ix_leaf, ixN15, ic) + defoli_n15(ic)

      !> 2.0 Calculate larva consumption and growth
      !>
      !> 2.1 Calculate eaten amounts with efficiency
      eaten_c(ic) = defoli_c(ic) * defoli_eat_efficiency
      eaten_n(ic) = defoli_n(ic) * defoli_eat_efficiency
      eaten_p(ic) = defoli_p(ic) * defoli_eat_efficiency
      eaten_c13(ic) = defoli_c13(ic) * defoli_eat_efficiency
      eaten_c14(ic) = defoli_c14(ic) * defoli_eat_efficiency
      eaten_n15(ic) = defoli_n15(ic) * defoli_eat_efficiency

      !> 2.2 Calculate leaf fragments (non-digested parts)
      leaf_fragment_c(ic) = defoli_c(ic) - eaten_c(ic)
      leaf_fragment_n(ic) = defoli_n(ic) - eaten_n(ic)
      leaf_fragment_p(ic) = defoli_p(ic) - eaten_p(ic)
      leaf_fragment_c13(ic) = defoli_c13(ic) - eaten_c13(ic)
      leaf_fragment_c14(ic) = defoli_c14(ic) - eaten_c14(ic)
      leaf_fragment_n15(ic) = defoli_n15(ic) - eaten_n15(ic)

      !> 2.3 Precaculate assimilation
      assim_c(ic) = eaten_c(ic) * defoli_assim_frac_C
      assim_c13(ic) = eaten_c13(ic) * defoli_assim_frac_C
      assim_c14(ic) = eaten_c14(ic) * defoli_assim_frac_C

      assim_n(ic) = eaten_n(ic) * defoli_assim_frac_N  ! now I just set defoli_assim_frac_N = 1.0
      assim_n15(ic) = eaten_n15(ic) * defoli_assim_frac_N  ! now I just set defoli_assim_frac_N15 = 1.0
      assim_p(ic) = eaten_p(ic) * defoli_assim_frac_P  ! now I just set defoli_assim_frac_P = 1.0

      !> 2.4 Calculate unassimilated feces
      feces_c(ic) = eaten_c(ic) - assim_c(ic)
      feces_c13(ic) = eaten_c13(ic) - assim_c13(ic)
      feces_c14(ic) = eaten_c14(ic) - assim_c14(ic)
      feces_n(ic) = eaten_n(ic) - assim_n(ic)
      feces_n15(ic) = eaten_n15(ic) - assim_n15(ic)
      feces_p(ic) = eaten_p(ic) - assim_p(ic)

      !> 2.5 Calculate growth respiration
      growth_resp_c(ic) = assim_c(ic) * defoli_growth_resp_frac
      growth_resp_c13(ic) = assim_c13(ic) * defoli_growth_resp_frac
      growth_resp_c14(ic) = assim_c14(ic) * defoli_growth_resp_frac

      c_avail(ic) = MAX(assim_c(ic) - growth_resp_c(ic), 0.0_wp)

      !> 2.6 Calculate mass gained

      p_from_c(ic) =  (c_avail(ic) * C_mole_mass) / defoli_CP_ratio / P_mole_mass
      p_from_n(ic) =  (assim_n(ic) * N_mole_mass) / defoli_NP_ratio / P_mole_mass
      p_from_p(ic) =  assim_p(ic)

      growth_p(ic) = MIN( p_from_c(ic), p_from_n(ic), p_from_p(ic) )
      growth_c(ic) = (growth_p(ic) * P_mole_mass) * defoli_CP_ratio / C_mole_mass
      growth_n(ic) = (growth_p(ic) * P_mole_mass) * defoli_NP_ratio / N_mole_mass

      IF (growth_c(ic) > eps12) THEN
        larva_gain_c13(ic) = growth_c(ic) * ( assim_c13(ic) / MAX(assim_c(ic), eps12) )
        larva_gain_c14(ic) = growth_c(ic) * ( assim_c14(ic) / MAX(assim_c(ic), eps12) )
      ELSE
        larva_gain_c13(ic) = 0.0_wp
        larva_gain_c14(ic) = 0.0_wp
      END IF
      IF (growth_n(ic) > eps12) THEN
        larva_gain_n15(ic) = growth_n(ic) * ( assim_n15(ic) / MAX(assim_n(ic), eps12) )
      ELSE
        larva_gain_n15(ic) = 0.0_wp
      END IF

      larva_gain_c(ic) = growth_c(ic)
      larva_gain_n(ic) = growth_n(ic)
      larva_gain_p(ic) = growth_p(ic)


      !> 3.0 Calculate mortality and maintenance respiration
      !>
      !> 3.1 Calculate mortality fluxes
      !IF (MAXVAL(dinsect_larva_mt(:,ic)) > eps12) THEN
        die_c(ic) = 1.0_wp / (tau_insect * one_day / dtime) * dinsect_larva_mt(ixC,ic)!    MAX(dinsect_larva_mt(ixC,ic) * (1.0_wp - EXP(-1.0_wp / (tau_insect * one_day / dtime))), 0.0_wp)
        die_n(ic) = 1.0_wp / (tau_insect * one_day / dtime) * dinsect_larva_mt(ixN,ic)!    MAX(dinsect_larva_mt(ixN,ic) * (1.0_wp - EXP(-1.0_wp / (tau_insect * one_day / dtime))), 0.0_wp)
        die_p(ic) = 1.0_wp / (tau_insect * one_day / dtime) * dinsect_larva_mt(ixP,ic)!    MAX(dinsect_larva_mt(ixP,ic) * (1.0_wp - EXP(-1.0_wp / (tau_insect * one_day / dtime))), 0.0_wp)
        die_c13(ic) = 1.0_wp / (tau_insect * one_day / dtime) * dinsect_larva_mt(ixC13,ic)!    MAX(dinsect_larva_mt(ixC13,ic) * (1.0_wp - EXP(-1.0_wp / (tau_insect * one_day / dtime))), 0.0_wp)
        die_c14(ic) = 1.0_wp / (tau_insect * one_day / dtime) * dinsect_larva_mt(ixC14,ic)!    MAX(dinsect_larva_mt(ixC14,ic) * (1.0_wp - EXP(-1.0_wp / (tau_insect * one_day / dtime))), 0.0_wp)
        die_n15(ic) = 1.0_wp / (tau_insect * one_day / dtime) * dinsect_larva_mt(ixN15,ic)!    MAX(dinsect_larva_mt(ixN15,ic) * (1.0_wp - EXP(-1.0_wp / (tau_insect * one_day / dtime))), 0.0_wp)
      !END IF

      !> 3.2 Calculate maintenance respiration
      maint_resp_c(ic) = (dinsect_larva_mt(ixC,ic) + larva_gain_c(ic) - die_c(ic)) * defoli_maint_resp_frac
      maint_resp_c13(ic) = (dinsect_larva_mt(ixC13,ic) + larva_gain_c13(ic) - die_c13(ic)) * defoli_maint_resp_frac
      maint_resp_c14(ic) = (dinsect_larva_mt(ixC14,ic) + larva_gain_c14(ic) - die_c14(ic)) * defoli_maint_resp_frac

      !> 3.3 Update larva pool
      dinsect_larva_mt(ixC,ic) = dinsect_larva_mt(ixC,ic) + larva_gain_c(ic) - die_c(ic) - maint_resp_c(ic)
      dinsect_larva_mt(ixN,ic) = dinsect_larva_mt(ixN,ic) + larva_gain_n(ic) - die_n(ic)
      dinsect_larva_mt(ixP,ic) = dinsect_larva_mt(ixP,ic) + larva_gain_p(ic) - die_p(ic)
      dinsect_larva_mt(ixC13,ic) = dinsect_larva_mt(ixC13,ic) + larva_gain_c13(ic) - die_c13(ic) - maint_resp_c13(ic)
      dinsect_larva_mt(ixC14,ic) = dinsect_larva_mt(ixC14,ic) + larva_gain_c14(ic) - die_c14(ic) - maint_resp_c14(ic)
      dinsect_larva_mt(ixN15,ic) = dinsect_larva_mt(ixN15,ic) + larva_gain_n15(ic) - die_n15(ic)

      !> 4.0 Calculate excretion and litter production
      !>
      !> 4.1 Calculate excretion (surplus)
      excrete_c(ic) = MAX(c_avail(ic) - larva_gain_c(ic), 0.0_wp)
      excrete_c13(ic) = MAX(assim_c13(ic) - larva_gain_c13(ic) - growth_resp_c13(ic), 0.0_wp)
      excrete_c14(ic) = MAX(assim_c14(ic) - larva_gain_c14(ic) - growth_resp_c14(ic), 0.0_wp)
      excrete_n(ic) = MAX(assim_n(ic) - larva_gain_n(ic), 0.0_wp)
      excrete_p(ic) = MAX(assim_p(ic) - larva_gain_p(ic), 0.0_wp)
      excrete_n15(ic) = MAX(assim_n15(ic) - larva_gain_n15(ic), 0.0_wp)

      !> 4.2 Update dropping pool (excretion + mortality)
      dinsect_dropping_mt(ixC, ic) = feces_c(ic) + excrete_c(ic) + die_c(ic)
      dinsect_dropping_mt(ixC13, ic) = feces_c13(ic) + excrete_c13(ic) + die_c13(ic)
      dinsect_dropping_mt(ixC14, ic) = feces_c14(ic) + excrete_c14(ic) + die_c14(ic)
      dinsect_dropping_mt(ixN, ic) = feces_n(ic) + excrete_n(ic) + die_n(ic)
      dinsect_dropping_mt(ixP, ic) = feces_p(ic) + excrete_p(ic) + die_p(ic)
      dinsect_dropping_mt(ixN15, ic) = feces_n15(ic) + excrete_n15(ic) + die_n15(ic)

      !> 4.3 Update litter pool with leaf fragments
      dinsect_slow_litter_mt(ix_leaf, ixC, ic) = dinsect_slow_litter_mt(ix_leaf, ixC, ic) + leaf_fragment_c(ic)
      dinsect_slow_litter_mt(ix_leaf, ixN, ic) = dinsect_slow_litter_mt(ix_leaf, ixN, ic) + leaf_fragment_n(ic)
      dinsect_slow_litter_mt(ix_leaf, ixP, ic) = dinsect_slow_litter_mt(ix_leaf, ixP, ic) + leaf_fragment_p(ic)
      dinsect_slow_litter_mt(ix_leaf, ixC13, ic) = dinsect_slow_litter_mt(ix_leaf, ixC13, ic) + leaf_fragment_c13(ic)
      dinsect_slow_litter_mt(ix_leaf, ixC14, ic) = dinsect_slow_litter_mt(ix_leaf, ixC14, ic) + leaf_fragment_c14(ic)
      dinsect_slow_litter_mt(ix_leaf, ixN15, ic) = dinsect_slow_litter_mt(ix_leaf, ixN15, ic) + leaf_fragment_n15(ic)

      !> 5.0 calculate respiration
      dinsect_defoliator_resp_c(ic) = growth_resp_c(ic) + maint_resp_c(ic)
      dinsect_defoliator_resp_c13(ic) = growth_resp_c13(ic) + maint_resp_c13(ic)
      dinsect_defoliator_resp_c14(ic) = growth_resp_c14(ic) + maint_resp_c14(ic)

      ! print*, 'larvae before', dinsect_larva_mt(ixC,ic)
      ! IF (doy > 246) THEN
      !     print*, 'droppings before', dinsect_dropping_mt(ixC,ic)
      !     ! move larvae to adult or to droppings with an emergence yield

      !     dinsect_dropping_mt(ixC,ic) = dinsect_dropping_mt(ixC,ic) + dinsect_larva_mt(ixC,ic)
      !     print*, 'droppings after', dinsect_dropping_mt(ixC,ic)
      !     dinsect_dropping_mt(ixC13,ic) = dinsect_dropping_mt(ixC13,ic) + dinsect_larva_mt(ixC13,ic)
      !     dinsect_dropping_mt(ixC14,ic) = dinsect_dropping_mt(ixC14,ic) + dinsect_larva_mt(ixC14,ic)
      !     dinsect_dropping_mt(ixN,ic) = dinsect_dropping_mt(ixN,ic) + dinsect_larva_mt(ixN,ic)
      !     dinsect_dropping_mt(ixP,ic) = dinsect_dropping_mt(ixP,ic) + dinsect_larva_mt(ixP,ic)
      !     dinsect_dropping_mt(ixN15,ic) = dinsect_dropping_mt(ixN15,ic) + dinsect_larva_mt(ixN15,ic)
      !     ! ... same for N,P and isotopes
      !     dinsect_larva_mt(ixC,ic) = 0.0_wp     ! or keep a pupae/adult pool instead
      !     dinsect_larva_mt(ixC13,ic) = 0.0_wp
      !     dinsect_larva_mt(ixC14,ic) = 0.0_wp
      !     dinsect_larva_mt(ixN,ic) = 0.0_wp
      !     dinsect_larva_mt(ixP,ic) = 0.0_wp
      !     dinsect_larva_mt(ixN15,ic) = 0.0_wp
      ! END IF
      ! print*, 'larvae out', dinsect_larva_mt(ixC,ic)

    END DO

    ! ----------------------------------------------------------------------------------------------------- !
  END SUBROUTINE calc_defoliator_effects


#endif
END MODULE mo_q_dinsect_defoliator_effect