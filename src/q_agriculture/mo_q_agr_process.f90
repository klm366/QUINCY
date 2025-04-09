!> QUINCY agriculture process calculations
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### routines for agriculture processes
!>
MODULE mo_q_agr_process
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_math_constants,  ONLY: eps8

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calc_crop_phenology, calc_planting_flux, set_crop_n_fixation_status
  PUBLIC :: update_crop_growth_phase, calc_crop_allocation_factors, calc_fertiliser_application

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_agr_process'

  CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate current gdd and its long-term average gdd_mavg for crops,
  !!  and determine inital values of growing demand for N and P at beginning of growing season
  !!  Also initialises crop growth phase to 0 at beginning of growing season
  !! This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!  temperature thresholds for gdd vary with crop type
  !!  deviations from CLM5:
  !!    - growing season is additionally constrained by soil moisture (same for all crop types)
  !!    - frost is assessed as frozen soil temperature, no other minimum temperatures are used
  !!    - weekly instead of 10day temperatures are used
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_crop_phenology(  &
    & nc                              , &
    & dtime                           , &
    & beta_soil_flush                 , &
    & gdd_req_max                     , &
    & beta_soil_senescence            , &
    & t_air_senescence                , &
    & cn_leaf                         , &
    & np_leaf                         , &
    & crop_type_index                 , &
    & t_air                           , &
    & t_air_week_mavg                 , &
    & t_air_month_mavg                , &
    & beta_soil_gs_tphen_mavg         , &
    & t_air_tphen_mavg                , &
    & t_soil_srf_tphen_mavg           , &
    & gpp_tlabile_mavg                , &
    & maint_respiration_tlabile_mavg  , &
    & growing_season                  , &
    & growth_req_n_tlabile_mavg       , &
    & growth_req_p_tlabile_mavg       , &
    & gdd                             , &
    & gdd_mavg                        , &
    & growth_phase                    , &
    & nd_dormance)

    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_jsb_impl_constants,      ONLY: true, false, test_false_true
    USE mo_q_agr_constants,         ONLY: mavg_period_cropseason, gdd_t_thres_dormseason, ctlib_gdd_max,   &
      &                                   ctlib_t_thres_gdd_gs, ctlib_t_thres_planting, ctlib_nd_dorm_max, &
      &                                   ctlib_gdd_mat_max,ix_planting, ix_harvest
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,      INTENT(in)    :: nc                                  !< dimensions
    REAL(wp),     INTENT(in)    :: dtime                               !< timestep length
    REAL(wp),     INTENT(in)    :: beta_soil_flush                     !< lctlib parameter
    REAL(wp),     INTENT(in)    :: gdd_req_max                         !< lctlib parameter
    REAL(wp),     INTENT(in)    :: beta_soil_senescence                !< lctlib parameter
    REAL(wp),     INTENT(in)    :: t_air_senescence                    !< lctlib parameter
    REAL(wp),     INTENT(in)    :: cn_leaf                             !< lctlib parameter
    REAL(wp),     INTENT(in)    :: np_leaf                             !< lctlib parameter
    REAL(wp),     INTENT(in)    :: crop_type_index(:)                  !< crop_type_index
    REAL(wp),     INTENT(in)    :: t_air(:)                            !< air temp
    REAL(wp),     INTENT(in)    :: t_air_week_mavg(:)                  !< moving average air temperature over a week
    REAL(wp),     INTENT(in)    :: t_air_month_mavg(:)                 !< moving average air temperature over a month
    REAL(wp),     INTENT(in)    :: beta_soil_gs_tphen_mavg(:)          !< soil water stress indicator
    REAL(wp),     INTENT(in)    :: t_air_tphen_mavg(:)                 !< moving average air temperature over phenology-relevant period
    REAL(wp),     INTENT(in)    :: t_soil_srf_tphen_mavg(:)            !< moving average soil surface temperature
    REAL(wp),     INTENT(in)    :: gpp_tlabile_mavg(:)                 !< moving average gpp over tlabile
    REAL(wp),     INTENT(in)    :: maint_respiration_tlabile_mavg(:)   !< moving average maintenance respiration over tlabile
    REAL(wp),     INTENT(inout) :: growing_season(:)                   !< growing season Yes/No
    REAL(wp),     INTENT(inout) :: growth_req_n_tlabile_mavg(:)        !< moles N required for a unit of C growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: growth_req_p_tlabile_mavg(:)        !< moles P required for a unit of N growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: gdd(:)                              !< growing degree days
    REAL(wp),     INTENT(inout) :: gdd_mavg(:)                         !< long-term average growing degree days
    REAL(wp),     INTENT(inout) :: growth_phase(:)                     !< category for growth phase of crops
    REAL(wp),     INTENT(inout) :: nd_dormance(:)                      !< number of days of dormancy
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ic                                  !< loop over grid cells
    INTEGER                     :: ix_ct                               !< index of crop type
    REAL(wp)                    :: hlp1                                !< helper variable for tmp values
    REAL(wp)                    :: growing_season_old                  !< growing season state of last time step
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_crop_phenology'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc
      !>  0.9 memorise old growing season state
      !>
      ix_ct = INT(crop_type_index(ic))
      growing_season_old = growing_season(ic)

      !>  1.0 check growing season has started or ended
      !>
      ! If not in growing season, check whether growing season has started
      IF (growing_season(ic) < test_false_true) THEN
        ! if temperature and soil moisture are not limiting growth
        IF (gdd(ic) > gdd_req_max .AND. beta_soil_gs_tphen_mavg(ic) > beta_soil_flush &
            & .AND. t_air_week_mavg(ic) > ctlib_t_thres_planting(ix_ct) .AND. t_soil_srf_tphen_mavg(ic) > TZero &
            & .AND. nd_dormance(ic) > ctlib_nd_dorm_max(ix_ct)) THEN
          growing_season(ic) = true
          growth_req_n_tlabile_mavg(ic) = 1._wp / cn_leaf
          growth_req_p_tlabile_mavg(ic) = 1._wp / np_leaf
        END IF
      ! If in the growing season, check whether growing season has ended
      ELSE
        ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND
        ! harvesting has occurred (ie no end of growing season during crop growth)
        IF ((t_air_tphen_mavg(ic) < t_air_senescence .OR. beta_soil_gs_tphen_mavg(ic) < beta_soil_senescence &
            & .OR. gpp_tlabile_mavg(ic) < maint_respiration_tlabile_mavg(ic) &
            & .OR. gdd(ic) > ctlib_gdd_mat_max(ix_ct)) .AND. &
            & .NOT.(growth_phase(ic) >= ix_planting .AND. growth_phase(ic) <= ix_harvest) ) THEN
          growing_season(ic) = false
        END IF
      END IF

      !>  2.0 update long-term average growing season length if now is the end of the season
      !>
      IF (growing_season(ic) <= test_false_true .AND. growing_season_old > test_false_true) THEN
        gdd_mavg(ic) = gdd_mavg(ic) * (1._wp - 1._wp / mavg_period_cropseason) &
          &                    + gdd(ic) * 1._wp / mavg_period_cropseason
      END IF

      !>  2.1 calculate current GDD increment
      !>
      !>    follows AgroIBIS in that there is a difference between the temperature
      !>    threshold for growing (crop-type specific) and dormant season
      !>    in the growing season, a maximum increment is applied
      IF (growing_season(ic) >= test_false_true) THEN
        hlp1 = MIN(MAX(t_air(ic) - ctlib_t_thres_gdd_gs(ix_ct),0.0_wp),ctlib_gdd_max(ix_ct))
      ELSE
        hlp1 = MAX(t_air(ic) - gdd_t_thres_dormseason,0.0_wp)
      END IF
      gdd(ic) = gdd(ic) + hlp1 * dtime / one_day

      !>  2.2 reset GDD counter at start and end of growing season
      !>       set growth phase to zero at start of growing season
      !>
      IF (growing_season(ic) <= test_false_true .AND. growing_season_old > test_false_true) THEN
        gdd(ic) = 0.0_wp
      END IF
      IF (growing_season(ic) > test_false_true .AND. growing_season_old <= test_false_true) THEN
        gdd(ic) = 0.0_wp
        growth_phase(ic) = 0.0_wp
      END IF
      ! reset gdd when there is snow coverage and soil is frozen
      IF (t_air_tphen_mavg(ic) < Tzero .OR. t_soil_srf_tphen_mavg(ic) <= Tzero) THEN
        gdd(ic) = 0.0_wp
      END IF

      !>  2.3 calculate current dormant days increment
      !>
      IF (growing_season(ic) < test_false_true) THEN
        nd_dormance(ic) = nd_dormance(ic) + dtime / one_day
      ELSE
        nd_dormance(ic) = 0.0_wp
      END IF

   END DO

  END SUBROUTINE calc_crop_phenology

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate flux from seed-bed to establishment if this is a planting time step
  !!  (i.e. crop growth phase is 0)
  !! This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_planting_flux( &
      & nc, &
      & growth_phase, &
      & veg_pool_seed_bed, &
      & veg_establishment_mt)

    USE mo_lnd_bgcm_idx
    USE mo_q_agr_constants,         ONLY: crop_planting_mass, ix_planting
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                        !< dimensions
    REAL(wp), INTENT(in)    :: growth_phase(:)           !< growth phase (0 = planting)
    REAL(wp), INTENT(in)    :: veg_pool_seed_bed(:,:)    !< bgcm vegetation seed-bed pool
    REAL(wp), INTENT(inout) :: veg_establishment_mt(:,:) !< bgcm flux of establishment
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                !< loop over grid cells
    REAL(wp)       :: planting_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_planting_flux'
    ! ----------------------------------------------------------------------------------------------------- !

    ! if growth phase equal zero and there is a seed pool, add transfer initial planting mass to establishment
    ! flux
    DO ic = 1,nc
      IF (INT(growth_phase(ic)) == ix_planting .AND. veg_pool_seed_bed(ixC, ic) > eps8 ) THEN
        planting_rate = MIN(crop_planting_mass, veg_pool_seed_bed(ixC, ic)) / veg_pool_seed_bed(ixC, ic)
        veg_establishment_mt(:,ic) = veg_establishment_mt(:,ic) + planting_rate * veg_pool_seed_bed( :, ic)
      END IF
    END DO

  END SUBROUTINE calc_planting_flux

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to set N fixation status to true (1) or false (0) depending on croptype
  !!  only soybeans are N fixers
  !! This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION set_crop_n_fixation_status(crop_type_index) RESULT (active_n_fixation)

    USE mo_q_agr_constants,         ONLY: ctlib_active_n_fixation

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: crop_type_index       !< crop type of current tile/gridcell
    REAL(wp)             :: active_n_fixation
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_crop_n_fixation_status'

    active_n_fixation = ctlib_active_n_fixation(INT(crop_type_index))

  END FUNCTION set_crop_n_fixation_status

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate current crop growth phase
  !!  This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation
  !!   Different to CLM5: crop growth phase 1 starts immediately after planting. This is because
  !!     inital growth has a larger fraction of roots, hence LAI development is slower than in CLM5
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_crop_growth_phase( &
      & nc, &
      & dtime, &
      & crop_type_index, &
      & gdd, &
      & gdd_mavg, &
      & lai, &
      & nd_crop_season, &
      & nd_crop_season_mavg, &
      & growth_phase)

    USE mo_jsb_math_constants,      ONLY: one_day
    USE mo_q_agr_constants
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                       !< dimensions
    REAL(wp), INTENT(in)    :: dtime                    !< timestep length
    REAL(wp), INTENT(in)    :: crop_type_index(:)       !< index of crop-type occupying this tile
    REAL(wp), INTENT(in)    :: gdd(:)                   !< current growing degree days (degC day)
    REAL(wp), INTENT(in)    :: gdd_mavg(:)              !< long-term average gdd (degC day)
    REAL(wp), INTENT(in)    :: lai(:)                   !< current leaf area index
    REAL(wp), INTENT(inout) :: nd_crop_season(:)        !< number of days in this growing season
    REAL(wp), INTENT(inout) :: nd_crop_season_mavg(:)   !< long-term average length of growing season
    REAL(wp), INTENT(inout) :: growth_phase(:)          !< growth phase (0 = planting,
                                                        !<                1 = emergence,
                                                        !<                2 = grain fill,
                                                        !<                3 = harvest
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                     !< loop over grid cells
    INTEGER        :: ix_ct                  !< crop type index
    REAL(wp)       :: hlp1
    REAL(wp)       :: gdd_mat                !< GDD at maturity
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_crop_growth_phase'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc

      !> 0.9  cropindex and GDD requirement for maturity
      !>
      ix_ct = INT(crop_type_index(ic))
      gdd_mat = MIN(MAX(gdd_mavg(ic),ctlib_gdd_mat_min(ix_ct)),ctlib_gdd_mat_max(ix_ct))

      !> 1.0 advance growth phase to next stage if conditions are met
      SELECT CASE (INT(growth_phase(ic)))
        !> 1.1  when the crop is in planting mode
        !>
        CASE (ix_planting)
          ! this only happens at the first time step of the growing season - move to emergence phase
          ! note this differs from AgroIBIS, which requires 5% of GDD_mat to be reached.
          ! this is a result of the labile pool in QUINCY, which ensures that crop growth starts only in
          ! favourable conditions.
          growth_phase(ic) = growth_phase(ic) + 1.0_wp

        !> 1.2  when the crop is in emergence mode
        !>
        CASE (ix_emergence)
          ! if GDD has accumulated to a given fraction of GDD maturity or
          ! LAI has reached crop specific LAI maximum, move to phase 2 (grain filling)
          IF (gdd(ic) > ctlib_fract_gdd_mat_gp2to3(ix_ct) * gdd_mat &
              & .OR. lai(ic) > ctlib_lai_max(ix_ct)) THEN
            growth_phase(ic) = growth_phase(ic) + 1.0_wp
          END IF

        !> 1.3  when the crop is in grainfilling mode
        !>
        CASE (ix_grainfilling)
          ! if GDD maturity is reached, move to phase 3 (harvest)
          IF (gdd(ic) > gdd_mat .OR. nd_crop_season(ic) >= ctlib_nd_gs_max(ix_ct)) THEN
            growth_phase(ic) = growth_phase(ic) + 1.0_wp
          END IF

        !> 1.4  when the crop is in harvesting mode
        !>
        CASE (ix_harvest)
          ! if LAI reaches zero during harvest, enter phase 4 (currently undefined - can be intercrop later)
          ! and update the long-term average crop growing season length
          IF (lai(ic) < eps8) THEN
            growth_phase(ic)        = growth_phase(ic) + 1.0_wp
            nd_crop_season_mavg(ic) = nd_crop_season_mavg(ic) * (1._wp - 1._wp / mavg_period_cropseason) &
              &                        + nd_crop_season(ic) * 1._wp / mavg_period_cropseason
          END IF

        CASE DEFAULT ! nothing to do, wait for reintialisation by crop phenology routine

      END SELECT

      !> 2.0  advance nd_crop_season counter, or set to zero outside crop growing season
      !>
      IF (growth_phase(ic) <= ix_harvest) THEN
        nd_crop_season(ic) = nd_crop_season(ic) + dtime / one_day
      ELSE
        nd_crop_season(ic) = 0.0_wp
      END IF

    END DO

  END SUBROUTINE update_crop_growth_phase

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate allometry of crops as well as fraction of carbon allocation to fruit
  !!  production. This follows the CLM5 implementation of AgroIBIS, see CLM5 documentation, but has been
  !!  adjusted to match QUINCY's allocation routine, which is based on allometry, not fixed fractions
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_crop_allocation_factors( &
      & nc, &
      & crop_type_index, &
      & growth_phase, &
      & gdd, &
      & gdd_mavg, &
      & lai, &
      & leaf2sapwood_mass_ratio, &
      & leaf2root_mass_ratio, &
      & falloc_fruit_crop)

    USE mo_q_agr_constants
    USE mo_veg_constants,      ONLY: leaf2root_min_ratio
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                          !< dimensions
    REAL(wp), INTENT(in)    :: crop_type_index(:)          !< index of crop-type occupying this tile
    REAL(wp), INTENT(in)    :: growth_phase(:)             !< growth phase (0 = planting, ..., 3 = harvest)
    REAL(wp), INTENT(in)    :: gdd(:)                      !< current growing degree days (degC day)
    REAL(wp), INTENT(in)    :: gdd_mavg(:)                 !< long-term average gdd (degC day)
    REAL(wp), INTENT(in)    :: lai(:)                      !< current leaf area index
    REAL(wp), INTENT(out)   :: leaf2sapwood_mass_ratio(:)  !< leaf to halm mass ratio (unitless)
    REAL(wp), INTENT(out)   :: leaf2root_mass_ratio(:)     !< leaf to fine root mass ratio (unitless)
    REAL(wp), INTENT(out)   :: falloc_fruit_crop(:)        !< fraction of allocation to fruits
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                     !< loop over grid cells
    INTEGER        :: ix_ct                  !< crop type index
    REAL(wp)       :: hlp1, hlp2
    REAL(wp)       :: gdd_mat                !< GDD at maturity
    REAL(wp)       :: heat_index             !< GDD based heat index
    REAL(wp)       :: fheat_index            !< GDD based heat index
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_crop_allocation_factors'

    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc

      !> 0.9  crop type index and other helper variables
      !>
      ix_ct = INT(crop_type_index(ic))
      ! GDD needed to reach maturity
      gdd_mat = MIN(MAX(gdd_mavg(ic),ctlib_gdd_mat_min(ix_ct)),ctlib_gdd_mat_max(ix_ct))
      ! Heat index (defined as the fraction of GDD for maturity that induces grain filling
      heat_index = ctlib_fract_gdd_mat_gp2to3(ix_ct) * gdd_mat
      ! Normalised heat index
      fheat_index = MIN(gdd(ic) / heat_index, 1._wp)

      !> 1.0  allometry of leaf, fine root and stem
      !>
      !> This is adapted from AgroIBIS but reformulated to transform from a fraction based to an
      !> allometry based allocation (consistent with the assumptions in mo_veg_growth:calc_allocation_fraction.
      !> This means that the parameters have been adjusted to match implied outcome
      !> Allocation is responding to cummulated heat sums during emergence (declining allocation to leaves
      !> with rising maturity). Once the heat requirement for grainfilling is met, allocation is more strongly
      !> towards roots.
      !>
      hlp1 = (exp(-k_heatindex) - exp(-k_heatindex*fheat_index))/(exp(-k_heatindex) - 1._wp)
      leaf2root_mass_ratio(ic) = ctlib_leaf2root_mass_fin(ix_ct) &
        &                        + (ctlib_leaf2root_mass_init(ix_ct) - ctlib_leaf2root_mass_fin(ix_ct)) * hlp1
      leaf2sapwood_mass_ratio(ic) = 1._wp / (ctlib_sapwood2leaf_mass_init(ix_ct) &
        &                                    + ctlib_sapwood2leaf_mass_fin(ix_ct) * (1._wp - hlp1))
      ! when grainfilling occurs, modify allocation fractions further
      IF (gdd(ic) > heat_index) THEN
        IF (gdd_mat - gdd(ic) > eps8) THEN
          hlp1 = (1._wp - ((gdd(ic) - heat_index) / (gdd_mat - heat_index)))**ctlib_k_hi_leaf(ix_ct)
          hlp2 = (1._wp - ((gdd(ic) - heat_index) / (gdd_mat - heat_index)))**ctlib_k_hi_stem(ix_ct)
        ELSE
          hlp1 = 0.0_wp
          hlp2 = 0.0_wp
        END IF

        leaf2root_mass_ratio(ic) = MAX(leaf2root_mass_ratio(ic) * hlp1, leaf2root_min_ratio)
        leaf2sapwood_mass_ratio(ic) = MAX(leaf2sapwood_mass_ratio(ic) * hlp2,leaf2sapwood_min_ratio)
      END IF
      ! if maximum LAI is reached, reduce leaf allocation to minimum
      IF (lai(ic) > ctlib_lai_max(ix_ct)) leaf2root_mass_ratio(ic) = leaf2root_min_ratio

      !> 2.0  calculate grainfill ratio for crop allocation
      !>
      !> this is adapted from AgroIBIS in that instead of decreasing the leaf and shoot allocation fraction
      !> by (1-hlp1)^d, we increase the fruit allocation fraction as accordingly. This is needed to make
      !> the calculation consistent with the assumptions in calc_allocation_fraction is mo_veg_growth
      !>
      IF (INT(growth_phase(ic)) == ix_grainfilling .AND. gdd(ic) > heat_index) THEN
        hlp1 = (gdd(ic) - heat_index) / ((1._wp - ctlib_fract_gdd_mat_gp2to3(ix_ct)) * gdd_mat)
        falloc_fruit_crop(ic) = ctlib_falloc_fruit_max(ix_ct) * (1._wp - (1._wp - hlp1)**ctlib_k_hi_leaf(ix_ct))
      ELSE
        falloc_fruit_crop(ic) = 0.0_wp
      END IF

    END DO

  END SUBROUTINE calc_crop_allocation_factors

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_cropland_dynamics
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to deal with fertiliser application
  !!  This follows OCN in that at specific days during the crop growing season
  !!  (expressed as a fraction of the number of days the crop growing season lasts) a given rate
  !!  of fertiliser is applied.
  !!  P fertiliser is assumed a time invariant constant fraction of N fertiliser. This is only a first
  !!  order guess
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_fertiliser_application( &
      & nc, &
      & lctlib_ps_pathway, &
      & nd_crop_season, &
      & nd_crop_season_mavg, &
      & n_fertiliser_c3, &
      & n_fertiliser_c4, &
      & fertiliser_nh4, &
      & fertiliser_nh4_n15, &
      & fertiliser_no3, &
      & fertiliser_no3_n15, &
      & fertiliser_po4)

    USE mo_jsb_math_constants,     ONLY: one_day
    USE mo_jsb_impl_constants,     ONLY: test_false_true
    USE mo_q_agr_constants,        ONLY: fcs_fertappl_day1, fcs_fertappl_day2, fcs_fertappl_day3, &
                                         fcs_fertappl_day4, frac_fertil_nh4, frac_fertil_po4, eta_fertiliser
    USE mo_q_assimi_constants,     ONLY: ic3phot, ic4phot
    USE mo_isotope_util,           ONLY: calc_mixing_ratio_N15N14
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,  INTENT(in)    :: nc                         !< dimensions
    INTEGER,  INTENT(in)    :: lctlib_ps_pathway          !< photosynthesis type (ic3phot|ic4phot)
    REAL(wp), INTENT(in)    :: nd_crop_season(:)       !< number of days this growing season has lasted
    REAL(wp), INTENT(in)    :: nd_crop_season_mavg(:)  !< number of days this growing season has lasted
    REAL(wp), INTENT(in)    :: n_fertiliser_c3(:)         !< annual N fertiliser rate for C3 crops (mol/m2/yr)
    REAL(wp), INTENT(in)    :: n_fertiliser_c4(:)         !< annual N fertiliser rate for C4 crops (mol/m2/yr)
    REAL(wp), INTENT(inout) :: fertiliser_nh4(:),       & !< NH4 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_nh4_n15(:),   & !< 15NH4 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_no3(:),       & !< NO3 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_no3_n15(:),   & !< 15NO3 fertiliser application [micro mol / m2 / s]
      &                        fertiliser_po4(:)          !< PO4 fertiliser application [micro mol / m2 / s]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER        :: ic                !< loop over grid cells
    REAL(wp)       :: hlp1
    LOGICAL        :: is_fertiliser_day
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fertiliser_application'
    ! ----------------------------------------------------------------------------------------------------- !

    DO ic = 1,nc
      !> 1.0 determine if today is a fertiliser day (given that this is a FLOOR, this is true for 24h)
      !>     the first application day is forced to be the first day after planting to avoid fertilisation
      !>     outside the growing season
      !>
      IF (FLOOR(nd_crop_season(ic)) == MAX(1,FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day1)) &
          .OR. FLOOR(nd_crop_season(ic)) == FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day2)   &
          .OR. FLOOR(nd_crop_season(ic)) == FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day3)   &
          .OR. FLOOR(nd_crop_season(ic)) == FLOOR(nd_crop_season_mavg(ic)*fcs_fertappl_day4)) THEN
        is_fertiliser_day = .TRUE.
      ELSE
        is_fertiliser_day = .FALSE.
      END IF
      !> 2.0 Mineral fertiliser application for croplands
      !>
      IF (is_fertiliser_day) THEN
        IF (lctlib_ps_pathway == ic3phot) THEN
          hlp1 = n_fertiliser_c3(ic) / one_day * 1.e6_wp / 4._wp ! conversion and split into 4 doses
        ELSE
          hlp1 = n_fertiliser_c4(ic) / one_day * 1.e6_wp / 4._wp ! conversion and split into 4 doses
        END IF
        fertiliser_nh4(ic) = hlp1 * frac_fertil_nh4
        fertiliser_nh4_n15(ic) = fertiliser_nh4(ic) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(eta_fertiliser))
        fertiliser_no3(ic) = hlp1 * (1._wp - frac_fertil_nh4)
        fertiliser_no3_n15(ic) = fertiliser_no3(ic) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(eta_fertiliser))
        fertiliser_po4(ic) = hlp1 * frac_fertil_po4
      END IF
    END DO

  END SUBROUTINE calc_fertiliser_application

#endif
END MODULE mo_q_agr_process
