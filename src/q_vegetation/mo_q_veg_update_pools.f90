!> QUINCY update vegetation pools
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
!>#### routines for calculating the vegetation pools using the fluxes that have been calculated
!>
MODULE mo_q_veg_update_pools
#ifndef __NO_QUINCY__

  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_TOTAL_BIO_ID, VEG_BGCM_GROWTH_ID, &
    &                                   VEG_BGCM_LITTERFALL_ID, VEG_BGCM_EXUDATION_ID, VEG_BGCM_ESTABLISHMENT_ID, &
    &                                   VEG_BGCM_RESERVE_USE_ID, VEG_BGCM_PP_FUEL_ID, VEG_BGCM_PP_PAPER_ID, &
    &                                   VEG_BGCM_PP_FIBERBOARD_ID, VEG_BGCM_PP_OIRW_ID, VEG_BGCM_PP_PV_ID,  &
    &                                   VEG_BGCM_PP_SAWNWOOD_ID, VEG_BGCM_FPROD_DECAY_ID
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_MYCO_EXPORT_ID

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_veg_pools

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_update_pools'

CONTAINS

  ! ======================================================================================================= !
  !>update vegetation pools
  !>
  SUBROUTINE update_veg_pools(tile, options)
    USE mo_kind,                            ONLY: wp
    USE mo_jsb_impl_constants,              ONLY: test_false_true
    USE mo_jsb_class,                       ONLY: Get_model
    USE mo_jsb_tile_class,                  ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,                  ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,                ONLY: t_lctlib_element
    USE mo_jsb_model_class,                 ONLY: t_jsb_model
    USE mo_quincy_model_config,             ONLY: QLAND, QPLANT, QSOIL, QCANOPY, Q_TEST_CANOPY, Q_TEST_RADIATION
    USE mo_jsb_process_class,               ONLY: SPQ_, VEG_, RAD_, Q_ASSIMI_, Q_PHENO_, SB_
    USE mo_jsb_grid_class,                  ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                        ONLY: Get_vgrid
    USE mo_jsb_physical_constants,          ONLY: lambda_C14, rhoh2o, grav
    USE mo_jsb_math_constants,              ONLY: one_day, one_year, eps4, eps8
    USE mo_isotope_util,                    ONLY: calc_mixing_ratio_N15N14, calc_fractionation, calc_mixing_ratio_C14C
    USE mo_q_veg_canopy,                    ONLY: calc_canopy_layers
    USE mo_q_veg_growth,                    ONLY: calc_diameter_from_woody_biomass, calc_height_from_diameter
    USE mo_veg_constants,                   ONLY: itree, eta_nfixation, max_leaf_shedding_rate, k_leafon_canopy, &
      &                                           fresp_growth, k_fpc, k_sai2lai, rfr_ratio_toc
    USE mo_q_assimi_parameters,             ONLY: CiCa_default_C3, CiCa_default_C4
    USE mo_q_assimi_constants,              ONLY: ic3phot, ic4phot
    USE mo_atmland_constants,               ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
      &                                           def_co2_deltaC13, def_co2_deltaC14
    USE mo_q_assimi_process,                ONLY: discrimination_ps
    USE mo_q_pheno_constants,               ONLY: ievergreen
    USE mo_exception,                       ONLY: finish
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(Q_ASSIMI_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(RAD_)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_jsb_vgrid),        POINTER :: vgrid_canopy_q_assimi  !< Vertical grid
    TYPE(t_jsb_vgrid),        POINTER :: vgrid_soil_sb          !< Vertical grid
    INTEGER                           :: ncanopy                !< number of canopy layers
    INTEGER                           :: nsoil_sb               !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(options%nc)   :: leaf_nitrogen          !< helper array for leaf nitrogen
    REAL(wp), ALLOCATABLE             :: vgrid_canopy_dz(:)     !< canopy grid helper variable dz
    REAL(wp), ALLOCATABLE             :: vgrid_canopy_lbounds(:)!< canopy grid helper variable lbounds
    REAL(wp), ALLOCATABLE             :: vgrid_canopy_ubounds(:)!< canopy grid helper variable ubounds
    REAL(wp)                          :: lctlib_g1              !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    REAL(wp)                          :: dtime                  !< timestep length
    INTEGER                           :: iblk, ics, ice, nc     !< grid dimensions
    INTEGER                           :: ic, is, icanopy        !< looping indices
    INTEGER                           :: id_elem, ix_elem       !< id and index of bgcm elements
    INTEGER                           :: ix_comp                !< currently id and index of bgcm compartments
    INTEGER                           :: nr_of_veg_bgcm_comp    !< number of compartments in veg bgcms
    CHARACTER(len=1024)               :: canopy_cond_scheme     !< canopy_conductance_scheme: medlyn / ballberry - q_assimi config
    LOGICAL                           :: flag_optimal_Nfraction !< on/off optimise leaf internal N allocation - q_assimi config
    REAL(wp) :: lctlib_sla, lctlib_tau_leaf, lctlib_cn_leaf, lctlib_np_leaf, lctlib_allom_k1, lctlib_allom_k2
    REAL(wp) :: lctlib_wood_density, lctlib_phi_leaf_min, lctlib_k_latosa, lctlib_k0_fn_struc, lctlib_fn_oth_min
    REAL(wp) :: lctlib_gmin, lctlib_g0, lctlib_t_jmax_omega
    REAL(wp) :: hlp
    INTEGER  :: lctlib_phenology_type, lctlib_growthform, lctlib_ps_pathway
    LOGICAL  :: config_l_use_product_pools
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_pools'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L2D :: veg_growth_mt
    dsl4jsb_Def_mt1L2D :: veg_exudation_mt
    dsl4jsb_Def_mt1L2D :: veg_establishment_mt
    dsl4jsb_Def_mt1L2D :: veg_reserve_use_mt
    dsl4jsb_Def_mt1L2D :: veg_total_biomass_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt

    dsl4jsb_Def_mt1L2D :: veg_pp_fuel_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_paper_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_fiberboard_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_oirw_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_pv_mt
    dsl4jsb_Def_mt1L2D :: veg_pp_sawnwood_mt
    dsl4jsb_Def_mt1L2D :: fprod_decay_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(Q_ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: lai_max
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: w_soil_root_pot
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    ! RAD_ 3D
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: do_cohort_harvest
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: sai
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: fract_fpc
    dsl4jsb_Real2D_onChunk      :: blended_height
    dsl4jsb_Real2D_onChunk      :: dphi
    dsl4jsb_Real2D_onChunk      :: npp
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration
    dsl4jsb_Real2D_onChunk      :: net_growth
    dsl4jsb_Real2D_onChunk      :: npp_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c13
    dsl4jsb_Real2D_onChunk      :: npp_c14
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c14
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: uptake_nh4_n15
    dsl4jsb_Real2D_onChunk      :: uptake_no3_n15
    dsl4jsb_Real2D_onChunk      :: n_fixation_n15
    dsl4jsb_Real2D_onChunk      :: uptake_po4
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: t_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: net_biosphere_production
    dsl4jsb_Real2D_onChunk      :: biological_n_fixation
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_p
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_pool_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_wood_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_wood_n
    dsl4jsb_Real2D_onChunk      :: veg_pool_fine_root_c
    dsl4jsb_Real2D_onChunk      :: veg_pool_fine_root_n
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_n
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_p
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_growth_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_n
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_p
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c13
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_c14
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_total_n15
    dsl4jsb_Real2D_onChunk      :: veg_products_total_c
    dsl4jsb_Real2D_onChunk      :: veg_products_total_n
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: root_fraction_sl
    dsl4jsb_Real3D_onChunk      :: delta_root_fraction_sl
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    dtime     = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(tile%owner_model_id)
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    vgrid_soil_sb         => Get_vgrid('soil_layer_sb')
    ncanopy               =  vgrid_canopy_q_assimi%n_levels
    nsoil_sb              =  vgrid_soil_sb%n_levels

    ! ----------------------------------------------------------------------------------------------------- !
    IF (dsl4jsb_Lctlib_param(BareSoilFlag)) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(Q_ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    config_l_use_product_pools = dsl4jsb_Config(VEG_)%l_use_product_pools

    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_GROWTH_ID, veg_growth_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_EXUDATION_ID, veg_exudation_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_ESTABLISHMENT_ID, veg_establishment_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_RESERVE_USE_ID, veg_reserve_use_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_TOTAL_BIO_ID, veg_total_biomass_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)

    IF (config_l_use_product_pools) THEN
      dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FUEL_ID, veg_pp_fuel_mt)
      dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PAPER_ID, veg_pp_paper_mt)
      dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_FIBERBOARD_ID, veg_pp_fiberboard_mt)
      dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_OIRW_ID, veg_pp_oirw_mt)
      dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_PV_ID, veg_pp_pv_mt)
      dsl4jsb_Get_mt1L2D(VEG_BGCM_PP_SAWNWOOD_ID, veg_pp_sawnwood_mt)
      dsl4jsb_Get_mt1L2D(VEG_BGCM_FPROD_DECAY_ID, fprod_decay_mt)
    ENDIF
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation)        ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C13)    ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C14)    ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air_tfrac_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soa_tphen_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps_tfrac_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tfrac_mavg)   ! in
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)             ! in
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, lai_max)                    ! in
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_pot)                ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)                  ! in
    ! RAD_ 3D
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)      ! in
    ! VEG 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, do_cohort_harvest)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                               ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, lai)                                  ! CASE "plant" & "land": out / other Cases: in
    dsl4jsb_Get_var2D_onChunk(VEG_, sai)                                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, fract_fpc)                            ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, blended_height)                       ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi)                                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)                                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth)                           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c13)                              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c13)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c14)                              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c14)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_n15)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_po4)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tfrac_mavg)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tfrac_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tfrac_mavg)                        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg)           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)                      ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_biosphere_production)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, biological_n_fixation)                ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_p)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c13)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_c14)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_total_n15)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_c)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_n)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_wood_c)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_wood_n)                      ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_fine_root_c)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_fine_root_n)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_n)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_p)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c13)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_c14)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_total_n15)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_n)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_p)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c13)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_c14)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_total_n15)             ! out
    IF (config_l_use_product_pools) THEN
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_c)               ! out
      dsl4jsb_Get_var2D_onChunk(VEG_, veg_products_total_n)               ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(VEG_, rfr_ratio_boc_tvegdyn_mavg)           ! in
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)     ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)               ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)                       ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)                     ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)                      ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, root_fraction_sl)               ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, delta_root_fraction_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                         ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)                    ! out
    ! ----------------------------------------------------------------------------------------------------- !

    ALLOCATE(vgrid_canopy_dz(ncanopy), vgrid_canopy_lbounds(ncanopy), vgrid_canopy_ubounds(ncanopy))
    !$ACC ENTER DATA CREATE(vgrid_canopy_dz(:), vgrid_canopy_lbounds(:), vgrid_canopy_ubounds(:))

    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO icanopy = 1, ncanopy
      vgrid_canopy_dz(icanopy) = vgrid_canopy_q_assimi%dz(icanopy)
      vgrid_canopy_lbounds(icanopy) = vgrid_canopy_q_assimi%lbounds(icanopy)
      vgrid_canopy_ubounds(icanopy) = vgrid_canopy_q_assimi%ubounds(icanopy)
    END DO
    !$ACC END PARALLEL LOOP

    !> 0.8 get number of bgcm compartments
    !>
    nr_of_veg_bgcm_comp = SIZE(veg_pool_mt,1)

    lctlib_sla = dsl4jsb_Lctlib_param(sla)
    lctlib_phenology_type = dsl4jsb_Lctlib_param(phenology_type)
    lctlib_tau_leaf = dsl4jsb_Lctlib_param(tau_leaf)
    lctlib_cn_leaf = dsl4jsb_Lctlib_param(cn_leaf)
    lctlib_np_leaf = dsl4jsb_Lctlib_param(np_leaf)
    lctlib_growthform = dsl4jsb_Lctlib_param(growthform)
    lctlib_allom_k1 = dsl4jsb_Lctlib_param(allom_k1)
    lctlib_allom_k2 = dsl4jsb_Lctlib_param(allom_k2)
    lctlib_wood_density = dsl4jsb_Lctlib_param(wood_density)
    lctlib_phi_leaf_min = dsl4jsb_Lctlib_param(phi_leaf_min)
    lctlib_k_latosa = dsl4jsb_Lctlib_param(k_latosa)
    lctlib_ps_pathway = dsl4jsb_Lctlib_param(ps_pathway)

    lctlib_k0_fn_struc = dsl4jsb_Lctlib_param(k0_fn_struc)
    lctlib_fn_oth_min = dsl4jsb_Lctlib_param(fn_oth_min)
    lctlib_gmin = dsl4jsb_Lctlib_param(gmin)
    lctlib_g0 = dsl4jsb_Lctlib_param(g0)
    lctlib_t_jmax_omega = dsl4jsb_Lctlib_param(t_jmax_omega)

    canopy_cond_scheme = TRIM(dsl4jsb_Config(Q_ASSIMI_)%canopy_conductance_scheme)
    flag_optimal_Nfraction = dsl4jsb_Config(Q_ASSIMI_)%flag_optimal_Nfraction

    !>0.9 set g1 according to canopy_conductance_scheme
    !>
    SELECT CASE(canopy_cond_scheme)
      CASE ("medlyn")
        lctlib_g1 = dsl4jsb_Lctlib_param(g1_medlyn)
      CASE ("ballberry")
        lctlib_g1 = dsl4jsb_Lctlib_param(g1_bberry)
    END SELECT

    !>1.0 for plant and land model, update vegetation pools with fluxes
    !>
    SELECT CASE(model%config%qmodel_id)
    CASE(QPLANT, QLAND)
#ifdef _OPENACC
      CALL finish(routine, 'Code blocks for QPLANT & QLAND modes not ported to GPU, yet. Stop.')
#endif

      DO ic = 1, nc

       !>  1.1 net balance of the labile carbon pool due to photosynthesis and respiration
       !>
       npp(ic)        = gross_assimilation(ic)     - &
                        growth_respiration(ic)     - maint_respiration(ic)     - n_processing_respiration(ic)
       npp_c13(ic)    = gross_assimilation_C13(ic) - &
                        growth_respiration_c13(ic) - maint_respiration_c13(ic) - n_processing_respiration_c13(ic)
       npp_c14(ic)    = gross_assimilation_C14(ic) - &
                        growth_respiration_c14(ic) - maint_respiration_c14(ic) - n_processing_respiration_c14(ic)
       net_growth(ic) = net_growth(ic) + npp(ic)

       veg_growth_mt(ix_labile, ixC, ic)   = veg_growth_mt(ix_labile, ixC, ic)   + npp(ic)     * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixC13, ic) = veg_growth_mt(ix_labile, ixC13, ic) + npp_c13(ic) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixC14, ic) = veg_growth_mt(ix_labile, ixC14, ic) + npp_c14(ic) * dtime / 1000000.0_wp

       !>  1.2 Add nutrients from uptake calculations
       !>
       n_fixation_n15(ic)                  = n_fixation(ic) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(-eta_nfixation))
       veg_growth_mt(ix_labile, ixN, ic)   = veg_growth_mt(ix_labile, ixN, ic) + &
                                             ( uptake_nh4(ic) + uptake_no3(ic) + n_fixation(ic)) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixN15, ic) = veg_growth_mt(ix_labile, ixN15, ic) + &
                                             ( uptake_nh4_n15(ic) + uptake_no3_n15(ic) + n_fixation_n15(ic)) * dtime / 1000000.0_wp
       veg_growth_mt(ix_labile, ixP, ic)   = veg_growth_mt(ix_labile, ixP, ic) + &
                                             uptake_po4(ic) * dtime / 1000000.0_wp

       !>  1.3 Add mycorrhiza export
       !>
       DO is = 1, nsoil_sb
        veg_growth_mt(ix_labile, ixC, ic)   = veg_growth_mt(ix_labile, ixC, ic)   + &
          &                                   sb_mycorrhiza_export_mt(ixC, ic, is)    * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixC13, ic) = veg_growth_mt(ix_labile, ixC13, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixC13, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixC14, ic) = veg_growth_mt(ix_labile, ixC14, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixC14, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixN, ic)   = veg_growth_mt(ix_labile, ixN, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixN, ic, is)    * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixN15, ic) = veg_growth_mt(ix_labile, ixN15, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixN15, ic, is)  * soil_depth_sl(ic,is)
        veg_growth_mt(ix_labile, ixP, ic)   = veg_growth_mt(ix_labile, ixP, ic) + &
          &                                   sb_mycorrhiza_export_mt(ixP, ic, is)    * soil_depth_sl(ic,is)
       END DO

      END DO

       !>  1.4 update labile pool and reserve
       !>
       veg_pool_mt(ix_labile, :, :)  = veg_pool_mt(ix_labile, :, :)  + veg_reserve_use_mt(:,:)
       veg_pool_mt(ix_reserve, :, :) = veg_pool_mt(ix_reserve, :, :) - veg_reserve_use_mt(:,:)

      DO ic = 1, nc
        !>  1.5 update fluxes from continuous nutrient turnover and resorption
       !>
       ! N
       veg_pool_mt(ix_labile, ixN, ic)     = veg_pool_mt(ix_labile, ixN, ic) &
         &                                   + recycling_leaf_n(ic) + recycling_fine_root_n(ic) + recycling_heart_wood_n(ic)
       veg_pool_mt(ix_leaf, ixN, ic)       = veg_pool_mt(ix_leaf, ixN, ic)       - recycling_leaf_n(ic)
       veg_pool_mt(ix_fine_root, ixN, ic)  = veg_pool_mt(ix_fine_root, ixN, ic)  - recycling_fine_root_n(ic)
       veg_pool_mt(ix_heart_wood, ixN, ic) = veg_pool_mt(ix_heart_wood, ixN, ic) - recycling_heart_wood_n(ic)
       ! P
       veg_pool_mt(ix_labile, ixP, ic)     = veg_pool_mt(ix_labile, ixP, ic) &
         &                                   + recycling_leaf_p(ic) + recycling_fine_root_p(ic) + recycling_heart_wood_p(ic)
       veg_pool_mt(ix_leaf, ixP, ic)       = veg_pool_mt(ix_leaf, ixP, ic)       - recycling_leaf_p(ic)
       veg_pool_mt(ix_fine_root, ixP, ic)  = veg_pool_mt(ix_fine_root, ixP, ic)  - recycling_fine_root_p(ic)
       veg_pool_mt(ix_heart_wood, ixP, ic) = veg_pool_mt(ix_heart_wood, ixP, ic) - recycling_heart_wood_p(ic)
       ! N15
       veg_pool_mt(ix_labile, ixN15, ic)     = veg_pool_mt(ix_labile, ixN15, ic) &
         &                                    + recycling_leaf_n15(ic) + recycling_fine_root_n15(ic) + recycling_heart_wood_n15(ic)
       veg_pool_mt(ix_leaf, ixN15, ic)       = veg_pool_mt(ix_leaf, ixN15, ic)       - recycling_leaf_n15(ic)
       veg_pool_mt(ix_fine_root, ixN15, ic)  = veg_pool_mt(ix_fine_root, ixN15, ic)  - recycling_fine_root_n15(ic)
       veg_pool_mt(ix_heart_wood, ixN15, ic) = veg_pool_mt(ix_heart_wood, ixN15, ic) - recycling_heart_wood_n15(ic)

       !>  1.6 update labile pool with fluxes to mycorrhiza
       !>
       veg_pool_mt(ix_labile, ixC, ic)   = veg_pool_mt(ix_labile, ixC, ic)   - veg_exudation_mt(ixC, ic)
       veg_pool_mt(ix_labile, ixN, ic)   = veg_pool_mt(ix_labile, ixN, ic)   - veg_exudation_mt(ixN, ic)
       veg_pool_mt(ix_labile, ixP, ic)   = veg_pool_mt(ix_labile, ixP, ic)   - veg_exudation_mt(ixP, ic)
       veg_pool_mt(ix_labile, ixC13, ic) = veg_pool_mt(ix_labile, ixC13, ic) - veg_exudation_mt(ixC13, ic)
       veg_pool_mt(ix_labile, ixC14, ic) = veg_pool_mt(ix_labile, ixC14, ic) - veg_exudation_mt(ixC14, ic)
       veg_pool_mt(ix_labile, ixN15, ic) = veg_pool_mt(ix_labile, ixN15, ic) - veg_exudation_mt(ixN15, ic)
      END DO


       !>  1.7 Update plant biogeochemical pools with growth, litterfall due to natural processes and establishment
       !>
      DO ic = 1, nc
        DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
          IF (model%config%is_element_used(id_elem)) THEN
            ix_elem = model%config%elements_index_map(id_elem)    ! get element index in bgcm
            DO ix_comp = 1, nr_of_veg_bgcm_comp
              veg_pool_mt(ix_comp, ix_elem, ic) = veg_pool_mt(ix_comp, ix_elem, ic) - veg_litterfall_mt(ix_comp, ix_elem, ic)
              veg_pool_mt(ix_comp, ix_elem, ic) = veg_pool_mt(ix_comp, ix_elem, ic) + veg_growth_mt(ix_comp, ix_elem, ic)
            END DO
            ! set woody biomass to zero after stand replacing harvest
            ! NOTE: the IF may never be TRUE in IQ, and it will be replaced with future work on "simple harvest"
            IF (model%config%l_do_stand_replacing_harvest .OR. do_cohort_harvest(ic) > test_false_true) THEN
              veg_pool_mt(ix_sap_wood, ix_elem, ic)   = 0.0_wp
              veg_pool_mt(ix_heart_wood, ix_elem, ic) = 0.0_wp
            END IF
          END IF
        END DO
      END DO

      DO ic = 1, nc
        DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
          IF (model%config%is_element_used(id_elem)) THEN
            ix_elem = model%config%elements_index_map(id_elem)    ! get element index in bgcm
            veg_pool_mt(ix_labile,   ix_elem, ic) = veg_pool_mt(ix_labile,   ix_elem, ic) + veg_establishment_mt(ix_elem, ic)
            veg_pool_mt(ix_seed_bed, ix_elem, ic) = veg_pool_mt(ix_seed_bed, ix_elem, ic) - veg_establishment_mt(ix_elem, ic)
          END IF
        END DO
      END DO

      ! @TODO
      !>    1.7.1 quick-fix ensure veg_pool does not become negative
      !>          apply only with IQ as it changes results with QS
      !>
      !   ==> need additional debug code
      !       to output litterfall & growth fluxes + lat/lon/PFT for grid cells where veg_pool gets negative
#ifdef __QUINCY_STANDALONE__
#else
      DO ic = 1, nc
        DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
          IF (model%config%is_element_used(id_elem)) THEN
            ix_elem = model%config%elements_index_map(id_elem)    ! get element index in bgcm
            DO ix_comp = 1, nr_of_veg_bgcm_comp
              veg_pool_mt(ix_comp, ix_elem, ic) = MAX(0.0_wp, veg_pool_mt(ix_comp, ix_elem, ic))
            END DO
          END IF
        END DO
      END DO
#endif

      !>  1.8 Radioactive decay of C14
      !>
      DO ic = 1, nc
        ! in living plant tissue
        DO ix_comp = 1, nr_of_veg_bgcm_comp
          veg_pool_mt(ix_comp, ixC14, ic) = veg_pool_mt(ix_comp, ixC14, ic) * (1._wp - lambda_C14 * dtime)
        END DO

        ! and if simulating with product pools also within these
        IF (config_l_use_product_pools) THEN
          veg_pp_fuel_mt(ixC14, ic) = veg_pp_fuel_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
          veg_pp_paper_mt(ixC14, ic) = veg_pp_paper_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
          veg_pp_fiberboard_mt(ixC14, ic) = veg_pp_fiberboard_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
          veg_pp_oirw_mt(ixC14, ic) = veg_pp_oirw_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
          veg_pp_pv_mt(ixC14, ic) = veg_pp_pv_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
          veg_pp_sawnwood_mt(ixC14, ic) = veg_pp_sawnwood_mt(ixC14, ic) * (1._wp - lambda_C14 * dtime)
        END IF
      END DO

    CASE (QCANOPY)

      !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
      DO ic = 1, nc

      ! diagnose required change in leaf C given phenology to attain prescribed LAI max
       IF (growing_season(ic) > test_false_true) THEN
         veg_growth_mt(ix_leaf, ixC, ic) = lai_max(ic) / lctlib_sla &
           &                              * (1._wp - exp(-k_leafon_canopy &
           &                              * (lai_max(ic) - lai(ic)) / lai_max(ic))) &
           &                              * dtime / one_day
       ELSE
         veg_litterfall_mt(ix_leaf, ixC, ic) = MIN(veg_pool_mt(ix_leaf, ixC, ic), &
           &                                  lai(ic)/lctlib_sla * max_leaf_shedding_rate * dtime / one_day &
           &                                  * lai_max(ic) &
           &                                  / MAX(eps4, lai(ic)))
       END IF

       ! calc leaf turnover for evergreen trees (needed, for e.g. mean_leaf_age, as evergreens are always within growing_season)
       IF (lctlib_phenology_type == ievergreen) THEN
         veg_litterfall_mt(ix_leaf, ixC, ic) = veg_pool_mt(ix_leaf, ixC, ic) / lctlib_tau_leaf / one_day / one_year * dtime
       ENDIF

       ! update leaf C:N:P pools
       veg_pool_mt(ix_leaf, ixC, ic) = veg_pool_mt(ix_leaf, ixC, ic) + veg_growth_mt(ix_leaf, ixC, ic) &
         &                            - veg_litterfall_mt(ix_leaf, ixC, ic)
       veg_pool_mt(ix_leaf, ixN, ic) = veg_pool_mt(ix_leaf, ixC, ic) / lctlib_cn_leaf
       veg_pool_mt(ix_leaf, ixP, ic) = veg_pool_mt(ix_leaf, ixN, ic) / lctlib_np_leaf

      END DO
      !$ACC END PARALLEL LOOP

    END SELECT

    !>2.0 Update plant diagnostics given the above fluxes
    !>
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      lai(ic)      = veg_pool_mt(ix_leaf, ixC, ic) * lctlib_sla
      dens_ind(ic) = dens_ind(ic) + delta_dens_ind(ic)

      IF (dens_ind(ic) > eps8 .AND. lctlib_growthform == itree) THEN
        diameter(ic)  = calc_diameter_from_woody_biomass(lctlib_allom_k1                     , &
          &                                              lctlib_allom_k2                     , &
          &                                              lctlib_wood_density                 , &
          &                                              veg_pool_mt(ix_sap_wood, ixC, ic)   , &
          &                                              veg_pool_mt(ix_heart_wood, ixC, ic) , &
          &                                              dens_ind(ic))
        height(ic)    = calc_height_from_diameter(lctlib_allom_k1, lctlib_allom_k2, diameter(ic))
      END IF

      IF (veg_pool_mt(ix_leaf, ixC, ic) > eps8) THEN
        mean_leaf_age(ic) = (mean_leaf_age(ic) + dtime / one_day) &
          &                * (1._wp - veg_growth_mt(ix_leaf, ixC, ic) / veg_pool_mt(ix_leaf, ixC, ic)) &
          &                - mean_leaf_age(ic) * veg_litterfall_mt(ix_leaf, ixC, ic) / veg_pool_mt(ix_leaf, ixC, ic)
      ELSE
        mean_leaf_age(ic) = 0.0_wp
      END IF

      dphi(ic) = w_soil_root_pot(ic) - lctlib_phi_leaf_min - grav * rhoh2o * height(ic) * 1.e-6_wp
    END DO
    !$ACC END PARALLEL LOOP

    IF (lctlib_growthform == itree) THEN
      SELECT CASE(model%config%qmodel_id)
      CASE(QPLANT, QLAND)
        DO ic = 1, nc
          IF (height(ic) > eps8) THEN
            sai(ic) = k_sai2lai * veg_pool_mt(ix_sap_wood, ixC, ic) &
              &       * lctlib_k_latosa / lctlib_wood_density / height(ic)
          ELSE
            sai(ic) = 0.0_wp
          END IF
        END DO
      CASE(QCANOPY)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
        DO ic = 1, nc
          sai(ic) = k_sai2lai * lai_max(ic)
        END DO
        !$ACC END PARALLEL LOOP
      END SELECT
    END IF

    !>  2.1 calculate foliage projected cover fraction and blended_height for turbulence calculations
    !>
    ! note: the '* 1.0_wp' for calc of blended_height(:) is the minimum/default height in meter
    ! using the VEG_ min_height parameter (= 0.1_wp) instead would change simulation results
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      fract_fpc(ic) = 1.0_wp - EXP(-k_fpc * (lai(ic) + sai(ic)))
      blended_height(ic) = height(ic) * fract_fpc(ic) + (1._wp - fract_fpc(ic)) * 1.0_wp
    END DO
    !$ACC END PARALLEL LOOP

    !>  2.2 implied change of root fraction given root growth
    !>
    IF (dsl4jsb_Config(VEG_)%flag_dynamic_roots) THEN
#ifdef _OPENACC
      CALL finish(routine, 'Code block ported to GPU but not tested, yet. Stop.')
#endif

      !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) DEFAULT(PRESENT)
      DO is = 1, nsoil_sb
        DO ic = 1, nc
          root_fraction_sl(ic,is) = root_fraction_sl(ic,is) + delta_root_fraction_sl(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !> 3.0 Given current vegetation pools, update canopy layers
    !>
    !$ACC DATA CREATE(leaf_nitrogen)
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      leaf_nitrogen(ic) = veg_pool_mt(ix_leaf, ixN, ic)
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    CALL calc_canopy_layers( nc                                                       , & ! in
                             ncanopy                                                  , &
                             dtime                                                    , &
                             vgrid_canopy_dz(:)                                       , &
                             vgrid_canopy_lbounds(:)                                  , &
                             vgrid_canopy_ubounds(:)                                  , &
                             lctlib_ps_pathway                                        , &
                             lctlib_k0_fn_struc                                       , &
                             lctlib_fn_oth_min                                        , &
                             lctlib_sla                                               , &
                             lctlib_np_leaf                                           , &
                             lctlib_gmin                                              , &
                             lctlib_g0                                                , &
                             lctlib_g1                                                , &
                             lctlib_t_jmax_omega                                      , &
                             flag_optimal_Nfraction                                   , & ! Q_ASSIMI_ config
                             canopy_cond_scheme                                       , & ! Q_ASSIMI_ config (medlyn/ballberry)
                             leaf_nitrogen(:)                                         , &
                             lai(:)                                                   , &
                             ppfd_sunlit_tfrac_mavg_cl(:,:)                           , &
                             ppfd_shaded_tfrac_mavg_cl(:,:)                           , &
                             fleaf_sunlit_tfrac_mavg_cl(:,:)                          , &
                             fn_rub_cl(:,:)                                           , &
                             fn_et_cl(:,:)                                            , &
                             fn_pepc_cl(:,:)                                          , &
                             fn_chl_cl(:,:)                                           , &
                             fn_oth_cl(:,:)                                           , & ! in
                             lai_cl(:,:)                                              , & ! out
                             cumm_lai_cl(:,:)                                         , & ! out
                             leaf_nitrogen_cl(:,:)                                    , & ! out
                             t_air            = t_air_tfrac_mavg(:), &                    ! optional in
                             t_acclim         = t_air_tacclim_mavg(:), &
                             press_srf        = press_srf_tfrac_mavg(:), &
                             co2_mixing_ratio = co2_mixing_ratio_tfrac_mavg(:), &
                             aerodyn_cond     = ga_tfrac_mavg(:), &
                             beta_air         = beta_air_tfrac_mavg(:), &
                             beta_soa         = beta_soa_tphen_mavg(:), &
                             beta_soil_ps     = beta_soil_ps_tfrac_mavg(:), &
                             beta_sinklim_ps  = beta_sinklim_ps_tfrac_mavg(:), &
                             beta_soil_gs     = beta_soil_gs_tfrac_mavg(:), &
                             t_jmax_opt       = t_jmax_opt_mavg(:) )                      ! optional in
    !$ACC WAIT(1)
    !$ACC END DATA
    !$ACC EXIT DATA DELETE(vgrid_canopy_dz, vgrid_canopy_lbounds, vgrid_canopy_ubounds)
    DEALLOCATE(vgrid_canopy_dz, vgrid_canopy_lbounds, vgrid_canopy_ubounds)

    !> 4.0 analysis of mass conservation: sum compartments of veg_pool into veg_total_biomass
    !>
    DO id_elem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (model%config%is_element_used(id_elem)) THEN
        ix_elem = model%config%elements_index_map(id_elem)    ! get element index in bgcm
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic = 1, nc
          veg_total_biomass_mt(ix_elem,ic) = 0.0_wp
        END DO
        !$ACC END PARALLEL

        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
        DO ix_comp = 1, nr_of_veg_bgcm_comp
          DO ic = 1, nc
            veg_total_biomass_mt(ix_elem,ic) = veg_total_biomass_mt(ix_elem,ic) + veg_pool_mt(ix_comp, ix_elem, ic)
          END DO
        END DO
        !$ACC END PARALLEL
      END IF
    END DO

    !> 5.0 bgc_material diagnostics
    !>     sum up veg_pool, veg_growth and veg_litterfall across components for each particular element
    !>     or simply add value (of one/two/few pool variables) to a memory variable for aggregation and output
    !>

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ! init pool totals
      veg_pool_total_c(ic)         = 0.0_wp
      veg_pool_total_n(ic)         = 0.0_wp
      veg_pool_total_p(ic)         = 0.0_wp
      veg_pool_total_c13(ic)       = 0.0_wp
      veg_pool_total_c14(ic)       = 0.0_wp
      veg_pool_total_n15(ic)       = 0.0_wp

      IF (config_l_use_product_pools) THEN
        veg_products_total_c(ic) = 0.0_wp
        veg_products_total_n(ic) = 0.0_wp
      END IF

      ! init flux totals
      veg_growth_total_c(ic) = 0.0_wp
      veg_growth_total_n(ic) = 0.0_wp
      veg_growth_total_p(ic) = 0.0_wp
      veg_growth_total_c13(ic) = 0.0_wp
      veg_growth_total_c14(ic) = 0.0_wp
      veg_growth_total_n15(ic) = 0.0_wp

      veg_litterfall_total_c(ic) = 0.0_wp
      veg_litterfall_total_n(ic) = 0.0_wp
      veg_litterfall_total_p(ic) = 0.0_wp
      veg_litterfall_total_c13(ic) = 0.0_wp
      veg_litterfall_total_c14(ic) = 0.0_wp
      veg_litterfall_total_n15(ic) = 0.0_wp

      ! store component totals for aggregation
      veg_pool_leaf_c(ic)          = veg_pool_mt(ix_leaf, ixC, ic)
      veg_pool_leaf_n(ic)          = veg_pool_mt(ix_leaf, ixN, ic)
      veg_pool_wood_c(ic)          = veg_pool_mt(ix_sap_wood, ixC, ic) + veg_pool_mt(ix_heart_wood, ixC, ic)
      veg_pool_wood_n(ic)          = veg_pool_mt(ix_sap_wood, ixN, ic) + veg_pool_mt(ix_heart_wood, ixN, ic)
      veg_pool_fine_root_c(ic)     = veg_pool_mt(ix_fine_root, ixC, ic)
      veg_pool_fine_root_n(ic)     = veg_pool_mt(ix_fine_root, ixN, ic)
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP GANG VECTOR PRIVATE(hlp)
    DO ic = 1, nc
      hlp = 0._wp
      !$ACC LOOP REDUCTION(+: hlp)
      DO ix_comp = 1, nr_of_veg_bgcm_comp
        hlp = hlp + veg_pool_mt(ix_comp,ixC, ic)
      END DO
      veg_pool_total_c(ic) = hlp
    END DO
    !$ACC END PARALLEL

    !TODO: Consider resolving "Complex loop carried dependence of" as above?
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO ix_comp = 1, nr_of_veg_bgcm_comp
      !$ACC LOOP GANG VECTOR
      DO ic = 1, nc
        ! calculate pool totals
!        veg_pool_total_c(ic)         = veg_pool_total_c(ic) + veg_pool_mt(ix_comp,ixC, ic)
        veg_pool_total_n(ic)         = veg_pool_total_n(ic) + veg_pool_mt(ix_comp,ixN, ic)
        veg_pool_total_p(ic)         = veg_pool_total_p(ic) + veg_pool_mt(ix_comp,ixP, ic)
        veg_pool_total_c13(ic)       = veg_pool_total_c13(ic) + veg_pool_mt(ix_comp,ixC13, ic)
        veg_pool_total_c14(ic)       = veg_pool_total_c14(ic) + veg_pool_mt(ix_comp,ixC14, ic)
        veg_pool_total_n15(ic)       = veg_pool_total_n15(ic) + veg_pool_mt(ix_comp,ixN15, ic)

        IF (config_l_use_product_pools) THEN
          veg_products_total_c(ic) = veg_pp_fuel_mt(ixC,ic) + veg_pp_paper_mt(ixC,ic)        &
            &                       + veg_pp_fiberboard_mt(ixC,ic) + veg_pp_oirw_mt(ixC,ic)  &
            &                       + veg_pp_pv_mt(ixC,ic) + veg_pp_sawnwood_mt(ixC,ic)
          veg_products_total_n(ic) = veg_pp_fuel_mt(ixN,ic) + veg_pp_paper_mt(ixN,ic)        &
            &                       + veg_pp_fiberboard_mt(ixN,ic) + veg_pp_oirw_mt(ixN,ic)  &
            &                       + veg_pp_pv_mt(ixN,ic) + veg_pp_sawnwood_mt(ixN,ic)
        END IF

        ! calculate flux totals
        veg_growth_total_c(ic) = veg_growth_total_c(ic) + veg_growth_mt(ix_comp, ixC, ic)
        veg_growth_total_n(ic) = veg_growth_total_n(ic) + veg_growth_mt(ix_comp, ixN, ic)
        veg_growth_total_p(ic) = veg_growth_total_p(ic) + veg_growth_mt(ix_comp, ixP, ic)
        veg_growth_total_c13(ic) = veg_growth_total_c13(ic) + veg_growth_mt(ix_comp, ixC13, ic)
        veg_growth_total_c14(ic) = veg_growth_total_c14(ic) + veg_growth_mt(ix_comp, ixC14, ic)
        veg_growth_total_n15(ic) = veg_growth_total_n15(ic) + veg_growth_mt(ix_comp, ixN15, ic)

        veg_litterfall_total_c(ic) = veg_litterfall_total_c(ic) + veg_litterfall_mt(ix_comp, ixC, ic)
        veg_litterfall_total_n(ic) = veg_litterfall_total_n(ic) + veg_litterfall_mt(ix_comp, ixN, ic)
        veg_litterfall_total_p(ic) = veg_litterfall_total_p(ic) + veg_litterfall_mt(ix_comp, ixP, ic)
        veg_litterfall_total_c13(ic) = veg_litterfall_total_c13(ic) + veg_litterfall_mt(ix_comp, ixC13, ic)
        veg_litterfall_total_c14(ic) = veg_litterfall_total_c14(ic) + veg_litterfall_mt(ix_comp, ixC14, ic)
        veg_litterfall_total_n15(ic) = veg_litterfall_total_n15(ic) + veg_litterfall_mt(ix_comp, ixN15, ic)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END PARALLEL

    !> 6.0 biosphere-level diagnostics across multiple processes
    !>
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1)
    DO ic = 1, nc
      ! see also 'update_sb_pools()'
      !
      ! TODO calculation may include 'fFire' and 'fLUC' (land-use change) once available
      !
      ! could also be calculated based on NPP
      ! (n_processing_respiration = n_fixation_respiration + n_transform_respiration)
      net_biosphere_production(ic) = net_biosphere_production(ic) &
        & + gross_assimilation(ic) - maint_respiration(ic) - growth_respiration(ic) - n_processing_respiration(ic)

      ! biological N fixation (VEG_ and SB_ processes)
      biological_n_fixation(ic)    = biological_n_fixation(ic) + n_fixation(ic)

      IF (config_l_use_product_pools) THEN
        ! The product pool decay is a loss to the atmosphere (fprod_decay_mt unit is mol m-2 timestep-1)
        net_biosphere_production(ic) = net_biosphere_production(ic) - (fprod_decay_mt(ixC,ic) * 1000000.0_wp / dtime)
      ENDIF
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE update_veg_pools

#endif
END MODULE mo_q_veg_update_pools
