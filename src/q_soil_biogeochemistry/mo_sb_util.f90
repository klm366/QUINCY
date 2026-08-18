!> QUINCY helper routines for soil-beogeochemistry
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
!>#### various helper routines for the soil-beogeochemistry process
!>
MODULE mo_sb_util
#ifndef __NO_JSBACH__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: zero

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, &
    & SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: reset_sb_fluxes, calculate_time_average_soilbiogeochemistry, calc_qmax_texture

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_util'

CONTAINS

  ! ======================================================================================================= !
  !>reset sb fluxes (to zero)
  !>
  SUBROUTINE reset_sb_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: SB_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)    :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                !< the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_jsb_vgrid),      POINTER   :: vgrid_soil_w         !< Vertical grid
    INTEGER                           :: nsoil_w              !< number of soil layers as used/defined by the SB_ process
    INTEGER                           :: iblk, ics, ice, nc   !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_sb_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_transport_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_TRANSPORT_ID, sb_transport_mt)
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 bgcm fluxes
    !>
    ! set the bgcm matrices to zero
    sb_formation_mt(:,:,:,:)       = zero
    sb_loss_mt(:,:,:,:)            = zero
    sb_transport_mt(:,:,:,:)       = zero
    sb_mycorrhiza_export_mt(:,:,:) = zero

    !>2.0 variables
    !>
    ! 2D
    dsl4jsb_var2D_onChunk(SB_, total_flux_carbon)              = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_nitrogen)            = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_phosphorus)          = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_carbon13)            = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_carbon14)            = zero
    dsl4jsb_var2D_onChunk(SB_, total_flux_nitrogen15)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_nh4_solute)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_no3_solute)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_po4_solute)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_carbon)            = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_nitrogen)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_phosphorus)        = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_carbon13)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_carbon14)          = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_dom_nitrogen15)        = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_nh4_n15_solute)        = zero
    dsl4jsb_var2D_onChunk(SB_, leaching_no3_n15_solute)        = zero
    dsl4jsb_var2D_onChunk(SB_, emission_noy)                   = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2o)                   = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2)                    = zero
    dsl4jsb_var2D_onChunk(SB_, emission_noy_n15)               = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2o_n15)               = zero
    dsl4jsb_var2D_onChunk(SB_, emission_n2_n15)                = zero
    dsl4jsb_var2D_onChunk(SB_, fertiliser_nh4)                 = zero
    dsl4jsb_var2D_onChunk(SB_, fertiliser_nh4_n15)             = zero
    dsl4jsb_var2D_onChunk(SB_, fertiliser_no3)                 = zero
    dsl4jsb_var2D_onChunk(SB_, fertiliser_no3_n15)             = zero
    dsl4jsb_var2D_onChunk(SB_, fertiliser_po4)                 = zero
    ! 3D
    dsl4jsb_var3D_onChunk(SB_, het_respiration)                  = zero
    dsl4jsb_var3D_onChunk(SB_, het_respiration_c13)              = zero
    dsl4jsb_var3D_onChunk(SB_, het_respiration_c14)              = zero
    dsl4jsb_var3D_onChunk(SB_, myc_respiration)                  = zero
    dsl4jsb_var3D_onChunk(SB_, myc_respiration_c13)              = zero
    dsl4jsb_var3D_onChunk(SB_, myc_respiration_c14)              = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_nh4)           = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_nh4_n15)       = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_no3)           = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_no3_n15)       = zero
    dsl4jsb_var3D_onChunk(SB_, net_mineralisation_po4)           = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_nh4_sl)              = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_nh4_n15_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_no3_sl)              = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_no3_n15_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, plant_uptake_po4_sl)              = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_nh4_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_nh4_n15_sl)      = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_no3_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_no3_n15_sl)      = zero
    dsl4jsb_var3D_onChunk(SB_, microbial_uptake_po4_sl)          = zero
    dsl4jsb_var3D_onChunk(SB_, weathering_po4)                   = zero
    dsl4jsb_var3D_onChunk(SB_, occlusion_po4)                    = zero
    dsl4jsb_var3D_onChunk(SB_, slow_exchange_po4)                = zero
    dsl4jsb_var3D_onChunk(SB_, fast_exchange_po4)                = zero
    dsl4jsb_var3D_onChunk(SB_, fast_exchange_nh4)                = zero
    dsl4jsb_var3D_onChunk(SB_, fast_adsorpt_po4)                 = zero
    dsl4jsb_var3D_onChunk(SB_, fast_adsorpt_nh4)                 = zero
    dsl4jsb_var3D_onChunk(SB_, biochem_mineralisation_po4)       = zero
    dsl4jsb_var3D_onChunk(SB_, gross_mineralisation_po4)         = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_nh4_sl)         = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_nh4_n15_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_no3_sl)         = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_no3_n15_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_norg_sl)        = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_norg_n15_sl)    = zero
    dsl4jsb_var3D_onChunk(SB_, mycorrhiza_uptake_po4_sl)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_solute)             = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_assoc)              = zero
    dsl4jsb_var3D_onChunk(SB_, transport_no3_solute)             = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_solute)             = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_assoc_fast)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_assoc_slow)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_occluded)           = zero
    dsl4jsb_var3D_onChunk(SB_, transport_po4_primary)            = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_n15_solute)         = zero
    dsl4jsb_var3D_onChunk(SB_, transport_nh4_n15_assoc)          = zero
    dsl4jsb_var3D_onChunk(SB_, transport_no3_n15_solute)         = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_nh4_solute_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_no3_solute_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_po4_solute_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_carbon_sl)       = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_nitrogen_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_phosphorus_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_carbon13_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_carbon14_sl)     = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_dom_nitrogen15_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_nh4_n15_solute_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, lateral_loss_no3_n15_solute_sl)   = zero
    dsl4jsb_var3D_onChunk(SB_, transport_noy)                    = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2o)                    = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2)                     = zero
    dsl4jsb_var3D_onChunk(SB_, transport_noy_n15)                = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2o_n15)                = zero
    dsl4jsb_var3D_onChunk(SB_, transport_n2_n15)                 = zero
    dsl4jsb_var3D_onChunk(SB_, volatilisation_nh4)               = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_no3)                = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_noy)                = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_n2o)                = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_noy)              = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2o)              = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2)               = zero
    dsl4jsb_var3D_onChunk(SB_, asymb_n_fixation)                 = zero
    dsl4jsb_var3D_onChunk(SB_, volatilisation_nh4_n15)           = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_no3_n15)            = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_noy_n15)            = zero
    dsl4jsb_var3D_onChunk(SB_, nitrification_n2o_n15)            = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_noy_n15)          = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2o_n15)          = zero
    dsl4jsb_var3D_onChunk(SB_, denitrification_n2_n15)           = zero
    dsl4jsb_var3D_onChunk(SB_, asymb_n_fixation_n15)             = zero

  END SUBROUTINE reset_sb_fluxes

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for SB_
  !>
  SUBROUTINE calculate_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_time,                ONLY: get_year
    USE mo_lnd_time_averages        ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    dsl4jsb_Use_processes SB_, HYDRO_
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                !< instance of the model
    TYPE(t_lnd_bgcm_store), POINTER   :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER   :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER   :: vgrid_soil_w         !< Vertical grid
    INTEGER                           :: nsoil_w              !< number of soil layers as used/defined by the SB_ process
    REAL(wp)                          :: dtime                !< timestep length
    INTEGER                           :: iblk, ics, ice, nc   !< dimensions
    INTEGER                           :: current_year, next_jump_year, i_jump, jump_year
    LOGICAL                           :: collect_analytical_forcing
    CHARACTER(len=*), PARAMETER       :: routine = modname//':calculate_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! 3D SB_
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff
    dsl4jsb_Real3D_onChunk      :: dom_cn
    dsl4jsb_Real3D_onChunk      :: dom_cp
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_AP
    dsl4jsb_Real3D_onChunk      :: rtm_depolymerisation
    dsl4jsb_Real3D_onChunk      :: rmm_depolymerisation
    dsl4jsb_Real3D_onChunk      :: rtm_mic_uptake
    dsl4jsb_Real3D_onChunk      :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk      :: rtm_sorption
    dsl4jsb_Real3D_onChunk      :: rmm_sorption
    dsl4jsb_Real3D_onChunk      :: rtm_desorption
    dsl4jsb_Real3D_onChunk      :: rmm_desorption
    dsl4jsb_Real3D_onChunk      :: rtm_hsc
    dsl4jsb_Real3D_onChunk      :: rmm_hsc
    dsl4jsb_Real3D_onChunk      :: k_bioturb
    dsl4jsb_Real3D_onChunk      :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk      :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl
    dsl4jsb_Real2D_onChunk      :: num_sl_above_bedrock
    dsl4jsb_Real3D_onChunk      :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onChunk      :: myc_export_c_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_p_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: dom_cn_mavg
    dsl4jsb_Real3D_onChunk      :: dom_cp_mavg
    dsl4jsb_Real3D_onChunk      :: residue_som_c_form_mavg_sl
    dsl4jsb_Real3D_onChunk      :: residue_som_c14_form_mavg_sl
    dsl4jsb_Real3D_onChunk      :: residue_assoc_som_c_form_mavg_sl
    dsl4jsb_Real3D_onChunk      :: residue_assoc_som_c14_form_mavg_sl
    dsl4jsb_Real3D_onChunk      :: assoc_dom_c_form_mavg_sl
    dsl4jsb_Real3D_onChunk      :: assoc_dom_c14_form_mavg_sl
    dsl4jsb_Real3D_onChunk      :: residue_som_c_loss_mavg_sl
    dsl4jsb_Real3D_onChunk      :: residue_assoc_som_c_loss_mavg_sl
    dsl4jsb_Real3D_onChunk      :: assoc_dom_c_loss_mavg_sl
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onChunk      :: enzyme_frac_poly_p_mavg
    dsl4jsb_Real2D_onChunk      :: analy_ss_forcing_time
    dsl4jsb_Real3D_onChunk      :: analy_ss_input_sol_c_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_input_poly_c_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_input_wood_c_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_input_sol_c14_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_input_poly_c14_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_input_wood_c14_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_doc_transp_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_depoly_mod_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_uptake_mod_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_sorption_mod_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_desorption_mod_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_hsc_mod_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_enzyme_ap_mean_sl
    dsl4jsb_Real3D_onChunk      :: analy_ss_bioturb_numerator_mean_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN !< only if the present tile is a pft
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(SB_)
    IF (model%config%flag_analytical_ss_spinup_accelerator) THEN
      dsl4jsb_Get_memory(HYDRO_)
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    ! 3D SB_
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff)                  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn)                             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp)                             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p)                 ! in
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_cue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_nue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  microbial_pue_eff_tmic_mavg)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_c_tmyc_mavg_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tmyc_mavg_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_n_tlabile_mavg_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  myc_export_p_tlabile_mavg_sl)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cn_mavg)                        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  dom_cp_mavg)                        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_c_mavg)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_n_mavg)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_poly_p_mavg)            ! out
    IF (model%config%flag_analytical_ss_spinup_accelerator) THEN
      dsl4jsb_Get_var3D_onChunk(SB_,  enzyme_frac_AP)                    ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rtm_depolymerisation)              ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rmm_depolymerisation)              ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rtm_mic_uptake)                    ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rmm_mic_uptake)                    ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rtm_sorption)                      ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rmm_sorption)                      ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rtm_desorption)                    ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rmm_desorption)                    ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rtm_hsc)                           ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  rmm_hsc)                           ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  k_bioturb)                         ! in
      dsl4jsb_Get_var3D_onChunk(SB_,  bulk_dens_corr_sl)                 ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, frac_wtr_transp_down_sl)         ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)                   ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, num_sl_above_bedrock)            ! in
      dsl4jsb_Get_var2D_onChunk(SB_,  analy_ss_forcing_time)       ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_input_sol_c_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_input_poly_c_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_input_wood_c_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_input_sol_c14_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_input_poly_c14_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_input_wood_c14_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_doc_transp_mean_sl)   ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_depoly_mod_mean_sl)   ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_uptake_mod_mean_sl)   ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_sorption_mod_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_desorption_mod_mean_sl) ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_hsc_mod_mean_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_enzyme_ap_mean_sl)    ! inout
      dsl4jsb_Get_var3D_onChunk(SB_,  analy_ss_bioturb_numerator_mean_sl) ! inout
    END IF
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      dsl4jsb_Get_var3D_onChunk(SB_,  residue_som_c_form_mavg_sl)         ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  residue_som_c14_form_mavg_sl)       ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  residue_assoc_som_c_form_mavg_sl)   ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  residue_assoc_som_c14_form_mavg_sl) ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  assoc_dom_c_form_mavg_sl)           ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  assoc_dom_c14_form_mavg_sl)         ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  residue_som_c_loss_mavg_sl)         ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  residue_assoc_som_c_loss_mavg_sl)   ! out
      dsl4jsb_Get_var3D_onChunk(SB_,  assoc_dom_c_loss_mavg_sl)           ! out
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    IF (model%config%flag_slow_sb_pool_spinup_accelerator .OR. &
      & model%config%flag_analytical_ss_spinup_accelerator) THEN
      dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
    END IF
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)
    END IF
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 daytime averages - not used here
    !>

    !>1.0 moving averages
    !>
    ! docu:
    ! calc_time_mavg(dtime, current average, new value, length of avg_period,  !
    !                do_calc=LOGICAL, avg_period_unit='day')            ! OPTIONAL
    !                RETURN(new current average)
    ! the unit of the averaging period is 'day' by default, but can also be 'week' or 'year'

    !>  1.1 tlabile (averages at the timescale of the labile pool)
    !>
    myc_export_n_tlabile_mavg_sl(:,:) = calc_time_mavg(dtime, myc_export_n_tlabile_mavg_sl(:,:), &
      &                                                sb_mycorrhiza_export_mt(ixN,:,:), mavg_period_tlabile)
    myc_export_p_tlabile_mavg_sl(:,:) = calc_time_mavg(dtime, myc_export_p_tlabile_mavg_sl(:,:), &
      &                                                sb_mycorrhiza_export_mt(ixP,:,:), mavg_period_tlabile)

    !>  1.2 lctlib tau_mycorrhiza
    !>
    myc_export_c_tmyc_mavg_sl(:,:) = calc_time_mavg(dtime, myc_export_c_tmyc_mavg_sl(:,:),                &
      &                                             sb_mycorrhiza_export_mt(ixC,:,:),        &
                                                    lctlib%tau_mycorrhiza, avg_period_unit='year')
    myc_export_n_tmyc_mavg_sl(:,:) = calc_time_mavg(dtime, myc_export_n_tmyc_mavg_sl(:,:),                &
      &                                             sb_mycorrhiza_export_mt(ixN,:,:),        &
                                                    lctlib%tau_mycorrhiza, avg_period_unit='year')

    !>  1.3 tmic (microbial community acclimation)
    !>
    microbial_cue_eff_tmic_mavg(:,:) = calc_time_mavg(dtime, microbial_cue_eff_tmic_mavg(:,:), microbial_cue_eff(:,:), &
                                                      mavg_period_tmic)
    microbial_nue_eff_tmic_mavg(:,:) = calc_time_mavg(dtime, microbial_nue_eff_tmic_mavg(:,:), microbial_nue_eff(:,:), &
                                                      mavg_period_tmic)
    microbial_pue_eff_tmic_mavg(:,:) = calc_time_mavg(dtime, microbial_pue_eff_tmic_mavg(:,:), microbial_pue_eff(:,:), &
                                                      mavg_period_tmic)

    !>  1.4 tenzyme (memory time-scale for enzyme allocation)
    !>
    dom_cn_mavg(:,:)             = calc_time_mavg(dtime, dom_cn_mavg(:,:), dom_cn(:,:), &
                                                  mavg_period_tenzyme)
    dom_cp_mavg(:,:)             = calc_time_mavg(dtime, dom_cp_mavg(:,:), dom_cp(:,:), &
                                                  mavg_period_tenzyme)
    enzyme_frac_poly_c_mavg(:,:) = calc_time_mavg(dtime, enzyme_frac_poly_c_mavg(:,:), enzyme_frac_poly_c(:,:), &
                                                  mavg_period_tenzyme)
    enzyme_frac_poly_n_mavg(:,:) = calc_time_mavg(dtime, enzyme_frac_poly_n_mavg(:,:), enzyme_frac_poly_n(:,:), &
                                                  mavg_period_tenzyme)
    enzyme_frac_poly_p_mavg(:,:) = calc_time_mavg(dtime, enzyme_frac_poly_p_mavg(:,:), enzyme_frac_poly_p(:,:), &
                                                  mavg_period_tenzyme)

    !>  1.5 tresidual (memory time-scale for som spinup accelarator)
    !>
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      residue_som_c_form_mavg_sl(:,:)         = calc_time_mavg(dtime, residue_som_c_form_mavg_sl(:,:),        &
                                                                      sb_formation_mt(ix_residue, ixC, :, :), &
                                                                      mavg_period_tresidual)
      residue_som_c14_form_mavg_sl(:,:)       = calc_time_mavg(dtime, residue_som_c14_form_mavg_sl(:,:),        &
                                                                      sb_formation_mt(ix_residue, ixC14, :, :), &
                                                                      mavg_period_tresidual)
      residue_assoc_som_c_form_mavg_sl(:,:)   = calc_time_mavg(dtime, residue_assoc_som_c_form_mavg_sl(:,:),        &
                                                                      sb_formation_mt(ix_residue_assoc, ixC, :, :), &
                                                                      mavg_period_tresidual)
      residue_assoc_som_c14_form_mavg_sl(:,:) = calc_time_mavg(dtime, residue_assoc_som_c14_form_mavg_sl(:,:),        &
                                                                      sb_formation_mt(ix_residue_assoc, ixC14, :, :), &
                                                                      mavg_period_tresidual)
      assoc_dom_c_form_mavg_sl(:,:)           = calc_time_mavg(dtime, assoc_dom_c_form_mavg_sl(:,:),            &
                                                                      sb_formation_mt(ix_dom_assoc, ixC, :, :), &
                                                                      mavg_period_tresidual)
      assoc_dom_c14_form_mavg_sl(:,:)         = calc_time_mavg(dtime, assoc_dom_c14_form_mavg_sl(:,:),            &
                                                                      sb_formation_mt(ix_dom_assoc, ixC14, :, :), &
                                                                      mavg_period_tresidual)
      residue_som_c_loss_mavg_sl(:,:)         = calc_time_mavg(dtime, residue_som_c_loss_mavg_sl(:,:),   &
                                                                      sb_loss_mt(ix_residue, ixC, :, :), &
                                                                      mavg_period_tresidual)
      residue_assoc_som_c_loss_mavg_sl(:,:)   = calc_time_mavg(dtime, residue_assoc_som_c_loss_mavg_sl(:,:),   &
                                                                      sb_loss_mt(ix_residue_assoc, ixC, :, :), &
                                                                      mavg_period_tresidual)
      assoc_dom_c_loss_mavg_sl(:,:)           = calc_time_mavg(dtime, assoc_dom_c_loss_mavg_sl(:,:),       &
                                                                      sb_loss_mt(ix_dom_assoc, ixC, :, :), &
                                                                      mavg_period_tresidual)
    END IF

    !> 1.6 exact time-weighted means for the analytical steady-state accelerator
    !> Litter formation in the BGCM store is timestep-integrated; all analytical
    !> forcing means use rate units so that they are independent of model timestep.
    collect_analytical_forcing = .FALSE.
    IF (model%config%flag_analytical_ss_spinup_accelerator .AND. model%config%l_spinup .AND. dtime > 0.0_wp) THEN
      current_year = get_year(options%current_datetime)
      next_jump_year = HUGE(next_jump_year)
      DO i_jump = 0, model%config%analytical_ss_spinup_accelerator_max_executions - 1
        jump_year = model%config%analytical_ss_spinup_accelerator_start_year + &
          & i_jump * model%config%analytical_ss_spinup_accelerator_frequency
        IF (jump_year > current_year) THEN
          next_jump_year = jump_year
          EXIT
        END IF
      END DO
      IF (next_jump_year < HUGE(next_jump_year)) THEN
        collect_analytical_forcing = current_year >= next_jump_year - &
          & model%config%analytical_ss_spinup_accelerator_forcing_years
      END IF
    END IF

    IF (collect_analytical_forcing) THEN
      analy_ss_forcing_time(:) = analy_ss_forcing_time(:) + dtime
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & sb_formation_mt(ix_soluable_litter, ixC, :, :) / dtime, analy_ss_input_sol_c_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & sb_formation_mt(ix_polymeric_litter, ixC, :, :) / dtime, analy_ss_input_poly_c_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & sb_formation_mt(ix_woody_litter, ixC, :, :) / dtime, analy_ss_input_wood_c_mean_sl)
      IF (model%config%include_carbon14) THEN
        CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
          & sb_formation_mt(ix_soluable_litter, ixC14, :, :) / dtime, analy_ss_input_sol_c14_mean_sl)
        CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
          & sb_formation_mt(ix_polymeric_litter, ixC14, :, :) / dtime, analy_ss_input_poly_c14_mean_sl)
        CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
          & sb_formation_mt(ix_woody_litter, ixC14, :, :) / dtime, analy_ss_input_wood_c14_mean_sl)
      END IF
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & frac_wtr_transp_down_sl, analy_ss_doc_transp_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & rtm_depolymerisation * rmm_depolymerisation, analy_ss_depoly_mod_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & rtm_mic_uptake * rmm_mic_uptake, analy_ss_uptake_mod_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & rtm_sorption * rmm_sorption, analy_ss_sorption_mod_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & rtm_desorption * rmm_desorption, analy_ss_desorption_mod_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & rtm_hsc * rmm_hsc, analy_ss_hsc_mod_mean_sl)
      CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
        & enzyme_frac_AP, analy_ss_enzyme_ap_mean_sl)
      IF (model%config%flag_analytical_ss_bioturbation) THEN
        CALL update_analytical_forcing_mean(dtime, analy_ss_forcing_time, &
          & k_bioturb * bulk_dens_corr_sl, analy_ss_bioturb_numerator_mean_sl)
      END IF
    END IF
  END SUBROUTINE calculate_time_average_soilbiogeochemistry

  ! ======================================================================================================= !
  !> update a layer-wise running mean using the elapsed averaging time
  SUBROUTINE update_analytical_forcing_mean(dtime, elapsed_time, sample, running_mean)

    REAL(wp), INTENT(in)    :: dtime
    REAL(wp), INTENT(in)    :: elapsed_time(:)
    REAL(wp), INTENT(in)    :: sample(:,:)
    REAL(wp), INTENT(inout) :: running_mean(:,:)

    REAL(wp) :: weight(SIZE(elapsed_time))
    INTEGER  :: is

    weight(:) = dtime / MAX(elapsed_time(:), dtime)
    DO is = 1, SIZE(sample, 2)
      running_mean(:,is) = running_mean(:,is) + weight(:) * (sample(:,is) - running_mean(:,is))
    END DO

  END SUBROUTINE update_analytical_forcing_mean

  ! ======================================================================================================= !
  !>calculate sorption capacity for organic material
  !>
  ELEMENTAL FUNCTION calc_qmax_texture(qmax_fine_particle, silt_sl, clay_sl) RESULT(qmax_texture)

    REAL(wp), INTENT(in) :: qmax_fine_particle    !< maximum sorption capacity of fine soil particle [mol / kg fine particle (silt or clay)]
    REAL(wp), INTENT(in) :: silt_sl               !< silt content of the mineral soil [kg silt / kg mineral soil]
    REAL(wp), INTENT(in) :: clay_sl               !< clay content of the mineral soil [kg clay / kg mineral soil]
    REAL(wp)             :: qmax_texture          !< maximum sorption capacity of the soil [mol / kg]

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_qmax_texture'

    ! maximum sorption capacity per kg mineral soil (mol/kg)
    qmax_texture = qmax_fine_particle * (silt_sl + clay_sl)

  END FUNCTION calc_qmax_texture


#endif
END MODULE mo_sb_util
