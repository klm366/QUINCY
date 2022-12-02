!>
!!  @par Copyright
!!    QUINCY terrestrial biosphere model with coupled carbon, nitrogen, and phosphorus cycles
!!    Copyright (C) 2019 Max Planck Institute for Biogeochemistry (MPI-BGC), Jena, Germany 
!!    Contributors: Caldararu S., Engel J., Kern M., Pallandt M., Schnur R., Thum T., Yu L., Zaehle S.
!!
!!    This file is part of the QUINCY model.
!!    The QUINCY model is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License version 3,
!!    as published by the Free Software Foundation.
!!
!!    The QUINCY model is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with the QUINCY model.  If not, see <https://www.gnu.org/licenses/>.
!! 
!! @brief Contains the routines for the plant growth 
!!
!! @author
!!  Soenke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Soenke Zaehle (2015-12-04)
!!
MODULE mo_veg_growth
#ifndef __NO_JSBACH__

#ifdef __QUINCY_STANDALONE__
#include "dsl_jsbach.inc"
#endif

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message
  USE mo_exception,             ONLY: message_text, finish


  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_veg_growth
  PUBLIC :: calc_diameter_from_woody_biomass, calc_height_from_diameter

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_growth'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update plant growth (flux\%growth)
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_veg_growth(tile, options)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    ! USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SSE_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_process_class,       ONLY: A2L_, L2A_, SEB_, TURB_, SPQ_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_veg_respiration,         ONLY: temperature_response_respiration, calc_maintenance_respiration 
    USE mo_veg_canopy,              ONLY: calc_dir_optimal_cn_leaf, calc_marginal_canopy_flux_increment
    USE mo_veg_constants,           ONLY: fmaint_rate_base

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(ASSIMI_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SB_)
    
    ! Use of process memories
    dsl4jsb_Use_memory(A2L_)        !  USE mo_a2l_memory_class,      ONLY: t_a2l_memory
    dsl4jsb_Use_memory(ASSIMI_)     !  USE mo_assimi_memory_class,   ONLY: t_assimi_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(PHENO_)      !  USE mo_pheno_memory_class,    ONLY: t_pheno_memory
    dsl4jsb_Use_memory(RAD_)        !  USE mo_rad_memory_class,      ONLY: t_rad_memory
    ! dsl4jsb_Use_memory(SSE_)       !  USE mo_sse_memory_class,     ONLY: t_sse_memory
    ! dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    ! Use of pools
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_veg
    USE mo_veg_memory_class,        ONLY: veg_pool, veg_growth, veg_falloc, veg_exudation, veg_reserve_use
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,         ONLY: sb_pool, sb_mycorrhiza_export
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool, t_veg_elements_pool
    dsl4jsb_Use_pools(SB_) t_sb_components_pool, t_sb_elements_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER         :: model            !< model 
    TYPE(t_lctlib_element), POINTER         :: lctlib           !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_canopy     !< Vertical grid
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_soil_sb    !< Vertical grid
    INTEGER                                 :: ncanopy          !< number of canopy layers
    INTEGER                                 :: nsoil_sb         !< number of soil layers as used/defined by the SB_ process
    REAL(wp), DIMENSION(options%nc)         :: dheight, target_labile_pool_carbon, target_labile_pool_nitrogen, &
                                               target_labile_pool_phosphorus, target_reserve_pool_carbon, &
                                               target_reserve_pool_nitrogen, target_reserve_pool_phosphorus
    REAL(wp), DIMENSION(options%nc)         :: hlp1
    REAL(wp)                                :: lctlib_g1              !< set to g1_medlyn or g1_bberry depending on canopy_conductance_scheme
    REAL(wp), ALLOCATABLE                   :: canopy_layer_depth(:)  !< thickness of each single canopy layer [m2 m-2 LAI]
    INTEGER                                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_growth'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration, process memory and pools
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(A2L_)                        !     TYPE(t_a2l_memory),      POINTER :: a2l__mem
    dsl4jsb_Def_memory(ASSIMI_)                     !     TYPE(t_assimi_memory),   POINTER :: assimi__mem
    dsl4jsb_Def_memory(VEG_)                        !     TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(PHENO_)                      !     TYPE(t_pheno_memory),    POINTER :: pheno__mem
    dsl4jsb_Def_memory(RAD_)                        !     TYPE(t_rad_memory),      POINTER :: rad__mem
    ! dsl4jsb_Def_memory(SSE_)                       !     TYPE(t_sse_memory),     POINTER :: sse__mem
    ! dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(SB_)                         !     TYPE(t_sb_memory),       POINTER :: sb__mem
#ifdef __QUINCY_STANDALONE__
#else
    dsl4jsb_Def_pool(t_veg_tissues_pool)    :: veg_pool
    dsl4jsb_Def_pool(t_veg_tissues_pool)    :: veg_growth
    dsl4jsb_Def_pool(t_veg_tissues_pool)    :: veg_falloc
    dsl4jsb_Def_pool(t_veg_elements_pool)   :: veg_exudation
    dsl4jsb_Def_pool(t_veg_elements_pool)   :: veg_reserve_use
    dsl4jsb_Def_pool(t_sb_components_pool)  :: sb_pool
    dsl4jsb_Def_pool(t_sb_elements_pool)    :: sb_mycorrhiza_export
#endif
    ! Declare pointers to variables in memory
    ! PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: lai
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: q_air
    dsl4jsb_Real2D_onChunk      :: press_srf
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio
    dsl4jsb_Real2D_onChunk      :: ga
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: swpar_fdiffuse
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_air
    dsl4jsb_Real2D_onChunk      :: beta_soa
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: maint_respiration_leaf
    dsl4jsb_Real2D_onChunk      :: net_assimilation
    dsl4jsb_Real2D_onChunk      :: canopy_cond
    dsl4jsb_Real2D_onChunk      :: beta_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tcnl_mavg
    ! SB_
    dsl4jsb_Real3D_onChunk      :: nh4_solute
    dsl4jsb_Real3D_onChunk      :: no3_solute
    dsl4jsb_Real3D_onChunk      :: plant_uptake_nh4_sl 
    dsl4jsb_Real3D_onChunk      :: plant_uptake_nh4_n15_sl 
    dsl4jsb_Real3D_onChunk      :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk      :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk      :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_p_tlabile_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_c_tmyc_mavg_sl
    dsl4jsb_Real3D_onChunk      :: myc_export_n_tmyc_mavg_sl
    ! ! SSE_ 2D
    ! dsl4jsb_Real2D_onChunk      :: t_soil_root
    ! ! HYDRO_ 2D
    ! dsl4jsb_Real2D_onChunk      :: w_soil_root_theta
    ! ! HYDRO_ 3D
    ! dsl4jsb_Real3D_onChunk      :: soil_depth_sl 
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: t_soil_root
    dsl4jsb_Real2D_onChunk      :: w_soil_root_theta
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk      :: soil_depth_sl 
    ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tcnl_mavg_cl
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: target_cn_leaf
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: press_srf_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: co2_mixing_ratio_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: ga_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: t_soil_root_tacclim_mavg
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: leaf_cn_direction
    dsl4jsb_Real2D_onChunk      :: beta_sinklim_ps
    dsl4jsb_Real2D_onChunk      :: unit_npp
    dsl4jsb_Real2D_onChunk      :: unit_transpiration
    dsl4jsb_Real2D_onChunk      :: kstar_labile
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_root
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_pot_troot_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_pot
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_act
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_act
    dsl4jsb_Real2D_onChunk      :: cost_n_uptake_root
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation
    dsl4jsb_Real2D_onChunk      :: target_np_leaf
    dsl4jsb_Real2D_onChunk      :: target_cn_fine_root
    dsl4jsb_Real2D_onChunk      :: target_cn_coarse_root
    dsl4jsb_Real2D_onChunk      :: target_cn_sap_wood
    dsl4jsb_Real2D_onChunk      :: target_cn_fruit
    dsl4jsb_Real2D_onChunk      :: target_np_fine_root
    dsl4jsb_Real2D_onChunk      :: target_np_coarse_root
    dsl4jsb_Real2D_onChunk      :: target_np_sap_wood
    dsl4jsb_Real2D_onChunk      :: target_np_fruit
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: uptake_n_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: growth_np_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: npp_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: fmaint_rate_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: labile_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: labile_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_carbon_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_nitrogen_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: reserve_phosphorus_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: w_root_lim_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cn_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_cp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: growth_np_talloc_mavg 
    dsl4jsb_Real2D_onChunk      :: npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_npp_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_transpiration_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_n_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: unit_uptake_p_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: dphi_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: leaf2sapwood_mass_ratio 
    dsl4jsb_Real2D_onChunk      :: leaf2root_mass_ratio
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: leaf2root_troot_mavg
    dsl4jsb_Real2D_onChunk      :: target_lai
    dsl4jsb_Real2D_onChunk      :: k1_opt
    dsl4jsb_Real2D_onChunk      :: k2_opt
    dsl4jsb_Real2D_onChunk      :: growth_req_p
    dsl4jsb_Real2D_onChunk      :: growth_req_n
    dsl4jsb_Real2D_onChunk      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: demand_uptake_p_tuptake_mavg
    dsl4jsb_Real2D_onChunk      :: cost_n_fixation_talloc_mavg
    dsl4jsb_Real2D_onChunk      :: exudation_c_tmyc_mavg
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_nh4_n15
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: uptake_no3_n15
    dsl4jsb_Real2D_onChunk      :: uptake_po4
    dsl4jsb_Real2D_onChunk      :: f_n_demand
    dsl4jsb_Real2D_onChunk      :: f_p_demand
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: n_fixation_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c13
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c14
    dsl4jsb_Real2D_onChunk      :: maint_respiration_pot
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tcnl_mavg_cl
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
#ifdef __QUINCY_STANDALONE__
    vgrid_canopy  => Get_vgrid('canopy_layer')
#else
    vgrid_canopy  => Get_vgrid('canopy_layer_q_')
#endif
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')    
    ncanopy       =  vgrid_canopy%n_levels
    nsoil_sb      =  vgrid_soil_sb%n_levels           ! also needed for sb_pool & sb_mycorrhiza_export
    ALLOCATE(canopy_layer_depth(ncanopy))
    canopy_layer_depth(:) = vgrid_canopy%dz(:)
    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(RAD_)
    ! dsl4jsb_Get_memory(SSE_)
    ! dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(SB_)
#ifdef __QUINCY_STANDALONE__
#else
    ! Get pools
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_pool,            veg_pool)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_growth,          veg_growth)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_falloc,          veg_falloc)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_exudation,       veg_exudation)
    dsl4jsb_Get_pool(VEG_, veg_main_pool%veg_reserve_use,     veg_reserve_use)
    dsl4jsb_Get_pool(SB_,  sb_main_pool%sb_pool,              sb_pool)
    dsl4jsb_Get_pool(SB_,  sb_main_pool%sb_mycorrhiza_export, sb_mycorrhiza_export)
    ! associate element-variables with chunk
    CALL veg_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_growth%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_falloc%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_exudation%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL veg_reserve_use%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL sb_pool%Associate_chunk_var_ptr(ics, ice, iblk)
    CALL sb_mycorrhiza_export%Associate_chunk_var_ptr(ics, ice, iblk)
#endif
    ! Get process variables (Set pointers to variables in memory)
    ! PHENO_ 
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)                   ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                              ! in
    ! A2L_ 
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)                              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)                          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, co2_mixing_ratio)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, ga)                                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)                     ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_fdiffuse)                     ! in
    ! ASSIMI_ 
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air)                        ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soa)                        ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps)                    ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs)                    ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, maint_respiration_leaf)          ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, net_assimilation)                ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, canopy_cond)                     ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_air_tcnl_mavg)              ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_ps_tcnl_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tcnl_mavg)          ! in
    ! SB_ 3D
    dsl4jsb_Get_var3D_onChunk(SB_, nh4_solute)                          ! in
    dsl4jsb_Get_var3D_onChunk(SB_, no3_solute)                          ! in
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_nh4_sl)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_nh4_n15_sl)             ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_no3_sl)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_no3_n15_sl)             ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, plant_uptake_po4_sl)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, myc_export_n_tlabile_mavg_sl)        ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, myc_export_p_tlabile_mavg_sl)        ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, myc_export_c_tmyc_mavg_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_, myc_export_n_tmyc_mavg_sl)           ! inout
    ! ! SSE_ 2D
    ! dsl4jsb_Get_var2D_onChunk(SSE_, t_soil_root)                        ! in
    ! ! HYDRO_ 2D
    ! dsl4jsb_Get_var2D_onChunk(HYDRO_, w_soil_root_theta)                ! in
    ! ! HYDRO_ 3D
    ! dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)                    ! in
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_, t_soil_root)                        ! in
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_theta)                  ! in
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_, soil_depth_sl)                      ! in
    ! RAD_ 3D
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)                     ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)                     ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tcnl_mavg_cl)           ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tcnl_mavg_cl)           ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_leaf)                     ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tcnl_mavg)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, press_srf_tcnl_mavg)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, ga_tcnl_mavg)                       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tacclim_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_soil_root_tacclim_mavg)           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps_tcnl_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf_cn_direction)                  ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, beta_sinklim_ps)                    ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp)                           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, kstar_labile)                       ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_root)                   ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_troot_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_troot_mavg)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot_troot_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_pot)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_pot)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_act)                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_act)                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_uptake_root)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation)                    ! out 
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_leaf)                     ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fine_root)                ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_coarse_root)              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_sap_wood)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, target_cn_fruit)                    ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fine_root)                ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_coarse_root)              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_sap_wood)                 ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, target_np_fruit)                    ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                           ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_n_tcnl_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_tcnl_mavg)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_np_tcnl_mavg)                ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_tcnl_mavg)                      ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, fmaint_rate_tcnl_mavg)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_tcnl_mavg)            ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_tcnl_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_tcnl_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_carbon_talloc_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_nitrogen_talloc_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, labile_phosphorus_talloc_mavg)      ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_carbon_talloc_mavg)         ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_nitrogen_talloc_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, reserve_phosphorus_talloc_mavg)     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, w_root_lim_talloc_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cn_talloc_mavg)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_cp_talloc_mavg)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_np_talloc_mavg)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_talloc_mavg)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_npp_talloc_mavg)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, transpiration_talloc_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration_talloc_mavg)     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_talloc_mavg)             ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_n_talloc_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_uptake_p_talloc_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, dphi_talloc_mavg)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2sapwood_mass_ratio)            ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_mass_ratio)               ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)          ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, leaf2root_troot_mavg)               ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, target_lai)                         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, k1_opt)                             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, k2_opt)                             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p)                       ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n)                       ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_n_tuptake_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, demand_uptake_p_tuptake_mavg)       ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, cost_n_fixation_talloc_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, exudation_c_tmyc_mavg)              ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)                         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4_n15)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)                         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3_n15)                     ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_po4)                         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, f_n_demand)                         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, f_p_demand)                         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)                         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)            ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_respiration)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration)           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)              ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)              ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c13)       ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c14)       ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_pot)              ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)                  ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)                 ! out
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                          ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)                           ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)                          ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)                         ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)                          ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)                    ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl)          ! in
#ifdef __QUINCY_STANDALONE__
    ! Get pools
#define Get_veg_pool
#define Get_veg_growth
#define Get_veg_falloc
#define Get_veg_exudation
#define Get_veg_reserve_use
#define Get_sb_pool
#define Get_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Get_veg_pool
#undef Get_veg_growth
#undef Get_veg_falloc
#undef Get_veg_exudation
#undef Get_veg_reserve_use
#undef Get_sb_pool
#undef Get_sb_mycorrhiza_export
#endif


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 0.9 set g1 according to canopy_conductance_scheme
    !! 
    SELECT CASE(TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme))
      CASE ("medlyn")
        lctlib_g1 = lctlib%g1_medlyn
      CASE ("ballberry")
        lctlib_g1 = lctlib%g1_bberry
    END SELECT

    !> 1.0 calculate potential growth rate relative to labile pool size
    !! 
    !! 
    kstar_labile(:) = calc_meristem_activity(t_air(:), w_soil_root_theta(:), growing_season(:))

    !> 1.1 calculate N-based maintenance respiration rate (as potential maitenance rate) and 
    !!     update N-specific maintenance respiration rate 
    !!
    maint_respiration_pot(:) = maint_respiration_leaf(:) + calc_maintenance_respiration( &
                               veg_pool%leaf%carbon(:)                     , &
                               t_air(:), t_air_tacclim_mavg(:)             , & 
                               t_soil_root(:), t_soil_root_tacclim_mavg(:) , &  
                               veg_pool%sap_wood%nitrogen(:)               , & 
                               veg_pool%fine_root%nitrogen(:)              , &
                               veg_pool%coarse_root%nitrogen(:)            )
                                
    fmaint_rate_root(:) = fmaint_rate_base * temperature_response_respiration(t_soil_root(:), t_soil_root_tacclim_mavg(:))

    !> 1.2 calculate utrient demand scalars
    !! 
    CALL calc_nutrient_demand_scalar( &
      & lctlib%cn_leaf                    , & ! lctlib
      & lctlib%np_leaf                    , & ! lctlib
      & leaf2root_troot_mavg(:)           , & ! in
      & veg_pool%labile%carbon(:)         , & ! in
      & veg_pool%labile%nitrogen(:)       , & ! in
      & veg_pool%labile%phosphorus(:)     , & ! in
      & f_n_demand(:)                     , & ! out
      & f_p_demand(:)                     )   ! out

    !> 1.3 calculate actual nutrient uptake and N fixation, as well as associated respiration
    !!
    CALL calc_actual_nutrient_uptake_rate( & 
      & nc, &                                             ! dimension
      & nsoil_sb, &                                       ! dimension
      & lctlib%cn_leaf, &                                 ! lctlib
      & lctlib%np_leaf, &                                 ! lctlib
      & lctlib%tau_fine_root, &                           ! lctlib
      & TRIM(model%config%quincy_model), &                ! config
      & TRIM(dsl4jsb_Config(VEG_)%bnf_scheme), &          ! config
      & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &      ! config
      & dsl4jsb_Config(SB_)%flag_mycorrhiza, &            ! config
      & kstar_labile(:), &                                ! in
      & t_air(:), &
      & veg_pool%fine_root%carbon(:), & 
      & veg_pool%fine_root%nitrogen(:), & 
      & veg_pool%labile%carbon(:), & 
      & veg_pool%labile%nitrogen(:), & 
      & veg_pool%labile%phosphorus(:), & 
      & sb_pool%mycorrhiza%carbon(:,:), & 
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & soil_depth_sl(:,:), &
      & growth_req_n_tlabile_mavg(:), &
      & demand_uptake_n_tuptake_mavg(:), &
      & demand_uptake_p_tuptake_mavg(:), &
      & fmaint_rate_troot_mavg(:), &
      & unit_npp_troot_mavg(:), &
      & unit_uptake_n_pot_troot_mavg(:), &
      & cost_n_fixation_talloc_mavg(:), &
      & exudation_c_tmyc_mavg(:), &
      & unit_uptake_n_pot(:), &
      & unit_uptake_p_pot(:), &                           ! in
      & plant_uptake_nh4_sl(:,:), &                       ! inout
      & plant_uptake_nh4_n15_sl(:,:), &                   ! inout
      & plant_uptake_no3_sl(:,:), &                       ! inout
      & plant_uptake_no3_n15_sl(:,:), &                   ! inout
      & plant_uptake_po4_sl(:,:), &                       ! inout
      & sb_mycorrhiza_export, &                           ! inout
      & myc_export_n_tlabile_mavg_sl(:,:), &              ! inout
      & myc_export_p_tlabile_mavg_sl(:,:), &              ! inout
      & myc_export_c_tmyc_mavg_sl(:,:), &                 ! inout
      & myc_export_n_tmyc_mavg_sl(:,:), &                 ! inout
      & uptake_nh4(:), &                                  ! out
      & uptake_nh4_n15(:), &                              ! out
      & uptake_no3(:), &                                  ! out
      & uptake_no3_n15(:), &                              ! out
      & uptake_po4(:), &                                  ! out
      & n_fixation(:), &                                  ! out
      & unit_uptake_n_act(:), &                           ! out
      & unit_uptake_p_act(:), &                           ! out
      & cost_n_uptake_root(:), &                          ! out
      & cost_n_fixation(:), &                             ! out
      & n_transform_respiration(:), &                     ! out
      & n_fixation_respiration(:), &                      ! out
      & n_processing_respiration(:), &                    ! out
      & veg_exudation%carbon(:) )                         ! out


    !> 2.0 Calculate target C:N:P for each tissue type
    !! 
    !! this function loops over:                    \n
    !!   points:          DO ipoint=1,nc            \n
    !!   direction:       DO dd = 1,3               \n
    !!   canopy layers:   DO icanopy =1, ncanopy
    !! 
    IF (TRIM(dsl4jsb_Config(VEG_)%leaf_stoichom_scheme) == "optimal") THEN
      CALL calc_dir_optimal_cn_leaf(    nc                                                        , & ! dimension
                                        ncanopy                                                   , &
                                        vgrid_canopy%dz(:)                                        , &
                                        vgrid_canopy%lbounds(:)                                   , &
                                        vgrid_canopy%ubounds(:)                                   , & ! dimension 
                                        lctlib%ps_pathway                                         , & ! lctlib
                                        lctlib%k0_fn_struc                                        , &
                                        lctlib%fn_oth_min                                         , &
                                        lctlib%sla                                                , &
                                        lctlib%np_leaf                                            , &
                                        lctlib%gmin                                               , &
                                        lctlib%g0                                                 , &
                                        lctlib_g1                                                 , &
                                        lctlib%t_jmax_omega                                       , &
                                        lctlib%tau_leaf                                           , & ! lctlib
                                        dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction            , & ! ASSIMI_ config
                                        TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme)   , & ! ASSIMI_ config
                                        canopy_layer_depth(:)                                     , & ! in
                                        growing_season(:)                                         , &
                                        veg_pool%leaf%carbon(:)                                   , &  
                                        veg_pool%leaf%nitrogen(:)                                 , &
                                        veg_pool%fine_root%carbon(:)                              , &
                                        lai(:)                                                    , &
                                        t_jmax_opt_mavg(:)                                        , &
                                        target_cn_leaf(:)                                         , &
                                        target_cn_fine_root(:)                                    , &
                                        t_air_tacclim_mavg(:)                                     , &
                                        t_air_tcnl_mavg(:)                                        , &
                                        ga_tcnl_mavg(:)                                           , &
                                        press_srf_tcnl_mavg(:)                                    , &
                                        co2_mixing_ratio_tcnl_mavg(:)                             , &
                                        beta_sinklim_ps_tcnl_mavg(:)                              , &
                                        beta_air_tcnl_mavg(:)                                     , &
                                        beta_soa(:)                                               , &
                                        beta_soil_ps_tcnl_mavg(:)                                 , &
                                        beta_soil_gs_tcnl_mavg(:)                                 , &
                                        fn_chl_cl(:,:)                                            , &
                                        fn_et_cl(:,:)                                             , &
                                        fn_rub_cl(:,:)                                            , &
                                        fn_pepc_cl(:,:)                                           , &
                                        fn_oth_cl(:,:)                                            , &
                                        fleaf_sunlit_tcnl_mavg_cl(:,:)                            , &
                                        ppfd_sunlit_tcnl_mavg_cl(:,:)                             , &
                                        ppfd_shaded_tcnl_mavg_cl(:,:)                             , &
                                        uptake_n_tcnl_mavg(:)                                     , &
                                        growth_cn_tcnl_mavg(:)                                    , &
                                        npp_tcnl_mavg(:)                                          , &
                                        fmaint_rate_tcnl_mavg(:)                                  , &
                                        labile_carbon_tcnl_mavg(:)                                , &
                                        labile_nitrogen_tcnl_mavg(:)                              , & ! in
                                        leaf_cn_direction(:)                                      )   ! inout
    END IF 

    CALL calc_stoichiometry_changes(  lctlib%growthform                               , & ! lctlib
                                      lctlib%phenology_type                           , &
                                      lctlib%cn_leaf                                  , &
                                      lctlib%np_leaf                                  , &
                                      lctlib%cn_leaf_min                              , &
                                      lctlib%cn_leaf_max                              , &
                                      lctlib%np_leaf_min                              , &
                                      lctlib%np_leaf_max                              , & ! lctlib
                                      TRIM(dsl4jsb_Config(VEG_)%leaf_stoichom_scheme) , & ! config
                                      growing_season(:)                               , & ! in 
                                      leaf_cn_direction(:)                            , &
                                      labile_carbon_tcnl_mavg(:)                      , &
                                      labile_nitrogen_tcnl_mavg(:)                    , &
                                      labile_phosphorus_tcnl_mavg(:)                  , &
                                      growth_cn_tcnl_mavg(:)                          , &
                                      growth_np_tcnl_mavg(:)                          , &
                                      lai(:)                                          , &
                                      target_lai(:)                                   , & ! in
                                      target_cn_leaf(:)                               , & ! inout
                                      target_np_leaf(:)                               , & ! inout
                                      target_cn_fine_root(:)                          , & ! out
                                      target_cn_coarse_root(:)                        , & ! out
                                      target_cn_sap_wood(:)                           , & ! out
                                      target_cn_fruit(:)                              , & ! out
                                      target_np_fine_root(:)                          , & ! out
                                      target_np_coarse_root(:)                        , & ! out
                                      target_np_sap_wood(:)                           , & ! out
                                      target_np_fruit(:)                              )   ! out
             

    !> 3.0 calculate current allometry 
    !! 
    !! 
    !> 3.1 calculate marginal NPP increment per unit C growth, needed for optimal allocation scheme
    !!
    CALL calc_marginal_canopy_flux_increment( nc                                                      , & ! dimension
                                              ncanopy                                                 , &
                                              vgrid_canopy%dz(:)                                      , &
                                              vgrid_canopy%lbounds(:)                                 , &
                                              vgrid_canopy%ubounds(:)                                 , & ! dimension
                                              lctlib%ps_pathway                                       , & ! lctlib
                                              lctlib%sla                                              , &
                                              lctlib%k0_fn_struc                                      , &
                                              lctlib%fn_oth_min                                       , &
                                              lctlib%np_leaf                                          , &
                                              lctlib%cn_leaf                                          , &
                                              lctlib%gmin                                             , &
                                              lctlib%g0                                               , &
                                              lctlib_g1                                               , &
                                              lctlib%t_jmax_omega                                     , &
                                              lctlib%sigma_vis                                        , & ! lctlib
                                              dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction          , & ! ASSIMI_ config
                                              TRIM(dsl4jsb_Config(ASSIMI_)%canopy_conductance_scheme) , & ! ASSIMI_ config
                                              canopy_layer_depth(:)                                   , & ! in
                                              t_air(:)                                                , &
                                              press_srf(:)                                            , &
                                              q_air(:)                                                , &
                                              co2_mixing_ratio(:)                                     , &
                                              ga(:)                                                   , & 
                                              net_assimilation(:)                                     , & 
                                              canopy_cond(:)                                          , & 
                                              beta_air(:)                                             , &
                                              beta_soa(:)                                             , & 
                                              beta_soil_ps(:)                                         , & 
                                              beta_soil_gs(:)                                         , & 
                                              veg_pool%leaf%carbon(:)                                 , &  
                                              veg_pool%leaf%nitrogen(:)                               , &
                                              lai(:)                                                  , &
                                              beta_sinklim_ps(:)                                      , &
                                              t_air_tacclim_mavg(:)                                   , &
                                              t_jmax_opt_mavg(:)                                      , &
                                              fn_chl_cl(:,:)                                          , &
                                              fn_et_cl(:,:)                                           , &
                                              fn_rub_cl(:,:)                                          , &
                                              fn_pepc_cl(:,:)                                         , &
                                              fn_oth_cl(:,:)                                          , &
                                              fleaf_sunlit_cl(:,:)                                    , &
                                              ppfd_sunlit_cl(:,:)                                     , &
                                              ppfd_shaded_cl(:,:)                                     , &
                                              swpar_srf_down(:)                                       , &
                                              swpar_fdiffuse(:)                                       , & ! in
                                              unit_npp(:)                                             , & ! out
                                              unit_transpiration(:)                                   )   ! out 
    
    !> 3.2 calculate allometry    
    !!
    CALL calc_allometry(lctlib%growthform                                 , & ! lctlib
                        lctlib%allom_k1                                   , &
                        lctlib%allom_k2                                   , &
                        lctlib%sla                                        , &
                        lctlib%wood_density                               , &
                        lctlib%tau_leaf                                   , &
                        lctlib%tau_fine_root                              , &
                        lctlib%k_latosa                                   , &
                        lctlib%k_root                                     , &
                        lctlib%k_rtos                                     , &
                        lctlib%c0_allom                                   , & ! lctlib
                        TRIM(dsl4jsb_Config(VEG_)%biomass_alloc_scheme)   , & ! config
                        TRIM(dsl4jsb_Config(VEG_)%bnf_scheme)             , & ! config
                        height(:)                                         , & ! in
                        veg_pool%sap_wood%carbon(:)                       , &
                        veg_pool%heart_wood%carbon(:)                     , &
                        dens_ind(:)                                       , & 
                        labile_carbon_talloc_mavg(:)                      , &
                        labile_nitrogen_talloc_mavg(:)                    , & 
                        labile_phosphorus_talloc_mavg(:)                  , &
                        reserve_carbon_talloc_mavg(:)                     , &
                        reserve_nitrogen_talloc_mavg(:)                   , &
                        reserve_phosphorus_talloc_mavg(:)                 , &
                        growth_cn_talloc_mavg(:)                          , &
                        growth_cp_talloc_mavg(:)                          , &
                        growth_np_talloc_mavg(:)                          , & 
                        npp_talloc_mavg(:)                                , &
                        unit_npp_talloc_mavg(:)                           , &
                        transpiration_talloc_mavg(:)                      , &
                        unit_transpiration_talloc_mavg(:)                 , &
                        n_fixation_talloc_mavg(:)                         , &
                        unit_uptake_n_talloc_mavg(:)                      , &
                        unit_uptake_p_talloc_mavg(:)                      , &
                        dphi_talloc_mavg(:)                               , &
                        kstar_labile(:)                                   , &
                        w_root_lim_talloc_mavg                            , & ! in
                        dheight(:)                                        , & ! out  (local variable)
                        leaf2sapwood_mass_ratio(:)                        , & ! out
                        leaf2root_mass_ratio(:)                           , & ! out
                        k1_opt(:)                                         , & ! out
                        k2_opt(:)                                         )   ! out

   
    !> 4.0 calculate reserve & labile pool dynamics
    !! 
    !> 4.1 calculate desired sizes of LAI, reserve and labile pool 
    !! 
    target_lai(:) = calc_target_lai(lctlib%growthform, lctlib%sla, lctlib%k_latosa, & 
                                    lctlib%wood_density, leaf2root_troot_mavg(:), &
                                    height(:), veg_pool%fine_root%carbon(:), veg_pool%sap_wood%carbon(:))

    CALL calc_target_labile_reserve_pool_size(nc                                    , & ! in
                                              lctlib%growthform                     , &
                                              lctlib%fstore_target                  , &
                                              lctlib%sla                            , &
                                              lctlib%tau_leaf                       , &
                                              lctlib%tau_fine_root                  , &
                                              gpp_tlabile_mavg(:)                   , &
                                              maint_respiration_tlabile_mavg(:)     , & 
                                              growth_req_n_tlabile_mavg(:)          , &
                                              growth_req_p_tlabile_mavg(:)          , &
                                              leaf2root_troot_mavg(:)               , &
                                              target_cn_leaf(:)                     , &
                                              target_cn_fine_root(:)                , &
                                              target_np_leaf(:)                     , &
                                              target_np_fine_root(:)                , & 
                                              target_lai(:)                         , &
                                              veg_pool%leaf%carbon(:)               , & 
                                              veg_pool%fine_root%carbon(:)          , & 
                                              veg_pool%sap_wood%carbon(:)           , & ! in
                                              target_labile_pool_carbon(:)          , & ! out  (local variable)
                                              target_labile_pool_nitrogen(:)        , & ! out  (local variable) 
                                              target_labile_pool_phosphorus(:)      , & ! out  (local variable) 
                                              target_reserve_pool_carbon(:)         , & ! out  (local variable) 
                                              target_reserve_pool_nitrogen(:)       , & ! out  (local variable) 
                                              target_reserve_pool_phosphorus(:)     )   ! out  (local variable)

    !> 4.2 calculate net flux between reserve and labile pool 
    !! 
    !! Info: veg_reserve_use is a pool structure
    CALL calc_reserve_useage(nc                                   , & ! in
                             growing_season(:)                    , &
                             kstar_labile(:)                      , &
                             target_labile_pool_carbon(:)         , & 
                             target_labile_pool_nitrogen(:)       , & 
                             target_labile_pool_phosphorus(:)     , & 
                             target_reserve_pool_carbon(:)        , & 
                             target_reserve_pool_nitrogen(:)      , & 
                             target_reserve_pool_phosphorus(:)    , &
                             target_lai(:)                        , &
                             lai(:)                               , &
                             veg_pool%labile                      , &
                             veg_pool%reserve                     , & ! in
                             veg_reserve_use                      )   ! inout 

    !> 4.3 signal sink limitation, identified by an accummulation of labile carbon, 
    !!     to photosynthesis routine                             
    !!
    beta_sinklim_ps(:) = calc_beta_sinklimitation(veg_pool%labile%carbon(:)  , &
                                                  veg_pool%labile%nitrogen(:), &
                                                  veg_pool%labile%phosphorus(:), &
                                                  target_labile_pool_carbon(:), &
                                                  target_labile_pool_nitrogen(:), &  
                                                  target_labile_pool_phosphorus(:))

    !> 5.0 calculate allocation coefficients for structural plant pools (leaves, wood, fine and coarse root, fruits)
    !! 
    !! 
    CALL calc_allocation_fraction(nc                                                , & ! in
                                  lctlib%growthform                                 , & 
                                  lctlib%k2_fruit_alloc                             , & 
                                  lctlib%k_crtos                                    , & 
                                  TRIM(dsl4jsb_Config(VEG_)%biomass_alloc_scheme)   , &
                                  dheight(:)                                        , &
                                  kstar_labile(:)                                   , &
                                  veg_pool%leaf%carbon(:)                           , &
                                  veg_pool%sap_wood%carbon(:)                       , &
                                  veg_pool%coarse_root%carbon(:)                    , &
                                  veg_pool%fine_root%carbon(:)                      , &
                                  veg_pool%labile%carbon(:)                         , &
                                  height(:)                                         , &
                                  leaf2sapwood_mass_ratio(:)                        , &
                                  leaf2root_mass_ratio(:)                           , & 
                                  k1_opt(:)                                         , &
                                  k2_opt(:)                                         , &
                                  target_cn_leaf(:)                                 , &
                                  target_cn_fine_root(:)                            , &
                                  target_cn_coarse_root(:)                          , &
                                  target_cn_sap_wood(:)                             , &
                                  target_cn_fruit(:)                                , &
                                  target_np_leaf(:)                                 , &
                                  target_np_fine_root(:)                            , &
                                  target_np_coarse_root(:)                          , &
                                  target_np_sap_wood(:)                             , &
                                  target_np_fruit(:)                                , &
                                  veg_reserve_use%carbon(:)                         , & ! in
                                  veg_falloc                                        , & ! inout
                                  growth_req_n(:)                                   , & ! out
                                  growth_req_p(:)                                   )   ! out

    !> 6.0 partition labile carbon into growth, growth_respiration, and maintenance 
    !! 
    !! 
    CALL calc_partitioning(nc                                     , & ! in
                           kstar_labile(:)                        , &
                           maint_respiration_pot(:)               , &
                           maint_respiration_leaf(:)              , &
                           growth_req_p(:)                        , &
                           growth_req_n(:)                        , &
                           veg_falloc                             , &
                           veg_pool%labile                        , & ! in
                           n_processing_respiration(:)            , & ! inout
                           maint_respiration_c13(:)               , & ! inout
                           maint_respiration_c14(:)               , & ! inout
                           growth_respiration_c13(:)              , & ! inout
                           growth_respiration_c14(:)              , & ! inout
                           n_processing_respiration_c13(:)        , & ! inout
                           n_processing_respiration_c14(:)        , & ! inout
                           veg_growth                             , & ! inout
                           veg_exudation                          , & ! inout
                           maint_respiration(:)                   , & ! out
                           growth_respiration(:)                  )   ! out

#ifdef __QUINCY_STANDALONE__
! Pass back values to process memory
#define Passback_veg_pool
#define Passback_veg_growth
#define Passback_veg_falloc
#define Passback_veg_exudation
#define Passback_veg_reserve_use
#define Passback_sb_pool
#define Passback_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Passback_veg_pool
#undef Passback_veg_growth
#undef Passback_veg_falloc
#undef Passback_veg_exudation
#undef Passback_veg_reserve_use
#undef Passback_sb_pool
#undef Passback_sb_mycorrhiza_export
#endif

  END SUBROUTINE update_veg_growth 


  !-----------------------------------------------------------------------------------------------------
  ! 1.0 Sub Task of update_veg_growth
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculates the activity of the meristem for growth 
  !!
  !! As this is proportional to the size of the labile pool, there is no need for a direct temperature 
  !! function of the turnover, as this will result from the dynamics of the labile pool. \nThis function
  !!  therefore only takes account of the impact of low temperatures (<5degC) and low moisture 
  !! (<20\% saturation). \nAdditionally, phenological dormancy is accounted for.
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_meristem_activity(t_air,w_soil_root_theta,growing_season) RESULT(kstar_labile)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_veg_constants,           ONLY: lambda_labile_temp  , &
                                          lambda_labile_theta , &
                                          k_labile_temp       , & 
                                          k_labile_theta      , &
                                          tau_labile
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),   INTENT(in)    :: t_air                  !< air temperature (K) 
    REAL(wp),   INTENT(in)    :: w_soil_root_theta      !< soil moisture scalar (fraction of available water holding capacity) 
    REAL(wp),   INTENT(in)    :: growing_season         !< is the plant growing? 
    REAL(wp)                  :: kstar_labile           !< current potential labile pool turnover (1/day)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                  :: ftemp                  ! temperature limit on cell construction rate
    REAL(wp)                  :: fmoist                 ! moisture limit on cell construction rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_meristem_activity'

  
    IF(t_air > Tzero) THEN 
        ftemp = 1.0_wp - exp(-(lambda_labile_temp * (t_air - Tzero))**k_labile_temp)
      ELSE 
        ftemp = 0.0_wp
    ENDIF
    fmoist = 1.0_wp - exp(-(lambda_labile_theta * w_soil_root_theta)**k_labile_theta) 
    IF(growing_season > test_false_true)THEN
        kstar_labile = ftemp * fmoist / tau_labile
      ELSE 
        kstar_labile = 0.0_wp
    ENDIF

  END FUNCTION calc_meristem_activity 


  !----------------------------------------------------------------------------------------------------------
  ! 1.1 Sub Task of update_veg_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Routine to calculate N and P demand to regulate N and P uptake 
  !!   
  !! output: N and P demand scalar for root uptake 
  !----------------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_nutrient_demand_scalar(lctlib_cn_leaf, &
                                                   lctlib_np_leaf, & 
                                                   leaf2root_mass_ratio, & 
                                                   labile_carbon, & 
                                                   labile_nitrogen, & 
                                                   labile_phosphorus, & 
                                                   f_n_demand, &
                                                   f_p_demand)

    USE mo_kind,                        ONLY: wp
    USE mo_jsb_math_constants,          ONLY: eps8
    USE mo_veg_constants,               ONLY: kappa_f_demand_nc, kappa_f_demand_pn, k_f_demand, &
                                              root2leaf_cn, root2leaf_np, fresp_growth

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),               INTENT(in)  :: lctlib_cn_leaf            !< PFT-specific foliar C:N (mol/mol) 
    REAL(wp),               INTENT(in)  :: lctlib_np_leaf            !< PFT-specific foliar N:P (mol/mol)
    REAL(wp),               INTENT(in)  :: leaf2root_mass_ratio      !< growing season mean leaf 2 fine root ratio 
    REAL(wp),               INTENT(in)  :: labile_carbon             !< labile pool carbon (mol C / m2)
    REAL(wp),               INTENT(in)  :: labile_nitrogen           !< labile pool nitrogen (mol N / m2)
    REAL(wp),               INTENT(in)  :: labile_phosphorus         !< labile pool phosphorus (mol P / m2)
    REAL(wp),               INTENT(out) :: f_n_demand                !< root N uptake demand scalar (unitless)
    REAL(wp),               INTENT(out) :: f_p_demand                !< root P uptake demand scalar (unitless)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: nc_max, pn_max, hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_nutrient_demand_scalar'


    !> 1.0 calculate reference level of N and P in labile pool that would induce zero uptake 
    !! The highest useful level corresponds to the C:N:P expected for leaves and fine root, 
    !! given the average root-lifetime ratio of leaf to root mass
    !! and the lctlib levels of leaf and fine root C:N:P \n
    !! i.e. CN = CN_leaf,pft * (1 + 1 / leaf2root_mass) / (1 + cn_root2leaf / leaf2root_mass) \n
    !! Since growth respiration takes some of the C away from growth, this is accounted in the C:N:P ratio
    nc_max = 1._wp / (lctlib_cn_leaf * & 
                      ((1._wp + 1._wp / leaf2root_mass_ratio) / &
                       (1._wp + root2leaf_cn / leaf2root_mass_ratio ) * (1._wp + fresp_growth)))
    pn_max = 1._wp / (lctlib_np_leaf / & 
                      ((1._wp + root2leaf_cn / leaf2root_mass_ratio) / &
                       (1._wp + root2leaf_cn * root2leaf_np / leaf2root_mass_ratio ))) 

    !> 1.1 N demand 
    !!
    !! lambda_f_demand_nc acconts for the range of C:N at which uptake becomes reduces
    IF(labile_carbon > eps8) THEN
      hlp1       = MAX(nc_max - labile_nitrogen / labile_carbon,0.0_wp) / MAX(nc_max * (1._wp - kappa_f_demand_nc),eps8)
      f_n_demand = 1._wp - exp(-hlp1 ** k_f_demand)
    ELSE
      f_n_demand = 1.0_wp
    ENDIF

    !> 1.2 P demand
    !!
    IF(labile_nitrogen > eps8) THEN
      hlp1       = MAX(pn_max - labile_phosphorus / labile_nitrogen,0.0_wp) / MAX(pn_max * (1._wp - kappa_f_demand_pn),eps8)
      f_p_demand = 1._wp - exp(-hlp1 ** k_f_demand)
    ELSE
      f_p_demand = 1.0_wp
    ENDIF

  END SUBROUTINE calc_nutrient_demand_scalar


  !----------------------------------------------------------------------------------------------------------
  ! 1.1 Sub Task of update_veg_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Routine to calculate actual uptake rates of nitrogen and phosphorus from soil given the 
  !! vegetation demands and resource availability 
  !!   
  !! output: actual uptake rates and associated respiration costs 
  !----------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_actual_nutrient_uptake_rate(nc, & 
                                              nsoil_sb, &
                                              lctlib_cn_leaf, &
                                              lctlib_np_leaf, &
                                              lctlib_tau_fine_root, &
                                              quincy_model, &
                                              bnf_scheme, &
                                              sb_model_scheme, &
                                              flag_mycorrhiza  , &
                                              kstar_labile, &
                                              t_air, &                                              
                                              fine_root_carbon, & 
                                              fine_root_nitrogen, & 
                                              labile_carbon, & 
                                              labile_nitrogen, & 
                                              labile_phosphorus, & 
                                              sb_pool_mycorrhiza_carbon, & 
                                              nh4_solute, & 
                                              no3_solute, & 
                                              soil_depth_sl, &
                                              growth_req_n_tlabile_mavg, &
                                              demand_uptake_n_tuptake_mavg, &
                                              demand_uptake_p_tuptake_mavg, &
                                              fmaint_rate_troot_mavg, &
                                              unit_npp_troot_mavg, &
                                              unit_uptake_n_troot_mavg, &
                                              cost_n_fixation_talloc_mavg, &
                                              exudation_c_tmyc_mavg, &
                                              unit_uptake_n_pot, &
                                              unit_uptake_p_pot, &
                                              plant_uptake_nh4_sl, &
                                              plant_uptake_nh4_n15_sl, &
                                              plant_uptake_no3_sl, &
                                              plant_uptake_no3_n15_sl, &
                                              plant_uptake_po4_sl, &
                                              sb_mycorrhiza_export, &
                                              myc_export_n_tlabile_mavg_sl, &
                                              myc_export_p_tlabile_mavg_sl, &
                                              myc_export_c_tmyc_mavg_sl, &
                                              myc_export_n_tmyc_mavg_sl, &
                                              uptake_nh4, &
                                              uptake_nh4_n15, &
                                              uptake_no3, &
                                              uptake_no3_n15, &
                                              uptake_po4, &
                                              n_fixation, & 
                                              unit_uptake_n_act, &
                                              unit_uptake_p_act, &
                                              cost_n_uptake_root, &
                                              cost_n_fixation, &
                                              n_transform_respiration, &
                                              n_fixation_respiration, &
                                              n_processing_respiration, &
                                              veg_exudation_carbon)
                       
    USE mo_kind,                    ONLY: wp
    USE mo_jsb_physical_constants,  ONLY: Tzero, molar_mass_C, molar_mass_N
    USE mo_isotope_util,            ONLY: calc_mixing_ratio_N15N14
    USE mo_jsb_math_constants,      ONLY: dtime, one_day, one_year, eps8
    USE mo_veg_constants,           ONLY: vmax_bnf, k_bnf, eta_nfixation, transform_cost_nh4, transform_cost_no3, &
                                          fresp_growth, cost_n_fixation_max, cost_n_fixation_min, & 
                                          k1_temp_bnf, k2_temp_bnf,  k3_temp_bnf,  k4_temp_bnf, tau_labile, &
                                          f_myc_exudation_max, f_mycorrhization_min, f_mycorrhization_max
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_sb
#else
    dsl4jsb_Use_pools(SB_) t_sb_elements_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                              INTENT(in)    :: nc, &                        !< dimensions
                                                           nsoil_sb
    REAL(wp),                             INTENT(in)    :: lctlib_cn_leaf               !< lctlib parameter
    REAL(wp),                             INTENT(in)    :: lctlib_np_leaf               !< lctlib parameter
    REAL(wp),                             INTENT(in)    :: lctlib_tau_fine_root         !< lctlib parameter
    CHARACTER(len=*),                     INTENT(in)    :: quincy_model                 !< model name: land plant soil canopy
    CHARACTER(len=*),                     INTENT(in)    :: bnf_scheme                   !< from veg config 
    CHARACTER(len=*),                     INTENT(in)    :: sb_model_scheme              !< from sb config
    LOGICAL,                              INTENT(in)    :: flag_mycorrhiza              !< from sb config
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: kstar_labile                 !< meristem turnover rate  
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: t_air                        !< air temperature (K) !! should be soil temperature!
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: fine_root_carbon             !< fine root carbon (mol C / m2)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: fine_root_nitrogen           !< fine root nitrogen (mol N / m2)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: labile_carbon                !< labile pool carbon (mol C / m2)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: labile_nitrogen              !< labile pool nitrogen (mol N / m2)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: labile_phosphorus            !< labile pool phosphorus (mol P / m2)
    REAL(wp), DIMENSION(nc, nsoil_sb),    INTENT(in)    :: sb_pool_mycorrhiza_carbon, & !< mycorrhiza biomass (mol C /m3)
                                                           nh4_solute         , &       !< soil solution ammonium (mol N /m3)
                                                           no3_solute         , &       !< soil solution nitrate (mol N /m3)
                                                           soil_depth_sl                !< depth of soil layers (m)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: growth_req_n_tlabile_mavg    !< time-averaged growth requirement for N 
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: demand_uptake_n_tuptake_mavg !< time-averaged N demand scalar (unitless) 
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: demand_uptake_p_tuptake_mavg !< time-averaged P demand scalar (unitless)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: fmaint_rate_troot_mavg       !< time-averaged N-speficic root maintentance respiration rate (micro-mol C/mol N / s)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: unit_npp_troot_mavg          !< time-averaged increase of NPP with increased leaf C investment (mol C/mol C)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: unit_uptake_n_troot_mavg     !< time-averaged increase of root N uptake with increased root C inventment (mol N / mol C)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: cost_n_fixation_talloc_mavg  !< time-averaged C costs for N fixation (mol C / mol N)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: exudation_c_tmyc_mavg        !< time-averaged C allocation to mycorrhizae (mol C / mol N)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: unit_uptake_n_pot            !< potential N uptake per unit root mass (mu mol N / mol C / s)
    REAL(wp), DIMENSION(nc),              INTENT(in)    :: unit_uptake_p_pot            !< potential P uptake per unit root mass (mu mol P / mol C / s)
    REAL(wp), DIMENSION(nc, nsoil_sb),    INTENT(inout) :: plant_uptake_nh4_sl     , &  !< plant NH4 uptake per soil layer (mu mol N / m3 / s)
                                                           plant_uptake_nh4_n15_sl , &  !< plant NH4-N15 uptake per soil layer (mu mol N / m3 / s)
                                                           plant_uptake_no3_sl     , &  !< plant NO3 uptake per soil layer (mu mol N / m3 / s)
                                                           plant_uptake_no3_n15_sl , &  !< plant NO3-N15 uptake per soil layer (mu mol N / m3 / s)
                                                           plant_uptake_po4_sl, &       !< plant PO4 uptake per soil layer (mu mol P / m3 / s)
                                                           myc_export_n_tlabile_mavg_sl  , & !< time averaged N export from myc2plant (mu mol N / m3)
                                                           myc_export_p_tlabile_mavg_sl  , & !< time averaged P export from myc2plant (mu mol P / m3)
                                                           myc_export_c_tmyc_mavg_sl  , & !< time averaged C export from myc2plant (mu mol N / m3)
                                                           myc_export_n_tmyc_mavg_sl  !< time averaged N export from myc2plant (mu mol P / m3)
#ifdef __QUINCY_STANDALONE__
    TYPE(t_sb_element_a),                 INTENT(inout) :: sb_mycorrhiza_export         !< total plant uptake from mycorrhiza (mol / m2 / timestep)  
#else
    TYPE(t_sb_elements_pool),             INTENT(inout) :: sb_mycorrhiza_export         !< total plant uptake from mycorrhiza (mol / m2 / timestep)  
#endif
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: uptake_nh4                   !< total plant NH4 uptake (mu mol N / m2 / s)  
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: uptake_nh4_n15               !< total plant NH4-N15 uptake (mu mol N / m2 / s)  
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: uptake_no3                   !< total plant NO3 uptake (mu mol N / m2 / s)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: uptake_no3_n15               !< total plant NO3-N15 uptake (mu mol N / m2 / s)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: uptake_po4                   !< total plant PO4 uptake (mu mol N / m2 / s)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: n_fixation                   !< symbiontic N fixation (mu mol N / m2 / s)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: unit_uptake_n_act            !< N uptake per unit root mass (mu mol N / mol C / s)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: unit_uptake_p_act            !< P uptake per unit root mass (mu mol P / mol C / s)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: cost_n_uptake_root           !< C costs for N uptake by roots (mol C / mol N)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: cost_n_fixation              !< C costs for N fixation (mol C / mol N)
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: n_transform_respiration      !< respiration associated with N transformation (mu mol C / m2 / s)     
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: n_fixation_respiration       !< respiration associated with N fixation (mu mol C / m2 / s)     
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: n_processing_respiration     !< respiration associated with N transformation and fixation (mu mol C / m2 / s)  
    REAL(wp), DIMENSION(nc),              INTENT(out)   :: veg_exudation_carbon         !< C investment into mycorrhizal fungi (mol C / m2 / timestep) 
    ! --------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)         :: dcost, &          ! cost_n_uptake_root - cost_n_fixation (mol C / mol N)
                                       n_demand, &       ! N demand of plants (mu mol N /m2 / s) 
                                       f_nacq_bnf, &     ! fraction of N uptake as BFN
                                       f_nacq_nup, &     ! fraction of N uptake as Nuptake
                                       uptake_nh4_pot, & ! potential NH4 uptake (mu mol N / m2 /s) 
                                       uptake_no3_pot, & ! potential NO3 uptake (mu mol N / m2 /s)
                                       uptake_po4_pot, & ! potential PO4 uptake (mu mol P / m2 /s)
                                       fresp_maint, &    ! fractional cost of maintenance respiration
                                       vmax_bnf_act, &   ! actual T-dependent vmax_bnf
                                       gain_c, &         ! marginal carbon gain per leaf investment
                                       gain_n, &         ! marginal nitrogen gain per root investment
                                       gain_mycorrhizae, &
                                       support, &
                                       nsup, &
                                       psup, &
                                       f_exudate_max, &
                                       f_exudate, &
                                       f_myc_opt, &
                                       fexu, &
                                       scal, &  
                                       hlp1, &
                                       hlp2 
    INTEGER                         :: isoil      ! looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_actual_nutrient_uptake_rate'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !>0.9 init output variables
    !>
    uptake_nh4(:)               = 0.0_wp
    uptake_nh4_n15(:)           = 0.0_wp
    uptake_no3(:)               = 0.0_wp
    uptake_no3_n15(:)           = 0.0_wp
    uptake_po4(:)               = 0.0_wp
    n_fixation(:)               = 0.0_wp
    unit_uptake_n_act(:)        = 0.0_wp
    unit_uptake_p_act(:)        = 0.0_wp
    cost_n_uptake_root(:)       = 0.0_wp
    cost_n_fixation(:)          = 0.0_wp
    n_transform_respiration(:)  = 0.0_wp
    n_fixation_respiration(:)   = 0.0_wp
    n_processing_respiration(:) = 0.0_wp
    veg_exudation_carbon(:)     = 0.0_wp

    !> 1.0 determine root nutrient uptake from soil  
    !!
    !> 1.1 in case of the plant-only model, prescribe amount of N and P taken up
    !!
    IF(quincy_model == 'plant')THEN
       ! in case of the plant-only model, prescribe N and P uptake to match demand
       plant_uptake_nh4_sl(:,1)     = MAX(labile_carbon(:) / lctlib_cn_leaf - labile_nitrogen(:), &
                                          0.0_wp) / soil_depth_sl(:,1)
       plant_uptake_po4_sl(:,1)     = MAX(labile_carbon(:) / lctlib_cn_leaf / lctlib_np_leaf - labile_phosphorus(:), &
                                          0.0_wp)  / soil_depth_sl(:,1)
       plant_uptake_nh4_n15_sl(:,1) = plant_uptake_nh4_sl(:,1) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
       plant_uptake_no3_n15_sl(:,1) = plant_uptake_no3_sl(:,1) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    ENDIF

    !> 1.2 integrate root uptake over soil layers
    !!
    uptake_nh4_pot(:) = 0.0_wp
    uptake_no3_pot(:) = 0.0_wp 
    uptake_po4_pot(:) = 0.0_wp
    DO isoil = 1,nsoil_sb
       uptake_nh4_pot(:) = uptake_nh4_pot(:) + plant_uptake_nh4_sl(:,isoil) * soil_depth_sl(:,isoil)
       uptake_no3_pot(:) = uptake_no3_pot(:) + plant_uptake_no3_sl(:,isoil) * soil_depth_sl(:,isoil)
       uptake_po4_pot(:) = uptake_po4_pot(:) + plant_uptake_po4_sl(:,isoil) * soil_depth_sl(:,isoil)
    ENDDO
   
    !> 1.3 apply demand driven reduction in uptake
    !! 
    uptake_nh4(:) = demand_uptake_n_tuptake_mavg(:) * uptake_nh4_pot(:)
    uptake_no3(:) = demand_uptake_n_tuptake_mavg(:) * uptake_no3_pot(:)
    uptake_po4(:) = demand_uptake_p_tuptake_mavg(:) * uptake_po4_pot(:)


    !> 2.0 C costs for N acquisition
    !!
    !> 2.1 C costs for (inorganic) N uptake, based on Rastetter (2001) and Meyerholt (2016, OPT)
    !!
    ! photosynthesis gain relative to carbon investment into leafs
    gain_c(:) = MAX(unit_npp_troot_mavg(:),0.0_wp)
 
    ! n-uptake gain relative to carbon investment into roots
    ! multiplier is Croot / Croot * (1 + fresp_growth + fresp_maint
    ! maintenance respiration cost of a unit fine root carbon over its lifetime (mol C/mol C)
    WHERE(fine_root_carbon(:) > eps8)
      fresp_maint(:) = fmaint_rate_troot_mavg(:) * fine_root_nitrogen(:) / fine_root_carbon(:)  * &
                       lctlib_tau_fine_root * one_day * one_year / 1.e6_wp
    ELSEWHERE
      fresp_maint(:) = 0.0_wp
    END WHERE
    gain_n(:) = unit_uptake_n_troot_mavg(:) / (1._wp + fresp_growth + fresp_maint(:))

    ! costs for mycorrhizae
    gain_mycorrhizae(:) = 0.0_wp
    IF(flag_mycorrhiza) THEN
      hlp1(:) = exudation_c_tmyc_mavg(:) - SUM(myc_export_c_tmyc_mavg_sl(:,:) * soil_depth_sl(:,:),DIM=2)
      hlp2(:) = SUM(myc_export_n_tmyc_mavg_sl(:,:) * soil_depth_sl(:,:),DIM=2)
      WHERE(hlp1(:) > eps8)
          gain_mycorrhizae(:) = hlp2(:) / hlp1(:)
      ENDWHERE
    ENDIF

    gain_n(:) = gain_n(:) + gain_mycorrhizae(:)
 
    ! transformation costs depend on NH4 and NO3 uptake rates   
    WHERE((uptake_nh4(:) + uptake_no3(:)) > eps8)
       hlp1(:) = uptake_nh4(:) / (uptake_nh4(:) + uptake_no3(:))
       hlp2(:) = uptake_no3(:) / (uptake_nh4(:) + uptake_no3(:))
    ELSEWHERE
       hlp1(:) = 1.0_wp
       hlp2(:) = 0.0_wp
    END WHERE

    ! root uptake costs
    WHERE(gain_n(:) > eps8)
      cost_n_uptake_root(:) = gain_c(:) / gain_n(:) + &
                              (hlp1(:) * transform_cost_nh4 + hlp2(:) * transform_cost_no3)
    ELSEWHERE
      cost_n_uptake_root(:) = 1000._wp ! an arbitrarily large value !
    END WHERE
                       
    !> 2.2 C costs for BNF
    !! Temperature function based on Houlton (2008)
    hlp1(:)            = k1_temp_bnf * &
                         exp(k2_temp_bnf + k3_temp_bnf * (t_air(:) - Tzero) * (1.0_wp - (t_air(:) - Tzero) / k4_temp_bnf))
    cost_n_fixation(:) = (cost_n_fixation_max - (cost_n_fixation_max - cost_n_fixation_min) * hlp1(:)) 


    !> 3.0 N fixation
    !!
    SELECT CASE (bnf_scheme)
      
      !> 3.1 Case 1: optimal N fixation
      !! BNF calculation with Michaelis Menten function (Rastetter, 2001) (MM)
      CASE ("optimal")
        hlp1(:)            = 1.0_wp /(k1_temp_bnf * &
                             exp(k2_temp_bnf + k3_temp_bnf * (t_air(:) - Tzero) * (1.0_wp - (t_air(:) - Tzero) / k4_temp_bnf)))
        vmax_bnf_act(:)    = vmax_bnf * (0.25_wp + 0.25_wp * hlp1(:))
        WHERE(cost_n_fixation_talloc_mavg(:) <= cost_n_uptake_root(:))
           dcost(:)      = cost_n_uptake_root(:) - cost_n_fixation_talloc_mavg(:)
           n_fixation(:) = fine_root_carbon(:) * vmax_bnf_act(:) * dcost(:) / (k_bnf + dcost(:))
        ELSEWHERE
           n_fixation(:) = 0.0_wp
        END WHERE

      !> 3.2 Case 2: Give enough N to maintain N requirement in the labile pool 
      CASE("unlimited")
        ! N fixation in micro-mol N / m2 / s
        n_fixation(:) = MAX(labile_carbon(:) * growth_req_n_tlabile_mavg(:) - labile_nitrogen(:), 0.0_wp) / & 
                        dtime * 1000000.0_wp
    ENDSELECT 


    !> 4.0 Determine C use for Root Exudates
    !! 
    IF(flag_mycorrhiza) THEN

       !> relative support of N/P acquisition by mycorrhizae
       ! scaling factor for plant C exudation
       ! high support = high exudation // low/no support = low exudation
       !!only sensitive to N at the moment (last line)
       nsup(:) = 0.0_wp
       hlp1(:) = SUM(myc_export_n_tlabile_mavg_sl(:,:) * soil_depth_sl(:,:),DIM=2) * 1000000.0_wp !coversion of mols to micro-mols
       nsup(:) = hlp1(:) / (hlp1(:) + 0.0001)
       psup(:) = 0.0_wp
       hlp1(:) = SUM(myc_export_p_tlabile_mavg_sl(:,:) * soil_depth_sl(:,:),DIM=2) * 1000000.0_wp
       psup(:) = hlp1(:) / (hlp1(:) + 0.000001)
       support(:) = nsup(:) !MIN(nsup(:), psup(:))

       !> assuming that plants will give max XX% of freshly assimilated carbon to mycorrhiza
       ! scaled with demand and support
       !!only sensitive to N at the moment (last line)
       f_myc_opt(:) = f_mycorrhization_min + (f_mycorrhization_max - f_mycorrhization_min) * support * &
                                              demand_uptake_n_tuptake_mavg(:)
                                             !MAX(demand_uptake_n_tuptake_mavg(:),demand_uptake_p_tuptake_mavg(:)) 

       !> exudation rate (mol / m2 / timestep)
       hlp1(:) = SUM(sb_pool_mycorrhiza_carbon(:,:) * soil_depth_sl(:,:),DIM=2)
       hlp2(:) = MIN(f_myc_exudation_max * labile_carbon(:) * tau_labile / one_day, &
                     MAX(f_myc_opt(:) * fine_root_carbon(:) - hlp1(:),0.0_wp))
       veg_exudation_carbon(:) = MAX(f_mycorrhization_min * fine_root_carbon(:) - hlp1(:), hlp2(:))
                      

       !> scale down uptake from mycorrhiza according to demand
       ! scale mycorrhizal export with demand
       ! to avoid DOWN-regulation with demand (between 0 and 1), when demand is high, demand for scaling is intensified by 50%
       ! limited to 100% potential export
       sb_mycorrhiza_export%carbon(:,:)     = sb_mycorrhiza_export%carbon(:,:)     * &
                                              MIN( SPREAD(demand_uptake_n_tuptake_mavg(:), 2, nsoil_sb) * 1.5_wp, 1.0_wp)
       sb_mycorrhiza_export%nitrogen(:,:)   = sb_mycorrhiza_export%nitrogen(:,:)   * &
                                              MIN( SPREAD(demand_uptake_n_tuptake_mavg(:), 2, nsoil_sb) * 1.5_wp, 1.0_wp)
       sb_mycorrhiza_export%phosphorus(:,:) = sb_mycorrhiza_export%phosphorus(:,:) * &
                                              MIN( SPREAD(demand_uptake_p_tuptake_mavg(:), 2, nsoil_sb) * 1.5_wp, 1.0_wp)
       sb_mycorrhiza_export%carbon13(:,:)   = sb_mycorrhiza_export%carbon13(:,:)   * &
                                              MIN( SPREAD(demand_uptake_n_tuptake_mavg(:), 2, nsoil_sb) * 1.5_wp, 1.0_wp)
       sb_mycorrhiza_export%carbon14(:,:)   = sb_mycorrhiza_export%carbon14(:,:)   * &
                                              MIN( SPREAD(demand_uptake_n_tuptake_mavg(:), 2, nsoil_sb) * 1.5_wp, 1.0_wp)
       sb_mycorrhiza_export%nitrogen15(:,:) = sb_mycorrhiza_export%nitrogen15(:,:) * &
                                              MIN( SPREAD(demand_uptake_n_tuptake_mavg(:), 2, nsoil_sb) * 1.5_wp, 1.0_wp)
       veg_exudation_carbon(:)              = veg_exudation_carbon(:) !+ SUM(sb_mycorrhiza_export_carbon(:,:), DIM=2)
    ENDIF ! flag_mycorrhiza


    !> 5.0 Determine C use for N uptake, and constrain this to available carbon
    !!
    IF(bnf_scheme == 'optimal' .OR. bnf_scheme == 'standard') THEN
      hlp1(:) = uptake_nh4(:) * transform_cost_nh4 + & 
                uptake_no3(:) * transform_cost_no3 + & 
                n_fixation(:) * cost_n_fixation(:)
      ! adjust N uptake to possible C expenditure
      hlp2(:) = MIN(1.0_wp / tau_labile / one_day * labile_carbon(:) * 1.e6_wp, hlp1(:))
      WHERE(hlp1(:) > eps8 .AND. hlp1(:) > hlp2(:))
        scal(:)       = hlp2(:) / hlp1(:) 
        uptake_nh4(:) = uptake_nh4(:) * scal(:)
        uptake_no3(:) = uptake_no3(:) * scal(:)
        n_fixation(:) = n_fixation(:) * scal(:) 
      ENDWHERE
      n_transform_respiration(:)  = uptake_nh4(:) * transform_cost_nh4 + uptake_no3(:) * transform_cost_no3
      n_fixation_respiration(:)   = n_fixation(:) * cost_n_fixation(:)
      n_processing_respiration(:) = n_fixation_respiration(:) + n_transform_respiration(:)
    ENDIF
  

    !> 6.0 update diagnostics and soil related variables, as well as N15 fluxes
    !!
    !> 6.1 Calculate actual N and P uptake per unit fine root carbon (micro-mol {N/P} / mol C / s)
    !!
    WHERE(fine_root_carbon(:) > eps8)
      unit_uptake_n_act(:) = (uptake_nh4(:) + uptake_no3(:)) / fine_root_carbon(:)
      unit_uptake_p_act(:) = uptake_po4(:) / fine_root_carbon(:)
    ELSEWHERE
      unit_uptake_n_act(:) = unit_uptake_n_pot(:) * demand_uptake_n_tuptake_mavg(:)
      unit_uptake_p_act(:) = unit_uptake_p_pot(:) * demand_uptake_p_tuptake_mavg(:)
    END WHERE

    !> 6.2 apply actual uptake rates to soil layers
    !!
    DO isoil = 1,nsoil_sb
      WHERE(uptake_nh4_pot(:) > eps8 .AND. uptake_nh4_pot(:) > uptake_nh4(:))
        plant_uptake_nh4_sl(:,isoil)      = plant_uptake_nh4_sl(:,isoil) * uptake_nh4(:) / uptake_nh4_pot(:) 
        plant_uptake_nh4_n15_sl(:,isoil)  = plant_uptake_nh4_n15_sl(:,isoil) * uptake_nh4(:) / uptake_nh4_pot(:) 
      ENDWHERE
      WHERE(uptake_no3_pot(:) > eps8 .AND. uptake_no3_pot(:) > uptake_no3(:))
        plant_uptake_no3_sl(:,isoil)      = plant_uptake_no3_sl(:,isoil) * uptake_no3(:) / uptake_no3_pot(:) 
        plant_uptake_no3_n15_sl(:,isoil)  = plant_uptake_no3_n15_sl(:,isoil) * uptake_no3(:) / uptake_no3_pot(:) 
      ENDWHERE
      WHERE(uptake_po4_pot(:) > eps8 .AND. uptake_po4_pot(:) > uptake_po4(:))
        plant_uptake_po4_sl(:,isoil)      = plant_uptake_po4_sl(:,isoil) * uptake_po4(:) / uptake_po4_pot(:) 
      ENDWHERE
    END DO

    !> 6.3 calculate actual total plant N15 uptake
    !!
    uptake_nh4_n15(:) = 0.0_wp
    uptake_no3_n15(:) = 0.0_wp 
    DO isoil = 1,nsoil_sb
      uptake_nh4_n15(:) = uptake_nh4_n15(:) + plant_uptake_nh4_n15_sl(:,isoil) * soil_depth_sl(:,isoil)
      uptake_no3_n15(:) = uptake_no3_n15(:) + plant_uptake_no3_n15_sl(:,isoil) * soil_depth_sl(:,isoil)
    ENDDO
 
  END SUBROUTINE calc_actual_nutrient_uptake_rate


  !-----------------------------------------------------------------------------------------------------
  ! 2.0 Sub Task of update_veg_growth
  ! 
  ! ------------------------------------------------------------------------------------------
  !> Routine to calculate the target stoichiometry of the newly growth tissue
  !! 
  !! Includes optimal CN ratio calculation \n
  !! output: target stoichiometries
  ! --------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_stoichiometry_changes(growthform              , &   ! lctlib
                                                  phenology_type          , &
                                                  cn_leaf                 , &   ! lctlib
                                                  np_leaf                 , &
                                                  cn_leaf_min             , &
                                                  cn_leaf_max             , &
                                                  np_leaf_min             , &
                                                  np_leaf_max             , &
                                                  leaf_stoichom_scheme    , &   ! in
                                                  growing_season          , &
                                                  leaf_cn_direction       , &
                                                  labile_carbon           , &
                                                  labile_nitrogen         , &
                                                  labile_phosphorus       , &
                                                  cn_growth               , &
                                                  np_growth               , &
                                                  lai                     , &
                                                  target_lai              , &
                                                  target_cn_leaf          , &   ! inout
                                                  target_np_leaf          , &
                                                  target_cn_fine_root     , &   ! out
                                                  target_cn_coarse_root   , &
                                                  target_cn_sap_wood      , &
                                                  target_cn_fruit         , &
                                                  target_np_fine_root     , &
                                                  target_np_coarse_root   , &
                                                  target_np_sap_wood      , &
                                                  target_np_fruit           )

    USE mo_jsb_math_constants,      ONLY: dtime, one_day, eps8
    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_veg_constants          ! e.g. delta_n_leaf
    USE mo_pheno_constants,         ONLY: ievergreen, iraingreen, isummergreen, iperennial

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                  INTENT(in)    :: growthform                !< land-cover-type library parameter
    INTEGER,                  INTENT(in)    :: phenology_type            !< land-cover-type library parameter
    REAL(wp),                 INTENT(in)    :: cn_leaf               , & !< land-cover-type library parameter
                                               cn_leaf_min           , & !< 
                                               cn_leaf_max           , &
                                               np_leaf               , &
                                               np_leaf_min           , & !<
                                               np_leaf_max           
    CHARACTER(len=*),         INTENT(in)    :: leaf_stoichom_scheme      !< from veg config
    REAL(wp),                 INTENT(in)    :: growing_season            !< 
    REAL(wp),                 INTENT(in)    :: leaf_cn_direction     
    REAL(Wp),                 INTENT(in)    :: labile_carbon         , & !< average labile carbon
                                               labile_nitrogen       , & !< average labile nitrogen
                                               labile_phosphorus     , & !< average labile phpsphorus
                                               cn_growth             , & !< average CN ratio of new growth
                                               np_growth
    REAL(wp),                 INTENT(in)    :: lai                   , &
                                               target_lai
    REAL(wp),                 INTENT(inout) :: target_cn_leaf        , &    !<
                                               target_np_leaf
    REAL(wp),                 INTENT(out)   :: target_cn_fine_root   , & !< fine root target C:N
                                               target_cn_coarse_root , & !< coarse root target C:N
                                               target_cn_sap_wood    , & !< sap wood target C:N
                                               target_cn_fruit       , & !< fruit target C:N
                                               target_np_fine_root   , & !< fine root target N:P
                                               target_np_coarse_root , & !< coarse root target N:P
                                               target_np_sap_wood    , & !< sap wood target N:P
                                               target_np_fruit           !< fruit target N:P
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: delta_cn, delta_np, f_growing_leaf
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_stoichiometry_changes'

    !> 1.0 calculate leaf C:N
    !!

    ! factor to limit changes in stoichiometry at the start of the growing season
    SELECT CASE (phenology_type)
    CASE (iraingreen, isummergreen)
       IF (lai/MAX(eps8, target_lai) < 0.5_wp) THEN
          f_growing_leaf = 0._wp
       ELSE
          f_growing_leaf = 1._wp
       ENDIF
    CASE DEFAULT
       f_growing_leaf = 1._wp
    END SELECT
    
    IF (growing_season > test_false_true) THEN
       SELECT CASE (leaf_stoichom_scheme)
       CASE ("optimal")
          ! update leaf CN target. leaf_cn_direction has values of -1, 0 or 1
          target_cn_leaf = target_cn_leaf * (1.0_wp - leaf_cn_direction * delta_n_leaf * dtime/one_day)

       CASE ("dynamic")
          
             IF (labile_carbon / MAX(eps8, labile_nitrogen) >= cn_growth) THEN
             ! If plants are N limited, i.e. if the CN of the labile pool is higher than the growth CN
             ! increase leaf CN (decrease leaf N concentration)
             delta_cn = (exp(-(2._wp * target_cn_leaf/(cn_leaf_min + cn_leaf_max))**8._wp))
          ELSE
             ! If plants are not N limited and have more N than they need
             ! decrease leaf CN (increase leaf N concentration)
             delta_cn = - (1._wp - exp(-(2._wp*target_cn_leaf/(cn_leaf_min + cn_leaf_max))**8._wp))
          END IF
          ! Apply CN modifier calculated above to actual change in leaf CN target, given by the maximum change per timestep
          target_cn_leaf = target_cn_leaf * (1._wp +  delta_cn  * f_growing_leaf * delta_n_leaf * dtime/one_day)

 
       END SELECT
    ENDIF

    !> 2.0 determine C:N of other plant tissues
    !!
    target_cn_fine_root = target_cn_leaf / root2leaf_cn
    SELECT CASE(growthform)
      CASE (itree) ! trees
         target_cn_coarse_root  = cn_leaf / wood2leaf_cn
         target_cn_sap_wood     = cn_leaf / wood2leaf_cn 
      CASE (igrass) ! grasses 
         target_cn_coarse_root  = cn_leaf / root2leaf_cn
         target_cn_sap_wood     = cn_leaf / root2leaf_cn 
    END SELECT
    target_cn_fruit = cn_leaf / root2leaf_cn

    !> 3.0 determine N:P of other plant tissues
    !!

    ! Leaf N:P
    IF (growing_season > test_false_true) THEN
       SELECT CASE (leaf_stoichom_scheme)
       CASE ("dynamic")
          
          IF (labile_nitrogen / MAX(eps8, labile_phosphorus) >= np_growth) THEN
             ! If plants are P limited, i.e. if the NP of the labile pool is higher than the growth NP
             ! increase leaf NP (decrease leaf P concentration)
             delta_np = (exp(-(2._wp * target_np_leaf/(np_leaf_min + np_leaf_max))**8._wp))
          ELSE
             ! If plants are not P limited and have more P than they need
             ! decrease leaf NP (increase leaf P concentration)
             delta_np = - (1._wp - exp(-(2._wp*target_np_leaf/(np_leaf_min + np_leaf_max))**8._wp))
          END IF
          ! Apply NP modifier calculated above to actual change in leaf NP target, given by the maximum change per timestep
             target_np_leaf = target_np_leaf * (1._wp +  delta_np  * f_growing_leaf * delta_n_leaf * dtime/one_day)
       END SELECT
    ENDIF

    target_np_fine_root = target_np_leaf / root2leaf_np
    SELECT CASE(growthform)
      CASE (itree) ! trees
         target_np_coarse_root  = np_leaf / wood2leaf_np
         target_np_sap_wood     = np_leaf / wood2leaf_np
      CASE (igrass) ! grasses 
         target_np_coarse_root  = np_leaf / root2leaf_np
         target_np_sap_wood     = np_leaf / root2leaf_np
    END SELECT
    target_np_fruit = np_leaf / root2leaf_np

  END SUBROUTINE calc_stoichiometry_changes


  !-----------------------------------------------------------------------------------------------------
  ! 3.2 Sub Task of update_veg_growth  (3.1 calc_marginal_canopy_flux_increment:mo_veg_canopy)
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculates the current allometry of a tree 
  !!
  !! diagnoses diameter and height from woody biomass by assuming that the trunk is a cylinder following
  !!  a given allometric relationship: H = k_allom1 * D ** k_allom2
  !!
  !! calculates leaf to sapwood and root mass ratios
  !!
  !! calculates maximum LAI given current root and sapwood mass
  !-----------------------------------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE calc_allometry( growthform                    , &  ! lctlib
                                            allom_k1                      , &
                                            allom_k2                      , &
                                            sla                           , &
                                            wood_density                  , &
                                            tau_leaf                      , &
                                            tau_fine_root                 , &
                                            k_latosa                      , &
                                            k_root                        , &
                                            k_rtos                        , &
                                            c0_allom                      , &  ! lctlib
                                            biomass_alloc_scheme          , &
                                            bnf_scheme                    , &
                                            height                        , &
                                            veg_pool_sap_wood_carbon      , &
                                            veg_pool_heart_wood_carbon    , &
                                            dens_ind                      , &
                                            labile_carbon                 , &
                                            labile_nitrogen               , &
                                            labile_phosphorus             , &
                                            reserve_carbon                , &
                                            reserve_nitrogen              , &
                                            reserve_phosphorus            , &
                                            growth_cn                     , &
                                            growth_cp                     , &
                                            growth_np                     , &
                                            npp                           , &
                                            unit_npp                      , &
                                            transpiration                 , &
                                            unit_transpiration            , &
                                            n_fixation                    , &
                                            unit_uptake_n                 , &
                                            unit_uptake_p                 , &
                                            dphi                          , &
                                            kstar_labile                  , &
                                            w_root_lim                    , &
                                            dheight                       , &
                                            leaf2sapwood_mass_ratio       , &
                                            leaf2root_mass_ratio          , &
                                            k1_opt                        , &
                                            k2_opt )                    

    USE mo_jsb_math_constants,      ONLY: one_day, one_year, eps4, eps8
    USE mo_jsb_physical_constants,  ONLY: rhoh2o
    USE mo_veg_constants,           ONLY: min_height, itree, igrass, sm2lm_grass, w_root_lim_max, &
                                          leaf2root_min_ratio, leaf2root_max_ratio
    USE mo_pheno_constants,         ONLY: ievergreen, iraingreen, isummergreen, iperennial

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,               INTENT(in)  :: growthform                   !< land-cover-type library parameter
    REAL(wp),              INTENT(in)  :: allom_k1               , &   !< land-cover-type library parameter
                                          allom_k2               , &   !< land-cover-type library parameter
                                          sla                    , &   !< land-cover-type library parameter
                                          wood_density           , &   !< land-cover-type library parameter
                                          tau_leaf               , &   !< land-cover-type library parameter
                                          tau_fine_root          , &   !< land-cover-type library parameter
                                          k_latosa               , &   !< land-cover-type library parameter
                                          k_root                 , &   !< land-cover-type library parameter
                                          k_rtos                 , &   !< land-cover-type library parameter
                                          c0_allom                     !< land-cover-type library parameter
    CHARACTER(len=*),      INTENT(in)  :: biomass_alloc_scheme         !< from veg config
    CHARACTER(len=*),      INTENT(in)  :: bnf_scheme                   !< from veg config
    REAL(wp),              INTENT(in)  :: height                       !< height (m)
    REAL(wp),              INTENT(in)  :: veg_pool_sap_wood_carbon     !< woody biomass of an individual tree (mol C/m2)
    REAL(wp),              INTENT(in)  :: veg_pool_heart_wood_carbon   !< woody biomass of an individual tree (mol C/m2)
    REAL(wp),              INTENT(in)  :: dens_ind                     !< density of populations (number of indivuduals /m2) 
    REAl(wp),              INTENT(in)  :: labile_carbon                !< labile carbon (mol C/m2) 
    REAL(wp),              INTENT(in)  :: labile_nitrogen              !< labile nitrogen (mol N/m2)
    REAL(wp),              INTENT(in)  :: labile_phosphorus            !< labile phosphorus (mol P/m2)
    REAL(wp),              INTENT(in)  :: reserve_carbon               !< reserve carbon (mol C)
    REAL(wp),              INTENT(in)  :: reserve_nitrogen             !< reserve nitrogen (mol N)
    REAL(wp),              INTENT(in)  :: reserve_phosphorus           !< reserve phosphorus (mol P)
    REAL(wp),              INTENT(in)  :: growth_cn                    !< C:N of whole plant growth increment
    REAL(wp),              INTENT(in)  :: growth_cp                    !< C:P of whole plant growth increment
    REAL(wp),              INTENT(in)  :: growth_np                    !< N:P of whole plant growth increment
    REAL(wp),              INTENT(in)  :: npp                          !< current average NPP (micro-mol C m-2 s-1)
    REAL(wp),              INTENT(in)  :: unit_npp                     !< increase in NPP for unit increase in leaf mass (micro-mol C m-2 s-1 mol-1 leaf C)
    REAL(wp),              INTENT(in)  :: transpiration                !< average transpiration flux (kg m-2 s-1)
    REAL(wp),              INTENT(in)  :: unit_transpiration           !< increase in transpiration with unit increase in leaf mass (kg m-1 s-2 mol-1 leaf C)
    REAL(wp),              INTENT(in)  :: n_fixation                   !< N fixation (micro-mol N m-2 s-1)
    REAL(wp),              INTENT(in)  :: unit_uptake_n                !< N uptake per unit root mass (micro-mol N m-2 s-1 mol-1 root C)
    REAL(wp),              INTENT(in)  :: unit_uptake_p                !< P uptake per unit root mass (micro-mol P m-2 s-1 mol-1 root C)
    REAL(wp),              INTENT(in)  :: dphi                         !< leaf-soil water potential gradient (MPa) 
    REAL(wp),              INTENT(in)  :: kstar_labile                 !< labile pool turnover rate
    REAL(wp),              INTENT(in)  :: w_root_lim
    REAL(wp),              INTENT(out) :: dheight                      !< change in height (m) resulting from a small change in woody biomass
    REAL(wp),              INTENT(out) :: leaf2sapwood_mass_ratio      !< leaf to sapwood mass ratio 
    REAL(wp),              INTENT(out) :: leaf2root_mass_ratio         !< leaf to root mass ratio
    REAL(wp),              INTENT(out) :: k1_opt                       !< coefficient in optimized leaf:root ratio
    REAL(wp),              INTENT(out) :: k2_opt                       !< coefficient in optimized leaf:root ratio
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)             :: diameter2        !< diameter (m) resulting from a infitisemal small increase in woody biomas 
    REAL(wp)             :: k1_opt_w  , &
                            k2_opt_w
    REAL(wp)             :: k1_opt_n  , &
                            k1_opt_p
    REAL(wp)             :: n_l2r_scal, &
                            p_l2r_scal, &
                            w_l2r_scal
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_allometry'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !>0.9 init out and local variables
    !>
    dheight                 = 0.0_wp
    leaf2sapwood_mass_ratio = 0.0_wp
    leaf2root_mass_ratio    = 0.0_wp
    k1_opt                  = 0.0_wp
    k2_opt                  = 0.0_wp    
    diameter2               = 0.0_wp
    k1_opt_w                = 0.0_wp
    k2_opt_w                = 0.0_wp
    k1_opt_n                = 0.0_wp
    k1_opt_p                = 0.0_wp
    n_l2r_scal              = 0.0_wp
    p_l2r_scal              = 0.0_wp
    w_l2r_scal              = 0.0_wp

    !>1.0 case growthform
    !>
    SELECT CASE (growthform)

      !> 1.1 trees
      !! 
      CASE (itree)
        ! Diameter and Height calculated assuming that 
        !   all woody biomass forms a cylinder: diameter = SQRT(4 * wood volume / height / pi)
        !   height = allom_k1 * diameter ^ allom_k2
        IF(dens_ind > eps8) THEN
           ! Calculate change of height given a small change in woody biomass
           diameter2 = calc_diameter_from_woody_biomass(allom_k1                      , &  ! lctlib
                                                        allom_k2                      , &
                                                        wood_density                  , &  ! lctlib
                                                        veg_pool_sap_wood_carbon+eps4 , & 
                                                        veg_pool_heart_wood_carbon    , &
                                                        dens_ind)
           dheight   = calc_height_from_diameter(allom_k1, allom_k2, diameter2) - height
        ELSE 
           dheight   = 0.0_wp
        ENDIF
 
        ! scaling factor for leaf to wood mass, derived from pipe-model (Shinozaki et al. 1964, J J Ecol)
        ! leaf area = k_latosa * sap_wood area
        ! => Cl * sla = k * Cs / wood_density / height
        ! => Cl = k_latosa / sla / wood_density / height * Cs = leaf2sapwood_mass_ratio * Cs / height
        leaf2sapwood_mass_ratio = k_latosa / sla / wood_density

        SELECT CASE (biomass_alloc_scheme)
        CASE("optimal")
          !> Calculate allocation fractions assuming optimal relationship between leaf and root allocation
          !! potential carbon growth rate = potential nittrogen growth rate * plant cn ratio
          !! NPP = gN * R * growth_cn
          !! Or similar relationships for P
          !! Difference between BNF schemes - for optimal scheme fixation depends on root mass
          SELECT CASE (bnf_scheme)
          CASE("optimal")
            k1_opt_n = (npp * tau_leaf * one_year * one_day + (labile_carbon - (labile_nitrogen) * growth_cn) * 1.e6_wp) / & 
              & MAX((unit_uptake_n + n_fixation) * growth_cn * tau_fine_root * one_year * one_day, eps8)
            k1_opt_p = (npp * tau_leaf * one_year * one_day + (labile_carbon - (labile_phosphorus) * growth_cp) * 1.e6_wp) / & 
              & MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            IF ( k1_opt_n >= k1_opt_p) THEN 
              ! if N is more limiting than P
              k1_opt = k1_opt_n 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX((unit_uptake_n + n_fixation) * growth_cn * tau_fine_root * one_year * one_day, eps8)
            ELSE
              k1_opt = k1_opt_p 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            ENDIF
          CASE DEFAULT 
            k1_opt_n = (npp*tau_leaf*one_year*one_day + &
                        (labile_carbon - (labile_nitrogen) * growth_cn) * 1.e6_wp - &
                        n_fixation * growth_cn * tau_fine_root * one_year * one_day) / & 
                       MAX(unit_uptake_n * growth_cn * tau_fine_root * one_year * one_day, eps8)
             
            k1_opt_p = (npp*tau_leaf*one_year*one_day + &
                        (labile_carbon - (labile_phosphorus) * growth_cp) * 1.e6_wp) / & 
                       MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            IF ( k1_opt_n >= k1_opt_p) THEN 
              ! if N is more limiting than P
              k1_opt = k1_opt_n 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX(unit_uptake_n * growth_cn * tau_fine_root * one_year * one_day, eps8)
            ELSE
              k1_opt = k1_opt_p 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            ENDIF
          END SELECT

          !> water based optimality, following Magnani 2000 
          !! R = k1 + k2*dL +k3*dS
          IF (height > 0.5_wp) THEN
            k1_opt_w = transpiration      / rhoh2o * (1._wp + c0_allom * height) / (k_root * dphi) 
            k2_opt_w = unit_transpiration / rhoh2o * (1._wp + c0_allom * height) / (k_root * dphi)
          ELSE
            ! all resistance is in the root for small trees
            k1_opt_w = transpiration      / rhoh2o * (1._wp) / (k_root * dphi) 
            k2_opt_w = unit_transpiration / rhoh2o * (1._wp) / (k_root * dphi)
          ENDIF
            
          ! Compare root increments for a hypothetical unit of leaf growth 
          IF (k1_opt + k2_opt < k1_opt_w + k2_opt_w) THEN 
            ! water is limiting 
            k1_opt = k1_opt_w 
            k2_opt = k2_opt_w 
          ENDIF
        CASE("dynamic")
          ! empirical scalling of leaf to root ratio based on nitrogen, phosphorus and water availability
          n_l2r_scal = growth_cn / MAX(eps8, labile_carbon / MAX(eps8, labile_nitrogen))
          p_l2r_scal = growth_np / MAX(eps8, labile_nitrogen / MAX(eps8, labile_phosphorus))
          w_l2r_scal = w_root_lim / w_root_lim_max
          ! actual ratio is given by the minimum of the three factors
          leaf2root_mass_ratio = MAX(leaf2root_min_ratio, MIN(leaf2root_max_ratio, &
            & k_rtos * leaf2sapwood_mass_ratio * MIN(MIN(n_l2r_scal, p_l2r_scal), w_l2r_scal)))
        CASE("fixed")
          ! scaling factor for root to leaf mass
          ! ad 1: derived from hydraulic homeostatis of root and wood hydraulic resistance (Magnani et al. 2000, PCE)
          ! Note: currently still lacks a moisture response!
          ! sapwood mass / height =  k_rtos * fine_root_mass
          ! Cr = Cs / k_rtos / height
          ! let Cs = Cl / leaf2sapwood_mass_ratio * height
          ! => Cl = Cr * leaf2sapwood_mass_ratio * k_rtos = Cl = Cr * leaf2root_mass_ratio 
          ! ad 2: nitrogen/phosphorus sensitivity still needed
          leaf2root_mass_ratio = k_rtos * leaf2sapwood_mass_ratio
        END SELECT  ! biomass_alloc_scheme

      !> 1.2 grass
      !! 
      CASE (igrass)
        ! leaf to stem mass ratio
        leaf2sapwood_mass_ratio = 1._wp / sm2lm_grass
        ! leaf:root ratio
        SELECT CASE (biomass_alloc_scheme)
        CASE("optimal")
          SELECT CASE (bnf_scheme)
          CASE("optimal")
            k1_opt_n = (npp * tau_leaf * one_year * one_day + (labile_carbon - (labile_nitrogen) * growth_cn) * 1.e6_wp) / & 
              & MAX((unit_uptake_n + n_fixation) * growth_cn * tau_fine_root * one_year * one_day, eps8)
            k1_opt_p = (npp * tau_leaf * one_year * one_day + (labile_carbon - (labile_phosphorus) * growth_cp) * 1.e6_wp) / & 
              & MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            IF ( k1_opt_n >= k1_opt_p) THEN 
              ! if N is more limiting than P
              k1_opt = k1_opt_n 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX((unit_uptake_n + n_fixation) * growth_cn * tau_fine_root * one_year * one_day, eps8)
            ELSE
              k1_opt = k1_opt_p 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            ENDIF
          CASE DEFAULT 
            k1_opt_n = (npp * tau_leaf * one_year * one_day + (labile_carbon - (labile_nitrogen) * growth_cn) * 1.e6_wp - &
                        n_fixation * growth_cn * tau_fine_root * one_year * one_day) / & 
                       MAX(MAX(unit_uptake_n, 0.05_wp) * growth_cn*tau_fine_root*one_year*one_day, eps8)
            k1_opt_p = (npp * tau_leaf * one_year * one_day + (labile_carbon - (labile_phosphorus) * growth_cp) * 1.e6_wp) / & 
                       MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            IF ( k1_opt_n >= k1_opt_p) THEN 
              ! if N is more limiting than P
              k1_opt = k1_opt_n 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX(unit_uptake_n * growth_cn * tau_fine_root * one_year * one_day, eps8)
            ELSE
              k1_opt = k1_opt_p 
              k2_opt = (unit_npp * tau_leaf * one_year * one_day) / &
                & MAX(unit_uptake_p * growth_cp * tau_fine_root * one_year * one_day, eps8)
            ENDIF
          END SELECT ! bnf_scheme

          !> water based optimality, following Magnani 2000 
          k1_opt_w = transpiration      / (k_root * dphi) 
          k2_opt_w = unit_transpiration / (k_root * dphi) 

          ! Compare root increments for a hypothetical unit of leaf growth 
          IF (k1_opt + k2_opt < k1_opt_w + k2_opt_w) THEN 
            ! water is limiting 
            k1_opt = k1_opt_w 
            k2_opt = k2_opt_w 
          ENDIF 
        CASE("dynamic")
            ! empirical scalling of leaf to root ratio based on nitrogen, phosphorus and water availability
           n_l2r_scal = growth_cn / MAX(eps8, labile_carbon / MAX(eps8, labile_nitrogen))
           p_l2r_scal = growth_np / MAX(eps8, labile_nitrogen / MAX(eps8, labile_phosphorus))
           w_l2r_scal = w_root_lim / w_root_lim_max
           ! actual ratio is given by the minimum of the three factors
           leaf2root_mass_ratio = MAX(leaf2root_min_ratio, MIN(leaf2root_max_ratio, &
                sm2lm_grass * k_rtos * MIN(MIN(n_l2r_scal, p_l2r_scal), w_l2r_scal)))
        CASE("fixed")
          leaf2root_mass_ratio = sm2lm_grass * k_rtos 
        END SELECT ! biomass_alloc_scheme

    ENDSELECT ! SELECT CASE (growthform)

  END SUBROUTINE calc_allometry 


  !----------------------------------------------------------------------------------------------------------
  ! 4.1 Sub Task of update_veg_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Calculate target leaf area index 
  !! 
  !! Input:\n
  !!    (1) for trees: height and sap wood mass, as well as lctlib parameters \n
  !!        for grasses: root mass and growing season leaf to root ratio \n 
  !!  
  !! output: target size of leaf area index 
  !----------------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_target_lai(growthform               , &  ! lctlib
                                     sla                      , &
                                     k_latosa                 , &  
                                     wood_density             , &  ! lctlib
                                     leaf2root_mass_ratio     , &
                                     height                   , & 
                                     veg_pool_fine_root_carbon, &
                                     veg_pool_sap_wood_carbon)    RESULT (target_lai)

    USE mo_veg_constants,           ONLY: min_lai, igrass, itree
    USE mo_jsb_math_constants,      ONLY: eps8

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,  INTENT(in) :: growthform                          !< plant's growth form
    REAL(wp), INTENT(in) :: sla                       , &       !< specific leaf area (m2 / mol C)
                            k_latosa                  , &       !< leaf to sapwood area ratio
                            wood_density                        !< wood density (mol C / m3)
    REAL(wp), INTENT(in) :: leaf2root_mass_ratio      , &       !< leaf to fine root mass ratio
                            height                    , &       !< height (m)
                            veg_pool_fine_root_carbon , &       !< fine root mass (mol C / m2)
                            veg_pool_sap_wood_carbon            !< sap wood mass (mol C / m2
    REAL(wp)             :: target_lai 
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_target_lai'


    SELECT CASE(growthform)
      CASE (itree)
        IF (height > eps8) THEN
          ! full leaf area implied by current sapwood area
          target_lai = MAX(veg_pool_sap_wood_carbon * k_latosa / wood_density / height, min_lai)
        ELSE
          target_lai = min_lai
        END IF
 
      CASE (igrass)
          ! full leaf area implied by current fine_root mass and leaf to root ratio 
          target_lai = MAX(veg_pool_fine_root_carbon * leaf2root_mass_ratio * sla, min_lai)

    END SELECT

  END FUNCTION calc_target_lai


  !-----------------------------------------------------------------------------------------------------
  ! 4.1 Sub Task of update_veg_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Calculate target size of the labile and reserve pools for each element 
  !! 
  !! Input:\n
  !!    (1) weekly GPP and maintenance respiration, and N,P growth requirements to calculate
  !!        target labile pool size;\n
  !!    (2) target LAI, leaf-to-root mass ratio and target C:N:P of foliage and fine roots to calculate
  !!        target reserve pool size;\n
  !!    (3) current leaf, fine root and sapwood mass to calculate physical limit of reserve pool size.
  !!  
  !! Output: target sizes of labile and reserve pools 
  !----------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_target_labile_reserve_pool_size(nc                            , &
                                                  lctlib_growthform             , &
                                                  lctlib_fstore_target          , &
                                                  lctlib_sla                    , &
                                                  lctlib_tau_leaf               , &
                                                  lctlib_tau_fine_root          , &
                                                  gpp_mavg                      , &
                                                  maint_respiration_mavg        , & 
                                                  growth_req_n_mavg             , &
                                                  growth_req_p_mavg             , &
                                                  leaf2root_mass_ratio_mavg     , &
                                                  target_cn_leaf                , &
                                                  target_cn_fine_root           , &
                                                  target_np_leaf                , &
                                                  target_np_fine_root           , & 
                                                  target_lai                    , &
                                                  veg_pool_leaf_carbon          , & 
                                                  veg_pool_fine_root_carbon     , & 
                                                  veg_pool_sap_wood_carbon      , & 
                                                  target_labile_pool_carbon     , & 
                                                  target_labile_pool_nitrogen   , & 
                                                  target_labile_pool_phosphorus , & 
                                                  target_reserve_pool_carbon    , & 
                                                  target_reserve_pool_nitrogen  , & 
                                                  target_reserve_pool_phosphorus )  


    USE mo_jsb_math_constants,    ONLY: one_day, eps8
    USE mo_veg_constants,         ONLY: tau_labile, fresp_growth, resorp_fract_leaf, igrass, itree, &
                                        fstore_leaf_max, fstore_fine_root_max, fstore_sap_wood_max

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                         INTENT(in)     :: nc                                !< dimensions
    INTEGER,                         INTENT(in)     :: lctlib_growthform                 !< lctlib parameter
    REAL(wp),                        INTENT(in)     :: lctlib_fstore_target              !< lctlib parameter
    REAL(wp),                        INTENT(in)     :: lctlib_sla                        !< lctlib parameter
    REAL(wp),                        INTENT(in)     :: lctlib_tau_leaf                   !< lctlib parameter
    REAL(wp),                        INTENT(in)     :: lctlib_tau_fine_root              !< lctlib parameter
    REAL(wp), DIMENSION(nc),         INTENT(in)     :: gpp_mavg                      , & !< last week's GPP (micro-mol CO2/m2/s) 
                                                       maint_respiration_mavg        , & !< last week's maintenance respiration (micro-mol CO2/m2/s)
                                                       growth_req_n_mavg             , & !< last week's requirement for growth (mol N / mol C)
                                                       growth_req_p_mavg             , & !< last week's requirement for growth (mol P / mol N)
                                                       leaf2root_mass_ratio_mavg     , & !< leaf to root mass ratio over the lifetime of fine roots
                                                       target_cn_leaf                , & !< target C:N ratio for leaves (mol C / mol N) 
                                                       target_cn_fine_root           , & !< target C:N ratio for fine roots (mol C / mol N)
                                                       target_np_leaf                , & !< target N:P ratio for leaves (mol N / mol P)
                                                       target_np_fine_root           , & !< target N:P ratio for fine roots (mol N / mol P)
                                                       target_lai                    , & !< target LAI 
                                                       veg_pool_leaf_carbon          , & !< current leaf carbon pool (mol C m-2) 
                                                       veg_pool_fine_root_carbon     , & !< current fine root carbon pool (mol C m-2)
                                                       veg_pool_sap_wood_carbon          !< current sap wood carbon pool (mol C m-2)
    REAL(wp), DIMENSION(nc),         INTENT(out)    :: target_labile_pool_carbon     , & !< desirable size (mol C/m2) of the labile pool
                                                       target_labile_pool_nitrogen   , & !< desirable size (mol N/m2) of the labile pool
                                                       target_labile_pool_phosphorus , & !< desirable size (mol P/m2) of the labile pool
                                                       target_reserve_pool_carbon    , & !< desirable size (mol C/m2) of the reserve pool
                                                       target_reserve_pool_nitrogen  , & !< desirable size (mol N/m2) of the reserve pool
                                                       target_reserve_pool_phosphorus    !< desirable size (mol P/m2) of the reserve pool
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)          :: reserve_pool_carbon_max                      , & ! physical limit of the reserve pool (mol C/m2) 
                                        hlp1                                             ! helper variable
    CHARACTER(len=*), PARAMETER      :: routine = TRIM(modname)//':calc_target_labile_reserve_pool_size'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    
    !> 1.0 target size of labile pool
    !! 
    !! This is taken as the average daily input or maintenance requirements of the last week 
    !!  scaled to the default turnover time of the labile pool
    target_labile_pool_carbon(:) = MAX(gpp_mavg(:) , maint_respiration_mavg(:)) * one_day / 1.e6_wp * tau_labile

    !! for N and P this corresponds to the current growth C:N:P in the target labile pool 
    target_labile_pool_nitrogen(:)   = target_labile_pool_carbon(:)    * growth_req_n_mavg(:)
    target_labile_pool_phosphorus(:) = target_labile_pool_nitrogen(:)  * growth_req_p_mavg(:)


    !> 2.0 target size reserve pool
    !!
    !! For the reserve pool, this is the C:N:P needed to regrow one year worth of folige and roots,
    !! accounting for losses to respiration and the nutrients currently required to grow these
    !!
    target_reserve_pool_carbon(:)     = lctlib_fstore_target * target_lai(:) / lctlib_sla * & 
                                        ( MIN(1.0_wp,1.0_wp/lctlib_tau_leaf) + &
                                         1._wp/leaf2root_mass_ratio_mavg(:)/lctlib_tau_fine_root ) * &
                                        (1._wp + fresp_growth)
    target_reserve_pool_nitrogen(:)   = lctlib_fstore_target * target_lai / lctlib_sla * & 
                                        ( MIN(1.0_wp, 1.0_wp/lctlib_tau_leaf) / &
                                         target_cn_leaf(:) + 1._wp/leaf2root_mass_ratio_mavg(:) / &
                                         lctlib_tau_fine_root / target_cn_fine_root(:)) * & 
                                        (1._wp + fresp_growth)
    target_reserve_pool_phosphorus(:) = lctlib_fstore_target * target_lai(:) / lctlib_sla * & 
                                        ( MIN(1.0_wp, 1.0_wp/lctlib_tau_leaf) / &
                                         target_cn_leaf(:) / target_np_leaf(:) + & 
                                         1._wp / leaf2root_mass_ratio_mavg(:)/ lctlib_tau_fine_root / &
                                         target_cn_fine_root(:) / target_np_fine_root(:)) * & 
                                        (1._wp + fresp_growth)

    !> 2.1 limit target size of the reserve pool to the size of the storage organs
    !!
    reserve_pool_carbon_max(:) = fstore_leaf_max      * veg_pool_leaf_carbon(:)      + & 
                                 fstore_fine_root_max * veg_pool_fine_root_carbon(:) + & 
                                 fstore_sap_wood_max  * veg_pool_sap_wood_carbon(:) 

    hlp1(:)                           = MIN(1.0_wp, reserve_pool_carbon_max(:) / target_reserve_pool_carbon(:))
    target_reserve_pool_carbon(:)     = target_reserve_pool_carbon(:)     * hlp1(:)
    target_reserve_pool_nitrogen(:)   = target_reserve_pool_nitrogen(:)   * hlp1(:)
    target_reserve_pool_phosphorus(:) = target_reserve_pool_phosphorus(:) * hlp1(:)

  END SUBROUTINE calc_target_labile_reserve_pool_size


  !----------------------------------------------------------------------------------------------------------
  ! 4.2 Sub Task of update_veg_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Calculate flux between reserve and labile pools for each element 
  !! 
  !! Input:\n
  !!    (1) target sizes of the labile and reserve pools
  !!    (2) current LAI, labile and reserve pools, growing season status and potential turnover rate
  !!          of the labile pool;\n    
  !!  
  !! output: net flux between reserve and labile pools 
  !----------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_reserve_useage( nc                            , &
                                  growing_season                , &
                                  kstar_labile                  , &
                                  target_labile_pool_carbon     , & 
                                  target_labile_pool_nitrogen   , & 
                                  target_labile_pool_phosphorus , & 
                                  target_reserve_pool_carbon    , & 
                                  target_reserve_pool_nitrogen  , & 
                                  target_reserve_pool_phosphorus, & 
                                  target_lai                    , &
                                  lai                           , &
                                  veg_pool_labile               , &
                                  veg_pool_reserve              , &
                                  veg_reserve_use               ) 

    USE mo_kind,                  ONLY: wp
    USE mo_jsb_math_constants,    ONLY: one_day, dtime, eps8
    USE mo_jsb_impl_constants,    ONLY: test_false_true
    USE mo_veg_constants,         ONLY: tau_labile, &
                                        lambda_phiphen,  k_phiphen, &
                                        lambda_phimaint, k_phimaint, &
                                        lambda_phimaint_nut, k_phimaint_nut, &
                                        lambda_phistore, k_phistore, &
                                        k_phi_interact, &
                                        min_lai, target_lai_max
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_veg
#else
    dsl4jsb_Use_pools(VEG_) t_veg_elements_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                          INTENT(in)    :: nc                                !< dimensions
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: growing_season                    !< does the plant intend to grow?
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: kstar_labile                  , & !< potential labile pool turnover (1/day) 
                                                       target_labile_pool_carbon     , & ! desirable size (mol C/m2) of the labile pool
                                                       target_labile_pool_nitrogen   , & ! desirable size (mol N/m2) of the labile pool
                                                       target_labile_pool_phosphorus , & ! desirable size (mol P/m2) of the labile pool
                                                       target_reserve_pool_carbon    , & ! desirable size (mol C/m2) of the reserve pool
                                                       target_reserve_pool_nitrogen  , & ! desirable size (mol N/m2) of the reserve pool
                                                       target_reserve_pool_phosphorus, & ! desirable size (mol P/m2) of the reserve pool
                                                       target_lai                    , & !< target LAI 
                                                       lai                               !< current LAI 
#ifdef __QUINCY_STANDALONE__
    TYPE(t_veg_element_a),            INTENT(in)    :: veg_pool_labile                   !< the plant's labile pool size
    TYPE(t_veg_element_a),            INTENT(in)    :: veg_pool_reserve                  !< the plant's reserve pool size
    TYPE(t_veg_element_a),            INTENT(inout) :: veg_reserve_use                   !< current resource use rate (mol/time step)
                                                                                         !!  (+ve: net flux to labile / -ve: net flux to reserve)
#else
    TYPE(t_veg_elements_pool),        INTENT(in)    :: veg_pool_labile                   !< the plant's labile pool size
    TYPE(t_veg_elements_pool),        INTENT(in)    :: veg_pool_reserve                  !< the plant's reserve pool size
    TYPE(t_veg_elements_pool),        INTENT(inout) :: veg_reserve_use                   !< current resource use rate (mol/time step)
                                                                                         !!  (+ve: net flux to labile / -ve: net flux to reserve)
#endif
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                         :: phi_phenology                 , & ! relative pull from reserve due to phenology
                                                       phi_maintenance               , & ! relative pull from reserve due to maintenance resp. 
                                                       phi_storage                   , & ! relative pull from labile due to low reserve
                                                       phi_phenology_base            , & ! helper variable for using for all elements
                                                       hlp1                              ! helper variable
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_reserve_useage'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 carbon reserve flow
    !! 
    !! calculate the pull from reserve to labile (phi_phenology,phi_maintenance) or
    !! 
    !! push from labile to reserve (phi_storage) depending on the difference of the LAI/labile/reserve from their targets
    !! 
    !! phi_storage is adjusted by phi_phenology to avoid storage building to interfere with flushing (push phi_storage towards 
    !!  reserve at times of low LAI or labile pool to give priority to growth and maintentance)
    hlp1(:) = MIN( MAX(min_lai, target_lai(:)), target_lai_max)
    WHERE(lai(:) < hlp1(:) .AND. growing_season(:) > test_false_true)
      phi_phenology(:) = exp(-(lambda_phiphen * lai(:) / hlp1(:) )**k_phiphen) - exp(-(lambda_phiphen)**k_phiphen) 
    ELSEWHERE
      phi_phenology(:) = 0._wp
    END WHERE

    ! Save value for use for N and P calculations
    phi_phenology_base(:) = phi_phenology(:)

    hlp1(:) = MAX(0.0_wp, veg_pool_labile%carbon(:))    ! for safety, catch potential negative values
    WHERE(target_labile_pool_carbon(:) > eps8 .AND. hlp1(:) < target_labile_pool_carbon(:))
      phi_maintenance(:) = exp(-(lambda_phimaint * hlp1(:) / target_labile_pool_carbon(:))**k_phimaint) & 
                           - exp(-(lambda_phimaint)**k_phimaint)
    ELSEWHERE
      phi_maintenance(:) = 0._wp
    END WHERE

    hlp1(:) = MAX(0.0_wp, veg_pool_reserve%carbon(:))   ! for safety, catch potential negative values
    WHERE(target_reserve_pool_carbon(:) > eps8 .AND. hlp1(:) < target_reserve_pool_carbon(:))
      phi_storage(:) = exp(-(lambda_phistore * hlp1(:) / target_reserve_pool_carbon(:))**k_phistore) & 
                       - exp(-(lambda_phistore)**k_phistore)
    ELSEWHERE
      phi_storage(:) = 0._wp
    END WHERE

    phi_storage(:) = phi_storage(:) * ( 1.0_wp - phi_phenology(:)) 
    WHERE(phi_maintenance(:) > k_phi_interact) 
      phi_storage(:) = phi_storage(:) * 1._wp / (1._wp - k_phi_interact) * ( 1.0_wp - phi_maintenance(:)) 
    END WHERE

    !> 1.1 convert phi to turnover times (1/timestep)
    !! 
    !! taking account that flushing can only happen if the meristem is active
    phi_phenology(:)    = phi_phenology(:)    * kstar_labile(:) / one_day * dtime
    phi_maintenance(:)  = phi_maintenance(:)  / tau_labile      / one_day * dtime
    phi_storage(:)      = phi_storage(:)      / tau_labile      / one_day * dtime

    veg_reserve_use%carbon(:)   = ((phi_phenology(:) + phi_maintenance(:)) * veg_pool_reserve%carbon(:)   - & 
                                    phi_storage(:)   * veg_pool_labile%carbon(:)) 
    veg_reserve_use%carbon13(:) = ((phi_phenology(:) + phi_maintenance(:)) * veg_pool_reserve%carbon13(:) - & 
                                    phi_storage(:)   * veg_pool_labile%carbon13(:)) 
    veg_reserve_use%carbon14(:) = ((phi_phenology(:) + phi_maintenance(:)) * veg_pool_reserve%carbon14(:) - & 
                                    phi_storage(:)   * veg_pool_labile%carbon14(:)) 


    !> 2.0 net transfer of N from labile to reserve
    !!
    !> 2.1 supply of reserve N to labile pool during growing season 
    !!
    phi_phenology(:) = phi_phenology_base(:)
    hlp1(:)          = MAX(0.0_wp, veg_pool_labile%nitrogen(:))        ! for safety, catch potential negative values
    WHERE(target_labile_pool_nitrogen(:) > eps8    .AND. & 
          hlp1(:) < target_labile_pool_nitrogen(:) .AND. & 
          growing_season(:) > test_false_true)
      phi_maintenance(:) = exp(-(lambda_phimaint_nut * hlp1(:) / target_labile_pool_nitrogen(:))**k_phimaint_nut) & 
                           - exp(-(lambda_phimaint_nut)**k_phimaint_nut)
    ELSEWHERE
      phi_maintenance(:) = 0._wp
    ENDWHERE

    !> 2.2 build-up of reserve N (composed of a pull of N to the reserve and a push from the labile to avoid 
    !!     too strong accumulation of N in the labile)
    !! 
    WHERE(target_reserve_pool_nitrogen(:) > eps8 .AND. veg_pool_reserve%nitrogen(:) < target_reserve_pool_nitrogen(:))
      phi_storage(:) = exp(-(lambda_phistore * veg_pool_reserve%nitrogen(:) / target_reserve_pool_nitrogen(:))**k_phistore) & 
                       - exp(-(lambda_phistore)**k_phistore)
    ELSEWHERE
      phi_storage(:) = 0._wp      
    ENDWHERE
    WHERE(target_labile_pool_nitrogen(:) > eps8 .AND. veg_pool_labile%nitrogen(:) > target_labile_pool_nitrogen(:))
      phi_storage(:) = MAX(phi_storage(:), &
                           MIN((veg_pool_labile%nitrogen(:) - target_labile_pool_nitrogen(:)) / target_labile_pool_nitrogen(:), &
                               1.0_wp)) 
    ENDWHERE

    !> 2.3 Limit push from labile to reserve to maximum storage capacity, and account for preference of bud-burst over storage
    !!
    phi_storage(:) = phi_storage(:) * ( 1.0_wp - phi_phenology(:)) 
    hlp1(:)        = veg_pool_reserve%nitrogen(:) + phi_storage(:) * veg_pool_labile%nitrogen(:) / tau_labile / one_day * dtime
    WHERE(veg_pool_labile%nitrogen(:) > eps8 .AND. hlp1(:) > target_reserve_pool_nitrogen(:))
       phi_storage(:) = MAX(0.0_wp,(target_reserve_pool_nitrogen(:) - veg_pool_reserve%nitrogen(:))) / &
                            (veg_pool_labile%nitrogen(:) / tau_labile / one_day * dtime)
    ENDWHERE

    !> 2.4 Limit push from reserve to labile to maintain a minuimum N concentration in the reserve
    !! 
    WHERE(phi_storage(:) > k_phi_interact) 
      phi_maintenance(:) = phi_maintenance(:) * 1._wp / (1._wp - k_phi_interact) * ( 1.0_wp - phi_storage(:)) 
    END WHERE

    veg_reserve_use%nitrogen(:)    = ((phi_maintenance(:) + phi_phenology(:)) * veg_pool_reserve%nitrogen(:) - & 
                                      phi_storage(:) * veg_pool_labile%nitrogen(:)) / tau_labile / one_day * dtime
    veg_reserve_use%nitrogen15(:)  = ((phi_maintenance(:) + phi_phenology(:)) * veg_pool_reserve%nitrogen15(:) - & 
                                      phi_storage(:) * veg_pool_labile%nitrogen15(:)) / tau_labile / one_day * dtime


    !> 3.0 net transfer of P from labile to reserve
    !! 
    !> 3.1 supply of reserve P to labile pool during growing season 
    !!
    phi_phenology(:) = phi_phenology_base(:)
    hlp1(:)          = MAX(0.0_wp, veg_pool_labile%phosphorus(:))      ! for safety, catch potential negative values
    WHERE(target_labile_pool_phosphorus(:) > eps8    .AND. & 
          hlp1(:) < target_labile_pool_phosphorus(:) .AND. & 
          growing_season(:) > test_false_true)
      phi_maintenance(:) = exp(-(lambda_phimaint_nut * hlp1(:) / target_labile_pool_phosphorus(:))**k_phimaint_nut) & 
                           - exp(-(lambda_phimaint_nut)**k_phimaint_nut)
    ELSEWHERE
      phi_maintenance(:) = 0._wp
    END WHERE

    !> 3.2 build-up of reserve P 
    !! 
    !! (composed of a pull of N to the reserve and a push from the labile to avoid too strong accumulation of P in the labile)
    WHERE(target_reserve_pool_phosphorus(:) > eps8 .AND. veg_pool_reserve%phosphorus(:) < target_reserve_pool_phosphorus(:))
      phi_storage(:) = exp(-(lambda_phistore * veg_pool_reserve%phosphorus(:) / target_reserve_pool_phosphorus(:))**k_phistore) & 
                       - exp(-(lambda_phistore)**k_phistore)
    ELSEWHERE
      phi_storage(:) = 0._wp      
    END WHERE
    WHERE(target_labile_pool_phosphorus(:) > eps8 .AND. veg_pool_labile%phosphorus(:) > target_labile_pool_phosphorus(:))
      phi_storage(:) = MAX(phi_storage(:), &
                           MIN((veg_pool_labile%phosphorus(:) - target_labile_pool_phosphorus(:)) / &
                               target_labile_pool_phosphorus(:), 1.0_wp)) 
    END WHERE

    !> 3.3 Limit push from labile to reserve to maximum storage capacity, and account for preference of bud-burst over storage
    !!
    phi_storage(:) = phi_storage(:) * ( 1.0_wp - phi_phenology(:)) 
    hlp1(:)        = veg_pool_reserve%phosphorus(:) + phi_storage(:) * veg_pool_labile%phosphorus(:) / tau_labile / one_day * dtime
    WHERE(veg_pool_labile%phosphorus(:) > eps8 .AND. hlp1(:) > target_reserve_pool_phosphorus(:))
      phi_storage(:) = MAX(0.0_wp,(target_reserve_pool_phosphorus(:) - veg_pool_reserve%phosphorus(:))) / &
                           (veg_pool_labile%phosphorus(:) / tau_labile / one_day * dtime)
    ENDWHERE

    !> 3.4 Limit push from reserve to labile to maintain a minuimum P concentration in the reserve
    !! 
    WHERE(phi_storage(:) > k_phi_interact) 
      phi_maintenance(:) = phi_maintenance(:) * 1._wp / (1._wp - k_phi_interact) * ( 1.0_wp - phi_storage(:)) 
    END WHERE

    veg_reserve_use%phosphorus(:) = ((phi_maintenance(:)+phi_phenology(:)) * veg_pool_reserve%phosphorus(:) - & 
                                     phi_storage(:) * veg_pool_labile%phosphorus(:)) / tau_labile / one_day * dtime

  END SUBROUTINE calc_reserve_useage


  !----------------------------------------------------------------------------------------------------------
  ! 4.3 Sub Task of update_veg_growth
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Calculate sinklimitation reduction factor for photosynthesis if reduced meristm activity or 
  !! stoichiometric constraints lead to an accumulation of labile C or a significant undersupply of 
  !! N and P 
  !! 
  !! Input:\n
  !!    (1) target and current size of the labile pool
  !!  
  !! output: beta factor for photosynthesis calculation 
  !----------------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_beta_sinklimitation(veg_pool_labile_carbon, &
                                              veg_pool_labile_nitrogen, &
                                              veg_pool_labile_phosphorus, &
                                              target_labile_pool_carbon, &
                                              target_labile_pool_nitrogen, &
                                              target_labile_pool_phosphorus) RESULT(beta_sinklim_ps)

    USE mo_kind,                          ONLY: wp
    USE mo_jsb_math_constants,            ONLY: eps8
    USE mo_veg_constants,                 ONLY: beta_sinklim_ps_min, k_sinklim_ps, lambda_sinklim_ps, k_cnp_sinklim_ps

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp), INTENT(in) :: veg_pool_labile_carbon, &
                            veg_pool_labile_nitrogen, &
                            veg_pool_labile_phosphorus, &
                            target_labile_pool_carbon, & 
                            target_labile_pool_nitrogen, &
                            target_labile_pool_phosphorus
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)             :: hlp1, hlp2 
    REAL(wp)             :: beta_sinklim_ps
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_beta_sinklimitation'


    !> 1.0 Sink limitation due to labile C accumulation
    !!
    IF(target_labile_pool_carbon > eps8 .AND. veg_pool_labile_carbon > target_labile_pool_carbon)THEN
        hlp1            = (veg_pool_labile_carbon - target_labile_pool_carbon) / target_labile_pool_carbon
        beta_sinklim_ps = beta_sinklim_ps_min + (1._wp - beta_sinklim_ps_min) * &
                             exp(-(lambda_sinklim_ps * hlp1)**k_sinklim_ps)
    ELSE 
        beta_sinklim_ps = 1.0_wp
    ENDIF

    !> 1.1 Sink limitation due to low NP content of labile pool
    !!
    IF(veg_pool_labile_carbon > eps8 .AND. target_labile_pool_carbon > eps8)THEN
       hlp1 = (veg_pool_labile_nitrogen / veg_pool_labile_carbon) / & 
              ( target_labile_pool_nitrogen / (k_cnp_sinklim_ps * target_labile_pool_carbon))
       hlp2 = (veg_pool_labile_phosphorus / veg_pool_labile_carbon) / &
              ( target_labile_pool_phosphorus / (k_cnp_sinklim_ps * target_labile_pool_carbon)) 

       ! take the minimum of the limit imposed by either N or P
       beta_sinklim_ps = MAX(beta_sinklim_ps_min,beta_sinklim_ps * MIN(MIN(hlp1,hlp2), 1.0_wp))
    ENDIF

  END FUNCTION calc_beta_sinklimitation


  !-----------------------------------------------------------------------------------------------------
  ! 5.0 Sub Task of update_veg_growth
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculates allocation fractions for each tissue type and element
  !! 
  !! Input: target stoichiometry of the tissue types, plant allometry, and current pool sizes
  !! 
  !! output:\n 
  !!         (1) Allocation coefficients for the distribution of biomass for growth;\n
  !!         (2) total N and P required to sustain a unit of growth (in mol); \n
  !!
  !! @todo: calc_allocation_fraction: check calculation of b_inc. Currently a fixed value for stability,
  !!        Needs to be adjusted to better match prescriebed allometry and to account for falloc_fruit
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_allocation_fraction(nc                          , &
                                      lctlib_growthform           , &
                                      lctlib_k2_fruit_alloc       , &
                                      lctlib_k_crtos              , &
                                      biomass_alloc_scheme        , &
                                      dheight                     , &
                                      kstar_labile                , &
                                      leaf_carbon                 , &
                                      sap_wood_carbon             , &
                                      coarse_root_carbon          , &
                                      root_carbon                 , &
                                      labile_carbon               , &
                                      height                      , &
                                      leaf2sapwood_mass_ratio     , &
                                      leaf2root_mass_ratio        , &
                                      k1_opt                      , &
                                      k2_opt                      , &
                                      target_cn_leaf              , &
                                      target_cn_fine_root         , &
                                      target_cn_coarse_root       , & 
                                      target_cn_sap_wood          , &
                                      target_cn_fruit             , &
                                      target_np_leaf              , &
                                      target_np_fine_root         , & 
                                      target_np_coarse_root       , &
                                      target_np_sap_wood          , &
                                      target_np_fruit             , &
                                      veg_reserve_use_carbon      , &
                                      veg_falloc                  , &       
                                      growth_req_n                , &     
                                      growth_req_p )       

    USE mo_jsb_math_constants,      ONLY: dtime, one_day, eps8
    USE mo_veg_constants,           ONLY: itree, igrass, &
                                          lambda_fruit_alloc, k1_fruit_alloc, k3_fruit_alloc, k4_fruit_alloc
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_veg
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                          INTENT(in)    :: nc                         !< dimensions
    INTEGER,                          INTENT(in)    :: lctlib_growthform          !< lctlib parameter
    REAL(wp),                         INTENT(in)    :: lctlib_k2_fruit_alloc      !< lctlib parameter
    REAL(wp),                         INTENT(in)    :: lctlib_k_crtos             !< lctlib parameter
    CHARACTER(len=*),                 INTENT(in)    :: biomass_alloc_scheme       !< from veg config
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: dheight                 ,& !< derivative of height per unit sapwood mass change 
                                                       kstar_labile            ,& !< potential labile pool turnover rate (1/day) without stoichiometric limit
                                                       leaf_carbon             ,& !< carbon in leaves (mol/individual)
                                                       sap_wood_carbon         ,& !< carbon in sap wood (mol/individual)
                                                       coarse_root_carbon      ,& !< carbon in heart wood (mol/individual)
                                                       root_carbon             ,& !< carbon in roots (mol/individual) 
                                                       labile_carbon           ,& !< carbon in labile pool (mol/individual) 
                                                       height                  ,& !< plant height  
                                                       leaf2sapwood_mass_ratio ,& !< leaf to shoot mass ratio
                                                       leaf2root_mass_ratio    ,& !< leaf to root mass ratio
                                                       k1_opt                  ,& !< coefficient for optimal leaf:root ratio
                                                       k2_opt                  ,& !< coefficient for optimal leaf:root ratio
                                                       target_cn_leaf          ,& !< target C:N ratio for tissue type 
                                                       target_cn_fine_root     ,& !< target C:N ratio for tissue type 
                                                       target_cn_coarse_root   ,& !< target C:N ratio for tissue type 
                                                       target_cn_sap_wood      ,& !< target C:N ratio for tissue type 
                                                       target_cn_fruit         ,& !< target C:N ratio for tissue type 
                                                       target_np_leaf          ,& !< target N:P ratio for tissue type 
                                                       target_np_fine_root     ,& !< target N:P ratio for tissue type 
                                                       target_np_coarse_root   ,& !< target N:P ratio for tissue type 
                                                       target_np_sap_wood      ,& !< target N:P ratio for tissue type 
                                                       target_np_fruit         ,& !< target N:P ratio for tissue type 
                                                       veg_reserve_use_carbon     !< current reserve usage (mol C / m2 / timestep)
#ifdef __QUINCY_STANDALONE__
    TYPE(t_veg_pool_a),               INTENT(inout) :: veg_falloc                 !< allocation coefficients 
#else
    TYPE(t_veg_tissues_pool),         INTENT(inout) :: veg_falloc                 !< allocation coefficients 
#endif
    REAL(wp), DIMENSION(nc),          INTENT(out)   :: growth_req_n            ,& !< N requirement for unit growth in C and N
                                                       growth_req_p               !< P requirement for unit growth in C and N
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                     :: b_inc,                   & ! total and compartment-wise carbon increment (mol/dtime)
                                                   all_inc,                 &
                                                   leaf_carbon_inc,         &
                                                   root_carbon_inc,         &
                                                   sap_wood_carbon_inc,     &
                                                   coarse_root_carbon_inc,  & 
                                                   leaf2wood_mass_ratio,    &
                                                   wood_carbon,             &
                                                   hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_allocation_fraction'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !>0.9 init output variables
    !>
    growth_req_n(:) = 0.0_wp
    growth_req_p(:) = 0.0_wp

    !> 1.0 Allocation factors for C (this should be calculated based on various principles of plant growth and allometry)
    !! 
    !! 

    ! SZ: uncomment this to revert to fixed constant allocation
    ! veg_falloc%leaf%carbon        = 0.3_wp
    ! veg_falloc%fine_root%carbon   = 0.3_wp
    ! veg_falloc%coarse_root%carbon = 0.1_wp
    ! veg_falloc%sap_wood%carbon    = 0.2_wp
    ! veg_falloc%fruit%carbon       = 0.1_wp 

    ! 
    ! seed for allocation coefficients is the potential growth rate per time step.
    ! given by the average turnover of the labile pool 
    b_inc(:) = kstar_labile(:) * labile_carbon(:) * dtime / one_day

    !> 1.1 Calculate allocation fraction to reproduction
    !!  Fruit allocation depends on whether or not the plant is able to build reserves
    !!  as soon as the carbon balance exceeds a PFT specific offset, the maximal fruit growth rate is achieved
    hlp1(:)                    = MAX(0.0_wp, veg_reserve_use_carbon(:) / dtime * 1000000.0_wp + k3_fruit_alloc)
    veg_falloc%fruit%carbon(:) = k1_fruit_alloc + ( lctlib_k2_fruit_alloc - k1_fruit_alloc ) * & 
                                 exp(-(lambda_fruit_alloc * hlp1(:))**k4_fruit_alloc)

    !> 1.2 Calculate allocation fractions to functional organs
    !!
    !! 
    SELECT CASE(lctlib_growthform)
      !> 1.2.1 trees
      !!
      CASE (itree)
        ! carbon increment
        leaf2wood_mass_ratio(:) = leaf2sapwood_mass_ratio(:) / (1._wp + lctlib_k_crtos)
        wood_carbon(:)          = sap_wood_carbon(:) + coarse_root_carbon(:)

        !! @todo calc_allocation_fraction() -> 1.2 Calculate allocation fractions -> 
        !! rewrite "CASE (biomass_alloc_scheme == "optimal")" such that instead of kopt1 and 2 leaf2root_mass ratio can be given
        SELECT CASE(biomass_alloc_scheme)
          CASE("optimal")
            WHERE(b_inc(:) > eps8)
              !! needed for stability - to be checked BLARPP! 
              !! @todo calc_allocation_fraction() -> 1.2 Calculate allocation fractions -> check stability issues with b_inc(:)
              b_inc(:) = 0.2_wp * dtime / one_day 

              !> 1.2.1a Calculate allocation fractions assuming optimal relationship between leaf and root allocation
              !! 
              leaf_carbon_inc(:) = MAX((b_inc(:) - k1_opt(:) + root_carbon(:) + wood_carbon(:) - dheight(:) / &
                                        leaf2wood_mass_ratio(:) * (b_inc(:)-k1_opt(:) + root_carbon(:)) - height(:) / &
                                        leaf2wood_mass_ratio(:)*leaf_carbon(:)) / &
                                       (1.0_wp + k2_opt(:) + dheight(:) / leaf2wood_mass_ratio(:) * &
                                        (b_inc(:) - (1.0_wp + k2_opt(:)) * leaf_carbon(:) - k1_opt(:) + root_carbon(:)) + &
                                        height(:)/leaf2wood_mass_ratio(:)), &
                                       0.0_wp)
              sap_wood_carbon_inc(:) = MAX((( leaf_carbon(:) + leaf_carbon_inc(:) ) / leaf2sapwood_mass_ratio(:) * &
                                            height(:) - sap_wood_carbon(:)) / & 
                                           (1._wp - ( leaf_carbon(:) + leaf_carbon_inc(:) ) / &
                                            leaf2sapwood_mass_ratio(:) * dheight(:) ), &
                                           0.0_wp)        
              WHERE(sap_wood_carbon_inc(:) > 0.0_wp)
                coarse_root_carbon_inc(:) = MAX((sap_wood_carbon(:) + sap_wood_carbon_inc(:)) * &
                                                lctlib_k_crtos - coarse_root_carbon(:), 0.0_wp)
                root_carbon_inc(:)        = MAX(k1_opt(:) + k2_opt(:) * leaf_carbon_inc(:) - root_carbon(:), 0.0_wp)
              ELSEWHERE
                leaf_carbon_inc(:)        = MAX((b_inc(:) - k1_opt(:) + root_carbon(:)) / (1.0_wp + k2_opt(:)), 0.0_wp)
                root_carbon_inc(:)        = MAX(k1_opt(:) + k2_opt(:) * leaf_carbon_inc(:) - root_carbon(:), 0.0_wp)
                sap_wood_carbon_inc(:)    = 0.0_wp
                coarse_root_carbon_inc(:) = 0.0_wp
              END WHERE
            END WHERE
          CASE("fixed", "dynamic")
            WHERE(b_inc(:) > eps8)
              !! needed for stability - to be checked BLARPP! 
              !! @todo calc_allocation_fraction() -> 1.2 Calculate allocation fractions -> check stability issues with b_inc(:)
              b_inc(:) = 0.2_wp * dtime / one_day 

              !> 1.2.1b solving the allometric rations between leaf and root, leaf and sapwood and sapwood and coarse roots 
              !! simulataneously with the height-to-diameter relationship 
              leaf_carbon_inc(:)     = MAX((b_inc(:) + leaf_carbon(:) * dheight(:) / leaf2wood_mass_ratio(:) * & 
                                            (leaf_carbon(:) / leaf2root_mass_ratio(:) - b_inc(:) - root_carbon(:) ) + &
                                            wood_carbon(:) + root_carbon(:) - leaf_carbon(:) * &
                                            (height(:) / leaf2wood_mass_ratio(:) + 1.0_wp / leaf2root_mass_ratio(:))) / & 
                                           (1.0_wp + ( b_inc(:) - (1.0_wp + 2.0_wp / leaf2root_mass_ratio(:)) * &
                                            leaf_carbon(:) + root_carbon(:)) * dheight(:) / &
                                            leaf2wood_mass_ratio(:) + height(:) / &
                                            leaf2wood_mass_ratio(:) + 1._wp / leaf2root_mass_ratio(:) ), & 
                                           0.0_wp) 

              sap_wood_carbon_inc(:) = MAX((( leaf_carbon(:) + leaf_carbon_inc(:) ) / leaf2sapwood_mass_ratio(:) * &
                                            height(:) - sap_wood_carbon(:)) / & 
                                           (1._wp - ( leaf_carbon(:) + leaf_carbon_inc(:) ) / leaf2sapwood_mass_ratio(:) * &
                                            dheight(:) ), &
                                           0.0_wp)          
              WHERE(sap_wood_carbon_inc(:) > 0.0_wp)
                coarse_root_carbon_inc(:) = MAX((sap_wood_carbon(:) + sap_wood_carbon_inc(:)) * &
                                                lctlib_k_crtos - coarse_root_carbon(:), 0.0_wp)
                root_carbon_inc(:)        = MAX((leaf_carbon(:) + leaf_carbon_inc(:)) / &
                                                leaf2root_mass_ratio(:) - root_carbon(:), 0.0_wp)
              ELSEWHERE ! special case of no wood growth (phenology etc.): redistribute available carbon 
                leaf_carbon_inc(:)        = MAX((b_inc(:) - leaf_carbon(:) / leaf2root_mass_ratio(:) + root_carbon(:)) / & 
                                                (1.0_wp + 1._wp / leaf2root_mass_ratio(:)), 0.0_wp)
                root_carbon_inc(:)        = MAX(0.0_wp, &
                                                MIN(b_inc(:) - leaf_carbon_inc(:), & 
                                                    (leaf_carbon_inc(:) + leaf_carbon(:)) / &
                                                    leaf2root_mass_ratio(:) - root_carbon(:)))
                sap_wood_carbon_inc(:)    = 0.0_wp
                coarse_root_carbon_inc(:) = 0.0_wp
              END WHERE
            END WHERE
        END SELECT ! biomass_alloc_scheme

      !> 1.2.2 grass
      !!
      CASE(igrass)
        SELECT CASE(biomass_alloc_scheme)
          CASE("optimal")
            leaf_carbon_inc(:)        = MAX((b_inc(:) - k1_opt(:) - leaf_carbon(:) / &
                                             leaf2sapwood_mass_ratio(:) + root_carbon(:) + sap_wood_carbon(:)) / &
                                            (1.0_wp + k2_opt(:) + 1.0_wp / leaf2sapwood_mass_ratio(:)), 0.0_wp)
            root_carbon_inc(:)        = MAX(k1_opt(:) + k2_opt(:) * leaf_carbon_inc(:) - root_carbon(:), 0.0_wp)
            sap_wood_carbon_inc(:)    = ( leaf_carbon_inc(:) + leaf_carbon(:) ) / leaf2sapwood_mass_ratio(:) - sap_wood_carbon(:)
            coarse_root_carbon_inc(:) = 0.0_wp
          CASE("fixed", "dynamic")
            leaf_carbon_inc(:)        = MAX((b_inc(:) - leaf_carbon(:) / leaf2root_mass_ratio(:) + &
                                             sap_wood_carbon(:) + root_carbon(:)) / & 
                                            (1.0_wp + 1._wp / leaf2sapwood_mass_ratio(:) + 1._wp / leaf2root_mass_ratio(:)), &
                                            0.0_wp) 
            root_carbon_inc(:)        = MAX(MIN(b_inc(:), ( leaf_carbon(:) + leaf_carbon_inc(:)) / &
                                                leaf2root_mass_ratio(:) - root_carbon(:)), &
                                            0.0_wp)
            sap_wood_carbon_inc(:)    = ( leaf_carbon_inc(:) + leaf_carbon(:) ) / leaf2sapwood_mass_ratio(:) - sap_wood_carbon(:)
            coarse_root_carbon_inc(:) = 0.0_wp
         END SELECT ! biomass_alloc_scheme

      !> 1.2.3 default - write and exit
      !!
      CASE DEFAULT
       !! Write statement
       !! 
       WRITE(message_text,'(a,i4)') 'invalid growthform type ',lctlib_growthform
       CALL finish("mo_veg_growth",message_text)
    END SELECT

    !> 1.3 Convert biomass increments to allocation fractions
    !! taking account of the previously calculated fraction of growth that goes to fruits
    !!
    WHERE(b_inc(:) > eps8)
       all_inc(:) = leaf_carbon_inc(:) + sap_wood_carbon_inc(:) + root_carbon_inc(:) + coarse_root_carbon_inc(:)
       WHERE(all_inc(:) > eps8)
          veg_falloc%leaf%carbon(:)        = (1._wp - veg_falloc%fruit%carbon(:) ) * leaf_carbon_inc(:)        / all_inc(:)
          veg_falloc%fine_root%carbon(:)   = (1._wp - veg_falloc%fruit%carbon(:) ) * root_carbon_inc(:)        / all_inc(:)
          veg_falloc%coarse_root%carbon(:) = (1._wp - veg_falloc%fruit%carbon(:) ) * coarse_root_carbon_inc(:) / all_inc(:)
          veg_falloc%sap_wood%carbon(:)    = (1._wp - veg_falloc%fruit%carbon(:) ) * sap_wood_carbon_inc(:)    / all_inc(:)
       END WHERE
    END WHERE


    !> 2.0 Allocation factors for N, derived from the target C:N ratios and the C partitioning coefficients
    !!  
    WHERE(b_inc(:) > eps8)
       !> 2.1 N required to sustain a unit of C growth
       !!
       growth_req_n(:) = veg_falloc%leaf%carbon(:)        / target_cn_leaf        + & 
                         veg_falloc%fine_root%carbon(:)   / target_cn_fine_root   + &
                         veg_falloc%coarse_root%carbon(:) / target_cn_coarse_root + &
                         veg_falloc%sap_wood%carbon(:)    / target_cn_sap_wood    + & 
                         veg_falloc%fruit%carbon(:)       / target_cn_fruit 
   
       !> 2.2 Adjustment of N allocation factors
       !!
       WHERE(growth_req_n(:) > eps8)
          veg_falloc%leaf%nitrogen(:)        = veg_falloc%leaf%carbon(:)        / target_cn_leaf        / growth_req_n(:)
          veg_falloc%fine_root%nitrogen(:)   = veg_falloc%fine_root%carbon(:)   / target_cn_fine_root   / growth_req_n(:) 
          veg_falloc%coarse_root%nitrogen(:) = veg_falloc%coarse_root%carbon(:) / target_cn_coarse_root / growth_req_n(:)
          veg_falloc%sap_wood%nitrogen(:)    = veg_falloc%sap_wood%carbon(:)    / target_cn_sap_wood    / growth_req_n(:)
          veg_falloc%fruit%nitrogen(:)       = veg_falloc%fruit%carbon(:)       / target_cn_fruit       / growth_req_n(:)
       END WHERE


       !> 3.0 Allocation factors for P, derived from the target N:P ratios and the N partitioning coefficients
       !! 
       !! 
       ! P required to sustain a unit of N growth
       growth_req_p(:) = veg_falloc%leaf%nitrogen(:)        / target_np_leaf        + & 
                         veg_falloc%fine_root%nitrogen(:)   / target_np_fine_root   + &
                         veg_falloc%coarse_root%nitrogen(:) / target_np_coarse_root + &
                         veg_falloc%sap_wood%nitrogen(:)    / target_np_sap_wood    + & 
                         veg_falloc%fruit%nitrogen(:)       / target_np_fruit 
   
       ! Adjustment of P allocation factors
       WHERE(growth_req_p(:) > eps8)
          veg_falloc%leaf%phosphorus(:)        = veg_falloc%leaf%nitrogen(:)        / target_np_leaf        / growth_req_p(:)
          veg_falloc%fine_root%phosphorus(:)   = veg_falloc%fine_root%nitrogen(:)   / target_np_fine_root   / growth_req_p(:) 
          veg_falloc%coarse_root%phosphorus(:) = veg_falloc%coarse_root%nitrogen(:) / target_np_coarse_root / growth_req_p(:)
          veg_falloc%sap_wood%phosphorus(:)    = veg_falloc%sap_wood%nitrogen(:)    / target_np_sap_wood    / growth_req_p(:)
          veg_falloc%fruit%phosphorus(:)       = veg_falloc%fruit%nitrogen(:)       / target_np_fruit       / growth_req_p(:)
       END WHERE

   
       !> 4.0 Allocation factors for isotopes follow the main isotope 
       !! 
       !! 
       veg_falloc%leaf%carbon13(:)          = veg_falloc%leaf%carbon(:)
       veg_falloc%fine_root%carbon13(:)     = veg_falloc%fine_root%carbon(:)
       veg_falloc%coarse_root%carbon13(:)   = veg_falloc%coarse_root%carbon(:)
       veg_falloc%sap_wood%carbon13(:)      = veg_falloc%sap_wood%carbon(:)
       veg_falloc%fruit%carbon13(:)         = veg_falloc%fruit%carbon(:)
   
       veg_falloc%leaf%carbon14(:)          = veg_falloc%leaf%carbon(:)
       veg_falloc%fine_root%carbon14(:)     = veg_falloc%fine_root%carbon(:)
       veg_falloc%coarse_root%carbon14(:)   = veg_falloc%coarse_root%carbon(:)
       veg_falloc%sap_wood%carbon14(:)      = veg_falloc%sap_wood%carbon(:)
       veg_falloc%fruit%carbon14(:)         = veg_falloc%fruit%carbon(:)
   
       veg_falloc%leaf%nitrogen15(:)        = veg_falloc%leaf%nitrogen(:)
       veg_falloc%fine_root%nitrogen15(:)   = veg_falloc%fine_root%nitrogen(:)
       veg_falloc%coarse_root%nitrogen15(:) = veg_falloc%coarse_root%nitrogen(:)
       veg_falloc%sap_wood%nitrogen15(:)    = veg_falloc%sap_wood%nitrogen(:)
       veg_falloc%fruit%nitrogen15(:)       = veg_falloc%fruit%nitrogen(:)
    END WHERE ! b_inc(:) > eps8

  END SUBROUTINE calc_allocation_fraction


  !-----------------------------------------------------------------------------------------------------
  ! 6.0 Sub Task of update_veg_growth
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> Calculate partitioning of labile pool outflow into maintenance, growth respiration and growth 
  !! Follows the concept of Amthor: Growth = Yg * (P - Rm), where\n
  !!   Yg = growth efficiency of biomass production (1-fresp_growth)\n
  !!   P  = current turnover of labile pool (Amthor's production)\n
  !!   Rm = maintenance respiration
  !!
  !! In the current implementation, maintenance is given preference over growth, but 
  !! nevertheless, maintenance is limited to the available C   
  !! 
  !! Note on veg_growth: it was an explicit decision making all veg_growth\%POOLs INTENT(inout) 
  !!       but overwrite veg_growth%labile (making it essentially INTENT(in)), see Ticket #162 and r504
  !! 
  !! Input: labile pool size and current turnover rate, N demand for respiration, & 
  !!        N & P requirements for growth, allocation fractions
  !! 
  !! output: growth of tissue types, maintenance, growth and N processing respiration, exudation 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_partitioning( nc                            , &
                                kstar_labile                  , &
                                maint_respiration_pot         , &
                                maint_respiration_leaf        , &
                                growth_req_p                  , &
                                growth_req_n                  , &
                                veg_falloc                    , &
                                veg_pool_labile               , &
                                n_processing_respiration      , &
                                maint_respiration_c13         , &
                                maint_respiration_c14         , &
                                growth_respiration_c13        , &
                                growth_respiration_c14        , &
                                n_processing_respiration_c13  , &
                                n_processing_respiration_c14  , &
                                veg_growth                    , &
                                veg_exudation                 , &
                                maint_respiration             , &
                                growth_respiration  )  

    USE mo_kind,                  ONLY: wp
    USE mo_jsb_math_constants,    ONLY: one_day, dtime, eps8
    USE mo_veg_constants,         ONLY: knut_labile, tau_labile, fresp_growth
#ifdef __QUINCY_STANDALONE__
    dsl4slm_Use_memory_a_veg
#else
    dsl4jsb_Use_pools(VEG_) t_veg_tissues_pool, t_veg_elements_pool
#endif

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                          INTENT(in)    :: nc                               !< dimensions
    REAL(wp), DIMENSION(nc),          INTENT(in)    :: kstar_labile                 , & !< potential labile pool turnover (1/day) 
                                                       maint_respiration_pot        , & !< potential (demand-based) plant maintenance respiration
                                                       maint_respiration_leaf       , & !< foliar maintenance respiration
                                                       growth_req_p                 , & !< moles P required for a unit of N growth under current allocation fractions
                                                       growth_req_n                     !< moles N required for a unit of C growth under current allocation fractions
#ifdef __QUINCY_STANDALONE__
    TYPE(t_veg_pool_a),               INTENT(in)    :: veg_falloc                       !< fractional allocation to pools
    TYPE(t_veg_element_a),            INTENT(in)    :: veg_pool_labile                  !< labile pool size of the plant 
#else
    TYPE(t_veg_tissues_pool),         INTENT(in)    :: veg_falloc                       !< fractional allocation to pools
    TYPE(t_veg_elements_pool),        INTENT(in)    :: veg_pool_labile                  !< labile pool size of the plant 
#endif
    REAL(wp), DIMENSION(nc),          INTENT(inout) :: n_processing_respiration     , & !< respiration associated with N processing
                                                       maint_respiration_c13        , & !< maintenance respiration
                                                       maint_respiration_c14        , & !< maintenance respiration
                                                       growth_respiration_c13       , & !< growth respiration
                                                       growth_respiration_c14       , & !< growth respiration
                                                       n_processing_respiration_c13 , & !< respiration associated with N processing
                                                       n_processing_respiration_c14     !< respiration associated with N processing
#ifdef __QUINCY_STANDALONE__
    TYPE(t_veg_pool_a),               INTENT(inout) :: veg_growth                       !< fluxes: growth; see section 4.1: veg_growth%labile is overwriten
    TYPE(t_veg_element_a),            INTENT(inout) :: veg_exudation                    !< vegetation flux to mycorrhiza (mol / m2 / timestep) 
#else
    TYPE(t_veg_tissues_pool),         INTENT(inout) :: veg_growth                       !< fluxes: growth; see section 4.1: veg_growth%labile is overwriten
    TYPE(t_veg_elements_pool),        INTENT(inout) :: veg_exudation                    !< vegetation flux to mycorrhiza (mol / m2 / timestep) 
#endif
    REAL(wp), DIMENSION(nc),          INTENT(out)   :: maint_respiration            , & !< maintenance respiration (micro-mol C / m2 / s)
                                                       growth_respiration               !< growth respiration (micro-mol C / m2 / s)
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc)                         :: hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_partitioning' 


    !> 1.0 Determine C use for maintenance (micro-mol C / m2 / s)
    !! 
    !! based on tissue N and Temperature, constrained by default turnover of the labile pool
    !! to avoid that maitenance respiration consumes all C in times of low GPP 
    !! constrained to support leaf respiration, but not more than the current labile pool turnover
    maint_respiration(:) = MAX(maint_respiration_leaf(:), & 
                               MIN(1.0_wp / tau_labile / one_day * veg_pool_labile%carbon(:) * 1.e6_wp, &
                                   maint_respiration_pot(:))) 


    !> 2.0 Determine C, N, P use for growth
    !!
    !> 2.1 potential C growth rate (mol C / m2 / timestep) 
    !! is determined from the labile turnover minus the carbon used for maintenance, i
    !! accounting for the fraction of growth respiration required to achieve this growth, to ensure that
    !! growth does not empty the labile pool
    veg_growth%labile%carbon(:) = (veg_pool_labile%carbon(:) * kstar_labile(:) / one_day - &
                                   (maint_respiration(:) + n_processing_respiration(:)) / 1.e6_wp ) * & 
                                  dtime /  (1._wp + fresp_growth)

    WHERE(veg_growth%labile%carbon(:) < eps8) 
      veg_growth%labile%carbon(:) = 0.0_wp
    END WHERE

    !> 2.2 Elemental limitation of growth by N & P (mol [N|P] / m2 / timestep)
    !! is determined by the minimum turnover of the labile N and P pools and the stoichiometric growth requirement
    !! given the potential growth rate of C
    !! The use of knut_labile allows plants to draw down the C:N:P of the labile pool, which will cause 
    !! a) an increased uptake per unit root mass and b) an increased cause 
    veg_growth%labile%nitrogen(:)   = MIN(growth_req_n(:) * veg_growth%labile%carbon(:), & 
                                          veg_pool_labile%nitrogen(:) * knut_labile * kstar_labile(:) * &
                                          dtime / one_day) 
    veg_growth%labile%phosphorus(:) = MIN(growth_req_p(:) * veg_growth%labile%nitrogen(:), & 
                                          veg_pool_labile%phosphorus(:) * knut_labile * knut_labile * kstar_labile(:) * &
                                          dtime / one_day) 
    ! adjust N growth rate to P growth rate if P is limiting 
    WHERE(growth_req_p(:) > eps8)
      hlp1(:) = veg_growth%labile%phosphorus(:) / growth_req_p(:)
      WHERE(veg_growth%labile%nitrogen(:) > hlp1(:))
        veg_growth%labile%nitrogen(:) = veg_growth%labile%phosphorus(:) / growth_req_p(:)
      END WHERE
    END WHERE
    ! adjust C growth rate to N growth rate if N is limiting 
    WHERE(growth_req_n(:) > eps8)
      hlp1(:) = veg_growth%labile%nitrogen(:) / growth_req_n(:)
      WHERE(veg_growth%labile%carbon(:) > hlp1(:))
        veg_growth%labile%carbon(:) = veg_growth%labile%nitrogen(:) / growth_req_n(:)
      END WHERE
    END WHERE  
    ! Diagnose growth respiration (micro-mol C / m2 / s)
    growth_respiration(:)   = fresp_growth * veg_growth%labile%carbon(:) * 1.e6_wp / dtime


    !> 3.0 Adjust exudation rate to available carbon (mol C / m2 / timestep)
    veg_exudation%carbon(:) = MIN(veg_exudation%carbon(:), &
                                  MAX(0.0_wp, & 
                                      0.2_wp * veg_pool_labile%carbon(:) - veg_growth%labile%carbon(:) - &
                                      (growth_respiration(:) + maint_respiration(:) + n_processing_respiration(:)) * &
                                      dtime / 1.e6_wp))


    !> 4.0 Calculate isotopic signature of respiration fluxes, and growth
    !!
    WHERE(veg_pool_labile%carbon(:) > eps8)
       ! growth (in mol C / m2 / timestep)
       veg_growth%labile%carbon13(:)   = veg_growth%labile%carbon(:)  * veg_pool_labile%carbon13(:) / veg_pool_labile%carbon(:)
       veg_growth%labile%carbon14(:)   = veg_growth%labile%carbon(:)  * veg_pool_labile%carbon14(:) / veg_pool_labile%carbon(:)
       ! exudation (in mol C / m2 / timestep)
       veg_exudation%carbon13(:)       = veg_exudation%carbon(:)      * veg_pool_labile%carbon13(:) / veg_pool_labile%carbon(:)
       veg_exudation%carbon14(:)       = veg_exudation%carbon(:)      * veg_pool_labile%carbon14(:) / veg_pool_labile%carbon(:)
       ! respiration (in micro-mol C / m2 / s)
       maint_respiration_c13(:)        = maint_respiration(:)         * veg_pool_labile%carbon13(:) / veg_pool_labile%carbon(:)
       maint_respiration_c14(:)        = maint_respiration(:)         * veg_pool_labile%carbon14(:) / veg_pool_labile%carbon(:)
       growth_respiration_c13(:)       = growth_respiration(:)        * veg_pool_labile%carbon13(:) / veg_pool_labile%carbon(:)
       growth_respiration_c14(:)       = growth_respiration(:)        * veg_pool_labile%carbon14(:) / veg_pool_labile%carbon(:)
       n_processing_respiration_c13(:) = n_processing_respiration(:)  * veg_pool_labile%carbon13(:) / veg_pool_labile%carbon(:) 
       n_processing_respiration_c14(:) = n_processing_respiration(:)  * veg_pool_labile%carbon14(:) / veg_pool_labile%carbon(:) 
    END WHERE

    WHERE(veg_pool_labile%nitrogen(:) > eps8)
      ! growth (in mol  / m2 / timestep)
      veg_growth%labile%nitrogen15(:) = veg_growth%labile%nitrogen(:) * veg_pool_labile%nitrogen15(:) / veg_pool_labile%nitrogen(:)
    END WHERE


    !> 5.0 partition growth to tissue types
    !! 
    veg_growth%leaf        = veg_growth%leaf        + veg_falloc%leaf        * veg_growth%labile
    veg_growth%fine_root   = veg_growth%fine_root   + veg_falloc%fine_root   * veg_growth%labile
    veg_growth%coarse_root = veg_growth%coarse_root + veg_falloc%coarse_root * veg_growth%labile
    veg_growth%sap_wood    = veg_growth%sap_wood    + veg_falloc%sap_wood    * veg_growth%labile
    veg_growth%fruit       = veg_growth%fruit       + veg_falloc%fruit       * veg_growth%labile

    !> 5.1 invert flux of labile to be out of the labile pool
    !! 
    !! essentially, this calculation makes veg_growth%labile INTENT(in)
    veg_growth%labile%carbon(:)     = veg_growth%labile%carbon(:)     * (-1.0_wp)
    veg_growth%labile%nitrogen(:)   = veg_growth%labile%nitrogen(:)   * (-1.0_wp)
    veg_growth%labile%phosphorus(:) = veg_growth%labile%phosphorus(:) * (-1.0_wp)
    veg_growth%labile%carbon13(:)   = veg_growth%labile%carbon13(:)   * (-1.0_wp)
    veg_growth%labile%carbon14(:)   = veg_growth%labile%carbon14(:)   * (-1.0_wp)
    veg_growth%labile%nitrogen15(:) = veg_growth%labile%nitrogen15(:) * (-1.0_wp)

  END SUBROUTINE calc_partitioning 


  !-----------------------------------------------------------------------------------------------------
  ! Additonal function used in:
  !   calc_allometry
  !   update_veg_pools
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Function to calculate diameter of a tree given its biomass
  !! 
  !! Diameter and Height calculated assuming that: \n
  !!   all woody biomass forms a cylinder: diameter = SQRT(4 * wood volume / height / pi) \n
  !!   height = allom_k1 * diameter ^ allom_k2
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_diameter_from_woody_biomass( allom_k1          , &  ! lctlib
                                                            allom_k2          , &
                                                            wood_density      , &  ! lctlib
                                                            sap_wood_carbon   , &
                                                            heart_wood_carbon , &
                                                            dens_ind)           &
                                                            RESULT(diameter)

    USE mo_jsb_math_constants,      ONLY: pi

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),               INTENT(in)  :: allom_k1         , &    !< land-cover-type library parameter
                                           allom_k2         , &    !< land-cover-type library parameter
                                           wood_density            !< land-cover-type library parameter
    REAL(wp),               INTENT(in)  :: sap_wood_carbon         !< woody biomass of an individual tree (mol C/m2)
    REAL(wp),               INTENT(in)  :: heart_wood_carbon       !< woody biomass of an individual tree (mol C/m2)
    REAL(wp),               INTENT(in)  :: dens_ind                !< density of populations (#indivuduals /m2) 
    REAL(wp)                            :: diameter                !< diameter of an individual (m2)
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_diameter_from_woody_biomass'


    diameter = ( ( sap_wood_carbon + heart_wood_carbon ) / dens_ind / & 
               (wood_density * pi/4._wp * allom_k1)) ** (1._wp/(2._wp + allom_k2))

  END FUNCTION calc_diameter_from_woody_biomass


  !-----------------------------------------------------------------------------------------------------
  ! Additonal function used in:
  !   calc_allometry
  !   update_veg_pools
  ! 
  !----------------------------------------------------------------------------------------------------------
  !> Function to calculate height of a tree given its diameter 
  !!
  !! height = allom_k1 * diameter ^ allom_k2
  !! min_height < height
  !----------------------------------------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION calc_height_from_diameter(allom_k1, allom_k2, diameter) RESULT(height)

    USE mo_veg_constants,           ONLY: min_height

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),               INTENT(in)  :: allom_k1 , & !< land-cover-type library parameter
                                           allom_k2 , & !< land-cover-type library parameter
                                           diameter     !< diameter of an individual (m2) 
    REAL(wp)                            :: height       !< height of an individual (m)
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_height_from_diameter'


    height = MAX(min_height,allom_k1 * (diameter ** allom_k2))

  END FUNCTION calc_height_from_diameter

#endif
END MODULE mo_veg_growth
