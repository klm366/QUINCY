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
!! @brief Contains the routines for the Task sb_simple_model 
!!
!! @author
!!  Silvia Caldararu, Soenke Zaehle 
!!
!! @par Revision History
!! First version by Soenke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!! Second version by Soenke Zaehle (20170331): added a 4 pool model 
!! 3rd version: Jan Engel (20170627): split calc_sb_fluxes_simple() into multiple calc_ routines and rename to update_sb_simple_model()
!!
MODULE mo_sb_simple_model

  USE mo_kind,                    ONLY: wp
  USE mo_jsb_control,             ONLY: debug_on
  USE mo_exception,               ONLY: message, message_text, finish
  
  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, SB_BGCM_TRANSPORT_ID, &
    &                                   SB_BGCM_MYCO_EXPORT_ID
  
  USE mo_jsb_math_constants,      ONLY: dtime, one_day, one_year, eps8, eps4

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_simple_model
  PUBLIC :: calc_sb_inorganic_n15_fluxes  ! used by mo_sb_jsm_main

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_simple_model'

CONTAINS

  ! ======================================================================================================= !
  !>update soil biogeochemical pools 
  !>
  !>  routine to update soil litter and mineral pools (simple version for vegetation model testing)
  !>
  SUBROUTINE update_sb_simple_model(tile, options)
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_jsb_process_class,     ONLY: SB_, VEG_, SPQ_, A2L_
    USE mo_jsb_grid_class,        ONLY: t_jsb_vgrid
    USE mo_jsb_grid,              ONLY: Get_vgrid
    USE mo_veg_constants,         ONLY: ITREE
    USE mo_sb_jsm_transport,      ONLY: calc_bioturbation_rate
    USE mo_sb_litter_processes,   ONLY: calc_litter_partitioning
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(SPQ_)
    dsl4jsb_Use_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model                !< the model
    TYPE(t_lnd_bgcm_store), POINTER       :: bgcm_store           !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_soil_sb        !< Vertical grid
    INTEGER                               :: nsoil_sb             !< number of soil layers as used/defined by the SB_ process
    INTEGER                               :: isoil                !< loop over nsoil_sb
    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: apparent_fast_som_cn      , &
                                             apparent_fast_som_np      , &
                                             fine_root_carbon_sl, &
                                             km1_uptake_nh4_act, &
                                             km2_uptake_nh4_act, &
                                             km1_uptake_no3_act, &
                                             km2_uptake_no3_act, &
                                             km1_uptake_po4_act, &
                                             km2_uptake_po4_act, &
                                             km_uptake_org_act, &
                                             km_mic_uptake_nh4_act, &
                                             km_mic_uptake_no3_act, &
                                             km_denitrification_c_act, &
                                             km_denitrification_no3_act, &
                                             p_deposition_act, &
                                             vmax_weath_mineral_act_sl, &
                                             k_adsorpt_po4_act, &
                                             k_desorpt_po4_act, &
                                             km_dip_biochem_po4_act, &
                                             qmax_po4_mineral_act
    INTEGER                               :: iblk, ics, ice, nc             !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_simple_model'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt2L3D :: sb_pool_mt
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt2L3D :: sb_transport_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SPQ_)
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk :: unit_uptake_p_pot
    ! VEG_ 3D
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    ! A2L_ 2D
    dsl4jsb_Real2D_onChunk :: p_deposition
    ! SPQ_ 3D
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl 
    dsl4jsb_Real3D_onChunk :: percolation_sl
    dsl4jsb_Real3D_onChunk :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onChunk :: w_soil_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: rtm_decomposition 
    dsl4jsb_Real3D_onChunk :: rtm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rtm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: rmm_decomposition 
    dsl4jsb_Real3D_onChunk :: rmm_plant_uptake
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_asymb_bnf
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl 
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: myc_respiration
    dsl4jsb_Real3D_onChunk :: myc_respiration_c13
    dsl4jsb_Real3D_onChunk :: myc_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nitrification_no3 
    dsl4jsb_Real3D_onChunk :: nitrification_noy 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_noy 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o 
    dsl4jsb_Real3D_onChunk :: denitrification_n2 
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_exchange_po4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow 
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2 
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2 
    dsl4jsb_Real3D_onChunk :: nh4_assoc 
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: k_bioturb 
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation 
    dsl4jsb_Real3D_onChunk :: qmax_po4
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2_n15 
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_solute_sl 
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_solute_sl 
    dsl4jsb_Real3D_onChunk :: lateral_loss_po4_solute_sl 
    dsl4jsb_Real3D_onChunk :: lateral_loss_nh4_n15_solute_sl 
    dsl4jsb_Real3D_onChunk :: lateral_loss_no3_n15_solute_sl 
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2_n15 
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15 
    dsl4jsb_Real3D_onChunk :: vmax_weath_mineral_sl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_soil_sb => Get_vgrid('soil_layer_sb')
    nsoil_sb      =  vgrid_soil_sb%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SPQ_)
    dsl4jsb_Get_memory(A2L_) 
    ! ----------------------------------------------------------------------------------------------------- !
    ALLOCATE(apparent_fast_som_cn(nc, nsoil_sb))
    ALLOCATE(apparent_fast_som_np(nc, nsoil_sb))
    ALLOCATE(fine_root_carbon_sl(nc, nsoil_sb))
    ALLOCATE(km1_uptake_nh4_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_nh4_act(nc, nsoil_sb))
    ALLOCATE(km1_uptake_no3_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_no3_act(nc, nsoil_sb))
    ALLOCATE(km1_uptake_po4_act(nc, nsoil_sb))
    ALLOCATE(km2_uptake_po4_act(nc, nsoil_sb))
    ALLOCATE(km_uptake_org_act(nc, nsoil_sb))
    ALLOCATE(km_mic_uptake_nh4_act(nc, nsoil_sb))
    ALLOCATE(km_mic_uptake_no3_act(nc, nsoil_sb))
    ALLOCATE(km_denitrification_c_act(nc, nsoil_sb))
    ALLOCATE(km_denitrification_no3_act(nc, nsoil_sb))
    ALLOCATE(p_deposition_act(nc, nsoil_sb))
    ALLOCATE(vmax_weath_mineral_act_sl(nc, nsoil_sb))
    ALLOCATE(k_adsorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(k_desorpt_po4_act(nc, nsoil_sb))
    ALLOCATE(km_dip_biochem_po4_act(nc, nsoil_sb))
    ALLOCATE(qmax_po4_mineral_act(nc, nsoil_sb))
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_POOL_ID, sb_pool_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_LOSS_ID, sb_loss_mt)
    dsl4jsb_Get_mt2L3D(SB_BGCM_TRANSPORT_ID, sb_transport_mt)
    dsl4jsb_Get_mt1L3D(SB_BGCM_MYCO_EXPORT_ID, sb_mycorrhiza_export_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Get_var2D_onChunk(A2L_,      p_deposition)                ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_n_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_p_pot)
    ! VEG_ 3D
    dsl4jsb_Get_var3D_onChunk(VEG_,      root_fraction_sl)
    ! SPQ_ 3D
    dsl4jsb_Get_var3D_onChunk(SPQ_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SPQ_,     bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SPQ_,     percolation_sl)
    dsl4jsb_Get_var3D_onChunk(SPQ_,     frac_w_lat_loss_sl) 
    dsl4jsb_Get_var3D_onChunk(SPQ_,     w_soil_sl)                   ! in
    ! SB_ 2D
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)          
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy_n15)          
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o_n15)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2_n15)         
    ! SB_ 3D
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_plant_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_asymb_bnf)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_plant_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_asymb_bnf)
    dsl4jsb_Get_var3D_onChunk(SB_,       anaerobic_volume_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       myc_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_exchange_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_nh4_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_nh4_n15_solute_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_no3_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_po4_solute_sl) 
    dsl4jsb_Get_var3D_onChunk(SB_,       lateral_loss_no3_n15_solute_sl)  
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       noy_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation_n15)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       vmax_weath_mineral_sl)
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 provide soil profile values
    !>  calculate local variables
    !>  
    CALL calc_local_soil_profile_values( &
      & nc, &                                   ! in
      & nsoil_sb, &
      & soil_depth_sl(:,:), &
      & bulk_dens_sl(:,:), &
      & w_soil_sl(:,:), &
      & nh4_solute(:,:), &
      & rtm_hsc(:,:), &
      & rmm_hsc(:,:), &
      & rtm_sorption(:,:), &
      & rtm_desorption(:,:), &
      & rmm_sorption(:,:), &
      & rmm_desorption(:,:), &
      & bulk_dens_corr_sl(:,:), &
      & p_deposition(:), &
      & vmax_weath_mineral_sl(:,:), &           ! in
      & km1_uptake_nh4_act(:,:), &              ! out
      & km2_uptake_nh4_act(:,:), &
      & km1_uptake_no3_act(:,:), &
      & km2_uptake_no3_act(:,:), &
      & km1_uptake_po4_act(:,:), &
      & km2_uptake_po4_act(:,:), &
      & km_uptake_org_act(:,:), &
      & km_mic_uptake_nh4_act, &
      & km_mic_uptake_no3_act, &
      & km_denitrification_c_act(:,:), &
      & km_denitrification_no3_act(:,:), &
      & apparent_fast_som_cn(:,:), &
      & apparent_fast_som_np(:,:), &
      & p_deposition_act(:,:), &
      & vmax_weath_mineral_act_sl(:,:), &
      & k_adsorpt_po4_act(:,:), &
      & k_desorpt_po4_act(:,:), &
      & km_dip_biochem_po4_act(:,:) )           ! out

    !>2.0 Calculate potential rates of uptake and losses
    !>
    
    !>  2.1 Potential plant nutrient uptake rates
    !>
    CALL calc_potential_plant_uptake_rate( &
      & nc, &                                     ! in
      & nsoil_sb, &
      & lctlib%vmax_uptake_n, &
      & lctlib%vmax_uptake_p, &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza, &
      & soil_depth_sl(:,:), &
      & veg_pool_mt(ix_fine_root, ixC, :), &
      & sb_pool_mt(ix_mycorrhiza, ixC, :, :), &
      & root_fraction_sl(:,:), &
      & rtm_plant_uptake(:,:), &
      & rmm_plant_uptake(:,:), &
      & km1_uptake_nh4_act(:,:), &
      & km2_uptake_nh4_act(:,:), &
      & km1_uptake_no3_act(:,:), &
      & km2_uptake_no3_act(:,:), &
      & km1_uptake_po4_act(:,:), &
      & km2_uptake_po4_act(:,:), &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & po4_solute(:,:), &                        ! in
      & plant_uptake_nh4_sl(:,:), &               ! out
      & plant_uptake_no3_sl(:,:), &
      & plant_uptake_po4_sl(:,:), &
      & unit_uptake_n_pot(:), &
      & unit_uptake_p_pot(:) )                    ! out

    !>  2.2 Potential mycorrhiza nutrient uptake rates
    !>
    IF (dsl4jsb_Config(SB_)%flag_mycorrhiza) THEN
      CALL calc_potential_mycorrhiza_uptake_rate( &
        & nc, &                                     ! in
        & nsoil_sb, &
        & lctlib%vmax_uptake_n, &
        & lctlib%vmax_uptake_p, &
        & soil_depth_sl(:,:), &
        & veg_pool_mt(ix_fine_root, ixC, :), &
        & sb_pool_mt(ix_mycorrhiza, ixC, :, :), &
        & sb_pool_mt(ix_mycorrhiza, ixN, :, :), &
        & rtm_plant_uptake(:,:), &
        & rmm_plant_uptake(:,:), &
        & km1_uptake_nh4_act(:,:), &
        & km2_uptake_nh4_act(:,:), &
        & km1_uptake_no3_act(:,:), &
        & km2_uptake_no3_act(:,:), &
        & km1_uptake_po4_act(:,:), &
        & km2_uptake_po4_act(:,:), &
        & nh4_solute(:,:), &
        & no3_solute(:,:), &
        & po4_solute(:,:), &                        ! in
        & mycorrhiza_uptake_nh4_sl(:,:), &          ! out
        & mycorrhiza_uptake_no3_sl(:,:), &
        & mycorrhiza_uptake_po4_sl(:,:) )           ! out
    END IF

    !>  2.3 Calculate the potential decomposition rate of soil organic matter and litter pools 
    !>
    CALL calc_potential_decomposition_rate( &
      & nc, &                                       ! in
      & nsoil_sb, &
      & soil_depth_sl(:,:), &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza_prim, &   
      & dsl4jsb_Config(SB_)%flag_mycorrhiza_org, &                    
      & rtm_decomposition(:,:), &
      & rmm_decomposition(:,:), &
      & sb_pool_mt(:,:,:,:), &                      ! in
      & het_respiration(:,:), &                     ! inout
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & myc_respiration(:,:), &
      & myc_respiration_c13(:,:), &
      & myc_respiration_c14(:,:), &
      & sb_loss_mt(:,:,:,:) )                       ! inout

    !>  2.4 Loss pathways for the dynamic loss calculation scheme 
    !>
    IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme) == 'dynamic')THEN
      
      CALL calc_potential_nitrification_denitrification_rate( &
        & nc, &                                     ! in
        & nsoil_sb, &
        & rtm_nitrification(:,:), &
        & rmm_nitrification(:,:), &
        & rtm_denitrification(:,:), &
        & km_denitrification_c_act(:,:), &
        & km_denitrification_no3_act(:,:), &
        & anaerobic_volume_fraction_sl(:,:), &
        & nh4_solute(:,:), &
        & no3_solute(:,:), &
        & sb_pool_mt(ix_microbial, ixC, :, :), &    ! in
        & nitrification_no3(:,:), &                 ! out
        & denitrification_n2(:,:) )                 ! out
        
      CALL calc_inorganic_transport_rate( &
        & nc, &                                   ! in
        & nsoil_sb, &
        & soil_depth_sl(:,:), &
        & percolation_sl(:,:), &
        & frac_w_lat_loss_sl(:,:), & 
        & rtm_gasdiffusion(:,:), &
        & rmm_gasdiffusion(:,:), &
        & nh4_solute(:,:), &
        & nh4_n15_solute(:,:), &
        & no3_solute(:,:), &
        & no3_n15_solute(:,:), &
        & po4_solute(:,:), &
        & noy(:,:), &
        & noy_n15(:,:), &
        & n2o(:,:), &
        & n2o_n15(:,:), &
        & n2(:,:), &
        & n2_n15(:,:), &                          ! in
        & transport_nh4_solute(:,:), &            ! inout
        & transport_nh4_n15_solute(:,:), &
        & transport_no3_solute(:,:), &
        & transport_no3_n15_solute(:,:), &
        & transport_po4_solute(:,:), &            ! inout
        & transport_noy(:,:), &                   ! out
        & transport_noy_n15(:,:), &   
        & transport_n2o(:,:), &     
        & transport_n2o_n15(:,:), & 
        & transport_n2(:,:), &      
        & transport_n2_n15(:,:), & 
        & leaching_nh4_solute(:), &
        & leaching_nh4_n15_solute(:), &
        & leaching_no3_solute(:), &
        & leaching_no3_n15_solute(:), &
        & leaching_po4_solute(:), &
        & lateral_loss_nh4_solute_sl(:,:), & 
        & lateral_loss_nh4_n15_solute_sl(:,:), & 
        & lateral_loss_no3_solute_sl(:,:), & 
        & lateral_loss_no3_n15_solute_sl(:,:), & 
        & lateral_loss_po4_solute_sl(:,:), & 
        & emission_noy(:), &    
        & emission_noy_n15(:), &
        & emission_n2o(:), &    
        & emission_n2o_n15(:), & 
        & emission_n2(:), &  
        & emission_n2_n15(:) )                    ! out
    
    ENDIF  ! IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme)=='dynamic')

    !>  2.5 transport due to bioturbation
    !>
    ! calculate the bioturbation rate per soil layer (vertical transport rate [1/s] )
    !   from bulk density and root fraction and k_diff_org
    k_bioturb(:,:) = calc_bioturbation_rate(nc, nsoil_sb, soil_depth_sl(:,:), bulk_dens_corr_sl(:,:), root_fraction_sl(:,:))
    ! bioturbation transport for the 'simple_1d' soil model
    CALL calc_bioturbation_transport_wrapper_s1d( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & soil_depth_sl(:,:), &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & k_bioturb(:,:), &
      & nh4_assoc, &
      & nh4_n15_assoc, & 
      & po4_assoc_fast, &
      & po4_assoc_slow, &
      & po4_occluded, &
      & sb_pool_mt(:,:,:,:), &                ! in
      & transport_nh4_assoc(:,:), &           ! inout
      & transport_nh4_n15_assoc(:,:), &
      & transport_po4_assoc_fast(:,:), &
      & transport_po4_assoc_slow(:,:), & 
      & transport_po4_occluded(:,:), &
      & sb_transport_mt(:,:,:,:) )            ! inout

    !>  2.6 Calculate the potential PO4 flux rate: biomineralisation, weathering, slow sorption, occlusion 
    !>
    ! calculate the fine root C of each soil layer, fine_root_carbon_sl [mol C m-3] 
    DO isoil = 1,nsoil_sb
      WHERE(soil_depth_sl(:,isoil) > eps8)
        fine_root_carbon_sl(:,isoil) = root_fraction_sl(:,isoil) / soil_depth_sl(:,isoil) * veg_pool_mt(ix_fine_root, ixC, :)
      ELSEWHERE
        fine_root_carbon_sl(:,isoil) = 0.0_wp
      ENDWHERE
    ENDDO
    ! calculate the potential PO4 flux rate
    CALL calc_potential_p_flux_rate( &
      & nc, &                                         ! in
      & nsoil_sb, &                  
      & dsl4jsb_Config(SB_)%flag_sb_prescribe_po4, &
      & rtm_decomposition(:,:), &
      & rmm_decomposition(:,:), &
      & km_dip_biochem_po4_act(:,:), &
      & bulk_dens_corr_sl(:,:), & 
      & fine_root_carbon_sl(:,:), &
      & po4_solute(:,:), &
      & po4_assoc_fast(:,:), &
      & po4_assoc_slow(:,:), &
      & po4_occluded(:,:), &
      & po4_primary(:,:), &
      & vmax_weath_mineral_act_sl(:,:), &
      & k_adsorpt_po4_act(:,:), &
      & k_desorpt_po4_act(:,:), &
      & sb_pool_mt(ix_residue, ixP, :, :), &          ! in
      & weathering_po4(:,:), &                        ! out
      & occlusion_po4(:,:), &
      & slow_exchange_po4(:,:), &
      & biochem_mineralisation_po4(:,:) )             ! out

    !>3.0 Actual fluxes of in and out of SOM 
    !>

    !>  3.1 Calculation of asymbiotic N fixation rate
    !>
    asymb_n_fixation(:,:) = calc_asymb_n_fixation( &
      &                         TRIM(dsl4jsb_Config(SB_)%bnf_scheme), & 
      &                         sb_loss_mt(ix_soluable_litter, ixC, :, :) + sb_loss_mt(ix_polymeric_litter, ixC, :, :), &
      &                         sb_loss_mt(ix_soluable_litter, ixN, :, :) + sb_loss_mt(ix_polymeric_litter, ixN, :, :), &
      &                         apparent_fast_som_cn(:,:), &
      &                         rtm_asymb_bnf(:,:), &
      &                         rmm_asymb_bnf(:,:), &
      &                         nh4_solute(:,:), &
      &                         no3_solute(:,:) )               

    !>  3.2 Calculation of actual plant, mycorrhizal, microbial uptake rate as well as 
    !>      nitrification, denitrification and biomineralisation rates
    !>
    CALL calc_actual_flux_rates( &
      & nc, &                                     ! in
      & nsoil_sb, &
      & TRIM(dsl4jsb_Config(SB_)%bnf_scheme), &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & po4_solute(:,:), &
      & po4_assoc_fast(:,:), &
      & transport_nh4_solute(:,:), &
      & transport_no3_solute(:,:), &
      & transport_po4_solute(:,:), &
      & transport_po4_assoc_fast(:,:), &
      & weathering_po4(:,:), &
      & p_deposition_act(:,:), &
      & slow_exchange_po4(:,:), &
      & apparent_fast_som_cn(:,:), &
      & km_mic_uptake_nh4_act(:,:), &
      & km_mic_uptake_no3_act(:,:), &             ! in
      & plant_uptake_nh4_sl(:,:), &               ! inout
      & plant_uptake_no3_sl(:,:), &
      & plant_uptake_po4_sl(:,:), &
      & mycorrhiza_uptake_nh4_sl(:,:), &
      & mycorrhiza_uptake_no3_sl(:,:), &
      & mycorrhiza_uptake_po4_sl(:,:), &
      & nitrification_no3(:,:), &
      & denitrification_n2(:,:), &
      & asymb_n_fixation(:,:), &
      & biochem_mineralisation_po4(:,:), &
      & sb_loss_mt(:,:,:,:), &                    ! inout
      & microbial_uptake_nh4_sl(:,:), &           ! out
      & microbial_uptake_no3_sl(:,:), &           ! out
      & microbial_uptake_po4_sl(:,:) )            ! out

    !>  3.3 ...
    !>
    IF (TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme) == 'dynamic') THEN
      CALL calc_partitioning_nitrification_denitrification_rate( &
        & nc, &
        & nsoil_sb, &
        & rtm_denitrification(:,:), &
        & nitrification_no3(:,:), &
        & denitrification_n2(:,:), &
        & nitrification_noy(:,:), &
        & nitrification_n2o(:,:), &
        & denitrification_noy(:,:), &
        & denitrification_n2o(:,:))
    ENDIF

    !>  3.4 update n15 for inorganic fluxes 
    !>      here the function is called for the simple_1d soil-model: flux rates in [mol/m3/timestep]
    !>
    !>    NOTE: JSM is using the same function but uses flux rates in [micro-mol/m3/s]; unit conversion is applied within the routine
    !>
    CALL calc_sb_inorganic_n15_fluxes( &
      & nc, &                                         ! in
      & nsoil_sb, &
      & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
      & nh4_solute(:,:), &
      & nh4_n15_solute(:,:), &
      & no3_solute(:,:), &
      & no3_n15_solute(:,:), &
      & plant_uptake_nh4_sl(:,:), & 
      & plant_uptake_no3_sl(:,:), &      
      & mycorrhiza_uptake_nh4_sl(:,:), & 
      & mycorrhiza_uptake_no3_sl(:,:), & 
      & microbial_uptake_nh4_sl(:,:), &        
      & microbial_uptake_no3_sl(:,:), &        
      & asymb_n_fixation(:,:), &
      & nitrification_no3(:,:), &       
      & nitrification_noy(:,:), &        
      & nitrification_n2o(:,:), &        
      & denitrification_noy(:,:), &      
      & denitrification_n2o(:,:), &      
      & denitrification_n2(:,:), &       
      & transport_nh4_n15_solute(:,:), & 
      & transport_no3_n15_solute(:,:), &              ! in
      & plant_uptake_nh4_n15_sl(:,:), &               ! out
      & plant_uptake_no3_n15_sl(:,:), &      
      & mycorrhiza_uptake_nh4_n15_sl(:,:), & 
      & mycorrhiza_uptake_no3_n15_sl(:,:), & 
      & microbial_uptake_nh4_n15_sl(:,:), &        
      & microbial_uptake_no3_n15_sl(:,:), &        
      & asymb_n_fixation_n15(:,:), &
      & nitrification_no3_n15(:,:), &        
      & nitrification_noy_n15(:,:), &        
      & nitrification_n2o_n15(:,:), &        
      & denitrification_noy_n15(:,:), &      
      & denitrification_n2o_n15(:,:), &      
      & denitrification_n2_n15(:,:) )                 ! out
    
    !>  3.5 formation of SOM and associated respiration, net mineralisation fluxes
    !>
    CALL calc_SOM_formation( &
      & nc, &                                 ! in
      & nsoil_sb, &
      & apparent_fast_som_cn(:,:), &
      & apparent_fast_som_np(:,:), &
      & microbial_uptake_nh4_sl(:,:), &
      & microbial_uptake_nh4_n15_sl(:,:), &       
      & microbial_uptake_no3_sl(:,:), &
      & microbial_uptake_no3_n15_sl(:,:), &        
      & asymb_n_fixation(:,:), &
      & asymb_n_fixation_n15(:,:), &        
      & microbial_uptake_po4_sl(:,:), &       ! in
      & het_respiration(:,:), &               ! inout
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & net_mineralisation_nh4(:,:), &
      & net_mineralisation_no3(:,:), &
      & net_mineralisation_po4(:,:), &
      & net_mineralisation_nh4_n15(:,:), &
      & net_mineralisation_no3_n15(:,:), &
      & biochem_mineralisation_po4(:,:), &
      & sb_loss_mt(:,:,:,:), &
      & sb_formation_mt(:,:,:,:) )            ! inout

    !>  3.6 formation of metabolic and strucutral litter by woody liter decay and 
    !>     associated respiration, net mineralisation fluxes    
    !>
    IF (lctlib%growthform == ITREE) THEN 
      CALL calc_litter_formation_from_woody_decay( &
        & nc, &                               ! in
        & nsoil_sb, &
        & sb_loss_mt(:,:,:,:), &              ! in
        & het_respiration(:,:), &             ! inout
        & het_respiration_c13(:,:), &
        & het_respiration_c14(:,:), &
        & net_mineralisation_nh4(:,:), &
        & net_mineralisation_po4(:,:), &
        & net_mineralisation_nh4_n15(:,:), &
        & sb_formation_mt(:,:,:,:) )          ! inout
    ENDIF

    !>4.0 Mycorrhiza dynamics 
    !>
    CALL calc_mycorrhiza_dynamics( &
      & nc, &                                       ! in
      & nsoil_sb, &
      & soil_depth_sl(:,:), &
      & model%config%elements_index_map(:), &
      & model%config%is_element_used(:), &
      & lctlib%tau_mycorrhiza, &
      & dsl4jsb_Config(SB_)%flag_mycorrhiza_org, &
      & veg_pool_mt(ix_fine_root, ixC, :), &
      & sb_pool_mt(:,:,:,:), &                      ! in
      & mycorrhiza_uptake_nh4_sl(:,:), &            ! inout
      & mycorrhiza_uptake_nh4_n15_sl(:,:), &
      & mycorrhiza_uptake_no3_sl(:,:), &
      & mycorrhiza_uptake_no3_n15_sl(:,:), &
      & mycorrhiza_uptake_norg_sl(:,:), &
      & mycorrhiza_uptake_norg_n15_sl(:,:), &
      & mycorrhiza_uptake_po4_sl(:,:), &
      & nh4_solute(:,:), &
      & no3_solute(:,:), &
      & po4_solute(:,:), &
      & km2_uptake_nh4_act(:,:), &
      & km2_uptake_no3_act(:,:), &
      & km2_uptake_po4_act(:,:), &
      & km_uptake_org_act(:,:), &
      & het_respiration(:,:), &
      & het_respiration_c13(:,:), &
      & het_respiration_c14(:,:), &
      & myc_respiration(:,:), &
      & myc_respiration_c13(:,:), &
      & myc_respiration_c14(:,:), &
      & sb_formation_mt(:,:,:,:), &
      & sb_loss_mt(:,:,:,:), &
      & sb_mycorrhiza_export_mt(:,:,:) )            ! inout

    !>5.0 Update soil litter pools with litter fall
    !>
    CALL calc_litter_partitioning( &
      & nc, &                                         ! in
      & nsoil_sb, &
      & lctlib%sla, &
      & lctlib%growthform, &
      & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
      & soil_depth_sl(:,:), &
      & root_fraction_sl(:,:), &
      & veg_litterfall_mt(:,:,:), &                   ! in
      & sb_formation_mt(:,:,:,:) )                    ! inout
                
    !>6.0 for 'sb_model = fixed': estimate N losses
    !> 
    IF (TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme) == 'fixed') THEN
      CALL calc_fixed_nloss_rate( &
        & nc, &
        & nsoil_sb, &
        & soil_depth_sl(:,:), &
        & net_mineralisation_nh4(:,:), & 
        & net_mineralisation_nh4_n15(:,:), & 
        & nitrification_no3(:,:), &
        & nitrification_no3_n15(:,:), &
        & volatilisation_nh4(:,:), &
        & volatilisation_nh4_n15(:,:), &
        & emission_n2(:), &
        & emission_n2_n15(:))
    ENDIF

    !>7.0 unit conversion of rates and fluxes
    !>
    !>    convert rates to [micro-mol/ m3 / s]
    !>    covert surface fluxes to [micro-mol / m2 / s]
    !>
    het_respiration(:,:)     = het_respiration(:,:) / dtime * 1000000._wp
    het_respiration_c13(:,:) = het_respiration_c13(:,:) / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) / dtime * 1000000._wp
    myc_respiration(:,:)     = myc_respiration(:,:) / dtime * 1000000._wp
    myc_respiration_c13(:,:) = myc_respiration_c13(:,:) / dtime * 1000000._wp
    myc_respiration_c14(:,:) = myc_respiration_c14(:,:) / dtime * 1000000._wp

    net_mineralisation_nh4(:,:)    = net_mineralisation_nh4(:,:) / dtime * 1000000._wp
    net_mineralisation_no3(:,:)    = net_mineralisation_no3(:,:) / dtime * 1000000._wp
    net_mineralisation_po4(:,:)    = net_mineralisation_po4(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_sl(:,:)       = plant_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_no3_sl(:,:)       = plant_uptake_no3_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_po4_sl(:,:)       = plant_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_nh4_sl(:,:)  = mycorrhiza_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_no3_sl(:,:)  = mycorrhiza_uptake_no3_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_norg_sl(:,:) = mycorrhiza_uptake_norg_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_po4_sl(:,:)  = mycorrhiza_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    microbial_uptake_nh4_sl(:,:)   = microbial_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    microbial_uptake_no3_sl(:,:)   = microbial_uptake_no3_sl(:,:) / dtime * 1000000._wp 
    microbial_uptake_po4_sl(:,:)   = microbial_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    asymb_n_fixation(:,:)          = asymb_n_fixation(:,:) / dtime * 1000000._wp 

    volatilisation_nh4(:,:)  = volatilisation_nh4(:,:) / dtime * 1000000._wp
    nitrification_no3(:,:)   = nitrification_no3(:,:) / dtime * 1000000._wp 
    nitrification_noy(:,:)   = nitrification_noy(:,:) / dtime * 1000000._wp 
    nitrification_n2o(:,:)   = nitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_noy(:,:) = denitrification_noy(:,:) / dtime * 1000000._wp 
    denitrification_n2o(:,:) = denitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_n2(:,:)  = denitrification_n2(:,:) / dtime * 1000000._wp 

    transport_nh4_solute(:,:)     = transport_nh4_solute(:,:) / dtime * 1000000._wp
    transport_nh4_assoc(:,:)      = transport_nh4_assoc(:,:) / dtime * 1000000._wp
    transport_no3_solute(:,:)     = transport_no3_solute(:,:) / dtime * 1000000._wp
    transport_po4_solute(:,:)     = transport_po4_solute(:,:) / dtime * 1000000._wp
    transport_po4_assoc_fast(:,:) = transport_po4_assoc_fast(:,:) / dtime * 1000000._wp
    transport_po4_assoc_slow(:,:) = transport_po4_assoc_slow(:,:) / dtime * 1000000._wp
    transport_po4_occluded(:,:)   = transport_po4_occluded(:,:) / dtime * 1000000._wp
    transport_noy(:,:)            = transport_noy(:,:) / dtime * 1000000._wp
    transport_n2o(:,:)            = transport_n2o(:,:) / dtime * 1000000._wp
    transport_n2(:,:)             = transport_n2(:,:) / dtime * 1000000._wp

    lateral_loss_nh4_solute_sl(:,:) = lateral_loss_nh4_solute_sl(:,:) / dtime * 1000000._wp 
    lateral_loss_no3_solute_sl(:,:) = lateral_loss_no3_solute_sl(:,:) / dtime * 1000000._wp 
    lateral_loss_po4_solute_sl(:,:) = lateral_loss_po4_solute_sl(:,:) / dtime * 1000000._wp 

    leaching_nh4_solute(:) = leaching_nh4_solute(:) / dtime * 1000000._wp
    leaching_no3_solute(:) = leaching_no3_solute(:) / dtime * 1000000._wp
    leaching_po4_solute(:) = leaching_po4_solute(:) / dtime * 1000000._wp
    emission_noy(:)        = emission_noy(:) / dtime * 1000000._wp
    emission_n2o(:)        = emission_n2o(:) / dtime * 1000000._wp
    emission_n2(:)         = emission_n2(:) / dtime * 1000000._wp

    net_mineralisation_nh4_n15(:,:)   = net_mineralisation_nh4_n15(:,:) / dtime * 1000000._wp
    net_mineralisation_no3_n15(:,:)   = net_mineralisation_no3_n15(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_nh4_n15_sl(:,:)  = microbial_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_no3_n15_sl(:,:)  = microbial_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp
    asymb_n_fixation_n15(:,:)         = asymb_n_fixation_n15(:,:) / dtime * 1000000._wp 

    volatilisation_nh4_n15(:,:)  = volatilisation_nh4_n15(:,:) / dtime * 1000000._wp
    nitrification_no3_n15(:,:)   = nitrification_no3_n15(:,:) / dtime * 1000000._wp 
    nitrification_noy_n15(:,:)   = nitrification_noy_n15(:,:) / dtime * 1000000._wp 
    nitrification_n2o_n15(:,:)   = nitrification_n2o_n15(:,:) / dtime * 1000000._wp 
    denitrification_noy_n15(:,:) = denitrification_noy_n15(:,:) / dtime * 1000000._wp 
    denitrification_n2o_n15(:,:) = denitrification_n2o_n15(:,:) / dtime * 1000000._wp 
    denitrification_n2_n15(:,:)  = denitrification_n2_n15(:,:) / dtime * 1000000._wp 

    transport_nh4_n15_solute(:,:) = transport_nh4_n15_solute(:,:) / dtime * 1000000._wp
    transport_nh4_n15_assoc(:,:)  = transport_nh4_n15_assoc(:,:) / dtime * 1000000._wp
    transport_no3_n15_solute(:,:) = transport_no3_n15_solute(:,:) / dtime * 1000000._wp
    transport_noy_n15(:,:)        = transport_noy_n15(:,:) / dtime * 1000000._wp
    transport_n2o_n15(:,:)        = transport_n2o_n15(:,:) / dtime * 1000000._wp
    transport_n2_n15(:,:)         = transport_n2_n15(:,:) / dtime * 1000000._wp

    leaching_nh4_n15_solute(:)          = leaching_nh4_n15_solute(:) / dtime * 1000000._wp
    leaching_no3_n15_solute(:)          = leaching_no3_n15_solute(:) / dtime * 1000000._wp
    lateral_loss_nh4_n15_solute_sl(:,:) = lateral_loss_nh4_n15_solute_sl(:,:) / dtime * 1000000._wp 
    lateral_loss_no3_n15_solute_sl(:,:) = lateral_loss_no3_n15_solute_sl(:,:) / dtime * 1000000._wp 
    
    emission_noy_n15(:)        = emission_noy_n15(:) / dtime * 1000000._wp
    emission_n2o_n15(:)        = emission_n2o_n15(:) / dtime * 1000000._wp
    emission_n2_n15(:)         = emission_n2_n15(:) / dtime * 1000000._wp

    weathering_po4(:,:)             = weathering_po4(:,:) / dtime * 1000000._wp
    occlusion_po4(:,:)              = occlusion_po4(:,:) / dtime * 1000000._wp
    slow_exchange_po4(:,:)          = slow_exchange_po4(:,:) / dtime * 1000000._wp
    biochem_mineralisation_po4(:,:) = biochem_mineralisation_po4(:,:) / dtime * 1000000._wp

    !------------------------------------------------------------------------------------------------------ !
    ! de-allocate local allocatable variables
    DEALLOCATE(apparent_fast_som_cn)
    DEALLOCATE(apparent_fast_som_np)
    DEALLOCATE(fine_root_carbon_sl)
    DEALLOCATE(km1_uptake_nh4_act)
    DEALLOCATE(km2_uptake_nh4_act)
    DEALLOCATE(km1_uptake_no3_act)
    DEALLOCATE(km2_uptake_no3_act)
    DEALLOCATE(km1_uptake_po4_act)
    DEALLOCATE(km2_uptake_po4_act)
    DEALLOCATE(km_uptake_org_act)
    DEALLOCATE(km_mic_uptake_nh4_act)
    DEALLOCATE(km_mic_uptake_no3_act)
    DEALLOCATE(km_denitrification_c_act)
    DEALLOCATE(km_denitrification_no3_act)
    DEALLOCATE(p_deposition_act)
    DEALLOCATE(vmax_weath_mineral_act_sl)
    DEALLOCATE(k_adsorpt_po4_act)
    DEALLOCATE(k_desorpt_po4_act)
    DEALLOCATE(km_dip_biochem_po4_act)
    DEALLOCATE(qmax_po4_mineral_act)

  END SUBROUTINE update_sb_simple_model 

  ! ======================================================================================================= !
  !>provide soil profile values, i.e., local variables
  !>
  SUBROUTINE calc_local_soil_profile_values( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & bulk_dens_sl, &
    & w_soil_sl, &
    & nh4_solute, &
    & rtm_hsc, &
    & rmm_hsc, &
    & rtm_sorption, &
    & rtm_desorption, &
    & rmm_sorption, &
    & rmm_desorption, &
    & bulk_dens_corr_sl, &
    & p_deposition, &
    & vmax_weath_mineral_sl, &
    & km1_uptake_nh4_act, &
    & km2_uptake_nh4_act, &
    & km1_uptake_no3_act, &
    & km2_uptake_no3_act, &
    & km1_uptake_po4_act, &
    & km2_uptake_po4_act, &
    & km_uptake_org_act, &
    & km_mic_uptake_nh4_act, &
    & km_mic_uptake_no3_act, &
    & km_denitrification_c_act, &
    & km_denitrification_no3_act, &
    & apparent_fast_som_cn, &
    & apparent_fast_som_np, &
    & p_deposition_act, &
    & vmax_weath_mineral_act_sl, &
    & k_adsorpt_po4_act, &
    & k_desorpt_po4_act, &
    & km_dip_biochem_po4_act )

    USE mo_sb_constants,            ONLY: km_denitrification_c, km_denitrification_no3, km1_up_nh4, km1_up_no3, & 
                                          km1_up_po4, km2_up_nh4, km2_up_no3, km2_up_po4, km_up_org, &
                                          k_fast_som_cn_min, k_fast_som_cn_max, km_mic_uptake_nh4, km_mic_uptake_no3, &
                                          f_fast_som_cn, k_fast_som_np, kref_fast_som_np, k_desorpt_po4, &
                                          k_adsorpt_po4, km_dip_biochem_po4
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                                     INTENT(in)  :: nc                           !< dimensions
    INTEGER,                                     INTENT(in)  :: nsoil_sb                     !< number of soil layers
    REAL(wp), DIMENSION(nc),                     INTENT(in)  :: p_deposition                 !< p deposition, [micro-mol m-2 s-1]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: &
                                                                soil_depth_sl, &             !< soil depth of a layer [m]
                                                                bulk_dens_sl, &              !< bulk density of a layer [kg m-3] 
                                                                w_soil_sl, &                 !< water content of a layer [m] 
                                                                nh4_solute, &                !< NH4 content of a layer [mol m-3] 
                                                                rtm_hsc, &                   !< temperature response scalar for half-saturation constants 
                                                                rmm_hsc, &                   !< moisture response scalar for half-saturation constants
                                                                rtm_sorption, &              !< temperature response scalar for po4 sorption
                                                                rtm_desorption, &            !< temperature response scalar for po4 desorption
                                                                rmm_sorption, &              !< moisture response scalar for po4 sorption
                                                                rmm_desorption, &            !< moisture response scalar for po4 desorption
                                                                bulk_dens_corr_sl, &         !< bulk density of the soil layer, [kg/m3]
                                                                vmax_weath_mineral_sl        !< weathering coefficient of mineral soil [mol P m-3 s-1]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(out) :: &
                                                                km1_uptake_nh4_act, &        !< actual km1 for NH4 uptake [m3 / mol]
                                                                km2_uptake_nh4_act, &        !< actual km2 for NH4 uptake [mol / m3]
                                                                km1_uptake_no3_act, &        !< actual km1 for NO3 uptake [m3 / mol]
                                                                km2_uptake_no3_act, &        !< actual km2 for NO3 uptake [mol / m3]
                                                                km1_uptake_po4_act, &        !< actual km1 for PO4 uptake [m3 / mol]
                                                                km2_uptake_po4_act, &        !< actual km2 for PO4 uptake [mol / m3]
                                                                km_uptake_org_act, &         !< actual km for organic uptake [mol / m3]
                                                                km_mic_uptake_nh4_act, &     !< actual km for microbial NH4 uptake [mol / m3]
                                                                km_mic_uptake_no3_act, &     !< actual km for microbial NO3 uptake [m3 / mol]
                                                                km_denitrification_c_act, &  !< actual km on C for denitrification [mol / m3]
                                                                km_denitrification_no3_act, & !< actual km on NO3 for denitrification [mol / m3]
                                                                apparent_fast_som_cn, &           !< apparent microbial molar carbon-nitrogen ratio 
                                                                apparent_fast_som_np, &           !< apparent microbial molar carbon-nitrogen ratio 
                                                                p_deposition_act, &          !< actual p deposition, [mol m-3 timestep-1]
                                                                vmax_weath_mineral_act_sl, & !< actual weathering rate for primary P [mol P kg-1 s-1]
                                                                k_adsorpt_po4_act, &         !< actual adsorption rate for fast_assoc_po4 [s-1]
                                                                k_desorpt_po4_act, &         !< actual desorption rate for slow_assoc_po4 [s-1]
                                                                km_dip_biochem_po4_act    !< actual km for soluble P control on biochem_mineralisation_po4
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':calc_local_soil_profile_values'
    ! ----------------------------------------------------------------------------------------------------- !
    
    km1_uptake_nh4_act(:,:) = km1_up_nh4 / (rtm_hsc(:,:) * rmm_hsc(:,:))
    km1_uptake_no3_act(:,:) = km1_up_no3 / (rtm_hsc(:,:) * rmm_hsc(:,:))
    km1_uptake_po4_act(:,:) = km1_up_po4 / (rtm_hsc(:,:) * rmm_hsc(:,:)) / 1000._wp / (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km2_uptake_nh4_act(:,:) = km2_up_nh4 * rtm_hsc(:,:) * rmm_hsc(:,:)
    km2_uptake_no3_act(:,:) = km2_up_no3 * rtm_hsc(:,:) * rmm_hsc(:,:)
    km2_uptake_po4_act(:,:) = km2_up_po4 * rtm_hsc(:,:) * rmm_hsc(:,:) * 1000._wp * (w_soil_sl(:,:) / soil_depth_sl(:,:))
    km_uptake_org_act(:,:)          = km_up_org  * rtm_hsc(:,:) * rmm_hsc(:,:)
    km_mic_uptake_nh4_act(:,:)      = km_mic_uptake_nh4 / (rtm_hsc(:,:) * rmm_hsc(:,:))
    km_mic_uptake_no3_act(:,:)      = km_mic_uptake_no3 / (rtm_hsc(:,:) * rmm_hsc(:,:))
    km_dip_biochem_po4_act(:,:)     = km_dip_biochem_po4 * rtm_hsc(:,:) * rmm_hsc(:,:)
    km_denitrification_c_act(:,:)   = km_denitrification_c
    km_denitrification_no3_act(:,:) = km_denitrification_no3
    apparent_fast_som_cn(:,:) = MAX(k_fast_som_cn_max - f_fast_som_cn * nh4_solute(:,:) / bulk_dens_sl(:,:), &
      &                         k_fast_som_cn_min)
    apparent_fast_som_np(:,:) = k_fast_som_np * MAX(k_fast_som_cn_max - f_fast_som_cn * kref_fast_som_np / bulk_dens_sl(:,:), &
      &                                             k_fast_som_cn_min) / apparent_fast_som_cn(:,:)
    vmax_weath_mineral_act_sl(:,:) = rtm_hsc(:,:) / rmm_hsc(:,:) * vmax_weath_mineral_sl(:,:)
    k_adsorpt_po4_act(:,:) = rtm_sorption(:,:) * k_adsorpt_po4
    k_desorpt_po4_act(:,:) = rtm_desorption(:,:) * k_desorpt_po4
    p_deposition_act(:,:)  = 0.0_wp
    p_deposition_act(:,1)  = p_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
  END SUBROUTINE calc_local_soil_profile_values

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates potential plant uptake of NH4, NO3 and PO4, given soil concentrations,
  !! and fine root mass per layer 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_plant_uptake_rate(nc, &                         ! in
                                              nsoil_sb, &
                                              lctlib_vmax_uptake_n, & 
                                              lctlib_vmax_uptake_p, & 
                                              flag_mycorrhiza, &
                                              soil_depth_sl, &
                                              veg_pool_fine_root_carbon, &
                                              sb_pool_mycorrhiza_carbon, &
                                              root_fraction_sl, &
                                              rtm_root_uptake_act, &
                                              rmm_root_uptake_act, &
                                              km1_uptake_nh4_act, &
                                              km2_uptake_nh4_act, &
                                              km1_uptake_no3_act, &
                                              km2_uptake_no3_act, &
                                              km1_uptake_po4_act, &
                                              km2_uptake_po4_act, &
                                              nh4_solute, &
                                              no3_solute, &
                                              po4_solute, &
                                              plant_uptake_nh4_sl, &        ! out
                                              plant_uptake_no3_sl, &
                                              plant_uptake_po4_sl, &
                                              unit_uptake_n_pot, &
                                              unit_uptake_p_pot  )

    USE mo_sb_constants,            ONLY: f_surface_myc2fineroot, myc_root_coverage_max
    
    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                     INTENT(in)  :: nc, &                      !< dimensions
                                                                nsoil_sb
    REAL(wp),                                    INTENT(in)  :: lctlib_vmax_uptake_n       !< lctlib parameter
    REAL(wp),                                    INTENT(in)  :: lctlib_vmax_uptake_p       !< lctlib parameter
    LOGICAL,                                     INTENT(in)  :: flag_mycorrhiza            !< from sb config
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: soil_depth_sl              !< soil depth per layer [m]
    REAL(wp), DIMENSION(nc),                     INTENT(in)  :: veg_pool_fine_root_carbon  !< fine root carbon [mol / m2]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: sb_pool_mycorrhiza_carbon  !< mycorrhiza carbon [mol / m2]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(in)  :: root_fraction_sl, &        !< fraction of fine roots per soil layer
                                                                rtm_root_uptake_act, &     !< temperature response of root uptake 
                                                                rmm_root_uptake_act, &     !< moisture response of root uptake 
                                                                km1_uptake_nh4_act, &      !< actual km1 for NH4 uptake [m3 / mol]
                                                                km2_uptake_nh4_act, &      !< actual km2 for NH4 uptake [mol / m3]
                                                                km1_uptake_no3_act, &      !< actual km1 for NO3 uptake [m3 / mol]
                                                                km2_uptake_no3_act, &      !< actual km2 for NO3 uptake [mol / m3]
                                                                km1_uptake_po4_act, &      !< actual km1 for PO4 uptake [m3 / mol]
                                                                km2_uptake_po4_act, &      !< actual km2 for PO4 uptake [mol / m3]
                                                                nh4_solute, &              !< NH4 solute concentration [mol / m3] 
                                                                no3_solute, &              !< NO3 solute concentration [mol / m3]
                                                                po4_solute                 !< PO4 solute concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),           INTENT(out) :: &
                                                                plant_uptake_nh4_sl, &     !< potential plant uptake of NH4 [mol / m3 / timestep]
                                                                plant_uptake_no3_sl, &     !< potential plant uptake of NO3 [mol / m3 / timestep]
                                                                plant_uptake_po4_sl        !< potential plant uptake of PO4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc),                     INTENT(out) :: unit_uptake_n_pot, &       !< potential N uptake by fine roots per unit mass [micro-mol / mol / s]
                                                                unit_uptake_p_pot          !< potential P uptake by fine roots per unit mass [micro-mol / mol / s]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: isoil        !< loop over nsoil_sb
    REAL(wp), DIMENSION(nc)                          :: myc_root_coverage_act, f_vmax_uptake_act
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                      :: routine = TRIM(modname)//':calc_potential_plant_uptake_rate'



    !> 0.9 zero output variables
    !!
    unit_uptake_n_pot(:) = 0.0_wp
    unit_uptake_p_pot(:) = 0.0_wp


    !> 1.0 inorganic plant N uptake 
    !!
    !! only root tips without mycorrhiza are allowed to take up N // max 30% of roots are infected! 'longer' myco if limited
    IF(flag_mycorrhiza)THEN
       WHERE(veg_pool_fine_root_carbon(:) > eps8) 
         ! myc_root_coverage_act(:) = myc_root_coverage_max
         ! dynamical root coverage to account for explorative mycorrhizae (idea for later...)  
         myc_root_coverage_act(:) = MIN(myc_root_coverage_max, &
                                    SUM(sb_pool_mycorrhiza_carbon(:,:) * soil_depth_sl(:,:), DIM=2) / veg_pool_fine_root_carbon(:))
         f_vmax_uptake_act(:)     =  1.0_wp / (1.0_wp - myc_root_coverage_max + myc_root_coverage_max * f_surface_myc2fineroot)
       ELSEWHERE
          myc_root_coverage_act(:) = 0.0_wp
          f_vmax_uptake_act(:)     = 1.0_wp
       ENDWHERE
    ELSE  
       f_vmax_uptake_act(:)     = 1.0_wp 
       myc_root_coverage_act(:) = 0.0_wp
    ENDIF

    DO isoil = 1,nsoil_sb
       !> 2.1 inorganic plant N uptake 
       !!
       !! potential N uptake rate (micro-mol N / mol C fine-root / s), weighted by root distribution
       hlp1(:,isoil) = root_fraction_sl(:,isoil) * rtm_root_uptake_act(:,isoil) * rmm_root_uptake_act(:,isoil) * &
                       (1.0_wp - myc_root_coverage_act(:)) * lctlib_vmax_uptake_n * f_vmax_uptake_act(:) * &
                       nh4_solute(:,isoil) * ( km1_uptake_nh4_act(:,isoil) + & 
                       1._wp / ( km2_uptake_nh4_act(:,isoil) + nh4_solute(:,isoil) )) 
       hlp2(:,isoil) = root_fraction_sl(:,isoil) * rtm_root_uptake_act(:,isoil) * rmm_root_uptake_act(:,isoil) * &
                       (1.0_wp - myc_root_coverage_act(:)) * lctlib_vmax_uptake_n * f_vmax_uptake_act(:) * &
                       no3_solute(:,isoil) * ( km1_uptake_no3_act(:,isoil) + & 
                       1._wp / ( km2_uptake_no3_act(:,isoil) + no3_solute(:,isoil) )) 
       unit_uptake_n_pot(:) = unit_uptake_n_pot(:) + hlp1(:,isoil) + hlp2(:,isoil) 
 
       ! potential uptake rate in mol N / m3 / timestep
       plant_uptake_nh4_sl(:,isoil) = hlp1(:,isoil) * &
                                      veg_pool_fine_root_carbon(:) / & 
                                      soil_depth_sl(:,isoil) * dtime / 1.e6_wp
       plant_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * &
                                      veg_pool_fine_root_carbon(:) / & 
                                      soil_depth_sl(:,isoil) * dtime / 1.e6_wp
 
       !> 2.2 inorganic plant P uptake 
       !!
       !! potential P uptake rate (micro-mol P / mol C fine-root / s), weighted by root distribution
       hlp1(:,isoil) = root_fraction_sl(:,isoil)  * rtm_root_uptake_act(:,isoil) * rmm_root_uptake_act(:,isoil) * &
                       (1.0_wp - myc_root_coverage_act(:)) * lctlib_vmax_uptake_p * f_vmax_uptake_act(:) * &
                        po4_solute(:,isoil) * ( km1_uptake_po4_act(:,isoil) + &
                        1._wp / ( km2_uptake_po4_act(:,isoil) + po4_solute(:,isoil) ))
       unit_uptake_p_pot(:) = unit_uptake_p_pot(:) + hlp1(:,isoil)
       ! potential uptake rate in mol P / m3 / timestep
       plant_uptake_po4_sl(:,isoil) = hlp1(:,isoil) * veg_pool_fine_root_carbon(:) / soil_depth_sl(:,isoil) * dtime / 1.e6_wp
   
    END DO    

  END SUBROUTINE calc_potential_plant_uptake_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates potential mycorrhiza uptake of NH4, NO3 and PO4, given soil concentrations,
  !! and mycorrhiza mass per layer 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_mycorrhiza_uptake_rate(nc, &                         ! in
                                                   nsoil_sb, &
                                                   lctlib_vmax_uptake_n, & 
                                                   lctlib_vmax_uptake_p, & 
                                                   soil_depth_sl, &
                                                   veg_pool_fine_root_carbon, &
                                                   sb_pool_mycorrhiza_carbon, &
                                                   sb_pool_mycorrhiza_nitrogen, &
                                                   rtm_mycorrhiza_uptake_act, &
                                                   rmm_mycorrhiza_uptake_act, &
                                                   km1_uptake_nh4_act, &
                                                   km2_uptake_nh4_act, &
                                                   km1_uptake_no3_act, &
                                                   km2_uptake_no3_act, &
                                                   km1_uptake_po4_act, &
                                                   km2_uptake_po4_act, &
                                                   nh4_solute, &
                                                   no3_solute, &
                                                   po4_solute, &
                                                   mycorrhiza_uptake_nh4_sl, &   ! out
                                                   mycorrhiza_uptake_no3_sl, &
                                                   mycorrhiza_uptake_po4_sl )
  
    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: f_maxresp_myc, f_surface_myc2fineroot, myc_root_coverage_max

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)  :: nc                           !< dimensions
    INTEGER,                                          INTENT(in)  :: nsoil_sb                     !< number of soil layers
    REAL(wp),                                         INTENT(in)  :: lctlib_vmax_uptake_n         !< lctlib parameter
    REAL(wp),                                         INTENT(in)  :: lctlib_vmax_uptake_p         !< lctlib parameter
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: soil_depth_sl                !< soil depth per layer [m]    
    REAL(wp), DIMENSION(nc),                          INTENT(in)  :: veg_pool_fine_root_carbon    !< fine root carbon [mol / m2]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: sb_pool_mycorrhiza_carbon    !< mycorrhizal biomass [mol C / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: sb_pool_mycorrhiza_nitrogen  !< mycorrhizal biomass [mol C / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: rtm_mycorrhiza_uptake_act, & !< temperature response of mycorrhiza uptake 
                                                                     rmm_mycorrhiza_uptake_act, & !< moisture response of mycorrhiza uptake 
                                                                     km1_uptake_nh4_act, &        !< actual km1 for NH4 uptake [m3 / mol]
                                                                     km2_uptake_nh4_act, &        !< actual km2 for NH4 uptake [mol / m3]
                                                                     km1_uptake_no3_act, &        !< actual km1 for NO3 uptake [m3 / mol]
                                                                     km2_uptake_no3_act, &        !< actual km2 for NO3 uptake [mol / m3]
                                                                     km1_uptake_po4_act, &        !< actual km1 for PO4 uptake [m3 / mol]
                                                                     km2_uptake_po4_act, &        !< actual km2 for PO4 uptake [mol / m3]
                                                                     nh4_solute, &                !< NH4 solute concentration [mol / m3] 
                                                                     no3_solute, &                !< NO3 solute concentration [mol / m3]
                                                                     po4_solute                   !< PO4 solute concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out) :: mycorrhiza_uptake_nh4_sl, &  !< potential mycorrhiza uptake of NH4 [mol / m3 / timestep]
                                                                     mycorrhiza_uptake_no3_sl, &  !< potential mycorrhiza uptake of NO3 [mol / m3 / timestep]
                                                                     mycorrhiza_uptake_po4_sl     !< potential mycorrhiza uptake of PO4 [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    INTEGER                                          :: isoil        !< loop over nsoil_sb
    REAL(wp), DIMENSION(nc)                          :: f_vmax_uptake_act
    REAL(wp), DIMENSION(nc,nsoil_sb)                 :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                      :: routine = TRIM(modname)//':calc_potential_mycorrhiza_uptake_rate'

    !> 0.9.1 init out arguments
    !>
    mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
    mycorrhiza_uptake_no3_sl(:,:) = 0.0_wp
    mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp
    
    !> 0.9.2 init local var
    !>
    f_vmax_uptake_act(:) = 0.0_wp

    !> 1.0 pot. N uptake by mycorrhizal fungi (mol / m3 / timestep)
    !>
    WHERE(veg_pool_fine_root_carbon(:) > eps8) 
       f_vmax_uptake_act(:) =  1.0_wp / (1.0_wp - myc_root_coverage_max + myc_root_coverage_max * f_surface_myc2fineroot)    
    ELSEWHERE
       f_vmax_uptake_act(:) = 1.0_wp
    ENDWHERE
    DO isoil = 1,nsoil_sb
      mycorrhiza_uptake_nh4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * rmm_mycorrhiza_uptake_act(:,isoil) * &
                                          f_surface_myc2fineroot * lctlib_vmax_uptake_n * f_vmax_uptake_act(:) * &
                                          nh4_solute(:,isoil) * ( km1_uptake_nh4_act(:,isoil) + & 
                                          1._wp / ( km2_uptake_nh4_act(:,isoil) + nh4_solute(:,isoil) )) * & 
                                          dtime / 1.e6_wp
      mycorrhiza_uptake_no3_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * rmm_mycorrhiza_uptake_act(:,isoil) * & 
                                          f_surface_myc2fineroot * lctlib_vmax_uptake_n * f_vmax_uptake_act(:) * &
                                          no3_solute(:,isoil) * ( km1_uptake_no3_act(:,isoil) + & 
                                          1._wp / ( km2_uptake_no3_act(:,isoil) + no3_solute(:,isoil) )) * & 
                                          dtime / 1.e6_wp
      !> 1.1 constrain potential N uptake by avaiable C for N transformation to organic
      !!
      hlp1(:,isoil) = mycorrhiza_uptake_nh4_sl(:,isoil) * transform_cost_nh4 + & 
                      mycorrhiza_uptake_no3_sl(:,isoil) * transform_cost_no3 
      WHERE(hlp1(:,isoil)  > (sb_pool_mycorrhiza_carbon(:,isoil) * f_maxresp_myc / one_day * dtime))
        WHERE(hlp1(:,isoil) > eps8)
          hlp2(:,isoil)                     = sb_pool_mycorrhiza_carbon(:,isoil) * f_maxresp_myc / one_day * dtime / & 
                                              (hlp1(:,isoil)) 
          mycorrhiza_uptake_nh4_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_nh4_sl(:,isoil)
          mycorrhiza_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_no3_sl(:,isoil)
        ELSEWHERE
          mycorrhiza_uptake_nh4_sl(:,isoil) = 0.0_wp
          mycorrhiza_uptake_no3_sl(:,isoil) = 0.0_wp
        ENDWHERE
      ENDWHERE 
      !> 2.0 pot. P uptake by mycorrhizal fungi (mol / m3 / timestep)
      !!
      mycorrhiza_uptake_po4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * rmm_mycorrhiza_uptake_act(:,isoil) *  & 
                                          f_surface_myc2fineroot * lctlib_vmax_uptake_p * f_vmax_uptake_act(:) * & 
                                          po4_solute(:,isoil) * ( km1_uptake_po4_act(:,isoil) + & 
                                          1._wp / ( km2_uptake_po4_act(:,isoil) + po4_solute(:,isoil) )) * & 
                                          dtime / 1.e6_wp
    END DO
  END SUBROUTINE calc_potential_mycorrhiza_uptake_rate
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the potential rates of nitrification and denitrification 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_nitrification_denitrification_rate(nc, &                            ! in
                                                               nsoil_sb, &
                                                               rtm_nitrification, &
                                                               rmm_nitrification, &
                                                               rtm_denitrification, &
                                                               km_denitrification_c_act, &
                                                               km_denitrification_no3_act, &
                                                               anaerobic_volume_fraction, &
                                                               nh4_solute, &
                                                               no3_solute, &
                                                               sb_pool_microbial_carbon, &
                                                               nitrification_no3, &               ! out
                                                               denitrification_n2)

    USE mo_sb_constants,     ONLY: vmax_nitrification, vmax_denitrification

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                                          INTENT(in)  :: nc                             !< dimensions
    INTEGER,                                          INTENT(in)  :: nsoil_sb                       !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: rtm_nitrification, &           !< actual temperature response for nitrification
                                                                     rmm_nitrification, &           !< actual moisture response for nitrification
                                                                     rtm_denitrification, &         !< actual temperature response for denitrification 
                                                                     km_denitrification_c_act, &    !< actual km on C for denitrification
                                                                     km_denitrification_no3_act, &  !< actual km on NO3 for denitrification
                                                                     anaerobic_volume_fraction, &   !< actual anaerobic volume fraction
                                                                     nh4_solute, &                  !< NH4 in solution [mol / m3]
                                                                     no3_solute                     !< NO3 in solution [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: sb_pool_microbial_carbon       !< microbial C [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out) :: nitrification_no3, &           !< potential nitrification rate [mol / m3 / timestep]
                                                                     denitrification_n2             !< potential denitrification rate [mol / m3 / timestep] 
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_nitrification_denitrification_rate'

    ! potential nitrification rate from NH4 to NO3 [mol/m3/timestep]
    nitrification_no3(:,:)   = vmax_nitrification * rtm_nitrification(:,:) * rmm_nitrification(:,:) * & 
                               sb_pool_microbial_carbon(:,:) / &
                               (km_denitrification_c_act(:,:) + sb_pool_microbial_carbon(:,:))  * &
                               (1._wp - anaerobic_volume_fraction(:,:)) * nh4_solute(:,:) * dtime / 1.e6_wp 

    ! potential denitrification rate from NO3 to N2 [mol/m3/timestep]
    denitrification_n2(:,:)  = vmax_denitrification * rtm_denitrification(:,:) * &
                               sb_pool_microbial_carbon(:,:) / &
                               (km_denitrification_c_act(:,:) + sb_pool_microbial_carbon(:,:))  * &
                               anaerobic_volume_fraction(:,:) * no3_solute(:,:) / &
                               (km_denitrification_no3_act(:,:) + no3_solute(:,:)) * dtime / 1.e6_wp 

  END SUBROUTINE calc_potential_nitrification_denitrification_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the leaching and gas transport within the soil column, 
  !! given current soil N and P concentrations, soil temperature, moisture and water percolation rate 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_inorganic_transport_rate(nc, &                              ! in
                                           nsoil_sb, &
                                           soil_depth_sl, &
                                           percolation_sl, &
                                           frac_w_lat_loss_sl, & 
                                           rtm_gasdiffusion, &  
                                           rmm_gasdiffusion, & 
                                           nh4_solute, &
                                           nh4_n15_solute, &
                                           no3_solute, &
                                           no3_n15_solute, &
                                           po4_solute, &
                                           noy, & 
                                           noy_n15, & 
                                           n2o, &
                                           n2o_n15, &
                                           n2, &
                                           n2_n15, &
                                           transport_nh4_solute, &            ! out
                                           transport_nh4_n15_solute, &
                                           transport_no3_solute, &
                                           transport_no3_n15_solute, &
                                           transport_po4_solute, &
                                           transport_noy, &     
                                           transport_noy_n15, &   
                                           transport_n2o, &     
                                           transport_n2o_n15, & 
                                           transport_n2, &      
                                           transport_n2_n15, &     
                                           leaching_nh4_solute, &
                                           leaching_nh4_n15_solute, &
                                           leaching_no3_solute, &
                                           leaching_no3_n15_solute, &
                                           leaching_po4_solute, &
                                           lateral_loss_nh4_solute_sl, & 
                                           lateral_loss_nh4_n15_solute_sl, & 
                                           lateral_loss_no3_solute_sl, & 
                                           lateral_loss_no3_n15_solute_sl, & 
                                           lateral_loss_po4_solute_sl, & 
                                           emission_noy, &    
                                           emission_noy_n15, &
                                           emission_n2o, &    
                                           emission_n2o_n15, & 
                                           emission_n2, &  
                                           emission_n2_n15)

    USE mo_sb_constants,     ONLY: fleach_nh4, fleach_po4
    USE mo_sb_jsm_transport, ONLY: calc_liquid_phase_transport

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                           INTENT(in)    :: nc                                 !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                           !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl                      !< soil layer depth [m]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: percolation_sl, &                  !< fraction of water lost for 1d model [s-1]
                                                        frac_w_lat_loss_sl, &              !< constrained fraction of lateral (horizontal) water loss of 'w_soil_sl_old' (prev. timestep)
                                                        rtm_gasdiffusion, &                !< temperature modifier for gas diffusion
                                                        rmm_gasdiffusion, &                !< moisture modifier for gas diffusion       
                                                        nh4_solute, &                      !< NH4 in solution [mol / m3]
                                                        nh4_n15_solute, &                  !< 15NH4 in solution [mol / m3]
                                                        no3_solute, &                      !< NO3 in solution [mol / m3]
                                                        no3_n15_solute, &                  !< 15NO3 in solution [mol / m3]
                                                        po4_solute, &                      !< PO4 in solution [mol / m3]
                                                        noy, &                             !< NOy concentration [mol / m3]
                                                        noy_n15, &                         !< 15NOy concentration [mol / m3]
                                                        n2o, &                             !< N2O concentration [mol / m3]
                                                        n2o_n15, &                         !< 15N2O concentration [mol / m3]
                                                        n2, &                              !< N2 concentration [mol / m3]
                                                        n2_n15                             !< 15N2 concentration [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_solute, &            !< vertical transport of NH4 inside soil column [mol / m3 / timestep]
                                                        transport_nh4_n15_solute, &        !< vertical transport of 15NH4 inside soil column [mol / m3 / timestep]
                                                        transport_no3_solute, &            !< vertical transport of NO3 inside soil column [mol / m3 / timestep]
                                                        transport_no3_n15_solute, &        !< vertical transport of 15NO3 inside soil column [mol / m3 / timestep]
                                                        transport_po4_solute               !< vertical transport of PO4 inside soil column [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: transport_noy, &                   !< vertical transport of NOy [mol / m3 / timestep]
                                                        transport_noy_n15, &               !< vertical transport of 15NOy [mol / m3 / timestep]
                                                        transport_n2o, &                   !< vertical transport of N2O [mol / m3 / timestep]
                                                        transport_n2o_n15, &               !< vertical transport of 15N2O [mol / m3 / timestep]
                                                        transport_n2, &                    !< vertical transport of N2 [mol / m3 / timestep]
                                                        transport_n2_n15, &                !< vertical transport of 15N2 [mol / m3 / timestep]
                                                        lateral_loss_nh4_solute_sl, &      !< lateral loss of NH4 [mol/ m2 / timestep]  
                                                        lateral_loss_nh4_n15_solute_sl, &  !< lateral loss of 15NH4 [mol / m2 / timestep] 
                                                        lateral_loss_no3_solute_sl, &      !< lateral loss of NO3 [mol/ m2 / timestep] 
                                                        lateral_loss_no3_n15_solute_sl, &  !< lateral loss of 15NO3 [mol / m2 / timestep] 
                                                        lateral_loss_po4_solute_sl         !< lateral loss of PO4 [mol / m2 / timestep]
    REAL(wp), DIMENSION(nc),           INTENT(out)   :: leaching_nh4_solute, &             !< leaching of NH4 [mol / m2 / timestep]
                                                        leaching_nh4_n15_solute, &         !< leaching of 15NH4 [mol / m2 / timestep]
                                                        leaching_no3_solute, &             !< leaching of NO3 [mol / m2 / timestep]
                                                        leaching_no3_n15_solute, &         !< leaching of 15NO3 [mol / m2 / timestep]
                                                        leaching_po4_solute, &             !< leaching of PO4 [mol / m2 / timestep]
                                                        emission_noy, &                    !< soil efflux of NOy [mol / m2 / timestep]
                                                        emission_noy_n15, &                !< soil efflux of 15NOy [mol / m2 / timestep]
                                                        emission_n2o, &                    !< soil efflux of N2O [mol / m2 / timestep]
                                                        emission_n2o_n15, &                !< soil efflux of 15N2O [mol / m2 / timestep]
                                                        emission_n2, &                     !< soil efflux of N2 [mol / m2 / timestep]
                                                        emission_n2_n15                    !< soil efflux of 15N2 [mol / m2 / timestep]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_inorganic_transport_rate'



    !>1.0. Leaching rate
    !>
    ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
    ! reduction of layer concentration in mol m-3 timestep-1, leaching rate in mol m-2 timestep-1
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     nh4_solute(:,:), transport_nh4_solute(:,:), leaching_nh4_solute(:), &
                                     lateral_loss_nh4_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     nh4_n15_solute(:,:), transport_nh4_n15_solute(:,:), leaching_nh4_n15_solute(:), &
                                     lateral_loss_nh4_n15_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     no3_solute(:,:), transport_no3_solute(:,:), leaching_no3_solute(:), &
                                     lateral_loss_no3_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     no3_n15_solute(:,:), transport_no3_n15_solute(:,:), leaching_no3_n15_solute(:), &
                                     lateral_loss_no3_n15_solute_sl(:,:))
    CALL calc_liquid_phase_transport(nc, nsoil_sb, percolation_sl(:,:), frac_w_lat_loss_sl(:,:), soil_depth_sl(:,:), &
                                     po4_solute(:,:), transport_po4_solute(:,:), leaching_po4_solute(:), &
                                     lateral_loss_po4_solute_sl(:,:))
  
    !> 2.0 gaseous loss of NOy, N2O, and N2, for now following Xu-Ri and Prentice 2008, 
    !!thus ignoring any actual gas-transport and instead assume that from each layer a certain fraction excapes directly to the atmosphere.
    !! [mol m-3 timestep-1] 
    transport_noy(:,:)     = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * noy(:,:)     / one_day * dtime
    transport_noy_n15(:,:) = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * noy_n15(:,:) / one_day * dtime
    transport_n2o(:,:)     = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2o(:,:)     / one_day * dtime 
    transport_n2o_n15(:,:) = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2o_n15(:,:) / one_day * dtime 
    transport_n2(:,:)      = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2(:,:)      / one_day * dtime 
    transport_n2_n15(:,:)  = -rtm_gasdiffusion(:,:) * rmm_gasdiffusion(:,:) * n2_n15(:,:)  / one_day * dtime 

    ! emission rates as positive flux to the atmosphere, integrated over the soil column [mol m-2 timestep-1]
    emission_noy(:)     = SUM(-transport_noy(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2o(:)     = SUM(-transport_n2o(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2(:)      = SUM(-transport_n2(:,:)      * soil_depth_sl(:,:), DIM=2)
    emission_noy_n15(:) = SUM(-transport_noy_n15(:,:) * soil_depth_sl(:,:), DIM=2)
    emission_n2o_n15(:) = SUM(-transport_n2o_n15(:,:) * soil_depth_sl(:,:), DIM=2)
    emission_n2_n15(:)  = SUM(-transport_n2_n15(:,:)  * soil_depth_sl(:,:), DIM=2)

  END SUBROUTINE calc_inorganic_transport_rate
  
  ! ======================================================================================================= !
  !>Calculate potential decomposition rate of SOM and litter pools according to:
  !>  `dX/dt = - f(T) * g(M) / tau_pool * X`
  !>
  SUBROUTINE calc_potential_decomposition_rate( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & flag_mycorrhiza_prim, &
    & flag_mycorrhiza_org, &
    & rtm_decomposition, &
    & rmm_decomposition, &
    & sb_pool_mt, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & myc_respiration, &
    & myc_respiration_c13, &
    & myc_respiration_c14, &
    & sb_loss_mt )

    USE mo_jsb_math_constants,      ONLY: eps8  
    USE mo_sb_constants,            ONLY: tau_metabolic_litter, tau_structural_litter, tau_woody_litter, tau_fast, tau_slow, &
      &                                   mycorrhiza_cn, vmax_priming, km_priming, resp_priming, f_tau_myc
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                            INTENT(in)    :: nc                       !< dimensions
    INTEGER,                            INTENT(in)    :: nsoil_sb                 !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: soil_depth_sl            !< soil layer depth [m]
    INTEGER,                            INTENT(in)    :: elements_index_map(:)    !< map bgcm element ID -> IDX
    LOGICAL,                            INTENT(in)    :: is_element_used(:)       !< is element in 'elements_index_map' used 
    LOGICAL,                            INTENT(in)    :: flag_mycorrhiza_prim, &  !< from sb config
                                                         flag_mycorrhiza_org      !< from sb config
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: rtm_decomposition, &     !< temperature response of decomposition
                                                         rmm_decomposition        !< moisture response of decomposition 
    REAL(wp),                           INTENT(in)    :: sb_pool_mt(:,:,:,:)      !< sb_pool [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(inout) :: het_respiration, &       !< heterotrophic respiration [mol C / m3 / timestep]
                                                         het_respiration_c13, &   !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                         het_respiration_c14, &   !< heterotrophic respiration [micro-mol 14C / m3 / timestep]
                                                         myc_respiration, &       !< mycorrhizae respiration [mol C / m3 / timestep]
                                                         myc_respiration_c13, &   !< mycorrhizae respiration [mol 13C / m3 / timestep]
                                                         myc_respiration_c14      !< mycorrhizae respiration [micro-mol 14C / m3 / timestep]
    REAL(wp),                           INTENT(inout) :: sb_loss_mt(:,:,:,:)      !< sb_loss flux [mol / m3  ??]
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                           :: ielem                !< loop over bgcm elements
    INTEGER                           :: ix_elem              !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc,nsoil_sb)  :: tau_slow_myc, &
                                         km_priming_act, &
                                         resp_priming_act, &
                                         tau_slow_act
    REAL(wp), DIMENSION(nc,nsoil_sb)  :: hlp1, hlp2, hlp3
    REAL(wp), DIMENSION(nc,nsoil_sb)  :: hlp_cn
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_decomposition_rate'
    ! ----------------------------------------------------------------------------------------------------- !

    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
    
        !>1.0 metabolic litter (= soluable pool)
        !>
        sb_loss_mt(ix_soluable_litter, ix_elem, :, :) = sb_loss_mt(ix_soluable_litter, ix_elem, :, :) &
          &                                             + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                             / tau_metabolic_litter &
          &                                             * dtime * sb_pool_mt(ix_soluable_litter, ix_elem, :, :)
    
        !>2.0 structural litter (= polymeric pool)
        !>
        sb_loss_mt(ix_polymeric_litter, ix_elem, :, :) = sb_loss_mt(ix_polymeric_litter, ix_elem, :, :) &
          &                                              + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                              / tau_structural_litter &
          &                                              * dtime * sb_pool_mt(ix_polymeric_litter, ix_elem, :, :)
        
        !>3.0 woody litter
        !>
        sb_loss_mt(ix_woody_litter, ix_elem, :, :) = sb_loss_mt(ix_woody_litter, ix_elem, :, :) &
          &                                          + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                          / tau_woody_litter &
          &                                          * dtime * sb_pool_mt(ix_woody_litter, ix_elem, :, :)
        
        !>4.0 fast SOM pool (= microbial pool)
        !>
        sb_loss_mt(ix_microbial, ix_elem, :, :) = sb_loss_mt(ix_microbial, ix_elem, :, :) &
          &                                       + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                       / tau_fast &
          &                                       * dtime * sb_pool_mt(ix_microbial, ix_elem, :, :)
      END IF  ! IF element is used
    END DO    ! loop over elements

    !>5.0 slow SOM pool (= residue pool)
    !>
    ! mycorrhiza_prim
    IF(flag_mycorrhiza_prim)THEN 
      km_priming_act(:,:)   = km_priming
      resp_priming_act(:,:) = resp_priming
      WHERE(sb_pool_mt(ix_mycorrhiza, ixC, :, :) > eps8)
          tau_slow_act(:,:) = tau_slow / ((sb_pool_mt(ix_mycorrhiza, ixC, :, :) * vmax_priming) / &
                              (km_priming_act(:,:) + sb_pool_mt(ix_mycorrhiza, ixC, :, :)))
      ELSEWHERE
          hlp1(:,:)         = eps8
          tau_slow_act(:,:) = tau_slow /((hlp1(:,:) * vmax_priming) / (km_priming_act(:,:) + hlp1(:,:)))
      ENDWHERE

      sb_loss_mt(ix_mycorrhiza, ixC, :, :)   = sb_loss_mt(ix_mycorrhiza, ixC, :, :)   + resp_priming_act(:,:) &
        &                                      * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      sb_loss_mt(ix_mycorrhiza, ixC13, :, :) = sb_loss_mt(ix_mycorrhiza, ixC13, :, :) + resp_priming_act(:,:) &
        &                                      * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
      sb_loss_mt(ix_mycorrhiza, ixC14, :, :) = sb_loss_mt(ix_mycorrhiza, ixC14, :, :) + resp_priming_act(:,:) & 
        &                                      * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

      het_respiration(:,:)     = het_respiration(:,:)     + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      het_respiration_c13(:,:) = het_respiration_c13(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
      het_respiration_c14(:,:) = het_respiration_c14(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

      myc_respiration(:,:)     = myc_respiration(:,:)     + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      myc_respiration_c13(:,:) = myc_respiration_c13(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
      myc_respiration_c14(:,:) = myc_respiration_c14(:,:) + resp_priming_act(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)
    ! mycorrhiza_org
    ELSEIF(flag_mycorrhiza_org)THEN
      tau_slow_act(:,:) = tau_slow * f_tau_myc 
    ! other
    ELSE
      tau_slow_act(:,:) = tau_slow
    END IF 

    ! ...
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        sb_loss_mt(ix_residue, ix_elem, :, :) = sb_loss_mt(ix_residue, ix_elem, :, :) &
          &                                     + rtm_decomposition(:,:) * rmm_decomposition(:,:) &
          &                                     / tau_slow_act(:,:) &
          &                                     * dtime * sb_pool_mt(ix_residue, ix_elem, :, :)
      END IF  ! IF element is used
    END DO    ! loop over elements
  END SUBROUTINE calc_potential_decomposition_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculate potential biomineralisation rate 
  !!  dP/dt = - f(T)*g(M)/tau_pool * P
  !> Calculate weathering rate, slow adsorption rate and occlusion rate based on Yang et al. 2014
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_p_flux_rate(nc, &
                                        nsoil_sb, &
                                        flag_sb_prescribe_po4, &
                                        rtm_decomposition, &
                                        rmm_decomposition, &
                                        km_dip_biochem_po4_act, &
                                        bulk_dens_corr_sl, &
                                        fine_root_carbon_sl, &
                                        po4_solute, &
                                        po4_assoc_fast, &
                                        po4_assoc_slow, &
                                        po4_occluded, &
                                        po4_primary, &
                                        vmax_weath_mineral_act_sl, &
                                        k_adsorpt_po4_act, &
                                        k_desorpt_po4_act, &
                                        sb_pool_residue_phosphorus, & 
                                        weathering_po4, &
                                        occlusion_po4, &
                                        slow_exchange_po4, &
                                        biochem_mineralisation_po4 )

    USE mo_sb_constants,            ONLY: tau_biomineralisation, k_occlude_po4, km_rootc_biochem_po4, km_rootc_weath_po4, &
      &                                   frac_background_weath_mineral

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                           INTENT(in)  :: nc                            !< dimensions
    INTEGER,                           INTENT(in)  :: nsoil_sb                      !< number of soil layers
    LOGICAL,                           INTENT(in)  :: flag_sb_prescribe_po4         !< whether or not solute PO4 is prescribed  
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)  :: rtm_decomposition, &          !< temperature response of decomposition
                                                      rmm_decomposition, &          !< moisture response of decomposition
                                                      km_dip_biochem_po4_act, &     !< actual half-saturation soluble P for biochem_mineralisation_po4 [mol P / m3]
                                                      bulk_dens_corr_sl, &          !< bulk density of the soil layer, [kg/m3]
                                                      fine_root_carbon_sl, &        !< carbon of fine roots per soil layer, [mol/m3]
                                                      po4_solute, &                 !< soluble PO4 pool [mol/m3]
                                                      po4_assoc_fast, &             !< fast minerally associated PO4 pool [mol/m3]
                                                      po4_assoc_slow, &             !< slow minerally associated PO4 pool [mol/m3]
                                                      po4_occluded, &               !< occluded PO4 pool [mol/m3]
                                                      po4_primary, &                !< primary PO4 pool [mol/m3]
                                                      vmax_weath_mineral_act_sl, &  !< actual weathering rate for primary P [mol P kg-1 s-1]
                                                      k_adsorpt_po4_act, &          !< actual adsorption rate for fast_assoc_po4 [s-1]
                                                      k_desorpt_po4_act, &          !< actual desorption rate for slow_assoc_po4 [s-1]
                                                      sb_pool_residue_phosphorus    !< slow SOM P pool [mol / m3] 
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out) :: weathering_po4, &             !< phosphorus weathering rate [mol P / m3 / timestep]
                                                      occlusion_po4, &              !< phosphorus occlusion rate, assoc_slow --> occluded [mol P / m3 / timestep]
                                                      slow_exchange_po4, &          !< phosphorus desorption rate, assoc_fast --> assoc_slow [mol P / m3 / timestep]
                                                      biochem_mineralisation_po4    !< potential biomineralisation rate [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, hlp2, hlp3
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_p_flux_rate'



    !> 1.0 Calculate P weathering rate from Primary P pool (mol /m3 / timestep)
    !! weathering rate is composed of a background rate (frac_background_weath_mineral) and pH/root biomass regulated rate 
    weathering_po4(:,:) = vmax_weath_mineral_act_sl(:,:) * po4_primary &
      &                   * (frac_background_weath_mineral + fine_root_carbon_sl(:,:) &
      &                   / (km_rootc_weath_po4 + fine_root_carbon_sl(:,:))) * dtime                     
    WHERE (weathering_po4(:,:) > po4_primary(:,:)) 
      weathering_po4(:,:) = po4_primary(:,:)
    ENDWHERE
    !> 1.1 Calculate P occlusion rate from occluded PO4 pool (mol /m3 / timestep)
    occlusion_po4(:,:)     = po4_assoc_slow(:,:) * k_occlude_po4 * dtime
    !> 1.2 Calculate P slow adsorption rate (mol /m3 / timestep)
    hlp1(:,:)              = k_desorpt_po4_act(:,:) * po4_assoc_slow(:,:) * dtime
    hlp2(:,:)              = k_adsorpt_po4_act(:,:) * po4_assoc_fast(:,:) * dtime
    slow_exchange_po4(:,:) = hlp2(:,:) - hlp1(:,:)  

    !> 2.0 Calculate potential P biomineralisation from slow organic P pool (mol /m3 / timestep)
    !!
    IF(.NOT. flag_sb_prescribe_po4) THEN
      hlp1(:,:)                       = km_dip_biochem_po4_act(:,:) / (km_dip_biochem_po4_act(:,:) + po4_solute(:,:))
      WHERE (fine_root_carbon_sl(:,:) > eps8)
        hlp2(:,:)                   = fine_root_carbon_sl(:,:) / (km_rootc_biochem_po4 + fine_root_carbon_sl(:,:))
      ELSEWHERE
        hlp2(:,:)                   = 0._wp
      ENDWHERE
      biochem_mineralisation_po4(:,:) = rtm_decomposition(:,:) * rmm_decomposition(:,:) / tau_biomineralisation * & 
        &                               sb_pool_residue_phosphorus(:,:) * hlp1(:,:) * hlp2(:,:) * dtime
    ELSE
      biochem_mineralisation_po4(:,:) = 0.0_wp
    ENDIF 

  END SUBROUTINE calc_potential_p_flux_rate
 
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates asymbiotic N fixation as a function of N deficit during litter decomposition
  !! and limiting functions from soil temperature, moisture and soil nutrient status 
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL FUNCTION calc_asymb_n_fixation(bnf_scheme, &
                                           sb_loss_litter_carbon, &
                                           sb_loss_litter_nitrogen, &
                                           apparent_fast_som_cn, &
                                           rtm_asymb_bnf, &
                                           rmm_asymb_bnf, &
                                           nh4_solute, &
                                           no3_solute) RESULT (asymb_n_fixation) 

    USE mo_sb_constants,          ONLY: frac_litter2fast_som, microbial_nue, km_asymb_bnf, vmax_asymb_bnf

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CHARACTER(len=*), INTENT(in)     :: bnf_scheme                    !< SB_ asymbiotic N fixation scheme  
    REAL(wp), INTENT(in)             :: sb_loss_litter_carbon,  &     !< C flux from litter to fast pool [mol / m3 / timestep]
                                        sb_loss_litter_nitrogen,  &   !< N flux from litter to fast pool [mol / m3 / timestep]
                                        apparent_fast_som_cn, &       !< apparent microbial carbon-nitrogen ratio [mol/mol]
                                        rtm_asymb_bnf, &              !< temperature rate modifier []
                                        rmm_asymb_bnf, &              !< moisture rate modifier []
                                        nh4_solute, &                 !< NH4 in solution [mol / m3] 
                                        no3_solute                    !< NO3 in solution [mol / m3]
    REAL(wp)                         :: asymb_n_fixation              !< asymbiotic N fixation [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                         :: f_nav,n_demand 
    CHARACTER(len=*), PARAMETER      :: routine = TRIM(modname)//':calc_asymb_n_fixation'
    
    SELECT CASE(TRIM(bnf_scheme))
    ! for the SB_ asymbiotic N fixation no difference is made between "dynamic" and "fixed" 
    !  ("fixed" is only available in the VEG_ symbiotic N fixation)
    CASE("dynamic", "fixed")
      !>1.0 N deficit of litter decomposition 
      !>
      
      !>1.1 calculate potential transfer from litter decay to fast SOM pool, and associated potential microbial
      !>immobilisation rate as difference between maximal growth rate and N in litter
      !>
      n_demand = frac_litter2fast_som / apparent_fast_som_cn * sb_loss_litter_carbon - & 
        &        microbial_nue * sb_loss_litter_nitrogen

      !>1.2 rate modifier for mineral soil nitrogen 
      !>
      IF(nh4_solute + no3_solute > eps8) THEN
        f_nav = MIN(1._wp, 1.0_wp - (nh4_solute + no3_solute) / (nh4_solute + no3_solute + km_asymb_bnf))
      ELSE
        f_nav = 1._wp 
      ENDIF

      asymb_n_fixation = vmax_asymb_bnf * rtm_asymb_bnf * rmm_asymb_bnf * f_nav * n_demand
    CASE("none", "unlimited")
      !>1.3 asymb_n_fixation is calculated during the decomposition process for these cases
      !>
      asymb_n_fixation = 0.0_wp
    END SELECT

  END FUNCTION calc_asymb_n_fixation 

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the immobilisation rate of NH4 and PO4 of microbial biomass, given
  !! litter decay and stoichiometry, the co-occurring rates of plant, mycorrhizal uptake,
  !! nitrification, denitrification and transport, biomineralisation of PO4
  !! such that the newly formed fast SOM has the prescribed C:N:P stoichiometry; the carbon-
  !! use-efficiency of microbes is within the observed bound and mass-balance constraints are 
  !! maintained. The result of this calculation are updates ("actual") rates of uptake, and processing
  !! as well as litter turnover rates.  
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_actual_flux_rates( nc, &                         ! in
                                     nsoil_sb, &
                                     bnf_scheme, & 
                                     nh4_solute, &
                                     no3_solute, &
                                     po4_solute, &
                                     po4_assoc_fast, &
                                     transport_nh4_solute, &
                                     transport_no3_solute, &
                                     transport_po4_solute, &
                                     transport_po4_assoc_fast, &
                                     weathering_po4, &
                                     p_deposition_act, &
                                     slow_exchange_po4, &
                                     apparent_fast_som_cn , &
                                     km_mic_uptake_nh4_act, & 
                                     km_mic_uptake_no3_act, & 
                                     plant_uptake_nh4_sl, &                 ! inout
                                     plant_uptake_no3_sl, &
                                     plant_uptake_po4_sl, &
                                     mycorrhiza_uptake_nh4_sl, &
                                     mycorrhiza_uptake_no3_sl, &
                                     mycorrhiza_uptake_po4_sl, &
                                     nitrification_no3, &
                                     denitrification_n2, &
                                     asymb_n_fixation, &
                                     biochem_mineralisation_po4, &
                                     sb_loss_mt, &                          ! inout
                                     microbial_uptake_nh4_sl, &             ! out
                                     microbial_uptake_no3_sl, &             ! out
                                     microbial_uptake_po4_sl)               ! out 

    USE mo_sb_constants,            ONLY: frac_litter2fast_som, k_fast_som_np, microbial_nue, microbial_pue, &
                                          vmax_mic_uptake_nh4, vmax_mic_uptake_no3

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                           INTENT(in)    :: nc, &                         !< dimensions
                                                        nsoil_sb
    CHARACTER(len=*),                  INTENT(in)    :: bnf_scheme                    !< SB_ asymbiotic N fixation scheme 
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_solute, &                 !< NH4 in solution [mol / m3] 
                                                        no3_solute, &                 !< NO3 in solution [mol / m3]
                                                        po4_solute, &                 !< PO4 in solution [mol / m3]
                                                        po4_assoc_fast, &             !< labile mineral P [mol / m3]
                                                        transport_nh4_solute, &       !< transport of NH4 [mol / m3 / timestep]
                                                        transport_no3_solute, &       !< transport of NO3 [mol / m3 / timestep]
                                                        transport_po4_solute, &       !< transport of soluble PO4 [mol / m3 / timestep]
                                                        transport_po4_assoc_fast, &   !< transport of labile PO4 [mol / m3 / timestep]
                                                        weathering_po4, &             !< weathering of po4 [mol / m3 / timestep]
                                                        slow_exchange_po4, &          !< po4 desorption rate, assoc_fast --> assoc_slow [mol P / m3 / timestep]
                                                        p_deposition_act, &           !< actual P deposition rate, [mol / m3 / timestep]
                                                        apparent_fast_som_cn, &       !< apparent microbial carbon-nitrogen ratio [mol/mol]
                                                        km_mic_uptake_nh4_act, &      !< actual km for microbial NH4 uptake [mol / m3]
                                                        km_mic_uptake_no3_act         !< actual km for microbial NO3 uptake [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: plant_uptake_nh4_sl, &        !< plant NH4 uptake [mol / m3 / timestep]
                                                        plant_uptake_no3_sl, &        !< plant NO3 uptake [mol / m3 / timestep]
                                                        plant_uptake_po4_sl, &        !< plant PO4 uptake [mol / m3 / timestep]
                                                        mycorrhiza_uptake_nh4_sl, &   !< mycorrhizal NH4 uptake [mol / m3 / timestep]
                                                        mycorrhiza_uptake_no3_sl, &   !< mycorrhizal NO3 uptake [mol / m3 / timestep]
                                                        mycorrhiza_uptake_po4_sl, &   !< mycorrhizal PO4 uptake [mol / m3 / timestep]
                                                        nitrification_no3, &          !< nitrification rate [mol / m3 / timestep]
                                                        denitrification_n2, &         !< denitrification rate [mol / m3 / timestep]
                                                        asymb_n_fixation, &           !< asymbiotic N fixation (unlimited BNF only) [mol / m3 / timestep]
                                                        biochem_mineralisation_po4    !< biomineralisation rate [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_loss_mt(:,:,:,:)           !< bgcm sb_loss flux: decomposition of SOM and litter pools [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: microbial_uptake_nh4_sl, &    !< microbial immobilisation of NH4 [mol / m3 / timestep]
                                                        microbial_uptake_no3_sl, &    !< microbial immobilisation of NO3 [mol / m3 / timestep]
                                                        microbial_uptake_po4_sl       !< microbial immobilisation of PO4 [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(nc, nsoil_sb) :: hlp1, hlp2
    REAL(wp), DIMENSION(nc, nsoil_sb) :: transfer_c
    REAL(wp), DIMENSION(nc, nsoil_sb) :: f_uptake_nh4_sl
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_actual_flux_rates'


    !> 1.0 formation of fast SOM from litter 
    !!

    !> 1.1 calculate potential transfer from litter decay to fast SOM pool, and associated potential microbial 
    !!     immobilisation rate as difference between maximal growth rate and N in litter
    transfer_c(:,:) = frac_litter2fast_som * & 
                      (sb_loss_mt(ix_soluable_litter, ixC, :, :) + & 
                       sb_loss_mt(ix_polymeric_litter, ixC, :, :)) 
   
    !> 1.2 calculate partitioning of N immobilisation to NH4 and NO3 according to their abundance and the 
    !!     affinities of the uptake channels to either species

    WHERE(nh4_solute(:,:) + no3_solute(:,:) > eps8)
      f_uptake_nh4_sl(:,:) = vmax_mic_uptake_nh4 * nh4_solute(:,:) / (nh4_solute(:,:) + km_mic_uptake_nh4_act) / &
                             (vmax_mic_uptake_nh4 * nh4_solute(:,:) / (nh4_solute(:,:) + km_mic_uptake_nh4_act) + &
                              vmax_mic_uptake_no3 * no3_solute(:,:) / (no3_solute(:,:) + km_mic_uptake_no3_act))
    ELSEWHERE
      f_uptake_nh4_sl(:,:) = 1.0_wp
    ENDWHERE
 
    microbial_uptake_nh4_sl(:,:) = f_uptake_nh4_sl(:,:) * & 
                                             (transfer_c(:,:) / apparent_fast_som_cn(:,:) - & 
                                              microbial_nue * (sb_loss_mt(ix_soluable_litter, ixN, :, :) + & 
                                                              sb_loss_mt(ix_polymeric_litter, ixN, :, :)) - &
                                              asymb_n_fixation(:,:))

    microbial_uptake_no3_sl(:,:) = (1.0_wp - f_uptake_nh4_sl(:,:)) * & 
                                             (transfer_c(:,:) / apparent_fast_som_cn(:,:) - & 
                                              microbial_nue * (sb_loss_mt(ix_soluable_litter, ixN, :, :) + & 
                                                              sb_loss_mt(ix_polymeric_litter, ixN, :, :)) - &
                                              asymb_n_fixation(:,:))

    !! limit uptake rate at low nitrogen assuming a Michaelis-Menten Kinetic 
    hlp1(:,:)                    = nh4_solute(:,:) / ( km_mic_uptake_nh4_act + nh4_solute(:,:) )
    microbial_uptake_nh4_sl(:,:) = microbial_uptake_nh4_sl(:,:) * MIN(hlp1(:,:), 1.0_wp)

    hlp1(:,:)                    = no3_solute(:,:) / ( km_mic_uptake_no3_act + no3_solute(:,:) )
    microbial_uptake_no3_sl(:,:) = microbial_uptake_no3_sl(:,:) * MIN(hlp1(:,:), 1.0_wp)

    !>1.3  Limit microbial, mycorrhizal and plant uptake, as well as nitrification to available NH4 minus transport losses
    !>
    SELECT CASE(TRIM(bnf_scheme))
    CASE("unlimited")
      asymb_n_fixation(:,:)        = MAX(microbial_uptake_nh4_sl(:,:) + plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + & 
                                         mycorrhiza_uptake_nh4_sl(:,:) - (nh4_solute(:,:) + transport_nh4_solute(:,:)), 0.0_wp)
    ! for the SB_ asymbiotic N fixation no difference is made between "dynamic" and "fixed" 
    !  ("fixed" is only available in the VEG_ symbiotic N fixation)
    CASE("none", "dynamic", "fixed")
      WHERE(microbial_uptake_nh4_sl(:,:) >= 0.0_wp)
        hlp1(:,:) = microbial_uptake_nh4_sl(:,:) + plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + &
                    mycorrhiza_uptake_nh4_sl(:,:)
        hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:)
        WHERE(hlp2(:,:) < hlp1(:,:)) 
          WHERE(hlp1(:,:) > eps8)
            microbial_uptake_nh4_sl(:,:)  = hlp2(:,:) * microbial_uptake_nh4_sl(:,:) / hlp1(:,:) 
            plant_uptake_nh4_sl(:,:)      = hlp2(:,:) * plant_uptake_nh4_sl(:,:) / hlp1(:,:)
            nitrification_no3(:,:)        = hlp2(:,:) * nitrification_no3(:,:) / hlp1(:,:)
            mycorrhiza_uptake_nh4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_nh4_sl(:,:) / hlp1(:,:)
          ELSEWHERE
            microbial_uptake_nh4_sl(:,:)  = 0.0_wp
            plant_uptake_nh4_sl(:,:)      = 0.0_wp 
            nitrification_no3(:,:)        = 0.0_wp 
            mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
          ENDWHERE 
        ENDWHERE
      ELSEWHERE
        hlp1(:,:) = plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + mycorrhiza_uptake_nh4_sl(:,:)
        hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:) - microbial_uptake_nh4_sl(:,:)
        WHERE(hlp2(:,:) < hlp1(:,:)) 
          WHERE(hlp1(:,:) > eps8)
            plant_uptake_nh4_sl(:,:)      = hlp2(:,:) * plant_uptake_nh4_sl(:,:) / hlp1(:,:)
            nitrification_no3(:,:)        = hlp2(:,:) * nitrification_no3(:,:) / hlp1(:,:)
            mycorrhiza_uptake_nh4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_nh4_sl(:,:) / hlp1(:,:)
          ELSEWHERE
            plant_uptake_nh4_sl(:,:)      = 0.0_wp 
            nitrification_no3(:,:)        = 0.0_wp 
            mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
          ENDWHERE 
        ENDWHERE
      ENDWHERE 
    CASE DEFAULT
      CALL finish(TRIM(routine),'No such bnf_scheme available here.')
    END SELECT
 
    !> 1.4  limit plant uptake and denitrification to available NO3 minus transport losses 
    !!
    WHERE(microbial_uptake_no3_sl(:,:) >= 0.0_wp)
      hlp1(:,:) = microbial_uptake_no3_sl(:,:) + plant_uptake_no3_sl(:,:) + denitrification_n2(:,:) + mycorrhiza_uptake_no3_sl(:,:)
      hlp2(:,:) = no3_solute(:,:) + transport_no3_solute(:,:)
      WHERE(hlp2(:,:) < hlp1(:,:))
        WHERE(hlp1(:,:) > eps8)
          microbial_uptake_no3_sl(:,:)  = hlp2(:,:) * microbial_uptake_no3_sl(:,:) / hlp1(:,:) 
          plant_uptake_no3_sl(:,:)      = hlp2(:,:) * plant_uptake_no3_sl(:,:) / hlp1(:,:)
          denitrification_n2(:,:)       = hlp2(:,:) * denitrification_n2(:,:) / hlp1(:,:)
          mycorrhiza_uptake_no3_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_no3_sl(:,:) / hlp1(:,:)
        ELSEWHERE
          microbial_uptake_no3_sl(:,:)  = 0.0_wp
          plant_uptake_no3_sl(:,:)      = 0.0_wp 
          denitrification_n2(:,:)       = 0.0_wp 
          mycorrhiza_uptake_no3_sl(:,:) = 0.0_wp
        ENDWHERE
      ENDWHERE 
    ELSEWHERE
      hlp1(:,:) = plant_uptake_no3_sl(:,:) + denitrification_n2(:,:) + mycorrhiza_uptake_no3_sl(:,:)
      hlp2(:,:) = no3_solute(:,:) + transport_no3_solute(:,:) - microbial_uptake_no3_sl(:,:)
      WHERE(hlp2(:,:) < hlp1(:,:))
        WHERE(hlp1(:,:) > eps8)
          plant_uptake_no3_sl(:,:)      = hlp2(:,:) * plant_uptake_no3_sl(:,:) / hlp1(:,:)
          denitrification_n2(:,:)       = hlp2(:,:) * denitrification_n2(:,:) / hlp1(:,:)
          mycorrhiza_uptake_no3_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_no3_sl(:,:) / hlp1(:,:)
        ELSEWHERE
          plant_uptake_no3_sl(:,:)      = 0.0_wp 
          denitrification_n2(:,:)       = 0.0_wp 
          mycorrhiza_uptake_no3_sl(:,:) = 0.0_wp
        ENDWHERE
      ENDWHERE
    ENDWHERE
           
    !> 2.0 determine actual growth rate of fast SOM pool
    !!
   
    !> 2.1 Calculate actual litter loss rate, calculated such that litter decomposition
    !! matches the actual immobilisation flux and CN uptake from litter given actual CUE
    hlp1(:,:) = frac_litter2fast_som / apparent_fast_som_cn(:,:) &
      &         * (sb_loss_mt(ix_soluable_litter, ixC, :, :) + sb_loss_mt(ix_polymeric_litter, ixC, :, :)) &
      &         - (asymb_n_fixation(:,:) + microbial_nue &
      &         * (sb_loss_mt(ix_soluable_litter, ixN, :, :) + sb_loss_mt(ix_polymeric_litter, ixN, :, :))) 
    hlp2(:,:) = 1.0_wp
    WHERE((microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) < hlp1(:,:))
      WHERE(hlp1(:,:) > eps8)
        hlp2(:,:) = (microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:)) / hlp1(:,:)
      ELSEWHERE
        hlp2(:,:) = 0.0_wp
      ENDWHERE
      sb_loss_mt(ix_soluable_litter, ixC, :, :)    = hlp2(:,:) * sb_loss_mt(ix_soluable_litter, ixC, :, :)
      sb_loss_mt(ix_soluable_litter, ixN, :, :)    = hlp2(:,:) * sb_loss_mt(ix_soluable_litter, ixN, :, :)
      sb_loss_mt(ix_soluable_litter, ixP, :, :)    = hlp2(:,:) * sb_loss_mt(ix_soluable_litter, ixP, :, :)
      sb_loss_mt(ix_soluable_litter, ixC13, :, :)  = hlp2(:,:) * sb_loss_mt(ix_soluable_litter, ixC13, :, :)
      sb_loss_mt(ix_soluable_litter, ixC14, :, :)  = hlp2(:,:) * sb_loss_mt(ix_soluable_litter, ixC14, :, :)
      sb_loss_mt(ix_soluable_litter, ixN15, :, :)  = hlp2(:,:) * sb_loss_mt(ix_soluable_litter, ixN15, :, :)

      sb_loss_mt(ix_polymeric_litter, ixC, :, :)   = hlp2(:,:) * sb_loss_mt(ix_polymeric_litter, ixC, :, :)
      sb_loss_mt(ix_polymeric_litter, ixN, :, :)   = hlp2(:,:) * sb_loss_mt(ix_polymeric_litter, ixN, :, :)
      sb_loss_mt(ix_polymeric_litter, ixP, :, :)   = hlp2(:,:) * sb_loss_mt(ix_polymeric_litter, ixP, :, :)
      sb_loss_mt(ix_polymeric_litter, ixC13, :, :) = hlp2(:,:) * sb_loss_mt(ix_polymeric_litter, ixC13, :, :)
      sb_loss_mt(ix_polymeric_litter, ixC14, :, :) = hlp2(:,:) * sb_loss_mt(ix_polymeric_litter, ixC14, :, :)
      sb_loss_mt(ix_polymeric_litter, ixN15, :, :) = hlp2(:,:) * sb_loss_mt(ix_polymeric_litter, ixN15, :, :)
    ENDWHERE 

    !>3.0 calculate potential microbial P uptake from soluable pool (mol /m3 / timestep)
    !>
    microbial_uptake_po4_sl(:,:) = frac_litter2fast_som / apparent_fast_som_cn(:,:) / k_fast_som_np &
      &                            * (sb_loss_mt(ix_soluable_litter, ixC, :, :) + sb_loss_mt(ix_polymeric_litter, ixC, :, :)) &
      &                            - microbial_pue &
      &                            * (sb_loss_mt(ix_soluable_litter, ixP, :, :) + sb_loss_mt(ix_polymeric_litter, ixP, :, :))

    !>  3.1 constrain P uptake by microbes and plants to the available sum of P from P_soluable, 
    !>      biomineralisation and transport (mol /m3 / timestep)
    !>
    WHERE(microbial_uptake_po4_sl(:,:) > 0.0_wp) ! if microbes take up P, constrain both microbial and plant uptake to available P 
      hlp1(:,:) = microbial_uptake_po4_sl(:,:) + plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:)
      hlp2(:,:) = po4_solute(:,:) + transport_po4_solute(:,:) 
      WHERE(hlp2(:,:) < hlp1(:,:))
        WHERE(hlp1(:,:) > eps8)
          microbial_uptake_po4_sl(:,:)  = hlp2(:,:) * microbial_uptake_po4_sl(:,:) / hlp1(:,:) 
          plant_uptake_po4_sl(:,:)      = hlp2(:,:) * plant_uptake_po4_sl(:,:) / hlp1(:,:) 
          mycorrhiza_uptake_po4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_po4_sl(:,:) / hlp1(:,:) 
        ELSEWHERE
          microbial_uptake_po4_sl(:,:)  = 0.0_wp
          plant_uptake_po4_sl(:,:)      = 0.0_wp
          mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp 
        ENDWHERE
      ENDWHERE
    ELSEWHERE ! if microbes release P, constrain plant uptake to available P, note the negative sign of microbial_uptake_po4_sl(:,:)
      hlp1(:,:) = plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) 
      hlp2(:,:) = po4_solute(:,:) + transport_po4_solute(:,:) - microbial_uptake_po4_sl(:,:)
      WHERE(hlp2(:,:) < hlp1(:,:))
        WHERE(hlp1(:,:) > eps8)
          plant_uptake_po4_sl(:,:)      = hlp2(:,:) * plant_uptake_po4_sl(:,:) / hlp1(:,:) 
          mycorrhiza_uptake_po4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_po4_sl(:,:) / hlp1(:,:) 
        ELSEWHERE
          plant_uptake_po4_sl(:,:)      = 0.0_wp 
          mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp 
        ENDWHERE
      ENDWHERE
    ENDWHERE
  
  END SUBROUTINE calc_actual_flux_rates
  
  ! ======================================================================================================= !
  !>For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to 
  !>  net mineralisation \n
  !>
  !>For the sb_nloss_scheme = dynamic: calculates the partitioning of nitrification and denitrification
  !>  to NO3, NOy, N2O and N2, and diagnoses the amount of leaching below the rooting zone
  !>
  SUBROUTINE calc_partitioning_nitrification_denitrification_rate( &
                                    nc, &                         ! in
                                    nsoil_sb, &
                                    rtm_denitrification, &
                                    nitrification_no3, &
                                    denitrification_n2, &
                                    nitrification_noy, &
                                    nitrification_n2o, &
                                    denitrification_noy, &
                                    denitrification_n2o)

    USE mo_sb_constants,            ONLY: f_nit_noy, f_nit_n2o, f_denit_noy, f_denit_n2o
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                      !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: rtm_denitrification     !< temperature modifier for denitrification
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: nitrification_no3, &    !< nitrification to NO3 [mol / m3 / timestep]
                                                        denitrification_n2      !< denitrification to N2 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(out)   :: nitrification_noy, &    !< nitrification to NOy [mol / m3 / timestep]
                                                        nitrification_n2o, &    !< nitrification to N2O [mol / m3 / timestep]
                                                        denitrification_noy, &  !< denitrification to NOy [mol / m3 / timestep]
                                                        denitrification_n2o     !< denitrification to N2O [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_partitioning_nitrification_denitrification_rate'
    ! ----------------------------------------------------------------------------------------------------- !
    
    !>1.0 case of dynamic N loss calculations, following Zaehle et al. 2011
    !>
    ! nitrification rate to produce NO3, NOy and N2O from NH4 [mol m-3 timestep-1]
    nitrification_noy(:,:)   = f_nit_noy * nitrification_no3(:,:)
    nitrification_n2o(:,:)   = f_nit_n2o * nitrification_no3(:,:)
    nitrification_no3(:,:)   = nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:)
    ! denitrification rate to produce NOy, N2O and N2 from NO3 [mol m-3 timestep-1]
    denitrification_noy(:,:) = f_denit_noy * rtm_denitrification(:,:) * denitrification_n2(:,:)
    denitrification_n2o(:,:) = f_denit_n2o * rtm_denitrification(:,:) * denitrification_n2(:,:)
    denitrification_n2(:,:)  = denitrification_n2(:,:) - denitrification_noy(:,:) - denitrification_n2o(:,:)
  END SUBROUTINE calc_partitioning_nitrification_denitrification_rate

  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to both: update_sb_simple_model and update_sb_jsm
  !   units for flux rates differ between the two models: \n
  !     in this routine they are assumed to be in [mol/m3/timestep] \n
  !     simple soil model: [mol/m3/timestep]
  !     JSM:               [micro-mol/m3/s] (unit conversion is applied within the routine)
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the N15 fluxes based on uptake or transformation of the inorganic NH4 / NO3 pools
  !! according to the framework proposed by Robinson, 2001, TREE. Includes a representation of the 
  !! Rayleigh model to account for the change in fractionation with increasing consumption  
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sb_inorganic_n15_fluxes(nc, & 
                                          nsoil_sb, &
                                          sb_model_scheme, &
                                          nh4_solute, &
                                          nh4_n15_solute, &
                                          no3_solute, &
                                          no3_n15_solute, &
                                          plant_uptake_nh4_sl, & 
                                          plant_uptake_no3_sl, &      
                                          mycorrhiza_uptake_nh4_sl, & 
                                          mycorrhiza_uptake_no3_sl, & 
                                          microbial_uptake_nh4_sl, &        
                                          microbial_uptake_no3_sl, &        
                                          asymb_n_fixation, &
                                          nitrification_no3, &        
                                          nitrification_noy, &        
                                          nitrification_n2o, &        
                                          denitrification_noy, &      
                                          denitrification_n2o, &      
                                          denitrification_n2, &       
                                          transport_nh4_n15_solute, & 
                                          transport_no3_n15_solute, &    
                                          plant_uptake_nh4_n15_sl, &      
                                          plant_uptake_no3_n15_sl, &      
                                          mycorrhiza_uptake_nh4_n15_sl, & 
                                          mycorrhiza_uptake_no3_n15_sl, & 
                                          microbial_uptake_nh4_n15_sl, &        
                                          microbial_uptake_no3_n15_sl, &        
                                          asymb_n_fixation_n15, &
                                          nitrification_no3_n15, &        
                                          nitrification_noy_n15, &        
                                          nitrification_n2o_n15, &        
                                          denitrification_noy_n15, &      
                                          denitrification_n2o_n15, &      
                                          denitrification_n2_n15) 
 
    USE mo_isotope_util,   ONLY: calc_fractionation, calc_mixing_ratio_N15N14 
    USE mo_sb_constants,   ONLY: eta_plant_uptake_nh4, eta_plant_uptake_no3, eta_mic_uptake_nh4, eta_mic_uptake_no3, & 
      &                          eta_nitrate_production, eta_nitrification, eta_denitrification
    USE mo_veg_constants,  ONLY: eta_nfixation
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                            INTENT(in)    :: nc, &                       !< dimensions
                                                         nsoil_sb
    CHARACTER(len=*),                   INTENT(in)    :: sb_model_scheme             !< simple_1d / jsm from sb config
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(in)    :: nh4_solute, &               !< soil solute NH4 concentration [mol m-3]
                                                         nh4_n15_solute, &           !< soil solute 15NH4 concentration [mol m-3]
                                                         no3_solute, &               !< soil solute NO3 concentration [mol m-3]
                                                         no3_n15_solute, &           !< soil solute 15NO3 concentration [mol m-3]
                                                         plant_uptake_nh4_sl, &      !< NH4 uptake by plants [mol / m3 / timestep]
                                                         plant_uptake_no3_sl, &      !< NO3 uptake by plants [mol / m3 / timestep]
                                                         mycorrhiza_uptake_nh4_sl, & !< NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                         mycorrhiza_uptake_no3_sl, & !< NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                         microbial_uptake_nh4_sl, &  !< NH4 uptake by microbes [mol / m3 / timestep]
                                                         microbial_uptake_no3_sl, &  !< NO3 uptake by microbes [mol / m3 / timestep]
                                                         asymb_n_fixation, &         !< asymbiotic N fixation [mol / m3 / timestep]
                                                         nitrification_no3, &        !< nitrification to NO3 [mol / m3 / timestep]
                                                         nitrification_noy, &        !< nitrification to NOy [mol / m3 / timestep]
                                                         nitrification_n2o, &        !< nitrification to N2O [mol / m3 / timestep]
                                                         denitrification_noy, &      !< denitrification to NOy [mol / m3 / timestep]
                                                         denitrification_n2o, &      !< denitrification to N2O [mol / m3 / timestep]
                                                         denitrification_n2, &       !< denitrification to N2 [mol / m3 / timestep]
                                                         transport_nh4_n15_solute, & !< transport of 15NH4 [mol / m3 / timestep]
                                                         transport_no3_n15_solute    !< transport of 15NO3 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),  INTENT(out)   :: plant_uptake_nh4_n15_sl, &      !< 15NH4 uptake by plants [mol / m3 / timestep]
                                                         plant_uptake_no3_n15_sl, &      !< 15NO3 uptake by plants [mol / m3 / timestep]
                                                         mycorrhiza_uptake_nh4_n15_sl, & !< 15NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                         mycorrhiza_uptake_no3_n15_sl, & !< 15NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                         microbial_uptake_nh4_n15_sl, &  !< 15NH4 uptake by microbes [mol / m3 / timestep]
                                                         microbial_uptake_no3_n15_sl, &  !< 15NH4 uptake by microbes [mol / m3 / timestep]
                                                         asymb_n_fixation_n15, &         !< asymbiotic 15N fixation [mol / m3 / timestep]
                                                         nitrification_no3_n15, &        !< nitrification to 15NO3 [mol / m3 / timestep]
                                                         nitrification_noy_n15, &        !< nitrification to 15NOy [mol / m3 / timestep]
                                                         nitrification_n2o_n15, &        !< nitrification to 15N2O [mol / m3 / timestep]
                                                         denitrification_noy_n15, &      !< denitrification to 15NOy [mol / m3 / timestep]
                                                         denitrification_n2o_n15, &      !< denitrification to 15N2O [mol / m3 / timestep]
                                                         denitrification_n2_n15          !< denitrification to 15N2 [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_sb)         :: hlp1, hlp2 
    REAL(wp), DIMENSION(nc, nsoil_sb)         :: &                ! flux rates with unit conversion for JSM
                                                 plant_uptake_nh4_local_sl, &
                                                 plant_uptake_no3_local_sl, &
                                                 mycorrhiza_uptake_nh4_local_sl, &
                                                 mycorrhiza_uptake_no3_local_sl, &
                                                 microbial_uptake_nh4_local_sl, &
                                                 microbial_uptake_no3_local_sl, &
                                                 asymb_n_fixation_local, &
                                                 nitrification_no3_local, &
                                                 nitrification_noy_local, &
                                                 nitrification_n2o_local, &
                                                 denitrification_noy_local, &
                                                 denitrification_n2o_local, &
                                                 denitrification_n2_local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sb_inorganic_n15_fluxes'
    ! ----------------------------------------------------------------------------------------------------- !
    
    
    !>0.9 'INTENT(in) valiables': flux rates unit conversion when this routines is called from JSM
    !>  unit conversion: [micro-mol/m3/s] -> [mol/m3/timestep]
    !>
    !> @NOTE @TODO ? not all INTENT(in) variables are considered ?
    IF(sb_model_scheme == 'jsm') THEN
      plant_uptake_nh4_local_sl(:,:)      = plant_uptake_nh4_sl(:,:)      * dtime / 1.e6_wp
      plant_uptake_no3_local_sl(:,:)      = plant_uptake_no3_sl(:,:)      * dtime / 1.e6_wp
      mycorrhiza_uptake_nh4_local_sl(:,:) = mycorrhiza_uptake_nh4_sl(:,:) * dtime / 1.e6_wp
      mycorrhiza_uptake_no3_local_sl(:,:) = mycorrhiza_uptake_no3_sl(:,:) * dtime / 1.e6_wp
      microbial_uptake_nh4_local_sl(:,:)  = microbial_uptake_nh4_sl(:,:)  * dtime / 1.e6_wp
      microbial_uptake_no3_local_sl(:,:)  = microbial_uptake_no3_sl(:,:)  * dtime / 1.e6_wp
      asymb_n_fixation_local(:,:)         = asymb_n_fixation(:,:)         * dtime / 1.e6_wp
      nitrification_no3_local(:,:)        = nitrification_no3(:,:)        * dtime / 1.e6_wp
      nitrification_noy_local(:,:)        = nitrification_noy(:,:)        * dtime / 1.e6_wp
      nitrification_n2o_local(:,:)        = nitrification_n2o(:,:)        * dtime / 1.e6_wp
      denitrification_noy_local(:,:)      = denitrification_noy(:,:)      * dtime / 1.e6_wp
      denitrification_n2o_local(:,:)      = denitrification_n2o(:,:)      * dtime / 1.e6_wp
      denitrification_n2_local(:,:)       = denitrification_n2(:,:)       * dtime / 1.e6_wp
    ELSE
      plant_uptake_nh4_local_sl(:,:)      = plant_uptake_nh4_sl(:,:)     
      plant_uptake_no3_local_sl(:,:)      = plant_uptake_no3_sl(:,:)     
      mycorrhiza_uptake_nh4_local_sl(:,:) = mycorrhiza_uptake_nh4_sl(:,:)
      mycorrhiza_uptake_no3_local_sl(:,:) = mycorrhiza_uptake_no3_sl(:,:)
      microbial_uptake_nh4_local_sl(:,:)  = microbial_uptake_nh4_sl(:,:) 
      microbial_uptake_no3_local_sl(:,:)  = microbial_uptake_no3_sl(:,:) 
      asymb_n_fixation_local(:,:)         = asymb_n_fixation(:,:) 
      nitrification_no3_local(:,:)        = nitrification_no3(:,:)       
      nitrification_noy_local(:,:)        = nitrification_noy(:,:)       
      nitrification_n2o_local(:,:)        = nitrification_n2o(:,:)       
      denitrification_noy_local(:,:)      = denitrification_noy(:,:)     
      denitrification_n2o_local(:,:)      = denitrification_n2o(:,:)     
      denitrification_n2_local(:,:)       = denitrification_n2(:,:)      
    ENDIF

    !> 1.0 calculate the actual discrimination for NH4, 
    !!     which depends on the fraction of the source consumed f = (sum of fluxes)/stock = hlp1. 
    !!     The discrimination is then reduced by a reduction factor r = - (1-f) ln (1-f) / f = hlp2
    ! fraction of NH4 that is consumed
    WHERE(nh4_solute(:,:) > eps8)
      hlp1(:,:) = MAX((plant_uptake_nh4_local_sl(:,:) + mycorrhiza_uptake_nh4_local_sl(:,:) + &
        &              microbial_uptake_nh4_local_sl(:,:) + nitrification_no3_local(:,:) + &
        &              nitrification_noy_local(:,:) + nitrification_n2o_local(:,:)) / nh4_solute(:,:), eps8)
    ELSEWHERE
      hlp1(:,:) = 1.0_wp
    ENDWHERE

    ! calculate reduction factor for discrimination 
    WHERE(hlp1(:,:) < (1.0_wp - eps8))
      hlp2(:,:) = (hlp1(:,:) - 1.0_wp) * LOG(1.0_wp - hlp1(:,:)) / hlp1(:,:)
    ELSEWHERE
      hlp2(:,:) = 0.0_wp
    ENDWHERE

    !> 1.1 15N fluxes involving ammonium
    !!
    plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_local_sl(:,:) * & 
      &                                 calc_fractionation(nh4_solute(:,:), nh4_n15_solute(:,:), &
      &                                                    eta_plant_uptake_nh4 * hlp2(:,:)) 
    mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_local_sl(:,:) * & 
      &                                 calc_fractionation(nh4_solute(:,:), nh4_n15_solute(:,:), &
      &                                                    eta_plant_uptake_nh4 * hlp2(:,:)) 
    microbial_uptake_nh4_n15_sl(:,:)  = microbial_uptake_nh4_local_sl(:,:) * & 
      &                                 calc_fractionation(nh4_solute(:,:), nh4_n15_solute(:,:), &
      &                                                    eta_mic_uptake_nh4 * hlp2(:,:))
    nitrification_no3_n15(:,:)        = nitrification_no3_local(:,:) * &
      &                                 calc_fractionation(nh4_solute(:,:), nh4_n15_solute(:,:), &
      &                                                    eta_nitrate_production * hlp2(:,:))
    nitrification_noy_n15(:,:)        = nitrification_noy_local(:,:) * & 
      &                                 calc_fractionation(nh4_solute(:,:), nh4_n15_solute(:,:), &
      &                                                    eta_nitrification * hlp2(:,:))
    nitrification_n2o_n15(:,:)        = nitrification_n2o_local(:,:) * &
      &                                 calc_fractionation(nh4_solute(:,:), nh4_n15_solute(:,:), &
      &                                                    eta_nitrification * hlp2(:,:)) 

    !> 1.2 contrain 15NH4 fluxes to available 15NH4
    !!
    hlp1(:,:) = plant_uptake_nh4_n15_sl(:,:) + mycorrhiza_uptake_nh4_n15_sl(:,:) + microbial_uptake_nh4_n15_sl(:,:) + & 
                nitrification_no3_n15(:,:) + nitrification_noy_n15(:,:) + nitrification_n2o_n15(:,:)
    hlp2(:,:) = nh4_n15_solute(:,:) + transport_nh4_n15_solute(:,:) 
    WHERE(hlp2(:,:) < hlp1(:,:))
       WHERE(hlp1(:,:) > eps8)
          plant_uptake_nh4_n15_sl(:,:)      = hlp2(:,:) * plant_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
          mycorrhiza_uptake_nh4_n15_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
          microbial_uptake_nh4_n15_sl(:,:)  = hlp2(:,:) * microbial_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
          nitrification_no3_n15(:,:)        = hlp2(:,:) * nitrification_no3_n15(:,:) / hlp1(:,:)
          nitrification_noy_n15(:,:)        = hlp2(:,:) * nitrification_noy_n15(:,:) / hlp1(:,:)
          nitrification_n2o_n15(:,:)        = hlp2(:,:) * nitrification_n2o_n15(:,:) / hlp1(:,:)
       ELSEWHERE
          plant_uptake_nh4_n15_sl(:,:)      = 0.0_wp 
          mycorrhiza_uptake_nh4_n15_sl(:,:) = 0.0_wp 
          microbial_uptake_nh4_n15_sl(:,:)  = 0.0_wp
          nitrification_no3_n15(:,:)        = 0.0_wp
          nitrification_noy_n15(:,:)        = 0.0_wp
          nitrification_n2o_n15(:,:)        = 0.0_wp
       ENDWHERE
    ENDWHERE

    !> 2.0 calculate the actual discrimination for NO3, 
    !!     which depends on the fraction of the source consumed f = (sum of fluxes)/stock = hlp1. 
    !!     The discrimination is then reduced by a reduction factor r = - (1-f) ln (1-f) / f = hlp2
    ! fraction of NO3 that is consumed
    WHERE(no3_solute(:,:) > eps8)
       hlp1(:,:) = MAX((plant_uptake_no3_local_sl(:,:) + mycorrhiza_uptake_no3_local_sl(:,:) + &
                        microbial_uptake_no3_local_sl(:,:) + denitrification_noy_local(:,:) + &
                        denitrification_n2o_local(:,:) + denitrification_n2_local(:,:)) / no3_solute(:,:), eps8)
    ELSEWHERE
       hlp1(:,:) = 1.0_wp
    ENDWHERE
    WHERE(hlp1(:,:) < (1.0_wp - eps8))
       hlp2(:,:) = (hlp1(:,:) - 1.0_wp) * LOG(1.0_wp - hlp1(:,:)) / hlp1(:,:)
    ELSEWHERE
       hlp2(:,:) = 0.0_wp
    ENDWHERE

    !> 2.1 N15 fluxes involving NO3
    !!
    plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_local_sl(:,:) * & 
      &                                 calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), &
      &                                                    eta_plant_uptake_no3 * hlp2(:,:))
    mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_local_sl(:,:) * & 
      &                                 calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), &
      &                                                    eta_plant_uptake_no3 * hlp2(:,:))
    microbial_uptake_no3_n15_sl(:,:)  = microbial_uptake_no3_local_sl(:,:) * & 
      &                                 calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), &
      &                                                    eta_mic_uptake_no3 * hlp2(:,:)) 
    denitrification_noy_n15(:,:)      = denitrification_noy_local(:,:) * & 
      &                                 calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), &
      &                                                    eta_denitrification * hlp2(:,:))
    denitrification_n2o_n15(:,:)      = denitrification_n2o_local(:,:) * & 
      &                                 calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), &
      &                                                    eta_denitrification * hlp2(:,:))
    denitrification_n2_n15(:,:)       = denitrification_n2_local(:,:) * & 
      &                                 calc_fractionation(no3_solute(:,:), no3_n15_solute(:,:), &
      &                                                    eta_denitrification * hlp2(:,:)) 

    !> 2.2 contrain 15NO3 fluxes to available 15NO3
    !!
    hlp1(:,:) = plant_uptake_no3_n15_sl(:,:) + mycorrhiza_uptake_no3_n15_sl(:,:) + microbial_uptake_no3_n15_sl(:,:) + & 
                denitrification_noy_n15(:,:) + denitrification_n2o_n15(:,:) + denitrification_n2_n15(:,:)
    hlp2(:,:) = no3_n15_solute(:,:) + transport_no3_n15_solute(:,:)
    WHERE(hlp2(:,:) < hlp1(:,:))
       WHERE(hlp1(:,:) > eps8)
          plant_uptake_no3_n15_sl(:,:)      = hlp2(:,:) * plant_uptake_no3_n15_sl(:,:) / hlp1(:,:) 
          mycorrhiza_uptake_no3_n15_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_no3_n15_sl(:,:) / hlp1(:,:) 
          microbial_uptake_no3_n15_sl(:,:)  = hlp2(:,:) * microbial_uptake_no3_n15_sl(:,:) / hlp1(:,:) 
          denitrification_noy_n15(:,:)      = hlp2(:,:) * denitrification_noy_n15(:,:) / hlp1(:,:)
          denitrification_n2o_n15(:,:)      = hlp2(:,:) * denitrification_n2o_n15(:,:) / hlp1(:,:)
          denitrification_n2_n15(:,:)       = hlp2(:,:) * denitrification_n2_n15(:,:) / hlp1(:,:)
       ELSEWHERE
          plant_uptake_no3_n15_sl(:,:)      = 0.0_wp 
          mycorrhiza_uptake_no3_n15_sl(:,:) = 0.0_wp 
          microbial_uptake_no3_n15_sl(:,:)  = 0.0_wp
          denitrification_noy_n15(:,:)      = 0.0_wp
          denitrification_n2o_n15(:,:)      = 0.0_wp
          denitrification_n2_n15(:,:)       = 0.0_wp
       ENDWHERE
    ENDWHERE

    !> 3.0 15N signal of asymbiotic fixation
    !! 
    asymb_n_fixation_n15(:,:) = asymb_n_fixation_local(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(-eta_nfixation))

    !> 4.0 unit conversion: [mol/m3/timestep] -> [micro-mol/m3/s]
    !!  
    !! 'INTENT(out)' variables unit conversion when this routines is called from JSM
    !!
    IF(sb_model_scheme == 'jsm') THEN
      plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_n15_sl(:,:)      / dtime * 1.e6_wp
      plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_n15_sl(:,:)      / dtime * 1.e6_wp
      mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_n15_sl(:,:) / dtime * 1.e6_wp
      mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_n15_sl(:,:) / dtime * 1.e6_wp
      microbial_uptake_nh4_n15_sl(:,:)  = microbial_uptake_nh4_n15_sl(:,:)  / dtime * 1.e6_wp
      microbial_uptake_no3_n15_sl(:,:)  = microbial_uptake_no3_n15_sl(:,:)  / dtime * 1.e6_wp
      asymb_n_fixation_n15(:,:)         = asymb_n_fixation_n15(:,:)         / dtime * 1.e6_wp
      nitrification_no3_n15(:,:)        = nitrification_no3_n15(:,:)        / dtime * 1.e6_wp
      nitrification_noy_n15(:,:)        = nitrification_noy_n15(:,:)        / dtime * 1.e6_wp
      nitrification_n2o_n15(:,:)        = nitrification_n2o_n15(:,:)        / dtime * 1.e6_wp
      denitrification_noy_n15(:,:)      = denitrification_noy_n15(:,:)      / dtime * 1.e6_wp
      denitrification_n2o_n15(:,:)      = denitrification_n2o_n15(:,:)      / dtime * 1.e6_wp
      denitrification_n2_n15(:,:)       = denitrification_n2_n15(:,:)       / dtime * 1.e6_wp
    ENDIF

  END SUBROUTINE calc_sb_inorganic_n15_fluxes 

  ! ======================================================================================================= !
  !>Calculates the formation of SOM resulting from the decomposition: dPrec/dt = CUE_Pdon * dPdon/dt
  !>
  !>  gross mineralisation follows the stoichiometric difference of the donor and recipient pool, taking account for CUE
  !>  for fast -> slow -> fast only causes gross mineralisation because of the C:N:P and CUE
  !>  for litter -> fast gross mineralisation is typcially less than microbial uptake from 
  !>                mineral pools (or organic matter for P)
  !>
  SUBROUTINE calc_SOM_formation( &
    & nc, &
    & nsoil_sb, &
    & apparent_fast_som_cn, &
    & apparent_fast_som_np, &
    & microbial_uptake_nh4_sl, &
    & microbial_uptake_nh4_n15_sl, &
    & microbial_uptake_no3_sl, &
    & microbial_uptake_no3_n15_sl, &
    & asymb_n_fixation, &
    & asymb_n_fixation_n15, &
    & microbial_uptake_po4_sl, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & net_mineralisation_nh4, &
    & net_mineralisation_no3, &
    & net_mineralisation_po4, &
    & net_mineralisation_nh4_n15, &
    & net_mineralisation_no3_n15, &
    & biochem_mineralisation_po4, &
    & sb_loss_mt, &
    & sb_formation_mt )
 
    USE mo_isotope_util,        ONLY: calc_fractionation
    USE mo_sb_constants,        ONLY: frac_litter2fast_som, frac_fast2slow_som, frac_slow2fast_som, &
      &                               k_slow_som_cn, k_slow_som_np, microbial_nue, microbial_pue, & 
      &                               eta_ammonification, eta_nitrification
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                              !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                        !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: apparent_fast_som_cn, &         !< apparent carbon-nitrogen ratio of microbial uptake [mol/mol] 
                                                        apparent_fast_som_np, &         !< apparent carbon-nitrogen ratio of microbial uptake [mol/mol]
                                                        microbial_uptake_nh4_sl, &      !< microbial NH4 uptake [mol / m3 / timestep]
                                                        microbial_uptake_nh4_n15_sl, &  !< microbial 15NH4 uptake [mol / m3 / timestep]
                                                        microbial_uptake_no3_sl, &      !< microbial NO3 uptake [mol / m3 / timestep]
                                                        microbial_uptake_no3_n15_sl, &  !< microbial 15NO3 uptake [mol / m3 / timestep]
                                                        asymb_n_fixation, &             !< asymbiotic N fixation [mol / m3 / timestep]
                                                        asymb_n_fixation_n15, &         !< asymbiotic 15N fixation [mol / m3 / timestep]
                                                        microbial_uptake_po4_sl         !< microbial PO4 uptake [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: het_respiration, &              !< heterotrophic respiration [mol C / m3 / timestep]
                                                        het_respiration_c13, &          !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                        het_respiration_c14, &          !< heterotrophic respiration [micro-mol 14C / m3 / timestep]
                                                        net_mineralisation_nh4, &       !< net mineralisation of NH4 [mol / m3 / timestep]
                                                        net_mineralisation_no3, &       !< net mineralisation of NO3 [mol / m3 / timestep]
                                                        net_mineralisation_po4, &       !< net mineralisation of PO4 [mol / m3 / timestep]
                                                        net_mineralisation_nh4_n15, &   !< net mineralisation of 15NH4 [mol / m3 / timestep]
                                                        net_mineralisation_no3_n15, &   !< net mineralisation of 15NO3 [mol / m3 / timestep]
                                                        biochem_mineralisation_po4      !< biomineralisation rate [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_loss_mt(:,:,:,:)             !< bgcm sb_loss flux: turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)        !< bgcm sb_formation flux: formation of SOM and litter pools  [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp2
    CHARACTER(len=*), PARAMETER                      :: routine = TRIM(modname)//':calc_SOM_formation'
    ! ----------------------------------------------------------------------------------------------------- !
    
    !>1.0 formation of fast SOM from litter 
    !>
    ! calculate 15N mineralised during litter decay, remainder enters fast pool
    hlp1(:,:) = (1.0_wp - microbial_nue) * sb_loss_mt(ix_soluable_litter, ixN, :, :) &
      &         * calc_fractionation(sb_loss_mt(ix_soluable_litter, ixN, :, :), &
      &                              sb_loss_mt(ix_soluable_litter, ixN15, :, :), &
      &                              eta_ammonification)
    hlp2(:,:) = (1.0_wp - microbial_nue) * sb_loss_mt(ix_polymeric_litter, ixN, :, :) &
      &         * calc_fractionation(sb_loss_mt(ix_polymeric_litter, ixN, :, :), &
      &                              sb_loss_mt(ix_polymeric_litter, ixN15, :, :), &
      &                              eta_ammonification)
    ! all elements
    sb_formation_mt(ix_microbial, ixC, :, :)   = sb_formation_mt(ix_microbial, ixC, :, :) + frac_litter2fast_som &
      &                                          * (sb_loss_mt(ix_soluable_litter, ixC, :, :) &
      &                                          + sb_loss_mt(ix_polymeric_litter, ixC, :, :)) 
    sb_formation_mt(ix_microbial, ixN, :, :)   = sb_formation_mt(ix_microbial, ixN, :, :) + microbial_nue &
      &                                          * (sb_loss_mt(ix_soluable_litter, ixN, :, :) &
      &                                          + sb_loss_mt(ix_polymeric_litter, ixN, :, :)) &
      &                                          + microbial_uptake_nh4_sl(:,:) + microbial_uptake_no3_sl(:,:) &
      &                                          + asymb_n_fixation(:,:)
    sb_formation_mt(ix_microbial, ixP, :, :)   = sb_formation_mt(ix_microbial, ixP, :, :) + microbial_pue &
      &                                          * (sb_loss_mt(ix_soluable_litter, ixP, :, :) &
      &                                          + sb_loss_mt(ix_polymeric_litter, ixP, :, :)) &
      &                                          + microbial_uptake_po4_sl(:,:)
    sb_formation_mt(ix_microbial, ixC13, :, :) = sb_formation_mt(ix_microbial, ixC13, :, :) + frac_litter2fast_som &
      &                                          * (sb_loss_mt(ix_soluable_litter, ixC13, :, :) &
      &                                          + sb_loss_mt(ix_polymeric_litter, ixC13, :, :))
    sb_formation_mt(ix_microbial, ixC14, :, :) = sb_formation_mt(ix_microbial, ixC14, :, :) + frac_litter2fast_som &
      &                                          * (sb_loss_mt(ix_soluable_litter, ixC14, :, :) &
      &                                          + sb_loss_mt(ix_polymeric_litter, ixC14, :, :))
    sb_formation_mt(ix_microbial, ixN15, :, :) = sb_formation_mt(ix_microbial, ixN15, :, :) & 
      &                                          + sb_loss_mt(ix_soluable_litter, ixN15, :, :)  - hlp1(:,:) &
      &                                          + sb_loss_mt(ix_polymeric_litter, ixN15, :, :) - hlp2(:,:) &
      &                                          + microbial_uptake_nh4_n15_sl(:,:) + microbial_uptake_no3_n15_sl(:,:) &
      &                                          + asymb_n_fixation_n15(:,:)

    het_respiration(:,:)            = het_respiration(:,:) + (1._wp - frac_litter2fast_som) * & 
      &                                                      (sb_loss_mt(ix_soluable_litter, ixC, :, :) + &  
      &                                                       sb_loss_mt(ix_polymeric_litter, ixC, :, :))
    het_respiration_c13(:,:)        = het_respiration_c13(:,:) + (1._wp - frac_litter2fast_som) * & 
      &                                                          (sb_loss_mt(ix_soluable_litter, ixC13, :, :) + &  
      &                                                           sb_loss_mt(ix_polymeric_litter, ixC13, :, :)) 
    het_respiration_c14(:,:)        = het_respiration_c14(:,:) + (1._wp - frac_litter2fast_som) * & 
      &                                                          (sb_loss_mt(ix_soluable_litter, ixC14, :, :) + &  
      &                                                           sb_loss_mt(ix_polymeric_litter, ixC14, :, :))
    net_mineralisation_nh4(:,:)     = net_mineralisation_nh4(:,:) - microbial_uptake_nh4_sl(:,:) + &  
      &                                                             (1.0_wp - microbial_nue) * &
      &                                                             (sb_loss_mt(ix_soluable_litter, ixN, :, :) + &  
      &                                                              sb_loss_mt(ix_polymeric_litter, ixN, :, :))
    net_mineralisation_no3(:,:)     = net_mineralisation_no3(:,:) - microbial_uptake_no3_sl(:,:) 
    net_mineralisation_po4(:,:)     = net_mineralisation_po4(:,:) - microbial_uptake_po4_sl(:,:) + & 
      &                                                             (1.0_wp - microbial_pue) * &
      &                                                             (sb_loss_mt(ix_soluable_litter, ixP, :, :) + &  
      &                                                              sb_loss_mt(ix_polymeric_litter, ixP, :, :))
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) - microbial_uptake_nh4_n15_sl(:,:) + hlp1(:,:) + hlp2(:,:)
    net_mineralisation_no3_n15(:,:) = net_mineralisation_no3_n15(:,:) - microbial_uptake_no3_n15_sl(:,:) 

    !>2.0 formation of slow SOM (= residue pool) from fast SOM decay
    !>
    sb_formation_mt(ix_residue, ixC, :, :)   = sb_formation_mt(ix_residue, ixC, :, :) &
      &                                        + frac_fast2slow_som &
      &                                        * sb_loss_mt(ix_microbial, ixC, :, :)
    sb_formation_mt(ix_residue, ixN, :, :)   = sb_formation_mt(ix_residue, ixN, :, :) &
      &                                        + frac_fast2slow_som / k_slow_som_cn &
      &                                        * sb_loss_mt(ix_microbial, ixC, :, :)
    sb_formation_mt(ix_residue, ixP, :, :)   = sb_formation_mt(ix_residue, ixP, :, :) &
      &                                        + frac_fast2slow_som / k_slow_som_cn / k_slow_som_np &
      &                                        * sb_loss_mt(ix_microbial, ixC, :, :)
    sb_formation_mt(ix_residue, ixC13, :, :) = sb_formation_mt(ix_residue, ixC13, :, :) &
      &                                        + frac_fast2slow_som &
      &                                        * sb_loss_mt(ix_microbial, ixC13, :, :)
    sb_formation_mt(ix_residue, ixC14, :, :) = sb_formation_mt(ix_residue, ixC14, :, :) &
      &                                        + frac_fast2slow_som &
      &                                        * sb_loss_mt(ix_microbial, ixC14, :, :)
    ! store formation of slow 15N for later use 
    hlp1(:,:) = frac_fast2slow_som / k_slow_som_cn * sb_loss_mt(ix_microbial, ixC, :, :) &
      &         * calc_fractionation(sb_loss_mt(ix_microbial, ixN, :, :), &
      &                              sb_loss_mt(ix_microbial, ixN15, :, :), &
      &                              0.0_wp)
    sb_formation_mt(ix_residue, ixN15, :, :) = sb_formation_mt(ix_residue, ixN15, :, :) + hlp1(:,:) 

    het_respiration(:,:)        = het_respiration(:,:)     + (1._wp - frac_fast2slow_som) * sb_loss_mt(ix_microbial, ixC, :, :)
    het_respiration_c13(:,:)    = het_respiration_c13(:,:) + (1._wp - frac_fast2slow_som) * sb_loss_mt(ix_microbial, ixC13, :, :)
    het_respiration_c14(:,:)    = het_respiration_c14(:,:) + (1._wp - frac_fast2slow_som) * sb_loss_mt(ix_microbial, ixC14, :, :)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb_loss_mt(ix_microbial, ixN, :, :) &
      &                           - frac_fast2slow_som / k_slow_som_cn * sb_loss_mt(ix_microbial, ixC, :, :)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb_loss_mt(ix_microbial, ixP, :, :) &
      &                           - frac_fast2slow_som / k_slow_som_cn / k_slow_som_np &
      &                           * sb_loss_mt(ix_microbial, ixC, :, :)
    ! store net mineralisation of 15N for later use
    hlp2(:,:) = (sb_loss_mt(ix_microbial, ixN, :, :) - frac_fast2slow_som / k_slow_som_cn &
      &         * sb_loss_mt(ix_microbial, ixC, :, :)) &
      &         * calc_fractionation(sb_loss_mt(ix_microbial, ixN, :, :), &
      &                              sb_loss_mt(ix_microbial, ixN15, :, :), &
      &                              eta_ammonification)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + hlp2(:,:)
    net_mineralisation_no3_n15(:,:) = net_mineralisation_no3_n15(:,:)

    ! adjust fast 15N loss rate to account for discrimination during ammonification
    sb_loss_mt(ix_microbial, ixN15, :, :) = sb_loss_mt(ix_microbial, ixN15, :, :) &
      &                                     - (sb_loss_mt(ix_microbial, ixN15, :, :) - hlp1(:,:) - hlp2(:,:))

    !>3.0 formation of fast SOM (= microbial pool) from slow SOM decay 
    !>
    sb_formation_mt(ix_microbial, ixC, :, :)   = sb_formation_mt(ix_microbial, ixC, :, :) &
      &                                          + frac_slow2fast_som &
      &                                          * sb_loss_mt(ix_residue, ixC, :, :)
    sb_formation_mt(ix_microbial, ixN, :, :)   = sb_formation_mt(ix_microbial, ixN, :, :) &
      &                                          + frac_slow2fast_som / apparent_fast_som_cn(:,:) &
      &                                          * sb_loss_mt(ix_residue, ixC, :, :)
    sb_formation_mt(ix_microbial, ixP, :, :)   = sb_formation_mt(ix_microbial, ixP, :, :) & 
      &                                          + frac_slow2fast_som / apparent_fast_som_cn(:,:) / apparent_fast_som_np(:,:) &
      &                                          * sb_loss_mt(ix_residue, ixC, :, :)
    sb_formation_mt(ix_microbial, ixC13, :, :) = sb_formation_mt(ix_microbial, ixC13, :, :) &
      &                                          + frac_slow2fast_som &
      &                                          * sb_loss_mt(ix_residue, ixC13, :, :)
    sb_formation_mt(ix_microbial, ixC14, :, :) = sb_formation_mt(ix_microbial, ixC14, :, :) &
      &                                          + frac_slow2fast_som &
      &                                          * sb_loss_mt(ix_residue, ixC14, :, :)
    ! store formation of fast 15N for later use 
    hlp1(:,:) = frac_slow2fast_som / apparent_fast_som_cn(:,:) * sb_loss_mt(ix_residue, ixC, :, :) &
      &         * calc_fractionation(sb_loss_mt(ix_residue, ixN, :, :), &
      &                              sb_loss_mt(ix_residue, ixN15, :, :), &
      &                              0.0_wp)
    sb_formation_mt(ix_microbial, ixN15, :, :) = sb_formation_mt(ix_microbial, ixN15, :, :) + hlp1(:,:) 

    het_respiration(:,:)        = het_respiration(:,:)     + (1._wp - frac_slow2fast_som) * sb_loss_mt(ix_residue, ixC, :, :)
    het_respiration_c13(:,:)    = het_respiration_c13(:,:) + (1._wp - frac_slow2fast_som) * sb_loss_mt(ix_residue, ixC13, :, :)
    het_respiration_c14(:,:)    = het_respiration_c14(:,:) + (1._wp - frac_slow2fast_som) * sb_loss_mt(ix_residue, ixC14, :, :)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb_loss_mt(ix_residue, ixN, :, :) &
      &                           - frac_slow2fast_som / apparent_fast_som_cn(:,:) * sb_loss_mt(ix_residue, ixC, :, :)
    ! update the loss of residue phosphorus due to biochem_mineralisation_po4, restricted by the CNP ratio of slow SOM
    biochem_mineralisation_po4(:,:) = MIN(biochem_mineralisation_po4(:,:), sb_loss_mt(ix_residue, ixP, :, :) &
      &                               - frac_slow2fast_som / apparent_fast_som_cn(:,:) / apparent_fast_som_np(:,:) &
      &                               * sb_loss_mt(ix_residue, ixC, :, :))
    net_mineralisation_po4(:,:)     = net_mineralisation_po4(:,:) + sb_loss_mt(ix_residue, ixP, :, :) &
      &                               - biochem_mineralisation_po4(:,:) &
      &                               - frac_slow2fast_som / apparent_fast_som_cn(:,:) / apparent_fast_som_np(:,:) &
      &                               * sb_loss_mt(ix_residue, ixC, :, :)
    ! store net mineralisation of 15N for later use
    hlp2(:,:) = (sb_loss_mt(ix_residue, ixN, :, :) - frac_slow2fast_som / apparent_fast_som_cn(:,:) &
      &         * sb_loss_mt(ix_residue, ixC, :, :)) &
      &         * calc_fractionation(sb_loss_mt(ix_residue, ixN, :, :), &
      &                              sb_loss_mt(ix_residue, ixN15, :, :), &
      &                              eta_ammonification)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + hlp2(:,:)

    ! adjust slow 15N loss rate to account for discrimination during ammonification
    sb_loss_mt(ix_residue, ixN15, :, :) = sb_loss_mt(ix_residue, ixN15, :, :) &
      &                                   - (sb_loss_mt(ix_residue, ixN15, :, :) - hlp1(:,:) - hlp2(:,:))
  END SUBROUTINE calc_SOM_formation  
  
  ! ======================================================================================================= !
  !>Calculates the formation of metabolic and structural litter from woody litter decay. 
  !>
  !>  Some C,N,P becomes lost during processing prescribed by the carbon and nutrient use
  !>  efficiencies. The remainder becomes incorporated into metabolic and structural litter
  !>  given the remaining litter C:N and the prescribed lignin content.
  !>
  SUBROUTINE calc_litter_formation_from_woody_decay( &
    & nc, &
    & nsoil_sb, &
    & sb_loss_mt, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & net_mineralisation_nh4, &
    & net_mineralisation_po4, &
    & net_mineralisation_nh4_n15, &
    & sb_formation_mt )
  
    USE mo_sb_constants,        ONLY: frac_woody2dec_litter, microbial_nue, microbial_pue, & 
      &                               fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, lc_woody_litter
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc, &                         !< dimensions
                                                        nsoil_sb
    REAL(wp),                          INTENT(in)    :: sb_loss_mt(:,:,:,:)           !< bgcm sb_loss flux: turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: het_respiration, &            !< heterotrophic respiration [mol C / m3 / timestep]
                                                        het_respiration_c13, &        !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                        het_respiration_c14, &        !< heterotrophic respiration [micro-mol 14C / m3 / timestep]
                                                        net_mineralisation_nh4, &     !< net mineralisation of NH4 [mol / m3 / timestep]
                                                        net_mineralisation_po4, &     !< net mineralisation of PO4 [mol / m3 / timestep]
                                                        net_mineralisation_nh4_n15    !< net mineralisation of 15NH4 [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)      !< bgcm sb_formation flux: formation of SOM and litter pools  [mol / m3 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_sb) :: fc_soluable_woody_litter, & !< fraction of non-respired woody litter becoming metabolic litter 
                                         fn_soluable_woody_litter, & 
                                         fp_soluable_woody_litter
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_formation_from_woody_decay'
    ! ----------------------------------------------------------------------------------------------------- !
    
    !>1.0 Calculate partitioning of decaying woody litter into metabolic and structural litter
    !>
    fc_soluable_woody_litter(:,:) = 0.0_wp
    fn_soluable_woody_litter(:,:) = 0.0_wp
    fp_soluable_woody_litter(:,:) = 0.0_wp
    WHERE(sb_loss_mt(ix_woody_litter, ixN, :, :) > eps8) 
      fc_soluable_woody_litter(:,:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_woody_litter &
        &                             * sb_loss_mt(ix_woody_litter, ixC, :, :) / sb_loss_mt(ix_woody_litter, ixN, :, :))
      WHERE(fc_soluable_woody_litter(:,:) > eps8)
        fn_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:)) &
          &                             / (k_fn_soluable * fc_soluable_woody_litter(:,:)))
        fp_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:)) &
          &                             / (k_fp_soluable * fc_soluable_woody_litter(:,:)))
      ENDWHERE
    ENDWHERE

    !>2.0 formation of metabolic and structural litter
    !>     taking account of the processing losses of C, N and P
    !>
    sb_formation_mt(ix_soluable_litter, ixC, :, :)   = sb_formation_mt(ix_soluable_litter, ixC, :, :) &
      &                                                + fc_soluable_woody_litter(:,:) * frac_woody2dec_litter &
      &                                                * sb_loss_mt(ix_woody_litter, ixC, :, :) 
    sb_formation_mt(ix_soluable_litter, ixN, :, :)   = sb_formation_mt(ix_soluable_litter, ixN, :, :) &
      &                                                + fn_soluable_woody_litter(:,:) * microbial_nue &
      &                                                * sb_loss_mt(ix_woody_litter, ixN, :, :) 
    sb_formation_mt(ix_soluable_litter, ixP, :, :)   = sb_formation_mt(ix_soluable_litter, ixP, :, :) &
      &                                                + fp_soluable_woody_litter(:,:) * microbial_pue &
      &                                                * sb_loss_mt(ix_woody_litter, ixP, :, :) 
    sb_formation_mt(ix_soluable_litter, ixC13, :, :) = sb_formation_mt(ix_soluable_litter, ixC13, :, :) &
      &                                                + fc_soluable_woody_litter(:,:) * frac_woody2dec_litter &
      &                                                * sb_loss_mt(ix_woody_litter, ixC13, :, :) 
    sb_formation_mt(ix_soluable_litter, ixC14, :, :) = sb_formation_mt(ix_soluable_litter, ixC14, :, :) &
      &                                                + fc_soluable_woody_litter(:,:) * frac_woody2dec_litter &
      &                                                * sb_loss_mt(ix_woody_litter, ixC14, :, :) 
    sb_formation_mt(ix_soluable_litter, ixN15, :, :) = sb_formation_mt(ix_soluable_litter, ixN15, :, :) &
      &                                                + fn_soluable_woody_litter(:,:) * microbial_nue &
      &                                                * sb_loss_mt(ix_woody_litter, ixN15, :, :) 

    sb_formation_mt(ix_polymeric_litter, ixC, :, :)   = sb_formation_mt(ix_polymeric_litter, ixC, :, :) &
      &                                                 + (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter &
      &                                                 * sb_loss_mt(ix_woody_litter, ixC, :, :) 
    sb_formation_mt(ix_polymeric_litter, ixN, :, :)   = sb_formation_mt(ix_polymeric_litter, ixN, :, :) &
      &                                                 + (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue &
      &                                                 * sb_loss_mt(ix_woody_litter, ixN, :, :) 
    sb_formation_mt(ix_polymeric_litter, ixP, :, :)   = sb_formation_mt(ix_polymeric_litter, ixP, :, :) &
      &                                                 + (1.0_wp - fp_soluable_woody_litter(:,:)) * microbial_pue &
      &                                                 * sb_loss_mt(ix_woody_litter, ixP, :, :) 
    sb_formation_mt(ix_polymeric_litter, ixC13, :, :) = sb_formation_mt(ix_polymeric_litter, ixC13, :, :) &
      &                                                 + (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter &
      &                                                 * sb_loss_mt(ix_woody_litter, ixC13, :, :) 
    sb_formation_mt(ix_polymeric_litter, ixC14, :, :) = sb_formation_mt(ix_polymeric_litter, ixC14, :, :) &
      &                                                 + (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter &
      &                                                 * sb_loss_mt(ix_woody_litter, ixC14, :, :) 
    sb_formation_mt(ix_polymeric_litter, ixN15, :, :) = sb_formation_mt(ix_polymeric_litter, ixN15, :, :) &
      &                                                 + (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue &
      &                                                 * sb_loss_mt(ix_woody_litter, ixN15, :, :) 

    !>3.0 record losses to heterotrophic respiration and mineralisation
    !>
    het_respiration(:,:)            = het_respiration(:,:) & 
      &                               + (1.0_wp - frac_woody2dec_litter) &
      &                               * sb_loss_mt(ix_woody_litter, ixC, :, :)
    net_mineralisation_nh4(:,:)     = net_mineralisation_nh4(:,:) &
      &                               + (1.0_wp - microbial_nue) &
      &                               * sb_loss_mt(ix_woody_litter, ixN, :, :)
    net_mineralisation_po4(:,:)     = net_mineralisation_po4(:,:) &
      &                               + (1.0_wp - microbial_pue) &
      &                               * sb_loss_mt(ix_woody_litter, ixP, :, :)
    het_respiration_c13(:,:)        = het_respiration_c13(:,:) & 
      &                               + (1.0_wp - frac_woody2dec_litter) &
      &                               * sb_loss_mt(ix_woody_litter, ixC13, :, :)
    het_respiration_c14(:,:)        = het_respiration_c14(:,:) & 
      &                               + (1.0_wp - frac_woody2dec_litter) &
      &                               * sb_loss_mt(ix_woody_litter, ixC14, :, :)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) &
      &                               + (1.0_wp - microbial_nue) &
      &                               * sb_loss_mt(ix_woody_litter, ixN15, :, :)
  END SUBROUTINE calc_litter_formation_from_woody_decay

  ! ======================================================================================================= !
  !>calculates matter transport by bioturbation, by calling calc_bioturbation_transport
  !>
  !> Input:  pools subject to bioturbation, percolation rate, and soil depth
  !>         This explicitly excludes mycorrhiza, as these die when eaten
  !>
  !> Output: vertical transport rate (mol/m2/timestep)
  !>
  SUBROUTINE calc_bioturbation_transport_wrapper_s1d( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & k_bioturb, &
    & nh4_assoc, &
    & nh4_n15_assoc, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & sb_pool_mt, &
    & transport_nh4_assoc, &
    & transport_nh4_n15_assoc, &
    & transport_po4_assoc_fast, &
    & transport_po4_assoc_slow, & 
    & transport_po4_occluded, &
    & sb_transport_mt)

    USE mo_sb_jsm_transport,        ONLY: calc_bioturbation_transport
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                          !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                    !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl               !< depth of each soil layer [m] 
    INTEGER,                           INTENT(in)    :: elements_index_map(:)       !< map bgcm element ID -> IDX
    LOGICAL,                           INTENT(in)    :: is_element_used(:)          !< is element in 'elements_index_map' used 
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: k_bioturb                   !< diffusion factor for bioturbation 
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_assoc                   !< adsorbed NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: nh4_n15_assoc               !< adsorbed 15NH4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: po4_assoc_fast              !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: po4_assoc_slow              !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: po4_occluded                !< occluded PO4 pool [mol/m3]
    REAL(wp),                          INTENT(in)    :: sb_pool_mt(:, :, :, :)      !< sb_pool
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_assoc         !< rate of associated NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_nh4_n15_assoc     !< rate of associated 15NH4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_po4_assoc_fast    !< rate of fast associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_po4_assoc_slow    !< rate of slow associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: transport_po4_occluded      !< rate of occluded PO4 transport [mol/m3/timestep]
    REAL(wp),                          INTENT(inout) :: sb_transport_mt(:, :, :, :) !< sb_transport flux
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ielem             !< loop over bgcm elements
    INTEGER                     :: ix_elem           !< index of element in bgcm, used for looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport_wrapper_s1d'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 organic pools
    !>
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        !>  1.1 fast SOM pool
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), & 
          &                              sb_pool_mt(ix_soluable_litter, ix_elem, :, :), &
          &                              sb_transport_mt(ix_soluable_litter, ix_elem, :, :))
        CALL calc_bioturbation_transport(nc, nsoil_sb, &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), & 
          &                              sb_pool_mt(ix_polymeric_litter, ix_elem, :, :), &
          &                              sb_transport_mt(ix_polymeric_litter, ix_elem, :, :))
        CALL calc_bioturbation_transport(nc, nsoil_sb, &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), & 
          &                              sb_pool_mt(ix_microbial, ix_elem, :, :), &
          &                              sb_transport_mt(ix_microbial, ix_elem, :, :))
        !>  1.1 slow SOM pool
        !>
        CALL calc_bioturbation_transport(nc, nsoil_sb, &
          &                              k_bioturb(:,:), &
          &                              soil_depth_sl(:,:), & 
          &                              sb_pool_mt(ix_residue, ix_elem, :, :), &
          &                              sb_transport_mt(ix_residue, ix_elem, :, :))
      END IF
    END DO

    !>2.0 inorganic pools
    !>
    CALL calc_bioturbation_transport(nc, nsoil_sb, &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), & 
      &                              nh4_assoc(:,:), transport_nh4_assoc(:,:)) 
    CALL calc_bioturbation_transport(nc, nsoil_sb, &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), & 
      &                              nh4_n15_assoc(:,:), transport_nh4_n15_assoc(:,:)) 
    CALL calc_bioturbation_transport(nc, nsoil_sb, &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), & 
      &                              po4_assoc_fast(:,:), transport_po4_assoc_fast(:,:)) 
    CALL calc_bioturbation_transport(nc, nsoil_sb, &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), & 
      &                              po4_assoc_slow(:,:), transport_po4_assoc_slow(:,:)) 
    CALL calc_bioturbation_transport(nc, nsoil_sb, &
      &                              k_bioturb(:,:),&
      &                              soil_depth_sl(:,:), & 
      &                              po4_occluded(:,:), transport_po4_occluded(:,:)) 
  END SUBROUTINE calc_bioturbation_transport_wrapper_s1d

  ! ======================================================================================================= !
  !>For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to net mineralisation
  !>
  SUBROUTINE calc_fixed_nloss_rate( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & net_mineralisation_nh4, & 
    & net_mineralisation_nh4_n15, & 
    & nitrification_no3, &
    & nitrification_no3_n15, &
    & volatilisation_nh4, &
    & volatilisation_nh4_n15, &
    & emission_n2, &
    & emission_n2_n15)

    USE mo_sb_constants,           ONLY: floss_nmin, fnit_nmin
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                                          INTENT(in)  :: nc                            !< dimensions
    INTEGER,                                          INTENT(in)  :: nsoil_sb                      !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(in)  :: soil_depth_sl, &              !< soil depth per layer [m]
                                                                     net_mineralisation_nh4, &     !< net mineralisation rate of NH4 [mol / m3 / timestep]
                                                                     net_mineralisation_nh4_n15    !< net mineralisation rate of 15NH4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc, nsoil_sb),                INTENT(out) :: nitrification_no3, &          !< nitrification to NO3 [mol / m3 / timestep]
                                                                     nitrification_no3_n15, &      !< nitrification to 15NO3 [mol / m3 / timestep]
                                                                     volatilisation_nh4, &         !< volatilisation of NH4 [mol / m3 / timestep]
                                                                     volatilisation_nh4_n15        !< volatilisation of 15NH4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(nc),                          INTENT(out) :: emission_n2, &                !< soil efflux of N2 [mol / m2 / timestep]
                                                                     emission_n2_n15               !< soil efflux of 15N2 [mol / m2 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fixed_nloss_rate'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 case of fixed N losses: N loss proportional to net mineralisation
    !>
    WHERE(net_mineralisation_nh4(:,:) > 0.0_wp)
      ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
      volatilisation_nh4(:,:)      = net_mineralisation_nh4(:,:)     * floss_nmin 
      volatilisation_nh4_n15(:,:)  = net_mineralisation_nh4_n15(:,:) * floss_nmin 
      ! NH4 transformed to NO3 [mol m-3 timestep-1]
      nitrification_no3(:,:)       = net_mineralisation_nh4(:,:)     * (1._wp - floss_nmin) * fnit_nmin 
      nitrification_no3_n15(:,:)   = net_mineralisation_nh4_n15(:,:) * (1._wp - floss_nmin) * fnit_nmin 
    ENDWHERE
    ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
    emission_n2(:)     = SUM(volatilisation_nh4(:,:)     * soil_depth_sl(:,:), DIM=2)
    emission_n2_n15(:) = SUM(volatilisation_nh4_n15(:,:) * soil_depth_sl(:,:), DIM=2)
  END SUBROUTINE calc_fixed_nloss_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to deal with mycorrhiza 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_mycorrhiza_dynamics( &
    & nc, &
    & nsoil_sb, &
    & soil_depth_sl, &
    & elements_index_map, &
    & is_element_used, &
    & lctlib_tau_mycorrhiza, &
    & flag_mycorrhiza_org, &
    & veg_pool_fine_root_carbon , &
    & sb_pool_mt, &
    & mycorrhiza_uptake_nh4_sl, &
    & mycorrhiza_uptake_nh4_n15_sl, &
    & mycorrhiza_uptake_no3_sl, &
    & mycorrhiza_uptake_no3_n15_sl, &
    & mycorrhiza_uptake_norg_sl, &
    & mycorrhiza_uptake_norg_n15_sl, &
    & mycorrhiza_uptake_po4_sl, &
    & nh4_solute, &
    & no3_solute, &
    & po4_solute, &
    & km2_uptake_nh4_act, &
    & km2_uptake_no3_act, &
    & km2_uptake_po4_act, &
    & km_uptake_org_act, &
    & het_respiration, &
    & het_respiration_c13, &
    & het_respiration_c14, &
    & myc_respiration, &
    & myc_respiration_c13, &
    & myc_respiration_c14, &
    & sb_formation_mt, &
    & sb_loss_mt, &
    & sb_mycorrhiza_export_mt) 

    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: mycorrhiza_cn, mycorrhiza_np, mycorrhiza_cue, &
                                          tau_slow, vmax_uptake_norg, f_nin2norg, f_norg_limit, organic_cn
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                           INTENT(in)    :: nc                                !< dimensions
    INTEGER,                           INTENT(in)    :: nsoil_sb                          !< number of soil layers
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(in)    :: soil_depth_sl                     !< depth per soil layer
    INTEGER,                           INTENT(in)    :: elements_index_map(:)             !< map bgcm element ID -> IDX
    LOGICAL,                           INTENT(in)    :: is_element_used(:)                !< is element in 'elements_index_map' used 
    REAL(wp),                          INTENT(in)    :: lctlib_tau_mycorrhiza             !< lctlib parameter
    LOGICAL,                           INTENT(in)    :: flag_mycorrhiza_org               !< from sb config
    REAL(wp), DIMENSION(nc),           INTENT(in)    :: veg_pool_fine_root_carbon         !< fine root carbon [mol / m2]
    REAL(wp),                          INTENT(in)    :: sb_pool_mt(:,:,:,:)               !< bgcm sb_pool: soil biogeochemistry pools [mol / m3]
    REAL(wp), DIMENSION(nc, nsoil_sb), INTENT(inout) :: mycorrhiza_uptake_nh4_sl, &       !< NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_nh4_n15_sl, &   !< 15NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_no3_sl, &       !< NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_no3_n15_sl, &   !< 15NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                        mycorrhiza_uptake_norg_sl, &      !< N uptake by mycorrhiza from organic (residue) pool [mol / m3 / timestep]
                                                        mycorrhiza_uptake_norg_n15_sl, &  !< N uptake by mycorrhiza from organic (residue) pool [mol / m3 / timestep]
                                                        mycorrhiza_uptake_po4_sl, &       !< PO4 uptake by mycorrhiza [mol / m3 / timestep]
                                                        nh4_solute, &                     !< soil solute NH4 concentration [mol m-3]
                                                        no3_solute, &                     !< soil solute NO3 concentration [mol m-3]
                                                        po4_solute, &                     !< soil solute PO4 concentration [mol m-3]
                                                        km2_uptake_nh4_act, &             !< ..
                                                        km2_uptake_no3_act, &             !< ..
                                                        km2_uptake_po4_act, &             !< ..
                                                        km_uptake_org_act, &              !< ..
                                                        het_respiration, &                !< heterotrophic respiration [mol / m3 / timestep]
                                                        het_respiration_c13, &            !< heterotrophic respiration [mol / m3 / timestep]
                                                        het_respiration_c14, &            !< heterotrophic respiration [mol / m3 / timestep]
                                                        myc_respiration, &                !< mycorrhizae respiration [mol / m3 / timestep]
                                                        myc_respiration_c13, &            !< mycorrhizae respiration [mol / m3 / timestep]
                                                        myc_respiration_c14               !< mycorrhizae respiration [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_formation_mt(:,:,:,:)          !< bgcm sb_formation flux: formation rate  [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_loss_mt(:,:,:,:)               !< bgcm sb_loss flux: loss rate  [mol / m3 / timestep]
    REAL(wp),                          INTENT(inout) :: sb_mycorrhiza_export_mt(:,:,:)    !< bgcm sb_mycorrhiza_export flux: export from mycorrhiza to plants [mol / m2 / timestep]
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: transform_respiration, &
                                                        fresp, &
                                                        fturn_mycorrhiza, &
                                                        fvmax, &
                                                        vmax_norg_up, &
                                                        k_norg_up, &
                                                        norg_act, &
                                                        fnorg, &
                                                        f_myc_opt
    REAL(wp), DIMENSION(nc, nsoil_sb)                :: hlp1, hlp2, hlp3
    INTEGER                                          :: ielem               !< loop over bgcm elements
    INTEGER                                          :: ix_elem             !< index of element in bgcm, used for looping
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_mycorrhiza_dynamics'
    ! ----------------------------------------------------------------------------------------------------- !
    
    !>1.0 Mycorrhiza dynamics
    !>

    !>  1.1 uptake of nutrients by mycorrhiza and associated respiration losses to heterotrophic respiration
    !>
    WHERE(sb_pool_mt(ix_mycorrhiza, ixN, :, :) > eps8)
      hlp1(:,:) = MIN(1.0_wp, MAX(0.0_wp, &
        & (((sb_pool_mt(ix_mycorrhiza, ixC, :, :) / sb_pool_mt(ix_mycorrhiza, ixN, :, :)) &
        & - mycorrhiza_cn * f_norg_limit) &
        & / (mycorrhiza_cn * (1._wp - f_norg_limit)))))
      mycorrhiza_uptake_nh4_sl(:,:)     = mycorrhiza_uptake_nh4_sl(:,:)     * hlp1(:,:) 
      mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_n15_sl(:,:) * hlp1(:,:)  
      mycorrhiza_uptake_no3_sl(:,:)     = mycorrhiza_uptake_no3_sl(:,:)     * hlp1(:,:)
      mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_n15_sl(:,:) * hlp1(:,:)
    ENDWHERE

    WHERE(sb_pool_mt(ix_mycorrhiza, ixP, :, :) > eps8)
      hlp1(:,:) = MIN(1.0_wp,MAX(0.0_wp, &
        & (((sb_pool_mt(ix_mycorrhiza, ixN, :, :) / sb_pool_mt(ix_mycorrhiza, ixP, :, :)) &
        & - mycorrhiza_np * f_norg_limit) &
        & / (mycorrhiza_np * (1._wp - f_norg_limit)))))
      mycorrhiza_uptake_po4_sl(:,:) = mycorrhiza_uptake_po4_sl(:,:) * hlp1(:,:) 
    ENDWHERE

    ! calculate respiration costs associated with N transformation (mol C / m3 / timestep)
    transform_respiration(:,:) = mycorrhiza_uptake_nh4_sl(:,:) * transform_cost_nh4 &
      &                          + mycorrhiza_uptake_no3_sl(:,:) * transform_cost_no3

    fresp(:,:) = 0.0_wp
    WHERE(sb_pool_mt(ix_mycorrhiza, ixC, :, :) > eps8)
      fresp(:,:) = transform_respiration(:,:) / sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    ENDWHERE

    sb_loss_mt(ix_mycorrhiza, ixC, :, :)   = sb_loss_mt(ix_mycorrhiza, ixC, :, :)   + fresp(:,:) &
      &                                      * sb_pool_mt(ix_mycorrhiza, ixC, :, :)   !! respiration
    sb_loss_mt(ix_mycorrhiza, ixC13, :, :) = sb_loss_mt(ix_mycorrhiza, ixC13, :, :) + fresp(:,:) &
      &                                      * sb_pool_mt(ix_mycorrhiza, ixC13, :, :) !! respiration
    sb_loss_mt(ix_mycorrhiza, ixC14, :, :) = sb_loss_mt(ix_mycorrhiza, ixC14, :, :) + fresp(:,:) &
      &                                      * sb_pool_mt(ix_mycorrhiza, ixC14, :, :) !! respiration

    het_respiration(:,:)     = het_respiration(:,:)     + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

    myc_respiration(:,:)     = myc_respiration(:,:)     + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    myc_respiration_c13(:,:) = myc_respiration_c13(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC13, :, :)
    myc_respiration_c14(:,:) = myc_respiration_c14(:,:) + fresp(:,:) * sb_pool_mt(ix_mycorrhiza, ixC14, :, :)

    !>  1.2 organic N uptake by mycorrhiza and associated respiration losses to heterotrophic respiration
    !>
    IF (flag_mycorrhiza_org) THEN 
      ! limit organic uptake when inorganic nitrogen is abundant
      fvmax(:,:) = MIN(1.0_wp, MAX(0.0_wp, &
        & 0.5_wp * (2.0_wp - nh4_solute(:,:) / (f_nin2norg * km2_uptake_nh4_act(:,:) + nh4_solute(:,:)) &
        & - no3_solute(:,:) / (f_nin2norg * km2_uptake_no3_act(:,:) + no3_solute(:,:)))))
      hlp1(:,:) = fvmax(:,:)

      ! limit organic uptake according to accessibility to substrate
      fvmax(:,:) = fvmax(:,:) * sb_pool_mt(ix_residue, ixN, :, :) / ( km_uptake_org_act(:,:) + sb_pool_mt(ix_residue, ixN, :, :))
      mycorrhiza_uptake_norg_sl(:,:) = fvmax(:,:) * vmax_uptake_norg * sb_pool_mt(ix_mycorrhiza, ixC, :, :) &
        &                              * dtime / 1000000._wp

      WHERE(sb_pool_mt(ix_residue, ixN, :, :) > eps8) 
        mycorrhiza_uptake_norg_n15_sl(:,:) = mycorrhiza_uptake_norg_sl(:,:) &
          &                                  * sb_pool_mt(ix_residue, ixN15, :, :) / sb_pool_mt(ix_residue, ixN, :, :) 
      ENDWHERE

      ! reduce organic uptake at low mycorrhiza CN (should be obsolete with above function)
      WHERE(sb_pool_mt(ix_mycorrhiza, ixN, :, :) > eps8)
        hlp1(:,:) = MIN(1.0_wp, MAX(0.0_wp, &
          & (((sb_pool_mt(ix_mycorrhiza, ixC, :, :) / sb_pool_mt(ix_mycorrhiza, ixN, :, :)) &
          & - mycorrhiza_cn * f_norg_limit) &
          & / (mycorrhiza_cn * (1._wp - f_norg_limit)))))
        mycorrhiza_uptake_norg_sl(:,:)     = mycorrhiza_uptake_norg_sl(:,:)     * hlp1(:,:)
        mycorrhiza_uptake_norg_n15_sl(:,:) = mycorrhiza_uptake_norg_n15_sl(:,:) * hlp1(:,:) 
      ENDWHERE

      ! calculate N uptake fraction
      fnorg(:,:) = 0.0_wp
      WHERE(sb_pool_mt(ix_residue, ixN, :, :) > eps8)
        fnorg(:,:) = mycorrhiza_uptake_norg_sl(:,:) / sb_pool_mt(ix_residue, ixN, :, :)
      ENDWHERE
  
      sb_loss_mt(ix_residue, ixC, :, :)   = sb_loss_mt(ix_residue, ixC, :, :)   + fnorg(:,:) * sb_pool_mt(ix_residue, ixC, :, :)
      sb_loss_mt(ix_residue, ixN, :, :)   = sb_loss_mt(ix_residue, ixN, :, :)   + fnorg(:,:) * sb_pool_mt(ix_residue, ixN, :, :)
      sb_loss_mt(ix_residue, ixP, :, :)   = sb_loss_mt(ix_residue, ixP, :, :)   + fnorg(:,:) * sb_pool_mt(ix_residue, ixP, :, :)
      sb_loss_mt(ix_residue, ixC13, :, :) = sb_loss_mt(ix_residue, ixC13, :, :) + fnorg(:,:) * sb_pool_mt(ix_residue, ixC13, :, :)
      sb_loss_mt(ix_residue, ixC14, :, :) = sb_loss_mt(ix_residue, ixC14, :, :) + fnorg(:,:) * sb_pool_mt(ix_residue, ixC14, :, :)
      sb_loss_mt(ix_residue, ixN15, :, :) = sb_loss_mt(ix_residue, ixN15, :, :) + fnorg(:,:) * sb_pool_mt(ix_residue, ixN15, :, :)
 
      !new mycorrhiza formation out of uptake from slow N pool 
      sb_formation_mt(ix_mycorrhiza, ixC, :, :)   = sb_formation_mt(ix_mycorrhiza, ixC, :, :) &
        &                                           + mycorrhiza_cue * fnorg(:,:) * sb_pool_mt(ix_residue, ixC, :, :) 
      sb_formation_mt(ix_mycorrhiza, ixN, :, :)   = sb_formation_mt(ix_mycorrhiza, ixN, :, :) &
        &                                           + fnorg(:,:) * sb_pool_mt(ix_residue, ixN, :, :)
      sb_formation_mt(ix_mycorrhiza, ixP, :, :)   = sb_formation_mt(ix_mycorrhiza, ixP, :, :) &
        &                                           + fnorg(:,:) * sb_pool_mt(ix_residue, ixP, :, :)
      sb_formation_mt(ix_mycorrhiza, ixC13, :, :) = sb_formation_mt(ix_mycorrhiza, ixC13, :, :) &
        &                                           + mycorrhiza_cue * fnorg(:,:) * sb_pool_mt(ix_residue, ixC13, :, :) 
      sb_formation_mt(ix_mycorrhiza, ixC14, :, :) = sb_formation_mt(ix_mycorrhiza, ixC14, :, :) &
        &                                           + mycorrhiza_cue * fnorg(:,:) * sb_pool_mt(ix_residue, ixC14, :, :) 
      sb_formation_mt(ix_mycorrhiza, ixN15, :, :) = sb_formation_mt(ix_mycorrhiza, ixN15, :, :) &
        &                                           + fnorg(:,:) * sb_pool_mt(ix_residue, ixN15, :, :)
      
      het_respiration(:,:)     = het_respiration(:,:)     + (1.0_wp - mycorrhiza_cue) &
        &                        * fnorg(:,:) * sb_pool_mt(ix_residue, ixC, :, :)   
      het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1.0_wp - mycorrhiza_cue) &
        &                        * fnorg(:,:) * sb_pool_mt(ix_residue, ixC13, :, :) 
      het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1.0_wp - mycorrhiza_cue) & 
        &                        * fnorg(:,:) * sb_pool_mt(ix_residue, ixC14, :, :)
    ENDIF

    !>  1.3 nutrient flux to vegetation (surplus uptake goes to plants) [mol (N/P) / m3 / timestep]
    !>
    ! maintain C:N of mycorrhizal fungi
    ! N supply for plants = N uptake by mycorrhiza - N demand by mycorrhiza
    ! N demand can be negative => N supply can be greater then current N uptake by plants!

    ! N demand of mycorrhiza in mol N / m3 / timestep
    hlp1(:,:) = MAX(0.0_wp, & 
      &         sb_pool_mt(ix_mycorrhiza, ixC, :, :) / mycorrhiza_cn - sb_pool_mt(ix_mycorrhiza, ixN, :, :))
    ! N uptake incorporated into mycorrhiza (minimum of potential uptake and demand)
    hlp2(:,:) = MIN(hlp1(:,:), &
      &         mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) + mycorrhiza_uptake_norg_sl(:,:))
    ! remainder is exported to plants (mol N / m3 / timestep)
    sb_mycorrhiza_export_mt( ixN, :, :) = sb_mycorrhiza_export_mt( ixN, :, :) &
      &                                   + mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) &
      &                                   + mycorrhiza_uptake_norg_sl(:,:) - hlp2(:,:)
    WHERE(mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) + mycorrhiza_uptake_norg_sl(:,:) > eps8)
      sb_mycorrhiza_export_mt( ixN15, :, :) = sb_mycorrhiza_export_mt( ixN15, :, :) &
        &                                     + sb_mycorrhiza_export_mt( ixN, :, :) * (mycorrhiza_uptake_nh4_n15_sl(:,:) &
        &                                     + mycorrhiza_uptake_no3_n15_sl(:,:) &
        &                                     + mycorrhiza_uptake_norg_n15_sl(:,:)) / (mycorrhiza_uptake_nh4_sl(:,:) &
        &                                     + mycorrhiza_uptake_no3_sl(:,:) &
        &                                     + mycorrhiza_uptake_norg_sl(:,:))
    ENDWHERE
    ! assuming that mycorrhiza export organic N (C:N = 3:1)
    sb_mycorrhiza_export_mt( ixC, :, :) = sb_mycorrhiza_export_mt( ixC, :, :) &
      &                                   + sb_mycorrhiza_export_mt( ixN, :, :) * organic_cn
    WHERE(sb_pool_mt(ix_mycorrhiza, ixC, :, :) > eps8)
      sb_mycorrhiza_export_mt( ixC13, :, :) = sb_mycorrhiza_export_mt( ixC13, :, :) &
        &                                     + sb_mycorrhiza_export_mt( ixC, :, :) &
        &                                     * sb_pool_mt(ix_mycorrhiza, ixC13, :, :) / sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      sb_mycorrhiza_export_mt( ixC14, :, :) = sb_mycorrhiza_export_mt( ixC14, :, :) &
        &                                     + sb_mycorrhiza_export_mt( ixC, :, :) &
        &                                     * sb_pool_mt(ix_mycorrhiza, ixC14, :, :) / sb_pool_mt(ix_mycorrhiza, ixC, :, :)
    ENDWHERE

    ! P demand of mycorrhiza in mol P / m3 / timestep 
    hlp1(:,:) = MAX(sb_pool_mt(ix_mycorrhiza, ixN, :, :) / mycorrhiza_np &
      &         - sb_pool_mt(ix_mycorrhiza, ixP, :, :), 0.0_wp)
    ! incorporation of P into mycorrhiza (mol P / m3 / timestep), constrain flux to minimum of potential uptake and demand 
    hlp2(:,:) = MIN(hlp1(:,:), mycorrhiza_uptake_po4_sl(:,:))
    ! remainder is exported to plants (mol P / m3 / timestep)
    hlp3(:,:) = mycorrhiza_uptake_po4_sl(:,:) - hlp2(:,:)
    sb_mycorrhiza_export_mt( ixP, :, :) = sb_mycorrhiza_export_mt( ixP, :, :) + hlp3(:,:) 

    !>  1.4 Formation of litter from mycorrhiza turnover (mol / m3 / timestep)
    !>
    ! calc turnover 
    fturn_mycorrhiza(:,:) = 1.0_wp / lctlib_tau_mycorrhiza / one_day / one_year * dtime
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        ! ...
        sb_loss_mt(ix_mycorrhiza, ixC, :, :) = sb_loss_mt(ix_mycorrhiza, ixC, :, :) &
          &                                    + fturn_mycorrhiza(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)     
        ! ...
        sb_formation_mt(ix_polymeric_litter, ixC, :, :) = sb_formation_mt(ix_polymeric_litter, ixC, :, :) &
          &                                               + fturn_mycorrhiza(:,:) * sb_pool_mt(ix_mycorrhiza, ixC, :, :)
      END IF
    END DO
  END SUBROUTINE calc_mycorrhiza_dynamics

END MODULE mo_sb_simple_model
