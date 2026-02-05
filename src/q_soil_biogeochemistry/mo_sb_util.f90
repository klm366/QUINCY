!> QUINCY helper routines for soil-beogeochemistry
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
!>#### various helper routines for the soil-beogeochemistry process
!>
MODULE mo_sb_util
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: zero

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: reset_sb_fluxes, calculate_time_average_soilbiogeochemistry, update_read_soil_forcing

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
    USE mo_lnd_time_averages        ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    dsl4jsb_Use_processes SB_
    dsl4jsb_Use_memory(SB_)
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
    CHARACTER(len=*), PARAMETER       :: routine = modname//':calculate_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L3D :: sb_formation_mt
    dsl4jsb_Def_mt2L3D :: sb_loss_mt
    dsl4jsb_Def_mt1L3D :: sb_mycorrhiza_export_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
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
    IF (model%config%flag_slow_sb_pool_spinup_accelerator) THEN
      dsl4jsb_Get_mt2L3D(SB_BGCM_FORMATION_ID, sb_formation_mt)
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
  END SUBROUTINE calculate_time_average_soilbiogeochemistry

  ! ======================================================================================================= !
  !> read soil forcing from input files
  !>
#ifndef __NO_JSBACH__
#ifndef  __QUINCY_STANDALONE__
  SUBROUTINE update_read_soil_forcing(tile, options)

    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_process_class,       ONLY: SB_, HYDRO_, SSE_, VEG_, A2L_

    USE mo_jsb_io_netcdf,           ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx, read_bgcm_from_file
    USE mo_lnd_bgcm_idx
    USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID
    USE mo_jsb_time,                ONLY: t_datetime, get_date_components, get_year_day

    dsl4jsb_Use_processes SB_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model            !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store       !< the bgcm store of this tile
    TYPE(t_input_file)              :: input_file       !< file with to be read veg state variables
    TYPE(t_datetime),       POINTER :: current_datetime !< date and time at the beginning of this time step
    REAL(wp)                        :: dtime, year_instant
    INTEGER                         :: iblk    !< Dimensions of the horizontal grid
    INTEGER                         :: veg_pool_idx, litterfall_pool_idx
    INTEGER                         :: year, month, day, hour, minute, second, yearday, &
                                    &  ith_timestep     !< the row number of current timestep in soil forcing
    REAL(wp), POINTER               :: veg_pool_mt_domain(:,:,:,:)      !< dim: compartments, elements, nc, nblks -- on domain!
    REAL(wp), POINTER               :: veg_litterfall_mt_domain(:,:,:,:)      !< dim: compartments, elements, nc, nblks -- on domain!
    REAL(wp), POINTER :: &
      & ptr_1D(:    ),   &                      !< Pointer to 1d variable
      & ptr_2D(:,  :),   &                      !< Pointer to 2d variable
      & ptr_3D(:,:,:),   &                      !< Pointer to 3d variable
      & ptr_4D(:,:,:,:)                         !< Pointer to 4d variable

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_read_soil_forcing'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(SB_)
    ! A2L_
    dsl4jsb_Real2D_onDomain      :: nhx_deposition
    dsl4jsb_Real2D_onDomain      :: noy_deposition
    dsl4jsb_Real2D_onDomain      :: nhx_n15_deposition
    dsl4jsb_Real2D_onDomain      :: noy_n15_deposition
    dsl4jsb_Real2D_onDomain      :: p_deposition
    ! VEG_ 2D
    dsl4jsb_Real2D_onDomain      :: demand_uptake_n_tuptake_mavg
    dsl4jsb_Real2D_onDomain      :: demand_uptake_p_tuptake_mavg
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain      :: root_fraction_sl
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain      :: frac_wtr_transp_down_sl
    dsl4jsb_Real3D_onDomain      :: frac_w_lat_loss_sl
    dsl4jsb_Real3D_onDomain      :: wtr_soil_sl
    dsl4jsb_Real3D_onDomain      :: vol_porosity_sl
    ! SSE_ 3D
    dsl4jsb_Real3D_onDomain      :: t_soil_sl
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    iblk = options%iblk
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var2D_onDomain(A2L_,  nhx_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  noy_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  nhx_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  noy_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_,  p_deposition)
    dsl4jsb_Get_var2D_onDomain(VEG_,  demand_uptake_n_tuptake_mavg)
    dsl4jsb_Get_var2D_onDomain(VEG_,  demand_uptake_p_tuptake_mavg)
    dsl4jsb_Get_var3D_onDomain(VEG_,  root_fraction_sl)
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onDomain(HYDRO_, frac_wtr_transp_down_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, frac_w_lat_loss_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    ! SSE_ 3D
    dsl4jsb_Get_var3D_onDomain(SSE_, t_soil_sl)

    !> get the timestep information of current timestep
    current_datetime = options%current_datetime
    dtime   = options%dtime
    CALL get_date_components(current_datetime, year, month, day, hour, minute, second)
    year_instant = get_year_day(current_datetime)
    yearday = INT(year_instant)
    ith_timestep = IMerge_TIME2Sec(yearday, hour, minute, second)/int(dtime) + 1 ! the row number of current timestep in forcing file

    !> 1.0 the first soil forcing file
    ! ----------------------------------------------------------------------------------------------------- !
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%soilphysics_forcing_filename), model%grid_id)

    !> read the 3D forcing variables
    !> - Fraction of fineroot
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='veg_root_fraction_sl_'//tile%name,          &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    root_fraction_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Fraction of vertical transport
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_frac_wtr_transp_down_sl_'//tile%name,   &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    frac_wtr_transp_down_sl(:,:,:) = ptr_4D(:,:,:,1)

    !> - Fraction of lateral water loss
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_frac_w_lat_loss_sl_'//tile%name,        &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    frac_w_lat_loss_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Water content
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_wtr_soil_sl_'//tile%name,              &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    wtr_soil_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Volumetric water content
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_vol_porosity_sl_'//tile%name,          &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    vol_porosity_sl(:,:,:) = MERGE(ptr_4D(:,:,:,1), 0._wp, ptr_4D(:,:,:,1) >= 0._wp)

    !> - Soil temperature
    ptr_4D => input_file%Read_3d_time(                  &
        & variable_name='sse_t_soil_sl_'//tile%name,                &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    t_soil_sl(:,:,:) = ptr_4D(:,:,:,1)

    !> read the 2D forcing variables
    !> - nhx deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_nhx_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    nhx_deposition(:,:) = ptr_3D(:,:,1)

    !> - noy deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_noy_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    noy_deposition(:,:) = ptr_3D(:,:,1)

    !> - nhx n15 deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_nhx_n15_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    nhx_n15_deposition(:,:) = ptr_3D(:,:,1)

    !> - noy n15 deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_noy_n15_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    noy_n15_deposition(:,:) = ptr_3D(:,:,1)

    !> - p deposition
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="a2l_p_deposition_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    p_deposition(:,:) = ptr_3D(:,:,1)

    !> - plant N uptake demand
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_demand_uptake_n_tuptake_mavg_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    demand_uptake_n_tuptake_mavg(:,:) = ptr_3D(:,:,1)

    !> - plant P uptake demand
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_demand_uptake_p_tuptake_mavg_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    demand_uptake_p_tuptake_mavg(:,:) = ptr_3D(:,:,1)

    CALL input_file%Close()

    !> 2.0 the second soil forcing file
    ! ----------------------------------------------------------------------------------------------------- !
    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%veglitterfall_forcing_filename), model%grid_id)
    IF (ASSOCIATED(tile%bgcm_store)) THEN
      bgcm_store => tile%bgcm_store
      ! litterfall_pool_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile%name, routine)
      ! veg_litterfall_mt_domain => bgcm_store%store_2l_2d_bgcms(bgcm_store%idx_in_store(litterfall_pool_idx))%mt_2l_2d_bgcm(:,:,:,:)
      veg_pool_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_POOL_ID, tile%name, routine)
      veg_pool_mt_domain => bgcm_store%store_2l_2d_bgcms(bgcm_store%idx_in_store(veg_pool_idx))%mt_2l_2d_bgcm(:,:,:,:)
    END IF

    !> - veg fineroot C
    ptr_3D => input_file%Read_2d_time(                  &
        & variable_name="veg_pool_fine_root_C_"//tile%name,            &
        & start_time_step=ith_timestep, end_time_step=ith_timestep)
    veg_pool_mt_domain(ix_fine_root, ixC, :,:) = ptr_3D(:,:,1)

    !> - veg_litterfall_mt_domain using read_bgcm_from_file
    CALL read_bgcm_from_file(bgcm_store, VEG_BGCM_LITTERFALL_ID, "veg", tile%name, routine, input_file, dt_step = ith_timestep)
    ! !> - litterfall leaf
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_leaf_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_leaf, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_leaf_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_leaf, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_leaf_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_leaf, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_leaf_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_leaf, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_leaf_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_leaf, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_leaf_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_leaf, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall fine root
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fine_root_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fine_root, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fine_root_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fine_root, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fine_root_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fine_root, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fine_root_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fine_root, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fine_root_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fine_root, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fine_root_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fine_root, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall coarse root
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_coarse_root_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_coarse_root, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_coarse_root_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_coarse_root, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_coarse_root_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_coarse_root, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_coarse_root_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_coarse_root, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_coarse_root_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_coarse_root, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_coarse_root_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_coarse_root, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall sap wood
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_sap_wood_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_sap_wood, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_sap_wood_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_sap_wood, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_sap_wood_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_sap_wood, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_sap_wood_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_sap_wood, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_sap_wood_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_sap_wood, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_sap_wood_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_sap_wood, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall heart wood
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_heart_wood_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_heart_wood, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_heart_wood_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_heart_wood, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_heart_wood_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_heart_wood, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_heart_wood_C13_"//tile%name,           &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_heart_wood, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_heart_wood_C14_"//tile%name,           &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_heart_wood, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_heart_wood_N15_"//tile%name,           &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_heart_wood, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall fruit
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fruit_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fruit, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fruit_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fruit, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fruit_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fruit, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fruit_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fruit, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fruit_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fruit, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_fruit_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_fruit, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall labile pool
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_labile_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_labile, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_labile_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_labile, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_labile_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_labile, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_labile_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_labile, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_labile_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_labile, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_labile_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_labile, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall reserve pool
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_reserve_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_reserve, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_reserve_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_reserve, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_reserve_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_reserve, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_reserve_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_reserve, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_reserve_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_reserve, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_reserve_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_reserve, ixN15, :,:) = ptr_3D(:,:,1)

    ! !> - litterfall seed bed
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_seed_bed_C_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_seed_bed, ixC, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_seed_bed_N_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_seed_bed, ixN, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_seed_bed_P_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_seed_bed, ixP, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_seed_bed_C13_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_seed_bed, ixC13, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_seed_bed_C14_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_seed_bed, ixC14, :,:) = ptr_3D(:,:,1)
    ! ptr_3D => input_file%Read_2d_time(                  &
        ! & variable_name="veg_litterfall_seed_bed_N15_"//tile%name,            &
        ! & start_time_step=ith_timestep, end_time_step=ith_timestep)
    ! veg_litterfall_mt_domain(ix_seed_bed, ixN15, :,:) = ptr_3D(:,:,1)

    CALL input_file%Close()

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_read_soil_forcing

  ! ======================================================================================================= !
  !>
  !> @TODO add docu
  !>
  FUNCTION IMerge_TIME2Sec (kd, kh, km, ks) RESULT (ix)
    !+
    !
    ! IMerge_TIME2Sec [function, integer]
    !    merge day, hour, minute and second to internal time
    !    (
    !    day    [integer] input (day of year, 0...364)
    !    hour   [integer] input (hour of day, 0...23)
    !    minute [integer] input (minute of hour, 0...59)
    !    second [integer] input (second of minute, 0...59)
    !    )
    !
    !-
    INTEGER, INTENT(in) :: kd, kh, km, ks
    INTEGER             :: ix

    ix = kd*86400 + kh*3600 + km*60 + ks

  END FUNCTION IMerge_TIME2Sec
#endif
#endif

#endif
END MODULE mo_sb_util
