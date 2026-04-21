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
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, finish
  USE mo_jsb_math_constants,    ONLY: zero
  USE mo_util,                  ONLY: int2string

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID, &
                                  VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_SEED_BED_LITTERFALL_ID
#ifndef  __QUINCY_STANDALONE__
  USE mo_jsb_io_netcdf,           ONLY: t_input_file
  USE mo_read_interface,          ONLY: read_2D_extdim, read_3D_extdim, on_cells
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: reset_sb_fluxes, calculate_time_average_soilbiogeochemistry
#ifndef  __QUINCY_STANDALONE__
  PUBLIC :: update_read_soil_forcing
#endif

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_util'

#ifndef  __QUINCY_STANDALONE__
  TYPE t_soil_forcing_cache
    TYPE(t_input_file) :: soilphysics_input_file
    TYPE(t_input_file) :: veglitterfall_input_file
    LOGICAL :: soilphysics_open = .FALSE.
    LOGICAL :: veglitterfall_open = .FALSE.
    LOGICAL :: arrays_initialized = .FALSE.
    LOGICAL :: preload_loaded = .FALSE.
    INTEGER :: nproma = 0
    INTEGER :: nsoil = 0
    INTEGER :: nblks = 0
    INTEGER :: npft = 0
    INTEGER :: veg_litter_nr_parts = 0
    INTEGER :: veg_litter_nr_elems = 0
    INTEGER :: seed_bed_nr_elems = 0
    INTEGER :: preload_nsteps = 0
    INTEGER :: cached_year = -HUGE(0)
    INTEGER :: cached_yearday = -HUGE(0)
    INTEGER :: cached_timestep = -HUGE(0)
    CHARACTER(len=:), ALLOCATABLE :: soilphysics_filename
    CHARACTER(len=:), ALLOCATABLE :: veglitterfall_filename
    REAL(wp), ALLOCATABLE :: root_fraction_sl_preload(:,:,:,:,:)      ! nproma, nsoil, nblks, npft, ntime
    REAL(wp), ALLOCATABLE :: demand_uptake_n_preload(:,:,:,:)         ! nproma, nblks, npft, ntime
    REAL(wp), ALLOCATABLE :: demand_uptake_p_preload(:,:,:,:)         ! nproma, nblks, npft, ntime
    REAL(wp), ALLOCATABLE :: frac_wtr_transp_down_sl_preload(:,:,:,:) ! nproma, nsoil, nblks, ntime
    REAL(wp), ALLOCATABLE :: frac_w_lat_loss_sl_preload(:,:,:,:)      ! nproma, nsoil, nblks, ntime
    REAL(wp), ALLOCATABLE :: wtr_soil_sl_preload(:,:,:,:)             ! nproma, nsoil, nblks, ntime
    REAL(wp), ALLOCATABLE :: vol_porosity_sl_preload(:,:,:,:)         ! nproma, nsoil, nblks, ntime
    REAL(wp), ALLOCATABLE :: t_soil_sl_preload(:,:,:,:)               ! nproma, nsoil, nblks, ntime
    REAL(wp), ALLOCATABLE :: nhx_deposition_preload(:,:,:)            ! nproma, nblks, ntime
    REAL(wp), ALLOCATABLE :: noy_deposition_preload(:,:,:)            ! nproma, nblks, ntime
    REAL(wp), ALLOCATABLE :: nhx_n15_deposition_preload(:,:,:)        ! nproma, nblks, ntime
    REAL(wp), ALLOCATABLE :: noy_n15_deposition_preload(:,:,:)        ! nproma, nblks, ntime
    REAL(wp), ALLOCATABLE :: p_deposition_preload(:,:,:)              ! nproma, nblks, ntime
    REAL(wp), ALLOCATABLE :: veg_pool_fine_root_c_preload(:,:,:,:)    ! nproma, nblks, npft, ntime
    REAL(wp), ALLOCATABLE :: veg_litterfall_preload(:,:,:,:,:,:)      ! part, elem, nproma, nblks, npft, ntime
    REAL(wp), ALLOCATABLE :: seed_bed_litterfall_preload(:,:,:,:,:)   ! elem, nproma, nblks, npft, ntime
    REAL(wp), ALLOCATABLE :: scratch_3d(:,:,:)                        ! nproma, nblks, 1
    REAL(wp), ALLOCATABLE :: scratch_4d(:,:,:,:)                      ! nproma, nsoil, nblks, 1
  END TYPE t_soil_forcing_cache

  TYPE(t_soil_forcing_cache), SAVE :: soil_forcing_cache
#endif

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
  !> read soil forcing from input files (ONLY for IQ)
  !>
#ifndef  __QUINCY_STANDALONE__
  SUBROUTINE update_read_soil_forcing(tile, options)

    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_process_class,       ONLY: SB_, HYDRO_, SSE_, VEG_, A2L_

    USE mo_lnd_bgcm_idx
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
    REAL(wp)                        :: dtime, year_instant
    INTEGER                         :: iblk             !< Dimensions of the horizontal grid
    INTEGER                         :: pft_idx, npft
    INTEGER                         :: year, month, day, hour, minute, second, yearday, forcing_yearday, &
                                    &  ith_timestep, cache_itimestep     !< forcing record index in file / cache
    CHARACTER(len=:), ALLOCATABLE   :: soilphysics_forcing_file, veglitterfall_forcing_file

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
    dtime   = options%dtime
    CALL get_date_components(options%current_datetime, year, month, day, hour, minute, second)
    year_instant = get_year_day(options%current_datetime)
    yearday = INT(year_instant)
    forcing_yearday = MODULO(yearday - 1, 365) + 1
    CALL get_soil_forcing_info(TRIM(dsl4jsb_Config(SB_)%soilphysics_forcing_filename), &
      & forcing_yearday, hour, minute, second, dtime, soilphysics_forcing_file, ith_timestep)
    CALL get_soil_forcing_info(TRIM(dsl4jsb_Config(SB_)%veglitterfall_forcing_filename), &
      & forcing_yearday, hour, minute, second, dtime, veglitterfall_forcing_file, ith_timestep)
    IF (debug_on()) CALL message(TRIM(routine), 'Reading the forcing record '//int2string(ith_timestep)// &
      & ' at forcing day '//int2string(forcing_yearday)//' (model day '//int2string(yearday)//') from ' // &
      & TRIM(soilphysics_forcing_file)//' ...')
    bgcm_store => tile%bgcm_store
    npft = model%config%npft
    pft_idx = get_pft_tile_index(tile%name)
    CALL ensure_soil_forcing_cache_loaded(model, bgcm_store, root_fraction_sl, soilphysics_forcing_file, &
      & veglitterfall_forcing_file, year, forcing_yearday, ith_timestep, dtime, npft, routine)
    cache_itimestep = ith_timestep - soil_forcing_cache%cached_timestep + 1
    IF (cache_itimestep < 1 .OR. cache_itimestep > soil_forcing_cache%preload_nsteps) THEN
      CALL finish(TRIM(routine), 'Invalid soil forcing cache index '//int2string(cache_itimestep)// &
        & ' for forcing timestep '//int2string(ith_timestep)//' with cache start '// &
        & int2string(soil_forcing_cache%cached_timestep)//' and preload_nsteps '// &
        & int2string(soil_forcing_cache%preload_nsteps))
    END IF

    IF (pft_idx > 0) THEN
      root_fraction_sl(:,:,:) = soil_forcing_cache%root_fraction_sl_preload(:,:,:,pft_idx,cache_itimestep)
      demand_uptake_n_tuptake_mavg(:,:) = soil_forcing_cache%demand_uptake_n_preload(:,:,pft_idx,cache_itimestep)
      demand_uptake_p_tuptake_mavg(:,:) = soil_forcing_cache%demand_uptake_p_preload(:,:,pft_idx,cache_itimestep)
      CALL apply_litterfall_cache_to_tile(bgcm_store, tile%name, pft_idx, cache_itimestep, routine)
    ELSE
      frac_wtr_transp_down_sl(:,:,:) = soil_forcing_cache%frac_wtr_transp_down_sl_preload(:,:,:,cache_itimestep)
      frac_w_lat_loss_sl(:,:,:)      = soil_forcing_cache%frac_w_lat_loss_sl_preload(:,:,:,cache_itimestep)
      wtr_soil_sl(:,:,:)             = soil_forcing_cache%wtr_soil_sl_preload(:,:,:,cache_itimestep)
      vol_porosity_sl(:,:,:)         = soil_forcing_cache%vol_porosity_sl_preload(:,:,:,cache_itimestep)
      t_soil_sl(:,:,:)               = soil_forcing_cache%t_soil_sl_preload(:,:,:,cache_itimestep)
      nhx_deposition(:,:)            = soil_forcing_cache%nhx_deposition_preload(:,:,cache_itimestep)
      noy_deposition(:,:)            = soil_forcing_cache%noy_deposition_preload(:,:,cache_itimestep)
      nhx_n15_deposition(:,:)        = soil_forcing_cache%nhx_n15_deposition_preload(:,:,cache_itimestep)
      noy_n15_deposition(:,:)        = soil_forcing_cache%noy_n15_deposition_preload(:,:,cache_itimestep)
      p_deposition(:,:)              = soil_forcing_cache%p_deposition_preload(:,:,cache_itimestep)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_read_soil_forcing

  ! ======================================================================================================= !
  !> Return the number of records to preload for the current forcing file.
  !>
  INTEGER FUNCTION get_soil_forcing_preload_nsteps(dtime)

    REAL(wp), INTENT(in) :: dtime

    IF (NINT(dtime) == 86400) THEN
      get_soil_forcing_preload_nsteps = 1
    ELSE
      get_soil_forcing_preload_nsteps = 86400 / INT(dtime)
    END IF

  END FUNCTION get_soil_forcing_preload_nsteps

  ! ======================================================================================================= !
  ! ======================================================================================================= !
  !> Return the pft index encoded in a tile name like pft01.
  !>
  INTEGER FUNCTION get_pft_tile_index(tile_name)

    CHARACTER(len=*), INTENT(in) :: tile_name
    INTEGER :: ios

    get_pft_tile_index = 0
    IF (LEN_TRIM(tile_name) /= 5) RETURN
    IF (tile_name(1:3) /= 'pft') RETURN

    READ(tile_name(4:5), '(I2)', IOSTAT=ios) get_pft_tile_index
    IF (ios /= 0) get_pft_tile_index = 0

  END FUNCTION get_pft_tile_index

  ! ======================================================================================================= !
  !> Return the canonical pft tile name.
  !>
  FUNCTION get_pft_tile_name(ipft) RESULT(tile_name)

    INTEGER, INTENT(in) :: ipft
    CHARACTER(len=5)    :: tile_name

    WRITE(tile_name, '("pft",I2.2)') ipft

  END FUNCTION get_pft_tile_name

  ! ======================================================================================================= !
  !> Return the litterfall compartment name used in the forcing variable names.
  !>
  FUNCTION get_veg_litter_part_name(i_part) RESULT(part_name)

    INTEGER, INTENT(in)    :: i_part
    CHARACTER(len=:), ALLOCATABLE :: part_name

    SELECT CASE(i_part)
    CASE (1)
      part_name = 'leaf'
    CASE (2)
      part_name = 'fine_root'
    CASE (3)
      part_name = 'coarse_root'
    CASE (4)
      part_name = 'sap_wood'
    CASE (5)
      part_name = 'heart_wood'
    CASE (6)
      part_name = 'labile'
    CASE (7)
      part_name = 'reserve'
    CASE (8)
      part_name = 'fruit'
    CASE DEFAULT
      CALL finish(modname//':get_veg_litter_part_name', 'Unsupported litterfall part index '//int2string(i_part))
    END SELECT

  END FUNCTION get_veg_litter_part_name

  ! ======================================================================================================= !
  !> Return the litterfall element name used in the forcing variable names.
  !>
  FUNCTION get_veg_litter_elem_name(i_elem) RESULT(elem_name)

    INTEGER, INTENT(in)    :: i_elem
    CHARACTER(len=:), ALLOCATABLE :: elem_name

    SELECT CASE(i_elem)
    CASE (1)
      elem_name = 'C'
    CASE (2)
      elem_name = 'N'
    CASE (3)
      elem_name = 'P'
    CASE (4)
      elem_name = 'C13'
    CASE (5)
      elem_name = 'C14'
    CASE (6)
      elem_name = 'N15'
    CASE DEFAULT
      CALL finish(modname//':get_veg_litter_elem_name', 'Unsupported litterfall element index '//int2string(i_elem))
    END SELECT

  END FUNCTION get_veg_litter_elem_name

  ! ======================================================================================================= !
  !> Keep the forcing files open and cache the currently loaded timestep.
  !>
  SUBROUTINE ensure_soil_forcing_cache_loaded(model, bgcm_store, root_fraction_sl, soilphysics_filename, &
    & veglitterfall_filename, year, yearday, itimestep, dtime, npft, caller)

    USE mo_jsb_model_class, ONLY: t_jsb_model

    TYPE(t_jsb_model),       POINTER, INTENT(in) :: model
    TYPE(t_lnd_bgcm_store),  POINTER, INTENT(in) :: bgcm_store
    REAL(wp),                        INTENT(in) :: root_fraction_sl(:,:,:)
    CHARACTER(len=*),                INTENT(in) :: soilphysics_filename, veglitterfall_filename, caller
    INTEGER,                         INTENT(in) :: year, yearday, itimestep, npft
    REAL(wp),                        INTENT(in) :: dtime

    LOGICAL :: cache_hit
    INTEGER :: preload_nsteps, cache_start_timestep

    preload_nsteps = get_soil_forcing_preload_nsteps(dtime)
    IF (preload_nsteps == 1) THEN
      cache_start_timestep = itimestep
    ELSE
      cache_start_timestep = 1
    END IF
    !$OMP CRITICAL (soil_forcing_netcdf_io)
    CALL ensure_soil_forcing_cache_shape(bgcm_store, root_fraction_sl, npft, preload_nsteps, caller)
    CALL open_cached_soil_forcing_file(model, soil_forcing_cache%soilphysics_input_file, &
      & soil_forcing_cache%soilphysics_open, soil_forcing_cache%soilphysics_filename, soilphysics_filename, caller)
    CALL open_cached_soil_forcing_file(model, soil_forcing_cache%veglitterfall_input_file, &
      & soil_forcing_cache%veglitterfall_open, soil_forcing_cache%veglitterfall_filename, veglitterfall_filename, caller)

    cache_hit = .FALSE.
    IF (ALLOCATED(soil_forcing_cache%soilphysics_filename) .AND. ALLOCATED(soil_forcing_cache%veglitterfall_filename)) THEN
      cache_hit = soil_forcing_cache%preload_loaded .AND. soil_forcing_cache%preload_nsteps == preload_nsteps .AND. &
        & TRIM(soil_forcing_cache%soilphysics_filename) == TRIM(soilphysics_filename) .AND. &
        & TRIM(soil_forcing_cache%veglitterfall_filename) == TRIM(veglitterfall_filename) .AND. &
        & soil_forcing_cache%cached_timestep == cache_start_timestep
    END IF

    IF (.NOT. cache_hit) THEN
      CALL load_soilphysics_cache(bgcm_store, preload_nsteps, cache_start_timestep, caller)
      CALL load_litterfall_cache(bgcm_store, preload_nsteps, cache_start_timestep, caller)
      soil_forcing_cache%preload_loaded = .TRUE.
      soil_forcing_cache%cached_year = year
      soil_forcing_cache%cached_yearday = yearday
      soil_forcing_cache%cached_timestep = cache_start_timestep
    END IF
    !$OMP END CRITICAL (soil_forcing_netcdf_io)

  END SUBROUTINE ensure_soil_forcing_cache_loaded

  ! ======================================================================================================= !
  !> Ensure that the cache arrays match the current model dimensions.
  !>
  SUBROUTINE ensure_soil_forcing_cache_shape(bgcm_store, root_fraction_sl, npft, preload_nsteps, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(in) :: bgcm_store
    REAL(wp),                      INTENT(in) :: root_fraction_sl(:,:,:)
    INTEGER,                       INTENT(in) :: npft
    INTEGER,                       INTENT(in) :: preload_nsteps
    CHARACTER(len=*),              INTENT(in) :: caller

    INTEGER :: nproma, nsoil, nblks
    INTEGER :: veg_litter_idx, seed_bed_idx
    INTEGER :: veg_litter_store_idx, seed_bed_store_idx
    LOGICAL :: needs_realloc
    CHARACTER(len=5) :: tile_name

    nproma = SIZE(root_fraction_sl, 1)
    nsoil  = SIZE(root_fraction_sl, 2)
    nblks  = SIZE(root_fraction_sl, 3)
    tile_name = get_pft_tile_name(1)

    veg_litter_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile_name, caller)
    seed_bed_idx   = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_LITTERFALL_ID, tile_name, caller)
    veg_litter_store_idx = bgcm_store%idx_in_store(veg_litter_idx)
    seed_bed_store_idx   = bgcm_store%idx_in_store(seed_bed_idx)

    needs_realloc = .NOT. soil_forcing_cache%arrays_initialized
    IF (.NOT. needs_realloc) THEN
      needs_realloc = soil_forcing_cache%nproma /= nproma .OR. soil_forcing_cache%nsoil /= nsoil .OR. &
        & soil_forcing_cache%nblks /= nblks .OR. soil_forcing_cache%npft /= npft .OR. &
        & soil_forcing_cache%preload_nsteps /= preload_nsteps .OR. &
        & soil_forcing_cache%veg_litter_nr_parts /= bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_parts .OR. &
        & soil_forcing_cache%veg_litter_nr_elems /= bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_elems .OR. &
        & soil_forcing_cache%seed_bed_nr_elems /= bgcm_store%store_1l_2d_bgcms(seed_bed_store_idx)%nr_of_elems
    END IF

    IF (.NOT. needs_realloc) RETURN

    CALL reset_soil_forcing_cache_arrays()

    soil_forcing_cache%nproma = nproma
    soil_forcing_cache%nsoil  = nsoil
    soil_forcing_cache%nblks  = nblks
    soil_forcing_cache%npft   = npft
    soil_forcing_cache%preload_nsteps = preload_nsteps
    soil_forcing_cache%veg_litter_nr_parts = bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_parts
    soil_forcing_cache%veg_litter_nr_elems = bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%nr_of_elems
    soil_forcing_cache%seed_bed_nr_elems   = bgcm_store%store_1l_2d_bgcms(seed_bed_store_idx)%nr_of_elems
    ALLOCATE(soil_forcing_cache%root_fraction_sl_preload(nproma, nsoil, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%demand_uptake_n_preload(nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%demand_uptake_p_preload(nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%frac_wtr_transp_down_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%frac_w_lat_loss_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%wtr_soil_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%vol_porosity_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%t_soil_sl_preload(nproma, nsoil, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%nhx_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%noy_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%nhx_n15_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%noy_n15_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%p_deposition_preload(nproma, nblks, preload_nsteps))
    ALLOCATE(soil_forcing_cache%veg_pool_fine_root_c_preload(nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%veg_litterfall_preload(soil_forcing_cache%veg_litter_nr_parts, soil_forcing_cache%veg_litter_nr_elems, &
      & nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%seed_bed_litterfall_preload(soil_forcing_cache%seed_bed_nr_elems, nproma, nblks, npft, preload_nsteps))
    ALLOCATE(soil_forcing_cache%scratch_3d(nproma, nblks, 1))
    ALLOCATE(soil_forcing_cache%scratch_4d(nproma, nsoil, nblks, 1))

    soil_forcing_cache%root_fraction_sl_preload = zero
    soil_forcing_cache%demand_uptake_n_preload = zero
    soil_forcing_cache%demand_uptake_p_preload = zero
    soil_forcing_cache%frac_wtr_transp_down_sl_preload = zero
    soil_forcing_cache%frac_w_lat_loss_sl_preload = zero
    soil_forcing_cache%wtr_soil_sl_preload = zero
    soil_forcing_cache%vol_porosity_sl_preload = zero
    soil_forcing_cache%t_soil_sl_preload = zero
    soil_forcing_cache%nhx_deposition_preload = zero
    soil_forcing_cache%noy_deposition_preload = zero
    soil_forcing_cache%nhx_n15_deposition_preload = zero
    soil_forcing_cache%noy_n15_deposition_preload = zero
    soil_forcing_cache%p_deposition_preload = zero
    soil_forcing_cache%veg_pool_fine_root_c_preload = zero
    soil_forcing_cache%veg_litterfall_preload = zero
    soil_forcing_cache%seed_bed_litterfall_preload = zero
    soil_forcing_cache%scratch_3d = zero
    soil_forcing_cache%scratch_4d = zero

    soil_forcing_cache%arrays_initialized = .TRUE.
    soil_forcing_cache%preload_loaded = .FALSE.
    soil_forcing_cache%cached_year = -HUGE(0)
    soil_forcing_cache%cached_yearday = -HUGE(0)
    soil_forcing_cache%cached_timestep = -HUGE(0)
  END SUBROUTINE ensure_soil_forcing_cache_shape

  ! ======================================================================================================= !
  !> Deallocate the cache arrays while keeping any open files untouched.
  !>
  SUBROUTINE reset_soil_forcing_cache_arrays()

    IF (ALLOCATED(soil_forcing_cache%root_fraction_sl_preload)) DEALLOCATE(soil_forcing_cache%root_fraction_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%demand_uptake_n_preload)) DEALLOCATE(soil_forcing_cache%demand_uptake_n_preload)
    IF (ALLOCATED(soil_forcing_cache%demand_uptake_p_preload)) DEALLOCATE(soil_forcing_cache%demand_uptake_p_preload)
    IF (ALLOCATED(soil_forcing_cache%frac_wtr_transp_down_sl_preload)) DEALLOCATE(soil_forcing_cache%frac_wtr_transp_down_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%frac_w_lat_loss_sl_preload)) DEALLOCATE(soil_forcing_cache%frac_w_lat_loss_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%wtr_soil_sl_preload)) DEALLOCATE(soil_forcing_cache%wtr_soil_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%vol_porosity_sl_preload)) DEALLOCATE(soil_forcing_cache%vol_porosity_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%t_soil_sl_preload)) DEALLOCATE(soil_forcing_cache%t_soil_sl_preload)
    IF (ALLOCATED(soil_forcing_cache%nhx_deposition_preload)) DEALLOCATE(soil_forcing_cache%nhx_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%noy_deposition_preload)) DEALLOCATE(soil_forcing_cache%noy_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%nhx_n15_deposition_preload)) DEALLOCATE(soil_forcing_cache%nhx_n15_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%noy_n15_deposition_preload)) DEALLOCATE(soil_forcing_cache%noy_n15_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%p_deposition_preload)) DEALLOCATE(soil_forcing_cache%p_deposition_preload)
    IF (ALLOCATED(soil_forcing_cache%veg_pool_fine_root_c_preload)) DEALLOCATE(soil_forcing_cache%veg_pool_fine_root_c_preload)
    IF (ALLOCATED(soil_forcing_cache%veg_litterfall_preload)) DEALLOCATE(soil_forcing_cache%veg_litterfall_preload)
    IF (ALLOCATED(soil_forcing_cache%seed_bed_litterfall_preload)) DEALLOCATE(soil_forcing_cache%seed_bed_litterfall_preload)
    IF (ALLOCATED(soil_forcing_cache%scratch_3d)) DEALLOCATE(soil_forcing_cache%scratch_3d)
    IF (ALLOCATED(soil_forcing_cache%scratch_4d)) DEALLOCATE(soil_forcing_cache%scratch_4d)

    soil_forcing_cache%arrays_initialized = .FALSE.
    soil_forcing_cache%preload_loaded = .FALSE.
    soil_forcing_cache%nproma = 0
    soil_forcing_cache%nsoil = 0
    soil_forcing_cache%nblks = 0
    soil_forcing_cache%npft = 0
    soil_forcing_cache%preload_nsteps = 0
    soil_forcing_cache%veg_litter_nr_parts = 0
    soil_forcing_cache%veg_litter_nr_elems = 0
    soil_forcing_cache%seed_bed_nr_elems = 0
    soil_forcing_cache%cached_year = -HUGE(0)
    soil_forcing_cache%cached_yearday = -HUGE(0)
    soil_forcing_cache%cached_timestep = -HUGE(0)

  END SUBROUTINE reset_soil_forcing_cache_arrays

  ! ======================================================================================================= !
  !> Reuse an open forcing file when possible, otherwise reopen it.
  !>
  SUBROUTINE open_cached_soil_forcing_file(model, input_file, is_open, cached_filename, target_filename, caller)

    USE mo_jsb_model_class, ONLY: t_jsb_model
    USE mo_jsb_io_netcdf,   ONLY: jsb_netcdf_open_input
    USE mo_jsb_parallel,    ONLY: my_process_is_stdio

    TYPE(t_jsb_model), POINTER, INTENT(in) :: model
    TYPE(t_input_file),        INTENT(inout) :: input_file
    LOGICAL,                   INTENT(inout) :: is_open
    CHARACTER(len=:), ALLOCATABLE, INTENT(inout) :: cached_filename
    CHARACTER(len=*),          INTENT(in) :: target_filename, caller

    IF (is_open .AND. ALLOCATED(cached_filename)) THEN
      IF (TRIM(cached_filename) == TRIM(target_filename)) THEN
        RETURN
      END IF
      CALL input_file%Close()
      is_open = .FALSE.
    END IF

    input_file = jsb_netcdf_open_input(TRIM(target_filename), model%grid_id)

    ! The standard JSBACH input wrapper can report is_open=.FALSE. on non-stdio ranks
    ! even though the distributed read path is still valid. Only the stdio rank should
    ! treat this as a hard open failure.
    IF (my_process_is_stdio() .AND. .NOT. input_file%is_open) THEN
      CALL finish(TRIM(caller), 'Failed to open soil forcing file '//TRIM(target_filename))
    END IF

    is_open = .TRUE.
    cached_filename = TRIM(target_filename)
    soil_forcing_cache%preload_loaded = .FALSE.
    soil_forcing_cache%cached_year = -HUGE(0)
    soil_forcing_cache%cached_yearday = -HUGE(0)
    soil_forcing_cache%cached_timestep = -HUGE(0)

  END SUBROUTINE open_cached_soil_forcing_file

  ! ======================================================================================================= !
  !> Read and cache the forcing fields from bc_qsoil_physics.
  !>
  SUBROUTINE load_soilphysics_cache(bgcm_store, nsteps, start_timestep, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(in) :: bgcm_store
    INTEGER,                       INTENT(in) :: nsteps
    INTEGER,                       INTENT(in) :: start_timestep
    CHARACTER(len=*),              INTENT(in) :: caller

    INTEGER :: ipft, itimestep, source_timestep
    CHARACTER(len=5) :: tile_name

    DO itimestep = 1, nsteps
      source_timestep = start_timestep + itimestep - 1
      DO ipft = 1, soil_forcing_cache%npft
        tile_name = get_pft_tile_name(ipft)

        CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'veg_root_fraction_sl_'//tile_name, soil_forcing_cache%scratch_4d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%root_fraction_sl_preload(:,:,:,ipft,itimestep) = MERGE( &
          & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
          & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

        CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'veg_demand_uptake_n_tuptake_mavg_'//tile_name, soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%demand_uptake_n_preload(:,:,ipft,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

        CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
          & 'veg_demand_uptake_p_tuptake_mavg_'//tile_name, soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%demand_uptake_p_preload(:,:,ipft,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)
      END DO

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_frac_wtr_transp_down_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%frac_wtr_transp_down_sl_preload(:,:,:,itimestep) = soil_forcing_cache%scratch_4d(:,:,:,1)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_frac_w_lat_loss_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%frac_w_lat_loss_sl_preload(:,:,:,itimestep) = MERGE( &
        & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
        & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_wtr_soil_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%wtr_soil_sl_preload(:,:,:,itimestep) = MERGE( &
        & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
        & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'hydro_vol_porosity_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%vol_porosity_sl_preload(:,:,:,itimestep) = MERGE( &
        & soil_forcing_cache%scratch_4d(:,:,:,1), zero, &
        & soil_forcing_cache%scratch_4d(:,:,:,1) >= zero)

      CALL read_3D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'sse_t_soil_sl_box', soil_forcing_cache%scratch_4d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%t_soil_sl_preload(:,:,:,itimestep) = soil_forcing_cache%scratch_4d(:,:,:,1)

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_nhx_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%nhx_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_noy_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%noy_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_nhx_n15_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%nhx_n15_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_noy_n15_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%noy_n15_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

      CALL read_2D_extdim(soil_forcing_cache%soilphysics_input_file%stream_id, on_cells, &
        & 'a2l_p_deposition_box', soil_forcing_cache%scratch_3d, &
        & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
      soil_forcing_cache%p_deposition_preload(:,:,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)
    END DO

  END SUBROUTINE load_soilphysics_cache

  ! ======================================================================================================= !
  !> Read and cache the forcing fields from bc_qsoil_litterfall.
  !>
  SUBROUTINE load_litterfall_cache(bgcm_store, nsteps, start_timestep, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(in) :: bgcm_store
    INTEGER,                       INTENT(in) :: nsteps
    INTEGER,                       INTENT(in) :: start_timestep
    CHARACTER(len=*),              INTENT(in) :: caller

    INTEGER :: ipft, i_part, i_elem, itimestep, source_timestep
    INTEGER :: veg_litter_idx, seed_bed_idx
    INTEGER :: veg_litter_store_idx, seed_bed_store_idx
    CHARACTER(len=5) :: tile_name

    tile_name = get_pft_tile_name(1)
    veg_litter_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile_name, caller)
    seed_bed_idx   = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_LITTERFALL_ID, tile_name, caller)
    veg_litter_store_idx = bgcm_store%idx_in_store(veg_litter_idx)
    seed_bed_store_idx   = bgcm_store%idx_in_store(seed_bed_idx)

    DO itimestep = 1, nsteps
      source_timestep = start_timestep + itimestep - 1
      DO ipft = 1, soil_forcing_cache%npft
        tile_name = get_pft_tile_name(ipft)

        CALL read_2D_extdim(soil_forcing_cache%veglitterfall_input_file%stream_id, on_cells, &
          & 'veg_veg_bgcm_pool_fine_root_C_'//tile_name, soil_forcing_cache%scratch_3d, &
          & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
        soil_forcing_cache%veg_pool_fine_root_c_preload(:,:,ipft,itimestep) = soil_forcing_cache%scratch_3d(:,:,1)

        DO i_part = 1, soil_forcing_cache%veg_litter_nr_parts
          DO i_elem = 1, soil_forcing_cache%veg_litter_nr_elems
            CALL read_2D_extdim(soil_forcing_cache%veglitterfall_input_file%stream_id, on_cells, &
              & 'veg_veg_bgcm_litterfall_'//get_veg_litter_part_name(i_part)//'_'// &
              & get_veg_litter_elem_name(i_elem)//'_'//tile_name, soil_forcing_cache%scratch_3d, &
              & start_extdim=source_timestep, end_extdim=source_timestep, extdim_name='time')
            soil_forcing_cache%veg_litterfall_preload(i_part, i_elem, :, :, ipft, itimestep) = &
              & soil_forcing_cache%scratch_3d(:,:,1)
          END DO
        END DO

        DO i_elem = 1, soil_forcing_cache%seed_bed_nr_elems
          CALL read_2D_extdim(soil_forcing_cache%veglitterfall_input_file%stream_id, on_cells, &
            & 'veg_veg_bgcm_seed_bed_litterfall_'//get_veg_litter_elem_name(i_elem)//'_'//tile_name, &
            & soil_forcing_cache%scratch_3d, start_extdim=source_timestep, end_extdim=source_timestep, &
            & extdim_name='time')
          soil_forcing_cache%seed_bed_litterfall_preload(i_elem, :, :, ipft, itimestep) = &
            & soil_forcing_cache%scratch_3d(:,:,1)
        END DO
      END DO
    END DO

  END SUBROUTINE load_litterfall_cache

  ! ======================================================================================================= !
  !> Copy cached litterfall and fine-root C into the tile bgcm storage.
  !>
  SUBROUTINE apply_litterfall_cache_to_tile(bgcm_store, tile_name, pft_idx, itimestep, caller)

    TYPE(t_lnd_bgcm_store), POINTER, INTENT(inout) :: bgcm_store
    CHARACTER(len=*),              INTENT(in)    :: tile_name, caller
    INTEGER,                       INTENT(in)    :: pft_idx, itimestep

    INTEGER :: veg_pool_idx, veg_litter_idx, seed_bed_idx
    INTEGER :: veg_pool_store_idx, veg_litter_store_idx, seed_bed_store_idx
    INTEGER :: i_part, i_elem

    veg_pool_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_POOL_ID, tile_name, caller)
    veg_litter_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_LITTERFALL_ID, tile_name, caller)
    seed_bed_idx = get_bgcm_idx(bgcm_store, VEG_BGCM_SEED_BED_LITTERFALL_ID, tile_name, caller)

    veg_pool_store_idx = bgcm_store%idx_in_store(veg_pool_idx)
    veg_litter_store_idx = bgcm_store%idx_in_store(veg_litter_idx)
    seed_bed_store_idx = bgcm_store%idx_in_store(seed_bed_idx)

    bgcm_store%store_2l_2d_bgcms(veg_pool_store_idx)%collection_2l_p(ix_fine_root, ixC)%p%ptr2d(:,:) = &
      & soil_forcing_cache%veg_pool_fine_root_c_preload(:,:,pft_idx,itimestep)

    DO i_part = 1, soil_forcing_cache%veg_litter_nr_parts
      DO i_elem = 1, soil_forcing_cache%veg_litter_nr_elems
        bgcm_store%store_2l_2d_bgcms(veg_litter_store_idx)%collection_2l_p(i_part, i_elem)%p%ptr2d(:,:) = &
          & soil_forcing_cache%veg_litterfall_preload(i_part, i_elem, :, :, pft_idx, itimestep)
      END DO
    END DO

    DO i_elem = 1, soil_forcing_cache%seed_bed_nr_elems
      bgcm_store%store_1l_2d_bgcms(seed_bed_store_idx)%collection_1l_p(i_elem)%p%ptr2d(:,:) = &
        & soil_forcing_cache%seed_bed_litterfall_preload(i_elem, :, :, pft_idx, itimestep)
    END DO

  END SUBROUTINE apply_litterfall_cache_to_tile

  ! ======================================================================================================= !
  !> Resolve soil forcing filename and record index for the current model time.
  !> The filename argument is interpreted as a full filename for yearly forcing
  !> and as a prefix for daily forcing, yielding <prefix>_<doy>.nc.
  !>
  SUBROUTINE get_soil_forcing_info(base_filename, yearday, hour, minute, second, dtime, &
      & resolved_filename, itimestep)

      CHARACTER(len=*),           INTENT(in)  :: base_filename
      INTEGER,                    INTENT(in)  :: yearday, hour, minute, second
      REAL(wp),                   INTENT(in)  :: dtime
      CHARACTER(len=:), ALLOCATABLE, INTENT(out) :: resolved_filename
      INTEGER,                    INTENT(out) :: itimestep

    INTEGER :: seconds_of_day
    CHARACTER(len=8) :: doy_string
    LOGICAL :: is_daily_timestep

    seconds_of_day = hour * 3600 + minute * 60 + second
    is_daily_timestep = NINT(dtime) == 86400

    IF (is_daily_timestep) THEN
      resolved_filename = TRIM(base_filename)
      itimestep = yearday
    ELSE
      WRITE(doy_string,'(I3.3)') yearday
      resolved_filename = trim_nc_filename(base_filename)//'_'//TRIM(doy_string)//'.nc'
      itimestep = seconds_of_day / INT(dtime) + 1
    END IF

  END SUBROUTINE get_soil_forcing_info

  ! ======================================================================================================= !
  !> Return the filename prefix without a trailing .nc suffix.
  !>
  FUNCTION trim_nc_filename(filename) RESULT(prefix)

    CHARACTER(len=*), INTENT(in) :: filename
    CHARACTER(len=:), ALLOCATABLE :: prefix
    INTEGER :: ilen

    ilen = LEN_TRIM(filename)
    IF (ilen >= 3 .AND. filename(ilen-2:ilen) == '.nc') THEN
      prefix = filename(:ilen-3)
    ELSE
      prefix = TRIM(filename)
    END IF

  END FUNCTION trim_nc_filename

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

    ix = (kd-1)*86400 + kh*3600 + km*60 + ks

  END FUNCTION IMerge_TIME2Sec
#endif

#endif
END MODULE mo_sb_util
