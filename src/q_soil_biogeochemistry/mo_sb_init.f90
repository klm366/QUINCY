!> QUINCY soil-biogeochemistry variables init
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
!>#### initialization of soil-biogeochemistry memory variables using, e.g., ic & bc input files
!>
MODULE mo_sb_init
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_exception,             ONLY: finish, message
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_jsb_math_constants,    ONLY: eps4, eps8, one_day
  USE mo_jsb_impl_constants,    ONLY: def_parameters

  USE mo_jsb_class,             ONLY: Get_model
  USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
  USE mo_jsb_model_class,       ONLY: t_jsb_model
  USE mo_jsb_grid_class,        ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,              ONLY: Get_grid, Get_vgrid

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store, get_bgcm_idx
  USE mo_lnd_bgcm_store_class,    ONLY: SB_BGCM_POOL_ID

  USE mo_jsb_process_class,       ONLY: SB_, HYDRO_, SSE_, VEG_, A2L_
  USE mo_sb_config_class,         ONLY: CONST_DEP, REFYEAR_DEP, TRANS_DEP

#ifndef __QUINCY_STANDALONE__
  USE mo_jsb_time,              ONLY: get_year, get_month
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sb_init
#ifndef __QUINCY_STANDALONE__
  PUBLIC :: provide_n_and_p_deposition, sb_read_states
#endif

  TYPE t_sb_init_vars
    REAL(wp), POINTER :: &
      & usda_taxonomy_class     (:,:) => NULL(), &
      & nwrb_taxonomy_class     (:,:) => NULL(), &
      & qmax_org_fine_particle  (:,:) => NULL(), &
      & soil_ph                 (:,:) => NULL(), &
      & soil_p_labile           (:,:) => NULL(), &
      & soil_p_slow             (:,:) => NULL(), &
      & soil_p_occluded         (:,:) => NULL(), &
      & soil_p_primary          (:,:) => NULL()
    INTEGER, POINTER :: &
      & usda_taxonomy_class_int (:,:) => NULL(), &
      & nwrb_taxonomy_class_int (:,:) => NULL()
  END TYPE t_sb_init_vars

  TYPE(t_sb_init_vars) :: sb_init_vars

  TYPE t_sb_analytical_init_vars
    REAL(wp), POINTER :: &
      & frac_wtr_transp_down_sl      (:,:,:) => NULL(), & !< HYDRO downward transport rate used in the analytical closure [s-1]
      & frac_w_lat_loss_sl           (:,:,:) => NULL(), & !< HYDRO lateral loss rate reserved for later analytical closures [s-1]
      & rtm_decomposition            (:,:,:) => NULL(), & !< temperature modifier for simple-model decomposition read from analytical status file [--]
      & rmm_decomposition            (:,:,:) => NULL(), & !< moisture modifier for simple-model decomposition read from analytical status file [--]
      & rtm_depolymerisation         (:,:,:) => NULL(), & !< temperature modifier for depolymerisation read from analytical status file [--]
      & rmm_depolymerisation         (:,:,:) => NULL(), & !< moisture modifier for depolymerisation read from analytical status file [--]
      & rtm_sorption                 (:,:,:) => NULL(), & !< temperature modifier for sorption read from analytical status file [--]
      & rmm_sorption                 (:,:,:) => NULL(), & !< moisture modifier for sorption read from analytical status file [--]
      & rtm_desorption               (:,:,:) => NULL(), & !< temperature modifier for desorption read from analytical status file [--]
      & rmm_desorption               (:,:,:) => NULL(), & !< moisture modifier for desorption read from analytical status file [--]
      & rtm_mic_uptake               (:,:,:) => NULL(), & !< temperature modifier for microbial uptake read from analytical status file [--]
      & rmm_mic_uptake               (:,:,:) => NULL(), & !< moisture modifier for microbial uptake read from analytical status file [--]
      & rtm_hsc                      (:,:,:) => NULL(), & !< temperature modifier on half-saturation terms read from analytical status file [--]
      & rmm_hsc                      (:,:,:) => NULL(), & !< moisture modifier on half-saturation terms read from analytical status file [--]
      & bulk_dens_sl                 (:,:,:) => NULL(), & !< SSE bulk density from analytical status file [kg m-3]
      & silt_sl                      (:,:,:) => NULL(), & !< SSE silt fraction from analytical status file [--]
      & clay_sl                      (:,:,:) => NULL(), & !< SSE clay fraction from analytical status file [--]
      & formation_soluable_litter_c (:,:,:) => NULL(), & !< external carbon input to soluable litter [mol C m-3 s-1]
      & formation_polymeric_litter_c(:,:,:) => NULL(), & !< external carbon input to polymeric litter [mol C m-3 s-1]
      & formation_woody_litter_c    (:,:,:) => NULL(), & !< external carbon input to woody litter [mol C m-3 s-1]
      & enzyme_frac_ap              (:,:,:) => NULL()    !< prescribed phosphatase fraction for the C-only analytical draft [--]
    LOGICAL :: is_loaded = .FALSE.                        !< analytical steady-state init data has been read from file and is ready to be consumed
  END TYPE t_sb_analytical_init_vars

  TYPE(t_sb_analytical_init_vars) :: sb_analytical_init_vars

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_init'

CONTAINS

  ! ======================================================================================================= !tile
  !> Intialize SB_ process
  !>
  SUBROUTINE sb_init(tile)

    USE mo_jsb_parallel,       ONLY: Get_omp_thread
    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    INTEGER :: no_omp_thread
    TYPE(t_jsb_model),  POINTER :: model
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_init'

    dsl4jsb_Def_config(SB_)

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)

#ifdef __QUINCY_STANDALONE__
    CALL sb_qs_read_init_vars(tile)
    IF (TRIM(dsl4jsb_Config(SB_)%sb_init_scheme) == 'analytical_ss') THEN
      CALL sb_qs_read_analytical_init_vars(tile)
    END IF
    CALL sb_init_ic_bc(tile)
    IF (TRIM(dsl4jsb_Config(SB_)%sb_init_scheme) == 'analytical_ss') THEN
      CALL sb_finalize_analytical_init_vars()
    END IF
    CALL sb_finalize_init_vars()
#else
    no_omp_thread = Get_omp_thread()

    ! call only at root tile, i.e., the only tile without associated parent tile (to avoid unnessary i/o)
    IF (.NOT. ASSOCIATED(tile%parent_tile)) THEN
      CALL sb_read_init_vars(tile)
      IF (TRIM(dsl4jsb_Config(SB_)%sb_init_scheme) == 'analytical_ss') THEN
        CALL sb_read_analytical_init_vars(tile)
      END IF
    END IF

    CALL sb_init_ic_bc(tile)

    IF (tile%Is_last_process_tile(SB_)) THEN
      IF (TRIM(dsl4jsb_Config(SB_)%sb_init_scheme) == 'analytical_ss') THEN
        CALL sb_finalize_analytical_init_vars()
      END IF
      CALL sb_finalize_init_vars()
    END IF
#endif

  END SUBROUTINE sb_init

  ! ======================================================================================================= !
  !> Intialize SB_ process from ic and bc input files
  !>
  SUBROUTINE sb_init_ic_bc(tile)
    USE mo_jsb_lctlib_class,              ONLY: t_lctlib_element
    USE mo_q_assimi_constants,            ONLY: IC3PHOT
    USE mo_q_assimi_process,              ONLY: discrimination_ps_c13, discrimination_ps_c14
    USE mo_q_assimi_parameters,           ONLY: CiCa_default_C3, CiCa_default_C4
    USE mo_atmland_constants,             ONLY: def_co2_mixing_ratio, def_co2_mixing_ratio_C13, def_co2_mixing_ratio_C14, &
      &                                         def_co2_deltaC13, def_co2_deltaC14
    USE mo_sb_constants,                  ONLY: microbial_cue_max, microbial_nue, microbial_pue, microbial_cn, microbial_np, &
      &   k_weath_mineral, &
      &   temp_freeze_thres_bio, temp_freeze_window_bio, frac_woody2dec_litter, frac_litter2fast_som, f_nit_noy, f_nit_n2o, &
      &   vmax_nitrification, vmax_denitrification, t_ref_decomposition, ea_decomposition, ed_decomposition
    USE mo_isotope_util,                  ONLY: calc_fractionation, calc_mixing_ratio_C14C
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_memory(HYDRO_)
    dsl4jsb_Use_memory(SSE_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model                 !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store            !< the bgcm store of this tile
    TYPE(t_lctlib_element), POINTER :: lctlib                !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_grid),       POINTER :: hgrid                 !< Horizontal grid
    TYPE(t_jsb_vgrid),      POINTER :: vgrid_soil_w          !< Vertical grid
    INTEGER                         :: nsoil_w               !< number of soil layers as used/defined by the SB_ process
    REAL(wp)                        :: hlp1                  !< helper var
    REAL(wp)                        :: hlp_frac_c13          !< helper var
    REAL(wp)                        :: hlp_frac_c14          !< helper var
    INTEGER                         :: iblk, ic, is          !< loop over dimensions
    INTEGER                         :: nblks                 !< number of blocks
    INTEGER                         :: nproma                !< number of grid points per iblk
    INTEGER                         :: sb_pool_idx           !< index of the sb pool within the bgcm store
    REAL(wp),           POINTER :: sb_pool_mt_domain(:,:,:,:,:) !< dim: compartments, elements, nc, nsoil, nblks
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init_ic_bc'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! VEG_ 3D
    dsl4jsb_Real3D_onDomain :: root_fraction_sl
    ! HYDRO_ 2D
    dsl4jsb_Real2D_onDomain :: num_sl_above_bedrock
    ! HYDRO_ 3D
    dsl4jsb_Real3D_onDomain :: soil_depth_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_width_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_center_sl
    dsl4jsb_Real3D_onDomain :: soil_lay_depth_ubound_sl
    dsl4jsb_Real3D_onDomain :: vol_porosity_sl
    ! SSE_ 3D
    dsl4jsb_Real3D_onDomain :: bulk_dens_sl
    dsl4jsb_Real3D_onDomain :: silt_sl
    dsl4jsb_Real3D_onDomain :: clay_sl
    ! SB_ 3D
    dsl4jsb_Real3D_onDomain :: nh4_assoc
    dsl4jsb_Real3D_onDomain :: nh4_solute
    dsl4jsb_Real3D_onDomain :: no3_solute
    dsl4jsb_Real3D_onDomain :: po4_solute
    dsl4jsb_Real3D_onDomain :: nh4_n15_assoc
    dsl4jsb_Real3D_onDomain :: nh4_n15_solute
    dsl4jsb_Real3D_onDomain :: no3_n15_solute
    dsl4jsb_Real3D_onDomain :: bulk_soil_carbon_sl
    dsl4jsb_Real3D_onDomain :: soil_litter_carbon_sl
    dsl4jsb_Real3D_onDomain :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onDomain :: qmax_org
    dsl4jsb_Real3D_onDomain :: surface_burial_frac
    dsl4jsb_Real3D_onDomain :: bulk_carbon_max
    dsl4jsb_Real3D_onDomain :: qmax_po4
    dsl4jsb_Real3D_onDomain :: qmax_nh4
    dsl4jsb_Real3D_onDomain :: qmax_fast_po4
    dsl4jsb_Real3D_onDomain :: qmax_slow_po4
    dsl4jsb_Real3D_onDomain :: km_fast_po4
    dsl4jsb_Real3D_onDomain :: km_slow_po4
    dsl4jsb_Real3D_onDomain :: km_adsorpt_po4_sl
    dsl4jsb_Real3D_onDomain :: km_adsorpt_nh4_sl
    dsl4jsb_Real3D_onDomain :: k_bioturb
    dsl4jsb_Real3D_onDomain :: po4_assoc_fast
    dsl4jsb_Real3D_onDomain :: po4_assoc_slow
    dsl4jsb_Real3D_onDomain :: po4_occluded
    dsl4jsb_Real3D_onDomain :: ph_sl
    dsl4jsb_Real3D_onDomain :: Qmax_AlFe_cor
    dsl4jsb_Real3D_onDomain :: po4_primary
    dsl4jsb_Real3D_onDomain :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: vmax_weath_mineral_sl
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly
    dsl4jsb_Real3D_onDomain :: enzyme_frac_residue
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_c
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_n
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_p
    dsl4jsb_Real3D_onDomain :: dom_cn
    dsl4jsb_Real3D_onDomain :: dom_cp
    dsl4jsb_Real3D_onDomain :: fact_n_status_mic_c_growth
    dsl4jsb_Real3D_onDomain :: fact_p_status_mic_c_growth
    dsl4jsb_Real3D_onDomain :: enzyme_frac_AP
    dsl4jsb_Real3D_onDomain :: qmax_org_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_po4_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_nh4_min_sl
    dsl4jsb_Real3D_onDomain :: qmax_po4_om_sl
    dsl4jsb_Real3D_onDomain :: volume_min_sl
    dsl4jsb_Real3D_onDomain :: rtm_depolymerisation
    dsl4jsb_Real3D_onDomain :: rtm_mic_uptake
    dsl4jsb_Real3D_onDomain :: rtm_sorption
    dsl4jsb_Real3D_onDomain :: rtm_desorption
    dsl4jsb_Real3D_onDomain :: rtm_hsc
    dsl4jsb_Real3D_onDomain :: rmm_depolymerisation
    dsl4jsb_Real3D_onDomain :: rmm_mic_uptake
    dsl4jsb_Real3D_onDomain :: rmm_sorption
    dsl4jsb_Real3D_onDomain :: rmm_desorption
    dsl4jsb_Real3D_onDomain :: rmm_hsc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN                !< run this init only if the present tile is a pft
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    hgrid         => Get_grid(model%grid_id)
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels
    nblks         =  hgrid%nblks
    nproma        =  hgrid%nproma
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    sb_pool_idx = get_bgcm_idx(bgcm_store, SB_BGCM_POOL_ID, tile%name, routine)
    sb_pool_mt_domain => bgcm_store%store_2l_3d_bgcms(bgcm_store%idx_in_store(sb_pool_idx))%mt_2l_3d_bgcm(:,:,:,:,:)
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_var3D_onDomain(VEG_,   root_fraction_sl)
    ! HYDRO_ 2D
    dsl4jsb_Get_var2D_onDomain(HYDRO_, num_sl_above_bedrock)     ! in
    ! HYDRO_ 3D
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_depth_sl)            ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_lay_width_sl)        ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_lay_depth_center_sl) ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_, soil_lay_depth_ubound_sl) ! in
    dsl4jsb_Get_var3D_onDomain(HYDRO_, vol_porosity_sl)
    ! SSE_ 3D
    dsl4jsb_Get_var3D_onDomain(SSE_,   bulk_dens_sl)             ! in
    dsl4jsb_Get_var3D_onDomain(SSE_,   silt_sl)
    dsl4jsb_Get_var3D_onDomain(SSE_,   clay_sl)
    ! ---------------------------
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_assoc)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    no3_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_n15_assoc)
    dsl4jsb_Get_var3D_onDomain(SB_,    nh4_n15_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    no3_n15_solute)
    dsl4jsb_Get_var3D_onDomain(SB_,    bulk_soil_carbon_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    soil_litter_carbon_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    bulk_dens_corr_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_org)
    dsl4jsb_Get_var3D_onDomain(SB_,    surface_burial_frac)
    dsl4jsb_Get_var3D_onDomain(SB_,    bulk_carbon_max)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_nh4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_fast_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_slow_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_fast_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_slow_po4)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_adsorpt_po4_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    km_adsorpt_nh4_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    k_bioturb)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_assoc_fast)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_assoc_slow)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_occluded)
    dsl4jsb_Get_var3D_onDomain(SB_,    ph_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    Qmax_AlFe_cor)
    dsl4jsb_Get_var3D_onDomain(SB_,    po4_primary)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_,    vmax_weath_mineral_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_residue)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_c)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_n)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_poly_p)
    dsl4jsb_Get_var3D_onDomain(SB_,    dom_cn)
    dsl4jsb_Get_var3D_onDomain(SB_,    dom_cp)
    dsl4jsb_Get_var3D_onDomain(SB_,    fact_n_status_mic_c_growth)
    dsl4jsb_Get_var3D_onDomain(SB_,    fact_p_status_mic_c_growth)
    dsl4jsb_Get_var3D_onDomain(SB_,    enzyme_frac_AP)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_org_min_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_po4_min_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_nh4_min_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    qmax_po4_om_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    volume_min_sl)
    dsl4jsb_Get_var3D_onDomain(SB_,    rtm_depolymerisation)
    dsl4jsb_Get_var3D_onDomain(SB_,    rtm_mic_uptake)
    dsl4jsb_Get_var3D_onDomain(SB_,    rtm_sorption)
    dsl4jsb_Get_var3D_onDomain(SB_,    rtm_desorption)
    dsl4jsb_Get_var3D_onDomain(SB_,    rtm_hsc)
    dsl4jsb_Get_var3D_onDomain(SB_,    rmm_depolymerisation)
    dsl4jsb_Get_var3D_onDomain(SB_,    rmm_mic_uptake)
    dsl4jsb_Get_var3D_onDomain(SB_,    rmm_sorption)
    dsl4jsb_Get_var3D_onDomain(SB_,    rmm_desorption)
    dsl4jsb_Get_var3D_onDomain(SB_,    rmm_hsc)
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.8 apply calibration of soil-biogeochemistry constants if enabled
    !>
    IF (dsl4jsb_Config(SB_)%flag_apply_nml_parameters) THEN
      ! temp_freeze_thres_bio
      IF (dsl4jsb_Config(SB_)%nml_temp_freeze_thres_bio > def_parameters) THEN
        temp_freeze_thres_bio = dsl4jsb_Config(SB_)%nml_temp_freeze_thres_bio
        !$ACC UPDATE DEVICE(temp_freeze_thres_bio) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: temp_freeze_thres_bio')
      END IF
      ! temp_freeze_window_bio
      IF (dsl4jsb_Config(SB_)%nml_temp_freeze_window_bio > def_parameters) THEN
        temp_freeze_window_bio = dsl4jsb_Config(SB_)%nml_temp_freeze_window_bio
        !$ACC UPDATE DEVICE(temp_freeze_window_bio) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: temp_freeze_window_bio')
      END IF
      ! frac_woody2dec_litter
      IF (dsl4jsb_Config(SB_)%nml_frac_woody2dec_litter > def_parameters) THEN
        frac_woody2dec_litter = dsl4jsb_Config(SB_)%nml_frac_woody2dec_litter
        !$ACC UPDATE DEVICE(frac_woody2dec_litter) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: frac_woody2dec_litter')
      END IF
      ! frac_litter2fast_som
      IF (dsl4jsb_Config(SB_)%nml_frac_litter2fast_som > def_parameters) THEN
        frac_litter2fast_som = dsl4jsb_Config(SB_)%nml_frac_litter2fast_som
        !$ACC UPDATE DEVICE(frac_litter2fast_som) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: frac_litter2fast_som')
      END IF
      ! f_nit_noy
      IF (dsl4jsb_Config(SB_)%nml_f_nit_noy > def_parameters) THEN
        f_nit_noy = dsl4jsb_Config(SB_)%nml_f_nit_noy
        !$ACC UPDATE DEVICE(f_nit_noy) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: f_nit_noy')
      END IF
      ! f_nit_n2o
      IF (dsl4jsb_Config(SB_)%nml_f_nit_n2o > def_parameters) THEN
        f_nit_n2o = dsl4jsb_Config(SB_)%nml_f_nit_n2o
        !$ACC UPDATE DEVICE(f_nit_n2o) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: f_nit_n2o')
      END IF
      ! vmax_nitrification
      IF (dsl4jsb_Config(SB_)%nml_vmax_nitrification > def_parameters) THEN
        vmax_nitrification = dsl4jsb_Config(SB_)%nml_vmax_nitrification / one_day * 1.E6_wp
        !$ACC UPDATE DEVICE(vmax_nitrification) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: vmax_nitrification')
      END IF
      ! vmax_denitrification
      IF (dsl4jsb_Config(SB_)%nml_vmax_denitrification > def_parameters) THEN
        vmax_denitrification = dsl4jsb_Config(SB_)%nml_vmax_denitrification / one_day * 1.E6_wp
        !$ACC UPDATE DEVICE(vmax_denitrification) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: vmax_denitrification')
      END IF
      ! t_ref_decomposition
      IF (dsl4jsb_Config(SB_)%nml_t_ref_decomposition > def_parameters) THEN
        t_ref_decomposition = dsl4jsb_Config(SB_)%nml_t_ref_decomposition
        !$ACC UPDATE DEVICE(t_ref_decomposition) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: t_ref_decomposition')
      END IF
      ! ea_decomposition
      IF (dsl4jsb_Config(SB_)%nml_ea_decomposition > def_parameters) THEN
        ea_decomposition = dsl4jsb_Config(SB_)%nml_ea_decomposition
        !$ACC UPDATE DEVICE(ea_decomposition) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: ea_decomposition')
      END IF
      ! ed_decomposition
      IF (dsl4jsb_Config(SB_)%nml_ed_decomposition > def_parameters) THEN
        ed_decomposition = dsl4jsb_Config(SB_)%nml_ed_decomposition
        !$ACC UPDATE DEVICE(ed_decomposition) ASYNC(1)
        CALL message(TRIM(routine), ' Modified value of vegetation parameter: ed_decomposition')
      END IF
    END IF

    !>0.9 init bgcm with zero
    !>  NOTE: reading the restart files is done after init
    !>  NOTE: not using a "Get" function in process init for bgcm
    !>
    sb_pool_mt_domain(:,:,:,:,:) = 0.0_wp

    !>1.0 all soil models
    !>  simple_sm & jsm
    !>

    !>  1.1 calc hlp_frac_c13 & hlp_frac_c14
    !>
    ! C13
    IF (lctlib%ps_pathway == IC3PHOT) THEN
      hlp1 = discrimination_ps_c13(CiCa_default_C3 * def_co2_mixing_ratio, &
        &                          def_co2_mixing_ratio, &
        &                          lctlib%ps_pathway)
    ELSE
      hlp1 = discrimination_ps_c13(CiCa_default_C4 * def_co2_mixing_ratio, &
        &                          def_co2_mixing_ratio, &
        &                          lctlib%ps_pathway)
    END IF
    hlp_frac_c13 = calc_fractionation(def_co2_mixing_ratio, &
      &                               def_co2_mixing_ratio_C13, &
      &                               hlp1 )
    ! C14
    IF (lctlib%ps_pathway == IC3PHOT) THEN
      hlp1 = discrimination_ps_c14(CiCa_default_C3 * def_co2_mixing_ratio, &
        &                          def_co2_mixing_ratio, &
        &                          lctlib%ps_pathway)
    ELSE
      hlp1 = discrimination_ps_c14(CiCa_default_C4 * def_co2_mixing_ratio, &
        &                          def_co2_mixing_ratio, &
        &                          lctlib%ps_pathway)
    END IF
    hlp_frac_c14 = calc_mixing_ratio_C14C(def_co2_deltaC13, &
      &                                   def_co2_deltaC14 - hlp1)

    !>  1.2 init site-specific soil properties
    !>
    ph_sl(:,:,:)              = SPREAD(sb_init_vars%soil_ph(:,:), DIM = 2, ncopies = nsoil_w)
    Qmax_AlFe_cor(:,:,:)      = 1._wp                           !< ... docu ...

    !>  1.3 init SB_ variables from constants
    !>
    microbial_cue_eff_tmic_mavg(:,:,:)  = microbial_cue_max
    microbial_nue_eff_tmic_mavg(:,:,:)  = microbial_nue
    microbial_pue_eff_tmic_mavg(:,:,:)  = microbial_pue
    vmax_weath_mineral_sl(:,:,:)        = k_weath_mineral
    dom_cn(:,:,:)                       = microbial_cn
    dom_cp(:,:,:)                       = microbial_cn * microbial_np
    fact_n_status_mic_c_growth(:,:,:)   = 1._wp              !< [0,1] microbial N status, to scale up/down some microbial processes
    fact_p_status_mic_c_growth(:,:,:)   = 1._wp              !< [0,1] microbial P status, to scale up/down some microbial processes

    SELECT CASE(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme))
    !>2.0 ssm - a simple four pool soil biogeochemical model
    !>
    CASE("simple_1d")
      CALL sb_init_simple_sm( &
        & nproma, &                                   ! in
        & nsoil_w, &
        & nblks, &
        & num_sl_above_bedrock(:,:), &
        & lctlib%growthform, &
        & lctlib%k_som_fast_init, &
        & lctlib%k_som_slow_init, &
        & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
        & TRIM(dsl4jsb_Config(SB_)%sb_init_scheme), &
        & dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
        & SPREAD(sb_init_vars%qmax_org_fine_particle(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_labile(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_slow(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_occluded(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_primary(:,:), DIM = 2, ncopies = nsoil_w), &
        & dsl4jsb_Config(SB_)%soil_p_depth, &
        & soil_depth_sl(:,:,:), &
        & soil_lay_width_sl(:,:,:), &
        & soil_lay_depth_ubound_sl(:,:,:), &
        & soil_lay_depth_center_sl(:,:,:), &
        & hlp_frac_c13, &
        & hlp_frac_c14, &
        & root_fraction_sl(:,:,:), &
        & clay_sl(:,:,:), &
        & silt_sl(:,:,:), &
        & vol_porosity_sl(:,:,:), &
        & bulk_dens_sl(:,:,:), &                      ! in
        & ph_sl(:,:,:), &                             ! inout
        & bulk_carbon_max(:,:,:), &
        & surface_burial_frac(:,:,:), &
        & bulk_soil_carbon_sl(:,:,:), &
        & soil_litter_carbon_sl(:,:,:), &
        & po4_solute(:,:,:), &
        & volume_min_sl(:,:,:), &
        & qmax_po4(:,:,:), &
        & qmax_po4_min_sl(:,:,:), &
        & qmax_po4_om_sl(:,:,:), &
        & qmax_fast_po4(:,:,:), &
        & qmax_slow_po4(:,:,:), &
        & km_fast_po4(:,:,:), &
        & km_slow_po4(:,:,:), &
        & km_adsorpt_po4_sl(:,:,:), &
        & po4_assoc_fast(:,:,:), &
        & po4_assoc_slow(:,:,:), &
        & po4_occluded(:,:,:), &
        & po4_primary(:,:,:), &
        & Qmax_AlFe_cor(:,:,:), &
        & sb_pool_mt_domain(:,:,:,:,:), &                    ! inout
        & nh4_solute(:,:,:), &                        ! out
        & no3_solute(:,:,:), &
        & nh4_n15_solute(:,:,:), &
        & no3_n15_solute(:,:,:), &
        & bulk_dens_corr_sl(:,:,:), &
        & qmax_org(:,:,:), &
        & qmax_org_min_sl(:,:,:), &
        & qmax_nh4(:,:,:), &
        & qmax_nh4_min_sl(:,:,:), &
        & km_adsorpt_nh4_sl(:,:,:), &
        & nh4_assoc(:,:,:), &
        & nh4_n15_assoc(:,:,:) &                      ! out
        & )
    !>3.0 jsm - Jena Soil Model (described in https://doi.org/10.5194/gmd-13-783-2020)
    !>
    CASE("jsm")
      CALL sb_init_jsm( &
        & nproma, &                                   ! in
        & nsoil_w, &
        & nblks, &
        & num_sl_above_bedrock(:,:), &
        & lctlib%growthform, &
        & TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
        & TRIM(dsl4jsb_Config(SB_)%sb_init_scheme), &
        & dsl4jsb_Config(SB_)%flag_sb_double_langmuir, &
        & SPREAD(sb_init_vars%qmax_org_fine_particle(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_labile(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_slow(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_occluded(:,:), DIM = 2, ncopies = nsoil_w), &
        & SPREAD(sb_init_vars%soil_p_primary(:,:), DIM = 2, ncopies = nsoil_w), &
        & dsl4jsb_Config(SB_)%soil_p_depth, &
        & soil_depth_sl(:,:,:), &
        & soil_lay_width_sl(:,:,:), &
        & soil_lay_depth_ubound_sl(:,:,:), &
        & soil_lay_depth_center_sl(:,:,:), &
        & hlp_frac_c13, &
        & hlp_frac_c14, &
        & root_fraction_sl(:,:,:), &
        & clay_sl(:,:,:), &
        & silt_sl(:,:,:), &
        & vol_porosity_sl(:,:,:), &
        & bulk_dens_sl(:,:,:), &                      ! in
        & ph_sl(:,:,:), &                             ! inout
        & qmax_po4_min_sl(:,:,:), &
        & qmax_po4_om_sl(:,:,:), &
        & volume_min_sl(:,:,:), &
        & po4_solute(:,:,:), &
        & bulk_carbon_max(:,:,:), &
        & surface_burial_frac(:,:,:), &
        & bulk_soil_carbon_sl(:,:,:), &
        & soil_litter_carbon_sl(:,:,:), &
        & qmax_po4(:,:,:), &
        & qmax_nh4(:,:,:), &
        & qmax_fast_po4(:,:,:), &
        & qmax_slow_po4(:,:,:), &
        & km_fast_po4(:,:,:), &
        & km_slow_po4(:,:,:), &
        & km_adsorpt_po4_sl(:,:,:), &
        & po4_assoc_fast(:,:,:), &
        & po4_assoc_slow(:,:,:), &
        & po4_occluded(:,:,:), &
        & Qmax_AlFe_cor(:,:,:), &
        & po4_primary(:,:,:), &
        & sb_pool_mt_domain(:,:,:,:,:), &                    ! inout
        & qmax_org_min_sl(:,:,:), &                   ! out
        & qmax_nh4_min_sl(:,:,:), &
        & nh4_assoc(:,:,:), &
        & nh4_solute(:,:,:), &
        & no3_solute(:,:,:), &
        & nh4_n15_assoc(:,:,:), &
        & nh4_n15_solute(:,:,:), &
        & no3_n15_solute(:,:,:), &
        & bulk_dens_corr_sl(:,:,:), &
        & qmax_org(:,:,:), &
        & km_adsorpt_nh4_sl(:,:,:), &
        & rtm_depolymerisation(:,:,:), &
        & rmm_depolymerisation(:,:,:), &
        & rtm_mic_uptake(:,:,:), &
        & rmm_mic_uptake(:,:,:), &
        & rtm_sorption(:,:,:), &
        & rmm_sorption(:,:,:), &
        & rtm_desorption(:,:,:), &
        & rmm_desorption(:,:,:), &
        & rtm_hsc(:,:,:), &
        & rmm_hsc(:,:,:), &
        & k_bioturb(:,:,:), &
        & enzyme_frac_poly(:,:,:), &
        & enzyme_frac_residue(:,:,:), &
        & enzyme_frac_poly_c(:,:,:), &
        & enzyme_frac_poly_n(:,:,:), &
        & enzyme_frac_poly_p(:,:,:), &
        & enzyme_frac_AP(:,:,:) &                     ! out
        & )
    CASE DEFAULT
      CALL finish(routine, 'Invalid setting of sb_model_scheme (simple_sm / jsm).')
    END SELECT  ! simple_sm | jsm

    !> 3.0 set sb_pool bgcm to zero for tiles such as bare soil and urban area
    !>
    IF (lctlib%BareSoilFlag) THEN
      sb_pool_mt_domain(:,:,:,:,:) = 0.0_wp
    END IF

    ! IN
    !$ACC UPDATE DEVICE(root_fraction_sl(:,:,:), num_sl_above_bedrock(:,:), bulk_dens_sl(:,:,:), vol_porosity_sl(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(soil_depth_sl(:,:,:), soil_lay_width_sl(:,:,:), soil_lay_depth_center_sl(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(soil_lay_depth_ubound_sl(:,:,:), silt_sl(:,:,:), clay_sl(:,:,:)) ASYNC(1)

    ! IN-OUT
    !$ACC UPDATE DEVICE(qmax_org_min_sl(:,:,:), qmax_po4_min_sl(:,:,:), qmax_nh4_min_sl(:,:,:), qmax_po4_om_sl(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(volume_min_sl(:,:,:)) ASYNC(1)

    ! SB_
    !$ACC UPDATE DEVICE(nh4_assoc(:,:,:), nh4_solute(:,:,:), no3_solute(:,:,:), po4_solute(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(nh4_n15_assoc(:,:,:), nh4_n15_solute(:,:,:), no3_n15_solute(:,:,:), bulk_soil_carbon_sl(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(surface_burial_frac(:,:,:), bulk_carbon_max(:,:,:), soil_litter_carbon_sl(:,:,:), bulk_dens_corr_sl(:,:,:), qmax_org(:,:,:), qmax_po4(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(qmax_nh4(:,:,:), qmax_fast_po4(:,:,:), qmax_slow_po4(:,:,:), km_fast_po4(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(km_slow_po4(:,:,:), km_adsorpt_po4_sl(:,:,:), km_adsorpt_nh4_sl(:,:,:), k_bioturb(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(po4_assoc_fast(:,:,:), po4_assoc_slow(:,:,:), po4_occluded(:,:,:), ph_sl(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(Qmax_AlFe_cor(:,:,:), po4_primary(:,:,:), microbial_cue_eff_tmic_mavg(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(microbial_nue_eff_tmic_mavg(:,:,:), microbial_pue_eff_tmic_mavg(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(vmax_weath_mineral_sl(:,:,:), enzyme_frac_poly(:,:,:), enzyme_frac_residue(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(enzyme_frac_poly_c(:,:,:), enzyme_frac_poly_n(:,:,:), enzyme_frac_poly_p(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(dom_cn(:,:,:), dom_cp(:,:,:), fact_n_status_mic_c_growth(:,:,:)) ASYNC(1)
    !$ACC UPDATE DEVICE(fact_p_status_mic_c_growth(:,:,:), enzyme_frac_AP(:,:,:)) ASYNC(1)

    !$ACC UPDATE DEVICE(sb_pool_mt_domain(:,:,:,:,:)) ASYNC(1)

  END SUBROUTINE sb_init_ic_bc

  ! ======================================================================================================= !
  !> de-allocate SB_ init vars
  !>
  SUBROUTINE sb_finalize_init_vars
    DEALLOCATE( &
      & sb_init_vars%usda_taxonomy_class      , &
      & sb_init_vars%nwrb_taxonomy_class      , &
      & sb_init_vars%qmax_org_fine_particle   , &
      & sb_init_vars%soil_ph                  , &
      & sb_init_vars%soil_p_labile            , &
      & sb_init_vars%soil_p_slow              , &
      & sb_init_vars%soil_p_occluded          , &
      & sb_init_vars%soil_p_primary           , &
      & sb_init_vars%usda_taxonomy_class_int  , &
      & sb_init_vars%nwrb_taxonomy_class_int)
  END SUBROUTINE sb_finalize_init_vars

  ! ======================================================================================================= !
  !> de-allocate analytical steady-state SB_ init vars
  !>
  SUBROUTINE sb_finalize_analytical_init_vars
    IF (ASSOCIATED(sb_analytical_init_vars%frac_wtr_transp_down_sl)) THEN
      DEALLOCATE(sb_analytical_init_vars%frac_wtr_transp_down_sl)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%frac_w_lat_loss_sl)) THEN
      DEALLOCATE(sb_analytical_init_vars%frac_w_lat_loss_sl)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rtm_decomposition)) THEN
      DEALLOCATE(sb_analytical_init_vars%rtm_decomposition)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rmm_decomposition)) THEN
      DEALLOCATE(sb_analytical_init_vars%rmm_decomposition)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rtm_depolymerisation)) THEN
      DEALLOCATE(sb_analytical_init_vars%rtm_depolymerisation)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rmm_depolymerisation)) THEN
      DEALLOCATE(sb_analytical_init_vars%rmm_depolymerisation)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rtm_sorption)) THEN
      DEALLOCATE(sb_analytical_init_vars%rtm_sorption)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rmm_sorption)) THEN
      DEALLOCATE(sb_analytical_init_vars%rmm_sorption)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rtm_desorption)) THEN
      DEALLOCATE(sb_analytical_init_vars%rtm_desorption)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rmm_desorption)) THEN
      DEALLOCATE(sb_analytical_init_vars%rmm_desorption)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rtm_mic_uptake)) THEN
      DEALLOCATE(sb_analytical_init_vars%rtm_mic_uptake)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rmm_mic_uptake)) THEN
      DEALLOCATE(sb_analytical_init_vars%rmm_mic_uptake)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rtm_hsc)) THEN
      DEALLOCATE(sb_analytical_init_vars%rtm_hsc)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%rmm_hsc)) THEN
      DEALLOCATE(sb_analytical_init_vars%rmm_hsc)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%bulk_dens_sl)) THEN
      DEALLOCATE(sb_analytical_init_vars%bulk_dens_sl)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%silt_sl)) THEN
      DEALLOCATE(sb_analytical_init_vars%silt_sl)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%clay_sl)) THEN
      DEALLOCATE(sb_analytical_init_vars%clay_sl)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%formation_soluable_litter_c)) THEN
      DEALLOCATE(sb_analytical_init_vars%formation_soluable_litter_c)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%formation_polymeric_litter_c)) THEN
      DEALLOCATE(sb_analytical_init_vars%formation_polymeric_litter_c)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%formation_woody_litter_c)) THEN
      DEALLOCATE(sb_analytical_init_vars%formation_woody_litter_c)
    END IF
    IF (ASSOCIATED(sb_analytical_init_vars%enzyme_frac_ap)) THEN
      DEALLOCATE(sb_analytical_init_vars%enzyme_frac_ap)
    END IF
    sb_analytical_init_vars%is_loaded = .FALSE.
  END SUBROUTINE sb_finalize_analytical_init_vars

  ! ======================================================================================================= !
  !> allocate analytical steady-state SB_ init vars
  !>
  SUBROUTINE sb_allocate_analytical_init_vars(nproma, nsoil_w, nblks)
    INTEGER, INTENT(in) :: nproma
    INTEGER, INTENT(in) :: nsoil_w
    INTEGER, INTENT(in) :: nblks

    ALLOCATE( &
      & sb_analytical_init_vars%frac_wtr_transp_down_sl      (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%frac_w_lat_loss_sl           (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rtm_decomposition            (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rmm_decomposition            (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rtm_depolymerisation         (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rmm_depolymerisation         (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rtm_sorption                 (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rmm_sorption                 (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rtm_desorption               (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rmm_desorption               (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rtm_mic_uptake               (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rmm_mic_uptake               (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rtm_hsc                      (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%rmm_hsc                      (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%bulk_dens_sl                 (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%silt_sl                      (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%clay_sl                      (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%formation_soluable_litter_c (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%formation_polymeric_litter_c(nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%formation_woody_litter_c    (nproma, nsoil_w, nblks), &
      & sb_analytical_init_vars%enzyme_frac_ap              (nproma, nsoil_w, nblks))
  END SUBROUTINE sb_allocate_analytical_init_vars

  ! ======================================================================================================= !
  !> set defaults for analytical steady-state SB_ init vars
  !>
  SUBROUTINE sb_set_analytical_init_vars_defaults()
    USE mo_sb_constants, ONLY: min_enzyme_fraction

    sb_analytical_init_vars%frac_wtr_transp_down_sl      (:,:,:) = 0.0_wp
    sb_analytical_init_vars%frac_w_lat_loss_sl           (:,:,:) = 0.0_wp
    sb_analytical_init_vars%rtm_decomposition            (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rmm_decomposition            (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rtm_depolymerisation         (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rmm_depolymerisation         (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rtm_sorption                 (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rmm_sorption                 (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rtm_desorption               (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rmm_desorption               (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rtm_mic_uptake               (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rmm_mic_uptake               (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rtm_hsc                      (:,:,:) = 1.0_wp
    sb_analytical_init_vars%rmm_hsc                      (:,:,:) = 1.0_wp
    sb_analytical_init_vars%bulk_dens_sl                 (:,:,:) = 0.0_wp
    sb_analytical_init_vars%silt_sl                      (:,:,:) = 0.0_wp
    sb_analytical_init_vars%clay_sl                      (:,:,:) = 0.0_wp
    sb_analytical_init_vars%formation_soluable_litter_c (:,:,:) = 0.0_wp
    sb_analytical_init_vars%formation_polymeric_litter_c(:,:,:) = 0.0_wp
    sb_analytical_init_vars%formation_woody_litter_c    (:,:,:) = 0.0_wp
    sb_analytical_init_vars%enzyme_frac_ap              (:,:,:) = min_enzyme_fraction
  END SUBROUTINE sb_set_analytical_init_vars_defaults

  ! ======================================================================================================= !
  !> assign one scalar value to one analytical steady-state init variable
  !>
  SUBROUTINE sb_assign_analytical_init_scalar(var_name, scalar_value)
    CHARACTER(len=*), INTENT(in) :: var_name
    REAL(wp),         INTENT(in) :: scalar_value

    SELECT CASE(TRIM(var_name))
    CASE('frac_wtr_transp_down_sl', 'hydro_frac_wtr_transp_down_sl')
      sb_analytical_init_vars%frac_wtr_transp_down_sl(:,:,:) = scalar_value
    CASE('frac_w_lat_loss_sl', 'hydro_frac_w_lat_loss_sl')
      sb_analytical_init_vars%frac_w_lat_loss_sl(:,:,:) = scalar_value
    CASE('rtm_decomposition', 'sb_rtm_decomposition')
      sb_analytical_init_vars%rtm_decomposition(:,:,:) = scalar_value
    CASE('rmm_decomposition', 'sb_rmm_decomposition')
      sb_analytical_init_vars%rmm_decomposition(:,:,:) = scalar_value
    CASE('rtm_depolymerisation', 'sb_rtm_depolymerisation')
      sb_analytical_init_vars%rtm_depolymerisation(:,:,:) = scalar_value
    CASE('rmm_depolymerisation', 'sb_rmm_depolymerisation')
      sb_analytical_init_vars%rmm_depolymerisation(:,:,:) = scalar_value
    CASE('rtm_sorption', 'sb_rtm_sorption')
      sb_analytical_init_vars%rtm_sorption(:,:,:) = scalar_value
    CASE('rmm_sorption', 'sb_rmm_sorption')
      sb_analytical_init_vars%rmm_sorption(:,:,:) = scalar_value
    CASE('rtm_desorption', 'sb_rtm_desorption')
      sb_analytical_init_vars%rtm_desorption(:,:,:) = scalar_value
    CASE('rmm_desorption', 'sb_rmm_desorption')
      sb_analytical_init_vars%rmm_desorption(:,:,:) = scalar_value
    CASE('rtm_mic_uptake', 'sb_rtm_mic_uptake')
      sb_analytical_init_vars%rtm_mic_uptake(:,:,:) = scalar_value
    CASE('rmm_mic_uptake', 'sb_rmm_mic_uptake')
      sb_analytical_init_vars%rmm_mic_uptake(:,:,:) = scalar_value
    CASE('rtm_hsc', 'sb_rtm_hsc')
      sb_analytical_init_vars%rtm_hsc(:,:,:) = scalar_value
    CASE('rmm_hsc', 'sb_rmm_hsc')
      sb_analytical_init_vars%rmm_hsc(:,:,:) = scalar_value
    CASE('bulk_dens_sl', 'sse_bulk_dens_sl')
      sb_analytical_init_vars%bulk_dens_sl(:,:,:) = scalar_value
    CASE('silt_sl', 'sse_silt_sl')
      sb_analytical_init_vars%silt_sl(:,:,:) = scalar_value
    CASE('clay_sl', 'sse_clay_sl')
      sb_analytical_init_vars%clay_sl(:,:,:) = scalar_value
    CASE('formation_soluable_litter_c')
      sb_analytical_init_vars%formation_soluable_litter_c(:,:,:) = scalar_value
    CASE('formation_polymeric_litter_c')
      sb_analytical_init_vars%formation_polymeric_litter_c(:,:,:) = scalar_value
    CASE('formation_woody_litter_c')
      sb_analytical_init_vars%formation_woody_litter_c(:,:,:) = scalar_value
    CASE DEFAULT
      CALL finish(TRIM(modname)//':sb_assign_analytical_init_scalar', &
        & 'Unknown analytical steady-state init variable in standalone text file: '//TRIM(var_name))
    END SELECT
  END SUBROUTINE sb_assign_analytical_init_scalar

  ! ======================================================================================================= !
  !> assign one scalar value to one analytical steady-state init variable and one specific soil layer
  !>
  SUBROUTINE sb_assign_analytical_init_layer_scalar(var_name, scalar_value, isoil_layer)
    CHARACTER(len=*), INTENT(in) :: var_name
    REAL(wp),         INTENT(in) :: scalar_value
    INTEGER,          INTENT(in) :: isoil_layer

    SELECT CASE(TRIM(var_name))
    CASE('frac_wtr_transp_down_sl', 'hydro_frac_wtr_transp_down_sl')
      sb_analytical_init_vars%frac_wtr_transp_down_sl(:,isoil_layer,:) = scalar_value
    CASE('frac_w_lat_loss_sl', 'hydro_frac_w_lat_loss_sl')
      sb_analytical_init_vars%frac_w_lat_loss_sl(:,isoil_layer,:) = scalar_value
    CASE('rtm_decomposition', 'sb_rtm_decomposition')
      sb_analytical_init_vars%rtm_decomposition(:,isoil_layer,:) = scalar_value
    CASE('rmm_decomposition', 'sb_rmm_decomposition')
      sb_analytical_init_vars%rmm_decomposition(:,isoil_layer,:) = scalar_value
    CASE('rtm_depolymerisation', 'sb_rtm_depolymerisation')
      sb_analytical_init_vars%rtm_depolymerisation(:,isoil_layer,:) = scalar_value
    CASE('rmm_depolymerisation', 'sb_rmm_depolymerisation')
      sb_analytical_init_vars%rmm_depolymerisation(:,isoil_layer,:) = scalar_value
    CASE('rtm_sorption', 'sb_rtm_sorption')
      sb_analytical_init_vars%rtm_sorption(:,isoil_layer,:) = scalar_value
    CASE('rmm_sorption', 'sb_rmm_sorption')
      sb_analytical_init_vars%rmm_sorption(:,isoil_layer,:) = scalar_value
    CASE('rtm_desorption', 'sb_rtm_desorption')
      sb_analytical_init_vars%rtm_desorption(:,isoil_layer,:) = scalar_value
    CASE('rmm_desorption', 'sb_rmm_desorption')
      sb_analytical_init_vars%rmm_desorption(:,isoil_layer,:) = scalar_value
    CASE('rtm_mic_uptake', 'sb_rtm_mic_uptake')
      sb_analytical_init_vars%rtm_mic_uptake(:,isoil_layer,:) = scalar_value
    CASE('rmm_mic_uptake', 'sb_rmm_mic_uptake')
      sb_analytical_init_vars%rmm_mic_uptake(:,isoil_layer,:) = scalar_value
    CASE('rtm_hsc', 'sb_rtm_hsc')
      sb_analytical_init_vars%rtm_hsc(:,isoil_layer,:) = scalar_value
    CASE('rmm_hsc', 'sb_rmm_hsc')
      sb_analytical_init_vars%rmm_hsc(:,isoil_layer,:) = scalar_value
    CASE('bulk_dens_sl', 'sse_bulk_dens_sl')
      sb_analytical_init_vars%bulk_dens_sl(:,isoil_layer,:) = scalar_value
    CASE('silt_sl', 'sse_silt_sl')
      sb_analytical_init_vars%silt_sl(:,isoil_layer,:) = scalar_value
    CASE('clay_sl', 'sse_clay_sl')
      sb_analytical_init_vars%clay_sl(:,isoil_layer,:) = scalar_value
    CASE('formation_soluable_litter_c')
      sb_analytical_init_vars%formation_soluable_litter_c(:,isoil_layer,:) = scalar_value
    CASE('formation_polymeric_litter_c')
      sb_analytical_init_vars%formation_polymeric_litter_c(:,isoil_layer,:) = scalar_value
    CASE('formation_woody_litter_c')
      sb_analytical_init_vars%formation_woody_litter_c(:,isoil_layer,:) = scalar_value
    CASE DEFAULT
      CALL finish(TRIM(modname)//':sb_assign_analytical_init_layer_scalar', &
        & 'Unknown analytical steady-state init variable in standalone text file: '//TRIM(var_name))
    END SELECT
  END SUBROUTINE sb_assign_analytical_init_layer_scalar

  ! ======================================================================================================= !
  !> assign one standalone analytical steady-state record to the analytical init storage
  !>
  !> The standalone text file now stores one site per row with the fixed column order:
  !>   site_name,
  !>   formation_soluable_litter_c_L1..Ln,
  !>   formation_polymeric_litter_c_L1..Ln,
  !>   formation_woody_litter_c_L1..Ln,
  !>   frac_wtr_transp_down_sl_L1..Ln,
  !>   frac_w_lat_loss_sl_L1..Ln,
  !>   rtm_decomposition_L1..Ln,
  !>   rmm_decomposition_L1..Ln,
  !>   rtm_depolymerisation_L1..Ln,
  !>   rmm_depolymerisation_L1..Ln,
  !>   rtm_mic_uptake_L1..Ln,
  !>   rmm_mic_uptake_L1..Ln,
  !>   rtm_sorption_L1..Ln,
  !>   rmm_sorption_L1..Ln,
  !>   rtm_desorption_L1..Ln,
  !>   rmm_desorption_L1..Ln,
  !>   rtm_hsc_L1..Ln,
  !>   rmm_hsc_L1..Ln,
  !>   bulk_dens_sl_L1..Ln,
  !>   silt_sl_L1..Ln,
  !>   clay_sl_L1..Ln
  !>
  SUBROUTINE sb_assign_analytical_init_record_values(record_values, nsoil_w)
    REAL(wp), INTENT(in) :: record_values(:)
    INTEGER,  INTENT(in) :: nsoil_w

    INTEGER :: idx
    INTEGER :: isoil
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_assign_analytical_init_record_values'

    IF (SIZE(record_values) /= 20 * nsoil_w) THEN
      CALL finish(TRIM(routine), 'Unexpected number of values in standalone analytical init record.')
    END IF

    idx = 0

    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%formation_soluable_litter_c(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%formation_polymeric_litter_c(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%formation_woody_litter_c(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%frac_wtr_transp_down_sl(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%frac_w_lat_loss_sl(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rtm_decomposition(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rmm_decomposition(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rtm_depolymerisation(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rmm_depolymerisation(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rtm_mic_uptake(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rmm_mic_uptake(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rtm_sorption(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rmm_sorption(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rtm_desorption(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rmm_desorption(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rtm_hsc(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%rmm_hsc(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%bulk_dens_sl(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%silt_sl(:, isoil, :) = record_values(idx)
    END DO
    DO isoil = 1, nsoil_w
      idx = idx + 1
      sb_analytical_init_vars%clay_sl(:, isoil, :) = record_values(idx)
    END DO
  END SUBROUTINE sb_assign_analytical_init_record_values

#ifdef __QUINCY_STANDALONE__
  ! ======================================================================================================= !
  !> read analytical steady-state SB_ init vars in QUINCY standalone mode
  !>
  SUBROUTINE sb_qs_read_analytical_init_vars(tile)
    USE mo_sb_constants, ONLY: min_enzyme_fraction

    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_w
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_w
    INTEGER                     :: analytical_unit
    INTEGER                     :: io_status
    INTEGER                     :: nline
    INTEGER                     :: values_per_site
    INTEGER                     :: data_lines_read
    INTEGER                     :: ivalue
    CHARACTER(len=32768)        :: line
    CHARACTER(len=32)           :: line_no_str
    CHARACTER(len=128)          :: site_name
    REAL(wp), ALLOCATABLE       :: record_values(:)
    LOGICAL                     :: header_read
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_qs_read_analytical_init_vars'

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       = vgrid_soil_w%n_levels

    dsl4jsb_Get_config(SB_)

    IF (TRIM(dsl4jsb_Config(SB_)%sb_init_scheme) /= 'analytical_ss') RETURN

    CALL sb_allocate_analytical_init_vars(nproma, nsoil_w, nblks)
    CALL sb_set_analytical_init_vars_defaults()

    !> Temporary standalone path: read analytical initialisation values from a simple text file.
    !> The same sb_analytical_init_filename config entry is used in both code paths:
    !> ICON-QUINCY points it to the gridded NetCDF file, whereas QUINCY standalone points it
    !> to a text table with one site record per row. The site-level run script extracts the
    !> header and the selected site row into the local run directory before model start, so
    !> the standalone reader here expects exactly one header row plus exactly one data row.
    !> The fixed column order is documented in sb_assign_analytical_init_record_values.
    values_per_site = 20 * nsoil_w
    ALLOCATE(record_values(values_per_site))

    OPEN(NEWUNIT=analytical_unit, FILE=TRIM(dsl4jsb_Config(SB_)%sb_analytical_init_filename), &
      & STATUS='OLD', ACTION='READ', IOSTAT=io_status)
    IF (io_status /= 0) THEN
      CALL finish(TRIM(routine), 'Cannot open standalone analytical init file: ' // &
        & TRIM(dsl4jsb_Config(SB_)%sb_analytical_init_filename))
    END IF

    nline = 0
    data_lines_read = 0
    header_read = .FALSE.
    DO
      READ(analytical_unit, '(A)', IOSTAT=io_status) line
      IF (io_status /= 0) EXIT
      nline = nline + 1

      IF (LEN_TRIM(line) == 0) CYCLE
      IF (line(1:1) == '#') CYCLE
      IF (line(1:1) == '!') CYCLE

      IF (.NOT. header_read) THEN
        READ(line, *, IOSTAT=io_status) site_name
        IF (io_status /= 0 .OR. TRIM(site_name) /= 'site_name') THEN
          CALL finish(TRIM(routine), &
          &  'Standalone analytical init file must start with a header line beginning with site_name.')
        END IF
        header_read = .TRUE.
      ELSE
        IF (data_lines_read >= 1) THEN
          CALL finish(TRIM(routine), &
          &  'Standalone analytical init file in the site run directory must contain exactly one site record.')
        END IF

        READ(line, *, IOSTAT=io_status) site_name, (record_values(ivalue), ivalue = 1, values_per_site)
        IF (io_status /= 0) THEN
          WRITE(line_no_str, '(I0)') nline
          CALL finish(TRIM(routine), 'Invalid data record in standalone analytical init file at line ' // &
            & TRIM(line_no_str) // ': ' // TRIM(line))
        END IF

        !> The selected site row already matches the current run directory, so only the
        !> ordered values remain to be unpacked into the analytical-init storage.
        CALL sb_assign_analytical_init_record_values(record_values, nsoil_w)
        data_lines_read = data_lines_read + 1
      END IF
    END DO

    CLOSE(analytical_unit)
    IF (.NOT. header_read) THEN
      CALL finish(TRIM(routine), &
      &  'Standalone analytical init file is empty or misses the header row.')
    END IF
    IF (data_lines_read /= 1) THEN
      CALL finish(TRIM(routine), &
      &  'Standalone analytical init file must contain exactly one selected site record.')
    END IF

    sb_analytical_init_vars%enzyme_frac_ap(:,:,:) = MIN(MAX(sb_analytical_init_vars%enzyme_frac_ap(:,:,:), &
      & min_enzyme_fraction), 1._wp - min_enzyme_fraction)
    sb_analytical_init_vars%is_loaded = .TRUE.
    DEALLOCATE(record_values)

    IF (debug_on()) THEN
      CALL message(TRIM(routine), 'Read standalone analytical steady-state SB init vars from: ' // &
        & TRIM(dsl4jsb_Config(SB_)%sb_analytical_init_filename))
    END IF
  END SUBROUTINE sb_qs_read_analytical_init_vars

  SUBROUTINE sb_qs_read_init_vars(tile)
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_w
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_w
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_qs_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting SB_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels

    ALLOCATE( &
      & sb_init_vars%usda_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%qmax_org_fine_particle  (nproma, nblks),       &
      & sb_init_vars%soil_ph                 (nproma, nblks),       &
      & sb_init_vars%soil_p_labile           (nproma, nblks),       &
      & sb_init_vars%soil_p_slow             (nproma, nblks),       &
      & sb_init_vars%soil_p_occluded         (nproma, nblks),       &
      & sb_init_vars%soil_p_primary          (nproma, nblks),       &
      & sb_init_vars%usda_taxonomy_class_int (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class_int (nproma, nblks)        &
      & )

    sb_init_vars%usda_taxonomy_class     (:,:) = missval
    sb_init_vars%nwrb_taxonomy_class     (:,:) = missval
    sb_init_vars%qmax_org_fine_particle  (:,:) = missval
    sb_init_vars%soil_ph                 (:,:) = missval
    sb_init_vars%soil_p_labile           (:,:) = missval
    sb_init_vars%soil_p_slow             (:,:) = missval
    sb_init_vars%soil_p_occluded         (:,:) = missval
    sb_init_vars%soil_p_primary          (:,:) = missval
    sb_init_vars%usda_taxonomy_class_int (:,:) = -9999
    sb_init_vars%nwrb_taxonomy_class_int (:,:) = -9999

    dsl4jsb_Get_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !

    !> USDA taxonomy soil class
    !>
    sb_init_vars%usda_taxonomy_class_int(:,:) = dsl4jsb_Config(SB_)%usda_taxonomy_class

    !> NWRB taxonomy soil class
    !>
    sb_init_vars%nwrb_taxonomy_class_int(:,:) = dsl4jsb_Config(SB_)%nwrb_taxonomy_class

    !> qmax organic fine particle
    !>
    sb_init_vars%qmax_org_fine_particle(:,:) = dsl4jsb_Config(SB_)%qmax_org_fine_particle

    !> soil ph
    !>
    sb_init_vars%soil_ph(:,:) = dsl4jsb_Config(SB_)%soil_ph

    !> P labile
    !>
    sb_init_vars%soil_p_labile(:,:) = dsl4jsb_Config(SB_)%soil_p_labile

    !> P slow
    !>
    sb_init_vars%soil_p_slow(:,:) = dsl4jsb_Config(SB_)%soil_p_slow

    !> P occluded
    !>
    sb_init_vars%soil_p_occluded(:,:) = dsl4jsb_Config(SB_)%soil_p_occluded

    !> P primary (apatite)
    !>
    sb_init_vars%soil_p_primary(:,:) = dsl4jsb_Config(SB_)%soil_p_primary
  END SUBROUTINE sb_qs_read_init_vars
#else
  ! ======================================================================================================= !
  !> read analytical steady-state SB_ init vars from input file
  !>
  SUBROUTINE sb_read_analytical_init_vars(tile)
    USE mo_jsb_io_netcdf,       ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_sb_constants,        ONLY: min_enzyme_fraction, max_enzyme_fraction, &
      &                               temp_freeze_window_bio, t_ref_decomposition, ea_decomposition, ed_decomposition, &
      &                               ea_depolymerisation, ea_mic_uptake, ea_sorption, ea_desorption, ea_hsc, &
      &                               kmr_decomposition, kmr_depolymerisation, kmr_mic_uptake, kmr_sorption, kmr_hsc
    USE mo_q_sb_rate_modifier,  ONLY: calc_temperature_rate_modifier, calc_oxygen_rate_modifier, &
      &                               calc_convert_bulk_to_liquid, calc_substrate_rate_modifier
    USE mo_phy_schemes,         ONLY: calc_peaked_arrhenius_function
    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SB_)

    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_w
    TYPE(t_input_file)          :: input_file
    REAL(wp), ALLOCATABLE       :: t_soil_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: wtr_soil_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: wtr_soil_pot_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: ice_soil_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: vol_field_cap_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: soil_depth_sl(:,:,:)
    REAL(wp)                    :: wtr_soil_fc_local
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_w
    INTEGER                     :: iblk
    INTEGER                     :: is
    INTEGER                     :: ic
    LOGICAL                     :: have_t_soil
    LOGICAL                     :: have_wtr_soil
    LOGICAL                     :: have_wtr_soil_pot
    LOGICAL                     :: have_ice_soil
    LOGICAL                     :: have_vol_field_cap
    LOGICAL                     :: have_soil_depth
    LOGICAL                     :: have_rtm_decomposition
    LOGICAL                     :: have_rmm_decomposition
    LOGICAL                     :: have_rtm_depolymerisation
    LOGICAL                     :: have_rmm_depolymerisation
    LOGICAL                     :: have_rtm_sorption
    LOGICAL                     :: have_rmm_sorption
    LOGICAL                     :: have_rtm_desorption
    LOGICAL                     :: have_rmm_desorption
    LOGICAL                     :: have_rtm_mic_uptake
    LOGICAL                     :: have_rmm_mic_uptake
    LOGICAL                     :: have_rtm_hsc
    LOGICAL                     :: have_rmm_hsc
    LOGICAL                     :: have_enzyme_frac_ap
    LOGICAL                     :: need_jsm_modifiers
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_read_analytical_init_vars'

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SB_)

    IF (TRIM(dsl4jsb_Config(SB_)%sb_init_scheme) /= 'analytical_ss') RETURN

    hgrid        => Get_grid(model%grid_id)
    nproma       =  hgrid%Get_nproma()
    nblks        =  hgrid%Get_nblks()
    vgrid_soil_w => Get_vgrid('soil_depth_water')
    nsoil_w      =  vgrid_soil_w%n_levels

    CALL sb_allocate_analytical_init_vars(nproma, nsoil_w, nblks)
    CALL sb_set_analytical_init_vars_defaults()

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%sb_analytical_init_filename), model%grid_id)

    !> Read the preprocessed analytical-init fields from the static NetCDF file.
    CALL sb_read_analytical_var(input_file, sb_analytical_init_vars%frac_wtr_transp_down_sl, 'frac_wtr_transp_down_sl', nsoil_w)
    CALL sb_read_analytical_var(input_file, sb_analytical_init_vars%frac_w_lat_loss_sl,    'frac_w_lat_loss_sl',    nsoil_w)
    CALL sb_read_analytical_var(input_file, sb_analytical_init_vars%formation_soluable_litter_c,  &
      & 'formation_soluable_litter_c', nsoil_w)
    CALL sb_read_analytical_var(input_file, sb_analytical_init_vars%formation_polymeric_litter_c, &
      & 'formation_polymeric_litter_c', nsoil_w)
    CALL sb_read_analytical_var(input_file, sb_analytical_init_vars%formation_woody_litter_c,     &
      & 'formation_woody_litter_c', nsoil_w)

    !> Optional direct fields: if the analytical output stream already contains them, use them directly.
    CALL sb_try_read_analytical_var(input_file, 'rtm_decomposition', &
      & sb_analytical_init_vars%rtm_decomposition, have_rtm_decomposition)
    CALL sb_try_read_analytical_var(input_file, 'rmm_decomposition', &
      & sb_analytical_init_vars%rmm_decomposition, have_rmm_decomposition)
    CALL sb_try_read_analytical_var(input_file, 'rtm_depolymerisation', &
      & sb_analytical_init_vars%rtm_depolymerisation, have_rtm_depolymerisation)
    CALL sb_try_read_analytical_var(input_file, 'rmm_depolymerisation', &
      & sb_analytical_init_vars%rmm_depolymerisation, have_rmm_depolymerisation)
    CALL sb_try_read_analytical_var(input_file, 'rtm_sorption', &
      & sb_analytical_init_vars%rtm_sorption, have_rtm_sorption)
    CALL sb_try_read_analytical_var(input_file, 'rmm_sorption', &
      & sb_analytical_init_vars%rmm_sorption, have_rmm_sorption)
    CALL sb_try_read_analytical_var(input_file, 'rtm_desorption', &
      & sb_analytical_init_vars%rtm_desorption, have_rtm_desorption)
    CALL sb_try_read_analytical_var(input_file, 'rmm_desorption', &
      & sb_analytical_init_vars%rmm_desorption, have_rmm_desorption)
    CALL sb_try_read_analytical_var(input_file, 'rtm_mic_uptake', &
      & sb_analytical_init_vars%rtm_mic_uptake, have_rtm_mic_uptake)
    CALL sb_try_read_analytical_var(input_file, 'rmm_mic_uptake', &
      & sb_analytical_init_vars%rmm_mic_uptake, have_rmm_mic_uptake)
    CALL sb_try_read_analytical_var(input_file, 'rtm_hsc', &
      & sb_analytical_init_vars%rtm_hsc, have_rtm_hsc)
    CALL sb_try_read_analytical_var(input_file, 'rmm_hsc', &
      & sb_analytical_init_vars%rmm_hsc, have_rmm_hsc)
    CALL sb_try_read_analytical_var(input_file, 'enzyme_frac_ap', &
      & sb_analytical_init_vars%enzyme_frac_ap, have_enzyme_frac_ap)

    need_jsm_modifiers = TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == 'jsm'

    !> Fall back to reconstructing the modifiers from the raw box-tile hydro and soil-temperature state.
    IF ((.NOT. have_rtm_decomposition) .OR. (.NOT. have_rmm_decomposition) .OR. &
      & (need_jsm_modifiers .AND. ((.NOT. have_rtm_depolymerisation) .OR. (.NOT. have_rmm_depolymerisation) .OR. &
      &                            (.NOT. have_rtm_sorption)         .OR. (.NOT. have_rmm_sorption)         .OR. &
      &                            (.NOT. have_rtm_desorption)       .OR. (.NOT. have_rmm_desorption)       .OR. &
      &                            (.NOT. have_rtm_mic_uptake)       .OR. (.NOT. have_rmm_mic_uptake)       .OR. &
      &                            (.NOT. have_rtm_hsc)              .OR. (.NOT. have_rmm_hsc))))) THEN

      ALLOCATE(t_soil_sl(nproma, nsoil_w, nblks))
      ALLOCATE(wtr_soil_sl(nproma, nsoil_w, nblks))
      ALLOCATE(wtr_soil_pot_sl(nproma, nsoil_w, nblks))
      ALLOCATE(ice_soil_sl(nproma, nsoil_w, nblks))
      ALLOCATE(vol_field_cap_sl(nproma, nsoil_w, nblks))
      ALLOCATE(soil_depth_sl(nproma, nsoil_w, nblks))

      CALL sb_try_read_analytical_var(input_file, 't_soil_sl', t_soil_sl, have_t_soil)
      CALL sb_try_read_analytical_var(input_file, 'wtr_soil_sl', wtr_soil_sl, have_wtr_soil)
      CALL sb_try_read_analytical_var(input_file, 'wtr_soil_pot_sl', wtr_soil_pot_sl, have_wtr_soil_pot)
      CALL sb_try_read_analytical_var(input_file, 'ice_soil_sl', ice_soil_sl, have_ice_soil)
      CALL sb_try_read_analytical_var(input_file, 'vol_field_cap_sl', vol_field_cap_sl, have_vol_field_cap)
      CALL sb_try_read_analytical_var(input_file, 'soil_depth_sl', soil_depth_sl, have_soil_depth)

      IF ((.NOT. have_rtm_decomposition .OR. .NOT. have_rmm_decomposition) .AND. &
        & (.NOT. (have_t_soil .AND. have_wtr_soil_pot))) THEN
        CALL finish(TRIM(routine), 'Need either direct analytical decomposition modifiers ' // &
          & '(sb_rtm_decomposition_box, sb_rmm_decomposition_box) or raw state variables ' // &
          & '(sse_t_soil_sl_box, hydro_wtr_soil_pot_sl_box) in ' // &
          & TRIM(dsl4jsb_Config(SB_)%sb_analytical_init_filename))
      END IF

      IF (need_jsm_modifiers .AND. ((.NOT. have_rtm_depolymerisation) .OR. (.NOT. have_rmm_depolymerisation) .OR. &
        &                            (.NOT. have_rtm_sorption)         .OR. (.NOT. have_rmm_sorption)         .OR. &
        &                            (.NOT. have_rtm_desorption)       .OR. (.NOT. have_rmm_desorption)       .OR. &
        &                            (.NOT. have_rtm_mic_uptake)       .OR. (.NOT. have_rmm_mic_uptake)       .OR. &
        &                            (.NOT. have_rtm_hsc)              .OR. (.NOT. have_rmm_hsc)) .AND. &
        &      (.NOT. (have_t_soil .AND. have_wtr_soil .AND. have_ice_soil .AND. have_vol_field_cap .AND. have_soil_depth)))) THEN
        CALL finish(TRIM(routine), 'Need either direct analytical JSM modifiers ' // &
          & '(sb_rtm_*/sb_rmm_*_box) or raw state variables ' // &
          & '(sse_t_soil_sl_box, hydro_wtr_soil_sl_box, hydro_ice_soil_sl_box, ' // &
          & 'hydro_vol_field_cap_sl_box, hydro_soil_depth_sl_box) in ' // &
          & TRIM(dsl4jsb_Config(SB_)%sb_analytical_init_filename))
      END IF

      DO iblk = 1, nblks
        DO ic = 1, nproma
          DO is = 1, nsoil_w
            IF (.NOT. have_rtm_decomposition) THEN
              sb_analytical_init_vars%rtm_decomposition(ic, is, iblk) = calc_peaked_arrhenius_function( &
                & t_soil_sl(ic, is, iblk), ea_decomposition, ed_decomposition, t_ref_decomposition)
            END IF
            IF (.NOT. have_rmm_decomposition) THEN
              sb_analytical_init_vars%rmm_decomposition(ic, is, iblk) = MAX(1._wp - &
                & wtr_soil_pot_sl(ic, is, iblk) / kmr_decomposition, 0.0_wp)
            END IF

            IF (need_jsm_modifiers) THEN
              wtr_soil_fc_local = vol_field_cap_sl(ic, is, iblk) * soil_depth_sl(ic, is, iblk)

              IF (.NOT. have_rtm_depolymerisation) THEN
                sb_analytical_init_vars%rtm_depolymerisation(ic, is, iblk) = calc_temperature_rate_modifier( &
                  & t_soil_sl(ic, is, iblk), ea_depolymerisation, temp_freeze_window_bio)
              END IF
              IF (.NOT. have_rtm_mic_uptake) THEN
                sb_analytical_init_vars%rtm_mic_uptake(ic, is, iblk) = calc_temperature_rate_modifier( &
                  & t_soil_sl(ic, is, iblk), ea_mic_uptake, temp_freeze_window_bio)
              END IF
              IF (.NOT. have_rtm_sorption) THEN
                sb_analytical_init_vars%rtm_sorption(ic, is, iblk) = calc_temperature_rate_modifier( &
                  & t_soil_sl(ic, is, iblk), ea_sorption)
              END IF
              IF (.NOT. have_rtm_desorption) THEN
                sb_analytical_init_vars%rtm_desorption(ic, is, iblk) = calc_temperature_rate_modifier( &
                  & t_soil_sl(ic, is, iblk), ea_desorption)
              END IF
              IF (.NOT. have_rtm_hsc) THEN
                sb_analytical_init_vars%rtm_hsc(ic, is, iblk) = calc_temperature_rate_modifier( &
                  & t_soil_sl(ic, is, iblk), ea_hsc)
              END IF

              IF (.NOT. have_rmm_depolymerisation) THEN
                sb_analytical_init_vars%rmm_depolymerisation(ic, is, iblk) = calc_oxygen_rate_modifier( &
                  & wtr_soil_sl(ic, is, iblk), ice_soil_sl(ic, is, iblk), wtr_soil_fc_local, kmr_depolymerisation)
              END IF
              IF (.NOT. have_rmm_mic_uptake) THEN
                sb_analytical_init_vars%rmm_mic_uptake(ic, is, iblk) = calc_oxygen_rate_modifier( &
                  & wtr_soil_sl(ic, is, iblk), ice_soil_sl(ic, is, iblk), wtr_soil_fc_local, kmr_mic_uptake)
              END IF
              IF (.NOT. have_rmm_sorption) THEN
                sb_analytical_init_vars%rmm_sorption(ic, is, iblk) = calc_convert_bulk_to_liquid( &
                  & wtr_soil_sl(ic, is, iblk), soil_depth_sl(ic, is, iblk), kmr_sorption)
              END IF
              IF (.NOT. have_rmm_desorption) THEN
                sb_analytical_init_vars%rmm_desorption(ic, is, iblk) = calc_convert_bulk_to_liquid( &
                  & wtr_soil_sl(ic, is, iblk), soil_depth_sl(ic, is, iblk), kmr_sorption)
              END IF
              IF (.NOT. have_rmm_hsc) THEN
                sb_analytical_init_vars%rmm_hsc(ic, is, iblk) = calc_substrate_rate_modifier( &
                  & wtr_soil_sl(ic, is, iblk), wtr_soil_fc_local, kmr_hsc)
              END IF
            END IF
          END DO
        END DO
      END DO

      DEALLOCATE(t_soil_sl, wtr_soil_sl, wtr_soil_pot_sl, ice_soil_sl, vol_field_cap_sl, soil_depth_sl)
    END IF

    IF (need_jsm_modifiers .AND. (.NOT. have_enzyme_frac_ap)) THEN
      CALL message(TRIM(routine), 'enzyme_frac_ap not found in analytical init file; using min_enzyme_fraction as fallback.')
      sb_analytical_init_vars%enzyme_frac_ap(:,:,:) = min_enzyme_fraction
    ELSE
      sb_analytical_init_vars%enzyme_frac_ap(:,:,:) = MIN(MAX(sb_analytical_init_vars%enzyme_frac_ap(:,:,:), &
        & min_enzyme_fraction), max_enzyme_fraction)
    END IF

    sb_analytical_init_vars%is_loaded = .TRUE.
    CALL input_file%Close()

    CALL message(TRIM(routine), 'Prepared analytical steady-state SB initialisation from preprocessed file: ' // &
      & TRIM(dsl4jsb_Config(SB_)%sb_analytical_init_filename))
  END SUBROUTINE sb_read_analytical_init_vars

  ! ======================================================================================================= !
  !> read a layer-resolved analytical-init variable from a static file
  !>
  SUBROUTINE sb_read_analytical_var(input_file, ptr3d, var_name, nsoil_w)
    USE mo_jsb_io_netcdf, ONLY: t_input_file

    TYPE(t_input_file), INTENT(inout) :: input_file
    REAL(wp),           INTENT(inout) :: ptr3d(:,:,:)
    CHARACTER(*),       INTENT(in)    :: var_name
    INTEGER,            INTENT(in)    :: nsoil_w
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_read_analytical_var'
    INTEGER                     :: i_soil
    REAL(wp), POINTER           :: return_3d(:,:,:)

    IF (.NOT. input_file%Has_var(var_name)) THEN
      CALL finish(TRIM(routine), 'Missing required analytical-init variable '//TRIM(var_name)//' in file.')
    END IF
    DO i_soil = 1, nsoil_w
      return_3d => input_file%Read_2d_extdim(variable_name=var_name, &
        & start_extdim=i_soil, end_extdim=i_soil, extdim_name='soillev')
      ptr3d(:, i_soil, :) = return_3d(:,:,1)
    END DO
  END SUBROUTINE sb_read_analytical_var

  ! ======================================================================================================= !
  !> read a layer-resolved output variable if it exists
  !>
  SUBROUTINE sb_try_read_analytical_var(input_file, var_name, ptr3d, found)
    USE mo_jsb_io_netcdf, ONLY: t_input_file

    TYPE(t_input_file), INTENT(inout) :: input_file
    CHARACTER(*),       INTENT(in)    :: var_name
    REAL(wp),           INTENT(inout) :: ptr3d(:,:,:)
    LOGICAL,            INTENT(out)   :: found
    INTEGER                         :: i_soil
    REAL(wp), POINTER               :: return_3d(:,:,:)

    found = input_file%Has_var(var_name)
    IF (found) THEN
      DO i_soil = 1, SIZE(ptr3d, 2)
        return_3d => input_file%Read_2d_extdim(variable_name=var_name, &
          & start_extdim=i_soil, end_extdim=i_soil, extdim_name='soillev')
        ptr3d(:, i_soil, :) = return_3d(:,:,1)
      END DO
    END IF
  END SUBROUTINE sb_try_read_analytical_var

  ! ======================================================================================================= !
  !> read SB_ init vars from input file
  !>
  SUBROUTINE sb_read_init_vars(tile)
    USE mo_jsb_io_netcdf,      ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_io,             ONLY: missval

    dsl4jsb_Use_config(SB_)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile

    dsl4jsb_Def_config(SB_)

    REAL(wp), POINTER           :: ptr_2D(:, :)
    TYPE(t_jsb_model), POINTER  :: model
    TYPE(t_jsb_grid),  POINTER  :: hgrid
    TYPE(t_jsb_vgrid), POINTER  :: vgrid_soil_w
    TYPE(t_input_file)          :: input_file
    INTEGER                     :: nproma
    INTEGER                     :: nblks
    INTEGER                     :: nsoil_w
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_read_init_vars'

    IF (debug_on()) CALL message(routine, 'Reading/setting SB_ init vars')

    model  => Get_model(tile%owner_model_id)
    hgrid  => Get_grid(model%grid_id)
    nproma = hgrid%Get_nproma()
    nblks  = hgrid%Get_nblks()
    vgrid_soil_w  => Get_vgrid('soil_depth_water')
    nsoil_w       =  vgrid_soil_w%n_levels

    ALLOCATE( &
      & sb_init_vars%usda_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class     (nproma, nblks),       &
      & sb_init_vars%qmax_org_fine_particle  (nproma, nblks),       &
      & sb_init_vars%soil_ph                 (nproma, nblks),       &
      & sb_init_vars%soil_p_labile           (nproma, nblks),       &
      & sb_init_vars%soil_p_slow             (nproma, nblks),       &
      & sb_init_vars%soil_p_occluded         (nproma, nblks),       &
      & sb_init_vars%soil_p_primary          (nproma, nblks),       &
      & sb_init_vars%usda_taxonomy_class_int (nproma, nblks),       &
      & sb_init_vars%nwrb_taxonomy_class_int (nproma, nblks)        &
      & )

    sb_init_vars%usda_taxonomy_class     (:,:) = missval
    sb_init_vars%nwrb_taxonomy_class     (:,:) = missval
    sb_init_vars%qmax_org_fine_particle  (:,:) = missval
    sb_init_vars%soil_ph                 (:,:) = missval
    sb_init_vars%soil_p_labile           (:,:) = missval
    sb_init_vars%soil_p_slow             (:,:) = missval
    sb_init_vars%soil_p_occluded         (:,:) = missval
    sb_init_vars%soil_p_primary          (:,:) = missval
    sb_init_vars%usda_taxonomy_class_int (:,:) = -9999
    sb_init_vars%nwrb_taxonomy_class_int (:,:) = -9999

    dsl4jsb_Get_config(SB_)
    ! ----------------------------------------------------------------------------------------------------- !

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%bc_quincy_soil_filename), model%grid_id)


    !> USDA taxonomy soil class
    !>
    ptr_2D => input_file%Read_2d(       &
      & variable_name='TAXOUSDA_10km',  &
      & fill_array = sb_init_vars%usda_taxonomy_class)
    sb_init_vars%usda_taxonomy_class(:,:)     = MERGE(ptr_2D, 1.0_wp, ptr_2D > 1.0_wp)
    sb_init_vars%usda_taxonomy_class_int(:,:) = INT(sb_init_vars%usda_taxonomy_class(:,:))

    !> NWRB taxonomy soil class
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='TAXNWRB',   &
      & fill_array = sb_init_vars%nwrb_taxonomy_class)
    sb_init_vars%nwrb_taxonomy_class(:,:)     = MERGE(ptr_2D, 1.0_wp, ptr_2D > 1.0_wp)
    sb_init_vars%nwrb_taxonomy_class_int(:,:) = INT(sb_init_vars%nwrb_taxonomy_class(:,:))

    !> qmax organic fine particle
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='qmax_org',  &
      & fill_array = sb_init_vars%qmax_org_fine_particle)
    sb_init_vars%qmax_org_fine_particle(:,:) = MERGE(ptr_2D, 6.5537_wp, ptr_2D > 0.0_wp)

    !> soil ph
    !>
    ptr_2D => input_file%Read_2d(   &
      & variable_name='PHIHOX',     &
      & fill_array = sb_init_vars%soil_ph)
    sb_init_vars%soil_ph(:,:) = MERGE(ptr_2D, 65.0_wp, ptr_2D > 0.0_wp)
    ! divide by 10 to get the actual values within the correct range
    sb_init_vars%soil_ph(:,:) = sb_init_vars%soil_ph(:,:) / 10.0_wp

    !> P labile
    !>
    ptr_2D => input_file%Read_2d(           &
      & variable_name='Labile_Inorganic_P', &
      & fill_array = sb_init_vars%soil_p_labile)
    sb_init_vars%soil_p_labile(:,:) = MERGE(ptr_2D, 50._wp, ptr_2D > 0.0_wp)

    !> P slow
    !>
    ptr_2D => input_file%Read_2d(             &
      & variable_name='Seconday_Mineral__P',  &
      & fill_array = sb_init_vars%soil_p_slow)
    sb_init_vars%soil_p_slow(:,:) = MERGE(ptr_2D, 70._wp, ptr_2D > 0.0_wp)

    !> P occluded
    !>
    ptr_2D => input_file%Read_2d(     &
      & variable_name='Occluded_P',   &
      & fill_array = sb_init_vars%soil_p_occluded)
    sb_init_vars%soil_p_occluded(:,:) = MERGE(ptr_2D, 150._wp, ptr_2D > 0.0_wp)

    !> P primary (apatite)
    !>
    ptr_2D => input_file%Read_2d(  &
      & variable_name='Apatite_P',  &
      & fill_array = sb_init_vars%soil_p_primary)
    sb_init_vars%soil_p_primary(:,:) = MERGE(ptr_2D, 100._wp, ptr_2D > 0.0_wp)

    CALL input_file%Close()
  END SUBROUTINE sb_read_init_vars
#endif

  ! ======================================================================================================= !
  !> calculate analytical steady-state carbon pools for the simple soil model
  !> The shared SB pool indices map to the simple-model notation as
  !>   C_met  <-> ix_soluable_litter
  !>   C_str  <-> ix_polymeric_litter
  !>   C_fast <-> ix_microbial
  !>   C_slow <-> ix_residue
  !> while the other carbon pools remain unused in the simple_1d analytical initialisation.
  !> The current C-only analytical draft uses
  !>   C_wl*   = I_wl * tau_wl
  !>   C_met*  = tau_met * (I_met + eta_met^wl * I_wl)
  !>   C_str*  = tau_str * (I_str + eta_str^wl * I_wl)
  !>   C_fast* = (tau_fast * eta_fast^l * A) / (1 - eta_fast^slow * eta_slow^fast)
  !>   C_slow* = (tau_slow * eta_slow^fast * eta_fast^l * A) / (1 - eta_fast^slow * eta_slow^fast)
  !>   A       = I_met + I_str + I_wl * (eta_met^wl + eta_str^wl)
  !> with the same decomposition modifiers as the prognostic simple-model path, i.e.
  !> every turnover time is divided by rtm_decomposition * rmm_decomposition.
  !> using the simple-model routing constants from mo_sb_constants.
  !>
  !> Unit convention of the analytical status file / sb_analytical_init_vars used here:
  !>   I_met, I_str, I_wl         [mol C m-3 s-1]
  !>   rtm_*, rmm_*               [-]
  !>   tau_*                      [s]
  !>   C_*                        [mol C m-3]
  !> Therefore the standalone converter must preserve litter-input rates when it aggregates
  !> timestep output into sb_analytical_init.dat.
  !>
  SUBROUTINE calc_analytical_ss_simple_sm( &
    & nproma, &
    & nsoil_w, &
    & nblks, &
    & num_sl_above_bedrock, &
    & sb_pool_mt_domain)

    USE mo_sb_constants, ONLY: tau_metabolic_litter, tau_structural_litter, tau_woody_litter, tau_fast, tau_slow, &
      &                        frac_woody2dec_litter, frac_litter2fast_som, frac_fast2slow_som, frac_slow2fast_som, &
      &                        fc_soluable_max

    INTEGER,  INTENT(in)    :: nproma
    INTEGER,  INTENT(in)    :: nsoil_w
    INTEGER,  INTENT(in)    :: nblks
    REAL(wp), INTENT(in)    :: num_sl_above_bedrock(:,:)
    REAL(wp), INTENT(inout) :: sb_pool_mt_domain(:,:,:,:,:)

    REAL(wp)                    :: input_met
    REAL(wp)                    :: input_str
    REAL(wp)                    :: input_wood
    REAL(wp)                    :: eta_met_wood
    REAL(wp)                    :: eta_str_wood
    REAL(wp)                    :: litter_input_total
    REAL(wp)                    :: denominator
    REAL(wp)                    :: decomp_modifier
    INTEGER                     :: iblk, is, ic

    !> In the C-only analytical simple-model draft, the woody split into metabolic versus structural
    !> litter cannot yet use the N-dependent partitioning from calc_litter_formation_from_woody_decay.
    !> Use the simple-model carbon partition implied by fc_soluable_max together with the prescribed
    !> woody-to-decomposable-litter fraction until the richer status-file inputs are added.
    eta_met_wood = frac_woody2dec_litter * fc_soluable_max
    eta_str_wood = frac_woody2dec_litter * (1._wp - fc_soluable_max)

    sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_microbial, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_residue, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_dom, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) = 0.0_wp

    denominator = 1._wp - frac_fast2slow_som * frac_slow2fast_som
    IF (denominator <= eps8) THEN
      CALL finish(TRIM(modname)//':calc_analytical_ss_simple_sm', &
        & 'Analytical simple_1d initialisation requires 1 - frac_fast2slow_som * frac_slow2fast_som > 0.')
    END IF

    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          input_met  = MAX(0.0_wp, sb_analytical_init_vars%formation_soluable_litter_c(ic, is, iblk))
          input_str  = MAX(0.0_wp, sb_analytical_init_vars%formation_polymeric_litter_c(ic, is, iblk))
          input_wood = MAX(0.0_wp, sb_analytical_init_vars%formation_woody_litter_c(ic, is, iblk))
          decomp_modifier = MAX(eps8, sb_analytical_init_vars%rtm_decomposition(ic, is, iblk) * &
            &                         sb_analytical_init_vars%rmm_decomposition(ic, is, iblk))

          !> C_wl* = I_wl [mol C m-3 s-1] * tau_wl [s] / (rtm_decomposition * rmm_decomposition) [-]
          sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk) = input_wood * tau_woody_litter / decomp_modifier

          !> C_met* = tau_met [s] * (I_met + eta_met^wl * I_wl) [mol C m-3 s-1] / (rtm_decomposition * rmm_decomposition) [-]
          sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) = &
            & tau_metabolic_litter * (input_met + eta_met_wood * input_wood) / decomp_modifier

          !> C_str* = tau_str [s] * (I_str + eta_str^wl * I_wl) [mol C m-3 s-1] / (rtm_decomposition * rmm_decomposition) [-]
          sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = &
            & tau_structural_litter * (input_str + eta_str_wood * input_wood) / decomp_modifier

          !> A = I_met + I_str + I_wl * (eta_met^wl + eta_str^wl) [mol C m-3 s-1]
          litter_input_total = input_met + input_str + input_wood * (eta_met_wood + eta_str_wood)

          !> C_fast* = (tau_fast * eta_fast^l * A) &
          !>           / ((rtm_decomposition * rmm_decomposition) * (1 - eta_fast^slow * eta_slow^fast))
          sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) = &
            & tau_fast * frac_litter2fast_som * litter_input_total / (decomp_modifier * denominator)

          !> C_slow* = (tau_slow * eta_slow^fast * eta_fast^l * A) &
          !>           / ((rtm_decomposition * rmm_decomposition) * (1 - eta_fast^slow * eta_slow^fast))
          sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) = &
            & tau_slow * frac_slow2fast_som * frac_litter2fast_som * litter_input_total / (decomp_modifier * denominator)
        END DO
      END DO
    END DO

  END SUBROUTINE calc_analytical_ss_simple_sm

  ! ======================================================================================================= !
  !> temporary debug dump of analytical steady-state carbon pools
  !>
  SUBROUTINE dump_analytical_ss_pools(filename, nproma, nsoil_w, nblks, num_sl_above_bedrock, sb_pool_mt_domain)

    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER,          INTENT(in) :: nproma
    INTEGER,          INTENT(in) :: nsoil_w
    INTEGER,          INTENT(in) :: nblks
    REAL(wp),         INTENT(in) :: num_sl_above_bedrock(:,:)
    REAL(wp),         INTENT(in) :: sb_pool_mt_domain(:,:,:,:,:)

    INTEGER :: lun
    INTEGER :: ios
    INTEGER :: iblk, ic, is

    OPEN(newunit=lun, file=TRIM(filename), status='replace', action='write', iostat=ios)
    IF (ios /= 0) THEN
      CALL finish(TRIM(modname)//':dump_analytical_ss_pools', &
        & 'Could not open analytical steady-state debug output file: '//TRIM(filename))
    END IF

    WRITE(lun,'(A)') 'iblk ic is c_soluable_litter c_polymeric_litter c_woody_litter &
      &                 c_microbial c_residue c_dom c_dom_assoc c_residue_assoc &
      &                 n_soluable_litter n_polymeric_litter n_woody_litter &
      &                 n_microbial n_residue n_dom n_dom_assoc n_residue_assoc &
      &                 p_soluable_litter p_polymeric_litter p_woody_litter &
      &                 p_microbial p_residue p_dom p_dom_assoc p_residue_assoc'

    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          WRITE(lun,'(3(I6,1X),24(ES24.16E3,1X))') iblk, ic, is, &
            & sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk), &
            & sb_pool_mt_domain(ix_soluable_litter, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_polymeric_litter, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_woody_litter, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_microbial, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_residue, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_dom, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_dom_assoc, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_residue_assoc, ixN, ic, is, iblk), &
            & sb_pool_mt_domain(ix_soluable_litter, ixP, ic, is, iblk), &
            & sb_pool_mt_domain(ix_polymeric_litter, ixP, ic, is, iblk), &
            & sb_pool_mt_domain(ix_woody_litter, ixP, ic, is, iblk), &
            & sb_pool_mt_domain(ix_microbial, ixP, ic, is, iblk), &
            & sb_pool_mt_domain(ix_residue, ixP, ic, is, iblk), &
            & sb_pool_mt_domain(ix_dom, ixP, ic, is, iblk), &
            & sb_pool_mt_domain(ix_dom_assoc, ixP, ic, is, iblk), &
            & sb_pool_mt_domain(ix_residue_assoc, ixP, ic, is, iblk)
        END DO
      END DO
    END DO

    CLOSE(lun)
  END SUBROUTINE dump_analytical_ss_pools

  ! ======================================================================================================= !
  !> calculate the layer-wise maximum bulk carbon allowed by the prescribed minimum mineral-soil volume
  !>
  SUBROUTINE calc_bulk_carbon_max_profile(nproma, nsoil_w, nblks, num_sl_above_bedrock, soil_lay_depth_center_sl, bulk_carbon_max)

    USE mo_sb_constants,           ONLY: min_volume_min_surface_analytical_ss, min_volume_min_gradient_analytical_ss, &
      &                                  carbon_per_dryweight_SOM, rho_bulk_org
    USE mo_jsb_physical_constants, ONLY: molar_mass_C

    INTEGER,  INTENT(in)  :: nproma
    INTEGER,  INTENT(in)  :: nsoil_w
    INTEGER,  INTENT(in)  :: nblks
    REAL(wp), INTENT(in)  :: num_sl_above_bedrock(:,:)
    REAL(wp), INTENT(in)  :: soil_lay_depth_center_sl(:,:,:)
    REAL(wp), INTENT(out) :: bulk_carbon_max(:,:,:)

    REAL(wp) :: qmax_b
    REAL(wp) :: min_volume_min_layer
    REAL(wp) :: min_volume_min_prev
    REAL(wp) :: layer_center_distance
    INTEGER  :: iblk, ic, is

    qmax_b = carbon_per_dryweight_SOM * rho_bulk_org * 1000._wp / molar_mass_C
    bulk_carbon_max(:,:,:) = 0.0_wp

    DO iblk = 1,nblks
      DO ic = 1,nproma
        min_volume_min_prev = min_volume_min_surface_analytical_ss
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          IF (is == 1) THEN
            min_volume_min_layer = min_volume_min_surface_analytical_ss
          ELSE
            layer_center_distance = MAX(0.0_wp, soil_lay_depth_center_sl(ic, is, iblk) - &
              &                              soil_lay_depth_center_sl(ic, is - 1, iblk))
            min_volume_min_layer = min_volume_min_prev + min_volume_min_gradient_analytical_ss * layer_center_distance
          END IF
          bulk_carbon_max(ic, is, iblk) = qmax_b * MAX(0.0_wp, 1.0_wp - min_volume_min_layer)
          min_volume_min_prev = min_volume_min_layer
        END DO
      END DO
    END DO

  END SUBROUTINE calc_bulk_carbon_max_profile

  ! ======================================================================================================= !
  !> evaluate one JSM analytical steady-state layer for a prescribed microbial biomass
  !> The helper is used by the steady-state-only physical correction to project the
  !> biogeochemical solution onto an admissible soil-volume state without changing the
  !> forward model equations.
  !>
  SUBROUTINE calc_analytical_ss_jsm_layer_state( &
    & l_surface_layer, &
    & input_sol, &
    & input_poly, &
    & input_wood, &
    & dom_carbon, &
    & microbial_carbon, &
    & enzyme_frac_ap_layer, &
    & enzyme_frac_total_scale, &
    & km_depoly_eff, &
    & vmax_poly_eff, &
    & vmax_res_eff, &
    & qmax_org_min_local, &
    & bulk_dens_local, &
    & rtm_sorb_eff, &
    & rmm_sorb_eff, &
    & rtm_desorb_eff, &
    & rmm_desorb_eff, &
    & alpha_poly, &
    & alpha_residue, &
    & polymeric_carbon, &
    & residue_carbon, &
    & dom_assoc_carbon, &
    & residue_assoc_carbon, &
    & qmax_org_layer, &
    & bulk_carbon_layer)

    USE mo_q_sb_jsm_processes,      ONLY: calc_enzyme_allocation
    USE mo_sb_constants
    USE mo_jsb_physical_constants,  ONLY: molar_mass_C

    LOGICAL,  INTENT(in)  :: l_surface_layer
    REAL(wp), INTENT(in)  :: input_sol
    REAL(wp), INTENT(in)  :: input_poly
    REAL(wp), INTENT(in)  :: input_wood
    REAL(wp), INTENT(in)  :: dom_carbon
    REAL(wp), INTENT(in)  :: microbial_carbon
    REAL(wp), INTENT(in)  :: enzyme_frac_ap_layer
    REAL(wp), INTENT(in)  :: enzyme_frac_total_scale
    REAL(wp), INTENT(in)  :: km_depoly_eff
    REAL(wp), INTENT(in)  :: vmax_poly_eff
    REAL(wp), INTENT(in)  :: vmax_res_eff
    REAL(wp), INTENT(in)  :: qmax_org_min_local
    REAL(wp), INTENT(in)  :: bulk_dens_local
    REAL(wp), INTENT(in)  :: rtm_sorb_eff
    REAL(wp), INTENT(in)  :: rmm_sorb_eff
    REAL(wp), INTENT(in)  :: rtm_desorb_eff
    REAL(wp), INTENT(in)  :: rmm_desorb_eff
    REAL(wp), INTENT(out) :: alpha_poly
    REAL(wp), INTENT(out) :: alpha_residue
    REAL(wp), INTENT(out) :: polymeric_carbon
    REAL(wp), INTENT(out) :: residue_carbon
    REAL(wp), INTENT(out) :: dom_assoc_carbon
    REAL(wp), INTENT(out) :: residue_assoc_carbon
    REAL(wp), INTENT(out) :: qmax_org_layer
    REAL(wp), INTENT(out) :: bulk_carbon_layer

    REAL(wp) :: soluable_carbon
    REAL(wp) :: woody_carbon
    REAL(wp) :: enzyme_depoly
    REAL(wp) :: enzyme_frac_total_local
    REAL(wp) :: d_litter
    REAL(wp) :: d_residue
    REAL(wp) :: qmax_a
    REAL(wp) :: qmax_b
    REAL(wp) :: c0_star
    REAL(wp) :: ratio_dom
    REAL(wp) :: ratio_res
    REAL(wp) :: qmax_f
    REAL(wp) :: sorption_denom

    !> C_wl*  = I_wl [mol C m-3 s-1] * tau_wl [s]
    woody_carbon    = input_wood * tau_woody_litter
    !> C_sol* = I_sol [mol C m-3 s-1] * tau_sol [s]
    soluable_carbon = input_sol  * tau_soluable_litter

    enzyme_frac_total_local = MAX(0.0_wp, enzyme_frac_total * enzyme_frac_total_scale)
    !> E_depoly* = C_B* [mol C m-3] * enzyme_frac_total [-] * (1 - enzyme_frac_AP) [-]
    enzyme_depoly   = microbial_carbon * enzyme_frac_total_local * (1._wp - enzyme_frac_ap_layer)
    !> d_L = I_poly + eta_poly^wl * I_wl  [mol C m-3 s-1]
    d_litter        = input_poly + fc_woody2polymeric * input_wood
    !> d_R = (1 - f_soluable_residue) * k_B * C_B*  [mol C m-3 s-1]
    d_residue       = (1._wp - f_soluable_residue) * k_mort_microbial * microbial_carbon

    IF (enzyme_depoly <= enzyme_frac_total_local * min_mic_biomass .OR. (d_litter <= eps8 .AND. d_residue <= eps8)) THEN
      alpha_poly = init_enzyme_fraction
    ELSE IF (ABS(d_litter - d_residue) <= eps8 * MAX(1._wp, d_litter, d_residue)) THEN
      alpha_poly = 0.5_wp
    ELSE
      alpha_poly = calc_enzyme_allocation(d_litter, d_residue, km_depoly_eff, enzyme_depoly)
    END IF

    IF (alpha_poly /= alpha_poly) alpha_poly = init_enzyme_fraction
    alpha_poly    = MIN(MAX(alpha_poly, min_enzyme_fraction), 1._wp - min_enzyme_fraction)
    alpha_residue = 1._wp - alpha_poly

    !> C_poly*=(I_poly + eta_poly^wl * I_wl) [mol C m-3 s-1]
    !>         / (vmax_poly [s-1] * alpha_poly [-] * E_depoly* [mol C m-3])
    !>         * (Km_depoly [mol C m-3] + alpha_poly * E_depoly*)
    IF (d_litter > 0.0_wp .AND. alpha_poly * enzyme_depoly > 0.0_wp .AND. vmax_poly_eff > 0.0_wp) THEN
      polymeric_carbon = d_litter * (km_depoly_eff + alpha_poly * enzyme_depoly) / &
        &                (vmax_poly_eff * alpha_poly * enzyme_depoly)
    ELSE
      polymeric_carbon = eps8
    END IF

    !> C_res*=((1 - f_soluable_residue) * k_B * C_B*) [mol C m-3 s-1]
    !>        / (vmax_res [s-1] * alpha_residue [-] * E_depoly* [mol C m-3])
    !>        * (Km_depoly [mol C m-3] + alpha_residue * E_depoly*)
    IF (d_residue > 0.0_wp .AND. alpha_residue * enzyme_depoly > 0.0_wp .AND. vmax_res_eff > 0.0_wp) THEN
      residue_carbon = d_residue * (km_depoly_eff + alpha_residue * enzyme_depoly) / &
        &             (vmax_res_eff * alpha_residue * enzyme_depoly)
    ELSE
      residue_carbon = eps8
    END IF

    !> A = qmax_min [mol C kg-1] * rho_bulk [kg m-3]
    qmax_a = qmax_org_min_local * bulk_dens_local
    !> B = f_C^SOM [-] * rho_org [kg SOM m-3], converted to [mol C m-3]
    qmax_b = carbon_per_dryweight_SOM * rho_bulk_org * 1000._wp / molar_mass_C

    IF (l_surface_layer) THEN
      !> Keep the first-layer physical constraint aligned with the mineral-volume treatment:
      !> woody litter stays above the mineral soil matrix, while soluable litter, polymeric litter,
      !> microbial biomass, DOM, and residue are the non-associated surface SOM that can be buried into
      !> the first layer.
      !> C0*_surface = C_sol* + C_poly* + C_D* + C_B* + C_res*
      c0_star = soluable_carbon + polymeric_carbon + dom_carbon + microbial_carbon + residue_carbon
    ELSE
      !> C0* = C_sol* + C_poly* + C_D* + C_B* + C_res*
      c0_star = soluable_carbon + polymeric_carbon + dom_carbon + microbial_carbon + residue_carbon
    END IF

    !> r_D = k_ads,D / k_des,D
    ratio_dom = k_adsorpt_dom * rtm_sorb_eff * rmm_sorb_eff / (k_desorpt_dom * rtm_desorb_eff * rmm_desorb_eff)
    !> r_R = k_ads,R / k_des,R
    ratio_res = k_adsorpt_det * rtm_sorb_eff * rmm_sorb_eff / (k_desorpt_det * rtm_desorb_eff * rmm_desorb_eff)
    !> F = (r_D * C_D* + r_R * C_res*) / (1 + r_D * C_D* + r_R * C_res*)
    qmax_f = (ratio_dom * dom_carbon + ratio_res * residue_carbon) / &
      &      (1._wp + ratio_dom * dom_carbon + ratio_res * residue_carbon)

    !> Qmax* = (A * (1 - C0* / B)) / (1 + A * F / B)
    IF (qmax_b > eps8) THEN
      qmax_org_layer = MAX(0.0_wp, qmax_a * (1._wp - c0_star / qmax_b) / (1._wp + qmax_a * qmax_f / qmax_b))
    ELSE
      qmax_org_layer = eps8
    END IF

    sorption_denom = 1._wp + ratio_res * residue_carbon + ratio_dom * dom_carbon
    IF (sorption_denom > eps8) THEN
      !> C_aDOM* = (r_D * C_D*) * Qmax* / (1 + r_R * C_res* + r_D * C_D*)
      dom_assoc_carbon     = ratio_dom * dom_carbon     * qmax_org_layer / sorption_denom
      !> C_aRes* = (r_R * C_res*) * Qmax* / (1 + r_R * C_res* + r_D * C_D*)
      residue_assoc_carbon = ratio_res * residue_carbon * qmax_org_layer / sorption_denom
    ELSE
      dom_assoc_carbon     = eps8
      residue_assoc_carbon = eps8
    END IF

    IF (l_surface_layer) THEN
      !> C_bulk* for the first mineral layer includes minerally associated SOM fully,
      !> and all other carbon pools as potentially buried surface SOM.
      bulk_carbon_layer = soluable_carbon + polymeric_carbon + microbial_carbon + dom_carbon &
                          + residue_carbon + dom_assoc_carbon + residue_assoc_carbon
    ELSE
      !> C_bulk* = C_wl* + C_sol* + C_poly* + C_D* + C_B* + C_res* + C_aDOM* + C_aRes*
      bulk_carbon_layer = woody_carbon + soluable_carbon + polymeric_carbon + microbial_carbon + dom_carbon &
        &                 + residue_carbon + dom_assoc_carbon + residue_assoc_carbon
    END IF
  END SUBROUTINE calc_analytical_ss_jsm_layer_state

  ! ======================================================================================================= !
  !> calculate analytical steady-state carbon pools for the Jena Soil Model
  !> The current C-only analytical draft uses
  !>   C_wl*   = I_wl tau_wl
  !>   C_sol*  = I_sol tau_sol
  !>   C_D*    = (k_B * K_m,U) / (CUE * vmax_U - k_B)
  !> and the layer-wise DOC transport closure
  !>   C_B,sl*  = CUE_sl / (k_B * (1 - CUE_sl)) * C_eff_input_sl
  !>   C_eff_input_sl = (eta_DOM^sol * I_sol,sl + I_poly,sl + eta_poly^wl * I_wl,sl
  !>                - lambda_sl * C_D,sl* + lambda_(sl-1) * C_D,sl-1*)
  !> together with the associated-pool and sorption-capacity closure
  !>   C_aDOM* = ((k_ads,D [s-1] * C_D* [mol C m-3]) / k_des,D [s-1]) * Qmax* [mol C m-3]
  !>             / (1 + (k_ads,R / k_des,R) * C_res* + (k_ads,D / k_des,D) * C_D*)
  !>   C_aRes* = ((k_ads,R [s-1] * C_res* [mol C m-3]) / k_des,R [s-1]) * Qmax* [mol C m-3]
  !>             / (1 + (k_ads,R / k_des,R) * C_res* + (k_ads,D / k_des,D) * C_D*)
  !>   Qmax*   = (A * (1 - C0* / B)) / (1 + A * F / B)
  !>   C0*     = C_sol* + C_poly* + C_D* + C_B* + C_res*
  !>   A       = qmax_min [mol C kg-1] * rho_bulk [kg m-3]
  !>   B       = f_C^SOM [-] * rho_org [kg SOM m-3], converted below to [mol C m-3]
  !>   F       = (r_D * C_D* + r_R * C_res*) / (1 + r_D * C_D* + r_R * C_res*)
  !> before solving the polymeric-litter and residue pools from the enzyme-allocation closure.
  !>   C_poly*=(I_poly + eta_poly^wl * I_wl)/(vmax_poly * enzyme_frac_poly * C_B* )(Km_depoly + enzyme_frac_poly * C_B*)
  !>   C_res*=((1-f_soluable_residue) * k_B * C_B*)/(vmax_res * enzyme_frac_residue * C_B*)(Km_depoly + enzyme_frac_residue * C_B*)
  !>
  !> Unit convention of the analytical status file / sb_analytical_init_vars used here:
  !>   I_sol, I_poly, I_wl        [mol C m-3 s-1]
  !>   lambda_doc                 [s-1]
  !>   tau_*                      [s]
  !>   k_B, vmax_*, k_ads, k_des  [s-1]
  !>   K_m,*                      [mol C m-3]
  !>   qmax_org_fine_particle     [mol C kg-1]
  !>   rho_bulk                   [kg m-3]
  !>   rtm_*, rmm_*, enzyme_frac* [-]
  !>   C_*, Qmax*                 [mol C m-3]
  !> This means the analytical-init file must store mean rate fields, not timestep sums,
  !> for the litter-formation and DOC-transport terms.
  SUBROUTINE calc_analytical_ss_jsm( &
    & nproma, &
    & nsoil_w, &
    & nblks, &
    & num_sl_above_bedrock, &
    & qmax_org_fine_particle, &
    & soil_lay_depth_center_sl, &
    & soil_lay_depth_ubound_sl, &
    & silt_sl, &
    & clay_sl, &
    & bulk_dens_sl, &
    & bulk_carbon_max, &
    & qmax_org, &
    & rtm_depolymerisation, &
    & rmm_depolymerisation, &
    & rtm_mic_uptake, &
    & rmm_mic_uptake, &
    & rtm_sorption, &
    & rmm_sorption, &
    & rtm_desorption, &
    & rmm_desorption, &
    & rtm_hsc, &
    & rmm_hsc, &
    & sb_pool_mt_domain, &
    & surface_burial_frac, &
    & enzyme_frac_poly, &
    & enzyme_frac_residue, &
    & enzyme_frac_poly_c, &
    & enzyme_frac_poly_n, &
    & enzyme_frac_poly_p, &
    & enzyme_frac_AP)

    USE mo_sb_util,                       ONLY: calc_qmax_texture
    USE mo_sb_constants
    USE mo_jsb_physical_constants, ONLY: molar_mass_C

    INTEGER,  INTENT(in)    :: nproma
    INTEGER,  INTENT(in)    :: nsoil_w
    INTEGER,  INTENT(in)    :: nblks
    REAL(wp), INTENT(in)    :: num_sl_above_bedrock(:,:)
    REAL(wp), INTENT(in)    :: qmax_org_fine_particle(:,:,:)
    REAL(wp), INTENT(in)    :: soil_lay_depth_center_sl(:,:,:)
    REAL(wp), INTENT(in)    :: soil_lay_depth_ubound_sl(:,:,:)
    REAL(wp), INTENT(in)    :: silt_sl(:,:,:)
    REAL(wp), INTENT(in)    :: clay_sl(:,:,:)
    REAL(wp), INTENT(in)    :: bulk_dens_sl(:,:,:)
    REAL(wp), INTENT(inout) :: bulk_carbon_max(:,:,:)
    REAL(wp), INTENT(inout) :: qmax_org(:,:,:)
    REAL(wp), INTENT(in)    :: rtm_depolymerisation(:,:,:)
    REAL(wp), INTENT(in)    :: rmm_depolymerisation(:,:,:)
    REAL(wp), INTENT(in)    :: rtm_mic_uptake(:,:,:)
    REAL(wp), INTENT(in)    :: rmm_mic_uptake(:,:,:)
    REAL(wp), INTENT(in)    :: rtm_sorption(:,:,:)
    REAL(wp), INTENT(in)    :: rmm_sorption(:,:,:)
    REAL(wp), INTENT(in)    :: rtm_desorption(:,:,:)
    REAL(wp), INTENT(in)    :: rmm_desorption(:,:,:)
    REAL(wp), INTENT(in)    :: rtm_hsc(:,:,:)
    REAL(wp), INTENT(in)    :: rmm_hsc(:,:,:)
    REAL(wp), INTENT(inout) :: sb_pool_mt_domain(:,:,:,:,:)
    REAL(wp), INTENT(inout) :: surface_burial_frac(:,:,:)
    REAL(wp), INTENT(out)   :: enzyme_frac_poly(:,:,:)
    REAL(wp), INTENT(out)   :: enzyme_frac_residue(:,:,:)
    REAL(wp), INTENT(out)   :: enzyme_frac_poly_c(:,:,:)
    REAL(wp), INTENT(out)   :: enzyme_frac_poly_n(:,:,:)
    REAL(wp), INTENT(out)   :: enzyme_frac_poly_p(:,:,:)
    REAL(wp), INTENT(out)   :: enzyme_frac_AP(:,:,:)

    REAL(wp)                    :: input_sol
    REAL(wp)                    :: input_poly
    REAL(wp)                    :: input_wood
    REAL(wp)                    :: eff_input_c
    REAL(wp)                    :: alpha_poly
    REAL(wp)                    :: alpha_residue
    REAL(wp)                    :: denominator
    REAL(wp)                    :: km_uptake_eff
    REAL(wp)                    :: vmax_uptake_eff
    REAL(wp)                    :: km_depoly_eff
    REAL(wp)                    :: vmax_poly_eff
    REAL(wp)                    :: vmax_res_eff
    REAL(wp)                    :: qmax_org_min_local
    REAL(wp)                    :: qmax_b
    REAL(wp)                    :: lambda_doc
    REAL(wp)                    :: lambda_doc_above
    REAL(wp)                    :: rtm_dep_eff
    REAL(wp)                    :: rmm_dep_eff
    REAL(wp)                    :: rtm_upt_eff
    REAL(wp)                    :: rmm_upt_eff
    REAL(wp)                    :: rtm_sorb_eff
    REAL(wp)                    :: rmm_sorb_eff
    REAL(wp)                    :: rtm_desorb_eff
    REAL(wp)                    :: rmm_desorb_eff
    REAL(wp)                    :: rtm_hsc_eff
    REAL(wp)                    :: rmm_hsc_eff
    REAL(wp)                    :: bulk_dens_local
    REAL(wp)                    :: silt_local
    REAL(wp)                    :: clay_local
    REAL(wp)                    :: microbial_c_bio
    REAL(wp)                    :: polymeric_carbon
    REAL(wp)                    :: residue_carbon
    REAL(wp)                    :: dom_assoc_carbon
    REAL(wp)                    :: residue_assoc_carbon
    REAL(wp)                    :: bulk_carbon_bio
    REAL(wp)                    :: bulk_carbon_mid
    REAL(wp)                    :: min_volume_min_layer
    REAL(wp)                    :: cb_bracket_tol
    REAL(wp)                    :: surface_bulk_nonsorbed
    REAL(wp)                    :: surface_bulk_assoc
    REAL(wp)                    :: qmax_org_target
    REAL(wp)                    :: assoc_scale
    REAL(wp)                    :: enzyme_frac_total_scale
    REAL(wp)                    :: enzyme_scale_lower
    REAL(wp)                    :: enzyme_scale_upper
    REAL(wp)                    :: enzyme_scale_mid
    REAL(wp)                    :: bulk_carbon_scale_upper
    REAL(wp), PARAMETER         :: enzyme_frac_total_phys_scale = 10.0_wp
    LOGICAL                     :: l_surface_layer
    LOGICAL                     :: l_zero_deep_profile
    INTEGER                     :: iblk, is, ic
    INTEGER                     :: iter_phys_constraint
    INTEGER, PARAMETER          :: max_iter_phys_constraint = 200
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_analytical_ss_jsm'
    !> init out var
    !> The current analytical steady-state solution only solves for the carbon pools, so we can safely
    !> initialise the other elements with zeros here. The actual NetCDF reader will need to fill
    sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_microbial, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_dom, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_residue, ixC, :, :, :) = 0.0_wp
    sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) = 0.0_wp
    qmax_org(:,:,:) = 0.0_wp
    surface_burial_frac(:,:,:) = 1.0_wp
    enzyme_frac_poly(:,:,:) = 0.0_wp
    enzyme_frac_residue(:,:,:) = 0.0_wp
    enzyme_frac_poly_c(:,:,:) = 0.0_wp
    enzyme_frac_poly_n(:,:,:) = 0.0_wp
    enzyme_frac_poly_p(:,:,:) = 0.0_wp
    enzyme_frac_AP(:,:,:) = 0.0_wp

    DO iblk = 1,nblks
      DO ic = 1,nproma
        l_zero_deep_profile = .FALSE.
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))

          input_sol  = MAX(0.0_wp, sb_analytical_init_vars%formation_soluable_litter_c(ic, is, iblk))
          input_poly = MAX(0.0_wp, sb_analytical_init_vars%formation_polymeric_litter_c(ic, is, iblk))
          input_wood = MAX(0.0_wp, sb_analytical_init_vars%formation_woody_litter_c(ic, is, iblk))

          !> Prefer the modifier and soil-property fields carried by the analytical status file.
          !> Until the dedicated NetCDF reader is implemented, these temporary arrays hold neutral
          !> defaults and we safely fall back to the currently available runtime fields.
          bulk_dens_local = MERGE(sb_analytical_init_vars%bulk_dens_sl(ic, is, iblk), &
            &                     bulk_dens_sl(ic, is, iblk), sb_analytical_init_vars%bulk_dens_sl(ic, is, iblk) > eps8)
          silt_local = MERGE(sb_analytical_init_vars%silt_sl(ic, is, iblk), &
            &                silt_sl(ic, is, iblk), sb_analytical_init_vars%silt_sl(ic, is, iblk) > eps8)
          clay_local = MERGE(sb_analytical_init_vars%clay_sl(ic, is, iblk), &
            &                clay_sl(ic, is, iblk), sb_analytical_init_vars%clay_sl(ic, is, iblk) > eps8)

          rtm_dep_eff    = MERGE(sb_analytical_init_vars%rtm_depolymerisation(ic, is, iblk), &
            &                    rtm_depolymerisation(ic, is, iblk), &
            &                    sb_analytical_init_vars%rtm_depolymerisation(ic, is, iblk) > eps8)
          rmm_dep_eff    = MERGE(sb_analytical_init_vars%rmm_depolymerisation(ic, is, iblk), &
            &                    rmm_depolymerisation(ic, is, iblk), &
            &                    sb_analytical_init_vars%rmm_depolymerisation(ic, is, iblk) > eps8)
          rtm_upt_eff    = MERGE(sb_analytical_init_vars%rtm_mic_uptake(ic, is, iblk), &
            &                    rtm_mic_uptake(ic, is, iblk), sb_analytical_init_vars%rtm_mic_uptake(ic, is, iblk) > eps8)
          rmm_upt_eff    = MERGE(sb_analytical_init_vars%rmm_mic_uptake(ic, is, iblk), &
            &                    rmm_mic_uptake(ic, is, iblk), sb_analytical_init_vars%rmm_mic_uptake(ic, is, iblk) > eps8)
          rtm_sorb_eff   = MERGE(sb_analytical_init_vars%rtm_sorption(ic, is, iblk), &
            &                    rtm_sorption(ic, is, iblk), sb_analytical_init_vars%rtm_sorption(ic, is, iblk) > eps8)
          rmm_sorb_eff   = MERGE(sb_analytical_init_vars%rmm_sorption(ic, is, iblk), &
            &                    rmm_sorption(ic, is, iblk), sb_analytical_init_vars%rmm_sorption(ic, is, iblk) > eps8)
          rtm_desorb_eff = MERGE(sb_analytical_init_vars%rtm_desorption(ic, is, iblk), &
            &                    rtm_desorption(ic, is, iblk), sb_analytical_init_vars%rtm_desorption(ic, is, iblk) > eps8)
          rmm_desorb_eff = MERGE(sb_analytical_init_vars%rmm_desorption(ic, is, iblk), &
            &                    rmm_desorption(ic, is, iblk), sb_analytical_init_vars%rmm_desorption(ic, is, iblk) > eps8)
          rtm_hsc_eff    = MERGE(sb_analytical_init_vars%rtm_hsc(ic, is, iblk), &
            &                    rtm_hsc(ic, is, iblk), sb_analytical_init_vars%rtm_hsc(ic, is, iblk) > eps8)
          rmm_hsc_eff    = MERGE(sb_analytical_init_vars%rmm_hsc(ic, is, iblk), &
            &                    rmm_hsc(ic, is, iblk), sb_analytical_init_vars%rmm_hsc(ic, is, iblk) > eps8)

          rtm_dep_eff    = MERGE(rtm_dep_eff,    1.0_wp, rtm_dep_eff    > eps8)
          rmm_dep_eff    = MERGE(rmm_dep_eff,    1.0_wp, rmm_dep_eff    > eps8)
          rtm_upt_eff    = MERGE(rtm_upt_eff,    1.0_wp, rtm_upt_eff    > eps8)
          rmm_upt_eff    = MERGE(rmm_upt_eff,    1.0_wp, rmm_upt_eff    > eps8)
          rtm_sorb_eff   = MERGE(rtm_sorb_eff,   1.0_wp, rtm_sorb_eff   > eps8)
          rmm_sorb_eff   = MERGE(rmm_sorb_eff,   1.0_wp, rmm_sorb_eff   > eps8)
          rtm_desorb_eff = MERGE(rtm_desorb_eff, 1.0_wp, rtm_desorb_eff > eps8)
          rmm_desorb_eff = MERGE(rmm_desorb_eff, 1.0_wp, rmm_desorb_eff > eps8)
          rtm_hsc_eff    = MERGE(rtm_hsc_eff,    1.0_wp, rtm_hsc_eff    > eps8)
          rmm_hsc_eff    = MERGE(rmm_hsc_eff,    1.0_wp, rmm_hsc_eff    > eps8)

          km_uptake_eff = km_mic_uptake * rtm_hsc_eff * rmm_hsc_eff
          vmax_uptake_eff = vmax_mic_uptake * rtm_upt_eff * rmm_upt_eff
          km_depoly_eff = km_depolymerisation * rtm_hsc_eff * rmm_hsc_eff
          vmax_poly_eff = vmax_depolymerisation_litter * rtm_dep_eff * rmm_dep_eff
          vmax_res_eff = vmax_depolymerisation_residue * rtm_dep_eff * rmm_dep_eff
          l_surface_layer = (is == 1)

          !> C_wl*   = I_wl [mol C m-3 s-1] * tau_wl [s]
          sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk)    = input_wood * tau_woody_litter
          !> C_sol*  = I_sol [mol C m-3 s-1] * tau_sol [s]
          sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) = input_sol * tau_soluable_litter

          denominator = microbial_cue_max * vmax_uptake_eff - k_mort_microbial
          IF (denominator <= eps8) THEN
            CALL finish(TRIM(routine), &
              & 'Analytical JSM initialisation requires modified microbial uptake to exceed microbial mortality.')
          END IF
          !> C_D*    = (k_B [s-1] * K_m,U [mol C m-3]) / (CUE * vmax_U [s-1] - k_B [s-1])
          sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) = k_mort_microbial * km_uptake_eff / denominator

          lambda_doc = MAX(0.0_wp, sb_analytical_init_vars%frac_wtr_transp_down_sl(ic, is, iblk))
          IF (is > 1) THEN
            lambda_doc_above = MAX(0.0_wp, sb_analytical_init_vars%frac_wtr_transp_down_sl(ic, is - 1, iblk))
          ELSE
            lambda_doc_above = 0.0_wp
          END IF

          eff_input_c = fc_soluable2dom * input_sol + input_poly + fc_woody2polymeric * input_wood &
            &            - lambda_doc * sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)
          IF (is > 1) THEN
            eff_input_c = eff_input_c + lambda_doc_above * sb_pool_mt_domain(ix_dom, ixC, ic, is - 1, iblk)
          END IF

          !> C_B,sl* = CUE [-] / (k_B [s-1] * (1 - CUE)) * C_eff_input_sl [mol C m-3 s-1]
          IF (l_zero_deep_profile) THEN
            microbial_c_bio = 0.0_wp
          ELSE IF (eff_input_c > eps8) THEN
            microbial_c_bio = microbial_cue_max / &
              & (k_mort_microbial * (1._wp - microbial_cue_max)) * eff_input_c
          ELSE
            microbial_c_bio = 0.0_wp
            IF (is > 1) THEN
              l_zero_deep_profile = .TRUE.
            END IF
          END IF

          enzyme_frac_AP(ic, is, iblk) = MIN(MAX(sb_analytical_init_vars%enzyme_frac_ap(ic, is, iblk), &
            & min_enzyme_fraction), max_enzyme_fraction)
          qmax_org_min_local = calc_qmax_texture(qmax_org_fine_particle(ic, is, iblk), silt_local, clay_local)
          qmax_org_min_local = qmax_org_min_local * EXP(-is / (INT(num_sl_above_bedrock(ic, iblk)) * k_sl_conv_moc_corr)) * 0.5_wp
          IF (soil_lay_depth_ubound_sl(ic, is, iblk) > z_deep_soil_moc_corr) THEN
            qmax_org_min_local = qmax_org_min_local * MAX(0.0001_wp, 10._wp**(-soil_lay_depth_ubound_sl(ic, is, iblk)))
          END IF

          !> Steady-state-only physical correction:
          !> phi_min^min(1) = 0.05
          !> phi_min^min(k) = phi_min^min(k-1) + 0.01 * (z_mid,k - z_mid,k-1)
          !> C_bulk^max(k) = B * (1 - phi_min^min(k))
          !> For the first mineral layer, if C_bulk^bio > C_bulk^max, solve the geometric
          !> surface-soil partition directly with
          !>   f_burial = (C_bulk^max - (C_aDOM* + C_aRes*))
          !>              / (C_sol* + C_poly* + C_B* + C_D* + C_res*)
          !> so only the buried fraction of the non-woody, non-associated surface SOM occupies soil volume.
          !> For deeper layers, keep the current bisection search on C_B.
          qmax_b = carbon_per_dryweight_SOM * rho_bulk_org * 1000._wp / molar_mass_C
          bulk_carbon_max(ic, is, iblk) = MAX(0.0_wp, bulk_carbon_max(ic, is, iblk))
          min_volume_min_layer = MAX(0.0_wp, 1.0_wp - bulk_carbon_max(ic, is, iblk) / qmax_b)
          enzyme_frac_total_scale = 1.0_wp

          CALL calc_analytical_ss_jsm_layer_state( &
            & l_surface_layer, &
            & input_sol, &
            & input_poly, &
            & input_wood, &
            & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
            & microbial_c_bio, &
            & enzyme_frac_AP(ic, is, iblk), &
            & enzyme_frac_total_scale, &
            & km_depoly_eff, &
            & vmax_poly_eff, &
            & vmax_res_eff, &
            & qmax_org_min_local, &
            & bulk_dens_local, &
            & rtm_sorb_eff, &
            & rmm_sorb_eff, &
            & rtm_desorb_eff, &
            & rmm_desorb_eff, &
            & alpha_poly, &
            & alpha_residue, &
            & polymeric_carbon, &
            & residue_carbon, &
            & dom_assoc_carbon, &
            & residue_assoc_carbon, &
            & qmax_org(ic, is, iblk), &
            & bulk_carbon_bio)

          surface_burial_frac(ic, is, iblk) = 1.0_wp
          IF (l_surface_layer .AND. bulk_carbon_bio > bulk_carbon_max(ic, is, iblk)) THEN
            !> For the capped surface layer, qmax_org is reset by the prescribed minimum
            !> mineral-soil volume rather than by the unconstrained analytical C_B solution:
            !>   Qmax_org = phi_min^min * qmax_org_min * rho_bulk
            qmax_org_target = min_volume_min_layer * qmax_org_min_local * bulk_dens_local
            IF (qmax_org(ic, is, iblk) > eps8) THEN
              assoc_scale = qmax_org_target / qmax_org(ic, is, iblk)
            ELSE
              assoc_scale = 0.0_wp
            END IF
            qmax_org(ic, is, iblk) = qmax_org_target
            dom_assoc_carbon       = dom_assoc_carbon * assoc_scale
            residue_assoc_carbon   = residue_assoc_carbon * assoc_scale
            !> C_bulk,1^soil = f_burial * (C_sol* + C_poly* + C_B* + C_D* + C_res*) + C_aDOM* + C_aRes*
            !> f_burial = (C_bulk,max - (C_aDOM* + C_aRes*))
            !>              / (C_sol* + C_poly* + C_B* + C_D* + C_res*)
            surface_bulk_nonsorbed = input_sol * tau_soluable_litter + polymeric_carbon + microbial_c_bio &
              &                     + sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) + residue_carbon
            surface_bulk_assoc     = dom_assoc_carbon + residue_assoc_carbon
            IF (surface_bulk_nonsorbed > eps8) THEN
              surface_burial_frac(ic, is, iblk) = (bulk_carbon_max(ic, is, iblk) - surface_bulk_assoc) / surface_bulk_nonsorbed
              surface_burial_frac(ic, is, iblk) = MIN(1.0_wp, MAX(0.0_wp, surface_burial_frac(ic, is, iblk)))
            ELSE
              surface_burial_frac(ic, is, iblk) = 0.0_wp
            END IF
          ELSE IF (bulk_carbon_bio > bulk_carbon_max(ic, is, iblk) .AND. microbial_c_bio > min_mic_biomass) THEN
            enzyme_scale_lower = 1.0_wp
            enzyme_scale_upper = MAX(enzyme_frac_total_phys_scale, bulk_carbon_bio / bulk_carbon_max(ic, is, iblk))
            enzyme_frac_total_scale = enzyme_scale_upper
            CALL calc_analytical_ss_jsm_layer_state( &
              & l_surface_layer, &
              & input_sol, &
              & input_poly, &
              & input_wood, &
              & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
              & microbial_c_bio, &
              & enzyme_frac_AP(ic, is, iblk), &
              & enzyme_frac_total_scale, &
              & km_depoly_eff, &
              & vmax_poly_eff, &
              & vmax_res_eff, &
              & qmax_org_min_local, &
              & bulk_dens_local, &
              & rtm_sorb_eff, &
              & rmm_sorb_eff, &
              & rtm_desorb_eff, &
              & rmm_desorb_eff, &
              & alpha_poly, &
              & alpha_residue, &
              & polymeric_carbon, &
              & residue_carbon, &
              & dom_assoc_carbon, &
              & residue_assoc_carbon, &
              & qmax_org(ic, is, iblk), &
              & bulk_carbon_scale_upper)
            IF (bulk_carbon_scale_upper <= bulk_carbon_max(ic, is, iblk)) THEN
              cb_bracket_tol = eps4 * MAX(1.0_wp, enzyme_scale_upper)
              !> Bisection solve for enzyme_frac_total_scale such that
              !> C_bulk(C_B^bio, enzyme_frac_total_scale) = C_bulk^max
              DO iter_phys_constraint = 1, max_iter_phys_constraint
                enzyme_scale_mid = 0.5_wp * (enzyme_scale_lower + enzyme_scale_upper)
                CALL calc_analytical_ss_jsm_layer_state( &
                  & l_surface_layer, &
                  & input_sol, &
                  & input_poly, &
                  & input_wood, &
                  & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
                  & microbial_c_bio, &
                  & enzyme_frac_AP(ic, is, iblk), &
                  & enzyme_scale_mid, &
                  & km_depoly_eff, &
                  & vmax_poly_eff, &
                  & vmax_res_eff, &
                  & qmax_org_min_local, &
                  & bulk_dens_local, &
                  & rtm_sorb_eff, &
                  & rmm_sorb_eff, &
                  & rtm_desorb_eff, &
                  & rmm_desorb_eff, &
                  & alpha_poly, &
                  & alpha_residue, &
                  & polymeric_carbon, &
                  & residue_carbon, &
                  & dom_assoc_carbon, &
                  & residue_assoc_carbon, &
                  & qmax_org(ic, is, iblk), &
                  & bulk_carbon_mid)
                IF (bulk_carbon_mid > bulk_carbon_max(ic, is, iblk)) THEN
                  enzyme_scale_lower = enzyme_scale_mid
                ELSE
                  enzyme_scale_upper = enzyme_scale_mid
                  bulk_carbon_scale_upper = bulk_carbon_mid
                  IF ((bulk_carbon_max(ic, is, iblk) - bulk_carbon_mid) <= 0.01_wp * bulk_carbon_max(ic, is, iblk)) EXIT
                END IF
                IF ((enzyme_scale_upper - enzyme_scale_lower) <= cb_bracket_tol) EXIT
              END DO
              enzyme_frac_total_scale = enzyme_scale_upper
            ELSE IF (bulk_carbon_bio > bulk_carbon_max(ic, is, iblk)) THEN
              enzyme_frac_total_scale = enzyme_scale_upper
            END IF

            CALL calc_analytical_ss_jsm_layer_state( &
              & l_surface_layer, &
              & input_sol, &
              & input_poly, &
              & input_wood, &
              & sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk), &
              & microbial_c_bio, &
              & enzyme_frac_AP(ic, is, iblk), &
              & enzyme_frac_total_scale, &
              & km_depoly_eff, &
              & vmax_poly_eff, &
              & vmax_res_eff, &
              & qmax_org_min_local, &
              & bulk_dens_local, &
              & rtm_sorb_eff, &
              & rmm_sorb_eff, &
              & rtm_desorb_eff, &
              & rmm_desorb_eff, &
              & alpha_poly, &
              & alpha_residue, &
              & polymeric_carbon, &
              & residue_carbon, &
              & dom_assoc_carbon, &
              & residue_assoc_carbon, &
              & qmax_org(ic, is, iblk), &
              & bulk_carbon_bio)
          END IF

          sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) = microbial_c_bio
          enzyme_frac_poly(ic, is, iblk)    = alpha_poly
          enzyme_frac_residue(ic, is, iblk) = alpha_residue
          !> In the C-only draft, the carbon-based analytical allocation is the only closed solution.
          !> Keep the N and P diagnostic fractions aligned with the carbon solution until the
          !> nutrient-limited analytical system is added.
          enzyme_frac_poly_c(ic, is, iblk)  = alpha_poly
          enzyme_frac_poly_n(ic, is, iblk)  = alpha_poly
          enzyme_frac_poly_p(ic, is, iblk)  = alpha_poly
          sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = polymeric_carbon
          sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk)          = residue_carbon
          sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)        = dom_assoc_carbon
          sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)    = residue_assoc_carbon
        END DO
      END DO
    END DO

  END SUBROUTINE calc_analytical_ss_jsm

  ! ======================================================================================================= !
  !>init SB_ with the simple soil model
  !>
  SUBROUTINE sb_init_simple_sm( &
    & nproma, &
    & nsoil_w, &
    & nblks, &
    & num_sl_above_bedrock, &
    & lctlib_growthform, &
    & lctlib_k_som_fast_init, &
    & lctlib_k_som_slow_init, &
    & sb_model_scheme, &
    & sb_init_scheme, &
    & flag_sb_double_langmuir, &
    & qmax_org_fine_particle, &
    & soil_p_labile, &
    & soil_p_slow, &
    & soil_p_occluded, &
    & soil_p_primary, &
    & soil_p_depth, &
    & soil_depth_sl, &
    & soil_lay_width_sl, &
    & soil_lay_depth_ubound_sl, &
    & soil_lay_depth_center_sl, &
    & hlp_frac_c13, &
    & hlp_frac_c14, &
    & root_fraction_sl, &
    & clay_sl, &
    & silt_sl, &
    & vol_porosity_sl, &
    & bulk_dens_sl, &
    & ph_sl, &
    & bulk_carbon_max, &
    & surface_burial_frac, &
    & bulk_soil_carbon_sl, &
    & soil_litter_carbon_sl, &
    & po4_solute, &
    & volume_min_sl, &
    & qmax_po4, &
    & qmax_po4_min_sl, &
    & qmax_po4_om_sl, &
    & qmax_fast_po4, &
    & qmax_slow_po4, &
    & km_fast_po4, &
    & km_slow_po4, &
    & km_adsorpt_po4_sl, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & po4_primary, &
    & Qmax_AlFe_cor, &
    & sb_pool_mt_domain, &
    & nh4_solute, &
    & no3_solute, &
    & nh4_n15_solute, &
    & no3_n15_solute, &
    & bulk_dens_corr_sl, &
    & qmax_org, &
    & qmax_org_min_sl, &
    & qmax_nh4, &
    & qmax_nh4_min_sl, &
    & km_adsorpt_nh4_sl, &
    & nh4_assoc, &
    & nh4_n15_assoc)

    USE mo_q_sb_jsm_processes,            ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, &
      &                                         calc_qmax_bulk_density_correction, calc_Psorption_parameter, &
      &                                         calc_fast_po4
    USE mo_sb_util,                       ONLY: calc_qmax_texture
    USE mo_isotope_util,                  ONLY: calc_mixing_ratio_N15N14
    USE mo_veg_constants,                 ONLY: ITREE
    USE mo_sb_constants,                  ONLY: def_sb_pool_metabol_litter, def_sb_pool_simple_struct_litter, &
      &                                         def_sb_pool_woody_litter, def_sb_pool_fast_som, def_sb_pool_slow_som, &
      &                                         nh4_solute_prescribe, no3_solute_prescribe, po4_solute_prescribe, &
      &                                         p_labile_slow_global_avg, frac_prod_not_root, qmax_nh4_clay, &
      &                                         km_adsorpt_OM_nh4, km_adsorpt_mineral_nh4, &
      &                                         k_init_soluable_litter_cn, k_init_litter_np, k_init_polymeric_litter_cn, &
      &                                         k_init_woody_litter_cn, k_fast_som_np, k_fast_som_cn_min, k_fast_som_cn_max, &
      &                                         k_slow_som_cn, k_slow_som_np, &
      &                                         reference_depth_som_init, k_som_init, &
      &                                         sb_pool_total_som_at_ref_depth_init, fract_som_fast
    USE mo_jsb_physical_constants,        ONLY: molar_mass_C
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                    INTENT(in)    :: nproma
    INTEGER,                    INTENT(in)    :: nsoil_w
    INTEGER,                    INTENT(in)    :: nblks
    REAL(wp),                   INTENT(in)    :: num_sl_above_bedrock(:,:)
    INTEGER,                    INTENT(in)    :: lctlib_growthform
    REAL(wp),                   INTENT(in)    :: lctlib_k_som_fast_init
    REAL(wp),                   INTENT(in)    :: lctlib_k_som_slow_init
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme                  !< SB_ config
    CHARACTER(len=*),           INTENT(in)    :: sb_init_scheme                   !< SB_ config
    LOGICAL,                    INTENT(in)    :: flag_sb_double_langmuir          !< SB_ config
    REAL(wp),                   INTENT(in)    :: qmax_org_fine_particle(:,:,:)    !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_labile(:,:,:)             !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_slow(:,:,:)               !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_occluded(:,:,:)           !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_primary(:,:,:)            !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_depth                     !< SB_ config
    REAL(wp),                   INTENT(in)    :: soil_depth_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_width_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_ubound_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_center_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: hlp_frac_c13
    REAL(wp),                   INTENT(in)    :: hlp_frac_c14
    REAL(wp),                   INTENT(in)    :: root_fraction_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: clay_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: silt_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: vol_porosity_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: bulk_dens_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: ph_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: bulk_carbon_max(:,:,:)
    REAL(wp),                   INTENT(inout) :: surface_burial_frac(:,:,:)
    REAL(wp),                   INTENT(inout) :: bulk_soil_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: soil_litter_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_solute(:,:,:)
    REAL(wp),                   INTENT(inout) :: volume_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_om_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_adsorpt_po4_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_fast(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_slow(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_occluded(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_primary(:,:,:)
    REAL(wp),                   INTENT(inout) :: Qmax_AlFe_cor(:,:,:)
    REAL(wp),                   INTENT(inout) :: sb_pool_mt_domain(:,:,:,:,:)          !< bgc_material sb_pool
    REAL(wp),                   INTENT(out)   :: nh4_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: bulk_dens_corr_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_org(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_org_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_nh4(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_nh4_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: km_adsorpt_nh4_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_assoc(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_assoc(:,:,:)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), ALLOCATABLE       :: hlp1(:,:)
    REAL(wp), ALLOCATABLE       :: hlp2(:,:)
    REAL(wp), ALLOCATABLE       :: weight_fast_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: weight_slow_sl(:,:,:)
    REAL(wp), ALLOCATABLE       :: init_soluable_litter(:,:)
    REAL(wp), ALLOCATABLE       :: init_polymeric_litter(:,:)
    REAL(wp), ALLOCATABLE       :: init_fast(:,:)
    REAL(wp), ALLOCATABLE       :: init_slow(:,:)
    REAL(wp), ALLOCATABLE       :: p_labile(:,:,:)
    REAL(wp), ALLOCATABLE       :: p_labile_slow_site(:,:,:)
    INTEGER                     :: iblk, is, ic
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init_simple_sm'
    IF (debug_on()) CALL message(TRIM(routine), 'Starting ...')
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    nh4_solute(:,:,:)         = 0.0_wp
    no3_solute(:,:,:)         = 0.0_wp
    nh4_n15_solute(:,:,:)     = 0.0_wp
    no3_n15_solute(:,:,:)     = 0.0_wp
    bulk_dens_corr_sl(:,:,:)  = 0.0_wp
    qmax_org(:,:,:)           = 0.0_wp
    qmax_org_min_sl(:,:,:)    = 0.0_wp
    surface_burial_frac(:,:,:) = 1.0_wp
    qmax_nh4(:,:,:)           = 0.0_wp
    qmax_nh4_min_sl(:,:,:)    = 0.0_wp
    km_adsorpt_nh4_sl(:,:,:)  = 0.0_wp
    nh4_assoc(:,:,:)          = 0.0_wp
    nh4_n15_assoc(:,:,:)      = 0.0_wp

    !>0.9 allocate local var
    !>
    ALLOCATE(hlp1(nproma, nblks))
    ALLOCATE(hlp2(nproma, nblks))
    ALLOCATE(weight_fast_sl(nproma, nsoil_w, nblks))
    ALLOCATE(weight_slow_sl(nproma, nsoil_w, nblks))
    ALLOCATE(init_soluable_litter(nproma, nblks))
    ALLOCATE(init_polymeric_litter(nproma, nblks))
    ALLOCATE(init_fast(nproma, nblks))
    ALLOCATE(init_slow(nproma, nblks))
    ALLOCATE(p_labile(nproma, nsoil_w, nblks))
    ALLOCATE(p_labile_slow_site(nproma, nsoil_w, nblks))

    SELECT CASE(TRIM(sb_init_scheme))
    CASE('legacy')
    !>1.0 init
    !>
    hlp1(:,:)             = 0.0_wp
    hlp2(:,:)             = 0.0_wp
    weight_fast_sl(:,:,:) = 0.0_wp
    weight_slow_sl(:,:,:) = 0.0_wp
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          weight_fast_sl(ic, is, iblk) = EXP(-(lctlib_k_som_fast_init * soil_lay_depth_center_sl(ic, is, iblk)) ** k_som_init)
          weight_slow_sl(ic, is, iblk) = EXP(-(lctlib_k_som_slow_init * soil_lay_depth_center_sl(ic, is, iblk)) ** k_som_init)
          ! only count up to reference depth for scaling
          IF (soil_lay_depth_center_sl(ic, is, iblk) < reference_depth_som_init) THEN
            hlp1(ic, iblk) = hlp1(ic, iblk) + weight_fast_sl(ic, is, iblk) * soil_lay_width_sl(ic, is, iblk)
            hlp2(ic, iblk) = hlp2(ic, iblk) + weight_slow_sl(ic, is, iblk) * soil_lay_width_sl(ic, is, iblk)
          END IF
        END DO
      END DO
    END DO

    ! calculate top-soil SOM concentration and convert from default gC/m2 to mol C / m3
    WHERE(hlp1 > eps8)
      init_soluable_litter(:,:)   = def_sb_pool_metabol_litter * fract_som_fast / molar_mass_C / hlp1(:,:)
      init_polymeric_litter(:,:)  = def_sb_pool_simple_struct_litter * fract_som_fast / molar_mass_C / hlp1(:,:)
      init_fast(:,:)              = sb_pool_total_som_at_ref_depth_init * fract_som_fast / molar_mass_C / hlp1(:,:)
    END WHERE
    WHERE(hlp2 > eps8)
      init_slow(:,:) = sb_pool_total_som_at_ref_depth_init * (1._wp - fract_som_fast) / molar_mass_C / hlp2(:,:)
    END WHERE

    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)  = init_soluable_litter(ic, iblk) &
            &                                                         * weight_fast_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = init_polymeric_litter(ic, iblk) &
            &                                                         * weight_fast_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        = init_fast(ic, iblk) * weight_fast_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk)          = init_slow(ic, iblk) * weight_slow_sl(ic, is, iblk)
        END DO
      END DO
    END DO
    ! re-calculate the first layer of the wood litter pool for trees only
    IF (lctlib_growthform == ITREE) THEN
      sb_pool_mt_domain(ix_woody_litter, ixC, :, 1, :) = &
        &   sb_pool_total_som_at_ref_depth_init / molar_mass_C / soil_lay_width_sl(:,1,:)
    END IF
    CASE('analytical_ss')
      CALL calc_analytical_ss_simple_sm( &
        & nproma, &
        & nsoil_w, &
        & nblks, &
        & num_sl_above_bedrock, &
        & sb_pool_mt_domain)
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid setting of sb_init_scheme (legacy / analytical_ss).')
    END SELECT

    ! NOTE: this differs on purpose from the calculation in mo_q_sb_update_pools for the prescribed option
    !       this is because here I'm trying to ensure that the models are initialised with the same amount of nutrient
    !       whereas in the prescribed option the target is unlimiting concentrations across the profile. The latter leads
    !       to more nutrient in the soil column in 1d compared to 0d.
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          nh4_solute(ic, is, iblk) = nh4_solute_prescribe &
            &                           / soil_lay_depth_ubound_sl(ic, NINT(num_sl_above_bedrock(ic, iblk)), iblk)
          no3_solute(ic, is, iblk) = no3_solute_prescribe &
            &                           / soil_lay_depth_ubound_sl(ic, NINT(num_sl_above_bedrock(ic, iblk)), iblk)
        END DO
      END DO
    END DO

    ! inital bulk density and qmax_po4 based on current soil carbon content
    DO iblk = 1,nblks
      CALL calc_bulk_soil_carbon( &
        & nproma, &                           ! in
        & nsoil_w, &
        & num_sl_above_bedrock(:,iblk), &
        & TRIM(sb_model_scheme), &
        & sb_pool_mt_domain(:,:,:,:,iblk), &  ! in
        & bulk_carbon_max(:,:,iblk), &        ! in
        & surface_burial_frac(:,:,iblk), &    ! in
        & bulk_soil_carbon_sl(:,:,iblk), &    ! inout
        & soil_litter_carbon_sl(:,:,iblk), &
        & volume_min_sl(:,:,iblk) )           ! inout

      bulk_dens_corr_sl(:,:,iblk) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:,iblk), &
        &                                                        soil_litter_carbon_sl(:,:,iblk), &
        &                                                        bulk_dens_sl(:,:,iblk))
      qmax_org_min_sl(:,:,iblk) = calc_qmax_texture(qmax_org_fine_particle(:,:,iblk), &
        &                                           silt_sl(:,:,iblk), &
        &                                           clay_sl(:,:,iblk))

      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &  ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0._wp, &
        & qmax_org_min_sl(:,:,iblk), &      ! in
        & qmax_org(:,:,iblk) )              ! out
    END DO

    ! Correct the 1st soil layer Qmax_AlFe_cor(:,:,:) with the woody litter
    Qmax_AlFe_cor(:, 1, :) = Qmax_AlFe_cor(:, 1, :) * 1._wp / ((1._wp - volume_min_sl(:, 1, :)) &
      &                      * (bulk_soil_carbon_sl(:, 1, :) + soil_litter_carbon_sl(:, 1, :)  &
      &                      + sb_pool_mt_domain(ix_woody_litter, ixC, :,  1, :)) &
      &                      / (bulk_soil_carbon_sl(:, 1, :) + soil_litter_carbon_sl(:, 1, :)) &
      &                      + volume_min_sl(:, 1, :))
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          CALL calc_Psorption_parameter( &
            & flag_sb_double_langmuir, &          ! in
            & clay_sl(ic, is, iblk), &
            & silt_sl(ic, is, iblk), &
            & vol_porosity_sl(ic, is, iblk) * 1000._wp, &
            & ph_sl(ic, is, iblk), &
            & bulk_soil_carbon_sl(ic, is, iblk), &
            & soil_litter_carbon_sl(ic, is, iblk), &
            & sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk), &
            & bulk_dens_sl(ic, is, iblk), &
            & volume_min_sl(ic, is, iblk), &
            & po4_solute_prescribe, &             ! in - po4_solute(ic, is, iblk)
            & qmax_po4_min_sl(ic, is, iblk), &        ! inout
            & qmax_po4_om_sl(ic, is, iblk), &
            & qmax_po4(ic, is, iblk), &
            & km_adsorpt_po4_sl(ic, is, iblk), &
            & qmax_fast_po4(ic, is, iblk), &
            & qmax_slow_po4(ic, is, iblk), &
            & km_fast_po4(ic, is, iblk), &
            & km_slow_po4(ic, is, iblk), &
            & Qmax_AlFe_cor(ic, is, iblk) )           ! inout
        END DO
      END DO
    END DO
    ! calculate the partition coefficient between labile and slow P pool based on
    ! their Qmax (double-Langmuir) or global average values (single-Langmuir)
    IF (flag_sb_double_langmuir) THEN
      p_labile_slow_site(:,:,:) = eps4
      ! if 'qmax_po4(:,:,:) > eps8'
      DO iblk = 1,nblks
        DO is = 1,nsoil_w
          DO ic = 1,nproma
            IF (qmax_po4(ic, is, iblk) > eps8) THEN
              p_labile_slow_site(ic, is, iblk) = MAX(qmax_fast_po4(ic, is, iblk) / qmax_po4(ic, is, iblk), eps4)
            END IF
          END DO
        END DO
      END DO
    ELSE
      p_labile_slow_site(:,:,:) = p_labile_slow_global_avg
    END IF

    ! nh4 adsorption parameter intialisation
    qmax_nh4_min_sl(:,:,:) = calc_qmax_texture(qmax_nh4_clay, &
      &                                        0.0_wp, &
      &                                        clay_sl(:,:,:))
    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &  ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0.0_wp, &
        & qmax_nh4_min_sl(:,:,iblk), &      ! in
        & qmax_nh4(:,:,iblk) )              ! out
    END DO
    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &  ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & km_adsorpt_OM_nh4, &
        & km_adsorpt_mineral_nh4, &         ! in
        & km_adsorpt_nh4_sl(:,:,iblk) )     ! out
    END DO
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          nh4_assoc(ic, is, iblk)     = qmax_nh4(ic, is, iblk) * nh4_solute(ic, is, iblk) &
              &                           / (nh4_solute(ic, is, iblk) + km_adsorpt_nh4_sl(ic, is, iblk))
          nh4_n15_assoc(ic, is, iblk) = nh4_assoc(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! metabolic litter all other elements
          sb_pool_mt_domain(ix_soluable_litter, ixN, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) / k_init_soluable_litter_cn
          sb_pool_mt_domain(ix_soluable_litter, ixP, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixN, ic, is, iblk) / k_init_litter_np
          sb_pool_mt_domain(ix_soluable_litter, ixC13, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) * hlp_frac_c13
          sb_pool_mt_domain(ix_soluable_litter, ixC14, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) * hlp_frac_c14
          sb_pool_mt_domain(ix_soluable_litter, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_soluable_litter, ixN, ic, is, iblk) &
            &                                                          / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! non-woody litter all other elements
          sb_pool_mt_domain(ix_polymeric_litter, ixN, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) / k_init_polymeric_litter_cn
          sb_pool_mt_domain(ix_polymeric_litter, ixP, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_polymeric_litter, ixN, ic, is, iblk) / k_init_litter_np
          sb_pool_mt_domain(ix_polymeric_litter, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
            &                                                           * hlp_frac_c13
          sb_pool_mt_domain(ix_polymeric_litter, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
            &                                                           * hlp_frac_c14
          sb_pool_mt_domain(ix_polymeric_litter, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_polymeric_litter, ixN, ic, is, iblk) &
            &                                                           / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! woody litter all other elements
          sb_pool_mt_domain(ix_woody_litter, ixN, ic, is, iblk) &
            &  = sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk) / k_init_woody_litter_cn
          sb_pool_mt_domain(ix_woody_litter, ixP, ic, is, iblk)   = sb_pool_mt_domain(ix_woody_litter, ixN, ic, is, iblk) &
            &                                                       / k_init_litter_np
          sb_pool_mt_domain(ix_woody_litter, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk) &
            &                                                       * hlp_frac_c13
          sb_pool_mt_domain(ix_woody_litter, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk) &
            &                                                       * hlp_frac_c14
          sb_pool_mt_domain(ix_woody_litter, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_woody_litter, ixN, ic, is, iblk) &
            &                                                       / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! fast pool (= microbial) all other elements
          sb_pool_mt_domain(ix_microbial, ixN, ic, is, iblk) &
            & = sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) / (k_fast_som_cn_min + k_fast_som_cn_max) * 2.0_wp
          sb_pool_mt_domain(ix_microbial, ixP, ic, is, iblk)   = sb_pool_mt_domain(ix_microbial, ixN, ic, is, iblk) / k_fast_som_np
          sb_pool_mt_domain(ix_microbial, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) * hlp_frac_c13
          sb_pool_mt_domain(ix_microbial, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) * hlp_frac_c14
          sb_pool_mt_domain(ix_microbial, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_microbial, ixN, ic, is, iblk) &
            &                                                    / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! ..
          sb_pool_mt_domain(ix_residue, ixN, ic, is, iblk)   = sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) / k_slow_som_cn
          sb_pool_mt_domain(ix_residue, ixP, ic, is, iblk)   = sb_pool_mt_domain(ix_residue, ixN, ic, is, iblk) / k_slow_som_np
          sb_pool_mt_domain(ix_residue, ixC13, ic, is, iblk) = sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) * hlp_frac_c13
          sb_pool_mt_domain(ix_residue, ixC14, ic, is, iblk) = sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) * hlp_frac_c14
          sb_pool_mt_domain(ix_residue, ixN15, ic, is, iblk) = sb_pool_mt_domain(ix_residue, ixN, ic, is, iblk) &
            &                                                  / (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          ! ..
          nh4_n15_solute(ic, is, iblk) = nh4_solute(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          no3_n15_solute(ic, is, iblk) = no3_solute(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        END DO
      END DO
    END DO

    !> Soil P pools intialisation
    !>
    ! initialise mineral soil P pools
    DO iblk = 1,nblks
      CALL init_soil_p_pools( &
        & nproma, &                         ! in
        & nsoil_w, &
        & num_sl_above_bedrock(:, iblk), &
        & soil_lay_depth_center_sl(:,:,iblk), &
        & flag_sb_double_langmuir, &
        & soil_depth_sl(:,:,iblk), &
        & bulk_dens_corr_sl(:,:,iblk), &
        & root_fraction_sl(:,:,iblk), &
        & soil_p_labile(:,:,iblk), &
        & soil_p_slow(:,:,iblk), &
        & soil_p_occluded(:,:,iblk), &
        & soil_p_primary(:,:,iblk), &
        & soil_p_depth, &
        & p_labile_slow_site(:,:,iblk), &   ! in
        & po4_solute(:,:,iblk), &           ! inout
        & po4_assoc_fast(:,:,iblk), &
        & po4_assoc_slow(:,:,iblk), &
        & po4_occluded(:,:,iblk), &
        & po4_primary(:,:,iblk), &
        & qmax_po4(:,:,iblk), &
        & km_adsorpt_po4_sl(:,:,iblk), &
        & qmax_fast_po4(:,:,iblk), &
        & qmax_slow_po4(:,:,iblk), &
        & km_fast_po4(:,:,iblk), &
        & km_slow_po4(:,:,iblk) )           ! inout
    END DO

    !> Temporary debug dump after the full simple_1d initialisation is complete so the output
    !> contains the final C, N, and P pool sizes for both legacy and analytical_ss options.
    CALL dump_analytical_ss_pools('sb_analytical_init_debug_simple_sm.txt', nproma, nsoil_w, nblks, &
      & num_sl_above_bedrock, sb_pool_mt_domain)

    ! de-allocate local allocatable var
    DEALLOCATE(p_labile_slow_site)
  END SUBROUTINE sb_init_simple_sm

  ! ======================================================================================================= !
  !>init SB_ with the jsm - jena soil model
  !>
  SUBROUTINE sb_init_jsm( &
    & nproma, &
    & nsoil_w, &
    & nblks, &
    & num_sl_above_bedrock, &
    & lctlib_growthform, &
    & sb_model_scheme, &
    & sb_init_scheme, &
    & flag_sb_double_langmuir, &
    & qmax_org_fine_particle, &
    & soil_p_labile, &
    & soil_p_slow, &
    & soil_p_occluded, &
    & soil_p_primary, &
    & soil_p_depth, &
    & soil_depth_sl, &
    & soil_lay_width_sl, &
    & soil_lay_depth_ubound_sl, &
    & soil_lay_depth_center_sl, &
    & hlp_frac_c13, &
    & hlp_frac_c14, &
    & root_fraction_sl, &
    & clay_sl, &
    & silt_sl, &
    & vol_porosity_sl, &
    & bulk_dens_sl, &
    & ph_sl, &
    & qmax_po4_min_sl, &
    & qmax_po4_om_sl, &
    & volume_min_sl, &
    & po4_solute, &
    & bulk_carbon_max, &
    & surface_burial_frac, &
    & bulk_soil_carbon_sl, &
    & soil_litter_carbon_sl, &
    & qmax_po4, &
    & qmax_nh4, &
    & qmax_fast_po4, &
    & qmax_slow_po4, &
    & km_fast_po4, &
    & km_slow_po4, &
    & km_adsorpt_po4_sl, &
    & po4_assoc_fast, &
    & po4_assoc_slow, &
    & po4_occluded, &
    & Qmax_AlFe_cor, &
    & po4_primary, &
    & sb_pool_mt_domain, &
    & qmax_org_min_sl, &
    & qmax_nh4_min_sl, &
    & nh4_assoc, &
    & nh4_solute, &
    & no3_solute, &
    & nh4_n15_assoc, &
    & nh4_n15_solute, &
    & no3_n15_solute, &
    & bulk_dens_corr_sl, &
    & qmax_org, &
    & km_adsorpt_nh4_sl, &
    & rtm_depolymerisation, &
    & rmm_depolymerisation, &
    & rtm_mic_uptake, &
    & rmm_mic_uptake, &
    & rtm_sorption, &
    & rmm_sorption, &
    & rtm_desorption, &
    & rmm_desorption, &
    & rtm_hsc, &
    & rmm_hsc, &
    & k_bioturb, &
    & enzyme_frac_poly, &
    & enzyme_frac_residue, &
    & enzyme_frac_poly_c, &
    & enzyme_frac_poly_n, &
    & enzyme_frac_poly_p, &
    & enzyme_frac_AP)

    USE mo_q_sb_jsm_processes,            ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, &
      &                                         calc_qmax_bulk_density_correction, calc_Psorption_parameter, &
      &                                         calc_fast_po4
    USE mo_q_sb_jsm_transport,            ONLY: calc_bioturbation_rate
    USE mo_sb_util,                       ONLY: calc_qmax_texture
    USE mo_isotope_util,                  ONLY: calc_mixing_ratio_N15N14
    USE mo_veg_constants,                 ONLY: ITREE
    USE mo_sb_constants                   ! ..
    USE mo_jsb_physical_constants,        ONLY: molar_mass_C
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                    INTENT(in)    :: nproma
    INTEGER,                    INTENT(in)    :: nsoil_w
    INTEGER,                    INTENT(in)    :: nblks
    REAL(wp),                   INTENT(in)    :: num_sl_above_bedrock(:,:)
    INTEGER,                    INTENT(in)    :: lctlib_growthform
    CHARACTER(len=*),           INTENT(in)    :: sb_model_scheme                  !< SB_ config
    CHARACTER(len=*),           INTENT(in)    :: sb_init_scheme                   !< SB_ config
    LOGICAL,                    INTENT(in)    :: flag_sb_double_langmuir          !< SB_ config
    REAL(wp),                   INTENT(in)    :: qmax_org_fine_particle(:,:,:)    !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_labile(:,:,:)             !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_slow(:,:,:)               !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_occluded(:,:,:)           !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_primary(:,:,:)            !< bc_quincy_soil input file
    REAL(wp),                   INTENT(in)    :: soil_p_depth                     !< SB_ config
    REAL(wp),                   INTENT(in)    :: soil_depth_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_width_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_ubound_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: soil_lay_depth_center_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: hlp_frac_c13
    REAL(wp),                   INTENT(in)    :: hlp_frac_c14
    REAL(wp),                   INTENT(in)    :: root_fraction_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: clay_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: silt_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: vol_porosity_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: bulk_dens_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: ph_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4_om_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: volume_min_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_solute(:,:,:)
    REAL(wp),                   INTENT(inout) :: bulk_carbon_max(:,:,:)
    REAL(wp),                   INTENT(inout) :: surface_burial_frac(:,:,:)
    REAL(wp),                   INTENT(inout) :: bulk_soil_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: soil_litter_carbon_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_nh4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: qmax_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_fast_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_slow_po4(:,:,:)
    REAL(wp),                   INTENT(inout) :: km_adsorpt_po4_sl(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_fast(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_assoc_slow(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_occluded(:,:,:)
    REAL(wp),                   INTENT(inout) :: Qmax_AlFe_cor(:,:,:)
    REAL(wp),                   INTENT(inout) :: po4_primary(:,:,:)
    REAL(wp),                   INTENT(inout) :: sb_pool_mt_domain(:,:,:,:,:)          !< bgc_material sb_pool
    REAL(wp),                   INTENT(out)   :: qmax_org_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_nh4_min_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_assoc(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_assoc(:,:,:)
    REAL(wp),                   INTENT(out)   :: nh4_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: no3_n15_solute(:,:,:)
    REAL(wp),                   INTENT(out)   :: bulk_dens_corr_sl(:,:,:)
    REAL(wp),                   INTENT(out)   :: qmax_org(:,:,:)
    REAL(wp),                   INTENT(out)   :: km_adsorpt_nh4_sl(:,:,:)
    REAL(wp),                   INTENT(in)    :: rtm_depolymerisation(:,:,:)
    REAL(wp),                   INTENT(in)    :: rmm_depolymerisation(:,:,:)
    REAL(wp),                   INTENT(in)    :: rtm_mic_uptake(:,:,:)
    REAL(wp),                   INTENT(in)    :: rmm_mic_uptake(:,:,:)
    REAL(wp),                   INTENT(in)    :: rtm_sorption(:,:,:)
    REAL(wp),                   INTENT(in)    :: rmm_sorption(:,:,:)
    REAL(wp),                   INTENT(in)    :: rtm_desorption(:,:,:)
    REAL(wp),                   INTENT(in)    :: rmm_desorption(:,:,:)
    REAL(wp),                   INTENT(in)    :: rtm_hsc(:,:,:)
    REAL(wp),                   INTENT(in)    :: rmm_hsc(:,:,:)
    REAL(wp),                   INTENT(out)   :: k_bioturb(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_residue(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly_c(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly_n(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_poly_p(:,:,:)
    REAL(wp),                   INTENT(out)   :: enzyme_frac_AP(:,:,:)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), ALLOCATABLE       :: p_labile(:,:,:)
    REAL(wp), ALLOCATABLE       :: p_labile_slow_site(:,:,:)
    REAL(wp), ALLOCATABLE       :: arr_hlp7(:,:,:)
    REAL(wp), ALLOCATABLE       :: arr_hlp8(:,:,:)
    REAL(wp), ALLOCATABLE       :: arr_hlp9(:,:,:)
    REAL(wp)                    :: frac_mineral_assoc_soil_c  !< fraction of mineral-associated C in the soil
    INTEGER                     :: iblk, is, ic
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':sb_init_jsm'
    IF (debug_on()) CALL message(TRIM(routine), 'Starting ...')
    ! ----------------------------------------------------------------------------------------------------- !

    !> init out var
    !>
    qmax_org_min_sl(:,:,:)      = 0.0_wp
    qmax_nh4_min_sl(:,:,:)      = 0.0_wp
    nh4_assoc(:,:,:)            = 0.0_wp
    nh4_solute(:,:,:)           = 0.0_wp
    no3_solute(:,:,:)           = 0.0_wp
    nh4_n15_assoc(:,:,:)        = 0.0_wp
    nh4_n15_solute(:,:,:)       = 0.0_wp
    no3_n15_solute(:,:,:)       = 0.0_wp
    bulk_dens_corr_sl(:,:,:)    = 0.0_wp
    qmax_org(:,:,:)             = 0.0_wp
    bulk_carbon_max(:,:,:)      = 0.0_wp
    surface_burial_frac(:,:,:)  = 1.0_wp
    km_adsorpt_nh4_sl(:,:,:)    = 0.0_wp
    k_bioturb(:,:,:)            = 0.0_wp
    enzyme_frac_poly(:,:,:)     = 0.0_wp
    enzyme_frac_residue(:,:,:)  = 0.0_wp
    enzyme_frac_poly_c(:,:,:)   = 0.0_wp
    enzyme_frac_poly_n(:,:,:)   = 0.0_wp
    enzyme_frac_poly_p(:,:,:)   = 0.0_wp
    enzyme_frac_AP(:,:,:)       = 0.0_wp

    !>0.9 allocate local var
    !>
    ALLOCATE(p_labile(nproma, nsoil_w, nblks))
    ALLOCATE(p_labile_slow_site(nproma, nsoil_w, nblks))
    ALLOCATE(arr_hlp7(nproma, nsoil_w, nblks))
    ALLOCATE(arr_hlp8(nproma, nsoil_w, nblks))
    ALLOCATE(arr_hlp9(nproma, nsoil_w, nblks))

    CALL calc_bulk_carbon_max_profile(nproma, nsoil_w, nblks, num_sl_above_bedrock, soil_lay_depth_center_sl, bulk_carbon_max)

    SELECT CASE(TRIM(sb_init_scheme))
    CASE('legacy')
    !>1.0 init
    !>
    ! initialisation of litter pools
    sb_pool_mt_domain(ix_soluable_litter, ixC, :, 1, :) &
      & = def_sb_pool_metabol_litter /  molar_mass_C / soil_lay_width_sl(:,1,:) &
      &                                             * (frac_prod_not_root + (1.0_wp - frac_prod_not_root) &
      &                                             * root_fraction_sl(:,1,:))
    sb_pool_mt_domain(ix_polymeric_litter, ixC, :, 1, :) &
      & = def_sb_pool_jsm_struct_litter / molar_mass_C / soil_lay_width_sl(:,1,:) &
      &                                             * (frac_prod_not_root + (1.0_wp - frac_prod_not_root) &
      &                                             * root_fraction_sl(:,1,:))
    IF (lctlib_growthform == ITREE) THEN
      sb_pool_mt_domain(ix_woody_litter, ixC, :, 1, :)  = def_sb_pool_woody_litter / molar_mass_C / soil_lay_width_sl(:,1,:)
    END IF
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)  = def_sb_pool_metabol_litter /  molar_mass_C &
            &                                                 / soil_lay_width_sl(ic, is, iblk) &
            &                                                 * (1.0_wp - frac_prod_not_root) * root_fraction_sl(ic, is, iblk)
          sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = def_sb_pool_jsm_struct_litter / molar_mass_C &
            &                                                 / soil_lay_width_sl(ic, is, iblk) &
            &                                                 * (1.0_wp - frac_prod_not_root) * root_fraction_sl(ic, is, iblk)
        END DO
      END DO
    END DO

    ! for each soil layer
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          IF (soil_lay_width_sl(ic, is, iblk) > 0._wp) THEN
            sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) = MAX(def_sb_pool_microbial_biomass &
              &                                          / soil_lay_width_sl(ic, is, iblk) &
              &                                          * root_fraction_sl(ic, is, iblk), min_mic_biomass)
            sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) = def_sb_pool_dom / soil_lay_width_sl(ic, is, iblk) &
              &                                          * root_fraction_sl(ic, is, iblk)
            sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk) &
                & = k_sb_pool_dom_assoc_c_corr * sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)
              sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) = &
                &                               def_sb_pool_microbial_necromass / soil_lay_width_sl(ic, is, iblk) &
              &                                              * root_fraction_sl(ic, is, iblk)
            ! calc fraction of mineral-associated C in the soil
            ! initialize the mineral-associated organic carbon (MOC) of soil layer by a fraction coefficient and the maximum sortion capacity
            ! assuming the fraction coefficient frac_mineral_assoc_soil_c decreases with soil depth
              frac_mineral_assoc_soil_c  = EXP(-is / (NINT(num_sl_above_bedrock(ic, iblk)) * k_sl_conv_moc_corr) ) * 0.5_wp
            ! deep soil is initialised with very low MOC, corrected by k_deep_soil_moc_corr
            IF (soil_lay_depth_ubound_sl(ic, is, iblk) > z_deep_soil_moc_corr) THEN
                frac_mineral_assoc_soil_c = frac_mineral_assoc_soil_c &
                &                           * MAX(0.0001_wp, 10._wp**(-soil_lay_depth_ubound_sl(ic, is, iblk)))
            END IF
            sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk) = qmax_org_fine_particle(ic, is, iblk) &
              &                                                      * (clay_sl(ic, is, iblk) + silt_sl(ic, is, iblk)) &
              &                                                      * bulk_dens_sl(ic, is, iblk) * frac_mineral_assoc_soil_c
            ! To avoid over-saturation of organic matter in soil layer
            ! calculate the fast C and stable C content (saturation at 1._wp) for each layer
            ! arr_hlp7: fast C including litter, microbial biomass, dissolved organic matter (DOM), and mineral-associated DOM
            ! arr_hlp8: stable C including microbial residue and mineral-associated residue
            IF (is == 1) THEN
              arr_hlp7(ic, is, iblk) = (sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) + &
              &                     sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)) * &
              &                     molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
              arr_hlp8(ic, is, iblk) = (sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)) * &
              &                     molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
            ELSE
              arr_hlp7(ic, is, iblk) = (sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)   &
              &                     + sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        &
              &                     + sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)              &
              &                     + sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk))       &
              &                     * molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
              arr_hlp8(ic, is, iblk) = (sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
              &                     + sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)) * &
              &                     molar_mass_C / carbon_per_dryweight_SOM / 1000._wp / rho_bulk_org
            ENDIF
            ! when soil is over over-saturated with fast C:
            !     downscale each pool in fast C
            !     set the stable C pools to zero
            IF (arr_hlp7(ic, is, iblk) > 1._wp ) THEN
              sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)  = &
                &                 sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk)  &
                &                 * 1._wp / arr_hlp7(ic, is, iblk)
              sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = &
                &                 sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
                &                 * 1._wp / arr_hlp7(ic, is, iblk)
              sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        = &
                &                 sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk)        &
                &                 * 1._wp / arr_hlp7(ic, is, iblk)
              sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)              = &
                &                 sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk)              &
                &                 * 1._wp / arr_hlp7(ic, is, iblk)
              sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)        = &
                &                 sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk)        &
                &                 * 1._wp / arr_hlp7(ic, is, iblk)
              sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk)          = 0._wp
              sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)    = 0._wp
            ! when soil is over-saturated with fast plus stable C, redistribute the stable C fraction (1 - arr_hlp7) to
            ! microbial residue and mineral-associated residue assuming:
            !     maximum sorption capacity (arr_hlp9) defined by fact_repartition_om_init_cond1 of soil mineral content
            !     total stable C pool (arr_hlp8)
            ELSE IF ((arr_hlp7(ic, is, iblk) + arr_hlp8(ic, is, iblk)) > 1._wp) THEN
              arr_hlp9(ic, is, iblk) = (1._wp - arr_hlp7(ic, is, iblk)) * fact_repartition_om_init_cond1 * &
                  &                     qmax_org_fine_particle(ic, is, iblk) * frac_mineral_assoc_soil_c * &
                &                     (clay_sl(ic, is, iblk) + silt_sl(ic, is, iblk)) * bulk_dens_sl(ic, is, iblk)
              arr_hlp8(ic, is, iblk) = (1._wp - arr_hlp7(ic, is, iblk)) * (1._wp - fact_repartition_om_init_cond1) * &
                &                     rho_bulk_org * 1000._wp * carbon_per_dryweight_SOM / molar_mass_C
              IF (arr_hlp8(ic, is, iblk) > arr_hlp9(ic, is, iblk)) THEN
                sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) = &
                  &                                                 sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) &
                  &                                                 + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                  &                                                 * fact_repartition_om_init_cond2
                sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = &
                  &                                                 sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
                  &                                                 + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                  &                                                 * (1._wp - fact_repartition_om_init_cond2)
                arr_hlp8(ic, is, iblk) = arr_hlp9(ic, is, iblk)
              ENDIF
              sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk) = calc_fast_po4(arr_hlp8(ic, is, iblk), &
                &                                                      arr_hlp9(ic, is, iblk), &
                &                                                      k_desorpt_det / k_adsorpt_det)
              sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
                &   = arr_hlp8(ic, is, iblk) - sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)
            ! when soil is not over-saturated with OM, check if the void soil space is enough to adsorb the OM
            !     arr_hlp9: the void space in soil layer, assuming as soil mineral
            !     arr_hlp8: the predefined mineral-associated OM content
            ELSE
              arr_hlp9(ic, is, iblk) = (1._wp - arr_hlp7(ic, is, iblk) - arr_hlp8(ic, is, iblk)) &
                &                     * qmax_org_fine_particle(ic, is, iblk) * frac_mineral_assoc_soil_c &
              &                     * (clay_sl(ic, is, iblk) + silt_sl(ic, is, iblk)) * bulk_dens_sl(ic, is, iblk)
              arr_hlp8(ic, is, iblk) &
              & = (sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) + sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk))
              IF (arr_hlp8(ic, is, iblk) > arr_hlp9(ic, is, iblk)) THEN
                sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) = &
                  &                                               sb_pool_mt_domain(ix_soluable_litter, ixC, ic, is, iblk) &
                  &                                               + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                  &                                               * fact_repartition_om_init_cond2
                sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) = &
                  &                                               sb_pool_mt_domain(ix_polymeric_litter, ixC, ic, is, iblk) &
                  &                                               + (arr_hlp8(ic, is, iblk) - arr_hlp9(ic, is, iblk)) &
                  &                                               * (1._wp - fact_repartition_om_init_cond2)
                arr_hlp8(ic, is, iblk) = arr_hlp9(ic, is, iblk)
              ENDIF
              sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk) = calc_fast_po4(arr_hlp8(ic, is, iblk), &
                &                                                      arr_hlp9(ic, is, iblk), &
                &                                                      k_desorpt_det / k_adsorpt_det)
              sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) &
                &   = arr_hlp8(ic, is, iblk) - sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk)
            ENDIF
          ELSE
            sb_pool_mt_domain(ix_microbial, ixC, ic, is, iblk) = 0._wp
            sb_pool_mt_domain(ix_dom, ixC, ic, is, iblk) = 0._wp
            sb_pool_mt_domain(ix_dom_assoc, ixC, ic, is, iblk) = 0._wp
            sb_pool_mt_domain(ix_residue, ixC, ic, is, iblk) = 0._wp
            sb_pool_mt_domain(ix_residue_assoc, ixC, ic, is, iblk) = 0._wp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    CASE('analytical_ss')
      CALL calc_analytical_ss_jsm( &
        & nproma, &
        & nsoil_w, &
        & nblks, &
        & num_sl_above_bedrock, &
        & qmax_org_fine_particle, &
        & soil_lay_depth_center_sl, &
        & soil_lay_depth_ubound_sl, &
        & silt_sl, &
        & clay_sl, &
        & bulk_dens_sl, &
        & bulk_carbon_max, &
        & qmax_org, &
        & rtm_depolymerisation, &
        & rmm_depolymerisation, &
        & rtm_mic_uptake, &
        & rmm_mic_uptake, &
        & rtm_sorption, &
        & rmm_sorption, &
        & rtm_desorption, &
        & rmm_desorption, &
        & rtm_hsc, &
        & rmm_hsc, &
        & sb_pool_mt_domain, &
        & surface_burial_frac, &
        & enzyme_frac_poly, &
        & enzyme_frac_residue, &
        & enzyme_frac_poly_c, &
        & enzyme_frac_poly_n, &
        & enzyme_frac_poly_p, &
        & enzyme_frac_AP)
    CASE DEFAULT
      CALL finish(TRIM(routine), 'Invalid setting of sb_init_scheme (legacy / analytical_ss).')
    END SELECT

    ! initial values of soluble nh4 and no3
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          nh4_solute(ic, is, iblk)      = nh4_solute_prescribe / soil_lay_depth_ubound_sl(ic, nsoil_w, iblk)
          nh4_n15_solute(ic, is, iblk)  = nh4_solute(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
          no3_solute(ic, is, iblk)      = no3_solute_prescribe / soil_lay_depth_ubound_sl(ic, nsoil_w, iblk)
          no3_n15_solute(ic, is, iblk)  = no3_solute(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        END DO
      END DO
    END DO
    ! inital bulk density, Qmax_org & Qmax_po4 based on current soil carbon content
    DO iblk = 1,nblks
      CALL calc_bulk_soil_carbon( &
        & nproma, &                           ! in
        & nsoil_w, &
        & num_sl_above_bedrock(:,iblk), &
        & TRIM(sb_model_scheme), &
        & sb_pool_mt_domain(:,:,:,:,iblk), &  ! in
        & bulk_carbon_max(:,:,iblk), &        ! in
        & surface_burial_frac(:,:,iblk), &    ! in
        & bulk_soil_carbon_sl(:,:,iblk), &    ! inout
        & soil_litter_carbon_sl(:,:,iblk), &
        & volume_min_sl(:,:,iblk) )           ! inout
      bulk_dens_corr_sl(:,:,iblk) = calc_bulk_density_correction(bulk_soil_carbon_sl(:,:,iblk), &
        &                                                        soil_litter_carbon_sl(:,:,iblk), &
        &                                                        bulk_dens_sl(:,:,iblk))
      IF (TRIM(sb_init_scheme) /= 'analytical_ss') THEN
      !> po4/OC adsorption parameter intialisation
      !! OC adsorption
      qmax_org_min_sl(:,:,iblk) = calc_qmax_texture(qmax_org_fine_particle(:,:,iblk), &
        &                                           silt_sl(:,:,iblk), &
        &                                           clay_sl(:,:,iblk))
      !! deep soil qmax correction for OM
      !! if we need to fix it also for NH4, and PO4, we need to modify calc_qmax_texture with the following code
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
            frac_mineral_assoc_soil_c  = EXP(-is / (INT(num_sl_above_bedrock(ic, iblk)) * k_sl_conv_moc_corr) ) * 0.5_wp
          qmax_org_min_sl(ic, is, iblk) = qmax_org_min_sl(ic, is, iblk) * frac_mineral_assoc_soil_c
          IF (soil_lay_depth_ubound_sl(ic, is, iblk) > z_deep_soil_moc_corr) THEN
              qmax_org_min_sl(ic, is, iblk) = qmax_org_min_sl(ic, is, iblk) &
              &                               * MAX(0.0001_wp, 10._wp**(-soil_lay_depth_ubound_sl(ic, is, iblk)))
          ENDIF
        END DO
      END DO

      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &    ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0._wp, &
        & qmax_org_min_sl(:,:,iblk), &        ! in
        & qmax_org(:,:,iblk) )                ! out
      END IF
    END DO

    !! po4 adsorption
    ! Correct the first layer Qmax_AlFe_cor(:,:,:) with the woody litter
    Qmax_AlFe_cor(:,1,:) = Qmax_AlFe_cor(:,1,:) * 1._wp / ((1._wp - volume_min_sl(:,1,:)) &
      &                    * (bulk_soil_carbon_sl(:,1,:) + soil_litter_carbon_sl(:,1,:) &
      &                    + sb_pool_mt_domain(ix_woody_litter, ixC, :, 1, :)) &
      &                    / (bulk_soil_carbon_sl(:,1,:) + soil_litter_carbon_sl(:,1,:)) &
      &                    + volume_min_sl(:,1,:))
    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          CALL calc_Psorption_parameter( &
            & flag_sb_double_langmuir, &          ! in
            & clay_sl(ic, is, iblk), &
            & silt_sl(ic, is, iblk), &
            & vol_porosity_sl(ic, is, iblk) * 1000._wp, &
            & ph_sl(ic, is, iblk), &
            & bulk_soil_carbon_sl(ic, is, iblk), &
            & soil_litter_carbon_sl(ic, is, iblk), &
            & sb_pool_mt_domain(ix_woody_litter, ixC, ic, is, iblk), &
            & bulk_dens_sl(ic, is, iblk), &
            & volume_min_sl(ic, is, iblk), &
            & po4_solute_prescribe, &             ! in - po4_solute(ic, is, iblk)
            & qmax_po4_min_sl(ic, is, iblk), &        ! inout
            & qmax_po4_om_sl(ic, is, iblk), &
            & qmax_po4(ic, is, iblk), &
            & km_adsorpt_po4_sl(ic, is, iblk), &
            & qmax_fast_po4(ic, is, iblk), &
            & qmax_slow_po4(ic, is, iblk), &
            & km_fast_po4(ic, is, iblk), &
            & km_slow_po4(ic, is, iblk), &
            & Qmax_AlFe_cor(ic, is, iblk) )           ! inout
        END DO
      END DO
    END DO

    ! calculate the partition coefficient between labile and slow P pool based on
    ! their Qmax (double-Langmuir) or global average values (single-Langmuir)
    IF (flag_sb_double_langmuir) THEN
      p_labile_slow_site(:,:,:) = eps4
      DO iblk = 1,nblks
        DO ic = 1,nproma
          DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
            IF (qmax_po4(ic, is, iblk) > eps8) THEN
              p_labile_slow_site(ic, is, iblk) = MAX(qmax_fast_po4(ic, is, iblk) / qmax_po4(ic, is, iblk), eps4)
            END IF
          END DO
        END DO
      END DO
    ELSE
      p_labile_slow_site(:,:,:) = p_labile_slow_global_avg
    END IF

    !> nh4 adsorption parameter intialisation
    !!
    qmax_nh4_min_sl(:,:,:) = calc_qmax_texture(qmax_nh4_clay, &
      &                                        0.0_wp, &
      &                                        clay_sl(:,:,:))
    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &   ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & 0.0_wp, &
        & qmax_nh4_min_sl(:,:,iblk), &       ! in
        & qmax_nh4(:,:,iblk) )               ! out
    END DO

    DO iblk = 1,nblks
      CALL calc_qmax_bulk_density_correction( &
        & bulk_soil_carbon_sl(:,:,iblk), &   ! in
        & volume_min_sl(:,:,iblk), &
        & bulk_dens_sl(:,:,iblk), &
        & km_adsorpt_OM_nh4, &
        & km_adsorpt_mineral_nh4, &          ! in
        & km_adsorpt_nh4_sl(:,:,iblk) )      ! out
    END DO

    DO iblk = 1,nblks
      DO ic = 1,nproma
        DO is = 1, NINT(num_sl_above_bedrock(ic, iblk))
          nh4_assoc(ic, is, iblk)     = qmax_nh4(ic, is, iblk) * nh4_solute (ic, is, iblk) &
            &                           / (nh4_solute(ic, is, iblk) + km_adsorpt_nh4_sl(ic, is, iblk))
          nh4_n15_assoc(ic, is, iblk) = nh4_assoc(ic, is, iblk) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        END DO
      END DO
    END DO

    !> Soil P pools intialisation
    !>
    ! initialise mineral soil P pools
    DO iblk = 1,nblks
      CALL init_soil_p_pools( &
        & nproma, &                               ! in
        & nsoil_w, &
        & num_sl_above_bedrock(:, iblk), &
        & soil_lay_depth_center_sl(:,:,iblk), &
        & flag_sb_double_langmuir, &
        & soil_depth_sl(:,:,iblk), &
        & bulk_dens_corr_sl(:,:,iblk), &
        & root_fraction_sl(:,:,iblk), &
        & soil_p_labile(:,:,iblk), &
        & soil_p_slow(:,:,iblk), &
        & soil_p_occluded(:,:,iblk), &
        & soil_p_primary(:,:,iblk), &
        & soil_p_depth, &
        & p_labile_slow_site(:,:,iblk), &         ! in
        & po4_solute(:,:,iblk), &                 ! inout
        & po4_assoc_fast(:,:,iblk), &
        & po4_assoc_slow(:,:,iblk), &
        & po4_occluded(:,:,iblk), &
        & po4_primary(:,:,iblk), &
        & qmax_po4(:,:,iblk), &
        & km_adsorpt_po4_sl(:,:,iblk), &
        & qmax_fast_po4(:,:,iblk), &
        & qmax_slow_po4(:,:,iblk), &
        & km_fast_po4(:,:,iblk), &
        & km_slow_po4(:,:,iblk) )                 ! inout
    END DO

    ! initialise diffusion by bioturbation
    DO iblk = 1,nblks
      k_bioturb(:,:,iblk) = calc_bioturbation_rate(nproma, &
        &                                          nsoil_w, &
        &                                          num_sl_above_bedrock(:,iblk), &
        &                                          soil_depth_sl(:,:,iblk), &
        &                                          bulk_dens_corr_sl(:,:,iblk), &
        &                                          root_fraction_sl(:,:,iblk))
    END DO

    ! initialisation of fractional enzyme allocation (range [0,1])
    IF (TRIM(sb_init_scheme) == 'legacy') THEN
    enzyme_frac_poly(:,:,:)     = init_enzyme_fraction
    enzyme_frac_residue(:,:,:)  = init_enzyme_fraction
    enzyme_frac_poly_c(:,:,:)   = init_enzyme_fraction
    enzyme_frac_poly_n(:,:,:)   = init_enzyme_fraction
    enzyme_frac_poly_p(:,:,:)   = init_enzyme_fraction
    enzyme_frac_AP(:,:,:)       = min_enzyme_fraction
    END IF

    ! isotopic signals in all pools
    ! microbial, DOM, aDOM, residue, and aRes pools
    sb_pool_mt_domain(ix_soluable_litter, ixN, :, :, :) &
      & = sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) / k_init_soluable_litter_cn
    sb_pool_mt_domain(ix_soluable_litter, ixP, :, :, :) &
      &  = sb_pool_mt_domain(ix_soluable_litter, ixN, :, :, :) / k_init_litter_np
    sb_pool_mt_domain(ix_soluable_litter, ixC13, :, :, :) = sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_soluable_litter, ixC14, :, :, :) = sb_pool_mt_domain(ix_soluable_litter, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_soluable_litter, ixN15, :, :, :) = sb_pool_mt_domain(ix_soluable_litter, ixN, :, :, :) / &
      &                                              (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_polymeric_litter, ixN, :, :, :) &
      & = sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) / k_init_polymeric_litter_cn
    sb_pool_mt_domain(ix_polymeric_litter, ixP, :, :, :) &
      & = sb_pool_mt_domain(ix_polymeric_litter, ixN, :, :, :) / k_init_litter_np
    sb_pool_mt_domain(ix_polymeric_litter, ixC13, :, :, :) = sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_polymeric_litter, ixC14, :, :, :) = sb_pool_mt_domain(ix_polymeric_litter, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_polymeric_litter, ixN15, :, :, :) = sb_pool_mt_domain(ix_polymeric_litter, ixN, :, :, :) / &
      &                                               (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_woody_litter, ixN, :, :, :) &
      & = sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) / k_init_woody_litter_cn
    sb_pool_mt_domain(ix_woody_litter, ixP, :, :, :)   = sb_pool_mt_domain(ix_woody_litter, ixN, :, :, :) / k_init_litter_np
    sb_pool_mt_domain(ix_woody_litter, ixC13, :, :, :) = sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_woody_litter, ixC14, :, :, :) = sb_pool_mt_domain(ix_woody_litter, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_woody_litter, ixN15, :, :, :) = sb_pool_mt_domain(ix_woody_litter, ixN, :, :, :) / &
      &                                           (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_microbial, ixN, :, :, :)   = sb_pool_mt_domain(ix_microbial, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_microbial, ixP, :, :, :)   = sb_pool_mt_domain(ix_microbial, ixN, :, :, :) / microbial_np
    sb_pool_mt_domain(ix_microbial, ixC13, :, :, :) = sb_pool_mt_domain(ix_microbial, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_microbial, ixC14, :, :, :) = sb_pool_mt_domain(ix_microbial, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_microbial, ixN15, :, :, :) = sb_pool_mt_domain(ix_microbial, ixN, :, :, :) / &
      &                                        (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_residue, ixN, :, :, :)   = sb_pool_mt_domain(ix_residue, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_residue, ixP, :, :, :)   = sb_pool_mt_domain(ix_residue, ixN, :, :, :) / k_slow_som_np
    sb_pool_mt_domain(ix_residue, ixC13, :, :, :) = sb_pool_mt_domain(ix_residue, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_residue, ixC14, :, :, :) = sb_pool_mt_domain(ix_residue, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_residue, ixN15, :, :, :) = sb_pool_mt_domain(ix_residue, ixN, :, :, :) / &
      &                                      (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_residue_assoc, ixN, :, :, :)   = sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_residue_assoc, ixP, :, :, :)   = sb_pool_mt_domain(ix_residue_assoc, ixN, :, :, :) / k_slow_som_np
    sb_pool_mt_domain(ix_residue_assoc, ixC13, :, :, :) = sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_residue_assoc, ixC14, :, :, :) = sb_pool_mt_domain(ix_residue_assoc, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_residue_assoc, ixN15, :, :, :) = sb_pool_mt_domain(ix_residue_assoc, ixN, :, :, :) / &
      &                                            (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_dom, ixN, :, :, :)   = sb_pool_mt_domain(ix_dom, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_dom, ixP, :, :, :)   = sb_pool_mt_domain(ix_dom, ixN, :, :, :) / microbial_np
    sb_pool_mt_domain(ix_dom, ixC13, :, :, :) = sb_pool_mt_domain(ix_dom, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_dom, ixC14, :, :, :) = sb_pool_mt_domain(ix_dom, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_dom, ixN15, :, :, :) = sb_pool_mt_domain(ix_dom, ixN, :, :, :) / &
      &                                  (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    sb_pool_mt_domain(ix_dom_assoc, ixN, :, :, :)   = sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) / microbial_cn
    sb_pool_mt_domain(ix_dom_assoc, ixP, :, :, :)   = sb_pool_mt_domain(ix_dom_assoc, ixN, :, :, :) / microbial_np
    sb_pool_mt_domain(ix_dom_assoc, ixC13, :, :, :) = sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_dom_assoc, ixC14, :, :, :) = sb_pool_mt_domain(ix_dom_assoc, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_dom_assoc, ixN15, :, :, :) = sb_pool_mt_domain(ix_dom_assoc, ixN, :, :, :) / &
      &                                        (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    ! mycorrhiza
    sb_pool_mt_domain(ix_mycorrhiza, ixC13, :, :, :) = sb_pool_mt_domain(ix_mycorrhiza, ixC, :, :, :) * hlp_frac_c13
    sb_pool_mt_domain(ix_mycorrhiza, ixC14, :, :, :) = sb_pool_mt_domain(ix_mycorrhiza, ixC, :, :, :) * hlp_frac_c14
    sb_pool_mt_domain(ix_mycorrhiza, ixN15, :, :, :) = sb_pool_mt_domain(ix_mycorrhiza, ixN, :, :, :) / &
      &                                         (1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))

    !> Temporary debug dump after the full JSM initialisation is complete so the output
    !> contains the final C, N, and P pool sizes for both legacy and analytical_ss options.
    CALL dump_analytical_ss_pools('sb_analytical_init_debug_jsm.txt', nproma, nsoil_w, nblks, &
      & num_sl_above_bedrock, sb_pool_mt_domain)

    ! de-allocate local allocatable var
    DEALLOCATE(p_labile)
    DEALLOCATE(p_labile_slow_site)
    DEALLOCATE(arr_hlp7)
    DEALLOCATE(arr_hlp8)
    DEALLOCATE(arr_hlp9)
  END SUBROUTINE sb_init_jsm

  ! ======================================================================================================= !
  !>Intialize soil P pools for all layers from the soil P maps (Yang et al. 2013)
  !> assuming that the soil P pool sizes are correlated with soil weight
  !>
  SUBROUTINE init_soil_p_pools(nc, &
                               nsoil_w, &
                               num_sl_above_bedrock, &
                               soil_lay_depth_center_sl, &
                               flag_sb_double_langmuir, &
                               soil_depth_sl, bulk_dens_sl, root_fraction_sl,&
                               soil_p_labile, soil_p_slow, soil_p_occluded, soil_p_primary, &
                               soil_p_depth, p_labile_slow_ratio, &
                               po4_solute, po4_assoc_fast, po4_assoc_slow, po4_occluded, po4_primary, &
                               qmax_po4,km_po4, qmax_fast_po4, km_fast_po4, qmax_slow_po4, km_slow_po4)

    USE mo_jsb_physical_constants, ONLY: molar_mass_P
    USE mo_sb_constants           ! e.g.: k1_p_pool_depth_corr, k2_p_pool_depth_corr, k3_p_pool_depth_corr, &
                                  !       k4_p_pool_depth_corr, k_p_pool_depth_corr_min_occlud, p_labile_slow_global_avg
    USE mo_q_sb_jsm_processes,     ONLY: calc_fast_po4

    INTEGER,                          INTENT(in)     :: nc                        !< dimensions
    INTEGER                         , INTENT(in)     :: nsoil_w                   !< number of soil layers (water)
    REAL(wp), DIMENSION(nc)         , INTENT(in)     :: num_sl_above_bedrock      !< number of soil layers with thickness > eps8
    REAL(wp), DIMENSION(nc,nsoil_w) , INTENT(in)     :: soil_lay_depth_center_sl  !< depth at the center of each soil layer
    LOGICAL                         , INTENT(in)     :: flag_sb_double_langmuir   !< T: double Langmuir; F: traditional Langmuir
    REAL(wp), DIMENSION(nc,nsoil_w) , INTENT(in)     :: soil_depth_sl, &      !< soil depth per soil layer, [m]
                                                        root_fraction_sl, &   !<
                                                        bulk_dens_sl, &       !< corrected soil bulk density per soil layer, [kg m-3]
                                                        p_labile_slow_ratio   !< ratio between labile and labile+slow P pools
    REAL(wp)                        , INTENT(in)     :: soil_p_depth          !< soil depth of the soil P maps
    REAL(wp), DIMENSION(nc,nsoil_w) , INTENT(in)     :: soil_p_labile, &      !< total labile P pool (solute + assoc_fast) of the soil P maps, [g m-2]
                                                        soil_p_slow, &        !< total slow P pool (secondary P) of the soil P maps, [g m-2]
                                                        soil_p_occluded, &    !< total occluded P pool of the soil P maps, [g m-2]
                                                        soil_p_primary        !< total primary P pool of the soil P maps, [g m-2]
    REAL(wp), DIMENSION(nc,nsoil_w) , INTENT(inout)  :: po4_solute, &         !< soluble Pi pool per soil layer, [mol m-3]
                                                        po4_assoc_fast, &     !< fast minerally associated PO4 pool per soil layer, [mol m-3]
                                                        po4_assoc_slow, &     !< slow minerally associated PO4 pool per soil layer, [mol m-3]
                                                        po4_occluded, &       !< occluded PO4 pool per soil layer, [mol m-3]
                                                        po4_primary, &        !< primary PO4 pool per soil layer, [mol m-3]
                                                        qmax_po4, &           !< maximum PO4 sorption capacity per soil layer, [mol m-3]
                                                        km_po4, &             !< half-saturation PO4 concentration for sorption per soil layer, [mol m-3]
                                                        qmax_fast_po4, &      !< maximum PO4 sorption capacity for po4_assoc_fast, [mol m-3]
                                                        km_fast_po4, &        !< half-saturation PO4 concentration for po4_assoc_fast, [mol m-3]
                                                        qmax_slow_po4, &      !< maximum PO4 sorption capacity for po4_assoc_slow, [mol m-3]
                                                        km_slow_po4           !< half-saturation PO4 concentration for po4_assoc_slow, [mol m-3]
    INTEGER                                          :: ic, is
    REAL(wp), DIMENSION(nc)                          :: rhlp1
    REAL(wp), DIMENSION(nc)                          :: soil_p_weight         !< total soil weight up to the soil_p_depth [kg]
    REAL(wp), DIMENSION(nc)                          :: soil_p_layer          !< up to soil_p_layer keep the P init map info
    REAL(wp), DIMENSION(nc)                          :: hlp1, hlp2
    REAL(wp), DIMENSION(nc)                          :: hlp3, hlp4
    LOGICAL,  DIMENSION(nc)                          :: aflag
    REAL(wp), DIMENSION(nc,nsoil_w)                  :: soil_weight_sl, &                 !< soil weight per soil layer
                                                        arr_hlp1
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_soil_p_pools'


    !> 0.8 init local variables
    !!
    soil_p_weight(:)    = eps8
    soil_p_layer(:)     = 1.0_wp
    hlp1(:)             = 0.0_wp
    hlp2(:)             = 0.0_wp
    aflag(:)            = .TRUE.
    soil_weight_sl(:,:) = 0.0_wp

    !> 0.9 init output variables
    !!
    po4_solute(:,:)     = 0._wp
    po4_assoc_fast(:,:) = 0._wp
    po4_assoc_slow(:,:) = 0._wp
    po4_occluded(:,:)   = 0._wp
    po4_primary(:,:)    = 0._wp

    !> 1.0 calculate the total soil weight up to the soil_p_depth,
    !! and the soil weight per soil layer
    !!
    !> NOTE at the moment the soil_p_* values are identical across soil layers, hence, using the value from sl 1
    rhlp1(:) = soil_p_slow(:,1) + soil_p_labile(:,1) + soil_p_occluded(:,1) + soil_p_primary(:,1)
    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        soil_weight_sl(ic,is) = soil_depth_sl(ic,is) * bulk_dens_sl(ic,is)
        IF (aflag(ic)) THEN
          hlp1(ic) = hlp1(ic) + soil_depth_sl(ic,is)
          hlp2(ic) = hlp2(ic) + soil_weight_sl(ic,is)
          IF (hlp1(ic) > soil_p_depth) THEN
            soil_p_weight(ic) = hlp2(ic) - soil_weight_sl(ic,is) * (hlp1(ic) - soil_p_depth) / soil_depth_sl(ic,is)
            aflag(ic)         = .FALSE.
            soil_p_layer(ic)  = is
          END IF
        END IF
      END DO
    END DO
    hlp4(:) = rhlp1(:) / soil_p_weight(:) ! the P density from init map, [g P / kg soil]

    !> 2.0 calculate the initial slow P, occluded P, primary P pool, and p_labile
    !! based on the assumption that soil P pool are linearly corrected with soil weight
    !! the depth correction of initial P pool distributions are fitted against the observed pattern of a young cambisol soil
    !! with basalt as parent material (BBR, Lang et al. 2017) that,
    !! the proportion of: labile+sorbed P decreases with increasing depth, and primary P increases with increasing depth (Yu et al., 2020)
    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        ! primary P, corrected with depth, [g P]
        hlp1(ic) = MIN(k2_p_pool_depth_corr * soil_lay_depth_center_sl(ic,is) ** k4_p_pool_depth_corr, &
                      (k1_p_pool_depth_corr - k_p_pool_depth_corr_min_occlud)) / k1_p_pool_depth_corr
        hlp2(ic) = k3_p_pool_depth_corr ** (k2_p_pool_depth_corr * soil_lay_depth_center_sl(ic,is)) / k1_p_pool_depth_corr
        ! occluded P, corrected with depth
        hlp3(ic) = MAX((1.0_wp - hlp1(ic) - hlp2(ic)), k_p_pool_depth_corr_min_occlud / k1_p_pool_depth_corr)

        IF ((hlp1(ic) + hlp2(ic) + hlp3(ic)) > 1._wp) THEN
          hlp1(ic) = hlp1(ic) / (hlp1(ic) + hlp2(ic) + hlp3(ic))
          hlp2(ic) = hlp2(ic) / (hlp1(ic) + hlp2(ic) + hlp3(ic))
          hlp3(ic) = hlp3(ic) / (hlp1(ic) + hlp2(ic) + hlp3(ic))
        END IF
        IF (soil_p_layer(ic) > is) THEN
          po4_assoc_slow(ic,is) = bulk_dens_sl(ic,is) * hlp4(ic) * (soil_p_slow(ic,is) + soil_p_labile(ic,is)) &
            &                     / rhlp1(ic) / molar_mass_P * (1._wp - p_labile_slow_ratio(ic,is))
          po4_occluded(ic,is)   = bulk_dens_sl(ic,is) * hlp4(ic) * soil_p_occluded(ic,is) &
            &                     / rhlp1(ic) / molar_mass_P
          po4_primary(ic,is)    = bulk_dens_sl(ic,is) * hlp4(ic) * soil_p_primary(ic,is) &
            &                     / rhlp1(ic) / molar_mass_P
          arr_hlp1(ic,is)       = bulk_dens_sl(ic,is) * hlp4(ic) * (soil_p_slow(ic,is) + soil_p_labile(ic,is)) &
            &                     / rhlp1(ic) / molar_mass_P * p_labile_slow_ratio(ic,is)
        ELSE
          po4_assoc_slow(ic,is) = bulk_dens_sl(ic,is) * hlp4(ic) * hlp2(ic) / molar_mass_P &
            &                     * (1._wp - p_labile_slow_ratio(ic,is))
          po4_occluded(ic,is)   = bulk_dens_sl(ic,is) * hlp4(ic) * hlp3(ic) / molar_mass_P
          po4_primary(ic,is)    = bulk_dens_sl(ic,is) * hlp4(ic) * hlp1(ic) / molar_mass_P
          arr_hlp1(ic,is)       = bulk_dens_sl(ic,is) * hlp4(ic) * hlp2(ic) / molar_mass_P &
            &                     * p_labile_slow_ratio(ic,is)
        END IF
      END DO
    END DO

    !> 2.1 calculate the initial solute P and fast P,
    !! based on the assumption that labile P is the sum of both,
    !! and they follow the Langmuir equilibrium
    !! po4_assoc_fast = (p_labile_sl + Smax + K - sqrt((p_labile_sl + Smax + K)**2._wp - 4*p_labile_sl*Smax)) / 2._wp
    DO ic = 1,nc
      DO is = 1, NINT(num_sl_above_bedrock(ic))
        IF (flag_sb_double_langmuir) THEN
          po4_assoc_fast(ic, is) = calc_fast_po4(arr_hlp1(ic, is), qmax_fast_po4(ic, is), km_fast_po4(ic, is))
          po4_solute(ic, is)     = arr_hlp1(ic, is) - po4_assoc_fast(ic, is)
        ELSE
          po4_assoc_fast(ic, is) = calc_fast_po4(arr_hlp1(ic, is), qmax_po4(ic, is), km_po4(ic, is))
          po4_solute(ic, is)     = arr_hlp1(ic, is) - po4_assoc_fast(ic, is)
        END IF
      END DO
    END DO
  END SUBROUTINE init_soil_p_pools

#ifndef __QUINCY_STANDALONE__
  ! ====================================================================================================== !
  !
  !> Read soil biogeochemistry state from file (after init/restart)
  !> (variables should be in sync with those listed in ./scripts/quincy_sb_and_veg_state_from_restart_file.sh)
  !
  ! ======================================================================================================= !
  SUBROUTINE sb_read_states(tile)
    USE mo_jsb_io_netcdf,           ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_util,                    ONLY: read_3D_var
    USE mo_lnd_bgcm_store,          ONLY: read_bgcm_from_file
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER :: model            !< the model
    TYPE(t_lnd_bgcm_store), POINTER :: bgcm_store       !< the bgcm store of this tile
    TYPE(t_input_file)              :: input_file       !< file with to be read sb state variables
    TYPE(t_jsb_vgrid),      POINTER :: vgrid_soil_w     !< vertical grid for sl variables
    INTEGER                         :: nsoil_w          !< number of soil layers as used/defined by the SB_ process
    CHARACTER(len=*), PARAMETER :: routine = modname//':sb_read_states'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Real3D_onDomain :: microbial_cue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_nue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: microbial_pue_eff_tmic_mavg
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_c_mavg
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_n_mavg
    dsl4jsb_Real3D_onDomain :: enzyme_frac_poly_p_mavg
    dsl4jsb_Real3D_onDomain :: nh4_assoc
    dsl4jsb_Real3D_onDomain :: po4_assoc_slow
    dsl4jsb_Real3D_onDomain :: po4_primary
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    vgrid_soil_w => Get_vgrid('soil_depth_water')
    nsoil_w = vgrid_soil_w%n_levels
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_var3D_onDomain(SB_, microbial_cue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_, microbial_nue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_, microbial_pue_eff_tmic_mavg)
    dsl4jsb_Get_var3D_onDomain(SB_, nh4_assoc)
    dsl4jsb_Get_var3D_onDomain(SB_, po4_assoc_slow)
    dsl4jsb_Get_var3D_onDomain(SB_, po4_primary)
    IF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == 'jsm') THEN
      dsl4jsb_Get_var3D_onDomain(SB_, enzyme_frac_poly_c_mavg)
      dsl4jsb_Get_var3D_onDomain(SB_, enzyme_frac_poly_n_mavg)
      dsl4jsb_Get_var3D_onDomain(SB_, enzyme_frac_poly_p_mavg)
    END IF
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on()) CALL message(TRIM(routine), ' for tile '//TRIM(tile%name))

    input_file = jsb_netcdf_open_input(TRIM(dsl4jsb_Config(SB_)%ic_filename), model%grid_id)

    IF (ASSOCIATED(tile%bgcm_store)) THEN
      bgcm_store => tile%bgcm_store
      CALL read_bgcm_from_file(bgcm_store, SB_BGCM_POOL_ID, "sb", tile%name, routine, input_file, 'soillev', nsoil_w)
    END IF

    CALL read_3D_var(input_file, microbial_cue_eff_tmic_mavg, 'sb_microbial_cue_eff_tmic_mavg_'//tile%name, 'soillev', nsoil_w)
    CALL read_3D_var(input_file, microbial_nue_eff_tmic_mavg, 'sb_microbial_nue_eff_tmic_mavg_'//tile%name, 'soillev', nsoil_w)
    CALL read_3D_var(input_file, microbial_pue_eff_tmic_mavg, 'sb_microbial_pue_eff_tmic_mavg_'//tile%name, 'soillev', nsoil_w)

    CALL read_3D_var(input_file, nh4_assoc, 'sb_nh4_assoc_'//tile%name, 'soillev', nsoil_w)
    CALL read_3D_var(input_file, po4_assoc_slow, 'sb_po4_primary_'//tile%name, 'soillev', nsoil_w)
    CALL read_3D_var(input_file, po4_primary, 'sb_po4_assoc_slow_'//tile%name, 'soillev', nsoil_w)

    IF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == 'jsm') THEN
      CALL read_3D_var(input_file, enzyme_frac_poly_c_mavg, 'sb_enzyme_frac_poly_c_mavg_'//tile%name, 'soillev', nsoil_w)
      CALL read_3D_var(input_file, enzyme_frac_poly_n_mavg, 'sb_enzyme_frac_poly_n_mavg_'//tile%name, 'soillev', nsoil_w)
      CALL read_3D_var(input_file, enzyme_frac_poly_p_mavg, 'sb_enzyme_frac_poly_p_mavg_'//tile%name, 'soillev', nsoil_w)
    END IF

    CALL input_file%Close()

  END SUBROUTINE sb_read_states

  ! ====================================================================================================== !
  !
  !> Provide n and p deposition depending on the deposition configuration and the time of the year - for IQ
  !
  SUBROUTINE provide_n_and_p_deposition(model_id, current_datetime, dtime)

    USE mo_jsb_time_iface,         ONLY: t_datetime
    USE mo_jsb_time,               ONLY: is_newyear, get_year, get_month, is_time_experiment_start
    USE mo_io_units,               ONLY: filename_max
    USE mo_jsb_parallel,           ONLY: Get_omp_thread
    USE mo_jsb_io_netcdf,          ONLY: t_input_file, jsb_netcdf_open_input
    USE mo_jsb_physical_constants, ONLY: molar_mass_N, molar_mass_P
    USE mo_isotope_util,           ONLY: calc_mixing_ratio_N15N14
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(SB_)
    dsl4jsb_Use_config(SB_)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id
    TYPE(t_datetime), POINTER, INTENT(in) :: current_datetime
    REAL(wp),                  INTENT(in) :: dtime
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),           POINTER :: model
    TYPE(t_jsb_grid),            POINTER :: hgrid
    CLASS(t_jsb_tile_abstract),  POINTER :: tile

    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(SB_)

    dsl4jsb_Real2D_onDomain :: nhx_deposition
    dsl4jsb_Real2D_onDomain :: noy_deposition
    dsl4jsb_Real2D_onDomain :: nhx_n15_deposition
    dsl4jsb_Real2D_onDomain :: noy_n15_deposition
    dsl4jsb_Real2D_onDomain :: p_deposition

    dsl4jsb_Real3D_onDomain :: nhx_deposition_monthly
    dsl4jsb_Real3D_onDomain :: noy_deposition_monthly
    dsl4jsb_Real3D_onDomain :: p_deposition_monthly

    REAL(wp),       POINTER :: ptr_3D(:,:,:)  ! tmp pointer

    INTEGER :: current_year, current_month, no_omp_thread, nproma, nblks, i_month
    LOGICAL :: is_experiment_start, read_deposition_data
    TYPE(t_input_file)          :: input_file
    CHARACTER(len=filename_max) :: filename_deposition_data
    CHARACTER(len=*), PARAMETER :: routine = modname//':provide_n_and_p_deposition'
    ! -------------------------------------------------------------------------------------------------- !
    model => Get_model(model_id)
    CALL model%Get_top_tile(tile)
    no_omp_thread = Get_omp_thread()
    hgrid => Get_grid(model%grid_id)

    nproma = hgrid%nproma
    nblks  = hgrid%nblks

    is_experiment_start = is_time_experiment_start(current_datetime)
    ! -------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_config(SB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_var2D_onDomain(A2L_, nhx_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_, noy_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_, nhx_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_, noy_n15_deposition)
    dsl4jsb_Get_var2D_onDomain(A2L_, p_deposition)
    dsl4jsb_Get_var3D_onDomain(SB_, nhx_deposition_monthly)
    dsl4jsb_Get_var3D_onDomain(SB_, noy_deposition_monthly)
    dsl4jsb_Get_var3D_onDomain(SB_, p_deposition_monthly)
    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on()) CALL message( TRIM(routine), 'Starting routine')

    ! In case that constants are to be used for the n and p deposition they are initialised here
    IF (dsl4jsb_Config(SB_)%deposition_scheme == CONST_DEP) THEN
      IF (is_experiment_start) THEN
        ! set constants ...
        ! SZ: "background N deposition is 2 kg N / ha / yr = 6.3419e-12 kg / m2 / s (2/365/24/3600/10000)
        !      the background N deposition is the sum of NHx and NOy, therefore "/ 2._wp"
        ! ...  and convert from kg/m2/s -> mumol/m2/s
        ! -- 1 kg/m2/s = 1 * 1000 / molar_mass_N / 1e-06 [mu mol/m2/s]
        nhx_deposition(:,:) = 6.3419E-12_wp / 2._wp * 1000._wp / molar_mass_N / 1.e-6_wp
        noy_deposition(:,:) = 6.3419E-12_wp / 2._wp * 1000._wp / molar_mass_N / 1.e-6_wp
        ! SZ: "Default P deposition is assumed to be stoichmetrically balanced (which would be higher than typical background values.)
        !      Thus, total N deposition to total P deposition = 14 g N/ g P -> 14 / molar_mass_N * molar_mass_P"
        !      (Sterner & Elser, Ecological Stoichmetry, 2002, Princeton University Press)
        p_deposition(:,:) =  6.3419E-12 / (14.0_wp / molar_mass_N * molar_mass_P) * 1000._wp / molar_mass_N / 1.e-6_wp
        ! SZ: NHX_deposition_N15 = NOY_deposition_N15 = NOY_deposition f(delta_n15=0.0)
        nhx_n15_deposition(:,:) = nhx_deposition(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
        noy_n15_deposition(:,:) = noy_deposition(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
      END IF

    ELSE
      ! Check if data needs to be read
      IF (is_experiment_start &
          & .OR. (is_newyear(current_datetime, dtime) .AND. (dsl4jsb_Config(SB_)%deposition_scheme == TRANS_DEP))) THEN
        read_deposition_data = .TRUE.
      ELSE
        read_deposition_data = .FALSE.
      END IF

      IF (read_deposition_data) THEN
        IF (.NOT. is_experiment_start) THEN
          ! In this case we search a file ending on the current year
          current_year  = get_year(current_datetime)

          !>
          !> Assertion: routine currently expects filenames with 4 digits
          !>
          IF (( current_year > 9999) .OR. (current_year < 1000)) THEN
            CALL finish(TRIM(routine), 'Violation of assertion: this routine currently expects filenames with 4 digits.')
          END IF
          WRITE (filename_deposition_data,'(a,a,I4.4,a)') &
            & TRIM(dsl4jsb_Config(SB_)%deposition_filename_prefix), '_', current_year, ".nc"
        ELSE
          ! else we expect a file name without a year
          WRITE (filename_deposition_data,'(a,a)') TRIM(dsl4jsb_Config(SB_)%deposition_filename_prefix), ".nc"
        END IF

        input_file = jsb_netcdf_open_input(TRIM(filename_deposition_data), model%grid_id)
        ptr_3D => input_file%Read_2d_time(variable_name=TRIM("NHx_deposition")) !, start_time_step=1, end_time_step=12)
        DO i_month = 1,12
          nhx_deposition_monthly(:,i_month,:) = ptr_3D(:,:,i_month)
        END DO
        ptr_3D => input_file%Read_2d_time(variable_name=TRIM("NOy_deposition")) !, start_time_step=1, end_time_step=12)
        DO i_month = 1,12
          noy_deposition_monthly(:,i_month,:) = ptr_3D(:,:,i_month)
        END DO
        ptr_3D => input_file%Read_2d_time(variable_name=TRIM("pdep")) !, start_time_step=1, end_time_step=12)
        DO i_month = 1,12
          p_deposition_monthly(:,i_month,:) = ptr_3D(:,:,i_month)
        END DO

        CALL input_file%Close()
      END IF

      ! Assign the deposition values of this month to the A2L variables
      current_month = get_month(current_datetime)

      ! ... convert from kg/m2/s -> mumol/m2/s
      ! -- 1 kg/m2/s = 1 * 1000 / molar_mass_N / 1e-06 [mu mol/m2/s]
      nhx_deposition(:,:) = nhx_deposition_monthly(:,current_month,:) * 1000._wp / molar_mass_N / 1.e-6_wp
      noy_deposition(:,:) = noy_deposition_monthly(:,current_month,:) * 1000._wp / molar_mass_N / 1.e-6_wp
      p_deposition  (:,:) = p_deposition_monthly  (:,current_month,:) * 1000._wp / molar_mass_P / 1.e-6_wp

      ! SZ: NHX_deposition_N15 = NOY_deposition_N15 = NOY_deposition f(delta_n15=0.0)
      nhx_n15_deposition(:,:) = nhx_deposition(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
      noy_n15_deposition(:,:) = noy_deposition(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp))
    END IF

    IF (debug_on()) CALL message(TRIM(routine), 'Finishing routine')

  END SUBROUTINE provide_n_and_p_deposition
#endif
#endif
END MODULE mo_sb_init
