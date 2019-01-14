!>
!! @brief Contains the routines for reading parameters from parameter_slm_run.list and replace parameter-values as defined in mo_PROC_constants
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (2016-10-04)
!!
MODULE mo_shared_set_parameters

  USE mo_kind,                ONLY: wp, eps8
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_jsb_io,              ONLY: missval  ! = -9.E+15_wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: parameters_slm_from_file                        ! it's an object
  PUBLIC :: set_parameter_lctlib_values_with_list_values    ! it's a subroutine

  !-----------------------------------------------------------------------------------------------------
  !> SLM process parameter values from file parameter_slm_run.list - replacing the values from mo_PROC_constants
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_parameters_slm_from_file
    !LOGICAL            :: active   ! de/activate this TYPE (this is non-functional for now 20161004)
    CHARACTER(LEN=30)  :: filename = 'parameter_slm_run.list'                 ! filename 
    ! ASSIMILATION_CTL
    REAL(wp)           :: jmax2n
    REAL(wp)           :: vcmax2n
    REAL(wp)           :: pepc2n
    REAL(wp)           :: chl2n
    REAL(wp)           :: k1_fn_struc
    REAL(wp)           :: E0kc
    REAL(wp)           :: E1kc
    REAL(wp)           :: E0ko
    REAL(wp)           :: E1ko
    REAL(wp)           :: E0pcp
    REAL(wp)           :: E1pcp
    REAL(wp)           :: E0v
    REAL(wp)           :: E1v
    REAL(wp)           :: t_jmax_offset
    REAL(wp)           :: t_jmax_slope
    REAL(wp)           :: t_jmax_opt_min
    REAL(wp)           :: t_jmax_opt_max
    REAL(wp)           :: Tref_pepc
    REAL(wp)           :: Tbase_pepc
    REAL(wp)           :: alpha_i
    REAL(wp)           :: pO2
    REAL(wp)           :: CiCa_default_C3
    REAL(wp)           :: CiCa_default_C4
    REAL(wp)           :: ps_it_max
    REAL(wp)           :: ci_max
    REAL(wp)           :: Dwv2co2_air
    REAL(wp)           :: Dwv2co2_turb
    REAL(wp)           :: discr_ps_a_C13
    REAL(wp)           :: discr_ps_b_C13
    REAL(wp)           :: discr_ps_c_C13
    REAL(wp)           :: discr_ps_a_C14
    REAL(wp)           :: discr_ps_b_C14
    REAL(wp)           :: discr_ps_c_C14
    REAL(wp)           :: discr_ps_phi
    REAL(wp)           :: ka
    REAL(wp)           :: kn
    REAL(wp)           :: soa_min_temp
    REAL(wp)           :: soa_max_temp
    ! PHENOLOGY_CTL
    REAL(wp)           :: gdd_t_air_threshold
    ! RADIATION_CTL
    REAL(wp)           :: kbl0_vis
    REAL(wp)           :: kdf0_vis
    REAL(wp)           :: kbl0_nir
    REAL(wp)           :: kdf0_nir
    REAL(wp)           :: rho2sbeta
    REAL(wp)           :: k_csf
    REAL(wp)           :: min_cos_zenith_angle_rad
    REAL(wp)           :: def_alb_vis_soil
    REAL(wp)           :: def_alb_nir_soil
    REAL(wp)           :: rad2ppfd
    REAL(wp)           :: rfr_ratio_toc
    REAL(wp)           :: k_r2fr_chl
    ! SHARED_CTL
    REAL(wp)           :: mavg_period_weekly
    REAL(wp)           :: mavg_period_monthly
    REAL(wp)           :: mavg_period_yearly
    REAL(wp)           :: mavg_period_tsoa
    REAL(wp)           :: mavg_period_tuptake
    REAL(wp)           :: mavg_period_tlabile
    REAL(wp)           :: mavg_period_tphen
    REAL(wp)           :: mavg_period_tenzyme
    REAL(wp)           :: mavg_period_tfrac
    REAL(wp)           :: mavg_period_tcnl
    REAL(wp)           :: mavg_period_talloc
    REAL(wp)           :: mavg_period_talloc_dynamic
    REAL(wp)           :: mavg_period_tacclim
    REAL(wp)           :: mavg_period_tmic
    REAL(wp)           :: mavg_period_tgrowth
    REAL(wp)           :: mavg_period_tvegdyn
    ! SSE_CTL
    REAL(wp)           :: soil_heat_cap
    REAL(wp)           :: soil_therm_cond
    REAL(wp)           :: InterceptionEfficiency
    REAL(wp)           :: w_skin_max
    REAL(wp)           :: drain_frac
    REAL(wp)           :: fdrain_runoff
    REAL(wp)           :: wpot_pwp
    REAL(wp)           :: wpot_fc
    REAL(wp)           :: wpot_min
    REAL(wp)           :: k_pwp_a
    REAL(wp)           :: k_pwp_c
    REAL(wp)           :: k_pwp_s
    REAL(wp)           :: k_pwp_sc
    REAL(wp)           :: k_pwp_at
    REAL(wp)           :: k_pwp_bt
    REAL(wp)           :: k_fc_a
    REAL(wp)           :: k_fc_c
    REAL(wp)           :: k_fc_s
    REAL(wp)           :: k_fc_sc
    REAL(wp)           :: k_fc_at
    REAL(wp)           :: k_fc_bt
    REAL(wp)           :: k_fc_ct
    REAL(wp)           :: k_sat_a
    REAL(wp)           :: k_sat_c
    REAL(wp)           :: k_sat_s
    REAL(wp)           :: k_sat_sc
    REAL(wp)           :: k_sat_at
    REAL(wp)           :: k_sat_bt
    REAL(wp)           :: k_sat_ct
    REAL(wp)           :: k_sat_dt
    REAL(wp)           :: kdiff_sat_max
    ! VEGETATION_CTL
    REAL(wp)           :: fresp_growth
    REAL(wp)           :: fmaint_rate_base
    REAL(wp)           :: fmaint_rate_wood
    REAL(wp)           :: t_maint_k1
    REAL(wp)           :: t_maint_k2
    REAL(wp)           :: t_maint_k3
    REAL(wp)           :: f_resp_acclim
    REAL(wp)           :: t_acclim_zero
    REAL(wp)           :: k0_fn_pepc_C3
    REAL(wp)           :: k0_fn_pepc_C4
    REAL(wp)           :: jmax2vcmax_C3
    REAL(wp)           :: jmax2vcmax_C4
    REAL(wp)           :: tau_labile
    REAL(wp)           :: lambda_labile_temp
    REAL(wp)           :: k_labile_temp
    REAL(wp)           :: lambda_labile_theta
    REAL(wp)           :: k_labile_theta
    REAL(wp)           :: lambda_phiphen
    REAL(wp)           :: k_phiphen
    REAL(wp)           :: lambda_phimaint
    REAL(wp)           :: k_phimaint
    REAL(wp)           :: lambda_phimaint_nut
    REAL(wp)           :: k_phimaint_nut
    REAL(wp)           :: lambda_phistore
    REAL(wp)           :: k_phistore
    REAL(wp)           :: knut_labile
    REAL(wp)           :: lambda_sinklim_ps
    REAL(wp)           :: k_sinklim_ps
    REAL(wp)           :: k_cnp_sinklim_ps
    REAL(wp)           :: beta_sinklim_ps_min
    REAL(wp)           :: fstore_leaf_max
    REAL(wp)           :: fstore_fine_root_max
    REAL(wp)           :: fstore_sap_wood_max
    REAL(wp)           :: k_phi_interact
    REAL(wp)           :: resorp_fract_leaf
    REAL(wp)           :: resorp_fract_wood
    REAL(wp)           :: root2leaf_cn
    REAL(wp)           :: wood2leaf_cn
    REAL(wp)           :: root2leaf_np
    REAL(wp)           :: wood2leaf_np
    REAL(wp)           :: tau_nutrient_recycling
    REAL(wp)           :: kappa_f_demand_nc
    REAL(wp)           :: kappa_f_demand_pn
    REAL(wp)           :: k_f_demand
    REAL(wp)           :: cost_n_fixation_max
    REAL(wp)           :: cost_n_fixation_min
    REAL(wp)           :: k1_temp_bnf
    REAL(wp)           :: k2_temp_bnf
    REAL(wp)           :: k3_temp_bnf
    REAL(wp)           :: k4_temp_bnf
    REAL(wp)           :: transform_cost_nh4
    REAL(wp)           :: transform_cost_no3
    REAL(wp)           :: transform_cost_org
    REAL(wp)           :: f_nacq_nh4
    REAL(wp)           :: f_nacq_no3
    REAL(wp)           :: vmax_bnf
    REAL(wp)           :: k_bnf
    REAL(wp)           :: offset_bnf
    REAL(wp)           :: eta_nfixation
    REAL(wp)           :: k0_fn_chl_C3
    REAL(wp)           :: k0_fn_chl_C4
    REAL(wp)           :: k1_fn_chl_C3
    REAL(wp)           :: k1_fn_chl_C4
    REAL(wp)           :: kfn_chl
    REAL(wp)           :: min_lai
    REAL(wp)           :: target_lai_max
    REAL(wp)           :: max_leaf_shedding_rate
    REAL(wp)           :: k_leafon_canopy
    REAL(wp)           :: background_mort_rate_tree
    REAL(wp)           :: background_mort_rate_grass
    REAL(wp)           :: lambda_mort_light
    REAL(wp)           :: k_mort_light
    REAL(wp)           :: lambda_est_temp
    REAL(wp)           :: k_est_temp
    REAL(wp)           :: lambda_est_moist
    REAL(wp)           :: k_est_moist
    REAL(wp)           :: wood_extraction_rate
    REAL(wp)           :: k1_mort_greff
    REAL(wp)           :: k2_mort_greff
    REAL(wp)           :: k3_mort_greff
    REAL(wp)           :: min_greff
    REAL(wp)           :: max_crown_area
    REAL(wp)           :: min_diameter
    REAL(wp)           :: k_crown_area
    REAL(wp)           :: k_rp
    REAL(wp)           :: k_fpc
    REAL(wp)           :: fpc_max
    REAL(wp)           :: carbon_per_dryweight_leaf
    REAL(wp)           :: min_height
    REAL(wp)           :: sm2lm_grass
    REAL(wp)           :: lambda_fruit_alloc
    REAL(wp)           :: k1_fruit_alloc
    REAL(wp)           :: k3_fruit_alloc
    REAL(wp)           :: k4_fruit_alloc
    REAL(wp)           :: w_root_lim_max
    REAL(wp)           :: leaf2root_min_ratio
    REAL(wp)           :: leaf2root_max_ratio
    REAL(wp)           :: delta_n_fraction
    REAL(wp)           :: delta_n_leaf
    ! SB_CTL
    REAL(wp)           :: tau_soluable_litter
    REAL(wp)           :: tau_metabolic_litter
    REAL(wp)           :: tau_structural_litter
    REAL(wp)           :: tau_woody_litter
    REAL(wp)           :: tau_fast
    REAL(wp)           :: tau_slow
    REAL(wp)           :: tau_biomineralisation
    REAL(wp)           :: t_ref_tresponse
    REAL(wp)           :: ea_depolymerisation
    REAL(wp)           :: ea_mic_uptake
    REAL(wp)           :: ea_sorption
    REAL(wp)           :: ea_desorption
    REAL(wp)           :: ea_hsc
    REAL(wp)           :: ea_decomposition
    REAL(wp)           :: ed_decomposition
    REAL(wp)           :: t_ref_decomposition
    REAL(wp)           :: k1_afps
    REAL(wp)           :: kmr_depolymerisation
    REAL(wp)           :: kmr_mic_uptake
    REAL(wp)           :: kmr_hsc
    REAL(wp)           :: kmr_sorption
    REAL(wp)           :: kmr_decomposition
    REAL(wp)           :: k_diff_org
    REAL(wp)           :: rho_bulk_org
    REAL(wp)           :: carbon_per_dryweight_SOM
    REAL(wp)           :: qmax_po4_fine_particle
    REAL(wp)           :: k_adsorpt_dom
    REAL(wp)           :: k_desorpt_dom
    REAL(wp)           :: k_adsorpt_det
    REAL(wp)           :: k_desorpt_det
    REAL(wp)           :: k_adsorpt_nh4
    REAL(wp)           :: k_desorpt_nh4
    REAL(wp)           :: k_adsorpt_po4
    REAL(wp)           :: k_desorpt_po4
    REAL(wp)           :: k_adsorpt_f_po4
    REAL(wp)           :: k_desorp_f_po4
    REAL(wp)           :: k_desorp_s_po4
    REAL(wp)           :: k_adsorpt_s_po4
    REAL(wp)           :: microbial_cn
    REAL(wp)           :: microbial_np
    REAL(wp)           :: SOM_np
    REAL(wp)           :: SOM_cn
    REAL(wp)           :: microbial_nue
    REAL(wp)           :: microbial_pue
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: microbial_cue_max
    REAL(wp)           :: microbial_cue_min
    REAL(wp)           :: microbial_cn_max
    REAL(wp)           :: microbial_cn_min
    REAL(wp)           :: f_mic_cn
    REAL(wp)           :: frac_woody2dec_litter
    REAL(wp)           :: frac_litter2fast_som
    REAL(wp)           :: frac_fast2slow_som
    REAL(wp)           :: frac_slow2fast_som
    REAL(wp)           :: fc_woody2polymeric
    REAL(wp)           :: fc_soluable2dom
    REAL(wp)           :: k_mort_microbial
    REAL(wp)           :: min_mic_biomass
    REAL(wp)           :: vmax_mic_uptake
    REAL(wp)           :: km_mic_uptake
    REAL(wp)           :: km_depolymerisation
    REAL(wp)           :: vmax_depolymerisation_residue
    REAL(wp)           :: vmax_depolymerisation_litter
    REAL(wp)           :: f_soluable_residue
    REAL(wp)           :: f_soluable_leaf
    REAL(wp)           :: f_soluable_fine_root
    REAL(wp)           :: f_soluable_fruit
    REAL(wp)           :: f_soluable_seed_bed
    REAL(wp)           :: fc_soluable_max
    REAL(wp)           :: k_fc_soluable
    REAL(wp)           :: k_fn_soluable
    REAL(wp)           :: k_fp_soluable
    REAL(wp)           :: lc_leaf_max
    REAL(wp)           :: lc_leaf2sla
    REAL(wp)           :: lc_fine_root
    REAL(wp)           :: lc_coarse_root
    REAL(wp)           :: lc_woody_litter
    REAL(wp)           :: lc_fruit
    REAL(wp)           :: lc_seed_bed
    REAL(wp)           :: mycorrhiza_cue
    REAL(wp)           :: mycorrhiza_cn
    REAL(wp)           :: mycorrhiza_np
    REAL(wp)           :: f_maxresp_myc
    REAL(wp)           :: f_surface_myc2fineroot
    REAL(wp)           :: km1_up_nh4
    REAL(wp)           :: km2_up_nh4
    REAL(wp)           :: km1_up_no3
    REAL(wp)           :: km2_up_no3
    REAL(wp)           :: km1_up_po4
    REAL(wp)           :: km2_up_po4
    REAL(wp)           :: vmax_mic_uptake_nh4
    REAL(wp)           :: km_mic_uptake_nh4
    REAL(wp)           :: vmax_mic_uptake_no3
    REAL(wp)           :: km_mic_uptake_no3
    REAL(wp)           :: f_mic_carrier_n
    REAL(wp)           :: vmax_plant_uptake_nh4
    REAL(wp)           :: km_plant_uptake_nh4
    REAL(wp)           :: vmax_plant_uptake_no3
    REAL(wp)           :: km_plant_uptake_no3
    REAL(wp)           :: f_plant_carrier_n
    REAL(wp)           :: vmax_mic_uptake_po4
    REAL(wp)           :: km_mic_uptake_po4
    REAL(wp)           :: vmax_plant_uptake_po4
    REAL(wp)           :: km_plant_uptake_po4
    REAL(wp)           :: k_mic_uptake_nh4
    REAL(wp)           :: k_plant_uptake_nh4
    REAL(wp)           :: k_mic_uptake_no3
    REAL(wp)           :: k_plant_uptake_no3
    REAL(wp)           :: k_mic_uptake_po4
    REAL(wp)           :: k_plant_uptake_po4
    REAL(wp)           :: f_mic_carrier_po4
    REAL(wp)           :: f_plant_carrier_po4
    REAL(wp)           :: qmax_po4_mineral
    REAL(wp)           :: km_adsorpt_mineral_po4
    REAL(wp)           :: km_adsorpt_OM_po4
    REAL(wp)           :: ea_nitrification
    REAL(wp)           :: ed_nitrification
    REAL(wp)           :: t_opt_nitrification
    REAL(wp)           :: f_nit_noy
    REAL(wp)           :: f_nit_n2o
    REAL(wp)           :: ea_denitrification
    REAL(wp)           :: vmax_nitrification
    REAL(wp)           :: vmax_denitrification
    REAL(wp)           :: km_denitrification_c
    REAL(wp)           :: km_denitrification_no3
    REAL(wp)           :: f_denit_noy
    REAL(wp)           :: f_denit_n2o
    REAL(wp)           :: ea_gasdiffusion
    REAL(wp)           :: lambda_anvf
    REAL(wp)           :: k_anvf
    REAL(wp)           :: floss_nmin
    REAL(wp)           :: fnit_nmin
    REAL(wp)           :: fleach_nh4
    REAL(wp)           :: fleach_po4
    REAL(wp)           :: soil_n_threshold_bnf
    REAL(wp)           :: standard_bnf_max_rate
    REAL(wp)           :: eta_mic_uptake_nh4
    REAL(wp)           :: eta_mic_uptake_no3
    REAL(wp)           :: eta_plant_uptake_nh4
    REAL(wp)           :: eta_plant_uptake_no3
    REAL(wp)           :: eta_nitrification
    REAL(wp)           :: eta_nitrate_production
    REAL(wp)           :: eta_denitrification
    REAL(wp)           :: eta_ammonification
    REAL(wp)           :: eta_volatilisation_nh3
    REAL(wp)           :: k_weath_mineral
    REAL(wp)           :: eq_sorp_mineral_po4
    REAL(wp)           :: k_occlude_po4
    REAL(wp)           :: vmax_biochem_po4
    REAL(wp)           :: km_mic_biochem_po4
    REAL(wp)           :: km_dip_biochem_po4
    REAL(wp)           :: km_rootc_biochem_po4
    REAL(wp)           :: km_rootc_weath_po4
    REAL(wp)           :: biomin_thres_PC
    REAL(wp)           :: fc_p_res2dom
    REAL(wp)           :: fc_n_res2dom
  CONTAINS
    PROCEDURE :: read_param => read_parameter_slm_list
    PROCEDURE :: set_param  => set_parameter_slm_values_with_list_values
  END TYPE t_parameters_slm_from_file

  TYPE(t_parameters_slm_from_file), SAVE :: parameters_slm_from_file


  !-----------------------------------------------------------------------------------------------------
  !> lctlib parameter values from file parameter_slm_run.list - replacing the values from lctlib_slm_*.def file
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_parameters_lctlib_from_file
    !LOGICAL            :: active   
    CHARACTER(LEN=30)   :: filename = 'parameter_slm_run.list'                 ! filename 

    REAL(wp)            :: sla
    REAL(wp)            :: sigma_vis
    REAL(wp)            :: sigma_nir
    REAL(wp)            :: omega_clumping
    REAL(wp)            :: crown_shape_factor
    REAL(wp)            :: cn_leaf
    REAL(wp)            :: cn_leaf_min
    REAL(wp)            :: cn_leaf_max
    REAL(wp)            :: np_leaf
    REAL(wp)            :: np_leaf_min
    REAL(wp)            :: np_leaf_max
    REAL(wp)            :: k0_fn_struc
    REAL(wp)            :: fn_oth_min
    REAL(wp)            :: t_jmax_opt
    REAL(wp)            :: t_jmax_omega
    REAL(wp)            :: g0
    REAL(wp)            :: g1
    REAL(wp)            :: gmin
    ! turnover times
    REAL(wp)            :: tau_leaf
    REAL(wp)            :: tau_fine_root
    REAL(wp)            :: tau_coarse_root
    REAL(wp)            :: tau_branch
    REAL(wp)            :: tau_sap_wood
    REAL(wp)            :: tau_fruit
    REAL(wp)            :: tau_seed_litter
    REAL(wp)            :: tau_seed_est
    REAL(wp)            :: tau_mycorrhiza
    ! N uptake parameters
    REAL(wp)            :: vmax_uptake_n
    REAL(wp)            :: vmax_uptake_p
    REAL(wp)            :: bnf_base
    ! Vegetation dynamics parameters
    REAL(wp)            :: lambda_est_light
    REAL(wp)            :: k_est_light
    REAL(wp)            :: seed_size
    ! Phenology parameters
    REAL(wp)            :: beta_soil_flush
    REAL(wp)            :: beta_soil_senescence
    REAL(wp)            :: gdd_req_max
    REAL(wp)            :: k_gdd_dormance
    REAL(wp)            :: t_air_senescence
    REAL(wp)            :: min_leaf_age
    ! Allocation parameters
    REAL(wp)            :: frac_sapwood_branch
    REAL(wp)            :: wood_density
    REAL(wp)            :: k_latosa
    REAL(wp)            :: k_crtos
    REAL(wp)            :: k_rtos
    REAL(wp)            :: k2_fruit_alloc
    REAL(wp)            :: allom_k1
    REAL(wp)            :: allom_k2
    REAL(wp)            :: phi_leaf_min
    REAL(wp)            :: k_root
    REAL(wp)            :: k_sapwood
    REAL(wp)            :: c0_allom
    REAL(wp)            :: fstore_target
    ! Soil
    REAL(wp)            :: k_root_dist

  END TYPE t_parameters_lctlib_from_file

  TYPE(t_parameters_lctlib_from_file) :: param_pft1_from_file,  &
                                         param_pft2_from_file,  &
                                         param_pft3_from_file,  &
                                         param_pft4_from_file,  &
                                         param_pft5_from_file,  &
                                         param_pft6_from_file,  &
                                         param_pft7_from_file,  &
                                         param_pft8_from_file


  CHARACTER(len=*), PARAMETER :: modname = 'mo_shared_set_parameters'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> read SLM parameter values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_parameter_slm_list(this)

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_parameters_slm_from_file), INTENT(inout) :: this
    ! ---------------------------
    ! 0.2 Local
    ! ASSIMILATION_CTL
    REAL(wp)           :: jmax2n
    REAL(wp)           :: vcmax2n
    REAL(wp)           :: pepc2n
    REAL(wp)           :: chl2n
    REAL(wp)           :: k1_fn_struc
    REAL(wp)           :: E0kc
    REAL(wp)           :: E1kc
    REAL(wp)           :: E0ko
    REAL(wp)           :: E1ko
    REAL(wp)           :: E0pcp
    REAL(wp)           :: E1pcp
    REAL(wp)           :: E0v
    REAL(wp)           :: E1v
    REAL(wp)           :: t_jmax_offset
    REAL(wp)           :: t_jmax_slope
    REAL(wp)           :: t_jmax_opt_min
    REAL(wp)           :: t_jmax_opt_max
    REAL(wp)           :: Tref_pepc
    REAL(wp)           :: Tbase_pepc
    REAL(wp)           :: alpha_i
    REAL(wp)           :: pO2
    REAL(wp)           :: CiCa_default_C3
    REAL(wp)           :: CiCa_default_C4
    REAL(wp)           :: ps_it_max
    REAL(wp)           :: ci_max
    REAL(wp)           :: Dwv2co2_air
    REAL(wp)           :: Dwv2co2_turb
    REAL(wp)           :: discr_ps_a_C13
    REAL(wp)           :: discr_ps_b_C13
    REAL(wp)           :: discr_ps_c_C13
    REAL(wp)           :: discr_ps_a_C14
    REAL(wp)           :: discr_ps_b_C14
    REAL(wp)           :: discr_ps_c_C14
    REAL(wp)           :: discr_ps_phi
    REAL(wp)           :: ka
    REAL(wp)           :: kn
    REAL(wp)           :: soa_min_temp
    REAL(wp)           :: soa_max_temp
    ! PHENOLOGY_CTL
    REAL(wp)           :: gdd_t_air_threshold
    ! RADIATION_CTL
    REAL(wp)           :: kbl0_vis
    REAL(wp)           :: kdf0_vis
    REAL(wp)           :: kbl0_nir
    REAL(wp)           :: kdf0_nir
    REAL(wp)           :: rho2sbeta
    REAL(wp)           :: k_csf
    REAL(wp)           :: min_cos_zenith_angle_rad
    REAL(wp)           :: def_alb_vis_soil
    REAL(wp)           :: def_alb_nir_soil
    REAL(wp)           :: rad2ppfd
    REAL(wp)           :: rfr_ratio_toc
    REAL(wp)           :: k_r2fr_chl
    ! SHARED_CTL
    REAL(wp)           :: mavg_period_weekly
    REAL(wp)           :: mavg_period_monthly
    REAL(wp)           :: mavg_period_yearly
    REAL(wp)           :: mavg_period_tsoa
    REAL(wp)           :: mavg_period_tuptake
    REAL(wp)           :: mavg_period_tlabile
    REAL(wp)           :: mavg_period_tphen
    REAL(wp)           :: mavg_period_tenzyme
    REAL(wp)           :: mavg_period_tfrac
    REAL(wp)           :: mavg_period_tcnl
    REAL(wp)           :: mavg_period_talloc
    REAL(wp)           :: mavg_period_talloc_dynamic
    REAL(wp)           :: mavg_period_tacclim
    REAL(wp)           :: mavg_period_tmic
    REAL(wp)           :: mavg_period_tgrowth
    REAL(wp)           :: mavg_period_tvegdyn
    ! SSE_CTL
    REAL(wp)           :: soil_heat_cap
    REAL(wp)           :: soil_therm_cond
    REAL(wp)           :: InterceptionEfficiency
    REAL(wp)           :: w_skin_max
    REAL(wp)           :: drain_frac
    REAL(wp)           :: fdrain_runoff
    REAL(wp)           :: wpot_pwp
    REAL(wp)           :: wpot_fc
    REAL(wp)           :: wpot_min
    REAL(wp)           :: k_pwp_a
    REAL(wp)           :: k_pwp_c
    REAL(wp)           :: k_pwp_s
    REAL(wp)           :: k_pwp_sc
    REAL(wp)           :: k_pwp_at
    REAL(wp)           :: k_pwp_bt
    REAL(wp)           :: k_fc_a
    REAL(wp)           :: k_fc_c
    REAL(wp)           :: k_fc_s
    REAL(wp)           :: k_fc_sc
    REAL(wp)           :: k_fc_at
    REAL(wp)           :: k_fc_bt
    REAL(wp)           :: k_fc_ct
    REAL(wp)           :: k_sat_a
    REAL(wp)           :: k_sat_c
    REAL(wp)           :: k_sat_s
    REAL(wp)           :: k_sat_sc
    REAL(wp)           :: k_sat_at
    REAL(wp)           :: k_sat_bt
    REAL(wp)           :: k_sat_ct
    REAL(wp)           :: k_sat_dt
    REAL(wp)           :: kdiff_sat_max
    ! VEGETATION_CTL
    REAL(wp)           :: fresp_growth
    REAL(wp)           :: fmaint_rate_base
    REAL(wp)           :: fmaint_rate_wood
    REAL(wp)           :: t_maint_k1
    REAL(wp)           :: t_maint_k2
    REAL(wp)           :: t_maint_k3
    REAL(wp)           :: f_resp_acclim
    REAL(wp)           :: t_acclim_zero
    REAL(wp)           :: k0_fn_pepc_C3
    REAL(wp)           :: k0_fn_pepc_C4
    REAL(wp)           :: jmax2vcmax_C3
    REAL(wp)           :: jmax2vcmax_C4
    REAL(wp)           :: tau_labile
    REAL(wp)           :: lambda_labile_temp
    REAL(wp)           :: k_labile_temp
    REAL(wp)           :: lambda_labile_theta
    REAL(wp)           :: k_labile_theta
    REAL(wp)           :: lambda_phiphen
    REAL(wp)           :: k_phiphen
    REAL(wp)           :: lambda_phimaint
    REAL(wp)           :: k_phimaint
    REAL(wp)           :: lambda_phimaint_nut
    REAL(wp)           :: k_phimaint_nut
    REAL(wp)           :: lambda_phistore
    REAL(wp)           :: k_phistore
    REAL(wp)           :: knut_labile
    REAL(wp)           :: lambda_sinklim_ps
    REAL(wp)           :: k_sinklim_ps
    REAL(wp)           :: k_cnp_sinklim_ps
    REAL(wp)           :: beta_sinklim_ps_min
    REAL(wp)           :: fstore_leaf_max
    REAL(wp)           :: fstore_fine_root_max
    REAL(wp)           :: fstore_sap_wood_max
    REAL(wp)           :: k_phi_interact
    REAL(wp)           :: resorp_fract_leaf
    REAL(wp)           :: resorp_fract_wood
    REAL(wp)           :: root2leaf_cn
    REAL(wp)           :: wood2leaf_cn
    REAL(wp)           :: root2leaf_np
    REAL(wp)           :: wood2leaf_np
    REAL(wp)           :: tau_nutrient_recycling
    REAL(wp)           :: kappa_f_demand_nc
    REAL(wp)           :: kappa_f_demand_pn
    REAL(wp)           :: k_f_demand
    REAL(wp)           :: cost_n_fixation_max
    REAL(wp)           :: cost_n_fixation_min
    REAL(wp)           :: k1_temp_bnf
    REAL(wp)           :: k2_temp_bnf
    REAL(wp)           :: k3_temp_bnf
    REAL(wp)           :: k4_temp_bnf
    REAL(wp)           :: transform_cost_nh4
    REAL(wp)           :: transform_cost_no3
    REAL(wp)           :: transform_cost_org
    REAL(wp)           :: f_nacq_nh4
    REAL(wp)           :: f_nacq_no3
    REAL(wp)           :: vmax_bnf
    REAL(wp)           :: k_bnf
    REAL(wp)           :: offset_bnf
    REAL(wp)           :: eta_nfixation
    REAL(wp)           :: k0_fn_chl_C3
    REAL(wp)           :: k0_fn_chl_C4
    REAL(wp)           :: k1_fn_chl_C3
    REAL(wp)           :: k1_fn_chl_C4
    REAL(wp)           :: kfn_chl
    REAL(wp)           :: min_lai
    REAL(wp)           :: target_lai_max
    REAL(wp)           :: max_leaf_shedding_rate
    REAL(wp)           :: k_leafon_canopy
    REAL(wp)           :: background_mort_rate_tree
    REAL(wp)           :: background_mort_rate_grass
    REAL(wp)           :: lambda_mort_light
    REAL(wp)           :: k_mort_light
    REAL(wp)           :: lambda_est_temp
    REAL(wp)           :: k_est_temp
    REAL(wp)           :: lambda_est_moist
    REAL(wp)           :: k_est_moist
    REAL(wp)           :: wood_extraction_rate
    REAL(wp)           :: k1_mort_greff
    REAL(wp)           :: k2_mort_greff
    REAL(wp)           :: k3_mort_greff
    REAL(wp)           :: min_greff
    REAL(wp)           :: max_crown_area
    REAL(wp)           :: min_diameter
    REAL(wp)           :: k_crown_area
    REAL(wp)           :: k_rp
    REAL(wp)           :: k_fpc
    REAL(wp)           :: fpc_max
    REAL(wp)           :: carbon_per_dryweight_leaf
    REAL(wp)           :: min_height
    REAL(wp)           :: sm2lm_grass
    REAL(wp)           :: lambda_fruit_alloc
    REAL(wp)           :: k1_fruit_alloc
    REAL(wp)           :: k3_fruit_alloc
    REAL(wp)           :: k4_fruit_alloc
    REAL(wp)           :: w_root_lim_max
    REAL(wp)           :: leaf2root_min_ratio
    REAL(wp)           :: leaf2root_max_ratio
    REAL(wp)           :: delta_n_fraction
    REAL(wp)           :: delta_n_leaf
    ! SB_CTL
    REAL(wp)           :: tau_soluable_litter
    REAL(wp)           :: tau_metabolic_litter
    REAL(wp)           :: tau_structural_litter
    REAL(wp)           :: tau_woody_litter
    REAL(wp)           :: tau_fast
    REAL(wp)           :: tau_slow
    REAL(wp)           :: tau_biomineralisation
    REAL(wp)           :: t_ref_tresponse
    REAL(wp)           :: ea_depolymerisation
    REAL(wp)           :: ea_mic_uptake
    REAL(wp)           :: ea_sorption
    REAL(wp)           :: ea_desorption
    REAL(wp)           :: ea_hsc
    REAL(wp)           :: ea_decomposition
    REAL(wp)           :: ed_decomposition
    REAL(wp)           :: t_ref_decomposition
    REAL(wp)           :: k1_afps
    REAL(wp)           :: kmr_depolymerisation
    REAL(wp)           :: kmr_mic_uptake
    REAL(wp)           :: kmr_hsc
    REAL(wp)           :: kmr_sorption
    REAL(wp)           :: kmr_decomposition
    REAL(wp)           :: k_diff_org
    REAL(wp)           :: rho_bulk_org
    REAL(wp)           :: carbon_per_dryweight_SOM
    REAL(wp)           :: qmax_po4_fine_particle
    REAL(wp)           :: k_adsorpt_dom
    REAL(wp)           :: k_desorpt_dom
    REAL(wp)           :: k_adsorpt_det
    REAL(wp)           :: k_desorpt_det
    REAL(wp)           :: k_adsorpt_nh4
    REAL(wp)           :: k_desorpt_nh4
    REAL(wp)           :: k_adsorpt_po4
    REAL(wp)           :: k_desorpt_po4
    REAL(wp)           :: k_adsorpt_f_po4
    REAL(wp)           :: k_desorp_f_po4
    REAL(wp)           :: k_desorp_s_po4
    REAL(wp)           :: k_adsorpt_s_po4
    REAL(wp)           :: microbial_cn
    REAL(wp)           :: microbial_np
    REAL(wp)           :: SOM_np
    REAL(wp)           :: SOM_cn
    REAL(wp)           :: microbial_nue
    REAL(wp)           :: microbial_pue
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: microbial_cue_max
    REAL(wp)           :: microbial_cue_min
    REAL(wp)           :: microbial_cn_max
    REAL(wp)           :: microbial_cn_min
    REAL(wp)           :: f_mic_cn
    REAL(wp)           :: frac_woody2dec_litter
    REAL(wp)           :: frac_litter2fast_som
    REAL(wp)           :: frac_fast2slow_som
    REAL(wp)           :: frac_slow2fast_som
    REAL(wp)           :: fc_woody2polymeric
    REAL(wp)           :: fc_soluable2dom
    REAL(wp)           :: k_mort_microbial
    REAL(wp)           :: min_mic_biomass
    REAL(wp)           :: vmax_mic_uptake
    REAL(wp)           :: km_mic_uptake
    REAL(wp)           :: km_depolymerisation
    REAL(wp)           :: vmax_depolymerisation_residue
    REAL(wp)           :: vmax_depolymerisation_litter
    REAL(wp)           :: f_soluable_residue
    REAL(wp)           :: f_soluable_leaf
    REAL(wp)           :: f_soluable_fine_root
    REAL(wp)           :: f_soluable_fruit
    REAL(wp)           :: f_soluable_seed_bed
    REAL(wp)           :: fc_soluable_max
    REAL(wp)           :: k_fc_soluable
    REAL(wp)           :: k_fn_soluable
    REAL(wp)           :: k_fp_soluable
    REAL(wp)           :: lc_leaf_max
    REAL(wp)           :: lc_leaf2sla
    REAL(wp)           :: lc_fine_root
    REAL(wp)           :: lc_coarse_root
    REAL(wp)           :: lc_woody_litter
    REAL(wp)           :: lc_fruit
    REAL(wp)           :: lc_seed_bed
    REAL(wp)           :: mycorrhiza_cue
    REAL(wp)           :: mycorrhiza_cn
    REAL(wp)           :: mycorrhiza_np
    REAL(wp)           :: f_maxresp_myc
    REAL(wp)           :: f_surface_myc2fineroot
    REAL(wp)           :: km1_up_nh4
    REAL(wp)           :: km2_up_nh4
    REAL(wp)           :: km1_up_no3
    REAL(wp)           :: km2_up_no3
    REAL(wp)           :: km1_up_po4
    REAL(wp)           :: km2_up_po4
    REAL(wp)           :: vmax_mic_uptake_nh4
    REAL(wp)           :: km_mic_uptake_nh4
    REAL(wp)           :: vmax_mic_uptake_no3
    REAL(wp)           :: km_mic_uptake_no3
    REAL(wp)           :: f_mic_carrier_n
    REAL(wp)           :: vmax_plant_uptake_nh4
    REAL(wp)           :: km_plant_uptake_nh4
    REAL(wp)           :: vmax_plant_uptake_no3
    REAL(wp)           :: km_plant_uptake_no3
    REAL(wp)           :: f_plant_carrier_n
    REAL(wp)           :: vmax_mic_uptake_po4
    REAL(wp)           :: km_mic_uptake_po4
    REAL(wp)           :: vmax_plant_uptake_po4
    REAL(wp)           :: km_plant_uptake_po4
    REAL(wp)           :: k_mic_uptake_nh4
    REAL(wp)           :: k_plant_uptake_nh4
    REAL(wp)           :: k_mic_uptake_no3
    REAL(wp)           :: k_plant_uptake_no3
    REAL(wp)           :: k_mic_uptake_po4
    REAL(wp)           :: k_plant_uptake_po4
    REAL(wp)           :: f_mic_carrier_po4
    REAL(wp)           :: f_plant_carrier_po4
    REAL(wp)           :: qmax_po4_mineral
    REAL(wp)           :: km_adsorpt_mineral_po4
    REAL(wp)           :: km_adsorpt_OM_po4
    REAL(wp)           :: ea_nitrification
    REAL(wp)           :: ed_nitrification
    REAL(wp)           :: t_opt_nitrification
    REAL(wp)           :: f_nit_noy
    REAL(wp)           :: f_nit_n2o
    REAL(wp)           :: ea_denitrification
    REAL(wp)           :: vmax_nitrification
    REAL(wp)           :: vmax_denitrification
    REAL(wp)           :: km_denitrification_c
    REAL(wp)           :: km_denitrification_no3
    REAL(wp)           :: f_denit_noy
    REAL(wp)           :: f_denit_n2o
    REAL(wp)           :: ea_gasdiffusion
    REAL(wp)           :: lambda_anvf
    REAL(wp)           :: k_anvf
    REAL(wp)           :: floss_nmin
    REAL(wp)           :: fnit_nmin
    REAL(wp)           :: fleach_nh4
    REAL(wp)           :: fleach_po4
    REAL(wp)           :: soil_n_threshold_bnf
    REAL(wp)           :: standard_bnf_max_rate
    REAL(wp)           :: eta_mic_uptake_nh4
    REAL(wp)           :: eta_mic_uptake_no3
    REAL(wp)           :: eta_plant_uptake_nh4
    REAL(wp)           :: eta_plant_uptake_no3
    REAL(wp)           :: eta_nitrification
    REAL(wp)           :: eta_nitrate_production
    REAL(wp)           :: eta_denitrification
    REAL(wp)           :: eta_ammonification
    REAL(wp)           :: eta_volatilisation_nh3
    REAL(wp)           :: k_weath_mineral
    REAL(wp)           :: eq_sorp_mineral_po4
    REAL(wp)           :: k_occlude_po4
    REAL(wp)           :: vmax_biochem_po4
    REAL(wp)           :: km_mic_biochem_po4
    REAL(wp)           :: km_dip_biochem_po4
    REAL(wp)           :: km_rootc_biochem_po4
    REAL(wp)           :: km_rootc_weath_po4
    REAL(wp)           :: biomin_thres_PC
    REAL(wp)           :: fc_p_res2dom
    REAL(wp)           :: fc_n_res2dom

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_parameter_slm_list'

    ! Namelist definition of section parameters
    ! 
    NAMELIST /assimilation_ctl/   &
        jmax2n                    , &
        vcmax2n                   , &
        pepc2n                    , &
        chl2n                     , &
        k1_fn_struc               , &
        E0kc                      , &
        E1kc                      , &
        E0ko                      , &
        E1ko                      , &
        E0pcp                     , &
        E1pcp                     , &
        E0v                       , &
        E1v                       , &
        t_jmax_offset             , &
        t_jmax_slope              , &
        t_jmax_opt_min            , &
        t_jmax_opt_max            , &
        Tref_pepc                 , &
        Tbase_pepc                , &
        alpha_i                   , &
        pO2                       , &
        CiCa_default_C3           , &
        CiCa_default_C4           , &
        ps_it_max                 , &
        ci_max                    , &
        Dwv2co2_air               , &
        Dwv2co2_turb              , &
        discr_ps_a_C13            , &
        discr_ps_b_C13            , &
        discr_ps_c_C13            , &
        discr_ps_a_C14            , &
        discr_ps_b_C14            , &
        discr_ps_c_C14            , &
        discr_ps_phi              , &
        ka                        , &
        kn                        , &
        soa_min_temp              , &
        soa_max_temp              

    NAMELIST /phenology_ctl/      &
        gdd_t_air_threshold

    NAMELIST /radiation_ctl/      &
        kbl0_vis                  , &
        kdf0_vis                  , &
        kbl0_nir                  , &
        kdf0_nir                  , &
        rho2sbeta                 , &
        k_csf                     , &
        min_cos_zenith_angle_rad  , &
        def_alb_vis_soil          , &
        def_alb_nir_soil          , &
        rad2ppfd                  , &
        rfr_ratio_toc             , &
        k_r2fr_chl              

    NAMELIST /shared_ctl/         &
        mavg_period_weekly        , &
        mavg_period_monthly       , &
        mavg_period_yearly        , &
        mavg_period_tsoa          , &
        mavg_period_tuptake       , &
        mavg_period_tlabile       , &
        mavg_period_tphen         , &
        mavg_period_tenzyme       , &
        mavg_period_tfrac         , &
        mavg_period_tcnl          , &
        mavg_period_talloc        , &
        mavg_period_talloc_dynamic, &
        mavg_period_tacclim       , &
        mavg_period_tmic          , &
        mavg_period_tgrowth       , &
        mavg_period_tvegdyn       

    NAMELIST /sse_ctl/            &  
        soil_heat_cap             , &
        soil_therm_cond           , &
        InterceptionEfficiency    , &
        w_skin_max                , &
        drain_frac                , &
        fdrain_runoff             , &
        wpot_pwp                  , &
        wpot_fc                   , &
        wpot_min                  , &
        k_pwp_a                   , &
        k_pwp_c                   , &
        k_pwp_s                   , &
        k_pwp_sc                  , &
        k_pwp_at                  , &
        k_pwp_bt                  , &
        k_fc_a                    , &
        k_fc_c                    , &
        k_fc_s                    , &
        k_fc_sc                   , &
        k_fc_at                   , &
        k_fc_bt                   , &
        k_fc_ct                   , &
        k_sat_a                   , &
        k_sat_c                   , &
        k_sat_s                   , &
        k_sat_sc                  , &
        k_sat_at                  , &
        k_sat_bt                  , &
        k_sat_ct                  , &
        k_sat_dt                  , &
        kdiff_sat_max             

    NAMELIST /vegetation_ctl/     &  
        fresp_growth              , &
        fmaint_rate_base          , &
        fmaint_rate_wood          , &
        t_maint_k1                , &
        t_maint_k2                , &
        t_maint_k3                , &
        f_resp_acclim             , &
        t_acclim_zero             , &
        k0_fn_pepc_C3             , &
        k0_fn_pepc_C4             , &
        jmax2vcmax_C3             , &
        jmax2vcmax_C4             , &
        tau_labile                , &
        lambda_labile_temp        , &
        k_labile_temp             , &
        lambda_labile_theta       , &
        k_labile_theta            , &
        lambda_phiphen            , &
        k_phiphen                 , &
        lambda_phimaint           , &
        k_phimaint                , &
        lambda_phimaint_nut       , &
        k_phimaint_nut            , &
        lambda_phistore           , &
        k_phistore                , &
        knut_labile               , &
        lambda_sinklim_ps         , &
        k_sinklim_ps              , &
        k_cnp_sinklim_ps          , &
        beta_sinklim_ps_min       , &
        fstore_leaf_max           , &
        fstore_fine_root_max      , &
        fstore_sap_wood_max       , &
        k_phi_interact            , &
        resorp_fract_leaf         , &
        resorp_fract_wood         , &
        root2leaf_cn              , &
        wood2leaf_cn              , &
        root2leaf_np              , &
        wood2leaf_np              , &
        tau_nutrient_recycling    , &
        kappa_f_demand_nc         , &
        kappa_f_demand_pn         , &
        k_f_demand                , &
        cost_n_fixation_max       , &
        cost_n_fixation_min       , &
        k1_temp_bnf               , &
        k2_temp_bnf               , &
        k3_temp_bnf               , &
        k4_temp_bnf               , &
        transform_cost_nh4        , &
        transform_cost_no3        , &
        transform_cost_org        , &
        f_nacq_nh4                , &
        f_nacq_no3                , &
        vmax_bnf                  , &
        k_bnf                     , &
        offset_bnf                , &
        eta_nfixation             , &
        k0_fn_chl_C3              , &
        k0_fn_chl_C4              , &
        k1_fn_chl_C3              , &
        k1_fn_chl_C4              , &
        kfn_chl                   , &
        min_lai                   , &
        target_lai_max            , &
        max_leaf_shedding_rate    , &
        k_leafon_canopy           , &
        background_mort_rate_tree , &
        background_mort_rate_grass, &
        lambda_mort_light         , &
        k_mort_light              , &
        lambda_est_temp           , &
        k_est_temp                , &
        lambda_est_moist          , &
        k_est_moist               , &
        wood_extraction_rate      , &
        k1_mort_greff             , &
        k2_mort_greff             , &
        k3_mort_greff             , &
        min_greff                 , &
        max_crown_area            , &
        min_diameter              , &
        k_crown_area              , &
        k_rp                      , &
        k_fpc                     , &
        fpc_max                   , &
        carbon_per_dryweight_leaf , &
        min_height                , &
        sm2lm_grass               , &
        lambda_fruit_alloc        , &
        k1_fruit_alloc            , &
        k3_fruit_alloc            , &
        k4_fruit_alloc            , &
        w_root_lim_max            , &
        leaf2root_min_ratio       , &
        leaf2root_max_ratio       , &
        delta_n_fraction          , &
        delta_n_leaf              

    NAMELIST /sb_ctl/             &
        tau_soluable_litter           , &
        tau_metabolic_litter          , &
        tau_structural_litter         , &
        tau_woody_litter              , &
        tau_fast                      , &
        tau_slow                      , &
        tau_biomineralisation         , &
        t_ref_tresponse               , &
        ea_depolymerisation           , &
        ea_mic_uptake                 , &
        ea_sorption                   , &
        ea_desorption                 , &
        ea_hsc                        , &
        ea_decomposition              , &
        ed_decomposition              , &
        t_ref_decomposition           , &
        k1_afps                       , &
        kmr_depolymerisation          , &
        kmr_mic_uptake                , &
        kmr_hsc                       , &
        kmr_sorption                  , &
        kmr_decomposition             , &
        k_diff_org                    , &
        rho_bulk_org                  , &
        carbon_per_dryweight_SOM      , &
        qmax_po4_fine_particle        , &
        k_adsorpt_dom                 , &
        k_desorpt_dom                 , &
        k_adsorpt_det                 , &
        k_desorpt_det                 , &
        k_adsorpt_nh4                 , &
        k_desorpt_nh4                 , &
        k_adsorpt_po4                 , &
        k_desorpt_po4                 , &
        k_adsorpt_f_po4               , &
        k_desorp_f_po4                , &
        k_desorp_s_po4                , &
        k_adsorpt_s_po4               , &
        microbial_cn                  , &
        microbial_np                  , &
        SOM_np                        , &
        SOM_cn                        , &
        microbial_nue                 , &
        microbial_pue                 , &
        microbial_cue                 , &
        microbial_cue_max             , &
        microbial_cue_min             , &
        microbial_cn_max              , &
        microbial_cn_min              , &
        f_mic_cn                      , &
        frac_woody2dec_litter         , &
        frac_litter2fast_som          , &
        frac_fast2slow_som            , &
        frac_slow2fast_som            , &
        fc_woody2polymeric            , &
        fc_soluable2dom               , &
        k_mort_microbial              , &
        min_mic_biomass               , &
        vmax_mic_uptake               , &
        km_mic_uptake                 , &
        km_depolymerisation           , &
        vmax_depolymerisation_residue , &
        vmax_depolymerisation_litter  , &
        f_soluable_residue            , &
        f_soluable_leaf               , &
        f_soluable_fine_root          , &
        f_soluable_fruit              , &
        f_soluable_seed_bed           , &
        fc_soluable_max               , &
        k_fc_soluable                 , &
        k_fn_soluable                 , &
        k_fp_soluable                 , &
        lc_leaf_max                   , &
        lc_leaf2sla                   , &
        lc_fine_root                  , &
        lc_coarse_root                , &
        lc_woody_litter               , &
        lc_fruit                      , &
        lc_seed_bed                   , &
        mycorrhiza_cue                , &
        mycorrhiza_cn                 , &
        mycorrhiza_np                 , &
        f_maxresp_myc                 , &
        f_surface_myc2fineroot        , &
        km1_up_nh4                    , &
        km2_up_nh4                    , &
        km1_up_no3                    , &
        km2_up_no3                    , &
        km1_up_po4                    , &
        km2_up_po4                    , &
        vmax_mic_uptake_nh4           , &
        km_mic_uptake_nh4             , &
        vmax_mic_uptake_no3           , &
        km_mic_uptake_no3             , &
        f_mic_carrier_n               , &
        vmax_plant_uptake_nh4         , &
        km_plant_uptake_nh4           , &
        vmax_plant_uptake_no3         , &
        km_plant_uptake_no3           , &
        f_plant_carrier_n             , &
        vmax_mic_uptake_po4           , &
        km_mic_uptake_po4             , &
        vmax_plant_uptake_po4         , &
        km_plant_uptake_po4           , &
        k_mic_uptake_nh4              , &
        k_plant_uptake_nh4            , &
        k_mic_uptake_no3              , &
        k_plant_uptake_no3            , &
        k_mic_uptake_po4              , &
        k_plant_uptake_po4            , &
        f_mic_carrier_po4             , &
        f_plant_carrier_po4           , &
        qmax_po4_mineral              , &
        km_adsorpt_mineral_po4        , &
        km_adsorpt_OM_po4             , &
        ea_nitrification              , &
        ed_nitrification              , &
        t_opt_nitrification           , &
        f_nit_noy                     , &
        f_nit_n2o                     , &
        ea_denitrification            , &
        vmax_nitrification            , &
        vmax_denitrification          , &
        km_denitrification_c          , &
        km_denitrification_no3        , &
        f_denit_noy                   , &
        f_denit_n2o                   , &
        ea_gasdiffusion               , &
        lambda_anvf                   , &
        k_anvf                        , &
        floss_nmin                    , &
        fnit_nmin                     , &
        fleach_nh4                    , &
        fleach_po4                    , &
        soil_n_threshold_bnf          , &
        standard_bnf_max_rate         , &
        eta_mic_uptake_nh4            , &
        eta_mic_uptake_no3            , &
        eta_plant_uptake_nh4          , &
        eta_plant_uptake_no3          , &
        eta_nitrification             , &
        eta_nitrate_production        , &
        eta_denitrification           , &
        eta_ammonification            , &
        eta_volatilisation_nh3        , &
        k_weath_mineral               , &
        eq_sorp_mineral_po4           , &
        k_occlude_po4                 , &
        vmax_biochem_po4              , &
        km_mic_biochem_po4            , &
        km_dip_biochem_po4            , &
        km_rootc_biochem_po4          , &
        km_rootc_weath_po4            , &
        biomin_thres_PC               , &
        fc_p_res2dom                  , &
        fc_n_res2dom                  


    ! variables for reading model-options from namelist
    ! 
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting to read parameter list')

    ! Set defaults
    ! 
    ! ASSIMILATION_CTL
    jmax2n                      = missval
    vcmax2n                     = missval
    pepc2n                      = missval
    chl2n                       = missval
    k1_fn_struc                 = missval
    E0kc                        = missval
    E1kc                        = missval
    E0ko                        = missval
    E1ko                        = missval
    E0pcp                       = missval
    E1pcp                       = missval
    E0v                         = missval
    E1v                         = missval
    t_jmax_offset               = missval
    t_jmax_slope                = missval
    t_jmax_opt_min              = missval
    t_jmax_opt_max              = missval
    Tref_pepc                   = missval
    Tbase_pepc                  = missval
    alpha_i                     = missval
    pO2                         = missval
    CiCa_default_C3             = missval
    CiCa_default_C4             = missval
    ps_it_max                   = missval
    ci_max                      = missval
    Dwv2co2_air                 = missval
    Dwv2co2_turb                = missval
    discr_ps_a_C13              = missval
    discr_ps_b_C13              = missval
    discr_ps_c_C13              = missval
    discr_ps_a_C14              = missval
    discr_ps_b_C14              = missval
    discr_ps_c_C14              = missval
    discr_ps_phi                = missval
    ka                          = missval
    kn                          = missval
    soa_min_temp                = missval
    soa_max_temp                = missval
    ! PHENOLOGY_CTL
    gdd_t_air_threshold         = missval
    ! RADIATION_CTL
    kbl0_vis                    = missval
    kdf0_vis                    = missval
    kbl0_nir                    = missval
    kdf0_nir                    = missval
    rho2sbeta                   = missval
    k_csf                       = missval
    min_cos_zenith_angle_rad    = missval
    def_alb_vis_soil            = missval
    def_alb_nir_soil            = missval
    rad2ppfd                    = missval
    rfr_ratio_toc               = missval
    k_r2fr_chl                  = missval
    ! SHARED_CTL
    mavg_period_weekly          = missval
    mavg_period_monthly         = missval
    mavg_period_yearly          = missval
    mavg_period_tsoa            = missval
    mavg_period_tuptake         = missval
    mavg_period_tlabile         = missval
    mavg_period_tphen           = missval
    mavg_period_tenzyme         = missval
    mavg_period_tfrac           = missval
    mavg_period_tcnl            = missval
    mavg_period_talloc          = missval
    mavg_period_talloc_dynamic  = missval
    mavg_period_tacclim         = missval
    mavg_period_tmic            = missval
    mavg_period_tgrowth         = missval
    mavg_period_tvegdyn         = missval
    ! SSE_CTL
    soil_heat_cap               = missval
    soil_therm_cond             = missval
    InterceptionEfficiency      = missval
    w_skin_max                  = missval
    drain_frac                  = missval
    fdrain_runoff               = missval
    wpot_pwp                    = missval
    wpot_fc                     = missval
    wpot_min                    = missval
    k_pwp_a                     = missval
    k_pwp_c                     = missval
    k_pwp_s                     = missval
    k_pwp_sc                    = missval
    k_pwp_at                    = missval
    k_pwp_bt                    = missval
    k_fc_a                      = missval
    k_fc_c                      = missval
    k_fc_s                      = missval
    k_fc_sc                     = missval
    k_fc_at                     = missval
    k_fc_bt                     = missval
    k_fc_ct                     = missval
    k_sat_a                     = missval
    k_sat_c                     = missval
    k_sat_s                     = missval
    k_sat_sc                    = missval
    k_sat_at                    = missval
    k_sat_bt                    = missval
    k_sat_ct                    = missval
    k_sat_dt                    = missval
    kdiff_sat_max               = missval
    ! VEGETATION_CTL
    fresp_growth                = missval
    fmaint_rate_base            = missval
    fmaint_rate_wood            = missval
    t_maint_k1                  = missval
    t_maint_k2                  = missval
    t_maint_k3                  = missval
    f_resp_acclim               = missval
    t_acclim_zero               = missval
    k0_fn_pepc_C3               = missval
    k0_fn_pepc_C4               = missval
    jmax2vcmax_C3               = missval
    jmax2vcmax_C4               = missval
    tau_labile                  = missval
    lambda_labile_temp          = missval
    k_labile_temp               = missval
    lambda_labile_theta         = missval
    k_labile_theta              = missval
    lambda_phiphen              = missval
    k_phiphen                   = missval
    lambda_phimaint             = missval
    k_phimaint                  = missval
    lambda_phimaint_nut         = missval
    k_phimaint_nut              = missval
    lambda_phistore             = missval
    k_phistore                  = missval
    knut_labile                 = missval
    lambda_sinklim_ps           = missval
    k_sinklim_ps                = missval
    k_cnp_sinklim_ps            = missval
    beta_sinklim_ps_min         = missval
    fstore_leaf_max             = missval
    fstore_fine_root_max        = missval
    fstore_sap_wood_max         = missval
    k_phi_interact              = missval
    resorp_fract_leaf           = missval
    resorp_fract_wood           = missval
    root2leaf_cn                = missval
    wood2leaf_cn                = missval
    root2leaf_np                = missval
    wood2leaf_np                = missval
    tau_nutrient_recycling      = missval
    kappa_f_demand_nc           = missval
    kappa_f_demand_pn           = missval
    k_f_demand                  = missval
    cost_n_fixation_max         = missval
    cost_n_fixation_min         = missval
    k1_temp_bnf                 = missval
    k2_temp_bnf                 = missval
    k3_temp_bnf                 = missval
    k4_temp_bnf                 = missval
    transform_cost_nh4          = missval
    transform_cost_no3          = missval
    transform_cost_org          = missval
    f_nacq_nh4                  = missval
    f_nacq_no3                  = missval
    vmax_bnf                    = missval
    k_bnf                       = missval
    offset_bnf                  = missval
    eta_nfixation               = missval
    k0_fn_chl_C3                = missval
    k0_fn_chl_C4                = missval
    k1_fn_chl_C3                = missval
    k1_fn_chl_C4                = missval
    kfn_chl                     = missval
    min_lai                     = missval
    target_lai_max              = missval
    max_leaf_shedding_rate      = missval
    k_leafon_canopy             = missval
    background_mort_rate_tree   = missval
    background_mort_rate_grass  = missval
    lambda_mort_light           = missval
    k_mort_light                = missval
    lambda_est_temp             = missval
    k_est_temp                  = missval
    lambda_est_moist            = missval
    k_est_moist                 = missval
    wood_extraction_rate        = missval
    k1_mort_greff               = missval
    k2_mort_greff               = missval
    k3_mort_greff               = missval
    min_greff                   = missval
    max_crown_area              = missval
    min_diameter                = missval
    k_crown_area                = missval
    k_rp                        = missval
    k_fpc                       = missval
    fpc_max                     = missval
    carbon_per_dryweight_leaf   = missval
    min_height                  = missval
    sm2lm_grass                 = missval
    lambda_fruit_alloc          = missval
    k1_fruit_alloc              = missval
    k3_fruit_alloc              = missval
    k4_fruit_alloc              = missval
    w_root_lim_max              = missval
    leaf2root_min_ratio         = missval
    leaf2root_max_ratio         = missval
    delta_n_fraction            = missval
    delta_n_leaf                = missval
    ! SB_CTL
    tau_soluable_litter           = missval
    tau_metabolic_litter          = missval
    tau_structural_litter         = missval
    tau_woody_litter              = missval
    tau_fast                      = missval
    tau_slow                      = missval
    tau_biomineralisation         = missval
    t_ref_tresponse               = missval
    ea_depolymerisation           = missval
    ea_mic_uptake                 = missval
    ea_sorption                   = missval
    ea_desorption                 = missval
    ea_hsc                        = missval
    ea_decomposition              = missval
    ed_decomposition              = missval
    t_ref_decomposition           = missval
    k1_afps                       = missval
    kmr_depolymerisation          = missval
    kmr_mic_uptake                = missval
    kmr_hsc                       = missval
    kmr_sorption                  = missval
    kmr_decomposition             = missval
    k_diff_org                    = missval
    rho_bulk_org                  = missval
    carbon_per_dryweight_SOM      = missval
    qmax_po4_fine_particle        = missval
    k_adsorpt_dom                 = missval
    k_desorpt_dom                 = missval
    k_adsorpt_det                 = missval
    k_desorpt_det                 = missval
    k_adsorpt_nh4                 = missval
    k_desorpt_nh4                 = missval
    k_adsorpt_po4                 = missval
    k_desorpt_po4                 = missval
    k_adsorpt_f_po4               = missval
    k_desorp_f_po4                = missval
    k_desorp_s_po4                = missval
    k_adsorpt_s_po4               = missval
    microbial_cn                  = missval
    microbial_np                  = missval
    SOM_np                        = missval
    SOM_cn                        = missval
    microbial_nue                 = missval
    microbial_pue                 = missval
    microbial_cue                 = missval
    microbial_cue_max             = missval
    microbial_cue_min             = missval
    microbial_cn_max              = missval
    microbial_cn_min              = missval
    f_mic_cn                      = missval
    frac_woody2dec_litter         = missval
    frac_litter2fast_som          = missval
    frac_fast2slow_som            = missval
    frac_slow2fast_som            = missval
    fc_woody2polymeric            = missval
    fc_soluable2dom               = missval
    k_mort_microbial              = missval
    min_mic_biomass               = missval
    vmax_mic_uptake               = missval
    km_mic_uptake                 = missval
    km_depolymerisation           = missval
    vmax_depolymerisation_residue = missval
    vmax_depolymerisation_litter  = missval
    f_soluable_residue            = missval
    f_soluable_leaf               = missval
    f_soluable_fine_root          = missval
    f_soluable_fruit              = missval
    f_soluable_seed_bed           = missval
    fc_soluable_max               = missval
    k_fc_soluable                 = missval
    k_fn_soluable                 = missval
    k_fp_soluable                 = missval
    lc_leaf_max                   = missval
    lc_leaf2sla                   = missval
    lc_fine_root                  = missval
    lc_coarse_root                = missval
    lc_woody_litter               = missval
    lc_fruit                      = missval
    lc_seed_bed                   = missval
    mycorrhiza_cue                = missval
    mycorrhiza_cn                 = missval
    mycorrhiza_np                 = missval
    f_maxresp_myc                 = missval
    f_surface_myc2fineroot        = missval
    km1_up_nh4                    = missval
    km2_up_nh4                    = missval
    km1_up_no3                    = missval
    km2_up_no3                    = missval
    km1_up_po4                    = missval
    km2_up_po4                    = missval
    vmax_mic_uptake_nh4           = missval
    km_mic_uptake_nh4             = missval
    vmax_mic_uptake_no3           = missval
    km_mic_uptake_no3             = missval
    f_mic_carrier_n               = missval
    vmax_plant_uptake_nh4         = missval
    km_plant_uptake_nh4           = missval
    vmax_plant_uptake_no3         = missval
    km_plant_uptake_no3           = missval
    f_plant_carrier_n             = missval
    vmax_mic_uptake_po4           = missval
    km_mic_uptake_po4             = missval
    vmax_plant_uptake_po4         = missval
    km_plant_uptake_po4           = missval
    k_mic_uptake_nh4              = missval
    k_plant_uptake_nh4            = missval
    k_mic_uptake_no3              = missval
    k_plant_uptake_no3            = missval
    k_mic_uptake_po4              = missval
    k_plant_uptake_po4            = missval
    f_mic_carrier_po4             = missval
    f_plant_carrier_po4           = missval
    qmax_po4_mineral              = missval
    km_adsorpt_mineral_po4        = missval
    km_adsorpt_OM_po4             = missval
    ea_nitrification              = missval
    ed_nitrification              = missval
    t_opt_nitrification           = missval
    f_nit_noy                     = missval
    f_nit_n2o                     = missval
    ea_denitrification            = missval
    vmax_nitrification            = missval
    vmax_denitrification          = missval
    km_denitrification_c          = missval
    km_denitrification_no3        = missval
    f_denit_noy                   = missval
    f_denit_n2o                   = missval
    ea_gasdiffusion               = missval
    lambda_anvf                   = missval
    k_anvf                        = missval
    floss_nmin                    = missval
    fnit_nmin                     = missval
    fleach_nh4                    = missval
    fleach_po4                    = missval
    soil_n_threshold_bnf          = missval
    standard_bnf_max_rate         = missval
    eta_mic_uptake_nh4            = missval
    eta_mic_uptake_no3            = missval
    eta_plant_uptake_nh4          = missval
    eta_plant_uptake_no3          = missval
    eta_nitrification             = missval
    eta_nitrate_production        = missval
    eta_denitrification           = missval
    eta_ammonification            = missval
    eta_volatilisation_nh3        = missval
    k_weath_mineral               = missval
    eq_sorp_mineral_po4           = missval
    k_occlude_po4                 = missval
    vmax_biochem_po4              = missval
    km_mic_biochem_po4            = missval
    km_dip_biochem_po4            = missval
    km_rootc_biochem_po4          = missval
    km_rootc_weath_po4            = missval
    biomin_thres_PC               = missval
    fc_p_res2dom                  = missval
    fc_n_res2dom                  = missval

    ! read the namelist
    ! 
    ! &ASSIMILATION_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('ASSIMILATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, assimilation_ctl)
    CALL close_nml(nml_handler)
    ! &PHENOLOGY_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('PHENOLOGY_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, phenology_ctl)
    CALL close_nml(nml_handler)
    ! &RADIATION_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('RADIATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, radiation_ctl)
    CALL close_nml(nml_handler)
    ! &SHARED_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('SHARED_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, shared_ctl)
    CALL close_nml(nml_handler)
    ! &SSE_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('SSE_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, sse_ctl)
    CALL close_nml(nml_handler)
    ! &VEGETATION_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('VEGETATION_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, vegetation_ctl)
    CALL close_nml(nml_handler)
    ! &SB_CTL
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('SB_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, sb_ctl)
    CALL close_nml(nml_handler)

    ! pass values as read from file
    ! 
    ! ASSIMILATION_CTL
    this%jmax2n                      = jmax2n          
    this%vcmax2n                     = vcmax2n         
    this%pepc2n                      = pepc2n          
    this%chl2n                       = chl2n           
    this%k1_fn_struc                 = k1_fn_struc     
    this%E0kc                        = E0kc            
    this%E1kc                        = E1kc            
    this%E0ko                        = E0ko            
    this%E1ko                        = E1ko            
    this%E0pcp                       = E0pcp           
    this%E1pcp                       = E1pcp           
    this%E0v                         = E0v             
    this%E1v                         = E1v             
    this%t_jmax_offset               = t_jmax_offset   
    this%t_jmax_slope                = t_jmax_slope    
    this%t_jmax_opt_min              = t_jmax_opt_min  
    this%t_jmax_opt_max              = t_jmax_opt_max  
    this%Tref_pepc                   = Tref_pepc       
    this%Tbase_pepc                  = Tbase_pepc      
    this%alpha_i                     = alpha_i         
    this%pO2                         = pO2             
    this%CiCa_default_C3             = CiCa_default_C3 
    this%CiCa_default_C4             = CiCa_default_C4 
    this%ps_it_max                   = ps_it_max       
    this%ci_max                      = ci_max          
    this%Dwv2co2_air                 = Dwv2co2_air     
    this%Dwv2co2_turb                = Dwv2co2_turb    
    this%discr_ps_a_C13              = discr_ps_a_C13  
    this%discr_ps_b_C13              = discr_ps_b_C13  
    this%discr_ps_c_C13              = discr_ps_c_C13  
    this%discr_ps_a_C14              = discr_ps_a_C14  
    this%discr_ps_b_C14              = discr_ps_b_C14  
    this%discr_ps_c_C14              = discr_ps_c_C14  
    this%discr_ps_phi                = discr_ps_phi    
    this%ka                          = ka              
    this%kn                          = kn              
    this%soa_min_temp                = soa_min_temp    
    this%soa_max_temp                = soa_max_temp    
    ! PHENOLOGY_CTL
    this%gdd_t_air_threshold         = gdd_t_air_threshold
    ! RADIATION_CTL
    this%kbl0_vis                    = kbl0_vis
    this%kdf0_vis                    = kdf0_vis
    this%kbl0_nir                    = kbl0_nir
    this%kdf0_nir                    = kdf0_nir
    this%rho2sbeta                   = rho2sbeta
    this%k_csf                       = k_csf
    this%min_cos_zenith_angle_rad    = min_cos_zenith_angle_rad
    this%def_alb_vis_soil            = def_alb_vis_soil
    this%def_alb_nir_soil            = def_alb_nir_soil
    this%rad2ppfd                    = rad2ppfd
    this%rfr_ratio_toc               = rfr_ratio_toc
    this%k_r2fr_chl                  = k_r2fr_chl
    ! SHARED_CTL
    this%mavg_period_weekly          = mavg_period_weekly
    this%mavg_period_monthly         = mavg_period_monthly
    this%mavg_period_yearly          = mavg_period_yearly
    this%mavg_period_tsoa            = mavg_period_tsoa
    this%mavg_period_tuptake         = mavg_period_tuptake
    this%mavg_period_tlabile         = mavg_period_tlabile
    this%mavg_period_tphen           = mavg_period_tphen
    this%mavg_period_tenzyme         = mavg_period_tenzyme
    this%mavg_period_tfrac           = mavg_period_tfrac
    this%mavg_period_tcnl            = mavg_period_tcnl
    this%mavg_period_talloc          = mavg_period_talloc
    this%mavg_period_talloc_dynamic  = mavg_period_talloc_dynamic
    this%mavg_period_tacclim         = mavg_period_tacclim
    this%mavg_period_tmic            = mavg_period_tmic
    this%mavg_period_tgrowth         = mavg_period_tgrowth
    this%mavg_period_tvegdyn         = mavg_period_tvegdyn
    ! SSE_CTL
    this%soil_heat_cap               = soil_heat_cap
    this%soil_therm_cond             = soil_therm_cond
    this%InterceptionEfficiency      = InterceptionEfficiency
    this%w_skin_max                  = w_skin_max
    this%drain_frac                  = drain_frac
    this%fdrain_runoff               = fdrain_runoff
    this%wpot_pwp                    = wpot_pwp
    this%wpot_fc                     = wpot_fc
    this%wpot_min                    = wpot_min
    this%k_pwp_a                     = k_pwp_a
    this%k_pwp_c                     = k_pwp_c
    this%k_pwp_s                     = k_pwp_s
    this%k_pwp_sc                    = k_pwp_sc
    this%k_pwp_at                    = k_pwp_at
    this%k_pwp_bt                    = k_pwp_bt
    this%k_fc_a                      = k_fc_a
    this%k_fc_c                      = k_fc_c
    this%k_fc_s                      = k_fc_s
    this%k_fc_sc                     = k_fc_sc
    this%k_fc_at                     = k_fc_at
    this%k_fc_bt                     = k_fc_bt
    this%k_fc_ct                     = k_fc_ct
    this%k_sat_a                     = k_sat_a
    this%k_sat_c                     = k_sat_c
    this%k_sat_s                     = k_sat_s
    this%k_sat_sc                    = k_sat_sc
    this%k_sat_at                    = k_sat_at
    this%k_sat_bt                    = k_sat_bt
    this%k_sat_ct                    = k_sat_ct
    this%k_sat_dt                    = k_sat_dt
    this%kdiff_sat_max               = kdiff_sat_max
    ! VEGETATION_CTL
    this%fresp_growth                = fresp_growth
    this%fmaint_rate_base            = fmaint_rate_base
    this%fmaint_rate_wood            = fmaint_rate_wood
    this%t_maint_k1                  = t_maint_k1
    this%t_maint_k2                  = t_maint_k2
    this%t_maint_k3                  = t_maint_k3
    this%f_resp_acclim               = f_resp_acclim
    this%t_acclim_zero               = t_acclim_zero
    this%k0_fn_pepc_C3               = k0_fn_pepc_C3
    this%k0_fn_pepc_C4               = k0_fn_pepc_C4
    this%jmax2vcmax_C3               = jmax2vcmax_C3
    this%jmax2vcmax_C4               = jmax2vcmax_C4
    this%tau_labile                  = tau_labile
    this%lambda_labile_temp          = lambda_labile_temp
    this%k_labile_temp               = k_labile_temp
    this%lambda_labile_theta         = lambda_labile_theta
    this%k_labile_theta              = k_labile_theta
    this%lambda_phiphen              = lambda_phiphen
    this%k_phiphen                   = k_phiphen
    this%lambda_phimaint             = lambda_phimaint
    this%k_phimaint                  = k_phimaint
    this%lambda_phimaint_nut         = lambda_phimaint_nut
    this%k_phimaint_nut              = k_phimaint_nut
    this%lambda_phistore             = lambda_phistore
    this%k_phistore                  = k_phistore
    this%knut_labile                 = knut_labile
    this%lambda_sinklim_ps           = lambda_sinklim_ps
    this%k_sinklim_ps                = k_sinklim_ps
    this%k_cnp_sinklim_ps            = k_cnp_sinklim_ps
    this%beta_sinklim_ps_min         = beta_sinklim_ps_min
    this%fstore_leaf_max             = fstore_leaf_max
    this%fstore_fine_root_max        = fstore_fine_root_max
    this%fstore_sap_wood_max         = fstore_sap_wood_max
    this%k_phi_interact              = k_phi_interact
    this%resorp_fract_leaf           = resorp_fract_leaf
    this%resorp_fract_wood           = resorp_fract_wood
    this%root2leaf_cn                = root2leaf_cn
    this%wood2leaf_cn                = wood2leaf_cn
    this%root2leaf_np                = root2leaf_np
    this%wood2leaf_np                = wood2leaf_np
    this%tau_nutrient_recycling      = tau_nutrient_recycling
    this%kappa_f_demand_nc           = kappa_f_demand_nc
    this%kappa_f_demand_pn           = kappa_f_demand_pn
    this%k_f_demand                  = k_f_demand
    this%cost_n_fixation_max         = cost_n_fixation_max
    this%cost_n_fixation_min         = cost_n_fixation_min
    this%k1_temp_bnf                 = k1_temp_bnf
    this%k2_temp_bnf                 = k2_temp_bnf
    this%k3_temp_bnf                 = k3_temp_bnf
    this%k4_temp_bnf                 = k4_temp_bnf
    this%transform_cost_nh4          = transform_cost_nh4
    this%transform_cost_no3          = transform_cost_no3
    this%transform_cost_org          = transform_cost_org
    this%f_nacq_nh4                  = f_nacq_nh4
    this%f_nacq_no3                  = f_nacq_no3
    this%vmax_bnf                    = vmax_bnf
    this%k_bnf                       = k_bnf
    this%offset_bnf                  = offset_bnf
    this%eta_nfixation               = eta_nfixation
    this%k0_fn_chl_C3                = k0_fn_chl_C3
    this%k0_fn_chl_C4                = k0_fn_chl_C4
    this%k1_fn_chl_C3                = k1_fn_chl_C3
    this%k1_fn_chl_C4                = k1_fn_chl_C4
    this%kfn_chl                     = kfn_chl
    this%min_lai                     = min_lai
    this%target_lai_max              = target_lai_max
    this%max_leaf_shedding_rate      = max_leaf_shedding_rate
    this%k_leafon_canopy             = k_leafon_canopy
    this%background_mort_rate_tree   = background_mort_rate_tree
    this%background_mort_rate_grass  = background_mort_rate_grass
    this%lambda_mort_light           = lambda_mort_light
    this%k_mort_light                = k_mort_light
    this%lambda_est_temp             = lambda_est_temp
    this%k_est_temp                  = k_est_temp
    this%lambda_est_moist            = lambda_est_moist
    this%k_est_moist                 = k_est_moist
    this%wood_extraction_rate        = wood_extraction_rate
    this%k1_mort_greff               = k1_mort_greff
    this%k2_mort_greff               = k2_mort_greff
    this%k3_mort_greff               = k3_mort_greff
    this%min_greff                   = min_greff
    this%max_crown_area              = max_crown_area
    this%min_diameter                = min_diameter
    this%k_crown_area                = k_crown_area
    this%k_rp                        = k_rp
    this%k_fpc                       = k_fpc
    this%fpc_max                     = fpc_max
    this%carbon_per_dryweight_leaf   = carbon_per_dryweight_leaf
    this%min_height                  = min_height
    this%sm2lm_grass                 = sm2lm_grass
    this%lambda_fruit_alloc          = lambda_fruit_alloc
    this%k1_fruit_alloc              = k1_fruit_alloc
    this%k3_fruit_alloc              = k3_fruit_alloc
    this%k4_fruit_alloc              = k4_fruit_alloc
    this%w_root_lim_max              = w_root_lim_max
    this%leaf2root_min_ratio         = leaf2root_min_ratio
    this%leaf2root_max_ratio         = leaf2root_max_ratio
    this%delta_n_fraction            = delta_n_fraction
    this%delta_n_leaf                = delta_n_leaf
    ! SB_CTL
    this%tau_soluable_litter           = tau_soluable_litter
    this%tau_metabolic_litter          = tau_metabolic_litter
    this%tau_structural_litter         = tau_structural_litter
    this%tau_woody_litter              = tau_woody_litter
    this%tau_fast                      = tau_fast
    this%tau_slow                      = tau_slow
    this%tau_biomineralisation         = tau_biomineralisation
    this%t_ref_tresponse               = t_ref_tresponse
    this%ea_depolymerisation           = ea_depolymerisation
    this%ea_mic_uptake                 = ea_mic_uptake
    this%ea_sorption                   = ea_sorption
    this%ea_desorption                 = ea_desorption
    this%ea_hsc                        = ea_hsc
    this%ea_decomposition              = ea_decomposition
    this%ed_decomposition              = ed_decomposition
    this%t_ref_decomposition           = t_ref_decomposition
    this%k1_afps                       = k1_afps
    this%kmr_depolymerisation          = kmr_depolymerisation
    this%kmr_mic_uptake                = kmr_mic_uptake
    this%kmr_hsc                       = kmr_hsc
    this%kmr_sorption                  = kmr_sorption
    this%kmr_decomposition             = kmr_decomposition
    this%k_diff_org                    = k_diff_org
    this%rho_bulk_org                  = rho_bulk_org
    this%carbon_per_dryweight_SOM      = carbon_per_dryweight_SOM
    this%qmax_po4_fine_particle        = qmax_po4_fine_particle
    this%k_adsorpt_dom                 = k_adsorpt_dom
    this%k_desorpt_dom                 = k_desorpt_dom
    this%k_adsorpt_det                 = k_adsorpt_det
    this%k_desorpt_det                 = k_desorpt_det
    this%k_adsorpt_nh4                 = k_adsorpt_nh4
    this%k_desorpt_nh4                 = k_desorpt_nh4
    this%k_adsorpt_po4                 = k_adsorpt_po4
    this%k_desorpt_po4                 = k_desorpt_po4
    this%k_adsorpt_f_po4               = k_adsorpt_f_po4
    this%k_desorp_f_po4                = k_desorp_f_po4
    this%k_desorp_s_po4                = k_desorp_s_po4
    this%k_adsorpt_s_po4               = k_adsorpt_s_po4
    this%microbial_cn                  = microbial_cn
    this%microbial_np                  = microbial_np
    this%SOM_np                        = SOM_np
    this%SOM_cn                        = SOM_cn
    this%microbial_nue                 = microbial_nue
    this%microbial_pue                 = microbial_pue
    this%microbial_cue                 = microbial_cue
    this%microbial_cue_max             = microbial_cue_max
    this%microbial_cue_min             = microbial_cue_min
    this%microbial_cn_max              = microbial_cn_max
    this%microbial_cn_min              = microbial_cn_min
    this%f_mic_cn                      = f_mic_cn
    this%frac_woody2dec_litter         = frac_woody2dec_litter
    this%frac_litter2fast_som          = frac_litter2fast_som
    this%frac_fast2slow_som            = frac_fast2slow_som
    this%frac_slow2fast_som            = frac_slow2fast_som
    this%fc_woody2polymeric            = fc_woody2polymeric
    this%fc_soluable2dom               = fc_soluable2dom
    this%k_mort_microbial              = k_mort_microbial
    this%min_mic_biomass               = min_mic_biomass
    this%vmax_mic_uptake               = vmax_mic_uptake
    this%km_mic_uptake                 = km_mic_uptake
    this%km_depolymerisation           = km_depolymerisation
    this%vmax_depolymerisation_residue = vmax_depolymerisation_residue
    this%vmax_depolymerisation_litter  = vmax_depolymerisation_litter
    this%f_soluable_residue            = f_soluable_residue
    this%f_soluable_leaf               = f_soluable_leaf
    this%f_soluable_fine_root          = f_soluable_fine_root
    this%f_soluable_fruit              = f_soluable_fruit
    this%f_soluable_seed_bed           = f_soluable_seed_bed
    this%fc_soluable_max               = fc_soluable_max
    this%k_fc_soluable                 = k_fc_soluable
    this%k_fn_soluable                 = k_fn_soluable
    this%k_fp_soluable                 = k_fp_soluable
    this%lc_leaf_max                   = lc_leaf_max
    this%lc_leaf2sla                   = lc_leaf2sla
    this%lc_fine_root                  = lc_fine_root
    this%lc_coarse_root                = lc_coarse_root
    this%lc_woody_litter               = lc_woody_litter
    this%lc_fruit                      = lc_fruit
    this%lc_seed_bed                   = lc_seed_bed
    this%mycorrhiza_cue                = mycorrhiza_cue
    this%mycorrhiza_cn                 = mycorrhiza_cn
    this%mycorrhiza_np                 = mycorrhiza_np
    this%f_maxresp_myc                 = f_maxresp_myc
    this%f_surface_myc2fineroot        = f_surface_myc2fineroot
    this%km1_up_nh4                    = km1_up_nh4
    this%km2_up_nh4                    = km2_up_nh4
    this%km1_up_no3                    = km1_up_no3
    this%km2_up_no3                    = km2_up_no3
    this%km1_up_po4                    = km1_up_po4
    this%km2_up_po4                    = km2_up_po4
    this%vmax_mic_uptake_nh4           = vmax_mic_uptake_nh4
    this%km_mic_uptake_nh4             = km_mic_uptake_nh4
    this%vmax_mic_uptake_no3           = vmax_mic_uptake_no3
    this%km_mic_uptake_no3             = km_mic_uptake_no3
    this%f_mic_carrier_n               = f_mic_carrier_n
    this%vmax_plant_uptake_nh4         = vmax_plant_uptake_nh4
    this%km_plant_uptake_nh4           = km_plant_uptake_nh4
    this%vmax_plant_uptake_no3         = vmax_plant_uptake_no3
    this%km_plant_uptake_no3           = km_plant_uptake_no3
    this%f_plant_carrier_n             = f_plant_carrier_n
    this%vmax_mic_uptake_po4           = vmax_mic_uptake_po4
    this%km_mic_uptake_po4             = km_mic_uptake_po4
    this%vmax_plant_uptake_po4         = vmax_plant_uptake_po4
    this%km_plant_uptake_po4           = km_plant_uptake_po4
    this%k_mic_uptake_nh4              = k_mic_uptake_nh4
    this%k_plant_uptake_nh4            = k_plant_uptake_nh4
    this%k_mic_uptake_no3              = k_mic_uptake_no3
    this%k_plant_uptake_no3            = k_plant_uptake_no3
    this%k_mic_uptake_po4              = k_mic_uptake_po4
    this%k_plant_uptake_po4            = k_plant_uptake_po4
    this%f_mic_carrier_po4             = f_mic_carrier_po4
    this%f_plant_carrier_po4           = f_plant_carrier_po4
    this%qmax_po4_mineral              = qmax_po4_mineral
    this%km_adsorpt_mineral_po4        = km_adsorpt_mineral_po4
    this%km_adsorpt_OM_po4             = km_adsorpt_OM_po4
    this%ea_nitrification              = ea_nitrification
    this%ed_nitrification              = ed_nitrification
    this%t_opt_nitrification           = t_opt_nitrification
    this%f_nit_noy                     = f_nit_noy
    this%f_nit_n2o                     = f_nit_n2o
    this%ea_denitrification            = ea_denitrification
    this%vmax_nitrification            = vmax_nitrification
    this%vmax_denitrification          = vmax_denitrification
    this%km_denitrification_c          = km_denitrification_c
    this%km_denitrification_no3        = km_denitrification_no3
    this%f_denit_noy                   = f_denit_noy
    this%f_denit_n2o                   = f_denit_n2o
    this%ea_gasdiffusion               = ea_gasdiffusion
    this%lambda_anvf                   = lambda_anvf
    this%k_anvf                        = k_anvf
    this%floss_nmin                    = floss_nmin
    this%fnit_nmin                     = fnit_nmin
    this%fleach_nh4                    = fleach_nh4
    this%fleach_po4                    = fleach_po4
    this%soil_n_threshold_bnf          = soil_n_threshold_bnf
    this%standard_bnf_max_rate         = standard_bnf_max_rate
    this%eta_mic_uptake_nh4            = eta_mic_uptake_nh4
    this%eta_mic_uptake_no3            = eta_mic_uptake_no3
    this%eta_plant_uptake_nh4          = eta_plant_uptake_nh4
    this%eta_plant_uptake_no3          = eta_plant_uptake_no3
    this%eta_nitrification             = eta_nitrification
    this%eta_nitrate_production        = eta_nitrate_production
    this%eta_denitrification           = eta_denitrification
    this%eta_ammonification            = eta_ammonification
    this%eta_volatilisation_nh3        = eta_volatilisation_nh3
    this%k_weath_mineral               = k_weath_mineral
    this%eq_sorp_mineral_po4           = eq_sorp_mineral_po4
    this%k_occlude_po4                 = k_occlude_po4
    this%vmax_biochem_po4              = vmax_biochem_po4
    this%km_mic_biochem_po4            = km_mic_biochem_po4
    this%km_dip_biochem_po4            = km_dip_biochem_po4
    this%km_rootc_biochem_po4          = km_rootc_biochem_po4
    this%km_rootc_weath_po4            = km_rootc_weath_po4
    this%biomin_thres_PC               = biomin_thres_PC
    this%fc_p_res2dom                  = fc_p_res2dom
    this%fc_n_res2dom                  = fc_n_res2dom

  END SUBROUTINE read_parameter_slm_list


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded SLM parameter values with values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_parameter_slm_values_with_list_values(this, proportional_parameter_changes)

    USE mo_base_constants,                  ONLY: one_day, one_year, dtime
    USE mo_assimi_constants
    USE mo_pheno_constants
    !USE mo_atmland_constants
    USE mo_rad_constants
    USE mo_jsb_phys_constants,              ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero
    USE mo_quincy_update_time_averages      ! includes the mavg_period_* parameters
    USE mo_sse_constants
    USE mo_sb_constants
    USE mo_veg_constants

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_parameters_slm_from_file), INTENT(inout) :: this
    LOGICAL                          , INTENT(in)    :: proportional_parameter_changes
    ! ---------------------------
    ! 0.2 Local
    LOGICAL                     :: dependency_recalculation   ! set TRUE for parameters that need to be re-calculated because a parameter it depends on got a new value
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_parameter_slm_values_with_list_values'


    !> 1.0 read parameter values from file
    !! 
    CALL this%read_param

    !> 2.0 replace parameter values from mo_PROC_constants for parameters provided by parameter_slm_run.list
    !! 
    !! dependency check: IF ((ABS(this%x    - missval) > eps8) .OR. &
    !!                       (ABS(this%y    - missval) > eps8)) THEN
    !!                       dependency_recalculation = .TRUE.
    !!                   ENDIF
    !! 
    !! routine:  SUBROUTINE Set_quincy_parameter_value(nml_abs_val, nml_prop_chng, proportional, dependency, def_val)
    !!    CALL Set_quincy_parameter_value(this%xyz, &         ! or "this%xyz */+- any other parameter"
    !!                                    this%xyz, &         ! just the parameter as defined in the namelist parameter_slm_run.list
    !!                                    "parameter name", proportional_parameter_changes, dependency_recalculation, &
    !!                                    xyz)                ! the original/default parameter variable as declared in mo_PROC_constants
    !! 
    ! CALL Set_quincy_parameter_value(this%xyz, &
    !                                 this%xyz, &
    !                                 "xyz", proportional_parameter_changes, dependency_recalculation, &
    !                                 xyz)
    ! 
    ! CALL Set_quincy_parameter_value(this%xyz, &
    !                                 this%xyz, &
    !                                 "xyz", proportional_parameter_changes, dependency_recalculation, &
    !                                 xyz)
    ! 
    ! CALL Set_quincy_parameter_value(this%xyz, &
    !                                 this%xyz, &
    !                                 "xyz", proportional_parameter_changes, dependency_recalculation, &
    !                                 xyz)
    !! 
    ! ASSIMILATION_CTL
    CALL Set_quincy_parameter_value(this%jmax2n / 4.0_wp, &
                                    this%jmax2n, &
                                    "jmax2n", proportional_parameter_changes, dependency_recalculation, &
                                    jmax2n)
                                    
    CALL Set_quincy_parameter_value(this%vcmax2n, &
                                    this%vcmax2n, &
                                    "vcmax2n", proportional_parameter_changes, dependency_recalculation, &
                                    vcmax2n)
                                    
    CALL Set_quincy_parameter_value(this%pepc2n, &
                                    this%pepc2n, &
                                    "pepc2n", proportional_parameter_changes, dependency_recalculation, &
                                    pepc2n)
                                    
    CALL Set_quincy_parameter_value(this%chl2n / 39.8_wp, &
                                    this%chl2n, &
                                    "chl2n", proportional_parameter_changes, dependency_recalculation, &
                                    chl2n)
                                    
    CALL Set_quincy_parameter_value(this%k1_fn_struc * molar_mass_N / 1000000._wp, &
                                    this%k1_fn_struc, &
                                    "k1_fn_struc", proportional_parameter_changes, dependency_recalculation, &
                                    k1_fn_struc)
                                    
    CALL Set_quincy_parameter_value(this%E0kc, &
                                    this%E0kc, &
                                    "E0kc", proportional_parameter_changes, dependency_recalculation, &
                                    E0kc)
                                    
    CALL Set_quincy_parameter_value(this%E1kc, &
                                    this%E1kc, &
                                    "E1kc", proportional_parameter_changes, dependency_recalculation, &
                                    E1kc)
                                    
    CALL Set_quincy_parameter_value(this%E0ko, &
                                    this%E0ko, &
                                    "E0ko", proportional_parameter_changes, dependency_recalculation, &
                                    E0ko)
                                    
    CALL Set_quincy_parameter_value(this%E1ko, &
                                    this%E1ko, &
                                    "E1ko", proportional_parameter_changes, dependency_recalculation, &
                                    E1ko)
                                    
    CALL Set_quincy_parameter_value(this%E0pcp, &
                                    this%E0pcp, &
                                    "E0pcp", proportional_parameter_changes, dependency_recalculation, &
                                    E0pcp)
                                    
    CALL Set_quincy_parameter_value(this%E1pcp, &
                                    this%E1pcp, &
                                    "E1pcp", proportional_parameter_changes, dependency_recalculation, &
                                    E1pcp)
                                    
    CALL Set_quincy_parameter_value(this%E0v, &
                                    this%E0v, &
                                    "E0v", proportional_parameter_changes, dependency_recalculation, &
                                    E0v)
                                    
    CALL Set_quincy_parameter_value(this%E1v, &
                                    this%E1v, &
                                    "E1v", proportional_parameter_changes, dependency_recalculation, &
                                    E1v)
                                    
    CALL Set_quincy_parameter_value(this%t_jmax_offset, &
                                    this%t_jmax_offset, &
                                    "t_jmax_offset", proportional_parameter_changes, dependency_recalculation, &
                                    t_jmax_offset)
                                    
    CALL Set_quincy_parameter_value(this%t_jmax_slope, &
                                    this%t_jmax_slope, &
                                    "t_jmax_slope", proportional_parameter_changes, dependency_recalculation, &
                                    t_jmax_slope)
                                    
    IF (ABS(this%t_jmax_offset - missval) > eps8) THEN
        dependency_recalculation = .TRUE.
        t_jmax_opt_min           = t_jmax_offset   ! copied from mo_assimi_constants
    ENDIF
    CALL Set_quincy_parameter_value(t_jmax_offset, &   ! default: t_jmax_opt_min = t_jmax_offset | hence, only proportional changes may have an effect
                                    this%t_jmax_opt_min, &
                                    "t_jmax_opt_min", proportional_parameter_changes, dependency_recalculation, &
                                    t_jmax_opt_min)
                                    
    CALL Set_quincy_parameter_value(this%t_jmax_opt_max, &
                                    this%t_jmax_opt_max, &
                                    "t_jmax_opt_max", proportional_parameter_changes, dependency_recalculation, &
                                    t_jmax_opt_max)
                                    
    CALL Set_quincy_parameter_value(this%Tref_pepc, &
                                    this%Tref_pepc, &
                                    "Tref_pepc", proportional_parameter_changes, dependency_recalculation, &
                                    Tref_pepc)
                                    
    CALL Set_quincy_parameter_value(this%Tbase_pepc, &
                                    this%Tbase_pepc, &
                                    "Tbase_pepc", proportional_parameter_changes, dependency_recalculation, &
                                    Tbase_pepc)
                                    
    CALL Set_quincy_parameter_value(this%alpha_i, &
                                    this%alpha_i, &
                                    "alpha_i", proportional_parameter_changes, dependency_recalculation, &
                                    alpha_i)
                                    
    CALL Set_quincy_parameter_value(this%pO2, &
                                    this%pO2, &
                                    "pO2", proportional_parameter_changes, dependency_recalculation, &
                                    pO2)
                                    
    CALL Set_quincy_parameter_value(this%CiCa_default_C3, &
                                    this%CiCa_default_C3, &
                                    "CiCa_default_C3", proportional_parameter_changes, dependency_recalculation, &
                                    CiCa_default_C3)
                                    
    CALL Set_quincy_parameter_value(this%CiCa_default_C4, &
                                    this%CiCa_default_C4, &
                                    "CiCa_default_C4", proportional_parameter_changes, dependency_recalculation, &
                                    CiCa_default_C4)
                                    
    CALL Set_quincy_parameter_value(this%ps_it_max, &
                                    this%ps_it_max, &
                                    "ps_it_max", proportional_parameter_changes, dependency_recalculation, &
                                    ps_it_max)
                                    
    CALL Set_quincy_parameter_value(this%ci_max, &
                                    this%ci_max, &
                                    "ci_max", proportional_parameter_changes, dependency_recalculation, &
                                    ci_max)
                                    
    ! current assumption: parameters from mo_jsb_phys_constants are not modified
    ! IF ((ABS(this%Dwv  - missval) > eps8) .OR. &
    !     (ABS(this%Dco2 - missval) > eps8)) THEN
    !     dependency_recalculation = .TRUE.
    ! ENDIF
    CALL Set_quincy_parameter_value(Dwv / Dco2, &                       ! default: Dwv2co2_air = Dwv / Dco2 | hence, only proportional changes may have an effect
                                    this%Dwv2co2_air, &
                                    "Dwv2co2_air", proportional_parameter_changes, dependency_recalculation, &
                                    Dwv2co2_air)

    ! current assumption: parameters from mo_jsb_phys_constants are not modified
    ! IF ((ABS(this%Dwv  - missval) > eps8) .OR. &
    !     (ABS(this%Dco2 - missval) > eps8)) THEN
    !     dependency_recalculation = .TRUE.
    ! ENDIF
    CALL Set_quincy_parameter_value(( Dwv / Dco2  ) ** (2._wp/3._wp), & ! default: Dwv2co2_turb = ( Dwv / Dco2  ) ** (2._wp/3._wp) | hence, only proportional changes may have an effect
                                    this%Dwv2co2_turb, &
                                    "Dwv2co2_turb", proportional_parameter_changes, dependency_recalculation, &
                                    Dwv2co2_turb)
                                    
    CALL Set_quincy_parameter_value(this%discr_ps_a_C13, &
                                    this%discr_ps_a_C13, &
                                    "discr_ps_a_C13", proportional_parameter_changes, dependency_recalculation, &
                                    discr_ps_a_C13)
                                    
    CALL Set_quincy_parameter_value(this%discr_ps_b_C13, &
                                    this%discr_ps_b_C13, &
                                    "discr_ps_b_C13", proportional_parameter_changes, dependency_recalculation, &
                                    discr_ps_b_C13)
                                    
    CALL Set_quincy_parameter_value(this%discr_ps_c_C13, &
                                    this%discr_ps_c_C13, &
                                    "discr_ps_c_C13", proportional_parameter_changes, dependency_recalculation, &
                                    discr_ps_c_C13)

    IF (ABS(this%discr_ps_a_C13 - missval) > eps8) THEN
        dependency_recalculation = .TRUE.
        discr_ps_a_C14           = 1.97_wp * discr_ps_a_C13   ! copied from mo_assimi_constants
    ENDIF
    CALL Set_quincy_parameter_value(this%discr_ps_a_C14 * discr_ps_a_C13, &
                                    this%discr_ps_a_C14, &
                                    "discr_ps_a_C14", proportional_parameter_changes, dependency_recalculation, &
                                    discr_ps_a_C14)

    IF (ABS(this%discr_ps_b_C13 - missval) > eps8) THEN
        dependency_recalculation = .TRUE.
        discr_ps_b_C14           = 1.89_wp * discr_ps_b_C13   ! copied from mo_assimi_constants
    ENDIF
    CALL Set_quincy_parameter_value(this%discr_ps_b_C14 * discr_ps_b_C13, &
                                    this%discr_ps_b_C14, &
                                    "discr_ps_b_C14", proportional_parameter_changes, dependency_recalculation, &
                                    discr_ps_b_C14)
                                    
    IF (ABS(this%discr_ps_c_C13 - missval) > eps8) THEN
        dependency_recalculation = .TRUE.
        discr_ps_c_C14           = 1.89_wp * discr_ps_c_C13   ! copied from mo_assimi_constants
    ENDIF
    CALL Set_quincy_parameter_value(this%discr_ps_c_C14 * discr_ps_c_C13, &
                                    this%discr_ps_c_C14, &
                                    "discr_ps_c_C14", proportional_parameter_changes, dependency_recalculation, &
                                    discr_ps_c_C14)
                                    
    CALL Set_quincy_parameter_value(this%discr_ps_phi, &
                                    this%discr_ps_phi, &
                                    "discr_ps_phi", proportional_parameter_changes, dependency_recalculation, &
                                    discr_ps_phi)
                                    
    CALL Set_quincy_parameter_value(this%ka, &
                                    this%ka, &
                                    "ka", proportional_parameter_changes, dependency_recalculation, &
                                    ka)
                                    
    CALL Set_quincy_parameter_value(this%kn, &
                                    this%kn, &
                                    "kn", proportional_parameter_changes, dependency_recalculation, &
                                    kn)
                                    
    CALL Set_quincy_parameter_value(this%soa_min_temp, &
                                    this%soa_min_temp, &
                                    "soa_min_temp", proportional_parameter_changes, dependency_recalculation, &
                                    soa_min_temp)
                                    
    CALL Set_quincy_parameter_value(this%soa_max_temp, &
                                    this%soa_max_temp, &
                                    "soa_max_temp", proportional_parameter_changes, dependency_recalculation, &
                                    soa_max_temp)

    ! PHENOLOGY_CTL
    CALL Set_quincy_parameter_value(this%gdd_t_air_threshold + Tzero, &
                                    this%gdd_t_air_threshold, &
                                    "gdd_t_air_threshold", proportional_parameter_changes, dependency_recalculation, &
                                    gdd_t_air_threshold)

    ! RADIATION_CTL
    CALL Set_quincy_parameter_value(this%kbl0_vis, &
                                    this%kbl0_vis, &
                                    "kbl0_vis", proportional_parameter_changes, dependency_recalculation, &
                                    kbl0_vis)
    
    CALL Set_quincy_parameter_value(this%kdf0_vis, &
                                    this%kdf0_vis, &
                                    "kdf0_vis", proportional_parameter_changes, dependency_recalculation, &
                                    kdf0_vis)
    
    CALL Set_quincy_parameter_value(this%kbl0_nir, &
                                    this%kbl0_nir, &
                                    "kbl0_nir", proportional_parameter_changes, dependency_recalculation, &
                                    kbl0_nir)

    CALL Set_quincy_parameter_value(this%kdf0_nir, &
                                    this%kdf0_nir, &
                                    "kdf0_nir", proportional_parameter_changes, dependency_recalculation, &
                                    kdf0_nir)
    
    CALL Set_quincy_parameter_value(this%rho2sbeta, &
                                    this%rho2sbeta, &
                                    "rho2sbeta", proportional_parameter_changes, dependency_recalculation, &
                                    rho2sbeta)
    
    CALL Set_quincy_parameter_value(this%k_csf, &
                                    this%k_csf, &
                                    "k_csf", proportional_parameter_changes, dependency_recalculation, &
                                    k_csf)

    CALL Set_quincy_parameter_value(this%min_cos_zenith_angle_rad, &
                                    this%min_cos_zenith_angle_rad, &
                                    "min_cos_zenith_angle_rad", proportional_parameter_changes, dependency_recalculation, &
                                    min_cos_zenith_angle_rad)
    
    CALL Set_quincy_parameter_value(this%def_alb_vis_soil, &
                                    this%def_alb_vis_soil, &
                                    "def_alb_vis_soil", proportional_parameter_changes, dependency_recalculation, &
                                    def_alb_vis_soil)
    
    CALL Set_quincy_parameter_value(this%def_alb_nir_soil, &
                                    this%def_alb_nir_soil, &
                                    "def_alb_nir_soil", proportional_parameter_changes, dependency_recalculation, &
                                    def_alb_nir_soil)
    
    CALL Set_quincy_parameter_value(this%rad2ppfd, &
                                    this%rad2ppfd, &
                                    "rad2ppfd", proportional_parameter_changes, dependency_recalculation, &
                                    rad2ppfd)
    
    CALL Set_quincy_parameter_value(this%rfr_ratio_toc, &
                                    this%rfr_ratio_toc, &
                                    "rfr_ratio_toc", proportional_parameter_changes, dependency_recalculation, &
                                    rfr_ratio_toc)
    
    CALL Set_quincy_parameter_value(this%k_r2fr_chl - 0.012_wp, &
                                    this%k_r2fr_chl, &
                                    "k_r2fr_chl", proportional_parameter_changes, dependency_recalculation, &
                                    k_r2fr_chl)
    
    ! SHARED_CTL
    CALL Set_quincy_parameter_value(this%mavg_period_weekly, &
                                    this%mavg_period_weekly, &
                                    "mavg_period_weekly", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_weekly)
    
    CALL Set_quincy_parameter_value(this%mavg_period_monthly, &
                                    this%mavg_period_monthly, &
                                    "mavg_period_monthly", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_monthly)
    
    CALL Set_quincy_parameter_value(this%mavg_period_yearly, &
                                    this%mavg_period_yearly, &
                                    "mavg_period_yearly", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_yearly)

    CALL Set_quincy_parameter_value(this%mavg_period_tsoa, &
                                    this%mavg_period_tsoa, &
                                    "mavg_period_tsoa", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tsoa)
    
    CALL Set_quincy_parameter_value(this%mavg_period_tuptake, &
                                    this%mavg_period_tuptake, &
                                    "mavg_period_tuptake", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tuptake)
    
    CALL Set_quincy_parameter_value(this%mavg_period_tlabile, &
                                    this%mavg_period_tlabile, &
                                    "mavg_period_tlabile", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tlabile)

    CALL Set_quincy_parameter_value(this%mavg_period_tphen, &
                                    this%mavg_period_tphen, &
                                    "mavg_period_tphen", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tphen)

    CALL Set_quincy_parameter_value(this%mavg_period_tenzyme, &
                                    this%mavg_period_tenzyme, &
                                    "mavg_period_tenzyme", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tenzyme)
    
    CALL Set_quincy_parameter_value(this%mavg_period_tfrac, &
                                    this%mavg_period_tfrac, &
                                    "mavg_period_tfrac", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tfrac)
    
    CALL Set_quincy_parameter_value(this%mavg_period_tcnl, &
                                    this%mavg_period_tcnl, &
                                    "mavg_period_tcnl", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tcnl)

    CALL Set_quincy_parameter_value(this%mavg_period_talloc, &
                                    this%mavg_period_talloc, &
                                    "mavg_period_talloc", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_talloc)
    
    CALL Set_quincy_parameter_value(this%mavg_period_talloc_dynamic * 365.0_wp, &  ! default value: 30.0_wp * 365.0_wp (i.e., 30 years)
                                    this%mavg_period_talloc_dynamic, &
                                    "mavg_period_talloc_dynamic", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_talloc_dynamic)
    
    CALL Set_quincy_parameter_value(this%mavg_period_tacclim, &
                                    this%mavg_period_tacclim, &
                                    "mavg_period_tacclim", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tacclim)

    CALL Set_quincy_parameter_value(this%mavg_period_tmic, &
                                    this%mavg_period_tmic, &
                                    "mavg_period_tmic", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tmic)
    
    CALL Set_quincy_parameter_value(this%mavg_period_tgrowth, &
                                    this%mavg_period_tgrowth, &
                                    "mavg_period_tgrowth", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tgrowth)
    
    CALL Set_quincy_parameter_value(this%mavg_period_tvegdyn, &
                                    this%mavg_period_tvegdyn, &
                                    "mavg_period_tvegdyn", proportional_parameter_changes, dependency_recalculation, &
                                    mavg_period_tvegdyn)

    ! SSE_CTL
    CALL Set_quincy_parameter_value(this%soil_heat_cap, &
                                    this%soil_heat_cap, &
                                    "soil_heat_cap", proportional_parameter_changes, dependency_recalculation, &
                                    soil_heat_cap)

    CALL Set_quincy_parameter_value(this%soil_therm_cond, &
                                    this%soil_therm_cond, &
                                    "soil_therm_cond", proportional_parameter_changes, dependency_recalculation, &
                                    soil_therm_cond)
    
    CALL Set_quincy_parameter_value(this%InterceptionEfficiency, &
                                    this%InterceptionEfficiency, &
                                    "InterceptionEfficiency", proportional_parameter_changes, dependency_recalculation, &
                                    InterceptionEfficiency)
    
    CALL Set_quincy_parameter_value(this%w_skin_max, &
                                    this%w_skin_max, &
                                    "w_skin_max", proportional_parameter_changes, dependency_recalculation, &
                                    w_skin_max)

    CALL Set_quincy_parameter_value(this%drain_frac, &
                                    this%drain_frac, &
                                    "drain_frac", proportional_parameter_changes, dependency_recalculation, &
                                    drain_frac)

    CALL Set_quincy_parameter_value(this%fdrain_runoff, &
                                    this%fdrain_runoff, &
                                    "fdrain_runoff", proportional_parameter_changes, dependency_recalculation, &
                                    fdrain_runoff)
    
    CALL Set_quincy_parameter_value(this%wpot_pwp, &
                                    this%wpot_pwp, &
                                    "wpot_pwp", proportional_parameter_changes, dependency_recalculation, &
                                    wpot_pwp)

    CALL Set_quincy_parameter_value(this%wpot_fc, &
                                    this%wpot_fc, &
                                    "wpot_fc", proportional_parameter_changes, dependency_recalculation, &
                                    wpot_fc)

    CALL Set_quincy_parameter_value(this%wpot_min, &
                                    this%wpot_min, &
                                    "wpot_min", proportional_parameter_changes, dependency_recalculation, &
                                    wpot_min)
    
    CALL Set_quincy_parameter_value(this%k_pwp_a, &
                                    this%k_pwp_a, &
                                    "k_pwp_a", proportional_parameter_changes, dependency_recalculation, &
                                    k_pwp_a)

    CALL Set_quincy_parameter_value(this%k_pwp_c, &
                                    this%k_pwp_c, &
                                    "k_pwp_c", proportional_parameter_changes, dependency_recalculation, &
                                    k_pwp_c)

    CALL Set_quincy_parameter_value(this%k_pwp_s, &
                                    this%k_pwp_s, &
                                    "k_pwp_s", proportional_parameter_changes, dependency_recalculation, &
                                    k_pwp_s)
    
    CALL Set_quincy_parameter_value(this%k_pwp_sc, &
                                    this%k_pwp_sc, &
                                    "k_pwp_sc", proportional_parameter_changes, dependency_recalculation, &
                                    k_pwp_sc)

    CALL Set_quincy_parameter_value(this%k_pwp_at, &
                                    this%k_pwp_at, &
                                    "k_pwp_at", proportional_parameter_changes, dependency_recalculation, &
                                    k_pwp_at)

    CALL Set_quincy_parameter_value(this%k_pwp_bt, &
                                    this%k_pwp_bt, &
                                    "k_pwp_bt", proportional_parameter_changes, dependency_recalculation, &
                                    k_pwp_bt)
    
    CALL Set_quincy_parameter_value(this%k_fc_a, &
                                    this%k_fc_a, &
                                    "k_fc_a", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_a)

    CALL Set_quincy_parameter_value(this%k_fc_c, &
                                    this%k_fc_c, &
                                    "k_fc_c", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_c)

    CALL Set_quincy_parameter_value(this%k_fc_s, &
                                    this%k_fc_s, &
                                    "k_fc_s", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_s)
    
    CALL Set_quincy_parameter_value(this%k_fc_sc, &
                                    this%k_fc_sc, &
                                    "k_fc_sc", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_sc)

    CALL Set_quincy_parameter_value(this%k_fc_at, &
                                    this%k_fc_at, &
                                    "k_fc_at", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_at)

    CALL Set_quincy_parameter_value(this%k_fc_bt, &
                                    this%k_fc_bt, &
                                    "k_fc_bt", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_bt)
    
    CALL Set_quincy_parameter_value(this%k_fc_ct, &
                                    this%k_fc_ct, &
                                    "k_fc_ct", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_ct)

    CALL Set_quincy_parameter_value(this%k_sat_a, &
                                    this%k_sat_a, &
                                    "k_sat_a", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_a)

    CALL Set_quincy_parameter_value(this%k_sat_c, &
                                    this%k_sat_c, &
                                    "k_sat_c", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_c)
    
    CALL Set_quincy_parameter_value(this%k_sat_s, &
                                    this%k_sat_s, &
                                    "k_sat_s", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_s)

    CALL Set_quincy_parameter_value(this%k_sat_sc, &
                                    this%k_sat_sc, &
                                    "k_sat_sc", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_sc)

    CALL Set_quincy_parameter_value(this%k_sat_at, &
                                    this%k_sat_at, &
                                    "k_sat_at", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_at)
    
    CALL Set_quincy_parameter_value(this%k_sat_bt, &
                                    this%k_sat_bt, &
                                    "k_sat_bt", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_bt)

    CALL Set_quincy_parameter_value(this%k_sat_ct, &
                                    this%k_sat_ct, &
                                    "k_sat_ct", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_ct)

    CALL Set_quincy_parameter_value(this%k_sat_dt, &
                                    this%k_sat_dt, &
                                    "k_sat_dt", proportional_parameter_changes, dependency_recalculation, &
                                    k_sat_dt)

    CALL Set_quincy_parameter_value(this%kdiff_sat_max  / 1000.0_wp / 3600._wp, &
                                    this%kdiff_sat_max, &
                                    "kdiff_sat_max", proportional_parameter_changes, dependency_recalculation, &
                                    kdiff_sat_max)

    ! VEGETATION_CTL
    CALL Set_quincy_parameter_value(this%fresp_growth, &
                                    this%fresp_growth, &
                                    "fresp_growth", proportional_parameter_changes, dependency_recalculation, &
                                    fresp_growth)
    
    CALL Set_quincy_parameter_value(this%fmaint_rate_base, &
                                    this%fmaint_rate_base, &
                                    "fmaint_rate_base", proportional_parameter_changes, dependency_recalculation, &
                                    fmaint_rate_base)
    
    IF (ABS(this%fmaint_rate_base - missval) > eps8) THEN
        dependency_recalculation = .TRUE.
        fmaint_rate_wood         = 0.25_wp * fmaint_rate_base   ! copied from mo_veg_constants
    ENDIF
    CALL Set_quincy_parameter_value(this%fmaint_rate_wood  * fmaint_rate_base, &
                                    this%fmaint_rate_wood, &
                                    "fmaint_rate_wood", proportional_parameter_changes, dependency_recalculation, &
                                    fmaint_rate_wood)

    CALL Set_quincy_parameter_value(this%t_maint_k1, &
                                    this%t_maint_k1, &
                                    "t_maint_k1", proportional_parameter_changes, dependency_recalculation, &
                                    t_maint_k1)
    
    CALL Set_quincy_parameter_value(this%t_maint_k2, &
                                    this%t_maint_k2, &
                                    "t_maint_k2", proportional_parameter_changes, dependency_recalculation, &
                                    t_maint_k2)
    
    CALL Set_quincy_parameter_value(this%t_maint_k3, &
                                    this%t_maint_k3, &
                                    "t_maint_k3", proportional_parameter_changes, dependency_recalculation, &
                                    t_maint_k3)

    CALL Set_quincy_parameter_value(this%f_resp_acclim, &
                                    this%f_resp_acclim, &
                                    "f_resp_acclim", proportional_parameter_changes, dependency_recalculation, &
                                    f_resp_acclim)
    
    CALL Set_quincy_parameter_value(this%t_acclim_zero, &
                                    this%t_acclim_zero, &
                                    "t_acclim_zero", proportional_parameter_changes, dependency_recalculation, &
                                    t_acclim_zero)
    
    CALL Set_quincy_parameter_value(this%k0_fn_pepc_C3, &
                                    this%k0_fn_pepc_C3, &
                                    "k0_fn_pepc_C3", proportional_parameter_changes, dependency_recalculation, &
                                    k0_fn_pepc_C3)

    CALL Set_quincy_parameter_value(this%k0_fn_pepc_C4, &
                                    this%k0_fn_pepc_C4, &
                                    "k0_fn_pepc_C4", proportional_parameter_changes, dependency_recalculation, &
                                    k0_fn_pepc_C4)

    CALL Set_quincy_parameter_value(this%jmax2vcmax_C3 / 4._wp, &
                                    this%jmax2vcmax_C3, &
                                    "jmax2vcmax_C3", proportional_parameter_changes, dependency_recalculation, &
                                    jmax2vcmax_C3)
    
    CALL Set_quincy_parameter_value(this%jmax2vcmax_C4, &
                                    this%jmax2vcmax_C4, &
                                    "jmax2vcmax_C4", proportional_parameter_changes, dependency_recalculation, &
                                    jmax2vcmax_C4)
    
    CALL Set_quincy_parameter_value(this%tau_labile, &
                                    this%tau_labile, &
                                    "tau_labile", proportional_parameter_changes, dependency_recalculation, &
                                    tau_labile)

    CALL Set_quincy_parameter_value(this%lambda_labile_temp, &
                                    this%lambda_labile_temp, &
                                    "lambda_labile_temp", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_labile_temp)
    
    CALL Set_quincy_parameter_value(this%k_labile_temp, &
                                    this%k_labile_temp, &
                                    "k_labile_temp", proportional_parameter_changes, dependency_recalculation, &
                                    k_labile_temp)
    
    CALL Set_quincy_parameter_value(this%lambda_labile_theta, &
                                    this%lambda_labile_theta, &
                                    "lambda_labile_theta", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_labile_theta)

    CALL Set_quincy_parameter_value(this%k_labile_theta, &
                                    this%k_labile_theta, &
                                    "k_labile_theta", proportional_parameter_changes, dependency_recalculation, &
                                    k_labile_theta)
    
    CALL Set_quincy_parameter_value(this%lambda_phiphen, &
                                    this%lambda_phiphen, &
                                    "lambda_phiphen", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_phiphen)
    
    CALL Set_quincy_parameter_value(this%k_phiphen, &
                                    this%k_phiphen, &
                                    "k_phiphen", proportional_parameter_changes, dependency_recalculation, &
                                    k_phiphen)

    CALL Set_quincy_parameter_value(this%lambda_phimaint, &
                                    this%lambda_phimaint, &
                                    "lambda_phimaint", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_phimaint)

    CALL Set_quincy_parameter_value(this%k_phimaint, &
                                    this%k_phimaint, &
                                    "k_phimaint", proportional_parameter_changes, dependency_recalculation, &
                                    k_phimaint)
    
    CALL Set_quincy_parameter_value(this%lambda_phimaint_nut, &
                                    this%lambda_phimaint_nut, &
                                    "lambda_phimaint_nut", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_phimaint_nut)
    
    CALL Set_quincy_parameter_value(this%k_phimaint_nut, &
                                    this%k_phimaint_nut, &
                                    "k_phimaint_nut", proportional_parameter_changes, dependency_recalculation, &
                                    k_phimaint_nut)

    CALL Set_quincy_parameter_value(this%lambda_phistore, &
                                    this%lambda_phistore, &
                                    "lambda_phistore", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_phistore)
    
    CALL Set_quincy_parameter_value(this%k_phistore, &
                                    this%k_phistore, &
                                    "k_phistore", proportional_parameter_changes, dependency_recalculation, &
                                    k_phistore)
    
    CALL Set_quincy_parameter_value(this%knut_labile, &
                                    this%knut_labile, &
                                    "knut_labile", proportional_parameter_changes, dependency_recalculation, &
                                    knut_labile)

    CALL Set_quincy_parameter_value(this%lambda_sinklim_ps, &
                                    this%lambda_sinklim_ps, &
                                    "lambda_sinklim_ps", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_sinklim_ps)
    
    CALL Set_quincy_parameter_value(this%k_sinklim_ps, &
                                    this%k_sinklim_ps, &
                                    "k_sinklim_ps", proportional_parameter_changes, dependency_recalculation, &
                                    k_sinklim_ps)
    
    CALL Set_quincy_parameter_value(this%k_cnp_sinklim_ps, &
                                    this%k_cnp_sinklim_ps, &
                                    "k_cnp_sinklim_ps", proportional_parameter_changes, dependency_recalculation, &
                                    k_cnp_sinklim_ps)

    CALL Set_quincy_parameter_value(this%beta_sinklim_ps_min, &
                                    this%beta_sinklim_ps_min, &
                                    "beta_sinklim_ps_min", proportional_parameter_changes, dependency_recalculation, &
                                    beta_sinklim_ps_min)

    CALL Set_quincy_parameter_value(this%fstore_leaf_max, &
                                    this%fstore_leaf_max, &
                                    "fstore_leaf_max", proportional_parameter_changes, dependency_recalculation, &
                                    fstore_leaf_max)
    
    CALL Set_quincy_parameter_value(this%fstore_fine_root_max, &
                                    this%fstore_fine_root_max, &
                                    "fstore_fine_root_max", proportional_parameter_changes, dependency_recalculation, &
                                    fstore_fine_root_max)
    
    CALL Set_quincy_parameter_value(this%fstore_sap_wood_max, &
                                    this%fstore_sap_wood_max, &
                                    "fstore_sap_wood_max", proportional_parameter_changes, dependency_recalculation, &
                                    fstore_sap_wood_max)

    CALL Set_quincy_parameter_value(this%k_phi_interact, &
                                    this%k_phi_interact, &
                                    "k_phi_interact", proportional_parameter_changes, dependency_recalculation, &
                                    k_phi_interact)
    
    CALL Set_quincy_parameter_value(this%resorp_fract_leaf, &
                                    this%resorp_fract_leaf, &
                                    "resorp_fract_leaf", proportional_parameter_changes, dependency_recalculation, &
                                    resorp_fract_leaf)
    
    CALL Set_quincy_parameter_value(this%resorp_fract_wood, &
                                    this%resorp_fract_wood, &
                                    "resorp_fract_wood", proportional_parameter_changes, dependency_recalculation, &
                                    resorp_fract_wood)

    CALL Set_quincy_parameter_value(this%root2leaf_cn, &
                                    this%root2leaf_cn, &
                                    "root2leaf_cn", proportional_parameter_changes, dependency_recalculation, &
                                    root2leaf_cn)
    
    CALL Set_quincy_parameter_value(this%wood2leaf_cn, &
                                    this%wood2leaf_cn, &
                                    "wood2leaf_cn", proportional_parameter_changes, dependency_recalculation, &
                                    wood2leaf_cn)
    
    CALL Set_quincy_parameter_value(this%root2leaf_np, &
                                    this%root2leaf_np, &
                                    "root2leaf_np", proportional_parameter_changes, dependency_recalculation, &
                                    root2leaf_np)

    CALL Set_quincy_parameter_value(this%wood2leaf_np, &
                                    this%wood2leaf_np, &
                                    "wood2leaf_np", proportional_parameter_changes, dependency_recalculation, &
                                    wood2leaf_np)

    CALL Set_quincy_parameter_value(this%tau_nutrient_recycling, &
                                    this%tau_nutrient_recycling, &
                                    "tau_nutrient_recycling", proportional_parameter_changes, dependency_recalculation, &
                                    tau_nutrient_recycling)
    
    CALL Set_quincy_parameter_value(this%kappa_f_demand_nc, &
                                    this%kappa_f_demand_nc, &
                                    "kappa_f_demand_nc", proportional_parameter_changes, dependency_recalculation, &
                                    kappa_f_demand_nc)
    
    CALL Set_quincy_parameter_value(this%kappa_f_demand_pn, &
                                    this%kappa_f_demand_pn, &
                                    "kappa_f_demand_pn", proportional_parameter_changes, dependency_recalculation, &
                                    kappa_f_demand_pn)

    CALL Set_quincy_parameter_value(this%k_f_demand, &
                                    this%k_f_demand, &
                                    "k_f_demand", proportional_parameter_changes, dependency_recalculation, &
                                    k_f_demand)
    
    CALL Set_quincy_parameter_value(this%cost_n_fixation_max / (molar_mass_C/molar_mass_N), &
                                    this%cost_n_fixation_max, &
                                    "cost_n_fixation_max", proportional_parameter_changes, dependency_recalculation, &
                                    cost_n_fixation_max)
    
    CALL Set_quincy_parameter_value(this%cost_n_fixation_min / (molar_mass_C/molar_mass_N), &
                                    this%cost_n_fixation_min, &
                                    "cost_n_fixation_min", proportional_parameter_changes, dependency_recalculation, &
                                    cost_n_fixation_min)

    CALL Set_quincy_parameter_value(this%k1_temp_bnf, &
                                    this%k1_temp_bnf, &
                                    "k1_temp_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    k1_temp_bnf)
    
    CALL Set_quincy_parameter_value(this%k2_temp_bnf, &
                                    this%k2_temp_bnf, &
                                    "k2_temp_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    k2_temp_bnf)
    
    CALL Set_quincy_parameter_value(this%k3_temp_bnf, &
                                    this%k3_temp_bnf, &
                                    "k3_temp_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    k3_temp_bnf)

    CALL Set_quincy_parameter_value(this%k4_temp_bnf, &
                                    this%k4_temp_bnf, &
                                    "k4_temp_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    k4_temp_bnf)

    CALL Set_quincy_parameter_value(this%transform_cost_nh4 / (molar_mass_C/molar_mass_N), &
                                    this%transform_cost_nh4, &
                                    "transform_cost_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    transform_cost_nh4)
    
    CALL Set_quincy_parameter_value(this%transform_cost_no3 / (molar_mass_C/molar_mass_N), &
                                    this%transform_cost_no3, &
                                    "transform_cost_no3", proportional_parameter_changes, dependency_recalculation, &
                                    transform_cost_no3)
    
    CALL Set_quincy_parameter_value(this%transform_cost_org / (molar_mass_C/molar_mass_N), &
                                    this%transform_cost_org, &
                                    "transform_cost_org", proportional_parameter_changes, dependency_recalculation, &
                                    transform_cost_org)

    IF ((ABS(this%transform_cost_no3 - missval) > eps8) .OR. & 
        (ABS(this%transform_cost_nh4 - missval) > eps8)) THEN
        dependency_recalculation = .TRUE.
        f_nacq_nh4               = transform_cost_no3 / (transform_cost_no3 + transform_cost_nh4)   ! copied from mo_veg_constants
    ENDIF
    CALL Set_quincy_parameter_value(transform_cost_no3 / (transform_cost_no3 + transform_cost_nh4), &   ! only proportional changes may have an effect
                                    this%f_nacq_nh4, &
                                    "f_nacq_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    f_nacq_nh4)

    IF ((ABS(this%transform_cost_no3 - missval) > eps8) .OR. & 
        (ABS(this%transform_cost_nh4 - missval) > eps8)) THEN
        dependency_recalculation = .TRUE.
        f_nacq_no3               = transform_cost_nh4 / (transform_cost_no3 + transform_cost_nh4)   ! copied from mo_veg_constants
    ENDIF    
    CALL Set_quincy_parameter_value(transform_cost_nh4 / (transform_cost_no3 + transform_cost_nh4), &   ! only proportional changes may have an effect
                                    this%f_nacq_no3, &
                                    "f_nacq_no3", proportional_parameter_changes, dependency_recalculation, &
                                    f_nacq_no3)
    
    CALL Set_quincy_parameter_value(this%vmax_bnf * 0.0225_wp / (molar_mass_N/molar_mass_C) / (one_day*one_year) * 1.e6_wp, &
                                    this%vmax_bnf, &
                                    "vmax_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_bnf)

    CALL Set_quincy_parameter_value(this%k_bnf / (molar_mass_C/molar_Mass_N), &
                                    this%k_bnf, &
                                    "k_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    k_bnf)

    IF ((ABS(this%transform_cost_no3 - missval) > eps8) .OR. & 
        (ABS(this%transform_cost_nh4 - missval) > eps8) .OR. & 
        (ABS(this%f_nacq_no3         - missval) > eps8) .OR. & 
        (ABS(this%f_nacq_nh4         - missval) > eps8)) THEN
        dependency_recalculation = .TRUE.
        offset_bnf               = 6.0_wp * (f_nacq_nh4 * transform_cost_nh4 + f_nacq_no3 * transform_cost_no3)   ! copied from mo_veg_constants
    ENDIF    
    CALL Set_quincy_parameter_value(this%offset_bnf * (f_nacq_nh4 * transform_cost_nh4 + f_nacq_no3 * transform_cost_no3), &
                                    this%offset_bnf, &
                                    "offset_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    offset_bnf)
    
    CALL Set_quincy_parameter_value(this%eta_nfixation, &
                                    this%eta_nfixation, &
                                    "eta_nfixation", proportional_parameter_changes, dependency_recalculation, &
                                    eta_nfixation)

    CALL Set_quincy_parameter_value(this%k0_fn_chl_C3, &
                                    this%k0_fn_chl_C3, &
                                    "k0_fn_chl_C3", proportional_parameter_changes, dependency_recalculation, &
                                    k0_fn_chl_C3)

    CALL Set_quincy_parameter_value(this%k0_fn_chl_C4, &
                                    this%k0_fn_chl_C4, &
                                    "k0_fn_chl_C4", proportional_parameter_changes, dependency_recalculation, &
                                    k0_fn_chl_C4)
    
    CALL Set_quincy_parameter_value(this%k1_fn_chl_C3, &
                                    this%k1_fn_chl_C3, &
                                    "k1_fn_chl_C3", proportional_parameter_changes, dependency_recalculation, &
                                    k1_fn_chl_C3)
    
    CALL Set_quincy_parameter_value(this%k1_fn_chl_C4, &
                                    this%k1_fn_chl_C4, &
                                    "k1_fn_chl_C4", proportional_parameter_changes, dependency_recalculation, &
                                    k1_fn_chl_C4)

    CALL Set_quincy_parameter_value(this%kfn_chl, &
                                    this%kfn_chl, &
                                    "kfn_chl", proportional_parameter_changes, dependency_recalculation, &
                                    kfn_chl)
    
    CALL Set_quincy_parameter_value(this%min_lai, &
                                    this%min_lai, &
                                    "min_lai", proportional_parameter_changes, dependency_recalculation, &
                                    min_lai)
    
    CALL Set_quincy_parameter_value(this%target_lai_max, &
                                    this%target_lai_max, &
                                    "target_lai_max", proportional_parameter_changes, dependency_recalculation, &
                                    target_lai_max)

    CALL Set_quincy_parameter_value(this%max_leaf_shedding_rate, &
                                    this%max_leaf_shedding_rate, &
                                    "max_leaf_shedding_rate", proportional_parameter_changes, dependency_recalculation, &
                                    max_leaf_shedding_rate)
    
    CALL Set_quincy_parameter_value(this%k_leafon_canopy, &
                                    this%k_leafon_canopy, &
                                    "k_leafon_canopy", proportional_parameter_changes, dependency_recalculation, &
                                    k_leafon_canopy)
    
    CALL Set_quincy_parameter_value(this%background_mort_rate_tree, &
                                    this%background_mort_rate_tree, &
                                    "background_mort_rate_tree", proportional_parameter_changes, dependency_recalculation, &
                                    background_mort_rate_tree)

    CALL Set_quincy_parameter_value(this%background_mort_rate_grass, &
                                    this%background_mort_rate_grass, &
                                    "background_mort_rate_grass", proportional_parameter_changes, dependency_recalculation, &
                                    background_mort_rate_grass)

    CALL Set_quincy_parameter_value(this%lambda_mort_light, &
                                    this%lambda_mort_light, &
                                    "lambda_mort_light", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_mort_light)
    
    CALL Set_quincy_parameter_value(this%k_mort_light, &
                                    this%k_mort_light, &
                                    "k_mort_light", proportional_parameter_changes, dependency_recalculation, &
                                    k_mort_light)
    
    CALL Set_quincy_parameter_value(this%lambda_est_temp, &
                                    this%lambda_est_temp, &
                                    "lambda_est_temp", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_est_temp)

    CALL Set_quincy_parameter_value(this%k_est_temp, &
                                    this%k_est_temp, &
                                    "k_est_temp", proportional_parameter_changes, dependency_recalculation, &
                                    k_est_temp)
    
    CALL Set_quincy_parameter_value(this%lambda_est_moist, &
                                    this%lambda_est_moist, &
                                    "lambda_est_moist", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_est_moist)
    
    CALL Set_quincy_parameter_value(this%k_est_moist, &
                                    this%k_est_moist, &
                                    "k_est_moist", proportional_parameter_changes, dependency_recalculation, &
                                    k_est_moist)

    CALL Set_quincy_parameter_value(this%wood_extraction_rate, &
                                    this%wood_extraction_rate, &
                                    "wood_extraction_rate", proportional_parameter_changes, dependency_recalculation, &
                                    wood_extraction_rate)
    
    CALL Set_quincy_parameter_value(this%k1_mort_greff, &
                                    this%k1_mort_greff, &
                                    "k1_mort_greff", proportional_parameter_changes, dependency_recalculation, &
                                    k1_mort_greff)
    
    CALL Set_quincy_parameter_value(this%k2_mort_greff, &
                                    this%k2_mort_greff, &
                                    "k2_mort_greff", proportional_parameter_changes, dependency_recalculation, &
                                    k2_mort_greff)

    CALL Set_quincy_parameter_value(this%k3_mort_greff, &
                                    this%k3_mort_greff, &
                                    "k3_mort_greff", proportional_parameter_changes, dependency_recalculation, &
                                    k3_mort_greff)

    CALL Set_quincy_parameter_value(this%min_greff, &
                                    this%min_greff, &
                                    "min_greff", proportional_parameter_changes, dependency_recalculation, &
                                    min_greff)
    
    CALL Set_quincy_parameter_value(this%max_crown_area, &
                                    this%max_crown_area, &
                                    "max_crown_area", proportional_parameter_changes, dependency_recalculation, &
                                    max_crown_area)
    
    CALL Set_quincy_parameter_value(this%min_diameter, &
                                    this%min_diameter, &
                                    "min_diameter", proportional_parameter_changes, dependency_recalculation, &
                                    min_diameter)

    CALL Set_quincy_parameter_value(this%k_crown_area, &
                                    this%k_crown_area, &
                                    "k_crown_area", proportional_parameter_changes, dependency_recalculation, &
                                    k_crown_area)
    
    CALL Set_quincy_parameter_value(this%k_rp, &
                                    this%k_rp, &
                                    "k_rp", proportional_parameter_changes, dependency_recalculation, &
                                    k_rp)
    
    CALL Set_quincy_parameter_value(this%k_fpc, &
                                    this%k_fpc, &
                                    "k_fpc", proportional_parameter_changes, dependency_recalculation, &
                                    k_fpc)

    CALL Set_quincy_parameter_value(this%fpc_max, &
                                    this%fpc_max, &
                                    "fpc_max", proportional_parameter_changes, dependency_recalculation, &
                                    fpc_max)
    
    CALL Set_quincy_parameter_value(this%carbon_per_dryweight_leaf, &
                                    this%carbon_per_dryweight_leaf, &
                                    "carbon_per_dryweight_leaf", proportional_parameter_changes, dependency_recalculation, &
                                    carbon_per_dryweight_leaf)
    
    CALL Set_quincy_parameter_value(this%min_height, &
                                    this%min_height, &
                                    "min_height", proportional_parameter_changes, dependency_recalculation, &
                                    min_height)

    CALL Set_quincy_parameter_value(this%sm2lm_grass, &
                                    this%sm2lm_grass, &
                                    "sm2lm_grass", proportional_parameter_changes, dependency_recalculation, &
                                    sm2lm_grass)

    CALL Set_quincy_parameter_value(this%lambda_fruit_alloc, &
                                    this%lambda_fruit_alloc, &
                                    "lambda_fruit_alloc", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_fruit_alloc)
    
    CALL Set_quincy_parameter_value(this%k1_fruit_alloc, &
                                    this%k1_fruit_alloc, &
                                    "k1_fruit_alloc", proportional_parameter_changes, dependency_recalculation, &
                                    k1_fruit_alloc)
    
    CALL Set_quincy_parameter_value(this%k3_fruit_alloc, &
                                    this%k3_fruit_alloc, &
                                    "k3_fruit_alloc", proportional_parameter_changes, dependency_recalculation, &
                                    k3_fruit_alloc)

    CALL Set_quincy_parameter_value(this%k4_fruit_alloc, &
                                    this%k4_fruit_alloc, &
                                    "k4_fruit_alloc", proportional_parameter_changes, dependency_recalculation, &
                                    k4_fruit_alloc)
    
    CALL Set_quincy_parameter_value(this%w_root_lim_max, &
                                    this%w_root_lim_max, &
                                    "w_root_lim_max", proportional_parameter_changes, dependency_recalculation, &
                                    w_root_lim_max)
    
    CALL Set_quincy_parameter_value(this%leaf2root_min_ratio, &
                                    this%leaf2root_min_ratio, &
                                    "leaf2root_min_ratio", proportional_parameter_changes, dependency_recalculation, &
                                    leaf2root_min_ratio)

    CALL Set_quincy_parameter_value(this%leaf2root_max_ratio, &
                                    this%leaf2root_max_ratio, &
                                    "leaf2root_max_ratio", proportional_parameter_changes, dependency_recalculation, &
                                    leaf2root_max_ratio)
    
    CALL Set_quincy_parameter_value(this%delta_n_fraction, &
                                    this%delta_n_fraction, &
                                    "delta_n_fraction", proportional_parameter_changes, dependency_recalculation, &
                                    delta_n_fraction)
    
    CALL Set_quincy_parameter_value(this%delta_n_leaf, &
                                    this%delta_n_leaf, &
                                    "delta_n_leaf", proportional_parameter_changes, dependency_recalculation, &
                                    delta_n_leaf)

    ! SB_CTL
    CALL Set_quincy_parameter_value(this%tau_soluable_litter * one_day, &
                                    this%tau_soluable_litter, &
                                    "tau_soluable_litter", proportional_parameter_changes, dependency_recalculation, &
                                    tau_soluable_litter)
    
    CALL Set_quincy_parameter_value(this%tau_metabolic_litter * one_year * one_day, &
                                    this%tau_metabolic_litter, &
                                    "tau_metabolic_litter", proportional_parameter_changes, dependency_recalculation, &
                                    tau_metabolic_litter)
    
    CALL Set_quincy_parameter_value(this%tau_structural_litter * one_year * one_day, &
                                    this%tau_structural_litter, &
                                    "tau_structural_litter", proportional_parameter_changes, dependency_recalculation, &
                                    tau_structural_litter)

    CALL Set_quincy_parameter_value(this%tau_woody_litter * one_year * one_day, &
                                    this%tau_woody_litter, &
                                    "tau_woody_litter", proportional_parameter_changes, dependency_recalculation, &
                                    tau_woody_litter)
    
    CALL Set_quincy_parameter_value(this%tau_fast * one_year * one_day, &
                                    this%tau_fast, &
                                    "tau_fast", proportional_parameter_changes, dependency_recalculation, &
                                    tau_fast)
    
    CALL Set_quincy_parameter_value(this%tau_slow * one_year * one_day, &
                                    this%tau_slow, &
                                    "tau_slow", proportional_parameter_changes, dependency_recalculation, &
                                    tau_slow)

    CALL Set_quincy_parameter_value(this%tau_biomineralisation * one_year * one_day, &
                                    this%tau_biomineralisation, &
                                    "tau_biomineralisation", proportional_parameter_changes, dependency_recalculation, &
                                    tau_biomineralisation)
    
    CALL Set_quincy_parameter_value(this%t_ref_tresponse, &
                                    this%t_ref_tresponse, &
                                    "t_ref_tresponse", proportional_parameter_changes, dependency_recalculation, &
                                    t_ref_tresponse)

    CALL Set_quincy_parameter_value(this%ea_depolymerisation, &
                                    this%ea_depolymerisation, &
                                    "ea_depolymerisation", proportional_parameter_changes, dependency_recalculation, &
                                    ea_depolymerisation)
    
    CALL Set_quincy_parameter_value(this%ea_mic_uptake, &
                                    this%ea_mic_uptake, &
                                    "ea_mic_uptake", proportional_parameter_changes, dependency_recalculation, &
                                    ea_mic_uptake)
    
    CALL Set_quincy_parameter_value(this%ea_sorption, &
                                    this%ea_sorption, &
                                    "ea_sorption", proportional_parameter_changes, dependency_recalculation, &
                                    ea_sorption)

    CALL Set_quincy_parameter_value(this%ea_desorption, &
                                    this%ea_desorption, &
                                    "ea_desorption", proportional_parameter_changes, dependency_recalculation, &
                                    ea_desorption)
    
    CALL Set_quincy_parameter_value(this%ea_hsc, &
                                    this%ea_hsc, &
                                    "ea_hsc", proportional_parameter_changes, dependency_recalculation, &
                                    ea_hsc)

    CALL Set_quincy_parameter_value(this%ea_decomposition, &
                                    this%ea_decomposition, &
                                    "ea_decomposition", proportional_parameter_changes, dependency_recalculation, &
                                    ea_decomposition)
    
    CALL Set_quincy_parameter_value(this%ed_decomposition, &
                                    this%ed_decomposition, &
                                    "ed_decomposition", proportional_parameter_changes, dependency_recalculation, &
                                    ed_decomposition)
    
    CALL Set_quincy_parameter_value(this%t_ref_decomposition, &   ! in function calculated for proportional changes as:
                                    this%t_ref_decomposition, &   ! def_val = 273.15 + (def_val - 273.15) * nml_prop_chng
                                    "t_ref_decomposition", proportional_parameter_changes, dependency_recalculation, &
                                    t_ref_decomposition)

    CALL Set_quincy_parameter_value(this%k1_afps / 3._wp, &
                                    this%k1_afps, &
                                    "k1_afps", proportional_parameter_changes, dependency_recalculation, &
                                    k1_afps)
    
    CALL Set_quincy_parameter_value(this%kmr_depolymerisation, &
                                    this%kmr_depolymerisation, &
                                    "kmr_depolymerisation", proportional_parameter_changes, dependency_recalculation, &
                                    kmr_depolymerisation)

    CALL Set_quincy_parameter_value(this%kmr_mic_uptake, &
                                    this%kmr_mic_uptake, &
                                    "kmr_mic_uptake", proportional_parameter_changes, dependency_recalculation, &
                                    kmr_mic_uptake)
    
    CALL Set_quincy_parameter_value(this%kmr_hsc, &
                                    this%kmr_hsc, &
                                    "kmr_hsc", proportional_parameter_changes, dependency_recalculation, &
                                    kmr_hsc)
    
    CALL Set_quincy_parameter_value(this%kmr_sorption, &
                                    this%kmr_sorption, &
                                    "kmr_sorption", proportional_parameter_changes, dependency_recalculation, &
                                    kmr_sorption)

    CALL Set_quincy_parameter_value(this%kmr_decomposition, &
                                    this%kmr_decomposition, &
                                    "kmr_decomposition", proportional_parameter_changes, dependency_recalculation, &
                                    kmr_decomposition)
    
    CALL Set_quincy_parameter_value(this%k_diff_org / one_year / one_day, &
                                    this%k_diff_org, &
                                    "k_diff_org", proportional_parameter_changes, dependency_recalculation, &
                                    k_diff_org)

    CALL Set_quincy_parameter_value(this%rho_bulk_org, &
                                    this%rho_bulk_org, &
                                    "rho_bulk_org", proportional_parameter_changes, dependency_recalculation, &
                                    rho_bulk_org)
    
    CALL Set_quincy_parameter_value(this%carbon_per_dryweight_SOM, &
                                    this%carbon_per_dryweight_SOM, &
                                    "carbon_per_dryweight_SOM", proportional_parameter_changes, dependency_recalculation, &
                                    carbon_per_dryweight_SOM)
    
    CALL Set_quincy_parameter_value(this%qmax_po4_fine_particle / molar_mass_P, &
                                    this%qmax_po4_fine_particle, &
                                    "qmax_po4_fine_particle", proportional_parameter_changes, dependency_recalculation, &
                                    qmax_po4_fine_particle)

    CALL Set_quincy_parameter_value(this%k_adsorpt_dom / one_year / one_day * molar_mass_C / 1000._wp, &
                                    this%k_adsorpt_dom, &
                                    "k_adsorpt_dom", proportional_parameter_changes, dependency_recalculation, &
                                    k_adsorpt_dom)
    
    CALL Set_quincy_parameter_value(this%k_desorpt_dom / one_year / one_day, &
                                    this%k_desorpt_dom, &
                                    "k_desorpt_dom", proportional_parameter_changes, dependency_recalculation, &
                                    k_desorpt_dom)

    CALL Set_quincy_parameter_value(this%k_adsorpt_det / one_year / one_day * molar_mass_C / 100000._wp, &
                                    this%k_adsorpt_det, &
                                    "k_adsorpt_det", proportional_parameter_changes, dependency_recalculation, &
                                    k_adsorpt_det)
    
    CALL Set_quincy_parameter_value(this%k_desorpt_det / one_year / one_day / 100._wp, &
                                    this%k_desorpt_det, &
                                    "k_desorpt_det", proportional_parameter_changes, dependency_recalculation, &
                                    k_desorpt_det)
    
    CALL Set_quincy_parameter_value(this%k_adsorpt_nh4, &
                                    this%k_adsorpt_nh4, &
                                    "k_adsorpt_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    k_adsorpt_nh4)

    CALL Set_quincy_parameter_value(this%k_desorpt_nh4, &
                                    this%k_desorpt_nh4, &
                                    "k_desorpt_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    k_desorpt_nh4)
    
    CALL Set_quincy_parameter_value(this%k_adsorpt_po4 * 8._wp / 9._wp / 30.0_wp / one_day, &
                                    this%k_adsorpt_po4, &
                                    "k_adsorpt_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_adsorpt_po4)

    CALL Set_quincy_parameter_value(this%k_desorpt_po4 / 30.0_wp / one_day, &
                                    this%k_desorpt_po4, &
                                    "k_desorpt_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_desorpt_po4)
    
    CALL Set_quincy_parameter_value(this%k_adsorpt_f_po4 / 3600._wp, &
                                    this%k_adsorpt_f_po4, &
                                    "k_adsorpt_f_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_adsorpt_f_po4)
    
    CALL Set_quincy_parameter_value(this%k_desorp_f_po4 / 3600._wp / 10._wp, &
                                    this%k_desorp_f_po4, &
                                    "k_desorp_f_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_desorp_f_po4)

    CALL Set_quincy_parameter_value(this%k_desorp_s_po4, &
                                    this%k_desorp_s_po4, &
                                    "k_desorp_s_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_desorp_s_po4)
    
    CALL Set_quincy_parameter_value(this%k_adsorpt_s_po4, &
                                    this%k_adsorpt_s_po4, &
                                    "k_adsorpt_s_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_adsorpt_s_po4)

    CALL Set_quincy_parameter_value(this%microbial_cn * molar_mass_N / molar_mass_C, &
                                    this%microbial_cn, &
                                    "microbial_cn", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_cn)
    
    CALL Set_quincy_parameter_value(this%microbial_np * molar_mass_P / molar_mass_N, &
                                    this%microbial_np, &
                                    "microbial_np", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_np)
    
    CALL Set_quincy_parameter_value(this%SOM_np * molar_mass_P / molar_mass_N, &
                                    this%SOM_np, &
                                    "SOM_np", proportional_parameter_changes, dependency_recalculation, &
                                    SOM_np)
    
    CALL Set_quincy_parameter_value(this%SOM_cn * molar_mass_N / molar_mass_C, &
                                    this%SOM_cn, &
                                    "SOM_cn", proportional_parameter_changes, dependency_recalculation, &
                                    SOM_cn)
    
    CALL Set_quincy_parameter_value(this%microbial_nue, &
                                    this%microbial_nue, &
                                    "microbial_nue", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_nue)

    CALL Set_quincy_parameter_value(this%microbial_pue, &
                                    this%microbial_pue, &
                                    "microbial_pue", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_pue)
    
    CALL Set_quincy_parameter_value(this%microbial_cue, &
                                    this%microbial_cue, &
                                    "microbial_cue", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_cue)

    CALL Set_quincy_parameter_value(this%microbial_cue_max, &
                                    this%microbial_cue_max, &
                                    "microbial_cue_max", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_cue_max)
    
    CALL Set_quincy_parameter_value(this%microbial_cue_min, &
                                    this%microbial_cue_min, &
                                    "microbial_cue_min", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_cue_min)
    
    CALL Set_quincy_parameter_value(this%microbial_cn_max * molar_mass_N / molar_mass_C, &
                                    this%microbial_cn_max, &
                                    "microbial_cn_max", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_cn_max)

    CALL Set_quincy_parameter_value(this%microbial_cn_min * molar_mass_N / molar_mass_C, &
                                    this%microbial_cn_min, &
                                    "microbial_cn_min", proportional_parameter_changes, dependency_recalculation, &
                                    microbial_cn_min)
    
    CALL Set_quincy_parameter_value(this%f_mic_cn, &
                                    this%f_mic_cn, &
                                    "f_mic_cn", proportional_parameter_changes, dependency_recalculation, &
                                    f_mic_cn)

    CALL Set_quincy_parameter_value(this%frac_woody2dec_litter, &
                                    this%frac_woody2dec_litter, &
                                    "frac_woody2dec_litter", proportional_parameter_changes, dependency_recalculation, &
                                    frac_woody2dec_litter)
    
    CALL Set_quincy_parameter_value(this%frac_litter2fast_som, &
                                    this%frac_litter2fast_som, &
                                    "frac_litter2fast_som", proportional_parameter_changes, dependency_recalculation, &
                                    frac_litter2fast_som)
    
    CALL Set_quincy_parameter_value(this%frac_fast2slow_som, &
                                    this%frac_fast2slow_som, &
                                    "frac_fast2slow_som", proportional_parameter_changes, dependency_recalculation, &
                                    frac_fast2slow_som)

    CALL Set_quincy_parameter_value(this%frac_slow2fast_som, &
                                    this%frac_slow2fast_som, &
                                    "frac_slow2fast_som", proportional_parameter_changes, dependency_recalculation, &
                                    frac_slow2fast_som)
    
    CALL Set_quincy_parameter_value(this%fc_woody2polymeric, &
                                    this%fc_woody2polymeric, &
                                    "fc_woody2polymeric", proportional_parameter_changes, dependency_recalculation, &
                                    fc_woody2polymeric)

    CALL Set_quincy_parameter_value(this%fc_soluable2dom, &
                                    this%fc_soluable2dom, &
                                    "fc_soluable2dom", proportional_parameter_changes, dependency_recalculation, &
                                    fc_soluable2dom)
    
    CALL Set_quincy_parameter_value(this%k_mort_microbial * 1.e-8_wp, &
                                    this%k_mort_microbial, &
                                    "k_mort_microbial", proportional_parameter_changes, dependency_recalculation, &
                                    k_mort_microbial)
    
    CALL Set_quincy_parameter_value(this%min_mic_biomass, &
                                    this%min_mic_biomass, &
                                    "min_mic_biomass", proportional_parameter_changes, dependency_recalculation, &
                                    min_mic_biomass)

    CALL Set_quincy_parameter_value(this%vmax_mic_uptake / one_year / one_day * 1000000._wp * 10._wp, &
                                    this%vmax_mic_uptake, &
                                    "vmax_mic_uptake", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_mic_uptake)
    
    CALL Set_quincy_parameter_value(this%km_mic_uptake / molar_mass_C * 1000_wp, &
                                    this%km_mic_uptake, &
                                    "km_mic_uptake", proportional_parameter_changes, dependency_recalculation, &
                                    km_mic_uptake)

    CALL Set_quincy_parameter_value(this%km_depolymerisation / 2._wp / molar_mass_C * 1000_wp, &
                                    this%km_depolymerisation, &
                                    "km_depolymerisation", proportional_parameter_changes, dependency_recalculation, &
                                    km_depolymerisation)
    
    CALL Set_quincy_parameter_value(this%vmax_depolymerisation_residue / one_year / one_day, &
                                    this%vmax_depolymerisation_residue, &
                                    "vmax_depolymerisation_residue", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_depolymerisation_residue)
    
    CALL Set_quincy_parameter_value(this%vmax_depolymerisation_litter / one_year / one_day, &
                                    this%vmax_depolymerisation_litter, &
                                    "vmax_depolymerisation_litter", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_depolymerisation_litter)

    CALL Set_quincy_parameter_value(this%f_soluable_residue, &
                                    this%f_soluable_residue, &
                                    "f_soluable_residue", proportional_parameter_changes, dependency_recalculation, &
                                    f_soluable_residue)
    
    CALL Set_quincy_parameter_value(this%f_soluable_leaf, &
                                    this%f_soluable_leaf, &
                                    "f_soluable_leaf", proportional_parameter_changes, dependency_recalculation, &
                                    f_soluable_leaf)

    CALL Set_quincy_parameter_value(this%f_soluable_fine_root, &
                                    this%f_soluable_fine_root, &
                                    "f_soluable_fine_root", proportional_parameter_changes, dependency_recalculation, &
                                    f_soluable_fine_root)
    
    CALL Set_quincy_parameter_value(this%f_soluable_fruit, &
                                    this%f_soluable_fruit, &
                                    "f_soluable_fruit", proportional_parameter_changes, dependency_recalculation, &
                                    f_soluable_fruit)
    
    CALL Set_quincy_parameter_value(this%f_soluable_seed_bed, &
                                    this%f_soluable_seed_bed, &
                                    "f_soluable_seed_bed", proportional_parameter_changes, dependency_recalculation, &
                                    f_soluable_seed_bed)

    CALL Set_quincy_parameter_value(this%fc_soluable_max, &
                                    this%fc_soluable_max, &
                                    "fc_soluable_max", proportional_parameter_changes, dependency_recalculation, &
                                    fc_soluable_max)
    
    CALL Set_quincy_parameter_value(this%k_fc_soluable, &
                                    this%k_fc_soluable, &
                                    "k_fc_soluable", proportional_parameter_changes, dependency_recalculation, &
                                    k_fc_soluable)

    CALL Set_quincy_parameter_value(this%k_fn_soluable, &
                                    this%k_fn_soluable, &
                                    "k_fn_soluable", proportional_parameter_changes, dependency_recalculation, &
                                    k_fn_soluable)
    
    CALL Set_quincy_parameter_value(this%k_fp_soluable, &
                                    this%k_fp_soluable, &
                                    "k_fp_soluable", proportional_parameter_changes, dependency_recalculation, &
                                    k_fp_soluable)
    
    CALL Set_quincy_parameter_value(this%lc_leaf_max / molar_mass_C * molar_mass_N, &
                                    this%lc_leaf_max, &
                                    "lc_leaf_max", proportional_parameter_changes, dependency_recalculation, &
                                    lc_leaf_max)

    CALL Set_quincy_parameter_value(this%lc_leaf2sla / molar_mass_C * molar_mass_N, &
                                    this%lc_leaf2sla, &
                                    "lc_leaf2sla", proportional_parameter_changes, dependency_recalculation, &
                                    lc_leaf2sla)
    
    CALL Set_quincy_parameter_value(this%lc_fine_root / molar_mass_C * molar_mass_N, &
                                    this%lc_fine_root, &
                                    "lc_fine_root", proportional_parameter_changes, dependency_recalculation, &
                                    lc_fine_root)

    CALL Set_quincy_parameter_value(this%lc_coarse_root / molar_mass_C * molar_mass_N, &
                                    this%lc_coarse_root, &
                                    "lc_coarse_root", proportional_parameter_changes, dependency_recalculation, &
                                    lc_coarse_root)
    
    CALL Set_quincy_parameter_value(this%lc_woody_litter / molar_mass_C * molar_mass_N, &
                                    this%lc_woody_litter, &
                                    "lc_woody_litter", proportional_parameter_changes, dependency_recalculation, &
                                    lc_woody_litter)
    
    CALL Set_quincy_parameter_value(this%lc_fruit / molar_mass_C * molar_mass_N, &
                                    this%lc_fruit, &
                                    "lc_fruit", proportional_parameter_changes, dependency_recalculation, &
                                    lc_fruit)

    CALL Set_quincy_parameter_value(this%lc_seed_bed / molar_mass_C * molar_mass_N, &
                                    this%lc_seed_bed, &
                                    "lc_seed_bed", proportional_parameter_changes, dependency_recalculation, &
                                    lc_seed_bed)
    
    CALL Set_quincy_parameter_value(this%mycorrhiza_cue, &
                                    this%mycorrhiza_cue, &
                                    "mycorrhiza_cue", proportional_parameter_changes, dependency_recalculation, &
                                    mycorrhiza_cue)

    CALL Set_quincy_parameter_value(this%mycorrhiza_cn * molar_mass_N / molar_mass_C, &
                                    this%mycorrhiza_cn, &
                                    "mycorrhiza_cn", proportional_parameter_changes, dependency_recalculation, &
                                    mycorrhiza_cn)
    
    CALL Set_quincy_parameter_value(this%mycorrhiza_np * molar_mass_P / molar_mass_N, &
                                    this%mycorrhiza_np, &
                                    "mycorrhiza_np", proportional_parameter_changes, dependency_recalculation, &
                                    mycorrhiza_np)
    
    CALL Set_quincy_parameter_value(this%f_maxresp_myc, &
                                    this%f_maxresp_myc, &
                                    "f_maxresp_myc", proportional_parameter_changes, dependency_recalculation, &
                                    f_maxresp_myc)

    CALL Set_quincy_parameter_value(this%f_surface_myc2fineroot, &
                                    this%f_surface_myc2fineroot, &
                                    "f_surface_myc2fineroot", proportional_parameter_changes, dependency_recalculation, &
                                    f_surface_myc2fineroot)
    
    CALL Set_quincy_parameter_value(this%km1_up_nh4, &
                                    this%km1_up_nh4, &
                                    "km1_up_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    km1_up_nh4)

    CALL Set_quincy_parameter_value(this%km2_up_nh4, &
                                    this%km2_up_nh4, &
                                    "km2_up_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    km2_up_nh4)
    
    CALL Set_quincy_parameter_value(this%km1_up_no3, &
                                    this%km1_up_no3, &
                                    "km1_up_no3", proportional_parameter_changes, dependency_recalculation, &
                                    km1_up_no3)
    
    CALL Set_quincy_parameter_value(this%km2_up_no3, &
                                    this%km2_up_no3, &
                                    "km2_up_no3", proportional_parameter_changes, dependency_recalculation, &
                                    km2_up_no3)

    CALL Set_quincy_parameter_value(this%km1_up_po4 / 0.0003_wp, &
                                    this%km1_up_po4, &
                                    "km1_up_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km1_up_po4)
    
    CALL Set_quincy_parameter_value(this%km2_up_po4 / 1000000._wp, &
                                    this%km2_up_po4, &
                                    "km2_up_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km2_up_po4)

    CALL Set_quincy_parameter_value(this%vmax_mic_uptake_nh4 / 3600._wp * molar_mass_C, &
                                    this%vmax_mic_uptake_nh4, &
                                    "vmax_mic_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_mic_uptake_nh4)
    
    CALL Set_quincy_parameter_value(this%km_mic_uptake_nh4 / molar_mass_N, &
                                    this%km_mic_uptake_nh4, &
                                    "km_mic_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    km_mic_uptake_nh4)
    
    CALL Set_quincy_parameter_value(this%vmax_mic_uptake_no3 / 3600._wp * molar_mass_C, &
                                    this%vmax_mic_uptake_no3, &
                                    "vmax_mic_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_mic_uptake_no3)

    CALL Set_quincy_parameter_value(this%km_mic_uptake_no3 / molar_mass_N, &
                                    this%km_mic_uptake_no3, &
                                    "km_mic_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    km_mic_uptake_no3)
    
    CALL Set_quincy_parameter_value(this%f_mic_carrier_n / 1000._wp, &
                                    this%f_mic_carrier_n, &
                                    "f_mic_carrier_n", proportional_parameter_changes, dependency_recalculation, &
                                    f_mic_carrier_n)

    CALL Set_quincy_parameter_value(this%vmax_plant_uptake_nh4 / 3600._wp * molar_mass_C, &
                                    this%vmax_plant_uptake_nh4, &
                                    "vmax_plant_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_plant_uptake_nh4)
    
    CALL Set_quincy_parameter_value(this%km_plant_uptake_nh4 / molar_mass_N, &
                                    this%km_plant_uptake_nh4, &
                                    "km_plant_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    km_plant_uptake_nh4)
    
    CALL Set_quincy_parameter_value(this%vmax_plant_uptake_no3 / 3600._wp * molar_mass_C, &
                                    this%vmax_plant_uptake_no3, &
                                    "vmax_plant_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_plant_uptake_no3)

    CALL Set_quincy_parameter_value(this%km_plant_uptake_no3 / molar_mass_N, &
                                    this%km_plant_uptake_no3, &
                                    "km_plant_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    km_plant_uptake_no3)
    
    CALL Set_quincy_parameter_value(this%f_plant_carrier_n / 10000._wp, &
                                    this%f_plant_carrier_n, &
                                    "f_plant_carrier_n", proportional_parameter_changes, dependency_recalculation, &
                                    f_plant_carrier_n)

    CALL Set_quincy_parameter_value(this%vmax_mic_uptake_po4 / 3600._wp * molar_mass_C, &
                                    this%vmax_mic_uptake_po4, &
                                    "vmax_mic_uptake_po4", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_mic_uptake_po4)
    
    CALL Set_quincy_parameter_value(this%km_mic_uptake_po4 / molar_mass_P, &
                                    this%km_mic_uptake_po4, &
                                    "km_mic_uptake_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_mic_uptake_po4)
    
    CALL Set_quincy_parameter_value(this%vmax_plant_uptake_po4 / 60._wp * molar_mass_C, &
                                    this%vmax_plant_uptake_po4, &
                                    "vmax_plant_uptake_po4", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_plant_uptake_po4)

    CALL Set_quincy_parameter_value(this%km_plant_uptake_po4 / molar_mass_P, &
                                    this%km_plant_uptake_po4, &
                                    "km_plant_uptake_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_plant_uptake_po4)
    
    CALL Set_quincy_parameter_value(this%k_mic_uptake_nh4 * 24._wp / one_day, &
                                    this%k_mic_uptake_nh4, &
                                    "k_mic_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    k_mic_uptake_nh4)

    CALL Set_quincy_parameter_value(this%k_plant_uptake_nh4 * 24._wp / one_day, &
                                    this%k_plant_uptake_nh4, &
                                    "k_plant_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    k_plant_uptake_nh4)
    
    CALL Set_quincy_parameter_value(this%k_mic_uptake_no3 * 24._wp / one_day, &
                                    this%k_mic_uptake_no3, &
                                    "k_mic_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    k_mic_uptake_no3)
    
    CALL Set_quincy_parameter_value(this%k_plant_uptake_no3 * 24._wp / one_day, &
                                    this%k_plant_uptake_no3, &
                                    "k_plant_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    k_plant_uptake_no3)

    CALL Set_quincy_parameter_value(this%k_mic_uptake_po4 / one_day, &
                                    this%k_mic_uptake_po4, &
                                    "k_mic_uptake_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_mic_uptake_po4)
    
    CALL Set_quincy_parameter_value(this%k_plant_uptake_po4 / one_day, &
                                    this%k_plant_uptake_po4, &
                                    "k_plant_uptake_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_plant_uptake_po4)

    CALL Set_quincy_parameter_value(this%f_mic_carrier_po4 / 1000._wp, &
                                    this%f_mic_carrier_po4, &
                                    "f_mic_carrier_po4", proportional_parameter_changes, dependency_recalculation, &
                                    f_mic_carrier_po4)
    
    CALL Set_quincy_parameter_value(this%f_plant_carrier_po4 / 10000._wp, &
                                    this%f_plant_carrier_po4, &
                                    "f_plant_carrier_po4", proportional_parameter_changes, dependency_recalculation, &
                                    f_plant_carrier_po4)
    
    CALL Set_quincy_parameter_value(this%qmax_po4_mineral / molar_mass_P, &
                                    this%qmax_po4_mineral, &
                                    "qmax_po4_mineral", proportional_parameter_changes, dependency_recalculation, &
                                    qmax_po4_mineral)

    CALL Set_quincy_parameter_value(this%km_adsorpt_mineral_po4 / molar_mass_P, &
                                    this%km_adsorpt_mineral_po4, &
                                    "km_adsorpt_mineral_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_adsorpt_mineral_po4)
    
    CALL Set_quincy_parameter_value(this%km_adsorpt_OM_po4 / molar_mass_P, &
                                    this%km_adsorpt_OM_po4, &
                                    "km_adsorpt_OM_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_adsorpt_OM_po4)

    CALL Set_quincy_parameter_value(this%ea_nitrification, &
                                    this%ea_nitrification, &
                                    "ea_nitrification", proportional_parameter_changes, dependency_recalculation, &
                                    ea_nitrification)
    
    CALL Set_quincy_parameter_value(this%ed_nitrification, &
                                    this%ed_nitrification, &
                                    "ed_nitrification", proportional_parameter_changes, dependency_recalculation, &
                                    ed_nitrification)
    
    CALL Set_quincy_parameter_value(this%t_opt_nitrification, &   ! in function calculated for proportional changes as:
                                    this%t_opt_nitrification, &   ! def_val = 273.15 + (def_val - 273.15) * nml_prop_chng
                                    "t_opt_nitrification", proportional_parameter_changes, dependency_recalculation, &
                                    t_opt_nitrification)

    CALL Set_quincy_parameter_value(this%f_nit_noy, &
                                    this%f_nit_noy, &
                                    "f_nit_noy", proportional_parameter_changes, dependency_recalculation, &
                                    f_nit_noy)
    
    CALL Set_quincy_parameter_value(this%f_nit_n2o, &
                                    this%f_nit_n2o, &
                                    "f_nit_n2o", proportional_parameter_changes, dependency_recalculation, &
                                    f_nit_n2o)

    CALL Set_quincy_parameter_value(this%ea_denitrification, &
                                    this%ea_denitrification, &
                                    "ea_denitrification", proportional_parameter_changes, dependency_recalculation, &
                                    ea_denitrification)
    
    CALL Set_quincy_parameter_value(this%vmax_nitrification / one_day * 1.e6_wp, &
                                    this%vmax_nitrification, &
                                    "vmax_nitrification", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_nitrification)
    
    CALL Set_quincy_parameter_value(this%vmax_denitrification / one_day * 1.e6_wp, &
                                    this%vmax_denitrification, &
                                    "vmax_denitrification", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_denitrification)

    CALL Set_quincy_parameter_value(this%km_denitrification_c, &
                                    this%km_denitrification_c, &
                                    "km_denitrification_c", proportional_parameter_changes, dependency_recalculation, &
                                    km_denitrification_c)
    
    CALL Set_quincy_parameter_value(this%km_denitrification_no3 * 1000.0_wp / molar_mass_N, &
                                    this%km_denitrification_no3, &
                                    "km_denitrification_no3", proportional_parameter_changes, dependency_recalculation, &
                                    km_denitrification_no3)

    CALL Set_quincy_parameter_value(this%f_denit_noy, &
                                    this%f_denit_noy, &
                                    "f_denit_noy", proportional_parameter_changes, dependency_recalculation, &
                                    f_denit_noy)
    
    CALL Set_quincy_parameter_value(this%f_denit_n2o, &
                                    this%f_denit_n2o, &
                                    "f_denit_n2o", proportional_parameter_changes, dependency_recalculation, &
                                    f_denit_n2o)
    
    CALL Set_quincy_parameter_value(this%ea_gasdiffusion, &
                                    this%ea_gasdiffusion, &
                                    "ea_gasdiffusion", proportional_parameter_changes, dependency_recalculation, &
                                    ea_gasdiffusion)

    CALL Set_quincy_parameter_value(this%lambda_anvf, &
                                    this%lambda_anvf, &
                                    "lambda_anvf", proportional_parameter_changes, dependency_recalculation, &
                                    lambda_anvf)
    
    CALL Set_quincy_parameter_value(this%k_anvf, &
                                    this%k_anvf, &
                                    "k_anvf", proportional_parameter_changes, dependency_recalculation, &
                                    k_anvf)

    CALL Set_quincy_parameter_value(this%floss_nmin, &
                                    this%floss_nmin, &
                                    "floss_nmin", proportional_parameter_changes, dependency_recalculation, &
                                    floss_nmin)
    
    CALL Set_quincy_parameter_value(this%fnit_nmin, &
                                    this%fnit_nmin, &
                                    "fnit_nmin", proportional_parameter_changes, dependency_recalculation, &
                                    fnit_nmin)
    
    CALL Set_quincy_parameter_value(this%fleach_nh4, &
                                    this%fleach_nh4, &
                                    "fleach_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    fleach_nh4)

    CALL Set_quincy_parameter_value(this%fleach_po4, &
                                    this%fleach_po4, &
                                    "fleach_po4", proportional_parameter_changes, dependency_recalculation, &
                                    fleach_po4)
    
    CALL Set_quincy_parameter_value(this%soil_n_threshold_bnf, &
                                    this%soil_n_threshold_bnf, &
                                    "soil_n_threshold_bnf", proportional_parameter_changes, dependency_recalculation, &
                                    soil_n_threshold_bnf)

    CALL Set_quincy_parameter_value(this%standard_bnf_max_rate / molar_mass_N / one_year / one_day * 1000000.0_wp, &
                                    this%standard_bnf_max_rate, &
                                    "standard_bnf_max_rate", proportional_parameter_changes, dependency_recalculation, &
                                    standard_bnf_max_rate)
    
    CALL Set_quincy_parameter_value(this%eta_mic_uptake_nh4, &
                                    this%eta_mic_uptake_nh4, &
                                    "eta_mic_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    eta_mic_uptake_nh4)
    
    CALL Set_quincy_parameter_value(this%eta_mic_uptake_no3, &
                                    this%eta_mic_uptake_no3, &
                                    "eta_mic_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    eta_mic_uptake_no3)

    CALL Set_quincy_parameter_value(this%eta_plant_uptake_nh4, &
                                    this%eta_plant_uptake_nh4, &
                                    "eta_plant_uptake_nh4", proportional_parameter_changes, dependency_recalculation, &
                                    eta_plant_uptake_nh4)
    
    CALL Set_quincy_parameter_value(this%eta_plant_uptake_no3, &
                                    this%eta_plant_uptake_no3, &
                                    "eta_plant_uptake_no3", proportional_parameter_changes, dependency_recalculation, &
                                    eta_plant_uptake_no3)

    CALL Set_quincy_parameter_value(this%eta_nitrification, &
                                    this%eta_nitrification, &
                                    "eta_nitrification", proportional_parameter_changes, dependency_recalculation, &
                                    eta_nitrification)
    
    CALL Set_quincy_parameter_value(this%eta_nitrate_production, &
                                    this%eta_nitrate_production, &
                                    "eta_nitrate_production", proportional_parameter_changes, dependency_recalculation, &
                                    eta_nitrate_production)
    
    CALL Set_quincy_parameter_value(this%eta_denitrification, &
                                    this%eta_denitrification, &
                                    "eta_denitrification", proportional_parameter_changes, dependency_recalculation, &
                                    eta_denitrification)

    CALL Set_quincy_parameter_value(this%eta_ammonification, &
                                    this%eta_ammonification, &
                                    "eta_ammonification", proportional_parameter_changes, dependency_recalculation, &
                                    eta_ammonification)
    
    CALL Set_quincy_parameter_value(this%eta_volatilisation_nh3, &
                                    this%eta_volatilisation_nh3, &
                                    "eta_volatilisation_nh3", proportional_parameter_changes, dependency_recalculation, &
                                    eta_volatilisation_nh3)

    CALL Set_quincy_parameter_value(this%k_weath_mineral / molar_mass_P / 1254._wp / one_year /one_day, &
                                    this%k_weath_mineral, &
                                    "k_weath_mineral", proportional_parameter_changes, dependency_recalculation, &
                                    k_weath_mineral)
    
    CALL Set_quincy_parameter_value(this%eq_sorp_mineral_po4 / molar_mass_P, &
                                    this%eq_sorp_mineral_po4, &
                                    "eq_sorp_mineral_po4", proportional_parameter_changes, dependency_recalculation, &
                                    eq_sorp_mineral_po4)
    
    CALL Set_quincy_parameter_value(this%k_occlude_po4 / 1000000._wp / 30._wp / one_day, &
                                    this%k_occlude_po4, &
                                    "k_occlude_po4", proportional_parameter_changes, dependency_recalculation, &
                                    k_occlude_po4)

    CALL Set_quincy_parameter_value(this%vmax_biochem_po4 / 40.0_wp / one_day, &
                                    this%vmax_biochem_po4, &
                                    "vmax_biochem_po4", proportional_parameter_changes, dependency_recalculation, &
                                    vmax_biochem_po4)
    
    CALL Set_quincy_parameter_value(this%km_mic_biochem_po4 / molar_mass_C, &
                                    this%km_mic_biochem_po4, &
                                    "km_mic_biochem_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_mic_biochem_po4)

    CALL Set_quincy_parameter_value(this%km_dip_biochem_po4, &
                                    this%km_dip_biochem_po4, &
                                    "km_dip_biochem_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_dip_biochem_po4)
    
    CALL Set_quincy_parameter_value(this%km_rootc_biochem_po4, &
                                    this%km_rootc_biochem_po4, &
                                    "km_rootc_biochem_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_rootc_biochem_po4)
    
    CALL Set_quincy_parameter_value(this%km_rootc_weath_po4, &
                                    this%km_rootc_weath_po4, &
                                    "km_rootc_weath_po4", proportional_parameter_changes, dependency_recalculation, &
                                    km_rootc_weath_po4)
    
    CALL Set_quincy_parameter_value(this%biomin_thres_PC / 2400._wp, &
                                    this%biomin_thres_PC, &
                                    "biomin_thres_PC", proportional_parameter_changes, dependency_recalculation, &
                                    biomin_thres_PC)
    
    CALL Set_quincy_parameter_value(this%fc_p_res2dom, &
                                    this%fc_p_res2dom, &
                                    "fc_p_res2dom", proportional_parameter_changes, dependency_recalculation, &
                                    fc_p_res2dom)
    
    CALL Set_quincy_parameter_value(this%fc_n_res2dom, &
                                    this%fc_n_res2dom, &
                                    "fc_n_res2dom", proportional_parameter_changes, dependency_recalculation, &
                                    fc_n_res2dom)
    
  END SUBROUTINE set_parameter_slm_values_with_list_values


  !-----------------------------------------------------------------------------------------------------
  !> read lctlib parameter values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_parameter_lctlib_list

    USE mo_namelist,           ONLY: open_nml, POSITIONED, position_nml, close_nml  ! read model options from file

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    !  just working on the above defined (private) objects 'param_pft1_from_file' etc.
    ! ---------------------------
    ! 0.2 Local
    ! PFT 1
    REAL(wp)            :: sla_p1
    REAL(wp)            :: sigma_vis_p1
    REAL(wp)            :: sigma_nir_p1
    REAL(wp)            :: omega_clumping_p1
    REAL(wp)            :: crown_shape_factor_p1
    REAL(wp)            :: cn_leaf_p1
    REAL(wp)            :: cn_leaf_min_p1
    REAL(wp)            :: cn_leaf_max_p1
    REAL(wp)            :: np_leaf_p1
    REAL(wp)            :: np_leaf_min_p1
    REAL(wp)            :: np_leaf_max_p1
    REAL(wp)            :: k0_fn_struc_p1
    REAL(wp)            :: fn_oth_min_p1
    REAL(wp)            :: t_jmax_opt_p1
    REAL(wp)            :: t_jmax_omega_p1
    REAL(wp)            :: g0_p1
    REAL(wp)            :: g1_p1
    REAL(wp)            :: gmin_p1
    REAL(wp)            :: tau_leaf_p1
    REAL(wp)            :: tau_fine_root_p1
    REAL(wp)            :: tau_coarse_root_p1
    REAL(wp)            :: tau_branch_p1
    REAL(wp)            :: tau_sap_wood_p1
    REAL(wp)            :: tau_fruit_p1
    REAL(wp)            :: tau_seed_litter_p1
    REAL(wp)            :: tau_seed_est_p1
    REAL(wp)            :: tau_mycorrhiza_p1
    REAL(wp)            :: vmax_uptake_n_p1
    REAL(wp)            :: vmax_uptake_p_p1
    REAL(wp)            :: bnf_base_p1
    REAL(wp)            :: lambda_est_light_p1
    REAL(wp)            :: k_est_light_p1
    REAL(wp)            :: seed_size_p1
    REAL(wp)            :: beta_soil_flush_p1
    REAL(wp)            :: beta_soil_senescence_p1
    REAL(wp)            :: gdd_req_max_p1
    REAL(wp)            :: k_gdd_dormance_p1
    REAL(wp)            :: t_air_senescence_p1
    REAL(wp)            :: min_leaf_age_p1
    REAL(wp)            :: frac_sapwood_branch_p1
    REAL(wp)            :: wood_density_p1
    REAL(wp)            :: k_latosa_p1
    REAL(wp)            :: k_crtos_p1
    REAL(wp)            :: k_rtos_p1
    REAL(wp)            :: k2_fruit_alloc_p1
    REAL(wp)            :: allom_k1_p1
    REAL(wp)            :: allom_k2_p1
    REAL(wp)            :: phi_leaf_min_p1
    REAL(wp)            :: k_root_p1
    REAL(wp)            :: k_sapwood_p1
    REAL(wp)            :: c0_allom_p1
    REAL(wp)            :: fstore_target_p1
    REAL(wp)            :: k_root_dist_p1
    ! PFT 2
    REAL(wp)            :: sla_p2
    REAL(wp)            :: sigma_vis_p2
    REAL(wp)            :: sigma_nir_p2
    REAL(wp)            :: omega_clumping_p2
    REAL(wp)            :: crown_shape_factor_p2
    REAL(wp)            :: cn_leaf_p2
    REAL(wp)            :: cn_leaf_min_p2
    REAL(wp)            :: cn_leaf_max_p2
    REAL(wp)            :: np_leaf_p2
    REAL(wp)            :: np_leaf_min_p2
    REAL(wp)            :: np_leaf_max_p2
    REAL(wp)            :: k0_fn_struc_p2
    REAL(wp)            :: fn_oth_min_p2
    REAL(wp)            :: t_jmax_opt_p2
    REAL(wp)            :: t_jmax_omega_p2
    REAL(wp)            :: g0_p2
    REAL(wp)            :: g1_p2
    REAL(wp)            :: gmin_p2
    REAL(wp)            :: tau_leaf_p2
    REAL(wp)            :: tau_fine_root_p2
    REAL(wp)            :: tau_coarse_root_p2
    REAL(wp)            :: tau_branch_p2
    REAL(wp)            :: tau_sap_wood_p2
    REAL(wp)            :: tau_fruit_p2
    REAL(wp)            :: tau_seed_litter_p2
    REAL(wp)            :: tau_seed_est_p2
    REAL(wp)            :: tau_mycorrhiza_p2
    REAL(wp)            :: vmax_uptake_n_p2
    REAL(wp)            :: vmax_uptake_p_p2
    REAL(wp)            :: bnf_base_p2
    REAL(wp)            :: lambda_est_light_p2
    REAL(wp)            :: k_est_light_p2
    REAL(wp)            :: seed_size_p2
    REAL(wp)            :: beta_soil_flush_p2
    REAL(wp)            :: beta_soil_senescence_p2
    REAL(wp)            :: gdd_req_max_p2
    REAL(wp)            :: k_gdd_dormance_p2
    REAL(wp)            :: t_air_senescence_p2
    REAL(wp)            :: min_leaf_age_p2
    REAL(wp)            :: frac_sapwood_branch_p2
    REAL(wp)            :: wood_density_p2
    REAL(wp)            :: k_latosa_p2
    REAL(wp)            :: k_crtos_p2
    REAL(wp)            :: k_rtos_p2
    REAL(wp)            :: k2_fruit_alloc_p2
    REAL(wp)            :: allom_k1_p2
    REAL(wp)            :: allom_k2_p2
    REAL(wp)            :: phi_leaf_min_p2
    REAL(wp)            :: k_root_p2
    REAL(wp)            :: k_sapwood_p2
    REAL(wp)            :: c0_allom_p2
    REAL(wp)            :: fstore_target_p2
    REAL(wp)            :: k_root_dist_p2
    ! PFT 3
    REAL(wp)            :: sla_p3
    REAL(wp)            :: sigma_vis_p3
    REAL(wp)            :: sigma_nir_p3
    REAL(wp)            :: omega_clumping_p3
    REAL(wp)            :: crown_shape_factor_p3
    REAL(wp)            :: cn_leaf_p3
    REAL(wp)            :: cn_leaf_min_p3
    REAL(wp)            :: cn_leaf_max_p3
    REAL(wp)            :: np_leaf_p3
    REAL(wp)            :: np_leaf_min_p3
    REAL(wp)            :: np_leaf_max_p3
    REAL(wp)            :: k0_fn_struc_p3
    REAL(wp)            :: fn_oth_min_p3
    REAL(wp)            :: t_jmax_opt_p3
    REAL(wp)            :: t_jmax_omega_p3
    REAL(wp)            :: g0_p3
    REAL(wp)            :: g1_p3
    REAL(wp)            :: gmin_p3
    REAL(wp)            :: tau_leaf_p3
    REAL(wp)            :: tau_fine_root_p3
    REAL(wp)            :: tau_coarse_root_p3
    REAL(wp)            :: tau_branch_p3
    REAL(wp)            :: tau_sap_wood_p3
    REAL(wp)            :: tau_fruit_p3
    REAL(wp)            :: tau_seed_litter_p3
    REAL(wp)            :: tau_seed_est_p3
    REAL(wp)            :: tau_mycorrhiza_p3
    REAL(wp)            :: vmax_uptake_n_p3
    REAL(wp)            :: vmax_uptake_p_p3
    REAL(wp)            :: bnf_base_p3
    REAL(wp)            :: lambda_est_light_p3
    REAL(wp)            :: k_est_light_p3
    REAL(wp)            :: seed_size_p3
    REAL(wp)            :: beta_soil_flush_p3
    REAL(wp)            :: beta_soil_senescence_p3
    REAL(wp)            :: gdd_req_max_p3
    REAL(wp)            :: k_gdd_dormance_p3
    REAL(wp)            :: t_air_senescence_p3
    REAL(wp)            :: min_leaf_age_p3
    REAL(wp)            :: frac_sapwood_branch_p3
    REAL(wp)            :: wood_density_p3
    REAL(wp)            :: k_latosa_p3
    REAL(wp)            :: k_crtos_p3
    REAL(wp)            :: k_rtos_p3
    REAL(wp)            :: k2_fruit_alloc_p3
    REAL(wp)            :: allom_k1_p3
    REAL(wp)            :: allom_k2_p3
    REAL(wp)            :: phi_leaf_min_p3
    REAL(wp)            :: k_root_p3
    REAL(wp)            :: k_sapwood_p3
    REAL(wp)            :: c0_allom_p3
    REAL(wp)            :: fstore_target_p3
    REAL(wp)            :: k_root_dist_p3
    ! PFT 4
    REAL(wp)            :: sla_p4
    REAL(wp)            :: sigma_vis_p4
    REAL(wp)            :: sigma_nir_p4
    REAL(wp)            :: omega_clumping_p4
    REAL(wp)            :: crown_shape_factor_p4
    REAL(wp)            :: cn_leaf_p4
    REAL(wp)            :: cn_leaf_min_p4
    REAL(wp)            :: cn_leaf_max_p4
    REAL(wp)            :: np_leaf_p4
    REAL(wp)            :: np_leaf_min_p4
    REAL(wp)            :: np_leaf_max_p4
    REAL(wp)            :: k0_fn_struc_p4
    REAL(wp)            :: fn_oth_min_p4
    REAL(wp)            :: t_jmax_opt_p4
    REAL(wp)            :: t_jmax_omega_p4
    REAL(wp)            :: g0_p4
    REAL(wp)            :: g1_p4
    REAL(wp)            :: gmin_p4
    REAL(wp)            :: tau_leaf_p4
    REAL(wp)            :: tau_fine_root_p4
    REAL(wp)            :: tau_coarse_root_p4
    REAL(wp)            :: tau_branch_p4
    REAL(wp)            :: tau_sap_wood_p4
    REAL(wp)            :: tau_fruit_p4
    REAL(wp)            :: tau_seed_litter_p4
    REAL(wp)            :: tau_seed_est_p4
    REAL(wp)            :: tau_mycorrhiza_p4
    REAL(wp)            :: vmax_uptake_n_p4
    REAL(wp)            :: vmax_uptake_p_p4
    REAL(wp)            :: bnf_base_p4
    REAL(wp)            :: lambda_est_light_p4
    REAL(wp)            :: k_est_light_p4
    REAL(wp)            :: seed_size_p4
    REAL(wp)            :: beta_soil_flush_p4
    REAL(wp)            :: beta_soil_senescence_p4
    REAL(wp)            :: gdd_req_max_p4
    REAL(wp)            :: k_gdd_dormance_p4
    REAL(wp)            :: t_air_senescence_p4
    REAL(wp)            :: min_leaf_age_p4
    REAL(wp)            :: frac_sapwood_branch_p4
    REAL(wp)            :: wood_density_p4
    REAL(wp)            :: k_latosa_p4
    REAL(wp)            :: k_crtos_p4
    REAL(wp)            :: k_rtos_p4
    REAL(wp)            :: k2_fruit_alloc_p4
    REAL(wp)            :: allom_k1_p4
    REAL(wp)            :: allom_k2_p4
    REAL(wp)            :: phi_leaf_min_p4
    REAL(wp)            :: k_root_p4
    REAL(wp)            :: k_sapwood_p4
    REAL(wp)            :: c0_allom_p4
    REAL(wp)            :: fstore_target_p4
    REAL(wp)            :: k_root_dist_p4
    ! PFT 5
    REAL(wp)            :: sla_p5
    REAL(wp)            :: sigma_vis_p5
    REAL(wp)            :: sigma_nir_p5
    REAL(wp)            :: omega_clumping_p5
    REAL(wp)            :: crown_shape_factor_p5
    REAL(wp)            :: cn_leaf_p5
    REAL(wp)            :: cn_leaf_min_p5
    REAL(wp)            :: cn_leaf_max_p5
    REAL(wp)            :: np_leaf_p5
    REAL(wp)            :: np_leaf_min_p5
    REAL(wp)            :: np_leaf_max_p5
    REAL(wp)            :: k0_fn_struc_p5
    REAL(wp)            :: fn_oth_min_p5
    REAL(wp)            :: t_jmax_opt_p5
    REAL(wp)            :: t_jmax_omega_p5
    REAL(wp)            :: g0_p5
    REAL(wp)            :: g1_p5
    REAL(wp)            :: gmin_p5
    REAL(wp)            :: tau_leaf_p5
    REAL(wp)            :: tau_fine_root_p5
    REAL(wp)            :: tau_coarse_root_p5
    REAL(wp)            :: tau_branch_p5
    REAL(wp)            :: tau_sap_wood_p5
    REAL(wp)            :: tau_fruit_p5
    REAL(wp)            :: tau_seed_litter_p5
    REAL(wp)            :: tau_seed_est_p5
    REAL(wp)            :: tau_mycorrhiza_p5
    REAL(wp)            :: vmax_uptake_n_p5
    REAL(wp)            :: vmax_uptake_p_p5
    REAL(wp)            :: bnf_base_p5
    REAL(wp)            :: lambda_est_light_p5
    REAL(wp)            :: k_est_light_p5
    REAL(wp)            :: seed_size_p5
    REAL(wp)            :: beta_soil_flush_p5
    REAL(wp)            :: beta_soil_senescence_p5
    REAL(wp)            :: gdd_req_max_p5
    REAL(wp)            :: k_gdd_dormance_p5
    REAL(wp)            :: t_air_senescence_p5
    REAL(wp)            :: min_leaf_age_p5
    REAL(wp)            :: frac_sapwood_branch_p5
    REAL(wp)            :: wood_density_p5
    REAL(wp)            :: k_latosa_p5
    REAL(wp)            :: k_crtos_p5
    REAL(wp)            :: k_rtos_p5
    REAL(wp)            :: k2_fruit_alloc_p5
    REAL(wp)            :: allom_k1_p5
    REAL(wp)            :: allom_k2_p5
    REAL(wp)            :: phi_leaf_min_p5
    REAL(wp)            :: k_root_p5
    REAL(wp)            :: k_sapwood_p5
    REAL(wp)            :: c0_allom_p5
    REAL(wp)            :: fstore_target_p5
    REAL(wp)            :: k_root_dist_p5
    ! PFT 6
    REAL(wp)            :: sla_p6
    REAL(wp)            :: sigma_vis_p6
    REAL(wp)            :: sigma_nir_p6
    REAL(wp)            :: omega_clumping_p6
    REAL(wp)            :: crown_shape_factor_p6
    REAL(wp)            :: cn_leaf_p6
    REAL(wp)            :: cn_leaf_min_p6
    REAL(wp)            :: cn_leaf_max_p6
    REAL(wp)            :: np_leaf_p6
    REAL(wp)            :: np_leaf_min_p6
    REAL(wp)            :: np_leaf_max_p6
    REAL(wp)            :: k0_fn_struc_p6
    REAL(wp)            :: fn_oth_min_p6
    REAL(wp)            :: t_jmax_opt_p6
    REAL(wp)            :: t_jmax_omega_p6
    REAL(wp)            :: g0_p6
    REAL(wp)            :: g1_p6
    REAL(wp)            :: gmin_p6
    REAL(wp)            :: tau_leaf_p6
    REAL(wp)            :: tau_fine_root_p6
    REAL(wp)            :: tau_coarse_root_p6
    REAL(wp)            :: tau_branch_p6
    REAL(wp)            :: tau_sap_wood_p6
    REAL(wp)            :: tau_fruit_p6
    REAL(wp)            :: tau_seed_litter_p6
    REAL(wp)            :: tau_seed_est_p6
    REAL(wp)            :: tau_mycorrhiza_p6
    REAL(wp)            :: vmax_uptake_n_p6
    REAL(wp)            :: vmax_uptake_p_p6
    REAL(wp)            :: bnf_base_p6
    REAL(wp)            :: lambda_est_light_p6
    REAL(wp)            :: k_est_light_p6
    REAL(wp)            :: seed_size_p6
    REAL(wp)            :: beta_soil_flush_p6
    REAL(wp)            :: beta_soil_senescence_p6
    REAL(wp)            :: gdd_req_max_p6
    REAL(wp)            :: k_gdd_dormance_p6
    REAL(wp)            :: t_air_senescence_p6
    REAL(wp)            :: min_leaf_age_p6
    REAL(wp)            :: frac_sapwood_branch_p6
    REAL(wp)            :: wood_density_p6
    REAL(wp)            :: k_latosa_p6
    REAL(wp)            :: k_crtos_p6
    REAL(wp)            :: k_rtos_p6
    REAL(wp)            :: k2_fruit_alloc_p6
    REAL(wp)            :: allom_k1_p6
    REAL(wp)            :: allom_k2_p6
    REAL(wp)            :: phi_leaf_min_p6
    REAL(wp)            :: k_root_p6
    REAL(wp)            :: k_sapwood_p6
    REAL(wp)            :: c0_allom_p6
    REAL(wp)            :: fstore_target_p6
    REAL(wp)            :: k_root_dist_p6
    ! PFT 7
    REAL(wp)            :: sla_p7
    REAL(wp)            :: sigma_vis_p7
    REAL(wp)            :: sigma_nir_p7
    REAL(wp)            :: omega_clumping_p7
    REAL(wp)            :: crown_shape_factor_p7
    REAL(wp)            :: cn_leaf_p7
    REAL(wp)            :: cn_leaf_min_p7
    REAL(wp)            :: cn_leaf_max_p7
    REAL(wp)            :: np_leaf_p7
    REAL(wp)            :: np_leaf_min_p7
    REAL(wp)            :: np_leaf_max_p7
    REAL(wp)            :: k0_fn_struc_p7
    REAL(wp)            :: fn_oth_min_p7
    REAL(wp)            :: t_jmax_opt_p7
    REAL(wp)            :: t_jmax_omega_p7
    REAL(wp)            :: g0_p7
    REAL(wp)            :: g1_p7
    REAL(wp)            :: gmin_p7
    REAL(wp)            :: tau_leaf_p7
    REAL(wp)            :: tau_fine_root_p7
    REAL(wp)            :: tau_coarse_root_p7
    REAL(wp)            :: tau_branch_p7
    REAL(wp)            :: tau_sap_wood_p7
    REAL(wp)            :: tau_fruit_p7
    REAL(wp)            :: tau_seed_litter_p7
    REAL(wp)            :: tau_seed_est_p7
    REAL(wp)            :: tau_mycorrhiza_p7
    REAL(wp)            :: vmax_uptake_n_p7
    REAL(wp)            :: vmax_uptake_p_p7
    REAL(wp)            :: bnf_base_p7
    REAL(wp)            :: lambda_est_light_p7
    REAL(wp)            :: k_est_light_p7
    REAL(wp)            :: seed_size_p7
    REAL(wp)            :: beta_soil_flush_p7
    REAL(wp)            :: beta_soil_senescence_p7
    REAL(wp)            :: gdd_req_max_p7
    REAL(wp)            :: k_gdd_dormance_p7
    REAL(wp)            :: t_air_senescence_p7
    REAL(wp)            :: min_leaf_age_p7
    REAL(wp)            :: frac_sapwood_branch_p7
    REAL(wp)            :: wood_density_p7
    REAL(wp)            :: k_latosa_p7
    REAL(wp)            :: k_crtos_p7
    REAL(wp)            :: k_rtos_p7
    REAL(wp)            :: k2_fruit_alloc_p7
    REAL(wp)            :: allom_k1_p7
    REAL(wp)            :: allom_k2_p7
    REAL(wp)            :: phi_leaf_min_p7
    REAL(wp)            :: k_root_p7
    REAL(wp)            :: k_sapwood_p7
    REAL(wp)            :: c0_allom_p7
    REAL(wp)            :: fstore_target_p7
    REAL(wp)            :: k_root_dist_p7
    ! PFT 8
    REAL(wp)            :: sla_p8
    REAL(wp)            :: sigma_vis_p8
    REAL(wp)            :: sigma_nir_p8
    REAL(wp)            :: omega_clumping_p8
    REAL(wp)            :: crown_shape_factor_p8
    REAL(wp)            :: cn_leaf_p8
    REAL(wp)            :: cn_leaf_min_p8
    REAL(wp)            :: cn_leaf_max_p8
    REAL(wp)            :: np_leaf_p8
    REAL(wp)            :: np_leaf_min_p8
    REAL(wp)            :: np_leaf_max_p8
    REAL(wp)            :: k0_fn_struc_p8
    REAL(wp)            :: fn_oth_min_p8
    REAL(wp)            :: t_jmax_opt_p8
    REAL(wp)            :: t_jmax_omega_p8
    REAL(wp)            :: g0_p8
    REAL(wp)            :: g1_p8
    REAL(wp)            :: gmin_p8
    REAL(wp)            :: tau_leaf_p8
    REAL(wp)            :: tau_fine_root_p8
    REAL(wp)            :: tau_coarse_root_p8
    REAL(wp)            :: tau_branch_p8
    REAL(wp)            :: tau_sap_wood_p8
    REAL(wp)            :: tau_fruit_p8
    REAL(wp)            :: tau_seed_litter_p8
    REAL(wp)            :: tau_seed_est_p8
    REAL(wp)            :: tau_mycorrhiza_p8
    REAL(wp)            :: vmax_uptake_n_p8
    REAL(wp)            :: vmax_uptake_p_p8
    REAL(wp)            :: bnf_base_p8
    REAL(wp)            :: lambda_est_light_p8
    REAL(wp)            :: k_est_light_p8
    REAL(wp)            :: seed_size_p8
    REAL(wp)            :: beta_soil_flush_p8
    REAL(wp)            :: beta_soil_senescence_p8
    REAL(wp)            :: gdd_req_max_p8
    REAL(wp)            :: k_gdd_dormance_p8
    REAL(wp)            :: t_air_senescence_p8
    REAL(wp)            :: min_leaf_age_p8
    REAL(wp)            :: frac_sapwood_branch_p8
    REAL(wp)            :: wood_density_p8
    REAL(wp)            :: k_latosa_p8
    REAL(wp)            :: k_crtos_p8
    REAL(wp)            :: k_rtos_p8
    REAL(wp)            :: k2_fruit_alloc_p8
    REAL(wp)            :: allom_k1_p8
    REAL(wp)            :: allom_k2_p8
    REAL(wp)            :: phi_leaf_min_p8
    REAL(wp)            :: k_root_p8
    REAL(wp)            :: k_sapwood_p8
    REAL(wp)            :: c0_allom_p8
    REAL(wp)            :: fstore_target_p8
    REAL(wp)            :: k_root_dist_p8


    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_parameter_lctlib_list'

    ! Namelist definitions for each PFT
    ! 
    NAMELIST /LCTLIB_PFT1_CTL/  &
          sla_p1, &
          sigma_vis_p1, &
          sigma_nir_p1, &
          omega_clumping_p1, &
          crown_shape_factor_p1, &
          cn_leaf_p1, &
          cn_leaf_min_p1, &
          cn_leaf_max_p1, &
          np_leaf_p1, &
          np_leaf_min_p1, &
          np_leaf_max_p1, &
          k0_fn_struc_p1, &
          fn_oth_min_p1, &
          t_jmax_opt_p1, &
          t_jmax_omega_p1, &
          g0_p1, &
          g1_p1, &
          gmin_p1, &
          tau_leaf_p1, &
          tau_fine_root_p1, &
          tau_coarse_root_p1, &
          tau_branch_p1, &
          tau_sap_wood_p1, &
          tau_fruit_p1, &
          tau_seed_litter_p1, &
          tau_seed_est_p1, &
          tau_mycorrhiza_p1, &
          vmax_uptake_n_p1, &
          vmax_uptake_p_p1, &
          bnf_base_p1, &
          lambda_est_light_p1, &
          k_est_light_p1, &
          seed_size_p1, &
          beta_soil_flush_p1, &
          beta_soil_senescence_p1, &
          gdd_req_max_p1, &
          k_gdd_dormance_p1, &
          t_air_senescence_p1, &
          min_leaf_age_p1, &
          frac_sapwood_branch_p1, &
          wood_density_p1, &
          k_latosa_p1, &
          k_crtos_p1, &
          k_rtos_p1, &
          k2_fruit_alloc_p1, &
          allom_k1_p1, &
          allom_k2_p1, &
          phi_leaf_min_p1, &
          k_root_p1, &
          k_sapwood_p1, &
          c0_allom_p1, &
          fstore_target_p1, &
          k_root_dist_p1

    NAMELIST /LCTLIB_PFT2_CTL/  &
          sla_p2, &
          sigma_vis_p2, &
          sigma_nir_p2, &
          omega_clumping_p2, &
          crown_shape_factor_p2, &
          cn_leaf_p2, &
          cn_leaf_min_p2, &
          cn_leaf_max_p2, &
          np_leaf_p2, &
          np_leaf_min_p2, &
          np_leaf_max_p2, &
          k0_fn_struc_p2, &
          fn_oth_min_p2, &
          t_jmax_opt_p2, &
          t_jmax_omega_p2, &
          g0_p2, &
          g1_p2, &
          gmin_p2, &
          tau_leaf_p2, &
          tau_fine_root_p2, &
          tau_coarse_root_p2, &
          tau_branch_p2, &
          tau_sap_wood_p2, &
          tau_fruit_p2, &
          tau_seed_litter_p2, &
          tau_seed_est_p2, &
          tau_mycorrhiza_p2, &
          vmax_uptake_n_p2, &
          vmax_uptake_p_p2, &
          bnf_base_p2, &
          lambda_est_light_p2, &
          k_est_light_p2, &
          seed_size_p2, &
          beta_soil_flush_p2, &
          beta_soil_senescence_p2, &
          gdd_req_max_p2, &
          k_gdd_dormance_p2, &
          t_air_senescence_p2, &
          min_leaf_age_p2, &
          frac_sapwood_branch_p2, &
          wood_density_p2, &
          k_latosa_p2, &
          k_crtos_p2, &
          k_rtos_p2, &
          k2_fruit_alloc_p2, &
          allom_k1_p2, &
          allom_k2_p2, &
          phi_leaf_min_p2, &
          k_root_p2, &
          k_sapwood_p2, &
          c0_allom_p2, &
          fstore_target_p2, &
          k_root_dist_p2

    NAMELIST /LCTLIB_PFT3_CTL/  &
          sla_p3, &
          sigma_vis_p3, &
          sigma_nir_p3, &
          omega_clumping_p3, &
          crown_shape_factor_p3, &
          cn_leaf_p3, &
          cn_leaf_min_p3, &
          cn_leaf_max_p3, &
          np_leaf_p3, &
          np_leaf_min_p3, &
          np_leaf_max_p3, &
          k0_fn_struc_p3, &
          fn_oth_min_p3, &
          t_jmax_opt_p3, &
          t_jmax_omega_p3, &
          g0_p3, &
          g1_p3, &
          gmin_p3, &
          tau_leaf_p3, &
          tau_fine_root_p3, &
          tau_coarse_root_p3, &
          tau_branch_p3, &
          tau_sap_wood_p3, &
          tau_fruit_p3, &
          tau_seed_litter_p3, &
          tau_seed_est_p3, &
          tau_mycorrhiza_p3, &
          vmax_uptake_n_p3, &
          vmax_uptake_p_p3, &
          bnf_base_p3, &
          lambda_est_light_p3, &
          k_est_light_p3, &
          seed_size_p3, &
          beta_soil_flush_p3, &
          beta_soil_senescence_p3, &
          gdd_req_max_p3, &
          k_gdd_dormance_p3, &
          t_air_senescence_p3, &
          min_leaf_age_p3, &
          frac_sapwood_branch_p3, &
          wood_density_p3, &
          k_latosa_p3, &
          k_crtos_p3, &
          k_rtos_p3, &
          k2_fruit_alloc_p3, &
          allom_k1_p3, &
          allom_k2_p3, &
          phi_leaf_min_p3, &
          k_root_p3, &
          k_sapwood_p3, &
          c0_allom_p3, &
          fstore_target_p3, &
          k_root_dist_p3

    NAMELIST /LCTLIB_PFT4_CTL/  &
          sla_p4, &
          sigma_vis_p4, &
          sigma_nir_p4, &
          omega_clumping_p4, &
          crown_shape_factor_p4, &
          cn_leaf_p4, &
          cn_leaf_min_p4, &
          cn_leaf_max_p4, &
          np_leaf_p4, &
          np_leaf_min_p4, &
          np_leaf_max_p4, &
          k0_fn_struc_p4, &
          fn_oth_min_p4, &
          t_jmax_opt_p4, &
          t_jmax_omega_p4, &
          g0_p4, &
          g1_p4, &
          gmin_p4, &
          tau_leaf_p4, &
          tau_fine_root_p4, &
          tau_coarse_root_p4, &
          tau_branch_p4, &
          tau_sap_wood_p4, &
          tau_fruit_p4, &
          tau_seed_litter_p4, &
          tau_seed_est_p4, &
          tau_mycorrhiza_p4, &
          vmax_uptake_n_p4, &
          vmax_uptake_p_p4, &
          bnf_base_p4, &
          lambda_est_light_p4, &
          k_est_light_p4, &
          seed_size_p4, &
          beta_soil_flush_p4, &
          beta_soil_senescence_p4, &
          gdd_req_max_p4, &
          k_gdd_dormance_p4, &
          t_air_senescence_p4, &
          min_leaf_age_p4, &
          frac_sapwood_branch_p4, &
          wood_density_p4, &
          k_latosa_p4, &
          k_crtos_p4, &
          k_rtos_p4, &
          k2_fruit_alloc_p4, &
          allom_k1_p4, &
          allom_k2_p4, &
          phi_leaf_min_p4, &
          k_root_p4, &
          k_sapwood_p4, &
          c0_allom_p4, &
          fstore_target_p4, &
          k_root_dist_p4

    NAMELIST /LCTLIB_PFT5_CTL/  &
          sla_p5, &
          sigma_vis_p5, &
          sigma_nir_p5, &
          omega_clumping_p5, &
          crown_shape_factor_p5, &
          cn_leaf_p5, &
          cn_leaf_min_p5, &
          cn_leaf_max_p5, &
          np_leaf_p5, &
          np_leaf_min_p5, &
          np_leaf_max_p5, &
          k0_fn_struc_p5, &
          fn_oth_min_p5, &
          t_jmax_opt_p5, &
          t_jmax_omega_p5, &
          g0_p5, &
          g1_p5, &
          gmin_p5, &
          tau_leaf_p5, &
          tau_fine_root_p5, &
          tau_coarse_root_p5, &
          tau_branch_p5, &
          tau_sap_wood_p5, &
          tau_fruit_p5, &
          tau_seed_litter_p5, &
          tau_seed_est_p5, &
          tau_mycorrhiza_p5, &
          vmax_uptake_n_p5, &
          vmax_uptake_p_p5, &
          bnf_base_p5, &
          lambda_est_light_p5, &
          k_est_light_p5, &
          seed_size_p5, &
          beta_soil_flush_p5, &
          beta_soil_senescence_p5, &
          gdd_req_max_p5, &
          k_gdd_dormance_p5, &
          t_air_senescence_p5, &
          min_leaf_age_p5, &
          frac_sapwood_branch_p5, &
          wood_density_p5, &
          k_latosa_p5, &
          k_crtos_p5, &
          k_rtos_p5, &
          k2_fruit_alloc_p5, &
          allom_k1_p5, &
          allom_k2_p5, &
          phi_leaf_min_p5, &
          k_root_p5, &
          k_sapwood_p5, &
          c0_allom_p5, &
          fstore_target_p5, &
          k_root_dist_p5

    NAMELIST /LCTLIB_PFT6_CTL/  &
          sla_p6, &
          sigma_vis_p6, &
          sigma_nir_p6, &
          omega_clumping_p6, &
          crown_shape_factor_p6, &
          cn_leaf_p6, &
          cn_leaf_min_p6, &
          cn_leaf_max_p6, &
          np_leaf_p6, &
          np_leaf_min_p6, &
          np_leaf_max_p6, &
          k0_fn_struc_p6, &
          fn_oth_min_p6, &
          t_jmax_opt_p6, &
          t_jmax_omega_p6, &
          g0_p6, &
          g1_p6, &
          gmin_p6, &
          tau_leaf_p6, &
          tau_fine_root_p6, &
          tau_coarse_root_p6, &
          tau_branch_p6, &
          tau_sap_wood_p6, &
          tau_fruit_p6, &
          tau_seed_litter_p6, &
          tau_seed_est_p6, &
          tau_mycorrhiza_p6, &
          vmax_uptake_n_p6, &
          vmax_uptake_p_p6, &
          bnf_base_p6, &
          lambda_est_light_p6, &
          k_est_light_p6, &
          seed_size_p6, &
          beta_soil_flush_p6, &
          beta_soil_senescence_p6, &
          gdd_req_max_p6, &
          k_gdd_dormance_p6, &
          t_air_senescence_p6, &
          min_leaf_age_p6, &
          frac_sapwood_branch_p6, &
          wood_density_p6, &
          k_latosa_p6, &
          k_crtos_p6, &
          k_rtos_p6, &
          k2_fruit_alloc_p6, &
          allom_k1_p6, &
          allom_k2_p6, &
          phi_leaf_min_p6, &
          k_root_p6, &
          k_sapwood_p6, &
          c0_allom_p6, &
          fstore_target_p6, &
          k_root_dist_p6

    NAMELIST /LCTLIB_PFT7_CTL/  &
          sla_p7, &
          sigma_vis_p7, &
          sigma_nir_p7, &
          omega_clumping_p7, &
          crown_shape_factor_p7, &
          cn_leaf_p7, &
          cn_leaf_min_p7, &
          cn_leaf_max_p7, &
          np_leaf_p7, &
          np_leaf_min_p7, &
          np_leaf_max_p7, &
          k0_fn_struc_p7, &
          fn_oth_min_p7, &
          t_jmax_opt_p7, &
          t_jmax_omega_p7, &
          g0_p7, &
          g1_p7, &
          gmin_p7, &
          tau_leaf_p7, &
          tau_fine_root_p7, &
          tau_coarse_root_p7, &
          tau_branch_p7, &
          tau_sap_wood_p7, &
          tau_fruit_p7, &
          tau_seed_litter_p7, &
          tau_seed_est_p7, &
          tau_mycorrhiza_p7, &
          vmax_uptake_n_p7, &
          vmax_uptake_p_p7, &
          bnf_base_p7, &
          lambda_est_light_p7, &
          k_est_light_p7, &
          seed_size_p7, &
          beta_soil_flush_p7, &
          beta_soil_senescence_p7, &
          gdd_req_max_p7, &
          k_gdd_dormance_p7, &
          t_air_senescence_p7, &
          min_leaf_age_p7, &
          frac_sapwood_branch_p7, &
          wood_density_p7, &
          k_latosa_p7, &
          k_crtos_p7, &
          k_rtos_p7, &
          k2_fruit_alloc_p7, &
          allom_k1_p7, &
          allom_k2_p7, &
          phi_leaf_min_p7, &
          k_root_p7, &
          k_sapwood_p7, &
          c0_allom_p7, &
          fstore_target_p7, &
          k_root_dist_p7

    NAMELIST /LCTLIB_PFT8_CTL/  &
          sla_p8, &
          sigma_vis_p8, &
          sigma_nir_p8, &
          omega_clumping_p8, &
          crown_shape_factor_p8, &
          cn_leaf_p8, &
          cn_leaf_min_p8, &
          cn_leaf_max_p8, &
          np_leaf_p8, &
          np_leaf_min_p8, &
          np_leaf_max_p8, &
          k0_fn_struc_p8, &
          fn_oth_min_p8, &
          t_jmax_opt_p8, &
          t_jmax_omega_p8, &
          g0_p8, &
          g1_p8, &
          gmin_p8, &
          tau_leaf_p8, &
          tau_fine_root_p8, &
          tau_coarse_root_p8, &
          tau_branch_p8, &
          tau_sap_wood_p8, &
          tau_fruit_p8, &
          tau_seed_litter_p8, &
          tau_seed_est_p8, &
          tau_mycorrhiza_p8, &
          vmax_uptake_n_p8, &
          vmax_uptake_p_p8, &
          bnf_base_p8, &
          lambda_est_light_p8, &
          k_est_light_p8, &
          seed_size_p8, &
          beta_soil_flush_p8, &
          beta_soil_senescence_p8, &
          gdd_req_max_p8, &
          k_gdd_dormance_p8, &
          t_air_senescence_p8, &
          min_leaf_age_p8, &
          frac_sapwood_branch_p8, &
          wood_density_p8, &
          k_latosa_p8, &
          k_crtos_p8, &
          k_rtos_p8, &
          k2_fruit_alloc_p8, &
          allom_k1_p8, &
          allom_k2_p8, &
          phi_leaf_min_p8, &
          k_root_p8, &
          k_sapwood_p8, &
          c0_allom_p8, &
          fstore_target_p8, &
          k_root_dist_p8


    ! variables for reading model-options from namelist
    ! 
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting to read parameter list')

    ! Set defaults
    ! 
    ! PFT1
    sla_p1      = missval
    sigma_vis_p1      = missval
    sigma_nir_p1      = missval
    omega_clumping_p1      = missval
    crown_shape_factor_p1      = missval
    cn_leaf_p1      = missval
    cn_leaf_min_p1      = missval
    cn_leaf_max_p1      = missval
    np_leaf_p1      = missval
    np_leaf_min_p1      = missval
    np_leaf_max_p1      = missval
    k0_fn_struc_p1      = missval
    fn_oth_min_p1      = missval
    t_jmax_opt_p1      = missval
    t_jmax_omega_p1      = missval
    g0_p1      = missval
    g1_p1      = missval
    gmin_p1      = missval
    tau_leaf_p1      = missval
    tau_fine_root_p1      = missval
    tau_coarse_root_p1      = missval
    tau_branch_p1      = missval
    tau_sap_wood_p1      = missval
    tau_fruit_p1      = missval
    tau_seed_litter_p1      = missval
    tau_seed_est_p1      = missval
    tau_mycorrhiza_p1      = missval
    vmax_uptake_n_p1      = missval
    vmax_uptake_p_p1      = missval
    bnf_base_p1      = missval
    lambda_est_light_p1      = missval
    k_est_light_p1      = missval
    seed_size_p1      = missval
    beta_soil_flush_p1      = missval
    beta_soil_senescence_p1      = missval
    gdd_req_max_p1      = missval
    k_gdd_dormance_p1      = missval
    t_air_senescence_p1      = missval
    min_leaf_age_p1      = missval
    frac_sapwood_branch_p1      = missval
    wood_density_p1      = missval
    k_latosa_p1      = missval
    k_crtos_p1      = missval
    k_rtos_p1      = missval
    k2_fruit_alloc_p1      = missval
    allom_k1_p1      = missval
    allom_k2_p1      = missval
    phi_leaf_min_p1      = missval
    k_root_p1      = missval
    k_sapwood_p1      = missval
    c0_allom_p1      = missval
    fstore_target_p1      = missval
    k_root_dist_p1      = missval
    ! PFT2
    sla_p2      = missval
    sigma_vis_p2      = missval
    sigma_nir_p2      = missval
    omega_clumping_p2      = missval
    crown_shape_factor_p2      = missval
    cn_leaf_p2      = missval
    cn_leaf_min_p2      = missval
    cn_leaf_max_p2      = missval
    np_leaf_p2      = missval
    np_leaf_min_p2      = missval
    np_leaf_max_p2      = missval
    k0_fn_struc_p2      = missval
    fn_oth_min_p2      = missval
    t_jmax_opt_p2      = missval
    t_jmax_omega_p2      = missval
    g0_p2      = missval
    g1_p2      = missval
    gmin_p2      = missval
    tau_leaf_p2      = missval
    tau_fine_root_p2      = missval
    tau_coarse_root_p2      = missval
    tau_branch_p2      = missval
    tau_sap_wood_p2      = missval
    tau_fruit_p2      = missval
    tau_seed_litter_p2      = missval
    tau_seed_est_p2      = missval
    tau_mycorrhiza_p2      = missval
    vmax_uptake_n_p2      = missval
    vmax_uptake_p_p2      = missval
    bnf_base_p2      = missval
    lambda_est_light_p2      = missval
    k_est_light_p2      = missval
    seed_size_p2      = missval
    beta_soil_flush_p2      = missval
    beta_soil_senescence_p2      = missval
    gdd_req_max_p2      = missval
    k_gdd_dormance_p2      = missval
    t_air_senescence_p2      = missval
    min_leaf_age_p2      = missval
    frac_sapwood_branch_p2      = missval
    wood_density_p2      = missval
    k_latosa_p2      = missval
    k_crtos_p2      = missval
    k_rtos_p2      = missval
    k2_fruit_alloc_p2      = missval
    allom_k1_p2      = missval
    allom_k2_p2      = missval
    phi_leaf_min_p2      = missval
    k_root_p2      = missval
    k_sapwood_p2      = missval
    c0_allom_p2      = missval
    fstore_target_p2      = missval
    k_root_dist_p2      = missval
    ! PFT3
    sla_p3      = missval
    sigma_vis_p3      = missval
    sigma_nir_p3      = missval
    omega_clumping_p3      = missval
    crown_shape_factor_p3      = missval
    cn_leaf_p3      = missval
    cn_leaf_min_p3      = missval
    cn_leaf_max_p3      = missval
    np_leaf_p3      = missval
    np_leaf_min_p3      = missval
    np_leaf_max_p3      = missval
    k0_fn_struc_p3      = missval
    fn_oth_min_p3      = missval
    t_jmax_opt_p3      = missval
    t_jmax_omega_p3      = missval
    g0_p3      = missval
    g1_p3      = missval
    gmin_p3      = missval
    tau_leaf_p3      = missval
    tau_fine_root_p3      = missval
    tau_coarse_root_p3      = missval
    tau_branch_p3      = missval
    tau_sap_wood_p3      = missval
    tau_fruit_p3      = missval
    tau_seed_litter_p3      = missval
    tau_seed_est_p3      = missval
    tau_mycorrhiza_p3      = missval
    vmax_uptake_n_p3      = missval
    vmax_uptake_p_p3      = missval
    bnf_base_p3      = missval
    lambda_est_light_p3      = missval
    k_est_light_p3      = missval
    seed_size_p3      = missval
    beta_soil_flush_p3      = missval
    beta_soil_senescence_p3      = missval
    gdd_req_max_p3      = missval
    k_gdd_dormance_p3      = missval
    t_air_senescence_p3      = missval
    min_leaf_age_p3      = missval
    frac_sapwood_branch_p3      = missval
    wood_density_p3      = missval
    k_latosa_p3      = missval
    k_crtos_p3      = missval
    k_rtos_p3      = missval
    k2_fruit_alloc_p3      = missval
    allom_k1_p3      = missval
    allom_k2_p3      = missval
    phi_leaf_min_p3      = missval
    k_root_p3      = missval
    k_sapwood_p3      = missval
    c0_allom_p3      = missval
    fstore_target_p3      = missval
    k_root_dist_p3      = missval
    ! PFT4
    sla_p4      = missval
    sigma_vis_p4      = missval
    sigma_nir_p4      = missval
    omega_clumping_p4      = missval
    crown_shape_factor_p4      = missval
    cn_leaf_p4      = missval
    cn_leaf_min_p4      = missval
    cn_leaf_max_p4      = missval
    np_leaf_p4      = missval
    np_leaf_min_p4      = missval
    np_leaf_max_p4      = missval
    k0_fn_struc_p4      = missval
    fn_oth_min_p4      = missval
    t_jmax_opt_p4      = missval
    t_jmax_omega_p4      = missval
    g0_p4      = missval
    g1_p4      = missval
    gmin_p4      = missval
    tau_leaf_p4      = missval
    tau_fine_root_p4      = missval
    tau_coarse_root_p4      = missval
    tau_branch_p4      = missval
    tau_sap_wood_p4      = missval
    tau_fruit_p4      = missval
    tau_seed_litter_p4      = missval
    tau_seed_est_p4      = missval
    tau_mycorrhiza_p4      = missval
    vmax_uptake_n_p4      = missval
    vmax_uptake_p_p4      = missval
    bnf_base_p4      = missval
    lambda_est_light_p4      = missval
    k_est_light_p4      = missval
    seed_size_p4      = missval
    beta_soil_flush_p4      = missval
    beta_soil_senescence_p4      = missval
    gdd_req_max_p4      = missval
    k_gdd_dormance_p4      = missval
    t_air_senescence_p4      = missval
    min_leaf_age_p4      = missval
    frac_sapwood_branch_p4      = missval
    wood_density_p4      = missval
    k_latosa_p4      = missval
    k_crtos_p4      = missval
    k_rtos_p4      = missval
    k2_fruit_alloc_p4      = missval
    allom_k1_p4      = missval
    allom_k2_p4      = missval
    phi_leaf_min_p4      = missval
    k_root_p4      = missval
    k_sapwood_p4      = missval
    c0_allom_p4      = missval
    fstore_target_p4      = missval
    k_root_dist_p4      = missval
    ! PFT5
    sla_p5      = missval
    sigma_vis_p5      = missval
    sigma_nir_p5      = missval
    omega_clumping_p5      = missval
    crown_shape_factor_p5      = missval
    cn_leaf_p5      = missval
    cn_leaf_min_p5      = missval
    cn_leaf_max_p5      = missval
    np_leaf_p5      = missval
    np_leaf_min_p5      = missval
    np_leaf_max_p5      = missval
    k0_fn_struc_p5      = missval
    fn_oth_min_p5      = missval
    t_jmax_opt_p5      = missval
    t_jmax_omega_p5      = missval
    g0_p5      = missval
    g1_p5      = missval
    gmin_p5      = missval
    tau_leaf_p5      = missval
    tau_fine_root_p5      = missval
    tau_coarse_root_p5      = missval
    tau_branch_p5      = missval
    tau_sap_wood_p5      = missval
    tau_fruit_p5      = missval
    tau_seed_litter_p5      = missval
    tau_seed_est_p5      = missval
    tau_mycorrhiza_p5      = missval
    vmax_uptake_n_p5      = missval
    vmax_uptake_p_p5      = missval
    bnf_base_p5      = missval
    lambda_est_light_p5      = missval
    k_est_light_p5      = missval
    seed_size_p5      = missval
    beta_soil_flush_p5      = missval
    beta_soil_senescence_p5      = missval
    gdd_req_max_p5      = missval
    k_gdd_dormance_p5      = missval
    t_air_senescence_p5      = missval
    min_leaf_age_p5      = missval
    frac_sapwood_branch_p5      = missval
    wood_density_p5      = missval
    k_latosa_p5      = missval
    k_crtos_p5      = missval
    k_rtos_p5      = missval
    k2_fruit_alloc_p5      = missval
    allom_k1_p5      = missval
    allom_k2_p5      = missval
    phi_leaf_min_p5      = missval
    k_root_p5      = missval
    k_sapwood_p5      = missval
    c0_allom_p5      = missval
    fstore_target_p5      = missval
    k_root_dist_p5      = missval
    ! PFT6
    sla_p6      = missval
    sigma_vis_p6      = missval
    sigma_nir_p6      = missval
    omega_clumping_p6      = missval
    crown_shape_factor_p6      = missval
    cn_leaf_p6      = missval
    cn_leaf_min_p6      = missval
    cn_leaf_max_p6      = missval
    np_leaf_p6      = missval
    np_leaf_min_p6      = missval
    np_leaf_max_p6      = missval
    k0_fn_struc_p6      = missval
    fn_oth_min_p6      = missval
    t_jmax_opt_p6      = missval
    t_jmax_omega_p6      = missval
    g0_p6      = missval
    g1_p6      = missval
    gmin_p6      = missval
    tau_leaf_p6      = missval
    tau_fine_root_p6      = missval
    tau_coarse_root_p6      = missval
    tau_branch_p6      = missval
    tau_sap_wood_p6      = missval
    tau_fruit_p6      = missval
    tau_seed_litter_p6      = missval
    tau_seed_est_p6      = missval
    tau_mycorrhiza_p6      = missval
    vmax_uptake_n_p6      = missval
    vmax_uptake_p_p6      = missval
    bnf_base_p6      = missval
    lambda_est_light_p6      = missval
    k_est_light_p6      = missval
    seed_size_p6      = missval
    beta_soil_flush_p6      = missval
    beta_soil_senescence_p6      = missval
    gdd_req_max_p6      = missval
    k_gdd_dormance_p6      = missval
    t_air_senescence_p6      = missval
    min_leaf_age_p6      = missval
    frac_sapwood_branch_p6      = missval
    wood_density_p6      = missval
    k_latosa_p6      = missval
    k_crtos_p6      = missval
    k_rtos_p6      = missval
    k2_fruit_alloc_p6      = missval
    allom_k1_p6      = missval
    allom_k2_p6      = missval
    phi_leaf_min_p6      = missval
    k_root_p6      = missval
    k_sapwood_p6      = missval
    c0_allom_p6      = missval
    fstore_target_p6      = missval
    k_root_dist_p6      = missval
    ! PFT7
    sla_p7      = missval
    sigma_vis_p7      = missval
    sigma_nir_p7      = missval
    omega_clumping_p7      = missval
    crown_shape_factor_p7      = missval
    cn_leaf_p7      = missval
    cn_leaf_min_p7      = missval
    cn_leaf_max_p7      = missval
    np_leaf_p7      = missval
    np_leaf_min_p7      = missval
    np_leaf_max_p7      = missval
    k0_fn_struc_p7      = missval
    fn_oth_min_p7      = missval
    t_jmax_opt_p7      = missval
    t_jmax_omega_p7      = missval
    g0_p7      = missval
    g1_p7      = missval
    gmin_p7      = missval
    tau_leaf_p7      = missval
    tau_fine_root_p7      = missval
    tau_coarse_root_p7      = missval
    tau_branch_p7      = missval
    tau_sap_wood_p7      = missval
    tau_fruit_p7      = missval
    tau_seed_litter_p7      = missval
    tau_seed_est_p7      = missval
    tau_mycorrhiza_p7      = missval
    vmax_uptake_n_p7      = missval
    vmax_uptake_p_p7      = missval
    bnf_base_p7      = missval
    lambda_est_light_p7      = missval
    k_est_light_p7      = missval
    seed_size_p7      = missval
    beta_soil_flush_p7      = missval
    beta_soil_senescence_p7      = missval
    gdd_req_max_p7      = missval
    k_gdd_dormance_p7      = missval
    t_air_senescence_p7      = missval
    min_leaf_age_p7      = missval
    frac_sapwood_branch_p7      = missval
    wood_density_p7      = missval
    k_latosa_p7      = missval
    k_crtos_p7      = missval
    k_rtos_p7      = missval
    k2_fruit_alloc_p7      = missval
    allom_k1_p7      = missval
    allom_k2_p7      = missval
    phi_leaf_min_p7      = missval
    k_root_p7      = missval
    k_sapwood_p7      = missval
    c0_allom_p7      = missval
    fstore_target_p7      = missval
    k_root_dist_p7      = missval
    ! PFT8
    sla_p8      = missval
    sigma_vis_p8      = missval
    sigma_nir_p8      = missval
    omega_clumping_p8      = missval
    crown_shape_factor_p8      = missval
    cn_leaf_p8      = missval
    cn_leaf_min_p8      = missval
    cn_leaf_max_p8      = missval
    np_leaf_p8      = missval
    np_leaf_min_p8      = missval
    np_leaf_max_p8      = missval
    k0_fn_struc_p8      = missval
    fn_oth_min_p8      = missval
    t_jmax_opt_p8      = missval
    t_jmax_omega_p8      = missval
    g0_p8      = missval
    g1_p8      = missval
    gmin_p8      = missval
    tau_leaf_p8      = missval
    tau_fine_root_p8      = missval
    tau_coarse_root_p8      = missval
    tau_branch_p8      = missval
    tau_sap_wood_p8      = missval
    tau_fruit_p8      = missval
    tau_seed_litter_p8      = missval
    tau_seed_est_p8      = missval
    tau_mycorrhiza_p8      = missval
    vmax_uptake_n_p8      = missval
    vmax_uptake_p_p8      = missval
    bnf_base_p8      = missval
    lambda_est_light_p8      = missval
    k_est_light_p8      = missval
    seed_size_p8      = missval
    beta_soil_flush_p8      = missval
    beta_soil_senescence_p8      = missval
    gdd_req_max_p8      = missval
    k_gdd_dormance_p8      = missval
    t_air_senescence_p8      = missval
    min_leaf_age_p8      = missval
    frac_sapwood_branch_p8      = missval
    wood_density_p8      = missval
    k_latosa_p8      = missval
    k_crtos_p8      = missval
    k_rtos_p8      = missval
    k2_fruit_alloc_p8      = missval
    allom_k1_p8      = missval
    allom_k2_p8      = missval
    phi_leaf_min_p8      = missval
    k_root_p8      = missval
    k_sapwood_p8      = missval
    c0_allom_p8      = missval
    fstore_target_p8      = missval
    k_root_dist_p8      = missval


    ! read the namelist (parameters of all PFT are stored in the same namelist file using param_pft1_from_file%filename here)
    ! 
    ! &LCTLIB_PFT1_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT1_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT1_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT2_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT2_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT2_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT3_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT3_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT3_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT4_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT4_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT4_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT5_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT5_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT5_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT6_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT6_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT6_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT7_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT7_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT7_CTL)
    CALL close_nml(nml_handler)
    ! &LCTLIB_PFT8_CTL
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('LCTLIB_PFT8_CTL', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, LCTLIB_PFT8_CTL)
    CALL close_nml(nml_handler)


    ! pass values as read from file
    ! 
    ! LCTLIB_PFT1_CTL
    param_pft1_from_file%sla   = sla_p1
    param_pft1_from_file%sigma_vis   = sigma_vis_p1
    param_pft1_from_file%sigma_nir   = sigma_nir_p1
    param_pft1_from_file%omega_clumping   = omega_clumping_p1
    param_pft1_from_file%crown_shape_factor   = crown_shape_factor_p1
    param_pft1_from_file%cn_leaf   = cn_leaf_p1
    param_pft1_from_file%cn_leaf_min   = cn_leaf_min_p1
    param_pft1_from_file%cn_leaf_max   = cn_leaf_max_p1
    param_pft1_from_file%np_leaf   = np_leaf_p1
    param_pft1_from_file%np_leaf_min   = np_leaf_min_p1
    param_pft1_from_file%np_leaf_max   = np_leaf_max_p1
    param_pft1_from_file%k0_fn_struc   = k0_fn_struc_p1
    param_pft1_from_file%fn_oth_min   = fn_oth_min_p1
    param_pft1_from_file%t_jmax_opt   = t_jmax_opt_p1
    param_pft1_from_file%t_jmax_omega   = t_jmax_omega_p1
    param_pft1_from_file%g0   = g0_p1
    param_pft1_from_file%g1   = g1_p1
    param_pft1_from_file%gmin   = gmin_p1
    param_pft1_from_file%tau_leaf   = tau_leaf_p1
    param_pft1_from_file%tau_fine_root   = tau_fine_root_p1
    param_pft1_from_file%tau_coarse_root   = tau_coarse_root_p1
    param_pft1_from_file%tau_branch   = tau_branch_p1
    param_pft1_from_file%tau_sap_wood   = tau_sap_wood_p1
    param_pft1_from_file%tau_fruit   = tau_fruit_p1
    param_pft1_from_file%tau_seed_litter   = tau_seed_litter_p1
    param_pft1_from_file%tau_seed_est   = tau_seed_est_p1
    param_pft1_from_file%tau_mycorrhiza   = tau_mycorrhiza_p1
    param_pft1_from_file%vmax_uptake_n   = vmax_uptake_n_p1
    param_pft1_from_file%vmax_uptake_p   = vmax_uptake_p_p1
    param_pft1_from_file%bnf_base   = bnf_base_p1
    param_pft1_from_file%lambda_est_light   = lambda_est_light_p1
    param_pft1_from_file%k_est_light   = k_est_light_p1
    param_pft1_from_file%seed_size   = seed_size_p1
    param_pft1_from_file%beta_soil_flush   = beta_soil_flush_p1
    param_pft1_from_file%beta_soil_senescence   = beta_soil_senescence_p1
    param_pft1_from_file%gdd_req_max   = gdd_req_max_p1
    param_pft1_from_file%k_gdd_dormance   = k_gdd_dormance_p1
    param_pft1_from_file%t_air_senescence   = t_air_senescence_p1
    param_pft1_from_file%min_leaf_age   = min_leaf_age_p1
    param_pft1_from_file%frac_sapwood_branch   = frac_sapwood_branch_p1
    param_pft1_from_file%wood_density   = wood_density_p1
    param_pft1_from_file%k_latosa   = k_latosa_p1
    param_pft1_from_file%k_crtos   = k_crtos_p1
    param_pft1_from_file%k_rtos   = k_rtos_p1
    param_pft1_from_file%k2_fruit_alloc   = k2_fruit_alloc_p1
    param_pft1_from_file%allom_k1   = allom_k1_p1
    param_pft1_from_file%allom_k2   = allom_k2_p1
    param_pft1_from_file%phi_leaf_min   = phi_leaf_min_p1
    param_pft1_from_file%k_root   = k_root_p1
    param_pft1_from_file%k_sapwood   = k_sapwood_p1
    param_pft1_from_file%c0_allom   = c0_allom_p1
    param_pft1_from_file%fstore_target   = fstore_target_p1
    param_pft1_from_file%k_root_dist   = k_root_dist_p1
    ! LCTLIB_PFT2_CTL
    param_pft2_from_file%sla   = sla_p2
    param_pft2_from_file%sigma_vis   = sigma_vis_p2
    param_pft2_from_file%sigma_nir   = sigma_nir_p2
    param_pft2_from_file%omega_clumping   = omega_clumping_p2
    param_pft2_from_file%crown_shape_factor   = crown_shape_factor_p2
    param_pft2_from_file%cn_leaf   = cn_leaf_p2
    param_pft2_from_file%cn_leaf_min   = cn_leaf_min_p2
    param_pft2_from_file%cn_leaf_max   = cn_leaf_max_p2
    param_pft2_from_file%np_leaf   = np_leaf_p2
    param_pft2_from_file%np_leaf_min   = np_leaf_min_p2
    param_pft2_from_file%np_leaf_max   = np_leaf_max_p2
    param_pft2_from_file%k0_fn_struc   = k0_fn_struc_p2
    param_pft2_from_file%fn_oth_min   = fn_oth_min_p2
    param_pft2_from_file%t_jmax_opt   = t_jmax_opt_p2
    param_pft2_from_file%t_jmax_omega   = t_jmax_omega_p2
    param_pft2_from_file%g0   = g0_p2
    param_pft2_from_file%g1   = g1_p2
    param_pft2_from_file%gmin   = gmin_p2
    param_pft2_from_file%tau_leaf   = tau_leaf_p2
    param_pft2_from_file%tau_fine_root   = tau_fine_root_p2
    param_pft2_from_file%tau_coarse_root   = tau_coarse_root_p2
    param_pft2_from_file%tau_branch   = tau_branch_p2
    param_pft2_from_file%tau_sap_wood   = tau_sap_wood_p2
    param_pft2_from_file%tau_fruit   = tau_fruit_p2
    param_pft2_from_file%tau_seed_litter   = tau_seed_litter_p2
    param_pft2_from_file%tau_seed_est   = tau_seed_est_p2
    param_pft2_from_file%tau_mycorrhiza   = tau_mycorrhiza_p2
    param_pft2_from_file%vmax_uptake_n   = vmax_uptake_n_p2
    param_pft2_from_file%vmax_uptake_p   = vmax_uptake_p_p2
    param_pft2_from_file%bnf_base   = bnf_base_p2
    param_pft2_from_file%lambda_est_light   = lambda_est_light_p2
    param_pft2_from_file%k_est_light   = k_est_light_p2
    param_pft2_from_file%seed_size   = seed_size_p2
    param_pft2_from_file%beta_soil_flush   = beta_soil_flush_p2
    param_pft2_from_file%beta_soil_senescence   = beta_soil_senescence_p2
    param_pft2_from_file%gdd_req_max   = gdd_req_max_p2
    param_pft2_from_file%k_gdd_dormance   = k_gdd_dormance_p2
    param_pft2_from_file%t_air_senescence   = t_air_senescence_p2
    param_pft2_from_file%min_leaf_age   = min_leaf_age_p2
    param_pft2_from_file%frac_sapwood_branch   = frac_sapwood_branch_p2
    param_pft2_from_file%wood_density   = wood_density_p2
    param_pft2_from_file%k_latosa   = k_latosa_p2
    param_pft2_from_file%k_crtos   = k_crtos_p2
    param_pft2_from_file%k_rtos   = k_rtos_p2
    param_pft2_from_file%k2_fruit_alloc   = k2_fruit_alloc_p2
    param_pft2_from_file%allom_k1   = allom_k1_p2
    param_pft2_from_file%allom_k2   = allom_k2_p2
    param_pft2_from_file%phi_leaf_min   = phi_leaf_min_p2
    param_pft2_from_file%k_root   = k_root_p2
    param_pft2_from_file%k_sapwood   = k_sapwood_p2
    param_pft2_from_file%c0_allom   = c0_allom_p2
    param_pft2_from_file%fstore_target   = fstore_target_p2
    param_pft2_from_file%k_root_dist   = k_root_dist_p2
    ! LCTLIB_PFT3_CTL
    param_pft3_from_file%sla   = sla_p3
    param_pft3_from_file%sigma_vis   = sigma_vis_p3
    param_pft3_from_file%sigma_nir   = sigma_nir_p3
    param_pft3_from_file%omega_clumping   = omega_clumping_p3
    param_pft3_from_file%crown_shape_factor   = crown_shape_factor_p3
    param_pft3_from_file%cn_leaf   = cn_leaf_p3
    param_pft3_from_file%cn_leaf_min   = cn_leaf_min_p3
    param_pft3_from_file%cn_leaf_max   = cn_leaf_max_p3
    param_pft3_from_file%np_leaf   = np_leaf_p3
    param_pft3_from_file%np_leaf_min   = np_leaf_min_p3
    param_pft3_from_file%np_leaf_max   = np_leaf_max_p3
    param_pft3_from_file%k0_fn_struc   = k0_fn_struc_p3
    param_pft3_from_file%fn_oth_min   = fn_oth_min_p3
    param_pft3_from_file%t_jmax_opt   = t_jmax_opt_p3
    param_pft3_from_file%t_jmax_omega   = t_jmax_omega_p3
    param_pft3_from_file%g0   = g0_p3
    param_pft3_from_file%g1   = g1_p3
    param_pft3_from_file%gmin   = gmin_p3
    param_pft3_from_file%tau_leaf   = tau_leaf_p3
    param_pft3_from_file%tau_fine_root   = tau_fine_root_p3
    param_pft3_from_file%tau_coarse_root   = tau_coarse_root_p3
    param_pft3_from_file%tau_branch   = tau_branch_p3
    param_pft3_from_file%tau_sap_wood   = tau_sap_wood_p3
    param_pft3_from_file%tau_fruit   = tau_fruit_p3
    param_pft3_from_file%tau_seed_litter   = tau_seed_litter_p3
    param_pft3_from_file%tau_seed_est   = tau_seed_est_p3
    param_pft3_from_file%tau_mycorrhiza   = tau_mycorrhiza_p3
    param_pft3_from_file%vmax_uptake_n   = vmax_uptake_n_p3
    param_pft3_from_file%vmax_uptake_p   = vmax_uptake_p_p3
    param_pft3_from_file%bnf_base   = bnf_base_p3
    param_pft3_from_file%lambda_est_light   = lambda_est_light_p3
    param_pft3_from_file%k_est_light   = k_est_light_p3
    param_pft3_from_file%seed_size   = seed_size_p3
    param_pft3_from_file%beta_soil_flush   = beta_soil_flush_p3
    param_pft3_from_file%beta_soil_senescence   = beta_soil_senescence_p3
    param_pft3_from_file%gdd_req_max   = gdd_req_max_p3
    param_pft3_from_file%k_gdd_dormance   = k_gdd_dormance_p3
    param_pft3_from_file%t_air_senescence   = t_air_senescence_p3
    param_pft3_from_file%min_leaf_age   = min_leaf_age_p3
    param_pft3_from_file%frac_sapwood_branch   = frac_sapwood_branch_p3
    param_pft3_from_file%wood_density   = wood_density_p3
    param_pft3_from_file%k_latosa   = k_latosa_p3
    param_pft3_from_file%k_crtos   = k_crtos_p3
    param_pft3_from_file%k_rtos   = k_rtos_p3
    param_pft3_from_file%k2_fruit_alloc   = k2_fruit_alloc_p3
    param_pft3_from_file%allom_k1   = allom_k1_p3
    param_pft3_from_file%allom_k2   = allom_k2_p3
    param_pft3_from_file%phi_leaf_min   = phi_leaf_min_p3
    param_pft3_from_file%k_root   = k_root_p3
    param_pft3_from_file%k_sapwood   = k_sapwood_p3
    param_pft3_from_file%c0_allom   = c0_allom_p3
    param_pft3_from_file%fstore_target   = fstore_target_p3
    param_pft3_from_file%k_root_dist   = k_root_dist_p3
    ! LCTLIB_PFT4_CTL
    param_pft4_from_file%sla   = sla_p4
    param_pft4_from_file%sigma_vis   = sigma_vis_p4
    param_pft4_from_file%sigma_nir   = sigma_nir_p4
    param_pft4_from_file%omega_clumping   = omega_clumping_p4
    param_pft4_from_file%crown_shape_factor   = crown_shape_factor_p4
    param_pft4_from_file%cn_leaf   = cn_leaf_p4
    param_pft4_from_file%cn_leaf_min   = cn_leaf_min_p4
    param_pft4_from_file%cn_leaf_max   = cn_leaf_max_p4
    param_pft4_from_file%np_leaf   = np_leaf_p4
    param_pft4_from_file%np_leaf_min   = np_leaf_min_p4
    param_pft4_from_file%np_leaf_max   = np_leaf_max_p4
    param_pft4_from_file%k0_fn_struc   = k0_fn_struc_p4
    param_pft4_from_file%fn_oth_min   = fn_oth_min_p4
    param_pft4_from_file%t_jmax_opt   = t_jmax_opt_p4
    param_pft4_from_file%t_jmax_omega   = t_jmax_omega_p4
    param_pft4_from_file%g0   = g0_p4
    param_pft4_from_file%g1   = g1_p4
    param_pft4_from_file%gmin   = gmin_p4
    param_pft4_from_file%tau_leaf   = tau_leaf_p4
    param_pft4_from_file%tau_fine_root   = tau_fine_root_p4
    param_pft4_from_file%tau_coarse_root   = tau_coarse_root_p4
    param_pft4_from_file%tau_branch   = tau_branch_p4
    param_pft4_from_file%tau_sap_wood   = tau_sap_wood_p4
    param_pft4_from_file%tau_fruit   = tau_fruit_p4
    param_pft4_from_file%tau_seed_litter   = tau_seed_litter_p4
    param_pft4_from_file%tau_seed_est   = tau_seed_est_p4
    param_pft4_from_file%tau_mycorrhiza   = tau_mycorrhiza_p4
    param_pft4_from_file%vmax_uptake_n   = vmax_uptake_n_p4
    param_pft4_from_file%vmax_uptake_p   = vmax_uptake_p_p4
    param_pft4_from_file%bnf_base   = bnf_base_p4
    param_pft4_from_file%lambda_est_light   = lambda_est_light_p4
    param_pft4_from_file%k_est_light   = k_est_light_p4
    param_pft4_from_file%seed_size   = seed_size_p4
    param_pft4_from_file%beta_soil_flush   = beta_soil_flush_p4
    param_pft4_from_file%beta_soil_senescence   = beta_soil_senescence_p4
    param_pft4_from_file%gdd_req_max   = gdd_req_max_p4
    param_pft4_from_file%k_gdd_dormance   = k_gdd_dormance_p4
    param_pft4_from_file%t_air_senescence   = t_air_senescence_p4
    param_pft4_from_file%min_leaf_age   = min_leaf_age_p4
    param_pft4_from_file%frac_sapwood_branch   = frac_sapwood_branch_p4
    param_pft4_from_file%wood_density   = wood_density_p4
    param_pft4_from_file%k_latosa   = k_latosa_p4
    param_pft4_from_file%k_crtos   = k_crtos_p4
    param_pft4_from_file%k_rtos   = k_rtos_p4
    param_pft4_from_file%k2_fruit_alloc   = k2_fruit_alloc_p4
    param_pft4_from_file%allom_k1   = allom_k1_p4
    param_pft4_from_file%allom_k2   = allom_k2_p4
    param_pft4_from_file%phi_leaf_min   = phi_leaf_min_p4
    param_pft4_from_file%k_root   = k_root_p4
    param_pft4_from_file%k_sapwood   = k_sapwood_p4
    param_pft4_from_file%c0_allom   = c0_allom_p4
    param_pft4_from_file%fstore_target   = fstore_target_p4
    param_pft4_from_file%k_root_dist   = k_root_dist_p4
    ! LCTLIB_PFT5_CTL
    param_pft5_from_file%sla   = sla_p5
    param_pft5_from_file%sigma_vis   = sigma_vis_p5
    param_pft5_from_file%sigma_nir   = sigma_nir_p5
    param_pft5_from_file%omega_clumping   = omega_clumping_p5
    param_pft5_from_file%crown_shape_factor   = crown_shape_factor_p5
    param_pft5_from_file%cn_leaf   = cn_leaf_p5
    param_pft5_from_file%cn_leaf_min   = cn_leaf_min_p5
    param_pft5_from_file%cn_leaf_max   = cn_leaf_max_p5
    param_pft5_from_file%np_leaf   = np_leaf_p5
    param_pft5_from_file%np_leaf_min   = np_leaf_min_p5
    param_pft5_from_file%np_leaf_max   = np_leaf_max_p5
    param_pft5_from_file%k0_fn_struc   = k0_fn_struc_p5
    param_pft5_from_file%fn_oth_min   = fn_oth_min_p5
    param_pft5_from_file%t_jmax_opt   = t_jmax_opt_p5
    param_pft5_from_file%t_jmax_omega   = t_jmax_omega_p5
    param_pft5_from_file%g0   = g0_p5
    param_pft5_from_file%g1   = g1_p5
    param_pft5_from_file%gmin   = gmin_p5
    param_pft5_from_file%tau_leaf   = tau_leaf_p5
    param_pft5_from_file%tau_fine_root   = tau_fine_root_p5
    param_pft5_from_file%tau_coarse_root   = tau_coarse_root_p5
    param_pft5_from_file%tau_branch   = tau_branch_p5
    param_pft5_from_file%tau_sap_wood   = tau_sap_wood_p5
    param_pft5_from_file%tau_fruit   = tau_fruit_p5
    param_pft5_from_file%tau_seed_litter   = tau_seed_litter_p5
    param_pft5_from_file%tau_seed_est   = tau_seed_est_p5
    param_pft5_from_file%tau_mycorrhiza   = tau_mycorrhiza_p5
    param_pft5_from_file%vmax_uptake_n   = vmax_uptake_n_p5
    param_pft5_from_file%vmax_uptake_p   = vmax_uptake_p_p5
    param_pft5_from_file%bnf_base   = bnf_base_p5
    param_pft5_from_file%lambda_est_light   = lambda_est_light_p5
    param_pft5_from_file%k_est_light   = k_est_light_p5
    param_pft5_from_file%seed_size   = seed_size_p5
    param_pft5_from_file%beta_soil_flush   = beta_soil_flush_p5
    param_pft5_from_file%beta_soil_senescence   = beta_soil_senescence_p5
    param_pft5_from_file%gdd_req_max   = gdd_req_max_p5
    param_pft5_from_file%k_gdd_dormance   = k_gdd_dormance_p5
    param_pft5_from_file%t_air_senescence   = t_air_senescence_p5
    param_pft5_from_file%min_leaf_age   = min_leaf_age_p5
    param_pft5_from_file%frac_sapwood_branch   = frac_sapwood_branch_p5
    param_pft5_from_file%wood_density   = wood_density_p5
    param_pft5_from_file%k_latosa   = k_latosa_p5
    param_pft5_from_file%k_crtos   = k_crtos_p5
    param_pft5_from_file%k_rtos   = k_rtos_p5
    param_pft5_from_file%k2_fruit_alloc   = k2_fruit_alloc_p5
    param_pft5_from_file%allom_k1   = allom_k1_p5
    param_pft5_from_file%allom_k2   = allom_k2_p5
    param_pft5_from_file%phi_leaf_min   = phi_leaf_min_p5
    param_pft5_from_file%k_root   = k_root_p5
    param_pft5_from_file%k_sapwood   = k_sapwood_p5
    param_pft5_from_file%c0_allom   = c0_allom_p5
    param_pft5_from_file%fstore_target   = fstore_target_p5
    param_pft5_from_file%k_root_dist   = k_root_dist_p5
    ! LCTLIB_PFT6_CTL
    param_pft6_from_file%sla   = sla_p6
    param_pft6_from_file%sigma_vis   = sigma_vis_p6
    param_pft6_from_file%sigma_nir   = sigma_nir_p6
    param_pft6_from_file%omega_clumping   = omega_clumping_p6
    param_pft6_from_file%crown_shape_factor   = crown_shape_factor_p6
    param_pft6_from_file%cn_leaf   = cn_leaf_p6
    param_pft6_from_file%cn_leaf_min   = cn_leaf_min_p6
    param_pft6_from_file%cn_leaf_max   = cn_leaf_max_p6
    param_pft6_from_file%np_leaf   = np_leaf_p6
    param_pft6_from_file%np_leaf_min   = np_leaf_min_p6
    param_pft6_from_file%np_leaf_max   = np_leaf_max_p6
    param_pft6_from_file%k0_fn_struc   = k0_fn_struc_p6
    param_pft6_from_file%fn_oth_min   = fn_oth_min_p6
    param_pft6_from_file%t_jmax_opt   = t_jmax_opt_p6
    param_pft6_from_file%t_jmax_omega   = t_jmax_omega_p6
    param_pft6_from_file%g0   = g0_p6
    param_pft6_from_file%g1   = g1_p6
    param_pft6_from_file%gmin   = gmin_p6
    param_pft6_from_file%tau_leaf   = tau_leaf_p6
    param_pft6_from_file%tau_fine_root   = tau_fine_root_p6
    param_pft6_from_file%tau_coarse_root   = tau_coarse_root_p6
    param_pft6_from_file%tau_branch   = tau_branch_p6
    param_pft6_from_file%tau_sap_wood   = tau_sap_wood_p6
    param_pft6_from_file%tau_fruit   = tau_fruit_p6
    param_pft6_from_file%tau_seed_litter   = tau_seed_litter_p6
    param_pft6_from_file%tau_seed_est   = tau_seed_est_p6
    param_pft6_from_file%tau_mycorrhiza   = tau_mycorrhiza_p6
    param_pft6_from_file%vmax_uptake_n   = vmax_uptake_n_p6
    param_pft6_from_file%vmax_uptake_p   = vmax_uptake_p_p6
    param_pft6_from_file%bnf_base   = bnf_base_p6
    param_pft6_from_file%lambda_est_light   = lambda_est_light_p6
    param_pft6_from_file%k_est_light   = k_est_light_p6
    param_pft6_from_file%seed_size   = seed_size_p6
    param_pft6_from_file%beta_soil_flush   = beta_soil_flush_p6
    param_pft6_from_file%beta_soil_senescence   = beta_soil_senescence_p6
    param_pft6_from_file%gdd_req_max   = gdd_req_max_p6
    param_pft6_from_file%k_gdd_dormance   = k_gdd_dormance_p6
    param_pft6_from_file%t_air_senescence   = t_air_senescence_p6
    param_pft6_from_file%min_leaf_age   = min_leaf_age_p6
    param_pft6_from_file%frac_sapwood_branch   = frac_sapwood_branch_p6
    param_pft6_from_file%wood_density   = wood_density_p6
    param_pft6_from_file%k_latosa   = k_latosa_p6
    param_pft6_from_file%k_crtos   = k_crtos_p6
    param_pft6_from_file%k_rtos   = k_rtos_p6
    param_pft6_from_file%k2_fruit_alloc   = k2_fruit_alloc_p6
    param_pft6_from_file%allom_k1   = allom_k1_p6
    param_pft6_from_file%allom_k2   = allom_k2_p6
    param_pft6_from_file%phi_leaf_min   = phi_leaf_min_p6
    param_pft6_from_file%k_root   = k_root_p6
    param_pft6_from_file%k_sapwood   = k_sapwood_p6
    param_pft6_from_file%c0_allom   = c0_allom_p6
    param_pft6_from_file%fstore_target   = fstore_target_p6
    param_pft6_from_file%k_root_dist   = k_root_dist_p6
    ! LCTLIB_PFT7_CTL
    param_pft7_from_file%sla   = sla_p7
    param_pft7_from_file%sigma_vis   = sigma_vis_p7
    param_pft7_from_file%sigma_nir   = sigma_nir_p7
    param_pft7_from_file%omega_clumping   = omega_clumping_p7
    param_pft7_from_file%crown_shape_factor   = crown_shape_factor_p7
    param_pft7_from_file%cn_leaf   = cn_leaf_p7
    param_pft7_from_file%cn_leaf_min   = cn_leaf_min_p7
    param_pft7_from_file%cn_leaf_max   = cn_leaf_max_p7
    param_pft7_from_file%np_leaf   = np_leaf_p7
    param_pft7_from_file%np_leaf_min   = np_leaf_min_p7
    param_pft7_from_file%np_leaf_max   = np_leaf_max_p7
    param_pft7_from_file%k0_fn_struc   = k0_fn_struc_p7
    param_pft7_from_file%fn_oth_min   = fn_oth_min_p7
    param_pft7_from_file%t_jmax_opt   = t_jmax_opt_p7
    param_pft7_from_file%t_jmax_omega   = t_jmax_omega_p7
    param_pft7_from_file%g0   = g0_p7
    param_pft7_from_file%g1   = g1_p7
    param_pft7_from_file%gmin   = gmin_p7
    param_pft7_from_file%tau_leaf   = tau_leaf_p7
    param_pft7_from_file%tau_fine_root   = tau_fine_root_p7
    param_pft7_from_file%tau_coarse_root   = tau_coarse_root_p7
    param_pft7_from_file%tau_branch   = tau_branch_p7
    param_pft7_from_file%tau_sap_wood   = tau_sap_wood_p7
    param_pft7_from_file%tau_fruit   = tau_fruit_p7
    param_pft7_from_file%tau_seed_litter   = tau_seed_litter_p7
    param_pft7_from_file%tau_seed_est   = tau_seed_est_p7
    param_pft7_from_file%tau_mycorrhiza   = tau_mycorrhiza_p7
    param_pft7_from_file%vmax_uptake_n   = vmax_uptake_n_p7
    param_pft7_from_file%vmax_uptake_p   = vmax_uptake_p_p7
    param_pft7_from_file%bnf_base   = bnf_base_p7
    param_pft7_from_file%lambda_est_light   = lambda_est_light_p7
    param_pft7_from_file%k_est_light   = k_est_light_p7
    param_pft7_from_file%seed_size   = seed_size_p7
    param_pft7_from_file%beta_soil_flush   = beta_soil_flush_p7
    param_pft7_from_file%beta_soil_senescence   = beta_soil_senescence_p7
    param_pft7_from_file%gdd_req_max   = gdd_req_max_p7
    param_pft7_from_file%k_gdd_dormance   = k_gdd_dormance_p7
    param_pft7_from_file%t_air_senescence   = t_air_senescence_p7
    param_pft7_from_file%min_leaf_age   = min_leaf_age_p7
    param_pft7_from_file%frac_sapwood_branch   = frac_sapwood_branch_p7
    param_pft7_from_file%wood_density   = wood_density_p7
    param_pft7_from_file%k_latosa   = k_latosa_p7
    param_pft7_from_file%k_crtos   = k_crtos_p7
    param_pft7_from_file%k_rtos   = k_rtos_p7
    param_pft7_from_file%k2_fruit_alloc   = k2_fruit_alloc_p7
    param_pft7_from_file%allom_k1   = allom_k1_p7
    param_pft7_from_file%allom_k2   = allom_k2_p7
    param_pft7_from_file%phi_leaf_min   = phi_leaf_min_p7
    param_pft7_from_file%k_root   = k_root_p7
    param_pft7_from_file%k_sapwood   = k_sapwood_p7
    param_pft7_from_file%c0_allom   = c0_allom_p7
    param_pft7_from_file%fstore_target   = fstore_target_p7
    param_pft7_from_file%k_root_dist   = k_root_dist_p7
    ! LCTLIB_PFT8_CTL
    param_pft8_from_file%sla   = sla_p8
    param_pft8_from_file%sigma_vis   = sigma_vis_p8
    param_pft8_from_file%sigma_nir   = sigma_nir_p8
    param_pft8_from_file%omega_clumping   = omega_clumping_p8
    param_pft8_from_file%crown_shape_factor   = crown_shape_factor_p8
    param_pft8_from_file%cn_leaf   = cn_leaf_p8
    param_pft8_from_file%cn_leaf_min   = cn_leaf_min_p8
    param_pft8_from_file%cn_leaf_max   = cn_leaf_max_p8
    param_pft8_from_file%np_leaf   = np_leaf_p8
    param_pft8_from_file%np_leaf_min   = np_leaf_min_p8
    param_pft8_from_file%np_leaf_max   = np_leaf_max_p8
    param_pft8_from_file%k0_fn_struc   = k0_fn_struc_p8
    param_pft8_from_file%fn_oth_min   = fn_oth_min_p8
    param_pft8_from_file%t_jmax_opt   = t_jmax_opt_p8
    param_pft8_from_file%t_jmax_omega   = t_jmax_omega_p8
    param_pft8_from_file%g0   = g0_p8
    param_pft8_from_file%g1   = g1_p8
    param_pft8_from_file%gmin   = gmin_p8
    param_pft8_from_file%tau_leaf   = tau_leaf_p8
    param_pft8_from_file%tau_fine_root   = tau_fine_root_p8
    param_pft8_from_file%tau_coarse_root   = tau_coarse_root_p8
    param_pft8_from_file%tau_branch   = tau_branch_p8
    param_pft8_from_file%tau_sap_wood   = tau_sap_wood_p8
    param_pft8_from_file%tau_fruit   = tau_fruit_p8
    param_pft8_from_file%tau_seed_litter   = tau_seed_litter_p8
    param_pft8_from_file%tau_seed_est   = tau_seed_est_p8
    param_pft8_from_file%tau_mycorrhiza   = tau_mycorrhiza_p8
    param_pft8_from_file%vmax_uptake_n   = vmax_uptake_n_p8
    param_pft8_from_file%vmax_uptake_p   = vmax_uptake_p_p8
    param_pft8_from_file%bnf_base   = bnf_base_p8
    param_pft8_from_file%lambda_est_light   = lambda_est_light_p8
    param_pft8_from_file%k_est_light   = k_est_light_p8
    param_pft8_from_file%seed_size   = seed_size_p8
    param_pft8_from_file%beta_soil_flush   = beta_soil_flush_p8
    param_pft8_from_file%beta_soil_senescence   = beta_soil_senescence_p8
    param_pft8_from_file%gdd_req_max   = gdd_req_max_p8
    param_pft8_from_file%k_gdd_dormance   = k_gdd_dormance_p8
    param_pft8_from_file%t_air_senescence   = t_air_senescence_p8
    param_pft8_from_file%min_leaf_age   = min_leaf_age_p8
    param_pft8_from_file%frac_sapwood_branch   = frac_sapwood_branch_p8
    param_pft8_from_file%wood_density   = wood_density_p8
    param_pft8_from_file%k_latosa   = k_latosa_p8
    param_pft8_from_file%k_crtos   = k_crtos_p8
    param_pft8_from_file%k_rtos   = k_rtos_p8
    param_pft8_from_file%k2_fruit_alloc   = k2_fruit_alloc_p8
    param_pft8_from_file%allom_k1   = allom_k1_p8
    param_pft8_from_file%allom_k2   = allom_k2_p8
    param_pft8_from_file%phi_leaf_min   = phi_leaf_min_p8
    param_pft8_from_file%k_root   = k_root_p8
    param_pft8_from_file%k_sapwood   = k_sapwood_p8
    param_pft8_from_file%c0_allom   = c0_allom_p8
    param_pft8_from_file%fstore_target   = fstore_target_p8
    param_pft8_from_file%k_root_dist   = k_root_dist_p8

  END SUBROUTINE read_parameter_lctlib_list


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded lctlib parameter values with values from file
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_parameter_lctlib_values_with_list_values(model)

    USE mo_jsb_model_class,                 ONLY: t_jsb_model

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_model), INTENT(inout) :: model
    ! ---------------------------
    ! 0.2 Local
    INTEGER                     :: pft_id 
    LOGICAL                     :: prop_parameter_changes
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_parameter_lctlib_values_with_list_values'



    !> 0.9 init local variables
    !! 
    prop_parameter_changes = model%config%set_parameter_values_proportional


    !> 1.0 read parameter lctlib values from file
    !! 
    CALL read_parameter_lctlib_list

    !> 2.0 replace parameter lctlib values for the particular PFT X that is simulated
    !! 
    !! lctlib(4)%sla is needed for calculation of k_rtos
    CALL change_lctlib_param_values_pft_x(model%lctlib(1), param_pft1_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    CALL change_lctlib_param_values_pft_x(model%lctlib(2), param_pft2_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    CALL change_lctlib_param_values_pft_x(model%lctlib(3), param_pft3_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    CALL change_lctlib_param_values_pft_x(model%lctlib(4), param_pft4_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    CALL change_lctlib_param_values_pft_x(model%lctlib(5), param_pft5_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    CALL change_lctlib_param_values_pft_x(model%lctlib(6), param_pft6_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    CALL change_lctlib_param_values_pft_x(model%lctlib(7), param_pft7_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    CALL change_lctlib_param_values_pft_x(model%lctlib(8), param_pft8_from_file, model%lctlib(4)%sla, prop_parameter_changes)

  END SUBROUTINE set_parameter_lctlib_values_with_list_values


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded lctlib parameter values for one PFT with values from file
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE change_lctlib_param_values_pft_x(lctlib, param_pft_from_file, pft4_sla, proportional_parameter_changes)

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_base_constants,          ONLY: one_day, one_year
    USE mo_jsb_phys_constants,      ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero 
    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, sm2lm_grass, igrass, itree

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_lctlib_element),               INTENT(inout) :: lctlib
    TYPE(t_parameters_lctlib_from_file),  INTENT(in)    :: param_pft_from_file
    REAL(wp),                             INTENT(in)    :: pft4_sla
    LOGICAL,                              INTENT(in)    :: proportional_parameter_changes
    ! ---------------------------
    ! 0.2 Local
    REAL(wp)                    :: pft4_sla_hlp
    LOGICAL                     :: dependency_recalculation   ! set TRUE for parameters that need to be re-calculated because a parameter it depends on got a new value
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':change_lctlib_param_values_pft_x'


    !> 0.9 set sla_helper to the sla value of pft4 always
    !! 
    pft4_sla_hlp = pft4_sla


    !> 1.0 replace lctlib parameter-value for the particular PFT if a new value is provided
    !! 
    CALL Set_quincy_parameter_value(1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * param_pft_from_file%sla, &
                                    param_pft_from_file%sla, &
                                    "lctlib_sla", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%sla)
    ! if the PFT of this simulation is 4 re-set the sla_helper variable (just in case the lctlib%sla was changed)
    IF (lctlib%lctnumber == 4) pft4_sla_hlp = lctlib%sla
    
    CALL Set_quincy_parameter_value(param_pft_from_file%sigma_vis, &
                                    param_pft_from_file%sigma_vis, &
                                    "lctlib_sigma_vis", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%sigma_vis)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%sigma_nir, &
                                    param_pft_from_file%sigma_nir, &
                                    "lctlib_sigma_nir", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%sigma_nir)

    CALL Set_quincy_parameter_value(param_pft_from_file%omega_clumping, &
                                    param_pft_from_file%omega_clumping, &
                                    "lctlib_omega_clumping", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%omega_clumping)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%crown_shape_factor, &
                                    param_pft_from_file%crown_shape_factor, &
                                    "lctlib_crown_shape_factor", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%crown_shape_factor)
    
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / &
                                      param_pft_from_file%cn_leaf, &
                                    param_pft_from_file%cn_leaf, &
                                    "lctlib_cn_leaf", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%cn_leaf)

    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / &
                                      param_pft_from_file%cn_leaf_min, &
                                    param_pft_from_file%cn_leaf_min, &
                                    "lctlib_cn_leaf_min", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%cn_leaf_min)
    
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / &
                                      param_pft_from_file%cn_leaf_max, &
                                    param_pft_from_file%cn_leaf_max, &
                                    "lctlib_cn_leaf_max", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%cn_leaf_max)
    
    ! not resolving the dependency to cn_leaf (as re-calc is not so easy, in case proportional changes are made, as the original input lctlib parameter is not available here)
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib%cn_leaf / &
                                      param_pft_from_file%np_leaf, &
                                    param_pft_from_file%np_leaf, &
                                    "lctlib_np_leaf", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%np_leaf)

    ! not resolving the dependency to cn_leaf (as re-calc is not so easy, in case proportional changes are made, as the original input lctlib parameter is not available here)
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib%cn_leaf / &
                                      param_pft_from_file%np_leaf_min, &
                                    param_pft_from_file%np_leaf_min, &
                                    "lctlib_np_leaf_min", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%np_leaf_min)

    ! not resolving the dependency to cn_leaf (as re-calc is not so easy, in case proportional changes are made, as the original input lctlib parameter is not available here)
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib%cn_leaf / &
                                      param_pft_from_file%np_leaf_max, &
                                    param_pft_from_file%np_leaf_max, &
                                    "lctlib_np_leaf_max", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%np_leaf_max)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%k0_fn_struc, &
                                    param_pft_from_file%k0_fn_struc, &
                                    "lctlib_k0_fn_struc", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k0_fn_struc)

    ! is always recalculated in case any of its dependencies has changed (this generality is fine as fn_oth_min is not varied itself; otherwise dependency_recalculation would need to be set correctly)
    lctlib%fn_oth_min = lctlib%k0_fn_struc - 5.0_wp / molar_mass_N * 0.5_wp * 1.e3_wp * 71.4_wp * molar_mass_N / 1.e6_wp
    ! no additional parameter changes to fn_oth_min
    ! CALL Set_quincy_parameter_value(param_pft_from_file%k0_fn_struc - 5.0_wp / molar_mass_N * 0.5_wp * 1.e3_wp * 71.4_wp * &
    !                                   molar_mass_N / 1.e6_wp, &
    !                                 param_pft_from_file%fn_oth_min, &
    !                                 "lctlib_fn_oth_min", proportional_parameter_changes, dependency_recalculation, &
    !                                 lctlib%fn_oth_min)

    CALL Set_quincy_parameter_value(param_pft_from_file%t_jmax_opt, &
                                    param_pft_from_file%t_jmax_opt, &
                                    "lctlib_t_jmax_opt", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%t_jmax_opt)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%t_jmax_omega, &
                                    param_pft_from_file%t_jmax_omega, &
                                    "lctlib_t_jmax_omega", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%t_jmax_omega)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%g0, &
                                    param_pft_from_file%g0, &
                                    "lctlib_g0", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%g0)

    CALL Set_quincy_parameter_value(param_pft_from_file%g1, &
                                    param_pft_from_file%g1, &
                                    "lctlib_g1", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%g1)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%gmin, &
                                    param_pft_from_file%gmin, &
                                    "lctlib_gmin", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%gmin)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%tau_leaf / 12.0_wp, &
                                    param_pft_from_file%tau_leaf, &
                                    "lctlib_tau_leaf", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_leaf)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_fine_root, &
                                    param_pft_from_file%tau_fine_root, &
                                    "lctlib_tau_fine_root", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_fine_root)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_coarse_root, &
                                    param_pft_from_file%tau_coarse_root, &
                                    "lctlib_tau_coarse_root", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_coarse_root)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%tau_branch, &
                                    param_pft_from_file%tau_branch, &
                                    "lctlib_tau_branch", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_branch)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%tau_sap_wood, &
                                    param_pft_from_file%tau_sap_wood, &
                                    "lctlib_tau_sap_wood", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_sap_wood)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_fruit, &
                                    param_pft_from_file%tau_fruit, &
                                    "lctlib_tau_fruit", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_fruit)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%tau_seed_litter, &
                                    param_pft_from_file%tau_seed_litter, &
                                    "lctlib_tau_seed_litter", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_seed_litter)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%tau_seed_est, &
                                    param_pft_from_file%tau_seed_est, &
                                    "lctlib_tau_seed_est", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_seed_est)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_mycorrhiza, &
                                    param_pft_from_file%tau_mycorrhiza, &
                                    "lctlib_tau_mycorrhiza", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%tau_mycorrhiza)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%vmax_uptake_n, &
                                    param_pft_from_file%vmax_uptake_n, &
                                    "lctlib_vmax_uptake_n", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%vmax_uptake_n)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%vmax_uptake_p, &
                                    param_pft_from_file%vmax_uptake_p, &
                                    "lctlib_vmax_uptake_p", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%vmax_uptake_p)

    CALL Set_quincy_parameter_value(param_pft_from_file%bnf_base / molar_mass_N / one_day / one_year * 1.e6_wp, &
                                    param_pft_from_file%bnf_base, &
                                    "lctlib_bnf_base", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%bnf_base)

    CALL Set_quincy_parameter_value(param_pft_from_file%lambda_est_light, &
                                    param_pft_from_file%lambda_est_light, &
                                    "lctlib_lambda_est_light", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%lambda_est_light)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%k_est_light, &
                                    param_pft_from_file%k_est_light, &
                                    "lctlib_k_est_light", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k_est_light)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%seed_size, &
                                    param_pft_from_file%seed_size, &
                                    "lctlib_seed_size", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%seed_size)

    CALL Set_quincy_parameter_value(param_pft_from_file%beta_soil_flush, &
                                    param_pft_from_file%beta_soil_flush, &
                                    "lctlib_beta_soil_flush", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%beta_soil_flush)
    ! beta_soil_flush may not be > 1
    IF (lctlib%beta_soil_flush > 1.0_wp) lctlib%beta_soil_flush = 1.0_wp
    
    CALL Set_quincy_parameter_value(param_pft_from_file%beta_soil_senescence, &
                                    param_pft_from_file%beta_soil_senescence, &
                                    "lctlib_beta_soil_senescence", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%beta_soil_senescence)
    ! beta_soil_senescence may not be > 1
    IF (lctlib%beta_soil_senescence > 1.0_wp) lctlib%beta_soil_senescence = 1.0_wp
    
    CALL Set_quincy_parameter_value(param_pft_from_file%gdd_req_max, &
                                    param_pft_from_file%gdd_req_max, &
                                    "lctlib_gdd_req_max", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%gdd_req_max)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_gdd_dormance, &
                                    param_pft_from_file%k_gdd_dormance, &
                                    "lctlib_k_gdd_dormance", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k_gdd_dormance)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%t_air_senescence + Tzero, &
                                    param_pft_from_file%t_air_senescence, &
                                    "lctlib_t_air_senescence", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%t_air_senescence)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%min_leaf_age, &
                                    param_pft_from_file%min_leaf_age, &
                                    "lctlib_min_leaf_age", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%min_leaf_age)

    CALL Set_quincy_parameter_value(param_pft_from_file%frac_sapwood_branch, &
                                    param_pft_from_file%frac_sapwood_branch, &
                                    "lctlib_frac_sapwood_branch", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%frac_sapwood_branch)
    ! frac_sapwood_branch may not be > 1
    IF (lctlib%frac_sapwood_branch > 1.0_wp) lctlib%frac_sapwood_branch = 1.0_wp

    CALL Set_quincy_parameter_value(param_pft_from_file%wood_density * 1.e6_wp / molar_mass_C, &
                                    param_pft_from_file%wood_density, &
                                    "lctlib_wood_density", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%wood_density)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%k_latosa, &
                                    param_pft_from_file%k_latosa, &
                                    "lctlib_k_latosa", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k_latosa)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%k_crtos, &
                                    param_pft_from_file%k_crtos, &
                                    "lctlib_k_crtos", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k_crtos)

    ! is always recalculated in case any of its dependencies has changed (this generality is fine as k_rtos is not varied itself; otherwise dependency_recalculation would need to be set correctly)
    IF(lctlib%growthform /= igrass) THEN
      lctlib%k_rtos = 1.75_wp * (pft4_sla_hlp * lctlib%wood_density) / lctlib%k_latosa
    ELSE
      lctlib%k_rtos = 0.5_wp / sm2lm_grass
    ENDIF
    ! no additional parameter changes to k_rtos
    ! CALL Set_quincy_parameter_value(k_rtos_hlp, &
    !                                 param_pft_from_file%k_rtos, &
    !                                 "lctlib_k_rtos", proportional_parameter_changes, dependency_recalculation, &
    !                                 lctlib%k_rtos)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%k2_fruit_alloc, &
                                    param_pft_from_file%k2_fruit_alloc, &
                                    "lctlib_k2_fruit_alloc", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k2_fruit_alloc)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%allom_k1, &
                                    param_pft_from_file%allom_k1, &
                                    "lctlib_allom_k1", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%allom_k1)

    CALL Set_quincy_parameter_value(param_pft_from_file%allom_k2, &
                                    param_pft_from_file%allom_k2, &
                                    "lctlib_allom_k2", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%allom_k2)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%phi_leaf_min, &
                                    param_pft_from_file%phi_leaf_min, &
                                    "lctlib_phi_leaf_min", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%phi_leaf_min)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%k_root * molar_mass_C * 1.e-10_wp, &
                                    param_pft_from_file%k_root, &
                                    "lctlib_k_root", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k_root)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_sapwood * 1.e-3_wp, &
                                    param_pft_from_file%k_sapwood, &
                                    "lctlib_k_sapwood", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k_sapwood)

    ! is always recalculated in case any of its dependencies has changed (this generality is fine as c0_allom is not varied itself; otherwise dependency_recalculation would need to be set correctly)
    IF(lctlib%growthform == itree) THEN
      lctlib%c0_allom = SQRT(lctlib%k_root / lctlib%k_sapwood * lctlib%tau_sap_wood / lctlib%tau_fine_root * &
                             lctlib%wood_density) * & 
                        lctlib%tau_fine_root / lctlib%tau_sap_wood
    ELSE
      lctlib%c0_allom = 0.0_wp
    ENDIF
    ! no additional parameter changes to c0_allom
    ! CALL Set_quincy_parameter_value(c0_allom_hlp, &
    !                                 param_pft_from_file%c0_allom, &
    !                                 "lctlib_c0_allom", proportional_parameter_changes, dependency_recalculation, &
    !                                 lctlib%c0_allom)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%fstore_target, &
                                    param_pft_from_file%fstore_target, &
                                    "lctlib_fstore_target", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%fstore_target)
    
    CALL Set_quincy_parameter_value(param_pft_from_file%k_root_dist, &
                                    param_pft_from_file%k_root_dist, &
                                    "lctlib_k_root_dist", proportional_parameter_changes, dependency_recalculation, &
                                    lctlib%k_root_dist)

  END SUBROUTINE change_lctlib_param_values_pft_x



  !-----------------------------------------------------------------------------------------------------
  !> return new QUINCY parameter value as defined by namellist \n
  !!  namelist created from generate_parameter_table.sh: parameter_slm_run.list \n
  !!  i.e., replace hard coded parameter value
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Set_quincy_parameter_value(nml_abs_val, nml_prop_chng, param_name, proportional, dependency, def_val)

    USE mo_kind,                ONLY: wp, eps8
    USE mo_jsb_io,              ONLY: missval  ! = -9.E+15_wp

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),         INTENT(in)    :: nml_abs_val        !< new value from namelist, potentially re-calculated with other parameters consistently with mo_PROC_constants
    REAL(wp),         INTENT(in)    :: nml_prop_chng      !< new value from namelist, not changed, used here for proportional change of def_val x from namelist [0 < x < 2]
    CHARACTER(LEN=*), INTENT(in)    :: param_name         !< name of the parameter for writing the new value to file
    LOGICAL,          INTENT(in)    :: proportional       !< logical: calulate new value using proportional change
    LOGICAL,          INTENT(inout) :: dependency         !< logical: this parameter value may be re-calculated because of a dependency to another param that has been changed 
    REAL(wp),         INTENT(inout) :: def_val            !< default parameter value as hardcoded in the mo_PROC_constants

    ! ---------------------------
    ! 0.2 Local
    REAL(wp)             :: param_value_diff    ! for negative def_val the absolute change-value is calculated
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_quincy_parameter_value'


    !> 1.0 re-calc def_val if namelist value differs from missval (i.e., a new namelist value is set) \n 
    !!      i.e., if a new value was defined in the namelist (parameter_slm_run.list) for this parameter \n
    !!      need to use the nml_prop_chng variable here, because the nml_abs_val might be changed because it has always "potential other parameters already included"
    IF ((ABS(nml_prop_chng - missval) > eps8) .OR. dependency) THEN
      ! differ between the two cases: just use the new value / change the (hardcoded) default value by x%
      IF (proportional) THEN
        ! special cases for some parameters, such as, temperatures given in K but modified based on the °C value
        SELECT CASE(TRIM(param_name))
        CASE('t_ref_decomposition','t_opt_nitrification','lctlib_t_air_senescence')
            ! values assumed to be in K but changed based on their °C value
            def_val = 273.15 + (def_val - 273.15) * nml_prop_chng
          CASE DEFAULT
            ! differ between negative and positive parameter values (because, e.g., 0.9 * -3.0 = -2.7 but should be -3.3)
            IF (def_val < 0.0_wp) THEN
              param_value_diff = def_val * ABS(1.0_wp - nml_prop_chng)
              IF (nml_prop_chng < 1.0) THEN
                def_val = def_val - param_value_diff
              ELSE
                def_val = def_val + param_value_diff
              ENDIF
            ELSE
              def_val = def_val * nml_prop_chng
            ENDIF
        END SELECT
      ELSE
        def_val = nml_abs_val
      ENDIF
    ENDIF
    
    !> 2.0 set parameter dependency back to the default value
    !!
    dependency = .FALSE.
    
    !> 3.0 write new value into file
    !!
    ! ...
    
  END SUBROUTINE Set_quincy_parameter_value

END MODULE mo_shared_set_parameters
