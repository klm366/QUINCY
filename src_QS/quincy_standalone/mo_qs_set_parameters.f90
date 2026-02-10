!> Contains the routines for reading parameters from parameter_slm_run.list and replace parameter-values as defined in mo_PROC_constants
!>
!> QUINCY standalone
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: https://doi.org/10.17871/quincy-model-2019
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_qs_set_parameters

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_jsb_io,              ONLY: missval  ! = -9.E+15_wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: parameters_slm_from_file                        ! it's an object
  PUBLIC :: set_parameter_lctlib_values_with_list_values    ! it's a subroutine

  !-----------------------------------------------------------------------------------------------------
  !> QUINCY process parameter values from file parameter_slm_run.list - replacing the values from mo_PROC_constants
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE t_parameters_slm_from_file
    !LOGICAL            :: active   ! de/activate this TYPE (this is non-functional for now 20161004)
    CHARACTER(LEN=30)  :: filename = 'parameter_slm_run.list'                 ! filename
    ! lnd_q_assimi_nml
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
    REAL(wp)           :: soa_b
    REAL(wp)           :: soa_t_s
    ! lnd_q_pheno_nml
    REAL(wp)           :: gdd_t_air_threshold
    REAL(wp)           :: k_leafon_canopy
    ! jsb_rad_nml
    REAL(wp)           :: kbl0_vis
    REAL(wp)           :: kdf0_vis
    REAL(wp)           :: kbl0_nir
    REAL(wp)           :: kdf0_nir
    REAL(wp)           :: rho2sbeta
    REAL(wp)           :: k_csf
    REAL(wp)           :: min_cos_zenith_angle_rad
    REAL(wp)           :: rad2ppfd
    ! qs_shared_nml
    REAL(wp)           :: mavg_period_weekly
    REAL(wp)           :: mavg_period_monthly
    REAL(wp)           :: mavg_period_yearly
    REAL(wp)           :: mavg_period_tsoa
    REAL(wp)           :: mavg_period_tuptake
    REAL(wp)           :: mavg_period_tlabile
    REAL(wp)           :: mavg_period_tphen
    REAL(wp)           :: mavg_period_tenzyme
    REAL(wp)           :: mavg_period_tresidual
    REAL(wp)           :: mavg_period_tfrac
    REAL(wp)           :: mavg_period_tcnl
    REAL(wp)           :: mavg_period_talloc
    REAL(wp)           :: mavg_period_talloc_dynamic
    REAL(wp)           :: mavg_period_tacclim
    REAL(wp)           :: mavg_period_tmic
    REAL(wp)           :: mavg_period_tgrowth
    REAL(wp)           :: mavg_period_tvegdyn
    ! lnd_spq_nml
    REAL(wp)           :: soil_heat_cap
    REAL(wp)           :: soil_therm_cond
    REAL(wp)           :: InterceptionEfficiency
    REAL(wp)           :: w_skin_max
    REAL(wp)           :: drain_frac
    REAL(wp)           :: fdrain_srf_runoff
    REAL(wp)           :: height_wind
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
    REAL(wp)           :: f_psensi_soil_sand
    REAL(wp)           :: f_psensi_soil_silt
    REAL(wp)           :: f_psensi_soil_awc_prescribe
    REAL(wp)           :: f_psensi_bulk_density
    ! lnd_veg_nml
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
    REAL(wp)           :: k_cost_symb_bnf
    REAL(wp)           :: k1_temp_bnf
    REAL(wp)           :: k2_temp_bnf
    REAL(wp)           :: k3_temp_bnf
    REAL(wp)           :: k4_temp_bnf
    REAL(wp)           :: transform_cost_nh4
    REAL(wp)           :: transform_cost_no3
    REAL(wp)           :: transform_cost_org
    REAL(wp)           :: f_nacq_nh4
    REAL(wp)           :: f_nacq_no3
    REAL(wp)           :: eta_nfixation
    REAL(wp)           :: k0_fn_chl_C3
    REAL(wp)           :: k0_fn_chl_C4
    REAL(wp)           :: k1_fn_chl_C3
    REAL(wp)           :: k1_fn_chl_C4
    REAL(wp)           :: kfn_chl
    REAL(wp)           :: min_lai
    REAL(wp)           :: target_lai_max
    REAL(wp)           :: max_leaf_shedding_rate
    REAL(wp)           :: background_mort_rate_tree
    REAL(wp)           :: background_mort_rate_grass
    REAL(wp)           :: lambda_mort_light
    REAL(wp)           :: k_mort_light
    REAL(wp)           :: lambda_est_temp
    REAL(wp)           :: k_est_temp
    REAL(wp)           :: lambda_est_moist
    REAL(wp)           :: k_est_moist
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
    REAL(wp)           :: k1_root_alloc
    REAL(wp)           :: leaf2root_min_ratio
    REAL(wp)           :: leaf2root_max_ratio
    REAL(wp)           :: delta_n_fraction
    REAL(wp)           :: delta_n_leaf
    REAL(wp)           :: f_myc_exudation_max
    REAL(wp)           :: f_mycorrhization_max
    REAL(wp)           :: rfr_ratio_toc
    REAL(wp)           :: k_r2fr_chl
    ! lnd_sb_nml
    REAL(wp)           :: tau_soluable_litter
    REAL(wp)           :: tau_metabolic_litter
    REAL(wp)           :: tau_structural_litter
    REAL(wp)           :: tau_woody_litter
    REAL(wp)           :: tau_fast
    REAL(wp)           :: tau_slow
    REAL(wp)           :: tau_biomineralisation
    REAL(wp)           :: temp_ref_tresponse
    REAL(wp)           :: ea_depolymerisation
    REAL(wp)           :: ea_mic_uptake
    REAL(wp)           :: ea_sorption
    REAL(wp)           :: ea_desorption
    REAL(wp)           :: ea_hsc
    REAL(wp)           :: ea_decomposition
    REAL(wp)           :: ed_decomposition
    REAL(wp)           :: t_ref_decomposition
    REAL(wp)           :: k1_afps
    REAL(wp)           :: k2_afps
    REAL(wp)           :: kmr_depolymerisation
    REAL(wp)           :: kmr_mic_uptake
    REAL(wp)           :: kmr_hsc
    REAL(wp)           :: kmr_sorption
    REAL(wp)           :: kmr_decomposition
    REAL(wp)           :: k_diff_org
    REAL(wp)           :: rho_bulk_org
    REAL(wp)           :: carbon_per_dryweight_SOM
    REAL(wp)           :: k_adsorpt_dom
    REAL(wp)           :: k_desorpt_dom
    REAL(wp)           :: k_adsorpt_det
    REAL(wp)           :: k_desorpt_det
    REAL(wp)           :: k_adsorpt_nh4
    REAL(wp)           :: k_desorpt_nh4
    REAL(wp)           :: k_adsorpt_po4
    REAL(wp)           :: k_desorpt_po4
    REAL(wp)           :: k_adsorpt_f_po4
    REAL(wp)           :: k_desorpt_f_po4
    REAL(wp)           :: microbial_cn
    REAL(wp)           :: microbial_np
    REAL(wp)           :: microbial_nue
    REAL(wp)           :: microbial_pue
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: microbial_cue_max
    REAL(wp)           :: microbial_cue_min
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
    REAL(wp)           :: vmax_nitrification_nh4
    REAL(wp)           :: km_nitrification_nh4
    REAL(wp)           :: f_nit_carrier_nh4
    REAL(wp)           :: vmax_denitrification_no3
    REAL(wp)           :: km_denit_no3
    REAL(wp)           :: f_denit_carrier_no3
    REAL(wp)           :: qmax_po4_mineral
    REAL(wp)           :: qmax_nh4_clay
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
    REAL(wp)           :: km_asymb_bnf
    REAL(wp)           :: vmax_asymb_bnf
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
    REAL(wp)           :: km_root_enzyme_weath_po4
    REAL(wp)           :: biomin_threshold_pc_ratio_mic_res
    REAL(wp)           :: fc_p_res2dom
    REAL(wp)           :: fc_n_res2dom
    REAL(wp)           :: f_norg_limit
    REAL(wp)           :: f_tau_myc
    REAL(wp)           :: km_priming
    REAL(wp)           :: km_up_org
    REAL(wp)           :: myc_root_coverage_max
    REAL(wp)           :: resp_priming
    REAL(wp)           :: vmax_uptake_norg
    REAL(wp)           :: k_fast_som_cn_max
    REAL(wp)           :: k_fast_som_cn_min
    REAL(wp)           :: f_fast_som_cn
    REAL(wp)           :: k_slow_som_cn
    REAL(wp)           :: k_fast_som_np
    REAL(wp)           :: k_slow_som_np
    REAL(wp)           :: f_psensi_soil_p_depth
    REAL(wp)           :: f_psensi_soil_p_labile
    REAL(wp)           :: f_psensi_soil_p_slow
    REAL(wp)           :: f_psensi_soil_p_occluded
    REAL(wp)           :: f_psensi_soil_p_primary
    REAL(wp)           :: f_psensi_qmax_org_fine_particle
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
    REAL(wp)            :: g1_medlyn
    REAL(wp)            :: g1_bberry
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
    REAL(wp)            :: vmax_symb_bnf
    REAL(wp)            :: km_symb_bnf
    ! Vegetation dynamics parameters
    REAL(wp)            :: lambda_est_light
    REAL(wp)            :: k_est_light
    REAL(wp)            :: seed_size
    REAL(wp)            :: k1_mort_greff
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

  INTEGER, PARAMETER          :: lun_param_out = 711

  CHARACTER(len=*), PARAMETER :: modname = 'mo_qs_set_parameters'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> read SLM parameter values from file
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_parameter_slm_list(this)

    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_parameters_slm_from_file), INTENT(inout) :: this
    ! ---------------------------
    ! 0.2 Local
    ! lnd_q_assimi_nml
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
    REAL(wp)           :: soa_b
    REAL(wp)           :: soa_t_s
    ! lnd_q_pheno_nml
    REAL(wp)           :: gdd_t_air_threshold
    REAL(wp)           :: k_leafon_canopy
    ! jsb_rad_nml
    REAL(wp)           :: kbl0_vis
    REAL(wp)           :: kdf0_vis
    REAL(wp)           :: kbl0_nir
    REAL(wp)           :: kdf0_nir
    REAL(wp)           :: rho2sbeta
    REAL(wp)           :: k_csf
    REAL(wp)           :: min_cos_zenith_angle_rad
    REAL(wp)           :: rad2ppfd
    ! qs_shared_nml
    REAL(wp)           :: mavg_period_weekly
    REAL(wp)           :: mavg_period_monthly
    REAL(wp)           :: mavg_period_yearly
    REAL(wp)           :: mavg_period_tsoa
    REAL(wp)           :: mavg_period_tuptake
    REAL(wp)           :: mavg_period_tlabile
    REAL(wp)           :: mavg_period_tphen
    REAL(wp)           :: mavg_period_tenzyme
    REAL(wp)           :: mavg_period_tresidual
    REAL(wp)           :: mavg_period_tfrac
    REAL(wp)           :: mavg_period_tcnl
    REAL(wp)           :: mavg_period_talloc
    REAL(wp)           :: mavg_period_talloc_dynamic
    REAL(wp)           :: mavg_period_tacclim
    REAL(wp)           :: mavg_period_tmic
    REAL(wp)           :: mavg_period_tgrowth
    REAL(wp)           :: mavg_period_tvegdyn
    ! lnd_spq_nml
    REAL(wp)           :: soil_heat_cap
    REAL(wp)           :: soil_therm_cond
    REAL(wp)           :: InterceptionEfficiency
    REAL(wp)           :: w_skin_max
    REAL(wp)           :: drain_frac
    REAL(wp)           :: fdrain_srf_runoff
    REAL(wp)           :: height_wind
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
    REAL(wp)           :: f_psensi_soil_sand
    REAL(wp)           :: f_psensi_soil_silt
    REAL(wp)           :: f_psensi_soil_awc_prescribe
    REAL(wp)           :: f_psensi_bulk_density
    ! lnd_veg_nml
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
    REAL(wp)           :: k_cost_symb_bnf
    REAL(wp)           :: k1_temp_bnf
    REAL(wp)           :: k2_temp_bnf
    REAL(wp)           :: k3_temp_bnf
    REAL(wp)           :: k4_temp_bnf
    REAL(wp)           :: transform_cost_nh4
    REAL(wp)           :: transform_cost_no3
    REAL(wp)           :: transform_cost_org
    REAL(wp)           :: f_nacq_nh4
    REAL(wp)           :: f_nacq_no3
    REAL(wp)           :: eta_nfixation
    REAL(wp)           :: k0_fn_chl_C3
    REAL(wp)           :: k0_fn_chl_C4
    REAL(wp)           :: k1_fn_chl_C3
    REAL(wp)           :: k1_fn_chl_C4
    REAL(wp)           :: kfn_chl
    REAL(wp)           :: min_lai
    REAL(wp)           :: target_lai_max
    REAL(wp)           :: max_leaf_shedding_rate
    REAL(wp)           :: background_mort_rate_tree
    REAL(wp)           :: background_mort_rate_grass
    REAL(wp)           :: lambda_mort_light
    REAL(wp)           :: k_mort_light
    REAL(wp)           :: lambda_est_temp
    REAL(wp)           :: k_est_temp
    REAL(wp)           :: lambda_est_moist
    REAL(wp)           :: k_est_moist
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
    REAL(wp)           :: k1_root_alloc
    REAL(wp)           :: leaf2root_min_ratio
    REAL(wp)           :: leaf2root_max_ratio
    REAL(wp)           :: delta_n_fraction
    REAL(wp)           :: delta_n_leaf
    REAL(wp)           :: f_myc_exudation_max
    REAL(wp)           :: f_mycorrhization_max
    REAL(wp)           :: rfr_ratio_toc
    REAL(wp)           :: k_r2fr_chl
    ! lnd_sb_nml
    REAL(wp)           :: tau_soluable_litter
    REAL(wp)           :: tau_metabolic_litter
    REAL(wp)           :: tau_structural_litter
    REAL(wp)           :: tau_woody_litter
    REAL(wp)           :: tau_fast
    REAL(wp)           :: tau_slow
    REAL(wp)           :: tau_biomineralisation
    REAL(wp)           :: temp_ref_tresponse
    REAL(wp)           :: ea_depolymerisation
    REAL(wp)           :: ea_mic_uptake
    REAL(wp)           :: ea_sorption
    REAL(wp)           :: ea_desorption
    REAL(wp)           :: ea_hsc
    REAL(wp)           :: ea_decomposition
    REAL(wp)           :: ed_decomposition
    REAL(wp)           :: t_ref_decomposition
    REAL(wp)           :: k1_afps
    REAL(wp)           :: k2_afps
    REAL(wp)           :: kmr_depolymerisation
    REAL(wp)           :: kmr_mic_uptake
    REAL(wp)           :: kmr_hsc
    REAL(wp)           :: kmr_sorption
    REAL(wp)           :: kmr_decomposition
    REAL(wp)           :: k_diff_org
    REAL(wp)           :: rho_bulk_org
    REAL(wp)           :: carbon_per_dryweight_SOM
    REAL(wp)           :: k_adsorpt_dom
    REAL(wp)           :: k_desorpt_dom
    REAL(wp)           :: k_adsorpt_det
    REAL(wp)           :: k_desorpt_det
    REAL(wp)           :: k_adsorpt_nh4
    REAL(wp)           :: k_desorpt_nh4
    REAL(wp)           :: k_adsorpt_po4
    REAL(wp)           :: k_desorpt_po4
    REAL(wp)           :: k_adsorpt_f_po4
    REAL(wp)           :: k_desorpt_f_po4
    REAL(wp)           :: microbial_cn
    REAL(wp)           :: microbial_np
    REAL(wp)           :: microbial_nue
    REAL(wp)           :: microbial_pue
    REAL(wp)           :: microbial_cue
    REAL(wp)           :: microbial_cue_max
    REAL(wp)           :: microbial_cue_min
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
    REAL(wp)           :: vmax_nitrification_nh4
    REAL(wp)           :: km_nitrification_nh4
    REAL(wp)           :: f_nit_carrier_nh4
    REAL(wp)           :: vmax_denitrification_no3
    REAL(wp)           :: km_denit_no3
    REAL(wp)           :: f_denit_carrier_no3
    REAL(wp)           :: qmax_po4_mineral
    REAL(wp)           :: qmax_nh4_clay
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
    REAL(wp)           :: km_asymb_bnf
    REAL(wp)           :: vmax_asymb_bnf
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
    REAL(wp)           :: km_root_enzyme_weath_po4
    REAL(wp)           :: biomin_threshold_pc_ratio_mic_res
    REAL(wp)           :: fc_p_res2dom
    REAL(wp)           :: fc_n_res2dom
    REAL(wp)           :: f_norg_limit
    REAL(wp)           :: f_tau_myc
    REAL(wp)           :: km_priming
    REAL(wp)           :: km_up_org
    REAL(wp)           :: myc_root_coverage_max
    REAL(wp)           :: resp_priming
    REAL(wp)           :: vmax_uptake_norg
    REAL(wp)           :: k_fast_som_cn_max
    REAL(wp)           :: k_fast_som_cn_min
    REAL(wp)           :: f_fast_som_cn
    REAL(wp)           :: k_slow_som_cn
    REAL(wp)           :: k_fast_som_np
    REAL(wp)           :: k_slow_som_np
    REAL(wp)           :: f_psensi_soil_p_depth
    REAL(wp)           :: f_psensi_soil_p_labile
    REAL(wp)           :: f_psensi_soil_p_slow
    REAL(wp)           :: f_psensi_soil_p_occluded
    REAL(wp)           :: f_psensi_soil_p_primary
    REAL(wp)           :: f_psensi_qmax_org_fine_particle

    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':read_parameter_slm_list'

    ! Namelist definition of section parameters
    !
    NAMELIST /lnd_q_assimi_nml/     &
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
        soa_b                     , &
        soa_t_s

    NAMELIST /lnd_q_pheno_nml/      &
        gdd_t_air_threshold       , &
        k_leafon_canopy

    NAMELIST /jsb_rad_nml/          &
        kbl0_vis                  , &
        kdf0_vis                  , &
        kbl0_nir                  , &
        kdf0_nir                  , &
        rho2sbeta                 , &
        k_csf                     , &
        min_cos_zenith_angle_rad  , &
        rad2ppfd

    NAMELIST /qs_shared_nml/        &
        mavg_period_weekly        , &
        mavg_period_monthly       , &
        mavg_period_yearly        , &
        mavg_period_tsoa          , &
        mavg_period_tuptake       , &
        mavg_period_tlabile       , &
        mavg_period_tphen         , &
        mavg_period_tenzyme       , &
        mavg_period_tresidual     , &
        mavg_period_tfrac         , &
        mavg_period_tcnl          , &
        mavg_period_talloc        , &
        mavg_period_talloc_dynamic, &
        mavg_period_tacclim       , &
        mavg_period_tmic          , &
        mavg_period_tgrowth       , &
        mavg_period_tvegdyn

    NAMELIST /lnd_spq_nml/          &
        soil_heat_cap             , &
        soil_therm_cond           , &
        InterceptionEfficiency    , &
        w_skin_max                , &
        drain_frac                , &
        fdrain_srf_runoff         , &
        height_wind               , &
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
        kdiff_sat_max             , &
        f_psensi_soil_sand        , &
        f_psensi_soil_silt        , &
        f_psensi_soil_awc_prescribe, &
        f_psensi_bulk_density

    NAMELIST /lnd_veg_nml/          &
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
        k_cost_symb_bnf           , &
        k1_temp_bnf               , &
        k2_temp_bnf               , &
        k3_temp_bnf               , &
        k4_temp_bnf               , &
        transform_cost_nh4        , &
        transform_cost_no3        , &
        transform_cost_org        , &
        f_nacq_nh4                , &
        f_nacq_no3                , &
        eta_nfixation             , &
        k0_fn_chl_C3              , &
        k0_fn_chl_C4              , &
        k1_fn_chl_C3              , &
        k1_fn_chl_C4              , &
        kfn_chl                   , &
        min_lai                   , &
        target_lai_max            , &
        max_leaf_shedding_rate    , &
        background_mort_rate_tree , &
        background_mort_rate_grass, &
        lambda_mort_light         , &
        k_mort_light              , &
        lambda_est_temp           , &
        k_est_temp                , &
        lambda_est_moist          , &
        k_est_moist               , &
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
        k1_root_alloc             , &
        leaf2root_min_ratio       , &
        leaf2root_max_ratio       , &
        delta_n_fraction          , &
        delta_n_leaf              , &
        f_myc_exudation_max       , &
        f_mycorrhization_max      , &
        rfr_ratio_toc             , &
        k_r2fr_chl

    NAMELIST /lnd_sb_nml/               &
        tau_soluable_litter           , &
        tau_metabolic_litter          , &
        tau_structural_litter         , &
        tau_woody_litter              , &
        tau_fast                      , &
        tau_slow                      , &
        tau_biomineralisation         , &
        temp_ref_tresponse            , &
        ea_depolymerisation           , &
        ea_mic_uptake                 , &
        ea_sorption                   , &
        ea_desorption                 , &
        ea_hsc                        , &
        ea_decomposition              , &
        ed_decomposition              , &
        t_ref_decomposition           , &
        k1_afps                       , &
        k2_afps                       , &
        kmr_depolymerisation          , &
        kmr_mic_uptake                , &
        kmr_hsc                       , &
        kmr_sorption                  , &
        kmr_decomposition             , &
        k_diff_org                    , &
        rho_bulk_org                  , &
        carbon_per_dryweight_SOM      , &
        k_adsorpt_dom                 , &
        k_desorpt_dom                 , &
        k_adsorpt_det                 , &
        k_desorpt_det                 , &
        k_adsorpt_nh4                 , &
        k_desorpt_nh4                 , &
        k_adsorpt_po4                 , &
        k_desorpt_po4                 , &
        k_adsorpt_f_po4               , &
        k_desorpt_f_po4               , &
        microbial_cn                  , &
        microbial_np                  , &
        microbial_nue                 , &
        microbial_pue                 , &
        microbial_cue                 , &
        microbial_cue_max             , &
        microbial_cue_min             , &
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
        vmax_nitrification_nh4        , &
        km_nitrification_nh4          , &
        f_nit_carrier_nh4             , &
        vmax_denitrification_no3      , &
        km_denit_no3                  , &
        f_denit_carrier_no3           , &
        qmax_po4_mineral              , &
        qmax_nh4_clay                 , &
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
        km_asymb_bnf                  , &
        vmax_asymb_bnf                , &
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
        km_root_enzyme_weath_po4      , &
        biomin_threshold_pc_ratio_mic_res , &
        fc_p_res2dom                  , &
        fc_n_res2dom                  , &
        f_norg_limit                  , &
        f_tau_myc                     , &
        km_priming                    , &
        km_up_org                     , &
        myc_root_coverage_max         , &
        resp_priming                  , &
        vmax_uptake_norg              , &
        k_fast_som_cn_max             , &
        k_fast_som_cn_min             , &
        f_fast_som_cn                 , &
        k_slow_som_cn                 , &
        k_fast_som_np                 , &
        k_slow_som_np                 , &
        f_psensi_soil_p_depth         , &
        f_psensi_soil_p_labile        , &
        f_psensi_soil_p_slow          , &
        f_psensi_soil_p_occluded      , &
        f_psensi_soil_p_primary       , &
        f_psensi_qmax_org_fine_particle


    ! variables for reading model-options from namelist
    !
    INTEGER  :: nml_handler, nml_unit, istat

    CALL message(TRIM(routine), 'Starting to read parameter list')

    ! Set defaults
    !
    ! lnd_q_assimi_nml
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
    soa_b                       = missval
    soa_t_s                     = missval
    ! lnd_q_pheno_nml
    gdd_t_air_threshold = missval
    k_leafon_canopy             = missval
    ! jsb_rad_nml
    kbl0_vis                    = missval
    kdf0_vis                    = missval
    kbl0_nir                    = missval
    kdf0_nir                    = missval
    rho2sbeta                   = missval
    k_csf                       = missval
    min_cos_zenith_angle_rad    = missval
    rad2ppfd                    = missval
    ! qs_shared_nml
    mavg_period_weekly          = missval
    mavg_period_monthly         = missval
    mavg_period_yearly          = missval
    mavg_period_tsoa            = missval
    mavg_period_tuptake         = missval
    mavg_period_tlabile         = missval
    mavg_period_tphen           = missval
    mavg_period_tenzyme         = missval
    mavg_period_tresidual       = missval
    mavg_period_tfrac           = missval
    mavg_period_tcnl            = missval
    mavg_period_talloc          = missval
    mavg_period_talloc_dynamic  = missval
    mavg_period_tacclim         = missval
    mavg_period_tmic            = missval
    mavg_period_tgrowth         = missval
    mavg_period_tvegdyn         = missval
    ! lnd_spq_nml
    soil_heat_cap               = missval
    soil_therm_cond             = missval
    InterceptionEfficiency      = missval
    w_skin_max                  = missval
    drain_frac                  = missval
    fdrain_srf_runoff           = missval
    height_wind                 = missval
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
    f_psensi_soil_sand          = missval
    f_psensi_soil_silt          = missval
    f_psensi_soil_awc_prescribe = missval
    f_psensi_bulk_density       = missval
    ! lnd_veg_nml
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
    k_cost_symb_bnf             = missval
    k1_temp_bnf                 = missval
    k2_temp_bnf                 = missval
    k3_temp_bnf                 = missval
    k4_temp_bnf                 = missval
    transform_cost_nh4          = missval
    transform_cost_no3          = missval
    transform_cost_org          = missval
    f_nacq_nh4                  = missval
    f_nacq_no3                  = missval
    eta_nfixation               = missval
    k0_fn_chl_C3                = missval
    k0_fn_chl_C4                = missval
    k1_fn_chl_C3                = missval
    k1_fn_chl_C4                = missval
    kfn_chl                     = missval
    min_lai                     = missval
    target_lai_max              = missval
    max_leaf_shedding_rate      = missval
    background_mort_rate_tree   = missval
    background_mort_rate_grass  = missval
    lambda_mort_light           = missval
    k_mort_light                = missval
    lambda_est_temp             = missval
    k_est_temp                  = missval
    lambda_est_moist            = missval
    k_est_moist                 = missval
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
    k1_root_alloc               = missval
    leaf2root_min_ratio         = missval
    leaf2root_max_ratio         = missval
    delta_n_fraction            = missval
    delta_n_leaf                = missval
    f_myc_exudation_max         = missval
    f_mycorrhization_max        = missval
    rfr_ratio_toc               = missval
    k_r2fr_chl                  = missval
    ! lnd_sb_nml
    tau_soluable_litter           = missval
    tau_metabolic_litter          = missval
    tau_structural_litter         = missval
    tau_woody_litter              = missval
    tau_fast                      = missval
    tau_slow                      = missval
    tau_biomineralisation         = missval
    temp_ref_tresponse            = missval
    ea_depolymerisation           = missval
    ea_mic_uptake                 = missval
    ea_sorption                   = missval
    ea_desorption                 = missval
    ea_hsc                        = missval
    ea_decomposition              = missval
    ed_decomposition              = missval
    t_ref_decomposition           = missval
    k1_afps                       = missval
    k2_afps                       = missval
    kmr_depolymerisation          = missval
    kmr_mic_uptake                = missval
    kmr_hsc                       = missval
    kmr_sorption                  = missval
    kmr_decomposition             = missval
    k_diff_org                    = missval
    rho_bulk_org                  = missval
    carbon_per_dryweight_SOM      = missval
    k_adsorpt_dom                 = missval
    k_desorpt_dom                 = missval
    k_adsorpt_det                 = missval
    k_desorpt_det                 = missval
    k_adsorpt_nh4                 = missval
    k_desorpt_nh4                 = missval
    k_adsorpt_po4                 = missval
    k_desorpt_po4                 = missval
    k_adsorpt_f_po4               = missval
    k_desorpt_f_po4               = missval
    microbial_cn                  = missval
    microbial_np                  = missval
    microbial_nue                 = missval
    microbial_pue                 = missval
    microbial_cue                 = missval
    microbial_cue_max             = missval
    microbial_cue_min             = missval
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
    vmax_nitrification_nh4        = missval
    km_nitrification_nh4          = missval
    f_nit_carrier_nh4             = missval
    vmax_denitrification_no3      = missval
    km_denit_no3                  = missval
    f_denit_carrier_no3           = missval
    qmax_po4_mineral              = missval
    qmax_nh4_clay                 = missval
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
    km_asymb_bnf                  = missval
    vmax_asymb_bnf                = missval
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
    km_root_enzyme_weath_po4      = missval
    biomin_threshold_pc_ratio_mic_res = missval
    fc_p_res2dom                  = missval
    fc_n_res2dom                  = missval
    f_norg_limit                  = missval
    f_tau_myc                     = missval
    km_priming                    = missval
    km_up_org                     = missval
    myc_root_coverage_max         = missval
    resp_priming                  = missval
    vmax_uptake_norg              = missval
    k_fast_som_cn_max             = missval
    k_fast_som_cn_min             = missval
    f_fast_som_cn                 = missval
    k_slow_som_cn                 = missval
    k_fast_som_np                 = missval
    k_slow_som_np                 = missval
    f_psensi_soil_p_depth         = missval
    f_psensi_soil_p_labile        = missval
    f_psensi_soil_p_slow          = missval
    f_psensi_soil_p_occluded      = missval
    f_psensi_soil_p_primary       = missval
    f_psensi_qmax_org_fine_particle = missval

    ! read the namelist
    !
    ! &lnd_q_assimi_nml
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('lnd_q_assimi_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_q_assimi_nml)
    CALL close_nml(nml_handler)
    ! &lnd_q_pheno_nml
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('lnd_q_pheno_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_q_pheno_nml)
    CALL close_nml(nml_handler)
    ! &jsb_rad_nml
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('jsb_rad_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, jsb_rad_nml)
    CALL close_nml(nml_handler)
    ! &qs_shared_nml
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('qs_shared_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, qs_shared_nml)
    CALL close_nml(nml_handler)
    ! &lnd_spq_nml
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('lnd_spq_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_spq_nml)
    CALL close_nml(nml_handler)
    ! &lnd_veg_nml
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('lnd_veg_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_veg_nml)
    CALL close_nml(nml_handler)
    ! &lnd_sb_nml
    nml_handler = open_nml(TRIM(this%filename))
    nml_unit = position_nml('lnd_sb_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lnd_sb_nml)
    CALL close_nml(nml_handler)

    ! pass values as read from file
    !
    ! lnd_q_assimi_nml
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
    this%soa_b                       = soa_b
    this%soa_t_s                     = soa_t_s
    ! lnd_q_pheno_nml
    this%gdd_t_air_threshold = gdd_t_air_threshold
    this%k_leafon_canopy             = k_leafon_canopy
    ! jsb_rad_nml
    this%kbl0_vis                    = kbl0_vis
    this%kdf0_vis                    = kdf0_vis
    this%kbl0_nir                    = kbl0_nir
    this%kdf0_nir                    = kdf0_nir
    this%rho2sbeta                   = rho2sbeta
    this%k_csf                       = k_csf
    this%min_cos_zenith_angle_rad    = min_cos_zenith_angle_rad
    this%rad2ppfd                    = rad2ppfd
    ! qs_shared_nml
    this%mavg_period_weekly          = mavg_period_weekly
    this%mavg_period_monthly         = mavg_period_monthly
    this%mavg_period_yearly          = mavg_period_yearly
    this%mavg_period_tsoa            = mavg_period_tsoa
    this%mavg_period_tuptake         = mavg_period_tuptake
    this%mavg_period_tlabile         = mavg_period_tlabile
    this%mavg_period_tphen           = mavg_period_tphen
    this%mavg_period_tenzyme         = mavg_period_tenzyme
    this%mavg_period_tresidual       = mavg_period_tresidual
    this%mavg_period_tfrac           = mavg_period_tfrac
    this%mavg_period_tcnl            = mavg_period_tcnl
    this%mavg_period_talloc          = mavg_period_talloc
    this%mavg_period_talloc_dynamic  = mavg_period_talloc_dynamic
    this%mavg_period_tacclim         = mavg_period_tacclim
    this%mavg_period_tmic            = mavg_period_tmic
    this%mavg_period_tgrowth         = mavg_period_tgrowth
    this%mavg_period_tvegdyn         = mavg_period_tvegdyn
    ! lnd_spq_nml
    this%soil_heat_cap               = soil_heat_cap
    this%soil_therm_cond             = soil_therm_cond
    this%InterceptionEfficiency      = InterceptionEfficiency
    this%w_skin_max                  = w_skin_max
    this%drain_frac                  = drain_frac
    this%fdrain_srf_runoff           = fdrain_srf_runoff
    this%height_wind                 = height_wind
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
    this%f_psensi_soil_sand          = f_psensi_soil_sand
    this%f_psensi_soil_silt          = f_psensi_soil_silt
    this%f_psensi_soil_awc_prescribe = f_psensi_soil_awc_prescribe
    this%f_psensi_bulk_density       = f_psensi_bulk_density
    ! lnd_veg_nml
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
    this%k_cost_symb_bnf             = k_cost_symb_bnf
    this%k1_temp_bnf                 = k1_temp_bnf
    this%k2_temp_bnf                 = k2_temp_bnf
    this%k3_temp_bnf                 = k3_temp_bnf
    this%k4_temp_bnf                 = k4_temp_bnf
    this%transform_cost_nh4          = transform_cost_nh4
    this%transform_cost_no3          = transform_cost_no3
    this%transform_cost_org          = transform_cost_org
    this%f_nacq_nh4                  = f_nacq_nh4
    this%f_nacq_no3                  = f_nacq_no3
    this%eta_nfixation               = eta_nfixation
    this%k0_fn_chl_C3                = k0_fn_chl_C3
    this%k0_fn_chl_C4                = k0_fn_chl_C4
    this%k1_fn_chl_C3                = k1_fn_chl_C3
    this%k1_fn_chl_C4                = k1_fn_chl_C4
    this%kfn_chl                     = kfn_chl
    this%min_lai                     = min_lai
    this%target_lai_max              = target_lai_max
    this%max_leaf_shedding_rate      = max_leaf_shedding_rate
    this%background_mort_rate_tree   = background_mort_rate_tree
    this%background_mort_rate_grass  = background_mort_rate_grass
    this%lambda_mort_light           = lambda_mort_light
    this%k_mort_light                = k_mort_light
    this%lambda_est_temp             = lambda_est_temp
    this%k_est_temp                  = k_est_temp
    this%lambda_est_moist            = lambda_est_moist
    this%k_est_moist                 = k_est_moist
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
    this%k1_root_alloc               = k1_root_alloc
    this%leaf2root_min_ratio         = leaf2root_min_ratio
    this%leaf2root_max_ratio         = leaf2root_max_ratio
    this%delta_n_fraction            = delta_n_fraction
    this%delta_n_leaf                = delta_n_leaf
    this%f_myc_exudation_max         = f_myc_exudation_max
    this%f_mycorrhization_max        = f_mycorrhization_max
    this%rfr_ratio_toc               = rfr_ratio_toc
    this%k_r2fr_chl                  = k_r2fr_chl
    ! lnd_sb_nml
    this%tau_soluable_litter           = tau_soluable_litter
    this%tau_metabolic_litter          = tau_metabolic_litter
    this%tau_structural_litter         = tau_structural_litter
    this%tau_woody_litter              = tau_woody_litter
    this%tau_fast                      = tau_fast
    this%tau_slow                      = tau_slow
    this%tau_biomineralisation         = tau_biomineralisation
    this%temp_ref_tresponse            = temp_ref_tresponse
    this%ea_depolymerisation           = ea_depolymerisation
    this%ea_mic_uptake                 = ea_mic_uptake
    this%ea_sorption                   = ea_sorption
    this%ea_desorption                 = ea_desorption
    this%ea_hsc                        = ea_hsc
    this%ea_decomposition              = ea_decomposition
    this%ed_decomposition              = ed_decomposition
    this%t_ref_decomposition           = t_ref_decomposition
    this%k1_afps                       = k1_afps
    this%k2_afps                       = k2_afps
    this%kmr_depolymerisation          = kmr_depolymerisation
    this%kmr_mic_uptake                = kmr_mic_uptake
    this%kmr_hsc                       = kmr_hsc
    this%kmr_sorption                  = kmr_sorption
    this%kmr_decomposition             = kmr_decomposition
    this%k_diff_org                    = k_diff_org
    this%rho_bulk_org                  = rho_bulk_org
    this%carbon_per_dryweight_SOM      = carbon_per_dryweight_SOM
    this%k_adsorpt_dom                 = k_adsorpt_dom
    this%k_desorpt_dom                 = k_desorpt_dom
    this%k_adsorpt_det                 = k_adsorpt_det
    this%k_desorpt_det                 = k_desorpt_det
    this%k_adsorpt_nh4                 = k_adsorpt_nh4
    this%k_desorpt_nh4                 = k_desorpt_nh4
    this%k_adsorpt_po4                 = k_adsorpt_po4
    this%k_desorpt_po4                 = k_desorpt_po4
    this%k_adsorpt_f_po4               = k_adsorpt_f_po4
    this%k_desorpt_f_po4               = k_desorpt_f_po4
    this%microbial_cn                  = microbial_cn
    this%microbial_np                  = microbial_np
    this%microbial_nue                 = microbial_nue
    this%microbial_pue                 = microbial_pue
    this%microbial_cue                 = microbial_cue
    this%microbial_cue_max             = microbial_cue_max
    this%microbial_cue_min             = microbial_cue_min
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
    this%vmax_nitrification_nh4        = vmax_nitrification_nh4
    this%km_nitrification_nh4          = km_nitrification_nh4
    this%f_nit_carrier_nh4             = f_nit_carrier_nh4
    this%vmax_denitrification_no3      = vmax_denitrification_no3
    this%km_denit_no3                  = km_denit_no3
    this%f_denit_carrier_no3           = f_denit_carrier_no3
    this%qmax_po4_mineral              = qmax_po4_mineral
    this%qmax_nh4_clay                 = qmax_nh4_clay
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
    this%km_asymb_bnf                  = km_asymb_bnf
    this%vmax_asymb_bnf                = vmax_asymb_bnf
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
    this%km_root_enzyme_weath_po4      = km_root_enzyme_weath_po4
    this%biomin_threshold_pc_ratio_mic_res = biomin_threshold_pc_ratio_mic_res
    this%fc_p_res2dom                  = fc_p_res2dom
    this%fc_n_res2dom                  = fc_n_res2dom
    this%f_norg_limit                  = f_norg_limit
    this%f_tau_myc                     = f_tau_myc
    this%km_priming                    = km_priming
    this%km_up_org                     = km_up_org
    this%myc_root_coverage_max         = myc_root_coverage_max
    this%resp_priming                  = resp_priming
    this%vmax_uptake_norg              = vmax_uptake_norg
    this%k_fast_som_cn_max             = k_fast_som_cn_max
    this%k_fast_som_cn_min             = k_fast_som_cn_min
    this%f_fast_som_cn                 = f_fast_som_cn
    this%k_slow_som_cn                 = k_slow_som_cn
    this%k_fast_som_np                 = k_fast_som_np
    this%k_slow_som_np                 = k_slow_som_np
    this%f_psensi_soil_p_depth         = f_psensi_soil_p_depth
    this%f_psensi_soil_p_labile        = f_psensi_soil_p_labile
    this%f_psensi_soil_p_slow          = f_psensi_soil_p_slow
    this%f_psensi_soil_p_occluded      = f_psensi_soil_p_occluded
    this%f_psensi_soil_p_primary       = f_psensi_soil_p_primary
    this%f_psensi_qmax_org_fine_particle = f_psensi_qmax_org_fine_particle

  END SUBROUTINE read_parameter_slm_list


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded SLM parameter values with values from file
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_parameter_slm_values_with_list_values(this, model)

    USE mo_jsb_math_constants,              ONLY: one_day, one_year, eps8
    USE mo_q_assimi_parameters
    USE mo_q_assimi_constants,              ONLY: Dwv2co2_air, Dwv2co2_turb
    USE mo_q_pheno_parameters
    USE mo_rad_parameters
    USE mo_jsb_physical_constants,          ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero
    USE mo_lnd_time_averages                ! includes the mavg_period_* parameters
    USE mo_spq_constants
    USE mo_sb_constants
    USE mo_veg_constants
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: SPQ_, SB_

    ! Use of process configurations >> needed for modifications of namelist options (which are stored in the process config)
    dsl4jsb_Use_config(SPQ_)
    dsl4jsb_Use_config(SB_)

    IMPLICIT NONE

    CLASS(t_parameters_slm_from_file), INTENT(inout) :: this
    TYPE(t_jsb_model),                 INTENT(inout) :: model

    CHARACTER(32)               :: fname                      ! filename
    LOGICAL                     :: proportional_parameter_changes
    dsl4jsb_Def_config(SPQ_)
    dsl4jsb_Def_config(SB_)
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_parameter_slm_values_with_list_values'

    proportional_parameter_changes = model%config%set_parameter_values_proportional
    dsl4jsb_Get_config(SPQ_)
    dsl4jsb_Get_config(SB_)


    !> 1.0 read parameter values from file
    !!
    CALL this%read_param


    !> 2.0 open and init output file for parameter values
    !!     file is closed in section 4.0
    fname = 'parameter_sensi_param_values.txt'
    OPEN (UNIT=lun_param_out, FILE=fname, ACTION="write", STATUS="unknown")
    WRITE(lun_param_out,'(a)') 'parName parVal parValChanged'


    !> 3.0 replace parameter values from mo_PROC_constants for parameters provided by parameter_slm_run.list
    !!
    !!
    !! routine:  SUBROUTINE Set_quincy_parameter_value(nml_abs_val, nml_prop_chng, proportional, def_val)
    !!    CALL Set_quincy_parameter_value(this%xyz, &         ! or "this%xyz */+- any other parameter"
    !!                                    this%xyz, &         ! just the parameter as defined in the namelist parameter_slm_run.list
    !!                                    "parameter name", proportional_parameter_changes, &
    !!                                    xyz)                ! the original/default parameter variable as declared in mo_PROC_constants
    !!
    !!
    !! routine: SUBROUTINE Modify_quincy_namelist_parameter(nml_param_factor, param_name, nml_param)
    !!    CALL Modify_quincy_namelist_parameter(this%f_psensi_nmlp, &       ! factor for multiplication of the namelist parameter defined in the namelist parameter_slm_run.list
    !!                                          "name of nmlp", &           ! name of the namelist parameter
    !!                                          dsl4jsb_Config(PROC_)%nmlp) ! namelist parameter from the process config it is stored in
    !!

    ! lnd_q_assimi_nml
    CALL Set_quincy_parameter_value(this%jmax2n / 4.0_wp, &
                                    this%jmax2n, &
                                    "jmax2n", proportional_parameter_changes, &
                                    jmax2n)

    CALL Set_quincy_parameter_value(this%vcmax2n, &
                                    this%vcmax2n, &
                                    "vcmax2n", proportional_parameter_changes, &
                                    vcmax2n)

    CALL Set_quincy_parameter_value(this%pepc2n, &
                                    this%pepc2n, &
                                    "pepc2n", proportional_parameter_changes, &
                                    pepc2n)

    CALL Set_quincy_parameter_value(this%chl2n / 39.8_wp, &
                                    this%chl2n, &
                                    "chl2n", proportional_parameter_changes, &
                                    chl2n)

    CALL Set_quincy_parameter_value(this%k1_fn_struc * molar_mass_N / 1000000._wp, &
                                    this%k1_fn_struc, &
                                    "k1_fn_struc", proportional_parameter_changes, &
                                    k1_fn_struc)

    CALL Set_quincy_parameter_value(this%E0kc, &
                                    this%E0kc, &
                                    "E0kc", proportional_parameter_changes, &
                                    E0kc)

    CALL Set_quincy_parameter_value(this%E1kc, &
                                    this%E1kc, &
                                    "E1kc", proportional_parameter_changes, &
                                    E1kc)

    CALL Set_quincy_parameter_value(this%E0ko, &
                                    this%E0ko, &
                                    "E0ko", proportional_parameter_changes, &
                                    E0ko)

    CALL Set_quincy_parameter_value(this%E1ko, &
                                    this%E1ko, &
                                    "E1ko", proportional_parameter_changes, &
                                    E1ko)

    CALL Set_quincy_parameter_value(this%E0pcp, &
                                    this%E0pcp, &
                                    "E0pcp", proportional_parameter_changes, &
                                    E0pcp)

    CALL Set_quincy_parameter_value(this%E1pcp, &
                                    this%E1pcp, &
                                    "E1pcp", proportional_parameter_changes, &
                                    E1pcp)

    CALL Set_quincy_parameter_value(this%E0v, &
                                    this%E0v, &
                                    "E0v", proportional_parameter_changes, &
                                    E0v)

    CALL Set_quincy_parameter_value(this%E1v, &
                                    this%E1v, &
                                    "E1v", proportional_parameter_changes, &
                                    E1v)

    CALL Set_quincy_parameter_value(this%t_jmax_offset, &
                                    this%t_jmax_offset, &
                                    "t_jmax_offset", proportional_parameter_changes, &
                                    t_jmax_offset)

    CALL Set_quincy_parameter_value(this%t_jmax_slope, &
                                    this%t_jmax_slope, &
                                    "t_jmax_slope", proportional_parameter_changes, &
                                    t_jmax_slope)

    ! is always recalculated in case any of its dependencies has changed
    t_jmax_opt_min = t_jmax_offset   ! mo_q_assimi_parameters
    CALL Set_quincy_parameter_value(t_jmax_offset, &   ! default: t_jmax_opt_min = t_jmax_offset | t_jmax_opt_min may be changed only via proportional changes
                                    this%t_jmax_opt_min, &
                                    "t_jmax_opt_min", proportional_parameter_changes, &
                                    t_jmax_opt_min)

    CALL Set_quincy_parameter_value(this%t_jmax_opt_max, &
                                    this%t_jmax_opt_max, &
                                    "t_jmax_opt_max", proportional_parameter_changes, &
                                    t_jmax_opt_max)

    CALL Set_quincy_parameter_value(this%Tref_pepc, &
                                    this%Tref_pepc, &
                                    "Tref_pepc", proportional_parameter_changes, &
                                    Tref_pepc)

    CALL Set_quincy_parameter_value(this%Tbase_pepc, &
                                    this%Tbase_pepc, &
                                    "Tbase_pepc", proportional_parameter_changes, &
                                    Tbase_pepc)

    CALL Set_quincy_parameter_value(this%alpha_i, &
                                    this%alpha_i, &
                                    "alpha_i", proportional_parameter_changes, &
                                    alpha_i)

    CALL Set_quincy_parameter_value(this%pO2, &
                                    this%pO2, &
                                    "pO2", proportional_parameter_changes, &
                                    pO2)

    CALL Set_quincy_parameter_value(this%CiCa_default_C3, &
                                    this%CiCa_default_C3, &
                                    "CiCa_default_C3", proportional_parameter_changes, &
                                    CiCa_default_C3)

    CALL Set_quincy_parameter_value(this%CiCa_default_C4, &
                                    this%CiCa_default_C4, &
                                    "CiCa_default_C4", proportional_parameter_changes, &
                                    CiCa_default_C4)

    CALL Set_quincy_parameter_value(this%ps_it_max, &
                                    this%ps_it_max, &
                                    "ps_it_max", proportional_parameter_changes, &
                                    ps_it_max)

    CALL Set_quincy_parameter_value(this%ci_max, &
                                    this%ci_max, &
                                    "ci_max", proportional_parameter_changes, &
                                    ci_max)

    ! current assumption: parameters from mo_jsb_physical_constants are not modified (this%Dwv & this%Dco2)
    CALL Set_quincy_parameter_value(Dwv / Dco2, &                       ! default: Dwv2co2_air = Dwv / Dco2 | Dwv2co2_air may be changed only via proportional changes
                                    this%Dwv2co2_air, &
                                    "Dwv2co2_air", proportional_parameter_changes, &
                                    Dwv2co2_air)

    ! current assumption: parameters from mo_jsb_physical_constants are not modified (this%Dwv & this%Dco2)
    CALL Set_quincy_parameter_value(( Dwv / Dco2  ) ** (2._wp/3._wp), & ! default: Dwv2co2_turb = ( Dwv / Dco2  ) ** (2._wp/3._wp) | Dwv2co2_turb may be changed only via proportional changes
                                    this%Dwv2co2_turb, &
                                    "Dwv2co2_turb", proportional_parameter_changes, &
                                    Dwv2co2_turb)

    CALL Set_quincy_parameter_value(this%discr_ps_a_C13, &
                                    this%discr_ps_a_C13, &
                                    "discr_ps_a_C13", proportional_parameter_changes, &
                                    discr_ps_a_C13)

    CALL Set_quincy_parameter_value(this%discr_ps_b_C13, &
                                    this%discr_ps_b_C13, &
                                    "discr_ps_b_C13", proportional_parameter_changes, &
                                    discr_ps_b_C13)

    CALL Set_quincy_parameter_value(this%discr_ps_c_C13, &
                                    this%discr_ps_c_C13, &
                                    "discr_ps_c_C13", proportional_parameter_changes, &
                                    discr_ps_c_C13)

    ! is always recalculated in case any of its dependencies has changed
    discr_ps_a_C14 = 1.97_wp * discr_ps_a_C13   ! mo_q_assimi_parameters
    CALL Set_quincy_parameter_value(this%discr_ps_a_C14 * discr_ps_a_C13, &
                                    this%discr_ps_a_C14, &
                                    "discr_ps_a_C14", proportional_parameter_changes, &
                                    discr_ps_a_C14)

    ! is always recalculated in case any of its dependencies has changed
    discr_ps_b_C14 = 1.89_wp * discr_ps_b_C13   ! mo_q_assimi_parameters
    CALL Set_quincy_parameter_value(this%discr_ps_b_C14 * discr_ps_b_C13, &
                                    this%discr_ps_b_C14, &
                                    "discr_ps_b_C14", proportional_parameter_changes, &
                                    discr_ps_b_C14)

    ! is always recalculated in case any of its dependencies has changed
    discr_ps_c_C14 = 1.89_wp * discr_ps_c_C13   ! mo_q_assimi_parameters
    CALL Set_quincy_parameter_value(this%discr_ps_c_C14 * discr_ps_c_C13, &
                                    this%discr_ps_c_C14, &
                                    "discr_ps_c_C14", proportional_parameter_changes, &
                                    discr_ps_c_C14)

    CALL Set_quincy_parameter_value(this%discr_ps_phi, &
                                    this%discr_ps_phi, &
                                    "discr_ps_phi", proportional_parameter_changes, &
                                    discr_ps_phi)

    CALL Set_quincy_parameter_value(this%ka, &
                                    this%ka, &
                                    "ka", proportional_parameter_changes, &
                                    ka)

    CALL Set_quincy_parameter_value(this%kn, &
                                    this%kn, &
                                    "kn", proportional_parameter_changes, &
                                    kn)

    CALL Set_quincy_parameter_value(this%soa_b, &
                                    this%soa_b, &
                                    "soa_b", proportional_parameter_changes, &
                                    soa_b)

    CALL Set_quincy_parameter_value(this%soa_t_s, &
                                    this%soa_t_s, &
                                    "soa_t_s", proportional_parameter_changes, &
                                    soa_t_s)

    ! lnd_q_pheno_nml
    CALL Set_quincy_parameter_value(this%gdd_t_air_threshold + Tzero, &
                                    this%gdd_t_air_threshold, &
                                    "gdd_t_air_threshold", proportional_parameter_changes, &
                                    gdd_t_air_threshold)

    CALL Set_quincy_parameter_value(this%k_leafon_canopy, &
                                    this%k_leafon_canopy, &
                                    "k_leafon_canopy", proportional_parameter_changes, &
                                    k_leafon_canopy)

    ! jsb_rad_nml
    CALL Set_quincy_parameter_value(this%kbl0_vis, &
                                    this%kbl0_vis, &
                                    "kbl0_vis", proportional_parameter_changes, &
                                    kbl0_vis)

    CALL Set_quincy_parameter_value(this%kdf0_vis, &
                                    this%kdf0_vis, &
                                    "kdf0_vis", proportional_parameter_changes, &
                                    kdf0_vis)

    CALL Set_quincy_parameter_value(this%kbl0_nir, &
                                    this%kbl0_nir, &
                                    "kbl0_nir", proportional_parameter_changes, &
                                    kbl0_nir)

    CALL Set_quincy_parameter_value(this%kdf0_nir, &
                                    this%kdf0_nir, &
                                    "kdf0_nir", proportional_parameter_changes, &
                                    kdf0_nir)

    CALL Set_quincy_parameter_value(this%rho2sbeta, &
                                    this%rho2sbeta, &
                                    "rho2sbeta", proportional_parameter_changes, &
                                    rho2sbeta)

    CALL Set_quincy_parameter_value(this%k_csf, &
                                    this%k_csf, &
                                    "k_csf", proportional_parameter_changes, &
                                    k_csf)

    CALL Set_quincy_parameter_value(this%min_cos_zenith_angle_rad, &
                                    this%min_cos_zenith_angle_rad, &
                                    "min_cos_zenith_angle_rad", proportional_parameter_changes, &
                                    min_cos_zenith_angle_rad)

    CALL Set_quincy_parameter_value(this%rad2ppfd, &
                                    this%rad2ppfd, &
                                    "rad2ppfd", proportional_parameter_changes, &
                                    rad2ppfd)

    ! qs_shared_nml
    CALL Set_quincy_parameter_value(this%mavg_period_weekly, &
                                    this%mavg_period_weekly, &
                                    "mavg_period_weekly", proportional_parameter_changes, &
                                    mavg_period_weekly)

    CALL Set_quincy_parameter_value(this%mavg_period_monthly, &
                                    this%mavg_period_monthly, &
                                    "mavg_period_monthly", proportional_parameter_changes, &
                                    mavg_period_monthly)

    CALL Set_quincy_parameter_value(this%mavg_period_yearly, &
                                    this%mavg_period_yearly, &
                                    "mavg_period_yearly", proportional_parameter_changes, &
                                    mavg_period_yearly)

    CALL Set_quincy_parameter_value(this%mavg_period_tsoa, &
                                    this%mavg_period_tsoa, &
                                    "mavg_period_tsoa", proportional_parameter_changes, &
                                    mavg_period_tsoa)

    CALL Set_quincy_parameter_value(this%mavg_period_tuptake, &
                                    this%mavg_period_tuptake, &
                                    "mavg_period_tuptake", proportional_parameter_changes, &
                                    mavg_period_tuptake)

    CALL Set_quincy_parameter_value(this%mavg_period_tlabile, &
                                    this%mavg_period_tlabile, &
                                    "mavg_period_tlabile", proportional_parameter_changes, &
                                    mavg_period_tlabile)

    CALL Set_quincy_parameter_value(this%mavg_period_tphen, &
                                    this%mavg_period_tphen, &
                                    "mavg_period_tphen", proportional_parameter_changes, &
                                    mavg_period_tphen)

    CALL Set_quincy_parameter_value(this%mavg_period_tenzyme, &
                                    this%mavg_period_tenzyme, &
                                    "mavg_period_tenzyme", proportional_parameter_changes, &
                                    mavg_period_tenzyme)

    CALL Set_quincy_parameter_value(this%mavg_period_tresidual, &
                                    this%mavg_period_tresidual, &
                                    "mavg_period_tresidual", proportional_parameter_changes, &
                                    mavg_period_tresidual)

    CALL Set_quincy_parameter_value(this%mavg_period_tfrac, &
                                    this%mavg_period_tfrac, &
                                    "mavg_period_tfrac", proportional_parameter_changes, &
                                    mavg_period_tfrac)

    CALL Set_quincy_parameter_value(this%mavg_period_tcnl, &
                                    this%mavg_period_tcnl, &
                                    "mavg_period_tcnl", proportional_parameter_changes, &
                                    mavg_period_tcnl)

    CALL Set_quincy_parameter_value(this%mavg_period_talloc, &
                                    this%mavg_period_talloc, &
                                    "mavg_period_talloc", proportional_parameter_changes, &
                                    mavg_period_talloc)

    CALL Set_quincy_parameter_value(this%mavg_period_talloc_dynamic * 365.0_wp, &  ! default value: 30.0_wp * 365.0_wp (i.e., 30 years)
                                    this%mavg_period_talloc_dynamic, &
                                    "mavg_period_talloc_dynamic", proportional_parameter_changes, &
                                    mavg_period_talloc_dynamic)

    CALL Set_quincy_parameter_value(this%mavg_period_tacclim, &
                                    this%mavg_period_tacclim, &
                                    "mavg_period_tacclim", proportional_parameter_changes, &
                                    mavg_period_tacclim)

    CALL Set_quincy_parameter_value(this%mavg_period_tmic, &
                                    this%mavg_period_tmic, &
                                    "mavg_period_tmic", proportional_parameter_changes, &
                                    mavg_period_tmic)

    CALL Set_quincy_parameter_value(this%mavg_period_tgrowth, &
                                    this%mavg_period_tgrowth, &
                                    "mavg_period_tgrowth", proportional_parameter_changes, &
                                    mavg_period_tgrowth)

    CALL Set_quincy_parameter_value(this%mavg_period_tvegdyn, &
                                    this%mavg_period_tvegdyn, &
                                    "mavg_period_tvegdyn", proportional_parameter_changes, &
                                    mavg_period_tvegdyn)

    ! lnd_spq_nml
    CALL Set_quincy_parameter_value(this%soil_heat_cap, &
                                    this%soil_heat_cap, &
                                    "soil_heat_cap", proportional_parameter_changes, &
                                    soil_heat_cap)

    CALL Set_quincy_parameter_value(this%soil_therm_cond, &
                                    this%soil_therm_cond, &
                                    "soil_therm_cond", proportional_parameter_changes, &
                                    soil_therm_cond)

    CALL Set_quincy_parameter_value(this%InterceptionEfficiency, &
                                    this%InterceptionEfficiency, &
                                    "InterceptionEfficiency", proportional_parameter_changes, &
                                    InterceptionEfficiency)

    CALL Set_quincy_parameter_value(this%w_skin_max, &
                                    this%w_skin_max, &
                                    "w_skin_max", proportional_parameter_changes, &
                                    w_skin_max)

    CALL Set_quincy_parameter_value(this%drain_frac, &
                                    this%drain_frac, &
                                    "drain_frac", proportional_parameter_changes, &
                                    drain_frac)

    CALL Set_quincy_parameter_value(this%fdrain_srf_runoff, &
                                    this%fdrain_srf_runoff, &
                                    "fdrain_srf_runoff", proportional_parameter_changes, &
                                    fdrain_srf_runoff)

    CALL Set_quincy_parameter_value(this%height_wind, &
                                    this%height_wind, &
                                    "height_wind", proportional_parameter_changes, &
                                    height_wind)

    CALL Set_quincy_parameter_value(this%wpot_pwp, &
                                    this%wpot_pwp, &
                                    "wpot_pwp", proportional_parameter_changes, &
                                    wpot_pwp)

    CALL Set_quincy_parameter_value(this%wpot_fc, &
                                    this%wpot_fc, &
                                    "wpot_fc", proportional_parameter_changes, &
                                    wpot_fc)

    CALL Set_quincy_parameter_value(this%wpot_min, &
                                    this%wpot_min, &
                                    "wpot_min", proportional_parameter_changes, &
                                    wpot_min)

    CALL Set_quincy_parameter_value(this%k_pwp_a, &
                                    this%k_pwp_a, &
                                    "k_pwp_a", proportional_parameter_changes, &
                                    k_pwp_a)

    CALL Set_quincy_parameter_value(this%k_pwp_c, &
                                    this%k_pwp_c, &
                                    "k_pwp_c", proportional_parameter_changes, &
                                    k_pwp_c)

    CALL Set_quincy_parameter_value(this%k_pwp_s, &
                                    this%k_pwp_s, &
                                    "k_pwp_s", proportional_parameter_changes, &
                                    k_pwp_s)

    CALL Set_quincy_parameter_value(this%k_pwp_sc, &
                                    this%k_pwp_sc, &
                                    "k_pwp_sc", proportional_parameter_changes, &
                                    k_pwp_sc)

    CALL Set_quincy_parameter_value(this%k_pwp_at, &
                                    this%k_pwp_at, &
                                    "k_pwp_at", proportional_parameter_changes, &
                                    k_pwp_at)

    CALL Set_quincy_parameter_value(this%k_pwp_bt, &
                                    this%k_pwp_bt, &
                                    "k_pwp_bt", proportional_parameter_changes, &
                                    k_pwp_bt)

    CALL Set_quincy_parameter_value(this%k_fc_a, &
                                    this%k_fc_a, &
                                    "k_fc_a", proportional_parameter_changes, &
                                    k_fc_a)

    CALL Set_quincy_parameter_value(this%k_fc_c, &
                                    this%k_fc_c, &
                                    "k_fc_c", proportional_parameter_changes, &
                                    k_fc_c)

    CALL Set_quincy_parameter_value(this%k_fc_s, &
                                    this%k_fc_s, &
                                    "k_fc_s", proportional_parameter_changes, &
                                    k_fc_s)

    CALL Set_quincy_parameter_value(this%k_fc_sc, &
                                    this%k_fc_sc, &
                                    "k_fc_sc", proportional_parameter_changes, &
                                    k_fc_sc)

    CALL Set_quincy_parameter_value(this%k_fc_at, &
                                    this%k_fc_at, &
                                    "k_fc_at", proportional_parameter_changes, &
                                    k_fc_at)

    CALL Set_quincy_parameter_value(this%k_fc_bt, &
                                    this%k_fc_bt, &
                                    "k_fc_bt", proportional_parameter_changes, &
                                    k_fc_bt)

    CALL Set_quincy_parameter_value(this%k_fc_ct, &
                                    this%k_fc_ct, &
                                    "k_fc_ct", proportional_parameter_changes, &
                                    k_fc_ct)

    CALL Set_quincy_parameter_value(this%k_sat_a, &
                                    this%k_sat_a, &
                                    "k_sat_a", proportional_parameter_changes, &
                                    k_sat_a)

    CALL Set_quincy_parameter_value(this%k_sat_c, &
                                    this%k_sat_c, &
                                    "k_sat_c", proportional_parameter_changes, &
                                    k_sat_c)

    CALL Set_quincy_parameter_value(this%k_sat_s, &
                                    this%k_sat_s, &
                                    "k_sat_s", proportional_parameter_changes, &
                                    k_sat_s)

    CALL Set_quincy_parameter_value(this%k_sat_sc, &
                                    this%k_sat_sc, &
                                    "k_sat_sc", proportional_parameter_changes, &
                                    k_sat_sc)

    CALL Set_quincy_parameter_value(this%k_sat_at, &
                                    this%k_sat_at, &
                                    "k_sat_at", proportional_parameter_changes, &
                                    k_sat_at)

    CALL Set_quincy_parameter_value(this%k_sat_bt, &
                                    this%k_sat_bt, &
                                    "k_sat_bt", proportional_parameter_changes, &
                                    k_sat_bt)

    CALL Set_quincy_parameter_value(this%k_sat_ct, &
                                    this%k_sat_ct, &
                                    "k_sat_ct", proportional_parameter_changes, &
                                    k_sat_ct)

    CALL Set_quincy_parameter_value(this%k_sat_dt, &
                                    this%k_sat_dt, &
                                    "k_sat_dt", proportional_parameter_changes, &
                                    k_sat_dt)

    CALL Set_quincy_parameter_value(this%kdiff_sat_max  / 1000.0_wp / 3600._wp, &
                                    this%kdiff_sat_max, &
                                    "kdiff_sat_max", proportional_parameter_changes, &
                                    kdiff_sat_max)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_sand, &
                                          "soil_sand", &
                                          dsl4jsb_Config(SPQ_)%soil_sand)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_silt, &
                                          "soil_silt", &
                                          dsl4jsb_Config(SPQ_)%soil_silt)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_awc_prescribe, &
                                          "soil_awc_prescribe", &
                                          dsl4jsb_Config(SPQ_)%soil_awc_prescribe)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_bulk_density, &
                                          "bulk_density", &
                                          dsl4jsb_Config(SPQ_)%bulk_density)

    ! lnd_veg_nml
    CALL Set_quincy_parameter_value(this%fresp_growth, &
                                    this%fresp_growth, &
                                    "fresp_growth", proportional_parameter_changes, &
                                    fresp_growth)

    CALL Set_quincy_parameter_value(this%fmaint_rate_base, &
                                    this%fmaint_rate_base, &
                                    "fmaint_rate_base", proportional_parameter_changes, &
                                    fmaint_rate_base)

    ! is always recalculated in case any of its dependencies has changed
    fmaint_rate_wood = 0.25_wp * fmaint_rate_base   ! copied from mo_veg_constants
    CALL Set_quincy_parameter_value(this%fmaint_rate_wood * fmaint_rate_base, &
                                    this%fmaint_rate_wood, &
                                    "fmaint_rate_wood", proportional_parameter_changes, &
                                    fmaint_rate_wood)

    CALL Set_quincy_parameter_value(this%t_maint_k1, &
                                    this%t_maint_k1, &
                                    "t_maint_k1", proportional_parameter_changes, &
                                    t_maint_k1)

    CALL Set_quincy_parameter_value(this%t_maint_k2, &
                                    this%t_maint_k2, &
                                    "t_maint_k2", proportional_parameter_changes, &
                                    t_maint_k2)

    CALL Set_quincy_parameter_value(this%t_maint_k3, &
                                    this%t_maint_k3, &
                                    "t_maint_k3", proportional_parameter_changes, &
                                    t_maint_k3)

    CALL Set_quincy_parameter_value(this%f_resp_acclim, &
                                    this%f_resp_acclim, &
                                    "f_resp_acclim", proportional_parameter_changes, &
                                    f_resp_acclim)

    CALL Set_quincy_parameter_value(this%t_acclim_zero, &
                                    this%t_acclim_zero, &
                                    "t_acclim_zero", proportional_parameter_changes, &
                                    t_acclim_zero)

    CALL Set_quincy_parameter_value(this%k0_fn_pepc_C3, &
                                    this%k0_fn_pepc_C3, &
                                    "k0_fn_pepc_C3", proportional_parameter_changes, &
                                    k0_fn_pepc_C3)

    CALL Set_quincy_parameter_value(this%k0_fn_pepc_C4, &
                                    this%k0_fn_pepc_C4, &
                                    "k0_fn_pepc_C4", proportional_parameter_changes, &
                                    k0_fn_pepc_C4)

    CALL Set_quincy_parameter_value(this%jmax2vcmax_C3 / 4.0_wp, &
                                    this%jmax2vcmax_C3, &
                                    "jmax2vcmax_C3", proportional_parameter_changes, &
                                    jmax2vcmax_C3)

    CALL Set_quincy_parameter_value(this%jmax2vcmax_C4, &
                                    this%jmax2vcmax_C4, &
                                    "jmax2vcmax_C4", proportional_parameter_changes, &
                                    jmax2vcmax_C4)

    CALL Set_quincy_parameter_value(this%tau_labile, &
                                    this%tau_labile, &
                                    "tau_labile", proportional_parameter_changes, &
                                    tau_labile)

    CALL Set_quincy_parameter_value(this%lambda_labile_temp, &
                                    this%lambda_labile_temp, &
                                    "lambda_labile_temp", proportional_parameter_changes, &
                                    lambda_labile_temp)

    CALL Set_quincy_parameter_value(this%k_labile_temp, &
                                    this%k_labile_temp, &
                                    "k_labile_temp", proportional_parameter_changes, &
                                    k_labile_temp)

    CALL Set_quincy_parameter_value(this%lambda_labile_theta, &
                                    this%lambda_labile_theta, &
                                    "lambda_labile_theta", proportional_parameter_changes, &
                                    lambda_labile_theta)

    CALL Set_quincy_parameter_value(this%k_labile_theta, &
                                    this%k_labile_theta, &
                                    "k_labile_theta", proportional_parameter_changes, &
                                    k_labile_theta)

    CALL Set_quincy_parameter_value(this%lambda_phiphen, &
                                    this%lambda_phiphen, &
                                    "lambda_phiphen", proportional_parameter_changes, &
                                    lambda_phiphen)

    CALL Set_quincy_parameter_value(this%k_phiphen, &
                                    this%k_phiphen, &
                                    "k_phiphen", proportional_parameter_changes, &
                                    k_phiphen)

    CALL Set_quincy_parameter_value(this%lambda_phimaint, &
                                    this%lambda_phimaint, &
                                    "lambda_phimaint", proportional_parameter_changes, &
                                    lambda_phimaint)

    CALL Set_quincy_parameter_value(this%k_phimaint, &
                                    this%k_phimaint, &
                                    "k_phimaint", proportional_parameter_changes, &
                                    k_phimaint)

    CALL Set_quincy_parameter_value(this%lambda_phimaint_nut, &
                                    this%lambda_phimaint_nut, &
                                    "lambda_phimaint_nut", proportional_parameter_changes, &
                                    lambda_phimaint_nut)

    CALL Set_quincy_parameter_value(this%k_phimaint_nut, &
                                    this%k_phimaint_nut, &
                                    "k_phimaint_nut", proportional_parameter_changes, &
                                    k_phimaint_nut)

    CALL Set_quincy_parameter_value(this%lambda_phistore, &
                                    this%lambda_phistore, &
                                    "lambda_phistore", proportional_parameter_changes, &
                                    lambda_phistore)

    CALL Set_quincy_parameter_value(this%k_phistore, &
                                    this%k_phistore, &
                                    "k_phistore", proportional_parameter_changes, &
                                    k_phistore)

    CALL Set_quincy_parameter_value(this%knut_labile, &
                                    this%knut_labile, &
                                    "knut_labile", proportional_parameter_changes, &
                                    knut_labile)

    CALL Set_quincy_parameter_value(this%lambda_sinklim_ps, &
                                    this%lambda_sinklim_ps, &
                                    "lambda_sinklim_ps", proportional_parameter_changes, &
                                    lambda_sinklim_ps)

    CALL Set_quincy_parameter_value(this%k_sinklim_ps, &
                                    this%k_sinklim_ps, &
                                    "k_sinklim_ps", proportional_parameter_changes, &
                                    k_sinklim_ps)

    CALL Set_quincy_parameter_value(this%k_cnp_sinklim_ps, &
                                    this%k_cnp_sinklim_ps, &
                                    "k_cnp_sinklim_ps", proportional_parameter_changes, &
                                    k_cnp_sinklim_ps)

    CALL Set_quincy_parameter_value(this%beta_sinklim_ps_min, &
                                    this%beta_sinklim_ps_min, &
                                    "beta_sinklim_ps_min", proportional_parameter_changes, &
                                    beta_sinklim_ps_min)

    CALL Set_quincy_parameter_value(this%fstore_leaf_max, &
                                    this%fstore_leaf_max, &
                                    "fstore_leaf_max", proportional_parameter_changes, &
                                    fstore_leaf_max)

    CALL Set_quincy_parameter_value(this%fstore_fine_root_max, &
                                    this%fstore_fine_root_max, &
                                    "fstore_fine_root_max", proportional_parameter_changes, &
                                    fstore_fine_root_max)

    CALL Set_quincy_parameter_value(this%fstore_sap_wood_max, &
                                    this%fstore_sap_wood_max, &
                                    "fstore_sap_wood_max", proportional_parameter_changes, &
                                    fstore_sap_wood_max)

    CALL Set_quincy_parameter_value(this%k_phi_interact, &
                                    this%k_phi_interact, &
                                    "k_phi_interact", proportional_parameter_changes, &
                                    k_phi_interact)

    CALL Set_quincy_parameter_value(this%resorp_fract_leaf, &
                                    this%resorp_fract_leaf, &
                                    "resorp_fract_leaf", proportional_parameter_changes, &
                                    resorp_fract_leaf)

    CALL Set_quincy_parameter_value(this%resorp_fract_wood, &
                                    this%resorp_fract_wood, &
                                    "resorp_fract_wood", proportional_parameter_changes, &
                                    resorp_fract_wood)

    CALL Set_quincy_parameter_value(this%root2leaf_cn, &
                                    this%root2leaf_cn, &
                                    "root2leaf_cn", proportional_parameter_changes, &
                                    root2leaf_cn)

    CALL Set_quincy_parameter_value(this%wood2leaf_cn, &
                                    this%wood2leaf_cn, &
                                    "wood2leaf_cn", proportional_parameter_changes, &
                                    wood2leaf_cn)

    CALL Set_quincy_parameter_value(this%root2leaf_np, &
                                    this%root2leaf_np, &
                                    "root2leaf_np", proportional_parameter_changes, &
                                    root2leaf_np)

    CALL Set_quincy_parameter_value(this%wood2leaf_np, &
                                    this%wood2leaf_np, &
                                    "wood2leaf_np", proportional_parameter_changes, &
                                    wood2leaf_np)

    CALL Set_quincy_parameter_value(this%tau_nutrient_recycling, &
                                    this%tau_nutrient_recycling, &
                                    "tau_nutrient_recycling", proportional_parameter_changes, &
                                    tau_nutrient_recycling)

    CALL Set_quincy_parameter_value(this%kappa_f_demand_nc, &
                                    this%kappa_f_demand_nc, &
                                    "kappa_f_demand_nc", proportional_parameter_changes, &
                                    kappa_f_demand_nc)

    CALL Set_quincy_parameter_value(this%kappa_f_demand_pn, &
                                    this%kappa_f_demand_pn, &
                                    "kappa_f_demand_pn", proportional_parameter_changes, &
                                    kappa_f_demand_pn)

    CALL Set_quincy_parameter_value(this%k_f_demand, &
                                    this%k_f_demand, &
                                    "k_f_demand", proportional_parameter_changes, &
                                    k_f_demand)

    CALL Set_quincy_parameter_value(this%k_cost_symb_bnf / (molar_mass_C / molar_mass_N), &
                                    this%k_cost_symb_bnf, &
                                    "k_cost_symb_bnf", proportional_parameter_changes, &
                                    k_cost_symb_bnf)

    CALL Set_quincy_parameter_value(this%k1_temp_bnf, &
                                    this%k1_temp_bnf, &
                                    "k1_temp_bnf", proportional_parameter_changes, &
                                    k1_temp_bnf)

    CALL Set_quincy_parameter_value(this%k2_temp_bnf, &
                                    this%k2_temp_bnf, &
                                    "k2_temp_bnf", proportional_parameter_changes, &
                                    k2_temp_bnf)

    CALL Set_quincy_parameter_value(this%k3_temp_bnf, &
                                    this%k3_temp_bnf, &
                                    "k3_temp_bnf", proportional_parameter_changes, &
                                    k3_temp_bnf)

    CALL Set_quincy_parameter_value(this%k4_temp_bnf, &
                                    this%k4_temp_bnf, &
                                    "k4_temp_bnf", proportional_parameter_changes, &
                                    k4_temp_bnf)

    CALL Set_quincy_parameter_value(this%transform_cost_nh4 / (molar_mass_C / molar_mass_N), &
                                    this%transform_cost_nh4, &
                                    "transform_cost_nh4", proportional_parameter_changes, &
                                    transform_cost_nh4)

    CALL Set_quincy_parameter_value(this%transform_cost_no3 / (molar_mass_C / molar_mass_N), &
                                    this%transform_cost_no3, &
                                    "transform_cost_no3", proportional_parameter_changes, &
                                    transform_cost_no3)

    CALL Set_quincy_parameter_value(this%transform_cost_org / (molar_mass_C / molar_mass_N), &
                                    this%transform_cost_org, &
                                    "transform_cost_org", proportional_parameter_changes, &
                                    transform_cost_org)

    ! is always recalculated in case any of its dependencies has changed
    f_nacq_nh4 = transform_cost_no3 / (transform_cost_no3 + transform_cost_nh4)   ! copied from mo_veg_constants
    CALL Set_quincy_parameter_value(transform_cost_no3 / (transform_cost_no3 + transform_cost_nh4), &   ! only proportional changes may have an effect
                                    this%f_nacq_nh4, &
                                    "f_nacq_nh4", proportional_parameter_changes, &
                                    f_nacq_nh4)

    ! is always recalculated in case any of its dependencies has changed
    f_nacq_no3 = transform_cost_nh4 / (transform_cost_no3 + transform_cost_nh4)   ! copied from mo_veg_constants
    CALL Set_quincy_parameter_value(transform_cost_nh4 / (transform_cost_no3 + transform_cost_nh4), &   ! only proportional changes may have an effect
                                    this%f_nacq_no3, &
                                    "f_nacq_no3", proportional_parameter_changes, &
                                    f_nacq_no3)

    CALL Set_quincy_parameter_value(this%eta_nfixation, &
                                    this%eta_nfixation, &
                                    "eta_nfixation", proportional_parameter_changes, &
                                    eta_nfixation)

    CALL Set_quincy_parameter_value(this%k0_fn_chl_C3, &
                                    this%k0_fn_chl_C3, &
                                    "k0_fn_chl_C3", proportional_parameter_changes, &
                                    k0_fn_chl_C3)

    CALL Set_quincy_parameter_value(this%k0_fn_chl_C4, &
                                    this%k0_fn_chl_C4, &
                                    "k0_fn_chl_C4", proportional_parameter_changes, &
                                    k0_fn_chl_C4)

    CALL Set_quincy_parameter_value(this%k1_fn_chl_C3, &
                                    this%k1_fn_chl_C3, &
                                    "k1_fn_chl_C3", proportional_parameter_changes, &
                                    k1_fn_chl_C3)

    CALL Set_quincy_parameter_value(this%k1_fn_chl_C4, &
                                    this%k1_fn_chl_C4, &
                                    "k1_fn_chl_C4", proportional_parameter_changes, &
                                    k1_fn_chl_C4)

    CALL Set_quincy_parameter_value(this%kfn_chl, &
                                    this%kfn_chl, &
                                    "kfn_chl", proportional_parameter_changes, &
                                    kfn_chl)

    CALL Set_quincy_parameter_value(this%min_lai, &
                                    this%min_lai, &
                                    "min_lai", proportional_parameter_changes, &
                                    min_lai)

    CALL Set_quincy_parameter_value(this%target_lai_max, &
                                    this%target_lai_max, &
                                    "target_lai_max", proportional_parameter_changes, &
                                    target_lai_max)

    CALL Set_quincy_parameter_value(this%max_leaf_shedding_rate, &
                                    this%max_leaf_shedding_rate, &
                                    "max_leaf_shedding_rate", proportional_parameter_changes, &
                                    max_leaf_shedding_rate)

    CALL Set_quincy_parameter_value(this%background_mort_rate_tree, &
                                    this%background_mort_rate_tree, &
                                    "background_mort_rate_tree", proportional_parameter_changes, &
                                    background_mort_rate_tree)

    CALL Set_quincy_parameter_value(this%background_mort_rate_grass, &
                                    this%background_mort_rate_grass, &
                                    "background_mort_rate_grass", proportional_parameter_changes, &
                                    background_mort_rate_grass)

    CALL Set_quincy_parameter_value(this%lambda_mort_light, &
                                    this%lambda_mort_light, &
                                    "lambda_mort_light", proportional_parameter_changes, &
                                    lambda_mort_light)

    CALL Set_quincy_parameter_value(this%k_mort_light, &
                                    this%k_mort_light, &
                                    "k_mort_light", proportional_parameter_changes, &
                                    k_mort_light)

    CALL Set_quincy_parameter_value(this%lambda_est_temp, &
                                    this%lambda_est_temp, &
                                    "lambda_est_temp", proportional_parameter_changes, &
                                    lambda_est_temp)

    CALL Set_quincy_parameter_value(this%k_est_temp, &
                                    this%k_est_temp, &
                                    "k_est_temp", proportional_parameter_changes, &
                                    k_est_temp)

    CALL Set_quincy_parameter_value(this%lambda_est_moist, &
                                    this%lambda_est_moist, &
                                    "lambda_est_moist", proportional_parameter_changes, &
                                    lambda_est_moist)

    CALL Set_quincy_parameter_value(this%k_est_moist, &
                                    this%k_est_moist, &
                                    "k_est_moist", proportional_parameter_changes, &
                                    k_est_moist)

    CALL Set_quincy_parameter_value(this%k2_mort_greff, &
                                    this%k2_mort_greff, &
                                    "k2_mort_greff", proportional_parameter_changes, &
                                    k2_mort_greff)

    CALL Set_quincy_parameter_value(this%k3_mort_greff, &
                                    this%k3_mort_greff, &
                                    "k3_mort_greff", proportional_parameter_changes, &
                                    k3_mort_greff)

    CALL Set_quincy_parameter_value(this%min_greff, &
                                    this%min_greff, &
                                    "min_greff", proportional_parameter_changes, &
                                    min_greff)

    CALL Set_quincy_parameter_value(this%max_crown_area, &
                                    this%max_crown_area, &
                                    "max_crown_area", proportional_parameter_changes, &
                                    max_crown_area)

    CALL Set_quincy_parameter_value(this%min_diameter, &
                                    this%min_diameter, &
                                    "min_diameter", proportional_parameter_changes, &
                                    min_diameter)

    CALL Set_quincy_parameter_value(this%k_crown_area, &
                                    this%k_crown_area, &
                                    "k_crown_area", proportional_parameter_changes, &
                                    k_crown_area)

    CALL Set_quincy_parameter_value(this%k_rp, &
                                    this%k_rp, &
                                    "k_rp", proportional_parameter_changes, &
                                    k_rp)

    CALL Set_quincy_parameter_value(this%k_fpc, &
                                    this%k_fpc, &
                                    "k_fpc", proportional_parameter_changes, &
                                    k_fpc)

    CALL Set_quincy_parameter_value(this%fpc_max, &
                                    this%fpc_max, &
                                    "fpc_max", proportional_parameter_changes, &
                                    fpc_max)

    CALL Set_quincy_parameter_value(this%carbon_per_dryweight_leaf, &
                                    this%carbon_per_dryweight_leaf, &
                                    "carbon_per_dryweight_leaf", proportional_parameter_changes, &
                                    carbon_per_dryweight_leaf)

    CALL Set_quincy_parameter_value(this%min_height, &
                                    this%min_height, &
                                    "min_height", proportional_parameter_changes, &
                                    min_height)

    CALL Set_quincy_parameter_value(this%sm2lm_grass, &
                                    this%sm2lm_grass, &
                                    "sm2lm_grass", proportional_parameter_changes, &
                                    sm2lm_grass)

    CALL Set_quincy_parameter_value(this%lambda_fruit_alloc, &
                                    this%lambda_fruit_alloc, &
                                    "lambda_fruit_alloc", proportional_parameter_changes, &
                                    lambda_fruit_alloc)

    CALL Set_quincy_parameter_value(this%k1_fruit_alloc, &
                                    this%k1_fruit_alloc, &
                                    "k1_fruit_alloc", proportional_parameter_changes, &
                                    k1_fruit_alloc)

    CALL Set_quincy_parameter_value(this%k3_fruit_alloc, &
                                    this%k3_fruit_alloc, &
                                    "k3_fruit_alloc", proportional_parameter_changes, &
                                    k3_fruit_alloc)

    CALL Set_quincy_parameter_value(this%k4_fruit_alloc, &
                                    this%k4_fruit_alloc, &
                                    "k4_fruit_alloc", proportional_parameter_changes, &
                                    k4_fruit_alloc)

    CALL Set_quincy_parameter_value(this%k1_root_alloc, &
                                    this%k1_root_alloc, &
                                    "k1_root_alloc", proportional_parameter_changes, &
                                    k1_root_alloc)

    CALL Set_quincy_parameter_value(this%leaf2root_min_ratio, &
                                    this%leaf2root_min_ratio, &
                                    "leaf2root_min_ratio", proportional_parameter_changes, &
                                    leaf2root_min_ratio)

    CALL Set_quincy_parameter_value(this%leaf2root_max_ratio, &
                                    this%leaf2root_max_ratio, &
                                    "leaf2root_max_ratio", proportional_parameter_changes, &
                                    leaf2root_max_ratio)

    CALL Set_quincy_parameter_value(this%delta_n_fraction, &
                                    this%delta_n_fraction, &
                                    "delta_n_fraction", proportional_parameter_changes, &
                                    delta_n_fraction)

    CALL Set_quincy_parameter_value(this%delta_n_leaf, &
                                    this%delta_n_leaf, &
                                    "delta_n_leaf", proportional_parameter_changes, &
                                    delta_n_leaf)

    CALL Set_quincy_parameter_value(this%f_myc_exudation_max, &
                                    this%f_myc_exudation_max, &
                                    "f_myc_exudation_max", proportional_parameter_changes, &
                                    f_myc_exudation_max)

    CALL Set_quincy_parameter_value(this%f_mycorrhization_max, &
                                    this%f_mycorrhization_max, &
                                    "f_mycorrhization_max", proportional_parameter_changes, &
                                    f_mycorrhization_max)

    CALL Set_quincy_parameter_value(this%rfr_ratio_toc, &
                                    this%rfr_ratio_toc, &
                                    "rfr_ratio_toc", proportional_parameter_changes, &
                                    rfr_ratio_toc)

    CALL Set_quincy_parameter_value(this%k_r2fr_chl - 0.012_wp, &
                                    this%k_r2fr_chl, &
                                    "k_r2fr_chl", proportional_parameter_changes, &
                                    k_r2fr_chl)

    ! lnd_sb_nml
    CALL Set_quincy_parameter_value(this%tau_soluable_litter * one_day, &
                                    this%tau_soluable_litter, &
                                    "tau_soluable_litter", proportional_parameter_changes, &
                                    tau_soluable_litter)

    CALL Set_quincy_parameter_value(this%tau_metabolic_litter * one_year * one_day, &
                                    this%tau_metabolic_litter, &
                                    "tau_metabolic_litter", proportional_parameter_changes, &
                                    tau_metabolic_litter)

    CALL Set_quincy_parameter_value(this%tau_structural_litter * one_year * one_day, &
                                    this%tau_structural_litter, &
                                    "tau_structural_litter", proportional_parameter_changes, &
                                    tau_structural_litter)

    CALL Set_quincy_parameter_value(this%tau_woody_litter * one_year * one_day, &
                                    this%tau_woody_litter, &
                                    "tau_woody_litter", proportional_parameter_changes, &
                                    tau_woody_litter)

    CALL Set_quincy_parameter_value(this%tau_fast * one_year * one_day, &
                                    this%tau_fast, &
                                    "tau_fast", proportional_parameter_changes, &
                                    tau_fast)

    CALL Set_quincy_parameter_value(this%tau_slow * one_year * one_day, &
                                    this%tau_slow, &
                                    "tau_slow", proportional_parameter_changes, &
                                    tau_slow)

    CALL Set_quincy_parameter_value(this%tau_biomineralisation * one_year * one_day, &
                                    this%tau_biomineralisation, &
                                    "tau_biomineralisation", proportional_parameter_changes, &
                                    tau_biomineralisation)

    CALL Set_quincy_parameter_value(this%temp_ref_tresponse, &
                                    this%temp_ref_tresponse, &
                                    "temp_ref_tresponse", proportional_parameter_changes, &
                                    temp_ref_tresponse)

    CALL Set_quincy_parameter_value(this%ea_depolymerisation, &
                                    this%ea_depolymerisation, &
                                    "ea_depolymerisation", proportional_parameter_changes, &
                                    ea_depolymerisation)

    CALL Set_quincy_parameter_value(this%ea_mic_uptake, &
                                    this%ea_mic_uptake, &
                                    "ea_mic_uptake", proportional_parameter_changes, &
                                    ea_mic_uptake)

    CALL Set_quincy_parameter_value(this%ea_sorption, &
                                    this%ea_sorption, &
                                    "ea_sorption", proportional_parameter_changes, &
                                    ea_sorption)

    CALL Set_quincy_parameter_value(this%ea_desorption, &
                                    this%ea_desorption, &
                                    "ea_desorption", proportional_parameter_changes, &
                                    ea_desorption)

    CALL Set_quincy_parameter_value(this%ea_hsc, &
                                    this%ea_hsc, &
                                    "ea_hsc", proportional_parameter_changes, &
                                    ea_hsc)

    CALL Set_quincy_parameter_value(this%ea_decomposition, &
                                    this%ea_decomposition, &
                                    "ea_decomposition", proportional_parameter_changes, &
                                    ea_decomposition)

    CALL Set_quincy_parameter_value(this%ed_decomposition, &
                                    this%ed_decomposition, &
                                    "ed_decomposition", proportional_parameter_changes, &
                                    ed_decomposition)

    CALL Set_quincy_parameter_value(this%t_ref_decomposition, &   ! in function calculated for proportional changes as:
                                    this%t_ref_decomposition, &   ! def_val = 273.15 + (def_val - 273.15) * nml_prop_chng
                                    "t_ref_decomposition", proportional_parameter_changes, &
                                    t_ref_decomposition)

    CALL Set_quincy_parameter_value(this%k1_afps / 3.0_wp, &
                                    this%k1_afps, &
                                    "k1_afps", proportional_parameter_changes, &
                                    k1_afps)

    CALL Set_quincy_parameter_value(this%k2_afps, &
                                    this%k2_afps, &
                                    "k2_afps", proportional_parameter_changes, &
                                    k2_afps)

    CALL Set_quincy_parameter_value(this%kmr_depolymerisation, &
                                    this%kmr_depolymerisation, &
                                    "kmr_depolymerisation", proportional_parameter_changes, &
                                    kmr_depolymerisation)

    CALL Set_quincy_parameter_value(this%kmr_mic_uptake, &
                                    this%kmr_mic_uptake, &
                                    "kmr_mic_uptake", proportional_parameter_changes, &
                                    kmr_mic_uptake)

    CALL Set_quincy_parameter_value(this%kmr_hsc, &
                                    this%kmr_hsc, &
                                    "kmr_hsc", proportional_parameter_changes, &
                                    kmr_hsc)

    CALL Set_quincy_parameter_value(this%kmr_sorption, &
                                    this%kmr_sorption, &
                                    "kmr_sorption", proportional_parameter_changes, &
                                    kmr_sorption)

    CALL Set_quincy_parameter_value(this%kmr_decomposition, &
                                    this%kmr_decomposition, &
                                    "kmr_decomposition", proportional_parameter_changes, &
                                    kmr_decomposition)

    CALL Set_quincy_parameter_value(this%k_diff_org / one_year / one_day, &
                                    this%k_diff_org, &
                                    "k_diff_org", proportional_parameter_changes, &
                                    k_diff_org)

    CALL Set_quincy_parameter_value(this%rho_bulk_org, &
                                    this%rho_bulk_org, &
                                    "rho_bulk_org", proportional_parameter_changes, &
                                    rho_bulk_org)

    CALL Set_quincy_parameter_value(this%carbon_per_dryweight_SOM, &
                                    this%carbon_per_dryweight_SOM, &
                                    "carbon_per_dryweight_SOM", proportional_parameter_changes, &
                                    carbon_per_dryweight_SOM)

    CALL Set_quincy_parameter_value(this%k_adsorpt_dom / one_year / one_day * molar_mass_C / 1000._wp, &
                                    this%k_adsorpt_dom, &
                                    "k_adsorpt_dom", proportional_parameter_changes, &
                                    k_adsorpt_dom)

    CALL Set_quincy_parameter_value(this%k_desorpt_dom / one_year / one_day / 1000._wp, &
                                    this%k_desorpt_dom, &
                                    "k_desorpt_dom", proportional_parameter_changes, &
                                    k_desorpt_dom)

    CALL Set_quincy_parameter_value(this%k_adsorpt_det / one_year / one_day * molar_mass_C / 1000._wp, &
                                    this%k_adsorpt_det, &
                                    "k_adsorpt_det", proportional_parameter_changes, &
                                    k_adsorpt_det)

    CALL Set_quincy_parameter_value(this%k_desorpt_det / one_year / one_day / 1000._wp, &
                                    this%k_desorpt_det, &
                                    "k_desorpt_det", proportional_parameter_changes, &
                                    k_desorpt_det)

    CALL Set_quincy_parameter_value(this%k_adsorpt_nh4, &
                                    this%k_adsorpt_nh4, &
                                    "k_adsorpt_nh4", proportional_parameter_changes, &
                                    k_adsorpt_nh4)

    CALL Set_quincy_parameter_value(this%k_desorpt_nh4 / 30.0_wp / one_day, &
                                    this%k_desorpt_nh4, &
                                    "k_desorpt_nh4", proportional_parameter_changes, &
                                    k_desorpt_nh4)

    CALL Set_quincy_parameter_value(this%k_adsorpt_po4 * 8.0_wp / 9.0_wp / 30.0_wp / one_day, &
                                    this%k_adsorpt_po4, &
                                    "k_adsorpt_po4", proportional_parameter_changes, &
                                    k_adsorpt_po4)

    CALL Set_quincy_parameter_value(this%k_desorpt_po4 / 30.0_wp / one_day, &
                                    this%k_desorpt_po4, &
                                    "k_desorpt_po4", proportional_parameter_changes, &
                                    k_desorpt_po4)

    CALL Set_quincy_parameter_value(this%k_adsorpt_f_po4 / 3600._wp, &
                                    this%k_adsorpt_f_po4, &
                                    "k_adsorpt_f_po4", proportional_parameter_changes, &
                                    k_adsorpt_f_po4)

    CALL Set_quincy_parameter_value(this%k_desorpt_f_po4 / 3600._wp / 10._wp, &
                                    this%k_desorpt_f_po4, &
                                    "k_desorpt_f_po4", proportional_parameter_changes, &
                                    k_desorpt_f_po4)

    CALL Set_quincy_parameter_value(this%microbial_cn * molar_mass_N / molar_mass_C, &
                                    this%microbial_cn, &
                                    "microbial_cn", proportional_parameter_changes, &
                                    microbial_cn)

    CALL Set_quincy_parameter_value(this%microbial_np * molar_mass_P / molar_mass_N, &
                                    this%microbial_np, &
                                    "microbial_np", proportional_parameter_changes, &
                                    microbial_np)

    CALL Set_quincy_parameter_value(this%microbial_nue, &
                                    this%microbial_nue, &
                                    "microbial_nue", proportional_parameter_changes, &
                                    microbial_nue)

    CALL Set_quincy_parameter_value(this%microbial_pue, &
                                    this%microbial_pue, &
                                    "microbial_pue", proportional_parameter_changes, &
                                    microbial_pue)

    CALL Set_quincy_parameter_value(this%microbial_cue, &
                                    this%microbial_cue, &
                                    "microbial_cue", proportional_parameter_changes, &
                                    microbial_cue)

    CALL Set_quincy_parameter_value(this%microbial_cue_max, &
                                    this%microbial_cue_max, &
                                    "microbial_cue_max", proportional_parameter_changes, &
                                    microbial_cue_max)

    CALL Set_quincy_parameter_value(this%microbial_cue_min, &
                                    this%microbial_cue_min, &
                                    "microbial_cue_min", proportional_parameter_changes, &
                                    microbial_cue_min)

    CALL Set_quincy_parameter_value(this%frac_woody2dec_litter, &
                                    this%frac_woody2dec_litter, &
                                    "frac_woody2dec_litter", proportional_parameter_changes, &
                                    frac_woody2dec_litter)

    CALL Set_quincy_parameter_value(this%frac_litter2fast_som, &
                                    this%frac_litter2fast_som, &
                                    "frac_litter2fast_som", proportional_parameter_changes, &
                                    frac_litter2fast_som)

    CALL Set_quincy_parameter_value(this%frac_fast2slow_som, &
                                    this%frac_fast2slow_som, &
                                    "frac_fast2slow_som", proportional_parameter_changes, &
                                    frac_fast2slow_som)

    CALL Set_quincy_parameter_value(this%frac_slow2fast_som, &
                                    this%frac_slow2fast_som, &
                                    "frac_slow2fast_som", proportional_parameter_changes, &
                                    frac_slow2fast_som)

    CALL Set_quincy_parameter_value(this%fc_woody2polymeric, &
                                    this%fc_woody2polymeric, &
                                    "fc_woody2polymeric", proportional_parameter_changes, &
                                    fc_woody2polymeric)

    CALL Set_quincy_parameter_value(this%fc_soluable2dom, &
                                    this%fc_soluable2dom, &
                                    "fc_soluable2dom", proportional_parameter_changes, &
                                    fc_soluable2dom)

    CALL Set_quincy_parameter_value(this%k_mort_microbial * 1.0e-8_wp, &
                                    this%k_mort_microbial, &
                                    "k_mort_microbial", proportional_parameter_changes, &
                                    k_mort_microbial)

    CALL Set_quincy_parameter_value(this%min_mic_biomass, &
                                    this%min_mic_biomass, &
                                    "min_mic_biomass", proportional_parameter_changes, &
                                    min_mic_biomass)

    CALL Set_quincy_parameter_value(this%vmax_mic_uptake / one_year / one_day * 1000000.0_wp * 10._wp, &
                                    this%vmax_mic_uptake, &
                                    "vmax_mic_uptake", proportional_parameter_changes, &
                                    vmax_mic_uptake)

    CALL Set_quincy_parameter_value(this%km_mic_uptake / molar_mass_C * 1000.0_wp, &
                                    this%km_mic_uptake, &
                                    "km_mic_uptake", proportional_parameter_changes, &
                                    km_mic_uptake)

    CALL Set_quincy_parameter_value(this%km_depolymerisation / 2._wp / molar_mass_C * 1000._wp, &
                                    this%km_depolymerisation, &
                                    "km_depolymerisation", proportional_parameter_changes, &
                                    km_depolymerisation)

    CALL Set_quincy_parameter_value(this%vmax_depolymerisation_residue / one_year / one_day, &
                                    this%vmax_depolymerisation_residue, &
                                    "vmax_depolymerisation_residue", proportional_parameter_changes, &
                                    vmax_depolymerisation_residue)

    CALL Set_quincy_parameter_value(this%vmax_depolymerisation_litter / one_year / one_day, &
                                    this%vmax_depolymerisation_litter, &
                                    "vmax_depolymerisation_litter", proportional_parameter_changes, &
                                    vmax_depolymerisation_litter)

    CALL Set_quincy_parameter_value(this%f_soluable_residue, &
                                    this%f_soluable_residue, &
                                    "f_soluable_residue", proportional_parameter_changes, &
                                    f_soluable_residue)

    CALL Set_quincy_parameter_value(this%f_soluable_leaf, &
                                    this%f_soluable_leaf, &
                                    "f_soluable_leaf", proportional_parameter_changes, &
                                    f_soluable_leaf)

    CALL Set_quincy_parameter_value(this%f_soluable_fine_root, &
                                    this%f_soluable_fine_root, &
                                    "f_soluable_fine_root", proportional_parameter_changes, &
                                    f_soluable_fine_root)

    CALL Set_quincy_parameter_value(this%f_soluable_fruit, &
                                    this%f_soluable_fruit, &
                                    "f_soluable_fruit", proportional_parameter_changes, &
                                    f_soluable_fruit)

    CALL Set_quincy_parameter_value(this%f_soluable_seed_bed, &
                                    this%f_soluable_seed_bed, &
                                    "f_soluable_seed_bed", proportional_parameter_changes, &
                                    f_soluable_seed_bed)

    CALL Set_quincy_parameter_value(this%fc_soluable_max, &
                                    this%fc_soluable_max, &
                                    "fc_soluable_max", proportional_parameter_changes, &
                                    fc_soluable_max)

    CALL Set_quincy_parameter_value(this%k_fc_soluable, &
                                    this%k_fc_soluable, &
                                    "k_fc_soluable", proportional_parameter_changes, &
                                    k_fc_soluable)

    CALL Set_quincy_parameter_value(this%k_fn_soluable, &
                                    this%k_fn_soluable, &
                                    "k_fn_soluable", proportional_parameter_changes, &
                                    k_fn_soluable)

    CALL Set_quincy_parameter_value(this%k_fp_soluable, &
                                    this%k_fp_soluable, &
                                    "k_fp_soluable", proportional_parameter_changes, &
                                    k_fp_soluable)

    CALL Set_quincy_parameter_value(this%lc_leaf_max / molar_mass_C * molar_mass_N, &
                                    this%lc_leaf_max, &
                                    "lc_leaf_max", proportional_parameter_changes, &
                                    lc_leaf_max)

    CALL Set_quincy_parameter_value(this%lc_leaf2sla / molar_mass_C * molar_mass_N, &
                                    this%lc_leaf2sla, &
                                    "lc_leaf2sla", proportional_parameter_changes, &
                                    lc_leaf2sla)

    CALL Set_quincy_parameter_value(this%lc_fine_root / molar_mass_C * molar_mass_N, &
                                    this%lc_fine_root, &
                                    "lc_fine_root", proportional_parameter_changes, &
                                    lc_fine_root)

    CALL Set_quincy_parameter_value(this%lc_coarse_root / molar_mass_C * molar_mass_N, &
                                    this%lc_coarse_root, &
                                    "lc_coarse_root", proportional_parameter_changes, &
                                    lc_coarse_root)

    CALL Set_quincy_parameter_value(this%lc_woody_litter / molar_mass_C * molar_mass_N, &
                                    this%lc_woody_litter, &
                                    "lc_woody_litter", proportional_parameter_changes, &
                                    lc_woody_litter)

    CALL Set_quincy_parameter_value(this%lc_fruit / molar_mass_C * molar_mass_N, &
                                    this%lc_fruit, &
                                    "lc_fruit", proportional_parameter_changes, &
                                    lc_fruit)

    CALL Set_quincy_parameter_value(this%lc_seed_bed / molar_mass_C * molar_mass_N, &
                                    this%lc_seed_bed, &
                                    "lc_seed_bed", proportional_parameter_changes, &
                                    lc_seed_bed)

    CALL Set_quincy_parameter_value(this%mycorrhiza_cue, &
                                    this%mycorrhiza_cue, &
                                    "mycorrhiza_cue", proportional_parameter_changes, &
                                    mycorrhiza_cue)

    CALL Set_quincy_parameter_value(this%mycorrhiza_cn * molar_mass_N / molar_mass_C, &
                                    this%mycorrhiza_cn, &
                                    "mycorrhiza_cn", proportional_parameter_changes, &
                                    mycorrhiza_cn)

    CALL Set_quincy_parameter_value(this%mycorrhiza_np * molar_mass_P / molar_mass_N, &
                                    this%mycorrhiza_np, &
                                    "mycorrhiza_np", proportional_parameter_changes, &
                                    mycorrhiza_np)

    CALL Set_quincy_parameter_value(this%f_maxresp_myc, &
                                    this%f_maxresp_myc, &
                                    "f_maxresp_myc", proportional_parameter_changes, &
                                    f_maxresp_myc)

    CALL Set_quincy_parameter_value(this%f_surface_myc2fineroot, &
                                    this%f_surface_myc2fineroot, &
                                    "f_surface_myc2fineroot", proportional_parameter_changes, &
                                    f_surface_myc2fineroot)

    CALL Set_quincy_parameter_value(this%km1_up_nh4, &
                                    this%km1_up_nh4, &
                                    "km1_up_nh4", proportional_parameter_changes, &
                                    km1_up_nh4)

    CALL Set_quincy_parameter_value(this%km2_up_nh4, &
                                    this%km2_up_nh4, &
                                    "km2_up_nh4", proportional_parameter_changes, &
                                    km2_up_nh4)

    CALL Set_quincy_parameter_value(this%km1_up_no3, &
                                    this%km1_up_no3, &
                                    "km1_up_no3", proportional_parameter_changes, &
                                    km1_up_no3)

    CALL Set_quincy_parameter_value(this%km2_up_no3, &
                                    this%km2_up_no3, &
                                    "km2_up_no3", proportional_parameter_changes, &
                                    km2_up_no3)

    CALL Set_quincy_parameter_value(this%km1_up_po4 / 0.0003_wp, &
                                    this%km1_up_po4, &
                                    "km1_up_po4", proportional_parameter_changes, &
                                    km1_up_po4)

    CALL Set_quincy_parameter_value(this%km2_up_po4 / 1000000.0_wp, &
                                    this%km2_up_po4, &
                                    "km2_up_po4", proportional_parameter_changes, &
                                    km2_up_po4)

    CALL Set_quincy_parameter_value(this%vmax_mic_uptake_nh4 / 3600.0_wp * molar_mass_C, &
                                    this%vmax_mic_uptake_nh4, &
                                    "vmax_mic_uptake_nh4", proportional_parameter_changes, &
                                    vmax_mic_uptake_nh4)

    CALL Set_quincy_parameter_value(this%km_mic_uptake_nh4 / molar_mass_N, &
                                    this%km_mic_uptake_nh4, &
                                    "km_mic_uptake_nh4", proportional_parameter_changes, &
                                    km_mic_uptake_nh4)

    CALL Set_quincy_parameter_value(this%vmax_mic_uptake_no3 / 3600.0_wp * molar_mass_C, &
                                    this%vmax_mic_uptake_no3, &
                                    "vmax_mic_uptake_no3", proportional_parameter_changes, &
                                    vmax_mic_uptake_no3)

    CALL Set_quincy_parameter_value(this%km_mic_uptake_no3 / molar_mass_N, &
                                    this%km_mic_uptake_no3, &
                                    "km_mic_uptake_no3", proportional_parameter_changes, &
                                    km_mic_uptake_no3)

    CALL Set_quincy_parameter_value(this%f_mic_carrier_n / 1000.0_wp, &
                                    this%f_mic_carrier_n, &
                                    "f_mic_carrier_n", proportional_parameter_changes, &
                                    f_mic_carrier_n)

    CALL Set_quincy_parameter_value(this%vmax_plant_uptake_nh4 / 3600.0_wp * molar_mass_C, &
                                    this%vmax_plant_uptake_nh4, &
                                    "vmax_plant_uptake_nh4", proportional_parameter_changes, &
                                    vmax_plant_uptake_nh4)

    CALL Set_quincy_parameter_value(this%km_plant_uptake_nh4 / molar_mass_N, &
                                    this%km_plant_uptake_nh4, &
                                    "km_plant_uptake_nh4", proportional_parameter_changes, &
                                    km_plant_uptake_nh4)

    CALL Set_quincy_parameter_value(this%vmax_plant_uptake_no3 / 3600.0_wp * molar_mass_C, &
                                    this%vmax_plant_uptake_no3, &
                                    "vmax_plant_uptake_no3", proportional_parameter_changes, &
                                    vmax_plant_uptake_no3)

    CALL Set_quincy_parameter_value(this%km_plant_uptake_no3 / molar_mass_N, &
                                    this%km_plant_uptake_no3, &
                                    "km_plant_uptake_no3", proportional_parameter_changes, &
                                    km_plant_uptake_no3)

    CALL Set_quincy_parameter_value(this%f_plant_carrier_n / 10000.0_wp, &
                                    this%f_plant_carrier_n, &
                                    "f_plant_carrier_n", proportional_parameter_changes, &
                                    f_plant_carrier_n)

    CALL Set_quincy_parameter_value(this%vmax_mic_uptake_po4 / 3600.0_wp * molar_mass_C, &
                                    this%vmax_mic_uptake_po4, &
                                    "vmax_mic_uptake_po4", proportional_parameter_changes, &
                                    vmax_mic_uptake_po4)

    CALL Set_quincy_parameter_value(this%km_mic_uptake_po4 / molar_mass_P, &
                                    this%km_mic_uptake_po4, &
                                    "km_mic_uptake_po4", proportional_parameter_changes, &
                                    km_mic_uptake_po4)

    CALL Set_quincy_parameter_value(this%vmax_plant_uptake_po4 / 60.0_wp * molar_mass_C, &
                                    this%vmax_plant_uptake_po4, &
                                    "vmax_plant_uptake_po4", proportional_parameter_changes, &
                                    vmax_plant_uptake_po4)

    CALL Set_quincy_parameter_value(this%km_plant_uptake_po4 / molar_mass_P, &
                                    this%km_plant_uptake_po4, &
                                    "km_plant_uptake_po4", proportional_parameter_changes, &
                                    km_plant_uptake_po4)

    CALL Set_quincy_parameter_value(this%k_mic_uptake_nh4 * 24._wp / one_day, &
                                    this%k_mic_uptake_nh4, &
                                    "k_mic_uptake_nh4", proportional_parameter_changes, &
                                    k_mic_uptake_nh4)

    CALL Set_quincy_parameter_value(this%k_plant_uptake_nh4 * 24._wp / one_day, &
                                    this%k_plant_uptake_nh4, &
                                    "k_plant_uptake_nh4", proportional_parameter_changes, &
                                    k_plant_uptake_nh4)

    CALL Set_quincy_parameter_value(this%k_mic_uptake_no3 * 24._wp / one_day, &
                                    this%k_mic_uptake_no3, &
                                    "k_mic_uptake_no3", proportional_parameter_changes, &
                                    k_mic_uptake_no3)

    CALL Set_quincy_parameter_value(this%k_plant_uptake_no3 * 24._wp / one_day, &
                                    this%k_plant_uptake_no3, &
                                    "k_plant_uptake_no3", proportional_parameter_changes, &
                                    k_plant_uptake_no3)

    CALL Set_quincy_parameter_value(this%k_mic_uptake_po4 / one_day, &
                                    this%k_mic_uptake_po4, &
                                    "k_mic_uptake_po4", proportional_parameter_changes, &
                                    k_mic_uptake_po4)

    CALL Set_quincy_parameter_value(this%k_plant_uptake_po4 / one_day, &
                                    this%k_plant_uptake_po4, &
                                    "k_plant_uptake_po4", proportional_parameter_changes, &
                                    k_plant_uptake_po4)

    CALL Set_quincy_parameter_value(this%f_mic_carrier_po4 / 100.0_wp, &
                                    this%f_mic_carrier_po4, &
                                    "f_mic_carrier_po4", proportional_parameter_changes, &
                                    f_mic_carrier_po4)

    CALL Set_quincy_parameter_value(this%f_plant_carrier_po4 / 10000.0_wp, &
                                    this%f_plant_carrier_po4, &
                                    "f_plant_carrier_po4", proportional_parameter_changes, &
                                    f_plant_carrier_po4)

    CALL Set_quincy_parameter_value(this%vmax_nitrification_nh4 / molar_mass_N / one_day * molar_mass_C, &
                                    this%vmax_nitrification_nh4, &
                                    "vmax_nitrification_nh4", proportional_parameter_changes, &
                                    vmax_nitrification_nh4)

    CALL Set_quincy_parameter_value(this%vmax_denitrification_no3 / molar_mass_N / one_day * molar_mass_C, &
                                    this%vmax_denitrification_no3, &
                                    "vmax_denitrification_no3", proportional_parameter_changes, &
                                    vmax_denitrification_no3)

    CALL Set_quincy_parameter_value(this%km_nitrification_nh4 / molar_mass_N, &
                                    this%km_nitrification_nh4, &
                                    "km_nitrification_nh4", proportional_parameter_changes, &
                                    km_nitrification_nh4)

    CALL Set_quincy_parameter_value(this%km_denit_no3 / molar_mass_N, &
                                    this%km_denit_no3, &
                                    "km_denit_no3", proportional_parameter_changes, &
                                    km_denit_no3)

    CALL Set_quincy_parameter_value(this%f_nit_carrier_nh4 / 10000.0_wp, &
                                    this%f_nit_carrier_nh4, &
                                    "f_nit_carrier_nh4", proportional_parameter_changes, &
                                    f_nit_carrier_nh4)

    CALL Set_quincy_parameter_value(this%f_denit_carrier_no3 / 10000.0_wp, &
                                    this%f_denit_carrier_no3, &
                                    "f_denit_carrier_no3", proportional_parameter_changes, &
                                    f_denit_carrier_no3)

    CALL Set_quincy_parameter_value(this%qmax_po4_mineral / molar_mass_P, &
                                    this%qmax_po4_mineral, &
                                    "qmax_po4_mineral", proportional_parameter_changes, &
                                    qmax_po4_mineral)

    CALL Set_quincy_parameter_value(this%qmax_nh4_clay / molar_mass_N, &
                                    this%qmax_nh4_clay, &
                                    "qmax_nh4_clay", proportional_parameter_changes, &
                                    qmax_nh4_clay)

    CALL Set_quincy_parameter_value(this%km_adsorpt_mineral_po4 / 1000._wp / molar_mass_P, &
                                    this%km_adsorpt_mineral_po4, &
                                    "km_adsorpt_mineral_po4", proportional_parameter_changes, &
                                    km_adsorpt_mineral_po4)

    CALL Set_quincy_parameter_value(this%km_adsorpt_OM_po4 / 1000._wp / molar_mass_P, &
                                    this%km_adsorpt_OM_po4, &
                                    "km_adsorpt_OM_po4", proportional_parameter_changes, &
                                    km_adsorpt_OM_po4)

    CALL Set_quincy_parameter_value(this%ea_nitrification, &
                                    this%ea_nitrification, &
                                    "ea_nitrification", proportional_parameter_changes, &
                                    ea_nitrification)

    CALL Set_quincy_parameter_value(this%ed_nitrification, &
                                    this%ed_nitrification, &
                                    "ed_nitrification", proportional_parameter_changes, &
                                    ed_nitrification)

    CALL Set_quincy_parameter_value(this%t_opt_nitrification, &   ! in function calculated for proportional changes as:
                                    this%t_opt_nitrification, &   ! def_val = 273.15 + (def_val - 273.15) * nml_prop_chng
                                    "t_opt_nitrification", proportional_parameter_changes, &
                                    t_opt_nitrification)

    CALL Set_quincy_parameter_value(this%f_nit_noy, &
                                    this%f_nit_noy, &
                                    "f_nit_noy", proportional_parameter_changes, &
                                    f_nit_noy)

    CALL Set_quincy_parameter_value(this%f_nit_n2o, &
                                    this%f_nit_n2o, &
                                    "f_nit_n2o", proportional_parameter_changes, &
                                    f_nit_n2o)

    CALL Set_quincy_parameter_value(this%ea_denitrification, &
                                    this%ea_denitrification, &
                                    "ea_denitrification", proportional_parameter_changes, &
                                    ea_denitrification)

    CALL Set_quincy_parameter_value(this%vmax_nitrification / one_day * 1.0e6_wp, &
                                    this%vmax_nitrification, &
                                    "vmax_nitrification", proportional_parameter_changes, &
                                    vmax_nitrification)

    CALL Set_quincy_parameter_value(this%vmax_denitrification / one_day * 1.0e6_wp, &
                                    this%vmax_denitrification, &
                                    "vmax_denitrification", proportional_parameter_changes, &
                                    vmax_denitrification)

    CALL Set_quincy_parameter_value(this%km_denitrification_c, &
                                    this%km_denitrification_c, &
                                    "km_denitrification_c", proportional_parameter_changes, &
                                    km_denitrification_c)

    CALL Set_quincy_parameter_value(this%km_denitrification_no3 * 1000.0_wp / molar_mass_N, &
                                    this%km_denitrification_no3, &
                                    "km_denitrification_no3", proportional_parameter_changes, &
                                    km_denitrification_no3)

    CALL Set_quincy_parameter_value(this%f_denit_noy, &
                                    this%f_denit_noy, &
                                    "f_denit_noy", proportional_parameter_changes, &
                                    f_denit_noy)

    CALL Set_quincy_parameter_value(this%f_denit_n2o, &
                                    this%f_denit_n2o, &
                                    "f_denit_n2o", proportional_parameter_changes, &
                                    f_denit_n2o)

    CALL Set_quincy_parameter_value(this%ea_gasdiffusion, &
                                    this%ea_gasdiffusion, &
                                    "ea_gasdiffusion", proportional_parameter_changes, &
                                    ea_gasdiffusion)

    CALL Set_quincy_parameter_value(this%lambda_anvf, &
                                    this%lambda_anvf, &
                                    "lambda_anvf", proportional_parameter_changes, &
                                    lambda_anvf)

    CALL Set_quincy_parameter_value(this%k_anvf, &
                                    this%k_anvf, &
                                    "k_anvf", proportional_parameter_changes, &
                                    k_anvf)

    CALL Set_quincy_parameter_value(this%floss_nmin, &
                                    this%floss_nmin, &
                                    "floss_nmin", proportional_parameter_changes, &
                                    floss_nmin)

    CALL Set_quincy_parameter_value(this%fnit_nmin, &
                                    this%fnit_nmin, &
                                    "fnit_nmin", proportional_parameter_changes, &
                                    fnit_nmin)

    CALL Set_quincy_parameter_value(this%fleach_nh4, &
                                    this%fleach_nh4, &
                                    "fleach_nh4", proportional_parameter_changes, &
                                    fleach_nh4)

    CALL Set_quincy_parameter_value(this%fleach_po4, &
                                    this%fleach_po4, &
                                    "fleach_po4", proportional_parameter_changes, &
                                    fleach_po4)

    CALL Set_quincy_parameter_value(this%km_asymb_bnf, &
                                    this%km_asymb_bnf, &
                                    "km_asymb_bnf", proportional_parameter_changes, &
                                    km_asymb_bnf)

    CALL Set_quincy_parameter_value(this%vmax_asymb_bnf, &
                                    this%vmax_asymb_bnf, &
                                    "vmax_asymb_bnf", proportional_parameter_changes, &
                                    vmax_asymb_bnf)

    CALL Set_quincy_parameter_value(this%eta_mic_uptake_nh4, &
                                    this%eta_mic_uptake_nh4, &
                                    "eta_mic_uptake_nh4", proportional_parameter_changes, &
                                    eta_mic_uptake_nh4)

    CALL Set_quincy_parameter_value(this%eta_mic_uptake_no3, &
                                    this%eta_mic_uptake_no3, &
                                    "eta_mic_uptake_no3", proportional_parameter_changes, &
                                    eta_mic_uptake_no3)

    CALL Set_quincy_parameter_value(this%eta_plant_uptake_nh4, &
                                    this%eta_plant_uptake_nh4, &
                                    "eta_plant_uptake_nh4", proportional_parameter_changes, &
                                    eta_plant_uptake_nh4)

    CALL Set_quincy_parameter_value(this%eta_plant_uptake_no3, &
                                    this%eta_plant_uptake_no3, &
                                    "eta_plant_uptake_no3", proportional_parameter_changes, &
                                    eta_plant_uptake_no3)

    CALL Set_quincy_parameter_value(this%eta_nitrification, &
                                    this%eta_nitrification, &
                                    "eta_nitrification", proportional_parameter_changes, &
                                    eta_nitrification)

    CALL Set_quincy_parameter_value(this%eta_nitrate_production, &
                                    this%eta_nitrate_production, &
                                    "eta_nitrate_production", proportional_parameter_changes, &
                                    eta_nitrate_production)

    CALL Set_quincy_parameter_value(this%eta_denitrification, &
                                    this%eta_denitrification, &
                                    "eta_denitrification", proportional_parameter_changes, &
                                    eta_denitrification)

    CALL Set_quincy_parameter_value(this%eta_ammonification, &
                                    this%eta_ammonification, &
                                    "eta_ammonification", proportional_parameter_changes, &
                                    eta_ammonification)

    CALL Set_quincy_parameter_value(this%eta_volatilisation_nh3, &
                                    this%eta_volatilisation_nh3, &
                                    "eta_volatilisation_nh3", proportional_parameter_changes, &
                                    eta_volatilisation_nh3)

    CALL Set_quincy_parameter_value(this%k_weath_mineral / molar_mass_P / 1254.0_wp / one_year / one_day, &
                                    this%k_weath_mineral, &
                                    "k_weath_mineral", proportional_parameter_changes, &
                                    k_weath_mineral)

    CALL Set_quincy_parameter_value(this%eq_sorp_mineral_po4 / molar_mass_P, &
                                    this%eq_sorp_mineral_po4, &
                                    "eq_sorp_mineral_po4", proportional_parameter_changes, &
                                    eq_sorp_mineral_po4)

    CALL Set_quincy_parameter_value(this%k_occlude_po4 / 1000000._wp / 30._wp / one_day, &
                                    this%k_occlude_po4, &
                                    "k_occlude_po4", proportional_parameter_changes, &
                                    k_occlude_po4)

    CALL Set_quincy_parameter_value(this%vmax_biochem_po4 / 50.0_wp / one_day, &
                                    this%vmax_biochem_po4, &
                                    "vmax_biochem_po4", proportional_parameter_changes, &
                                    vmax_biochem_po4)

    CALL Set_quincy_parameter_value(this%km_mic_biochem_po4 / molar_mass_C, &
                                    this%km_mic_biochem_po4, &
                                    "km_mic_biochem_po4", proportional_parameter_changes, &
                                    km_mic_biochem_po4)

    CALL Set_quincy_parameter_value(this%km_dip_biochem_po4, &
                                    this%km_dip_biochem_po4, &
                                    "km_dip_biochem_po4", proportional_parameter_changes, &
                                    km_dip_biochem_po4)

    CALL Set_quincy_parameter_value(this%km_rootc_biochem_po4, &
                                    this%km_rootc_biochem_po4, &
                                    "km_rootc_biochem_po4", proportional_parameter_changes, &
                                    km_rootc_biochem_po4)

    CALL Set_quincy_parameter_value(this%km_rootc_weath_po4, &
                                    this%km_rootc_weath_po4, &
                                    "km_rootc_weath_po4", proportional_parameter_changes, &
                                    km_rootc_weath_po4)

    CALL Set_quincy_parameter_value(this%km_root_enzyme_weath_po4 / molar_mass_C, &
                                    this%km_root_enzyme_weath_po4, &
                                    "km_root_enzyme_weath_po4", proportional_parameter_changes, &
                                    km_root_enzyme_weath_po4)

    CALL Set_quincy_parameter_value(this%biomin_threshold_pc_ratio_mic_res / 5000.0_wp, &
                                    this%biomin_threshold_pc_ratio_mic_res, &
                                    "biomin_threshold_pc_ratio_mic_res", proportional_parameter_changes, &
                                    biomin_threshold_pc_ratio_mic_res)

    CALL Set_quincy_parameter_value(this%fc_p_res2dom, &
                                    this%fc_p_res2dom, &
                                    "fc_p_res2dom", proportional_parameter_changes, &
                                    fc_p_res2dom)

    CALL Set_quincy_parameter_value(this%fc_n_res2dom, &
                                    this%fc_n_res2dom, &
                                    "fc_n_res2dom", proportional_parameter_changes, &
                                    fc_n_res2dom)

    CALL Set_quincy_parameter_value(this%f_norg_limit, &
                                    this%f_norg_limit, &
                                    "f_norg_limit", proportional_parameter_changes, &
                                    f_norg_limit)

    CALL Set_quincy_parameter_value(this%f_tau_myc, &
                                    this%f_tau_myc, &
                                    "f_tau_myc", proportional_parameter_changes, &
                                    f_tau_myc)

    CALL Set_quincy_parameter_value(this%km_priming, &
                                    this%km_priming, &
                                    "km_priming", proportional_parameter_changes, &
                                    km_priming)

    CALL Set_quincy_parameter_value(this%km_up_org, &
                                    this%km_up_org, &
                                    "km_up_org", proportional_parameter_changes, &
                                    km_up_org)

    CALL Set_quincy_parameter_value(this%myc_root_coverage_max, &
                                    this%myc_root_coverage_max, &
                                    "myc_root_coverage_max", proportional_parameter_changes, &
                                    myc_root_coverage_max)

    CALL Set_quincy_parameter_value(this%resp_priming / one_day, &
                                    this%resp_priming, &
                                    "resp_priming", proportional_parameter_changes, &
                                    resp_priming)

    CALL Set_quincy_parameter_value(this%vmax_uptake_norg, &
                                    this%vmax_uptake_norg, &
                                    "vmax_uptake_norg", proportional_parameter_changes, &
                                    vmax_uptake_norg)

    CALL Set_quincy_parameter_value(this%k_fast_som_cn_max * molar_mass_N / molar_mass_C, &
                                    this%k_fast_som_cn_max, &
                                    "k_fast_som_cn_max", proportional_parameter_changes, &
                                    k_fast_som_cn_max)

    CALL Set_quincy_parameter_value(this%k_fast_som_cn_min * molar_mass_N / molar_mass_C, &
                                    this%k_fast_som_cn_min, &
                                    "k_fast_som_cn_min", proportional_parameter_changes, &
                                    k_fast_som_cn_min)

    CALL Set_quincy_parameter_value(this%f_fast_som_cn, &
                                    this%f_fast_som_cn, &
                                    "f_fast_som_cn", proportional_parameter_changes, &
                                    f_fast_som_cn)

    CALL Set_quincy_parameter_value(this%k_slow_som_cn * molar_mass_N / molar_mass_C, &
                                    this%k_slow_som_cn, &
                                    "k_slow_som_cn", proportional_parameter_changes, &
                                    k_slow_som_cn)

    CALL Set_quincy_parameter_value(this%k_fast_som_np * molar_mass_P / molar_mass_N, &
                                    this%k_fast_som_np, &
                                    "k_fast_som_np", proportional_parameter_changes, &
                                    k_fast_som_np)

    CALL Set_quincy_parameter_value(this%k_slow_som_np * molar_mass_P / molar_mass_N, &
                                    this%k_slow_som_np, &
                                    "k_slow_som_np", proportional_parameter_changes, &
                                    k_slow_som_np)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_p_depth, &
                                          "soil_p_depth", &
                                          dsl4jsb_Config(SB_)%soil_p_depth)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_p_labile, &
                                          "soil_p_labile", &
                                          dsl4jsb_Config(SB_)%soil_p_labile)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_p_slow, &
                                          "soil_p_slow", &
                                          dsl4jsb_Config(SB_)%soil_p_slow)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_p_occluded, &
                                          "soil_p_occluded", &
                                          dsl4jsb_Config(SB_)%soil_p_occluded)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_soil_p_primary, &
                                          "soil_p_primary", &
                                          dsl4jsb_Config(SB_)%soil_p_primary)

    CALL Modify_quincy_namelist_parameter(this%f_psensi_qmax_org_fine_particle, &
                                          "qmax_org_fine_particle", &
                                          dsl4jsb_Config(SB_)%qmax_org_fine_particle)

    !> 4.0 close output file for parameter values
    !!     init was done in section 2.0
    CLOSE(lun_param_out)

  END SUBROUTINE set_parameter_slm_values_with_list_values


  !-----------------------------------------------------------------------------------------------------
  !> read lctlib parameter values from file
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE read_parameter_lctlib_list

    USE mo_jsb_namelist_iface,    ONLY: open_nml, POSITIONED, position_nml, close_nml

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
    REAL(wp)            :: g1_medlyn_p1
    REAL(wp)            :: g1_bberry_p1
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
    REAL(wp)            :: vmax_symb_bnf_p1
    REAL(wp)            :: km_symb_bnf_p1
    REAL(wp)            :: lambda_est_light_p1
    REAL(wp)            :: k_est_light_p1
    REAL(wp)            :: seed_size_p1
    REAL(wp)            :: k1_mort_greff_p1
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
    REAL(wp)            :: g1_medlyn_p2
    REAL(wp)            :: g1_bberry_p2
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
    REAL(wp)            :: vmax_symb_bnf_p2
    REAL(wp)            :: km_symb_bnf_p2
    REAL(wp)            :: lambda_est_light_p2
    REAL(wp)            :: k_est_light_p2
    REAL(wp)            :: seed_size_p2
    REAL(wp)            :: k1_mort_greff_p2
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
    REAL(wp)            :: g1_medlyn_p3
    REAL(wp)            :: g1_bberry_p3
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
    REAL(wp)            :: vmax_symb_bnf_p3
    REAL(wp)            :: km_symb_bnf_p3
    REAL(wp)            :: lambda_est_light_p3
    REAL(wp)            :: k_est_light_p3
    REAL(wp)            :: seed_size_p3
    REAL(wp)            :: k1_mort_greff_p3
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
    REAL(wp)            :: g1_medlyn_p4
    REAL(wp)            :: g1_bberry_p4
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
    REAL(wp)            :: vmax_symb_bnf_p4
    REAL(wp)            :: km_symb_bnf_p4
    REAL(wp)            :: lambda_est_light_p4
    REAL(wp)            :: k_est_light_p4
    REAL(wp)            :: seed_size_p4
    REAL(wp)            :: k1_mort_greff_p4
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
    REAL(wp)            :: g1_medlyn_p5
    REAL(wp)            :: g1_bberry_p5
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
    REAL(wp)            :: vmax_symb_bnf_p5
    REAL(wp)            :: km_symb_bnf_p5
    REAL(wp)            :: lambda_est_light_p5
    REAL(wp)            :: k_est_light_p5
    REAL(wp)            :: seed_size_p5
    REAL(wp)            :: k1_mort_greff_p5
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
    REAL(wp)            :: g1_medlyn_p6
    REAL(wp)            :: g1_bberry_p6
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
    REAL(wp)            :: vmax_symb_bnf_p6
    REAL(wp)            :: km_symb_bnf_p6
    REAL(wp)            :: lambda_est_light_p6
    REAL(wp)            :: k_est_light_p6
    REAL(wp)            :: seed_size_p6
    REAL(wp)            :: k1_mort_greff_p6
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
    REAL(wp)            :: g1_medlyn_p7
    REAL(wp)            :: g1_bberry_p7
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
    REAL(wp)            :: vmax_symb_bnf_p7
    REAL(wp)            :: km_symb_bnf_p7
    REAL(wp)            :: lambda_est_light_p7
    REAL(wp)            :: k_est_light_p7
    REAL(wp)            :: seed_size_p7
    REAL(wp)            :: k1_mort_greff_p7
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
    REAL(wp)            :: g1_medlyn_p8
    REAL(wp)            :: g1_bberry_p8
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
    REAL(wp)            :: vmax_symb_bnf_p8
    REAL(wp)            :: km_symb_bnf_p8
    REAL(wp)            :: lambda_est_light_p8
    REAL(wp)            :: k_est_light_p8
    REAL(wp)            :: seed_size_p8
    REAL(wp)            :: k1_mort_greff_p8
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
    NAMELIST /lctlib_pft1_nml/  &
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
          g1_medlyn_p1, &
          g1_bberry_p1, &
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
          vmax_symb_bnf_p1, &
          km_symb_bnf_p1, &
          lambda_est_light_p1, &
          k_est_light_p1, &
          seed_size_p1, &
          k1_mort_greff_p1, &
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

    NAMELIST /lctlib_pft2_nml/  &
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
          g1_medlyn_p2, &
          g1_bberry_p2, &
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
          vmax_symb_bnf_p2, &
          km_symb_bnf_p2, &
          lambda_est_light_p2, &
          k_est_light_p2, &
          seed_size_p2, &
          k1_mort_greff_p2, &
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

    NAMELIST /lctlib_pft3_nml/  &
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
          g1_medlyn_p3, &
          g1_bberry_p3, &
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
          vmax_symb_bnf_p3, &
          km_symb_bnf_p3, &
          lambda_est_light_p3, &
          k_est_light_p3, &
          seed_size_p3, &
          k1_mort_greff_p3, &
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

    NAMELIST /lctlib_pft4_nml/  &
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
          g1_medlyn_p4, &
          g1_bberry_p4, &
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
          vmax_symb_bnf_p4, &
          km_symb_bnf_p4, &
          lambda_est_light_p4, &
          k_est_light_p4, &
          seed_size_p4, &
          k1_mort_greff_p4, &
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

    NAMELIST /lctlib_pft5_nml/  &
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
          g1_medlyn_p5, &
          g1_bberry_p5, &
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
          vmax_symb_bnf_p5, &
          km_symb_bnf_p5, &
          lambda_est_light_p5, &
          k_est_light_p5, &
          seed_size_p5, &
          k1_mort_greff_p5, &
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

    NAMELIST /lctlib_pft6_nml/  &
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
          g1_medlyn_p6, &
          g1_bberry_p6, &
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
          vmax_symb_bnf_p6, &
          km_symb_bnf_p6, &
          lambda_est_light_p6, &
          k_est_light_p6, &
          seed_size_p6, &
          k1_mort_greff_p6, &
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

    NAMELIST /lctlib_pft7_nml/  &
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
          g1_medlyn_p7, &
          g1_bberry_p7, &
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
          vmax_symb_bnf_p7, &
          km_symb_bnf_p7, &
          lambda_est_light_p7, &
          k_est_light_p7, &
          seed_size_p7, &
          k1_mort_greff_p7, &
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

    NAMELIST /lctlib_pft8_nml/  &
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
          g1_medlyn_p8, &
          g1_bberry_p8, &
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
          vmax_symb_bnf_p8, &
          km_symb_bnf_p8, &
          lambda_est_light_p8, &
          k_est_light_p8, &
          seed_size_p8, &
          k1_mort_greff_p8, &
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
    g1_medlyn_p1      = missval
    g1_bberry_p1      = missval
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
    vmax_symb_bnf_p1      = missval
    km_symb_bnf_p1      = missval
    lambda_est_light_p1      = missval
    k_est_light_p1      = missval
    seed_size_p1      = missval
    k1_mort_greff_p1      = missval
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
    g1_medlyn_p2      = missval
    g1_bberry_p2      = missval
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
    vmax_symb_bnf_p2      = missval
    km_symb_bnf_p2      = missval
    lambda_est_light_p2      = missval
    k_est_light_p2      = missval
    seed_size_p2      = missval
    k1_mort_greff_p2      = missval
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
    g1_medlyn_p3      = missval
    g1_bberry_p3      = missval
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
    vmax_symb_bnf_p3      = missval
    km_symb_bnf_p3      = missval
    lambda_est_light_p3      = missval
    k_est_light_p3      = missval
    seed_size_p3      = missval
    k1_mort_greff_p3      = missval
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
    g1_medlyn_p4      = missval
    g1_bberry_p4      = missval
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
    vmax_symb_bnf_p4      = missval
    km_symb_bnf_p4      = missval
    lambda_est_light_p4      = missval
    k_est_light_p4      = missval
    seed_size_p4      = missval
    k1_mort_greff_p4      = missval
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
    g1_medlyn_p5      = missval
    g1_bberry_p5      = missval
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
    vmax_symb_bnf_p5      = missval
    km_symb_bnf_p5      = missval
    lambda_est_light_p5      = missval
    k_est_light_p5      = missval
    seed_size_p5      = missval
    k1_mort_greff_p5      = missval
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
    g1_medlyn_p6      = missval
    g1_bberry_p6      = missval
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
    vmax_symb_bnf_p6      = missval
    km_symb_bnf_p6      = missval
    lambda_est_light_p6      = missval
    k_est_light_p6      = missval
    seed_size_p6      = missval
    k1_mort_greff_p6      = missval
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
    g1_medlyn_p7      = missval
    g1_bberry_p7      = missval
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
    vmax_symb_bnf_p7      = missval
    km_symb_bnf_p7      = missval
    lambda_est_light_p7      = missval
    k_est_light_p7      = missval
    seed_size_p7      = missval
    k1_mort_greff_p7      = missval
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
    g1_medlyn_p8      = missval
    g1_bberry_p8      = missval
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
    vmax_symb_bnf_p8      = missval
    km_symb_bnf_p8      = missval
    lambda_est_light_p8      = missval
    k_est_light_p8      = missval
    seed_size_p8      = missval
    k1_mort_greff_p8      = missval
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
    ! &lctlib_pft1_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft1_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft1_nml)
    CALL close_nml(nml_handler)
    ! &lctlib_pft2_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft2_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft2_nml)
    CALL close_nml(nml_handler)
    ! &lctlib_pft3_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft3_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft3_nml)
    CALL close_nml(nml_handler)
    ! &lctlib_pft4_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft4_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft4_nml)
    CALL close_nml(nml_handler)
    ! &lctlib_pft5_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft5_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft5_nml)
    CALL close_nml(nml_handler)
    ! &lctlib_pft6_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft6_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft6_nml)
    CALL close_nml(nml_handler)
    ! &lctlib_pft7_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft7_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft7_nml)
    CALL close_nml(nml_handler)
    ! &lctlib_pft8_nml
    nml_handler = open_nml(TRIM(param_pft1_from_file%filename))
    nml_unit = position_nml('lctlib_pft8_nml', nml_handler, STATUS=istat)
    IF (istat == POSITIONED) READ(nml_unit, lctlib_pft8_nml)
    CALL close_nml(nml_handler)


    ! pass values as read from file
    !
    ! lctlib_pft1_nml
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
    param_pft1_from_file%g1_medlyn   = g1_medlyn_p1
    param_pft1_from_file%g1_bberry   = g1_bberry_p1
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
    param_pft1_from_file%vmax_symb_bnf   = vmax_symb_bnf_p1
    param_pft1_from_file%km_symb_bnf   = km_symb_bnf_p1
    param_pft1_from_file%lambda_est_light   = lambda_est_light_p1
    param_pft1_from_file%k_est_light   = k_est_light_p1
    param_pft1_from_file%seed_size   = seed_size_p1
    param_pft1_from_file%k1_mort_greff   = k1_mort_greff_p1
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
    ! lctlib_pft2_nml
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
    param_pft2_from_file%g1_medlyn   = g1_medlyn_p2
    param_pft2_from_file%g1_bberry   = g1_bberry_p2
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
    param_pft2_from_file%vmax_symb_bnf   = vmax_symb_bnf_p2
    param_pft2_from_file%km_symb_bnf   = km_symb_bnf_p2
    param_pft2_from_file%lambda_est_light   = lambda_est_light_p2
    param_pft2_from_file%k_est_light   = k_est_light_p2
    param_pft2_from_file%seed_size   = seed_size_p2
    param_pft2_from_file%k1_mort_greff   = k1_mort_greff_p2
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
    ! lctlib_pft3_nml
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
    param_pft3_from_file%g1_medlyn   = g1_medlyn_p3
    param_pft3_from_file%g1_bberry   = g1_bberry_p3
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
    param_pft3_from_file%vmax_symb_bnf   = vmax_symb_bnf_p3
    param_pft3_from_file%km_symb_bnf   = km_symb_bnf_p3
    param_pft3_from_file%lambda_est_light   = lambda_est_light_p3
    param_pft3_from_file%k_est_light   = k_est_light_p3
    param_pft3_from_file%seed_size   = seed_size_p3
    param_pft3_from_file%k1_mort_greff   = k1_mort_greff_p3
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
    ! lctlib_pft4_nml
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
    param_pft4_from_file%g1_medlyn   = g1_medlyn_p4
    param_pft4_from_file%g1_bberry   = g1_bberry_p4
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
    param_pft4_from_file%vmax_symb_bnf   = vmax_symb_bnf_p4
    param_pft4_from_file%km_symb_bnf   = km_symb_bnf_p4
    param_pft4_from_file%lambda_est_light   = lambda_est_light_p4
    param_pft4_from_file%k_est_light   = k_est_light_p4
    param_pft4_from_file%seed_size   = seed_size_p4
    param_pft4_from_file%k1_mort_greff   = k1_mort_greff_p4
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
    ! lctlib_pft5_nml
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
    param_pft5_from_file%g1_medlyn   = g1_medlyn_p5
    param_pft5_from_file%g1_bberry   = g1_bberry_p5
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
    param_pft5_from_file%vmax_symb_bnf   = vmax_symb_bnf_p5
    param_pft5_from_file%km_symb_bnf   = km_symb_bnf_p5
    param_pft5_from_file%lambda_est_light   = lambda_est_light_p5
    param_pft5_from_file%k_est_light   = k_est_light_p5
    param_pft5_from_file%seed_size   = seed_size_p5
    param_pft5_from_file%k1_mort_greff   = k1_mort_greff_p5
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
    ! lctlib_pft6_nml
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
    param_pft6_from_file%g1_medlyn   = g1_medlyn_p6
    param_pft6_from_file%g1_bberry   = g1_bberry_p6
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
    param_pft6_from_file%vmax_symb_bnf   = vmax_symb_bnf_p6
    param_pft6_from_file%km_symb_bnf   = km_symb_bnf_p6
    param_pft6_from_file%lambda_est_light   = lambda_est_light_p6
    param_pft6_from_file%k_est_light   = k_est_light_p6
    param_pft6_from_file%seed_size   = seed_size_p6
    param_pft6_from_file%k1_mort_greff   = k1_mort_greff_p6
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
    ! lctlib_pft7_nml
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
    param_pft7_from_file%g1_medlyn   = g1_medlyn_p7
    param_pft7_from_file%g1_bberry   = g1_bberry_p7
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
    param_pft7_from_file%vmax_symb_bnf   = vmax_symb_bnf_p7
    param_pft7_from_file%km_symb_bnf   = km_symb_bnf_p7
    param_pft7_from_file%lambda_est_light   = lambda_est_light_p7
    param_pft7_from_file%k_est_light   = k_est_light_p7
    param_pft7_from_file%seed_size   = seed_size_p7
    param_pft7_from_file%k1_mort_greff   = k1_mort_greff_p7
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
    ! lctlib_pft8_nml
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
    param_pft8_from_file%g1_medlyn   = g1_medlyn_p8
    param_pft8_from_file%g1_bberry   = g1_bberry_p8
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
    param_pft8_from_file%vmax_symb_bnf   = vmax_symb_bnf_p8
    param_pft8_from_file%km_symb_bnf   = km_symb_bnf_p8
    param_pft8_from_file%lambda_est_light   = lambda_est_light_p8
    param_pft8_from_file%k_est_light   = k_est_light_p8
    param_pft8_from_file%seed_size   = seed_size_p8
    param_pft8_from_file%k1_mort_greff   = k1_mort_greff_p8
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
    USE mo_jsb_process_class,               ONLY: VEG_
    ! process configurations (t_PROC_config)
    dsl4jsb_Use_config(VEG_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_model), INTENT(inout) :: model
    ! ---------------------------
    ! 0.2 Local
    INTEGER                     :: pft_id
    LOGICAL                     :: prop_parameter_changes
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_parameter_lctlib_values_with_list_values'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(VEG_)
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    !
    ! ---------------------------
    ! 0.5 Get Memory
    !
    ! Get process config
    dsl4jsb_Get_config(VEG_)


    !> 0.9 init local variables
    !!
    prop_parameter_changes = model%config%set_parameter_values_proportional
    pft_id                 = dsl4jsb_Config(VEG_)%pft_id


    !> 1.0 read parameter lctlib values from file
    !!
    CALL read_parameter_lctlib_list

    !> 2.0 replace parameter lctlib values for the particular PFT X that is simulated
    !!
    !! lctlib(4)%sla is needed for calculation of k_rtos
    SELECT CASE(pft_id)
      CASE(1)
        CALL change_lctlib_param_values_pft_x(model%lctlib(1), param_pft1_from_file, model%lctlib(4)%sla, prop_parameter_changes)
      CASE(2)
        CALL change_lctlib_param_values_pft_x(model%lctlib(2), param_pft2_from_file, model%lctlib(4)%sla, prop_parameter_changes)
      CASE(3)
        CALL change_lctlib_param_values_pft_x(model%lctlib(3), param_pft3_from_file, model%lctlib(4)%sla, prop_parameter_changes)
      CASE(4)
        CALL change_lctlib_param_values_pft_x(model%lctlib(4), param_pft4_from_file, model%lctlib(4)%sla, prop_parameter_changes)
      CASE(5)
        CALL change_lctlib_param_values_pft_x(model%lctlib(5), param_pft5_from_file, model%lctlib(4)%sla, prop_parameter_changes)
      CASE(6)
        CALL change_lctlib_param_values_pft_x(model%lctlib(6), param_pft6_from_file, model%lctlib(4)%sla, prop_parameter_changes)
      CASE(7)
        CALL change_lctlib_param_values_pft_x(model%lctlib(7), param_pft7_from_file, model%lctlib(4)%sla, prop_parameter_changes)
      CASE(8)
        CALL change_lctlib_param_values_pft_x(model%lctlib(8), param_pft8_from_file, model%lctlib(4)%sla, prop_parameter_changes)
    END SELECT

  END SUBROUTINE set_parameter_lctlib_values_with_list_values


  !-----------------------------------------------------------------------------------------------------
  !> replace hard coded lctlib parameter values for one PFT with values from file
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE change_lctlib_param_values_pft_x(lctlib, param_pft_from_file, pft4_sla, proportional_parameter_changes)

    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_math_constants,      ONLY: one_day, one_year
    USE mo_jsb_physical_constants,  ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero
    USE mo_veg_constants,           ONLY: carbon_per_dryweight_leaf, sm2lm_grass, itree

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
    CHARACTER(33)               :: fname                      ! filename
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':change_lctlib_param_values_pft_x'



    !> 0.9 set sla_helper to the sla value of pft4 always
    !!
    pft4_sla_hlp = pft4_sla


    !> 1.0 open and init output file for parameter values
    !!     file is closed in section 3.0
    fname = 'parameter_sensi_lctlib_values.txt'
    OPEN (UNIT=lun_param_out, FILE=fname, ACTION="write", STATUS="unknown")
    WRITE(lun_param_out,'(a)') 'parName parVal parValChanged'


    !> 2.0 replace lctlib parameter-value for the particular PFT if a new value is provided
    !!
    CALL Set_quincy_parameter_value(1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * param_pft_from_file%sla, &
                                    param_pft_from_file%sla, &
                                    "lctlib_sla", proportional_parameter_changes, &
                                    lctlib%sla)
    ! if the PFT of this simulation is 4 re-set the sla_helper variable (just in case the lctlib%sla was changed)
    IF (lctlib%lctnumber == 4) pft4_sla_hlp = lctlib%sla

    CALL Set_quincy_parameter_value(param_pft_from_file%sigma_vis, &
                                    param_pft_from_file%sigma_vis, &
                                    "lctlib_sigma_vis", proportional_parameter_changes, &
                                    lctlib%sigma_vis)

    CALL Set_quincy_parameter_value(param_pft_from_file%sigma_nir, &
                                    param_pft_from_file%sigma_nir, &
                                    "lctlib_sigma_nir", proportional_parameter_changes, &
                                    lctlib%sigma_nir)

    CALL Set_quincy_parameter_value(param_pft_from_file%omega_clumping, &
                                    param_pft_from_file%omega_clumping, &
                                    "lctlib_omega_clumping", proportional_parameter_changes, &
                                    lctlib%omega_clumping)

    CALL Set_quincy_parameter_value(param_pft_from_file%crown_shape_factor, &
                                    param_pft_from_file%crown_shape_factor, &
                                    "lctlib_crown_shape_factor", proportional_parameter_changes, &
                                    lctlib%crown_shape_factor)

    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / &
                                      param_pft_from_file%cn_leaf, &   ! line break - this part belongs to prev. line!
                                    param_pft_from_file%cn_leaf, &
                                    "lctlib_cn_leaf", proportional_parameter_changes, &
                                    lctlib%cn_leaf)

    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / &
                                      param_pft_from_file%cn_leaf_min, &   ! line break - this part belongs to prev. line!
                                    param_pft_from_file%cn_leaf_min, &
                                    "lctlib_cn_leaf_min", proportional_parameter_changes, &
                                    lctlib%cn_leaf_min)

    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / &
                                      param_pft_from_file%cn_leaf_max, &   ! line break - this part belongs to prev. line!
                                    param_pft_from_file%cn_leaf_max, &
                                    "lctlib_cn_leaf_max", proportional_parameter_changes, &
                                    lctlib%cn_leaf_max)

    ! not resolving the dependency to cn_leaf (as re-calc is not so easy, in case proportional changes are made, as the original input lctlib parameter is not available here)
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib%cn_leaf / &
                                      param_pft_from_file%np_leaf, &   ! line break - this part belongs to prev. line!
                                    param_pft_from_file%np_leaf, &
                                    "lctlib_np_leaf", proportional_parameter_changes, &
                                    lctlib%np_leaf)

    ! not resolving the dependency to cn_leaf (as re-calc is not so easy, in case proportional changes are made, as the original input lctlib parameter is not available here)
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib%cn_leaf / &
                                      param_pft_from_file%np_leaf_min, &   ! line break - this part belongs to prev. line!
                                    param_pft_from_file%np_leaf_min, &
                                    "lctlib_np_leaf_min", proportional_parameter_changes, &
                                    lctlib%np_leaf_min)

    ! not resolving the dependency to cn_leaf (as re-calc is not so easy, in case proportional changes are made, as the original input lctlib parameter is not available here)
    CALL Set_quincy_parameter_value(carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib%cn_leaf / &
                                      param_pft_from_file%np_leaf_max, &   ! line break - this part belongs to prev. line!
                                    param_pft_from_file%np_leaf_max, &
                                    "lctlib_np_leaf_max", proportional_parameter_changes, &
                                    lctlib%np_leaf_max)

    CALL Set_quincy_parameter_value(param_pft_from_file%k0_fn_struc, &
                                    param_pft_from_file%k0_fn_struc, &
                                    "lctlib_k0_fn_struc", proportional_parameter_changes, &
                                    lctlib%k0_fn_struc)

    ! is always recalculated in case any of its dependencies has changed
    lctlib%fn_oth_min = lctlib%k0_fn_struc - 5.0_wp / molar_mass_N &
      &                 * 0.5_wp * 1.e3_wp * 71.4_wp * molar_mass_N / 1.e6_wp
    ! no additional parameter changes to fn_oth_min
    ! CALL Set_quincy_parameter_value(param_pft_from_file%k0_fn_struc - 5.0_wp / molar_mass_N * 0.5_wp * 1.e3_wp * 71.4_wp * &
    !                                   molar_mass_N / 1.e6_wp, &   ! line break - this part belongs to prev. line!
    !                                 param_pft_from_file%fn_oth_min, &
    !                                 "lctlib_fn_oth_min", proportional_parameter_changes, &
    !                                 lctlib%fn_oth_min)

    CALL Set_quincy_parameter_value(param_pft_from_file%t_jmax_opt, &
                                    param_pft_from_file%t_jmax_opt, &
                                    "lctlib_t_jmax_opt", proportional_parameter_changes, &
                                    lctlib%t_jmax_opt)

    CALL Set_quincy_parameter_value(param_pft_from_file%t_jmax_omega, &
                                    param_pft_from_file%t_jmax_omega, &
                                    "lctlib_t_jmax_omega", proportional_parameter_changes, &
                                    lctlib%t_jmax_omega)

    CALL Set_quincy_parameter_value(param_pft_from_file%g0, &
                                    param_pft_from_file%g0, &
                                    "lctlib_g0", proportional_parameter_changes, &
                                    lctlib%g0)

    CALL Set_quincy_parameter_value(param_pft_from_file%g1_medlyn, &
                                    param_pft_from_file%g1_medlyn, &
                                    "lctlib_g1_medlyn", proportional_parameter_changes, &
                                    lctlib%g1_medlyn)

    CALL Set_quincy_parameter_value(param_pft_from_file%g1_bberry, &
                                    param_pft_from_file%g1_bberry, &
                                    "lctlib_g1_bberry", proportional_parameter_changes, &
                                    lctlib%g1_bberry)

    CALL Set_quincy_parameter_value(param_pft_from_file%gmin, &
                                    param_pft_from_file%gmin, &
                                    "lctlib_gmin", proportional_parameter_changes, &
                                    lctlib%gmin)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_leaf / 12.0_wp, &
                                    param_pft_from_file%tau_leaf, &
                                    "lctlib_tau_leaf", proportional_parameter_changes, &
                                    lctlib%tau_leaf)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_fine_root, &
                                    param_pft_from_file%tau_fine_root, &
                                    "lctlib_tau_fine_root", proportional_parameter_changes, &
                                    lctlib%tau_fine_root)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_coarse_root, &
                                    param_pft_from_file%tau_coarse_root, &
                                    "lctlib_tau_coarse_root", proportional_parameter_changes, &
                                    lctlib%tau_coarse_root)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_branch, &
                                    param_pft_from_file%tau_branch, &
                                    "lctlib_tau_branch", proportional_parameter_changes, &
                                    lctlib%tau_branch)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_sap_wood, &
                                    param_pft_from_file%tau_sap_wood, &
                                    "lctlib_tau_sap_wood", proportional_parameter_changes, &
                                    lctlib%tau_sap_wood)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_fruit, &
                                    param_pft_from_file%tau_fruit, &
                                    "lctlib_tau_fruit", proportional_parameter_changes, &
                                    lctlib%tau_fruit)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_seed_litter, &
                                    param_pft_from_file%tau_seed_litter, &
                                    "lctlib_tau_seed_litter", proportional_parameter_changes, &
                                    lctlib%tau_seed_litter)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_seed_est, &
                                    param_pft_from_file%tau_seed_est, &
                                    "lctlib_tau_seed_est", proportional_parameter_changes, &
                                    lctlib%tau_seed_est)

    CALL Set_quincy_parameter_value(param_pft_from_file%tau_mycorrhiza, &
                                    param_pft_from_file%tau_mycorrhiza, &
                                    "lctlib_tau_mycorrhiza", proportional_parameter_changes, &
                                    lctlib%tau_mycorrhiza)

    CALL Set_quincy_parameter_value(param_pft_from_file%vmax_uptake_n, &
                                    param_pft_from_file%vmax_uptake_n, &
                                    "lctlib_vmax_uptake_n", proportional_parameter_changes, &
                                    lctlib%vmax_uptake_n)

    CALL Set_quincy_parameter_value(param_pft_from_file%vmax_uptake_p, &
                                    param_pft_from_file%vmax_uptake_p, &
                                    "lctlib_vmax_uptake_p", proportional_parameter_changes, &
                                    lctlib%vmax_uptake_p)

    CALL Set_quincy_parameter_value(param_pft_from_file%bnf_base / molar_mass_N / one_day / one_year * 1.e6_wp, &
                                    param_pft_from_file%bnf_base, &
                                    "lctlib_bnf_base", proportional_parameter_changes, &
                                    lctlib%bnf_base)

    CALL Set_quincy_parameter_value(param_pft_from_file%vmax_symb_bnf / (molar_mass_N / molar_mass_C) &
                                    / one_day / one_year * 1.e6_wp, &
                                    param_pft_from_file%vmax_symb_bnf, &
                                    "lctlib_vmax_symb_bnf", proportional_parameter_changes, &
                                    lctlib%vmax_symb_bnf)

    CALL Set_quincy_parameter_value(param_pft_from_file%km_symb_bnf / (molar_mass_C / molar_mass_N), &
                                    param_pft_from_file%km_symb_bnf, &
                                    "lctlib_km_symb_bnf", proportional_parameter_changes, &
                                    lctlib%km_symb_bnf)

    CALL Set_quincy_parameter_value(param_pft_from_file%lambda_est_light, &
                                    param_pft_from_file%lambda_est_light, &
                                    "lctlib_lambda_est_light", proportional_parameter_changes, &
                                    lctlib%lambda_est_light)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_est_light, &
                                    param_pft_from_file%k_est_light, &
                                    "lctlib_k_est_light", proportional_parameter_changes, &
                                    lctlib%k_est_light)

    CALL Set_quincy_parameter_value(param_pft_from_file%seed_size, &
                                    param_pft_from_file%seed_size, &
                                    "lctlib_seed_size", proportional_parameter_changes, &
                                    lctlib%seed_size)

    CALL Set_quincy_parameter_value(param_pft_from_file%k1_mort_greff, &
                                    param_pft_from_file%k1_mort_greff, &
                                    "lctlib_k1_mort_greff", proportional_parameter_changes, &
                                    lctlib%k1_mort_greff)

    CALL Set_quincy_parameter_value(param_pft_from_file%beta_soil_flush, &
                                    param_pft_from_file%beta_soil_flush, &
                                    "lctlib_beta_soil_flush", proportional_parameter_changes, &
                                    lctlib%beta_soil_flush)
    ! beta_soil_flush may not be > 1
    IF (lctlib%beta_soil_flush > 1.0_wp) lctlib%beta_soil_flush = 1.0_wp

    CALL Set_quincy_parameter_value(param_pft_from_file%beta_soil_senescence, &
                                    param_pft_from_file%beta_soil_senescence, &
                                    "lctlib_beta_soil_senescence", proportional_parameter_changes, &
                                    lctlib%beta_soil_senescence)
    ! beta_soil_senescence may not be > 1
    IF (lctlib%beta_soil_senescence > 1.0_wp) lctlib%beta_soil_senescence = 1.0_wp

    CALL Set_quincy_parameter_value(param_pft_from_file%gdd_req_max, &
                                    param_pft_from_file%gdd_req_max, &
                                    "lctlib_gdd_req_max", proportional_parameter_changes, &
                                    lctlib%gdd_req_max)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_gdd_dormance, &
                                    param_pft_from_file%k_gdd_dormance, &
                                    "lctlib_k_gdd_dormance", proportional_parameter_changes, &
                                    lctlib%k_gdd_dormance)

    CALL Set_quincy_parameter_value(param_pft_from_file%t_air_senescence + Tzero, &
                                    param_pft_from_file%t_air_senescence, &
                                    "lctlib_t_air_senescence", proportional_parameter_changes, &
                                    lctlib%t_air_senescence)

    CALL Set_quincy_parameter_value(param_pft_from_file%min_leaf_age, &
                                    param_pft_from_file%min_leaf_age, &
                                    "lctlib_min_leaf_age", proportional_parameter_changes, &
                                    lctlib%min_leaf_age)

    CALL Set_quincy_parameter_value(param_pft_from_file%frac_sapwood_branch, &
                                    param_pft_from_file%frac_sapwood_branch, &
                                    "lctlib_frac_sapwood_branch", proportional_parameter_changes, &
                                    lctlib%frac_sapwood_branch)
    ! frac_sapwood_branch may not be > 1
    IF (lctlib%frac_sapwood_branch > 1.0_wp) lctlib%frac_sapwood_branch = 1.0_wp

    CALL Set_quincy_parameter_value(param_pft_from_file%wood_density * 1.e6_wp / molar_mass_C, &
                                    param_pft_from_file%wood_density, &
                                    "lctlib_wood_density", proportional_parameter_changes, &
                                    lctlib%wood_density)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_latosa, &
                                    param_pft_from_file%k_latosa, &
                                    "lctlib_k_latosa", proportional_parameter_changes, &
                                    lctlib%k_latosa)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_crtos, &
                                    param_pft_from_file%k_crtos, &
                                    "lctlib_k_crtos", proportional_parameter_changes, &
                                    lctlib%k_crtos)

    ! is always recalculated in case any of its dependencies has changed
    IF(lctlib%growthform == itree) THEN
      lctlib%k_rtos = 1.75_wp * (pft4_sla_hlp * lctlib%wood_density) / lctlib%k_latosa
    ELSE
      lctlib%k_rtos = 0.5_wp / sm2lm_grass
    ENDIF
    ! no additional parameter changes to k_rtos
    ! CALL Set_quincy_parameter_value(k_rtos_hlp, &
    !                                 param_pft_from_file%k_rtos, &
    !                                 "lctlib_k_rtos", proportional_parameter_changes, &
    !                                 lctlib%k_rtos)

    CALL Set_quincy_parameter_value(param_pft_from_file%k2_fruit_alloc, &
                                    param_pft_from_file%k2_fruit_alloc, &
                                    "lctlib_k2_fruit_alloc", proportional_parameter_changes, &
                                    lctlib%k2_fruit_alloc)

    CALL Set_quincy_parameter_value(param_pft_from_file%allom_k1, &
                                    param_pft_from_file%allom_k1, &
                                    "lctlib_allom_k1", proportional_parameter_changes, &
                                    lctlib%allom_k1)

    CALL Set_quincy_parameter_value(param_pft_from_file%allom_k2, &
                                    param_pft_from_file%allom_k2, &
                                    "lctlib_allom_k2", proportional_parameter_changes, &
                                    lctlib%allom_k2)

    CALL Set_quincy_parameter_value(param_pft_from_file%phi_leaf_min, &
                                    param_pft_from_file%phi_leaf_min, &
                                    "lctlib_phi_leaf_min", proportional_parameter_changes, &
                                    lctlib%phi_leaf_min)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_root * molar_mass_C * 1.0e-10_wp, &
                                    param_pft_from_file%k_root, &
                                    "lctlib_k_root", proportional_parameter_changes, &
                                    lctlib%k_root)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_sapwood * 1.0e-3_wp, &
                                    param_pft_from_file%k_sapwood, &
                                    "lctlib_k_sapwood", proportional_parameter_changes, &
                                    lctlib%k_sapwood)

    ! is always recalculated in case any of its dependencies has changed
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
    !                                 "lctlib_c0_allom", proportional_parameter_changes, &
    !                                 lctlib%c0_allom)

    CALL Set_quincy_parameter_value(param_pft_from_file%fstore_target, &
                                    param_pft_from_file%fstore_target, &
                                    "lctlib_fstore_target", proportional_parameter_changes, &
                                    lctlib%fstore_target)

    CALL Set_quincy_parameter_value(param_pft_from_file%k_root_dist, &
                                    param_pft_from_file%k_root_dist, &
                                    "lctlib_k_root_dist", proportional_parameter_changes, &
                                    lctlib%k_root_dist)



    !> 3.0 close output file for parameter values
    !!     init was done in section 1.0
    CLOSE(lun_param_out)

  END SUBROUTINE change_lctlib_param_values_pft_x

  !-----------------------------------------------------------------------------------------------------
  !> return new QUINCY parameter value as defined by namellist \n
  !!  namelist created from generate_parameter_table.sh: parameter_slm_run.list \n
  !!  i.e., replace hard coded parameter value
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Set_quincy_parameter_value(nml_abs_val, nml_prop_chng, param_name, proportional, def_val)

    USE mo_kind,                ONLY: wp
    USE mo_jsb_math_constants,  ONLY: eps8
    USE mo_jsb_io,              ONLY: missval  ! = -9.E+15_wp

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    REAL(wp),         INTENT(in)    :: nml_abs_val        !< new value from namelist, potentially re-calculated with other parameters consistently with mo_PROC_constants
    REAL(wp),         INTENT(in)    :: nml_prop_chng      !< new value from namelist, not changed, used here for proportional change of def_val x from namelist [0 < x < 2]
    CHARACTER(LEN=*), INTENT(in)    :: param_name         !< name of the parameter for writing the new value to file
    LOGICAL,          INTENT(in)    :: proportional       !< logical: calulate new value using proportional change
    REAL(wp),         INTENT(inout) :: def_val            !< default parameter value as hardcoded in the mo_PROC_constants

    ! ---------------------------
    ! 0.2 Local
    REAL(wp)             :: param_value_diff    ! for negative def_val the absolute change-value is calculated
    REAL(wp)             :: def_val_save        ! safe the original value for file output
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Set_quincy_parameter_value'


    !> 0.9 safe the original value for file output
    !!
    def_val_save = def_val


    !> 1.0 re-calc def_val if namelist value differs from missval (i.e., a new namelist value is set or a proportional-change value is set) \n
    !!     using the nml_prop_chng variable here, because the nml_abs_val might be changed because it has always "potential other parameters already included"
    !!
    IF (ABS(nml_prop_chng - missval) > eps8) THEN
      ! differ between the two cases: just use the new value / change the (hardcoded) default value by x%
      IF (proportional) THEN
        ! special cases for some parameters, such as, temperatures given in K but modified based on the degC value
        SELECT CASE(TRIM(param_name))
        CASE('t_ref_decomposition','t_opt_nitrification','lctlib_t_air_senescence')
            ! values assumed to be in K but changed based on their degC value
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

    !> 2.0 write default and new value to file
    !!
    WRITE(lun_param_out,'(A,2F24.12)') TRIM(param_name), REAL(def_val_save, wp), REAL(def_val, wp)

  END SUBROUTINE Set_quincy_parameter_value

  !-----------------------------------------------------------------------------------------------------
  !> modify QUINCY namelist parameter by multiplication with a factor
  !!
  !!  this factor is associated with the particular namelist parameter and may range [eps4, 10000]
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Modify_quincy_namelist_parameter(nml_param_factor, param_name, nml_param)

    USE mo_kind,                ONLY: wp
    USE mo_jsb_math_constants,  ONLY: eps8, eps4
    USE mo_jsb_io,              ONLY: missval  ! = -9.E+15_wp

    IMPLICIT NONE

    REAL(wp),         INTENT(in)    :: nml_param_factor   !< factor (range [eps4, 10000]) to multiply with the namelist parameter
    CHARACTER(LEN=*), INTENT(in)    :: param_name         !< name of the namelist parameter for writing the new value to file
    REAL(wp),         INTENT(inout) :: nml_param          !< original value of the namelist parameter as provided by the namelist 'qs.namelist'

    REAL(wp)                    :: nml_param_save         ! save the original value for output
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':Modify_quincy_namelist_parameter'

    !> 0.9 safe the original value for file output
    !!
    nml_param_save = nml_param

    !> 1.0 multiply namelist parameter with factor
    !!
    ! if this namelist parameter should be changed (then the factor may differs from missval)
    IF (ABS(nml_param_factor - missval) > eps8) THEN
      !> 1.1 test if factor is within range
      IF (nml_param_factor < eps4 .OR. nml_param_factor > 10000.0_wp) THEN
         WRITE(message_text, *) 'Factor is out of range [eps4, 10000] for multiplying the namelist parameter', param_name
         CALL finish(routine, message_text)
      ENDIF
      !> 1.2 modify namelist parameter
      nml_param = nml_param_factor * nml_param
    ENDIF

    !> 2.0 write default and new value to file
    !!
    WRITE(lun_param_out,'(A,2F24.12)') TRIM(param_name), REAL(nml_param_save, wp), REAL(nml_param, wp)

  END SUBROUTINE Modify_quincy_namelist_parameter

END MODULE mo_qs_set_parameters
