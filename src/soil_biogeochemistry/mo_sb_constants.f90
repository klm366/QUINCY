!>
!! @brief Contains constants for the soil_biogeochemistry process
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170124)
!!
MODULE mo_sb_constants

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PUBLIC

  ! parameters for simple soil biogeochemical model
  REAL(wp), SAVE :: & 
      floss_nmin, &             !< fractional loss of N per unit N mineralisation 
      fnit_nmin, &              !< fractional nitrification per unit N mineralisation 
      tau_slow, &               !< ..
      tau_fast, &               !< ..
      tau_non_woody_litter, &   !< ..
      residue_cue               !< ..

  ! parameter for JSM implementation

  ! rate modifiers
  REAL(wp), SAVE :: &
      t_ref_tresponse, &      !< reference temperature for temperature response function [K]
      ea_depolymerisation, &  !< activation energy for depolymerisation [J mol-1]
      ea_mic_uptake, &        !< activation energy for microbial uptake [J mol-1]
      ea_sorption, &          !< activation energy for sorption to mineral surfaces [J mol-1]
      ea_desorption, &        !< activation energy for desorption to mineral surfaces [J mol-1]
      ea_hsc, &               !< activation energy for half-saturation constants [J mol-1]
      k1_afps, &              !< scaling factor in the air-filled pore space to modifier response
      kmr_depolymerisation, & !< scaling factor for the sensitivity of depolymerisation to oxygen limitation 
      kmr_mic_uptake, &       !< scaling factor for the sensitivity of microbial uptake to oxygen limitation 
      kmr_hsc                 !< scaling factor for the sensitivity of half-saturation constant to moisture limitation 

  ! transport and sorption parameters
  REAL(wp), SAVE :: &
      k_diff_org, &         !< diffusion coefficient for organic material due to bioturbation? [m2 (kg m-3) s-1]                       
      rho_bulk_org, &       !< bulk density of organic material [kg m-3]
      sorp_cap_org, &       !< sorption capacity of organic matter [kg C kg-1 OM]
      qmax_org_mineral_default, & !< sorption capacity of mineral soil [mol C kg-1 mineral soil]
      k_adsorpt_dom, & 
      k_desorpt_dom, &
      k_adsorpt_det, &
      k_desorpt_det, &
      k_adsorpt_nh4, &
      k_desorpt_nh4, &
      k_adsorpt_po4, &
      k_desorpt_po4

  ! microbial parameters
  REAL(wp), SAVE :: & 
      k_mort_microbial , &       !< mortality rate of microbes [s-1]
      f_soluable_residue, &      !< fraction of microbial residue entering the DOM pool
      min_mic_biomass, &         !< minimum microbial biomass density (stopping mortality) [mol C/m3]
      vmax_mic_uptake, &         !< maximum uptake rates of microbes [s-1] 
      km_mic_uptake, &           !< half-saturation DOM density for microbial DOM consumption 
      km_depolymerisation, &     !< half-saturation microbial biomass density for depolymerisation
      tau_soluable_litter, &     !< turnover time of soluable litter
      tau_woody_litter, &        !< turnover time of woody litter
      vmax_depolymerisation_residue, & !< maximum depolymerisation rate of residue
      vmax_depolymerisation_litter, &  !< maximum depolymerisation rate of polymeric litter
      microbial_cue              !< microbial carbon-use efficiency

   REAL(wp), SAVE :: & 
      f_soluable_leaf, &         !< fraction of leaf litter that becomes DOM
      f_soluable_fine_root, &    !< fraction of fine root litter that becomes DOM
      f_soluable_fruit, &        !< fraction of fruit litter that becomes DOM
      f_soluable_seed_bed        !< fraction of seed bed litter that becomes DOM

   REAL(wp), SAVE :: & 
     carbon_per_dryweight_SOM    !< C per unit dryweight SOM [gC / gOM] 

   REAL(wp), SAVE :: &
      microbial_cn, &            !< molar ratio of C:N in microbes (testing only!)
      microbial_np, &            !< molar ratio of N:P in microbes (testing only!)
      microbial_nue, & 
      microbial_pue

   ! nutrient uptake kinetics
   REAL(wp), SAVE :: &
      vmax_mic_uptake_nh4, &     !< maximum NH4 uptake capacity of microbes (mol N / mol biomass C / s)  
      km_mic_uptake_nh4, &       !< half-saturation concentration (mol N / m3)
      vmax_mic_uptake_no3, &     !< maximum NO3 uptake capacity of microbes (mol N / mol biomass C / s)
      km_mic_uptake_no3, &       !< half-saturation concentration (mol N / m3)
      vmax_plant_uptake_nh4, &   !< maximum NH4 uptake capacity of fine roots (mol N / mol biomass C / s)
      km_plant_uptake_nh4, &     ! < half-saturation concentration (mol N / m3)
      vmax_plant_uptake_no3, &   !< maximum NO3 uptake capacity of fine roots (mol N / mol biomass C / s)
      km_plant_uptake_no3, &     !< half-saturation concentration (mol N / m3)
      vmax_mic_uptake_po4, &     !< maximum PO4 uptake capacity of microbes (mol N / mol biomass C / s)  
      km_mic_uptake_po4, &       !< half-saturation concentration (mol N / m3)
      vmax_plant_uptake_po4, &   !< maximum PO4 uptake capacity of fine roots (mol N / mol biomass C / s)
      km_plant_uptake_po4, &     ! < half-saturation concentration (mol N / m3)
      vmax_nitrification, &      !< maximum rate of nitrification per N uptake
      vmax_denitirfication       !< maximum rate of denitrification per N uptake

   ! isoptop stuff
   REAL(wp), SAVE :: & 
      eta_mic_uptake_nh4, &      !< fractionation of N15 due to microbial NH4 assimilation (per mill)
      eta_mic_uptake_no3, &      !< fractionation of N15 due to microbial NO3 assimilation (per mill)
      eta_plant_uptake_nh4, &    !< fractionation of N15 due to plant NH4 assimilation (per mill)
      eta_plant_uptake_no3, &    !< fractionation of N15 due to plant NO3 assimilation (per mill)
      eta_nitrification, &       !< fractionation of N15 due to microbial nitrification to NO and N2O (per mill)
      eta_nitrate_production, &  !< fractionation of N15 due to microbial NO3 production during nitrification (per mill)
      eta_denitrification, &     !< fractionation of N15 due to microbial denitrification (per mill)
      eta_ammonification, &      !< fractionation of N15 due to microbial ammonification (per mill)
      eta_volatilisation_nh3     !< fractionation of N15 due to NH3 volatilisation (per mill)

   ! only temporarily
   REAL(wp), SAVE :: & 
      nh4_solute_prescribe, no3_solute_prescribe, po4_solute_prescribe

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil_biogeochemistry
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_sb_constants

    USE mo_base_constants, ONLY: one_day, one_year
    USE mo_shared_constants, ONLY: molar_mass_C, molar_mass_N, molar_mass_P
 
    IMPLICIT NONE
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sb_constants'

    ! parameters for simple soil biogeochemical model 
    floss_nmin = 0.05_wp
    fnit_nmin  = 0.1_wp
    tau_slow = 50._wp * one_year * one_day
    tau_fast = 5._wp * one_year * one_day
    tau_non_woody_litter = 0.7_wp * one_year * one_day
    residue_cue = 0.1_wp

    ! parameter for JSM implementation
    t_ref_tresponse = 293.15_wp   ! Wang, G., et al. 2012., Soil Biology and Biochemistry 48, 28-38.
    ea_depolymerisation = 53000._wp   ! B. Ahrens, pers. comm. 2016
    ea_mic_uptake       = 47000._wp   ! B. Ahrens, pers. comm. 2016
    ea_sorption         = 5000._wp   ! B. Ahrens, pers. comm. 2016
    ea_desorption       = 20000._wp   ! B. Ahrens, pers. comm. 2016
    ea_hsc              = 30000._wp   ! B. Ahrens, pers. comm. 2016
    k1_afps             = 4._wp / 3._wp ! B. Ahrens, pers. comm. 2016
    kmr_depolymerisation= 0.001_wp
    kmr_mic_uptake      = 0.001_wp
    kmr_hsc             = 0.0_wp

    k_diff_org          = 0.286933_wp / 8.0_wp / one_year / one_day
    rho_bulk_org        = 150.3935_wp 
    sorp_cap_org        = 0.0_wp
    qmax_org_mineral_default = 0.055_wp / molar_mass_C * 1000._wp 
    k_adsorpt_dom       = 843.391_wp / one_year / one_day * molar_mass_C / 1000._wp ! conversion to m3 / mol C / s
    k_desorpt_dom       = 4301.741_wp / one_year / one_day
    k_adsorpt_det       = 6733.843_wp / one_year / one_day * molar_mass_C / 1000._wp ! conversion to m3 / mol C / s
    k_desorpt_det       = 2065.491_wp / one_year / one_day
    k_adsorpt_nh4       = 0.0_wp
    k_desorpt_nh4       = 0.0_wp
    k_adsorpt_po4       = 0.0_wp
    k_desorpt_po4       = 0.0_wp

    k_mort_microbial    = 7.48_wp * 1E-8
    f_soluable_residue  = 0.172_wp
    min_mic_biomass     = 0.01_wp
    vmax_mic_uptake     = 3495.18053_wp / one_year / one_day * 1000000._wp
    km_mic_uptake       = 1.02404788_wp / molar_mass_C * 1000_wp ! conversion to mol C / m3
    tau_soluable_litter = 30._wp * one_day ! 30 days 
    tau_woody_litter    = 20._wp * one_year * one_day ! 20 years 
    vmax_depolymerisation_residue = 0.2316605_wp / one_year / one_day
    vmax_depolymerisation_litter = 0.369800811690762_wp / one_year / one_day 
    km_depolymerisation = 0.08876_wp / molar_mass_C * 1000_wp ! conversion to mol C / m3
    microbial_cue       = 0.392_wp !! COMISSION paper
   
    f_soluable_leaf      = 0.5_wp
    f_soluable_fine_root = 0.5_wp
    f_soluable_fruit     = 0.5_wp
    f_soluable_seed_bed  = 0.5_wp

    carbon_per_dryweight_SOM = 0.4271_wp   !! COMISSION for Heinich

    microbial_cn = 9._wp * molar_mass_N / molar_mass_C
    microbial_np = 14._wp * molar_mass_P / molar_mass_N
    microbial_nue = 0.8_wp
    microbial_pue = 0.8_wp
 
    vmax_mic_uptake_nh4   = 48.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein 
    km_mic_uptake_nh4     = 1.2_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein    
    vmax_mic_uptake_no3   = 39.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein 
    km_mic_uptake_no3     = 1.8_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein   
    vmax_plant_uptake_nh4 = 9.7_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein
    km_plant_uptake_nh4   = 0.18_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein 
    vmax_plant_uptake_no3 = 8.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein
    km_plant_uptake_no3   = 0.41_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein   
    vmax_mic_uptake_po4   = 48.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s 
    km_mic_uptake_po4     = 1.2_wp / molar_mass_N ! mol / m3 
    vmax_plant_uptake_po4 = 9.7_wp / 3600._wp * molar_mass_C ! µmol / mol / s 
    km_plant_uptake_po4   = 0.18_wp / molar_mass_N ! mol / m3 
    vmax_nitrification    = 0.01_wp    
    vmax_denitirfication  = 0.01_wp  

    ! All values from Robinson 2001, TREE
    eta_mic_uptake_nh4     = 17._wp
    eta_mic_uptake_no3     = 13._wp
    eta_plant_uptake_nh4   = 13.5_wp
    eta_plant_uptake_no3   = 9.5_wp
    eta_nitrification      = 47.5_wp
    eta_nitrate_production = 25._wp
    eta_denitrification    = 31._wp
    eta_ammonification     = 2.5_wp
    eta_volatilisation_nh3 = 50._wp

    ! BLARPP temporarily
    nh4_solute_prescribe = 0.1_wp
    no3_solute_prescribe = 0.1_wp
    po4_solute_prescribe = 0.1_wp
 
  END SUBROUTINE init_sb_constants

END MODULE mo_sb_constants
