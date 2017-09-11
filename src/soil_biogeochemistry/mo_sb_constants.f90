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

  ! turnover times (mostly for simple model) 
  REAL(wp), SAVE :: & 
      tau_soluable_litter, &    !< turnover time of soluable litter [s]
      tau_metabolic_litter, &   !< turnover time of metabolic litter pool [s]
      tau_structural_litter, &  !< turnover time of structural litter pool [s]
      tau_woody_litter, &       !< turnover time of woody litter [s]
      tau_fast, &               !< turnover time of fast SOM pool [s] 
      tau_slow, &               !< turnover time of slow SOM pool [s] 
      tau_biomineralisation     !< turnover time of slow SOM pool for biomineralisation [s]

  ! temperature and moisture rate modifiers
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
      k_diff_org, &                !< diffusion coefficient for organic material due to bioturbation? [m2 (kg m-3) s-1]                       
      rho_bulk_org, &              !< bulk density of organic material [kg m-3]
      carbon_per_dryweight_SOM, &  !< C per unit dryweight SOM [gC / gOM] 
      sorp_cap_org, &              !< sorption capacity of organic matter [kg C kg-1 OM]
      qmax_org_mineral_default, &  !< sorption capacity of mineral soil [mol C kg-1 mineral soil]
      k_adsorpt_dom, & 
      k_desorpt_dom, &
      k_adsorpt_det, &
      k_desorpt_det, &
      k_adsorpt_nh4, &
      k_desorpt_nh4, &
      k_desorpt_po4, &            ! < Slow (de)sorption rate (between sorbed Pi and labile Pi) of PO4, [mol kg-1 soil s-1]
      k_adsorpt_f_po4, &          ! < fast adsorption rate (from soluble Pi to labile Pi) of PO4, [[kg soil mol-1 P s-1]]
      k_desorp_f_po4, &           ! < fast desorption rate (from labile Pi to soluble Pi) of PO4, [s-1]
      k_adsorpt_s_po4, &          ! < slow adsorption rate (from soluble Pi to sorbed Pi) of PO4, [[kg soil mol-1 P s-1]]
      k_desorp_s_po4              ! < slow desorption rate (from sorbed Pi to soluble Pi) of PO4, [s-1]

  ! microbial parameters
   REAL(wp), SAVE :: &
      microbial_cn, &            !< molar ratio of C:N in microbes (testing only!)
      microbial_np, &            !< molar ratio of N:P in microbes (testing only!)
      microbial_nue, &           !< microbial nitrogen use efficiency
      microbial_pue, &           !< microibal phosphorus use efficiency
      microbial_cue, &           !< microbial carbon-use efficiency
      microbial_cue_max, &       !< maximum microbial carbon-use efficiency
      microbial_cue_min, &       !< minimum microbial carbon-use efficiency
      residue_cue                !< carbon-use efficiency of residue formation

  ! additional microbial parameters for JSM implementation
  REAL(wp), SAVE :: & 
      k_mort_microbial , &       !< mortality rate of microbes [s-1]
      min_mic_biomass, &         !< minimum microbial biomass density (stopping mortality) [mol C/m3]
      vmax_mic_uptake, &         !< maximum uptake rates of microbes [s-1] 
      km_mic_uptake, &           !< half-saturation DOM density for microbial DOM consumption 
      km_depolymerisation, &     !< half-saturation microbial biomass density for depolymerisation
      vmax_depolymerisation_residue, & !< maximum depolymerisation rate of residue
      vmax_depolymerisation_litter     !< maximum depolymerisation rate of polymeric litter

   ! parameters controlliong litter partitioning
   REAL(wp), SAVE :: & 
      f_soluable_residue, &      !< fraction of microbial residue entering the DOM pool
      f_soluable_leaf, &         !< fraction of leaf litter that becomes DOM
      f_soluable_fine_root, &    !< fraction of fine root litter that becomes DOM
      f_soluable_fruit, &        !< fraction of fruit litter that becomes DOM
      f_soluable_seed_bed, &     !< fraction of seed bed litter that becomes DOM
      fc_soluable_max, &         !< maximum fraction of soluable litter formation
      k_fc_soluable, &           !< slope of soluable fraction with lignin to N ratio
      k_fn_soluable, &           !< proportionality factor controlling C:N of solubable vs. polymeric pool 
      k_fp_soluable, &           !< proportionality factor controlling C:P of solubable vs. polymeric pool
      lc_leaf_max, &             !< maximum lignin content of leaves [1/mol]
      lc_leaf2sla, &             !< slope of lignin to sla relationship [1/m2]
      lc_fine_root, &            !< lignin content of fine roots [1/mol]
      lc_coarse_root, &          !< lignin content of coarse roots [1/mol]
      lc_woody_litter, &         !< lignin content of woody litter [1/mol]
      lc_fruit, &                !< lignin content of fruits [1/mol]
      lc_seed_bed                !< lignin content of seed bed [1/mol]

   ! mycorrhizal parameters
   REAL(wp), SAVE :: & 
      mycorrhizal_cue, &         !< mycorrhizal carbon-use efficiency
      mycorrhiza_cn, &           !< molar ratio of C:N of mycorrhizal fung
      mycorrhiza_np, &           !< molar ratio of N:P of mycorrhizal fungi
      f_maxresp_myc, &           !< maximum fraction of mycorrhizal carbon used for N transformation [1/day]
      f_surface_myc2fineroot     !< ratio of surface to mass relationship between mycorrhizal fungi and fine roots

   ! nutrient uptake kinetics & ECA approach
   ! simple model
   REAL(wp), SAVE :: & 
      km1_up_nh4, &             !< low-affinity parameter for plant and mycorrhizal uptake [m3/mol]
      km1_up_no3, &             !< low-affinity parameter for plant and mycorrhizal uptake [m3/mol]
      km1_up_po4, &             !< low-affinity parameter for plant and mycorrhizal uptake [m3/mol]
      km2_up_nh4, &             !< high-affinity parameter for plant and mycorrhizal uptake [m3/mol]
      km2_up_no3, &             !< high-affinity parameter for plant and mycorrhizal uptake [m3/mol]
      km2_up_po4                !< high-affinity parameter for plant and mycorrhizal uptake [m3/mol]

   ! jsm parameters
   REAL(wp), SAVE :: &
      vmax_mic_uptake_nh4, &     !< maximum NH4 uptake capacity of microbes (mol N / mol biomass C / s)  
      km_mic_uptake_nh4, &       !< half-saturation concentration (mol N / m3)
      vmax_mic_uptake_no3, &     !< maximum NO3 uptake capacity of microbes (mol N / mol biomass C / s)
      km_mic_uptake_no3, &       !< half-saturation concentration (mol N / m3)
      vmax_plant_uptake_nh4, &   !< maximum NH4 uptake capacity of fine roots (mol N / mol biomass C / s)
      km_plant_uptake_nh4, &     ! < half-saturation concentration (mol N / m3)
      vmax_plant_uptake_no3, &   !< maximum NO3 uptake capacity of fine roots (mol N / mol biomass C / s)
      km_plant_uptake_no3, &     !< half-saturation concentration (mol N / m3)
      vmax_mic_uptake_po4, &     !< maximum PO4 uptake capacity of microbes (mol PO4 / mol biomass C / s)
      km_mic_uptake_po4, &       !< half-saturation concentration for decomposing microbe PO4 immobilization [mol PO4 m-3]
      vmax_plant_uptake_po4, &   !< maximum PO4 uptake capacity of fine roots (mol PO4 / mol biomass C / s)
      km_plant_uptake_po4, &     !< half-saturation concentration (mol PO4 / m3)
      vmax_denitirfication, &    !< maximum rate of denitrification per N uptake
      k_mic_uptake_po4, &        !< Reaction rate of microbe PO4 carrier enzyme [s-1], ECA specific
      k_plant_uptake_po4, &      !< Reaction rate of plant PO4 carrier enzyme [s-1], ECA specific
      f_mic_carrier_po4, &       !< ratio between microbial P transport and microbial biomass, []
      f_plant_carrier_po4, &     !< ratio between plant P transport and fine root biomass, []
      qmax_po4_mineral, &        !< Maximum PO4 sorption capacity of mineral soil [mol P kg-1 soil], ECA and Langmuir
      km_adsorpt_po4             !< Half-saturation concentration for PO4 mineral soil adsorption [mol P kg-1 soil], ECA and Langmuir

   ! nitrification, denitrification parameters (sb_nloss_scheme = dynamic)
   REAL(wp), SAVE :: & 
      ea_nitrification, &        !< activation energy for nitrification [J mol-1] 
      ed_nitrification, &        !< deactivation energy for nitrification [J mol-1]
      t_opt_nitrification, &     !< optimum temperature for nitrification [K]
      ea_denitrification, &      !< activation energy for denitrification [J mol-1]
      ea_gasdiffusion, &         !< activation energy for denitrification [J mol-1]
      lambda_anvf, &             !< lambda in the Weibull function relating soil moisture to the anaerobic volume fraction
      k_anvf, &                  !< k in the Weibull function relating soil moisture to the anaerobic volume fraction
      vmax_nitrification, &      !< maximum nitrification rate [µmol mol-1 s-1] 
      vmax_denitrification, &    !< maximum denitrification rate [µmol mol-1 s-1]
      km_denitrification_c, &    !< half-saturation constant for carbon-limited denitrification [mol C m-3]
      km_denitrification_no3, &  !< half-saturation constant for nitrate-limited denitrifcation [mol N m-3]
      f_nit_noy, &               !< fraction of nitrification lost as NOy [--]
      f_nit_n2o, &               !< fraction of nitrification lost as N2O [--]
      f_denit_noy, &             !< fraction of denitrification lost as NOy [--]
      f_denit_n2o                !< fraction of denitrification lost as N2O [--]

   ! Simple-model fixed loss and simple leaching parameterisation
   REAL(wp), SAVE :: &
      floss_nmin, &              !< fractional loss of N per unit N mineralisation 
      fnit_nmin, &               !< fractional nitrification per unit N mineralisation 
      fleach_nh4, &              !< fraction of NH4 amenable to leaching
      fleach_po4, &              !< fraction of PO4 amenable to leaching
      soil_n_threshold_bnf       !< maximum soil N at which BNF occurs in standard scheme (mol N / m2) 

   ! Isoptope stuff
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

	! soil Pi fluxes
	REAL(wp), SAVE :: &
      k_weath_mineral, &         !< Weathering coefficient of mineral soil [mol P kg-1 mineral soil s-1]
      eq_sorp_mineral_po4, &     !< Slow Pi sorption equilibrium content of mineral soil [mol P kg-1 mineral soil]
      k_occlude_po4, &           !< Occlusion coefficient of sorbed PO4, [mol P kg-1 mineral soil s-1]
      vmax_biochem_po4, &        !< Maximum PO4 biochemical mineralization rate [s-1], MM kinetics (potentially also for ECA)
      km_mic_biochem_po4         !< Half-saturation microbial biomass for PO4 biochemical mineralization [mol C m-3], MM kinetics (potentially also for ECA)

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

    ! turnover times
    tau_soluable_litter   = 30._wp * one_day ! 30 days 
    tau_metabolic_litter  = 0.066_wp * one_year * one_day
    tau_structural_litter = 0.245_wp * one_year * one_day
    tau_woody_litter      = 10._wp * one_year * one_day ! 10 years 
    tau_fast              = 5._wp * one_year * one_day
    tau_slow              = 50._wp * one_year * one_day
    tau_biomineralisation = 5._wp * one_year * one_day

    ! temperature and moisture rate modifiers 
    t_ref_tresponse     = 293.15_wp   ! Wang, G., et al. 2012., Soil Biology and Biochemistry 48, 28-38.
    ea_depolymerisation = 53000._wp   ! B. Ahrens, pers. comm. 2016
    ea_mic_uptake       = 47000._wp   ! B. Ahrens, pers. comm. 2016
    ea_sorption         = 5000._wp    ! B. Ahrens, pers. comm. 2016
    ea_desorption       = 20000._wp   ! B. Ahrens, pers. comm. 2016
    ea_hsc              = 30000._wp   ! B. Ahrens, pers. comm. 2016
    k1_afps             = 4._wp / 3._wp ! B. Ahrens, pers. comm. 2016
    kmr_depolymerisation= 0.001_wp
    kmr_mic_uptake      = 0.001_wp
    kmr_hsc             = 0.001_wp

    ! transport and sorption parameters (all B. Ahrens, pers. comm. 2016)
    k_diff_org          = 0.286933_wp / 8.0_wp / one_year / one_day
    rho_bulk_org        = 150.3935_wp 
    carbon_per_dryweight_SOM = 0.4271_wp   !! COMISSION for Heinich
    sorp_cap_org        = 0.0_wp
    qmax_org_mineral_default = 0.055_wp / molar_mass_C * 1000._wp 
    k_adsorpt_dom       = 843.391_wp / one_year / one_day * molar_mass_C / 1000._wp ! conversion to m3 / mol C / s
    k_desorpt_dom       = 4301.741_wp / one_year / one_day
    k_adsorpt_det       = 6733.843_wp / one_year / one_day * molar_mass_C / 1000._wp ! conversion to m3 / mol C / s
    k_desorpt_det       = 2065.491_wp / one_year / one_day
    k_adsorpt_nh4       = 0.0_wp
    k_desorpt_nh4       = 0.0_wp
    k_adsorpt_f_po4     = 2000._wp*57.9_wp/3600._wp   ! [kg soil mol-1 P s-1], based on Van der Zee et al. 1989
    k_desorpt_po4       = 0.0_wp                ! [mol kg-1 soil s-1], estimated from McGechan and Lewis 2002
    k_desorp_f_po4      = 0.14_wp / 3600._wp ! [s-1], Van der Zee et al. 1989
    k_desorp_s_po4      = 0.0_wp                ! [kg soil mol-1 P s-1], McGechan and Lewis 2002
    k_adsorpt_s_po4     = 0.0_wp               ! [s-1], McGechan and Lewis 2002

    ! microbial parameters
    microbial_cn        = 9._wp * molar_mass_N / molar_mass_C    ! textbook knowledge
    microbial_np        = 14._wp * molar_mass_P / molar_mass_N   ! textbook knowledge
    microbial_nue       = 0.8_wp                                 ! Manzoni et al. 2008
    microbial_pue       = 0.8_wp                                 ! Manzoni et al. 2008
    microbial_cue       = 0.392_wp                               ! COMISSION paper
    microbial_cue_max   = 0.6_wp                                 ! Manzoni et al. 2008
    microbial_cue_min   = 0.3_wp                                 ! Manzoni et al. 2008
    residue_cue         = 0.1_wp                                 ! made-up
   
    ! additional microbial parameters for JSM implementation
    k_mort_microbial    = 7.48_wp * 1E-8
    min_mic_biomass     = 0.01_wp
    vmax_mic_uptake     = 3495.18053_wp / one_year / one_day * 1000000._wp
    km_mic_uptake       = 1.02404788_wp / molar_mass_C * 1000_wp ! conversion to mol C / m3
    vmax_depolymerisation_residue = 0.2316605_wp / one_year / one_day
    vmax_depolymerisation_litter = 0.369800811690762_wp / one_year / one_day 
    km_depolymerisation = 0.08876_wp / molar_mass_C * 1000_wp ! conversion to mol C / m3

    ! parameters controlliong litter partitioning
    f_soluable_residue  = 0.172_wp
    f_soluable_leaf      = 0.5_wp
    f_soluable_fine_root = 0.5_wp
    f_soluable_fruit     = 0.5_wp
    f_soluable_seed_bed  = 0.5_wp
    fc_soluable_max      = 0.85_wp  ! Century
    k_fc_soluable        = 0.018_wp ! Century
    k_fn_soluable        = 5.0_wp   ! Century
    k_fp_soluable        = 5.0_wp   ! Century
    lc_leaf_max          = 0.295_wp / molar_mass_C * molar_mass_N
    lc_leaf2sla          = -0.3712_wp / molar_mass_C * molar_mass_N
    lc_fine_root         = 0.22_wp / molar_mass_C * molar_mass_N
    lc_coarse_root       = 0.7_wp / molar_mass_C * molar_mass_N
    lc_woody_litter      = 0.7_wp / molar_mass_C * molar_mass_N
    lc_fruit             = 0.22_wp / molar_mass_C * molar_mass_N 
    lc_seed_bed          = 0.22_wp / molar_mass_C * molar_mass_N

    ! mycorrhizal parameters
    mycorrhizal_cue     = 0.7_wp !! ??? SOURCE
    mycorrhiza_cn       = 10.0_wp * molar_mass_N / molar_mass_C
    mycorrhiza_np       = 14.0_wp * molar_mass_P / molar_mass_N
    f_maxresp_myc       = 0.4_wp !! ??? SOURCE
    f_surface_myc2fineroot = 10.0_wp !! ??? SOURCE

    ! nutrient uptake kinetics
    ! simple model
    km1_up_nh4 = 0.0416_wp
    km2_up_nh4 = 1.0_wp 
    km1_up_no3 = 0.0416_wp
    km2_up_no3 = 1.0_wp
    km1_up_po4 = 0.0416_wp
    km2_up_po4 = 1.0_wp

    ! additional jsm parameters 
    vmax_mic_uptake_nh4   = 48.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein 
    km_mic_uptake_nh4     = 1.2_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein    
    vmax_mic_uptake_no3   = 39.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein 
    km_mic_uptake_no3     = 1.8_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein   
    vmax_plant_uptake_nh4 = 9.7_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein
    km_plant_uptake_nh4   = 0.18_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein 
    vmax_plant_uptake_no3 = 8.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s from Zhu et al. 2017, and references therein
    km_plant_uptake_no3   = 0.41_wp / molar_mass_N ! mol / m3 from Zhu et al. 2017, and references therein
    vmax_mic_uptake_po4   = 48.0_wp / 3600._wp * molar_mass_C ! µmol / mol / s
    km_mic_uptake_po4     = 0.02_wp / molar_mass_P ! mol / m3
    vmax_plant_uptake_po4 = 9.7_wp / 3600._wp * molar_mass_C ! µmol / mol / s
    km_plant_uptake_po4   = 0.067_wp / molar_mass_P ! mol / m3
    vmax_denitirfication  = 0.01_wp
	k_mic_uptake_po4      = 400._wp / 10000._wp/ one_day ! [s-1], Zhu et al. 2016, and references therein
	k_plant_uptake_po4    = 12._wp / 100._wp/ one_day ! [s-1], Zhu et al. 2016, and references therein/ 100._wp
	f_mic_carrier_po4     = 0.0015_wp ! assuming 30% of the microbial enzymes are P transporter, based on Zhu et al. 2016
	f_plant_carrier_po4   = 0.1_wp / 80._wp ! assuming 10% of the plant enzymes are P transporter, based on Zhu et al. 2016
	qmax_po4_mineral      = 0.134_wp / molar_mass_P ! [mol P kg-1 soil], Olander and Vitousek 2005
	km_adsorpt_po4        = 0.075_wp / molar_mass_P ! [mol P kg-1 soil], Wang et al. 2010


	k_weath_mineral       = 0.025_wp / molar_mass_P / one_year / one_day ! [mol P m-3 s-1], Wang et al. 2010
	eq_sorp_mineral_po4   = 80._wp / molar_mass_P ![mol P m-3] ! random guess based on init value of assoc_slow
	k_occlude_po4         = 1E-6_wp / 30._wp / one_day ! s-1, Yang et al. 2014
	vmax_biochem_po4      = 1.0_wp / 40.0_wp / one_day ![s-1], estimated from Bunemann et al. 2016, LUS site
	km_mic_biochem_po4    = 0.1_wp / molar_mass_C ![mol C m-3] ! random number, equals to the usual microbe biomass (0.1 g C m-3)

    ! Nitrification, Denitrification parameters
    lambda_anvf            = 1.8_wp         ! fitted from Zaehle et al. 2011, NatGeo 
    k_anvf                 = 3.0_wp         ! fitted from Zaehle et al. 2011, NatGeo
    ea_nitrification       = 80000._wp      ! Xu-Ri & Prentice, 2008, to fit their function
    ed_nitrification       = 200000._wp     ! Xu-Ri & Prentice, 2008, to fit their function
    t_opt_nitrification    = 331.15_wp      ! Xu-Ri & Prentice, 2008, to fit their function 
    f_nit_noy              = 0.02           ! Xu-Ri & Prentice, 2008
    f_nit_n2o              = 0.002          ! Xu-Ri & Prentice, 2008
    ea_denitrification     = 47000._wp      ! Xu-Ri & Prentice, 2008, to fit their function
    vmax_nitrification     = 0.1_wp / one_day * 1.e6             ! Xu-Ri & Prentice, 2008, 10% / day
    vmax_denitrification   = 0.1_wp / one_day * 1.e6          ! Xu-Ri & Prentice, 2008, 10% / day
    km_denitrification_c   = 0.017 * 1000.0_wp / molar_mass_C ! Xu-Ri & Prentice, 2008   
    km_denitrification_no3 = 0.083 * 1000.0_wp / molar_mass_N ! Xu-Ri & Prentice, 2008
    f_denit_noy            = 0.002          ! Xu-Ri & Prentice, 2008
    f_denit_n2o            = 0.02           ! Xu-Ri & Prentice, 2008
    ea_gasdiffusion        = 47000._wp      ! Xu-Ri & Prentice, 2008, to fit their function

    ! parameters for simple soil biogeochemical model 
    floss_nmin = 0.05_wp
    fnit_nmin  = 0.1_wp
    ! simple leaching
    fleach_nh4 = 0.1_wp
    fleach_po4 = 0.001_wp
    soil_n_threshold_bnf = 0.04 ! mol N == 0.56 g N /m3

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
