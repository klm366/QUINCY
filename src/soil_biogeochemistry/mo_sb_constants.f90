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
      tau_litter, &        !< soil litter pool turnover time (years)
      floss_nmin, &        !< fractional loss of N per unit N mineralisation 
      fnit_nmin            !< fractional nitrification per unit N mineralisation 

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
      kmr_sorption, &         !< scaling factor for the sensitivity of sorption to moisture limitation 
      kmr_desorption, &       !< scaling factor for the sensitivity of desorption to moisture limitation 
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
      microbial_np               !< molar ratio of N:P in microbes (testing only!)

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants from soil_biogeochemistry
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_sb_constants

    USE mo_base_constants, ONLY: one_day, one_year
    USE mo_shared_constants, ONLY: molar_mass_C 
 
    IMPLICIT NONE
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_sb_constants'

    ! parameters for simple soil biogeochemical model 
    tau_litter = 5.0_wp  ! years
    floss_nmin = 0.05_wp
    fnit_nmin  = 0.1_wp

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
    kmr_sorption        = 0.50853_wp
    kmr_desorption      = 0.50853_wp
    kmr_hsc             = 0.0_wp

    k_diff_org          = 0.286933_wp / one_year / one_day
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
    vmax_mic_uptake     = 3495.18053_wp / one_year / one_day
    km_mic_uptake       = 1.02404788_wp / molar_mass_C * 1000_wp ! conversion to mol C / m3
    tau_soluable_litter = 10._wp * one_day ! 10 days 
    tau_woody_litter    = 20._wp * one_year * one_day ! 20 years 
    vmax_depolymerisation_residue = 0.2316605_wp / one_year / one_day
    vmax_depolymerisation_litter = 0.369800811690762_wp / one_year / one_day 
    km_depolymerisation = 0.08876_wp / molar_mass_C * 1000_wp ! conversion to mol C / m3
    microbial_cue       = 0.392_wp !! COMISSION paper
   
    f_soluable_leaf    = 0.5_wp
    f_soluable_fine_root    = 0.5_wp
    f_soluable_fruit    = 0.5_wp
    f_soluable_seed_bed    = 0.5_wp

    carbon_per_dryweight_SOM = 0.4271_wp   !! COMISSION for Heinich

    microbial_cn = 6._wp
    microbial_np = 14._wp

  END SUBROUTINE init_sb_constants

END MODULE mo_sb_constants
