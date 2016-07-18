!>
!! @brief Contains the TYPE that holds all parameters and information on the landcover data considering lct "vegetated"
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (2016-06-13)
!!
MODULE mo_base_slm_parameter_class

  USE mo_kind
  USE mo_base_constants, ONLY: npft

  IMPLICIT NONE

  !-----------------------------------------------------------------------------------------------------
  !> initial version of the TYPE that holds all parameters of SLM
  !!
  !! global parameters to describe the vegetation processes \n 
  !! each variable is a vector over all plant functional types (pft with number npft)
  !-----------------------------------------------------------------------------------------------------
  TYPE t_slm_parameter
  !renamed using: TYPE t_slm_veg_param
  ! Respiration parameters
  REAL(wp),DIMENSION(npft) :: fresp_growth                 !< growth respiration per unit new biomass 
  REAL(wp),DIMENSION(npft) :: fmaint_rate_base             !< N-specific maintenance respiration rate (µmol CO2 / mmol N / s) 
  REAL(wp),DIMENSION(npft) :: t_maint_k1                   !< Coefficient in the equation for temperature sensitivity of foliar respiration 
  REAL(wp),DIMENSION(npft) :: t_maint_k2                   !< Coefficient in the equation for temperature sensitivity of foliar respiration
  REAL(wp),DIMENSION(npft) :: t_maint_k3                   !< Coefficient in the equation for temperature sensitivity of foliar respiration 

  ! Labile pool parameters 
  REAL(wp),DIMENSION(npft) :: tau_labile_default           !< default turnover time of the labile carbo-hydrates (in days)
  REAL(wp),DIMENSION(npft) :: lambda_tau_labile_temp       !< lambda in temperature response function of labile pool 
  REAL(wp),DIMENSION(npft) :: weibull_k_tau_labile_temp    !< k in temperature response function of labile pool 
  REAL(wp),DIMENSION(npft) :: lambda_tau_labile_theta      !< lambda in moisture response function of labile pool 
  REAL(wp),DIMENSION(npft) :: weibull_k_tau_labile_theta   !< k in moisture response function of labile pool 
  REAL(wp),DIMENSION(npft) :: lambda_phiphen               !< lambda in storage response function to phenology 
  REAL(wp),DIMENSION(npft) :: weibull_k_phiphen            !< k in in storage response function to phenology 
  REAL(wp),DIMENSION(npft) :: lambda_phimaint              !< lambda in storage response function to maintenance needs 
  REAL(wp),DIMENSION(npft) :: weibull_k_phimaint           !< k in in storage response function to maintenance needs
  REAL(wp),DIMENSION(npft) :: lambda_phistore              !< lambda in storage response function to storage size
  REAL(wp),DIMENSION(npft) :: weibull_k_phistore           !< k in in storage response function to storage size 

  ! Resporption parameters
  REAL(wp),DIMENSION(npft) :: resorption_default           !< default fraction of nutrient resporption before leaf shedding
  REAL(wp),DIMENSION(npft) :: resorption_wood              !< default fraction of nutrient resporption before wood death
  REAL(wp),DIMENSION(npft) :: root2leaf_cn                 !< relative C:N of fine roots compared to leaves
  REAL(wp),DIMENSION(npft) :: wood2leaf_cn                 !< relative C:N of woody biomass compared to leaves
  REAL(wp),DIMENSION(npft) :: root2leaf_np                 !< relative N of fine roots compared to leaves
  REAL(wp),DIMENSION(npft) :: wood2leaf_np                 !< relative N:P of woody biomass compared to leaves

  ! Parameters describing photosynthetic activity 
  REAL(wp),DIMENSION(npft) :: jmax2n                       !< electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
  REAL(wp),DIMENSION(npft) :: vcmax2n                      !< Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
  REAL(wp),DIMENSION(npft) :: k0_fnl_pepc                  !< fraction of N in PEP and PPKD (C4 plants only) 
  REAL(wp),DIMENSION(npft) :: pepc2n                       !< PEP C limited carboxylation rate per µmol CO2 / mmol N in PEP C
  REAL(wp),DIMENSION(npft) :: chl2n                        !< moles chlorophyll per mol of N
  REAL(wp),DIMENSION(npft) :: k1_fnl_struc                 !< slope of the fraction of leaf N not associated with PS against leaf N (--/mmolN) 
  REAL(wp),DIMENSION(npft) :: jmax2vcmax_C3                !< ratio of Jmax25/Vcmax25 for prescribed leaf N (C3 plants)
  REAL(wp),DIMENSION(npft) :: jmax2vcmax_C4                !< ratio of Jmax25/Vcmax25 for prescribed leaf N (C4 plants)
  REAL(wp),DIMENSION(npft) :: tau_canopy_acclimation       !< time scale of acclimation of canopy processes (N, P, T) (days)
  ! temperature sensitivies of PS processes
  REAL(wp),DIMENSION(npft) :: E0kc                         !< Baseline of Michaelis-Menten constant of Rubisco for CO2 (unit?)
  REAL(wp),DIMENSION(npft) :: E1kc                         !< Temperature sensitivity of Michaelis-Menten constant of Rubisco for CO2 (unit?)
  REAL(wp),DIMENSION(npft) :: E0ko                         !< Baseline of Michealis-Menten Kintetics for Rubisco-oxigination
  REAL(wp),DIMENSION(npft) :: E1ko                         !< Temperature sensitivity of Michealis-Menten Kintetics for Rubisco-oxigination
  REAL(wp),DIMENSION(npft) :: E0pcp                        !< Baseline of photosynthetic compensation point
  REAL(wp),DIMENSION(npft) :: E1pcp                        !< Temperature sensitivity of photosynthetic compensation point
  REAL(wp),DIMENSION(npft) :: E0v                          !< Baseline of Rubisco
  REAL(wp),DIMENSION(npft) :: E1v                          !< Temperature sensitivity of Rubisco
  REAL(wp),DIMENSION(npft) :: Tref_pepc                    !< Temperature sensitivity of PeP C4 photosynthesis 
  REAL(wp),DIMENSION(npft) :: Tbase_pepc                   !< Temperature sensitivity of PeP C4 photosynthesis
  REAL(wp),DIMENSION(npft) :: alpha_i                      !< intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I))
                                                !< µmol CO2 / mol quanta) 
  REAL(wp),DIMENSION(npft) :: pO2                          !< Partial Pressure of O2 in kPa
  REAL(wp),DIMENSION(npft) :: CiCa_default_C3              !< default Ci:Ca for C3 plants
  REAL(wp),DIMENSION(npft) :: CiCa_default_C4              !< default Ci:Ca for C4 plants
  REAL(wp),DIMENSION(npft) :: ps_it_max                    !< maximum iteration of PS calculation
  REAL(wp),DIMENSION(npft) :: ci_max                       !< saturating Ci in Pa C4 plants
  REAL(wp),DIMENSION(npft) :: Dwv2co2_air                  !< ratio of diffusion coefficient for H2O and CO2 in air 
  REAL(wp),DIMENSION(npft) :: Dwv2co2_turb                 !< ratio of diffusion coefficient for H2O and CO2 in turbulent air 
  ! parameters describing the fractionation due to photosynthesis
  REAL(wp),DIMENSION(npft) :: discr_ps_a_C13               !< discrimination of C13 due to stomatal diffusion
  REAL(wp),DIMENSION(npft) :: discr_ps_b_C13               !< discrimination of C13 due to Rubisco
  REAL(wp),DIMENSION(npft) :: discr_ps_c_C13               !< discrimination of C13 due to PEP C 
  REAL(wp),DIMENSION(npft) :: discr_ps_a_C14               !< discrimination of C13 due to stomatal diffusion
  REAL(wp),DIMENSION(npft) :: discr_ps_b_C14               !< discrimination of C13 due to Rubisco
  REAL(wp),DIMENSION(npft) :: discr_ps_c_C14               !< discrimination of C13 due to PEP C
  REAL(wp),DIMENSION(npft) :: discr_ps_phi                 !< leakage rate of bundle sheath cells
  
  ! parameters for canopy profile and canopy light extinction
  REAL(wp),DIMENSION(npft) :: kn                           !< extinction coefficient to describe decline of N within the canopy 
  REAL(wp),DIMENSION(npft) :: k0_fnl_chl_C3                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp),DIMENSION(npft) :: k0_fnl_chl_C4                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp),DIMENSION(npft) :: k1_fnl_chl_C3                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp),DIMENSION(npft) :: k1_fnl_chl_C4                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp),DIMENSION(npft) :: kfnl_chl                     !< parameter to describe the varying fraction of N in chlorophyll through the canopy

  ! parameters for changes in canopy N fractions
  REAl(wp),DIMENSION(npft) :: delta_n_leaf                 !< maximum fractional change for optimal adjsutment of canopy N allocation 
  REAL(wp),DIMENSION(npft) :: delta_n_fraction             !< maximum fractional change for optimal adjustment of photosynthetic N fraction

  ! parameters used for phenology calculations
  REAL(wp),DIMENSION(npft) :: min_lai                      !< minimum LAI beyond which all leaves are shed
  REAL(wp),DIMENSION(npft) :: target_lai_max               !< maximum LAI target for reserve use calculations 
  REAL(wp),DIMENSION(npft) :: min_leaf_age                 !< minimum leaf age before senescence is permitted (days) 
  REAL(wp),DIMENSION(npft) :: gdd_t_air_threshold          !< temperature threshold for growing degree days 
  REAL(wp),DIMENSION(npft) :: gdd_req_max                  !< maximum GDD requirement
  REAL(wp),DIMENSION(npft) :: k_gdd_dormance               !< exponent factor in GDD - nd_dormance equation 
  REAL(wp),DIMENSION(npft) :: t_air_senescence             !< weekly air temperature initiating leaf senescence
  REAL(wp),DIMENSION(npft) :: beta_soil_flush_grass        !< value of beta_soil at leaf flushing for grass
  REAL(wp),DIMENSION(npft) :: beta_soil_flush_tree         !< value of beta_soil at leaf flushing for trees 
  REAL(wp),DIMENSION(npft) :: beta_soil_senescence_tree    !< value of beta_soil at leaf shedding for trees

  ! Parameters required for albedo and fapar calculation
  REAL(wp),DIMENSION(npft) :: ka                           !< Extinction coefficient for PAR on chlorophyll
  REAL(wp),DIMENSION(npft) :: sigma                        !< single leaf scattering albedo 
  REAL(wp),DIMENSION(npft) :: kbl0                         !< extinction coefficient over black leaves 
  REAL(wp),DIMENSION(npft) :: kdf0                         !< extinction coefficient for diffuse PAR (Spitters eq. 7) 
  REAL(wp),DIMENSION(npft) :: rho2sbeta                    !< scaling factor of solar angle in reflection calculation (mo_veg_fapar) 

  ! other
  REAL(wp),DIMENSION(npft) :: carbon_per_dryweight_leaf    !< carbon mass per unit dry weight of leaves (gC / g DW) 
  REAL(wp),DIMENSION(npft) :: min_height                   !< minimum permissable height of an individual
  REAL(wp),DIMENSION(npft) :: sm2lm_grass                  !< stem mass to leaf mass ratio of grasses
  REAL(wp),DIMENSION(npft) :: lambda_fruit_alloc           !< parameter in the Weibull function controlling fruit growth
  REAL(wp),DIMENSION(npft) :: weibull_k_fruit_alloc        !< parameter in the Weibull function controlling fruit growth
  REAL(wp),DIMENSION(npft) :: offset_fruit_alloc           !< parameter in the Weibull function controlling fruit growth
  REAL(wp),DIMENSION(npft) :: min_falloc_fruit             !< minimum fraction of allocation going to fruit (flowers etc.)

  ! vegetation dynamics
  REAL(wp),DIMENSION(npft) :: background_mort_rate_tree   !< background mortality rate of trees in static mortality case
  REAL(wp),DIMENSION(npft) :: background_mort_rate_grass  !< background mortality rate of grasses in static mortality case

  ! soil pools
  REAL(wp),DIMENSION(npft) :: tau_litter                  !< soil litter pool turnover time


  ! previously stored in TYPE(t_slm_veg_param),INTENT(inout) :: veg_param

     INTEGER,DIMENSION(npft) :: growthform 
     INTEGER,DIMENSION(npft) :: ps_pathway
     INTEGER,DIMENSION(npft) :: phenology_type
     REAL(wp),DIMENSION(npft) :: tau_leaf
     REAL(wp),DIMENSION(npft) :: sla
     REAL(wp),DIMENSION(npft) :: cn_leaf
     REAL(wp),DIMENSION(npft) :: np_leaf
     REAL(wp),DIMENSION(npft) :: k0_fnl_struc
     REAL(wp),DIMENSION(npft) :: tjmax_opt
     REAL(wp),DIMENSION(npft) :: tjmax_omega
     REAL(wp),DIMENSION(npft) :: g0
     REAL(wp),DIMENSION(npft) :: g1
     REAL(wp),DIMENSION(npft) :: lambda_soil_ps
     REAL(wp),DIMENSION(npft) :: weibull_k_soil_ps
     REAL(wp),DIMENSION(npft) :: lambda_soil_gs
     REAL(wp),DIMENSION(npft) :: weibull_k_soil_gs
     REAL(wp),DIMENSION(npft) :: tau_branches
     REAL(wp),DIMENSION(npft) :: tau_sap_wood
     REAL(wp),DIMENSION(npft) :: wood_density 
     REAL(wp),DIMENSION(npft) :: k_latosa 
     REAL(wp),DIMENSION(npft) :: k_rtos 
     REAL(wp),DIMENSION(npft) :: k_crtos 
     REAL(wp),DIMENSION(npft) :: k_max_falloc_fruit
     REAL(wp),DIMENSION(npft) :: tau_fine_root
     REAL(wp),DIMENSION(npft) :: tau_coarse_root
     REAL(wp),DIMENSION(npft) :: allom_k1
     REAL(wp),DIMENSION(npft) :: allom_k2
     REAL(wp),DIMENSION(npft) :: tau_fruit
     REAL(wp),DIMENSION(npft) :: tau_seed_litter
     REAL(wp),DIMENSION(npft) :: tau_seed_est
     REAL(wp),DIMENSION(npft) :: vmax_uptake_n
     REAL(wp),DIMENSION(npft) :: soil_n_aff1 
     REAL(wp),DIMENSION(npft) :: soil_n_aff2 
     REAL(wp),DIMENSION(npft) :: vmax_uptake_p
     REAL(wp),DIMENSION(npft) :: soil_p_aff1 
     REAL(wp),DIMENSION(npft) :: soil_p_aff2 
     REAL(wp),DIMENSION(npft) :: bnf_base
     REAL(wp),DIMENSION(npft) :: lambda_est_light
     REAL(wp),DIMENSION(npft) :: k_est_light
     REAL(wp),DIMENSION(npft) :: seed_size
  CONTAINS
    PROCEDURE :: init => init_slm_parameter
  END TYPE t_slm_parameter

  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  INTEGER,PARAMETER :: itree = 1 
  INTEGER,PARAMETER :: igrass = 2 
  INTEGER,PARAMETER :: ic3phot = 1
  INTEGER,PARAMETER :: ic4phot = 2 
  INTEGER,PARAMETER :: ievergreen = 1 
  INTEGER,PARAMETER :: iraingreen = 2 
  INTEGER,PARAMETER :: isummergreen = 3 
  INTEGER,PARAMETER :: iperennial = 4 


  !-----------------------------------------------------------------------------------------------------
  !> one dimensional version of t_slm_parameter
  !!
  !! includes all parameters of t_slm_parameter, with only one value \n
  !! the pointer lctlib(:) initializes and stores one t_slm_parameter for each PFT \n
  !! see Function lctlib_wrapper below
  !-----------------------------------------------------------------------------------------------------
  TYPE t_lctlib_element
    ! this is JSBACH4 stuff
    !INTEGER  :: LctNumber            !! Unique number for each landcover type (1, 2, 3, ..)
    !CHARACTER(LEN=16)  :: LctName    !! Names of landcover types
    !INTEGER  :: LandcoverClass       !! Landcover classes (not to be confused with landcover types!):
    ! ----

  ! Respiration parameters
  REAL(wp) :: fresp_growth                 !< growth respiration per unit new biomass 
  REAL(wp) :: fmaint_rate_base             !< N-specific maintenance respiration rate (µmol CO2 / mmol N / s) 
  REAL(wp) :: t_maint_k1                   !< Coefficient in the equation for temperature sensitivity of foliar respiration 
  REAL(wp) :: t_maint_k2                   !< Coefficient in the equation for temperature sensitivity of foliar respiration
  REAL(wp) :: t_maint_k3                   !< Coefficient in the equation for temperature sensitivity of foliar respiration 

  ! Labile pool parameters 
  REAL(wp) :: tau_labile_default           !< default turnover time of the labile carbo-hydrates (in days)
  REAL(wp) :: lambda_tau_labile_temp       !< lambda in temperature response function of labile pool 
  REAL(wp) :: weibull_k_tau_labile_temp    !< k in temperature response function of labile pool 
  REAL(wp) :: lambda_tau_labile_theta      !< lambda in moisture response function of labile pool 
  REAL(wp) :: weibull_k_tau_labile_theta   !< k in moisture response function of labile pool 
  REAL(wp) :: lambda_phiphen               !< lambda in storage response function to phenology 
  REAL(wp) :: weibull_k_phiphen            !< k in in storage response function to phenology 
  REAL(wp) :: lambda_phimaint              !< lambda in storage response function to maintenance needs 
  REAL(wp) :: weibull_k_phimaint           !< k in in storage response function to maintenance needs
  REAL(wp) :: lambda_phistore              !< lambda in storage response function to storage size
  REAL(wp) :: weibull_k_phistore           !< k in in storage response function to storage size 

  ! Resporption parameters
  REAL(wp) :: resorption_default           !< default fraction of nutrient resporption before leaf shedding
  REAL(wp) :: resorption_wood              !< default fraction of nutrient resporption before wood death
  REAL(wp) :: root2leaf_cn                 !< relative C:N of fine roots compared to leaves
  REAL(wp) :: wood2leaf_cn                 !< relative C:N of woody biomass compared to leaves
  REAL(wp) :: root2leaf_np                 !< relative N of fine roots compared to leaves
  REAL(wp) :: wood2leaf_np                 !< relative N:P of woody biomass compared to leaves

  ! Parameters describing photosynthetic activity 
  REAL(wp) :: jmax2n                       !< electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
  REAL(wp) :: vcmax2n                      !< Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
  REAL(wp) :: k0_fnl_pepc                  !< fraction of N in PEP and PPKD (C4 plants only) 
  REAL(wp) :: pepc2n                       !< PEP C limited carboxylation rate per µmol CO2 / mmol N in PEP C
  REAL(wp) :: chl2n                        !< moles chlorophyll per mol of N
  REAL(wp) :: k1_fnl_struc                 !< slope of the fraction of leaf N not associated with PS against leaf N (--/mmolN) 
  REAL(wp) :: jmax2vcmax_C3                !< ratio of Jmax25/Vcmax25 for prescribed leaf N (C3 plants)
  REAL(wp) :: jmax2vcmax_C4                !< ratio of Jmax25/Vcmax25 for prescribed leaf N (C4 plants)
  REAL(wp) :: tau_canopy_acclimation       !< time scale of acclimation of canopy processes (N, P, T) (days)
  ! temperature sensitivies of PS processes
  REAL(wp) :: E0kc                         !< Baseline of Michaelis-Menten constant of Rubisco for CO2 (unit?)
  REAL(wp) :: E1kc                         !< Temperature sensitivity of Michaelis-Menten constant of Rubisco for CO2 (unit?)
  REAL(wp) :: E0ko                         !< Baseline of Michealis-Menten Kintetics for Rubisco-oxigination
  REAL(wp) :: E1ko                         !< Temperature sensitivity of Michealis-Menten Kintetics for Rubisco-oxigination
  REAL(wp) :: E0pcp                        !< Baseline of photosynthetic compensation point
  REAL(wp) :: E1pcp                        !< Temperature sensitivity of photosynthetic compensation point
  REAL(wp) :: E0v                          !< Baseline of Rubisco
  REAL(wp) :: E1v                          !< Temperature sensitivity of Rubisco
  REAL(wp) :: Tref_pepc                    !< Temperature sensitivity of PeP C4 photosynthesis 
  REAL(wp) :: Tbase_pepc                   !< Temperature sensitivity of PeP C4 photosynthesis
  REAL(wp) :: alpha_i                      !< intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I))
                                                !< µmol CO2 / mol quanta) 
  REAL(wp) :: pO2                          !< Partial Pressure of O2 in kPa
  REAL(wp) :: CiCa_default_C3              !< default Ci:Ca for C3 plants
  REAL(wp) :: CiCa_default_C4              !< default Ci:Ca for C4 plants
  REAL(wp) :: ps_it_max                    !< maximum iteration of PS calculation
  REAL(wp) :: ci_max                       !< saturating Ci in Pa C4 plants
  REAL(wp) :: Dwv2co2_air                  !< ratio of diffusion coefficient for H2O and CO2 in air 
  REAL(wp) :: Dwv2co2_turb                 !< ratio of diffusion coefficient for H2O and CO2 in turbulent air 
  ! parameters describing the fractionation due to photosynthesis
  REAL(wp) :: discr_ps_a_C13               !< discrimination of C13 due to stomatal diffusion
  REAL(wp) :: discr_ps_b_C13               !< discrimination of C13 due to Rubisco
  REAL(wp) :: discr_ps_c_C13               !< discrimination of C13 due to PEP C 
  REAL(wp) :: discr_ps_a_C14               !< discrimination of C13 due to stomatal diffusion
  REAL(wp) :: discr_ps_b_C14               !< discrimination of C13 due to Rubisco
  REAL(wp) :: discr_ps_c_C14               !< discrimination of C13 due to PEP C
  REAL(wp) :: discr_ps_phi                 !< leakage rate of bundle sheath cells
  
  ! parameters for canopy profile and canopy light extinction
  REAL(wp) :: kn                           !< extinction coefficient to describe decline of N within the canopy 
  REAL(wp) :: k0_fnl_chl_C3                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp) :: k0_fnl_chl_C4                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp) :: k1_fnl_chl_C3                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp) :: k1_fnl_chl_C4                !< parameter to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp) :: kfnl_chl                     !< parameter to describe the varying fraction of N in chlorophyll through the canopy

  ! parameters for changes in canopy N fractions
  REAl(wp) :: delta_n_leaf                 !< maximum fractional change for optimal adjsutment of canopy N allocation 
  REAL(wp) :: delta_n_fraction             !< maximum fractional change for optimal adjustment of photosynthetic N fraction

  ! parameters used for phenology calculations
  REAL(wp) :: min_lai                      !< minimum LAI beyond which all leaves are shed
  REAL(wp) :: target_lai_max               !< maximum LAI target for reserve use calculations 
  REAL(wp) :: min_leaf_age                 !< minimum leaf age before senescence is permitted (days) 
  REAL(wp) :: gdd_t_air_threshold          !< temperature threshold for growing degree days 
  REAL(wp) :: gdd_req_max                  !< maximum GDD requirement
  REAL(wp) :: k_gdd_dormance               !< exponent factor in GDD - nd_dormance equation 
  REAL(wp) :: t_air_senescence             !< weekly air temperature initiating leaf senescence
  REAL(wp) :: beta_soil_flush_grass        !< value of beta_soil at leaf flushing for grass
  REAL(wp) :: beta_soil_flush_tree         !< value of beta_soil at leaf flushing for trees 
  REAL(wp) :: beta_soil_senescence_tree    !< value of beta_soil at leaf shedding for trees

  ! Parameters required for albedo and fapar calculation
  REAL(wp) :: ka                           !< Extinction coefficient for PAR on chlorophyll
  REAL(wp) :: sigma                        !< single leaf scattering albedo 
  REAL(wp) :: kbl0                         !< extinction coefficient over black leaves 
  REAL(wp) :: kdf0                         !< extinction coefficient for diffuse PAR (Spitters eq. 7) 
  REAL(wp) :: rho2sbeta                    !< scaling factor of solar angle in reflection calculation (mo_veg_fapar) 

  ! other
  REAL(wp) :: carbon_per_dryweight_leaf    !< carbon mass per unit dry weight of leaves (gC / g DW) 
  REAL(wp) :: min_height                   !< minimum permissable height of an individual
  REAL(wp) :: sm2lm_grass                  !< stem mass to leaf mass ratio of grasses
  REAL(wp) :: lambda_fruit_alloc           !< parameter in the Weibull function controlling fruit growth
  REAL(wp) :: weibull_k_fruit_alloc        !< parameter in the Weibull function controlling fruit growth
  REAL(wp) :: offset_fruit_alloc           !< parameter in the Weibull function controlling fruit growth
  REAL(wp) :: min_falloc_fruit             !< minimum fraction of allocation going to fruit (flowers etc.)

  ! vegetation dynamics
  REAL(wp) :: background_mort_rate_tree   !< background mortality rate of trees in static mortality case
  REAL(wp) :: background_mort_rate_grass  !< background mortality rate of grasses in static mortality case

  ! soil pools
  REAL(wp) :: tau_litter                  !< soil litter pool turnover time


  ! previously stored in TYPE(t_slm_veg_param),INTENT(inout) :: veg_param

     INTEGER :: growthform 
     INTEGER :: ps_pathway
     INTEGER :: phenology_type
     REAL(wp) :: tau_leaf
     REAL(wp) :: sla
     REAL(wp) :: cn_leaf
     REAL(wp) :: np_leaf
     REAL(wp) :: k0_fnl_struc
     REAL(wp) :: tjmax_opt
     REAL(wp) :: tjmax_omega
     REAL(wp) :: g0
     REAL(wp) :: g1
     REAL(wp) :: lambda_soil_ps
     REAL(wp) :: weibull_k_soil_ps
     REAL(wp) :: lambda_soil_gs
     REAL(wp) :: weibull_k_soil_gs
     REAL(wp) :: tau_branches
     REAL(wp) :: tau_sap_wood
     REAL(wp) :: wood_density 
     REAL(wp) :: k_latosa 
     REAL(wp) :: k_rtos 
     REAL(wp) :: k_crtos 
     REAL(wp) :: k_max_falloc_fruit
     REAL(wp) :: tau_fine_root
     REAL(wp) :: tau_coarse_root
     REAL(wp) :: allom_k1
     REAL(wp) :: allom_k2
     REAL(wp) :: tau_fruit
     REAL(wp) :: tau_seed_litter
     REAL(wp) :: tau_seed_est
     REAL(wp) :: vmax_uptake_n
     REAL(wp) :: soil_n_aff1 
     REAL(wp) :: soil_n_aff2 
     REAL(wp) :: vmax_uptake_p
     REAL(wp) :: soil_p_aff1 
     REAL(wp) :: soil_p_aff2 
     REAL(wp) :: bnf_base
     REAL(wp) :: lambda_est_light
     REAL(wp) :: k_est_light
     REAL(wp) :: seed_size
  END TYPE t_lctlib_element

!  PRIVATE
!  PUBLIC :: 

CONTAINS 

  !-----------------------------------------------------------------------------------------------------
  !> initialize parameter values
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_slm_parameter(this)

    USE mo_shared_constants,        ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2 
    USE mo_atmland_constants

    IMPLICIT NONE
      ! 0.1 InOut
      CLASS(t_slm_parameter), INTENT(out) :: this        !< the slm parameter TYPE
 
      !> plant functional types are: \n
      !! 1) broad leaved evergreen (rainforest) \n
      !! 2) broad leaved evergreen (xeric forest) \n
      !! 3) broad leaved deciduous (rain green) \n
      !! 4) broad leaved deciduous (summer green) \n
      !! 5) needle-leaved evergreen \n
      !! 6) needle-leaved deciduous \n
      !! 7) C3 grass \n
      !! 8) C4 grass

      ! growth respiration per unit new biomass
      this%fresp_growth(:) = 0.25_wp
      ! N-specific foliar maintenance respiration rate (µmol CO2 / mol N / s)
      ! Taken from Sprugel et al. Respiration from the organ level to the stand. 
      ! In: Smith WK, Hinckley, T.M., editor. Physiological Ecology of Coniferous Forests.
      ! San Diego, CA, USA: Academic Press; 1996.
      ! Note converted from 7.4 e-7 to 1.25 e-6 to account for the different values of the Lloyd&Taylor equation
      ! such that the applied rates a comparable at 15-20°C
      this%fmaint_rate_base(:) = 1._wp 
      ! parameters in the Lloyd and Taylor temperature response function (Lloyd & Taylor, 1996, Funct. Ecol.) 
      this%t_maint_k1(:) = 308.56_wp  
      this%t_maint_k2(:) = 56.02_wp
      this%t_maint_k3(:) = 227.13_wp

      ! default turnover time of the labile carbo-hydrates (in days)
      this%tau_labile_default(:) = 0.2_wp 
      ! temperature scalar on turnover time of labile pool (1-exp(-(lambda * T)^k))
      ! tuned to reduce growth between 5°C to zero at 0°C
      this%lambda_tau_labile_temp(:) = 0.5_wp
      this%weibull_k_tau_labile_temp(:) = 2.0_wp
      ! moisture scalar on turnover time of labile pool (1-exp(-(lambda * theta)^k))
      ! tuned to reduce growth from 20% saturation to zero at 0% saturation
      this%lambda_tau_labile_theta(:) = 0.1_wp
      this%weibull_k_tau_labile_theta(:) = 2.0_wp
      ! pull from reserve to labile pool due to phenology (1-exp(-(lambda * lai/target_lai)^k))
      this%lambda_phiphen(:) = 1.8_wp
      this%weibull_k_phiphen(:) = 4.0_wp
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      this%lambda_phimaint(:) = 4.0_wp
      this%weibull_k_phimaint(:) = 1.2_wp
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      this%lambda_phistore(:) = 2._wp
      this%weibull_k_phistore(:) = 3._wp

      ! default fraction of nutrient resporption before leaf shedding
      this%resorption_default(:) = 0.5_wp 
      ! default fraction of nutrient resporption before wood death 
      this%resorption_wood(:) = 0.2_wp 
      ! ratio of root C:N to leaf C:N
      this%root2leaf_cn(:) = 0.85_wp 
      ! ratio of root C:N to leaf C:N
      this%wood2leaf_cn(:) = 0.145_wp 
      ! ratio of root N:P to leaf N:P
      this%root2leaf_np(:) = 1.0_wp 
      ! ratio of root N:P to leaf N:P
      this%wood2leaf_np(:) = 1.0_wp 

      ! Parameters used in the photosynthesis routine
      ! Partial Pressure of O2 in kPa
      this%pO2(:) = 20.9_wp
      ! Temperature sensitivity and baseline of Michaelis-Menten constant of Rubisco for CO2 (unit?)
      this%E0kc(:) = 38.05_wp
      this%E1kc(:) = 79430.0_wp
      ! Temperature sensitivity and baseline of Michealis-Menten Kintetics for Rubisco-oxigination     
      this%E0ko(:) = 20.3_wp
      this%E1ko(:) = 36380.0_wp
      ! Temperature sensitivity and baseline of photosynthetic compensation point
      this%E0pcp(:) = 19.02_wp
      this%E1pcp(:) = 37830.0_wp
      ! Temperature sensitivity and baseline of Rubisco
      this%E0v(:) = 26.35_wp
      this%E1v(:) = 65330.0_wp
   ! Temperature sensitivity and baseline value of foliar respiration
   !this%E0kt(:) = 18.72_wp
   !this%E1kt(:) = 46390.0_wp
      ! Temperature sensitivity of PeP C4 photosynthesis
      this%Tref_pepc(:) = 25._wp
      this%Tbase_pepc(:) = 10._wp
      ! intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I))) (µmol CO2 / mol quanta)
      ! accounting for the absorptance of leaves (0.85), the maximum quantum yield of PS II (0.7)
      ! the fraction of total light that reaches PS II (0.5) and the NADPH limited RubGeneration (jointly leading to 0.066)
      this%alpha_i(:) = 0.85_wp * 0.066_wp 
      ! initial guess of the leaf Ci to Ca ratio
      this%CiCa_default_C3(:) = 0.7_wp
      this%CiCa_default_C4(:) = 0.4_wp
      this%ci_max(:) = 7800.0_wp ! saturating Ci in Pa
      ! maximum number of iteration in photosynthesis calculation (should probably be in mo_base_constants)
      this%ps_it_max(:) = 10.0_wp
      ! ratio of the molecular diffusion constants of water and CO2 in air (Nobel, Plant Phys. 2005, p 380)
      this%Dwv2co2_air(:) = Dwv / Dco2 
      ! ratio of the diffusion constant of water and CO2 in air, accounting for turbulent transfer (Nobel, Plant Phys. 2005, p 380) 
      this%Dwv2co2_turb(:) = ( Dwv / Dco2  ) ** (2._wp/3._wp)

      ! carbon content of leaves per dry matter
      this%carbon_per_dryweight_leaf(:) = 0.48_wp ! g C / g DW
      ! Parameters describing the PS ~ An relationship
      ! Values from Niinemets 1997 (SC 01/16)
      this%jmax2n(:) = 17.6_wp/4.0_wp   ! electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
      this%vcmax2n(:) = 1.8_wp          ! Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
      this%k0_fnl_pepc(:) = 0.045_wp    ! fraction of N in PEP and PPDK (Makino et al. 2003 / Tazoe et al. PCE 2006, 29, 691-700)
      ! PEP C limited carboxylation rate per µmol CO2 / mmol N in PEP C, derived from k = 0.7 mol m-2 s-1 (Collatz et al. 1992)
      ! assuming 157 mmol N / m-2 leaf N and 4.5 % of leaf N in PS 
      this%pepc2n(:) = 98777.97_wp  ! PEP C limited carboxylation rate per µmol CO2 / mmol N in PEP C
      ! Value from Friend et al. 2009 
      ! values from Evans 1989, Oecologia, 78: 9-19 (Table 2)
      this%chl2n(:) = 1000._wp/39.8_wp  ! µmol Chloroplast per mmol N (light harvesting complex)
      ! slope of the relationship of structural leaf N against leaf N in fraction/mmol (Friend 1997)
      this%k1_fnl_struc(:) = 71.4_wp * molar_mass_N / 1000000._wp
      ! in static leaf N mode: prescribed ratio of Jmax25/Vcmax25
      this%jmax2vcmax_C3(:) = 1.92_wp      ! text book
      this%jmax2vcmax_C4(:) = 1.4_wp       ! calibrated to have Jmax=vcmax at 25 deg
      ! time scale of canopy acclimation (days)
      this%tau_canopy_acclimation(:) = 10._wp 
 
     ! parameters describing the isotopic discrimination of PS, derived from Drake 2014, Radiocarbon, 56, 29-38
      this%discr_ps_a_C13(:) = 4.4_wp  
      this%discr_ps_b_C13(:) = 27.0_wp
      this%discr_ps_c_C13(:) = 5.7_wp
      this%discr_ps_a_C14(:) = this%discr_ps_a_C13(:) * 1.97_wp
      this%discr_ps_b_C14(:) = this%discr_ps_b_C13(:) * 1.97_wp
      this%discr_ps_c_C14(:) = this%discr_ps_c_C13(:) * 1.89_wp 
      this%discr_ps_phi(:) = 0.16_wp ! a typical value cf Hatch et al. 1995, Plant Phys. 108, 173-181

      ! parameters describing the canopy profile of N 
      !  extinction coefficient to describe decline of N within the canopy 
      this%kn(:) = 0.11_wp
      ! parameters to describe the varying fraction of N in chlorophyll through the canopy
      this%k0_fnl_chl_C3(:) = 6.0_wp
      this%k0_fnl_chl_C4(:) = 15.0_wp
      this%k1_fnl_chl_C3(:) = 3.6_wp
      this%k1_fnl_chl_C4(:) = 4.4_wp
      this%kfnl_chl(:) = 0.7_wp

      ! Parameters required for albedo and fapar calculation
      ! single leaf scattering albedo
      this%sigma(:) = 0.2_wp
      ! extinction coefficient over black leaves
      this%kbl0(:) = 0.5_wp
      ! extinction coefficient for diffuse PAR (Spitters eq. 7)
      this%kdf0(:) = 0.8_wp
      ! scaling factor of solar angle in reflection calculation (mo_veg_fapar)
      this%rho2sbeta(:) = 1.6_wp
      ! Extinction coefficient for PAR on chlorophyll
      this%ka(:) = 0.005_wp

      ! maximum change in CN ratio of canopy
      this%delta_n_leaf(:) = 0.001_wp
      ! maximum change in photosynthetic N fractions
      this%delta_n_fraction(:) = 0.01_wp
  
      ! minimum LAI for phenology calculations
      this%min_lai(:) = 0.1_wp
      ! minimum LAI for phenology calculations
      this%target_lai_max(:) = 5.0_wp
      ! temperature threshold for growing degree days 
      this%gdd_t_air_threshold(:)= 5.0_wp + Tzero
      ! maximal GDD requirement (assuming no dormance) 
      this%gdd_req_max(:) = 800._wp
      ! temperature threshold for growing degree days 
      this%k_gdd_dormance(:) = 0.0098_wp
      ! weekly temperature initiating leaf senescence
      this%t_air_senescence(:) = 5.0_wp + Tzero
      ! value of beta_soil for leaf flushing for trees
      this%beta_soil_flush_grass(:) = 0.9_wp
      ! value of beta_soil for leaf flushing for trees
      this%beta_soil_flush_tree(:) = 0.98_wp
      ! value of beta_soil at shedding of leaves for trees
      this%beta_soil_senescence_tree(:)  = 0.5_wp
      ! minimum leaf age beyond which phenological senescence can set in (days) 
      this%min_leaf_age(:) = 10.0_wp

      ! background mortality rate for trees in static mortality case (1/year)
      this%background_mort_rate_tree(:) = 0.01_wp
      ! background mortality rate for grass in static mortality case (1/year)
      this%background_mort_rate_grass(:) = 0.01_wp

      ! minimum permissable height of an individual
      this%min_height(:) = 0.1_wp
      ! stem mass to leaf mass ratio of grasses
      this%sm2lm_grass(:) = 0.05_wp
      ! Parameters controlling fruit allocation in: exp(-(lambda * reserve_use+offset)^k) 
      this%lambda_fruit_alloc(:) = 10.0_wp
      this%weibull_k_fruit_alloc(:) = 2.0_wp
      this%offset_fruit_alloc(:) = 0.1_wp ! ~µmol C / m2 / s
      ! minimum fraction of allocation going to fruit (flowers etc.)
      this%min_falloc_fruit(:) = 0.01_wp
      ! soil litter turnover pool, s equivalent to 1 year
      this%tau_litter(:) = 1.0_wp * 3600.0_wp * 24.0_wp * 365.0_wp


      ! previously stored in TYPE(t_slm_veg_param),INTENT(inout) :: veg_param

      ! vegetation growth form
      this%growthform(:) = & 
      (/itree,itree,itree,itree,itree,itree,igrass,igrass/)
      ! photosynthetic path way 
      this%ps_pathway(:) = & 
      (/ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic4phot/)
      ! photosynthetic path way 
      this%phenology_type(:) = & 
      (/ievergreen,ievergreen,iraingreen,isummergreen,ievergreen,isummergreen,iperennial,iperennial/)
      ! average turnover (1/lifetime) of a leaf (1/years) 
      ! From Kattge et al. 2011, Table 5, including conversion from month lifetime to annual turnover
      this%tau_leaf(:) = 12._wp / & 
      (/16.83_wp,15.88_wp,5.83_wp,5.83_wp,39.7_wp,6.08_wp,3.85_wp,3.85_wp/)
      ! specific leaf area of a plant (m2 / mol C) 
      ! From Kattge et al. 2011, Table 5, including conversion from mm/mgDW to m2/mol
      this%sla(:) = 1._wp / this%carbon_per_dryweight_leaf(:) / 1000.0_wp * molar_mass_C * &
      (/9.46_wp,8.99_wp,15.4_wp,15.4_wp,5.0_wp,10.09_wp,20.12_wp,19.23_wp/) 
      ! default leaf C:N (molC/molN) if not calculated dynamically 
      ! From Kattge et al. 2011, Table 5, including conversion from mg/g to mol/mol
      this%cn_leaf(:) = this%carbon_per_dryweight_leaf(:) * 1000._wp * molar_mass_N / molar_mass_C / & 
      (/16.89_wp,13.73_wp,21.32_wp,21.32_wp,12.09_wp,19.37_wp,17.84_wp,14.14_wp/) 
      ! default leaf N:P (molN/molP) if not calculated dynamically 
      ! From Kattge et al. 2011, Table 5, including conversion from mg/g to mol/mol
      this%np_leaf(:) = this%carbon_per_dryweight_leaf(:) * 1000._wp * molar_mass_P / molar_mass_C / this%cn_leaf(:) / &
      (/0.86_wp,0.84_wp,1.44_wp,1.44_wp,1.23_wp,1.83_wp,1.43_wp,1.36_wp/)
      ! fraction of leaf N not used for photosynthesis (Friend 1997) 
      ! temporarily tuned to match Amax/N in Kattge etal. 2012, Table 5
      this%k0_fnl_struc(:) = & 
      !(/0.67_wp,0.67_wp,0.67_wp,0.67_wp,0.83_wp,0.83_wp,0.67_wp,0.67_wp/) 
      (/0.75_wp,0.70_wp,0.63_wp,0.63_wp,0.83_wp,0.79_wp,0.3_wp,0.3_wp/) 
      !Le Maire et al. 2012 (derived from grasses, does not give plausble results for needles)
      !0.012_wp / 0.48_wp * molar_mass_C / molar_mass_N / this%sla(:)/&
      !(1._wp / this%sla(:) / this%cn_leaf(:)) 
      ! initial value of the temperature optimum of electron transport and its shape parameter (deg C)
      this%tjmax_opt(:) = & 
      (/17._wp,17._wp,17._wp,17._wp,17._wp,17._wp,17._wp,35._wp/)
      this%tjmax_omega(:) = 18._wp
      ! intercept and slope of the An ~gs relationship (mmol/µmol)
      this%g0(:) = & ! 0.01_wp
      (/0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.03_wp,0.03_wp/)
      this%g1(:) = & !9.3_wp
      ! g1 in Ball & Berry relationship, tuned to get similar Amax/N values as with Lin et al. g1 
      (/9.3_wp,8.3_wp,7.0_wp,10.9_wp,5.5_wp,7.0_wp,9.3_wp,2.0_wp/)
      ! g1 in kPa^0.5 as in Lin et al. 2015, Nat CC. OBS: no values for Larches - assigned broad leaved value!
      !(/3.77_wp,3.37_wp,2.98_wp,4.64_wp,2.35_wp,4.64_wp,4.16_wp,1.62_wp/)
      ! Weibull parameters describing the response of leaf-level photosynthesis to 
      ! fractional water availability
      this%lambda_soil_ps(:) = 3.0_wp
      this%weibull_k_soil_ps(:) = 3.0_wp
      ! Weibull parameters describing the response of stomatal conductance to 
      ! fractional water availability
      this%lambda_soil_gs(:) = this%lambda_soil_ps(:)
      this%weibull_k_soil_gs(:) = this%weibull_k_soil_ps(:)
      ! turnover time of the fraction of sapwood that is in branches (10%) (1/years)
      this%tau_branches(:) = 0.01_wp
      ! turnover time of the sapwood pool (1/years)
      this%tau_sap_wood(:) = 0.05_wp
      ! wood density (mol C / m3)
      this%wood_density(:) = & 
      (/30000.0_wp,30000.0_wp,30000.0_wp,30000.0_wp,30000.0_wp,30000.0_wp,22000._wp,22000._wp/) / molar_mass_C
      ! leaf area to sapwood area ratio 
      this%k_latosa(:) = 4000.0_wp 
      ! coarse root to sapwood mass ratio 
      this%k_crtos(:) = 0.2_wp
      ! trade-off parametr for hydraulic investment into sapwood or fine roots 
      WHERE(this%growthform /= igrass)
          this%k_rtos(:) = 1.0_wp / (this%k_latosa(:) / this%sla(:) / this%wood_density(:))
        ELSEWHERE
          this%k_rtos(:) = 0.5_wp / this%sm2lm_grass
      END WHERE
      ! maximum fraction of biomass growth going to fruits 
      WHERE(this%growthform /= igrass)
          this%k_max_falloc_fruit(:) = 0.1_wp 
        ELSEWHERE
          this%k_max_falloc_fruit(:) = 0.2_wp 
      END WHERE
      ! average turnover (1/lifetime) of a fine root (1/years) 
      this%tau_fine_root(:) = & 
      (/0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp,1.2_wp,1.2_wp/) 
      ! average turnover (1/lifetime) of a coarse root (1/years) 
      this%tau_coarse_root(:) = 0.1_wp 
      ! paramter in height diameter relationship 
      this%allom_k1(:) = & 
      (/55.0_wp,55.0_wp,55.0_wp,55.0_wp,55.0_wp,55.0_wp,20._wp,20._wp/) 
      ! paramter in height diameter relationship
      this%allom_k2(:) = & 
      (/0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.5_wp,0.5_wp/) 
      ! turnover time of the fruit pool (1/years)
      this%tau_fruit(:) = 1.0_wp
      ! turnover time of the seed bed to litter (1/years)
      this%tau_seed_litter(:) = 0.1_wp
      ! turnover time of the seed bed to establishment (1/years)
      this%tau_seed_est(:) = 0.3_wp
      ! vmax of N uptake (~µmol N / mol C / s)
      this%vmax_uptake_n(:) = 14._wp !714.28_wp   ! tuned from OCN's 1.5 gN/gC/1800s
      ! low-affinity uptake (L/mmol) 
      this%soil_n_aff1(:) = 0.105_wp      ! tuned from OCN's 0.05 m2/gN assuming 150 L per m2
      ! half saturation of high-affinity uptake (mmol/L) 
      this%soil_n_aff2(:) = 0.4_wp        ! tuned from OCN's 0.84 gN/m2 assuming 150 L per m2
      ! vmax of P uptake (~µmol P / mol C / s)
      this%vmax_uptake_p(:) = 14.28_wp   ! tuned from OCN's 1.5 gN/gC/1800s
      ! low-affinity uptake (L/mmol) 
      this%soil_p_aff1(:) = 0.105_wp      ! tuned from OCN's 0.05 m2/gN assuming 150 L per m2
      ! half saturation of high-affinity uptake (mmol/L) 
      this%soil_p_aff2(:) = 0.4_wp        ! tuned from OCN's 0.84 gN/m2 assuming 150 L per m2
      ! placeholder constant BNf rate, required to keep soil NH4 at constant levels mumol N s-1 
      this%bnf_base(:) = & 
      (/2.6_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9, 5.0_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9/) 
      ! parameter in the Weibull function controlling light-limited establishment 
      this%lambda_est_light(:) = 1.1_wp        
      ! parameter in the Weibull function controlling light-limited establishment
      this%k_est_light(:) = 4.0_wp
      ! seed size (mol C)
      this%seed_size(:) = 0.1_wp
    END SUBROUTINE init_slm_parameter 

  !-----------------------------------------------------------------------------------------------------
  !> pass parameter values from slm_parameter to lctlib
  !!
  !! this Function has NO functionality other than changing the STRUCTURE of how the parameter values are stored \n
  !! Goal: increasing compatibility with jsbach4
  !-----------------------------------------------------------------------------------------------------
    FUNCTION slm_lctlib_wrapper(slm_parameter) RESULT(lctlib)

      USE mo_base_constants, ONLY: npft
      IMPLICIT NONE

      TYPE(t_slm_parameter),INTENT(in)     :: slm_parameter
      TYPE(t_lctlib_element),POINTER       :: lctlib(:)

      ALLOCATE(lctlib(npft))    !    jsbach4: ALLOCATE(lctlib(nlct))

      !> pass parameters from slm_parameter to lctlib
      !! 
!      ! Respiration parameters
!      lctlib(1:npft)%fresp_growth               = slm_parameter%fresp_growth(1:npft)
!      lctlib(1:npft)%fmaint_rate_base           = slm_parameter%fmaint_rate_base(1:npft)
!      lctlib(1:npft)%t_maint_k1                 = slm_parameter%t_maint_k1(1:npft)
!      lctlib(1:npft)%t_maint_k2                 = slm_parameter%t_maint_k2(1:npft)
!      lctlib(1:npft)%t_maint_k3                 = slm_parameter%t_maint_k3(1:npft)

!      ! Labile pool parameters 
!      lctlib(1:npft)%tau_labile_default         = slm_parameter%tau_labile_default(1:npft)
!      lctlib(1:npft)%lambda_tau_labile_temp     = slm_parameter%lambda_tau_labile_temp(1:npft)
!      lctlib(1:npft)%weibull_k_tau_labile_temp  = slm_parameter%weibull_k_tau_labile_temp(1:npft)
!      lctlib(1:npft)%lambda_tau_labile_theta    = slm_parameter%lambda_tau_labile_theta(1:npft)
!      lctlib(1:npft)%weibull_k_tau_labile_theta = slm_parameter%weibull_k_tau_labile_theta(1:npft)
!      lctlib(1:npft)%lambda_phiphen             = slm_parameter%lambda_phiphen(1:npft)
!      lctlib(1:npft)%weibull_k_phiphen          = slm_parameter%weibull_k_phiphen(1:npft)
!      lctlib(1:npft)%lambda_phimaint            = slm_parameter%lambda_phimaint(1:npft)
!      lctlib(1:npft)%weibull_k_phimaint         = slm_parameter%weibull_k_phimaint(1:npft)
!      lctlib(1:npft)%lambda_phistore            = slm_parameter%lambda_phistore(1:npft)
!      lctlib(1:npft)%weibull_k_phistore         = slm_parameter%weibull_k_phistore(1:npft)

!      ! Resporption parameters
!      lctlib(1:npft)%resorption_default         = slm_parameter%resorption_default(1:npft)
!      lctlib(1:npft)%resorption_wood            = slm_parameter%resorption_wood(1:npft)
!      lctlib(1:npft)%root2leaf_cn               = slm_parameter%root2leaf_cn(1:npft)
!      lctlib(1:npft)%wood2leaf_cn               = slm_parameter%wood2leaf_cn(1:npft)
!      lctlib(1:npft)%root2leaf_np               = slm_parameter%root2leaf_np(1:npft)
!      lctlib(1:npft)%wood2leaf_np               = slm_parameter%wood2leaf_np(1:npft)

!      ! Parameters describing photosynthetic activity 
!      lctlib(1:npft)%jmax2n                     = slm_parameter%jmax2n(1:npft)
!      lctlib(1:npft)%vcmax2n                    = slm_parameter%vcmax2n(1:npft)
!      lctlib(1:npft)%k0_fnl_pepc                = slm_parameter%k0_fnl_pepc(1:npft)
!      lctlib(1:npft)%pepc2n                     = slm_parameter%pepc2n(1:npft)
!      lctlib(1:npft)%chl2n                      = slm_parameter%chl2n(1:npft)
!      lctlib(1:npft)%k1_fnl_struc               = slm_parameter%k1_fnl_struc(1:npft)
!      lctlib(1:npft)%jmax2vcmax_C3              = slm_parameter%jmax2vcmax_C3(1:npft)
!      lctlib(1:npft)%jmax2vcmax_C4              = slm_parameter%jmax2vcmax_C4(1:npft)
!      lctlib(1:npft)%tau_canopy_acclimation     = slm_parameter%tau_canopy_acclimation(1:npft)
!      ! temperature sensitivies of PS processes
!      lctlib(1:npft)%E0kc                       = slm_parameter%E0kc(1:npft)
!      lctlib(1:npft)%E1kc                       = slm_parameter%E1kc(1:npft)
!      lctlib(1:npft)%E0ko                       = slm_parameter%E0ko(1:npft)
!      lctlib(1:npft)%E1ko                       = slm_parameter%E1ko(1:npft)
!      lctlib(1:npft)%E0pcp                      = slm_parameter%E0pcp(1:npft)
!      lctlib(1:npft)%E1pcp                      = slm_parameter%E1pcp(1:npft)
!      lctlib(1:npft)%E0v                        = slm_parameter%E0v(1:npft)
!      lctlib(1:npft)%E1v                        = slm_parameter%E1v(1:npft)
!      lctlib(1:npft)%Tref_pepc                  = slm_parameter%Tref_pepc(1:npft)
!      lctlib(1:npft)%Tbase_pepc                 = slm_parameter%Tbase_pepc(1:npft)
!      lctlib(1:npft)%alpha_i                    = slm_parameter%alpha_i(1:npft)
!                            
!      lctlib(1:npft)%pO2                        = slm_parameter%pO2(1:npft)
!      lctlib(1:npft)%CiCa_default_C3            = slm_parameter%CiCa_default_C3(1:npft)
!      lctlib(1:npft)%CiCa_default_C4            = slm_parameter%CiCa_default_C4(1:npft)
!      lctlib(1:npft)%ps_it_max                  = slm_parameter%ps_it_max(1:npft)
!      lctlib(1:npft)%ci_max                     = slm_parameter%ci_max(1:npft)
!      lctlib(1:npft)%Dwv2co2_air                = slm_parameter%Dwv2co2_air(1:npft)
!      lctlib(1:npft)%Dwv2co2_turb               = slm_parameter%Dwv2co2_turb(1:npft)
!      ! parameters describing the fractionation due to photosynthesis
!      lctlib(1:npft)%discr_ps_a_C13             = slm_parameter%discr_ps_a_C13(1:npft)
!      lctlib(1:npft)%discr_ps_b_C13             = slm_parameter%discr_ps_b_C13(1:npft)
!      lctlib(1:npft)%discr_ps_c_C13             = slm_parameter%discr_ps_c_C13(1:npft)
!      lctlib(1:npft)%discr_ps_a_C14             = slm_parameter%discr_ps_a_C14(1:npft)
!      lctlib(1:npft)%discr_ps_b_C14             = slm_parameter%discr_ps_b_C14(1:npft)
!      lctlib(1:npft)%discr_ps_c_C14             = slm_parameter%discr_ps_c_C14(1:npft)
!      lctlib(1:npft)%discr_ps_phi               = slm_parameter%discr_ps_phi(1:npft)
!  
!      ! parameters for canopy profile and canopy light extinction
!      lctlib(1:npft)%kn                         = slm_parameter%kn(1:npft)
!      lctlib(1:npft)%k0_fnl_chl_C3              = slm_parameter%k0_fnl_chl_C3(1:npft)
!      lctlib(1:npft)%k0_fnl_chl_C4              = slm_parameter%k0_fnl_chl_C4(1:npft)
!      lctlib(1:npft)%k1_fnl_chl_C3              = slm_parameter%k1_fnl_chl_C3(1:npft)
!      lctlib(1:npft)%k1_fnl_chl_C4              = slm_parameter%k1_fnl_chl_C4(1:npft)
!      lctlib(1:npft)%kfnl_chl                   = slm_parameter%kfnl_chl(1:npft)

!      ! parameters for changes in canopy N fractions
!      lctlib(1:npft)%delta_n_leaf               = slm_parameter%delta_n_leaf(1:npft)
!      lctlib(1:npft)%delta_n_fraction           = slm_parameter%delta_n_fraction(1:npft)

!      ! parameters used for phenology calculations
!      lctlib(1:npft)%min_lai                    = slm_parameter%min_lai(1:npft)
!      lctlib(1:npft)%target_lai_max             = slm_parameter%target_lai_max(1:npft)
!      lctlib(1:npft)%min_leaf_age               = slm_parameter%min_leaf_age(1:npft)
!      lctlib(1:npft)%gdd_t_air_threshold        = slm_parameter%gdd_t_air_threshold(1:npft)
!      lctlib(1:npft)%gdd_req_max                = slm_parameter%gdd_req_max(1:npft)
!      lctlib(1:npft)%k_gdd_dormance             = slm_parameter%k_gdd_dormance(1:npft)
!      lctlib(1:npft)%t_air_senescence           = slm_parameter%t_air_senescence(1:npft)
!      lctlib(1:npft)%beta_soil_flush_grass      = slm_parameter%beta_soil_flush_grass(1:npft)
!      lctlib(1:npft)%beta_soil_flush_tree       = slm_parameter%beta_soil_flush_tree(1:npft)
!      lctlib(1:npft)%beta_soil_senescence_tree  = slm_parameter%beta_soil_senescence_tree(1:npft)

!      ! Parameters required for albedo and fapar calculation
!      lctlib(1:npft)%ka                         = slm_parameter%ka(1:npft)
!      lctlib(1:npft)%sigma                      = slm_parameter%sigma(1:npft)
!      lctlib(1:npft)%kbl0                       = slm_parameter%kbl0(1:npft)
!      lctlib(1:npft)%kdf0                       = slm_parameter%kdf0(1:npft)
!      lctlib(1:npft)%rho2sbeta                  = slm_parameter%rho2sbeta(1:npft)

!      ! other
!      lctlib(1:npft)%carbon_per_dryweight_leaf  = slm_parameter%carbon_per_dryweight_leaf(1:npft)
!      lctlib(1:npft)%min_height                 = slm_parameter%min_height(1:npft)
!      lctlib(1:npft)%sm2lm_grass                = slm_parameter%sm2lm_grass(1:npft)
!      lctlib(1:npft)%lambda_fruit_alloc         = slm_parameter%lambda_fruit_alloc(1:npft)
!      lctlib(1:npft)%weibull_k_fruit_alloc      = slm_parameter%weibull_k_fruit_alloc(1:npft)
!      lctlib(1:npft)%offset_fruit_alloc         = slm_parameter%offset_fruit_alloc(1:npft)
!      lctlib(1:npft)%min_falloc_fruit           = slm_parameter%min_falloc_fruit(1:npft)

!      ! vegetation dynamics
!      lctlib(1:npft)%background_mort_rate_tree  = slm_parameter%background_mort_rate_tree(1:npft)
!      lctlib(1:npft)%background_mort_rate_grass = slm_parameter%background_mort_rate_grass(1:npft)

!      ! soil pools
!      lctlib(1:npft)%tau_litter                 = slm_parameter%tau_litter(1:npft)

      ! previously stored in TYPE(t_slm_veg_param),INTENT(inout) :: veg_param
         lctlib(1:npft)%growthform        = slm_parameter%growthform(1:npft)
         lctlib(1:npft)%ps_pathway        = slm_parameter%ps_pathway(1:npft)
         lctlib(1:npft)%phenology_type    = slm_parameter%phenology_type(1:npft)
         lctlib(1:npft)%tau_leaf          = slm_parameter%tau_leaf(1:npft)
         lctlib(1:npft)%sla               = slm_parameter%sla(1:npft)
         lctlib(1:npft)%cn_leaf           = slm_parameter%cn_leaf(1:npft)
         lctlib(1:npft)%np_leaf           = slm_parameter%np_leaf(1:npft)
         lctlib(1:npft)%k0_fnl_struc      = slm_parameter%k0_fnl_struc(1:npft)
         lctlib(1:npft)%tjmax_opt         = slm_parameter%tjmax_opt(1:npft)
         lctlib(1:npft)%tjmax_omega       = slm_parameter%tjmax_omega(1:npft)
         lctlib(1:npft)%g0                = slm_parameter%g0(1:npft)
         lctlib(1:npft)%g1                = slm_parameter%g1(1:npft)
         lctlib(1:npft)%lambda_soil_ps    = slm_parameter%lambda_soil_ps(1:npft)
         lctlib(1:npft)%weibull_k_soil_ps = slm_parameter%weibull_k_soil_ps(1:npft)
         lctlib(1:npft)%lambda_soil_gs    = slm_parameter%lambda_soil_gs(1:npft)
         lctlib(1:npft)%weibull_k_soil_gs = slm_parameter%weibull_k_soil_gs(1:npft)
         lctlib(1:npft)%tau_branches      = slm_parameter%tau_branches(1:npft)
         lctlib(1:npft)%tau_sap_wood      = slm_parameter%tau_sap_wood(1:npft)
         lctlib(1:npft)%wood_density      = slm_parameter%wood_density(1:npft)
         lctlib(1:npft)%k_latosa          = slm_parameter%k_latosa(1:npft)
         lctlib(1:npft)%k_rtos            = slm_parameter%k_rtos(1:npft)
         lctlib(1:npft)%k_crtos           = slm_parameter%k_crtos(1:npft)
         lctlib(1:npft)%k_max_falloc_fruit = slm_parameter%k_max_falloc_fruit(1:npft)
         lctlib(1:npft)%tau_fine_root     = slm_parameter%tau_fine_root(1:npft)
         lctlib(1:npft)%tau_coarse_root   = slm_parameter%tau_coarse_root(1:npft)
         lctlib(1:npft)%allom_k1          = slm_parameter%allom_k1(1:npft)
         lctlib(1:npft)%allom_k2          = slm_parameter%allom_k2(1:npft)
         lctlib(1:npft)%tau_fruit         = slm_parameter%tau_fruit(1:npft)
         lctlib(1:npft)%tau_seed_litter   = slm_parameter%tau_seed_litter(1:npft)
         lctlib(1:npft)%tau_seed_est      = slm_parameter%tau_seed_est(1:npft)
         lctlib(1:npft)%vmax_uptake_n     = slm_parameter%vmax_uptake_n(1:npft)
         lctlib(1:npft)%soil_n_aff1       = slm_parameter%soil_n_aff1(1:npft)
         lctlib(1:npft)%soil_n_aff2       = slm_parameter%soil_n_aff2(1:npft)
         lctlib(1:npft)%vmax_uptake_p     = slm_parameter%vmax_uptake_p(1:npft)
         lctlib(1:npft)%soil_p_aff1       = slm_parameter%soil_p_aff1(1:npft)
         lctlib(1:npft)%soil_p_aff2       = slm_parameter%soil_p_aff2(1:npft)
         lctlib(1:npft)%bnf_base          = slm_parameter%bnf_base(1:npft)
         lctlib(1:npft)%lambda_est_light  = slm_parameter%lambda_est_light(1:npft)
         lctlib(1:npft)%k_est_light       = slm_parameter%k_est_light(1:npft)
         lctlib(1:npft)%seed_size         = slm_parameter%seed_size(1:npft)

    END FUNCTION slm_lctlib_wrapper


END MODULE mo_base_slm_parameter_class
