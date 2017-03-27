!>
!! @brief Contains constants for the vegetation process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-07-21)
!!
MODULE mo_veg_constants

  USE mo_kind,                ONLY: wp

  IMPLICIT NONE
  PUBLIC


  ! Respiration parameters
  REAL(wp), SAVE :: &      ! In Fortran 2003 and Fortran 2008, a maximum of 255 continuation lines are allowed
      ! growth respiration per unit new biomass
      fresp_growth  ,     & !< growth respiration per unit new biomass 
      ! N-specific foliar maintenance respiration rate (µmol CO2 / mol N / s)
      ! Taken from Sprugel et al. Respiration from the organ level to the stand. 
      ! In: Smith WK, Hinckley, T.M., editor. Physiological Ecology of Coniferous Forests.
      ! San Diego, CA, USA: Academic Press; 1996.
      ! Note converted from 7.4 e-7 to 1.25 e-6 to account for the different values of the Lloyd&Taylor equation
      ! such that the applied rates a comparable at 15-20°C
      fmaint_rate_base   ,     & !< N-specific maintenance respiration rate (µmol CO2 / mmol N / s) 
      ! parameters in the Lloyd and Taylor temperature response function (Lloyd & Taylor, 1996, Funct. Ecol.) 
      t_maint_k1    ,     & !< Coefficient in the equation for temperature sensitivity of foliar respiration
      t_maint_k2  ,     &  !< Coefficient in the equation for temperature sensitivity of foliar respiration
      t_maint_k3  , &      !< Coefficient in the equation for temperature sensitivity of foliar respiration
      f_resp_acclim ,   &    !< respiration temperature acclimation factor 
      t_acclim_zero        !< base temperature for respiration acclimation 

  ! Labile pool parameters 
  REAL(wp), SAVE :: &
      ! default turnover time of the labile carbo-hydrates (in days)
      tau_labile,           &  !< default turnover time of the labile carbo-hydrates (in days)
      ! temperature scalar on turnover time of labile pool (1-exp(-(lambda * T)^k))
      ! tuned to reduce growth between 5°C to zero at 0°C
      lambda_labile_temp,   & !< lambda in temperature response function of labile pool 
      k_labile_temp,        & !< k in temperature response function of labile pool 
      ! moisture scalar on turnover time of labile pool (1-exp(-(lambda * theta)^k))
      ! tuned to reduce growth from 20% saturation to zero at 0% saturation
      lambda_labile_theta,  & !< lambda in moisture response function of labile pool 
      k_labile_theta,       & !< k in moisture response function of labile pool 
      ! pull from reserve to labile pool due to phenology (1-exp(-(lambda * lai/target_lai)^k))
      lambda_phiphen,       & !< lambda in storage response function to phenology 
      k_phiphen,            & !< k in storage response function to phenology 
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phimaint,      & !< lambda in storage response function to maintenance needs 
      k_phimaint,           & !< k in storage response function to maintenance needs
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phistore,      & !< lambda in labile response function to storage size
      k_phistore,           & !< k in labile response function to storage size 
      ! sink limitation signal strength to photosynthesis (exp(-(lambda * labile/target_labile)^k))
      lambda_sinklim_ps,    & !< lambda in photosynthesis response function to labile pool size 
      k_sinklim_ps,         & !< k in photosynthesis response function to labile pool size 
      beta_sinklim_ps_min     !< maximum downacclimation of photosynthetic capacity due to sink limitation 

  ! Resorption parameters
  REAL(wp), SAVE :: & 
      resorp_fract_leaf,    & !< default fraction of nutrient resporption before leaf shedding
      resorp_fract_wood,    & !< default fraction of nutrient resporption before wood death
      root2leaf_cn,         & !< relative C:N of fine roots compared to leaves (mol C / mol N)
      wood2leaf_cn,         & !< relative C:N of sapwood biomass compared to leaves (mol C / mol N)
      root2leaf_np,         & !< relative N:P of fine roots compared to leaves (mol N / mol P)
      wood2leaf_np,         & !< relative N:P of sapwood biomass compared to leaves (mol N / mol P)
      tau_nutrient_recycling  !< time scale of enzyme, RNA etc. turnover (as in return to labile pool!, days) 


  ! Parameters concerning N uptake calculations
  REAL(wp), SAVE :: &
      ! maximum cost of N fixation at low / high temperatures (mol C / mol N) 
      cost_n_fixation_max  , &
      ! mimumum cost of N fixation at 25°C (mol C / mol N)
      cost_n_fixation_min  , &
      ! factors in temperature response function of nitrogenase activity
      k1_temp_bnf          , &
      k2_temp_bnf          , &
      k3_temp_bnf          , &
      k4_temp_bnf          , &
      ! N transformation costs of converting NH4 to organic N  (mol C / mol N)
      transform_cost_nh4   , &
      ! N transformation costs of converting NO3 to organic N (mol C / mol N)
      transform_cost_no3   , &
      ! N transformation costs of converting organic N to organic N (mol C / mol N)
      transform_cost_org   , &
      ! fractions of NO3 and NH4 uptake based on transformation costs
      f_nacq_no3           , &
      f_nacq_nh4           , &
      ! maximum rate of N fixation (µmol N / mol C fine root / s )
      vmax_bnf             , & 
      ! half_saturating cost of the response of N fixation (mol C / mol N) 
      k_bnf                , &
      ! factor for BNF response in resistance approach (mol C / mol N)
      offset_bnf           , &
      ! maximum soil N at which BNF occurs in standard scheme (mol N / m2)
      soil_n_threshold_bnf

  ! parameters used for phenology calculations
  REAL(wp), SAVE :: & 
      min_lai  ,                 & !< minimum LAI beyond which all leaves are shed
      target_lai_max  ,          & !< maximum LAI target for reserve use calculations 
      max_leaf_shedding_rate       !< maximum rate of leaf shedding (1/days) 

  ! vegetation dynamics
  REAL(wp), SAVE :: & 
      background_mort_rate_tree  ,     &  !< background mortality rate for trees in static mortality case (1/year)
      background_mort_rate_grass  ,     & !< background mortality rate for grass in static mortality case (1/year)
      lambda_mort_light   ,     & !< parameter in the Weibull function controlling self-thinning mortality 
      k_mort_light         ,     & !< parameter in the Weibull function controlling self-thinning mortality 
      lambda_est_temp      ,     & !< parameter in the Weibull function controlling temperature-limited establishment 
      k_est_temp           ,     & !< parameter in the Weibull function controlling temperature-limited establishment 
      lambda_est_moist     ,     & !< parameter in the Weibull function controlling moisture-limited establishment 
      k_est_moist          ,     & !< parameter in the Weibull function controlling moisture-limited establishment 
      ! these parameters are more canopy process related, but only useful for veg-dynamics?
      toc_rfr_ratio        ,     &               !< top of the canopy red:farred ratio  (Kull & Kruij 1998)
      k_r2fr_chl                   !< extinction of red:farred ratio on leaf chlorophyll (1/mmol N); 
                                                         !! (Kull & Krujit 1998 for the far red value, tuned to get a R:FR ratio
                                                          !! of 0.15 at full canopy closure (a common value) 

  ! other
  REAL(wp), SAVE :: & 
      ! carbon content of leaves per dry matter
      carbon_per_dryweight_leaf   ,     &   !< carbon mass per unit dry weight of leaves (gC / g DW) 
      ! minimum permissable height of an individual
      min_height  ,     & !< minimum permissable height of an individual
      ! stem mass to leaf mass ratio of grasses
      sm2lm_grass  ,     & !< stem mass to leaf mass ratio of grasses
      ! Parameters controlling fruit allocation in: exp(-(lambda * reserve_use+offset)^k) 
      lambda_fruit_alloc  ,     & !< parameter in the Weibull function controlling fruit growth
      k4_fruit_alloc  ,     & !< parameter in the Weibull function controlling fruit growth
      k3_fruit_alloc   ,     &  !< parameter in the Weibull function controlling fruit growth; ~µmol C / m2 / s
      ! minimum fraction of allocation going to fruit (flowers etc.)
      k1_fruit_alloc            !< minimum fraction of allocation going to fruit (flowers etc.)

  ! parameters for changes in canopy N fractions
  REAL(wp), SAVE :: &     
      ! maximum change in photosynthetic N fractions
      delta_n_fraction  ,     & !< maximum fractional change for optimal adjustment of photosynthetic N fraction
      ! maximum change in CN ratio of canopy
      delta_n_leaf      !< maximum fractional change for optimal adjsutment of canopy N allocation


  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  ! need to be parameters (i.e. constant across model runtime)
      INTEGER, PARAMETER ::    &
      itree         = 1  ,     &
      igrass        = 2  


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_veg_constants'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the vegetation process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_veg_constants

    USE mo_base_constants,      ONLY: one_day, one_year 
    USE mo_shared_constants,    ONLY: molar_mass_C, molar_mass_N

    IMPLICIT NONE
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_veg_constants'


  ! Respiration parameters
      fresp_growth      = 0.25_wp 
      fmaint_rate_base  = 1.0_wp 
      t_maint_k1        = 308.56_wp  
      t_maint_k2        = 56.02_wp 
      t_maint_k3        = 227.13_wp 
      f_resp_acclim     = -0.008_wp
      t_acclim_zero     = 283.15_wp 

  ! Labile pool parameters 
      tau_labile              = 5.0_wp  
      lambda_labile_temp      = 0.5_wp 
      k_labile_temp           = 2.0_wp 
      lambda_labile_theta     = 10.0_wp
      k_labile_theta          = 2.0_wp 
      lambda_phiphen          = 1.8_wp 
      k_phiphen               = 4.0_wp 
      lambda_phimaint         = 4.0_wp 
      k_phimaint              = 1.2_wp 
      lambda_phistore         = 2._wp
      k_phistore              = 3._wp 
      lambda_sinklim_ps       = 0.05_wp
      k_sinklim_ps            = 2.0_wp 
      beta_sinklim_ps_min     = 0.25_wp 

  ! Resorption parameters
      resorp_fract_leaf     = 0.5_wp  
      resorp_fract_wood     = 0.2_wp 
      root2leaf_cn          = 0.85_wp 
      wood2leaf_cn          = 0.145_wp 
      root2leaf_np          = 1.0_wp 
      wood2leaf_np          = 1.0_wp 

  ! N uptake parameters
      cost_n_fixation_max = 17.5_wp   / (molar_mass_C/molar_mass_N) 
      cost_n_fixation_min = 8.0_wp    / (molar_mass_C/molar_mass_N) 
      k1_temp_bnf         = 1.25_wp
      k2_temp_bnf         = -3.62_wp
      k3_temp_bnf         = 0.27_wp
      k4_temp_bnf         = 50.3_wp
      transform_cost_nh4  = 1.7_wp    / (molar_mass_C/molar_mass_N)
      transform_cost_no3  = 2.3_wp    / (molar_mass_C/molar_mass_N)
      transform_cost_org  = 1.0_wp    / (molar_mass_C/molar_mass_N)
      f_nacq_nh4          = transform_cost_no3 / (transform_cost_no3 + transform_cost_nh4)
      f_nacq_no3          = transform_cost_nh4 / (transform_cost_no3 + transform_cost_nh4)
      vmax_bnf            = 5.0_wp * 0.0225_wp / (molar_mass_N/molar_mass_C) / (one_day*one_year) *1.e6_wp
      k_bnf               = 50.0_wp   / (molar_mass_C/molar_Mass_N)
      offset_bnf          = 6.0_wp * (f_nacq_nh4 * transform_cost_nh4 + f_nacq_no3 * transform_cost_no3)
      soil_n_threshold_bnf = 0.04 ! mol N == 0.56 g N /m2

  ! parameters used for phenology calculations
      min_lai                     = 0.1_wp
      target_lai_max              = 5.0_wp 
      max_leaf_shedding_rate      = 0.05_wp

  ! vegetation dynamics
      background_mort_rate_tree   = 0.01_wp 
      background_mort_rate_grass  = 0.01_wp 
      lambda_mort_light  = 0.15_wp 
      k_mort_light       = 4.0_wp  
      lambda_est_temp    = 0.075_wp  
      k_est_temp         = 4.0_wp  
      lambda_est_moist   = 10._wp  
      k_est_moist        = 2.0_wp  
      toc_rfr_ratio      = 1.2_wp  
      k_r2fr_chl         = 0.0012_wp - 0.012_wp  

  ! other
      carbon_per_dryweight_leaf   = 0.48_wp  
      min_height                  = 0.1_wp 
      sm2lm_grass                 = 0.05_wp 
      lambda_fruit_alloc          = 10.0_wp 
      k1_fruit_alloc              = 0.01_wp 
      k3_fruit_alloc              = 0.1_wp  
      k4_fruit_alloc              = 2.0_wp 

  ! moved here from mo_assimi_constants
      tau_nutrient_recycling      = 10._wp 

  ! parameters for changes in canopy N fractions
!@todo:should this be 1/20days?
      delta_n_fraction  = 0.048_wp ! 0.001_wp (changed to account for multiplication with dtime/one_day in the code)
      delta_n_leaf      = 0.048_wp  ! 0.0001 (changed to account for multiplication with dtime/one_day in the code)   

  END SUBROUTINE init_veg_constants

END MODULE mo_veg_constants
