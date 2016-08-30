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
      t_maint_k3           !< Coefficient in the equation for temperature sensitivity of foliar respiration

  ! Labile pool parameters 
  REAL(wp), SAVE :: &
      ! default turnover time of the labile carbo-hydrates (in days)
      tau_labile_default   ,     &  !< default turnover time of the labile carbo-hydrates (in days)
      ! temperature scalar on turnover time of labile pool (1-exp(-(lambda * T)^k))
      ! tuned to reduce growth between 5°C to zero at 0°C
      lambda_tau_labile_temp  ,     & !< lambda in temperature response function of labile pool 
      weibull_k_tau_labile_temp  ,     & !< k in temperature response function of labile pool 
      ! moisture scalar on turnover time of labile pool (1-exp(-(lambda * theta)^k))
      ! tuned to reduce growth from 20% saturation to zero at 0% saturation
      lambda_tau_labile_theta  ,     & !< lambda in moisture response function of labile pool 
      weibull_k_tau_labile_theta  ,     & !< k in moisture response function of labile pool 
      ! pull from reserve to labile pool due to phenology (1-exp(-(lambda * lai/target_lai)^k))
      lambda_phiphen  ,     & !< lambda in storage response function to phenology 
      weibull_k_phiphen  ,     & !< k in in storage response function to phenology 
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phimaint  ,     & !< lambda in storage response function to maintenance needs 
      weibull_k_phimaint  ,     & !< k in in storage response function to maintenance needs
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phistore  ,     & !< lambda in storage response function to storage size
      weibull_k_phistore        !< k in in storage response function to storage size 

  ! Resporption parameters
  REAL(wp), SAVE :: & 
      ! default fraction of nutrient resporption before leaf shedding
      resorption_default   ,     & !< default fraction of nutrient resporption before leaf shedding
      ! default fraction of nutrient resporption before wood death 
      resorption_wood   ,     & !< default fraction of nutrient resporption before wood death
      ! ratio of root C:N to leaf C:N
      root2leaf_cn   ,     & !< relative C:N of fine roots compared to leaves
      ! ratio of root C:N to leaf C:N
      wood2leaf_cn   ,     & !< relative C:N of woody biomass compared to leaves
      ! ratio of root N:P to leaf N:P
      root2leaf_np   ,     & !< relative N of fine roots compared to leaves
      ! ratio of root N:P to leaf N:P
      wood2leaf_np           !< relative N:P of woody biomass compared to leaves

  ! parameters used for phenology calculations
  REAL(wp), SAVE :: & 
      ! minimum LAI for phenology calculations
      min_lai  ,     & !< minimum LAI beyond which all leaves are shed
      target_lai_max  ,     & !< maximum LAI target for reserve use calculations 
      ! minimum leaf age beyond which phenological senescence can set in (days) 
      min_leaf_age  ,     & !< minimum leaf age before senescence is permitted (days) 
      ! temperature threshold for growing degree days 
      gdd_t_air_threshold  ,     & !< temperature threshold for growing degree days 
      ! maximal GDD requirement (assuming no dormance) 
      gdd_req_max  ,     & !< maximum GDD requirement
      ! temperature threshold for growing degree days 
      k_gdd_dormance  ,     & !< exponent factor in GDD - nd_dormance equation 
      ! weekly temperature initiating leaf senescence
      t_air_senescence  ,     & !< weekly air temperature initiating leaf senescence
      ! value of beta_soil for leaf flushing for trees
      beta_soil_flush_grass  ,     & !< value of beta_soil at leaf flushing for grass
      ! value of beta_soil for leaf flushing for trees
      beta_soil_flush_tree  ,     & !< value of beta_soil at leaf flushing for trees 
      ! value of beta_soil at shedding of leaves for trees
      beta_soil_senescence_tree      !< value of beta_soil at leaf shedding for trees

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
      weibull_k_fruit_alloc  ,     & !< parameter in the Weibull function controlling fruit growth
      offset_fruit_alloc   ,     &  !< parameter in the Weibull function controlling fruit growth; ~µmol C / m2 / s
      ! minimum fraction of allocation going to fruit (flowers etc.)
      min_falloc_fruit            !< minimum fraction of allocation going to fruit (flowers etc.)

  ! moved here from mo_assimi_constants at about r204
  REAL(wp), SAVE :: & 
      tau_canopy_acclimation      !< time scale of acclimation of canopy processes (N, P, T) (days)

  ! soil pools  (may be moved to mo_soil_constants later ?)
  REAL(wp), SAVE :: & 
      ! soil litter turnover pool, s equivalent to 1 year
      tau_litter      !< soil litter pool turnover time

  ! parameters for changes in canopy N fractions
  REAL(wp), SAVE :: &     
      ! maximum change in photosynthetic N fractions
      delta_n_fraction  ,     & !< maximum fractional change for optimal adjustment of photosynthetic N fraction
      ! maximum change in CN ratio of canopy
      delta_n_leaf      !< maximum fractional change for optimal adjsutment of canopy N allocation


  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  ! need to be parameters (i.e. constant across model runtime)
	  INTEGER, PARAMETER ::      &
      itree         = 1  ,     &
      igrass        = 2  ,     &
      ic3phot       = 1  ,     &
      ic4phot       = 2  ,     &
      ievergreen    = 1  ,     &
      iraingreen    = 2  ,     &
      isummergreen  = 3  ,     &
      iperennial    = 4


CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize constants used in the vegetation process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_veg_constants

    USE mo_shared_constants,    ONLY: Tzero

    IMPLICIT NONE


  ! Respiration parameters
      fresp_growth      = 0.25_wp 
      fmaint_rate_base  = 1._wp 
      t_maint_k1        = 308.56_wp  
      t_maint_k2        = 56.02_wp 
      t_maint_k3        = 227.13_wp 

  ! Labile pool parameters 
      tau_labile_default          = 0.2_wp  
      lambda_tau_labile_temp      = 0.5_wp 
      weibull_k_tau_labile_temp   = 2.0_wp 
      lambda_tau_labile_theta     = 0.1_wp
      weibull_k_tau_labile_theta  = 2.0_wp 
      lambda_phiphen              = 1.8_wp 
      weibull_k_phiphen           = 4.0_wp 
      lambda_phimaint             = 4.0_wp 
      weibull_k_phimaint          = 1.2_wp 
      lambda_phistore             = 2._wp
      weibull_k_phistore          = 3._wp 

  ! Resporption parameters
      resorption_default    = 0.5_wp  
      resorption_wood       = 0.2_wp 
      root2leaf_cn          = 0.85_wp 
      wood2leaf_cn          = 0.145_wp 
      root2leaf_np          = 1.0_wp 
      wood2leaf_np          = 1.0_wp 

  ! parameters used for phenology calculations
      min_lai                     = 0.1_wp
      target_lai_max              = 5.0_wp 
      min_leaf_age                = 10.0_wp 
      gdd_t_air_threshold         = 5.0_wp + Tzero 
      gdd_req_max                 = 800._wp 
      k_gdd_dormance              = 0.0098_wp 
      t_air_senescence            = 5.0_wp + Tzero 
      beta_soil_flush_grass       = 0.9_wp 
      beta_soil_flush_tree        = 0.98_wp 
      beta_soil_senescence_tree   = 0.5_wp

  ! vegetation dynamics
      background_mort_rate_tree   = 0.01_wp 
      background_mort_rate_grass  = 0.01_wp 
      lambda_mort_light  = 0.15_wp 
      k_mort_light       = 4.0_wp  
      lambda_est_temp    = 0.1_wp  
      k_est_temp         = 4.0_wp  
      lambda_est_moist   = 0.1_wp  
      k_est_moist        = 2.0_wp  
      toc_rfr_ratio      = 1.2_wp  
      k_r2fr_chl         = 0.0012_wp - 0.012_wp  

  ! other
      carbon_per_dryweight_leaf   = 0.48_wp  
      min_height                  = 0.1_wp 
      sm2lm_grass                 = 0.05_wp 
      lambda_fruit_alloc          = 10.0_wp 
      weibull_k_fruit_alloc       = 2.0_wp 
      offset_fruit_alloc          = 0.1_wp  
      min_falloc_fruit            = 0.01_wp 

  ! moved here from mo_assimi_constants
      tau_canopy_acclimation      = 10._wp 

  ! soil pools (may move to mo_soil_constants later) !
      tau_litter                  = 1.0_wp * 3600.0_wp * 24.0_wp * 365.0_wp  

  ! parameters for changes in canopy N fractions
!@todo:should this be 1/20days?
      delta_n_fraction  = 0.048 ! 0.001_wp (changed to account for multiplication with dt/one_day in the code)
      delta_n_leaf      = 0.0048_wp  ! 0.0001 (changed to account for multiplication with dt/one_day in the code)   

  END SUBROUTINE init_veg_constants

END MODULE mo_veg_constants
