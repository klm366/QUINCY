!>
!! @brief Contains the routines for plant parameters 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_parameters

  USE mo_kind
  USE mo_base_constants, ONLY: npft

  ! Respiration parameters
  REAL(wp),SAVE :: fresp_growth                 ! growth respiration per unit new biomass 
  REAL(wp),SAVE :: fmaint_rate_base             ! N-specific maintenance respiration rate (µmol CO2 / mmol N / s) 
  REAL(wp),SAVE :: t_maint_k1                   ! Coefficient in the equation for temperature sensitivity of foliar respiration 
  REAL(wp),SAVE :: t_maint_k2                   ! Coefficient in the equation for temperature sensitivity of foliar respiration
  REAL(wp),SAVE :: t_maint_k3                   ! Coefficient in the equation for temperature sensitivity of foliar respiration 

  ! Labile pool parameters 
  REAL(wp),SAVE :: tau_labile_default           ! default turnover time of the labile carbo-hydrates (in days)
  REAL(wp),SAVE :: lambda_tau_labile_temp       ! lambda in temperature response function of labile pool 
  REAL(wp),SAVE :: weibull_k_tau_labile_temp    ! k in temperature response function of labile pool 
  REAL(wp),SAVE :: lambda_tau_labile_theta      ! lambda in moisture response function of labile pool 
  REAL(wp),SAVE :: weibull_k_tau_labile_theta   ! k in moisture response function of labile pool 
  REAL(wp),SAVE :: lambda_phiphen               ! lambda in storage response function to phenology 
  REAL(wp),SAVE :: weibull_k_phiphen            ! k in in storage response function to phenology 
  REAL(wp),SAVE :: lambda_phimaint              ! lambda in storage response function to maintenance needs 
  REAL(wp),SAVE :: weibull_k_phimaint           ! k in in storage response function to maintenance needs
  REAL(wp),SAVE :: lambda_phistore              ! lambda in storage response function to storage size
  REAL(wp),SAVE :: weibull_k_phistore           ! k in in storage response function to storage size 

  ! Resporption parameters
  REAL(wp),SAVE :: resorption_default           ! default fraction of nutrient resporption before leaf shedding
  REAL(wp),SAVE :: resorption_wood              ! default fraction of nutrient resporption before wood death
  REAL(wp),SAVE :: root2leaf_cn                 ! relative C:N of fine roots compared to leaves
  REAL(wp),SAVE :: wood2leaf_cn                 ! relative C:N of woody biomass compared to leaves
  REAL(wp),SAVE :: root2leaf_np                 ! relative N of fine roots compared to leaves
  REAL(wp),SAVE :: wood2leaf_np                 ! relative N:P of woody biomass compared to leaves

  ! Parameters used in the photosynthesis routine
  REAL(wp),SAVE :: carbon_per_dryweight_leaf    ! carbon mass per unit dry weight of leaves (gC / g DW) 
  REAL(wp),SAVE :: jmax2n                       ! electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
  REAL(wp),SAVE :: vcmax2n                      ! Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
  REAL(wp),SAVE :: chl2n                        ! moles chlorophyll per mol of N
  REAL(wp),SAVE :: k1_fnl_struc                 ! slope of the fraction of leaf N not associated with PS against leaf N (--/mmolN) 
  REAL(wp),SAVE :: jmax2vcmax                   ! ratio of Jmax25/Vcmax25 for prescribed leaf N 
  REAL(wp),SAVE :: pO2                          ! Partial Pressure of O2 in kPa
  REAL(wp),SAVE :: E0kc                         ! Baseline of Michaelis-Menten constant of Rubisco for CO2 (unit?)
  REAL(wp),SAVE :: E1kc                         ! Temperature sensitivity of Michaelis-Menten constant of Rubisco for CO2 (unit?)
  REAL(wp),SAVE :: E0ko                         ! Baseline of Michealis-Menten Kintetics for Rubisco-oxigination
  REAL(wp),SAVE :: E1ko                         ! Temperature sensitivity of Michealis-Menten Kintetics for Rubisco-oxigination
  REAL(wp),SAVE :: E0pcp                        ! Baseline of photosynthetic compensation point
  REAL(wp),SAVE :: E1pcp                        ! Temperature sensitivity of photosynthetic compensation point
  REAL(wp),SAVE :: E0v                          ! Baseline of Rubisco
  REAL(wp),SAVE :: E1v                          ! Temperature sensitivity of Rubisco
  REAL(wp),SAVE :: Tref_pep                     ! Temperature sensitivity of PeP C4 photosynthesis 
  REAL(wp),SAVE :: Tbase_pep                    ! Temperature sensitivity of PeP C4 photosynthesis
  REAL(wp),SAVE :: alpha_i                      ! intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I))
                                                ! µmol CO2 / mol quanta) 
  REAL(wp),SAVE :: CiCa_min                     ! allowable minimum leaf Ci:Ca for numerical stability
  REAL(wp),SAVE :: CiCa_max                     ! allowable maximum leaf Ci:Ca for numerical stability
  REAL(wp),SAVE :: CiCa_default_C3              ! default Ci:Ca for C3 plants
  REAL(wp),SAVE :: CiCa_default_C4              ! default Ci:Ca for C4 plants
  REAL(wp),SAVE :: ps_it_max                    ! maximum iteration of PS calculation
  
  ! parameters for canopy profile and canopy light extinction
  REAL(wp),SAVE :: kn                           ! extinction coefficient to describe decline of N within the canopy 
  REAL(wp),SAVE :: k0_fnl_chl                   ! parameters to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp),SAVE :: k1_fnl_chl                   ! parameters to describe the varying fraction of N in chlorophyll through the canopy
  REAL(wp),SAVE :: kfnl_chl                     ! parameters to describe the varying fraction of N in chlorophyll through the canopy

  ! Parameters required for albedo and fapar calculation
  REAL(wp),SAVE :: ka                           ! Extinction coefficient for PAR on chlorophyll
  REAL(wp),SAVE :: sigma                        ! single leaf scattering albedo 
  REAL(wp),SAVE :: kbl0                         ! extinction coefficient over black leaves 
  REAL(wp),SAVE :: kdf0                         ! extinction coefficient for diffuse PAR (Spitters eq. 7) 
  REAL(wp),SAVE :: rho2sbeta                    ! scaling factor of solar angle in reflection calculation (mo_veg_fapar) 

  ! some useful numbers for easy identification of photosynthetic pathways and growth forms
  INTEGER,SAVE :: itree = 1 
  INTEGER,SAVE :: igrass = 2 
  INTEGER,SAVE :: ic3phot = 1
  INTEGER,SAVE :: ic4phot = 2 

  TYPE VegetationParameterClass
     INTEGER,DIMENSION(npft) :: growthform 
     INTEGER,DIMENSION(npft) :: ps_pathway
     REAL(wp),DIMENSION(npft) :: tau_leaf
     REAL(wp),DIMENSION(npft) :: sla
     REAL(wp),DIMENSION(npft) :: cn_leaf
     REAL(wp),DIMENSION(npft) :: np_leaf
     REAL(wp),DIMENSION(npft) :: k0_fnl_struc
     REAL(wp),DIMENSION(npft) :: nf 
     REAL(wp),DIMENSION(npft) :: tjmax_opt
     REAL(wp),DIMENSION(npft) :: tjmax_omega
     REAL(wp),DIMENSION(npft) :: g0
     REAL(wp),DIMENSION(npft) :: g1
     REAL(wp),DIMENSION(npft) :: lambda_soil_ps
     REAL(wp),DIMENSION(npft) :: weibull_k_soil_ps
     REAL(wp),DIMENSION(npft) :: lambda_soil_gs
     REAL(wp),DIMENSION(npft) :: weibull_k_soil_gs
     REAL(wp),DIMENSION(npft) :: tau_sap_wood
     REAL(wp),DIMENSION(npft) :: wood_density 
     REAL(wp),DIMENSION(npft) :: k_latosa 
     REAL(wp),DIMENSION(npft) :: k_rtos 
     REAL(wp),DIMENSION(npft) :: tau_fine_root
     REAL(wp),DIMENSION(npft) :: tau_coarse_root
     REAL(wp),DIMENSION(npft) :: allom_k1
     REAL(wp),DIMENSION(npft) :: allom_k2
     REAL(wp),DIMENSION(npft) :: tau_fruit
     REAL(wp),DIMENSION(npft) :: tau_seed
  END TYPE VegetationParameterClass

  TYPE AccumulationIntervalPlantClass 
    ! interval over which the values are accumulated/ averaged
    ! interval in days !!
    ! see mo_base_memory for defintion of TYPEs that keep the starting/default value
    REAL(wp) :: par_sunlit   = 1._wp
    REAL(wp) :: par_shaded   = 1._wp
    REAL(wp) :: fleaf_sunlit = 1._wp
    REAL(wp) :: beta_air     = 1._wp
    REAL(wp) :: beta_soil_ps = 1._wp
    REAL(wp) :: beta_soil_gs = 1._wp

    REAL(wp) :: weekly       = 7._wp         ! 7 days; constant across the model
  END TYPE AccumulationIntervalPlantClass


  PUBLIC :: init_veg_parameters

  CONTAINS 

    SUBROUTINE init_veg_parameters(veg_param)
   
      USE mo_land_constants
      TYPE(VegetationParameterClass),INTENT(inout) :: veg_param

      !IMPLICIT NONE
      ! 
      ! this subroutine sets global constants to describe the vegetation processes
      ! defined in the module header
      !   
      ! growth respiration per unit new biomass
      fresp_growth = 0.25_wp
      ! N-specific foliar maintenance respiration rate (µmol CO2 / mol N / s)
      ! Taken from Sprugel et al. Respiration from the organ level to the stand. 
      ! In: Smith WK, Hinckley, T.M., editor. Physiological Ecology of Coniferous Forests.
      ! San Diego, CA, USA: Academic Press; 1996.
      ! Note converted from 7.4 e-7 to 1.25 e-6 to account for the different values of the Lloyd&Taylor equation
      ! such that the applied rates a comparable at 15-20°C
      fmaint_rate_base = 1.25_wp 
      ! parameters in the Lloyd and Taylor temperature response function (Lloyd & Taylor, 1996, Funct. Ecol.) 
      t_maint_k1 = 308.56_wp  
      t_maint_k2 = 56.02_wp
      t_maint_k3 = 227.13_wp

      ! default turnover time of the labile carbo-hydrates (in days)
      tau_labile_default = 0.1_wp 
      ! temperature scalar on turnover time of labile pool (1-exp(-(lambda * T)^k))
      ! tuned to reduce growth between 5°C to zero at 0°C
      lambda_tau_labile_temp = 0.5_wp
      weibull_k_tau_labile_temp = 2.0_wp
      ! moisture scalar on turnover time of labile pool (1-exp(-(lambda * theta)^k))
      ! tuned to reduce growth from 20% saturation to zero at 0% saturation
      lambda_tau_labile_theta = 0.1_wp
      weibull_k_tau_labile_theta = 2.0_wp
      ! pull from reserve to labile pool due to phenology (1-exp(-(lambda * lai/target_lai)^k))
      lambda_phiphen = 1.8_wp
      weibull_k_phiphen = 4.0_wp
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phimaint = 4.0_wp
      weibull_k_phimaint = 1.2_wp
      ! pull from reserve to labile pool due to maintenance of labile pool (1-exp(-(lambda * labile/target_labile)^k))
      lambda_phistore = 2._wp
      weibull_k_phistore = 3._wp

      ! default fraction of nutrient resporption before leaf shedding
      resorption_default = 0.5_wp 
      ! default fraction of nutrient resporption before wood death 
      resorption_wood = 0.2_wp 
      ! ratio of root C:N to leaf C:N
      root2leaf_cn = 0.85_wp 
      ! ratio of root C:N to leaf C:N
      wood2leaf_cn = 0.145_wp 
      ! ratio of root N:P to leaf N:P
      root2leaf_np = 1.0_wp 
      ! ratio of root N:P to leaf N:P
      wood2leaf_np = 1.0_wp 

      ! Parameters used in the photosynthesis routine
      ! Partial Pressure of O2 in kPa
      pO2 = 20.9_wp
      ! Temperature sensitivity and baseline of Michaelis-Menten constant of Rubisco for CO2 (unit?)
      E0kc = 38.05_wp
      E1kc = 79430.0_wp
      ! Temperature sensitivity and baseline of Michealis-Menten Kintetics for Rubisco-oxigination     
      E0ko = 20.3_wp
      E1ko = 36380.0_wp
      ! Temperature sensitivity and baseline of photosynthetic compensation point
      E0pcp = 19.02_wp
      E1pcp = 37830.0_wp
      ! Temperature sensitivity and baseline of Rubisco
      E0v = 26.35_wp
      E1v = 65330.0_wp
      ! Temperature sensitivity and baseline value of foliar respiration
      E0kt = 18.72_wp
      E1kt = 46390.0_wp
      ! Temperature sensitivity of PeP C4 photosynthesis
      Tref_pep = 25._wp
      Tbase_pep = 10._wp
      ! intrinsic quantum yield (efficiency of the quanta absorption of PS II (or PS I)) 
      ! µmol CO2 / mol quanta)
      alpha_i = 0.08_wp
      ! initial guess of the leaf Ci to Ca ratio
      CiCa_default_C3 = 0.7_wp
      CiCa_default_C4 = 0.95_wp
      ! maximum number of iteration in photosynthesis calculation (should probably be in mo_base_constants)
      ps_it_max = 10

      ! carbon content of leaves per dry matter
      carbon_per_dryweight_leaf = 0.48_wp ! g C / g DW
      ! Parameters describing the PS ~ An relationship
      ! Values from Niinemets 1997 (SC 01/16)
      jmax2n=17.6_wp/4.0_wp   ! electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
      vcmax2n=1.8_wp          ! Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
      ! values from Evans 1988
      chl2n=1000._wp/25.5_wp  ! moles chlorophyll per mmol N
      ! slope of the relationship of structural leaf N against leaf N in fraction/mmol (Friend 1997)
      k1_fnl_struc=71.4_wp * molar_mass_N / 1000000._wp
      ! in static leaf N mode: prescribed ratio of Jmax25/Vcmax25
      jmax2vcmax=1.92_wp      ! text book
      ! parameters describing the canopy profile of N 
      !  extinction coefficient to describe decline of N within the canopy 
      kn = 0.11_wp
      ! parameters to describe the varying fraction of N in chlorophyll through the canopy
      k0_fnl_chl = 6.0_wp
      k1_fnl_chl = 3.6_wp
      kfnl_chl = 0.7_wp

      ! Parameters required for albedo and fapar calculation
      ! single leaf scattering albedo
      sigma = 0.2_wp
      ! extinction coefficient over black leaves
      kbl0 = 0.5_wp
      ! extinction coefficient for diffuse PAR (Spitters eq. 7)
      kdf0 = 0.8_wp
      ! scaling factor of solar angle in reflection calculation (mo_veg_fapar)
      rho2sbeta = 1.6_wp
      ! Extinction coefficient for PAR on chlorophyll
      ka = 0.005_wp

      CALL init_pft_param_table(veg_param)

    END SUBROUTINE init_veg_parameters

    SUBROUTINE init_pft_param_table(veg_param)
      ! 
      ! This subroutine initalise plant type specific parameters
      ! This should obviously not be hard coded, but read from a file
      !
      USE mo_land_constants, ONLY: molar_mass_C,molar_mass_N,molar_mass_P
      IMPLICIT NONE
      TYPE(VegetationParameterClass),INTENT(inout) :: veg_param
 
      ! types are 
      ! 1) broad leaved evergreen (rainforest) 
      ! 2) broad leaved evergreen (xeric forest) 
      ! 3) broad leaved deciduous (rain green) 
      ! 4) broad leaved deciduous (summer green) 
      ! 5) needle-leaved evergreen
      ! 6) needle-leaved deciduous 
      ! 7) C3 grass
      ! 8) C4 grass
      ! vegetation growth form
      veg_param%growthform(:) = & 
      (/itree,itree,itree,itree,itree,itree,igrass,igrass/)
      ! photosynthetic path way 
      veg_param%ps_pathway(:) = & 
      (/ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic4phot/)
      ! average turnover (1/lifetime) of a leaf (1/years) 
      ! From Kattge et al. 2011, Table 5, including conversion from month lifetime to annual turnover
      veg_param%tau_leaf(:) = 12._wp / & 
      (/16.83_wp,15.88_wp,5.83_wp,5.83_wp,39.7_wp,6.08_wp,3.85_wp,3.85_wp/)
      ! specific leaf area of a plant (m2 / mol C) 
      ! From Kattge et al. 2011, Table 5, including conversion from mm/mgDW to m2/mol
      veg_param%sla(:) = 1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * &
      (/9.46_wp,8.99_wp,15.4_wp,15.4_wp,5.0_wp,10.09_wp,20.12_wp,19.23_wp/) 
      ! default leaf C:N (molC/molN) if not calculated dynamically 
      ! From Kattge et al. 2011, Table 5, including conversion from mg/g to mol/mol
      veg_param%cn_leaf(:) = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / & 
      (/16.89_wp,13.73_wp,21.32_wp,21.32_wp,12.09_wp,19.37_wp,17.84_wp,14.14_wp/) 
      ! default leaf N:P (molN/molP) if not calculated dynamically 
      ! From Kattge et al. 2011, Table 5, including conversion from mg/g to mol/mol
      veg_param%np_leaf(:) = carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / veg_param%cn_leaf(:) / &
      (/0.86_wp,0.84_wp,1.44_wp,1.44_wp,1.23_wp,1.83_wp,1.43_wp,1.36_wp/)
      ! fraction of leaf N not used for photosynthesis (Friend 1997) 
      ! temporarily tuned to match Amax/N in Kattge etal. 2012, Table 5
      veg_param%k0_fnl_struc(:) = & 
      !(/0.67_wp,0.67_wp,0.67_wp,0.67_wp,0.83_wp,0.83_wp,0.67_wp,0.67_wp/) 
      (/0.79_wp,0.73_wp,0.67_wp,0.67_wp,0.87_wp,0.83_wp,0.3_wp,0.67_wp/) 
      !Le Maire et al. 2012 (derived from grasses, does not give plausble results for needles)
      !0.012_wp / 0.48_wp * molar_mass_C / molar_mass_N / veg_param%sla(:)/&
      !(1._wp / veg_param%sla(:) / veg_param%cn_leaf(:)) 
      ! nitrogen scaling parameter for photosynthesis calculation 
      veg_param%nf(:) = 1.83_wp
      ! initial value of the temperature optimum of electron transport and its shape parameter (deg C)
      veg_param%tjmax_opt(:) = 17._wp
      veg_param%tjmax_omega(:) = 18._wp
      ! intercept and slope of the An ~gs relationship (mmol/µmol)
      veg_param%g0(:) = & ! 0.01_wp
      (/0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.03_wp,0.03_wp/)
      veg_param%g1(:) = & !9.3_wp
      ! g1 in Ball & Berry relationship, tuned to get similar Amax/N values as with Lin et al. g1 
      (/9.3_wp,8.3_wp,7.0_wp,10.9_wp,5.5_wp,7.0_wp,9.3_wp,1.62_wp/)
      ! g1 in kPa^0.5 as in Lin et al. 2015, Nat CC. OBS: no values for Larches - assigned broad leaved value!
      !(/3.77_wp,3.37_wp,2.98_wp,4.64_wp,2.35_wp,4.64_wp,4.16_wp,1.62_wp/)
      ! Weibull parameters describing the response of leaf-level photosynthesis to 
      ! fractional water availability
      veg_param%lambda_soil_ps(:) = 3.0_wp
      veg_param%weibull_k_soil_ps(:) = 3.0_wp
      ! Weibull parameters describing the response of stomatal conductance to 
      ! fractional water availability
      veg_param%lambda_soil_gs(:) = veg_param%lambda_soil_ps(:)
      veg_param%weibull_k_soil_gs(:) = veg_param%weibull_k_soil_ps(:)
      ! turnover time of the sapwood pool (1/years)
      veg_param%tau_sap_wood(:) = 0.05_wp
      ! wood density (mol C / m3)
      veg_param%wood_density(:) = 30000.0_wp / molar_mass_C
      ! leaf area to sapwood area ratio 
      veg_param%k_latosa(:) = 4000.0_wp 
      ! trade-off parametr for hydraulic investment into sapwood or fine roots 
      veg_param%k_rtos(:) = 1.0_wp / (veg_param%k_latosa(:) / veg_param%sla(:) / veg_param%wood_density(:))
      ! average turnover (1/lifetime) of a fine root (1/years) 
      veg_param%tau_fine_root(:) = 0.7_wp 
      ! average turnover (1/lifetime) of a coarse root (1/years) 
      veg_param%tau_coarse_root(:) = 0.1_wp 
      ! paramter in height diameter relationship 
      veg_param%allom_k1(:) = 55.0_wp 
      ! paramter in height diameter relationship
      veg_param%allom_k2(:) = 0.65 
      ! turnover time of the fruit pool (1/years)
      veg_param%tau_fruit(:) = 1.0_wp
      ! turnover time of the fruit pool (1/years)
      veg_param%tau_seed(:) = 1.0_wp

    END SUBROUTINE init_pft_param_table 

END MODULE mo_veg_parameters
