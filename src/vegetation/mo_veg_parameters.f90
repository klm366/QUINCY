!!
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
  USE mo_base_constants, ONLY : npft

  ! Respiration parameters
  REAL(wp),SAVE :: fresp_growth                 ! growth respiration per unit new biomass 
  REAL(wp),SAVE :: fmaint_rate_base             ! N-specific maintenance respiration rate (µmol CO2 / mmol N / s) 
  REAL(wp),SAVE :: t_maint_k1                   ! Coefficient in the equation for temperature sensitivity of foliar respiration 
  REAL(wp),SAVE :: t_maint_k2                   ! Coefficient in the equation for temperature sensitivity of foliar respiration
  REAL(wp),SAVE :: t_maint_k3                   ! Coefficient in the equation for temperature sensitivity of foliar respiration 

  ! Standard pool turnover times
  REAL(wp),SAVE :: tau_labile_default           ! default turnover time of the labile carbo-hydrates (in days)

  ! Resporption parameters
  REAL(wp),SAVE :: resorption_default           ! default fraction of nutrient resporption before leaf shedding
  REAL(wp),SAVE :: resorption_wood              ! default fraction of nutrient resporption before wood death

  ! Parameters used in the photosynthesis routine
  REAL(wp),SAVE :: jmax2n                       ! electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
  REAL(wp),SAVE :: vcmax2n                      ! Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
  REAL(wp),SAVE :: chl2n                        ! moles chlorophyll per mol of N
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
     REAL(wp),DIMENSION(npft) :: sla
     REAL(wp),DIMENSION(npft) :: nf 
     REAL(wp),DIMENSION(npft) :: tjmax_opt
     REAL(wp),DIMENSION(npft) :: tjmax_omega
     REAL(wp),DIMENSION(npft) :: g0
     REAL(wp),DIMENSION(npft) :: g1
     REAL(wp),DIMENSION(npft) :: lambda_soil_ps
     REAL(wp),DIMENSION(npft) :: weibull_k_soil_ps
     REAL(wp),DIMENSION(npft) :: lambda_soil_gs
     REAL(wp),DIMENSION(npft) :: weibull_k_soil_gs
     REAL(wp),DIMENSION(npft) :: tau_sapwood
     REAL(wp),DIMENSION(npft) :: wood_density 
     REAL(wp),DIMENSION(npft) :: k_latosa 
     REAL(wp),DIMENSION(npft) :: k_rtos 
     REAL(wp),DIMENSION(npft) :: allom_k1
     REAL(wp),DIMENSION(npft) :: allom_k2
     REAL(wp),DIMENSION(npft) :: tau_fruit_pool 
  END TYPE VegetationParameterClass

  TYPE AccumulationPeriodClass   ! new
    ! length of period over which the values are accumulated
    REAL(wp) :: light = 1.
    REAL(wp) :: npp = 1.
    REAL(wp) :: gpp = 1.
  END TYPE AccumulationPeriodClass


  PUBLIC :: init_veg_parameters,init_pft_param_table

  CONTAINS 

    SUBROUTINE init_veg_parameters
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
      ! Note converted from 7.4 e-7 to 1.48 e-6 to account for the different values of the Lloyd&Taylor equation
      fmaint_rate_base = 1.48_wp 
      ! parameters in the Lloyd and Taylor temperature response function (Lloyd & Taylor, 1996, Funct. Ecol.) 
      t_maint_k1 = 308.56_wp  
      t_maint_k2 = 56.02_wp
      t_maint_k3 = 227.13_wp

      ! default turnover time of the labile carbo-hydrates (in days)
      tau_labile_default = 0.1_wp 

      ! default fraction of nutrient resporption before leaf shedding
      resorption_default = 0.5_wp 
      ! default fraction of nutrient resporption before wood death 
      resorption_wood = 0.2_wp 

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

      ! Parameters describing the PS ~ An relationship
      ! OBS: CURRENTLY THESE VALUES ARE MULTIPLIED BY 10 FOR BACKWARDS COMPATABILITY
      ! NEEDS TO BE REPLACED BY REAL VALUES
      jmax2n=1.2_wp  ! electron-transport limited carboxylation rate per µmol CO2 / mmol N in ET
      vcmax2n=2.3_wp ! Rubisco limited carboxylation rate per µmol CO2 / mmol N in RUB
      chl2n=1000._wp/39.5_wp  ! moles chlorophyll per mol of N
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

    END SUBROUTINE init_veg_parameters

    SUBROUTINE init_pft_param_table(veg_param)
      ! 
      ! This subroutine initalise plant type specific parameters
      ! This should obviously not be hard coded, but read from a file
      !
      USE mo_land_constants, ONLY: molar_mass_C
      IMPLICIT NONE
      TYPE(VegetationParameterClass),INTENT(inout) :: veg_param
 
      ! vegetation growth form
      veg_param%growthform(:) = itree
      ! specific leaf area of a plant (m2 / mol C) 
      veg_param%sla(:) = 10.0_wp / 1000.0_wp * molar_mass_C 
      ! nitrogen scaling parameter for photosynthesis calculation 
      veg_param%nf(:) = 1.83_wp
      ! initial value of the temperature optimum of electron transport and its shape parameter (deg C)
      veg_param%tjmax_opt(:) = 17._wp
      veg_param%tjmax_omega(:) = 18._wp
      ! intercept and slope of the An ~gs relationship (mmol/µmol)
      veg_param%g0(:) = 0.01_wp
      veg_param%g1(:) = 9.3_wp
      ! Weibull parameters describing the response of leaf-level photosynthesis to 
      ! fractional water availability
      veg_param%lambda_soil_ps(:) = 3.0_wp
      veg_param%weibull_k_soil_ps(:) = 3.0_wp
      ! Weibull parameters describing the response of stomatal conductance to 
      ! fractional water availability
      veg_param%lambda_soil_gs(:) = veg_param%lambda_soil_ps(:)
      veg_param%weibull_k_soil_gs(:) = veg_param%weibull_k_soil_ps(:)
      ! turnover time of the sapwood pool (1/years)
      veg_param%tau_sapwood(:) = 0.05_wp
      ! wood density (mol C / m3)
      veg_param%wood_density(:) = 30000.0_wp / molar_mass_C
      ! leaf area to sapwood area ratio 
      veg_param%k_latosa(:) = 4000.0_wp 
      ! trade-off parametr for hydraulic investment into sapwood or fine roots 
      veg_param%k_rtos(:) = 1.0_wp / (veg_param%k_latosa(:) / veg_param%sla(:) / veg_param%wood_density(:))
      ! paramter in height diameter relationship 
      veg_param%allom_k1(:) = 55.0_wp 
      ! paramter in height diameter relationship
      veg_param%allom_k2(:) = 0.65 
      ! turnover time of the fruit pool (1/years)
      veg_param%tau_fruit_pool(:) = 1.0_wp

    END SUBROUTINE init_pft_param_table 

END MODULE mo_veg_parameters
