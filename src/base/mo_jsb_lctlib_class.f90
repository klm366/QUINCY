!>
!! @brief 
!!   Contains the TYPE that holds all parameters and information on the landcover data considering lct "vegetated" \n
!!   Contains the lctlib and a wrapper for t_slm_parameter -> t_lctlib_element
!!
!! @author
!!  Jan Engel, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Jan Engel (2016-06-13)
!! 2nd version by Jan Engel (2016-07-27): renamed mo_base_slm_parameter_class.f90 -> mo_jsb_lctlib_class.f90
!!
MODULE mo_jsb_lctlib_class

  USE mo_kind


  IMPLICIT NONE
  PRIVATE
  PUBLIC :: init_lctlib, t_lctlib_element

  !-----------------------------------------------------------------------------------------------------
  !> vector holding all parameters necessary for one pft
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE t_lctlib_element
    ! this is JSBACH4 stuff
    !INTEGER  :: LctNumber            !! Unique number for each landcover type (1, 2, 3, ..)
    !CHARACTER(LEN=16)  :: LctName    !! Names of landcover types
    !INTEGER  :: LandcoverClass       !! Landcover classes (not to be confused with landcover types!):
    ! ----

     INTEGER :: growthform 
     INTEGER :: ps_pathway
     INTEGER :: phenology_type
     REAL(wp) :: sla                         !< specific leaf area (m2 / mol C)
     REAL(wp) :: sigma_vis 
     REAL(wp) :: sigma_nir 
     REAL(wp) :: omega_clumping 
     REAL(wp) :: crown_shape_factor 
     REAL(wp) :: cn_leaf                    !< default C:N of leaves (mol C / mol N)
     REAL(wp) :: np_leaf                    !< default N:P of leaves (mol N / mol P)
     REAL(wp) :: k0_fn_struc
     REAL(wp) :: fn_oth_min
     REAL(wp) :: t_jmax_opt
     REAL(wp) :: t_jmax_omega
     REAL(wp) :: g0
     REAL(wp) :: g1
     REAL(wp) :: lambda_soil_ps
     REAL(wp) :: weibull_k_soil_ps
     REAL(wp) :: lambda_soil_gs
     REAL(wp) :: weibull_k_soil_gs
! Turnover 
     REAL(wp) :: tau_leaf                    !< leaf turnover time (years)
     REAL(wp) :: tau_branch                  !< branch turnover time (years)
     REAL(wp) :: tau_sap_wood                !< sap wood turnover time (years)
     REAL(wp) :: tau_fine_root               !< fine root turnover time (years)
     REAL(wp) :: tau_coarse_root             !< coarse root turnover time (years)
     REAL(wp) :: tau_fruit                   !< fruit turnover time (years)
     REAL(wp) :: tau_seed_litter             !< seed bed turnover time to litter (years)
     REAL(wp) :: tau_seed_est                !< seed bed turnover time to establishment (years)
     REAL(wp) :: wood_density 
     REAL(wp) :: k_latosa 
     REAL(wp) :: k_rtos 
     REAL(wp) :: k_crtos 
     REAL(wp) :: frac_sapwood_branch
     REAL(wp) :: phi_leaf_min 
     REAL(wp) :: k_root 
     REAL(wp) :: k_sapwood 
     REAL(wp) :: k2_fruit_alloc
     REAL(wp) :: allom_k1
     REAL(wp) :: allom_k2
     REAL(wp) :: c0_allom
     REAL(wp) :: vmax_uptake_n
     REAL(wp) :: soil_n_aff1 
     REAL(wp) :: soil_n_aff2
     REAL(wp) :: labile_cn_max
     REAL(wp) :: labile_cn_min
     REAL(wp) :: vmax_uptake_p
     REAL(wp) :: soil_p_aff1 
     REAL(wp) :: soil_p_aff2 
     REAL(wp) :: bnf_base
     REAL(wp) :: lambda_est_light
     REAL(wp) :: k_est_light
     REAL(wp) :: seed_size
! Phenology
     REAL(wp) :: beta_soil_flush             !< soil moisture limitation factor on stomatal conductance inducing leaf flushing
     REAL(wp) :: beta_soil_senescence        !< soil moisture limitation factor on stomatal conductance inducing leaf senescence
     REAL(wp) :: gdd_t_air_threshold         !< temperature threshold for the accumulation of growing degree days (K)
     REAL(wp) :: gdd_req_max                 !< maximum growing degree days requirements without chilling (degree days)
     REAL(wp) :: k_gdd_dormance              !< scaling factor in the GDD to number of dormant days relationship (days -1)
     REAL(wp) :: t_air_senescence            !< air temperature inducing senescence (K)
! Soil
     REAL(wp) :: k_root_dist                 !< scaling factor describing root distribution 

  END TYPE t_lctlib_element


  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_lctlib_class'

CONTAINS 

  !-----------------------------------------------------------------------------------------------------
  !> initialize parameter values across all pft
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION init_lctlib() RESULT(lctlib)

    USE mo_base_constants,          ONLY: npft, one_day, one_year
    USE mo_shared_constants,        ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2, Tzero 
    USE mo_atmland_constants
    USE mo_veg_constants
    USE mo_assimi_constants,        ONLY: ic3phot, ic4phot
    USE mo_pheno_constants,         ONLY: ievergreen, iraingreen, isummergreen, iperennial
    USE mo_assimi_config_class,     ONLY: assimi_config

    IMPLICIT NONE

      ! 0.1 InOut
      TYPE(t_lctlib_element),POINTER      :: lctlib(:)
      ! 0.2 Local
      INTEGER                             :: ipft         ! loop over pft
      CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_lctlib'

      ALLOCATE(lctlib(npft))    !    jsbach4: ALLOCATE(lctlib(nlct))

      !> plant functional types are: \n
      !! 1) broad leaved evergreen (rainforest) \n
      !! 2) broad leaved evergreen (xeric forest) \n
      !! 3) broad leaved deciduous (rain green) \n
      !! 4) broad leaved deciduous (summer green) \n
      !! 5) needle-leaved evergreen \n
      !! 6) needle-leaved deciduous \n
      !! 7) C3 grass \n
      !! 8) C4 grass

      ! vegetation growth form
      lctlib(1:npft)%growthform = & 
      (/itree,itree,itree,itree,itree,itree,igrass,igrass/)
      ! photosynthetic path way 
      lctlib(1:npft)%ps_pathway = & 
      (/ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic3phot,ic4phot/)
      ! photosynthetic path way 
      lctlib(1:npft)%phenology_type = & 
      (/ievergreen,ievergreen,iraingreen,isummergreen,ievergreen,isummergreen,iperennial,iperennial/)
      ! specific leaf area of a plant (m2 / mol C) 
      ! From Kattge et al. 2011, Table 5, including conversion from mm/mgDW to m2/mol
      lctlib(1:npft)%sla = 1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * &
      (/9.46_wp,8.99_wp,15.4_wp,15.4_wp,5.0_wp,10.09_wp,20.12_wp,19.23_wp/) 
      ! single leaf scattering albedo in the visible range
      !  extrapolated from Otto, 2014, BG for trees and Spitter et al. 1986 for grasses
      lctlib(1:npft)%sigma_vis = & 
      (/0.17_wp,0.17_wp,0.17_wp,0.17_wp,0.15_wp,0.15_wp,0.2_wp,0.2_wp/)
      ! single leaf scattering albedo in the near infrared range
      !  extrapolated from Otto, 2014, BG for trees and Spitter et al. 1986 for grasses
      lctlib(1:npft)%sigma_nir = & 
      (/0.76_wp,0.76_wp,0.76_wp,0.76_wp,0.73_wp,0.73_wp,0.8_wp,0.8_wp/)
      ! canopy clumping factor
      !  based on Campell and Norman 1999
      lctlib(1:npft)%omega_clumping = & 
      (/0.9_wp,0.9_wp,0.9_wp,0.9_wp,0.5_wp,0.7_wp,1.0_wp,1.0_wp/)
      IF(.NOT.assimi_config%flag_heterogeneous_canopy) lctlib(1:npft)%omega_clumping = 1.0_wp
    
      ! crown shape factor = 3.8 - 0.46 * crown depth / crown diameter, limited to range 1-3.34 
      !  based on Campell and Norman 1998, eq. 15.35
      lctlib(1:npft)%crown_shape_factor = & 
      (/3.34_wp,3.34_wp,3.34_wp,3.34_wp,2.19_wp,2.88_wp,3.34_wp,3.34_wp/)
      ! default leaf C:N (molC/molN) if not calculated dynamically 
      ! From Kattge et al. 2011, Table 5, including conversion from mg/g to mol/mol
      lctlib(1:npft)%cn_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_N / molar_mass_C / & 
      (/16.89_wp,13.73_wp,21.32_wp,21.32_wp,12.09_wp,19.37_wp,17.84_wp,14.14_wp/) 
      ! default leaf N:P (molN/molP) if not calculated dynamically 
      ! From Kattge et al. 2011, Table 5, including conversion from mg/g to mol/mol
      lctlib(1:npft)%np_leaf = carbon_per_dryweight_leaf * 1000._wp * molar_mass_P / molar_mass_C / lctlib(1:npft)%cn_leaf / &
      (/0.86_wp,0.84_wp,1.44_wp,1.44_wp,1.23_wp,1.83_wp,1.43_wp,1.36_wp/)
      ! fraction of leaf N not used for photosynthesis (Friend 1997) 
      ! temporarily tuned to match Amax/N in Kattge etal. 2012, Table 5
      lctlib(1:npft)%k0_fn_struc = & 
      !(/0.67_wp,0.67_wp,0.67_wp,0.67_wp,0.83_wp,0.83_wp,0.67_wp,0.67_wp/) 
      (/0.75_wp,0.70_wp,0.63_wp,0.63_wp,0.83_wp,0.79_wp,0.3_wp,0.3_wp/) 
      ! minimum fraction of non-photosynthetic leaf N. Calculated using range of leaf N values from GLOPNET 
      lctlib(1:npft)%fn_oth_min= lctlib(1:npft)%k0_fn_struc - 5._wp/molar_mass_n * 0.5_wp * 1.e3_wp * 71.4_wp * &
                                  molar_mass_N / 1000000._wp
      !Le Maire et al. 2012 (derived from grasses, does not give plausble results for needles)
      !0.012_wp / 0.48_wp * molar_mass_C / molar_mass_N / lctlib(1:npft)%sla/&
      !(1._wp / lctlib(1:npft)%sla / lctlib(1:npft)%cn_leaf) 
      ! initial value of the temperature optimum of electron transport and its shape parameter (deg C)
      lctlib(1:npft)%t_jmax_opt = & 
      (/30._wp,30._wp,25._wp,25._wp,25._wp,25._wp,25._wp,35._wp/)
      lctlib(1:npft)%t_jmax_omega = 18._wp
      ! intercept and slope of the An ~gs relationship (mmol/µmol)
      lctlib(1:npft)%g0 = & ! 0.01_wp
      (/0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.01_wp,0.03_wp,0.03_wp/)
      lctlib(1:npft)%g1 = & !9.3_wp
      ! g1 in Ball & Berry relationship, tuned to get similar Amax/N values as with Lin et al. g1 
      (/9.3_wp,8.3_wp,7.0_wp,10.9_wp,5.5_wp,7.0_wp,9.3_wp,2.0_wp/)
      ! g1 in kPa^0.5 as in Lin et al. 2015, Nat CC. OBS: no values for Larches - assigned broad leaved value!
      !(/3.77_wp,3.37_wp,2.98_wp,4.64_wp,2.35_wp,4.64_wp,4.16_wp,1.62_wp/)
      ! Weibull parameters describing the response of leaf-level photosynthesis to 
      ! fractional water availability
      lctlib(1:npft)%lambda_soil_ps = 3.0_wp
      lctlib(1:npft)%weibull_k_soil_ps = 3.0_wp
      ! Weibull parameters describing the response of stomatal conductance to 
      ! fractional water availability
      lctlib(1:npft)%lambda_soil_gs = lctlib(1:npft)%lambda_soil_ps
      lctlib(1:npft)%weibull_k_soil_gs = lctlib(1:npft)%weibull_k_soil_ps

! turnover times
      ! average turnover time of a leaf (years) 
      ! From Kattge et al. 2011, Table 5, including conversion from month lifetime to annual turnover
      lctlib(1:npft)%tau_leaf = & 
      (/16.83_wp,15.88_wp,5.83_wp,5.83_wp,39.7_wp,6.08_wp,3.85_wp,3.85_wp/) / 12.0_wp
      ! average turnover time of a fine root (years), Ahrens et al. 2014, NP 
      lctlib(1:npft)%tau_fine_root = & 
      (/0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp/) 
      ! average turnover time of a coarse root (years), Ahrens et al. 2014, NP 
      lctlib(1:npft)%tau_coarse_root = 8.0_wp 
      ! turnover time of the fraction of sapwood that is in branches (years)
      lctlib(1:npft)%tau_branch = 10._wp
      ! turnover time of the sapwood pool (years)
      lctlib(1:npft)%tau_sap_wood = 40.0_wp 
      ! turnover time of the fruit pool (years)
      lctlib(1:npft)%tau_fruit = 1.0_wp
      ! turnover time of the seed bed to litter (years)
      lctlib(1:npft)%tau_seed_litter = 2.0_wp
      ! turnover time of the seed bed to establishment (years)
      lctlib(1:npft)%tau_seed_est = 10.0_wp

! N uptake parameters
      ! vmax of N uptake (~µmol N / mol C / s)
      lctlib(1:npft)%vmax_uptake_n = 0.14_wp !714.28_wp   ! tuned from OCN's 1.5 gN/gC/1800s
      ! low-affinity uptake (L/mmol) 
      lctlib(1:npft)%soil_n_aff1 = 0.105_wp      ! tuned from OCN's 0.05 m2/gN assuming 150 L per m2
      ! half saturation of high-affinity uptake (mmol/L) 
      lctlib(1:npft)%soil_n_aff2 = 0.4_wp        ! tuned from OCN's 0.84 gN/m2 assuming 150 L per m2
      ! Maximum CN of labile pool, needed for N demand
      lctlib(1:npft)%labile_cn_max = 75._wp
      ! Minimum CN of labile pool, needed for N demand
      lctlib(1:npft)%labile_cn_min = 16._wp
      ! vmax of P uptake (~µmol P / mol C / s)
      lctlib(1:npft)%vmax_uptake_p = 0.1428_wp   ! tuned from OCN's 1.5 gN/gC/1800s
      ! low-affinity uptake (L/mmol) 
      lctlib(1:npft)%soil_p_aff1 = 0.105_wp      ! tuned from OCN's 0.05 m2/gN assuming 150 L per m2
      ! half saturation of high-affinity uptake (mmol/L) 
      lctlib(1:npft)%soil_p_aff2 = 0.4_wp        ! tuned from OCN's 0.84 gN/m2 assuming 150 L per m2
      ! placeholder constant BNf rate, required to keep soil NH4 at constant levels mumol N s-1 
      lctlib(1:npft)%bnf_base = 2.0_wp / molar_mass_N / one_day / one_year * 1.e6
      !(/2.6_wp, 2.6_wp, 2.6_wp, 5.0_wp, 2.6_wp, 2.6_wp, 2.6_wp, 2.6_wp/) * 1.e-2_wp *0.25_wp

! Vegetation dynamics parameters
      ! parameter in the Weibull function controlling light-limited establishment 
      lctlib(1:npft)%lambda_est_light = 1.1_wp        
      ! parameter in the Weibull function controlling light-limited establishment
      lctlib(1:npft)%k_est_light = 8.0_wp
      ! seed size (mol C)
      lctlib(1:npft)%seed_size = 0.5_wp 

! Phenology parameters
      ! soil moisture threshold for leaf flushing
      lctlib(1:npft)%beta_soil_flush = &
      (/0.0_wp,0.0_wp,0.0_wp,0.98_wp,0.0_wp,0.0_wp,0.9_wp,0.9_wp/)
      ! soil moisture threshold for leaf shedding 
      lctlib(1:npft)%beta_soil_senescence = &
      (/0.0_wp,0.0_wp,0.0_wp,0.5_wp,0.0_wp,0.0_wp,0.01_wp,0.01_wp/)
      ! GDD temperature threshold (above which GDDs are accumulated)
      lctlib(1:npft)%gdd_t_air_threshold = 5.0_wp + Tzero
      ! maximum GDD requirement (degC days) in the absence of chilling  
      lctlib(1:npft)%gdd_req_max = &
      (/0._wp,0._wp,0._wp,800._wp,0._wp,800._wp,10._wp,10._wp/)
      ! scaling coefficient in the decline of GDD requirements with increased length of the dormance period 
      lctlib(1:npft)%k_gdd_dormance = &
      (/0._wp,0._wp,0._wp,0.0098_wp,0._wp,0.0098_wp,0.1_wp,0.1_wp/)
      ! Weekly air temperature threshold inducing senecence of leaves 
      lctlib(1:npft)%t_air_senescence = & 
      (/0._wp,0._wp,0._wp,5.0_wp,0._wp,5.0_wp,0._wp,0._wp/) + Tzero

! Allocation paramters
      ! fraction of sapwood that is in branches 
      lctlib(1:npft)%frac_sapwood_branch = & 
      (/0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp,0.1_wp,1.0_wp,1.0_wp/) 
      ! wood density (mol C / m3)
      lctlib(1:npft)%wood_density = & 
      (/0.3_wp,0.3_wp,0.3_wp,0.3_wp,0.3_wp,0.3_wp,0.22_wp,0.22_wp/) * 1000000._wp/ molar_mass_C
      ! leaf area to sapwood area ratio 
      lctlib(1:npft)%k_latosa = &
       (/4000.0_wp, 4000.0_wp, 4000.0_wp, 4000.0_wp, 4000.0_wp, 4000.0_wp, 0.0_wp, 0.0_wp/)  
      ! coarse root to sapwood mass ratio 
      lctlib(1:npft)%k_crtos = 0.2_wp
      ! trade-off parametr for hydraulic investment into sapwood or fine roots 
      WHERE(lctlib(1:npft)%growthform /= igrass)
          lctlib(1:npft)%k_rtos = 1.0_wp / (lctlib(1:npft)%k_latosa / lctlib(1:npft)%sla / lctlib(1:npft)%wood_density)
        ELSEWHERE
          lctlib(1:npft)%k_rtos = 0.5_wp / sm2lm_grass
      END WHERE
      ! maximum fraction of biomass growth going to fruits 
      WHERE(lctlib(1:npft)%growthform /= igrass)
          lctlib(1:npft)%k2_fruit_alloc = 0.1_wp 
        ELSEWHERE
          lctlib(1:npft)%k2_fruit_alloc = 0.2_wp 
      END WHERE
      ! paramter in height diameter relationship 
      lctlib(1:npft)%allom_k1 = & 
      (/55.0_wp,55.0_wp,55.0_wp,55.0_wp,55.0_wp,55.0_wp,20._wp,20._wp/) 
      ! paramter in height diameter relationship
      lctlib(1:npft)%allom_k2 = & 
      (/0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.5_wp,0.5_wp/) 
      ! minimum leaf water potential, MPa 
      lctlib(1:npft)%phi_leaf_min = &  
      (/-1.0_wp,-1.0_wp,-1.0_wp,-1.0_wp,-1.0_wp,-1.0_wp,-0.75_wp,-0.75_wp/) 
      ! fine root hydraulic conductance, m3 mol-1 s-1 MPa-1 
      lctlib(1:npft)%k_root = 2.3_wp * molar_mass_C * 1.e-10 
      ! sapwood hydraulic conductance, m2 s-1 MPa-1 
      ! @todo check tropics check grasses 
      lctlib(1:npft)%k_sapwood = 1.e-3 * & 
      (/3.0_wp,3.0_wp,3.0_wp,3.0_wp,1.3_wp,1.3_wp,0._wp,0._wp/) 
      ! hydraulic trade-off parameter for fine root growth
      !! c0_allom: differentiate between trees & grasses
      DO ipft=1,npft
        IF(lctlib(ipft)%growthform == itree) THEN
          lctlib(ipft)%c0_allom = SQRT(lctlib(ipft)%k_root/lctlib(ipft)%k_sapwood * & 
                                lctlib(ipft)%tau_sap_wood/lctlib(ipft)%tau_fine_root * lctlib(ipft)%wood_density) * & 
                                lctlib(ipft)%tau_fine_root/lctlib(ipft)%tau_sap_wood
        ELSE
          lctlib(ipft)%c0_allom = 0.0_wp
        END IF
      END DO
      ! tuned according to Jackson et al. 1996, Oecologia, to reproduce their cummulative root distribution profiles
      ! only making use of the contrast sclerophyllic, tree, grass
      lctlib(1:npft)%k_root_dist = &
      (/3.5_wp,2.8_wp,2.5_wp,3.5_wp,3.5_wp,3.5_wp,5.5_wp,5.5_wp/)


    END FUNCTION init_lctlib 


END MODULE mo_jsb_lctlib_class
