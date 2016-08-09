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
     REAL(wp) :: tau_leaf
     REAL(wp) :: sla
     REAL(wp) :: cn_leaf
     REAL(wp) :: np_leaf
     REAL(wp) :: k0_fnl_struc
     REAL(wp) :: t_jmax_opt
     REAL(wp) :: t_jmax_omega
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


  PRIVATE
  PUBLIC :: init_lctlib, t_lctlib_element

CONTAINS 

  !-----------------------------------------------------------------------------------------------------
  !> initialize parameter values across all pft
  !!
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION init_lctlib() RESULT(lctlib)

    USE mo_base_constants,          ONLY: npft
    USE mo_shared_constants,        ONLY: molar_mass_C, molar_mass_N, molar_mass_P, Dwv, Dco2 
    USE mo_atmland_constants
    USE mo_veg_constants

    IMPLICIT NONE

      TYPE(t_lctlib_element),POINTER      :: lctlib(:)
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
      ! average turnover (1/lifetime) of a leaf (1/years) 
      ! From Kattge et al. 2011, Table 5, including conversion from month lifetime to annual turnover
      lctlib(1:npft)%tau_leaf = 12._wp / & 
      (/16.83_wp,15.88_wp,5.83_wp,5.83_wp,39.7_wp,6.08_wp,3.85_wp,3.85_wp/)
      ! specific leaf area of a plant (m2 / mol C) 
      ! From Kattge et al. 2011, Table 5, including conversion from mm/mgDW to m2/mol
      lctlib(1:npft)%sla = 1._wp / carbon_per_dryweight_leaf / 1000.0_wp * molar_mass_C * &
      (/9.46_wp,8.99_wp,15.4_wp,15.4_wp,5.0_wp,10.09_wp,20.12_wp,19.23_wp/) 
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
      lctlib(1:npft)%k0_fnl_struc = & 
      !(/0.67_wp,0.67_wp,0.67_wp,0.67_wp,0.83_wp,0.83_wp,0.67_wp,0.67_wp/) 
      (/0.75_wp,0.70_wp,0.63_wp,0.63_wp,0.83_wp,0.79_wp,0.3_wp,0.3_wp/) 
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
      ! turnover time of the fraction of sapwood that is in branches (10%) (1/years)
      lctlib(1:npft)%tau_branches = 0.01_wp
      ! turnover time of the sapwood pool (1/years)
      lctlib(1:npft)%tau_sap_wood = 0.05_wp
      ! wood density (mol C / m3)
      lctlib(1:npft)%wood_density = & 
      (/30000.0_wp,30000.0_wp,30000.0_wp,30000.0_wp,30000.0_wp,30000.0_wp,22000._wp,22000._wp/) / molar_mass_C
      ! leaf area to sapwood area ratio 
      lctlib(1:npft)%k_latosa = 4000.0_wp 
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
          lctlib(1:npft)%k_max_falloc_fruit = 0.1_wp 
        ELSEWHERE
          lctlib(1:npft)%k_max_falloc_fruit = 0.2_wp 
      END WHERE
      ! average turnover (1/lifetime) of a fine root (1/years) 
      lctlib(1:npft)%tau_fine_root = & 
      (/0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp,0.7_wp,1.2_wp,1.2_wp/) 
      ! average turnover (1/lifetime) of a coarse root (1/years) 
      lctlib(1:npft)%tau_coarse_root = 0.1_wp 
      ! paramter in height diameter relationship 
      lctlib(1:npft)%allom_k1 = & 
      (/55.0_wp,55.0_wp,55.0_wp,55.0_wp,55.0_wp,55.0_wp,20._wp,20._wp/) 
      ! paramter in height diameter relationship
      lctlib(1:npft)%allom_k2 = & 
      (/0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.65_wp,0.5_wp,0.5_wp/) 
      ! turnover time of the fruit pool (1/years)
      lctlib(1:npft)%tau_fruit = 1.0_wp
      ! turnover time of the seed bed to litter (1/years)
      lctlib(1:npft)%tau_seed_litter = 0.1_wp
      ! turnover time of the seed bed to establishment (1/years)
      lctlib(1:npft)%tau_seed_est = 0.3_wp
      ! vmax of N uptake (~µmol N / mol C / s)
      lctlib(1:npft)%vmax_uptake_n = 14._wp !714.28_wp   ! tuned from OCN's 1.5 gN/gC/1800s
      ! low-affinity uptake (L/mmol) 
      lctlib(1:npft)%soil_n_aff1 = 0.105_wp      ! tuned from OCN's 0.05 m2/gN assuming 150 L per m2
      ! half saturation of high-affinity uptake (mmol/L) 
      lctlib(1:npft)%soil_n_aff2 = 0.4_wp        ! tuned from OCN's 0.84 gN/m2 assuming 150 L per m2
      ! vmax of P uptake (~µmol P / mol C / s)
      lctlib(1:npft)%vmax_uptake_p = 14.28_wp   ! tuned from OCN's 1.5 gN/gC/1800s
      ! low-affinity uptake (L/mmol) 
      lctlib(1:npft)%soil_p_aff1 = 0.105_wp      ! tuned from OCN's 0.05 m2/gN assuming 150 L per m2
      ! half saturation of high-affinity uptake (mmol/L) 
      lctlib(1:npft)%soil_p_aff2 = 0.4_wp        ! tuned from OCN's 0.84 gN/m2 assuming 150 L per m2
      ! placeholder constant BNf rate, required to keep soil NH4 at constant levels mumol N s-1 
      lctlib(1:npft)%bnf_base = & 
      (/2.6_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9, 5.0_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9, 2.6_wp*1.e-9/) 
      ! parameter in the Weibull function controlling light-limited establishment 
      lctlib(1:npft)%lambda_est_light = 1.1_wp        
      ! parameter in the Weibull function controlling light-limited establishment
      lctlib(1:npft)%k_est_light = 4.0_wp
      ! seed size (mol C)
      lctlib(1:npft)%seed_size = 0.1_wp

    END FUNCTION init_lctlib 


END MODULE mo_jsb_lctlib_class
