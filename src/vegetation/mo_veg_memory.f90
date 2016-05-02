!>
!! @brief Contains the routines for plant memory operations 
!!
!! @author
!!  Sönke Zaehle, MPI-BGC Jena 
!!
!! @par Revision History
!! First version by Sönke Zaehle (2015-12-04)
!!
MODULE mo_veg_memory

  IMPLICIT NONE

  PUBLIC :: VegPoolAdd,VegPoolSub,VegPoolMul,VegPoolDiv,set_zero_fluxes_of_plant_pools

CONTAINS

 FUNCTION VegPoolAdd(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(t_slm_veg_pool), INTENT(in) :: p1,p2
     TYPE(t_slm_veg_pool) :: p3
     p3%leaf          = ElementAdd(p1%leaf,p2%leaf)
     p3%fine_root     = ElementAdd(p1%fine_root,p2%fine_root)
     p3%coarse_root   = ElementAdd(p1%coarse_root,p2%coarse_root)
     p3%sap_wood      = ElementAdd(p1%sap_wood,p2%sap_wood)
     p3%heart_wood    = ElementAdd(p1%heart_wood,p2%heart_wood)
     p3%labile        = ElementAdd(p1%labile,p2%labile)
     p3%reserve       = ElementAdd(p1%reserve,p2%reserve)
     p3%fruit         = ElementAdd(p1%fruit,p2%fruit)
     p3%seed_bed      = ElementAdd(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolAdd  

 FUNCTION VegPoolSub(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(t_slm_veg_pool), INTENT(in) :: p1,p2
     TYPE(t_slm_veg_pool) :: p3
     p3%leaf = ElementSub(p1%leaf,p2%leaf)
     p3%fine_root = ElementSub(p1%fine_root,p2%fine_root)
     p3%coarse_root = ElementSub(p1%coarse_root,p2%coarse_root)
     p3%sap_wood = ElementSub(p1%sap_wood,p2%sap_wood)
     p3%heart_wood = ElementSub(p1%heart_wood,p2%heart_wood)
     p3%labile = ElementSub(p1%labile,p2%labile)
     p3%reserve = ElementSub(p1%reserve,p2%reserve)
     p3%fruit = ElementSub(p1%fruit,p2%fruit)
     p3%seed_bed = ElementSub(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolSub  

 FUNCTION VegPoolMul(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(t_slm_veg_pool), INTENT(in) :: p1,p2
     TYPE(t_slm_veg_pool) :: p3
     p3%leaf = ElementMul(p1%leaf,p2%leaf)
     p3%fine_root = ElementMul(p1%fine_root,p2%fine_root)
     p3%coarse_root = ElementMul(p1%coarse_root,p2%coarse_root)
     p3%sap_wood = ElementMul(p1%sap_wood,p2%sap_wood)
     p3%heart_wood = ElementMul(p1%heart_wood,p2%heart_wood)
     p3%labile = ElementMul(p1%labile,p2%labile)
     p3%reserve = ElementMul(p1%reserve,p2%reserve)
     p3%fruit = ElementMul(p1%fruit,p2%fruit)
     p3%seed_bed = ElementMul(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolMul  

 FUNCTION VegPoolDiv(p1,p2) RESULT(p3)
     USE mo_base_memory
     IMPLICIT NONE
     TYPE(t_slm_veg_pool), INTENT(in) :: p1,p2
     TYPE(t_slm_veg_pool) :: p3
     p3%leaf = ElementDiv(p1%leaf,p2%leaf)
     p3%fine_root = ElementDiv(p1%fine_root,p2%fine_root)
     p3%coarse_root = ElementDiv(p1%coarse_root,p2%coarse_root)
     p3%sap_wood = ElementDiv(p1%sap_wood,p2%sap_wood)
     p3%heart_wood = ElementDiv(p1%heart_wood,p2%heart_wood)
     p3%labile = ElementDiv(p1%labile,p2%labile)
     p3%reserve = ElementDiv(p1%reserve,p2%reserve)
     p3%fruit = ElementDiv(p1%fruit,p2%fruit)
     p3%seed_bed = ElementDiv(p1%seed_bed,p2%seed_bed)
 END FUNCTION VegPoolDiv  


  !-----------------------------------------------------------------------------------------------------
  !> initalize and set zero elements of all tissues 
  !!
  !! is called for flux%litterfall & flux%growth each timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_zero_flux_of_plant_pools(one_flux_of_pools)
    USE mo_base_memory
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_pool), INTENT(inout) :: one_flux_of_pools

     one_flux_of_pools%leaf         = element_zero   ! TYPE(t_element) :: element_zero; see mo_base_memory
     one_flux_of_pools%fine_root    = element_zero
     one_flux_of_pools%coarse_root  = element_zero
     one_flux_of_pools%sap_wood     = element_zero
     one_flux_of_pools%heart_wood   = element_zero
     one_flux_of_pools%labile       = element_zero
     one_flux_of_pools%reserve      = element_zero
     one_flux_of_pools%fruit        = element_zero
     one_flux_of_pools%seed_bed     = element_zero
  END SUBROUTINE set_zero_flux_of_plant_pools

  !-----------------------------------------------------------------------------------------------------
  !> initalize and set zero all vegetation memory 
  !!
  !! is initially called in plant_model_driver, 
  !! and then called every timestep in beginning of mo_veg_interface
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_zero_fluxes_of_plant_pools(flux)
    USE mo_base_memory
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_flux), INTENT(inout) :: flux

    CALL set_zero_flux_of_plant_pools(flux%growth)
    CALL set_zero_flux_of_plant_pools(flux%litterfall)

    flux%reserve_use             = element_zero    ! TYPE(t_element) :: reserve_use

    flux%gpp                     = 0.0_wp  !< gross photosynthesis (µmol CO2 m-2 s-1)
    flux%an                      = 0.0_wp  !< net photosynthesis (µmol CO2 m-2 s-1)
    flux%gpp_c13                 = 0.0_wp  !< gross photosynthesis (µmol 13CO2 m-2 s-1)
    flux%gpp_c14                 = 0.0_wp  !< gross photosynthesis (µmol 14CO2 m-2 s-1)
    flux%maint_respiration_leaf  = 0.0_wp  !< foliar maintenance respiration (µmol CO2 m-2 s-1)
    flux%maint_respiration       = 0.0_wp  !< maintenance respiration (µmol CO2 m-2 s-1)
    flux%maint_respiration_c13   = 0.0_wp  !< maintenance respiration (µmol 13CO2 m-2 s-1)
    flux%maint_respiration_c14   = 0.0_wp  !< maintenance respiration (µmol 14CO2 m-2 s-1)
    flux%growth_respiration      = 0.0_wp  !< growth respiration (µmol CO2 m-2 s-1)
    flux%growth_respiration_c13  = 0.0_wp  !< growth respiration (µmol 13CO2 m-2 s-1)
    flux%growth_respiration_c14  = 0.0_wp  !< growth respiration (µmol 14CO2 m-2 s-1)
    flux%npp                     = 0.0_wp  !< net primary production (µmol CO2 m-2 s-1)
    flux%npp_c13                 = 0.0_wp  !< net primary production (µmol 13CO2 m-2 s-1)
    flux%npp_c14                 = 0.0_wp  !< net primary production (µmol 14CO2 m-2 s-1)
    flux%uptake_n                = 0.0_wp  !< ??
    flux%uptake_n15              = 0.0_wp  !< ??
    flux%uptake_p                = 0.0_wp  !< ??
    flux%leaf_resorption_n       = 0.0_wp  !< resorption of N upon leaf senescence (mol m-2 timestep-1)
    flux%leaf_resorption_n15     = 0.0_wp  !< resorption of 15N upon leaf senescence (mol m-2 timestep-1)
    flux%leaf_resorption_p       = 0.0_wp  !< resorption of P upon leaf senescence (mol m-2 timestep-1)
    flux%transpiration           = 0.0_wp  !< transpiration in kg m-2 s-1
    flux%ci                      = 0.0_wp  !< ??
  END SUBROUTINE set_zero_fluxes_of_plant_pools

  !-----------------------------------------------------------------------------------------------------
  !> set intial values for parameters that accumulate values (average) across a certain period; 
  !!
  !! called in plant_model_driver section 1.0; 
  !! the interval is defined in mo_veg_parameters
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_plant_variable_accumulation(accumulated)
    USE mo_base_memory
    IMPLICIT NONE
    ! 0.2 Output
    TYPE(AccumulationPlantClass), INTENT(inout) :: accumulated 
    ! 
    ! 1.0 initialize default values; tile(:)%plant%accumulated
    ! 
    accumulated%par_sunlit        = 10._wp
    accumulated%par_shaded        = 10._wp
    accumulated%fleaf_sunlit      = 10._wp
    accumulated%beta_air          = 1._wp
    accumulated%beta_soil_ps      = 1._wp
    accumulated%beta_soil_gs      = 1._wp
    accumulated%gpp_week               = 1._wp
    accumulated%maint_respiration_week = 1._wp
  END SUBROUTINE init_plant_variable_accumulation

END MODULE mo_veg_memory
