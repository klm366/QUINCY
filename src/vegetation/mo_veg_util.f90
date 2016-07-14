!>
!! @brief diverse subroutines and functions for the plant-model
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-17); was deleted because empty at Rev. 126
!! 2nd version by Jan Engel (2016-05-11)
!!
MODULE mo_veg_util

  USE mo_kind, ONLY: wp

  IMPLICIT NONE

  PUBLIC :: set_zero_fluxes_of_plant_pools

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> initalize and set zero elements of all tissues 
  !!
  !! is called for flux%litterfall & flux%growth each timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_zero_flux_of_plant_pools(one_flux_of_pools)

    USE mo_veg_memory,  ONLY: t_slm_veg_pool
    USE mo_shared_memory, ONLY: element_zero

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_pool), INTENT(inout) :: one_flux_of_pools

     one_flux_of_pools%leaf         = element_zero   ! TYPE(t_element) :: element_zero; see mo_shared_memory
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
  SUBROUTINE set_zero_fluxes_of_plant_pools(tile)

    USE mo_base_memory,       ONLY: t_slm_tile
    USE mo_veg_memory,        ONLY: t_slm_veg_flux
    USE mo_shared_memory,       ONLY: element_zero

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET :: tile    !< one tile
    ! 0.2 Local
    TYPE(t_slm_veg_flux),POINTER           :: flux    ! fluxes of plant of the certain tile

    !  bi-directional pointer
    flux => tile%plant%flux

    CALL set_zero_flux_of_plant_pools(flux%growth)
    CALL set_zero_flux_of_plant_pools(flux%litterfall)

    flux%reserve_use             = element_zero    ! TYPE(t_element) :: reserve_use

    flux%gpp                     = 0.0_wp  !< gross photosynthesis (µmol CO2 m-2 s-1)
    flux%an                      = 0.0_wp  !< net photosynthesis (µmol CO2 m-2 s-1)
    flux%gpp_c13                 = 0.0_wp  !< gross photosynthesis (µmol 13CO2 m-2 s-1)
    flux%gpp_c14                 = 0.0_wp  !< gross photosynthesis (µmol 14CO2 m-2 s-1)
    flux%maint_respiration_leaf  = 0.0_wp  !< foliar maintenance respiration (µmol CO2 m-2 s-1)
    flux%maint_respiration_pot   = 0.0_wp  !< maintenance respiration (µmol CO2 m-2 s-1)
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
    flux%recycling_leaf_n        = 0.0_wp  !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    flux%recycling_leaf_n15      = 0.0_wp  !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
    flux%recycling_leaf_p        = 0.0_wp  !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    flux%recycling_fine_root_n   = 0.0_wp  !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    flux%recycling_fine_root_n15 = 0.0_wp  !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
    flux%recycling_fine_root_p   = 0.0_wp  !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    flux%recycling_sap_wood_n    = 0.0_wp  !< net flux of N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
    flux%recycling_sap_wood_n15  = 0.0_wp  !< net flux of 15N from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1) 
    flux%recycling_sap_wood_p    = 0.0_wp  !< net flux of P from sapwood to labile (senescence and maintenance) (mol m-2 timestep-1)
    flux%transpiration           = 0.0_wp  !< transpiration in kg m-2 s-1
    flux%ci                      = 0.0_wp  !< ??
  END SUBROUTINE set_zero_fluxes_of_plant_pools



END MODULE mo_veg_util
