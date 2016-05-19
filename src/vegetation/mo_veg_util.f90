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

!  PUBLIC :: 

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialize all plant variables for pools and so, not fluxes
  !! 
  !! called in the driver; for fluxes see mo_veg_memory:set_zero_fluxes_of_plant_pools
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_plant(tile)

    USE mo_base_memory,       ONLY: t_slm_tile
    USE mo_veg_memory,        ONLY: t_slm_plant, element_zero

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout), TARGET :: tile    !< one tile
    ! 0.2 Local
    TYPE(t_slm_plant),POINTER              :: plant   ! plant of one tile

    !  bi-directional pointer
    plant => tile%plant

    ! pools
    plant%pool%leaf         = element_zero   ! TYPE(t_element) :: element_zero
    plant%pool%fine_root    = element_zero
    plant%pool%coarse_root  = element_zero
    plant%pool%sap_wood     = element_zero
    plant%pool%heart_wood   = element_zero
    plant%pool%labile       = element_zero
    plant%pool%reserve      = element_zero
    plant%pool%fruit        = element_zero
    plant%pool%seed_bed     = element_zero
    ! canopy layers
    plant%canopy_layer(:)%nl_tot       = 0.0_wp
    plant%canopy_layer(:)%fnl_chl      = 0.0_wp
    plant%canopy_layer(:)%fnl_et       = 0.0_wp
    plant%canopy_layer(:)%fnl_rub      = 0.0_wp
    plant%canopy_layer(:)%par_sunlit   = 0.0_wp
    plant%canopy_layer(:)%par_shaded   = 0.0_wp
    plant%canopy_layer(:)%fleaf_sunlit = 0.0_wp
    ! plant variables
    plant%pftid            = 1          ! INTEGER
    plant%height           = 0.0_wp
    plant%diameter         = 0.0_wp
    plant%lai              = 0.0_wp
    plant%target_cn_leaf   = 0.0_wp
    plant%target_np_leaf   = 0.0_wp
    plant%mean_leaf_age    = 0.0_wp
    plant%gc               = 0.0_wp
    plant%f_n_demand       = 0.0_wp
    plant%f_p_demand       = 0.0_wp
    plant%beta_air         = 0.0_wp
    plant%beta_soil_ps     = 0.0_wp
    plant%beta_soil_gs     = 0.0_wp
    plant%growing_season   = .TRUE.     ! LOGICAL
  END SUBROUTINE init_plant

  !-----------------------------------------------------------------------------------------------------
  !> initialize the derived TYPE t_element element_zero with zero
  !! 
  !! this is done only one time, this element_zero is used to set other derived TYPE's of t_element to zero
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_element_zero

    USE mo_veg_memory,  ONLY: element_zero   ! TYPE(t_element) :: element_zero | declared in mo_veg_memory

    element_zero%carbon     = 0.0_wp
    element_zero%nitrogen   = 0.0_wp
    element_zero%phosphorus = 0.0_wp
    element_zero%carbon13   = 0.0_wp
    element_zero%carbon14   = 0.0_wp
    element_zero%nitrogen15 = 0.0_wp
  END SUBROUTINE init_element_zero

  !-----------------------------------------------------------------------------------------------------
  !> initalize and set zero elements of all tissues 
  !!
  !! is called for flux%litterfall & flux%growth each timestep
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_zero_flux_of_plant_pools(one_flux_of_pools)

    USE mo_veg_memory,  ONLY: t_slm_veg_pool, element_zero

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_veg_pool), INTENT(inout) :: one_flux_of_pools

     one_flux_of_pools%leaf         = element_zero   ! TYPE(t_element) :: element_zero; see mo_land_memory and mo_veg_memory
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
    USE mo_veg_memory,        ONLY: t_slm_veg_flux, element_zero

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
    flux%leaf_resorption_n       = 0.0_wp  !< resorption of N upon leaf senescence (mol m-2 timestep-1)
    flux%leaf_resorption_n15     = 0.0_wp  !< resorption of 15N upon leaf senescence (mol m-2 timestep-1)
    flux%leaf_resorption_p       = 0.0_wp  !< resorption of P upon leaf senescence (mol m-2 timestep-1)
    flux%transpiration           = 0.0_wp  !< transpiration in kg m-2 s-1
    flux%ci                      = 0.0_wp  !< ??
  END SUBROUTINE set_zero_fluxes_of_plant_pools



END MODULE mo_veg_util
