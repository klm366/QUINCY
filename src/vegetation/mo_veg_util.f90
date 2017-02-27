!>
!! @brief diverse subroutines and functions for the plant-model
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (2016-03-17); was deleted because empty at Rev. 126
!! 2nd version by Jan Engel (2016-05-11)
!! 3rd version by Jan Engel (2016-07-15); all content removed, see SVN r198 & trac ticket #63
!!
MODULE mo_veg_util

#include "dsl_jsbach.inc"

  USE mo_kind, ONLY: wp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: set_zero_vegetation_fluxes    ! called in the land model interface prior to all Tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_util'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> set all vegetation fluxes to zero
  !! 
  !! called in the land model interface prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE set_zero_vegetation_fluxes(tile)

    USE mo_base_memory,           ONLY: t_slm_tile
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_slm_tile),INTENT(inout)      :: tile           !< one tile
    ! 0.2 Local
    REAL(wp) :: flux_zero
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_zero_vegetation_fluxes'
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory), POINTER      :: veg__mem  
    TYPE(t_veg_memory), POINTER         :: veg        ! local pointer to veg__mem for SLM usability
    dsl4jsb_Get_memory(VEG_)                          ! veg__mem => tile%mem(VEG_)%p

    veg => veg__mem 

    flux_zero = 0.0_wp


    ! t_slm_struct_veg_element
    CALL veg%reserve_use%init(flux_zero)
    CALL veg%establishment%init(flux_zero)
    ! t_slm_struct_veg_pool
    CALL veg%growth%init(flux_zero)
    CALL veg%litterfall%init(flux_zero)
    ! variables
    CALL veg%maint_respiration_pot%init(flux_zero)         !< maintenance respiration (µmol CO2 m-2 s-1)
    CALL veg%maint_respiration%init(flux_zero)             !< maintenance respiration (µmol CO2 m-2 s-1)
    CALL veg%maint_respiration_c13%init(flux_zero)         !< maintenance respiration (µmol 13CO2 m-2 s-1)
    CALL veg%maint_respiration_c14%init(flux_zero)         !< maintenance respiration (µmol 14CO2 m-2 s-1)
    CALL veg%growth_respiration%init(flux_zero)            !< growth respiration (µmol CO2 m-2 s-1)
    CALL veg%growth_respiration_c13%init(flux_zero)        !< growth respiration (µmol 13CO2 m-2 s-1)
    CALL veg%growth_respiration_c14%init(flux_zero)        !< growth respiration (µmol 14CO2 m-2 s-1)
    CALL veg%n_transform_respiration%init(flux_zero)       !< respiration associated with N transformation (µmol CO2 m-2 s-1)
    CALL veg%n_fixation_respiration%init(flux_zero)        !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    CALL veg%n_processing_respiration%init(flux_zero)      !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    CALL veg%n_processing_respiration_c13%init(flux_zero)  !< respiration associated with N fixation (µmol 13CO2 m-2 s-1)
    CALL veg%n_processing_respiration_c14%init(flux_zero)  !< respiration associated with N fixation (µmol 14CO2 m-2 s-1)
    CALL veg%npp%init(flux_zero)                           !< net primary production (µmol CO2 m-2 s-1)
    CALL veg%npp_c13%init(flux_zero)                       !< net primary production (µmol 13CO2 m-2 s-1)
    CALL veg%npp_c14%init(flux_zero)                       !< net primary production (µmol 14CO2 m-2 s-1)
    CALL veg%uptake_nh4%init(flux_zero)                    !< NH4 uptake (µmol N m-2 s-1) 
    CALL veg%uptake_nh4_n15%init(flux_zero)                !< 15NH4 uptake (µmol N m-2 s-1)
    CALL veg%uptake_no3%init(flux_zero)                    !< NO3 uptake (µmol N m-2 s-1)
    CALL veg%uptake_no3_n15%init(flux_zero)                !< 15NO3 uptake (µmol N m-2 s-1)
    CALL veg%n_fixation%init(flux_zero)                    !< ??
    CALL veg%n_fixation_n15%init(flux_zero)                !< ??
    CALL veg%uptake_p%init(flux_zero)                      !< ??
    CALL veg%recycling_leaf_n%init(flux_zero)              !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_leaf_n15%init(flux_zero)            !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
    CALL veg%recycling_leaf_p%init(flux_zero)              !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_fine_root_n%init(flux_zero)         !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_fine_root_n15%init(flux_zero)       !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
    CALL veg%recycling_fine_root_p%init(flux_zero)         !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_heart_wood_n%init(flux_zero)        !< net flux of N from heart wood to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%recycling_heart_wood_n15%init(flux_zero)      !< net flux of 15N from heart wood to labile (senescence and maintenance) (mol m-2 timestep-1) 
    CALL veg%recycling_heart_wood_p%init(flux_zero)        !< net flux of P from heart wood to labile (senescence and maintenance) (mol m-2 timestep-1)
    CALL veg%transpiration%init(flux_zero)                 !< transpiration in kg m-2 s-1
    CALL veg%delta_dens_ind%init(flux_zero)                !< ??
    CALL veg%unit_transpiration%init(flux_zero)            !< ??

  END SUBROUTINE set_zero_vegetation_fluxes

END MODULE mo_veg_util
