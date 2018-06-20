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
!! 4th version by Jan Engel (2017-09-25); apply JSBACH4 memory usage via dsl4jsb macros, see r759 & #277
!!
MODULE mo_veg_util

#include "dsl_jsbach.inc"

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: reset_veg_fluxes    ! called in the land model interface prior to all Tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> init/reset vegetation fluxes (with/to zero)
  !! 
  !! called in the land model interface prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_veg_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_constants,        ONLY: zero
    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local 
    TYPE(t_jsb_model),        POINTER       :: model
    INTEGER                                 :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_veg_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    ! Declare pointers to variables in memory
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: maint_respiration_pot
    dsl4jsb_Real2D_onChunk      :: maint_respiration
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c13
    dsl4jsb_Real2D_onChunk      :: maint_respiration_c14
    dsl4jsb_Real2D_onChunk      :: growth_respiration
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c13
    dsl4jsb_Real2D_onChunk      :: growth_respiration_c14
    dsl4jsb_Real2D_onChunk      :: n_transform_respiration
    dsl4jsb_Real2D_onChunk      :: n_fixation_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c13
    dsl4jsb_Real2D_onChunk      :: n_processing_respiration_c14
    dsl4jsb_Real2D_onChunk      :: npp
    dsl4jsb_Real2D_onChunk      :: npp_c13
    dsl4jsb_Real2D_onChunk      :: npp_c14
    dsl4jsb_Real2D_onChunk      :: net_growth 
    dsl4jsb_Real2D_onChunk      :: uptake_nh4
    dsl4jsb_Real2D_onChunk      :: uptake_nh4_n15
    dsl4jsb_Real2D_onChunk      :: uptake_no3
    dsl4jsb_Real2D_onChunk      :: uptake_no3_n15
    dsl4jsb_Real2D_onChunk      :: n_fixation
    dsl4jsb_Real2D_onChunk      :: n_fixation_n15
    dsl4jsb_Real2D_onChunk      :: uptake_po4
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_n15
    dsl4jsb_Real2D_onChunk      :: recycling_leaf_p
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_n15
    dsl4jsb_Real2D_onChunk      :: recycling_fine_root_p
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_n15
    dsl4jsb_Real2D_onChunk      :: recycling_heart_wood_p
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: unit_transpiration
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(VEG_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)           ! veg__mem => veg_memory
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    !lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_transform_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_processing_respiration_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c13)
    dsl4jsb_Get_var2D_onChunk(VEG_, npp_c14)
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_nh4_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_no3_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation)
    dsl4jsb_Get_var2D_onChunk(VEG_, n_fixation_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, uptake_po4)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_leaf_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_fine_root_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_, recycling_heart_wood_p)
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)
    dsl4jsb_Get_var2D_onChunk(VEG_, unit_transpiration)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! t_slm_struct_veg_pool
    CALL veg__mem%growth%init(zero)
    CALL veg__mem%litterfall%init(zero)

    ! t_slm_struct_veg_element
    CALL veg__mem%reserve_use%init(zero)
    CALL veg__mem%establishment%init(zero)
    CALL veg__mem%exudation%init(zero)

    ! flux variables
    maint_respiration_pot(:) = zero          !< maintenance respiration (µmol CO2 m-2 s-1)
    maint_respiration(:) = zero              !< maintenance respiration (µmol CO2 m-2 s-1)
    maint_respiration_c13(:) = zero          !< maintenance respiration (µmol 13CO2 m-2 s-1)
    maint_respiration_c14(:) = zero          !< maintenance respiration (µmol 14CO2 m-2 s-1)
    growth_respiration(:) = zero             !< growth respiration (µmol CO2 m-2 s-1)
    growth_respiration_c13(:) = zero         !< growth respiration (µmol 13CO2 m-2 s-1)
    growth_respiration_c14(:) = zero         !< growth respiration (µmol 14CO2 m-2 s-1)
    n_transform_respiration(:) = zero        !< respiration associated with N transformation (µmol CO2 m-2 s-1)
    n_fixation_respiration(:) = zero         !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    n_processing_respiration(:) = zero       !< respiration associated with N fixation (µmol CO2 m-2 s-1)
    n_processing_respiration_c13(:) = zero   !< respiration associated with N fixation (µmol 13CO2 m-2 s-1)
    n_processing_respiration_c14(:) = zero   !< respiration associated with N fixation (µmol 14CO2 m-2 s-1)
    npp(:) = zero                            !< net primary production (µmol CO2 m-2 s-1)
    npp_c13(:) = zero                        !< net primary production (µmol 13CO2 m-2 s-1)
    npp_c14(:) = zero                        !< net primary production (µmol 14CO2 m-2 s-1)
    net_growth(:) = zero                     !< net growth (µmol CO2 m-2 s-1)
    uptake_nh4(:) = zero                     !< NH4 uptake (µmol N m-2 s-1) 
    uptake_nh4_n15(:) = zero                 !< 15NH4 uptake (µmol N m-2 s-1)
    uptake_no3(:) = zero                     !< NO3 uptake (µmol N m-2 s-1)
    uptake_no3_n15(:) = zero                 !< 15NO3 uptake (µmol N m-2 s-1)
    n_fixation(:) = zero                     !< ??
    n_fixation_n15(:) = zero                 !< ??
    uptake_po4(:) = zero                     !< ??
    recycling_leaf_n(:) = zero               !< net flux of N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    recycling_leaf_n15(:) = zero             !< net flux of 15N from leaf to labile (senescence and maintenance) (mol m-2 timestep-1) 
    recycling_leaf_p(:) = zero               !< net flux of P from leaf to labile (senescence and maintenance) (mol m-2 timestep-1)
    recycling_fine_root_n(:) = zero          !< net flux of N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    recycling_fine_root_n15(:) = zero        !< net flux of 15N from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1) 
    recycling_fine_root_p(:) = zero          !< net flux of P from fine roots to labile (senescence and maintenance) (mol m-2 timestep-1)
    recycling_heart_wood_n(:) = zero         !< net flux of N from heart wood to labile (senescence and maintenance) (mol m-2 timestep-1)
    recycling_heart_wood_n15(:) = zero       !< net flux of 15N from heart wood to labile (senescence and maintenance) (mol m-2 timestep-1) 
    recycling_heart_wood_p(:) = zero         !< net flux of P from heart wood to labile (senescence and maintenance) (mol m-2 timestep-1)

    ! other variables
    delta_dens_ind(:) = zero                 !< ??
    unit_transpiration(:) = zero             !< ??

  END SUBROUTINE reset_veg_fluxes

END MODULE mo_veg_util
