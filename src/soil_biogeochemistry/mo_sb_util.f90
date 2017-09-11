!>
!! @brief Contains the utilities for the soil_biogeochemistry process
!!
!! @author
!!  Jan Engel, MPI-BGC Jena
!!
!! @par Revision History
!! First version by Jan Engel (20170104)
!!
MODULE mo_sb_util
#include "dsl_jsbach.inc"
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp
  USE mo_exception, ONLY: finish

  IMPLICIT NONE
  PRIVATE
  PUBLIC reset_sb_fluxes

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_sb_util'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  !> init/reset sb fluxes (with/to zero)
  !! 
  !! called in the land model interface prior to all Tasks
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_sb_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_constants,        ONLY: zero

    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET  :: tile         !< one tile
    TYPE(t_jsb_task_options), INTENT(in)             :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    !TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_sb_fluxes'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    ! 2D SB_
    dsl4jsb_Real2D_onChunk :: total_flux_carbon
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen15
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: emission_nh3
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2
    ! 3D SB_
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_sorption_po4
    dsl4jsb_Real3D_onChunk :: fast_adsorption_po4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: gross_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow 
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2_n15
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: emission_nh3_n15
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4
    dsl4jsb_Real3D_onChunk :: nitrification_no3
    dsl4jsb_Real3D_onChunk :: nitrification_noy
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_noy
    dsl4jsb_Real3D_onChunk :: denitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_n2
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation 
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15
! 2 be removed later: pointer to sb   *slmtemp
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process config
!    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    !lctlib => model%lctlib(tile%lcts%lib_id) 
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_nh3)                ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)                ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)                ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)                 ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)              ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_sorption_po4)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_adsorption_po4)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       gross_mineralisation_po4)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_n15_sl)! out
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)        ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o)               ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2)                ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)         ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)          ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)  ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)! out
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)! out
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4_n15)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3_n15)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy_n15)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o_n15)       ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy_n15)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o_n15)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2_n15)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)     ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy_n15)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o_n15)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2_n15)            ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute)     ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute)     ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_nh3_n15)            ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy_n15)            ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o_n15)            ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2_n15)             ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation_n15)        ! out


! 2 be removed  *slmtemp
    sb   => sb__mem


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 initialise/reset fluxes
    !!
    ! 2D fields  
    total_flux_carbon(:)         = zero
    total_flux_nitrogen(:)       = zero
    total_flux_phosphorus(:)     = zero
    total_flux_carbon13(:)       = zero
    total_flux_carbon14(:)       = zero
    total_flux_nitrogen15(:)     = zero
    leaching_nh4_solute(:)       = zero
    leaching_no3_solute(:)       = zero
    leaching_po4_solute(:)       = zero
    emission_nh3(:)              = zero
    emission_noy(:)              = zero
    emission_n2o(:)              = zero
    emission_n2(:)               = zero

    ! pools
    CALL sb%loss%init(zero)
    CALL sb%formation%init(zero)
    CALL sb%transport%init(zero)

    ! t_slm_struct_sb_element
    CALL sb%mycorrhiza_export%init(zero)

    ! 3D fields
    het_respiration(:,:)              = zero
    het_respiration_c13(:,:)          = zero
    het_respiration_c14(:,:)          = zero
    net_mineralisation_nh4(:,:)       = zero
    net_mineralisation_nh4_n15(:,:)   = zero
    net_mineralisation_no3(:,:)       = zero
    net_mineralisation_no3_n15(:,:)   = zero
    net_mineralisation_po4(:,:)       = zero
    plant_uptake_nh4_sl(:,:)          = zero
    plant_uptake_nh4_n15_sl(:,:)      = zero
    plant_uptake_no3_sl(:,:)          = zero
    plant_uptake_no3_n15_sl(:,:)      = zero
    plant_uptake_po4_sl(:,:)          = zero
    weathering_po4(:,:)               = zero
    occlusion_po4(:,:)                = zero
    slow_sorption_po4(:,:)            = zero
    fast_adsorption_po4(:,:)          = zero
    biochem_mineralisation_po4(:,:)   = zero
    gross_mineralisation_po4(:,:)     = zero
    mycorrhiza_uptake_nh4_sl(:,:)     = zero
    mycorrhiza_uptake_nh4_n15_sl(:,:) = zero
    mycorrhiza_uptake_no3_sl(:,:)     = zero
    mycorrhiza_uptake_no3_n15_sl(:,:) = zero
    mycorrhiza_uptake_norg_sl(:,:)    = zero
    mycorrhiza_uptake_norg_n15_sl(:,:)= zero
    mycorrhiza_uptake_po4_sl(:,:)     = zero
    transport_nh4_solute(:,:)         = zero 
    transport_nh4_assoc(:,:)          = zero 
    transport_no3_solute(:,:)         = zero 
    transport_po4_solute(:,:)         = zero 
    transport_po4_assoc_fast(:,:)     = zero 
    transport_po4_assoc_slow(:,:)     = zero 
    transport_po4_occluded(:,:)       = zero 
    transport_nh4_n15_solute(:,:)     = zero 
    transport_nh4_n15_assoc(:,:)      = zero 
    transport_no3_n15_solute(:,:)     = zero 
    transport_noy(:,:)                = zero 
    transport_n2o(:,:)                = zero 
    transport_n2(:,:)                 = zero 
    transport_noy_n15(:,:)            = zero 
    transport_n2o_n15(:,:)            = zero 
    transport_n2_n15(:,:)             = zero 
    volatilisation_nh4(:,:)           = zero
    nitrification_no3(:,:)            = zero
    nitrification_noy(:,:)            = zero
    nitrification_n2o(:,:)            = zero
    denitrification_noy(:,:)          = zero
    denitrification_n2o(:,:)          = zero
    denitrification_n2(:,:)           = zero
    asymb_n_fixation(:,:)             = zero
    volatilisation_nh4_n15(:,:)       = zero
    nitrification_no3_n15(:,:)        = zero
    nitrification_noy_n15(:,:)        = zero
    nitrification_n2o_n15(:,:)        = zero
    denitrification_noy_n15(:,:)      = zero
    denitrification_n2o_n15(:,:)      = zero
    denitrification_n2_n15(:,:)       = zero
    leaching_nh4_n15_solute(:)        = zero
    leaching_no3_n15_solute(:)        = zero
    emission_nh3_n15(:)               = zero
    emission_noy_n15(:)               = zero
    emission_n2o_n15(:)               = zero
    emission_n2_n15(:)                = zero
    asymb_n_fixation_n15(:,:)         = zero
  
  END SUBROUTINE reset_sb_fluxes




#endif
END MODULE mo_sb_util
