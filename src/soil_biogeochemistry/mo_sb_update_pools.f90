!>
!! @brief Contains the routines for the Task update_sb_pools
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!!
MODULE mo_sb_update_pools

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp, eps8, eps4
  USE mo_exception,       ONLY: message, message_text, finish
  USE mo_base_memory,     ONLY: t_jsb_task_options

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_pools


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_update_pools'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_pools(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_base_config_class,     ONLY: base_config
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
    USE mo_shared_constants,      ONLY: lambda_C14
    USE mo_sb_constants
    USE mo_sb_jsm_processes,      ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, calc_qmax
    USE mo_shared_isotope,        ONLY: calc_mixing_ratio_N15N14

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
 
    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_soil_memory
    dsl4jsb_Use_memory(A2L_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET  :: tile         !< one tile
    TYPE(t_jsb_task_options), INTENT(in)             :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    INTEGER                               :: iblk, ics, ice, nc
    INTEGER                               :: isoil
!!
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: fine_root_carbon,hlp1
!!
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_memory(SB_)                               ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(SOIL_)                             ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(VEG_)                              ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(A2L_)
    ! Declare pointers to variables in memory
    ! A2L
    dsl4jsb_Real2D_onChunk :: nhx_deposition
    dsl4jsb_Real2D_onChunk :: noy_deposition
    dsl4jsb_Real2D_onChunk :: p_deposition
    ! SOIL
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl 
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: total_flux_carbon
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen
    dsl4jsb_Real2D_onChunk :: total_flux_phosphorus
    dsl4jsb_Real2D_onChunk :: total_flux_carbon13
    dsl4jsb_Real2D_onChunk :: total_flux_carbon14
    dsl4jsb_Real2D_onChunk :: total_flux_nitrogen15
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: emission_nh3
    dsl4jsb_Real2D_onChunk :: emission_nh3_n15
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4 
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3 
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_noy 
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy 
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2 
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15 
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_sorption_po4
    dsl4jsb_Real3D_onChunk :: fast_adsorption_po4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: gross_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_nh4_assoc 
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_solute 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow 
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2
    dsl4jsb_Real3D_onChunk :: transport_n2_n15 
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: nh4_assoc 
    dsl4jsb_Real3D_onChunk :: nh4_n15_assoc 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: noy 
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o 
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2
    dsl4jsb_Real3D_onChunk :: n2_n15
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl

! 2 be removed?: pointers to sb *slmtemp
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    TYPE(t_veg_memory), POINTER             :: veg        ! local pointer to veg__mem for SLM usability
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
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id) 
    ! Get process config
    dsl4jsb_Get_config(SB_)
    ! Get process memories
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(A2L_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_,      nhx_deposition)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      noy_deposition)              ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      p_deposition)                ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)               ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     bulk_dens_sl)                ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_phosphorus)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon13)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_carbon14)         ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       total_flux_nitrogen15)       ! out
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute)     ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute)         ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_nh3)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_nh3_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)                ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o_n15)            ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)                 ! in
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2_n15)             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)             ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o_n15)       ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)         ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy_n15)     ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o_n15)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)          ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2_n15)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4)          ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4_n15)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation_n15)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_sl)   ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_n15_sl) ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_assoc)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_assoc)     ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute)        ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy_n15)           ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o_n15)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2_n15)            ! in 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_assoc)                   ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_assoc)               ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute)                  ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       noy)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       noy_n15)                     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o)                         ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o_n15)                     ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2)                          ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       n2_n15)                      ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)               ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_sorption_po4)           ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       fast_adsorption_po4)         ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)  ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       gross_mineralisation_po4)    ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)                    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)                    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)         ! out



! 2 be removed  *slmtemp
    sb   => sb__mem
    veg  => veg__mem


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 0.9 infrastructure for mass balance checks
    !!
    !! fluxes are calculated at the end of this routine (to include het_respiration)
    !!   
    !! changes in pools per time step are calculated here AND at the end of this routine
    sb%delta_pools%carbon%ptr(ics:ice,:,iblk) = &
       sb%pool%dom%carbon%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) + sb%pool%fungi%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%carbon%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%carbon%ptr(ics:ice,:,iblk)
    sb%delta_pools%nitrogen%ptr(ics:ice,:,iblk) = &
       sb%pool%dom%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%fungi%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + &
       nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:) + noy(:,:) + n2o(:,:) + n2(:,:)
    sb%delta_pools%phosphorus%ptr(ics:ice,:,iblk) = &
       sb%pool%dom%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%fungi%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + &
       po4_solute(:,:) + po4_occluded(:,:) + po4_assoc_slow(:,:) + po4_assoc_fast(:,:)
    sb%delta_pools%carbon13%ptr(ics:ice,:,iblk) = &
       sb%pool%dom%carbon13%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) + sb%pool%fungi%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%carbon13%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk)
    ! for mass balance purposes, already reduce by radioactive decay, since this is not a reported flux
    sb%delta_pools%carbon14%ptr(ics:ice,:,iblk) = &
       (sb%pool%dom%carbon14%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) + sb%pool%fungi%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%carbon14%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk)) * &
       (1._wp - lambda_C14 * dtime)
    sb%delta_pools%nitrogen15%ptr(ics:ice,:,iblk) = &
       sb%pool%dom%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%fungi%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + &
       nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + no3_n15_solute(:,:) + noy_n15(:,:) + n2o_n15(:,:) + n2_n15(:,:)


    !> 1.0 loss of C14 due to radioactive decay
    !!
    sb%loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk)  = sb%loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
       sb%pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk)  * lambda_C14 * dtime
    sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb%loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) * lambda_C14 * dtime
    sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)     = sb%loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)     * lambda_C14 * dtime
    sb%loss%dom%carbon14%ptr(ics:ice,:,iblk)              = sb%loss%dom%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%dom%carbon14%ptr(ics:ice,:,iblk)              * lambda_C14 * dtime
    sb%loss%dom_assoc%carbon14%ptr(ics:ice,:,iblk)        = sb%loss%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk)        * lambda_C14 * dtime
    sb%loss%fungi%carbon14%ptr(ics:ice,:,iblk)            = sb%loss%fungi%carbon14%ptr(ics:ice,:,iblk) + & 
       sb%pool%fungi%carbon14%ptr(ics:ice,:,iblk)            * lambda_C14 * dtime
    sb%loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)       = sb%loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + & 
       sb%pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)       * lambda_C14 * dtime
    sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk)        = sb%loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
       sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk)        * lambda_C14 * dtime
    sb%loss%residue%carbon14%ptr(ics:ice,:,iblk)          = sb%loss%residue%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%carbon14%ptr(ics:ice,:,iblk)          * lambda_C14 * dtime
    sb%loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk)    = sb%loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
       sb%pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk)    * lambda_C14 * dtime 


    !> 2.0 update mycorrhiza with flux due to exudation from vegetation and N uptake from soil AND correct flux due to export to plants
    !!
    sb%formation%mycorrhiza%carbon%ptr(ics:ice,1,iblk)     = sb%formation%mycorrhiza%carbon%ptr(ics:ice,1,iblk) + &
                                                             veg%exudation%carbon%ptr * mycorrhizal_cue / soil_depth_sl(:,1)
    sb%formation%mycorrhiza%nitrogen%ptr(ics:ice,1,iblk)   = sb%formation%mycorrhiza%nitrogen%ptr(ics:ice,1,iblk) + &
                                                             veg%exudation%nitrogen%ptr / soil_depth_sl(:,1) + &
                                                             (mycorrhiza_uptake_nh4_sl(:,1) + mycorrhiza_uptake_no3_sl(:,1)) * & !+ mycorrhiza_uptake_norg_sl(:,1)) * &
                                                             dtime / 1000000._wp
    sb%formation%mycorrhiza%phosphorus%ptr(ics:ice,1,iblk) = sb%formation%mycorrhiza%phosphorus%ptr(ics:ice,1,iblk) + &
                                                             veg%exudation%phosphorus%ptr / soil_depth_sl(:,1) + &
                                                             mycorrhiza_uptake_po4_sl(:,1) * dtime / 1000000._wp
    sb%formation%mycorrhiza%carbon13%ptr(ics:ice,1,iblk)   = sb%formation%mycorrhiza%carbon13%ptr(ics:ice,1,iblk) + &
                                                             veg%exudation%carbon13%ptr * mycorrhizal_cue / soil_depth_sl(:,1)
    sb%formation%mycorrhiza%carbon14%ptr(ics:ice,1,iblk)   = sb%formation%mycorrhiza%carbon14%ptr(ics:ice,1,iblk) + &
                                                             veg%exudation%carbon14%ptr * mycorrhizal_cue / soil_depth_sl(:,1)
    sb%formation%mycorrhiza%nitrogen15%ptr(ics:ice,1,iblk) = sb%formation%mycorrhiza%nitrogen15%ptr(ics:ice,1,iblk) + &
                                                             veg%exudation%nitrogen15%ptr / soil_depth_sl(:,1) + &
                                                             (mycorrhiza_uptake_nh4_n15_sl(:,1) + mycorrhiza_uptake_no3_n15_sl(:,1)) * & !+ mycorrhiza_uptake_norg_n15_sl(:,1)) * &
                                                             dtime / 1000000._wp

    sb%loss%mycorrhiza%carbon%ptr(ics:ice,1,iblk)     = sb%loss%mycorrhiza%carbon%ptr(ics:ice,1,iblk) + &
                                                        sb%mycorrhiza_export%carbon%ptr(ics:ice,1,iblk) 
    sb%loss%mycorrhiza%nitrogen%ptr(ics:ice,1,iblk)   = sb%loss%mycorrhiza%nitrogen%ptr(ics:ice,1,iblk) + &
                                                        sb%mycorrhiza_export%nitrogen%ptr(ics:ice,1,iblk)
    sb%loss%mycorrhiza%phosphorus%ptr(ics:ice,1,iblk) = sb%loss%mycorrhiza%phosphorus%ptr(ics:ice,1,iblk) + &
                                                        sb%mycorrhiza_export%phosphorus%ptr(ics:ice,1,iblk)
    sb%loss%mycorrhiza%carbon13%ptr(ics:ice,1,iblk)   = sb%loss%mycorrhiza%carbon13%ptr(ics:ice,1,iblk) + &
                                                        sb%mycorrhiza_export%carbon13%ptr(ics:ice,1,iblk)
    sb%loss%mycorrhiza%carbon14%ptr(ics:ice,1,iblk)   = sb%loss%mycorrhiza%carbon14%ptr(ics:ice,1,iblk) + &
                                                        sb%mycorrhiza_export%carbon14%ptr(ics:ice,1,iblk)
    sb%loss%mycorrhiza%nitrogen15%ptr(ics:ice,1,iblk) = sb%loss%mycorrhiza%nitrogen15%ptr(ics:ice,1,iblk) + &
                                                        sb%mycorrhiza_export%nitrogen15%ptr(ics:ice,1,iblk)

    het_respiration(:,:)     = het_respiration(:,:)     + (1._wp - mycorrhizal_cue) * veg%exudation%carbon%ptr / dtime * 1000000._wp                          
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - mycorrhizal_cue) * veg%exudation%carbon13%ptr / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - mycorrhizal_cue) * veg%exudation%carbon14%ptr / dtime * 1000000._wp

    !> 3.0 update of organic pools
    !!
    CALL calc_organic_pool_change(options,sb%pool,sb%formation,sb%loss,sb%transport)

    !> 3.1 update diagnostics and other state variables depending on organic matter content 
    !!
    bulk_soil_carbon_sl = calc_bulk_soil_carbon(options, sb%pool)
    bulk_dens_corr_sl = calc_bulk_density_correction(bulk_soil_carbon_sl,bulk_dens_sl)
    qmax_org = calc_qmax(bulk_dens_corr_sl,qmax_org_mineral_default)
    qmax_po4 = calc_qmax(bulk_dens_corr_sl,qmax_po4_mineral)

    !> 4.0 update of inorganic pools
    !!
    !> 4.1 update of N pools  
    !!
    nh4_solute(:,:)     = nh4_solute(:,:)     + (net_mineralisation_nh4(:,:) + asymb_n_fixation(:,:) + transport_nh4_solute(:,:) - & 
                                                 plant_uptake_nh4_sl(:,:) - mycorrhiza_uptake_nh4_sl(:,:) - & 
                                                 nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:) - &
                                                 volatilisation_nh4(:,:)) * dtime/1.e6_wp
    nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + (net_mineralisation_nh4_n15(:,:) + asymb_n_fixation_n15(:,:) + transport_nh4_n15_solute(:,:) - & 
                                                 plant_uptake_nh4_n15_sl(:,:) - mycorrhiza_uptake_nh4_n15_sl(:,:) - &
                                                 nitrification_no3_n15(:,:) - nitrification_noy_n15(:,:) - nitrification_n2o_n15(:,:) - &
                                                 volatilisation_nh4_n15(:,:)) * dtime/1.e6_wp
    no3_solute(:,:)     = no3_solute(:,:)     + (net_mineralisation_no3(:,:) + nitrification_no3(:,:) + transport_no3_solute(:,:) - & 
                                                 plant_uptake_no3_sl(:,:) - mycorrhiza_uptake_no3_sl(:,:) - &
                                                 denitrification_noy(:,:) - denitrification_n2o(:,:) - denitrification_n2(:,:)) * dtime/1.e6_wp
    no3_n15_solute(:,:) = no3_n15_solute(:,:) + (net_mineralisation_no3_n15(:,:) + nitrification_no3_n15(:,:) + transport_no3_n15_solute(:,:) - & 
                                                 plant_uptake_no3_n15_sl(:,:) - mycorrhiza_uptake_no3_n15_sl(:,:) + &
                                                 denitrification_noy_n15(:,:) - denitrification_n2o_n15(:,:) - denitrification_n2_n15(:,:)) * dtime/1.e6_wp

    noy(:,:)            = noy(:,:)            + (nitrification_noy(:,:) + denitrification_noy(:,:) + & 
                                                 transport_noy(:,:)) * dtime/1.e6_wp
    noy_n15(:,:)        = noy_n15(:,:)        + (nitrification_noy_n15(:,:) + denitrification_noy_n15(:,:) + & 
                                                 transport_noy_n15(:,:)) * dtime/1.e6_wp
    n2o(:,:)            = n2o(:,:)            + (nitrification_n2o(:,:) + denitrification_n2o(:,:) + & 
                                                 transport_n2o(:,:)) * dtime/1.e6_wp
    n2o_n15(:,:)        = n2o_n15(:,:)        + (nitrification_n2o_n15(:,:) + denitrification_n2o_n15(:,:) + & 
                                                 transport_n2o_n15(:,:)) * dtime/1.e6_wp
    n2(:,:)             = n2(:,:)             + (denitrification_n2(:,:) + transport_n2(:,:)) * dtime/1.e6_wp
    n2_n15(:,:)         = n2_n15(:,:)         + (denitrification_n2_n15(:,:) + transport_n2_n15(:,:)) * dtime/1.e6_wp

    !nh4_assoc(:,:) = nh4_assoc(:,:) + transport_nh4_assoc(:,:)
    !po4_assoc(:,:) = po4_assoc(:,:) + transport_po4_assoc(:,:)
    !nh4_n15_assoc(:,:) = nh4_n15_assoc(:,:) + transport_nh4_n15_assoc(:,:)

    !> 4.2 update of mineral P pools
    !!
    po4_solute(:,:)     = po4_solute(:,:)     + (net_mineralisation_po4(:,:) - plant_uptake_po4_sl(:,:) - mycorrhiza_uptake_po4_sl(:,:) + weathering_po4(:,:) &
                          + biochem_mineralisation_po4(:,:) + transport_po4_solute(:,:) - fast_adsorption_po4(:,:)) * dtime/1.e6_wp
    po4_assoc_fast(:,:) = po4_assoc_fast(:,:) + (fast_adsorption_po4 - slow_sorption_po4 + transport_po4_assoc_fast) * dtime/1.e6_wp
    po4_assoc_slow(:,:) = po4_assoc_slow(:,:) + (slow_sorption_po4 - occlusion_po4 + transport_po4_assoc_slow) * dtime/1.e6_wp
    po4_occluded(:,:)   = po4_occluded(:,:)   + (occlusion_po4 + transport_po4_occluded) * dtime/1.e6_wp


    !> 4.3 atmospheric deposition
    !!
    nh4_solute(:,1)     = nh4_solute(:,1)     + nhx_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
    no3_solute(:,1)     = no3_solute(:,1)     + noy_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
    po4_solute(:,1)     = po4_solute(:,1)     + p_deposition(:)   / 1.e6_wp / soil_depth_sl(:,1) * dtime
    nh4_n15_solute(:,1) = nh4_n15_solute(:,1) + nhx_deposition(:) / 1.e6_wp / & 
                                                ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) / soil_depth_sl(:,1) * dtime
    no3_n15_solute(:,1) = no3_n15_solute(:,1) + noy_deposition(:) / 1.e6_wp / & 
                                                ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) / soil_depth_sl(:,1) * dtime

    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_nh4)THEN
       hlp1(:,:) = nh4_solute_prescribe - nh4_solute(:,:)
       nh4_solute(:,:) = nh4_solute(:,:) + hlp1(:,:) 
       nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    ENDIF
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_no3)THEN
       hlp1(:,:) = no3_solute_prescribe - no3_solute(:,:)
       no3_solute(:,:) = no3_solute(:,:) + hlp1(:,:) 
       no3_n15_solute(:,:) = no3_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    ENDIF
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_po4)THEN
       po4_solute(:,:) = po4_solute_prescribe 
    ENDIF

    !> 5.0 infrastructure for mass balance checks (see also section 0.9)
    !!
    
    !! pools: calc the difference to what was calculated at the beginning of this routine
    sb%delta_pools%carbon%ptr(ics:ice,:,iblk) = &
       (sb%pool%dom%carbon%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%carbon%ptr(ics:ice,:,iblk) + sb%pool%fungi%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + sb%pool%microbial%carbon%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%carbon%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%carbon%ptr(ics:ice,:,iblk)) - & 
       sb%delta_pools%carbon%ptr(ics:ice,:,iblk)
    sb%delta_pools%nitrogen%ptr(ics:ice,:,iblk) = &
       (sb%pool%dom%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%fungi%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%microbial%nitrogen%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%nitrogen%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + &
       nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:) + noy(:,:) + n2o(:,:) + n2(:,:)) - &
       sb%delta_pools%nitrogen%ptr(ics:ice,:,iblk)
    sb%delta_pools%phosphorus%ptr(ics:ice,:,iblk) = &
       (sb%pool%dom%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%fungi%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%microbial%phosphorus%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + &
       po4_solute(:,:) + po4_occluded(:,:) + po4_assoc_slow(:,:) + po4_assoc_fast(:,:)) - &
       sb%delta_pools%phosphorus%ptr(ics:ice,:,iblk)
    sb%delta_pools%carbon13%ptr(ics:ice,:,iblk) = &
       (sb%pool%dom%carbon13%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) + sb%pool%fungi%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + sb%pool%microbial%carbon13%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%carbon13%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk)) - &
       sb%delta_pools%carbon13%ptr(ics:ice,:,iblk)
    sb%delta_pools%carbon14%ptr(ics:ice,:,iblk) = &
       (sb%pool%dom%carbon14%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) + sb%pool%fungi%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + sb%pool%microbial%carbon14%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%carbon14%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk)) - &
       sb%delta_pools%carbon14%ptr(ics:ice,:,iblk)
    sb%delta_pools%nitrogen15%ptr(ics:ice,:,iblk) = &
       (sb%pool%dom%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%fungi%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) + &
       sb%pool%residue%nitrogen15%ptr(ics:ice,:,iblk) + sb%pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + &
       nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + no3_n15_solute(:,:) + noy_n15(:,:) + n2o_n15(:,:) + n2_n15(:,:)) - &
       sb%delta_pools%nitrogen15%ptr(ics:ice,:,iblk)

    !! fluxes

    total_flux_carbon(:) =                  veg%litterfall%leaf%carbon%ptr        + &
      veg%litterfall%fine_root%carbon%ptr + veg%litterfall%coarse_root%carbon%ptr + &
      veg%litterfall%sap_wood%carbon%ptr  + veg%litterfall%heart_wood%carbon%ptr  + &
      veg%litterfall%labile%carbon%ptr    + veg%litterfall%reserve%carbon%ptr     + &
      veg%litterfall%fruit%carbon%ptr     + veg%litterfall%seed_bed%carbon%ptr    + &
      veg%exudation%carbon%ptr                                                    - &
      SUM(sb%mycorrhiza_export%carbon%ptr(:,:,iblk) * soil_depth_sl(:,:), DIM=2)  - &
      SUM(het_respiration(:,:) * soil_depth_sl(:,:), DIM=2) * dtime/1.e6_wp                          ! DIM=2: sum only over dimension 2, i.e., soil layers
    total_flux_nitrogen(:) =                  veg%litterfall%leaf%nitrogen%ptr        + &
      veg%litterfall%fine_root%nitrogen%ptr + veg%litterfall%coarse_root%nitrogen%ptr + &
      veg%litterfall%sap_wood%nitrogen%ptr  + veg%litterfall%heart_wood%nitrogen%ptr  + &
      veg%litterfall%labile%nitrogen%ptr    + veg%litterfall%reserve%nitrogen%ptr     + &
      veg%litterfall%fruit%nitrogen%ptr     + veg%litterfall%seed_bed%nitrogen%ptr    + &
      SUM(asymb_n_fixation(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp           + &
      veg%exudation%nitrogen%ptr                                                      - &
      SUM(sb%mycorrhiza_export%nitrogen%ptr(:,:,iblk) * soil_depth_sl(:,:), DIM=2)    - &
      (leaching_nh4_solute(:) + leaching_no3_solute(:)                                + &
       emission_nh3(:) + emission_noy(:) + emission_n2o(:) + emission_n2(:)            + &
       SUM(plant_uptake_nh4_sl(:,:) * soil_depth_sl(:,:), DIM=2)                       + &
       SUM(plant_uptake_no3_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp                       
    total_flux_phosphorus(:) =                  veg%litterfall%leaf%phosphorus%ptr        + &
      veg%litterfall%fine_root%phosphorus%ptr + veg%litterfall%coarse_root%phosphorus%ptr + &
      veg%litterfall%sap_wood%phosphorus%ptr  + veg%litterfall%heart_wood%phosphorus%ptr  + &
      veg%litterfall%labile%phosphorus%ptr    + veg%litterfall%reserve%phosphorus%ptr     + &
      veg%litterfall%fruit%phosphorus%ptr     + veg%litterfall%seed_bed%phosphorus%ptr    + &
      veg%exudation%phosphorus%ptr                                                        - &
      SUM(sb%mycorrhiza_export%phosphorus%ptr(:,:,iblk) * soil_depth_sl(:,:), DIM=2)      - &
      (leaching_po4_solute(:)                                                             + &
       SUM(plant_uptake_po4_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp       
    total_flux_carbon13(:) =                  veg%litterfall%leaf%carbon13%ptr        + &
      veg%litterfall%fine_root%carbon13%ptr + veg%litterfall%coarse_root%carbon13%ptr + &
      veg%litterfall%sap_wood%carbon13%ptr  + veg%litterfall%heart_wood%carbon13%ptr  + &
      veg%litterfall%labile%carbon13%ptr    + veg%litterfall%reserve%carbon13%ptr     + &
      veg%litterfall%fruit%carbon13%ptr     + veg%litterfall%seed_bed%carbon13%ptr    + &
      veg%exudation%carbon13%ptr                                                      - &
      SUM(sb%mycorrhiza_export%carbon13%ptr(:,:,iblk) * soil_depth_sl(:,:), DIM=2)    - &
      SUM(het_respiration_c13(:,:) * soil_depth_sl(:,:), DIM=2) * dtime/1.e6_wp
    total_flux_carbon14(:) =                  veg%litterfall%leaf%carbon14%ptr        + &
      veg%litterfall%fine_root%carbon14%ptr + veg%litterfall%coarse_root%carbon14%ptr + &
      veg%litterfall%sap_wood%carbon14%ptr  + veg%litterfall%heart_wood%carbon14%ptr  + &
      veg%litterfall%labile%carbon14%ptr    + veg%litterfall%reserve%carbon14%ptr     + &
      veg%litterfall%fruit%carbon14%ptr     + veg%litterfall%seed_bed%carbon14%ptr    + &
      veg%exudation%carbon14%ptr                                                      - &
      SUM(sb%mycorrhiza_export%carbon14%ptr(:,:,iblk) * soil_depth_sl(:,:), DIM=2)    - &
      SUM(het_respiration_c14(:,:) * soil_depth_sl(:,:), DIM=2) * dtime/1.e6_wp
    total_flux_nitrogen15(:) =                  veg%litterfall%leaf%nitrogen15%ptr          + &
      veg%litterfall%fine_root%nitrogen15%ptr + veg%litterfall%coarse_root%nitrogen15%ptr   + &
      veg%litterfall%sap_wood%nitrogen15%ptr  + veg%litterfall%heart_wood%nitrogen15%ptr    + &
      veg%litterfall%labile%nitrogen15%ptr    + veg%litterfall%reserve%nitrogen15%ptr       + &
      veg%litterfall%fruit%nitrogen15%ptr     + veg%litterfall%seed_bed%nitrogen15%ptr      + &
      SUM(asymb_n_fixation_n15(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp             + &
      veg%exudation%nitrogen15%ptr                                                          - &
      SUM(sb%mycorrhiza_export%nitrogen15%ptr(:,:,iblk) * soil_depth_sl(:,:), DIM=2)        - &
      (leaching_nh4_n15_solute(:) - leaching_no3_n15_solute(:)                              - &
       emission_nh3_n15(:) - emission_noy_n15(:) - emission_n2o_n15(:) - emission_n2_n15(:) + &
       SUM(plant_uptake_nh4_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)                        + &
       SUM(plant_uptake_no3_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp  


  END SUBROUTINE update_sb_pools


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> updates all SB pools given the JSM equations for formation, loss and transport
  !! should be a simple pool operator overload function, but that doesn't work yet
  !!
  !! Input: formation, loss and transport rate, pool sizes
  !!
  !! Output: updated pool sizes
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_organic_pool_change(options,pool,formation,loss,transport)

    USE mo_sb_memory_class,       ONLY: t_slm_struct_sb_pool, t_slm_struct_sb_element

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options), INTENT(in)          :: options      !< model options
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: formation    !< rate of formation of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: loss         !< rate of loss of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(in)       :: transport    !< rate of transport of a pool [mol/m3/timestep]
    TYPE(t_slm_struct_sb_pool),  INTENT(inout)    :: pool         !< the pool [mol/m3]
    ! ---------------------------
    ! 0.2 Local 
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_organic_pool_change'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! carbon
    pool%dom%carbon%ptr(ics:ice,:,iblk) = pool%dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon%ptr(ics:ice,:,iblk)  
    pool%fungi%carbon%ptr(ics:ice,:,iblk) = pool%fungi%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%carbon%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) = pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%carbon%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon%ptr(ics:ice,:,iblk) = pool%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon%ptr(ics:ice,:,iblk)  
    pool%residue%carbon%ptr(ics:ice,:,iblk) = pool%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon%ptr(ics:ice,:,iblk)  

    ! nitrogen
    pool%dom%nitrogen%ptr(ics:ice,:,iblk) = pool%dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%nitrogen%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) = pool%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%nitrogen%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%nitrogen%ptr(ics:ice,:,iblk)  
    pool%fungi%nitrogen%ptr(ics:ice,:,iblk) = pool%fungi%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%nitrogen%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) = pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk)  
    pool%microbial%nitrogen%ptr(ics:ice,:,iblk) = pool%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%nitrogen%ptr(ics:ice,:,iblk)  
    pool%residue%nitrogen%ptr(ics:ice,:,iblk) = pool%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%nitrogen%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) = pool%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%nitrogen%ptr(ics:ice,:,iblk)  

    ! phosphorus
    pool%dom%phosphorus%ptr(ics:ice,:,iblk) = pool%dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%phosphorus%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) = pool%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%phosphorus%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%phosphorus%ptr(ics:ice,:,iblk)  
    pool%fungi%phosphorus%ptr(ics:ice,:,iblk) = pool%fungi%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%phosphorus%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) = pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk)  
    pool%microbial%phosphorus%ptr(ics:ice,:,iblk) = pool%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%phosphorus%ptr(ics:ice,:,iblk)  
    pool%residue%phosphorus%ptr(ics:ice,:,iblk) = pool%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%phosphorus%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) = pool%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%phosphorus%ptr(ics:ice,:,iblk)  

    ! C13
    pool%dom%carbon13%ptr(ics:ice,:,iblk) = pool%dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon13%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon13%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon13%ptr(ics:ice,:,iblk)  
    pool%fungi%carbon13%ptr(ics:ice,:,iblk) = pool%fungi%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%carbon13%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) = pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%carbon13%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon13%ptr(ics:ice,:,iblk) = pool%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon13%ptr(ics:ice,:,iblk)  
    pool%residue%carbon13%ptr(ics:ice,:,iblk) = pool%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon13%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon13%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon13%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon13%ptr(ics:ice,:,iblk)  
    ! C14
    pool%dom%carbon14%ptr(ics:ice,:,iblk) = pool%dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%carbon14%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) = pool%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%carbon14%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) = pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) = pool%woody_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%carbon14%ptr(ics:ice,:,iblk)  
    pool%fungi%carbon14%ptr(ics:ice,:,iblk) = pool%fungi%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%carbon14%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) = pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)  
    pool%microbial%carbon14%ptr(ics:ice,:,iblk) = pool%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%carbon14%ptr(ics:ice,:,iblk)  
    pool%residue%carbon14%ptr(ics:ice,:,iblk) = pool%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%carbon14%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk) = pool%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%carbon14%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%carbon14%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%carbon14%ptr(ics:ice,:,iblk)  
    ! N15
    pool%dom%nitrogen15%ptr(ics:ice,:,iblk) = pool%dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%dom%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%dom%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%dom%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) = pool%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%dom_assoc%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%fungi%nitrogen15%ptr(ics:ice,:,iblk) = pool%fungi%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%fungi%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%fungi%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%fungi%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) = pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) = pool%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%microbial%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%residue%nitrogen15%ptr(ics:ice,:,iblk) = pool%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%residue%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%residue%nitrogen15%ptr(ics:ice,:,iblk)  
    pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) = pool%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                formation%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) - & 
                                                loss%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                transport%residue_assoc%nitrogen15%ptr(ics:ice,:,iblk)  

  END SUBROUTINE calc_organic_pool_change


END MODULE mo_sb_update_pools
