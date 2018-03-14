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
  USE mo_jsb_control,     ONLY: debug_on
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
    USE mo_sb_jsm_processes,      ONLY: calc_bulk_soil_carbon, calc_bulk_density_correction, calc_qmax, calc_adsorption_po4
    USE mo_shared_isotope,        ONLY: calc_mixing_ratio_N15N14

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
 
    ! Use of process memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool, sb_formation, sb_loss, sb_transport, sb_delta_pools, sb_mycorrhiza_export
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_soil_memory
    dsl4slm_Use_memory_a_veg
    USE mo_veg_memory_class,      ONLY: veg_litterfall, veg_exudation, veg_pool
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
    INTEGER                               :: nsoil_sb
    INTEGER                               :: isoil
!!
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: fine_root_carbon, p_deposition_sl, hlp1
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
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl
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
    dsl4jsb_Real3D_onChunk :: microbial_uptake_po4_sl
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
    dsl4jsb_Real3D_onChunk :: transport_po4_primary
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
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: noy 
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o 
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2
    dsl4jsb_Real3D_onChunk :: n2_n15
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: qmax_org
    dsl4jsb_Real3D_onChunk :: qmax_po4
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl
    dsl4jsb_Real3D_onChunk :: bulk_soil_carbon_sl
    ! Get local variables from options argument
    iblk      = options%iblk
    ics       = options%ics
    ice       = options%ice
    nc        = options%nc
    nsoil_sb  = options%nsoil_sb
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(SB_)) RETURN
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
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
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)              ! in
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
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)     ! in
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
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_primary)       ! in
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
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)                 ! inout
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
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)                ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)              ! in
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_org)                    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)                    ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_soil_carbon_sl)         ! out
    ! Get process memory pool structures
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)
#define Get_veg_litterfall
#define Get_veg_exudation
#define Get_sb_pool
#define Get_sb_formation
#define Get_sb_loss
#define Get_sb_transport
#define Get_sb_delta_pools
#define Get_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Get_veg_litterfall
#undef Get_veg_exudation
#undef Get_sb_pool
#undef Get_sb_formation
#undef Get_sb_loss
#undef Get_sb_transport
#undef Get_sb_delta_pools
#undef Get_sb_mycorrhiza_export


    p_deposition_sl(:,:) = 0.0_wp
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 0.9 infrastructure for mass balance checks
    !!
    !! fluxes are calculated at the end of this routine (to include het_respiration)
    !!   
    !! changes in pools per time step are calculated here AND at the end of this routine
    sb_delta_pools%carbon(:,:) = &
       sb_pool%dom%carbon(:,:) + sb_pool%dom_assoc%carbon(:,:) + &
       sb_pool%soluable_litter%carbon(:,:) + sb_pool%polymeric_litter%carbon(:,:) + &
       sb_pool%woody_litter%carbon(:,:) + sb_pool%fungi%carbon(:,:) + &
       sb_pool%mycorrhiza%carbon(:,:) + sb_pool%microbial%carbon(:,:) + &
       sb_pool%residue%carbon(:,:) + sb_pool%residue_assoc%carbon(:,:)
    sb_delta_pools%nitrogen(:,:) = &
       sb_pool%dom%nitrogen(:,:) + sb_pool%dom_assoc%nitrogen(:,:) + &
       sb_pool%soluable_litter%nitrogen(:,:) + sb_pool%polymeric_litter%nitrogen(:,:) + &
       sb_pool%woody_litter%nitrogen(:,:) + sb_pool%fungi%nitrogen(:,:) + &
       sb_pool%mycorrhiza%nitrogen(:,:) + sb_pool%microbial%nitrogen(:,:) + &
       sb_pool%residue%nitrogen(:,:) + sb_pool%residue_assoc%nitrogen(:,:) + &
       nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:) + noy(:,:) + n2o(:,:) + n2(:,:)
    sb_delta_pools%phosphorus(:,:) = &
       sb_pool%dom%phosphorus(:,:) + sb_pool%dom_assoc%phosphorus(:,:) + &
       sb_pool%soluable_litter%phosphorus(:,:) + sb_pool%polymeric_litter%phosphorus(:,:) + &
       sb_pool%woody_litter%phosphorus(:,:) + sb_pool%fungi%phosphorus(:,:) + &
       sb_pool%mycorrhiza%phosphorus(:,:) + sb_pool%microbial%phosphorus(:,:) + &
       sb_pool%residue%phosphorus(:,:) + sb_pool%residue_assoc%phosphorus(:,:) + &
       po4_solute(:,:) + po4_occluded(:,:) + po4_assoc_slow(:,:) + po4_assoc_fast(:,:) + po4_primary(:,:)
    sb_delta_pools%carbon13(:,:) = &
       sb_pool%dom%carbon13(:,:) + sb_pool%dom_assoc%carbon13(:,:) + &
       sb_pool%soluable_litter%carbon13(:,:) + sb_pool%polymeric_litter%carbon13(:,:) + &
       sb_pool%woody_litter%carbon13(:,:) + sb_pool%fungi%carbon13(:,:) + &
       sb_pool%mycorrhiza%carbon13(:,:) + sb_pool%microbial%carbon13(:,:) + &
       sb_pool%residue%carbon13(:,:) + sb_pool%residue_assoc%carbon13(:,:)
    ! for mass balance purposes, already reduce by radioactive decay, since this is not a reported flux
    sb_delta_pools%carbon14(:,:) = &
       (sb_pool%dom%carbon14(:,:) + sb_pool%dom_assoc%carbon14(:,:) + &
       sb_pool%soluable_litter%carbon14(:,:) + sb_pool%polymeric_litter%carbon14(:,:) + &
       sb_pool%woody_litter%carbon14(:,:) + sb_pool%fungi%carbon14(:,:) + &
       sb_pool%mycorrhiza%carbon14(:,:) + sb_pool%microbial%carbon14(:,:) + &
       sb_pool%residue%carbon14(:,:) + sb_pool%residue_assoc%carbon14(:,:)) * &
       (1._wp - lambda_C14 * dtime)
    sb_delta_pools%nitrogen15(:,:) = &
       sb_pool%dom%nitrogen15(:,:) + sb_pool%dom_assoc%nitrogen15(:,:) + &
       sb_pool%soluable_litter%nitrogen15(:,:) + sb_pool%polymeric_litter%nitrogen15(:,:) + &
       sb_pool%woody_litter%nitrogen15(:,:) + sb_pool%fungi%nitrogen15(:,:) + &
       sb_pool%mycorrhiza%nitrogen15(:,:) + sb_pool%microbial%nitrogen15(:,:) + &
       sb_pool%residue%nitrogen15(:,:) + sb_pool%residue_assoc%nitrogen15(:,:) + &
       nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + no3_n15_solute(:,:) + noy_n15(:,:) + n2o_n15(:,:) + n2_n15(:,:)


    !> 1.0 loss of C14 due to radioactive decay
    !!
    sb_loss%soluable_litter%carbon14(:,:)  = sb_loss%soluable_litter%carbon14(:,:) + & 
       sb_pool%soluable_litter%carbon14(:,:)  * lambda_C14 * dtime
    sb_loss%polymeric_litter%carbon14(:,:) = sb_loss%polymeric_litter%carbon14(:,:) + &
       sb_pool%polymeric_litter%carbon14(:,:) * lambda_C14 * dtime
    sb_loss%woody_litter%carbon14(:,:)     = sb_loss%woody_litter%carbon14(:,:) + &
       sb_pool%woody_litter%carbon14(:,:)     * lambda_C14 * dtime
    sb_loss%dom%carbon14(:,:)              = sb_loss%dom%carbon14(:,:) + &
       sb_pool%dom%carbon14(:,:)              * lambda_C14 * dtime
    sb_loss%dom_assoc%carbon14(:,:)        = sb_loss%dom_assoc%carbon14(:,:) + &
       sb_pool%dom_assoc%carbon14(:,:)        * lambda_C14 * dtime
    sb_loss%fungi%carbon14(:,:)            = sb_loss%fungi%carbon14(:,:) + & 
       sb_pool%fungi%carbon14(:,:)            * lambda_C14 * dtime
    sb_loss%mycorrhiza%carbon14(:,:)       = sb_loss%mycorrhiza%carbon14(:,:) + & 
       sb_pool%mycorrhiza%carbon14(:,:)       * lambda_C14 * dtime
    sb_loss%microbial%carbon14(:,:)        = sb_loss%microbial%carbon14(:,:) + & 
       sb_pool%microbial%carbon14(:,:)        * lambda_C14 * dtime
    sb_loss%residue%carbon14(:,:)          = sb_loss%residue%carbon14(:,:) + &
       sb_pool%residue%carbon14(:,:)          * lambda_C14 * dtime
    sb_loss%residue_assoc%carbon14(:,:)    = sb_loss%residue_assoc%carbon14(:,:) + & 
       sb_pool%residue_assoc%carbon14(:,:)    * lambda_C14 * dtime 


    !> 2.0 update mycorrhiza with flux due to exudation from vegetation and N uptake from soil AND correct flux due to export to plants
    !!
    sb_formation%mycorrhiza%carbon(:,1)     = sb_formation%mycorrhiza%carbon(:,1) + &
                                                             veg_exudation%carbon * mycorrhizal_cue / soil_depth_sl(:,1)
    sb_formation%mycorrhiza%nitrogen(:,1)   = sb_formation%mycorrhiza%nitrogen(:,1) + &
                                                             veg_exudation%nitrogen / soil_depth_sl(:,1) + &
                                                             (mycorrhiza_uptake_nh4_sl(:,1) + mycorrhiza_uptake_no3_sl(:,1)) * & !+ mycorrhiza_uptake_norg_sl(:,1)) * &
                                                             dtime / 1000000._wp
    sb_formation%mycorrhiza%phosphorus(:,1) = sb_formation%mycorrhiza%phosphorus(:,1) + &
                                                             veg_exudation%phosphorus / soil_depth_sl(:,1) + &
                                                             mycorrhiza_uptake_po4_sl(:,1) * dtime / 1000000._wp
    sb_formation%mycorrhiza%carbon13(:,1)   = sb_formation%mycorrhiza%carbon13(:,1) + &
                                                             veg_exudation%carbon13 * mycorrhizal_cue / soil_depth_sl(:,1)
    sb_formation%mycorrhiza%carbon14(:,1)   = sb_formation%mycorrhiza%carbon14(:,1) + &
                                                             veg_exudation%carbon14 * mycorrhizal_cue / soil_depth_sl(:,1)
    sb_formation%mycorrhiza%nitrogen15(:,1) = sb_formation%mycorrhiza%nitrogen15(:,1) + &
                                                             veg_exudation%nitrogen15 / soil_depth_sl(:,1) + &
                                                             (mycorrhiza_uptake_nh4_n15_sl(:,1) + mycorrhiza_uptake_no3_n15_sl(:,1)) * & !+ mycorrhiza_uptake_norg_n15_sl(:,1)) * &
                                                             dtime / 1000000._wp

    sb_loss%mycorrhiza%carbon(:,1)     = sb_loss%mycorrhiza%carbon(:,1) + &
                                                        sb_mycorrhiza_export%carbon(:,1) 
    sb_loss%mycorrhiza%nitrogen(:,1)   = sb_loss%mycorrhiza%nitrogen(:,1) + &
                                                        sb_mycorrhiza_export%nitrogen(:,1)
    sb_loss%mycorrhiza%phosphorus(:,1) = sb_loss%mycorrhiza%phosphorus(:,1) + &
                                                        sb_mycorrhiza_export%phosphorus(:,1)
    sb_loss%mycorrhiza%carbon13(:,1)   = sb_loss%mycorrhiza%carbon13(:,1) + &
                                                        sb_mycorrhiza_export%carbon13(:,1)
    sb_loss%mycorrhiza%carbon14(:,1)   = sb_loss%mycorrhiza%carbon14(:,1) + &
                                                        sb_mycorrhiza_export%carbon14(:,1)
    sb_loss%mycorrhiza%nitrogen15(:,1) = sb_loss%mycorrhiza%nitrogen15(:,1) + &
                                                        sb_mycorrhiza_export%nitrogen15(:,1)

    DO isoil = 1,nsoil_sb
      het_respiration(:,isoil)     = het_respiration(:,isoil)     + (1._wp - mycorrhizal_cue) * veg_exudation%carbon(:)   / &
                                      dtime * 1000000._wp                          
      het_respiration_c13(:,isoil) = het_respiration_c13(:,isoil) + (1._wp - mycorrhizal_cue) * veg_exudation%carbon13(:) / &
                                      dtime * 1000000._wp
      het_respiration_c14(:,isoil) = het_respiration_c14(:,isoil) + (1._wp - mycorrhizal_cue) * veg_exudation%carbon14(:) / &
                                      dtime * 1000000._wp
    END DO

    !> 3.0 update of organic pools
    !!
    !! updates all SB pools given the JSM equations for formation, loss and transport
    !! Output: updated pool sizes
    sb_pool = sb_pool + sb_formation - sb_loss + sb_transport 


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

    !> 4.2 atmospheric deposition
    !!
    nh4_solute(:,1)      = nh4_solute(:,1)     + nhx_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
    no3_solute(:,1)      = no3_solute(:,1)     + noy_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime
    po4_solute(:,1)      = po4_solute(:,1)     + p_deposition(:)   / 1.e6_wp / soil_depth_sl(:,1) * dtime
    p_deposition_sl(:,1) = p_deposition(:) / soil_depth_sl(:,1) !µmol P / m3 / s
    nh4_n15_solute(:,1)  = nh4_n15_solute(:,1) + nhx_deposition(:) / 1.e6_wp / & 
                                                ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) / soil_depth_sl(:,1) * dtime
    no3_n15_solute(:,1)  = no3_n15_solute(:,1) + noy_deposition(:) / 1.e6_wp / & 
                                                ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) / soil_depth_sl(:,1) * dtime

    !> 4.3 update of mineral P pools
    !!
    IF(.NOT.dsl4jsb_Config(SB_)%flag_sb_prescribe_po4)THEN
       IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d" .OR. TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_1d" ) THEN
            CALL calc_pinorg_dynamics( options, &
                                       dtime, &
                                       TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                                       po4_solute(:,:), &
                                       po4_assoc_fast(:,:), &
                                       po4_assoc_slow(:,:), &
                                       root_fraction_sl(:,:), &
                                       bulk_dens_corr_sl(:,:), &
                                       qmax_po4(:,:), &
                                       p_deposition_sl(:,:), &
                                       rtm_sorption(:,:), &
                                       rmm_sorption(:,:), &
                                       rtm_desorption(:,:), &
                                       rmm_desorption(:,:), &
                                       biochem_mineralisation_po4(:,:), &
                                       mycorrhiza_uptake_po4_sl(:,:), &
                                       plant_uptake_po4_sl(:,:), &
                                       microbial_uptake_po4_sl(:,:), &
                                       net_mineralisation_po4(:,:), &
                                       transport_po4_solute(:,:), &
                                       transport_po4_assoc_fast(:,:), &
                                       transport_po4_assoc_slow(:,:), &
                                       weathering_po4(:,:), &
                                       slow_sorption_po4(:,:), &
                                       occlusion_po4(:,:), &
                                       fast_adsorption_po4(:,:))

       ELSEIF (TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "jsm" ) THEN
            DO isoil = 1,nsoil_sb
              fine_root_carbon(:,isoil) = root_fraction_sl(:,isoil) * veg_pool%fine_root%carbon(:)
            ENDDO
            !fine_root_carbon = root_fraction_sl * 600._wp !temporary solution
            CALL calc_adsorption_po4(options, fine_root_carbon(:,:),  qmax_po4(:,:), p_deposition_sl(:,:), bulk_dens_corr_sl(:,:), &
                                       po4_solute(:,:), po4_assoc_fast(:,:), po4_assoc_slow(:,:), &
                                       weathering_po4(:,:), occlusion_po4(:,:), slow_sorption_po4(:,:), &
                                       transport_po4_solute(:,:), transport_po4_assoc_fast(:,:), transport_po4_assoc_slow(:,:), &
                                       gross_mineralisation_po4(:,:), biochem_mineralisation_po4(:,:), sb_pool%microbial%carbon(:,:),  &
                                       microbial_uptake_po4_sl(:,:), plant_uptake_po4_sl(:,:), fast_adsorption_po4(:,:), &
                                       TRIM(dsl4jsb_Config(SB_)%sb_adsorp_scheme), dsl4jsb_Config(SB_)%flag_sb_CLE_po4, dsl4jsb_Config(SB_)%flag_sb_LK_po4)
            po4_solute(:,:)     = po4_solute(:,:)     + (net_mineralisation_po4(:,:) - plant_uptake_po4_sl(:,:) - mycorrhiza_uptake_po4_sl(:,:) + weathering_po4(:,:) &
                                  + biochem_mineralisation_po4(:,:) + transport_po4_solute(:,:) - fast_adsorption_po4(:,:)) * dtime/1.e6_wp
            po4_assoc_fast(:,:) = po4_assoc_fast(:,:) + (fast_adsorption_po4 - slow_sorption_po4 + transport_po4_assoc_fast) * dtime/1.e6_wp
    
       ENDIF

       po4_assoc_slow(:,:) = po4_assoc_slow(:,:) + (slow_sorption_po4(:,:) - occlusion_po4(:,:) + transport_po4_assoc_slow(:,:)) * dtime/1.e6_wp
       po4_occluded(:,:)   = po4_occluded(:,:)   + (occlusion_po4(:,:) + transport_po4_occluded(:,:)) * dtime/1.e6_wp
       po4_primary(:,:)    = po4_primary(:,:)    + ( - weathering_po4(:,:) + transport_po4_primary(:,:)) * dtime/1.e6_wp

    ENDIF


    !> 4.4 update diagnostics and other state variables depending on organic matter content 
    !!
    bulk_soil_carbon_sl = calc_bulk_soil_carbon(options, sb_pool)
    bulk_dens_corr_sl = calc_bulk_density_correction(bulk_soil_carbon_sl,bulk_dens_sl)
    !DO isoil = 1,nsoil_sb
    !   hlp1(:,isoil) = 0.9_wp * exp(-0.9_wp * SUM(soil_depth_sl(:,1:isoil),DIM=2)) + 0.1_wp
    !ENDDO
    !qmax_org = calc_qmax(bulk_soil_carbon_sl,bulk_dens_sl,qmax_org_om,qmax_org_mineral*hlp1)
    qmax_org = calc_qmax(bulk_soil_carbon_sl,bulk_dens_sl,qmax_org_om,qmax_org_mineral)
    qmax_po4 = calc_qmax(bulk_soil_carbon_sl,bulk_dens_sl,qmax_po4_mineral,qmax_po4_mineral)

    ! NH4
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_nh4)THEN
      IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
        hlp1(:,:) = 0.0_wp
        hlp1(:,1) = nh4_solute_prescribe / soil_depth_sl(:,1) - nh4_solute(:,1)
      ELSE 
        hlp1(:,:) = nh4_solute_prescribe - nh4_solute(:,:)
      END IF
      nh4_solute(:,:)     = nh4_solute(:,:) + hlp1(:,:) 
      nh4_n15_solute(:,:) = nh4_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    ENDIF
    ! NO3
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_no3)THEN
      IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
        hlp1(:,:) = 0.0_wp
        hlp1(:,1) = no3_solute_prescribe / soil_depth_sl(:,1) - no3_solute(:,1)
      ELSE
        hlp1(:,:) = no3_solute_prescribe - no3_solute(:,:)
      END IF
      no3_solute(:,:)     = no3_solute(:,:) + hlp1(:,:) 
      no3_n15_solute(:,:) = no3_n15_solute(:,:) + hlp1(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) 
    ENDIF
    ! PO4
    IF(dsl4jsb_Config(SB_)%flag_sb_prescribe_po4)THEN
      IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == "simple_0d") THEN
        po4_solute(:,:) = 0.0_wp
        po4_solute(:,1) = po4_solute_prescribe / soil_depth_sl(:,1)
      ELSE
        po4_solute(:,:) = po4_solute_prescribe 
      END IF
    ENDIF


    !> 5.0 infrastructure for mass balance checks (see also section 0.9)
    !!
    
    !! pools: calc the difference to what was calculated at the beginning of this routine
    sb_delta_pools%carbon(:,:) = &
       (sb_pool%dom%carbon(:,:) + sb_pool%dom_assoc%carbon(:,:) + &
       sb_pool%soluable_litter%carbon(:,:) + sb_pool%polymeric_litter%carbon(:,:) + &
       sb_pool%woody_litter%carbon(:,:) + sb_pool%fungi%carbon(:,:) + &
       sb_pool%mycorrhiza%carbon(:,:) + sb_pool%microbial%carbon(:,:) + &
       sb_pool%residue%carbon(:,:) + sb_pool%residue_assoc%carbon(:,:)) - & 
       sb_delta_pools%carbon(:,:)
    sb_delta_pools%nitrogen(:,:) = &
       (sb_pool%dom%nitrogen(:,:) + sb_pool%dom_assoc%nitrogen(:,:) + &
       sb_pool%soluable_litter%nitrogen(:,:) + sb_pool%polymeric_litter%nitrogen(:,:) + &
       sb_pool%woody_litter%nitrogen(:,:) + sb_pool%fungi%nitrogen(:,:) + &
       sb_pool%mycorrhiza%nitrogen(:,:) + sb_pool%microbial%nitrogen(:,:) + &
       sb_pool%residue%nitrogen(:,:) + sb_pool%residue_assoc%nitrogen(:,:) + &
       nh4_solute(:,:) + nh4_assoc(:,:) + no3_solute(:,:) + noy(:,:) + n2o(:,:) + n2(:,:)) - &
       sb_delta_pools%nitrogen(:,:)
    sb_delta_pools%phosphorus(:,:) = &
       (sb_pool%dom%phosphorus(:,:) + sb_pool%dom_assoc%phosphorus(:,:) + &
       sb_pool%soluable_litter%phosphorus(:,:) + sb_pool%polymeric_litter%phosphorus(:,:) + &
       sb_pool%woody_litter%phosphorus(:,:) + sb_pool%fungi%phosphorus(:,:) + &
       sb_pool%mycorrhiza%phosphorus(:,:) + sb_pool%microbial%phosphorus(:,:) + &
       sb_pool%residue%phosphorus(:,:) + sb_pool%residue_assoc%phosphorus(:,:) + &
       po4_solute(:,:) + po4_occluded(:,:) + po4_assoc_slow(:,:) + po4_assoc_fast(:,:) + po4_primary(:,:)) - &
       sb_delta_pools%phosphorus(:,:)
    sb_delta_pools%carbon13(:,:) = &
       (sb_pool%dom%carbon13(:,:) + sb_pool%dom_assoc%carbon13(:,:) + &
       sb_pool%soluable_litter%carbon13(:,:) + sb_pool%polymeric_litter%carbon13(:,:) + &
       sb_pool%woody_litter%carbon13(:,:) + sb_pool%fungi%carbon13(:,:) + &
       sb_pool%mycorrhiza%carbon13(:,:) + sb_pool%microbial%carbon13(:,:) + &
       sb_pool%residue%carbon13(:,:) + sb_pool%residue_assoc%carbon13(:,:)) - &
       sb_delta_pools%carbon13(:,:)
    sb_delta_pools%carbon14(:,:) = &
       (sb_pool%dom%carbon14(:,:) + sb_pool%dom_assoc%carbon14(:,:) + &
       sb_pool%soluable_litter%carbon14(:,:) + sb_pool%polymeric_litter%carbon14(:,:) + &
       sb_pool%woody_litter%carbon14(:,:) + sb_pool%fungi%carbon14(:,:) + &
       sb_pool%mycorrhiza%carbon14(:,:) + sb_pool%microbial%carbon14(:,:) + &
       sb_pool%residue%carbon14(:,:) + sb_pool%residue_assoc%carbon14(:,:)) - &
       sb_delta_pools%carbon14(:,:)
    sb_delta_pools%nitrogen15(:,:) = &
       (sb_pool%dom%nitrogen15(:,:) + sb_pool%dom_assoc%nitrogen15(:,:) + &
       sb_pool%soluable_litter%nitrogen15(:,:) + sb_pool%polymeric_litter%nitrogen15(:,:) + &
       sb_pool%woody_litter%nitrogen15(:,:) + sb_pool%fungi%nitrogen15(:,:) + &
       sb_pool%mycorrhiza%nitrogen15(:,:) + sb_pool%microbial%nitrogen15(:,:) + &
       sb_pool%residue%nitrogen15(:,:) + sb_pool%residue_assoc%nitrogen15(:,:) + &
       nh4_n15_solute(:,:) + nh4_n15_assoc(:,:) + no3_n15_solute(:,:) + noy_n15(:,:) + n2o_n15(:,:) + n2_n15(:,:)) - &
       sb_delta_pools%nitrogen15(:,:)

    !! fluxes
    total_flux_carbon(:) =                  veg_litterfall%leaf%carbon        + &
      veg_litterfall%fine_root%carbon + veg_litterfall%coarse_root%carbon + &
      veg_litterfall%sap_wood%carbon  + veg_litterfall%heart_wood%carbon  + &
      veg_litterfall%labile%carbon    + veg_litterfall%reserve%carbon     + &
      veg_litterfall%fruit%carbon     + veg_litterfall%seed_bed%carbon    + &
      veg_exudation%carbon                                                    - &
      SUM(sb_mycorrhiza_export%carbon(:,:) * soil_depth_sl(:,:), DIM=2)  - &
      SUM(het_respiration(:,:) * soil_depth_sl(:,:), DIM=2) * dtime/1.e6_wp                          ! DIM=2: sum only over dimension 2, i.e., soil layers
    total_flux_nitrogen(:) =                  veg_litterfall%leaf%nitrogen    + &
      veg_litterfall%fine_root%nitrogen + veg_litterfall%coarse_root%nitrogen + &
      veg_litterfall%sap_wood%nitrogen  + veg_litterfall%heart_wood%nitrogen  + &
      veg_litterfall%labile%nitrogen    + veg_litterfall%reserve%nitrogen     + &
      veg_litterfall%fruit%nitrogen     + veg_litterfall%seed_bed%nitrogen    + &
      SUM(asymb_n_fixation(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp   + &
      (nhx_deposition(:) + noy_deposition(:)) * dtime/1.e6_wp                 + & 
      veg_exudation%nitrogen                                                  - &
      SUM(sb_mycorrhiza_export%nitrogen(:,:) * soil_depth_sl(:,:), DIM=2)     - &
      (leaching_nh4_solute(:) + leaching_no3_solute(:)                        + &
       emission_nh3(:) + emission_noy(:) + emission_n2o(:) + emission_n2(:)   + &
       SUM(plant_uptake_nh4_sl(:,:) * soil_depth_sl(:,:), DIM=2)              + &
       SUM(plant_uptake_no3_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp                       
    total_flux_phosphorus(:) =                  veg_litterfall%leaf%phosphorus    + &
      veg_litterfall%fine_root%phosphorus + veg_litterfall%coarse_root%phosphorus + &
      veg_litterfall%sap_wood%phosphorus  + veg_litterfall%heart_wood%phosphorus  + &
      veg_litterfall%labile%phosphorus    + veg_litterfall%reserve%phosphorus     + &
      veg_litterfall%fruit%phosphorus     + veg_litterfall%seed_bed%phosphorus    + &
      p_deposition(:) * dtime/1.e6_wp                                             + &
      veg_exudation%phosphorus                                                    - &
      SUM(sb_mycorrhiza_export%phosphorus(:,:) * soil_depth_sl(:,:), DIM=2)       - &
      (leaching_po4_solute(:)                                                     + &
       SUM(plant_uptake_po4_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp       
    total_flux_carbon13(:) =                  veg_litterfall%leaf%carbon13        + &
      veg_litterfall%fine_root%carbon13 + veg_litterfall%coarse_root%carbon13 + &
      veg_litterfall%sap_wood%carbon13  + veg_litterfall%heart_wood%carbon13  + &
      veg_litterfall%labile%carbon13    + veg_litterfall%reserve%carbon13     + &
      veg_litterfall%fruit%carbon13     + veg_litterfall%seed_bed%carbon13    + &
      veg_exudation%carbon13                                                      - &
      SUM(sb_mycorrhiza_export%carbon13(:,:) * soil_depth_sl(:,:), DIM=2)    - &
      SUM(het_respiration_c13(:,:) * soil_depth_sl(:,:), DIM=2) * dtime/1.e6_wp
    total_flux_carbon14(:) =                  veg_litterfall%leaf%carbon14        + &
      veg_litterfall%fine_root%carbon14 + veg_litterfall%coarse_root%carbon14 + &
      veg_litterfall%sap_wood%carbon14  + veg_litterfall%heart_wood%carbon14  + &
      veg_litterfall%labile%carbon14    + veg_litterfall%reserve%carbon14     + &
      veg_litterfall%fruit%carbon14     + veg_litterfall%seed_bed%carbon14    + &
      veg_exudation%carbon14                                                      - &
      SUM(sb_mycorrhiza_export%carbon14(:,:) * soil_depth_sl(:,:), DIM=2)    - &
      SUM(het_respiration_c14(:,:) * soil_depth_sl(:,:), DIM=2) * dtime/1.e6_wp
    total_flux_nitrogen15(:) =                  veg_litterfall%leaf%nitrogen15      + &
      veg_litterfall%fine_root%nitrogen15 + veg_litterfall%coarse_root%nitrogen15   + &
      veg_litterfall%sap_wood%nitrogen15  + veg_litterfall%heart_wood%nitrogen15    + &
      veg_litterfall%labile%nitrogen15    + veg_litterfall%reserve%nitrogen15       + &
      veg_litterfall%fruit%nitrogen15     + veg_litterfall%seed_bed%nitrogen15      + &
      SUM(asymb_n_fixation_n15(:,:) * soil_depth_sl(:,:), DIM=2) *dtime/1.e6_wp     + &
      (nhx_deposition(:) + noy_deposition(:)) /                                       & 
         ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(0.0_wp)) *dtime/1.e6_wp         + &
      veg_exudation%nitrogen15                                                      - &
      SUM(sb_mycorrhiza_export%nitrogen15(:,:) * soil_depth_sl(:,:), DIM=2)         - &
      (leaching_nh4_n15_solute(:) + leaching_no3_n15_solute(:)                      + &
       emission_nh3_n15(:) + emission_noy_n15(:) + emission_n2o_n15(:)              + & 
       emission_n2_n15(:)                                                           + &
       SUM(plant_uptake_nh4_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)                + &
       SUM(plant_uptake_no3_n15_sl(:,:) * soil_depth_sl(:,:), DIM=2)) * dtime/1.e6_wp  


! Pass back values to process memory
#define Passback_veg_litterfall
#define Passback_veg_exudation
#define Passback_sb_pool
#define Passback_sb_formation
#define Passback_sb_loss
#define Passback_sb_transport
#define Passback_sb_delta_pools
#define Passback_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Passback_veg_litterfall
#undef Passback_veg_exudation
#undef Passback_sb_pool
#undef Passback_sb_formation
#undef Passback_sb_loss
#undef Passback_sb_transport
#undef Passback_sb_delta_pools
#undef Passback_sb_mycorrhiza_export

  END SUBROUTINE update_sb_pools



  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_pools
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to calculate fast_adsorption_po4 based on Langmuir equilibrium(Yang et al. 2014)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_pinorg_dynamics(options, &
                                  dtime, &
                                  sb_model_scheme, &
                                  po4_solute, &
                                  po4_assoc_fast, &
                                  po4_assoc_slow, &
                                  root_fraction_sl, &
                                  bulk_dens_corr_sl, &
                                  qmax_po4_corr_sl, &
                                  p_deposition_sl, &
                                  rtm_sorption, &
                                  rmm_sorption, &
                                  rtm_desorption, &
                                  rmm_desorption, &
                                  biochem_mineralisation_po4, &
                                  mycorrhiza_uptake_po4_sl, &
                                  plant_uptake_po4_sl, &
                                  microbial_uptake_po4_sl, &
                                  net_mineralisation_po4, &
                                  transport_po4_solute, &
                                  transport_po4_assoc_fast, &
                                  transport_po4_assoc_slow, &
                                  weathering_po4, &
                                  occlusion_po4, &
                                  slow_sorption_po4, &
                                  fast_adsorption_po4)

    USE mo_sb_constants

    IMPLICIT NONE
    ! 0.1 inout
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options                       !< dimensions
    REAL(wp),                                         INTENT(in)    :: dtime                         !< length of one time step in seconds
    CHARACTER(len=*),                                 INTENT(in)    :: sb_model_scheme               !< from sb config
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: po4_assoc_slow, &             !< slow minerally associated PO4 pool [mol/m3]
                                                                       root_fraction_sl, &           !< fraction of root mass in layer
                                                                       bulk_dens_corr_sl, &          !< bulk density of the mineral soil [kg / m3]
                                                                       qmax_po4_corr_sl, &           !< maximum sorption capacity of po4, [mol P/ m3]
                                                                       p_deposition_sl, &            !< actual P deposition rate, [µmol P / m3 / s]
                                                                       rtm_sorption, &               !< temperature response scalar for po4 sorption
                                                                       rtm_desorption, &             !< temperature response scalar for po4 desorption
                                                                       rmm_sorption, &               !< moisture response scalar for po4 sorption
                                                                       rmm_desorption, &             !< moisture response scalar for po4 desorption
                                                                       biochem_mineralisation_po4, & !< biomineralisation rate [mol / m3 / timestep]
                                                                       transport_po4_solute, &       !< transport of soluble PO4 [mol / m3 / timestep]
                                                                       transport_po4_assoc_fast, &   !< transport of fast_assoc_po4 [mol / m3 / timestep]
                                                                       transport_po4_assoc_slow, &   !< transport of slow_assoc_po4 [mol / m3 / timestep]
                                                                       weathering_po4                !< phosphorus weathering rate [mol P / m3 / timestep]   
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: po4_solute, &                 !< PO4 in solution [mol / m3]
                                                                       po4_assoc_fast, &             !< fast minerally associated PO4 pool [mol/m3]
                                                                       mycorrhiza_uptake_po4_sl, &   !< mycorrhizal PO4 uptake [mol / m3 / timestep]
                                                                       plant_uptake_po4_sl, &        !< plant PO4 uptake [mol / m3 / timestep]
                                                                       microbial_uptake_po4_sl, &    !< microbial PO4 uptake [mol / m3 / timestep]
                                                                       occlusion_po4, &              !< phosphorus occlusion rate [mol P / m3 / timestep]
                                                                       net_mineralisation_po4, &     !< net mineralisation of PO4 [mol / m3 / timestep]
                                                                       slow_sorption_po4             !< phosphorus desorption rate, assoc_fast --> assoc_slow [mol P / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out)   :: fast_adsorption_po4           !< phosphorus adsorption rate, solute --> assoc_fast [mol P / m3 / timestep]


    ! 0.2 local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1, hlp2, hlp3, hlp4
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: km_adsorpt_po4_act, qmax_po4_mineral_act, k_desorpt_po4_act, k_adsorpt_po4_act, &
                                                        rtm_sorption_act, rmm_sorption_act, rtm_desorption_act, rmm_desorption_act, bulk_dens_corr_sl_act 
    INTEGER                                          :: isoil                !< loop over options%nsoil_sb                        
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_pinorg_dynamics'

    ! 0.3 rate responsers for 1d and 0d soil model
    IF (sb_model_scheme=="simple_0d") THEN
      rtm_sorption_act(:,:) = 0.0_wp 
      rmm_sorption_act(:,:) = 0.0_wp
      rtm_desorption_act(:,:) = 0.0_wp 
      rmm_desorption_act(:,:) = 0.0_wp
      k_adsorpt_po4_act(:,:) = 0.0_wp 
      k_desorpt_po4_act(:,:) = 0.0_wp
      km_adsorpt_po4_act(:,:) = 0.0_wp
      qmax_po4_mineral_act(:,:) = 0.0_wp
      bulk_dens_corr_sl_act(:,:) = 0.0_wp
      DO isoil = 1,options%nsoil_sb
         rtm_sorption_act(:,1) = rtm_sorption_act(:,1) + & 
                                         rtm_sorption(:,isoil) * root_fraction_sl(:,isoil) 
         rmm_sorption_act(:,1) = rmm_sorption_act(:,1) + & 
                                         rmm_sorption(:,isoil) * root_fraction_sl(:,isoil) 
         rtm_desorption_act(:,1) = rtm_desorption_act(:,1) + & 
                                         rtm_desorption(:,isoil) * root_fraction_sl(:,isoil) 
         rmm_desorption_act(:,1) = rmm_desorption_act(:,1) + & 
                                         rmm_desorption(:,isoil) * root_fraction_sl(:,isoil) 
         bulk_dens_corr_sl_act(:,1) = bulk_dens_corr_sl_act(:,1) + & 
                                         bulk_dens_corr_sl(:,isoil) * root_fraction_sl(:,isoil) 
         qmax_po4_mineral_act(:,1) = qmax_po4_mineral_act(:,1) + &
                                         qmax_po4_corr_sl(:,isoil) * root_fraction_sl(:,isoil) 
      ENDDO
      k_adsorpt_po4_act(:,1) = rtm_sorption_act(:,1) * rmm_sorption_act(:,1) * k_adsorpt_po4
      k_desorpt_po4_act(:,1) = rtm_desorption_act(:,1) * rmm_desorption_act(:,1) * k_desorpt_po4
      km_adsorpt_po4_act(:,1) = bulk_dens_corr_sl_act(:,1) * km_adsorpt_po4
    ELSEIF (sb_model_scheme=="simple_1d") THEN
      k_adsorpt_po4_act(:,:) = rtm_sorption(:,:) * rmm_sorption(:,:) * k_adsorpt_po4
      k_desorpt_po4_act(:,:) = rtm_desorption(:,:) * rmm_desorption(:,:) * k_desorpt_po4
      km_adsorpt_po4_act(:,:) = bulk_dens_corr_sl(:,:) * km_adsorpt_po4
      qmax_po4_mineral_act(:,:) = qmax_po4_corr_sl(:,:)
    ENDIF

    !> 1.0 limit to possible assoc_slow, 0d model
    !!
    IF (sb_model_scheme=="simple_0d") THEN
      hlp1(:,1) = k_desorpt_po4_act(:,1) * po4_assoc_slow(:,1) * dtime
      hlp2(:,1) = k_adsorpt_po4_act(:,1) * po4_assoc_fast(:,1) * dtime
      hlp3(:,1) = po4_assoc_slow(:,1) + hlp2(:,1) - hlp1(:,1)  - occlusion_po4(:,1) * dtime / 1000000._wp
      WHERE (hlp3(:,1) < eps8)
        hlp3(:,1) = po4_assoc_slow(:,1) + hlp2(:,1)
        hlp4(:,1) = hlp1(:,1) + occlusion_po4(:,1) * dtime / 1000000._wp
        hlp1(:,1) = hlp1(:,1) * hlp3(:,1) / hlp4(:,1)
        occlusion_po4(:,1) = occlusion_po4(:,1) * hlp3(:,1) / hlp4(:,1)
        slow_sorption_po4(:,1) = (hlp2(:,1) - hlp1(:,1)) / dtime * 1000000._wp
      ENDWHERE
    !> 1.1 Langmuir equilibrium between soluble P and labile P, 0d model
    !!
      fast_adsorption_po4(:,:) = 0.0_wp
      hlp1(:,:) = 0.0_wp; hlp2(:,:) = 0.0_wp; hlp4(:,:) = 0.0_wp; hlp3(:,:) = 0.0_wp; 
      hlp1(:,1) = (po4_solute(:,1) + km_adsorpt_po4_act(:,1))**2.0_wp / &
                  ((po4_solute(:,1) + km_adsorpt_po4_act(:,1))**2.0_wp + qmax_po4_mineral_act(:,1) * km_adsorpt_po4_act(:,1))
      hlp2(:,1) = (weathering_po4(:,1) + biochem_mineralisation_po4(:,1) + net_mineralisation_po4(:,1) + transport_po4_solute(:,1) + p_deposition_sl(:,1) - &
            (plant_uptake_po4_sl(:,1) + mycorrhiza_uptake_po4_sl(:,1) + slow_sorption_po4(:,1))) * dtime / 1000000._wp

      WHERE (po4_solute(:,1) < (-hlp2(:,1)) .OR. po4_assoc_fast(:,1) < (-hlp2(:,1)) )
        hlp1(:,1) = po4_solute(:,1)/(po4_solute(:,1)+po4_assoc_fast(:,1))
      ENDWHERE

      fast_adsorption_po4(:,1) = hlp2(:,1) + slow_sorption_po4(:,1) * dtime / 1000000._wp - hlp1(:,1) * hlp2(:,1)
      fast_adsorption_po4(:,1) = fast_adsorption_po4(:,1) / dtime * 1000000._wp 

    !> 2.0 limit to possible assoc_slow, 1d model
    !!
    ELSEIF (sb_model_scheme=="simple_1d") THEN
      hlp1(:,:) = k_desorpt_po4_act(:,:) * po4_assoc_slow(:,:) * dtime
      hlp2(:,:) = k_adsorpt_po4_act(:,:) * po4_assoc_fast(:,:) * dtime
      hlp3(:,:) = po4_assoc_slow(:,:) + hlp2(:,:) - hlp1(:,:) + (transport_po4_assoc_slow(:,:) - occlusion_po4(:,:)) * dtime / 1000000._wp
      WHERE (hlp3(:,:) < eps8)
        hlp3(:,:) = po4_assoc_slow(:,:) + transport_po4_assoc_slow(:,:) * dtime / 1000000._wp + hlp2(:,:)
        hlp4(:,:) = hlp1(:,:) + occlusion_po4(:,:) * dtime / 1000000._wp
        hlp1(:,:) = hlp1(:,:) * hlp3(:,:) / hlp4(:,:)
        occlusion_po4(:,:) = occlusion_po4(:,:) * hlp3(:,:) / hlp4(:,:)
        slow_sorption_po4(:,:) = (hlp2(:,:) - hlp1(:,:)) / dtime * 1000000._wp
      ENDWHERE
    !> 2.1 Langmuir equilibrium between soluble P and labile P, 1d model
    !!
      hlp1(:,:) = (po4_solute(:,:) + km_adsorpt_po4_act(:,:))**2.0_wp / &
                  ((po4_solute(:,:) + km_adsorpt_po4_act(:,:))**2.0_wp + qmax_po4_mineral_act(:,:) * km_adsorpt_po4_act(:,:))
      hlp2(:,:) = (weathering_po4(:,:) + biochem_mineralisation_po4(:,:) + net_mineralisation_po4(:,:) + transport_po4_solute(:,:) + p_deposition_sl(:,:) - &
            (plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) + slow_sorption_po4(:,:)) + transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp
      WHERE (po4_solute(:,:) < (-hlp2(:,:)) .OR. po4_assoc_fast(:,:) < (-hlp2(:,:)) )
        hlp1(:,:) = po4_solute(:,:)/(po4_solute(:,:)+po4_assoc_fast(:,:))
      ENDWHERE

      fast_adsorption_po4(:,:) = hlp2(:,:) + (slow_sorption_po4(:,:) - transport_po4_assoc_fast(:,:)) * dtime / 1000000._wp - hlp1(:,:) * hlp2(:,:)
      fast_adsorption_po4(:,:) = fast_adsorption_po4(:,:) / dtime * 1000000._wp 
    ENDIF

    !> 3.0 update pools, both 1d and 0d
    !!
    po4_solute(:,:) = po4_solute(:,:) + hlp1(:,:)*hlp2(:,:)
    po4_assoc_fast(:,:) = po4_assoc_fast(:,:) + (1 - hlp1(:,:))*hlp2(:,:)


  END SUBROUTINE calc_pinorg_dynamics

END MODULE mo_sb_update_pools
