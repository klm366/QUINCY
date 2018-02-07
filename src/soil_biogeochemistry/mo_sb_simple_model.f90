!>
!! @brief Contains the routines for the Task sb_simple_model 
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!! Second version by Sönke Zaehle (20170331): added a 4 pool model 
!! 3rd version: Jan Engel (20170627): split calc_sb_fluxes_simple() into multiple calc_ routines and rename to update_sb_simple_model()
!!
MODULE mo_sb_simple_model

#include "dsl_jsbach.inc"

  USE mo_kind,                  ONLY: wp, eps8, eps4
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_jsb_class,             ONLY: Get_model
  USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
  USE mo_jsb_model_class,       ONLY: t_jsb_model
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
  USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
  USE mo_base_constants,        ONLY: dtime, one_day, one_year

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_sb_simple_model 


  CHARACTER(len=*), PARAMETER :: modname = 'mo_sb_simple_model'

CONTAINS


  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> update soil biogeochemical pools 
  !! 
  !! routine to update soil litter and mineral pools (simple version for vegetation model testing)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_simple_model(tile, options)

    USE mo_veg_constants,         ONLY: itree
    USE mo_sb_jsm_transport,      ONLY: calc_bioturbation_rate

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_config(SOIL_)

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4slm_Use_memory_a_veg        !  USE mo_veg_memory_class,      ONLY: t_veg_element_a, t_veg_pool_a
    USE mo_veg_memory_class,      ONLY: veg_pool, veg_litterfall
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
    dsl4jsb_Use_memory(A2L_)    
    dsl4slm_Use_memory_a_sb
    USE mo_sb_memory_class,       ONLY: sb_pool, sb_formation, sb_loss, sb_transport, sb_mycorrhiza_export

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET    :: tile          !< one tile
    TYPE(t_jsb_task_options), INTENT(in)               :: options       !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    INTEGER                               :: isoil                !< loop over options%nsoil_sb
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) ::      &
                                             apparent_mic_cn      , &
                                             root_fraction_act_sl , &      ! and below needed for switch between 0d and 1d
                                             rtm_decomposition_act, &
                                             rmm_decomposition_act, & 
                                             km1_uptake_nh4_act, &
                                             km2_uptake_nh4_act, &
                                             km1_uptake_no3_act, &
                                             km2_uptake_no3_act, &
                                             km1_uptake_po4_act, &
                                             km2_uptake_po4_act, &
                                             km_denitrification_c_act, &
                                             km_denitrification_no3_act, &
                                             rtm_nitrification_act, &
                                             rtm_denitrification_act, &
                                             rmm_nitrification_act, &
                                             rtm_gasdiffusion_act, &
                                             rmm_gasdiffusion_act, &
                                             p_deposition_act, &
                                             k_weath_mineral_act, &
                                             k_adsorpt_po4_act, &
                                             k_desorpt_po4_act, &
                                             km_dip_biochem_po4_act, &
                                             qmax_po4_mineral_act, &
                                             anaerobic_volume_fraction_act
    INTEGER                               :: nsoil_sb
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_simple_model'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_config(SOIL_)
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    dsl4jsb_Def_memory(A2L_)
    ! Declare pointers to variables in memory
    ! VEG_
    dsl4jsb_Real2D_onChunk :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk :: unit_uptake_p_pot
    ! Declare pointers to variables in memory
    ! A2L
    dsl4jsb_Real2D_onChunk :: p_deposition
    ! SOIL_ 2D
    dsl4jsb_Real2D_onChunk :: drainage_fraction 
    ! SOIL_ 3D
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: bulk_dens_sl 
    dsl4jsb_Real3D_onChunk :: root_fraction_sl
    dsl4jsb_Real3D_onChunk :: percolation_sl
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
    dsl4jsb_Real2D_onChunk :: emission_noy
    dsl4jsb_Real2D_onChunk :: emission_n2o
    dsl4jsb_Real2D_onChunk :: emission_n2
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: rtm_decomposition 
    dsl4jsb_Real3D_onChunk :: rtm_hsc
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rtm_sorption
    dsl4jsb_Real3D_onChunk :: rtm_desorption
    dsl4jsb_Real3D_onChunk :: rmm_decomposition 
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl 
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nitrification_no3 
    dsl4jsb_Real3D_onChunk :: nitrification_noy 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_noy 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o 
    dsl4jsb_Real3D_onChunk :: denitrification_n2 
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: microbial_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_norg_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: weathering_po4
    dsl4jsb_Real3D_onChunk :: occlusion_po4
    dsl4jsb_Real3D_onChunk :: slow_sorption_po4
    dsl4jsb_Real3D_onChunk :: biochem_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: transport_po4_assoc_slow 
    dsl4jsb_Real3D_onChunk :: transport_po4_occluded
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2 
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2 
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: po4_assoc_fast 
    dsl4jsb_Real3D_onChunk :: po4_assoc_slow
    dsl4jsb_Real3D_onChunk :: po4_occluded
    dsl4jsb_Real3D_onChunk :: po4_primary
    dsl4jsb_Real3D_onChunk :: k_bioturb 
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation 
    dsl4jsb_Real3D_onChunk :: qmax_po4
    ! SB_ 2D
    dsl4jsb_Real2D_onChunk :: leaching_nh4_n15_solute
    dsl4jsb_Real2D_onChunk :: leaching_no3_n15_solute
    dsl4jsb_Real2D_onChunk :: emission_noy_n15
    dsl4jsb_Real2D_onChunk :: emission_n2o_n15
    dsl4jsb_Real2D_onChunk :: emission_n2_n15
    ! SB_ 3D
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: noy_n15
    dsl4jsb_Real3D_onChunk :: n2o_n15
    dsl4jsb_Real3D_onChunk :: n2_n15 
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: transport_nh4_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_n15_solute
    dsl4jsb_Real3D_onChunk :: transport_noy_n15
    dsl4jsb_Real3D_onChunk :: transport_n2o_n15
    dsl4jsb_Real3D_onChunk :: transport_n2_n15 
    dsl4jsb_Real3D_onChunk :: volatilisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: nitrification_no3_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o_n15
    dsl4jsb_Real3D_onChunk :: denitrification_noy_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o_n15 
    dsl4jsb_Real3D_onChunk :: denitrification_n2_n15 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation_n15 
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
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_config(SOIL_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(A2L_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(A2L_,      p_deposition)                ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_n_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_p_pot)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_,     drainage_fraction)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     bulk_dens_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_n15_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_n15_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)          
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy_n15)          
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o_n15)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2_n15)         
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       anaerobic_volume_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       microbial_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_norg_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       weathering_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       occlusion_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       slow_sorption_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       biochem_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_fast)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_assoc_slow)    ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_occluded)      ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_fast)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_assoc_slow)              ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_occluded)                ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_primary)                 ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       qmax_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       noy_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2_n15) 
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation_n15)            ! out
    ! Get process memory pool structures or individual variables from pools
    dsl4slm_get_vegpool_var(pool,fine_root,carbon)
#define Get_veg_litterfall
#define Get_sb_pool
#define Get_sb_formation
#define Get_sb_loss
#define Get_sb_transport
#define Get_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Get_veg_litterfall
#undef Get_sb_pool
#undef Get_sb_formation
#undef Get_sb_loss
#undef Get_sb_transport
#undef Get_sb_mycorrhiza_export



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    microbial_uptake_no3_sl = 0.0_wp

    !> 1.0 Set values of soil profile variables to allow switching between simple_0d and simple1d in the 
    !!     following routines 
    !! 
    CALL calc_vertical_profile(options, &
                               lctlib, &
                               TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                               soil_depth_sl(:,:), &
                               bulk_dens_sl(:,:), &
                               nh4_solute(:,:), &
                               root_fraction_sl(:,:), &
                               rtm_decomposition(:,:), &
                               rmm_decomposition(:,:), &
                               rtm_hsc(:,:), &
                               rmm_hsc(:,:), &
                               rtm_nitrification(:,:), &
                               rmm_nitrification(:,:), &
                               rtm_denitrification(:,:), &
                               rtm_gasdiffusion(:,:), &
                               rmm_gasdiffusion(:,:), &
                               rtm_sorption(:,:), &
                               rtm_desorption(:,:), &
                               rmm_sorption(:,:), &
                               rmm_desorption(:,:), &
                               bulk_dens_corr_sl(:,:), &
                               p_deposition(:), &
                               anaerobic_volume_fraction_sl(:,:), &
                               root_fraction_act_sl(:,:), &
                               km1_uptake_nh4_act(:,:), &
                               km2_uptake_nh4_act(:,:), &
                               km1_uptake_no3_act(:,:), &
                               km2_uptake_no3_act(:,:), &
                               km1_uptake_po4_act(:,:), &
                               km2_uptake_po4_act(:,:), &
                               km_denitrification_c_act(:,:), &
                               km_denitrification_no3_act(:,:), &
                               apparent_mic_cn(:,:), &
                               rtm_decomposition_act(:,:), &
                               rmm_decomposition_act(:,:), &
                               rtm_nitrification_act(:,:), &
                               rmm_nitrification_act(:,:), &
                               rtm_denitrification_act(:,:), &
                               rtm_gasdiffusion_act(:,:), &
                               rmm_gasdiffusion_act(:,:), &
                               p_deposition_act(:,:), &
                               k_weath_mineral_act(:,:), &
                               k_adsorpt_po4_act(:,:), &
                               k_desorpt_po4_act(:,:), &
                               km_dip_biochem_po4_act(:,:), &
                               anaerobic_volume_fraction_act(:,:) )

    !> 2.0 Calculate potential rates of uptake and losses
    !!
    !> 2.1 Potential plant nutrient uptake rates
    !!
    CALL calc_potential_plant_uptake_rate(options, &
                                          lctlib, &
                                          soil_depth_sl(:,:), &
                                          veg_pool%fine_root%carbon(:), &
                                          root_fraction_act_sl(:,:), &
                                          rtm_decomposition_act(:,:), &
                                          km1_uptake_nh4_act(:,:), &
                                          km2_uptake_nh4_act(:,:), &
                                          km1_uptake_no3_act(:,:), &
                                          km2_uptake_no3_act(:,:), &
                                          km1_uptake_po4_act(:,:), &
                                          km2_uptake_po4_act(:,:), &
                                          nh4_solute(:,:), &
                                          no3_solute(:,:), &
                                          po4_solute(:,:), &
                                          plant_uptake_nh4_sl(:,:), &
                                          plant_uptake_no3_sl(:,:), &
                                          plant_uptake_po4_sl(:,:), &
                                          unit_uptake_n_pot(:), &
                                          unit_uptake_p_pot(:) )

    !> 2.2 Potential plant nutrient uptake rates
    !!
    CALL calc_potential_mycorrhiza_uptake_rate(options, &
                                               lctlib, &
                                               dsl4jsb_Config(SB_)%flag_myco_inorganic, &
                                               soil_depth_sl(:,:), &
                                               sb_pool%mycorrhiza%carbon(:,:), &
                                               sb_pool%mycorrhiza%nitrogen(:,:), &
                                               rtm_decomposition_act(:,:), &
                                               km1_uptake_nh4_act(:,:), &
                                               km2_uptake_nh4_act(:,:), &
                                               km1_uptake_no3_act(:,:), &
                                               km2_uptake_no3_act(:,:), &
                                               km1_uptake_po4_act(:,:), &
                                               km2_uptake_po4_act(:,:), &
                                               nh4_solute(:,:), &
                                               no3_solute(:,:), &
                                               po4_solute(:,:), &
                                               mycorrhiza_uptake_nh4_sl(:,:), &
                                               mycorrhiza_uptake_no3_sl(:,:), &
                                               mycorrhiza_uptake_po4_sl(:,:) ) 

    !> 2.4 Calculate the potential decomposition rate of soil organic matter and litter pools 
    !!
    CALL calc_potential_decomposition_rate(options, &  
                                           dsl4jsb_Config(SB_)%flag_myco_priming, &                  
                                           rtm_decomposition_act(:,:), &
                                           rmm_decomposition_act(:,:), &
                                           sb_pool, &
                                           sb_loss )


    !> 2.5 Loss pathways for the dynamic loss calculation scheme 
    !!
    IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme)=='dynamic')THEN
       CALL calc_potential_nitrification_denitrification_rate(options, &
                                                               rtm_nitrification_act(:,:), &
                                                               rmm_nitrification_act(:,:), &
                                                               rtm_denitrification_act(:,:), &
                                                               km_denitrification_c_act(:,:), &
                                                               km_denitrification_no3_act(:,:), &
                                                               anaerobic_volume_fraction_act(:,:), &
                                                               nh4_solute(:,:), &
                                                               no3_solute(:,:), &
                                                               sb_pool%microbial%carbon(:,:), &
                                                               nitrification_no3(:,:), &
                                                               denitrification_n2(:,:))

       CALL calc_inorganic_transport_rate(options, &
                                          TRIM(dsl4jsb_Config(SOIL_)%soil_moisture_scheme), &
                                          TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                                          soil_depth_sl(:,:), &
                                          drainage_fraction(:), &
                                          percolation_sl(:,:), &
                                          rtm_gasdiffusion_act(:,:), &
                                          rmm_gasdiffusion_act(:,:), &
                                          nh4_solute(:,:), &
                                          nh4_n15_solute(:,:), &
                                          no3_solute(:,:), &
                                          no3_n15_solute(:,:), &
                                          po4_solute(:,:), &
                                          noy(:,:), &
                                          noy_n15(:,:), &
                                          n2o(:,:), &
                                          n2o_n15(:,:), &
                                          n2(:,:), &
                                          n2_n15(:,:), &
                                          transport_nh4_solute(:,:), &
                                          transport_nh4_n15_solute(:,:), &
                                          transport_no3_solute(:,:), &
                                          transport_no3_n15_solute(:,:), &
                                          transport_po4_solute(:,:), &
                                          transport_noy(:,:), &     
                                          transport_noy_n15(:,:), &   
                                          transport_n2o(:,:), &     
                                          transport_n2o_n15(:,:), & 
                                          transport_n2(:,:), &      
                                          transport_n2_n15(:,:), & 
                                          leaching_nh4_solute(:), &
                                          leaching_nh4_n15_solute(:), &
                                          leaching_no3_solute(:), &
                                          leaching_no3_n15_solute(:), &
                                          leaching_po4_solute(:), &
                                          emission_noy(:), &    
                                          emission_noy_n15(:), &
                                          emission_n2o(:), &    
                                          emission_n2o_n15(:), & 
                                          emission_n2(:), &  
                                          emission_n2_n15(:))

    ENDIF  ! IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme)=='dynamic')

    !> 2.6 transport due to bioturbation
    !!
    IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme)=='simple_1d')THEN
       k_bioturb(:,:) = calc_bioturbation_rate(options, bulk_dens_corr_sl(:,:))
       CALL calc_bioturbation_transport_wrapper(options, k_bioturb(:,:), soil_depth_sl(:,:), &
                                                sb_pool%microbial, sb_pool%residue, & 
                                                sb_transport%microbial, sb_transport%residue, &
                                                po4_assoc_fast, po4_assoc_slow, po4_occluded, &
                                                transport_po4_assoc_fast, transport_po4_assoc_slow, & 
                                                transport_po4_occluded)
    ENDIF

    !> 2.7 Calculate the potential PO4 flux rate: biomineralisation, weathering, slow sorption, occlusion 
    !!
    CALL calc_potential_p_flux_rate(options, &                  
                                    rtm_decomposition_act(:,:), &
                                    rmm_decomposition_act(:,:), &
                                    km_dip_biochem_po4_act(:,:), &
                                    po4_solute(:,:), &
                                    po4_assoc_fast(:,:), &
                                    po4_assoc_slow(:,:), &
                                    po4_occluded(:,:), &
                                    po4_primary(:,:), &
                                    k_weath_mineral_act(:,:), &
                                    k_adsorpt_po4_act(:,:), &
                                    k_desorpt_po4_act(:,:), &
                                    sb_pool%residue%phosphorus(:,:), & 
                                    weathering_po4(:,:), &
                                    occlusion_po4(:,:), &
                                    slow_sorption_po4(:,:), &
                                    biochem_mineralisation_po4(:,:))

    !> 3.0 Actual fluxes of in and out of SOM 
    !!
 
    !> 3.1 Calculation of actual plant, mycorrhizal, microbial uptake rate as well as 
    !!     nitrification, denitrification and biomineralisation rates
    !!
    CALL calc_actual_flux_rates(options, &
                                nh4_solute(:,:), &
                                no3_solute(:,:), &
                                po4_solute(:,:), &
                                po4_assoc_fast(:,:), &
                                transport_nh4_solute(:,:), &
                                transport_no3_solute(:,:), &
                                transport_po4_solute(:,:), &
                                transport_po4_assoc_fast(:,:), &
                                weathering_po4(:,:), &
                                p_deposition_act(:,:), &
                                slow_sorption_po4(:,:), &
                                apparent_mic_cn(:,:), &
                                km2_uptake_nh4_act(:,:), &
                                plant_uptake_nh4_sl(:,:), &
                                plant_uptake_no3_sl(:,:), &
                                plant_uptake_po4_sl(:,:), &
                                mycorrhiza_uptake_nh4_sl(:,:), &
                                mycorrhiza_uptake_no3_sl(:,:), &
                                mycorrhiza_uptake_po4_sl(:,:), &
                                nitrification_no3(:,:), &
                                denitrification_n2(:,:), &
                                biochem_mineralisation_po4(:,:), &
                                sb_loss, &
                                microbial_uptake_nh4_sl(:,:), &
                                microbial_uptake_po4_sl(:,:))

    IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme)=='dynamic')THEN
       CALL calc_partitioning_nitrification_denitrification_rate(options, &
                                                                 rtm_denitrification_act(:,:), &
                                                                 nitrification_no3(:,:), &
                                                                 denitrification_n2(:,:), &
                                                                 nitrification_noy(:,:), &
                                                                 nitrification_n2o(:,:), &
                                                                 denitrification_noy(:,:), &
                                                                 denitrification_n2o(:,:))
    ENDIF

    CALL calc_sb_inorganic_n15_fluxes(options, & 
                                          nh4_solute, &
                                          nh4_n15_solute, &
                                          no3_solute, &
                                          no3_n15_solute, &
                                          plant_uptake_nh4_sl, & 
                                          plant_uptake_no3_sl, &      
                                          mycorrhiza_uptake_nh4_sl, & 
                                          mycorrhiza_uptake_no3_sl, & 
                                          microbial_uptake_nh4_sl, &        
                                          microbial_uptake_no3_sl, &        
                                          nitrification_no3, &        
                                          nitrification_noy, &        
                                          nitrification_n2o, &        
                                          denitrification_noy, &      
                                          denitrification_n2o, &      
                                          denitrification_n2, &       
                                          transport_nh4_n15_solute, & 
                                          transport_no3_n15_solute, &    
                                          plant_uptake_nh4_n15_sl, &      
                                          plant_uptake_no3_n15_sl, &      
                                          mycorrhiza_uptake_nh4_n15_sl, & 
                                          mycorrhiza_uptake_no3_n15_sl, & 
                                          microbial_uptake_nh4_n15_sl, &        
                                          microbial_uptake_no3_n15_sl, &        
                                          nitrification_no3_n15, &        
                                          nitrification_noy_n15, &        
                                          nitrification_n2o_n15, &        
                                          denitrification_noy_n15, &      
                                          denitrification_n2o_n15, &      
                                          denitrification_n2_n15) 
    
    !> 3.2 formation of SOM and associated respiration, net mineralisation fluxes
    !!          
    CALL calc_SOM_formation(options, &
                            apparent_mic_cn(:,:), &
                            microbial_uptake_nh4_sl(:,:), &
                            microbial_uptake_nh4_n15_sl, &        
                            microbial_uptake_po4_sl(:,:), &
                            sb_loss, &
                            het_respiration(:,:), &
                            het_respiration_c13(:,:), &
                            het_respiration_c14(:,:), &
                            net_mineralisation_nh4(:,:), &
                            net_mineralisation_po4(:,:), &
                            net_mineralisation_nh4_n15(:,:), &
                            sb_formation )

    !> 3.3 formation of metabolic and strucutral litter by woody liter decay and 
    !!     associated respiration, net mineralisation fluxes    
    IF(lctlib%growthform==itree)THEN 
       CALL calc_litter_formation_from_woody_decay(options, &
                                                   sb_loss, &
                                                   het_respiration(:,:), &
                                                   het_respiration_c13(:,:), &
                                                   het_respiration_c14(:,:), &
                                                   net_mineralisation_nh4(:,:), &
                                                   net_mineralisation_po4(:,:), &
                                                   net_mineralisation_nh4_n15(:,:), &
                                                   sb_formation )
    ENDIF


    !> 4.0 Mycorrhiza dynamics 
    !! 
    CALL calc_mycorrhiza_dynamics(options, &
                                  lctlib, &
                                  dsl4jsb_Config(SB_)%flag_myco_organic, &
                                  soil_depth_sl(:,:), &
                                  mycorrhiza_uptake_nh4_sl(:,:), &
                                  mycorrhiza_uptake_no3_sl(:,:), &
                                  mycorrhiza_uptake_norg_sl(:,:), &
                                  mycorrhiza_uptake_po4_sl(:,:), &
                                  sb_pool, &
                                  sb_formation, &
                                  sb_loss, &
                                  het_respiration(:,:), &
                                  het_respiration_c13(:,:), &
                                  het_respiration_c14(:,:), &
                                  sb_mycorrhiza_export)


    !> 5.0 Update soil litter pools with litter fall
    !!
    CALL calc_litter_partitioning(options, &
                                  lctlib, &  
                                  soil_depth_sl(:,:), &
                                  root_fraction_act_sl(:,:), &
                                  veg_litterfall, &
                                  sb_formation )

    !> 6.0 Add asymbiontic N fixation to NH4 pool 
    !! 
    CALL calc_asymbiontic_n_fixation(options, &
                                     lctlib, &
                                     TRIM(dsl4jsb_Config(VEG_)%bnf_scheme), & 
                                     soil_depth_sl(:,:), &
                                     nh4_solute(:,:), &
                                     no3_solute(:,:), &
                                     sb_pool%microbial%carbon(:,:), &
                                     root_fraction_act_sl(:,:), &
                                     asymb_n_fixation(:,:), &
                                     asymb_n_fixation_n15(:,:))
                    
    !> 6.1 For sb_model = fixed, estimate N losses
    !! 
    IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme)=='fixed')THEN
       CALL calc_fixed_nloss_rate(options, &
                                  soil_depth_sl(:,:), &
                                  net_mineralisation_nh4(:,:), & 
                                  net_mineralisation_nh4_n15(:,:), & 
                                  nitrification_no3(:,:), &
                                  nitrification_no3_n15(:,:), &
                                  volatilisation_nh4(:,:), &
                                  volatilisation_nh4_n15(:,:), &
                                  emission_n2(:), &
                                  emission_n2_n15(:))
    ENDIF

    !> 7.0 convert rates to µmol/ m3 / s (µmol / m2 / s for surface fluxes)
    !!
    het_respiration(:,:) = het_respiration(:,:) / dtime * 1000000._wp
    het_respiration_c13(:,:) = het_respiration_c13(:,:) / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) / dtime * 1000000._wp

    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) / dtime * 1000000._wp
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_sl(:,:) = plant_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_no3_sl(:,:) = plant_uptake_no3_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_po4_sl(:,:) = plant_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_nh4_sl(:,:) = mycorrhiza_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_no3_sl(:,:) = mycorrhiza_uptake_no3_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_norg_sl(:,:) = mycorrhiza_uptake_norg_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_po4_sl(:,:) = mycorrhiza_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    microbial_uptake_nh4_sl(:,:) = microbial_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    microbial_uptake_no3_sl(:,:) = microbial_uptake_no3_sl(:,:) / dtime * 1000000._wp 
    microbial_uptake_po4_sl(:,:) = microbial_uptake_po4_sl(:,:) / dtime * 1000000._wp 

    volatilisation_nh4(:,:) = volatilisation_nh4(:,:) / dtime * 1000000._wp
    nitrification_no3(:,:) = nitrification_no3(:,:) / dtime * 1000000._wp 
    nitrification_noy(:,:) = nitrification_noy(:,:) / dtime * 1000000._wp 
    nitrification_n2o(:,:) = nitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_noy(:,:) = denitrification_noy(:,:) / dtime * 1000000._wp 
    denitrification_n2o(:,:) = denitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_n2(:,:) = denitrification_n2(:,:) / dtime * 1000000._wp 

    transport_nh4_solute(:,:) = transport_nh4_solute(:,:) / dtime * 1000000._wp
    transport_no3_solute(:,:) = transport_no3_solute(:,:) / dtime * 1000000._wp
    transport_po4_solute(:,:) = transport_po4_solute(:,:) / dtime * 1000000._wp
    transport_po4_assoc_fast(:,:) = transport_po4_assoc_fast(:,:) / dtime * 1000000._wp
    transport_po4_assoc_slow(:,:) = transport_po4_assoc_slow(:,:) / dtime * 1000000._wp
    transport_po4_occluded(:,:) = transport_po4_occluded(:,:) / dtime * 1000000._wp
    transport_noy(:,:) = transport_noy(:,:) / dtime * 1000000._wp
    transport_n2o(:,:) = transport_n2o(:,:) / dtime * 1000000._wp
    transport_n2(:,:) = transport_n2(:,:) / dtime * 1000000._wp

    leaching_nh4_solute(:) = leaching_nh4_solute(:) / dtime * 1000000._wp
    leaching_no3_solute(:) = leaching_no3_solute(:) / dtime * 1000000._wp
    leaching_po4_solute(:) = leaching_po4_solute(:) / dtime * 1000000._wp
    emission_noy(:) = emission_noy(:) / dtime * 1000000._wp
    emission_n2o(:) = emission_n2o(:) / dtime * 1000000._wp
    emission_n2(:) = emission_n2(:) / dtime * 1000000._wp

    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_n15_sl(:,:) = plant_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_no3_n15_sl(:,:) = plant_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_nh4_n15_sl(:,:) = microbial_uptake_nh4_n15_sl(:,:) / dtime * 1000000._wp
    microbial_uptake_no3_n15_sl(:,:) = microbial_uptake_no3_n15_sl(:,:) / dtime * 1000000._wp

    volatilisation_nh4_n15(:,:) = volatilisation_nh4_n15(:,:) / dtime * 1000000._wp
    nitrification_no3_n15(:,:) = nitrification_no3_n15(:,:) / dtime * 1000000._wp 
    nitrification_noy_n15(:,:) = nitrification_noy_n15(:,:) / dtime * 1000000._wp 
    nitrification_n2o_n15(:,:) = nitrification_n2o_n15(:,:) / dtime * 1000000._wp 
    denitrification_noy_n15(:,:) = denitrification_noy_n15(:,:) / dtime * 1000000._wp 
    denitrification_n2o_n15(:,:) = denitrification_n2o_n15(:,:) / dtime * 1000000._wp 
    denitrification_n2_n15(:,:) = denitrification_n2_n15(:,:) / dtime * 1000000._wp 

    transport_nh4_n15_solute(:,:) = transport_nh4_n15_solute(:,:) / dtime * 1000000._wp
    transport_no3_n15_solute(:,:) = transport_no3_n15_solute(:,:) / dtime * 1000000._wp
    transport_noy_n15(:,:) = transport_noy_n15(:,:) / dtime * 1000000._wp
    transport_n2o_n15(:,:) = transport_n2o_n15(:,:) / dtime * 1000000._wp
    transport_n2_n15(:,:) = transport_n2_n15(:,:) / dtime * 1000000._wp

    leaching_nh4_n15_solute(:) = leaching_nh4_n15_solute(:) / dtime * 1000000._wp
    leaching_no3_n15_solute(:) = leaching_no3_n15_solute(:) / dtime * 1000000._wp
    emission_noy_n15(:) = emission_noy_n15(:) / dtime * 1000000._wp
    emission_n2o_n15(:) = emission_n2o_n15(:) / dtime * 1000000._wp
    emission_n2_n15(:) = emission_n2_n15(:) / dtime * 1000000._wp

    weathering_po4(:,:) = weathering_po4(:,:) / dtime * 1000000._wp
    occlusion_po4(:,:) = occlusion_po4(:,:) / dtime * 1000000._wp
    slow_sorption_po4(:,:) = slow_sorption_po4(:,:) / dtime * 1000000._wp
    biochem_mineralisation_po4(:,:) = biochem_mineralisation_po4(:,:) / dtime * 1000000._wp

! Pass back values to process memory
dsl4slm_passback_vegpool_var(pool,fine_root,carbon)
#define Passback_veg_litterfall
#define Passback_sb_pool
#define Passback_sb_formation
#define Passback_sb_loss
#define Passback_sb_transport
#define Passback_sb_mycorrhiza_export
#include "slm_pools_a_incl.f90"
#undef Passback_veg_litterfall
#undef Passback_sb_pool
#undef Passback_sb_formation
#undef Passback_sb_loss
#undef Passback_sb_transport
#undef Passback_sb_mycorrhiza_export

  END SUBROUTINE update_sb_simple_model 


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Routine to provide soil profile values for "simple_0d" and "simple_1d" (1 and many soil layers)
  !!
  !! For the 0D model, all processes happen only in the first layer of the soil \n
  !! therefore the soil temperature and moisture modifiers are aggregated 
  !! across the profile (assuming that the weighting follows the prescribed root profile). \n
  !! The kms are expressed as mol/m3 or m3/mol and therefore need to be adjusted to avoid
  !! the implicit assumption that all material would be concentrated in the top layer. \n
  !! This is done by assuming the material is distributed within the average rooting depth \n
  !!
  !! For the 1D model, values are simply copied.
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_vertical_profile( options, & 
                                    lctlib, &                             ! in 
                                    sb_model_scheme, &
                                    soil_depth_sl, &
                                    bulk_dens_sl, &
                                    nh4_solute, &
                                    root_fraction_sl, &
                                    rtm_decomposition, &
                                    rmm_decomposition, &
                                    rtm_hsc, &
                                    rmm_hsc, &
                                    rtm_nitrification, &
                                    rmm_nitrification, &
                                    rtm_denitrification, &
                                    rtm_gasdiffusion, &
                                    rmm_gasdiffusion, &
                                    rtm_sorption, &
                                    rtm_desorption, &
                                    rmm_sorption, &
                                    rmm_desorption, &
                                    bulk_dens_corr_sl, &
                                    p_deposition, &
                                    anaerobic_volume_fraction_sl, &
                                    root_fraction_act_sl, &              ! out
                                    km1_uptake_nh4_act, &
                                    km2_uptake_nh4_act, &
                                    km1_uptake_no3_act, &
                                    km2_uptake_no3_act, &
                                    km1_uptake_po4_act, &
                                    km2_uptake_po4_act, &
                                    km_denitrification_c_act, &
                                    km_denitrification_no3_act, &
                                    apparent_mic_cn, &
                                    rtm_decomposition_act, &
                                    rmm_decomposition_act, &
                                    rtm_nitrification_act, &
                                    rmm_nitrification_act, &
                                    rtm_denitrification_act, &
                                    rtm_gasdiffusion_act, &
                                    rmm_gasdiffusion_act, &
                                    p_deposition_act, &
                                    k_weath_mineral_act, &
                                    k_adsorpt_po4_act, &
                                    k_desorpt_po4_act, &
                                    km_dip_biochem_po4_act, &
                                    anaerobic_volume_fraction_act )

    USE mo_sb_constants,            ONLY: km_denitrification_c, km_denitrification_no3, km1_up_nh4, km1_up_no3, & 
                                          km1_up_po4, km2_up_nh4, km2_up_no3, km2_up_po4,microbial_cn_max,microbial_cn_min, &
                                          f_mic_cn, k_weath_mineral, k_desorpt_po4, k_adsorpt_po4, km_dip_biochem_po4

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)  :: options
    TYPE(t_lctlib_element),                      INTENT(in)  :: lctlib                       !< land-cover-type library - parameter across pft's
    CHARACTER(len=*),                            INTENT(in)  :: sb_model_scheme              !< from sb config
    REAL(wp), DIMENSION(options%nc),             INTENT(in)  :: p_deposition                 !< p deposition, [µmol m-2 s-1]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: &
                                                                soil_depth_sl, &             !< soil depth of a layer [m]
                                                                bulk_dens_sl, &              !< bulk density of a layer [kg m-3] 
                                                                nh4_solute, &                !< NH4 content of a layer [mol m-3] 
                                                                root_fraction_sl, &          !< fraction of root mass in layer
                                                                rtm_decomposition, &         !< temperature response scalar for decomposition
                                                                rmm_decomposition, &         !< moisture response scalar for decomposition
                                                                rtm_hsc, &                   !< temperature response scalar for half-saturation constants 
                                                                rmm_hsc, &                   !< moisture response scalar for half-saturation constants
                                                                rtm_nitrification, &         !< temperature response scalar for nitrification
                                                                rmm_nitrification, &         !< moisture response scalar for nitrification
                                                                rtm_denitrification, &       !< temperature response scalar for denitrification
                                                                rtm_gasdiffusion, &          !< temperature response scalar for gas diffusion
                                                                rmm_gasdiffusion, &          !< moisture response scalar for gas diffusion
                                                                rtm_sorption, &              !< temperature response scalar for po4 sorption
                                                                rtm_desorption, &            !< temperature response scalar for po4 desorption
                                                                rmm_sorption, &              !< moisture response scalar for po4 sorption
                                                                rmm_desorption, &            !< moisture response scalar for po4 desorption
                                                                bulk_dens_corr_sl, &         !< bulk density of the soil layer, [kg/m3]
                                                                anaerobic_volume_fraction_sl !< anaerobic volume fraction
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out) :: &
                                                                root_fraction_act_sl, &      !< actual fraction of root mass in layer
                                                                km1_uptake_nh4_act, &        !< actual km1 for NH4 uptake [m3 / mol]
                                                                km2_uptake_nh4_act, &        !< actual km2 for NH4 uptake [mol / m3]
                                                                km1_uptake_no3_act, &        !< actual km1 for NO3 uptake [m3 / mol]
                                                                km2_uptake_no3_act, &        !< actual km2 for NO3 uptake [mol / m3]
                                                                km1_uptake_po4_act, &        !< actual km1 for PO4 uptake [m3 / mol]
                                                                km2_uptake_po4_act, &        !< actual km2 for PO4 uptake [mol / m3]
                                                                km_denitrification_c_act, &  !< actual km on C for denitrification [mol / m3]
                                                                km_denitrification_no3_act, & !< actual km on NO3 for denitrification [mol / m3]
                                                                apparent_mic_cn, &           !< apparent microbial molar carbon-nitrogen ratio 
                                                                rtm_decomposition_act, &     !< actual temperature response scalar for decomposition
                                                                rmm_decomposition_act, &     !< actual moisture response scalar for decomposition
                                                                rtm_nitrification_act, &     !< actual temperature response scalar for nitrification
                                                                rmm_nitrification_act, &     !< actual moisture response scalar for nitrification
                                                                rtm_denitrification_act, &   !< actual temperature response scalar for denitrification
                                                                rtm_gasdiffusion_act, &      !< actual temperature response scalar for gas diffusion
                                                                rmm_gasdiffusion_act, &      !< actual moisture response scalar for gas diffusion
                                                                p_deposition_act, &          !< actual p deposition, [mol m-3 timestep-1]
                                                                k_weath_mineral_act, &       !< actual weathering rate for primary P [s-1]
                                                                k_adsorpt_po4_act, &         !< actual adsorption rate for fast_assoc_po4 [s-1]
                                                                k_desorpt_po4_act, &         !< actual desorption rate for slow_assoc_po4 [s-1]
                                                                km_dip_biochem_po4_act, &    !< actual km for soluble P control on biochem_mineralisation_po4
                                                                anaerobic_volume_fraction_act!< actual anaerobic volume fraction
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc) :: root_depth_avg, &       !< average rooting depth [m]
                                       rtm_hsc_act, rmm_hsc_act, & 
                                       rtm_sorption_act, rmm_sorption_act, &
                                       rtm_desorption_act, rmm_desorption_act
    INTEGER                         :: isoil                !< loop over options%nsoil_sb
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':calc_vertical_profile'

    !> 1.0 Case of 0d model
    !!
    IF(sb_model_scheme == 'simple_0d') THEN
      root_fraction_act_sl(:,1)                  = 1.0_wp
      root_fraction_act_sl(:,2:options%nsoil_sb) = 0.0_wp
      ! set to default values to avoid division by zero below (2:options%nsoil_sb are unused in 0D)
      km1_uptake_nh4_act(:,:) = km1_up_nh4 
      km1_uptake_no3_act(:,:) = km1_up_no3
      km1_uptake_po4_act(:,:) = km1_up_po4
      km2_uptake_nh4_act(:,:) = km2_up_nh4
      km2_uptake_no3_act(:,:) = km2_up_no3
      km2_uptake_po4_act(:,:) = km2_up_po4
      km_dip_biochem_po4_act (:,:) = km_dip_biochem_po4
      km_denitrification_c_act(:,:) = km_denitrification_c
      km_denitrification_no3_act(:,:) = km_denitrification_no3
      root_depth_avg(:) = sum(root_fraction_sl(:,:)*soil_depth_sl(:,:),DIM=2)
      WHERE(root_depth_avg(:)>eps8) 
         km1_uptake_nh4_act(:,1) = km1_up_nh4 / (rtm_hsc(:,1) * rmm_hsc(:,1)) * soil_depth_sl(:,1) / root_depth_avg(:)
         km1_uptake_no3_act(:,1) = km1_up_no3 / (rtm_hsc(:,1) * rmm_hsc(:,1)) * soil_depth_sl(:,1) / root_depth_avg(:)
         km1_uptake_po4_act(:,1) = km1_up_po4 / (rtm_hsc(:,1) * rmm_hsc(:,1)) * soil_depth_sl(:,1) / root_depth_avg(:)
         km2_uptake_nh4_act(:,1) = km2_up_nh4 * rtm_hsc(:,1) * rmm_hsc(:,1) / soil_depth_sl(:,1) * root_depth_avg(:)
         km2_uptake_no3_act(:,1) = km2_up_no3 * rtm_hsc(:,1) * rmm_hsc(:,1) / soil_depth_sl(:,1) * root_depth_avg(:)
         km2_uptake_po4_act(:,1) = km2_up_po4 * rtm_hsc(:,1) * rmm_hsc(:,1) / soil_depth_sl(:,1) * root_depth_avg(:)
         km_dip_biochem_po4_act(:,1) = km_dip_biochem_po4 * rtm_hsc(:,1) * rmm_hsc(:,1) / soil_depth_sl(:,1) * root_depth_avg(:)
         km_denitrification_c_act(:,1) = km_denitrification_c / soil_depth_sl(:,1) * root_depth_avg(:)
         km_denitrification_no3_act(:,1) = km_denitrification_no3 / soil_depth_sl(:,1) * root_depth_avg(:)
      ENDWHERE
      apparent_mic_cn(:,:) = microbial_cn_max
      apparent_mic_cn(:,1) = MAX(microbial_cn_max - f_mic_cn * nh4_solute(:,1) * soil_depth_sl(:,1) / root_depth_avg(:) / bulk_dens_sl(:,1), &
                                 microbial_cn_min)

      rtm_decomposition_act(:,:) = 0.0_wp 
      rmm_decomposition_act(:,:) = 0.0_wp 
      rtm_nitrification_act(:,:) = 0.0_wp 
      rmm_nitrification_act(:,:) = 0.0_wp 
      rtm_denitrification_act(:,:) = 0.0_wp 
      rtm_gasdiffusion_act(:,:) = 0.0_wp 
      rmm_gasdiffusion_act(:,:) = 0.0_wp 
      anaerobic_volume_fraction_act(:,:) = 0.0_wp 
      rtm_hsc_act(:) = 0.0_wp 
      rmm_hsc_act(:) = 0.0_wp 
      rtm_sorption_act(:) = 0.0_wp 
      rmm_sorption_act(:) = 0.0_wp
      rtm_desorption_act(:) = 0.0_wp 
      rmm_desorption_act(:) = 0.0_wp
      k_weath_mineral_act(:,:) = 0.0_wp 
      k_adsorpt_po4_act(:,:) = 0.0_wp 
      k_desorpt_po4_act(:,:) = 0.0_wp 
      DO isoil = 1,options%nsoil_sb
         rtm_decomposition_act(:,1) = rtm_decomposition_act(:,1) + & 
                                         rtm_decomposition(:,isoil) * root_fraction_sl(:,isoil)
         rmm_decomposition_act(:,1) = rmm_decomposition_act(:,1) + & 
                                         rmm_decomposition(:,isoil) * root_fraction_sl(:,isoil)
         rtm_nitrification_act(:,1) = rtm_nitrification_act(:,1) + & 
                                         rtm_nitrification(:,isoil) * root_fraction_sl(:,isoil)
         rmm_nitrification_act(:,1) = rmm_nitrification_act(:,1) + & 
                                         rmm_nitrification(:,isoil) * root_fraction_sl(:,isoil)
         rtm_denitrification_act(:,1) = rtm_denitrification_act(:,1) + & 
                                         rtm_denitrification(:,isoil) * root_fraction_sl(:,isoil)
         rtm_gasdiffusion_act(:,1) = rtm_gasdiffusion_act(:,1) + & 
                                         rtm_gasdiffusion(:,isoil) * root_fraction_sl(:,isoil)
         rmm_gasdiffusion_act(:,1) = rmm_gasdiffusion_act(:,1) + & 
                                         rmm_gasdiffusion(:,isoil) * root_fraction_sl(:,isoil)
         anaerobic_volume_fraction_act(:,1) = anaerobic_volume_fraction_act(:,1) + & 
                                         anaerobic_volume_fraction_sl(:,isoil) * root_fraction_sl(:,isoil)
         rtm_hsc_act(:) = rtm_hsc_act(:) + & 
                                         rtm_hsc(:,isoil) * root_fraction_sl(:,isoil)
         rmm_hsc_act(:) = rmm_hsc_act(:) + & 
                                         rmm_hsc(:,isoil) * root_fraction_sl(:,isoil)
         rtm_sorption_act(:) = rtm_sorption_act(:) + & 
                                         rtm_sorption(:,isoil) * root_fraction_sl(:,isoil) 
         rmm_sorption_act(:) = rmm_sorption_act(:) + & 
                                         rmm_sorption(:,isoil) * root_fraction_sl(:,isoil) 
         rtm_desorption_act(:) = rtm_desorption_act(:) + & 
                                         rtm_desorption(:,isoil) * root_fraction_sl(:,isoil) 
         rmm_desorption_act(:) = rmm_desorption_act(:) + & 
                                         rmm_desorption(:,isoil) * root_fraction_sl(:,isoil) 
      ENDDO
      k_weath_mineral_act(:,1) = rtm_hsc_act(:) / rmm_hsc_act(:) * k_weath_mineral
      k_adsorpt_po4_act(:,1) = rtm_sorption_act(:) * rmm_sorption_act(:) * k_adsorpt_po4
      k_desorpt_po4_act(:,1) = rtm_desorption_act(:) * rmm_desorption_act(:) * k_desorpt_po4
    !> 2.0 Case of 1d model, simply copy values
    !!
    ELSE 
      root_fraction_act_sl(:,:) = root_fraction_sl(:,:)
      km1_uptake_nh4_act(:,:) = km1_up_nh4 / (rtm_hsc(:,:) * rmm_hsc(:,:))
      km1_uptake_no3_act(:,:) = km1_up_no3 / (rtm_hsc(:,:) * rmm_hsc(:,:))
      km1_uptake_po4_act(:,:) = km1_up_po4 / (rtm_hsc(:,:) * rmm_hsc(:,:))
      km2_uptake_nh4_act(:,:) = km2_up_nh4 * rtm_hsc(:,:) * rmm_hsc(:,:)
      km2_uptake_no3_act(:,:) = km2_up_no3 * rtm_hsc(:,:) * rmm_hsc(:,:)
      km2_uptake_po4_act(:,:) = km2_up_po4 * rtm_hsc(:,:) * rmm_hsc(:,:)
      km_dip_biochem_po4_act(:,:) = km_dip_biochem_po4 * rtm_hsc(:,:) * rmm_hsc(:,:)
      km_denitrification_c_act(:,:) = km_denitrification_c
      km_denitrification_no3_act(:,:) = km_denitrification_no3
      apparent_mic_cn(:,:) = MAX(microbial_cn_max - f_mic_cn * nh4_solute(:,:)/bulk_dens_sl(:,:), &
                                 microbial_cn_min)

      rtm_decomposition_act(:,:) = rtm_decomposition(:,:)
      rmm_decomposition_act(:,:) = rmm_decomposition(:,:) 
      rtm_nitrification_act(:,:) = rtm_nitrification(:,:)
      rmm_nitrification_act(:,:) = rmm_nitrification(:,:) 
      rtm_denitrification_act(:,:) = rtm_denitrification(:,:)
      rtm_gasdiffusion_act(:,:) = rtm_gasdiffusion(:,:)
      rmm_gasdiffusion_act(:,:) = rmm_gasdiffusion(:,:) 
      k_weath_mineral_act(:,:) = rtm_hsc(:,:) / rmm_hsc(:,:) * k_weath_mineral
      k_adsorpt_po4_act(:,:) = rtm_sorption(:,:) * rmm_sorption(:,:) * k_adsorpt_po4
      k_desorpt_po4_act(:,:) = rtm_desorption(:,:) * rmm_desorption(:,:) * k_desorpt_po4
      anaerobic_volume_fraction_act(:,:) = anaerobic_volume_fraction_sl(:,:) 

    ENDIF

    p_deposition_act(:,:) = 0.0_wp
    p_deposition_act(:,1) = p_deposition(:) / 1.e6_wp / soil_depth_sl(:,1) * dtime

  END SUBROUTINE calc_vertical_profile


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates potential plant uptake of NH4, NO3 and PO4, given soil concentrations,
  !! and fine root mass per layer 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_plant_uptake_rate(options, &
                                              lctlib, &                     ! in
                                              soil_depth_sl, &
                                              veg_pool_fine_root_carbon, &
                                              root_fraction_act_sl, &
                                              rtm_root_uptake_act, &
                                              km1_uptake_nh4_act, &
                                              km2_uptake_nh4_act, &
                                              km1_uptake_no3_act, &
                                              km2_uptake_no3_act, &
                                              km1_uptake_po4_act, &
                                              km2_uptake_po4_act, &
                                              nh4_solute, &
                                              no3_solute, &
                                              po4_solute, &
                                              plant_uptake_nh4_sl, &        ! out
                                              plant_uptake_no3_sl, &
                                              plant_uptake_po4_sl, &
                                              unit_uptake_n_pot, &
                                              unit_uptake_p_pot  )
  
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)  :: options
    TYPE(t_lctlib_element),                      INTENT(in)  :: lctlib                     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc),             INTENT(in)  :: veg_pool_fine_root_carbon  !< fine root carbon [mol / m2]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: &
                                                                soil_depth_sl, &           !< soil depth per layer [m]
                                                                root_fraction_act_sl, &    !< fraction of fine roots per soil layer
                                                                rtm_root_uptake_act, &     !< temperature response of root uptake (taken identical to that of decomposition) 
                                                                km1_uptake_nh4_act, &      !< actual km1 for NH4 uptake [m3 / mol]
                                                                km2_uptake_nh4_act, &      !< actual km2 for NH4 uptake [mol / m3]
                                                                km1_uptake_no3_act, &      !< actual km1 for NO3 uptake [m3 / mol]
                                                                km2_uptake_no3_act, &      !< actual km2 for NO3 uptake [mol / m3]
                                                                km1_uptake_po4_act, &      !< actual km1 for PO4 uptake [m3 / mol]
                                                                km2_uptake_po4_act, &      !< actual km2 for PO4 uptake [mol / m3]
                                                                nh4_solute, &              !< NH4 solute concentration [mol / m3] 
                                                                no3_solute, &              !< NO3 solute concentration [mol / m3]
                                                                po4_solute                 !< PO4 solute concentration [mol / m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out) :: &
                                                                plant_uptake_nh4_sl, &     !< potential plant uptake of NH4 [mol / m3 / timestep]
                                                                plant_uptake_no3_sl, &     !< potential plant uptake of NO3 [mol / m3 / timestep]
                                                                plant_uptake_po4_sl        !< potential plant uptake of PO4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc),             INTENT(out) :: unit_uptake_n_pot, &       !< potential N uptake by fine roots per unit mass [µmol / mol / s]
                                                                unit_uptake_p_pot          !< potential P uptake by fine roots per unit mass [µmol / mol / s]
    ! 0.2 Local
    INTEGER                                          :: isoil        !< loop over options%nsoil_sb
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_potential_plant_uptake_rate'

    unit_uptake_n_pot(:) = 0.0_wp
    unit_uptake_p_pot(:) = 0.0_wp

    DO isoil = 1,options%nsoil_sb
      !> 1.0 inorganic plant N uptake 
      !!
      !! potential N uptake rate (µmol N / mol C fine-root / s), weighted by root distribution
      hlp1(:,isoil) = root_fraction_act_sl(:,isoil) * rtm_root_uptake_act(:,isoil) * lctlib%vmax_uptake_n * & 
                  nh4_solute(:,isoil) * ( km1_uptake_nh4_act(:,isoil) + & 
                  1._wp / ( km2_uptake_nh4_act(:,isoil) + nh4_solute(:,isoil) )) 
      hlp2(:,isoil) = root_fraction_act_sl(:,isoil) * rtm_root_uptake_act(:,isoil) * lctlib%vmax_uptake_n * & 
                  no3_solute(:,isoil) * ( km1_uptake_no3_act(:,isoil) + & 
                  1._wp / ( km2_uptake_no3_act(:,isoil) + no3_solute(:,isoil) )) 
      unit_uptake_n_pot(:) = unit_uptake_n_pot(:) + hlp1(:,isoil) + hlp2(:,isoil) 
  
      ! potential uptake rate in mol N / m3 / timiestep
      plant_uptake_nh4_sl(:,isoil) = hlp1(:,isoil) * veg_pool_fine_root_carbon(:) / & 
                                     soil_depth_sl(:,isoil) * dtime/1.e6_wp
      plant_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * veg_pool_fine_root_carbon(:) / & 
                                     soil_depth_sl(:,isoil) * dtime/1.e6_wp
  
      !> 2.0 inorganic plant P uptake 
      !!
      !! potential P uptake rate (µmol P / mol C fine-root / s), weighted by root distribution
      hlp1(:,isoil) = root_fraction_act_sl(:,isoil)  * rtm_root_uptake_act(:,isoil) * lctlib%vmax_uptake_p * &
                             po4_solute(:,isoil) * ( km1_uptake_po4_act(:,isoil) + &
                             1._wp / ( km2_uptake_po4_act(:,isoil) + po4_solute(:,isoil) ))
      unit_uptake_p_pot(:) = unit_uptake_p_pot(:) + hlp1(:,isoil)
      ! potential uptake rate in mol P / m3 / timestep
      plant_uptake_po4_sl(:,isoil) = hlp1(:,isoil) * & 
                                     veg_pool_fine_root_carbon(:) / soil_depth_sl(:,isoil) * dtime/1.e6_wp
   
    END DO    

  END SUBROUTINE calc_potential_plant_uptake_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates potential mycorrhiza uptake of NH4, NO3 and PO4, given soil concentrations,
  !! and mycorrhiza mass per layer 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_mycorrhiza_uptake_rate(options, &
                                                   lctlib, &                     ! in
                                                   flag_myco_inorganic, &
                                                   soil_depth_sl, &
                                                   sb_pool_mycorrhiza_carbon, &
                                                   sb_pool_mycorrhiza_nitrogen, &
                                                   rtm_mycorrhiza_uptake_act, &
                                                   km1_uptake_nh4_act, &
                                                   km2_uptake_nh4_act, &
                                                   km1_uptake_no3_act, &
                                                   km2_uptake_no3_act, &
                                                   km1_uptake_po4_act, &
                                                   km2_uptake_po4_act, &
                                                   nh4_solute, &
                                                   no3_solute, &
                                                   po4_solute, &
                                                   mycorrhiza_uptake_nh4_sl, &   ! out
                                                   mycorrhiza_uptake_no3_sl, &
                                                   mycorrhiza_uptake_po4_sl )
  
    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: f_maxresp_myc, f_surface_myc2fineroot

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)  :: options                     !< dimensions...
    TYPE(t_lctlib_element),                           INTENT(in)  :: lctlib                      !< land-cover-type library - parameter across pft's
    LOGICAL,                                          INTENT(in)  :: flag_myco_inorganic           !< from sb config
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: soil_depth_sl, &            !< soil depth per layer [m]
                                                                     sb_pool_mycorrhiza_carbon, &!< mycorrhizal biomass [mol C / m3]
                                                                     sb_pool_mycorrhiza_nitrogen, &!< mycorrhizal biomass [mol C / m3]
                                                                     rtm_mycorrhiza_uptake_act, &!< temperature response of mycorrhiza uptake (taken identical to that of decomposition) 
                                                                     km1_uptake_nh4_act, &       !< actual km1 for NH4 uptake [m3 / mol]
                                                                     km2_uptake_nh4_act, &       !< actual km2 for NH4 uptake [mol / m3]
                                                                     km1_uptake_no3_act, &       !< actual km1 for NO3 uptake [m3 / mol]
                                                                     km2_uptake_no3_act, &       !< actual km2 for NO3 uptake [mol / m3]
                                                                     km1_uptake_po4_act, &       !< actual km1 for PO4 uptake [m3 / mol]
                                                                     km2_uptake_po4_act, &       !< actual km2 for PO4 uptake [mol / m3]
                                                                     nh4_solute, &               !< NH4 solute concentration [mol / m3] 
                                                                     no3_solute, &               !< NO3 solute concentration [mol / m3]
                                                                     po4_solute                  !< PO4 solute concentration [mol / m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out) :: mycorrhiza_uptake_nh4_sl, & !< potential mycorrhiza uptake of NH4 [mol / m3 / timestep]
                                                                     mycorrhiza_uptake_no3_sl, & !< potential mycorrhiza uptake of NO3 [mol / m3 / timestep]
                                                                     mycorrhiza_uptake_po4_sl    !< potential mycorrhiza uptake of PO4 [mol / m3 / timestep]

    ! 0.2 Local
    INTEGER                                          :: isoil        !< loop over options%nsoil_sb
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_potential_mycorrhiza_uptake_rate'

    

    !> 1.0 pot. N uptake by mycorrhizal fungi (mol / m3 / timestep)
    !!
    IF(flag_myco_inorganic)THEN
      DO isoil = 1,options%nsoil_sb
        mycorrhiza_uptake_nh4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * f_surface_myc2fineroot * lctlib%vmax_uptake_n * & 
                                          nh4_solute(:,isoil) * ( km1_uptake_nh4_act(:,isoil) + & 
                                          1._wp / ( km2_uptake_nh4_act(:,isoil) + nh4_solute(:,isoil) )) * & 
                                          dtime/1.e6_wp
        mycorrhiza_uptake_no3_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * f_surface_myc2fineroot * lctlib%vmax_uptake_n * & 
                                          no3_solute(:,isoil) * ( km1_uptake_no3_act(:,isoil) + & 
                                          1._wp / ( km2_uptake_no3_act(:,isoil) + no3_solute(:,isoil) )) * & 
                                          dtime/1.e6_wp

        !> 1.1 constrain potential N uptake by avaiable C for N transformation to organic
        !!
        hlp1(:,isoil) = mycorrhiza_uptake_nh4_sl(:,isoil) * transform_cost_nh4 + & 
                      mycorrhiza_uptake_no3_sl(:,isoil) * transform_cost_no3 

        WHERE(hlp1(:,isoil)  > sb_pool_mycorrhiza_carbon(:,isoil) * f_maxresp_myc / one_day * dtime )
          WHERE(hlp1(:,isoil) > eps8)
            hlp2(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * f_maxresp_myc / one_day * dtime / & 
                            (hlp1(:,isoil)) 
            mycorrhiza_uptake_nh4_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_nh4_sl(:,isoil)
            mycorrhiza_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_no3_sl(:,isoil)
          ELSEWHERE
            mycorrhiza_uptake_nh4_sl(:,isoil) = 0.0_wp
            mycorrhiza_uptake_no3_sl(:,isoil) = 0.0_wp
          ENDWHERE
        ENDWHERE 

        !> 2.0 pot. P uptake by mycorrhizal fungi (mol / m3 / timestep)
        !!
        mycorrhiza_uptake_po4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          rtm_mycorrhiza_uptake_act(:,isoil) * f_surface_myc2fineroot * lctlib%vmax_uptake_p * & 
                                          po4_solute(:,isoil) * ( km1_uptake_po4_act(:,isoil) + & 
                                          1._wp / ( km2_uptake_po4_act(:,isoil) + po4_solute(:,isoil) )) * & 
                                          dtime/1.e6_wp
      END DO
    ENDIF 
  
  END SUBROUTINE calc_potential_mycorrhiza_uptake_rate
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the potential rates of nitrification and denitrification 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_nitrification_denitrification_rate(options, &
                                                               rtm_nitrification_act, &
                                                               rmm_nitrification_act, &
                                                               rtm_denitrification_act, &
                                                               km_denitrification_c_act, &
                                                               km_denitrification_no3_act, &
                                                               anaerobic_volume_fraction_act, &
                                                               nh4_solute, &
                                                               no3_solute, &
                                                               sb_pool_microbial_carbon, &
                                                               nitrification_no3, &
                                                               denitrification_n2)

    USE mo_sb_constants,     ONLY: vmax_nitrification, vmax_denitrification

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)  :: options
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: rtm_nitrification_act, &           !< actual temperature response for nitrification
                                                                     rmm_nitrification_act, &           !< actual moisture response for nitrification
                                                                     rtm_denitrification_act, &         !< actual temperature response for denitrification 
                                                                     km_denitrification_c_act, &        !< actual km on C for denitrification
                                                                     km_denitrification_no3_act, &      !< actual km on NO3 for denitrification
                                                                     anaerobic_volume_fraction_act, &   !< actual anaerobic volume fraction
                                                                     nh4_solute, &                      !< NH4 in solution [mol / m3]
                                                                     no3_solute, &                      !< NO3 in solution [mol / m3]
                                                                     sb_pool_microbial_carbon           !< microbial C [mol / m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out) :: nitrification_no3, &               !< potential nitrification rate [mol / m3 / timestep]
                                                                     denitrification_n2                 !< potential denitrification rate [mol / m3 / timestep] 
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_nitrification_denitrification_rate'

    ! potential nitrification rate from NH4 to NO3 [mol/m3/timestep]
    nitrification_no3(:,:)   = vmax_nitrification * rtm_nitrification_act(:,:) * rmm_nitrification_act(:,:) * & 
                               (1._wp - anaerobic_volume_fraction_act(:,:)) * nh4_solute(:,:) * dtime/1.e6_wp 

    ! potential denitrification rate from NO3 to N2 [mol/m3/timestep]
    denitrification_n2(:,:)  = vmax_denitrification * rtm_denitrification_act(:,:) * &
                               sb_pool_microbial_carbon(:,:) / &
                               ( km_denitrification_c_act(:,:) + sb_pool_microbial_carbon(:,:))  * &
                               anaerobic_volume_fraction_act(:,:) * no3_solute(:,:) / &
                               ( km_denitrification_no3_act(:,:) + no3_solute(:,:)) * dtime/1.e6_wp 

  END SUBROUTINE calc_potential_nitrification_denitrification_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the leaching and gas transport within the soil column, 
  !! given current soil N and P concentrations, soil temperature, moisture and water percolation rate 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_inorganic_transport_rate(options, &
                                           soil_moisture_scheme, &
                                           sb_model_scheme, &
                                           soil_depth_sl, &
                                           drainage_fraction, &
                                           percolation_sl, &
                                           rtm_gasdiffusion_act, &  
                                           rmm_gasdiffusion_act, & 
                                           nh4_solute, &
                                           nh4_n15_solute, &
                                           no3_solute, &
                                           no3_n15_solute, &
                                           po4_solute, &
                                           noy, & 
                                           noy_n15, & 
                                           n2o, &
                                           n2o_n15, &
                                           n2, &
                                           n2_n15, &
                                           transport_nh4_solute, &
                                           transport_nh4_n15_solute, &
                                           transport_no3_solute, &
                                           transport_no3_n15_solute, &
                                           transport_po4_solute, &
                                           transport_noy, &     
                                           transport_noy_n15, &   
                                           transport_n2o, &     
                                           transport_n2o_n15, & 
                                           transport_n2, &      
                                           transport_n2_n15, &     
                                           leaching_nh4_solute, &
                                           leaching_nh4_n15_solute, &
                                           leaching_no3_solute, &
                                           leaching_no3_n15_solute, &
                                           leaching_po4_solute, &
                                           emission_noy, &    
                                           emission_noy_n15, &
                                           emission_n2o, &    
                                           emission_n2o_n15, & 
                                           emission_n2, &  
                                           emission_n2_n15)

    USE mo_sb_constants,     ONLY: fleach_nh4, fleach_po4
    USE mo_sb_jsm_transport, ONLY: calc_liquid_phase_transport

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)  :: options
    CHARACTER(len=*),                                 INTENT(in)  :: soil_moisture_scheme               !< type of soil model scheme: bucket / layer
    CHARACTER(len=*),                                 INTENT(in)  :: sb_model_scheme                    !< type of sb model scheme: simple0d / simple1d 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: soil_depth_sl                      !< soil layer depth [m]
    REAL(wp), DIMENSION(options%nc),                  INTENT(in)  :: drainage_fraction                  !< fraction of water lost for 0d model [s-1]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: percolation_sl, &                  !< fraction of water lost for 1d model [s-1]
                                                                     rtm_gasdiffusion_act, &            !< temperature modifier for gas diffusion
                                                                     rmm_gasdiffusion_act, &            !< moisture modifier for gas diffusion       
                                                                     nh4_solute, &                      !< NH4 in solution [mol / m3]
                                                                     nh4_n15_solute, &                  !< 15NH4 in solution [mol / m3]
                                                                     no3_solute, &                      !< NO3 in solution [mol / m3]
                                                                     no3_n15_solute, &                  !< 15NO3 in solution [mol / m3]
                                                                     po4_solute, &                      !< PO4 in solution [mol / m3]
                                                                     noy, &                             !< NOy concentration [mol / m3]
                                                                     noy_n15, &                         !< 15NOy concentration [mol / m3]
                                                                     n2o, &                             !< N2O concentration [mol / m3]
                                                                     n2o_n15, &                         !< 15N2O concentration [mol / m3]
                                                                     n2, &                              !< N2 concentration [mol / m3]
                                                                     n2_n15                             !< 15N2 concentration [mol / m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out) :: transport_nh4_solute, &            !< vertical transport of NH4 inside soil column [mol / m3 / timestep]
                                                                     transport_nh4_n15_solute, &        !< vertical transport of 15NH4 inside soil column [mol / m3 / timestep]
                                                                     transport_no3_solute, &            !< vertical transport of NO3 inside soil column [mol / m3 / timestep]
                                                                     transport_no3_n15_solute, &        !< vertical transport of 15NO3 inside soil column [mol / m3 / timestep]
                                                                     transport_po4_solute, &            !< vertical transport of PO4 inside soil column [mol / m3 / timestep]
                                                                     transport_noy, &                   !< vertical transport of NOy [mol / m3 / timestep]
                                                                     transport_noy_n15, &               !< vertical transport of 15NOy [mol / m3 / timestep]
                                                                     transport_n2o, &                   !< vertical transport of N2O [mol / m3 / timestep]
                                                                     transport_n2o_n15, &               !< vertical transport of 15N2O [mol / m3 / timestep]
                                                                     transport_n2, &                    !< vertical transport of N2 [mol / m3 / timestep]
                                                                     transport_n2_n15                   !< vertical transport of 15N2 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc),                  INTENT(out) :: leaching_nh4_solute, &             !< leaching of NH4 [mol / m2 / timestep]
                                                                     leaching_nh4_n15_solute, &         !< leaching of 15NH4 [mol / m2 / timestep]
                                                                     leaching_no3_solute, &             !< leaching of NO3 [mol / m2 / timestep]
                                                                     leaching_no3_n15_solute, &         !< leaching of 15NO3 [mol / m2 / timestep]
                                                                     leaching_po4_solute, &             !< leaching of PO4 [mol / m2 / timestep]
                                                                     emission_noy, &                    !< soil efflux of NOy [mol / m2 / timestep]
                                                                     emission_noy_n15, &                !< soil efflux of 15NOy [mol / m2 / timestep]
                                                                     emission_n2o, &                    !< soil efflux of N2O [mol / m2 / timestep]
                                                                     emission_n2o_n15, &                !< soil efflux of 15N2O [mol / m2 / timestep]
                                                                     emission_n2, &                     !< soil efflux of N2 [mol / m2 / timestep]
                                                                     emission_n2_n15                    !< soil efflux of 15N2 [mol / m2 / timestep]
    ! 0.2 Local
    ! 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_inorganic_transport_rate'

    !> 1.0. Leaching rate
    !!
    IF(sb_model_scheme == 'simple_0d')THEN
       ! For the 0D soil model, loss is proportional to the fractional waterloss of the soil column
       ! reduction of layer concentration in mol m-3 timestep-1
       transport_nh4_solute(:,1)     = - fleach_nh4 * nh4_solute(:,1) * drainage_fraction(:) * dtime
       transport_nh4_n15_solute(:,1) = - fleach_nh4 * nh4_n15_solute(:,1) * drainage_fraction(:) * dtime
       transport_no3_solute(:,1)     = - no3_solute(:,1) * drainage_fraction(:) * dtime
       transport_no3_n15_solute(:,1) = - no3_n15_solute(:,1) * drainage_fraction(:) * dtime
       transport_po4_solute(:,1)     = - po4_solute(:,1) * drainage_fraction(:) * dtime

       ! record leaching losses as [mol m-2 timestep-1] 
       leaching_nh4_solute(:)     = -transport_nh4_solute(:,1) * soil_depth_sl(:,1)
       leaching_nh4_n15_solute(:) = -transport_nh4_n15_solute(:,1) * soil_depth_sl(:,1)
       leaching_no3_solute(:)     = -transport_no3_solute(:,1) * soil_depth_sl(:,1)
       leaching_no3_n15_solute(:) = -transport_no3_n15_solute(:,1) * soil_depth_sl(:,1)
       leaching_po4_solute(:)     = -transport_po4_solute(:,1) * soil_depth_sl(:,1)

    ELSE 
       IF(soil_moisture_scheme=='layer')THEN
         ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
         ! reduction of layer concentration in mol m-3 timestep-1, leaching rate in mol m-2 timestep-1
         CALL calc_liquid_phase_transport(options, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          fleach_nh4 * nh4_solute(:,:),transport_nh4_solute(:,:),leaching_nh4_solute(:))
         CALL calc_liquid_phase_transport(options, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          fleach_nh4 * nh4_n15_solute(:,:),transport_nh4_n15_solute(:,:),leaching_nh4_n15_solute(:))
         CALL calc_liquid_phase_transport(options, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          no3_solute(:,:),transport_no3_solute(:,:),leaching_no3_solute(:))
         CALL calc_liquid_phase_transport(options, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          no3_n15_solute(:,:),transport_no3_n15_solute(:,:),leaching_no3_n15_solute(:))
         CALL calc_liquid_phase_transport(options, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          po4_solute(:,:),transport_po4_solute(:,:),leaching_po4_solute(:))
       ENDIF
    ENDIF
  
    !> 2.0 gaseous loss of NOy, N2O, and N2, for now following Xu-Ri and Prentice 2008, 
    !!thus ignoring any actual gas-transport and instead assume that from each layer a certain fraction excapes directly to the atmosphere.
    !! [mol m-3 timestep-1] 
    transport_noy(:,:)     = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * noy(:,:) / one_day * dtime
    transport_noy_n15(:,:) = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * noy_n15(:,:) / one_day * dtime
    transport_n2o(:,:)     = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2o(:,:) / one_day * dtime 
    transport_n2o_n15(:,:) = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2o_n15(:,:) / one_day * dtime 
    transport_n2(:,:)      = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2(:,:)  / one_day * dtime 
    transport_n2_n15(:,:)  = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2_n15(:,:)  / one_day * dtime 

    ! emission rates as positive flux to the atmosphere, integrated over the soil column [mol m-2 timestep-1]
    emission_noy(:)     = SUM(-transport_noy(:,:) * soil_depth_sl(:,:),DIM=2)
    emission_n2o(:)     = SUM(-transport_n2o(:,:) * soil_depth_sl(:,:),DIM=2)
    emission_n2(:)      = SUM(-transport_n2(:,:)  * soil_depth_sl(:,:),DIM=2)
    emission_noy_n15(:) = SUM(-transport_noy_n15(:,:) * soil_depth_sl(:,:),DIM=2)
    emission_n2o_n15(:) = SUM(-transport_n2o_n15(:,:) * soil_depth_sl(:,:),DIM=2)
    emission_n2_n15(:)  = SUM(-transport_n2_n15(:,:)  * soil_depth_sl(:,:),DIM=2)

  END SUBROUTINE calc_inorganic_transport_rate
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculate potential decomposition rate of SOM and litter pools according to 
  !!  dX/dt = - f(T)*g(M)/tau_pool * X
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_decomposition_rate(options, &                ! in
                                               flag_myco_priming, &
                                               rtm_decomposition_act, &
                                               rmm_decomposition_act, &
                                               sb_pool, &
                                               sb_loss )                ! out

    USE mo_sb_constants,            ONLY: tau_metabolic_litter, tau_structural_litter, tau_woody_litter, tau_fast, tau_slow
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options                  !< dimensions
    LOGICAL,                                          INTENT(in)    :: flag_myco_priming        !< from sb config
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_decomposition_act, & !< temperature response of decomposition
                                                                       rmm_decomposition_act    !< moisture response of decomposition 
    TYPE(t_sb_pool_a),                                INTENT(in)    :: sb_pool                  !< soil organic matter and litter pools [mol / m3]
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_loss                  !< soil organic matter and litter loss [mol / m3]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: k_priming, &
                                                        tau_slow_myc, &
                                                        hlp1
    REAL(wp)                                         :: vmax_priming
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_decomposition_rate'



    !> 1.0 metabolic litter (= soluable pool)
    !! 
    sb_loss%soluable_litter%carbon(:,:) = sb_loss%soluable_litter%carbon(:,:) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                          sb_pool%soluable_litter%carbon(:,:)
    sb_loss%soluable_litter%nitrogen(:,:) = sb_loss%soluable_litter%nitrogen(:,:) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                            sb_pool%soluable_litter%nitrogen(:,:)
    sb_loss%soluable_litter%phosphorus(:,:) = sb_loss%soluable_litter%phosphorus(:,:) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                              sb_pool%soluable_litter%phosphorus(:,:)
    sb_loss%soluable_litter%carbon13(:,:) = sb_loss%soluable_litter%carbon13(:,:) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                            sb_pool%soluable_litter%carbon13(:,:)
    sb_loss%soluable_litter%carbon14(:,:) = sb_loss%soluable_litter%carbon14(:,:) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                            sb_pool%soluable_litter%carbon14(:,:)
    sb_loss%soluable_litter%nitrogen15(:,:) = sb_loss%soluable_litter%nitrogen15(:,:) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                              sb_pool%soluable_litter%nitrogen15(:,:)

    !> 2.0 structural litter (= polymeric pool)
    !! 
    sb_loss%polymeric_litter%carbon(:,:) = sb_loss%polymeric_litter%carbon(:,:) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                          sb_pool%polymeric_litter%carbon(:,:)
    sb_loss%polymeric_litter%nitrogen(:,:) = sb_loss%polymeric_litter%nitrogen(:,:) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                            sb_pool%polymeric_litter%nitrogen(:,:)
    sb_loss%polymeric_litter%phosphorus(:,:) = sb_loss%polymeric_litter%phosphorus(:,:) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                              sb_pool%polymeric_litter%phosphorus(:,:)
    sb_loss%polymeric_litter%carbon13(:,:) = sb_loss%polymeric_litter%carbon13(:,:) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                            sb_pool%polymeric_litter%carbon13(:,:)
    sb_loss%polymeric_litter%carbon14(:,:) = sb_loss%polymeric_litter%carbon14(:,:) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                            sb_pool%polymeric_litter%carbon14(:,:)
    sb_loss%polymeric_litter%nitrogen15(:,:) = sb_loss%polymeric_litter%nitrogen15(:,:) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                              sb_pool%polymeric_litter%nitrogen15(:,:)

    !> 3.0 woody litter
    !! 
    sb_loss%woody_litter%carbon(:,:) = sb_loss%woody_litter%carbon(:,:) + &
                                                      rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                      sb_pool%woody_litter%carbon(:,:)
    sb_loss%woody_litter%nitrogen(:,:) = sb_loss%woody_litter%nitrogen(:,:) + &
                                                        rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%nitrogen(:,:)
    sb_loss%woody_litter%phosphorus(:,:) = sb_loss%woody_litter%phosphorus(:,:) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                          sb_pool%woody_litter%phosphorus(:,:)
    sb_loss%woody_litter%carbon13(:,:) = sb_loss%woody_litter%carbon13(:,:) + &
                                                        rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%carbon13(:,:)
    sb_loss%woody_litter%carbon14(:,:) = sb_loss%woody_litter%carbon14(:,:) + &
                                                        rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%carbon14(:,:)
    sb_loss%woody_litter%nitrogen15(:,:) = sb_loss%woody_litter%nitrogen15(:,:) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                          sb_pool%woody_litter%nitrogen15(:,:)

    !> 4.0 fast SOM pool (=microbial pool)
    !!
    sb_loss%microbial%carbon(:,:) = sb_loss%microbial%carbon(:,:) + & 
                                                   rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                   sb_pool%microbial%carbon(:,:)
    sb_loss%microbial%nitrogen(:,:) = sb_loss%microbial%nitrogen(:,:) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%nitrogen(:,:)
    sb_loss%microbial%phosphorus(:,:) = sb_loss%microbial%phosphorus(:,:) + & 
                                                       rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                       sb_pool%microbial%phosphorus(:,:)
    sb_loss%microbial%carbon13(:,:) = sb_loss%microbial%carbon13(:,:) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%carbon13(:,:)
    sb_loss%microbial%carbon14(:,:) = sb_loss%microbial%carbon14(:,:) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%carbon14(:,:)
    sb_loss%microbial%nitrogen15(:,:) = sb_loss%microbial%nitrogen15(:,:) + & 
                                                       rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                       sb_pool%microbial%nitrogen15(:,:)

    !> 5.0 slow SOM pool (=residue pool)
    !!
    sb_loss%residue%carbon(:,:) = sb_loss%residue%carbon(:,:) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/ &
                                                     tau_slow * dtime * sb_pool%residue%carbon(:,:)
    sb_loss%residue%nitrogen(:,:) = sb_loss%residue%nitrogen(:,:) + & 
                                                       rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/ &
                                                       tau_slow * dtime * sb_pool%residue%nitrogen(:,:)
    sb_loss%residue%phosphorus(:,:) = sb_loss%residue%phosphorus(:,:) + & 
                                                         rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/ &
                                                         tau_slow * dtime * sb_pool%residue%phosphorus(:,:)
    sb_loss%residue%carbon13(:,:) = sb_loss%residue%carbon13(:,:) + & 
                                                       rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/ &
                                                       tau_slow * dtime * sb_pool%residue%carbon13(:,:)
    sb_loss%residue%carbon14(:,:) = sb_loss%residue%carbon14(:,:) + & 
                                                       rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/ &
                                                       tau_slow * dtime * sb_pool%residue%carbon14(:,:)
    sb_loss%residue%nitrogen15(:,:) = sb_loss%residue%nitrogen15(:,:) + & 
                                                         rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/ &
                                                         tau_slow * dtime * sb_pool%residue%nitrogen15(:,:) 


  END SUBROUTINE calc_potential_decomposition_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculate potential biomineralisation rate 
  !!  dP/dt = - f(T)*g(M)/tau_pool * P
  !> Calculate weathering rate, slow adsorption rate and occlusion rate based on Yang et al. 2014
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_p_flux_rate(options, &                      ! in 
                                                   rtm_decomposition_act, &
                                                   rmm_decomposition_act, &
                                                   km_dip_biochem_po4_act, &
                                                   po4_solute, &
                                                   po4_assoc_fast, &
                                                   po4_assoc_slow, &
                                                   po4_occluded, &
                                                   po4_primary, &
                                                   k_weath_mineral_act, &
                                                   k_adsorpt_po4_act, &
                                                   k_desorpt_po4_act, &
                                                   sb_pool_residue_phosphorus, & 
                                                   weathering_po4, &
                                                   occlusion_po4, &
                                                   slow_sorption_po4, &
                                                   biochem_mineralisation_po4 )         ! out 

    USE mo_sb_constants,            ONLY: tau_biomineralisation, k_occlude_po4

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)  :: options                    !< dimensions
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: rtm_decomposition_act, &   !< temperature response of decomposition
                                                                     rmm_decomposition_act, &   !< moisture response of decomposition
                                                                     km_dip_biochem_po4_act, &  !< actual half-saturation soluble P for biochem_mineralisation_po4 [mol P / m3]
                                                                     po4_solute, &              !< soluble PO4 pool [mol/m3]
                                                                     po4_assoc_fast, &          !< fast minerally associated PO4 pool [mol/m3]
                                                                     po4_assoc_slow, &          !< slow minerally associated PO4 pool [mol/m3]
                                                                     po4_occluded, &            !< occluded PO4 pool [mol/m3]
                                                                     po4_primary, &             !< primary PO4 pool [mol/m3]
                                                                     k_weath_mineral_act, &     !< actual weathering rate for primary P [s-1]
                                                                     k_adsorpt_po4_act, &       !< actual adsorption rate for fast_assoc_po4 [s-1]
                                                                     k_desorpt_po4_act, &       !< actual desorption rate for slow_assoc_po4 [s-1]
                                                                     sb_pool_residue_phosphorus !< slow SOM P pool [mol / m3] 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out) :: weathering_po4, &          !< phosphorus weathering rate [mol P / m3 / timestep]
                                                                     occlusion_po4, &           !< phosphorus occlusion rate, assoc_slow --> occluded [mol P / m3 / timestep]
                                                                     slow_sorption_po4, &       !< phosphorus desorption rate, assoc_fast --> assoc_slow [mol P / m3 / timestep]
                                                                     biochem_mineralisation_po4 !< potential biomineralisation rate [mol / m3 / timestep]
    ! 0.2 Local
    ! 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1, hlp2
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_p_flux_rate'

    !> 1.0 Calculate P weathering rate from Primary P pool (mol /m3 / timestep)
    !!
    weathering_po4 = k_weath_mineral_act * po4_primary * dtime
    !> 1.1 Calculate P occlusion rate from occluded PO4 pool (mol /m3 / timestep)
    occlusion_po4 = po4_assoc_slow * k_occlude_po4 * dtime
    !> 1.2 Calculate P slow adsorption rate (mol /m3 / timestep)
    hlp1 = k_desorpt_po4_act * po4_assoc_slow * dtime
    hlp2 = k_adsorpt_po4_act * po4_assoc_fast * dtime
    slow_sorption_po4 = hlp2 - hlp1  

    !> 2.0 Calculate potential P biomineralisation from slow organic P pool (mol /m3 / timestep)
    !!
    hlp1 = km_dip_biochem_po4_act / (km_dip_biochem_po4_act + po4_solute)
    biochem_mineralisation_po4(:,:) = rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_biomineralisation * & 
                                 dtime * sb_pool_residue_phosphorus(:,:) * hlp1


  END SUBROUTINE calc_potential_p_flux_rate
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the immobilisation rate of NH4 and PO4 of microbial biomass, given
  !! litter decay and stoichiometry, the co-occurring rates of plant, mycorrhizal uptake,
  !! nitrification, denitrification and transport, biomineralisation of PO4
  !! such that the newly formed fast SOM has the prescribed C:N:P stoichiometry; the carbon-
  !! use-efficiency of microbes is within the observed bound and mass-balance constraints are 
  !! maintained. The result of this calculation are updates ("actual") rates of uptake, and processing
  !! as well as litter turnover rates.  
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_actual_flux_rates( options, &
                                     nh4_solute, &
                                     no3_solute, &
                                     po4_solute, &
                                     po4_assoc_fast, &
                                     transport_nh4_solute, &
                                     transport_no3_solute, &
                                     transport_po4_solute, &
                                     transport_po4_assoc_fast, &
                                     weathering_po4, &
                                     p_deposition_act, &
                                     slow_sorption_po4, &
                                     apparent_mic_cn , &
                                     km2_uptake_nh4_act, & 
                                     plant_uptake_nh4_sl, &
                                     plant_uptake_no3_sl, &
                                     plant_uptake_po4_sl, &
                                     mycorrhiza_uptake_nh4_sl, &
                                     mycorrhiza_uptake_no3_sl, &
                                     mycorrhiza_uptake_po4_sl, &
                                     nitrification_no3, &
                                     denitrification_n2, &
                                     biochem_mineralisation_po4, &
                                     sb_loss, &
                                     microbial_uptake_nh4_sl, &
                                     microbial_uptake_po4_sl)

    USE mo_sb_constants,            ONLY: frac_litter2fast_som, microbial_cn, microbial_np, &
                                          microbial_nue, microbial_pue
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options                       !< dimensions
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: nh4_solute, &                 !< NH4 in solution [mol / m3] 
                                                                       no3_solute, &                 !< NO3 in solution [mol / m3]
                                                                       po4_solute, &                 !< PO4 in solution [mol / m3]
                                                                       po4_assoc_fast, &             !< labile mineral P [mol / m3]
                                                                       transport_nh4_solute, &       !< transport of NH4 [mol / m3 / timestep]
                                                                       transport_no3_solute, &       !< transport of NO3 [mol / m3 / timestep]
                                                                       transport_po4_solute, &       !< transport of soluble PO4 [mol / m3 / timestep]
                                                                       transport_po4_assoc_fast, &   !< transport of labile PO4 [mol / m3 / timestep]
                                                                       weathering_po4, &             !< weathering of po4 [mol / m3 / timestep]
                                                                       slow_sorption_po4, &          !< po4 desorption rate, assoc_fast --> assoc_slow [mol P / m3 / timestep]
                                                                       p_deposition_act, &           !< actual P deposition rate, [mol / m3 / timestep]
                                                                       apparent_mic_cn, &            !< apparent microbial carbon-nitrogen ratio [mol/mol]
                                                                       km2_uptake_nh4_act            !< actual km2 for NH4 uptake [mol / m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: plant_uptake_nh4_sl, &        !< plant NH4 uptake [mol / m3 / timestep]
                                                                       plant_uptake_no3_sl, &        !< plant NO3 uptake [mol / m3 / timestep]
                                                                       plant_uptake_po4_sl, &        !< plant PO4 uptake [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_nh4_sl, &   !< mycorrhizal NH4 uptake [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_no3_sl, &   !< mycorrhizal NO3 uptake [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_po4_sl, &   !< mycorrhizal PO4 uptake [mol / m3 / timestep]
                                                                       nitrification_no3, &          !< nitrification rate [mol / m3 / timestep]
                                                                       denitrification_n2, &         !< denitrification rate [mol / m3 / timestep]
                                                                       biochem_mineralisation_po4    !< biomineralisation rate [mol / m3 / timestep]
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_loss                       !< decomposition of SOM and litter pools [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out)   :: microbial_uptake_nh4_sl, &    !< microbial immobilisation of NH4 [mol / m3 / timestep]
                                                                       microbial_uptake_po4_sl       !< microbial immobilisation of PO4 [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb)                :: hlp1, hlp2, transfer_c
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_actual_flux_rates'


    !> 1.0 formation of fast SOM from litter 
    !!

    !> 1.1 calculate potential transfer from litter decay to fast SOM pool, and associated potential microbial 
    !!     immobilisation rate as difference between maximal growth rate and N in litter
    transfer_c(:,:) = frac_litter2fast_som * & 
                      (sb_loss%soluable_litter%carbon(:,:) + & 
                       sb_loss%polymeric_litter%carbon(:,:)) 
    
    microbial_uptake_nh4_sl(:,:) = transfer_c(:,:) / apparent_mic_cn(:,:) - & 
                                   microbial_nue * (sb_loss%soluable_litter%nitrogen(:,:) + & 
                                                    sb_loss%polymeric_litter%nitrogen(:,:))
    !! limit uptake rate at low nitrogen assuming a Michaelis-Menten Kinetic 
    hlp1(:,:) = 1.25_wp * nh4_solute(:,:) / ( km2_uptake_nh4_act(:,:) + nh4_solute(:,:) )
    microbial_uptake_nh4_sl(:,:) = microbial_uptake_nh4_sl(:,:) * MIN(hlp1(:,:),1.0_wp)

    !> 1.2  Limit microbial, mycorrhizal and plant uptake, as well as nitrification to available NH4 minus transport losses
    !!
    WHERE(microbial_uptake_nh4_sl(:,:)>=0.0_wp)
       hlp1(:,:) = microbial_uptake_nh4_sl(:,:) + plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + mycorrhiza_uptake_nh4_sl(:,:)
       hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:)
       WHERE(hlp2(:,:) < hlp1(:,:)) 
           WHERE(hlp1(:,:) > eps8)
              microbial_uptake_nh4_sl(:,:)   = hlp2(:,:) * microbial_uptake_nh4_sl(:,:) / hlp1(:,:) 
              plant_uptake_nh4_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_sl(:,:) / hlp1(:,:)
              nitrification_no3(:,:)   = hlp2(:,:) * nitrification_no3(:,:) / hlp1(:,:)
              mycorrhiza_uptake_nh4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_nh4_sl(:,:) / hlp1(:,:)
           ELSEWHERE
              microbial_uptake_nh4_sl(:,:) = 0.0_wp
              plant_uptake_nh4_sl(:,:) = 0.0_wp 
              nitrification_no3(:,:) = 0.0_wp 
              mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
           ENDWHERE 
        ENDWHERE
    ELSEWHERE
       hlp1(:,:) = plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + mycorrhiza_uptake_nh4_sl(:,:)
       hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:) - microbial_uptake_nh4_sl(:,:)
       WHERE(hlp2(:,:) < hlp1(:,:)) 
           WHERE(hlp1(:,:) > eps8)
              plant_uptake_nh4_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_sl(:,:) / hlp1(:,:)
              nitrification_no3(:,:)   = hlp2(:,:) * nitrification_no3(:,:) / hlp1(:,:)
              mycorrhiza_uptake_nh4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_nh4_sl(:,:) / hlp1(:,:)
           ELSEWHERE
              plant_uptake_nh4_sl(:,:) = 0.0_wp 
              nitrification_no3(:,:) = 0.0_wp 
              mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
           ENDWHERE 
        ENDWHERE
    ENDWHERE 
   
    !> 1.3  limit plant uptake and denitrification to available NO3 minus transport losses 
    !!
    hlp1(:,:) = plant_uptake_no3_sl(:,:) + denitrification_n2(:,:) + mycorrhiza_uptake_no3_sl(:,:)
    hlp2(:,:) = no3_solute(:,:) + transport_no3_solute(:,:)
    WHERE(hlp2(:,:) < hlp1(:,:))
        WHERE(hlp1(:,:) > eps8)
           plant_uptake_no3_sl(:,:) = hlp2(:,:) * plant_uptake_no3_sl(:,:) / hlp1(:,:)
           denitrification_n2(:,:) = hlp2(:,:) * denitrification_n2(:,:) / hlp1(:,:)
           mycorrhiza_uptake_no3_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_no3_sl(:,:) / hlp1(:,:)
        ELSEWHERE
           plant_uptake_no3_sl(:,:) = 0.0_wp 
           denitrification_n2(:,:) = 0.0_wp 
           mycorrhiza_uptake_no3_sl(:,:) = 0.0_wp
        ENDWHERE
    ENDWHERE 

    !> 2.0 determine actual growth rate of fast SOM pool
    !!
   
    !> 2.1 Calculate actual litter loss rate, calculated such that litter decomposition
    !! matches the actual immobilisation flux and CN uptake from litter given actual CUE
    hlp1(:,:) = frac_litter2fast_som / apparent_mic_cn(:,:) * &
                      (sb_loss%soluable_litter%carbon(:,:) + & 
                       sb_loss%polymeric_litter%carbon(:,:)) - & 
                microbial_nue * & 
                       (sb_loss%soluable_litter%nitrogen(:,:) + & 
                        sb_loss%polymeric_litter%nitrogen(:,:)) 
    hlp2(:,:) = 1.0_wp
    WHERE(microbial_uptake_nh4_sl(:,:)<hlp1(:,:))
       WHERE(hlp1(:,:)>eps8)
           hlp2(:,:) = microbial_uptake_nh4_sl(:,:) / hlp1(:,:)
       ELSEWHERE
           hlp2(:,:) = 0.0_wp
       ENDWHERE
       sb_loss%soluable_litter%carbon(:,:) = hlp2(:,:) * sb_loss%soluable_litter%carbon(:,:)
       sb_loss%soluable_litter%nitrogen(:,:) = hlp2(:,:) * sb_loss%soluable_litter%nitrogen(:,:)
       sb_loss%soluable_litter%phosphorus(:,:) = hlp2(:,:) * sb_loss%soluable_litter%phosphorus(:,:)
       sb_loss%soluable_litter%carbon13(:,:) = hlp2(:,:) * sb_loss%soluable_litter%carbon13(:,:)
       sb_loss%soluable_litter%carbon14(:,:) = hlp2(:,:) * sb_loss%soluable_litter%carbon14(:,:)
       sb_loss%soluable_litter%nitrogen15(:,:) = hlp2(:,:) * sb_loss%soluable_litter%nitrogen15(:,:)

       sb_loss%polymeric_litter%carbon(:,:) = hlp2(:,:) * sb_loss%polymeric_litter%carbon(:,:)
       sb_loss%polymeric_litter%nitrogen(:,:) = hlp2(:,:) * sb_loss%polymeric_litter%nitrogen(:,:)
       sb_loss%polymeric_litter%phosphorus(:,:) = hlp2(:,:) * sb_loss%polymeric_litter%phosphorus(:,:)
       sb_loss%polymeric_litter%carbon13(:,:) = hlp2(:,:) * sb_loss%polymeric_litter%carbon13(:,:)
       sb_loss%polymeric_litter%carbon14(:,:) = hlp2(:,:) * sb_loss%polymeric_litter%carbon14(:,:)
       sb_loss%polymeric_litter%nitrogen15(:,:) = hlp2(:,:) * sb_loss%polymeric_litter%nitrogen15(:,:)
    ENDWHERE 

    !> 3.0 calculate potential microbial P uptake from soluable pool (mol /m3 / timestep) 
    microbial_uptake_po4_sl(:,:) = frac_litter2fast_som / apparent_mic_cn(:,:) / microbial_np * & 
                                    (sb_loss%soluable_litter%carbon(:,:) + &
                                     sb_loss%polymeric_litter%carbon(:,:)) - &
                                   microbial_pue * & 
                                    (sb_loss%soluable_litter%phosphorus(:,:) + & 
                                     sb_loss%polymeric_litter%phosphorus(:,:)) 

    !> 3.1 constrain P uptake by microbes and plants to the available sum of P from P_soluable, 
    !!     biomineralisation and transport (mol /m3 / timestep)   
    WHERE(microbial_uptake_po4_sl(:,:)>0.0_wp) ! if microbes take up P, constrain both microbial and plant uptake to available P 
      hlp1(:,:) = microbial_uptake_po4_sl(:,:) + plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:)
      hlp2(:,:) = po4_solute(:,:) + transport_po4_solute(:,:) 
      WHERE(hlp2(:,:) < hlp1(:,:))
          WHERE(hlp1(:,:) > eps8)
             microbial_uptake_po4_sl(:,:)  = hlp2(:,:) * microbial_uptake_po4_sl(:,:) / hlp1(:,:) 
             plant_uptake_po4_sl(:,:)      = hlp2(:,:) * plant_uptake_po4_sl(:,:) / hlp1(:,:) 
             mycorrhiza_uptake_po4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_po4_sl(:,:) / hlp1(:,:) 
          ELSEWHERE
             microbial_uptake_po4_sl(:,:)  = 0.0_wp
             plant_uptake_po4_sl(:,:)      = 0.0_wp
             mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp 
          ENDWHERE
       ENDWHERE
    ELSEWHERE ! if microbes release P, constrain plant uptake to available P, note the negative sign of microbial_uptake_po4_sl(:,:)
       hlp1(:,:) = plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) 
       hlp2(:,:) = po4_solute(:,:) + transport_po4_solute(:,:) - microbial_uptake_po4_sl(:,:)
       WHERE(hlp2(:,:) < hlp1(:,:))
          WHERE(hlp1(:,:) > eps8)
             plant_uptake_po4_sl(:,:)      = hlp2(:,:) * plant_uptake_po4_sl(:,:) / hlp1(:,:) 
             mycorrhiza_uptake_po4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_po4_sl(:,:) / hlp1(:,:) 
          ELSEWHERE
             plant_uptake_po4_sl(:,:)      = 0.0_wp 
             mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp 
          ENDWHERE
       ENDWHERE
    ENDWHERE

    !> 3.2 Limit biomineralisation to the actual rate of biomineralisation required to meet microbial and plant demands
    !!      and substract from slow=residue pool 
    sb_loss%residue%phosphorus(:,:) = sb_loss%residue%phosphorus(:,:) + & 
                                                     biochem_mineralisation_po4(:,:)
  
  END SUBROUTINE calc_actual_flux_rates
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to 
  !! net mineralisation \n
  !! For the sb_nloss_scheme = dynamic: calculates the partitioning of nitrification and denitrification
  !! to NO3, NOy, N2O and N2, and diagnoses the amount of leaching below the rooting zone
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_partitioning_nitrification_denitrification_rate(options, &
                                    rtm_denitrification_act, &
                                    nitrification_no3, &
                                    denitrification_n2, &
                                    nitrification_noy, &
                                    nitrification_n2o, &
                                    denitrification_noy, &
                                    denitrification_n2o)

    USE mo_sb_constants,            ONLY: f_nit_noy, f_nit_n2o, f_denit_noy, f_denit_n2o
  
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: rtm_denitrification_act    !< temperature modifier for denitrification
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: nitrification_no3, &       !< nitrification to NO3 [mol / m3 / timestep]
                                                                       denitrification_n2         !< denitrification to N2 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out)   :: nitrification_noy, &       !< nitrification to NOy [mol / m3 / timestep]
                                                                       nitrification_n2o, &       !< nitrification to N2O [mol / m3 / timestep]
                                                                       denitrification_noy, &     !< denitrification to NOy [mol / m3 / timestep]
                                                                       denitrification_n2o        !< denitrification to N2O [mol / m3 / timestep]
    ! 0.2 Local
    ! 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_partitioning_nitrification_denitrification_rate'
    
    !> 1.2 case of dynamic N loss calculations, following Zaehle et al. 2011
    !!
    ! nitrification rate to produce NO3, NOy and N2O from NH4 [mol m-3 timestep-1]
    nitrification_noy(:,:) = f_nit_noy * nitrification_no3(:,:)
    nitrification_n2o(:,:) = f_nit_n2o * nitrification_no3(:,:)
    nitrification_no3(:,:) = nitrification_no3(:,:) - nitrification_noy(:,:) - nitrification_n2o(:,:)
 
    ! denitrification rate to produce NOy, N2O and N2 from NO3 [mol m-3 timestep-1]
    denitrification_noy(:,:) = f_denit_noy * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
    denitrification_n2o(:,:) = f_denit_n2o * rtm_denitrification_act(:,:) * denitrification_n2(:,:)
    denitrification_n2(:,:)  = denitrification_n2(:,:) - denitrification_noy(:,:) - denitrification_n2o(:,:)

  END SUBROUTINE calc_partitioning_nitrification_denitrification_rate

  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the N15 fluxes based on uptake or transformation of the inorganic NH4 / NO3 pools
  !! according to the framework proposed by Robinson, 2001, TREE. Includes a representation of the 
  !! Rayleigh model to account for the change in fractionation with increasing consumption  
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sb_inorganic_n15_fluxes(options, & 
                                          nh4_solute, &
                                          nh4_n15_solute, &
                                          no3_solute, &
                                          no3_n15_solute, &
                                          plant_uptake_nh4_sl, & 
                                          plant_uptake_no3_sl, &      
                                          mycorrhiza_uptake_nh4_sl, & 
                                          mycorrhiza_uptake_no3_sl, & 
                                          microbial_uptake_nh4_sl, &        
                                          microbial_uptake_no3_sl, &        
                                          nitrification_no3, &        
                                          nitrification_noy, &        
                                          nitrification_n2o, &        
                                          denitrification_noy, &      
                                          denitrification_n2o, &      
                                          denitrification_n2, &       
                                          transport_nh4_n15_solute, & 
                                          transport_no3_n15_solute, &    
                                          plant_uptake_nh4_n15_sl, &      
                                          plant_uptake_no3_n15_sl, &      
                                          mycorrhiza_uptake_nh4_n15_sl, & 
                                          mycorrhiza_uptake_no3_n15_sl, & 
                                          microbial_uptake_nh4_n15_sl, &        
                                          microbial_uptake_no3_n15_sl, &        
                                          nitrification_no3_n15, &        
                                          nitrification_noy_n15, &        
                                          nitrification_n2o_n15, &        
                                          denitrification_noy_n15, &      
                                          denitrification_n2o_n15, &      
                                          denitrification_n2_n15) 
 
    USE mo_shared_isotope, ONLY: calc_fractionation 
    USE mo_sb_constants,   ONLY: eta_plant_uptake_nh4,eta_plant_uptake_no3,eta_mic_uptake_nh4,eta_mic_uptake_no3, & 
                                 eta_nitrate_production,eta_nitrification,eta_denitrification

    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options                     !< model options
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: nh4_solute, &               !< soil solute NH4 concentration [mol m-3]
                                                                       nh4_n15_solute, &           !< soil solute 15NH4 concentration [mol m-3]
                                                                       no3_solute, &               !< soil solute NO3 concentration [mol m-3]
                                                                       no3_n15_solute, &           !< soil solute 15NO3 concentration [mol m-3]
                                                                       plant_uptake_nh4_sl, &      !< NH4 uptake by plants [mol / m3 / timestep]
                                                                       plant_uptake_no3_sl, &      !< NO3 uptake by plants [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_nh4_sl, & !< NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_no3_sl, & !< NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                                       microbial_uptake_nh4_sl, &        !< NH4 uptake by microbes [mol / m3 / timestep]
                                                                       microbial_uptake_no3_sl, &        !< NH4 uptake by microbes [mol / m3 / timestep]
                                                                       nitrification_no3, &        !< nitrification to NO3 [mol / m3 / timestep]
                                                                       nitrification_noy, &        !< nitrification to NOy [mol / m3 / timestep]
                                                                       nitrification_n2o, &        !< nitrification to N2O [mol / m3 / timestep]
                                                                       denitrification_noy, &      !< denitrification to NOy [mol / m3 / timestep]
                                                                       denitrification_n2o, &      !< denitrification to N2O [mol / m3 / timestep]
                                                                       denitrification_n2, &       !< denitrification to N2 [mol / m3 / timestep]
                                                                       transport_nh4_n15_solute, & !< transport of 15NH4 [mol / m3 / timestep]
                                                                       transport_no3_n15_solute    !< transport of 15NO3 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out)   :: plant_uptake_nh4_n15_sl, &      !< 15NH4 uptake by plants [mol / m3 / timestep]
                                                                       plant_uptake_no3_n15_sl, &      !< 15NO3 uptake by plants [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_nh4_n15_sl, & !< 15NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_no3_n15_sl, & !< 15NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                                       microbial_uptake_nh4_n15_sl, &        !< 15NH4 uptake by microbes [mol / m3 / timestep]
                                                                       microbial_uptake_no3_n15_sl, &        !< 15NH4 uptake by microbes [mol / m3 / timestep]
                                                                       nitrification_no3_n15, &        !< nitrification to 15NO3 [mol / m3 / timestep]
                                                                       nitrification_noy_n15, &        !< nitrification to 15NOy [mol / m3 / timestep]
                                                                       nitrification_n2o_n15, &        !< nitrification to 15N2O [mol / m3 / timestep]
                                                                       denitrification_noy_n15, &      !< denitrification to 15NOy [mol / m3 / timestep]
                                                                       denitrification_n2o_n15, &      !< denitrification to 15N2O [mol / m3 / timestep]
                                                                       denitrification_n2_n15          !< denitrification to 15N2 [mol / m3 / timestep]
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1,hlp2 


    !> 2.0 calculate the actual discrimination for NH4, 
    !!     which depends on the fraction of the source consumed f = (sum of fluxes)/stock = hlp1. 
    !!     The discrimination is then reduced by a reduction factor r = - (1-f) ln (1-f) / f = hlp2
    ! fraction of NH4 that is consumed
    WHERE(nh4_solute(:,:)>eps8)
       hlp1(:,:) = MAX((plant_uptake_nh4_sl(:,:) + mycorrhiza_uptake_nh4_sl(:,:) + microbial_uptake_nh4_sl(:,:) + & 
                    nitrification_no3(:,:) + nitrification_noy(:,:) + nitrification_n2o(:,:)) / nh4_solute(:,:),eps8)
    ELSEWHERE
       hlp1(:,:) = 1.0_wp
    ENDWHERE

    ! calculate reduction factor for discrimination 
    WHERE(hlp1(:,:)<1.0_wp-eps8)
       hlp2(:,:) = (hlp1(:,:) - 1.0_wp) * log ( 1.0_wp - hlp1(:,:)) / hlp1(:,:)
    ELSEWHERE
       hlp2(:,:) = 0.0_wp
    ENDWHERE

    !> 1.1 15N fluxes involving ammonium
    !!
    plant_uptake_nh4_n15_sl(:,:)      = plant_uptake_nh4_sl(:,:) * & 
                                        calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_plant_uptake_nh4*hlp2(:,:)) 
    mycorrhiza_uptake_nh4_n15_sl(:,:) = mycorrhiza_uptake_nh4_sl(:,:) * & 
                                        calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_plant_uptake_nh4*hlp2(:,:)) 
    microbial_uptake_nh4_n15_sl(:,:)   = microbial_uptake_nh4_sl(:,:) * & 
                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_mic_uptake_nh4*hlp2(:,:))
    nitrification_no3_n15(:,:)   = nitrification_no3(:,:) * &
                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_nitrate_production*hlp2(:,:))
    nitrification_noy_n15(:,:)   = nitrification_noy(:,:) * & 
                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_nitrification*hlp2(:,:))
    nitrification_n2o_n15(:,:)   = nitrification_n2o(:,:) * &
                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_nitrification*hlp2(:,:)) 

    !> 1.2 contrain 15NH4 fluxes to available 15NH4
    !!
    hlp1(:,:) = plant_uptake_nh4_n15_sl(:,:) + mycorrhiza_uptake_nh4_n15_sl(:,:) + microbial_uptake_nh4_n15_sl(:,:) + & 
                nitrification_no3_n15(:,:) + nitrification_noy_n15(:,:) + nitrification_n2o_n15(:,:)
    hlp2(:,:) = nh4_n15_solute + transport_nh4_n15_solute(:,:) 
    WHERE(hlp2 < hlp1(:,:))
       WHERE(hlp1>eps8)
          plant_uptake_nh4_n15_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
          mycorrhiza_uptake_nh4_n15_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
          microbial_uptake_nh4_n15_sl(:,:)   = hlp2(:,:) * microbial_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
          nitrification_no3_n15(:,:)   = hlp2(:,:) * nitrification_no3_n15(:,:) / hlp1(:,:)
          nitrification_noy_n15(:,:)   = hlp2(:,:) * nitrification_noy_n15(:,:) / hlp1(:,:)
          nitrification_n2o_n15(:,:)   = hlp2(:,:) * nitrification_n2o_n15(:,:) / hlp1(:,:)
       ELSEWHERE
          plant_uptake_nh4_n15_sl(:,:) = 0.0_wp 
          mycorrhiza_uptake_nh4_n15_sl(:,:) = 0.0_wp 
          microbial_uptake_nh4_n15_sl(:,:)   = 0.0_wp
          nitrification_no3_n15(:,:)   = 0.0_wp
          nitrification_noy_n15(:,:)   = 0.0_wp
          nitrification_n2o_n15(:,:)   = 0.0_wp
       ENDWHERE
    ENDWHERE

    !> 2.0 calculate the actual discrimination for NO3, 
    !!     which depends on the fraction of the source consumed f = (sum of fluxes)/stock = hlp1. 
    !!     The discrimination is then reduced by a reduction factor r = - (1-f) ln (1-f) / f = hlp2
    ! fraction of NO3 that is consumed
    WHERE(no3_solute(:,:)>eps8)
       hlp1(:,:) = MAX((plant_uptake_no3_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) + microbial_uptake_no3_sl(:,:) + & 
                    denitrification_noy(:,:) + denitrification_n2o(:,:) + denitrification_n2(:,:)) / no3_solute(:,:),eps8)
    ELSEWHERE
       hlp1(:,:) = 1.0_wp
    ENDWHERE
    WHERE(hlp1(:,:)<1.0_wp-eps8)
       hlp2(:,:) = (hlp1(:,:) - 1.0_wp) * log ( 1.0_wp - hlp1(:,:)) / hlp1(:,:)
    ELSEWHERE
       hlp2(:,:) = 0.0_wp
    ENDWHERE

    !> 2.1 N15 fluxes involving NO3
    !!
    plant_uptake_no3_n15_sl(:,:)      = plant_uptake_no3_sl(:,:) * & 
                                        calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_plant_uptake_no3*hlp2(:,:))
    mycorrhiza_uptake_no3_n15_sl(:,:) = mycorrhiza_uptake_no3_sl(:,:) * & 
                                        calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_plant_uptake_no3*hlp2(:,:))
    microbial_uptake_no3_n15_sl(:,:)   = microbial_uptake_no3_sl(:,:) * & 
                                   calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_mic_uptake_no3*hlp2(:,:)) 
    denitrification_noy_n15(:,:) = denitrification_noy(:,:) * & 
                                   calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_denitrification*hlp2(:,:))
    denitrification_n2o_n15(:,:) = denitrification_n2o(:,:) * & 
                                   calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_denitrification*hlp2(:,:))
    denitrification_n2_n15(:,:)  = denitrification_n2(:,:) * & 
                                   calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_denitrification*hlp2(:,:)) 

    !> 2.2 contrain 15NO3 fluxes to available 15NO3
    !!
    hlp1(:,:) = plant_uptake_no3_n15_sl(:,:) + mycorrhiza_uptake_no3_n15_sl(:,:) + microbial_uptake_no3_n15_sl(:,:) + & 
                denitrification_noy_n15(:,:) + denitrification_n2o_n15(:,:) + denitrification_n2_n15(:,:)
    hlp2(:,:) = no3_n15_solute(:,:) + transport_no3_n15_solute(:,:)
    WHERE(hlp2 < hlp1(:,:))
       WHERE(hlp1>eps8)
          plant_uptake_no3_n15_sl(:,:) = hlp2(:,:) * plant_uptake_no3_n15_sl(:,:) / hlp1(:,:) 
          mycorrhiza_uptake_no3_n15_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_no3_n15_sl(:,:) / hlp1(:,:) 
          microbial_uptake_no3_n15_sl(:,:)   = hlp2(:,:) * microbial_uptake_no3_n15_sl(:,:) / hlp1(:,:) 
          denitrification_noy_n15(:,:)   = hlp2(:,:) * denitrification_noy_n15(:,:) / hlp1(:,:)
          denitrification_n2o_n15(:,:)   = hlp2(:,:) * denitrification_n2o_n15(:,:) / hlp1(:,:)
          denitrification_n2_n15(:,:)   = hlp2(:,:) * denitrification_n2_n15(:,:) / hlp1(:,:)
       ELSEWHERE
          plant_uptake_no3_n15_sl(:,:) = 0.0_wp 
          mycorrhiza_uptake_no3_n15_sl(:,:) = 0.0_wp 
          microbial_uptake_no3_n15_sl(:,:)   = 0.0_wp
          denitrification_noy_n15(:,:)   = 0.0_wp
          denitrification_n2o_n15(:,:)   = 0.0_wp
          denitrification_n2_n15(:,:)   = 0.0_wp
       ENDWHERE
    ENDWHERE

  END SUBROUTINE calc_sb_inorganic_n15_fluxes 


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the formation of SOM resulting from the decomposition dPrec/dt = CUE_Pdon * dPdon/dt
  !! gross mineralisation follows the stoichiometric difference of the donor and recipient pool, taking account for CUE
  !! for fast -> slow -> fast only causes gross mineralisation because of the C:N:P and CUE
  !! for litter -> fast gross mineralisation is typcially less than microbial uptake from 
  !!               mineral pools (or organic matter for P)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_SOM_formation(options, &                       ! in 
                                apparent_mic_cn, &
                                microbial_uptake_nh4_sl, &
                                microbial_uptake_nh4_n15_sl, &
                                microbial_uptake_po4_sl, &
                                sb_loss, &
                                het_respiration, &               ! out
                                het_respiration_c13, &
                                het_respiration_c14, &
                                net_mineralisation_nh4, &
                                net_mineralisation_po4, &
                                net_mineralisation_nh4_n15, &
                                sb_formation  )
 
    USE mo_shared_isotope,      ONLY: calc_fractionation
    USE mo_sb_constants,        ONLY: frac_litter2fast_som, frac_fast2slow_som, frac_slow2fast_som, microbial_np, microbial_cn, microbial_nue, microbial_pue, & 
                                      eta_ammonification
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options                       !< dimensions
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: apparent_mic_cn, &            !< apparent carbon-nitrogen ratio of microbial uptake [mol/mol] 
                                                                       microbial_uptake_nh4_sl, &          !< microbial NH4 uptake [mol / m3 / timestep]
                                                                       microbial_uptake_nh4_n15_sl, &      !< microbial 15NH4 uptake [mol / m3 / timestep]
                                                                       microbial_uptake_po4_sl             !< microbial PO4 uptake [mol / m3 / timestep]
    TYPE(t_sb_pool_a),                                INTENT(in)    :: sb_loss                       !< turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: het_respiration, &            !< heterotrophic respiration [mol C / m3 / timestep]
                                                                       het_respiration_c13, &        !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                                       het_respiration_c14, &        !< heterotrophic respiration [µmol 14C / m3 / timestep]
                                                                       net_mineralisation_nh4, &     !< net mineralisation of NH4 [mol / m3 / timestep]
                                                                       net_mineralisation_po4, &     !< net mineralisation of PO4 [mol / m3 / timestep]
                                                                       net_mineralisation_nh4_n15    !< net mineralisation of 15NH4 [mol / m3 / timestep]
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_formation                  !< formation of SOM and litter pools  [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1
    CHARACTER(len=*), PARAMETER                      :: routine = TRIM(modname)//':calc_SOM_formation'

    
    !> 1.0 formation of fast SOM from litter 
    !! 
    sb_formation%microbial%carbon(:,:) = sb_formation%microbial%carbon(:,:) + & 
                                                        frac_litter2fast_som * & 
                                                        (sb_loss%soluable_litter%carbon(:,:) + &  
                                                         sb_loss%polymeric_litter%carbon(:,:)) 
    sb_formation%microbial%nitrogen(:,:) = sb_formation%microbial%nitrogen(:,:) + & 
                                                          microbial_nue * &
                                                          (sb_loss%soluable_litter%nitrogen(:,:) + &  
                                                           sb_loss%polymeric_litter%nitrogen(:,:)) + & 
                                                          microbial_uptake_nh4_sl(:,:) 
    sb_formation%microbial%phosphorus(:,:) = sb_formation%microbial%phosphorus(:,:) + & 
                                                            microbial_pue * &
                                                            (sb_loss%soluable_litter%phosphorus(:,:) + &  
                                                             sb_loss%polymeric_litter%phosphorus(:,:)) + & 
                                                            microbial_uptake_po4_sl(:,:) 
    sb_formation%microbial%carbon13(:,:) = sb_formation%microbial%carbon13(:,:) + & 
                                                        frac_litter2fast_som * & 
                                                        (sb_loss%soluable_litter%carbon13(:,:) + &  
                                                         sb_loss%polymeric_litter%carbon13(:,:)) 
    sb_formation%microbial%carbon14(:,:) = sb_formation%microbial%carbon14(:,:) + & 
                                                        frac_litter2fast_som * & 
                                                        (sb_loss%soluable_litter%carbon14(:,:) + &  
                                                         sb_loss%polymeric_litter%carbon14(:,:))
    sb_formation%microbial%nitrogen15(:,:) = sb_formation%microbial%nitrogen15(:,:) + & 
                                                          microbial_nue * &
                                                          (sb_loss%soluable_litter%nitrogen15(:,:) + &  
                                                           sb_loss%polymeric_litter%nitrogen15(:,:)) + &
                                                           microbial_uptake_nh4_n15_sl(:,:) 

    het_respiration(:,:)     = het_respiration(:,:) + (1._wp - frac_litter2fast_som) * & 
                                                      (sb_loss%soluable_litter%carbon(:,:) + &  
                                                       sb_loss%polymeric_litter%carbon(:,:))
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - frac_litter2fast_som) * & 
                                                          (sb_loss%soluable_litter%carbon13(:,:) + &  
                                                           sb_loss%polymeric_litter%carbon13(:,:)) 
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - frac_litter2fast_som) * & 
                                                          (sb_loss%soluable_litter%carbon14(:,:) + &  
                                                           sb_loss%polymeric_litter%carbon14(:,:))
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) - microbial_uptake_nh4_sl(:,:) + &  
                                                                (1.0_wp - microbial_nue) * &
                                                                (sb_loss%soluable_litter%nitrogen(:,:) + &  
                                                                 sb_loss%polymeric_litter%nitrogen(:,:))
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) - microbial_uptake_po4_sl(:,:) + & 
                                                                (1.0_wp - microbial_pue) * &
                                                                (sb_loss%soluable_litter%phosphorus(:,:) + &  
                                                                 sb_loss%polymeric_litter%phosphorus(:,:))
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) - microbial_uptake_nh4_n15_sl(:,:) + & 
                                                                (1.0_wp - microbial_nue) * &
                                                                (sb_loss%soluable_litter%nitrogen15(:,:) + &  
                                                                 sb_loss%polymeric_litter%nitrogen15(:,:))


    !> 2.0 formation of slow SOM (=residue pool) from fast SOM decay
    !! 
    sb_formation%residue%carbon(:,:)     = sb_formation%residue%carbon(:,:) + & 
                                              frac_fast2slow_som * sb_loss%microbial%carbon(:,:)
    sb_formation%residue%nitrogen(:,:)   = sb_formation%residue%nitrogen(:,:) + & 
                                              frac_fast2slow_som / microbial_cn * sb_loss%microbial%carbon(:,:)
    sb_formation%residue%phosphorus(:,:) = sb_formation%residue%phosphorus(:,:) + & 
                                              frac_fast2slow_som / microbial_cn / microbial_np * sb_loss%microbial%carbon(:,:)
    sb_formation%residue%carbon13(:,:)   = sb_formation%residue%carbon13(:,:) + & 
                                              frac_fast2slow_som * sb_loss%microbial%carbon13(:,:)
    sb_formation%residue%carbon14(:,:)   = sb_formation%residue%carbon14(:,:) + & 
                                              frac_fast2slow_som * sb_loss%microbial%carbon14(:,:)
    hlp1(:,:) = (sb_loss%microbial%nitrogen(:,:) - frac_fast2slow_som / microbial_cn * sb_loss%microbial%carbon(:,:)) * & 
                calc_fractionation(sb_loss%microbial%nitrogen(:,:), &
                                   sb_loss%microbial%nitrogen15(:,:),eta_ammonification)
    sb_formation%residue%nitrogen15(:,:) = sb_formation%residue%nitrogen15(:,:) + &
                                             sb_loss%microbial%nitrogen15(:,:) - hlp1(:,:) 

    het_respiration(:,:)     = het_respiration(:,:)     + (1._wp - frac_fast2slow_som) * sb_loss%microbial%carbon(:,:)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - frac_fast2slow_som) * sb_loss%microbial%carbon13(:,:)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - frac_fast2slow_som) * sb_loss%microbial%carbon14(:,:)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb_loss%microbial%nitrogen(:,:) - &
                                                                frac_fast2slow_som / microbial_cn * sb_loss%microbial%carbon(:,:)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb_loss%microbial%phosphorus(:,:) - &
                                                                frac_fast2slow_som / microbial_cn / microbial_np * sb_loss%microbial%carbon(:,:)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + hlp1(:,:) 

    !> 3.0 formation of fast SOM (=microbial pool) from slow SOM decay 
    !! 
    sb_formation%microbial%carbon(:,:)     = sb_formation%microbial%carbon(:,:) + & 
                                                frac_slow2fast_som * sb_loss%residue%carbon(:,:)
    sb_formation%microbial%nitrogen(:,:)   = sb_formation%microbial%nitrogen(:,:) + & 
                                                frac_slow2fast_som / apparent_mic_cn(:,:) * sb_loss%residue%carbon(:,:)
    sb_formation%microbial%phosphorus(:,:) = sb_formation%microbial%phosphorus(:,:) + & 
                                                frac_slow2fast_som / apparent_mic_cn(:,:) / microbial_np * sb_loss%residue%carbon(:,:)
    sb_formation%microbial%carbon13(:,:) = sb_formation%microbial%carbon13(:,:) + & 
                                                frac_slow2fast_som * sb_loss%residue%carbon13(:,:)
    sb_formation%microbial%carbon14(:,:) = sb_formation%microbial%carbon14(:,:) + & 
                                                frac_slow2fast_som * sb_loss%residue%carbon14(:,:)
    hlp1(:,:) = (sb_loss%residue%nitrogen(:,:) - frac_slow2fast_som / apparent_mic_cn(:,:) * sb_loss%residue%carbon(:,:)) * & 
                calc_fractionation(sb_loss%residue%nitrogen(:,:), &
                                   sb_loss%residue%nitrogen15(:,:),eta_ammonification)
    sb_formation%microbial%nitrogen15(:,:) = sb_formation%microbial%nitrogen15(:,:) + &
                                                sb_loss%residue%nitrogen15(:,:) - hlp1(:,:) 

    het_respiration(:,:)      = het_respiration(:,:)     + (1._wp - frac_slow2fast_som) * sb_loss%residue%carbon(:,:)
    het_respiration_c13(:,:)  = het_respiration_c13(:,:) + (1._wp - frac_slow2fast_som) * sb_loss%residue%carbon13(:,:)
    het_respiration_c14(:,:)  = het_respiration_c14(:,:) + (1._wp - frac_slow2fast_som) * sb_loss%residue%carbon14(:,:)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb_loss%residue%nitrogen(:,:) - &
                                                                frac_slow2fast_som / apparent_mic_cn(:,:) * sb_loss%residue%carbon(:,:)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb_loss%residue%phosphorus(:,:) - &
                                                                frac_fast2slow_som / apparent_mic_cn(:,:) / microbial_np * sb_loss%residue%carbon(:,:)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + hlp1(:,:)

  END SUBROUTINE calc_SOM_formation  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the formation of metabolic and structural litter from woody litter decay. 
  !! Some C,N,P becomes lost during processing prescribed by the carbon and nutrient use
  !! efficiencies. The remainder becomes incorporated into metabolic and structural litter
  !! given the remaining litter C:N and the prescribed lignin content.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_litter_formation_from_woody_decay(options, &                       ! in 
                                                    sb_loss, &
                                                    het_respiration, &               ! out
                                                    het_respiration_c13, &
                                                    het_respiration_c14, &
                                                    net_mineralisation_nh4, &
                                                    net_mineralisation_po4, &
                                                    net_mineralisation_nh4_n15, &
                                                    sb_formation  )
  
    USE mo_sb_constants,        ONLY: frac_woody2dec_litter, microbial_nue, microbial_pue, & 
                                      fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, lc_woody_litter
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options                       !< dimensions
    TYPE(t_sb_pool_a),                                INTENT(in)    :: sb_loss                       !< turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: het_respiration, &            !< heterotrophic respiration [mol C / m3 / timestep]
                                                                       het_respiration_c13, &        !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                                       het_respiration_c14, &        !< heterotrophic respiration [µmol 14C / m3 / timestep]
                                                                       net_mineralisation_nh4, &     !< net mineralisation of NH4 [mol / m3 / timestep]
                                                                       net_mineralisation_po4, &     !< net mineralisation of PO4 [mol / m3 / timestep]
                                                                       net_mineralisation_nh4_n15    !< net mineralisation of 15NH4 [mol / m3 / timestep]
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_formation                  !< formation of SOM and litter pools  [mol / m3 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: fc_soluable_woody_litter, & !< fraction of non-respired woody litter becoming metabolic litter 
                                                        fn_soluable_woody_litter, & 
                                                        fp_soluable_woody_litter
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_litter_formation_from_woody_decay'

    
    !> 1.0 Calculate partitioning of decaying woody litter into metabolic and structural litter
    !!
    fc_soluable_woody_litter(:,:) = 0.0_wp
    fn_soluable_woody_litter(:,:) = 0.0_wp
    fp_soluable_woody_litter(:,:) = 0.0_wp
    WHERE(sb_loss%woody_litter%nitrogen(:,:) > eps8) 
       fc_soluable_woody_litter(:,:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_woody_litter * &
           sb_loss%woody_litter%carbon(:,:) / sb_loss%woody_litter%nitrogen(:,:))
       WHERE(fc_soluable_woody_litter(:,:)>eps8)
          fn_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:)) / &
                                          (k_fn_soluable*fc_soluable_woody_litter(:,:)))
          fp_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:)) / &
                                          (k_fp_soluable*fc_soluable_woody_litter(:,:)))
       ENDWHERE
    ENDWHERE

    !> 2.0 formation of metabolic and structural litter 
    !!     taking account of the processing losses of C, N and P
    sb_formation%soluable_litter%carbon(:,:) = sb_formation%soluable_litter%carbon(:,:) + & 
       fc_soluable_woody_litter(:,:) * frac_woody2dec_litter * sb_loss%woody_litter%carbon(:,:) 
    sb_formation%soluable_litter%nitrogen(:,:) = sb_formation%soluable_litter%nitrogen(:,:) + & 
       fn_soluable_woody_litter(:,:) * microbial_nue * sb_loss%woody_litter%nitrogen(:,:) 
    sb_formation%soluable_litter%phosphorus(:,:) = sb_formation%soluable_litter%phosphorus(:,:) + & 
       fp_soluable_woody_litter(:,:) * microbial_pue * sb_loss%woody_litter%phosphorus(:,:) 
    sb_formation%soluable_litter%carbon13(:,:) = sb_formation%soluable_litter%carbon13(:,:) + & 
       fc_soluable_woody_litter(:,:) * frac_woody2dec_litter * sb_loss%woody_litter%carbon13(:,:) 
    sb_formation%soluable_litter%carbon14(:,:) = sb_formation%soluable_litter%carbon14(:,:) + & 
       fc_soluable_woody_litter(:,:) * frac_woody2dec_litter * sb_loss%woody_litter%carbon14(:,:) 
    sb_formation%soluable_litter%nitrogen15(:,:) = sb_formation%soluable_litter%nitrogen15(:,:) + & 
       fn_soluable_woody_litter(:,:) * microbial_nue * sb_loss%woody_litter%nitrogen15(:,:) 

    sb_formation%polymeric_litter%carbon(:,:) = sb_formation%polymeric_litter%carbon(:,:) + & 
       (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter * sb_loss%woody_litter%carbon(:,:) 
    sb_formation%polymeric_litter%nitrogen(:,:) = sb_formation%polymeric_litter%nitrogen(:,:) + & 
       (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue * sb_loss%woody_litter%nitrogen(:,:) 
    sb_formation%polymeric_litter%phosphorus(:,:) = sb_formation%polymeric_litter%phosphorus(:,:) + & 
       (1.0_wp - fp_soluable_woody_litter(:,:)) * microbial_pue * sb_loss%woody_litter%phosphorus(:,:) 
    sb_formation%polymeric_litter%carbon13(:,:) = sb_formation%polymeric_litter%carbon13(:,:) + & 
       (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter * sb_loss%woody_litter%carbon13(:,:) 
    sb_formation%polymeric_litter%carbon14(:,:) = sb_formation%polymeric_litter%carbon14(:,:) + & 
       (1.0_wp - fc_soluable_woody_litter(:,:)) * frac_woody2dec_litter * sb_loss%woody_litter%carbon14(:,:) 
    sb_formation%polymeric_litter%nitrogen15(:,:) = sb_formation%polymeric_litter%nitrogen15(:,:) + & 
       (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue * sb_loss%woody_litter%nitrogen15(:,:) 

    !> 3.0 record losses to heterotrophic respiration and mineralisation
    !!
    het_respiration(:,:)            = het_respiration(:,:) + & 
                                      (1.0_wp - frac_woody2dec_litter) * sb_loss%woody_litter%carbon(:,:)
    net_mineralisation_nh4(:,:)     = net_mineralisation_nh4(:,:) + &
                                      (1.0_wp - microbial_nue) * sb_loss%woody_litter%nitrogen(:,:)
    net_mineralisation_po4(:,:)     = net_mineralisation_po4(:,:) + &
                                      (1.0_wp - microbial_pue) * sb_loss%woody_litter%phosphorus(:,:)
    het_respiration_c13(:,:)        = het_respiration_c13(:,:) + & 
                                      (1.0_wp - frac_woody2dec_litter) * sb_loss%woody_litter%carbon13(:,:)
    het_respiration_c14(:,:)        = het_respiration_c14(:,:) + & 
                                      (1.0_wp - frac_woody2dec_litter) * sb_loss%woody_litter%carbon14(:,:)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + &
                                      (1.0_wp - microbial_nue) * sb_loss%woody_litter%nitrogen15(:,:)
 
  END SUBROUTINE calc_litter_formation_from_woody_decay

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task called from update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> calculates matter transport by bioturbation, by calling calc_bioturbation_transport
  !! partly necessary because of a missing operator overload...
  !!
  !! Input:  pools subject to bioturbation, percolation rate, and soil depth 
  !!         This explicitly excludes mycorrhiza, as these die when eaten 
  !!
  !! Output: vertical transport rate (mol/m2/timestep) 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_bioturbation_transport_wrapper(options, k_bioturb, soil_depth_sl, &
                                                 pool_microbial, pool_residue, &
                                                 transport_microbial, transport_residue, & 
                                                 po4_assoc_fast, po4_assoc_slow, po4_occluded, &
                                                 transport_po4_assoc_fast, transport_po4_assoc_slow, & 
                                                 transport_po4_occluded)

    USE mo_sb_jsm_transport,        ONLY: calc_bioturbation_transport
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: k_bioturb                  !< diffusion factor for bioturbation 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: soil_depth_sl              !< depth of each soil layer [m] 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: po4_assoc_fast             !< fast minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: po4_assoc_slow             !< slow minerally associated PO4 pool [mol/m3]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: po4_occluded               !< occluded PO4 pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_microbial             !< microbial pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(in)    :: pool_residue               !< microbial residue pool [mol/m3]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_microbial        !< rate of microbial transport [mol/m3/timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: transport_residue          !< rate of microbial residue transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: transport_po4_assoc_fast   !< rate of fast associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: transport_po4_assoc_slow   !< rate of slow associated PO4 transport [mol/m3/timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: transport_po4_occluded     !< rate of occluded PO4 transport [mol/m3/timestep]
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_bioturbation_transport_wrapper'



    !> 1.0 fast SOM pool
    !!
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%carbon(:,:),transport_microbial%carbon(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%nitrogen(:,:),transport_microbial%nitrogen(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%phosphorus(:,:),transport_microbial%phosphorus(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%carbon13(:,:),transport_microbial%carbon13(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%carbon14(:,:),transport_microbial%carbon14(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%nitrogen15(:,:),transport_microbial%nitrogen15(:,:))

    !> 1.1 slow SOM pool
    !!
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_residue%carbon(:,:),transport_residue%carbon(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_residue%nitrogen(:,:),transport_residue%nitrogen(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_residue%phosphorus(:,:),transport_residue%phosphorus(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_residue%carbon13(:,:),transport_residue%carbon13(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_residue%carbon14(:,:),transport_residue%carbon14(:,:))
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        pool_residue%nitrogen15(:,:),transport_residue%nitrogen15(:,:))

    !> 2.0 inorganic pools
    !!
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        po4_assoc_fast,transport_po4_assoc_fast) 
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        po4_assoc_slow,transport_po4_assoc_slow) 
    CALL calc_bioturbation_transport(options,k_bioturb,soil_depth_sl, & 
                                        po4_occluded,transport_po4_occluded) 

  END SUBROUTINE calc_bioturbation_transport_wrapper

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> For the sb_nloss_scheme = fixed: calculates the amount of N lost in proportion to 
  !! net mineralisation \n
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_fixed_nloss_rate(options, &
                                   soil_depth_sl, &
                                   net_mineralisation_nh4, & 
                                   net_mineralisation_nh4_n15, & 
                                   nitrification_no3, &
                                   nitrification_no3_n15, &
                                   volatilisation_nh4, &
                                   volatilisation_nh4_n15, &
                                   emission_n2, &
                                   emission_n2_n15)

    USE mo_sb_constants,           ONLY: floss_nmin,fnit_nmin
  
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)  :: options
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)  :: soil_depth_sl, &              !< soil depth per layer [m]
                                                                     net_mineralisation_nh4, &     !< net mineralisation rate of NH4 [mol / m3 / timestep]
                                                                     net_mineralisation_nh4_n15    !< net mineralisation rate of 15NH4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out) :: nitrification_no3, &          !< nitrification to NO3 [mol / m3 / timestep]
                                                                     nitrification_no3_n15, &      !< nitrification to 15NO3 [mol / m3 / timestep]
                                                                     volatilisation_nh4, &         !< volatilisation of NH4 [mol / m3 / timestep]
                                                                     volatilisation_nh4_n15        !< volatilisation of 15NH4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc),                  INTENT(out) :: emission_n2, &                !< soil efflux of N2 [mol / m2 / timestep]
                                                                     emission_n2_n15               !< soil efflux of 15N2 [mol / m2 / timestep]
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_fixed_nloss_rate'


    !> 1.0 case of fixed N losses: N loss proportional to net mineralisation
    !!
    WHERE(net_mineralisation_nh4(:,:)>0.0_wp)
        ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
        volatilisation_nh4(:,:)      = net_mineralisation_nh4(:,:)     * floss_nmin 
        volatilisation_nh4_n15(:,:)  = net_mineralisation_nh4_n15(:,:) * floss_nmin 
        ! NH4 transformed to NO3 [mol m-3 timestep-1]
        nitrification_no3(:,:)     = net_mineralisation_nh4(:,:)     * (1._wp - floss_nmin) * fnit_nmin 
        nitrification_no3_n15(:,:) = net_mineralisation_nh4_n15(:,:) * (1._wp - floss_nmin) * fnit_nmin 
    ENDWHERE
    ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
    emission_n2(:)     = SUM(volatilisation_nh4(:,:)     * soil_depth_sl(:,:),DIM=2)
    emission_n2_n15(:) = SUM(volatilisation_nh4_n15(:,:) * soil_depth_sl(:,:),DIM=2)

  END SUBROUTINE calc_fixed_nloss_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the partitioning of litter fall from vegetation into metabolic, structural, 
  !! and woody litter given a) tissue specifc preferences (labile + reserve = metabolic; 
  !! sapwood and heartwood = woody); b) the tissue-specific lignin concentration and the actual
  !! tissue C:N ratio.
  !! All above-ground litter enters the first layer, only fine and coarse roots are distributed
  !! according to their root profile 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_litter_partitioning(options, &
                                      lctlib, & 
                                      soil_depth_sl, &
                                      root_fraction_act_sl, &
                                      veg_litterfall, &
                                      sb_formation )

    USE mo_veg_constants,           ONLY: itree, igrass
    USE mo_sb_constants,            ONLY: fc_soluable_max,k_fc_soluable,k_fn_soluable,k_fp_soluable, &
                                          lc_leaf_max,lc_leaf2sla,lc_fine_root,lc_coarse_root,lc_fruit,lc_seed_bed 
    dsl4slm_Use_memory_a_veg
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options
    TYPE(t_lctlib_element),                           INTENT(in)    :: lctlib                     !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: soil_depth_sl, &           !< depth per soil layer
                                                                       root_fraction_act_sl       !< root fraction per soil layer
    TYPE(t_veg_pool_a),                               INTENT(in)    :: veg_litterfall             !< vegetation litterfall pool
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_formation               !< soil biogeochemistry formation pool
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc) :: fc_soluable_leaf, fc_soluable_fine_root, fc_soluable_coarse_root, fc_soluable_fruit, & 
                                       fc_soluable_seed_bed, fc_soluable_labile, fc_soluable_reserve 
    REAL(wp), DIMENSION(options%nc) :: fn_soluable_leaf, fn_soluable_fine_root, fn_soluable_coarse_root, fn_soluable_fruit, & 
                                       fn_soluable_seed_bed, fn_soluable_labile, fn_soluable_reserve 
    REAL(wp), DIMENSION(options%nc) :: fp_soluable_leaf, fp_soluable_fine_root, fp_soluable_coarse_root, fp_soluable_fruit, & 
                                       fp_soluable_seed_bed, fp_soluable_labile, fp_soluable_reserve 
    INTEGER                     :: isoil                      !< loop over soil layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'



    !> 1.0 determine allocation of fresh litter fall to litter pools
    !!
    fc_soluable_leaf(:) = 0.0_wp
    fn_soluable_leaf(:) = 0.0_wp
    fp_soluable_leaf(:) = 0.0_wp
    WHERE(veg_litterfall%leaf%nitrogen(:) > eps8)
       fc_soluable_leaf(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib%sla) * & 
                                 veg_litterfall%leaf%carbon(:) / veg_litterfall%leaf%nitrogen(:))
       WHERE(fc_soluable_leaf(:)>eps8)
          fn_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fn_soluable*fc_soluable_leaf(:)))
          fp_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fp_soluable*fc_soluable_leaf(:)))
       END WHERE
    END WHERE

    fc_soluable_fine_root(:) = 0.0_wp 
    fn_soluable_fine_root(:) = 0.0_wp 
    fp_soluable_fine_root(:) = 0.0_wp 
    WHERE(veg_litterfall%fine_root%nitrogen(:) > eps8)
       fc_soluable_fine_root(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fine_root * & 
                                      veg_litterfall%fine_root%carbon(:) / veg_litterfall%fine_root%nitrogen(:))
       WHERE(fc_soluable_fine_root(:)>eps8)
          fn_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fn_soluable*fc_soluable_fine_root(:)))
          fp_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fp_soluable*fc_soluable_fine_root(:)))
       END WHERE
    END WHERE

    fc_soluable_coarse_root(:) = 0.0_wp 
    fn_soluable_coarse_root(:) = 0.0_wp 
    fp_soluable_coarse_root(:) = 0.0_wp 
    WHERE(veg_litterfall%coarse_root%nitrogen(:) > eps8)
       fc_soluable_coarse_root(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_coarse_root * & 
                                        veg_litterfall%coarse_root%carbon(:) / veg_litterfall%coarse_root%nitrogen(:))
       WHERE(fc_soluable_coarse_root(:)>eps8)
          fn_soluable_coarse_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(:))/(k_fn_soluable*fc_soluable_coarse_root(:)))
          fp_soluable_coarse_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(:))/(k_fp_soluable*fc_soluable_coarse_root(:)))
        END WHERE
    END WHERE

    fc_soluable_fruit(:) = 0.0_wp
    fn_soluable_fruit(:) = 0.0_wp
    fp_soluable_fruit(:) = 0.0_wp
    WHERE(veg_litterfall%fruit%nitrogen(:) > eps8)
       fc_soluable_fruit(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fruit * & 
                                  veg_litterfall%fruit%carbon(:) / veg_litterfall%fruit%nitrogen(:))
       WHERE(fc_soluable_fruit(:)>eps8)
          fn_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fn_soluable*fc_soluable_fruit(:)))
          fp_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fp_soluable*fc_soluable_fruit(:)))
       END WHERE
    END WHERE

    fc_soluable_seed_bed(:) = 0.0_wp
    fn_soluable_seed_bed(:) = 0.0_wp
    fp_soluable_seed_bed(:) = 0.0_wp
    WHERE(veg_litterfall%seed_bed%nitrogen(:) > eps8)
       fc_soluable_seed_bed(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_seed_bed * & 
                                     veg_litterfall%seed_bed%carbon(:) / veg_litterfall%seed_bed%nitrogen(:))
       WHERE(fc_soluable_seed_bed(:)>eps8)
          fn_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fn_soluable*fc_soluable_seed_bed(:)))
          fp_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fp_soluable*fc_soluable_seed_bed(:)))
       END WHERE
    END WHERE

    fc_soluable_labile(:)  = 1.0_wp 
    fn_soluable_labile(:)  = 1.0_wp 
    fp_soluable_labile(:)  = 1.0_wp 
    fc_soluable_reserve(:) = 1.0_wp
    fn_soluable_reserve(:) = 1.0_wp
    fp_soluable_reserve(:) = 1.0_wp

    !> 2.0 Update surface metabolic and soluable litter pools with litter fall
    !!
    sb_formation%soluable_litter%carbon(:,1) = sb_formation%soluable_litter%carbon(:,1) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon(:)    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon(:)   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon(:)  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon(:) + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%nitrogen(:,1) = sb_formation%soluable_litter%nitrogen(:,1) + & 
       (fn_soluable_leaf(:)     * veg_litterfall%leaf%nitrogen(:)    + & 
        fn_soluable_fruit(:)    * veg_litterfall%fruit%nitrogen(:)   + & 
        fn_soluable_labile(:)   * veg_litterfall%labile%nitrogen(:)  + & 
        fn_soluable_reserve(:)  * veg_litterfall%reserve%nitrogen(:) + & 
        fn_soluable_seed_bed(:) * veg_litterfall%seed_bed%nitrogen(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%phosphorus(:,1) = sb_formation%soluable_litter%phosphorus(:,1) + & 
       (fp_soluable_leaf(:)     * veg_litterfall%leaf%phosphorus(:)    + & 
        fp_soluable_fruit(:)    * veg_litterfall%fruit%phosphorus(:)   + & 
        fp_soluable_labile(:)   * veg_litterfall%labile%phosphorus(:)  + & 
        fp_soluable_reserve(:)  * veg_litterfall%reserve%phosphorus(:) + & 
        fp_soluable_seed_bed(:) * veg_litterfall%seed_bed%phosphorus(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%carbon13(:,1) = sb_formation%soluable_litter%carbon13(:,1) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon13(:)    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon13(:)   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon13(:)  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon13(:) + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon13(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%carbon14(:,1) = sb_formation%soluable_litter%carbon14(:,1) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon14(:)    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon14(:)   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon14(:)  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon14(:) + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon14(:) ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%nitrogen15(:,1) = sb_formation%soluable_litter%nitrogen15(:,1) + & 
       (fn_soluable_leaf(:)     * veg_litterfall%leaf%nitrogen15(:)    + & 
        fn_soluable_fruit(:)    * veg_litterfall%fruit%nitrogen15(:)   + & 
        fn_soluable_labile(:)   * veg_litterfall%labile%nitrogen15(:)  + & 
        fn_soluable_reserve(:)  * veg_litterfall%reserve%nitrogen15(:) + & 
        fn_soluable_seed_bed(:) * veg_litterfall%seed_bed%nitrogen15(:) ) / soil_depth_sl(:,1)

    sb_formation%polymeric_litter%carbon(:,1) = sb_formation%polymeric_litter%carbon(:,1) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon(:)    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon(:)   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon(:)  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon(:) + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen(:,1) = sb_formation%polymeric_litter%nitrogen(:,1) + & 
       ((1.0_wp - fn_soluable_leaf(:))     * veg_litterfall%leaf%nitrogen(:)    + & 
        (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall%fruit%nitrogen(:)   + & 
        (1.0_wp - fn_soluable_labile(:))   * veg_litterfall%labile%nitrogen(:)  + & 
        (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall%reserve%nitrogen(:) + & 
        (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall%seed_bed%nitrogen(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%phosphorus(:,1) = sb_formation%polymeric_litter%phosphorus(:,1) + & 
       ((1.0_wp - fp_soluable_leaf(:))     * veg_litterfall%leaf%phosphorus(:)    + & 
        (1.0_wp - fp_soluable_fruit(:))    * veg_litterfall%fruit%phosphorus(:)   + & 
        (1.0_wp - fp_soluable_labile(:))   * veg_litterfall%labile%phosphorus(:)  + & 
        (1.0_wp - fp_soluable_reserve(:))  * veg_litterfall%reserve%phosphorus(:) + & 
        (1.0_wp - fp_soluable_seed_bed(:)) * veg_litterfall%seed_bed%phosphorus(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon13(:,1) = sb_formation%polymeric_litter%carbon13(:,1) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon13(:)    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon13(:)   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon13(:)  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon13(:) + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon13(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon14(:,1) = sb_formation%polymeric_litter%carbon14(:,1) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon14(:)    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon14(:)   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon14(:)  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon14(:) + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon14(:) ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen15(:,1) = sb_formation%polymeric_litter%nitrogen15(:,1) + & 
       ((1.0_wp - fn_soluable_leaf(:))     * veg_litterfall%leaf%nitrogen15(:)    + & 
        (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall%fruit%nitrogen15(:)   + & 
        (1.0_wp - fn_soluable_labile(:))   * veg_litterfall%labile%nitrogen15(:)  + & 
        (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall%reserve%nitrogen15(:) + & 
        (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall%seed_bed%nitrogen15(:) ) / soil_depth_sl(:,1)

    !> 3.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
    !!     for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to 
    !!     metabolic and structural litter as if their were fine roots
    SELECT CASE(lctlib%growthform)
       CASE (itree)
          sb_formation%woody_litter%carbon(:,1) = sb_formation%woody_litter%carbon(:,1) + & 
             (veg_litterfall%sap_wood%carbon(:) + veg_litterfall%heart_wood%carbon) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%nitrogen(:,1) = sb_formation%woody_litter%nitrogen(:,1) +  & 
             (veg_litterfall%sap_wood%nitrogen(:) + veg_litterfall%heart_wood%nitrogen) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%phosphorus(:,1) = sb_formation%woody_litter%phosphorus(:,1) + & 
             (veg_litterfall%sap_wood%phosphorus(:) + veg_litterfall%heart_wood%phosphorus) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%carbon13(:,1) = sb_formation%woody_litter%carbon13(:,1) + & 
             (veg_litterfall%sap_wood%carbon13(:) + veg_litterfall%heart_wood%carbon13) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%carbon14(:,1) = sb_formation%woody_litter%carbon14(:,1) + & 
             (veg_litterfall%sap_wood%carbon14(:) + veg_litterfall%heart_wood%carbon14) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%nitrogen15(:,1) = sb_formation%woody_litter%nitrogen15(:,1) + &
             (veg_litterfall%sap_wood%nitrogen15 + veg_litterfall%heart_wood%nitrogen15) / &
             soil_depth_sl(:,1)
       CASE (igrass)
          sb_formation%soluable_litter%carbon(:,1) = sb_formation%soluable_litter%carbon(:,1) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon(:) + veg_litterfall%heart_wood%carbon) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%nitrogen(:,1) = sb_formation%soluable_litter%nitrogen(:,1) + & 
             fn_soluable_fine_root(:) * (veg_litterfall%sap_wood%nitrogen(:) + veg_litterfall%heart_wood%nitrogen) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%phosphorus(:,1) = sb_formation%soluable_litter%phosphorus(:,1) + & 
             fp_soluable_fine_root(:) * (veg_litterfall%sap_wood%phosphorus(:) + veg_litterfall%heart_wood%phosphorus) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%carbon13(:,1) = sb_formation%soluable_litter%carbon13(:,1) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon13(:) + veg_litterfall%heart_wood%carbon13) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%carbon14(:,1) = sb_formation%soluable_litter%carbon14(:,1) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon14(:) + veg_litterfall%heart_wood%carbon14) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%nitrogen15(:,1) = sb_formation%soluable_litter%nitrogen15(:,1) + & 
             fn_soluable_fine_root(:) * (veg_litterfall%sap_wood%nitrogen15 + veg_litterfall%heart_wood%nitrogen15) / &
             soil_depth_sl(:,1)
      
          sb_formation%polymeric_litter%carbon(:,1) = sb_formation%polymeric_litter%carbon(:,1) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * (veg_litterfall%sap_wood%carbon(:) + veg_litterfall%heart_wood%carbon) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%nitrogen(:,1) = sb_formation%polymeric_litter%nitrogen(:,1) + & 
             (1.0_wp - fn_soluable_fine_root(:)) * (veg_litterfall%sap_wood%nitrogen(:) + veg_litterfall%heart_wood%nitrogen) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%phosphorus(:,1) = sb_formation%polymeric_litter%phosphorus(:,1) + & 
             (1.0_wp - fp_soluable_fine_root(:)) * (veg_litterfall%sap_wood%phosphorus(:) + veg_litterfall%heart_wood%phosphorus) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%carbon13(:,1) = sb_formation%polymeric_litter%carbon13(:,1) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * (veg_litterfall%sap_wood%carbon13(:) + veg_litterfall%heart_wood%carbon13) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%carbon14(:,1) = sb_formation%polymeric_litter%carbon14(:,1) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * (veg_litterfall%sap_wood%carbon14(:) + veg_litterfall%heart_wood%carbon14) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%nitrogen15(:,1) = sb_formation%polymeric_litter%nitrogen15(:,1) + & 
             (1.0_wp - fn_soluable_fine_root(:)) * (veg_litterfall%sap_wood%nitrogen15 + veg_litterfall%heart_wood%nitrogen15) / &
             soil_depth_sl(:,1)

    END SELECT

    !> 4.0 Distribute fine and coarse root litter to soil depth according to the root profile 
    !!
    DO isoil = 1,options%nsoil_sb
       sb_formation%soluable_litter%carbon(:,isoil)     = sb_formation%soluable_litter%carbon(:,isoil) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon(:)     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%nitrogen(:,isoil)   = sb_formation%soluable_litter%nitrogen(:,isoil) + & 
          (fn_soluable_fine_root(:)   * veg_litterfall%fine_root%nitrogen(:)     + &
           fn_soluable_coarse_root(:) * veg_litterfall%coarse_root%nitrogen(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%phosphorus(:,isoil) = sb_formation%soluable_litter%phosphorus(:,isoil) + & 
          (fp_soluable_fine_root(:)   * veg_litterfall%fine_root%phosphorus(:)     + &
           fp_soluable_coarse_root(:) * veg_litterfall%coarse_root%phosphorus(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%carbon13(:,isoil)   = sb_formation%soluable_litter%carbon13(:,isoil) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon13(:)     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon13(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%carbon14(:,isoil)   = sb_formation%soluable_litter%carbon14(:,isoil) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon14(:)     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon14(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%nitrogen15(:,isoil) = sb_formation%soluable_litter%nitrogen15(:,isoil) + & 
          (fn_soluable_fine_root(:)   * veg_litterfall%fine_root%nitrogen15     + &
           fn_soluable_coarse_root(:) * veg_litterfall%coarse_root%nitrogen15 ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 

       sb_formation%polymeric_litter%carbon(:,isoil)    = sb_formation%polymeric_litter%carbon(:,isoil) + & 
          ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon(:)     + &
           (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%nitrogen(:,isoil)    = sb_formation%polymeric_litter%nitrogen(:,isoil) + & 
          ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen(:)     + &
           (1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%phosphorus(:,isoil)    = sb_formation%polymeric_litter%phosphorus(:,isoil) + & 
          ((1.0_wp - fp_soluable_fine_root(:))  * veg_litterfall%fine_root%phosphorus(:)     + &
           (1.0_wp - fp_soluable_coarse_root(:)) * veg_litterfall%coarse_root%phosphorus(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%carbon13(:,isoil)    = sb_formation%polymeric_litter%carbon13(:,isoil) + & 
          ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon13(:)     + &
           (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon13(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%carbon14(:,isoil)    = sb_formation%polymeric_litter%carbon14(:,isoil) + & 
          ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon14(:)     + &
           (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon14(:) ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%nitrogen15(:,isoil)    = sb_formation%polymeric_litter%nitrogen15(:,isoil) + & 
          ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen15     + &
           (1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen15 ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
    END DO  

  END SUBROUTINE calc_litter_partitioning


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Subroutine to deal with mycorrhiza 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_mycorrhiza_dynamics(options, &
                                      lctlib, &
                                      flag_myco_organic  , &
                                      soil_depth_sl, &
                                      mycorrhiza_uptake_nh4_sl, &
                                      mycorrhiza_uptake_no3_sl, &
                                      mycorrhiza_uptake_norg_sl, &
                                      mycorrhiza_uptake_po4_sl, &
                                      sb_pool, &
                                      sb_formation, &
                                      sb_loss, &
                                      het_respiration, &
                                      het_respiration_c13, &
                                      het_respiration_c14, &
                                      sb_mycorrhiza_export) 

    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: mycorrhiza_cn, mycorrhiza_np, tau_slow
    dsl4slm_Use_memory_a_sb

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options
    TYPE(t_lctlib_element),                           INTENT(in)    :: lctlib                       !< land-cover-type library - parameter across pft's
    LOGICAL,                                          INTENT(in)    :: flag_myco_organic            !< from sb config
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: soil_depth_sl                !< depth per soil layer
    TYPE(t_sb_pool_a),                                INTENT(in)    :: sb_pool                      !< soil biogeochemistry pools [mol / m3]
    TYPE(t_sb_pool_a),                                INTENT(inout) :: sb_formation, &              !< pool formation rate  [mol / m3 / timestep]
                                                                       sb_loss                      !< pool loss rate  [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: mycorrhiza_uptake_nh4_sl, &  !< NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_no3_sl, &  !< NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_norg_sl, & !< N uptake by mycorrhiza from organic (residue) pool [mol / m3 / timestep]
                                                                       mycorrhiza_uptake_po4_sl     !< PO4 uptake by mycorrhiza [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(out)   :: het_respiration, &           !< heterotrophic respiration [mol / m3 / timestep]
                                                                       het_respiration_c13, &       !< heterotrophic respiration [mol / m3 / timestep]
                                                                       het_respiration_c14          !< heterotrophic respiration [mol / m3 / timestep]
    TYPE(t_sb_element_a),                             INTENT(inout) :: sb_mycorrhiza_export         !< export from mycorrhiza to plants [mol / m2 / timestep]
    ! ---------------------------
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: transform_respiration, &
                                                        fresp, &
                                                        fturn_mycorrhiza, &
                                                        vmax_norg_up, &
                                                        k_norg_up, &
                                                        fnorg, &
                                                        mycorrhiza_cue, &
                                                        f_myc_opt
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb) :: hlp1, hlp2, hlp3
    INTEGER                     :: isoil                                            !< loop over soil layers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_mycorrhiza_dynamics'


    
    !> 1.0 Mycorrhiza dynamics 
    !! 

    !> 1.1 uptake of nutrients by mycorrhiza and associated respiration losses to heterotrophic respiration
    !!
    WHERE(sb_pool%mycorrhiza%nitrogen(:,:) > eps8 .AND. &
          sb_pool%mycorrhiza%carbon(:,:) / sb_pool%mycorrhiza%nitrogen(:,:) < mycorrhiza_cn)
        mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
        mycorrhiza_uptake_no3_sl(:,:) = 0.0_wp
        mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp
    ENDWHERE

    ! calculate respiration costs associated with N transformation (mol C / m3 / timestep)
    transform_respiration(:,:) = mycorrhiza_uptake_nh4_sl(:,:) * transform_cost_nh4 + & 
                                 mycorrhiza_uptake_no3_sl(:,:) * transform_cost_no3

    fresp(:,:) = 0.0_wp
    WHERE(sb_pool%mycorrhiza%carbon(:,:) > eps8)
       fresp(:,:) = transform_respiration(:,:) / sb_pool%mycorrhiza%carbon(:,:)
    ENDWHERE

    sb_loss%mycorrhiza%carbon(:,:)   = sb_loss%mycorrhiza%carbon(:,:)   + fresp * sb_pool%mycorrhiza%carbon(:,:)   !! respiration
    sb_loss%mycorrhiza%carbon13(:,:) = sb_loss%mycorrhiza%carbon13(:,:) + fresp * sb_pool%mycorrhiza%carbon13(:,:) !! respiration
    sb_loss%mycorrhiza%carbon14(:,:) = sb_loss%mycorrhiza%carbon14(:,:) + fresp * sb_pool%mycorrhiza%carbon14(:,:) !! respiration

    het_respiration(:,:)     = het_respiration(:,:)     + fresp * sb_pool%mycorrhiza%carbon(:,:)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + fresp * sb_pool%mycorrhiza%carbon13(:,:)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + fresp * sb_pool%mycorrhiza%carbon14(:,:)

    !> 1.2 organic N uptake by mycorrhiza and associated respiration losses to heterotrophic respiration
    !!
    mycorrhiza_uptake_norg_sl(:,:) = 0.0_wp

    fnorg(:,:) = 0.0_wp
    WHERE(sb_pool%residue%nitrogen(:,:) > eps8)
        fnorg(:,:) = mycorrhiza_uptake_norg_sl(:,:) / sb_pool%residue%nitrogen(:,:)
    ENDWHERE

    sb_loss%residue%carbon(:,:)     = sb_loss%residue%carbon(:,:)     + fnorg * sb_pool%residue%carbon(:,:)
    sb_loss%residue%nitrogen(:,:)   = sb_loss%residue%nitrogen(:,:)   + fnorg * sb_pool%residue%nitrogen(:,:)
    sb_loss%residue%phosphorus(:,:) = sb_loss%residue%phosphorus(:,:) + fnorg * sb_pool%residue%phosphorus(:,:)
    sb_loss%residue%carbon13(:,:)   = sb_loss%residue%carbon13(:,:)   + fnorg * sb_pool%residue%carbon13(:,:)
    sb_loss%residue%carbon14(:,:)   = sb_loss%residue%carbon14(:,:)   + fnorg * sb_pool%residue%carbon14(:,:)
    sb_loss%residue%nitrogen15(:,:) = sb_loss%residue%nitrogen15(:,:) + fnorg * sb_pool%residue%nitrogen15(:,:)

    !> 1.3 nutrient flux to vegetation (surplus uptake goes to plants) (mol (N/P) / m3 / timestep)
    !!
!@todo: N15 uptake
 
    ! N demand of mycorrhiza in mol N / m3 / timestep
    hlp1(:,:) = sb_pool%mycorrhiza%carbon(:,:) / mycorrhiza_cn - sb_pool%mycorrhiza%nitrogen(:,:)

    ! constrain flux to minimum of potential uptake and demand
    hlp2(:,:) = MIN(hlp1(:,:), &
                mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) + mycorrhiza_uptake_norg_sl(:,:))

    ! maintain C:N of mycorrhizal fungi
    ! N supply for plants = N uptake by mycorrhiza - N demand by mycorrhiza
    ! N demand can be negative => N supply can be greater then current N uptake by plants!

    ! remainder is exported to plants (mol N / m3 / timestep)
    hlp3(:,:) = mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) + mycorrhiza_uptake_norg_sl(:,:) - hlp2(:,:)

    sb_mycorrhiza_export%nitrogen(:,:) = 0.0_wp
    sb_mycorrhiza_export%nitrogen(:,:) = sb_mycorrhiza_export%nitrogen(:,:) + hlp3(:,:)

    ! P demand of mycorrhiza in mol P / m3 / timestep ! BLARPP SZ: constrained to positive solutions
    hlp1(:,:) = MAX(sb_pool%mycorrhiza%carbon(:,:) / mycorrhiza_cn / mycorrhiza_np - &
                    sb_pool%mycorrhiza%phosphorus(:,:), 0.0_wp)
    ! incorporation of P into mycorrhiza (mol P / m3 / timestep)
    ! constrain flux to minimum of potential uptake and demand 
    hlp2(:,:) = MIN(hlp1(:,:), mycorrhiza_uptake_po4_sl(:,:))

    ! remainder is exported to plants (mol P / m3 / timestep)
    hlp3(:,:) = mycorrhiza_uptake_po4_sl(:,:) - hlp2(:,:)
    sb_mycorrhiza_export%phosphorus(:,:) = 0.0_wp
    sb_mycorrhiza_export%phosphorus(:,:) = sb_mycorrhiza_export%phosphorus(:,:) + hlp3(:,:) 

    !new mycorrhiza formation out of uptake from slow N pool 
    mycorrhiza_cue = 0.5_wp

    sb_formation%mycorrhiza%carbon(:,:)     = sb_formation%mycorrhiza%carbon(:,:) + &
                                                             fnorg * sb_pool%residue%carbon(:,:) * mycorrhiza_cue
    sb_formation%mycorrhiza%nitrogen(:,:)   = sb_formation%mycorrhiza%nitrogen(:,:) + &
                                                             fnorg * sb_pool%residue%nitrogen(:,:)
    sb_formation%mycorrhiza%phosphorus(:,:) = sb_formation%mycorrhiza%phosphorus(:,:) + &
                                                             fnorg * sb_pool%residue%phosphorus(:,:)
    sb_formation%mycorrhiza%carbon13(:,:)   = sb_formation%mycorrhiza%carbon13(:,:) + &
                                                             fnorg * sb_pool%residue%carbon13(:,:) * mycorrhiza_cue
    sb_formation%mycorrhiza%carbon14(:,:)   = sb_formation%mycorrhiza%carbon14(:,:) + &
                                                             fnorg * sb_pool%residue%carbon14(:,:) * mycorrhiza_cue
    sb_formation%mycorrhiza%nitrogen15(:,:) = sb_formation%mycorrhiza%nitrogen15(:,:) + &
                                                             fnorg * sb_pool%residue%nitrogen15(:,:)

    het_respiration(:,:)     = het_respiration(:,:)     + fnorg * sb_pool%residue%carbon(:,:)   * (1.0_wp - mycorrhiza_cue)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + fnorg * sb_pool%residue%carbon13(:,:) * (1.0_wp - mycorrhiza_cue)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + fnorg * sb_pool%residue%carbon14(:,:) * (1.0_wp - mycorrhiza_cue)

    !> 1.4 Formation of litter from mycorrhiza turnover (mol / m3 / timestep)
    !!
    fturn_mycorrhiza(:,:) = 1.0_wp/lctlib%tau_mycorrhiza/one_day/one_year * dtime

    sb_loss%mycorrhiza%carbon(:,:)   = sb_loss%mycorrhiza%carbon(:,:) + &
                                                      fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon(:,:)     !! turnover
    sb_loss%mycorrhiza%nitrogen(:,:) = sb_loss%mycorrhiza%nitrogen(:,:) + &
                                                      fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen(:,:)   !! turnover
    sb_loss%mycorrhiza%phosphorus(:,:) = sb_loss%mycorrhiza%phosphorus(:,:) + &
                                                        fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%phosphorus(:,:) !! turnover
    sb_loss%mycorrhiza%carbon13(:,:) = sb_loss%mycorrhiza%carbon13(:,:) + &
                                                      fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon13(:,:)   !! turnover
    sb_loss%mycorrhiza%carbon14(:,:) = sb_loss%mycorrhiza%carbon14(:,:) + &
                                                        fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon14(:,:)   !! turnover
    sb_loss%mycorrhiza%nitrogen15(:,:) = sb_loss%mycorrhiza%nitrogen15(:,:) + &
                                                        fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen15(:,:) !! turnover

    sb_formation%polymeric_litter%carbon(:,:) = sb_formation%polymeric_litter%carbon(:,:) + &
                                                               fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon(:,:)
    sb_formation%polymeric_litter%carbon13(:,:) = sb_formation%polymeric_litter%carbon13(:,:) + &
                                                                 fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon13(:,:)
    sb_formation%polymeric_litter%carbon14(:,:) = sb_formation%polymeric_litter%carbon14(:,:) + &
                                                                 fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon14(:,:)
    sb_formation%polymeric_litter%nitrogen(:,:) = sb_formation%polymeric_litter%nitrogen(:,:) + &
                                                                 fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen(:,:)
    sb_formation%polymeric_litter%nitrogen15(:,:) = sb_formation%polymeric_litter%nitrogen15(:,:) + &
                                                                   fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen15(:,:)
    sb_formation%polymeric_litter%phosphorus(:,:) = sb_formation%polymeric_litter%phosphorus(:,:) + &
                                                                   fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%phosphorus(:,:)

  END SUBROUTINE calc_mycorrhiza_dynamics

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the rate of asymbiontic N fixation, wherever there is a bit of soil carbon 
  !! declines within the soil profile according to the distribution of roots as a proxy for litter input
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_asymbiontic_n_fixation(options, &
                                         lctlib, &                                      ! in
                                         veg_bnf_scheme, &
                                         soil_depth_sl, &
                                         nh4_solute, &
                                         no3_solute, &
                                         sb_pool_microbial_carbon, &
                                         root_fraction_act_sl, & 
                                         asymb_n_fixation, &                                ! out
                                         asymb_n_fixation_n15)                              ! out

    USE mo_shared_isotope,     ONLY: calc_mixing_ratio_N15N14
    USE mo_sb_constants,       ONLY: soil_n_threshold_bnf
    USe mo_veg_constants,      ONLY: eta_nfixation

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                         INTENT(in)    :: options
    TYPE(t_lctlib_element),                           INTENT(in)    :: lctlib                      !< land-cover-type library - parameter across pft's
    CHARACTER(len=*),                                 INTENT(in)    :: veg_bnf_scheme              !< veg n fixation scheme  
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(in)    :: soil_depth_sl, &            !< soil depth per layer [m] 
                                                                       nh4_solute, &               !< NH4 concentration [mol m-3]
                                                                       no3_solute, &               !< NO3 concentration [mol m-3]
                                                                       sb_pool_microbial_carbon, & !< microbial carbon pool [mol m-3]
                                                                       root_fraction_act_sl        !< fraction of roots per root layer [m-3] 
    REAL(wp), DIMENSION(options%nc,options%nsoil_sb), INTENT(inout) :: asymb_n_fixation, &         !< asymbiontic rate of N fixation [umol m-3 s-1]
                                                                       asymb_n_fixation_n15        !< asymbiontic rate of 15N fixation [umol m-3 s-1]
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_asymbiontic_n_fixation'

 
    !> 1.0 set asymbiontic N fixation to value in lctlib for the fixed and standard schemes
    !!
    SELECT CASE(TRIM(veg_bnf_scheme))
      CASE("fixed", "standard")
        WHERE((nh4_solute(:,:) + no3_solute) * soil_depth_sl(:,:) < soil_n_threshold_bnf .AND. & 
               sb_pool_microbial_carbon(:,:) > eps4 .AND. soil_depth_sl(:,:)>eps8)
                 asymb_n_fixation(:,:) = root_fraction_act_sl(:,:) * lctlib%bnf_base / soil_depth_sl(:,:) 
        ENDWHERE
    END SELECT

    !> 1.1 calculate 15N flux associated with asymbiontic N fixation, which is identical to the 
    !!     atmospheric dN15 minus the discrimination of the nitrogenase (Robinson, 2001, TREE)
    !! 
    asymb_n_fixation_n15(:,:) = asymb_n_fixation(:,:) / ( 1._wp + 1._wp / calc_mixing_ratio_N15N14(-eta_nfixation))

  END SUBROUTINE calc_asymbiontic_n_fixation

END MODULE mo_sb_simple_model
