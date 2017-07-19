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
  USE mo_exception,             ONLY: message, message_text, finish
  USE mo_jsb_class,             ONLY: Get_model
  USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
  USE mo_jsb_model_class,       ONLY: t_jsb_model
  USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element    
  USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
  USE mo_base_constants,        ONLY: dtime, one_day, one_year
  USE mo_sb_memory_class,       ONLY: t_slm_struct_sb_pool,t_slm_struct_sb_element

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

    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_veg_constants,         ONLY: itree
    USE mo_sb_jsm_transport,      ONLY: calc_bioturbation_rate

    ! Use of process configurations (t_PROC_config)
    dsl4jsb_Use_config(SB_)
    dsl4jsb_Use_config(VEG_)

    ! Use of process memories
    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET    :: tile          !< one tile
    TYPE(t_jsb_task_options), INTENT(in)               :: options       !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
    INTEGER                               :: isoil                !< loop over soil_config%nsoil_water
    REAL(wp), DIMENSION(options%nc,soil_config%nsoil_water) ::      &
                                             apparent_mic_cue     , &
                                             mic_uptake_nh4_sl    , &
                                             mic_uptake_po4_sl    , &
                                             biomineralisation_po4, &
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
                                             anaerobic_volume_fraction_act
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_simple_model'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(SB_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    ! VEG_
    dsl4jsb_Real2D_onChunk :: unit_uptake_n_pot
    dsl4jsb_Real2D_onChunk :: unit_uptake_p_pot
    ! SOIL_ 2D
    dsl4jsb_Real2D_onChunk :: drainage_fraction 
    ! SOIL_ 3D
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
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
    dsl4jsb_Real3D_onChunk :: rmm_decomposition 
    dsl4jsb_Real3D_onChunk :: rmm_hsc
    dsl4jsb_Real3D_onChunk :: rmm_nitrification
    dsl4jsb_Real3D_onChunk :: rmm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: anaerobic_volume_fraction_sl 
    dsl4jsb_Real3D_onChunk :: het_respiration
    dsl4jsb_Real3D_onChunk :: het_respiration_c13
    dsl4jsb_Real3D_onChunk :: het_respiration_c14
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
    dsl4jsb_Real3D_onChunk :: net_mineralisation_po4
    dsl4jsb_Real3D_onChunk :: nitrification_no3 
    dsl4jsb_Real3D_onChunk :: nitrification_noy 
    dsl4jsb_Real3D_onChunk :: nitrification_n2o
    dsl4jsb_Real3D_onChunk :: denitrification_noy 
    dsl4jsb_Real3D_onChunk :: denitrification_n2o 
    dsl4jsb_Real3D_onChunk :: denitrification_n2 
    dsl4jsb_Real3D_onChunk :: volatilisation_n
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_nh4_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_no3_n15_sl
    dsl4jsb_Real3D_onChunk :: mycorrhiza_uptake_po4_sl
    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
    dsl4jsb_Real3D_onChunk :: transport_no3_solute
    dsl4jsb_Real3D_onChunk :: transport_po4_solute
    dsl4jsb_Real3D_onChunk :: transport_noy
    dsl4jsb_Real3D_onChunk :: transport_n2o
    dsl4jsb_Real3D_onChunk :: transport_n2 
    dsl4jsb_Real3D_onChunk :: noy
    dsl4jsb_Real3D_onChunk :: n2o
    dsl4jsb_Real3D_onChunk :: n2 
    dsl4jsb_Real3D_onChunk :: nh4_solute 
    dsl4jsb_Real3D_onChunk :: no3_solute
    dsl4jsb_Real3D_onChunk :: po4_solute 
    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
    dsl4jsb_Real3D_onChunk :: k_bioturb 
    dsl4jsb_Real3D_onChunk :: bulk_dens_corr_sl 
    dsl4jsb_Real3D_onChunk :: asymb_n_fixation 
!! *slmtemp    
    TYPE(t_sb_memory), POINTER              :: sb         ! local pointer to sb__mem for SLM usability
    TYPE(t_veg_memory), POINTER             :: veg       ! local pointer to soil__mem for SLM usability
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
    dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SB_)
    ! Get process variables (Set pointers to variables in memory)
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_n_pot)
    dsl4jsb_Get_var2D_onChunk(VEG_,      unit_uptake_p_pot)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SOIL_,     drainage_fraction)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     soil_depth_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     root_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute) 
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_noy)          
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2o)         
    dsl4jsb_Get_var2D_onChunk(SB_,       emission_n2)         
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_decomposition)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_hsc)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       anaerobic_volume_fraction_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c13)
    dsl4jsb_Get_var3D_onChunk(SB_,       het_respiration_c14)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_po4)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)
    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)
    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_n)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_nh4_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_no3_n15_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       mycorrhiza_uptake_po4_sl)
    dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o) 
    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       po4_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       nh4_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       no3_n15_solute) 
    dsl4jsb_Get_var3D_onChunk(SB_,       k_bioturb)                   ! out
    dsl4jsb_Get_var3D_onChunk(SB_,       bulk_dens_corr_sl)           ! inout
    dsl4jsb_Get_var3D_onChunk(SB_,       asymb_n_fixation)            ! out


!! *slmtemp
    veg  => veg__mem
    sb   => sb__mem

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 Set values of soil profile variables to allow switching between simple_0d and simple1d in the 
    !!     following routines 
    !! 
    CALL calc_vertical_profile(lctlib, &
                               options, &
                               soil_config%nsoil_water, &
                               TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                               soil_depth_sl(:,:), &
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
                               rtm_decomposition_act(:,:), &
                               rmm_decomposition_act(:,:), &
                               rtm_nitrification_act(:,:), &
                               rmm_nitrification_act(:,:), &
                               rtm_denitrification_act(:,:), &
                               rtm_gasdiffusion_act(:,:), &
                               rmm_gasdiffusion_act(:,:), &
                               anaerobic_volume_fraction_act(:,:) )

    !> 2.0 Calculate potential rates of uptake and losses
    !!
    !> 2.1 Potential plant nutrient uptake rates
    !!
    CALL calc_potential_plant_uptake_rate(lctlib, &
                                          options, &
                                          soil_config%nsoil_water, &
                                          soil_depth_sl(:,:), &
                                          veg%pool%fine_root%carbon%ptr, &
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
    CALL calc_potential_mycorrhiza_uptake_rate(lctlib, &
                                               options, &
                                               soil_config%nsoil_water, &
                                               soil_depth_sl(:,:), &
                                               sb%pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk), &
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
                                           soil_config%nsoil_water, &
                                           rtm_decomposition_act(:,:), &
                                           rmm_decomposition_act(:,:), &
                                           sb%pool, &
                                           sb%loss )

    !> 2.5 Calculate the potential PO4 biomineralisation rate 
    !!
    CALL calc_potential_biomineralisation_rate(options, &                  
                                               soil_config%nsoil_water, &
                                               rtm_decomposition_act(:,:), &
                                               rmm_decomposition_act(:,:), &
                                               sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk), & 
                                               biomineralisation_po4(:,:) )

    !> 2.6 Loss pathways for the dynamic loss calculation scheme 
    !!
    IF(TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme)=='dynamic')THEN
       CALL calc_potential_loss_rates(options, &
                                   soil_config%nsoil_water, &
                                   TRIM(soil_config%soil_moisture_scheme), &
                                   TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                                   soil_depth_sl(:,:), &
                                   drainage_fraction(:), &
                                   percolation_sl(:,:), &
                                   rtm_nitrification_act(:,:), &
                                   rmm_nitrification_act(:,:), &
                                   rtm_denitrification_act(:,:), &
                                   km_denitrification_c_act(:,:), &
                                   km_denitrification_no3_act(:,:), &
                                   anaerobic_volume_fraction_act(:,:), &
                                   nh4_solute(:,:), &
                                   no3_solute(:,:), &
                                   po4_solute(:,:), &
                                   sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), &
                                   nitrification_no3(:,:), &
                                   denitrification_n2(:,:), &
                                   transport_nh4_solute(:,:), &
                                   transport_no3_solute(:,:), &
                                   transport_po4_solute(:,:), &
                                   leaching_nh4_solute(:), &
                                   leaching_no3_solute(:), &
                                   leaching_po4_solute(:) )  
    ENDIF

    !> 3.0 Actual fluxes of in and out of SOM 
    !!
 
    !> 3.1 Calculation of actual plant, mycorrhizal, microbial uptake rate as well as 
    !!     nitrification, denitrification and biomineralisation rates
    !!
    CALL calc_actual_flux_rates(options, &
                                soil_config%nsoil_water, &
                                nh4_solute(:,:), &
                                no3_solute(:,:), &
                                po4_solute(:,:), &
                                transport_nh4_solute(:,:), &
                                transport_no3_solute(:,:), &
                                transport_po4_solute(:,:), &
                                plant_uptake_nh4_sl(:,:), &
                                plant_uptake_no3_sl(:,:), &
                                plant_uptake_po4_sl(:,:), &
                                mycorrhiza_uptake_nh4_sl(:,:), &
                                mycorrhiza_uptake_no3_sl(:,:), &
                                mycorrhiza_uptake_po4_sl(:,:), &
                                nitrification_no3(:,:), &
                                denitrification_n2(:,:), &
                                biomineralisation_po4(:,:), &
                                sb%loss, &
                                mic_uptake_nh4_sl(:,:), &
                                mic_uptake_po4_sl(:,:), &
                                apparent_mic_cue(:,:) )
    
    !> 3.2 formation of SOM and associated respiration, net mineralisation fluxes
    !!          
    CALL calc_SOM_formation(options, &
                            soil_config%nsoil_water, &
                            apparent_mic_cue(:,:), &
                            mic_uptake_nh4_sl(:,:), &
                            mic_uptake_po4_sl(:,:), &
                            sb%loss, &
                            het_respiration(:,:), &
                            het_respiration_c13(:,:), &
                            het_respiration_c14(:,:), &
                            net_mineralisation_nh4(:,:), &
                            net_mineralisation_po4(:,:), &
                            net_mineralisation_nh4_n15(:,:), &
                            sb%formation )

    !> 3.3 formation of metabolic and strucutral litter by woody liter decay and 
    !!     associated respiration, net mineralisation fluxes    
    IF(lctlib%growthform==itree)THEN 
       CALL calc_litter_formation_from_woody_decay(options, &
                                                   soil_config%nsoil_water, &
                                                   sb%loss, &
                                                   het_respiration(:,:), &
                                                   het_respiration_c13(:,:), &
                                                   het_respiration_c14(:,:), &
                                                   net_mineralisation_nh4(:,:), &
                                                   net_mineralisation_po4(:,:), &
                                                   net_mineralisation_nh4_n15(:,:), &
                                                   sb%formation )
    ENDIF

    !> 3.4 transport due to bioturbation
    !!
    IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme)=='simple_1d')THEN
       k_bioturb(:,:) = calc_bioturbation_rate(options, soil_config%nsoil_water, bulk_dens_corr_sl(:,:))
       CALL calc_bioturbation_transport_wrapper(options, soil_config%nsoil_water, k_bioturb(:,:), soil_depth_sl(:,:), &
                                                sb%pool%microbial, sb%pool%residue, & 
                                                sb%transport%microbial, sb%transport%residue)
    ENDIF

    !> 4.0 calculate N losses and nitrification
    !! 
    CALL calc_actual_nloss_rate(options, &
                                soil_config%nsoil_water, &
                                TRIM(dsl4jsb_Config(SB_)%sb_model_scheme), &
                                TRIM(dsl4jsb_Config(SB_)%sb_nloss_scheme), &
                                soil_depth_sl(:,:), &
                                rtm_denitrification_act(:,:), &
                                rtm_gasdiffusion_act(:,:), &
                                rmm_gasdiffusion_act(:,:), &
                                noy(:,:), &
                                n2o(:,:), &
                                n2(:,:), &
                                net_mineralisation_nh4(:,:), &
                                transport_nh4_solute(:,:), &
                                transport_no3_solute(:,:), &
                                transport_po4_solute(:,:), &
                                nitrification_no3(:,:), &
                                denitrification_n2(:,:), &
                                nitrification_noy(:,:), &
                                nitrification_n2o(:,:), &
                                denitrification_noy(:,:), &
                                denitrification_n2o(:,:), &
                                volatilisation_n(:,:), &
                                transport_noy(:,:), &
                                transport_n2o(:,:), &
                                transport_n2(:,:), &
                                emission_noy(:), &
                                emission_n2o(:), &
                                emission_n2(:), &
                                leaching_nh4_solute(:), &
                                leaching_no3_solute(:), &
                                leaching_po4_solute(:) )

    !> 5.0 Mycorrhiza dynamics 
    !! 
    CALL calc_mycorrhiza_dynamics(lctlib, &
                                  options, &
                                  soil_config%nsoil_water, &
                                  soil_depth_sl(:,:), &
                                  mycorrhiza_uptake_nh4_sl(:,:), &
                                  mycorrhiza_uptake_no3_sl(:,:), &
                                  mycorrhiza_uptake_po4_sl(:,:), &
                                  sb%pool, &
                                  sb%formation, &
                                  sb%loss, &
                                  het_respiration(:,:), &
                                  het_respiration_c13(:,:), &
                                  het_respiration_c14(:,:), &
                                  veg%mycorrhiza_export)


    !> 6.0 Update soil litter pools with litter fall
    !!
    CALL calc_litter_partitioning(lctlib, &  
                                  options, &
                                  soil_config%nsoil_water, &
                                  soil_depth_sl(:,:), &
                                  root_fraction_act_sl(:,:), &
                                  veg%litterfall, &
                                  sb%formation )
  
    CALL calc_asymbiontic_n_fixation(lctlib, &
                                     options, &
                                     soil_config%nsoil_water, &
                                     TRIM(dsl4jsb_Config(VEG_)%bnf_scheme), & 
                                     nh4_solute(:,:), &
                                     no3_solute(:,:), &
                                     sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), &
                                     root_fraction_act_sl(:,:), &
                                     asymb_n_fixation(:,:))
                    
    !> 7.0 convert rates to µmol/ m3 / s (µmol / m2 / s for surface fluxes)
    !!
    het_respiration(:,:) = het_respiration(:,:) / dtime * 1000000._wp
    het_respiration_c13(:,:) = het_respiration_c13(:,:) / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) / dtime * 1000000._wp
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_sl(:,:) = plant_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_no3_sl(:,:) = plant_uptake_no3_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_po4_sl(:,:) = plant_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_nh4_sl(:,:) = mycorrhiza_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    mycorrhiza_uptake_no3_sl(:,:) = mycorrhiza_uptake_no3_sl(:,:) / dtime * 1000000._wp
    mycorrhiza_uptake_po4_sl(:,:) = mycorrhiza_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    volatilisation_n(:,:) = volatilisation_n(:,:) / dtime * 1000000._wp
    nitrification_no3(:,:) = nitrification_no3(:,:) / dtime * 1000000._wp 
    nitrification_noy(:,:) = nitrification_noy(:,:) / dtime * 1000000._wp 
    nitrification_n2o(:,:) = nitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_noy(:,:) = denitrification_noy(:,:) / dtime * 1000000._wp 
    denitrification_n2o(:,:) = denitrification_n2o(:,:) / dtime * 1000000._wp 
    denitrification_n2(:,:) = denitrification_n2(:,:) / dtime * 1000000._wp 
    transport_nh4_solute(:,:) = transport_nh4_solute(:,:) / dtime * 1000000._wp
    transport_no3_solute(:,:) = transport_no3_solute(:,:) / dtime * 1000000._wp
    transport_po4_solute(:,:) = transport_po4_solute(:,:) / dtime * 1000000._wp
    leaching_nh4_solute(:) = leaching_nh4_solute(:) / dtime * 1000000._wp
    leaching_no3_solute(:) = leaching_no3_solute(:) / dtime * 1000000._wp
    leaching_po4_solute(:) = leaching_po4_solute(:) / dtime * 1000000._wp
    transport_noy(:,:) = transport_noy(:,:) / dtime * 1000000._wp
    transport_n2o(:,:) = transport_n2o(:,:) / dtime * 1000000._wp
    transport_n2(:,:) = transport_n2(:,:) / dtime * 1000000._wp
    emission_noy(:) = emission_noy(:) / dtime * 1000000._wp
    emission_n2o(:) = emission_n2o(:) / dtime * 1000000._wp
    emission_n2(:) = emission_n2(:) / dtime * 1000000._wp

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
  SUBROUTINE calc_vertical_profile(lctlib, &                             ! in 
                                    options, &                         
                                    nsoil_water, &
                                    sb_model_scheme, &
                                    soil_depth_sl, &
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
                                    rtm_decomposition_act, &
                                    rmm_decomposition_act, &
                                    rtm_nitrification_act, &
                                    rmm_nitrification_act, &
                                    rtm_denitrification_act, &
                                    rtm_gasdiffusion_act, &
                                    rmm_gasdiffusion_act, &
                                    anaerobic_volume_fraction_act )

    USE mo_sb_constants,            ONLY: km_denitrification_c, km_denitrification_no3, km1_up_nh4, km1_up_no3, & 
                                          km1_up_po4, km2_up_nh4, km2_up_no3, km2_up_po4

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),                      INTENT(in)  :: lctlib                       !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_task_options),                    INTENT(in)  :: options
    INTEGER,                                     INTENT(in)  :: nsoil_water                  !< number of soil layers
    CHARACTER(len=*),                            INTENT(in)  :: sb_model_scheme              !< from sb config
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)  :: soil_depth_sl, &             !< soil depth of a layer [m]
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
                                                                anaerobic_volume_fraction_sl !< anaerobic volume fraction
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(out) :: root_fraction_act_sl, &      !< actual fraction of root mass in layer
                                                                km1_uptake_nh4_act, &        !< actual km1 for NH4 uptake [m3 / mol]
                                                                km2_uptake_nh4_act, &        !< actual km2 for NH4 uptake [mol / m3]
                                                                km1_uptake_no3_act, &        !< actual km1 for NO3 uptake [m3 / mol]
                                                                km2_uptake_no3_act, &        !< actual km2 for NO3 uptake [mol / m3]
                                                                km1_uptake_po4_act, &        !< actual km1 for PO4 uptake [m3 / mol]
                                                                km2_uptake_po4_act, &        !< actual km2 for PO4 uptake [mol / m3]
                                                                km_denitrification_c_act, &  !< actual km on C for denitrification [mol / m3]
                                                                km_denitrification_no3_act, & !< actual km on NO3 for denitrification [mol / m3]
                                                                rtm_decomposition_act, &     !< actual temperature response scalar for decomposition
                                                                rmm_decomposition_act, &     !< actual moisture response scalar for decomposition
                                                                rtm_nitrification_act, &     !< actual temperature response scalar for nitrification
                                                                rmm_nitrification_act, &     !< actual moisture response scalar for nitrification
                                                                rtm_denitrification_act, &   !< actual temperature response scalar for denitrification
                                                                rtm_gasdiffusion_act, &      !< actual temperature response scalar for gas diffusion
                                                                rmm_gasdiffusion_act, &      !< actual moisture response scalar for gas diffusion
                                                                anaerobic_volume_fraction_act !< actual anaerobic volume fraction
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc) :: root_depth_avg       !< average rooting depth [m]
    INTEGER                         :: isoil                !< loop over soil_config%nsoil_water
    CHARACTER(len=*), PARAMETER     :: routine = TRIM(modname)//':calc_vertical_profile'

    !> 1.0 Case of 0d model
    !!
    IF(sb_model_scheme == 'simple_0d') THEN
      root_fraction_act_sl(:,1)             = 1.0_wp
      root_fraction_act_sl(:,2:nsoil_water) = 0.0_wp
      ! set to default values to avoid division by zero below (2:nsoil_water are unused in 0D)
      km1_uptake_nh4_act(:,:) = km1_up_nh4 
      km1_uptake_no3_act(:,:) = km1_up_no3
      km1_uptake_po4_act(:,:) = km1_up_po4
      km2_uptake_nh4_act(:,:) = km2_up_nh4
      km2_uptake_no3_act(:,:) = km2_up_no3
      km2_uptake_po4_act(:,:) = km2_up_po4
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
         km_denitrification_c_act(:,1) = km_denitrification_c / soil_depth_sl(:,1) * root_depth_avg(:)
         km_denitrification_no3_act(:,1) = km_denitrification_no3 / soil_depth_sl(:,1) * root_depth_avg(:)
      ENDWHERE
      rtm_decomposition_act(:,:) = 0.0_wp 
      rmm_decomposition_act(:,:) = 0.0_wp 
      rtm_nitrification_act(:,:) = 0.0_wp 
      rmm_nitrification_act(:,:) = 0.0_wp 
      rtm_denitrification_act(:,:) = 0.0_wp 
      rtm_gasdiffusion_act(:,:) = 0.0_wp 
      rmm_gasdiffusion_act(:,:) = 0.0_wp 
      anaerobic_volume_fraction_act(:,:) = 0.0_wp 
      DO isoil = 1,nsoil_water
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
      ENDDO
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
      km_denitrification_c_act(:,:) = km_denitrification_c
      km_denitrification_no3_act(:,:) = km_denitrification_no3
      rtm_decomposition_act(:,:) = rtm_decomposition(:,:)
      rmm_decomposition_act(:,:) = rmm_decomposition(:,:) 
      rtm_nitrification_act(:,:) = rtm_nitrification(:,:)
      rmm_nitrification_act(:,:) = rmm_nitrification(:,:) 
      rtm_denitrification_act(:,:) = rtm_denitrification(:,:)
      rtm_gasdiffusion_act(:,:) = rtm_gasdiffusion(:,:)
      rmm_gasdiffusion_act(:,:) = rmm_gasdiffusion(:,:) 
      anaerobic_volume_fraction_act(:,:) = anaerobic_volume_fraction_sl(:,:) 
    ENDIF

  END SUBROUTINE calc_vertical_profile


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates potential plant uptake of NH4, NO3 and PO4, given soil concentrations,
  !! and fine root mass per layer 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_plant_uptake_rate(lctlib, &                     ! in
                                              options, &
                                              nsoil_water, &
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
    TYPE(t_lctlib_element),                      INTENT(in)  :: lctlib                     !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_task_options),                    INTENT(in)  :: options
    INTEGER,                                     INTENT(in)  :: nsoil_water                !< number of soil layers
    REAL(wp),                                    INTENT(in)  :: veg_pool_fine_root_carbon  !< fine root carbon [mol / m2]
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)  :: soil_depth_sl, &           !< soil depth per layer [m]
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
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(out) :: plant_uptake_nh4_sl, &     !< potential plant uptake of NH4 [mol / m3 / timestep]
                                                                plant_uptake_no3_sl, &     !< potential plant uptake of NO3 [mol / m3 / timestep]
                                                                plant_uptake_po4_sl        !< potential plant uptake of PO4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc),             INTENT(out) :: unit_uptake_n_pot, &       !< potential N uptake by fine roots per unit mass [µmol / mol / s]
                                                                unit_uptake_p_pot          !< potential P uptake by fine roots per unit mass [µmol / mol / s]
    ! 0.2 Local
    INTEGER                                     :: isoil        !< loop over nsoil_water
    REAL(wp), DIMENSION(options%nc,nsoil_water) :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_potential_plant_uptake_rate'

    unit_uptake_n_pot(:) = 0.0_wp
    unit_uptake_p_pot(:) = 0.0_wp

    DO isoil = 1,nsoil_water
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
      plant_uptake_nh4_sl(:,isoil) = hlp1(:,isoil) * veg_pool_fine_root_carbon / & 
                                     soil_depth_sl(:,isoil) * dtime/1.e6_wp
      plant_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * veg_pool_fine_root_carbon / & 
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
                                     veg_pool_fine_root_carbon / soil_depth_sl(:,isoil) * dtime/1.e6_wp
   
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
  SUBROUTINE calc_potential_mycorrhiza_uptake_rate(lctlib, &                     ! in
                                                   options, &
                                                   nsoil_water, &
                                                   soil_depth_sl, &
                                                   sb_pool_mycorrhiza_carbon, &
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
    TYPE(t_lctlib_element),                      INTENT(in)  :: lctlib                      !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_task_options),                    INTENT(in)  :: options                     !< dimensions...
    INTEGER,                                     INTENT(in)  :: nsoil_water                 !< number of soil layers
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)  :: soil_depth_sl, &            !< soil depth per layer [m]
                                                                sb_pool_mycorrhiza_carbon, &!< mycorrhizal biomass [mol C / m3]
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
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(out) :: mycorrhiza_uptake_nh4_sl, & !< potential mycorrhiza uptake of NH4 [mol / m3 / timestep]
                                                                mycorrhiza_uptake_no3_sl, & !< potential mycorrhiza uptake of NO3 [mol / m3 / timestep]
                                                                mycorrhiza_uptake_po4_sl    !< potential mycorrhiza uptake of PO4 [mol / m3 / timestep]
    ! 0.2 Local
    INTEGER                                     :: isoil        !< loop over nsoil_water
    REAL(wp), DIMENSION(options%nc,nsoil_water) :: hlp1, hlp2   !< helpers
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_potential_mycorrhiza_uptake_rate'

    
    DO isoil = 1,nsoil_water

      !> 1.0 pot. N uptake by mycorrhizal fungi (mol / m3 / timestep)
      !!
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
  
  END SUBROUTINE calc_potential_mycorrhiza_uptake_rate
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the potential rates of nitrification, denitrification and advection within the 
  !! soil column, given current soil N and P concentrations, soil temperature, moisture and 
  !! water percolation rate 
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_loss_rates(options, &
                                       nsoil_water, &
                                       soil_moisture_scheme, &
                                       sb_model_scheme, &
                                       soil_depth_sl, &
                                       drainage_fraction, &
                                       percolation_sl, &
                                       rtm_nitrification_act, &
                                       rmm_nitrification_act, &
                                       rtm_denitrification_act, &
                                       km_denitrification_c_act, &
                                       km_denitrification_no3_act, &
                                       anaerobic_volume_fraction_act, &
                                       nh4_solute, &
                                       no3_solute, &
                                       po4_solute, &
                                       sb_pool_microbial_carbon, &
                                       nitrification_no3, &
                                       denitrification_n2, &
                                       transport_nh4_solute, &
                                       transport_no3_solute, &
                                       transport_po4_solute, &
                                       leaching_nh4_solute, &
                                       leaching_no3_solute, &
                                       leaching_po4_solute  )

    USE mo_sb_constants,     ONLY: vmax_nitrification, vmax_denitrification, fleach_nh4, fleach_po4
    USE mo_sb_jsm_transport, ONLY: calc_liquid_phase_transport

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)  :: options
    INTEGER,                                     INTENT(in)  :: nsoil_water                        !< number of soil layers
    CHARACTER(len=*),                            INTENT(in)  :: soil_moisture_scheme               !< type of soil model scheme: bucket / layer
    CHARACTER(len=*),                            INTENT(in)  :: sb_model_scheme                    !< type of sb model scheme: simple0d / simple1d 
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)  :: soil_depth_sl                      !< soil layer depth [m]
    REAL(wp), DIMENSION(options%nc),             INTENT(in)  :: drainage_fraction                  !< fraction of water lost for 0d model [s-1]
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)  :: percolation_sl, &                  !< fraction of water lost for 1d model [s-1]
                                                                rtm_nitrification_act, &           !< actual temperature response for nitrification
                                                                rmm_nitrification_act, &           !< actual moisture response for nitrification
                                                                rtm_denitrification_act, &         !< actual temperature response for denitrification 
                                                                km_denitrification_c_act, &        !< actual km on C for denitrification
                                                                km_denitrification_no3_act, &      !< actual km on NO3 for denitrification
                                                                anaerobic_volume_fraction_act, &   !< actual anaerobic volume fraction
                                                                nh4_solute, &                      !< NH4 in solution [mol / m3]
                                                                no3_solute, &                      !< NO3 in solution [mol / m3]
                                                                po4_solute, &                      !< PO4 in solution [mol / m3]
                                                                sb_pool_microbial_carbon           !< microbial C [mol / m3]
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(out) :: nitrification_no3, &               !< nitrification rate [mol / m3 / timestep]
                                                                denitrification_n2, &              !< denitrification rate [mol / m3 / timestep] 
                                                                transport_nh4_solute, &            !< vertical transport of NH4 inside soil column [mol / m3 / timestep]
                                                                transport_no3_solute, &            !< vertical transport of NO3 inside soil column [mol / m3 / timestep]
                                                                transport_po4_solute               !< vertical transport of PO4 inside soil column [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc),             INTENT(out) :: leaching_nh4_solute, &             !< leaching of NH4 [mol / m2 / timestep]
                                                                leaching_no3_solute, &             !< leaching of NO3 [mol / m2 / timestep]
                                                                leaching_po4_solute                !< leaching of PO4 [mol / m2 / timestep]
    ! 0.2 Local
    ! 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_nloss_rate'

    !> 1.0 Potential nitrification and denitrification rates (according to Zaehle et al. 2011)
    !! potential nitrification rate from NH4 to NO3 [mol/m3/timestep]
    nitrification_no3(:,:)   = vmax_nitrification * rtm_nitrification_act(:,:) * rmm_nitrification_act(:,:) * & 
                               (1._wp - anaerobic_volume_fraction_act(:,:)) * nh4_solute(:,:) * dtime/1.e6_wp 
    !! potential denitrification rate from NO3 to N2 [mol/m3/timestep]
    denitrification_n2(:,:)  = vmax_denitrification * rtm_denitrification_act(:,:) * &
                               sb_pool_microbial_carbon(:,:) / &
                               ( km_denitrification_c_act(:,:) + sb_pool_microbial_carbon(:,:))  * &
                               anaerobic_volume_fraction_act(:,:) * no3_solute(:,:) / &
                               ( km_denitrification_no3_act(:,:) + no3_solute(:,:)) * dtime/1.e6_wp 

    !> 2.0. Leaching rate
    !!
    IF(sb_model_scheme == 'simple_0d')THEN
       ! For the 0D soil model, loss is proportional to the fractional waterloss of the soil column
       ! reduction of layer concentration in mol m-3 timestep-1
       transport_nh4_solute(:,1) = - fleach_nh4 * nh4_solute(:,1) * drainage_fraction(:) * dtime
       transport_no3_solute(:,1) = - no3_solute(:,1) * drainage_fraction(:) * dtime
       transport_po4_solute(:,1) = - fleach_po4 * po4_solute(:,1) * drainage_fraction(:) * dtime

    ELSE 
       IF(soil_moisture_scheme=='layer')THEN
         ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
         ! reduction of layer concentration in mol m-3 timestep-1, leaching rate in mol m-2 timestep-1
         CALL calc_liquid_phase_transport(options, nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          fleach_nh4 * nh4_solute(:,:),transport_nh4_solute(:,:),leaching_nh4_solute(:))
         CALL calc_liquid_phase_transport(options, nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          no3_solute(:,:),transport_no3_solute(:,:),leaching_no3_solute(:))
         CALL calc_liquid_phase_transport(options, nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                          fleach_po4 * po4_solute(:,:),transport_po4_solute(:,:),leaching_po4_solute(:))
       ENDIF
    ENDIF
  
  END SUBROUTINE calc_potential_loss_rates
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculate potential decomposition rate of SOM and litter pools according to 
  !!  dX/dt = - f(T)*g(M)/tau_pool * X
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_decomposition_rate(options, &                ! in
                                                nsoil_water, &
                                                rtm_decomposition_act, &
                                                rmm_decomposition_act, &
                                                sb_pool, &
                                                sb_loss )                ! out

    USE mo_sb_constants,            ONLY: tau_metabolic_litter, tau_structural_litter, tau_woody_litter, tau_fast, tau_slow

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)    :: options                  !< dimensions
    INTEGER,                                     INTENT(in)    :: nsoil_water              !< number of soil layers
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)    :: rtm_decomposition_act, & !< temperature response of decomposition
                                                                  rmm_decomposition_act    !< moisture response of decomposition 
    TYPE(t_slm_struct_sb_pool),                  INTENT(in)    :: sb_pool                  !< soil organic matter and litter pools [mol / m3]
    TYPE(t_slm_struct_sb_pool),                  INTENT(inout) :: sb_loss                  !< soil organic matter and litter loss [mol / m3]
    ! 0.2 Local
    ! 
    INTEGER                     :: ics,ice,iblk
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_decomposition_rate'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    !> 1.0 metabolic litter (= soluable pool)
    !! 
    sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) = sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                          sb_pool%soluable_litter%carbon%ptr(ics:ice,:,iblk)
    sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                            sb_pool%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                              sb_pool%soluable_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) = sb_loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                            sb_pool%soluable_litter%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) = sb_loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                            sb_pool%soluable_litter%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_metabolic_litter * dtime * &
                                                              sb_pool%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 2.0 structural litter (= polymeric pool)
    !! 
    sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                          sb_pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                            sb_pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                              sb_pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                            sb_pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                            rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                            sb_pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                              rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_structural_litter * dtime * &
                                                              sb_pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.0 woody litter
    !! 
    sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                      rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                      sb_pool%woody_litter%carbon%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                        rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                          sb_pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                        rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                        rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                          rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_woody_litter * dtime * &
                                                          sb_pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 4.0 fast SOM pool (=microbial pool)
    !!
    sb_loss%microbial%carbon%ptr(ics:ice,:,iblk) = sb_loss%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                   rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                   sb_pool%microbial%carbon%ptr(ics:ice,:,iblk)
    sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                       rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                       sb_pool%microbial%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%microbial%carbon13%ptr(ics:ice,:,iblk) = sb_loss%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%microbial%carbon14%ptr(ics:ice,:,iblk) = sb_loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                       rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_fast * dtime * &
                                                       sb_pool%microbial%nitrogen15%ptr(ics:ice,:,iblk)

    !> 5.0 slow SOM pool (=residue pool)
    !!
    sb_loss%residue%carbon%ptr(ics:ice,:,iblk) = sb_loss%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                 rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_slow * dtime * &
                                                 sb_pool%residue%carbon%ptr(ics:ice,:,iblk)
    sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                   rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_slow * dtime * &
                                                   sb_pool%residue%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_slow * dtime * &
                                                     sb_pool%residue%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%residue%carbon13%ptr(ics:ice,:,iblk) = sb_loss%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                   rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_slow * dtime * &
                                                   sb_pool%residue%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%residue%carbon14%ptr(ics:ice,:,iblk) = sb_loss%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                   rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_slow * dtime * &
                                                   sb_pool%residue%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                     rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_slow * dtime * &
                                                     sb_pool%residue%nitrogen15%ptr(ics:ice,:,iblk)  
  
  END SUBROUTINE calc_potential_decomposition_rate

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculate potential biomineralisation rate 
  !!  dP/dt = - f(T)*g(M)/tau_pool * P
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_biomineralisation_rate(options, &                      ! in 
                                                   nsoil_water, &
                                                   rtm_decomposition_act, &
                                                   rmm_decomposition_act, &
                                                   sb_pool_residue_phosphorus, & 
                                                   biomineralisation_po4 )         ! out 

    USE mo_sb_constants,            ONLY: tau_biomineralisation 

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)  :: options                    !< dimensions
    INTEGER,                                     INTENT(in)  :: nsoil_water                !< number of soil layers
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)  :: rtm_decomposition_act, &   !< temperature response of decomposition
                                                                rmm_decomposition_act, &   !< moisture response of decomposition
                                                                sb_pool_residue_phosphorus !< slow SOM P pool [mol / m3] 
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(out) :: biomineralisation_po4      !< potential biomineralisation rate [mol / m3 / timestep] 
    ! 0.2 Local
    ! 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_potential_biomineralisation_rate'

    !> 1.0 Calculate potential P biomineralisation from slow organic P pool (mol /m3 / timestep)
    biomineralisation_po4(:,:) = rtm_decomposition_act(:,:)*rmm_decomposition_act(:,:)/tau_biomineralisation * & 
                                 dtime * sb_pool_residue_phosphorus(:,:) 

  END SUBROUTINE calc_potential_biomineralisation_rate
  
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
  SUBROUTINE calc_actual_flux_rates(options, &
                                     nsoil_water, &
                                     nh4_solute, &
                                     no3_solute, &
                                     po4_solute, &
                                     transport_nh4_solute, &
                                     transport_no3_solute, &
                                     transport_po4_solute, &
                                     plant_uptake_nh4_sl, &
                                     plant_uptake_no3_sl, &
                                     plant_uptake_po4_sl, &
                                     mycorrhiza_uptake_nh4_sl, &
                                     mycorrhiza_uptake_no3_sl, &
                                     mycorrhiza_uptake_po4_sl, &
                                     nitrification_no3, &
                                     denitrification_n2, &
                                     biomineralisation_po4, &
                                     sb_loss, &
                                     mic_uptake_nh4_sl, &
                                     mic_uptake_po4_sl, &
                                     apparent_mic_cue )

    USE mo_sb_constants,            ONLY: microbial_cue_max, microbial_cue_min, microbial_cn, microbial_np, &
                                          microbial_nue, microbial_pue

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)    :: options                       !< dimensions
    INTEGER,                                     INTENT(in)    :: nsoil_water                   !< number of soil layers
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)    :: nh4_solute, &                 !< NH4 in solution [mol / m3] 
                                                                  no3_solute, &                 !< NO3 in solution [mol / m3]
                                                                  po4_solute, &                 !< PO4 in solution [mol / m3]
                                                                  transport_nh4_solute, &       !< transport of NH4 [mol / m3 / timestep]
                                                                  transport_no3_solute, &       !< transport of NO3 [mol / m3 / timestep]
                                                                  transport_po4_solute          !< transport of PO4 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(inout) :: plant_uptake_nh4_sl, &        !< plant NH4 uptake [mol / m3 / timestep]
                                                                  plant_uptake_no3_sl, &        !< plant NO3 uptake [mol / m3 / timestep]
                                                                  plant_uptake_po4_sl, &        !< plant PO4 uptake [mol / m3 / timestep]
                                                                  mycorrhiza_uptake_nh4_sl, &   !< mycorrhizal NH4 uptake [mol / m3 / timestep]
                                                                  mycorrhiza_uptake_no3_sl, &   !< mycorrhizal NO3 uptake [mol / m3 / timestep]
                                                                  mycorrhiza_uptake_po4_sl, &   !< mycorrhizal PO4 uptake [mol / m3 / timestep]
                                                                  nitrification_no3, &          !< nitrification rate [mol / m3 / timestep]
                                                                  denitrification_n2, &         !< denitrification rate [mol / m3 / timestep]
                                                                  biomineralisation_po4         !< biomineralisation rate [mol / m3 / timestep]
    TYPE(t_slm_struct_sb_pool),                  INTENT(inout) :: sb_loss                       !< decomposition of SOM and litter pools [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(out)   :: mic_uptake_nh4_sl, &          !< microbial immobilisation of NH4 [mol / m3 / timestep]
                                                                  mic_uptake_po4_sl, &          !< microbial immobilisation of PO4 [mol / m3 / timestep]
                                                                  apparent_mic_cue              !< apparent microbial carbon-use-efficiency
    ! 0.2 Local
    INTEGER                                     :: ics, ice, iblk         !< chunks [ics:ice] and block numbers
    REAL(wp), DIMENSION(options%nc,nsoil_water) :: hlp1, hlp2, growth_c, growth_n
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_actual_flux_rates'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    !> 1.0 formation of fast SOM from litter 
    !!

    !> 1.1 calculate potential growth rate of fast SOM pool, and associated potential microbial 
    !!     immobilisation rate as difference between maximal growth rate and N in litter
    growth_c(:,:) = microbial_cue_max * & 
                    (sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                     sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk)) 
    
    mic_uptake_nh4_sl(:,:) = growth_c(:,:) / microbial_cn - & 
                             microbial_nue * (sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                              sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk))

    !> 1.2  Limit microbial, mycorrhizal and plant uptake, as well as nitrification to available NH4 minus transport losses
    !!
    WHERE(mic_uptake_nh4_sl(:,:)>=0.0_wp)
       hlp1(:,:) = mic_uptake_nh4_sl(:,:) + plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + mycorrhiza_uptake_nh4_sl(:,:)
       hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:)
       WHERE(hlp2(:,:) < hlp1(:,:)) 
           WHERE(hlp1(:,:) > eps8)
              mic_uptake_nh4_sl(:,:)   = hlp2(:,:) * mic_uptake_nh4_sl(:,:) / hlp1(:,:) 
              plant_uptake_nh4_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_sl(:,:) / hlp1(:,:)
              nitrification_no3(:,:)   = hlp2(:,:) * nitrification_no3(:,:) / hlp1(:,:)
              mycorrhiza_uptake_nh4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_nh4_sl(:,:) / hlp1(:,:)
           ELSEWHERE
              mic_uptake_nh4_sl(:,:) = 0.0_wp
              plant_uptake_nh4_sl(:,:) = 0.0_wp 
              nitrification_no3(:,:) = 0.0_wp 
              mycorrhiza_uptake_nh4_sl(:,:) = 0.0_wp
           ENDWHERE 
        ENDWHERE
    ELSEWHERE
       hlp1(:,:) = plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) + mycorrhiza_uptake_nh4_sl(:,:)
       hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:) - mic_uptake_nh4_sl(:,:)
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
   
    !> 2.1 calculate C growth based on N availability
    !!
    growth_n(:,:) = microbial_cn * ( mic_uptake_nh4_sl(:,:) + & 
                                     microbial_nue * (sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                      sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk))) 

    !> 2.2 Adjust carbon-use efficiency to match N availability by unchanged litter decomposition
    !!
    WHERE(growth_c(:,:) <= growth_n(:,:)) ! enough N is available -> constant CUE
       apparent_mic_cue(:,:) = microbial_cue_max
    ELSEWHERE(growth_c(:,:) > eps8) ! not enough N is available -> changed CUE (within min/max CUE)
       apparent_mic_cue(:,:) = MAX(microbial_cue_min, & 
                                   growth_n(:,:) / ( growth_c(:,:) / microbial_cue_max ))
    ELSEWHERE ! growth C is too small -> CUE does not matter
       apparent_mic_cue(:,:) = microbial_cue_min 
    ENDWHERE

    !> 2.3 Calculate actual litter loss rate, calculated such that litter decomposition
    !! matches the actual immobilisation flux and CN uptake from litter given actual CUE
    hlp1(:,:) = apparent_mic_cue(:,:) / microbial_cn * &
                      (sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
                       sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk)) - & 
                microbial_nue * & 
                       (sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                        sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)) 
    hlp2(:,:) = 1.0_wp
    WHERE(mic_uptake_nh4_sl(:,:)<hlp1(:,:))
       WHERE(hlp1(:,:)>eps8)
           hlp2(:,:) = mic_uptake_nh4_sl(:,:) / hlp1(:,:)
       ELSEWHERE
           hlp2(:,:) = 0.0_wp
       ENDWHERE
       sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk)
       sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)
       sb_loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk)
       sb_loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk)
       sb_loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk)
       sb_loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk)

       sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)
    ENDWHERE 

    !> 3.0 calculate potential microbial P uptake from soluable pool (mol /m3 / timestep) 
    mic_uptake_po4_sl(:,:) =  apparent_mic_cue(:,:) / microbial_cn / microbial_np * & 
                                (sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + &
                                  sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk)) - &
                              microbial_pue * & 
                                (sb_loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                                 sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)) 

    !> 3.1 constrain P uptake by microbes and plants to the available sum of P from P_soluable, 
    !!     biomineralisation and transport (mol /m3 / timestep)   
    WHERE(mic_uptake_po4_sl(:,:)>0.0_wp) ! if microbes take up P, constrain both microbial and plant uptake to available P 
       hlp1(:,:) = mic_uptake_po4_sl(:,:) + plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:)
       hlp2(:,:) = po4_solute(:,:) + transport_po4_solute(:,:) + biomineralisation_po4(:,:)
       WHERE(hlp2(:,:) < hlp1(:,:))
          WHERE(hlp1(:,:) > eps8)
             mic_uptake_po4_sl(:,:)        = hlp2(:,:) * mic_uptake_po4_sl(:,:) / hlp1(:,:) 
             plant_uptake_po4_sl(:,:)      = hlp2(:,:) * plant_uptake_po4_sl(:,:) / hlp1(:,:) 
             mycorrhiza_uptake_po4_sl(:,:) = hlp2(:,:) * mycorrhiza_uptake_po4_sl(:,:) / hlp1(:,:) 
          ELSEWHERE
             mic_uptake_po4_sl(:,:)        = 0.0_wp
             plant_uptake_po4_sl(:,:)      = 0.0_wp
             mycorrhiza_uptake_po4_sl(:,:) = 0.0_wp 
          ENDWHERE
       ENDWHERE
    ELSEWHERE ! if microbes release P, constrain plant uptake to available P, note the negative sign of mic_uptake_po4_sl(:,:)
       hlp1(:,:) = plant_uptake_po4_sl(:,:) + mycorrhiza_uptake_po4_sl(:,:) 
       hlp2(:,:) = po4_solute(:,:) + transport_po4_solute(:,:) + biomineralisation_po4(:,:) - mic_uptake_po4_sl(:,:)
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
    biomineralisation_po4(:,:) = MAX(mic_uptake_po4_sl(:,:) + plant_uptake_po4_sl(:,:) - & 
                                     (po4_solute(:,:) + transport_po4_solute(:,:)),0.0_wp)
    sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                     biomineralisation_po4(:,:)
  
  END SUBROUTINE calc_actual_flux_rates
  
  
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
                                nsoil_water, &
                                apparent_mic_cue, &
                                mic_uptake_nh4_sl, &
                                mic_uptake_po4_sl, &
                                sb_loss, &
                                het_respiration, &               ! out
                                het_respiration_c13, &
                                het_respiration_c14, &
                                net_mineralisation_nh4, &
                                net_mineralisation_po4, &
                                net_mineralisation_nh4_n15, &
                                sb_formation  )
  
    USE mo_sb_constants,        ONLY: residue_cue, microbial_cue_min, microbial_np, microbial_cn, microbial_nue, microbial_pue

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)    :: options                       !< dimensions
    INTEGER,                                     INTENT(in)    :: nsoil_water                   !< number of soil layers 
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(in)    :: apparent_mic_cue, &           !< apparent carbon-use efficiency of microbial uptake
                                                                  mic_uptake_nh4_sl, &          !< microbial NH4 uptake [mol / m3 / timestep]
                                                                  mic_uptake_po4_sl             !< microbial PO4 uptake [mol / m3 / timestep]
    TYPE(t_slm_struct_sb_pool),                  INTENT(in)    :: sb_loss                       !< turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(inout) :: het_respiration, &            !< heterotrophic respiration [mol C / m3 / timestep]
                                                                  het_respiration_c13, &        !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                                  het_respiration_c14, &        !< heterotrophic respiration [µmol 14C / m3 / timestep]
                                                                  net_mineralisation_nh4, &     !< net mineralisation of NH4 [mol / m3 / timestep]
                                                                  net_mineralisation_po4, &     !< net mineralisation of PO4 [mol / m3 / timestep]
                                                                  net_mineralisation_nh4_n15    !< net mineralisation of 15NH4 [mol / m3 / timestep]
    TYPE(t_slm_struct_sb_pool),                  INTENT(inout) :: sb_formation                  !< formation of SOM and litter pools  [mol / m3 / timestep]
    
    ! 0.2 Local
    ! 
    INTEGER                     :: ics, ice, iblk         !< chunks [ics:ice] and block numbers
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_SOM_formation'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    
    !> 1.0 formation of fast SOM from litter 
    !! 
    sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                        apparent_mic_cue(:,:) * & 
                                                        (sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + &  
                                                         sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk)) 
    sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                          microbial_nue * &
                                                          (sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + &  
                                                           sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)) + & 
                                                          mic_uptake_nh4_sl(:,:) 
    sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                            microbial_pue * &
                                                            (sb_loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + &  
                                                             sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)) + & 
                                                            mic_uptake_po4_sl(:,:) 
    sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                        apparent_mic_cue(:,:) * & 
                                                        (sb_loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + &  
                                                         sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)) 
    sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                        apparent_mic_cue(:,:) * & 
                                                        (sb_loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + &  
                                                         sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk))
!! to be moved to N15 calculation
    WHERE(sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)+sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                              (sb_loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                               sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)) / & 
                              (sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                               sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)) *  & 
                              apparent_mic_cue(:,:) / microbial_cn * & 
                              (sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + &  
                               sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk))
    ENDWHERE

    het_respiration(:,:)     = het_respiration(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                                                      (sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + &  
                                                       sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk))
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                                                          (sb_loss%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + &  
                                                           sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)) 
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                                                          (sb_loss%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + &  
                                                           sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk))
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) - mic_uptake_nh4_sl(:,:) + &  
                                                                (1.0_wp - microbial_nue) * &
                                                                (sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + &  
                                                                 sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk))
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) - mic_uptake_po4_sl(:,:) + & 
                                                                (1.0_wp - microbial_pue) * &
                                                                (sb_loss%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + &  
                                                                 sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk))
 
    WHERE(sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk)+sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       net_mineralisation_nh4_n15(:,:)     = net_mineralisation_nh4_n15(:,:) + & 
                              sb_loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                              sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                              (sb_loss%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
                               sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)) / & 
                              (sb_loss%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                               sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)) *  & 
                              apparent_mic_cue(:,:) / microbial_cn * & 
                              (sb_loss%soluable_litter%carbon%ptr(ics:ice,:,iblk) + &  
                               sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk))
    ENDWHERE

    !> 2.0 formation of slow SOM (=residue pool) from fast SOM decay
    !! 
    sb_formation%residue%carbon%ptr(ics:ice,:,iblk) = sb_formation%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                      residue_cue * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    sb_formation%residue%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue / microbial_cn * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    sb_formation%residue%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue / microbial_cn / microbial_np * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    sb_formation%residue%carbon13%ptr(ics:ice,:,iblk) = sb_formation%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue * sb_loss%microbial%carbon13%ptr(ics:ice,:,iblk)
    sb_formation%residue%carbon14%ptr(ics:ice,:,iblk) = sb_formation%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue * sb_loss%microbial%carbon14%ptr(ics:ice,:,iblk)
    WHERE(sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) > eps8)
       sb_formation%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%residue%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                             sb_loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) / &  
                                                             sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) * & 
                                                             residue_cue / microbial_cn * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE
    het_respiration(:,:)     = het_respiration(:,:)     + (1._wp - residue_cue) * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - residue_cue) * sb_loss%microbial%carbon13%ptr(ics:ice,:,iblk)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - residue_cue) * sb_loss%microbial%carbon14%ptr(ics:ice,:,iblk)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb_loss%microbial%phosphorus%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn / microbial_np * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    WHERE(sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) > eps8)
       net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + & 
                                    sb_loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) - &
                                    sb_loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) / sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) * & 
                                    residue_cue / microbial_cn * sb_loss%microbial%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE

    !> 3.0 formation of fast SOM (=microbial pool) from slow SOM decay 
    !! 
    sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                        microbial_cue_min * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                          microbial_cue_min / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                            microbial_cue_min / microbial_cn / microbial_np * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                          microbial_cue_min * sb_loss%residue%carbon13%ptr(ics:ice,:,iblk)
    sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                          microbial_cue_min * sb_loss%residue%carbon14%ptr(ics:ice,:,iblk)
    WHERE(sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                               sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) / &  
                                                               sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) * &  
                                                               microbial_cue_min / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE

    het_respiration(:,:)      = het_respiration(:,:)     + (1._wp - microbial_cue_min) * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    het_respiration_c13(:,:)  = het_respiration_c13(:,:) + (1._wp - microbial_cue_min) * sb_loss%residue%carbon13%ptr(ics:ice,:,iblk)
    het_respiration_c14(:,:)  = het_respiration_c14(:,:) + (1._wp - microbial_cue_min) * sb_loss%residue%carbon14%ptr(ics:ice,:,iblk)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) - &
                                                                microbial_cue_min / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) - &
                                                                microbial_cue_min / microbial_cn / microbial_np * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    WHERE(sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + & 
                                    sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) - &
                                    sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) / sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) * &  
                                    microbial_cue_min / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE  
    
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
                                                    nsoil_water, &
                                                    sb_loss, &
                                                    het_respiration, &               ! out
                                                    het_respiration_c13, &
                                                    het_respiration_c14, &
                                                    net_mineralisation_nh4, &
                                                    net_mineralisation_po4, &
                                                    net_mineralisation_nh4_n15, &
                                                    sb_formation  )
  
    USE mo_sb_constants,        ONLY: microbial_cue_min, microbial_nue, microbial_pue, & 
                                      fc_soluable_max, k_fc_soluable, k_fn_soluable, k_fp_soluable, lc_woody_litter

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),                    INTENT(in)    :: options                       !< dimensions
    INTEGER,                                     INTENT(in)    :: nsoil_water                   !< number of soil layers 
    TYPE(t_slm_struct_sb_pool),                  INTENT(in)    :: sb_loss                       !< turnover of SOM and litter pools  [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,nsoil_water), INTENT(inout) :: het_respiration, &            !< heterotrophic respiration [mol C / m3 / timestep]
                                                                  het_respiration_c13, &        !< heterotrophic respiration [mol 13C / m3 / timestep]
                                                                  het_respiration_c14, &        !< heterotrophic respiration [µmol 14C / m3 / timestep]
                                                                  net_mineralisation_nh4, &     !< net mineralisation of NH4 [mol / m3 / timestep]
                                                                  net_mineralisation_po4, &     !< net mineralisation of PO4 [mol / m3 / timestep]
                                                                  net_mineralisation_nh4_n15    !< net mineralisation of 15NH4 [mol / m3 / timestep]
    TYPE(t_slm_struct_sb_pool),                  INTENT(inout) :: sb_formation                  !< formation of SOM and litter pools  [mol / m3 / timestep]
    
    ! 0.2 Local
    ! 
    REAL(wp), DIMENSION(options%nc,nsoil_water) :: fc_soluable_woody_litter, & !< fraction of non-respired woody litter becoming metabolic litter 
                                                   fn_soluable_woody_litter, & 
                                                   fp_soluable_woody_litter
    INTEGER                                     :: ics, ice, iblk              !< chunks [ics:ice] and block numbers
    CHARACTER(len=*), PARAMETER                 :: routine = TRIM(modname)//':calc_litter_formation_from_woody_decay'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    
    !> 1.0 Calculate partitioning of decaying woody litter into metabolic and structural litter
    !!
    fc_soluable_woody_litter(:,:) = 0.0_wp
    fn_soluable_woody_litter(:,:) = 0.0_wp
    fp_soluable_woody_litter(:,:) = 0.0_wp
    WHERE(sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) > eps8) 
       fc_soluable_woody_litter(:,:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_woody_litter * &
           microbial_cue_min / microbial_nue * & 
           sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) / sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk))
       WHERE(fc_soluable_woody_litter(:,:)>eps8)
          fn_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:))/(k_fn_soluable*fc_soluable_woody_litter(:,:)))
          fp_soluable_woody_litter(:,:) = 1._wp / (1._wp + (1._wp - fc_soluable_woody_litter(:,:))/(k_fp_soluable*fc_soluable_woody_litter(:,:)))
       ENDWHERE
    ENDWHERE

    !> 2.0 formation of metabolic and structural litter 
    !!     taking account of the processing losses of C, N and P
    sb_formation%soluable_litter%carbon%ptr(ics:ice,:,iblk) = sb_formation%soluable_litter%carbon%ptr(ics:ice,:,iblk) + & 
       fc_soluable_woody_litter(:,:) * microbial_cue_min * sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) 
    sb_formation%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%soluable_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
       fn_soluable_woody_litter(:,:) * microbial_nue * sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) 
    sb_formation%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%soluable_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
       fp_soluable_woody_litter(:,:) * microbial_pue * sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) 
    sb_formation%soluable_litter%carbon13%ptr(ics:ice,:,iblk) = sb_formation%soluable_litter%carbon13%ptr(ics:ice,:,iblk) + & 
       fc_soluable_woody_litter(:,:) * microbial_cue_min * sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) 
    sb_formation%soluable_litter%carbon14%ptr(ics:ice,:,iblk) = sb_formation%soluable_litter%carbon14%ptr(ics:ice,:,iblk) + & 
       fc_soluable_woody_litter(:,:) * microbial_cue_min * sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) 
    sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
       fn_soluable_woody_litter(:,:) * microbial_nue * sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) 

    sb_formation%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
       (1.0_wp - fc_soluable_woody_litter(:,:)) * microbial_cue_min * sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) 
    sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
       (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue * sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) 
    sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
       (1.0_wp - fp_soluable_woody_litter(:,:)) * microbial_pue * sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) 
    sb_formation%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + & 
       (1.0_wp - fc_soluable_woody_litter(:,:)) * microbial_cue_min * sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) 
    sb_formation%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + & 
       (1.0_wp - fc_soluable_woody_litter(:,:)) * microbial_cue_min * sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) 
    sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + & 
       (1.0_wp - fn_soluable_woody_litter(:,:)) * microbial_nue * sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) 

    !> 3.0 record losses to heterotrophic respiration and mineralisation
    !!
    het_respiration(:,:)            = het_respiration(:,:) + & 
                                      (1.0_wp - microbial_cue_min) * sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)
    net_mineralisation_nh4(:,:)     = net_mineralisation_nh4(:,:) + &
                                      (1.0_wp - microbial_nue) * sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)
    net_mineralisation_po4(:,:)     = net_mineralisation_po4(:,:) + &
                                      (1.0_wp - microbial_pue) * sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk)
    het_respiration_c13(:,:)        = het_respiration_c13(:,:) + & 
                                      (1.0_wp - microbial_cue_min) * sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk)
    het_respiration_c14(:,:)        = het_respiration_c14(:,:) + & 
                                      (1.0_wp - microbial_cue_min) * sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + &
                                      (1.0_wp - microbial_nue) * sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)
 
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
  SUBROUTINE calc_bioturbation_transport_wrapper(options, nsoil, k_bioturb, soil_depth_sl, &
                                                 pool_microbial, pool_residue, &
                                                 transport_microbial, transport_residue)

    USE mo_sb_jsm_transport,        ONLY: calc_bioturbation_transport

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: k_bioturb              !< diffusion factor for bioturbation 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_sl          !< depth of each soil layer [m] 
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_microbial         !< microbial pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(in)    :: pool_residue           !< microbial residue pool [mol/m3]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_microbial        !< rate of microbial transport [mol/m3/timestep]
    TYPE(t_slm_struct_sb_element),         INTENT(inout) :: transport_residue          !< rate of microbial residue transport [mol/m3/timestep]
    ! 0.2 Local
    ! 
    INTEGER                               :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_liquid_phase_transport_wrapper'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice

    !> 1.0 fast SOM pool
    !!
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%carbon%ptr(ics:ice,:,iblk), &
                                        transport_microbial%carbon%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%nitrogen%ptr(ics:ice,:,iblk), &
                                        transport_microbial%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%phosphorus%ptr(ics:ice,:,iblk), &
                                        transport_microbial%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%carbon13%ptr(ics:ice,:,iblk), &
                                        transport_microbial%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%carbon14%ptr(ics:ice,:,iblk), &
                                        transport_microbial%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_microbial%nitrogen15%ptr(ics:ice,:,iblk), &
                                        transport_microbial%nitrogen15%ptr(ics:ice,:,iblk))

    !> 1.1 slow SOM pool
    !!
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_residue%carbon%ptr(ics:ice,:,iblk), &
                                        transport_residue%carbon%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_residue%nitrogen%ptr(ics:ice,:,iblk), &
                                        transport_residue%nitrogen%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_residue%phosphorus%ptr(ics:ice,:,iblk), &
                                        transport_residue%phosphorus%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_residue%carbon13%ptr(ics:ice,:,iblk), &
                                        transport_residue%carbon13%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_residue%carbon14%ptr(ics:ice,:,iblk), &
                                        transport_residue%carbon14%ptr(ics:ice,:,iblk))
    CALL calc_bioturbation_transport(options,nsoil,k_bioturb,soil_depth_sl, & 
                                        pool_residue%nitrogen15%ptr(ics:ice,:,iblk), &
                                        transport_residue%nitrogen15%ptr(ics:ice,:,iblk))

  END SUBROUTINE calc_bioturbation_transport_wrapper

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
  SUBROUTINE calc_actual_nloss_rate(options, &
                                    nsoil, &
                                    sb_model_scheme, &
                                    sb_nloss_scheme, &
                                    soil_depth_sl, &
                                    rtm_denitrification_act, &
                                    rtm_gasdiffusion_act, &
                                    rmm_gasdiffusion_act, &
                                    noy, &
                                    n2o, &
                                    n2, &
                                    net_mineralisation_nh4, &
                                    transport_nh4_solute, &
                                    transport_no3_solute, &
                                    transport_po4_solute, &
                                    nitrification_no3, &
                                    denitrification_n2, &
                                    nitrification_noy, &
                                    nitrification_n2o, &
                                    denitrification_noy, &
                                    denitrification_n2o, &
                                    volatilisation_n, &
                                    transport_noy, &
                                    transport_n2o, &
                                    transport_n2, &
                                    emission_noy, &
                                    emission_n2o, &
                                    emission_n2, &
                                    leaching_nh4_solute, &
                                    leaching_no3_solute, &
                                    leaching_po4_solute )

    USE mo_sb_constants,            ONLY: floss_nmin, fnit_nmin, f_nit_noy, f_nit_n2o, f_denit_noy, f_denit_n2o
  
    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),              INTENT(in)     :: options
    INTEGER,                               INTENT(in)     :: nsoil                      !< number of soil layers
    CHARACTER(len=*),                      INTENT(in)     :: sb_model_scheme            !< sb model type 
    CHARACTER(len=*),                      INTENT(in)     :: sb_nloss_scheme            !< sb Nloss scheme  
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)     :: soil_depth_sl, &           !< soil depth per layer [m]
                                                             rtm_denitrification_act, & !< temperature modifier for denitrification
                                                             rtm_gasdiffusion_act, &    !< temperature modifier for gas diffusion
                                                             rmm_gasdiffusion_act, &    !< moisture modifier for gas diffusion       
                                                             noy, &                     !< NOy concentration [mol / m3]
                                                             n2o, &                     !< N2O concentration [mol / m3]
                                                             n2, &                      !< N2 concentration [mol / m3]
                                                             net_mineralisation_nh4, &  !< net mineralisation of NH4 [mol / m3 / timestep]
                                                             transport_nh4_solute, &    !< vertical transport of NH4 [mol / m3 / timestep]
                                                             transport_no3_solute, &    !< vertical transport of NO3 [mol / m3 / timestep]
                                                             transport_po4_solute       !< vertical transport of PO4 [mol / m3 / timestep]
     REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: nitrification_no3, &       !< nitrification to NO3 [mol / m3 / timestep]
                                                             denitrification_n2         !< denitrification to N2 [mol / m3 / timestep]
     REAL(wp), DIMENSION(options%nc,nsoil), INTENT(out)   :: nitrification_noy, &       !< nitrification to NOy [mol / m3 / timestep]
                                                             nitrification_n2o, &       !< nitrification to N2O [mol / m3 / timestep]
                                                             denitrification_noy, &     !< denitrification to NOy [mol / m3 / timestep]
                                                             denitrification_n2o, &     !< denitrification to N2O [mol / m3 / timestep]
                                                             volatilisation_n, &        !< volatilisation of NH4 [mol / m3 / timestep]
                                                             transport_noy, &           !< vertical transport of NOy [mol / m3 / timestep]
                                                             transport_n2o, &           !< vertical transport of N2O [mol / m3 / timestep]
                                                             transport_n2               !< vertical transport of N2 [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc),        INTENT(out)   :: emission_noy, &            !< soil efflux of NOy [mol / m2 / timestep]
                                                             emission_n2o, &            !< soil efflux of N2O [mol / m2 / timestep]
                                                             emission_n2, &             !< soil efflux of N2 [mol / m2 / timestep]
                                                             leaching_nh4_solute, &     !< leaching of NH4 [mol / m2 / timestep]
                                                             leaching_no3_solute, &     !< leaching of NO3 [mol / m2 / timestep]
                                                             leaching_po4_solute        !< leaching of PO4 [mol / m2 / timestep]
    ! 0.2 Local
    ! 
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_actual_nloss_rate'
    
    
    !> 1.0 calculate N losses and nitrification
    !! 
    SELECT CASE(TRIM(sb_nloss_scheme))
    CASE ('fixed') 
       !> 1.1 case of fixed N losses: N loss proportional to net mineralisation
       !!
       WHERE(net_mineralisation_nh4(:,:)>0.0_wp)
          ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
          volatilisation_n(:,:)  = net_mineralisation_nh4(:,:) * floss_nmin 
          ! NH4 transformed to NO3 [mol m-3 timestep-1]
          nitrification_no3(:,:) = net_mineralisation_nh4(:,:) * (1._wp - floss_nmin) * fnit_nmin 
       ENDWHERE
       ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
       emission_n2(:) = SUM(volatilisation_n(:,:) * soil_depth_sl(:,:),DIM=2)
    CASE ('dynamic') 
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

       ! gaseous loss of NOy, N2O, and N2, for now following Xu-Ri and Prentice 2008, thus ignoring any actual gas-transport
       ! and instead assume that from each layer a certain fraction excapes directly to the atmosphere.
       ! [mol m-3 timestep-1] 
       transport_noy(:,:) = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * noy(:,:) / one_day * dtime
       transport_n2o(:,:) = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2o(:,:) / one_day * dtime 
       transport_n2(:,:)  = - rtm_gasdiffusion_act(:,:) * rmm_gasdiffusion_act(:,:) * n2(:,:)  / one_day * dtime 

       ! emission rates as positive flux to the atmosphere, integrated over the soil column [mol m-2 timestep-1]
       emission_noy(:) = SUM(-transport_noy(:,:) * soil_depth_sl(:,:),DIM=2)
       emission_n2o(:) = SUM(-transport_n2o(:,:) * soil_depth_sl(:,:),DIM=2)
       emission_n2(:)  = SUM(-transport_n2(:,:)  * soil_depth_sl(:,:),DIM=2)

       ! record leaching losses for simple_0d model as [mol m-2 timestep-1] (simple_1d has already been calculated in 2.)
       IF(sb_model_scheme == 'simple_0d') THEN
          leaching_nh4_solute(:) = -transport_nh4_solute(:,1) * soil_depth_sl(:,1)
          leaching_no3_solute(:) = -transport_no3_solute(:,1) * soil_depth_sl(:,1)
          leaching_po4_solute(:) = -transport_po4_solute(:,1) * soil_depth_sl(:,1)
       ENDIF
    END SELECT
 
  END SUBROUTINE calc_actual_nloss_rate


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
  SUBROUTINE calc_litter_partitioning(lctlib, & 
                                      options, &
                                      nsoil, &
                                      soil_depth_sl, &
                                      root_fraction_act_sl, &
                                      veg_litterfall, &
                                      sb_formation )

    USE mo_veg_memory_class,        ONLY: t_slm_struct_veg_pool
    USE mo_veg_constants,           ONLY: itree, igrass
    USE mo_sb_constants,            ONLY: fc_soluable_max,k_fc_soluable,k_fn_soluable,k_fp_soluable, &
                                          lc_leaf_max,lc_leaf2sla,lc_fine_root,lc_coarse_root,lc_fruit,lc_seed_bed 

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),                INTENT(in)    :: lctlib                     !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil                      !< number of soil layers
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_sl, &           !< depth per soil layer
                                                            root_fraction_act_sl       !< root fraction per soil layer
    TYPE(t_slm_struct_veg_pool),           INTENT(in)    :: veg_litterfall             !< vegetation litterfall pool
    TYPE(t_slm_struct_sb_pool),            INTENT(inout) :: sb_formation               !< soil biogeochemistry formation pool
    ! 0.2 Local
    REAL(wp), DIMENSION(options%nc) :: fc_soluable_leaf, fc_soluable_fine_root, fc_soluable_coarse_root, fc_soluable_fruit, & 
                                       fc_soluable_seed_bed, fc_soluable_labile, fc_soluable_reserve 
    REAL(wp), DIMENSION(options%nc) :: fn_soluable_leaf, fn_soluable_fine_root, fn_soluable_coarse_root, fn_soluable_fruit, & 
                                       fn_soluable_seed_bed, fn_soluable_labile, fn_soluable_reserve 
    REAL(wp), DIMENSION(options%nc) :: fp_soluable_leaf, fp_soluable_fine_root, fp_soluable_coarse_root, fp_soluable_fruit, & 
                                       fp_soluable_seed_bed, fp_soluable_labile, fp_soluable_reserve 
    INTEGER                     :: isoil                      !< loop over soil layers
    INTEGER                     :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_litter_partitioning'
    ! ---------------------------
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
   
    !> 1.0 determine allocation of fresh litter fall to litter pools
    !!
    fc_soluable_leaf(:) = 0.0_wp
    fn_soluable_leaf(:) = 0.0_wp
    fp_soluable_leaf(:) = 0.0_wp
    IF(veg_litterfall%leaf%nitrogen%ptr > eps8) THEN
       fc_soluable_leaf(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * (lc_leaf_max + lc_leaf2sla * lctlib%sla) * & 
                                 veg_litterfall%leaf%carbon%ptr / veg_litterfall%leaf%nitrogen%ptr)
       WHERE(fc_soluable_leaf(:)>eps8)
          fn_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fn_soluable*fc_soluable_leaf(:)))
          fp_soluable_leaf(:) = 1._wp / (1._wp + (1._wp - fc_soluable_leaf(:))/(k_fp_soluable*fc_soluable_leaf(:)))
       ENDWHERE
    ENDIF

    fc_soluable_fine_root(:) = 0.0_wp 
    fn_soluable_fine_root(:) = 0.0_wp 
    fp_soluable_fine_root(:) = 0.0_wp 
    IF(veg_litterfall%fine_root%nitrogen%ptr > eps8) THEN
       fc_soluable_fine_root(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fine_root * & 
                                      veg_litterfall%fine_root%carbon%ptr / veg_litterfall%fine_root%nitrogen%ptr)
       WHERE(fc_soluable_fine_root(:)>eps8)
          fn_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fn_soluable*fc_soluable_fine_root(:)))
          fp_soluable_fine_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fine_root(:))/(k_fp_soluable*fc_soluable_fine_root(:)))
       ENDWHERE
    ENDIF

    fc_soluable_coarse_root(:) = 0.0_wp 
    fn_soluable_coarse_root(:) = 0.0_wp 
    fp_soluable_coarse_root(:) = 0.0_wp 
    IF(veg_litterfall%coarse_root%nitrogen%ptr > eps8) THEN
       fc_soluable_coarse_root(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_coarse_root * & 
                                        veg_litterfall%coarse_root%carbon%ptr / veg_litterfall%coarse_root%nitrogen%ptr)
       WHERE(fc_soluable_coarse_root(:)>eps8)
          fn_soluable_coarse_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(:))/(k_fn_soluable*fc_soluable_coarse_root(:)))
          fp_soluable_coarse_root(:) = 1._wp / (1._wp + (1._wp - fc_soluable_coarse_root(:))/(k_fp_soluable*fc_soluable_coarse_root(:)))
        ENDWHERE
    ENDIF

    fc_soluable_fruit(:) = 0.0_wp
    fn_soluable_fruit(:) = 0.0_wp
    fp_soluable_fruit(:) = 0.0_wp
    IF(veg_litterfall%fruit%nitrogen%ptr > eps8) THEN
       fc_soluable_fruit(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_fruit * & 
                                  veg_litterfall%fruit%carbon%ptr / veg_litterfall%fruit%nitrogen%ptr)
       WHERE(fc_soluable_fruit(:)>eps8)
          fn_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fn_soluable*fc_soluable_fruit(:)))
          fp_soluable_fruit(:) = 1._wp / (1._wp + (1._wp - fc_soluable_fruit(:))/(k_fp_soluable*fc_soluable_fruit(:)))
       ENDWHERE
    ENDIF

    fc_soluable_seed_bed(:) = 0.0_wp
    fn_soluable_seed_bed(:) = 0.0_wp
    fp_soluable_seed_bed(:) = 0.0_wp
    IF(veg_litterfall%seed_bed%nitrogen%ptr > eps8) THEN
       fc_soluable_seed_bed(:) = MAX(0.0_wp,fc_soluable_max - k_fc_soluable * lc_seed_bed * & 
                                     veg_litterfall%seed_bed%carbon%ptr / veg_litterfall%seed_bed%nitrogen%ptr)
       WHERE(fc_soluable_seed_bed(:)>eps8)
          fn_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fn_soluable*fc_soluable_seed_bed(:)))
          fp_soluable_seed_bed(:) = 1._wp / (1._wp + (1._wp - fc_soluable_seed_bed(:))/(k_fp_soluable*fc_soluable_seed_bed(:)))
       ENDWHERE
    ENDIF

    fc_soluable_labile(:)  = 1.0_wp 
    fn_soluable_labile(:)  = 1.0_wp 
    fp_soluable_labile(:)  = 1.0_wp 
    fc_soluable_reserve(:) = 1.0_wp
    fn_soluable_reserve(:) = 1.0_wp
    fp_soluable_reserve(:) = 1.0_wp

    !> 2.0 Update surface metabolic and soluable litter pools with litter fall
    !!
    sb_formation%soluable_litter%carbon%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%carbon%ptr(ics:ice,1,iblk) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon%ptr    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon%ptr   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon%ptr  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon%ptr + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon%ptr ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%nitrogen%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%nitrogen%ptr(ics:ice,1,iblk) + & 
       (fn_soluable_leaf(:)     * veg_litterfall%leaf%nitrogen%ptr    + & 
        fn_soluable_fruit(:)    * veg_litterfall%fruit%nitrogen%ptr   + & 
        fn_soluable_labile(:)   * veg_litterfall%labile%nitrogen%ptr  + & 
        fn_soluable_reserve(:)  * veg_litterfall%reserve%nitrogen%ptr + & 
        fn_soluable_seed_bed(:) * veg_litterfall%seed_bed%nitrogen%ptr ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%phosphorus%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
       (fp_soluable_leaf(:)     * veg_litterfall%leaf%phosphorus%ptr    + & 
        fp_soluable_fruit(:)    * veg_litterfall%fruit%phosphorus%ptr   + & 
        fp_soluable_labile(:)   * veg_litterfall%labile%phosphorus%ptr  + & 
        fp_soluable_reserve(:)  * veg_litterfall%reserve%phosphorus%ptr + & 
        fp_soluable_seed_bed(:) * veg_litterfall%seed_bed%phosphorus%ptr ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%carbon13%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%carbon13%ptr(ics:ice,1,iblk) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon13%ptr    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon13%ptr   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon13%ptr  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon13%ptr + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon13%ptr ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%carbon14%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%carbon14%ptr(ics:ice,1,iblk) + & 
       (fc_soluable_leaf(:)     * veg_litterfall%leaf%carbon14%ptr    + & 
        fc_soluable_fruit(:)    * veg_litterfall%fruit%carbon14%ptr   + & 
        fc_soluable_labile(:)   * veg_litterfall%labile%carbon14%ptr  + & 
        fc_soluable_reserve(:)  * veg_litterfall%reserve%carbon14%ptr + & 
        fc_soluable_seed_bed(:) * veg_litterfall%seed_bed%carbon14%ptr ) / soil_depth_sl(:,1)
    sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) + & 
       (fn_soluable_leaf(:)     * veg_litterfall%leaf%nitrogen15%ptr    + & 
        fn_soluable_fruit(:)    * veg_litterfall%fruit%nitrogen15%ptr   + & 
        fn_soluable_labile(:)   * veg_litterfall%labile%nitrogen15%ptr  + & 
        fn_soluable_reserve(:)  * veg_litterfall%reserve%nitrogen15%ptr + & 
        fn_soluable_seed_bed(:) * veg_litterfall%seed_bed%nitrogen15%ptr ) / soil_depth_sl(:,1)

    sb_formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon%ptr    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon%ptr   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon%ptr  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon%ptr + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) + & 
       ((1.0_wp - fn_soluable_leaf(:))     * veg_litterfall%leaf%nitrogen%ptr    + & 
        (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall%fruit%nitrogen%ptr   + & 
        (1.0_wp - fn_soluable_labile(:))   * veg_litterfall%labile%nitrogen%ptr  + & 
        (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall%reserve%nitrogen%ptr + & 
        (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall%seed_bed%nitrogen%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
       ((1.0_wp - fp_soluable_leaf(:))     * veg_litterfall%leaf%phosphorus%ptr    + & 
        (1.0_wp - fp_soluable_fruit(:))    * veg_litterfall%fruit%phosphorus%ptr   + & 
        (1.0_wp - fp_soluable_labile(:))   * veg_litterfall%labile%phosphorus%ptr  + & 
        (1.0_wp - fp_soluable_reserve(:))  * veg_litterfall%reserve%phosphorus%ptr + & 
        (1.0_wp - fp_soluable_seed_bed(:)) * veg_litterfall%seed_bed%phosphorus%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon13%ptr    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon13%ptr   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon13%ptr  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon13%ptr + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon13%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + & 
       ((1.0_wp - fc_soluable_leaf(:))     * veg_litterfall%leaf%carbon14%ptr    + & 
        (1.0_wp - fc_soluable_fruit(:))    * veg_litterfall%fruit%carbon14%ptr   + & 
        (1.0_wp - fc_soluable_labile(:))   * veg_litterfall%labile%carbon14%ptr  + & 
        (1.0_wp - fc_soluable_reserve(:))  * veg_litterfall%reserve%carbon14%ptr + & 
        (1.0_wp - fc_soluable_seed_bed(:)) * veg_litterfall%seed_bed%carbon14%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + & 
       ((1.0_wp - fn_soluable_leaf(:))     * veg_litterfall%leaf%nitrogen15%ptr    + & 
        (1.0_wp - fn_soluable_fruit(:))    * veg_litterfall%fruit%nitrogen15%ptr   + & 
        (1.0_wp - fn_soluable_labile(:))   * veg_litterfall%labile%nitrogen15%ptr  + & 
        (1.0_wp - fn_soluable_reserve(:))  * veg_litterfall%reserve%nitrogen15%ptr + & 
        (1.0_wp - fn_soluable_seed_bed(:)) * veg_litterfall%seed_bed%nitrogen15%ptr ) / soil_depth_sl(:,1)

    !> 3.0 Add surface litter fall from sapwood and heartwood pools to surface woody litter pool for trees
    !!     for grasses, stems and halms (=sapwood) are not actually wood, and therefore paritioned to 
    !!     metabolic and structural litter as if their were fine roots
    SELECT CASE(lctlib%growthform)
       CASE (itree)
          sb_formation%woody_litter%carbon%ptr(ics:ice,1,iblk) = sb_formation%woody_litter%carbon%ptr(ics:ice,1,iblk) + & 
             (veg_litterfall%sap_wood%carbon%ptr + veg_litterfall%heart_wood%carbon%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%nitrogen%ptr(ics:ice,1,iblk) = sb_formation%woody_litter%nitrogen%ptr(ics:ice,1,iblk) +  & 
             (veg_litterfall%sap_wood%nitrogen%ptr + veg_litterfall%heart_wood%nitrogen%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%phosphorus%ptr(ics:ice,1,iblk) = sb_formation%woody_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
             (veg_litterfall%sap_wood%phosphorus%ptr + veg_litterfall%heart_wood%phosphorus%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%carbon13%ptr(ics:ice,1,iblk) = sb_formation%woody_litter%carbon13%ptr(ics:ice,1,iblk) + & 
             (veg_litterfall%sap_wood%carbon13%ptr + veg_litterfall%heart_wood%carbon13%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%carbon14%ptr(ics:ice,1,iblk) = sb_formation%woody_litter%carbon14%ptr(ics:ice,1,iblk) + & 
             (veg_litterfall%sap_wood%carbon14%ptr + veg_litterfall%heart_wood%carbon14%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%woody_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb_formation%woody_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
             (veg_litterfall%sap_wood%nitrogen15%ptr + veg_litterfall%heart_wood%nitrogen15%ptr) / &
             soil_depth_sl(:,1)
       CASE (igrass)
          sb_formation%soluable_litter%carbon%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%carbon%ptr(ics:ice,1,iblk) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon%ptr + veg_litterfall%heart_wood%carbon%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%nitrogen%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%nitrogen%ptr(ics:ice,1,iblk) + & 
             fn_soluable_fine_root(:) * (veg_litterfall%sap_wood%nitrogen%ptr + veg_litterfall%heart_wood%nitrogen%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%phosphorus%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
             fp_soluable_fine_root(:) * (veg_litterfall%sap_wood%phosphorus%ptr + veg_litterfall%heart_wood%phosphorus%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%carbon13%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%carbon13%ptr(ics:ice,1,iblk) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon13%ptr + veg_litterfall%heart_wood%carbon13%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%carbon14%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%carbon14%ptr(ics:ice,1,iblk) + & 
             fc_soluable_fine_root(:) * (veg_litterfall%sap_wood%carbon14%ptr + veg_litterfall%heart_wood%carbon14%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,1,iblk) + & 
             fn_soluable_fine_root(:) * (veg_litterfall%sap_wood%nitrogen15%ptr + veg_litterfall%heart_wood%nitrogen15%ptr) / &
             soil_depth_sl(:,1)
      
          sb_formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * (veg_litterfall%sap_wood%carbon%ptr + veg_litterfall%heart_wood%carbon%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) + & 
             (1.0_wp - fn_soluable_fine_root(:)) * (veg_litterfall%sap_wood%nitrogen%ptr + veg_litterfall%heart_wood%nitrogen%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
             (1.0_wp - fp_soluable_fine_root(:)) * (veg_litterfall%sap_wood%phosphorus%ptr + veg_litterfall%heart_wood%phosphorus%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * (veg_litterfall%sap_wood%carbon13%ptr + veg_litterfall%heart_wood%carbon13%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + & 
             (1.0_wp - fc_soluable_fine_root(:)) * (veg_litterfall%sap_wood%carbon14%ptr + veg_litterfall%heart_wood%carbon14%ptr) / & 
             soil_depth_sl(:,1)
          sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + & 
             (1.0_wp - fn_soluable_fine_root(:)) * (veg_litterfall%sap_wood%nitrogen15%ptr + veg_litterfall%heart_wood%nitrogen15%ptr) / &
             soil_depth_sl(:,1)

    END SELECT

    !> 4.0 Distribute fine and coarse root litter to soil depth according to the root profile 
    !!
    DO isoil = 1,nsoil
       sb_formation%soluable_litter%carbon%ptr(ics:ice,isoil,iblk)     = sb_formation%soluable_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon%ptr     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%nitrogen%ptr(ics:ice,isoil,iblk)   = sb_formation%soluable_litter%nitrogen%ptr(ics:ice,isoil,iblk) + & 
          (fn_soluable_fine_root(:)   * veg_litterfall%fine_root%nitrogen%ptr     + &
           fn_soluable_coarse_root(:) * veg_litterfall%coarse_root%nitrogen%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%phosphorus%ptr(ics:ice,isoil,iblk) = sb_formation%soluable_litter%phosphorus%ptr(ics:ice,isoil,iblk) + & 
          (fp_soluable_fine_root(:)   * veg_litterfall%fine_root%phosphorus%ptr     + &
           fp_soluable_coarse_root(:) * veg_litterfall%coarse_root%phosphorus%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%carbon13%ptr(ics:ice,isoil,iblk)   = sb_formation%soluable_litter%carbon13%ptr(ics:ice,isoil,iblk) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon13%ptr     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon13%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%carbon14%ptr(ics:ice,isoil,iblk)   = sb_formation%soluable_litter%carbon14%ptr(ics:ice,isoil,iblk) + & 
          (fc_soluable_fine_root(:)   * veg_litterfall%fine_root%carbon14%ptr     + &
           fc_soluable_coarse_root(:) * veg_litterfall%coarse_root%carbon14%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = sb_formation%soluable_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + & 
          (fn_soluable_fine_root(:)   * veg_litterfall%fine_root%nitrogen15%ptr     + &
           fn_soluable_coarse_root(:) * veg_litterfall%coarse_root%nitrogen15%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 

       sb_formation%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk)    = sb_formation%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
          ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon%ptr     + &
           (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk)    = sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) + & 
          ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen%ptr     + &
           (1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk)    = sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) + & 
          ((1.0_wp - fp_soluable_fine_root(:))  * veg_litterfall%fine_root%phosphorus%ptr     + &
           (1.0_wp - fp_soluable_coarse_root(:)) * veg_litterfall%coarse_root%phosphorus%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk)    = sb_formation%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) + & 
          ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon13%ptr     + &
           (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon13%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk)    = sb_formation%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) + & 
          ((1.0_wp - fc_soluable_fine_root(:))  * veg_litterfall%fine_root%carbon14%ptr     + &
           (1.0_wp - fc_soluable_coarse_root(:)) * veg_litterfall%coarse_root%carbon14%ptr ) * & 
           root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk)    = sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + & 
          ((1.0_wp - fn_soluable_fine_root(:))  * veg_litterfall%fine_root%nitrogen15%ptr     + &
           (1.0_wp - fn_soluable_coarse_root(:)) * veg_litterfall%coarse_root%nitrogen15%ptr ) * & 
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
  SUBROUTINE calc_mycorrhiza_dynamics(lctlib, &
                                      options, &
                                      nsoil, &
                                      soil_depth_sl, &
                                      mycorrhiza_uptake_nh4_sl, &
                                      mycorrhiza_uptake_no3_sl, &
                                      mycorrhiza_uptake_po4_sl, &
                                      sb_pool, &
                                      sb_formation, &
                                      sb_loss, &
                                      het_respiration, &
                                      het_respiration_c13, &
                                      het_respiration_c14, &
                                      veg_mycorrhiza_export) 

    USE mo_veg_memory_class
    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: mycorrhiza_cn, mycorrhiza_np

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),                INTENT(in)    :: lctlib                       !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil                        !< number of soil layers
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: soil_depth_sl, &             !< depth per soil layer
                                                            mycorrhiza_uptake_nh4_sl, &  !< NH4 uptake by mycorrhiza [mol / m3 / timestep]
                                                            mycorrhiza_uptake_no3_sl, &  !< NO3 uptake by mycorrhiza [mol / m3 / timestep]
                                                            mycorrhiza_uptake_po4_sl     !< PO4 uptake by mycorrhiza [mol / m3 / timestep]
    TYPE(t_slm_struct_sb_pool),            INTENT(in)    :: sb_pool                      !< soil biogeochemistry pools [mol / m3]
    TYPE(t_slm_struct_sb_pool),            INTENT(inout) :: sb_formation, &              !< pool formation rate  [mol / m3 / timestep]
                                                            sb_loss                      !< pool loss rate  [mol / m3 / timestep]
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(out)   :: het_respiration, &           !< heterotrophic respiration [mol / m3 / timestep]
                                                            het_respiration_c13, &       !< heterotrophic respiration [mol / m3 / timestep]
                                                            het_respiration_c14          !< heterotrophic respiration [mol / m3 / timestep]
    TYPE(t_slm_struct_veg_element),        INTENT(inout) :: veg_mycorrhiza_export        !< export from mycorrhiza to plants [mol / m2 / timestep]
    ! 0.2 Local
    INTEGER                               :: isoil                      !< loop over soil layers
    INTEGER                               :: iblk, ics, ice
    REAL(wp), DIMENSION(options%nc,nsoil) :: transform_respiration,fresp,fturn_mycorrhiza
    REAL(wp), DIMENSION(options%nc,nsoil) :: hlp1, hlp2, hlp3
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc__mycorrhiza_dynamics'
    ! ---------------------------
    ! Get local variables from options argument
    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    
    !> 8.0 Mycorrhiza dynamics 
    !! 
    !> 8.1 uptake of nutrients by mycorrhiza and associated respiration losses to heterotrophic respiration
    !!
    ! calculate respiration costs associated with N transformation (mol C / m3 / timestep)
    transform_respiration(:,:) = mycorrhiza_uptake_nh4_sl(:,:) * transform_cost_nh4 + & 
                                 mycorrhiza_uptake_no3_sl(:,:) * transform_cost_no3
    fresp(:,:) = 0.0_wp
    WHERE(sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) > eps8)
       fresp(:,:) = transform_respiration(:,:) / sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE
    sb_loss%mycorrhiza%carbon%ptr(ics:ice,:,iblk)   = sb_loss%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + &
                                                      fresp * sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk)   !! respiration
    sb_loss%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + &
                                                      fresp * sb_pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) !! respiration
    sb_loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + &
                                                      fresp * sb_pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) !! respiration

    het_respiration(:,:)     = het_respiration(:,:)     + fresp * sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk)
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + fresp * sb_pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk)
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + fresp * sb_pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)

    !> 8.2 nutrient flux to vegetation (surplus uptake goes to plants) (mol (N/P) / m2 / timestep)
    !!
!@todo: N15 uptake
 
    ! N demand of mycorrhiza in mol N / m3 / timestep
    hlp1(:,:) = sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) / mycorrhiza_cn - &
                sb_pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk)

    ! constrain flux to minimum of potential uptake and demand
    hlp2(:,:) = MIN(hlp1(:,:), &
                mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:))

    ! maintain C:N of mycorrhizal fungi
    ! N supply for plants = N uptake by mycorrhiza - N demand by mycorrhiza
    ! N demand can be negative => N supply can be greater then current N uptake by plants!
    WHERE(hlp1(:,:) > 0.0_wp)
        ! incorporation of N into mycorrhiza (mol N / m3 / timestep)
        sb_formation%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + &
                                                               hlp2(:,:)
    ELSEWHERE
        !surplus N removal from myorrhiza (mol N / m3 / timestep, note the negative sign
        sb_loss%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) - &
                                                          hlp2(:,:)
    ENDWHERE

    ! remainder is exported to plants (mol N / m2 / timestep)
    hlp3(:,:) = mycorrhiza_uptake_nh4_sl(:,:) + mycorrhiza_uptake_no3_sl(:,:) - hlp2(:,:)

    veg_mycorrhiza_export%nitrogen%ptr = 0.0_wp
    DO isoil = 1, nsoil
         veg_mycorrhiza_export%nitrogen%ptr = veg_mycorrhiza_export%nitrogen%ptr + hlp3(1,isoil) * soil_depth_sl(1,isoil)
    ENDDO 

    ! P demand of mycorrhiza in mol P / m3 / timestep ! BLARPP SZ: constrained to positive solutions
    hlp1(:,:) = MAX(sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) / mycorrhiza_cn / mycorrhiza_np - &
                    sb_pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk), 0.0_wp)
    ! incorporation of P into mycorrhiza (mol P / m3 / timestep)
    ! constrain flux to minimum of potential uptake and demand 
    hlp2(:,:) = MIN(hlp1(:,:), mycorrhiza_uptake_po4_sl(:,:))
    sb_formation%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + &
                                                             hlp2(:,:)

    ! remainder is exported to plants (mol P / m2 / timestep)
    hlp3(:,:) = mycorrhiza_uptake_po4_sl(:,:) - hlp2(:,:)
    veg_mycorrhiza_export%phosphorus%ptr = 0.0_wp
    DO isoil = 1, nsoil
         veg_mycorrhiza_export%phosphorus%ptr = veg_mycorrhiza_export%phosphorus%ptr + hlp3(1,isoil) * soil_depth_sl(1,isoil) 
    ENDDO

    !> 8.3 Formation of litter from mycorrhiza turnover (mol / m3 / timestep)
    !!
    fturn_mycorrhiza(:,:) = 1.0_wp/lctlib%tau_mycorrhiza/one_day/one_year * dtime

    sb_loss%mycorrhiza%carbon%ptr(ics:ice,:,iblk)   = sb_loss%mycorrhiza%carbon%ptr(ics:ice,:,iblk) + &
                                                      fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk)     !! turnover
    sb_loss%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk) + &
                                                      fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk)   !! turnover
    sb_loss%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + &
                                                        fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) !! turnover
    sb_loss%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%carbon13%ptr(ics:ice,:,iblk) + &
                                                      fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk)   !! turnover
    sb_loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%carbon14%ptr(ics:ice,:,iblk) + &
                                                        fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)   !! turnover
    sb_loss%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                        fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk) !! turnover

    sb_formation%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                               fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk)
    sb_formation%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                                 fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon13%ptr(ics:ice,:,iblk)
    sb_formation%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                                 fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%carbon14%ptr(ics:ice,:,iblk)
    sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                                 fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen%ptr(ics:ice,:,iblk)
    sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                                   fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%nitrogen15%ptr(ics:ice,:,iblk)
    sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                                   fturn_mycorrhiza(:,:) * sb_pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk)

  END SUBROUTINE calc_mycorrhiza_dynamics

  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to update_sb_simple_model
  !
  !-----------------------------------------------------------------------------------------------------
  !> Calculates the rate of asymbiontic N fixation, wherever there is a bit of soil carbon 
  !! declines within the soil profile according to the distribution of roots as a proxy for litter input
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_asymbiontic_n_fixation(lctlib, &                                      ! in
                                         options, &
                                         nsoil, &
                                         veg_bnf_scheme, &
                                         nh4_solute, &
                                         no3_solute, &
                                         sb_pool_microbial_carbon, &
                                         root_fraction_act_sl, & 
                                         asymb_n_fixation)                              ! out

    USE mo_sb_constants,       ONLY: soil_n_threshold_bnf

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),                INTENT(in)    :: lctlib                      !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_task_options),              INTENT(in)    :: options
    INTEGER,                               INTENT(in)    :: nsoil                       !< number of soil layers
    CHARACTER(len=*),                      INTENT(in)    :: veg_bnf_scheme              !< veg n fixation scheme  
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(in)    :: nh4_solute, &               !< NH4 concentration [mol m-3] 
                                                            no3_solute, &               !< NO3 concentration [mol m-3]
                                                            sb_pool_microbial_carbon, & !< microbial carbon pool [mol m-3]
                                                            root_fraction_act_sl        !< fraction of roots per root layer [m-3] 
    REAL(wp), DIMENSION(options%nc,nsoil), INTENT(inout) :: asymb_n_fixation            !< asymbiontic rate of N fixation [umol m-3 s-1]
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_asymbiontic_n_fixation'
  
    SELECT CASE(TRIM(veg_bnf_scheme))
      CASE("fixed", "standard")
        WHERE(nh4_solute(:,:) + no3_solute < soil_n_threshold_bnf .AND. & 
          sb_pool_microbial_carbon(:,:) > eps4)
          asymb_n_fixation(:,:) = root_fraction_act_sl(:,:) * lctlib%bnf_base 
        END WHERE
    END SELECT

  END SUBROUTINE calc_asymbiontic_n_fixation

!!!!  !-----------------------------------------------------------------------------------------------------
!!  ! Sub Task
!!  !
!!  !-----------------------------------------------------------------------------------------------------
!!  !> routine to calculate 15N fluxes, given N fluxes calculated by the simple soil model 
!!  !!
!!  !-----------------------------------------------------------------------------------------------------
!!  SUBROUTINE calc_sb_simple_n15_fluxes(tile, options)
!!
!!    USE mo_jsb_class,             ONLY: Get_model
!!    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
!!    USE mo_jsb_model_class,       ONLY: t_jsb_model
!!    USE mo_base_constants,        ONLY: dtime, one_day
!!    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
!!    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
!!    USE mo_soil_config_class,     ONLY: soil_config
!!    USE mo_sb_constants
!!    USE mo_sb_jsm_transport
!!
!!    ! Use of process configurations (t_PROC_config)
!!    dsl4jsb_Use_config(SB_)
!!
!!    ! Use of process memories
!!    dsl4jsb_Use_memory(VEG_)        !  USE mo_veg_memory_class,      ONLY: t_veg_memory
!!    dsl4jsb_Use_memory(SOIL_)       !  USE mo_soil_memory_class,     ONLY: t_soil_memory
!!    dsl4jsb_Use_memory(SB_)         !  USE mo_sb_memory_class,       ONLY: t_sb_memory
!!
!!    IMPLICIT NONE
!!    ! ---------------------------
!!    ! 0.1 InOut
!!    TYPE(t_slm_tile),         INTENT(inout), TARGET    :: tile          !< one tile
!!    TYPE(t_jsb_task_options), INTENT(in)               :: options       !< model options
!!    ! ---------------------------
!!    ! 0.2 Local
!!    TYPE(t_jsb_model),      POINTER       :: model
!!    TYPE(t_lctlib_element), POINTER       :: lctlib               !< land-cover-type library - parameter across pft's
!!    INTEGER                               :: isoil                !< loop over soil_config%nsoil_water
!!    INTEGER                               :: iblk, ics, ice, nc
!!    REAL(wp), DIMENSION(options%nc,soil_config%nsoil_water) ::      &
!!                                             apparent_mic_cue     , &
!!                                             growth_c             , &
!!                                             growth_n             , &
!!                                             hlp1                 , &
!!                                             hlp2                 , &
!!                                             hlp3                 , &
!!                                             mic_uptake_nh4_sl        , &
!!                                             root_fraction_act_sl , &      ! and below needed for switch between 0d and 1d
!!                                             rtm_decomposition_act, &
!!                                             rmm_decomposition_act, & 
!!                                             km1_plant_uptake_n, &
!!                                             km2_plant_uptake_n, &
!!                                             km1_plant_uptake_p, &
!!                                             km2_plant_uptake_p, &
!!                                             km_denitrification_c_act, &
!!                                             km_denitrification_no3_act, &
!!                                             rtm_nitrification_act, &
!!                                             rtm_denitrification_act, &
!!                                             rmm_nitrification_act, &
!!                                             rtm_gasdiffusion_act, &
!!                                             rmm_gasdiffusion_act, &
!!                                             anaerobic_volume_fraction_act
!!    REAL(wp), DIMENSION(options%nc) :: root_depth_avg
!!    
!!    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_pools_lesssimple'
!!    ! ---------------------------
!!    ! 0.3 Memory allocation
!!    ! Declare process configuration and memory Pointers
!!    dsl4jsb_Def_config(SB_)
!!    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
!!    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
!!    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
!!    ! Declare pointers to variables in memory
!!    dsl4jsb_Real3D_onChunk :: percolation_sl
!!    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4
!!    dsl4jsb_Real3D_onChunk :: net_mineralisation_nh4_n15
!!    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3
!!    dsl4jsb_Real3D_onChunk :: net_mineralisation_no3_n15
!!    dsl4jsb_Real3D_onChunk :: nitrification_no3 
!!    dsl4jsb_Real3D_onChunk :: nitrification_noy 
!!    dsl4jsb_Real3D_onChunk :: nitrification_n2o
!!    dsl4jsb_Real3D_onChunk :: denitrification_noy 
!!    dsl4jsb_Real3D_onChunk :: denitrification_n2o 
!!    dsl4jsb_Real3D_onChunk :: denitrification_n2 
!!    dsl4jsb_Real3D_onChunk :: volatilisation_n
!!    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_sl
!!    dsl4jsb_Real3D_onChunk :: plant_uptake_nh4_n15_sl
!!    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_sl
!!    dsl4jsb_Real3D_onChunk :: plant_uptake_no3_n15_sl
!!    dsl4jsb_Real3D_onChunk :: plant_uptake_po4_sl
!!    dsl4jsb_Real3D_onChunk :: transport_nh4_solute
!!    dsl4jsb_Real3D_onChunk :: transport_no3_solute
!!    dsl4jsb_Real3D_onChunk :: transport_po4_solute
!!    dsl4jsb_Real2D_onChunk :: leaching_nh4_solute
!!    dsl4jsb_Real2D_onChunk :: leaching_no3_solute
!!    dsl4jsb_Real2D_onChunk :: leaching_po4_solute
!!    dsl4jsb_Real3D_onChunk :: transport_noy
!!    dsl4jsb_Real3D_onChunk :: transport_n2o
!!    dsl4jsb_Real3D_onChunk :: transport_n2 
!!    dsl4jsb_Real3D_onChunk :: noy
!!    dsl4jsb_Real3D_onChunk :: n2o
!!    dsl4jsb_Real3D_onChunk :: n2 
!!    dsl4jsb_Real3D_onChunk :: nh4_solute 
!!    dsl4jsb_Real3D_onChunk :: no3_solute
!!    dsl4jsb_Real3D_onChunk :: po4_solute 
!!    dsl4jsb_Real3D_onChunk :: nh4_n15_solute 
!!    dsl4jsb_Real3D_onChunk :: no3_n15_solute 
!!    dsl4jsb_Real2D_onChunk :: emission_noy
!!    dsl4jsb_Real2D_onChunk :: emission_n2o
!!    dsl4jsb_Real2D_onChunk :: emission_n2
!!    ! Get local variables from options argument
!!    iblk  = options%iblk
!!    ics   = options%ics
!!    ice   = options%ice
!!    nc    = options%nc
!!    ! ---------------------------
!!    ! 0.4 Process Activity, Debug Option
!!!    IF (.NOT. tile%Is_process_active(xx_)) RETURN
!!!    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
!!    ! ---------------------------
!!    ! 0.5 Get Memory
!!    model  => Get_model(tile%owner_model_id)
!!    lctlib => model%lctlib(tile%lcts%lib_id)
!!    ! Get process config
!!    dsl4jsb_Get_config(SB_)
!!    ! Get process memories
!!    dsl4jsb_Get_memory(SOIL_)
!!    dsl4jsb_Get_memory(SB_)
!!    ! Get process variables (Set pointers to variables in memory)
!!    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)
!!    ! ---------------------------
!!    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_nh4_n15)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       net_mineralisation_no3_n15)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_no3)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_noy)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       nitrification_n2o)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_noy)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2o)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       denitrification_n2)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       volatilisation_n)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_sl)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_nh4_n15_sl)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_sl)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_no3_n15_sl)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       plant_uptake_po4_sl)
!!    dsl4jsb_Get_var3D_onChunk(SB_,       noy) 
!!    dsl4jsb_Get_var3D_onChunk(SB_,       n2o) 
!!    dsl4jsb_Get_var3D_onChunk(SB_,       n2) 
!!    dsl4jsb_Get_var3D_onChunk(SB_,       transport_nh4_solute) 
!!    dsl4jsb_Get_var3D_onChunk(SB_,       transport_no3_solute) 
!!    dsl4jsb_Get_var3D_onChunk(SB_,       transport_po4_solute) 
!!    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_nh4_solute) 
!!    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_no3_solute) 
!!    dsl4jsb_Get_var2D_onChunk(SB_,       leaching_po4_solute) 
!!    dsl4jsb_Get_var3D_onChunk(SB_,       transport_noy) 
!!    dsl4jsb_Get_var3D_onChunk(SB_,       transport_n2o) 
!!
!!
!!
!!
!!    !> 1.0 potential 15N fluxes involving ammonium
!!    !!
!!    IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == 'simple_0d')THEN
!!       ! For the 0D soil model, loss is proportional to the fractional waterloss of the soil column
!!       ! reduction of layer concentration in µmol m-3 s-1
!!       transport_nh4_n15_solute(:,1) = - fleach_nh4 * nh4_n15_solute(:,1) * drainage_fraction(:) *1.e6_wp
!!    ELSE 
!!       IF(soil_config%soil_moisture_scheme=='layer')THEN
!!         ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
!!         ! reduction of layer concentration in µmol m-3 s-1, leaching rate in µmol m-2 s-1
!!         CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
!!                                          fleach_nh4 * nh4_n15_solute(:,:),transport_nh4_n15_solute(:,:),leaching_nh4_n15_solute(:))
!!       ENDIF
!!       transport_nh4_n15_solute(:,:) = transport_nh4_n15_solute(:,:) / dtime * 1.e6_wp
!!       leaching_nh4_n15_solute(:,:) = leaching_nh4_n15_solute(:,:) / dtime * 1.e6_wp
!!    ENDIF
!!    plant_uptake_nh4_n15_sl(:,:) = plant_uptake_nh4_sl(:,:) * & 
!!                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_plant_uptake_nh4) 
!!    mic_uptake_nh4_n15_sl(:,:)   = mic_uptake_nh4_sl(:,:) * & 
!!                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_mic_uptake_nh4)
!!    nitrification_no3_n15(:,:)   = nitrification_no3(:,:) * &
!!                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_nitrate_production)
!!    nitrification_noy_n15(:,:)   = nitrification_noy(:,:) * & 
!!                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_nitrification)
!!    nitrification_n2o_n15(:,:)   = nitrification_n2o(:,:) * &
!!                                   calc_fractionation(nh4_solute(:,:),nh4_n15_solute(:,:),eta_nitrification) 
!!
!!    !> 1.1 contrain 15NH4 fluxes to available 15NH4
!!    !!
!!    hlp1(:,:) = plant_uptake_nh4_n15_sl(:,:) + mic_uptake_nh4_n15_sl(:,:) + nitrification_no3_n15(:,:) + & 
!!                nitrification_noy_n15(:,:) + nitrification_n2o_n15(:,:)
!!    hlp2(:,:) = nh14_n15_solute / dtime * 1.e6_wp
!!    WHERE(hlp2 < hlp1(:,:))
!!       WHERE(hlp1>eps8)
!!          plant_uptake_nh4_n15_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
!!          mic_uptake_nh4_n15_sl(:,:)   = hlp2(:,:) * mic_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
!!          nitrification_no3_n15(:,:)   = hlp2(:,:) * nitrification_no3_n15(:,:) / hlp1(:,:)
!!          nitrification_noy_n15(:,:)   = hlp2(:,:) * nitrification_noy_n15(:,:) / hlp1(:,:)
!!          nitrification_n2o_n15(:,:)   = hlp2(:,:) * nitrification_n2o_n15(:,:) / hlp1(:,:)
!!       ELSEWHERE
!!          plant_uptake_nh4_n15_sl(:,:) = 0.0_wp 
!!          mic_uptake_nh4_n15_sl(:,:)   = 0.0_wp
!!          nitrification_no3_n15(:,:)   = 0.0_wp
!!          nitrification_noy_n15(:,:)   = 0.0_wp
!!          nitrification_n2o_n15(:,:)   = 0.0_wp
!!       ENDWHERE
!!    ENDWHERE
!!
!!    !> 2.0 potential 15N fluxes involving nitrate
!!    !!
!!    IF(TRIM(dsl4jsb_Config(SB_)%sb_model_scheme) == 'simple_0d')THEN
!!       ! For the 0D soil model, loss is proportional to the fractional waterloss of the soil column
!!       ! reduction of layer concentration in mol m-3 timestep-1
!!       transport_nh4_solute(:,1) = - fleach_nh4 * nh4_solute(:,1) * drainage_fraction(:) * dtime
!!       transport_no3_solute(:,1) = - no3_solute(:,1) * drainage_fraction(:) * dtime
!!    ELSE 
!!       IF(soil_config%soil_moisture_scheme=='layer')THEN
!!         ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
!!         ! reduction of layer concentration in mol m-3 timestep-1, leaching rate in mol m-2 timestep-1
!!         CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
!!                                          fleach_nh4 * nh4_solute(:,:),transport_nh4_solute(:,:),leaching_nh4_solute(:))
!!         CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
!!                                          no3_solute(:,:),transport_no3_solute(:,:),leaching_no3_solute(:))
!!       ENDIF
!!    ENDIF
!!
!!    plant_uptake_no3_n15_sl(:,:)    = plant_uptake_no3_sl(:,:) * & 
!!                                      calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_plant_uptake_no3)
!!    mic_uptake_no3_n15_sl(:,:)      = mic_uptake_no3_sl(:,:) * & 
!!                                      calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_mic_uptake_no3) 
!!    denitrification_noy_n15_sl(:,:) = denitrification_noy(:,:) * & 
!!                                      calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_denitrification)
!!    denitrification_n2o_n15_sl(:,:) = denitrification_n2o(:,:) * & 
!!                                      calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_denitrification)
!!    denitrification_n2_n15_sl(:,:)  = denitrification_n2(:,:) * & 
!!                                      calc_fractionation(no3_solute(:,:),no3_n15_solute(:,:),eta_denitrification) 
!!
!!    !> 2.1 contrain 15NO3 fluxes to available 15NO3
!!    !!
!!    hlp1(:,:) = plant_uptake_no3_n15_sl(:,:) + mic_uptake_no3_n15_sl(:,:) + & 
!!                denitrification_noy_n15_sl(:,:) + denitrification_n2o_n15_sl(:,:) + denitrification_n2_n15_sl(:,:)
!!    hlp2(:,:) = no3_n15_solute(:,:) / dtime * 1.e6_wp + transport_no3_n15(:,:)
!!    WHERE(hlp2 < hlp1(:,:))
!!       WHERE(hlp1>eps8)
!!          plant_uptake_nh4_n15_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
!!          mic_uptake_nh4_n15_sl(:,:)   = hlp2(:,:) * mic_uptake_nh4_n15_sl(:,:) / hlp1(:,:) 
!!          nitrification_no3_n15(:,:)   = hlp2(:,:) * nitrification_no3_n15(:,:) / hlp1(:,:)
!!          nitrification_noy_n15(:,:)   = hlp2(:,:) * nitrification_noy_n15(:,:) / hlp1(:,:)
!!          nitrification_n2o_n15(:,:)   = hlp2(:,:) * nitrification_n2o_n15(:,:) / hlp1(:,:)
!!       ELSEWHERE
!!          plant_uptake_nh4_n15_sl(:,:) = 0.0_wp 
!!          mic_uptake_nh4_n15_sl(:,:)   = 0.0_wp
!!          nitrification_no3_n15(:,:)   = 0.0_wp
!!          nitrification_noy_n15(:,:)   = 0.0_wp
!!          nitrification_n2o_n15(:,:)   = 0.0_wp
!!       ENDWHERE
!!    ENDWHERE
!!
!!
!!  END SUBROUTINE calc_sb_simple_n15_fluxes

END MODULE mo_sb_simple_model
