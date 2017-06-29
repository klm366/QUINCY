!>
!! @brief Contains the routines for the Task sb_simple_model 
!!
!! @author
!!  Silvia Caldararu, Sönke Zaehle 
!!
!! @par Revision History
!! First version by Sönke Zaehle (20170105): copyied from vegetation/mo_veg_soil.f90
!! Second version by Sönke Zaehle (20170331): added a 4 pool model 
!! 3rd version: Jan Engel (20170627): split calc_sb_fluxes_simple() into multiple calc_ routines
!!
MODULE mo_sb_simple_model

#include "dsl_jsbach.inc"

  USE mo_kind,            ONLY: wp, eps8, eps4
  USE mo_exception,       ONLY: message, message_text, finish

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
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_sb_simple_model(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    TYPE(t_slm_tile),         INTENT(inout), TARGET  :: tile         !< one tile
    TYPE(t_jsb_task_options), INTENT(in)             :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    !TYPE(t_jsb_model),      POINTER       :: model
    !TYPE(t_lctlib_element), POINTER       :: lctlib     !< land-cover-type library - parameter across pft's
    !INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_sb_simple_model'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! ...





    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    CALL calc_sb_fluxes_simple(tile, options)


!!    CALL calc_n15_fluxes(tile, options)

  END SUBROUTINE update_sb_simple_model


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task
  !
  !-----------------------------------------------------------------------------------------------------
  !> routine to update soil litter and mineral pools (simple version for vegetation model testing)
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_sb_fluxes_simple(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_base_memory,           ONLY: t_slm_tile, t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_base_constants,        ONLY: dtime
    USE mo_jsb_process_class,     ONLY: A2L_, L2A_, SEB_, TURB_, SOIL_, HYDRO_, VEG_, HD_, RAD_, ASSIMI_, PHENO_, SB_
    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_soil_config_class,     ONLY: soil_config
    USE mo_sb_config_class,       ONLY: sb_config

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
                                             growth_c             , &
                                             growth_n             , &
                                             hlp1                 , &
                                             hlp2                 , &
                                             hlp3                 , &
                                             mic_uptake_nh4_sl    , &
                                             fturn_mycorrhiza     , &
                                             transform_respiration , &
                                             fresp                , &
                                             n_demand_mycorrhiza  , &
                                             n2mycorrhiza         , &
                                             n2plants             , &
                                             mic_uptake_po4_sl    , &
                                             biomineralisation_po4, &
                                             root_fraction_act_sl , &      ! and below needed for switch between 0d and 1d
                                             rtm_depolymerisation_act, &
                                             rmm_depolymerisation_act, & 
                                             km1_plant_uptake_n, &
                                             km2_plant_uptake_n, &
                                             km1_plant_uptake_p, &
                                             km2_plant_uptake_p, &
                                             km_denitrification_c_act, &
                                             km_denitrification_no3_act, &
                                             rtm_nitrification_act, &
                                             rtm_denitrification_act, &
                                             rmm_nitrification_act, &
                                             rtm_gasdiffusion_act, &
                                             rmm_gasdiffusion_act, &
                                             anaerobic_volume_fraction_act
    REAL(wp), DIMENSION(options%nc)       :: root_depth_avg
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_sb_fluxes_simple'
    ! ---------------------------
    ! 0.3 Memory allocation
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(VEG_)                          ! TYPE(t_veg_memory),      POINTER :: veg__mem
    dsl4jsb_Def_memory(SOIL_)                         ! TYPE(t_soil_memory),     POINTER :: soil__mem
    dsl4jsb_Def_memory(SB_)                           ! TYPE(t_sb_memory),       POINTER :: sb__mem
    ! Declare pointers to variables in memory
    ! VEG_
    dsl4jsb_Real2D_onChunk :: uptake_norg
    dsl4jsb_Real2D_onChunk :: uptake_norg_n15
    dsl4jsb_Real2D_onChunk :: uptake_porg
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
    dsl4jsb_Real3D_onChunk :: rtm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rtm_nitrification
    dsl4jsb_Real3D_onChunk :: rtm_denitrification
    dsl4jsb_Real3D_onChunk :: rtm_gasdiffusion
    dsl4jsb_Real3D_onChunk :: rmm_depolymerisation 
    dsl4jsb_Real3D_onChunk :: rmm_mic_uptake
    dsl4jsb_Real3D_onChunk :: rmm_sorption
    dsl4jsb_Real3D_onChunk :: rmm_desorption
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
    ! Get process config
!    dsl4jsb_Get_config(xx_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SOIL_)
    dsl4jsb_Get_memory(SB_)
    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    lctlib => model%lctlib(tile%lcts%lib_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_,      uptake_norg)
    dsl4jsb_Get_var2D_onChunk(VEG_,      uptake_norg_n15)
    dsl4jsb_Get_var2D_onChunk(VEG_,      uptake_porg)
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
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_depolymerisation)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_nitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_denitrification)
    dsl4jsb_Get_var3D_onChunk(SB_,       rtm_gasdiffusion)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_depolymerisation)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_mic_uptake)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_sorption)
    dsl4jsb_Get_var3D_onChunk(SB_,       rmm_desorption)
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


!! *slmtemp
    veg  => veg__mem
    sb   => sb__mem

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.9 
    !! 
    CALL calc_vertical_profile(lctlib, &
                               nc, &
                               soil_config%nsoil_water, &
                               root_fraction_act_sl(:,:), &
                               root_depth_avg(:), &
                               root_fraction_sl(:,:), &
                               soil_depth_sl(:,:), &
                               km1_plant_uptake_n(:,:), &
                               km1_plant_uptake_p(:,:), &
                               km2_plant_uptake_n(:,:), &
                               km2_plant_uptake_p(:,:), &
                               km_denitrification_c_act(:,:), &
                               km_denitrification_no3_act(:,:), &
                               rtm_depolymerisation(:,:), &
                               rmm_depolymerisation(:,:), &
                               rtm_nitrification(:,:), &
                               rmm_nitrification(:,:), &
                               rtm_denitrification(:,:), &
                               rtm_gasdiffusion(:,:), &
                               rmm_gasdiffusion(:,:), &
                               anaerobic_volume_fraction_sl(:,:), &
                               rtm_depolymerisation_act(:,:), &
                               rmm_depolymerisation_act(:,:), &
                               rtm_nitrification_act(:,:), &
                               rmm_nitrification_act(:,:), &
                               rtm_denitrification_act(:,:), &
                               rtm_gasdiffusion_act(:,:), &
                               rmm_gasdiffusion_act(:,:), &
                               anaerobic_volume_fraction_act(:,:) )

    !> 1.0 Plant nutrient uptake rates
    !!
    CALL calc_potential_plant_uptake_rate(lctlib, &
                                          nc, &
                                          soil_config%nsoil_water, &
                                          veg%pool%fine_root%carbon%ptr, &
                                          sb%pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk), &
                                          root_fraction_act_sl(:,:), &
                                          soil_depth_sl(:,:), &
                                          nh4_solute(:,:), &
                                          no3_solute(:,:), &
                                          po4_solute(:,:), &
                                          km1_plant_uptake_n(:,:), &
                                          km2_plant_uptake_n(:,:), &
                                          km1_plant_uptake_p(:,:), &
                                          km2_plant_uptake_p(:,:), &
                                          plant_uptake_nh4_sl(:,:), &
                                          plant_uptake_no3_sl(:,:), &
                                          plant_uptake_po4_sl(:,:), &
                                          mycorrhiza_uptake_nh4_sl(:,:), &
                                          mycorrhiza_uptake_no3_sl(:,:), &
                                          mycorrhiza_uptake_po4_sl(:,:), &
                                          unit_uptake_n_pot(:), &
                                          unit_uptake_p_pot(:) )

    !> 2.0 Loss pathways for the dynamic loss calculation scheme 
    !!
    CALL calc_potential_nloss_rate(options, &
                                   nc, &
                                   soil_config%nsoil_water, &
                                   sb%pool%microbial%carbon%ptr(ics:ice,:,iblk), &
                                   nitrification_no3(:,:), &
                                   rtm_nitrification_act(:,:), &
                                   rmm_nitrification_act(:,:), &
                                   denitrification_n2(:,:), &
                                   rtm_denitrification_act(:,:), &
                                   km_denitrification_c_act(:,:), &
                                   km_denitrification_no3_act(:,:), &
                                   anaerobic_volume_fraction_act(:,:), &
                                   nh4_solute(:,:), &
                                   no3_solute(:,:), &
                                   po4_solute(:,:), &
                                   transport_nh4_solute(:,:), &
                                   transport_no3_solute(:,:), &
                                   transport_po4_solute(:,:), &
                                   percolation_sl(:,:), &
                                   soil_depth_sl(:,:), &
                                   drainage_fraction(:), &
                                   leaching_nh4_solute(:), &
                                   leaching_no3_solute(:), &
                                   leaching_po4_solute(:) )  

    !> 3.0 Calculate decomposition of pools: dP/dt = - f(T)*f(M)/tau_pool * P
    !!
    CALL calc_potential_decomposition_rate(ics, &
                                           ice, &
                                           nc, &
                                           iblk, &
                                           soil_config%nsoil_water, &
                                           rtm_depolymerisation_act(:,:), &
                                           rmm_depolymerisation_act(:,:), &
                                           sb%loss, &
                                           sb%pool )

    !> 4.0 formation of SOM resulting from the decomposition dPrec/dt = CUE_Pdon * dPdon/dt
    !! net mineralisation follows the stoichiometric difference of the donor and recipient pool, taking account for CUE
    !! for fast -> slow -> fast this only causes positive mineralisation because of the C:N:P and CUE
    !! for litter -> fast this causes in most circumstances negative mineralisation, which is then taken from the 
    !!               mineral pools (or organic matter for P)
    CALL calc_actual_flux_rate(ics, &
                               ice, &
                               nc, &
                               iblk, &
                               soil_config%nsoil_water, &
                               sb%pool%residue%phosphorus%ptr(ics:ice,:,iblk), &
                               growth_c(:,:), &
                               growth_n(:,:), &
                               mic_uptake_nh4_sl(:,:), &
                               plant_uptake_nh4_sl(:,:), &
                               mycorrhiza_uptake_nh4_sl(:,:), &
                               plant_uptake_no3_sl(:,:), &
                               mycorrhiza_uptake_no3_sl(:,:), &
                               mic_uptake_po4_sl(:,:), &
                               plant_uptake_po4_sl(:,:), &
                               mycorrhiza_uptake_po4_sl(:,:), &
                               nitrification_no3(:,:), &
                               denitrification_n2(:,:), &
                               nh4_solute(:,:), &
                               no3_solute(:,:), &
                               po4_solute(:,:), &
                               transport_nh4_solute(:,:), &
                               transport_no3_solute(:,:), &
                               transport_po4_solute(:,:), &
                               apparent_mic_cue(:,:), &
                               biomineralisation_po4(:,:), &
                               rtm_depolymerisation_act(:,:), &
                               rmm_depolymerisation_act(:,:), &
                               sb%loss )
    
    !> 5.0 formation of SOM and associated respiration, net mineralisation fluxes
    !!          
    CALL calc_actual_SOM_formation(ics, &
                                   ice, &
                                   nc, &
                                   iblk, &
                                   soil_config%nsoil_water, &
                                   apparent_mic_cue(:,:), &
                                   mic_uptake_nh4_sl(:,:), &
                                   mic_uptake_po4_sl(:,:), &
                                   het_respiration(:,:), &
                                   het_respiration_c13(:,:), &
                                   het_respiration_c14(:,:), &
                                   net_mineralisation_nh4(:,:), &
                                   net_mineralisation_po4(:,:), &
                                   net_mineralisation_nh4_n15(:,:), &
                                   sb%loss, &
                                   sb%formation )

    !> 6.0 calculate N losses and nitrification
    !! 
    CALL calc_actual_nloss_rate(nc, &
                                soil_config%nsoil_water, &
                                soil_depth_sl(:,:), &
                                net_mineralisation_nh4(:,:), &
                                volatilisation_n(:,:), &
                                nitrification_no3(:,:), &
                                nitrification_noy(:,:), &
                                nitrification_n2o(:,:), &
                                denitrification_noy(:,:), &
                                denitrification_n2o(:,:), &
                                denitrification_n2(:,:), &
                                noy(:,:), &
                                n2o(:,:), &
                                n2(:,:), &
                                transport_noy(:,:), &
                                transport_n2o(:,:), &
                                transport_n2(:,:), &
                                transport_nh4_solute(:,:), &
                                transport_no3_solute(:,:), &
                                transport_po4_solute(:,:), &
                                rtm_denitrification_act(:,:), &
                                rtm_gasdiffusion_act(:,:), &
                                rmm_gasdiffusion_act(:,:), &
                                emission_noy(:), &
                                emission_n2o(:), &
                                emission_n2(:), &
                                leaching_nh4_solute(:), &
                                leaching_no3_solute(:), &
                                leaching_po4_solute(:) )

    !> 7.0 Update soil litter pools with litter fall
    !!
    CALL calc_litter_partitioning(ics, &
                                  ice, &
                                  nc, &
                                  iblk, &
                                  soil_config%nsoil_water, &
                                  soil_depth_sl(:,:), &
                                  root_fraction_act_sl(:,:), &
                                  veg%litterfall, &
                                  sb%formation )

    !> 8.0 Mycorrhiza dynamics 
    !! 
    CALL calc_mycorrhiza_dynamics(lctlib, &
                                  ics, &
                                  ice, &
                                  nc, &
                                  iblk, &
                                  soil_config%nsoil_water, &
                                  soil_depth_sl(:,:), &
                                  mycorrhiza_uptake_nh4_sl(:,:), &
                                  mycorrhiza_uptake_no3_sl(:,:), &
                                  mycorrhiza_uptake_po4_sl(:,:), &
                                  fresp(:,:), &
                                  transform_respiration(:,:), &
                                  het_respiration(:,:), &
                                  het_respiration_c13(:,:), &
                                  het_respiration_c14(:,:), &
                                  fturn_mycorrhiza(:,:), &
                                  uptake_norg(:), &
                                  uptake_porg(:), &
                                  sb%formation, &
                                  sb%pool, &
                                  sb%loss )

    !> 9.0 convert rates to µmol/ m3 / s (µmol / m2 / s for surface fluxes)
    !!
    het_respiration(:,:) = het_respiration(:,:) / dtime * 1000000._wp
    het_respiration_c13(:,:) = het_respiration_c13(:,:) / dtime * 1000000._wp
    het_respiration_c14(:,:) = het_respiration_c14(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) / dtime * 1000000._wp
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) / dtime * 1000000._wp
    net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) / dtime * 1000000._wp
    plant_uptake_nh4_sl(:,:) = plant_uptake_nh4_sl(:,:) / dtime * 1000000._wp 
    plant_uptake_no3_sl(:,:) = plant_uptake_no3_sl(:,:) / dtime * 1000000._wp 
    uptake_norg(:) = uptake_norg(:) / dtime * 1000000._wp
    plant_uptake_po4_sl(:,:) = plant_uptake_po4_sl(:,:) / dtime * 1000000._wp 
    uptake_porg(:) = uptake_porg(:) / dtime * 1000000._wp
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

  END SUBROUTINE calc_sb_fluxes_simple 


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_vertical_profile  differentiate between simple_0d and simple_1d (1 and many soil layers)
  !!
  !! for the 0D model, all processes happen only in the first layer of the soil \n
  !! therefore the soil temperature and moisture modifiers are aggregated 
  !! across the profile (assuming that the weighting follows the prescribed root profile) \n
  !! 
  !! The km's are expressed as mol/m3 or m3/mol and therefore need to be adjusted to avoid
  !! the implicit assumption that all material would be concentrated in the top layer 
  !! this is done here by assuming the material is distributed within the average rooting depth
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_vertical_profile(lctlib, &
                                    nc, &
                                    nsoil_water, &
                                    root_fraction_act_sl, &
                                    root_depth_avg, &
                                    root_fraction_sl, &
                                    soil_depth_sl, &
                                    km1_plant_uptake_n, &
                                    km1_plant_uptake_p, &
                                    km2_plant_uptake_n, &
                                    km2_plant_uptake_p, &
                                    km_denitrification_c_act, &
                                    km_denitrification_no3_act, &
                                    rtm_depolymerisation, &
                                    rmm_depolymerisation, &
                                    rtm_nitrification, &
                                    rmm_nitrification, &
                                    rtm_denitrification, &
                                    rtm_gasdiffusion, &
                                    rmm_gasdiffusion, &
                                    anaerobic_volume_fraction_sl, &
                                    rtm_depolymerisation_act, &
                                    rmm_depolymerisation_act, &
                                    rtm_nitrification_act, &
                                    rmm_nitrification_act, &
                                    rtm_denitrification_act, &
                                    rtm_gasdiffusion_act, &
                                    rmm_gasdiffusion_act, &
                                    anaerobic_volume_fraction_act )

    USE mo_sb_config_class,         ONLY: sb_config
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element    
    USE mo_sb_constants,            ONLY: km_denitrification_c, km_denitrification_no3

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),               INTENT(in)    :: lctlib                     !< land-cover-type library - parameter across pft's
    INTEGER,                              INTENT(in)    :: nc                         !< number of chunks [ics:ice]
    INTEGER,                              INTENT(in)    :: nsoil_water                !< number of soil layers
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: root_fraction_act_sl, &
                                                           root_fraction_sl, &
                                                           soil_depth_sl, &
                                                           km1_plant_uptake_n, &
                                                           km1_plant_uptake_p, &
                                                           km2_plant_uptake_n, &
                                                           km2_plant_uptake_p, &
                                                           km_denitrification_c_act, &
                                                           km_denitrification_no3_act, &
                                                           rtm_depolymerisation, &
                                                           rmm_depolymerisation, &
                                                           rtm_nitrification, &
                                                           rmm_nitrification, &
                                                           rtm_denitrification, &
                                                           rtm_gasdiffusion, &
                                                           rmm_gasdiffusion, &
                                                           anaerobic_volume_fraction_sl, &
                                                           rtm_depolymerisation_act, &
                                                           rmm_depolymerisation_act, &
                                                           rtm_nitrification_act, &
                                                           rmm_nitrification_act, &
                                                           rtm_denitrification_act, &
                                                           rtm_gasdiffusion_act, &
                                                           rmm_gasdiffusion_act, &
                                                           anaerobic_volume_fraction_act
    REAL(wp), DIMENSION(nc),              INTENT(inout) :: root_depth_avg
    ! 0.2 Local
    INTEGER                               :: isoil                !< loop over soil_config%nsoil_water


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
    

    IF(sb_config%sb_model_scheme == 'simple_0d') THEN
      root_fraction_act_sl(:,1)                         = 1.0_wp
      root_fraction_act_sl(:,2:nsoil_water) = 0.0_wp
      ! set to default values to avoid division by zero below (2:nsoil_water are unused in 0D)
      km1_plant_uptake_n(:,:) = lctlib%soil_n_aff1 
      km1_plant_uptake_p(:,:) = lctlib%soil_p_aff1 
      km2_plant_uptake_n(:,:) = lctlib%soil_n_aff2 
      km2_plant_uptake_p(:,:) = lctlib%soil_p_aff2
      km_denitrification_c_act(:,:) = km_denitrification_c
      km_denitrification_no3_act(:,:) = km_denitrification_no3
      root_depth_avg(:) = sum(root_fraction_sl(:,:)*soil_depth_sl(:,:),DIM=2)
      WHERE(root_depth_avg(:)>eps8) 
         km1_plant_uptake_n(:,1) = lctlib%soil_n_aff1 * soil_depth_sl(:,1) / root_depth_avg(:)
         km1_plant_uptake_p(:,1) = lctlib%soil_p_aff1 * soil_depth_sl(:,1) / root_depth_avg(:)
         km2_plant_uptake_n(:,1) = lctlib%soil_n_aff2 / soil_depth_sl(:,1) * root_depth_avg(:)
         km2_plant_uptake_p(:,1) = lctlib%soil_p_aff2 / soil_depth_sl(:,1) * root_depth_avg(:)
         km_denitrification_c_act(:,1) = km_denitrification_c / soil_depth_sl(:,1) * root_depth_avg(:)
         km_denitrification_no3_act(:,1) = km_denitrification_no3 / soil_depth_sl(:,1) * root_depth_avg(:)
      ENDWHERE
      rtm_depolymerisation_act(:,:) = 0.0_wp 
      rmm_depolymerisation_act(:,:) = 0.0_wp 
      rtm_nitrification_act(:,:) = 0.0_wp 
      rmm_nitrification_act(:,:) = 0.0_wp 
      rtm_denitrification_act(:,:) = 0.0_wp 
      rtm_gasdiffusion_act(:,:) = 0.0_wp 
      rmm_gasdiffusion_act(:,:) = 0.0_wp 
      anaerobic_volume_fraction_act(:,:) = 0.0_wp 
      DO isoil = 1,nsoil_water
         rtm_depolymerisation_act(:,1) = rtm_depolymerisation_act(:,1) + & 
                                         rtm_depolymerisation(:,isoil) * root_fraction_sl(:,isoil)
         rmm_depolymerisation_act(:,1) = rmm_depolymerisation_act(:,1) + & 
                                         rmm_depolymerisation(:,isoil) * root_fraction_sl(:,isoil)
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
    ELSE ! simple_1d, simply copy values
      root_fraction_act_sl(:,:) = root_fraction_sl(:,:)
      rtm_depolymerisation_act(:,:) = rtm_depolymerisation(:,:)
      rmm_depolymerisation_act(:,:) = rmm_depolymerisation(:,:) 
      rtm_nitrification_act(:,:) = rtm_nitrification(:,:)
      rmm_nitrification_act(:,:) = rmm_nitrification(:,:) 
      rtm_denitrification_act(:,:) = rtm_denitrification(:,:)
      rtm_gasdiffusion_act(:,:) = rtm_gasdiffusion(:,:)
      rmm_gasdiffusion_act(:,:) = rmm_gasdiffusion(:,:) 
      anaerobic_volume_fraction_act(:,:) = anaerobic_volume_fraction_sl(:,:) 
      km1_plant_uptake_n(:,:) = lctlib%soil_n_aff1 
      km1_plant_uptake_p(:,:) = lctlib%soil_p_aff1 
      km2_plant_uptake_n(:,:) = lctlib%soil_n_aff2 
      km2_plant_uptake_p(:,:) = lctlib%soil_p_aff2 
      km_denitrification_c_act(:,:) = km_denitrification_c
      km_denitrification_no3_act(:,:) = km_denitrification_no3
    ENDIF

  END SUBROUTINE calc_vertical_profile


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_potential_plant_uptake_rate, Plant nutrient uptake rates
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_plant_uptake_rate(lctlib, &
                                              nc, &
                                              nsoil_water, &
                                              veg_pool_fine_root_carbon, &
                                              sb_pool_mycorrhiza_carbon, &
                                              root_fraction_act_sl, &
                                              soil_depth_sl, &
                                              nh4_solute, &
                                              no3_solute, &
                                              po4_solute, &
                                              km1_plant_uptake_n, &
                                              km2_plant_uptake_n, &
                                              km1_plant_uptake_p, &
                                              km2_plant_uptake_p, &
                                              plant_uptake_nh4_sl, &
                                              plant_uptake_no3_sl, &
                                              plant_uptake_po4_sl, &
                                              mycorrhiza_uptake_nh4_sl, &
                                              mycorrhiza_uptake_no3_sl, &
                                              mycorrhiza_uptake_po4_sl, &
                                              unit_uptake_n_pot, &
                                              unit_uptake_p_pot  )
  
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element    
    USE mo_base_constants,          ONLY: dtime, one_day
    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),               INTENT(in)    :: lctlib                     !< land-cover-type library - parameter across pft's
    INTEGER,                              INTENT(in)    :: nc                         !< number of chunks [ics:ice]
    INTEGER,                              INTENT(in)    :: nsoil_water                !< number of soil layers
    REAL(wp),                             INTENT(in)    :: veg_pool_fine_root_carbon
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(in)    :: sb_pool_mycorrhiza_carbon
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: root_fraction_act_sl, &
                                                           soil_depth_sl, &
                                                           nh4_solute, &
                                                           no3_solute, &
                                                           po4_solute, &
                                                           km1_plant_uptake_n, &
                                                           km2_plant_uptake_n, &
                                                           km1_plant_uptake_p, &
                                                           km2_plant_uptake_p, &
                                                           plant_uptake_nh4_sl, &
                                                           plant_uptake_no3_sl, &
                                                           plant_uptake_po4_sl, &
                                                           mycorrhiza_uptake_nh4_sl, &
                                                           mycorrhiza_uptake_no3_sl, &
                                                           mycorrhiza_uptake_po4_sl
    REAL(wp), DIMENSION(nc),              INTENT(inout) :: unit_uptake_n_pot, &
                                                           unit_uptake_p_pot
    ! 0.2 Local
    INTEGER                               :: isoil                !< loop over soil_config%nsoil_water
    REAL(wp), DIMENSION(nc,nsoil_water)   :: hlp1, hlp2

    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
    

    !> 1.0 Plant nutrient uptake rates
    !!
    unit_uptake_n_pot(:) = 0.0_wp
    unit_uptake_p_pot(:) = 0.0_wp

    DO isoil = 1,nsoil_water
      !> 1.1 inorganic plant N uptake 
      !!
      !! potential N uptake rate (µmol N / mol C fine-root / s), weighted by root distribution
      hlp1(:,isoil) = root_fraction_act_sl(:,isoil) * lctlib%vmax_uptake_n * & 
                  nh4_solute(:,isoil) * ( km1_plant_uptake_n(:,isoil) + & 
                  1._wp / ( km2_plant_uptake_n(:,isoil) + nh4_solute(:,isoil) )) 
      hlp2(:,isoil) = root_fraction_act_sl(:,isoil) * lctlib%vmax_uptake_n * & 
                  no3_solute(:,isoil) * ( km1_plant_uptake_n(:,isoil) + & 
                  1._wp / ( km2_plant_uptake_n(:,isoil) + no3_solute(:,isoil) )) 
      unit_uptake_n_pot(:) = unit_uptake_n_pot(:) + hlp1(:,isoil) + hlp2(:,isoil) 
  
      ! potential uptake rate in mol N / m3 / timiestep
      plant_uptake_nh4_sl(:,isoil) = hlp1(:,isoil) * veg_pool_fine_root_carbon / & 
                                     soil_depth_sl(:,isoil) * dtime/1.e6_wp
      plant_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * veg_pool_fine_root_carbon / & 
                                     soil_depth_sl(:,isoil) * dtime/1.e6_wp
  
      !> 1.2 inorganic plant P uptake 
      !!
      !! potential P uptake rate (µmol P / mol C fine-root / s), weighted by root distribution
      hlp1(:,isoil) = root_fraction_act_sl(:,isoil) * lctlib%vmax_uptake_p * &
                             po4_solute(:,isoil) * ( km1_plant_uptake_p(:,isoil) + &
                             1._wp / ( km2_plant_uptake_p(:,isoil) + po4_solute(:,isoil) ))
      unit_uptake_p_pot(:) = unit_uptake_p_pot(:) + hlp1(:,isoil)
      ! potential uptake rate in mol P / m3 / timestep
      plant_uptake_po4_sl(:,isoil) = hlp1(:,isoil) * & 
                                     veg_pool_fine_root_carbon / soil_depth_sl(:,isoil) * dtime/1.e6_wp
   
      !> 1.3 pot. N uptake by mycorrhizal fungi (mol / m3 / timestep)
      !!
      mycorrhiza_uptake_nh4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          10.0_wp * lctlib%vmax_uptake_n * & 
                                          nh4_solute(:,isoil) * ( km1_plant_uptake_n(:,isoil) + & 
                                          1._wp / ( km2_plant_uptake_n(:,isoil) + nh4_solute(:,isoil) )) * & 
                                          dtime/1.e6_wp
      mycorrhiza_uptake_no3_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          10.0_wp * lctlib%vmax_uptake_n * & 
                                          no3_solute(:,isoil) * ( km1_plant_uptake_n(:,isoil) + & 
                                          1._wp / ( km2_plant_uptake_n(:,isoil) + no3_solute(:,isoil) )) * & 
                                          dtime/1.e6_wp
      hlp1(:,isoil) = mycorrhiza_uptake_nh4_sl(:,isoil) * transform_cost_nh4 + & 
                      mycorrhiza_uptake_no3_sl(:,isoil) * transform_cost_no3 

      WHERE(hlp1(:,isoil)  > sb_pool_mycorrhiza_carbon(:,isoil) * 0.4_wp / one_day * dtime )
         WHERE(hlp1(:,isoil) > eps8)
            hlp2(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * 0.4_wp / one_day * dtime / & 
                            (hlp1(:,isoil)) 
            mycorrhiza_uptake_nh4_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_nh4_sl(:,isoil)
            mycorrhiza_uptake_no3_sl(:,isoil) = hlp2(:,isoil) * mycorrhiza_uptake_no3_sl(:,isoil)
         ELSEWHERE
            mycorrhiza_uptake_nh4_sl(:,isoil) = 0.0_wp
            mycorrhiza_uptake_no3_sl(:,isoil) = 0.0_wp
         ENDWHERE
      ENDWHERE 

      mycorrhiza_uptake_po4_sl(:,isoil) = sb_pool_mycorrhiza_carbon(:,isoil) * &
                                          10.0_wp * lctlib%vmax_uptake_p * & 
                                          po4_solute(:,isoil) * ( km1_plant_uptake_p(:,isoil) + & 
                                          1._wp / ( km2_plant_uptake_p(:,isoil) + po4_solute(:,isoil) )) * & 
                                          dtime/1.e6_wp
    END DO    
  
  END SUBROUTINE calc_potential_plant_uptake_rate
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_potential_nloss_rate, Loss pathways for the dynamic loss calculation scheme
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_nloss_rate(options, &
                                       nc, &
                                       nsoil_water, &
                                       sb_pool_microbial_carbon, &
                                       nitrification_no3, &
                                       rtm_nitrification_act, &
                                       rmm_nitrification_act, &
                                       denitrification_n2, &
                                       rtm_denitrification_act, &
                                       km_denitrification_c_act, &
                                       km_denitrification_no3_act, &
                                       anaerobic_volume_fraction_act, &
                                       nh4_solute, &
                                       no3_solute, &
                                       po4_solute, &
                                       transport_nh4_solute, &
                                       transport_no3_solute, &
                                       transport_po4_solute, &
                                       drainage_fraction, &
                                       percolation_sl, &
                                       soil_depth_sl, &
                                       leaching_nh4_solute, &
                                       leaching_no3_solute, &
                                       leaching_po4_solute  )

    USE mo_base_memory,             ONLY: t_jsb_task_options
    USE mo_sb_config_class,         ONLY: sb_config
    USE mo_soil_config_class,       ONLY: soil_config
    USE mo_base_constants,          ONLY: dtime
    USE mo_sb_constants,            ONLY: vmax_nitrification, vmax_denitrification, fleach_nh4, fleach_po4
    USE mo_sb_jsm_transport,        ONLY: calc_liquid_phase_transport

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_jsb_task_options),             INTENT(in)    :: options
    INTEGER,                              INTENT(in)    :: nc                         !< number of chunks [ics:ice]
    INTEGER,                              INTENT(in)    :: nsoil_water                !< number of soil layers
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(in)    :: sb_pool_microbial_carbon
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: nitrification_no3, &
                                                           rtm_nitrification_act, &
                                                           rmm_nitrification_act, &
                                                           denitrification_n2, &
                                                           rtm_denitrification_act, &
                                                           km_denitrification_c_act, &
                                                           km_denitrification_no3_act, &
                                                           anaerobic_volume_fraction_act, &
                                                           nh4_solute, &
                                                           no3_solute, &
                                                           po4_solute, &
                                                           transport_nh4_solute, &
                                                           transport_no3_solute, &
                                                           transport_po4_solute, &
                                                           percolation_sl, &
                                                           soil_depth_sl
    REAL(wp), DIMENSION(nc),              INTENT(inout) :: drainage_fraction, &
                                                           leaching_nh4_solute, &
                                                           leaching_no3_solute, &
                                                           leaching_po4_solute
    ! 0.2 Local
    ! 


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 2.0 Loss pathways for the dynamic loss calculation scheme 
    !!
    IF(sb_config%sb_nloss_scheme == 'dynamic') THEN
       !> 2.1 Potential nitrification and denitrification rates (according to Zaehle et al. 2011)
       !! potential nitrification rate from NH4 to NO3 [mol/m3/timestep]
       nitrification_no3(:,:)   = vmax_nitrification * rtm_nitrification_act(:,:) * rmm_nitrification_act(:,:) * & 
                                  (1._wp - anaerobic_volume_fraction_act(:,:)) * nh4_solute(:,:) * dtime/1.e6_wp 
       !! potential denitrification rate from NO3 to N2 [mol/m3/timestep]
       denitrification_n2(:,:)  = vmax_denitrification * rtm_denitrification_act(:,:) * &
                                  sb_pool_microbial_carbon(:,:) / &
                                  ( km_denitrification_c_act(:,:) + sb_pool_microbial_carbon(:,:))  * &
                                  anaerobic_volume_fraction_act(:,:) * no3_solute(:,:) / &
                                  ( km_denitrification_no3_act(:,:) + no3_solute(:,:)) * dtime/1.e6_wp 

       !> 2.2. Leaching rate
       !!
       IF(sb_config%sb_model_scheme == 'simple_0d')THEN
          ! For the 0D soil model, loss is proportional to the fractional waterloss of the soil column
          ! reduction of layer concentration in mol m-3 timestep-1
          transport_nh4_solute(:,1) = - fleach_nh4 * nh4_solute(:,1) * drainage_fraction(:) * dtime
          transport_no3_solute(:,1) = - no3_solute(:,1) * drainage_fraction(:) * dtime
          transport_po4_solute(:,1) = - fleach_po4 * po4_solute(:,1) * drainage_fraction(:) * dtime

       ELSE 
          IF(soil_config%soil_moisture_scheme=='layer')THEN
            ! For the 1D scheme with multiple soil hydrology layers, use adjection-diffusion scheme of JSM
            ! reduction of layer concentration in mol m-3 timestep-1, leaching rate in mol m-2 timestep-1
            CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                             fleach_nh4 * nh4_solute(:,:),transport_nh4_solute(:,:),leaching_nh4_solute(:))
            CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                             no3_solute(:,:),transport_no3_solute(:,:),leaching_no3_solute(:))
            CALL calc_liquid_phase_transport(options, soil_config%nsoil_water, percolation_sl(:,:), soil_depth_sl(:,:), &
                                             fleach_po4 * po4_solute(:,:),transport_po4_solute(:,:),leaching_po4_solute(:))
          ENDIF
       ENDIF
    ENDIF 
  
  END SUBROUTINE calc_potential_nloss_rate
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_potential_decomposition_rate, Calculate decomposition of pools
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_potential_decomposition_rate(ics, &
                                                ice, &
                                                nc, &
                                                iblk, &
                                                nsoil_water, &
                                                rtm_depolymerisation_act, &
                                                rmm_depolymerisation_act, &
                                                sb_loss, &
                                                sb_pool )

    USE mo_base_constants,          ONLY: dtime
    USE mo_sb_memory_class,         ONLY: t_slm_struct_sb_pool
    USE mo_sb_constants,            ONLY: tau_non_woody_litter, tau_woody_litter, tau_fast, tau_slow

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER,                              INTENT(in)    :: ics, ice, nc, iblk         !< chunks [ics:ice] and block numbers
    INTEGER,                              INTENT(in)    :: nsoil_water
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: rtm_depolymerisation_act, &
                                                           rmm_depolymerisation_act
    TYPE(t_slm_struct_sb_pool),           INTENT(inout) :: sb_loss, &
                                                           sb_pool
    ! 0.2 Local
    ! 
    
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 3.0 Calculate decomposition of pools: dP/dt = - f(T)*f(M)/tau_pool * P
    !!
    !> 3.1 non-woody litter (= polymeric pool)
    !! 
    sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                          sb_pool%polymeric_litter%carbon%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                            sb_pool%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                              rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                              sb_pool%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                            sb_pool%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                            rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                            sb_pool%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                              rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_non_woody_litter * dtime * &
                                                              sb_pool%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.2 woody litter
    !! 
    sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) + &
                                                      rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                      sb_pool%woody_litter%carbon%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                          sb_pool%woody_litter%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) + &
                                                        rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                        sb_pool%woody_litter%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) + &
                                                          rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_woody_litter * dtime * &
                                                          sb_pool%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.3 fast SOM pool (=microbial pool)
    !!
    sb_loss%microbial%carbon%ptr(ics:ice,:,iblk) = sb_loss%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                   sb_pool%microbial%carbon%ptr(ics:ice,:,iblk)
    sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                       rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                       sb_pool%microbial%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%microbial%carbon13%ptr(ics:ice,:,iblk) = sb_loss%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%microbial%carbon14%ptr(ics:ice,:,iblk) = sb_loss%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                     sb_pool%microbial%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                       rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * &
                                                       sb_pool%microbial%nitrogen15%ptr(ics:ice,:,iblk)

    !> 3.4 slow SOM pool (=residue pool)
    !!
    sb_loss%residue%carbon%ptr(ics:ice,:,iblk) = sb_loss%residue%carbon%ptr(ics:ice,:,iblk) + & 
                                                 rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                 sb_pool%residue%carbon%ptr(ics:ice,:,iblk)
    sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) = sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                   sb_pool%residue%nitrogen%ptr(ics:ice,:,iblk)
    sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                     sb_pool%residue%phosphorus%ptr(ics:ice,:,iblk)
    sb_loss%residue%carbon13%ptr(ics:ice,:,iblk) = sb_loss%residue%carbon13%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                   sb_pool%residue%carbon13%ptr(ics:ice,:,iblk)
    sb_loss%residue%carbon14%ptr(ics:ice,:,iblk) = sb_loss%residue%carbon14%ptr(ics:ice,:,iblk) + & 
                                                   rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                   sb_pool%residue%carbon14%ptr(ics:ice,:,iblk)
    sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) = sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                     rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_slow * dtime * &
                                                     sb_pool%residue%nitrogen15%ptr(ics:ice,:,iblk)  
  
  END SUBROUTINE calc_potential_decomposition_rate

  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_actual_flux_rate, formation of SOM resulting from the decomposition
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_actual_flux_rate(ics, &
                                   ice, &
                                   nc, &
                                   iblk, &
                                   nsoil_water, &
                                   sb_pool_residue_phosphorus, &
                                   growth_c, &
                                   growth_n, &
                                   mic_uptake_nh4_sl, &
                                   plant_uptake_nh4_sl, &
                                   mycorrhiza_uptake_nh4_sl, &
                                   plant_uptake_no3_sl, &
                                   mycorrhiza_uptake_no3_sl, &
                                   mic_uptake_po4_sl, &
                                   plant_uptake_po4_sl, &
                                   mycorrhiza_uptake_po4_sl, &
                                   nitrification_no3, &
                                   denitrification_n2, &
                                   nh4_solute, &
                                   no3_solute, &
                                   po4_solute, &
                                   transport_nh4_solute, &
                                   transport_no3_solute, &
                                   transport_po4_solute, &
                                   apparent_mic_cue, &
                                   biomineralisation_po4, &
                                   rtm_depolymerisation_act, &
                                   rmm_depolymerisation_act, &
                                   sb_loss  )
                                                
    USE mo_base_constants,          ONLY: dtime
    USE mo_sb_memory_class,         ONLY: t_slm_struct_sb_pool
    USE mo_sb_constants,            ONLY: microbial_cue_max, microbial_cue_min, microbial_cn, microbial_np, tau_fast

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER,                              INTENT(in)    :: ics, ice, nc, iblk         !< chunks [ics:ice] and block numbers
    INTEGER,                              INTENT(in)    :: nsoil_water
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: sb_pool_residue_phosphorus, &
                                                           growth_c, &
                                                           growth_n, &
                                                           mic_uptake_nh4_sl, &
                                                           plant_uptake_nh4_sl, &
                                                           mycorrhiza_uptake_nh4_sl, &
                                                           plant_uptake_no3_sl, &
                                                           mycorrhiza_uptake_no3_sl, &
                                                           mic_uptake_po4_sl, &
                                                           plant_uptake_po4_sl, &
                                                           mycorrhiza_uptake_po4_sl, &
                                                           nitrification_no3, &
                                                           denitrification_n2, &
                                                           nh4_solute, &
                                                           no3_solute, &
                                                           po4_solute, &
                                                           transport_nh4_solute, &
                                                           transport_no3_solute, &
                                                           transport_po4_solute, &
                                                           apparent_mic_cue, &
                                                           biomineralisation_po4, &
                                                           rtm_depolymerisation_act, &
                                                           rmm_depolymerisation_act
    TYPE(t_slm_struct_sb_pool),           INTENT(inout) :: sb_loss
    ! 0.2 Local
    REAL(wp), DIMENSION(nc,nsoil_water)                 :: hlp1, hlp2
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------
  
    !> 4.0 formation of SOM resulting from the decomposition dPrec/dt = CUE_Pdon * dPdon/dt
    !! net mineralisation follows the stoichiometric difference of the donor and recipient pool, taking account for CUE
    !! for fast -> slow -> fast this only causes positive mineralisation because of the C:N:P and CUE
    !! for litter -> fast this causes in most circumstances negative mineralisation, which is then taken from the 
    !!               mineral pools (or organic matter for P)

    !> 4.1 formation of fast SOM from litter 
    !! 
    ! calculate C growth of fast pool
    growth_c(:,:) = microbial_cue_max * & 
                    (sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + & 
                     sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk))
    
    ! calculate potential microbial immobilisation rate as difference between maximal growth rate and N in litter 
    mic_uptake_nh4_sl(:,:) = growth_c(:,:) / microbial_cn - & 
                             (sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                             sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk))

    ! limit microbial and plant uptake, as well as nitrification to available NH4 minus transport losses
    WHERE(mic_uptake_nh4_sl(:,:)>=0.0_wp)
       hlp1(:,:) = mic_uptake_nh4_sl(:,:) + plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:)  + mycorrhiza_uptake_nh4_sl(:,:)
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
       hlp1(:,:) = plant_uptake_nh4_sl(:,:) + nitrification_no3(:,:) 
       hlp2(:,:) = nh4_solute(:,:) + transport_nh4_solute(:,:) - mic_uptake_nh4_sl(:,:)
       WHERE(hlp2(:,:) < hlp1(:,:)) 
           WHERE(hlp1(:,:) > eps8)
              plant_uptake_nh4_sl(:,:) = hlp2(:,:) * plant_uptake_nh4_sl(:,:) / hlp1(:,:)
              nitrification_no3(:,:)   = hlp2(:,:) * nitrification_no3(:,:) / hlp1(:,:)
           ELSEWHERE
              plant_uptake_nh4_sl(:,:) = 0.0_wp 
              nitrification_no3(:,:) = 0.0_wp 
           ENDWHERE 
        ENDWHERE
    ENDWHERE 
   
    ! limit plant uptake and denitrification to available NO3 minus transport losses
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

    ! calculate C growth based on N availability
    growth_n(:,:) = microbial_cn * (sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + & 
                                    sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) + mic_uptake_nh4_sl(:,:))

    ! adjust carbon-use efficiency to match N availability by unchanged litter decomposition
    WHERE(growth_c(:,:) <= growth_n(:,:)) ! enough N is available -> constant CUE
       apparent_mic_cue(:,:) = microbial_cue_max
    ELSEWHERE(growth_c(:,:) > eps8) ! not enough N is available -> changed CUE (within min/max CUE)
       apparent_mic_cue(:,:) = MAX(microbial_cue_min, & 
                                   growth_n(:,:) / ( growth_c(:,:) / microbial_cue_max ))
    ELSEWHERE ! growth C is too small -> CUE does not matter
       apparent_mic_cue(:,:) = microbial_cue_min 
    ENDWHERE

    ! calculate actual litter loss rate, calculated such that litter decomposition
    ! matches the actual immobilisation flux and CN uptake from litter given actual CUE
    hlp1(:,:) = (apparent_mic_cue(:,:) / microbial_cn * &
                      (sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) - &
                      (sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)))
    hlp2(:,:) = 1.0_wp
    WHERE(mic_uptake_nh4_sl(:,:)<hlp1(:,:))
       WHERE(hlp1(:,:)>eps8)
           hlp2(:,:) = mic_uptake_nh4_sl(:,:) / hlp1(:,:)
       ELSEWHERE
           hlp2(:,:) = 0.0_wp
       ENDWHERE
       sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk)
       sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)

       sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)
       sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)
       sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk)
       sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk)
       sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)
       sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) = hlp2(:,:) * sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)
    ENDWHERE 

    ! calculate potential P biomineralisation from slow organic P pool (mol /m3 / timestep)
    biomineralisation_po4(:,:) = rtm_depolymerisation_act(:,:)*rmm_depolymerisation_act(:,:)/tau_fast * dtime * & 
                                 sb_pool_residue_phosphorus(:,:) 

    ! calculate potential microbial P uptake from soluable pool (mol /m3 / timestep) 
    mic_uptake_po4_sl(:,:) =  apparent_mic_cue(:,:) / microbial_cn / microbial_np * & 
                              (sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + &
                               sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) - & 
                              (sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + & 
                               sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk))

    ! constrain P uptake by microbes and plants to the available sum of P from P_soluable, biomineralisation and transport (mol /m3 / timestep)   
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

    ! limit biomineralisation to the actual rate of biomineralisation required to meet microbial and plant demands
    ! and substract from slow=residue pool 
    biomineralisation_po4(:,:) = MAX(mic_uptake_po4_sl(:,:) + plant_uptake_po4_sl(:,:) - & 
                                     (po4_solute(:,:) + transport_po4_solute(:,:)),0.0_wp)
    sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) = sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                     biomineralisation_po4(:,:)
  
  END SUBROUTINE calc_actual_flux_rate
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_actual_SOM_formation, formation of SOM and associated respiration, net mineralisation fluxes
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_actual_SOM_formation(ics, &
                                       ice, &
                                       nc, &
                                       iblk, &
                                       nsoil_water, &
                                       apparent_mic_cue, &
                                       mic_uptake_nh4_sl, &
                                       mic_uptake_po4_sl, &
                                       het_respiration, &
                                       het_respiration_c13, &
                                       het_respiration_c14, &
                                       net_mineralisation_nh4, &
                                       net_mineralisation_po4, &
                                       net_mineralisation_nh4_n15, &
                                       sb_loss, &
                                       sb_formation  )
  
    USE mo_sb_memory_class,         ONLY: t_slm_struct_sb_pool
    USE mo_sb_constants,            ONLY: residue_cue, microbial_cue_min, microbial_np, microbial_cn

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER,                              INTENT(in)    :: ics, ice, nc, iblk         !< chunks [ics:ice] and block numbers
    INTEGER,                              INTENT(in)    :: nsoil_water
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: apparent_mic_cue, &
                                                           mic_uptake_nh4_sl, &
                                                           mic_uptake_po4_sl, &
                                                           het_respiration, &
                                                           het_respiration_c13, &
                                                           het_respiration_c14, &
                                                           net_mineralisation_nh4, &
                                                           net_mineralisation_po4, &
                                                           net_mineralisation_nh4_n15
    TYPE(t_slm_struct_sb_pool),           INTENT(inout) :: sb_loss, &
                                                           sb_formation
    ! 0.2 Local
    ! 
    
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

    !> 5.0 formation of SOM and associated respiration, net mineralisation fluxes
    !!
    !> 5.1 formation of fast SOM from litter 
    !! 
    sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                        apparent_mic_cue(:,:) * & 
                                                        (sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                          mic_uptake_nh4_sl(:,:) 
    sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                            (sb_loss%polymeric_litter%phosphorus%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%phosphorus%ptr(ics:ice,:,iblk)) + & 
                                                            mic_uptake_po4_sl(:,:)
    sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                          apparent_mic_cue(:,:) * & 
                                                          (sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk)) 
    sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                          apparent_mic_cue(:,:) * & 
                                                          (sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)) 
!! to be moved to N15 calculation
    WHERE(sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                              (sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)+sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)) / & 
                              (sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)) * & 
                               apparent_mic_cue(:,:) / microbial_cn * & 
                              (sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    ENDWHERE

    het_respiration(:,:)     = het_respiration(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                               (sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    het_respiration_c13(:,:) = het_respiration_c13(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                               (sb_loss%polymeric_litter%carbon13%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon13%ptr(ics:ice,:,iblk)) 
    het_respiration_c14(:,:) = het_respiration_c14(:,:) + (1._wp - apparent_mic_cue(:,:)) * & 
                               (sb_loss%polymeric_litter%carbon14%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon14%ptr(ics:ice,:,iblk)) 
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) - mic_uptake_nh4_sl(:,:)  
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) - mic_uptake_po4_sl(:,:)  
    WHERE(sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       net_mineralisation_nh4_n15(:,:)     = net_mineralisation_nh4_n15(:,:) + & 
                              sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk) - & 
                              (sb_loss%polymeric_litter%nitrogen15%ptr(ics:ice,:,iblk)+sb_loss%woody_litter%nitrogen15%ptr(ics:ice,:,iblk)) / & 
                              (sb_loss%polymeric_litter%nitrogen%ptr(ics:ice,:,iblk)+sb_loss%woody_litter%nitrogen%ptr(ics:ice,:,iblk)) * & 
                               apparent_mic_cue(:,:) / microbial_cn * & 
                              (sb_loss%polymeric_litter%carbon%ptr(ics:ice,:,iblk) + sb_loss%woody_litter%carbon%ptr(ics:ice,:,iblk)) 
    ENDWHERE

    !> 5.2 formation of slow SOM (=residue pool) from fast SOM decay
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

    !> 5.3 formation of fast SOM (=microbial pool) from slow SOM decay 
    !! 
    sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon%ptr(ics:ice,:,iblk) + & 
                                                        residue_cue * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%microbial%phosphorus%ptr(ics:ice,:,iblk) + & 
                                                            residue_cue / microbial_cn / microbial_np * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon13%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue * sb_loss%residue%carbon13%ptr(ics:ice,:,iblk)
    sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) = sb_formation%microbial%carbon14%ptr(ics:ice,:,iblk) + & 
                                                          residue_cue * sb_loss%residue%carbon14%ptr(ics:ice,:,iblk)
    WHERE(sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) = sb_formation%microbial%nitrogen15%ptr(ics:ice,:,iblk) + & 
                                                               sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) / &  
                                                               sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) * &  
                                                               residue_cue / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE

    het_respiration(:,:)      = het_respiration(:,:)     + (1._wp - residue_cue) * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    het_respiration_c13(:,:)  = het_respiration_c13(:,:) + (1._wp - residue_cue) * sb_loss%residue%carbon13%ptr(ics:ice,:,iblk)
    het_respiration_c14(:,:)  = het_respiration_c14(:,:) + (1._wp - residue_cue) * sb_loss%residue%carbon14%ptr(ics:ice,:,iblk)
    net_mineralisation_nh4(:,:) = net_mineralisation_nh4(:,:) + sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    net_mineralisation_po4(:,:) = net_mineralisation_po4(:,:) + sb_loss%residue%phosphorus%ptr(ics:ice,:,iblk) - &
                                                                residue_cue / microbial_cn / microbial_np * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    WHERE(sb_loss%microbial%nitrogen%ptr(ics:ice,:,iblk)>eps8)
       net_mineralisation_nh4_n15(:,:) = net_mineralisation_nh4_n15(:,:) + & 
                                    sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) - &
                                    sb_loss%residue%nitrogen15%ptr(ics:ice,:,iblk) / sb_loss%residue%nitrogen%ptr(ics:ice,:,iblk) * &  
                                    residue_cue / microbial_cn * sb_loss%residue%carbon%ptr(ics:ice,:,iblk)
    ENDWHERE  
    
  END SUBROUTINE calc_actual_SOM_formation  
  
  
  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_actual_nloss_rate, calculate N losses and nitrification
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_actual_nloss_rate(nc, &
                                    nsoil_water, &
                                    soil_depth_sl, &
                                    net_mineralisation_nh4, &
                                    volatilisation_n, &
                                    nitrification_no3, &
                                    nitrification_noy, &
                                    nitrification_n2o, &
                                    denitrification_noy, &
                                    denitrification_n2o, &
                                    denitrification_n2, &
                                    noy, &
                                    n2o, &
                                    n2, &
                                    transport_noy, &
                                    transport_n2o, &
                                    transport_n2, &
                                    transport_nh4_solute, &
                                    transport_no3_solute, &
                                    transport_po4_solute, &
                                    rtm_denitrification_act, &
                                    rtm_gasdiffusion_act, &
                                    rmm_gasdiffusion_act, &
                                    emission_noy, &
                                    emission_n2o, &
                                    emission_n2, &
                                    leaching_nh4_solute, &
                                    leaching_no3_solute, &
                                    leaching_po4_solute )

    USE mo_sb_config_class,         ONLY: sb_config
    USE mo_base_constants,          ONLY: dtime, one_day
    USE mo_sb_constants,            ONLY: floss_nmin, fnit_nmin, f_nit_noy, f_nit_n2o, f_denit_noy, f_denit_n2o
  
    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER,                              INTENT(in)    :: nc                         !< number of chunks [ics:ice]
    INTEGER,                              INTENT(in)    :: nsoil_water                !< number of soil layers
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: net_mineralisation_nh4, &
                                                           volatilisation_n, &
                                                           nitrification_no3, &
                                                           nitrification_noy, &
                                                           nitrification_n2o, &
                                                           denitrification_noy, &
                                                           denitrification_n2o, &
                                                           denitrification_n2, &
                                                           noy, &
                                                           n2o, &
                                                           n2, &
                                                           transport_noy, &
                                                           transport_n2o, &
                                                           transport_n2, &
                                                           transport_nh4_solute, &
                                                           transport_no3_solute, &
                                                           transport_po4_solute, &
                                                           rtm_denitrification_act, &
                                                           rtm_gasdiffusion_act, &
                                                           rmm_gasdiffusion_act, &
                                                           soil_depth_sl
    REAL(wp), DIMENSION(nc),              INTENT(inout) :: emission_noy, &
                                                           emission_n2o, &
                                                           emission_n2, &
                                                           leaching_nh4_solute, &
                                                           leaching_no3_solute, &
                                                           leaching_po4_solute
    ! 0.2 Local
    ! 
    
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------

  
    !> 6.0 calculate N losses and nitrification
    !! 
    IF(sb_config%sb_nloss_scheme == 'fixed') THEN
       !> 6.1 case of fixed N losses: N loss proportional to net mineralisation
       WHERE(net_mineralisation_nh4(:,:)>0.0_wp)
          ! N lost to atmosphere during processing of net mineralisation [mol m-3 timestep-1]
          volatilisation_n(:,:)  = net_mineralisation_nh4(:,:) * floss_nmin 
          ! NH4 transformed to NO3 [mol m-3 timestep-1]
          nitrification_no3(:,:) = net_mineralisation_nh4(:,:) * (1._wp - floss_nmin) * fnit_nmin 
       ENDWHERE
       ! report volatilisation losses as N2 emission [mol m-2 timestep-1]
       emission_n2(:) = SUM(volatilisation_n(:,:) * soil_depth_sl(:,:),DIM=2)
    ELSE
       !> 6.2 case of dynamic N loss calculations, following Zaehle et al. 2011
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
       IF(sb_config%sb_model_scheme=='simple_0d')THEN
          leaching_nh4_solute(:) = -transport_nh4_solute(:,1) * soil_depth_sl(:,1)
          leaching_no3_solute(:) = -transport_no3_solute(:,1) * soil_depth_sl(:,1)
          leaching_po4_solute(:) = -transport_po4_solute(:,1) * soil_depth_sl(:,1)
       ENDIF
    ENDIF
 
  END SUBROUTINE calc_actual_nloss_rate


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_litter_partitioning, Update soil litter pools with litter fall
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_litter_partitioning(ics, &
                                      ice, &
                                      nc, &
                                      iblk, &
                                      nsoil_water, &
                                      soil_depth_sl, &
                                      root_fraction_act_sl, &
                                      veg_litterfall, &
                                      sb_formation )

    USE mo_sb_memory_class,         ONLY: t_slm_struct_sb_pool
    USE mo_veg_memory_class,        ONLY: t_slm_struct_veg_pool

    IMPLICIT NONE
    ! 0.1 InOut
    INTEGER,                              INTENT(in)    :: ics, ice, nc, iblk         !< number of chunks and block number
    INTEGER,                              INTENT(in)    :: nsoil_water                !< number of soil layers
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: soil_depth_sl, &           !< depth per soil layer
                                                           root_fraction_act_sl       !< root fraction per soil layer
    TYPE(t_slm_struct_veg_pool),          INTENT(inout) :: veg_litterfall             !< vegetation litterfall pool
    TYPE(t_slm_struct_sb_pool),           INTENT(inout) :: sb_formation               !< soil biogeochemistry formation pool
    ! 0.2 Local
    INTEGER                                             :: isoil                      !< loop over soil layers
    
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 7.0 Update soil litter pools with litter fall
    !!
    sb_formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon%ptr(ics:ice,1,iblk) + & 
       (veg_litterfall%leaf%carbon%ptr + veg_litterfall%fruit%carbon%ptr + & 
        veg_litterfall%labile%carbon%ptr + veg_litterfall%reserve%carbon%ptr + & 
        veg_litterfall%seed_bed%carbon%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,1,iblk) +  & 
       (veg_litterfall%leaf%nitrogen%ptr + veg_litterfall%fruit%nitrogen%ptr + & 
        veg_litterfall%labile%nitrogen%ptr + veg_litterfall%reserve%nitrogen%ptr + & 
        veg_litterfall%seed_bed%nitrogen%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,1,iblk) + & 
       (veg_litterfall%leaf%phosphorus%ptr + veg_litterfall%fruit%phosphorus%ptr + & 
        veg_litterfall%labile%phosphorus%ptr + veg_litterfall%reserve%phosphorus%ptr + & 
        veg_litterfall%seed_bed%phosphorus%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon13%ptr(ics:ice,1,iblk) + & 
       (veg_litterfall%leaf%carbon13%ptr + veg_litterfall%fruit%carbon13%ptr + & 
        veg_litterfall%labile%carbon13%ptr + veg_litterfall%reserve%carbon13%ptr + & 
        veg_litterfall%seed_bed%carbon13%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%carbon14%ptr(ics:ice,1,iblk) + & 
       (veg_litterfall%leaf%carbon14%ptr + veg_litterfall%fruit%carbon14%ptr + & 
        veg_litterfall%labile%carbon14%ptr + veg_litterfall%reserve%carbon14%ptr + & 
        veg_litterfall%seed_bed%carbon14%ptr ) / soil_depth_sl(:,1)
    sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) = sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,1,iblk) + &
       (veg_litterfall%leaf%nitrogen15%ptr + veg_litterfall%fruit%nitrogen15%ptr + & 
        veg_litterfall%labile%nitrogen15%ptr + veg_litterfall%reserve%nitrogen15%ptr + & 
        veg_litterfall%seed_bed%nitrogen15%ptr ) / soil_depth_sl(:,1)

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

    DO isoil = 1,nsoil_water
       sb_formation%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk)     = sb_formation%polymeric_litter%carbon%ptr(ics:ice,isoil,iblk) + & 
          ( veg_litterfall%fine_root%carbon%ptr     + veg_litterfall%coarse_root%carbon%ptr )     * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk)   = sb_formation%polymeric_litter%nitrogen%ptr(ics:ice,isoil,iblk) + & 
          ( veg_litterfall%fine_root%nitrogen%ptr   + veg_litterfall%coarse_root%nitrogen%ptr )   * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) = sb_formation%polymeric_litter%phosphorus%ptr(ics:ice,isoil,iblk) + & 
          ( veg_litterfall%fine_root%phosphorus%ptr + veg_litterfall%coarse_root%phosphorus%ptr ) * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk)   = sb_formation%polymeric_litter%carbon13%ptr(ics:ice,isoil,iblk) + & 
          ( veg_litterfall%fine_root%carbon13%ptr   + veg_litterfall%coarse_root%carbon13%ptr )   * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk)   = sb_formation%polymeric_litter%carbon14%ptr(ics:ice,isoil,iblk) + & 
          ( veg_litterfall%fine_root%carbon14%ptr   + veg_litterfall%coarse_root%carbon14%ptr )   * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
       sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) = sb_formation%polymeric_litter%nitrogen15%ptr(ics:ice,isoil,iblk) + & 
          ( veg_litterfall%fine_root%nitrogen15%ptr + veg_litterfall%coarse_root%nitrogen15%ptr ) * & 
          root_fraction_act_sl(:,isoil) / soil_depth_sl(:,isoil) 
    END DO  

  END SUBROUTINE calc_litter_partitioning


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task to calc_sb_fluxes_simple
  !
  !-----------------------------------------------------------------------------------------------------
  !> calc_mycorrhiza_dynamics, nutrient uptake and flux, and litter formation
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE calc_mycorrhiza_dynamics(lctlib, &
                                      ics, &
                                      ice, &
                                      nc, &
                                      iblk, &
                                      nsoil_water, &
                                      soil_depth_sl, &
                                      mycorrhiza_uptake_nh4_sl, &
                                      mycorrhiza_uptake_no3_sl, &
                                      mycorrhiza_uptake_po4_sl, &
                                      fresp, &
                                      transform_respiration, &
                                      het_respiration, &
                                      het_respiration_c13, &
                                      het_respiration_c14, &
                                      fturn_mycorrhiza, &
                                      uptake_norg, &
                                      uptake_porg, &
                                      sb_formation, &
                                      sb_pool, &
                                      sb_loss )

    USE mo_base_constants,          ONLY: dtime, one_day, one_year
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element    
    USE mo_sb_memory_class,         ONLY: t_slm_struct_sb_pool
    USE mo_veg_constants,           ONLY: transform_cost_nh4, transform_cost_no3
    USE mo_sb_constants,            ONLY: cn_mycorrhiza, np_mycorrhiza

    IMPLICIT NONE
    ! 0.1 InOut
    TYPE(t_lctlib_element),               INTENT(in)    :: lctlib                       !< land-cover-type library - parameter across pft's
    INTEGER,                              INTENT(in)    :: ics, ice, nc, iblk           !< number of chunks and block number
    INTEGER,                              INTENT(in)    :: nsoil_water                  !< number of soil layers
    REAL(wp), DIMENSION(nc,nsoil_water),  INTENT(inout) :: soil_depth_sl, &             !< depth per soil layer
                                                           mycorrhiza_uptake_nh4_sl, &
                                                           mycorrhiza_uptake_no3_sl, &
                                                           mycorrhiza_uptake_po4_sl, &
                                                           fresp, &
                                                           transform_respiration, &
                                                           het_respiration, &
                                                           het_respiration_c13, &
                                                           het_respiration_c14, &
                                                           fturn_mycorrhiza
    REAL(wp), DIMENSION(nc),              INTENT(inout) :: uptake_norg, &
                                                           uptake_porg
    TYPE(t_slm_struct_sb_pool),           INTENT(inout) :: sb_formation, &            !< soil biogeochemistry formation pool
                                                           sb_pool, &                 !< soil biogeochemistry pool
                                                           sb_loss                    !< soil biogeochemistry loss pool
    ! 0.2 Local
    INTEGER                                             :: isoil                      !< loop over soil layers
    REAL(wp), DIMENSION(nc,nsoil_water)                 :: hlp1, hlp2, hlp3
    
    
    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


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
    hlp1(:,:) = sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) / cn_mycorrhiza - &
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

    uptake_norg(:) = 0.0_wp
    DO isoil = 1, nsoil_water
         uptake_norg(:) = uptake_norg(:) + hlp3(:,isoil) * soil_depth_sl(:,isoil)
    ENDDO 

    ! P demand of mycorrhiza in mol P / m3 / timestep ! BLARPP SZ: constrained to positive solutions
    hlp1(:,:) = MAX(sb_pool%mycorrhiza%carbon%ptr(ics:ice,:,iblk) / cn_mycorrhiza / np_mycorrhiza - &
                    sb_pool%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk), 0.0_wp)
    ! incorporation of P into mycorrhiza (mol P / m3 / timestep)
    ! constrain flux to minimum of potential uptake and demand 
    hlp2(:,:) = MIN(hlp1(:,:), mycorrhiza_uptake_po4_sl(:,:))
    sb_formation%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) = sb_formation%mycorrhiza%phosphorus%ptr(ics:ice,:,iblk) + &
                                                             hlp2(:,:)

    ! remainder is exported to plants (mol P / m2 / timestep)
    hlp3(:,:) = mycorrhiza_uptake_po4_sl(:,:) - hlp2(:,:)
    uptake_porg(:) = 0.0_wp
    DO isoil = 1, nsoil_water
         uptake_porg(:) = uptake_porg(:) + hlp3(:,isoil) * soil_depth_sl(:,isoil) 
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
!!    USE mo_sb_config_class,       ONLY: sb_config
!!    USE mo_sb_constants
!!    USE mo_sb_jsm_transport
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
!!                                             rtm_depolymerisation_act, &
!!                                             rmm_depolymerisation_act, & 
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
!!    ! Get process config
!!!    dsl4jsb_Get_config(xx_)
!!    ! Get process memories
!!    dsl4jsb_Get_memory(SOIL_)
!!    dsl4jsb_Get_memory(SB_)
!!    ! Get process variables (Set pointers to variables in memory)
!!    model  => Get_model(tile%owner_model_id)
!!    lctlib => model%lctlib(tile%lcts%lib_id)
!!    ! ---------------------------
!!    dsl4jsb_Get_var3D_onChunk(SOIL_,     percolation_sl)
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
!!    !> 1.0 potential 15N fluxes involving ammonium
!!    !!
!!    IF(sb_config%sb_model_scheme == 'simple_0d')THEN
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
!!    IF(sb_config%sb_model_scheme == 'simple_0d')THEN
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
